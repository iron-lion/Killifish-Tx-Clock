"""
BayesAge 2.0 aging clock implementation.

Normalization : frequency count normalization (raw / total reads per sample).
Training      : LOWESS regression per gene; Spearman correlation with age.
Prediction    : Poisson log-likelihood maximization across age range.
Validation    : Leave-One-Sample-Out Cross-Validation (LOSO-CV).

Optimization  : Vectorized Spearman (rank-correlation via matrix ops);
                LOWESS computed only on top-ranked genes (up to `lowess_top_n`).
"""

import math
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import rankdata
from statsmodels.nonparametric.smoothers_lowess import lowess

from normalization import FrequencyNormalizer


# Age prediction range for killifish (days), per method.txt
AGE_MIN = 47
AGE_MAX = 163
AGE_STEP = 1
AGE_RANGE = np.arange(AGE_MIN, AGE_MAX + AGE_STEP, AGE_STEP)


class BayesAge2Clock:
    """BayesAge 2.0 aging clock.

    Parameters
    ----------
    age_range : array-like
        Candidate ages (in days) to evaluate during prediction.
        Default: 47–163 days at 1-day resolution.
    lowess_frac : float
        LOWESS smoothing fraction (frac parameter).
    lowess_top_n : int
        Only compute LOWESS for the top N genes by |Spearman r|.
        Reduces computation without affecting predictions (predictions
        already select the top M << lowess_top_n genes).
    """

    def __init__(
        self,
        age_range: np.ndarray = AGE_RANGE,
        lowess_frac: float = 0.3,
        lowess_top_n: int = 250,
    ):
        self.age_range = np.array(age_range)
        self.lowess_frac = lowess_frac
        self.lowess_top_n = lowess_top_n
        self._reference: pd.DataFrame | None = None
        self._freq_normalizer = FrequencyNormalizer()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _vectorized_spearman(mat: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Compute Spearman correlation of each row in mat with vector y.

        Parameters
        ----------
        mat : np.ndarray, shape (n_genes, n_samples)
        y   : np.ndarray, shape (n_samples,)

        Returns
        -------
        np.ndarray, shape (n_genes,) — Spearman r values (NaN for constant rows).
        """
        # Rank each gene across samples
        ranked_mat = np.apply_along_axis(rankdata, 1, mat)
        ranked_y = rankdata(y)

        rm = ranked_mat - ranked_mat.mean(axis=1, keepdims=True)
        ry = ranked_y - ranked_y.mean()

        num = rm @ ry
        denom = np.sqrt((rm ** 2).sum(axis=1) * (ry ** 2).sum())

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            rs = np.where(denom > 0, num / denom, np.nan)
        return rs

    def _build_single_reference(
        self,
        freq: pd.DataFrame,
        ages: pd.Series,
    ) -> pd.DataFrame:
        """Build reference matrix from frequency-normalized training data.

        Parameters
        ----------
        freq : pd.DataFrame
            Frequency-normalized counts, genes × samples.
        ages : pd.Series
            Age in days indexed by sample ID.

        Returns
        -------
        pd.DataFrame
            Genes × [age_47, ..., age_163, spearman_r, p_value]
            Only rows for the top `lowess_top_n` genes (by |Spearman r|)
            have meaningful LOWESS values; remaining genes have NaN LOWESS
            but valid spearman_r for reference sorting.
        """
        ages_arr = ages.loc[freq.columns].values.astype(float)
        mat = freq.values  # (n_genes, n_samples)
        age_cols = [f"age_{int(a)}" for a in self.age_range]

        # Vectorized Spearman for all genes
        spearman_r = self._vectorized_spearman(mat, ages_arr)

        # Select top genes by |r| for LOWESS (saves time)
        abs_r = np.where(np.isnan(spearman_r), 0.0, np.abs(spearman_r))
        top_idx = np.argsort(abs_r)[::-1][: self.lowess_top_n]
        top_mask = np.zeros(len(freq), dtype=bool)
        top_mask[top_idx] = True

        # LOWESS only on top genes
        lowess_matrix = np.full((len(freq), len(self.age_range)), np.nan)
        for i in np.where(top_mask)[0]:
            expr = mat[i]
            smoothed = lowess(expr, ages_arr, frac=self.lowess_frac, return_sorted=True)
            lowess_matrix[i] = np.interp(self.age_range, smoothed[:, 0], smoothed[:, 1])

        ref_data = np.hstack([
            lowess_matrix,
            spearman_r[:, None],
            np.zeros((len(freq), 1)),  # p-value placeholder (expensive to compute)
        ])
        ref_df = pd.DataFrame(
            ref_data,
            index=freq.index,
            columns=age_cols + ["spearman_r", "p_value"],
        )
        return ref_df

    def _predict_one_sample(
        self,
        raw_counts: pd.Series,
        reference: pd.DataFrame,
        n_genes: int,
    ) -> float:
        """Predict tAge for a single sample using Poisson log-likelihood.

        Parameters
        ----------
        raw_counts : pd.Series
            Observed raw counts (index = gene names).
        reference : pd.DataFrame
            Reference matrix with LOWESS-fit expected frequencies.
        n_genes : int
            Number of top |Spearman r| genes to use.

        Returns
        -------
        float
            Predicted transcriptomic age (tAge) in days.
        """
        # Select top M genes by |Spearman r| that also have LOWESS values
        ref_valid = reference.dropna(subset=["spearman_r"])
        ref_valid = ref_valid[ref_valid.iloc[:, :-2].notna().all(axis=1)]  # has LOWESS

        top_genes = (
            ref_valid["spearman_r"].abs().nlargest(n_genes).index
        )
        ref_sub = reference.loc[top_genes]
        age_cols = [f"age_{int(a)}" for a in self.age_range]
        expected = ref_sub[age_cols].values  # shape: (n_genes, n_ages)

        obs = raw_counts.reindex(top_genes).fillna(0).values.astype(float)

        # Replace zero/NaN expected with small value to avoid log(0)
        expected = np.where((expected <= 0) | np.isnan(expected), 1e-6, expected)

        # Scale expected frequencies by total observed reads to get expected counts
        total_obs = raw_counts.sum()
        lambda_mat = expected * total_obs  # shape: (n_genes, n_ages)

        # Poisson log-PMF: log P(obs_k | lambda) = k*log(lambda) - lambda - log(k!)
        # Floor per the original BayesAge2 implementation: when poisson.pmf underflows
        # to 0 (obs far from lambda), the original replaces PMF with 1e-3, equivalent
        # to clipping logP at log(1e-3). Without this, extreme outlier genes dominate
        # the argmax and force predictions to the ceiling age.
        _LOG_PMF_FLOOR = np.log(1e-3)  # = -6.908, matches original code behaviour
        obs_col = obs[:, None]
        log_factorial = np.array([math.lgamma(int(k) + 1) for k in obs])[:, None]
        log_probs = obs_col * np.log(lambda_mat) - lambda_mat - log_factorial
        log_probs = np.maximum(log_probs, _LOG_PMF_FLOOR)
        log_likelihood = log_probs.sum(axis=0)

        predicted_age = self.age_range[np.argmax(log_likelihood)]
        return float(predicted_age)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build_reference(
        self,
        raw_counts: pd.DataFrame,
        metadata: pd.DataFrame,
        save_path: Path | None = None,
    ) -> pd.DataFrame:
        """Build reference matrix from ALL training samples (for query prediction).

        Parameters
        ----------
        raw_counts : pd.DataFrame
            Raw counts, genes × samples.
        metadata : pd.DataFrame
            Must contain 'age_days' column, indexed by sample ID.
        save_path : Path, optional
            Save reference as tab-separated .tsv.

        Returns
        -------
        pd.DataFrame
            Reference matrix (stored in self._reference).
        """
        freq = self._freq_normalizer.normalize(raw_counts)
        ages = metadata["age_days"].astype(float)
        ref = self._build_single_reference(freq, ages)
        self._reference = ref

        if save_path is not None:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            ref.to_csv(save_path, sep="\t")
            print(f"Reference saved to {save_path}")

        return ref

    def loso_cv(
        self,
        raw_counts: pd.DataFrame,
        metadata: pd.DataFrame,
        m_values: list[int] | None = None,
    ) -> pd.DataFrame:
        """Leave-One-Sample-Out Cross-Validation.

        For each sample, builds a reference on the remaining N-1 samples,
        then predicts tAge for the left-out sample.

        Parameters
        ----------
        raw_counts : pd.DataFrame
            Raw counts, genes × samples.
        metadata : pd.DataFrame
            Must contain 'age_days', indexed by sample ID.
        m_values : list of int, optional
            Gene-set sizes to test. Default: 5, 10, ..., 50.

        Returns
        -------
        pd.DataFrame
            Columns: age_days, tAge_M5, tAge_M10, ..., indexed by sample_id.
        """
        if m_values is None:
            m_values = list(range(5, 55, 5))

        samples = raw_counts.columns.tolist()
        ages = metadata["age_days"].astype(float)
        freq_all = self._freq_normalizer.normalize(raw_counts)

        records = []
        for i, test_id in enumerate(samples):
            train_ids = [s for s in samples if s != test_id]
            freq_train = freq_all[train_ids]
            ages_train = ages.loc[train_ids]

            ref = self._build_single_reference(freq_train, ages_train)

            raw_test = raw_counts[test_id]
            row = {"sample_id": test_id, "age_days": ages[test_id]}
            for m in m_values:
                row[f"tAge_M{m}"] = self._predict_one_sample(raw_test, ref, n_genes=m)

            records.append(row)
            if (i + 1) % 5 == 0 or (i + 1) == len(samples):
                print(f"  LOSO-CV: {i+1}/{len(samples)} done")

        return pd.DataFrame(records).set_index("sample_id")

    def predict(
        self,
        raw_counts: pd.DataFrame,
        n_genes: int = 25,
        reference: pd.DataFrame | None = None,
    ) -> pd.Series:
        """Predict tAge for query samples using a pre-built reference.

        Parameters
        ----------
        raw_counts : pd.DataFrame
            Raw counts of query samples, genes × samples.
        n_genes : int
            Number of top Spearman-correlated genes to use.
        reference : pd.DataFrame, optional
            Overrides self._reference if provided.

        Returns
        -------
        pd.Series
            Predicted tAge per sample.
        """
        ref = reference if reference is not None else self._reference
        if ref is None:
            raise RuntimeError("No reference. Call build_reference() first.")

        return pd.Series(
            {sid: self._predict_one_sample(raw_counts[sid], ref, n_genes)
             for sid in raw_counts.columns},
            name=f"tAge_M{n_genes}",
        )

    @classmethod
    def load_reference(cls, path: Path) -> "BayesAge2Clock":
        """Load a saved reference matrix and return a ready-to-predict clock."""
        clock = cls()
        clock._reference = pd.read_csv(path, sep="\t", index_col=0)
        return clock
