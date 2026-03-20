"""
Preprocessor: gene filtering, tissue/sex stratification, PCA-based outlier detection.

Used before any clock training or prediction.
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class Preprocessor:
    """Filter genes and stratify Atlas samples by tissue and sex.

    Parameters
    ----------
    min_count : int
        Minimum total raw count across all samples for a gene to be kept.
        Genes with rowSum <= min_count are removed.
    """

    def __init__(self, min_count: int = 1):
        self.min_count = min_count

    # ------------------------------------------------------------------
    # Gene filtering
    # ------------------------------------------------------------------

    def filter_genes(self, counts: pd.DataFrame) -> pd.DataFrame:
        """Remove genes with very low/zero total expression.

        Parameters
        ----------
        counts : pd.DataFrame
            Raw counts, genes × samples.

        Returns
        -------
        pd.DataFrame
            Filtered counts.
        """
        mask = counts.sum(axis=1) > self.min_count
        filtered = counts.loc[mask]
        n_removed = counts.shape[0] - filtered.shape[0]
        print(f"Gene filter: kept {filtered.shape[0]:,} / {counts.shape[0]:,} genes "
              f"(removed {n_removed:,} with total count ≤ {self.min_count})")
        return filtered

    # ------------------------------------------------------------------
    # Stratification
    # ------------------------------------------------------------------

    def stratify(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        tissue: str,
        sex: str | None = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Subset counts and metadata to a specific tissue (and optionally sex).

        Parameters
        ----------
        counts : pd.DataFrame
            Expression matrix, genes × samples.
        metadata : pd.DataFrame
            Sample metadata indexed by sample ID.
        tissue : str
            Tissue name to select (case-insensitive match to metadata 'tissue' column).
        sex : str or None
            'M', 'F', or None (combined).

        Returns
        -------
        (counts_sub, meta_sub) : tuple of DataFrames
            Filtered expression matrix and metadata.
        """
        mask = metadata["tissue"].str.lower() == tissue.lower()
        if sex is not None:
            mask &= metadata["sex"].str.upper() == sex.upper()

        meta_sub = metadata.loc[mask]
        shared = counts.columns.intersection(meta_sub.index)
        meta_sub = meta_sub.loc[shared]
        counts_sub = counts[shared]

        label = f"{tissue} / {sex}" if sex else tissue
        print(f"Stratified [{label}]: {counts_sub.shape[1]} samples, "
              f"age range {meta_sub['age_days'].min()}–{meta_sub['age_days'].max()} days")
        return counts_sub, meta_sub

    # ------------------------------------------------------------------
    # Outlier detection (PCA-based, per method.txt)
    # ------------------------------------------------------------------

    def detect_outliers(
        self,
        counts: pd.DataFrame,
        n_sd: float = 2.0,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Remove samples whose PC1 score is > n_sd standard deviations from the mean.

        Used before BayesAge2 and PCR modeling (not EN).

        Parameters
        ----------
        counts : pd.DataFrame
            Expression matrix, genes × samples. Already filtered/normalized.
        n_sd : float
            Number of standard deviations to use as the outlier threshold.

        Returns
        -------
        (counts_clean, outlier_ids) : tuple
            Expression matrix with outliers removed, and list of removed sample IDs.
        """
        # PCA on samples (samples × genes)
        mat = counts.T.values
        scaler = StandardScaler()
        mat_scaled = scaler.fit_transform(mat)

        pca = PCA(n_components=1, random_state=1)
        pc1 = pca.fit_transform(mat_scaled)[:, 0]

        mean_pc1 = pc1.mean()
        sd_pc1 = pc1.std()
        threshold = n_sd * sd_pc1

        outlier_mask = np.abs(pc1 - mean_pc1) > threshold
        sample_ids = counts.columns.tolist()
        outlier_ids = [sid for sid, is_out in zip(sample_ids, outlier_mask) if is_out]

        counts_clean = counts.drop(columns=outlier_ids)

        pct_var = pca.explained_variance_ratio_[0] * 100
        print(f"Outlier detection: PC1 explains {pct_var:.1f}% variance | "
              f"threshold = ±{threshold:.3f} | "
              f"removed {len(outlier_ids)} outlier(s): {outlier_ids}")
        return counts_clean, outlier_ids
