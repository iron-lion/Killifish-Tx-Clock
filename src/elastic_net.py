"""
Elastic Net aging clock implementation.

Normalization : StandardScaler (z-scaling).
Training      : GridSearchCV over alpha × l1_ratio with LeaveOneOut CV.
Validation    : LOSO-CV via scikit-learn.
random_state  : 42.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import pearsonr


# Hyperparameter grids (per method.txt)
ALPHA_GRID = [1e-2, 1e-1, ]
#ALPHA_GRID = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0]
L1_RATIO_GRID = list(np.arange(0.0, 1.1, 0.2))

# max_iter per tissue (method.txt: brain uses 30,000; others 10,000)
MAX_ITER_DEFAULT = 10_000
MAX_ITER_BRAIN = 30_000
MAX_ITER_FINAL = 100_000

RANDOM_STATE = 42


class ElasticNetClock:
    """Elastic Net aging clock.

    Parameters
    ----------
    tissue : str
        Tissue name; used to set max_iter (brain gets higher limit).
    random_state : int
        Random seed for reproducibility.
    top_n_var_genes : int or None
        If set, pre-filter to the top N most variable genes (by variance
        across samples) before training. Dramatically reduces compute time
        for large gene sets (e.g., 5000 out of 25K). Set to None to use
        all genes (default; matches method.txt exactly but is much slower).
    """

    def __init__(
        self,
        tissue: str = "",
        random_state: int = RANDOM_STATE,
        top_n_var_genes: int | None = None,
    ):
        self.tissue = tissue
        self.random_state = random_state
        self.top_n_var_genes = top_n_var_genes
        self._scaler: StandardScaler | None = None
        self._model: ElasticNet | None = None
        self._best_params: dict | None = None
        self._nonzero_genes: pd.Series | None = None  # gene → coefficient
        self._selected_genes: list[str] | None = None  # genes kept after var filter

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @property
    def _max_iter_cv(self) -> int:
        return MAX_ITER_BRAIN if self.tissue.lower() == "brain" else MAX_ITER_DEFAULT

    def _select_genes(self, norm_counts: pd.DataFrame) -> pd.DataFrame:
        """Optionally filter to top N most variable genes.

        If top_n_var_genes is None, returns norm_counts unchanged.
        The selected gene list is stored in self._selected_genes so that
        predict() can apply the same filter to query data.
        """
        if self.top_n_var_genes is None:
            self._selected_genes = norm_counts.index.tolist()
            return norm_counts

        var = norm_counts.var(axis=1)
        top_genes = var.nlargest(self.top_n_var_genes).index.tolist()
        self._selected_genes = top_genes
        n_total = len(norm_counts)
        print(f"  Var-filter: using top {len(top_genes):,} / {n_total:,} genes by variance")
        return norm_counts.loc[top_genes]

    def _prepare_matrices(
        self,
        norm_counts: pd.DataFrame,
        metadata: pd.DataFrame,
        scaler: StandardScaler | None = None,
    ) -> tuple[np.ndarray, np.ndarray, StandardScaler]:
        """Scale expression matrix and extract age vector.

        Parameters
        ----------
        norm_counts : pd.DataFrame
            Normalized counts, genes × samples (already gene-filtered).
        metadata : pd.DataFrame
            Indexed by sample ID; must contain 'age_days'.
        scaler : StandardScaler or None
            If provided, use existing scaler (for test/query sets).

        Returns
        -------
        (X, y, scaler)
        """
        shared = norm_counts.columns.intersection(metadata.index)
        X_raw = norm_counts[shared].T.values  # samples × genes
        y = metadata.loc[shared, "age_days"].astype(float).values

        if scaler is None:
            scaler = StandardScaler()
            X = scaler.fit_transform(X_raw)
        else:
            X = scaler.transform(X_raw)

        return X, y, scaler

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def tune_and_train(
        self,
        norm_counts: pd.DataFrame,
        metadata: pd.DataFrame,
    ) -> dict:
        """Step 1: GridSearchCV with LOO-CV to find optimal alpha and l1_ratio.
        Step 2: Refit final model with max_iter=100,000.

        Parameters
        ----------
        norm_counts : pd.DataFrame
            Normalized counts, genes × samples.
        metadata : pd.DataFrame
            Must contain 'age_days', indexed by sample ID.

        Returns
        -------
        dict
            Best hyperparameters found by GridSearchCV.
        """
        norm_counts = self._select_genes(norm_counts)
        X, y, scaler = self._prepare_matrices(norm_counts, metadata)
        self._scaler = scaler

        gene_names = norm_counts.index.tolist()
        shared = norm_counts.columns.intersection(metadata.index)

        param_grid = {
            "alpha": ALPHA_GRID,
            "l1_ratio": L1_RATIO_GRID,
        }
        base_model = ElasticNet(
            max_iter=self._max_iter_cv,
            random_state=self.random_state,
        )
        gs = GridSearchCV(
            estimator=base_model,
            param_grid=param_grid,
            cv=LeaveOneOut(),
            scoring="neg_mean_absolute_error",
            n_jobs=-1,
            refit=False,
        )
        print(f"  GridSearchCV ({len(X)} samples, "
              f"{len(ALPHA_GRID)}×{len(L1_RATIO_GRID)} grid)...")
        gs.fit(X, y)
        self._best_params = gs.best_params_
        print(f"  Best params: {self._best_params}  "
              f"(MAE={-gs.best_score_:.2f})")

        # Final model with high max_iter for convergence
        self._model = ElasticNet(
            alpha=self._best_params["alpha"],
            l1_ratio=self._best_params["l1_ratio"],
            max_iter=MAX_ITER_FINAL,
            random_state=self.random_state,
        )
        self._model.fit(X, y)

        # Extract non-zero coefficient genes
        coefs = pd.Series(self._model.coef_, index=gene_names)
        self._nonzero_genes = coefs[coefs != 0].sort_values(key=abs, ascending=False)
        print(f"  Non-zero genes: {len(self._nonzero_genes)}")

        return self._best_params

    def loso_cv(
        self,
        norm_counts: pd.DataFrame,
        metadata: pd.DataFrame,
        alpha: float | None = None,
        l1_ratio: float | None = None,
    ) -> pd.DataFrame:
        """Leave-One-Sample-Out CV predictions using fixed hyperparameters.

        If alpha / l1_ratio are not given, uses self._best_params.
        Fits a fresh scaler + model for each fold.

        Parameters
        ----------
        norm_counts : pd.DataFrame
            Normalized counts, genes × samples.
        metadata : pd.DataFrame
            Must contain 'age_days'.
        alpha : float, optional
        l1_ratio : float, optional

        Returns
        -------
        pd.DataFrame
            Columns: age_days, predicted_age — indexed by sample_id.
        """
        if alpha is None:
            if self._best_params is None:
                raise RuntimeError("Run tune_and_train() first or pass alpha/l1_ratio.")
            alpha = self._best_params["alpha"]
            l1_ratio = self._best_params["l1_ratio"]

        # Apply same gene filter used during training
        if self._selected_genes is not None:
            norm_counts = norm_counts.loc[
                norm_counts.index.intersection(self._selected_genes)
            ]

        shared = norm_counts.columns.intersection(metadata.index)
        X_raw = norm_counts[shared].T.values  # samples × genes
        y = metadata.loc[shared, "age_days"].astype(float).values
        sample_ids = shared.tolist()

        predictions = np.empty(len(sample_ids))
        loo = LeaveOneOut()
        for train_idx, test_idx in loo.split(X_raw):
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_raw[train_idx])
            X_test = scaler.transform(X_raw[test_idx])

            model = ElasticNet(
                alpha=alpha,
                l1_ratio=l1_ratio,
                max_iter=MAX_ITER_FINAL,
                random_state=self.random_state,
            )
            model.fit(X_train, y[train_idx])
            predictions[test_idx] = model.predict(X_test)

        results = pd.DataFrame(
            {"age_days": y, "predicted_age": predictions},
            index=sample_ids,
        )
        results.index.name = "sample_id"
        return results

    def predict(self, norm_counts: pd.DataFrame) -> pd.Series:
        """Predict ages for new samples using the trained final model.

        Parameters
        ----------
        norm_counts : pd.DataFrame
            Normalized counts, genes × samples.

        Returns
        -------
        pd.Series
            Predicted ages, indexed by sample ID.
        """
        if self._model is None or self._scaler is None:
            raise RuntimeError("Model not trained. Run tune_and_train() first.")

        # Apply same gene filter used during training
        if self._selected_genes is not None:
            norm_counts = norm_counts.loc[
                norm_counts.index.intersection(self._selected_genes)
            ]

        X = self._scaler.transform(norm_counts.T.values)
        preds = self._model.predict(X)
        return pd.Series(preds, index=norm_counts.columns, name="predicted_age")

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, out_dir: Path, prefix: str = "EN") -> None:
        """Save best params, non-zero genes, and model info."""
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        if self._best_params is not None:
            with open(out_dir / f"{prefix}_best_params.json", "w") as f:
                json.dump(self._best_params, f, indent=2)

        if self._nonzero_genes is not None:
            self._nonzero_genes.to_csv(
                out_dir / f"{prefix}_non_zero_genes.csv", header=True
            )
        print(f"Saved EN outputs to {out_dir}")

    @property
    def nonzero_genes(self) -> pd.Series | None:
        return self._nonzero_genes

    @property
    def best_params(self) -> dict | None:
        return self._best_params
