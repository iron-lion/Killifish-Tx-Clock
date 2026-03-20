"""
Principal Component Regression (PCR) aging clock implementation.

Normalization : StandardScaler (z-scaling).
Training      : Pipeline(StandardScaler → PCA → LinearRegression).
Validation    : LOSO-CV via cross_val_predict for each n_components in range.
Optimal       : n_components that maximises R^2 / minimises MAE.
random_state  : 1. 
"""

import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import pearsonr

RANDOM_STATE = 1


class PCRClock:
    """Principal Component Regression aging clock.

    Parameters
    ----------
    n_components_range : list of int
        PCA component counts to test in LOSO-CV. Default: [5, 10, 15, 20].
    random_state : int
        Random seed for PCA reproducibility.
    top_n_var_genes : int or None
        If set, pre-filter to the top N most variable genes before training.
        None = use all genes (matches method.txt, but slower for large gene sets).
    """

    def __init__(
        self,
        n_components_range: list[int] = None,
        random_state: int = RANDOM_STATE,
        top_n_var_genes: int | None = None,
    ):
        self.n_components_range = n_components_range or [5, 10, 15, 20]
        self.random_state = random_state
        self.top_n_var_genes = top_n_var_genes
        self._pipeline: Pipeline | None = None      # fitted final model
        self._optimal_n: int | None = None          # best n_components from CV
        self._selected_genes: list[str] | None = None
        self._cv_metrics: pd.DataFrame | None = None  # MAE/R² per n_components

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _select_genes(self, norm_counts: pd.DataFrame) -> pd.DataFrame:
        """Optionally filter to top N most variable genes (same logic as ElasticNetClock)."""
        if self.top_n_var_genes is None:
            self._selected_genes = norm_counts.index.tolist()
            return norm_counts
        var = norm_counts.var(axis=1)
        top_genes = var.nlargest(self.top_n_var_genes).index.tolist()
        self._selected_genes = top_genes
        print(f"  Var-filter: using top {len(top_genes):,} / {len(norm_counts):,} genes")
        return norm_counts.loc[top_genes]

    def _make_pipeline(self, n_components: int) -> Pipeline:
        return Pipeline([
            ("scaler", StandardScaler()),
            ("pca", PCA(n_components=n_components, random_state=self.random_state)),
            ("lr", LinearRegression()),
        ])

    def _get_X_y(self, norm_counts: pd.DataFrame, metadata: pd.DataFrame):
        shared = norm_counts.columns.intersection(metadata.index)
        X = norm_counts[shared].T.values  # samples × genes
        y = metadata.loc[shared, "age_days"].astype(float).values
        sample_ids = shared.tolist()
        return X, y, sample_ids

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def loso_cv(
        self,
        norm_counts: pd.DataFrame,
        metadata: pd.DataFrame,
    ) -> pd.DataFrame:
        """LOSO-CV for each n_components; select and store optimal.

        Parameters
        ----------
        norm_counts : pd.DataFrame
            Normalized counts, genes × samples.
        metadata : pd.DataFrame
            Must contain 'age_days', indexed by sample ID.

        Returns
        -------
        pd.DataFrame
            Columns: age_days, pred_ncomp_5, pred_ncomp_10, ...
            Indexed by sample_id. Also stores CV metrics in self._cv_metrics.
        """
        norm_counts = self._select_genes(norm_counts)
        X, y, sample_ids = self._get_X_y(norm_counts, metadata)

        records = {"age_days": y}
        metrics = []

        for n in self.n_components_range:
            pipe = self._make_pipeline(n)
            preds = cross_val_predict(pipe, X, y, cv=LeaveOneOut())
            records[f"pred_ncomp_{n}"] = preds
            mae = mean_absolute_error(y, preds)
            r2 = r2_score(y, preds)
            r, _ = pearsonr(y, preds)
            metrics.append({"n_components": n, "MAE": mae, "R2": r2, "pearson_r": r})
            print(f"  n_components={n:2d}  MAE={mae:.2f}  R²={r2:.3f}  r={r:.3f}")

        self._cv_metrics = pd.DataFrame(metrics).set_index("n_components")

        # Optimal: highest R² (ties broken by lowest MAE)
        self._optimal_n = int(
            self._cv_metrics.sort_values(["R2", "MAE"], ascending=[False, True]).index[0]
        )
        print(f"  Optimal n_components: {self._optimal_n}")

        return pd.DataFrame(records, index=sample_ids).rename_axis("sample_id")

    def fit(
        self,
        norm_counts: pd.DataFrame,
        metadata: pd.DataFrame,
        n_components: int | None = None,
    ) -> None:
        """Fit final PCR model on all data.

        Parameters
        ----------
        norm_counts : pd.DataFrame
            Normalized counts, genes × samples.
        metadata : pd.DataFrame
            Must contain 'age_days'.
        n_components : int, optional
            Use this value; defaults to self._optimal_n (set by loso_cv).
        """
        if n_components is None:
            if self._optimal_n is None:
                raise RuntimeError("Run loso_cv() first or pass n_components explicitly.")
            n_components = self._optimal_n

        # Apply same gene filter
        if self._selected_genes is not None:
            norm_counts = norm_counts.loc[
                norm_counts.index.intersection(self._selected_genes)
            ]

        X, y, _ = self._get_X_y(norm_counts, metadata)
        self._pipeline = self._make_pipeline(n_components)
        self._pipeline.fit(X, y)
        print(f"  Fitted final PCR model  n_components={n_components}")

    def predict(self, norm_counts: pd.DataFrame) -> pd.Series:
        """Predict ages for new samples using the fitted final model.

        Parameters
        ----------
        norm_counts : pd.DataFrame
            Normalized counts, genes × samples.

        Returns
        -------
        pd.Series
            Predicted ages indexed by sample ID.
        """
        if self._pipeline is None:
            raise RuntimeError("Model not fitted. Run fit() first.")

        if self._selected_genes is not None:
            norm_counts = norm_counts.loc[
                norm_counts.index.intersection(self._selected_genes)
            ]

        preds = self._pipeline.predict(norm_counts.T.values)
        return pd.Series(preds, index=norm_counts.columns, name="predicted_age")

    def get_feature_importance(self) -> pd.Series:
        """Composite feature importance per gene (equation 7 from paper).

        FI_g = sum_j ( loading_gj × coeff_j )

        where loading_gj = PCA component j loading for gene g,
              coeff_j    = LinearRegression coefficient for PC j.

        Returns
        -------
        pd.Series
            FI scores indexed by gene name, sorted by absolute value descending.
        """
        if self._pipeline is None:
            raise RuntimeError("Run fit() first.")

        pca: PCA = self._pipeline.named_steps["pca"]
        lr: LinearRegression = self._pipeline.named_steps["lr"]

        # components_ shape: (n_components, n_genes); coef_ shape: (n_components,)
        # FI_g = loadings (n_genes, n_components) @ coef (n_components,)
        fi = pca.components_.T @ lr.coef_  # shape: (n_genes,)

        genes = self._selected_genes or []
        fi_series = pd.Series(fi, index=genes, name="feature_importance")
        return fi_series.reindex(fi_series.abs().sort_values(ascending=False).index)

    def save_loadings(self, out_dir: Path, prefix: str = "PCR", top_n: int = 50) -> None:
        """Save per-component PCA loadings and feature importance scores.

        Parameters
        ----------
        out_dir : Path
            Output directory.
        prefix : str
            Filename prefix.
        top_n : int
            Number of top genes (by |loading|) to save per component.
        """
        if self._pipeline is None:
            raise RuntimeError("Run fit() first.")

        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        pca: PCA = self._pipeline.named_steps["pca"]
        genes = self._selected_genes or []
        n_comp = pca.n_components_

        # Full loadings matrix (genes × components)
        loadings = pd.DataFrame(
            pca.components_.T,
            index=genes,
            columns=[f"PC{j+1}" for j in range(n_comp)],
        )
        loadings.to_csv(out_dir / f"{prefix}_loadings.tsv", sep="\t")

        # Top genes per component
        for j in range(n_comp):
            col = f"PC{j+1}"
            top = loadings[col].abs().nlargest(top_n).index
            top_df = loadings.loc[top, [col]].sort_values(col, ascending=False)
            top_df.to_csv(out_dir / f"{prefix}_top_genes_{col}.csv")

        # Feature importance
        fi = self.get_feature_importance()
        fi.to_csv(out_dir / f"{prefix}_feature_importance.csv", header=True)

        # CV metrics
        if self._cv_metrics is not None:
            self._cv_metrics.to_csv(out_dir / f"{prefix}_cv_metrics.csv")

        print(f"Saved PCR outputs to {out_dir}")

    @property
    def optimal_n_components(self) -> int | None:
        return self._optimal_n

    @property
    def cv_metrics(self) -> pd.DataFrame | None:
        return self._cv_metrics
