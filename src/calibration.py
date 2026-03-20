"""
Apply Atlas aging clocks to query datasets (Atlas = train, Query = test).

Query data: per-tissue fasting/refeeding experiment.
Count values in query xlsx files are DESeq2-normalized (decimal), NOT raw integers.

Strategy:
  BayesAge2 : build reference on Atlas only; predict all query samples.
  PCR       : fit Pipeline on Atlas only; predict all query samples.
  EN        : tune + train on Atlas only; predict all query samples.

"""

import re
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from scipy.stats import mannwhitneyu
from inmoose.pycombat import pycombat_seq

from bayesage2 import BayesAge2Clock
from elastic_net import ElasticNetClock
from gene_mapping import GeneMapper

REPO_ROOT = Path(__file__).resolve().parent.parent
QUERY_DIR = REPO_ROOT / "query_data"

# ── USER-CONFIGURABLE ──────────────────────────────────────────────────────
YOUNG_AGE_DAYS: int = 56     # ~8 weeks 
OLD_AGE_DAYS: int = 126      # ~18 weeks
CONTROL_CONDITION: str = ""  # NOTE: untreated control
# ──────────────────────────────────────────────────────────────────────────

_TISSUE_MAP = {"FAT": "Fat", "LIVER": "Liver", "MUSCLE": "Muscle"}

# Columns in query xlsx that are not raw/norm counts
_NON_COUNT_COLS = frozenset({
    "ensembl_gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE",
    "pvalue", "padj", "gene_biotype", "GO_id", "GO_term", "log2FCrev",
})


def _normalize_condition(raw: str) -> str:
    """Canonicalize condition name (handles mixed case / abbreviations)."""
    c = raw.upper()
    if "72" in c:
        return "72"
    if "24" in c:
        return "24"
    if "6" in c:
        return "6"
    return c


def _parse_sample_name(col: str) -> dict | None:
    """Parse sample column name to (tissue, age_group, condition, rep).

    Expected format: Tissue_AgeGroup_Condition.Rep_N
    """
    m = re.match(r"^([A-Za-z]+)_(Young|Old)_(.+)\.Rep_(\d+)$", col, re.IGNORECASE)
    if not m:
        return None
    raw_tissue, age_group, cond_raw, rep = m.groups()
    tissue = _TISSUE_MAP.get(raw_tissue.upper(), raw_tissue.capitalize())
    return dict(
        tissue=tissue,
        age_group=age_group.capitalize(),
        condition=_normalize_condition(cond_raw),
        rep=int(rep),
    )


class QueryCountExtractor:
    """Extract per-tissue count matrices from query DE result xlsx files.

    The xlsx files contain per-sample count columns alongside DE statistics.
    This class identifies count columns, deduplicates samples appearing
    in multiple comparison files, converts ENSNFUG gene IDs to Atlas
    gene names via GeneMapper, and returns a count matrix + metadata.

    Notes
    -----
    - Count values are DESeq2-normalized (decimal). Use as-is for EN/PCR.
      Round to int for BayesAge2 (approximation).
    - Genes are intersected across all tissue files to ensure consistency.
    """

    def __init__(self, query_dir: Path = QUERY_DIR, mapper: GeneMapper | None = None):
        self.query_dir = Path(query_dir)
        self._mapper = mapper or GeneMapper()

    def extract_tissue(
        self,
        tissue: str,
        young_age_days: int = YOUNG_AGE_DAYS,
        old_age_days: int = OLD_AGE_DAYS,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Extract count matrix and metadata for one tissue.

        Parameters
        ----------
        tissue : str
            Atlas tissue label.
        young_age_days : int
            Age (days) assigned to Young samples.
        old_age_days : int
            Age (days) assigned to Old samples.

        Returns
        -------
        counts : pd.DataFrame
            Atlas gene names × samples. Values are DESeq2-normalized.
        metadata : pd.DataFrame
            Index = sample_id. Columns: tissue, age_group, condition, rep, age_days.
        """
        tissue_key = tissue.upper()

        # Collect per-file gene sets and sample columns
        file_gene_sets: list[set] = []
        sample_cols_per_file: dict[str, list[str]] = {}  # file_path → sample col names
        file_dfs: dict[str, pd.DataFrame] = {}

        for xlsx_path in sorted(self.query_dir.glob("*.xlsx")):
            if tissue_key not in xlsx_path.name.upper():
                continue

            df = pd.read_excel(xlsx_path)
            if "ensembl_gene_id" in df.columns:
                df = df.set_index("ensembl_gene_id")

            count_cols = [c for c in df.columns if c not in _NON_COUNT_COLS and _parse_sample_name(c) is not None]
            # Filter to this tissue's columns
            count_cols = [c for c in count_cols if (_parse_sample_name(c) or {}).get("tissue", "").upper() == tissue_key]

            if not count_cols:
                continue

            file_dfs[str(xlsx_path)] = df
            sample_cols_per_file[str(xlsx_path)] = count_cols
            file_gene_sets.append(set(df.index))

        if not file_dfs:
            raise ValueError(f"No query xlsx files found for tissue='{tissue}'.")

        # Intersect genes across files for consistency
        common_genes = file_gene_sets[0]
        for gs in file_gene_sets[1:]:
            common_genes &= gs
        common_genes = sorted(common_genes)

        # Collect unique samples (deduplicate across files)
        sample_frames: dict[str, pd.Series] = {}
        sample_info: dict[str, dict] = {}

        for fpath, df in file_dfs.items():
            for col in sample_cols_per_file[fpath]:
                if col in sample_frames:
                    continue  # already collected from another file
                info = _parse_sample_name(col)
                if info is None:
                    continue
                sample_frames[col] = df.loc[common_genes, col].astype(float)
                sample_info[col] = info

        # Build raw count matrix (ENSNFUG × samples)
        counts_ensnfug = pd.DataFrame(sample_frames, index=common_genes)

        # Convert gene IDs to Atlas names
        counts = self._mapper.convert(counts_ensnfug)

        # Build metadata
        rows = []
        for sid, info in sample_info.items():
            age_days = young_age_days if info["age_group"] == "Young" else old_age_days
            rows.append({"sample_id": sid, **info, "age_days": age_days})
        metadata = pd.DataFrame(rows).set_index("sample_id")

        return counts, metadata

    def correct_batch(
        self,
        query_counts: pd.DataFrame,
        atlas_raw: pd.DataFrame,
    ) -> pd.DataFrame:
        """Batch-correct query counts toward the Atlas distribution via Combat-seq.

        Atlas is used as the reference batch so that the corrected query counts
        are on the same count scale as the Atlas raw data.

        Parameters
        ----------
        query_counts : pd.DataFrame
            Query counts, genes × samples.
        atlas_raw : pd.DataFrame
            Raw integer Atlas counts, genes × samples.

        Returns
        -------
        pd.DataFrame
            Batch-corrected query counts (genes × samples, float) on Atlas scale.
        """
        shared_genes = atlas_raw.index.intersection(query_counts.index)
        atlas_sub  = atlas_raw.loc[shared_genes]
        query_int  = query_counts.loc[shared_genes].round().astype(int)

        # Concatenate into one genes × (atlas + query) matrix
        combined = pd.concat([atlas_sub, query_int], axis=1)
        batch    = ["Atlas"] * atlas_sub.shape[1] + ["Query"] * query_int.shape[1]

        corrected = pycombat_seq(combined, batch=batch, ref_batch="Atlas")

        corrected_df = pd.DataFrame(
            corrected,
            index=shared_genes,
            columns=combined.columns,
        )
        return corrected_df[query_int.columns].astype(float)


class CalibrationManager:
    """Apply Atlas aging clocks to query datasets (Atlas = train, Query = test)."""

    # ------------------------------------------------------------------
    # BayesAge2
    # ------------------------------------------------------------------

    def run_bayesage2(
        self,
        atlas_raw: pd.DataFrame,
        atlas_meta: pd.DataFrame,
        query_counts: pd.DataFrame,
        query_meta: pd.DataFrame,
        m_values: list[int] | None = None,
        ref_save_path: Path | None = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Build reference on Atlas only; predict all query samples.

        Query count values are DESeq2-normalized (decimal) — rounded to int
        as a raw-count approximation for the Poisson model.

        Returns
        -------
        result : pd.DataFrame
            tAge_M{m} columns + age_group + condition, indexed by sample_id.
        gene_importance : pd.DataFrame
            All genes with valid LOWESS ranked by |spearman_r|.
            Columns: spearman_r, abs_spearman_r, rank.
        """
        if m_values is None:
            m_values = list(range(5, 205, 5))

        shared_genes = atlas_raw.index.intersection(query_counts.index)

        clock = BayesAge2Clock(lowess_top_n=250)
        ref = clock.build_reference(
            atlas_raw.loc[shared_genes], atlas_meta, save_path=ref_save_path
        )

        test_raw = query_counts.loc[shared_genes].round().astype(int)
        test_ids = test_raw.columns.tolist()

        records = {sid: {} for sid in test_ids}
        for m in m_values:
            preds = clock.predict(test_raw, n_genes=m, reference=ref)
            for sid, tage in preds.items():
                records[sid][f"tAge_M{m}"] = tage

        result = pd.DataFrame(records).T
        result.index.name = "sample_id"
        result["age_group"] = query_meta["age_group"]
        result["condition"] = query_meta["condition"]

        # Feature importance: genes ranked by |Spearman r| with LOWESS fits
        age_cols = [c for c in ref.columns if c.startswith("age_")]
        fi = ref[["spearman_r"]].dropna()
        fi = fi[ref[age_cols].notna().all(axis=1)]   # only genes with LOWESS
        fi = fi.copy()
        fi["abs_spearman_r"] = fi["spearman_r"].abs()
        fi = fi.sort_values("abs_spearman_r", ascending=False)
        fi["rank"] = range(1, len(fi) + 1)
        fi.index.name = "gene"

        return result, fi

    # ------------------------------------------------------------------
    # PCR
    # ------------------------------------------------------------------

    def run_pcr(
        self,
        atlas_norm: pd.DataFrame,
        atlas_meta: pd.DataFrame,
        query_norm: pd.DataFrame,
        query_meta: pd.DataFrame,
        n_components_range: list[int] | None = None,
        top_n_var_genes: int | None = None,
    ) -> tuple[pd.DataFrame, dict, dict]:
        """Fit Pipeline on Atlas only; predict all query samples.

        Returns
        -------
        predictions : pd.DataFrame
            tAge_n{n} columns + age_group + condition, indexed by sample_id.
        mw_pvals : dict
            Mann-Whitney U p-values (Young vs Old) per n_components.
        gene_importance : dict[int, pd.DataFrame]
            Per n_components: DataFrame with columns pca_loading_score,
            regression_coef_contribution, importance_score, rank.
            importance_score = |loadings.T @ reg_coefs| — the absolute
            contribution of each gene to the final age prediction.
        """
        if n_components_range is None:
            n_components_range = [5, 10, 15, 20]

        shared_genes = atlas_norm.index.intersection(query_norm.index)
        train_norm = atlas_norm.loc[shared_genes]

        if top_n_var_genes:
            top_genes = train_norm.var(axis=1).nlargest(top_n_var_genes).index
            shared_genes = top_genes
            train_norm = train_norm.loc[shared_genes]

        gene_names = shared_genes.tolist()
        X_train = train_norm.T.values
        y_train = atlas_meta["age_days"].values

        test_norm = query_norm.loc[shared_genes]
        test_ids = test_norm.columns.tolist()

        records = {sid: {} for sid in test_ids}
        mw_pvals: dict[int, float] = {}
        gene_importance: dict[int, pd.DataFrame] = {}

        for n in n_components_range:
            pipe = Pipeline([
                ("scaler", StandardScaler()),
                ("pca", PCA(n_components=n, random_state=1)),
                ("reg", LinearRegression()),
            ])
            pipe.fit(X_train, y_train)
            preds = pipe.predict(test_norm.T.values)

            for sid, p in zip(test_ids, preds):
                records[sid][f"tAge_n{n}"] = p

            young_p = [records[s][f"tAge_n{n}"] for s in test_ids if query_meta.loc[s, "age_group"] == "Young"]
            old_p   = [records[s][f"tAge_n{n}"] for s in test_ids if query_meta.loc[s, "age_group"] == "Old"]
            if young_p and old_p:
                _, pval = mannwhitneyu(young_p, old_p, alternative="two-sided")
                mw_pvals[n] = pval

            # Feature importance: |loadings.T @ reg_coefs|
            # loadings: (n_components, n_genes); reg_coefs: (n_components,)
            loadings = pipe["pca"].components_          # (n, n_genes)
            reg_coefs = pipe["reg"].coef_               # (n,)
            importance = np.abs(loadings.T @ reg_coefs) # (n_genes,)
            # Per-PC contribution for inspection: |loading| * |coef| summed
            pc_contrib = (np.abs(loadings.T) * np.abs(reg_coefs)).sum(axis=1)

            fi = pd.DataFrame({
                "importance_score": importance,
                "pc_contribution":  pc_contrib,
            }, index=gene_names)
            fi = fi.sort_values("importance_score", ascending=False)
            fi["rank"] = range(1, len(fi) + 1)
            fi.index.name = "gene"
            gene_importance[n] = fi

        result = pd.DataFrame(records).T
        result.index.name = "sample_id"
        result["age_group"] = query_meta["age_group"]
        result["condition"] = query_meta["condition"]
        return result, mw_pvals, gene_importance

    # ------------------------------------------------------------------
    # Elastic Net
    # ------------------------------------------------------------------

    def run_en(
        self,
        atlas_norm: pd.DataFrame,
        atlas_meta: pd.DataFrame,
        query_norm: pd.DataFrame,
        query_meta: pd.DataFrame,
        tissue: str = "",
        top_n_var_genes: int | None = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Tune + train on Atlas only; predict all query samples.

        Returns
        -------
        result : pd.DataFrame
            Columns: age_days, predicted_age, source, age_group, condition.
            source = 'Atlas' for atlas LOSO-CV rows, 'Query' for query predictions.
        gene_importance : pd.DataFrame
            Non-zero EN coefficients ranked by |coefficient|.
            Columns: coefficient, abs_coefficient, rank.
        """
        shared_genes = atlas_norm.index.intersection(query_norm.index)

        en = ElasticNetClock(tissue=tissue, top_n_var_genes=top_n_var_genes)
        en.tune_and_train(atlas_norm.loc[shared_genes], atlas_meta)

        # Atlas: LOSO-CV predictions (for reference/evaluation)
        atlas_cv = en.loso_cv(atlas_norm.loc[shared_genes], atlas_meta)
        atlas_cv["source"] = "Atlas"

        # Query: predict using the model trained on Atlas
        query_preds = en.predict(query_norm.loc[shared_genes])
        query_rows = pd.DataFrame({
            "age_days": query_meta["age_days"],
            "predicted_age": query_preds,
            "source": "Query",
            "age_group": query_meta["age_group"],
            "condition": query_meta["condition"],
        }, index=query_meta.index)

        result = pd.concat([atlas_cv, query_rows])

        # Feature importance: non-zero ElasticNet coefficients
        fi = en.nonzero_genes.to_frame(name="coefficient")
        fi["abs_coefficient"] = fi["coefficient"].abs()
        fi = fi.sort_values("abs_coefficient", ascending=False)
        fi["rank"] = range(1, len(fi) + 1)
        fi.index.name = "gene"

        return result, fi
