"""
run_clocks.py — Apply Atlas aging clocks to a user-provided count matrix.

Input
-----
A genes × samples CSV (or TSV) where sample columns follow:

    TISSUE_repN        e.g.  Liver_rep1,  SpinalCord_rep2,  Muscle_rep3

Gene rows may use Ensembl ENSNFUG IDs (auto-mapped via GeneMapper) or Atlas
gene names (pass --gene-id-type atlas to skip mapping).

Optional --metadata CSV adds per-sample info:
    sample_id, tissue [, age_days, ...]
age_days is only required for EN LOSO-CV on query; pure prediction works
without it.

Clocks
------
  bayesage2  — build reference on Atlas; predict query at M=25..200
  pcr        — LOSO-CV on Atlas to select optimal n_components;
               fit final model; predict query
  en         — GridSearchCV + LOO-CV on Atlas; predict query

Prerequisite
------------
Run once before PCR / EN:
    python src/normalize_reference.py

Usage examples
--------------
    # Minimal (auto-parse tissue from column names, all three clocks):
    python run_clocks.py --counts my_counts.csv

    # With metadata, select tissues, skip batch correction:
    python run_clocks.py --counts my_counts.csv \\
        --metadata my_meta.csv \\
        --tissues Liver Muscle \\
        --clocks bayesage2 pcr \\
        --no-batch-correct

    # Already have Atlas gene names, custom output directory:
    python run_clocks.py --counts my_counts.csv \\
        --gene-id-type atlas \\
        --out-dir results/

Outputs (under --out-dir)
-------------------------
  bayesage2/{tissue}_BayesAge2_predictions.csv
  bayesage2/{tissue}_BayesAge2_feature_importance.csv
  bayesage2/references/{tissue}_reference.tsv
  pcr/{tissue}_PCR_predictions.csv
  pcr/{tissue}_PCR_cv_metrics.csv
  pcr/{tissue}_PCR_feature_importance.csv
  pcr/{tissue}_PCR_loadings.tsv
  elastic_net/{tissue}_EN_predictions.csv
  elastic_net/{tissue}_EN_feature_importance.csv
"""

import argparse
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from bayesage2 import BayesAge2Clock
from data_loader import DataLoader
from elastic_net import ElasticNetClock
from gene_mapping import GeneMapper
from pcr import PCRClock
from preprocessing import Preprocessor

# ── Defaults ──────────────────────────────────────────────────────────────────
_REPO_ROOT       = Path(__file__).resolve().parent
_NORM_DIR        = _REPO_ROOT / "outputs" / "normalized"
_ATLAS_NORM_FILE = _NORM_DIR / "Atlas_DESeq2_normalized.csv"
_DEFAULT_OUT     = _REPO_ROOT / "outputs"

DEFAULT_M_VALUES   = list(range(25, 205, 5))
DEFAULT_N_COMP     = [5, 10, 15, 20]


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_counts(path: Path) -> pd.DataFrame:
    """Read a genes × samples CSV or TSV; coerce all values to float."""
    sep = "\t" if path.suffix in (".tsv", ".txt") else ","
    df = pd.read_csv(path, index_col=0, sep=sep)
    return df.apply(pd.to_numeric, errors="coerce").dropna(how="all")


def parse_sample_cols(columns) -> dict[str, dict]:
    """Parse tissue and replicate number from column names.

    Accepted formats (case-insensitive _rep / _Rep / _REP delimiter):
        TISSUE_repN        Liver_rep1, SpinalCord_Rep2
        TISSUE_N           Liver_1, Muscle_2  (bare number suffix)

    Returns
    -------
    dict  {col_name: {"tissue": str, "rep": int}}
    Columns that do not match either pattern are silently skipped.
    """
    result = {}
    for col in columns:
        m = re.match(r"^(.+?)_rep(\d+)$", col, re.IGNORECASE)
        if m:
            result[col] = {"tissue": m.group(1), "rep": int(m.group(2))}
            continue
        m2 = re.match(r"^(.+?)_(\d+)$", col)
        if m2:
            result[col] = {"tissue": m2.group(1), "rep": int(m2.group(2))}
    return result


def detect_gene_id_type(index) -> str:
    """Return 'ensembl' if index looks like ENSNFUG IDs, else 'atlas'."""
    sample = [g for g in list(index[:30]) if isinstance(g, str)]
    if any(g.startswith("ENSNFUG") for g in sample):
        return "ensembl"
    return "atlas"


# ── Batch correction ──────────────────────────────────────────────────────────

def batch_correct(
    query_counts: pd.DataFrame,
    atlas_raw: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ComBat-seq batch correction applied jointly to Atlas and query.

    Both datasets are corrected in a single pass so that training (Atlas)
    and test (query) data share the same corrected feature space.
    Atlas is used as the reference batch, preserving its overall scale.

    Parameters
    ----------
    query_counts : genes × samples  (float, any scale)
    atlas_raw    : genes × samples  (raw integer Atlas counts)

    Returns
    -------
    (atlas_corrected, query_corrected) : tuple of DataFrames
        Both restricted to the shared gene set, corrected floats.
    """
    from inmoose.pycombat import pycombat_seq

    shared    = atlas_raw.index.intersection(query_counts.index)
    atlas_sub = atlas_raw.loc[shared]
    query_int = query_counts.loc[shared].round().astype(int)

    combined = pd.concat([atlas_sub, query_int], axis=1)
    batch    = ["Atlas"] * atlas_sub.shape[1] + ["Query"] * query_int.shape[1]

    corrected    = pycombat_seq(combined, batch=batch, ref_batch="Atlas")
    corrected_df = pd.DataFrame(corrected, index=shared, columns=combined.columns)

    return (
        corrected_df[atlas_sub.columns].astype(float),
        corrected_df[query_int.columns].astype(float),
    )


# ── Clock runners ─────────────────────────────────────────────────────────────

def run_bayesage2(
    atlas_raw: pd.DataFrame,
    atlas_meta: pd.DataFrame,
    query_counts: pd.DataFrame,
    out_dir: Path,
    m_values: list[int],
    tissue: str,
) -> None:
    """Build reference on Atlas; predict query at each M in m_values."""
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # BayesAge2 expects raw integer counts for the Poisson model
    query_int = query_counts.round().astype(int)

    clock = BayesAge2Clock(lowess_top_n=250)
    ref = clock.build_reference(
        atlas_raw,
        atlas_meta,
        save_path=out_dir.parent / "references" / f"{tissue}_reference.tsv",
    )

    # Predict for each M; collect into a wide DataFrame
    shared = atlas_raw.index.intersection(query_int.index)
    query_int = query_int.loc[shared]

    records: dict[str, dict] = {sid: {} for sid in query_int.columns}
    for m in m_values:
        preds = clock.predict(query_int, n_genes=m, reference=ref)
        for sid, tage in preds.items():
            records[sid][f"tAge_M{m}"] = tage

    result = pd.DataFrame(records).T
    result.index.name = "sample_id"
    result.to_csv(out_dir / f"{tissue}_BayesAge2_predictions.csv")

    # Feature importance: genes ranked by |Spearman r|
    age_cols = [c for c in ref.columns if c.startswith("age_")]
    fi = ref[["spearman_r"]].dropna()
    fi = fi[ref.loc[fi.index, age_cols].notna().all(axis=1)].copy()
    fi["abs_spearman_r"] = fi["spearman_r"].abs()
    fi = fi.sort_values("abs_spearman_r", ascending=False)
    fi["rank"] = range(1, len(fi) + 1)
    fi.index.name = "gene"
    fi.to_csv(out_dir / f"{tissue}_BayesAge2_feature_importance.csv")

    print(f"  [BayesAge2] {tissue}  n_query={len(result)}  "
          f"top_gene={fi.index[0]}  {time.time()-t0:.0f}s")


def run_pcr(
    atlas_norm: pd.DataFrame,
    atlas_meta: pd.DataFrame,
    query_counts: pd.DataFrame,
    out_dir: Path,
    n_components: list[int],
    top_n_var: int | None,
    tissue: str,
) -> None:
    """LOSO-CV on Atlas to select optimal n_components; fit; predict query."""
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    shared = atlas_norm.index.intersection(query_counts.index)

    clock = PCRClock(
        n_components_range=n_components,
        top_n_var_genes=top_n_var,
    )
    clock.loso_cv(atlas_norm.loc[shared], atlas_meta)
    clock.fit(atlas_norm.loc[shared], atlas_meta)

    preds = clock.predict(query_counts)   # PCRClock.predict filters genes internally
    result = preds.to_frame(name=f"predicted_age_n{clock.optimal_n_components}")
    result.index.name = "sample_id"
    result.to_csv(out_dir / f"{tissue}_PCR_predictions.csv")

    # Save CV metrics, loadings, feature importance
    clock.save_loadings(out_dir, prefix=tissue)

    print(f"  [PCR]       {tissue}  n_query={len(result)}  "
          f"optimal_n={clock.optimal_n_components}  {time.time()-t0:.0f}s")


def run_en(
    atlas_norm: pd.DataFrame,
    atlas_meta: pd.DataFrame,
    query_counts: pd.DataFrame,
    out_dir: Path,
    top_n_var: int | None,
    tissue: str,
) -> None:
    """GridSearchCV + LOO-CV on Atlas; predict query samples."""
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    shared = atlas_norm.index.intersection(query_counts.index)

    clock = ElasticNetClock(tissue=tissue, top_n_var_genes=top_n_var)
    clock.tune_and_train(atlas_norm.loc[shared], atlas_meta)

    preds = clock.predict(query_counts.loc[shared])
    result = preds.to_frame(name="predicted_age")
    result.index.name = "sample_id"
    result.to_csv(out_dir / f"{tissue}_EN_predictions.csv")

    if clock.nonzero_genes is not None:
        fi = clock.nonzero_genes.to_frame(name="coefficient")
        fi["abs_coefficient"] = fi["coefficient"].abs()
        fi = fi.sort_values("abs_coefficient", ascending=False)
        fi["rank"] = range(1, len(fi) + 1)
        fi.index.name = "gene"
        fi.to_csv(out_dir / f"{tissue}_EN_feature_importance.csv")
        clock.save(out_dir, prefix=tissue)

    print(f"  [EN]        {tissue}  n_query={len(result)}  "
          f"n_genes={len(clock.nonzero_genes) if clock.nonzero_genes is not None else 0}  "
          f"{time.time()-t0:.0f}s")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Apply Atlas aging clocks to a query count matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--counts", required=True,
        help="genes × samples CSV/TSV. Column names: TISSUE_repN  (e.g. Liver_rep1)",
    )
    parser.add_argument(
        "--metadata", default=None,
        help="Optional sample metadata CSV. Required columns: sample_id, tissue. "
             "Optional: age_days (used for EN LOSO-CV on query).",
    )
    parser.add_argument(
        "--tissues", nargs="+", default=None,
        help="Tissues to run. Must match Atlas tissue labels "
             "(Brain, Liver, Muscle, Fat, Gut, Heart, Kidney, Skin, Bone, "
             "Eye, Spleen, Gonad, SpinalCord). "
             "Default: all tissues detected in column names.",
    )
    parser.add_argument(
        "--clocks", nargs="+", default=["bayesage2", "pcr", "en"],
        choices=["bayesage2", "pcr", "en"],
        help="Clocks to run.",
    )
    parser.add_argument(
        "--gene-id-type", default="auto",
        choices=["auto", "ensembl", "atlas"],
        help="Gene ID type in the count matrix. "
             "'ensembl': ENSNFUG IDs mapped via GeneMapper; "
             "'atlas': Atlas gene names used as-is; "
             "'auto': detected from the index.",
    )
    parser.add_argument(
        "--no-batch-correct", action="store_true",
        help="Skip ComBat-seq batch correction. "
             "Use when query counts are already on the Atlas scale.",
    )
    parser.add_argument(
        "--m-values", nargs="+", type=int, default=DEFAULT_M_VALUES,
        help="BayesAge2: gene-set sizes M to evaluate.",
    )
    parser.add_argument(
        "--n-components", nargs="+", type=int, default=DEFAULT_N_COMP,
        help="PCR: n_components values tested via LOSO-CV on Atlas.",
    )
    parser.add_argument(
        "--top-n-var", type=int, default=None,
        help="Pre-filter to top N most variable genes before PCR/EN training. "
             "Default: use all shared genes.",
    )
    parser.add_argument(
        "--out-dir", default=str(_DEFAULT_OUT),
        help="Output base directory.",
    )
    args = parser.parse_args()

    out_base = Path(args.out_dir)
    clocks   = [c.lower() for c in args.clocks]

    # ── Load query counts ─────────────────────────────────────────────────────
    counts_path = Path(args.counts)
    print(f"Loading query counts: {counts_path}")
    query_all = load_counts(counts_path)
    print(f"  {query_all.shape[0]:,} genes × {query_all.shape[1]} samples")

    # ── Gene ID mapping ───────────────────────────────────────────────────────
    gene_id_type = args.gene_id_type
    if gene_id_type == "auto":
        gene_id_type = detect_gene_id_type(query_all.index)
        print(f"  Gene ID type detected: {gene_id_type}")

    if gene_id_type == "ensembl":
        print("  Mapping ENSNFUG → Atlas gene names via GeneMapper...")
        query_all = GeneMapper().convert(query_all)
        print(f"  After mapping: {query_all.shape[0]:,} genes retained")

    # ── Build sample metadata ─────────────────────────────────────────────────
    if args.metadata:
        meta_all = pd.read_csv(args.metadata).set_index("sample_id")
        shared   = query_all.columns.intersection(meta_all.index)
        query_all = query_all[shared]
        meta_all  = meta_all.loc[shared]
    else:
        parsed = parse_sample_cols(query_all.columns)
        if not parsed:
            print(
                "ERROR: Could not parse tissue from column names.\n"
                "  Expected format: TISSUE_repN  (e.g. Liver_rep1)\n"
                "  Or supply --metadata with columns [sample_id, tissue]."
            )
            sys.exit(1)
        meta_all  = pd.DataFrame(
            [{"sample_id": col, **info} for col, info in parsed.items()]
        ).set_index("sample_id")
        query_all = query_all[list(parsed.keys())]
        n_skip = query_all.shape[1] - len(parsed)
        if n_skip:
            print(f"  Warning: {n_skip} column(s) skipped (unrecognised format)")

    tissues = args.tissues or sorted(meta_all["tissue"].unique())
    print(f"  Tissues: {tissues}")

    # ── Load Atlas data ───────────────────────────────────────────────────────
    print("\nLoading Atlas data...")
    dl = DataLoader()
    atlas_raw_all, atlas_meta_all = dl.load_atlas()
    pp = Preprocessor(min_count=1)
    atlas_raw_all = pp.filter_genes(atlas_raw_all)

    # DESeq2-normalized Atlas is only needed when batch correction is skipped
    # (batch-corrected raw counts are used directly for PCR/EN otherwise).
    atlas_norm_all = None
    if args.no_batch_correct and ("pcr" in clocks or "en" in clocks):
        if not _ATLAS_NORM_FILE.exists():
            print(f"\nERROR: DESeq2-normalized Atlas not found at {_ATLAS_NORM_FILE}")
            print("  Run first:  python src/normalize_reference.py")
            sys.exit(1)
        print("Loading Atlas DESeq2-normalized counts...")
        atlas_norm_all = pd.read_csv(_ATLAS_NORM_FILE, index_col=0)

    # ── Loop over tissues ─────────────────────────────────────────────────────
    total_t0 = time.time()
    for tissue in tissues:
        print(f"\n{'─'*55}")
        print(f"  Tissue: {tissue}")
        print(f"{'─'*55}")

        # Query subset
        t_mask = meta_all["tissue"].str.lower() == tissue.lower()
        if not t_mask.any():
            print(f"  No query samples found, skipping.")
            continue
        q_meta   = meta_all.loc[t_mask].copy()
        q_counts = query_all[q_meta.index]
        print(f"  Query samples: {len(q_meta)}")

        # Atlas subset
        atlas_raw_sub, atlas_meta_sub = pp.stratify(
            atlas_raw_all, atlas_meta_all, tissue=tissue, sex=None
        )
        if len(atlas_meta_sub) < 8:
            print(f"  Too few Atlas samples ({len(atlas_meta_sub)}), skipping.")
            print(f"  Valid Atlas tissues: {sorted(atlas_meta_all['tissue'].unique())}")
            continue
        print(f"  Atlas samples: {len(atlas_meta_sub)}")

        # ── Batch correction ─────────────────────────────────────────────────
        if args.no_batch_correct:
            # Use raw Atlas for BayesAge2; DESeq2-normalized for PCR/EN.
            atlas_corrected = atlas_raw_sub.astype(float)
            q_corrected     = q_counts
            print("  Batch correction: skipped")
        else:
            # Joint correction: Atlas and query corrected in the same pass.
            # The corrected Atlas is used for all clock training so that
            # training and test data share the same corrected feature space.
            print("  Running ComBat-seq batch correction (Atlas + query)...")
            t_bc = time.time()
            atlas_corrected, q_corrected = batch_correct(q_counts, atlas_raw_sub)
            print(f"  Done in {time.time()-t_bc:.0f}s  "
                  f"({len(q_corrected):,} shared genes)")

        # ── Run clocks ───────────────────────────────────────────────────────
        if "bayesage2" in clocks:
            run_bayesage2(
                atlas_corrected, atlas_meta_sub,
                q_corrected,
                out_base / "bayesage2", args.m_values, tissue,
            )

        if "pcr" in clocks:
            atlas_for_pcr = (
                pp.stratify(atlas_norm_all, atlas_meta_all, tissue=tissue, sex=None)[0]
                if args.no_batch_correct else atlas_corrected
            )
            run_pcr(
                atlas_for_pcr, atlas_meta_sub,
                q_corrected,
                out_base / "pcr", args.n_components, args.top_n_var, tissue,
            )

        if "en" in clocks:
            atlas_for_en = (
                pp.stratify(atlas_norm_all, atlas_meta_all, tissue=tissue, sex=None)[0]
                if args.no_batch_correct else atlas_corrected
            )
            run_en(
                atlas_for_en, atlas_meta_sub,
                q_corrected,
                out_base / "elastic_net", args.top_n_var, tissue,
            )

    elapsed = time.time() - total_t0
    print(f"\n{'─'*55}")
    print(f"Done in {elapsed/60:.1f} min  |  outputs in: {out_base}")


if __name__ == "__main__":
    main()
