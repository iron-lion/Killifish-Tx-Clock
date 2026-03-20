"""
run_query_clocks.py — Apply aging clocks to query dataset.

Transfer calibration strategy (method.txt):
  BayesAge2 / PCR : 2 young + 2 old untreated controls from query added to Atlas training.
  EN              : combined Atlas + query LOSO-CV (all query samples used).

"""

import sys
import time
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, './src')

from data_loader import DataLoader
from preprocessing import Preprocessor
from calibration import CalibrationManager, QueryCountExtractor, YOUNG_AGE_DAYS, OLD_AGE_DAYS

OUT_BASE = Path(__file__).resolve().parents[0] / "outputs"
NORM_DIR = OUT_BASE / "normalized"

# Tissues in query data (must match Atlas tissue labels)
TISSUES = [] #NOTE

# ── PARAMETERS ──────────────────────────────────────────────────────────────
YOUNG_AGE = YOUNG_AGE_DAYS   # days — set in calibration.py
OLD_AGE   = OLD_AGE_DAYS     # days — set in calibration.py

# Gene-set sizes for BayesAge2 (method.txt: 5..200 step 5 for external datasets)
M_VALUES = list(range(25, 205, 5))
# PCR n_components range
N_COMPONENTS = [5, 10, 15, 20]
# Optional gene pre-filter for EN/PCR speed (None = use all shared genes)
TOP_N_VAR = None#5000
# ──────────────────────────────────────────────────────────────────────────


def prefix(tissue: str) -> str:
    return f"{tissue}_sexcombined"   # sex-combined (sex unknown for query data)


def run_bayesage2_query(atlas_raw, atlas_meta, query_counts, query_meta, tissue, out_dir, force=False):
    tag = prefix(tissue)
    pred_path = out_dir / f"{tag}_BayesAge2_query.csv"
    fi_path   = out_dir / f"{tag}_BayesAge2_feature_importance.csv"
    ref_path  = out_dir.parent / "references" / f"{tag}_query_reference.tsv"

    if not force and pred_path.exists():
        print(f"  [BayesAge2] {tag} — already done, skipping.")
        return

    t0 = time.time()
    mgr = CalibrationManager()
    result, gene_importance = mgr.run_bayesage2(
        atlas_raw, atlas_meta,
        query_counts, query_meta,
        m_values=M_VALUES,
        ref_save_path=ref_path,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(pred_path)
    gene_importance.to_csv(fi_path)
    print(f"  [BayesAge2] {tag}  n_test={len(result)}  top_gene={gene_importance.index[0]}  "
          f"{time.time()-t0:.0f}s  → {pred_path.name}")


def run_pcr_query(atlas_norm, atlas_meta, query_counts, query_meta, tissue, out_dir, force=False):
    tag = prefix(tissue)
    pred_path = out_dir / f"{tag}_PCR_query.csv"
    mw_path   = out_dir / f"{tag}_PCR_query_mw_pvals.csv"

    if not force and pred_path.exists():
        print(f"  [PCR]       {tag} — already done, skipping.")
        return

    t0 = time.time()
    mgr = CalibrationManager()
    result, mw, gene_importance = mgr.run_pcr(
        atlas_norm, atlas_meta,
        query_counts, query_meta,
        n_components_range=N_COMPONENTS,
        top_n_var_genes=TOP_N_VAR,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(pred_path)
    pd.Series(mw, name="mw_pval").to_csv(mw_path, header=True)
    for n, fi in gene_importance.items():
        fi.to_csv(out_dir / f"{tag}_PCR_feature_importance_n{n}.csv")
    optimal_n = min(mw, key=lambda n: mw[n]) if mw else N_COMPONENTS[0]
    print(f"  [PCR]       {tag}  n_test={len(result)}  {time.time()-t0:.0f}s  "
          f"optimal_n={optimal_n}  mw_pvals={mw}")


def run_en_query(atlas_norm, atlas_meta, query_counts, query_meta, tissue, out_dir, force=False):
    tag = prefix(tissue)
    pred_path = out_dir / f"{tag}_EN_query_loso.csv"
    fi_path   = out_dir / f"{tag}_EN_feature_importance.csv"

    if not force and pred_path.exists():
        print(f"  [EN]        {tag} — already done, skipping.")
        return

    t0 = time.time()
    mgr = CalibrationManager()
    result, gene_importance = mgr.run_en(
        atlas_norm, atlas_meta,
        query_counts, query_meta,
        tissue=tissue,
        top_n_var_genes=TOP_N_VAR,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(pred_path)
    gene_importance.to_csv(fi_path)
    n_q = (result["source"] == "Query").sum()
    print(f"  [EN]        {tag}  n_atlas={len(atlas_meta)}  n_query={n_q}  "
          f"n_genes_selected={len(gene_importance)}  {time.time()-t0:.0f}s  → {pred_path.name}")


def main():
    print("=" * 60)
    print("Query Clock Runner  (Atlas=train, Query=test)")
    print(f"  YOUNG_AGE_DAYS = {YOUNG_AGE}")
    print(f"  OLD_AGE_DAYS   = {OLD_AGE}")
    print("=" * 60)

    # ── Load Atlas data ──────────────────────────────────────────────
    print("\nLoading Atlas data...")
    dl = DataLoader()
    atlas_raw_all, atlas_meta_all = dl.load_atlas()
    pp = Preprocessor(min_count=1)
    atlas_raw_all = pp.filter_genes(atlas_raw_all)

    print("Loading pre-saved Atlas DESeq2-normalized counts...")
    atlas_norm_all = pd.read_csv(NORM_DIR / "Atlas_DESeq2_normalized.csv", index_col=0)

    # ── Load query data ──────────────────────────────────────────────
    print("Loading query counts via GeneMapper...")
    extractor = QueryCountExtractor()

    # ── Loop over tissues ────────────────────────────────────────────
    total_t0 = time.time()
    for tissue in TISSUES:
        print(f"\n{'='*60}")
        print(f"  Tissue = {tissue}")
        print(f"{'='*60}")

        # Query counts for this tissue
        try:
            query_counts, query_meta = extractor.extract_tissue(
                tissue, young_age_days=YOUNG_AGE, old_age_days=OLD_AGE
            )
        except ValueError as e:
            print(f"  Skipping {tissue}: {e}")
            continue

        print(f"  Query samples: {len(query_meta)} | genes overlap: {len(query_counts)}")
        print(f"  Conditions: {query_meta['condition'].value_counts().to_dict()}")

        # Atlas subset for this tissue (sex-combined)
        atlas_raw_sub, atlas_meta_sub = pp.stratify(atlas_raw_all, atlas_meta_all, tissue=tissue, sex=None)
        atlas_norm_sub, _ = pp.stratify(atlas_norm_all, atlas_meta_all, tissue=tissue, sex=None)

        if len(atlas_meta_sub) < 8:
            print(f"  Too few Atlas samples ({len(atlas_meta_sub)}), skipping.")
            continue

        # ── Combat-seq batch correction (query → Atlas distribution) ─────
        print(f"  Running Combat-seq batch correction...")
        t_bc = time.time()
        query_corrected = extractor.correct_batch(query_counts, atlas_raw_sub)
        print(f"  Batch correction done in {time.time()-t_bc:.0f}s  "
              f"({len(query_corrected)} shared genes)")

        # Output directories
        ba2_dir = OUT_BASE / "bayesage2"
        pcr_dir = OUT_BASE / "pcr"
        en_dir  = OUT_BASE / "elastic_net"

        run_bayesage2_query(atlas_raw_sub, atlas_meta_sub, query_corrected, query_meta, tissue, ba2_dir, force=True)
        run_pcr_query(atlas_norm_sub, atlas_meta_sub, query_corrected, query_meta, tissue, pcr_dir, force=True)
        run_en_query(atlas_norm_sub, atlas_meta_sub, query_corrected, query_meta, tissue, en_dir, force=True)

    elapsed = time.time() - total_t0
    print(f"\n{'='*60}")
    print(f"Done in {elapsed/60:.1f} min")
    print(f"Outputs in: {OUT_BASE}")


if __name__ == "__main__":
    main()
