"""Unit tests for GeneMapper (src/gene_mapping.py).

Synthetic tests: verify convert() behaviour with controlled data.
Integration test: measure coverage against the Atlas count matrix
    (data/GSE308970_Counts_Atlas_allbatches_merged_v3.csv).
"""

import sys
import unittest
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from gene_mapping import GeneMapper, MAPPING_FILE, ATLAS_COUNTS_FILE


# ── helpers ───────────────────────────────────────────────────────────────────

def _make_counts(genes, n_samples=4, seed=0):
    """Return a genes × samples DataFrame with random integer counts."""
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        rng.integers(1, 100, size=(len(genes), n_samples)),
        index=genes,
        columns=[f"s{i}" for i in range(n_samples)],
    )


# ── synthetic unit tests ──────────────────────────────────────────────────────

class TestGeneMapperInit(unittest.TestCase):
    def setUp(self):
        self.mapper = GeneMapper()

    def test_map_non_empty(self):
        self.assertGreater(len(self.mapper._map), 0)

    def test_map_size(self):
        # Expected unique ENSNFUG → Atlas mappings after dedup
        self.assertGreaterEqual(len(self.mapper._map), 12_000)

    def test_map_index_is_ensembl(self):
        self.assertTrue(self.mapper._map.index.str.startswith("ENSNFUG").all())

    def test_no_nan_in_map(self):
        self.assertFalse(self.mapper._map.isna().any())

    def test_no_duplicate_ensembl_ids(self):
        self.assertFalse(self.mapper._map.index.duplicated().any())


class TestGeneMapperConvert(unittest.TestCase):
    def setUp(self):
        # Minimal synthetic mapping: 3 Ensembl IDs → Atlas names
        df = pd.DataFrame({
            "ensembl_gene_id": ["ENSNFUG00000000001", "ENSNFUG00000000002", "ENSNFUG00000000003"],
            "atlas_gene":      ["actb",               "gapdh",              "actb"],  # dup atlas name
        })
        self.mapper = GeneMapper.__new__(GeneMapper)
        df = df.dropna(subset=["ensembl_gene_id", "atlas_gene"])
        df = df.drop_duplicates(subset="ensembl_gene_id", keep="first")
        self.mapper._map = df.set_index("ensembl_gene_id")["atlas_gene"]

    def test_mapped_genes_renamed(self):
        counts = _make_counts(["ENSNFUG00000000001", "ENSNFUG00000000002"])
        out = self.mapper.convert(counts)
        self.assertIn("actb", out.index)
        self.assertIn("gapdh", out.index)

    def test_unmapped_genes_dropped(self):
        counts = _make_counts(["ENSNFUG00000000001", "ENSNFUG99999999999"])
        out = self.mapper.convert(counts)
        self.assertEqual(len(out), 1)
        self.assertIn("actb", out.index)

    def test_shape_preserved(self):
        counts = _make_counts(["ENSNFUG00000000001", "ENSNFUG00000000002"])
        out = self.mapper.convert(counts)
        self.assertEqual(out.shape, (2, 4))

    def test_duplicate_atlas_name_deduplicated(self):
        # ENSNFUG00000000001 and ENSNFUG00000000003 both map to "actb"
        counts = _make_counts([
            "ENSNFUG00000000001",
            "ENSNFUG00000000002",
            "ENSNFUG00000000003",
        ])
        out = self.mapper.convert(counts)
        self.assertEqual(out.index.tolist().count("actb"), 1)

    def test_all_unmapped_returns_empty(self):
        counts = _make_counts(["ENSNFUG99999999999", "ENSNFUG88888888888"])
        out = self.mapper.convert(counts)
        self.assertEqual(len(out), 0)

    def test_values_preserved(self):
        counts = _make_counts(["ENSNFUG00000000001"])
        out = self.mapper.convert(counts)
        np.testing.assert_array_equal(out.values, counts.values)


# ── integration test against Atlas matrix ─────────────────────────────────────

@unittest.skipUnless(
    ATLAS_COUNTS_FILE.exists(),
    f"Atlas count matrix not found at {ATLAS_COUNTS_FILE}; skipping integration test.",
)
class TestGeneMapperAtlasCoverage(unittest.TestCase):
    """Measure reverse coverage: how many Atlas genes can GeneMapper reach?

    Coverage stats (computed 2026-05-26 against GCF_043380555.1 / GSE308970):
        Atlas total      : 25,122 genes
        Covered          : 12,017  (47.8 %)
        LOC genes        : 4,311 / 14,589  (29.5 %)
        Symbol genes     : 7,706 / 10,533  (73.2 %)
    """

    @classmethod
    def setUpClass(cls):
        cls.mapper = GeneMapper()
        # Read only the gene index (nrows=0 loads header + index only)
        atlas = pd.read_csv(ATLAS_COUNTS_FILE, index_col=0, nrows=0)
        cls.atlas_genes = pd.Index(atlas.columns)  # genes are columns after index_col=0

        # Fall back: some versions store genes as row index
        if len(cls.atlas_genes) < 1000:
            atlas_t = pd.read_csv(ATLAS_COUNTS_FILE, index_col=0, nrows=1)
            cls.atlas_genes = pd.Index(
                pd.read_csv(ATLAS_COUNTS_FILE, usecols=[0]).iloc[:, 0].tolist()
            )

        cls.atlas_gene_set = set(cls.atlas_genes)
        cls.mapped_atlas = set(cls.mapper._map.values)

        cls.covered = cls.atlas_gene_set & cls.mapped_atlas
        cls.loc_genes = {g for g in cls.atlas_gene_set if g.startswith("LOC")}
        cls.sym_genes = cls.atlas_gene_set - cls.loc_genes
        cls.covered_loc = cls.loc_genes & cls.mapped_atlas
        cls.covered_sym = cls.sym_genes & cls.mapped_atlas

    def test_atlas_gene_count(self):
        self.assertGreater(len(self.atlas_gene_set), 20_000,
                           "Atlas should have >20k genes")

    def test_mapper_size(self):
        self.assertGreaterEqual(len(self.mapper._map), 12_000,
                                "GeneMapper should carry ≥12,000 ENSNFUG mappings")

    def test_overall_coverage_fraction(self):
        frac = len(self.covered) / len(self.atlas_gene_set)
        self.assertGreater(frac, 0.40, f"Overall coverage {frac:.1%} < 40 %")
        self.assertLess(frac, 0.65, f"Overall coverage {frac:.1%} suspiciously high")

    def test_symbol_coverage_higher_than_loc(self):
        sym_frac = len(self.covered_sym) / len(self.sym_genes) if self.sym_genes else 0
        loc_frac = len(self.covered_loc) / len(self.loc_genes) if self.loc_genes else 0
        self.assertGreater(sym_frac, loc_frac,
                           "Symbol gene coverage should exceed LOC gene coverage")

    def test_symbol_coverage_above_60pct(self):
        sym_frac = len(self.covered_sym) / len(self.sym_genes) if self.sym_genes else 0
        self.assertGreater(sym_frac, 0.60,
                           f"Symbol gene coverage {sym_frac:.1%} < 60 %")

    def test_loc_coverage_above_20pct(self):
        loc_frac = len(self.covered_loc) / len(self.loc_genes) if self.loc_genes else 0
        self.assertGreater(loc_frac, 0.20,
                           f"LOC gene coverage {loc_frac:.1%} < 20 %")

    def test_print_coverage_report(self):
        """Not an assertion — prints a human-readable coverage summary."""
        n_atlas = len(self.atlas_gene_set)
        n_cov = len(self.covered)
        n_loc = len(self.loc_genes)
        n_sym = len(self.sym_genes)
        n_cloc = len(self.covered_loc)
        n_csym = len(self.covered_sym)
        print(
            f"\n── GeneMapper coverage against Atlas ──\n"
            f"  Atlas total genes : {n_atlas:,}\n"
            f"  Covered (overall) : {n_cov:,}  ({n_cov/n_atlas:.1%})\n"
            f"  Symbol genes      : {n_csym:,} / {n_sym:,}  ({n_csym/n_sym:.1%})\n"
            f"  LOC genes         : {n_cloc:,} / {n_loc:,}  ({n_cloc/n_loc:.1%})\n"
            f"  Not covered       : {n_atlas - n_cov:,}\n"
        )


# ── integration test: PRJNA817434 direct Atlas overlap ───────────────────────

PRJNA817434_FILE = (
    Path(__file__).resolve().parents[1]
    / "raw_RNAseq_process/results/PRJNA817434/PRJNA817434_raw_count.csv"
)


@unittest.skipUnless(
    PRJNA817434_FILE.exists() and ATLAS_COUNTS_FILE.exists(),
    "PRJNA817434 count file or Atlas matrix not found; skipping.",
)
class TestPRJNA817434AtlasCoverage(unittest.TestCase):
    """Direct Atlas name overlap for PRJNA817434 (no ENSNFUG IDs — GeneMapper N/A).

    PRJNA817434 is produced by raw_RNAseq_process using GCF_043380555.1 GTF.
    Gene IDs are NCBI-style: named symbols, LOC107XXXXXX, LOC139XXXXXX, tRNA,
    and KEG92_tXX entries.  No ENSNFUG IDs → GeneMapper.convert() does not apply.
    Coverage is measured as direct set intersection with Atlas gene names.

    Coverage stats (computed 2026-05-26):
        PRJNA817434 total        : 36,530 genes
        Direct Atlas match       : 12,249  (33.5 %)
        Symbol genes             : 7,098 / 14,706  (48.3 %)
        LOC107 genes             : 5,103 / 6,203   (82.3 %)
        LOC139 genes             : 0     / 7,031   ( 0.0 %)  ← different annotation
        tRNA / KEG92             : not in Atlas
    """

    @classmethod
    def setUpClass(cls):
        prjna = pd.read_csv(PRJNA817434_FILE, usecols=[0]).iloc[:, 0]
        cls.query_genes = pd.Index(prjna.tolist())
        cls.atlas_set = set(
            pd.read_csv(ATLAS_COUNTS_FILE, usecols=[0]).iloc[:, 0].tolist()
        )

        cls.loc107 = {g for g in cls.query_genes if g.startswith("LOC107")}
        cls.loc139 = {g for g in cls.query_genes if g.startswith("LOC139")}
        cls.trna   = {g for g in cls.query_genes if g.startswith("trn")}
        cls.keg    = {g for g in cls.query_genes if g.startswith("KEG92")}
        cls.sym    = (
            set(cls.query_genes) - cls.loc107 - cls.loc139 - cls.trna - cls.keg
        )

        cls.covered       = set(cls.query_genes) & cls.atlas_set
        cls.covered_sym   = cls.sym   & cls.atlas_set
        cls.covered_loc107 = cls.loc107 & cls.atlas_set
        cls.covered_loc139 = cls.loc139 & cls.atlas_set

    def test_query_gene_count(self):
        self.assertGreater(len(self.query_genes), 30_000)

    def test_loc139_not_in_atlas(self):
        # LOC139XXXXXX come from a newer NCBI annotation; Atlas uses LOC107XXXXXX
        self.assertEqual(len(self.covered_loc139), 0,
                         "LOC139 genes should have 0 % Atlas coverage")

    def test_loc107_coverage_high(self):
        frac = len(self.covered_loc107) / len(self.loc107) if self.loc107 else 0
        self.assertGreater(frac, 0.70,
                           f"LOC107 coverage {frac:.1%} < 70 %")

    def test_symbol_coverage(self):
        frac = len(self.covered_sym) / len(self.sym) if self.sym else 0
        self.assertGreater(frac, 0.40,
                           f"Symbol coverage {frac:.1%} < 40 %")

    def test_overall_coverage(self):
        frac = len(self.covered) / len(self.query_genes)
        self.assertGreater(frac, 0.25, f"Overall coverage {frac:.1%} < 25 %")

    def test_print_coverage_report(self):
        n_total = len(self.query_genes)
        n_cov   = len(self.covered)
        n_sym   = len(self.sym);    n_csym  = len(self.covered_sym)
        n_107   = len(self.loc107); n_c107  = len(self.covered_loc107)
        n_139   = len(self.loc139); n_c139  = len(self.covered_loc139)
        n_trna  = len(self.trna);   n_keg   = len(self.keg)
        print(
            f"\n── PRJNA817434 direct Atlas coverage ──\n"
            f"  Total genes               : {n_total:,}\n"
            f"  Direct Atlas match        : {n_cov:,}  ({n_cov/n_total:.1%})\n"
            f"  Symbol genes              : {n_csym:,} / {n_sym:,}  ({n_csym/n_sym:.1%})\n"
            f"  LOC107 genes              : {n_c107:,} / {n_107:,}  ({n_c107/n_107:.1%})\n"
            f"  LOC139 genes (new annot.) : {n_c139:,} / {n_139:,}  ( 0.0 % — not in Atlas)\n"
            f"  tRNA genes                : {n_trna:,}  (excluded)\n"
            f"  KEG92 entries             : {n_keg:,}  (excluded)\n"
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
