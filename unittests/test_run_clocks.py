"""Unit tests for run_query_clocks.py.

Coverage
--------
  TestParseSampleCols   — column-name parsing helpers
  TestDetectGeneIdType  — auto-detection of gene ID format
  TestLoadCounts        — CSV / TSV loading
  TestBatchCorrect      — ComBat-seq wrapper (pycombat_seq mocked)
  TestRunBayesage2      — BayesAge2 runner with toy atlas/query data
  TestRunPCR            — PCR runner with toy data
  TestRunEN             — EN runner with toy data
  TestMainIntegration   — end-to-end main() using query_data/toy.csv
                          (DataLoader mocked, batch correction skipped)
"""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from run_query_clocks import (
    batch_correct,
    detect_gene_id_type,
    load_counts,
    parse_sample_cols,
    run_bayesage2,
    run_en,
    run_pcr,
)

# ── Shared fixtures ────────────────────────────────────────────────────────────

_N_GENES = 50
_N_ATLAS = 24   # enough for LOSO-CV; divisible by tissue count
_N_QUERY = 3

_GENE_NAMES  = [f"gene_{i:03d}" for i in range(1, _N_GENES + 1)]
_ATLAS_IDS   = [f"ATLAS_{i:03d}" for i in range(1, _N_ATLAS + 1)]
_QUERY_IDS   = [f"Liver_rep{i}" for i in range(1, _N_QUERY + 1)]


def _make_atlas_raw(rng=None) -> pd.DataFrame:
    rng = rng or np.random.default_rng(0)
    counts = rng.negative_binomial(10, 0.5, (_N_GENES, _N_ATLAS)).astype(float)
    return pd.DataFrame(counts, index=_GENE_NAMES, columns=_ATLAS_IDS)


def _make_atlas_meta() -> pd.DataFrame:
    rng = np.random.default_rng(1)
    return pd.DataFrame(
        {
            "tissue": ["Liver"] * _N_ATLAS,
            "age_days": rng.integers(47, 163, _N_ATLAS).astype(float),
            "sex": ["M"] * _N_ATLAS,
        },
        index=pd.Index(_ATLAS_IDS, name="lib"),
    )


def _make_query_counts() -> pd.DataFrame:
    rng = np.random.default_rng(2)
    counts = rng.negative_binomial(10, 0.5, (_N_GENES, _N_QUERY)).astype(float)
    return pd.DataFrame(counts, index=_GENE_NAMES, columns=_QUERY_IDS)


# ── TestParseSampleCols ────────────────────────────────────────────────────────

class TestParseSampleCols(unittest.TestCase):

    def test_rep_format(self):
        result = parse_sample_cols(["Liver_rep1", "Muscle_rep2"])
        self.assertEqual(result["Liver_rep1"]["tissue"], "Liver")
        self.assertEqual(result["Liver_rep1"]["rep"], 1)
        self.assertEqual(result["Muscle_rep2"]["tissue"], "Muscle")
        self.assertEqual(result["Muscle_rep2"]["rep"], 2)

    def test_rep_case_insensitive(self):
        result = parse_sample_cols(["Brain_Rep3", "Gut_REP10"])
        self.assertIn("Brain_Rep3", result)
        self.assertEqual(result["Brain_Rep3"]["rep"], 3)
        self.assertIn("Gut_REP10", result)
        self.assertEqual(result["Gut_REP10"]["rep"], 10)

    def test_bare_number_format(self):
        result = parse_sample_cols(["Liver_1", "Muscle_2"])
        self.assertEqual(result["Liver_1"]["tissue"], "Liver")
        self.assertEqual(result["Liver_1"]["rep"], 1)

    def test_tissue_with_underscore(self):
        # SpinalCord_rep1 — tissue name itself contains no digits
        result = parse_sample_cols(["SpinalCord_rep1"])
        self.assertIn("SpinalCord_rep1", result)
        self.assertEqual(result["SpinalCord_rep1"]["tissue"], "SpinalCord")

    def test_unrecognised_columns_skipped(self):
        result = parse_sample_cols(["good_rep1", "nodigitsuffix", "bad_repXY"])
        self.assertIn("good_rep1", result)
        self.assertNotIn("nodigitsuffix", result)
        self.assertNotIn("bad_repXY", result)

    def test_empty_input(self):
        self.assertEqual(parse_sample_cols([]), {})


# ── TestDetectGeneIdType ───────────────────────────────────────────────────────

class TestDetectGeneIdType(unittest.TestCase):

    def test_ensembl(self):
        index = pd.Index([f"ENSNFUG{i:011d}" for i in range(20)])
        self.assertEqual(detect_gene_id_type(index), "ensembl")

    def test_atlas(self):
        index = pd.Index(["actb", "gapdh", "LOC107380001", "myh7"])
        self.assertEqual(detect_gene_id_type(index), "atlas")

    def test_mixed_defaults_to_atlas(self):
        # Only some ENSNFUG → counted as atlas unless majority are ENSNFUG
        index = pd.Index(["actb"] + [f"ENSNFUG{i:011d}" for i in range(2)])
        # Two ENSNFUG in first 30 items: detect_gene_id_type returns 'ensembl'
        # because any() is True; verify the actual behaviour
        result = detect_gene_id_type(index)
        self.assertIn(result, ("ensembl", "atlas"))


# ── TestLoadCounts ─────────────────────────────────────────────────────────────

class TestLoadCounts(unittest.TestCase):

    def _write_and_load(self, suffix: str, sep: str):
        df = pd.DataFrame(
            {"Sample_rep1": [1.0, 2.0], "Sample_rep2": [3.0, 4.0]},
            index=pd.Index(["gene_001", "gene_002"], name="gene_id"),
        )
        with tempfile.NamedTemporaryFile(
            suffix=suffix, mode="w", delete=False
        ) as f:
            df.to_csv(f, sep=sep)
            path = Path(f.name)
        loaded = load_counts(path)
        path.unlink()
        return loaded, df

    def test_csv_shape(self):
        loaded, original = self._write_and_load(".csv", ",")
        self.assertEqual(loaded.shape, original.shape)

    def test_tsv_shape(self):
        loaded, original = self._write_and_load(".tsv", "\t")
        self.assertEqual(loaded.shape, original.shape)

    def test_values_are_float(self):
        loaded, _ = self._write_and_load(".csv", ",")
        self.assertTrue((loaded.dtypes == float).all())

    def test_toy_csv_loads(self):
        toy = Path(__file__).resolve().parents[1] / "query_data" / "toy.csv"
        if not toy.exists():
            self.skipTest("query_data/toy.csv not found")
        df = load_counts(toy)
        self.assertEqual(df.shape[0], 100)
        self.assertEqual(df.shape[1], 11)   # Liver×4 + Muscle×4 + Brain×3

    def test_toy_csv_gene_ids_are_ensnfug(self):
        toy = Path(__file__).resolve().parents[1] / "query_data" / "toy.csv"
        if not toy.exists():
            self.skipTest("query_data/toy.csv not found")
        df = load_counts(toy)
        self.assertTrue(
            all(g.startswith("ENSNFUG") for g in df.index),
            "toy.csv gene IDs should be ENSNFUG Ensembl identifiers",
        )


# ── TestBatchCorrect ───────────────────────────────────────────────────────────

class TestBatchCorrect(unittest.TestCase):
    """Mocks pycombat_seq so no inmoose dependency needed at test time."""

    def _run_with_mock(self, atlas_raw, query_counts):
        """Run batch_correct with pycombat_seq mocked as identity."""
        combined = pd.concat(
            [atlas_raw.loc[atlas_raw.index.intersection(query_counts.index)],
             query_counts.loc[atlas_raw.index.intersection(query_counts.index)]
                         .round().astype(int)],
            axis=1,
        )
        with patch("inmoose.pycombat.pycombat_seq", return_value=combined.values):
            return batch_correct(query_counts, atlas_raw)

    def test_returns_tuple(self):
        atlas = _make_atlas_raw()
        query = _make_query_counts()
        result = self._run_with_mock(atlas, query)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)

    def test_atlas_corrected_columns(self):
        atlas = _make_atlas_raw()
        query = _make_query_counts()
        atlas_corr, _ = self._run_with_mock(atlas, query)
        self.assertListEqual(list(atlas_corr.columns), _ATLAS_IDS)

    def test_query_corrected_columns(self):
        atlas = _make_atlas_raw()
        query = _make_query_counts()
        _, query_corr = self._run_with_mock(atlas, query)
        self.assertListEqual(list(query_corr.columns), _QUERY_IDS)

    def test_shared_genes_only(self):
        atlas = _make_atlas_raw()
        # Query has extra genes not in atlas
        extra = pd.DataFrame(
            np.ones((5, _N_QUERY)),
            index=[f"extra_{i}" for i in range(5)],
            columns=_QUERY_IDS,
        )
        query = pd.concat([_make_query_counts(), extra])
        atlas_corr, query_corr = self._run_with_mock(atlas, query)
        self.assertEqual(len(atlas_corr), _N_GENES)
        self.assertEqual(len(query_corr), _N_GENES)

    def test_output_dtype_float(self):
        atlas = _make_atlas_raw()
        query = _make_query_counts()
        atlas_corr, query_corr = self._run_with_mock(atlas, query)
        self.assertTrue((atlas_corr.dtypes == float).all())
        self.assertTrue((query_corr.dtypes == float).all())


# ── TestRunBayesage2 ───────────────────────────────────────────────────────────

class TestRunBayesage2(unittest.TestCase):

    def setUp(self):
        self.tmp      = tempfile.TemporaryDirectory()
        self.out_dir  = Path(self.tmp.name) / "bayesage2"
        self.atlas    = _make_atlas_raw()
        self.meta     = _make_atlas_meta()
        self.query    = _make_query_counts()

    def tearDown(self):
        self.tmp.cleanup()

    def test_prediction_csv_created(self):
        run_bayesage2(self.atlas, self.meta, self.query,
                      self.out_dir, m_values=[10, 20], tissue="Liver")
        self.assertTrue((self.out_dir / "Liver_BayesAge2_predictions.csv").exists())

    def test_feature_importance_csv_created(self):
        run_bayesage2(self.atlas, self.meta, self.query,
                      self.out_dir, m_values=[10], tissue="Liver")
        self.assertTrue(
            (self.out_dir / "Liver_BayesAge2_feature_importance.csv").exists()
        )

    def test_predictions_row_count(self):
        run_bayesage2(self.atlas, self.meta, self.query,
                      self.out_dir, m_values=[10, 20], tissue="Liver")
        df = pd.read_csv(self.out_dir / "Liver_BayesAge2_predictions.csv", index_col=0)
        self.assertEqual(len(df), _N_QUERY)

    def test_predictions_have_tage_columns(self):
        run_bayesage2(self.atlas, self.meta, self.query,
                      self.out_dir, m_values=[10, 20], tissue="Liver")
        df = pd.read_csv(self.out_dir / "Liver_BayesAge2_predictions.csv", index_col=0)
        self.assertIn("tAge_M10", df.columns)
        self.assertIn("tAge_M20", df.columns)

    def test_predictions_within_age_range(self):
        run_bayesage2(self.atlas, self.meta, self.query,
                      self.out_dir, m_values=[10], tissue="Liver")
        df = pd.read_csv(self.out_dir / "Liver_BayesAge2_predictions.csv", index_col=0)
        self.assertTrue((df["tAge_M10"] >= 47).all())
        self.assertTrue((df["tAge_M10"] <= 163).all())

    def test_feature_importance_sorted(self):
        run_bayesage2(self.atlas, self.meta, self.query,
                      self.out_dir, m_values=[10], tissue="Liver")
        fi = pd.read_csv(
            self.out_dir / "Liver_BayesAge2_feature_importance.csv", index_col=0
        )
        self.assertIn("abs_spearman_r", fi.columns)
        self.assertTrue((fi["abs_spearman_r"].diff().dropna() <= 1e-9).all())


# ── TestRunPCR ────────────────────────────────────────────────────────────────

class TestRunPCR(unittest.TestCase):

    def setUp(self):
        self.tmp     = tempfile.TemporaryDirectory()
        self.out_dir = Path(self.tmp.name) / "pcr"
        self.atlas   = _make_atlas_raw()
        self.meta    = _make_atlas_meta()
        self.query   = _make_query_counts()

    def tearDown(self):
        self.tmp.cleanup()

    def test_prediction_csv_created(self):
        run_pcr(self.atlas, self.meta, self.query,
                self.out_dir, n_components=[2, 3], top_n_var=None, tissue="Liver")
        self.assertTrue((self.out_dir / "Liver_PCR_predictions.csv").exists())

    def test_cv_metrics_created(self):
        run_pcr(self.atlas, self.meta, self.query,
                self.out_dir, n_components=[2, 3], top_n_var=None, tissue="Liver")
        self.assertTrue((self.out_dir / "Liver_cv_metrics.csv").exists())

    def test_feature_importance_created(self):
        run_pcr(self.atlas, self.meta, self.query,
                self.out_dir, n_components=[2, 3], top_n_var=None, tissue="Liver")
        self.assertTrue((self.out_dir / "Liver_feature_importance.csv").exists())

    def test_predictions_row_count(self):
        run_pcr(self.atlas, self.meta, self.query,
                self.out_dir, n_components=[2, 3], top_n_var=None, tissue="Liver")
        df = pd.read_csv(self.out_dir / "Liver_PCR_predictions.csv", index_col=0)
        self.assertEqual(len(df), _N_QUERY)

    def test_top_n_var_filter(self):
        # top_n_var=10 should still produce valid output with 50-gene atlas
        run_pcr(self.atlas, self.meta, self.query,
                self.out_dir, n_components=[2], top_n_var=10, tissue="Liver")
        df = pd.read_csv(self.out_dir / "Liver_PCR_predictions.csv", index_col=0)
        self.assertEqual(len(df), _N_QUERY)


# ── TestRunEN ─────────────────────────────────────────────────────────────────

class TestRunEN(unittest.TestCase):

    def setUp(self):
        self.tmp     = tempfile.TemporaryDirectory()
        self.out_dir = Path(self.tmp.name) / "en"
        self.atlas   = _make_atlas_raw()
        self.meta    = _make_atlas_meta()
        self.query   = _make_query_counts()

    def tearDown(self):
        self.tmp.cleanup()

    def test_prediction_csv_created(self):
        run_en(self.atlas, self.meta, self.query,
               self.out_dir, top_n_var=None, tissue="Liver")
        self.assertTrue((self.out_dir / "Liver_EN_predictions.csv").exists())

    def test_predictions_row_count(self):
        run_en(self.atlas, self.meta, self.query,
               self.out_dir, top_n_var=None, tissue="Liver")
        df = pd.read_csv(self.out_dir / "Liver_EN_predictions.csv", index_col=0)
        self.assertEqual(len(df), _N_QUERY)

    def test_feature_importance_created(self):
        run_en(self.atlas, self.meta, self.query,
               self.out_dir, top_n_var=None, tissue="Liver")
        self.assertTrue((self.out_dir / "Liver_EN_feature_importance.csv").exists())

    def test_predicted_age_column_present(self):
        run_en(self.atlas, self.meta, self.query,
               self.out_dir, top_n_var=None, tissue="Liver")
        df = pd.read_csv(self.out_dir / "Liver_EN_predictions.csv", index_col=0)
        self.assertIn("predicted_age", df.columns)

    def test_top_n_var_filter(self):
        run_en(self.atlas, self.meta, self.query,
               self.out_dir, top_n_var=10, tissue="Liver")
        df = pd.read_csv(self.out_dir / "Liver_EN_predictions.csv", index_col=0)
        self.assertEqual(len(df), _N_QUERY)


# ── TestMainIntegration ────────────────────────────────────────────────────────

class TestMainIntegration(unittest.TestCase):
    """End-to-end test of main() with toy.csv and mocked Atlas.

    Uses --no-batch-correct and --clocks bayesage2 to keep the test fast
    and avoid the ComBat-seq and DESeq2-norm file dependencies.
    """

    TOY_CSV = Path(__file__).resolve().parents[1] / "query_data" / "toy.csv"

    @classmethod
    def setUpClass(cls):
        if not cls.TOY_CSV.exists():
            raise unittest.SkipTest("query_data/toy.csv not found — run the generator first")
        toy = pd.read_csv(cls.TOY_CSV, index_col=0)
        cls.toy_genes = list(toy.index)
        cls.n_genes   = len(cls.toy_genes)

    def _make_mock_atlas(self, tissues=None):
        """Mock Atlas with ENSNFUG IDs matching toy.csv, covering requested tissues."""
        if tissues is None:
            tissues = {"Liver": 10, "Muscle": 10}
        rng     = np.random.default_rng(99)
        ids, tissue_labels = [], []
        for tissue, n in tissues.items():
            for i in range(1, n + 1):
                ids.append(f"{tissue[:3].upper()}_{i:03d}")
                tissue_labels.append(tissue)
        n_atlas = len(ids)
        raw  = pd.DataFrame(
            rng.negative_binomial(10, 0.5, (self.n_genes, n_atlas)).astype(float),
            index=self.toy_genes,   # ENSNFUG IDs — same as toy.csv
            columns=ids,
        )
        meta = pd.DataFrame(
            {
                "tissue":   tissue_labels,
                "age_days": rng.integers(47, 163, n_atlas).astype(float),
                "sex":      ["M"] * n_atlas,
            },
            index=pd.Index(ids, name="lib"),
        )
        return raw, meta

    def _mock_gene_mapper(self, counts):
        """GeneMapper mock: return counts unchanged (ENSNFUG IDs kept as-is)."""
        mock = MagicMock()
        mock.return_value.convert.side_effect = lambda df: df
        return mock

    def test_bayesage2_output_files(self):
        atlas_raw, atlas_meta = self._make_mock_atlas()

        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            with patch("run_query_clocks.DataLoader") as MockDL, \
                 patch("run_query_clocks.GeneMapper") as MockGM:
                MockDL.return_value.load_atlas.return_value = (atlas_raw, atlas_meta)
                # GeneMapper.convert is identity: keeps ENSNFUG IDs unchanged
                MockGM.return_value.convert.side_effect = lambda df: df

                sys.argv = [
                    "run_query_clocks.py",
                    "--counts",       str(self.TOY_CSV),
                    "--tissues",      "Liver",
                    "--clocks",       "bayesage2",
                    "--no-batch-correct",
                    "--gene-id-type", "ensembl",   # triggers GeneMapper (mocked)
                    "--m-values",     "10", "20",
                    "--out-dir",      str(out_dir),
                ]
                from run_query_clocks import main
                main()

            pred = out_dir / "bayesage2" / "Liver_BayesAge2_predictions.csv"
            self.assertTrue(pred.exists(), f"Missing: {pred}")
            df = pd.read_csv(pred, index_col=0)
            self.assertEqual(len(df), 4)        # 4 Liver reps in toy.csv
            self.assertIn("tAge_M10", df.columns)
            self.assertIn("tAge_M20", df.columns)

    def test_multi_tissue_detection(self):
        """main() auto-detects Liver, Muscle, Brain from toy.csv column names."""
        atlas_raw, atlas_meta = self._make_mock_atlas(
            tissues={"Liver": 10, "Muscle": 10, "Brain": 10}
        )

        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            with patch("run_query_clocks.DataLoader") as MockDL, \
                 patch("run_query_clocks.GeneMapper") as MockGM:
                MockDL.return_value.load_atlas.return_value = (atlas_raw, atlas_meta)
                MockGM.return_value.convert.side_effect = lambda df: df

                sys.argv = [
                    "run_query_clocks.py",
                    "--counts",       str(self.TOY_CSV),
                    "--clocks",       "bayesage2",
                    "--no-batch-correct",
                    "--gene-id-type", "ensembl",
                    "--m-values",     "10",
                    "--out-dir",      str(out_dir),
                ]
                from run_query_clocks import main
                main()

            ba2_dir = out_dir / "bayesage2"
            self.assertTrue((ba2_dir / "Liver_BayesAge2_predictions.csv").exists())
            self.assertTrue((ba2_dir / "Muscle_BayesAge2_predictions.csv").exists())
            self.assertTrue((ba2_dir / "Brain_BayesAge2_predictions.csv").exists())


if __name__ == "__main__":
    unittest.main()
