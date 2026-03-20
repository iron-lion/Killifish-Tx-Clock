"""Unit tests for pcr.py using toy data.

n_components_range uses [2, 5] during dev/test (per CLAUDE.md).
"""

import sys
import unittest
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pcr import PCRClock


def make_toy_data(n_genes=60, n_samples=15, seed=0):
    """Toy normalized counts with a linear age signal."""
    rng = np.random.default_rng(seed)
    ages = np.linspace(47, 162, n_samples).astype(float)
    counts = rng.uniform(1, 100, size=(n_genes, n_samples))
    # 10 genes correlated with age
    for i in range(10):
        counts[i] = ages * rng.uniform(0.3, 1.5) + rng.uniform(5, 20)
    genes = [f"gene_{i}" for i in range(n_genes)]
    samples = [f"s{i}" for i in range(n_samples)]
    counts_df = pd.DataFrame(counts, index=genes, columns=samples)
    meta_df = pd.DataFrame({"age_days": ages}, index=samples)
    return counts_df, meta_df


class TestPCRClockLosoCv(unittest.TestCase):

    def setUp(self):
        self.counts, self.meta = make_toy_data()
        # 2 values per param per CLAUDE.md
        self.clock = PCRClock(n_components_range=[2, 5])

    def test_loso_cv_returns_correct_columns(self):
        results = self.clock.loso_cv(self.counts, self.meta)
        self.assertIn("age_days", results.columns)
        self.assertIn("pred_ncomp_2", results.columns)
        self.assertIn("pred_ncomp_5", results.columns)

    def test_loso_cv_index_matches_samples(self):
        results = self.clock.loso_cv(self.counts, self.meta)
        self.assertListEqual(
            sorted(results.index.tolist()),
            sorted(self.counts.columns.tolist()),
        )

    def test_loso_cv_sets_optimal_n(self):
        self.clock.loso_cv(self.counts, self.meta)
        self.assertIn(self.clock.optimal_n_components, [2, 5])

    def test_cv_metrics_shape(self):
        self.clock.loso_cv(self.counts, self.meta)
        metrics = self.clock.cv_metrics
        self.assertIsNotNone(metrics)
        self.assertEqual(len(metrics), 2)
        self.assertIn("MAE", metrics.columns)
        self.assertIn("R2", metrics.columns)
        self.assertIn("pearson_r", metrics.columns)


class TestPCRClockFitPredict(unittest.TestCase):

    def setUp(self):
        self.counts, self.meta = make_toy_data()
        self.clock = PCRClock(n_components_range=[2, 5])
        self.clock.loso_cv(self.counts, self.meta)

    def test_fit_with_optimal_n(self):
        self.clock.fit(self.counts, self.meta)
        self.assertIsNotNone(self.clock._pipeline)

    def test_fit_with_explicit_n(self):
        self.clock.fit(self.counts, self.meta, n_components=2)
        self.assertIsNotNone(self.clock._pipeline)

    def test_predict_shape_and_index(self):
        self.clock.fit(self.counts, self.meta)
        preds = self.clock.predict(self.counts)
        self.assertEqual(len(preds), self.counts.shape[1])
        self.assertListEqual(preds.index.tolist(), self.counts.columns.tolist())

    def test_predict_without_fit_raises(self):
        clock = PCRClock()
        with self.assertRaises(RuntimeError):
            clock.predict(self.counts)

    def test_fit_without_loso_cv_raises(self):
        clock = PCRClock()
        with self.assertRaises(RuntimeError):
            clock.fit(self.counts, self.meta)


class TestPCRClockFeatureImportance(unittest.TestCase):

    def setUp(self):
        self.counts, self.meta = make_toy_data()
        self.clock = PCRClock(n_components_range=[2, 5])
        self.clock.loso_cv(self.counts, self.meta)
        self.clock.fit(self.counts, self.meta)

    def test_fi_length_matches_genes(self):
        fi = self.clock.get_feature_importance()
        self.assertEqual(len(fi), self.counts.shape[0])

    def test_fi_index_matches_genes(self):
        fi = self.clock.get_feature_importance()
        self.assertSetEqual(set(fi.index), set(self.counts.index))

    def test_fi_sorted_by_abs_descending(self):
        fi = self.clock.get_feature_importance()
        abs_vals = fi.abs().values
        self.assertTrue(all(abs_vals[i] >= abs_vals[i + 1] for i in range(len(abs_vals) - 1)))

    def test_fi_without_fit_raises(self):
        clock = PCRClock()
        with self.assertRaises(RuntimeError):
            clock.get_feature_importance()


class TestPCRClockVarFilter(unittest.TestCase):

    def test_top_n_var_genes_filter(self):
        counts, meta = make_toy_data(n_genes=60)
        clock = PCRClock(n_components_range=[2], top_n_var_genes=20)
        clock.loso_cv(counts, meta)
        self.assertEqual(len(clock._selected_genes), 20)

    def test_save_loadings_creates_files(self):
        import tempfile
        counts, meta = make_toy_data()
        clock = PCRClock(n_components_range=[2, 5])
        clock.loso_cv(counts, meta)
        clock.fit(counts, meta)
        with tempfile.TemporaryDirectory() as tmpdir:
            clock.save_loadings(Path(tmpdir), prefix="test")
            self.assertTrue((Path(tmpdir) / "test_loadings.tsv").exists())
            self.assertTrue((Path(tmpdir) / "test_feature_importance.csv").exists())
            self.assertTrue((Path(tmpdir) / "test_cv_metrics.csv").exists())


if __name__ == "__main__":
    unittest.main()
