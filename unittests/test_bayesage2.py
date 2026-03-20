"""Unit tests for bayesage2.py using toy data."""

import sys
import unittest
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from bayesage2 import BayesAge2Clock, AGE_MIN, AGE_MAX


def make_toy_data(n_genes=30, n_samples=12, seed=42):
    """
    Toy counts with a weak age signal: a few genes scale linearly with age.
    Ages span 47–162 days.
    """
    rng = np.random.default_rng(seed)
    ages = np.linspace(AGE_MIN, AGE_MAX, n_samples).astype(int)
    counts = rng.integers(50, 300, size=(n_genes, n_samples)).astype(int)

    # Make first 5 genes correlated with age (positive)
    for i in range(5):
        counts[i] = (ages * rng.uniform(0.5, 2.0) + rng.integers(10, 50)).astype(int)
    # Make next 5 genes anti-correlated with age
    for i in range(5, 10):
        counts[i] = (ages[::-1] * rng.uniform(0.5, 2.0) + rng.integers(10, 50)).astype(int)

    genes = [f"gene_{i}" for i in range(n_genes)]
    samples = [f"s{i}" for i in range(n_samples)]
    counts_df = pd.DataFrame(counts, index=genes, columns=samples)
    meta_df = pd.DataFrame({"age_days": ages.astype(float)}, index=samples)
    return counts_df, meta_df


class TestBayesAge2BuildReference(unittest.TestCase):

    def setUp(self):
        self.counts, self.meta = make_toy_data()
        self.clock = BayesAge2Clock(lowess_top_n=15)

    def test_reference_has_correct_columns(self):
        ref = self.clock.build_reference(self.counts, self.meta)
        self.assertIn("spearman_r", ref.columns)
        self.assertIn("p_value", ref.columns)
        age_cols = [f"age_{a}" for a in range(AGE_MIN, AGE_MAX + 1)]
        for col in age_cols:
            self.assertIn(col, ref.columns)

    def test_reference_index_matches_genes(self):
        ref = self.clock.build_reference(self.counts, self.meta)
        self.assertListEqual(ref.index.tolist(), self.counts.index.tolist())

    def test_spearman_values_in_valid_range(self):
        ref = self.clock.build_reference(self.counts, self.meta)
        valid = ref["spearman_r"].dropna()
        self.assertTrue((valid >= -1).all() and (valid <= 1).all())

    def test_top_genes_have_highest_correlation(self):
        ref = self.clock.build_reference(self.counts, self.meta)
        top = ref["spearman_r"].abs().nlargest(10).index.tolist()
        # The correlated genes (gene_0..9) should dominate the top
        n_signal = sum(g in top for g in [f"gene_{i}" for i in range(10)])
        self.assertGreaterEqual(n_signal, 5)

    def test_save_and_load_reference(self):
        import tempfile, os
        ref = self.clock.build_reference(self.counts, self.meta)
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
            path = Path(f.name)
        try:
            self.clock.build_reference(self.counts, self.meta, save_path=path)
            loaded_clock = BayesAge2Clock.load_reference(path)
            pd.testing.assert_frame_equal(ref, loaded_clock._reference)
        finally:
            path.unlink(missing_ok=True)


class TestBayesAge2Predict(unittest.TestCase):

    def setUp(self):
        self.counts, self.meta = make_toy_data()
        self.clock = BayesAge2Clock(lowess_top_n=15)
        self.clock.build_reference(self.counts, self.meta)

    def test_predict_returns_series_with_correct_index(self):
        preds = self.clock.predict(self.counts, n_genes=5)
        self.assertIsInstance(preds, pd.Series)
        self.assertListEqual(preds.index.tolist(), self.counts.columns.tolist())

    def test_predicted_ages_within_range(self):
        preds = self.clock.predict(self.counts, n_genes=5)
        self.assertTrue((preds >= AGE_MIN).all())
        self.assertTrue((preds <= AGE_MAX).all())

    def test_predict_without_reference_raises(self):
        clock = BayesAge2Clock()
        with self.assertRaises(RuntimeError):
            clock.predict(self.counts, n_genes=5)


class TestBayesAge2LosoCv(unittest.TestCase):

    def setUp(self):
        self.counts, self.meta = make_toy_data(n_samples=8)
        self.clock = BayesAge2Clock(lowess_top_n=15)

    def test_loso_cv_output_shape(self):
        results = self.clock.loso_cv(self.counts, self.meta, m_values=[5, 10])
        self.assertEqual(len(results), self.counts.shape[1])
        self.assertIn("tAge_M5", results.columns)
        self.assertIn("tAge_M10", results.columns)

    def test_loso_cv_index_matches_samples(self):
        results = self.clock.loso_cv(self.counts, self.meta, m_values=[5])
        self.assertListEqual(
            sorted(results.index.tolist()),
            sorted(self.counts.columns.tolist()),
        )

    def test_loso_cv_ages_within_range(self):
        results = self.clock.loso_cv(self.counts, self.meta, m_values=[5])
        self.assertTrue((results["tAge_M5"] >= AGE_MIN).all())
        self.assertTrue((results["tAge_M5"] <= AGE_MAX).all())


if __name__ == "__main__":
    unittest.main()
