"""Unit tests for normalization.py using toy data."""

import sys
import unittest
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from normalization import FrequencyNormalizer, DESeq2Normalizer


def make_toy_counts(n_genes=20, n_samples=10, seed=0):
    """Small integer count matrix with known properties."""
    rng = np.random.default_rng(seed)
    counts = rng.integers(1, 500, size=(n_genes, n_samples)).astype(int)
    genes = [f"gene_{i}" for i in range(n_genes)]
    samples = [f"s{i}" for i in range(n_samples)]
    return pd.DataFrame(counts, index=genes, columns=samples)


def make_toy_metadata(n_samples=10, seed=0):
    rng = np.random.default_rng(seed)
    ages = rng.integers(47, 163, size=n_samples).astype(float)
    samples = [f"s{i}" for i in range(n_samples)]
    return pd.DataFrame({"age_days": ages}, index=samples)


class TestFrequencyNormalizer(unittest.TestCase):

    def setUp(self):
        self.counts = make_toy_counts()
        self.fn = FrequencyNormalizer()

    def test_output_shape(self):
        freq = self.fn.normalize(self.counts)
        self.assertEqual(freq.shape, self.counts.shape)

    def test_column_sums_to_one(self):
        freq = self.fn.normalize(self.counts)
        col_sums = freq.sum(axis=0)
        np.testing.assert_allclose(col_sums.values, np.ones(self.counts.shape[1]), atol=1e-10)

    def test_values_between_zero_and_one(self):
        freq = self.fn.normalize(self.counts)
        self.assertTrue((freq.values >= 0).all())
        self.assertTrue((freq.values <= 1).all())

    def test_preserves_index_and_columns(self):
        freq = self.fn.normalize(self.counts)
        self.assertListEqual(freq.index.tolist(), self.counts.index.tolist())
        self.assertListEqual(freq.columns.tolist(), self.counts.columns.tolist())

    def test_zero_column_handled(self):
        counts = self.counts.copy()
        counts["s0"] = 0  # all-zero sample
        # Should not raise; result will be NaN for that column
        freq = self.fn.normalize(counts)
        self.assertTrue(np.isnan(freq["s0"]).all())


class TestDESeq2Normalizer(unittest.TestCase):

    def setUp(self):
        self.counts = make_toy_counts(n_genes=20, n_samples=10)
        self.meta = make_toy_metadata(n_samples=10)
        self.dn = DESeq2Normalizer(design_factors=["age_days"])

    def test_output_shape(self):
        norm = self.dn.normalize(self.counts, self.meta)
        self.assertEqual(norm.shape, self.counts.shape)

    def test_output_index_and_columns(self):
        norm = self.dn.normalize(self.counts, self.meta)
        self.assertListEqual(norm.index.tolist(), self.counts.index.tolist())
        self.assertListEqual(norm.columns.tolist(), self.counts.columns.tolist())

    def test_size_factors_available(self):
        self.dn.normalize(self.counts, self.meta)
        sf = self.dn.size_factors
        self.assertIsNotNone(sf)
        self.assertEqual(len(sf), self.counts.shape[1])
        self.assertTrue((sf > 0).all())

    def test_normalized_values_are_positive(self):
        norm = self.dn.normalize(self.counts, self.meta)
        self.assertTrue((norm.values >= 0).all())


if __name__ == "__main__":
    unittest.main()
