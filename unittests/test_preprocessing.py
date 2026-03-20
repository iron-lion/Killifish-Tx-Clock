"""Unit tests for preprocessing.py using toy data."""

import sys
import unittest
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from preprocessing import Preprocessor


def make_toy_counts(n_genes=20, n_samples=12, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.integers(1, 300, size=(n_genes, n_samples)).astype(int)
    genes = [f"gene_{i}" for i in range(n_genes)]
    samples = [f"s{i}" for i in range(n_samples)]
    return pd.DataFrame(counts, index=genes, columns=samples)


def make_toy_metadata(n_samples=12):
    """6 Liver/M, 3 Liver/F, 3 Brain/M samples with ages."""
    ages = [50, 70, 90, 110, 130, 150, 60, 80, 100, 55, 75, 95]
    tissues = ["Liver"] * 6 + ["Liver"] * 3 + ["Brain"] * 3
    sexes = ["M"] * 6 + ["F"] * 3 + ["M"] * 3
    samples = [f"s{i}" for i in range(n_samples)]
    return pd.DataFrame({"age_days": ages, "tissue": tissues, "sex": sexes}, index=samples)


class TestFilterGenes(unittest.TestCase):

    def setUp(self):
        self.pp = Preprocessor(min_count=1)

    def test_removes_zero_genes(self):
        counts = make_toy_counts()
        # Force two genes to all-zero
        counts.iloc[0] = 0
        counts.iloc[3] = 0
        filtered = self.pp.filter_genes(counts)
        self.assertEqual(filtered.shape[0], counts.shape[0] - 2)

    def test_keeps_all_when_expressed(self):
        counts = make_toy_counts()
        filtered = self.pp.filter_genes(counts)
        self.assertEqual(filtered.shape[0], counts.shape[0])

    def test_output_columns_unchanged(self):
        counts = make_toy_counts()
        filtered = self.pp.filter_genes(counts)
        self.assertListEqual(filtered.columns.tolist(), counts.columns.tolist())


class TestStratify(unittest.TestCase):

    def setUp(self):
        self.pp = Preprocessor()
        self.counts = make_toy_counts()
        self.meta = make_toy_metadata()

    def test_stratify_tissue_only(self):
        c, m = self.pp.stratify(self.counts, self.meta, tissue="Liver")
        self.assertEqual(c.shape[1], 9)  # 6 M + 3 F Liver
        self.assertTrue((m["tissue"] == "Liver").all())

    def test_stratify_tissue_and_sex(self):
        c, m = self.pp.stratify(self.counts, self.meta, tissue="Liver", sex="M")
        self.assertEqual(c.shape[1], 6)
        self.assertTrue((m["sex"] == "M").all())

    def test_case_insensitive(self):
        c1, _ = self.pp.stratify(self.counts, self.meta, tissue="liver", sex="m")
        c2, _ = self.pp.stratify(self.counts, self.meta, tissue="Liver", sex="M")
        self.assertEqual(c1.shape, c2.shape)

    def test_counts_and_meta_aligned(self):
        c, m = self.pp.stratify(self.counts, self.meta, tissue="Brain")
        self.assertListEqual(sorted(c.columns.tolist()), sorted(m.index.tolist()))


class TestDetectOutliers(unittest.TestCase):

    def setUp(self):
        self.pp = Preprocessor()

    def test_no_outliers_normal_data(self):
        rng = np.random.default_rng(1)
        counts = pd.DataFrame(
            rng.integers(100, 200, size=(30, 15)),
            index=[f"g{i}" for i in range(30)],
            columns=[f"s{i}" for i in range(15)],
        )
        # Use n_sd=3: at 2 SD, random small samples can have false positives by chance
        clean, outliers = self.pp.detect_outliers(counts, n_sd=3)
        self.assertEqual(len(outliers), 0)
        self.assertEqual(clean.shape[1], 15)

    def test_detects_extreme_outlier(self):
        rng = np.random.default_rng(2)
        counts = pd.DataFrame(
            rng.integers(100, 200, size=(50, 15)),
            index=[f"g{i}" for i in range(50)],
            columns=[f"s{i}" for i in range(15)],
        )
        # Make one sample an extreme outlier
        counts["s0"] = counts["s0"] * 1000
        clean, outliers = self.pp.detect_outliers(counts, n_sd=2)
        self.assertIn("s0", outliers)

    def test_output_shape_consistent(self):
        rng = np.random.default_rng(3)
        counts = pd.DataFrame(
            rng.integers(10, 100, size=(20, 10)),
            index=[f"g{i}" for i in range(20)],
            columns=[f"s{i}" for i in range(10)],
        )
        clean, outliers = self.pp.detect_outliers(counts)
        self.assertEqual(clean.shape[1], 10 - len(outliers))
        self.assertEqual(clean.shape[0], counts.shape[0])


if __name__ == "__main__":
    unittest.main()
