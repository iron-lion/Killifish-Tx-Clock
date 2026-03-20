"""Unit tests for elastic_net.py using toy data.

GridSearch uses 2 values per param during dev/test (per CLAUDE.md).
"""

import sys
import unittest
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from elastic_net import ElasticNetClock


# 2-value dev grids (per CLAUDE.md)
DEV_ALPHA_GRID = [0.1, 1.0]
DEV_L1_RATIO_GRID = [0.5, 1.0]


def make_toy_data(n_genes=50, n_samples=15, seed=0):
    """Toy normalized counts with a linear age signal in a few genes."""
    rng = np.random.default_rng(seed)
    ages = np.linspace(47, 162, n_samples).astype(float)
    counts = rng.uniform(1, 100, size=(n_genes, n_samples))

    # Make first 10 genes correlated with age
    for i in range(10):
        counts[i] = ages * rng.uniform(0.3, 1.5) + rng.uniform(5, 20)

    genes = [f"gene_{i}" for i in range(n_genes)]
    samples = [f"s{i}" for i in range(n_samples)]
    counts_df = pd.DataFrame(counts, index=genes, columns=samples)
    meta_df = pd.DataFrame({"age_days": ages}, index=samples)
    return counts_df, meta_df


class TestElasticNetClock(unittest.TestCase):

    def setUp(self):
        self.counts, self.meta = make_toy_data()
        self.en = ElasticNetClock(tissue="Liver", top_n_var_genes=None)
        # Patch grids to dev-size for speed
        import elastic_net as en_mod
        self._orig_alpha = en_mod.ALPHA_GRID
        self._orig_l1 = en_mod.L1_RATIO_GRID
        en_mod.ALPHA_GRID = DEV_ALPHA_GRID
        en_mod.L1_RATIO_GRID = DEV_L1_RATIO_GRID

    def tearDown(self):
        import elastic_net as en_mod
        en_mod.ALPHA_GRID = self._orig_alpha
        en_mod.L1_RATIO_GRID = self._orig_l1

    def test_tune_and_train_returns_best_params(self):
        params = self.en.tune_and_train(self.counts, self.meta)
        self.assertIn("alpha", params)
        self.assertIn("l1_ratio", params)
        self.assertIn(params["alpha"], DEV_ALPHA_GRID)
        self.assertIn(params["l1_ratio"], DEV_L1_RATIO_GRID)

    def test_nonzero_genes_are_subset_of_input(self):
        self.en.tune_and_train(self.counts, self.meta)
        nonzero = self.en.nonzero_genes
        self.assertIsNotNone(nonzero)
        for gene in nonzero.index:
            self.assertIn(gene, self.counts.index)

    def test_predict_shape(self):
        self.en.tune_and_train(self.counts, self.meta)
        preds = self.en.predict(self.counts)
        self.assertEqual(len(preds), self.counts.shape[1])
        self.assertListEqual(preds.index.tolist(), self.counts.columns.tolist())

    def test_loso_cv_output_shape(self):
        self.en.tune_and_train(self.counts, self.meta)
        cv = self.en.loso_cv(self.counts, self.meta)
        self.assertEqual(len(cv), self.counts.shape[1])
        self.assertIn("age_days", cv.columns)
        self.assertIn("predicted_age", cv.columns)

    def test_top_n_var_genes_filter(self):
        en_filtered = ElasticNetClock(tissue="Liver", top_n_var_genes=20)
        import elastic_net as en_mod
        en_mod.ALPHA_GRID = DEV_ALPHA_GRID
        en_mod.L1_RATIO_GRID = DEV_L1_RATIO_GRID
        en_filtered.tune_and_train(self.counts, self.meta)
        self.assertIsNotNone(en_filtered._selected_genes)
        self.assertEqual(len(en_filtered._selected_genes), 20)

    def test_predict_without_training_raises(self):
        en = ElasticNetClock()
        with self.assertRaises(RuntimeError):
            en.predict(self.counts)

    def test_save_creates_files(self):
        import tempfile, os
        self.en.tune_and_train(self.counts, self.meta)
        with tempfile.TemporaryDirectory() as tmpdir:
            self.en.save(Path(tmpdir), prefix="test")
            self.assertTrue((Path(tmpdir) / "test_best_params.json").exists())
            self.assertTrue((Path(tmpdir) / "test_non_zero_genes.csv").exists())


if __name__ == "__main__":
    unittest.main()
