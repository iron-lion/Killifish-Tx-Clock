"""Unit tests for data_loader.py using real files (smoke tests only)."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from data_loader import DataLoader


class TestDataLoader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dl = DataLoader()

    def test_load_counts_shape(self):
        counts = self.dl.load_counts()
        self.assertGreater(counts.shape[0], 0, "No genes loaded")
        self.assertGreater(counts.shape[1], 0, "No samples loaded")

    def test_load_metadata_index(self):
        meta = self.dl.load_metadata()
        self.assertIn("age_days", meta.columns)
        self.assertIn("sex", meta.columns)
        self.assertIn("tissue", meta.columns)

    def test_load_atlas_alignment(self):
        counts, meta = self.dl.load_atlas()
        # Columns of counts must match index of meta exactly
        self.assertListEqual(sorted(counts.columns.tolist()), sorted(meta.index.tolist()))

    def test_gonad_harmonization(self):
        meta = self.dl.load_metadata()
        # Ovary and Testis should be replaced with Gonad
        self.assertNotIn("Ovary", meta["tissue"].values)
        self.assertNotIn("Testis", meta["tissue"].values)

    def test_load_query_files_returns_dict(self):
        query = self.dl.load_query_files()
        self.assertIsInstance(query, dict)
        self.assertGreater(len(query), 0, "No query files loaded")
        for df in query.values():
            self.assertGreater(df.shape[0], 0)


if __name__ == "__main__":
    unittest.main()
