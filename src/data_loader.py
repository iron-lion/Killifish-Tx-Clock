"""
DataLoader: load Atlas counts, metadata, and query DE result files.
All paths are relative to the repository root (KillifishAtlas/).
"""

import os
import pandas as pd
from pathlib import Path

# Repository root = one level up from this file
REPO_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = REPO_ROOT / "data_matrices"
QUERY_DIR = REPO_ROOT / "query_data"


class DataLoader:
    """Load and merge Atlas expression data and metadata."""

    COUNTS_FILE = DATA_DIR / "GSE308970_Counts_Atlas_allbatches_merged_v3.csv"
    TPM_FILE = DATA_DIR / "GSE308970_TPM_Atlas_allbatches_merged_v3.csv"
    META_FILE = DATA_DIR / "ExperimentDesign_allbatches_combined_v7.csv"

    def __init__(self):
        self._counts = None
        self._tpm = None
        self._metadata = None

    # ------------------------------------------------------------------
    # Atlas data
    # ------------------------------------------------------------------

    def load_counts(self) -> pd.DataFrame:
        """Load raw count matrix (genes × samples).

        Returns
        -------
        pd.DataFrame
            Index = gene names, columns = sample IDs (lib codes).
        """
        if self._counts is None:
            self._counts = pd.read_csv(self.COUNTS_FILE, index_col=0)
        return self._counts.copy()

    def load_tpm(self) -> pd.DataFrame:
        """Load TPM matrix (genes × samples)."""
        if self._tpm is None:
            self._tpm = pd.read_csv(self.TPM_FILE, index_col=0)
        return self._tpm.copy()

    def load_metadata(self) -> pd.DataFrame:
        """Load sample metadata.

        Returns
        -------
        pd.DataFrame
            Index = lib (sample ID). Key columns: sex, age_days, tissue, cohort.
        """
        if self._metadata is None:
            meta = pd.read_csv(self.META_FILE, index_col=0)
            # Harmonize gonad labels (same as original DESeq2 scripts)
            meta["tissue"] = meta["tissue"].replace({"Ovary": "Gonad", "Testis": "Gonad"})
            # Use lib as the index to match count matrix columns
            meta = meta.set_index("lib")
            self._metadata = meta
        return self._metadata.copy()

    def load_atlas(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return (counts, metadata) aligned to shared samples.

        Keeps only samples present in both the count matrix and metadata.
        """
        counts = self.load_counts()
        meta = self.load_metadata()

        shared = counts.columns.intersection(meta.index)
        counts = counts[shared]
        meta = meta.loc[shared]
        return counts, meta

    # ------------------------------------------------------------------
    # Query data
    # ------------------------------------------------------------------

    def load_query_files(self) -> dict[str, pd.DataFrame]:
        """Load all query .xlsx DE result files.

        Returns
        -------
        dict
            Keys = comparison name (stem of filename), values = DataFrames.
        """
        results = {}
        for path in sorted(QUERY_DIR.glob("*.xlsx")):
            key = path.stem
            df = pd.read_excel(path, index_col=0)
            results[key] = df
        return results

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self) -> None:
        """Print a concise summary of the loaded Atlas data."""
        counts, meta = self.load_atlas()

        print(f"Counts matrix : {counts.shape[0]:,} genes × {counts.shape[1]:,} samples")
        print(f"Metadata      : {meta.shape[0]:,} samples, {meta.shape[1]} columns")
        print()

        print("Samples per tissue:")
        print(meta["tissue"].value_counts().to_string())
        print()

        print("Samples per sex:")
        print(meta["sex"].value_counts().to_string())
        print()

        print(f"Age range (days): {meta['age_days'].min()} – {meta['age_days'].max()}")

        query = self.load_query_files()
        print(f"\nQuery files loaded: {len(query)}")
        for k in query:
            print(f"  {k}")


if __name__ == "__main__":
    dl = DataLoader()
    dl.summary()
