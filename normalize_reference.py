"""
normalize_reference.py — Run normalization on the full Atlas reference dataset.

Produces two output files in outputs/normalized/:
  Atlas_freq_normalized.csv   — frequency normalization (input for BayesAge 2.0)
  Atlas_DESeq2_normalized.csv — DESeq2 size-factor normalization (input for EN / PCR)

"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))

from data_loader import DataLoader
from normalization import FrequencyNormalizer, DESeq2Normalizer
from preprocessing import Preprocessor

OUT_DIR = Path(__file__).resolve().parent / "outputs" / "normalized"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FREQ_OUT = OUT_DIR / "Atlas_freq_normalized.csv"
DESEQ2_OUT = OUT_DIR / "Atlas_DESeq2_normalized.csv"


def main():
    # --- Load ---
    print("=== Loading Atlas data ===")
    dl = DataLoader()
    counts, meta = dl.load_atlas()
    print(f"Raw counts: {counts.shape[0]:,} genes × {counts.shape[1]:,} samples")

    # --- Gene filter (remove zero-count genes) ---
    pp = Preprocessor(min_count=1)
    counts = pp.filter_genes(counts)

    # --- Frequency normalization (BayesAge 2.0) ---
    print("\n=== Frequency normalization ===")
    freq = FrequencyNormalizer().normalize(counts)
    freq.to_csv(FREQ_OUT)
    print(f"Saved → {FREQ_OUT}")
    print(f"Shape: {freq.shape[0]:,} genes × {freq.shape[1]:,} samples")
    print(f"Column sum check (first 3): {freq.sum(axis=0).iloc[:3].round(6).tolist()}")

    # --- DESeq2 normalization (EN / PCR) ---
    # Design: tissue (categorical), matching original R script design intent.
    # Size factors are estimated from median-of-ratios across all samples.
    print("\n=== DESeq2 normalization (design: tissue) ===")
    dn = DESeq2Normalizer(design_factors=["tissue"])
    norm = dn.normalize(counts, meta, save_path=DESEQ2_OUT)
    print(f"Shape: {norm.shape[0]:,} genes × {norm.shape[1]:,} samples")
    print(f"Size factor range: {dn.size_factors.min():.3f} – {dn.size_factors.max():.3f}")

    print("\n=== Done ===")
    print(f"Freq normalized : {FREQ_OUT}")
    print(f"DESeq2 normalized: {DESEQ2_OUT}")


if __name__ == "__main__":
    main()
