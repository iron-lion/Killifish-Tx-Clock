"""
plot_gene_expression.py — Plot TPM expression of a gene across Atlas tissues.

Usage:
    python plot_gene_expression.py <gene_id> [--log] [--out PATH] [--timecourse] [--outdir DIR]

Arguments:
    gene_id         Gene name/ID as it appears in the TPM matrix row index (case-insensitive).
    --log           Plot log2(TPM + 1) instead of raw TPM.
    --out PATH      Output file for the across-tissues summary plot (default: <gene_id>_expression.svg).
    --timecourse    Also generate per-tissue scatter plots (age vs expression, sex-colored).
    --outdir DIR    Directory for per-tissue timecourse plots (default: <gene_id>_timecourse/).

Examples:
    python plot_gene_expression.py myc
    python plot_gene_expression.py tp53 --log --timecourse
    python plot_gene_expression.py actb --timecourse --outdir figures/actb/
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from statsmodels.nonparametric.smoothers_lowess import lowess

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR  = Path(__file__).resolve().parents[1] / "data"
TPM_PATH  = DATA_DIR / "GSE308970_TPM_Atlas_allbatches_merged_v3.csv"
META_PATH = DATA_DIR / "ExperimentDesign_allbatches_combined_v7.csv"

# Tissue display order (roughly body regions)
TISSUE_ORDER = [
    "Brain", "Eye", "SpinalCord",
    "Heart", "Lung",
    "Liver", "Gut", "Fat",
    "Kidney", "Spleen",
    "Muscle", "Skin", "Bone",
    "Testis", "Ovary",
]
SEX_COLORS = {"F": "#E07B8A", "M": "#5B8DB8"}
SEX_LABELS = {"F": "Female", "M": "Male"}
# ─────────────────────────────────────────────────────────────────────────────


def load_data(gene_id: str) -> tuple[pd.DataFrame, str]:
    """Load TPM row for gene_id and merge with tissue metadata.

    Returns (merged_df, matched_gene_name).
    Raises SystemExit if the gene is not found.
    """
    # Load metadata
    meta = pd.read_csv(META_PATH, index_col=0, usecols=["Unnamed: 0", "tissue", "age_days", "sex"])

    # Load only the gene row from TPM (scan index first for speed)
    tpm_index = pd.read_csv(TPM_PATH, index_col=0, usecols=[0]).index

    # Case-insensitive lookup
    matches = [g for g in tpm_index if g.lower() == gene_id.lower()]
    if not matches:
        # Partial match fallback
        matches = [g for g in tpm_index if gene_id.lower() in g.lower()]
    if not matches:
        print(f"ERROR: Gene '{gene_id}' not found in TPM matrix.", file=sys.stderr)
        print(f"  Matrix has {len(tpm_index):,} genes.", file=sys.stderr)
        sys.exit(1)
    if len(matches) > 1:
        print(f"WARNING: Multiple matches for '{gene_id}': {matches[:5]}. Using '{matches[0]}'.")

    matched = matches[0]

    # Read just that row
    tpm_full = pd.read_csv(TPM_PATH, index_col=0)
    expr = tpm_full.loc[matched]  # Series: sample_id → TPM

    # Merge with metadata
    df = meta.copy()
    df["TPM"] = expr.reindex(df.index)
    df = df.dropna(subset=["TPM", "tissue"])

    return df, matched


def make_plot(df: pd.DataFrame, gene_name: str, log_scale: bool, out_path: Path) -> None:
    # Keep only tissues present in the data; respect TISSUE_ORDER
    present = df["tissue"].unique()
    order = [t for t in TISSUE_ORDER if t in present] + \
            sorted(t for t in present if t not in TISSUE_ORDER)

    y_col = "TPM"
    y_label = "TPM"
    if log_scale:
        df = df.copy()
        df["log2TPM"] = np.log2(df["TPM"] + 1)
        y_col = "log2TPM"
        y_label = "log₂(TPM + 1)"

    n_tissues = len(order)
    fig, ax = plt.subplots(figsize=(max(8, n_tissues * 0.9), 5))

    # Tissue-specific colors (palette cycles if >14 tissues)
    palette = sns.color_palette("tab20", n_tissues)
    color_map = dict(zip(order, palette))

    df = df.copy()
    df["tissue"] = pd.Categorical(df["tissue"], categories=order, ordered=True)

    sns.violinplot(
        data=df, x="tissue", y=y_col, hue="tissue",
        order=order, hue_order=order, palette=color_map, fill=False,
        width=0.55, linewidth=1.0, ax=ax, legend=False, #fliersize=0, 
    )
    sns.stripplot(
        data=df, x="tissue", y=y_col, hue="tissue",
        order=order, hue_order=order, palette=color_map,
        size=4, jitter=True, alpha=0.6,
        edgecolor="white", linewidth=0.3, ax=ax, legend=False,
    )

    # Sample-count annotation below each tick
    for i, tissue in enumerate(order):
        n = (df["tissue"] == tissue).sum()
        ax.text(i, ax.get_ylim()[0], f"n={n}", ha="center", va="top",
                fontsize=7, color="#666666")

    ax.set_title(f"{gene_name}  —  expression across tissues", fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel(y_label, fontsize=10)
    ax.set_ylim((0,350))
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=40, ha="right", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)

    csv_path = out_path.with_suffix(".csv")
    df[["tissue", "sex", "age_days", y_col]].to_csv(csv_path)
    print(f"Saved → {out_path}  +  {csv_path.name}")


def make_timecourse_plots(
    df: pd.DataFrame,
    gene_name: str,
    log_scale: bool,
    out_dir: Path,
) -> None:
    """One scatter plot per tissue: age_days (x) vs expression (y), sex-colored.

    A LOWESS trend line is drawn per sex when ≥ 5 data points are available.
    """
    present = df["tissue"].unique()
    tissues = [t for t in TISSUE_ORDER if t in present] + \
              sorted(t for t in present if t not in TISSUE_ORDER)

    y_col = "TPM"
    y_label = "TPM"
    if log_scale:
        df = df.copy()
        df["log2TPM"] = np.log2(df["TPM"] + 1)
        y_col = "log2TPM"
        y_label = "log₂(TPM + 1)"

    out_dir.mkdir(parents=True, exist_ok=True)

    for tissue in tissues:
        sub = df[df["tissue"] == tissue].dropna(subset=["age_days", y_col, "sex"])
        if sub.empty:
            continue

        sexes_present = [s for s in ["F", "M"] if s in sub["sex"].values]

        fig, ax = plt.subplots(figsize=(5.5, 4))

        for sex in sexes_present:
            grp = sub[sub["sex"] == sex].sort_values("age_days")
            color = SEX_COLORS[sex]
            label = SEX_LABELS[sex]

            ax.scatter(
                grp["age_days"], grp[y_col],
                color=color, label=label,
                s=30, alpha=0.75, edgecolors="white", linewidths=0.4, zorder=3,
            )

            # LOWESS trend line (need ≥ 5 points for a meaningful smooth)
            if len(grp) >= 5:
                smoothed = lowess(
                    grp[y_col].values, grp["age_days"].values,
                    frac=0.6, return_sorted=True,
                )
                ax.plot(
                    smoothed[:, 0], smoothed[:, 1],
                    color=color, linewidth=1.8, alpha=0.9, zorder=4,
                )

        ax.set_title(f"{gene_name}  —  {tissue}", fontsize=11)
        ax.set_xlabel("Age (days)", fontsize=10)
        ax.set_ylabel(y_label, fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if len(sexes_present) > 1:
            ax.legend(title="Sex", fontsize=9, title_fontsize=9, frameon=False)

        plt.tight_layout()
        out_path = out_dir / f"{gene_name}_{tissue}_timecourse.svg"
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)

        csv_path = out_path.with_suffix(".csv")
        sub[["age_days", "sex", y_col]].to_csv(csv_path)
        print(f"  Saved → {out_path.name}  +  {csv_path.name}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot TPM expression of a gene across Atlas tissues."
    )
    parser.add_argument("gene_id", help="Gene name (case-insensitive, partial match supported)")
    parser.add_argument("--log",       action="store_true", help="Plot log2(TPM+1)")
    parser.add_argument("--out",       type=Path, default=None, help="Output path for summary plot (.svg/.png)")
    parser.add_argument("--timecourse", action="store_true", help="Also generate per-tissue age vs expression plots")
    parser.add_argument("--outdir",    type=Path, default=None, help="Directory for timecourse plots")
    args = parser.parse_args()

    df, matched = load_data(args.gene_id)
    print(f"Gene: {matched}  |  {len(df)} samples across {df['tissue'].nunique()} tissues")

    out_path = args.out or Path(f"{matched}_expression.svg")
    make_plot(df, matched, log_scale=args.log, out_path=out_path)

    if args.timecourse:
        out_dir = args.outdir or Path(f"{matched}_timecourse")
        print(f"Generating timecourse plots → {out_dir}/")
        make_timecourse_plots(df, matched, log_scale=args.log, out_dir=out_dir)


if __name__ == "__main__":
    main()
