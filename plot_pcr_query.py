"""
plot_pcr_query.py — Scatter-box plots of query age predictions (BayesAge2, PCR, EN).

Usage:
    python plot_pcr_query.py

Plots only Query samples. For each tissue × model:
  - X-axis: conditions in biological order
  - Y-axis: predicted tAge (days)
  - Boxes split by Young / Old, individual dots labeled by replicate

Optimal parameter selection:
  BayesAge2 : first M where mean predictions stabilise (< 3 day change over 3 steps)
  PCR       : n_components with lowest Mann-Whitney p-value (from mw_pvals.csv)
  EN        : single predicted_age column (no selection needed)

Outputs: outputs/figures/{Tissue}_combined_{Model}_query_<param>.svg
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ── Config ───────────────────────────────────────────────────────────────────
TISSUES = [] #NOTE
MODELS  = ["BayesAge2", "PCR", "EN"]

OUT_BASE = Path(__file__).resolve().parent / "outputs"
FIG_DIR  = OUT_BASE / "figures"

CONDITION_ORDER  = [] #NOTE
CONDITION_LABELS = {} #NOTE
AGE_ORDER  = ["Young", "Old"]
AGE_COLORS = {"Young": "#4C72B0", "Old": "#DD8452"}
# ─────────────────────────────────────────────────────────────────────────────


# ── File paths ────────────────────────────────────────────────────────────────

def _paths(tissue: str, model: str) -> tuple[Path, Path | None]:
    """Return (pred_path, optional_param_path) for the given tissue + model."""
    prefix = f"{tissue}_sexcombined"
    T = ""
    if model == "BayesAge2":
        return OUT_BASE / "bayesage2" / f"{prefix}_BayesAge2_{T}query.csv", None
    if model == "PCR":
        return (
            OUT_BASE / "pcr" / f"{prefix}_PCR_{T}query.csv",
            OUT_BASE / "pcr" / f"{prefix}_PCR_{T}query_mw_pvals.csv",
        )
    if model == "EN":
        return OUT_BASE / "elastic_net" / f"{prefix}_EN_{T}query_loso.csv", None
    raise ValueError(f"Unknown model: {model}")


# ── Optimal parameter selection ───────────────────────────────────────────────

def _pick_stable_m(df: pd.DataFrame) -> int:
    """BayesAge2: first M where mean predictions stabilise.

    Stability = mean absolute change across all samples between consecutive
    M values drops below 3 days for 3 consecutive steps.

    Degenerate M values (std < 1 day across all samples) are excluded first,
    as they indicate a collapsed prediction where all samples return the same
    integer age due to a degenerate argmax in the Poisson likelihood.
    """
    m_cols = [c for c in df.columns if c.startswith("tAge_M")]
    m_vals = sorted(int(c.replace("tAge_M", "")) for c in m_cols)

    # Exclude M values where all samples collapse to the same prediction
    valid_m = [m for m in m_vals if df[f"tAge_M{m}"].std() >= 1.0]
    if not valid_m:
        return m_vals[-1]

    means = [df[f"tAge_M{m}"].mean() for m in valid_m]
    threshold = 3.0
    stable_streak = 0

    for i in range(1, len(valid_m)):
        change = abs(means[i] - means[i - 1])
        stable_streak = stable_streak + 1 if change < threshold else 0
        if stable_streak >= 3:
            return valid_m[i - 2]   # first M of the stable run

    return valid_m[-1]   # fallback: largest valid M


def _pick_optimal_n(mw_path: Path) -> int:
    """PCR: n_components with lowest Mann-Whitney p-value."""
    mw = pd.read_csv(mw_path, index_col=0)["mw_pval"]
    return int(mw.idxmin())


def _select_column(df: pd.DataFrame, model: str, param_path: Path | None) -> tuple[str, str]:
    """Return (tage_column_name, subtitle_string) for the model."""
    if model == "BayesAge2":
        m = _pick_stable_m(df)
        return f"tAge_M{m}", f"optimal M = {m}"
    if model == "PCR":
        n = _pick_optimal_n(param_path)
        mw = pd.read_csv(param_path, index_col=0)["mw_pval"]
        p  = mw[n]
        return f"tAge_n{n}", f"n_components = {n}  |  MW p = {p:.3g}"
    if model == "EN":
        return "predicted_age", ""
    raise ValueError(model)


# ── Core plot function ────────────────────────────────────────────────────────

def _parse_rep(sample_id: str) -> str:
    """Extract 'Rep_N' from sample_id for dot labels."""
    return sample_id.split(".")[-1]



def plot_model(tissue: str, model: str) -> None:
    """Load data, select optimal parameter, draw and save figure."""
    pred_path, param_path = _paths(tissue, model)
    if not pred_path.exists():
        print(f"  [SKIP] {tissue} {model}: {pred_path.name} not found.")
        return

    df = pd.read_csv(pred_path)

    # EN: filter to Query samples only
    if model == "EN":
        df = df[df["source"] == "Query"].copy()

    # Drop rows with missing metadata
    df = df.dropna(subset=["age_group", "condition"])

    tage_col, subtitle = _select_column(df, model, param_path)
    conditions = [c for c in CONDITION_ORDER if c in df["condition"].unique()]
    n_cond = len(conditions)

    age_groups = [g for g in AGE_ORDER if g in df["age_group"].unique()]
    n_age = len(age_groups)

    fig, axes = plt.subplots(1, n_age, figsize=(2.8 * n_cond, 4.5), sharey=False)
    if n_age == 1:
        axes = [axes]

    title = f"{model} Predicted Age — {tissue}"
    if subtitle:
        title += f"\n{subtitle}"
    fig.suptitle(title, fontsize=11)

    cond_labels = [CONDITION_LABELS.get(c, c) for c in conditions]

    for ax, age_group in zip(axes, age_groups):
        sub = df[df["age_group"] == age_group].copy()
        sub["condition"] = pd.Categorical(sub["condition"], categories=conditions, ordered=True)
        sub = sub.sort_values("condition")

        sns.boxplot(
            data=sub, x="condition", y=tage_col,
            order=conditions, color=AGE_COLORS[age_group],
            width=0.45, fliersize=0, linewidth=1.2, ax=ax,
        )
        sns.stripplot(
            data=sub, x="condition", y=tage_col,
            order=conditions, color=AGE_COLORS[age_group],
            size=6, jitter=True,
            edgecolor="white", linewidth=0.6, ax=ax,
        )
        for _, row in sub.iterrows():
            x_pos = conditions.index(row["condition"])
            ax.text(
                x_pos, row[tage_col] + 0.3, _parse_rep(row["sample_id"]),
                ha="center", va="bottom", fontsize=6, color="#555555",
            )

        ax.set_title(age_group, fontsize=10, color=AGE_COLORS[age_group], fontweight="bold")
        ax.set_xlabel("")
        ax.set_xticks(range(len(conditions)))
        ax.set_xticklabels(cond_labels, fontsize=9)
        ax.set_ylabel("Predicted tAge (days)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    param_str = tage_col.replace("tAge_", "").replace("predicted_age", "en")
    out_path = FIG_DIR / f"{tissue}_combined_{model}_query_{param_str}.svg"
    fig.savefig(out_path, bbox_inches="tight")
    print(f"  Saved → {out_path.name}")
    plt.close(fig)


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    for tissue in TISSUES:
        for model in MODELS:
            print(f"{tissue} / {model}")
            plot_model(tissue, model)
