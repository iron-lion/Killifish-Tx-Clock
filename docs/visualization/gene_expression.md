# plot_gene_expression.py

Plot TPM expression of a single gene across all Atlas tissues, with optional
per-tissue age timecourse plots.

**File:** `plot_gene_expression.py`

## Usage

```bash
# Basic: violin + strip plot across all tissues
python plot_gene_expression.py actb

# Log-scale
python plot_gene_expression.py tp53 --log

# Also generate per-tissue age timecourse plots
python plot_gene_expression.py mb21d1 --timecourse

# Full options
python plot_gene_expression.py actb \
    --log \
    --out figures/actb_expression.svg \
    --timecourse \
    --outdir figures/actb_timecourse/
```

## Arguments

| Argument | Default | Description |
|---|---|---|
| `gene_id` | *(required)* | Gene name or ID (case-insensitive, partial match supported) |
| `--log` | off | Plot log₂(TPM + 1) instead of raw TPM |
| `--out PATH` | `<gene>_expression.svg` | Output path for the across-tissues summary plot |
| `--timecourse` | off | Also generate per-tissue age × expression scatter plots |
| `--outdir DIR` | `<gene>_timecourse/` | Directory for timecourse SVGs |

## Data sources

| File | Description |
|---|---|
| `data/GSE308970_TPM_Atlas_allbatches_merged_v3.csv` | TPM expression matrix (genes × samples) |
| `data/ExperimentDesign_allbatches_combined_v7.csv` | Sample metadata (tissue, age_days, sex) |

## Gene lookup

- **Exact match** — case-insensitive comparison against the full gene index (~25K genes).
- **Partial match fallback** — substring search if no exact match found.
- **Multiple matches** — uses the first match and prints a warning listing alternatives.

## Summary plot (across tissues)

One violin + strip plot per tissue, ordered by body region:

```
Brain → Eye → SpinalCord → Heart → Lung →
Liver → Gut → Fat → Kidney → Spleen →
Muscle → Skin → Bone → Testis → Ovary
```

Each violin shows the expression distribution; dots show individual samples.
Tissues not present in the data are skipped automatically.
Sample count (`n=`) is annotated below each violin.

## Timecourse plots (`--timecourse`)

One SVG per tissue. Each plot shows:

- **X-axis:** Age (days)
- **Y-axis:** TPM or log₂(TPM + 1)
- **Colors:** Female (#E07B8A) / Male (#5B8DB8)
- **LOWESS trend line** — drawn per sex when ≥ 5 data points are available (`frac=0.6`)

Single-sex tissues (Testis, Ovary) omit the sex legend automatically.

## Example output

```
Gene: actb  |  650 samples across 15 tissues
Saved → actb_expression.svg
Generating timecourse plots → actb_timecourse/
  Saved → actb_Brain_timecourse.svg
  Saved → actb_Eye_timecourse.svg
  ...
  Saved → actb_Ovary_timecourse.svg
```
