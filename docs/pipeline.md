# Pipeline Overview

All run scripts share the same internal flow. They loop over tissues and apply the steps below.

## Architecture

### 1. Load Atlas reference

- **`DataLoader.load_atlas()`** — loads raw counts + metadata from `data/`
- **`Preprocessor.filter_genes(min_count=1)`** — removes all-zero genes
- **`Preprocessor.stratify(tissue, sex=None)`** — subsets to the target tissue

### 2. Load query data

Two paths depending on dataset format:

**Generic CSV/TSV** (`run_query_clocks.py`, `run_PRJNA817434_clocks.py`)

- `load_counts(path)` — reads genes × samples matrix
- Gene ID handling (auto-detected from index):
    - ENSNFUG IDs → `GeneMapper.convert()` renames rows to Atlas gene names
    - Atlas-style names (symbols / LOC107) → no mapping; intersected by name directly
- Column names parsed for tissue + replicate, or supplied via `--metadata`

**AAlab xlsx** (`run_AAlab_clocks.py`, `run_Eugen_clocks.py`)

- `QueryCountExtractor.extract_tissue(tissue)` — parses xlsx DE result files, deduplicates samples, applies `GeneMapper`

### 3. Batch correction (default; skip with `--no-batch-correct`)

- `ComBat-seq` (`inmoose.pycombat_seq`) applied jointly to Atlas + query counts
- Atlas is the reference batch — preserves its scale
- Returns batch-corrected counts for both Atlas and query, restricted to shared genes
- Corrected Atlas counts are used for all clock training so that training and test data share the same corrected feature space

### 4. Run clocks

Three clocks run in parallel over the same corrected gene set:

**BayesAge2**

- Input: raw integer counts (rounds batch-corrected floats)
- `BayesAge2Clock(lowess_top_n=250).build_reference(atlas_raw, atlas_meta)` — LOWESS fits + Spearman ranking internally apply frequency normalization
- `.predict(query, n_genes=M)` for each M in 25…200 (step 5)
- Outputs: `tAge_M{m}` columns per sample

**PCR**

- Input: batch-corrected counts (float)
- `PCRClock.loso_cv(atlas, meta)` — LOSO-CV across n_components candidates; selects optimal n by R²
- `.fit(atlas, meta)` — refit on full Atlas at optimal n
- `.predict(query)` — single prediction column `predicted_age_n{optimal_n}`
- `.save_loadings()` — saves CV metrics, gene loadings, feature importance

**Elastic Net**

- Input: batch-corrected counts (float)
- `ElasticNetClock.tune_and_train(atlas, meta)` — GridSearchCV + LOO-CV on Atlas
- `.predict(query)` — single `predicted_age` column

When `--no-batch-correct` is used, PCR and EN use the pre-saved DESeq2-normalized Atlas
(`outputs/normalized/Atlas_DESeq2_normalized.csv`) instead of the batch-corrected counts.

### 5. Save outputs

Each runner saves predictions + feature importance to disk under the output directory.

## Output files

| Sub-directory | File | Contents |
|---|---|---|
| `bayesage2/` | `{tissue}_BayesAge2_predictions.csv` | `tAge_M{m}` columns, index=`sample_id` |
| `bayesage2/` | `{tissue}_BayesAge2_feature_importance.csv` | Genes ranked by \|Spearman r\| |
| `bayesage2/references/` | `{tissue}_reference.tsv` | Full BayesAge2 reference table |
| `pcr/` | `{tissue}_PCR_predictions.csv` | `predicted_age_n{optimal_n}` column |
| `pcr/` | `{tissue}_PCR_cv_metrics.csv` | LOSO-CV R² and MAE per n_components |
| `pcr/` | `{tissue}_PCR_feature_importance.csv` | Gene importance scores |
| `pcr/` | `{tissue}_PCR_loadings.tsv` | Per-component gene loadings |
| `elastic_net/` | `{tissue}_EN_predictions.csv` | `predicted_age` column |
| `elastic_net/` | `{tissue}_EN_feature_importance.csv` | Non-zero EN coefficients |

## Run scripts

| Script | Dataset type | Gene ID format |
|---|---|---|
| `scripts/run_query_clocks.py` | Any genes × samples CSV/TSV | ENSNFUG or Atlas-style (auto-detected) |
| `scripts/run_PRJNA817434_clocks.py` | PRJNA817434 Fat tissue (NCBI GTF output) | Atlas-style (LOC107 / symbol); no mapping |
| `scripts/run_Eugen_clocks.py` | Eugene killifish xlsx (Gut/Kidney/Spleen) | ENSNFUG via `GeneMapper` |
| `scripts/run_AAlab_clocks.py` | AAlab xlsx datasets | ENSNFUG via `GeneMapper` |

## Step-by-step

### Step 1 — Normalize the Atlas reference (run once)

```bash
python scripts/normalize_reference.py
```

Outputs:

- `outputs/normalized/Atlas_freq_normalized.csv` — frequency-normalized (BayesAge 2.0 reference)
- `outputs/normalized/Atlas_DESeq2_normalized.csv` — DESeq2 size-factor normalized (PCR / EN fallback)

Only needed when running PCR or EN with `--no-batch-correct`.

### Step 2 — Apply clocks to query data

```bash
# Generic CSV — ENSNFUG IDs, column names: TISSUE_repN
python scripts/run_query_clocks.py --counts my_counts.csv

# Generic CSV — Atlas gene names already, skip mapping
python scripts/run_query_clocks.py --counts my_counts.csv --gene-id-type atlas

# PRJNA817434 (column names don't follow TISSUE_repN convention)
python scripts/run_PRJNA817434_clocks.py
```

See [run_query_clocks.py](scripts/run_query_clocks.md) for the full option reference.

## Key design decisions

### ComBat-seq batch correction

ComBat-seq is applied jointly to Atlas + query before clock training so that training and test
data share the same corrected feature space. Atlas is the reference batch, preserving its scale.
Disable with `--no-batch-correct` if query counts are already on the Atlas scale.

### M-value range (BayesAge 2.0)

The default range 25–200 (step 5) provides broader exploration than the original (5–100, step 5).
Configurable via `--m-values` in `run_query_clocks.py` or `M_VALUES` at the top of the
dataset-specific scripts.

### Gene ID handling

`GeneMapper` is only invoked for ENSNFUG datasets. Datasets from `raw_RNAseq_process`
(NCBI GTF pipeline) already use Atlas gene names and are intersected directly. See
[Gene ID Mapping](data/gene_mapping.md) for coverage statistics.
