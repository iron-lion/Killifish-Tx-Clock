# run_query_clocks.py

General-purpose CLI for applying all three aging clocks to any genes × samples count matrix.

## Input format

A CSV or TSV where:

- Rows are genes (Ensembl `ENSNFUG...` IDs or Atlas gene names)
- Columns are samples named `TISSUE_repN` (e.g. `Liver_rep1`, `Muscle_Rep3`, `SpinalCord_rep2`)

Column names are parsed to extract the tissue label automatically.
An optional `--metadata` CSV can override tissue assignments and provide `age_days`.

## Usage

Scripts are in `scripts/`. Run from the repo root:

```bash
# Minimal — all clocks, auto-detect tissues, all available genes:
python scripts/run_query_clocks.py --counts query_data/toy.csv

# Select tissues and clocks, skip batch correction:
python scripts/run_query_clocks.py --counts my_counts.csv \
    --tissues Liver Muscle \
    --clocks bayesage2 pcr \
    --no-batch-correct

# Provide metadata explicitly:
python scripts/run_query_clocks.py --counts my_counts.csv \
    --metadata my_meta.csv   # columns: sample_id, tissue [, age_days]

# Atlas gene names already (skip Ensembl mapping):
python scripts/run_query_clocks.py --counts my_counts.csv \
    --gene-id-type atlas

# Custom output directory:
python scripts/run_query_clocks.py --counts my_counts.csv --out-dir results/
```

## Options

| Flag | Default | Description |
|---|---|---|
| `--counts` | *(required)* | genes × samples CSV or TSV |
| `--metadata` | — | optional sample metadata CSV (`sample_id`, `tissue`, `age_days`) |
| `--tissues` | auto-detected | Atlas tissue labels to include |
| `--clocks` | `bayesage2 pcr en` | subset of clocks to run |
| `--gene-id-type` | `auto` | `ensembl` (ENSNFUG IDs) / `atlas` (gene names) / `auto` |
| `--no-batch-correct` | off | skip ComBat-seq batch correction |
| `--m-values` | `25..200 step 5` | BayesAge2 gene-set sizes |
| `--n-components` | `5 10 15 20` | PCR components tested via LOSO-CV |
| `--top-n-var` | all genes | pre-filter to top-N most variable genes |
| `--out-dir` | `outputs/` | output base directory |

## Outputs

```
outputs/
  bayesage2/{tissue}_BayesAge2_predictions.csv
  bayesage2/{tissue}_BayesAge2_feature_importance.csv
  bayesage2/references/{tissue}_reference.tsv
  pcr/{tissue}_PCR_predictions.csv
  pcr/{tissue}_PCR_cv_metrics.csv
  pcr/{tissue}_PCR_feature_importance.csv
  pcr/{tissue}_PCR_loadings.tsv
  elastic_net/{tissue}_EN_predictions.csv
  elastic_net/{tissue}_EN_feature_importance.csv
```

## Gene ID handling

`run_query_clocks.py` supports two gene ID formats:

| Gene ID type | Example | What happens |
|---|---|---|
| `ensembl` | `ENSNFUG00015001234` | `GeneMapper.convert()` renames rows to Atlas names; unmapped genes dropped |
| `atlas` | `actb`, `LOC107374091` | No mapping — genes intersected with Atlas by name directly |

When `--gene-id-type auto` (default), the script inspects the first 30 index entries:

- Any `ENSNFUG`-prefixed entry → `ensembl` mode
- Otherwise → `atlas` mode

Pass `--gene-id-type ensembl` or `--gene-id-type atlas` to force one mode.

### When to use `--gene-id-type atlas`

Datasets produced by `raw_RNAseq_process/run_rnaseq.sh` (NCBI GTF pipeline) already
use Atlas-compatible NCBI gene names. Pass `--gene-id-type atlas` (or rely on
auto-detection) — `GeneMapper` is skipped and genes are matched by name.

Note: `LOC139XXXXXX` genes produced by the NfurGRZ-RIMD1 GTF have **0 % Atlas
coverage** (newer NCBI annotation not present in the Atlas). They are silently dropped
at the intersection step. See [Gene ID Mapping](../data/gene_mapping.md) for coverage details.

### PRJNA817434

PRJNA817434 uses Atlas-style gene IDs but its column names (`Old_wt_fed_fat_1`) do not
follow the `TISSUE_repN` convention expected by `run_query_clocks.py`. Use the dedicated
wrapper instead:

```bash
python scripts/run_PRJNA817434_clocks.py
```

## Prerequisites

Run once before using PCR or EN:

```bash
python scripts/normalize_reference.py
```
