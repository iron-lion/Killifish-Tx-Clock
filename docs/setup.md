# Setup

## Prerequisites

- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/) package manager
- Python 3.11

## Create the environment

```bash
conda env create -f environment.yml
conda activate killifish-tx-clock
```

The environment installs:

| Package | Purpose |
|---|---|
| `numpy`, `pandas`, `scipy` | Core numerics |
| `scikit-learn` | PCA, ElasticNet, GridSearchCV, LOSO-CV |
| `statsmodels` | LOWESS smoothing (BayesAge 2.0 + visualization) |
| `matplotlib`, `seaborn` | Plotting |
| `openpyxl` | Reading `.xlsx` DE result files |
| `pydeseq2` | Python DESeq2 normalization for Atlas reference |
| `inmoose` | ComBat-seq batch correction (`pycombat_seq`) |
| `pybiomart` | (optional) Rebuild Ensembl gene ID mapping |

## Data matrices

The large Atlas data matrices are **not tracked in git**. Place them in `data/`:

| File | Description |
|---|---|
| `GSE308970_TPM_Atlas_allbatches_merged_v3.csv` | TPM matrix (genes × samples) |
| `GSE308970_rawcount_Atlas_allbatches_merged_v3.csv` | Raw count matrix |
| `ExperimentDesign_allbatches_combined_v7.csv` | Sample metadata (tissue, age_days, sex, batch) |

These files are available from GEO accession **GSE308970**.

## One-time normalization

Before running PCR or Elastic Net clocks, pre-compute and cache the DESeq2-normalized
Atlas matrix:

```bash
python src/normalize_reference.py
# → outputs/normalized/Atlas_DESeq2_normalized.csv
# → outputs/normalized/Atlas_freq_normalized.csv
```

This step takes ~5–10 minutes. The output files are loaded automatically by all run scripts.

## Verify installation

```bash
pytest unittests/ -v
```

All 8 test files should pass.
