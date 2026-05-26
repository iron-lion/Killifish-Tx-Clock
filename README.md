# Killifish Transcriptome Aging Clock

An adaptation of the KillifishAtlas aging clock pipeline for [AAlab](https://www.age.mpg.de/antebi). Applies three transcriptomic aging clocks — **BayesAge 2.0**, **Elastic Net (EN)**, and **Principal Component Regression (PCR)** — to *N. furzeri* RNA-seq data and predicts transcriptomic age (tAge) in query samples.

---

## Repository Structure

```
src/                        # Core pipeline modules
├── data_loader.py          # Load Atlas counts, metadata, and query files
├── normalization.py        # Frequency and DESeq2 normalization
├── preprocessing.py        # Gene filtering, stratification, outlier detection
├── bayesage2.py            # BayesAge 2.0 clock
├── elastic_net.py          # Elastic Net clock
├── pcr.py                  # Principal Component Regression clock
├── calibration.py          # Apply clocks to query datasets
├── gene_mapping.py         # Map query Ensembl IDs to Atlas gene names
└── normalize_reference.py  # Run normalization on the full Atlas dataset

run_query_clocks.py         # General-purpose CLI: apply clocks to any count matrix
run_AAlab_clocks.py         # AAlab fasting/refeeding experiment (xlsx files, Fat/Liver/Muscle)
run_Eugen_clocks.py         # Eugene experiment (CSV, Gut/Kidney/Spleen, WT vs KO vs IR)

unittests/                  # Unit tests (pytest)
Costa_et_al/                # Original KillifishAtlas analysis scripts
data/                       # Atlas count/TPM matrices, metadata, and gene mapping
query_data/                 # Query input files (xlsx DE results, CSV count matrices)
outputs/                    # Generated outputs
environment.yml             # Conda environment definition
```

---

## Setup
Download KilifishAtlas data from [GEO PRJNA1274512](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308970)

```bash
conda env create -f environment.yml
conda activate killifish-tx-clock
```

---


## Running the Pipeline

### Step 1 — Normalize the Atlas reference (run once)

```bash
python src/normalize_reference.py
# → outputs/normalized/Atlas_DESeq2_normalized.csv
```

### Step 2 — Apply clocks to query data

Three run scripts are provided. All three share the same underlying clocks and
output format; they differ only in how the query data is loaded and parsed.

---

#### `run_query_clocks.py` — general-purpose CLI

Accepts any genes × samples CSV or TSV. Column names must follow
`TISSUE_repN` (e.g. `Liver_rep1`, `Muscle_Rep3`).

```bash
# All clocks, auto-detect tissues and gene ID type:
python run_query_clocks.py --counts query_data/toy.csv

# Select tissues, skip batch correction, BayesAge2 + PCR only:
python run_query_clocks.py --counts my_counts.csv \
    --tissues Liver Muscle \
    --clocks bayesage2 pcr \
    --no-batch-correct

# Provide sample metadata explicitly:
python run_query_clocks.py --counts my_counts.csv \
    --metadata my_meta.csv          # columns: sample_id, tissue [, age_days]

# Custom output directory:
python run_query_clocks.py --counts my_counts.csv --out-dir results/
```

Key options:

| Flag | Default | Description |
|---|---|---|
| `--counts` | *(required)* | genes × samples CSV/TSV |
| `--metadata` | — | optional sample metadata CSV (`sample_id`, `tissue`, `age_days`) |
| `--tissues` | auto-detected | Atlas tissue labels to run |
| `--clocks` | `bayesage2 pcr en` | subset of clocks to run |
| `--gene-id-type` | `auto` | `ensembl` (ENSNFUG IDs) / `atlas` / `auto` |
| `--no-batch-correct` | off | skip ComBat-seq correction |
| `--m-values` | 25..200 step 5 | BayesAge2 gene-set sizes |
| `--n-components` | 5 10 15 20 | PCR components tested via LOSO-CV |
| `--top-n-var` | all genes | pre-filter to top-N variable genes |
| `--out-dir` | `outputs/` | output base directory |

---

#### `run_AAlab_clocks.py` — AAlab fasting/refeeding experiment

Reads the DE result xlsx files in `query_data/` (Fat, Liver, Muscle tissues;
conditions: 72H\_FASTED, 6H\_REFED, 24H\_REFED; Young vs Old).
Set `TISSUES` inside the script before running.

```bash
python run_AAlab_clocks.py
# → outputs/{bayesage2,pcr,elastic_net}/{tissue}_sexcombined_*.csv
```

---

#### `run_Eugen_clocks.py` — Eugene experiment

Reads `query_data/eugene_killifish.csv` (Gut, Kidney, Spleen tissues;
conditions: WT, cGAS\_KO, STINGg\_KO; age groups: Old, Young, Young\_IR).
Column format: `Tissue_Condition_AgeGroup.Rep_N`.

```bash
python run_Eugen_clocks.py
# → outputs/eugen/{bayesage2,pcr,elastic_net}/{tissue}_sexcombined_*.csv
```

### Step 3 — Plot results

```bash
python src/plot_pcr_query.py
```

---

## Pipeline Overview

All three run scripts share the same internal flow. They loop over tissues and execute the following steps:

```
  Atlas (train)                       Query (test)          
                                                            
  DataLoader                          QueryCountExtractor
  └─ raw counts + metadata            └─ parse xlsx/csv files  
     filter_genes (min_count=1)          GeneMapper (ENSNFUG→Atlas)
     Preprocessor.stratify()             extract per-tissue counts 
           │                                    │                 
           └──── ComBat-seq batch correction ───┘
                     (inmoose.pycombat_seq, ref=Atlas)
                                    │
              batch-corrected query counts
                       ┌────────────┤
                       │            │
              ┌────────┘            └───────────────┐
              ▼                                     ▼
   FrequencyNormalize(Atlas raw)       Atlas DESeq2-normalized
              │                        (pre-saved from step 1)
              ▼                                     │
   BayesAge2Clock                       ┌───────────┴──────────┐
   .build_reference(Atlas)              ▼                       ▼
   .predict(query, M=25..200)    PCRClock                ElasticNetClock
              │                  .loso_cv(Atlas)         .tune_and_train(Atlas)
              │                  .fit(Atlas)             .loso_cv(Atlas)
              │                  .predict(query)         .predict(query)
              │                  + Mann-Whitney U        (currently disabled)
              │                  per n_components
              ▼                          ▼                       ▼
   outputs/bayesage2/            outputs/pcr/           outputs/elastic_net/
   *_BayesAge2_query.csv         *_PCR_query.csv        *_EN_query_loso.csv
   *_BayesAge2_feature_          *_PCR_mw_pvals.csv     *_EN_feature_
     importance.csv              *_PCR_feature_           importance.csv
                                   importance_n*.csv
```

### Output Files

AAlab outputs are written to `outputs/{bayesage2,pcr,elastic_net}/`.
Eugene outputs are written to `outputs/eugen/{bayesage2,pcr,elastic_net}/`.
The file naming convention is the same for both.

| Sub-directory | File pattern | Contents |
|---|---|---|
| `bayesage2/` | `{tissue}_sexcombined_BayesAge2_query.csv` | tAge predictions at each M value |
| `bayesage2/` | `{tissue}_sexcombined_BayesAge2_feature_importance.csv` | Genes ranked by \|Spearman r\| |
| `bayesage2/references/` | `{tissue}_sexcombined_query_reference.tsv` | Full BayesAge2 reference table |
| `pcr/` | `{tissue}_sexcombined_PCR_query.csv` | tAge predictions per n_components |
| `pcr/` | `{tissue}_sexcombined_PCR_query_mw_pvals.csv` | Mann-Whitney U p-values (Young vs Old) |
| `pcr/` | `{tissue}_sexcombined_PCR_feature_importance_n{n}.csv` | Gene importance per n_components |
| `elastic_net/` | `{tissue}_sexcombined_EN_query_loso.csv` | Atlas LOSO-CV + query predictions |
| `elastic_net/` | `{tissue}_sexcombined_EN_feature_importance.csv` | Non-zero EN coefficients |
| `figures/` | `{tissue}_combined_{Model}_query_{param}.svg` | Scatter-box plots per tissue/model |

---

## Note
The refactoring has been done with the assistance of CLAUDE.
