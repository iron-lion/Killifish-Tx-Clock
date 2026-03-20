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


unittests/                  # Unit tests (pytest)
Costa_et_al/                # Original KillifishAtlas analysis scripts
data_matrices/              # Atlas count/TPM matrices and metadata
query_data/                 # Query xlsx DE result files
outputs/                    # Generated outputs
environment.yml             # Conda environment definition
```

---

## Setup

```bash
conda env create -f environment.yml
conda activate killifish-tx-clock
```

---


## Running the Pipeline

```bash
# 1. Normalize the full Atlas reference (run once)
python src/normalize_reference.py

# 2. Apply clocks to query data
#    Check query_data/toy.csv
python run_query_clocks.py

# 3. Plot results
python src/plot_pcr_query.py
```

---

## Pipeline Overview

`run_query_clocks.py` loops over tissues and executes the following steps:

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

| Directory | File pattern | Contents |
|---|---|---|
| `outputs/bayesage2/` | `{tissue}_sexcombined_BayesAge2_query.csv` | tAge predictions at each M value |
| `outputs/bayesage2/` | `{tissue}_sexcombined_BayesAge2_feature_importance.csv` | Genes ranked by \|Spearman r\| |
| `outputs/pcr/` | `{tissue}_sexcombined_PCR_query.csv` | tAge predictions per n_components |
| `outputs/pcr/` | `{tissue}_sexcombined_PCR_query_mw_pvals.csv` | Mann-Whitney U p-values (Young vs Old) |
| `outputs/pcr/` | `{tissue}_sexcombined_PCR_feature_importance_n{n}.csv` | Gene importance per n_components |
| `outputs/elastic_net/` | `{tissue}_sexcombined_EN_query_loso.csv` | Atlas LOSO-CV + query predictions |
| `outputs/elastic_net/` | `{tissue}_sexcombined_EN_feature_importance.csv` | Non-zero EN coefficients |
| `outputs/figures/` | `{tissue}_combined_{Model}_query_{param}.svg` | Scatter-box plots per tissue/model |

---
