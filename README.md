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


## Note
The refactoring has been done with the assistance of CLAUDE.
