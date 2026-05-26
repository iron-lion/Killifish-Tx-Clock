# Killifish Transcriptome Aging Clock

Python adaptation of the KillifishAtlas aging clock pipeline for [AAlab](https://www.age.mpg.de/antebi).
Applies three transcriptomic aging clocks — **BayesAge 2.0**, **Elastic Net (EN)**, and **Principal
Component Regression (PCR)** — to *Nothobranchius furzeri* RNA-seq data and predicts transcriptomic
age (tAge) in query samples.

---

## What this repo does

1. **Trains** three aging clocks on the published KillifishAtlas reference dataset
   (677 samples, 15 tissues, ages 47–163 days).
2. **Transfers** the trained models to novel query datasets via ComBat-seq batch correction.
3. **Predicts** transcriptomic age (tAge) for each query sample.
4. **Visualizes** predictions as scatter-box plots and generates per-gene expression
   profiles across tissues and ages.
5. **Preprocessor of RNAseq** for Killifish process other public dataset.

---

## Repository layout

```
src/                        # Core pipeline modules
├── data_loader.py          # Atlas counts, metadata, query file loading
├── normalization.py        # Frequency and DESeq2 normalization
├── preprocessing.py        # Gene filtering, stratification, outlier detection
├── bayesage2.py            # BayesAge 2.0 clock
├── elastic_net.py          # Elastic Net clock
├── pcr.py                  # Principal Component Regression clock
├── calibration.py          # Apply clocks to query datasets
└── gene_mapping.py         # Map Ensembl IDs → Atlas gene names

run_query_clocks.py         # General-purpose CLI: apply clocks to any count matrix

plot_gene_expression.py     # Per-gene TPM expression across tissues + timecourse

data/                       # Gene mapping tables and reference GTF
raw_RNAseq_process/         # Total RNA-seq: download → QC → align → count matrix
unittests/                  # pytest test suite (8 test files)
Costa_et_al/                # Original KillifishAtlas analysis scripts
data/                       # Atlas TPM/count matrices, metadata, and gene mapping
query_data/                 # Query input files
outputs/                    # Generated prediction and figure outputs
```

---

## Quick start

```bash
conda env create -f environment.yml
conda activate killifish-tx-clock

# 1. Normalize Atlas reference (once)
python src/normalize_reference.py

# 2. Run all three clocks on any count matrix
python run_query_clocks.py --counts query_data/toy.csv

# 3. Plot a gene's expression across tissues
python plot_gene_expression.py actb --log --timecourse
```

See [Setup](setup.md) and [Pipeline Overview](pipeline.md) for full details.
