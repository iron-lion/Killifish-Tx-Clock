# Module Descriptions

## `data_loader.py` — `DataLoader`
Loads Atlas expression data and query DE result files.

- `load_counts()` — raw count matrix (25,122 genes × 677 samples)
- `load_metadata()` — sample metadata; harmonizes Ovary/Testis → Gonad
- `load_atlas()` — returns `(counts, metadata)` aligned to shared samples
- `load_query_files()` — loads all `query_data/*.xlsx` DE result files
- `summary()` — prints tissue/sex breakdown and age range

---

## `normalization.py` — `FrequencyNormalizer`, `DESeq2Normalizer`
Two normalization strategies required by the different clocks.

- **`FrequencyNormalizer`** — divides each gene by total reads per sample → relative frequency in [0, 1]. Used as input for BayesAge 2.0.
- **`DESeq2Normalizer`** — median-of-ratios size-factor normalization via `pydeseq2`. Used as input for EN and PCR clocks.

---

## `preprocessing.py` — `Preprocessor`
Gene filtering and sample stratification before clock training.

- `filter_genes(counts)` — removes genes with total count ≤ `min_count` (default 1)
- `stratify(counts, metadata, tissue, sex)` — subsets to a tissue (and optionally sex)
- `detect_outliers(counts, n_sd=2.0)` — PCA PC1-based outlier removal; applied before BayesAge2 and PCR (not EN)

---

## `bayesage2.py` — `BayesAge2Clock`
BayesAge 2.0 implementation. Predicts transcriptomic age via Poisson log-likelihood maximization over a LOWESS reference.

- `build_reference(raw_counts, metadata)` — computes Spearman correlation per gene with age; fits LOWESS on top 250 genes; stores reference matrix (genes × age grid)
- `loso_cv(raw_counts, metadata, m_values)` — leave-one-sample-out CV for multiple gene-set sizes M
- `predict(raw_counts, n_genes, reference)` — predicts tAge for query samples using top M genes
- `load_reference(path)` — loads a pre-saved reference `.tsv`

Key details:
- Age range: 47–163 days at 1-day resolution (killifish-specific)
- Spearman correlation is vectorized (rank-correlation via matrix ops) for speed (~0.85 s for 25K genes)
- Poisson log-PMF floor at log(1e-3) matches original BayesAge2 behavior

---

## `elastic_net.py` — `ElasticNetClock`
Elastic Net regularized regression clock.

- `tune_and_train(norm_counts, metadata)` — GridSearchCV over `alpha × l1_ratio` with leave-one-out CV; refits final model at `max_iter=100,000`
- `loso_cv(norm_counts, metadata)` — leave-one-sample-out CV predictions at best hyperparameters
- `predict(norm_counts)` — predicts ages for new samples using the trained model
- `save(out_dir)` — saves best params and non-zero gene coefficients

Key details:
- Input: DESeq2-normalized counts, z-scaled per gene (`StandardScaler`)
- Default `alpha` grid: `[1e-2, 1e-1]`; `l1_ratio` grid: 0.0–1.0 in steps of 0.2
- `top_n_var_genes` optional pre-filter (e.g., 5000) for speed; default=None uses all genes
- Brain tissue uses `max_iter=30,000` during grid search (vs 10,000 for other tissues)

---

## `pcr.py` — `PCRClock`
Principal Component Regression clock using a scikit-learn Pipeline.

- `loso_cv(norm_counts, metadata)` — LOSO-CV across a range of `n_components`; selects optimal by highest R² (ties broken by lowest MAE)
- `fit(norm_counts, metadata, n_components)` — fits final `Pipeline(StandardScaler → PCA → LinearRegression)` on all data
- `predict(norm_counts)` — predicts ages for new samples
- `get_feature_importance()` — gene importance scores as `loadings.T @ coef` (paper eq. 7)
- `save_loadings(out_dir)` — saves per-component loadings, top genes, feature importance, and CV metrics

Key details:
- Default `n_components_range`: `[5, 10, 15, 20]`
- `top_n_var_genes` optional pre-filter (same pattern as EN)

---

## `gene_mapping.py` — `GeneMapper`
Maps query Ensembl gene IDs (ENSNFUG...) to Atlas NCBI gene names.

Three-layer mapping strategy:
1. **Direct** — lowercase `gene_name` → Atlas gene name (7,775 genes)
2. **BioMart fallback** — Ensembl `external_gene_name` for unmatched genes (+393)
3. **LOC107 via NCBI gene2ensembl** — ENSNFUG → GeneID 107XXXXXX → `LOC107XXXXXX` (+4,685)

Total: **12,482 / 23,991 query genes mapped (52%)**

- `convert(counts)` — maps ENSNFUG index to Atlas gene names; drops unmapped genes
- `build_and_save(...)` — rebuilds mapping from scratch (requires internet for BioMart); saves to `data_matrices/query_to_atlas_gene_mapping.csv`

---

## `calibration.py` — `QueryCountExtractor`, `CalibrationManager`
Applies trained Atlas clocks to query datasets.

**`QueryCountExtractor`**
- Parses per-tissue count columns from query xlsx files (DESeq2-normalized values)
- Deduplicates samples appearing across multiple comparison files
- Intersects genes across files for consistency
- Applies `GeneMapper` to convert ENSNFUG IDs to Atlas names
- `correct_batch(query_counts, atlas_raw)` — optional ComBat-seq batch correction via `inmoose.pycombat`

**`CalibrationManager`**
- `run_bayesage2(atlas_raw, atlas_meta, query_counts, query_meta)` — builds reference on Atlas; predicts query samples at M=5..200
- `run_pcr(atlas_norm, atlas_meta, query_norm, query_meta)` — fits Pipeline on Atlas; predicts query; computes Mann-Whitney U (Young vs Old) per n_components
- `run_en(atlas_norm, atlas_meta, query_norm, query_meta)` — tunes/trains EN on Atlas; runs LOSO-CV; predicts query

> **Note:** Verify `YOUNG_AGE_DAYS`, `OLD_AGE_DAYS`, and `CONTROL_CONDITION` constants in `calibration.py` before running against your experiment.

---

## `normalize_reference.py`
Convenience script to normalize the full Atlas dataset and save outputs.

Produces:
- `outputs/normalized/Atlas_freq_normalized.csv` — input for BayesAge 2.0
- `outputs/normalized/Atlas_DESeq2_normalized.csv` — input for EN / PCR

---

## `plot_pcr_query.py`
Scatter-box plots of query age predictions for each tissue × model combination.

- X-axis: conditions in biological order (72h Fasted → 6h Refed → 24h Refed)
- Y-axis: predicted tAge (days), split by Young / Old
- Optimal parameter selection: stable-M heuristic for BayesAge2; lowest Mann-Whitney p for PCR
- Outputs: `outputs/figures/{Tissue}_combined_{Model}_query_<param>.svg`

---
