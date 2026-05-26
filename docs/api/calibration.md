# calibration — `QueryCountExtractor`, `CalibrationManager`

End-to-end application of trained Atlas clocks to query datasets.

**File:** `src/calibration.py`  
---

## `QueryCountExtractor`

Parses and prepares query count data from AAlab-style xlsx DE result files.

### Constructor

```python
QueryCountExtractor(query_dir=QUERY_DIR, mapper=None)
```

`mapper` defaults to a freshly constructed `GeneMapper()`.

### `extract_tissue(tissue, young_age_days=56, old_age_days=126)`

Loads all `query_data/*.xlsx` files matching the given tissue name and returns
`(counts, metadata)`:

1. Load xlsx files whose filename contains the tissue name.
2. Deduplicate samples appearing across multiple comparison files.
3. Intersect gene sets across files for consistency.
4. Apply `GeneMapper` to convert `ENSNFUG` IDs → Atlas gene names.

`counts` values are DESeq2-normalized decimals. Use as-is for EN/PCR; round to int for BayesAge 2.0.

### `correct_batch(query_counts, atlas_raw)`

Applies ComBat-seq batch correction (`inmoose.pycombat_seq`) to the
concatenated Atlas + query count matrix, with Atlas as the reference batch.

Returns batch-corrected query counts (float) aligned to the Atlas gene set.

```python
from src.calibration import QueryCountExtractor

extractor = QueryCountExtractor()
query_counts, query_meta = extractor.extract_tissue("Liver")
corrected = extractor.correct_batch(query_counts, atlas_raw)
```

---

## `CalibrationManager`

Orchestrates training and prediction for all three clocks. Has no constructor
arguments — instantiate directly and call the run methods.

Each method returns results as DataFrames; saving to disk is done by the caller
(run scripts).

### `run_bayesage2(atlas_raw, atlas_meta, query_counts, query_meta, m_values=None, ref_save_path=None)`

1. Intersects genes between Atlas and query.
2. Builds `BayesAge2Clock` reference on Atlas raw counts (frequency normalization happens inside the clock; `lowess_top_n=250`).
3. Predicts tAge for query at M = 5, 10, …, 200 (step 5).
4. Returns `(result_df, feature_importance_df)`.

`result_df` columns: `tAge_M5`, `tAge_M10`, …, `age_group`, `condition`.  
`feature_importance_df`: genes ranked by `|spearman_r|` with LOWESS fits.

### `run_pcr(atlas_norm, atlas_meta, query_norm, query_meta, n_components_range=None, top_n_var_genes=None)`

1. Intersects genes; optionally pre-filters to top-N variable genes.
2. For each n in `n_components_range` (default `[5, 10, 15, 20]`): fits `Pipeline(StandardScaler → PCA → LinearRegression)` on Atlas; predicts query.
3. Computes Mann-Whitney U (Young vs Old) per n_components on query predictions.
4. Returns `(result_df, mw_pvals_dict, gene_importance_dict)`.

`result_df` columns: `tAge_n5`, `tAge_n10`, …, `age_group`, `condition`.

### `run_en(atlas_norm, atlas_meta, query_norm, query_meta, tissue="", top_n_var_genes=None)`

1. Intersects genes between Atlas and query.
2. Calls `ElasticNetClock.tune_and_train()` (GridSearchCV + LOO-CV) on Atlas.
3. Runs `ElasticNetClock.loso_cv()` on Atlas.
4. Predicts query samples using the trained model.
5. Returns `(result_df, feature_importance_df)`.

`result_df` columns: `age_days`, `predicted_age`, `source` (`Atlas` / `Query`), `age_group`, `condition`.

### Usage

```python
from src.calibration import CalibrationManager, QueryCountExtractor

extractor = QueryCountExtractor()
query_counts, query_meta = extractor.extract_tissue("Liver")
corrected = extractor.correct_batch(query_counts, atlas_raw)

mgr = CalibrationManager()
result, fi = mgr.run_bayesage2(atlas_raw, atlas_meta, corrected, query_meta)
result, mw, gi = mgr.run_pcr(atlas_norm, atlas_meta, corrected.astype(float), query_meta)
result, fi = mgr.run_en(atlas_norm, atlas_meta, corrected.astype(float), query_meta, tissue="Liver")
```
