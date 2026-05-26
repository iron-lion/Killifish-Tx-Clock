# pcr — `PCRClock`

Principal Component Regression clock using a scikit-learn Pipeline.

**File:** `src/pcr.py`  
**Test:** `unittests/test_pcr.py`

## Algorithm

1. **LOSO-CV model selection** — fit `StandardScaler → PCA(n) → LinearRegression`
   for each candidate `n_components`. Select the n with highest LOSO-CV R²
   (ties broken by lowest MAE).

2. **Final model** — refit the selected Pipeline on all Atlas samples.

3. **Prediction** — transform query counts through the fitted Pipeline.

4. **Feature importance** — `gene_importance = loadings.T @ coef`
   (equation 7 from the BayesAge2 paper).

## Constructor parameters

| Parameter | Default | Description |
|---|---|---|
| `n_components_range` | `[5, 10, 15, 20]` | Candidate component counts for LOSO-CV |
| `top_n_var_genes` | `None` | Pre-filter to top-N most variable genes |

## Methods

### `loso_cv(norm_counts, metadata)`

LOSO-CV across all `n_components_range` values using `cross_val_predict` with `LeaveOneOut`.
Selects and stores the optimal n_components (highest R², ties broken by lowest MAE).

Returns a DataFrame with columns `age_days, pred_ncomp_5, pred_ncomp_10, …` indexed by `sample_id`.
CV metrics (MAE, R², Pearson r) per n_components are stored in `self._cv_metrics`.

Note: Mann-Whitney U (Young vs Old) is computed in `CalibrationManager.run_pcr()`, not here.

### `fit(norm_counts, metadata, n_components)`

Fits the final `Pipeline(StandardScaler → PCA → LinearRegression)` on all data.

### `predict(norm_counts)`

Transforms query counts through the fitted Pipeline and returns predicted tAge.

### `get_feature_importance()`

Returns per-gene importance as `loadings.T @ coef`. Higher magnitude = stronger
age-predictive signal.

### `save_loadings(out_dir)`

Saves to `out_dir`:

- per-component gene loadings
- top genes per component
- feature importance scores
- LOSO-CV R² and MAE metrics

## Output columns

**`PCRClock.loso_cv()` returns** (Atlas self-evaluation):

| Column | Description |
|---|---|
| `age_days` | True age |
| `pred_ncomp_5`, `pred_ncomp_10`, … | LOSO-CV predicted tAge per n_components |

**`CalibrationManager.run_pcr()` output** (`{tissue}_sexcombined_PCR_query.csv`):

| Column | Description |
|---|---|
| `sample_id` | Sample identifier |
| `tAge_n5`, `tAge_n10`, … | Predicted tAge at each n_components |
| `age_group` | Young / Old |
| `condition` | Experimental condition |

## Design choice: data-driven model selection

The original notebook selects `n_components` by visual inspection of
Mann-Whitney U results on the query set — involving the query data in model
selection. This repo uses LOSO-CV R² on the Atlas training set only, which
is a proper cross-validation procedure.
