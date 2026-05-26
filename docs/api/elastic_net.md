# elastic_net — `ElasticNetClock`

Elastic Net regularized linear regression clock.

**File:** `src/elastic_net.py`  
**Test:** `unittests/test_elastic_net.py`

## Algorithm

1. **Hyperparameter search** — GridSearchCV with leave-one-out CV over an
   `alpha × l1_ratio` grid. Each fold independently z-scales the training data
   (no leakage from held-out sample).

2. **Final model** — refit `ElasticNet` on all Atlas samples with best
   `(alpha, l1_ratio)` at `max_iter=100,000`.

3. **LOSO-CV** — leave-one-sample-out predictions on the Atlas training set,
   each fold with an independent `StandardScaler`.

4. **Prediction** — z-scale query counts using Atlas statistics; apply trained model.

## Hyperparameter grid

| Parameter | Values |
|---|---|
| `alpha` | `[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100]` |
| `l1_ratio` | `0.0, 0.1, 0.2, …, 1.0` (11 values) |

!!! warning "Computation time"
    GridSearchCV with LOO-CV and 88 parameter combinations is computationally
    intensive (~minutes per tissue on a single CPU).

## Constructor parameters

| Parameter | Default | Description |
|---|---|---|
| `top_n_var_genes` | `None` | Pre-filter to top-N most variable genes before fitting |

## Methods

### `tune_and_train(norm_counts, metadata)`

GridSearchCV over `alpha × l1_ratio`; fits final model on all data.

### `loso_cv(norm_counts, metadata)`

Leave-one-sample-out CV predictions at best hyperparameters.
Returns DataFrame with columns `[sample_id, age_days, predicted_age]`.

### `predict(norm_counts)`

Predicts tAge for new samples. Input must be DESeq2-normalized.

### `save(out_dir)`

Saves best hyperparameters and non-zero gene coefficients to `out_dir`.

## Output columns

`{tissue}_sexcombined_EN_query_loso.csv`:

| Column | Description |
|---|---|
| `sample_id` | Sample identifier |
| `age_days` | True age (Atlas LOSO samples only) |
| `predicted_age` | Predicted tAge |
| `source` | `Atlas` or `Query` |

## Design choices

### No scaler leakage

The original Costa et al. notebook fits `StandardScaler` on the full Atlas before
calling `GridSearchCV`, leaking held-out statistics. This repo fits the scaler
independently per fold inside each CV split.

### `top_n_var_genes`

For Brain tissue, the large number of samples makes GridSearchCV slow.
Setting `top_n_var_genes=5000` provides a practical speed-up with minimal
accuracy loss. The default (`None`) uses all genes.
