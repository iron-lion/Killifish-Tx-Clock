# bayesage2 — `BayesAge2Clock`

BayesAge 2.0 implementation. Predicts transcriptomic age via
Poisson log-likelihood maximization over a LOWESS-smoothed reference.

**File:** `src/bayesage2.py`  
**Test:** `unittests/test_bayesage2.py`

## Algorithm

1. **Spearman correlation** — rank-correlate every gene's frequency-normalized
   expression with `age_days` across Atlas training samples.
   Vectorized: all ~25K genes at once via rank-matrix operations.

2. **LOWESS smoothing** — fit a smoothed age–expression curve (LOWESS, `frac=0.7`)
   for the top `lowess_top_n` genes ranked by |Spearman r|.
   Produces a reference matrix of shape `(genes × age_grid)`.

3. **Prediction** — for a query sample, select the top M genes by |Spearman r|
   and find the age that maximizes the Poisson log-PMF between the query's
   frequency-normalized counts and the reference curves.

The age grid spans **47–163 days** at 1-day resolution (killifish-specific).

## Constructor parameters

| Parameter | Default | Description |
|---|---|---|
| `lowess_frac` | `0.7` | Fraction of points used per local regression |
| `lowess_top_n` | `500` | Gene pool size for LOWESS fitting |

!!! note "LOWESS fraction"
    The original BayesAge 2.0 notebook uses `frac=0.7`. An earlier version of
    this repo used `frac=0.3`, which was the single largest driver of differences
    from the original. Fixed in this implementation.

## Methods

### `build_reference(raw_counts, metadata)`

Builds the age-expression reference matrix on Atlas training data.

Stores internally: Spearman correlations, LOWESS-smoothed reference, age grid.

### `predict(raw_counts, n_genes, reference=None)`

Predicts tAge for all samples in `raw_counts` using the top `n_genes` (M) genes.
Optionally accepts a pre-loaded reference (from `load_reference()`).

Returns a `pd.Series` mapping sample ID → predicted tAge (days).

### `loso_cv(raw_counts, metadata, m_values)`

Leave-one-sample-out cross-validation over multiple M values.
Returns a DataFrame with columns `[sample_id, age_days, tAge_M{m}, ...]`.

### `load_reference(path)`

Loads a pre-saved reference `.tsv` instead of rebuilding from scratch.

### `save_reference(path)`

Saves the reference matrix to a TSV file.

## Example

```python
from src.data_loader import DataLoader
from src.preprocessing import Preprocessor
from src.normalization import FrequencyNormalizer
from src.bayesage2 import BayesAge2Clock

loader = DataLoader()
counts, meta = loader.load_atlas()

pp = Preprocessor()
tissue_counts, tissue_meta = pp.stratify(counts, meta, tissue="Liver")
clean_counts, clean_meta, _ = pp.detect_outliers(tissue_counts)

freq = FrequencyNormalizer().normalize(clean_counts)

clock = BayesAge2Clock()
clock.build_reference(freq, clean_meta)
predictions = clock.predict(freq, n_genes=100)
```

## Performance notes

- Spearman correlation: ~0.85 s for 25K genes (vectorized)
- LOWESS fitting: ~2–5 s for 500 genes
- Prediction: ~0.1 s per sample
