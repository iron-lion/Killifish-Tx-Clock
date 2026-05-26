# preprocessing — `Preprocessor`

Gene filtering, sample stratification, and outlier removal before clock training.

**File:** `src/preprocessing.py`  
**Test:** `unittests/test_preprocessing.py`

## Methods

### `filter_genes(counts, min_count=1)`

Removes genes whose total count across all samples is ≤ `min_count`.
Default `min_count=1` removes genes that are completely absent.

```python
from src.preprocessing import Preprocessor

pp = Preprocessor()
filtered = pp.filter_genes(raw_counts)
```

### `stratify(counts, metadata, tissue, sex=None)`

Subsets the count matrix and metadata to a specific tissue, and optionally to
a single sex.

```python
# All samples for Liver, both sexes (sex-combined mode used in this repo)
liver_counts, liver_meta = pp.stratify(counts, meta, tissue="Liver", sex=None)

# Female Liver samples only
liver_f_counts, liver_f_meta = pp.stratify(counts, meta, tissue="Liver", sex="F")
```

### `detect_outliers(counts, n_sd=2.0)`

PCA-based outlier detection. Computes PC1 scores across all samples and flags
any sample whose PC1 score deviates more than `n_sd` standard deviations from
the mean.

Applied before BayesAge 2.0 and PCR. Not applied before Elastic Net.

Returns a **2-tuple** `(counts_clean, outlier_ids)` — metadata is not modified by this method; the caller must filter it separately using the returned `outlier_ids`.

```python
clean_counts, outlier_ids = pp.detect_outliers(counts, n_sd=2.0)
clean_meta = meta.drop(index=outlier_ids)
```

## Order of operations

```
raw_counts
    └─► filter_genes()         # remove zero-count genes
    └─► stratify()             # subset to tissue (and sex)
    └─► detect_outliers()      # PCA-based outlier removal
    └─► normalize()            # FrequencyNormalizer or DESeq2Normalizer
    └─► clock.build_reference() / clock.fit()
```
