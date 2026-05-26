# normalization — `FrequencyNormalizer`, `DESeq2Normalizer`

Two normalization strategies, each required by a different subset of clocks.

**File:** `src/normalization.py`  
**Test:** `unittests/test_normalization.py`

## `FrequencyNormalizer`

Divides each gene's count by the total read count for that sample, producing
relative frequencies in the range [0, 1].

```
freq_ij = count_ij / sum_j(count_ij)
```

**Used by:** BayesAge 2.0

**Why raw counts:** BayesAge 2.0 models gene expression as a Poisson process.
Relative frequencies are the natural input; applying DESeq2 first would
double-normalize.

### Usage

```python
from src.normalization import FrequencyNormalizer

norm = FrequencyNormalizer()
freq_counts = norm.normalize(raw_counts)   # DataFrame, same shape
```

---

## `DESeq2Normalizer`

Implements median-of-ratios size-factor normalization via
[`pydeseq2`](https://github.com/owkin/PyDESeq2).

Default design factors: `["age_days"]`

**Used by:** Elastic Net clock, PCR clock

### Usage

```python
from src.normalization import DESeq2Normalizer

norm = DESeq2Normalizer(design_factors=["age_days"])
norm_counts = norm.normalize(raw_counts, metadata)

# Optionally save to disk at the same time
norm_counts = norm.normalize(raw_counts, metadata, save_path=Path("outputs/normalized/Atlas_DESeq2_normalized.csv"))
```

### Caching

The pre-normalized Atlas matrix is cached by `src/normalize_reference.py`:

```bash
python src/normalize_reference.py
# → outputs/normalized/Atlas_DESeq2_normalized.csv
```

All run scripts load this cached file rather than re-normalizing each time.

---

## Normalization choice per clock

| Clock | Input | Normalizer |
|---|---|---|
| BayesAge 2.0 | Raw counts | `FrequencyNormalizer` |
| PCR | DESeq2-normalized | `DESeq2Normalizer` |
| Elastic Net | DESeq2-normalized | `DESeq2Normalizer` |
