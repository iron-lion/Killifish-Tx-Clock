# data_loader — `DataLoader`

Loads the KillifishAtlas expression data and query DE result files.

**File:** `src/data_loader.py`  
**Test:** `unittests/test_data_loader.py`

## Atlas dataset

The Atlas contains **677 samples** across **15 tissues** at ages 47–163 days.

| Property | Value |
|---|---|
| Genes | 25,122 |
| Samples | 677 |
| Tissues | Brain, Eye, SpinalCord, Heart, Lung, Liver, Gut, Fat, Kidney, Spleen, Muscle, Skin, Bone, Testis, Ovary |
| Age range | 47–163 days |

## Methods

### `load_counts()`

Loads the raw count matrix from `data/GSE308970_rawcount_Atlas_allbatches_merged_v3.csv`.

Returns a `pd.DataFrame` (genes × samples).

### `load_metadata()`

Loads sample metadata from `data/ExperimentDesign_allbatches_combined_v7.csv`.

Columns: `tissue`, `age_days`, `sex`, `batch`.  
Note: `Ovary` and `Testis` are harmonized to `Gonad` in some internal operations.

### `load_atlas()`

Returns `(counts, metadata)` aligned to shared samples.

### `load_query_files()`

Loads all `query_data/*.xlsx` DE result files. Returns a dict of DataFrames keyed
by tissue + comparison.

### `summary()`

Prints tissue/sex breakdown and age range to stdout.

## Example

```python
from src.data_loader import DataLoader

loader = DataLoader()
counts, meta = loader.load_atlas()
print(counts.shape)   # (25122, 677)
print(meta.dtypes)
```
