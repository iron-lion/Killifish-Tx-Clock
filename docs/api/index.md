# Source Modules (`src/`)

The `src/` package contains all reusable pipeline components. Each module
is independently testable and has a corresponding test file under `unittests/`.

## Module summary

| Module | Class(es) | Role |
|---|---|---|
| [`data_loader`](data_loader.md) | `DataLoader` | Load Atlas counts, metadata, and query files |
| [`normalization`](normalization.md) | `FrequencyNormalizer`, `DESeq2Normalizer` | Two normalization strategies |
| [`preprocessing`](preprocessing.md) | `Preprocessor` | Gene filtering, stratification, outlier detection |
| [`bayesage2`](bayesage2.md) | `BayesAge2Clock` | BayesAge 2.0 — Poisson log-likelihood clock |
| [`elastic_net`](elastic_net.md) | `ElasticNetClock` | Elastic Net regularized regression |
| [`pcr`](pcr.md) | `PCRClock` | Principal Component Regression |
| [`gene_mapping`](gene_mapping.md) | `GeneMapper` | Ensembl ENSNFUG IDs → Atlas gene names |
| [`calibration`](calibration.md) | `QueryCountExtractor`, `CalibrationManager` | End-to-end clock application to query data |

## Data flow

```
DataLoader ──► Preprocessor ──► (FrequencyNormalizer) ──► BayesAge2Clock
                     │
                     └──────► (DESeq2Normalizer) ──► PCRClock
                                                  └──► ElasticNetClock

QueryCountExtractor ──► GeneMapper ──► correct_batch ──► CalibrationManager
```

## Tests

```bash
pytest unittests/ -v
```

Eight test files cover all modules:

```
unittests/
  test_data_loader.py
  test_normalization.py
  test_preprocessing.py
  test_bayesage2.py
  test_elastic_net.py
  test_pcr.py
  test_run_clocks.py
```
