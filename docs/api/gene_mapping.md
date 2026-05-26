# gene_mapping — `GeneMapper`

Maps query Ensembl gene IDs (`ENSNFUG...`) to Atlas NCBI gene names.

**File:** `src/gene_mapping.py`  
**Test:** `unittests/test_gene_mapping.py`

## When mapping is needed

| Dataset source | Gene ID format | Mapping needed |
|---|---|---|
| Ensembl pipeline | `ENSNFUG00000001234` | **Yes — `GeneMapper.convert()`** |
| `raw_RNAseq_process` output (NCBI GTF) | `actb`, `LOC107XXXXXX` | **No — Atlas names directly** |

`GeneMapper` is only invoked when the input index contains `ENSNFUG` IDs.
Datasets already using Atlas NCBI gene names (symbols + LOC107) skip mapping entirely
and are intersected with the Atlas gene set by name.

## The mapping problem (ENSNFUG datasets)

The KillifishAtlas uses NCBI gene names (e.g. `actb`, `mb21d1`, `LOC107374091`).
Query datasets produced from newer NCBI genome assemblies use Ensembl IDs
(`ENSNFUG00000001234`). Direct string matching fails for ~48% of genes.

## Three-layer mapping strategy

| Layer | Source | Genes mapped |
|---|---|---|
| 1 — Direct | lowercase `gene_name` → Atlas gene name | 7,775 genes |
| 2 — BioMart | Ensembl `external_gene_name` for unmatched genes | +393 genes |
| 3 — LOC107 via NCBI gene2ensembl | ENSNFUG → GeneID 107XXXXXX → `LOC107XXXXXX` | +4,685 genes |

**Total: 12,482 / 23,991 query genes mapped (52%)**

## Methods

### `convert(counts)`

Maps the ENSNFUG row index to Atlas gene names. Drops unmapped genes and
deduplicates any many-to-one mappings.

```python
from src.gene_mapping import GeneMapper

mapper = GeneMapper()
atlas_counts = mapper.convert(ensembl_counts)
```

### `build_and_save(...)`

Rebuilds the mapping table from scratch.
Requires internet access for BioMart queries.
Saves to `data/gene_id_map.csv`.

Not needed for normal use — the pre-built CSV is included in the repo.

## Pre-built mapping files

Located in `data/`:

| File | Description |
|---|---|
| `gene_id_map.csv` | Full unified mapping table (GTF + NCBI gene2ensembl + Atlas names) |
| `query_to_atlas_gene_mapping.csv` | Ensembl ID → Atlas gene name (used by `GeneMapper`) |
| `ncbi_gene2ensembl_nfurzeri.csv` | NCBI gene2ensembl mapping for *N. furzeri* |

See [Gene ID Mapping](../data/gene_mapping.md) for details on how these files were built.
