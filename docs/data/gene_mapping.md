# Reference genome

**Species:** *Nothobranchius furzeri* (African turquoise killifish)  
**Strain:** GRZ  
**Assembly:** GCF_043380555.1 NfurGRZ-RIMD1  
**Source:** NCBI RefSeq, December 2024  
**Genome size:** ~1.3 Gb  
**GTF file:** `data/GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf`

The GTF is also used by `raw_RNAseq_process/run_rnaseq.sh` to build the STAR index
and quantify gene counts. See [Raw RNA-seq Pipeline](../raw_rnaseq/pipeline.md).

# Gene ID Mapping

The KillifishAtlas uses NCBI-style gene names (e.g. `actb`, `mb21d1`, `LOC107374091`).
Query datasets produced from the NfurGRZ-RIMD1 genome assembly (NCBI RefSeq GCF_043380555.1)
use Ensembl IDs (`ENSNFUG00000001234`). Bridging these two namespaces requires a
multi-source mapping pipeline.

## Building the unified map (`data/build_gene_map.py`)

**Script:** `data/build_gene_map.py`

Merges three sources into `data/gene_id_map.csv`:

```
Source 1: GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf
          → gtf_gene_id, gtf_gene_name, ncbi_gene_id

Source 2: ncbi_gene2ensembl_nfurzeri.csv
          → ncbi_gene_id ↔ ensembl_gene_id

Source 3: query_to_atlas_gene_mapping.csv
          → ensembl_gene_id → atlas_gene
```

Merge strategy:

1. Parse `gene` records from the GTF, extracting `gene_id`, `gene` name, and `GeneID` db_xref.
2. Join with `gene2ensembl` on `ncbi_gene_id` to add `ensembl_gene_id`.
3. Join with the Atlas mapping on `ensembl_gene_id` to add `atlas_gene`.
4. Supplement with Ensembl IDs in `gene2ensembl` not reachable via the GTF.
5. Supplement with Atlas entries not reachable via either route.

```bash
cd data/
python build_gene_map.py
# → gene_id_map.csv
```

## Output columns

| Column | Description |
|---|---|
| `gtf_gene_id` | RefSeq gene ID from GTF (e.g. `gene1234`) |
| `gtf_gene_name` | Gene symbol from GTF (e.g. `actb`) |
| `ncbi_gene_id` | NCBI GeneID integer |
| `ensembl_gene_id` | Ensembl ID (`ENSNFUG...`) |
| `atlas_gene` | Atlas gene name used as row index in count matrices |

## Three-layer mapping in `GeneMapper`

`src/gene_mapping.py` applies three layers in order:

| Layer | Method | 
|---|---|
| 1 | Direct lowercase `gene_name` → Atlas |
| 2 | BioMart `external_gene_name` fallback |
| 3 | ENSNFUG → GeneID 107XXXXXX → `LOC107XXXXXX` |

### Coverage against the Atlas gene symbols -> ENSNFUGxxx

Measured by `unittests/test_gene_mapping.py` against `GSE308970_Counts_Atlas_allbatches_merged_v3.csv`

If a query experiment measured all Atlas gene symbols by their ENSNFUG IDs, how many could GeneMapper successfully translate?
(25,122 Atlas genes × 677 samples; run 2026-05-26) :

| Gene type | Atlas total | Covered by GeneMapper | Coverage |
|---|---|---|---|
| Named symbols (e.g. `actb`) | 10,533 | 7,706 | **73.2 %** |
| LOC genes (e.g. `LOC107374091`) | 14,589 | 4,311 | **29.5 %** |
| **All genes** | **25,122** | **12,017** | **47.8 %** |

The lower LOC-gene coverage reflects a biological reality: most `LOC107XXXXXX` entries
in the Atlas are unannotated loci that lack Ensembl cross-references in either BioMart or
the NCBI gene2ensembl table. The 13,105 uncovered genes are dropped before clock training
and have negligible impact on clock accuracy because unannotated loci carry little
age-predictive signal.

### Coverage: raw_RNAseq_process output (PRJNA817434)

Data produced by `raw_RNAseq_process/run_rnaseq.sh` uses NCBI-style gene names directly
from the GTF — **no ENSNFUG IDs, so `GeneMapper` does not apply**. Coverage is measured
as a direct set intersection with Atlas gene names.

Measured against `raw_RNAseq_process/results/PRJNA817434/PRJNA817434_raw_count.csv`
(36,530 genes × 9 samples; run 2026-05-26):

| Gene type | PRJNA817434 total | In Atlas | Coverage |
|---|---|---|---|
| Named symbols (e.g. `acsl4a`) | 17,657 | 7,098 | **40.2 %** |
| LOC107 genes (e.g. `LOC107374091`) | 6,203 | 5,103 | **82.3 %** |
| LOC139 genes (e.g. `LOC139071432`) | 7,031 | 0 | **0.0 %** |
| tRNA / KEG92 entries | 5,639 | — | excluded |
| **Usable total** | **36,530** | **12,249** | **33.5 %** |

`LOC139XXXXXX` genes come from a newer NCBI annotation (GeneID range 139M) that does not
exist in the Atlas (which uses the older 107M range). They are unmappable without rebuilding
the Atlas with the updated GTF.
