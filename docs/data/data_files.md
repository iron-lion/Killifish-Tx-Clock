# Data Files

## Atlas data matrices (`data/`)

Large files not tracked in git. Obtain from GEO accession **GSE308970**.

| File | Dimensions | Description |
|---|---|---|
| `GSE308970_TPM_Atlas_allbatches_merged_v3.csv` | ~25K genes × 677 samples | TPM-normalized expression matrix |
| `GSE308970_rawcount_Atlas_allbatches_merged_v3.csv` | ~25K genes × 677 samples | Raw count matrix |
| `ExperimentDesign_allbatches_combined_v7.csv` | 677 rows | Sample metadata: `tissue`, `age_days`, `sex`, `batch` |

### Metadata columns

| Column | Values | Notes |
|---|---|---|
| `tissue` | Brain, Eye, SpinalCord, Heart, Lung, Liver, Gut, Fat, Kidney, Spleen, Muscle, Skin, Bone, Testis, Ovary | 15 tissues |
| `age_days` | 47–163 | Killifish lifespan: ~3–6 months |
| `sex` | F, M | Some tissues are single-sex (Testis=M, Ovary=F) |
| `batch` | integer | Sequencing batch number (used by ComBat-seq) |

---

## Gene mapping files

Files tracked in git. Used to map query Ensembl IDs → Atlas gene names.

| File | Description |
|---|---|
| `gene_id_map.csv` | Full unified mapping table (GTF + NCBI + Atlas) |
| `query_to_atlas_gene_mapping.csv` | Ensembl ID → Atlas gene name (direct input for `GeneMapper`) |
| `ncbi_gene2ensembl_nfurzeri.csv` | NCBI gene2ensembl table for *N. furzeri* |
| `GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf` | GTF annotation for NfurGRZ-RIMD1 assembly |
| `build_gene_map.py` | Script that produced `gene_id_map.csv` |

See [Gene ID Mapping](gene_mapping.md) for the full mapping pipeline.

---

## Query data (`query_data/`)

| File | Description |
|---|---|
| `toy.csv` | Minimal toy count matrix for testing `run_query_clocks.py` |

