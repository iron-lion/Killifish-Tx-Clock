# Raw RNA-seq Processing Pipeline

**Location:** `raw_RNAseq_process/`

Standalone pipeline for total RNA-seq data: GEO/SRA download → fastp QC →
STAR alignment → count matrix. Self-contained; runs inside an Apptainer container.

## Files

| File | Purpose |
|---|---|
| `TotalRNAseq.def` | Apptainer container definition: fastp, STAR 2.7.11b, samtools 1.22.1, Python + pandas |
| `setup_genome.sh` | Download killifish genome + GTF from NCBI; build STAR index |
| `run_rnaseq.sh` | Pipeline: GEO/SRA download → QC → alignment → count matrix |

---

## Quick start

```bash
cd raw_RNAseq_process/

# 1. Build container (once, requires Apptainer)
apptainer build TotalRNAseq.sif TotalRNAseq.def

# 2. Download genome + build STAR index (once, ~30–60 min, ~30 GB RAM)
./setup_genome.sh -o /path/to/ref

# 3. Run pipeline from a GEO accession
./run_rnaseq.sh -g GSE123456 -i /path/to/ref/star_index

# 4. Or from explicit SRR accessions
./run_rnaseq.sh -i /path/to/ref/star_index SRR12345678 SRR12345679

# 5. Or from a TSV mapping file
./run_rnaseq.sh -m samples.tsv -i /path/to/ref/star_index
```

---

## Pipeline steps

| Step | Tool | Input → Output |
|---|---|---|
| 1 | `prefetch` + `fasterq-dump` | SRR accession → `tmp/<s>_R1.fastq.gz`, `_R2.fastq.gz` |
| 2 | fastp | paired/single FASTQ → trimmed FASTQ + QC JSON/HTML |
| 3 | STAR | trimmed FASTQ → sorted BAM + `ReadsPerGene.out.tab` |
| 3b | samtools | BAM → `.bam.bai` index |
| 4 | Python | all `ReadsPerGene.out.tab` → `count_matrix.csv` |

Trimmed FASTQs are deleted after STAR alignment (use `-k` to keep).
Steps are idempotent — completed samples are skipped on re-runs.

---

## `run_rnaseq.sh` options

| Flag | Default | Description |
|---|---|---|
| `-g GSE_ID` | — | GEO series ID; fetches SRR list via `esearch`/`efetch` |
| `-m FILE` | — | TSV: `sample_id <tab> SRR` or `sample_id <tab> R1.fastq.gz,R2.fastq.gz` |
| `-i DIR` | required | STAR genome index directory (must contain `SA` file) |
| `-o DIR` | `./results` | Output directory |
| `-s STRAND` | `reverse` | `reverse` / `forward` / `unstranded` |
| `-e END` | `paired` | `paired` / `single` |
| `-t INT` | `8` | Threads |
| `-k` | off | Keep trimmed FASTQs after alignment |
| `--skip-merge` | off | Process samples but skip final count matrix merge (for SLURM arrays) |
| `--merge-only` | off | Skip sample processing; only merge existing tab files |

---

## Strandedness

| Library type | Flag | STAR column used |
|---|---|---|
| dUTP/TruSeq Stranded | `-s reverse` (default) | column 4 |
| Forward-stranded | `-s forward` | column 3 |
| Unstranded | `-s unstranded` | column 2 |

---

## Outputs

```
results/
  count_matrix.csv                            ← genes × samples raw count matrix
  qc_reports/<sample>_fastp.html/.json
  STAR_out_<sample>_Aligned.sortedByCoord.out.bam
  STAR_out_<sample>_Aligned.sortedByCoord.out.bam.bai
  STAR_out_<sample>_ReadsPerGene.out.tab
```

The `count_matrix.csv` rows are gene IDs from the GTF; summary rows prefixed `N_`
(unmapped, multimapping, noFeature, ambiguous) are excluded automatically.

---

## Setting up the genome reference (`setup_genome.sh`)

Downloads and indexes the *N. furzeri* GRZ reference genome for STAR alignment.

```bash
./setup_genome.sh                         # download + index into ./ref/
./setup_genome.sh -o /data/killifish      # custom output directory
./setup_genome.sh -t 16                   # more threads for STAR index build
./setup_genome.sh --skip-download         # index only (files already present)
```

| Flag | Default | Description |
|---|---|---|
| `-o DIR` | `./ref` | Output directory |
| `-t INT` | `8` | Threads for STAR `--runMode genomeGenerate` |
| `--skip-download` | off | Skip wget/gunzip; go straight to index build |

**Reference genome details:**

| Property | Value |
|---|---|
| Species | *Nothobranchius furzeri* (GRZ strain) |
| Assembly | GCF_043380555.1 NfurGRZ-RIMD1 |
| Source | NCBI RefSeq FTP, December 2024 |
| FASTA size | ~1.3 GB |
| STAR index size | ~27 GB |
| `genomeSAindexNbases` | 14 (correct for genome > 1 Gb) |

---

## Prerequisites

- Apptainer installed (`apptainer build` for container)
- SRA Toolkit (`prefetch`, `fasterq-dump`) for SRA/GEO input
- entrez-direct (`esearch`, `efetch`) for `-g GSE_ID` mode only
