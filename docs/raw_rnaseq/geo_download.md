# Total RNA-seq — GEO Download Pipeline

**Location:** `raw_RNAseq_process/`

Standalone pipeline for total RNA-seq: GEO/SRA download → fastp QC → STAR alignment → count matrix.
Self-contained; does not depend on any external shell library.

---

## Files

| File | Purpose |
|---|---|
| `TotalRNAseq.def` | Apptainer definition: fastp, STAR 2.7.11b, samtools 1.22.1, pandas |
| `setup_genome.sh` | Download killifish genome + GTF from NCBI; build STAR index |
| `run_rnaseq.sh` | Pipeline script: GEO download → alignment → count matrix |

---

## Container Build

```bash
cd raw_RNAseq_process/
apptainer build TotalRNAseq.sif TotalRNAseq.def
```

Tools included (no bowtie2/rRNA filter — not needed for total RNA-seq):

- fastp (opengene binary)
- STAR 2.7.11b (compiled from source)
- samtools 1.22.1 (compiled from source)
- Python 3 + pandas

---

## Pipeline Steps

| Step | Tool | Input → Output |
|---|---|---|
| 1 | prefetch + fasterq-dump | SRR accession → `tmp/<s>_R1.fastq.gz`, `_R2.fastq.gz` |
| 2 | fastp | paired FASTQ → trimmed FASTQ + QC JSON/HTML |
| 3 | STAR | trimmed FASTQ → sorted BAM + `ReadsPerGene.out.tab` |
| 4 | Python | all `ReadsPerGene.out.tab` → `count_matrix.csv` |

Trimmed FASTQs are deleted after STAR alignment (use `-k` to keep).

---

## Usage

```bash
# Build container (once)
apptainer build TotalRNAseq.sif TotalRNAseq.def

# From GEO series ID (requires entrez-direct)
./run_rnaseq.sh -g GSE123456 -i /path/to/star_genome_index

# From explicit SRR accessions
./run_rnaseq.sh -i /path/to/star_index SRR12345678 SRR12345679

# From TSV mapping file (sample_id <tab> SRR  OR  sample_id <tab> R1.fastq.gz,R2.fastq.gz)
./run_rnaseq.sh -m samples.tsv -i /path/to/star_index

# Override strandedness (default: reverse = dUTP/TruSeq Stranded)
./run_rnaseq.sh -g GSE123456 -i /path/to/star_index -s unstranded
```

**Options:**

| Flag | Default | Description |
|---|---|---|
| `-g GSE_ID` | — | GEO series ID; fetches SRR list via esearch/efetch |
| `-m FILE` | — | TSV: sample_id ↔ SRR or R1,R2 paths |
| `-i DIR` | required | STAR genome index directory |
| `-o DIR` | `./results` | Output directory |
| `-s STRAND` | `reverse` | `reverse` / `forward` / `unstranded` |
| `-t INT` | `8` | Threads |
| `-k` | off | Keep trimmed FASTQs after alignment |

---

## Outputs

```
results/
  count_matrix.csv                              ← gene × sample raw count matrix
  qc_reports/<sample>_fastp.html/.json
  STAR_out_<sample>_Aligned.sortedByCoord.out.bam
  STAR_out_<sample>_ReadsPerGene.out.tab
```

The count matrix uses STAR GeneCounts column selection based on strandedness:

- `reverse` → column 4 (reverse-strand)
- `forward` → column 3 (forward-strand)
- `unstranded` → column 2

---

## Killifish Reference Genome

**Species:** *Nothobranchius furzeri* (African turquoise killifish)  
**Assembly:** GCF_043380555.1 NfurGRZ-RIMD1 (NCBI RefSeq, Dec 2024)

```bash
# One-time setup: download genome + GTF + build STAR index
./setup_genome.sh                        # outputs to raw_RNAseq_process/ref/
./setup_genome.sh -o /data/killifish     # custom location
./setup_genome.sh -t 16                  # more threads for STAR index build

# Then run the pipeline pointing at the built index
./run_rnaseq.sh -i ref/star_index -g GSE123456
```

Downloaded files (in `ref/`):

| File | Source | Size |
|---|---|---|
| `GCF_043380555.1_NfurGRZ-RIMD1_genomic.fna` | NCBI FTP (354 MB gz) | ~1.3 GB |
| `GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf` | NCBI FTP (17 MB gz) | ~65 MB |
| `ref/star_index/` | Built by STAR | ~27 GB |

`genomeSAindexNbases 14` is correct for this genome size (~1.3 GB).
`--skip-download` flag skips wget/gunzip if files are already present.

---

## Prerequisites

- Apptainer installed
- STAR genome index (build separately with `./setup_genome.sh`)
- SRA Toolkit (`prefetch`, `fasterq-dump`) for SRA/GEO input
- entrez-direct (`esearch`, `efetch`) for `-g GSE_ID` mode only
