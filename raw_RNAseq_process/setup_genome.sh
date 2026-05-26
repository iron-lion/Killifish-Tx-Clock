#!/usr/bin/env bash
# setup_genome.sh
#
# Download and index the Nothobranchius furzeri (killifish) reference genome
# for STAR alignment in the total RNA-seq pipeline.
#
# Source:
#   NCBI RefSeq GCF_043380555.1 NfurGRZ-RIMD1 (Dec 2024)
#   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/043/380/555/GCF_043380555.1_NfurGRZ-RIMD1/
#
# Steps:
#   1. Download genomic FASTA + GTF from NCBI FTP
#   2. Decompress both files
#   3. Build STAR genome index (via TotalRNAseq.sif)
#
# Outputs (under -o DIR, default: ./ref):
#   ref/GCF_043380555.1_NfurGRZ-RIMD1_genomic.fna    (~1.3 GB)
#   ref/GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf
#   ref/star_index/                                   (STAR index, ~27 GB)
#
# Usage:
#   ./setup_genome.sh                       # download + index into ./ref/
#   ./setup_genome.sh -o /data/killifish    # custom output directory
#   ./setup_genome.sh -t 16                 # override thread count
#   ./setup_genome.sh --skip-download       # index only (files already present)
#
# Options:
#   -o DIR           Output directory (default: ./ref)
#   -t INT           Threads for STAR (default: 8)
#   --skip-download  Skip download/decompress; go straight to STAR index build
#   -h               Show this help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER="/ptmp/parky/apptainer_images/TotalRNAseq.sif"

# ── Defaults ──────────────────────────────────────────────────────────────────
REF_DIR="${SCRIPT_DIR}/ref2"
THREADS=8
SKIP_DOWNLOAD=false

BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/043/380/555/GCF_043380555.1_NfurGRZ-RIMD1"
GENOME_GZ="GCF_043380555.1_NfurGRZ-RIMD1_genomic.fna.gz"
GTF_GZ="GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf.gz"
GENOME_FA="GCF_043380555.1_NfurGRZ-RIMD1_genomic.fna"
GTF="GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf"

# ── Logging ───────────────────────────────────────────────────────────────────
log_info()  { echo "[INFO]  $(date '+%H:%M:%S') $*"; }
log_error() { echo "[ERROR] $(date '+%H:%M:%S') $*" >&2; }
log_step()  { echo ""; echo "=== $* ==="; }

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) REF_DIR="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=true; shift ;;
        -h)
            sed -n '2,/^set /p' "$0" | grep '^#' | sed 's/^# \?//'; exit 0 ;;
        *) log_error "Unknown argument: $1"; exit 1 ;;
    esac
done

STAR_INDEX="${REF_DIR}/star_index"
MNT_REF="/mnt/ref"

# ── Container wrapper ─────────────────────────────────────────────────────────
run_container() {
    apptainer exec --bind "${REF_DIR}:${MNT_REF}" "${CONTAINER}" bash -c "$1"
}

# ── Pre-flight ────────────────────────────────────────────────────────────────
[[ -f "${CONTAINER}" ]] || {
    log_error "Container not found: ${CONTAINER}"
    log_error "Build with: apptainer build TotalRNAseq.sif TotalRNAseq.def"
    exit 1; }

mkdir -p "${REF_DIR}"
log_info "Output dir : ${REF_DIR}"
log_info "Threads    : ${THREADS}"

# ── Step 1: Download ──────────────────────────────────────────────────────────
if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    log_step "Downloading genome FASTA and GTF"

    if [[ ! -f "${REF_DIR}/${GENOME_FA}" ]]; then
        if [[ ! -f "${REF_DIR}/${GENOME_GZ}" ]]; then
            log_info "Downloading ${GENOME_GZ} (~354 MB)..."
            wget -c -P "${REF_DIR}" "${BASE_URL}/${GENOME_GZ}"
        else
            log_info "${GENOME_GZ} already downloaded, skipping"
        fi
        log_info "Decompressing ${GENOME_GZ}..."
        gunzip "${REF_DIR}/${GENOME_GZ}"
    else
        log_info "Genome FASTA already present: ${GENOME_FA}"
    fi

    if [[ ! -f "${REF_DIR}/${GTF}" ]]; then
        if [[ ! -f "${REF_DIR}/${GTF_GZ}" ]]; then
            log_info "Downloading ${GTF_GZ} (~17 MB)..."
            wget -c -P "${REF_DIR}" "${BASE_URL}/${GTF_GZ}"
        else
            log_info "${GTF_GZ} already downloaded, skipping"
        fi
        log_info "Decompressing ${GTF_GZ}..."
        gunzip "${REF_DIR}/${GTF_GZ}"
    else
        log_info "GTF already present: ${GTF}"
    fi
else
    log_info "--skip-download: assuming files already in ${REF_DIR}"
fi

# ── Verify reference files exist ──────────────────────────────────────────────
[[ -f "${REF_DIR}/${GENOME_FA}" ]] || {
    log_error "Genome FASTA not found: ${REF_DIR}/${GENOME_FA}"; exit 1; }
[[ -f "${REF_DIR}/${GTF}" ]] || {
    log_error "GTF not found: ${REF_DIR}/${GTF}"; exit 1; }

log_info "Genome : ${REF_DIR}/${GENOME_FA} ($(du -sh "${REF_DIR}/${GENOME_FA}" | cut -f1))"
log_info "GTF    : ${REF_DIR}/${GTF} ($(du -sh "${REF_DIR}/${GTF}" | cut -f1))"

# ── Step 2: Build STAR genome index ──────────────────────────────────────────
log_step "Building STAR genome index"

if [[ -f "${STAR_INDEX}/SA" ]]; then
    log_info "STAR index already exists: ${STAR_INDEX}"
    log_info "Delete ${STAR_INDEX}/SA to force rebuild."
else
    mkdir -p "${STAR_INDEX}"
    log_info "Running STAR --runMode genomeGenerate (this takes ~30-60 min and ~30 GB RAM)..."
    run_container \
        "STAR \
            --runMode            genomeGenerate \
            --runThreadN         ${THREADS} \
            --genomeDir          ${MNT_REF}/star_index \
            --genomeFastaFiles   ${MNT_REF}/${GENOME_FA} \
            --sjdbGTFfile        ${MNT_REF}/${GTF} \
            --genomeSAindexNbases 14"

    log_info "STAR index complete: ${STAR_INDEX}"
fi

# ── Done ──────────────────────────────────────────────────────────────────────
log_step "Done"
echo ""
echo "Genome FASTA : ${REF_DIR}/${GENOME_FA}"
echo "GTF          : ${REF_DIR}/${GTF}"
echo "STAR index   : ${STAR_INDEX}/"
echo ""
echo "Run the RNA-seq pipeline with:"
echo "  ./run_rnaseq.sh -i ${STAR_INDEX} -g <GSE_ID>"
