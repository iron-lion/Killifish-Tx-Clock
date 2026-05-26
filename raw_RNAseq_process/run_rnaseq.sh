#!/usr/bin/env bash
# run_rnaseq.sh
#
# Total RNA-seq pipeline: GEO/SRA download → fastp QC → STAR align → count matrix
#
# Pipeline steps per sample:
#   1. [GEO/SRA] prefetch + fasterq-dump  SRR accession  → tmp/<s>_R1.fastq.gz, _R2.fastq.gz
#   2. fastp     paired-end adapter trim + QC             → tmp/<s>_R1_trimmed.fastq.gz, _R2_trimmed.fastq.gz
#   3. STAR      genome alignment + GeneCounts            → results/STAR_out_<s>_*.bam + ReadsPerGene.out.tab
#      [delete trimmed FASTQs]
# Final:
#   4. Merge ReadsPerGene.out.tab files                   → results/count_matrix.csv
#
# Input modes (choose one):
#   -g GSE_ID    GEO series ID (requires entrez-direct: esearch/efetch)
#   -m FILE      TSV: sample_id <tab> SRR_accession (or paired R1,R2 FASTQ paths)
#   positional   bare SRR/ERR/DRR accessions
#
# Options:
#   -g GSE_ID      GEO series accession (auto-fetches SRR list via E-utilities)
#   -m FILE        Sample → SRR/FASTQ mapping TSV
#   -o DIR         Output directory (default: ./results)
#   -i DIR         STAR genome index directory (default: ./genome_index)
#   -s STRAND      Strandedness: reverse|forward|unstranded (default: reverse)
#   -e END         Library layout: paired|single (default: paired)
#   -t INT         Threads (default: 8)
#   -k             Keep trimmed FASTQs after alignment
#   --skip-merge   Process samples but skip final count matrix merge (SLURM array tasks)
#   --merge-only   Skip sample processing; only merge existing ReadsPerGene.out.tab files
#   -h             Show this help
#
# Column mapping for count matrix (STAR GeneCounts):
#   reverse     → column 4 (reverse-strand counts)
#   forward     → column 3 (forward-strand counts)
#   unstranded  → column 2 (unstranded counts)
#
# Examples:
#   # From GEO series (fetches SRR list automatically)
#   ./run_rnaseq.sh -g GSE123456 -i /path/to/star_index
#
#   # From explicit SRR accessions
#   ./run_rnaseq.sh -i /path/to/star_index SRR12345678 SRR12345679
#
#   # From mapping TSV (sample_id <tab> SRR or R1.fastq.gz,R2.fastq.gz)
#   ./run_rnaseq.sh -m samples.tsv -i /path/to/star_index --strand unstranded
#
# Prerequisites:
#   - RNA_seq_process/TotalRNAseq.sif   (build: apptainer build TotalRNAseq.sif TotalRNAseq.def)
#   - STAR genome index (-i flag)
#   - SRA Toolkit (prefetch, fasterq-dump) for SRA/GEO input
#   - entrez-direct (esearch, efetch) for -g GSE_ID mode

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER="/ptmp/parky/apptainer_images/TotalRNAseq.sif"

# ── Defaults ──────────────────────────────────────────────────────────────────
GSE_ID=""
MAP_FILE=""
OUT_DIR="${SCRIPT_DIR}/results"
GENOME_INDEX=""
STRAND="reverse"
ENDEDNESS="paired"
THREADS=8
KEEP_TMP=false
SKIP_MERGE=false
MERGE_ONLY=false

# ── Logging ───────────────────────────────────────────────────────────────────
log_info()  { echo "[INFO]  $(date '+%H:%M:%S') $*"; }
log_error() { echo "[ERROR] $(date '+%H:%M:%S') $*" >&2; }
log_step()  { echo ""; echo "=== $* ==="; }

# ── Parse options ─────────────────────────────────────────────────────────────
# getopts handles short flags; long flags (--skip-merge, --merge-only) handled below
_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-merge)  SKIP_MERGE=true;  shift ;;
        --merge-only)  MERGE_ONLY=true;  shift ;;
        *) _ARGS+=("$1"); shift ;;
    esac
done
set -- "${_ARGS[@]}"

while getopts ":g:m:o:i:s:e:t:kh" opt; do
    case $opt in
        g) GSE_ID="$OPTARG" ;;
        m) MAP_FILE="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        i) GENOME_INDEX="$OPTARG" ;;
        s) STRAND="$OPTARG" ;;
        e) ENDEDNESS="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        k) KEEP_TMP=true ;;
        h) sed -n '2,/^set /p' "$0" | grep '^#' | sed 's/^# \?//'; exit 0 ;;
        :) log_error "-$OPTARG requires an argument"; exit 1 ;;
        \?) log_error "Unknown option: -$OPTARG"; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

[[ "${STRAND}" =~ ^(reverse|forward|unstranded)$ ]] || {
    log_error "-s must be reverse, forward, or unstranded"; exit 1; }
[[ "${ENDEDNESS}" =~ ^(paired|single)$ ]] || {
    log_error "-e must be paired or single"; exit 1; }

# ── Strand → STAR column index ────────────────────────────────────────────────
case "${STRAND}" in
    reverse)    STRAND_COL=4 ;;
    forward)    STRAND_COL=3 ;;
    unstranded) STRAND_COL=2 ;;
esac

# ── Derived paths ─────────────────────────────────────────────────────────────
TMP_DIR="${OUT_DIR}/tmp"
QC_DIR="${OUT_DIR}/qc_reports"
SRA_CACHE="${OUT_DIR}/sra_cache"
MNT_OUT="/mnt/results"
MNT_INDEX="/mnt/index"

# ── Container wrapper ─────────────────────────────────────────────────────────
run_container() {
    local extra_binds="$1"; shift
    local cmd="$*"
    local binds="${OUT_DIR}:${MNT_OUT}"
    [[ -n "${GENOME_INDEX}" ]] && binds="${binds},${GENOME_INDEX}:${MNT_INDEX}"
    [[ -n "${extra_binds}" ]]  && binds="${binds},${extra_binds}"
    apptainer exec --bind "${binds}" "${CONTAINER}" bash -c "${cmd}"
}

# ── Build sample map from all input modes ─────────────────────────────────────
# SAMPLE_MAP[sample_id] = "SRR_acc"  or  "R1.fastq.gz,R2.fastq.gz"
declare -A SAMPLE_MAP

# Mode 1: GSE ID (requires entrez-direct)
if [[ -n "${GSE_ID}" ]]; then
    log_step "Fetching SRR list for ${GSE_ID}"
    command -v esearch &>/dev/null || { log_error "esearch not found — install entrez-direct"; exit 1; }
    while IFS= read -r acc; do
        [[ -z "${acc}" || "${acc}" == "Run" ]] && continue
        SAMPLE_MAP["${acc}"]="${acc}"
    done < <(esearch -db sra -query "${GSE_ID}" | efetch -format runinfo | cut -d',' -f1 | tail -n +2)
    log_info "Found ${#SAMPLE_MAP[@]} SRR accessions for ${GSE_ID}"
fi

# Mode 2: TSV mapping file (sample_id <tab> SRR or R1,R2)
if [[ -n "${MAP_FILE}" ]]; then
    [[ -f "${MAP_FILE}" ]] || { log_error "Map file not found: ${MAP_FILE}"; exit 1; }
    while IFS=$'\t' read -r sid src || [[ -n "$sid" ]]; do
        sid="${sid%%#*}"; sid="${sid//[$'\t\r\n ']/}"
        src="${src%%#*}"; src="${src//[$'\r\n']/}"
        [[ -z "${sid}" ]] && continue
        SAMPLE_MAP["${sid}"]="${src}"
    done < "${MAP_FILE}"
fi

# Mode 3: positional SRR accessions
for arg in "$@"; do
    if [[ "${arg}" =~ ^[SED]RR[0-9]+ ]]; then
        SAMPLE_MAP["${arg}"]="${arg}"
    else
        log_error "Unrecognized argument: ${arg} (expected SRR/ERR/DRR accession)"; exit 1
    fi
done

[[ ${#SAMPLE_MAP[@]} -gt 0 ]] || {
    log_error "No samples specified. Use -g GSE_ID, -m FILE, or positional SRR accessions."
    log_error "Run with -h for usage."
    exit 1; }

# ── Pre-flight ────────────────────────────────────────────────────────────────
preflight() {
    log_step "Pre-flight checks"
    [[ -f "${CONTAINER}" ]] || {
        log_error "Container not found: ${CONTAINER}"
        log_error "Build with: apptainer build TotalRNAseq.sif TotalRNAseq.def"
        exit 1; }
    # Genome index not needed for merge-only mode
    if [[ "${MERGE_ONLY}" == false ]]; then
        [[ -n "${GENOME_INDEX}" && -f "${GENOME_INDEX}/SA" ]] || {
            log_error "STAR genome index not found. Provide with -i DIR (must contain SA file)."
            log_error "Build with STAR --runMode genomeGenerate ..."
            exit 1; }
    fi
    mkdir -p "${OUT_DIR}" "${TMP_DIR}" "${QC_DIR}" "${SRA_CACHE}"
    log_info "Output dir   : ${OUT_DIR}"
    log_info "Genome index : ${GENOME_INDEX}"
    log_info "Strandedness : ${STRAND} (column ${STRAND_COL})"
    log_info "Endedness    : ${ENDEDNESS}"
    log_info "Threads      : ${THREADS}"
    log_info "Samples      : ${!SAMPLE_MAP[*]}"
}

# ── Step 1: SRA download ──────────────────────────────────────────────────────
step_sra_download() {
    local sample="$1" acc="$2"
    local r1="${TMP_DIR}/${sample}_R1.fastq.gz"
    local r2="${TMP_DIR}/${sample}_R2.fastq.gz"
    if [[ -f "${r1}" && -f "${r2}" ]]; then
        log_info "[${sample}] step1 SRA: FASTQs exist, skipping"
        return 0
    fi
    command -v prefetch     &>/dev/null || { log_error "prefetch not found — install SRA Toolkit"; return 1; }
    command -v fasterq-dump &>/dev/null || { log_error "fasterq-dump not found — install SRA Toolkit"; return 1; }

    local sra_file="${SRA_CACHE}/${acc}/${acc}.sra"
    if [[ ! -f "${sra_file}" ]]; then
        log_info "[${sample}] step1a prefetch ${acc}..."
        prefetch --progress --max-size 50G --output-directory "${SRA_CACHE}" "${acc}" \
            || { log_error "[${sample}] prefetch failed for ${acc}"; return 1; }
    else
        log_info "[${sample}] step1a .sra cached: ${sra_file}"
    fi

    # Use a separate temp dir — fasterq-dump conflicts when --temp == --outdir
    local FQ_TMP="${TMP_DIR}/fqtmp_${acc}"
    mkdir -p "${FQ_TMP}"

    log_info "[${sample}] step1b fasterq-dump → ${TMP_DIR}..."
    if [[ "${ENDEDNESS}" == "paired" ]]; then
        fasterq-dump \
            --split-files \
            --outdir  "${TMP_DIR}" \
            --temp    "${FQ_TMP}" \
            --threads "${THREADS}" \
            "${sra_file}" \
            || { log_error "[${sample}] fasterq-dump failed for ${acc}"; rm -rf "${FQ_TMP}"; return 1; }
        rm -rf "${FQ_TMP}"

        local fq1="${TMP_DIR}/${acc}_1.fastq"
        local fq2="${TMP_DIR}/${acc}_2.fastq"
        [[ -f "${fq1}" ]] || { log_error "[${sample}] ${acc}_1.fastq missing — is this paired-end? Use -e single if not"; return 1; }
        [[ -f "${fq2}" ]] || { log_error "[${sample}] ${acc}_2.fastq missing — is this paired-end? Use -e single if not"; return 1; }

        gzip -f "${fq1}" \
            || { log_error "[${sample}] gzip failed on ${fq1}"; return 1; }
        mv "${TMP_DIR}/${acc}_1.fastq.gz" "${r1}" \
            || { log_error "[${sample}] mv failed: ${acc}_1.fastq.gz → ${r1}"; return 1; }

        gzip -f "${fq2}" \
            || { log_error "[${sample}] gzip failed on ${fq2}"; return 1; }
        mv "${TMP_DIR}/${acc}_2.fastq.gz" "${r2}" \
            || { log_error "[${sample}] mv failed: ${acc}_2.fastq.gz → ${r2}"; return 1; }

        log_info "[${sample}] step1 SRA: done → $(basename "${r1}"), $(basename "${r2}")"
    else
        fasterq-dump \
            --outdir  "${TMP_DIR}" \
            --temp    "${FQ_TMP}" \
            --threads "${THREADS}" \
            "${sra_file}" \
            || { log_error "[${sample}] fasterq-dump failed for ${acc}"; rm -rf "${FQ_TMP}"; return 1; }
        rm -rf "${FQ_TMP}"

        local fq1="${TMP_DIR}/${acc}.fastq"
        [[ -f "${fq1}" ]] || { log_error "[${sample}] ${acc}.fastq missing after fasterq-dump"; return 1; }

        gzip -f "${fq1}" \
            || { log_error "[${sample}] gzip failed on ${fq1}"; return 1; }
        mv "${TMP_DIR}/${acc}.fastq.gz" "${r1}" \
            || { log_error "[${sample}] mv failed: ${acc}.fastq.gz → ${r1}"; return 1; }

        log_info "[${sample}] step1 SRA: done → $(basename "${r1}")"
    fi
}

# ── Step 2: fastp — QC + trim (paired or single) ──────────────────────────────
step_fastp() {
    local sample="$1" r1="$2" r2="$3"
    local out_r1="${TMP_DIR}/${sample}_R1_trimmed.fastq.gz"
    [[ -f "${r1}" ]] || { log_error "[${sample}] R1 not found: ${r1}"; return 1; }

    if [[ "${ENDEDNESS}" == "paired" ]]; then
        local out_r2="${TMP_DIR}/${sample}_R2_trimmed.fastq.gz"
        if [[ -f "${out_r1}" && -f "${out_r2}" ]]; then
            log_info "[${sample}] step2 fastp: trimmed FASTQs exist, skipping"
            return 0
        fi
        [[ -f "${r2}" ]] || { log_error "[${sample}] R2 not found: ${r2}"; return 1; }
        log_info "[${sample}] step2 fastp (paired-end, length_required=30)..."
        run_container "" \
            "fastp \
                --in1  ${MNT_OUT}/tmp/${sample}_R1.fastq.gz \
                --in2  ${MNT_OUT}/tmp/${sample}_R2.fastq.gz \
                --out1 ${MNT_OUT}/tmp/${sample}_R1_trimmed.fastq.gz \
                --out2 ${MNT_OUT}/tmp/${sample}_R2_trimmed.fastq.gz \
                --json ${MNT_OUT}/qc_reports/${sample}_fastp.json \
                --html ${MNT_OUT}/qc_reports/${sample}_fastp.html \
                --detect_adapter_for_pe \
                --correction \
                --length_required 30 \
                --thread ${THREADS}"
    else
        if [[ -f "${out_r1}" ]]; then
            log_info "[${sample}] step2 fastp: trimmed FASTQ exists, skipping"
            return 0
        fi
        log_info "[${sample}] step2 fastp (single-end, length_required=30)..."
        run_container "" \
            "fastp \
                --in1  ${MNT_OUT}/tmp/${sample}_R1.fastq.gz \
                --out1 ${MNT_OUT}/tmp/${sample}_R1_trimmed.fastq.gz \
                --json ${MNT_OUT}/qc_reports/${sample}_fastp.json \
                --html ${MNT_OUT}/qc_reports/${sample}_fastp.html \
                --length_required 30 \
                --thread ${THREADS}"
    fi

    log_info "[${sample}] step2 fastp: done"
}

# ── Step 3: STAR — paired-end alignment + GeneCounts ─────────────────────────
step_star_align() {
    local sample="$1"
    local out_bam="${OUT_DIR}/STAR_out_${sample}_Aligned.sortedByCoord.out.bam"
    if [[ -f "${out_bam}.bai" ]]; then
        log_info "[${sample}] step3 STAR: BAM exists, skipping"
        return 0
    fi
    [[ -f "${TMP_DIR}/${sample}_R1_trimmed.fastq.gz" ]] || {
        log_error "[${sample}] trimmed R1 not found (run step2 first)"; return 1; }

    log_info "[${sample}] step3 STAR alignment (${ENDEDNESS}-end)..."
    local reads_in="${MNT_OUT}/tmp/${sample}_R1_trimmed.fastq.gz"
    [[ "${ENDEDNESS}" == "paired" ]] && \
        reads_in="${reads_in} ${MNT_OUT}/tmp/${sample}_R2_trimmed.fastq.gz"
    run_container "" \
        "STAR \
            --runThreadN         ${THREADS} \
            --genomeDir          ${MNT_INDEX} \
            --readFilesIn        ${reads_in} \
            --readFilesCommand   zcat \
            --outFileNamePrefix  ${MNT_OUT}/STAR_out_${sample}_ \
            --outSAMtype         BAM SortedByCoordinate \
            --outSAMattributes   NH HI AS NM \
            --quantMode          GeneCounts && \
         samtools index ${MNT_OUT}/STAR_out_${sample}_Aligned.sortedByCoord.out.bam"

    local n; n=$(run_container "" \
        "samtools view -c -F 4 ${MNT_OUT}/STAR_out_${sample}_Aligned.sortedByCoord.out.bam")
    log_info "[${sample}] step3 STAR: ${n} aligned reads"
}

# ── Cleanup tmp FASTQs ────────────────────────────────────────────────────────
cleanup_tmp() {
    local sample="$1"
    [[ "${KEEP_TMP}" == true ]] && return 0
    rm -f \
        "${TMP_DIR}/${sample}_R1.fastq.gz" \
        "${TMP_DIR}/${sample}_R1_trimmed.fastq.gz"
    if [[ "${ENDEDNESS}" == "paired" ]]; then
        rm -f \
            "${TMP_DIR}/${sample}_R2.fastq.gz" \
            "${TMP_DIR}/${sample}_R2_trimmed.fastq.gz"
    fi
    log_info "[${sample}] tmp FASTQs removed"
}

# ── Per-sample orchestration ──────────────────────────────────────────────────
run_sample() {
    local sample="$1" src="$2"
    log_step "Sample: ${sample}"

    local out_bam="${OUT_DIR}/STAR_out_${sample}_Aligned.sortedByCoord.out.bam"
    if [[ -f "${out_bam}.bai" ]]; then
        log_info "[${sample}] already complete — skipping"
        return 0
    fi

    # Resolve R1/R2 paths
    local r1 r2=""
    if [[ "${src}" =~ ^[SED]RR[0-9]+ ]]; then
        step_sra_download "${sample}" "${src}" || return 1
        r1="${TMP_DIR}/${sample}_R1.fastq.gz"
        [[ "${ENDEDNESS}" == "paired" ]] && r2="${TMP_DIR}/${sample}_R2.fastq.gz"
    elif [[ "${src}" == *","* ]]; then
        # TSV row with explicit R1,R2 paths (paired only)
        [[ "${ENDEDNESS}" == "paired" ]] || {
            log_error "[${sample}] comma-separated paths require -e paired"; return 1; }
        r1="${src%%,*}"; r2="${src##*,}"
        [[ -f "${r1}" ]] || { log_error "[${sample}] R1 not found: ${r1}"; return 1; }
        [[ -f "${r2}" ]] || { log_error "[${sample}] R2 not found: ${r2}"; return 1; }
        ln -sf "$(realpath "${r1}")" "${TMP_DIR}/${sample}_R1.fastq.gz"
        ln -sf "$(realpath "${r2}")" "${TMP_DIR}/${sample}_R2.fastq.gz"
        r1="${TMP_DIR}/${sample}_R1.fastq.gz"
        r2="${TMP_DIR}/${sample}_R2.fastq.gz"
    else
        # Single FASTQ path (single-end only)
        [[ "${ENDEDNESS}" == "single" ]] || {
            log_error "[${sample}] single path '${src}' requires -e single"; return 1; }
        [[ -f "${src}" ]] || { log_error "[${sample}] R1 not found: ${src}"; return 1; }
        ln -sf "$(realpath "${src}")" "${TMP_DIR}/${sample}_R1.fastq.gz"
        r1="${TMP_DIR}/${sample}_R1.fastq.gz"
    fi

    step_fastp      "${sample}" "${r1}" "${r2}" || return 1
    step_star_align "${sample}"                 || return 1
    cleanup_tmp     "${sample}"
}

# ── Step 4: Merge ReadsPerGene.out.tab → count_matrix.csv ────────────────────
merge_counts() {
    log_step "Merging count tables → ${OUT_DIR}/count_matrix.csv"

    local tab_files=()
    for f in "${OUT_DIR}"/STAR_out_*_ReadsPerGene.out.tab; do
        [[ -f "${f}" ]] && tab_files+=("${f}")
    done
    [[ ${#tab_files[@]} -gt 0 ]] || { log_error "No ReadsPerGene.out.tab files found in ${OUT_DIR}"; return 1; }

    # Build sample→file list for Python
    local py_args=""
    for f in "${tab_files[@]}"; do
        local sid; sid=$(basename "${f}" | sed 's/^STAR_out_//;s/_ReadsPerGene\.out\.tab$//')
        py_args="${py_args} ${sid}:${f}"
    done

    python3 - "${STRAND_COL}" ${py_args} <<'PYEOF'
import sys, csv
col = int(sys.argv[1])
entries = sys.argv[2:]

samples = {}
for entry in entries:
    sid, path = entry.split(":", 1)
    samples[sid] = path

gene_ids = None
counts = {}

for sid, path in samples.items():
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t")
                if not r[0].startswith("N_")]
    if gene_ids is None:
        gene_ids = [r[0] for r in rows]
    counts[sid] = {r[0]: int(r[col - 1]) for r in rows}

sample_ids = list(samples.keys())
with open("count_matrix.csv", "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["gene_id"] + sample_ids)
    for g in gene_ids:
        w.writerow([g] + [counts[s].get(g, 0) for s in sample_ids])

print(f"Wrote count_matrix.csv: {len(gene_ids)} genes × {len(sample_ids)} samples")
PYEOF

    mv -f count_matrix.csv "${OUT_DIR}/count_matrix.csv"
    log_info "Count matrix: ${OUT_DIR}/count_matrix.csv"
}

# ── Main ──────────────────────────────────────────────────────────────────────
main() {
    preflight

    if [[ "${MERGE_ONLY}" == true ]]; then
        merge_counts
        log_step "Done (merge only)"
        echo "Count matrix : ${OUT_DIR}/count_matrix.csv"
        return 0
    fi

    local n_ok=0 n_fail=0
    for sample in $(echo "${!SAMPLE_MAP[@]}" | tr ' ' '\n' | sort); do
        src="${SAMPLE_MAP[$sample]}"
        if run_sample "${sample}" "${src}"; then
            (( n_ok  += 1 ))
        else
            (( n_fail += 1 ))
            log_error "Sample ${sample} failed — continuing with remaining samples"
        fi
    done

    if [[ "${SKIP_MERGE}" == false ]]; then
        merge_counts
        echo "Count matrix : ${OUT_DIR}/count_matrix.csv"
    else
        log_info "--skip-merge set: skipping count matrix merge"
    fi

    log_step "Done"
    log_info "Completed: ${n_ok} | Failed: ${n_fail}"
    echo ""
    echo "QC reports   : ${QC_DIR}/"
    echo "BAM files    : ${OUT_DIR}/STAR_out_*_Aligned.sortedByCoord.out.bam"
}

main
