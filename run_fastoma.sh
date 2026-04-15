#!/usr/bin/env bash
# =============================================================================
# run_fastoma.sh
# Wrapper around the FastOMA Nextflow pipeline (DessimozLab/FastOMA)
# =============================================================================
# No other repositories required. Dependencies managed via environment.yml:
# nextflow, java, python3, biopython, pandas, wget, unzip, ncbi-datasets-cli.
# See README.md for full setup instructions.
#
# Usage:
#   ./run_fastoma.sh --config my_config.json --output-dir results/
#   ./run_fastoma.sh --species "human:Homo sapiens,chimp:Pan troglodytes" \
#                    --tree "((human:6,chimp:6):75,mouse:81);" \
#                    --proteins-dir /data/proteomes \
#                    --gff-dir /data/annotations \
#                    --output-dir results/
# =============================================================================

set -euo pipefail

SCRIPT_VERSION="1.0.0"

# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_THREADS="$(nproc 2>/dev/null || echo 4)"
DEFAULT_MEMORY="128G"
DEFAULT_MIN_PROTEIN_LENGTH=50
DEFAULT_MIN_SEQUENCE_LENGTH=60
DEFAULT_FILTER_GAP_RATIO_ROW=0.4
DEFAULT_FILTER_GAP_RATIO_COL=0.6
DEFAULT_NR_REPR_PER_HOG=3

LUCA_DB_URL="https://omabrowser.org/All/LUCA.h5"
FASTOMA_REPO="DessimozLab/FastOMA"
FASTOMA_REVISION="main"

# =============================================================================
# RUNTIME VARIABLES
# =============================================================================

ANALYSIS_NAME=""
CONFIG_FILE=""
SPECIES_LIST_STR=""
SPECIES_TREE=""

OUTPUT_DIR=""
WORK_DIR=""
PROTEINS_DIR=""
GFF_DIR=""
LUCA_DB_PATH="${HOME}/.fastoma/omamerdb.h5"
LOG_FILE=""

THREADS="${DEFAULT_THREADS}"
MEMORY="${DEFAULT_MEMORY}"
MIN_PROTEIN_LENGTH="${DEFAULT_MIN_PROTEIN_LENGTH}"
MIN_SEQUENCE_LENGTH="${DEFAULT_MIN_SEQUENCE_LENGTH}"
FILTER_GAP_RATIO_ROW="${DEFAULT_FILTER_GAP_RATIO_ROW}"
FILTER_GAP_RATIO_COL="${DEFAULT_FILTER_GAP_RATIO_COL}"
NR_REPR_PER_HOG="${DEFAULT_NR_REPR_PER_HOG}"

AUTO_DOWNLOAD=false
DRY_RUN=false
FORCE_REPROCESS=false
VERBOSE=false

# =============================================================================
# UTILITIES
# =============================================================================

log() {
    local ts; ts="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo "${ts} $*"
    [[ -n "${LOG_FILE:-}" ]] && echo "${ts} $*" >> "${LOG_FILE}" || true
}

error_exit() { echo "[ERROR] $*" >&2; exit 1; }
ensure_dir()  { mkdir -p "$@" || error_exit "Cannot create directory: $*"; }

# =============================================================================
# USAGE
# =============================================================================

usage() {
    cat << 'EOF'
run_fastoma.sh — FastOMA orthology inference wrapper

USAGE:
    ./run_fastoma.sh --config CONFIG.json --output-dir DIR [options]
    ./run_fastoma.sh --species "key:Name,..." --tree NEWICK \
        --proteins-dir DIR --gff-dir DIR --output-dir DIR [options]

INPUT (choose one):
    --config FILE             JSON config file
    --species "k:Name,..."    Comma-separated "key:Scientific name" pairs
    --tree NEWICK             Newick tree with branch lengths

DATA SOURCES:
    --proteins-dir DIR        Directory of per-species protein FASTAs
                              Files must be named <species_key>.fa
    --gff-dir DIR             Directory of per-species GFF3 files
                              Files must be named <species_key>.gff
    --auto-download           Download genome/protein data from NCBI

OUTPUT:
    --output-dir DIR          Required. Where results are written.
    --work-dir DIR            Temp/intermediate files
                              Default: OUTPUT_DIR/.work
    --luca-db PATH            Path to LUCA database (omamerdb.h5)
                              Default: ~/.fastoma/omamerdb.h5
                              Downloaded automatically on first run (~7 GB)
    --log-file PATH           Log file
                              Default: OUTPUT_DIR/run.log

RESOURCES:
    --threads N               CPU threads (default: nproc)
    --memory SIZE             Memory allocation, e.g. 64G (default: 128G)
    --min-protein-length N    Discard proteins shorter than N aa (default: 50)

FASTOMA / NEXTFLOW:
    --min-sequence-length N   default: 60
    --filter-gap-row F        default: 0.4
    --filter-gap-col F        default: 0.6
    --nr-repr-per-hog N       default: 3
    --fastoma-revision REV    FastOMA git tag or branch (default: main)

EXECUTION:
    --analysis-name NAME      Label for this run (default: timestamp)
    --dry-run                 Validate inputs, exit without running
    --force-reprocess         Re-run even if output already exists
    --verbose                 Print species tree and proteome list before run
    --help, -h                Show this message
    --version                 Print version

CONFIG FILE FORMAT:
    {
        "ANALYSIS_NAME": "primates_v1",
        "SPECIES_LIST": {
            "human":  "Homo sapiens",
            "chimp":  "Pan troglodytes",
            "gorilla":"Gorilla gorilla",
            "mouse":  "Mus musculus"
        },
        "SPECIES_TREE": "(((human:6,chimp:6):2,gorilla:8):72,mouse:80);"
    }

    Species keys must match the filenames in --proteins-dir and --gff-dir.

OUTPUT:
    OUTPUT_DIR/
    ├── OrthologousGroups.tsv     Main result
    ├── run_summary.txt           Run statistics
    ├── run.log                   Full log
    └── logs/
        ├── timeline.html
        ├── report.html
        └── trace.txt

EXAMPLES:
    # Pre-existing data
    ./run_fastoma.sh \
        --config configs/primates.json \
        --proteins-dir /data/proteomes \
        --gff-dir /data/annotations \
        --output-dir results/primates_v1 \
        --threads 32 --memory 256G

    # Auto-download from NCBI
    ./run_fastoma.sh \
        --config configs/primates.json \
        --auto-download \
        --output-dir results/primates_v1

    # Dry run
    ./run_fastoma.sh --config configs/primates.json \
        --proteins-dir /data/proteomes --gff-dir /data/annotations \
        --output-dir results/primates_v1 --dry-run

TROUBLESHOOTING:
    Java heap:       export NXF_OPTS="-Xms2G -Xmx8G"
    Resume failed:   rm -rf OUTPUT_DIR/.work/work
    Low proteins:    use --min-protein-length 30
    LUCA missing:    wget -O ~/.fastoma/omamerdb.h5 https://omabrowser.org/All/LUCA.h5
    Tree error:      Newick must end with ; and include branch lengths
EOF
}

version() { echo "run_fastoma.sh v${SCRIPT_VERSION}"; }

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================

parse_arguments() {
    [[ $# -eq 0 ]] && { usage; exit 1; }

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --config)               CONFIG_FILE="$2";            shift 2 ;;
            --species)              SPECIES_LIST_STR="$2";       shift 2 ;;
            --tree)                 SPECIES_TREE="$2";           shift 2 ;;
            --output-dir)           OUTPUT_DIR="$2";             shift 2 ;;
            --work-dir)             WORK_DIR="$2";               shift 2 ;;
            --proteins-dir)         PROTEINS_DIR="$2";           shift 2 ;;
            --gff-dir)              GFF_DIR="$2";                shift 2 ;;
            --luca-db)              LUCA_DB_PATH="$2";           shift 2 ;;
            --log-file)             LOG_FILE="$2";               shift 2 ;;
            --analysis-name)        ANALYSIS_NAME="$2";          shift 2 ;;
            --threads)              THREADS="$2";                shift 2 ;;
            --memory)               MEMORY="$2";                 shift 2 ;;
            --min-protein-length)   MIN_PROTEIN_LENGTH="$2";     shift 2 ;;
            --min-sequence-length)  MIN_SEQUENCE_LENGTH="$2";    shift 2 ;;
            --filter-gap-row)       FILTER_GAP_RATIO_ROW="$2";  shift 2 ;;
            --filter-gap-col)       FILTER_GAP_RATIO_COL="$2";  shift 2 ;;
            --nr-repr-per-hog)      NR_REPR_PER_HOG="$2";       shift 2 ;;
            --fastoma-revision)     FASTOMA_REVISION="$2";       shift 2 ;;
            --auto-download)        AUTO_DOWNLOAD=true;          shift ;;
            --dry-run)              DRY_RUN=true;                shift ;;
            --force-reprocess)      FORCE_REPROCESS=true;        shift ;;
            --verbose)              VERBOSE=true;                shift ;;
            --help|-h)              usage; exit 0 ;;
            --version)              version; exit 0 ;;
            *)                      echo "Unknown option: $1"; usage; exit 1 ;;
        esac
    done

    [[ -z "${OUTPUT_DIR}" ]] && error_exit "--output-dir is required"
    [[ -z "${CONFIG_FILE}" && -z "${SPECIES_LIST_STR}" ]] && \
        error_exit "Provide either --config or --species"
    [[ -n "${SPECIES_LIST_STR}" && -z "${SPECIES_TREE}" ]] && \
        error_exit "--tree is required with --species"
    if [[ "${AUTO_DOWNLOAD}" == false ]]; then
        [[ -z "${PROTEINS_DIR}" ]] && error_exit "--proteins-dir required (or use --auto-download)"
        [[ -z "${GFF_DIR}" ]]      && error_exit "--gff-dir required (or use --auto-download)"
    fi

    [[ -z "${WORK_DIR}" ]]      && WORK_DIR="${OUTPUT_DIR}/.work"
    [[ -z "${LOG_FILE}" ]]      && LOG_FILE="${OUTPUT_DIR}/run.log"
    [[ -z "${ANALYSIS_NAME}" ]] && ANALYSIS_NAME="fastoma_$(date +%Y%m%d_%H%M%S)"
}

# =============================================================================
# CONFIG
# =============================================================================

load_config() {
    # Build a temp config from CLI args if --species was used
    if [[ -n "${SPECIES_LIST_STR}" ]]; then
        CONFIG_FILE="$(mktemp /tmp/fastoma_config.XXXXXX.json)"
        local entries="" first=true
        IFS=',' read -ra pairs <<< "${SPECIES_LIST_STR}"
        for pair in "${pairs[@]}"; do
            IFS=':' read -ra kv <<< "${pair}"
            [[ "${first}" == true ]] && first=false || entries+=","
            entries+="\"${kv[0]// /}\": \"${kv[1]:-}\""
        done
        printf '{"ANALYSIS_NAME":"%s","SPECIES_LIST":{%s},"SPECIES_TREE":"%s"}\n' \
            "${ANALYSIS_NAME}" "${entries}" "${SPECIES_TREE}" > "${CONFIG_FILE}"
    fi

    # Validate
    python3 << PYEOF || error_exit "Invalid config: ${CONFIG_FILE}"
import json, sys
with open('${CONFIG_FILE}') as f:
    c = json.load(f)
for k in ('ANALYSIS_NAME', 'SPECIES_LIST', 'SPECIES_TREE'):
    if k not in c:
        print(f"Missing required field: {k}", file=sys.stderr); sys.exit(1)
if not c['SPECIES_LIST']:
    print("SPECIES_LIST is empty", file=sys.stderr); sys.exit(1)
if not c['SPECIES_TREE'].strip().endswith(';'):
    print("SPECIES_TREE must end with ;", file=sys.stderr); sys.exit(1)
PYEOF

    ANALYSIS_NAME="$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['ANALYSIS_NAME'])")"
    SPECIES_TREE="$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['SPECIES_TREE'])")"
    local species_list
    species_list="$(python3 -c "import json; print(', '.join(json.load(open('${CONFIG_FILE}'))['SPECIES_LIST']))")"
    log "Analysis: ${ANALYSIS_NAME} | Species: ${species_list}"
}

# =============================================================================
# ENVIRONMENT
# =============================================================================

validate_environment() {
    log "Checking dependencies..."
    for cmd in nextflow python3 wget awk grep sed; do
        command -v "${cmd}" &>/dev/null || error_exit "Required tool not found: ${cmd}. See README.md."
    done
    local jver; jver="$(java -version 2>&1 | head -1 | cut -d'"' -f2 | cut -d'.' -f1)"
    [[ "${jver}" -ge 11 ]] || error_exit "Java 11+ required (found ${jver})"
    python3 -c "import Bio, pandas" 2>/dev/null || \
        error_exit "Missing Python packages. Activate the conda environment: mamba activate fastoma"
    if [[ "${AUTO_DOWNLOAD}" == true ]]; then
        command -v datasets &>/dev/null || \
            error_exit "ncbi-datasets-cli not found. Activate the conda environment: mamba activate fastoma"
    fi
    local nver; nver="$(nextflow -version 2>&1 | grep -oP 'version \K[0-9.]+' || echo unknown)"
    log "Nextflow ${nver} | Java ${jver} | ${THREADS} threads | ${MEMORY} memory"
}

# =============================================================================
# LUCA DATABASE
# =============================================================================

ensure_luca_database() {
    if [[ -f "${LUCA_DB_PATH}" ]] && [[ "${FORCE_REPROCESS}" == false ]]; then
        log "LUCA database: ${LUCA_DB_PATH}"; return 0
    fi
    log "Downloading LUCA database (~7 GB) → ${LUCA_DB_PATH}"
    ensure_dir "$(dirname "${LUCA_DB_PATH}")"
    wget --show-progress -O "${LUCA_DB_PATH}.tmp" "${LUCA_DB_URL}" \
        || { rm -f "${LUCA_DB_PATH}.tmp"; error_exit "LUCA download failed"; }
    mv "${LUCA_DB_PATH}.tmp" "${LUCA_DB_PATH}"
    log "LUCA database ready"
}

# =============================================================================
# NCBI DATA DOWNLOAD  (--auto-download only)
# =============================================================================

download_species_data() {
    local species="$1" taxon="$2"
    local dl_dir="${WORK_DIR}/ncbi_downloads/${species}"
    ensure_dir "${dl_dir}" "${WORK_DIR}/raw_proteins" "${WORK_DIR}/raw_gff"

    local prot_out="${WORK_DIR}/raw_proteins/${species}.fa"
    local gff_out="${WORK_DIR}/raw_gff/${species}.gff"

    # Check that a reference assembly with annotation exists before downloading
    local summary
    summary="$(datasets summary genome taxon "${taxon}" \
        --reference --assembly-level complete,chromosome \
        --assembly-source RefSeq --annotated 2>/dev/null || true)"
    local count
    count="$(echo "${summary}" | python3 -c \
        "import json,sys; d=json.load(sys.stdin); print(d.get('total_count',0))" 2>/dev/null || echo 0)"
    if [[ "${count}" -eq 0 ]]; then
        error_exit "No annotated RefSeq reference assembly found for '${taxon}' (${species}).
  Try a more specific taxon name or provide protein/GFF files directly with --proteins-dir / --gff-dir."
    fi
    log "  Found ${count} assembly match(es) for ${species} (${taxon})"

    if [[ ! -f "${gff_out}" ]] || [[ "${FORCE_REPROCESS}" == true ]]; then
        log "  Downloading genome/GFF for ${species}..."
        datasets download genome taxon "${taxon}" \
            --reference --assembly-level complete,chromosome \
            --assembly-source RefSeq --annotated \
            --include genome,gff3 --filename "${dl_dir}/genome.zip" \
            || error_exit "Download failed for ${species}"
        unzip -q "${dl_dir}/genome.zip" -d "${dl_dir}/genome_data"
        find "${dl_dir}/genome_data/ncbi_dataset/data" \( -name "*.gff" -o -name "*.gff3" \) \
            | head -1 | xargs -I{} cp {} "${gff_out}"
        [[ -f "${gff_out}" ]] || error_exit "No GFF extracted for ${species} — assembly may lack annotation"
    fi

    if [[ ! -f "${prot_out}" ]] || [[ "${FORCE_REPROCESS}" == true ]]; then
        log "  Downloading proteins for ${species}..."
        datasets download genome taxon "${taxon}" \
            --reference --assembly-level complete,chromosome \
            --assembly-source RefSeq --annotated \
            --include protein --filename "${dl_dir}/protein.zip" \
            || error_exit "Protein download failed for ${species}"
        unzip -q "${dl_dir}/protein.zip" -d "${dl_dir}/protein_data"
        find "${dl_dir}/protein_data/ncbi_dataset/data" -name "*.faa" \
            | head -1 | xargs -I{} cp {} "${prot_out}"
        [[ -s "${prot_out}" ]] || error_exit "No proteins extracted for ${species} — assembly may lack annotation"
    fi

    rm -rf "${dl_dir}"
    PROTEINS_DIR="${WORK_DIR}/raw_proteins"
    GFF_DIR="${WORK_DIR}/raw_gff"
}

# =============================================================================
# GFF PROCESSING
# =============================================================================

process_gff() {
    local species="$1"
    local gff_in="${GFF_DIR}/${species}.gff"
    local gff_out="${WORK_DIR}/input/gff/${species}.gff"

    [[ -f "${gff_in}" ]] || error_exit "GFF not found: ${gff_in}"
    ensure_dir "$(dirname "${gff_out}")"

    # Prefix feature IDs with species key to prevent cross-species ID collisions
    awk -v sp="${species}" '
        BEGIN { OFS="\t" }
        /^#/ { print; next }
        $3 ~ /^(gene|mRNA|CDS)$/ {
            gsub(/;$/, "", $9)
            gsub(/ID=/,     "ID="     sp "_", $9)
            gsub(/Parent=/, "Parent=" sp "_", $9)
            print
        }
    ' "${gff_in}" > "${gff_out}"

    [[ -s "${gff_out}" ]] || error_exit "GFF output is empty for ${species}"
    log "  GFF processed: ${species}"
}

# =============================================================================
# PROTEIN FILTERING  (logic inlined — no external script needed)
# =============================================================================

filter_proteins() {
    local species="$1"
    local prot_in="${PROTEINS_DIR}/${species}.fa"
    local prot_out="${WORK_DIR}/input/proteome/${species}.fa"

    [[ -f "${prot_in}" ]] || error_exit "Protein FASTA not found: ${prot_in}"
    ensure_dir "$(dirname "${prot_out}")"

    python3 - "${prot_in}" "${prot_out}" "${MIN_PROTEIN_LENGTH}" << 'PYEOF'
import sys
from Bio import SeqIO
infile, outfile, min_len = sys.argv[1], sys.argv[2], int(sys.argv[3])
kept = [r for r in SeqIO.parse(infile, "fasta") if len(r.seq) >= min_len]
SeqIO.write(kept, outfile, "fasta")
print(f"  {len(kept)} proteins >= {min_len} aa retained", flush=True)
PYEOF

    local count; count="$(grep -c '^>' "${prot_out}" 2>/dev/null || echo 0)"
    [[ "${count}" -lt 5 ]] && \
        log "  WARNING: only ${count} proteins for ${species} after filtering — consider --min-protein-length 30"
    log "  Proteins: ${species} — ${count} sequences retained"
}

# =============================================================================
# SPECIES TREE
# =============================================================================

write_species_tree() {
    local tree_out="${WORK_DIR}/input/species_tree.nwk"

    # FastOMA requires internal nodes to be named
    python3 - "${SPECIES_TREE}" "${tree_out}" << 'PYEOF'
import re, sys

def name_internal_nodes(newick):
    counter = [1]
    def replace(m):
        label = f")internal_{counter[0]}"
        counter[0] += 1
        return label
    return re.sub(r'\)(?![a-zA-Z0-9_])(?=[:,);])', replace, newick)

tree, outfile = sys.argv[1], sys.argv[2]
with open(outfile, 'w') as f:
    f.write(name_internal_nodes(tree) + '\n')
PYEOF

    log "Species tree written"
    [[ "${VERBOSE}" == true ]] && { echo "  Tree:"; cat "${WORK_DIR}/input/species_tree.nwk"; }
}

# =============================================================================
# NEXTFLOW CONFIG
# =============================================================================

write_nextflow_config() {
    local logs_dir="${OUTPUT_DIR}/logs"
    ensure_dir "${logs_dir}"

    # Auto-size Nextflow heap at 1/4 of available RAM, capped at 8 GB
    local avail_gb; avail_gb="$(free -g 2>/dev/null | awk '/^Mem:/ {print $2}' || echo 64)"
    local nxf_heap=$(( avail_gb / 4 ))
    [[ "${nxf_heap}" -lt 2 ]] && nxf_heap=2
    [[ "${nxf_heap}" -gt 8 ]] && nxf_heap=8
    export NXF_OPTS="-Xms1G -Xmx${nxf_heap}G"

    cat > "${WORK_DIR}/nextflow.config" << EOF
process {
    executor = 'local'
    cpus     = ${THREADS}
    memory   = '${MEMORY}'
}

env {
    JAVA_OPTS = "-Xms1G -Xmx${MEMORY} -XX:+UseG1GC -XX:MaxGCPauseMillis=200"
}

timeline { enabled = true; file = '${logs_dir}/timeline.html' }
report   { enabled = true; file = '${logs_dir}/report.html'   }
trace    { enabled = true; file = '${logs_dir}/trace.txt'     }
dag      { enabled = true; file = '${logs_dir}/dag.dot'; overwrite = true }

workDir = '${WORK_DIR}/work'
EOF
}

# =============================================================================
# RUN FASTOMA
# =============================================================================

run_fastoma() {
    log "Running FastOMA (${FASTOMA_REPO} @ ${FASTOMA_REVISION})..."
    [[ "${FORCE_REPROCESS}" == true ]] && rm -rf "${WORK_DIR}/work" 2>/dev/null || true

    # Symlink LUCA DB into the input folder
    ln -sf "${LUCA_DB_PATH}" "${WORK_DIR}/input/omamerdb.h5"

    if [[ "${VERBOSE}" == true ]]; then
        log "Proteome files:"
        ls -lh "${WORK_DIR}/input/proteome/"
    fi

    nextflow run "${FASTOMA_REPO}" \
        -r             "${FASTOMA_REVISION}" \
        -c             "${WORK_DIR}/nextflow.config" \
        --input_folder         "${WORK_DIR}/input" \
        --output_folder        "${OUTPUT_DIR}" \
        --species_tree         "${WORK_DIR}/input/species_tree.nwk" \
        --database             "${WORK_DIR}/input/omamerdb.h5" \
        --min_sequence_length  "${MIN_SEQUENCE_LENGTH}" \
        --filter_gap_ratio_row "${FILTER_GAP_RATIO_ROW}" \
        --filter_gap_ratio_col "${FILTER_GAP_RATIO_COL}" \
        --nr_repr_per_hog      "${NR_REPR_PER_HOG}" \
        --report false \
        -profile standard \
        -resume \
        || error_exit "FastOMA pipeline failed. Check ${LOG_FILE} for details."

    log "FastOMA complete"
}

# =============================================================================
# SUMMARY
# =============================================================================

write_summary() {
    local orth="${OUTPUT_DIR}/OrthologousGroups.tsv"
    [[ -f "${orth}" ]] || { log "WARNING: OrthologousGroups.tsv not found in output"; return; }
    local total; total="$(( $(wc -l < "${orth}") - 1 ))"  # subtract header

    {
        echo "FastOMA Run Summary"
        echo "==================="
        echo "Date:             $(date)"
        echo "Analysis:         ${ANALYSIS_NAME}"
        echo "FastOMA revision: ${FASTOMA_REVISION}"
        echo "Threads:          ${THREADS}"
        echo "Memory:           ${MEMORY}"
        echo "Min protein len:  ${MIN_PROTEIN_LENGTH} aa"
        echo ""
        echo "Orthology groups: ${total}"
        echo "Output:           ${OUTPUT_DIR}/OrthologousGroups.tsv"
        echo "Log:              ${LOG_FILE}"
    } | tee "${OUTPUT_DIR}/run_summary.txt"
}

# =============================================================================
# MAIN
# =============================================================================

main() {
    parse_arguments "$@"

    ensure_dir "${OUTPUT_DIR}" "${WORK_DIR}"
    ensure_dir "${WORK_DIR}/input/proteome" "${WORK_DIR}/input/gff"

    exec > >(tee -a "${LOG_FILE}") 2>&1

    echo "========================================"
    echo " fastoma-wrapper v${SCRIPT_VERSION}"
    echo "========================================"

    load_config
    validate_environment

    if [[ "${DRY_RUN}" == true ]]; then
        log "Dry run — all inputs valid. Exiting without execution."; return 0
    fi

    ensure_luca_database

    # Read species list from config
    mapfile -t species_keys < <(
        python3 -c "import json; [print(k) for k in json.load(open('${CONFIG_FILE}'))['SPECIES_LIST']]"
    )

    # Download data if requested
    if [[ "${AUTO_DOWNLOAD}" == true ]]; then
        log "Downloading species data from NCBI..."
        for sp in "${species_keys[@]}"; do
            local taxon; taxon="$(python3 -c \
                "import json; print(json.load(open('${CONFIG_FILE}'))['SPECIES_LIST']['${sp}'])")"
            download_species_data "${sp}" "${taxon}"
        done
    fi

    # Pre-process inputs
    log "Pre-processing inputs..."
    for sp in "${species_keys[@]}"; do
        process_gff    "${sp}"
        filter_proteins "${sp}"
    done

    write_species_tree
    write_nextflow_config
    run_fastoma
    write_summary

    echo ""
    echo "========================================"
    echo " Complete."
    echo " Results: ${OUTPUT_DIR}/OrthologousGroups.tsv"
    echo " Log:     ${LOG_FILE}"
    echo "========================================"
}

main "$@"
