#!/usr/bin/env bash
#
# fetch_host_genome.sh - download a host reference genome and build the
#                        minimap2 indexes bioThreatFinder uses for host
#                        depletion.
#
# Host depletion is the first real filter in the pipeline. If it is weak, host
# reads leak into the classifier and read out as low-abundance pathogens - which
# in a biosurveillance context means false alarms. So the indexes are built
# once, deliberately, and pinned by path in packs/<pack>/pack.yaml.
#
# Paths written here MUST match packs/human/pack.yaml exactly:
#   host.fasta              ${db_root}/bioThreatFinder/host/human/chm13v2.0.fa.gz
#   host.minimap2_index_sr  ${db_root}/bioThreatFinder/host/human/chm13v2.0.sr.mmi
#   host.minimap2_index_lr  ${db_root}/bioThreatFinder/host/human/chm13v2.0.map-ont.mmi
#
set -euo pipefail

# An ERR trap, because `set -e` alone exits SILENTLY.
#
# This build died mid-way through the multi-hour bacteria download and left a log whose
# last line was ordinary progress output - no error, no exit code, nothing to diagnose.
# `set -e` had done its job and aborted, but said nothing. For a job measured in hours
# that is the difference between "resume from the failing step" and "discover tomorrow
# that nothing happened".
_btf_on_err() {
    local rc=$? line=${1:-?}
    printf '[%s] ERROR: aborted at line %s with exit %s: %s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$line" "$rc" "${BASH_COMMAND:-unknown}" >&2
    printf '[%s] ERROR: completed steps are stamped, so re-running resumes rather than restarting.\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >&2
    exit "$rc"
}
trap '_btf_on_err $LINENO' ERR
# --------------------------------------------------------------------------
# Defaults
# --------------------------------------------------------------------------
DB_ROOT="/media/36tb/db"
# minimap2 indexing is dominated by a poorly-parallel sort; past ~16 threads the
# return is negligible while the memory high-water mark keeps climbing. The box
# has 64 cores but there is no reason to use them all here.
THREADS=16
SPECIES_LIST=""
FORCE=0
KEEP_UNCOMPRESSED=0

# Minimum free space, in GiB, required under DB_ROOT before we start.
# Human: 1 GB .fa.gz + ~8 GB per minimap2 index + working headroom.
MIN_FREE_GIB_HUMAN=40
# Non-human hosts are looped, so the requirement scales with the list; this is
# the per-genome allowance used to compute the total.
MIN_FREE_GIB_PER_GENOME=40

# --------------------------------------------------------------------------
# Source URLs
# --------------------------------------------------------------------------
# T2T-CHM13v2.0 = GCF_009914755.1 (includes the HG002 chrY).
#
# Why CHM13 and not GRCh38: GRCh38 leaves centromeric/satellite and rDNA
# regions as gaps or models. Reads from those regions fail to map, survive
# depletion, and then hit bacterial/fungal repeats in the classifier. CHM13
# closes them, which measurably reduces host carry-over.
#
# VERIFY: NCBI FTP layout is stable but the exact directory name embeds the
# assembly name; re-check this path if the download 404s.
CHM13_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0"
CHM13_FASTA_URL="${CHM13_BASE}/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
CHM13_MD5_URL="${CHM13_BASE}/md5checksums.txt"
# NCBI accession, used when the `datasets` CLI route is taken instead of FTP.
CHM13_ACCESSION="GCF_009914755.1"

# Documented fallback: the T2T consortium's own "analysis set" copy. Same
# assembly, chrEBV/decoy-free, UCSC-style chr* names rather than RefSeq NC_*.
# VERIFY: bucket path; also note the different sequence naming, which matters if
# anything downstream joins on contig names.
CHM13_FALLBACK_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"

usage() {
    cat <<'EOF'
Usage:
  fetch_host_genome.sh <human|animal|plant> [options]

Downloads a host reference genome and builds the two minimap2 indexes used for
host depletion (short-read and ONT).

Positional:
  human | animal | plant   Which domain pack's host(s) to fetch.
                           `human` is fully implemented (T2T-CHM13v2.0).
                           `animal` and `plant` require --species-list, because
                           those packs have many hosts and no single reference.

Options:
  --db-root DIR       Shared database root. Artifacts are written under
                      DIR/bioThreatFinder/host/<domain>/.
                      (default: /media/36tb/db)
  --threads N         Threads for minimap2 indexing (default: 16)
  --species-list FILE TSV of hosts for animal/plant. See "Species list" below.
  --keep-uncompressed Keep the decompressed FASTA (some tools want plain .fa).
                      Off by default: it costs ~3.1 GB for human and minimap2
                      reads the .gz directly.
  --force             Rebuild even if outputs already exist.
  -h, --help          This message.

Species list (animal/plant), tab-separated, '#' comments allowed:
  label <TAB> taxid <TAB> accession
e.g.
  Bos taurus      9913    GCF_002263795.3
  Gallus gallus   9031    GCF_016699485.2

Behaviour:
  * Idempotent - any step whose output already exists and is non-empty is
    skipped with a message. Use --force to override.
  * Checksums are verified against the NCBI md5checksums.txt when available.
  * Downloads land on a .part file and are renamed only on success, so an
    interrupted run never leaves a truncated FASTA that looks complete.

Examples:
  fetch_host_genome.sh human
  fetch_host_genome.sh human --threads 32
  fetch_host_genome.sh animal --species-list packs/animal/hosts.tsv
EOF
}

# --------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------
log()  { printf '[%s] %s\n'        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
warn() { printf '[%s] WARNING: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n'   "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; exit 1; }
skip() { printf '[%s] SKIP: %s\n'    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }

# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------
[[ $# -ge 1 ]] || { usage; exit 2; }

DOMAIN=""
case "$1" in
    human|animal|plant) DOMAIN="$1"; shift ;;
    -h|--help)          usage; exit 0 ;;
    *)                  usage; die "first argument must be human, animal or plant (got '$1')" ;;
esac

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db-root)           DB_ROOT="${2:?--db-root needs a value}"; shift 2 ;;
        --threads)           THREADS="${2:?--threads needs a value}"; shift 2 ;;
        --species-list)      SPECIES_LIST="${2:?--species-list needs a value}"; shift 2 ;;
        --keep-uncompressed) KEEP_UNCOMPRESSED=1; shift ;;
        --force)             FORCE=1; shift ;;
        -h|--help)           usage; exit 0 ;;
        *)                   usage; die "unknown option: $1" ;;
    esac
done

[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || die "--threads must be a positive integer"

HOST_DIR="${DB_ROOT}/bioThreatFinder/host/${DOMAIN}"

# --------------------------------------------------------------------------
# Tool resolution
# --------------------------------------------------------------------------
# The shared conda prefix is not necessarily on PATH for non-interactive runs
# (cron, Nextflow), so look there explicitly before giving up.
CONDA_BIN="/home/bioinfo/miniconda3/bin"

resolve_tool() {
    # resolve_tool <name> -> prints absolute path, or returns 1
    local name="$1" p
    if p="$(command -v "$name" 2>/dev/null)"; then printf '%s\n' "$p"; return 0; fi
    if [[ -x "${CONDA_BIN}/${name}" ]]; then printf '%s\n' "${CONDA_BIN}/${name}"; return 0; fi
    return 1
}

MINIMAP2="$(resolve_tool minimap2 || true)"
if [[ -z "$MINIMAP2" ]]; then
    # Last resort: the copy bundled with the ONT dorado server. It works, but it
    # may be an older release than the one the pipeline was validated against,
    # so say so loudly rather than silently using it.
    if [[ -x /home/bioinfo/prog/ont-dorado-server/bin/minimap2 ]]; then
        MINIMAP2=/home/bioinfo/prog/ont-dorado-server/bin/minimap2
        warn "using minimap2 bundled with ont-dorado-server: ${MINIMAP2}"
        warn "this may be an old build; prefer: ${CONDA_BIN}/mamba install -c bioconda minimap2"
    else
        die "minimap2 not found. Install it, e.g.: ${CONDA_BIN}/mamba install -c bioconda minimap2"
    fi
fi

CURL="$(resolve_tool curl || true)"
[[ -n "$CURL" ]] || die "curl not found"

DATASETS="$(resolve_tool datasets || true)"   # optional; FTP is the primary route here

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

# check_free_space <dir> <required_gib>
# The shared volume sits at 97% capacity, so a build that silently runs the
# filesystem dry would take other databases down with it. Refuse up front.
check_free_space() {
    local dir="$1" need_gib="$2" probe avail_gib
    probe="$dir"
    while [[ ! -d "$probe" && "$probe" != "/" ]]; do probe="$(dirname "$probe")"; done
    avail_gib="$(df -BG --output=avail "$probe" | tail -1 | tr -dc '0-9')"
    [[ -n "$avail_gib" ]] || die "could not determine free space for ${probe}"
    log "free space on $(df --output=target "$probe" | tail -1 | tr -d ' '): ${avail_gib} GiB (need ${need_gib} GiB)"
    if (( avail_gib < need_gib )); then
        die "insufficient free space: ${avail_gib} GiB available under ${probe}, need ${need_gib} GiB.
       Free space or pass --db-root pointing at a larger volume.
       Nothing was downloaded."
    fi
}

# exists_nonempty <path>
exists_nonempty() { [[ -s "$1" ]]; }

# download <url> <dest>
# Atomic: writes <dest>.part, renames on success. --retry survives the flaky
# stretches that NCBI FTP-over-HTTPS has under load.
download() {
    local url="$1" dest="$2"
    log "downloading ${url}"
    log "        -> ${dest}"
    "$CURL" --fail --location --show-error --silent \
            --retry 5 --retry-delay 10 --retry-connrefused \
            --connect-timeout 30 \
            --continue-at - \
            --output "${dest}.part" "$url" \
        || { rm -f "${dest}.part"; return 1; }
    [[ -s "${dest}.part" ]] || { rm -f "${dest}.part"; return 1; }
    mv -f "${dest}.part" "$dest"
}

# unpack_zip <zipfile> <destdir>
# `unzip` is not installed everywhere; python3 is, so fall back to it rather
# than failing after a multi-GB download has already succeeded.
unpack_zip() {
    local zf="$1" dest="$2"
    if command -v unzip >/dev/null 2>&1; then
        unzip -q -o "$zf" -d "$dest"
    elif command -v python3 >/dev/null 2>&1; then
        python3 -m zipfile -e "$zf" "$dest"
    else
        die "neither unzip nor python3 is available to unpack ${zf}"
    fi
}

# verify_md5 <file> <expected_md5>
verify_md5() {
    local file="$1" expected="$2" got
    got="$(md5sum "$file" | awk '{print $1}')"
    if [[ "$got" != "$expected" ]]; then
        warn "md5 MISMATCH for ${file}"
        warn "  expected ${expected}"
        warn "  got      ${got}"
        return 1
    fi
    log "md5 OK for $(basename "$file")"
}

# ncbi_md5_for <md5checksums_url> <basename> -> prints md5 or nothing
# NCBI's md5checksums.txt lines look like: "<md5>  ./<filename>"
ncbi_md5_for() {
    local md5url="$1" fname="$2"
    "$CURL" --fail --location --silent --retry 3 --connect-timeout 20 "$md5url" 2>/dev/null \
        | awk -v f="$fname" '{ sub(/^\.\//, "", $2); if ($2 == f) { print $1; exit } }'
}

# build_mmi <fasta> <out.mmi> <preset>
# The k-mer and window size are baked into a minimap2 index at build time, so an
# index built with -x sr is NOT usable for map-ont and vice versa. Hence two
# separate indexes rather than one shared .mmi - this is the single most common
# way to get silently terrible host depletion.
build_mmi() {
    local fasta="$1" out="$2" preset="$3"
    if exists_nonempty "$out" && (( FORCE == 0 )); then
        skip "minimap2 index exists: ${out}"
        return 0
    fi
    log "building minimap2 index (-x ${preset}, ${THREADS} threads): ${out}"
    "$MINIMAP2" -x "$preset" -t "$THREADS" -d "${out}.part" "$fasta"
    [[ -s "${out}.part" ]] || { rm -f "${out}.part"; die "minimap2 produced an empty index for ${out}"; }
    mv -f "${out}.part" "$out"
    log "done: ${out} ($(du -h "$out" | cut -f1))"
}

# --------------------------------------------------------------------------
# human
# --------------------------------------------------------------------------
fetch_human() {
    local fasta="${HOST_DIR}/chm13v2.0.fa.gz"
    local mmi_sr="${HOST_DIR}/chm13v2.0.sr.mmi"
    local mmi_ont="${HOST_DIR}/chm13v2.0.map-ont.mmi"

    check_free_space "$DB_ROOT" "$MIN_FREE_GIB_HUMAN"
    mkdir -p "$HOST_DIR"

    if exists_nonempty "$fasta" && (( FORCE == 0 )); then
        skip "host FASTA exists: ${fasta} ($(du -h "$fasta" | cut -f1))"
    else
        local want_md5
        want_md5="$(ncbi_md5_for "$CHM13_MD5_URL" "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" || true)"
        if [[ -n "$want_md5" ]]; then
            log "expected md5 from NCBI md5checksums.txt: ${want_md5}"
        else
            warn "could not read md5checksums.txt; download will not be checksum-verified"
        fi

        if download "$CHM13_FASTA_URL" "$fasta"; then
            :
        elif [[ -n "$DATASETS" ]]; then
            # datasets emits a zip with a nested directory layout, so this route
            # needs an extra unpack step and yields no md5 to check against.
            warn "FTP download failed; trying the datasets CLI for ${CHM13_ACCESSION}"
            local tmpd; tmpd="$(mktemp -d "${HOST_DIR}/.datasets.XXXXXX")"
            "$DATASETS" download genome accession "$CHM13_ACCESSION" \
                --include genome --filename "${tmpd}/chm13.zip" \
                || die "datasets download failed for ${CHM13_ACCESSION}"
            unpack_zip "${tmpd}/chm13.zip" "$tmpd"
            local found
            found="$(find "$tmpd" -name '*_genomic.fna' -o -name '*.fna' | head -1)"
            [[ -n "$found" ]] || die "no FASTA inside the datasets archive"
            pigz -p "$THREADS" -c "$found" > "${fasta}.part" 2>/dev/null \
                || gzip -c "$found" > "${fasta}.part"
            mv -f "${fasta}.part" "$fasta"
            rm -rf "$tmpd"
            want_md5=""   # checksum from NCBI does not apply to a re-compressed file
        else
            warn "FTP download failed and datasets CLI is unavailable."
            warn "trying the T2T consortium fallback (different contig naming - see comments)"
            download "$CHM13_FALLBACK_URL" "$fasta" \
                || die "all download routes for T2T-CHM13v2.0 failed. Verify the URLs in this script."
            want_md5=""
        fi

        if [[ -n "$want_md5" ]]; then
            verify_md5 "$fasta" "$want_md5" || die "checksum verification failed; ${fasta} left in place for inspection"
        fi
        # A sanity check that costs nothing and catches a truncated or
        # HTML-error-page download that happens to be non-empty.
        gzip -t "$fasta" || die "${fasta} is not a valid gzip stream"
        log "host FASTA ready: ${fasta} ($(du -h "$fasta" | cut -f1))"
    fi

    if (( KEEP_UNCOMPRESSED == 1 )); then
        local plain="${HOST_DIR}/chm13v2.0.fa"
        if exists_nonempty "$plain" && (( FORCE == 0 )); then
            skip "uncompressed FASTA exists: ${plain}"
        else
            log "decompressing to ${plain} (--keep-uncompressed)"
            if command -v pigz >/dev/null 2>&1; then
                pigz -dc -p "$THREADS" "$fasta" > "${plain}.part"
            else
                gzip -dc "$fasta" > "${plain}.part"
            fi
            mv -f "${plain}.part" "$plain"
        fi
    fi

    # minimap2 indexes both read straight from the gzipped FASTA.
    build_mmi "$fasta" "$mmi_sr"  sr
    build_mmi "$fasta" "$mmi_ont" map-ont

    # Cross-check against pack.yaml so a path drift is caught here rather than
    # three stages downstream when depletion silently finds no index.
    local pack_yaml
    pack_yaml="$(dirname "$0")/../packs/human/pack.yaml"
    if [[ -f "$pack_yaml" ]]; then
        local p
        for p in "chm13v2.0.fa.gz" "chm13v2.0.sr.mmi" "chm13v2.0.map-ont.mmi"; do
            grep -q "$p" "$pack_yaml" \
                || warn "packs/human/pack.yaml does not mention ${p}; paths may have drifted"
        done
    fi
}

# --------------------------------------------------------------------------
# animal / plant
# --------------------------------------------------------------------------
# These packs have MANY hosts (livestock and poultry species for animal; dozens
# of crops for plant), so there is no single reference to hard-code. The host
# set is therefore data, supplied as a species list, exactly like watchlists.
#
# TODO: curate packs/animal/hosts.tsv and packs/plant/hosts.tsv, choosing one
#       assembly accession per host species. Deliberately NOT guessed here:
#       picking the wrong assembly (a scaffold-level submission, or the wrong
#       breed/cultivar) degrades depletion in ways that are hard to notice.
#       For each host decide and record:
#         - RefSeq accession (prefer GCF_ over GCA_, chromosome-level)
#         - whether organelle sequences must be present (plant chloroplast and
#           mitochondrial reads are a large fraction of plant total RNA/DNA
#           libraries; leaving them out means they survive depletion)
#       Once those files exist, this function needs no further changes.
fetch_from_species_list() {
    [[ -n "$SPECIES_LIST" ]] || die "${DOMAIN} requires --species-list FILE (see --help). \
No hosts are hard-coded for ${DOMAIN}: see the TODO in this script."
    [[ -f "$SPECIES_LIST" ]] || die "species list not found: ${SPECIES_LIST}"

    local n_hosts
    n_hosts="$(grep -cvE '^[[:space:]]*(#|$)' "$SPECIES_LIST" || true)"
    (( n_hosts > 0 )) || die "species list ${SPECIES_LIST} has no data rows"
    log "${DOMAIN}: ${n_hosts} host(s) listed in ${SPECIES_LIST}"

    check_free_space "$DB_ROOT" "$(( n_hosts * MIN_FREE_GIB_PER_GENOME ))"
    mkdir -p "$HOST_DIR"

    [[ -n "$DATASETS" ]] || die "the ${DOMAIN} route downloads by accession and needs the NCBI datasets CLI.
       Install it:  ${CONDA_BIN}/mamba install -c conda-forge ncbi-datasets-cli
       (or see scripts/build_pack_refs.py for the documented curl/FTP fallback)"

    local n_ok=0 n_fail=0
    local label taxid accession
    while IFS=$'\t' read -r label taxid accession _rest; do
        # Note the explicit if/then: with `set -e`, a bare `[[ ... ]] && continue`
        # whose test is false would abort the whole script.
        if [[ -z "${label:-}" || "$label" == \#* ]]; then continue; fi
        if [[ -z "${accession:-}" ]]; then
            warn "row '${label}' has no accession; skipping"
            n_fail=$((n_fail+1)); continue
        fi

        # Filesystem-safe stem: "Bos taurus" -> Bos_taurus
        local stem; stem="$(printf '%s' "$label" | tr -c 'A-Za-z0-9._-' '_' | sed 's/__*/_/g; s/^_//; s/_$//')"
        local sub="${HOST_DIR}/${stem}"
        local fasta="${sub}/${stem}.fa.gz"
        local mmi_sr="${sub}/${stem}.sr.mmi"
        local mmi_ont="${sub}/${stem}.map-ont.mmi"
        mkdir -p "$sub"

        log "--- ${label} (taxid ${taxid:-?}, ${accession}) ---"
        if exists_nonempty "$fasta" && (( FORCE == 0 )); then
            skip "FASTA exists: ${fasta}"
        else
            local tmpd; tmpd="$(mktemp -d "${sub}/.datasets.XXXXXX")"
            if ! "$DATASETS" download genome accession "$accession" \
                    --include genome --filename "${tmpd}/g.zip" >/dev/null; then
                warn "datasets download failed for ${label} (${accession})"
                rm -rf "$tmpd"; n_fail=$((n_fail+1)); continue
            fi
            unpack_zip "${tmpd}/g.zip" "$tmpd" || { warn "unzip failed for ${label}"; rm -rf "$tmpd"; n_fail=$((n_fail+1)); continue; }
            local found; found="$(find "$tmpd" -name '*.fna' | head -1)"
            if [[ -z "$found" ]]; then
                warn "no FASTA in archive for ${label}"; rm -rf "$tmpd"; n_fail=$((n_fail+1)); continue
            fi
            if command -v pigz >/dev/null 2>&1; then
                pigz -p "$THREADS" -c "$found" > "${fasta}.part"
            else
                gzip -c "$found" > "${fasta}.part"
            fi
            mv -f "${fasta}.part" "$fasta"
            rm -rf "$tmpd"
        fi

        build_mmi "$fasta" "$mmi_sr"  sr
        build_mmi "$fasta" "$mmi_ont" map-ont
        n_ok=$((n_ok+1))
    done < "$SPECIES_LIST"

    log "${DOMAIN}: ${n_ok} host(s) built, ${n_fail} failed"
    (( n_fail == 0 )) || warn "${n_fail} host(s) failed - rerun to retry (completed hosts are skipped)"
}

# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
log "domain=${DOMAIN} db_root=${DB_ROOT} threads=${THREADS}"
log "output directory: ${HOST_DIR}"

case "$DOMAIN" in
    human)        fetch_human ;;
    animal|plant) fetch_from_species_list ;;
esac

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
echo >&2
log "=============================== SUMMARY ==============================="
if [[ -d "$HOST_DIR" ]]; then
    while IFS= read -r line; do log "  ${line}"; done < <(
        find "$HOST_DIR" -maxdepth 2 -type f \( -name '*.fa.gz' -o -name '*.fa' -o -name '*.mmi' \) \
            -printf '%s\t%p\n' | sort -k2 | awk -F'\t' '{
                sz=$1; unit="B";
                if (sz>=1073741824) { sz=sz/1073741824; unit="GiB" }
                else if (sz>=1048576) { sz=sz/1048576; unit="MiB" }
                printf "%8.1f %-3s  %s\n", sz, unit, $2
            }'
    )
    log "  total for ${HOST_DIR}: $(du -sh "$HOST_DIR" | cut -f1)"
fi
log "free space remaining: $(df -h --output=avail "$DB_ROOT" | tail -1 | tr -d ' ')"
log "======================================================================="
