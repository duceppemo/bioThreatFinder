#!/usr/bin/env bash
#
# build_sylph_db.sh - build the sylph ANI/containment database for a pack.
#
# Target path (must match screen.sylph_db in packs/*/pack.yaml):
#   ${BTF_DB}/sylph/<pack>_pack_refs.syldb     (default BTF_DB=/media/8tb/bioThreatFinder)
#
# WHY THIS SCRIPT EXISTS, LATE
# ----------------------------
# Every other database in this project has a builder. This one did not. pack.yaml credited
# "scripts/build_sylph_db.sh" for the index while no such file existed, and the 184 MB
# .syldb on disk was produced by a hand-run sylph sketch on 2026-09-29 with no recorded
# command - so the only arm with no builder was also the only one whose provenance could
# not be checked. A reader following that comment went looking for a file that was not
# there. This is that file.
#
# THE FAILURE MODE THIS BUILD HAS TO PREVENT IS SILENT, NOT LOUD
# --------------------------------------------------------------
# sylph reports GENOME FILENAMES, never taxids. bin/parse_screen.py recovers a taxid by
# tokenising the reported path on "/", "." and "_" and taking any all-digit token of two or
# more characters (parse_sylph). That works because this pack's references are named
# <taxid>_<name>.fna.gz.
#
# Build the index from anything else - the GTDB distribution is the obvious temptation,
# since it is prebuilt and large - and every filename is a GTDB accession, no token is an
# NCBI taxid, and EVERY lookup misses. The arm then contributes nothing while the report
# still reads "no sylph/ANI support", which is indistinguishable from a clean negative.
# That is a false negative on every watchlist organism, reported as a result.
#
# So this script refuses inputs it cannot map, rather than building an index that would be
# quietly useless (--check-names, on by default).
#
# PARAMETERS ARE NOT A FREE CHOICE
# --------------------------------
# modules/local/detect.nf runs `sylph profile ${db} ...` and passes neither -c nor -k, so
# the query side uses sylph's defaults (c=200, k=31, min-spacing=30). A database built at
# any other -c or -k does not match the query sketch. This script builds at the defaults and
# then VERIFIES the written database reports them back, because the cost of a mismatch is
# again silence rather than an error.
#
#   scripts/build_sylph_db.sh                      # human pack, default paths
#   scripts/build_sylph_db.sh --pack human --force # rebuild over an existing index
#   scripts/build_sylph_db.sh --dry-run            # print the plan and the command
#
set -eEuo pipefail

# An ERR trap, because `set -e` alone exits SILENTLY - and it needs `set -E` (errtrace) or
# the trap is not inherited by shell functions, which is where every step here runs. That
# exact omission made three other builds in this repo die without printing anything.
_btf_on_err() {
    local rc=$? line=${1:-?}
    printf '[%s] ERROR: aborted at line %s with exit %s: %s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$line" "$rc" "${BASH_COMMAND:-unknown}" >&2
    exit "$rc"
}
trap '_btf_on_err $LINENO' ERR

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
BTF_DB="${BTF_DB:-/media/8tb/bioThreatFinder}"
PACK="human"
REFS_DIR=""
OUT=""
SYLPH="${SYLPH:-sylph}"
THREADS="$(nproc 2>/dev/null || echo 8)"
# Defaults, deliberately: `sylph profile` in modules/local/detect.nf passes neither.
SYLPH_C=200
SYLPH_K=31
SYLPH_MIN_SPACING=30
DRY_RUN=0
FORCE=0
CHECK_NAMES=1

log()  { printf '[%s] %s\n'        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; exit 1; }

usage() {
    sed -n '2,42p' "$0" | sed 's/^# \{0,1\}//'
    cat >&2 <<EOF

Options:
  --pack NAME        pack whose refs are sketched        [${PACK}]
  --btf-db DIR       project database root               [${BTF_DB}]
  --refs-dir DIR     override the pack's refs directory
  --out FILE         override the output .syldb path
  --threads N        threads for sylph sketch            [${THREADS}]
  --sylph PATH       sylph binary                        [${SYLPH}]
  --no-check-names   build even if filenames carry no taxid (see above; you are
                     asking for an index whose hits cannot be mapped)
  --force            overwrite an existing database
  --dry-run          print the plan and the command, build nothing
  -h, --help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pack)           PACK="${2:?--pack needs a value}"; shift 2 ;;
        --btf-db)         BTF_DB="${2:?}"; shift 2 ;;
        --refs-dir)       REFS_DIR="${2:?}"; shift 2 ;;
        --out)            OUT="${2:?}"; shift 2 ;;
        --threads)        THREADS="${2:?}"; shift 2 ;;
        --sylph)          SYLPH="${2:?}"; shift 2 ;;
        --no-check-names) CHECK_NAMES=0; shift ;;
        --force)          FORCE=1; shift ;;
        --dry-run)        DRY_RUN=1; shift ;;
        -h|--help)        usage; exit 0 ;;
        *)                usage; die "unknown argument: $1" ;;
    esac
done

[[ "$THREADS" =~ ^[0-9]+$ ]] && (( THREADS > 0 )) || die "--threads must be a positive integer, got '${THREADS}'"

REFS_DIR="${REFS_DIR:-${BTF_DB}/packs/${PACK}/refs}"
OUT="${OUT:-${BTF_DB}/sylph/${PACK}_pack_refs.syldb}"

command -v "$SYLPH" >/dev/null 2>&1 || die "sylph not found as '${SYLPH}'. Install it (conda: bioconda::sylph=1.0.0) or pass --sylph."
[[ -d "$REFS_DIR" ]] || die "refs directory not found: ${REFS_DIR}. Run scripts/build_pack_refs.py first."

SYLPH_VERSION="$("$SYLPH" --version 2>&1 | head -1 | sed 's/^sylph //')"

# ---------------------------------------------------------------------------
# Collect inputs
# ---------------------------------------------------------------------------
mapfile -t FASTAS < <(find "$REFS_DIR" -maxdepth 1 -type f \
    \( -name '*.fna.gz' -o -name '*.fa.gz' -o -name '*.fasta.gz' \
       -o -name '*.fna' -o -name '*.fa' -o -name '*.fasta' \) | sort)
(( ${#FASTAS[@]} > 0 )) || die "no FASTA files in ${REFS_DIR}"

# GUARD 1: every basename must carry a taxid parse_screen can recover.
# Checked against the SAME rule parse_sylph uses - an all-digit token of >= 2 characters -
# rather than a prettier rule that might accept a name the parser then rejects.
BAD=()
for f in "${FASTAS[@]}"; do
    base="$(basename "$f")"
    [[ "$base" =~ ^[0-9]{2,}_ ]] || BAD+=("$base")
done
if (( ${#BAD[@]} > 0 )); then
    printf '[%s] %s of %s filenames do not start with <taxid>_ :\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${#BAD[@]}" "${#FASTAS[@]}" >&2
    printf '    %s\n' "${BAD[@]:0:10}" >&2
    (( ${#BAD[@]} > 10 )) && printf '    ... and %s more\n' "$(( ${#BAD[@]} - 10 ))" >&2
    if (( CHECK_NAMES == 1 )); then
        die "refusing to build. sylph reports FILENAMES and bin/parse_screen.py recovers the
       taxid from them, so these genomes could never map to a watchlist taxon: the arm
       would contribute nothing while the report still said 'no sylph/ANI support', which
       reads as a clean negative. Rename them <taxid>_<name>.<ext>, or pass
       --no-check-names if you genuinely want an unmappable index."
    fi
    log "WARNING: --no-check-names given; ${#BAD[@]} genomes will be unmappable."
fi

# GUARD 2: the stored PATH is tokenised too, not only the basename.
# parse_sylph splits the whole Genome_file string, so a numeric directory component
# injects a token indistinguishable from a taxid - "/media/36tb/..." contributes "36",
# which is a real NCBI taxid. Every sylph row would then carry phantom support for it.
PATH_TOKENS="$(printf '%s\n' "$REFS_DIR" | tr '/._' '\n\n\n' | grep -E '^[0-9]{2,}$' || true)"
if [[ -n "$PATH_TOKENS" ]]; then
    log "WARNING: the refs path '${REFS_DIR}' contains numeric component(s) that"
    log "         parse_sylph cannot tell from a taxid: $(printf '%s ' $PATH_TOKENS)"
    log "         Every hit in this database will appear to support those taxids."
    log "         Move the refs under a non-numeric path if any of them is on a watchlist."
fi

log "pack        : ${PACK}"
log "refs dir    : ${REFS_DIR}"
log "genomes     : ${#FASTAS[@]}"
log "output      : ${OUT}"
log "sylph       : ${SYLPH} (${SYLPH_VERSION})"
log "parameters  : c=${SYLPH_C} k=${SYLPH_K} min-spacing=${SYLPH_MIN_SPACING} (sylph defaults; MUST match 'sylph profile')"

if [[ -e "$OUT" ]] && (( FORCE == 0 )); then
    die "${OUT} exists. Pass --force to rebuild, or --out to write elsewhere."
fi

if (( DRY_RUN == 1 )); then
    printf '[dry-run] %q sketch -t %q -c %q -k %q --min-spacing %q -o %q <%s genomes>\n' \
        "$SYLPH" "$THREADS" "$SYLPH_C" "$SYLPH_K" "$SYLPH_MIN_SPACING" \
        "${OUT%.syldb}" "${#FASTAS[@]}" >&2
    log "[dry-run] nothing written"
    exit 0
fi

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
mkdir -p "$(dirname "$OUT")"
# Build to a temporary name and move it into place, so an interrupted run cannot leave a
# half-written .syldb where pack.yaml expects a usable one - the pipeline would load it and
# report an empty arm.
TMP_PREFIX="${OUT%.syldb}.building.$$"
cleanup() { rm -f "${TMP_PREFIX}.syldb"; }
trap 'cleanup' EXIT

log "sketching ${#FASTAS[@]} genomes ..."
"$SYLPH" sketch -t "$THREADS" -c "$SYLPH_C" -k "$SYLPH_K" \
    --min-spacing "$SYLPH_MIN_SPACING" -o "$TMP_PREFIX" "${FASTAS[@]}"

[[ -s "${TMP_PREFIX}.syldb" ]] || die "sylph sketch produced no ${TMP_PREFIX}.syldb"

# ---------------------------------------------------------------------------
# GUARD 3: verify the database we just wrote, before trusting it
# ---------------------------------------------------------------------------
INSPECT="$("$SYLPH" inspect "${TMP_PREFIX}.syldb" 2>/dev/null || true)"
got_c="$(printf '%s\n'  "$INSPECT" | sed -n 's/^[[:space:]]*c:[[:space:]]*\([0-9]\+\).*/\1/p' | head -1)"
got_k="$(printf '%s\n'  "$INSPECT" | sed -n 's/^[[:space:]]*k:[[:space:]]*\([0-9]\+\).*/\1/p' | head -1)"
got_ms="$(printf '%s\n' "$INSPECT" | sed -n 's/^[[:space:]]*min_spacing_parameter:[[:space:]]*\([0-9]\+\).*/\1/p' | head -1)"
got_n="$(printf '%s\n'  "$INSPECT" | grep -c '^[[:space:]]*- file_name:' || true)"

[[ -n "$got_c" ]] || die "could not read c from 'sylph inspect'; refusing to install an unverified database."
(( got_c == SYLPH_C ))   || die "built database reports c=${got_c}, expected ${SYLPH_C}. A database whose c differs from the query sketch matches nothing, silently."
(( got_k == SYLPH_K ))   || die "built database reports k=${got_k}, expected ${SYLPH_K}."
(( got_ms == SYLPH_MIN_SPACING )) || die "built database reports min_spacing=${got_ms}, expected ${SYLPH_MIN_SPACING}."
(( got_n == ${#FASTAS[@]} )) || die "built database holds ${got_n} genomes but ${#FASTAS[@]} were given. Some input was skipped; investigate rather than shipping a short index."

mv -f "${TMP_PREFIX}.syldb" "$OUT"
trap - EXIT

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
MANIFEST="$(dirname "$OUT")/BTF_BUILD_MANIFEST.txt"
{
    printf 'tool=sylph\n'
    printf 'built_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'db_file=%s\n' "$OUT"
    printf 'pack=%s\n' "$PACK"
    printf 'refs_dir=%s\n' "$REFS_DIR"
    printf 'genomes=%s\n' "${#FASTAS[@]}"
    printf 'sylph=%s\n' "$SYLPH"
    printf 'sylph_version=%s\n' "$SYLPH_VERSION"
    printf 'c=%s\nk=%s\nmin_spacing=%s\n' "$got_c" "$got_k" "$got_ms"
    printf 'threads=%s\n' "$THREADS"
    printf 'filename_taxid_check=%s\n' "$( (( CHECK_NAMES == 1 )) && echo enforced || echo SKIPPED )"
    # The namespace note matters for the same reason it does for Metabuli: a hit can only
    # be mapped back if the filename carries an NCBI taxid.
    printf 'taxid_namespace=ncbi (recovered from <taxid>_<name> filenames by parse_screen.py)\n'
    printf 'built_by=scripts/build_sylph_db.sh\n'
} > "$MANIFEST"

log "wrote ${OUT} ($(du -h "$OUT" | cut -f1), ${got_n} genomes)"
log "wrote ${MANIFEST}"
log "pack.yaml must point screen.sylph_db at this path."
