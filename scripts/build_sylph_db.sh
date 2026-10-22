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
# EXISTENCE IS NOT FRESHNESS (--check)
# ------------------------------------
# A sketch is built once and then goes on existing while the pack it describes keeps
# growing. Nothing noticed: on 2026-10-21 the installed sketch was from 2026-09-29 and held
# 359 genomes against 375 references, so the sixteen added since - including all twelve
# Filoviridae, whose whole purpose (D25) is catching a divergent filovirus no species-keyed
# row will catch - had been invisible to an ON-BY-DEFAULT arm for 22 days. No run failed.
# screen_summary.tsv recorded sylph as `ok` throughout, because the arm ran, the database
# loaded, and a genome that was never sketched simply produces no row - which is what a
# clean negative also looks like.
#
# D45 had already named this predicate for the competitive bundles ("existence is not
# freshness"); it had not travelled to the sketch. --check is that lesson here: it compares
# the installed database against the refs directory and exits non-zero when they disagree.
#
# WHERE THE NOTICE BELONGS, AND WHY THE SIDECAR EXISTS
# ----------------------------------------------------
# The only consumer of this sketch is a pipeline run, so a run is where the operator has to
# be told - not a schedule, which would report staleness to a log on days when nobody was
# going to use the database anyway. main.nf warns at launch.
#
# But a run cannot open a .syldb: sylph lives in the btf-sylph env and is used inside the
# process, so the head node has no sylph and can only compare timestamps. mtime cannot see
# a reference added, removed or swapped with every timestamp settled - the state that
# actually shipped. So every build writes <db>.genomes.txt, the names read back out of
# `sylph inspect`, and the run makes the exact comparison against that. --stamp writes it
# for a database that already exists.
#
#   scripts/build_sylph_db.sh                      # human pack, default paths
#   scripts/build_sylph_db.sh --pack human --force # rebuild over an existing index
#   scripts/build_sylph_db.sh --check              # read-only; is the sketch current?
#   scripts/build_sylph_db.sh --stamp              # record what the installed db holds
#   scripts/build_sylph_db.sh --dry-run            # print the plan and the command
#
# --check exit codes are distinct, so a caller can tell a verdict from a failure:
#   0 fresh    3 stale, rebuild with --force    1 could not tell (no sylph, bad db, ...)
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
CHECK_ONLY=0
STAMP_ONLY=0

log()  { printf '[%s] %s\n'        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; exit 1; }
mtime() { date -u -r "$1" +%Y-%m-%dT%H:%M:%SZ; }

# One parser for `sylph inspect`, used BOTH to verify a database this script just built and
# to judge one that is already installed (--check). Two copies of this sed would be two
# chances for the post-build verification and the freshness check to disagree about what the
# database says - and the whole point of either is that a parameter mismatch is silent.
inspect_num() {
    printf '%s\n' "$1" | sed -n "s/^[[:space:]]*$2:[[:space:]]*\([0-9]\+\).*/\1/p" | head -1
}
inspect_basenames() {
    printf '%s\n' "$1" |
        sed -n 's|^[[:space:]]*- file_name:[[:space:]]*\(.*\)$|\1|p' |
        sed 's|.*/||'
}

# THE SIDECAR EXISTS SO A RUN CAN ASK THE EXACT QUESTION WITHOUT SYLPH.
#
# sylph lives in the btf-sylph env and is used INSIDE the process, so the Nextflow head
# node cannot read the genome list out of the database. Without a record of what the sketch
# holds, a run can only compare timestamps - and mtime cannot see a reference added,
# removed or swapped with every timestamp settled, which is precisely the state that
# shipped for 22 days.
#
# This is D45's sidecar, one artifact over: the names are read back out of `sylph inspect`,
# so the file records what the database ACTUALLY HOLDS rather than what the build intended
# to put in it. A sidecar derived from the refs directory would be a claim nobody verified -
# the exact error this line of work exists to avoid - and it would agree with refs/ by
# construction, which is the one thing it must be able to disagree with.
sidecar_path() { printf '%s.genomes.txt\n' "$1"; }
write_sidecar() {
    local inspect="$1"
    local path tmp
    path="$(sidecar_path "$OUT")"
    tmp="${path}.tmp.$$"
    {
        printf '# genomes held by %s\n' "$(basename "$OUT")"
        printf '# read back from `sylph inspect`, not from the refs directory\n'
        printf '# written %s by scripts/build_sylph_db.sh\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        inspect_basenames "$inspect" | sort
    } > "$tmp"
    # Atomic, because main.nf reads this during a run: a reader must never see it half
    # written and conclude the sketch is missing genomes.
    mv -f "$tmp" "$path"
    log "wrote $(basename "$path") ($(inspect_basenames "$inspect" | grep -c . || true) genomes)"
}

# Prints up to eight members of a category and a count of the rest. Eight because a
# warning nobody can read to the end is a warning nobody reads.
show_some() {
    local label="$1"; shift
    printf '    %s (%s):\n' "$label" "$#" >&2
    printf '      %s\n' "${@:1:8}" >&2
    if (( $# > 8 )); then printf '      ... and %s more\n' "$(( $# - 8 ))" >&2; fi
}

usage() {
    # Bounded by the `set -` line rather than a line number, which silently truncated the
    # header as it grew: the range used to stop at 42 and had already swallowed the
    # --dry-run example.
    sed -n '2,/^set -/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'
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
  --check            compare the installed database against the refs directory and
                     exit 0 fresh / 3 stale / 1 undetermined. Writes nothing.
  --stamp            record what the installed database holds, read back from
                     'sylph inspect', so a RUN can make the exact comparison without
                     sylph. Writes only the sidecar; does not rebuild.
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
        --check)          CHECK_ONLY=1; shift ;;
        --stamp)          STAMP_ONLY=1; shift ;;
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

# ---------------------------------------------------------------------------
# --check: is the installed database a sketch of this pack, AS IT IS NOW?
# ---------------------------------------------------------------------------
# Three comparisons, because they fail in three different ways:
#
#   missing  a reference the sketch has never seen. The silent one, and the one that
#            actually happened for 22 days. The arm emits no row for that genome, which
#            is also what finding nothing looks like.
#   newer    a reference whose FILE is newer than the sketch. Its name IS in the sketch,
#            so a set comparison alone calls it present while the sequence sketched is
#            one the pack no longer carries: a re-download of one accession, a restore
#            from backup, anything that rewrites a file under a name already known.
#   extra    a genome in the sketch that is no longer in refs/. sylph reports filenames,
#            so this one is not silent but wrong in the other direction: screen support
#            for a reference stage 2 cannot competitively confirm, because it is not
#            there to map against.
#
# Parameters are checked too. A database at the wrong -c or -k matches nothing while
# loading cleanly, so for the purposes of this verdict that is a stale database as well:
# the fix is the same rebuild.
check_freshness() {
    if [[ ! -e "$OUT" ]]; then
        log "STALE: no database at ${OUT}"
        log "verdict: STALE. Build it: scripts/build_sylph_db.sh --pack ${PACK}"
        exit 3
    fi

    local inspect
    inspect="$("$SYLPH" inspect "$OUT" 2>/dev/null || true)"
    [[ -n "$inspect" ]] || die "'${SYLPH} inspect ${OUT}' produced no output, so freshness
       cannot be judged. sylph refuses any path not ending in .syldb, .sylsp or .syl2db, so
       a backup named <name>.syldb.bak-<date> is unreadable BY SYLPH and useless as a
       rollback; name it <name>.bak-<date>.syldb instead."

    local -a sketched=()
    local p
    while IFS= read -r p; do
        if [[ -n "$p" ]]; then sketched+=("$p"); fi
    done < <(inspect_basenames "$inspect")
    (( ${#sketched[@]} > 0 )) || die "${OUT} reports no genome files at all. It is damaged
       rather than stale; rebuild with --force."

    local got_c got_k got_ms
    got_c="$(inspect_num  "$inspect" 'c')"
    got_k="$(inspect_num  "$inspect" 'k')"
    got_ms="$(inspect_num "$inspect" 'min_spacing_parameter')"

    # Basenames on both sides. The stored file_name is an absolute path from build time, so
    # comparing paths would call the whole database stale after the pack moved volumes -
    # and a check that cries wolf is a check people learn to skip.
    local -A in_sketch=() in_refs=()
    local b f
    for b in "${sketched[@]}"; do in_sketch["$b"]=1; done

    local -a missing=() newer=() extra=()
    for f in "${FASTAS[@]}"; do
        b="$(basename "$f")"
        in_refs["$b"]=1
        if [[ -z "${in_sketch[$b]:-}" ]]; then
            missing+=("$b")
        elif [[ "$f" -nt "$OUT" ]]; then
            newer+=("$b")
        fi
    done
    for b in "${sketched[@]}"; do
        if [[ -z "${in_refs[$b]:-}" ]]; then extra+=("$b"); fi
    done

    log "sketch      : ${#sketched[@]} genomes, c=${got_c:-?} k=${got_k:-?} min-spacing=${got_ms:-?}, written $(mtime "$OUT")"
    log "refs        : ${#FASTAS[@]} FASTA files in ${REFS_DIR}"

    local stale=0
    if (( ${#missing[@]} > 0 )); then
        stale=1
        log "STALE: ${#missing[@]} reference(s) were never sketched. The sylph arm is BLIND"
        log "       to these genomes and reports nothing for them, which is indistinguishable"
        log "       from a clean negative in screen_summary.tsv."
        show_some "not in the sketch" "${missing[@]}"
    fi
    if (( ${#newer[@]} > 0 )); then
        stale=1
        log "STALE: ${#newer[@]} reference file(s) are NEWER than the sketch. The name is in"
        log "       the database but the sequence behind it has changed since."
        show_some "newer than the sketch" "${newer[@]}"
    fi
    if (( ${#extra[@]} > 0 )); then
        stale=1
        log "STALE: ${#extra[@]} sketched genome(s) are no longer in refs/. A hit on one is"
        log "       screen support for a reference stage 2 has nothing to map against."
        show_some "in the sketch, not in refs" "${extra[@]}"
    fi
    # FOURTH COMPARISON: the sidecar a RUN reads, against the database it describes.
    # A run cannot open the .syldb, so it trusts this file. A sidecar that disagrees with
    # the database is therefore worse than none - it would answer the exact question
    # wrongly and confidently - and the fix is --stamp, not a rebuild.
    local side
    side="$(sidecar_path "$OUT")"
    if [[ ! -e "$side" ]]; then
        log "NOTE: no $(basename "$side") beside the database, so a RUN can only compare"
        log "      timestamps and will not see a reference added, removed or swapped with"
        log "      every mtime settled. Write it from the database: --stamp"
    else
        local -a recorded=()
        while IFS= read -r b; do
            if [[ -n "$b" && "${b:0:1}" != "#" ]]; then recorded+=("$b"); fi
        done < "$side"
        local -A in_rec=()
        for b in "${recorded[@]}"; do in_rec["$b"]=1; done
        local -a drift=()
        for b in "${sketched[@]}"; do
            if [[ -z "${in_rec[$b]:-}" ]]; then drift+=("$b (in db, not recorded)"); fi
        done
        for b in "${recorded[@]}"; do
            if [[ -z "${in_sketch[$b]:-}" ]]; then drift+=("$b (recorded, not in db)"); fi
        done
        if (( ${#drift[@]} > 0 )); then
            stale=1
            log "STALE RECORD: $(basename "$side") disagrees with the database it describes."
            log "              A run reads that file and cannot open the .syldb, so it would"
            log "              answer the exact question wrongly. Fix with --stamp."
            show_some "disagreements" "${drift[@]}"
        fi
    fi

    if [[ "${got_c}" != "$SYLPH_C" || "${got_k}" != "$SYLPH_K" || "${got_ms}" != "$SYLPH_MIN_SPACING" ]]; then
        stale=1
        log "STALE: parameters are c=${got_c:-?} k=${got_k:-?} min-spacing=${got_ms:-?}, but"
        log "       'sylph profile' in modules/local/detect.nf queries at the defaults"
        log "       c=${SYLPH_C} k=${SYLPH_K} min-spacing=${SYLPH_MIN_SPACING}. A mismatch matches nothing, silently."
    fi

    if (( stale == 1 )); then
        log "verdict: STALE. Rebuild: scripts/build_sylph_db.sh --pack ${PACK} --force"
        exit 3
    fi
    log "verdict: FRESH. Every reference is sketched, none is newer than the sketch, no"
    log "         sketched genome has been removed, and the parameters match the query."
    exit 0
}

# --stamp: record what an ALREADY INSTALLED database holds, without rebuilding it.
#
# This is how the shipped sketch got a sidecar. It is legitimate where a stamp derived from
# the refs manifest would not be (D45 refused exactly that for the bundles): the names are
# read out of the artifact with `sylph inspect`, so this records a measurement of the
# database rather than an assumption about it.
stamp_sidecar() {
    [[ -e "$OUT" ]] || die "no database at ${OUT}; there is nothing to stamp. Build it first."
    local inspect
    inspect="$("$SYLPH" inspect "$OUT" 2>/dev/null || true)"
    [[ -n "$inspect" ]] || die "'${SYLPH} inspect ${OUT}' produced no output; refusing to
       write a record of a database that cannot be read."
    local n
    n="$(inspect_basenames "$inspect" | grep -c . || true)"
    (( n > 0 )) || die "${OUT} reports no genome files; refusing to record an empty set."
    write_sidecar "$inspect"
    exit 0
}

if (( STAMP_ONLY == 1 )); then
    stamp_sidecar       # always exits
fi

if (( CHECK_ONLY == 1 )); then
    check_freshness     # always exits
fi

if [[ -e "$OUT" ]] && (( FORCE == 0 )); then
    die "${OUT} exists. Pass --force to rebuild, or --out to write elsewhere.
       To ask whether it is still current rather than replacing it: --check"
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
got_c="$(inspect_num  "$INSPECT" 'c')"
got_k="$(inspect_num  "$INSPECT" 'k')"
got_ms="$(inspect_num "$INSPECT" 'min_spacing_parameter')"
got_n="$(inspect_basenames "$INSPECT" | grep -c . || true)"

[[ -n "$got_c" ]] || die "could not read c from 'sylph inspect'; refusing to install an unverified database."
(( got_c == SYLPH_C ))   || die "built database reports c=${got_c}, expected ${SYLPH_C}. A database whose c differs from the query sketch matches nothing, silently."
(( got_k == SYLPH_K ))   || die "built database reports k=${got_k}, expected ${SYLPH_K}."
(( got_ms == SYLPH_MIN_SPACING )) || die "built database reports min_spacing=${got_ms}, expected ${SYLPH_MIN_SPACING}."
(( got_n == ${#FASTAS[@]} )) || die "built database holds ${got_n} genomes but ${#FASTAS[@]} were given. Some input was skipped; investigate rather than shipping a short index."

mv -f "${TMP_PREFIX}.syldb" "$OUT"
trap - EXIT

write_sidecar "$INSPECT"

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
# NOTE: every sketch under ${BTF_DB}/sylph/ shares this one manifest, so building a second
# pack's sketch overwrites the first pack's provenance. `db_file=` below is what makes that
# recoverable, and anything reading `genomes=` must check it names the sketch it is judging.
MANIFEST="$(dirname "$OUT")/BTF_BUILD_MANIFEST.txt"
# The newest input, recorded so staleness has a watermark that does not need sylph to read:
# this timestamp against the refs directory's own newest file answers "has the pack moved
# on since?" with two stat calls.
NEWEST_REF="${FASTAS[0]}"
for f in "${FASTAS[@]}"; do
    if [[ "$f" -nt "$NEWEST_REF" ]]; then NEWEST_REF="$f"; fi
done
{
    printf 'tool=sylph\n'
    printf 'built_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'db_file=%s\n' "$OUT"
    printf 'pack=%s\n' "$PACK"
    printf 'refs_dir=%s\n' "$REFS_DIR"
    printf 'genomes=%s\n' "${#FASTAS[@]}"
    printf 'refs_newest=%s\n' "$(basename "$NEWEST_REF")"
    printf 'refs_newest_mtime=%s\n' "$(mtime "$NEWEST_REF")"
    printf 'freshness_check=scripts/build_sylph_db.sh --pack %s --check\n' "$PACK"
    printf 'genome_record=%s\n' "$(basename "$(sidecar_path "$OUT")")"
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
