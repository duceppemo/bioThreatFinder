#!/usr/bin/env bash
#
# build_metabuli_db.sh - obtain the Metabuli index that gives bioThreatFinder
#                        stage 1 a JOINT DNA + AMINO-ACID screen.
#
# Target path (must match screen.metabuli_db in packs/*/pack.yaml):
#   ${BTF_DB}/metabuli/<mode>              (default BTF_DB=/media/8tb/bioThreatFinder)
#
# Why Metabuli, given that Kraken2 is already the primary screen
# --------------------------------------------------------------
# Kraken2 fixes the OTHER half of the stage-1 problem. Its HyperLogLog
# distinct-minimizer counts (--report-minimizer-data) tell you whether a taxon's
# reads are spread over a genome or stacked on one conserved locus - that is a specificity fix, and it only
# helps for taxa whose exact nucleotide k-mers are in the index.
#
# The half it cannot fix: a divergent or genuinely novel agent shares almost no
# exact 31-mers with any reference, so it is not merely under-counted, it is
# invisible. Nucleotide k-mer matching fails silently there. Today the only tier
# in this pipeline that sees through divergence is the marker tier (protein and
# HMM evidence, stage 2b) - and that only works for the 163 curated markers, on
# reads that some earlier stage already flagged. A novel agent that stage 1 never
# raises never reaches it.
#
# Metabuli indexes "metamers": a DNA k-mer paired with its amino-acid translation,
# and matches on both. The amino-acid half is degenerate with respect to
# synonymous substitution, so a divergent orthologue still matches while the
# DNA half keeps the specificity that a protein-only search would lose. That is
# what puts a divergent agent on the stage-1 candidate list in the first place.
#
# So the three nucleotide/protein screens are complementary, not redundant:
#   kraken2     PRIMARY nucleotide screen; distinct-minimizer specificity,
#               abundance (Bracken)
#   krakenuniq  OPT-IN second opinion, off by default (docs/decisions.md D2)
#   metabuli    divergence and novelty; joint DNA + AA
#
# NVMe placement matters
# ----------------------
# Metabuli indexes are large (see the size table below) and querying them is
# dominated by reads into the index, not by CPU. The project DB root is on NVMe
# (/media/8tb) deliberately - putting a 400 GB Metabuli index on the 36 TB
# spinning volume would make stage 1 disk-bound and is the single easiest way to
# make this pipeline unusably slow. Metabuli also streams the index in chunks
# governed by --max-ram, so on a machine that cannot hold the index the seek
# pattern is exactly what NVMe is good at.
#
set -eEuo pipefail

# An ERR trap, because `set -e` alone exits SILENTLY.
#
# This build died mid-way through the multi-hour bacteria download and left a log whose
# last line was ordinary progress output - no error, no exit code, nothing to diagnose.
# `set -e` had done its job and aborted, but said nothing. For a job measured in hours
# that is the difference between "resume from the failing step" and "discover tomorrow
# that nothing happened".
#
# The trap needs `set -E` (errtrace), not just `set -e`: without it a bash ERR trap is
# NOT inherited by shell functions, and every step here runs inside one. The trap was
# added without -E, so the build kept dying silently anyway - kraken2-build --build
# failed, the script exited 1, and it printed nothing of its own.
_btf_on_err() {
    local rc=$? line=${1:-?}
    printf '[%s] ERROR: aborted at line %s with exit %s: %s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$line" "$rc" "${BASH_COMMAND:-unknown}" >&2
    printf '[%s] ERROR: completed steps are stamped, so re-running resumes rather than restarting.\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >&2
    exit "$rc"
}
trap '_btf_on_err $LINENO' ERR
# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
BTF_DB="/media/8tb/bioThreatFinder"
# Pinned dated snapshot, NOT the shared /media/36tb/db/taxdump that D8 measured
# at 201 days stale and decided to stop using.
TAXDUMP_SRC="/media/8tb/bioThreatFinder/taxdump/2026-09-29"
MAX_TAXDUMP_AGE_DAYS=90
THREADS=32

# gtdb  = GTDB species representatives. Best bacterial/archaeal resolution and
#         consistent taxonomy, but GTDB taxids are NOT NCBI taxids (see the long
#         comment at prepare_taxonomy) and the watchlist is resolved against NCBI.
# refseq= RefSeq-based, NCBI taxonomy, comparable in scope to Kraken2 PlusPF.
#         The DEFAULT, because it is the only mode whose taxids line up with
#         packs/*/watchlist.resolved.tsv without a translation layer.
# custom= build from your own FASTA list (pack reference bundles, curated
#         high-consequence genomes, a pathogen-only index for fast triage).
MODE="refseq"

ACCESSION2TAXID=""      # required for --mode custom (NCBI-style 4-column map)
FASTA_LIST=""           # required for --mode custom (one absolute path per line)
PREBUILT=""             # URL, or a known name, or "auto"
LOCAL_BUILD=0           # force a local build even for gtdb/refseq
MAX_RAM_GIB=0           # 0 = auto from /proc/meminfo
MIN_FREE_GIB=0          # 0 = derive from the route/mode, see check_free_space()
DRY_RUN=0
FORCE=0

# The conda prefix is not necessarily on PATH for non-interactive runs (cron, Nextflow),
# so look there explicitly before giving up. Overridable: this was a hardcoded personal
# miniconda prefix, and $CONDA_PREFIX covers the common case with no configuration.
CONDA_BIN="${BTF_CONDA_BIN:-${CONDA_PREFIX:-/opt/conda}/bin}"

# ---------------------------------------------------------------------------
# Prebuilt index catalogue
# ---------------------------------------------------------------------------
# PREBUILT IS THE RECOMMENDED ROUTE, and by a wide margin. Building a
# GTDB-scale Metabuli index locally means downloading ~100k+ genomes and then
# running the metamer extraction and indexing over all of them: days of wall
# time, a peak disk footprint several times the final index, and a result that is
# not reproducible against anyone else's because your genome set is your own.
# Downloading a published index is network-bound only, is the same artifact other
# groups cite, and is what the upstream authors validated their accuracy claims
# against. Build locally only for --mode custom, or when you specifically need a
# collection date or a genome set the published indexes do not cover.
#
# Published indexes and their stated sizes (VERIFY - these change with every
# release, and the size column is the figure upstream publishes for the archive,
# not necessarily the unpacked size):
#
#   name             stated size   taxonomy
#   ---------------- ------------  ------------------------------------------
#   refseq_standard  ~110-115 GB   NCBI   archaea, bacteria, viral, plasmid,
#                                         eukaryotes - scoped like Kraken2 PlusPF
#   gtdb226          ~375-380 GB   GTDB   species representatives
#   gtdb232          ~740-750 GB   GTDB   ~200k species reps + human + RefSeq viral
#   hrgm2            ~85 GB        NCBI   human gut reference genomes
#   hrom             ~40 GB        NCBI   human oral microbiome
#
# VERIFY: the URLs below were read from the Metabuli documentation's prebuilt
# database page and are recorded here as a STARTING POINT, not as fact. Upstream
# has already moved these between an S3 bucket and opendata.mmseqs.org at least
# once. Confirm the current location and the archive's checksum on the Metabuli
# docs before trusting a multi-hundred-GB download, and pass the confirmed URL
# with --prebuilt. A wrong URL costs you a 404; a stale URL silently gives you a
# months-old index, which is the failure mode that matters for biosurveillance.
#   index listing : https://opendata.mmseqs.org/metabuli
#   docs          : https://jaebeom-kim.github.io/metabuli-doc/
prebuilt_url_for() {
    # VERIFY every one of these. Also VERIFY the required minimum Metabuli
    # version - the newer indexes are documented as needing >= 1.2.0, and an
    # older binary reading a newer index is exactly the kind of silent
    # semantic mismatch this script refuses to gamble on.
    case "$1" in
        refseq|refseq_standard) printf '%s\n' "https://opendata.mmseqs.org/metabuli/refseq_standard.tar.gz" ;;
        gtdb|gtdb226)           printf '%s\n' "https://opendata.mmseqs.org/metabuli/gtdb226.tar.gz" ;;
        gtdb232)                printf '%s\n' "https://steineggerlab.s3.amazonaws.com/metabuli/gtdb232.tar.gz" ;;
        hrgm2)                  printf '%s\n' "https://opendata.mmseqs.org/metabuli/hrgm2.tar.gz" ;;
        hrom)                   printf '%s\n' "https://opendata.mmseqs.org/metabuli/hrom.tar.gz" ;;
        *)                      return 1 ;;
    esac
}

usage() {
    cat <<'EOF'
Usage:
  build_metabuli_db.sh [options]

Obtains the Metabuli index that gives bioThreatFinder stage 1 a joint
DNA + amino-acid screen, at
  <btf-db>/metabuli/<mode>

Metabuli matches "metamers" - a DNA k-mer together with its amino-acid
translation - so a divergent or novel agent that shares no exact nucleotide
k-mers with any reference is still seen at stage 1. Today only the marker tier
sees through divergence, and only for the 163 curated markers.

DOWNLOADING A PREBUILT INDEX IS THE RECOMMENDED ROUTE. A local GTDB-scale build
takes days and a peak disk footprint several times the final index; a prebuilt
index is network-bound only and is the same artifact upstream validated. Local
build is the fallback, and the only route for --mode custom.

Options:
  --btf-db DIR         Project-owned database root
                       (default: /media/8tb/bioThreatFinder)
  --taxdump DIR        Unpacked NCBI taxdump to COPY nodes/names/merged/
                       delnodes.dmp from (default: /media/36tb/db/taxdump,
                       shared and read-only). Local build only.
  --threads N          Threads (default: 32)
  --mode gtdb|refseq|custom
                       Which index to produce, and the subdirectory name.
                       refseq (DEFAULT) - NCBI taxonomy, PlusPF-like scope.
                         The only mode whose taxids line up with
                         packs/*/watchlist.resolved.tsv directly.
                       gtdb   - GTDB species representatives. Better bacterial
                         resolution, but GTDB taxids are NOT NCBI taxids and the
                         watchlist needs a translation layer first.
                       custom - build from your own FASTA list. Requires
                         --fasta-list and --accession2taxid.
  --accession2taxid FILE
                       NCBI-style accession2taxid map (custom mode; required).
  --fasta-list FILE    Text file, one ABSOLUTE path to a reference FASTA per
                       line (custom mode; required).
  --prebuilt URL|NAME|auto
                       Download a prebuilt index instead of building.
                       A full https:// URL is used as given (preferred - you
                       have verified it). A NAME (refseq_standard, gtdb226,
                       gtdb232, hrgm2, hrom) or "auto" uses this script's
                       recorded URL for that name, which is marked VERIFY in
                       the source and must be confirmed against the Metabuli
                       docs before you trust it.
                       For --mode refseq|gtdb this is the DEFAULT route.
  --local-build        Build locally even for --mode refseq|gtdb (the fallback).
  --max-ram GIB        Cap Metabuli's memory use (default: ~70% of this box's
                       RAM). Metabuli's own default is 128 GB, which leaves most
                       of a 503 GB box idle during a build.
  --min-free-gib N     Override the pre-flight free-space requirement. Default is
                       derived from the chosen route and mode; the check prints
                       its breakdown.
  --force              Redo steps even when their outputs already exist
  --dry-run            Print the plan and the exact commands, download nothing
  -h, --help           This message

What a local build consumes
---------------------------
  metabuli build <DBDIR> <FASTA_LIST> <accession2taxid> --taxonomy-path <TAXDUMP>
GTDB variant:
  metabuli build --gtdb 1 <DBDIR> <FASTA_LIST> <GTDB_TAXDUMP/taxid.map> \
      --taxonomy-path <GTDB_TAXDUMP>
  FASTA_LIST       plain text, one ABSOLUTE path to a reference FASTA per line
  accession2taxid  NCBI 4-column format: accession, accession.version, taxid, gi
                   (tab-separated, with the NCBI header line). Sequences whose
                   accession is absent from this map are SILENTLY SKIPPED, which
                   is why this script counts what it fed in and what came out.
  --taxonomy-path  directory holding nodes.dmp / names.dmp (+ merged/delnodes)
VERIFY: argument order, the positional-vs-flag split, and the flag spellings
(--taxonomy-path, --gtdb, --make-library, --max-ram, --accession-level) against
`metabuli build -h` for YOUR version before a multi-day run. A flag that has
been renamed fails loudly; a flag whose MEANING changed does not, and this
script would rather stop and ask.

Nothing under /media/36tb/db is ever written.
EOF
}

# ---------------------------------------------------------------------------
# Logging - same vocabulary as build_kraken2_db.sh / build_krakenuniq_db.sh
# ---------------------------------------------------------------------------
log()  { printf '[%s] %s\n'          "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
warn() { printf '[%s] WARNING: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n'   "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; exit 1; }
skip() { printf '[%s] SKIP: %s\n'    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2; }
run()  {
    if (( DRY_RUN == 1 )); then
        printf '[dry-run] %s\n' "$*" >&2
    else
        log "+ $*"
        "$@"
    fi
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --btf-db)          BTF_DB="${2:?--btf-db needs a value}"; shift 2 ;;
        --max-taxdump-age-days) MAX_TAXDUMP_AGE_DAYS="${2:?}"; shift 2 ;;
        --taxdump)         TAXDUMP_SRC="${2:?--taxdump needs a value}"; shift 2 ;;
        --threads)         THREADS="${2:?--threads needs a value}"; shift 2 ;;
        --mode)            MODE="${2:?--mode needs a value}"; shift 2 ;;
        --accession2taxid) ACCESSION2TAXID="${2:?--accession2taxid needs a value}"; shift 2 ;;
        --fasta-list)      FASTA_LIST="${2:?--fasta-list needs a value}"; shift 2 ;;
        --prebuilt)        PREBUILT="${2:?--prebuilt needs a value}"; shift 2 ;;
        --local-build)     LOCAL_BUILD=1; shift ;;
        --max-ram)         MAX_RAM_GIB="${2:?--max-ram needs a value}"; shift 2 ;;
        --min-free-gib)    MIN_FREE_GIB="${2:?--min-free-gib needs a value}"; shift 2 ;;
        --force)           FORCE=1; shift ;;
        --dry-run)         DRY_RUN=1; shift ;;
        -h|--help)         usage; exit 0 ;;
        *)                 usage; die "unknown option: $1" ;;
    esac
done

case "$MODE" in
    gtdb|refseq|custom) : ;;
    *) die "--mode must be one of gtdb, refseq, custom (got '${MODE}')" ;;
esac
[[ "$THREADS"      =~ ^[0-9]+$ ]] || die "--threads must be an integer"
[[ "$MAX_RAM_GIB"  =~ ^[0-9]+$ ]] || die "--max-ram must be an integer (GiB)"
[[ "$MIN_FREE_GIB" =~ ^[0-9]+$ ]] || die "--min-free-gib must be an integer"
(( THREADS >= 1 )) || die "--threads must be >= 1"

DB_DIR="${BTF_DB}/metabuli/${MODE}"
STAMP_DIR="${DB_DIR}/.build_stamps"
TAX_DIR="${DB_DIR}/taxonomy"
DL_DIR="${DB_DIR}/.download"

# Hard guard, same posture as build_kraken2_db.sh: a mistyped --btf-db must not
# be able to write anywhere else, and must never reach the shared volume.
case "$DB_DIR" in
    */metabuli/gtdb|*/metabuli/refseq|*/metabuli/custom) : ;;
    *) die "refusing to build at '${DB_DIR}' - expected a path ending in /metabuli/{gtdb,refseq,custom}" ;;
esac
[[ "$DB_DIR" == "${BTF_DB}/metabuli/"* ]] \
    || die "refusing to build at '${DB_DIR}' - it is not under ${BTF_DB}/metabuli"
case "$DB_DIR" in
    /media/36tb/db/*) die "refusing to write under /media/36tb/db - that volume is shared and read-only for this project" ;;
esac
for protected in "$TAXDUMP_SRC" \
                 "/media/36tb/db/taxdump" \
                 "/media/36tb/db/kraken2" \
                 "/media/36tb/db/nt" \
                 "/media/36tb/db/sylph"; do
    if [[ "$DB_DIR" == "$protected" || "$DB_DIR" == "$protected"/* ]]; then
        die "refusing to write inside the pre-existing shared resource ${protected}"
    fi
done

# ---------------------------------------------------------------------------
# Route selection: prebuilt (recommended) vs local build (fallback)
# ---------------------------------------------------------------------------
ROUTE="prebuilt"
PREBUILT_URL=""
if [[ "$MODE" == "custom" ]]; then
    # There is nothing published to download for a private genome set.
    if [[ -n "$PREBUILT" ]]; then
        die "--prebuilt makes no sense with --mode custom (there is no published index for your own genome set)"
    fi
    ROUTE="build"
elif (( LOCAL_BUILD == 1 )); then
    ROUTE="build"
    warn "--local-build: taking the FALLBACK route for --mode ${MODE}."
    warn "  A prebuilt index would be network-bound only, days faster, and is the"
    warn "  artifact upstream validated. Continue only if you specifically need a"
    warn "  different collection date or genome set."
else
    case "$PREBUILT" in
        ""|auto) PREBUILT_URL="$(prebuilt_url_for "$MODE" || true)" ;;
        http://*|https://*) PREBUILT_URL="$PREBUILT" ;;
        *) PREBUILT_URL="$(prebuilt_url_for "$PREBUILT" || true)"
           [[ -n "$PREBUILT_URL" ]] \
               || die "unknown prebuilt name '${PREBUILT}'. Known names: refseq_standard gtdb226 gtdb232 hrgm2 hrom.
       Or pass the full https:// URL you have verified against
         https://jaebeom-kim.github.io/metabuli-doc/" ;;
    esac
    [[ -n "$PREBUILT_URL" ]] \
        || die "no prebuilt URL recorded for --mode ${MODE}. Pass --prebuilt <URL> or --local-build."
fi

if [[ "$ROUTE" == "build" ]]; then
    if [[ "$MODE" == "custom" ]]; then
        [[ -n "$FASTA_LIST" ]] \
            || die "--mode custom needs --fasta-list FILE (one ABSOLUTE path to a reference FASTA per line)"
        [[ -n "$ACCESSION2TAXID" ]] \
            || die "--mode custom needs --accession2taxid FILE (NCBI 4-column accession2taxid map).
       Metabuli SILENTLY SKIPS any sequence whose accession is not in that map,
       so building without it does not fail - it quietly produces a smaller index."
        [[ -s "$FASTA_LIST" ]]      || die "--fasta-list not found or empty: ${FASTA_LIST}"
        [[ -s "$ACCESSION2TAXID" ]] || die "--accession2taxid not found or empty: ${ACCESSION2TAXID}"
    else
        # A local refseq/gtdb build still needs a genome set and a map; this script
        # does not fetch 100k+ genomes for you, and pretending otherwise would be
        # the "silently skip a step" failure this project forbids.
        [[ -n "$FASTA_LIST" && -n "$ACCESSION2TAXID" ]] || die \
"a local --mode ${MODE} build needs both --fasta-list and --accession2taxid.
       This script deliberately does NOT fetch the genome set for you: assembling
       a reproducible RefSeq/GTDB genome collection is a separate, auditable step
       (ncbi-datasets-cli, or the GTDB representative genome tarball), and burying
       it here would make the index's provenance unrecoverable.
       Either:
         * take the recommended route and drop --local-build (downloads a
           published index; see --help), or
         * assemble the genome set, write the FASTA list and the
           accession2taxid map, and pass both.
       For GTDB, the taxonkit-derived GTDB taxdump and its taxid.map are the
       accepted source: https://github.com/shenwei356/gtdb-taxdump/releases
       VERIFY that URL and the release you pin."
        [[ -s "$FASTA_LIST" ]]      || die "--fasta-list not found or empty: ${FASTA_LIST}"
        [[ -s "$ACCESSION2TAXID" ]] || die "--accession2taxid not found or empty: ${ACCESSION2TAXID}"
    fi
fi

# ---------------------------------------------------------------------------
# Tool resolution - PATH first, then the conda bin, never a silent skip
# ---------------------------------------------------------------------------
resolve_tool() {
    local name="$1" p
    if p="$(command -v "$name" 2>/dev/null)"; then printf '%s\n' "$p"; return 0; fi
    if [[ -x "${CONDA_BIN}/${name}" ]]; then printf '%s\n' "${CONDA_BIN}/${name}"; return 0; fi
    return 1
}

INSTALL_HINT="Install Metabuli into its own environment:
         ${CONDA_BIN}/mamba create -y -n metabuli -c bioconda -c conda-forge metabuli
         conda activate metabuli
       VERIFY the version you get: the current prebuilt indexes are documented as
       requiring Metabuli >= 1.2.0, and an older binary reading a newer index is a
       silent semantic mismatch, not a clean error."

METABULI="$(resolve_tool metabuli || true)"
if [[ -z "$METABULI" ]]; then
    msg="metabuli not found on PATH or in ${CONDA_BIN}.
       ${INSTALL_HINT}"
    # Under --dry-run the point is to show the plan, so carry on with the bare
    # command name instead of aborting - same behaviour as build_kraken2_db.sh.
    if (( DRY_RUN == 1 )); then warn "$msg"; METABULI="metabuli"; else die "$msg"; fi
fi

DOWNLOADER=""
if [[ "$ROUTE" == "prebuilt" ]]; then
    # curl is preferred over wget only because -C - resume semantics are the ones
    # this script relies on for a multi-hundred-GB archive.
    if resolve_tool curl >/dev/null; then DOWNLOADER="curl"
    elif resolve_tool wget >/dev/null; then DOWNLOADER="wget"
    else
        msg="neither curl nor wget found; cannot download a prebuilt index.
       Install one, or take the local-build route with --local-build."
        if (( DRY_RUN == 1 )); then warn "$msg"; DOWNLOADER="curl"; else die "$msg"; fi
    fi
    command -v tar >/dev/null 2>&1 || die "tar not found; cannot unpack a prebuilt index"
fi

if (( MAX_RAM_GIB == 0 )); then
    # Metabuli's own default is 128 GB. On a 503 GB box that throttles the build
    # and the query for no reason. 70% leaves room for the rest of the pipeline,
    # which on this shared machine is not hypothetical.
    MAX_RAM_GIB="$(awk '/^MemTotal:/ {printf "%d", ($2/1048576)*0.70}' /proc/meminfo)"
fi

# ---------------------------------------------------------------------------
# Pre-flight free space
# ---------------------------------------------------------------------------
# The arithmetic, so the number is arguable instead of magic. All GiB, and all
# based on the STATED archive sizes above, which are marked VERIFY.
#
# PREBUILT route, peak = archive + unpacked index:
#   refseq_standard  ~112 archive + ~130 unpacked  -> ~250, resident ~130
#   gtdb226          ~378 archive + ~420 unpacked  -> ~800, resident ~420
#   gtdb232          ~745 archive + ~820 unpacked  -> ~1600, resident ~820
# The archive is deleted after a verified unpack, so steady state is the second
# number. A resumable download plus a separate unpack is chosen over streaming
# `curl | tar -xz` on purpose: streaming halves the peak, but a network blip at
# 600 GB then costs the whole download again. Disk is cheap here (5.9 TB free);
# a repeated 745 GB transfer is not. To trade the other way, stream it by hand:
#   curl -fSL <url> | tar -xzf - -C <db-dir>
#
# LOCAL BUILD route, peak:
#   genome set FASTA               ~150 (refseq scope) to ~700 (GTDB scope)
#   metamer extraction scratch     ~1-2x the final index
#   final index                    ~130 (refseq scope) to ~800 (GTDB scope)
# -> ~600 for a refseq-scope build, ~2600 for a GTDB-scope build. This is the
# concrete reason prebuilt is the recommendation, not a stylistic preference.
default_min_free() {
    if [[ "$ROUTE" == "prebuilt" ]]; then
        case "$PREBUILT_URL" in
            *gtdb232*) printf '1600\n' ;;
            *gtdb*)    printf '800\n'  ;;
            *hrgm2*)   printf '200\n'  ;;
            *hrom*)    printf '120\n'  ;;
            *)         printf '300\n'  ;;   # refseq_standard, or an unrecognised URL
        esac
    else
        case "$MODE" in
            gtdb)   printf '2600\n' ;;
            refseq) printf '600\n'  ;;
            *)      printf '400\n'  ;;   # custom: usually far smaller, but not always
        esac
    fi
}

check_free_space() {
    local probe="$BTF_DB" avail_gib
    while [[ ! -d "$probe" && "$probe" != "/" ]]; do probe="$(dirname "$probe")"; done
    avail_gib="$(df -BG --output=avail "$probe" | tail -1 | tr -dc '0-9')"
    [[ -n "$avail_gib" ]] || die "could not determine free space for ${probe}"
    log "free space under ${probe}: ${avail_gib} GiB"
    log "required (transient peak): ${MIN_FREE_GIB} GiB   [route: ${ROUTE}, mode: ${MODE}]"
    if [[ "$ROUTE" == "prebuilt" ]]; then
        log "  breakdown: archive .tar.gz + the unpacked index coexist until the"
        log "             archive is verified and deleted, so the peak is roughly"
        log "             2.2x the published archive size. Steady state is ~1.15x."
        log "             Published sizes (VERIFY): refseq_standard ~112 GB,"
        log "             gtdb226 ~378 GB, gtdb232 ~745 GB, hrgm2 ~85 GB, hrom ~42 GB."
    else
        log "  breakdown: ~150-700 GiB genome-set FASTA,"
        log "             ~1-2x the final index in metamer-extraction scratch,"
        log "             ~130-800 GiB final index, ~1 GiB copied taxonomy."
        log "             A GTDB-scope local build peaks near 2600 GiB and runs for"
        log "             days - which is why the prebuilt route is the default."
    fi
    if (( avail_gib < MIN_FREE_GIB )); then
        die "insufficient free space: ${avail_gib} GiB available, ${MIN_FREE_GIB} GiB required.
       Options:
         * choose a smaller published index (refseq_standard is ~1/7 of gtdb232
           and, unlike the GTDB variants, uses NCBI taxids the watchlist already
           resolves against), or
         * stream the archive instead of staging it (halves the peak, loses
           resume):  curl -fSL <url> | tar -xzf - -C ${DB_DIR}
         * free space on the volume, or
         * lower the requirement with --min-free-gib if you are certain.
       Nothing was downloaded."
    fi
}

# ---------------------------------------------------------------------------
# Stamps - same resume pattern as build_kraken2_db.sh
# ---------------------------------------------------------------------------
stamp()      { printf '%s\n' "${STAMP_DIR}/$1.done"; }
have_stamp() { [[ -f "$(stamp "$1")" ]]; }
set_stamp()  { (( DRY_RUN == 1 )) || { mkdir -p "$STAMP_DIR"; date -u +%Y-%m-%dT%H:%M:%SZ > "$(stamp "$1")"; }; }

# ---------------------------------------------------------------------------
# Taxonomy (local build only)
# ---------------------------------------------------------------------------
# `metabuli build --taxonomy-path DIR` expects a directory holding the NCBI dmp
# files. Copy them in rather than pointing at /media/36tb/db/taxdump directly,
# for the same reason build_kraken2_db.sh does: a build tool that writes anything
# into its taxonomy directory - a derived map, a merged copy, an unpacked
# refresh - would corrupt the one taxdump that resolve_watchlist.py pins for
# every pack, and a corrupted taxdump does not announce itself. It shows up later
# as watchlist entries that stop resolving, i.e. taxa that silently stop being
# screened for. Four files, about a gigabyte, and that entire class of accident
# is gone. Copy, never symlink - a symlink would be followed straight back to the
# shared volume.
#
# VERIFY: the flag spelling (--taxonomy-path), and whether the version in use
# also wants merged.dmp/delnodes.dmp or only nodes.dmp/names.dmp. They are cheap,
# so all four are copied regardless.
#
# GTDB WARNING that outlives this script: a GTDB-based index does NOT use NCBI
# taxids. packs/*/watchlist.resolved.tsv is resolved against the NCBI taxdump, so
# a GTDB index's taxids will not match the watchlist and the screen will appear to
# find nothing. --mode gtdb therefore needs the taxonkit GTDB taxdump AND a
# GTDB->NCBI translation step before its output can be joined to a watchlist.
# That translation does not exist in this pipeline yet, which is why --mode refseq
# is the default.
# The taxonomy a screen is built against decides whether its taxids join the
# watchlist. decisions.md D8 measured the shared /media/36tb/db/taxdump at 201 days
# stale (+288,764 taxids of drift) and decided to stop using it; this script defaulted
# to it anyway and had no age check, so a build could bake a 201-day-old taxonomy into
# an index while the watchlist resolved against the pinned snapshot. Taxa added or
# merged in between then silently fail to join - partial, invisible false negatives.
require_fresh_taxdump() {
    local src="$1" max="${2:-90}" age
    [[ -s "${src}/nodes.dmp" ]] || return 0
    age=$(( ( $(date -u +%s) - $(stat -c %Y "${src}/nodes.dmp") ) / 86400 ))
    log "taxdump ${src} is ${age} day(s) old"
    if (( age > max )); then
        die "taxdump ${src} is ${age} days old (limit ${max}).
  Reference sequences are current, so an old taxonomy drops or misassigns every
  organism added since the snapshot, and the index then disagrees with the
  watchlist - silently, as zero matches rather than an error.
  Refresh from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/ and pass
  --taxdump, or pass --max-taxdump-age-days to accept the risk deliberately."
    fi
}

prepare_taxonomy() {
    if have_stamp taxonomy && (( FORCE == 0 )); then
        skip "taxonomy already prepared (${TAX_DIR})"
        return 0
    fi
    [[ -d "$TAXDUMP_SRC" ]] || die "--taxdump directory does not exist: ${TAXDUMP_SRC}"
    for f in nodes.dmp names.dmp; do
        [[ -s "${TAXDUMP_SRC}/${f}" ]] \
            || die "required taxonomy file missing or empty: ${TAXDUMP_SRC}/${f}"
    done
    require_fresh_taxdump "$TAXDUMP_SRC" "$MAX_TAXDUMP_AGE_DAYS"
    log "copying taxonomy from ${TAXDUMP_SRC} (copy, never symlink - see comment)"
    run mkdir -p "$TAX_DIR"
    run cp -f "${TAXDUMP_SRC}/nodes.dmp" "${TAX_DIR}/nodes.dmp"
    run cp -f "${TAXDUMP_SRC}/names.dmp" "${TAX_DIR}/names.dmp"
    for extra in merged.dmp delnodes.dmp; do
        if [[ -s "${TAXDUMP_SRC}/${extra}" ]]; then
            run cp -f "${TAXDUMP_SRC}/${extra}" "${TAX_DIR}/${extra}"
        else
            warn "${TAXDUMP_SRC}/${extra} not found; metabuli may want it, and the provenance record will be incomplete"
        fi
    done
    if (( DRY_RUN == 0 )); then
        printf 'taxdump_source=%s\ntaxdump_copied_utc=%s\n' \
            "$TAXDUMP_SRC" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "${TAX_DIR}/BTF_PROVENANCE.txt"
    fi
    set_stamp taxonomy
}

# ---------------------------------------------------------------------------
# Prebuilt route
# ---------------------------------------------------------------------------
fetch_prebuilt() {
    local archive="${DL_DIR}/$(basename "${PREBUILT_URL%%\?*}")"
    if have_stamp "unpack" && (( FORCE == 0 )); then
        skip "prebuilt index already unpacked into ${DB_DIR}"
        return 0
    fi

    warn "about to download a prebuilt Metabuli index. This is HUNDREDS OF GIGABYTES"
    warn "  and takes many hours. Run under nohup, tmux or screen."
    warn "  URL: ${PREBUILT_URL}"
    warn "  VERIFY this URL against https://jaebeom-kim.github.io/metabuli-doc/"
    warn "  before trusting it - upstream has moved these archives at least once,"
    warn "  and a stale URL gives you a months-old index without any error."

    run mkdir -p "$DL_DIR"
    if ! have_stamp "download" || (( FORCE == 1 )); then
        # -C - resumes a partial transfer, which is the entire reason the archive
        # is staged rather than streamed. --retry rides out transient 5xx/timeouts.
        case "$DOWNLOADER" in
            curl) run curl -fSL --retry 10 --retry-delay 30 --retry-all-errors \
                          -C - -o "$archive" "$PREBUILT_URL" ;;
            wget) run wget --continue --tries=10 --waitretry=30 \
                          -O "$archive" "$PREBUILT_URL" ;;
        esac
        if (( DRY_RUN == 0 )); then
            [[ -s "$archive" ]] || die "download produced no file at ${archive}"
            log "downloaded $(du -h "$archive" | cut -f1)"
        fi
        set_stamp "download"
    else
        skip "archive already downloaded: ${archive}"
    fi

    # VERIFY: upstream does not consistently publish checksums beside these
    # archives. If a .md5 or .sha256 exists for the URL you used, check it here -
    # a truncated multi-hundred-GB index that unpacks "successfully" is a real
    # possibility and it would degrade the screen silently.
    if (( DRY_RUN == 0 )); then
        log "verifying the archive is intact before unpacking (gzip -t over $(du -h "$archive" | cut -f1); this reads the whole file)"
        gzip -t "$archive" 2>/dev/null \
            || die "archive at ${archive} is not a valid gzip stream - the download is truncated or corrupt.
       Delete ${STAMP_DIR}/download.done and the archive, then rerun to re-fetch."
    fi

    log "unpacking into ${DB_DIR}"
    # --strip-components is deliberately NOT guessed: these archives may or may
    # not carry a top-level directory, and stripping the wrong number of levels
    # scatters the index. Unpack as published, then report what landed where.
    run tar -xzf "$archive" -C "$DB_DIR"
    set_stamp "unpack"

    if (( DRY_RUN == 0 )); then
        log "removing the archive to reclaim $(du -h "$archive" | cut -f1)"
        run rm -f "$archive"
        run rmdir --ignore-fail-on-non-empty "$DL_DIR"
        log "contents of ${DB_DIR}:"
        find "$DB_DIR" -maxdepth 2 -mindepth 1 -not -name '.build_stamps' \
            -printf '    %y %10s  %P\n' 2>/dev/null | head -40 >&2 || true
        warn "If the index landed inside a nested subdirectory, point"
        warn "  screen.metabuli_db at THAT directory, not at ${DB_DIR}."
    fi
}

# ---------------------------------------------------------------------------
# Local build route
# ---------------------------------------------------------------------------
build_local() {
    if have_stamp "build" && (( FORCE == 0 )); then
        skip "index already built in ${DB_DIR}"
        return 0
    fi

    local n_fasta missing=0 line
    if (( DRY_RUN == 0 )); then
        n_fasta="$(grep -cve '^[[:space:]]*$' "$FASTA_LIST" || true)"
        (( n_fasta > 0 )) || die "--fasta-list ${FASTA_LIST} has no entries"
        # Metabuli requires ABSOLUTE paths in the list, and a path it cannot open
        # is one genome silently missing from the index. Check before spending days.
        while IFS= read -r line; do
            [[ -z "${line// }" ]] && continue
            case "$line" in
                /*) [[ -s "$line" ]] || { warn "listed FASTA missing or empty: ${line}"; missing=$((missing+1)); } ;;
                *)  warn "listed FASTA path is not absolute: ${line}"; missing=$((missing+1)) ;;
            esac
        done < "$FASTA_LIST"
        (( missing == 0 )) || die "${missing} of ${n_fasta} entries in ${FASTA_LIST} are unusable.
       Fix the list rather than building an index that is quietly incomplete."
        log "FASTA list: ${n_fasta} genomes, all present and absolute"
    fi

    local args=(build)
    if [[ "$MODE" == "gtdb" ]]; then
        # VERIFY: --gtdb 1 plus a GTDB taxid.map as the third positional, per the
        # Metabuli docs' GTDB recipe. Confirm against `metabuli build -h`.
        args+=(--gtdb 1)
    fi
    args+=("$DB_DIR" "$FASTA_LIST" "$ACCESSION2TAXID"
           --taxonomy-path "$TAX_DIR"
           --threads "$THREADS"
           --max-ram "$MAX_RAM_GIB")
    # VERIFY: --make-library builds a per-species genome library before indexing
    # and is documented as strongly recommended when one FASTA holds several
    # species. It is NOT enabled by default here because it costs another full
    # copy of the sequence on disk, and for one-genome-per-file reference bundles
    # (how build_pack_refs.py writes them) it buys nothing. Add it by hand for a
    # multi-species-per-file collection.

    warn "the metamer extraction and indexing phases now run for many hours (days at GTDB scale)"
    run "$METABULI" "${args[@]}"

    if (( DRY_RUN == 0 )); then
        # Do not assert specific filenames - they have changed between Metabuli
        # versions and a hardcoded list would fail on a working index. Assert that
        # something substantial exists instead.
        local produced
        produced="$(find "$DB_DIR" -maxdepth 2 -type f -size +1M -not -path "*/.build_stamps/*" -not -path "*/taxonomy/*" | wc -l)"
        (( produced > 0 )) || die "build finished but ${DB_DIR} holds no index files - do not point the pipeline at this directory"
        log "build produced ${produced} index file(s) >1 MiB"
    fi
    set_stamp "build"
}

detect_taxonomy_namespace() {
    # ASK THE ARTIFACT, NOT THE FLAG.
    #
    # This used to be `[[ "$MODE" == gtdb ]] && echo GTDB || echo NCBI`: an assertion
    # about how the script was invoked, not a check on what was written. The prebuilt
    # route accepts `--prebuilt gtdb232` INDEPENDENTLY of --mode, and DB_DIR is
    # ${BTF_DB}/metabuli/${MODE} - so `build_metabuli_db.sh --prebuilt gtdb232`
    # unpacked a GTDB index into .../metabuli/refseq, exactly the path pack.yaml points
    # at, and stamped it taxonomy_namespace=NCBI. main.nf's assertNcbiTaxonomy() then
    # read that line, saw NCBI and let the run proceed; every Metabuli taxid failed to
    # join the watchlist and the entire watchlist came back a clean negative. That is
    # the very failure the guard exists to catch, defeated by trusting a flag.
    #
    # GTDB taxonomy carries rank-prefixed names (d__Bacteria, s__Escherichia coli),
    # which NCBI never does, so the names file itself is the evidence.
    local d="${1:-}" f names="" requested_gtdb=0
    case "${MODE,,}|${PREBUILT,,}" in *gtdb*) requested_gtdb=1 ;; esac
    # The prebuilt archives ship their taxonomy under ncbi-taxdump/, which this list did
    # not know about - so a correctly unpacked NCBI index was stamped
    # "NCBI (UNVERIFIED - no taxonomy names file found)" and main.nf:171 then refused to
    # run on it. Searching for the FILE rather than guessing its directory also covers
    # whatever the next archive layout calls it.
    for f in "${d}/taxonomy/names.dmp" "${d}/names.dmp" "${d}/ncbi-taxdump/names.dmp" \
             "${d}/taxonomy/nodes.dmp"; do
        [[ -s "$f" ]] && { names="$f"; break; }
    done
    if [[ -z "$names" ]]; then
        names="$(find "$d" -maxdepth 2 -name names.dmp -size +0 -print -quit 2>/dev/null || true)"
    fi
    if [[ -n "$names" ]]; then
        if head -c 2000000 "$names" 2>/dev/null \
                | grep -qE '(^|[|[:space:]])[dpcofgs]__'; then
            echo "GTDB (NOT NCBI - needs translation before joining watchlist.resolved.tsv)"
            return 0
        fi
        # Absence of rank prefixes in a real names file is positive evidence of NCBI,
        # and it OUTRANKS what was requested. If the two disagree, say so loudly: one
        # of them is wrong and a reader needs to find out which before trusting a run.
        if (( requested_gtdb )); then
            echo "NCBI (CONTRADICTION: a GTDB collection was requested but ${names} carries no rank-prefixed names. Verify the index before use.)"
        else
            echo "NCBI (verified against ${names})"
        fi
        return 0
    fi
    # No artifact evidence either way: fall back to what was requested, and say the
    # answer is unverified rather than implying it was checked.
    if (( requested_gtdb )); then
        echo "GTDB (NOT NCBI - needs translation before joining watchlist.resolved.tsv)"
    else
        echo "NCBI (UNVERIFIED - no taxonomy names file found under ${d} to confirm)"
    fi
}

write_manifest() {
    local out="${DB_DIR}/BTF_BUILD_MANIFEST.txt"
    if (( DRY_RUN == 1 )); then
        printf '[dry-run] write %s\n' "$out" >&2
        return 0
    fi
    {
        printf 'tool=metabuli\n'
        printf 'built_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'db_dir=%s\n' "$DB_DIR"
        printf 'mode=%s\nroute=%s\n' "$MODE" "$ROUTE"
        printf 'metabuli=%s\n' "$METABULI"
        printf 'metabuli_version=%s\n' "$("$METABULI" version 2>&1 | head -1 || echo unknown)"
        printf 'threads=%s\nmax_ram_gib=%s\n' "$THREADS" "$MAX_RAM_GIB"
        if [[ "$ROUTE" == "prebuilt" ]]; then
            printf 'prebuilt_url=%s\n' "$PREBUILT_URL"
            printf 'prebuilt_url_verified_by_operator=UNKNOWN\n'
        else
            printf 'fasta_list=%s\n' "$FASTA_LIST"
            printf 'accession2taxid=%s\n' "$ACCESSION2TAXID"
            printf 'taxdump_source=%s\n' "$TAXDUMP_SRC"
        fi
        # The taxonomy a screen was built against decides whether its taxids can be
        # joined to the watchlist at all. Record it where a reader will trip over it.
        printf 'taxonomy_namespace=%s\n' "$(detect_taxonomy_namespace "$DB_DIR")"
        printf 'host=%s\n' "$(hostname)"
    } > "$out"
    log "wrote build manifest ${out}"
}

# ---------------------------------------------------------------------------
# Plan / main
# ---------------------------------------------------------------------------
# Written as an if, not `(( ... )) && ...`: under `set -e` a trailing false
# arithmetic test at the top level would abort the script.
if (( MIN_FREE_GIB == 0 )); then
    MIN_FREE_GIB="$(default_min_free)"
fi

log "==================== bioThreatFinder metabuli build ===================="
log "database dir : ${DB_DIR}"
log "mode         : ${MODE}"
log "route        : ${ROUTE}$( [[ "$ROUTE" == "prebuilt" ]] && echo "  (RECOMMENDED)" || echo "  (fallback - a local build is days of work)" )"
if [[ "$ROUTE" == "prebuilt" ]]; then
    log "prebuilt url : ${PREBUILT_URL}"
    log "               ^ VERIFY this against the Metabuli docs before trusting it"
else
    log "fasta list   : ${FASTA_LIST}"
    log "acc2taxid    : ${ACCESSION2TAXID}"
    log "taxonomy     : copy nodes/names/merged/delnodes.dmp from ${TAXDUMP_SRC}"
fi
log "threads      : ${THREADS}"
log "max-ram      : ${MAX_RAM_GIB} GiB (metabuli's own default is 128)"
log "taxids       : $( [[ "$MODE" == "gtdb" ]] && echo "GTDB - NOT NCBI; needs translation before the watchlist can be joined" || echo "NCBI - joins packs/*/watchlist.resolved.tsv directly" )"
log "placement    : NVMe (${BTF_DB}) on purpose - query speed is index-I/O bound"
log "untouched    : /media/36tb/db/** (shared, read-only for this project)"
log "EXPECT MANY HOURS (a local GTDB-scale build is DAYS). Run under nohup, tmux or screen."
log "======================================================================="

if (( DRY_RUN == 1 )); then
    log "--dry-run: printing the plan, downloading and building nothing"
fi

check_free_space
run mkdir -p "$DB_DIR" "$STAMP_DIR"

# ONE BUILD PER DATABASE DIRECTORY (see build_kraken2_db.sh for the incident that
# prompted this). Stamps make resuming safe SEQUENTIALLY and say nothing about
# concurrency; two builders sharing a directory can interleave partial downloads into
# a reference file, and a corrupt reference in an index does not error - it
# misclassifies. The lock releases on exit however the shell exits.
if (( DRY_RUN == 0 )); then
    exec {BTF_LOCK_FD}> "${DB_DIR}/.btf_build.lock"
    if ! flock -n "$BTF_LOCK_FD"; then
        holder="$(cat "${DB_DIR}/.btf_build.owner" 2>/dev/null || echo unknown)"
        die "another build already holds the lock on ${DB_DIR} (${holder}).
  Wait for it, or stop it before resuming. Completed steps are stamped, so resuming
  afterwards restarts nothing."
    fi
    printf 'pid=%s started=%s host=%s\n' "$$" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        "$(hostname)" > "${DB_DIR}/.btf_build.owner"
fi

if [[ "$ROUTE" == "prebuilt" ]]; then
    fetch_prebuilt
else
    prepare_taxonomy
    build_local
fi
write_manifest

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo >&2
log "=============================== SUMMARY ==============================="
if (( DRY_RUN == 1 )); then
    log "dry run complete - nothing was downloaded or built"
else
    log "  total: $(du -sh "$DB_DIR" 2>/dev/null | cut -f1) in ${DB_DIR}"
    log "free space remaining: $(df -h --output=avail "$BTF_DB" | tail -1 | tr -d ' ')"
    log "Point screen.metabuli_db at: ${DB_DIR}"
    log "  (or at the nested directory the archive created - check the listing above)"
    log "Query with: ${METABULI} classify --threads ${THREADS} --max-ram ${MAX_RAM_GIB} ..."
    log "  VERIFY the classify subcommand's flags for your version."
    if [[ "$MODE" == "gtdb" ]]; then
        warn "This index uses GTDB taxids. Do NOT join its report to"
        warn "  packs/*/watchlist.resolved.tsv without a GTDB->NCBI translation:"
        warn "  the taxids will not match and the screen will look like a clean negative."
    fi
fi
log "======================================================================="
