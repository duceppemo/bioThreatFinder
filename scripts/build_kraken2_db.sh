#!/usr/bin/env bash
#
# build_kraken2_db.sh - build a fresh Kraken2 "PlusPF"-style database for the
#                       bioThreatFinder stage-1 broad screen, plus Bracken
#                       k-mer distributions for short-read abundance.
#
# Target path (must match screen.kraken2_db in packs/*/pack.yaml):
#   ${db_root}/bioThreatFinder/kraken2/pluspf
#
# Why a fresh build at all
# ------------------------
# The Kraken2 databases already on this box are from 2022. For biosurveillance
# that is not merely "a bit old": SARS-CoV-2 lineages, the mpox clade II
# outbreak genomes, the ICTV renaming of the bunyavirus/arenavirus families and
# large swathes of Candida auris and Burkholderia genomes all postdate them. A
# taxon that does not exist in the database cannot be detected, and a taxon that
# was renamed reads out under a taxid the watchlist no longer matches. So this
# script builds into a NEW directory and never touches
# ${db_root}/kraken2/standard_2022-03-15 or minikraken2_v2_8GB_201904_UPDATE.
#
# FASTER ALTERNATIVE - consider this first
# ----------------------------------------
# The Langmead lab publishes prebuilt Kraken2/Bracken indexes (the "k2" AWS
# indexes page, https://benlangmead.github.io/aws-indexes/k2). A PlusPF index is
# roughly 70-80 GB compressed and downloads in under an hour, versus many hours
# of CPU here, and it ships the Bracken .kmer_distrib files already built.
# Build locally only when you need a specific collection date, a custom library
# set, or a non-default k/minimizer. VERIFY: check the current sizes and dates on
# that page; both change with each release.
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
DB_ROOT="/media/36tb/db"
# PlusPF = Kraken2 "standard" (archaea, bacteria, viral, plasmid, human,
# UniVec_Core) + protozoa + fungi. The protozoa and fungi libraries are the
# whole point for this pipeline: the human watchlist carries 20 parasites and 13
# fungi, and standard would miss every one of them at stage 1.
LIBRARIES="archaea bacteria viral plasmid fungi protozoa human UniVec_Core"
# 64 cores are available, but the build's memory-bound phases stop scaling well
# before that and this is a shared machine.
THREADS=32
KMER=35          # kraken2 default; also what Bracken's default models assume
MINIMIZER=31     # kraken2 default minimizer length (must be <= KMER)
MAX_DB_SIZE=""   # empty = no cap (see the comment at build_db)
SKIP_BRACKEN=0
BRACKEN_READ_LENGTHS=(100 150 250)
DRY_RUN=0
CLEAN=0
REUSE_TAXDUMP=1
# rsync (port 873) is blocked on many networks, including this host, while HTTPS is
# open. kraken2-build defaults to rsync and dies with "Error downloading assembly
# summary file", which reads like an NCBI outage rather than a local firewall.
# 'auto' probes port 873 and falls back to --use-ftp (wget over HTTPS/FTP).
TRANSPORT="auto"
SKIP_MAPS=1
FORCE=0

# Existing, current, unpacked NCBI taxdump on the shared volume.
# Default is the PINNED, dated project snapshot, not the shared one on /media/36tb.
# The index and the watchlist must resolve against the SAME taxonomy snapshot; if they
# differ, taxids disagree and watchlist lookups miss, giving an empty result that is
# indistinguishable from a clean negative. The shared /media/36tb/db/taxdump was
# 6.5 months stale and is used by other work on this machine, so it is left alone.
TAXDUMP_SRC="/media/8tb/bioThreatFinder/taxdump/2026-09-29"

# Pre-flight free-space requirement, GiB.
#   library/    ~300-400 GiB of downloaded FASTA (bacteria dominates)
#   taxonomy/   ~1 GiB with --skip-maps, ~25 GiB with the accession2taxid maps
#   build       ~100 GiB for the final hash.k2d plus transient working files
# 700 GiB is the conservative total for a PlusPF build from scratch. The final
# database is only ~100-120 GiB; the rest is recoverable with --clean afterwards.
MIN_FREE_GIB=700

usage() {
    cat <<'EOF'
Usage:
  build_kraken2_db.sh [options]

Builds a fresh Kraken2 PlusPF-style database at
  <db-root>/bioThreatFinder/kraken2/pluspf
and Bracken k-mer distributions for 100/150/250 bp reads.

Options:
  --db-root DIR        Shared database root (default: /media/36tb/db)
  --threads N          Threads for download and build (default: 32)
  --libraries "A B C"  Space-separated kraken2-build library names
                       (default: "archaea bacteria viral plasmid fungi protozoa human UniVec_Core")
                       Valid names: archaea bacteria plasmid viral human fungi
                       plant protozoa UniVec UniVec_Core nt
  --kmer N             k-mer length (default: 35)
  --minimizer N        minimizer length, must be <= --kmer (default: 31)
  --max-db-size SIZE   Cap the hash table (accepts e.g. 400G, or bytes).
                       Omitted by default - see WHY in the script.
  --read-lengths "100 150"  Bracken read lengths (default: "100 150 250")
  --skip-bracken       Do not build Bracken distributions
  --transport MODE     auto | rsync | ftp (default: auto). kraken2-build uses
                       rsync by default; many networks block port 873 while
                       allowing HTTPS, and the resulting error looks like an NCBI
                       outage. 'auto' probes and falls back to ftp.
  --use-ftp            Shorthand for --transport ftp.
  --taxdump DIR        Pinned NCBI taxonomy snapshot to reuse
                       (default: /media/8tb/bioThreatFinder/taxdump/2026-09-29).
                       MUST be the same snapshot the watchlist was resolved
                       against, or taxids will silently disagree.
  --download-taxonomy  Force kraken2-build --download-taxonomy instead of
                       reusing /media/36tb/db/taxdump
  --with-maps          Also download the accession2taxid maps (~25 GB). Needed
                       only if you will --add-to-library your own sequences
                       whose headers carry no taxid.
  --min-free-gib N     Override the pre-flight free-space requirement
                       (default: 700 GiB transient; see the breakdown it prints)
  --clean              After a successful build, run kraken2-build --clean to
                       delete library/ and taxonomy/ (recovers several hundred
                       GB, but any future rebuild must re-download everything)
  --force              Rebuild steps even when their outputs already exist
  --dry-run            Print the plan and the exact commands, download nothing
  -h, --help           This message

This takes MANY HOURS (commonly 6-16 h on this hardware: the RefSeq library
download is network-bound and the hash build is CPU- and memory-bound). Run it
under nohup/tmux/screen. Every step is skipped if its output already exists, so
an interrupted run can simply be restarted.

Existing 2022 databases under <db-root>/kraken2/ are never touched.
EOF
}

# ---------------------------------------------------------------------------
# Logging
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
        --db-root)           DB_ROOT="${2:?--db-root needs a value}"; shift 2 ;;
        --threads)           THREADS="${2:?--threads needs a value}"; shift 2 ;;
        --libraries)         LIBRARIES="${2:?--libraries needs a value}"; shift 2 ;;
        --kmer)              KMER="${2:?--kmer needs a value}"; shift 2 ;;
        --minimizer)         MINIMIZER="${2:?--minimizer needs a value}"; shift 2 ;;
        --max-db-size)       MAX_DB_SIZE="${2:?--max-db-size needs a value}"; shift 2 ;;
        --read-lengths)      read -r -a BRACKEN_READ_LENGTHS <<< "${2:?--read-lengths needs a value}"; shift 2 ;;
        --skip-bracken)      SKIP_BRACKEN=1; shift ;;
        --taxdump)           TAXDUMP_SRC="${2:?--taxdump needs a value}"; shift 2 ;;
        --transport)         TRANSPORT="${2:?--transport needs auto|rsync|ftp}"; shift 2 ;;
        --use-ftp)           TRANSPORT="ftp"; shift ;;
        --download-taxonomy) REUSE_TAXDUMP=0; shift ;;
        --with-maps)         SKIP_MAPS=0; shift ;;
        --clean)             CLEAN=1; shift ;;
        --force)             FORCE=1; shift ;;
        --dry-run)           DRY_RUN=1; shift ;;
        --min-free-gib)      MIN_FREE_GIB="${2:?--min-free-gib needs a value}"; shift 2 ;;
        -h|--help)           usage; exit 0 ;;
        *)                   usage; die "unknown option: $1" ;;
    esac
done

[[ "$THREADS"   =~ ^[0-9]+$ ]] || die "--threads must be an integer"
[[ "$KMER"      =~ ^[0-9]+$ ]] || die "--kmer must be an integer"
[[ "$MINIMIZER" =~ ^[0-9]+$ ]] || die "--minimizer must be an integer"
(( MINIMIZER <= KMER )) || die "--minimizer ($MINIMIZER) must be <= --kmer ($KMER)"

DB_DIR="${DB_ROOT}/bioThreatFinder/kraken2/pluspf"
STAMP_DIR="${DB_DIR}/.build_stamps"

# Hard guard: never let a mistyped --db-root make us build on top of the stale
# 2022 databases. They are the only thing currently able to classify anything on
# this box and are referenced by other projects.
case "$DB_DIR" in
    */bioThreatFinder/kraken2/pluspf) : ;;
    *) die "refusing to build at '${DB_DIR}' - expected a path ending in /bioThreatFinder/kraken2/pluspf" ;;
esac
for protected in "${DB_ROOT}/kraken2/standard_2022-03-15" \
                 "${DB_ROOT}/kraken2/minikraken2_v2_8GB_201904_UPDATE" \
                 "${DB_ROOT}/kraken2/taxonomy"; do
    if [[ "$DB_DIR" == "$protected" || "$DB_DIR" == "$protected"/* ]]; then
        die "refusing to write inside the pre-existing database ${protected}"
    fi
done

# ---------------------------------------------------------------------------
# Tool resolution
# ---------------------------------------------------------------------------
CONDA_BIN="/home/bioinfo/miniconda3/bin"
resolve_tool() {
    local name="$1" p
    if p="$(command -v "$name" 2>/dev/null)"; then printf '%s\n' "$p"; return 0; fi
    if [[ -x "${CONDA_BIN}/${name}" ]]; then printf '%s\n' "${CONDA_BIN}/${name}"; return 0; fi
    return 1
}

KRAKEN2_BUILD="$(resolve_tool kraken2-build || true)"
if [[ -z "$KRAKEN2_BUILD" ]]; then
    msg="kraken2-build not found.
       Install into a dedicated environment, e.g.:
         ${CONDA_BIN}/mamba create -y -n kraken2 -c bioconda -c conda-forge kraken2 bracken
         conda activate kraken2
       kraken2-build also needs rsync, wget/curl and dustmasker (blast+) on PATH."
    # Under --dry-run the point is to show the plan, so carry on with the bare
    # command name instead of aborting.
    if (( DRY_RUN == 1 )); then warn "$msg"; KRAKEN2_BUILD="kraken2-build"; else die "$msg"; fi
fi
BRACKEN_BUILD="$(resolve_tool bracken-build || true)"
if [[ -z "$BRACKEN_BUILD" && $SKIP_BRACKEN -eq 0 ]]; then
    msg="bracken-build not found (needed unless --skip-bracken).
       ${CONDA_BIN}/mamba install -c bioconda bracken"
    if (( DRY_RUN == 1 )); then warn "$msg"; BRACKEN_BUILD="bracken-build"; else die "$msg"; fi
fi
# kraken2-build shells out to a downloader; failing here beats failing four hours in.
#
# It requires EITHER rsync OR wget, depending on the transport chosen below - and this
# check used to hard-fail on a missing rsync unconditionally. That is wrong on any host
# where rsync's port 873 is blocked (this one), because the build then legitimately runs
# over wget and never touches rsync. Warn about whichever is absent and let
# choose_transport() decide; ensure_working_wget() enforces the wget side.
if (( DRY_RUN == 0 )); then
    _have_rsync=0; _have_wget=0
    command -v rsync >/dev/null 2>&1 && _have_rsync=1
    command -v wget  >/dev/null 2>&1 && _have_wget=1
    if (( _have_rsync == 0 && _have_wget == 0 )); then
        die "kraken2-build needs rsync or wget to download libraries and neither is on
  PATH. Install one:  mamba install -n <env> -c conda-forge wget"
    fi
    (( _have_rsync == 0 )) && warn "rsync not found - only the ftp/wget transport is available"
    (( _have_wget  == 0 )) && warn "wget not found - only the rsync transport is available"
    command -v dustmasker >/dev/null 2>&1 \
        || warn "dustmasker (from blast+) not found; it is only needed if you pass --no-masking's opposite, i.e. default low-complexity masking. Install blast+ if the build complains."
fi

# ---------------------------------------------------------------------------
# Pre-flight
# ---------------------------------------------------------------------------
check_free_space() {
    local probe="$DB_ROOT" avail_gib
    while [[ ! -d "$probe" && "$probe" != "/" ]]; do probe="$(dirname "$probe")"; done
    avail_gib="$(df -BG --output=avail "$probe" | tail -1 | tr -dc '0-9')"
    [[ -n "$avail_gib" ]] || die "could not determine free space for ${probe}"
    log "free space under ${probe}: ${avail_gib} GiB"
    log "required (transient) : ${MIN_FREE_GIB} GiB"
    log "  breakdown: ~300-400 GiB downloaded RefSeq FASTA in library/,"
    log "             ~100 GiB for hash.k2d + transient build files,"
    log "             ~1 GiB taxonomy (~25 GiB with --with-maps)."
    log "             The finished DB is ~100-120 GiB; --clean recovers the rest."
    if (( avail_gib < MIN_FREE_GIB )); then
        die "insufficient free space: ${avail_gib} GiB available, ${MIN_FREE_GIB} GiB required.
       This volume is already ~97% full. Options:
         * free space elsewhere on the volume, or
         * download a prebuilt PlusPF index instead (see the header comment:
           ~70-80 GB, no transient library/ at all), or
         * lower the requirement with --min-free-gib if you are certain.
       Nothing was downloaded."
    fi
}

stamp()      { printf '%s\n' "${STAMP_DIR}/$1.done"; }
have_stamp() { [[ -f "$(stamp "$1")" ]]; }
set_stamp()  { (( DRY_RUN == 1 )) || { mkdir -p "$STAMP_DIR"; date -u +%Y-%m-%dT%H:%M:%SZ > "$(stamp "$1")"; }; }

# ---------------------------------------------------------------------------
# Taxonomy
# ---------------------------------------------------------------------------
# Can kraken2-build consume the existing /media/36tb/db/taxdump?
#
# Partly, and it is worth doing. There is NO kraken2-build flag that points at an
# external taxonomy directory: `--download-taxonomy` always re-fetches. But
# `--build` only reads ${DB}/taxonomy/nodes.dmp and ${DB}/taxonomy/names.dmp, so
# pre-populating ${DB}/taxonomy with those two files from the existing unpacked
# taxdump and then NOT calling --download-taxonomy works, and saves both the
# taxdump download and (much more importantly) the ~25 GB accession2taxid maps.
#
# The catch: the accession2taxid maps are what let kraken2-build assign a taxid
# to a sequence whose header carries none. Sequences fetched by
# `--download-library` do carry one - download_genomic_library.sh rewrites their
# headers to `>kraken:taxid|<taxid>|...` using assembly_summary.txt, and writes a
# prelim_map.txt beside them - so for a pure RefSeq library build the maps are
# not needed. They ARE needed if you later `--add-to-library` your own FASTA
# without a kraken:taxid header. Hence --with-maps.
#
# VERIFY: this depends on kraken2-build internals (checked against kraken2
# 2.1.x). If --build complains about a missing taxonomy file, rerun with
# --download-taxonomy and let kraken2 fetch its own copy.
prepare_taxonomy() {
    local taxdir="${DB_DIR}/taxonomy"
    if have_stamp taxonomy && (( FORCE == 0 )); then
        skip "taxonomy already prepared (${taxdir})"
        return 0
    fi

    if (( REUSE_TAXDUMP == 1 )); then
        if [[ -s "${TAXDUMP_SRC}/nodes.dmp" && -s "${TAXDUMP_SRC}/names.dmp" ]]; then
            # Refuse to bake a stale taxonomy into a 100+ GB index. kraken2-build
            # downloads CURRENT library sequences whose kraken:taxid| headers carry
            # today's taxids; pairing those with an old taxonomy drops or misassigns
            # every organism added since the snapshot, silently.
            _td_age_days=$(( ( $(date -u +%s) - $(stat -c %Y "${TAXDUMP_SRC}/nodes.dmp") ) / 86400 ))
            log "taxdump ${TAXDUMP_SRC} is ${_td_age_days} day(s) old"
            if (( _td_age_days > 90 )); then
                die "taxdump is ${_td_age_days} days old. kraken2-build downloads CURRENT
  library sequences, so an old taxonomy drops or misassigns every organism added
  since the snapshot - and the index would then disagree with the watchlist.
  Refresh from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/ and pass
  --taxdump, or pass --download-taxonomy to let kraken2-build fetch its own
  (which then will NOT match the watchlist)."
            fi
            log "reusing existing taxdump from ${TAXDUMP_SRC}"
            run mkdir -p "$taxdir"
            # Copy rather than symlink: kraken2-build has historically rewritten
            # files in ${DB}/taxonomy, and a symlink would then corrupt the
            # shared taxdump that resolve_watchlist.py and other projects pin.
            run cp -f "${TAXDUMP_SRC}/nodes.dmp" "${taxdir}/nodes.dmp"
            run cp -f "${TAXDUMP_SRC}/names.dmp" "${taxdir}/names.dmp"
            # merged.dmp/delnodes.dmp are not required by --build, but having them
            # beside the DB records which taxonomy snapshot it was built from.
            for extra in merged.dmp delnodes.dmp; do
                if [[ -s "${TAXDUMP_SRC}/${extra}" ]]; then
                    run cp -f "${TAXDUMP_SRC}/${extra}" "${taxdir}/${extra}"
                fi
            done
            if (( SKIP_MAPS == 0 )); then
                log "--with-maps: fetching accession2taxid maps only (~25 GB)"
                run "$KRAKEN2_BUILD" --db "$DB_DIR" --download-taxonomy --threads "$THREADS"
                # --download-taxonomy FETCHES AND UNPACKS NCBI'S CURRENT TAXDUMP INTO
                # THIS SAME DIRECTORY, replacing the pinned nodes.dmp/names.dmp copied
                # above - and then the step was stamped complete, so re-running skipped
                # it. The result was a ~100 GB index built against current taxonomy
                # while the watchlist is pinned to a dated snapshot: the two disagree,
                # lookups miss for every drifted taxon, and nothing says so. Restore
                # the pinned files and VERIFY, rather than assuming the copy held.
                log "restoring pinned taxonomy over kraken2-build's downloaded copy"
                run cp -f "${TAXDUMP_SRC}/nodes.dmp" "${taxdir}/nodes.dmp"
                run cp -f "${TAXDUMP_SRC}/names.dmp" "${taxdir}/names.dmp"
                if (( DRY_RUN == 0 )); then
                    for f in nodes.dmp names.dmp; do
                        a="$(md5sum < "${TAXDUMP_SRC}/${f}" | cut -d" " -f1)"
                        b="$(md5sum < "${taxdir}/${f}" | cut -d" " -f1)"
                        [[ "$a" == "$b" ]] || die "after --with-maps, ${taxdir}/${f} does
  not match the pinned ${TAXDUMP_SRC}/${f} (${b} vs ${a}). The index would be built
  against a different taxonomy than the watchlist resolves against, which produces
  silent lookup misses rather than an error. Refusing to stamp this step."
                    done
                    log "pinned taxonomy verified against ${TAXDUMP_SRC}"
                fi
            fi
            set_stamp taxonomy
            return 0
        fi
        warn "${TAXDUMP_SRC} is missing nodes.dmp/names.dmp; falling back to --download-taxonomy"
    fi

    local extra_args=()
    if (( SKIP_MAPS == 1 )); then extra_args+=(--skip-maps); fi
    run "$KRAKEN2_BUILD" --db "$DB_DIR" --download-taxonomy --threads "$THREADS" ${extra_args[@]+"${extra_args[@]}"}
    set_stamp taxonomy
}

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------
# One kraken2-build call per library, each with its own stamp file, so an
# interrupted run resumes at the library that failed instead of re-downloading
# the ~200 GB bacteria set.
# MINIMUM VERSION GUARD.
#
# The stage-1 design rests on Kraken2's --report-minimizer-data columns (distinct
# minimizer count and duplicity), which is what lets Kraken2 substitute for a separate
# KrakenUniq index. That flag arrived in 2.1.0. Building with anything older produces an
# index that works but cannot supply the evidence the pipeline scores on - a silent
# capability loss, not an error.
#
# This is not hypothetical: an unconstrained `mamba install kraken2` into an environment
# that also held mmseqs2, hmmer and mafft resolved to 2.0.8-beta (2019), which has no
# such flag. Installing kraken2 in its OWN environment is what avoids that.
require_kraken2_version() {
    local v
    v="$("$KRAKEN2_BUILD" --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
    if [[ -z "$v" ]]; then
        warn "could not determine the kraken2 version; proceeding, but verify that"
        warn "  --report-minimizer-data is supported or stage 1 loses its k-mer evidence."
        return
    fi
    local maj min
    maj="${v%%.*}"; min="$(printf '%s' "$v" | cut -d. -f2)"
    log "kraken2 version: ${v}"
    if (( maj < 2 )) || { (( maj == 2 )) && (( min < 1 )); }; then
        die "kraken2 ${v} is too old. --report-minimizer-data (added in 2.1.0) is what
  supplies the distinct-minimizer and duplicity columns this pipeline scores stage 1
  on; an index built here would work but silently carry no k-mer evidence.
  Install kraken2 in its OWN environment so the solver is not constrained:
    mamba create -n btf-kraken2 -c conda-forge -c bioconda kraken2=2.17.2 bracken"
    fi
    if ! kraken2 --help 2>&1 | grep -q 'report-minimizer-data'; then
        die "this kraken2 (${v}) does not advertise --report-minimizer-data. Stage 1
  would lose the distinct-k-mer evidence it scores on. Install a newer kraken2."
    fi
}

# kraken2's --use-ftp path shells out to `wget`, and a broken wget produces the same
# "Error downloading assembly summary file" as a network failure - which sent this
# build chasing NCBI and firewalls for a while. The conda-forge wget in this env is
# linked against libunistring.so.2 while the env ships .so.5, so it cannot start at
# all. Detect that and fall back to a working wget rather than mis-diagnosing it.
#
# A shim directory is used instead of removing the conda package, because kraken2
# DECLARES wget as a dependency - removing it uninstalls kraken2 and bracken too.
ensure_working_wget() {
    if wget --version >/dev/null 2>&1; then
        log "wget: $(command -v wget) (working)"
        return
    fi
    local broken; broken="$(command -v wget 2>/dev/null || echo none)"
    local cand
    for cand in /usr/bin/wget /bin/wget /usr/local/bin/wget; do
        if [[ -x "$cand" ]] && "$cand" --version >/dev/null 2>&1; then
            local shim="${BTF_SHIM_DIR:-${DB_ROOT}/bioThreatFinder/toolshim}"
            mkdir -p "$shim"
            ln -sf "$cand" "${shim}/wget"
            export PATH="${shim}:${PATH}"
            warn "wget at ${broken} is broken (cannot load its shared libraries)."
            warn "  Shimming ${cand} ahead of it via ${shim}."
            warn "  Left the conda package in place on purpose: kraken2 declares wget as"
            warn "  a dependency, so removing it uninstalls kraken2 and bracken as well."
            return
        fi
    done
    die "wget at ${broken} cannot run and no working wget was found.
  kraken2-build needs it for --use-ftp, and rsync (port 873) is unavailable here.
  Install a working wget, or pass --transport rsync if port 873 is open for you."
}

# Decide the transport once, before spending hours discovering it the hard way.
choose_transport() {
    DL_OPTS=""
    case "$TRANSPORT" in
        rsync) log "transport: rsync (forced)" ;;
        ftp)   DL_OPTS="--use-ftp"; log "transport: ftp/https (forced)" ;;
        auto)
            if timeout 12 bash -c 'cat < /dev/null > /dev/tcp/ftp.ncbi.nlm.nih.gov/873' 2>/dev/null; then
                log "transport: rsync (port 873 reachable)"
            else
                DL_OPTS="--use-ftp"
                log "transport: ftp/https - rsync port 873 is NOT reachable from this host."
                log "  Without this fallback kraken2-build fails with 'Error downloading"
                log "  assembly summary file', which looks like an NCBI problem but is a"
                log "  local firewall blocking the rsync port."
            fi
            ;;
        *) die "--transport must be auto, rsync or ftp (got '${TRANSPORT}')" ;;
    esac
}

download_libraries() {
    local lib
    require_kraken2_version
    choose_transport
    [[ -n "$DL_OPTS" ]] && ensure_working_wget
    for lib in $LIBRARIES; do
        if have_stamp "library_${lib}" && (( FORCE == 0 )); then
            skip "library already downloaded: ${lib}"
            continue
        fi
        log "downloading library '${lib}' (this is the slow part; bacteria alone is hours)"
        run "$KRAKEN2_BUILD" --db "$DB_DIR" --download-library "$lib" --threads "$THREADS" ${DL_OPTS}
        set_stamp "library_${lib}"
    done
}

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
# WHY --max-db-size is unset by default: it does not merely truncate, it makes
# kraken2-build subsample minimizers to fit, which lowers sensitivity in a way
# that is invisible at run time - you simply stop seeing some low-abundance
# hits. For a biosurveillance screen that is the wrong trade. The box has 503 GB
# of RAM and the full PlusPF hash is ~100 GB, so there is no need to shrink it.
# Pass --max-db-size only when the DB must run on a smaller machine, and record
# that fact in the run metadata.
build_db() {
    if [[ -s "${DB_DIR}/hash.k2d" && -s "${DB_DIR}/taxo.k2d" && -s "${DB_DIR}/opts.k2d" ]] \
       && (( FORCE == 0 )); then
        skip "database already built ($(du -sh "$DB_DIR" 2>/dev/null | cut -f1) in ${DB_DIR})"
        return 0
    fi
    local args=(--db "$DB_DIR" --build --threads "$THREADS"
                --kmer-len "$KMER" --minimizer-len "$MINIMIZER")
    if [[ -n "$MAX_DB_SIZE" ]]; then args+=(--max-db-size "$MAX_DB_SIZE"); fi
    warn "the hash build is the CPU/RAM-heavy phase and can run for several hours"
    run "$KRAKEN2_BUILD" "${args[@]}"
    if (( DRY_RUN == 0 )); then
        [[ -s "${DB_DIR}/hash.k2d" ]] || die "build finished but ${DB_DIR}/hash.k2d is missing or empty"
    fi
}

# ---------------------------------------------------------------------------
# Bracken
# ---------------------------------------------------------------------------
# Bracken redistributes reads assigned to internal nodes down to species level,
# using a k-mer distribution that is specific to BOTH the database and the read
# length. Hence one distribution per read length: 100, 150 and 250 bp cover the
# Illumina/Element configurations in use.
#
# ONT: do NOT use Bracken for nanopore data. Its model assumes reads of a single
# fixed length; ONT read lengths span orders of magnitude within one run, so the
# redistribution is not just approximate but biased. For ONT, bioThreatFinder
# uses the Kraken2 report (and the per-read assignments) directly, and relies on
# stage-2 competitive alignment for quantitative claims.
build_bracken() {
    if (( SKIP_BRACKEN == 1 )); then
        log "--skip-bracken: not building k-mer distributions"
        return 0
    fi
    local rl
    for rl in "${BRACKEN_READ_LENGTHS[@]}"; do
        local out="${DB_DIR}/database${rl}mers.kmer_distrib"
        if [[ -s "$out" ]] && (( FORCE == 0 )); then
            skip "bracken distribution exists: $(basename "$out")"
            continue
        fi
        log "building bracken distribution for read length ${rl}"
        # -k must match the k used for the kraken2 DB, or the distribution is
        # silently wrong.
        run "$BRACKEN_BUILD" -d "$DB_DIR" -t "$THREADS" -k "$KMER" -l "$rl"
    done
}

# ---------------------------------------------------------------------------
# Plan / main
# ---------------------------------------------------------------------------
log "===================== bioThreatFinder kraken2 build ====================="
log "database dir : ${DB_DIR}"
log "libraries    : ${LIBRARIES}"
log "k / minimizer: ${KMER} / ${MINIMIZER}"
log "threads      : ${THREADS}"
log "taxonomy     : $( (( REUSE_TAXDUMP == 1 )) && echo "reuse ${TAXDUMP_SRC}" || echo "kraken2-build --download-taxonomy" )$( (( SKIP_MAPS == 1 )) && echo " (no accession2taxid maps)" )"
log "bracken      : $( (( SKIP_BRACKEN == 1 )) && echo "skipped" || echo "read lengths ${BRACKEN_READ_LENGTHS[*]}" )"
log "max-db-size  : ${MAX_DB_SIZE:-<unset: full-sensitivity hash>}"
log "untouched    : ${DB_ROOT}/kraken2/* (the stale 2022 databases)"
log "EXPECT THIS TO TAKE MANY HOURS. Run under nohup, tmux or screen."
log "========================================================================="

if (( DRY_RUN == 1 )); then
    log "--dry-run: printing the plan, downloading nothing"
fi

check_free_space
run mkdir -p "$DB_DIR" "$STAMP_DIR"

# ---------------------------------------------------------------------------
# ONE BUILD PER DATABASE DIRECTORY, ENFORCED.
#
# This was learned the hard way during the PlusPF build. A resume was started while
# the original build was still running, both invoked `kraken2-build
# --download-library bacteria` against the same directory, and the second one's
# internal cleanup failed with `rm: cannot remove 'all/': Directory not empty`. Two
# concurrent FTP transfers into one staging directory can interleave partial files
# into library.fna, and a truncated or duplicated reference in a classification index
# is invisible: it does not error, it just misclassifies or silently fails to match.
#
# The stamps make resuming safe SEQUENTIALLY; they say nothing about concurrency.
# flock closes that gap. It is released automatically when this shell exits, however
# it exits, so a killed build does not leave the directory permanently locked.
# ---------------------------------------------------------------------------
if (( DRY_RUN == 0 )); then
    exec {BTF_LOCK_FD}> "${DB_DIR}/.btf_build.lock"
    if ! flock -n "$BTF_LOCK_FD"; then
        holder="$(cat "${DB_DIR}/.btf_build.owner" 2>/dev/null || echo unknown)"
        die "another build already holds the lock on ${DB_DIR} (${holder}).
  Two kraken2-build processes sharing one database directory can interleave partial
  downloads into library.fna, and a corrupt reference in the index does not announce
  itself - it misclassifies. Wait for that build, or stop it before resuming.
  Completed steps are stamped, so resuming afterwards restarts nothing."
    fi
    printf 'pid=%s started=%s host=%s\n' "$$" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        "$(hostname)" > "${DB_DIR}/.btf_build.owner"
fi

prepare_taxonomy
download_libraries
build_db
build_bracken

if (( CLEAN == 1 )); then
    # kraken2-build --clean removes library/ and taxonomy/, keeping only the
    # *.k2d files the classifier actually loads.
    warn "--clean: removing library/ and taxonomy/ from ${DB_DIR}; a future rebuild will need a full re-download"
    run "$KRAKEN2_BUILD" --db "$DB_DIR" --clean
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo >&2
log "=============================== SUMMARY ==============================="
if (( DRY_RUN == 1 )); then
    log "dry run complete - nothing was downloaded or built"
else
    for f in hash.k2d opts.k2d taxo.k2d seqid2taxid.map; do
        if [[ -e "${DB_DIR}/${f}" ]]; then
            log "  $(du -h "${DB_DIR}/${f}" | cut -f1)	${f}"
        else
            warn "  MISSING	${f}"
        fi
    done
    if (( SKIP_BRACKEN == 0 )); then
        for rl in "${BRACKEN_READ_LENGTHS[@]}"; do
            if [[ -s "${DB_DIR}/database${rl}mers.kmer_distrib" ]]; then
                log "  bracken	database${rl}mers.kmer_distrib"
            else
                warn "  bracken	database${rl}mers.kmer_distrib MISSING"
            fi
        done
    fi
    log "  total: $(du -sh "$DB_DIR" 2>/dev/null | cut -f1) in ${DB_DIR}"
    if [[ -d "${DB_DIR}/library" && $CLEAN -eq 0 ]]; then
        log "  library/ still holds $(du -sh "${DB_DIR}/library" 2>/dev/null | cut -f1)."
        log "  Recover it with: ${KRAKEN2_BUILD:-kraken2-build} --db ${DB_DIR} --clean"
    fi
    log "free space remaining: $(df -h --output=avail "$DB_ROOT" | tail -1 | tr -d ' ')"
    log "Point screen.kraken2_db at: ${DB_DIR}"
fi
log "======================================================================="
