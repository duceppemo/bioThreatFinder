#!/usr/bin/env bash
#
# build_krakenuniq_db.sh - build the KrakenUniq index for bioThreatFinder's
#                          OPT-IN second-opinion stage-1 screen.
#
# READ docs/decisions.md D2 BEFORE RUNNING THIS. KrakenUniq is off by default
# (--skip_krakenuniq defaults true) and Kraken2 is the primary nucleotide arm. It
# is only switched on if the bake-off (scripts/run_bakeoff.sh, variant
# plus_krakenuniq) shows a confirm-PPV gain worth the cost - and even then D2 says
# use a PREBUILT index rather than this build (~1400 GiB transient disk, a
# Jellyfish 1.x dependency, 12+ hours).
#
# Target path (must match screen.krakenuniq_db in packs/*/pack.yaml):
#   ${BTF_DB}/krakenuniq/standard          (default BTF_DB=/media/8tb/bioThreatFinder)
#
# Why KrakenUniq and not just Kraken2
# -----------------------------------
# Kraken2 reports READ counts per taxon. It never asks how many DISTINCT k-mers
# of that taxon were hit, so 500 reads piled on one conserved 16S/insertion-
# sequence k-mer and 500 reads spread evenly over a 5 Mb chromosome produce the
# same line in the report. That is the documented mechanism behind Kraken2's
# species-level false positives, and it is exactly the failure mode that matters
# here: a watchlist name appearing for the wrong reason.
#
# KrakenUniq answers the missing question directly. For every taxon it maintains
# a HyperLogLog sketch of the distinct database k-mers observed, and emits a
# `kmers` (cardinality) and `dup` column beside the read count. Reads stacked on
# one conserved locus give a large read count with a tiny distinct-k-mer count -
# visible, and filterable, at stage 1 instead of at stage 2.
#
# That argument was the original case for making KrakenUniq primary, and it was
# wrong for Kraken2 >= 2.1.0: --report-minimizer-data emits distinct-minimizer
# counts and duplicity from the same HyperLogLog code (see docs/decisions.md D2).
# Kraken2 is therefore the primary screen; KrakenUniq remains an opt-in second
# opinion, counting distinct k-mers rather than minimizers.
#
# CAN THE EXISTING KRAKEN2 DATABASE BE REUSED? NO.
# ------------------------------------------------
# KrakenUniq descends from Kraken 1, not Kraken 2, and uses Kraken 1's on-disk
# format: database.kdb + database.idx + taxDB. Kraken2's hash.k2d/opts.k2d/
# taxo.k2d is a different, incompatible structure (compact hashed minimizers with
# a spaced seed, not a sorted k-mer list). KrakenUniq reads a Kraken 1 database
# natively; it cannot read a Kraken 2 database at all, and there is no converter
# - the minimizer/spaced-seed transform in Kraken2 is lossy with respect to the
# k-mers KrakenUniq needs to count.
#
# Consequence for disk, and why the pre-flight below quotes a COMBINED figure:
# the two indexes are two separate, additive costs on the same NVMe volume. You
# do not get the KrakenUniq screen "for free" once Kraken2 PlusPF exists; budget
# for both, plus the KrakenUniq build's transient k-mer files, which are the
# largest single item of the three.
#
# The downloaded FASTA library IS shareable in principle (both tools consume
# plain RefSeq FASTA), but kraken2-build rewrites headers to `>kraken:taxid|...`
# in its own library/ layout while krakenuniq-download writes a different layout
# plus its own *.map files, so this script keeps its own library/ rather than
# pointing --library-dir at Kraken2's. Cheaper in engineering risk than in disk.
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
# Project-owned databases live on the FAST NVMe volume. KrakenUniq classification
# is dominated by random reads into database.kdb, so the difference between NVMe
# and spinning disk here is not a nicety.
BTF_DB="/media/8tb/bioThreatFinder"

# Library specs are passed straight to krakenuniq-download, whose grammar is
#   refseq/DOMAIN[/ASS_LEVEL[/COLUMN=value]*]  |  genbank/...  |  contaminants
#   |  viral-neighbors  |  nt  |  microbial-nt
# VERIFY: checked against krakenuniq-download in krakenuniq master (0.7.x).
# Valid DOMAIN values there are exactly: bacteria viral archaea fungi protozoa
# invertebrate plant vertebrate_mammalian vertebrate_other. Note that this is a
# DIFFERENT vocabulary from kraken2-build --download-library (no "human", no
# "plasmid", no "UniVec_Core"), which is why this list looks nothing like the one
# in build_kraken2_db.sh.
#
# WHY these libraries:
#   archaea/bacteria at Complete_Genome  - draft assemblies are the main source
#       of cross-contamination in reference sets, and a contaminated contig
#       becomes a permanent false positive for a watchlist taxon.
#   viral/fungi/protozoa at Any          - the human watchlist carries 20
#       parasites and 13 fungi; restricting these to Complete_Genome would drop
#       most of them, and a taxon absent from the index cannot be detected.
#   vertebrate_mammalian species_taxid=9606 - human, so host reads are absorbed
#       by the human node instead of being mis-assigned to microbial taxa.
#   contaminants                         - UniVec/EmVec vector and adapter
#       sequence; the Kraken2 equivalent of UniVec_Core.
#
# Deliberately NOT in the default list: "viral-neighbors" (NCBI Viral Genome
# Resource strain sequences). It is the single best addition for divergent viral
# strains, but it is fetched one Entrez query at a time and can run for days or
# stall on rate limiting. Metabuli (build_metabuli_db.sh) is this pipeline's
# answer to viral divergence instead. Add it explicitly when you want it:
#   --libraries "... viral-neighbors"
# Its per-library stamp file means adding it later does not re-pull anything else.
LIBRARIES="refseq/archaea/Complete_Genome refseq/bacteria/Complete_Genome refseq/viral/Any refseq/fungi/Any refseq/protozoa/Any refseq/vertebrate_mammalian/Any/species_taxid=9606 contaminants"

# 64 cores exist, but the jellyfish count and the LCA-setting phase are
# memory-bandwidth bound long before 64 threads, and this is a shared machine.
THREADS=32

# KrakenUniq (Kraken 1 lineage) caps k at 31 and requires minimizer < k strictly.
KMER=31
MINIMIZER=15

# HyperLogLog precision. See the long comment at write_classify_defaults().
HLL_PRECISION=12

# Build-time memory policy: "auto" decides between the RAM build and the on-disk
# build from the measured library size. See decide_work_on_disk().
WORK_ON_DISK="auto"     # auto | yes | no
JELLYFISH_BIN=""        # empty = let krakenuniq-build find its own jellyfish 1

# Pre-flight free-space requirement, GiB - see check_free_space() for the
# breakdown and the arithmetic.
MIN_FREE_GIB=1400

WORK_DIR=""             # optional: put the transient k-mer files elsewhere
DRY_RUN=0
FORCE=0
CLEAN=0

# Existing, current, unpacked NCBI taxdump on the shared READ-ONLY volume.
# Pinned dated snapshot, NOT the shared /media/36tb/db/taxdump that D8 measured
# at 201 days stale and decided to stop using.
TAXDUMP_SRC="/media/8tb/bioThreatFinder/taxdump/2026-09-29"
MAX_TAXDUMP_AGE_DAYS=90

CONDA_BIN="/home/bioinfo/miniconda3/bin"

usage() {
    cat <<'EOF'
Usage:
  build_krakenuniq_db.sh [options]

Builds the KrakenUniq index for bioThreatFinder's OPT-IN second-opinion
stage-1 screen (off by default; Kraken2 is primary - see docs/decisions.md D2,
which prefers a prebuilt index over this build), at
  <btf-db>/krakenuniq/standard

KrakenUniq reports DISTINCT k-mer counts per taxon (HyperLogLog cardinality)
alongside read counts. That distinction is the specificity fix Kraken2 cannot
provide: reads piled on one conserved k-mer are separable from reads spread
across a genome.

Options:
  --btf-db DIR         Project-owned database root
                       (default: /media/8tb/bioThreatFinder)
  --taxdump DIR        Unpacked NCBI taxdump to COPY nodes.dmp/names.dmp from
                       (default: /media/36tb/db/taxdump, read-only and shared)
  --threads N          Threads for download and build (default: 32)
  --libraries "A B C"  Space-separated krakenuniq-download library specs.
                       Grammar: refseq/DOMAIN[/ASS_LEVEL[/COL=value]*],
                       genbank/..., contaminants, viral-neighbors, nt,
                       microbial-nt. DOMAIN is one of bacteria viral archaea
                       fungi protozoa invertebrate plant vertebrate_mammalian
                       vertebrate_other. ASS_LEVEL is one of Any
                       Complete_Genome Chromosome Scaffold Contig.
                       NOTE: this vocabulary is NOT kraken2-build's.
  --kmer N             k-mer length, max 31 (default: 31)
  --minimizer N        minimizer length, must be < --kmer (default: 15)
  --hll-precision N    HyperLogLog precision, 10-18 (default: 12). This is a
                       CLASSIFY-time parameter, not a build-time one - see
                       below. The value is recorded in the database directory
                       for the pipeline to use at query time.
  --min-free-gib N     Override the pre-flight free-space requirement
                       (default: 1400 GiB transient; the check prints its
                       breakdown and the combined Kraken2+KrakenUniq figure)
  --work-dir DIR       Hold the transient jellyfish/sorted k-mer files here
                       instead of inside the database directory. Use when the
                       database volume is tight but another fast volume is not.
  --work-on-disk       Force the low-RAM, on-disk build (slower, survives)
  --in-ram             Force the in-RAM build (faster, can be OOM-killed)
                       Default is automatic: see the note below.
  --jellyfish-bin PATH Path to a Jellyfish *version 1* binary
  --clean              After a successful build, remove the transient k-mer
                       files and the downloaded library/ (recovers the bulk of
                       the space; any future rebuild re-downloads everything)
  --force              Redo steps even when their outputs already exist
  --dry-run            Print the plan and the exact commands, download nothing
  -h, --help           This message

--hll-precision is not a krakenuniq-build option
------------------------------------------------
Verified against krakenuniq master: --hll-precision belongs to the `krakenuniq`
CLASSIFIER (range 10-18, default 12), not to `krakenuniq-build`. The sketches
are built while reading a sample, not baked into the index. This script
therefore validates the value and writes it to
  <db>/btf_classify_defaults.env
so the Nextflow module and this build agree, instead of silently accepting a
flag it cannot pass on.

Memory, and why the default is automatic
----------------------------------------
KrakenUniq's build is the RAM-hungry step. By default build_db.sh passes -M
("minimize disk writes") and holds the k-mer database in memory; --work-on-disk
inverts that ("minimize RAM usage") at a significant speed cost. An in-RAM build
of a PlusPF-scale library wants roughly as much RAM as the finished
database.kdb, which can approach or exceed this box's 503 GB. Being OOM-killed
twelve hours in is the worst outcome, so this script estimates the index size
from the downloaded library and chooses; override with --work-on-disk/--in-ram.

Jellyfish 1 is required
-----------------------
The build shells out to Jellyfish and needs version 1.x. The jellyfish on this
box is 2.3.1 and WILL NOT WORK. The bioconda krakenuniq package (>=0.7.3) ships
its own Jellyfish 1; install KrakenUniq into its own environment and let it use
that, or pass --jellyfish-bin.

This takes MANY HOURS - commonly 12-30 h here: the RefSeq download is
network-bound, the jellyfish k-mer count and sort are I/O- and RAM-bound, and
the LCA pass is CPU-bound. Run it under nohup, tmux or screen. Every step
carries a stamp file, so an interrupted run resumes instead of re-downloading
hundreds of GB.

Nothing under /media/36tb/db is ever written, and the stale 2022 Kraken2
databases at /media/36tb/db/kraken2/ are never touched.
EOF
}

# ---------------------------------------------------------------------------
# Logging - same vocabulary as build_kraken2_db.sh so log greps work on both
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
        --btf-db)        BTF_DB="${2:?--btf-db needs a value}"; shift 2 ;;
        --max-taxdump-age-days) MAX_TAXDUMP_AGE_DAYS="${2:?}"; shift 2 ;;
        --taxdump)       TAXDUMP_SRC="${2:?--taxdump needs a value}"; shift 2 ;;
        --threads)       THREADS="${2:?--threads needs a value}"; shift 2 ;;
        --libraries)     LIBRARIES="${2:?--libraries needs a value}"; shift 2 ;;
        --kmer)          KMER="${2:?--kmer needs a value}"; shift 2 ;;
        --minimizer)     MINIMIZER="${2:?--minimizer needs a value}"; shift 2 ;;
        --hll-precision) HLL_PRECISION="${2:?--hll-precision needs a value}"; shift 2 ;;
        --min-free-gib)  MIN_FREE_GIB="${2:?--min-free-gib needs a value}"; shift 2 ;;
        --work-dir)      WORK_DIR="${2:?--work-dir needs a value}"; shift 2 ;;
        --work-on-disk)  WORK_ON_DISK="yes"; shift ;;
        --in-ram)        WORK_ON_DISK="no"; shift ;;
        --jellyfish-bin) JELLYFISH_BIN="${2:?--jellyfish-bin needs a value}"; shift 2 ;;
        --clean)         CLEAN=1; shift ;;
        --force)         FORCE=1; shift ;;
        --dry-run)       DRY_RUN=1; shift ;;
        -h|--help)       usage; exit 0 ;;
        *)               usage; die "unknown option: $1" ;;
    esac
done

[[ "$THREADS"       =~ ^[0-9]+$ ]] || die "--threads must be an integer"
[[ "$KMER"          =~ ^[0-9]+$ ]] || die "--kmer must be an integer"
[[ "$MINIMIZER"     =~ ^[0-9]+$ ]] || die "--minimizer must be an integer"
[[ "$HLL_PRECISION" =~ ^[0-9]+$ ]] || die "--hll-precision must be an integer"
[[ "$MIN_FREE_GIB"  =~ ^[0-9]+$ ]] || die "--min-free-gib must be an integer"
(( THREADS >= 1 ))      || die "--threads must be >= 1"
(( KMER >= 2 ))         || die "--kmer must be >= 2"
# Hard ceiling from Kraken 1's 2-bit-packed 64-bit k-mer representation.
(( KMER <= 31 ))        || die "--kmer must be <= 31 (KrakenUniq inherits Kraken 1's 64-bit packed k-mer limit)"
# Note STRICTLY less-than, unlike kraken2-build which allows minimizer == k.
(( MINIMIZER < KMER ))  || die "--minimizer ($MINIMIZER) must be STRICTLY less than --kmer ($KMER)"
(( MINIMIZER >= 1 ))    || die "--minimizer must be >= 1"

DB_DIR="${BTF_DB}/krakenuniq/standard"
STAMP_DIR="${DB_DIR}/.build_stamps"
TAX_DIR="${DB_DIR}/taxonomy"
LIB_DIR="${DB_DIR}/library"

# Hard guard: a mistyped --btf-db must not be able to write anywhere else. The
# same posture as build_kraken2_db.sh, and for a stronger reason here: the shared
# volume holds the only databases on this box that currently classify anything.
case "$DB_DIR" in
    */krakenuniq/standard) : ;;
    *) die "refusing to build at '${DB_DIR}' - expected a path ending in /krakenuniq/standard" ;;
esac
[[ "$DB_DIR" == "${BTF_DB}/krakenuniq/"* ]] \
    || die "refusing to build at '${DB_DIR}' - it is not under ${BTF_DB}/krakenuniq"
case "$DB_DIR" in
    /media/36tb/db/*) die "refusing to write under /media/36tb/db - that volume is shared and read-only for this project" ;;
esac
# Never write inside the shared taxdump, and never near the stale 2022 Kraken2
# databases, whatever --btf-db was set to.
for protected in "$TAXDUMP_SRC" \
                 "/media/36tb/db/taxdump" \
                 "/media/36tb/db/kraken2" \
                 "/media/36tb/db/nt" \
                 "/media/36tb/db/sylph"; do
    if [[ "$DB_DIR" == "$protected" || "$DB_DIR" == "$protected"/* ]]; then
        die "refusing to write inside the pre-existing shared resource ${protected}"
    fi
done

if (( HLL_PRECISION < 10 || HLL_PRECISION > 18 )); then
    die "--hll-precision must be between 10 and 18 (krakenuniq classifier accepts 10-18; got ${HLL_PRECISION})"
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

INSTALL_HINT="Install KrakenUniq into its own environment:
         ${CONDA_BIN}/mamba create -y -n krakenuniq -c bioconda -c conda-forge krakenuniq
         conda activate krakenuniq
       Use a dedicated environment on purpose: krakenuniq needs Jellyfish 1.x for
       the build, while other tools in a shared env will want Jellyfish 2.x, and
       the two cannot both be first on PATH."

KU_BUILD="$(resolve_tool krakenuniq-build || true)"
if [[ -z "$KU_BUILD" ]]; then
    msg="krakenuniq-build not found on PATH or in ${CONDA_BIN}.
       ${INSTALL_HINT}"
    if (( DRY_RUN == 1 )); then warn "$msg"; KU_BUILD="krakenuniq-build"; else die "$msg"; fi
fi
KU_DOWNLOAD="$(resolve_tool krakenuniq-download || true)"
if [[ -z "$KU_DOWNLOAD" ]]; then
    msg="krakenuniq-download not found on PATH or in ${CONDA_BIN}.
       It ships with krakenuniq and is what actually fetches the libraries.
       ${INSTALL_HINT}"
    if (( DRY_RUN == 1 )); then warn "$msg"; KU_DOWNLOAD="krakenuniq-download"; else die "$msg"; fi
fi
# The classifier itself is not needed to build, but if it is missing the index
# has no consumer, and finding that out now is cheaper than after 20 hours.
if ! resolve_tool krakenuniq >/dev/null; then
    warn "krakenuniq (the classifier) is not on PATH. The index will build, but nothing can query it yet.
       ${INSTALL_HINT}"
fi

# Jellyfish 1 check. krakenuniq-build fails at the k-mer counting step, well into
# the run, if it gets Jellyfish 2 - so check up front.
check_jellyfish() {
    local jf="$JELLYFISH_BIN" ver
    if [[ -z "$jf" ]]; then
        # krakenuniq-build prefers its own bundled jellyfish-install/bin/jellyfish
        # next to the script; if that exists we trust it and say nothing.
        local kudir
        kudir="$(dirname "$(readlink -f "$KU_BUILD" 2>/dev/null || printf '%s' "$KU_BUILD")")"
        if [[ -x "${kudir}/jellyfish-install/bin/jellyfish" ]]; then
            log "jellyfish: using krakenuniq's bundled ${kudir}/jellyfish-install/bin/jellyfish"
            return 0
        fi
        jf="$(resolve_tool jellyfish || true)"
    fi
    if [[ -z "$jf" ]]; then
        warn "no jellyfish found and krakenuniq ships none next to krakenuniq-build.
       The build WILL fail at the k-mer counting step. Install krakenuniq >=0.7.3
       (which bundles Jellyfish 1) or pass --jellyfish-bin /path/to/jellyfish1."
        return 0
    fi
    ver="$("$jf" --version 2>&1 | head -1 || true)"
    case "$ver" in
        *" 1."*|*jellyfish\ 1*) log "jellyfish: ${jf} (${ver}) - version 1, correct for krakenuniq-build" ;;
        *) warn "jellyfish at ${jf} reports '${ver}'.
       krakenuniq-build REQUIRES Jellyfish version 1.x and will fail with 2.x.
       Install krakenuniq >=0.7.3 (bundles Jellyfish 1) or pass
       --jellyfish-bin /path/to/jellyfish1." ;;
    esac
}

if (( DRY_RUN == 0 )); then
    # krakenuniq-download uses these for the RefSeq pull.
    for dep in perl wget; do
        command -v "$dep" >/dev/null 2>&1 \
            || warn "krakenuniq-download normally needs '${dep}' on PATH; the download may fail"
    done
    command -v rsync >/dev/null 2>&1 \
        || warn "rsync not found; krakenuniq-download --rsync is unavailable (ftp/http fallback will be used)"
fi

# ---------------------------------------------------------------------------
# Pre-flight free space
# ---------------------------------------------------------------------------
# The arithmetic, so the number is arguable instead of magic. All GiB.
#
#   library/ (the default library set, Complete_Genome archaea+bacteria plus all
#            viral/fungi/protozoa plus human plus contaminants)
#            ......................................................  200-350
#   transient k-mer files: jellyfish counts the whole library at k=31, then the
#            counts are sorted into database0/database.kdb. KrakenUniq stores a
#            2-bit-packed 64-bit k-mer (8 B) plus a 4 B taxid per DISTINCT
#            k-mer, i.e. ~12 B. A PlusPF-scale library holds on the order of
#            2-4e10 distinct 31-mers -> 220-450 GiB, and the unsorted jellyfish
#            output and the sorted copy coexist briefly.
#            ......................................................  450-800
#   final database.kdb + database.idx + taxDB .................  250-450
#
# Peak is not the sum of all three (the jellyfish output is dropped once sorted),
# but library + sorted k-mers + final index do overlap. 1400 GiB is the
# conservative peak for the default library set.
#
# COMBINED with the other stage-1 nucleotide index on the same NVMe: because a
# Kraken2 database CANNOT be reused by KrakenUniq (see the header), the finished
# Kraken2 PlusPF index (~100-120 GiB) is an ADDITIONAL resident cost, and its own
# build needs ~700 GiB transient. Steady state for both: ~400-600 GiB. Peak if
# both are built in the same window: ~2100 GiB. There is 5.9 TB free here, so
# both fit; the point of stating it is that "we already have Kraken2" is not a
# disk argument for skipping this build.
check_free_space() {
    local probe="$BTF_DB" avail_gib
    while [[ ! -d "$probe" && "$probe" != "/" ]]; do probe="$(dirname "$probe")"; done
    avail_gib="$(df -BG --output=avail "$probe" | tail -1 | tr -dc '0-9')"
    [[ -n "$avail_gib" ]] || die "could not determine free space for ${probe}"
    log "free space under ${probe}: ${avail_gib} GiB"
    log "required (transient peak): ${MIN_FREE_GIB} GiB"
    log "  breakdown: ~200-350 GiB RefSeq FASTA in library/,"
    log "             ~450-800 GiB transient jellyfish + sorted k-mer files,"
    log "             ~250-450 GiB final database.kdb + database.idx + taxDB,"
    log "             ~1 GiB copied taxonomy (nodes/names/merged/delnodes.dmp)."
    log "  COMBINED with Kraken2: a Kraken2 DB cannot be reused by KrakenUniq,"
    log "             so PlusPF's ~100-120 GiB index is an additional resident"
    log "             cost (+~700 GiB transient while IT builds). Both resident:"
    log "             ~400-600 GiB. Both building at once: ~2100 GiB."
    if [[ -n "$WORK_DIR" ]]; then
        local wprobe="$WORK_DIR" wavail
        while [[ ! -d "$wprobe" && "$wprobe" != "/" ]]; do wprobe="$(dirname "$wprobe")"; done
        wavail="$(df -BG --output=avail "$wprobe" | tail -1 | tr -dc '0-9')"
        log "  --work-dir ${WORK_DIR} lives on ${wprobe} with ${wavail:-?} GiB free"
        log "  (the transient 450-800 GiB lands there instead of on ${probe})"
    fi
    if (( avail_gib < MIN_FREE_GIB )); then
        die "insufficient free space: ${avail_gib} GiB available, ${MIN_FREE_GIB} GiB required.
       Options:
         * put the transient files on another fast volume with --work-dir, or
         * shrink the library set with --libraries (dropping
           refseq/bacteria/Complete_Genome is by far the biggest saving, at the
           cost of exactly the sensitivity this screen exists for), or
         * build with --clean so library/ is reclaimed immediately afterwards, or
         * lower the requirement with --min-free-gib if you are certain.
       Nothing was downloaded."
    fi
}

# ---------------------------------------------------------------------------
# Stamps - one per library spec, so a restart resumes instead of re-pulling
# hundreds of GB. Same pattern as build_kraken2_db.sh.
# ---------------------------------------------------------------------------
stamp_name() {
    # Library specs contain '/' and '=' (refseq/vertebrate_mammalian/Any/species_taxid=9606)
    # which cannot appear in a filename, so flatten them reversibly enough to read.
    printf '%s\n' "$1" | tr '/=,' '___'
}
stamp()      { printf '%s\n' "${STAMP_DIR}/$1.done"; }
have_stamp() { [[ -f "$(stamp "$1")" ]]; }
set_stamp()  { (( DRY_RUN == 1 )) || { mkdir -p "$STAMP_DIR"; date -u +%Y-%m-%dT%H:%M:%SZ > "$(stamp "$1")"; }; }

# ---------------------------------------------------------------------------
# Taxonomy
# ---------------------------------------------------------------------------
# Does krakenuniq-build accept an external taxonomy?
#
# Yes - and unlike kraken2-build it has a real flag for it: `--taxonomy-dir DIR`
# (it sets KRAKEN_TAXONOMY_DIR, which build_db.sh uses in place of DBDIR/taxonomy).
# The build reads only names.dmp and nodes.dmp from that directory, and writes
# its derived `taxDB` into the DATABASE directory, not into the taxonomy one.
#
# So why COPY instead of pointing --taxonomy-dir at /media/36tb/db/taxdump?
# Because build_db.sh contains this fallback: if names.dmp or nodes.dmp is
# missing from TAXONOMY_DIR, it cd's INTO that directory and downloads a fresh
# taxdump there. One bad path, one partial file, and the shared taxdump that
# resolve_watchlist.py pins for every pack gets overwritten - a silent
# watchlist-resolution change across the whole project. Copying four files costs
# about a gigabyte and removes that class of accident entirely. Same reasoning,
# and the same four files, as build_kraken2_db.sh.
#
# VERIFY: --taxonomy-dir and the build_db.sh download fallback were read from
# krakenuniq master (0.7.x). If a future version rejects an external taxonomy,
# the fix is `krakenuniq-download --db <db> taxonomy`, which writes into
# <db>/taxonomy - still never into the shared volume.
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
            || die "required taxonomy file missing or empty: ${TAXDUMP_SRC}/${f}
       krakenuniq-build needs both nodes.dmp and names.dmp. Either point
       --taxdump at a complete unpacked taxdump, or fetch a private copy with
         ${KU_DOWNLOAD} --db ${DB_DIR} taxonomy
       (which writes into ${DB_DIR}/taxonomy and never touches the shared volume)."
    done

    require_fresh_taxdump "$TAXDUMP_SRC" "$MAX_TAXDUMP_AGE_DAYS"
    log "copying taxonomy from ${TAXDUMP_SRC} (copy, never symlink - see comment)"
    run mkdir -p "$TAX_DIR"
    run cp -f "${TAXDUMP_SRC}/nodes.dmp" "${TAX_DIR}/nodes.dmp"
    run cp -f "${TAXDUMP_SRC}/names.dmp" "${TAX_DIR}/names.dmp"
    # merged.dmp/delnodes.dmp are not read by the build, but keeping them beside
    # the index records which taxonomy snapshot it was built from - which is what
    # lets a later watchlist mismatch be diagnosed instead of guessed at.
    for extra in merged.dmp delnodes.dmp; do
        if [[ -s "${TAXDUMP_SRC}/${extra}" ]]; then
            run cp -f "${TAXDUMP_SRC}/${extra}" "${TAX_DIR}/${extra}"
        else
            warn "${TAXDUMP_SRC}/${extra} not found; provenance record will be incomplete"
        fi
    done
    if (( DRY_RUN == 0 )); then
        printf 'taxdump_source=%s\ntaxdump_copied_utc=%s\n' \
            "$TAXDUMP_SRC" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "${TAX_DIR}/BTF_PROVENANCE.txt"
    fi
    set_stamp taxonomy
}

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------
# One krakenuniq-download call per spec, each stamped. krakenuniq-download --db
# writes to <db>/library/<domain>/<assembly_level>/ and drops the *.map files
# that build_db.sh later concatenates into seqid2taxid.map, so the layout must be
# left exactly as the downloader made it.
download_libraries() {
    local spec key
    for spec in $LIBRARIES; do
        key="library_$(stamp_name "$spec")"
        if have_stamp "$key" && (( FORCE == 0 )); then
            skip "library already downloaded: ${spec}"
            continue
        fi
        case "$spec" in
            *bacteria*) log "downloading library '${spec}' - the slowest single step, hours on its own" ;;
            *)          log "downloading library '${spec}'" ;;
        esac
        # --dust masks low-complexity regions at download time. KrakenUniq has no
        # separate masking step at build, and unmasked low-complexity k-mers are a
        # direct source of the spurious hits this index exists to suppress.
        # VERIFY: --dust/-D is a krakenuniq-download flag in 0.7.x; confirm it is
        # still present, and that it has not moved to krakenuniq-build.
        run "$KU_DOWNLOAD" --db "$DB_DIR" --threads "$THREADS" --dust "$spec"
        set_stamp "$key"
    done
    if (( DRY_RUN == 0 )); then
        [[ -d "$LIB_DIR" ]] || die "no library/ directory after download - krakenuniq-download wrote nothing to ${DB_DIR}"
        # A silently empty library produces a syntactically valid but useless
        # index, which is the worst failure mode for a screen.
        if [[ -z "$(find -L "$LIB_DIR" -type f \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' -o -name '*.fna.gz' \) -print -quit 2>/dev/null)" ]]; then
            die "library/ under ${DB_DIR} contains no FASTA. Refusing to build an empty index."
        fi
        log "library/ now holds $(du -sh "$LIB_DIR" 2>/dev/null | cut -f1)"
    fi
}

# ---------------------------------------------------------------------------
# RAM vs disk
# ---------------------------------------------------------------------------
# build_db.sh's default is the RAM build (it passes -M, "minimize disk writes").
# --work-on-disk inverts that. The in-RAM build wants roughly as much memory as
# the finished index; if that exceeds what the box has, the run dies late and
# expensively, so estimate and choose.
#
# The estimate: distinct 31-mers scale with the non-redundant sequence content of
# the library, and KrakenUniq spends ~12 B per distinct k-mer. Measured across
# RefSeq-scale builds the finished index lands near 1.0-1.3x the size of the
# (uncompressed-equivalent) library. Using 1.0x is deliberately optimistic about
# the index and therefore CONSERVATIVE about switching to disk.
# VERIFY: this ratio is an empirical rule of thumb, not a documented figure.
# Recalibrate it once this box has built one index: the honest number is
# du(database.kdb) / du(library).
decide_work_on_disk() {
    local ram_gib lib_gib est_gib budget_gib
    ram_gib="$(awk '/^MemTotal:/ {printf "%d", $2/1048576}' /proc/meminfo)"
    # Leave headroom: an in-RAM build that uses every last page thrashes and the
    # OOM killer picks the build, not the estimate.
    budget_gib=$(( ram_gib * 80 / 100 ))

    if [[ "$WORK_ON_DISK" == "yes" ]]; then
        log "memory policy: ON-DISK build (forced with --work-on-disk). Slower, but cannot be OOM-killed."
        return 0
    fi
    if [[ "$WORK_ON_DISK" == "no" ]]; then
        warn "memory policy: IN-RAM build (forced with --in-ram). If the index exceeds ~${budget_gib} GiB this run will be OOM-killed hours in."
        return 0
    fi

    if [[ ! -d "$LIB_DIR" ]]; then
        # Nothing measured yet (dry run, or download skipped). Do not guess.
        WORK_ON_DISK="no"
        log "memory policy: library/ not present, cannot estimate; defaulting to the IN-RAM build (${ram_gib} GiB RAM, ~${budget_gib} GiB usable)"
        return 0
    fi
    lib_gib="$(du -s --block-size=1G -L "$LIB_DIR" 2>/dev/null | cut -f1)"
    lib_gib="${lib_gib:-0}"
    est_gib="$lib_gib"
    log "memory policy: library/ is ${lib_gib} GiB -> estimated index ~${est_gib} GiB; usable RAM ~${budget_gib} of ${ram_gib} GiB"
    if (( est_gib > budget_gib )); then
        WORK_ON_DISK="yes"
        warn "estimated index (~${est_gib} GiB) exceeds the ~${budget_gib} GiB RAM budget."
        warn "DEGRADING to the on-disk build (krakenuniq-build --work-on-disk)."
        warn "This is slower - plan on roughly 1.5-3x the wall time - but it completes"
        warn "instead of being OOM-killed. Override with --in-ram if you disagree with"
        warn "the estimate, or reduce --libraries. Recorded in the build manifest."
    else
        WORK_ON_DISK="no"
        log "estimate fits in RAM; using the faster in-RAM build. Force the safe route with --work-on-disk."
    fi
}

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
build_db() {
    if [[ -s "${DB_DIR}/database.kdb" && -s "${DB_DIR}/database.idx" && -s "${DB_DIR}/taxDB" ]] \
       && (( FORCE == 0 )); then
        skip "database already built ($(du -sh "$DB_DIR" 2>/dev/null | cut -f1) in ${DB_DIR})"
        return 0
    fi

    decide_work_on_disk

    local args=(--db "$DB_DIR" --build --threads "$THREADS"
                --kmer-len "$KMER" --minimizer-len "$MINIMIZER"
                --taxonomy-dir "$TAX_DIR")
    [[ "$WORK_ON_DISK" == "yes" ]] && args+=(--work-on-disk)
    [[ -n "$JELLYFISH_BIN" ]] && args+=(--jellyfish-bin "$JELLYFISH_BIN")
    # Deliberately NOT passing --max-db-size. It does not truncate, it SHRINKS the
    # k-mer set to fit, discarding k-mers - and a discarded k-mer is a detection
    # that silently never happens. The same reasoning as build_kraken2_db.sh.
    # Deliberately NOT passing --uid-database: it roughly doubles build time and
    # size, and the LCA database plus HLL cardinalities already answer the
    # question this pipeline asks.

    warn "the k-mer count, sort and LCA phases now run for many hours"
    if [[ -n "$WORK_DIR" ]]; then
        # build_db.sh cd's into the DB directory and writes its transient files
        # there; TMPDIR is the only lever that moves the large external sorts.
        # VERIFY: this relocates the sort spill files, not database0.*; if the
        # database volume is genuinely too tight, put the whole DB_DIR on the
        # roomier volume and symlink instead of relying on this.
        log "exporting TMPDIR=${WORK_DIR} for the external sorts"
        run mkdir -p "$WORK_DIR"
        (( DRY_RUN == 1 )) || export TMPDIR="$WORK_DIR"
    fi
    run "$KU_BUILD" "${args[@]}"

    if (( DRY_RUN == 0 )); then
        for f in database.kdb database.idx taxDB; do
            [[ -s "${DB_DIR}/${f}" ]] \
                || die "build finished but ${DB_DIR}/${f} is missing or empty - do not point the pipeline at this directory"
        done
    fi
}

# ---------------------------------------------------------------------------
# Classify-time defaults
# ---------------------------------------------------------------------------
# --hll-precision is the knob that decides how accurate the distinct-k-mer
# counts - the whole reason this index exists - actually are. It belongs to the
# classifier, so it cannot be baked in; recording it here keeps the build and the
# Nextflow module from drifting apart.
#
# The tradeoff: HyperLogLog uses 2^p registers, so p is a straight trade of
# memory and per-taxon bookkeeping against relative error, which falls roughly as
# 1.04/sqrt(2^p):
#     p=10   1024 registers   ~3.3% relative error   cheapest
#     p=12   4096 registers   ~1.6% relative error   <- default
#     p=14  16384 registers   ~0.8% relative error
#     p=18  262144 registers  ~0.2% relative error   costliest
# Registers are allocated PER TAXON observed in a sample, so the cost is
# 2^p x (register bytes) x (taxa seen), and a metagenome can touch tens of
# thousands of taxa - which is why the high end is not free.
#
# Why 12 is the right operating point here: the decision this pipeline makes from
# cardinality is a threshold comparison - "did this taxon's distinct-k-mer count
# clear the gate, or are these reads stacked on one locus?" - and the real signal
# separation is an order of magnitude, not a few percent. A 1.6% error on a
# quantity whose decision boundary is 10x away changes no call, while p=18 would
# multiply the sketch memory by 64 for accuracy nothing reads. Raise p only if a
# future calibration step (roadmap item 2, empirical null / calibrated FDR) starts
# using cardinality as a continuous statistic rather than a gate.
write_classify_defaults() {
    local out="${DB_DIR}/btf_classify_defaults.env"
    if (( DRY_RUN == 1 )); then
        printf '[dry-run] write %s (KRAKENUNIQ_HLL_PRECISION=%s)\n' "$out" "$HLL_PRECISION" >&2
        return 0
    fi
    cat > "$out" <<EOF
# Written by build_krakenuniq_db.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ).
# These are CLASSIFY-time parameters. krakenuniq-build accepts none of them;
# --hll-precision belongs to the krakenuniq classifier (valid range 10-18).
# Source this in the Nextflow krakenuniq process so query and build agree.
KRAKENUNIQ_DB=${DB_DIR}
KRAKENUNIQ_HLL_PRECISION=${HLL_PRECISION}
# Built with k=${KMER}, minimizer=${MINIMIZER}. The classifier reads k from the
# index, so these are recorded for provenance, not to be passed back in.
KRAKENUNIQ_BUILD_KMER=${KMER}
KRAKENUNIQ_BUILD_MINIMIZER=${MINIMIZER}
# If RAM at query time is tight, prefer --preload-size (e.g. 32G) over --preload:
# it streams the index in chunks and is what makes a ~400 GB database queryable
# on a machine that cannot hold it. VERIFY the current flag name for your version.
EOF
    log "recorded classify-time defaults in ${out}"
}

write_manifest() {
    local out="${DB_DIR}/BTF_BUILD_MANIFEST.txt"
    if (( DRY_RUN == 1 )); then
        printf '[dry-run] write %s\n' "$out" >&2
        return 0
    fi
    {
        printf 'tool=krakenuniq\n'
        printf 'built_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'db_dir=%s\n' "$DB_DIR"
        printf 'krakenuniq_build=%s\n' "$KU_BUILD"
        printf 'krakenuniq_version=%s\n' "$("$KU_BUILD" --version 2>&1 | head -1 || echo unknown)"
        printf 'kmer=%s\nminimizer=%s\n' "$KMER" "$MINIMIZER"
        printf 'hll_precision_for_classify=%s\n' "$HLL_PRECISION"
        printf 'work_on_disk=%s\n' "$WORK_ON_DISK"
        printf 'threads=%s\n' "$THREADS"
        printf 'libraries=%s\n' "$LIBRARIES"
        printf 'taxdump_source=%s\n' "$TAXDUMP_SRC"
        printf 'host=%s\n' "$(hostname)"
    } > "$out"
    log "wrote build manifest ${out}"
}

# ---------------------------------------------------------------------------
# Plan / main
# ---------------------------------------------------------------------------
log "=================== bioThreatFinder krakenuniq build ==================="
log "database dir : ${DB_DIR}"
log "libraries    : ${LIBRARIES}"
log "k / minimizer: ${KMER} / ${MINIMIZER}   (KrakenUniq caps k at 31)"
log "threads      : ${THREADS}"
log "taxonomy     : copy nodes/names/merged/delnodes.dmp from ${TAXDUMP_SRC}"
log "hll-precision: ${HLL_PRECISION} (classify-time; recorded, not passed to the build)"
case "$WORK_ON_DISK" in
    auto) log "memory policy: auto - decided from the measured library size vs 80% of RAM" ;;
    yes)  log "memory policy: on-disk build (forced with --work-on-disk)" ;;
    no)   log "memory policy: in-RAM build (forced with --in-ram)" ;;
esac
log "work dir     : ${WORK_DIR:-<inside the database directory>}"
log "format       : Kraken 1 lineage - a Kraken2 DB CANNOT be reused here"
log "untouched    : /media/36tb/db/** (shared, read-only for this project)"
log "EXPECT THIS TO TAKE MANY HOURS (commonly 12-30 h). Run under nohup, tmux or screen."
log "======================================================================="

if (( DRY_RUN == 1 )); then
    log "--dry-run: printing the plan, downloading and building nothing"
fi

check_jellyfish
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

prepare_taxonomy
download_libraries
build_db
write_classify_defaults
write_manifest

if (( CLEAN == 1 )); then
    # krakenuniq-build --clean removes the intermediates the classifier does not
    # read. VERIFY whether the version in use also removes library/ - in 0.7.x it
    # does not, so remove it explicitly below.
    warn "--clean: removing build intermediates; a future rebuild will re-download everything"
    run "$KU_BUILD" --db "$DB_DIR" --clean
    if [[ -d "$LIB_DIR" ]]; then
        warn "--clean: also removing ${LIB_DIR} ($(du -sh "$LIB_DIR" 2>/dev/null | cut -f1))"
        run rm -rf "$LIB_DIR"
        # Drop the library stamps too, or a later --force run would "resume" past
        # a download whose output no longer exists.
        run find "$STAMP_DIR" -name 'library_*.done' -delete
    fi
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo >&2
log "=============================== SUMMARY ==============================="
if (( DRY_RUN == 1 )); then
    log "dry run complete - nothing was downloaded or built"
else
    for f in database.kdb database.idx taxDB seqid2taxid.map; do
        if [[ -e "${DB_DIR}/${f}" ]]; then
            log "  $(du -h "${DB_DIR}/${f}" | cut -f1)	${f}"
        else
            warn "  MISSING	${f}"
        fi
    done
    log "  total: $(du -sh "$DB_DIR" 2>/dev/null | cut -f1) in ${DB_DIR}"
    if [[ -d "$LIB_DIR" && $CLEAN -eq 0 ]]; then
        log "  library/ still holds $(du -sh "$LIB_DIR" 2>/dev/null | cut -f1); rerun with --clean to reclaim it"
    fi
    log "free space remaining: $(df -h --output=avail "$BTF_DB" | tail -1 | tr -d ' ')"
    log "Point screen.krakenuniq_db at: ${DB_DIR}"
    log "Classify with --report-file and --hll-precision ${HLL_PRECISION}; the report's"
    log "  'kmers' (distinct k-mer cardinality) and 'dup' columns are what stage 1 gates on."
fi
log "======================================================================="
