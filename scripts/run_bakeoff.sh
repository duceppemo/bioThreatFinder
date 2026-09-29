#!/usr/bin/env bash
#
# Tool bake-off driver: run the pipeline several times over the SAME simulated
# data, changing only the configuration, then compare accuracy AND cost.
#
# WHY THIS EXISTS
# Tool choice in metagenomic pathogen detection is usually argued from papers.
# Papers benchmark on their authors' data, with their authors' thresholds, for a
# community-profiling question rather than an "is organism X present" question.
# The only answer that binds for THIS pipeline, THIS watchlist and THIS specimen
# type is a measurement on identical input.
#
# The cost half matters as much as the accuracy half. KrakenUniq costs ~1400 GiB
# transient disk, a ~250-450 GiB resident index, a Jellyfish 1.x build dependency
# and 12+ hours to build. Kraken2's --report-minimizer-data already yields distinct
# minimizer counts and a duplicity ratio from an index we build anyway. Whether the
# HyperLogLog cardinality is worth that difference is an empirical question, and
# this script is how it gets answered rather than asserted.
#
# Usage:
#   run_bakeoff.sh --batch v1 [--simdir DIR] [--outdir DIR] [--profile P]
#                  [--variants "name1 name2 ..."] [--list] [--dry-run]
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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

BATCH=""
SIMDIR="/media/8tb/bioThreatFinder/simulated"
OUTDIR=""
PROFILE="conda,workstation"
PACK="human"
SCENARIOS="${REPO_DIR}/benchmarks/scenarios/human_core.yaml"
NEXTFLOW="${NEXTFLOW:-$(command -v nextflow || echo /home/bioinfo/prog/nextflow)}"
WANTED=""
DRY_RUN=0
LIST_ONLY=0

# ---------------------------------------------------------------------------
# The variant matrix.
#
# Each entry is "name|extra nextflow args|one-line rationale".
#
# Keep these MINIMAL PAIRS wherever possible. A variant that changes two things at
# once cannot attribute the difference to either, which is the commonest way a
# tool comparison ends up unable to answer the question it was run for.
# ---------------------------------------------------------------------------
# NOTE ON --skip_krakenuniq: it DEFAULTS TO TRUE (nextflow.config), so omitting it
# does not enable KrakenUniq - `--skip_krakenuniq false` does. The matrix originally
# read `plus_krakenuniq|--skip_metabuli`, which left the arm off and made that variant
# a byte-for-byte duplicate of the baseline: the comparison it exists to make returned
# a guaranteed null difference, readable as "KrakenUniq adds nothing". decisions.md D2
# names this exact variant as what would reverse the Kraken2-primary decision, so D2
# was unfalsifiable by its own stated test. `all_arms` had the same hole.
VARIANTS=(
  "kraken2_minimizer|--skip_krakenuniq true --skip_metabuli|Baseline. Kraken2 with --report-minimizer-data (distinct minimizers + duplicity) plus Bracken and sylph. The cheapest configuration that still has distinct-k-mer-style evidence."
  "plus_krakenuniq|--skip_krakenuniq false --skip_metabuli|Baseline + KrakenUniq. ISOLATES the question: does HyperLogLog cardinality beat Kraken2 minimizer counts by enough to justify a second index? Needs a built KrakenUniq index."
  "plus_metabuli|--skip_krakenuniq true|Baseline + Metabuli. Isolates the value of the DNA+amino-acid arm, which should show up on the divergent_* scenarios specifically."
  "all_arms|--skip_krakenuniq false|Everything on, KrakenUniq included. Upper bound on sensitivity and on cost."
  "no_sylph|--skip_krakenuniq true --skip_metabuli --skip_sylph|Baseline minus sylph, to check whether ANI/containment is contributing anything the k-mer arms do not."
  "no_markers|--skip_markers|All arms except KrakenUniq, marker gating OFF. Quantifies what marker gating buys - expect confirm PPV to fall and false confirmations to rise."
)

usage() {
  cat <<EOF
Tool bake-off driver for bioThreatFinder.

Runs the pipeline once per variant over the same simulated batch, then calls
bin/evaluate_benchmark.py with every variant's findings AND Nextflow trace, so the
comparison reports accuracy and cost together.

Required:
  --batch NAME          simulated batch produced by scripts/simulate_benchmark.py
                        (expects <simdir>/<batch>/samplesheet.csv and manifest.tsv)

Options:
  --simdir DIR          where batches live            [${SIMDIR}]
  --outdir DIR          results root                  [<simdir>/<batch>/bakeoff]
  --profile P           nextflow -profile             [${PROFILE}]
  --pack NAME           domain pack                   [${PACK}]
  --scenarios FILE      scenario definitions          [benchmarks/scenarios/human_core.yaml]
  --variants "a b c"    run only these named variants [all]
  --list                print the variant matrix and exit
  --dry-run             print the commands without running them
  -h, --help

Variants available: $(printf '%s ' "${VARIANTS[@]%%|*}")

NOTE: every variant must run over the SAME data. This script enforces that by
reusing one samplesheet; do not point variants at different batches.
EOF
}

log()  { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*" >&2; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch)     BATCH="${2:?--batch needs a value}"; shift 2 ;;
    --simdir)    SIMDIR="${2:?}"; shift 2 ;;
    --outdir)    OUTDIR="${2:?}"; shift 2 ;;
    --profile)   PROFILE="${2:?}"; shift 2 ;;
    --pack)      PACK="${2:?}"; shift 2 ;;
    --scenarios) SCENARIOS="${2:?}"; shift 2 ;;
    --variants)  WANTED="${2:?}"; shift 2 ;;
    --list)      LIST_ONLY=1; shift ;;
    --dry-run)   DRY_RUN=1; shift ;;
    -h|--help)   usage; exit 0 ;;
    *)           die "unknown argument: $1 (try --help)" ;;
  esac
done

if [[ "$LIST_ONLY" == 1 ]]; then
  printf '%-20s %-52s %s\n' "VARIANT" "EXTRA ARGS" "RATIONALE"
  for v in "${VARIANTS[@]}"; do
    IFS='|' read -r n a r <<<"$v"
    printf '%-20s %-52s %s\n' "$n" "$a" "$r"
  done
  exit 0
fi

[[ -n "$BATCH" ]] || { usage; die "--batch is required"; }

BATCH_DIR="${SIMDIR}/${BATCH}"
SAMPLESHEET="${BATCH_DIR}/samplesheet.csv"
MANIFEST="${BATCH_DIR}/manifest.tsv"
OUTDIR="${OUTDIR:-${BATCH_DIR}/bakeoff}"
WATCHLIST="${REPO_DIR}/packs/${PACK}/watchlist.resolved.tsv"

[[ -s "$SAMPLESHEET" ]] || die "no samplesheet at ${SAMPLESHEET}. Generate the batch first with scripts/simulate_benchmark.py --batch ${BATCH}"
[[ -s "$MANIFEST" ]]    || die "no manifest at ${MANIFEST} - the evaluator needs it for ground truth"
[[ -s "$SCENARIOS" ]]   || die "no scenario file at ${SCENARIOS}"
[[ -s "$WATCHLIST" ]]   || die "no resolved watchlist at ${WATCHLIST}. Run scripts/resolve_watchlist.py first."
[[ -x "$NEXTFLOW" ]] || command -v "$NEXTFLOW" >/dev/null 2>&1 || die "nextflow not found (set \$NEXTFLOW)"

selected=()
for v in "${VARIANTS[@]}"; do
  name="${v%%|*}"
  if [[ -z "$WANTED" ]] || [[ " $WANTED " == *" $name "* ]]; then
    selected+=("$v")
  fi
done
[[ ${#selected[@]} -gt 0 ]] || die "no variants selected (check --variants against --list)"

log "batch       : ${BATCH}  ($(grep -c . "$SAMPLESHEET" | awk '{print $1-1}') samples)"
log "variants    : ${#selected[@]}"
log "outdir      : ${OUTDIR}"
mkdir -p "$OUTDIR"

EVAL_ARGS=()
for v in "${selected[@]}"; do
  IFS='|' read -r name extra rationale <<<"$v"
  [[ "$extra" == "(none)" ]] && extra=""
  vdir="${OUTDIR}/${name}"
  findings="${vdir}/score_evidence/findings.tsv"

  # Resume-friendly: a completed variant is skipped, so a bake-off interrupted
  # after four of six runs does not restart from scratch.
  # `-s` is true for a header-only file, which would skip a variant that produced
  # no findings at all and then feed an empty table to the comparison.
  if [[ -s "$findings" ]] && (( $(grep -vc '^$' "$findings" 2>/dev/null || echo 0) > 1 )); then
    log "SKIP ${name} (findings already present)"
  else
    log "RUN  ${name}  ${extra:-<all arms>}"
    log "     ${rationale}"
    # Word splitting on $extra is INTENTIONAL: it carries several flags
    # ("--skip_krakenuniq --skip_metabuli") that must become separate array elements.
    # Quoting it would pass them to nextflow as one argument.
    # shellcheck disable=SC2086,SC2206
    cmd=("$NEXTFLOW" run "$REPO_DIR"
         -profile "$PROFILE"
         --input "$SAMPLESHEET"
         --pack "$PACK"
         --outdir "$vdir"
         -work-dir "${OUTDIR}/work/${name}"
         $extra)
    if [[ "$DRY_RUN" == 1 ]]; then
      printf '  %q ' "${cmd[@]}"; printf '\n'
    else
      # A single variant failing must not abandon the whole bake-off; record it and
      # carry on, so five good comparisons are not lost to one bad configuration.
      if ! "${cmd[@]}" >"${OUTDIR}/${name}.log" 2>&1; then
        log "WARNING: variant ${name} FAILED - see ${OUTDIR}/${name}.log"
        log "         it will be omitted from the comparison"
        continue
      fi
    fi
  fi

  [[ "$DRY_RUN" == 1 ]] && continue
  if [[ ! -s "$findings" ]] || (( $(grep -vc '^$' "$findings" 2>/dev/null || echo 0) <= 1 )); then
    log "WARNING: ${name} produced no findings (file absent or header only), omitting"
    continue
  fi

  EVAL_ARGS+=(--variant "${name}=${findings}")
  # Newest trace for this variant, for the cost columns.
  trace="$(ls -t "${vdir}"/pipeline_info/trace-*.txt 2>/dev/null | head -1 || true)"
  [[ -n "$trace" ]] && EVAL_ARGS+=(--variant-trace "${name}=${trace}")
done

if [[ "$DRY_RUN" == 1 ]]; then
  log "dry run - nothing executed"
  exit 0
fi

[[ ${#EVAL_ARGS[@]} -gt 0 ]] || die "every variant failed; nothing to compare (see ${OUTDIR}/*.log)"

# The reference findings for the expectation checks are the all-arms run when it is
# present, otherwise the first variant that succeeded.
REF="${OUTDIR}/all_arms/score_evidence/findings.tsv"
[[ -s "$REF" ]] || REF="$(echo "${EVAL_ARGS[1]}" | cut -d= -f2-)"

log "evaluating ..."
python3 "${REPO_DIR}/bin/evaluate_benchmark.py" \
  --scenarios "$SCENARIOS" \
  --manifest "$MANIFEST" \
  --watchlist "$WATCHLIST" \
  --findings "$REF" \
  "${EVAL_ARGS[@]}" \
  --out-dir "${OUTDIR}/evaluation" || log "evaluator reported failing expectations (exit non-zero)"

log "done"
log "  comparison : ${OUTDIR}/evaluation/bakeoff.tsv"
log "  per taxon  : ${OUTDIR}/evaluation/bakeoff_per_taxon.tsv"
log "  summary    : ${OUTDIR}/evaluation/summary.md"
log ""
log "Read bakeoff_per_taxon.tsv before ranking anything. An aggregate sensitivity"
log "figure hides the case where one arm finds the divergent filovirus and another"
log "finds the extra Plasmodium - that is a trade to decide, not a ranking."
