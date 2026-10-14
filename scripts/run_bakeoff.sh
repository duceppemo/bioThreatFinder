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
# NOTE ON THE SKIP FLAGS: EVERY VARIANT SPELLS OUT EVERY ARM.
# Two separate failures here, both of which made a variant silently equal the baseline:
#
#  1. DEFAULTS DRIFT. The matrix once read `plus_krakenuniq|--skip_metabuli`, relying on
#     skip_krakenuniq defaulting false. It defaults TRUE, so that variant was a
#     byte-for-byte duplicate of the baseline and its null result read as "KrakenUniq adds
#     nothing". skip_metabuli has since ALSO flipped to true (D15), which would have done
#     the same thing to plus_metabuli and all_arms.
#  2. `--flag false` WAS A NO-OP. The repair for (1) was `--skip_krakenuniq false`, and on
#     the Nextflow CLI that is the string "false", which is truthy - so the repair did not
#     work either and D2 stayed unfalsifiable by its own named test. Fixed in
#     lib/BtfBool.groovy; `--skip_X false` now genuinely enables the arm.
#
# So: no variant may rely on a default, and none may omit an arm it has an opinion about.
# If you add a boolean arm, add it to EVERY row here.
VARIANTS=(
  "kraken2_minimizer|--skip_krakenuniq true --skip_metabuli true|Baseline. Kraken2 with --report-minimizer-data (distinct minimizers + duplicity) plus Bracken and sylph. The cheapest configuration that still has distinct-k-mer-style evidence."
  "plus_krakenuniq|--skip_krakenuniq false --skip_metabuli true|Baseline + KrakenUniq. ISOLATES the question: does HyperLogLog cardinality beat Kraken2 minimizer counts by enough to justify a second index? Needs a built KrakenUniq index. This is the variant decisions.md D2 names as what would reverse the Kraken2-primary decision."
  "plus_metabuli|--skip_krakenuniq true --skip_metabuli false|Baseline + Metabuli. Isolates the value of the DNA+amino-acid arm. Run on divergent_v1 (D15) it changed no decision; its untested case is an agent absent from the index."
  "all_arms|--skip_krakenuniq false --skip_metabuli false|Everything on, KrakenUniq included. Upper bound on sensitivity and on cost."
  "no_sylph|--skip_krakenuniq true --skip_metabuli true --skip_sylph true|Baseline minus sylph, to check whether ANI/containment is contributing anything the k-mer arms do not."
  "no_markers|--skip_krakenuniq true --skip_metabuli false --skip_markers true|All arms except KrakenUniq, marker gating OFF. Quantifies what marker gating buys - expect confirm PPV to fall and false confirmations to rise."
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
[[ -s "$WATCHLIST" ]]   || die "no resolved watchlist at ${WATCHLIST}. Run bin/resolve_watchlist.py first."
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
  # THE SKIP MUST CHECK WHAT PRODUCED THE FILE, NOT JUST THAT IT EXISTS.
  # Existence alone reuses a table built under different flags, and that is not
  # hypothetical here: before lib/BtfBool.groovy, `--skip_metabuli false` was a no-op, so
  # plus_metabuli/findings.tsv was produced with the arm OFF and is byte-identical to the
  # baseline. Re-running after the fix would print SKIP, feed the stale table to the
  # comparison, and reproduce the very conclusion the fix exists to re-test - D15 for
  # Metabuli, D2 for KrakenUniq. The args are stamped beside the findings and compared.
  argstamp="${vdir}/.variant_args"
  # THE FLAGS ALONE ARE NOT WHAT MAKES A RESULT REUSABLE. The motivating case for this
  # stamp was a CODE fix (--skip_metabuli false being a no-op), after which $extra is
  # unchanged, the stamp matches, SKIP fires and the stale table is fed to the comparison
  # anyway - the bug intact. Every fix since has the same property: a breadth metric that
  # was 1000x wrong does not change a variant's flags. The pack and profile matter too,
  # because `--pack animal` over the same --batch resolves to the same OUTDIR.
  stamp_now="$(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || echo nogit)|${PACK}|${PROFILE}|${extra}"
  stamped=""
  [[ -r "$argstamp" ]] && stamped="$(cat "$argstamp")"
  if [[ -s "$findings" ]] && (( $(grep -vc '^$' "$findings" 2>/dev/null || echo 0) > 1 )) \
     && [[ "$stamped" == "$stamp_now" ]]; then
    log "SKIP ${name} (findings already present, same flags)"
  else
    if [[ -s "$findings" ]] && [[ "$stamped" != "$stamp_now" ]]; then
      log "RERUN ${name}: findings exist but were built as ${stamped:-<unrecorded>},"
      log "      not ${stamp_now}. Reusing them would answer a different question."
    fi
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
      # The stamp is written only AFTER the run succeeds. Writing it first looks
      # harmless and is not: nextflow fails (OOM, ctrl-C, a starved arm), `continue`
      # runs, and the PREVIOUS findings.tsv survives because nextflow does not delete
      # what it already published. The next invocation then sees findings + a matching
      # stamp and skips - feeding a table built under the old flags to a comparison
      # labelled with the new ones, which is the bug this stamp exists to prevent.
      mkdir -p "$vdir"
      rm -f "$argstamp"
      if ! "${cmd[@]}" >"${OUTDIR}/${name}.log" 2>&1; then
        log "WARNING: variant ${name} FAILED - see ${OUTDIR}/${name}.log"
        log "         it will be omitted from the comparison, and its flag stamp is NOT"
        log "         written, so a later run re-executes it rather than reusing whatever"
        log "         findings.tsv survived this failure"
        continue
      fi
      printf '%s' "$stamp_now" > "$argstamp"
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

# The descendants file is BUILT BY THE RUN, not shipped in the pack: it is derived
# from whichever taxdump the pipeline was pointed at, and a stale copy would join a
# group lead to the wrong species set. Taken from the first arm that published one;
# every arm resolves the same watchlist against the same taxdump, so they agree.
# Absent (an older run, or a pack with no group-rank entry) the evaluator simply
# scores no leads, which is the behaviour before this existed.
DESCENDANTS=""
for _arm in "${OUTDIR}"/*/watchlist/watchlist.descendants.tsv; do
  [[ -s "$_arm" ]] && { DESCENDANTS="$_arm"; break; }
done
DESC_ARGS=()
if [[ -n "$DESCENDANTS" ]]; then
  DESC_ARGS=(--descendants "$DESCENDANTS")
  log "descendants  : ${DESCENDANTS}"
else
  log "no watchlist.descendants.tsv published; group-rank leads will not be scored"
fi

log "evaluating ..."
python3 "${REPO_DIR}/bin/evaluate_benchmark.py" \
  --scenarios "$SCENARIOS" \
  --manifest "$MANIFEST" \
  --watchlist "$WATCHLIST" \
  "${DESC_ARGS[@]}" \
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
