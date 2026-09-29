#!/usr/bin/env bash
#
# Put every conda spec declared in modules/local/*.nf through the solver.
#
# WHY THIS EXISTS
# The pins were written from memory and none had ever been exercised. On the first
# end-to-end run MARKER_EVIDENCE died with:
#
#   samtools =1.24 requires htslib >=1.24,<1.25.0a0, which requires
#   libdeflate >=1.25,<1.27.0a0, which does not exist (perhaps a missing channel)
#
# An unsatisfiable pin is not a subtle failure - it stops the pipeline dead - but it
# only shows up when that particular process first runs, which on a long pipeline can
# be an hour in and after the expensive steps have already completed. Solving all of
# them up front costs minutes.
#
# NOT part of tests/test_pipeline.py: that suite is deliberately offline and runs in
# seconds. This needs the network and takes minutes, so it is a separate, deliberate
# check to run after changing a pin.
#
# Usage: bash scripts/check_conda_specs.sh [--solver mamba|conda]
set -euo pipefail

_on_err() {
    local rc=$? line=${1:-?}
    printf 'ERROR: aborted at line %s with exit %s: %s\n' \
        "$line" "$rc" "${BASH_COMMAND:-unknown}" >&2
    exit "$rc"
}
trap '_on_err $LINENO' ERR

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOLVER="mamba"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --solver) SOLVER="${2:?--solver needs a value}"; shift 2 ;;
        -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
        *) printf 'unknown argument: %s\n' "$1" >&2; exit 2 ;;
    esac
done
command -v "$SOLVER" >/dev/null 2>&1 \
    || { printf 'ERROR: %s not on PATH\n' "$SOLVER" >&2; exit 2; }

mapfile -t SPECS < <(
    grep -h 'conda "' "${REPO_DIR}"/modules/local/*.nf \
        | sed 's/.*conda "//; s/"$//' | sort -u
)
(( ${#SPECS[@]} > 0 )) || { printf 'ERROR: no conda specs found\n' >&2; exit 2; }

printf 'solver: %s\nspecs : %d distinct\n\n' "$SOLVER" "${#SPECS[@]}"
fails=0
i=0
for spec in "${SPECS[@]}"; do
    i=$((i + 1))
    printf '%2d  %-52s ' "$i" "$spec"
    log="$(mktemp)"
    # --dry-run so nothing is written; we only want the solve.
    if timeout 420 "$SOLVER" create --dry-run --yes --quiet \
            -n "_btf_specprobe_${i}" -c conda-forge -c bioconda $spec \
            >"$log" 2>&1; then
        printf 'OK\n'
    else
        printf 'UNSATISFIABLE\n'
        sed -n 's/^ *//; /could not be installed\|does not exist\|Could not solve/p' \
            "$log" | head -4 | sed 's/^/        /'
        fails=$((fails + 1))
    fi
    rm -f "$log"
done

printf '\n'
if (( fails )); then
    printf '%d of %d spec(s) cannot be solved. The pipeline will die when the first\n' \
        "$fails" "${#SPECS[@]}" >&2
    printf 'affected process runs - possibly an hour in, after the expensive steps.\n' >&2
    exit 1
fi
printf 'all %d spec(s) solve\n' "${#SPECS[@]}"
