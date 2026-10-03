#!/usr/bin/env python3
"""check_pedestal_model.py - test `pedestal = recruitment_rpm x source_reads` against a run.

WHY THIS IS A SCRIPT AND NOT ARITHMETIC IN A COMMIT MESSAGE

The D12/D14 pedestal figure has been published four times and was wrong three times:
7.74% (records divided by pairs), 3.87% (right unit, measured at 101 bp for a 150 bp
batch), 5.51% (a 50k-read spot check quoted while a 500k run was producing a better
number), and a "within 6%" claim that was 6.31% on its own two inputs. Every one of those
was hand-derived from numbers copied between a terminal and a document.

So the comparison is computed here, from primary sources, and printed with the conditions
attached. Hand-copying a coefficient is the failure mode this file exists to retire.

WHAT IT READS
  recruitment matrix   scripts/measure_recruitment.py output (carries read_len and mapq)
  component JSON       <batch>/logs/<sample>.json - how many PAIRS of each organism were
                       simulated. Note `reads_requested`/`reads_actual` are PAIRS; the
                       pipeline's .nreads and total_reads_screened are RECORDS, i.e. 2x.
  coverage TSV         <run>/confirm/<sample>.<taxid>.coverage.tsv - the observed
                       best_neighbour_reads and unique_reads, in records.

WHAT IT REFUSES TO DO
  Invert the coefficient to estimate surviving HOST reads. Survivors of host depletion are
  not a random sample of host reads - they are the reads that failed to map to the host
  reference, which is the same class that cross-maps to a divergent eukaryotic genome. A
  coefficient measured on whole-genome simulated host reads is measured on a population
  the pipeline never sees, and inverting it produced an estimate larger than the sample's
  entire non-host remainder. The host arm is reported as an OBSERVED rate, never inverted.

    python3 scripts/check_pedestal_model.py \\
        --matrix /media/8tb/bioThreatFinder/benchmarks/recruitment_matrix_150bp.tsv \\
        --batch  /media/8tb/bioThreatFinder/simulated/lod_v1 \\
        --run    /media/8tb/bioThreatFinder/runs/lod_v2_kraken2 \\
        --source 562_Escherichia_coli.fna.gz --bundle 632 --neighbour 630
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
from pathlib import Path


def die(msg: str) -> None:
    sys.exit(f"ERROR: {msg}")


def load_matrix(path: Path, source: str, bundle: str) -> dict[str, dict]:
    rows = [r for r in csv.DictReader(path.open(), delimiter="\t")
            if r["source"] == source and r["bundle_taxid"] == bundle]
    if not rows:
        die(f"no rows in {path} for source={source} bundle={bundle}")
    if "read_len" not in rows[0]:  # noqa: SIM102 - message is the point
        die(f"{path} has no read_len column, so the conditions it was measured under are "
            f"unknown. Regenerate with scripts/measure_recruitment.py - applying a "
            f"coefficient at the wrong read length is the error this check exists to stop.")
    return {r["ref_taxid"]: r for r in rows}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True, type=Path)
    ap.add_argument("--batch", required=True, type=Path, help="simulated batch dir")
    ap.add_argument("--run", required=True, type=Path, help="pipeline run dir")
    ap.add_argument("--source", required=True, help="matrix source filename")
    ap.add_argument("--organism", default="Escherichia coli",
                    help="component label in the batch JSON for --source")
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--neighbour", required=True)
    ap.add_argument("--attributed", help="comma-separated observed records on the TARGET "
                    "reference attributable to --organism, from read attribution. Without "
                    "it the target arm is reported as not compared.")
    a = ap.parse_args()

    refs = load_matrix(a.matrix, a.source, a.bundle)
    for t in (a.bundle, a.neighbour):
        if t not in refs:
            die(f"matrix has no row for ref {t} in bundle {a.bundle}")
    cond = refs[a.bundle]
    print(f"coefficients from {a.matrix.name}")
    print(f"  measured at read_len={cond['read_len']} mapq>={cond['mapq_min']} "
          f"insert_mean={cond.get('insert_mean', '?')} "
          f"over {int(cond['source_reads']):,} source reads")

    # Batch read length must match the coefficient's, or the comparison is meaningless.
    cmds = sorted(glob.glob(str(a.batch / "logs" / "*.cmds")))
    lens = set()
    for c in cmds:
        for line in open(c):
            if " -l " in line:
                parts = line.split()
                lens.add(parts[parts.index("-l") + 1])
    if not lens:
        die(f"no ART commands found under {a.batch / 'logs'}/*.cmds, so the batch's read "
            f"length is unknown. Previously this printed 'conditions match' and continued, "
            f"which is worse than no guard: it asserts a condition it could not check.")
    if {str(cond["read_len"])} != lens:
        die(f"coefficient measured at {cond['read_len']} bp but the batch was simulated at "
            f"{sorted(lens)} bp. Recruitment is read-length dependent; this comparison "
            f"would be invalid.")
    # Insert size is checked too: the two 150 bp matrices on disk differ ONLY by insert
    # (300 vs 350) and that moves the target-arm coefficient by 10% and ref 633 by 28%.
    inserts = set()
    for c in cmds:
        for line in open(c):
            parts = line.split()
            if "-m" in parts:
                inserts.add(parts[parts.index("-m") + 1])
    if inserts and str(cond.get("insert_mean", "")) not in inserts:
        die(f"coefficient measured at insert {cond.get('insert_mean')} but the batch used "
            f"{sorted(inserts)}. That alone moves coefficients by 10-28%.")
    print(f"  batch simulated at {sorted(lens)} bp, insert {sorted(inserts) or ['?']} "
          f"-> conditions match")

    rows = []
    for sample_json in sorted(glob.glob(str(a.batch / "logs" / "*.json"))):
        d = json.load(open(sample_json))
        sample = d["sample"]
        cov = a.run / "confirm" / f"{sample}.{a.bundle}.coverage.tsv"
        if not cov.exists():
            continue
        # reads_generated is what ART actually produced; reads_requested is the plan.
        # There is no reads_actual - an earlier version looked for it and so always
        # silently used the planned figure.
        pairs = sum((c.get("reads_generated") or c.get("reads_requested") or 0)
                    for c in d["components"] if c.get("label") == a.organism)
        if not pairs:
            continue
        src_reads = pairs * 2          # components are PAIRS; the pipeline counts RECORDS
        obs = next(csv.DictReader(cov.open(), delimiter="\t"))
        if obs.get("best_neighbour_taxid") != a.neighbour:
            # Comparing a coefficient for ref X against a count for ref Y is silent
            # nonsense, so it is refused rather than averaged in.
            die(f"{sample}: best_neighbour_taxid is {obs.get('best_neighbour_taxid')}, "
                f"not the --neighbour {a.neighbour} the coefficient was taken for.")
        rows.append((sample, src_reads,
                     float(refs[a.neighbour]["rpm_source_reads"]) * src_reads / 1e6,
                     int(obs["best_neighbour_reads"]),
                     float(refs[a.bundle]["rpm_source_reads"]) * src_reads / 1e6))
    skipped = len(glob.glob(str(a.batch / "logs" / "*.json"))) - len(rows)
    if not rows:
        die(f"no sample had both a component JSON and {a.bundle}.coverage.tsv")

    print(f"\n{'sample':<42} {'src reads':>10} {'pred nb':>8} {'obs nb':>7} {'err':>7}")
    errs = []
    for s, src, pred_nb, obs_nb, _ in rows:
        e = pred_nb / obs_nb - 1
        errs.append(e)
        print(f"{s:<42} {src:>10,} {pred_nb:>8.0f} {obs_nb:>7,} {e:>+6.1%}")
    mean = sum(errs) / len(errs)
    print(f"\nNEIGHBOUR ARM: mean error {mean:+.1%} over {len(errs)} sample(s)"
          f"  ({skipped} sample(s) had no {a.bundle} coverage TSV and were skipped - "
          f"stage 1 never nominated the taxon there)")

    # The target arm is reported too. Quoting only the arm that agrees is how a model gets
    # called calibrated when half of it is not.
    tgt_pred = sum(r[4] for r in rows) / len(rows)
    print(f"TARGET ARM:    predicted {tgt_pred:.0f} records on ref {a.bundle} per sample.")
    if a.attributed:
        obs_t = [float(x) for x in a.attributed.split(",")]
        m = sum(obs_t) / len(obs_t)
        print(f"               observed {min(obs_t):.0f}-{max(obs_t):.0f} "
              f"(mean {m:.0f}, n={len(obs_t)}, supplied via --attributed) "
              f"-> {tgt_pred / m - 1:+.1%}")
    else:
        print("               NOT COMPARED: the coverage TSV does not separate which of "
              "the target's reads came from the background organism, so the observed "
              "value must come from read attribution and be passed with --attributed. "
              "Do not quote a target-arm error without it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
