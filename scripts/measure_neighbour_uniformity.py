#!/usr/bin/env python3
"""measure_neighbour_uniformity.py - is a competitive-bundle reference carrying a REAL
organism, or a cross-mapping pedestal?

WHY THIS IS A SCRIPT

D12 records that `neighbour_ratio`'s denominator, for the *Y. pestis* bundle in
`lod_gram_negative`, is ~98% *E. coli* background and 0% *Yersinia*: there is no neighbour
in the specimen at all, so the "competition" the gate represents does not exist. Remedy
option 2 there is "require the neighbour to be credible before it can out-vote - gate
attribution on the neighbour's own breadth and occupancy, not its raw count", with the
stated risk that it "adds a second threshold to calibrate".

This measures the quantity that threshold would be set on, so the choice rests on data
rather than on the plausibility of the sentence. The coverage TSV the pipeline already
writes carries `best_neighbour_reads` and `neighbour_ratio` but NO breadth for the
neighbour, so the number cannot be read off an existing artifact - it has to come from the
BAM.

THE METRIC

    uniformity = observed_breadth / breadth_expected_if_reads_were_spread_uniformly
               = (covered_bp / reference_bp) / min(1, reads * read_len / reference_bp)

A genuinely present organism's reads land across its genome, so observed breadth
approaches the uniform expectation and uniformity approaches 1. Reads that cross-map from
another organism pile onto the conserved fraction the two genomes share, covering far less
than their bulk implies, so uniformity collapses. It is a SHAPE statistic: deliberately
near-independent of abundance, which is what lets one threshold serve every level.

WHY NOT RAW BREADTH. Breadth alone cannot separate the two - it scales with abundance, so a
trace true positive and a fat pedestal can have the same breadth. At 100 rpm in lod_v4 the
target holds 2,336 reads and the absent neighbour 2,731, counts close enough that the raw
ratio attributes the sample to the neighbour; their uniformities are 0.86 and 0.12.

    python3 scripts/measure_neighbour_uniformity.py \\
        --run /media/8tb/bioThreatFinder/runs/lod_v4 --bundle 632 --read-len 150

Read length is REQUIRED and not defaulted: the expectation is linear in it, so a wrong
value rescales every uniformity reported. D13 and D14 both record figures that were wrong
because a read length was assumed rather than taken from the batch.
"""
from __future__ import annotations

import argparse
import collections
import glob
import os
import re
import subprocess
import sys
from pathlib import Path


def die(msg: str) -> None:
    sys.exit(f"ERROR: {msg}")


def measure(bam: Path, samtools: str, mapq: int) -> dict:
    """Per-taxid covered bases, reference length and read count from one bundle BAM.

    `--ff UNMAP,SECONDARY,QCFAIL,DUP` plus `-q` is the same population the pipeline's own
    gates are computed over: primary, mapped, MAPQ-filtered. Mixing that with an -F 0x4
    count is the unit error D13/D14 record three times over.
    """
    proc = subprocess.run(
        [samtools, "coverage", "-q", str(mapq),
         "--ff", "UNMAP,SECONDARY,QCFAIL,DUP", str(bam)],
        capture_output=True, text=True)
    if proc.returncode != 0:
        die(f"samtools coverage failed on {bam.name}: {proc.stderr.strip()[:300]}")
    agg: dict = collections.defaultdict(lambda: {"bp": 0, "cov": 0, "reads": 0})
    for line in proc.stdout.splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 7:
            continue
        parts = f[0].split("|")
        if len(parts) < 2:
            die(f"{bam.name}: reference {f[0]!r} is not taxid|<taxid>|<accession>; this "
                f"script reads the taxid out of the reference name and cannot proceed.")
        a = agg[parts[1]]
        a["bp"] += int(f[2])
        a["reads"] += int(f[3])
        a["cov"] += int(f[4])
    if not agg:
        die(f"{bam.name}: no references passed the filters; nothing to measure.")
    return agg


def uniformity(d: dict, read_len: int) -> float:
    if d["bp"] == 0 or d["reads"] == 0:
        return float("nan")
    expected = min(1.0, d["reads"] * read_len / d["bp"])
    if expected <= 0:
        return float("nan")
    return (d["cov"] / d["bp"]) / expected


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, type=Path, help="pipeline outdir")
    ap.add_argument("--bundle", required=True,
                    help="bundle taxid, i.e. the TARGET whose competitive BAMs to read")
    ap.add_argument("--read-len", required=True, type=int,
                    help="REQUIRED. The batch's read length; the expectation is linear in it.")
    ap.add_argument("--mapq", type=int, default=10)
    ap.add_argument("--min-neighbour-reads", type=int, default=50,
                    help="ignore neighbours below this; a handful of reads has no shape")
    ap.add_argument("--pedestal-max", type=float, default=0.30,
                    help="at or below this, a neighbour's reads are pile-up, not coverage")
    ap.add_argument("--present-min", type=float, default=0.60,
                    help="at or above this, a neighbour is genuinely in the specimen")
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args()

    if a.read_len <= 0:
        die("--read-len must be positive")

    bams = sorted(glob.glob(str(a.run / "minimap2_competitive" / f"*.{a.bundle}.bam")))
    if not bams:
        die(f"no *.{a.bundle}.bam under {a.run}/minimap2_competitive")

    print(f"# run={a.run}  bundle={a.bundle}  read_len={a.read_len}  mapq>={a.mapq}")
    print(f"# uniformity = observed_breadth / uniform-expectation; ~1 present, <<1 pedestal")
    print(f"{'sample':<44}{'tgt_reads':>11}{'tgt_unif':>10}"
          f"{'nb_taxid':>10}{'nb_reads':>10}{'nb_unif':>9}{'unif_ratio':>11}")
    rows = []
    for b in bams:
        bam = Path(b)
        samp = bam.name[: -len(f".{a.bundle}.bam")]
        agg = measure(bam, a.samtools, a.mapq)
        tgt = agg.get(a.bundle)
        if tgt is None:
            die(f"{bam.name}: bundle taxid {a.bundle} has no reference in its own BAM.")
        nbs = {k: v for k, v in agg.items()
               if k != a.bundle and v["reads"] >= a.min_neighbour_reads}
        if not nbs:
            print(f"{samp:<44}{tgt['reads']:>11,}{uniformity(tgt, a.read_len):>10.3f}"
                  f"{'-':>10}{0:>10}{'-':>9}{'-':>11}")
            continue
        nk, nb = max(nbs.items(), key=lambda kv: kv[1]["reads"])
        ut, un = uniformity(tgt, a.read_len), uniformity(nb, a.read_len)
        ratio = ut / un if un and un == un and un > 0 else float("nan")
        rows.append((ut, un))
        print(f"{samp:<44}{tgt['reads']:>11,}{ut:>10.3f}"
              f"{nk:>10}{nb['reads']:>10,}{un:>9.3f}{ratio:>11.2f}")

    if rows:
        ts = [u for u, _ in rows]
        ns = [n for _, n in rows]
        print()
        print(f"target    uniformity: min {min(ts):.3f}  max {max(ts):.3f}")
        print(f"neighbour uniformity: min {min(ns):.3f}  max {max(ns):.3f}")
        print()
        # THE QUESTION IS ABOUT THE NEIGHBOUR ALONE, NOT ABOUT A COMPARISON.
        # An earlier version of this summary asked "are the neighbours below the targets?",
        # which silently assumed the defective case. On a bundle whose declared neighbour IS
        # in the specimen the ordering simply inverts - the neighbour is the credible one and
        # the unspiked target is the pedestal - and that version printed NOT SEPARATED for
        # three bundles that are in fact cleanly separated the other way round. The gate
        # being considered (D12 option 2) tests whether the NEIGHBOUR deserves to out-vote,
        # so the neighbour's own shape is the whole input; the target's is printed for
        # context only.
        lo, hi = min(ns), max(ns)
        if hi < a.pedestal_max:
            print(f"NEIGHBOUR IS PEDESTAL-LIKE on every sample ({lo:.3f}-{hi:.3f} < "
                  f"{a.pedestal_max:.2f}): its reads pile on shared sequence rather than "
                  f"covering its genome, so on this bundle it is not evidence of a "
                  f"competing organism and D12's defect applies.")
        elif lo > a.present_min:
            print(f"NEIGHBOUR IS PRESENT-LIKE on every sample ({lo:.3f}-{hi:.3f} > "
                  f"{a.present_min:.2f}): it really is in these specimens, so the "
                  f"attribution gate is doing the job it claims here and must keep "
                  f"out-voting.")
        else:
            print(f"MIXED on this bundle ({lo:.3f}-{hi:.3f} spans "
                  f"{a.pedestal_max:.2f}-{a.present_min:.2f}): some samples' neighbours look "
                  f"credible and others look like pedestals, so no single threshold "
                  f"classifies this bundle. Report it rather than tuning until it resolves.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
