#!/usr/bin/env python3
"""check_host_recruitment.py - what fraction of SURVIVING host reads lands on a reference?

WHY THIS IS A SCRIPT

This number has now been published five times across D13/D14 and was wrong every time, in
the same direction. The last version read "residual human reads are at most 20,236, so the
in-pipeline rate is ~41%, about 7.6x the whole-genome coefficient". Three errors in one
sentence:

  1. THE BOUND WAS INVERTED. residual_host = total_screened - nonhost_survivors, and
     survivors are at most what was GENERATED, so the expression is a FLOOR on residual
     host, not a ceiling. Writing "at most" inverted the inequality.
  2. The bound was refuted by the same BAM it described: more host-origin reads map to the
     one reference than the bound allowed for the entire sample.
  3. The numerator counted supplementary alignments (-F 0x4) while the coefficient it was
     compared against excludes them (-F 0x904), so the ratio mixed units.

Every one was hand-derived. So the derivation lives here, the inequality direction is in
the code, and the units are forced to match.

THE DENOMINATOR IS EXACT, NOT A BOUND
An earlier version reported a FLOOR on residual host (screened minus non-host generated) and
therefore a CEILING on the rate, and its docstring asserted the exact count was unknowable
because "the surviving reads are not labelled". They ARE labelled - this script relies on
that label for its own numerator - and the BAM already opened holds EVERY screened read as a
primary record (primary count == total_reads_screened). So one more pass over the same file
gives the exact denominator. The floor is still computed and printed, as a cross-check: if
the exact count ever falls below it, something upstream is miscounted.

    python3 scripts/check_host_recruitment.py \\
        --batch $BTF_DB/simulated/lod_v1 \\
        --run   $BTF_DB/runs/lod_v2_kraken2 \\
        --sample negative_host_only_illumina_r1 --bundle 5811 \\
        --matrix $BTF_DB/benchmarks/recruitment_150bp_key.tsv \\
        --source chm13v2.0.fa.gz
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def die(msg: str) -> None:
    sys.exit(f"ERROR: {msg}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--batch", required=True, type=Path)
    ap.add_argument("--run", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--matrix", required=True, type=Path)
    ap.add_argument("--source", required=True, help="host row in the matrix")
    ap.add_argument("--host-prefix", default="NC_0609",
                    help="read-name prefix identifying host-origin simulated reads")
    ap.add_argument("--mapq", type=int, default=10)
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args()

    cov = a.run / "confirm" / f"{a.sample}.{a.bundle}.coverage.tsv"
    bam = a.run / "minimap2_competitive" / f"{a.sample}.{a.bundle}.bam"
    js = a.batch / "logs" / f"{a.sample}.json"
    for pth in (cov, bam, js):
        if not pth.exists():
            die(f"missing {pth}")

    obs = next(csv.DictReader(cov.open(), delimiter="\t"))
    screened = int(obs["total_reads_screened"])
    if int(obs["mapq_min"]) != a.mapq:
        die(f"coverage TSV says mapq_min={obs['mapq_min']} but --mapq is {a.mapq}; the "
            f"numerator and the pipeline must be filtered alike.")

    d = json.load(open(js))
    nonhost_pairs = sum((c.get("reads_generated") or c.get("reads_requested") or 0)
                        for c in d["components"] if c["kind"] != "host")
    nonhost_reads = nonhost_pairs * 2

    # FLOOR, not ceiling. Survivors of fastp and depletion are at most what was generated,
    # so subtracting the generated count can only UNDERSTATE residual host.
    residual_floor = screened - nonhost_reads
    if residual_floor <= 0:
        die(f"non-host generated ({nonhost_reads:,}) >= screened ({screened:,}); the "
            f"floor is not informative for this sample.")

    refs = subprocess.run(
        [a.samtools, "idxstats", str(bam)], capture_output=True, text=True, check=True
    ).stdout
    target = [ln.split("\t")[0] for ln in refs.splitlines()
              if ln.startswith(f"taxid|{a.bundle}|")]
    if not target:
        die(f"no taxid|{a.bundle}| references in {bam.name}")

    # -F 0x904 = primary, mapped, non-supplementary - the SAME filter
    # measure_recruitment.py uses, so numerator and coefficient are like for like.
    out = subprocess.run(
        [a.samtools, "view", "-q", str(a.mapq), "-F", "0x904", str(bam)] + target,
        capture_output=True, text=True, check=True).stdout
    recs = out.splitlines()
    host_recs = sum(1 for ln in recs if ln.startswith(a.host_prefix))

    # EXACT denominator: every screened read appears once as a primary record in this BAM.
    allrec = subprocess.run(
        [a.samtools, "view", "-F", "0x900", str(bam)],
        capture_output=True, text=True, check=True).stdout.splitlines()
    if len(allrec) != screened:
        die(f"BAM holds {len(allrec):,} primary records but total_reads_screened is "
            f"{screened:,}; they must match for the exact denominator to be valid.")
    residual_exact = sum(1 for ln in allrec if ln.startswith(a.host_prefix))
    # The guard belongs HERE, on the denominator, not on the numerator. Testing host_recs
    # rejected the common correct answer - 7 of this sample's 8 bundles legitimately recruit
    # ZERO host reads, which is the whole point of D13's table. What is never legitimate is a
    # prefix that matches nothing in the sample at all.
    if residual_exact == 0:
        die(f"--host-prefix {a.host_prefix!r} matches none of the {len(allrec):,} records in "
            f"this sample, so it cannot be the host naming for this batch (CHM13v2.0 is "
            f"NC_0609*; a GRCh38 analysis set is NC_0000* plus chrM NC_012920.1). A "
            f"zero-recruitment ANSWER is fine; a prefix that matches nothing is not.")
    if residual_exact < residual_floor:
        die(f"exact residual {residual_exact:,} is below the floor {residual_floor:,}, "
            f"which is impossible - one of the two is miscounted.")

    coef = None
    for r in csv.DictReader(a.matrix.open(), delimiter="\t"):
        if r["source"] == a.source and r["bundle_taxid"] == a.bundle \
                and r["ref_taxid"] == a.bundle:
            coef = r
    if coef is None:
        die(f"{a.matrix} has no row source={a.source} bundle=ref={a.bundle}")
    if int(coef["mapq_min"]) != a.mapq:
        die(f"matrix mapq_min={coef['mapq_min']} but numerator uses {a.mapq}")

    whole = float(coef["rpm_source_reads"]) / 1e4            # rpm -> percent
    ceiling = 100.0 * host_recs / residual_floor
    exact = 100.0 * host_recs / residual_exact

    print(f"sample                     {a.sample}")
    print(f"bundle                     {a.bundle}")
    print(f"screened records           {screened:,}")
    print(f"non-host generated         {nonhost_reads:,}  ({nonhost_pairs:,} pairs x 2)")
    print(f"residual host FLOOR        {residual_floor:,}   (cross-check only)")
    print(f"residual host EXACT        {residual_exact:,}   <- host-origin primary records "
          f"in the BAM")
    print(f"records on ref (MAPQ>={a.mapq}, primary)  {len(recs):,}")
    print(f"  of which host-origin     {host_recs:,}")
    print()
    print(f"whole-genome coefficient   {whole:.2f}%  (read_len={coef['read_len']}, "
          f"insert={coef.get('insert_mean','?')}, mapq>={coef['mapq_min']})")
    print(f"in-pipeline rate EXACT     {exact:.1f}%   = {host_recs:,} / {residual_exact:,}"
          f"   (ceiling against the floor would be {ceiling:.1f}%)")
    print(f"enrichment EXACT           {exact / whole:.1f}x")
    print()
    print("The mechanism: survivors of depletion are the reads that failed to map to the "
          "host reference, which is the same class that cross-maps to a divergent "
          "eukaryotic genome - so the whole-genome coefficient must not be inverted to "
          "estimate them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
