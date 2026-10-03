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

WHAT IT CANNOT TELL YOU
The exact residual host count. The surviving reads are not labelled, and the classifier
undercounts them badly - host reads that survive depletion are the ones that failed to map
to the host reference, which is the same class that evades exact k-mer matching. So this
reports a FLOOR on residual host and therefore a CEILING on the recruitment rate, and says
so.

    python3 scripts/check_host_recruitment.py \\
        --batch /media/8tb/bioThreatFinder/simulated/lod_v1 \\
        --run   /media/8tb/bioThreatFinder/runs/lod_v2_kraken2 \\
        --sample negative_host_only_illumina_r1 --bundle 5811 \\
        --matrix /media/8tb/bioThreatFinder/benchmarks/recruitment_150bp_key.tsv \\
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

    print(f"sample                     {a.sample}")
    print(f"bundle                     {a.bundle}")
    print(f"screened records           {screened:,}")
    print(f"non-host generated         {nonhost_reads:,}  ({nonhost_pairs:,} pairs x 2)")
    print(f"residual host FLOOR        {residual_floor:,}   <- a LOWER bound; fastp and "
          f"depletion losses on the non-host side make the true value larger")
    print(f"records on ref (MAPQ>={a.mapq}, primary)  {len(recs):,}")
    print(f"  of which host-origin     {host_recs:,}")
    print()
    print(f"whole-genome coefficient   {whole:.2f}%  (read_len={coef['read_len']}, "
          f"insert={coef.get('insert_mean','?')}, mapq>={coef['mapq_min']})")
    print(f"in-pipeline rate CEILING   {ceiling:.1f}%   = {host_recs:,} / {residual_floor:,}")
    print(f"enrichment CEILING         {ceiling / whole:.1f}x")
    print()
    print("Both in-pipeline figures are CEILINGS because the denominator is a floor. The "
          "qualitative finding - survivors of depletion recruit far above the whole-genome "
          "rate - does not depend on the exact value; the mechanism is that survivors are "
          "the reads that failed to map to the host reference, which is the same class "
          "that cross-maps to a divergent eukaryotic genome.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
