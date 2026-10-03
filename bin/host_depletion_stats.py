#!/usr/bin/env python3
"""host_depletion_stats.py - turn samtools flagstat into a host-depletion figure.

WHY THIS IS A FIRST-CLASS NUMBER AND NOT A LINE IN A LOG

Every eukaryotic reference in the pack recruits host reads in proportion to how many
survive host depletion. Measured 2026-10-03 (`scripts/measure_recruitment.py`), per
million human reads reaching the competitive alignment:

    Toxoplasma gondii      77,354 rpm   (7.74% of human reads)
    Plasmodium knowlesi    44,340 rpm
    Entamoeba dispar        8,160 rpm

In a negative control with no parasite present, that produced 8,380 reads on the
T. gondii reference - 100% of them human - at 44% breadth and 96.5% window occupancy,
reported as tier `supported`, an assertion of presence for a notifiable agent. Only
`depth_cv` stopped it being `confirmed_in_silico`. See docs/decisions.md D14.

Recruitment scales LINEARLY with whatever survives depletion, so this fraction is the
real control on that entire class of false positive. It was already being computed by
`samtools flagstat` inside MINIMAP2_HOST_DEPLETE, published, and then read by nothing:
`ch_host_stats` was a dead channel.

WHAT THE DENOMINATOR IS, SAID PLAINLY
flagstat counts alignment RECORDS, not fragments: a pair contributes two records plus one
per supplementary alignment. So `host_fraction` here is the fraction of alignment records
that matched the host, which is the right quantity for "how much host sequence got
through" but is NOT the same as the fraction of input pairs removed. The column is named
for records so the two cannot be confused.

A missing or unparseable flagstat is an ERROR, not a zero. "0% host" and "could not tell"
look identical downstream, and the whole point of this file is that the number is
load-bearing.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

TOTAL_RE = re.compile(r"^(\d+)\s*\+\s*(\d+)\s+in total")
# "113 + 0 mapped (0.03% : N/A)" - and NOT "primary mapped", which is a subset.
MAPPED_RE = re.compile(r"^(\d+)\s*\+\s*(\d+)\s+mapped\s*\(")

COLS = ("sample", "platform", "alignment_records", "host_matching_records",
        "host_fraction_of_records", "surviving_records")


def parse(path: Path) -> tuple[int, int]:
    total = mapped = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if total is None:
            m = TOTAL_RE.match(line)
            if m:
                total = int(m.group(1)) + int(m.group(2))
                continue
        if mapped is None and " mapped (" in line and "primary" not in line:
            m = MAPPED_RE.match(line)
            if m:
                mapped = int(m.group(1)) + int(m.group(2))
    if total is None:
        sys.exit(f"ERROR: no 'in total' line in {path}. Refusing to report 0% host, "
                 f"which is indistinguishable from a depletion step that did nothing.")
    if mapped is None:
        sys.exit(f"ERROR: no 'mapped' line in {path}. Same reason: silence here would "
                 f"read as a clean specimen.")
    if mapped > total:
        sys.exit(f"ERROR: {path} reports {mapped} host-matching of {total} total records.")
    return total, mapped


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flagstat", type=Path)
    ap.add_argument("--sample")
    ap.add_argument("--platform", default="")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if not args.self_test:
        missing = [n for n, v in (("--flagstat", args.flagstat), ("--sample", args.sample),
                                  ("--out", args.out)) if not v]
        if missing:
            ap.error(f"the following arguments are required: {', '.join(missing)}")

    if args.self_test:
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            f = Path(td) / "fs.txt"
            f.write_text("399275 + 0 in total (QC-passed reads + QC-failed reads)\n"
                         "399246 + 0 primary\n"
                         "113 + 0 mapped (0.03% : N/A)\n"
                         "84 + 0 primary mapped (0.02% : N/A)\n")
            t, m = parse(f)
            assert (t, m) == (399275, 113), (t, m)
            # "primary mapped" must NOT be taken for "mapped"
            assert m != 84
        print("host_depletion_stats self-test OK", file=sys.stderr)
        return 0

    total, mapped = parse(args.flagstat)
    frac = mapped / total if total else 0.0
    args.out.write_text(
        "\t".join(COLS) + "\n"
        + f"{args.sample}\t{args.platform}\t{total}\t{mapped}\t{frac:.6f}\t{total - mapped}\n")
    print(f"[host] {args.sample}: {mapped}/{total} records host-matching "
          f"({frac * 100:.2f}%), {total - mapped} surviving", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
