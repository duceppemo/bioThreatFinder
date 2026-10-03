#!/usr/bin/env python3
"""measure_recruitment.py - how many reads from organism X land on reference Y?

WHAT THIS MEASURES AND WHY IT IS A DIFFERENT QUESTION FROM DETECTION

D12 and D13 both rest on one observation made in one simulated batch: a reference can
recruit reads from an organism that is NOT it, in quantity, at MAPQ>=10. Two numbers were
derived from that single batch and then used as if they were constants -

  * ~2,700 E. coli reads land on the Y. enterocolitica reference while only ~210 land on
    the Y. pestis reference, a ~13x asymmetry that biases neighbour_ratio (D12);
  * the T. gondii reference absorbs thousands of HUMAN reads, three orders of magnitude
    more than any bacterial reference (D13).

Both are properties of REFERENCE CONTENT, not of error profiles or of any specimen. So
they can be measured directly, per (source organism -> reference) pair, without a benchmark
run, without ground-truth bookkeeping and without downloading anything: simulate reads from
a single organism, align them against a competitive bundle, and count what sticks.

The output is a recruitment matrix in reads per million SOURCE reads, which is the unit
that makes bundles comparable. A pedestal predicted from it is then
`rpm_recruitment * background_abundance`, which is the quantity D12 needs and currently
guesses.

WHY SINGLE-ORGANISM INPUT IS THE WHOLE POINT
Every read is known to be foreign by construction, so no read-name trick, no truth table
and no classifier is involved. A read on the reference is recruitment, full stop. This is
also why real mock-community data cannot answer it while a single-organism library can.

WHAT IT DOES NOT MEASURE
Whether recruitment causes a false POSITIVE. Breadth, occupancy, identity and depth_cv all
still apply downstream, and in the T. gondii case they are what hold the line. Recruitment
is the input to those gates, not a verdict.

Reads are simulated with ART at the read length that matters for the pack, because mapping
is length-dependent: a 101 bp read maps on shorter and less identical stretches than the
200 bp / 95% identity alignments mask_pack_refs.py looks for, which is exactly the scale
mismatch D13 records.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

TAXID_RE = re.compile(r"^taxid\|(\d+)\|")


def log(msg: str) -> None:
    print(f"[recruit] {msg}", file=sys.stderr, flush=True)


def die(msg: str) -> None:
    print(f"[recruit] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


def need(tool: str) -> str:
    p = shutil.which(tool)
    if not p:
        die(f"{tool} not found on PATH. art_illumina lives in the btf-sim env and "
            f"minimap2/samtools in btf-refs, so one PATH rarely has all three - "
            f"prepend both, or pass --art/--minimap2/--samtools.")
    return p


def genome_bp(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    n = 0
    with opener(path, "rt") as fh:                             # type: ignore[operator]
        for line in fh:
            if not line.startswith(">"):
                n += len(line.strip())
    return n


def simulate(src: Path, n_reads: int, read_len: int, art: str, tmp: Path) -> Path:
    """ART reads from one genome. Coverage is derived from the requested read count so
    the number is comparable across genomes of wildly different size."""
    bp = genome_bp(src)
    if bp == 0:
        die(f"{src} contains no sequence")
    # Paired 2x{read_len}; fold coverage that yields roughly n_reads pairs.
    fold = (n_reads * read_len * 2) / bp
    plain = tmp / (src.name.replace(".gz", ""))
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as fi, plain.open("wb") as fo:
            shutil.copyfileobj(fi, fo)
    else:
        plain = src
    out = tmp / "sim"
    log(f"  simulating ~{n_reads:,} pairs from {src.name} "
        f"({bp/1e6:.1f} Mbp, {fold:.4f}x)")
    cmd = [art, "-ss", "HS25", "-i", str(plain), "-p", "-l", str(read_len),
           "-f", f"{fold:.6f}", "-m", "300", "-s", "50", "-o", str(out),
           "-na", "-rs", "42"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        die(f"art_illumina failed on {src.name}: {r.stderr.strip()[-400:]}")
    r1, r2 = tmp / "sim1.fq", tmp / "sim2.fq"
    if not r1.exists() or not r2.exists():
        die(f"art_illumina produced no reads for {src.name} at {fold:.6f}x - the genome "
            f"may be too small for the requested count.")
    got = sum(1 for _ in r1.open()) // 4
    if got == 0:
        die(f"0 reads simulated from {src.name}; refusing to report 0 recruitment, which "
            f"would be indistinguishable from a reference that recruits nothing.")
    log(f"  simulated {got:,} pairs")
    return r1


def recruit(r1: Path, r2: Path, index: Path, mapq: int,
            minimap2: str, samtools: str, threads: int) -> dict[str, int]:
    """MAPQ-passing read count per reference sequence in one bundle."""
    mm = subprocess.Popen([minimap2, "-ax", "sr", "-t", str(threads), "--secondary=no",
                           str(index), str(r1), str(r2)],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    st = subprocess.Popen([samtools, "view", "-q", str(mapq), "-F", "0x904", "-"],
                          stdin=mm.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, text=True)
    if mm.stdout:
        mm.stdout.close()
    per_ref: dict[str, int] = defaultdict(int)
    assert st.stdout is not None
    for line in st.stdout:
        f = line.split("\t", 3)
        if len(f) > 2:
            per_ref[f[2]] += 1
    st.wait()
    mm.wait()
    if mm.returncode not in (0, None):
        die(f"minimap2 exited {mm.returncode} for {index.name}")
    return dict(per_ref)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", action="append", required=True, type=Path,
                    help="single-organism genome FASTA; repeatable")
    ap.add_argument("--bundles", required=True, type=Path,
                    help="dir of competitive bundles (<taxid>.sr.mmi)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--only-bundle", action="append", default=[],
                    help="bundle taxid; repeatable. Default: all.")
    ap.add_argument("--reads", type=int, default=200_000, help="read pairs per source")
    ap.add_argument("--read-len", type=int, default=101,
                    help="must match the batch being explained; 101 for lod_v1")
    ap.add_argument("--mapq", type=int, default=10, help="match pack.yaml mapq_min")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--art", default="art_illumina")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    args = ap.parse_args()

    art, mm2, st = need(args.art), need(args.minimap2), need(args.samtools)
    for s in args.source:
        if not s.exists():
            die(f"source genome not found: {s}")

    idx = sorted(args.bundles.glob("*.sr.mmi"))
    if args.only_bundle:
        want = set(args.only_bundle)
        idx = [i for i in idx if i.name.split(".")[0] in want]
        missing = want - {i.name.split(".")[0] for i in idx}
        if missing:
            die(f"no bundle index for taxid(s): {', '.join(sorted(missing))}")
    if not idx:
        die(f"no *.sr.mmi bundles under {args.bundles}")
    log(f"{len(args.source)} source(s) x {len(idx)} bundle(s), "
        f"{args.reads:,} pairs each, MAPQ>={args.mapq}, {args.read_len} bp")

    rows: list[dict] = []
    for src in args.source:
        with tempfile.TemporaryDirectory(prefix="recruit.") as td:
            tmp = Path(td)
            r1 = simulate(src, args.reads, args.read_len, art, tmp)
            r2 = r1.with_name("sim2.fq")
            n_pairs = sum(1 for _ in r1.open()) // 4
            for i, index in enumerate(idx, 1):
                bundle = index.name.split(".")[0]
                per_ref = recruit(r1, r2, index, args.mapq, mm2, st, args.threads)
                by_taxid: dict[str, int] = defaultdict(int)
                for ref, n in per_ref.items():
                    m = TAXID_RE.match(ref)
                    by_taxid[m.group(1) if m else ref] += n
                total = sum(by_taxid.values())
                for ref_taxid, n in sorted(by_taxid.items(), key=lambda kv: -kv[1]):
                    rows.append({
                        "source": src.name,
                        "bundle_taxid": bundle,
                        "ref_taxid": ref_taxid,
                        "is_bundle_target": "yes" if ref_taxid == bundle else "no",
                        "reads": n,
                        "source_pairs": n_pairs,
                        # per million SOURCE reads - the unit that makes bundles comparable
                        # and that multiplies by background abundance to predict a pedestal
                        "rpm_source": f"{n / n_pairs * 1e6:.1f}",
                    })
                if i % 10 == 0 or total:
                    log(f"  [{i}/{len(idx)}] bundle {bundle}: {total:,} read(s) recruited")

    if not rows:
        die("no recruitment rows produced; refusing to write an empty matrix")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {args.out} ({len(rows)} row(s))")

    # Loudest line last: the references that absorb the most foreign reads.
    worst = sorted(rows, key=lambda r: -float(r["rpm_source"]))[:10]
    log("highest recruitment (reads per million source reads):")
    for r in worst:
        log(f"  {r['source']:<34} -> bundle {r['bundle_taxid']:<7} "
            f"ref {r['ref_taxid']:<7} {r['rpm_source']:>10} rpm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
