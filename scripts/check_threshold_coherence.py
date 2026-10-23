#!/usr/bin/env python3
"""Report which of a pack's two stage-2 read/breadth gates actually binds, per taxon.

WHY THIS EXISTS

`min_unique_reads` and `min_breadth` are ANDed in bin/score_evidence.py, and the
watchlist sets both per taxon. But breadth is bounded by read count: with reads of
length L on a genome of size G, the best breadth that N reads can possibly reach is
N*L/G, and that is only achieved if no two reads overlap. So whenever

    min_breadth * G / L  >  min_unique_reads

the read floor CANNOT be the binding constraint - it is decorative, the true limit of
detection is set by breadth alone, and the number written in the watchlist column
understates the real LoD by exactly that ratio.

Measured on the shipped human pack (150 bp reads, genome sizes from the built refs
manifest), that ratio reaches 243x for Schistosoma mansoni: the file says 50 reads,
the gate needs ~12,151. Balamuthia mandrillaris, whose own watchlist note says "low
read counts are expected and still meaningful", needs ~3,159 reads to clear its own
breadth floor.

WHY THIS SCRIPT REPORTS RATHER THAN FIXES

The obvious fix is to recompute one gate from the other. That means choosing a
spread factor, and choosing it from reasoning rather than measurement would replace
one arbitrary set of numbers with another - which is the thing this project is trying
not to do. These thresholds are due to be set from the benchmark. Until then the
incoherence is recorded here, and tests/test_pipeline.py pins the counts so that a
change to the pack shows up as coherence improving or regressing rather than
silently staying broken.

Exit status: 0 always, unless --strict, which returns 1 if any taxon is incoherent.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def read_watchlist(path: Path) -> list[dict]:
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines(keepends=True)
             if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def read_genome_sizes(manifest: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    if not manifest.is_file():
        return out
    with manifest.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("role") != "target":
                continue
            try:
                size = int(row.get("genome_size") or 0)
            except ValueError:
                continue
            if size > 0:
                out.setdefault(row["name"], size)
    return out


def analyse(rows: list[dict], sizes: dict[str, int], read_len: int) -> list[dict]:
    out = []
    for r in rows:
        genome = sizes.get(r.get("name", ""))
        if not genome:
            continue
        try:
            min_reads = float(r["min_unique_reads"])
            min_breadth = float(r["min_breadth"])
        except (KeyError, ValueError):
            continue
        if min_reads <= 0:
            continue
        reads_for_breadth = min_breadth * genome / read_len
        out.append({
            "taxon": r["name"], "group": r.get("group", ""),
            "category": r.get("category", ""), "genome_bp": genome,
            "min_unique_reads": int(min_reads), "min_breadth": min_breadth,
            "reads_needed_for_breadth": round(reads_for_breadth, 1),
            "ratio_breadth_over_reads": round(reads_for_breadth / min_reads, 2),
            "binding_gate": "breadth" if reads_for_breadth > min_reads else "reads",
        })
    out.sort(key=lambda d: -d["ratio_breadth_over_reads"])
    return out


def summarise(rows: list[dict], tol: float) -> dict:
    return {
        "taxa_assessed": len(rows),
        "breadth_binds": sum(1 for r in rows
                             if r["ratio_breadth_over_reads"] > tol),
        "reads_bind": sum(1 for r in rows
                          if r["ratio_breadth_over_reads"] < 1.0 / tol),
        "coherent": sum(1 for r in rows
                        if 1.0 / tol <= r["ratio_breadth_over_reads"] <= tol),
        "worst_ratio": max((r["ratio_breadth_over_reads"] for r in rows), default=0.0),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pack", default="human")
    ap.add_argument("--packs-dir", type=Path, default=REPO / "packs")
    # Derived from $BTF_DB and from --pack, rather than one machine's volume with the
    # pack name baked in: the old default named packs/human whatever --pack said.
    _env = os.environ.get("BTF_DB")
    ap.add_argument("--refs-manifest", type=Path, default=None,
                    help="built refs manifest, for real genome sizes (default: "
                         f"$BTF_DB/packs/<pack>/refs/manifest.tsv, here "
                         f"{_env + '/packs/<pack>/refs/manifest.tsv' if _env else 'unset'})")
    ap.add_argument("--read-length", type=int, default=150,
                    help="the read length the LoD claim is made at [150]")
    ap.add_argument("--tolerance", type=float, default=1.5,
                    help="a gate within this factor of the other counts as coherent")
    ap.add_argument("--out", type=Path, default=None, help="write the table as TSV")
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 if any taxon's gates are incoherent")
    args = ap.parse_args(argv)

    wl = args.packs_dir / args.pack / "watchlist.tsv"
    if not wl.is_file():
        sys.stderr.write(f"ERROR: no watchlist at {wl}\n")
        return 2
    refs_manifest = args.refs_manifest
    if refs_manifest is None:
        btf_db = os.environ.get("BTF_DB")
        if not btf_db:
            sys.stderr.write(
                "ERROR: no --refs-manifest and $BTF_DB is not in the environment.\n"
                "       export BTF_DB=/path/to/project/bioThreatFinder, or pass\n"
                "       --refs-manifest PATH.\n")
            return 2
        refs_manifest = (Path(btf_db) / "packs" / args.pack / "refs" / "manifest.tsv")
    sizes = read_genome_sizes(refs_manifest)
    if not sizes:
        sys.stderr.write(
            f"ERROR: no genome sizes from {refs_manifest}.\n"
            f"       Without real genome sizes nothing here can be computed, and\n"
            f"       guessing them would produce a confident wrong answer.\n"
            f"       Build the pack refs first, or pass --refs-manifest.\n")
        return 2

    rows = analyse(read_watchlist(wl), sizes, args.read_length)
    stats = summarise(rows, args.tolerance)

    print(f"pack={args.pack} read_length={args.read_length} "
          f"tolerance={args.tolerance}x")
    print(f"{'taxon':38} {'group':9} {'genome_bp':>12} {'min_rd':>7} "
          f"{'min_br':>8} {'rd_needed':>10} {'ratio':>7} {'binds':>8}")
    for r in rows:
        print(f"{r['taxon'][:37]:38} {r['group']:9} {r['genome_bp']:>12,} "
              f"{r['min_unique_reads']:>7} {r['min_breadth']:>8.4f} "
              f"{r['reads_needed_for_breadth']:>10,.0f} "
              f"{r['ratio_breadth_over_reads']:>7.1f} {r['binding_gate']:>8}")
    print()
    print(f"taxa assessed            : {stats['taxa_assessed']}")
    print(f"breadth is the real gate : {stats['breadth_binds']}  "
          f"(min_unique_reads is decorative for these)")
    print(f"reads are the real gate  : {stats['reads_bind']}  "
          f"(min_breadth is decorative for these)")
    print(f"roughly coherent         : {stats['coherent']}")
    print(f"worst ratio              : {stats['worst_ratio']:.1f}x")
    print()
    print("A ratio of N means the per-taxon read floor understates that taxon's true")
    print("limit of detection by N-fold, because breadth cannot be reached with the")
    print("declared number of reads. These numbers are the input to threshold")
    print("calibration; they are not themselves a calibration.")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {args.out}")

    if args.strict and stats["coherent"] != stats["taxa_assessed"]:
        sys.stderr.write(
            f"\nERROR (--strict): {stats['taxa_assessed'] - stats['coherent']} taxon/taxa "
            f"have mutually incoherent read and breadth gates.\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
