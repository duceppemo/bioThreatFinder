#!/usr/bin/env python3
"""
Intersect stage-1 screen output with the watchlist to produce stage-2 candidates.

The important detail: Kraken2 assigns many reads to INTERNAL nodes. A read that
is consistent with several species in a genus is placed at the genus. So matching
a watchlist species taxid against the report's taxid column alone misses most of
the relevant signal - which is the difference between a screen that works and one
that quietly reports nothing.

This script therefore rolls up every descendant of each watchlist taxon
(precomputed by scripts/resolve_watchlist.py) and also walks UP from the
watchlist taxon to record ancestor-level assignments separately, flagged as
``ancestor_only``. Ancestor-level reads cannot confirm a species, but a genus-
level pile-up with zero species-level reads is exactly the pattern a divergent or
under-represented strain produces, so it must reach stage 2 rather than be
dropped.

Stage 1 is deliberately permissive. Stage 2 does the discriminating.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    if not path or not path.exists() or path.name == "NO_FILE":
        return []
    with path.open() as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t")) if lines else []


def parse_kraken2_report(path: Path) -> dict[str, dict]:
    """taxid -> {clade_reads, direct_reads, rank, name}.

    Kraken2 report columns: pct, clade_reads, direct_reads, rank, taxid, name.
    A 6-column variant and an 8-column (--report-minimizer-data) variant both
    occur in the wild, so index from the right-hand end where the layout is
    stable rather than assuming a width.
    """
    out: dict[str, dict] = {}
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            name = f[-1].strip()
            taxid = f[-2].strip()
            rank = f[-3].strip()
            if not taxid.isdigit():
                continue
            try:
                clade = int(f[1])
                direct = int(f[2])
            except (ValueError, IndexError):
                continue
            out[taxid] = {"clade_reads": clade, "direct_reads": direct,
                          "rank": rank, "name": name}
    return out


def parse_sylph(path: Path) -> dict[str, float]:
    """Genome file -> taxonomic abundance. Keyed loosely; sylph reports genome
    paths, not taxids, so match on any taxid-looking token in the filename."""
    out: dict[str, float] = {}
    if not path.exists() or path.name == "NO_FILE":
        return out
    rows = read_tsv(path)
    for row in rows:
        ab = row.get("Taxonomic_abundance") or row.get("Sequence_abundance") or "0"
        gf = row.get("Genome_file") or row.get("Contig_name") or ""
        try:
            val = float(ab)
        except ValueError:
            continue
        for token in gf.replace("/", "_").replace(".", "_").split("_"):
            if token.isdigit() and len(token) >= 2:
                out[token] = max(out.get(token, 0.0), val)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kraken2-report", required=True, type=Path)
    ap.add_argument("--bracken-report", type=Path, default=None)
    ap.add_argument("--sylph-profile", type=Path, default=None)
    ap.add_argument("--watchlist", required=True, type=Path, help="watchlist.resolved.tsv")
    ap.add_argument("--descendants", required=True, type=Path,
                    help="watchlist_taxid -> descendant_taxid TSV")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--role", default="sample")
    ap.add_argument("--platform", default="illumina")
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    watch = read_tsv(args.watchlist)
    if not watch:
        sys.exit(f"ERROR: empty watchlist {args.watchlist}")

    desc_of: dict[str, set[str]] = defaultdict(set)
    with args.descendants.open() as fh:
        header = next(fh, "")
        if "descendant" not in header:
            fh.seek(0)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                desc_of[parts[0]].add(parts[2])

    report = parse_kraken2_report(args.kraken2_report)
    bracken = parse_kraken2_report(args.bracken_report) if args.bracken_report else {}
    sylph = parse_sylph(args.sylph_profile) if args.sylph_profile else {}

    if not report:
        sys.stderr.write(
            f"[screen] WARNING: no usable rows in {args.kraken2_report}. "
            "No candidates will be produced for this sample.\n"
        )

    rows = []
    for entry in watch:
        taxid = entry["taxid"]
        family = desc_of.get(taxid) or {taxid}

        assigned = 0
        direct = 0
        matched_nodes = []
        for d in family:
            hit = report.get(d)
            if hit and hit["direct_reads"] > 0:
                direct += hit["direct_reads"]
                matched_nodes.append(f"{d}:{hit['direct_reads']}")
        # Clade count at the watchlist node itself captures internal-node
        # assignments in one number; prefer it when present.
        self_hit = report.get(taxid)
        assigned = max(direct, self_hit["clade_reads"] if self_hit else 0)

        ancestor_only = bool(self_hit and self_hit["clade_reads"] > 0 and direct == 0)

        bracken_reads = 0
        for d in family:
            b = bracken.get(d)
            if b:
                bracken_reads += b["clade_reads"]

        sylph_ab = sylph.get(taxid, 0.0)

        if assigned < args.min_reads and sylph_ab <= 0:
            continue

        rows.append({
            "sample": args.sample,
            "role": args.role,
            "platform": args.platform,
            "taxid": taxid,
            "name": entry["name"],
            "group": entry.get("group", ""),
            "category": entry.get("category", ""),
            "assigned_reads": assigned,
            "direct_reads": direct,
            "bracken_reads": bracken_reads,
            "sylph_abundance": f"{sylph_ab:.6f}" if sylph_ab else "",
            "ancestor_only": "yes" if ancestor_only else "no",
            "matched_nodes": ";".join(sorted(matched_nodes)[:20]),
            "n_descendants_searched": len(family),
            "near_neighbour_taxids": entry.get("near_neighbour_taxids", ""),
            "min_unique_reads": entry.get("min_unique_reads", ""),
            "min_breadth": entry.get("min_breadth", ""),
            "notes": entry.get("notes", ""),
        })

    rows.sort(key=lambda r: -r["assigned_reads"])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0]) if rows else [
        "sample", "role", "platform", "taxid", "name", "group", "category",
        "assigned_reads", "direct_reads", "bracken_reads", "sylph_abundance",
        "ancestor_only", "matched_nodes", "n_descendants_searched",
        "near_neighbour_taxids", "min_unique_reads", "min_breadth", "notes",
    ]
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    sys.stderr.write(
        f"[screen] {args.sample}: {len(rows)} of {len(watch)} watchlist taxa raised as "
        f"stage-2 candidates (min_reads={args.min_reads})\n"
    )
    for r in rows[:15]:
        flag = " [ancestor-level only]" if r["ancestor_only"] == "yes" else ""
        sys.stderr.write(f"[screen]   {r['name']} (taxid {r['taxid']}): "
                         f"{r['assigned_reads']} reads{flag}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
