#!/usr/bin/env python3
"""Can the Metabuli database report every taxon on the watchlist?

The last arm with no check at all. D47 closed the marker database, D46 the sketch, D48 the
Kraken2 index; Metabuli stayed unchecked because it is opt-in and off by default
(`--skip_metabuli` defaults true) - which is a reason it was not urgent, not a reason it was
covered. An arm nobody checks is an arm that answers nothing when it is finally switched on
for the question it exists for: a divergent agent that the nucleotide arm would miss.

WHY NODE PRESENCE IS THE WRONG TEST HERE, UNLIKE KRAKEN2
--------------------------------------------------------
`check_index_coverage.py` asks whether a watchlist taxon has a node in the Kraken2 index's
PRUNED taxonomy, and that works because kraken2-build keeps only nodes spanning sequences it
loaded - so an absent node proves the arm cannot report the taxon.

A Metabuli database copies the **whole** NCBI taxdump into `ncbi-taxdump/`. Every watchlist
taxon has a node there, and so would every taxon that has ever existed. Testing node
presence would therefore pass unconditionally: a check that cannot fail is not a check.

The discriminating input is `taxID_list`, the taxids the database actually holds sequence
for. So the question is asked over the CLADE: does this taxon, or any descendant of it,
appear in `taxID_list`? A taxon with no sequenced descendant can never be reported, exactly
as with an absent Kraken2 node - the arm runs, `screen_summary.tsv` records it `ok`, and the
row simply never exists.

Descendants are expanded from the database's OWN `ncbi-taxdump/nodes.dmp`, not from the
pinned snapshot. That is deliberate: the question is what THIS database can report, and it
will use its own copy of the taxonomy at query time. Using the pinned snapshot instead would
answer a question about a database that is not installed.

WHAT IT DOES NOT CLAIM
----------------------
A sequenced descendant is necessary, not sufficient. One genome under a genus does not mean
a given species can be distinguished, and Metabuli's amino-acid route can place a read on a
clade without resolving the species. So an absent clade PROVES the arm is blind to that
taxon, while a present one only says the arm has something to answer with.

    scripts/check_metabuli_coverage.py --pack human
      exit 0  every watchlist taxon has a sequenced clade in this database
      exit 3  at least one cannot be reported by this arm
      exit 1  could not tell (no database, unreadable, or a namespace mismatch)
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

EXIT_OK, EXIT_UNDETERMINED, EXIT_GAP = 0, 1, 3


def read_watchlist(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if (r.get("taxid") or "").strip().isdigit()]


def read_taxid_list(path: Path) -> set[int]:
    """The taxids this database holds sequence for."""
    out: set[int] = set()
    with path.open() as fh:
        for line in fh:
            tok = line.strip().split("\t")[0].strip()
            if tok.isdigit():
                out.add(int(tok))
    return out


def read_parents(nodes_dmp: Path) -> dict[int, int]:
    """child -> parent, from the database's own copy of the taxonomy."""
    parents: dict[int, int] = {}
    with nodes_dmp.open() as fh:
        for line in fh:
            # nodes.dmp is "taxid\t|\tparent\t|\trank\t|..."
            parts = line.split("\t|\t", 2)
            if len(parts) < 2:
                continue
            a, b = parts[0].strip(), parts[1].strip()
            if a.isdigit() and b.isdigit():
                parents[int(a)] = int(b)
    return parents


def sequenced_clades(seq_taxids: set[int], parents: dict[int, int]) -> set[int]:
    """Every taxid that has sequence at or below it.

    Walked UP from each sequenced taxid rather than down from each watchlist taxon: the
    taxonomy has ~3M nodes and 35k sequenced ones, so this is one pass over the small set
    instead of 90 descents through the large one. Ancestors are memoised, so a deep lineage
    is walked once however many sequenced genomes sit under it.
    """
    covered: set[int] = set()
    for t in seq_taxids:
        node, chain = t, []
        while node and node not in covered:
            chain.append(node)
            nxt = parents.get(node)
            if nxt is None or nxt == node:        # root, or a self-parent (taxdump does this at 1)
                break
            node = nxt
        covered.update(chain)
    return covered


def self_test() -> int:
    """Synthetic taxonomy: 1 -> 10 -> {11, 12}, 1 -> 20 -> 21. Sequence under 11 and 21."""
    parents = {1: 1, 10: 1, 11: 10, 12: 10, 20: 1, 21: 20}
    seq = {11, 21}
    covered = sequenced_clades(seq, parents)
    checks = [
        ("a sequenced taxon is covered", 11 in covered),
        ("its ancestor is covered", 10 in covered and 1 in covered),
        ("a SIBLING with no sequence is NOT covered", 12 not in covered),
        ("the other branch is covered through its own leaf", 20 in covered and 21 in covered),
        # The property that makes the walk correct: memoisation must not stop a later
        # lineage short. 21's walk hits 1, which 11's walk already added - if the loop
        # broke on `node in covered` before recording 20, this would fail.
        ("a shared ancestor does not truncate a later lineage", 20 in covered),
    ]
    bad = [n for n, ok in checks if not ok]
    for n, ok in checks:
        print(f"  {'ok  ' if ok else 'FAIL'} {n}")
    print("\nSELF-TEST " + ("PASSED" if not bad else f"FAILED ({len(bad)})"))
    return EXIT_OK if not bad else EXIT_GAP


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pack", default="human")
    ap.add_argument("--packs-dir", type=Path,
                    default=Path(__file__).resolve().parent.parent / "packs")
    _env = os.environ.get("BTF_DB")
    ap.add_argument("--btf-db", type=Path, default=Path(_env) if _env else None,
                    help="base for project-owned artifacts (default: $BTF_DB, here "
                         f"{_env or 'unset'})")
    ap.add_argument("--metabuli-db", type=Path, default=None,
                    help="the database directory holding taxID_list "
                         "[<btf-db>/metabuli/<mode>]")
    ap.add_argument("--mode", default="refseq", help="which database (default: refseq)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()

    if args.metabuli_db is None and args.btf_db is None:
        print("UNDETERMINED: no --metabuli-db or --btf-db, and $BTF_DB is not set.\n"
              "  export BTF_DB=/path/to/project/bioThreatFinder, or pass --metabuli-db DIR.")
        return EXIT_UNDETERMINED

    db = args.metabuli_db or (args.btf_db / "metabuli" / args.mode)
    watchlist = args.packs_dir / args.pack / "watchlist.resolved.tsv"
    if not watchlist.exists():
        print(f"UNDETERMINED: no resolved watchlist at {watchlist}")
        return EXIT_UNDETERMINED

    taxid_list = db / "taxID_list"
    nodes = db / "ncbi-taxdump" / "nodes.dmp"
    for p in (taxid_list, nodes):
        if not p.is_file():
            print(f"UNDETERMINED: no Metabuli database to check - missing {p}")
            print("  This arm is opt-in and off by default, so an absent database is a")
            print("  configuration, not a fault. Build it with scripts/build_metabuli_db.sh")
            print("  or leave --skip_metabuli at its default.")
            return EXIT_UNDETERMINED

    # THE NAMESPACE GUARD, RESTATED HERE. Metabuli publishes GTDB-based databases too, and
    # GTDB taxids are not NCBI taxids: point the pipeline at one and every watchlist lookup
    # misses, producing an empty result rather than an error. main.nf refuses such an index
    # at launch; this check must not quietly report coverage against the wrong namespace.
    manifest = db / "BTF_BUILD_MANIFEST.txt"
    if manifest.is_file():
        text = manifest.read_text(encoding="utf-8", errors="replace")
        line = next((l for l in text.splitlines()
                     if l.startswith("taxonomy_namespace=")), "")
        if line and "NCBI" not in line.upper():
            print(f"UNDETERMINED: {manifest.name} reports {line.strip()}")
            print("  Coverage against a non-NCBI namespace is not comparable to the")
            print("  watchlist, which carries NCBI taxids. Nothing is reported.")
            return EXIT_UNDETERMINED
    else:
        print(f"NOTE: no {manifest.name}; the namespace could not be confirmed as NCBI.")

    rows = read_watchlist(watchlist)
    seq = read_taxid_list(taxid_list)
    parents = read_parents(nodes)
    if not seq or not parents:
        print(f"UNDETERMINED: read {len(seq)} sequenced taxid(s) and {len(parents)} "
              f"taxonomy node(s); one of them is empty, so nothing can be concluded.")
        return EXIT_UNDETERMINED

    covered = sequenced_clades(seq, parents)
    print(f"database : {db} ({len(seq)} sequenced taxids, "
          f"{len(parents)} taxonomy nodes, {len(covered)} clades with sequence)")
    print(f"watchlist: {watchlist} ({len(rows)} taxa)\n")

    absent = [r for r in rows if int(r["taxid"]) not in covered]
    if not absent:
        print(f"verdict: COVERED. All {len(rows)} watchlist taxa have sequence at or below "
              f"them in this database.")
        print("         A sequenced clade is necessary, not sufficient: it says the arm has "
              "something to\n         answer with, not that a species can be resolved.")
        return EXIT_OK

    print(f"{len(absent)} of {len(rows)} watchlist taxa have NO sequence at or below them "
          f"and can never appear in this arm's output:\n")
    print(f"    {'taxid':>8s}  {'name':<34s} {'category':<17s} rank")
    for r in sorted(absent, key=lambda r: (r.get("category", ""), r.get("name", ""))):
        print(f"    {r['taxid']:>8s}  {(r.get('name') or '')[:34]:<34s} "
              f"{(r.get('category') or '')[:17]:<17s} {r.get('resolved_rank') or ''}")
    print()
    print("verdict: GAP. This arm cannot report the taxa above. It is opt-in and off by")
    print("         default, so that is not in itself a run-time defect - but it is what")
    print("         the arm would answer if switched on for them, and the reason to switch")
    print("         it on is a divergent agent the nucleotide arm would miss.")
    print("         This is a property of the database's contents, not staleness.")
    return EXIT_GAP


if __name__ == "__main__":
    sys.exit(main())
