#!/usr/bin/env python3
"""
Resolve a bioThreatFinder watchlist from names to NCBI taxonomy IDs.

Why this exists
---------------
Watchlists keyed by taxid rot. NCBI merges, splits and retires nodes; a stale
taxid silently stops matching and the organism drops out of screening with no
error. So watchlists are keyed by NAME and resolved against a pinned local
taxdump at build time. Anything that does not resolve to exactly one node is a
hard error with candidates printed, so a human fixes it deliberately.

Also emits, for each watchlist taxon, the full set of descendant taxids. Kraken2
assigns reads to internal nodes (often genus or above), so matching only the
species taxid misses most of the signal.

Usage
-----
  resolve_watchlist.py --watchlist packs/human/watchlist.tsv \\
                       --taxdump /media/36tb/db/taxdump \\
                       --out packs/human/watchlist.resolved.tsv
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
from collections import defaultdict
from pathlib import Path

REQUIRED_COLUMNS = [
    "name", "rank", "group", "category",
    "min_unique_reads", "min_breadth", "near_neighbours", "notes", "source",
]

# NCBI name classes we accept, in order of preference. A scientific name always
# wins; synonyms are accepted because virus and fungal taxonomy churns and the
# official disease lists lag NCBI by years.
NAME_CLASS_PREFERENCE = [
    "scientific name",
    "equivalent name",
    "synonym",
    "genbank synonym",
    "includes",
]


def read_watchlist(path: Path) -> list[dict]:
    rows = []
    with path.open() as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    if not lines:
        sys.exit(f"ERROR: {path} contains no non-comment lines")
    reader = csv.DictReader(lines, delimiter="\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
    if missing:
        sys.exit(f"ERROR: {path} is missing required columns: {', '.join(missing)}")
    for i, row in enumerate(reader, start=2):
        name = (row.get("name") or "").strip()
        if not name:
            continue
        row["name"] = name
        row["_line"] = i
        rows.append(row)
    if not rows:
        sys.exit(f"ERROR: {path} yielded no watchlist entries")
    return rows


def taxdump_age_days(taxdump: Path) -> tuple[float | None, str]:
    """Age of the taxonomy snapshot in days, and where that date came from.

    Prefers PROVENANCE.txt's recorded NCBI Last-Modified, because a directory's file
    mtimes say when it was unpacked locally, not when NCBI generated it.
    """
    prov = taxdump / "PROVENANCE.txt"
    if prov.exists():
        for line in prov.read_text().splitlines():
            if line.lower().startswith("ncbi_last_modified"):
                raw = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else ""
                for fmt in ("%a, %d %b %Y %H:%M:%S %Z", "%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%d"):
                    try:
                        when = datetime.datetime.strptime(raw, fmt)
                        return (datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None) - when).total_seconds() / 86400.0, \
                               f"PROVENANCE.txt ncbi_last_modified={raw}"
                    except ValueError:
                        continue
    newest = 0.0
    for f in ("nodes.dmp", "names.dmp"):
        fp = taxdump / f
        if fp.exists():
            newest = max(newest, fp.stat().st_mtime)
    if not newest:
        return None, "unknown"
    age = (datetime.datetime.now(datetime.timezone.utc).timestamp() - newest) / 86400.0
    return age, f"file mtime of nodes.dmp/names.dmp ({datetime.date.fromtimestamp(newest)})"


def load_names(taxdump: Path) -> dict[str, list[tuple[int, str]]]:
    """Map lowercased name -> [(taxid, name_class), ...]."""
    names_dmp = taxdump / "names.dmp"
    if not names_dmp.exists():
        sys.exit(f"ERROR: {names_dmp} not found. Point --taxdump at an unpacked NCBI taxdump.")
    index: dict[str, list[tuple[int, str]]] = defaultdict(list)
    with names_dmp.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            try:
                taxid = int(parts[0])
            except ValueError:
                continue
            index[parts[1].lower()].append((taxid, parts[3]))
    return index


def load_sci_names(taxdump: Path) -> dict[int, str]:
    """Map taxid -> scientific name, in the taxdump's own casing.

    load_names() lowercases its keys, which is right for lookup and wrong for output:
    the descendants file is joined against benchmark truth columns and watchlist
    entries, both of which carry NCBI casing. Reversing the lowercased index would
    emit "bundibugyo ebolavirus" and match nothing.
    """
    out: dict[int, str] = {}
    with (taxdump / "names.dmp").open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4 or parts[3] != "scientific name":
                continue
            try:
                out[int(parts[0])] = parts[1]
            except ValueError:
                continue
    return out


def load_nodes(taxdump: Path) -> tuple[dict[int, int], dict[int, str]]:
    """Return (child -> parent, taxid -> rank)."""
    nodes_dmp = taxdump / "nodes.dmp"
    if not nodes_dmp.exists():
        sys.exit(f"ERROR: {nodes_dmp} not found")
    parent: dict[int, int] = {}
    rank: dict[int, str] = {}
    with nodes_dmp.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            try:
                tid, pid = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            parent[tid] = pid
            rank[tid] = parts[2]
    return parent, rank


def load_merged(taxdump: Path) -> dict[int, int]:
    """Old taxid -> current taxid, so historical lists keep working."""
    merged = taxdump / "merged.dmp"
    out: dict[int, int] = {}
    if not merged.exists():
        return out
    with merged.open() as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 2:
                try:
                    out[int(parts[0])] = int(parts[1])
                except ValueError:
                    pass
    return out


def build_children(parent: dict[int, int]) -> dict[int, list[int]]:
    children: dict[int, list[int]] = defaultdict(list)
    for child, par in parent.items():
        if child != par:
            children[par].append(child)
    return children


def descendants(taxid: int, children: dict[int, list[int]]) -> set[int]:
    """Iterative DFS - the tree is deep enough that recursion is a real risk."""
    seen = {taxid}
    stack = [taxid]
    while stack:
        node = stack.pop()
        for kid in children.get(node, ()):
            if kid not in seen:
                seen.add(kid)
                stack.append(kid)
    return seen


def resolve_one(name: str, names_idx, rank_of, expected_rank: str):
    """Resolve a single name. Returns (taxid, note) or raises ValueError."""
    hits = names_idx.get(name.lower())
    if not hits:
        raise ValueError(f"no taxonomy node matches name {name!r}")

    # Collapse to unique taxids, keeping the best name class per taxid.
    best: dict[int, str] = {}
    for taxid, name_class in hits:
        cur = best.get(taxid)
        if cur is None:
            best[taxid] = name_class
            continue
        def pref(nc: str) -> int:
            return NAME_CLASS_PREFERENCE.index(nc) if nc in NAME_CLASS_PREFERENCE else 99
        if pref(name_class) < pref(cur):
            best[taxid] = name_class

    accepted = {t: nc for t, nc in best.items() if nc in NAME_CLASS_PREFERENCE}
    if not accepted:
        classes = sorted(set(best.values()))
        raise ValueError(
            f"name {name!r} exists but only as {classes}; not an accepted name class"
        )

    if len(accepted) == 1:
        taxid, name_class = next(iter(accepted.items()))
        note = "" if name_class == "scientific name" else f"matched via {name_class}"
        return taxid, note

    # Ambiguous. Try to break the tie on the expected rank.
    by_rank = {t: nc for t, nc in accepted.items()
               if rank_of.get(t, "no rank") == expected_rank}
    if len(by_rank) == 1:
        taxid, name_class = next(iter(by_rank.items()))
        note = f"ambiguous name disambiguated by rank={expected_rank}"
        return taxid, note

    # Then on name class - a scientific name beats a synonym.
    sci = {t: nc for t, nc in accepted.items() if nc == "scientific name"}
    if len(sci) == 1:
        taxid = next(iter(sci))
        return taxid, "ambiguous name disambiguated by scientific-name class"

    detail = ", ".join(
        f"{t} (rank={rank_of.get(t, '?')}, class={nc})"
        for t, nc in sorted(accepted.items())
    )
    raise ValueError(f"name {name!r} is ambiguous across {len(accepted)} nodes: {detail}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--watchlist", required=True, type=Path)
    ap.add_argument("--taxdump", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--descendants-out", type=Path,
                    help="optional TSV of watchlist_taxid -> descendant taxid")
    ap.add_argument("--max-taxdump-age-days", type=float, default=90.0,
                    help="warn when the taxonomy snapshot is older than this [90]")
    ap.add_argument("--require-fresh-taxdump", action="store_true",
                    help="fail instead of warning when the snapshot is stale")
    ap.add_argument("--allow-unresolved", action="store_true",
                    help="warn instead of failing on unresolved names (NOT for operational use)")
    args = ap.parse_args()

    rows = read_watchlist(args.watchlist)
    sys.stderr.write(f"[resolve] {len(rows)} watchlist entries from {args.watchlist}\n")

    age, how = taxdump_age_days(args.taxdump)
    if age is None:
        sys.stderr.write(f"[resolve] WARNING: cannot determine the age of {args.taxdump}\n")
    elif age > args.max_taxdump_age_days:
        msg = (
            f"taxonomy snapshot is {age:.0f} days old (from {how}), above the "
            f"{args.max_taxdump_age_days:.0f}-day limit.\n"
            "  NCBI regenerates the taxdump daily and adds roughly 40,000 taxids a\n"
            "  month. A stale snapshot does not fail loudly - it silently loses\n"
            "  descendant nodes inside watchlist clades, so reads a classifier assigns\n"
            "  to a newer node are never recognised as the watchlist organism. A\n"
            "  6.5-month-old snapshot was measured here as 288,764 taxids short of\n"
            "  current and 2,581 descendant nodes short inside watchlist clades.\n"
            "  It must also match the snapshot the k-mer index was built against.\n"
            "  Refresh: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/\n")
        if args.require_fresh_taxdump:
            sys.exit(f"ERROR: {msg}")
        sys.stderr.write(f"[resolve] WARNING: {msg}")
    else:
        sys.stderr.write(f"[resolve] taxonomy snapshot is {age:.0f} days old ({how})\n")

    sys.stderr.write(f"[resolve] loading taxonomy from {args.taxdump}\n")
    names_idx = load_names(args.taxdump)
    parent, rank_of = load_nodes(args.taxdump)
    merged = load_merged(args.taxdump)
    children = build_children(parent)
    sys.stderr.write(f"[resolve] {len(parent)} taxonomy nodes, {len(merged)} merged ids\n")

    resolved: list[dict] = []
    failures: list[str] = []

    for row in rows:
        try:
            taxid, note = resolve_one(row["name"], names_idx, rank_of, row["rank"].strip())
        except ValueError as exc:
            failures.append(f"  line {row['_line']}: {exc}")
            continue

        taxid = merged.get(taxid, taxid)
        actual_rank = rank_of.get(taxid, "no rank")
        if actual_rank != row["rank"].strip():
            note = (note + "; " if note else "") + \
                   f"rank is {actual_rank}, watchlist expected {row['rank'].strip()}"

        desc = descendants(taxid, children)
        # RETIRED TAXIDS BELONG IN THE CLADE TOO. merged.dmp was applied to the watchlist
        # taxid and to near neighbours, but never to the descendant set, which walks
        # nodes.dmp and therefore holds CURRENT ids only. When a classifier index predates
        # the taxdump - the normal direction, since indexes are rebuilt rarely - it assigns
        # reads to ids that have since been merged away, parse_screen matches those ids
        # against this set, finds nothing, and the taxon never becomes a candidate. A silent
        # missed detection with no trace in any output. Measured against the pinned
        # 2026-09-29 taxdump: 3,491 retired ids point into 50 of the 89 watchlist clades,
        # 2,728 of them inside Influenza A virus alone, plus 453 Salmonella enterica, 66
        # Listeria monocytogenes and 9 Enterovirus C.
        desc |= {old for old, new in merged.items() if new in desc or new == taxid}

        # Near neighbours are advisory - a missing one degrades competitive
        # mapping but must not block the run.
        nn_names = [n.strip() for n in (row.get("near_neighbours") or "").split(";")
                    if n.strip() and n.strip() != "-"]
        nn_taxids, nn_missing = [], []
        for nn in nn_names:
            try:
                # expected_rank was hardcoded "species", so a neighbour that is
                # legitimately a subspecies, strain or genus got no rank tie-break and an
                # ambiguous name landed in nn_missing - silently weakening competitive
                # mapping for exactly the target/neighbour pairs that matter. The live
                # watchlist already carries "Francisella tularensis subsp. novicida" (a
                # subspecies) and "Shigella" (a genus) as neighbours; they survive only
                # because those names happen to be unambiguous. Try species first, then
                # fall back to accepting any rank rather than discarding the neighbour.
                try:
                    nn_taxid, _ = resolve_one(nn, names_idx, rank_of, "species")
                except ValueError:
                    nn_taxid, _ = resolve_one(nn, names_idx, rank_of, "")
                nn_taxids.append(str(merged.get(nn_taxid, nn_taxid)))
            except ValueError:
                nn_missing.append(nn)
        if nn_missing:
            note = (note + "; " if note else "") + \
                   f"unresolved near neighbours: {', '.join(nn_missing)}"

        resolved.append({
            "taxid": taxid,
            "name": row["name"],
            "resolved_rank": actual_rank,
            "group": row["group"],
            "category": row["category"],
            "min_unique_reads": row["min_unique_reads"],
            "min_breadth": row["min_breadth"],
            "near_neighbour_taxids": ";".join(nn_taxids),
            "n_descendants": len(desc),
            "notes": row["notes"],
            "source": row["source"],
            "resolution_note": note,
            "_descendants": desc,
        })

    if failures:
        sys.stderr.write(
            f"\n{'WARNING' if args.allow_unresolved else 'ERROR'}: "
            f"{len(failures)} watchlist entries did not resolve:\n"
        )
        sys.stderr.write("\n".join(failures) + "\n")
        sys.stderr.write(
            "\nFix the names in the watchlist, or update the taxdump. Names must match\n"
            "NCBI taxonomy exactly (scientific name or accepted synonym).\n"
        )
        if not args.allow_unresolved:
            return 1

    # Warn on duplicates - two watchlist rows mapping to one node means
    # conflicting thresholds and a double-counted report.
    seen: dict[int, str] = {}
    for r in resolved:
        if r["taxid"] in seen:
            sys.stderr.write(
                f"WARNING: {r['name']!r} and {seen[r['taxid']]!r} both resolve to "
                f"taxid {r['taxid']}\n"
            )
        else:
            seen[r["taxid"]] = r["name"]

    # A HEADERLESS 2-BYTE FILE IS NOT AN OUTPUT.
    #
    # With --allow-unresolved and nothing resolvable, `cols` fell back to [] and
    # csv.DictWriter with empty fieldnames wrote only a line terminator - then exited 0.
    # scripts/run_bakeoff.sh gates on `[[ -s "$WATCHLIST" ]]`, and a 2-byte file is
    # non-empty, so it passed the check and the bake-off ran against an empty watchlist.
    if not resolved:
        sys.stderr.write(
            "[resolve] ERROR: nothing resolved, so no watchlist would be written.\n"
            "[resolve] Refusing to emit an empty file: it is non-empty to `test -s`, so "
            "downstream steps would accept it and screen for no organisms at all - a "
            "clean negative for the entire watchlist.\n")
        return 2
    cols = [c for c in resolved[0] if not c.startswith("_")]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in resolved:
            w.writerow(r)
    sys.stderr.write(f"[resolve] wrote {len(resolved)} resolved entries to {args.out}\n")

    if args.descendants_out:
        # The NAME goes out beside the taxid. parse_screen matches on taxids, but
        # evaluate_benchmark joins findings to benchmark truth on the taxon NAME, and
        # the manifest's truth_present column carries names only. Without this column
        # a group-rank finding cannot be connected to the species the truth names, so
        # the lead reads as a call on an unlisted taxon and the species it covers
        # still scores as a miss. A retired descendant with no current scientific
        # name keeps its row and an empty name - it is still a valid classifier
        # target, it just cannot participate in the name join.
        sci = load_sci_names(args.taxdump)
        with args.descendants_out.open("w") as fh:
            fh.write("watchlist_taxid\twatchlist_name\tdescendant_taxid"
                     "\tdescendant_name\n")
            for r in resolved:
                for d in sorted(r["_descendants"]):
                    fh.write(f"{r['taxid']}\t{r['name']}\t{d}\t{sci.get(d, '')}\n")
        total = sum(len(r["_descendants"]) for r in resolved)
        sys.stderr.write(
            f"[resolve] wrote {total} descendant mappings to {args.descendants_out}\n"
        )

    notes = [r for r in resolved if r["resolution_note"]]
    if notes:
        sys.stderr.write(f"\n[resolve] {len(notes)} entries resolved with notes:\n")
        for r in notes:
            sys.stderr.write(f"  {r['name']} -> {r['taxid']}: {r['resolution_note']}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
