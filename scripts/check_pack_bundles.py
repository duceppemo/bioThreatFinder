#!/usr/bin/env python3
"""Are the competitive bundles current? Asked without invoking the builder.

Every other database here can be interrogated read-only: D46 gave the sketch `--check`,
D47 the marker database, D48 the Kraken2 index, D60 the Metabuli one. The competitive
bundles - the thing stage 2 actually aligns against - were the exception. Their staleness
predicate (D45, `bundle_staleness` in build_pack_refs.py) existed and was correct, but it
lived halfway down a 43 GB download tool, after the watchlist has been resolved against
NCBI and the per-taxon subtree queries have gone out. So "are my bundles current?" was a
question an operator could only ask by starting the builder, and `--dry-run` prints the
plan without ever reaching the predicate. A check nobody can run without side effects is a
check nobody runs.

THE PREDICATE IS IMPORTED, NOT REIMPLEMENTED
--------------------------------------------
`bundle_staleness` and `index_staleness` come from build_pack_refs.py itself. A second copy
would be the worse failure of the two available: the check would answer "current" with the
builder's old rules while the builder rebuilt, or cry stale on a pack the builder is happy
with, and either way the operator would learn to ignore it. One predicate, two callers.

WHAT THE SIDECAR CANNOT DO, AND WHY --deep IS THE DEFAULT
---------------------------------------------------------
D45's predicate reads `<taxid>.inputs.tsv`, the record of what went into a bundle. It is
exact - but only 10 of the human pack's 90 bundles have one, because the other 80 were
built before the sidecar existed. For those, D45 falls back to mtime, which cannot see a
restore from backup, a copy that preserved timestamps, or a re-download of the same
accession. Metadata alone therefore answers "could not tell" for 89% of this pack, and a
check whose default answer is "could not tell" teaches nothing.

So the default reads the bundles. The manifest records `genome_size` and `n_contigs` per
accession, computed by the same `fasta_stats` counting rule that this re-applies, and
`append_with_taxid_prefix` copies sequence lines through untouched - so per-taxid-prefix
contig counts and base totals inside a bundle must equal the manifest's sums exactly. That
comparison needs no sidecar, and it is stronger than one: the sidecar records what the
builder believed it wrote, while this reads what is there. The bundles are 2.1 GB
compressed and the pass takes about a minute; `--quick` skips it for the metadata-only
answer.

Reading the bundles also re-runs the two post-write checks that currently only ever fire on
a REBUILD, which means no bundle built before they existed has ever been subjected to them:

  * identical sequence under two taxid prefixes - the MAPQ-0 inversion of D40. Both copies
    match every read equally, the MAPQ filter removes them, `unique_reads` collapses and
    the taxon reads as attributed_to_neighbour. A systematic false negative, worst on
    exactly the target/neighbour pairs the bundle exists to discriminate;
  * a header with no `taxid|` prefix - a sequence no read can be attributed to, which for a
    segmented virus means losing a segment silently.

WHAT IT DOES NOT ANSWER
-----------------------
Whether the manifest agrees with `refs.lock.tsv`. That is a different question with a
different authority - the lock outranks the disk (D53) - and it is answered by
`build_pack_refs.py --from-lock`, which must reach NCBI to answer it. This check takes the
manifest as given and asks only whether the bundles are consistent with it.

    scripts/check_pack_bundles.py --pack human
      exit 0  every bundle and index is current
      exit 3  at least one is stale, incomplete or internally wrong
      exit 1  could not tell (no pack on disk, unreadable manifest, or --quick with
              no input record to compare against)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import os
import sys
from pathlib import Path

EXIT_OK, EXIT_UNDETERMINED, EXIT_GAP = 0, 1, 3


def load_builder():
    """build_pack_refs.py as a module, for its staleness predicates.

    Loaded by path rather than by name: `scripts/` is not a package, and putting it on
    sys.path would make every stdlib import in this process resolve against a directory
    full of scripts.
    """
    src = Path(__file__).resolve().parent / "build_pack_refs.py"
    spec = importlib.util.spec_from_file_location("btf_build_pack_refs", src)
    if spec is None or spec.loader is None:
        raise OSError(f"cannot load {src}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def bundle_prefix_stats(comp: Path) -> tuple[dict[str, tuple[int, int]], list[str], set]:
    """One pass over a bundle: per-taxid-prefix (contigs, bases), unprefixed headers,
    and sequences that appear under two different prefixes.

    Counted with `len(line.strip())` on non-header lines, which is `fasta_stats`'s rule
    exactly - the manifest's `genome_size` was produced by it, so any other rule here
    would disagree by a few bytes per contig and report every bundle stale.

    Sequence identity is hashed incrementally rather than buffered: the builder's own
    post-write check joins the record into one string first, which is fine for a bacterium
    and 1.5 Gbp for the largest bundle here.
    """
    stats: dict[str, tuple[int, int]] = {}
    unprefixed: list[str] = []
    seen: dict[str, str] = {}
    collisions: set[tuple[str, str]] = set()
    cur: str | None = None
    h = None

    def flush() -> None:
        if cur is None or h is None:
            return
        digest = h.hexdigest()
        prev = seen.get(digest)
        if prev is not None and prev != cur:
            collisions.add(tuple(sorted((prev, cur))))
        else:
            seen.setdefault(digest, cur)

    with gzip.open(comp, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                bits = line[1:].split("|", 2)
                if len(bits) >= 2 and bits[0] == "taxid" and bits[1].strip().isdigit():
                    cur = bits[1].strip()
                else:
                    cur = "?"
                    if len(unprefixed) < 5:
                        unprefixed.append(line[1:].rstrip()[:60])
                n, b = stats.get(cur, (0, 0))
                stats[cur] = (n + 1, b)
                h = hashlib.md5()
            else:
                seq = line.strip()
                if cur is not None:
                    n, b = stats[cur]
                    stats[cur] = (n, b + len(seq))
                if h is not None:
                    h.update(seq.encode())
        flush()
    return stats, unprefixed, collisions


def plan(builder, watchlist: list[dict], manifest: list[dict]) -> tuple[list[dict], list[str]]:
    """What each bundle should contain, from the watchlist and the manifest alone.

    This mirrors the builder's own construction, including the order it writes in, because
    the order decides the de-duplication: target genomes go in first, so when a neighbour
    taxid resolved to the SAME assembly as its target - observed for F. tularensis 263 vs
    subsp. novicida 264, both GCF_000833355.1 - the target keeps it and the neighbour is
    dropped. Expecting the dropped copy would report every such bundle stale forever.
    """
    rows_by_taxid: dict[str, list[dict]] = {}
    notes: list[str] = []
    for r in manifest:
        tid = (r.get("taxid") or "").strip()
        p = (r.get("path") or "").strip()
        if not tid or not p:
            continue
        if not Path(p).is_file():
            notes.append(f"manifest names {Path(p).name} for taxid {tid}, which is not on "
                         f"disk")
            continue
        rows_by_taxid.setdefault(tid, []).append(r)

    accession_by_path = {str(Path(r["path"]).resolve()): (r.get("accession") or "").strip()
                         for rs in rows_by_taxid.values() for r in rs
                         if (r.get("accession") or "").strip()}

    out: list[dict] = []
    for w in watchlist:
        tid = (w.get("taxid") or "").strip()
        if not tid.isdigit():
            continue
        target_rows = rows_by_taxid.get(tid, [])
        neighbours = [str(n) for n in builder.parse_neighbours(
            w.get("near_neighbour_taxids", ""))]
        present = [n for n in neighbours if rows_by_taxid.get(n)]
        missing = [n for n in neighbours if not rows_by_taxid.get(n)]

        keys: list[str] = []
        sources: list[Path] = []
        expect: dict[str, tuple[int, int]] = {}
        dropped: list[tuple[str, str]] = []
        for owner in [tid] + present:
            for r in rows_by_taxid.get(owner, []):
                k = builder.assembly_key(Path(r["path"]), accession_by_path)
                if k in keys:
                    dropped.append((owner, k))
                    continue
                keys.append(k)
                sources.append(Path(r["path"]))
                n, b = expect.get(owner, (0, 0))
                expect[owner] = (n + int(r.get("n_contigs") or 0),
                                 b + int(r.get("genome_size") or 0))
        out.append({"taxid": tid, "name": (w.get("name") or "").strip(),
                    "keys": keys, "sources": sources, "expect": expect,
                    "missing_neighbours": missing, "dropped": dropped,
                    "has_target": bool(target_rows)})
    return out, notes


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
    ap.add_argument("--refs-dir", type=Path, default=None,
                    help="override [<btf-db>/packs/<pack>/refs]")
    ap.add_argument("--mm2-preset", default=None,
                    help="indexes expected per bundle, comma-separated. Default: whichever "
                         "presets the pack already holds, so a pack built with one preset "
                         "is not reported stale for lacking the other - while a bundle "
                         "missing an index its siblings have still is")
    ap.add_argument("--quick", action="store_true",
                    help="metadata only - do not read the bundles. Answers from the input "
                         "record and mtimes, which is UNDETERMINED for any bundle built "
                         "before the sidecar existed")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)

    try:
        builder = load_builder()
    except Exception as exc:                                   # noqa: BLE001
        print(f"UNDETERMINED: cannot load build_pack_refs.py to borrow its staleness "
              f"predicate ({exc}). This check deliberately has no copy of it.")
        return EXIT_UNDETERMINED

    if args.self_test:
        return self_test(builder)

    if args.refs_dir is None and args.btf_db is None:
        print("UNDETERMINED: no --refs-dir or --btf-db, and $BTF_DB is not set.\n"
              "  export BTF_DB=/path/to/project/bioThreatFinder, or pass --refs-dir DIR.")
        return EXIT_UNDETERMINED

    refs = args.refs_dir or (args.btf_db / "packs" / args.pack / "refs")
    comp_dir = refs / "competitive"
    manifest_path = refs / "manifest.tsv"
    watchlist_path = args.packs_dir / args.pack / "watchlist.resolved.tsv"

    for p in (watchlist_path, manifest_path):
        if not p.is_file():
            print(f"UNDETERMINED: no {p}")
            print("  Nothing can be said about bundles without both the watchlist that")
            print("  defines them and the manifest that records their inputs.")
            return EXIT_UNDETERMINED
    if not comp_dir.is_dir():
        print(f"UNDETERMINED: no competitive directory at {comp_dir}")
        print("  The pack's references may exist without its bundles; build them with")
        print("  scripts/build_pack_refs.py --pack " + args.pack)
        return EXIT_UNDETERMINED

    try:
        watchlist = [r for r in read_tsv(watchlist_path)
                     if (r.get("taxid") or "").strip().isdigit()]
        manifest = read_tsv(manifest_path)
    except OSError as exc:
        print(f"UNDETERMINED: cannot read the pack's tables ({exc})")
        return EXIT_UNDETERMINED

    if args.mm2_preset:
        presets = [p.strip() for p in args.mm2_preset.split(",") if p.strip()]
    else:
        # `<taxid>.<preset>.mmi`. Derived rather than assumed: `--mm2-preset` is a build
        # option, so hardcoding "sr,map-ont" here would report every pack built with one
        # of them stale for lacking the other, on a deployment that never wanted it.
        # Union over the directory, so a bundle missing an index that the other 89 have is
        # still stale. The limit: on a pack with ONE bundle the union is whatever that
        # bundle holds, and a preset that was never built cannot be told from one that
        # went missing - which is what `--mm2-preset` is for.
        presets = sorted({p.name.split(".")[1] for p in comp_dir.glob("*.*.mmi")
                          if p.name.count(".") == 2})
    if not presets:
        print(f"UNDETERMINED: no minimap2 index of any preset in {comp_dir}")
        print("  Either the pack was built with --skip-index, in which case stage 2 has")
        print("  nothing to align against, or this is not a bundle directory. Which of")
        print("  those it is cannot be told from here, so nothing is reported.")
        return EXIT_UNDETERMINED

    targets, notes = plan(builder, watchlist, manifest)
    if not targets:
        print(f"UNDETERMINED: {watchlist_path} resolved no targets")
        return EXIT_UNDETERMINED

    print(f"refs      : {refs}")
    print(f"manifest  : {manifest_path} ({len(manifest)} rows)")
    print(f"watchlist : {watchlist_path} ({len(targets)} targets)")
    print(f"mode      : {'metadata only (--quick)' if args.quick else 'reading every bundle'}"
          f"; indexes checked: {', '.join(presets)}\n")

    stale: list[tuple[str, str, str]] = []      # taxid, name, why
    wrong: list[tuple[str, str, str]] = []      # taxid, name, what is wrong inside
    no_record: list[str] = []
    unbuildable: list[tuple[str, str]] = []
    weakened: list[tuple[str, str, list[str]]] = []
    deduped: list[tuple[str, str, list[str]]] = []
    unreadable: list[tuple[str, str]] = []
    checked_deep = 0

    for t in targets:
        tid, name = t["taxid"], t["name"]
        comp = comp_dir / f"{tid}.fna.gz"
        if not t["has_target"]:
            unbuildable.append((tid, name))
            continue
        if t["missing_neighbours"]:
            weakened.append((tid, name, t["missing_neighbours"]))
        if t["dropped"]:
            deduped.append((tid, name, sorted({d[0] for d in t["dropped"]})))

        why = builder.bundle_staleness(comp, t["keys"], t["sources"])
        if why:
            stale.append((tid, name, why))
        elif builder.read_bundle_inputs(comp) is None:
            no_record.append(tid)

        for preset in presets:
            iwhy = builder.index_staleness(comp_dir / f"{tid}.{preset}.mmi", comp)
            if iwhy:
                stale.append((tid, name, f"{preset} index: {iwhy}"))

        if args.quick or not comp.is_file():
            continue
        try:
            got, unprefixed, collisions = bundle_prefix_stats(comp)
        except (OSError, EOFError, gzip.BadGzipFile) as exc:
            unreadable.append((tid, f"{comp.name}: {exc}"))
            continue
        checked_deep += 1
        for owner, (n_exp, b_exp) in sorted(t["expect"].items()):
            n_got, b_got = got.get(owner, (0, 0))
            if (n_got, b_got) != (n_exp, b_exp):
                wrong.append((tid, name,
                              f"taxid {owner}: bundle holds {n_got} contig(s)/{b_got} bp, "
                              f"the manifest says {n_exp}/{b_exp}"))
        for owner in sorted(set(got) - set(t["expect"]) - {"?"}):
            n_got, b_got = got[owner]
            wrong.append((tid, name, f"taxid {owner} is in the bundle ({n_got} "
                                     f"contig(s)/{b_got} bp) and in nothing that defines "
                                     f"it - not the target, not a current neighbour"))
        if unprefixed:
            wrong.append((tid, name, f"{got.get('?', (0, 0))[0]} header(s) with no taxid| "
                                     f"prefix, e.g. {unprefixed[0]!r} - no read over that "
                                     f"sequence can be attributed"))
        for a, b in sorted(collisions):
            wrong.append((tid, name, f"identical sequence under taxid {a} and taxid {b} - "
                                     f"every read over it gets MAPQ 0 and this taxon "
                                     f"cannot confirm"))

    if notes:
        print(f"--- {len(notes)} manifest row(s) point at a file that is gone ---")
        for n in notes[:10]:
            print(f"    {n}")
        if len(notes) > 10:
            print(f"    ... and {len(notes) - 10} more")
        print("    The builder would re-download these, so the bundles above them are not")
        print("    current whatever their timestamps say.")
        print()

    if stale:
        print(f"--- {len(stale)} stale artifact(s) ---")
        for tid, name, why in stale:
            print(f"    {tid:>8s}  {name[:30]:<30s} {why}")
        print()
    if wrong:
        print(f"--- {len(wrong)} bundle(s) disagree with the manifest, or are internally "
              f"wrong ---")
        for tid, name, what in wrong:
            print(f"    {tid:>8s}  {name[:30]:<30s} {what}")
        print()
    if unreadable:
        print(f"--- {len(unreadable)} bundle(s) could not be read ---")
        for tid, what in unreadable:
            print(f"    {tid:>8s}  {what}")
        print()
    if unbuildable:
        print(f"--- {len(unbuildable)} target(s) have no reference at all; stage 2 cannot "
              f"confirm them ---")
        for tid, name in unbuildable:
            print(f"    {tid:>8s}  {name}")
        print("    Not staleness: there is nothing to be stale. Reported because a bundle")
        print("    that was never built looks exactly like one that is current here.")
        print()
    if weakened:
        print(f"--- {len(weakened)} bundle(s) are missing a declared near neighbour ---")
        for tid, name, miss in weakened[:12]:
            print(f"    {tid:>8s}  {name[:30]:<30s} no genome for {', '.join(miss)}")
        if len(weakened) > 12:
            print(f"    ... and {len(weakened) - 12} more")
        print("    Competitive mapping is weakened, not stale: a relative that is absent")
        print("    from the bundle cannot take its own reads back, which is a")
        print("    false-POSITIVE risk on this target.")
        print()
    if deduped:
        print(f"--- {len(deduped)} bundle(s) lost a neighbour to de-duplication ---")
        for tid, name, ntids in deduped[:12]:
            print(f"    {tid:>8s}  {name[:30]:<30s} same assembly as {', '.join(ntids)}")
        if len(deduped) > 12:
            print(f"    ... and {len(deduped) - 12} more")
        print("    Correct behaviour - two copies would give every read MAPQ 0 - but those")
        print("    pairs cannot be discriminated until a distinct assembly is curated.")
        print()

    n = len(targets)
    if not args.quick:
        print(f"read {checked_deep} of {n} bundle(s) and compared every taxid prefix "
              f"against the manifest.")
    elif no_record:
        print(f"{len(no_record)} of {n} bundle(s) have no input record and were checked by "
              f"mtime alone.")

    # A vanished reference counts. The builder would re-download it, the fresh file would
    # then be newer than the bundle, and the bundle would be rebuilt - so a pack with a
    # hole in it is not current, however settled its timestamps look.
    if stale or wrong or unreadable or notes:
        print(f"\nverdict: STALE. "
              f"{len(stale) + len(wrong) + len(unreadable) + len(notes)} problem(s) "
              f"above. Rebuild with:")
        print(f"         scripts/build_pack_refs.py --pack {args.pack}")
        print("         Existing, current bundles are skipped; only what is named above is")
        print("         rebuilt, and nothing is re-downloaded that the manifest's md5s")
        print("         still match.")
        return EXIT_GAP

    if args.quick and no_record:
        print("\nverdict: UNDETERMINED. Nothing looks stale, but the bundles without an")
        print("         input record were compared by timestamp only, which cannot see a")
        print("         restored or re-downloaded genome. Drop --quick to compare their")
        print("         contents against the manifest instead.")
        return EXIT_UNDETERMINED

    print(f"\nverdict: CURRENT. All {n} bundle(s) and {n * len(presets)} index(es) match "
          f"the manifest.")
    print("         This says the bundles agree with the manifest, not that the manifest")
    print("         agrees with refs.lock.tsv - that needs build_pack_refs.py --from-lock,")
    print("         which must reach NCBI.")
    return EXIT_OK


def self_test(builder) -> int:
    """Build a tiny pack in a temporary directory and mutate it.

    A staleness check that has never been shown to FAIL is the thing it exists to catch:
    it would report `CURRENT` over a pack with nothing in it. Each case below changes one
    thing about a bundle that is otherwise current, and the check must notice that one
    thing - including the two cases metadata cannot see, which is the whole argument for
    reading the bundles.
    """
    import contextlib
    import io
    import tempfile
    import textwrap

    def write_fasta(p: Path, records: list[tuple[str, str]]) -> tuple[int, int]:
        with gzip.open(p, "wt") as fh:
            for hdr, seq in records:
                fh.write(f">{hdr}\n{seq}\n")
        return sum(len(s) for _, s in records), len(records)

    results: list[tuple[str, bool]] = []
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        packs, refs = root / "packs", root / "db" / "packs" / "human" / "refs"
        comp = refs / "competitive"
        comp.mkdir(parents=True)
        (packs / "human").mkdir(parents=True)

        (packs / "human" / "watchlist.resolved.tsv").write_text(textwrap.dedent("""\
            taxid\tname\tgroup\tnear_neighbour_taxids
            1000\tTarget one\tbacteria\t2000
            """))
        tgt = refs / "1000_Target_one.fna.gz"
        nbr = refs / "2000_Neighbour.fna.gz"
        t_bases, t_contigs = write_fasta(tgt, [("c1", "ACGTACGTAC"), ("c2", "TTTTGGGG")])
        n_bases, n_contigs = write_fasta(nbr, [("n1", "CCCCAAAAGG")])
        rows = [("1000", "target", "GCF_000000001.1", t_bases, t_contigs, tgt),
                ("2000", "near_neighbour", "GCF_000000002.1", n_bases, n_contigs, nbr)]
        bundle = comp / "1000.fna.gz"

        def write_manifest(extra: list | None = None) -> None:
            lines = ["\t".join(["taxid", "name", "group", "role", "accession",
                                "assembly_level", "genome_size", "n_contigs", "path",
                                "md5", "downloaded_utc"])]
            for tid, role, acc, bases, contigs, path in (rows + (extra or [])):
                lines.append("\t".join([tid, f"taxon {tid}", "bacteria", role, acc,
                                        "complete", str(bases), str(contigs), str(path),
                                        "", ""]))
            (refs / "manifest.tsv").write_text("\n".join(lines) + "\n")

        def build_bundle(members: list[tuple[str, Path, str]], sidecar: bool = True) -> None:
            """Write a bundle the way the builder does, then stamp and index it."""
            with gzip.open(bundle, "wt") as out:
                for tid, src, _ in members:
                    builder.append_with_taxid_prefix(src, int(tid), out)
            p = builder.bundle_inputs_path(bundle)
            p.unlink(missing_ok=True)
            if sidecar:
                builder.write_bundle_inputs(bundle, {acc for _, _, acc in members})
            for preset in ("sr", "map-ont"):
                (comp / f"1000.{preset}.mmi").write_bytes(b"x")

        both = [("1000", tgt, "GCF_000000001.1"), ("2000", nbr, "GCF_000000002.1")]

        def run(*extra: str) -> int:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                return main(["--pack", "human", "--packs-dir", str(packs),
                             "--refs-dir", str(refs), *extra])

        def case(label: str, want: int, *extra: str) -> None:
            got = run(*extra)
            results.append((f"{label} -> exit {got} (want {want})", got == want))

        write_manifest()
        build_bundle(both)
        case("a current pack is CURRENT", EXIT_OK)
        # The baseline must pass in --quick too, or every case below exercises the deep
        # pass alone and the sidecar path is never tested at all.
        case("a current pack is CURRENT in --quick too", EXIT_OK, "--quick")

        # Named explicitly: this pack has one bundle, so the preset set derived from the
        # directory would become whatever that one bundle still holds. On a real pack the
        # union over 90 bundles makes a single missing index stale on its own.
        (comp / "1000.sr.mmi").unlink()
        case("a missing index is caught", EXIT_GAP, "--mm2-preset", "sr,map-ont")
        # No index of ANY preset: either --skip-index was used or this is not a bundle
        # directory, and the two cannot be told apart from here. Undetermined, not stale.
        (comp / "1000.map-ont.mmi").unlink()
        case("...and a pack with no index of any preset is UNDETERMINED",
             EXIT_UNDETERMINED)
        build_bundle(both)

        build_bundle([("1000", tgt, "GCF_000000001.1")])
        case("a neighbour left out of the bundle is caught", EXIT_GAP)

        build_bundle(both)
        os.utime(bundle, (1, 1))
        case("a bundle older than its inputs is caught", EXIT_GAP)

        # THE CASE THAT MOTIVATES READING THE BUNDLES. The manifest gains a second
        # assembly for the target - D42's shape exactly: the manifest claims two
        # references while the aligner uses one. With a sidecar this is exact; without
        # one, and with the new file's mtime older than the bundle, metadata is blind.
        extra_ref = refs / "1000_Target_one.GCF_000000003.1.fna.gz"
        e_bases, e_contigs = write_fasta(extra_ref, [("x1", "GGGGCCCCTTTT")])
        extra_row = [("1000", "target", "GCF_000000003.1", e_bases, e_contigs, extra_ref)]

        build_bundle(both)
        write_manifest(extra_row)
        case("a reference the manifest added is caught by the sidecar", EXIT_GAP, "--quick")

        build_bundle(both, sidecar=False)
        os.utime(extra_ref, (1, 1))          # older than the bundle: no mtime signal
        case("...and by reading the bundle when there is no sidecar", EXIT_GAP)
        got_quick = run("--quick")
        results.append((
            f"...while --quick says UNDETERMINED rather than CURRENT -> exit {got_quick} "
            f"(want {EXIT_UNDETERMINED})", got_quick == EXIT_UNDETERMINED))
        write_manifest()

        # Identical sequence under two prefixes: the MAPQ-0 inversion. Written by hand,
        # because the builder de-duplicates by accession and cannot produce it - which is
        # also why no bundle built before that de-duplication has ever been tested for it.
        build_bundle([("1000", tgt, "GCF_000000001.1"), ("2000", tgt, "GCF_000000002.1")])
        case("identical sequence under two prefixes is caught", EXIT_GAP)

        build_bundle(both)
        with gzip.open(bundle, "at") as out:
            out.write(">stray contig\nACGT\n")
        case("a header with no taxid| prefix is caught", EXIT_GAP)

        build_bundle(both)
        case("the pack is CURRENT again once each mutation is undone", EXIT_OK)

        tgt_hidden = tgt.with_suffix(".hidden")
        tgt.rename(tgt_hidden)
        case("a manifest row whose file is gone is caught", EXIT_GAP)
        tgt_hidden.rename(tgt)

        (refs / "manifest.tsv").unlink()
        case("no manifest is UNDETERMINED, not stale", EXIT_UNDETERMINED)
        write_manifest()
        case("and CURRENT once more afterwards", EXIT_OK)

    print("\n--- self-test ---")
    for label, ok in results:
        print(f"  {'ok  ' if ok else 'FAIL'} {label}")
    bad = [l for l, ok in results if not ok]
    print("\nSELF-TEST " + ("PASSED" if not bad else f"FAILED ({len(bad)})"))
    return EXIT_OK if not bad else EXIT_GAP


if __name__ == "__main__":
    sys.exit(main())
