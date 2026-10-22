#!/usr/bin/env python3
"""Can the Kraken2 index report every taxon on the watchlist?

D47 left the Kraken2 index without a check and said why the bundles' freshness predicate
does not transfer: its inputs are NCBI library snapshots, not `packs/`. The honest question
for it is a different one, and this is that question - for every watchlist taxon, does the
index hold a node that could carry its reads?

A taxon absent from the index's taxonomy can never appear in a Kraken2 report. Nothing
errors, the arm runs, `screen_summary.tsv` records it as `ok`, and the row simply never
exists: the silent subtraction D34 found in the masking BEDs and D46 in the sketch.
`pack.yaml` makes it sharper than a general worry - it sets
`eukaryote_strategy: kraken2_then_align` and says to rely on Kraken2 PlusPF for parasites
and fungi, so for those taxa this is the arm the pack nominates.

READS taxo.k2d DIRECTLY. The head node has no kraken2 (it lives in btf-kraken2 and is used
inside the process), which is D46's constraint again - but unlike a `.syldb` this file has a
documented layout and needs no tool:

    8 bytes  "K2TAXDAT"
    3 x uint64 LE: node_count, name_data_len, rank_data_len
    node_count records of 7 x uint64: parent_id, first_child, child_count,
                                      name_offset, rank_offset, EXTERNAL_ID, godparent_id
    then name_data_len bytes, then rank_data_len bytes

The arithmetic is self-checking: 32 + node_count*56 + name_data_len + rank_data_len must
equal the file size exactly. If it does not, this refuses to report rather than guess - a
wrong stride reads garbage taxids and reports them confidently.

TWO THINGS IT DOES NOT CLAIM.

The pruned taxonomy is NOT the pinned taxdump. kraken2-build keeps only nodes spanning
sequences it loaded, so a correct index holds far fewer nodes than the snapshot. Demanding
equality would call every good index broken.

A node present is not SEQUENCE present: a node can exist because an ancestor or descendant
carried a genome. So an absent taxid PROVES the arm cannot report that taxon, while a
present one is necessary but not sufficient. The verdict is worded accordingly.

    scripts/check_index_coverage.py --pack human --btf-db /media/8tb/bioThreatFinder
      exit 0  every watchlist taxon has a node
      exit 3  at least one cannot be reported by this arm
      exit 1  could not tell (no index, unreadable, layout does not add up)
"""
from __future__ import annotations

import argparse
import csv
import struct
import sys
from pathlib import Path

K2_MAGIC = b"K2TAXDAT"
K2_HEADER_LEN = 32
K2_NODE_STRIDE = 56        # 7 x uint64
K2_EXTERNAL_ID_OFFSET = 40  # field index 5

EXIT_OK = 0
EXIT_UNDETERMINED = 1
EXIT_GAP = 3


def read_index_taxids(taxo: Path) -> set[int]:
    """Every NCBI taxid the index's pruned taxonomy holds.

    Raises ValueError if the file is not a taxo.k2d or its layout does not add up, so a
    wrong assumption cannot be reported as a confident set of taxids.
    """
    size = taxo.stat().st_size
    with taxo.open("rb") as fh:
        head = fh.read(K2_HEADER_LEN)
        if len(head) < K2_HEADER_LEN or head[:8] != K2_MAGIC:
            raise ValueError(f"{taxo.name} does not begin with {K2_MAGIC!r}")
        n_nodes, name_len, rank_len = struct.unpack("<QQQ", head[8:32])
        expect = K2_HEADER_LEN + n_nodes * K2_NODE_STRIDE + name_len + rank_len
        if expect != size:
            raise ValueError(
                f"{taxo.name} layout does not add up: {K2_HEADER_LEN} + "
                f"{n_nodes}*{K2_NODE_STRIDE} + {name_len} + {rank_len} = {expect}, "
                f"but the file is {size} bytes. Refusing to read taxids on a wrong stride.")
        buf = fh.read(n_nodes * K2_NODE_STRIDE)
    if len(buf) != n_nodes * K2_NODE_STRIDE:
        raise ValueError(f"{taxo.name} truncated while reading nodes")
    return {struct.unpack_from("<Q", buf, i * K2_NODE_STRIDE + K2_EXTERNAL_ID_OFFSET)[0]
            for i in range(n_nodes)}


def read_watchlist(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(
            (ln for ln in fh if not ln.startswith("#")), delimiter="\t"))
    return [r for r in rows if (r.get("taxid") or "").strip().isdigit()]


def read_sketched(sidecar: Path) -> set[str]:
    """Basenames in the sylph sketch's record, or empty if there is none."""
    if not sidecar.exists():
        return set()
    return {ln.strip() for ln in sidecar.read_text(encoding="utf-8").splitlines()
            if ln.strip() and not ln.startswith("#")}


def refs_by_taxid(manifest: Path) -> dict[int, list[dict]]:
    out: dict[int, list[dict]] = {}
    if not manifest.exists():
        return out
    with manifest.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = (r.get("taxid") or "").strip()
            if t.isdigit():
                out.setdefault(int(t), []).append(r)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pack", default="human")
    ap.add_argument("--packs-dir", type=Path,
                    default=Path(__file__).resolve().parent.parent / "packs")
    ap.add_argument("--btf-db", type=Path, default=Path("/media/8tb/bioThreatFinder"))
    ap.add_argument("--kraken2-db", type=Path, default=None,
                    help="the index directory holding taxo.k2d "
                         "[<btf-db>/kraken2/pluspf]")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()

    watchlist = args.packs_dir / args.pack / "watchlist.resolved.tsv"
    if not watchlist.exists():
        print(f"UNDETERMINED: no resolved watchlist at {watchlist}")
        return EXIT_UNDETERMINED
    k2_dir = args.kraken2_db or (args.btf_db / "kraken2" / "pluspf")
    taxo = k2_dir / "taxo.k2d"
    if not taxo.exists():
        print(f"UNDETERMINED: no taxo.k2d under {k2_dir}")
        return EXIT_UNDETERMINED

    try:
        in_index = read_index_taxids(taxo)
    except (OSError, ValueError) as exc:
        print(f"UNDETERMINED: {exc}")
        return EXIT_UNDETERMINED

    rows = read_watchlist(watchlist)
    if not rows:
        print(f"UNDETERMINED: no taxid rows in {watchlist}")
        return EXIT_UNDETERMINED

    missing = [r for r in rows if int(r["taxid"]) not in in_index]
    print(f"index   : {taxo} ({len(in_index)} distinct taxids in its pruned taxonomy)")
    print(f"watchlist: {watchlist} ({len(rows)} taxa)")

    if not missing:
        print(f"\nverdict: OK. All {len(rows)} watchlist taxa have a node in this index. "
              f"A node is necessary but not sufficient - it may exist because a relative "
              f"carried the genome - so this bounds what the arm CANNOT report, not what "
              f"it will.")
        return EXIT_OK

    # What else could raise these, so the consequence is bounded rather than alarming.
    # sylph profiles THIS PACK's refs and is on by default, so a taxon with a sketched
    # reference is still reachable by stage 1 through that arm.
    sketched = read_sketched(args.btf_db / "sylph" /
                             f"{args.pack}_pack_refs.syldb.genomes.txt")
    by_taxid = refs_by_taxid(args.btf_db / "packs" / args.pack / "refs" / "manifest.tsv")

    print(f"\n{len(missing)} of {len(rows)} watchlist taxa have NO node in this index "
          f"and can never appear in its report:\n")
    print(f"  {'taxid':>8}  {'name':<34} {'category':<17} other stage-1 cover")
    uncovered = []
    for r in sorted(missing, key=lambda r: r.get("category", "")):
        t = int(r["taxid"])
        with_seq = [x for x in by_taxid.get(t, []) if (x.get("path") or "").strip()]
        n_sk = sum(1 for x in with_seq if x["path"].rsplit("/", 1)[-1] in sketched)
        cover = f"sylph ({n_sk} ref)" if n_sk else "NONE"
        if not n_sk:
            uncovered.append(r)
        print(f"  {t:>8}  {r.get('name',''):<34} {r.get('category',''):<17} {cover}")

    print()
    if uncovered:
        print(f"verdict: GAP. {len(uncovered)} of those have NO other stage-1 arm either, "
              f"so they are invisible to stage 1 entirely:")
        for r in uncovered:
            print(f"           {r['taxid']} {r.get('name','')}")
    else:
        print(f"verdict: GAP, bounded. Every one of the {len(missing)} is sketched for "
              f"sylph, which is on by default and profiles this pack, so stage 1 can still "
              f"raise them through that arm - one arm is blind, not the stage. "
              f"`combining_rule: any_classifier_or` is what makes that true.")
    print("         This is a property of the index's library selection, not staleness: "
          "rebuilding against the same libraries will not change it.")
    return EXIT_GAP


def self_test() -> int:
    """Offline checks on the two things that could report confident nonsense."""
    import tempfile
    failures: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        print(f"  {'ok  ' if ok else 'FAIL'} {name}{'' if ok else ': ' + detail}")
        if not ok:
            failures.append(name)

    print("check_index_coverage.py --self-test")

    def make_taxo(path: Path, taxids: list[int], *, name_len: int = 4,
                  rank_len: int = 2, magic: bytes = K2_MAGIC,
                  pad: int = 0) -> None:
        body = b"".join(
            struct.pack("<QQQQQQQ", 1, 0, 0, 0, 0, t, 0) for t in taxids)
        with path.open("wb") as fh:
            fh.write(magic)
            fh.write(struct.pack("<QQQ", len(taxids), name_len, rank_len))
            fh.write(body)
            fh.write(b"\x00" * (name_len + rank_len + pad))

    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        good = d / "taxo.k2d"
        make_taxo(good, [1, 9606, 1392, 632])
        try:
            got = read_index_taxids(good)
            check("a well-formed taxonomy round-trips its external ids",
                  got == {1, 9606, 1392, 632}, str(sorted(got)))
        except ValueError as exc:
            check("a well-formed taxonomy round-trips its external ids", False, str(exc))

        # THE DECISIVE ONE. A wrong stride still "reads" - it just reads the wrong 8 bytes
        # of every record and reports them as taxids. The size arithmetic is the only thing
        # standing between that and a confident, wrong answer, so a file whose arithmetic
        # does not add up must REFUSE rather than return a set.
        bad = d / "bad.k2d"
        make_taxo(bad, [1, 2, 3], pad=17)
        raised = False
        try:
            read_index_taxids(bad)
        except ValueError:
            raised = True
        check("a layout that does not add up is refused, not read", raised,
              "a wrong stride would report garbage taxids confidently")

        notk2 = d / "other.k2d"
        make_taxo(notk2, [1], magic=b"NOTK2DAT")
        raised = False
        try:
            read_index_taxids(notk2)
        except ValueError:
            raised = True
        check("a file that is not a taxo.k2d is refused", raised)

        # A pruned taxonomy is SMALLER than the snapshot by design; nothing here may
        # compare the two. Guard the constant rather than the comparison.
        check("the node stride is 7 x uint64", K2_NODE_STRIDE == 7 * 8)
        check("external_id is field index 5", K2_EXTERNAL_ID_OFFSET == 5 * 8)

    print(f"\n{'SELF-TEST PASSED' if not failures else str(len(failures)) + ' FAILED'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
