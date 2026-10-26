#!/usr/bin/env python3
"""Convert stage-1 classifier reports into Pavian-readable Kraken reports.

WHAT PAVIAN ACTUALLY ACCEPTS
----------------------------
Verified against pavian/R/datainput-read_report.R, not assumed. Pavian has TWO
readers and they do not accept the same things:

  read_report()  - the Shiny app's entry point. Headerless files are dispatched
                   on the TAB COUNT OF THE FIRST LINE:
                       ntabs == 5 -> percentage, cladeReads, taxonReads,
                                     taxRank, taxID, name
                       ntabs == 7 -> percentage, cladeReads, taxonReads,
                                     n_unique_kmers, n_kmers, taxRank, taxID, name
                   Any other tab count leaves col_names undefined, the read
                   throws, and the FILE IS REJECTED WHOLESALE.
                   So Kraken2 --report-minimizer-data (8 fields) IS supported
                   here, as is KrakenUniq via its own header.

  read_report2() - headerless files are read with EXACTLY six column names.
                   The 8-field minimizer report does NOT work here.

Consequences this script is built around:

  * '#' IS NOT A COMMENT. read_report decides on
        has_header <- grepl("^[a-zA-Z#%\"]", first.line)
    so a leading '#' makes Pavian treat that line as the HEADER ROW and the
    parse collapses. Nothing written here starts a report file with '#'.
    Provenance goes to pavian_manifest.tsv and README.txt beside the reports.

  * 6 fields is the only layout BOTH readers accept, so that is what is written.
    The k-mer columns are therefore dropped on conversion; the untouched 8-field
    Kraken2 report stays published under kraken2/ and can be loaded directly
    into the app when the minimizer evidence is wanted.

  * Indentation of the name field is DATA, not cosmetics:
        report$depth <- nchar(gsub("\\S.*","",report$name))/2
    Tree depth is recovered from it. Passthrough preserves it byte-for-byte;
    synthesis rebuilds it at two spaces per level.

WHAT ACTUALLY NEEDED CONVERTING
-------------------------------
Kraken2, Bracken and KrakenUniq are already loadable by the app as published.
They are re-emitted anyway so that one directory holds every arm under one
naming scheme, which is what makes Pavian's sample-comparison view line up.
The two that genuinely do not work without this script are:

  metabuli  Its report opens with a '#'-prefixed header row
            (#clade_proportion, clade_count, taxon_count, rank, taxID, name).
            Per the rule above that '#' makes Pavian treat line 1 as a HEADER,
            and since those column names are not in its harmonisation table it
            finds no cladeReads/name and rejects the file. Dropping the header
            row fixes it. (This was first attributed here to the trailing tab
            of upstream issue #194; the two installed reports carry 0 trailing
            tabs and 6 fields on every line, so #194 is not what bites. The
            trailing-field strip below is kept anyway - it is correct, it costs
            nothing, and #194 is still open.)

  sylph     Not a Kraken report at all - no tree, no read counts. Synthesised
            from the NCBI taxdump; see synth_from_sylph() for what is real and
            what is estimated.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

RANK_CODE = {
    "superkingdom": "D", "domain": "D", "realm": "D",
    "kingdom": "K", "phylum": "P", "class": "C", "order": "O",
    "family": "F", "genus": "G", "species": "S",
    "subspecies": "S1", "strain": "S1", "serotype": "S1", "serogroup": "S1",
    "biotype": "S1", "genotype": "S1", "forma specialis": "S1", "pathogroup": "S1",
    "subgenus": "G1", "species group": "G1", "species subgroup": "G2",
    "subfamily": "F1", "tribe": "F2", "suborder": "O1", "subclass": "C1",
    "subphylum": "P1", "clade": "-", "no rank": "-",
    # Metabuli's spellings for the nodes above the domains.
    "root": "R", "cellular root": "R1", "acellular root": "R1",
}

# KrakenUniq header spellings -> canonical field. Same map as parse_screen.py.
KU_COLUMNS = {
    "%": "percentage", "reads": "clade_reads", "taxreads": "taxon_reads",
    "taxid": "taxid", "rank": "rank", "taxname": "name",
}


def normalise_rank(rank: str, taxid: str) -> str:
    """Give every arm the single-letter rank codes Pavian actually filters on.

    Kraken2 and Bracken already emit codes (U, R, R1, D, ... S, S1) and are
    passed through untouched. Metabuli and KrakenUniq emit FULL WORDS ("domain",
    "cellular root", "no rank"), and Pavian maps only a fixed set of those,
    rendering everything else as "-" - so Bacteria would sit at "-" for those
    arms while Kraken2 put it at "D", and the arms would not line up in the
    rank filters or the comparison view.
    """
    if taxid == "0":
        return "U"
    if taxid == "1":
        return "R"
    r = rank.strip()
    if re.fullmatch(r"[A-Z][0-9]*", r):
        return r
    return RANK_CODE.get(r.lower(), "-")


def as_int(v, default=0) -> int:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def as_float(v, default=0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def usable(path: Path | None) -> bool:
    return bool(path) and path.exists() and path.name != "NO_FILE" and path.stat().st_size > 0


def indent_of(raw_name: str) -> str:
    """Leading whitespace of a report name field, preserved verbatim."""
    return raw_name[: len(raw_name) - len(raw_name.lstrip(" "))]


# --------------------------------------------------------------------------
# Passthrough converters
# --------------------------------------------------------------------------
def looks_like_abundance_table(f: list[str]) -> bool:
    """Is this row a Bracken ABUNDANCE table line rather than a report line?

    `name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads,
    new_est_reads, fraction_total_reads` - seven fields, the taxid SECOND. The
    BRACKEN process writes one of these beside the `-w` report it writes for
    this script, so the wrong file is one argument away in the same directory.
    Read from the right it parses WITHOUT ERROR into nonsense: f[-2] is
    `new_est_reads`, all digits, so it becomes the taxid; f[-3] is
    `added_reads`, so it becomes the rank; f[-1] is the fraction, so it becomes
    the name. The result is a Pavian report of fabricated taxids with zero
    unparseable lines and exit 0 - this project's recurring shape, a degraded
    result indistinguishable from a good one.
    """
    return len(f) == 7 and f[1].strip().isdigit() and f[2].strip() in {
        "D", "K", "P", "C", "O", "F", "G", "S", "S1", "S2", "R", "R1", "U",
    }


def rows_from_kraken_style(path: Path) -> tuple[list[tuple], int]:
    """Kraken2 (6 or 8 col), Bracken, Metabuli -> canonical 6-field rows.

    Indexed from the RIGHT, which is what makes one reader cover all three
    layouts: they differ on the left (Kraken2's two minimizer columns) but keep
    ... rank, taxid, name stable at the end. Lifted from parse_screen.py so the
    two never disagree about what a column means.

    Right-indexing is why the shape has to be checked on the LEFT as well: see
    looks_like_abundance_table(), and the percentage check below. Every one of
    these layouts leads with a clade percentage, so field 0 is a number; the
    file this reader must not silently accept leads with a taxon name.
    """
    rows, bad = [], 0
    abundance_shaped = 0
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "%")):
                continue
            f = line.rstrip("\n").split("\t")
            # Belt and braces for upstream #194 (a trailing tab on unclassified
            # rows), which would shift every right-indexed field by one. The
            # installed Metabuli does not do this, but the issue is still open.
            while f and f[-1].strip() == "":
                f.pop()
            if len(f) < 6:
                bad += 1
                continue
            if looks_like_abundance_table(f):
                abundance_shaped += 1
                bad += 1
                continue
            # A report row's first field is the clade percentage. Refusing a
            # non-numeric one costs nothing on a real report and is what makes
            # a mis-passed file report zero usable rows instead of nonsense.
            try:
                float(f[0].strip() or "x")
            except ValueError:
                bad += 1
                continue
            taxid, rank, name = f[-2].strip(), f[-3].strip(), f[-1]
            if not taxid.isdigit():
                bad += 1
                continue
            rows.append((f[0].strip(), as_int(f[1]), as_int(f[2]),
                         normalise_rank(rank, taxid), taxid, name))
    if abundance_shaped:
        sys.stderr.write(
            f"[pavian] WARNING: {path.name}: {abundance_shaped} row(s) have the "
            f"shape of a Bracken ABUNDANCE table (name, taxonomy_id, "
            f"taxonomy_lvl, ...), not a Kraken-style report. Pass the file "
            f"written by `bracken -w`, not the abundance table beside it; these "
            f"rows were rejected rather than read from the right into "
            f"fabricated taxids.\n")
    return rows, bad


def rows_from_krakenuniq(path: Path) -> tuple[list[tuple], int]:
    """KrakenUniq -> canonical 6-field rows, header-aware.

    Column order is ... taxID, rank, taxName, i.e. NOT Kraken's, so the fields
    are located by header name instead of by position.
    """
    rows, bad, idx = [], 0, None
    with path.open() as fh:
        for line in fh:
            raw = line.rstrip("\n")
            if not raw.strip():
                continue
            fields = raw.split("\t")
            low = [c.strip().lower().lstrip("#").strip() for c in fields]
            if idx is None and any(c in ("kmers", "dup", "cov") for c in low):
                idx = {KU_COLUMNS[c]: i for i, c in enumerate(low) if c in KU_COLUMNS}
                continue
            if raw.startswith("#") or idx is None:
                continue
            if "taxid" not in idx or len(fields) <= max(idx.values()):
                bad += 1
                continue
            taxid = fields[idx["taxid"]].strip()
            if not taxid.isdigit():
                bad += 1
                continue
            rows.append((
                fields[idx["percentage"]].strip() if "percentage" in idx else "0.00",
                as_int(fields[idx["clade_reads"]]) if "clade_reads" in idx else 0,
                as_int(fields[idx["taxon_reads"]]) if "taxon_reads" in idx else 0,
                normalise_rank(fields[idx["rank"]] if "rank" in idx else "-", taxid),
                taxid,
                fields[idx["name"]] if "name" in idx else "",
            ))
    return rows, bad


# --------------------------------------------------------------------------
# Sylph synthesis
# --------------------------------------------------------------------------
def load_taxdump(taxdump: Path, wanted: set[str]) -> tuple[dict, dict, dict, dict]:
    """parent / rank for every taxid, scientific names for the needed closure,
    and the old->new remap for any wanted taxid that has since been merged.

    merged.dmp matters because pack references are NAMED with the taxid that was
    current when the pack was built. If NCBI later merges that taxid, the name on
    disk still carries the old one, nodes.dmp no longer has it, and without this
    the organism silently vanishes from the sylph view. resolve_watchlist.py
    already consults merged.dmp for the same reason; omitting it here was an
    inconsistency, not a decision.
    """
    parent: dict[str, str] = {}
    rank: dict[str, str] = {}
    with (taxdump / "nodes.dmp").open() as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) >= 3:
                parent[f[0]], rank[f[0]] = f[1], f[2]

    remap: dict[str, str] = {}
    missing = {t for t in wanted if t not in parent}
    merged = taxdump / "merged.dmp"
    if missing and merged.exists():
        with merged.open() as fh:
            for line in fh:
                f = [x.strip() for x in line.split("|")]
                if len(f) >= 2 and f[0] in missing and f[1] in parent:
                    remap[f[0]] = f[1]
    wanted = {remap.get(t, t) for t in wanted}

    needed: set[str] = set()
    for t in wanted:
        cur, guard = t, 0
        while cur and cur in parent and guard < 200:
            needed.add(cur)
            if cur == "1":
                break
            cur = parent[cur]
            guard += 1
    needed.add("1")

    names: dict[str, str] = {}
    with (taxdump / "names.dmp").open() as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) >= 4 and f[3] == "scientific name" and f[0] in needed:
                names[f[0]] = f[1]
    return parent, rank, names, remap


def taxid_from_genome_file(genome_file: str) -> str | None:
    """Pack references are named '<taxid>_<Name>...', so take the leading
    numeric token of the basename.

    Deliberately STRICTER than parse_screen.parse_sylph(), which accepts any
    digit token anywhere in the path. That looseness is defensible for raising a
    candidate (more evidence is better) but not for building a taxonomy tree,
    where one wrong taxid silently invents a lineage.
    """
    stem = Path(genome_file).name
    head = stem.split("_", 1)[0]
    return head if head.isdigit() and len(head) >= 2 else None


def synth_from_sylph(path: Path, taxdump: Path, total_reads: int) -> tuple[list[tuple], list[str]]:
    """Build a Kraken-shaped tree from a sylph profile.

    HONESTY NOTE, which is also written into every output file: sylph reports
    ANI/containment and an ABUNDANCE, not read assignments. The counts below are
    abundance x total_reads, i.e. ESTIMATES back-derived so that Pavian (which
    can only draw read counts) has something to draw. They are not reads that
    were assigned to these taxa, and they must not be quoted as such.
    """
    notes: list[str] = []
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return [], ["sylph profile has no rows (no genome passed its ANI threshold)"]

    leaf: dict[str, float] = {}
    unnamed = 0
    for r in rows:
        gf = r.get("Genome_file") or ""
        taxid = taxid_from_genome_file(gf)
        if not taxid:
            # COUNT, do not list. A sylph database built with --no-check-names
            # (e.g. a prebuilt GTDB sketch) fails every row, and listing them put
            # thousands of entries in a single manifest cell.
            unnamed += 1
            continue
        ab = as_float(r.get("Sequence_abundance"), -1.0)
        if ab < 0:
            ab = as_float(r.get("Taxonomic_abundance"), 0.0)
        leaf[taxid] = leaf.get(taxid, 0.0) + ab

    if unnamed:
        notes.append(f"{unnamed} sylph row(s) carried no taxid in the genome "
                     f"file name and were skipped (expected '<taxid>_<Name>'; a "
                     f"--no-check-names sylph database fails every row)")
    if not leaf:
        return [], notes or ["no usable sylph rows"]

    parent, rank, names, remap = load_taxdump(taxdump, set(leaf))
    if remap:
        merged_leaf: dict[str, float] = {}
        for taxid, ab in leaf.items():
            key = remap.get(taxid, taxid)
            merged_leaf[key] = merged_leaf.get(key, 0.0) + ab
        leaf = merged_leaf
        notes.append("resolved via merged.dmp: "
                     + ", ".join(f"{o}->{n}" for o, n in sorted(remap.items())))

    # Accumulate clade counts up every lineage.
    clade: dict[str, int] = {}
    direct: dict[str, int] = {}
    children: dict[str, set[str]] = {}
    for taxid, ab in leaf.items():
        reads = round(ab / 100.0 * total_reads)
        if taxid not in parent:
            notes.append(f"taxid {taxid} is not in the taxdump - skipped")
            continue
        direct[taxid] = direct.get(taxid, 0) + reads
        cur, guard = taxid, 0
        while cur and guard < 200:
            clade[cur] = clade.get(cur, 0) + reads
            if cur == "1":
                break
            par = parent.get(cur)
            if not par or par == cur:
                break
            children.setdefault(par, set()).add(cur)
            cur = par
            guard += 1

    if not clade:
        return [], notes
    # walk() indexes clade[] unguarded, so root MUST be present. It is absent only
    # if every lineage broke before taxid 1 (a parent missing from nodes.dmp, a
    # self-parent below root, or the guard tripping). Returning nothing beats a
    # KeyError: this step is presentational, and nextflow.config sets no
    # errorStrategy, so an exception here would kill the run before stage 2.
    if "1" not in clade:
        notes.append("no lineage reached root (taxid 1) - nothing emitted; the "
                     "taxdump may be incomplete or inconsistent with the pack")
        return [], notes

    total_clade = clade.get("1", 0)
    out: list[tuple] = []

    unassigned = max(total_reads - total_clade, 0)
    if unassigned:
        out.append((f"{unassigned / total_reads * 100:.2f}", unassigned, unassigned,
                    "U", "0", "unclassified"))

    def walk(taxid: str, depth: int) -> None:
        # Per-leaf round() can push the sum a few reads past total_reads, so
        # clamp rather than print 100.01%.
        pct = min(clade[taxid] / total_reads * 100, 100.0) if total_reads else 0.0
        code = "R" if taxid == "1" else RANK_CODE.get(rank.get(taxid, "no rank"), "-")
        name = names.get(taxid, f"taxid {taxid}")
        out.append((f"{pct:.2f}", clade[taxid], direct.get(taxid, 0),
                    code, taxid, "  " * depth + name))
        for kid in sorted(children.get(taxid, ()), key=lambda k: -clade.get(k, 0)):
            walk(kid, depth + 1)

    walk("1", 0)
    return out, notes


# --------------------------------------------------------------------------
def write_report(dest: Path, rows: list[tuple]) -> None:
    """Six tab-separated fields and nothing else.

    No header, no comment line, no trailing blank: the first line must carry
    exactly five tabs or Pavian's read_report() discards the whole file.
    """
    with dest.open("w") as fh:
        for pct, cladeR, taxonR, rank, taxid, name in rows:
            fh.write(f"{pct}\t{cladeR}\t{taxonR}\t{rank}\t{taxid}\t{name}\n")


def total_reads_from(rows: list[tuple]) -> int:
    """Total reads in a converted report: clade count at root plus unclassified."""
    root = next((r[1] for r in rows if r[4] == "1"), 0)
    uncl = next((r[1] for r in rows if r[4] == "0"), 0)
    return root + uncl


README = """bioThreatFinder - Pavian exports
=====================================

Load these in Pavian (https://github.com/fbreitwieser/pavian):

    install.packages("remotes")
    remotes::install_github("fbreitwieser/pavian")
    pavian::runApp(port = 5000)

then in "Data Selection" point the directory browser at THIS directory, or
upload the *.pavian.txt files. They are plain 6-column Kraken reports:

    percentage  cladeReads  taxonReads  taxRank  taxID  name

One file per sample per classifier, named <sample>.<classifier>.pavian.txt, so
Pavian's sample-comparison view lines the arms up side by side.

Provenance for every file, including what is measured and what is estimated,
is in <sample>.pavian_manifest.tsv next to this README.

Two things worth knowing
------------------------
1. The k-mer columns are not here. Pavian's two readers disagree about them:
   the app's read_report() accepts an 8-column Kraken2 --report-minimizer-data
   file, but read_report2() accepts only 6. These exports use 6 so that both
   work. When you want the minimizer evidence, load the untouched report from
   ../kraken2/<sample>.kraken2.report.txt directly - the app reads it fine.

2. The sylph export is SYNTHESISED, not converted. sylph reports ANI and
   containment, not read assignments, so its counts here are back-derived from
   abundance so that Pavian has something to draw. Do not quote them as read
   counts. See <sample>.pavian_manifest.tsv for the basis used.
"""


def self_test() -> int:
    """The converters, on synthetic reports. No network, no taxdump, no samtools.

    This script is presentational and was the one pipeline-invoked script with no
    coverage at all, which is the wrong combination: because it never dies (see
    the DEGRADE, DO NOT DIE note in main()), a reader that silently produces
    nonsense produces it under exit 0 and nobody looks.
    """
    import tempfile

    def w(td, name, text):
        p = Path(td) / name
        p.write_text(text)
        return p

    with tempfile.TemporaryDirectory() as td:
        # ---- 6-column Kraken-style report ----
        rows, bad = rows_from_kraken_style(w(td, "k6.txt",
            "100.00\t1000\t0\tR\t1\troot\n"
            " 50.00\t500\t500\tS\t632\t  Yersinia pestis\n"
            "  0.00\t0\t0\tU\t0\tunclassified\n"))
        assert bad == 0, bad
        assert [r[4] for r in rows] == ["1", "632", "0"], rows
        assert rows[1][5] == "  Yersinia pestis", "indentation must survive verbatim"
        assert total_reads_from(rows) == 1000, total_reads_from(rows)

        # ---- 8-column Kraken2 --report-minimizer-data: two EXTRA columns on
        # the left, which is the whole reason the fields are taken from the
        # right. The minimizer columns must not shift the taxid.
        rows8, bad8 = rows_from_kraken_style(w(td, "k8.txt",
            "100.00\t1000\t0\t900\t800\tR\t1\troot\n"
            " 50.00\t500\t500\t450\t400\tS\t632\t  Yersinia pestis\n"))
        assert bad8 == 0, bad8
        assert [r[4] for r in rows8] == ["1", "632"], rows8
        assert rows8[1][1] == 500 and rows8[1][2] == 500, rows8

        # ---- Metabuli: a '#'-prefixed header, dropped ----
        rowsm, badm = rows_from_kraken_style(w(td, "mb.txt",
            "#clade_proportion\tclade_count\tread_count\trank\ttaxID\tname\n"
            " 12.50\t125\t125\tspecies\t1392\tBacillus anthracis\n"))
        assert badm == 0 and len(rowsm) == 1, (badm, rowsm)
        assert rowsm[0][3] == "S", rowsm[0]      # rank normalised

        # ---- THE WRONG FILE. A Bracken ABUNDANCE table, which the BRACKEN
        # process writes beside the `-w` report this script wants. Before the
        # shape check it parsed with bad=0 into fabricated taxids taken from
        # new_est_reads. It must now yield NO rows.
        abund = w(td, "s.bracken.tsv",
            "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\t"
            "added_reads\tnew_est_reads\tfraction_total_reads\n"
            "Yersinia pestis\t632\tS\t400\t120\t520\t0.05200\n"
            "Bacillus anthracis\t1392\tS\t90\t15\t105\t0.01050\n")
        rows_a, bad_a = rows_from_kraken_style(abund)
        assert rows_a == [], f"an abundance table must give no rows, got {rows_a}"
        assert bad_a == 3, bad_a                 # header + both data rows
        # and the taxids it WOULD have fabricated are the read counts
        assert "520" not in {r[4] for r in rows_a}

        # ---- a report whose first field is not a percentage is refused ----
        rows_x, bad_x = rows_from_kraken_style(w(td, "x.txt",
            "sample_one\t10\t10\tS\t632\tYersinia pestis\n"))
        assert rows_x == [] and bad_x == 1, (rows_x, bad_x)

        # ---- KrakenUniq: header-aware, because its column order is not
        # Kraken's (taxID, rank, taxName at the end in a different order).
        rows_ku, bad_ku = rows_from_krakenuniq(w(td, "ku.txt",
            "%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName\n"
            "100.00\t1000\t0\t5000\t1.1\t0.01\t1\tno rank\troot\n"
            " 50.00\t500\t500\t2500\t1.2\t0.02\t632\tspecies\t  Yersinia pestis\n"))
        assert bad_ku == 0, bad_ku
        assert [r[4] for r in rows_ku] == ["1", "632"], rows_ku
        assert rows_ku[1][3] == "S", rows_ku[1]

        # ---- write_report: exactly five tabs per line, or Pavian's
        # read_report() discards the whole file.
        dest = Path(td) / "out.txt"
        write_report(dest, rows)
        lines = dest.read_text().splitlines()
        assert len(lines) == 3, lines
        assert all(ln.count("\t") == 5 for ln in lines), lines
        assert not dest.read_text().endswith("\n\n")

        # ---- normalise_rank ----
        assert normalise_rank("no rank", "1") == "R"
        assert normalise_rank("no rank", "0") == "U"
        assert normalise_rank("species", "632") == "S"

    print("to_pavian self-test OK", file=sys.stderr)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample")
    ap.add_argument("--outdir", type=Path, default=Path("."))
    ap.add_argument("--kraken2-report", type=Path)
    ap.add_argument("--bracken-report", type=Path)
    ap.add_argument("--metabuli-report", type=Path)
    ap.add_argument("--krakenuniq-report", type=Path)
    ap.add_argument("--sylph-profile", type=Path)
    ap.add_argument("--taxdump", type=Path,
                    help="taxdump directory (nodes.dmp, names.dmp); required for sylph")
    ap.add_argument("--total-reads", type=int, default=0,
                    help="read total for sylph scaling; derived from a Kraken-style "
                         "report when not given")
    ap.add_argument("--self-test", action="store_true",
                    help="exercise the converters on synthetic reports and exit")
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    if not args.sample:
        ap.error("the following arguments are required: --sample")

    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict] = []
    derived_total = args.total_reads

    passthrough = [
        ("kraken2", args.kraken2_report, rows_from_kraken_style,
         ("measured read counts; the two minimizer columns were dropped so that "
          "both of Pavian's readers accept the file (the untouched 8-column "
          "report under kraken2/ also loads in the app)")),
        ("bracken", args.bracken_report, rows_from_kraken_style,
         ("measured; Bracken-reestimated abundances, already 6-column, "
          "re-emitted unchanged for consistent naming")),
        ("metabuli", args.metabuli_report, rows_from_kraken_style,
         ("measured; its '#'-prefixed header row is dropped - Pavian would read "
          "that line as a header, fail to map the column names and reject the "
          "whole file")),
        ("krakenuniq", args.krakenuniq_report, rows_from_krakenuniq,
         ("measured; read header-aware because KrakenUniq's column order is not "
          "Kraken's; its kmers/dup/cov columns were dropped")),
    ]

    for label, path, reader, how in passthrough:
        if not usable(path):
            continue
        rows, bad = reader(path)
        if not rows:
            sys.stderr.write(f"[pavian] WARNING: {label} report {path.name} gave no "
                             f"usable rows ({bad} unparseable lines); skipped\n")
            continue
        if not derived_total:
            derived_total = total_reads_from(rows)
        dest = args.outdir / f"{args.sample}.{label}.pavian.txt"
        write_report(dest, rows)
        manifest.append({
            "sample": args.sample, "classifier": label, "file": dest.name,
            "source": str(path), "rows": len(rows), "unparseable_source_lines": bad,
            "counts": "measured", "basis": "source report", "notes": how,
        })

    if usable(args.sylph_profile):
        if not args.taxdump:
            sys.stderr.write("[pavian] WARNING: --sylph-profile given without "
                             "--taxdump; cannot build a tree, skipped\n")
        else:
            basis = derived_total or 1_000_000
            # DEGRADE, DO NOT DIE. This whole script is presentational and the
            # pipeline sets no errorStrategy on it, so an unhandled exception in
            # the sylph synthesis - a malformed profile, a surprising taxdump -
            # would terminate the run before stage 2 for the sake of a view.
            # The kraken-style exports above are already written by this point.
            try:
                rows, notes = synth_from_sylph(args.sylph_profile, args.taxdump, basis)
            except Exception as exc:                       # noqa: BLE001
                sys.stderr.write(f"[pavian] WARNING: sylph synthesis failed for "
                                 f"{args.sample}: {exc!r}; continuing without it\n")
                rows, notes = [], [f"synthesis FAILED: {exc!r}"]
            dest = args.outdir / f"{args.sample}.sylph.pavian.txt"
            # An EMPTY file is worse than no file: Pavian's readLines(n=1) returns
            # nothing, it logs "Could not read" and skips it, so all an empty
            # export does is put a dead entry in the user's file browser. sylph
            # legitimately finds nothing when no genome clears its ANI threshold.
            if rows:
                write_report(dest, rows)
            manifest.append({
                "sample": args.sample, "classifier": "sylph",
                "file": dest.name if rows else "(not written - no rows)",
                "source": str(args.sylph_profile), "rows": len(rows),
                "unparseable_source_lines": 0,
                "counts": "ESTIMATED - not read assignments",
                "basis": f"Sequence_abundance x {basis} reads"
                         + ("" if derived_total else " (NOMINAL - no read total available)"),
                "notes": "synthesised from the NCBI taxdump; sylph reports ANI and "
                         "containment, not reads. "
                         + ("; ".join(notes) if notes else "no warnings"),
            })

    # ONE MANIFEST PER SAMPLE, not one per directory. In the pipeline each
    # sample is a separate task writing into its own work directory, and
    # publishDir then copies them all into the same place - a shared filename
    # would mean every sample but the last silently disappeared.
    man_path = args.outdir / f"{args.sample}.pavian_manifest.tsv"
    cols = ["sample", "classifier", "file", "counts", "basis", "rows",
            "unparseable_source_lines", "source", "notes"]
    with man_path.open("w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in manifest:
            w.writerow(row)

    readme = args.outdir / "README.txt"
    if not readme.exists():
        readme.write_text(README)

    if not manifest:
        sys.stderr.write(f"[pavian] no classifier report was convertible for "
                         f"{args.sample}\n")
    else:
        sys.stderr.write(f"[pavian] {args.sample}: wrote "
                         f"{', '.join(m['file'] for m in manifest)}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
