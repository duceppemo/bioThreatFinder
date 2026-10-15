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

  metabuli  Metabuli emits a TRAILING TAB on unclassified rows (upstream issue
            #194, open). Unclassified is usually line 1, so the first line
            carries 6 tabs, matches neither 5 nor 7, and Pavian rejects the
            ENTIRE file. Stripping the trailing empty field fixes it.

  sylph     Not a Kraken report at all - no tree, no read counts. Synthesised
            from the NCBI taxdump; see synth_from_sylph() for what is real and
            what is estimated.
"""

import argparse
import csv
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
}

# KrakenUniq header spellings -> canonical field. Same map as parse_screen.py.
KU_COLUMNS = {
    "%": "percentage", "reads": "clade_reads", "taxreads": "taxon_reads",
    "taxid": "taxid", "rank": "rank", "taxname": "name",
}


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
def rows_from_kraken_style(path: Path) -> tuple[list[tuple], int]:
    """Kraken2 (6 or 8 col), Bracken, Metabuli -> canonical 6-field rows.

    Indexed from the RIGHT, which is what makes one reader cover all three
    layouts: they differ on the left (Kraken2's two minimizer columns) but keep
    ... rank, taxid, name stable at the end. Lifted from parse_screen.py so the
    two never disagree about what a column means.
    """
    rows, bad = [], 0
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "%")):
                continue
            f = line.rstrip("\n").split("\t")
            # Metabuli emits a trailing tab on unclassified rows (upstream #194),
            # which would shift every right-indexed field by one.
            while f and f[-1].strip() == "":
                f.pop()
            if len(f) < 6:
                bad += 1
                continue
            taxid, rank, name = f[-2].strip(), f[-3].strip(), f[-1]
            if not taxid.isdigit():
                bad += 1
                continue
            rows.append((f[0].strip(), as_int(f[1]), as_int(f[2]), rank, taxid, name))
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
                fields[idx["rank"]].strip() if "rank" in idx else "-",
                taxid,
                fields[idx["name"]] if "name" in idx else "",
            ))
    return rows, bad


# --------------------------------------------------------------------------
# Sylph synthesis
# --------------------------------------------------------------------------
def load_taxdump(taxdump: Path, wanted: set[str]) -> tuple[dict, dict, dict]:
    """parent / rank for every taxid, scientific names for the needed closure."""
    parent: dict[str, str] = {}
    rank: dict[str, str] = {}
    with (taxdump / "nodes.dmp").open() as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) >= 3:
                parent[f[0]], rank[f[0]] = f[1], f[2]

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
    return parent, rank, names


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
    for r in rows:
        gf = r.get("Genome_file") or ""
        taxid = taxid_from_genome_file(gf)
        if not taxid:
            notes.append(f"could not read a taxid from genome file name: {gf}")
            continue
        ab = as_float(r.get("Sequence_abundance"), -1.0)
        if ab < 0:
            ab = as_float(r.get("Taxonomic_abundance"), 0.0)
        leaf[taxid] = leaf.get(taxid, 0.0) + ab

    if not leaf:
        return [], notes or ["no usable sylph rows"]

    parent, rank, names = load_taxdump(taxdump, set(leaf))

    # Accumulate clade counts up every lineage.
    clade: dict[str, int] = {}
    direct: dict[str, int] = {}
    children: dict[str, set[str]] = {}
    for taxid, ab in leaf.items():
        reads = int(round(ab / 100.0 * total_reads))
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

    total_clade = clade.get("1", 0)
    out: list[tuple] = []

    unassigned = max(total_reads - total_clade, 0)
    if unassigned:
        out.append((f"{unassigned / total_reads * 100:.2f}", unassigned, unassigned,
                    "U", "0", "unclassified"))

    def walk(taxid: str, depth: int) -> None:
        pct = clade[taxid] / total_reads * 100 if total_reads else 0.0
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
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
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict] = []
    derived_total = args.total_reads

    passthrough = [
        ("kraken2", args.kraken2_report, rows_from_kraken_style,
         "measured read counts; the two minimizer columns were dropped so that both "
         "of Pavian's readers accept the file (the untouched 8-column report under "
         "kraken2/ also loads in the app)"),
        ("bracken", args.bracken_report, rows_from_kraken_style,
         "measured; Bracken-reestimated abundances, already 6-column, re-emitted "
         "unchanged for consistent naming"),
        ("metabuli", args.metabuli_report, rows_from_kraken_style,
         "measured; trailing empty field stripped from unclassified rows (upstream "
         "issue #194) - without this Pavian rejects the entire file"),
        ("krakenuniq", args.krakenuniq_report, rows_from_krakenuniq,
         "measured; read header-aware because KrakenUniq's column order is not "
         "Kraken's; its kmers/dup/cov columns were dropped"),
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
            rows, notes = synth_from_sylph(args.sylph_profile, args.taxdump, basis)
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
