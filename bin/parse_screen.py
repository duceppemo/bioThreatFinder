#!/usr/bin/env python3
"""
Intersect stage-1 screen output with the watchlist to produce stage-2 candidates.

MULTI-CLASSIFIER BY DESIGN

Stage 1 runs up to four independent classifiers, because they fail in different
directions and no single one is good enough for a "is organism X present" question:

  krakenuniq  PRIMARY nucleotide screen. Unlike Kraken2 it reports the number of
              DISTINCT k-mers hit per taxon (HyperLogLog cardinality) plus a
              duplication ratio and an estimated genome coverage. This is the
              specificity fix: Kraken2 counts READS, so 500 reads piled on one
              conserved k-mer are indistinguishable from 500 reads spread across a
              genome. KrakenUniq tells those apart at stage 1, before any
              alignment is spent on them.
  metabuli    Joint DNA + amino-acid k-mer screen. Sees divergent and novel agents
              that nucleotide k-mers miss entirely. Until now only the marker tier
              could see through divergence, and only for curated markers.
  kraken2     Retained for abundance estimation via Bracken and for continuity
              with existing runs. No longer the primary evidence.
  sylph       ANI/containment based. Very low false-positive rate for bacteria and
              fungi, and it disagrees with k-mer methods in exactly the cases that
              matter - a taxon assembled from scattered conserved k-mers gets no
              sylph support.

A candidate is raised if ANY classifier supports it. This is deliberate: stage 1
stays permissive and stage 2 does the discriminating. Requiring consensus here
would discard the sensitivity that running four classifiers is meant to buy. What
the consensus DOES do is get recorded - `n_classifiers` and `classifiers` travel
into scoring and into the report, because a taxon seen by krakenuniq, metabuli and
sylph is far stronger evidence than one seen by kraken2 alone.

HYPERLOGLOG CARDINALITIES ARE NOT ADDITIVE

A subtle but load-bearing detail. KrakenUniq's `kmers` column is a cardinality
ESTIMATE for the clade rooted at that node. Summing it across descendant nodes
would double-count shared k-mers and inflate the number, sometimes badly. So for
distinct-k-mer evidence this script uses the watchlist node's OWN clade row, and
falls back to the MAXIMUM over descendants rather than the sum. Read counts, being
plain integers over disjoint assignments, are summed as normal.

INTERNAL-NODE ASSIGNMENT

Classifiers assign many reads to internal nodes - a read consistent with several
species in a genus is placed at the genus. Matching only the species taxid misses
most of the signal, which is the difference between a screen that works and one
that quietly reports nothing. Every watchlist taxon is therefore expanded to all
its descendants (precomputed by scripts/resolve_watchlist.py), and ancestor-level
assignment is recorded separately as `ancestor_only` rather than dropped: a
genus-level pile-up with no species-level reads is exactly what a divergent or
under-represented strain produces.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

# KrakenUniq report columns, by header name. Parsed by name when a header is
# present so a column-order change cannot silently shift the k-mer counts into
# the wrong field.
# VERIFY: header spellings against the installed krakenuniq version on first run.
KU_COLUMNS = {
    "%": "pct", "reads": "clade_reads", "taxreads": "direct_reads",
    "kmers": "kmers", "dup": "dup", "cov": "cov",
    "taxid": "taxid", "rank": "rank", "taxname": "name",
}


def read_tsv(path: Path | None) -> list[dict]:
    if not path or not path.exists() or path.name == "NO_FILE":
        return []
    with path.open() as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t")) if lines else []


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


def unclassified_fraction(table: dict[str, dict]) -> tuple[int, int]:
    """(unclassified_reads, total_reads) from a kraken-style report.

    WHY THIS IS A FIRST-CLASS OUTPUT, NOT A DIAGNOSTIC

    Perseus (Nguyen & Schatz, Bioinformatics 2026, doi:10.1093/bioinformatics/btag687)
    measured what Kraken2 does with sequence from a family absent from its database:
    the false-assignment rate approaches 50% on long reads. Divergent reads are not
    left unclassified - they are CONFIDENTLY ASSIGNED TO THE WRONG KNOWN TAXON. One
    family-exclusion scenario produced misclassifications while 89% of k-mers went
    unclassified, and species assignments were frequently supported by a median of
    only five k-mers.

    For a surveillance tool the consequence is specific and severe: a novel agent is
    not missed, it is mislabelled as something familiar and then dismissed as
    background. The unclassified fraction is the only signal that survives that
    failure mode, so it is reported on every sample rather than thrown away.
    """
    unc = table.get("0")
    unclassified = unc["clade_reads"] if unc else 0
    root = table.get("1")
    classified = root["clade_reads"] if root else sum(
        v["direct_reads"] for v in table.values())
    return unclassified, unclassified + classified


def parse_kraken_style(path: Path | None, label: str) -> dict[str, dict]:
    """Kraken2 / Bracken / Metabuli report: pct, clade, direct, rank, taxid, name.

    Index from the right-hand end: a 6-column report, an 8-column
    --report-minimizer-data report and Metabuli's variant all keep the trailing
    layout stable while differing on the left.
    """
    out: dict[str, dict] = {}
    if not path or not path.exists() or path.name == "NO_FILE":
        return out
    bad = 0
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("%"):
                continue
            # Metabuli emits a TRAILING TAB on unclassified rows (upstream issue
            # #194, open), which shifts every right-indexed field by one and would
            # silently mis-read the taxid column. Strip trailing empties before
            # indexing from the right.
            f = line.rstrip("\n").split("\t")
            while f and f[-1].strip() == "":
                f.pop()
            if len(f) < 6:
                bad += 1
                continue
            taxid, rank, name = f[-2].strip(), f[-3].strip(), f[-1].strip()
            if not taxid.isdigit():
                bad += 1
                continue
            out[taxid] = {"clade_reads": as_int(f[1]), "direct_reads": as_int(f[2]),
                          "rank": rank, "name": name}
    if bad and not out:
        sys.stderr.write(
            f"[screen] WARNING: {label} report {path.name} yielded no usable rows "
            f"({bad} unparseable lines). Is it the right format?\n")
    return out


def parse_krakenuniq(path: Path | None) -> dict[str, dict]:
    """KrakenUniq report, header-aware. Carries distinct-k-mer evidence."""
    out: dict[str, dict] = {}
    if not path or not path.exists() or path.name == "NO_FILE":
        return out

    header: list[str] | None = None
    with path.open() as fh:
        for line in fh:
            raw = line.rstrip("\n")
            if not raw.strip():
                continue
            fields = raw.split("\t")
            low = [c.strip().lower().lstrip("#").strip() for c in fields]

            if header is None and any(c in ("kmers", "dup", "cov") for c in low):
                header = low
                continue
            if raw.startswith("#"):
                continue

            if header:
                idx = {}
                for i, col in enumerate(header):
                    key = KU_COLUMNS.get(col)
                    if key:
                        idx[key] = i
                if "taxid" not in idx:
                    continue
                if len(fields) <= max(idx.values()):
                    continue
                taxid = fields[idx["taxid"]].strip()
                if not taxid.isdigit():
                    continue
                out[taxid] = {
                    "clade_reads": as_int(fields[idx.get("clade_reads", -1)]
                                         if "clade_reads" in idx else 0),
                    "direct_reads": as_int(fields[idx.get("direct_reads", -1)]
                                          if "direct_reads" in idx else 0),
                    "kmers": as_int(fields[idx["kmers"]]) if "kmers" in idx else 0,
                    "dup": as_float(fields[idx["dup"]]) if "dup" in idx else 0.0,
                    "cov": as_float(fields[idx["cov"]]) if "cov" in idx else 0.0,
                    "rank": fields[idx["rank"]].strip() if "rank" in idx else "",
                    "name": fields[idx["name"]].strip() if "name" in idx else "",
                }
            else:
                # No recognisable header. Fall back to the documented positional
                # layout, but say so - silently guessing which column holds the
                # k-mer count is exactly the kind of quiet mis-parse that must not
                # happen in this file.
                if len(fields) < 9 or not fields[6].strip().isdigit():
                    continue
                out[fields[6].strip()] = {
                    "clade_reads": as_int(fields[1]), "direct_reads": as_int(fields[2]),
                    "kmers": as_int(fields[3]), "dup": as_float(fields[4]),
                    "cov": as_float(fields[5]), "rank": fields[7].strip(),
                    "name": fields[8].strip(),
                }

    if out and header is None:
        sys.stderr.write(
            f"[screen] WARNING: no header found in KrakenUniq report {path.name}; "
            "fell back to positional parsing (%, reads, taxReads, kmers, dup, cov, "
            "taxID, rank, taxName). Verify the column order.\n")
    return out


def parse_sylph(path: Path | None) -> dict[str, float]:
    """Sylph reports genome file paths, not taxids, so match any taxid-looking
    token in the filename. Loose by necessity; treated as supporting evidence
    only, never as the sole basis for a candidate's thresholds."""
    out: dict[str, float] = {}
    for row in read_tsv(path):
        ab = row.get("Taxonomic_abundance") or row.get("Sequence_abundance") or "0"
        gf = row.get("Genome_file") or row.get("Contig_name") or ""
        val = as_float(ab, -1.0)
        if val < 0:
            continue
        for token in gf.replace("/", "_").replace(".", "_").split("_"):
            if token.isdigit() and len(token) >= 2:
                out[token] = max(out.get(token, 0.0), val)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--krakenuniq-report", type=Path,
                    help="primary nucleotide screen; carries distinct-k-mer evidence")
    ap.add_argument("--kraken2-report", type=Path)
    ap.add_argument("--bracken-report", type=Path)
    ap.add_argument("--metabuli-report", type=Path,
                    help="joint DNA+AA screen; catches divergent agents")
    ap.add_argument("--sylph-profile", type=Path)
    ap.add_argument("--watchlist", required=True, type=Path)
    ap.add_argument("--descendants", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--role", default="sample")
    ap.add_argument("--platform", default="illumina")
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--min-distinct-kmers", type=int, default=150,
                    help="below this, distinct-k-mer evidence is flagged as weak "
                         "(FLAG ONLY - candidates are not dropped here) [150]")
    ap.add_argument("--max-kmer-dup", type=float, default=25.0,
                    help="reads-per-distinct-k-mer above this is flagged as piled "
                         "up on few k-mers (FLAG ONLY) [25.0]")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--out-summary", type=Path,
                    help="per-sample screen summary incl. the unclassified fraction")
    ap.add_argument("--high-unclassified", type=float, default=0.40,
                    help="warn above this unclassified fraction - in a surveillance "
                         "context it is the signature of sequence the databases cannot "
                         "name [0.40]")
    args = ap.parse_args()

    watch = read_tsv(args.watchlist)
    if not watch:
        sys.exit(f"ERROR: empty watchlist {args.watchlist}")

    desc_of: dict[str, set[str]] = defaultdict(set)
    with args.descendants.open() as fh:
        first = next(fh, "")
        if "descendant" not in first:
            fh.seek(0)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                desc_of[parts[0]].add(parts[2])

    ku = parse_krakenuniq(args.krakenuniq_report)
    k2 = parse_kraken_style(args.kraken2_report, "kraken2")
    br = parse_kraken_style(args.bracken_report, "bracken")
    mb = parse_kraken_style(args.metabuli_report, "metabuli")
    sy = parse_sylph(args.sylph_profile)

    supplied = [n for n, d in (("krakenuniq", ku), ("kraken2", k2),
                               ("metabuli", mb), ("sylph", sy)) if d]
    if not supplied:
        sys.stderr.write(
            "[screen] WARNING: no classifier produced usable output for "
            f"{args.sample}. No candidates will be raised - this is a pipeline "
            "failure, not a clean negative.\n")
    else:
        sys.stderr.write(f"[screen] {args.sample}: classifiers with output: "
                         f"{', '.join(supplied)}\n")

    rows = []
    for entry in watch:
        taxid = entry["taxid"]
        family = desc_of.get(taxid) or {taxid}

        def clade_or_sum(table: dict[str, dict],
                         taxid: str = taxid, family: set = family) -> tuple[int, int]:
            """(reads, direct_reads) for this watchlist clade.

            taxid/family are bound as defaults rather than captured: the closure is
            called inside this iteration so late binding is harmless today, but the
            pattern fails silently if it ever outlives the loop.
            """
            self_hit = table.get(taxid)
            direct = sum(table[d]["direct_reads"] for d in family
                         if d in table and table[d]["direct_reads"] > 0)
            clade = self_hit["clade_reads"] if self_hit else 0
            return max(direct, clade), direct

        ku_reads, ku_direct = clade_or_sum(ku)
        k2_reads, k2_direct = clade_or_sum(k2)
        mb_reads, _ = clade_or_sum(mb)
        br_reads = sum(br[d]["clade_reads"] for d in family if d in br)

        # HyperLogLog cardinalities are NOT additive - see module docstring.
        # Prefer the node's own clade estimate; fall back to the max over
        # descendants, never the sum.
        ku_self = ku.get(taxid)
        if ku_self:
            ku_kmers = ku_self["kmers"]
            ku_dup = ku_self["dup"]
            ku_cov = ku_self["cov"]
            kmer_basis = "clade"
        else:
            cands = [ku[d] for d in family if d in ku]
            if cands:
                best = max(cands, key=lambda r: r["kmers"])
                ku_kmers, ku_dup, ku_cov = best["kmers"], best["dup"], best["cov"]
                kmer_basis = "max_descendant"
            else:
                ku_kmers, ku_dup, ku_cov = 0, 0.0, 0.0
                kmer_basis = "none"

        sylph_ab = sy.get(taxid, 0.0)

        # Which classifiers support this taxon, and how many.
        support = []
        if ku_reads >= args.min_reads:
            support.append("krakenuniq")
        if k2_reads >= args.min_reads:
            support.append("kraken2")
        if mb_reads >= args.min_reads:
            support.append("metabuli")
        if sylph_ab > 0:
            support.append("sylph")

        # Permissive OR: any classifier raises the candidate. Stage 2 decides.
        if not support:
            continue

        # Distinct-k-mer quality. FLAGGED, NOT ENFORCED - dropping here would
        # undo the architecture's "stage 1 permissive" rule, and a flag that
        # reaches scoring is more useful than a candidate that silently vanished.
        flags = []
        if kmer_basis == "none":
            flags.append("no_kmer_evidence")
        else:
            if ku_kmers < args.min_distinct_kmers:
                flags.append(f"low_distinct_kmers({ku_kmers}<{args.min_distinct_kmers})")
            implied_dup = (ku_reads / ku_kmers) if ku_kmers else 0.0
            if ku_dup > args.max_kmer_dup or implied_dup > args.max_kmer_dup:
                flags.append(f"high_kmer_duplication(dup={ku_dup:.1f})")
        if not args.krakenuniq_report:
            flags = ["krakenuniq_not_run"]

        rows.append({
            "sample": args.sample,
            "role": args.role,
            "platform": args.platform,
            "taxid": taxid,
            "name": entry["name"],
            "group": entry.get("group", ""),
            "category": entry.get("category", ""),
            # `assigned_reads` stays the headline number downstream consumers use;
            # it is now the best evidence across classifiers rather than kraken2's.
            "assigned_reads": max(ku_reads, k2_reads, mb_reads),
            "krakenuniq_reads": ku_reads,
            "krakenuniq_distinct_kmers": ku_kmers,
            "krakenuniq_dup": f"{ku_dup:.3f}" if kmer_basis != "none" else "",
            "krakenuniq_cov": f"{ku_cov:.6f}" if kmer_basis != "none" else "",
            "kmer_basis": kmer_basis,
            "kraken2_reads": k2_reads,
            "bracken_reads": br_reads,
            "metabuli_reads": mb_reads,
            "sylph_abundance": f"{sylph_ab:.6f}" if sylph_ab else "",
            "n_classifiers": len(support),
            "classifiers": ";".join(support),
            "screen_flags": ";".join(flags),
            "direct_reads": max(ku_direct, k2_direct),
            "ancestor_only": "yes" if (max(ku_reads, k2_reads) > 0
                                       and max(ku_direct, k2_direct) == 0) else "no",
            "n_descendants_searched": len(family),
            "near_neighbour_taxids": entry.get("near_neighbour_taxids", ""),
            "min_unique_reads": entry.get("min_unique_reads", ""),
            "min_breadth": entry.get("min_breadth", ""),
            "notes": entry.get("notes", ""),
        })

    # Data-side detection of the same taxonomy-namespace mismatch the pipeline
    # refuses up front. If Metabuli produced plenty of rows but matched no watchlist
    # taxon at all while another arm did, the most likely cause by far is that its
    # index is keyed in a different namespace (GTDB taxids are not NCBI taxids). Say
    # so explicitly - the alternative is an arm that looks like it ran and found
    # nothing, which is indistinguishable from a clean negative.
    if mb and not any(r["metabuli_reads"] > 0 for r in rows):
        other = any(r["krakenuniq_reads"] > 0 or r["kraken2_reads"] > 0 for r in rows)
        if other:
            sys.stderr.write(
                f"[screen] WARNING: the metabuli report held {len(mb)} taxon rows but "
                "matched NO watchlist taxid, while another arm did match. The usual "
                "cause is a taxonomy-namespace mismatch - a GTDB-based Metabuli index "
                "uses GTDB taxids, which do not join an NCBI-keyed watchlist. This arm "
                "is contributing nothing. Verify with "
                "scripts/build_metabuli_db.sh --mode refseq, or pass --skip_metabuli.\n")

    # ---- unclassified fraction, per classifier ----
    summary = []
    for label, table in (("krakenuniq", ku), ("kraken2", k2), ("metabuli", mb)):
        if not table:
            continue
        unc, tot = unclassified_fraction(table)
        frac = (unc / tot) if tot else 0.0
        summary.append({
            "sample": args.sample, "role": args.role, "platform": args.platform,
            "classifier": label, "total_reads": tot, "unclassified_reads": unc,
            "unclassified_fraction": f"{frac:.6f}",
            "flag": ("HIGH_UNCLASSIFIED" if frac >= args.high_unclassified else ""),
        })
        if frac >= args.high_unclassified:
            sys.stderr.write(
                f"[screen] {args.sample}: {label} left {frac:.1%} of reads UNCLASSIFIED "
                f"({unc}/{tot}). In a surveillance context a high unclassified fraction "
                "is the signature of sequence the databases cannot name. It is also the "
                "only signal that survives the failure mode where a divergent agent is "
                "confidently mislabelled as a known taxon rather than missed - so "
                "triage this fraction, do not discard it.\n")

    if args.out_summary and summary:
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        with args.out_summary.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary[0]), delimiter="\t")
            w.writeheader()
            for r in summary:
                w.writerow(r)

    rows.sort(key=lambda r: (-r["n_classifiers"], -r["assigned_reads"]))

    cols = list(rows[0]) if rows else [
        "sample", "role", "platform", "taxid", "name", "group", "category",
        "assigned_reads", "krakenuniq_reads", "krakenuniq_distinct_kmers",
        "krakenuniq_dup", "krakenuniq_cov", "kmer_basis", "kraken2_reads",
        "bracken_reads", "metabuli_reads", "sylph_abundance", "n_classifiers",
        "classifiers", "screen_flags", "direct_reads", "ancestor_only",
        "n_descendants_searched", "near_neighbour_taxids", "min_unique_reads",
        "min_breadth", "notes",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    sys.stderr.write(
        f"[screen] {args.sample}: {len(rows)} of {len(watch)} watchlist taxa raised "
        f"as stage-2 candidates (min_reads={args.min_reads}, any-classifier OR)\n")
    for r in rows[:15]:
        extra = []
        if r["ancestor_only"] == "yes":
            extra.append("ancestor-level only")
        if r["screen_flags"]:
            extra.append(r["screen_flags"])
        tail = f" [{'; '.join(extra)}]" if extra else ""
        sys.stderr.write(
            f"[screen]   {r['name']} (taxid {r['taxid']}): {r['assigned_reads']} reads, "
            f"{r['krakenuniq_distinct_kmers']} distinct k-mers, "
            f"{r['n_classifiers']} classifier(s) [{r['classifiers']}]{tail}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
