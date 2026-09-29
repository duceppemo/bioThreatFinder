#!/usr/bin/env python3
"""
Turn marker search output into three-valued marker states.

Feeds bin/marker_rules.py. The contract with that engine is that PRESENT,
ABSENT and UNKNOWN are three distinct facts:

  PRESENT  the marker was searched for and found above its thresholds
  ABSENT   the marker was searched for and confidently not found
  UNKNOWN  the marker was never actually searched - it is missing from the
           marker database, retrieval failed at build time, or no input for its
           evidence tier was supplied to this run

The manifest written by scripts/build_marker_db.py is AUTHORITATIVE for that
distinction. A marker whose manifest status is not `built` or `partial` can only
ever be UNKNOWN here, no matter what the alignment files contain. This is the
whole point: a gap in the marker database must not read as a clean negative, or
the pipeline will cheerfully report "RD1 absent, therefore BCG" for a sample it
never tested, and "no pXO2, therefore not anthrax" for a marker it never had.

THREE EVIDENCE TIERS, most to least specific

  dna      read alignment to the marker nucleotide sequence. Precise, but blind
           to a divergent allele.
  protein  translated search (MMseqs2/DIAMOND blastx). Sees through synonymous
           divergence, so it catches an allele that nucleotide alignment misses.
  hmm      profile search. The most remote-homology-sensitive tier. This is what
           catches a novel BoNT or Stx subtype, or a divergent filovirus gene,
           which is the difference between a tool that detects a known strain and
           one that detects a new agent.

A marker is PRESENT if ANY applicable tier passes. Tiers are complementary, not
a consensus vote - requiring agreement would throw away exactly the sensitivity
the protein and HMM tiers exist to provide.

Needs only `samtools` on PATH. No third-party Python.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_OPS = {"M", "D", "N", "=", "X"}
ALN_OPS = {"M", "I", "=", "X"}

PRESENT, ABSENT, UNKNOWN = "PRESENT", "ABSENT", "UNKNOWN"
BUILT_STATUSES = {"built", "partial"}


def stream(cmd: list[str]):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    for line in proc.stdout:
        yield line
    proc.stdout.close()
    rc = proc.wait()
    if rc != 0:
        err = proc.stderr.read() if proc.stderr else ""
        sys.exit(f"ERROR: command failed ({rc}): {' '.join(cmd)}\n{err}")


def read_tsv(path: Path | None) -> list[dict]:
    if not path or not path.exists() or path.name == "NO_FILE":
        return []
    with path.open() as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t")) if lines else []


def marker_of(header: str) -> str | None:
    """`marker|<marker_id>|<accession>|<desc>` -> marker_id.

    Split with maxsplit so a description containing '|' cannot shift the field
    we care about.
    """
    if not header.startswith("marker|"):
        return None
    parts = header.split("|", 3)
    return parts[1] if len(parts) >= 2 and parts[1] else None


def union_length(intervals: list[tuple[int, int]]) -> int:
    """Total length covered by a set of possibly overlapping 1-based closed
    intervals. Summing raw lengths would double-count stacked hits and inflate
    coverage past 1.0."""
    if not intervals:
        return 0
    merged = []
    for start, end in sorted(intervals):
        if start > end:
            start, end = end, start
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return sum(e - s + 1 for s, e in merged)


def fasta_lengths(path: Path | None) -> dict[str, int]:
    """marker_id -> summed sequence length, from a combined marker FASTA."""
    out: dict[str, int] = defaultdict(int)
    if not path or not path.exists():
        return {}
    current = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                current = marker_of(line[1:].strip())
            elif current:
                out[current] += len(line.strip())
    return dict(out)


# ---------------------------------------------------------------------------
# Tier: nucleotide read alignment
# ---------------------------------------------------------------------------

def nt_evidence(bam: Path, mapq_min: int) -> dict[str, dict]:
    lengths: dict[str, int] = defaultdict(int)
    ref_to_marker: dict[str, str] = {}
    for line in stream(["samtools", "view", "-H", str(bam)]):
        if not line.startswith("@SQ"):
            continue
        name = ln = None
        for field in line.rstrip("\n").split("\t")[1:]:
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                ln = int(field[3:])
        if name and ln is not None:
            mid = marker_of(name)
            if mid:
                ref_to_marker[name] = mid
                lengths[mid] += ln

    if not ref_to_marker:
        sys.stderr.write(
            f"[marker] WARNING: no reference in {bam.name} uses the "
            "'marker|<id>|...' convention; nucleotide tier yields nothing\n"
        )
        return {}

    agg = defaultdict(lambda: {"reads": 0, "nm": 0, "aln": 0})
    for line in stream(["samtools", "view", "-F", "0x904", "-q", str(mapq_min), str(bam)]):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        mid = ref_to_marker.get(f[2])
        if not mid:
            continue
        aln = 0
        for num, op in CIGAR_RE.findall(f[5]):
            if op in ALN_OPS:
                aln += int(num)
        nm = 0
        for tag in f[11:]:
            if tag.startswith("NM:i:"):
                try:
                    nm = int(tag[5:])
                except ValueError:
                    nm = 0
                break
        rec = agg[mid]
        rec["reads"] += 1
        rec["nm"] += nm
        rec["aln"] += aln

    covered: dict[str, int] = defaultdict(int)
    for line in stream(["samtools", "depth", "-a", "-Q", str(mapq_min), str(bam)]):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        mid = ref_to_marker.get(f[0])
        if mid and int(f[2]) > 0:
            covered[mid] += 1

    out = {}
    for mid, length in lengths.items():
        rec = agg.get(mid, {"reads": 0, "nm": 0, "aln": 0})
        out[mid] = {
            "nt_reads": rec["reads"],
            "nt_breadth": (covered.get(mid, 0) / length) if length else 0.0,
            "nt_identity": (1.0 - rec["nm"] / rec["aln"]) if rec["aln"] else 0.0,
            "nt_length": length,
        }
    return out


# ---------------------------------------------------------------------------
# Tier: translated search (MMseqs2 / DIAMOND tabular m8)
# ---------------------------------------------------------------------------

def aa_evidence(hits: Path, aa_lengths: dict[str, int]) -> dict[str, dict]:
    """m8: query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bits

    The marker protein is the TARGET; reads are queries. Marker coverage is the
    union of target intervals, so many reads stacked on one domain do not read as
    full-length coverage.
    """
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    best: dict[str, dict] = {}
    with hits.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            mid = marker_of(f[1])
            if not mid:
                continue
            try:
                fident = float(f[2])
                tstart, tend = int(f[8]), int(f[9])
                evalue, bits = float(f[10]), float(f[11])
            except ValueError:
                continue
            # MMseqs2 reports fident as a fraction, DIAMOND as a percentage.
            if fident > 1.0:
                fident /= 100.0
            intervals[mid].append((tstart, tend))
            cur = best.get(mid)
            if cur is None or bits > cur["aa_bits"]:
                best[mid] = {"aa_identity": fident, "aa_bits": bits,
                             "aa_evalue": evalue, "aa_hits": 0}
    out = {}
    for mid, ivs in intervals.items():
        length = aa_lengths.get(mid, 0)
        rec = dict(best.get(mid, {"aa_identity": 0.0, "aa_bits": 0.0, "aa_evalue": 1.0}))
        rec["aa_hits"] = len(ivs)
        rec["aa_coverage"] = (union_length(ivs) / length) if length else 0.0
        rec["aa_length"] = length
        out[mid] = rec
    return out


# ---------------------------------------------------------------------------
# Tier: profile HMM (hmmsearch --domtblout)
# ---------------------------------------------------------------------------

def hmm_evidence(domtbl: Path) -> dict[str, dict]:
    """Profiles are the QUERY, translated reads are targets.

    Columns (0-based): 3 query name (= marker_id), 5 qlen, 13 domain score,
    15 hmm-from, 16 hmm-to.
    """
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    qlen: dict[str, int] = {}
    best: dict[str, float] = defaultdict(float)
    with domtbl.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 17:
                continue
            mid = marker_of(f[3]) or f[3]
            try:
                qlen[mid] = int(f[5])
                score = float(f[13])
                hfrom, hto = int(f[15]), int(f[16])
            except ValueError:
                continue
            intervals[mid].append((hfrom, hto))
            best[mid] = max(best[mid], score)
    out = {}
    for mid, ivs in intervals.items():
        length = qlen.get(mid, 0)
        out[mid] = {
            "hmm_bits": best[mid],
            "hmm_coverage": (union_length(ivs) / length) if length else 0.0,
            "hmm_domains": len(ivs),
        }
    return out


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--markers", required=True, type=Path, help="markers.tsv catalogue")
    ap.add_argument("--manifest", type=Path,
                    help="markers.manifest.tsv - AUTHORITATIVE for tested vs untested")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--nt-bam", type=Path, help="reads aligned to markers.nt.fna")
    ap.add_argument("--aa-hits", type=Path, help="MMseqs2/DIAMOND m8 output")
    ap.add_argument("--aa-fasta", type=Path,
                    help="markers.aa.faa, for marker protein lengths (needed with --aa-hits)")
    ap.add_argument("--hmm-hits", type=Path, help="hmmsearch --domtblout output")
    ap.add_argument("--mapq-min", type=int, default=10)
    ap.add_argument("--hmm-min-bits", type=float, default=25.0,
                    help="minimum domain bitscore for the HMM tier [25]")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    catalogue = read_tsv(args.markers)
    if not catalogue:
        sys.exit(f"ERROR: no markers read from {args.markers}")

    # Manifest status decides what CAN be evaluated at all.
    status: dict[str, str] = {}
    manifest_note: dict[str, str] = {}
    if args.manifest:
        rows = read_tsv(args.manifest)
        if not rows:
            sys.stderr.write(
                f"[marker] WARNING: manifest {args.manifest} is empty or unreadable; "
                "every marker will be reported UNKNOWN\n"
            )
        for row in rows:
            mid = row.get("marker_id")
            if mid:
                status[mid] = (row.get("status") or "").strip().lower()
                manifest_note[mid] = row.get("note", "")
    else:
        sys.stderr.write(
            "[marker] WARNING: no --manifest supplied. Cannot distinguish 'not in "
            "database' from 'absent from sample', so markers with no hits will be "
            "reported UNKNOWN rather than ABSENT. This is the safe direction, but "
            "supply the manifest for correct ABSENT calls.\n"
        )

    nt = nt_evidence(args.nt_bam, args.mapq_min) if args.nt_bam else {}
    aa = {}
    if args.aa_hits:
        if not args.aa_fasta:
            sys.stderr.write(
                "[marker] WARNING: --aa-hits without --aa-fasta; marker protein "
                "lengths are unknown so aa_coverage will be 0 and the protein tier "
                "will be judged on identity and bitscore alone\n"
            )
        aa = aa_evidence(args.aa_hits, fasta_lengths(args.aa_fasta))
    hmm = hmm_evidence(args.hmm_hits) if args.hmm_hits else {}

    tiers_supplied = {
        "dna": args.nt_bam is not None,
        "protein": args.aa_hits is not None,
        "hmm": args.hmm_hits is not None,
    }

    results = []
    counts = defaultdict(int)

    for row in catalogue:
        mid = row.get("marker_id")
        if not mid:
            continue
        requested = {t.strip() for t in (row.get("evidence_type") or "").split("+") if t.strip()}
        try:
            min_cov = float(row.get("min_coverage") or 0.0)
        except ValueError:
            min_cov = 0.0
        try:
            min_id = float(row.get("min_identity") or 0.0)
        except ValueError:
            min_id = 0.0

        rec = {
            "sample": args.sample, "marker_id": mid,
            "taxon": row.get("taxon", ""), "gene": row.get("gene", ""),
            "marker_class": row.get("marker_class", ""),
            "copy_expectation": row.get("copy_expectation", ""),
            "state": UNKNOWN, "tier_used": "", "reason": "",
            "nt_reads": "", "nt_breadth": "", "nt_identity": "",
            "aa_hits": "", "aa_coverage": "", "aa_identity": "", "aa_bits": "",
            "hmm_domains": "", "hmm_coverage": "", "hmm_bits": "",
        }

        # Gate 1: was this marker ever built?
        mstat = status.get(mid)
        if args.manifest and mstat not in BUILT_STATUSES:
            rec["reason"] = (
                f"not in marker database (manifest status "
                f"{mstat or 'missing'}{': ' + manifest_note[mid] if manifest_note.get(mid) else ''})"
            )
            results.append(rec)
            counts[UNKNOWN] += 1
            continue

        # Gate 2: were inputs for any requested tier supplied to this run?
        evaluable = {t for t in requested if tiers_supplied.get(t)}
        if not evaluable:
            rec["reason"] = (
                f"no input supplied for requested tier(s) {'+'.join(sorted(requested)) or 'none'}"
            )
            results.append(rec)
            counts[UNKNOWN] += 1
            continue

        passed, tried = [], []

        if "dna" in evaluable:
            e = nt.get(mid)
            if e:
                rec.update({k: (f"{v:.4f}" if isinstance(v, float) else v)
                            for k, v in e.items() if k != "nt_length"})
                tried.append("dna")
                # AT LEAST ONE READ, WHATEVER THE THRESHOLDS SAY.
                #
                # nt_evidence() emits a record for every marker in the BAM header,
                # including those with no alignments at all (nt_breadth 0.0,
                # nt_identity 0.0). With min_coverage and min_identity both 0.00 -
                # which bmallei_deletions carries - the gate passed on zero evidence
                # and the marker was reported PRESENT. Any marker with a 0 floor
                # became a permanent free PRESENT the moment a rule referenced it.
                if (e["nt_reads"] >= 1
                        and e["nt_breadth"] >= min_cov
                        and e["nt_identity"] >= min_id):
                    passed.append(
                        f"dna(breadth={e['nt_breadth']:.3f}>={min_cov:g}, "
                        f"id={e['nt_identity']:.3f}>={min_id:g}, reads={e['nt_reads']})")
            else:
                tried.append("dna")

        if "protein" in evaluable:
            e = aa.get(mid)
            if e:
                rec.update({
                    "aa_hits": e["aa_hits"],
                    "aa_coverage": f"{e['aa_coverage']:.4f}",
                    "aa_identity": f"{e['aa_identity']:.4f}",
                    "aa_bits": f"{e['aa_bits']:.1f}",
                })
                tried.append("protein")
                # Protein identity floors are lower than nucleotide by nature;
                # scale rather than reusing min_identity directly, since a
                # synonymous-divergent allele is the case this tier exists for.
                aa_min_id = max(0.50, min_id - 0.15)
                if e["aa_identity"] >= aa_min_id and (
                        e["aa_coverage"] >= min_cov or e["aa_length"] == 0):
                    passed.append(
                        f"protein(id={e['aa_identity']:.3f}>={aa_min_id:.2f}, "
                        f"cov={e['aa_coverage']:.3f}, bits={e['aa_bits']:.0f})")
            else:
                tried.append("protein")

        if "hmm" in evaluable:
            e = hmm.get(mid)
            if e:
                rec.update({
                    "hmm_domains": e["hmm_domains"],
                    "hmm_coverage": f"{e['hmm_coverage']:.4f}",
                    "hmm_bits": f"{e['hmm_bits']:.1f}",
                })
                tried.append("hmm")
                if e["hmm_bits"] >= args.hmm_min_bits:
                    passed.append(
                        f"hmm(bits={e['hmm_bits']:.0f}>={args.hmm_min_bits:g}, "
                        f"cov={e['hmm_coverage']:.3f})")
            else:
                tried.append("hmm")

        if passed:
            rec["state"] = PRESENT
            rec["tier_used"] = ",".join(p.split("(")[0] for p in passed)
            rec["reason"] = "; ".join(passed)
        elif not args.manifest:
            # THE WARNING ABOVE PROMISED THIS AND THE CODE DID THE OPPOSITE.
            #
            # Gate 1 is `if args.manifest and mstat not in BUILT_STATUSES`, so with no
            # manifest it was skipped entirely and a marker that was never built fell
            # through to ABSENT - while the warning printed at startup asserted such
            # markers "will be reported UNKNOWN rather than ABSENT. This is the safe
            # direction." The unsafe direction has teeth: genuine M. tuberculosis
            # reads make mtb_rpoB PRESENT, mtb_RD1 was never retrieved into the
            # database, so `vaccine_if: absent(mtb_RD1) and mtb_rpoB` evaluates TRUE
            # and a real TB detection is reported as a BCG stand-down - verbatim the
            # failure this module exists to prevent.
            rec["state"] = UNKNOWN
            rec["tier_used"] = ",".join(sorted(set(tried)))
            rec["reason"] = (
                f"searched via {rec['tier_used'] or 'none'} and not found, but no "
                f"--manifest was supplied, so 'never built' cannot be distinguished "
                f"from 'absent from sample'. Reported UNKNOWN, which cannot promote "
                f"or veto a finding."
            )
        else:
            rec["state"] = ABSENT
            rec["tier_used"] = ",".join(sorted(set(tried)))
            rec["reason"] = (
                f"searched via {rec['tier_used'] or 'none'}; no tier met thresholds "
                f"(min_coverage={min_cov:g}, min_identity={min_id:g})"
            )
            if mstat == "partial":
                rec["reason"] += " [marker only partially built - absence is weaker evidence]"

        counts[rec["state"]] += 1
        results.append(rec)

    cols = list(results[0]) if results else ["sample", "marker_id", "state", "reason"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow(r)

    sys.stderr.write(
        f"[marker] {args.sample}: {counts[PRESENT]} present, {counts[ABSENT]} absent, "
        f"{counts[UNKNOWN]} unknown (of {len(results)} catalogue markers) -> {args.out}\n"
    )
    for r in results:
        if r["state"] == PRESENT:
            sys.stderr.write(f"[marker]   PRESENT {r['marker_id']} ({r['taxon']}): {r['reason']}\n")
    if counts[UNKNOWN]:
        sys.stderr.write(
            f"[marker] {counts[UNKNOWN]} marker(s) UNKNOWN. Rules depending on these "
            "return indeterminate and will HOLD findings at 'supported' rather than "
            "confirming them.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
