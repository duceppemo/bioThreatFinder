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

The manifest is authoritative PER TIER as well, via `n_nt_seqs`, `n_aa_seqs` and
`has_hmm`. `status: partial` says some tier is missing; those three say which.
Without them a marker that declares `dna+protein+hmm` and has no HMM still
reported `tier_used = dna,hmm,protein`, claiming a search that never ran. A
declared-but-unbuilt tier is now excluded from `tier_used`, named in
`tiers_not_built` and in the reason, and - where it is the marker's ONLY tier -
forces UNKNOWN rather than ABSENT. An older manifest lacking the three columns
keeps the previous behaviour rather than reading as all-gaps.

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
import math
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
    """REFERENCE NAME -> sequence length, from a combined marker FASTA.

    PER SEQUENCE, NOT SUMMED PER MARKER. This returned `marker_id -> summed
    length`, and aa_evidence() divided covered residues by that sum - the exact
    defect nt_evidence() documents fixing, never carried across:

        "BREADTH IS PER SEQUENCE, NOT OVER THEIR SUMMED LENGTH. A marker often
         carries several ALTERNATIVE representations of one locus ... Dividing
         covered bases by the summed 1878 bp caps breadth at 0.5 even when one
         allele is covered end to end, so a min_coverage of 0.80 is structurally
         unreachable at ANY depth."

    Measured on the shipped pack before this fix: 107 of 125 protein markers
    could not reach their own `min_coverage` on a PERFECT full-length hit. An
    8-allele marker caps at 0.125 against floors of 0.50-0.60, so `ebola_NP`,
    `nipah_L`, `marburg_NP`, `lassa_NP` and most of the other divergent-agent
    markers wrote a confident ABSENT - which feeds `absent()` as TRUE and can
    fire a veto_if or vaccine_if rule.

    The key is the first whitespace-delimited token of the header, which is what
    MMseqs2/DIAMOND report in the m8 target column.
    """
    out: dict[str, int] = {}
    if not path or not path.exists():
        return {}
    current = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                current = name or None
                if current:
                    out.setdefault(current, 0)
            elif current:
                out[current] += len(line.strip())
    return out


# ---------------------------------------------------------------------------
# Tier: nucleotide read alignment
# ---------------------------------------------------------------------------

def expected_breadth_at(depth: float) -> float:
    """Breadth a locus can reach at this depth: Lander-Waterman, 1 - e^-d.

    THE REASON A FIXED BREADTH FLOOR IS NOT DEPTH-NEUTRAL. A catalogue min_coverage of
    0.80 silently requires depth >= -ln(0.20) = 1.6x on the locus. Below that the
    threshold is unreachable no matter how real the organism is, and the marker is
    recorded ABSENT - which is not a missing call but positive evidence against a true
    finding, and it propagates into support_requires and the tier.

    Measured on the plumbing positive control: Y. pestis at 0.34x mean depth, marker
    ypest_pmt1_caf1 with 13 reads at 99.8% identity, observed breadth 0.323 against an
    expectation of 1 - e^-0.34 = 0.29. Fully covered for the depth available, and called
    ABSENT against a 0.80 floor - so the four plasmid markers all failed, support_requires
    evaluated FALSE, and a correct Y. pestis detection was held at `screened_only` on the
    one scenario whose file says any outcome other than confirmed is a configuration fault.
    """
    if depth <= 0:
        return 0.0
    return 1.0 - math.exp(-depth)


def nt_evidence(bam: Path, mapq_min: int) -> dict[str, dict]:
    """Per-marker nucleotide evidence, judged on the BEST single reference sequence.

    TWO THINGS THIS GETS RIGHT THAT THE FIRST VERSION DID NOT, both found by running the
    plumbing positive control at 6.8x depth and still failing to confirm Y. pestis.

    (1) BREADTH IS PER SEQUENCE, NOT OVER THEIR SUMMED LENGTH. A marker often carries
        several ALTERNATIVE representations of one locus - ypest_ppcp1_pla holds two
        939 bp alleles, WP_002218234.1 and WP_002231174.1. Dividing covered bases by the
        summed 1878 bp caps breadth at 0.5 even when one allele is covered end to end, so
        a min_coverage of 0.80 is structurally unreachable at ANY depth. A marker is
        present if ANY of its representatives is well covered, so the best one decides.

    (2) THE MAPQ FLOOR APPLIES TO EVERY MARKER. An earlier version of this dropped the
        floor for multi-sequence markers, reasoning that MAPQ 0 there means "ambiguous
        between two alleles of this marker" and that the identity requirement would keep
        a mismapped read out. That reasoning is WRONG, and the run proved it: identity
        does not discriminate when the other source is a >99% identical homologue.
        ypseudo_chr_inv went from 2 reads to 271 at 99.4% identity in a sample containing
        no Y. pseudotuberculosis at all - those were Y. pestis reads on the homologous
        inv locus - and the false PRESENT vetoed a correct select-agent detection, the
        worst failure direction this tool has.
        The redundancy that motivated the relaxation is now removed at its source:
        build_marker_db.py de-duplicates near-identical sequences within a marker when it
        builds the combined index, so reads map uniquely and MAPQ is high without any
        relaxation here. Fix the cause, not the symptom.
    """
    ref_len: dict[str, int] = {}
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
                ref_len[name] = ln

    if not ref_to_marker:
        sys.stderr.write(
            f"[marker] WARNING: no reference in {bam.name} uses the "
            "'marker|<id>|...' convention; nucleotide tier yields nothing\n"
        )
        return {}

    n_seqs: dict[str, int] = defaultdict(int)
    for mid in ref_to_marker.values():
        n_seqs[mid] += 1

    # Per REFERENCE, so the best representative can be chosen later.
    agg = defaultdict(lambda: {"reads": 0, "nm": 0, "aln": 0})
    for line in stream(["samtools", "view", "-F", "0x904", str(bam)]):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        ref = f[2]
        mid = ref_to_marker.get(ref)
        if not mid:
            continue
        try:
            mapq = int(f[4])
        except ValueError:
            mapq = 0
        if mapq < mapq_min:
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
        rec = agg[ref]
        rec["reads"] += 1
        rec["nm"] += nm
        rec["aln"] += aln

    # Covered bases per reference. -Q 0 for multi-sequence markers, filtered otherwise;
    # one pass each, and the BAM is small.
    covered: dict[str, int] = defaultdict(int)
    # -G 0x800 EXCLUDES SUPPLEMENTARY, matching the `view -F 0x904` that counts marker
    # reads a few lines above. samtools depth's default filter drops
    # UNMAP/SECONDARY/QCFAIL/DUP but KEEPS supplementary, so a chimeric read's extra
    # segments lifted nt_breadth - the PRESENCE gate - while contributing nothing to
    # nt_reads. The same mismatch was fixed in coverage_stats.py (e69fd01); this file was
    # not in that review's scope. Direction of harm here is the worse one: an inflated
    # breadth makes a neighbour-specific marker read PRESENT, which fires veto_if and
    # vetoes a correct select-agent detection.
    for line in stream(["samtools", "depth", "-a", "-G", "0x800",
                        "-Q", str(mapq_min), str(bam)]):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        if f[0] in ref_to_marker and int(f[2]) > 0:
            covered[f[0]] += 1

    out: dict[str, dict] = {}
    for mid in sorted(n_seqs):
        refs = [r for r, m in ref_to_marker.items() if m == mid]
        best, best_frac = None, -1.0
        for r in refs:
            length = ref_len.get(r, 0)
            frac = (covered.get(r, 0) / length) if length else 0.0
            if frac > best_frac:
                best, best_frac = r, frac
        length = ref_len.get(best, 0) if best else 0
        rec = agg.get(best, {"reads": 0, "nm": 0, "aln": 0}) if best else {
            "reads": 0, "nm": 0, "aln": 0}
        # Reads and identity are summed over the whole marker: a read on any allele is
        # evidence for the marker. Breadth and depth are the BEST single sequence.
        tot = {"reads": 0, "nm": 0, "aln": 0}
        for r in refs:
            a = agg.get(r)
            if a:
                for k in tot:
                    tot[k] += a[k]
        depth = (rec["aln"] / length) if length else 0.0
        out[mid] = {
            "nt_reads": tot["reads"],
            "nt_breadth": max(0.0, best_frac),
            "nt_identity": (1.0 - tot["nm"] / tot["aln"]) if tot["aln"] else 0.0,
            "nt_length": length,
            "nt_depth": depth,
            "nt_expected_breadth": expected_breadth_at(depth),
            "nt_n_seqs": n_seqs[mid],
            "nt_best_ref": (best or ""),
        }
    return out


def aa_evidence(hits: Path, aa_lengths: dict[str, int]) -> dict[str, dict]:
    """m8: query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bits

    The marker protein is the TARGET; reads are queries. Coverage is the union of
    target intervals, so many reads stacked on one domain do not read as
    full-length coverage.

    COVERAGE IS PER TARGET SEQUENCE AND THE BEST ONE DECIDES, mirroring
    nt_evidence(). A marker holds ALTERNATIVE representations of one locus; a
    read matching any of them is evidence for the marker, and dividing by their
    summed length made the gate unreachable (see fasta_lengths()).

    Identity, bits and e-value are taken from the best-scoring hit ON THE CHOSEN
    SEQUENCE, so the two numbers the gate compares describe the same protein.
    They previously came from the best hit anywhere in the marker while coverage
    was unioned over every allele, so a marker could pass identity on one allele
    using coverage accumulated from another. `aa_hits` stays a whole-marker
    count, because it reports how much evidence there is, not where it landed.
    """
    # Per TARGET, so the best representative can be chosen below.
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    best_on: dict[str, dict] = {}
    n_hits: dict[str, int] = defaultdict(int)
    unknown_targets: set[str] = set()
    with hits.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            target = f[1]
            mid = marker_of(target)
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
            if aa_lengths and target not in aa_lengths:
                unknown_targets.add(target)
            intervals[target].append((tstart, tend))
            n_hits[mid] += 1
            cur = best_on.get(target)
            if cur is None or bits > cur["aa_bits"]:
                best_on[target] = {"aa_identity": fident, "aa_bits": bits,
                                   "aa_evalue": evalue}

    # A TARGET THAT IS IN THE HITS BUT NOT IN THE FASTA HAS NO LENGTH, so its
    # coverage silently computes as 0 and the marker reads ABSENT. nt_evidence()
    # warns on the same mismatch; this tier did not, so an aligner run against a
    # different or stale marker FASTA produced confident negatives in silence.
    if unknown_targets:
        sample = ", ".join(sorted(unknown_targets)[:3])
        sys.stderr.write(
            f"[marker] WARNING: {len(unknown_targets)} protein target(s) in "
            f"{hits.name} are absent from the marker FASTA (e.g. {sample}); "
            "their coverage cannot be computed and those markers will read "
            "ABSENT. Is --aa-fasta the same build the search used?\n")

    by_marker: dict[str, list[str]] = defaultdict(list)
    for target in intervals:
        by_marker[marker_of(target)].append(target)

    n_seqs: dict[str, int] = defaultdict(int)
    for name in aa_lengths:
        mid = marker_of(name)
        if mid:
            n_seqs[mid] += 1

    out = {}
    for mid, targets in by_marker.items():
        best, best_frac = None, -1.0
        for t in sorted(targets):
            length = aa_lengths.get(t, 0)
            frac = (union_length(intervals[t]) / length) if length else 0.0
            if frac > best_frac:
                best, best_frac = t, frac
        hit = dict(best_on.get(best, {"aa_identity": 0.0, "aa_bits": 0.0,
                                      "aa_evalue": 1.0}))
        hit["aa_hits"] = n_hits[mid]
        hit["aa_coverage"] = max(0.0, best_frac)
        hit["aa_length"] = aa_lengths.get(best, 0) if best else 0
        hit["aa_n_seqs"] = n_seqs.get(mid, 0)
        hit["aa_best_ref"] = best or ""
        out[mid] = hit
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
            # THE PROFILE LENGTH, because without it `hmm_bits` cannot be interpreted.
            # A bitscore is roughly linear in profile length, so 94 bits on marburg_L's
            # ~2,100 aa L polymerase and 106 bits on ebola_GP's ~676 aa glycoprotein are
            # not comparable numbers - per position they differ by more than 3x, and
            # that ratio is the only thing that distinguishes them. It was computed here
            # and discarded, so the question "is this hit strong for this profile?"
            # could not be asked of marker_states.tsv at all. D23.
            "hmm_qlen": length,
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
    ap.add_argument("--platform", default="",
                    help="illumina|ont. The catalogue's per-marker min_identity is an "
                         "Illumina number (76 markers sit at 0.95 and 21 at 0.97); ONT "
                         "reads cannot meet those on sequencing error alone, so a "
                         "genuine marker would read ABSENT on that platform - positive "
                         "evidence against a real finding. See --ont-identity-slack.")
    ap.add_argument("--min-low-depth-reads", type=int, default=5,
                    help="Minimum reads for the depth-limited presence route, used only "
                         "when the marker's achievable breadth is below its min_coverage. "
                         "Identity is never relaxed. [5]")
    ap.add_argument("--low-depth-breadth-fraction", type=float, default=0.60,
                    help="On the depth-limited route, observed breadth must be at least "
                         "this fraction of what the depth allows (1-e^-d). [0.60]")
    ap.add_argument("--ont-identity-slack", type=float, default=0.05,
                    help="Subtracted from each marker's nucleotide min_identity when "
                         "--platform ont, floored at 0.70. PROVISIONAL: reasoned from "
                         "R10.4.1 read accuracy (~97-98%%) and badread's simulated 95%%, "
                         "not measured on this pipeline. Set from the benchmark once it "
                         "has run. [0.05]")
    ap.add_argument("--hmm-min-bits", type=float, default=25.0,
                    help="minimum domain bitscore for the HMM tier [25]")
    ap.add_argument("--protein-identity-slack", type=float, default=0.0,
                    help="subtracted from the catalogue's min_identity for the PROTEIN "
                         "tier. Default 0: honour what the curator wrote. This was an "
                         "undeclared blanket 0.15, which made six housekeeping genes "
                         "PRESENT off unrelated background on a single-organism sample "
                         "(D23).")
    ap.add_argument("--ont-protein-identity-slack", type=float, default=0.05,
                    help="protein-tier slack on ONT, where translated identity still "
                         "carries some basecalling error - but far less than nucleotide, "
                         "so this is deliberately smaller than --ont-identity-slack and "
                         "does NOT compound with it [0.05]")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    catalogue = read_tsv(args.markers)
    if not catalogue:
        sys.exit(f"ERROR: no markers read from {args.markers}")

    # Manifest status decides what CAN be evaluated at all.
    status: dict[str, str] = {}
    manifest_note: dict[str, str] = {}
    # PER-MARKER TIER AVAILABILITY, which is a different question from the
    # run-level `tiers_supplied` below.
    #
    # `tiers_supplied` answers "did this run receive an input file for the tier".
    # It cannot answer "does this marker have anything in that file", and that gap
    # was being papered over: `tried` was appended unconditionally inside each
    # tier's branch, so a marker whose HMM was never built still reported
    # `tier_used = dna,hmm,protein` and a reason reading "searched via
    # dna,hmm,protein; no tier met thresholds". The row asserted that a tier looked
    # and found nothing when the model did not exist. On the shipped pack that is
    # 12 of 37 hmm-declaring markers - cbut_bont, ebola_NP, marburg_GP, nipah_L and
    # the rest - every one of them a divergence-catching marker on a
    # high-consequence agent, which is exactly the set where markers.tsv argues the
    # HMM tier is the one that matters.
    built_tiers: dict[str, set[str]] = {}
    manifest_has_tier_cols = False
    if args.manifest:
        rows = read_tsv(args.manifest)
        if not rows:
            sys.stderr.write(
                f"[marker] WARNING: manifest {args.manifest} is empty or unreadable; "
                "every marker will be reported UNKNOWN\n"
            )
        for row in rows:
            mid = row.get("marker_id")
            if not mid:
                continue
            status[mid] = (row.get("status") or "").strip().lower()
            manifest_note[mid] = row.get("note", "")
            # An older manifest predating these columns must not silently turn
            # every marker into a gap, so a column that is absent from the HEADER
            # means "not recorded, assume built" - the behaviour from before this
            # was checked at all. A column that is present but empty means no.
            have = set()
            for tier, col in (("dna", "n_nt_seqs"), ("protein", "n_aa_seqs")):
                if col in row:
                    manifest_has_tier_cols = True
                    try:
                        if int(float((row.get(col) or "0").strip() or 0)) > 0:
                            have.add(tier)
                    except ValueError:
                        pass
            if "has_hmm" in row:
                manifest_has_tier_cols = True
                if (row.get("has_hmm") or "").strip().lower() in ("yes", "true", "1"):
                    have.add("hmm")
            built_tiers[mid] = have
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
        # PLATFORM, for the same reason score_evidence.py needs it. A 0.95 or 0.97
        # nucleotide floor is a statement about Illumina error as much as about
        # sequence divergence, and ONT reads do not reach it; applying it unchanged
        # made a present marker read ABSENT, which is not a missing call but positive
        # evidence against a real finding. Only the nucleotide tier is adjusted - the
        # protein floor is already derived from min_id with its own allowance, and
        # translated identity is far less sensitive to basecalling error.
        id_note = ""
        # The catalogue value is KEPT, because the protein floor below is derived from it
        # and the comment above promises only the nucleotide tier moves. It did not: the
        # relaxation assigned to min_id in place and `aa_min_id = max(0.50, min_id - 0.15)`
        # then read the relaxed number, so an ONT run at a catalogue floor of 0.95 got a
        # protein floor of 0.75 instead of the documented 0.80 - widening a false-PRESENT
        # route on the tier that has no bitscore or coverage-independent floor of its own,
        # and false PRESENT on a marker is what fires a wrong veto or vaccine call.
        min_id_catalogue = min_id
        if (args.platform or "").strip().lower() == "ont" and min_id > 0:
            relaxed = max(0.70, min_id - max(0.0, args.ont_identity_slack))
            if relaxed < min_id:
                id_note = f" [ont: floor {min_id:g} -> {relaxed:g}, nucleotide tier only]"
                min_id = relaxed

        rec = {
            "sample": args.sample, "marker_id": mid,
            "taxon": row.get("taxon", ""), "gene": row.get("gene", ""),
            "marker_class": row.get("marker_class", ""),
            "copy_expectation": row.get("copy_expectation", ""),
            "state": UNKNOWN, "tier_used": "", "reason": "",
            # Declared in markers.tsv but absent from the marker database for this
            # marker. Initialised here so the column exists even on rows that
            # return at gate 1, since the writer takes its fieldnames from the
            # first record.
            "tiers_not_built": "",
            "nt_reads": "", "nt_breadth": "", "nt_identity": "",
            "aa_hits": "", "aa_coverage": "", "aa_identity": "", "aa_bits": "",
            "hmm_domains": "", "hmm_coverage": "", "hmm_bits": "", "hmm_qlen": "",
            # Emitted for the same reason as hmm_qlen: `aa_coverage` is a ratio whose
            # denominator was thrown away, so an aa_coverage of 0 could not be told
            # apart from a marker whose FASTA header did not parse and whose coverage
            # gate was therefore SKIPPED entirely (the aa_length == 0 bypass below).
            "aa_length": "",
            # The denominator's PROVENANCE. `aa_length` is now one sequence's
            # length, not the marker's summed length, so without these a reader
            # cannot tell which of several alleles the coverage describes - the
            # same gap `aa_length` itself was added to close.
            #
            # `aa_n_seqs` COUNTS THE SEARCH INDEX, NOT THE MANIFEST. The manifest's
            # `n_aa_seqs` counts the per-marker archive; concat_files() collapses
            # near-identical alleles when it builds markers.aa.faa, so the two
            # disagree by design and the gap is large (ebola_NP: 8 archived, 2
            # searched). This column is the one that describes what was actually
            # searched, which is what a coverage denominator needs.
            "aa_n_seqs": "", "aa_best_ref": "",
            # Carried so score_evidence.py can compute, from the target's observed
            # depth, whether a marker that was not found would have been SAMPLED at
            # all. Without the length, "searched and not found" cannot be told apart
            # from "never had a chance of being seen".
            "nt_length": "",
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
        supplied = {t for t in requested if tiers_supplied.get(t)}

        # Gate 2b: and does THIS marker actually have that tier built? A tier the
        # marker does not have was never searched, so it must not be reported as
        # having been searched, and its silence is not evidence of absence.
        unbuilt: set[str] = set()
        if args.manifest and manifest_has_tier_cols:
            unbuilt = {t for t in requested if t not in built_tiers.get(mid, set())}
        rec["tiers_not_built"] = ",".join(sorted(unbuilt))

        evaluable = supplied - unbuilt
        if not evaluable:
            if supplied and unbuilt >= supplied:
                # Distinct from "no input supplied": the inputs arrived, this marker
                # just has nothing in them. UNKNOWN either way - nothing was
                # searched, so absence cannot be asserted - but the reason has to
                # say which, or a database gap reads as a run configuration error.
                rec["reason"] = (
                    f"requested tier(s) {'+'.join(sorted(requested))} but none is "
                    f"built for this marker (not built: {','.join(sorted(unbuilt))}); "
                    f"nothing was searched"
                )
            else:
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
                            for k, v in e.items()})
                tried.append("dna")
                # AT LEAST ONE READ, WHATEVER THE THRESHOLDS SAY.
                #
                # nt_evidence() emits a record for every marker in the BAM header,
                # including those with no alignments at all (nt_breadth 0.0,
                # nt_identity 0.0). With min_coverage and min_identity both 0.00 -
                # which bmallei_deletions carries - the gate passed on zero evidence
                # and the marker was reported PRESENT. Any marker with a 0 floor
                # became a permanent free PRESENT the moment a rule referenced it.
                # DEPTH-AWARE PRESENCE.
                #
                # The strict floor still applies wherever the depth can actually reach
                # it. Where it cannot - expected breadth below the floor - a marker may
                # instead qualify on being as covered as the depth permits, backed by a
                # minimum read count and the full identity requirement. Identity is
                # never relaxed: that is what distinguishes a real locus from a
                # mismapped read, and it is the discriminator that survives low depth.
                exp_b = e.get("nt_expected_breadth", 0.0)
                depth_limited = exp_b < min_cov
                strict_ok = e["nt_breadth"] >= min_cov
                low_depth_ok = (
                    depth_limited
                    and e["nt_reads"] >= args.min_low_depth_reads
                    and exp_b > 0
                    and e["nt_breadth"] >= args.low_depth_breadth_fraction * exp_b)
                if (e["nt_reads"] >= 1
                        and (strict_ok or low_depth_ok)
                        and e["nt_identity"] >= min_id):
                    how = "" if strict_ok else (
                        f", DEPTH-LIMITED: {e['nt_depth']:.2f}x on this marker allows "
                        f"at most {exp_b:.3f} breadth vs a {min_cov:g} floor, "
                        f"observed {e['nt_breadth']:.3f}")
                    passed.append(
                        f"dna(breadth={e['nt_breadth']:.3f}>={min_cov:g}, "
                        f"id={e['nt_identity']:.3f}>={min_id:g}, "
                        f"reads={e['nt_reads']}){id_note}{how}")
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
                    "aa_length": e.get("aa_length", ""),
                    "aa_n_seqs": e.get("aa_n_seqs", ""),
                    "aa_best_ref": e.get("aa_best_ref", ""),
                })
                tried.append("protein")
                # THE CATALOGUE'S FLOOR, NOT A BLANKET SUBTRACTION FROM IT.
                #
                # This was `max(0.50, min_id_catalogue - 0.15)`: every marker lost 15
                # identity points regardless of what the curator wrote, on the reasoning
                # that "protein floors are lower than nucleotide by nature". The premise
                # is true and the implementation inverted the curation that encodes it.
                # The catalogue ALREADY carries the distinction: housekeeping genes that
                # cross-match across genera are floored at 0.95-0.98, and the divergent
                # virus block is floored at 0.78-0.90 on purpose, with a note saying
                # "identity floors here are deliberately LOW and the HMM tier carries
                # the detection". The blanket -0.15 overwrote exactly that.
                #
                # `[MEASURED]` on divergent_bundibugyo_filovirus_illumina_r1, a sample
                # containing one virus and no other watchlist pathogen, the old floor
                # made SIX markers PRESENT off unrelated background, every one a
                # conserved housekeeping gene: mtb_rpoB 0.877 (catalogue 0.98),
                # cdiph_rpoB 0.900 (0.97), bholm_recA 0.860 (0.97), cburn_icd 0.877
                # (0.95), rprow_gltA 0.833 (0.95), mtb_katG_inhA 0.840 (0.97). Honouring
                # the catalogue removes all six and keeps both real filovirus hits -
                # ebola_NP 0.880 and ebolaS_NP 1.000 against a deliberate 0.85 floor.
                # One of the six, mtb_katG_inhA, is `report_always`, so it was being
                # surfaced to analysts as isoniazid resistance on samples with no
                # mycobacteria in them.
                #
                # The slack is now DECLARED and defaults to zero. ONT keeps a small one
                # because translated identity still carries some basecalling error - but
                # far less than nucleotide, so it is not the nucleotide slack. Note this
                # reads min_id_catalogue, NOT min_id: the ONT nucleotide relaxation must
                # not compound with this one, which is the bug the comment above records.
                aa_slack = max(0.0, args.protein_identity_slack)
                if (args.platform or "").strip().lower() == "ont":
                    aa_slack = max(aa_slack, max(0.0, args.ont_protein_identity_slack))
                aa_min_id = max(0.50, min_id_catalogue - aa_slack)
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
            if unbuilt:
                rec["reason"] += (
                    f"; NOT SEARCHED - declared but not built for this marker: "
                    f"{','.join(sorted(unbuilt))}"
                )
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
            if unbuilt:
                # Named, not merely omitted. "searched via dna,protein" is true but
                # incomplete; a reader weighing an ABSENT on a divergent agent needs
                # to know the tier built to catch divergence never ran.
                rec["reason"] += (
                    f"; NOT SEARCHED - declared but not built for this marker: "
                    f"{','.join(sorted(unbuilt))}"
                )
            if mstat == "partial":
                rec["reason"] += " [marker only partially built - absence is weaker evidence]"

        counts[rec["state"]] += 1
        results.append(rec)

    # FIELDNAMES FROM THE UNION OF EVERY RECORD, not from results[0]. The keys on a
    # record depend on how far it got: a row that returns at gate 1 never reaches the
    # nt measurement block, so it carries no `nt_depth`, `nt_expected_breadth`,
    # `nt_n_seqs` or `nt_best_ref`. Taking the header from the FIRST record alone meant
    # the whole file's schema was decided by whichever catalogue row happened to sort
    # first, and `extrasaction="ignore"` then dropped those four columns from every
    # other row WITHOUT error. Measured 2026-10-20: a states file whose first marker was
    # UNKNOWN lost all four, including on the PRESENT rows that had measured them.
    # nt_expected_breadth is the input to the depth gate - the rule that marker absence
    # is only called when depth could have shown it - so the column that exists to stop
    # a false ABSENT can disappear silently, which is the failure this file is built to
    # prevent. Some keys are pre-initialised above for this very reason; that fix was
    # per-column and could not cover a block that adds four at once.
    cols: list[str] = []
    _seen: set[str] = set()
    for _r in results:
        for _k in _r:
            if _k not in _seen:
                _seen.add(_k)
                cols.append(_k)
    if not cols:
        cols = ["sample", "marker_id", "state", "reason"]
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

    # A DATABASE GAP IS NOT A RUN EVENT, so it is reported once per run rather than
    # per marker. These rows are still ABSENT or PRESENT on the tiers that exist -
    # the point is that one declared tier never ran, and on the shipped pack the
    # tier that is missing is almost always the HMM one, i.e. the tier carried
    # specifically to catch the divergent allele.
    gapped = [r for r in results if r.get("tiers_not_built")]
    if gapped:
        sys.stderr.write(
            f"[marker] {len(gapped)} marker(s) declare a tier that is not built in "
            f"this database; it was NOT searched and its silence is not evidence:\n"
        )
        for r in gapped:
            sys.stderr.write(
                f"[marker]   GAP {r['marker_id']} ({r['taxon']}): "
                f"missing {r['tiers_not_built']}, state {r['state']}\n"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
