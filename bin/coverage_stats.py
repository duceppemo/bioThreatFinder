#!/usr/bin/env python3
"""
Per-taxon alignment evidence from a competitive-mapping BAM.

This is the step that separates a detection from a k-mer artifact.

The BAM is produced by mapping reads against a bundle containing the target
genome AND its near neighbours, with every reference renamed
``taxid|<taxid>|<original header>``. Because minimap2 emits one primary
alignment per read at its best locus, restricting to primary alignments gives
competitive assignment for free: a read from an environmental relative lands on
the relative, not on the high-consequence target.

Metrics, and why each one is here
---------------------------------
unique_reads     Primary, MAPQ-filtered reads assigned to the target. MAPQ
                 matters: a read mapping equally well to target and neighbour
                 gets MAPQ 0 from minimap2 and must not count as evidence.

breadth          Fraction of target bases at depth >= 1. The single most useful
                 metric. A true low-level organism spreads thin across the whole
                 genome; an artifact stacks on one conserved gene.

window_occupancy Fraction of fixed-width windows containing >= 1 aligned base.
                 Breadth's more robust sibling. At very low depth, breadth is
                 small for both real and artifactual signal, but occupancy still
                 separates them: 40 reads scattered over 40 windows is an
                 organism, 40 reads in 2 windows is a repeat.

depth_cv         Coefficient of variation of per-base depth. High CV = piled up.
                 Read as a companion to occupancy, not alone: genuinely uneven
                 real coverage happens (rRNA operons, prophage, plasmids).

mean_identity    From NM over aligned length. A divergent-but-real relative
                 shows lower identity; so does a forced alignment of a
                 mis-assigned read. Interpreted per taxon - some watchlist
                 viruses are legitimately far from any available reference.

neighbour_ratio  unique_reads(target) / unique_reads(best neighbour). Below ~1
                 the evidence points at the neighbour, not the target.

Needs only `samtools` on PATH. No third-party Python.
"""

from __future__ import annotations

import argparse
import re
import os
import math
import statistics
import tempfile
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
# Reference-consuming CIGAR ops, for aligned reference span.
REF_OPS = {"M", "D", "N", "=", "X"}
# Query-aligned ops, for identity denominator.
ALN_OPS = {"M", "I", "=", "X"}


def run(cmd: list[str]) -> str:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"ERROR: command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr}")
    return proc.stdout


def stream(cmd: list[str]):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            bufsize=1 << 20)
    assert proc.stdout is not None
    for line in proc.stdout:
        yield line
    proc.stdout.close()
    rc = proc.wait()
    if rc != 0:
        err = proc.stderr.read() if proc.stderr else ""
        sys.exit(f"ERROR: command failed ({rc}): {' '.join(cmd)}\n{err}")


def load_mask(path: Path | None) -> dict[str, list[tuple[int, int]]]:
    """contig -> merged, sorted masked intervals, from a BED file or a directory of them.

    WHY MASKED REGIONS ARE HANDLED HERE RATHER THAN BY DELETING SEQUENCE

    scripts/mask_pack_refs.py soft-masks: it records intervals that are vector,
    host-derived, low-complexity, or shared near-identically with a phylogenetically
    distant reference (assembly cross-contamination). Hard-masking those to N would
    make the failure mode invisible - the reads would simply stop aligning and the
    taxon would quietly vanish with no record of why.

    Handling it here instead gives the auditable version:

      * masked bases are excluded from the breadth and occupancy DENOMINATORS, so a
        genome is judged on the portion of it we trust;
      * reads lying ENTIRELY inside a masked interval do not count toward
        unique_reads, because that is untrustworthy evidence;
      * but they are reported as reads_in_masked_regions, so the situation is visible.

    That combination is what catches the Wright et al. 2026 case directly: an
    H. pylori call whose support turned out to be contaminating phage sequence inside
    the reference genome now shows unique_reads near zero AND a high
    reads_in_masked_regions, instead of a confident false positive.
    """
    if not path:
        return {}
    beds: list[Path] = []
    if path.is_dir():
        beds = sorted(path.glob("*.mask.bed"))
        if not beds:
            sys.stderr.write(f"[coverage] WARNING: no *.mask.bed under {path}\n")
    elif path.exists():
        beds = [path]
    else:
        sys.stderr.write(f"[coverage] WARNING: mask path {path} not found; "
                         "proceeding UNMASKED\n")
        return {}

    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for bed in beds:
        with bed.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                try:
                    raw[f[0]].append((int(f[1]), int(f[2])))
                except ValueError:
                    continue
    out: dict[str, list[tuple[int, int]]] = {}
    for contig, ivs in raw.items():
        ivs.sort()
        merged: list[list[int]] = []
        for a, b in ivs:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        out[contig] = [(a, b) for a, b in merged]
    return out


def contig_key(refname: str) -> str:
    """`taxid|632|NZ_CP009973.1` -> `NZ_CP009973.1`.

    Competitive bundles prefix every reference with its taxid, but mask BEDs are
    written against the original contig names, so the prefix has to come off before
    lookup. Silently failing to match would leave everything unmasked while
    appearing to work.
    """
    if refname.startswith("taxid|"):
        parts = refname.split("|", 2)
        if len(parts) >= 3:
            return parts[2].split()[0]
    return refname.split()[0]


def in_masked(ivs: list[tuple[int, int]], start: int, end: int) -> bool:
    """True if [start, end) lies entirely within one masked interval."""
    lo, hi = 0, len(ivs) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        a, b = ivs[mid]
        if start < a:
            hi = mid - 1
        elif start >= b:
            lo = mid + 1
        else:
            return end <= b
    return False


def masked_len_within(ivs: list[tuple[int, int]], length: int) -> int:
    return sum(min(b, length) - a for a, b in ivs if a < length)


def taxid_of(refname: str) -> str | None:
    """`taxid|632|NZ_CP009973.1 Yersinia...` -> '632'."""
    if not refname.startswith("taxid|"):
        return None
    parts = refname.split("|", 2)
    return parts[1] if len(parts) >= 2 and parts[1].isdigit() else None


def read_header(bam: Path) -> tuple[dict[str, int], dict[str, list[str]]]:
    """Return (refname -> length, taxid -> [refnames])."""
    lengths: dict[str, int] = {}
    by_taxid: dict[str, list[str]] = defaultdict(list)
    for line in stream(["samtools", "view", "-H", str(bam)]):
        if not line.startswith("@SQ"):
            continue
        name = length = None
        for field in line.rstrip("\n").split("\t")[1:]:
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                length = int(field[3:])
        if name is not None and length is not None:
            lengths[name] = length
            tid = taxid_of(name)
            if tid:
                by_taxid[tid].append(name)
    if not lengths:
        sys.exit(f"ERROR: no @SQ lines in {bam}; is it a valid BAM with a header?")
    if not by_taxid:
        sys.exit(
            f"ERROR: no reference in {bam} uses the 'taxid|<id>|...' naming convention.\n"
            "Competitive attribution needs it. Rebuild the bundle with "
            "scripts/build_pack_refs.py."
        )
    return lengths, dict(by_taxid)


def cigar_spans(cigar: str) -> tuple[int, int]:
    """(reference span, query-aligned length)."""
    ref_span = aln_len = 0
    for num, op in CIGAR_RE.findall(cigar):
        n = int(num)
        if op in REF_OPS:
            ref_span += n
        if op in ALN_OPS:
            aln_len += n
    return ref_span, aln_len


def collect_alignments(bam: Path, mapq_min: int, mask=None):
    """Primary, non-supplementary, mapped alignments grouped by taxid.

    The extra counters here mirror the per-alignment metric set UKHSA's mSCAPE
    programme stores in production (evenness, coverage_1x/10x, mapped vs uniquely
    mapped reads, mean_read_identity, read_duplication_rate, forward_proportion,
    mean_alignment_length). Two of them catch artifacts the others miss:

      forward_proportion    a genuine shotgun alignment is near 0.5. A strong strand
                            imbalance points at a single amplicon, an adapter
                            artifact or a repeat, not an organism.
      duplication_rate      reads stacked at identical start positions inflate depth
                            without adding evidence. Without an upstream markdup
                            step this is a PROXY - the fraction of primary reads
                            sharing a start position and strand - so it reads high
                            for genuinely deep, narrow coverage too. Interpret it
                            with breadth, not alone.
    """
    mask = mask or {}
    per_taxid = defaultdict(lambda: {"reads": 0, "nm": 0, "aln": 0, "bases": 0,
                                     "mapped_any_mapq": 0, "forward": 0,
                                     "masked_reads": 0, "starts": None})
    # 0x904 = unmapped + secondary + supplementary. Excluding these leaves one
    # record per read at its best locus, which is the competitive assignment.
    # No -q here: we want BOTH the MAPQ-filtered count (evidence) and the raw mapped
    # count (denominator), because their ratio is itself informative - a taxon whose
    # reads nearly all fall below the MAPQ floor is one whose reads map equally well
    # somewhere else.
    cmd = ["samtools", "view", "-F", "0x904", str(bam)]
    for line in stream(cmd):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        tid = taxid_of(f[2])
        if tid is None:
            continue
        try:
            flag = int(f[1])
            mapq = int(f[4])
        except ValueError:
            continue
        rec0 = per_taxid[tid]
        rec0["mapped_any_mapq"] += 1
        if rec0["starts"] is None:
            rec0["starts"] = defaultdict(int)
        if mapq < mapq_min:
            continue
        ref_span, aln_len = cigar_spans(f[5])

        # A read lying ENTIRELY inside a masked interval is not evidence. It is
        # counted and reported, but excluded from unique_reads - otherwise a
        # reference contaminated with phage or vector sequence collects reads and
        # manufactures a confident false positive, which is exactly the documented
        # Wright et al. failure.
        ivs = mask.get(contig_key(f[2]))
        if ivs:
            try:
                start0 = int(f[3]) - 1
            except ValueError:
                start0 = -1
            if start0 >= 0 and in_masked(ivs, start0, start0 + ref_span):
                rec0["masked_reads"] += 1
                continue

        rec0["starts"][(f[2], f[3], flag & 0x10)] += 1
        if not (flag & 0x10):
            rec0["forward"] += 1
        nm = 0
        for tag in f[11:]:
            if tag.startswith("NM:i:"):
                try:
                    nm = int(tag[5:])
                except ValueError:
                    nm = 0
                break
        rec = per_taxid[tid]
        rec["reads"] += 1
        rec["nm"] += nm
        rec["aln"] += aln_len
        rec["bases"] += ref_span
    return per_taxid


def covered_by_taxid(bam: Path, mapq_min: int) -> dict[str, tuple[int, int, float]]:
    """Per-taxid (reference_bp, covered_bp, aligned_bp) in ONE pass, for the neighbour gate.

    D12 measured that `neighbour_ratio`'s denominator, for the case the gate exists to
    handle, is background cross-mapping and not a competing organism at all: in
    lod_gram_negative the declared neighbour holds ~2,700 reads of which ZERO are Yersinia.
    Raw read count cannot tell those apart from a real competitor - at 100 rpm the target
    had 2,336 reads and the absent neighbour 2,731 - but their SHAPE can. Cross-mapped
    reads pile onto the conserved fraction two genomes share; a present organism's reads
    spread over its genome.

        uniformity = (covered_bp / reference_bp) / min(1, aligned_bp / reference_bp)

    i.e. OBSERVED BREADTH DIVIDED BY THE BREADTH THOSE SAME BASES WOULD GIVE IF THEY WERE
    SPREAD UNIFORMLY. All three quantities come from one `samtools coverage` pass:
    reference_bp from endpos, covered_bp from covbases, aligned_bp from meandepth x endpos.
    No read length is assumed, which is what lets the same formula serve Illumina and ONT.

    THIS IS ALGEBRAICALLY THE PREVIOUS FORM, `covered_bp / min(reference_bp, span)`, AND
    THE CORRECTION THAT CLAIMED OTHERWISE WAS WRONG. `span` there is collect_alignments()'s
    summed ALIGNED reference span, i.e. the same normalising term, so the two differ only
    by indel and supplementary-filter handling: over all 467 neighbour rows on disk,
    max |old - new| = 0.0073, mean 0.00015. The old docstring's "the expectation needs no
    read-length assumption" was correct; it was misread as meaning there was no expectation.
    The present form is kept only because both terms now come from one samtools pass
    instead of two sources. It fixes nothing, because nothing was broken here.

    WHAT IS REAL IS THAT NORMALISING DOES NOT MAKE THIS ABUNDANCE-INDEPENDENT, and the
    docstring used to claim it did. Below ~1x depth the statistic is covered-per-aligned-
    base, which is ~1 for ANY taxon whose reads are too sparse to overlap - so a rare
    pedestal is indistinguishable from a present organism. Above ~1x it saturates to plain
    breadth, which for a sibling reference plateaus at the shared fraction. It is a
    DEPTH-GATED pile-up detector, informative only where there is enough depth for overlap,
    and it was described as a shape statistic that one threshold could serve at every level.

    Measured on ONE ABSENT neighbour
    (L. innocua, taxid 1642) across the Zymo libraries:

        14 reads 0.9676 | 21 reads 0.9234 | 1,083 reads 0.8513 | 11,593 reads 0.4176
        190,739 reads 0.2653 | 1,880,316 reads 0.3642 | 7,639,529 reads 0.3912

    A monotone decline across five orders of magnitude, crossing every threshold ever
    proposed for it. The same absent organism was "genuinely covered" when rare and a
    "pedestal" when abundant. See D20.

    WHAT THE RECALIBRATION FOUND (467 neighbour rows, all competitive BAMs on disk,
    truth labelled per specimen as PRESENT / ABSENT_SIBLING / ABSENT_DISTANT). Counting
    DISTINCT CONDITIONS rather than rows, because 288 of the 292 present rows are one
    simulated background replicated:

        PRESENT         9 conditions   0.524 - 0.991
        ABSENT_SIBLING 17 conditions   0.022 - 0.861
        ABSENT_DISTANT 26 conditions   0.000 - 0.604

    Three things follow. (1) Rows below ~0.05x neighbour depth carry no information and
    must not be judged - unfiltered, ABSENT_DISTANT reaches p90 0.932 purely from sparse
    non-overlapping reads. (2) Above that floor PRESENT and ABSENT_DISTANT nearly separate,
    but NOT cleanly: present E. coli in zymo_log scores 0.524 while absent Citrobacter in
    zymo_even ONT scores 0.604. (3) ABSENT_SIBLING spans the whole range and always will -
    a sibling reference IS covered by its relative's reads, which is the honest answer and
    the reason indistinguishable_complex exists.

    So this number can support "the neighbour is a pedestal" only with a depth
    precondition, and even then the present arm rests on ~5 independent conditions. It is
    reported for audit; params.contest_attribution is OFF and no gate reads it.

    SUPPLEMENTARY is in the filter list because samtools coverage does NOT exclude it by
    default, while the read counts this is compared against come from `view -F 0x904`.
    """
    out: dict[str, list[int]] = {}
    proc = subprocess.run(
        ["samtools", "coverage", "-q", str(mapq_min),
         "--ff", "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY", str(bam)],
        capture_output=True, text=True)
    if proc.returncode != 0:
        # Not fatal: the gate degrades to "no credibility information", which
        # score_evidence.py treats as the status quo rather than as a pass.
        sys.stderr.write(f"[coverage] WARNING: samtools coverage failed, neighbour "
                         f"credibility unavailable: {proc.stderr.strip()[:200]}\n")
        return {}
    for line in proc.stdout.splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 7:
            continue
        parts = f[0].split("|")
        if len(parts) < 2:
            continue
        rec = out.setdefault(parts[1], [0, 0, 0.0])
        try:
            ref_bp = int(f[2])
            rec[0] += ref_bp
            rec[1] += int(f[4])
            # meandepth is over the WHOLE reference, so meandepth x reference_bp is the
            # aligned base count - the normalising term, free in this same pass.
            rec[2] += float(f[6]) * ref_bp
        except ValueError:
            continue
    return {k: (v[0], v[1], v[2]) for k, v in out.items()}


def depth_profile(bam: Path, mapq_min: int, target_refs: set[str], window: int,
                  mask=None, lengths=None):
    """Per-base depth for the target's references only.

    `-aa`, NOT `-a`, AND THE DIFFERENCE IS ~1000x ON A FRAGMENTED GENOME.
    `samtools depth -a` emits positions only for references that received at least one
    counted read; a contig with zero coverage is omitted entirely. `total` is the
    denominator for breadth here and `window_possible` for occupancy, so on a multi-contig
    reference both were computed over only the contigs that recruited something. Measured
    on the T. gondii bundle of runs/lod_v4 (2,264 target contigs, 65,668,120 bp): `-a`
    emitted 57,470 positions - 0.087% of the reference - and the published row read
    breadth 0.469 / occupancy 0.976 / mean_depth 21.2x where the true values are
    0.000428 / 0.00125 / 0.0168x. So the two gates this module calls the most useful were
    never applied to a fragmented reference, and 177 of 363 pack references have more than
    one contig (Hammondia hammondi: 36,360). mean_depth crossing 20x also drove
    score_evidence.py's "genome_typing_supported" claim off a 0.017x genome.

    `-aa` alone would emit every neighbour base too, which is what the `-a` here was
    avoiding, so the target's contigs go to samtools as a BED and it does the restriction.
    Measured: 15.7 s for the 65.7 Mb target.

    Depth is accumulated with Welford rather than collected in a list, because one row per
    base is 65.7 M ints for that reference - about 2.6 GB.

    -G 0x800 excludes SUPPLEMENTARY alignments, matching collect_alignments()'s
    `view -F 0x904`. samtools depth's default filter drops UNMAP/SECONDARY/QCFAIL/DUP but
    KEEPS supplementary, so a chimeric ONT read's extra segments were counted toward
    breadth, coverage_10x and occupancy while contributing nothing to unique_reads.
    """
    mask = mask or {}
    lengths = lengths or {}
    covered = 0
    covered10 = 0
    total = 0
    masked_positions = 0
    _n = 0
    _mean = 0.0
    _m2 = 0.0
    window_hits: dict[tuple[str, int], bool] = {}
    window_possible: dict[tuple[str, int], bool] = {}
    cur_win = None
    cur_ref = None
    cur_ivs = None

    # samtools' BED reader splits on WHITESPACE, not tabs alone, so a reference whose
    # name contains a space cannot be addressed as a region at all - it parses as
    # name/start/end and fails. The SAM spec forbids spaces in @SQ SN for exactly this
    # reason and every reference this pack builds is space-free, but a hand-made or
    # third-party BAM can carry one, and the failure would otherwise be an opaque
    # "[bed_read] Parse error ... at line 1" from a subprocess.
    spacey = sorted(r for r in target_refs if any(c.isspace() for c in r))
    if spacey:
        sys.exit(f"ERROR: target reference name(s) contain whitespace: {spacey[:2]}. "
                 f"The SAM spec forbids spaces in @SQ SN, and samtools cannot address such "
                 f"a reference as a region, so per-base depth cannot be restricted to the "
                 f"target. Rebuild the bundle from FASTA headers whose first token is the "
                 f"whole name.")
    missing = sorted(r for r in target_refs if r not in lengths)
    if missing:
        sys.exit(f"ERROR: no header length for target reference(s) {missing[:3]}; without "
                 f"them the region BED cannot be written and `samtools depth -aa` would "
                 f"emit every neighbour base in the bundle.")

    bed_fh = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
    try:
        for r in sorted(target_refs):
            bed_fh.write(f"{r}\t0\t{lengths[r]}\n")
        bed_fh.close()
        cmd = ["samtools", "depth", "-aa", "-b", bed_fh.name,
               "-Q", str(mapq_min), "-G", "0x800", str(bam)]
        for line in stream(cmd):
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] not in target_refs:
                continue
            pos0 = int(f[1]) - 1
            # PER CONTIG, NOT PER BASE. `samtools depth` groups its output by reference, so
            # contig_key() and the mask lookup were being recomputed 65.7 M times for the
            # T. gondii reference to return the same answer. With masking on, that lookup -
            # not the samtools call, and not the window bookkeeping - was the dominant cost.
            if f[0] != cur_ref:
                cur_ref = f[0]
                cur_ivs = mask.get(contig_key(f[0]))
            ivs = cur_ivs
            if ivs and in_masked(ivs, pos0, pos0 + 1):
                masked_positions += 1
                continue
            d = int(f[2])
            total += 1
            _n += 1
            _delta = d - _mean
            _mean += _delta / _n
            _m2 += _delta * (d - _mean)
            # ONE DICT WRITE PER WINDOW, NOT PER BASE. `samtools depth` emits positions in
            # order, so the window key changes every `window` rows; writing it on every row
            # cost ~65.7 M dict operations on the T. gondii reference and dominated the
            # task (3m29s, against 7.3s for the samtools call feeding it). Measured before
            # changing anything: samtools coverage is ~1s and depth ~1-7s, so the loop was
            # the whole cost and the samtools calls were never the problem.
            win = (f[0], pos0 // window)
            if win != cur_win:
                cur_win = win
                window_possible[win] = True
            if d > 0:
                covered += 1
                if win not in window_hits:
                    window_hits[win] = True
            if d >= 10:
                covered10 += 1
    finally:
        os.unlink(bed_fh.name)

    mean_depth = _mean if _n else 0.0
    sd = math.sqrt(_m2 / _n) if _n else 0.0
    return (covered, total, (mean_depth, sd), len(window_hits), covered10,
            masked_positions, len(window_possible))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--taxid", required=True, help="target taxid (must appear in the bundle)")
    ap.add_argument("--taxon", default="", help="target taxon name, for the report")
    ap.add_argument("--group", default="", help="bacteria|virus|parasite|fungus")
    ap.add_argument("--role", default="sample", help="sample|control")
    ap.add_argument("--platform", default="",
                    help="illumina|ont. Recorded so score_evidence.py can apply a\nplatform-appropriate read-identity floor; a 0.98 floor built for Illumina cannot\nbe met by ONT reads at all.")
    ap.add_argument("--mapq-min", type=int, default=10)
    ap.add_argument("--window", type=int, default=1000,
                    help="window size in bp for occupancy [1000]")
    ap.add_argument("--total-reads", type=int, default=0,
                    help="host-depleted read count, for normalisation across samples")
    ap.add_argument("--assigned-reads", type=int, default=0,
                    help="reads the stage-1 classifier assigned to this taxon. Used for "
                         "mapped_fraction_of_assigned - if most of the reads a k-mer "
                         "classifier attributed to a taxon do not survive alignment, the "
                         "k-mer evidence did not hold up (Wright et al. 2026, "
                         "doi:10.1099/mgen.0.001739)")
    ap.add_argument("--mask-bed", type=Path, default=None,
                    help="BED file, or directory of *.mask.bed, from "
                         "scripts/mask_pack_refs.py. Masked bases are excluded from "
                         "the breadth/occupancy denominators and reads lying entirely "
                         "within masked intervals are reported separately instead of "
                         "counting as evidence")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    if not args.bam.exists():
        sys.exit(f"ERROR: {args.bam} not found")

    mask = load_mask(args.mask_bed)
    if mask:
        sys.stderr.write(f"[coverage] mask loaded for {len(mask)} contig(s)\n")

    lengths, by_taxid = read_header(args.bam)

    if args.taxid not in by_taxid:
        sys.exit(
            f"ERROR: target taxid {args.taxid} is not present in the competitive bundle "
            f"{args.bam.name}. Bundle contains: {sorted(by_taxid)[:10]}"
        )

    target_refs = set(by_taxid[args.taxid])
    target_len_raw = sum(lengths[r] for r in target_refs)
    masked_bp = sum(masked_len_within(mask.get(contig_key(r), []), lengths[r])
                    for r in target_refs) if mask else 0
    # Judge the genome on the portion we trust.
    target_len = max(1, target_len_raw - masked_bp)

    # Windows must not be larger than the genome, or occupancy is always 1/1.
    # Small viral genomes need a proportionally smaller window.
    window = min(args.window, max(100, target_len // 50)) if target_len else args.window

    per_taxid = collect_alignments(args.bam, args.mapq_min, mask)
    (covered, total_bases, depth_stats, occupied_windows, covered10, masked_positions,
     possible_windows) = depth_profile(args.bam, args.mapq_min, target_refs,
                                       window, mask, lengths)

    tgt = per_taxid.get(args.taxid, {"reads": 0, "nm": 0, "aln": 0, "bases": 0})
    unique_reads = tgt["reads"]

    # Best competing neighbour, for the attribution ratio.
    neighbours = {t: v["reads"] for t, v in per_taxid.items() if t != args.taxid}
    best_nb, best_nb_reads = ("", 0)
    if neighbours:
        best_nb, best_nb_reads = max(neighbours.items(), key=lambda kv: kv[1])

    # Is that neighbour CREDIBLE, or is it a cross-mapping pedestal? See covered_by_taxid.
    # Reported, never enforced here: coverage_stats.py measures, score_evidence.py decides.
    nb_uniformity = ""
    if best_nb:
        cov_by_tax = covered_by_taxid(args.bam, args.mapq_min)
        nb_bp, nb_cov, nb_aln = cov_by_tax.get(best_nb, (0, 0, 0.0))
        if nb_bp > 0 and nb_aln > 0:
            expected = min(1.0, nb_aln / nb_bp)
            if expected > 0:
                nb_uniformity = f"{(nb_cov / nb_bp) / expected:.4f}"

    breadth = (covered / total_bases) if total_bases else 0.0
    breadth10 = (covered10 / total_bases) if total_bases else 0.0
    mean_depth, _depth_sd = depth_stats
    # CV over the whole reference (zeros included) - an organism present at low
    # depth everywhere has a lower CV than one stacked on a single locus.
    depth_cv = (_depth_sd / mean_depth) if mean_depth > 0 else 0.0
    # OCCUPANCY IS A FRACTION, AND IT USED TO BE ABLE TO EXCEED 1.0.
    #
    # window_hits is keyed (contig, position // window) - PER CONTIG, over raw
    # coordinates. The denominator was ceil(target_len / window) where target_len is
    # the mask-SUBTRACTED sum over ALL contigs. Two independent inflations followed.
    # Multi-contig: each contig contributes ceil(len_i/window) windows, and their sum
    # exceeds ceil(sum(len_i)/window) - eight influenza segments of 1,400 bp gave 56
    # real windows against an expected 50, so a fully covered reference scored 1.12.
    # Masking: a 5 Mb genome with 1.5 Mb masked expected 3,500 windows while raw
    # coordinates span 5,000 indices, capping occupancy at 1.43. Both loosened the
    # `occupancy >= 0.15` gate in the false-positive direction, worst on segmented
    # viruses, fragmented drafts, and the heavily masked references D5 says to
    # distrust most.
    #
    # possible_windows counts, in the same pass and on the same key, the windows that
    # hold at least one unmasked base - so numerator and denominator are now the same
    # kind of thing and occupancy cannot exceed 1.0.
    expected_windows = max(1, possible_windows)
    occupancy = occupied_windows / expected_windows if expected_windows else 0.0
    mean_identity = (1.0 - tgt["nm"] / tgt["aln"]) if tgt["aln"] else 0.0
    # Reads per million host-depleted reads - comparable across samples and
    # against the negative control.
    rpm = (unique_reads / args.total_reads * 1e6) if args.total_reads else 0.0
    nb_ratio = (unique_reads / best_nb_reads) if best_nb_reads else float("inf")
    # Wright et al. criterion (ii). Capped at 1.0: a competitive bundle can recruit
    # reads the stage-1 classifier attributed elsewhere, which would otherwise give a
    # nonsensical fraction above one.
    mapped_frac = min(1.0, unique_reads / args.assigned_reads) if args.assigned_reads else -1.0

    starts = tgt.get("starts") or {}
    n_start_positions = len(starts)
    dup_rate = (1.0 - n_start_positions / unique_reads) if unique_reads else 0.0
    fwd_prop = (tgt.get("forward", 0) / unique_reads) if unique_reads else 0.0
    mean_aln_len = (tgt["aln"] / unique_reads) if unique_reads else 0.0
    mapped_any = tgt.get("mapped_any_mapq", 0)
    mapq_pass_frac = (unique_reads / mapped_any) if mapped_any else 0.0

    row = {
        "sample": args.sample,
        "role": args.role,
        "taxid": args.taxid,
        "taxon": args.taxon,
        "group": args.group,
        "platform": (args.platform or "").strip().lower(),
        "reference_bp": target_len,
        "reference_bp_raw": target_len_raw,
        "masked_bp": masked_bp,
        "masked_fraction_of_reference": (f"{masked_bp / target_len_raw:.6f}"
                                        if target_len_raw else "0"),
        "reads_in_masked_regions": tgt.get("masked_reads", 0),
        "reference_contigs": len(target_refs),
        "unique_reads": unique_reads,
        "reads_per_million": f"{rpm:.4f}",
        "aligned_bases": tgt["bases"],
        "breadth": f"{breadth:.6f}",
        "coverage_10x": f"{breadth10:.6f}",
        "mean_depth": f"{mean_depth:.4f}",
        "depth_cv": f"{depth_cv:.4f}",
        "window_bp": window,
        "windows_expected": expected_windows,
        "windows_occupied": occupied_windows,
        "window_occupancy": f"{occupancy:.6f}",
        "mean_identity": f"{mean_identity:.6f}",
        "mean_alignment_length": f"{mean_aln_len:.1f}",
        "mapped_reads_any_mapq": mapped_any,
        "mapq_pass_fraction": f"{mapq_pass_frac:.4f}",
        "forward_proportion": f"{fwd_prop:.4f}",
        "duplication_rate_proxy": f"{dup_rate:.4f}",
        "distinct_start_positions": n_start_positions,
        "best_neighbour_taxid": best_nb,
        "best_neighbour_reads": best_nb_reads,
        # Blank when it cannot be computed, NOT 0: an absent measurement must not read as
        # "maximally incredible" and silently disable the attribution gate.
        "best_neighbour_uniformity": nb_uniformity,
        "neighbour_ratio": ("inf" if nb_ratio == float("inf") else f"{nb_ratio:.4f}"),
        "screen_assigned_reads": args.assigned_reads,
        "mapped_fraction_of_assigned": (f"{mapped_frac:.4f}" if mapped_frac >= 0 else ""),
        "total_reads_screened": args.total_reads,
        "mapq_min": args.mapq_min,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        fh.write("\t".join(row) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")

    if masked_bp:
        sys.stderr.write(
            f"[coverage] {args.sample} taxid={args.taxid}: {masked_bp}/{target_len_raw} bp "
            f"({masked_bp / target_len_raw:.2%}) of the reference is masked; "
            f"{tgt.get('masked_reads', 0)} read(s) fell entirely within masked "
            "intervals and are NOT counted as evidence\n")
    sys.stderr.write(
        f"[coverage] {args.sample} taxid={args.taxid} ({args.taxon}): "
        f"{unique_reads} reads, breadth={breadth:.4f}, occupancy={occupancy:.4f}, "
        f"nb_ratio={row['neighbour_ratio']} (best neighbour {best_nb or 'none'}: "
        f"{best_nb_reads} reads), mapped_frac="
        f"{row['mapped_fraction_of_assigned'] or 'n/a'}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
