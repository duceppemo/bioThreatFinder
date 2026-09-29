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
import statistics
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


def collect_alignments(bam: Path, mapq_min: int):
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
    per_taxid = defaultdict(lambda: {"reads": 0, "nm": 0, "aln": 0, "bases": 0,
                                     "mapped_any_mapq": 0, "forward": 0,
                                     "starts": None})
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
        rec0["starts"][(f[2], f[3], flag & 0x10)] += 1
        if not (flag & 0x10):
            rec0["forward"] += 1
        ref_span, aln_len = cigar_spans(f[5])
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


def depth_profile(bam: Path, mapq_min: int, target_refs: set[str], window: int):
    """Per-base depth for the target's references only.

    ``samtools depth -a`` over the whole competitive bundle would emit a row per
    base of every neighbour genome too - tens of millions of lines for no
    benefit. Only the target needs a depth profile.
    """
    covered = 0
    covered10 = 0
    total = 0
    depths: list[int] = []
    window_hits: dict[tuple[str, int], bool] = {}

    cmd = ["samtools", "depth", "-a", "-Q", str(mapq_min), str(bam)]
    for line in stream(cmd):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3 or f[0] not in target_refs:
            continue
        d = int(f[2])
        total += 1
        depths.append(d)
        if d > 0:
            covered += 1
            window_hits[(f[0], int(f[1]) // window)] = True
        if d >= 10:
            covered10 += 1

    return covered, total, depths, len(window_hits), covered10


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--taxid", required=True, help="target taxid (must appear in the bundle)")
    ap.add_argument("--taxon", default="", help="target taxon name, for the report")
    ap.add_argument("--group", default="", help="bacteria|virus|parasite|fungus")
    ap.add_argument("--role", default="sample", help="sample|control")
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
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    if not args.bam.exists():
        sys.exit(f"ERROR: {args.bam} not found")

    lengths, by_taxid = read_header(args.bam)

    if args.taxid not in by_taxid:
        sys.exit(
            f"ERROR: target taxid {args.taxid} is not present in the competitive bundle "
            f"{args.bam.name}. Bundle contains: {sorted(by_taxid)[:10]}"
        )

    target_refs = set(by_taxid[args.taxid])
    target_len = sum(lengths[r] for r in target_refs)

    # Windows must not be larger than the genome, or occupancy is always 1/1.
    # Small viral genomes need a proportionally smaller window.
    window = min(args.window, max(100, target_len // 50)) if target_len else args.window

    per_taxid = collect_alignments(args.bam, args.mapq_min)
    covered, total_bases, depths, occupied_windows, covered10 = depth_profile(
        args.bam, args.mapq_min, target_refs, window
    )

    tgt = per_taxid.get(args.taxid, {"reads": 0, "nm": 0, "aln": 0, "bases": 0})
    unique_reads = tgt["reads"]

    # Best competing neighbour, for the attribution ratio.
    neighbours = {t: v["reads"] for t, v in per_taxid.items() if t != args.taxid}
    best_nb, best_nb_reads = ("", 0)
    if neighbours:
        best_nb, best_nb_reads = max(neighbours.items(), key=lambda kv: kv[1])

    breadth = (covered / total_bases) if total_bases else 0.0
    breadth10 = (covered10 / total_bases) if total_bases else 0.0
    mean_depth = statistics.fmean(depths) if depths else 0.0
    # CV over the whole reference (zeros included) - an organism present at low
    # depth everywhere has a lower CV than one stacked on a single locus.
    depth_cv = (statistics.pstdev(depths) / mean_depth) if depths and mean_depth > 0 else 0.0
    expected_windows = max(1, -(-target_len // window))  # ceil
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
        "reference_bp": target_len,
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
