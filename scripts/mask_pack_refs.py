#!/usr/bin/env python3
"""
mask_pack_refs.py - identify contaminant and untrustworthy regions in a pack's
reference genomes.

WHY THIS IS THE HIGHEST-VALUE STEP IN THE PIPELINE

Two independent 2025-2026 reviews reached the same conclusion: curating and masking
the reference database buys more specificity than anything available downstream,
including assembly, and it costs almost no sensitivity. VirMask-style removal of
recurrent contaminants cut false detections by 40-94% with 100% true-positive
preservation (J Clin Microbiol 2026). Conterminator-style cross-taxon screening is
described as table stakes.

The concrete failure mode, from Wright et al. 2026 (Microb Genom,
doi:10.1099/mgen.0.001739): a Helicobacter pylori call in stool had ZERO reads
mapping to its reference once checked, and the cause was traced to contaminating
PHAGE SEQUENCE INSIDE THE REFERENCE GENOME. Unmasked references do not merely fail to
help - they manufacture false positives, and no downstream threshold can tell the
difference because the reads genuinely do align.

SOFT-MASK, NOT HARD-MASK - the central design decision

This script does NOT overwrite sequence with Ns by default. It emits BED intervals
with a reason per interval, and downstream `coverage_stats.py --mask-bed`:

  * excludes masked bases from the breadth and occupancy DENOMINATORS, so a genome
    is judged on the part of it we trust;
  * counts reads falling only in masked regions SEPARATELY, as
    reads_in_masked_regions;
  * lets a finding supported mainly by masked-region reads be flagged rather than
    silently counted.

Hard-masking would make the Wright case invisible: the reads would simply stop
aligning and the taxon would quietly drop out, with no record of why. Soft-masking
makes it auditable, and auditable matters more than tidy for a surveillance tool.
`--hard-mask` is available for pipelines that need it, but it is not the default.

MASKING SOURCES (each independently switchable via --sources)

  lowcomplexity  Homopolymer runs, short tandem repeats, and windows of low
                 linguistic complexity. The classic spurious-mapping substrate.
                 Implemented in stdlib; approximate compared with dustmasker, and
                 documented as such rather than presented as equivalent.
  vector         NCBI UniVec_Core. Cloning-vector and adapter sequence left in
                 assemblies is common and is pure artifact.
  host           The pack's own host genome. Human sequence contaminating microbial
                 assemblies is well documented, and a host-derived region in a
                 pathogen reference collects host reads that survived depletion.
  crossgroup     Regions of a genome that also occur near-identically in a
                 PHYLOGENETICALLY DISTANT reference in this same pack (a different
                 group: bacteria / virus / parasite / fungus). This is the
                 Conterminator insight. Same-group hits are deliberately NOT masked,
                 because shared sequence between close relatives is real biology and
                 is exactly what competitive mapping exists to resolve.

A GUARD, because this step can quietly destroy a reference

If masking would exceed --max-masked-fraction of a genome, the genome is reported
as OVER_MASKED and its BED is still written but marked, rather than being applied
silently. A reference that is 60% masked is not a curated reference, it is a
misconfiguration - most likely a genuine phage or mobile-element-rich genome being
masked against something it legitimately resembles.

Needs `minimap2` on PATH for the alignment-based sources. Stdlib Python otherwise.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

DEFAULT_BTF_DB = Path("/media/8tb/bioThreatFinder")
# VERIFY: UniVec_Core lives here; confirmed reachable 2026-09-29 (HTTP 206).
UNIVEC_URL = "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core"
ALL_SOURCES = ("lowcomplexity", "vector", "host", "crossgroup")


def log(msg: str) -> None:
    sys.stderr.write(f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] INFO {msg}\n")


def warn(msg: str) -> None:
    sys.stderr.write(f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] WARN {msg}\n")


def die(msg: str) -> None:
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def which(name: str) -> str | None:
    p = shutil.which(name)
    if p:
        return p
    for cand in (Path("/home/bioinfo/miniconda3/envs/btf-refs/bin") / name,
                 Path("/home/bioinfo/miniconda3/bin") / name):
        if cand.exists() and os.access(cand, os.X_OK):
            return str(cand)
    return None


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}PB"


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------

def open_maybe_gz(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path: Path):
    """Yield (header_first_token, sequence). Sequences are uppercased."""
    name, chunks = None, []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip().upper())
    if name is not None:
        yield name, "".join(chunks)


# ---------------------------------------------------------------------------
# Interval algebra
# ---------------------------------------------------------------------------

def merge_intervals(ivs: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    """Merge overlapping half-open [start, end) intervals, joining reasons.

    Reasons are combined rather than one winning, because an interval that is both
    low-complexity AND matches a vector is more suspicious than either alone, and
    throwing half of that away would lose information for no benefit.
    """
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for start, end, reason in ivs[1:]:
        if start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
            existing = set(out[-1][2].split("+"))
            existing.add(reason)
            out[-1][2] = "+".join(sorted(existing))
        else:
            out.append([start, end, reason])
    return [(a, b, c) for a, b, c in out]


def total_len(ivs: list[tuple[int, int, str]]) -> int:
    return sum(b - a for a, b, _ in ivs)


# ---------------------------------------------------------------------------
# Source 1: low complexity (stdlib)
# ---------------------------------------------------------------------------

def low_complexity_intervals(seq: str, window: int = 64, step: int = 32,
                             k: int = 4, min_distinct: int = 18,
                             homopolymer: int = 20) -> list[tuple[int, int, str]]:
    """Homopolymers, short tandem repeats and low-linguistic-complexity windows.

    NOT equivalent to dustmasker or sdust - this is a deliberately simple stdlib
    approximation so the pipeline keeps no extra dependency for this step. It catches
    the cases that matter for spurious mapping (long single-base runs, di/tri-nucleotide
    tracts, and windows with few distinct k-mers) and will miss subtler ones.
    """
    n = len(seq)
    out: list[tuple[int, int, str]] = []

    # Homopolymer runs.
    run_start, run_char = 0, ""
    for i, c in enumerate(seq + "\0"):
        if c != run_char:
            if run_char and run_char in "ACGT" and i - run_start >= homopolymer:
                out.append((run_start, i, "lowcomplexity"))
            run_start, run_char = i, c

    # Windowed distinct-k-mer count.
    if n >= window:
        for start in range(0, n - window + 1, step):
            w = seq[start:start + window]
            if w.count("N") > window // 2:
                continue
            distinct = len({w[i:i + k] for i in range(len(w) - k + 1)})
            if distinct < min_distinct:
                out.append((start, start + window, "lowcomplexity"))
    return merge_intervals(out)


# ---------------------------------------------------------------------------
# Sources 2-4: minimap2 based
# ---------------------------------------------------------------------------

def minimap2_mask_intervals(mm2: str, subject_fa: Path, query_fa: Path, reason: str,
                            threads: int, min_len: int, min_identity: float,
                            preset: str = "asm5") -> dict[str, list[tuple[int, int, str]]]:
    """Intervals ON THE QUERY that align to the subject.

    The query is the genome we intend to mask and the subject is the contaminant
    source, so PAF query coordinates are what we want. Getting this backwards would
    mask the contaminant database instead, which would look like it worked.
    """
    if not subject_fa.exists() or subject_fa.stat().st_size == 0:
        return {}
    cmd = [mm2, "-x", preset, "-t", str(threads), "--secondary=no",
           str(subject_fa), str(query_fa)]
    per_seq: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        try:
            qname, qstart, qend = f[0], int(f[2]), int(f[3])
            nmatch, alnlen = int(f[9]), int(f[10])
        except ValueError:
            continue
        if alnlen < min_len:
            continue
        if alnlen and (nmatch / alnlen) < min_identity:
            continue
        per_seq[qname].append((qstart, qend, reason))
    proc.stdout.close()
    proc.wait()
    return {k: merge_intervals(v) for k, v in per_seq.items()}


def fetch_univec(dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        log(f"SKIP UniVec already present ({human(dest.stat().st_size)})")
        return dest
    log(f"downloading UniVec_Core from {UNIVEC_URL}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(UNIVEC_URL, timeout=120) as r, tmp.open("wb") as out:
        shutil.copyfileobj(r, out)
    tmp.rename(dest)
    log(f"UniVec_Core: {human(dest.stat().st_size)}")
    return dest


# ---------------------------------------------------------------------------

def read_manifest(path: Path) -> list[dict]:
    if not path.exists():
        die(f"no refs manifest at {path}. Run scripts/build_pack_refs.py first.")
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    usable = [r for r in rows if (r.get("path") or "").strip()]
    if not usable:
        die(f"{path} lists no genomes with a path - nothing to mask")
    return usable


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pack", required=True)
    ap.add_argument("--btf-db", type=Path, default=DEFAULT_BTF_DB)
    ap.add_argument("--packs-dir", type=Path,
                    default=Path(__file__).resolve().parent.parent / "packs")
    ap.add_argument("--sources", default=",".join(ALL_SOURCES),
                    help=f"comma-separated subset of {','.join(ALL_SOURCES)}")
    ap.add_argument("--univec", type=Path, default=None,
                    help="UniVec_Core FASTA; downloaded to the pack cache if omitted")
    ap.add_argument("--host-fasta", type=Path, default=None,
                    help="host genome; read from the pack's host.fasta if omitted")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-hit-len", type=int, default=200,
                    help="ignore alignments shorter than this [200]")
    ap.add_argument("--min-identity", type=float, default=0.95,
                    help="ignore alignments below this identity [0.95]")
    ap.add_argument("--max-masked-fraction", type=float, default=0.20,
                    help="flag a genome as OVER_MASKED above this fraction [0.20]")
    ap.add_argument("--hard-mask", action="store_true",
                    help="ALSO write N-masked FASTAs. Not the default: hard masking "
                         "makes a contaminated reference silently lose its reads "
                         "instead of reporting why")
    ap.add_argument("--only-taxid", action="append", default=[])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    sources = [s.strip() for s in args.sources.split(",") if s.strip()]
    bad = [s for s in sources if s not in ALL_SOURCES]
    if bad:
        die(f"unknown source(s): {', '.join(bad)}. Valid: {', '.join(ALL_SOURCES)}")

    refs_dir = args.btf_db / "packs" / args.pack / "refs"
    mask_dir = refs_dir / "masked"
    manifest_path = refs_dir / "manifest.tsv"
    cache = args.btf_db / "cache"

    all_rows = read_manifest(manifest_path)
    rows = all_rows
    if args.only_taxid:
        want = set(args.only_taxid)
        rows = [r for r in all_rows if r["taxid"] in want]
        if not rows:
            die(f"no genomes matched --only-taxid {sorted(want)}")
        log(f"--only-taxid: masking {len(rows)} of {len(all_rows)} genome(s)")

    log(f"pack={args.pack} genomes={len(rows)} sources={','.join(sources)}")
    log(f"masked output -> {mask_dir}")

    mm2 = which("minimap2")
    needs_mm2 = any(s in sources for s in ("vector", "host", "crossgroup"))
    if needs_mm2 and not mm2:
        die("minimap2 not found but required for vector/host/crossgroup masking.\n"
            "  mamba install -n btf-refs -c bioconda minimap2\n"
            "  (or run with --sources lowcomplexity)")

    # ---- prepare contaminant subjects -------------------------------------
    univec = None
    if "vector" in sources:
        univec = args.univec or (cache / "UniVec_Core.fna")
        if not args.dry_run:
            univec = fetch_univec(Path(univec))

    host_fa = args.host_fasta
    if "host" in sources and host_fa is None:
        # Read the host path out of pack.yaml, resolving ${btf_db}.
        pack_yaml = args.packs_dir / args.pack / "pack.yaml"
        if pack_yaml.exists():
            for line in pack_yaml.read_text().splitlines():
                t = line.strip()
                if t.startswith("fasta:"):
                    v = t.split(":", 1)[1].strip().strip('"\'')
                    host_fa = Path(v.replace("${btf_db}", str(args.btf_db)))
                    break
        if host_fa is None or not Path(host_fa).exists():
            warn(f"host genome not found (looked for {host_fa}); "
                 "dropping the 'host' source. Run scripts/fetch_host_genome.sh to "
                 "enable it - host sequence inside a pathogen reference collects the "
                 "host reads that survived depletion.")
            sources = [s for s in sources if s != "host"]

    # Group genomes so crossgroup can compare only across distant groups.
    #
    # Built from ALL manifest rows, not the --only-taxid subset. Restricting the
    # comparison set to the subset would make a subset run silently UNDER-MASK -
    # it would report a clean genome because it had almost nothing to compare
    # against, which is indistinguishable from a genome that really is clean.
    by_group: dict[str, list[dict]] = defaultdict(list)
    for r in all_rows:
        by_group[(r.get("group") or "unknown").strip().lower()].append(r)
    log("groups: " + ", ".join(f"{g}={len(v)}" for g, v in sorted(by_group.items())))

    if args.dry_run:
        log("[dry-run] would mask:")
        for g, items in sorted(by_group.items()):
            log(f"[dry-run]   {g}: {len(items)} genome(s)")
        if "crossgroup" in sources:
            for g in sorted(by_group):
                others = sum(len(v) for k, v in by_group.items() if k != g)
                log(f"[dry-run]   crossgroup for {g}: vs {others} genome(s) from other groups")
            log("[dry-run]   (crossgroup subjects are built from the full manifest even "
                "under --only-taxid, so a subset run does not under-mask)")
        log("[dry-run] nothing written")
        return 0

    mask_dir.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(tempfile.mkdtemp(prefix="btf_mask_", dir=str(cache)))

    # -----------------------------------------------------------------------
    # BATCHED BY SOURCE, NOT BY GENOME.
    #
    # The obvious implementation maps each genome against each contaminant source in
    # turn, which makes minimap2 re-index the 932 MB human genome once per genome -
    # 276 times here. Concatenating all queries into one run per source turns that
    # into a handful of index builds. Contig names are assembly accessions and so are
    # unique across the pack, which is what makes the results attributable afterwards.
    # -----------------------------------------------------------------------
    contig_owner: dict[str, str] = {}          # contig -> genome stem
    contig_len: dict[str, int] = {}
    genome_seqs: dict[str, dict[str, str]] = {}
    stem_meta: dict[str, dict] = {}
    group_query: dict[str, Path] = {}
    todo = []

    for r in rows:
        gpath = Path(r["path"])
        if not gpath.exists():
            warn(f"{r['taxid']} {r['name']}: {gpath} missing, skipping")
            continue
        stem = gpath.name.replace(".fna.gz", "").replace(".fna", "")
        bed_path = mask_dir / f"{stem}.mask.bed"
        if bed_path.exists() and not args.force:
            log(f"SKIP {r['taxid']} {r['name']} (bed present)")
            continue
        seqs = dict(read_fasta(gpath))
        if not seqs:
            warn(f"{r['taxid']} {r['name']}: empty genome, skipping")
            continue
        genome_seqs[stem] = seqs
        stem_meta[stem] = r
        for cname, seq in seqs.items():
            if cname in contig_owner:
                warn(f"duplicate contig name {cname} in {stem} and "
                     f"{contig_owner[cname]}; masking will be attributed to the first")
                continue
            contig_owner[cname] = stem
            contig_len[cname] = len(seq)
        todo.append((stem, r))

    if not todo:
        log("nothing to do (all BEDs present; use --force to redo)")
        shutil.rmtree(tmpdir, ignore_errors=True)
        return 0

    log(f"masking {len(todo)} genome(s), {sum(contig_len.values())} bp total")

    # One concatenated query for the whole pack, plus one per group for crossgroup.
    all_query = tmpdir / "all_queries.fna"
    per_group_handles: dict[str, object] = {}
    with all_query.open("w") as out_all:
        for stem, r in todo:
            g = (r.get("group") or "unknown").strip().lower()
            if g not in per_group_handles:
                gp = tmpdir / f"query_{g}.fna"
                group_query[g] = gp
                per_group_handles[g] = gp.open("w")
            for cname, seq in genome_seqs[stem].items():
                if contig_owner.get(cname) != stem:
                    continue
                rec = f">{cname}\n{seq}\n"
                out_all.write(rec)
                per_group_handles[g].write(rec)
    for fh in per_group_handles.values():
        fh.close()
    log(f"query set: {human(all_query.stat().st_size)}")

    # ---- accumulate intervals per contig ----
    per_contig: dict[str, list[tuple[int, int, str]]] = defaultdict(list)

    if "lowcomplexity" in sources:
        log("scanning low complexity (stdlib)...")
        for stem, _ in todo:
            for cname, seq in genome_seqs[stem].items():
                if contig_owner.get(cname) != stem:
                    continue
                per_contig[cname].extend(low_complexity_intervals(seq))
        log(f"  low complexity: {sum(total_len(v) for v in per_contig.values())} bp")

    if "vector" in sources and univec:
        log("aligning against UniVec_Core (one pass)...")
        hits = minimap2_mask_intervals(mm2, Path(univec), all_query, "vector",
                                      args.threads, args.min_hit_len, args.min_identity)
        for cname, ivs in hits.items():
            per_contig[cname].extend(ivs)
        log(f"  vector: {sum(total_len(v) for v in hits.values())} bp "
            f"across {len(hits)} contig(s)")

    if "host" in sources and host_fa:
        log(f"aligning against host {Path(host_fa).name} (one pass, this is the slow one)...")
        hits = minimap2_mask_intervals(mm2, Path(host_fa), all_query, "host",
                                      args.threads, args.min_hit_len, args.min_identity)
        for cname, ivs in hits.items():
            per_contig[cname].extend(ivs)
        log(f"  host: {sum(total_len(v) for v in hits.values())} bp "
            f"across {len(hits)} contig(s)")

    if "crossgroup" in sources:
        for g, qpath in sorted(group_query.items()):
            others = [r for k, v in by_group.items() if k != g for r in v]
            if not others:
                continue
            sub = tmpdir / f"notgroup_{g}.fna"
            n = 0
            with sub.open("w") as out:
                for r in others:
                    pth = Path(r["path"])
                    if not pth.exists():
                        continue
                    for cname, seq in read_fasta(pth):
                        out.write(f">{r['taxid']}|{cname}\n{seq}\n")
                        n += 1
            if not n:
                continue
            log(f"crossgroup: {g} ({human(qpath.stat().st_size)}) vs {n} sequences "
                f"from other groups ({human(sub.stat().st_size)})")
            hits = minimap2_mask_intervals(mm2, sub, qpath, "crossgroup", args.threads,
                                          args.min_hit_len, args.min_identity)
            for cname, ivs in hits.items():
                per_contig[cname].extend(ivs)
            log(f"  crossgroup {g}: {sum(total_len(v) for v in hits.values())} bp "
                f"across {len(hits)} contig(s)")
            sub.unlink(missing_ok=True)

    # ---- attribute back to genomes and write ----
    summary = []
    over_masked = []
    for stem, r in todo:
        taxid, name = r["taxid"], r["name"]
        group = (r.get("group") or "unknown").strip().lower()
        seqs = genome_seqs[stem]
        genome_bp = sum(len(v) for k, v in seqs.items() if contig_owner.get(k) == stem)
        merged = {}
        for cname in seqs:
            if contig_owner.get(cname) != stem:
                continue
            ivs = per_contig.get(cname)
            if ivs:
                merged[cname] = merge_intervals(ivs)

        masked_bp = sum(total_len(v) for v in merged.values())
        frac = masked_bp / genome_bp if genome_bp else 0.0
        by_reason: dict[str, int] = defaultdict(int)
        for ivs in merged.values():
            for a, b, reason in ivs:
                for part in reason.split("+"):
                    by_reason[part] += b - a

        flag = ""
        if frac > args.max_masked_fraction:
            flag = "OVER_MASKED"
            over_masked.append((taxid, name, frac))
            warn(f"{taxid} {name}: {frac:.1%} of the genome would be masked, above "
                 f"--max-masked-fraction {args.max_masked_fraction:.0%}. BED is written "
                 "and MARKED, not applied silently - a reference this heavily masked is "
                 "a misconfiguration, most likely a genuinely mobile-element-rich "
                 "genome masked against something it legitimately resembles.")

        bed_path = mask_dir / f"{stem}.mask.bed"
        tmp_bed = bed_path.with_suffix(".bed.part")
        with tmp_bed.open("w") as out:
            out.write(f"# taxid={taxid} name={name} group={group}\n")
            out.write(f"# genome_bp={genome_bp} masked_bp={masked_bp} "
                      f"masked_fraction={frac:.6f} flag={flag}\n")
            out.write(f"# sources={','.join(sources)} min_hit_len={args.min_hit_len} "
                      f"min_identity={args.min_identity}\n")
            for cname in sorted(merged):
                for a, b, reason in merged[cname]:
                    out.write(f"{cname}\t{a}\t{b}\t{reason}\n")
        tmp_bed.rename(bed_path)

        if args.hard_mask:
            hm = mask_dir / f"{stem}.masked.fna.gz"
            with gzip.open(hm, "wt") as out:
                for cname, seq in seqs.items():
                    if contig_owner.get(cname) != stem:
                        continue
                    lst = list(seq)
                    for a, b, _ in merged.get(cname, []):
                        for j in range(a, min(b, len(lst))):
                            lst[j] = "N"
                    out.write(f">{cname}\n{''.join(lst)}\n")

        summary.append({
            "taxid": taxid, "name": name, "group": group, "role": r.get("role", ""),
            "assembly_level": r.get("assembly_level", ""),
            "genome_bp": genome_bp, "masked_bp": masked_bp,
            "masked_fraction": f"{frac:.6f}",
            "lowcomplexity_bp": by_reason.get("lowcomplexity", 0),
            "vector_bp": by_reason.get("vector", 0),
            "host_bp": by_reason.get("host", 0),
            "crossgroup_bp": by_reason.get("crossgroup", 0),
            "n_masked_intervals": sum(len(v) for v in merged.values()),
            "bed_path": str(bed_path), "flag": flag,
        })

    shutil.rmtree(tmpdir, ignore_errors=True)

    if summary:
        mpath = mask_dir / "mask_manifest.tsv"
        with mpath.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary[0]), delimiter="\t")
            w.writeheader()
            for s in summary:
                w.writerow(s)
        log(f"wrote {mpath}")

        tot_g = sum(s["genome_bp"] for s in summary)
        tot_m = sum(s["masked_bp"] for s in summary)
        log("")
        log(f"SUMMARY: {len(summary)} genome(s), {human(tot_g)} total, "
            f"{tot_m} bp masked ({tot_m / tot_g:.3%} overall)")
        for src in ALL_SOURCES:
            if src in sources:
                v = sum(s[f"{src}_bp"] for s in summary)
                log(f"  {src:<14} {v} bp")
        if over_masked:
            log("")
            warn(f"{len(over_masked)} genome(s) flagged OVER_MASKED - review before use:")
            for taxid, name, frac in sorted(over_masked, key=lambda x: -x[2])[:20]:
                warn(f"  {taxid} {name}: {frac:.1%}")
        log("")
        log("Masking is SOFT: coverage_stats.py --mask-bed excludes these intervals "
            "from the breadth denominator and counts reads landing in them "
            "separately. Nothing is deleted.")
    else:
        warn("no genomes were masked")
    return 0


if __name__ == "__main__":
    sys.exit(main())
