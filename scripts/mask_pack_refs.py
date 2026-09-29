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
import os
import re
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

# One regex handles homopolymers and short tandem repeats: a unit of 1-6 bases
# followed by at least three more copies of itself. Matches are filtered by total span,
# so a homopolymer needs 20 bases and a multi-base unit needs 24.
#
# This replaces a windowed distinct-k-mer scan written in pure Python which, over the
# 5.7 Gbp of reference sequence in this pack, would have run for hours - and which also
# copied every sequence in order to append a sentinel character. Regex runs at C speed
# and catches the cases that actually drive spurious mapping.
TANDEM_RE = re.compile(r"([ACGT]{1,6}?)\1{3,}")


def low_complexity_intervals(seq: str, min_span: int = 24,
                             min_homopolymer: int = 20) -> list[tuple[int, int, str]]:
    """Homopolymer runs and short tandem repeats.

    NOT equivalent to dustmasker or sdust - deliberately a simple stdlib approximation
    so this step adds no dependency, and documented as such rather than presented as
    equivalent. It catches long single-base runs and di- to hexa-nucleotide tracts,
    which are the substrate for spurious mapping, and will miss subtler low-complexity
    sequence.
    """
    out: list[tuple[int, int, str]] = []
    for m in TANDEM_RE.finditer(seq):
        span = m.end() - m.start()
        threshold = min_homopolymer if len(m.group(1)) == 1 else min_span
        if span >= threshold:
            out.append((m.start(), m.end(), "lowcomplexity"))
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
        die(f"masking subject {subject_fa} is missing or empty, so source '{reason}' "
            f"would contribute 0 bp to every BED and then be cached as complete. "
            f"Refusing: fix the subject, or drop the source with --sources.")
    cmd = [mm2, "-x", preset, "-t", str(threads), "--secondary=no",
           str(subject_fa), str(query_fa)]
    per_seq: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
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
    rc = proc.wait()
    err = ""
    if proc.stderr is not None:
        try:
            err = proc.stderr.read()
        except (OSError, ValueError):
            # Reading the diagnostic must never replace the real failure below.
            err = ""
        proc.stderr.close()
    # A FAILED ALIGNMENT PASS MUST NOT LOOK LIKE A CLEAN GENOME.
    #
    # This previously ignored the return code, with stderr sent to DEVNULL. When
    # minimap2 died - OOM on the 3.1 Gbp host index, killed, a missing shared
    # library - the function returned {} and the caller logged
    # "host: 0 bp across 0 contig(s)". Every BED was then written without host
    # intervals, every manifest row said host_bp=0, and because the skip decision is
    # "does the BED exist", no later run would redo it. The masking step D5 calls the
    # highest-value step in the pipeline would have been silently off, leaving an
    # artifact that reads as "checked, none found".
    if rc != 0:
        tail = "\n    ".join((err or "").strip().splitlines()[-5:]) or "(no stderr)"
        die(f"minimap2 failed (exit {rc}) masking against {subject_fa.name} for "
            f"source '{reason}'. Refusing to write BEDs that would record 0 bp for "
            f"this source and then be cached as complete.\n"
            f"  minimap2 stderr:\n    {tail}")
    return {k: merge_intervals(v) for k, v in per_seq.items()}


# UniVec_Core carries ~2,300 sequences; anything near zero is a truncated download.
# Deliberately loose - this is a structural sanity check, not a version pin.
MIN_UNIVEC_SEQS = 100


def bed_sources(path: Path) -> set[str]:
    """The sources recorded in an existing BED's header.

    The skip decision used to be "does the BED file exist", which cannot tell a BED
    built from all four sources from one built by `--sources lowcomplexity` during a
    smoke test. The later full run then logged SKIP for every genome and exited 0
    with host and crossgroup masking never done - and the manifest's host_bp column
    read 0, indistinguishable from "checked, none found".
    """
    try:
        with path.open() as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                if "sources=" in line:
                    tok = line.split("sources=", 1)[1].split()[0]
                    return {x for x in tok.split(",") if x}
    except OSError:
        return set()
    return set()


def fetch_univec(dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    # A TRUNCATED DOWNLOAD IS CACHED FOREVER.
    #
    # The skip was "exists and size > 0". shutil.copyfileobj over a dropped HTTPS
    # connection can end short without raising, and the truncated UniVec would then be
    # reused by every later run - quietly reducing vector masking, with a manifest
    # column reading a smaller vector_bp that looks like a cleaner genome. UniVec_Core
    # is ~1 MB and every record is a FASTA entry, so a cheap structural check catches
    # the realistic failure.
    if dest.exists() and dest.stat().st_size > 0:
        n_seqs = sum(1 for _ in read_fasta(dest))
        if n_seqs >= MIN_UNIVEC_SEQS:
            log(f"SKIP UniVec already present ({human(dest.stat().st_size)}, "
                f"{n_seqs} sequences)")
            return dest
        warn(f"cached {dest.name} holds only {n_seqs} sequence(s), below the "
             f"{MIN_UNIVEC_SEQS} expected - treating it as a truncated download and "
             f"re-fetching rather than masking against a partial vector set")
        dest.unlink(missing_ok=True)
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
    # stem -> owner stem, for genomes that ARE the same sequence under another taxid.
    alias_of: dict[str, str] = {}

    for r in rows:
        gpath = Path(r["path"])
        if not gpath.exists():
            warn(f"{r['taxid']} {r['name']}: {gpath} missing, skipping")
            continue
        stem = gpath.name.replace(".fna.gz", "").replace(".fna", "")
        bed_path = mask_dir / f"{stem}.mask.bed"
        if bed_path.exists() and not args.force:
            have = bed_sources(bed_path)
            missing = [x for x in sources if x not in have]
            if not missing:
                log(f"SKIP {r['taxid']} {r['name']} (bed present, sources complete)")
                continue
            log(f"REDO {r['taxid']} {r['name']}: existing bed was built from "
                f"{','.join(sorted(have)) or 'unknown'}, missing "
                f"{','.join(missing)}")
        seqs = dict(read_fasta(gpath))
        if not seqs:
            warn(f"{r['taxid']} {r['name']}: empty genome, skipping")
            continue
        # TWO TAXIDS, ONE ASSEMBLY.
        #
        # F. tularensis (263) / F. novicida (264) and S. dysenteriae (622) /
        # Shigella (620) each resolve to a single NCBI assembly, so their FASTAs
        # carry identical contig names. The old code warned, skipped the second
        # genome's contigs, and then computed genome_bp as a sum filtered to
        # contigs this stem owns - which is zero. The result was a manifest row
        # reading genome_bp=0 masked_bp=0 flag= and an intervals-free BED: byte
        # identical in meaning to "examined and found clean". Both shipped BEDs
        # were in that state. Since the sequence is the SAME sequence, the honest
        # answer is the owner's intervals, recorded as shared.
        shared = [c for c in seqs if c in contig_owner]
        if shared:
            owners = {contig_owner[c] for c in shared}
            same = (len(shared) == len(seqs) and len(owners) == 1)
            owner = next(iter(owners)) if len(owners) == 1 else ""
            if same and all(seqs[c] == genome_seqs[owner].get(c) for c in seqs):
                alias_of[stem] = owner
                genome_seqs[stem] = seqs
                stem_meta[stem] = r
                todo.append((stem, r))
                log(f"SHARED SEQUENCE {r['taxid']} {r['name']}: identical to "
                    f"{owner}; masking it once and writing both BEDs from that "
                    f"result (flagged SHARED_SEQUENCE in the manifest)")
                continue
            die(f"{r['taxid']} {r['name']} ({gpath.name}) shares contig name(s) "
                f"{', '.join(shared[:3])} with {owner or 'several genomes'} but the "
                f"sequence differs. Masking cannot be attributed, and attributing it "
                f"to the first genome would silently leave this one unmasked. "
                f"Rename the contigs or drop one reference.")
        genome_seqs[stem] = seqs
        stem_meta[stem] = r
        for cname, seq in seqs.items():
            contig_owner[cname] = stem
            contig_len[cname] = len(seq)
        todo.append((stem, r))

    if not todo:
        log("nothing to do (all BEDs present and built from the requested sources; "
            "use --force to redo)")
        if not (mask_dir / "mask_manifest.tsv").exists():
            warn("mask_manifest.tsv is missing and cannot be rebuilt without "
                 "re-masking. The BEDs are the operative artifact, but the manifest "
                 "is the auditable record D5 cites - re-run with --force to restore it.")
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
        # For a shared-sequence alias the intervals live under the owner's name.
        owner = alias_of.get(stem, stem)
        genome_bp = sum(len(v) for k, v in seqs.items() if contig_owner.get(k) == owner)
        merged = {}
        for cname in seqs:
            if contig_owner.get(cname) != owner:
                continue
            ivs = per_contig.get(cname)
            if ivs:
                merged[cname] = merge_intervals(ivs)

        masked_bp = sum(total_len(v) for v in merged.values())
        frac = masked_bp / genome_bp if genome_bp else 0.0

        # PER-SOURCE FOOTPRINTS, NOT A PARTITION OF masked_bp.
        #
        # The old accounting walked the MERGED intervals and charged every base of a
        # merged span in full to each contributing reason, so an interval tagged
        # "lowcomplexity+vector" added its whole length to both columns. A 0-100
        # lowcomplexity and a 50-150 vector hit merged to 150 bp and reported 150 +
        # 150 = 300. Recomputing from the raw per-source intervals gives each source
        # its true footprint. The columns still need not sum to masked_bp, because
        # sources genuinely overlap - that is a fact about the data, not an error.
        by_reason: dict[str, int] = defaultdict(int)
        for cname in seqs:
            if contig_owner.get(cname) != owner:
                continue
            raw = per_contig.get(cname) or []
            for src in ALL_SOURCES:
                own = [(a, b, src) for a, b, reason in raw
                       if src in reason.split("+")]
                if own:
                    by_reason[src] += total_len(merge_intervals(own))

        flag = ""
        if stem in alias_of:
            flag = f"SHARED_SEQUENCE:{alias_of[stem]}"
        if frac > args.max_masked_fraction:
            flag = f"OVER_MASKED,{flag}" if flag else "OVER_MASKED"
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
            # .part + rename, like every other writer in this file. Written straight to
            # the final path, an interrupt left a truncated .masked.fna.gz that looks
            # like a valid gzip stream to anything that only checks existence.
            hm = mask_dir / f"{stem}.masked.fna.gz"
            hm_tmp = hm.with_name(hm.name + ".part")
            with gzip.open(hm_tmp, "wt") as out:
                for cname, seq in seqs.items():
                    if contig_owner.get(cname) != owner:
                        continue
                    lst = list(seq)
                    for a, b, _ in merged.get(cname, []):
                        for j in range(a, min(b, len(lst))):
                            lst[j] = "N"
                    out.write(f">{cname}\n{''.join(lst)}\n")
            hm_tmp.rename(hm)

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

        # MERGE, DO NOT REWRITE.
        #
        # `summary` holds only the genomes masked in THIS invocation, and `todo`
        # excludes every genome whose BED already exists. So a run interrupted at
        # genome 300 of 363 and then restarted rewrote the manifest with the 63
        # remaining rows, destroying the record of the 296 already masked - with all
        # 363 BEDs still on disk and exit code 0. `--only-taxid 263` cut it to one
        # row. This is the same truncation already fixed in build_pack_refs.py and
        # build_marker_db.py; it was never fixed here.
        fresh = {str(row["taxid"]) for row in summary}
        carried: list[dict] = []
        if mpath.exists():
            try:
                with mpath.open() as fh:
                    for old in csv.DictReader(fh, delimiter="\t"):
                        if str(old.get("taxid", "")).strip() not in fresh:
                            carried.append(old)
            except OSError as exc:
                die(f"cannot read {mpath} to merge ({exc}). Refusing to continue: "
                    f"rewriting it from this run alone would drop "
                    f"{'the other genomes' if fresh else 'everything'} masked earlier.")
        if carried:
            log(f"carrying forward {len(carried)} manifest row(s) from earlier runs")

        cols = list(summary[0])
        allrows = summary + [{c: old.get(c, "") for c in cols} for old in carried]
        allrows.sort(key=lambda row: str(row.get("name", "")))
        with mpath.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for row in allrows:
                w.writerow(row)
        log(f"wrote {mpath} ({len(allrows)} row(s): {len(summary)} this run, "
            f"{len(carried)} carried)")

        tot_g = sum(int(row["genome_bp"] or 0) for row in allrows)
        tot_m = sum(int(row["masked_bp"] or 0) for row in allrows)
        log("")
        pct = f"{tot_m / tot_g:.3%} overall" if tot_g else "no bases"
        log(f"SUMMARY: {len(allrows)} genome(s), {human(tot_g)} total, "
            f"{tot_m} bp masked ({pct})")
        for src in ALL_SOURCES:
            if src in sources:
                v = sum(int(row.get(f"{src}_bp") or 0) for row in allrows)
                log(f"  {src:<14} {v} bp")
        log("  (per-source figures are FOOTPRINTS and may overlap, so they need not")
        log("   sum to the masked total; only masked_bp is a partition of the genome)")
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
