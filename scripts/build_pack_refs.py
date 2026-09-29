#!/usr/bin/env python3
"""
build_pack_refs.py - fetch reference genomes for a bioThreatFinder domain pack
                     and assemble the competitive-mapping references used by
                     stage-2 confirmation.

Why this exists
---------------
Stage 1 (Kraken2 / sylph) says "something that looks like Bacillus anthracis is
here". That is not a finding. Stage 2 re-aligns the candidate reads against the
target genome AND against its nearest relatives simultaneously, then assigns
each read to its best hit. Without those relatives in the index, every read from
a common environmental cousin (B. cereus, F. novicida, Y. pseudotuberculosis,
any Brucella) piles onto the high-consequence target and produces a false
alarm. The near_neighbour_taxids column of the resolved watchlist is what makes
that discrimination possible, so this script treats neighbours as first-class
download targets, not as an optional extra.

What it writes (paths must match packs/<pack>/pack.yaml)
-------------------------------------------------------
  <db_root>/bioThreatFinder/packs/<pack>/refs/<taxid>_<name>.fna.gz
  <db_root>/bioThreatFinder/packs/<pack>/refs/manifest.tsv
  <db_root>/bioThreatFinder/packs/<pack>/refs/competitive/<taxid>.fna.gz
  <db_root>/bioThreatFinder/packs/<pack>/refs/competitive/<taxid>.<preset>.mmi

Every sequence in a competitive FASTA is renamed to

    >taxid|<taxid>|<original header>

so that a downstream aligner can attribute a read to a source organism from the
reference name alone, with no extra lookup table to fall out of sync.

Dependencies
------------
Python 3.8+ standard library only, plus these external command-line tools:
  * `datasets`  - NCBI datasets CLI v2 (preferred download route)
  * `curl`      - required for the documented FTP fallback
  * `minimap2`  - only needed for the competitive indexes (--skip-index to omit)

Usage
-----
  build_pack_refs.py --pack human
  build_pack_refs.py --pack human --dry-run
  build_pack_refs.py --pack human --threads 32 --max-per-taxon 2
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

LOG = logging.getLogger("build_pack_refs")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_DB_ROOT = Path("/media/36tb/db")
DEFAULT_TAXDUMP = Path("/media/36tb/db/taxdump")
CONDA_BIN = Path("/home/bioinfo/miniconda3/bin")

# The shared volume runs at ~97% capacity. Refuse to start rather than fill it.
# Bacteria and viruses are trivial; the eukaryotic parasites (Schistosoma,
# Echinococcus, Trichinella, Acanthamoeba) and their competitive sets plus the
# minimap2 indexes are what make this number as large as it is.
DEFAULT_MIN_FREE_GIB = 150

MANIFEST_COLUMNS = [
    "taxid", "name", "group", "role", "accession", "assembly_level",
    "genome_size", "n_contigs", "path", "md5", "downloaded_utc",
]

# Ranked best-first. Used both for --assembly-level parsing and for choosing
# among candidate assemblies.
ASSEMBLY_LEVELS = ["complete", "chromosome", "scaffold", "contig"]
ASSEMBLY_LEVEL_ALIASES = {
    "complete genome": "complete",
    "complete": "complete",
    "chromosome": "chromosome",
    "scaffold": "scaffold",
    "contig": "contig",
}

# RefSeq FTP division summaries. Far smaller than the monolithic
# assembly_summary_refseq.txt (which is ~0.5 GB), and the pack's `group` column
# tells us which ones can possibly contain a given taxon.
# VERIFY: these paths are long-standing but re-check if a fetch 404s.
NCBI_FTP = "https://ftp.ncbi.nlm.nih.gov"
DIVISIONS_BY_GROUP = {
    "bacteria": ["bacteria", "archaea"],
    "virus": ["viral"],
    "fungus": ["fungi"],
    # Helminths live under "invertebrate", single-celled parasites under
    # "protozoa"; the pack's `parasite` group spans both.
    "parasite": ["protozoa", "invertebrate"],
}
ALL_DIVISIONS = ["bacteria", "archaea", "viral", "fungi", "protozoa", "invertebrate"]

# assembly_summary.txt column indices (0-based), format version 1.2.
AS_COL_ACCESSION = 0
AS_COL_REFSEQ_CATEGORY = 4
AS_COL_TAXID = 5
AS_COL_SPECIES_TAXID = 6
AS_COL_ORGANISM = 7
AS_COL_VERSION_STATUS = 10
AS_COL_ASSEMBLY_LEVEL = 11
AS_COL_GENOME_REP = 13
AS_COL_FTP_PATH = 19

CACHE_MAX_AGE_DAYS = 7


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sanitise(name: str) -> str:
    """Filesystem-safe stem. 'Bacillus anthracis' -> 'Bacillus_anthracis'."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return re.sub(r"_+", "_", s).strip("_") or "unknown"


def md5_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def which(name: str) -> Optional[str]:
    """PATH first, then the shared conda prefix (not always on PATH under
    Nextflow/cron), then the minimap2 that ships with the ONT dorado server."""
    found = shutil.which(name)
    if found:
        return found
    candidate = CONDA_BIN / name
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    if name == "minimap2":
        dorado = Path("/home/bioinfo/prog/ont-dorado-server/bin/minimap2")
        if dorado.is_file() and os.access(dorado, os.X_OK):
            LOG.warning("using minimap2 bundled with ont-dorado-server (%s); "
                        "it may predate the validated release", dorado)
            return str(dorado)
    return None


def normalise_level(raw: str) -> str:
    return ASSEMBLY_LEVEL_ALIASES.get((raw or "").strip().lower(), (raw or "").strip().lower())


def dig(obj, *keys, default=None):
    """Tolerant nested lookup - the datasets CLI has reshuffled its JSON between
    minor versions, so never index blindly into it."""
    cur = obj
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur if cur is not None else default


def check_free_space(path: Path, need_gib: int) -> None:
    probe = path
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    usage = shutil.disk_usage(probe)
    free_gib = usage.free / (1024 ** 3)
    LOG.info("free space at %s: %.0f GiB (need %d GiB)", probe, free_gib, need_gib)
    if free_gib < need_gib:
        raise SystemExit(
            f"ERROR: insufficient free space: {free_gib:.0f} GiB available at {probe}, "
            f"need {need_gib} GiB.\n"
            f"       Free space, lower --min-free-gib if you know better, or point "
            f"--db-root at a larger volume.\n"
            f"       Nothing was downloaded."
        )


# ---------------------------------------------------------------------------
# Taxonomy names (from the pinned local taxdump)
# ---------------------------------------------------------------------------

class Taxdump:
    """Lazy reader for the unpacked NCBI taxdump already on the shared volume.

    Needed for two things:
      * naming near-neighbour taxids that are not themselves watchlist rows
        (the watchlist only carries names for its own targets);
      * expanding a taxid to its descendants for the FTP fallback, which can
        only match assembly_summary rows on an exact taxid.
    """

    def __init__(self, path: Path):
        self.path = path
        self._names: Optional[Dict[int, str]] = None
        self._children: Optional[Dict[int, List[int]]] = None

    @property
    def available(self) -> bool:
        return (self.path / "names.dmp").is_file() and (self.path / "nodes.dmp").is_file()

    def _load_names(self) -> Dict[int, str]:
        if self._names is None:
            self._names = {}
            if not (self.path / "names.dmp").is_file():
                return self._names
            LOG.info("loading scientific names from %s", self.path / "names.dmp")
            with (self.path / "names.dmp").open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    parts = [p.strip() for p in line.rstrip("\n").rstrip("\t|").split("\t|")]
                    if len(parts) >= 4 and parts[3] == "scientific name":
                        try:
                            self._names[int(parts[0])] = parts[1]
                        except ValueError:
                            continue
        return self._names

    def name(self, taxid: int) -> Optional[str]:
        return self._load_names().get(taxid)

    def _load_children(self) -> Dict[int, List[int]]:
        if self._children is None:
            self._children = {}
            if not (self.path / "nodes.dmp").is_file():
                return self._children
            LOG.info("loading taxonomy tree from %s (needed for the FTP fallback)",
                     self.path / "nodes.dmp")
            with (self.path / "nodes.dmp").open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    parts = line.split("\t|\t", 2)
                    if len(parts) < 2:
                        continue
                    try:
                        tid, parent = int(parts[0]), int(parts[1])
                    except ValueError:
                        continue
                    if tid != parent:
                        self._children.setdefault(parent, []).append(tid)
        return self._children

    def descendants(self, taxid: int) -> set:
        """Inclusive subtree of taxid."""
        children = self._load_children()
        out, stack = {taxid}, [taxid]
        while stack:
            for child in children.get(stack.pop(), ()):
                if child not in out:
                    out.add(child)
                    stack.append(child)
        return out


# ---------------------------------------------------------------------------
# Watchlist
# ---------------------------------------------------------------------------

REQUIRED_WATCHLIST_COLUMNS = [
    "taxid", "name", "group", "near_neighbour_taxids",
]


def read_resolved_watchlist(path: Path) -> List[dict]:
    if not path.is_file():
        raise SystemExit(
            f"ERROR: resolved watchlist not found: {path}\n"
            f"       Run scripts/resolve_watchlist.py first - this script deliberately "
            f"does not resolve names itself."
        )
    with path.open(encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    missing = [c for c in REQUIRED_WATCHLIST_COLUMNS if c not in (reader.fieldnames or [])]
    if missing:
        raise SystemExit(f"ERROR: {path} is missing required columns: {', '.join(missing)}")
    rows = []
    for row in reader:
        taxid = (row.get("taxid") or "").strip()
        if not taxid:
            continue
        rows.append(row)
    if not rows:
        raise SystemExit(f"ERROR: {path} has no data rows")
    return rows


def parse_neighbours(raw: str) -> List[int]:
    """near_neighbour_taxids is ';'-separated in the resolved watchlist."""
    out = []
    for tok in re.split(r"[;,\s]+", (raw or "").strip()):
        if tok.isdigit():
            tid = int(tok)
            if tid not in out:
                out.append(tid)
    return out


# ---------------------------------------------------------------------------
# Assembly candidate selection
# ---------------------------------------------------------------------------

class Candidate:
    __slots__ = ("accession", "organism", "taxid", "level", "refseq_category",
                 "ftp_path", "size_hint", "contig_hint", "source")

    def __init__(self, accession: str, organism: str = "", taxid: int = 0,
                 level: str = "", refseq_category: str = "", ftp_path: str = "",
                 size_hint: int = 0, contig_hint: int = 0, source: str = ""):
        self.accession = accession
        self.organism = organism
        self.taxid = taxid
        self.level = normalise_level(level)
        self.refseq_category = (refseq_category or "").strip().lower()
        self.ftp_path = ftp_path
        self.size_hint = size_hint
        self.contig_hint = contig_hint
        self.source = source

    def sort_key(self) -> tuple:
        cat_rank = {"reference genome": 0, "reference": 0,
                    "representative genome": 1, "representative": 1}.get(self.refseq_category, 2)
        level_rank = ASSEMBLY_LEVELS.index(self.level) if self.level in ASSEMBLY_LEVELS else len(ASSEMBLY_LEVELS)
        # GCF (RefSeq, curated) before GCA (GenBank, submitter-supplied).
        source_rank = 0 if self.accession.startswith("GCF_") else 1
        # Among equals, the larger assembly is the more complete one; negate so
        # that plain ascending sort puts it first.
        return (cat_rank, level_rank, source_rank, -self.size_hint, self.accession)

    def __repr__(self) -> str:
        return f"<Candidate {self.accession} {self.level or '?'} {self.refseq_category or 'na'} via {self.source}>"


def datasets_candidates(datasets_bin: str, taxid: int, levels: Sequence[str],
                        want: int, timeout: int = 600) -> List[Candidate]:
    """Ask the NCBI datasets CLI for assemblies in the taxid's subtree.

    Two passes: reference-only first (one curated assembly per species, which is
    what we want for bacteria), then unrestricted if that came back empty or we
    need more than one (viruses).
    """
    results: Dict[str, Candidate] = {}

    def run(extra: List[str]) -> None:
        cmd = [datasets_bin, "summary", "genome", "taxon", str(taxid),
               "--as-json-lines"] + extra
        LOG.debug("running: %s", " ".join(cmd))
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            LOG.warning("taxid %s: datasets summary timed out after %ss", taxid, timeout)
            return
        except OSError as exc:
            LOG.warning("taxid %s: could not run datasets: %s", taxid, exc)
            return
        if proc.returncode != 0:
            # An empty result set is a non-zero exit in some versions, so this is
            # informational rather than fatal - the FTP fallback still gets a go.
            msg = (proc.stderr or "").strip().splitlines()
            LOG.debug("taxid %s: datasets summary exit %s: %s",
                      taxid, proc.returncode, msg[-1] if msg else "")
            return
        # Cap the number of lines parsed: a genus-level neighbour can return
        # tens of thousands of assemblies and we only ever keep a handful.
        for i, line in enumerate(proc.stdout.splitlines()):
            if i >= 5000:
                LOG.debug("taxid %s: stopped parsing datasets output at 5000 records", taxid)
                break
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            acc = rec.get("accession") or dig(rec, "assembly_info", "assembly_accession", default="")
            if not acc:
                continue
            cand = Candidate(
                accession=acc,
                organism=dig(rec, "organism", "organism_name", default="") or "",
                taxid=int(dig(rec, "organism", "tax_id", default=0) or 0),
                level=dig(rec, "assembly_info", "assembly_level", default="") or "",
                refseq_category=dig(rec, "assembly_info", "refseq_category", default="") or "",
                size_hint=int(dig(rec, "assembly_stats", "total_sequence_length", default=0) or 0),
                contig_hint=int(dig(rec, "assembly_stats", "number_of_contigs", default=0) or 0),
                source="datasets",
            )
            results.setdefault(acc, cand)

    run(["--reference"])
    if len(results) < want:
        run([])

    keep = [c for c in results.values() if not levels or c.level in levels or not c.level]
    keep.sort(key=Candidate.sort_key)
    return keep


def fetch_assembly_summary(division: str, cache_dir: Path, curl_bin: str,
                           dry_run: bool) -> Optional[Path]:
    """Download (and cache) one RefSeq division's assembly_summary.txt."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / f"assembly_summary_{division}.txt"
    if dest.is_file() and dest.stat().st_size > 0:
        age_days = (time.time() - dest.stat().st_mtime) / 86400.0
        if age_days < CACHE_MAX_AGE_DAYS:
            return dest
        LOG.info("cached %s is %.1f days old; refreshing", dest.name, age_days)
    if dry_run:
        LOG.info("[dry-run] would fetch %s/genomes/refseq/%s/assembly_summary.txt",
                 NCBI_FTP, division)
        return dest if dest.is_file() else None
    url = f"{NCBI_FTP}/genomes/refseq/{division}/assembly_summary.txt"
    LOG.info("fetching %s", url)
    tmp = dest.with_suffix(".txt.part")
    cmd = [curl_bin, "--fail", "--location", "--silent", "--show-error",
           "--retry", "5", "--retry-delay", "10", "--connect-timeout", "30",
           "--output", str(tmp), url]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not tmp.exists() or tmp.stat().st_size == 0:
        LOG.warning("could not fetch assembly summary for division %s: %s",
                    division, (proc.stderr or "").strip())
        tmp.unlink(missing_ok=True)
        return dest if dest.is_file() else None
    tmp.replace(dest)
    return dest


def ftp_candidates(taxid: int, group: str, levels: Sequence[str], cache_dir: Path,
                   curl_bin: str, taxdump: Taxdump, dry_run: bool) -> List[Candidate]:
    """Documented fallback for when the `datasets` CLI is unavailable or fails.

    Route: RefSeq assembly_summary.txt per division -> filter on taxid /
    species_taxid -> take ftp_path and append the conventional
    '<basename>_genomic.fna.gz'. That naming convention is part of NCBI's
    published FTP layout, not a guess, but see VERIFY notes near NCBI_FTP.
    """
    divisions = DIVISIONS_BY_GROUP.get((group or "").strip().lower(), ALL_DIVISIONS)
    wanted = {taxid}
    out: List[Candidate] = []
    for division in divisions:
        summary = fetch_assembly_summary(division, cache_dir, curl_bin, dry_run)
        if summary is None or not summary.is_file():
            continue
        hits = _scan_assembly_summary(summary, wanted, levels)
        if not hits and taxdump.available:
            # Species-level queries usually match directly; strain/genus-level or
            # recently reorganised virus taxids need the subtree.
            wanted = taxdump.descendants(taxid)
            hits = _scan_assembly_summary(summary, wanted, levels)
        out.extend(hits)
    out.sort(key=Candidate.sort_key)
    return out


def _scan_assembly_summary(summary: Path, wanted: set, levels: Sequence[str]) -> List[Candidate]:
    hits: List[Candidate] = []
    with summary.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= AS_COL_FTP_PATH:
                continue
            try:
                tid = int(f[AS_COL_TAXID])
                sp_tid = int(f[AS_COL_SPECIES_TAXID])
            except ValueError:
                continue
            if tid not in wanted and sp_tid not in wanted:
                continue
            if f[AS_COL_VERSION_STATUS].strip().lower() not in ("latest", ""):
                continue
            if f[AS_COL_GENOME_REP].strip().lower() == "partial":
                continue
            level = normalise_level(f[AS_COL_ASSEMBLY_LEVEL])
            if levels and level and level not in levels:
                continue
            ftp_path = f[AS_COL_FTP_PATH].strip()
            if not ftp_path or ftp_path == "na":
                continue
            hits.append(Candidate(
                accession=f[AS_COL_ACCESSION].strip(),
                organism=f[AS_COL_ORGANISM].strip(),
                taxid=tid,
                level=level,
                refseq_category=f[AS_COL_REFSEQ_CATEGORY].strip(),
                # assembly_summary still advertises ftp://; https is the same tree
                # and behaves better through proxies.
                ftp_path=ftp_path.replace("ftp://", "https://", 1),
                source="ftp",
            ))
    return hits


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def download_via_datasets(datasets_bin: str, accession: str, dest: Path,
                          threads: int, timeout: int = 3600) -> bool:
    """`datasets download genome accession` yields a zip; unpack and recompress.

    Note: the resulting .fna.gz will NOT match NCBI's published md5, because we
    re-gzip the plain FASTA from the archive. The manifest therefore records the
    md5 of the file we wrote, which is what re-verification on the next run
    needs anyway.
    """
    import zipfile

    tmp_out = dest.with_suffix(dest.suffix + ".part")
    with tempfile.TemporaryDirectory(dir=str(dest.parent), prefix=".datasets.") as td:
        zip_path = Path(td) / "g.zip"
        cmd = [datasets_bin, "download", "genome", "accession", accession,
               "--include", "genome", "--no-progressbar",
               "--filename", str(zip_path)]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except (subprocess.TimeoutExpired, OSError) as exc:
            LOG.warning("%s: datasets download failed: %s", accession, exc)
            return False
        if proc.returncode != 0 or not zip_path.is_file():
            LOG.warning("%s: datasets download failed: %s", accession,
                        (proc.stderr or "").strip().splitlines()[-1:])
            return False
        try:
            with zipfile.ZipFile(zip_path) as zf:
                members = [m for m in zf.namelist()
                           if m.endswith((".fna", ".fna.gz", ".fa", ".fa.gz"))]
                if not members:
                    LOG.warning("%s: no FASTA inside datasets archive", accession)
                    return False
                # An accession can carry several FASTA members (e.g. unplaced
                # scaffolds); concatenate them all into one .fna.gz.
                with gzip.open(tmp_out, "wb", compresslevel=6) as out:
                    for member in sorted(members):
                        with zf.open(member) as src:
                            if member.endswith(".gz"):
                                with gzip.open(src, "rb") as gz:
                                    shutil.copyfileobj(gz, out, length=1 << 20)
                            else:
                                shutil.copyfileobj(src, out, length=1 << 20)
        except Exception as exc:  # zipfile / gzip corruption, truncated download
            LOG.warning("%s: could not unpack datasets archive: %s", accession, exc)
            tmp_out.unlink(missing_ok=True)
            return False
    tmp_out.replace(dest)
    return True


def download_via_ftp(curl_bin: str, cand: Candidate, dest: Path) -> bool:
    """Fetch <ftp_path>/<basename>_genomic.fna.gz, verifying against the
    md5checksums.txt that NCBI publishes beside it."""
    base = cand.ftp_path.rstrip("/").rsplit("/", 1)[-1]
    fname = f"{base}_genomic.fna.gz"
    url = f"{cand.ftp_path.rstrip('/')}/{fname}"
    tmp = dest.with_suffix(dest.suffix + ".part")
    cmd = [curl_bin, "--fail", "--location", "--silent", "--show-error",
           "--retry", "5", "--retry-delay", "10", "--retry-connrefused",
           "--connect-timeout", "30", "--output", str(tmp), url]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not tmp.exists() or tmp.stat().st_size == 0:
        LOG.warning("%s: FTP download failed (%s): %s", cand.accession, url,
                    (proc.stderr or "").strip())
        tmp.unlink(missing_ok=True)
        return False

    expected = _ncbi_md5(curl_bin, cand.ftp_path, fname)
    if expected:
        got = md5_of(tmp)
        if got != expected:
            LOG.warning("%s: md5 mismatch (expected %s, got %s) - discarding",
                        cand.accession, expected, got)
            tmp.unlink(missing_ok=True)
            return False
        LOG.debug("%s: md5 verified against NCBI md5checksums.txt", cand.accession)
    else:
        LOG.debug("%s: no md5checksums.txt entry; skipping remote verification",
                  cand.accession)
    tmp.replace(dest)
    return True


def _ncbi_md5(curl_bin: str, ftp_dir: str, fname: str) -> Optional[str]:
    url = f"{ftp_dir.rstrip('/')}/md5checksums.txt"
    proc = subprocess.run(
        [curl_bin, "--fail", "--location", "--silent", "--retry", "3",
         "--connect-timeout", "20", url],
        capture_output=True, text=True)
    if proc.returncode != 0:
        return None
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1].lstrip("./") == fname:
            return parts[0]
    return None


# ---------------------------------------------------------------------------
# FASTA statistics and competitive concatenation
# ---------------------------------------------------------------------------

def fasta_stats(path: Path) -> Tuple[int, int]:
    """(genome_size, n_contigs) from a gzipped FASTA, streamed."""
    total = 0
    contigs = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                contigs += 1
            else:
                total += len(line.strip())
    return total, contigs


def append_with_taxid_prefix(src: Path, taxid: int, out_handle) -> None:
    """Copy a gzipped FASTA into an open gzip handle, rewriting every header to
    'taxid|<taxid>|<original header>'.

    The prefix is what lets stage 2 attribute a read to an organism straight
    from the reference name in the BAM, so it must be present on every sequence,
    including plasmids and segments - a segmented virus with an unprefixed
    segment would be silently unattributable.
    """
    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                out_handle.write(f">taxid|{taxid}|{line[1:].rstrip()}\n")
            else:
                out_handle.write(line)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(
        prog="build_pack_refs.py",
        description="Download reference genomes for a bioThreatFinder pack and build "
                    "competitive-mapping references for stage-2 confirmation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  build_pack_refs.py --pack human --dry-run
  build_pack_refs.py --pack human --threads 32
  build_pack_refs.py --pack human --max-per-taxon 2 --assembly-level complete,chromosome,scaffold

notes:
  * Idempotent and resumable. A genome already present with a matching md5 is
    skipped; interrupt and rerun freely. The manifest is rewritten after every
    taxon so progress survives a crash.
  * Taxa with no downloadable genome are WARNINGs, are listed in the final
    summary, and are written to the manifest with an empty path. They are never
    silently dropped - a missing near neighbour means a target is mapped
    non-competitively, which is exactly the false-positive mode this pipeline
    exists to avoid.
""")
    p.add_argument("--pack", required=True,
                   help="pack name, e.g. human (must exist under --packs-dir)")
    p.add_argument("--packs-dir", type=Path, default=repo_root / "packs",
                   help="directory containing pack definitions (default: %(default)s)")
    p.add_argument("--btf-db", type=Path, default=None,
                   help="explicit base for project-owned artifacts, e.g. "
                        "/media/8tb/bioThreatFinder. Preferred over --db-root, and "
                        "checked against refs_dir in the pack's pack.yaml - a "
                        "mismatch makes the pipeline report no_reference for every "
                        "taxon, which reads as 'checked and clear'.")
    p.add_argument("--db-root", type=Path, default=DEFAULT_DB_ROOT,
                   help="shared database root; retains the implicit "
                        "/bioThreatFinder suffix (default: %(default)s)")
    p.add_argument("--taxdump", type=Path, default=DEFAULT_TAXDUMP,
                   help="unpacked NCBI taxdump, used to name neighbour taxids and to "
                        "expand subtrees for the FTP fallback (default: %(default)s)")
    p.add_argument("--threads", type=int, default=16,
                   help="threads for minimap2 indexing (default: %(default)s)")
    p.add_argument("--max-per-taxon", type=int, default=1,
                   help="assemblies per taxon (default: 1; see --max-per-taxon-virus)")
    p.add_argument("--max-per-taxon-virus", type=int, default=3,
                   help="assemblies per taxon for group=virus (default: %(default)s). "
                        "Viruses need more: intra-species diversity is high enough that "
                        "one RefSeq genome can lose reads from a divergent lineage "
                        "(e.g. Influenza A, Lassa, the hantaviruses).")
    p.add_argument("--assembly-level", default="complete,chromosome",
                   help="comma-separated, best-first, from: complete,chromosome,scaffold,contig "
                        "(default: %(default)s). Relaxed automatically per taxon if nothing "
                        "matches, because many parasite and fungal references never get past "
                        "scaffold level.")
    p.add_argument("--no-relax", action="store_true",
                   help="do not relax --assembly-level when a taxon has no match")
    p.add_argument("--skip-index", action="store_true",
                   help="build competitive FASTAs but not minimap2 indexes")
    p.add_argument("--mm2-preset", default="sr,map-ont",
                   help="minimap2 presets to index per competitive FASTA "
                        "(default: %(default)s). Two indexes because k-mer and window "
                        "size are baked in at index time: an -x sr index is not usable "
                        "for ONT reads and vice versa.")
    p.add_argument("--min-free-gib", type=int, default=DEFAULT_MIN_FREE_GIB,
                   help="abort if less than this much free space (default: %(default)s)")
    p.add_argument("--only-taxid", action="append", default=[],
                   help="restrict to these target taxids (repeatable); for debugging")
    p.add_argument("--force", action="store_true",
                   help="re-download and rebuild even when outputs exist")
    p.add_argument("--dry-run", action="store_true",
                   help="print the plan and exit without downloading anything")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    return p.parse_args(argv)


def load_manifest(path: Path) -> Dict[Tuple[str, str], dict]:
    """Existing manifest keyed by (taxid, accession) - the resume state."""
    out: Dict[Tuple[str, str], dict] = {}
    if not path.is_file():
        return out
    with path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row.get("taxid", ""), row.get("accession", ""))] = row
    LOG.info("resuming from existing manifest: %d row(s) in %s", len(out), path)
    return out


def write_manifest(path: Path, rows: Iterable[dict]) -> None:
    tmp = path.with_suffix(".tsv.part")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    tmp.replace(path)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%SZ",
        stream=sys.stderr,
    )
    logging.Formatter.converter = time.gmtime

    levels = [normalise_level(x) for x in args.assembly_level.split(",") if x.strip()]
    bad = [l for l in levels if l not in ASSEMBLY_LEVELS]
    if bad:
        raise SystemExit(f"ERROR: unknown --assembly-level value(s): {', '.join(bad)}")
    presets = [p.strip() for p in args.mm2_preset.split(",") if p.strip()]

    pack_dir = args.packs_dir / args.pack
    watchlist_path = pack_dir / "watchlist.resolved.tsv"
    # Project artifacts now live on NVMe under an explicit --btf-db base. --db-root
    # is retained for back-compat and keeps its implicit /bioThreatFinder suffix.
    # Getting this wrong is not a cosmetic error: refs built under the wrong base
    # make the pipeline report no_reference for every watchlist taxon, which reads
    # as "checked and clear" rather than "never built".
    btf_base = args.btf_db if args.btf_db else (args.db_root / "bioThreatFinder")
    refs_dir = btf_base / "packs" / args.pack / "refs"
    comp_dir = refs_dir / "competitive"
    manifest_path = refs_dir / "manifest.tsv"
    cache_dir = btf_base / "cache"

    LOG.info("pack=%s btf_db=%s threads=%d", args.pack, btf_base, args.threads)

    # Cross-check where the PIPELINE will look, before hours of downloading.
    #
    # An earlier version of this check substituted the base we were given into
    # pack.yaml's ${btf_db} template and compared - which is tautological and always
    # passes. The invariant that actually matters is different: this script must
    # write where the pipeline READS, and the pipeline's base is params.btf_db in
    # nextflow.config. Getting it wrong makes every watchlist taxon come back
    # no_reference, which an analyst reads as "checked and clear" rather than
    # "never built".
    pack_yaml = args.packs_dir / args.pack / "pack.yaml"
    nf_config = args.packs_dir.parent / "nextflow.config"

    refs_template = None
    if pack_yaml.exists():
        for line in pack_yaml.read_text().splitlines():
            t = line.strip()
            if t.startswith("refs_dir:"):
                refs_template = t.split(":", 1)[1].strip().strip('"\'')
                break

    pipeline_base = None
    if nf_config.exists():
        for line in nf_config.read_text().splitlines():
            t = line.strip()
            if t.startswith("btf_db") and "=" in t:
                val = t.split("=", 1)[1].split("//")[0].strip().strip('"\'')
                if val and not val.startswith("$") and val != "null":
                    pipeline_base = Path(val)
                break

    if refs_template is None:
        LOG.warning("no refs_dir in %s; cannot verify the pipeline will look here",
                    pack_yaml)
    elif "${btf_db}" in refs_template:
        if pipeline_base is None:
            LOG.warning("could not read params.btf_db from %s; cannot verify that the "
                        "pipeline will look in %s", nf_config, refs_dir)
        elif pipeline_base != btf_base:
            LOG.error("BASE MISMATCH - refusing to download")
            LOG.error("  this script would write to : %s", refs_dir)
            LOG.error("  the pipeline reads from    : %s",
                      pipeline_base / "packs" / args.pack / "refs")
            LOG.error("pack.yaml resolves refs_dir from params.btf_db in "
                      "nextflow.config, so genomes written here would never be found. "
                      "Every watchlist taxon would report no_reference, which reads "
                      "as 'checked and clear' rather than 'never built'.")
            LOG.error("Pass --btf-db %s, or change params.btf_db in nextflow.config.",
                      pipeline_base)
            return 2
        else:
            LOG.info("base agrees with params.btf_db in nextflow.config: %s", btf_base)
    else:
        if Path(refs_template) != refs_dir:
            LOG.error("BASE MISMATCH - refusing to download")
            LOG.error("  this script would write to : %s", refs_dir)
            LOG.error("  pack.yaml refs_dir is      : %s", refs_template)
            return 2
        LOG.info("refs_dir matches the absolute path in pack.yaml: %s", refs_dir)
    LOG.info("refs dir:  %s", refs_dir)
    LOG.info("manifest:  %s", manifest_path)

    rows = read_resolved_watchlist(watchlist_path)
    if args.only_taxid:
        keep = set(args.only_taxid)
        rows = [r for r in rows if r["taxid"].strip() in keep]
        if not rows:
            raise SystemExit("ERROR: --only-taxid matched no watchlist rows")

    taxdump = Taxdump(args.taxdump)
    if not taxdump.available:
        LOG.warning("taxdump not usable at %s; neighbour taxids will be named from "
                    "the downloaded assembly metadata instead, and the FTP fallback "
                    "cannot expand subtrees", args.taxdump)

    # ---- build the work plan -------------------------------------------
    # A taxid can be both a target and someone else's near neighbour (e.g. 235
    # Brucella abortus). Download it once; role=target wins in the manifest.
    known_names = {r["taxid"].strip(): (r.get("name") or "").strip() for r in rows}
    known_groups = {r["taxid"].strip(): (r.get("group") or "").strip() for r in rows}

    tasks: Dict[int, dict] = {}       # taxid -> {name, group, role, want}
    targets: List[dict] = []          # ordered target records

    for r in rows:
        tid = int(r["taxid"].strip())
        group = (r.get("group") or "").strip().lower()
        want = args.max_per_taxon_virus if group == "virus" else args.max_per_taxon
        tasks[tid] = {"taxid": tid, "name": (r.get("name") or "").strip(),
                      "group": group, "role": "target", "want": want}
        neighbours = parse_neighbours(r.get("near_neighbour_taxids", ""))
        targets.append({"taxid": tid, "name": (r.get("name") or "").strip(),
                        "group": group, "neighbours": neighbours})

    for t in targets:
        for ntid in t["neighbours"]:
            if ntid in tasks:
                continue
            key = str(ntid)
            name = known_names.get(key) or (taxdump.name(ntid) if taxdump.available else None) \
                or f"taxid_{ntid}"
            # A neighbour inherits its target's group: it is by construction a
            # close relative, so the group is the same, and the group drives both
            # --max-per-taxon and which RefSeq divisions the fallback searches.
            tasks[ntid] = {"taxid": ntid, "name": name,
                           "group": known_groups.get(key) or t["group"],
                           "role": "near_neighbour",
                           "want": args.max_per_taxon_virus if t["group"] == "virus"
                                   else args.max_per_taxon}

    n_targets = len(targets)
    n_neighbours = sum(1 for v in tasks.values() if v["role"] == "near_neighbour")
    LOG.info("plan: %d target taxa, %d additional near-neighbour taxa, %d unique downloads",
             n_targets, n_neighbours, len(tasks))
    LOG.info("plan: %d competitive FASTA(s) -> %s", n_targets, comp_dir)
    LOG.info("plan: assembly levels %s%s; %d genome(s) per taxon (%d for viruses)",
             ",".join(levels), "" if args.no_relax else " (relaxable)",
             args.max_per_taxon, args.max_per_taxon_virus)

    # ---- tools ----------------------------------------------------------
    datasets_bin = which("datasets")
    curl_bin = which("curl")
    minimap2_bin = which("minimap2")
    if datasets_bin:
        LOG.info("using datasets CLI: %s", datasets_bin)
    else:
        LOG.warning("NCBI datasets CLI not found - falling back to the curl/FTP "
                    "assembly_summary route. Install the preferred tool with: "
                    "%s/mamba install -c conda-forge ncbi-datasets-cli", CONDA_BIN)
    if not curl_bin:
        raise SystemExit("ERROR: curl not found and it is required for the FTP fallback")
    if not minimap2_bin and not args.skip_index and not args.dry_run:
        raise SystemExit(
            "ERROR: minimap2 not found; competitive indexes cannot be built.\n"
            f"       Install it ({CONDA_BIN}/mamba install -c bioconda minimap2) "
            "or pass --skip-index.")

    if args.dry_run:
        LOG.info("[dry-run] would create %s and %s", refs_dir, comp_dir)
        for t in targets:
            comp = comp_dir / f"{t['taxid']}.fna.gz"
            LOG.info("[dry-run] %s (%s, %s) + %d neighbour(s) -> %s%s",
                     t["taxid"], t["name"], t["group"], len(t["neighbours"]), comp,
                     "" if args.skip_index else " + " + ", ".join(
                         f"{t['taxid']}.{p}.mmi" for p in presets))
        LOG.info("[dry-run] manifest would be written to %s", manifest_path)
        LOG.info("[dry-run] nothing downloaded")
        return 0

    check_free_space(btf_base.parent if btf_base.parent.exists() else args.db_root,
                     args.min_free_gib)
    refs_dir.mkdir(parents=True, exist_ok=True)
    comp_dir.mkdir(parents=True, exist_ok=True)

    prior = load_manifest(manifest_path)
    manifest_rows: List[dict] = []
    # taxid -> list of local FASTA paths successfully present
    genomes_by_taxid: Dict[int, List[Path]] = {}
    failed_taxa: List[dict] = []

    # ---- download loop --------------------------------------------------
    ordered = sorted(tasks.values(), key=lambda t: (t["role"] != "target", t["taxid"]))
    for i, task in enumerate(ordered, start=1):
        tid, name, group, role, want = (task["taxid"], task["name"], task["group"],
                                        task["role"], task["want"])
        stem = f"{tid}_{sanitise(name)}"
        LOG.info("[%d/%d] taxid %s (%s, %s, %s) - want %d assembly/assemblies",
                 i, len(ordered), tid, name, group, role, want)

        # Reuse anything the previous run already fetched for this taxid.
        got: List[Path] = []
        if not args.force:
            for (mtid, macc), row in prior.items():
                if mtid != str(tid) or not row.get("path"):
                    continue
                p = Path(row["path"])
                if not p.is_file() or p.stat().st_size == 0:
                    continue
                if row.get("md5") and md5_of(p) != row["md5"]:
                    LOG.warning("taxid %s: md5 changed for %s - will re-download", tid, p)
                    continue
                LOG.info("  skip (already present, md5 OK): %s", p.name)
                got.append(p)
                manifest_rows.append({**row, "role": "target" if role == "target" else row.get("role", role)})
                if len(got) >= want:
                    break

        if len(got) < want:
            use_levels = list(levels)
            cands: List[Candidate] = []
            for attempt in range(2):
                if datasets_bin:
                    cands = datasets_candidates(datasets_bin, tid, use_levels, want)
                if not cands:
                    cands = ftp_candidates(tid, group, use_levels, cache_dir, curl_bin,
                                           taxdump, args.dry_run)
                if cands or args.no_relax or attempt == 1:
                    break
                # Relax rather than report a spurious zero: many eukaryotic
                # parasite references are scaffold- or contig-level only.
                use_levels = list(ASSEMBLY_LEVELS)
                LOG.info("  no assembly at level(s) %s; relaxing to any level",
                         ",".join(levels))

            already = {Path(p).name for p in got}
            for cand in cands:
                if len(got) >= want:
                    break
                # The primary genome keeps the plain <taxid>_<name>.fna.gz name
                # promised by pack.yaml; extras are disambiguated by accession.
                if not got:
                    dest = refs_dir / f"{stem}.fna.gz"
                else:
                    dest = refs_dir / f"{stem}.{cand.accession}.fna.gz"
                if dest.name in already:
                    continue
                if dest.is_file() and dest.stat().st_size > 0 and not args.force:
                    LOG.info("  skip (file exists): %s", dest.name)
                else:
                    LOG.info("  downloading %s (%s, %s) -> %s",
                             cand.accession, cand.level or "level?", cand.source, dest.name)
                    ok = False
                    if datasets_bin:
                        ok = download_via_datasets(datasets_bin, cand.accession, dest, args.threads)
                    if not ok and cand.ftp_path:
                        ok = download_via_ftp(curl_bin, cand, dest)
                    if not ok and not cand.ftp_path:
                        # Candidate came from `datasets` so we have no FTP path;
                        # look it up in the division summaries before giving up.
                        for alt in ftp_candidates(tid, group, [], cache_dir, curl_bin,
                                                  taxdump, args.dry_run):
                            if alt.accession == cand.accession:
                                ok = download_via_ftp(curl_bin, alt, dest)
                                break
                    if not ok:
                        LOG.warning("  taxid %s: could not download %s", tid, cand.accession)
                        continue
                try:
                    size, contigs = fasta_stats(dest)
                except OSError as exc:
                    LOG.warning("  taxid %s: %s is unreadable (%s); discarding", tid, dest, exc)
                    dest.unlink(missing_ok=True)
                    continue
                if contigs == 0:
                    LOG.warning("  taxid %s: %s contains no sequences; discarding", tid, dest.name)
                    dest.unlink(missing_ok=True)
                    continue
                got.append(dest)
                already.add(dest.name)
                manifest_rows.append({
                    "taxid": str(tid), "name": name, "group": group, "role": role,
                    "accession": cand.accession,
                    "assembly_level": cand.level or "unknown",
                    "genome_size": str(size), "n_contigs": str(contigs),
                    "path": str(dest), "md5": md5_of(dest),
                    "downloaded_utc": utc_now(),
                })

        if got:
            genomes_by_taxid[tid] = got
        else:
            # Never a silent skip: record the hole so it shows up in the manifest
            # and in the final summary.
            LOG.warning("taxid %s (%s, %s, %s): NO downloadable genome", tid, name, group, role)
            failed_taxa.append(task)
            manifest_rows.append({
                "taxid": str(tid), "name": name, "group": group, "role": role,
                "accession": "", "assembly_level": "", "genome_size": "0",
                "n_contigs": "0", "path": "", "md5": "", "downloaded_utc": utc_now(),
            })

        # Rewrite after each taxon so an interrupted run resumes cleanly.
        write_manifest(manifest_path, manifest_rows)

    # ---- competitive FASTAs + indexes -----------------------------------
    fully_built, partially_built, unbuildable = [], [], []
    duplicate_neighbours: list[tuple[str, str, list[str]]] = []

    # Path -> accession. The dedup below MUST key on the assembly accession, not on
    # the file path and not on the file's md5:
    #   * paths differ, because each taxid gets its own <taxid>_<name>.fna.gz even
    #     when `datasets` returned the same assembly for both;
    #   * gzip md5s differ too, because the container embeds the original filename
    #     and mtime, so two byte-identical genomes hash differently.
    # Only the accession identifies the underlying assembly. An earlier version of
    # this fix keyed on resolved path, silently changed nothing, and was caught only
    # by re-reading the written bundle - which is why the post-write check below now
    # exists.
    accession_by_path: dict[str, str] = {}
    for mrow in manifest_rows:
        mp = (mrow.get("path") or "").strip()
        acc = (mrow.get("accession") or "").strip()
        if mp and acc:
            accession_by_path[str(Path(mp).resolve())] = acc

    for t in targets:
        tid = t["taxid"]
        comp = comp_dir / f"{tid}.fna.gz"
        target_genomes = genomes_by_taxid.get(tid, [])
        if not target_genomes:
            LOG.warning("target %s (%s): no target genome, cannot build a competitive "
                        "reference - stage 2 cannot confirm this taxon", tid, t["name"])
            unbuildable.append(t)
            continue

        present_neighbours = [n for n in t["neighbours"] if genomes_by_taxid.get(n)]
        missing_neighbours = [n for n in t["neighbours"] if not genomes_by_taxid.get(n)]

        if comp.is_file() and comp.stat().st_size > 0 and not args.force:
            LOG.info("target %s: competitive FASTA exists, skipping build: %s", tid, comp.name)
        else:
            LOG.info("target %s (%s): building competitive FASTA from 1 target + "
                     "%d/%d neighbour(s)", tid, t["name"],
                     len(present_neighbours), len(t["neighbours"]))
            tmp = comp.with_suffix(".gz.part")
            try:
                # DE-DUPLICATE BY SOURCE FILE.
                #
                # A neighbour taxid can resolve to the SAME assembly as its target -
                # `datasets` returns the parent species' reference for a subspecies
                # taxid, and a genus taxid returns one of its species. Observed here
                # for F. tularensis (263) vs subsp. novicida (264), which both
                # resolved to GCF_000833355.1.
                #
                # Writing that assembly twice under two taxid prefixes does not just
                # waste space, it INVERTS the result: every read maps equally well to
                # both copies, minimap2 assigns MAPQ 0, the MAPQ filter removes them,
                # unique_reads collapses to ~0 and neighbour_ratio approaches 1. The
                # taxon then reads as attributed_to_neighbour and can never confirm -
                # a systematic FALSE NEGATIVE, and worst of all on exactly the
                # target/near-neighbour pairs the bundle exists to discriminate.
                def _key(src) -> str:
                    rp = str(Path(src).resolve())
                    return accession_by_path.get(rp) or rp

                written: set[str] = set()
                dropped: list[tuple[str, str]] = []
                with gzip.open(tmp, "wt", compresslevel=6, encoding="utf-8") as out:
                    for src in target_genomes:
                        append_with_taxid_prefix(src, tid, out)
                        written.add(_key(src))
                    for ntid in present_neighbours:
                        for src in genomes_by_taxid[ntid]:
                            k = _key(src)
                            if k in written:
                                dropped.append((ntid, k))
                                continue
                            append_with_taxid_prefix(src, ntid, out)
                            written.add(k)
                tmp.replace(comp)

                # Verify what was actually written, rather than trusting the logic
                # above. Hash each sequence and assert no identical sequence appears
                # under two taxid prefixes.
                seen_md5: dict[str, str] = {}
                collisions: set[tuple[str, str]] = set()
                cur_tax, buf = None, []

                def _flush(cur_tax, buf):
                    if not cur_tax or not buf:
                        return
                    h = hashlib.md5("".join(buf).encode()).hexdigest()
                    prev = seen_md5.get(h)
                    if prev is not None and prev != cur_tax:
                        collisions.add(tuple(sorted((prev, cur_tax))))
                    else:
                        seen_md5[h] = cur_tax

                with gzip.open(comp, "rt") as chk:
                    for line in chk:
                        if line.startswith(">"):
                            _flush(cur_tax, buf)
                            bits = line[1:].split("|", 2)
                            cur_tax = bits[1] if len(bits) >= 2 else "?"
                            buf = []
                        else:
                            buf.append(line.strip())
                _flush(cur_tax, buf)

                if collisions:
                    LOG.error(
                        "target %s (%s): POST-WRITE CHECK FAILED - the bundle still "
                        "contains identical sequence under two taxid prefixes %s. Every "
                        "read over that sequence will get MAPQ 0 and this taxon cannot "
                        "confirm. Do not use this bundle.",
                        tid, t["name"], sorted(collisions))
                    unbuildable.append(t)
                    continue
                for ntid, fname in dropped:
                    LOG.warning(
                        "target %s (%s): neighbour taxid %s resolved to the SAME "
                        "assembly as the target (%s) and was DROPPED from the "
                        "competitive bundle. Keeping it would give every read MAPQ 0 "
                        "and make this taxon impossible to confirm. Consequence: this "
                        "bundle cannot discriminate %s from %s - curate a distinct "
                        "assembly for the neighbour to restore that discrimination.",
                        tid, t["name"], ntid, fname, tid, ntid)
                if dropped:
                    duplicate_neighbours.append((tid, t["name"],
                                                 [d[0] for d in dropped]))
            except OSError as exc:
                tmp.unlink(missing_ok=True)
                LOG.error("target %s: failed to write %s: %s", tid, comp, exc)
                unbuildable.append(t)
                continue

        if not args.skip_index:
            for preset in presets:
                mmi = comp_dir / f"{tid}.{preset}.mmi"
                if mmi.is_file() and mmi.stat().st_size > 0 and not args.force:
                    LOG.info("target %s: index exists, skipping: %s", tid, mmi.name)
                    continue
                LOG.info("target %s: minimap2 -x %s index -> %s", tid, preset, mmi.name)
                tmp_mmi = mmi.with_suffix(".mmi.part")
                proc = subprocess.run(
                    [minimap2_bin, "-x", preset, "-t", str(args.threads),
                     "-d", str(tmp_mmi), str(comp)],
                    capture_output=True, text=True)
                if proc.returncode != 0 or not tmp_mmi.exists() or tmp_mmi.stat().st_size == 0:
                    tmp_mmi.unlink(missing_ok=True)
                    LOG.error("target %s: minimap2 indexing failed: %s", tid,
                              (proc.stderr or "").strip().splitlines()[-3:])
                    continue
                tmp_mmi.replace(mmi)

        if missing_neighbours:
            partially_built.append((t, missing_neighbours))
        else:
            fully_built.append(t)

    write_manifest(manifest_path, manifest_rows)

    # ---- summary --------------------------------------------------------
    print("", file=sys.stderr)
    LOG.info("=========================== SUMMARY ===========================")
    LOG.info("pack:                 %s", args.pack)
    LOG.info("refs dir:             %s", refs_dir)
    LOG.info("manifest:             %s (%d rows)", manifest_path, len(manifest_rows))
    LOG.info("unique taxa attempted: %d (%d targets, %d neighbours)",
             len(tasks), n_targets, n_neighbours)
    LOG.info("genomes on disk:      %d", sum(len(v) for v in genomes_by_taxid.values()))
    LOG.info("targets fully built:      %d / %d", len(fully_built), n_targets)
    LOG.info("targets partially built:  %d / %d (target present, some neighbours missing)",
             len(partially_built), n_targets)
    LOG.info("targets NOT built:        %d / %d", len(unbuildable), n_targets)

    if partially_built:
        LOG.warning("--- partially built targets (competitive mapping is weakened; "
                    "a missing neighbour is a false-positive risk) ---")
        for t, missing in partially_built:
            LOG.warning("  %s %s: missing neighbour taxid(s) %s",
                        t["taxid"], t["name"], ",".join(str(m) for m in missing))
    if unbuildable:
        LOG.warning("--- targets with NO genome (stage 2 cannot confirm these) ---")
        for t in unbuildable:
            LOG.warning("  %s %s (%s)", t["taxid"], t["name"], t["group"])
    if failed_taxa:
        LOG.warning("--- taxa with zero downloadable genomes (recorded in the "
                    "manifest with an empty path) ---")
        for task in failed_taxa:
            LOG.warning("  %s %s (%s, %s)", task["taxid"], task["name"],
                        task["group"], task["role"])
        LOG.warning("Investigate each: the taxid may be a parent with no assembly of "
                    "its own, newly renamed, or genuinely unsequenced. Rerunning is "
                    "cheap - completed taxa are skipped.")

    try:
        LOG.info("disk used by refs:    %s",
                 _human(sum(f.stat().st_size for f in refs_dir.rglob("*") if f.is_file())))
    except OSError:
        pass
    if duplicate_neighbours:
        LOG.warning("")
        LOG.warning("%d competitive bundle(s) lost a neighbour because it resolved to "
                    "the target's own assembly. Those pairs CANNOT be discriminated "
                    "until a distinct assembly is curated for the neighbour:",
                    len(duplicate_neighbours))
        for tid, name, ntids in duplicate_neighbours:
            LOG.warning("  %s (%s): lost neighbour(s) %s", tid, name, ", ".join(ntids))

    LOG.info("free space remaining: %s", _human(shutil.disk_usage(btf_base).free))
    LOG.info("===============================================================")

    # Non-zero if any target cannot be confirmed at all: that is a pack-level
    # defect a human must look at, and a CI/Nextflow wrapper should notice.
    return 1 if unbuildable else 0


def _human(nbytes: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if nbytes < 1024 or unit == "TiB":
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} TiB"


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        LOG.warning("interrupted - rerun to resume (completed downloads are skipped)")
        sys.exit(130)
