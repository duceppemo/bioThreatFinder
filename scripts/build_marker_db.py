#!/usr/bin/env python3
"""
build_marker_db.py - resolve the marker catalogue of a bioThreatFinder domain
                     pack to real sequences and build the marker search
                     databases used by the marker-evidence tier.

Why this exists
---------------
packs/<pack>/markers.tsv is keyed by GENE SYMBOL, not by accession, on purpose:
an accession written from memory is a fabrication risk in a file that gates a
high-consequence call. This script is the only place where a symbol becomes a
sequence. It resolves each catalogue row over NCBI E-utilities, records the
accession / query / retrieval date it actually used, and - crucially - records
the rows it could NOT resolve.

That last point is the whole design constraint. Downstream (bin/marker_evidence.py,
bin/marker_rules.py) uses three-valued logic:

    PRESENT   marker was searched for and found
    ABSENT    marker was searched for and not found
    UNKNOWN   marker was never testable, because no sequence was ever built

"ABSENT" and "UNKNOWN" must never be conflated. A pXO1 marker that is absent
from a sample argues against anthrax; a pXO1 marker that was never built tells
you nothing at all, and treating it as absence would silently downgrade a real
finding. Equally, treating a build failure as a clean absence for an
`absence_expected` marker would silently UPGRADE a false positive. So every row
of the catalogue appears in markers.manifest.tsv with an explicit status, and
rows with no usable sequence carry a non-`built` status and empty paths.

What it writes (paths must match packs/<pack>/pack.yaml)
-------------------------------------------------------
  <db_root>/bioThreatFinder/packs/<pack>/markers/
      nt/<marker_id>.fna            per-marker nucleotide sequences
      aa/<marker_id>.faa            per-marker protein sequences
      hmm/<marker_id>.hmm           per-marker profile HMM
      markers.nt.fna                all nucleotide sequences, concatenated
      markers.nt.sr.mmi             minimap2 index, short reads
      markers.nt.map-ont.mmi        minimap2 index, ONT reads
      markers.aa.faa                all protein sequences, concatenated
      mmseqs/markers_aa*            MMseqs2 database built from markers.aa.faa
      markers.hmm                   all profiles, concatenated + hmmpress'ed
      markers.manifest.tsv          THE record of tested vs untested
      curation_required.tsv         unmet curation tasks (see below)

THE HEADER CONVENTION IS LOAD-BEARING
-------------------------------------
Every sequence written by this script - in the per-marker files and in the
combined files alike - is renamed to

    >marker|<marker_id>|<accession>|<original description>

'|'-delimited, marker_id at field index 1. bin/marker_evidence.py attributes an
alignment or a translated hit to a marker by splitting the reference name on
'|' and taking field 1, with no side lookup table to fall out of sync. The
original description is kept verbatim as the LAST field and may itself contain
'|' characters (NCBI descriptions do), so the parse must split with maxsplit=3
and never with a plain split(). make_marker_header() / parse_marker_header()
below are the only two functions allowed to know this layout, and --self-test
proves the round trip including pipe-bearing descriptions.

Profile HMMs get NAME 'marker|<marker_id>' for the same reason: the HMM name is
what appears in the hmmsearch --domtblout query column, so the same
split('|')[1] rule works for the HMM tier.

Retrieval
---------
NCBI E-utilities over HTTPS via urllib - Python standard library only, no
Biopython, no third-party imports anywhere in this file, because this runs
inside a minimal container on the cluster. Raw responses are cached under
<db_root>/bioThreatFinder/cache/entrez/ so reruns are cheap and resumable;
see ENTREZ_CACHE_TTL_*_DAYS for the TTLs and the reasoning behind them.

Dependencies
------------
Python 3.8+ standard library only, plus these external command-line tools,
each resolved at runtime and reported with its install command if missing:
  * `hmmbuild`/`hmmpress` (HMMER 3)  - profile tier      (--skip-hmm to omit)
  * `mmseqs`             (MMseqs2)   - protein tier      (--skip-protein to omit)
  * `minimap2`                       - nucleotide indexes (warn-only if absent)
  * `mafft`/`muscle`/`famsa`/`clustalo` - MSA for hmmbuild; a stdlib
    centre-star fallback is used if none is present (see align_protein())
  * `seqkit`                         - optional; only used for a dedup fast path
  * `datasets` (NCBI datasets CLI)   - not used here; reported for parity with
    scripts/build_pack_refs.py, which does need it

Usage
-----
  build_marker_db.py --self-test
  build_marker_db.py --pack human --dry-run
  build_marker_db.py --pack human --email you@example.org --api-key $NCBI_API_KEY
  build_marker_db.py --pack human --only-taxon 'Bacillus anthracis' -v
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import http.client
import json
import logging
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

LOG = logging.getLogger("build_marker_db")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_DB_ROOT = Path("/media/36tb/db")
CONDA_BIN = Path("/home/bioinfo/miniconda3/bin")

# The shared volume runs at ~97% capacity. A marker database is small (tens of
# MB of sequence; the minimap2 indexes and the MMseqs2 database are the bulk),
# but refuse to start on a nearly-full volume rather than die halfway through
# and leave truncated indexes that look valid.
DEFAULT_MIN_FREE_GIB = 20

# VERIFY: E-utilities base URL and endpoint names. Long-standing and stable,
# but re-check against the current E-utilities documentation if a request
# starts returning HTML or 404.
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENTREZ_TOOL = "bioThreatFinder-build_marker_db"

# NCBI's published limits: 3 requests/second anonymously, 10/second with an API
# key. Exceeding them gets the whole host blocked, which would break every
# other pipeline on this box, so the limiter is unconditional and shared.
# VERIFY: rate limits are policy, not API, and NCBI has changed them before.
RATE_LIMIT_NO_KEY = 3.0
RATE_LIMIT_WITH_KEY = 10.0

# Cache TTLs. A sequence record fetched by versioned accession is immutable, so
# efetch responses can be cached for a year; esearch results change as GenBank
# grows, so searches expire in a month. Both are overridden by --force.
ENTREZ_CACHE_TTL_SEARCH_DAYS = 30
ENTREZ_CACHE_TTL_FETCH_DAYS = 365

HTTP_TIMEOUT = 120
HTTP_RETRIES = 5
# Retry on transient server-side conditions only. A 400 means the query itself
# is wrong and retrying just burns the rate limit.
HTTP_RETRY_STATUS = {429, 500, 502, 503, 504}

# How many sequences to keep per marker. A single-copy gene in one organism
# needs a couple of alleles; a divergent family (BoNT, Stx, filovirus GP, CCHFV
# L) needs enough members that the profile spans the family's real diversity.
MAX_SEQS_SINGLE = 3
MAX_SEQS_FAMILY = 8
# esearch retmax. Kept small deliberately: we fetch candidates and then filter
# hard on length and gene name, and a huge candidate list mostly adds latency.
ESEARCH_RETMAX = 40

# A marker must be a locus, not a chromosome. Anything longer than this is only
# accepted after CDS extraction; see fetch_nucleotide().
NT_MAX_DIRECT_LEN = 50_000
NT_MIN_LEN = 60
AA_MAX_LEN = 10_000
AA_MIN_LEN = 30

# Fewer than three sequences cannot express a family's conserved columns: a
# profile built from one or two sequences is just that sequence with a cosmetic
# weight matrix, and it would report high-confidence hits with no statistical
# basis. Such a marker is recorded `partial` instead - the protein and
# nucleotide tiers still work, and the manifest says the HMM tier does not.
MIN_SEQS_FOR_HMM = 3

# Cap on the stdlib fallback aligner. Needleman-Wunsch in pure Python is
# O(n*m) per pair; beyond this the run would stall for minutes per marker, and
# the honest answer is "install an aligner" rather than a silent hang.
MAX_NW_CELLS = 4_000_000

MANIFEST_COLUMNS = [
    "marker_id", "taxon", "gene", "marker_class", "evidence_type",
    "query_strategy", "status", "n_nt_seqs", "n_aa_seqs", "has_hmm",
    "accessions", "nt_path", "aa_path", "hmm_path", "md5", "retrieved_utc",
    "note",
]

CATALOGUE_COLUMNS = [
    "marker_id", "taxon", "gene", "aliases", "marker_class", "evidence_type",
    "locus", "min_coverage", "min_identity", "copy_expectation",
    "query_strategy", "source_db", "notes",
]

# The status taxonomy. This is the contract with bin/marker_evidence.py: only
# `built` (and the tiers a `partial` row actually produced) may be evaluated as
# PRESENT/ABSENT; everything else is UNKNOWN downstream.
STATUS_BUILT = "built"
STATUS_PARTIAL = "partial"
STATUS_CURATION = "curation_required"
STATUS_FAILED = "failed"
STATUS_MOTIF = "motif_todo"
ALL_STATUSES = [STATUS_BUILT, STATUS_PARTIAL, STATUS_CURATION,
                STATUS_FAILED, STATUS_MOTIF]

KNOWN_STRATEGIES = {"gene_in_organism", "gene_family", "named_region", "motif",
                    "curate_required"}
KNOWN_TIERS = {"dna", "protein", "hmm"}

MM2_PRESETS = ["sr", "map-ont"]

# Install hints. The cluster runs containers, so nothing is assumed to be on
# PATH. VERIFY: biocontainer image tags change constantly - pick a current tag
# from https://quay.io/organization/biocontainers rather than trusting a tag
# written here, which is why none is written here.
TOOL_INSTALL = {
    "hmmbuild": f"{CONDA_BIN}/mamba install -c bioconda hmmer"
                "   # or: singularity pull docker://quay.io/biocontainers/hmmer:<tag>",
    "hmmpress": f"{CONDA_BIN}/mamba install -c bioconda hmmer"
                "   # or: singularity pull docker://quay.io/biocontainers/hmmer:<tag>",
    "mmseqs": f"{CONDA_BIN}/mamba install -c bioconda mmseqs2"
              "   # or: singularity pull docker://quay.io/biocontainers/mmseqs2:<tag>",
    "minimap2": f"{CONDA_BIN}/mamba install -c bioconda minimap2"
                "   # or: singularity pull docker://quay.io/biocontainers/minimap2:<tag>",
    "seqkit": f"{CONDA_BIN}/mamba install -c bioconda seqkit"
              "   # or: singularity pull docker://quay.io/biocontainers/seqkit:<tag>",
    "mafft": f"{CONDA_BIN}/mamba install -c bioconda mafft"
             "   # or: singularity pull docker://quay.io/biocontainers/mafft:<tag>",
    "muscle": f"{CONDA_BIN}/mamba install -c bioconda muscle",
    "famsa": f"{CONDA_BIN}/mamba install -c bioconda famsa",
    "clustalo": f"{CONDA_BIN}/mamba install -c bioconda clustalo",
    "datasets": f"{CONDA_BIN}/mamba install -c conda-forge ncbi-datasets-cli",
}

ALIGNERS = ["mafft", "famsa", "muscle", "clustalo"]


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sanitise(name: str) -> str:
    """Filesystem-safe stem. Marker ids are already tame; taxa are not."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", (name or "").strip())
    return re.sub(r"_+", "_", s).strip("_") or "unknown"


def md5_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def combined_md5(paths: Sequence[Path]) -> str:
    """One digest over a marker's whole output set.

    Digest of "<basename>:<md5>" lines rather than of the concatenated bytes, so
    that adding a tier changes the digest even if the existing files do not -
    which is what makes the manifest md5 usable as a resume check.
    """
    h = hashlib.md5()
    for p in sorted(paths, key=lambda x: x.name):
        if p.is_file():
            h.update(f"{p.name}:{md5_of(p)}\n".encode("utf-8"))
    return h.hexdigest()


def which(name: str) -> Optional[str]:
    """PATH first, then the shared conda prefix (not always on PATH under
    Nextflow/cron), then the minimap2 bundled with the ONT dorado server.

    Same resolution order as scripts/build_pack_refs.py on purpose: two builders
    that disagree about which minimap2 they used produce indexes that cannot be
    reasoned about together.
    """
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


def human(nbytes: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if nbytes < 1024 or unit == "TiB":
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} TiB"


# ---------------------------------------------------------------------------
# THE HEADER CONVENTION
# ---------------------------------------------------------------------------
# Only these two functions may know the layout. Everything else calls them.

MARKER_TAG = "marker"
MARKER_HEADER_FIELDS = 4  # marker | marker_id | accession | description


def sanitise_accession(accession: str) -> str:
    """Make an accession safe for the pipe-delimited header without losing it.

    NCBI protein fetches return SwissProt-style names such as
    `sp|Q306W5.1|POLS_EEEV1`, whose pipes collide with the header's field delimiter.
    Raising on them, as this module first did, turned a routine identifier format into
    a hard retrieval failure for every affected marker.

    For `sp|ACC|ENTRY` and `tr|ACC|ENTRY` the middle field IS the accession, so it is
    extracted rather than mangled. Anything else has its pipes replaced, which keeps
    the identifier readable and the contract intact.
    """
    acc = (accession or "").strip()
    m = re.match(r"^(?:sp|tr|pdb|gb|emb|dbj|ref|pir|prf|tpg|tpe|tpd)\|([^|]+)\|", acc)
    if m:
        return m.group(1)
    return acc.replace("|", "_")


def make_marker_header(marker_id: str, accession: str, description: str) -> str:
    """Build the load-bearing sequence name (without the leading '>').

    marker_id and accession must not contain '|' - they are positional fields.
    The description is the final field and is kept verbatim, pipes included.
    """
    if "|" in marker_id:
        raise ValueError(f"marker_id must not contain '|': {marker_id!r}")
    accession = sanitise_accession(accession)
    # Newlines/CR would split a FASTA record in two; tabs would break any
    # downstream TSV that carries the reference name.
    desc = re.sub(r"[\r\n\t]+", " ", description or "").strip()
    return f"{MARKER_TAG}|{marker_id}|{accession}|{desc}"


def parse_marker_header(name: str) -> Optional[Tuple[str, str, str]]:
    """Inverse of make_marker_header(): -> (marker_id, accession, description).

    Returns None for anything that is not one of ours, so a mixed index cannot
    silently attribute a foreign reference to a marker.

    maxsplit is essential: NCBI descriptions contain '|' (e.g. the old
    gi|...|ref|... style and many strain names), and a plain split() would
    truncate the description or, worse, shift the fields.

    Accepts the 2-field form 'marker|<marker_id>' as well, because that is the
    NAME given to each profile HMM and it appears in hmmsearch output.
    """
    if not name:
        return None
    text = name[1:] if name.startswith(">") else name
    # A FASTA name ends at the first whitespace, but our description is part of
    # the header line, so take the whole line and strip trailing whitespace.
    text = text.rstrip()
    parts = text.split("|", MARKER_HEADER_FIELDS - 1)
    if len(parts) < 2 or parts[0] != MARKER_TAG:
        return None
    marker_id = parts[1]
    if not marker_id:
        return None
    accession = parts[2] if len(parts) > 2 else ""
    description = parts[3] if len(parts) > 3 else ""
    return marker_id, accession, description


def marker_id_of(name: str) -> Optional[str]:
    """Field index 1 of the '|'-delimited name - what marker_evidence.py does."""
    parsed = parse_marker_header(name)
    return parsed[0] if parsed else None


# ---------------------------------------------------------------------------
# Minimal FASTA I/O
# ---------------------------------------------------------------------------

class SeqRecord:
    __slots__ = ("accession", "description", "seq", "gene")

    def __init__(self, accession: str, description: str, seq: str, gene: str = ""):
        self.accession = accession
        self.description = description
        self.seq = seq
        self.gene = gene

    def __len__(self) -> int:
        return len(self.seq)

    def __repr__(self) -> str:
        return f"<SeqRecord {self.accession} {len(self.seq)}aa/nt>"


CDS_ATTR_RE = re.compile(r"\[([A-Za-z_]+)=([^\]]*)\]")


def parse_fasta(text: str) -> List[SeqRecord]:
    """Parse NCBI FASTA text, including the CDS flavours.

    Plain records look like '>NC_000962.3 Mycobacterium tuberculosis ...'.
    rettype=fasta_cds_na/_aa records look like
    '>lcl|NC_000962.3_cds_NP_214515.1_1 [gene=dnaA] [protein_id=NP_214515.1] ...'
    so the accession is taken from [protein_id=] where present and the gene name
    from [gene=]/[locus_tag=], which is what lets a CDS be matched to a marker.
    """
    records: List[SeqRecord] = []
    name: Optional[str] = None
    chunks: List[str] = []

    def flush() -> None:
        if name is None:
            return
        seq = "".join(chunks).replace(" ", "")
        if not seq:
            return
        first, _, rest = name.partition(" ")
        attrs = dict(CDS_ATTR_RE.findall(rest))
        gene = attrs.get("gene") or attrs.get("locus_tag") or ""
        accession = attrs.get("protein_id") or ""
        if not accession:
            token = first
            if token.startswith("lcl|"):
                token = token[4:]
                # 'NC_000962.3_cds_NP_214515.1_1' -> keep the source accession
                token = token.split("_cds_")[0]
            accession = token
        records.append(SeqRecord(accession=accession, description=rest.strip(),
                                 seq=seq.upper(), gene=gene))

    for line in text.splitlines():
        if line.startswith(">"):
            flush()
            name = line[1:].strip()
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    flush()
    return records


def write_marker_fasta(path: Path, marker_id: str, records: Sequence[SeqRecord],
                       wrap: int = 60) -> None:
    """Write records under the load-bearing header convention, atomically."""
    tmp = path.with_name(path.name + ".part")
    with tmp.open("w", encoding="utf-8", newline="\n") as fh:
        for rec in records:
            fh.write(">" + make_marker_header(marker_id, rec.accession,
                                              rec.description) + "\n")
            seq = rec.seq
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")
    tmp.replace(path)


def dedup_records(records: Iterable[SeqRecord]) -> List[SeqRecord]:
    """Drop exact duplicate sequences, keeping the first (best-ranked) one.

    Exact duplicates are the norm - a gene is deposited hundreds of times - and
    they would bias an HMM towards whichever allele happens to be popular.
    (`seqkit rmdup -s` does the same thing; it is not used here because the
    sequence sets are small enough that a dict is faster than a subprocess, and
    it keeps this step working in a container without seqkit.)
    """
    seen = set()
    out: List[SeqRecord] = []
    for rec in records:
        key = rec.seq
        if key in seen:
            continue
        seen.add(key)
        out.append(rec)
    return out


# ---------------------------------------------------------------------------
# Catalogue + watchlist
# ---------------------------------------------------------------------------

def read_marker_catalogue(path: Path) -> List[dict]:
    """Read packs/<pack>/markers.tsv.

    The file is heavily commented; the header is the last comment-free line
    starting with 'marker_id', and comment lines are interleaved with data rows
    (they are section banners), so filter them out before csv sees the file.
    """
    if not path.is_file():
        raise SystemExit(
            f"ERROR: marker catalogue not found: {path}\n"
            f"       Expected packs/<pack>/markers.tsv. A pack is data; this script "
            f"does not invent markers."
        )
    with path.open(encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    if not lines:
        raise SystemExit(f"ERROR: {path} contains no non-comment lines")
    reader = csv.DictReader(lines, delimiter="\t")
    missing = [c for c in CATALOGUE_COLUMNS if c not in (reader.fieldnames or [])]
    if missing:
        raise SystemExit(
            f"ERROR: {path} is missing required column(s): {', '.join(missing)}\n"
            f"       found: {', '.join(reader.fieldnames or [])}")
    rows: List[dict] = []
    seen_ids = set()
    for row in reader:
        marker_id = (row.get("marker_id") or "").strip()
        if not marker_id:
            continue
        if "|" in marker_id:
            raise SystemExit(
                f"ERROR: marker_id {marker_id!r} contains '|', which is the field "
                f"separator of the sequence-name convention. Rename it in {path}.")
        if marker_id in seen_ids:
            raise SystemExit(f"ERROR: duplicate marker_id in {path}: {marker_id}")
        seen_ids.add(marker_id)
        rows.append({k: (v or "").strip() for k, v in row.items() if k})
    if not rows:
        raise SystemExit(f"ERROR: {path} has no data rows")
    return rows


def read_taxon_categories(pack_dir: Path) -> Dict[str, str]:
    """taxon name (lowercased) -> watchlist category.

    Needed only for the exit code: a failed `required` marker for a
    high_consequence taxon is a build defect a human must see, whereas a failed
    lineage marker for a near neighbour is a nuisance. The resolved watchlist is
    preferred but the raw one carries the same `category` column, so either
    works and neither is a hard requirement.
    """
    out: Dict[str, str] = {}
    for candidate in ("watchlist.resolved.tsv", "watchlist.tsv"):
        path = pack_dir / candidate
        if not path.is_file():
            continue
        with path.open(encoding="utf-8") as fh:
            lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
        reader = csv.DictReader(lines, delimiter="\t")
        if "name" not in (reader.fieldnames or []) or "category" not in (reader.fieldnames or []):
            continue
        for row in reader:
            name = (row.get("name") or "").strip().lower()
            if name:
                out.setdefault(name, (row.get("category") or "").strip().lower())
        if out:
            LOG.info("read %d taxon categories from %s", len(out), path)
            break
    if not out:
        LOG.warning("no watchlist category information found under %s; every marker "
                    "failure will be treated as non-fatal", pack_dir)
    return out


def parse_tiers(evidence_type: str) -> List[str]:
    """'dna+protein+hmm' -> ['dna', 'protein', 'hmm'], order preserved."""
    tiers = []
    for tok in (evidence_type or "").split("+"):
        tok = tok.strip().lower()
        if tok and tok in KNOWN_TIERS and tok not in tiers:
            tiers.append(tok)
    return tiers


def parse_aliases(raw: str) -> List[str]:
    """The aliases column is ';'-separated, with '-' meaning none."""
    out = []
    for tok in re.split(r"[;,]", raw or ""):
        tok = tok.strip()
        if not tok or tok == "-":
            continue
        if tok not in out:
            out.append(tok)
    return out


# ---------------------------------------------------------------------------
# Entrez client
# ---------------------------------------------------------------------------

class EntrezError(Exception):
    pass


class Entrez:
    """Rate-limited, cached, backing-off E-utilities client (stdlib only).

    Caching is by request, not by marker, so that the many markers that search
    the same organism/gene space reuse each other's responses and an interrupted
    run costs nothing to resume.
    """

    def __init__(self, cache_dir: Path, email: str = "", api_key: str = "",
                 force: bool = False, offline: bool = False):
        self.cache_dir = cache_dir
        self.email = email or ""
        self.api_key = api_key or ""
        self.force = force
        self.offline = offline
        self.min_interval = 1.0 / (RATE_LIMIT_WITH_KEY if self.api_key else RATE_LIMIT_NO_KEY)
        self._last_request = 0.0
        self.n_requests = 0
        self.n_cache_hits = 0
        if not self.api_key:
            LOG.info("no NCBI API key supplied: limiting to %.0f request(s)/second. "
                     "An API key raises this to %.0f and makes a full pack build "
                     "roughly three times faster.",
                     RATE_LIMIT_NO_KEY, RATE_LIMIT_WITH_KEY)

    # -- cache ------------------------------------------------------------

    def _cache_path(self, endpoint: str, params: Dict[str, str]) -> Path:
        # The key deliberately excludes api_key/email/tool: they are credentials
        # and courtesy fields, not part of the query, and including them would
        # give every operator a private cache of identical responses.
        payload = json.dumps({k: v for k, v in sorted(params.items())
                              if k not in ("api_key", "email", "tool")},
                             sort_keys=True)
        key = hashlib.md5(f"{endpoint}?{payload}".encode("utf-8")).hexdigest()
        return self.cache_dir / endpoint / key[:2] / f"{key}.txt"

    def _read_cache(self, path: Path, ttl_days: int) -> Optional[str]:
        if self.force or not path.is_file() or path.stat().st_size == 0:
            return None
        age_days = (time.time() - path.stat().st_mtime) / 86400.0
        if age_days > ttl_days:
            LOG.debug("cache entry %s is %.1f days old (TTL %d); refetching",
                      path.name, age_days, ttl_days)
            return None
        self.n_cache_hits += 1
        return path.read_text(encoding="utf-8", errors="replace")

    def _write_cache(self, path: Path, url: str, body: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".part")
        tmp.write_text(body, encoding="utf-8")
        tmp.replace(path)
        # Sidecar provenance: which URL produced this blob and when. Without it a
        # cache directory is unauditable, and provenance is part of the result.
        meta = path.with_suffix(".meta.json")
        meta.write_text(json.dumps({"url": url, "retrieved_utc": utc_now()},
                                   indent=2), encoding="utf-8")

    # -- transport --------------------------------------------------------

    def _throttle(self) -> None:
        wait = self.min_interval - (time.monotonic() - self._last_request)
        if wait > 0:
            time.sleep(wait)
        self._last_request = time.monotonic()

    def _request(self, endpoint: str, params: Dict[str, str], ttl_days: int) -> str:
        cache_path = self._cache_path(endpoint, params)
        cached = self._read_cache(cache_path, ttl_days)
        if cached is not None:
            return cached
        if self.offline:
            raise EntrezError(f"offline and no cached response for {endpoint} "
                              f"{params.get('term') or params.get('id') or ''}")

        full = dict(params)
        full["tool"] = ENTREZ_TOOL
        if self.email:
            full["email"] = self.email
        if self.api_key:
            full["api_key"] = self.api_key
        url = f"{ENTREZ_BASE}/{endpoint}?" + urllib.parse.urlencode(full)

        last_err = ""
        for attempt in range(HTTP_RETRIES):
            self._throttle()
            self.n_requests += 1
            try:
                req = urllib.request.Request(url, headers={"User-Agent": ENTREZ_TOOL})
                with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as resp:
                    body = resp.read().decode("utf-8", errors="replace")
                if body.strip():
                    self._write_cache(cache_path, url, body)
                    return body
                last_err = "empty response body"
            except urllib.error.HTTPError as exc:
                last_err = f"HTTP {exc.code}"
                if exc.code not in HTTP_RETRY_STATUS:
                    # 400-class: the query is malformed or the id is unknown.
                    # Retrying is pointless and costs the rate-limit budget.
                    raise EntrezError(f"{endpoint}: HTTP {exc.code} (not retryable)")
            except (urllib.error.URLError, OSError, TimeoutError,
                    http.client.HTTPException) as exc:
                # http.client.IncompleteRead lives under HTTPException, not OSError, so
                # without it here a truncated efetch body escaped the retry loop
                # entirely and surfaced as "unexpected error". That cost 30+ markers on
                # the first full run - transient network truncation recorded as a
                # permanent retrieval failure.
                last_err = f"{type(exc).__name__}: {exc}"
            # Exponential backoff with jitter: a whole cluster retrying in
            # lockstep is how a transient 503 becomes a sustained one.
            delay = min(60.0, (2 ** attempt)) * (1.0 + random.random() * 0.3)
            LOG.warning("%s: %s - retrying in %.1fs (attempt %d/%d)",
                        endpoint, last_err, delay, attempt + 1, HTTP_RETRIES)
            time.sleep(delay)
        raise EntrezError(f"{endpoint}: giving up after {HTTP_RETRIES} attempts ({last_err})")

    # -- API --------------------------------------------------------------

    def esearch(self, db: str, term: str, retmax: int = ESEARCH_RETMAX) -> List[str]:
        """Return accession-agnostic UIDs, best-first as NCBI ranked them."""
        params = {"db": db, "term": term, "retmax": str(retmax), "retmode": "json"}
        try:
            body = self._request("esearch.fcgi", params, ENTREZ_CACHE_TTL_SEARCH_DAYS)
        except EntrezError as exc:
            LOG.debug("esearch failed (%s): %s", term, exc)
            return []
        try:
            data = json.loads(body)
        except json.JSONDecodeError:
            LOG.debug("esearch returned non-JSON for %r", term)
            return []
        result = data.get("esearchresult") or {}
        if result.get("ERROR"):
            LOG.debug("esearch error for %r: %s", term, result["ERROR"])
            return []
        return [str(u) for u in (result.get("idlist") or [])]

    def efetch(self, db: str, uids: Sequence[str], rettype: str) -> str:
        """efetch as text. POST is not used: GET keeps the cache key simple and
        our id lists are short by construction (ESEARCH_RETMAX).

        VERIFY: rettype values. 'fasta' is universal; 'fasta_cds_na' and
        'fasta_cds_aa' are nuccore-only and are what give per-CDS sequences with
        [gene=]/[protein_id=] attributes. Re-check the E-utilities rettype table
        if a CDS fetch starts coming back empty for annotated records.
        """
        if not uids:
            return ""
        params = {"db": db, "id": ",".join(uids), "rettype": rettype,
                  "retmode": "text"}
        try:
            return self._request("efetch.fcgi", params, ENTREZ_CACHE_TTL_FETCH_DAYS)
        except EntrezError as exc:
            LOG.debug("efetch failed (%s %s): %s", db, rettype, exc)
            return ""


# ---------------------------------------------------------------------------
# Query construction
# ---------------------------------------------------------------------------
# Queries are built as an ordered list of progressively looser attempts. The
# first attempt that yields usable sequence wins, and whichever attempt was used
# is recorded in the manifest note - "we found pagA" is not a result, "we found
# pagA with this query on this date" is.

MIN_QUERY_NAME_LEN = 2
# Amino-acid substitution designations (L98H, Y132F, S639F) and tandem-repeat
# labels (TR34, TR46) appear in the aliases column because they are what makes
# the marker actionable, but they are not searchable names: as an Entrez term
# they match unrelated records or nothing at all. Dropped for querying only -
# they stay in the catalogue, and the rows that hinge on them are motif rows.
SUBSTITUTION_RE = re.compile(r"^(TR\d{1,4}|[A-Z]\d{1,4}[A-Z*]|\d{1,4}[A-Z])$")


def query_names(row: dict) -> List[str]:
    """Searchable names for a row: the gene symbol plus usable aliases.

    A one-character alias ('t' for tul4) would match essentially everything in
    [All Fields] and flood the candidate list with unrelated records, so short
    tokens are dropped too.
    """
    gene = row.get("gene", "")
    raw = ([gene] if gene and gene != "-" else []) + parse_aliases(row.get("aliases", ""))
    out: List[str] = []
    for name in raw:
        if len(name) < MIN_QUERY_NAME_LEN:
            continue
        if SUBSTITUTION_RE.match(name):
            continue
        if name not in out:
            out.append(name)
    return out


def organism_clause(taxon: str) -> str:
    return f'"{taxon}"[Organism]'


def genus_of(taxon: str) -> str:
    return (taxon or "").split()[0] if taxon else ""


def nt_queries(row: dict) -> List[str]:
    """Nucleotide queries for one catalogue row, strictest first."""
    taxon, gene = row["taxon"], row["gene"]
    strategy = row["query_strategy"]
    names = query_names(row)
    org = organism_clause(taxon)
    out: List[str] = []

    if strategy == "named_region":
        # A named region (RD1, RD9, TTS1, wboA, LIPI-4, eri) is usually not a
        # gene symbol, so [Gene Name] would return nothing. Search titles and
        # all fields, and accept that some of these genuinely need curated
        # coordinates - hence the VERIFY notes in the catalogue itself.
        for name in names or [gene]:
            out.append(f'{name}[Title] AND {org}')
            out.append(f'{name}[All Fields] AND {org}')
        return dedup_str(out)

    for name in names:
        # RefSeq first: curated, one record per locus, stable accession.
        out.append(f'{name}[Gene Name] AND {org} AND srcdb_refseq[Properties]')
        out.append(f'{name}[Gene Name] AND {org}')
        out.append(f'{name}[All Fields] AND {org}')
    if strategy == "gene_family" and names:
        # Divergence is the point of a family marker: a member that is only
        # deposited under a sister species (BoNT/E in C. butyricum, Stx in
        # Shigella) still belongs in the family, so relax the organism to genus.
        genus = genus_of(taxon)
        if genus:
            for name in names:
                out.append(f'{name}[All Fields] AND "{genus}"[Organism]')
    return dedup_str(out)


def aa_queries(row: dict) -> List[str]:
    taxon, gene = row["taxon"], row["gene"]
    strategy = row["query_strategy"]
    names = query_names(row)
    org = organism_clause(taxon)
    out: List[str] = []
    for name in names:
        out.append(f'{name}[Gene Name] AND {org} AND srcdb_refseq[Properties]')
        out.append(f'{name}[Protein Name] AND {org} AND srcdb_refseq[Properties]')
        out.append(f'{name}[Gene Name] AND {org}')
        out.append(f'{name}[All Fields] AND {org}')
    if strategy == "gene_family" and names:
        genus = genus_of(taxon)
        if genus:
            for name in names:
                out.append(f'{name}[Protein Name] AND "{genus}"[Organism]')
                out.append(f'{name}[All Fields] AND "{genus}"[Organism]')
    if strategy == "named_region":
        for name in names or [gene]:
            out.append(f'{name}[All Fields] AND {org}')
    return dedup_str(out)


def dedup_str(items: Iterable[str]) -> List[str]:
    out: List[str] = []
    for i in items:
        if i and i not in out:
            out.append(i)
    return out


def gene_name_matches(cds_gene: str, wanted: Sequence[str]) -> bool:
    """Case-insensitive match of a CDS [gene=] value against gene + aliases.

    Deliberately permissive at the prefix level ('capB' matches 'capBCADE'
    operon annotations, 'stx2' matches 'stx2A'), because NCBI annotation of the
    same locus varies between submitters. The cost of a loose match here is a
    slightly larger marker set; the cost of a strict one is an empty one.
    """
    g = (cds_gene or "").strip().lower()
    if not g:
        return False
    for w in wanted:
        w = (w or "").strip().lower()
        if not w:
            continue
        if g == w or g.startswith(w) or w.startswith(g):
            return True
    return False


# ---------------------------------------------------------------------------
# Retrieval
# ---------------------------------------------------------------------------

def fetch_protein(ez: Entrez, row: dict, want: int) -> Tuple[List[SeqRecord], str]:
    """Protein records for a marker. Returns (records, query-used)."""
    for query in aa_queries(row):
        uids = ez.esearch("protein", query)
        if not uids:
            continue
        text = ez.efetch("protein", uids[:max(want * 4, want)], "fasta")
        records = [r for r in parse_fasta(text) if AA_MIN_LEN <= len(r) <= AA_MAX_LEN]
        if not records:
            continue
        records = dedup_records(records)
        # Prefer RefSeq accessions (NP_/YP_/WP_/XP_) - curated, and stable across
        # rebuilds, which matters because the accession goes in the manifest.
        records.sort(key=lambda r: (0 if re.match(r"^[NYWX]P_", r.accession) else 1,
                                    -len(r)))
        LOG.debug("  protein: %d record(s) via %s", len(records), query)
        return records[:want], query
    return [], (aa_queries(row)[0] if aa_queries(row) else "")


def fetch_nucleotide(ez: Entrez, row: dict, want: int) -> Tuple[List[SeqRecord], str]:
    """Nucleotide records for a marker. Returns (records, query-used).

    Two routes, because "the gene" is not a record in GenBank:
      1. CDS extraction (rettype=fasta_cds_na) from whichever records the search
         found, keeping only the CDSs whose [gene=] matches. This is what makes
         a marker a locus rather than a chromosome, and it works even when the
         only deposit is a 5 Mb complete genome.
      2. Plain FASTA, accepted only when the record is short enough to BE a
         locus (NT_MAX_DIRECT_LEN). This is the route for the things that have
         no CDS at all: 18S rRNA, IS elements, named regions, promoters.
    """
    # Same filtered name list as the queries use: matching a CDS on a
    # one-character alias would prefix-match a large fraction of every genome.
    wanted_names = query_names(row)
    strategy = row["query_strategy"]
    for query in nt_queries(row):
        uids = ez.esearch("nuccore", query)
        if not uids:
            continue
        probe = uids[:max(want * 3, want)]
        records: List[SeqRecord] = []

        # Route 1 - CDS extraction. Skipped for named_region/no-gene rows, where
        # there is no gene symbol to match a CDS against.
        if strategy != "named_region" and wanted_names:
            cds_text = ez.efetch("nuccore", probe, "fasta_cds_na")
            for rec in parse_fasta(cds_text):
                if not gene_name_matches(rec.gene, wanted_names):
                    continue
                if NT_MIN_LEN <= len(rec) <= NT_MAX_DIRECT_LEN:
                    records.append(rec)

        # Route 2 - whole record, only if it is locus-sized.
        if not records:
            text = ez.efetch("nuccore", probe, "fasta")
            for rec in parse_fasta(text):
                if NT_MIN_LEN <= len(rec) <= NT_MAX_DIRECT_LEN:
                    records.append(rec)
                else:
                    LOG.debug("  nt: rejecting %s (%d bp; not locus-sized)",
                              rec.accession, len(rec))

        if not records:
            continue
        records = dedup_records(records)
        records.sort(key=lambda r: (0 if re.match(r"^(N[CGMRTWZ]|X[MR])_", r.accession) else 1,
                                    -len(r)))
        LOG.debug("  nt: %d record(s) via %s", len(records), query)
        return records[:want], query
    return [], (nt_queries(row)[0] if nt_queries(row) else "")


# ---------------------------------------------------------------------------
# Alignment + HMM
# ---------------------------------------------------------------------------

def nw_align(a: str, b: str, match: int = 2, mismatch: int = -1,
             gap: int = -2) -> Tuple[str, str]:
    """Global Needleman-Wunsch with a linear gap penalty.

    Used only by the stdlib fallback aligner. Scoring is deliberately generic
    (not BLOSUM): the fallback exists to produce a usable MSA when no real
    aligner is installed, not to compete with one.
    """
    n, m = len(a), len(b)
    prev = [gap * j for j in range(m + 1)]
    # Traceback stored as one bytearray row per query position: 0=diag, 1=up, 2=left
    pointers: List[bytearray] = []
    for i in range(1, n + 1):
        cur = [gap * i] + [0] * m
        ptr = bytearray(m + 1)
        ptr[0] = 1
        ai = a[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (match if ai == b[j - 1] else mismatch)
            up = prev[j] + gap
            left = cur[j - 1] + gap
            best = diag
            code = 0
            if up > best:
                best, code = up, 1
            if left > best:
                best, code = left, 2
            cur[j] = best
            ptr[j] = code
        pointers.append(ptr)
        prev = cur

    out_a: List[str] = []
    out_b: List[str] = []
    i, j = n, m
    while i > 0 and j > 0:
        code = pointers[i - 1][j]
        if code == 0:
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); i -= 1; j -= 1
        elif code == 1:
            out_a.append(a[i - 1]); out_b.append("-"); i -= 1
        else:
            out_a.append("-"); out_b.append(b[j - 1]); j -= 1
    while i > 0:
        out_a.append(a[i - 1]); out_b.append("-"); i -= 1
    while j > 0:
        out_a.append("-"); out_b.append(b[j - 1]); j -= 1
    return "".join(reversed(out_a)), "".join(reversed(out_b))


def centre_star_msa(seqs: Sequence[str]) -> Optional[List[str]]:
    """Pure-Python centre-star multiple alignment; None if it would be too slow.

    Documented fallback for when no real aligner (mafft/famsa/muscle/clustalo)
    is resolvable - the cluster runs containers and cannot be assumed to have
    one. The centre is the longest sequence, each other sequence is aligned to
    it pairwise, and insertions are merged into the growing MSA. This is worse
    than a progressive aligner with a substitution matrix, which is exactly why
    the HMM built from it should be regarded as provisional; the manifest note
    says which aligner was used so that is visible after the fact.
    """
    if len(seqs) < 2:
        return None
    order = sorted(range(len(seqs)), key=lambda i: -len(seqs[i]))
    centre = seqs[order[0]]
    cells = sum(len(centre) * len(seqs[i]) for i in order[1:])
    if cells > MAX_NW_CELLS:
        LOG.warning("stdlib fallback alignment would need %d DP cells (> %d); "
                    "refusing. Install an aligner: %s",
                    cells, MAX_NW_CELLS, TOOL_INSTALL["mafft"])
        return None

    msa = [centre]
    for idx in order[1:]:
        ac, as_ = nw_align(msa[0].replace("-", ""), seqs[idx])
        new_rows = ["" for _ in msa]
        new_s = ""
        col = 0
        width = len(msa[0])
        for ca, cs in zip(ac, as_):
            if ca == "-":
                # Insertion relative to the centre: a brand-new column.
                for r in range(len(msa)):
                    new_rows[r] += "-"
                new_s += cs
            else:
                while col < width and msa[0][col] == "-":
                    for r in range(len(msa)):
                        new_rows[r] += msa[r][col]
                    new_s += "-"
                    col += 1
                for r in range(len(msa)):
                    new_rows[r] += msa[r][col]
                new_s += cs
                col += 1
        while col < width:
            for r in range(len(msa)):
                new_rows[r] += msa[r][col]
            new_s += "-"
            col += 1
        msa = new_rows + [new_s]

    # Restore the caller's order so row i still corresponds to seqs[i].
    out = [""] * len(seqs)
    for pos, idx in enumerate(order):
        out[idx] = msa[pos]
    return out


def align_protein(records: Sequence[SeqRecord], aligner: Optional[str],
                  aligner_name: str, workdir: Path,
                  threads: int) -> Tuple[Optional[Path], str]:
    """Write an aligned FASTA for hmmbuild. Returns (path, aligner-used)."""
    unaligned = workdir / "in.faa"
    with unaligned.open("w", encoding="utf-8") as fh:
        for i, rec in enumerate(records):
            fh.write(f">s{i}_{sanitise(rec.accession)}\n{rec.seq}\n")
    aligned = workdir / "aln.afa"

    if aligner:
        cmds = {
            "mafft": [aligner, "--quiet", "--auto", "--thread", str(threads), str(unaligned)],
            "famsa": [aligner, "-t", str(threads), str(unaligned), str(aligned)],
            "muscle": [aligner, "-align", str(unaligned), "-output", str(aligned)],
            "clustalo": [aligner, "-i", str(unaligned), "-o", str(aligned),
                         "--force", "--threads", str(threads)],
        }
        cmd = cmds.get(aligner_name)
        if cmd:
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
            except (subprocess.TimeoutExpired, OSError) as exc:
                LOG.warning("  aligner %s failed: %s", aligner_name, exc)
                proc = None
            if proc is not None and proc.returncode == 0:
                if aligner_name == "mafft":
                    aligned.write_text(proc.stdout, encoding="utf-8")
                if aligned.is_file() and aligned.stat().st_size > 0:
                    return aligned, aligner_name
            elif proc is not None:
                LOG.warning("  aligner %s exited %d: %s", aligner_name,
                            proc.returncode, (proc.stderr or "").strip().splitlines()[-1:])

    msa = centre_star_msa([r.seq for r in records])
    if msa is None:
        return None, "none"
    with aligned.open("w", encoding="utf-8") as fh:
        for i, rowseq in enumerate(msa):
            fh.write(f">s{i}_{sanitise(records[i].accession)}\n{rowseq}\n")
    return aligned, "centre-star (stdlib fallback)"


def build_hmm(hmmbuild_bin: str, marker_id: str, aligned: Path, dest: Path,
              threads: int) -> bool:
    """hmmbuild one profile, NAMEd 'marker|<marker_id>' so that hmmsearch output
    carries the marker id in the same field index as the FASTA convention."""
    tmp = dest.with_name(dest.name + ".part")
    cmd = [hmmbuild_bin, "--amino", "--informat", "afa",
           "-n", f"{MARKER_TAG}|{marker_id}", "--cpu", str(max(1, threads)),
           str(tmp), str(aligned)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except (subprocess.TimeoutExpired, OSError) as exc:
        LOG.warning("  %s: hmmbuild failed: %s", marker_id, exc)
        tmp.unlink(missing_ok=True)
        return False
    if proc.returncode != 0 or not tmp.is_file() or tmp.stat().st_size == 0:
        LOG.warning("  %s: hmmbuild exited %d: %s", marker_id, proc.returncode,
                    (proc.stderr or "").strip().splitlines()[-2:])
        tmp.unlink(missing_ok=True)
        return False
    tmp.replace(dest)
    return True


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def load_manifest(path: Path) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    if not path.is_file():
        return out
    with path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mid = (row.get("marker_id") or "").strip()
            if mid:
                out[mid] = row
    if out:
        LOG.info("resuming from existing manifest: %d row(s) in %s", len(out), path)
    return out


def write_manifest(path: Path, rows: Iterable[dict]) -> None:
    """Atomic rewrite. Called after EVERY marker, because an interrupted run
    that leaves a half-written manifest would make untested markers look tested
    or vice versa - the one failure mode this file exists to prevent."""
    tmp = path.with_name(path.name + ".part")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    tmp.replace(path)


def blank_manifest_row(row: dict, status: str, note: str) -> dict:
    """A manifest row with no sequence: explicit status, empty paths."""
    return {
        "marker_id": row["marker_id"], "taxon": row["taxon"], "gene": row["gene"],
        "marker_class": row["marker_class"], "evidence_type": row["evidence_type"],
        "query_strategy": row["query_strategy"], "status": status,
        "n_nt_seqs": "0", "n_aa_seqs": "0", "has_hmm": "no", "accessions": "",
        "nt_path": "", "aa_path": "", "hmm_path": "", "md5": "",
        "retrieved_utc": utc_now(), "note": note,
    }


def write_curation_todo(path: Path, rows: Sequence[dict]) -> None:
    """Unmet curation tasks, as a file rather than only a log line.

    curate_required and motif rows are not failures to retry - they are work
    items for a human. Writing them out means the gap is reviewable and
    diffable, instead of scrolling past in a build log.
    """
    tmp = path.with_name(path.name + ".part")
    cols = ["marker_id", "taxon", "gene", "marker_class", "evidence_type",
            "query_strategy", "status", "what_is_needed", "catalogue_notes"]
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        fh.write("# bioThreatFinder unmet marker-curation tasks\n")
        fh.write(f"# written {utc_now()} by build_marker_db.py\n")
        fh.write("# Every marker listed here evaluates as UNKNOWN downstream: it is\n")
        fh.write("# neither present nor absent, because it was never testable.\n")
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    tmp.replace(path)


# ---------------------------------------------------------------------------
# Per-marker build
# ---------------------------------------------------------------------------

class Context:
    """Everything a single marker build needs, resolved once."""

    def __init__(self, args, out_dir: Path, ez: Entrez, tools: Dict[str, Optional[str]]):
        self.args = args
        self.out_dir = out_dir
        self.nt_dir = out_dir / "nt"
        self.aa_dir = out_dir / "aa"
        self.hmm_dir = out_dir / "hmm"
        self.ez = ez
        self.tools = tools
        self.aligner_name = next((a for a in ALIGNERS if tools.get(a)), "")
        self.aligner = tools.get(self.aligner_name) if self.aligner_name else None


def process_marker(row: dict, ctx: Context, prior: Dict[str, dict]) -> dict:
    """Resolve one catalogue row to sequences and return its manifest row.

    Never raises for a retrieval problem: a failure is a `failed` row with the
    query in the note, because a crash halfway through a 163-marker build would
    leave the manifest silent about everything after it.
    """
    args = ctx.args
    marker_id = row["marker_id"]
    strategy = row["query_strategy"]
    tiers = parse_tiers(row["evidence_type"])

    if strategy not in KNOWN_STRATEGIES:
        return blank_manifest_row(
            row, STATUS_FAILED,
            f"unknown query_strategy {strategy!r}; this build does not know how to "
            f"resolve it - extend build_marker_db.py or fix the catalogue")

    # ---- strategies that are not retrievals ----------------------------
    if strategy == "curate_required":
        LOG.info("  %s: curation required (not automatable) - recorded as an unmet task",
                 marker_id)
        return blank_manifest_row(
            row, STATUS_CURATION,
            "query_strategy=curate_required: no automated retrieval is possible. "
            f"Catalogue note: {row.get('notes', '')}")

    if strategy == "motif":
        # A motif is a set of residues at defined positions (HA polybasic
        # cleavage site, Sabin attenuating positions, TC-83 substitutions), not
        # a sequence to align against. Fetching the parent gene would produce a
        # marker that matches every wild-type allele too, which is the exact
        # opposite of what these rows are for - so nothing is built and the row
        # is recorded as a TODO with the specification a human must supply.
        LOG.info("  %s: motif marker - needs a manual residue/position specification",
                 marker_id)
        return blank_manifest_row(
            row, STATUS_MOTIF,
            "query_strategy=motif: requires a manual motif specification "
            "(reference coordinates, residue positions, and wild-type vs variant "
            "states) plus a matching rule in marker_rules.yaml. A whole-gene fetch "
            f"would NOT be equivalent. Catalogue note: {row.get('notes', '')}")

    # ---- resume ---------------------------------------------------------
    nt_path = ctx.nt_dir / f"{marker_id}.fna"
    aa_path = ctx.aa_dir / f"{marker_id}.faa"
    hmm_path = ctx.hmm_dir / f"{marker_id}.hmm"

    if not args.force:
        old = prior.get(marker_id)
        if old and old.get("status") in (STATUS_BUILT, STATUS_PARTIAL):
            declared = [Path(old[k]) for k in ("nt_path", "aa_path", "hmm_path")
                        if old.get(k)]
            if declared and all(p.is_file() and p.stat().st_size > 0 for p in declared):
                got = combined_md5(declared)
                if not old.get("md5") or old["md5"] == got:
                    LOG.info("  SKIP %s: %s already built (%s), md5 OK",
                             marker_id, ", ".join(p.name for p in declared),
                             old.get("status"))
                    return dict(old)
                LOG.warning("  %s: md5 changed on disk (manifest %s, actual %s) - "
                            "rebuilding", marker_id, old["md5"], got)
            else:
                LOG.warning("  %s: manifest says %s but file(s) are missing or empty "
                            "- rebuilding", marker_id, old.get("status"))

    want = MAX_SEQS_FAMILY if strategy == "gene_family" else MAX_SEQS_SINGLE

    want_dna = "dna" in tiers
    want_aa = "protein" in tiers
    want_hmm = "hmm" in tiers

    notes: List[str] = []
    skipped_tiers: List[str] = []
    if want_aa and args.skip_protein:
        want_aa = False
        skipped_tiers.append("protein (--skip-protein)")
    if want_hmm and args.skip_hmm:
        want_hmm = False
        skipped_tiers.append("hmm (--skip-hmm)")
    if want_hmm and args.skip_protein:
        # An HMM here is a protein profile; without the protein tier there is
        # nothing to align.
        want_hmm = False
        skipped_tiers.append("hmm (needs the protein tier, which --skip-protein disabled)")

    nt_records: List[SeqRecord] = []
    aa_records: List[SeqRecord] = []
    nt_query = aa_query = ""
    accessions: List[str] = []

    if want_dna:
        try:
            nt_records, nt_query = fetch_nucleotide(ctx.ez, row, want)
        except EntrezError as exc:
            notes.append(f"nucleotide retrieval error: {exc}")
        if nt_records:
            write_marker_fasta(nt_path, marker_id, nt_records)
            LOG.info("  %s: nt %d seq(s) -> %s", marker_id, len(nt_records), nt_path.name)
        else:
            notes.append(f"no usable nucleotide sequence; last query: {nt_query or 'none built'}")

    # The protein tier is what catches a divergent or engineered allele that
    # nucleotide alignment would miss, so it is attempted even when the
    # nucleotide tier already succeeded.
    if want_aa or want_hmm:
        try:
            aa_records, aa_query = fetch_protein(ctx.ez, row, max(want, MIN_SEQS_FOR_HMM)
                                                 if want_hmm else want)
        except EntrezError as exc:
            notes.append(f"protein retrieval error: {exc}")
        if aa_records and want_aa:
            write_marker_fasta(aa_path, marker_id, aa_records)
            LOG.info("  %s: aa %d seq(s) -> %s", marker_id, len(aa_records), aa_path.name)
        elif want_aa:
            notes.append(f"no usable protein sequence; last query: {aa_query or 'none built'}")

    accessions = dedup_str([r.accession for r in nt_records] +
                           [r.accession for r in aa_records])

    # ---- HMM tier -------------------------------------------------------
    has_hmm = False
    if want_hmm:
        if not ctx.tools.get("hmmbuild"):
            notes.append("hmm tier not built: hmmbuild not found "
                         f"(install: {TOOL_INSTALL['hmmbuild']})")
        elif len(aa_records) < MIN_SEQS_FOR_HMM:
            # See MIN_SEQS_FOR_HMM: a profile from one or two sequences is not a
            # profile, and it would report confident hits with no basis.
            notes.append(f"hmm tier not built: only {len(aa_records)} protein "
                         f"sequence(s) retrieved, need >= {MIN_SEQS_FOR_HMM} for a "
                         f"meaningful profile")
        else:
            with tempfile.TemporaryDirectory(prefix=".hmm.", dir=str(ctx.out_dir)) as td:
                aligned, used = align_protein(aa_records, ctx.aligner, ctx.aligner_name,
                                              Path(td), args.threads)
                if aligned is None:
                    notes.append("hmm tier not built: no aligner available and the "
                                 "stdlib fallback declined (sequences too long)")
                elif build_hmm(ctx.tools["hmmbuild"], marker_id, aligned, hmm_path,
                               args.threads):
                    has_hmm = True
                    notes.append(f"hmm built from {len(aa_records)} sequence(s), "
                                 f"aligner: {used}")
                    LOG.info("  %s: hmm built (%d seqs, aligner %s)",
                             marker_id, len(aa_records), used)
                else:
                    notes.append("hmm tier not built: hmmbuild failed")

    # ---- status ---------------------------------------------------------
    requested = list(tiers)
    obtained = []
    if nt_records:
        obtained.append("dna")
    if aa_records and want_aa:
        obtained.append("protein")
    if has_hmm:
        obtained.append("hmm")

    if not obtained:
        status = STATUS_FAILED
        if skipped_tiers:
            notes.append("tier(s) skipped by flag: " + "; ".join(skipped_tiers))
    elif set(requested) <= set(obtained):
        status = STATUS_BUILT
    else:
        status = STATUS_PARTIAL
        missing = [t for t in requested if t not in obtained]
        notes.insert(0, "tiers built: " + "+".join(obtained) +
                     "; tiers MISSING: " + "+".join(missing))
        if skipped_tiers:
            notes.append("tier(s) skipped by flag: " + "; ".join(skipped_tiers))

    paths = [p for p in (nt_path if nt_records else None,
                         aa_path if (aa_records and want_aa) else None,
                         hmm_path if has_hmm else None) if p is not None]

    # Remove any file this run did NOT produce. Reaching this point means the
    # resume check declined to trust what was on disk, so a leftover file is
    # stale: leaving it would put sequence into markers.nt.fna/markers.aa.faa
    # that no manifest row vouches for, which is the same conflation of
    # "untested" and "tested" that this script exists to prevent.
    kept = set(paths)
    for p in (nt_path, aa_path, hmm_path):
        if p not in kept and p.is_file():
            LOG.debug("  %s: removing stale %s", marker_id, p.name)
            p.unlink(missing_ok=True)

    if status == STATUS_FAILED:
        return blank_manifest_row(row, STATUS_FAILED, " | ".join(notes) or
                                  "retrieval attempted, nothing usable returned")

    if nt_query:
        notes.append(f"nt query: {nt_query}")
    if aa_query:
        notes.append(f"aa query: {aa_query}")

    return {
        "marker_id": marker_id, "taxon": row["taxon"], "gene": row["gene"],
        "marker_class": row["marker_class"], "evidence_type": row["evidence_type"],
        "query_strategy": strategy, "status": status,
        "n_nt_seqs": str(len(nt_records)),
        "n_aa_seqs": str(len(aa_records) if want_aa else 0),
        "has_hmm": "yes" if has_hmm else "no",
        "accessions": ";".join(accessions),
        "nt_path": str(nt_path) if nt_records else "",
        "aa_path": str(aa_path) if (aa_records and want_aa) else "",
        "hmm_path": str(hmm_path) if has_hmm else "",
        "md5": combined_md5(paths),
        "retrieved_utc": utc_now(),
        "note": " | ".join(notes),
    }


# ---------------------------------------------------------------------------
# Combined outputs
# ---------------------------------------------------------------------------

def concat_files(sources: Sequence[Path], dest: Path) -> int:
    """Concatenate, atomically. Returns the number of sequences written.

    A plain concatenation is correct because the per-marker files were already
    written under the header convention - one renaming step, one place to get it
    wrong, and the per-marker file on disk is byte-identical to its slice of the
    combined file, which is what makes the manifest md5 meaningful.
    """
    n = 0
    tmp = dest.with_name(dest.name + ".part")
    with tmp.open("w", encoding="utf-8", newline="\n") as out:
        for src in sources:
            with src.open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
                    out.write(line if line.endswith("\n") else line + "\n")
    tmp.replace(dest)
    return n


def build_mm2_indexes(minimap2_bin: str, fasta: Path, out_dir: Path, stem: str,
                      threads: int, force: bool) -> List[Path]:
    """One index per preset: k-mer and window size are baked in at index time,
    so an -x sr index is unusable for ONT reads and vice versa."""
    built = []
    for preset in MM2_PRESETS:
        mmi = out_dir / f"{stem}.{preset}.mmi"
        if mmi.is_file() and mmi.stat().st_size > 0 and not force:
            LOG.info("SKIP minimap2 index (exists): %s", mmi.name)
            built.append(mmi)
            continue
        LOG.info("minimap2 -x %s index -> %s", preset, mmi.name)
        tmp = mmi.with_name(mmi.name + ".part")
        proc = subprocess.run([minimap2_bin, "-x", preset, "-t", str(threads),
                               "-d", str(tmp), str(fasta)],
                              capture_output=True, text=True)
        if proc.returncode != 0 or not tmp.is_file() or tmp.stat().st_size == 0:
            tmp.unlink(missing_ok=True)
            LOG.error("minimap2 indexing failed (%s): %s", preset,
                      (proc.stderr or "").strip().splitlines()[-3:])
            continue
        tmp.replace(mmi)
        built.append(mmi)
    return built


def build_mmseqs_db(mmseqs_bin: str, fasta: Path, db_path: Path,
                    force: bool) -> bool:
    """`mmseqs createdb` only.

    createindex is deliberately NOT run: it multiplies the on-disk footprint and
    this volume sits at ~97% full, while a marker database is small enough that
    an unindexed search is fast. Add it later if profiling says otherwise.
    """
    if db_path.is_file() and db_path.stat().st_size > 0 and not force:
        LOG.info("SKIP MMseqs2 database (exists): %s", db_path.name)
        return True
    db_path.parent.mkdir(parents=True, exist_ok=True)
    LOG.info("mmseqs createdb -> %s", db_path)
    proc = subprocess.run([mmseqs_bin, "createdb", str(fasta), str(db_path)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        LOG.error("mmseqs createdb failed: %s",
                  (proc.stderr or proc.stdout or "").strip().splitlines()[-3:])
        return False
    return True


def press_hmms(hmmpress_bin: str, combined: Path, force: bool) -> bool:
    """hmmpress writes <file>.h3{m,i,f,p} beside the input and refuses to
    overwrite them, so clear them first when rebuilding."""
    suffixes = [".h3m", ".h3i", ".h3f", ".h3p"]
    existing = [combined.with_name(combined.name + s) for s in suffixes]
    if all(p.is_file() and p.stat().st_size > 0 for p in existing) and not force:
        LOG.info("SKIP hmmpress (pressed files exist): %s.h3*", combined.name)
        return True
    for p in existing:
        p.unlink(missing_ok=True)
    LOG.info("hmmpress %s", combined.name)
    proc = subprocess.run([hmmpress_bin, str(combined)], capture_output=True, text=True)
    if proc.returncode != 0:
        LOG.error("hmmpress failed: %s",
                  (proc.stderr or proc.stdout or "").strip().splitlines()[-3:])
        return False
    return True


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def self_test() -> int:
    """Unit tests for the two load-bearing conventions, runnable offline.

    The header convention and the manifest contract are what the downstream
    three-valued logic rests on, so they are tested here rather than trusted.
    """
    failures: List[str] = []

    def check(name: str, cond: bool, detail: str = "") -> None:
        if cond:
            print(f"  ok   {name}")
        else:
            print(f"  FAIL {name}{': ' + detail if detail else ''}")
            failures.append(name)

    print("build_marker_db.py --self-test")
    print("-- header convention: rename and parse round-trip --")

    cases = [
        ("banth_pxo1_pagA", "NP_052803.1", "protective antigen [Bacillus anthracis]"),
        # Descriptions genuinely contain pipes: the legacy gi|...|ref|... style,
        # strain names, and submitter free text. The parse must keep them whole.
        ("cbot_bont", "WP_012704175.1", "botulinum neurotoxin type A|subtype A1|strain Hall"),
        ("ecoli_stx2", "AAB29202.1", "gi|123456|gb|AAB29202.1| Shiga toxin 2 subunit A"),
        ("mtb_RD9", "NC_000962.3", ""),
        ("lassa_GPC", "YP_694872.1", "  leading and trailing whitespace  "),
        ("flua_HA_H5", "ABC12345.1", "hemagglutinin [Influenza A virus (A/x/1/2004(H5N1))]"),
        ("pfal_18S", "XR_002273095.1", "18S ribosomal RNA | partial sequence"),
        ("tbg_TgsGP", "CAA_no_desc", "|"),
        ("marburg_L", "YP_001531159.1", "L protein\twith\ttabs"),
    ]
    for marker_id, acc, desc in cases:
        header = make_marker_header(marker_id, acc, desc)
        parsed = parse_marker_header(header)
        check(f"parse returns 3-tuple for {marker_id}", parsed is not None)
        if not parsed:
            continue
        got_id, got_acc, got_desc = parsed
        expect_desc = re.sub(r"[\r\n\t]+", " ", desc).strip()
        check(f"marker_id round-trips for {marker_id}", got_id == marker_id,
              f"{got_id!r} != {marker_id!r}")
        check(f"accession round-trips for {marker_id}", got_acc == acc,
              f"{got_acc!r} != {acc!r}")
        check(f"description round-trips for {marker_id} (pipes preserved)",
              got_desc == expect_desc, f"{got_desc!r} != {expect_desc!r}")
        # THE downstream contract: split on '|', take field index 1.
        naive_field_1 = header.split("|")[1]
        check(f"field index 1 is the marker_id for {marker_id}",
              naive_field_1 == marker_id, f"{naive_field_1!r} != {marker_id!r}")
        check(f"marker_id_of() agrees for {marker_id}",
              marker_id_of(header) == marker_id)
        check(f"'>' prefixed form parses for {marker_id}",
              marker_id_of(">" + header) == marker_id)

    print("-- header convention: rejections and the HMM NAME form --")
    check("non-marker header returns None", parse_marker_header(
        "NC_000962.3 Mycobacterium tuberculosis H37Rv complete genome") is None)
    check("taxid-prefixed reference (build_pack_refs style) returns None",
          parse_marker_header("taxid|1392|contig1 Bacillus anthracis") is None)
    check("empty name returns None", parse_marker_header("") is None)
    check("bare tag with empty id returns None", parse_marker_header("marker||x|y") is None)
    # 'marker|<id>' is the NAME given to each profile HMM.
    hmm_parsed = parse_marker_header(f"{MARKER_TAG}|cbot_bont")
    check("2-field HMM NAME form parses", hmm_parsed is not None and
          hmm_parsed[0] == "cbot_bont" and hmm_parsed[1] == "" and hmm_parsed[2] == "")
    try:
        make_marker_header("bad|id", "X", "y")
        check("marker_id containing '|' is rejected", False)
    except ValueError:
        check("marker_id containing '|' is rejected", True)
    try:
        make_marker_header("ok_id", "A|B", "y")
        check("accession containing '|' is rejected", False)
    except ValueError:
        check("accession containing '|' is rejected", True)

    print("-- FASTA write/read round-trip through a real file --")
    with tempfile.TemporaryDirectory(prefix="btf_selftest.") as td:
        tmpdir = Path(td)
        recs = [
            SeqRecord("NP_052803.1", "protective antigen|pXO1|Bacillus anthracis",
                      "MKKRKVLIPLMALSTILVSSTGNLEVIQAEVKQENRLLNESESSSQGLLGYYFSDLNFQAPMVVTSSTTGDLSIPS"),
            SeqRecord("WP_000728089.1", "protective antigen", "MKKRKVLIPLMALSTILVSST"),
        ]
        out = tmpdir / "banth_pxo1_pagA.faa"
        write_marker_fasta(out, "banth_pxo1_pagA", recs)
        text = out.read_text(encoding="utf-8")
        reread = parse_fasta(text)
        check("all records survive the write/read cycle", len(reread) == len(recs),
              f"{len(reread)} != {len(recs)}")
        header_lines = [ln for ln in text.splitlines() if ln.startswith(">")]
        check("every written header carries the convention",
              all(marker_id_of(ln) == "banth_pxo1_pagA" for ln in header_lines))
        # Sequences must be unchanged by wrapping.
        rejoined = []
        for ln in text.splitlines():
            if ln.startswith(">"):
                rejoined.append("")
            else:
                rejoined[-1] += ln
        check("sequences unchanged by line wrapping",
              rejoined == [r.seq for r in recs], f"{rejoined}")
        # A combined file is a plain concatenation, and attribution must survive.
        combined = tmpdir / "markers.aa.faa"
        second = tmpdir / "cbot_bont.faa"
        write_marker_fasta(second, "cbot_bont",
                           [SeqRecord("WP_012704175.1", "BoNT/A|subtype A1", "MPFVNKQFNYKDPV")])
        n = concat_files([out, second], combined)
        check("concat counts every sequence", n == 3, str(n))
        ids = [marker_id_of(ln) for ln in combined.read_text(encoding="utf-8").splitlines()
               if ln.startswith(">")]
        check("combined file attributes each sequence to its marker",
              ids == ["banth_pxo1_pagA", "banth_pxo1_pagA", "cbot_bont"], str(ids))

        print("-- manifest contract --")
        row = {"marker_id": "fnov_specific", "taxon": "Francisella novicida",
               "gene": "-", "marker_class": "absence_expected",
               "evidence_type": "dna+protein", "query_strategy": "curate_required",
               "notes": "CURATION NEEDED"}
        blank = blank_manifest_row(row, STATUS_CURATION, "because")
        check("blank row has every manifest column",
              set(blank) == set(MANIFEST_COLUMNS),
              str(set(MANIFEST_COLUMNS) ^ set(blank)))
        check("blank row has empty paths",
              blank["nt_path"] == "" and blank["aa_path"] == "" and blank["hmm_path"] == "")
        check("blank row status is not 'built'", blank["status"] != STATUS_BUILT)
        mpath = tmpdir / "markers.manifest.tsv"
        write_manifest(mpath, [blank])
        back = load_manifest(mpath)
        check("manifest round-trips by marker_id", list(back) == ["fnov_specific"])
        check("manifest preserves the status", back["fnov_specific"]["status"] == STATUS_CURATION)
        check("status taxonomy is exactly the documented five",
              set(ALL_STATUSES) == {"built", "partial", "curation_required",
                                    "failed", "motif_todo"})

        print("-- catalogue parsing --")
        cat = tmpdir / "markers.tsv"
        cat.write_text(
            "# a comment\n#\n" + "\t".join(CATALOGUE_COLUMNS) + "\n"
            "#---- section banner ----\n"
            "m1\tBacillus anthracis\tpagA\tpag;PA\trequired\tdna+protein\tpXO1\t0.8\t0.95\tsingle\tgene_in_organism\tNCBI\tnote\n"
            "m2\tClostridium botulinum\tbont\tbotA;botB\tpathotype\tdna+protein+hmm\tvariable\t0.7\t0.9\tsingle\tgene_family\tNCBI\tnote\n",
            encoding="utf-8")
        rows = read_marker_catalogue(cat)
        check("comment lines and section banners are skipped", len(rows) == 2, str(len(rows)))
        check("columns are read by name", rows[0]["gene"] == "pagA")
        check("tiers parse", parse_tiers(rows[1]["evidence_type"]) == ["dna", "protein", "hmm"])
        check("tiers ignore junk", parse_tiers("dna+DNA+nonsense+hmm") == ["dna", "hmm"])
        check("aliases parse", parse_aliases("pag;PA;protective antigen") ==
              ["pag", "PA", "protective antigen"])
        check("'-' aliases parse to nothing", parse_aliases("-") == [])

        print("-- combined md5 is sensitive to the tier set --")
        d1 = combined_md5([out])
        d2 = combined_md5([out, second])
        check("adding a file changes the digest", d1 != d2)
        check("digest is order-independent", combined_md5([second, out]) == d2)

    print("-- query construction --")
    grow = {"marker_id": "banth_pxo1_pagA", "taxon": "Bacillus anthracis",
            "gene": "pagA", "aliases": "pag;PA", "query_strategy": "gene_in_organism"}
    q = nt_queries(grow)
    check("gene_in_organism queries constrain the organism",
          all('"Bacillus anthracis"[Organism]' in x for x in q))
    check("gene_in_organism tries RefSeq first",
          "srcdb_refseq[Properties]" in q[0], q[0])
    check("aliases are tried too", any("pag[" in x for x in q))
    fam = {"marker_id": "cbot_bont", "taxon": "Clostridium botulinum", "gene": "bont",
           "aliases": "botA;botE", "query_strategy": "gene_family"}
    check("gene_family relaxes to the genus",
          any('"Clostridium"[Organism]' in x for x in nt_queries(fam)))
    reg = {"marker_id": "mtb_RD9", "taxon": "Mycobacterium tuberculosis", "gene": "RD9",
           "aliases": "region of difference 9", "query_strategy": "named_region"}
    check("named_region does not use [Gene Name]",
          all("[Gene Name]" not in x for x in nt_queries(reg)))
    check("one-character aliases are not searched (would match everything)",
          "t" not in query_names({"gene": "lpnA", "aliases": "tul4;t"}),
          str(query_names({"gene": "lpnA", "aliases": "tul4;t"})))
    check("substitution designations are not searched",
          query_names({"gene": "cyp51A", "aliases": "TR34;L98H;TR46;Y121F;T289A"})
          == ["cyp51A"],
          str(query_names({"gene": "cyp51A", "aliases": "TR34;L98H;TR46;Y121F;T289A"})))
    check("real subtype aliases survive the filter",
          query_names({"gene": "HA", "aliases": "H5"}) == ["HA", "H5"])
    check("rows with gene '-' yield no searchable name",
          query_names({"gene": "-", "aliases": "-"}) == [])
    check("gene name matching is case- and prefix-tolerant",
          gene_name_matches("capBCADE", ["capB"]) and gene_name_matches("STX2A", ["stx2"])
          and not gene_name_matches("recA", ["capB"]))

    print("-- CDS-flavoured FASTA parsing --")
    cds = (">lcl|NC_000962.3_cds_NP_214515.1_1 [gene=dnaA] [protein=chromosomal "
           "replication initiator protein DnaA] [protein_id=NP_214515.1] "
           "[location=1..1524]\nATGGCGGAT\n"
           ">lcl|NC_000962.3_cds_NP_214516.1_2 [locus_tag=Rv0002] [protein_id=NP_214516.1]\n"
           "ATGTAA\n")
    parsed_cds = parse_fasta(cds)
    check("CDS accession comes from [protein_id=]",
          [r.accession for r in parsed_cds] == ["NP_214515.1", "NP_214516.1"],
          str([r.accession for r in parsed_cds]))
    check("CDS gene comes from [gene=]/[locus_tag=]",
          [r.gene for r in parsed_cds] == ["dnaA", "Rv0002"])
    plain = parse_fasta(">NC_000962.3 Mycobacterium tuberculosis H37Rv\nACGT\nACGT\n")
    check("plain FASTA accession is the first token",
          len(plain) == 1 and plain[0].accession == "NC_000962.3" and plain[0].seq == "ACGTACGT")
    check("dedup drops identical sequences",
          len(dedup_records([SeqRecord("a", "", "ACGT"), SeqRecord("b", "", "ACGT"),
                             SeqRecord("c", "", "ACGA")])) == 2)

    print("-- stdlib fallback aligner --")
    msa = centre_star_msa(["MKKRKVLIPLMAL", "MKKRVLIPLMAL", "MKKRKVLIPLMALST"])
    check("centre-star returns one row per input", msa is not None and len(msa) == 3)
    if msa:
        check("all alignment rows are the same length",
              len({len(r) for r in msa}) == 1, str([len(r) for r in msa]))
        check("ungapped rows equal the inputs",
              [r.replace("-", "") for r in msa] ==
              ["MKKRKVLIPLMAL", "MKKRVLIPLMAL", "MKKRKVLIPLMALST"])
    a, b = nw_align("ACGT", "ACGT")
    check("identical sequences align without gaps", a == "ACGT" and b == "ACGT")

    print("-- Entrez cache key --")
    with tempfile.TemporaryDirectory(prefix="btf_selftest_cache.") as td:
        ez_a = Entrez(Path(td), email="a@example.org", api_key="KEY1", offline=True)
        ez_b = Entrez(Path(td), email="b@example.org", api_key="KEY2", offline=True)
        params = {"db": "protein", "term": "pagA[Gene Name]", "retmax": "40"}
        check("cache key ignores credentials",
              ez_a._cache_path("esearch.fcgi", params) ==
              ez_b._cache_path("esearch.fcgi", params))
        check("cache key separates different queries",
              ez_a._cache_path("esearch.fcgi", params) !=
              ez_a._cache_path("esearch.fcgi", dict(params, term="lef[Gene Name]")))
        check("no-key client throttles to 3/s",
              abs(Entrez(Path(td), offline=True).min_interval - 1 / 3.0) < 1e-9)
        check("api-key client throttles to 10/s",
              abs(ez_a.min_interval - 1 / 10.0) < 1e-9)
        try:
            ez_a.esearch("protein", "anything")
            offline_ok = True  # esearch swallows EntrezError and returns []
        except Exception:
            offline_ok = False
        check("offline client never touches the network", offline_ok)

    print("")
    if failures:
        print(f"SELF-TEST FAILED: {len(failures)} check(s) failed:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("SELF-TEST PASSED")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(
        prog="build_marker_db.py",
        description="Resolve a bioThreatFinder pack's marker catalogue to real "
                    "sequences and build the nucleotide, protein and profile-HMM "
                    "marker databases, recording for every catalogue row whether "
                    "usable sequence was actually obtained.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  build_marker_db.py --self-test
  build_marker_db.py --pack human --dry-run
  build_marker_db.py --pack human --email you@example.org --api-key $NCBI_API_KEY
  build_marker_db.py --pack human --only-marker banth_pxo1_pagA -v
  build_marker_db.py --pack human --only-taxon 'Yersinia pestis' --skip-hmm

manifest status values (markers.manifest.tsv):
  built              usable sequence for every evidence tier the catalogue asked for
  partial            some tiers only; the note says which are missing and why
  curation_required  query_strategy=curate_required - not automatable, a human task
  failed             retrieval was attempted and returned nothing usable; the note
                     records the query that was used
  motif_todo         a motif row needing a manual residue/position specification

  EVERY catalogue row gets a manifest row. A marker with no usable sequence has a
  non-'built' status and empty paths, and evaluates as UNKNOWN downstream - which
  is NOT the same verdict as 'absent from the sample'.

notes:
  * Idempotent and resumable. An already-built marker whose files exist and whose
    md5 matches the manifest is SKIPped; the manifest is rewritten atomically
    after every marker, so interrupt and rerun freely. --force overrides.
  * Retrieval goes through NCBI E-utilities and is cached under
    <db-root>/bioThreatFinder/cache/entrez/ (esearch TTL 30 d, efetch TTL 365 d).
  * Sequences per marker are fixed constants, not flags: %d for a single-locus
    marker, %d for a gene_family marker (divergence is the point of a family).
  * Exit status is non-zero only when a 'required'-class marker for a
    high_consequence taxon failed; everything else is reported and warned about.
""" % (MAX_SEQS_SINGLE, MAX_SEQS_FAMILY))
    p.add_argument("--pack", help="pack name, e.g. human (must exist under --packs-dir)")
    p.add_argument("--packs-dir", type=Path, default=repo_root / "packs",
                   help="directory containing pack definitions (default: %(default)s)")
    p.add_argument("--db-root", type=Path, default=DEFAULT_DB_ROOT,
                   help="shared database root (default: %(default)s)")
    p.add_argument("--threads", type=int, default=8,
                   help="threads for minimap2/hmmbuild/aligners (default: %(default)s)")
    p.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""),
                   help="contact address sent to NCBI E-utilities as required by "
                        "their usage policy (default: $NCBI_EMAIL)")
    p.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""),
                   help="NCBI API key; raises the rate limit from 3 to 10 req/s "
                        "(default: $NCBI_API_KEY)")
    p.add_argument("--only-marker", action="append", default=[], metavar="ID",
                   help="restrict to these marker_ids (repeatable); for debugging. "
                        "NOTE: the manifest then describes only these markers.")
    p.add_argument("--only-taxon", action="append", default=[], metavar="NAME",
                   help="restrict to catalogue rows whose taxon matches (repeatable, "
                        "case-insensitive, substring). Same manifest caveat.")
    p.add_argument("--dry-run", action="store_true",
                   help="print the plan, the resolved tools and the queries that "
                        "would be issued, then exit without any network access")
    p.add_argument("--force", action="store_true",
                   help="rebuild everything and bypass the Entrez response cache")
    p.add_argument("--min-free-gib", type=int, default=DEFAULT_MIN_FREE_GIB,
                   help="abort if less than this much free space (default: %(default)s)")
    p.add_argument("--skip-hmm", action="store_true",
                   help="do not build profile HMMs (markers needing them are recorded "
                        "as 'partial', never as 'built')")
    p.add_argument("--skip-protein", action="store_true",
                   help="do not build the protein tier (implies --skip-hmm, since the "
                        "profiles here are protein profiles)")
    p.add_argument("--self-test", action="store_true",
                   help="run the built-in unit tests (offline) and exit")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    return p.parse_args(argv)


def resolve_tools(args) -> Dict[str, Optional[str]]:
    """Resolve every external tool and print the exact install command for the
    missing ones. Nothing is assumed to be on PATH: this runs in containers.
    """
    needed = ["hmmbuild", "hmmpress", "mmseqs", "minimap2", "seqkit", "datasets"] + ALIGNERS
    tools = {name: which(name) for name in needed}

    LOG.info("--- external tools ---")
    for name in needed:
        path = tools[name]
        LOG.info("  %-9s %s", name, path or "NOT FOUND")

    def complain(name: str, why: str, fatal: bool) -> None:
        msg = (f"{name} not found - {why}\n"
               f"       install: {TOOL_INSTALL.get(name, 'see the tool documentation')}")
        if fatal:
            raise SystemExit(f"ERROR: {msg}")
        LOG.warning(msg)

    if not args.skip_hmm and not args.skip_protein:
        if not tools["hmmbuild"]:
            complain("hmmbuild", "the hmm tier cannot be built; markers requesting it "
                                 "will be recorded 'partial'. Pass --skip-hmm to make "
                                 "that explicit.", fatal=False)
        if not tools["hmmpress"]:
            complain("hmmpress", "markers.hmm cannot be pressed, so hmmscan will not "
                                 "use it", fatal=False)
    if not args.skip_protein and not tools["mmseqs"]:
        complain("mmseqs", "the MMseqs2 marker database cannot be built; the translated "
                           "tier is what catches a divergent or engineered allele, so "
                           "this is a real loss of sensitivity", fatal=False)
    if not tools["minimap2"]:
        complain("minimap2", "markers.nt.*.mmi cannot be built; the nucleotide FASTA is "
                             "still written and can be indexed later", fatal=False)
    if not any(tools[a] for a in ALIGNERS):
        LOG.warning("no MSA tool found (%s) - profile HMMs will use the pure-Python "
                    "centre-star fallback, which is adequate but not equivalent. "
                    "install: %s", ", ".join(ALIGNERS), TOOL_INSTALL["mafft"])
    return tools


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%SZ",
        stream=sys.stderr,
    )
    logging.Formatter.converter = time.gmtime

    if args.self_test:
        return self_test()
    if not args.pack:
        raise SystemExit("ERROR: --pack is required (or use --self-test)")

    pack_dir = args.packs_dir / args.pack
    catalogue_path = pack_dir / "markers.tsv"
    out_dir = args.db_root / "bioThreatFinder" / "packs" / args.pack / "markers"
    manifest_path = out_dir / "markers.manifest.tsv"
    todo_path = out_dir / "curation_required.tsv"
    cache_dir = args.db_root / "bioThreatFinder" / "cache" / "entrez"

    LOG.info("pack=%s db_root=%s threads=%d", args.pack, args.db_root, args.threads)
    LOG.info("catalogue: %s", catalogue_path)
    LOG.info("output:    %s", out_dir)
    LOG.info("manifest:  %s", manifest_path)
    LOG.info("cache:     %s (esearch TTL %d d, efetch TTL %d d)",
             cache_dir, ENTREZ_CACHE_TTL_SEARCH_DAYS, ENTREZ_CACHE_TTL_FETCH_DAYS)

    if not args.email:
        LOG.warning("no --email given. NCBI's usage policy asks for a contact address "
                    "and they will throttle or block anonymous bulk use; set --email "
                    "or $NCBI_EMAIL.")

    rows = read_marker_catalogue(catalogue_path)
    LOG.info("catalogue: %d marker row(s)", len(rows))

    selected = rows
    if args.only_marker:
        keep = set(args.only_marker)
        selected = [r for r in selected if r["marker_id"] in keep]
        unknown = keep - {r["marker_id"] for r in rows}
        if unknown:
            raise SystemExit(f"ERROR: --only-marker matched no catalogue row: "
                             f"{', '.join(sorted(unknown))}")
    if args.only_taxon:
        wanted = [t.strip().lower() for t in args.only_taxon]
        selected = [r for r in selected
                    if any(w == r["taxon"].lower() or w in r["taxon"].lower()
                           for w in wanted)]
        if not selected:
            raise SystemExit("ERROR: --only-taxon matched no catalogue rows")
    if len(selected) != len(rows):
        LOG.warning("restricted to %d of %d catalogue row(s): the manifest written by "
                    "this run describes ONLY these markers, so do not treat it as the "
                    "pack's full tested/untested record", len(selected), len(rows))

    # ---- plan -----------------------------------------------------------
    by_strategy: Dict[str, int] = {}
    by_tier: Dict[str, int] = {}
    for r in selected:
        by_strategy[r["query_strategy"]] = by_strategy.get(r["query_strategy"], 0) + 1
        for t in parse_tiers(r["evidence_type"]):
            by_tier[t] = by_tier.get(t, 0) + 1
    LOG.info("plan: query_strategy breakdown: %s",
             ", ".join(f"{k}={v}" for k, v in sorted(by_strategy.items())))
    LOG.info("plan: evidence tiers requested: %s",
             ", ".join(f"{k}={v}" for k, v in sorted(by_tier.items())))
    unknown_strategies = set(by_strategy) - KNOWN_STRATEGIES
    if unknown_strategies:
        LOG.warning("catalogue contains query_strategy value(s) this build cannot "
                    "handle: %s - they will be recorded 'failed'",
                    ", ".join(sorted(unknown_strategies)))

    tools = resolve_tools(args)
    categories = read_taxon_categories(pack_dir)

    if args.dry_run:
        LOG.info("[dry-run] would create %s (nt/, aa/, hmm/, mmseqs/)", out_dir)
        for r in selected:
            strategy = r["query_strategy"]
            if strategy == "curate_required":
                LOG.info("[dry-run] %-22s %-14s -> %s (no retrieval)",
                         r["marker_id"], strategy, STATUS_CURATION)
                continue
            if strategy == "motif":
                LOG.info("[dry-run] %-22s %-14s -> %s (no retrieval)",
                         r["marker_id"], strategy, STATUS_MOTIF)
                continue
            tiers = parse_tiers(r["evidence_type"])
            LOG.info("[dry-run] %-22s %-14s tiers=%s", r["marker_id"], strategy,
                     "+".join(tiers))
            nq, aq = nt_queries(r), aa_queries(r)
            if "dna" in tiers and nq:
                LOG.info("[dry-run]     nuccore: %s", nq[0])
            if ("protein" in tiers or "hmm" in tiers) and aq:
                LOG.info("[dry-run]     protein: %s", aq[0])
            if args.verbose:
                for q in nq[1:]:
                    LOG.debug("[dry-run]     nuccore fallback: %s", q)
                for q in aq[1:]:
                    LOG.debug("[dry-run]     protein fallback: %s", q)
        LOG.info("[dry-run] combined outputs: markers.nt.fna (+ %s), markers.aa.faa "
                 "(+ MMseqs2 db), markers.hmm (+ hmmpress)",
                 ", ".join(f"markers.nt.{p}.mmi" for p in MM2_PRESETS))
        LOG.info("[dry-run] manifest would be written to %s", manifest_path)
        LOG.info("[dry-run] curation tasks would be written to %s", todo_path)
        LOG.info("[dry-run] nothing fetched, nothing written")
        return 0

    check_free_space(args.db_root, args.min_free_gib)
    for d in (out_dir, out_dir / "nt", out_dir / "aa", out_dir / "hmm"):
        d.mkdir(parents=True, exist_ok=True)

    ez = Entrez(cache_dir, email=args.email, api_key=args.api_key, force=args.force)
    ctx = Context(args, out_dir, ez, tools)
    if ctx.aligner_name:
        LOG.info("MSA tool for hmmbuild: %s (%s)", ctx.aligner_name, ctx.aligner)

    prior = load_manifest(manifest_path)
    manifest_rows: List[dict] = []

    # ---- per-marker loop ------------------------------------------------
    for i, row in enumerate(selected, start=1):
        LOG.info("[%d/%d] %s (%s / %s, %s, %s)", i, len(selected), row["marker_id"],
                 row["taxon"], row["gene"] or "-", row["marker_class"],
                 row["query_strategy"])
        try:
            result = process_marker(row, ctx, prior)
        except KeyboardInterrupt:
            raise
        except Exception as exc:  # noqa: BLE001 - see the docstring: never crash
            # An unexpected error must still leave a truthful manifest row, or
            # the marker silently becomes "absent" instead of "untested".
            LOG.exception("%s: unexpected error during build", row["marker_id"])
            result = blank_manifest_row(row, STATUS_FAILED,
                                        f"unexpected error: {type(exc).__name__}: {exc}")
        manifest_rows.append(result)
        # Rewrite after each marker so an interrupted run resumes cleanly.
        write_manifest(manifest_path, manifest_rows)

    LOG.info("entrez: %d request(s), %d cache hit(s)", ez.n_requests, ez.n_cache_hits)

    # ---- combined outputs ----------------------------------------------
    nt_sources = [Path(r["nt_path"]) for r in manifest_rows
                  if r.get("nt_path") and Path(r["nt_path"]).is_file()]
    aa_sources = [Path(r["aa_path"]) for r in manifest_rows
                  if r.get("aa_path") and Path(r["aa_path"]).is_file()]
    hmm_sources = [Path(r["hmm_path"]) for r in manifest_rows
                   if r.get("hmm_path") and Path(r["hmm_path"]).is_file()]

    combined_nt = out_dir / "markers.nt.fna"
    combined_aa = out_dir / "markers.aa.faa"
    combined_hmm = out_dir / "markers.hmm"

    if nt_sources:
        n = concat_files(nt_sources, combined_nt)
        LOG.info("wrote %s (%d sequence(s) from %d marker(s))",
                 combined_nt.name, n, len(nt_sources))
        if tools.get("minimap2"):
            build_mm2_indexes(tools["minimap2"], combined_nt, out_dir, "markers.nt",
                              args.threads, args.force)
        else:
            LOG.warning("minimap2 absent: %s left unindexed (install: %s)",
                        combined_nt.name, TOOL_INSTALL["minimap2"])
    else:
        LOG.warning("no nucleotide sequences at all - %s not written", combined_nt.name)

    if aa_sources:
        n = concat_files(aa_sources, combined_aa)
        LOG.info("wrote %s (%d sequence(s) from %d marker(s))",
                 combined_aa.name, n, len(aa_sources))
        if tools.get("mmseqs"):
            build_mmseqs_db(tools["mmseqs"], combined_aa, out_dir / "mmseqs" / "markers_aa",
                            args.force)
        else:
            LOG.warning("mmseqs absent: no MMseqs2 database built (install: %s)",
                        TOOL_INSTALL["mmseqs"])
    elif not args.skip_protein:
        LOG.warning("no protein sequences at all - %s not written", combined_aa.name)

    if hmm_sources:
        tmp = combined_hmm.with_name(combined_hmm.name + ".part")
        with tmp.open("w", encoding="utf-8", newline="\n") as out:
            for src in hmm_sources:
                out.write(src.read_text(encoding="utf-8", errors="replace"))
        tmp.replace(combined_hmm)
        LOG.info("wrote %s (%d profile(s))", combined_hmm.name, len(hmm_sources))
        if tools.get("hmmpress"):
            press_hmms(tools["hmmpress"], combined_hmm, args.force)
        else:
            LOG.warning("hmmpress absent: %s not pressed (install: %s)",
                        combined_hmm.name, TOOL_INSTALL["hmmpress"])
    elif not args.skip_hmm:
        LOG.warning("no profile HMMs built - %s not written", combined_hmm.name)

    write_manifest(manifest_path, manifest_rows)

    # ---- curation TODO file ---------------------------------------------
    todo_rows = []
    for r in manifest_rows:
        if r["status"] not in (STATUS_CURATION, STATUS_MOTIF):
            continue
        src = next((x for x in selected if x["marker_id"] == r["marker_id"]), {})
        todo_rows.append({**r, "what_is_needed": r.get("note", ""),
                          "catalogue_notes": src.get("notes", "")})
    write_curation_todo(todo_path, todo_rows)

    # ---- summary --------------------------------------------------------
    counts = {s: 0 for s in ALL_STATUSES}
    for r in manifest_rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1

    by_status: Dict[str, List[dict]] = {s: [] for s in ALL_STATUSES}
    for r in manifest_rows:
        by_status.setdefault(r["status"], []).append(r)

    print("", file=sys.stderr)
    LOG.info("=========================== SUMMARY ===========================")
    LOG.info("pack:              %s", args.pack)
    LOG.info("markers dir:       %s", out_dir)
    LOG.info("manifest:          %s (%d row(s))", manifest_path, len(manifest_rows))
    LOG.info("catalogue rows:    %d selected of %d in %s",
             len(selected), len(rows), catalogue_path.name)
    for status in ALL_STATUSES:
        LOG.info("  %-18s %d", status + ":", counts.get(status, 0))
    LOG.info("sequences:         %d nt marker file(s), %d aa marker file(s), %d profile(s)",
             len(nt_sources), len(aa_sources), len(hmm_sources))

    for status, banner in (
        (STATUS_CURATION, "--- curation_required: NOT automatable, a human must curate these ---"),
        (STATUS_MOTIF, "--- motif_todo: need a manual residue/position specification ---"),
        (STATUS_FAILED, "--- failed: retrieval was attempted and returned nothing usable ---"),
    ):
        if by_status.get(status):
            LOG.warning(banner)
            for r in by_status[status]:
                LOG.warning("  %-24s %-20s %-16s %s", r["marker_id"],
                            r["taxon"][:20], r["marker_class"],
                            (r["note"] or "")[:160])
    if by_status.get(STATUS_PARTIAL):
        LOG.warning("--- partial: some evidence tiers only (the missing tier is UNKNOWN, "
                    "not negative) ---")
        for r in by_status[STATUS_PARTIAL]:
            LOG.warning("  %-24s %-20s %s", r["marker_id"], r["taxon"][:20],
                        (r["note"] or "")[:160])

    LOG.info("---------------------------------------------------------------")
    LOG.warning("REMINDER: every marker above that is not 'built' evaluates as UNKNOWN "
                "downstream - it is neither present nor absent. UNKNOWN markers cannot "
                "satisfy a required-marker rule, so any finding that depends on one is "
                "held at 'supported' and is NOT promoted to 'confirmed_in_silico'. "
                "Fix or curate them before operational use, and rerun: completed "
                "markers are skipped.")
    LOG.info("curation tasks:    %s (%d row(s))", todo_path, len(todo_rows))
    try:
        LOG.info("disk used here:    %s",
                 human(sum(f.stat().st_size for f in out_dir.rglob("*") if f.is_file())))
    except OSError:
        pass
    LOG.info("free space:         %s", human(shutil.disk_usage(args.db_root).free))
    LOG.info("===============================================================")

    # ---- exit code ------------------------------------------------------
    # Non-zero ONLY for the case a human must act on before the pack is usable:
    # a required-class marker for a high_consequence taxon that could not be
    # retrieved. Those taxa cannot reach confirmed_in_silico without it, and a
    # CI or Nextflow wrapper should fail loudly rather than ship a pack that
    # silently cannot confirm anthrax.
    blocking = [r for r in manifest_rows
                if r["status"] == STATUS_FAILED
                and r["marker_class"] == "required"
                and categories.get(r["taxon"].lower()) == "high_consequence"]
    if blocking:
        LOG.error("--- BLOCKING: required marker(s) for high_consequence taxa failed ---")
        for r in blocking:
            LOG.error("  %s (%s)", r["marker_id"], r["taxon"])
        LOG.error("exiting non-zero: these taxa cannot be confirmed at all.")
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        LOG.warning("interrupted - rerun to resume (built markers are skipped and the "
                    "manifest was rewritten after the last completed marker)")
        sys.exit(130)
