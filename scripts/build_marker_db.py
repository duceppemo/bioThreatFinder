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

EXISTENCE IS NOT FRESHNESS (--check)
------------------------------------
D45 established the predicate for the competitive bundles and D46 carried it to
the sylph sketch: a database is current when it covers its INPUTS, not when it
opens. This database's input is markers.tsv, and nothing checked that until
--check existed.

The consequence of a stale marker database is quiet and points the wrong way.
marker_evidence.py treats markers.manifest.tsv as AUTHORITATIVE for tested vs
untested, so a marker added to the catalogue after the last build has no
manifest row and evaluates UNKNOWN. UNKNOWN never promotes a finding, so the
cost is a lost confirmation - and it is indistinguishable from the 16 markers
that are legitimately unevaluable on this pack (3 failed, 4 curation_required,
9 motif_todo). An operator reading "UNKNOWN" attributes it to the known curation
gap, which is the one explanation that requires no action.

Four comparisons, because each one sees what the others cannot:

  catalogue vs manifest   EXACT, and the predicate proper. Catches a marker
                          added or removed with every mtime settled.
  mtime                   catches an EDITED marker, which set coverage cannot
                          see: change an accession or a gene under a marker_id
                          the manifest already knows and every set still agrees
                          while the sequence underneath is wrong. That is D41's
                          wrong-gene bug exactly, and it is the inverse of the
                          sylph case, where the sets caught what mtime missed.
  manifest vs artifact    the manifest claims a tier the combined file does not
                          hold. This is the one that matters most, because it
                          fails UNSAFE: marker_evidence.py believes the claim,
                          searches, finds nothing, and reports ABSENT rather
                          than UNKNOWN. ABSENT is a positive assertion that the
                          marker was looked for and is not there.
  index vs source         an index older than the FASTA or HMM it was built
                          from means the aligner searches the previous
                          sequences while every set above agrees. D45's
                          index check, same reasoning.

No sidecar record is needed here, which is the difference from D46. That record
exists because a run cannot open a .syldb - sylph lives in its own env and is
used inside the process, so the head node can only compare timestamps. These
artifacts are plain text: the combined FASTAs carry the marker_id at field 1 of
every header and each profile carries it in its NAME line, so the run reads the
genuine artifact directly and there is nothing to keep stamped. Presence is
compared, never sequence COUNTS - the combined files are de-duplicated within
each marker while the manifest counts the per-marker archive, so equal counts
are not expected and demanding them would report every de-duplicated marker as
stale.

  scripts/build_marker_db.py --pack human --check    # read-only
    exit 0  fresh      exit 3  stale, rebuild      exit 1  could not tell

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
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

LOG = logging.getLogger("build_marker_db")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# PINNED taxonomy snapshot. The date is deliberate - the watchlist resolves against this
# exact snapshot and NCBI overwrites new_taxdump.tar.gz in place, so it is not re-fetchable.
PINNED_TAXDUMP = "2026-09-29"

# WHERE THE DATABASES LIVE IS NOT THIS SCRIPT'S TO ASSUME.
#
# These two used to be fixed paths - two volumes on
# the machine this was written on, which nobody else has. nextflow.config stopped defaulting
# to them (it refuses to start with either root unset); the build scripts had not, so a
# clone's build still pointed at paths that do not exist. $BTF_DB is the same variable the
# pipeline reads, and --db-root is its PARENT because this script appends bioThreatFinder/
# itself. None means "not determined": main() refuses rather than guessing.
_ENV_BTF_DB = os.environ.get("BTF_DB")
DEFAULT_BTF_DB = Path(_ENV_BTF_DB) if _ENV_BTF_DB else None
DEFAULT_DB_ROOT = DEFAULT_BTF_DB.parent if DEFAULT_BTF_DB else None
DEFAULT_TAXDUMP = (DEFAULT_BTF_DB / "taxdump" / PINNED_TAXDUMP
                   if DEFAULT_BTF_DB else None)
# LAST-RESORT place to look for a tool, after $PATH. Overridable, because this used to be
# one developer's miniconda prefix hardcoded. $CONDA_PREFIX is set inside any activated
# env, so the common case needs no configuration.
CONDA_BIN = Path(os.environ.get("BTF_CONDA_BIN")
                 or (Path(os.environ["CONDA_PREFIX"]) / "bin"
                     if os.environ.get("CONDA_PREFIX") else "/opt/conda/bin"))

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
    # HOW MANY ENTREZ REQUESTS FAILED WHILE THIS ROW WAS BUILT.
    #
    # A `built` row is trusted and SKIPPED by the resume path, so a row built
    # while the service was refusing requests is pinned in with no way to tell
    # it apart from a clean one. The 2026-10-15 --force rebuild lost 14 requests
    # across 6 markers and all 6 recovered from other rungs, which was only
    # knowable by correlating the log against the manifest by hand. Without this
    # column "re-run the ones that were affected" is not an executable
    # instruction - which is the whole recovery plan STATUS_RETRIEVAL_ERROR
    # assumes.
    "retrieval_errors",
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
# A REQUEST THAT FAILED IS NOT A SEARCH THAT FOUND NOTHING.
#
# Entrez.esearch()/efetch() catch EntrezError and return []/"" at LOG.debug, so the
# `except EntrezError` handlers in process_marker() are unreachable and a rate-limit
# block, a 5xx exhaustion or a truncated efetch was written as `failed` - documented
# as "retrieval was attempted and returned nothing usable", i.e. evidence that NCBI
# has nothing. The note then asserted "N queries tried, all empty", a positive claim
# about a search that did not complete. This host was HTTP 429'd throughout the
# 2026-10-15 rebuild, so that was not hypothetical.
STATUS_RETRIEVAL_ERROR = "retrieval_error"
ALL_STATUSES = [STATUS_BUILT, STATUS_PARTIAL, STATUS_CURATION,
                STATUS_FAILED, STATUS_MOTIF, STATUS_RETRIEVAL_ERROR]

KNOWN_STRATEGIES = {"gene_in_organism", "gene_family", "named_region", "motif",
                    "curate_required", "subtype_in_organism"}
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


def md5_or_none(path: Path) -> Optional[str]:
    """md5 of a file, or None when it does not exist. Used to decide whether a
    downstream index is stale.

    WHY. concat_files() always rewrites markers.nt.fna / markers.aa.faa / markers.hmm
    from every manifest row, but build_mm2_indexes(), build_mmseqs_db() and
    press_hmms() all SKIP when their output exists and --force was not passed. So
    `--only-marker bmand_18S` put the marker's sequences into markers.nt.fna, logged
    "SKIP minimap2 index (exists)", and left an index that does not contain them. The
    manifest then said status=built, marker_evidence.py treated the marker as testable,
    and the nucleotide tier could never match it - producing state=ABSENT, which is
    positive evidence AGAINST a real finding, rather than UNKNOWN. bmand_18S is the
    `required` marker for Balamuthia mandrillaris.
    """
    try:
        return md5_of(path) if path.is_file() else None
    except OSError:
        return None


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
    # A vendored minimap2 was accepted here from a hardcoded path on the machine this was
    # written on. Removed rather than generalised: a copy that "may predate the validated
    # release" is a worse outcome than a clear failure. BTF_<TOOL> points at one explicitly.
    override = os.environ.get(f"BTF_{name.upper()}")
    if override:
        p = Path(override)
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
        LOG.warning("BTF_%s is set to %s, which is not an executable file; ignoring",
                    name.upper(), override)
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


def sam_name_of(header: str) -> str:
    """The part of a FASTA header an aligner uses as the reference NAME.

    Everything from the first whitespace onwards is a description and is discarded by
    minimap2/samtools. That distinction is load-bearing and it was missed: see
    write_marker_fasta().
    """
    return header.split(None, 1)[0]


def write_marker_fasta(path: Path, marker_id: str, records: Sequence[SeqRecord],
                       wrap: int = 60) -> None:
    """Write records under the load-bearing header convention, atomically.

    REFERENCE NAMES MUST BE UNIQUE, AND THE DESCRIPTION DOES NOT HELP.
    samtools rejects a SAM header carrying two @SQ lines with the same SN, and an
    aligner's reference name is only the text before the first whitespace. The header
    is `marker|<id>|<accession>|<description>`, so the name ends at the first space
    INSIDE the description - in practice `marker|banth_pxo1_cya|WP_000197748.1|[gene=cya]`.
    For a CDS-level fetch the same PROTEIN accession recurs for the same gene across
    different strains (protein_id WP_000197748.1 appeared for both CWH_RS28280 and
    DBA56_RS28000), and the only distinguishing text - locus_tag, location - sits after
    a space where no aligner looks.
    Result: 16 of 524 sequences collided, minimap2 emitted the duplicate @SQ lines,
    `samtools sort` failed with "Duplicate entry ... in sam header", and MARKER_ALIGN_NT
    died - taking the entire nucleotide marker tier with it. Found by running the
    pipeline; no amount of reading the FASTA found it, because the full headers ARE
    distinct and only the first tokens collide.
    """
    tmp = path.with_name(path.name + ".part")
    seen: dict[str, int] = {}
    with tmp.open("w", encoding="utf-8", newline="\n") as fh:
        for rec in records:
            # SANITISE FIRST, THEN DISAMBIGUATE.
            #
            # The first version of this appended the suffix to the RAW accession, and
            # sanitise_accession() then matched `pdb|7CKL|A.dup2` against its
            # `^pdb\|([^|]+)\|` rule and returned `7CKL` - silently discarding the
            # suffix, so the collision survived and concat_files' guard rejected the
            # file again. Two PDB chains of one entry (7CKL chain A and B for lassa_Z,
            # 4CO6 for nipah_N) both reduce to the same accession, and both
            # descriptions begin "Chain", so the reference name matched exactly.
            # sanitise_accession is idempotent on an already-clean accession, so
            # suffixing afterwards is safe.
            acc = sanitise_accession(rec.accession)
            header = make_marker_header(marker_id, acc, rec.description)
            name = sam_name_of(header)
            if name in seen:
                seen[name] += 1
                # Disambiguate in the ACCESSION field, which keeps marker_id at
                # index 1 for parse_marker_header() and marker_evidence.py, and keeps
                # the original accession readable as the prefix.
                header = make_marker_header(
                    marker_id, f"{acc}.dup{seen[name]}", rec.description)
                LOG.info("  %s: duplicate reference name %s -> %s",
                         marker_id, name, sam_name_of(header))
            else:
                seen[name] = 1
            fh.write(">" + header + "\n")
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


def missing_tiers(manifest_row: dict) -> List[str]:
    """Tiers the catalogue declares that this manifest row does not actually have.

    MIRRORS bin/marker_evidence.py's rule exactly, and must keep doing so. That
    module decides what `tier_used` may claim from the same three columns; if the
    two drift, the build's to-do list and the run's reported gaps disagree about
    which markers are instruments and which are claims.

    A `partial` status says SOME tier is missing. These three columns say which,
    and that is the part a curator needs.
    """
    have = set()
    try:
        if int(float(manifest_row.get("n_nt_seqs") or 0)) > 0:
            have.add("dna")
    except ValueError:
        pass
    try:
        if int(float(manifest_row.get("n_aa_seqs") or 0)) > 0:
            have.add("protein")
    except ValueError:
        pass
    if (manifest_row.get("has_hmm") or "").strip().lower() in ("yes", "true", "1"):
        have.add("hmm")
    return [t for t in parse_tiers(manifest_row.get("evidence_type", ""))
            if t not in have]


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
        # Requests that FAILED, as distinct from requests that returned nothing.
        self.n_errors = 0
        self.last_error = ""
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
            self.n_errors += 1
            self.last_error = f"esearch {term!r}: {exc}"
            LOG.warning("esearch FAILED (%s): %s - this is not an empty result",
                        term, exc)
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
            self.n_errors += 1
            self.last_error = f"efetch {db}/{rettype}: {exc}"
            LOG.warning("efetch FAILED (%s %s): %s - this is not an empty result",
                        db, rettype, exc)
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
    gene = (row.get("gene") or "").strip()
    out: List[str] = []

    # THE SUBSTITUTION FILTER APPLIES TO ALIASES ONLY, NEVER TO THE DECLARED GENE.
    #
    # SUBSTITUTION_RE exists for alias lists like cyp51A's `TR34;L98H;TR46;Y121F;
    # T289A` - resistance substitutions that are not searchable gene symbols. Every
    # value it was written for is an ALIAS. Applied to the gene column it ate real
    # symbols that happen to share the shape:
    #
    #   mpox_C3L   gene `C3L`  - poxvirus ORF names are <fragment><n><L|R>, which is
    #                            exactly [A-Z]\d{1,4}[A-Z]. The declared symbol was
    #                            dropped before any query existed, so the marker fell
    #                            through to the [All Fields] last resort and was built
    #                            from B5R and haemagglutinin - not the complement
    #                            control protein, and useless for the clade call it
    #                            exists to make.
    #   five rows  gene `18S`  - matches \d{1,4}[A-Z]. pfal_18S, ehist_ssu,
    #                            crypto_ssu, nfow_18S, bmand_18S, all `required`.
    #                            They survived on their `SSU rRNA` alias and the rRNA
    #                            title fallback, so they are degraded rather than
    #                            broken - searched without their most precise term.
    #
    # A curator who writes a symbol in the `gene` column means it is the gene. The
    # length filter still applies to both: a one-character name is unsearchable
    # whatever column it came from (see MIN_QUERY_NAME_LEN).
    if gene and gene != "-" and len(gene) >= MIN_QUERY_NAME_LEN:
        out.append(gene)
    for name in parse_aliases(row.get("aliases", "")):
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
    if strategy == "gene_family" and names and clade_relaxation_allowed(row):
        # Divergence is the point of a family marker: a member that is only
        # deposited under a sister species (BoNT/E in C. butyricum, Stx in
        # Shigella) still belongs in the family, so relax the organism to genus.
        # Restricted to `pathotype` rows - see clade_relaxation_allowed(). The
        # nucleotide tier needs the gate as much as the protein tier does:
        # henipa_cedar's nucleotide set came from the Paramyxoviridae rung and
        # was eight Measles N CDSs.
        genus = genus_of(taxon)
        lin = _LINEAGE.get(taxon.lower())
        clades = [genus] if genus else []
        if lin:
            clades += [c for c in (lin[1], lin[2]) if c]   # family, order
        for clade in dedup_str(clades):
            for name in names:
                out.append(f'{name}[All Fields] AND "{clade}"[Organism]')
    if looks_like_rrna_or_organelle(row):
        # rRNA and organelle loci are not indexed under [Gene Name]; they appear in
        # titles as e.g. "18S ribosomal RNA gene". Without these the queries come back
        # empty, which on the first build silently removed a `required` marker and made
        # Balamuthia mandrillaris unconfirmable.
        for name in (names or [gene]):
            out.append(f'{name}[Title] AND {org}')
        out.append(f'ribosomal RNA[Title] AND {org}')
        out.append(f'small subunit ribosomal RNA[Title] AND {org}')
        out.append(f'mitochondrion[Title] AND {org}')
        out.append(f'complete genome[Title] AND {org}')
    return dedup_str(out)


# taxon name (lowercased) -> (genus, family, order), from rankedlineage.dmp.
# Populated once in main(); empty if the file is absent.
_LINEAGE: dict = {}


def load_rankedlineage(taxdump: Path) -> dict:
    """Parse rankedlineage.dmp into name -> (genus, family, order).

    WHY THIS IS NEEDED

    The HMM tier is the remote-homology tier - the one that catches a divergent or
    novel allele, and the reason the marker catalogue lists `hmm` for the filovirus
    and BoNT families at all. A profile needs at least three sequences to mean
    anything, and the first full build produced only 1-2 for ebola_NP, ebolaS_NP,
    marburg_NP and cbut_bont, because queries relaxed no further than the GENUS. A
    single orthoebolavirus species has almost no distinct RefSeq NP proteins.

    Relaxing to the FAMILY (Filoviridae) is what makes a profile buildable, and family
    is precisely the level at which remote homology should be modelled anyway - a
    novel filovirus is a new member of the family, not of the species.

    rankedlineage.dmp carries the name in column 2, so no join against names.dmp is
    required.
    """
    path = taxdump / "rankedlineage.dmp"
    out: dict = {}
    if not path.exists():
        LOG.warning("no rankedlineage.dmp under %s - gene_family queries cannot relax "
                    "to family level, so HMM profiles for divergent families may not "
                    "reach the 3-sequence minimum", taxdump)
        return out

    # taxid -> (genus, family, order), plus the scientific name as a key.
    by_taxid: dict = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) < 7 or not f[0].isdigit():
                continue
            lin = (f[3], f[4], f[5])          # genus, family, order
            by_taxid[int(f[0])] = lin
            if f[1]:
                out[f[1].lower()] = lin

    # ALSO key on synonyms and equivalent names.
    #
    # rankedlineage.dmp carries only the current SCIENTIFIC name, and ICTV has renamed
    # most of the virus families this catalogue targets - "Marburg marburgvirus" is now
    # a synonym, so a scientific-name-only lookup found no lineage for it and silently
    # produced zero family-level queries. That is the same failure the watchlist avoids
    # by resolving names through every accepted name class, and the fix is the same.
    names_dmp = taxdump / "names.dmp"
    n_syn = 0
    if names_dmp.exists():
        with names_dmp.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                f = [x.strip() for x in line.split("|")]
                if len(f) < 4 or not f[0].isdigit():
                    continue
                key = f[1].lower()
                if key in out:
                    continue
                lin = by_taxid.get(int(f[0]))
                if lin:
                    out[key] = lin
                    n_syn += 1
    else:
        LOG.warning("no names.dmp under %s - synonym lookups unavailable, so renamed "
                    "taxa will get no family-level relaxation", taxdump)

    LOG.info("loaded lineage for %d taxon names (%d via synonyms)", len(out), n_syn)
    return out


# Markers whose gene term is an rRNA or organelle designation rather than a gene
# symbol. [Gene Name] does not index these - rRNA loci are annotated as e.g.
# "18S ribosomal RNA" - which is why bmand_18S, nfow_18S and pjiro_mtLSU queries
# came back empty on the first build. bmand_18S is a `required` marker, so its
# absence made Balamuthia mandrillaris impossible to confirm at all.
RRNA_TERMS = ("18s", "16s", "23s", "28s", "5.8s", "ssu", "lsu", "rrna", "its",
              "mtlsu", "mito", "mitochondrial", "rdna")


def looks_like_rrna_or_organelle(row: dict) -> bool:
    blob = f"{row.get('gene','')} {row.get('aliases','')} {row.get('locus','')}".lower()
    return any(t in blob for t in RRNA_TERMS)


SUBTYPE_STRATEGY = "subtype_in_organism"


def subtype_tokens(row: dict) -> List[str]:
    """Tokens a retrieved record MUST carry, for a subtype-scoped marker.

    WHY A THIRD KIND OF SCOPE EXISTS.

    clade_relaxation_allowed() asks "does this rung stay inside the declared
    organism". For flua_HA_H5 the organism is right and the SUBTYPE is the
    claim, so that question passes a reference set that is exactly as wrong as
    the clade rungs it blocks. Measured on the shipped pack: `flua_HA_H5` and
    `flua_HA_H7` held the SAME seven HA proteins - H5N1, H1N1, H3N2, H1N1,
    H2N2, H9N2, H7N9 - one of seven being H5 and one being H7, both `built`
    with an HMM over the mixture. They are the two markers for highly
    pathogenic avian influenza, and they are `report_always`, so any seasonal
    influenza A sample put "H5 present" and "H7 present" in front of an analyst.

    The winning rung for both was `HA[Gene Name] AND "Influenza A virus"
    [Organism] AND srcdb_refseq[Properties]`; the H5/H7 aliases contributed
    nothing because the gene symbol's rungs come first and already satisfied
    `want`.

    NCBI encodes the subtype in the organism string - "Influenza A virus
    (A/goose/Guangdong/1/1996(H5N1))" - so it is checkable on the record rather
    than guessable from the query. A row declares the requirement by setting
    query_strategy=subtype_in_organism; the aliases are the tokens.
    """
    if (row.get("query_strategy") or "").strip() != SUBTYPE_STRATEGY:
        return []
    return parse_aliases(row.get("aliases", ""))


def _subtype_re(token: str) -> "re.Pattern[str]":
    """H5 matches H5 and H5N1, but not H15N1, H51 or AH5."""
    return re.compile(r"(?<![A-Za-z0-9])" + re.escape(token) + r"(?:N\d{1,2})?(?![0-9])")


def record_has_subtype(rec: SeqRecord, tokens: Sequence[str]) -> bool:
    """True when the record's own text carries one of the required subtypes.

    No tokens means no requirement, so every non-subtype marker is unaffected.
    """
    if not tokens:
        return True
    text = f"{rec.description} {rec.accession}"
    return any(_subtype_re(t).search(text) for t in tokens)


def clade_relaxation_allowed(row: dict) -> bool:
    """May this marker be searched OUTSIDE its declared organism?

    Only for a `pathotype` marker, where the gene is the finding and the host
    species is not: BoNT in C. butyricum, Stx in Shigella. `cbot_bont` is built
    on exactly that - `marker_rules.yaml` pairs it with `veto_if: cbut_bont` so
    a cross-species hit is ATTRIBUTED rather than suppressed (D29).

    For every other class the organism IS the claim, and a sequence from a
    sister taxon turns the marker into a detector for the wrong thing. Measured
    on the shipped pack: relaxing `nipah_L` ("polymerase", class `required`) to
    Paramyxoviridae returned seven non-Nipah polymerases out of eight, and
    `henipa_cedar` ("nucleoprotein", class `absence_expected`) returned eight
    Measles nucleoproteins and no Cedar at all. Since `marker_evidence.py`
    scores a marker PRESENT if ANY tier passes, those are not diluted profiles -
    they are false PRESENT calls on a Nipah marker and on the veto that would
    stand down a Nipah finding.

    No marker in the shipped pack needs the relaxation to build a correct set:
    every `gene_family` row that is `built` filled from in-organism rungs alone,
    `cbot_bont` included. So this gate costs nothing that was working, and the
    two markers it leaves short are reported as gaps - which is the answer the
    rest of this module is built to give when evidence is absent.
    """
    return (row.get("marker_class") or "").strip().lower() == "pathotype"


def aa_query_plan(row: dict) -> List[Tuple[str, bool]]:
    """The protein ladder as (query, safe-to-accumulate) pairs.

    WHY A RUNG NEEDS A SAFETY BIT.

    fetch_protein() walks the ladder until it has enough sequences, and
    everything it collects is written to the marker's protein FASTA - the
    reference set reads are aligned against. `marker_evidence.py` states the
    rule plainly: "A marker is PRESENT if ANY applicable tier passes." So a
    sequence from the wrong organism in that file is not a diluted profile, it
    is a false PRESENT on a high-consequence marker.

    Two kinds of rung cannot be accumulated, both MEASURED on the shipped pack:

      * CLADE-RELAXED rungs. The genus/family/order relaxation was written to
        get a profile past the 3-sequence minimum, and for a generic product
        name it leaves the organism entirely. Accumulating them gave
        `henipa_cedar` (gene N, "nucleoprotein") EIGHT Measles virus
        nucleoproteins and zero Cedar, and gave `nipah_L` ("polymerase") seven
        non-Nipah paramyxovirus polymerases out of eight - Langya, Achimota,
        Respirovirus, two Parajeilongvirus, a bat morbillivirus. henipa_cedar
        is the `veto_if` for BOTH Nipah and Hendra, so that set would have let
        a measles-positive sample veto a genuine Nipah detection.

      * [All Fields] rungs, even inside the right organism. The term matches
        strain designations, not just products: `bont[All Fields] AND
        "Clostridium butyricum"[Organism]` matches every protein of strain
        "C. butyricum E4 str. BoNT E BL5262", so `cbut_bont` collected a
        thiazole biosynthesis protein, an MFS transporter and a
        3-hydroxybutyryl-CoA dehydrogenase alongside the real BoNT/E.

    Both stay in the ladder as a LAST RESORT - a marker that can only be served
    by a loose rung still gets what it gets today, one rung, no accumulation -
    but neither may contribute to a growing pool.
    """
    taxon, gene = row["taxon"], row["gene"]
    strategy = row["query_strategy"]
    names = query_names(row)
    org = organism_clause(taxon)
    out: List[Tuple[str, bool]] = []
    for name in names:
        # Own organism, matched on a NAME field: precise enough to pool.
        out.append((f'{name}[Gene Name] AND {org} AND srcdb_refseq[Properties]', True))
        out.append((f'{name}[Protein Name] AND {org} AND srcdb_refseq[Properties]', True))
        out.append((f'{name}[Gene Name] AND {org}', True))
        out.append((f'{name}[Protein Name] AND {org}', True))
        out.append((f'{name}[All Fields] AND {org}', False))
    if strategy == "gene_family" and names and clade_relaxation_allowed(row):
        genus = genus_of(taxon)
        lin = _LINEAGE.get(taxon.lower())
        clades = [genus] if genus else []
        if lin:
            clades += [c for c in (lin[1], lin[2]) if c]   # family, order
        for clade in dedup_str(clades):
            for name in names:
                out.append((f'{name}[Protein Name] AND "{clade}"[Organism]', False))
                out.append((f'{name}[All Fields] AND "{clade}"[Organism]', False))
    if strategy == "named_region":
        for name in names or [gene]:
            out.append((f'{name}[All Fields] AND {org}', False))
    seen, uniq = set(), []
    for q, safe in out:
        if q and q not in seen:
            seen.add(q)
            uniq.append((q, safe))
    return uniq


def aa_queries(row: dict) -> List[str]:
    """The protein ladder as plain strings (dry-run, diagnostics, self-tests)."""
    return [q for q, _ in aa_query_plan(row)]


def dedup_str(items: Iterable[str]) -> List[str]:
    out: List[str] = []
    for i in items:
        if i and i not in out:
            out.append(i)
    return out


# Minimum characters a CDS gene name must have before it may claim a LONGER wanted
# name as a prefix. See gene_match_rank().
MIN_PREFIX_CLAIM = 3


def gene_match_rank(cds_gene: str, wanted: Sequence[str]) -> Optional[int]:
    """How well a CDS [gene=] value matches gene + aliases. Lower is better.

        0  exact, case-insensitive
        1  the CDS name EXTENDS a wanted name  ('capBCADE' for a wanted 'capB')
        2  the CDS name is a PREFIX of a wanted name  ('pag' for a wanted 'pagA')
        None  no match

    Rank 1 is the permissiveness this function was written for and it is kept:
    NCBI annotation of one locus varies between submitters, so 'capB' must match
    'capBCADE' and 'stx2' must match 'stx2A'. The cost of a loose match there is a
    slightly larger marker set; the cost of a strict one is an empty one.

    RANK 2 IS WHERE IT WENT WRONG. The old code accepted `w.startswith(g)` with no
    length condition, so a one-character CDS gene name could claim any alias that
    merely begins with that letter - and this pack pairs single-letter viral gene
    names with word aliases throughout. Measured on the shipped marker database
    2026-10-20: `nipah_L` declares gene `L` with alias `polymerase`, and the Nipah
    **P** (phosphoprotein) CDS matched it because "polymerase".startswith("p"). The
    nucleotide reference for the L-polymerase marker was therefore 2130 bp of P, while
    the protein tier held the real 2244 aa L - two tiers, two different genes, and
    min_coverage measured against the wrong one. nipah_L sits in Nipah's
    confirm_requires. Six markers share the L/polymerase shape (cchf_L, marburg_L,
    lassa_L, junin_L, rvfv_L, nipah_L); nipah_L is the one that took the CDS route and
    so the one that was hit. MIN_PREFIX_CLAIM makes a genuine truncation ('pag', 'inv')
    still match while a single letter cannot.

    Returning a RANK rather than a bool also lets the caller prefer an exact hit over a
    prefix one, which is the other half of the problem: `ypest_pmt1_caf1` declares
    `caf1` and its nucleotide reference holds caf1A, caf1M and caf1R - the operon's
    accessory genes - and not caf1 itself, because every one of them is a legitimate
    rank-1 match and nothing preferred the gene actually asked for.
    """
    g = (cds_gene or "").strip().lower()
    if not g:
        return None
    best: Optional[int] = None
    for w in wanted:
        w = (w or "").strip().lower()
        if not w:
            continue
        if g == w:
            return 0
        if g.startswith(w):
            best = 1 if best is None else min(best, 1)
        elif w.startswith(g) and len(g) >= MIN_PREFIX_CLAIM:
            best = 2 if best is None else min(best, 2)
    return best


def gene_name_matches(cds_gene: str, wanted: Sequence[str]) -> bool:
    """Whether a CDS [gene=] value matches gene + aliases at all."""
    return gene_match_rank(cds_gene, wanted) is not None


# ---------------------------------------------------------------------------
# Retrieval
# ---------------------------------------------------------------------------

def _tried_phrase(ladder: Sequence[str], reported: str) -> str:
    """Describe a retrieval that came back empty.

    The note used to read "last query: <q>" while the fetchers returned
    ladder[0] - the FIRST rung. That wording is why these rows read as "one
    query was tried and it was a bad one" when in fact the whole ladder had
    been walked and returned nothing. Both facts now appear.
    """
    if not ladder:
        return "no query could be built from the catalogue row"
    return (f"{len(ladder)} quer{'y' if len(ladder) == 1 else 'ies'} tried, "
            f"all empty; first: {reported or ladder[0]}")


def fetch_protein(ez: Entrez, row: dict, want: int) -> Tuple[List[SeqRecord], str]:
    """Protein records for a marker. Returns (records, queries-used).

    THE LADDER IS WALKED UNTIL IT HAS ENOUGH, NOT UNTIL IT HAS ANY.

    aa_queries() is ordered strictest-first and, for a `gene_family` row, ends
    with genus/family/order-level rungs whose stated purpose is to "get a
    profile past the 3-sequence minimum". Returning at the first query that
    yielded ANY record made that tail unreachable: a virus species has exactly
    one RefSeq protein per locus, so rung 1 always returned 1 record and the
    search stopped there. Every `gene_family` virus marker in the shipped pack
    was therefore starved at 1-2 sequences and its HMM tier - the tier the
    catalogue says exists to catch a divergent allele - could never be built.
    The relaxation was written, committed and unreachable.

    So accumulate across rungs and stop when `want` distinct records are in
    hand - but ONLY across the rungs aa_query_plan() marks poolable, which are
    the ones that stay inside the declared organism and match on a name field.
    Pooling the clade-relaxed and [All Fields] rungs was tried, measured, and
    rejected: see aa_query_plan() for the three markers it corrupted. Those
    rungs remain as a single-rung last resort, so a marker that only they can
    serve behaves exactly as it did before pooling existed.

    `fetch_nucleotide` is deliberately NOT given the same treatment: the
    nucleotide tier has no minimum-count threshold to be starved against, and
    the two markers missing it (orthopox_ATI, mpox_C3L) get nothing from any
    rung - that is an empty ladder, not a short-circuited one, and no amount of
    accumulation fixes it.
    """
    plan = aa_query_plan(row)
    subtypes = subtype_tokens(row)
    pool: List[SeqRecord] = []
    used: List[str] = []
    n_wrong_subtype = 0

    def take(query: str) -> int:
        nonlocal n_wrong_subtype
        uids = ez.esearch("protein", query)
        if not uids:
            return 0
        text = ez.efetch("protein", uids[:max(want * 4, want)], "fasta")
        records = [r for r in parse_fasta(text) if AA_MIN_LEN <= len(r) <= AA_MAX_LEN]
        if subtypes:
            kept = [r for r in records if record_has_subtype(r, subtypes)]
            n_wrong_subtype += len(records) - len(kept)
            records = kept
        if not records:
            return 0
        before = len(pool)
        # Dedup against everything already collected, not just within this rung:
        # the later rungs re-find the earlier rungs' records by construction.
        pool[:] = dedup_records(pool + records)
        gained = len(pool) - before
        if gained:
            used.append(query)
            LOG.debug("  protein: +%d record(s) (%d total) via %s",
                      gained, len(pool), query)
        return gained

    # Pass 1: pool across the rungs that stay inside the declared organism and
    # match on a name field. This is the pass that fills an HMM profile.
    for query, poolable in plan:
        if len(pool) >= want:
            break
        if poolable:
            take(query)

    # Pass 2: LAST RESORT. Nothing precise existed, so fall back to the loose and
    # clade-relaxed rungs - and take only the FIRST that yields anything, which is
    # exactly what this function did for every marker before pooling existed. A
    # marker served only by a loose rung is therefore unchanged, and a loose rung
    # can never top up an otherwise-good set with off-target sequence.
    if not pool:
        for query, poolable in plan:
            if not poolable and take(query):
                break

    if not pool:
        if subtypes and n_wrong_subtype:
            LOG.warning("  %s: %d protein record(s) matched the gene and organism "
                        "but carried none of the required subtype(s) %s - not used",
                        row.get("marker_id"), n_wrong_subtype, "/".join(subtypes))
        return [], (plan[0][0] if plan else "")
    if subtypes and n_wrong_subtype:
        LOG.info("  %s: dropped %d off-subtype protein record(s), kept %d for %s",
                 row.get("marker_id"), n_wrong_subtype, len(pool), "/".join(subtypes))
    # Prefer RefSeq accessions (NP_/YP_/WP_/XP_) - curated, and stable across
    # rebuilds, which matters because the accession goes in the manifest.
    pool.sort(key=lambda r: (0 if re.match(r"^[NYWX]P_", r.accession) else 1,
                             -len(r)))
    return pool[:want], " + ".join(used)


def fetch_nucleotide(ez: Entrez, row: dict,
                     want: int) -> Tuple[List[SeqRecord], str, str]:
    """Nucleotide records for a marker.

    Returns (records, query-used, why-empty). The third element is non-empty
    only when the ladder MATCHED records and every one was filtered out; see
    the counters below for why that distinction is load-bearing.

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
    subtypes = subtype_tokens(row)
    n_wrong_subtype = 0
    # WHY THESE COUNTERS EXIST: "the search found nothing" and "the search found
    # forty complete genomes and every one was discarded" produced the SAME note.
    # orthopox_ATI and mpox_C3L are the second case - C3L[All Fields] returns 40
    # MPXV genomes of ~195 kb, all above NT_MAX_DIRECT_LEN, whose CDSs are named
    # MPXV-UK_P2-0NN rather than C3L - and they read as the first, which points a
    # curator at the query when the actual need is the ORF designation.
    n_uids = n_oversize = n_cds_seen = 0
    for query in nt_queries(row):
        uids = ez.esearch("nuccore", query)
        if not uids:
            continue
        n_uids += len(uids)
        probe = uids[:max(want * 3, want)]
        records: List[SeqRecord] = []

        # Route 1 - CDS extraction. Skipped for named_region/no-gene rows, where
        # there is no gene symbol to match a CDS against.
        if strategy != "named_region" and wanted_names:
            cds_text = ez.efetch("nuccore", probe, "fasta_cds_na")
            ranked: List[Tuple[int, SeqRecord]] = []
            for rec in parse_fasta(cds_text):
                n_cds_seen += 1
                rank = gene_match_rank(rec.gene, wanted_names)
                if rank is None:
                    continue
                if not record_has_subtype(rec, subtypes):
                    n_wrong_subtype += 1
                    continue
                if NT_MIN_LEN <= len(rec) <= NT_MAX_DIRECT_LEN:
                    ranked.append((rank, rec))
            # AN EXACT GENE-NAME MATCH WINS OVER A PREFIX ONE. Without this, every
            # member of an operon is an equally good answer to a query for one of its
            # genes, and which ones survive is decided by the later sort and truncation
            # to `want`. That is how ypest_pmt1_caf1 - the F1 capsular antigen marker in
            # Y. pestis's confirm_requires - ended up holding caf1A, caf1M and caf1R and
            # not caf1. The permissiveness is still there as a FALLBACK: when no CDS
            # names the gene exactly, the prefix matches are all that is on offer and
            # they are still used.
            if ranked:
                best_rank = min(r for r, _ in ranked)
                records = [rec for r, rec in ranked if r == best_rank]

        # Route 2 - whole record, only if it is locus-sized.
        if not records:
            text = ez.efetch("nuccore", probe, "fasta")
            for rec in parse_fasta(text):
                if not record_has_subtype(rec, subtypes):
                    n_wrong_subtype += 1
                    continue
                if NT_MIN_LEN <= len(rec) <= NT_MAX_DIRECT_LEN:
                    records.append(rec)
                else:
                    n_oversize += 1
                    LOG.debug("  nt: rejecting %s (%d bp; not locus-sized)",
                              rec.accession, len(rec))

        if not records:
            continue
        records = dedup_records(records)
        records.sort(key=lambda r: (0 if re.match(r"^(N[CGMRTWZ]|X[MR])_", r.accession) else 1,
                                    -len(r)))
        LOG.debug("  nt: %d record(s) via %s", len(records), query)
        return records[:want], query, ""
    ladder = nt_queries(row)
    if n_uids and (n_oversize or n_cds_seen):
        # Found records, discarded all of them. Say which filter did it.
        why = (f"{len(ladder)} quer{'y' if len(ladder) == 1 else 'ies'} tried; "
               f"they MATCHED {n_uids} record(s), but none was usable as a locus: "
               f"{n_oversize} record(s) exceeded {NT_MAX_DIRECT_LEN} bp (whole "
               f"genomes, not genes) and none of {n_cds_seen} annotated CDS "
               f"matched {wanted_names or 'any name'}. THIS IS NOT A MISSING "
               "SEQUENCE - it is a naming mismatch: supply the ORF/locus "
               "designation this organism's deposits actually use as an alias. "
               f"First query: {ladder[0] if ladder else ''}")
        return [], (ladder[0] if ladder else ""), why
    return [], (ladder[0] if ladder else ""), ""


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
        for ca, cs in zip(ac, as_, strict=True):
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


# ---------------------------------------------------------------------------
# Freshness: does the built database still cover markers.tsv?
# See "EXISTENCE IS NOT FRESHNESS" in the module docstring for why each of the
# four comparisons is here and what it catches that the others do not.
# ---------------------------------------------------------------------------

CHECK_FRESH = 0
CHECK_UNDETERMINED = 1
CHECK_STALE = 3

# Which manifest column claims which tier, and which combined artifact has to
# hold the marker if the claim is true. One table, so a tier cannot be checked
# in one direction and forgotten in the other.
TIER_CLAIMS = (
    ("dna", "n_nt_seqs", "markers.nt.fna"),
    ("protein", "n_aa_seqs", "markers.aa.faa"),
    ("hmm", "has_hmm", "markers.hmm"),
)

# An index and the file it was built from. An index merely EXISTING says nothing:
# D45's lesson on the competitive bundles, which is why mtime is compared.
INDEX_SOURCES = (
    ("markers.nt.sr.mmi", "markers.nt.fna"),
    ("markers.nt.map-ont.mmi", "markers.nt.fna"),
    ("markers.hmm.h3f", "markers.hmm"),
    ("markers.hmm.h3i", "markers.hmm"),
    ("markers.hmm.h3m", "markers.hmm"),
    ("markers.hmm.h3p", "markers.hmm"),
)


def _manifest_claims_tier(row: dict, column: str) -> bool:
    """Does this manifest row claim the tier that `column` records?

    `has_hmm` is a yes/no word while the other two are counts, so the two forms
    are read here rather than at each call site.
    """
    raw = (row.get(column) or "").strip()
    if column == "has_hmm":
        return raw.lower() in ("yes", "true", "1")
    try:
        return int(float(raw or 0)) > 0
    except ValueError:
        return False


def marker_ids_in_artifact(path: Path) -> Optional[Set[str]]:
    """The marker_ids a combined artifact actually holds, or None if unreadable.

    Reads FASTA headers and HMM NAME lines through marker_id_of(), the same
    parse marker_evidence.py attributes a hit with, so this cannot disagree with
    the pipeline about what a header means. A line that is not one of ours is
    skipped rather than guessed at.
    """
    if not path.exists():
        return None
    found: Set[str] = set()
    try:
        with path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    mid = marker_id_of(line)
                elif line.startswith("NAME"):
                    parts = line.split(None, 1)
                    mid = marker_id_of(parts[1]) if len(parts) == 2 else None
                else:
                    continue
                if mid:
                    found.add(mid)
    except OSError:
        return None
    return found


def check_marker_db_freshness(catalogue_path: Path, out_dir: Path,
                              manifest_path: Path,
                              rebuild: str = "scripts/build_marker_db.py --pack <pack>"
                              ) -> int:
    """Read-only. Is the built marker database current for this catalogue?

    Prints what disagrees and returns an exit code: CHECK_FRESH,
    CHECK_STALE, or CHECK_UNDETERMINED when the question cannot be answered.
    Nothing here writes, so it is safe to run against a live database.
    """
    def say(msg: str) -> None:
        print(msg, flush=True)

    if not catalogue_path.exists():
        say(f"UNDETERMINED: no catalogue at {catalogue_path}")
        return CHECK_UNDETERMINED
    if not out_dir.exists():
        say(f"STALE: no marker database at {out_dir}")
        say(f"verdict: STALE. Build it: {rebuild}")
        return CHECK_STALE
    if not manifest_path.exists():
        # The sequence files may well be intact, but the manifest is what
        # decides tested-vs-untested, so without it every marker is UNKNOWN.
        say(f"STALE: no {manifest_path.name} in {out_dir}. marker_evidence.py "
            f"treats it as authoritative, so every marker would be UNKNOWN.")
        say(f"verdict: STALE. Rebuild: {rebuild}")
        return CHECK_STALE

    try:
        catalogue = read_marker_catalogue(catalogue_path)
    except Exception as exc:                              # noqa: BLE001
        say(f"UNDETERMINED: cannot read catalogue {catalogue_path}: {exc}")
        return CHECK_UNDETERMINED

    try:
        with manifest_path.open(encoding="utf-8") as fh:
            man_rows = [r for r in csv.DictReader(fh, delimiter="\t")
                        if (r.get("marker_id") or "").strip()]
    except OSError as exc:
        say(f"UNDETERMINED: cannot read {manifest_path}: {exc}")
        return CHECK_UNDETERMINED

    if not man_rows:
        say(f"STALE: {manifest_path.name} holds no marker row.")
        say(f"verdict: STALE. Rebuild: {rebuild}")
        return CHECK_STALE

    cat_ids = [r["marker_id"] for r in catalogue if r.get("marker_id")]
    manifest = {r["marker_id"].strip(): r for r in man_rows}
    cat_set, man_set = set(cat_ids), set(manifest)

    proved: List[str] = []       # a set disagreement: the database IS stale
    suggested: List[str] = []    # only a timestamp: it MAY be stale

    # ---- 1. the predicate proper: does the database cover the catalogue? ----
    missing = [m for m in cat_ids if m not in man_set]
    if missing:
        proved.append(
            f"{len(missing)} catalogue marker(s) have NO manifest row "
            f"({', '.join(missing[:6])}{', ...' if len(missing) > 6 else ''}), "
            f"so they evaluate UNKNOWN and can never confirm a finding")
    orphans = sorted(man_set - cat_set)
    if orphans:
        proved.append(
            f"{len(orphans)} manifest row(s) are no longer in the catalogue "
            f"({', '.join(orphans[:6])}{', ...' if len(orphans) > 6 else ''})")

    # ---- 2. mtime: an EDITED marker, which the sets above cannot see --------
    try:
        cat_mtime = catalogue_path.stat().st_mtime
        man_mtime = manifest_path.stat().st_mtime
        if cat_mtime > man_mtime:
            suggested.append(
                f"the catalogue is NEWER than the manifest, so a marker may have "
                f"been edited under a marker_id the database already knows - "
                f"a changed accession or gene is invisible to the set comparison")
    except OSError:
        pass

    # ---- 3. the manifest's claims against the artifacts themselves ---------
    # Read an artifact only when some row claims that tier, so a pack with no
    # HMM tier does not pay to scan a profile file that is legitimately absent.
    unreadable: List[str] = []
    for tier, column, filename in TIER_CLAIMS:
        claimants = [m for m, r in manifest.items()
                     if _manifest_claims_tier(r, column)]
        if not claimants:
            continue
        present = marker_ids_in_artifact(out_dir / filename)
        if present is None:
            unreadable.append(filename)
            proved.append(
                f"{len(claimants)} manifest row(s) claim the {tier} tier but "
                f"{filename} is missing or unreadable")
            continue
        absent = sorted(m for m in claimants if m not in present)
        if absent:
            proved.append(
                f"{len(absent)} manifest row(s) claim the {tier} tier but "
                f"{filename} does not hold them "
                f"({', '.join(absent[:6])}{', ...' if len(absent) > 6 else ''}) - "
                f"these report ABSENT rather than UNKNOWN, which asserts the "
                f"marker was searched for")
        # The other direction: sequence in the index for a marker the manifest
        # records as having none. Nothing downstream reads it, but the manifest
        # and the database disagree, and the manifest is what is trusted.
        stray = sorted(m for m in present
                       if m in manifest and not _manifest_claims_tier(
                           manifest[m], column))
        if stray:
            proved.append(
                f"{filename} holds {len(stray)} marker(s) the manifest records "
                f"as having no {tier} sequence "
                f"({', '.join(stray[:6])}{', ...' if len(stray) > 6 else ''})")

    # ---- 4. indexes against the files they were built from -----------------
    for index_name, source_name in INDEX_SOURCES:
        index, source = out_dir / index_name, out_dir / source_name
        if not source.exists():
            continue
        if not index.exists():
            proved.append(f"{index_name} is missing, but {source_name} is present")
            continue
        try:
            if index.stat().st_mtime < source.stat().st_mtime:
                proved.append(
                    f"{index_name} is OLDER than {source_name}, so searches use "
                    f"the previous sequences")
        except OSError:
            pass

    if not proved and not suggested:
        say(f"verdict: FRESH. All {len(cat_ids)} catalogue marker(s) have a "
            f"manifest row, every claimed tier is present in the combined "
            f"artifact that has to hold it, no index predates its source, and "
            f"the catalogue is not newer than the build.")
        return CHECK_FRESH

    for line in proved + suggested:
        say(f"STALE: {line}")
    # "is" when a set proves it, "may be" when only a timestamp suggests it.
    # Overstating a timestamp as proof is how a guard loses its reader.
    if proved:
        say(f"verdict: STALE. Rebuild: {rebuild}")
        return CHECK_STALE
    say(f"verdict: MAY BE STALE on timestamps alone. Rebuild to be certain: "
        f"{rebuild}")
    return CHECK_STALE


def blank_manifest_row(row: dict, status: str, note: str,
                       retrieval_errors: int = 0) -> dict:
    """A manifest row with no sequence: explicit status, empty paths."""
    return {
        "marker_id": row["marker_id"], "taxon": row["taxon"], "gene": row["gene"],
        "marker_class": row["marker_class"], "evidence_type": row["evidence_type"],
        "query_strategy": row["query_strategy"], "status": status,
        "n_nt_seqs": "0", "n_aa_seqs": "0", "has_hmm": "no", "accessions": "",
        "nt_path": "", "aa_path": "", "hmm_path": "", "md5": "",
        "retrieved_utc": utc_now(), "retrieval_errors": str(retrieval_errors),
        "note": note,
    }


def write_curation_todo(path: Path, rows: Sequence[dict]) -> None:
    """Unmet curation tasks, as a file rather than only a log line.

    curate_required and motif rows are not failures to retry - they are work
    items for a human. Writing them out means the gap is reviewable and
    diffable, instead of scrolling past in a build log.

    THE FILTER WAS `status in (curation_required, motif_todo)`, WHICH MISSED TWO
    WHOLE CATEGORIES. On the shipped human pack it listed 13 markers out of 31.

    `failed` was never listed at all - 5 markers whose retrieval failed, which
    evaluate UNKNOWN downstream exactly like a curate_required row, which is the
    stated reason the file exists. Three of them are load-bearing in
    marker_rules.yaml: `bthai_ara` is the veto for both B. pseudomallei and
    B. mallei, `henipa_cedar` is the veto for both Nipah and Hendra, and
    `bruc_S19_eri` is half of the Brucella vaccine rule. A veto that can never
    fire is the quiet kind of gap - it costs specificity, not sensitivity, so
    nothing crashes and no call is obviously wrong.

    `partial` was excluded on the reasoning that those rows DO evaluate, on the
    tiers they have. They do - and on fewer tiers than markers.tsv declares, with
    the missing one recorded nowhere but a `note` string. That hid 13 markers,
    11 missing the HMM tier, which the catalogue says exists to catch a divergent
    allele: cbut_bont, ebola_NP, ebolaS_NP, marburg_NP, marburg_GP, junin_GPC,
    rvfv_L, rvfv_Gn, nipah_L, nipah_G, eeev_nsP. orthopox_ATI (Variola) and
    mpox_C3L (Monkeypox) are missing the nucleotide tier instead.

    A partially built marker is a weaker instrument than the catalogue claims,
    which is a curation task whether or not it also produces an UNKNOWN.
    """
    tmp = path.with_name(path.name + ".part")
    cols = ["marker_id", "taxon", "gene", "marker_class", "evidence_type",
            "query_strategy", "status", "missing_tiers", "what_is_needed",
            "catalogue_notes"]
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        fh.write("# bioThreatFinder unmet marker-curation tasks\n")
        fh.write(f"# written {utc_now()} by build_marker_db.py\n")
        fh.write("# Two kinds of row, told apart by `status`:\n")
        fh.write("#   curation_required / motif_todo / failed - NOT testable at all.\n")
        fh.write("#     Evaluates UNKNOWN downstream: neither present nor absent.\n")
        fh.write("#   partial - testable, but on fewer tiers than markers.tsv declares.\n")
        fh.write("#     Evaluates PRESENT/ABSENT on the tiers in the database, and an\n")
        fh.write("#     ABSENT is correspondingly weaker. `missing_tiers` says which\n")
        fh.write("#     were never searched.\n")
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
                # A RESUMED ROW MUST STILL DESCRIBE THE CURRENT CATALOGUE.
                #
                # This returned dict(old) verbatim - the prior row's evidence_type,
                # marker_class and taxon included - and never compared it against the
                # catalogue row being processed. Two consequences, both silent:
                #
                #  * Change a marker from `dna` to `dna+protein+hmm` in markers.tsv to
                #    add the remote-homology tier D3 relies on, rerun without --force,
                #    and the marker is SKIPped: no protein is fetched, the manifest
                #    keeps evidence_type=dna and status=built, while
                #    marker_evidence.py reads evidence_type from the CATALOGUE and so
                #    expects protein and hmm. The manifest says built, so the marker is
                #    treated as testable, and the two missing tiers evaluate as ABSENCE
                #    rather than UNKNOWN.
                #  * Promote a marker to `required` and rerun: the carried-forward row
                #    still says the old class, and the blocking exit code filters on the
                #    MANIFEST's marker_class - so the build exits 0 with a
                #    now-required marker unbuilt.
                #
                # A fresh run and a resumed run produced different manifests from
                # identical inputs. Anything the catalogue declares is compared here.
                # Only columns the manifest actually carries - copy_expectation is
                # catalogue-only, and comparing it would report drift on every marker.
                drifted = [k for k in ("taxon", "marker_class", "evidence_type",
                                       "query_strategy")
                           if (row.get(k) or "").strip()
                           != (old.get(k) or "").strip()]
                if drifted:
                    LOG.warning("  %s: catalogue changed since this was built (%s) - "
                                "rebuilding rather than carrying the stale row forward",
                                marker_id,
                                ", ".join(f"{k}: {old.get(k)!r} -> {row.get(k)!r}"
                                          for k in drifted))
                elif not old.get("md5") or old["md5"] == got:
                    LOG.info("  SKIP %s: %s already built (%s), md5 OK",
                             marker_id, ", ".join(p.name for p in declared),
                             old.get("status"))
                    return dict(old)
                elif old.get("md5"):
                    LOG.warning("  %s: md5 changed on disk (manifest %s, actual %s) - "
                                "rebuilding", marker_id, old["md5"], got)
            else:
                LOG.warning("  %s: manifest says %s but file(s) are missing or empty "
                            "- rebuilding", marker_id, old.get("status"))

    want = MAX_SEQS_FAMILY if strategy == "gene_family" else MAX_SEQS_SINGLE
    # Snapshot so an incomplete search cannot be written as a completed one.
    errors_before = ctx.ez.n_errors

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
    nt_query = aa_query = nt_why = ""
    accessions: List[str] = []

    if want_dna:
        try:
            nt_records, nt_query, nt_why = fetch_nucleotide(ctx.ez, row, want)
        except EntrezError as exc:
            notes.append(f"nucleotide retrieval error: {exc}")
        if nt_records:
            write_marker_fasta(nt_path, marker_id, nt_records)
            LOG.info("  %s: nt %d seq(s) -> %s", marker_id, len(nt_records), nt_path.name)
        else:
            notes.append("no usable nucleotide sequence; " +
                         (nt_why or _tried_phrase(nt_queries(row), nt_query)))

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
            notes.append("no usable protein sequence; "
                         f"{_tried_phrase(aa_queries(row), aa_query)}")

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

    n_errors = ctx.ez.n_errors - errors_before

    if not obtained and n_errors:
        # THE SEARCH DID NOT COMPLETE, SO ITS EMPTINESS MEANS NOTHING.
        #
        # Reporting `failed` here would record a transient 429 or a truncated
        # efetch as "NCBI has nothing", put the marker on the curation list, and
        # send a human looking for sequence that may well exist. Downstream this
        # is UNKNOWN either way; the difference is what the operator is told, and
        # whether a rebuild can be mistaken for a result.
        notes.insert(0, f"NOT ESTABLISHED - {n_errors} Entrez request(s) FAILED "
                        f"for this marker (last: {ctx.ez.last_error}). Nothing "
                        "here is evidence about what NCBI holds; re-run when the "
                        "service is reachable.")
        LOG.error("  %s: %d Entrez request(s) failed - recorded as %s, NOT as a "
                  "negative result", marker_id, n_errors, STATUS_RETRIEVAL_ERROR)
        return blank_manifest_row(row, STATUS_RETRIEVAL_ERROR, " | ".join(notes),
                                  retrieval_errors=n_errors)

    if not obtained:
        # A ROW THAT BUILT NO QUERY WAS NEVER SEARCHED, AND MUST NOT SAY `failed`.
        #
        # `failed` is documented as "retrieval was attempted and returned nothing
        # usable" - evidence that NCBI has nothing. A row whose gene symbol is one
        # character with no alias (henipa_cedar: gene 'N', aliases '-') is filtered
        # out by query_names(), so nt_queries()/aa_queries() return [] and not a
        # single request is issued. Recording that as `failed` is the same
        # conflation D26 removed from the tier reporting: silence from an
        # instrument that was never switched on, filed as a negative result.
        #
        # It is a catalogue authoring defect - the fix is an alias, not a retry -
        # so it is `curation_required`, which is where a human task belongs.
        unsearchable = (strategy not in ("curate_required", "motif")
                        and not nt_queries(row) and not aa_queries(row))
        if unsearchable:
            status = STATUS_CURATION
            notes.insert(0, "NOT SEARCHED - no query could be built from this "
                            f"catalogue row: gene {row.get('gene')!r} with aliases "
                            f"{row.get('aliases')!r} yields no searchable name "
                            f"(names under {MIN_QUERY_NAME_LEN} characters are "
                            "dropped). Add an alias, as every other single-letter "
                            "virus gene in the catalogue does.")
            LOG.error("  %s: NOT SEARCHED - catalogue row builds no query; add an "
                      "alias for gene %r", marker_id, row.get("gene"))
        else:
            status = STATUS_FAILED
        if skipped_tiers:
            notes.append("tier(s) skipped by flag: " + "; ".join(skipped_tiers))
    elif set(requested) <= set(obtained):
        status = STATUS_BUILT
    else:
        status = STATUS_PARTIAL
        if n_errors:
            notes.insert(0, f"CAUTION: {n_errors} Entrez request(s) failed during "
                            "this marker's retrieval, so a missing tier below may "
                            "be a transient service failure rather than absent "
                            "sequence")
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

    if status in (STATUS_FAILED, STATUS_CURATION):
        return blank_manifest_row(row, status, " | ".join(notes) or
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
        "retrieval_errors": str(n_errors),
        "note": " | ".join(notes),
    }


# ---------------------------------------------------------------------------
# Combined outputs
# ---------------------------------------------------------------------------

DEDUPE_MIN_IDENTITY = 0.99


def _read_fasta_records(path: Path) -> List[Tuple[str, str]]:
    """(header_without_gt, sequence) pairs, sequence upper-cased and unwrapped."""
    out: List[Tuple[str, str]] = []
    name: Optional[str] = None
    buf: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:], []
        else:
            buf.append(line.strip().upper())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def dedupe_marker_records(records: Sequence[Tuple[str, str]]
                          ) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    """Drop sequences that are redundant with one already kept. -> (kept, dropped).

    WHY REDUNDANCY INSIDE ONE MARKER IS HARMFUL, not merely wasteful. Two near-identical
    sequences under the same marker make every read that covers them map equally well to
    both, so minimap2 assigns MAPQ 0 and the `-q 10` filter discards exactly the reads
    that prove the marker is present. Measured on the plumbing control: ypest_ppcp1_pla
    holds two 939 bp alleles differing by ONE base (identity 0.9989) and kept 16 of ~85
    expected reads. It also inflates the breadth denominator when breadth is taken over
    the summed length.
    77 of 143 markers were in this state - mers_ORF1a held 8 sequences reducing to 1,
    one of them an exact duplicate - so the nucleotide tier has been systematically
    under-calling PRESENT across more than half the catalogue.

    Three redundancy tests, all chosen from measurement rather than assumed:
      * identical sequence;
      * one contained in another as a substring - ypest_pmt1_ymt's 1599 bp and 474 bp
        entries are both inside its 1764 bp entry;
      * equal length and identity >= DEDUPE_MIN_IDENTITY.
    The LONGEST is kept, so the representative is the most informative one.

    Genuine diversity is preserved: ypest_pmt1_caf1's three sequences are different genes
    of the F1 operon (2502/906/777 bp, no containment) and all survive, as do the
    gene-family markers whose spread is the point.
    """
    order = sorted(range(len(records)), key=lambda i: -len(records[i][1]))
    kept_idx: List[int] = []
    dropped: List[Tuple[str, str]] = []
    for i in order:
        seq = records[i][1]
        why = None
        for k in kept_idx:
            other = records[k][1]
            if seq == other:
                why = "identical"
            elif seq and other and (seq in other or other in seq):
                why = "substring"
            elif len(seq) == len(other) and seq:
                same = sum(1 for a, b in zip(seq, other) if a == b)
                if same / len(seq) >= DEDUPE_MIN_IDENTITY:
                    why = f"identity={same / len(seq):.4f}"
            if why:
                break
        if why:
            dropped.append((records[i][0], why))
        else:
            kept_idx.append(i)
    kept = [records[i] for i in sorted(kept_idx)]
    return kept, dropped


def concat_files(sources: Sequence[Path], dest: Path, wrap: int = 60,
                 dedupe: bool = True) -> int:
    """Concatenate per-marker FASTAs into the combined file, atomically.

    One source file is one marker, which is what makes de-duplication possible here:
    redundancy is only ever collapsed WITHIN a marker, never across markers.

    The per-marker files on disk are deliberately left untouched - they are the archive,
    their md5s are what the manifest records, and dropping a sequence from them would
    lose provenance. This combined file is the working index, so it is the right place to
    keep one representative per redundant group. See dedupe_marker_records() for why the
    redundancy is harmful rather than merely wasteful.
    """
    n = 0
    names: dict[str, str] = {}
    dupes: list[tuple[str, str, str]] = []
    total_dropped = 0
    tmp = dest.with_name(dest.name + ".part")
    with tmp.open("w", encoding="utf-8", newline="\n") as out:
        for src in sources:
            records = _read_fasta_records(src)
            if dedupe:
                records, dropped = dedupe_marker_records(records)
                if dropped:
                    total_dropped += len(dropped)
                    LOG.info("  dedupe %s: %d -> %d sequence(s) [%s]", src.stem,
                             len(records) + len(dropped), len(records),
                             ", ".join(f"{h.split('|')[2] if h.count('|') >= 2 else h}"
                                       f" {w}" for h, w in dropped[:3])
                             + (", ..." if len(dropped) > 3 else ""))
            for header, seq in records:
                nm = sam_name_of(header)
                if nm in names:
                    dupes.append((nm, names[nm], src.name))
                else:
                    names[nm] = src.name
                n += 1
                out.write(">" + header + "\n")
                for i in range(0, len(seq), wrap):
                    out.write(seq[i:i + wrap] + "\n")
    # REFUSE TO PUBLISH A FASTA AN ALIGNER CANNOT INDEX.
    #
    # A duplicate reference NAME (the text before the first whitespace) makes
    # `samtools sort` fail with "Duplicate entry ... in sam header", which killed
    # MARKER_ALIGN_NT and with it the whole nucleotide marker tier. Per-marker writes
    # disambiguate, but two DIFFERENT markers could still collide, and this is the
    # check against the written artifact rather than trust in the writer.
    if dupes:
        tmp.unlink(missing_ok=True)
        detail = "; ".join(f"{nm} (in {a} and {b})" for nm, a, b in dupes[:5])
        raise SystemExit(
            f"ERROR: {dest.name} would contain {len(dupes)} duplicate reference "
            f"name(s): {detail}\n"
            f"       An aligner's reference name is only the text before the first "
            f"whitespace, so samtools rejects the SAM header and the tier that reads "
            f"this file dies. Refusing to write it.")
    if total_dropped:
        LOG.info("de-duplicated %d redundant sequence(s) across the catalogue; the "
                 "per-marker archive files are unchanged", total_dropped)
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
    # An accession containing '|' is SANITISED, not rejected. Rejecting it was the
    # original behaviour and it was wrong: NCBI protein fetches routinely return
    # SwissProt-style names like sp|Q306W5.1|POLS_EEEV1, so raising turned an ordinary
    # identifier format into a hard retrieval failure for every affected marker.
    check("sp|ACC|ENTRY yields the real accession",
          sanitise_accession("sp|Q306W5.1|POLS_EEEV1") == "Q306W5.1")
    check("tr|ACC|ENTRY yields the real accession",
          sanitise_accession("tr|A0A1B0|SOMETHING") == "A0A1B0")
    check("a plain accession is untouched",
          sanitise_accession("NC_003131.1") == "NC_003131.1")
    check("other pipes are replaced, not dropped",
          sanitise_accession("weird|thing|here") == "weird_thing_here")
    _h = make_marker_header("ok_id", "sp|A1|B2", "desc | with pipe")
    _p = parse_marker_header(_h)
    check("header round-trips a sanitised accession and a piped description",
          _p is not None and _p[0] == "ok_id" and _p[1] == "A1"
          and _p[2] == "desc | with pipe", str(_p))

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
        # The second pagA record is a SUBSTRING of the first, so it is redundant and the
        # combined index keeps one representative. The per-marker archive above still
        # holds both, which the round-trip checks just asserted.
        n = concat_files([out, second], combined)
        check("concat collapses a redundant sequence within a marker", n == 2, str(n))
        ids = [marker_id_of(ln) for ln in combined.read_text(encoding="utf-8").splitlines()
               if ln.startswith(">")]
        check("combined file attributes each sequence to its marker",
              ids == ["banth_pxo1_pagA", "cbot_bont"], str(ids))
        check("the LONGEST representative is the one kept",
              "GNLEVIQAEV" in combined.read_text(encoding="utf-8"))
        # dedupe=False must preserve everything, for callers that want the raw archive.
        combined_raw = tmpdir / "markers.raw.faa"
        check("dedupe can be turned off",
              concat_files([out, second], combined_raw, dedupe=False) == 3)

        print("-- de-duplication rules --")
        kept, dropped = dedupe_marker_records([
            ("m|a|A1|x", "ACGTACGTAA"), ("m|a|A2|x", "ACGTACGTAA")])
        check("identical sequences collapse", len(kept) == 1 and len(dropped) == 1,
              f"{len(kept)}/{len(dropped)}")
        kept, dropped = dedupe_marker_records([
            ("m|a|A1|x", "ACGTACGTAACCGG"), ("m|a|A2|x", "GTACGTAA")])
        check("a contained sequence collapses",
              len(kept) == 1 and dropped and dropped[0][1] == "substring",
              f"{kept}/{dropped}")
        # 1 mismatch in 1000 = 0.999, above the 0.99 floor -> redundant.
        a = "A" * 1000
        b = "A" * 999 + "C"
        kept, dropped = dedupe_marker_records([("m|a|A1|x", a), ("m|a|A2|x", b)])
        check("equal-length near-identity collapses (the pla case)",
              len(kept) == 1 and dropped and dropped[0][1].startswith("identity="),
              f"{dropped}")
        # 5% divergence must be KEPT: that is real allelic diversity.
        c = "A" * 950 + "C" * 50
        kept, dropped = dedupe_marker_records([("m|a|A1|x", a), ("m|a|A2|x", c)])
        check("genuinely divergent sequences are kept (the caf1 case)",
              len(kept) == 2 and not dropped, f"{len(kept)}/{dropped}")
        # Different lengths, no containment -> different genes, keep both.
        kept, dropped = dedupe_marker_records([
            ("m|a|A1|x", "ACGT" * 100), ("m|a|A2|x", "TTGCAA" * 40)])
        check("different genes of one operon are kept", len(kept) == 2, str(len(kept)))

        print("-- missing tiers --")
        # The curation to-do list used to carry only the rows that cannot be
        # evaluated at all, so a marker declaring an HMM it does not have was a
        # silent degradation - 13 of them on the shipped human pack, 11 missing
        # exactly the tier the catalogue says catches a divergent allele.
        check("a declared tier with no sequences is reported missing",
              missing_tiers({"evidence_type": "dna+protein+hmm", "n_nt_seqs": "4",
                             "n_aa_seqs": "2", "has_hmm": "no"}) == ["hmm"],
              str(missing_tiers({"evidence_type": "dna+protein+hmm",
                                 "n_nt_seqs": "4", "n_aa_seqs": "2",
                                 "has_hmm": "no"})))
        check("the nucleotide gap is caught too (the orthopox_ATI case)",
              missing_tiers({"evidence_type": "dna+protein", "n_nt_seqs": "0",
                             "n_aa_seqs": "1", "has_hmm": "no"}) == ["dna"],
              str(missing_tiers({"evidence_type": "dna+protein", "n_nt_seqs": "0",
                                 "n_aa_seqs": "1", "has_hmm": "no"})))
        check("a fully built marker reports no gap",
              missing_tiers({"evidence_type": "dna+protein+hmm", "n_nt_seqs": "3",
                             "n_aa_seqs": "3", "has_hmm": "yes"}) == [])
        # An UNDECLARED tier that happens to exist is not a gap. has_hmm=yes on a
        # dna+protein marker is spare capacity, not a missing instrument.
        check("a tier that exists but is not declared is not a gap",
              missing_tiers({"evidence_type": "dna+protein", "n_nt_seqs": "3",
                             "n_aa_seqs": "3", "has_hmm": "yes"}) == [])
        check("order follows the catalogue declaration",
              missing_tiers({"evidence_type": "dna+protein+hmm", "n_nt_seqs": "0",
                             "n_aa_seqs": "0", "has_hmm": "no"})
              == ["dna", "protein", "hmm"])

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
        check("status taxonomy is exactly the documented six",
              set(ALL_STATUSES) == {"built", "partial", "curation_required",
                                    "failed", "motif_todo", "retrieval_error"})
        check("a failed request is a DIFFERENT status from an empty result",
              STATUS_RETRIEVAL_ERROR != STATUS_FAILED)

        print("-- a failed request is not an empty result --")

        class _FlakyEz:
            """esearch raises the way _request does when NCBI refuses."""

            def __init__(self):
                self.n_requests = 0
                self.n_cache_hits = 0
                self.n_errors = 0
                self.last_error = ""

            def esearch(self, db, term):
                # Mirrors the real client: swallow, count, return empty.
                self.n_errors += 1
                self.last_error = f"esearch {term!r}: HTTP 429"
                return []

            def efetch(self, db, uids, rettype):
                return ""

        flaky = _FlakyEz()
        famrow2 = {"marker_id": "ebola_NP", "taxon": "Zaire ebolavirus", "gene": "NP",
                   "aliases": "nucleoprotein", "query_strategy": "gene_family",
                   "marker_class": "required"}
        recs_f, _ = fetch_protein(flaky, famrow2, MIN_SEQS_FOR_HMM)
        check("a throttled ladder still returns no records", recs_f == [])
        check("every failed request is counted, not swallowed",
              flaky.n_errors == len(aa_queries(famrow2)),
              "%d errors for %d rungs" % (flaky.n_errors, len(aa_queries(famrow2))))
        check("the client records what failed", "429" in flaky.last_error,
              flaky.last_error)

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
           "aliases": "botA;botE", "query_strategy": "gene_family",
           "marker_class": "pathotype"}
    check("a pathotype gene_family row relaxes to the genus",
          any('"Clostridium"[Organism]' in x for x in nt_queries(fam)))
    # The gate. BoNT legitimately sits in a sister species; a nucleoprotein does
    # not, and relaxing one gave henipa_cedar eight Measles sequences.
    notfam = dict(fam, marker_id="nipah_L", taxon="Nipah henipavirus", gene="L",
                  aliases="polymerase", marker_class="required")
    check("a non-pathotype gene_family row does NOT leave its organism",
          all('"Nipah henipavirus"[Organism]' in x for x in nt_queries(notfam)),
          str([x for x in nt_queries(notfam)
               if '"Nipah henipavirus"[Organism]' not in x]))
    check("the same gate applies to the protein ladder",
          all('"Nipah henipavirus"[Organism]' in x for x in aa_queries(notfam)),
          str([x for x in aa_queries(notfam)
               if '"Nipah henipavirus"[Organism]' not in x]))
    check("clade relaxation is keyed on marker_class, not strategy",
          clade_relaxation_allowed(fam) and not clade_relaxation_allowed(notfam))
    reg = {"marker_id": "mtb_RD9", "taxon": "Mycobacterium tuberculosis", "gene": "RD9",
           "aliases": "region of difference 9", "query_strategy": "named_region"}
    check("named_region does not use [Gene Name]",
          all("[Gene Name]" not in x for x in nt_queries(reg)))
    check("one-character aliases are not searched (would match everything)",
          "t" not in query_names({"gene": "lpnA", "aliases": "tul4;t"}),
          str(query_names({"gene": "lpnA", "aliases": "tul4;t"})))
    check("a DECLARED GENE is never dropped as a substitution",
          query_names({"gene": "C3L", "aliases": "complement control protein"})
          == ["C3L", "complement control protein"],
          str(query_names({"gene": "C3L", "aliases": "complement control protein"})))
    check("the 18S gene symbol survives too (five required parasite markers)",
          query_names({"gene": "18S", "aliases": "SSU rRNA"}) == ["18S", "SSU rRNA"],
          str(query_names({"gene": "18S", "aliases": "SSU rRNA"})))
    check("poxvirus ORF names have the substitution SHAPE, which is why "
          "the filter may not see the gene column",
          SUBSTITUTION_RE.match("C3L") and SUBSTITUTION_RE.match("18S"))
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

    print("-- the protein ladder pools precise rungs only --")

    class _LadderEz:
        """Entrez stub: each query maps to a list of (accession, sequence)."""

        def __init__(self, responses):
            self.responses = responses
            self.searched = []

        def esearch(self, db, term):
            self.searched.append(term)
            return ["u%d" % i for i in range(len(self.responses.get(term, [])))]

        def efetch(self, db, uids, rettype):
            term = self.searched[-1]
            return "".join(">%s x\n%s\n" % (acc, seq)
                           for acc, seq in self.responses.get(term, []))

    famrow = {"marker_id": "ebola_NP", "taxon": "Zaire ebolavirus", "gene": "NP",
              "aliases": "nucleoprotein", "query_strategy": "gene_family"}
    plan = aa_query_plan(famrow)
    poolable = [q for q, ok in plan if ok]
    loose = [q for q, ok in plan if not ok]
    aa = "M" + "A" * 80

    check("the plan marks some rungs poolable and some not", poolable and loose,
          "%d poolable, %d loose" % (len(poolable), len(loose)))
    check("no poolable rung leaves the declared organism",
          all('"Zaire ebolavirus"[Organism]' in q for q in poolable),
          str([q for q in poolable if '"Zaire ebolavirus"[Organism]' not in q]))
    check("no poolable rung is an [All Fields] match",
          all("[All Fields]" not in q for q in poolable),
          str([q for q in poolable if "[All Fields]" in q]))
    check("aa_queries() still returns every rung as a plain string",
          aa_queries(famrow) == [q for q, _ in plan])

    # THE FIX: one RefSeq record on rung 1 must not end the search, because the
    # other precise in-organism rungs are where a virus's allele diversity is.
    ez = _LadderEz({poolable[0]: [("NP_000001.1", aa)],
                    poolable[2]: [("NP_000001.1", aa),        # re-found, deduped
                                  ("AAA00002.1", aa + "C"),
                                  ("AAA00003.1", aa + "DE")]})
    recs, used = fetch_protein(ez, famrow, MIN_SEQS_FOR_HMM)
    check("a 1-record first rung no longer ends the search",
          len(recs) >= MIN_SEQS_FOR_HMM, "%d record(s)" % len(recs))
    check("records re-found by a later rung are not duplicated",
          len({r.accession for r in recs}) == len(recs),
          str([r.accession for r in recs]))
    check("RefSeq still sorts first across the whole pool",
          recs[0].accession == "NP_000001.1", recs[0].accession)
    check("every contributing rung is named in the query note",
          poolable[0] in used and poolable[2] in used, used)

    # THE GUARD: henipa_cedar / nipah_L. A clade-relaxed rung is rich - it is
    # where measles and the bat paramyxoviruses live - and it must NOT be allowed
    # to top up a short in-organism set. Short is the correct answer; the marker
    # is then reported as a gap instead of being built from the wrong organism.
    ez2 = _LadderEz({poolable[0]: [("NP_000001.1", aa)],
                     loose[-1]: [("XXX1.1", aa + "C"), ("XXX2.1", aa + "DE"),
                                 ("XXX3.1", aa + "FG"), ("XXX4.1", aa + "HI")]})
    recs2, _ = fetch_protein(ez2, famrow, MIN_SEQS_FOR_HMM)
    check("a clade-relaxed rung may not top up a short in-organism set",
          len(recs2) == 1, "%d record(s): %s" % (len(recs2), [r.accession for r in recs2]))
    check("no off-organism record reaches a set that had an in-organism hit",
          all(r.accession == "NP_000001.1" for r in recs2),
          str([r.accession for r in recs2]))

    # The cbut_bont shape: [All Fields] inside the right organism still matches
    # the strain designation, so it may not top up either.
    allfields = [q for q in loose if '"Zaire ebolavirus"[Organism]' in q]
    check("an in-organism [All Fields] rung exists and is marked loose", allfields,
          str(loose[:2]))
    if allfields:
        ez3 = _LadderEz({poolable[0]: [("NP_000001.1", aa)],
                         allfields[0]: [("JUNK1.1", "MQQ" + "W" * 60),
                                        ("JUNK2.1", "MQQ" + "Y" * 60)]})
        recs3, _ = fetch_protein(ez3, famrow, MIN_SEQS_FOR_HMM)
        check("an [All Fields] rung may not top up a short set either",
              len(recs3) == 1, "%d record(s)" % len(recs3))

    # LAST RESORT unchanged: when nothing precise exists at all, the first loose
    # rung that yields anything is used, exactly as before pooling existed.
    ez4 = _LadderEz({loose[0]: [("L1.1", aa), ("L2.1", aa + "C")],
                     loose[-1]: [("L3.1", aa + "DE")]})
    recs4, used4 = fetch_protein(ez4, famrow, MIN_SEQS_FOR_HMM)
    check("with nothing precise, one loose rung is still used",
          len(recs4) == 2, "%d record(s)" % len(recs4))
    check("the last-resort pass stops at the first productive loose rung",
          loose[-1] not in used4, used4)

    ez5 = _LadderEz({poolable[0]: [("NP_1.1", aa), ("NP_2.1", aa + "C"),
                                   ("NP_3.1", aa + "DE")]})
    recs5, used5 = fetch_protein(ez5, famrow, MIN_SEQS_FOR_HMM)
    check("a satisfied first rung stops the ladder", ez5.searched == [poolable[0]],
          str(ez5.searched))
    check("a satisfied first rung is truncated to want", len(recs5) == MIN_SEQS_FOR_HMM)
    check("a single contributing rung yields a plain query note", used5 == poolable[0],
          used5)

    ez6 = _LadderEz({})
    recs6, used6 = fetch_protein(ez6, famrow, MIN_SEQS_FOR_HMM)
    check("an empty ladder returns nothing", recs6 == [])
    check("an empty ladder walks every rung", len(ez6.searched) == len(plan),
          "%d of %d" % (len(ez6.searched), len(plan)))
    check("the empty-retrieval note counts the rungs rather than claiming one",
          "quer" in _tried_phrase(aa_queries(famrow), used6)
          and str(len(plan)) in _tried_phrase(aa_queries(famrow), used6),
          _tried_phrase(aa_queries(famrow), used6))
    check("a row that builds no query says so, not 'none built' as a query",
          "no query could be built" in _tried_phrase([], ""))

    print("-- subtype scope: the organism is right and the subtype is the claim --")
    h5 = {"marker_id": "flua_HA_H5", "taxon": "Influenza A virus", "gene": "HA",
          "aliases": "H5", "marker_class": "pathotype",
          "query_strategy": "subtype_in_organism"}
    check("a subtype-scoped row declares its tokens", subtype_tokens(h5) == ["H5"],
          str(subtype_tokens(h5)))
    check("a row that is not subtype-scoped requires nothing",
          subtype_tokens(dict(h5, query_strategy="gene_family")) == [])
    check("no requirement means every record passes",
          record_has_subtype(SeqRecord("A", "anything at all", "M" * 50), []))
    def _rec(desc):
        return SeqRecord("X", desc, "M" * 50)
    check("the declared subtype is kept",
          record_has_subtype(_rec("HA [Influenza A virus (A/goose/1/96(H5N1))]"), ["H5"]))
    check("a different N is still the same H",
          record_has_subtype(_rec("HA [Influenza A virus (A/x/20(H5N8))]"), ["H5"]))
    check("another subtype is dropped - this is the seven-protein bug",
          not record_has_subtype(_rec("HA [Influenza A virus (A/x/09(H1N1))]"), ["H5"])
          and not record_has_subtype(_rec("HA [Influenza A virus (A/x/13(H7N9))]"), ["H5"]))
    check("H15 is not H5",
          not record_has_subtype(_rec("HA [Influenza A virus (A/y/20(H15N1))]"), ["H5"]))
    check("H5 is not matched inside a longer number",
          not record_has_subtype(_rec("HA [Influenza A virus (A/z/20(H51))]"), ["H5"]))
    check("a subtype-scoped row does NOT leave its organism",
          all('"Influenza A virus"[Organism]' in q for q in nt_queries(h5))
          and all('"Influenza A virus"[Organism]' in q for q in aa_queries(h5)),
          str([q for q in nt_queries(h5)
               if '"Influenza A virus"[Organism]' not in q]))
    check("subtype_in_organism is a known strategy",
          SUBTYPE_STRATEGY in KNOWN_STRATEGIES)

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
  curation_required  not automatable, a human task. Either query_strategy=
                     curate_required, or a catalogue row that builds no query at all
  failed             retrieval COMPLETED and returned nothing usable; the note records
                     the queries that were used. This is a statement about what NCBI
                     holds, so it is only written when the search actually finished
  retrieval_error    retrieval DID NOT COMPLETE - one or more Entrez requests failed
                     (rate limit, 5xx, truncated read). Says nothing about what NCBI
                     holds. Re-run when the service is reachable; the run exits
                     non-zero so an incomplete build cannot be mistaken for a result
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
                   help="the volume this script builds under; it appends "
                        "bioThreatFinder/packs/<pack>/markers itself, so this is the "
                        "PARENT of $BTF_DB (default: $BTF_DB's parent, "
                        f"here {DEFAULT_DB_ROOT or 'unset - $BTF_DB is not set'})")
    p.add_argument("--taxdump", type=Path, default=DEFAULT_TAXDUMP,
                   help="pinned NCBI taxonomy snapshot, used for rankedlineage.dmp so "
                        "gene_family queries can relax to family level. Should be the "
                        "same snapshot the watchlist was resolved against "
                        "(default: %(default)s)")
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
    p.add_argument("--check", action="store_true",
                   help="read-only: is the built database current for the pack's "
                        "markers.tsv? Exits 0 fresh, 3 stale, 1 undetermined. "
                        "Writes nothing, so it is safe against a live database")
    p.add_argument("--markers-dir", type=Path, default=None, metavar="DIR",
                   help="with --check: the marker directory to judge, instead of the one "
                        "derived from --db-root. The default --db-root is the SHARED db "
                        "root while a pack's marker database is built on the project "
                        "volume, so a bare --check would report a healthy pack as having "
                        "no database at all. Also the way to check a backup directly")
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

    # Refuse rather than guess, and name the two ways to answer. This used to default to a
    # volume on one particular machine, so a clone's build wrote to a path that did not
    # exist - or, worse for --check, read one and reported a healthy pack as having no
    # database. --markers-dir answers the question on its own, so it is accepted instead.
    if args.db_root is None and not (args.check and args.markers_dir is not None):
        raise SystemExit(
            "ERROR: --db-root is not set and $BTF_DB is not in the environment.\n"
            "  $BTF_DB is where the pipeline reads its databases (nextflow.config refuses\n"
            "  to start without it). --db-root is its PARENT, because this script appends\n"
            "  bioThreatFinder/packs/<pack>/markers itself.\n"
            "    export BTF_DB=/path/to/project/bioThreatFinder\n"
            "  or:\n"
            "    scripts/build_marker_db.py --pack PACK --db-root /path/to/project\n"
            "  To check a marker directory without either, pass --markers-dir DIR --check.")

    pack_dir = args.packs_dir / args.pack
    catalogue_path = pack_dir / "markers.tsv"
    # db_root is None only on the --markers-dir --check path, which overrides all three
    # below; the refusal above has already rejected every other way of getting here.
    base = args.db_root / "bioThreatFinder" if args.db_root else None
    out_dir = base / "packs" / args.pack / "markers" if base else None
    manifest_path = out_dir / "markers.manifest.tsv" if out_dir else None
    todo_path = out_dir / "curation_required.tsv" if out_dir else None
    cache_dir = base / "cache" / "entrez" if base else None

    # --markers-dir judges a directory directly, which is what the launch guard points at:
    # it knows exactly where the database is, so its advice should not depend on a default
    # that differs from it. Refused outside --check rather than silently ignored - a build
    # writes a dozen artifacts and honouring this for only some of them would be worse.
    if args.markers_dir is not None:
        if not args.check:
            raise SystemExit("ERROR: --markers-dir applies to --check only")
        out_dir = args.markers_dir
        manifest_path = out_dir / "markers.manifest.tsv"

    # Before the taxdump load, which is minutes and gigabytes and which a
    # read-only freshness question has no use for.
    if args.check:
        # The remediation must be RUNNABLE as printed, in a shell that may not have
        # $BTF_DB exported - so the resolved --db-root is always named rather than left to
        # a default. It used to be carried only when it differed from the default, which
        # was correct while the default was a fixed volume and is not now that it comes
        # from the environment.
        rebuild = f"scripts/build_marker_db.py --pack {args.pack}"
        if args.db_root is not None:
            rebuild += f" --db-root {args.db_root}"
        return check_marker_db_freshness(catalogue_path, out_dir, manifest_path, rebuild)

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

    global _LINEAGE
    _LINEAGE = load_rankedlineage(args.taxdump)

    rows = read_marker_catalogue(catalogue_path)
    LOG.info("catalogue: %d marker row(s)", len(rows))

    # MERGE, never truncate, when this is a subset run.
    #
    # write_manifest() rewrites the whole file. With --only-marker that replaced a
    # 163-row manifest with the 19 rows just processed, destroying the record of
    # everything built earlier - and because marker_evidence.py treats the manifest as
    # AUTHORITATIVE for tested-vs-untested, every absent row silently becomes UNKNOWN.
    # A subset rebuild would therefore have taken the entire marker tier offline while
    # the sequence files sat intact on disk.
    existing_manifest: dict = {}
    if args.only_marker or args.only_taxon:
        mpath = manifest_path
        if mpath.exists():
            try:
                with mpath.open(encoding="utf-8") as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        if r.get("marker_id"):
                            existing_manifest[r["marker_id"]] = r
                LOG.info("subset run: carrying forward %d existing manifest row(s)",
                         len(existing_manifest))
            except OSError as exc:
                LOG.warning("could not read %s to merge (%s); a subset run would "
                            "otherwise truncate it, so refusing to continue", mpath, exc)
                return 2

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
                for q in nq[1:]:
                    LOG.info("[dry-run]       then: %s", q)
            if ("protein" in tiers or "hmm" in tiers) and aq:
                # Every rung, not just the first. fetch_protein walks the ladder
                # until it has enough sequences, so a dry run that showed only
                # rung 1 was describing a retrieval that does not happen - and
                # it is why the family-level relaxation could sit unreachable
                # without anyone reading a plan and noticing.
                LOG.info("[dry-run]     protein: %s", aq[0])
                for q in aq[1:]:
                    LOG.info("[dry-run]       then: %s", q)
                if "hmm" in tiers:
                    LOG.info("[dry-run]       (hmm tier needs >= %d distinct protein "
                             "sequences from the rungs above)", MIN_SEQS_FOR_HMM)
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
    # Seed with the carried-forward rows so a subset run updates the manifest instead
    # of replacing it. Rows for markers processed in THIS run are overwritten below.
    manifest_rows: List[dict] = [existing_manifest[k] for k in existing_manifest
                                 if k not in {r["marker_id"] for r in selected}]
    if manifest_rows:
        LOG.info("manifest seeded with %d row(s) from the previous build",
                 len(manifest_rows))

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
        before = md5_or_none(combined_nt)
        n = concat_files(nt_sources, combined_nt)
        nt_changed = md5_or_none(combined_nt) != before
        LOG.info("wrote %s (%d sequence(s) from %d marker(s))%s",
                 combined_nt.name, n, len(nt_sources),
                 " [content changed - indexes will be rebuilt]" if nt_changed else "")
        if tools.get("minimap2"):
            build_mm2_indexes(tools["minimap2"], combined_nt, out_dir, "markers.nt",
                              args.threads, args.force or nt_changed)
        else:
            LOG.warning("minimap2 absent: %s left unindexed (install: %s)",
                        combined_nt.name, TOOL_INSTALL["minimap2"])
    else:
        LOG.warning("no nucleotide sequences at all - %s not written", combined_nt.name)

    if aa_sources:
        before = md5_or_none(combined_aa)
        n = concat_files(aa_sources, combined_aa)
        aa_changed = md5_or_none(combined_aa) != before
        LOG.info("wrote %s (%d sequence(s) from %d marker(s))%s",
                 combined_aa.name, n, len(aa_sources),
                 " [content changed - database will be rebuilt]" if aa_changed else "")
        if tools.get("mmseqs"):
            build_mmseqs_db(tools["mmseqs"], combined_aa, out_dir / "mmseqs" / "markers_aa",
                            args.force or aa_changed)
        else:
            LOG.warning("mmseqs absent: no MMseqs2 database built (install: %s)",
                        TOOL_INSTALL["mmseqs"])
    elif not args.skip_protein:
        LOG.warning("no protein sequences at all - %s not written", combined_aa.name)

    if hmm_sources:
        before = md5_or_none(combined_hmm)
        tmp = combined_hmm.with_name(combined_hmm.name + ".part")
        with tmp.open("w", encoding="utf-8", newline="\n") as out:
            for src in hmm_sources:
                out.write(src.read_text(encoding="utf-8", errors="replace"))
        tmp.replace(combined_hmm)
        hmm_changed = md5_or_none(combined_hmm) != before
        LOG.info("wrote %s (%d profile(s))%s", combined_hmm.name, len(hmm_sources),
                 " [content changed - will be re-pressed]" if hmm_changed else "")
        if tools.get("hmmpress"):
            press_hmms(tools["hmmpress"], combined_hmm, args.force or hmm_changed)
        else:
            LOG.warning("hmmpress absent: %s not pressed (install: %s)",
                        combined_hmm.name, TOOL_INSTALL["hmmpress"])
    elif not args.skip_hmm:
        LOG.warning("no profile HMMs built - %s not written", combined_hmm.name)

    write_manifest(manifest_path, manifest_rows)

    # ---- curation TODO file ---------------------------------------------
    todo_rows = []
    n_untestable = n_degraded = 0
    for r in manifest_rows:
        gaps = missing_tiers(r)
        untestable = r["status"] in (STATUS_CURATION, STATUS_MOTIF, STATUS_FAILED)
        # A `built` row with a gap would mean the status and the counts disagree;
        # list it anyway rather than trusting the status over the evidence.
        if not untestable and not gaps:
            continue
        src = next((x for x in selected if x["marker_id"] == r["marker_id"]), {})
        if untestable:
            n_untestable += 1
            need = r.get("note", "")
        else:
            n_degraded += 1
            need = (f"tier(s) not built: {', '.join(gaps)} - declared in markers.tsv "
                    f"but absent from the database, so never searched. "
                    f"{r.get('note', '')}").strip()
        todo_rows.append({**r, "missing_tiers": ",".join(gaps),
                          "what_is_needed": need,
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
    LOG.info("curation tasks:    %s (%d row(s): %d not testable, %d built on fewer "
             "tiers than declared)", todo_path, len(todo_rows), n_untestable, n_degraded)
    if n_degraded:
        # Named here because this is the line a build operator actually reads, and
        # the whole reason these rows were invisible is that the detail lived only
        # in a manifest `note`.
        for r in todo_rows:
            if r.get("missing_tiers") and r["status"] not in (
                    STATUS_CURATION, STATUS_MOTIF, STATUS_FAILED):
                LOG.warning("  DEGRADED %s (%s): missing %s",
                            r["marker_id"], r["taxon"], r["missing_tiers"])
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
    # STATUS_RETRIEVAL_ERROR counts here too. The whole point of separating it from
    # STATUS_FAILED is that it is NOT a verdict about NCBI - which makes shipping a
    # pack built from one worse than shipping a pack with a known gap, because the
    # gap is real and the error is an artefact of the network.
    blocking = [r for r in manifest_rows
                if r["status"] in (STATUS_FAILED, STATUS_RETRIEVAL_ERROR)
                and r["marker_class"] == "required"
                and categories.get(r["taxon"].lower()) == "high_consequence"]
    if blocking:
        LOG.error("--- BLOCKING: required marker(s) for high_consequence taxa failed ---")
        for r in blocking:
            LOG.error("  %s (%s) [%s]", r["marker_id"], r["taxon"], r["status"])
        LOG.error("exiting non-zero: these taxa cannot be confirmed at all.")
        return 1

    # A RUN THAT HIT THE NETWORK AND LOST MUST NOT LOOK LIKE A CLEAN RUN, even when
    # every affected marker happened to be low-consequence. Without this, a rebuild
    # throttled end to end exits 0 and the manifest it leaves behind is read as a
    # measurement of what NCBI holds.
    errored = [r for r in manifest_rows if r["status"] == STATUS_RETRIEVAL_ERROR]
    touched = [r for r in manifest_rows
               if (r.get("retrieval_errors") or "0") not in ("", "0")]
    if errored or ez.n_errors:
        LOG.error("--- %d Entrez request(s) FAILED during this run ---", ez.n_errors)
        # NAME THE ROWS, because "re-run the affected markers" is only an
        # instruction if the affected markers are identifiable. A `built` row is
        # skipped by the resume path, so one built under a failing service is
        # otherwise indistinguishable from a clean one.
        for r in errored:
            LOG.error("  NOT ESTABLISHED  %-20s (%s) - retry when reachable",
                      r["marker_id"], r["taxon"])
        for r in touched:
            if r["status"] == STATUS_RETRIEVAL_ERROR:
                continue
            LOG.error("  SUSPECT %-8s %-20s (%s) - %s request(s) failed while it was "
                      "built; it recovered from other rungs, but its set may be "
                      "short", r["status"], r["marker_id"], r["taxon"],
                      r["retrieval_errors"])
        if touched:
            LOG.error("  refresh just these: %s",
                      " ".join(f"--only-marker {r['marker_id']}" for r in touched))
        LOG.error("exiting non-zero: this build may be INCOMPLETE. `retrieval_error` "
                  "rows are retried by a plain re-run; a SUSPECT row is `built`, so "
                  "the resume path SKIPS it and refreshing it needs --force.")
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        LOG.warning("interrupted - rerun to resume (built markers are skipped and the "
                    "manifest was rewritten after the last completed marker)")
        sys.exit(130)
