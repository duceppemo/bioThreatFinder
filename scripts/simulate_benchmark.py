#!/usr/bin/env python3
"""
simulate_benchmark.py - generate synthetic metagenomes with KNOWN ground truth so
                        that bioThreatFinder's detection thresholds can be
                        measured instead of reasoned about.

Why this exists
---------------
Every threshold in this pipeline (min_unique_reads, min_breadth, neighbour_ratio,
the marker gates) is currently a defensible guess. A guess is not a performance
claim. This script builds the datasets that turn them into numbers:

  * sensitivity / limit of detection  - spike a target at a known reads-per-million
    level into a realistic host+background matrix and see where detection fails;
  * specificity                       - spike only the near neighbour (B. cereus,
    F. novicida, Y. pseudotuberculosis, B. thailandensis) and count how often the
    high-consequence cousin is called anyway;
  * vaccine discrimination            - spike BCG / vaccinia / RB51-like material
    and check it lands in `vaccine_strain`, not in disease;
  * marker-gate correctness           - spike a B. anthracis chromosome with
    pXO1/pXO2 REMOVED and confirm the marker gate refuses to promote it;
  * control behaviour                 - NTCs and index hopping, which is where
    low-biomass metagenomics actually generates its false positives.

The ground truth is the whole product. Consequently this script is paranoid in
one specific direction: it will refuse to build a dataset it cannot describe
truthfully. A spike whose genome cannot be resolved is a hard error, never a
silent skip, because a missing spike read out as "the pipeline failed to detect a
true positive" and would corrupt every sensitivity number computed from the batch.

What it writes
--------------
  <outdir>/<batch>/manifest.tsv        ground truth, one row per sample
  <outdir>/<batch>/samplesheet.csv     directly runnable by the Nextflow pipeline
  <outdir>/<batch>/simulate.log        full run log
  <outdir>/<batch>/logs/<sample>.cmds  every command used for that sample
  <outdir>/<batch>/logs/<sample>.json  per-component read counts (the resume state)
  <outdir>/<batch>/reads/<sample>*.fastq.gz
  <outdir>/_genome_cache/              staged (and locus-filtered) plain FASTAs

Dependencies
------------
Python 3.8+ standard library only (matching bin/*.py), plus external simulators:

  Illumina: `iss` (InSilicoSeq)  preferred, `art_illumina` (ART) fallback
  ONT:      `badread`            preferred, `simulator.py` (NanoSim) fallback

Scenario files are parsed by bin/btf_yaml.py - this script deliberately contains
no YAML parsing of its own.

Usage
-----
  simulate_benchmark.py --scenarios benchmarks/scenarios.yaml --batch v1 --dry-run
  simulate_benchmark.py --scenarios benchmarks/scenarios.yaml --batch v1 --threads 32
  simulate_benchmark.py --scenarios benchmarks/scenarios.yaml --batch v1 \
      --scenario banth_lod --platform illumina --replicates 1
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import hashlib
import json
import logging
import os
import random
import re
import shutil
import struct
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

LOG = logging.getLogger("simulate_benchmark")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_DB_ROOT = Path("/media/36tb/db")
CONDA_BIN = Path("/home/bioinfo/miniconda3/bin")

# Headroom required *on top of* the estimated output size. The shared volume runs
# near capacity and a half-written benchmark batch is worse than a refusal: its
# manifest would claim samples that are truncated.
DEFAULT_MIN_FREE_GIB = 20

# The consumer of this batch (the scoring script) reads these columns positionally
# as well as by name. Do not reorder without updating it.
MANIFEST_COLUMNS = [
    "sample", "scenario_id", "kind", "replicate", "platform", "level_rpm", "role",
    "truth_present", "truth_absent", "truth_vaccine", "host_fraction",
    "total_reads", "spiked_reads", "excluded_loci", "fastq_1", "fastq_2",
    "seed", "notes",
]
SAMPLESHEET_COLUMNS = ["sample", "platform", "fastq_1", "fastq_2", "role"]

PLATFORMS = ["illumina", "ont"]
ROLES = ["sample", "control"]

# `kind` drives the semantics of a scenario, so an unrecognised value is an error
# rather than a default: guessing would mislabel ground truth.
KINDS = ["lod_series", "specificity", "vaccine", "contamination", "negative",
         "marker_ablation"]
SPIKE_SOURCES = ["target", "neighbour", "vaccine", "background"]

DEFAULTS = {
    "replicates": 3,
    "platforms": list(PLATFORMS),
    "host": "human",
    "host_fraction": 0.95,
    "total_reads": 20_000_000,
    "ont_total_reads": 4_000_000,
    "read_length": 150,
    "ont_mean_length": 5000,
}
# Used when a spike carries no levels_rpm. 100 rpm is comfortably above every
# threshold in pack.yaml, so it behaves as "clearly present" rather than as an
# accidental limit-of-detection test.
DEFAULT_LEVEL_RPM = 100.0

# An NTC is a real library, just a tiny one. 1% of the sample depth keeps the
# control informative (it is what control subtraction is computed from) without
# doubling the cost of the batch.
NTC_READ_FRACTION = 0.01
NTC_MIN_READS = 2000

# InSilicoSeq's read length is baked into its error model, so the model is chosen
# from the requested read_length rather than passed alongside it.
# VERIFY: model names against `iss generate --help` for the installed version;
# older releases used --model_file and a different name set.
ISS_MODELS = [(125, "hiseq"), (150, "novaseq"), (300, "miseq")]

# Locus name -> regex matched (case-insensitively) against FASTA headers, for
# --exclude-loci style strain surgery. Headers in RefSeq genomic FASTA look like
#   >NC_007322.2 Bacillus anthracis str. Ames plasmid pXO2, complete sequence
# so the bare plasmid name is usually enough; the aliases exist for the loci whose
# published names differ from what the header says.
# VERIFY: header spellings for pMT1/pCD1 vary between assemblies ("plasmid pMT1",
# "plasmid MT1", "pMT1-like"); the patterns below are deliberately loose, and a
# zero-match request is a hard error precisely so a wrong pattern cannot pass
# unnoticed.
LOCUS_ALIASES = {
    "pxo1": r"pxo1",
    "pxo2": r"pxo2",
    "ppcp1": r"p?pcp1|pesticin|\b9\.5\s*kb\b",
    "pmt1": r"p?mt1|murine\s+toxin",
    "pcd1": r"p?cd1|\blcr\b",
    "rd1": r"\brd1\b|esx-?1",
    "bont": r"bont|botulinum\s+neurotoxin",
    "ctxab": r"ctxa|ctxb|cholera\s+toxin",
    "chromosome": r"chromosome",
    "plasmid": r"plasmid",
}

# Bytes per FASTQ record, uncompressed: header+sequence+'+'+quality. Used only for
# the pre-flight estimate.
FASTQ_OVERHEAD_BYTES = 50
GZIP_RATIO = 0.35          # measured typical for simulated FASTQ with quality
ESTIMATE_SAFETY = 1.25     # simulators write intermediates before we concatenate


# ---------------------------------------------------------------------------
# Small helpers (kept consistent with scripts/build_pack_refs.py)
# ---------------------------------------------------------------------------

def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sanitise(name: str) -> str:
    """Filesystem-safe stem. 'Bacillus anthracis' -> 'Bacillus_anthracis'."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name).strip())
    return re.sub(r"_+", "_", s).strip("_") or "unknown"


def md5_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def which(name: str) -> Optional[str]:
    """PATH first, then the shared conda prefix (not always on PATH under
    Nextflow/cron) - same resolution order as build_pack_refs.py."""
    found = shutil.which(name)
    if found:
        return found
    candidate = CONDA_BIN / name
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    return None


def human(nbytes: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if nbytes < 1024 or unit == "TiB":
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} TiB"


def check_free_space(path: Path, need_gib: float) -> None:
    probe = path
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    usage = shutil.disk_usage(probe)
    free_gib = usage.free / (1024 ** 3)
    LOG.info("free space at %s: %.0f GiB (need %.0f GiB)", probe, free_gib, need_gib)
    if free_gib < need_gib:
        raise SystemExit(
            f"ERROR: insufficient free space: {free_gib:.0f} GiB available at {probe}, "
            f"need {need_gib:.0f} GiB.\n"
            f"       Free space, lower --min-free-gib if you know better, reduce "
            f"total_reads/replicates in the scenario file, or point --outdir at a "
            f"larger volume.\n"
            f"       Nothing was generated."
        )


def derive_seed(master: int, *parts) -> int:
    """Per-sample seed = hash(master, scenario, replicate, platform, ...).

    Hashed rather than counted on purpose: a counter would renumber every
    downstream sample the moment a scenario is inserted, so an old batch could
    never be reproduced from a newer scenario file. Truncated to 31 bits because
    ART's -rs and badread's --seed both want a plain positive int.
    """
    key = "|".join([str(master)] + [str(p) for p in parts])
    return int(hashlib.md5(key.encode("utf-8")).hexdigest()[:8], 16) & 0x7FFFFFFF


def level_tag(level: float) -> str:
    """0.1 -> '0p1', 100 -> '100'. Sample names must be filesystem- and TSV-safe."""
    txt = ("%g" % float(level))
    return txt.replace(".", "p").replace("-", "m")


# ---------------------------------------------------------------------------
# The shared YAML parser lives in bin/. Import it rather than reimplementing it:
# scenario files gate what the measured numbers mean, and two parsers would
# eventually disagree about one of them.
# ---------------------------------------------------------------------------

def import_btf_yaml():
    bin_dir = Path(__file__).resolve().parent.parent / "bin"
    if str(bin_dir) not in sys.path:
        sys.path.insert(0, str(bin_dir))
    try:
        import btf_yaml  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            f"ERROR: could not import the shared YAML parser btf_yaml from {bin_dir}\n"
            f"       ({exc})\n"
            f"       This script does not parse YAML itself by design. Expected\n"
            f"       {bin_dir / 'btf_yaml.py'} exposing load(path) -> dict and BtfYamlError.\n"
            f"       Check out the full repository (bin/ must sit beside scripts/)."
        )
    for attr in ("load", "BtfYamlError"):
        if not hasattr(btf_yaml, attr):
            raise SystemExit(
                f"ERROR: {bin_dir / 'btf_yaml.py'} does not expose `{attr}`.\n"
                f"       This script needs load(path) -> dict and BtfYamlError from it."
            )
    return btf_yaml


# ---------------------------------------------------------------------------
# Scenario file: validation
#
# btf_yaml already types scalars, but a scenario file is hand-written and is the
# definition of ground truth, so every field is re-validated here with its own
# path in the error message. Coercion from str is tolerated (a quoted number is a
# plausible typo, not a reason to refuse); a value that cannot be a number at all
# is not.
# ---------------------------------------------------------------------------

class ScenarioError(SystemExit):
    pass


def _bad(where: str, msg: str) -> "ScenarioError":
    return ScenarioError(f"ERROR: scenario file: {where}: {msg}")


def as_mapping(value, where: str) -> dict:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise _bad(where, f"expected a mapping, got {type(value).__name__}")
    return value


def as_list(value, where: str) -> list:
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or isinstance(value, dict):
        raise _bad(where, f"expected a list, got {type(value).__name__}")
    if not isinstance(value, (list, tuple)):
        raise _bad(where, f"expected a list, got {type(value).__name__}")
    return list(value)


def as_str(value, where: str, default: str = "") -> str:
    if value is None:
        return default
    if isinstance(value, (dict, list, tuple)):
        raise _bad(where, f"expected a string, got {type(value).__name__}")
    return str(value).strip()


def as_float(value, where: str, default: Optional[float] = None) -> float:
    if value is None:
        if default is None:
            raise _bad(where, "is required")
        return float(default)
    try:
        return float(value)
    except (TypeError, ValueError):
        raise _bad(where, f"expected a number, got {value!r}")


def as_int(value, where: str, default: Optional[int] = None) -> int:
    if value is None:
        if default is None:
            raise _bad(where, "is required")
        return int(default)
    try:
        return int(float(value))
    except (TypeError, ValueError):
        raise _bad(where, f"expected an integer, got {value!r}")


def as_bool(value, where: str, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in ("true", "yes", "1", "on"):
        return True
    if text in ("false", "no", "0", "off"):
        return False
    raise _bad(where, f"expected a boolean, got {value!r}")


def as_name_list(value, where: str) -> List[str]:
    """Accept either a YAML list or a ';'/',' separated string of names."""
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [as_str(v, where) for v in value if as_str(v, where)]
    return [t.strip() for t in re.split(r"[;,]", str(value)) if t.strip()]


class Spike:
    __slots__ = ("taxon", "source", "genome", "levels_rpm", "exclude_loci")

    def __init__(self, taxon: str, source: str, genome: str,
                 levels_rpm: List[float], exclude_loci: List[str]):
        self.taxon = taxon
        self.source = source
        self.genome = genome
        self.levels_rpm = levels_rpm
        self.exclude_loci = exclude_loci

    def __repr__(self) -> str:
        return f"<Spike {self.taxon} {self.source} {self.levels_rpm}>"


class Scenario:
    __slots__ = ("sid", "description", "kind", "spikes", "present", "absent",
                 "vaccine", "include_ntc", "index_hop_fraction", "replicates",
                 "platforms", "host", "host_fraction", "total_reads",
                 "ont_total_reads", "read_length", "ont_mean_length", "background",
                 "downsample_total_reads")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    @property
    def primary_spike(self) -> Optional[Spike]:
        """The spike whose level the manifest's level_rpm column reports.

        A scenario can carry several spikes (target plus a competing neighbour, a
        vaccine strain alongside wild type). The primary is the declared target,
        because that is the organism whose limit of detection is being measured;
        secondary spikes are held at their first declared level.
        """
        for sp in self.spikes:
            if sp.source == "target":
                return sp
        return self.spikes[0] if self.spikes else None


# RECOGNISED KEYS AT EACH LEVEL OF THE SCENARIO FILE.
#
# The loader had no key-set validation, so an unrecognised key was silently ignored.
# The shipped file declared
#     downsample_total_reads: [40000000, 20000000, ...]
# on the depth-titration scenario, and `grep -rn downsample_total_reads` matched that
# YAML line and nothing else in the repository: the scenario generated ONE depth while
# the file said it measured six, and the run reported success. The same silence would
# swallow `host_fracton: 0.998`, leaving the sample at the default with no warning.
#
# `expect` is listed as a scenario key because evaluate_benchmark.py reads it even
# though this script does not - omitting it would make every real scenario fail.
DOC_KEYS = {"schema_version", "defaults", "scenarios"}
DEFAULTS_KEYS = {"replicates", "platforms", "host", "host_fraction", "total_reads",
                 "ont_total_reads", "read_length", "ont_mean_length", "background"}
SCENARIO_KEYS = {"id", "description", "kind", "spikes", "truth", "controls",
                 "background", "host", "host_fraction", "total_reads",
                 "ont_total_reads", "read_length", "ont_mean_length", "platforms",
                 "replicates", "expect", "downsample_total_reads",
                 # Documentation/metadata read elsewhere: tests/test_pipeline.py asserts
                 # all 7 NYSDOH pathogen categories are covered via nysdoh_category.
                 "nysdoh_category"}
SPIKE_KEYS = {"taxon", "source", "genome", "levels_rpm", "exclude_loci"}
TRUTH_KEYS = {"present", "absent", "vaccine"}
CONTROLS_KEYS = {"include_ntc", "index_hop_fraction"}


def reject_unknown_keys(mapping, allowed: set, where: str) -> None:
    """Fail on any key this loader does not read. Silence here is a missing axis."""
    unknown = sorted(k for k in (mapping or {}) if k not in allowed)
    if unknown:
        raise _bad(where,
                   f"unrecognised key(s) {', '.join(repr(k) for k in unknown)}. "
                   f"Nothing reads them, so whatever they declare would silently not "
                   f"happen while the run reported success. "
                   f"Recognised here: {', '.join(sorted(allowed))}")


def load_scenarios(path: Path) -> Tuple[Scenario, ...]:
    btf_yaml = import_btf_yaml()
    if not path.is_file():
        raise SystemExit(
            f"ERROR: scenario file not found: {path}\n"
            f"       Point --scenarios at the YAML that defines the benchmark "
            f"(see benchmarks/ in the repository)."
        )
    try:
        doc = btf_yaml.load(path)
    except btf_yaml.BtfYamlError as exc:
        raise SystemExit(
            f"ERROR: could not parse scenario file {path}:\n"
            f"       {exc}\n"
            f"       bin/btf_yaml.py supports a restricted YAML subset and refuses "
            f"anything it might mis-read; fix the line it names."
        )
    doc = as_mapping(doc, "<document>")
    reject_unknown_keys(doc, DOC_KEYS, "<document>")

    version = as_int(doc.get("schema_version"), "schema_version", 1)
    if version != 1:
        raise _bad("schema_version",
                   f"this script understands version 1, file declares {version}")

    d = as_mapping(doc.get("defaults"), "defaults")
    reject_unknown_keys(d, DEFAULTS_KEYS, "defaults")
    base = {
        "replicates": as_int(d.get("replicates"), "defaults.replicates",
                             DEFAULTS["replicates"]),
        "platforms": [p.lower() for p in as_name_list(d.get("platforms"),
                                                      "defaults.platforms")]
                     or list(DEFAULTS["platforms"]),
        "host": as_str(d.get("host"), "defaults.host", DEFAULTS["host"]),
        "host_fraction": as_float(d.get("host_fraction"), "defaults.host_fraction",
                                  DEFAULTS["host_fraction"]),
        "total_reads": as_int(d.get("total_reads"), "defaults.total_reads",
                              DEFAULTS["total_reads"]),
        # ONT depth is not the same number as Illumina depth. 20M reads at a 5 kb mean
        # is ~100 Gbp, which no real flow cell delivers - simulating it would measure a
        # machine that does not exist and flatter long-read sensitivity accordingly.
        "ont_total_reads": as_int(d.get("ont_total_reads"), "defaults.ont_total_reads",
                                  DEFAULTS["ont_total_reads"]),
        "read_length": as_int(d.get("read_length"), "defaults.read_length",
                              DEFAULTS["read_length"]),
        "ont_mean_length": as_int(d.get("ont_mean_length"), "defaults.ont_mean_length",
                                  DEFAULTS["ont_mean_length"]),
        "background": {},
    }
    for plat in base["platforms"]:
        if plat not in PLATFORMS:
            raise _bad("defaults.platforms",
                       f"unknown platform {plat!r}; expected one of {', '.join(PLATFORMS)}")
    if not 0.0 <= base["host_fraction"] < 1.0:
        raise _bad("defaults.host_fraction", "must be >= 0 and < 1")
    if base["total_reads"] <= 0:
        raise _bad("defaults.total_reads", "must be positive")
    if base["replicates"] <= 0:
        raise _bad("defaults.replicates", "must be positive")

    for key, val in as_mapping(d.get("background"), "defaults.background").items():
        frac = as_float(val, f"defaults.background[{key}]")
        if not 0.0 <= frac < 1.0:
            raise _bad(f"defaults.background[{key}]", "must be >= 0 and < 1")
        base["background"][as_str(key, "defaults.background key")] = frac

    raw_scenarios = as_list(doc.get("scenarios"), "scenarios")
    if not raw_scenarios:
        raise _bad("scenarios", "no scenarios defined - nothing to simulate")

    out: List[Scenario] = []
    seen: Dict[str, int] = {}
    for i, raw in enumerate(raw_scenarios):
        where = f"scenarios[{i}]"
        s = as_mapping(raw, where)
        reject_unknown_keys(s, SCENARIO_KEYS, where)
        sid = as_str(s.get("id"), f"{where}.id")
        if not sid:
            raise _bad(f"{where}.id", "is required (it names the samples on disk)")
        if sid != sanitise(sid):
            raise _bad(f"{where}.id",
                       f"{sid!r} is not filesystem-safe; use letters, digits, '.', '_', '-'")
        if sid in seen:
            raise _bad(f"{where}.id",
                       f"duplicate id {sid!r} (also at scenarios[{seen[sid]}]) - "
                       f"sample names and truth rows would collide")
        seen[sid] = i

        kind = as_str(s.get("kind"), f"{where}.kind", "lod_series").lower()
        if kind not in KINDS:
            raise _bad(f"{where}.kind",
                       f"unknown kind {kind!r}; expected one of {', '.join(KINDS)}")

        platforms = [p.lower() for p in as_name_list(s.get("platforms"),
                                                     f"{where}.platforms")] \
                    or list(base["platforms"])
        for plat in platforms:
            if plat not in PLATFORMS:
                raise _bad(f"{where}.platforms", f"unknown platform {plat!r}")

        background = dict(base["background"])
        if "background" in s:
            background = {}
            for key, val in as_mapping(s.get("background"),
                                       f"{where}.background").items():
                background[as_str(key, f"{where}.background key")] = as_float(
                    val, f"{where}.background[{key}]")

        spikes: List[Spike] = []
        for j, rawsp in enumerate(as_list(s.get("spikes"), f"{where}.spikes")):
            spwhere = f"{where}.spikes[{j}]"
            sp = as_mapping(rawsp, spwhere)
            reject_unknown_keys(sp, SPIKE_KEYS, spwhere)
            taxon = as_str(sp.get("taxon"), f"{spwhere}.taxon")
            if not taxon:
                raise _bad(f"{spwhere}.taxon",
                           "is required - an unnamed spike cannot be scored")
            source = as_str(sp.get("source"), f"{spwhere}.source", "target").lower()
            if source not in SPIKE_SOURCES:
                raise _bad(f"{spwhere}.source",
                           f"unknown source {source!r}; expected one of "
                           f"{', '.join(SPIKE_SOURCES)}")
            levels = [as_float(v, f"{spwhere}.levels_rpm[]")
                      for v in as_list(sp.get("levels_rpm"), f"{spwhere}.levels_rpm")]
            if not levels:
                levels = [DEFAULT_LEVEL_RPM]
                LOG.warning("%s (%s): spike %r has no levels_rpm; defaulting to "
                            "%g rpm", sid, spwhere, taxon, DEFAULT_LEVEL_RPM)
            for lv in levels:
                if lv <= 0:
                    raise _bad(f"{spwhere}.levels_rpm",
                               f"level {lv!r} must be > 0 (use kind: negative for a "
                               f"sample with nothing spiked in)")
                if lv > 1e6:
                    raise _bad(f"{spwhere}.levels_rpm",
                               f"level {lv!r} exceeds 1e6 rpm (the whole library)")
            if kind != "lod_series" and len(levels) > 1:
                raise _bad(f"{spwhere}.levels_rpm",
                           f"kind {kind!r} takes a single level; got {len(levels)}. "
                           f"Use kind: lod_series for a dose series so that each "
                           f"level becomes its own sample.")
            spikes.append(Spike(
                taxon=taxon,
                source=source,
                genome=as_str(sp.get("genome"), f"{spwhere}.genome"),
                levels_rpm=levels,
                exclude_loci=as_name_list(sp.get("exclude_loci"),
                                          f"{spwhere}.exclude_loci"),
            ))

        if kind == "negative" and spikes:
            raise _bad(f"{where}.spikes",
                       "kind: negative must not declare spikes - the point of the "
                       "row is that the answer is 'nothing'")
        if kind != "negative" and not spikes:
            raise _bad(f"{where}.spikes",
                       f"kind {kind!r} needs at least one spike (or use kind: negative)")
        if kind == "marker_ablation" and not any(sp.exclude_loci for sp in spikes):
            raise _bad(f"{where}.spikes",
                       "kind: marker_ablation must exclude at least one locus, "
                       "otherwise it is an ordinary positive control")

        # IN-SILICO DOWN-SAMPLING: one library, several depths.
        #
        # This key was declared in the shipped scenario file and read by nothing, so
        # downsample_accuracy_floor - whose entire purpose is to find the depth at which
        # the pipeline stops being accurate, and which cites NYSDOH CLEP's endorsement
        # of exactly this design - generated a SINGLE depth while the file said six,
        # and reported success. Implemented rather than deleted: removing the key would
        # have removed a validation axis the design deliberately included.
        depths = [as_int(v, f"{where}.downsample_total_reads[]")
                  for v in as_list(s.get("downsample_total_reads"),
                                   f"{where}.downsample_total_reads")]
        for dv in depths:
            if dv <= 0:
                raise _bad(f"{where}.downsample_total_reads",
                           "every depth must be positive")

        truth = as_mapping(s.get("truth"), f"{where}.truth")
        reject_unknown_keys(truth, TRUTH_KEYS, f"{where}.truth")
        controls = as_mapping(s.get("controls"), f"{where}.controls")
        reject_unknown_keys(controls, CONTROLS_KEYS, f"{where}.controls")
        hop = as_float(controls.get("index_hop_fraction"),
                       f"{where}.controls.index_hop_fraction", 0.0)
        if not 0.0 <= hop < 0.5:
            raise _bad(f"{where}.controls.index_hop_fraction",
                       "must be >= 0 and < 0.5 (it is a bleed fraction, not a mix)")

        host_fraction = as_float(s.get("host_fraction"), f"{where}.host_fraction",
                                 base["host_fraction"])
        if not 0.0 <= host_fraction < 1.0:
            raise _bad(f"{where}.host_fraction", "must be >= 0 and < 1")
        total_reads = as_int(s.get("total_reads"), f"{where}.total_reads",
                             base["total_reads"])
        if total_reads <= 0:
            raise _bad(f"{where}.total_reads", "must be positive")
        ont_total_reads = as_int(s.get("ont_total_reads"), f"{where}.ont_total_reads",
                                 base["ont_total_reads"])
        if ont_total_reads <= 0:
            raise _bad(f"{where}.ont_total_reads", "must be positive")

        out.append(Scenario(
            sid=sid,
            description=as_str(s.get("description"), f"{where}.description"),
            kind=kind,
            spikes=spikes,
            present=as_name_list(truth.get("present"), f"{where}.truth.present"),
            absent=as_name_list(truth.get("absent"), f"{where}.truth.absent"),
            vaccine=as_name_list(truth.get("vaccine"), f"{where}.truth.vaccine"),
            include_ntc=as_bool(controls.get("include_ntc"),
                                f"{where}.controls.include_ntc", False),
            index_hop_fraction=hop,
            replicates=as_int(s.get("replicates"), f"{where}.replicates",
                              base["replicates"]),
            platforms=platforms,
            host=as_str(s.get("host"), f"{where}.host", base["host"]),
            host_fraction=host_fraction,
            total_reads=total_reads,
            ont_total_reads=ont_total_reads,
            read_length=as_int(s.get("read_length"), f"{where}.read_length",
                               base["read_length"]),
            ont_mean_length=as_int(s.get("ont_mean_length"),
                                   f"{where}.ont_mean_length",
                                   base["ont_mean_length"]),
            background=background,
            downsample_total_reads=depths,
        ))

    # Truth cross-check. A scenario that spikes an organism and then lists it as
    # absent is usually a copy-paste error, and it would score every correct
    # detection as a false positive.
    #
    # Two kinds are exempt, because for them "spiked but absent" is the whole
    # point: an ablated strain (anthrax chromosome with no pXO1/pXO2) and a
    # vaccine strain (BCG) both put the taxon's reads in the sample while the
    # correct verdict is still "not the disease agent". Individual spikes carrying
    # exclude_loci are exempt for the same reason.
    for sc in out:
        if sc.kind in ("marker_ablation", "vaccine"):
            continue
        spiked = {sp.taxon for sp in sc.spikes if not sp.exclude_loci}
        contradiction = spiked & set(sc.absent)
        if contradiction:
            raise _bad(f"scenarios[{seen[sc.sid]}].truth.absent",
                       f"{', '.join(sorted(contradiction))} is spiked into this "
                       f"scenario but declared absent - one of the two is wrong")
    return tuple(out)


# ---------------------------------------------------------------------------
# Genome resolution: reuse exactly what build_pack_refs.py produced
# ---------------------------------------------------------------------------

class GenomeIndex:
    """Name -> reference genome, read out of the pack refs manifest.

    Resolution is by NAME, mirroring the watchlist: taxids rot silently when NCBI
    merges nodes, and a rotted taxid here would mean a spike quietly sourced from
    the wrong organism.
    """

    def __init__(self, manifest_path: Path, taxdump: Optional[Path] = None):
        self.manifest_path = manifest_path
        self.taxdump = taxdump
        self.by_name: Dict[str, List[dict]] = {}
        self.by_taxid: Dict[str, List[dict]] = {}
        self.by_accession: Dict[str, dict] = {}
        self._synonyms: Optional[Dict[str, str]] = None
        self._load()

    def _synonym_taxid(self, taxon: str) -> Optional[str]:
        """Resolve a name to a taxid through EVERY accepted NCBI name class.

        The manifest stores the CURRENT SCIENTIFIC name that build_pack_refs.py
        resolved, while scenario files and watchlists use the names people actually
        write. Those diverge constantly: "Vaccinia virus" is now
        "Orthopoxvirus vaccinia", "Francisella novicida" is
        "Francisella tularensis subsp. novicida", "Bundibugyo ebolavirus" is
        "Bundibugyo virus". An exact-name-only lookup declares all of them missing
        even though the genome is sitting on disk.

        Loaded lazily: it costs ~3M lines and is only needed when an exact match
        fails.
        """
        if self.taxdump is None:
            return None
        if self._synonyms is None:
            self._synonyms = {}
            names = Path(self.taxdump) / "names.dmp"
            if not names.is_file():
                LOG.warning("no names.dmp under %s; cannot resolve a scenario name that "
                            "differs from the manifest's scientific name", self.taxdump)
                return None
            with names.open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    f = [x.strip() for x in line.split("|")]
                    if len(f) < 2 or not f[0].isdigit():
                        continue
                    self._synonyms.setdefault(self._key(f[1]), f[0])
            LOG.info("loaded %d taxonomy name(s) for synonym resolution",
                     len(self._synonyms))
        return self._synonyms.get(self._key(taxon))

    @staticmethod
    def _key(name: str) -> str:
        return re.sub(r"[^a-z0-9]+", " ", str(name).lower()).strip()

    def _load(self) -> None:
        if not self.manifest_path.is_file():
            raise SystemExit(
                f"ERROR: reference manifest not found: {self.manifest_path}\n"
                f"       Spike genomes are taken from what build_pack_refs.py "
                f"already downloaded, never fetched here.\n"
                f"       Run: python3 scripts/build_pack_refs.py --pack <pack> "
                f"--db-root <db_root>"
            )
        with self.manifest_path.open(encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required = {"taxid", "name", "path"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise SystemExit(
                    f"ERROR: {self.manifest_path} is missing column(s): "
                    f"{', '.join(sorted(missing))}. Rebuild it with build_pack_refs.py."
                )
            rows = list(reader)
        usable = 0
        for row in rows:
            path = (row.get("path") or "").strip()
            if not path:
                # build_pack_refs.py records taxa it could not download with an
                # empty path. Keep them for the error message but never as a
                # genome source.
                continue
            usable += 1
            self.by_name.setdefault(self._key(row.get("name", "")), []).append(row)
            self.by_taxid.setdefault((row.get("taxid") or "").strip(), []).append(row)
            acc = (row.get("accession") or "").strip()
            if acc:
                self.by_accession.setdefault(acc, row)
        LOG.info("reference manifest: %d row(s), %d with a genome on disk (%s)",
                 len(rows), usable, self.manifest_path)
        if not usable:
            raise SystemExit(
                f"ERROR: {self.manifest_path} lists no genome with a path on disk.\n"
                f"       build_pack_refs.py has not successfully downloaded anything "
                f"yet; there is nothing to spike."
            )

    def resolve(self, taxon: str, explicit: str = "") -> dict:
        """Return {'name','taxid','path','genome_size','accession'} or raise.

        Never falls back to "something similar". A substituted genome would make
        the ground truth a lie, which is worse than a failed run.
        """
        if explicit:
            # An explicit path wins, and an explicit accession is looked up in the
            # same manifest so it still has to be a genome we actually hold.
            p = Path(explicit).expanduser()
            if p.is_file():
                return {"name": taxon, "taxid": "", "path": str(p),
                        "genome_size": "0", "accession": "", "source": "explicit-path"}
            if explicit in self.by_accession:
                row = dict(self.by_accession[explicit])
                row["source"] = "explicit-accession"
                return row
            raise SystemExit(
                f"ERROR: spike {taxon!r}: explicit genome {explicit!r} is neither an "
                f"existing file nor an accession in {self.manifest_path}.\n"
                f"       Fix the scenario file, or download the assembly with "
                f"build_pack_refs.py first. Refusing to substitute another genome."
            )

        key = self._key(taxon)
        hits = self.by_name.get(key)
        if not hits and str(taxon).strip().isdigit():
            hits = self.by_taxid.get(str(taxon).strip())
        resolved_via = "name"
        if not hits:
            # Synonym fallback: resolve the name to a taxid through the taxdump, then
            # match the manifest on taxid. This is matching the SAME ORGANISM under a
            # different accepted name, not guessing at a similar one.
            syn = self._synonym_taxid(taxon)
            if syn:
                hits = self.by_taxid.get(syn)
                if hits:
                    resolved_via = f"synonym -> taxid {syn}"
                    LOG.info("spike %r resolved via taxonomy synonym to %r (taxid %s)",
                             taxon, hits[0].get("name", "?"), syn)
        if not hits:
            # Genus-level fallback is deliberately NOT done. Suggest instead.
            suggestions = self._suggest(taxon)
            raise SystemExit(
                f"ERROR: spike taxon {taxon!r} could not be resolved to a reference "
                f"genome in {self.manifest_path}.\n"
                f"       This is fatal on purpose: a missing spike would read out as "
                f"a detection failure and corrupt every sensitivity number in the "
                f"batch.\n"
                + (f"       Did you mean: {', '.join(suggestions)}?\n" if suggestions else "")
                + "       Add the organism to the pack watchlist (or its "
                  "near_neighbour_taxids) and rerun build_pack_refs.py, or give the "
                  "spike an explicit `genome:` path."
            )
        # Prefer the row with the largest recorded genome: for a taxon with both a
        # chromosome-only and a chromosome+plasmid assembly, the complete one is
        # the honest choice, and exclude_loci exists to subtract from it.
        best = max(hits, key=lambda r: int(r.get("genome_size") or 0))
        row = dict(best)
        row["source"] = "manifest" if resolved_via == "name" else f"manifest ({resolved_via})"
        path = Path(row["path"])
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(
                f"ERROR: spike taxon {taxon!r} resolves to {path}, which is missing or "
                f"empty.\n"
                f"       The manifest and the refs directory have diverged; rerun "
                f"build_pack_refs.py (it re-verifies md5 and re-downloads)."
            )
        return row

    def _suggest(self, taxon: str, limit: int = 4) -> List[str]:
        """Cheap token-overlap suggestions - no third-party fuzzy matching."""
        want = set(self._key(taxon).split())
        if not want:
            return []
        scored = []
        for key, rows in self.by_name.items():
            overlap = len(want & set(key.split()))
            if overlap:
                scored.append((overlap, rows[0].get("name", key)))
        scored.sort(key=lambda t: (-t[0], t[1]))
        return [name for _, name in scored[:limit]]


def uncompressed_bases(path: Path) -> int:
    """Approximate base count of a FASTA, cheaply, without decompressing it.

    Needed because ART expresses depth as FOLD COVERAGE, so a wrong reference size
    silently produces a wrong number of reads. Getting it wrong is expensive: treating
    the 3.1 Gbp human genome as 1 bp produced -f 567600000 and would have run ART until
    the disk filled; treating it as its GZIPPED size (932 MB) overstates coverage 3.3x
    and quietly breaks the declared host fraction.

    For a .gz, gzip's trailing ISIZE field gives the uncompressed byte count in O(1).
    It is stored mod 2^32, so it is only trustworthy below 4 GiB - checked here, with a
    fall back to streaming the file when it is not. Byte count overstates the base count
    by roughly 1% for headers and newlines, which is immaterial for a coverage figure.
    """
    if str(path).endswith(".gz"):
        try:
            with path.open("rb") as fh:
                fh.seek(-4, os.SEEK_END)
                isize = struct.unpack("<I", fh.read(4))[0]
            # ISIZE wraps at 4 GiB. Accept it only when the compressed size makes a
            # sub-4-GiB original plausible (gzip rarely beats ~4x on FASTA).
            if isize > 0 and path.stat().st_size * 12 > isize:
                return isize
        except (OSError, struct.error):
            pass
        n = 0
        with gzip.open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 24), b""):
                n += len(chunk)
        return n
    try:
        return path.stat().st_size
    except OSError:
        return 0


def resolve_host_genome(db_root: Path, host: str) -> Optional[Path]:
    """Host genome as written by scripts/fetch_host_genome.sh.

    Globbed rather than hardcoded to chm13v2.0.fa.gz so that a re-hosted or
    updated assembly does not silently disable host simulation.
    """
    host_dir = db_root / "bioThreatFinder" / "host" / host
    if not host_dir.is_dir():
        return None
    for pattern in ("*.fa.gz", "*.fna.gz", "*.fasta.gz", "*.fa", "*.fna", "*.fasta"):
        for cand in sorted(host_dir.glob(pattern)):
            if cand.is_file() and cand.stat().st_size > 0:
                return cand
    return None


# ---------------------------------------------------------------------------
# Genome staging: decompress once, optionally remove loci
# ---------------------------------------------------------------------------

def locus_regex(locus: str) -> "re.Pattern":
    """Regex for a requested locus.

    A bare name is matched loosely (substring, case-insensitive) because header
    conventions vary; a value that looks like a regex is used as given so an
    operator can be precise when the alias table is not enough.
    """
    alias = LOCUS_ALIASES.get(locus.strip().lower())
    if alias:
        return re.compile(alias, re.IGNORECASE)
    if re.search(r"[\\^$.\[\]()|*+?{}]", locus):
        try:
            return re.compile(locus, re.IGNORECASE)
        except re.error as exc:
            raise SystemExit(
                f"ERROR: exclude_loci entry {locus!r} looks like a regex but does not "
                f"compile: {exc}"
            )
    return re.compile(re.escape(locus.strip()), re.IGNORECASE)


def open_fasta(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def staged_path(genome_cache: Path, label: str, genome_src: Path,
                exclude_loci: Sequence[str]) -> Path:
    """Where a staged plain FASTA lives in the shared cache.

    The source path is hashed into the name because the cache is shared between
    batches and packs: two spikes can carry the same organism label while pointing
    at different assemblies (one from the manifest, one an explicit `genome:`), and
    a collision there would silently simulate reads from the wrong genome.
    """
    suffix = ("_no_" + "_".join(sanitise(l) for l in exclude_loci)
              if exclude_loci else "")
    digest = hashlib.md5(str(genome_src).encode("utf-8")).hexdigest()[:8]
    return genome_cache / f"{sanitise(label)}{suffix}.{digest}.fna"


def stage_genome(src: Path, dest: Path, exclude_loci: Sequence[str],
                 taxon: str, force: bool) -> dict:
    """Write a plain (uncompressed) FASTA for the simulators, minus excluded loci.

    Plain FASTA because ART and NanoSim will not read gzip, and because handing
    every simulator the same input removes one variable from a comparison between
    them.

    Locus removal is how "B. anthracis chromosome WITHOUT pXO1/pXO2" is built. It
    is the case the marker gate exists for: chromosomal evidence looks exactly as
    strong as real anthrax, and the only correct answer is `supported`, not
    `confirmed_in_silico`. A silently-unmatched pattern would hand the marker gate
    the plasmids after all and the test would pass for the wrong reason - hence the
    hard error below.
    """
    meta_path = dest.with_suffix(dest.suffix + ".json")
    if (not force and dest.is_file() and dest.stat().st_size > 0
            and meta_path.is_file()):
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            LOG.info("  SKIP stage (exists): %s", dest.name)
            return meta
        except (OSError, ValueError):
            LOG.warning("  staged genome %s has an unreadable sidecar; restaging",
                        dest.name)

    patterns = [(loc, locus_regex(loc)) for loc in exclude_loci]
    matched: Dict[str, List[str]] = {loc: [] for loc, _ in patterns}
    excluded_headers: List[str] = []
    kept_headers: List[str] = []
    kept_bases = 0

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with open_fasta(src) as fh, tmp.open("w", encoding="utf-8") as out:
            keep = True
            for line in fh:
                if line.startswith(">"):
                    header = line[1:].rstrip()
                    hits = [loc for loc, rx in patterns if rx.search(header)]
                    keep = not hits
                    if hits:
                        for loc in hits:
                            matched[loc].append(header)
                        excluded_headers.append(header)
                    else:
                        kept_headers.append(header)
                    if keep:
                        out.write(line)
                elif keep:
                    out.write(line)
                    kept_bases += len(line.strip())
    except OSError as exc:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"ERROR: could not stage {src} -> {dest}: {exc}")

    unmatched = [loc for loc, hits in matched.items() if not hits]
    if unmatched:
        tmp.unlink(missing_ok=True)
        preview = "\n".join(f"         >{h}" for h in kept_headers[:12])
        raise SystemExit(
            f"ERROR: spike {taxon!r}: exclude_loci {', '.join(unmatched)} matched no "
            f"sequence in {src}.\n"
            f"       Refusing to continue: the sample would silently CONTAIN the locus "
            f"the scenario says it lacks, so a marker-gate test would pass for the "
            f"wrong reason.\n"
            f"       Headers available in this assembly:\n{preview}\n"
            f"       Fix the locus name, add an alias to LOCUS_ALIASES, or give the "
            f"spike an explicit `genome:` that has the locus named as expected."
        )
    if not kept_headers or kept_bases == 0:
        tmp.unlink(missing_ok=True)
        raise SystemExit(
            f"ERROR: spike {taxon!r}: exclude_loci "
            f"{', '.join(exclude_loci) or '(none)'} removed every sequence in {src}; "
            f"there is nothing left to simulate reads from."
        )
    tmp.replace(dest)

    meta = {
        "taxon": taxon,
        "source_fasta": str(src),
        "staged_fasta": str(dest),
        "exclude_loci": list(exclude_loci),
        "excluded_sequences": excluded_headers,
        "kept_sequences": len(kept_headers),
        "kept_bases": kept_bases,
        "md5": md5_of(dest),
        "staged_utc": utc_now(),
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    if exclude_loci:
        LOG.info("  staged %s minus %s: removed %d sequence(s) [%s], kept %d (%s)",
                 src.name, ",".join(exclude_loci), len(excluded_headers),
                 "; ".join(h.split()[0] for h in excluded_headers[:4]),
                 len(kept_headers), human(kept_bases))
    else:
        LOG.info("  staged %s: %d sequence(s), %s", src.name, len(kept_headers),
                 human(kept_bases))
    return meta


# ---------------------------------------------------------------------------
# Simulators
#
# Illumina, InSilicoSeq vs ART:
#   ISS models per-platform error and quality profiles empirically and takes an
#   abundance distribution over several genomes at once, which is exactly the
#   metagenomic mixture we need; it is slow and memory-hungry at 20M reads.
#   ART is much faster and is the long-standing validated choice for aligner
#   benchmarks, but its read count is expressed as fold coverage (or as reads PER
#   REFERENCE SEQUENCE), so exact library composition takes arithmetic on our side
#   and the realised count drifts. Realised counts are therefore always measured
#   from the FASTQ, never assumed - see count_fastq_reads().
#
# ONT, Badread vs NanoSim:
#   Badread needs no training data and deliberately generates the pathologies that
#   break metagenomic callers (chimeras, adapters, junk and random reads, glitches),
#   which is what a false-positive study wants. NanoSim is more faithful to a real
#   run's error profile but requires a trained model directory built from real
#   reads of the same chemistry; we do not ship one, so it is a fallback that must
#   be pointed at a model explicitly.
# ---------------------------------------------------------------------------

class Cmd:
    __slots__ = ("argv", "stdout_path")

    def __init__(self, argv: List[str], stdout_path: Optional[Path] = None):
        self.argv = [str(a) for a in argv]
        self.stdout_path = stdout_path

    def shell_repr(self) -> str:
        quoted = " ".join(
            a if re.fullmatch(r"[A-Za-z0-9_@%+=:,./-]+", a) else "'" + a.replace("'", "'\\''") + "'"
            for a in self.argv)
        return quoted + (f" > {self.stdout_path}" if self.stdout_path else "")


class Simulator:
    name = "abstract"
    platform = "abstract"
    paired = False

    def __init__(self, binary: str):
        self.binary = binary

    def commands(self, genome: Path, reads: int, seed: int, prefix: Path,
                 threads: int, opts: dict) -> List[Cmd]:
        raise NotImplementedError

    def output_globs(self, prefix: Path) -> Tuple[str, Optional[str]]:
        raise NotImplementedError


class IssSimulator(Simulator):
    name = "InSilicoSeq"
    platform = "illumina"
    paired = True

    @staticmethod
    def model_for(read_length: int) -> Tuple[str, int]:
        best = min(ISS_MODELS, key=lambda m: abs(m[0] - read_length))
        return best[1], best[0]

    def commands(self, genome, reads, seed, prefix, threads, opts):
        model, model_len = self.model_for(opts["read_length"])
        if model_len != opts["read_length"]:
            LOG.warning("  InSilicoSeq read length is fixed by its error model: "
                        "requested %d bp, model %s gives %d bp",
                        opts["read_length"], model, model_len)
        # VERIFY: --n_reads counts individual reads, so a paired run of N pairs is
        # requested as 2N. Confirm against `iss generate --help` for the installed
        # version before trusting the first batch's read counts; the manifest
        # records the measured count either way.
        return [Cmd([self.binary, "generate",
                     "--genomes", genome,
                     "--n_reads", max(2, reads * 2),
                     "--model", model,
                     "--abundance", "uniform",
                     "--cpus", threads,
                     "--seed", seed,
                     "--compress",
                     "--quiet",
                     "--output", prefix])]

    def output_globs(self, prefix):
        return (f"{prefix}_R1.fastq*", f"{prefix}_R2.fastq*")


class ArtSimulator(Simulator):
    name = "ART"
    platform = "illumina"
    paired = True

    def commands(self, genome, reads, seed, prefix, threads, opts):
        genome_bases = int(opts.get("genome_bases") or 0)
        if genome_bases <= 0:
            # Was `max(1, ...)`, which silently substituted 1 bp for "unknown" and
            # turned 1.89M host reads into -f 567600000 - 567 million-fold coverage of
            # the human genome. ART would have run until the disk filled. A missing
            # genome size is not a value to default; it is a bug to report.
            raise SystemExit(
                f"ERROR: ART needs the reference size to convert {reads} reads into fold "
                f"coverage, and genome_bases was not supplied for {genome}.\n"
                f"       Refusing to guess: defaulting it produces an absurd -f and ART "
                f"will run until the disk fills.")
        read_length = opts["read_length"]
        # art_illumina -c is reads PER REFERENCE SEQUENCE, which for a multi-contig
        # or plasmid-bearing assembly would over-generate wildly. Fold coverage is
        # the only knob that maps cleanly onto "N pairs from this whole genome".
        fold = (reads * 2.0 * read_length) / genome_bases
        if fold > 10000.0:
            raise SystemExit(
                f"ERROR: computed fold coverage {fold:,.0f}x for {genome} "
                f"({reads} reads over {genome_bases:,} bases) is implausible.\n"
                f"       Refusing to run: ART would generate until the disk fills. This "
                f"almost always means the reference size is wrong, not that the depth "
                f"was requested.")
        fold = max(fold, 1e-6)
        # VERIFY: -ss profile name (HS25 = HiSeq 2500) and the -m/-s fragment
        # parameters against `art_illumina --help`; -na suppresses the SAM output.
        return [Cmd([self.binary,
                     "-ss", "HS25",
                     "-i", genome,
                     "-p",
                     "-l", read_length,
                     "-f", f"{fold:.6f}",
                     "-m", max(read_length * 2 + 50, 300),
                     "-s", 50,
                     "-rs", seed,
                     "-na",
                     "-o", f"{prefix}_R"])]

    def output_globs(self, prefix):
        return (f"{prefix}_R1.fq*", f"{prefix}_R2.fq*")


class BadreadSimulator(Simulator):
    name = "Badread"
    platform = "ont"
    paired = False

    def commands(self, genome, reads, seed, prefix, threads, opts):
        mean_len = max(200, int(opts["ont_mean_length"]))
        bases = max(mean_len, reads * mean_len)
        # VERIFY: --quantity accepts a bare base count as well as "50x"/"250M";
        # --length takes MEAN,STDEV. Both per badread's documented CLI.
        return [Cmd([self.binary, "simulate",
                     "--reference", genome,
                     "--quantity", str(bases),
                     "--length", f"{mean_len},{max(100, mean_len // 2)}",
                     "--seed", seed],
                    stdout_path=Path(f"{prefix}.fastq"))]

    def output_globs(self, prefix):
        return (f"{prefix}.fastq", None)


class NanosimSimulator(Simulator):
    name = "NanoSim"
    platform = "ont"
    paired = False

    def __init__(self, binary: str, model_prefix: Optional[str]):
        super().__init__(binary)
        self.model_prefix = model_prefix

    def commands(self, genome, reads, seed, prefix, threads, opts):
        if not self.model_prefix:
            raise SystemExit(
                "ERROR: NanoSim needs a trained model and none was given.\n"
                "       Set NANOSIM_MODEL_PREFIX to a pre-trained model prefix (the\n"
                "       path stem of the *_model_profile files), or install the\n"
                "       preferred simulator instead:\n"
                f"       {CONDA_BIN}/mamba install -c bioconda badread\n"
                "       Badread needs no training data, which is why it is preferred.")
        # VERIFY: NanoSim's entry point has been `simulator.py genome` and, in
        # recent releases, `nanosim genome`; flags -rg/-c/-n/-o/-t and --fastq
        # against the installed version.
        return [Cmd([self.binary, "genome",
                     "-rg", genome,
                     "-c", self.model_prefix,
                     "-n", reads,
                     "-o", prefix,
                     "-t", threads,
                     "--seed", seed,
                     "--fastq"])]

    def output_globs(self, prefix):
        # genome mode writes <prefix>_aligned_reads.fastq (+ _unaligned_reads.fastq
        # when a chimera/unaligned fraction is modelled); take whatever appeared.
        return (f"{prefix}*_reads.fastq", None)


def resolve_simulators(platforms: Sequence[str], illumina_choice: str = "auto",
                       ont_choice: str = "auto") -> Dict[str, Optional[Simulator]]:
    """Pick the simulator per platform, and say how to install it.

    THROUGHPUT, measured on this pack rather than assumed: InSilicoSeq runs at roughly
    350k read-pairs/min on a single sample and parallelises poorly, which for the full
    359-sample matrix is ~124 h sequential. ART is C rather than Python and is
    dramatically faster, at the cost of expressing depth as FOLD COVERAGE rather than a
    read count - so requested counts are approximate. That is tolerable here because the
    manifest records the MEASURED count (`spiked_reads`), which is the denominator a
    sensitivity number should use anyway.

    'auto' keeps InSilicoSeq's more realistic error model when it is present; pass 'art'
    when wall time is the binding constraint.
    """
    out: Dict[str, Optional[Simulator]] = {}
    if "illumina" in platforms:
        iss, art = which("iss"), which("art_illumina")
        if illumina_choice == "iss" and not iss:
            raise SystemExit("ERROR: --illumina-simulator iss requested but `iss` is not "
                             "on PATH.\n       mamba install -n <env> -c bioconda insilicoseq")
        if illumina_choice == "art" and not art:
            raise SystemExit("ERROR: --illumina-simulator art requested but `art_illumina` "
                             "is not on PATH.\n       mamba install -n <env> -c bioconda art")
        if illumina_choice == "art":
            out["illumina"] = ArtSimulator(art)
            LOG.info("illumina simulator: ART (%s) - chosen explicitly. Depth is "
                     "expressed as fold coverage, so requested read counts are "
                     "approximate; the manifest records the measured count.", art)
        elif iss and illumina_choice in ("auto", "iss"):
            out["illumina"] = IssSimulator(iss)
            LOG.info("illumina simulator: InSilicoSeq (%s)", iss)
        elif art:
            out["illumina"] = ArtSimulator(art)
            LOG.warning("illumina simulator: InSilicoSeq not found, falling back to "
                        "ART (%s). Read counts become approximate (ART expresses "
                        "depth as fold coverage); measured counts are still recorded.",
                        art)
        else:
            out["illumina"] = None
            LOG.error("no Illumina simulator found. Install the preferred one:\n"
                      "    %s/mamba install -c bioconda insilicoseq\n"
                      "  or the fallback:\n"
                      "    %s/mamba install -c bioconda art", CONDA_BIN, CONDA_BIN)
    if "ont" in platforms:
        badread = which("badread")
        nanosim = which("simulator.py") or which("nanosim")
        if ont_choice == "nanosim" and not nanosim:
            raise SystemExit("ERROR: --ont-simulator nanosim requested but it is not on "
                             "PATH.\n       mamba install -n <env> -c bioconda nanosim")
        if ont_choice == "badread" and not badread:
            raise SystemExit("ERROR: --ont-simulator badread requested but it is not on "
                             "PATH.\n       mamba install -n <env> -c bioconda badread")
        if ont_choice == "nanosim":
            out["ont"] = NanosimSimulator(nanosim, os.environ.get("NANOSIM_MODEL_PREFIX"))
            LOG.info("ont simulator: NanoSim (%s) - chosen explicitly", nanosim)
        elif badread and ont_choice in ("auto", "badread"):
            out["ont"] = BadreadSimulator(badread)
            LOG.info("ont simulator: Badread (%s)", badread)
        elif nanosim:
            out["ont"] = NanosimSimulator(nanosim, os.environ.get("NANOSIM_MODEL_PREFIX"))
            LOG.warning("ont simulator: Badread not found, falling back to NanoSim "
                        "(%s). NanoSim needs a trained model; set "
                        "NANOSIM_MODEL_PREFIX.", nanosim)
        else:
            out["ont"] = None
            LOG.error("no ONT simulator found. Install the preferred one:\n"
                      "    %s/mamba install -c bioconda badread\n"
                      "  or the fallback:\n"
                      "    %s/mamba install -c bioconda nanosim", CONDA_BIN, CONDA_BIN)
    return out


# ---------------------------------------------------------------------------
# FASTQ plumbing
# ---------------------------------------------------------------------------

def count_fastq_reads(path: Path) -> int:
    """Records in a FASTQ, counted rather than assumed.

    Every simulator misses its requested count to some degree (ART works in fold
    coverage, Badread in bases, ISS rounds across genomes). The manifest's
    spiked_reads is the denominator of a sensitivity figure, so it has to be the
    number of reads that actually exist.
    """
    lines = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    # Counted in binary chunks rather than by iterating lines: a 20M-read library
    # is 80M lines, and the per-line Python loop makes the count cost as much as
    # the simulation.
    with opener(path, "rb") as fh:  # type: ignore[operator]
        while True:
            chunk = fh.read(1 << 22)
            if not chunk:
                break
            lines += chunk.count(b"\n")
    if lines % 4:
        LOG.warning("%s has %d lines, not a multiple of 4 - truncated FASTQ?",
                    path, lines)
    return lines // 4


def concat_to_gz(sources: Sequence[Path], dest: Path) -> None:
    """Concatenate FASTQ parts (plain or gzipped) into one .fastq.gz.

    Gzip members concatenate byte-wise into a valid stream, so an already-gzipped
    part is copied verbatim instead of being decompressed and recompressed. Plain
    parts get their own member. Every reader we care about (zcat, seqkit,
    fastp, minimap2) handles multi-member gzip.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with tmp.open("wb") as out:
        if not sources:
            # A genuinely empty library (a blank NTC with no background panel) is a
            # legitimate benchmark row, but it must still be a VALID empty gzip -
            # a zero-byte .fastq.gz would crash the pipeline instead of testing it.
            with gzip.GzipFile(fileobj=out, mode="wb", compresslevel=6, mtime=0):
                pass
        for src in sources:
            if str(src).endswith(".gz"):
                with src.open("rb") as fh:
                    shutil.copyfileobj(fh, out, length=1 << 20)
            else:
                with gzip.GzipFile(fileobj=out, mode="wb", compresslevel=6, mtime=0) as gz:
                    with src.open("rb") as fh:
                        shutil.copyfileobj(fh, gz, length=1 << 20)
    tmp.replace(dest)


def append_selected_reads(src: Path, dest: Path, indices: Sequence[int]) -> int:
    """Append the selected read records of src to an existing .fastq.gz.

    Written as a new gzip member (see concat_to_gz). Records are copied unmodified:
    a hopped read is a real read that carries the wrong barcode, so altering it
    would make the test easier than reality.
    """
    wanted = set(indices)
    if not wanted:
        return 0
    written = 0
    opener = gzip.open if str(src).endswith(".gz") else open
    with dest.open("ab") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6, mtime=0) as out:
            with opener(src, "rb") as fh:  # type: ignore[operator]
                idx = 0
                while True:
                    record = [fh.readline() for _ in range(4)]
                    if not record[0]:
                        break
                    if idx in wanted:
                        out.write(b"".join(record))
                        written += 1
                    idx += 1
    return written


# ---------------------------------------------------------------------------
# Sample planning
# ---------------------------------------------------------------------------

class Component:
    """One organism's contribution to one sample."""
    __slots__ = ("kind", "label", "genome_src", "staged", "reads", "seed",
                 "exclude_loci", "genome_bases", "meta")

    def __init__(self, kind: str, label: str, genome_src: Path, reads: int,
                 seed: int, exclude_loci: Sequence[str] = (),
                 genome_bases: int = 0):
        self.kind = kind                # host | spike | background
        self.label = label
        self.genome_src = genome_src
        self.staged: Optional[Path] = None
        self.reads = reads
        self.seed = seed
        self.exclude_loci = list(exclude_loci)
        self.genome_bases = genome_bases
        self.meta: dict = {}

    @property
    def tag(self) -> str:
        suffix = ("_no_" + "_".join(sanitise(l) for l in self.exclude_loci)
                  if self.exclude_loci else "")
        return f"{self.kind}_{sanitise(self.label)}{suffix}"


class SamplePlan:
    __slots__ = ("sample", "scenario", "platform", "replicate", "level_rpm",
                 "role", "seed", "components", "fastq_1", "fastq_2",
                 "host_fraction", "total_reads", "notes", "excluded_loci",
                 "truth_present", "truth_absent", "truth_vaccine",
                 "spiked_reads", "actual_total", "status", "failure")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))
        self.notes = list(kw.get("notes") or [])
        self.spiked_reads = 0
        self.actual_total = 0
        self.status = "planned"
        self.failure = ""

    def manifest_row(self) -> dict:
        return {
            "sample": self.sample,
            "scenario_id": self.scenario.sid,
            "kind": self.scenario.kind,
            "replicate": str(self.replicate),
            "platform": self.platform,
            "level_rpm": "" if self.level_rpm is None else ("%g" % self.level_rpm),
            "role": self.role,
            "truth_present": ";".join(self.truth_present),
            "truth_absent": ";".join(self.truth_absent),
            "truth_vaccine": ";".join(self.truth_vaccine),
            "host_fraction": "%.6g" % self.host_fraction,
            "total_reads": str(self.total_reads),
            "spiked_reads": str(self.spiked_reads),
            "excluded_loci": ";".join(self.excluded_loci),
            "fastq_1": str(self.fastq_1),
            "fastq_2": str(self.fastq_2 or ""),
            "seed": str(self.seed),
            "notes": ";".join(n for n in self.notes if n),
        }

    def samplesheet_row(self) -> dict:
        return {
            "sample": self.sample,
            "platform": self.platform,
            "fastq_1": str(self.fastq_1),
            "fastq_2": str(self.fastq_2 or ""),
            "role": self.role,
        }


def plan_batch(scenarios: Sequence[Scenario], genomes: GenomeIndex,
               db_root: Path, reads_dir: Path, master_seed: int,
               platform_filter: Optional[List[str]],
               replicates_override: Optional[int]) -> Tuple[List[SamplePlan], List[str]]:
    """Expand scenarios into concrete samples. Resolution failures are collected
    and reported together, then fatal: seeing all of them at once is the
    difference between one rerun and eight."""
    plans: List[SamplePlan] = []
    problems: List[str] = []
    host_cache: Dict[str, Optional[Path]] = {}

    for sc in scenarios:
        platforms = [p for p in sc.platforms
                     if not platform_filter or p in platform_filter]
        if not platforms:
            continue
        reps = replicates_override or sc.replicates

        if sc.host not in host_cache:
            host_cache[sc.host] = resolve_host_genome(db_root, sc.host)
        host_fasta = host_cache[sc.host]
        if sc.host_fraction > 0 and host_fasta is None:
            problems.append(
                f"{sc.sid}: host_fraction={sc.host_fraction} but no host genome under "
                f"{db_root / 'bioThreatFinder' / 'host' / sc.host}. Run "
                f"scripts/fetch_host_genome.sh {sc.host}, or set host_fraction: 0.")
            continue

        # Resolve every genome first; a scenario is either fully buildable or an
        # error, never partially spiked.
        try:
            spike_refs = [genomes.resolve(sp.taxon, sp.genome) for sp in sc.spikes]
            bg_refs = [(name, frac, resolve_background(name, frac, genomes))
                       for name, frac in sc.background.items()]
        except SystemExit as exc:
            problems.append(f"{sc.sid}: {exc}")
            continue

        primary = sc.primary_spike
        levels = primary.levels_rpm if primary else [None]
        if sc.kind != "lod_series":
            levels = levels[:1]

        # Staged in a local list so that a scenario which fails halfway through its
        # replicates contributes nothing: half a dose series is not a limit of
        # detection measurement, it is a misleading one.
        mine: List[SamplePlan] = []
        try:
            # None means "the scenario's own depth"; a list adds a depth axis.
            depth_axis = sc.downsample_total_reads or [None]
            for platform in platforms:
                for rep in range(1, reps + 1):
                    for level in levels:
                        for depth in depth_axis:
                            mine.append(build_sample_plan(
                                sc, platform, rep, level, spike_refs, bg_refs,
                                host_fasta, reads_dir, master_seed, primary,
                                role="sample", depth_override=depth))
                    if sc.include_ntc:
                        mine.append(build_ntc_plan(
                            sc, platform, rep, bg_refs, reads_dir, master_seed))
        except SystemExit as exc:
            problems.append(f"{sc.sid}: {exc}")
            continue
        plans.extend(mine)
    return plans, problems


def resolve_background(name: str, frac: float, genomes: GenomeIndex) -> dict:
    """A background key is either a path to a FASTA or an organism in the manifest.

    Background is not ground truth for the watchlist, but it is still not allowed
    to be silently dropped: a scenario whose background did not materialise has a
    different host/background ratio than its manifest claims, and every rpm figure
    computed from it is wrong.
    """
    p = Path(str(name)).expanduser()
    if p.is_file():
        return {"name": p.stem, "path": str(p), "taxid": "", "genome_size": "0",
                "source": "background-file"}
    return genomes.resolve(name)


def _compose(sc: Scenario, level: Optional[float],
             spike_refs: List[dict], bg_refs: List[Tuple[str, float, dict]],
             host_fasta: Optional[Path], total_reads: int,
             primary: Optional[Spike], master_seed: int,
             sample: str) -> Tuple[List[Component], float, List[str]]:
    """Turn fractions into per-component read counts.

    Rule: spikes and background take their declared share; the host absorbs the
    remainder. That is what a real clinical specimen looks like (overwhelmingly
    host, with the pathogen as a trace), and it keeps the realised library size
    equal to total_reads so that an rpm level means what it says.
    """
    notes: List[str] = []
    components: List[Component] = []
    claimed = 0.0

    # strict=True on purpose: spike_refs is built as a comprehension over sc.spikes so
    # the lengths match today, but a silent zip truncation here would DROP A SPIKE while
    # still writing its truth label - the exact failure this script refuses elsewhere.
    for sp, ref in zip(sc.spikes, spike_refs, strict=True):
        rpm = level if (primary is not None and sp is primary and level is not None) \
            else sp.levels_rpm[0]
        frac = rpm / 1e6
        reads = int(round(total_reads * frac))
        if reads <= 0:
            # Honesty guard: a level so low it rounds to nothing is not a limit of
            # detection measurement, it is an empty sample mislabelled as a spike.
            raise SystemExit(
                f"ERROR: {sc.sid}/{sample}: spike {sp.taxon!r} at {rpm:g} rpm of "
                f"{total_reads} reads rounds to 0 reads.\n"
                f"       Raise total_reads (>= {int(1e6 / max(rpm, 1e-9)):,} needed "
                f"for one read at this level) or raise the level.")
        claimed += frac
        components.append(Component(
            kind="spike", label=sp.taxon, genome_src=Path(ref["path"]),
            reads=reads,
            seed=derive_seed(master_seed, sample, "spike", sp.taxon,
                             ",".join(sp.exclude_loci), "%g" % rpm),
            exclude_loci=sp.exclude_loci,
            genome_bases=int(ref.get("genome_size") or 0)))

    for name, frac, ref in bg_refs:
        reads = int(round(total_reads * frac))
        if reads <= 0:
            LOG.warning("%s: background %s at fraction %g rounds to 0 reads; dropped",
                        sample, name, frac)
            notes.append(f"background_dropped={sanitise(name)}")
            continue
        claimed += frac
        components.append(Component(
            kind="background", label=ref.get("name") or name,
            genome_src=Path(ref["path"]), reads=reads,
            seed=derive_seed(master_seed, sample, "background", name),
            genome_bases=int(ref.get("genome_size") or 0)))

    if claimed > 1.0 + 1e-9:
        raise SystemExit(
            f"ERROR: {sc.sid}: spikes + background claim {claimed:.4f} of the library "
            f"(> 1.0). Lower the levels_rpm or the background fractions.")

    host_frac = sc.host_fraction
    remainder = 1.0 - claimed
    if host_fasta is not None and host_frac > 0:
        if host_frac > remainder + 1e-9:
            raise SystemExit(
                f"ERROR: {sc.sid}: host_fraction {host_frac:g} plus spikes/background "
                f"{claimed:.4f} exceeds the library. Lower host_fraction to "
                f"<= {remainder:.4f}.")
        if remainder - host_frac > 1e-9:
            # Host is the filler, so the realised host fraction is recorded rather
            # than the requested one - the manifest must describe the FASTQ.
            notes.append(f"host_filled={remainder - host_frac:.6g}")
        host_frac = remainder
        # genome_bases is required by the ART path to compute fold coverage. At this
        # point host_fasta is still the gzipped reference, so its file size is NOT the
        # base count - see uncompressed_bases().
        host_bases = uncompressed_bases(host_fasta)
        components.insert(0, Component(
            kind="host", label=sc.host, genome_src=host_fasta,
            reads=int(round(total_reads * host_frac)),
            seed=derive_seed(master_seed, sample, "host", sc.host),
            genome_bases=host_bases))
    else:
        host_frac = 0.0
        if remainder > 1e-9:
            has_spike = any(c.kind == "spike" for c in components)
            if components and has_spike:
                # THE OLD BEHAVIOUR HERE MISREPORTED CONCENTRATION BY ORDERS OF
                # MAGNITUDE, AND SAID THE OPPOSITE IN A COMMENT.
                #
                # It scaled every component by 1/claimed and noted "rpm keeps its
                # meaning". It does not. With no host there is nothing to dilute with,
                # so the spike's share of the DELIVERED library is spike_frac/claimed,
                # not spike_frac - scaling changes the library's size, not that ratio.
                # For a 100 rpm spike with 0.1% background and host_fraction 0,
                # claimed is 0.0011 and the realised concentration is ~90,900 rpm while
                # the manifest's level_rpm column still said 100. Every LoD point from
                # such a scenario would be wrong by ~3 orders of magnitude, in the
                # direction that flatters sensitivity, and evaluate_benchmark keys its
                # curves on level_rpm.
                #
                # Refused rather than recorded, because there is no honest number to
                # put in level_rpm: the operator asked for a concentration the mixture
                # cannot express.
                #
                # Only when a SPIKE is present. A no-template control and a negative
                # scenario are background-only BY DESIGN (build_ntc_plan sets
                # host_fraction to 0 deliberately) and carry level_rpm = None, so there
                # is no declared concentration to misreport and filling the library to
                # depth is correct. Gating on spikes rather than on host_fraction is
                # what separates the two cases.
                raise SystemExit(
                    f"ERROR: {sc.sid}: spikes and background claim only "
                    f"{claimed:.6f} of the library and there is no host to make up "
                    f"the remaining {remainder:.6f}.\n"
                    f"       The delivered library would then be {claimed:.6f} of the "
                    f"requested depth, so every component's concentration would be "
                    f"{1.0 / max(claimed, 1e-12):.4g}x what was declared - and "
                    f"level_rpm would report the declared value.\n"
                    f"       Fix by one of:\n"
                    f"         - set host_fraction (the host absorbs the remainder), or\n"
                    f"         - add background organisms summing to "
                    f"{remainder:.6f}, or\n"
                    f"         - raise levels_rpm so the spikes fill the library.")
            elif components:
                # Spike-free: ratios among background organisms are preserved and no
                # declared concentration is affected, so fill the library to depth.
                scale = 1.0 / max(claimed, 1e-12)
                for comp in components:
                    comp.reads = max(1, int(round(comp.reads * scale)))
                notes.append(f"scaled_to_fill={scale:.6g}")
            else:
                notes.append("empty_library")
    return components, host_frac, notes


def build_sample_plan(sc: Scenario, platform: str, rep: int, level: Optional[float],
                      spike_refs: List[dict], bg_refs, host_fasta: Optional[Path],
                      reads_dir: Path, master_seed: int, primary: Optional[Spike],
                      role: str, depth_override: Optional[int] = None) -> SamplePlan:
    parts = [sc.sid, platform, f"r{rep}"]
    if sc.kind == "lod_series" and level is not None:
        parts.append(f"rpm{level_tag(level)}")
    # The depth belongs in the NAME, or the six down-sampled libraries of one
    # scenario would all resolve to one filename and overwrite each other.
    if depth_override is not None:
        parts.append(f"n{depth_override}")
    sample = "_".join(parts)
    seed = derive_seed(master_seed, sc.sid, rep, platform,
                       "" if level is None else "%g" % level, role,
                       "" if depth_override is None else str(depth_override))
    # Per-platform depth. An ONT run at a 5 kb mean does not deliver an Illumina
    # read count, so using one number for both would simulate a flow cell that does
    # not exist and overstate long-read sensitivity.
    depth = sc.ont_total_reads if platform == "ont" else sc.total_reads
    depth_note = None
    if depth_override is not None:
        # downsample_total_reads is written in ILLUMINA read counts. Applying them
        # verbatim to ONT would simulate 40M reads at a 5 kb mean - 200 Gbp, which no
        # flow cell delivers - and would overstate long-read sensitivity for the same
        # reason ont_total_reads exists at all. Scale by the platform's own depth ratio
        # so the titration keeps its SHAPE at depths the platform can actually reach,
        # and record the scaling rather than leaving the manifest to imply otherwise.
        if platform == "ont" and sc.total_reads:
            ratio = sc.ont_total_reads / float(sc.total_reads)
            depth = max(1, int(round(depth_override * ratio)))
            depth_note = (f"downsample={depth_override}|ont_scaled_by={ratio:g}"
                          f"|realised={depth}")
        else:
            depth = depth_override
            depth_note = f"downsample={depth_override}"
    components, host_frac, notes = _compose(
        sc, level, spike_refs, bg_refs, host_fasta, depth,
        primary, master_seed, sample)

    if depth_note:
        notes.append(depth_note)

    excluded = []
    for sp in sc.spikes:
        for loc in sp.exclude_loci:
            excluded.append(f"{sp.taxon}:{loc}")

    # truth.present defaults to the spiked organisms: writing it out twice in the
    # scenario file is an opportunity to disagree with the mixture.
    present = sc.present or [sp.taxon for sp in sc.spikes]
    return SamplePlan(
        sample=sample, scenario=sc, platform=platform, replicate=rep,
        level_rpm=level if sc.kind != "negative" else None, role=role, seed=seed,
        components=components,
        fastq_1=reads_dir / (f"{sample}_R1.fastq.gz" if platform == "illumina"
                             else f"{sample}.fastq.gz"),
        fastq_2=(reads_dir / f"{sample}_R2.fastq.gz") if platform == "illumina" else None,
        host_fraction=host_frac, total_reads=depth,
        truth_present=present, truth_absent=sc.absent, truth_vaccine=sc.vaccine,
        excluded_loci=excluded, notes=notes)


def build_ntc_plan(sc: Scenario, platform: str, rep: int, bg_refs,
                   reads_dir: Path, master_seed: int) -> SamplePlan:
    """A no-template control for this scenario: background only, no host, no spike.

    The NTC is what control subtraction is computed against, so it must be built
    from the same background panel as its samples - a pristine NTC would make the
    `suppressed_control` tier untestable.
    """
    sample = f"{sc.sid}_{platform}_r{rep}_NTC"
    # PER-PLATFORM DEPTH, like build_sample_plan. This used sc.total_reads on both
    # platforms, so with the defaults (20M Illumina / 4M ONT) the ONT NTC was planned
    # at 200,000 reads against ONT samples of 4M - a control five times more
    # concentrated relative to its samples than intended, which is the wrong
    # denominator for control subtraction. At a 5 kb mean it is also ~1 Gbp of
    # simulated ONT reads per NTC, for nothing.
    sample_depth = sc.ont_total_reads if platform == "ont" else sc.total_reads
    total = max(NTC_MIN_READS, int(sample_depth * NTC_READ_FRACTION))
    seed = derive_seed(master_seed, sc.sid, rep, platform, "NTC")
    ntc_sc = Scenario(**{k: getattr(sc, k) for k in Scenario.__slots__})
    ntc_sc.spikes = []
    ntc_sc.host_fraction = 0.0
    ntc_sc.total_reads = total
    components, host_frac, notes = _compose(
        ntc_sc, None, [], bg_refs, None, total, None, master_seed, sample)
    notes.append("ntc")
    # The spiked targets belong in truth_absent even if index hopping later puts a
    # handful of their reads here: "present in trace amounts via barcode bleed" is
    # exactly the thing the pipeline must NOT call, so absent is the correct
    # expectation and the hop is recorded in notes.
    absent = list(dict.fromkeys(list(sc.absent) + [sp.taxon for sp in sc.spikes]))
    return SamplePlan(
        sample=sample, scenario=ntc_sc, platform=platform, replicate=rep,
        level_rpm=None, role="control", seed=seed, components=components,
        fastq_1=reads_dir / (f"{sample}_R1.fastq.gz" if platform == "illumina"
                             else f"{sample}.fastq.gz"),
        fastq_2=(reads_dir / f"{sample}_R2.fastq.gz") if platform == "illumina" else None,
        host_fraction=host_frac, total_reads=total,
        truth_present=[c.label for c in components if c.kind == "background"],
        truth_absent=absent, truth_vaccine=[], excluded_loci=[], notes=notes)


# ---------------------------------------------------------------------------
# Size estimate
# ---------------------------------------------------------------------------

def estimate_bytes(plans: Sequence[SamplePlan], genome_bytes: int) -> Tuple[int, int]:
    """(estimated output bytes, total reads).

    Per FASTQ record: header + sequence + '+' + quality, i.e. roughly
    2 * read_length + FASTQ_OVERHEAD_BYTES uncompressed, times GZIP_RATIO for
    gzip -6 on simulated data, times ESTIMATE_SAFETY because each simulator writes
    uncompressed per-component parts before we concatenate. Paired Illumina reads
    are counted once each (total_reads is the count of reads, not of pairs).
    Staged plain-FASTA genomes are added separately - the human host genome alone
    is ~3 GB decompressed.
    """
    total_reads = 0
    raw = 0
    for plan in plans:
        sc = plan.scenario
        length = sc.read_length if plan.platform == "illumina" else sc.ont_mean_length
        total_reads += plan.total_reads
        raw += plan.total_reads * (length + FASTQ_OVERHEAD_BYTES)
    return int(raw * GZIP_RATIO * ESTIMATE_SAFETY) + genome_bytes, total_reads


def estimate_genome_bytes(plans: Sequence[SamplePlan]) -> int:
    """Decompressed size of the distinct genomes we will stage (~3x the .gz)."""
    seen: Dict[str, int] = {}
    for plan in plans:
        for comp in plan.components:
            key = str(comp.genome_src) + "|" + ",".join(comp.exclude_loci)
            if key in seen:
                continue
            try:
                size = comp.genome_src.stat().st_size
            except OSError:
                size = 0
            seen[key] = size * 3 if str(comp.genome_src).endswith(".gz") else size
    return sum(seen.values())


# ---------------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------------

def run_cmd(cmd: Cmd, log_fh, cwd: Optional[Path] = None) -> bool:
    log_fh.write(f"# {utc_now()}\n{cmd.shell_repr()}\n")
    log_fh.flush()
    LOG.debug("running: %s", cmd.shell_repr())
    try:
        if cmd.stdout_path:
            cmd.stdout_path.parent.mkdir(parents=True, exist_ok=True)
            with cmd.stdout_path.open("wb") as out:
                proc = subprocess.run(cmd.argv, stdout=out, stderr=subprocess.PIPE,
                                      cwd=str(cwd) if cwd else None)
        else:
            proc = subprocess.run(cmd.argv, capture_output=True, cwd=str(cwd) if cwd else None)
    except OSError as exc:
        LOG.error("  could not run %s: %s", cmd.argv[0], exc)
        log_fh.write(f"# FAILED to launch: {exc}\n")
        return False
    stderr = (proc.stderr or b"").decode("utf-8", "replace").strip()
    if proc.returncode != 0:
        LOG.error("  %s exited %d: %s", Path(cmd.argv[0]).name, proc.returncode,
                  "; ".join(stderr.splitlines()[-3:]))
        log_fh.write(f"# exit {proc.returncode}\n{stderr}\n")
        return False
    if stderr:
        log_fh.write(f"# stderr:\n{stderr}\n")
    return True


def resolve_glob(pattern: Optional[str]) -> List[Path]:
    if not pattern:
        return []
    return [Path(p) for p in sorted(glob.glob(pattern)) if Path(p).is_file()]


def generate_sample(plan: SamplePlan, sim: Simulator, work_dir: Path,
                    genome_cache: Path, logs_dir: Path, threads: int,
                    force: bool, staged_seen: set) -> bool:
    """Build one sample. Returns True on success; sets plan.failure on failure."""
    sidecar = logs_dir / f"{plan.sample}.json"
    outputs = [plan.fastq_1] + ([plan.fastq_2] if plan.fastq_2 else [])

    if not force and all(p.is_file() and p.stat().st_size > 0 for p in outputs):
        if sidecar.is_file():
            try:
                data = json.loads(sidecar.read_text(encoding="utf-8"))

                # THE READS ON DISK MUST MATCH THE PLAN THAT WAS JUST COMPUTED.
                #
                # The skip only required the FASTQs to exist and the sidecar to parse.
                # It never compared what the sidecar recorded against the plan. So
                # editing a scenario - bumping host_fraction, changing a level, a new
                # seed - and re-running the same batch reused the old FASTQs for every
                # sample whose NAME was unchanged, while the manifest was rewritten
                # from the NEW plan. The manifest then described libraries that were
                # never generated, and the evaluator computed an LoD curve against
                # declared levels that did not match the reads. Ground truth that
                # silently does not describe the data is worse than no ground truth.
                mismatches = []
                if data.get("seed") is not None and int(data["seed"]) != int(plan.seed):
                    mismatches.append(f"seed {data['seed']} != {plan.seed}")
                old_lvl, new_lvl = data.get("level_rpm"), plan.level_rpm
                if (old_lvl is None) != (new_lvl is None) or (
                        old_lvl is not None and new_lvl is not None
                        and abs(float(old_lvl) - float(new_lvl)) > 1e-9):
                    mismatches.append(f"level_rpm {old_lvl} != {new_lvl}")
                old_hf = data.get("host_fraction")
                if old_hf is not None and abs(float(old_hf)
                                              - float(plan.host_fraction)) > 1e-6:
                    mismatches.append(
                        f"host_fraction {old_hf} != {plan.host_fraction}")
                old_tr = data.get("total_reads_intended")
                if old_tr is not None and int(old_tr) != int(plan.total_reads):
                    mismatches.append(f"total_reads {old_tr} != {plan.total_reads}")
                if mismatches:
                    LOG.warning("REGENERATE %s: existing reads do not match the "
                                "current plan (%s). Rebuilding rather than shipping a "
                                "manifest that describes data that was never made.",
                                plan.sample, "; ".join(mismatches))
                    for pth in outputs:
                        pth.unlink(missing_ok=True)
                    sidecar.unlink(missing_ok=True)
                else:
                    plan.spiked_reads = int(data.get("spiked_reads", 0))
                    plan.actual_total = int(data.get("actual_total_reads", 0))
                    plan.status = "skipped"
                    plan.notes.append("resumed")
                    # Re-add the notes a fresh run would have written, so a resumed
                    # batch's manifest is byte-comparable with the original.
                    if plan.actual_total != plan.total_reads:
                        plan.notes.append(
                            f"actual_total_reads={plan.actual_total}")
                    LOG.info("SKIP %s (outputs present, counts from %s)",
                             plan.sample, sidecar.name)
                    return True
            except (OSError, ValueError, TypeError):
                LOG.warning("%s: outputs exist but %s is unreadable", plan.sample,
                            sidecar.name)
        # Reads exist but the per-component counts do not. The FASTQ cannot be
        # re-derived into truth (spike reads are indistinguishable once merged), so
        # a dataset without its sidecar has no usable ground truth and is rebuilt.
        LOG.warning("%s: FASTQ present but per-component counts are missing; "
                    "regenerating (its ground truth is unknown otherwise)",
                    plan.sample)

    sample_work = work_dir / plan.sample
    if sample_work.exists():
        shutil.rmtree(sample_work, ignore_errors=True)
    sample_work.mkdir(parents=True, exist_ok=True)

    cmd_log = logs_dir / f"{plan.sample}.cmds"
    parts_1: List[Path] = []
    parts_2: List[Path] = []
    spiked = 0
    total = 0
    comp_records: List[dict] = []

    with cmd_log.open("w", encoding="utf-8") as log_fh:
        log_fh.write(f"# sample: {plan.sample}\n# scenario: {plan.scenario.sid}\n"
                     f"# platform: {plan.platform}\n# sample seed: {plan.seed}\n"
                     f"# simulator: {sim.name} ({sim.binary})\n")
        for comp in plan.components:
            staged = staged_path(genome_cache, comp.label, comp.genome_src,
                                 comp.exclude_loci)
            LOG.info("  %s: %s -> %d read(s), seed %d", plan.sample, comp.tag,
                     comp.reads, comp.seed)
            # --force must not re-stage the same genome once per sample: the
            # decompressed human host genome is ~3 GB, and the cache is shared.
            restage = force and str(staged) not in staged_seen
            comp.meta = stage_genome(comp.genome_src, staged, comp.exclude_loci,
                                     comp.label, restage)
            staged_seen.add(str(staged))
            comp.staged = staged
            if not comp.genome_bases:
                comp.genome_bases = int(comp.meta.get("kept_bases") or 0)

            prefix = sample_work / comp.tag
            opts = {"read_length": plan.scenario.read_length,
                    "ont_mean_length": plan.scenario.ont_mean_length,
                    "genome_bases": comp.meta.get("kept_bases") or comp.genome_bases}
            for cmd in sim.commands(staged, comp.reads, comp.seed, prefix,
                                    threads, opts):
                if not run_cmd(cmd, log_fh):
                    plan.failure = f"{sim.name} failed on component {comp.tag}"
                    return False

            g1, g2 = sim.output_globs(str(prefix))
            f1, f2 = resolve_glob(g1), resolve_glob(g2)
            if not f1:
                plan.failure = (f"{sim.name} produced no reads for {comp.tag} "
                                f"(expected {g1})")
                LOG.error("  %s", plan.failure)
                return False
            if sim.paired and not f2:
                plan.failure = (f"{sim.name} produced no mate file for {comp.tag} "
                                f"(expected {g2})")
                LOG.error("  %s", plan.failure)
                return False
            got = sum(count_fastq_reads(p) for p in f1)
            parts_1.extend(f1)
            parts_2.extend(f2)
            total += got
            if comp.kind == "spike":
                spiked += got
            if got != comp.reads:
                LOG.info("  %s: %s asked for %d read(s), got %d (recorded as-is)",
                         plan.sample, comp.tag, comp.reads, got)
            comp_records.append({
                "kind": comp.kind, "label": comp.label,
                "genome_src": str(comp.genome_src), "staged_fasta": str(staged),
                "exclude_loci": comp.exclude_loci,
                "excluded_sequences": comp.meta.get("excluded_sequences", []),
                "seed": comp.seed, "reads_requested": comp.reads,
                "reads_generated": got,
            })

        concat_to_gz(parts_1, plan.fastq_1)
        if plan.fastq_2:
            concat_to_gz(parts_2, plan.fastq_2)
        log_fh.write(f"# wrote {plan.fastq_1}"
                     + (f" and {plan.fastq_2}\n" if plan.fastq_2 else "\n"))

    # Intermediates are large; the command log plus the seeds are what makes the
    # sample reproducible, so the parts themselves are disposable.
    shutil.rmtree(sample_work, ignore_errors=True)

    plan.spiked_reads = spiked
    plan.actual_total = total
    plan.status = "generated"
    if total != plan.total_reads:
        plan.notes.append(f"actual_total_reads={total}")
    sidecar.write_text(json.dumps({
        "sample": plan.sample, "scenario_id": plan.scenario.sid,
        "platform": plan.platform, "replicate": plan.replicate,
        "role": plan.role, "seed": plan.seed,
        "level_rpm": plan.level_rpm, "host_fraction": plan.host_fraction,
        "total_reads_intended": plan.total_reads,
        "actual_total_reads": total, "spiked_reads": spiked,
        "simulator": {"name": sim.name, "binary": sim.binary},
        "components": comp_records,
        "fastq_1": str(plan.fastq_1), "fastq_2": str(plan.fastq_2 or ""),
        "generated_utc": utc_now(),
    }, indent=2) + "\n", encoding="utf-8")
    return True


# ---------------------------------------------------------------------------
# Index hopping
#
# Why this is simulated at all: on a patterned-flowcell Illumina run, free adapter
# lets a fragment from one library acquire another library's index, at rates around
# 0.1-2%. In low-biomass metagenomics that is routinely the dominant source of
# false positives - a high-titre positive next to an NTC will put real target reads
# into the NTC, and a caller with no control-subtraction step will confirm them.
# `suppressed_control` exists for exactly this, and it cannot be validated on clean
# data. The hopped reads are recorded in notes but deliberately NOT added to
# spiked_reads: they are contamination, not a deliberate spike, and the expected
# behaviour for the recipient remains "do not call it".
# ---------------------------------------------------------------------------

def spike_read_ranges(logs_dir: Optional[Path], plan: SamplePlan) -> List[Tuple[int, int]]:
    """Half-open [start, end) read-index ranges of the SPIKE components in plan's FASTQ.

    generate_sample() concatenates the per-component FASTQs in component order and
    records reads_generated for each, in that same order, so the cumulative sum locates
    any component's reads inside the merged file.

    WHY THIS IS NEEDED. The hop used `rng.sample(range(donor.actual_total), n)` - a
    uniform draw over the donor's WHOLE library, host and background included. The
    number of TARGET reads moved was therefore n x spike_fraction, not n. With the
    shipped human_core.yaml (host_fraction 0.94, background 4.9%, 20M reads) a 100 rpm
    donor is 1e-4 spike, an NTC recipient plans 200,000 reads, and at
    index_hop_fraction 0.001 that is n=200 reads carrying an EXPECTED 0.02 target
    reads. The NTC received 200 host and background reads, `notes` asserted a hop had
    happened, and the contamination scenarios were effectively clean with respect to
    the target - so a pipeline with no control-subtraction step at all would pass them.
    The module comment above says `suppressed_control` "cannot be validated on clean
    data", and that is exactly what the data was.
    """
    if not logs_dir:
        return []
    sidecar = logs_dir / f"{plan.sample}.json"
    if not sidecar.is_file():
        return []
    try:
        data = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    ranges: List[Tuple[int, int]] = []
    offset = 0
    for comp in data.get("components") or []:
        try:
            count = int(comp.get("reads_generated") or 0)
        except (TypeError, ValueError):
            count = 0
        if count <= 0:
            continue
        if comp.get("kind") == "spike":
            ranges.append((offset, offset + count))
        offset += count
    return ranges


def apply_index_hopping(plans: Sequence[SamplePlan], master_seed: int,
                        logs_dir: Optional[Path], dry_run: bool) -> int:
    groups: Dict[Tuple[str, str, int], List[SamplePlan]] = {}
    for plan in plans:
        if plan.scenario.index_hop_fraction > 0 and plan.status in ("generated", "skipped"):
            groups.setdefault((plan.scenario.sid, plan.platform, plan.replicate),
                              []).append(plan)
    hops = 0
    for (sid, platform, rep), members in sorted(groups.items()):
        frac = members[0].scenario.index_hop_fraction
        donors = [p for p in members if p.spiked_reads > 0]
        if not donors:
            LOG.warning("%s/%s/r%d: index_hop_fraction=%g but no spiked sample to "
                        "bleed from; skipped", sid, platform, rep, frac)
            continue
        # Bleed from the most concentrated library, which is where hopping is worst
        # and which is the worst case the pipeline has to survive.
        donor = max(donors, key=lambda p: p.spiked_reads)
        for recipient in members:
            if recipient is donor:
                continue
            n = int(round(recipient.total_reads * frac))
            if n <= 0:
                continue
            if dry_run:
                LOG.info("[dry-run] would bleed %d read(s) (%g) from %s into %s",
                         n, frac, donor.sample, recipient.sample)
                hops += 1
                continue

            # Appending reads is not idempotent, so the sidecar records which hops
            # have already been applied. Without this, a resumed run would bleed a
            # second dose into every recipient and quietly change the dataset that
            # earlier numbers were measured on.
            sidecar = (logs_dir / f"{recipient.sample}.json") if logs_dir else None
            data: dict = {}
            if sidecar and sidecar.is_file():
                try:
                    data = json.loads(sidecar.read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    data = {}
            done = data.get("index_hops") or []
            prior = next((h for h in done if h.get("from") == donor.sample), None)
            if prior:
                recipient.notes.append(
                    f"index_hop={frac:g}|from={donor.sample}|reads={prior.get('reads', 0)}")
                recipient.actual_total += int(prior.get("reads", 0))
                LOG.info("  SKIP index hop %s -> %s (already applied)",
                         donor.sample, recipient.sample)
                continue

            available = donor.actual_total or count_fastq_reads(donor.fastq_1)
            if available <= 0:
                LOG.warning("  donor %s has no reads; no hop into %s",
                            donor.sample, recipient.sample)
                continue

            # Draw from the donor's SPIKE reads, so the recipient really does receive
            # target sequence. This makes the simulation a deliberate WORST CASE rather
            # than a faithful uniform-bleed model: on a real flow cell the hopped reads
            # would be dominated by host, and at these spike fractions essentially no
            # target read would move - which is physically accurate and tests nothing.
            # The tier being validated (`suppressed_control`) is only exercised when
            # target reads are present in the control, so that is what is injected, and
            # the note records it as target-derived rather than implying a uniform bleed.
            spike_ranges = spike_read_ranges(logs_dir, donor)
            eligible = [i for a, b in spike_ranges for i in range(a, b)]
            if eligible:
                n = min(n, len(eligible))
                pool, drawn_from = eligible, "spike"
            else:
                # No per-component record (an older batch, or a donor whose sidecar is
                # gone): fall back to the whole library and SAY SO, because the hop
                # then probably carries no target reads and the scenario is not
                # measuring what it claims.
                LOG.warning("  donor %s has no per-component read counts; hopping "
                            "uniformly from the whole library, which at its spike "
                            "fraction will carry few or no target reads - the "
                            "contamination scenario will not exercise "
                            "suppressed_control. Regenerate the donor with --force.",
                            donor.sample)
                n = min(n, available)
                pool, drawn_from = list(range(available)), "library"
            rng = random.Random(derive_seed(master_seed, donor.sample,
                                           recipient.sample, "index_hop"))
            picks = sorted(rng.sample(pool, n))
            moved = append_selected_reads(donor.fastq_1, recipient.fastq_1, picks)
            if recipient.fastq_2 and donor.fastq_2:
                # Same record indices in both mates: the pair hops together,
                # because it is one physical fragment with one mis-assigned index.
                append_selected_reads(donor.fastq_2, recipient.fastq_2, picks)
            recipient.notes.append(
                f"index_hop={frac:g}|from={donor.sample}|reads={moved}"
                f"|drawn_from={drawn_from}")
            recipient.actual_total += moved
            if sidecar:
                done.append({"from": donor.sample, "fraction": frac, "reads": moved,
                             "drawn_from": drawn_from, "applied_utc": utc_now()})
                # actual_total_reads is deliberately NOT updated: it means "reads
                # this script generated for this sample", and hopped reads are
                # contamination from another library, tracked in index_hops.
                data["index_hops"] = done
                sidecar.write_text(json.dumps(data, indent=2) + "\n",
                                   encoding="utf-8")
            LOG.info("  index hop: %d read(s) %s -> %s (%g)", moved, donor.sample,
                     recipient.sample, frac)
            hops += 1
    return hops


# ---------------------------------------------------------------------------
# Output writing (atomic, after every sample)
# ---------------------------------------------------------------------------

def carry_forward(path: Path, produced: set, key: str, delimiter: str) -> list:
    """Rows of an existing manifest/samplesheet for samples NOT produced this run.

    WHY THIS EXISTS. Both writers rewrote the whole file from the current
    invocation's plans, and `plans` is built from the --scenario/--platform filtered
    scenario list. So generating a 359-sample batch and then re-running
    `--batch v1 --scenario lod_banthracis` to regenerate one series rewrote
    manifest.tsv with that scenario's ~30 rows. Every other sample's FASTQ was still
    in reads/ with its sidecar intact - only its TRUTH row was gone.

    That is worse than an obviously incomplete file, because evaluate_benchmark.py
    iterates the manifest and scores one observation per row: a sample missing from
    the manifest is not "untested" and not an error, it silently leaves the
    denominator. Sensitivity would then be computed over 30 samples and reported as
    the batch result.
    """
    if not path.exists():
        return []
    try:
        with path.open(encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter=delimiter))
    except OSError as exc:
        raise SystemExit(
            f"ERROR: cannot read {path} to merge ({exc}). Refusing to rewrite it from "
            f"this run alone: the other samples' reads would stay on disk while their "
            f"ground-truth rows vanished, and the evaluator would silently score a "
            f"smaller batch and report it as the batch result.") from exc
    return [r for r in rows if (r.get(key) or "").strip() not in produced]


def write_manifest(path: Path, plans: Sequence[SamplePlan]) -> None:
    fresh = [p for p in plans if p.status in ("generated", "skipped")]
    produced = {p.sample for p in fresh}
    carried = carry_forward(path, produced, "sample", "\t")
    tmp = path.with_suffix(".tsv.part")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for plan in fresh:
            w.writerow(plan.manifest_row())
        for row in carried:
            w.writerow({c: row.get(c, "") for c in MANIFEST_COLUMNS})
    tmp.replace(path)


def write_samplesheet(path: Path, plans: Sequence[SamplePlan]) -> None:
    fresh = [p for p in plans if p.status in ("generated", "skipped")]
    produced = {p.sample for p in fresh}
    carried = carry_forward(path, produced, "sample", ",")
    tmp = path.with_suffix(".csv.part")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SAMPLESHEET_COLUMNS, lineterminator="\n")
        w.writeheader()
        for plan in fresh:
            w.writerow(plan.samplesheet_row())
        for row in carried:
            w.writerow({c: row.get(c, "") for c in SAMPLESHEET_COLUMNS})
    tmp.replace(path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(
        prog="simulate_benchmark.py",
        description="Generate synthetic metagenomes with known ground truth so that "
                    "bioThreatFinder's sensitivity, specificity and limit of "
                    "detection can be measured per taxon.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  simulate_benchmark.py --scenarios benchmarks/scenarios.yaml --batch v1 --dry-run
  simulate_benchmark.py --scenarios benchmarks/scenarios.yaml --batch v1 --threads 32
  simulate_benchmark.py --scenarios benchmarks/scenarios.yaml --batch v1 \\
      --scenario banth_lod --scenario banth_ablation --platform illumina

notes:
  * Spike genomes come from what scripts/build_pack_refs.py already downloaded.
    A taxon that cannot be resolved there is a HARD ERROR: a missing spike would
    read out as a detection failure and corrupt the measurement.
  * Deterministic. Every sample's seed is derived from --seed plus its scenario,
    replicate, platform and level, and is recorded in the manifest, so adding a
    scenario does not change any existing sample.
  * Idempotent and resumable: existing non-empty FASTQ with its counts sidecar is
    skipped; --force rebuilds. The manifest is rewritten after every sample.
""")
    p.add_argument("--scenarios", required=True, type=Path,
                   help="scenario definitions (YAML, parsed by bin/btf_yaml.py)")
    p.add_argument("--batch", required=True,
                   help="batch name; everything lands under <outdir>/<batch>/")
    p.add_argument("--outdir", type=Path, default=None,
                   help="output root (default: <db-root>/bioThreatFinder/simulated). "
                        "Simulated FASTQ is large - keep it off the repo volume.")
    p.add_argument("--taxdump", type=Path,
                   default=Path("/media/8tb/bioThreatFinder/taxdump/2026-09-29"),
                   help="pinned taxonomy snapshot, used ONLY to resolve a scenario "
                        "taxon name that differs from the manifest's current "
                        "scientific name (default: %(default)s)")
    p.add_argument("--db-root", type=Path, default=DEFAULT_DB_ROOT,
                   help="shared database root, as used by build_pack_refs.py "
                        "(default: %(default)s)")
    p.add_argument("--packs-dir", type=Path, default=repo_root / "packs",
                   help="directory containing pack definitions (default: %(default)s)")
    p.add_argument("--pack", default="human",
                   help="pack whose reference genomes are spiked (default: %(default)s)")
    p.add_argument("--scenario", action="append", default=[], metavar="ID",
                   help="only build this scenario id (repeatable)")
    p.add_argument("--platform", choices=PLATFORMS + ["both"], default="both",
                   help="restrict platforms (default: %(default)s)")
    p.add_argument("--replicates", type=int, default=None,
                   help="override the replicate count from the scenario file")
    p.add_argument("--seed", type=int, default=1,
                   help="master seed; per-sample seeds are derived from it "
                        "(default: %(default)s)")
    p.add_argument("--illumina-simulator", choices=("auto", "iss", "art"), default="auto",
                   help="auto prefers InSilicoSeq for its error model; 'art' is far "
                        "faster (C, not Python) but expresses depth as fold coverage so "
                        "read counts are approximate - the manifest records the measured "
                        "count either way (default: %(default)s)")
    p.add_argument("--ont-simulator", choices=("auto", "badread", "nanosim"),
                   default="auto",
                   help="auto prefers Badread, which needs no training data "
                        "(default: %(default)s)")
    p.add_argument("--threads", type=int, default=8,
                   help="threads passed to the simulators (default: %(default)s)")
    p.add_argument("--dry-run", action="store_true",
                   help="print the plan and every command, generate nothing")
    p.add_argument("--force", action="store_true",
                   help="regenerate samples whose outputs already exist")
    p.add_argument("--min-free-gib", type=int, default=DEFAULT_MIN_FREE_GIB,
                   help="headroom required on top of the estimated output size "
                        "(default: %(default)s)")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%SZ",
        stream=sys.stderr,
    )
    logging.Formatter.converter = time.gmtime

    if args.batch != sanitise(args.batch):
        raise SystemExit(f"ERROR: --batch {args.batch!r} is not filesystem-safe; "
                         f"use letters, digits, '.', '_', '-'")
    if args.replicates is not None and args.replicates <= 0:
        raise SystemExit("ERROR: --replicates must be positive")
    if args.threads <= 0:
        raise SystemExit("ERROR: --threads must be positive")

    outdir = args.outdir or (args.db_root / "bioThreatFinder" / "simulated")
    batch_dir = outdir / args.batch
    reads_dir = batch_dir / "reads"
    logs_dir = batch_dir / "logs"
    work_dir = batch_dir / "work"
    genome_cache = outdir / "_genome_cache"
    manifest_path = batch_dir / "manifest.tsv"
    samplesheet_path = batch_dir / "samplesheet.csv"

    platform_filter = None if args.platform == "both" else [args.platform]

    LOG.info("batch=%s outdir=%s pack=%s seed=%d threads=%d",
             args.batch, outdir, args.pack, args.seed, args.threads)
    LOG.info("scenarios: %s", args.scenarios)

    pack_dir = args.packs_dir / args.pack
    if not pack_dir.is_dir():
        raise SystemExit(
            f"ERROR: pack {args.pack!r} not found at {pack_dir}\n"
            f"       Available: "
            f"{', '.join(sorted(d.name for d in args.packs_dir.glob('*') if d.is_dir())) or '(none)'}")
    refs_manifest = (args.db_root / "bioThreatFinder" / "packs" / args.pack /
                     "refs" / "manifest.tsv")

    scenarios = load_scenarios(args.scenarios)
    LOG.info("scenario file: %d scenario(s)", len(scenarios))
    if args.scenario:
        wanted = set(args.scenario)
        unknown = wanted - {sc.sid for sc in scenarios}
        if unknown:
            raise SystemExit(
                f"ERROR: --scenario {', '.join(sorted(unknown))} not in "
                f"{args.scenarios}\n"
                f"       Defined: {', '.join(sc.sid for sc in scenarios)}")
        scenarios = tuple(sc for sc in scenarios if sc.sid in wanted)

    genomes = GenomeIndex(refs_manifest, taxdump=args.taxdump)
    plans, problems = plan_batch(scenarios, genomes, args.db_root, reads_dir,
                                 args.seed, platform_filter, args.replicates)
    if problems:
        LOG.error("--- scenarios that cannot be built ---")
        for msg in problems:
            for line in str(msg).splitlines():
                LOG.error("  %s", line)
        raise SystemExit(
            "ERROR: refusing to generate a partial batch. Every scenario must "
            "resolve before any sample is written, because a batch with a silently "
            "missing spike produces sensitivity numbers that are wrong rather than "
            "absent.")
    if not plans:
        raise SystemExit("ERROR: nothing to do - no scenario matched the "
                         "--scenario/--platform filters")

    genome_bytes = estimate_genome_bytes(plans)
    est_bytes, est_reads = estimate_bytes(plans, genome_bytes)

    # ---- plan report ----------------------------------------------------
    LOG.info("plan: %d sample(s) across %d scenario(s)", len(plans), len(scenarios))
    for sc in scenarios:
        mine = [p for p in plans if p.scenario.sid == sc.sid]
        LOG.info("  %-24s kind=%-15s samples=%-4d platforms=%s%s", sc.sid, sc.kind,
                 len(mine),
                 ",".join(sorted({p.platform for p in mine})),
                 "  " + sc.description if sc.description else "")
    LOG.info("plan: %s intended reads total; estimated output %s "
             "(%s of that is staged genomes)",
             f"{est_reads:,}", human(est_bytes), human(genome_bytes))
    LOG.info("plan: reads -> %s", reads_dir)
    LOG.info("plan: manifest -> %s", manifest_path)
    LOG.info("plan: samplesheet -> %s", samplesheet_path)

    sims = resolve_simulators(sorted({p.platform for p in plans}),
                              illumina_choice=args.illumina_simulator,
                              ont_choice=args.ont_simulator)
    missing = [plat for plat, sim in sims.items() if sim is None]
    if missing and not args.dry_run:
        raise SystemExit(
            f"ERROR: no simulator available for: {', '.join(missing)}. "
            f"See the install commands above.")

    # ---- dry run ---------------------------------------------------------
    if args.dry_run:
        LOG.info("[dry-run] would create %s, %s, %s", reads_dir, logs_dir, genome_cache)
        for plan in plans:
            sim = sims.get(plan.platform)
            LOG.info("[dry-run] %s  role=%s seed=%d level_rpm=%s host_fraction=%.4g "
                     "total_reads=%d", plan.sample, plan.role, plan.seed,
                     "-" if plan.level_rpm is None else "%g" % plan.level_rpm,
                     plan.host_fraction, plan.total_reads)
            LOG.info("[dry-run]   truth present=[%s] absent=[%s] vaccine=[%s] "
                     "excluded_loci=[%s]",
                     ";".join(plan.truth_present), ";".join(plan.truth_absent),
                     ";".join(plan.truth_vaccine), ";".join(plan.excluded_loci))
            LOG.info("[dry-run]   fastq_1=%s fastq_2=%s", plan.fastq_1,
                     plan.fastq_2 or "")
            for comp in plan.components:
                staged = staged_path(genome_cache, comp.label, comp.genome_src,
                                     comp.exclude_loci)
                LOG.info("[dry-run]   %-40s %10d read(s) seed=%d <- %s",
                         comp.tag, comp.reads, comp.seed, comp.genome_src)
                if comp.exclude_loci:
                    LOG.info("[dry-run]     would drop sequences matching: %s",
                             ", ".join(comp.exclude_loci))
                if sim is None:
                    LOG.info("[dry-run]     (no simulator installed for %s)",
                             plan.platform)
                    continue
                opts = {"read_length": plan.scenario.read_length,
                        "ont_mean_length": plan.scenario.ont_mean_length,
                        "genome_bases": comp.genome_bases}
                try:
                    cmds = sim.commands(staged, comp.reads, comp.seed,
                                        Path(f"/dry-run/{plan.sample}/{comp.tag}"),
                                        args.threads, opts)
                except SystemExit as exc:
                    LOG.info("[dry-run]     %s", exc)
                    continue
                for cmd in cmds:
                    LOG.info("[dry-run]     %s", cmd.shell_repr())
        apply_index_hopping(plans, args.seed, None, dry_run=True)
        LOG.info("[dry-run] free-space check would require %s "
                 "(estimate %s + %d GiB headroom)",
                 human(est_bytes + args.min_free_gib * (1024 ** 3)),
                 human(est_bytes), args.min_free_gib)
        LOG.info("[dry-run] nothing generated")
        return 0

    # ---- for real --------------------------------------------------------
    need_gib = est_bytes / (1024 ** 3) + args.min_free_gib
    check_free_space(outdir, need_gib)

    for d in (reads_dir, logs_dir, work_dir, genome_cache):
        d.mkdir(parents=True, exist_ok=True)

    # simulate.log is part of the deliverable: a batch nobody can trace back to a
    # command is not reproducible.
    handler = logging.FileHandler(batch_dir / "simulate.log", encoding="utf-8")
    handler.setFormatter(logging.Formatter(
        "[%(asctime)s] %(levelname)s %(message)s", datefmt="%Y-%m-%dT%H:%M:%SZ"))
    logging.getLogger().addHandler(handler)
    LOG.info("run started %s: %s", utc_now(), " ".join(sys.argv))

    failed: List[SamplePlan] = []
    staged_seen: set = set()
    for i, plan in enumerate(plans, start=1):
        sim = sims[plan.platform]
        LOG.info("[%d/%d] %s (%s, %s, seed %d)", i, len(plans), plan.sample,
                 plan.scenario.kind, plan.platform, plan.seed)
        try:
            ok = generate_sample(plan, sim, work_dir, genome_cache, logs_dir,
                                 args.threads, args.force, staged_seen)
        except SystemExit:
            raise
        except Exception as exc:  # a broken sample must not lose the batch
            ok = False
            plan.failure = f"{type(exc).__name__}: {exc}"
            LOG.error("  %s failed: %s", plan.sample, plan.failure)
        if not ok:
            plan.status = "failed"
            failed.append(plan)
            for p in (plan.fastq_1, plan.fastq_2):
                # Never leave a truncated FASTQ that a later run would "resume".
                if p and p.exists():
                    p.unlink(missing_ok=True)
        # Rewrite after every sample so an interrupted run still has a valid,
        # self-consistent manifest for what it did finish.
        write_manifest(manifest_path, plans)
        write_samplesheet(samplesheet_path, plans)

    hops = apply_index_hopping(plans, args.seed, logs_dir, dry_run=False)
    write_manifest(manifest_path, plans)
    write_samplesheet(samplesheet_path, plans)
    shutil.rmtree(work_dir, ignore_errors=True)

    # ---- summary ---------------------------------------------------------
    built = [p for p in plans if p.status in ("generated", "skipped")]
    print("", file=sys.stderr)
    LOG.info("=========================== SUMMARY ===========================")
    LOG.info("batch:            %s", args.batch)
    LOG.info("batch dir:        %s", batch_dir)
    LOG.info("manifest:         %s (%d row(s))", manifest_path, len(built))
    LOG.info("samplesheet:      %s", samplesheet_path)
    LOG.info("samples built:    %d generated, %d skipped, %d FAILED",
             sum(1 for p in plans if p.status == "generated"),
             sum(1 for p in plans if p.status == "skipped"), len(failed))
    LOG.info("index hops applied: %d", hops)
    LOG.info("--- samples per scenario / platform ---")
    for sc in scenarios:
        for platform in PLATFORMS:
            mine = [p for p in built if p.scenario.sid == sc.sid
                    and p.platform == platform]
            if not mine:
                continue
            LOG.info("  %-24s %-9s %3d sample(s), %15s read(s), %12s spike read(s)",
                     sc.sid, platform, len(mine),
                     f"{sum(p.actual_total for p in mine):,}",
                     f"{sum(p.spiked_reads for p in mine):,}")
    LOG.info("total reads written:  %s", f"{sum(p.actual_total for p in built):,}")
    LOG.info("total spike reads:    %s", f"{sum(p.spiked_reads for p in built):,}")
    try:
        used = sum(f.stat().st_size for f in batch_dir.rglob("*") if f.is_file())
        LOG.info("disk used by batch:   %s (estimated %s)", human(used),
                 human(est_bytes))
    except OSError:
        pass
    try:
        cached = sum(f.stat().st_size for f in genome_cache.rglob("*") if f.is_file())
        LOG.info("staged genome cache:  %s (%s, shared between batches)",
                 human(cached), genome_cache)
    except OSError:
        pass
    LOG.info("free space remaining: %s", human(shutil.disk_usage(outdir).free))

    if failed:
        LOG.warning("--- samples that could not be built ---")
        for plan in failed:
            LOG.warning("  %s (%s/%s): %s", plan.sample, plan.scenario.sid,
                        plan.platform, plan.failure or "unknown")
        LOG.warning("These are absent from the manifest rather than present with a "
                    "guessed row: a benchmark must not contain a sample whose truth "
                    "it cannot state. Rerun to retry - finished samples are skipped.")
    LOG.info("===============================================================")

    # Non-zero if anything failed: a partially built benchmark batch must not pass
    # silently into a measurement.
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        LOG.warning("interrupted - rerun to resume (completed samples are skipped)")
        sys.exit(130)
