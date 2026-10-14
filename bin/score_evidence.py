#!/usr/bin/env python3
"""
Turn alignment metrics into tiered findings.

Design principle: every finding explains itself. The ``reasons`` column records
which gate decided the tier, so an analyst can see why something was promoted or
held back without re-running anything. A tool that emits a bare verdict on a
high-consequence organism is not usable in a public-health workflow.

Tiers (ascending)
-----------------
no_reference            Watchlist taxon had no reference genome in the pack.
                        NOT a negative - a database gap. Must be fixed and rerun.
screened_only           Stage-1 lead that failed stage-2. The normal outcome for
                        k-mer noise. Reported for audit, not for action.
attributed_to_neighbour Reads align better to a near neighbour. Usually the
                        correct answer for environmental relatives - but only
                        where the target has no claim of its own, which is what
                        --contest-attribution checks before letting it stand.
indistinguishable_complex
                        The target covers its own reference AND so does a near
                        neighbour that out-counts it only slightly. One of the
                        two IS present and this evidence does not say which.
                        Neither a detection nor a negative; needs a lineage
                        marker or a targeted assay to resolve.
suppressed_control      Present in a negative control at a comparable level.
                        Reagent/index contamination until proven otherwise.
control_background      What the negative control itself contains. Reported for
                        QC, never as a finding about a specimen.
marker_veto             A marker that argues AGAINST the target was found (a
                        neighbour-specific gene, or a region the target lacks).
vaccine_strain          Marker evidence indicates a vaccine or attenuated strain
                        rather than wild type. Actionable as a stand-down, so it
                        is reported prominently rather than suppressed.
supported               Real alignment evidence, below confirmation bar.
confirmed_in_silico     Meets every threshold. Still requires wet-lab confirmation.

MARKER GATING

For many taxa the chromosome is not the actionable signal, so alignment metrics
alone cannot reach `confirmed_in_silico`. Where packs/<pack>/marker_rules.yaml
defines a rule, it gates promotion under three-valued logic:

  confirm_requires TRUE     -> promotion allowed
  confirm_requires FALSE    -> capped at `supported`, with the failing rule named
  confirm_requires UNKNOWN  -> capped at `supported`, because the markers were
                               never tested. A database gap must not read as a
                               clean negative.

This is why a chromosomally convincing B. anthracis signal with no pXO1/pXO2
evidence stops at `supported`: without plasmid markers the call is not
defensible, whatever the breadth says.

Nothing here is a diagnosis.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

# marker_rules lives beside this script in bin/. Imported rather than
# reimplemented so the three-valued semantics have exactly one definition.
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    import marker_rules as mr
except ImportError:  # pragma: no cover - marker gating simply disabled
    mr = None

# Tiers that assert SOMETHING IS THERE, as opposed to explaining why nothing was
# called. A control row reaching any of them is control background, not a finding.
# Kept as a named set because the comment at the reclassification said the test was
# "did anything assert presence, not an enumeration" while the code was an inline
# enumeration - so adding a presence-asserting tier silently bypassed it, which is
# the exact failure that comment was written about.
ASSERTS_PRESENCE = {
    "supported", "confirmed_in_silico", "vaccine_strain", "marker_veto",
    "indistinguishable_complex",
}

TIER_ORDER = [
    "no_reference",
    "control_background",
    "screened_only",
    "attributed_to_neighbour",
    # Below `supported` because it asserts no single organism, above
    # `attributed_to_neighbour` because it asserts that one of a named pair IS
    # present - which a silent attribution does not.
    "indistinguishable_complex",
    "marker_veto",
    "suppressed_control",
    "vaccine_strain",
    "supported",
    "confirmed_in_silico",
]


def read_tsv(path: Path) -> list[dict]:
    if not path or not path.exists() or path.name == "NO_FILE":
        return []
    with path.open() as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    if not lines:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def apply_threshold_overrides(thresholds: dict, specs: list) -> dict:
    """Apply KEY=VALUE overrides onto the pack thresholds, loudly.

    Anything unparseable, or any key the pack does not define, is FATAL: a silently
    ignored override is how -profile test came to assert relaxed gates that were never
    applied. A typo must not look like a threshold that took effect.
    """
    out = dict(thresholds)
    for spec in specs or []:
        if "=" not in spec:
            sys.exit(f"[score] FATAL: --threshold expects KEY=VALUE, got {spec!r}")
        key, _, raw = spec.partition("=")
        key, raw = key.strip(), raw.strip()
        if key not in out:
            sys.exit(f"[score] FATAL: --threshold {key!r} is not a known threshold. "
                     f"Known: {', '.join(sorted(k for k in out))}")
        if isinstance(out[key], dict):
            sys.exit(f"[score] FATAL: --threshold {key!r} is a nested mapping and "
                     f"cannot be overridden from the command line; edit the pack.")
        try:
            out[key] = float(raw)
        except ValueError:
            sys.exit(f"[score] FATAL: --threshold {key}={raw!r} is not a number")
        sys.stderr.write(f"[score] threshold override: {key} = {out[key]} "
                         f"(pack.yaml said {thresholds[key]})\n")
    return out


TYPING_CLASSES = ("pathotype", "lineage", "vaccine_strain", "resistance")


def expected_marker_reads(mean_depth: float, marker_bp: float,
                          read_len: float) -> float:
    """Reads expected on a single-copy locus of `marker_bp` at `mean_depth`.

    depth = reads * read_len / length, so reads = depth * length / read_len. Using the
    CHROMOSOMAL depth for a locus that may sit on a multi-copy plasmid underestimates
    the expectation, which only ever makes the gate more conservative - it errs towards
    UNKNOWN rather than towards a confident absence.
    """
    if marker_bp <= 0 or read_len <= 0:
        return 0.0
    return mean_depth * marker_bp / read_len


def gate_absence_on_depth(states: dict, metas: dict, taxon: str, mean_depth: float,
                          read_len: float, thresholds: dict) -> tuple[dict, list]:
    """Demote ABSENT to UNKNOWN for this taxon's markers that depth could not sample.

    WHY THIS EXISTS. Marker typing survives low depth because it needs reads on a short
    locus rather than across a genome - that is what makes pathotype calling useful at
    detection-level titres. ABSENCE does not survive it. A single-copy 1.5 kb marker at
    0.012x mean depth expects 0.12 reads, so finding none is the expected outcome
    whether or not the locus is present, and recording ABSENT turns no data into
    positive evidence against a real finding.

    Returning UNKNOWN instead composes with the rest of the tiering: under Kleene rules
    UNKNOWN cannot satisfy confirm_requires, and an unevaluable veto or vaccine check
    already holds a finding at `supported` rather than confirming it.
    """
    floor = as_float(thresholds.get("min_expected_reads_for_absence"), 3.0)
    assumed_bp = as_float(thresholds.get("assumed_marker_bp"), 1000.0)
    adjusted = dict(states)
    demoted: list[tuple[str, float]] = []
    for mid, meta in metas.items():
        if meta.get("state") != "ABSENT":
            continue
        # Only this taxon's markers: depth is per taxon, and a neighbour's marker must
        # not be re-judged using this taxon's coverage.
        if taxon and meta.get("taxon") and meta["taxon"] != taxon:
            continue
        marker_bp = meta.get("nt_length") or assumed_bp
        exp = expected_marker_reads(mean_depth, marker_bp, read_len)
        if exp < floor:
            adjusted.pop(mid, None)          # absent from the dict == UNKNOWN
            demoted.append((mid, exp))
    return adjusted, sorted(demoted)


def describe_strain_claim(mean_depth: float, thresholds: dict) -> tuple[str, str]:
    """What kind of strain claim the observed depth can carry, and the reason.

    Says "not attempted" out loud rather than leaving a reader to assume genome-wide
    typing was checked. Nothing in the pipeline computes cgMLST, SNP typing or a
    consensus lineage yet; this records the depth finding either way so that when one
    is added the gate is already the thing deciding.
    """
    need = as_float(thresholds.get("min_depth_for_genome_typing"), 20.0)
    if mean_depth >= need:
        return ("genome_typing_supported", (
            f"mean depth {mean_depth:.2f}x >= {need:g}x, so genome-wide typing "
            f"(cgMLST/SNP/consensus lineage) is supportable - NOT YET IMPLEMENTED"))
    return ("marker_typing_only", (
        f"genome-wide strain typing NOT ATTEMPTED: mean depth {mean_depth:.2f}x "
        f"< {need:g}x required for confident SNP/cgMLST calling. Sub-species "
        f"evidence below is marker presence/absence only"))


def resolve_identity_floor(thresholds: dict, group: str,
                           platform: str) -> tuple[float, str]:
    """The read-identity floor for this group ON THIS PLATFORM, and why.

    Resolution order, most specific first:
        platform + group  ->  platform default  ->  group  ->  global default

    WHY PLATFORM IS IN HERE AT ALL. The floor was group-only, and `platform` did not
    appear anywhere in this file. The pack's default is 0.98, which is a statement
    about Illumina substitution error as much as about phylogenetic distance. ONT
    reads do not reach it: R10.4.1 simplex means ~97-98%, and the badread defaults
    the benchmark simulates with produce 95%. So every bacterial and fungal ONT
    finding failed the identity gate on sequencing error alone and capped at
    `supported`, at every spike level, with a reason string that reads exactly like
    genuinely weak evidence. Returning the reason alongside the number is deliberate:
    the findings row now says which rule was applied, so a floor that is too strict
    is visible in the output instead of having to be inferred.
    """
    by_plat = (thresholds.get("min_read_identity_by_platform") or {}).get(platform)
    if isinstance(by_plat, dict):
        if group and group in by_plat:
            return as_float(by_plat[group], thresholds["min_read_identity"]), \
                f"{platform}/{group}"
        if "default" in by_plat:
            return as_float(by_plat["default"], thresholds["min_read_identity"]), \
                f"{platform}/default"
    by_group = (thresholds.get("min_read_identity_by_group") or {})
    if group and group in by_group:
        return as_float(by_group[group], thresholds["min_read_identity"]), group
    # No platform recorded at all is worth saying out loud: it means the floor is the
    # Illumina-shaped global default, which on ONT data is unreachable.
    return thresholds["min_read_identity"], (
        f"{group or 'default'}, no platform recorded" if not platform
        else f"{group or 'default'}")


def resolve_alignment_floor(thresholds: dict, platform: str,
                            avg_read_len: float = 0.0) -> tuple[float, str]:
    """The minimum mean alignment length for this platform, and which rule applied.

    WHY THIS IS PER-PLATFORM AND WAS NOT. pack.yaml states the gate's rationale in its
    own comment - "a 50 bp match inside a 10 kb read says almost nothing about where
    that read came from" - which is a statement about ONT. It was then applied flat at
    100 bp to every platform, where it is not a weak-evidence filter but a read-length
    filter: a 2x101 library trims and soft-clips to a mean of 98-99, so the gate is
    UNREACHABLE and no finding on such a run can ever be confirmed.

    Measured on the Zymo isolate libraries (2x101, real reads): 81 of 82 rows failed it,
    including Salmonella enterica in its own single-organism library at 24,057,383 reads,
    breadth 0.9137, occupancy 0.9378, depth_cv 0.42 and neighbour_ratio 1928 - held at
    `supported` with the reason "failed: mean_alignment_length 99 >= 100". Across all
    nine runs on disk, every confirmed_in_silico call the project has ever produced is on
    150 bp simulated reads, and this is why.

    This is the fourth recurrence of the standing lesson that a threshold must be
    REACHABLE, not merely strict (D10, D11, D20). It stayed invisible because the only
    fixture in the suite sets mean_alignment_length to exactly 100.0 - the single value
    that passes.
    """
    # A FRACTION, BECAUSE THE PACK'S OWN RATIONALE IS A RATIO. "A 50 bp match inside a
    # 10 kb read says almost nothing about where that read came from" is a statement
    # about alignment length RELATIVE TO READ LENGTH, and a constant cannot express it.
    # Measured, the mismatch is total: the same 100 bp constant is 99% of a 2x101 read
    # (binding, and unreachable), 67% of a 150 bp read (slack), and 1.4% of the 7,276
    # bp median ONT read on disk - where the minimum observed mean_alignment_length is
    # 727.6, so the ONT branch could never fire on anything resembling real ONT data.
    #
    # It also fixes the defect ONE READ LENGTH DOWN. A flat 60 reproduces exactly what
    # a flat 100 did: a 2x50 library aligns ~49 bp and could never be confirmed, and
    # nothing on disk or in the suite would have shown it.
    frac = (thresholds.get("min_alignment_fraction_of_read") or {}).get(platform)
    if avg_read_len and avg_read_len > 0 and frac is not None:
        floor = as_float(thresholds.get("min_alignment_length_floor"), 40.0)
        want = max(floor, as_float(frac, 0.6) * avg_read_len)
        return want, (f"{platform} {as_float(frac, 0.6):.0%} of {avg_read_len:.0f} bp read"
                      + (f", floored at {floor:.0f}" if want == floor else ""))
    by_plat = (thresholds.get("min_alignment_length_by_platform") or {})
    if platform and platform in by_plat:
        return as_float(by_plat[platform], thresholds["min_alignment_length"]), (
            f"{platform}, read length unknown")
    return as_float(thresholds.get("min_alignment_length"), 100.0), (
        f"{platform or 'default'}, no platform rule")


def load_pack_thresholds(pack_yaml: Path) -> dict:
    """Read the pack's thresholds block.

    Uses bin/btf_yaml.py rather than the ad-hoc line scraper this function used to
    be. The scraper could not read a nested block, and min_read_identity_by_group is
    nested - a divergent virus needs a far lower identity floor than a bacterium
    whose near neighbour sits at 96%.
    """
    defaults = {
        "min_unique_reads": 20.0,
        "min_breadth": 0.005,
        "max_depth_cv": 2.5,
        "min_read_identity": 0.98,
        "min_mapped_fraction_of_assigned": 0.50,
        "min_alignment_length": 100.0,
        "min_alignment_length_by_platform": {},
        "min_alignment_fraction_of_read": {},
        "min_alignment_length_floor": 40.0,
        "min_read_identity_by_group": {},
        "min_read_identity_by_platform": {},
        "min_expected_reads_for_absence": 3.0,
        "assumed_marker_bp": 1000.0,
        "min_depth_for_genome_typing": 20.0,
    }
    if not pack_yaml or not pack_yaml.exists():
        return defaults
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import btf_yaml
        doc = btf_yaml.load(pack_yaml)
    except Exception as exc:
        # FATAL, NOT A FALLBACK.
        #
        # This used to warn and return the built-in defaults, whose
        # min_read_identity_by_group is EMPTY. A single editing slip in pack.yaml -
        # any construct btf_yaml rejects - therefore silently replaced the pack's
        # per-group identity floors with the flat 0.98 bacterial floor, downgrading
        # every virus confirmation in the run from confirmed_in_silico to supported.
        # The only trace was one stderr line in score.log, while the findings TSV
        # reason string read like a legitimate threshold decision. btf_yaml's own
        # contract is "fail loudly rather than mis-parse quietly"; honouring that
        # contract and then discarding the failure defeats it.
        sys.stderr.write(
            f"[score] FATAL: could not parse {pack_yaml}: {exc}\n"
            "[score] Refusing to score with built-in defaults: they carry no "
            "per-group or per-platform identity floors, so every virus and every "
            "ONT finding would be judged against the Illumina bacterial floor and "
            "silently downgraded. Fix the pack file.\n")
        sys.exit(2)
    block = (doc or {}).get("thresholds") or {}
    out = dict(defaults)
    for k, v in block.items():
        out[k] = v
    if not isinstance(out.get("min_read_identity_by_group"), dict):
        out["min_read_identity_by_group"] = {}
    if not isinstance(out.get("min_alignment_length_by_platform"), dict):
        out["min_alignment_length_by_platform"] = {}
    if not isinstance(out.get("min_alignment_fraction_of_read"), dict):
        out["min_alignment_fraction_of_read"] = {}
    if not isinstance(out.get("min_read_identity_by_platform"), dict):
        out["min_read_identity_by_platform"] = {}
    return out


def rank_order(tier: str) -> int:
    return TIER_ORDER.index(tier) if tier in TIER_ORDER else 0


def as_float(value, default=0.0) -> float:
    if value is None or value == "":
        return default
    if isinstance(value, str) and value.strip().lower() in ("inf", "+inf"):
        return math.inf
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def neighbour_is_credible(row: dict, args) -> bool:
    """Does the best neighbour look like an organism that is actually present?

    D12: for the case the attribution gate exists to handle, the denominator is background
    cross-mapping with ZERO reads from the declared neighbour - 2,700 reads on the
    Y. enterocolitica reference of which not one is Yersinia. Raw count cannot separate that
    from a real competitor (at 100 rpm: target 2,336, absent neighbour 2,731), so the gate
    turned min_neighbour_ratio into an undeclared read floor whose height is set by
    background composition. Measured: it put the Y. pestis cliff at 176.6 rpm, confirmed
    against a prediction of 171-184.

    Coverage SHAPE separates them, and is near-independent of abundance, which is what lets
    one threshold serve every spike level. Measured over runs/lod_v4 with the same
    ref-span denominator this reads:

        absent neighbour  (630)                     0.1532
        present neighbours (1642, 1428, 486, 562)   0.9720 - 0.9913

    ABSENT EVIDENCE IS NOT EVIDENCE OF INCREDIBILITY. A blank uniformity - no neighbour, or
    samtools coverage unavailable - returns True, leaving the gate exactly as it was. The
    change may only ever WITHHOLD an attribution on positive evidence that the neighbour is
    a pedestal; it must never manufacture one, and it must never silently disable the gate
    when the measurement is missing.
    """
    raw = (row.get("best_neighbour_uniformity") or "").strip()
    if not raw:
        return True
    try:
        return float(raw) >= args.min_neighbour_uniformity
    except ValueError:
        return True


def target_covers_its_reference(breadth: float, min_breadth: float,
                                occupancy: float, min_occupancy: float) -> bool:
    """Does the TARGET's own coverage geometry stand up, ignoring the neighbour?

    Deliberately the same two gates `--actionable-requires` defaults to, and for the
    same reason D14 chose them: a finding that does not cover its reference is not a
    detection, whatever else it passes. Those gates separated true from false
    positives on the Zymo isolates at median breadth 0.8661 against 0.0015.

    This is NOT a re-run of the attribution decision. It asks a prior question: is
    there anything here worth contesting? Where the answer is no - 84 of 84
    attributed rows in `lod_v4`, 39 of 45 in `zymo_rl100` - the attribution stands
    untouched, which is what keeps this from re-opening the false positives D18
    measured when the neighbour gate was applied to every row.
    """
    return breadth >= min_breadth and occupancy >= min_occupancy


def load_ref_names(manifest) -> dict:
    """taxid -> (name, role) from the pack's reference manifest.

    WHY A NAME IS NOT COSMETIC. `attributed_to_neighbour` is the only tier that
    reports a watchlist organism as absent on the strength of a DIFFERENT organism,
    and until now it named that organism by number: "reads favour neighbour taxid
    620". An analyst cannot audit that. The two calls this was written for are
    *E. coli* -> 620 and *C. neoformans* -> 40410; both read as opaque, and both are
    cases where the pack's own notes say the pair cannot be separated at species
    level from a metagenome. Spelled out as "Shigella" and "Cryptococcus
    deneoformans" the call is at least reviewable by the person holding it.

    The manifest's `role` column comes with it, because "the neighbour is itself a
    watchlist target" and "the neighbour is a near_neighbour reference that exists
    only to absorb cross-mapping" are different situations for a reader, and the
    distinction is already recorded.

    MISSING IS NOT FATAL. No manifest, or a taxid absent from it, falls back to the
    bare number - the behaviour before this function existed. Scoring must not
    depend on a cosmetic file.
    """
    out: dict = {}
    if not manifest or manifest.name == "NO_FILE" or not manifest.exists():
        return out
    try:
        with manifest.open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                tid = (row.get("taxid") or "").strip()
                if tid and tid not in out:
                    out[tid] = ((row.get("name") or "").strip(),
                                (row.get("role") or "").strip())
    except Exception as exc:
        sys.stderr.write(f"[score] WARN: could not read {manifest}: {exc}; "
                         "neighbours will be named by taxid only\n")
    return out


def neighbour_label(taxid: str, names: dict) -> str:
    """Human-readable identification of a neighbour, for a reason string."""
    tid = (taxid or "").strip()
    if not tid:
        return "an unidentified neighbour"
    name, role = names.get(tid, ("", ""))
    if not name:
        return f"taxid {tid}"
    qual = " (near-neighbour reference)" if role == "near_neighbour" else ""
    return f"{name} [taxid {tid}]{qual}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coverage", nargs="*", type=Path, default=[],
                    help="per-(sample,taxon) TSVs from coverage_stats.py")
    ap.add_argument("--candidates", nargs="*", type=Path, default=[],
                    help="per-sample stage-1 candidate TSVs from parse_screen.py")
    ap.add_argument("--missing-refs", type=Path, default=None)
    ap.add_argument("--samplesheet", type=Path,
                    help="the run's samplesheet. THE ONLY AUTHORITATIVE SOURCE for "
                         "which samples are controls: a control that raises nothing "
                         "produces a header-only candidates.tsv and no coverage rows, "
                         "so counting controls from either output misses the cleanest "
                         "case - the one where the reassurance matters most.")
    ap.add_argument("--watchlist", required=True, type=Path,
                    help="watchlist.resolved.tsv")
    ap.add_argument("--pack-yaml", type=Path, default=None)
    ap.add_argument("--read-stats", nargs="*", type=Path, default=[],
                    help="seqkit_stats TSVs. Supplies each sample's avg_len, which turns "
                         "the alignment-length floor from a constant into a fraction of "
                         "the read. SCORE_EVIDENCE was the only process not given these, "
                         "which is why that gate was a constant at all.")
    ap.add_argument("--refs-manifest", type=Path, default=None,
                    help="pack refs manifest.tsv; supplies the NAME of an "
                         "attributed-to neighbour. Absent = name by taxid only.")
    ap.add_argument("--min-occupancy", type=float, default=0.15,
                    help="minimum window occupancy for confirmation [0.15]")
    ap.add_argument("--min-neighbour-ratio", type=float, default=1.5,
                    help="target reads must exceed best neighbour by this factor [1.5]")
    ap.add_argument("--actionable-requires", default="breadth,occupancy",
                    help="comma-separated confirmation gates that MUST pass for a finding "
                         "to be actionable; failing one drops it to screened_only instead "
                         "of `supported`. Empty string restores the old behaviour, where "
                         "`supported` was the bin for any finding that cleared unique_reads "
                         "and neighbour_ratio however badly it then failed.")
    # DEFAULT 0, MATCHING nextflow.config, BECAUSE A DEFAULT THAT DISAGREES WITH WHAT
    # SHIPS IS A TRAP FOR EVERY MANUAL RE-SCORE. D18 measured this gate on real reads and
    # set `min_neighbour_uniformity = 0` in nextflow.config: on the Zymo isolates it ADDED
    # false positives (25 vs 4 in zymo_rl100) without rescuing the true positives it exists
    # for. The argparse default stayed 0.50, so every re-scoring of a published batch that
    # forgot the flag silently reproduced PRE-D18 behaviour. That is not hypothetical - it
    # happened during the D16a review (2026-10-13) and produced 7 tier changes where the
    # shipped configuration produces 5, which read as an error in docs/decisions.md that
    # was not there. 0.50 remains available; it is no longer what you get by accident.
    ap.add_argument("--min-neighbour-uniformity", type=float, default=0.0,
                    help="a neighbour below this coverage uniformity is a cross-mapping "
                         "pedestal, not a competing organism, and cannot out-vote the "
                         "target (D12). Calibrated gap: 0.1532 absent vs 0.9720-0.9913 "
                         "present. DEFAULT 0 (disabled), matching nextflow.config: D18 "
                         "measured the gate on real reads and it cost more than it "
                         "bought. Set 0.50 to restore the D12 behaviour.")
    ap.add_argument("--contest-attribution", action="store_true",
                    help="when a neighbour out-counts the target BUT the target's own "
                         "breadth and occupancy pass, do not let the attribution stand "
                         "silently: either the neighbour is a cross-mapping pedestal "
                         "(judge the target on its own evidence) or both references are "
                         "genuinely covered (report indistinguishable_complex). "
                         "--pedestal-max-uniformity decides which.")
    ap.add_argument("--pedestal-max-uniformity", type=float, default=0.30,
                    help="with --contest-attribution: a neighbour at or below this "
                         "coverage uniformity is pile-up on shared sequence, not a "
                         "present organism [0.30]")
    ap.add_argument("--control-ratio", type=float, default=5.0,
                    help="sample RPM must exceed control RPM by this factor [5.0]")
    ap.add_argument("--marker-evidence", nargs="*", type=Path, default=[],
                    help="per-sample marker state TSVs from marker_evidence.py")
    ap.add_argument("--marker-rules", type=Path, default=None,
                    help="packs/<pack>/marker_rules.yaml")
    ap.add_argument("--threshold", action="append", default=[], metavar="KEY=VALUE",
                    help="Override one pack.yaml threshold, repeatable. pack.yaml is "
                         "the source of truth; this exists so a profile or a CLI flag "
                         "can relax a gate for a smoke test WITHOUT editing the pack. "
                         "Seven params used to be declared in nextflow.config and "
                         "reach no process at all: -profile test set min_unique_reads=1 "
                         "and min_breadth=0.0 intending to relax stage 2, both were "
                         "silently ignored, and the tiny test input produced zero "
                         "confirmations that looked like a clean negative.")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.actionable_requires = [x.strip() for x in
                                (args.actionable_requires or "").split(",") if x.strip()]

    thresholds = apply_threshold_overrides(load_pack_thresholds(args.pack_yaml),
                                           args.threshold)

    # taxid -> (name, role). Used ONLY to name an attributed-to neighbour in prose;
    # no gate reads it, so a missing manifest degrades the wording and nothing else.
    ref_names = load_ref_names(args.refs_manifest)

    # sample -> mean read length, from the stats the pipeline already computes.
    read_len_of: dict = {}
    for path in args.read_stats:
        for r in read_tsv(path):
            samp, avg = (r.get("sample") or "").strip(), as_float(r.get("avg_len"))
            if samp and avg > 0:
                read_len_of[samp] = max(read_len_of.get(samp, 0.0), avg)
    if read_len_of:
        sys.stderr.write(f"[score] read lengths known for {len(read_len_of)} sample(s)\n")

    # ---- marker evidence and rules -------------------------------------
    # marker_state[sample][marker_id] = mr.TRUE / mr.FALSE. A marker that is
    # UNKNOWN is deliberately ABSENT FROM THE DICT, so the rule engine's default
    # of UNKNOWN applies and indeterminacy propagates correctly.
    marker_state: dict[str, dict[str, int]] = defaultdict(dict)
    # Metadata per (sample, marker), so the depth gate below can ask whether a marker
    # that was NOT found had any chance of being sampled at the observed depth.
    marker_meta: dict[str, dict[str, dict]] = defaultdict(dict)
    marker_seen = 0
    for path in args.marker_evidence:
        for row in read_tsv(path):
            state = (row.get("state") or "").strip().upper()
            sample = row.get("sample", "")
            mid = row.get("marker_id", "")
            if not mid or not sample or mr is None:
                continue
            marker_seen += 1
            marker_meta[sample][mid] = {
                "taxon": (row.get("taxon") or "").strip(),
                "marker_class": (row.get("marker_class") or "").strip(),
                "nt_length": as_float(row.get("nt_length"), 0.0),
                "nt_reads": as_float(row.get("nt_reads"), 0.0),
                "state": state,
            }
            if state == "PRESENT":
                marker_state[sample][mid] = mr.TRUE
            elif state == "ABSENT":
                marker_state[sample][mid] = mr.FALSE
            # UNKNOWN: intentionally not recorded.

    marker_rules: dict[str, dict] = {}
    if args.marker_rules and mr is not None:
        if not args.marker_rules.exists():
            sys.exit(f"ERROR: --marker-rules {args.marker_rules} not found")
        try:
            marker_rules = mr.load_rules(args.marker_rules)
        except mr.RuleError as exc:
            sys.exit(f"ERROR: cannot parse {args.marker_rules}: {exc}")
        sys.stderr.write(
            f"[score] marker gating active: {len(marker_rules)} rule set(s), "
            f"{marker_seen} marker observation(s)\n"
        )
    elif args.marker_evidence:
        sys.stderr.write(
            "[score] WARNING: marker evidence supplied without --marker-rules; "
            "marker gating is INACTIVE and taxa needing marker confirmation may be "
            "promoted on alignment evidence alone\n"
        )

    # Per-taxon threshold overrides. A single global breadth cutoff cannot serve
    # a 5 kb virus and a 23 Mb Plasmodium genome at once.
    watch = {}
    for row in read_tsv(args.watchlist):
        watch[row["taxid"]] = row
    if not watch:
        sys.exit(f"ERROR: no entries read from watchlist {args.watchlist}")

    # Stage-1 support, for concordance reporting.
    screen: dict[tuple[str, str], dict] = {}
    # Control samples seen ANYWHERE, stage 1 included. n_controls used to be counted from
    # coverage rows alone, which exist only for taxa that reached stage 2 - so the cleanest
    # possible NTC, the one that nominates nothing, counted as zero controls and every
    # finding in the run carried "NO NEGATIVE CONTROL IN RUN - contamination cannot be
    # excluded". The run HAD a control; it was clean. That is the one QC line an analyst
    # reads first, and it asserted the opposite of the truth in the best case.
    control_samples: set[str] = set()
    for path in args.candidates:
        for row in read_tsv(path):
            key = (row.get("sample", ""), row.get("taxid", ""))
            screen[key] = row
            if row.get("role") == "control" and row.get("sample"):
                control_samples.add(row["sample"])

    coverage: list[dict] = []
    for path in args.coverage:
        coverage.extend(read_tsv(path))

    # Highest control signal per taxon, used to subtract run-wide contamination.
    control_rpm: dict[str, float] = defaultdict(float)
    control_src: dict[str, str] = {}
    for row in coverage:
        if row.get("role") == "control":
            rpm = as_float(row.get("reads_per_million"))
            # .get(), NOT [] - control_rpm is a defaultdict, so a bare read INSERTS the
            # key even when the comparison is false. Every control with any coverage row,
            # including one at 0.00 RPM, then made `control_rpm` truthy and the
            # "controls all clean" branch below unreachable.
            if rpm > control_rpm.get(row["taxid"], 0.0):
                control_rpm[row["taxid"]] = rpm
                control_src[row["taxid"]] = row.get("sample", "")
    control_samples |= {r["sample"] for r in coverage if r.get("role") == "control"}
    # The samplesheet wins when present. Counting from pipeline OUTPUT asks "did a control
    # produce a row?", a different question from "did the run have a control?", and answers
    # it wrongly in exactly the good case: parse_screen writes a header-only candidates.tsv
    # for a sample that raises nothing, so a spotless NTC contributes no rows anywhere and
    # counted as zero controls - putting "NO NEGATIVE CONTROL IN RUN - contamination cannot
    # be excluded" on every finding of the cleanest possible run. An earlier fix moved the
    # count from coverage rows to candidate rows, which only widened it from "control
    # reached stage 2" to "control reached stage 1" and left the spotless case broken.
    if args.samplesheet and args.samplesheet.exists():
        with args.samplesheet.open(newline="") as fh:
            for r in csv.DictReader(fh):
                if (r.get("role") or "").strip().lower() == "control" and r.get("sample"):
                    control_samples.add(r["sample"].strip())
    n_controls = len(control_samples)

    # CONTROL SUBTRACTION HAS BEEN INVISIBLE. Across the six published batches it moved
    # 108 rows, and the only record of any of it was prose inside the `reasons` column -
    # every number in the 2026-10-13 audit had to be recovered with a regex over that
    # prose, which is precisely the artifact that breaks silently when a message is
    # reworded. Worse, a control that carries a watchlist organism but happens to
    # suppress nothing today left NO trace at all: lod_v3_threearm's NTC holds
    # B. anthracis at 32.58 RPM, which at control_ratio 5.0 is a standing 162.9 RPM
    # floor on a select agent, and nothing anywhere said so. These counters exist so the
    # subtraction can be read off a table instead of reconstructed.
    #
    # `rows` and `eligible` differ on purpose. Suppression lives in the final else of the
    # chain that starts at the read floor, so a row below min_unique_reads or already
    # taken by attributed_to_neighbour NEVER REACHES IT. Counting only what reached the
    # branch would understate every control's true reach; counting only the rows would
    # imply a suppression that the ladder had already pre-empted. Both are reported.
    control_reach: dict[str, dict] = defaultdict(
        lambda: {"rows": 0, "eligible": 0, "suppressed": 0, "exceeds": 0})

    findings: list[dict] = []

    for row in coverage:
        taxid = row["taxid"]
        meta = watch.get(taxid, {})
        # Machine-readable record of what control subtraction did to THIS row. The
        # prose in `reasons` is for the analyst; this column is so that
        # evaluate_benchmark.py can one day count the false negatives subtraction
        # causes, which it cannot do today.
        control_subtraction = ""

        # `as_float(x, default)` already falls back when the field is absent or
        # unparseable. The trailing `or thresholds[...]` that used to be here also
        # swallowed a LEGITIMATE ZERO: a watchlist row asking for
        # min_unique_reads: 0 ("report every read for this select agent") had its
        # request replaced by the pack default, and the reason string then cited a
        # threshold the operator never set.
        min_reads = as_float(meta.get("min_unique_reads"),
                             thresholds["min_unique_reads"])
        min_breadth = as_float(meta.get("min_breadth"), thresholds["min_breadth"])

        if row.get("role") != "control":
            control_reach[taxid]["rows"] += 1

        # HOISTED OUT OF THE TIER CHAIN ON PURPOSE. These used to be bound inside the
        # final `else`, which meant a row pre-empted by the read floor or by
        # attributed_to_neighbour never bound them at all - and Python leaks loop locals
        # across iterations, so reading them from outside the `else` would have handed
        # that row the PREVIOUS row's control name. A wrong attribution that looks right
        # is worse than a blank one, and this file already carries that exact bug's scar
        # for `platform`. Binding unconditionally, once per row, removes the class.
        ctrl = control_rpm.get(taxid, 0.0)
        ctrl_name = control_src.get(taxid, "?")

        unique_reads = as_float(row.get("unique_reads"))
        breadth = as_float(row.get("breadth"))
        occupancy = as_float(row.get("window_occupancy"))
        depth_cv = as_float(row.get("depth_cv"))
        identity = as_float(row.get("mean_identity"))
        nb_ratio = as_float(row.get("neighbour_ratio"), math.inf)
        rpm = as_float(row.get("reads_per_million"))

        # Read per-row BEFORE any branching. These were assigned only inside the
        # innermost else, but `platform` is written into every finding and Python
        # leaks loop locals across iterations, so a row that short-circuited on
        # unique_reads / attributed_to_neighbour / suppressed_control either crashed
        # with UnboundLocalError (if it was the FIRST row) or silently inherited the
        # PREVIOUS row's platform - into the very field a reader uses to tell an ONT
        # confirmation from an Illumina one.
        grp = (row.get("group") or meta.get("group") or "").strip().lower()
        plat = (row.get("platform") or "").strip().lower()

        # Computed BEFORE the branch because the attribution branch now has to ask
        # whether there is anything worth contesting, and that question is about the
        # target's own evidence - which the attribution branch previously never looked at.
        _nb_raw = (row.get("best_neighbour_uniformity") or "").strip()
        try:
            nb_unif = float(_nb_raw) if _nb_raw else None
        except ValueError:
            nb_unif = None
        # A ROW WITH NO MEASUREMENT IS NOT CONTESTED, AND THE FIRST VERSION OF THIS GOT
        # THAT BACKWARDS. It routed a blank to `indistinguishable_complex`, reasoning
        # from "absent evidence is not evidence" - but that principle, which
        # neighbour_is_credible() states in its own docstring, says a missing
        # measurement must leave the decision EXACTLY AS IT WAS, not manufacture a
        # different one. The symptom is concrete: runs/lod_v1 and runs/lod_v3_threearm
        # pre-date the column entirely, so byte-identical rows scored `supported` in one
        # run and `indistinguishable_complex` in another, and the latter ranks lower -
        # so the very expectation D19 cited as proof would have FAILED on them. The live
        # risk is not old runs, it is covered_by_taxid() returning {} when samtools
        # coverage fails, which blanks every row in that sample at once.
        contested = (nb_unif is not None
                     and target_covers_its_reference(breadth, min_breadth,
                                                     occupancy, args.min_occupancy))

        reasons: list[str] = []
        tier = "screened_only"

        if unique_reads < min_reads:
            reasons.append(f"unique_reads {unique_reads:.0f} < {min_reads:.0f}")
        elif (nb_ratio < args.min_neighbour_ratio and neighbour_is_credible(row, args)
              and not (args.contest_attribution and contested)):
            tier = "attributed_to_neighbour"
            reasons.append(
                f"reads favour "
                f"{neighbour_label(row.get('best_neighbour_taxid'), ref_names)} "
                f"(ratio {row.get('neighbour_ratio')} < {args.min_neighbour_ratio}; "
                f"neighbour uniformity {row.get('best_neighbour_uniformity') or 'n/a'})"
            )
        elif (args.contest_attribution and contested
              and nb_ratio < args.min_neighbour_ratio
              and nb_unif > args.pedestal_max_uniformity):
            # BOTH REFERENCES ARE GENUINELY COVERED, SO NEITHER CALL IS SAFE.
            #
            # The attribution machinery assumes one of the two organisms is real and
            # the other is cross-mapping. For a sibling pair that assumption fails:
            # Shigella IS a pathotype within E. coli, and C. deneoformans was
            # C. neoformans var. neoformans until 2015. Both references then cover
            # well, the read counts land within ~10% of each other, and whichever way
            # the ratio falls the pipeline printed a confident one-sided answer.
            #
            # Measured on the Zymo single-organism libraries, where the truth is
            # known by construction: E. coli at 3,349,860 reads and breadth 0.667 was
            # reported ABSENT because Shigella took the vote at ratio 0.9032, and
            # C. neoformans at 5,870,155 reads and breadth 0.642 the same way at
            # 0.7072. Neither neighbour is on the watchlist, so NO finding anywhere
            # carried those reads: the entire output for a library that is 100% a
            # notifiable pathogen was `attributed_to_neighbour` and `screened_only`.
            #
            # "One of these two is present and this method cannot say which" is a
            # true statement and a useful one. A one-sided attribution is neither.
            tier = "indistinguishable_complex"
            reasons.append(
                f"CANNOT SEPARATE: the target covers its own reference (breadth "
                f"{breadth:.4f} >= {min_breadth:.4f}, occupancy {occupancy:.4f} >= "
                f"{args.min_occupancy:.4f}) AND "
                f"{neighbour_label(row.get('best_neighbour_taxid'), ref_names)} "
                + (f"covers its own (uniformity {nb_unif:.4f} > "
                   f"{args.pedestal_max_uniformity:.2f})")
                + f" while out-counting it only "
                f"{nb_ratio:.4f}:1. Both references are genuinely covered, so the "
                f"reads cannot be assigned to one organism at this resolution. One of "
                f"the two IS present; which one is not determined by this evidence - "
                f"resolve with lineage markers or a targeted assay."
            )
        else:
            # WITHHELD ATTRIBUTION MUST FALL THROUGH, NOT SHORT-CIRCUIT. The first version
            # of this put the not-credible case in its own `elif`, which appended the
            # reason and then skipped the gate block below - leaving `tier` at its initial
            # "screened_only". The comment said it fell through to the confirmation gates
            # and it did not, so a 170 rpm sample was scored BELOW a 190 rpm one not
            # because its evidence was weaker but because it was never evaluated. Caught by
            # the d12fine revalidation: the cliff did not move, which is the symptom.
            if (args.contest_attribution and contested
                    and nb_ratio < args.min_neighbour_ratio
                    and nb_unif <= args.pedestal_max_uniformity):
                reasons.append(
                    f"ATTRIBUTION CONTESTED AND REJECTED: "
                    f"{neighbour_label(row.get('best_neighbour_taxid'), ref_names)} "
                    f"out-counts the target {nb_ratio:.4f}:1, but the target covers its "
                    f"own reference (breadth {breadth:.4f} >= {min_breadth:.4f}, "
                    f"occupancy {occupancy:.4f} >= {args.min_occupancy:.4f}) while the "
                    f"neighbour's uniformity is {nb_unif:.4f} <= "
                    f"{args.pedestal_max_uniformity:.2f} - its reads pile on shared "
                    f"sequence rather than covering its genome, so it is not in the "
                    f"specimen and cannot out-vote (D12). Judged on its own evidence below."
                )
            elif nb_ratio < args.min_neighbour_ratio:
                reasons.append(
                    f"neighbour "
                    f"{neighbour_label(row.get('best_neighbour_taxid'), ref_names)} "
                    f"out-counts "
                    f"the target (ratio {row.get('neighbour_ratio')}) but is NOT CREDIBLE: "
                    f"coverage uniformity {row.get('best_neighbour_uniformity')} < "
                    f"{args.min_neighbour_uniformity} means its reads pile on shared "
                    f"sequence rather than covering its genome, so this is background "
                    f"cross-mapping and not a competing organism (D12). Attribution "
                    f"withheld; judged on its own evidence below."
                )
            if row.get("role") != "control" and ctrl > 0:
                control_reach[taxid]["eligible"] += 1
            if row.get("role") != "control" and ctrl > 0 and rpm < args.control_ratio * ctrl:
                tier = "suppressed_control"
                control_reach[taxid]["suppressed"] += 1
                reasons.append(
                    f"control {ctrl_name} has {ctrl:.2f} RPM; "
                    f"sample {rpm:.2f} RPM does not exceed {args.control_ratio}x"
                )
                # THE SAMPLE HAS MORE OF IT THAN THE BLANK DOES, AND ONLY THE MARGIN IS
                # HOLDING IT DOWN. Measured 2026-10-13: 8 rows in lod_v4 and 3 in
                # d12fine are suppressed at a sample/control ratio of 1.01-1.03, i.e.
                # the contamination model did not decide them - control_ratio 5.0 did,
                # and nothing in docs/decisions.md derives that 5.0 from anything. A
                # reader seeing "does not exceed 5.0x" reasonably assumes the control
                # held MORE of the organism. In these rows it held less.
                if rpm > ctrl:
                    control_reach[taxid]["exceeds"] += 1
                    control_subtraction = f"exceeds_control:{ctrl_name}@{ctrl:.2f}"
                    reasons.append(
                        f"NOTE: this sample carries MORE of the organism than the "
                        f"control does ({rpm:.2f} vs {ctrl:.2f} RPM, "
                        f"{rpm / ctrl:.2f}x) - it is suppressed by the "
                        f"{args.control_ratio}x MARGIN, not by having less than the "
                        f"blank. Review before treating this as a negative."
                    )
                else:
                    control_subtraction = f"suppressed:{ctrl_name}@{ctrl:.2f}"
            else:
                if row.get("role") != "control" and ctrl > 0:
                    control_subtraction = f"cleared:{ctrl_name}@{ctrl:.2f}"
                # Identity floor is group-aware: 98% is the near-neighbour
                # discriminator for bacteria (Wright et al. found ~96% meant
                # genus-correct, species-wrong), but a divergent virus is
                # legitimately far from any available reference and the same floor
                # would manufacture false negatives where sensitivity matters most.
                id_floor, floor_src = resolve_identity_floor(thresholds, grp, plat)

                gates = {
                    f"breadth {breadth:.4f} >= {min_breadth:.4f}": breadth >= min_breadth,
                    f"occupancy {occupancy:.4f} >= {args.min_occupancy:.4f}":
                        occupancy >= args.min_occupancy,
                    f"identity {identity:.4f} >= {id_floor:.4f} ({floor_src})":
                        identity >= id_floor,
                    f"depth_cv {depth_cv:.2f} <= {thresholds['max_depth_cv']:.2f}":
                        depth_cv <= thresholds["max_depth_cv"],
                }
                # Minimum alignment length. Short alignments are a documented
                # false-positive source on ONT data specifically, where a 50 bp
                # match inside a 10 kb read says almost nothing about the read's
                # origin. Only applied when the metric is present.
                mean_aln = as_float(row.get("mean_alignment_length"), -1.0)
                if mean_aln >= 0:
                    aln_floor, aln_src = resolve_alignment_floor(
                        thresholds, plat, read_len_of.get(row.get("sample", ""), 0.0))
                    gates[f"mean_alignment_length {mean_aln:.0f} >= {aln_floor:.0f} "
                          f"({aln_src})"] = mean_aln >= aln_floor

                # Wright et al. criterion (ii). Only applied when the denominator is
                # known - absent stage-1 counts must not fabricate a failed gate.
                mapped_frac = as_float(row.get("mapped_fraction_of_assigned"), -1.0)
                if mapped_frac >= 0:
                    floor = as_float(thresholds.get("min_mapped_fraction_of_assigned"), 0.5)
                    gates[f"mapped_fraction {mapped_frac:.3f} >= {floor:.2f}"] = \
                        mapped_frac >= floor
                failed = [k for k, ok in gates.items() if not ok]
                if not failed:
                    tier = "confirmed_in_silico"
                    reasons.append("all confirmation gates passed: " + "; ".join(gates))
                else:
                    # `supported` IS AN ASSERTION OF PRESENCE, SO IT CANNOT BE THE BIN FOR
                    # "FAILED EVERYTHING". It used to be: any finding clearing unique_reads
                    # and neighbour_ratio landed here however badly it then failed, and the
                    # evaluator maps `supported` to POSITIVE. Measured on the Zymo
                    # single-organism libraries - where every watchlist call but the
                    # library's own organism is false by construction - that produced 45
                    # actionable false positives against 4 true ones, including B. anthracis,
                    # Y. pestis and B. pseudomallei, at a median breadth of 0.0015 and
                    # depth_cv of 78.9.
                    #
                    # The gates were not wrong; the ladder discarded their answer. True
                    # positives in that data fail a median of 1 gate, false positives 5.
                    # Requiring the coverage-geometry gates to PASS keeps 4 of 4 true
                    # positives and drops 45 false positives to 4:
                    #
                    #     fail <= 2 gates                 TP 4/4   FP 4
                    #     occupancy must pass             TP 4/4   FP 4
                    #     breadth AND occupancy must pass TP 4/4   FP 4   <- this
                    #     breadth must pass               TP 4/4   FP 6
                    #     depth_cv must pass              TP 3/4   FP 3
                    #
                    # Breadth and occupancy are chosen over a bare failure count because
                    # they say something: a finding that does not cover its reference is not
                    # a detection, whatever else it passes. coverage_stats.py's own
                    # docstring calls them the most useful metrics it computes.
                    # NOT FOR CONTROLS. This rule governs whether a finding is
                    # ACTIONABLE, and a control row never is: its content is QC about the
                    # control itself, reclassified to control_background below. Applying
                    # the rule to controls sent weak control rows to screened_only, which
                    # is not in that reclassification's list, so the row silently stopped
                    # being reported as control background - losing the QC signal without
                    # preventing any false positive. Caught by tests/test_pipeline.py
                    # `tier:NTC`, which asserts exactly this.
                    blocking = [] if row.get("role") == "control" else [
                        k for k in failed
                        if any(k.startswith(g) for g in args.actionable_requires)]
                    if blocking and args.actionable_requires:
                        tier = "screened_only"
                        reasons.append(
                            f"NOT ACTIONABLE: failed {'; '.join(blocking)} - a finding that "
                            f"does not cover its reference is not a detection. Also failed: "
                            f"{'; '.join(k for k in failed if k not in blocking) or 'nothing else'}")
                    else:
                        tier = "supported"
                        reasons.append(f"read threshold met; failed: {'; '.join(failed)}")
                if ctrl > 0:
                    reasons.append(
                        f"control present at {ctrl:.2f} RPM, sample {rpm:.2f} RPM "
                        f"({rpm / ctrl:.1f}x)"
                    )

        # ---- masked-region support ------------------------------------
        # A finding whose support sits in masked regions is the documented Wright
        # et al. failure: contaminating phage or vector sequence inside a reference
        # collecting reads and manufacturing a confident call. It usually fails the
        # read floor on its own once masked reads are excluded, but it must say so
        # rather than just quietly scoring low.
        masked_reads = as_float(row.get("reads_in_masked_regions"), 0.0)
        masked_frac_ref = as_float(row.get("masked_fraction_of_reference"), 0.0)
        if masked_reads > 0:
            if masked_reads >= max(1.0, unique_reads):
                reasons.append(
                    f"CAVEAT: {masked_reads:.0f} read(s) fell entirely inside MASKED "
                    f"reference regions vs {unique_reads:.0f} trustworthy read(s) - "
                    "support is dominated by vector, host, low-complexity or "
                    "cross-contaminated sequence in the reference itself")
                if tier == "confirmed_in_silico":
                    tier = "supported"
            else:
                reasons.append(
                    f"{masked_reads:.0f} read(s) in masked reference regions "
                    "(excluded from the read count)")
        if masked_frac_ref >= 0.20:
            reasons.append(
                f"CAVEAT: {masked_frac_ref:.0%} of this reference is masked, so breadth "
                "is computed over a small trusted fraction - treat the reference as "
                "suspect and consider recuration")

        # ---- depth-conditional typing gate ---------------------------
        # Runs BEFORE rule evaluation, so a marker absence that the depth could not
        # have established never reaches the rules as FALSE.
        this_taxon = row.get("taxon") or meta.get("name", "")
        this_sample = row.get("sample", "")
        mean_depth = as_float(row.get("mean_depth"), 0.0)
        read_len = as_float(row.get("mean_alignment_length"), 0.0) or 150.0
        metas = marker_meta.get(this_sample, {})
        eff_state, depth_demoted = gate_absence_on_depth(
            marker_state.get(this_sample, {}), metas, this_taxon,
            mean_depth, read_len, thresholds)
        strain_claim, strain_reason = describe_strain_claim(mean_depth, thresholds)
        if depth_demoted:
            shown = ", ".join(f"{m}~{e:.2f} reads" for m, e in depth_demoted[:4])
            more = ", ..." if len(depth_demoted) > 4 else ""
            need_reads = as_float(
                thresholds.get("min_expected_reads_for_absence"), 3.0)
            reasons.append(
                f"DEPTH GATE: {len(depth_demoted)} marker absence(s) for this taxon "
                f"reported UNKNOWN rather than ABSENT - at {mean_depth:.3f}x mean "
                f"depth fewer than {need_reads:g} reads were expected on them, so not "
                f"finding them is uninformative ({shown}{more})")

        # Sub-species typing evidence: the marker classes that carry a strain claim.
        typing_bits = []
        for mid, m in sorted(metas.items()):
            if m.get("marker_class") not in TYPING_CLASSES:
                continue
            if this_taxon and m.get("taxon") and m["taxon"] != this_taxon:
                continue
            st = ("PRESENT" if eff_state.get(mid) == mr.TRUE
                  else "ABSENT" if eff_state.get(mid) == mr.FALSE else "UNKNOWN")
            typing_bits.append(f"{m['marker_class']}:{mid}={st}")
        strain_typing = "; ".join(typing_bits)

        # ---- marker gating -------------------------------------------
        marker_summary = ""
        spec = marker_rules.get(this_taxon)
        # An explicit, auditable declaration that this taxon HAS no usable markers, as
        # opposed to nobody having written a block for it. Silence must not grant
        # permission; a key in the pack must.
        declared_no_markers = bool(spec) and bool(spec.get("no_markers_available"))
        if (spec and not declared_no_markers and mr is not None
                and tier in ("supported", "confirmed_in_silico")):
            verdict = mr.evaluate(row.get("taxon", ""), spec, eff_state)
            marker_summary = "; ".join(verdict["reasons"])

            if verdict["veto"] == mr.TRUE:
                tier = "marker_veto"
                reasons.append(f"MARKER VETO: {spec.get('veto_if')} is satisfied")
            elif verdict["vaccine"] == mr.TRUE:
                tier = "vaccine_strain"
                reasons.append(
                    f"VACCINE/ATTENUATED STRAIN indicated ({spec.get('vaccine_if')}). "
                    "Not wild type - do not report as disease.")
            elif verdict["support"] == mr.FALSE:
                tier = "screened_only"
                reasons.append(
                    f"minimum marker evidence not met ({spec.get('support_requires')})")
            # AN UNTESTED VETO OR VACCINE CHECK NOW CAPS THE TIER, NOT JUST THE PROSE.
            #
            # Kleene logic is right everywhere else: UNKNOWN in confirm_requires holds a
            # finding at `supported`, the conservative direction. But UNKNOWN in
            # veto_if/vaccine_if means NO DOWNGRADE, which is the unsafe direction, and
            # this used to disclose it in a caveat and confirm anyway. That is not a
            # hypothetical gap. Cross-checking every rule against the built marker
            # manifest, 18 referenced markers can only ever be UNKNOWN (5 failed, 4
            # curation_required, 9 motif_todo), so:
            #   * 8 of 9 vaccine_if rules can never fire - Variola, Enterovirus C,
            #     Measles, Yellow fever, Rabies, VEEV, Junin, RVFV. Only BCG and RB51 can.
            #   * 6 veto_if rules can never fire, including both near-neighbour vetoes
            #     the README highlights: B. thailandensis -> B. pseudomallei and
            #     F. novicida -> F. tularensis.
            # The concrete consequence: a Sabin vaccine strain satisfies Enterovirus C's
            # confirm_requires (polio_VP1 is built) while polio_sabin is UNKNOWN, so the
            # finding reached `confirmed_in_silico` - a vaccine strain reported as
            # confirmed poliovirus, which is the report that starts an outbreak response.
            #
            # The asymmetry decides it. "Supported, vaccine status untestable" costs an
            # analyst a follow-up; "confirmed wild-type poliovirus" costs a public health
            # investigation. So an unevaluable vaccine or veto check holds the finding at
            # `supported`. The finding is still reported, and prominently - it simply is
            # not stamped as confirmed wild type on evidence that was never gathered.
            if verdict["veto"] == mr.UNKNOWN:
                reasons.append(
                    f"CAVEAT: veto check could not be evaluated ({spec.get('veto_if')}) "
                    "- untested marker(s); an argument against this call may exist "
                    "but was not tested")
                if tier == "confirmed_in_silico":
                    tier = "supported"
                    reasons.append(
                        "HELD at supported: a near-neighbour veto that could not be "
                        "evaluated cannot be treated as having passed")
            if verdict["vaccine"] == mr.UNKNOWN:
                reasons.append(
                    f"CAVEAT: vaccine/attenuated-strain check could not be evaluated "
                    f"({spec.get('vaccine_if')}) - a vaccine strain CANNOT be excluded "
                    "for this finding")
                if tier == "confirmed_in_silico":
                    tier = "supported"
                    reasons.append(
                        "HELD at supported: confirming wild type requires actually "
                        "excluding the vaccine strain, not merely failing to test for it")

            # Only evaluate confirm_requires if veto/vaccine/support has not already
            # decided the tier. Without this guard the confirm branch OVERWRITES a
            # vaccine_strain determination: BCG (rpoB present, RD9 absent) sets
            # vaccine_strain, then confirm_requires = rpoB AND RD9 evaluates FALSE and
            # drops it back to `supported` - reporting a vaccine strain as a possible
            # tuberculosis detection, which is precisely the wrong direction. Introduced
            # when the veto/vaccine CAVEAT disclosures broke the original elif chain;
            # caught by tests/test_pipeline.py::marker:BCG.
            already_decided = tier in ("marker_veto", "vaccine_strain", "screened_only")
            if not already_decided and (tier == "confirmed_in_silico"
                                        or verdict["confirm"] is not None):
                if verdict["confirm"] == mr.TRUE:
                    reasons.append(f"marker rule satisfied: {spec.get('confirm_requires')}")
                elif verdict["confirm"] == mr.FALSE:
                    tier = "supported"
                    reasons.append(
                        f"HELD at supported: required markers not found "
                        f"({spec.get('confirm_requires')})")
                elif verdict["confirm"] == mr.UNKNOWN:
                    tier = "supported"
                    reasons.append(
                        "HELD at supported: marker evidence INCOMPLETE - untested "
                        f"marker(s) {', '.join(verdict['untested']) or 'unknown'}. "
                        "This is a database gap, not a negative result.")
                elif spec.get("confirm_requires") is None:
                    tier = min(tier, "supported", key=lambda t: rank_order(t))
                    reasons.append(
                        "no confirm_requires rule for this taxon - cannot reach "
                        "confirmed_in_silico by design (see marker_rules.yaml rationale)")

        if n_controls == 0 and row.get("role") != "control":
            reasons.append(
                "NO NEGATIVE CONTROL IN RUN - contamination cannot be excluded"
            )
        elif not control_rpm and row.get("role") != "control":
            # Controls exist but none of them raised anything for any taxon. That is the
            # good case and must not read like the bad one: subtraction had nothing to
            # subtract because the controls were clean, not because they were absent.
            reasons.append(
                f"{n_controls} negative control(s) in run, all clean - no taxon raised in "
                f"any control, so no control subtraction applied"
            )

        # A control's own content is QC information, not a finding about a
        # specimen. Keeping it in the actionable tiers would inflate the counts an
        # analyst triages on.
        # A NEGATIVE CONTROL'S OWN CONTENT IS NEVER AN ACTIONABLE FINDING.
        #
        # This only caught `supported` and `confirmed_in_silico`, but marker gating runs
        # BEFORE it and can set vaccine_strain or marker_veto. A control whose marker
        # evidence tripped vaccine_if therefore stayed `vaccine_strain` and was rendered
        # in the report's "Vaccine or attenuated strain" block - presented as an
        # actionable stand-down about a specimen, when it is the no-template control's
        # own background. The tier list is meant to be exhaustive for controls, so the
        # test is now "did anything assert presence", not an enumeration.
        # D17a: NO MARKER RULE BLOCK MEANS NO CONFIRMATION. 31 of the 89 watchlist
        # entries have no block in marker_rules.yaml - 29 of them `notifiable`,
        # including Dengue, Zika, West Nile, all three orthohantaviruses, Salmonella
        # enterica and C. neoformans. For those taxa the veto/vaccine/confirm block
        # above simply does not execute, and the finding reached `confirmed_in_silico`
        # on alignment geometry alone with an EMPTY marker_verdict and no caveat.
        #
        # This was dormant only because nothing could reach the top tier: until
        # min_alignment_length was made reachable, every confirmation in the project's
        # history was Y. pestis or E. coli, both of which have blocks. The first
        # real-data confirmation that fix produces is Salmonella enterica - notifiable,
        # marker_verdict empty, no block - which is why the two land together (D20).
        #
        # Capped at `supported`, matching the unevaluable-veto rule a few lines above
        # and for the same asymmetry: the finding is still reported and still
        # prominent, it is simply not stamped as confirmed on evidence nobody gathered.
        if (tier == "confirmed_in_silico" and mr is not None and marker_rules
                and not spec and not declared_no_markers):
            tier = "supported"
            reasons.append(
                "CAPPED AT SUPPORTED: no marker rule block exists for this taxon, so "
                "the veto, vaccine and confirmation checks never ran. Alignment "
                "geometry alone cannot distinguish wild type from a vaccine or "
                "attenuated strain, nor the target from a near neighbour carrying a "
                "diagnostic deletion. Add a block to marker_rules.yaml, or declare "
                "`no_markers_available: true` there to state that none exist.")

        # A ROW A CONTROL NEVER GOT TO VOTE ON IS NOT A ROW NO CONTROL TOUCHED.
        # Suppression lives in the final else of the tier chain, so a row below the read
        # floor or already taken by attributed_to_neighbour exits before the control is
        # consulted. Those rows used to carry "" - the same value as a taxon no control
        # ever raised - which made the column unable to express the very case the stderr
        # table was built to expose. Measured on lod_v3_threearm: 18 B. anthracis rows
        # sitting under a 40.09 RPM control (a standing 200.46 RPM floor) were
        # indistinguishable from 6 Y. pestis rows no control had raised at all.
        if row.get("role") != "control" and ctrl > 0 and not control_subtraction:
            control_subtraction = f"not_evaluated:{ctrl_name}@{ctrl:.2f}"

        if row.get("role") == "control" and tier in ASSERTS_PRESENCE:
            reasons.append(f"negative control background at {rpm:.2f} RPM "
                           f"(was {tier} before control reclassification)")
            tier = "control_background"

        scr = screen.get((row.get("sample", ""), taxid), {})
        sylph_support = scr.get("sylph_abundance", "")
        concordant = bool(sylph_support and as_float(sylph_support) > 0)
        n_clf = int(as_float(scr.get("n_classifiers"), 0))
        clf_list = scr.get("classifiers", "")
        distinct_kmers = scr.get("krakenuniq_distinct_kmers", "")
        screen_flags = scr.get("screen_flags", "")

        if tier in ("supported", "confirmed_in_silico"):
            # Stage-1 concordance is recorded, never required. Four classifiers are
            # run precisely because they fail in different directions; demanding
            # agreement would throw away the sensitivity that buys. But a taxon seen
            # by three independent methods is far stronger evidence than one seen by
            # a single k-mer counter, and the analyst should see which they have.
            if clf_list:
                reasons.append(
                    f"stage-1 concordance: {n_clf} classifier(s) [{clf_list}]")
            if not concordant:
                reasons.append(
                    "no sylph/ANI support (expected for viruses and eukaryotes, "
                    "which GTDB does not cover)")

            # Distinct-k-mer evidence. This is the stage-1 analogue of window
            # occupancy: reads piled on a few conserved k-mers look identical to
            # reads spread across a genome if you only count reads.
            if screen_flags and "krakenuniq_not_run" not in screen_flags:
                reasons.append(
                    f"CAVEAT: weak distinct-k-mer evidence at stage 1 ({screen_flags}) "
                    "- alignment evidence carried this finding, so treat the k-mer "
                    "disagreement as a reason to look closely")
            elif "krakenuniq_not_run" in screen_flags:
                reasons.append(
                    "CAVEAT: KrakenUniq not run, so no distinct-k-mer evidence was "
                    "available to corroborate or contradict the read counts")

        findings.append({
            "sample": row.get("sample", ""),
            "role": row.get("role", "sample"),
            # Carried through because it decides which read-identity floor was applied,
            # so a reader can tell an ONT confirmation from an Illumina one.
            "platform": plat,
            "tier": tier,
            "taxid": taxid,
            "taxon": row.get("taxon") or meta.get("name", ""),
            "group": row.get("group") or meta.get("group", ""),
            "category": meta.get("category", ""),
            "unique_reads": row.get("unique_reads", "0"),
            "reads_per_million": row.get("reads_per_million", "0"),
            "breadth": row.get("breadth", "0"),
            "window_occupancy": row.get("window_occupancy", "0"),
            "mean_depth": row.get("mean_depth", "0"),
            "depth_cv": row.get("depth_cv", "0"),
            "mean_identity": row.get("mean_identity", "0"),
            "neighbour_ratio": row.get("neighbour_ratio", ""),
            "mapped_fraction_of_assigned": row.get("mapped_fraction_of_assigned", ""),
            "reads_in_masked_regions": row.get("reads_in_masked_regions", ""),
            "masked_fraction_of_reference": row.get("masked_fraction_of_reference", ""),
            "best_neighbour_taxid": row.get("best_neighbour_taxid", ""),
            # THE NEIGHBOUR'S NAME, NOT ONLY ITS NUMBER. `attributed_to_neighbour`
            # reports a watchlist organism as absent because another organism won the
            # reads; a column holding "620" makes that unauditable without an external
            # lookup, and the two calls this was written for - Shigella out-voting
            # E. coli, C. deneoformans out-voting C. neoformans - are exactly the ones
            # a reader most needs to recognise on sight. Empty when the manifest is
            # absent or the taxid is not in it.
            "best_neighbour_name": ref_names.get(
                (row.get("best_neighbour_taxid") or "").strip(), ("", ""))[0],
            # The number that decides whether that neighbour may out-vote the
            # target travels with the finding. "Every finding explains itself" is
            # this module's stated principle, and a gate whose input is only in the
            # prose of `reasons` cannot be audited or re-thresholded from the TSV.
            "best_neighbour_uniformity": row.get("best_neighbour_uniformity", ""),
            "screen_reads": scr.get("assigned_reads", ""),
            "n_classifiers": n_clf or "",
            "classifiers": clf_list,
            "distinct_kmers": distinct_kmers,
            "screen_flags": screen_flags,
            "sylph_abundance": sylph_support,
            "orthogonal_concordance": "yes" if concordant else "no",
            "reference_bp": row.get("reference_bp", ""),
            "taxon_notes": meta.get("notes", ""),
            "marker_verdict": marker_summary,
            # Sub-species / strain evidence, and what the depth can actually carry.
            "strain_typing": strain_typing,
            "strain_claim": strain_claim,
            "strain_claim_reason": strain_reason,
            # "" | cleared:<c>@<rpm> | suppressed:<c>@<rpm> | exceeds_control:<c>@<rpm>
            # | not_evaluated:<c>@<rpm>. Empty means NO CONTROL RAISED THIS TAXON in
            # this run - distinct from `cleared:` (a control held it and the sample beat
            # the margin) and from `not_evaluated:` (a control held it, but an earlier
            # gate decided the row before subtraction could vote).
            "control_subtraction": control_subtraction,
            "reasons": " | ".join(reasons),
        })

    # Database gaps must survive into the report - a watchlist taxon that could
    # not be checked is not the same as one that was checked and cleared.
    if args.missing_refs and args.missing_refs.exists() and args.missing_refs.name != "NO_FILE":
        for line in args.missing_refs.read_text().splitlines():
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sample, taxid, taxon = parts[0], parts[1], parts[2]
            meta = watch.get(taxid, {})
            findings.append({
                # role was hardcoded "sample" here, so a CONTROL's database gap was
                # reported as a sample's. The screen row knows which it was.
                "sample": sample,
                "role": (screen.get((sample, taxid), {}).get("role") or "sample"),
                "platform": (screen.get((sample, taxid), {}).get("platform") or ""),
                "tier": "no_reference",
                "taxid": taxid, "taxon": taxon or meta.get("name", ""),
                "group": meta.get("group", ""), "category": meta.get("category", ""),
                "unique_reads": "", "reads_per_million": "", "breadth": "",
                "window_occupancy": "", "mean_depth": "", "depth_cv": "",
                "mean_identity": "", "neighbour_ratio": "",
                "mapped_fraction_of_assigned": "", "reads_in_masked_regions": "",
                "masked_fraction_of_reference": "", "best_neighbour_taxid": "",
                "best_neighbour_name": "",
                "screen_reads": screen.get((sample, taxid), {}).get("assigned_reads", ""),
                # These were left empty even though the screen row is available two
                # lines above (screen_reads uses it), so the report lost which arms had
                # raised a taxon whose reference is missing - exactly the information
                # needed to judge whether the database gap matters.
                "n_classifiers": screen.get((sample, taxid), {}).get("n_classifiers", ""),
                "classifiers": screen.get((sample, taxid), {}).get("classifiers", ""),
                "distinct_kmers": screen.get((sample, taxid), {}).get(
                    "krakenuniq_distinct_kmers", ""),
                "screen_flags": screen.get((sample, taxid), {}).get("screen_flags", ""),
                "sylph_abundance": "", "orthogonal_concordance": "",
                "reference_bp": "", "taxon_notes": meta.get("notes", ""),
                "marker_verdict": "",
                "strain_typing": "", "strain_claim": "not_attempted_no_reference",
                "strain_claim_reason": ("no reference genome, so neither depth nor "
                                        "marker typing could be assessed"),
                "control_subtraction": "",
                "reasons": "DATABASE GAP: no reference genome in pack; stage-2 "
                           "confirmation could not run. Build refs and rerun.",
            })

    rank = {t: i for i, t in enumerate(TIER_ORDER)}
    findings.sort(key=lambda f: (-rank.get(f["tier"], 0), f["sample"],
                                 -as_float(f["unique_reads"])))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(findings[0]) if findings else [
        "sample", "role", "tier", "taxid", "taxon", "group", "category",
        "control_subtraction", "reasons"
    ]
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for f in findings:
            w.writerow(f)

    counts = defaultdict(int)
    for f in findings:
        counts[f["tier"]] += 1
    sys.stderr.write(f"[score] {len(findings)} findings -> {args.out}\n")
    for tier in reversed(TIER_ORDER):
        if counts[tier]:
            sys.stderr.write(f"[score]   {tier}: {counts[tier]}\n")
    # THE SUBTRACTION AUDIT. Every taxon ANY control raised gets a line, including the
    # ones that suppressed nothing - those are the dangerous ones, because a control
    # holding a select agent at a level too low to bite today is a standing floor that
    # will bite the moment a sample arrives just under it, and today it is invisible.
    # `if control_rpm:` alone is truthy when every entry is 0.0, and the per-row
    # `ctrl <= 0` skip below then prints a header with no body.
    if any(v > 0 for v in control_rpm.values()):
        taxon_of = {r["taxid"]: r.get("taxon", "") for r in coverage}
        live = [t for t, v in control_rpm.items() if v > 0]

        # WIDTHS COME FROM THE DATA, AND TRUNCATION IS VISIBLE WHEN IT HAPPENS.
        # The first version of this table padded to constants and sliced `src` with
        # [-37:], which truncates from the LEFT and without a marker. d12fine's control
        # names are 38-40 characters, so every one of them printed as an identifier that
        # does not exist - `od_gram_negative_fine_illumina_r1_NTC` is not a sample, is
        # not in the samplesheet, and greps to nothing. In the one table whose whole
        # purpose is to name the control responsible, a silently altered identifier is
        # worse than no table: it reads as a data error rather than a display one.
        def fit(text: str, width: int) -> str:
            return text if len(text) <= width else text[:width - 1] + "\u2026"

        w_name = max(12, min(40, max(len(taxon_of.get(t, "") or f"taxid {t}")
                                     for t in live)))
        w_src = max(12, min(48, max(len(control_src.get(t, "?") or "?")
                                    for t in live)))
        sys.stderr.write(
            f"\n[score] CONTROL SUBTRACTION AUDIT (control_ratio {args.control_ratio}x)\n"
            f"[score]   {'taxon':<{w_name}}  {'control':<{w_src}}{'ctrl RPM':>12}"
            f"{'floor RPM':>12}{'suppr':>7}{'elig':>6}{'rows':>6}\n")
        for tid, ctrl in sorted(control_rpm.items(),
                                key=lambda kv: -kv[1]):
            if ctrl <= 0:
                continue
            reach = control_reach.get(tid, {"rows": 0, "eligible": 0,
                                            "suppressed": 0, "exceeds": 0})
            name = fit(taxon_of.get(tid, "") or f"taxid {tid}", w_name)
            src = fit(control_src.get(tid, "?") or "?", w_src)
            flags = []
            if reach["suppressed"] == 0:
                # Armed and unfired - but only when there was something to fire AT.
                # A taxon with no sample row in this run at all is merely unexercised,
                # and calling that LATENT would bury the dangerous case (rows present,
                # every one pre-empted) among harmless ones.
                flags.append("LATENT" if reach["rows"] else "UNEXERCISED")
            if reach["exceeds"]:
                flags.append(f"{reach['exceeds']} EXCEED-CONTROL")
            if reach["rows"] and not reach["eligible"]:
                # Every row for this taxon was taken by an earlier gate, so the
                # control never got a vote. Without this the table would imply the
                # control had been consulted and had cleared them.
                flags.append("all rows pre-empted by an earlier gate")
            sys.stderr.write(
                f"[score]   {name:<{w_name}}  {src:<{w_src}}{ctrl:>12.2f}"
                f"{args.control_ratio * ctrl:>12.2f}{reach['suppressed']:>7}"
                f"{reach['eligible']:>6}{reach['rows']:>6}"
                f"{('   <- ' + ', '.join(flags)) if flags else ''}\n")
        n_exc = sum(v["exceeds"] for v in control_reach.values())
        if n_exc:
            sys.stderr.write(
                f"[score]   {n_exc} suppressed row(s) carry MORE of the organism than "
                f"the control does and are held down only by the "
                f"{args.control_ratio}x margin.\n")

    if n_controls == 0:
        sys.stderr.write(
            "[score] WARNING: no negative control in this run. Contamination "
            "cannot be excluded for any finding.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
