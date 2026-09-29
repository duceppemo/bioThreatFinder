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
                        correct answer for environmental relatives.
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

TIER_ORDER = [
    "no_reference",
    "control_background",
    "screened_only",
    "attributed_to_neighbour",
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coverage", nargs="*", type=Path, default=[],
                    help="per-(sample,taxon) TSVs from coverage_stats.py")
    ap.add_argument("--candidates", nargs="*", type=Path, default=[],
                    help="per-sample stage-1 candidate TSVs from parse_screen.py")
    ap.add_argument("--missing-refs", type=Path, default=None)
    ap.add_argument("--watchlist", required=True, type=Path,
                    help="watchlist.resolved.tsv")
    ap.add_argument("--pack-yaml", type=Path, default=None)
    ap.add_argument("--min-occupancy", type=float, default=0.15,
                    help="minimum window occupancy for confirmation [0.15]")
    ap.add_argument("--min-neighbour-ratio", type=float, default=1.5,
                    help="target reads must exceed best neighbour by this factor [1.5]")
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

    thresholds = apply_threshold_overrides(load_pack_thresholds(args.pack_yaml),
                                           args.threshold)

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
    for path in args.candidates:
        for row in read_tsv(path):
            key = (row.get("sample", ""), row.get("taxid", ""))
            screen[key] = row

    coverage: list[dict] = []
    for path in args.coverage:
        coverage.extend(read_tsv(path))

    # Highest control signal per taxon, used to subtract run-wide contamination.
    control_rpm: dict[str, float] = defaultdict(float)
    control_src: dict[str, str] = {}
    for row in coverage:
        if row.get("role") == "control":
            rpm = as_float(row.get("reads_per_million"))
            if rpm > control_rpm[row["taxid"]]:
                control_rpm[row["taxid"]] = rpm
                control_src[row["taxid"]] = row.get("sample", "")
    n_controls = len({r["sample"] for r in coverage if r.get("role") == "control"})

    findings: list[dict] = []

    for row in coverage:
        taxid = row["taxid"]
        meta = watch.get(taxid, {})

        # `as_float(x, default)` already falls back when the field is absent or
        # unparseable. The trailing `or thresholds[...]` that used to be here also
        # swallowed a LEGITIMATE ZERO: a watchlist row asking for
        # min_unique_reads: 0 ("report every read for this select agent") had its
        # request replaced by the pack default, and the reason string then cited a
        # threshold the operator never set.
        min_reads = as_float(meta.get("min_unique_reads"),
                             thresholds["min_unique_reads"])
        min_breadth = as_float(meta.get("min_breadth"), thresholds["min_breadth"])

        unique_reads = as_float(row.get("unique_reads"))
        breadth = as_float(row.get("breadth"))
        occupancy = as_float(row.get("window_occupancy"))
        depth_cv = as_float(row.get("depth_cv"))
        identity = as_float(row.get("mean_identity"))
        nb_ratio = as_float(row.get("neighbour_ratio"), math.inf)
        rpm = as_float(row.get("reads_per_million"))

        reasons: list[str] = []
        tier = "screened_only"

        if unique_reads < min_reads:
            reasons.append(f"unique_reads {unique_reads:.0f} < {min_reads:.0f}")
        elif nb_ratio < args.min_neighbour_ratio:
            tier = "attributed_to_neighbour"
            reasons.append(
                f"reads favour neighbour taxid {row.get('best_neighbour_taxid') or '?'} "
                f"(ratio {row.get('neighbour_ratio')} < {args.min_neighbour_ratio})"
            )
        else:
            ctrl = control_rpm.get(taxid, 0.0)
            if row.get("role") != "control" and ctrl > 0 and rpm < args.control_ratio * ctrl:
                tier = "suppressed_control"
                reasons.append(
                    f"control {control_src.get(taxid, '?')} has {ctrl:.2f} RPM; "
                    f"sample {rpm:.2f} RPM does not exceed {args.control_ratio}x"
                )
            else:
                # Identity floor is group-aware: 98% is the near-neighbour
                # discriminator for bacteria (Wright et al. found ~96% meant
                # genus-correct, species-wrong), but a divergent virus is
                # legitimately far from any available reference and the same floor
                # would manufacture false negatives where sensitivity matters most.
                grp = (row.get("group") or meta.get("group") or "").strip().lower()
                plat = (row.get("platform") or "").strip().lower()
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
                    aln_floor = as_float(thresholds.get("min_alignment_length"), 100.0)
                    gates[f"mean_alignment_length {mean_aln:.0f} >= {aln_floor:.0f}"] = \
                        mean_aln >= aln_floor

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
        if spec and mr is not None and tier in ("supported", "confirmed_in_silico"):
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
        if row.get("role") == "control" and tier in (
                "supported", "confirmed_in_silico", "vaccine_strain", "marker_veto"):
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
                "reasons": "DATABASE GAP: no reference genome in pack; stage-2 "
                           "confirmation could not run. Build refs and rerun.",
            })

    rank = {t: i for i, t in enumerate(TIER_ORDER)}
    findings.sort(key=lambda f: (-rank.get(f["tier"], 0), f["sample"],
                                 -as_float(f["unique_reads"])))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(findings[0]) if findings else [
        "sample", "role", "tier", "taxid", "taxon", "group", "category", "reasons"
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
    if n_controls == 0:
        sys.stderr.write(
            "[score] WARNING: no negative control in this run. Contamination "
            "cannot be excluded for any finding.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
