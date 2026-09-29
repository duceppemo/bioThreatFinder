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


def load_pack_thresholds(pack_yaml: Path) -> dict:
    """Minimal YAML scrape for the thresholds block.

    Deliberately not a YAML parser: this script is stdlib-only so it runs in any
    conda env the pipeline happens to hand it, including ones without PyYAML.
    """
    defaults = {
        "min_unique_reads": 20.0,
        "min_breadth": 0.10,
        "max_depth_cv": 2.5,
        "min_read_identity": 0.90,
    }
    if not pack_yaml or not pack_yaml.exists():
        return defaults
    in_block = False
    for raw in pack_yaml.read_text().splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith((" ", "\t")):
            in_block = line.strip().rstrip(":") == "thresholds"
            continue
        if in_block and ":" in line:
            k, v = line.strip().split(":", 1)
            try:
                defaults[k.strip()] = float(v.strip())
            except ValueError:
                pass
    return defaults


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
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    thresholds = load_pack_thresholds(args.pack_yaml)

    # ---- marker evidence and rules -------------------------------------
    # marker_state[sample][marker_id] = mr.TRUE / mr.FALSE. A marker that is
    # UNKNOWN is deliberately ABSENT FROM THE DICT, so the rule engine's default
    # of UNKNOWN applies and indeterminacy propagates correctly.
    marker_state: dict[str, dict[str, int]] = defaultdict(dict)
    marker_seen = 0
    for path in args.marker_evidence:
        for row in read_tsv(path):
            state = (row.get("state") or "").strip().upper()
            sample = row.get("sample", "")
            mid = row.get("marker_id", "")
            if not mid or not sample or mr is None:
                continue
            marker_seen += 1
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

        min_reads = as_float(meta.get("min_unique_reads"), thresholds["min_unique_reads"]) \
            or thresholds["min_unique_reads"]
        min_breadth = as_float(meta.get("min_breadth"), thresholds["min_breadth"]) \
            or thresholds["min_breadth"]

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
                gates = {
                    f"breadth {breadth:.4f} >= {min_breadth:.4f}": breadth >= min_breadth,
                    f"occupancy {occupancy:.4f} >= {args.min_occupancy:.4f}":
                        occupancy >= args.min_occupancy,
                    f"identity {identity:.4f} >= {thresholds['min_read_identity']:.4f}":
                        identity >= thresholds["min_read_identity"],
                    f"depth_cv {depth_cv:.2f} <= {thresholds['max_depth_cv']:.2f}":
                        depth_cv <= thresholds["max_depth_cv"],
                }
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

        # ---- marker gating -------------------------------------------
        marker_summary = ""
        spec = marker_rules.get(row.get("taxon") or meta.get("name", ""))
        if spec and mr is not None and tier in ("supported", "confirmed_in_silico"):
            verdict = mr.evaluate(row.get("taxon", ""), spec, marker_state.get(row.get("sample", ""), {}))
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
            # An untested veto or vaccine check must not block a confirmation that
            # its own required markers support - but it must be DISCLOSED. A
            # confirmed finding whose BCG or vaccinia check could not run is
            # materially different from one where the check ran and passed, and
            # the analyst is the one who needs to know which they are looking at.
            if verdict["veto"] == mr.UNKNOWN:
                reasons.append(
                    f"CAVEAT: veto check could not be evaluated ({spec.get('veto_if')}) "
                    "- untested marker(s); an argument against this call may exist "
                    "but was not tested")
            if verdict["vaccine"] == mr.UNKNOWN:
                reasons.append(
                    f"CAVEAT: vaccine/attenuated-strain check could not be evaluated "
                    f"({spec.get('vaccine_if')}) - a vaccine strain CANNOT be excluded "
                    "for this finding")

            if tier == "confirmed_in_silico" or verdict["confirm"] is not None:
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
        if row.get("role") == "control" and tier in ("supported", "confirmed_in_silico"):
            reasons.append(f"negative control background at {rpm:.2f} RPM")
            tier = "control_background"

        scr = screen.get((row.get("sample", ""), taxid), {})
        sylph_support = scr.get("sylph_abundance", "")
        concordant = bool(sylph_support and as_float(sylph_support) > 0)
        if tier in ("supported", "confirmed_in_silico"):
            reasons.append(
                "orthogonal sylph support present" if concordant
                else "no orthogonal sylph support (expected for viruses and eukaryotes)"
            )

        findings.append({
            "sample": row.get("sample", ""),
            "role": row.get("role", "sample"),
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
            "best_neighbour_taxid": row.get("best_neighbour_taxid", ""),
            "screen_reads": scr.get("assigned_reads", ""),
            "sylph_abundance": sylph_support,
            "orthogonal_concordance": "yes" if concordant else "no",
            "reference_bp": row.get("reference_bp", ""),
            "taxon_notes": meta.get("notes", ""),
            "marker_verdict": marker_summary,
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
                "sample": sample, "role": "sample", "tier": "no_reference",
                "taxid": taxid, "taxon": taxon or meta.get("name", ""),
                "group": meta.get("group", ""), "category": meta.get("category", ""),
                "unique_reads": "", "reads_per_million": "", "breadth": "",
                "window_occupancy": "", "mean_depth": "", "depth_cv": "",
                "mean_identity": "", "neighbour_ratio": "", "best_neighbour_taxid": "",
                "screen_reads": screen.get((sample, taxid), {}).get("assigned_reads", ""),
                "sylph_abundance": "", "orthogonal_concordance": "",
                "reference_bp": "", "taxon_notes": meta.get("notes", ""),
                "marker_verdict": "",
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
