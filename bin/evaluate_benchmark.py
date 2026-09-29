#!/usr/bin/env python3
"""
Measure detection performance against benchmark ground truth.

This is what turns the thresholds in packs/*/pack.yaml from reasoned guesses into
measured operating points.

EXPLICIT SEMANTICS - the whole point is that these are not left implicit

A finding's tier is mapped to one of four outcomes:

  POSITIVE        supported, confirmed_in_silico
                  The tool is asserting the organism is there.
  VACCINE         vaccine_strain
                  Asserting a vaccine/attenuated strain, i.e. a stand-down.
  NEGATIVE        screened_only, attributed_to_neighbour, marker_veto,
                  suppressed_control, control_background, or NO finding at all
                  The tool is not asserting presence.
  NOT_ASSESSED    no_reference
                  A database gap. EXCLUDED from sensitivity and specificity
                  denominators and reported separately, because scoring it as a
                  negative would flatter the tool for a question it never asked.

That last exclusion matters. If a missing reference genome counted as a true
negative, adding watchlist taxa without building their references would appear to
improve specificity.

METRICS

  sensitivity      per taxon, platform and spike level: fraction of replicates
                   reaching POSITIVE
  LOD              lowest level at which every replicate is POSITIVE
                   (reported as lod_all; lod_half is the 50% point). With few
                   replicates this is a coarse estimate and is labelled as such.
  specificity      for taxa the truth says are ABSENT: rate of POSITIVE calls,
                   and separately of confirmed_in_silico calls
  tier PPV         of all calls at tier T, the fraction where truth says present.
                   This is the calibration number - if confirmed_in_silico does
                   not have markedly higher PPV than supported, the tier
                   structure is not earning anything.
  vaccine accuracy for vaccine-strain scenarios: fraction correctly tiered
                   vaccine_strain rather than reported as disease
  ablation         the same data scored with and without marker gating, to
                   quantify what the marker tier actually contributes

A specificity failure is weighted more heavily than a sensitivity failure in the
summary, because a false confirmation on a high-consequence organism triggers a
real-world response.

Exits non-zero if any scenario expectation fails, so it can gate CI.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    import btf_yaml
except ImportError:
    sys.exit("ERROR: bin/btf_yaml.py not found beside this script")

try:
    # Single source of truth for tier ordering - duplicating it here would drift.
    from score_evidence import TIER_ORDER
except ImportError:
    sys.exit("ERROR: bin/score_evidence.py not importable; tier ordering unavailable")

POSITIVE_TIERS = {"supported", "confirmed_in_silico"}
CONFIRMED_TIERS = {"confirmed_in_silico"}
VACCINE_TIERS = {"vaccine_strain"}
NOT_ASSESSED_TIERS = {"no_reference"}

POSITIVE, VACCINE, NEGATIVE, NOT_ASSESSED = "POSITIVE", "VACCINE", "NEGATIVE", "NOT_ASSESSED"

MANIFEST_REQUIRED = ["sample", "scenario_id", "platform", "truth_present", "truth_absent"]


def tier_rank(tier: str | None) -> int:
    """Absence of a finding ranks below every tier."""
    if not tier:
        return -1
    return TIER_ORDER.index(tier) if tier in TIER_ORDER else -1


def outcome_of(tier: str | None) -> str:
    if not tier:
        return NEGATIVE
    if tier in NOT_ASSESSED_TIERS:
        return NOT_ASSESSED
    if tier in VACCINE_TIERS:
        return VACCINE
    if tier in POSITIVE_TIERS:
        return POSITIVE
    return NEGATIVE


def read_tsv(path: Path) -> list[dict]:
    if not path or not path.exists():
        return []
    with path.open() as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t")) if lines else []


def semi(value: str | None) -> list[str]:
    return [v.strip() for v in (value or "").split(";") if v.strip()]


def as_float(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def load_scenarios(path: Path) -> dict[str, dict]:
    doc = btf_yaml.load(path)
    if not isinstance(doc, dict) or "scenarios" not in doc:
        sys.exit(f"ERROR: {path} has no top-level 'scenarios' list")
    out = {}
    for sc in doc["scenarios"] or []:
        sid = sc.get("id")
        if not sid:
            sys.exit(f"ERROR: {path} contains a scenario with no 'id'")
        if sid in out:
            sys.exit(f"ERROR: duplicate scenario id {sid!r} in {path}")
        out[sid] = sc
    if not out:
        sys.exit(f"ERROR: no scenarios defined in {path}")
    return out


def load_findings(paths: list[Path]) -> dict[tuple[str, str], dict]:
    """(sample, taxon) -> best finding row. Keeps the highest-ranked tier if a
    taxon somehow appears twice for one sample."""
    out: dict[tuple[str, str], dict] = {}
    for p in paths:
        for row in read_tsv(p):
            key = (row.get("sample", ""), row.get("taxon", ""))
            prev = out.get(key)
            if prev is None or tier_rank(row.get("tier")) > tier_rank(prev.get("tier")):
                out[key] = row
    return out


# ---------------------------------------------------------------------------

def evaluate(scenarios, manifest, findings, targeted: set[str] | None = None):
    """Build per-observation records: one per (sample, taxon-of-interest).

    `targeted` is the pack's watchlist. A spiked organism that is NOT on the
    watchlist cannot be reported by the pipeline at all - background commensals
    and near-neighbour confusers are spiked precisely so they are present, but
    the tool is not looking for them. Scoring those as detection failures would
    put "LOD not reached" in the report for an organism no LOD exists for. They
    are marked `present_untargeted` and excluded from sensitivity and LOD.
    """
    obs = []
    for m in manifest:
        sample = m["sample"]
        sid = m.get("scenario_id", "")
        sc = scenarios.get(sid, {})
        present = set(semi(m.get("truth_present")))
        absent = set(semi(m.get("truth_absent")))
        vaccine = set(semi(m.get("truth_vaccine")))

        # Taxa of interest: anything the truth mentions, plus anything the tool
        # called for this sample. The latter is essential - a false positive on a
        # taxon nobody listed still has to be counted.
        called = {t for (s, t) in findings if s == sample and t}
        for taxon in present | absent | vaccine | called:
            f = findings.get((sample, taxon))
            tier = f.get("tier") if f else None
            if taxon in vaccine:
                truth = "vaccine"
            elif taxon in present:
                truth = ("present" if (targeted is None or taxon in targeted)
                         else "present_untargeted")
            elif taxon in absent:
                truth = "absent"
            else:
                # Not mentioned by the truth set at all. In a scenario where the
                # spike list is exhaustive this is a false positive; treat it as
                # 'unlisted' and report separately rather than guessing.
                truth = "unlisted"
            obs.append({
                "sample": sample, "scenario_id": sid, "kind": sc.get("kind", ""),
                "platform": m.get("platform", ""), "role": m.get("role", "sample"),
                "level_rpm": m.get("level_rpm", ""),
                "replicate": m.get("replicate", ""),
                "taxon": taxon, "truth": truth, "tier": tier or "",
                "outcome": outcome_of(tier),
                "unique_reads": (f or {}).get("unique_reads", ""),
                "reasons": (f or {}).get("reasons", ""),
            })
    return obs


def sensitivity_table(obs):
    groups = defaultdict(list)
    for o in obs:
        if o["truth"] != "present" or o["role"] == "control":
            continue
        groups[(o["scenario_id"], o["taxon"], o["platform"], o["level_rpm"])].append(o)
    rows = []
    for (sid, taxon, plat, level), items in sorted(groups.items()):
        assessable = [i for i in items if i["outcome"] != NOT_ASSESSED]
        n = len(assessable)
        pos = sum(1 for i in assessable if i["outcome"] == POSITIVE)
        conf = sum(1 for i in assessable if i["tier"] in CONFIRMED_TIERS)
        rows.append({
            "scenario_id": sid, "taxon": taxon, "platform": plat, "level_rpm": level,
            "n_replicates": len(items), "n_assessable": n,
            "n_not_assessed": len(items) - n,
            "n_positive": pos, "n_confirmed": conf,
            "detection_rate": f"{pos / n:.3f}" if n else "",
            "confirm_rate": f"{conf / n:.3f}" if n else "",
        })
    return rows


def lod_table(sens_rows):
    """Lowest level at which every assessable replicate was POSITIVE."""
    by = defaultdict(list)
    for r in sens_rows:
        lvl = as_float(r["level_rpm"])
        if lvl is None or not r["n_assessable"]:
            continue
        by[(r["scenario_id"], r["taxon"], r["platform"])].append(
            (lvl, float(r["detection_rate"]), r["n_assessable"]))
    rows = []
    for (sid, taxon, plat), pts in sorted(by.items()):
        pts.sort()
        lod_all = next((l for l, rate, _ in pts if rate >= 0.999), None)
        lod_half = next((l for l, rate, _ in pts if rate >= 0.5), None)
        min_n = min(n for _, _, n in pts)
        rows.append({
            "scenario_id": sid, "taxon": taxon, "platform": plat,
            "lod_all_rpm": lod_all if lod_all is not None else "not reached",
            "lod_half_rpm": lod_half if lod_half is not None else "not reached",
            "levels_tested": ";".join(f"{l:g}" for l, _, _ in pts),
            "min_replicates": min_n,
            "precision_note": ("coarse - fewer than 5 replicates" if min_n < 5 else ""),
        })
    return rows


def specificity_table(obs):
    groups = defaultdict(list)
    for o in obs:
        if o["truth"] != "absent":
            continue
        groups[(o["scenario_id"], o["taxon"], o["platform"])].append(o)
    rows = []
    for (sid, taxon, plat), items in sorted(groups.items()):
        assessable = [i for i in items if i["outcome"] != NOT_ASSESSED]
        n = len(assessable)
        fp = sum(1 for i in assessable if i["outcome"] == POSITIVE)
        fc = sum(1 for i in assessable if i["tier"] in CONFIRMED_TIERS)
        rows.append({
            "scenario_id": sid, "taxon": taxon, "platform": plat,
            "n_samples": len(items), "n_assessable": n,
            "n_false_positive": fp, "n_false_confirmed": fc,
            "false_positive_rate": f"{fp / n:.3f}" if n else "",
            "false_confirm_rate": f"{fc / n:.3f}" if n else "",
            "worst_tier": max((i["tier"] for i in assessable),
                              key=tier_rank, default=""),
        })
    return rows


def tier_ppv_table(obs):
    """Calibration: of calls at tier T, how many were actually true?"""
    counts = defaultdict(lambda: {"n": 0, "true": 0, "false": 0, "unlisted": 0})
    for o in obs:
        if not o["tier"] or o["outcome"] == NOT_ASSESSED or o["role"] == "control":
            continue
        c = counts[o["tier"]]
        c["n"] += 1
        if o["truth"] in ("present", "vaccine"):
            c["true"] += 1
        elif o["truth"] == "absent":
            c["false"] += 1
        else:
            c["unlisted"] += 1
    rows = []
    for tier in sorted(counts, key=tier_rank, reverse=True):
        c = counts[tier]
        denom = c["true"] + c["false"]
        rows.append({
            "tier": tier, "n_calls": c["n"], "n_truth_present": c["true"],
            "n_truth_absent": c["false"], "n_unlisted": c["unlisted"],
            "ppv": f"{c['true'] / denom:.3f}" if denom else "",
        })
    return rows


def vaccine_table(obs):
    groups = defaultdict(list)
    for o in obs:
        if o["truth"] != "vaccine":
            continue
        groups[(o["scenario_id"], o["taxon"], o["platform"])].append(o)
    rows = []
    for (sid, taxon, plat), items in sorted(groups.items()):
        n = len(items)
        correct = sum(1 for i in items if i["outcome"] == VACCINE)
        as_disease = sum(1 for i in items if i["outcome"] == POSITIVE)
        rows.append({
            "scenario_id": sid, "taxon": taxon, "platform": plat, "n_samples": n,
            "n_tiered_vaccine": correct, "n_reported_as_disease": as_disease,
            "vaccine_accuracy": f"{correct / n:.3f}" if n else "",
        })
    return rows


def check_expectations(scenarios, obs):
    """Machine-checkable assertions from the scenario file."""
    by_sc = defaultdict(list)
    for o in obs:
        by_sc[(o["scenario_id"], o["taxon"])].append(o)

    results = []
    for sid, sc in scenarios.items():
        for exp in sc.get("expect") or []:
            taxon = exp.get("taxon")
            if not taxon:
                results.append({"scenario_id": sid, "taxon": "", "check": "malformed",
                                "expected": "", "observed": "", "result": "ERROR",
                                "detail": "expectation has no 'taxon'"})
                continue
            role = exp.get("sample_role")
            items = [o for o in by_sc.get((sid, taxon), [])
                     if role is None or o["role"] == role]

            if not items:
                results.append({
                    "scenario_id": sid, "taxon": taxon, "check": "coverage",
                    "expected": "at least one observation", "observed": "none",
                    "result": "SKIP",
                    "detail": "scenario produced no observations - was it simulated and run?",
                })
                continue

            if "max_tier" in exp:
                cap = exp["max_tier"]
                bad = [o for o in items if tier_rank(o["tier"]) > tier_rank(cap)]
                results.append({
                    "scenario_id": sid, "taxon": taxon, "check": "max_tier",
                    "expected": f"<= {cap}", "observed": _worst(items),
                    "result": "FAIL" if bad else "PASS",
                    "detail": _detail(bad),
                })

            if "min_tier" in exp:
                floor = exp["min_tier"]
                bad = [o for o in items if tier_rank(o["tier"]) < tier_rank(floor)]
                results.append({
                    "scenario_id": sid, "taxon": taxon, "check": "min_tier",
                    "expected": f">= {floor}", "observed": _best(items),
                    "result": "FAIL" if bad else "PASS",
                    "detail": _detail(bad),
                })

            for level, want in (exp.get("min_tier_at_rpm") or {}).items():
                lvl = as_float(level)
                at = [o for o in items if as_float(o["level_rpm"]) == lvl]
                if not at:
                    results.append({
                        "scenario_id": sid, "taxon": taxon,
                        "check": f"min_tier_at_rpm[{level}]",
                        "expected": f">= {want}", "observed": "level not run",
                        "result": "SKIP", "detail": "",
                    })
                    continue
                bad = [o for o in at if tier_rank(o["tier"]) < tier_rank(want)]
                results.append({
                    "scenario_id": sid, "taxon": taxon,
                    "check": f"min_tier_at_rpm[{level}]",
                    "expected": f">= {want}", "observed": _best(at),
                    "result": "FAIL" if bad else "PASS",
                    "detail": _detail(bad),
                })
    return results


def _worst(items):
    return max((o["tier"] or "no finding" for o in items),
               key=lambda t: tier_rank(None if t == "no finding" else t), default="")


def _best(items):
    return min((o["tier"] or "no finding" for o in items),
               key=lambda t: tier_rank(None if t == "no finding" else t), default="")


def _detail(bad):
    if not bad:
        return ""
    return "; ".join(f"{o['sample']}={o['tier'] or 'no finding'}" for o in bad[:5])


def variant_stats(obs) -> dict:
    """Headline accuracy for one pipeline configuration."""
    assessable = [o for o in obs if o["outcome"] != NOT_ASSESSED and o["role"] != "control"]
    tp = sum(1 for o in assessable if o["truth"] == "present" and o["outcome"] == POSITIVE)
    fn = sum(1 for o in assessable if o["truth"] == "present" and o["outcome"] != POSITIVE)
    fp = sum(1 for o in assessable if o["truth"] == "absent" and o["outcome"] == POSITIVE)
    fc = sum(1 for o in assessable if o["truth"] == "absent" and o["tier"] in CONFIRMED_TIERS)
    tc = sum(1 for o in assessable if o["truth"] == "present" and o["tier"] in CONFIRMED_TIERS)
    return {"true_positive": tp, "false_negative": fn, "false_positive": fp,
            "false_confirmed": fc, "true_confirmed": tc,
            "sensitivity": f"{tp / (tp + fn):.3f}" if (tp + fn) else "",
            "confirm_ppv": f"{tc / (tc + fc):.3f}" if (tc + fc) else ""}


def ablation(obs_with, obs_without):
    """What does marker gating actually buy? (A two-variant special case.)"""
    a, b = variant_stats(obs_with), variant_stats(obs_without)
    return [{"metric": k, "with_marker_gating": a[k], "without_marker_gating": b[k]}
            for k in a]


# ---------------------------------------------------------------------------
# Tool bake-off: compare pipeline configurations on IDENTICAL data
# ---------------------------------------------------------------------------

def _dur_min(v: str) -> float:
    """Nextflow durations look like '2h 31m 4s', '4.2s', '1.5m'."""
    if not v or v in ("-", "0"):
        return 0.0
    total = 0.0
    for tok in str(v).split():
        tok = tok.strip()
        try:
            if tok.endswith("ms"):
                total += float(tok[:-2]) / 60000.0
            elif tok.endswith("h"):
                total += float(tok[:-1]) * 60.0
            elif tok.endswith("m"):
                total += float(tok[:-1])
            elif tok.endswith("s"):
                total += float(tok[:-1]) / 60.0
            else:
                total += float(tok) / 60.0
        except ValueError:
            continue
    return total


def _mem_gib(v: str) -> float:
    if not v or v in ("-", "0"):
        return 0.0
    s = str(v).strip().replace(" ", "")
    mult = {"KIB": 1 / 1048576, "MIB": 1 / 1024, "GIB": 1.0, "TIB": 1024.0,
            "KB": 1 / 1048576, "MB": 1 / 1024, "GB": 1.0, "TB": 1024.0}
    for suf, m in mult.items():
        if s.upper().endswith(suf):
            try:
                return float(s[:-len(suf)]) * m
            except ValueError:
                return 0.0
    try:
        return float(s) / 1073741824.0
    except ValueError:
        return 0.0


def parse_trace(path: Path) -> dict:
    """Aggregate resource cost from a Nextflow trace file.

    Accuracy alone cannot decide between tools. A classifier that is one point
    better for 1400 GiB of index and twelve hours of build time is a bad trade, and
    a comparison that hides the cost side invites exactly that mistake. Nextflow
    already writes a trace per run, so the cost is free to collect.
    """
    if not path or not path.exists():
        return {}
    rows = read_tsv(path)
    if not rows:
        return {}

    dur_to_min = _dur_min
    mem_to_gib = _mem_gib
    per_proc = defaultdict(lambda: {"tasks": 0, "cpu_min": 0.0, "peak_gib": 0.0})
    total_cpu_min = 0.0
    peak_gib = 0.0
    for r in rows:
        name = (r.get("process") or r.get("name") or "").split(":")[-1].split(" ")[0]
        mins = dur_to_min(r.get("realtime") or r.get("duration") or "")
        gib = mem_to_gib(r.get("peak_rss") or r.get("rss") or "")
        p = per_proc[name]
        p["tasks"] += 1
        p["cpu_min"] += mins
        p["peak_gib"] = max(p["peak_gib"], gib)
        total_cpu_min += mins
        peak_gib = max(peak_gib, gib)
    return {"per_process": dict(per_proc), "total_cpu_min": total_cpu_min,
            "peak_gib": peak_gib, "n_tasks": len(rows)}


def bakeoff_table(variants: dict, traces: dict) -> list[dict]:
    """One row per configuration: accuracy AND cost, on identical input data."""
    rows = []
    for name, obs in variants.items():
        st = variant_stats(obs)
        tr = traces.get(name) or {}
        rows.append({
            "variant": name,
            "sensitivity": st["sensitivity"],
            "confirm_ppv": st["confirm_ppv"],
            "true_positive": st["true_positive"],
            "false_negative": st["false_negative"],
            "false_positive": st["false_positive"],
            "false_confirmed": st["false_confirmed"],
            "cpu_minutes": f"{tr.get('total_cpu_min', 0):.1f}" if tr else "",
            "peak_rss_gib": f"{tr.get('peak_gib', 0):.1f}" if tr else "",
            "n_tasks": tr.get("n_tasks", "") if tr else "",
        })
    # Fewest false confirmations first, then sensitivity: on a biothreat watchlist a
    # false confirmation costs more than a missed borderline lead.
    rows.sort(key=lambda r: (r["false_confirmed"], -as_float(r["sensitivity"], 0.0)))
    return rows


def bakeoff_per_taxon(variants: dict) -> list[dict]:
    """Per-taxon detection and LOD, one column per variant.

    This is where a tool choice actually shows itself: an aggregate sensitivity
    figure hides that classifier A finds the divergent filovirus and classifier B
    finds the extra Plasmodium, which is a real trade rather than a ranking.
    """
    names = list(variants)
    per: dict[tuple, dict] = {}
    for name, obs in variants.items():
        sens = sensitivity_table(obs)
        lod = {(r["taxon"], r["platform"]): r for r in lod_table(sens)}
        for r in sens:
            key = (r["taxon"], r["platform"], r["level_rpm"])
            e = per.setdefault(key, {"taxon": r["taxon"], "platform": r["platform"],
                                     "level_rpm": r["level_rpm"]})
            e[f"detect__{name}"] = r["detection_rate"]
            e[f"confirm__{name}"] = r["confirm_rate"]
        for (taxon, plat), r in lod.items():
            for e in per.values():
                if e["taxon"] == taxon and e["platform"] == plat:
                    e[f"lod__{name}"] = r["lod_all_rpm"]
    cols = ["taxon", "platform", "level_rpm"]
    for n in names:
        cols += [f"detect__{n}", f"confirm__{n}", f"lod__{n}"]
    out = []
    for e in sorted(per.values(), key=lambda x: (x["taxon"], x["platform"],
                                                 as_float(x["level_rpm"], 0))):
        out.append({c: e.get(c, "") for c in cols})
    return out


# ---------------------------------------------------------------------------

def write_tsv(path: Path, rows: list[dict], cols: list[str] | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = cols or (list(rows[0]) if rows else ["note"])
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_summary(path: Path, tables: dict, exp_results: list[dict], abl):
    fails = [r for r in exp_results if r["result"] == "FAIL"]
    skips = [r for r in exp_results if r["result"] == "SKIP"]
    passes = [r for r in exp_results if r["result"] == "PASS"]

    L = ["# bioThreatFinder benchmark results", ""]
    L.append(f"Expectations: **{len(passes)} passed, {len(fails)} failed, "
             f"{len(skips)} skipped**")
    L.append("")
    L.append("> Simulated reads flatter every tool. These numbers establish relative "
             "performance and catch regressions; they do NOT establish clinical "
             "performance. See docs/benchmarking.md.")
    L.append("")

    if fails:
        L.append("## Failed expectations")
        L.append("")
        L.append("| scenario | taxon | check | expected | observed | detail |")
        L.append("|---|---|---|---|---|---|")
        for r in fails:
            L.append(f"| {r['scenario_id']} | {r['taxon']} | {r['check']} | "
                     f"{r['expected']} | {r['observed']} | {r['detail']} |")
        L.append("")

    spec = tables.get("specificity", [])
    bad_spec = [r for r in spec if r["n_false_confirmed"]]
    L.append("## Specificity (weighted most heavily)")
    L.append("")
    if bad_spec:
        L.append("**FALSE CONFIRMATIONS on taxa the truth says are absent:**")
        L.append("")
        L.append("| scenario | taxon | platform | false confirms / n |")
        L.append("|---|---|---|---|")
        for r in bad_spec:
            L.append(f"| {r['scenario_id']} | {r['taxon']} | {r['platform']} | "
                     f"{r['n_false_confirmed']}/{r['n_assessable']} |")
    else:
        L.append("No false confirmations on absent taxa.")
    L.append("")

    ppv = tables.get("tier_ppv", [])
    if ppv:
        L.append("## Tier calibration")
        L.append("")
        L.append("If `confirmed_in_silico` does not show markedly higher PPV than "
                 "`supported`, the tier structure is not earning its complexity.")
        L.append("")
        L.append("| tier | calls | truth present | truth absent | PPV |")
        L.append("|---|---|---|---|---|")
        for r in ppv:
            L.append(f"| {r['tier']} | {r['n_calls']} | {r['n_truth_present']} | "
                     f"{r['n_truth_absent']} | {r['ppv']} |")
        L.append("")

    lod = tables.get("lod", [])
    if lod:
        L.append("## Limit of detection")
        L.append("")
        L.append("| scenario | taxon | platform | LOD (all reps) | LOD (50%) | note |")
        L.append("|---|---|---|---|---|---|")
        for r in lod:
            L.append(f"| {r['scenario_id']} | {r['taxon']} | {r['platform']} | "
                     f"{r['lod_all_rpm']} | {r['lod_half_rpm']} | {r['precision_note']} |")
        L.append("")

    vac = tables.get("vaccine", [])
    if vac:
        L.append("## Vaccine-strain discrimination")
        L.append("")
        L.append("| scenario | taxon | platform | tiered vaccine | reported as disease |")
        L.append("|---|---|---|---|---|")
        for r in vac:
            L.append(f"| {r['scenario_id']} | {r['taxon']} | {r['platform']} | "
                     f"{r['n_tiered_vaccine']}/{r['n_samples']} | "
                     f"{r['n_reported_as_disease']} |")
        L.append("")

    bo = tables.get("bakeoff", [])
    if bo:
        L.append("## Tool bake-off")
        L.append("")
        L.append("All variants run over the SAME simulated data, so only the "
                 "configuration differs. Ordered by fewest false confirmations, then "
                 "sensitivity: on a biothreat watchlist a false confirmation costs "
                 "more than a missed borderline lead.")
        L.append("")
        L.append("Read the cost columns as part of the decision, not a footnote. A "
                 "configuration that gains a point of sensitivity for several hundred "
                 "GiB of index and hours of build time has to justify that.")
        L.append("")
        L.append("| variant | sens | confirm PPV | TP | FN | FP | **false confirms** | CPU-min | peak RSS GiB |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for r in bo:
            L.append(f"| {r['variant']} | {r['sensitivity']} | {r['confirm_ppv']} | "
                     f"{r['true_positive']} | {r['false_negative']} | "
                     f"{r['false_positive']} | **{r['false_confirmed']}** | "
                     f"{r['cpu_minutes']} | {r['peak_rss_gib']} |")
        L.append("")
        L.append("Per-taxon detection and LOD per variant: `bakeoff_per_taxon.tsv`. "
                 "Check it before ranking - an aggregate figure hides that one "
                 "configuration finds the divergent filovirus while another finds the "
                 "extra Plasmodium, which is a trade rather than a ranking.")
        L.append("")

    if abl:
        L.append("## Marker-gating ablation")
        L.append("")
        L.append("| metric | with gating | without gating |")
        L.append("|---|---|---|")
        for r in abl:
            L.append(f"| {r['metric']} | {r['with_marker_gating']} | "
                     f"{r['without_marker_gating']} |")
        L.append("")

    untargeted = sum(1 for o in tables.get("_obs", []) or []
                     if o.get("truth") == "present_untargeted")
    if untargeted:
        L.append(f"## Untargeted spikes\n\n{untargeted} observation(s) are organisms "
                 "spiked as background or as near-neighbour confusers that are NOT on "
                 "the watchlist. The pipeline is not looking for them, so they are "
                 "excluded from sensitivity and LOD rather than counted as misses.\n")

    na = sum(r["n_not_assessed"] for r in tables.get("sensitivity", []))
    if na:
        L.append(f"## Database gaps\n\n{na} observation(s) could not be assessed "
                 "(`no_reference`). These are excluded from all rates above - a "
                 "missing reference is not a negative result.\n")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenarios", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path,
                    help="benchmark manifest.tsv from simulate_benchmark.py")
    ap.add_argument("--findings", required=True, nargs="+", type=Path,
                    help="findings.tsv from one or more pipeline runs")
    ap.add_argument("--ablation-findings", nargs="*", type=Path, default=[],
                    help="findings.tsv from the same data scored WITHOUT marker gating")
    ap.add_argument("--variant", action="append", default=[], metavar="NAME=FINDINGS[,FINDINGS]",
                    help="a pipeline configuration to compare, e.g. "
                         "--variant kraken2_only=runA/findings.tsv "
                         "--variant plus_krakenuniq=runB/findings.tsv . Repeatable. All "
                         "variants MUST be runs over the SAME simulated data, or the "
                         "comparison is meaningless.")
    ap.add_argument("--variant-trace", action="append", default=[], metavar="NAME=TRACEFILE",
                    help="Nextflow trace file for a variant, so the comparison reports "
                         "cost alongside accuracy. A tool that is one point better for "
                         "1400 GiB and 12 h is a bad trade.")
    ap.add_argument("--watchlist", type=Path,
                    help="packs/<pack>/watchlist.resolved.tsv - spiked organisms not on "
                         "the watchlist are excluded from sensitivity/LOD, since the "
                         "pipeline is not looking for them")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--fail-on-skip", action="store_true",
                    help="treat SKIPped expectations as failures (for CI)")
    ap.add_argument("--self-test", action="store_true")
    args, _ = ap.parse_known_args()

    if args.self_test:
        return self_test()

    scenarios = load_scenarios(args.scenarios)
    manifest = read_tsv(args.manifest)
    if not manifest:
        sys.exit(f"ERROR: empty or unreadable manifest {args.manifest}")
    missing = [c for c in MANIFEST_REQUIRED if c not in manifest[0]]
    if missing:
        sys.exit(f"ERROR: manifest is missing required column(s): {', '.join(missing)}")

    targeted = None
    if args.watchlist:
        wl = read_tsv(args.watchlist)
        if not wl:
            sys.exit(f"ERROR: empty or unreadable watchlist {args.watchlist}")
        targeted = {r["name"] for r in wl if r.get("name")}
        sys.stderr.write(f"[bench] {len(targeted)} targeted taxa from {args.watchlist}\n")
    else:
        sys.stderr.write(
            "[bench] WARNING: no --watchlist given; spiked background and confuser "
            "organisms will be scored as detection targets and may appear as "
            "spurious LOD failures\n")

    findings = load_findings(args.findings)
    sys.stderr.write(
        f"[bench] {len(scenarios)} scenarios, {len(manifest)} benchmark samples, "
        f"{len(findings)} findings\n")

    unknown = {m["scenario_id"] for m in manifest} - set(scenarios)
    if unknown:
        sys.stderr.write(
            f"[bench] WARNING: manifest references {len(unknown)} scenario id(s) not in "
            f"the scenario file: {', '.join(sorted(unknown))}\n")

    # ---- tool bake-off -------------------------------------------------
    def split_named(items, label):
        out = {}
        for item in items:
            if "=" not in item:
                sys.exit(f"ERROR: --{label} expects NAME=PATH, got {item!r}")
            name, paths = item.split("=", 1)
            name = name.strip()
            if not name:
                sys.exit(f"ERROR: --{label} has an empty name in {item!r}")
            if name in out:
                sys.exit(f"ERROR: duplicate --{label} name {name!r}")
            out[name] = [Path(p) for p in paths.split(",") if p.strip()]
        return out

    variant_paths = split_named(args.variant, "variant")
    trace_paths = split_named(args.variant_trace, "variant-trace")

    obs = evaluate(scenarios, manifest, findings, targeted)
    tables = {
        "sensitivity": sensitivity_table(obs),
        "specificity": specificity_table(obs),
        "tier_ppv": tier_ppv_table(obs),
        "vaccine": vaccine_table(obs),
    }
    tables["lod"] = lod_table(tables["sensitivity"])
    exp_results = check_expectations(scenarios, obs)

    variants = {}
    traces = {}
    if variant_paths:
        for name, paths in variant_paths.items():
            missing = [str(p) for p in paths if not p.exists()]
            if missing:
                sys.exit(f"ERROR: variant {name!r} findings not found: {', '.join(missing)}")
            variants[name] = evaluate(scenarios, manifest, load_findings(paths), targeted)
        for name, paths in trace_paths.items():
            if name not in variants:
                sys.stderr.write(
                    f"[bench] WARNING: --variant-trace {name!r} has no matching "
                    "--variant; its cost figures will be ignored\n")
                continue
            traces[name] = parse_trace(paths[0])
        sys.stderr.write(
            f"[bench] bake-off across {len(variants)} variant(s): "
            f"{', '.join(variants)}\n")

    abl = None
    if args.ablation_findings:
        obs_wo = evaluate(scenarios, manifest,
                          load_findings(args.ablation_findings), targeted)
        abl = ablation(obs, obs_wo)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out_dir / "observations.tsv", obs)
    for name, rows in tables.items():
        write_tsv(args.out_dir / f"{name}.tsv", rows)
    write_tsv(args.out_dir / "expectations.tsv", exp_results)
    if abl:
        write_tsv(args.out_dir / "ablation.tsv", abl)
    if variants:
        bo = bakeoff_table(variants, traces)
        write_tsv(args.out_dir / "bakeoff.tsv", bo)
        write_tsv(args.out_dir / "bakeoff_per_taxon.tsv", bakeoff_per_taxon(variants))
        tables["bakeoff"] = bo

    tables["_obs"] = obs
    write_summary(args.out_dir / "summary.md", tables, exp_results, abl)
    tables.pop("_obs", None)

    fails = [r for r in exp_results if r["result"] == "FAIL"]
    skips = [r for r in exp_results if r["result"] == "SKIP"]
    false_confirms = sum(r["n_false_confirmed"] for r in tables["specificity"])

    sys.stderr.write(
        f"[bench] expectations: {sum(1 for r in exp_results if r['result'] == 'PASS')} pass, "
        f"{len(fails)} fail, {len(skips)} skip\n")
    if false_confirms:
        sys.stderr.write(
            f"[bench] *** {false_confirms} FALSE CONFIRMATION(S) on absent taxa ***\n")
    for r in fails:
        sys.stderr.write(f"[bench] FAIL {r['scenario_id']}/{r['taxon']} {r['check']}: "
                         f"expected {r['expected']}, observed {r['observed']}\n")
    sys.stderr.write(f"[bench] wrote {args.out_dir}/summary.md\n")

    if fails:
        return 1
    if skips and args.fail_on_skip:
        return 1
    return 0


# ---------------------------------------------------------------------------

def self_test() -> int:
    fails = []

    scenarios = {
        "spec_x": {"id": "spec_x", "kind": "specificity",
                   "expect": [{"taxon": "Bacillus anthracis",
                               "max_tier": "attributed_to_neighbour"}]},
        "lod_y": {"id": "lod_y", "kind": "lod_series",
                  "expect": [{"taxon": "Yersinia pestis",
                              "min_tier_at_rpm": {50: "confirmed_in_silico"}}]},
        "vac_z": {"id": "vac_z", "kind": "vaccine",
                  "expect": [{"taxon": "Mycobacterium tuberculosis",
                              "max_tier": "vaccine_strain"}]},
    }

    manifest = [
        # Specificity: anthracis absent. One replicate wrongly confirms.
        {"sample": "s1", "scenario_id": "spec_x", "platform": "illumina",
         "replicate": "1", "level_rpm": "500", "role": "sample",
         "truth_present": "Bacillus cereus", "truth_absent": "Bacillus anthracis",
         "truth_vaccine": ""},
        {"sample": "s2", "scenario_id": "spec_x", "platform": "illumina",
         "replicate": "2", "level_rpm": "500", "role": "sample",
         "truth_present": "Bacillus cereus", "truth_absent": "Bacillus anthracis",
         "truth_vaccine": ""},
        # LOD: two replicates at 50 RPM, one confirms and one only supports.
        {"sample": "p1", "scenario_id": "lod_y", "platform": "illumina",
         "replicate": "1", "level_rpm": "50", "role": "sample",
         "truth_present": "Yersinia pestis", "truth_absent": "", "truth_vaccine": ""},
        {"sample": "p2", "scenario_id": "lod_y", "platform": "illumina",
         "replicate": "2", "level_rpm": "50", "role": "sample",
         "truth_present": "Yersinia pestis", "truth_absent": "", "truth_vaccine": ""},
        # Vaccine scenario, correctly tiered.
        {"sample": "v1", "scenario_id": "vac_z", "platform": "illumina",
         "replicate": "1", "level_rpm": "200", "role": "sample",
         "truth_present": "", "truth_absent": "",
         "truth_vaccine": "Mycobacterium tuberculosis"},
        # A no_reference observation that must NOT count as a true negative.
        {"sample": "n1", "scenario_id": "spec_x", "platform": "illumina",
         "replicate": "3", "level_rpm": "500", "role": "sample",
         "truth_present": "", "truth_absent": "Variola virus", "truth_vaccine": ""},
    ]

    findings = {
        ("s1", "Bacillus anthracis"): {"sample": "s1", "taxon": "Bacillus anthracis",
                                       "tier": "confirmed_in_silico"},
        ("s2", "Bacillus anthracis"): {"sample": "s2", "taxon": "Bacillus anthracis",
                                       "tier": "attributed_to_neighbour"},
        ("p1", "Yersinia pestis"): {"sample": "p1", "taxon": "Yersinia pestis",
                                    "tier": "confirmed_in_silico"},
        ("p2", "Yersinia pestis"): {"sample": "p2", "taxon": "Yersinia pestis",
                                    "tier": "supported"},
        ("v1", "Mycobacterium tuberculosis"): {"sample": "v1",
                                               "taxon": "Mycobacterium tuberculosis",
                                               "tier": "vaccine_strain"},
        ("n1", "Variola virus"): {"sample": "n1", "taxon": "Variola virus",
                                  "tier": "no_reference"},
    }

    obs = evaluate(scenarios, manifest, findings)

    # Untargeted spikes must not be scored as detection targets.
    obs_t = evaluate(scenarios, manifest, findings, targeted={"Yersinia pestis",
                                                             "Bacillus anthracis",
                                                             "Mycobacterium tuberculosis",
                                                             "Variola virus"})
    untgt = [o for o in obs_t if o["truth"] == "present_untargeted"]
    sens_t = sensitivity_table(obs_t)

    def check(cond, label):
        if not cond:
            fails.append(label)

    # Outcome mapping
    check(outcome_of("confirmed_in_silico") == POSITIVE, "confirmed -> POSITIVE")
    check(outcome_of("supported") == POSITIVE, "supported -> POSITIVE")
    check(outcome_of("vaccine_strain") == VACCINE, "vaccine -> VACCINE")
    check(outcome_of("marker_veto") == NEGATIVE, "veto -> NEGATIVE")
    check(outcome_of("suppressed_control") == NEGATIVE, "suppressed -> NEGATIVE")
    check(outcome_of("no_reference") == NOT_ASSESSED, "no_reference -> NOT_ASSESSED")
    check(outcome_of(None) == NEGATIVE, "no finding -> NEGATIVE")
    check(tier_rank(None) == -1, "no finding ranks below all tiers")

    spec = specificity_table(obs)
    anth = next((r for r in spec if r["taxon"] == "Bacillus anthracis"), None)
    check(anth is not None, "anthracis specificity row exists")
    if anth:
        check(anth["n_false_positive"] == 1, f"1 false positive, got {anth}")
        check(anth["n_false_confirmed"] == 1, f"1 false confirm, got {anth}")
        check(anth["false_positive_rate"] == "0.500", f"FP rate 0.5, got {anth}")

    # no_reference must be excluded from the denominator, not scored as a pass
    var = next((r for r in spec if r["taxon"] == "Variola virus"), None)
    check(var is not None and var["n_assessable"] == 0 and var["n_samples"] == 1,
          f"no_reference excluded from specificity denominator, got {var}")

    sens = sensitivity_table(obs)
    yp = next((r for r in sens if r["taxon"] == "Yersinia pestis"), None)
    check(yp is not None, "pestis sensitivity row exists")
    if yp:
        check(yp["n_positive"] == 2 and yp["detection_rate"] == "1.000",
              f"both replicates POSITIVE, got {yp}")
        check(yp["n_confirmed"] == 1 and yp["confirm_rate"] == "0.500",
              f"one confirmed, got {yp}")

    lod = lod_table(sens)
    ypl = next((r for r in lod if r["taxon"] == "Yersinia pestis"), None)
    check(ypl is not None and ypl["lod_all_rpm"] == 50.0,
          f"LOD all-replicates 50, got {ypl}")
    check(ypl is not None and "coarse" in ypl["precision_note"],
          "few-replicate precision caveat present")

    ppv = {r["tier"]: r for r in tier_ppv_table(obs)}
    check("confirmed_in_silico" in ppv, "confirmed tier in PPV table")
    if "confirmed_in_silico" in ppv:
        c = ppv["confirmed_in_silico"]
        check(c["n_truth_present"] == 1 and c["n_truth_absent"] == 1
              and c["ppv"] == "0.500", f"confirmed PPV 0.5, got {c}")

    vac = vaccine_table(obs)
    check(len(vac) == 1 and vac[0]["n_tiered_vaccine"] == 1
          and vac[0]["n_reported_as_disease"] == 0,
          f"vaccine correctly tiered, got {vac}")

    exp = check_expectations(scenarios, obs)
    by = {(r["scenario_id"], r["check"]): r for r in exp}
    check(by[("spec_x", "max_tier")]["result"] == "FAIL",
          "false confirmation must FAIL the max_tier expectation")
    check(by[("lod_y", "min_tier_at_rpm[50]")]["result"] == "FAIL",
          "one replicate below the floor must FAIL min_tier_at_rpm")
    check(by[("vac_z", "max_tier")]["result"] == "PASS",
          "correct vaccine tiering must PASS")

    # Ablation arithmetic
    obs_wo = evaluate(scenarios, manifest, {
        **findings,
        ("s2", "Bacillus anthracis"): {"sample": "s2", "taxon": "Bacillus anthracis",
                                       "tier": "confirmed_in_silico"},
    })
    abl = {r["metric"]: r for r in ablation(obs, obs_wo)}
    check(abl["false_confirmed"]["with_marker_gating"] == 1
          and abl["false_confirmed"]["without_marker_gating"] == 2,
          f"ablation shows gating prevents a false confirm, got {abl['false_confirmed']}")

    check(any(o["taxon"] == "Bacillus cereus" for o in untgt),
          "B. cereus spiked but not on the watchlist -> present_untargeted")
    check(not any(r["taxon"] == "Bacillus cereus" for r in sens_t),
          "untargeted taxon must NOT appear in the sensitivity/LOD tables")
    check(any(r["taxon"] == "Yersinia pestis" for r in sens_t),
          "targeted taxon still appears in the sensitivity table")

    # Bake-off: two configurations over identical data, one with an extra false
    # confirmation, plus cost figures.
    v_good = evaluate(scenarios, manifest, findings)
    v_bad = evaluate(scenarios, manifest, {
        **findings,
        ("s2", "Bacillus anthracis"): {"sample": "s2", "taxon": "Bacillus anthracis",
                                       "tier": "confirmed_in_silico"},
    })
    bo = {r["variant"]: r for r in bakeoff_table({"good": v_good, "bad": v_bad}, {})}
    check(bo["good"]["false_confirmed"] == 1 and bo["bad"]["false_confirmed"] == 2,
          f"bake-off counts false confirms per variant, got {bo}")
    order = [r["variant"] for r in bakeoff_table({"bad": v_bad, "good": v_good}, {})]
    check(order[0] == "good",
          f"bake-off must rank fewest false confirmations first, got {order}")
    pt = bakeoff_per_taxon({"good": v_good, "bad": v_bad})
    check(any("detect__good" in r and "detect__bad" in r for r in pt),
          "per-taxon bake-off emits a column per variant")

    # Duration and memory parsing, which decide the cost half of the comparison.
    check(abs(_dur_min("2h 30m") - 150.0) < 1e-6, f"2h 30m -> 150 min, got {_dur_min('2h 30m')}")
    check(abs(_dur_min("90s") - 1.5) < 1e-6, f"90s -> 1.5 min, got {_dur_min('90s')}")
    check(abs(_mem_gib("2 GB") - 2.0) < 1e-6, f"2 GB -> 2.0, got {_mem_gib('2 GB')}")
    check(abs(_mem_gib("512 MB") - 0.5) < 1e-6, f"512 MB -> 0.5, got {_mem_gib('512 MB')}")
    check(_mem_gib("-") == 0.0 and _dur_min("-") == 0.0, "missing trace values are 0")

    if fails:
        sys.stderr.write(f"evaluate_benchmark self-test FAILED ({len(fails)}):\n")
        for f in fails:
            sys.stderr.write(f"  {f}\n")
        return 1
    sys.stderr.write("evaluate_benchmark self-test: all checks passed\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
