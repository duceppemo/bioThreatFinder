#!/usr/bin/env python3
"""
Render findings as a self-contained HTML report plus a flat summary TSV.

Reporting choices that are deliberate:

* Findings are grouped by evidence tier, highest first, and every row carries the
  ``reasons`` string. An analyst should never have to guess why something landed
  where it did.
* ``no_reference`` findings are shown in their own block near the top. A
  watchlist organism that could not be checked is not a negative result, and
  burying it would misrepresent the run's coverage.
* Taxon-specific caveats from the watchlist (near-neighbour confusion, vaccine
  strains, species complexes) are printed with the finding, not in a footnote.
* Absence of a negative control is a run-level warning, not a silent default.

stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import html
import sys
from collections import defaultdict
from pathlib import Path

TIERS = [
    ("confirmed_in_silico", "Confirmed in silico",
     "Met every alignment threshold. Requires orthogonal laboratory confirmation."),
    ("supported", "Supported",
     "Real alignment evidence present but below the confirmation bar."),
    ("vaccine_strain", "Vaccine or attenuated strain",
     "Marker evidence indicates a vaccine or attenuated strain rather than wild "
     "type. Actionable as a stand-down, not as disease."),
    ("no_reference", "Not assessed - database gap",
     "On the watchlist and raised by the screen, but no reference genome was "
     "available, so confirmation could not run. This is NOT a negative result."),
    ("suppressed_control", "Suppressed by control",
     "Detected at a level comparable to a negative control. Treated as reagent "
     "or index contamination."),
    ("marker_veto", "Vetoed by marker evidence",
     "A marker arguing against the target was found - a neighbour-specific gene, "
     "or a region the target is known to lack."),
    ("attributed_to_neighbour", "Attributed to a near neighbour",
     "Reads aligned better to a related organism than to the watchlist target. "
     "The organism they were attributed to is named in 'Why this tier'."),
    ("indistinguishable_complex", "Cannot be separated from a near neighbour",
     "The target covers its own reference AND so does a near neighbour, which "
     "out-counts it only slightly. One of the two IS present; this evidence does not "
     "say which. Resolve with lineage markers or a targeted assay - do not read this "
     "as either a detection or a negative."),
    ("group_lead", "Lead - a member of this group, species unresolved",
     "The watchlist entry is a taxonomic group, not an organism, and the evidence "
     "says a member of it is present without saying which. This exists so a "
     "divergent member with no reference of its own still surfaces instead of "
     "vanishing. It is a lead to follow up, not a species call - and not a negative "
     "either."),
    ("control_background", "Negative-control background",
     "Organisms detected in the negative control itself. QC information about "
     "the run, not findings about any specimen."),
    ("screened_only", "Screen lead only",
     "Raised by k-mer screening but not supported by alignment. The expected "
     "outcome for classifier noise."),
]

TIER_TONE = {
    "confirmed_in_silico": "crit",
    "supported": "warn",
    "no_reference": "gap",
    "suppressed_control": "muted",
    "attributed_to_neighbour": "muted",
    # Not "muted": an unresolved complex needs an analyst to act on it, which is the
    # whole reason it exists as a tier rather than a silent attribution.
    "indistinguishable_complex": "warn",
    # Not "muted", for the same reason: a lead that nobody looks at is the Bundibugyo
    # miss happening again one step later.
    "group_lead": "warn",
    "screened_only": "muted",
    "control_background": "muted",
    "marker_veto": "muted",
    "vaccine_strain": "gap",
}

CSS = """
:root{--bg:#f7f7f8;--panel:#fff;--ink:#16181d;--dim:#5b6070;--line:#e3e4e8;
--crit:#b42318;--critbg:#fef3f2;--warn:#b54708;--warnbg:#fffaeb;
--gap:#5925dc;--gapbg:#f4f3ff;--ok:#067647;--okbg:#ecfdf3;--mutedbg:#f2f3f5;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
--bg:#0f1115;--panel:#16181d;--ink:#e7e8ec;--dim:#9aa0ae;--line:#2a2d35;
--critbg:#2d1513;--crit:#f97066;--warnbg:#2b1c07;--warn:#fdb022;
--gapbg:#1e1733;--gap:#b9a6ff;--okbg:#0d2419;--ok:#47cd89;--mutedbg:#1c1f26;}}
*{box-sizing:border-box}
body{margin:0;padding:0 16px 64px;background:var(--bg);color:var(--ink);
font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1180px;margin:0 auto}
h1{font-size:1.5rem;margin:32px 0 4px;letter-spacing:-.01em}
h2{font-size:1.12rem;margin:36px 0 2px;letter-spacing:-.01em}
.sub{color:var(--dim);font-size:.9rem;margin:0 0 8px}
.banner{background:var(--warnbg);border:1px solid var(--warn);border-left-width:4px;
border-radius:8px;padding:12px 14px;margin:18px 0;font-size:.9rem}
.banner.crit{background:var(--critbg);border-color:var(--crit)}
.meta{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px;margin:18px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:10px 12px}
.card .k{color:var(--dim);font-size:.74rem;text-transform:uppercase;letter-spacing:.04em}
.card .v{font-size:1.02rem;font-weight:600;margin-top:2px;word-break:break-word}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0 4px}
.tile{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:12px;text-align:center}
.tile .n{font-size:1.7rem;font-weight:650;line-height:1.1}
.tile .l{color:var(--dim);font-size:.76rem;margin-top:3px}
.tile.crit{background:var(--critbg);border-color:var(--crit)}.tile.crit .n{color:var(--crit)}
.tile.warn{background:var(--warnbg);border-color:var(--warn)}.tile.warn .n{color:var(--warn)}
.tile.gap{background:var(--gapbg);border-color:var(--gap)}.tile.gap .n{color:var(--gap)}
.scroll{overflow-x:auto;background:var(--panel);border:1px solid var(--line);border-radius:8px;margin-top:10px}
table{border-collapse:collapse;width:100%;font-size:.86rem}
th,td{text-align:left;padding:8px 10px;border-bottom:1px solid var(--line);white-space:nowrap}
th{background:var(--mutedbg);font-weight:600;font-size:.76rem;text-transform:uppercase;
letter-spacing:.03em;color:var(--dim);position:sticky;top:0}
tr:last-child td{border-bottom:none}
td.taxon{white-space:normal;min-width:190px;font-weight:600}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.why{white-space:normal;min-width:300px;color:var(--dim);font-size:.8rem}
.note{white-space:normal;min-width:240px;color:var(--dim);font-size:.78rem;font-style:italic}
.pill{display:inline-block;padding:1px 7px;border-radius:99px;font-size:.72rem;font-weight:600}
.pill.crit{background:var(--critbg);color:var(--crit)}
.pill.warn{background:var(--warnbg);color:var(--warn)}
.pill.gap{background:var(--gapbg);color:var(--gap)}
.pill.muted{background:var(--mutedbg);color:var(--dim)}
.empty{padding:14px;color:var(--dim);font-size:.88rem;background:var(--panel);
border:1px dashed var(--line);border-radius:8px;margin-top:10px}
footer{margin-top:44px;padding-top:14px;border-top:1px solid var(--line);
color:var(--dim);font-size:.78rem}
code{background:var(--mutedbg);padding:1px 5px;border-radius:4px;font-size:.82em;
word-break:break-all}
"""


def read_tsv(path: Path) -> list[dict]:
    if not path or not path.exists() or path.name == "NO_FILE":
        return []
    with path.open() as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t")) if lines else []


def esc(v) -> str:
    return html.escape("" if v is None else str(v))


def fmt(v, places=3) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return esc(v) if v else "-"
    return f"{f:.{places}f}"


def pack_disclaimer(pack_yaml: Path) -> str:
    """Pull the disclaimer out of pack.yaml without a YAML dependency."""
    default = ("In-silico screening result. Requires orthogonal laboratory "
               "confirmation and interpretation by a qualified laboratorian "
               "before any clinical or public-health action.")
    if not pack_yaml or not pack_yaml.exists():
        return default
    lines = pack_yaml.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("disclaimer:"):
            after = line.split("disclaimer:", 1)[1].strip()
            if after and after not in (">", "|", ">-", "|-"):
                return after.strip("'\"")
            body = []
            for nxt in lines[i + 1:]:
                if not nxt.strip():
                    break
                if not nxt.startswith(("    ", "\t")):
                    break
                body.append(nxt.strip())
            if body:
                return " ".join(body)
    return default


def render_table(rows: list[dict]) -> str:
    cols = [
        ("sample", "Sample", ""), ("taxon", "Taxon", "taxon"),
        ("group", "Group", ""), ("category", "Category", ""),
        ("unique_reads", "Reads", "num"), ("reads_per_million", "RPM", "num"),
        ("breadth", "Breadth", "num"), ("window_occupancy", "Occupancy", "num"),
        ("mean_depth", "Depth", "num"), ("depth_cv", "Depth CV", "num"),
        ("mean_identity", "Identity", "num"), ("neighbour_ratio", "Nb ratio", "num"),
        ("n_classifiers", "Clf", "num"),
        ("classifiers", "Stage-1 arms", ""),
        ("distinct_kmers", "Distinct k-mers", "num"),
        ("reads_in_masked_regions", "Masked reads", "num"),
        ("marker_verdict", "Marker rules", "why"),
        ("reasons", "Why this tier", "why"),
        ("taxon_notes", "Taxon caveats", "note"),
    ]
    head = "".join(f"<th>{esc(label)}</th>" for _, label, _ in cols)
    body = []
    for r in rows:
        tds = []
        for key, _, cls in cols:
            val = r.get(key, "")
            if key in ("breadth", "window_occupancy", "mean_identity"):
                val = fmt(val, 4)
            elif key in ("mean_depth", "depth_cv"):
                val = fmt(val, 2)
            elif key == "reads_per_million":
                val = fmt(val, 2)
            elif key == "taxon":
                val = f"{esc(r.get('taxon',''))} <span class='pill muted'>{esc(r.get('taxid',''))}</span>"
                tds.append(f"<td class='{cls}'>{val}</td>")
                continue
            tds.append(f"<td class='{cls}'>{esc(val) if val not in ('', None) else '-'}</td>")
        body.append("<tr>" + "".join(tds) + "</tr>")
    return (f"<div class='scroll'><table><thead><tr>{head}</tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--findings", required=True, type=Path)
    ap.add_argument("--read-stats", nargs="*", type=Path, default=[])
    ap.add_argument("--screen-summary", nargs="*", type=Path, default=[],
                    help="per-sample stage-1 coverage from parse_screen.py: the "
                         "unclassified fraction per arm, and any arm that was enabled "
                         "and produced nothing")
    ap.add_argument("--samplesheet", type=Path,
                    help="authoritative source of sample roles; see has_control")
    ap.add_argument("--pack-yaml", type=Path, default=None)
    ap.add_argument("--pack", default="unknown")
    ap.add_argument("--lab", default="Unspecified laboratory")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--pipeline-version", default="")
    ap.add_argument("--command", default="")
    ap.add_argument("--html", required=True, type=Path)
    ap.add_argument("--summary", required=True, type=Path)
    args = ap.parse_args()

    findings = read_tsv(args.findings)

    # STAGE-1 COVERAGE: what was searched, not what was found.
    #
    # Two things live only in this file, and until it was wired into this script both
    # reached the published screen/ directory and no further:
    #   * the unclassified fraction per arm. Perseus (Bioinformatics 2026) measured
    #     that sequence from a family absent from the database is not left
    #     unclassified but CONFIDENTLY ASSIGNED to a known taxon, approaching a 50%
    #     false-assignment rate on long reads. The novel agent is then not missed, it
    #     is mislabelled as background and dismissed. A high unclassified fraction is
    #     the only signal that survives that, so it belongs on the page an analyst
    #     reads rather than in a TSV they would have to know to open.
    #   * an arm that was ENABLED and produced nothing. METABULI, KRAKENUNIQ and
    #     SYLPH run under `errorStrategy 'ignore'` so that one dead arm cannot sink
    #     the run; the cost is that its failure looks exactly like a clean negative.
    #     An enabled arm that produced no report is a HOLE IN THE SCREEN'S COVERAGE,
    #     and a report that does not say so is claiming a search it did not perform.
    #
    # THIS STEP MUST NOT BE ABLE TO END THE RUN, AND IT IS NEW ENOUGH TO.
    #
    # MAKE_REPORT sets no errorStrategy, so an exception here kills the pipeline
    # AFTER stage 1, stage 2, the marker tier and scoring have all finished -
    # trading the actual result for the view of it. That is the exact failure
    # bf9744d fixed in PAVIAN_EXPORT, and adding a new input to the last step in
    # the run re-opens it: a truncated or half-written summary makes read_tsv
    # raise, and a non-numeric unclassified_fraction makes float() raise.
    #
    # So the whole thing degrades to a visible note. A report that cannot show the
    # coverage table is worth having; a run that produces no report because the
    # coverage table was malformed is not.
    screen_rows: list[dict] = []
    screen_error = ""
    for sp in args.screen_summary:
        try:
            screen_rows.extend(read_tsv(sp))
        except Exception as exc:                                  # noqa: BLE001
            screen_error = f"{getattr(sp, 'name', sp)}: {exc}"
            sys.stderr.write(f"[report] WARNING: cannot read screen summary {sp}: "
                             f"{exc}; stage-1 coverage will not be shown\n")
    dead_arms = [r for r in screen_rows
                 if r.get("flag") in ("ARM_PRODUCED_NOTHING", "ARM_EMPTY_REPORT")]
    high_unc = [r for r in screen_rows if r.get("flag") == "HIGH_UNCLASSIFIED"]

    control_samples: set[str] = set()
    if args.samplesheet and args.samplesheet.exists():
        with args.samplesheet.open(newline="") as fh:
            for r in csv.DictReader(fh):
                if (r.get("role") or "").strip().lower() == "control" and r.get("sample"):
                    control_samples.add(r["sample"].strip())
    by_tier: dict[str, list[dict]] = defaultdict(list)
    for f in findings:
        by_tier[f.get("tier", "screened_only")].append(f)

    samples = sorted({f["sample"] for f in findings if f.get("sample")})
    # findings come only from coverage rows, so a control that raised nothing is
    # invisible here too and the banner an analyst reads first - "No negative control
    # in this run" - fired precisely when the control was spotless. The samplesheet is
    # authoritative; the findings check stays as a fallback when it is not supplied.
    has_control = (bool(control_samples)
                   or any(f.get("role") == "control" for f in findings))
    n_crit = len(by_tier["confirmed_in_silico"])
    n_warn = len(by_tier["supported"])
    n_gap = len(by_tier["no_reference"])

    generated = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    disclaimer = pack_disclaimer(args.pack_yaml)

    parts = [f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>bioThreatFinder Report</title><style>{CSS}</style></head><body><div class="wrap">
<h1>bioThreatFinder screening report</h1>
<p class="sub">{esc(args.pack)} pack &middot; {len(samples)} sample(s) &middot; generated {esc(generated)}</p>
<div class="banner crit"><strong>Not a diagnostic result.</strong> {esc(disclaimer)}</div>"""]

    if not has_control:
        parts.append(
            "<div class='banner'><strong>No negative control in this run.</strong> "
            "Reagent, index-hopping and carry-over contamination cannot be excluded "
            "for any finding below. Include a no-template control in future runs.</div>"
        )

    if dead_arms:
        by_arm: dict[str, set] = defaultdict(set)
        for r in dead_arms:
            by_arm[r.get("classifier", "?")].add(r.get("sample", "?"))
        detail = "; ".join(f"{a} on {len(ss)} sample(s)"
                           for a, ss in sorted(by_arm.items()))
        parts.append(
            "<div class='banner'><strong>A stage-1 classifier was enabled and "
            f"produced nothing: {esc(detail)}.</strong> These arms run so that one "
            "failure cannot sink the run, which means a failure here is silent. "
            "Any organism that only the missing arm would have raised is absent "
            "from this report, and its absence is NOT a negative result. See "
            "&lsquo;Stage-1 screen coverage&rsquo; below and the screen logs.</div>")

    if high_unc:
        def _frac(v) -> float:
            try:
                return float(v or 0)
            except (TypeError, ValueError):
                return 0.0

        worst = max(_frac(r.get("unclassified_fraction")) for r in high_unc)
        n_s = len({r.get("sample") for r in high_unc})
        parts.append(
            f"<div class='banner'><strong>{n_s} sample(s) left a high fraction of "
            f"reads unclassified (up to {worst:.0%}).</strong> Sequence whose family "
            "is absent from the database is not reliably left unclassified - it is "
            "frequently assigned, with confidence, to the nearest known taxon, so a "
            "novel agent appears as familiar background rather than as a gap. The "
            "unclassified fraction is the signal that survives that, and it should "
            "be triaged (assembled, protein-searched) rather than discarded.</div>")

    parts.append(f"""<div class="meta">
<div class="card"><div class="k">Laboratory</div><div class="v">{esc(args.lab)}</div></div>
<div class="card"><div class="k">Run ID</div><div class="v">{esc(args.run_id) or '-'}</div></div>
<div class="card"><div class="k">Domain pack</div><div class="v">{esc(args.pack)}</div></div>
<div class="card"><div class="k">Pipeline</div><div class="v">v{esc(args.pipeline_version)}</div></div>
</div>
<div class="tiles">
<div class="tile crit"><div class="n">{n_crit}</div><div class="l">Confirmed in silico</div></div>
<div class="tile warn"><div class="n">{n_warn}</div><div class="l">Supported</div></div>
<div class="tile gap"><div class="n">{n_gap}</div><div class="l">Not assessed (DB gap)</div></div>
<div class="tile"><div class="n">{len(findings)}</div><div class="l">Total findings</div></div>
</div>""")

    for key, label, blurb in TIERS:
        rows = by_tier.get(key, [])
        tone = TIER_TONE[key]
        parts.append(f"<h2>{esc(label)} <span class='pill {tone}'>{len(rows)}</span></h2>")
        parts.append(f"<p class='sub'>{esc(blurb)}</p>")
        parts.append(render_table(rows) if rows
                     else "<div class='empty'>None.</div>")

    # WHAT CONTROL SUBTRACTION DID, IN THE REPORT RATHER THAN IN SCROLLBACK.
    # This module's own principle is that a caveat belongs on the finding and in the
    # HTML, not in terminal output nobody keeps. Control subtraction was the exception:
    # it silently moved 108 rows across six batches with no record outside the prose of
    # `reasons`. Two cases get their own emphasis because neither is what a reader
    # assumes "suppressed by control" means:
    #   - exceeds_control: the SAMPLE held more of the organism than the control did,
    #     and only the control_ratio margin held it down.
    #   - a control carrying a watchlist organism that suppressed nothing yet. That is
    #     a standing floor, not a clean control, and it has never been displayed.
    subtracted = [f for f in findings if (f.get("control_subtraction") or "").strip()]
    if subtracted:
        by_kind: dict = {}
        for f in subtracted:
            kind = (f.get("control_subtraction") or "").split(":", 1)[0]
            by_kind.setdefault(kind, []).append(f)
        n_exc = len(by_kind.get("exceeds_control", []))
        # exceeds_control rows ARE suppressed rows - the kind only records why. Counting
        # the two buckets separately reported "28 suppressed" on a batch that suppressed
        # 36, which is wrong on its own terms even though the next sentence named the
        # other 8.
        n_suppressed = len(by_kind.get("suppressed", [])) + n_exc
        n_notev = len(by_kind.get("not_evaluated", []))
        parts.append("<h2>Control subtraction <span class='pill "
                     f"{'warn' if n_exc else 'muted'}'>{len(subtracted)}</span></h2>")
        blurb = (f"{n_suppressed} finding(s) suppressed"
                 + (f" ({n_exc} of them carrying more of the organism than the control "
                    f"does)" if n_exc else "")
                 + f", {len(by_kind.get('cleared', []))} cleared the margin.")
        if n_exc:
            blurb += (f" Those {n_exc} are held down only by the control ratio, not by "
                      f"the control having more - review before treating them as "
                      f"negative.")
        if n_notev:
            blurb += (f" A further {n_notev} finding(s) sit under a control floor but "
                      f"were decided by an earlier gate, so control subtraction never "
                      f"voted on them.")
        parts.append(f"<p class='sub'>{esc(blurb)}</p>")
        keys = ["sample", "taxon", "tier", "reads_per_million", "control_subtraction"]
        head = "".join(f"<th>{esc(k)}</th>" for k in keys)
        body = "".join(
            "<tr>" + "".join(
                f"<td class=\"{'num' if k == 'reads_per_million' else ''}\">"
                f"{esc(r.get(k, '-'))}</td>" for k in keys) + "</tr>"
            for r in sorted(subtracted,
                            key=lambda r: (not (r.get("control_subtraction") or "")
                                           .startswith("exceeds_control"),
                                           r.get("taxon", ""))))
        parts.append(f"<div class='scroll'><table><thead><tr>{head}</tr></thead>"
                     f"<tbody>{body}</tbody></table></div>")

    # WHAT STAGE 1 SEARCHED, PER SAMPLE AND PER ARM.
    #
    # Every enabled arm gets a row even when it behaved, because the question this
    # table answers is "was this searched?", and an answer that is only printed when
    # something went wrong cannot be distinguished from a missing table. Rows that
    # cost the run coverage are sorted to the top; `not_run` rows stay, because an
    # arm the operator disabled is also a limit on what this report can claim.
    if screen_error and not screen_rows:
        parts.append(
            "<h2>Stage-1 screen coverage <span class='pill warn'>?</span></h2>"
            "<p class='sub'>The per-sample screen summary could not be read, so "
            "this report cannot say which classifier arms actually ran. Do not read "
            "the absence of this table as a clean screen &mdash; check "
            f"<code>screen/</code> directly. ({esc(screen_error)})</p>")
    if screen_rows:
        tone = "warn" if dead_arms else "muted"
        parts.append("<h2>Stage-1 screen coverage <span class='pill "
                     f"{tone}'>{len(dead_arms)}</span></h2>")
        blurb = ("One row per sample per classifier arm. `no_output` means the arm "
                 "was enabled and produced nothing - treat it as not searched, not "
                 "as negative. The unclassified fraction is reported only for the "
                 "arms that classify reads themselves.")
        parts.append(f"<p class='sub'>{esc(blurb)}</p>")
        keys = ["sample", "role", "platform", "classifier", "status", "total_reads",
                "unclassified_reads", "unclassified_fraction", "flag", "note"]
        keys = [k for k in keys if k in screen_rows[0]]
        head = "".join(f"<th>{esc(k.replace('_', ' '))}</th>" for k in keys)
        num = {"total_reads", "unclassified_reads", "unclassified_fraction"}

        def _rank(r: dict) -> tuple:
            return ({"ARM_PRODUCED_NOTHING": 0, "ARM_EMPTY_REPORT": 1,
                     "HIGH_UNCLASSIFIED": 2}.get(r.get("flag", ""), 3),
                    r.get("sample", ""), r.get("classifier", ""))

        body = ""
        for r in sorted(screen_rows, key=_rank):
            cells = ""
            for k in keys:
                cls = "num" if k in num else ("note" if k == "note" else "")
                v = r.get(k, "")
                if k == "flag" and v:
                    pill = "crit" if v != "HIGH_UNCLASSIFIED" else "warn"
                    cells += f"<td><span class='pill {pill}'>{esc(v)}</span></td>"
                else:
                    cells += f"<td class=\"{cls}\">{esc(v) if v != '' else '-'}</td>"
            body += f"<tr>{cells}</tr>"
        parts.append(f"<div class='scroll'><table><thead><tr>{head}</tr></thead>"
                     f"<tbody>{body}</tbody></table></div>")

    stats_rows = []
    for p in args.read_stats:
        stats_rows.extend(read_tsv(p))
    if stats_rows:
        parts.append("<h2>Read statistics after QC and host depletion</h2>")
        parts.append("<p class='sub'>Low post-depletion depth limits sensitivity; "
                     "interpret negatives in light of these numbers.</p>")
        keys = [k for k in ("sample", "role", "platform", "num_seqs", "sum_len",
                            "avg_len", "N50", "Q20(%)", "Q30(%)") if k in stats_rows[0]]
        head = "".join(f"<th>{esc(k)}</th>" for k in keys)
        body = "".join(
            "<tr>" + "".join(f"<td class=\"{'num' if k not in ('sample','role','platform') else ''}\">"
                             f"{esc(r.get(k,'-'))}</td>" for k in keys) + "</tr>"
            for r in stats_rows)
        parts.append(f"<div class='scroll'><table><thead><tr>{head}</tr></thead>"
                     f"<tbody>{body}</tbody></table></div>")

    parts.append(f"""<footer>
<p><strong>Interpretation limits.</strong> Sequence detection does not establish
viable organism, active infection, or exposure. Metagenomic sensitivity depends on
input depth, host background and database completeness; a taxon absent from this
report may still be present below the limit of detection. Species-level calls
within closely related complexes are frequently not defensible from metagenomic
data alone - read the taxon caveats column.</p>
<p><strong>Reproducibility.</strong> <code>{esc(args.command)}</code></p>
<p>bioThreatFinder v{esc(args.pipeline_version)} &middot; {esc(generated)}</p>
</footer></div></body></html>""")

    args.html.parent.mkdir(parents=True, exist_ok=True)
    args.html.write_text("\n".join(parts))

    # marker_verdict and taxon_notes were omitted here while extrasaction="ignore"
    # hid the omission, so a downstream consumer of the summary TSV could not see WHY a
    # finding was vetoed or held - and could not see the report_always markers
    # (V. cholerae serogroup, A. fumigatus cyp51A resistance) that ride in
    # marker_verdict. The HTML showed both; the machine-readable output did not.
    # platform is included because it decides which identity floor was applied.
    summary_cols = ["sample", "role", "platform", "tier", "taxid", "taxon", "group",
                    "category", "unique_reads", "reads_per_million", "breadth",
                    "window_occupancy", "mean_identity", "neighbour_ratio",
                    "n_classifiers", "classifiers", "distinct_kmers", "screen_flags",
                    "reads_in_masked_regions", "masked_fraction_of_reference",
                    "orthogonal_concordance", "marker_verdict",
                    # What control subtraction did to this row, machine-readable. The
                    # prose in `reasons` had to be regexed to audit it, and a regex over
                    # a message is what breaks silently when the message is reworded.
                    "control_subtraction",
                    # Sub-species typing and what the depth can carry. A reader must be
                    # able to see that genome-wide typing was NOT attempted, rather than
                    # assume it was checked and came back negative.
                    "strain_claim", "strain_typing", "strain_claim_reason",
                    "taxon_notes", "reasons"]
    with args.summary.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=summary_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for f in findings:
            w.writerow(f)

    sys.stderr.write(
        f"[report] {len(findings)} findings, {n_crit} confirmed / {n_warn} supported / "
        f"{n_gap} db-gap -> {args.html}\n"
    )
    if not has_control:
        sys.stderr.write("[report] WARNING: no negative control in this run\n")
    for r in sorted(dead_arms, key=lambda r: (r.get("sample", ""),
                                              r.get("classifier", ""))):
        sys.stderr.write(
            f"[report] WARNING: {r.get('sample')}: the {r.get('classifier')} arm was "
            f"enabled and produced {'an empty report' if r.get('flag') == 'ARM_EMPTY_REPORT' else 'nothing'}. "
            "Anything only that arm would have raised is missing from this report "
            "and is not a negative.\n")
    if high_unc:
        sys.stderr.write(
            f"[report] WARNING: {len({r.get('sample') for r in high_unc})} sample(s) "
            "flagged HIGH_UNCLASSIFIED - triage the unclassified fraction\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
