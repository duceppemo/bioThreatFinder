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
     "Reads aligned better to a related organism than to the watchlist target."),
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
    by_tier: dict[str, list[dict]] = defaultdict(list)
    for f in findings:
        by_tier[f.get("tier", "screened_only")].append(f)

    samples = sorted({f["sample"] for f in findings if f.get("sample")})
    has_control = any(f.get("role") == "control" for f in findings)
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

    summary_cols = ["sample", "role", "tier", "taxid", "taxon", "group", "category",
                    "unique_reads", "reads_per_million", "breadth", "window_occupancy",
                    "mean_identity", "neighbour_ratio", "n_classifiers", "classifiers",
                    "distinct_kmers", "screen_flags", "orthogonal_concordance", "reasons"]
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
