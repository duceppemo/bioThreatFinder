#!/usr/bin/env python3
"""
make_workflow_svg.py - regenerate the README workflow diagram in assets/workflow/.

Two files, light and dark, because a README `<img>` has no page context: `currentColor`
and `prefers-color-scheme` inside an SVG served as an image cannot be relied on through
GitHub's image proxy. The repo already solves this for the logo with a `<picture>` and two
files, so the diagram does the same.

WHAT THE PICTURE IS FOR, since an ASCII version already exists in the README: it shows the
two fan-outs and the competitive bundle, which are the parts prose keeps having to explain.
Stage 1 is an OR across arms - one arm nominating a taxon is enough - and stage 2 re-asks
the question per candidate against an index holding the target AND its near neighbours, so
a read that fits a harmless relative better is attributed to the relative. That split is
what `neighbour_ratio` measures, and it is the mechanism behind decisions D2 and D12.

Geometry is computed from a grid rather than eyeballed, so the two files cannot drift and
boxes cannot creep out of alignment.

Text is NOT outlined to paths (unlike make_logo.py, which does outline its wordmark). The
labels are short and the boxes are padded generously, so a substituted font shifts glyphs
without colliding; outlining ~60 labels would bloat the file for little gain.

    python3 scripts/make_workflow_svg.py
    python3 scripts/make_workflow_svg.py --out assets/workflow
"""
from __future__ import annotations

import argparse
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Brand palette, taken from assets/logo/bioThreatFinder-logo-{light,dark}.svg.
THEMES = {
    "light": {
        "ink": "#10231F",        # primary text and strokes
        "accent": "#0F766E",     # the teal that carries meaning
        "muted": "#5B6B66",      # secondary annotation
        "rule": "#C8D3CE",       # band dividers, box outlines
        "panel": "#F4F7F5",      # neutral box fill
        "accent_fill": "#E4F0EE",  # fill for the boxes that decide
        "warn": "#8A5A00",       # [opt] markers - off by default
    },
    "dark": {
        "ink": "#EEF2EF",
        "accent": "#5EC4B6",
        "muted": "#9FB3AD",
        "rule": "#2C4A43",
        "panel": "#16302B",
        "accent_fill": "#1B3B36",
        "warn": "#D9A54A",
    },
}

MONO = "ui-monospace, SFMono-Regular, Menlo, Consolas, 'Liberation Mono', monospace"
SANS = "ui-sans-serif, system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, sans-serif"

W, H = 1140, 790


def esc(t: str) -> str:
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class Canvas:
    def __init__(self, c: dict):
        self.c = c
        self.parts: list[str] = []

    def add(self, s: str) -> None:
        self.parts.append("  " + s)

    def box(self, x, y, w, h, *, fill=None, stroke=None, dash=None, r=6, sw=1.3):
        f = fill or self.c["panel"]
        s = stroke or self.c["rule"]
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" '
                 f'fill="{f}" stroke="{s}" stroke-width="{sw}"{d}/>')

    def text(self, x, y, t, *, size=13, fill=None, anchor="middle", font=MONO,
             weight="400", spacing=None):
        ls = f' letter-spacing="{spacing}"' if spacing else ""
        self.add(f'<text x="{x}" y="{y}" font-family="{font}" font-size="{size}" '
                 f'font-weight="{weight}" fill="{fill or self.c["ink"]}" '
                 f'text-anchor="{anchor}"{ls}>{esc(t)}</text>')

    def line(self, x1, y1, x2, y2, *, stroke=None, sw=1.4, dash=None, arrow=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        a = f' marker-end="url(#{arrow})"' if arrow else ""
        self.add(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
                 f'stroke="{stroke or self.c["ink"]}" stroke-width="{sw}"{d}{a}/>')

    def path(self, d, *, stroke=None, sw=1.4, dash=None, arrow=None):
        dd = f' stroke-dasharray="{dash}"' if dash else ""
        a = f' marker-end="url(#{arrow})"' if arrow else ""
        self.add(f'<path d="{d}" fill="none" stroke="{stroke or self.c["ink"]}" '
                 f'stroke-width="{sw}"{dd}{a}/>')

    def band(self, y, label, note):
        """A divider that names a stage and its posture."""
        self.add(f'<line x1="28" y1="{y}" x2="{W - 28}" y2="{y}" '
                 f'stroke="{self.c["rule"]}" stroke-width="1"/>')
        self.text(36, y - 8, label, size=11, fill=self.c["accent"], anchor="start",
                  font=SANS, weight="700", spacing="1.4")
        self.text(W - 36, y - 8, note, size=11, fill=self.c["muted"], anchor="end",
                  font=SANS)


def build(theme: str) -> str:
    c = THEMES[theme]
    g = Canvas(c)
    cx = 470                      # main spine
    rx = 930                      # right-hand lane: side outputs and the watchlist

    # ---- input ----------------------------------------------------------------
    # The trim boxes end exactly on the spine (x = cx), so the platform branches rejoin
    # with a straight vertical rather than a dogleg that reads as a closed loop.
    g.box(40, 44, 190, 40)
    g.text(135, 69, "samplesheet.csv", size=13)

    g.box(cx - 190, 26, 190, 36, fill=c["panel"])
    g.text(cx - 95, 49, "FASTP", size=13)
    g.text(cx - 16, 49, "illumina", size=10, fill=c["muted"], anchor="end", font=SANS)
    g.box(cx - 190, 72, 190, 36, fill=c["panel"])
    g.text(cx - 95, 95, "CHOPPER", size=13)
    g.text(cx - 16, 95, "ont", size=10, fill=c["muted"], anchor="end", font=SANS)

    g.path(f"M230 64 H258 V44 H{cx - 190}", arrow="ink")
    g.path(f"M230 64 H258 V90 H{cx - 190}", arrow="ink")
    g.text(254, 20, "by platform", size=10, fill=c["muted"], anchor="start", font=SANS)
    g.line(cx, 44, cx, 90, stroke=c["ink"], sw=1.4)

    # ---- host depletion + its QC artifact -------------------------------------
    g.box(cx - 130, 140, 260, 44, fill=c["panel"])
    g.text(cx, 160, "MINIMAP2_HOST_DEPLETE", size=12)
    g.text(cx, 176, "both mates unmapped kept", size=10, fill=c["muted"], font=SANS)
    g.path(f"M{cx} 90 V140", arrow="ink")

    g.box(rx - 150, 140, 300, 44, dash="4 3")
    g.text(rx, 160, "qc/host_depletion.tsv", size=12, fill=c["accent"])
    g.text(rx, 176, "what fraction of host got through", size=10, fill=c["muted"],
           font=SANS)
    g.line(cx + 130, 162, rx - 150, 162, dash="4 3", stroke=c["muted"], arrow="muted")

    # ---- stage 1 --------------------------------------------------------------
    g.band(214, "STAGE 1  ·  WHAT MIGHT BE HERE", "permissive — missing a lead is worse")

    arms = [
        (155, "KRAKEN2", "+ BRACKEN", None),
        (365, "SYLPH_PROFILE", "ANI-based", None),
        (575, "METABULI", "DNA + AA", None),
        (785, "KRAKENUNIQ", "distinct k-mers", "opt"),
    ]
    for ax, name, sub, opt in arms:
        g.box(ax, 248, 190, 48, fill=c["panel"])
        g.text(ax + 95, 268, name, size=12)
        g.text(ax + 95, 285, sub, size=10, fill=c["muted"], font=SANS)
        if opt:
            g.text(ax + 183, 263, "[opt]", size=9, fill=c["warn"], anchor="end", font=SANS)
        g.path(f"M{ax + 95} 296 V318", stroke=c["muted"], sw=1.2)

    g.path(f"M{cx} 184 V206", arrow="ink")
    g.line(250, 318, 880, 318, stroke=c["muted"], sw=1.2)
    g.text(cx + 26, 340, "OR  ·  one arm is enough to carry a taxon forward", size=11,
           fill=c["accent"], anchor="start", font=SANS, weight="600")

    g.box(cx - 130, 352, 260, 44, fill=c["accent_fill"], stroke=c["accent"])
    g.text(cx, 372, "PARSE_SCREEN", size=12)
    g.text(cx, 388, "candidates per sample", size=10, fill=c["muted"], font=SANS)
    g.path(f"M{cx} 318 V352", arrow="ink")

    g.box(rx - 150, 352, 300, 44)
    g.text(rx, 372, "RESOLVE_WATCHLIST", size=12)
    g.text(rx, 388, "name → taxid + descendants", size=10, fill=c["muted"], font=SANS)
    g.line(rx - 150, 374, cx + 130, 374, stroke=c["muted"], arrow="muted")

    # ---- stage 2 --------------------------------------------------------------
    g.band(424, "STAGE 2  ·  DOES THE ALIGNMENT HOLD UP", "strict — this is what decides")

    # Competitive bundle, drawn open so the mechanism is visible rather than named.
    g.box(60, 458, 480, 196, fill="none", stroke=c["rule"], dash="5 4")
    g.text(80, 478, "MINIMAP2_COMPETITIVE", size=12, anchor="start")
    g.text(520, 478, "× candidate", size=10, fill=c["accent"], anchor="end", font=SANS,
           weight="600")

    g.text(80, 500, "one index, target + near neighbours:", size=10, fill=c["muted"],
           anchor="start", font=SANS)
    for i, (lbl, fill, stroke) in enumerate([
            ("target", c["accent_fill"], c["accent"]),
            ("neighbour", c["panel"], c["rule"]),
            ("neighbour", c["panel"], c["rule"])]):
        bx = 80 + i * 152
        g.box(bx, 510, 140, 30, fill=fill, stroke=stroke, r=4)
        g.text(bx + 70, 530, lbl, size=11)

    g.text(80, 566, "a read that fits a relative better is scored to the relative", size=10,
           fill=c["muted"], anchor="start", font=SANS)

    g.box(80, 580, 440, 58, fill=c["panel"])
    g.text(100, 600, "COVERAGE_STATS", size=12, anchor="start")
    g.text(100, 618, "breadth · window occupancy · depth_cv · identity", size=10,
           fill=c["muted"], anchor="start", font=SANS)
    g.text(100, 632, "neighbour_ratio  ·  masked regions excluded", size=10,
           fill=c["muted"], anchor="start", font=SANS)

    # Marker tier.
    g.box(600, 458, 480, 196, fill="none", stroke=c["rule"], dash="5 4")
    g.text(620, 478, "2b  MARKER TIER", size=12, anchor="start")
    g.text(1060, 478, "virulence loci", size=10, fill=c["muted"], anchor="end", font=SANS)

    for i, (name, kind) in enumerate([("MARKER_ALIGN_NT", "nucleotide"),
                                      ("MARKER_SEARCH_AA", "protein")]):
        g.box(620, 498 + i * 44, 300, 36, fill=c["panel"])
        g.text(640, 521 + i * 44, name, size=11, anchor="start")
        g.text(900, 521 + i * 44, kind, size=10, fill=c["muted"], anchor="end", font=SANS)
    g.box(620, 586, 300, 36, fill=c["panel"])
    g.text(640, 609, "MARKER_SEARCH_HMM", size=11, anchor="start")
    g.text(900, 609, "profile", size=10, fill=c["muted"], anchor="end", font=SANS)
    g.path("M700 578 V586", arrow="muted", stroke=c["muted"], sw=1.2)
    g.text(712, 584, "hit reads", size=9, fill=c["muted"], anchor="start", font=SANS)

    g.box(950, 498, 110, 124, fill=c["accent_fill"], stroke=c["accent"])
    g.text(1005, 540, "PRESENT", size=10, fill=c["ink"])
    g.text(1005, 562, "ABSENT", size=10, fill=c["ink"])
    g.text(1005, 584, "UNKNOWN", size=10, fill=c["accent"])
    g.text(1005, 520, "three-valued", size=9, fill=c["muted"], font=SANS)
    g.line(920, 560, 950, 560, stroke=c["muted"], sw=1.2, arrow="muted")

    # Stage 2's two halves take DIFFERENT inputs. The competitive alignment is per
    # candidate and follows PARSE_SCREEN; the marker tier runs once per SAMPLE straight off
    # the host-depleted reads (main.nf:498-507, "Run per sample, not per candidate"). An
    # earlier version drew both as fed by PARSE_SCREEN, which misstated the fan-out.
    g.path(f"M{cx} 396 V432 H300 V458", arrow="ink")
    # Routed outside RESOLVE_WATCHLIST (which ends at x=1080) rather than through it.
    # Leaves the host-deplete box from its BOTTOM edge at x=560, not along y=184 where the
    # dashed qc/host_depletion.tsv box has its border - two coincident dashed strokes made
    # it readable as leaving the QC box. The vertical sits at x=1096, inside the band rules'
    # endpoint (W-28 = 1112) rather than exactly on it.
    # Three things to miss on the way: the dashed qc/ box border (y=184), the band notes
    # (baseline y=206, right-aligned to x=1104) and the band rules (which end at W-28=1112).
    # So: leave the host box bottom at x=560 - to the RIGHT of the main flow at x=470, so it
    # never crosses it - run at y=192 ABOVE the notes, and drop at x=1124 beyond every rule
    # endpoint. Verified by rendering; two earlier routes each traded one collision for
    # another.
    g.path("M560 184 V192 H1124 V446 H1040 V458", stroke=c["muted"], dash="4 3",
           arrow="muted")
    g.text(1118, 320, "per sample,", size=9, fill=c["muted"], anchor="end", font=SANS)
    g.text(1118, 332, "from clean reads", size=9, fill=c["muted"], anchor="end", font=SANS)

    # ---- decide ---------------------------------------------------------------
    g.band(684, "DECIDE", "gates + marker rules + control subtraction")

    g.box(cx - 160, 702, 320, 40, fill=c["accent_fill"], stroke=c["accent"], sw=1.6)
    g.text(cx, 727, "SCORE_EVIDENCE", size=13)
    g.path("M300 654 V676 H" + str(cx) + " V702", arrow="ink")
    g.path("M840 654 V676 H" + str(cx), stroke=c["ink"], sw=1.4)

    g.text(cx, 762, "findings.tsv  ·  9 tiers  ·  nothing above screened_only "
                    "without alignment evidence", size=11, fill=c["muted"], font=SANS)
    g.box(rx - 60, 702, 180, 40)
    g.text(rx + 30, 727, "MAKE_REPORT", size=12)
    g.line(cx + 160, 722, rx - 60, 722, arrow="ink")

    label = ("bioThreatFinder workflow: a permissive stage-1 OR across classifier arms "
             "feeds a strict, per-candidate stage 2 that maps reads against the target "
             "together with its near neighbours, so reads preferring a harmless relative "
             "are attributed to the relative.")
    head = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
            f'width="{W}" height="{H}" role="img" aria-label="{esc(label)}">\n'
            f'  <defs>\n'
            f'    <marker id="ink" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
            f'markerHeight="6" orient="auto-start-reverse">\n'
            f'      <path d="M0 1 L9 5 L0 9 z" fill="{c["ink"]}"/>\n    </marker>\n'
            f'    <marker id="muted" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
            f'markerHeight="6" orient="auto-start-reverse">\n'
            f'      <path d="M0 1 L9 5 L0 9 z" fill="{c["muted"]}"/>\n    </marker>\n'
            f'  </defs>\n')
    return head + "\n".join(g.parts) + "\n</svg>\n"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default=os.path.join(REPO, "assets", "workflow"))
    a = p.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    for theme in THEMES:
        path = os.path.join(a.out, f"bioThreatFinder-workflow-{theme}.svg")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(build(theme))
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
