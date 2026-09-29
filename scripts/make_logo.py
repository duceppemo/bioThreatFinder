#!/usr/bin/env python3
"""
make_logo.py - regenerate the bioThreatFinder logo files in assets/logo/.

The mark is a magnifying lens over a coverage track whose bars span the whole
reference: the breadth signal stage 2 relies on, rather than an alarm symbol.

Wordmark and tagline text are outlined to SVG paths, so the files render the
same everywhere without Space Grotesk or IBM Plex Mono installed. The fonts
(both SIL OFL) are downloaded once from Google Fonts into a cache directory.

Requires: fonttools, uharfbuzz  (pip install fonttools uharfbuzz)

    python3 scripts/make_logo.py                      # SVGs into assets/logo/
    python3 scripts/make_logo.py --chrome "$(which chromium)"
                                                      # also the social-preview PNG
"""
import argparse
import os
import subprocess
import sys
import urllib.request

try:
    import uharfbuzz as hb
    from fontTools.pens.svgPathPen import SVGPathPen
    from fontTools.pens.transformPen import TransformPen
    from fontTools.ttLib import TTFont
except ImportError:
    sys.exit("make_logo.py needs fonttools and uharfbuzz: pip install fonttools uharfbuzz")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Pinned Google Fonts builds, so a re-run reproduces the committed outlines.
FONT_URLS = {
    "SpaceGrotesk-Medium.ttf": "https://fonts.gstatic.com/s/spacegrotesk/v22/V8mQoQDjQSkFtoMM3T6r8E7mF71Q-gOoraIAEj7aUUsj.ttf",
    "SpaceGrotesk-Bold.ttf": "https://fonts.gstatic.com/s/spacegrotesk/v22/V8mQoQDjQSkFtoMM3T6r8E7mF71Q-gOoraIAEj4PVksj.ttf",
    "IBMPlexMono-Regular.ttf": "https://fonts.gstatic.com/s/ibmplexmono/v20/-F63fjptAgt5VM-kVkqdyU8n5ig.ttf",
}

LIGHT = dict(ink="#10231F", accent="#0F766E", dim="#3E4F4B", bg="#F4F1EA")
DARK = dict(ink="#EEF2EF", accent="#5EC4B6", dim="#A9BAB5", bg="#0E1B19")


def fetch_fonts(font_dir):
    os.makedirs(font_dir, exist_ok=True)
    for name, url in FONT_URLS.items():
        path = os.path.join(font_dir, name)
        if not os.path.exists(path):
            print(f"downloading {name}", file=sys.stderr)
            urllib.request.urlretrieve(url, path)


class Face:
    def __init__(self, path):
        self.tt = TTFont(path)
        self.gs = self.tt.getGlyphSet()
        self.upem = self.tt["head"].unitsPerEm
        self.cap = self.tt["OS/2"].sCapHeight / self.upem
        with open(path, "rb") as f:
            self.hbfont = hb.Font(hb.Face(f.read()))

    def run(self, text, size, x, y, tracking=0.0):
        """Return (path d, advance width) for text with baseline at (x, y)."""
        buf = hb.Buffer()
        buf.add_str(text)
        buf.guess_segment_properties()
        hb.shape(self.hbfont, buf, {"kern": True, "liga": True})
        s = size / self.upem
        order = self.tt.getGlyphOrder()
        pen = SVGPathPen(self.gs, ntos=lambda v: f"{v:.2f}".rstrip("0").rstrip("."))
        cx = x
        for info, pos in zip(buf.glyph_infos, buf.glyph_positions):
            gx = cx + pos.x_offset * s
            gy = y - pos.y_offset * s
            self.gs[order[info.codepoint]].draw(TransformPen(pen, (s, 0, 0, -s, gx, gy)))
            cx += pos.x_advance * s + tracking
        return pen.getCommands(), cx - x - tracking



# Coverage bars inside the lens: (x, top y), all bottoms at y=70.
BARS = [(28, 58), (34, 50), (40, 54), (46, 44), (52, 40), (58, 47), (64, 52), (70, 57)]


def mark(c, tx=0, ty=0, scale=1.0):
    """Coverage-lens mark, drawn in a 120x120 box."""
    bars = "".join(
        f'<rect x="{x}" y="{t}" width="4" height="{70 - t}" rx="1"/>' for x, t in BARS
    )
    return (
        f'<g transform="translate({tx} {ty}) scale({scale})">'
        f'<circle cx="52" cy="52" r="38" fill="none" stroke="{c["ink"]}" stroke-width="8"/>'
        f'<line x1="80" y1="80" x2="106" y2="106" stroke="{c["ink"]}" stroke-width="12" stroke-linecap="round"/>'
        f'<line x1="26" y1="72" x2="78" y2="72" stroke="{c["ink"]}" stroke-width="2" stroke-linecap="round" opacity="0.35"/>'
        f'<g fill="{c["accent"]}">{bars}</g>'
        "</g>"
    )


def small_mark(c):
    """Simplified mark for 16-32 px: fewer, fatter bars and heavier strokes."""
    return (
        f'<circle cx="52" cy="52" r="38" fill="none" stroke="{c["ink"]}" stroke-width="11"/>'
        f'<line x1="80" y1="80" x2="106" y2="106" stroke="{c["ink"]}" stroke-width="15" stroke-linecap="round"/>'
        f'<g fill="{c["accent"]}"><rect x="31" y="51" width="9" height="19"/>'
        f'<rect x="45" y="38" width="9" height="32"/><rect x="59" y="46" width="9" height="24"/></g>'
    )


def wordmark(c, x, baseline, size):
    track = -3 * size / 104
    d1, w1 = SG500.run("bio", size, x, baseline, track)
    d2, w2 = SG700.run("ThreatFinder", size, x + w1 + track, baseline, track)
    svg = f'<path fill="{c["accent"]}" d="{d1}"/><path fill="{c["ink"]}" d="{d2}"/>'
    return svg, w1 + track + w2


def doc(w, h, body, title):
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:g}" height="{h:g}" '
        f'viewBox="0 0 {w:g} {h:g}" role="img" aria-label="{title}">'
        f"<title>{title}</title>{body}</svg>\n"
    )


def lockup(c):
    """Mark + wordmark, transparent background, for README headers."""
    m = 200
    size = 104
    text_x = m + 40
    baseline = 52 / 120 * m + SG700.cap * size / 2
    words, ww = wordmark(c, text_x, baseline, size)
    return doc(text_x + ww + 4, m, mark(c, scale=m / 120) + words, "bioThreatFinder")


def social(c):
    """1280x640 GitHub social preview card."""
    W, H, pad = 1280, 640, 96
    m, gap, size = 180, 36, 92
    tag = "metagenomic screening · in-silico confirmation"
    my = (H - m) / 2
    baseline = my + 52 / 120 * m + SG700.cap * size / 2 - 18
    _, ww = wordmark(c, 0, baseline, size)
    pad = max(pad, (W - (m + gap + ww)) / 2)
    tx = pad + m + gap
    words, ww = wordmark(c, tx, baseline, size)
    assert tx + ww <= W - pad + 1, f"wordmark overflows: {tx + ww}"
    td, _ = MONO.run(tag, 24, tx + 2, baseline + 56, 0.5)
    body = (
        f'<rect width="{W}" height="{H}" fill="{c["bg"]}"/>'
        + mark(c, pad, my, m / 120)
        + words
        + f'<path fill="{c["dim"]}" d="{td}"/>'
    )
    return doc(W, H, body, "bioThreatFinder")


def main(argv=None):
    p = argparse.ArgumentParser(description="Regenerate the bioThreatFinder logo files.")
    p.add_argument("--out", default=os.path.join(REPO, "assets", "logo"),
                   help="output directory (default: assets/logo)")
    p.add_argument("--font-dir", default=os.path.expanduser("~/.cache/bioThreatFinder/fonts"),
                   help="where the downloaded fonts are cached")
    p.add_argument("--chrome", help="Chromium/Chrome binary; if given, also render the social-preview PNG")
    args = p.parse_args(argv)

    fetch_fonts(args.font_dir)
    global SG500, SG700, MONO
    SG500 = Face(os.path.join(args.font_dir, "SpaceGrotesk-Medium.ttf"))
    SG700 = Face(os.path.join(args.font_dir, "SpaceGrotesk-Bold.ttf"))
    MONO = Face(os.path.join(args.font_dir, "IBMPlexMono-Regular.ttf"))

    os.makedirs(args.out, exist_ok=True)
    files = {
        "bioThreatFinder-mark.svg": doc(120, 120, mark(LIGHT), "bioThreatFinder"),
        "bioThreatFinder-mark-dark.svg": doc(120, 120, mark(DARK), "bioThreatFinder"),
        "bioThreatFinder-favicon.svg": doc(120, 120, small_mark(LIGHT), "bioThreatFinder"),
        "bioThreatFinder-logo-light.svg": lockup(LIGHT),
        "bioThreatFinder-logo-dark.svg": lockup(DARK),
        "bioThreatFinder-social-preview.svg": social(LIGHT),
    }
    for name, text in files.items():
        with open(os.path.join(args.out, name), "w") as f:
            f.write(text)
        print(os.path.join(args.out, name))

    if args.chrome:
        svg = os.path.abspath(os.path.join(args.out, "bioThreatFinder-social-preview.svg"))
        png = os.path.abspath(os.path.join(args.out, "bioThreatFinder-social-preview.png"))
        subprocess.run([args.chrome, "--headless", "--no-sandbox", "--disable-gpu",
                        "--hide-scrollbars", f"--screenshot={png}", "--window-size=1280,640",
                        f"file://{svg}"], check=True, capture_output=True)
        print(png)


if __name__ == "__main__":
    main()
