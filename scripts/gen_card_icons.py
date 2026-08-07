#!/usr/bin/env python
"""Generate the small schematic SVG icons used on the landing page (index.html)
and mirrored in the repo README, one per model card.

These are NOT screenshots — they're minimal line-art diagrams of what each
model actually measures (a plane through a sphere, two dots and a distance
line, a scatter against a time axis, ...), so they read at thumbnail size and
survive the light/dark themes via `currentColor` (the card sets `color` to its
own `--accent`). Regenerate after adding/renaming a model:

    python scripts/gen_card_icons.py
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "card-icons")
os.makedirs(OUT, exist_ok=True)

W, H = 64, 44
CX, CY = 32, 22
SW = 2.4          # stroke width
DOT_R = 2.6

HEAD = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'fill="none" stroke="currentColor" stroke-width="{SW}" '
        f'stroke-linecap="round" stroke-linejoin="round">\n')
TAIL = "</svg>\n"


def dot(x, y, r=DOT_R, fill=True, opacity=None):
    op = f' fill-opacity="{opacity}"' if opacity is not None else ""
    return (f'<circle cx="{x}" cy="{y}" r="{r}" '
            f'{"fill=\"currentColor\" stroke=\"none\"" if fill else "fill=\"none\""}{op}/>\n')


def line(x1, y1, x2, y2, dash=None, opacity=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    op = f' stroke-opacity="{opacity}"' if opacity is not None else ""
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}"{d}{op}/>\n'


def circle(cx, cy, r, opacity=None, dash=None):
    op = f' stroke-opacity="{opacity}"' if opacity is not None else ""
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<circle cx="{cx}" cy="{cy}" r="{r}"{op}{d}/>\n'


def ellipse(cx, cy, rx, ry, rot=0, opacity=None, dash=None):
    op = f' stroke-opacity="{opacity}"' if opacity is not None else ""
    d = f' stroke-dasharray="{dash}"' if dash else ""
    t = f' transform="rotate({rot} {cx} {cy})"' if rot else ""
    return f'<ellipse cx="{cx}" cy="{cy}" rx="{rx}" ry="{ry}"{op}{d}{t}/>\n'


def path(d, opacity=None, dash=None, fill=None):
    op = f' stroke-opacity="{opacity}"' if opacity is not None else ""
    da = f' stroke-dasharray="{dash}"' if dash else ""
    f = f' fill="{fill}"' if fill else ""
    return f'<path d="{d}"{op}{da}{f}/>\n'


def arrow(x1, y1, x2, y2, opacity=None):
    """A line with a small open arrowhead at (x2,y2)."""
    import math
    ang = math.atan2(y2 - y1, x2 - x1)
    a1 = ang + 2.6; a2 = ang - 2.6
    hx1, hy1 = x2 + 5 * math.cos(a1), y2 + 5 * math.sin(a1)
    hx2, hy2 = x2 + 5 * math.cos(a2), y2 + 5 * math.sin(a2)
    return (line(x1, y1, x2, y2, opacity=opacity) +
            line(x2, y2, hx1, hy1, opacity=opacity) +
            line(x2, y2, hx2, hy2, opacity=opacity))


def body(cell_r=13, dashed_outline=False):
    """The base 'embryo' circle every geometry icon builds on."""
    return circle(CX, CY, cell_r, opacity=(.5 if dashed_outline else .85),
                  dash=("3 3" if dashed_outline else None))


def write(name, inner):
    with open(os.path.join(OUT, f"{name}.svg"), "w", encoding="utf-8", newline="\n") as f:
        f.write(HEAD + inner + TAIL)
    print("wrote", name + ".svg")


ICONS = {}

# ── Pseudotime: model development ──────────────────────────────────────
ICONS["pronuclear-pseudotime"] = (
    body(14) +
    dot(CX - 5, CY - 2, 4.6, fill=False) + dot(CX + 6, CY + 3, 4.2, fill=False) +
    path(f"M 10 {H - 6} Q 32 {H - 14} 54 {H - 6}", opacity=.35) +
    dot(34, H - 9, 2.2) +
    line(34, H - 9, 34, H - 15, opacity=.5) +
    path(f"M 27 {H - 5} Q 34 {H - 12} 41 {H - 5}", opacity=.22, dash="2 2")
)
ICONS["pseudotime-calibration"] = (
    dot(20, 14, 3.4) + dot(34, 20, 3.4) + line(20, 14, 34, 20) +
    line(8, H - 7, 56, H - 7, opacity=.4) +
    path(f"M 8 {H - 7} L 24 20 L 44 12 L 56 {H - 20}", opacity=.7) +
    dot(24, 20, 1.8) + dot(44, 12, 1.8)
)
ICONS["vision-pseudotime"] = (
    "".join(f'<rect x="{6+i*4}" y="{10+i*2}" width="18" height="14" rx="2" '
            f'stroke-opacity="{.9 - i*.25}"/>\n' for i in range(3)) +
    arrow(30, 20, 44, 20) +
    circle(52, 20, 7, opacity=.8, dash="2 2") +
    line(52, 15, 52, 25, opacity=.5) + line(47, 20, 57, 20, opacity=.5)
)

# ── Pseudotime: applied ────────────────────────────────────────────────
ICONS["pn3d-transcripts"] = (
    line(8, H - 7, 56, H - 7, opacity=.4) +
    dot(16, 15, 2) + dot(24, 24, 2) + dot(30, 12, 2) + dot(38, 20, 2) + dot(46, 16, 2) + dot(52, 26, 2) +
    line(20, 30, 44, 30, opacity=.55, dash="2 2") +
    line(20, 27, 20, 33, opacity=.55) + line(44, 27, 44, 33, opacity=.55)
)
ICONS["pronuclei"] = (
    body(14) + dot(CX - 5, CY - 1, 3.6, fill=False) + dot(CX + 5, CY + 2, 3.6, fill=False) +
    line(CX - 5 + 3.4, CY - 1, CX + 5 - 3.4, CY + 2, dash="2 2")
)
ICONS["extpt"] = (
    dot(12, 22, 4, fill=False) +
    dot(28, 19, 3) + dot(28, 25, 3) +
    dot(44, 16, 2.4) + dot(48, 20, 2.4) + dot(44, 25, 2.4) + dot(48, 28, 2.4) +
    arrow(18, 22, 22, 22) + arrow(34, 22, 38, 22) +
    line(6, H - 5, 58, H - 5, opacity=.35)
)

# ── Pronuclei identity ─────────────────────────────────────────────────
ICONS["pronuclei-assignments"] = (
    body(14) +
    dot(CX - 5, CY - 1, 3.6, fill=False) + dot(CX + 5, CY + 2, 3.6, fill=False) +
    dot(CX, CY - 15, 2, opacity=.7) +
    line(CX - 5, CY - 1, CX, CY - 15, opacity=.45, dash="1.5 2") +
    line(CX + 5, CY + 2, CX, CY - 15, opacity=.45, dash="1.5 2")
)

# ── Sperm: locating & aligning ─────────────────────────────────────────
ICONS["sperm-map"] = (
    body(15) + dot(CX + 6, CY - 4, 2.4) +
    circle(CX, CY, 8, opacity=.35, dash="2 2") + circle(CX, CY, 15, opacity=.2, dash="2 2") +
    dot(CX - 8, CY + 6, 1.6, opacity=.6) + dot(CX + 2, CY + 9, 1.6, opacity=.6)
)
ICONS["alignment"] = (
    "".join(circle(CX, CY, r, opacity=.16 + .1 * (3 - i), dash=None)
            for i, r in enumerate([16, 13, 10])) +
    arrow(CX, CY + 15, CX, CY - 15) + dot(CX + 5, CY - 6, 2.2)
)
ICONS["alphabeta"] = (
    dot(20, CY, 8, fill=False) + dot(44, CY, 8, fill=False) +
    arrow(28, CY - 10, 36, CY - 10) + arrow(36, CY + 10, 28, CY + 10)
)
ICONS["sperm-pca"] = (
    dot(14, 28, 1.8) + dot(20, 22, 1.8) + dot(26, 24, 1.8) + dot(30, 16, 1.8) +
    dot(38, 18, 1.8) + dot(44, 12, 1.8) + dot(48, 20, 1.8) +
    arrow(12, 30, 50, 10) +
    dot(50, 10, 2.6, opacity=.9)
)
ICONS["sperm-pseudotime"] = (
    dot(CX, CY, 2.6) +
    dot(CX - 12, CY - 6, 3, fill=False) + dot(CX + 12, CY + 4, 3, fill=False) +
    line(CX, CY, CX - 12, CY - 6, dash="1.5 2", opacity=.55) +
    line(CX, CY, CX + 12, CY + 4, dash="1.5 2", opacity=.55) +
    line(8, H - 6, 56, H - 6, opacity=.35)
)

# ── Division-plane geometry ─────────────────────────────────────────────
ICONS["zygote-planes"] = (
    circle(CX, CY, 14) +
    "".join(line(CX, CY, round(CX + 14 * __import__("math").cos(a), 1),
                 round(CY + 14 * __import__("math").sin(a), 1), opacity=.28)
            for a in [i * 3.14159 / 6 for i in range(6)])
)
ICONS["sperm-division"] = (
    circle(CX, CY, 14) + dot(CX, CY, 1.8) +
    dot(CX - 8, CY - 9, 2, opacity=.85) + dot(CX + 9, CY + 8, 2, opacity=.85) +
    line(CX - 8, CY - 9, CX + 9, CY + 8, opacity=.8)
)
ICONS["planes-all"] = (
    ellipse(CX, CY, 14, 14) +
    ellipse(CX, CY, 14, 5, opacity=.35) + ellipse(CX, CY, 14, 9, opacity=.3, rot=35) +
    ellipse(CX, CY, 14, 9, opacity=.3, rot=-35) + ellipse(CX, CY, 5, 14, opacity=.3)
)
ICONS["equatorial-planes"] = (
    circle(CX, CY, 14) + line(CX - 14, CY, CX + 14, CY, opacity=.9) +
    dot(CX, CY - 14, 1.8, opacity=.7)
)
ICONS["compare-planes"] = (
    circle(CX, CY, 14) +
    line(CX - 14, CY, CX + 14, CY, opacity=.85) +
    ellipse(CX, CY, 14, 5, opacity=.55, rot=25) +
    ellipse(CX, CY, 14, 5, opacity=.4, rot=-25) +
    ellipse(CX, CY, 5, 14, opacity=.3)
)
ICONS["sperm-sphere"] = (
    circle(CX, CY, 14) + dot(CX + 7, CY - 5, 2.2) +
    circle(CX + 7, CY - 5, 7, opacity=.45, dash="2 2")
)
ICONS["axes"] = (
    circle(CX, CY, 14) +
    arrow(CX - 10, CY + 10, CX + 10, CY - 10) +
    arrow(CX - 9, CY - 3, CX + 6, CY + 9)
)

# ── Gene expression & spatial pattern ───────────────────────────────────
ICONS["contact"] = (
    dot(20, CY, 9, fill=False) + dot(44, CY, 9, fill=False) +
    "".join(line(29, CY - 8 + i * 3.5, 35, CY - 8 + i * 3.5, opacity=.5) for i in range(6)) +
    line(32, CY - 10, 32, CY + 10, opacity=.5, dash="1 3")
)
ICONS["clustering"] = (
    dot(16, 14, 2) + dot(20, 18, 2) + dot(15, 20, 2) +
    dot(42, 12, 2) + dot(46, 17, 2) + dot(48, 10, 2) +
    dot(28, 30, 2) + dot(33, 33, 2) + dot(37, 29, 2) +
    circle(18, 17, 7, opacity=.3, dash="2 2") +
    circle(45, 13, 7, opacity=.3, dash="2 2") +
    circle(33, 31, 7, opacity=.3, dash="2 2")
)
ICONS["stage-expression"] = (
    dot(12, CY, 6, fill=False) +
    dot(27, CY - 3, 5, fill=False) + dot(27, CY + 4, 5, fill=False) +
    dot(42, CY - 4, 4, fill=False) + dot(48, CY - 1, 4, fill=False) +
    dot(42, CY + 5, 4, fill=False) + dot(48, CY + 8, 4, fill=False) +
    line(6, H - 4, 58, H - 4, opacity=.3)
)
ICONS["segments"] = (
    circle(CX, CY, 14) +
    path(f"M {CX} {CY-14} A 14 14 0 0 1 {CX+14} {CY} L {CX} {CY} Z", opacity=.5) +
    path(f"M {CX} {CY} L {CX+14} {CY} A 14 14 0 0 1 {CX+3} {CY+13.7} Z", opacity=.3) +
    line(CX, CY - 14, CX, CY, opacity=.5) + line(CX, CY, CX + 14, CY, opacity=.5) +
    line(CX, CY, CX + 3, CY + 13.7, opacity=.5)
)
ICONS["diffusion"] = (
    dot(CX, CY, 2) +
    circle(CX, CY, 7, opacity=.55, dash="2 2") +
    circle(CX, CY, 13, opacity=.3, dash="2 2") +
    dot(CX + 9, CY - 8, 1.6, opacity=.8) + dot(CX - 11, CY + 4, 1.6, opacity=.8)
)

missing = [k for k in ICONS if not k]
for name, inner in ICONS.items():
    write(name, inner)

print(f"\n{len(ICONS)} icons written to {os.path.abspath(OUT)}")
