#!/usr/bin/env python3
"""
Tests for the geometric time-axis aggregate (data/pronuclei_clocks.json) that adds two
extra clocks to the Transcripts-vs-Pronuclear-Distance page:
  * mat_polar — maternal ♀ pronucleus → polar body,
  * sperm_pat — sperm entry → paternal ♂ pronucleus.

Designed to FAIL on the ways this could silently drift:
  * a split (no ♀/♂ consensus) zygote getting a maternal/paternal distance;
  * the maternal/paternal identity not matching the assignments consensus;
  * a distance not matching a direct recomputation from the source geometry;
  * a distance defined for a zygote missing the required part (polar / sperm);
  * the developmental DIRECTION flipping (both distances must rise with τ);
  * the declared per-clock counts disagreeing with the rows;
  * absolute paths leaking into the committed artifact.

Run:  python3 scripts/test_pronuclei_clocks.py
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AGG = os.path.join(HERE, "data", "pronuclei_clocks.json")
ASSIGN = os.path.join(HERE, "data", "pronuclei_assignments.json")
PT = os.path.join(HERE, "data", "pronuclei_pseudotime.json")
MANIFEST = os.path.join(HERE, "data", "pronuclei_manifest.json")

FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def d3(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return 0.0
    mx = sum(xs) / n; my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs); syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    return sxy / math.sqrt(sxx * syy) if sxx and syy else 0.0


def main():
    print("pronuclei geometric clocks — tests\n")
    if not os.path.isfile(AGG):
        print("  FAIL  aggregate missing — run python3 build_pronuclei_clocks.py")
        return 1
    d = json.load(open(AGG))
    E = d["embryos"]
    assign = {r["id"]: r for r in json.load(open(ASSIGN))["embryos"]}
    scale = d["unit_um_per_plot"]

    print("[structure]")
    check("scale is 0.15 µm/plot-unit", abs(scale - 0.15) < 1e-9, str(scale))
    keys = {c["key"] for c in d["clocks"]}
    check("declares both migration clocks", keys == {"mat_polar", "sperm_pat"}, str(keys))
    check("every clock is flagged larger_is_later", all(c["larger_is_later"] for c in d["clocks"]))
    check("declared n matches row count", d["n"] == len(E))

    print("\n[consensus & identity]")
    # a split zygote has no ♀/♂ identity → cannot have EITHER distance, so must not appear at all
    bad_split = [e["id"] for e in E
                 if (assign[e["id"]].get("consensus") or {}).get("split")]
    check("no split zygote is emitted", not bad_split, str(bad_split))
    # every row's maternal/paternal identity is that of the assignments consensus
    id_ok = all((assign[e["id"]].get("consensus") or {}).get("female") is not None for e in E)
    check("every emitted zygote has a ♀/♂ consensus", id_ok)

    print("\n[distances reproduce the source geometry]")
    bad = 0
    for e in E:
        r = assign[e["id"]]
        female = (r.get("consensus") or {}).get("female")
        pron = r["pron"]
        polar = (r.get("polar") or {}).get("com_plot")
        sperm = r.get("sperm_plot")
        mat = pron[female]["com_plot"]
        pat = pron[1 - female]["com_plot"]
        exp_mp = round(d3(mat, polar) * scale, 2) if polar else None
        exp_sp = round(d3(sperm, pat) * scale, 2) if sperm else None
        if e["mat_polar"] != exp_mp:
            bad += 1
        if e["sperm_pat"] != exp_sp:
            bad += 1
    check("recomputed distances match the stored ones", bad == 0, f"{bad} mismatch")

    print("\n[definedness follows the required parts]")
    mp_bad = [e["id"] for e in E
              if (e["mat_polar"] is not None) != bool((assign[e["id"]].get("polar") or {}).get("com_plot"))]
    sp_bad = [e["id"] for e in E
              if (e["sperm_pat"] is not None) != bool(assign[e["id"]].get("sperm_plot"))]
    check("mat_polar defined iff a polar body exists", not mp_bad, str(mp_bad[:4]))
    check("sperm_pat defined iff a sperm is labelled", not sp_bad, str(sp_bad[:4]))
    check("per-clock n matches non-null rows",
          all(c["n"] == sum(1 for e in E if e[c["key"]] is not None) for c in d["clocks"]))

    print("\n[developmental direction — both rise with τ]")
    tau = {e["id"]: e for e in json.load(open(PT))["embryos"]}
    for key, floor in (("mat_polar", 0.0), ("sperm_pat", 0.3)):
        xy = [(tau[e["id"]]["tau"], e[key]) for e in E
              if e[key] is not None and tau.get(e["id"], {}).get("tau") is not None]
        r = pearson([p[0] for p in xy], [p[1] for p in xy])
        check(f"{key} rises with τ (r > {floor})", r > floor, f"r={r:.2f} (n={len(xy)})")

    print("\n[joins to the page's manifest]")
    mids = {e["id"] for e in json.load(open(MANIFEST))["embryos"]}
    missing = [e["id"] for e in E if e["id"] not in mids]
    check("every clock zygote is in the pronuclei manifest", not missing, str(missing[:4]))

    print("\n[provenance]")
    txt = open(AGG, "rb").read()
    leaks = [s for s in (b"/Users/", b"/Volumes/", b"C:\\") if s in txt]
    check("no absolute paths in the artifact", not leaks, str(leaks))

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"  FAILED: {f}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
