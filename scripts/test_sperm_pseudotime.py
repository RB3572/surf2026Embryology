#!/usr/bin/env python3
"""
Tests for the Sperm Location vs Pseudotime aggregate (data/sperm_pseudotime.json).

Designed to FAIL on the ways this could silently drift:
  * an embryo without a labelled sperm leaking in;
  * maternal/paternal distances defined for a split (no-consensus) zygote;
  * the sperm not sitting closest to the paternal pronucleus (the geometry sanity check);
  * the sperm→paternal-vs-τ trend not reproducing;
  * positions or distances not matching a direct recomputation from the source files;
  * absolute paths in the committed artifact.

Run: python3 scripts/test_sperm_pseudotime.py
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AGG = os.path.join(HERE, "data", "sperm_pseudotime.json")
ASSIGN = os.path.join(HERE, "data", "pronuclei_assignments.json")
PT = os.path.join(HERE, "data", "pronuclei_pseudotime.json")

FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def pearson(xs, ys):
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs); syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    if sxx == 0 or syy == 0:
        return 0.0, 0.0
    r = sxy / math.sqrt(sxx * syy)
    return r, sxy / sxx


def main():
    print("sperm location vs pseudotime — tests\n")
    if not os.path.isfile(AGG):
        print("  FAIL  aggregate missing — run python3 build_sperm_pseudotime.py")
        return 1
    d = json.load(open(AGG))
    E = d["embryos"]

    print("[sample set]")
    assign = {r["id"]: r for r in json.load(open(ASSIGN))["embryos"]}
    pt = {e["id"]: e for e in json.load(open(PT))["embryos"]}
    check("every embryo has a labelled sperm", all(assign[e["id"]].get("sperm_plot") for e in E))
    check("every embryo has a pseudotime τ", all(pt[e["id"]].get("tau") is not None for e in E))
    exp_n = sum(1 for r in assign.values() if r.get("sperm_plot") and pt.get(r["id"], {}).get("tau") is not None)
    check("count matches the sperm∩τ set", len(E) == exp_n, f"{len(E)} vs {exp_n}")
    check("declared n matches", d["n"] == len(E))

    print("\n[maternal/paternal consensus]")
    splits = [e for e in E if e["split"]]
    check("split zygotes have no maternal/paternal distance",
          all(e["dist_um"]["maternal"] is None and e["dist_um"]["paternal"] is None for e in splits))
    check("non-split zygotes have both maternal and paternal distances",
          all(e["dist_um"]["maternal"] is not None and e["dist_um"]["paternal"] is not None
              for e in E if not e["split"]))
    check("female index matches the assignments consensus",
          all(e["female"] == (assign[e["id"]].get("consensus") or {}).get("female")
              for e in E if not e["split"]))

    print("\n[distances reproduce the source]")
    scale = d["unit_um_per_plot"]
    bad = 0
    for e in E:
        sp = assign[e["id"]]["sperm_plot"]
        pol = (assign[e["id"]].get("polar") or {}).get("com_plot")
        if pol:
            exp = round(math.dist(sp, pol) * scale, 2)
            if abs(exp - e["dist_um"]["polar"]) > 0.05:
                bad += 1
        if not e["split"]:
            pat = e["pron"][1 - e["female"]]["com"]
            exp = round(math.dist(sp, pat) * scale, 2)
            if abs(exp - e["dist_um"]["paternal"]) > 0.05:
                bad += 1
    check("recomputed distances match the stored ones", bad == 0, f"{bad} mismatch")

    print("\n[geometry & trend sanity]")
    # sperm should be closest to the paternal pronucleus (it forms from the sperm)
    closer = sum(1 for e in E if not e["split"]
                 and e["dist_um"]["paternal"] < e["dist_um"]["maternal"])
    ns = sum(1 for e in E if not e["split"])
    check("sperm is closer to paternal than maternal in the large majority",
          closer >= 0.8 * ns, f"{closer}/{ns}")
    # sperm→paternal grows with tau (paternal PN migrates inward from the cortex)
    pat = [(e["tau"], e["dist_um"]["paternal"]) for e in E if not e["split"]]
    r, slope = pearson([p[0] for p in pat], [p[1] for p in pat])
    check("sperm→paternal distance rises with τ (r > 0.3)", r > 0.3, f"r={r:.2f}")
    check("positive slope of tens of µm per τ", slope > 5, f"{slope:.0f} µm/τ")

    print("\n[provenance]")
    txt = open(AGG, "rb").read()
    leaks = [s for s in (b"/Users/", b"/Volumes/", b"C:\\") if s in txt]
    check("no absolute paths in the artifact", not leaks, str(leaks))
    check("model version recorded", d.get("model_version", "").startswith("pnpt"))

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"  FAILED: {f}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
