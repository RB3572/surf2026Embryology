#!/usr/bin/env python3
"""
Tests for the Compare Division Planes aggregate (data/compare_planes.json.gz).
Guards the invariants the front-end relies on:
  * the four plane pathways are present and named;
  * the polar-body axis is a unit vector; fixed-plane volume splits sum to Vtot;
  * has_sperm ⟺ the sperm plane exists;
  * every stored gene clears the min-count floor and its per-plane side counts
    sum to the in-cell total and never exceed it; per-side volumes are positive
    and sum to Vtot;
  * no absolute paths leak into the artifact.
Run: python3 scripts/test_compare_planes.py
"""
from __future__ import annotations
import gzip, json, math, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AGG = os.path.join(HERE, "data", "compare_planes.json.gz")
FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def approx(a, b, tol=2.0):
    return abs(a - b) <= tol


def main():
    print("compare division planes — tests\n")
    if not os.path.isfile(AGG):
        print("  FAIL  aggregate missing — run python3 build_compare_planes.py"); return 1
    d = json.load(gzip.open(AGG, "rt"))
    E = d["embryos"]; minc = d["minCount"]

    print("[structure]")
    check("four plane pathways", [p["key"] for p in d["planes"]] == ["polar", "exhaustive", "equatorial", "sperm"],
          str([p.get("key") for p in d["planes"]]))
    check("enough zygotes", len(E) >= 20, str(len(E)))
    check("declared n matches", d["n"] == len(E))
    check("unit µm/plot is 0.15", abs(d["unit_um_per_plot"] - 0.15) < 1e-9)
    check("declared n_sperm matches", d["n_sperm"] == sum(1 for e in E if e["has_sperm"]))

    print("\n[per-zygote geometry]")
    bad_axis = bad_eqv = bad_sd = bad_sperm = 0
    for e in E:
        ax = e["axis_um"]
        if abs(math.sqrt(sum(x * x for x in ax)) - 1.0) > 1e-3:
            bad_axis += 1
        if not approx(e["eq"]["v"][0] + e["eq"]["v"][1], e["Vtot"], max(2.0, 0.001 * e["Vtot"])):
            bad_eqv += 1
        if (e["sd"] is not None) != e["has_sperm"]:
            bad_sperm += 1
        if e["sd"] is not None and not approx(e["sd"]["v"][0] + e["sd"]["v"][1], e["Vtot"], max(2.0, 0.001 * e["Vtot"])):
            bad_sd += 1
    check("polar axis is a unit vector", bad_axis == 0, f"{bad_axis} bad")
    check("equatorial vA+vB = Vtot", bad_eqv == 0, f"{bad_eqv} bad")
    check("has_sperm ⟺ sperm plane present", bad_sperm == 0, f"{bad_sperm} bad")
    check("sperm vA+vB = Vtot", bad_sd == 0, f"{bad_sd} bad")

    print("\n[per-gene splits]")
    bad_minc = bad_eqc = bad_pb = bad_ex = bad_vol = ngene = 0
    for e in E:
        Vtot = e["Vtot"]
        for g, r in e["genes"].items():
            ngene += 1
            nc = r["nc"]
            if nc < minc:
                bad_minc += 1
            if not (0 <= r["eq"] <= nc):
                bad_eqc += 1
            if r["sd"] is not None and not (0 <= r["sd"] <= nc):
                bad_eqc += 1
            pb = r["pb"]
            if not (pb["c"][0] + pb["c"][1] == nc and 0 <= pb["c"][0] <= nc and 0 <= pb["ang"] < 180):
                bad_pb += 1
            if not (pb["v"][0] > 0 and pb["v"][1] > 0 and approx(pb["v"][0] + pb["v"][1], Vtot, max(2.0, 0.001 * Vtot))):
                bad_vol += 1
            ex = r["ex"]
            if ex is not None:
                if not (ex["c"][0] + ex["c"][1] == nc and ex["v"][0] > 0 and ex["v"][1] > 0):
                    bad_ex += 1
                if not approx(ex["v"][0] + ex["v"][1], Vtot, max(2.0, 0.001 * Vtot)):
                    bad_vol += 1
    check("every gene clears the min-count floor", bad_minc == 0, f"{bad_minc}/{ngene}")
    check("equatorial/sperm side counts in [0,nc]", bad_eqc == 0, f"{bad_eqc} bad")
    check("polar-best counts sum to nc + angle in range", bad_pb == 0, f"{bad_pb} bad")
    check("exhaustive-best counts sum to nc", bad_ex == 0, f"{bad_ex} bad")
    check("per-side volumes positive and sum to Vtot", bad_vol == 0, f"{bad_vol} bad")

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
