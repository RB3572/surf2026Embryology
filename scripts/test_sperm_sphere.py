#!/usr/bin/env python3
"""
Tests for the Sperm-Entry-Site Enrichment aggregate (data/sperm_sphere.json.gz).
Fails on the ways this could silently drift:
  * the fold not equalling (nsph/Vsph)/(n/Vcell) for the stored counts/volumes;
  * n_sph or the sphere volume not increasing monotonically with radius;
  * p_null not equal to Vsph/Vcell, or not increasing with radius;
  * empirical spatial p's out of (0,1];
  * labels not in the canonical Z-P?-fov? form;
  * absolute paths leaking into the artifact.
Run: python3 scripts/test_sperm_sphere.py
"""
from __future__ import annotations
import gzip, json, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AGG = os.path.join(HERE, "data", "sperm_sphere.json.gz")
FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def main():
    print("sperm-entry-site enrichment — tests\n")
    if not os.path.isfile(AGG):
        print("  FAIL  aggregate missing — run python3 build_sperm_sphere.py"); return 1
    d = json.load(gzip.open(AGG, "rt"))
    R = d["radii"]; E = d["embryos"]

    print("[structure]")
    check("radii are increasing", all(R[i] < R[i + 1] for i in range(len(R) - 1)), str(R))
    check("has sperm zygotes", len(E) >= 20, str(len(E)))
    check("meta defaultRadiusIdx valid", 0 <= d["meta"]["defaultRadiusIdx"] < len(R))
    check("byGene present", isinstance(d.get("byGene"), dict) and len(d["byGene"]) > 0)

    print("\n[per-zygote sphere geometry]")
    bad_pnull = bad_mono = bad_vmono = 0
    for e in E:
        Vcell = e["V_cell"]
        vs = [e["v_sph"][str(r)] for r in R]
        pn = [e["p_null"][str(r)] for r in R]
        if any(abs(pn[i] - vs[i] / Vcell) > 1e-4 for i in range(len(R))):
            bad_pnull += 1
        if any(vs[i] > vs[i + 1] for i in range(len(R) - 1)):
            bad_vmono += 1
        # every gene: n_sph non-decreasing with radius (bigger sphere ⊇ smaller)
        for rec in e["genes"].values():
            if any(rec["nsph"][i] > rec["nsph"][i + 1] for i in range(len(R) - 1)):
                bad_mono += 1; break
    check("p_null == Vsph/Vcell", bad_pnull == 0, f"{bad_pnull} bad")
    check("sphere volume grows with radius", bad_vmono == 0, f"{bad_vmono} bad")
    check("n_sph grows with radius (nested spheres)", bad_mono == 0, f"{bad_mono} zygotes bad")

    print("\n[fold + p sanity]")
    bad_fold = bad_p = 0
    checked = 0
    for e in E:
        Vcell = e["V_cell"]
        for rec in e["genes"].values():
            n = rec["n"]
            for i, r in enumerate(R):
                vsph = e["v_sph"][str(r)]
                exp_fold = round((rec["nsph"][i] / vsph) / (n / Vcell), 3) if (vsph and n) else 0.0
                if abs(exp_fold - rec["fold"][i]) > 0.02:
                    bad_fold += 1
                for pk in ("pE", "pD", "pSE", "pSD"):
                    v = rec[pk][i]
                    if not (0 <= v <= 1.0001):
                        bad_p += 1
                checked += 1
    check("fold == (nsph/Vsph)/(n/Vcell)", bad_fold == 0, f"{bad_fold}/{checked}")
    check("all p-values in (0,1]", bad_p == 0, f"{bad_p} bad")

    print("\n[labels + provenance]")
    import re
    ok_lbl = sum(1 for e in E if re.match(r"^Z-P\d+-fov", e["label"]))
    check("labels are canonical Z-P?-fov?", ok_lbl == len(E), f"{ok_lbl}/{len(E)}")
    txt = open(AGG, "rb").read()
    leaks = [s for s in (b"/Users/", b"/Volumes/", b"C:\\") if s in txt]
    check("no absolute paths in the artifact", not leaks, str(leaks))

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"  FAILED: {f}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
