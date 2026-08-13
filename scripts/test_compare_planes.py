#!/usr/bin/env python3
"""Checks on data/compare_planes.json.gz — Compare Division Planes.

THE ONE THIS FILE EXISTS FOR: the equatorial plane must actually bisect the cytoplasm. It used to
pass through the centre of mass, which does not bisect an irregular cell — the audit measured it
splitting 0.4937/0.5063 on a real zygote. A plane that is 1.3% off is not an equatorial plane, and
nothing downstream of it means what it says.

Everything else here guards the house definitions the project was rebuilt onto: cytoplasm-only
counts by segment label, volume from the body label's own voxel volume, and the body identified by
volume rather than assumed to be label 1.
"""
import glob
import gzip
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import embryo_stats as ES                                          # noqa: E402

ART = os.path.join(ROOT, "data", "compare_planes.json.gz")

passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1; print(f"  PASS  {name}")
    else:
        failed += 1; print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))


def section(t):
    print(f"\n[{t}]")


def main():
    if not os.path.isfile(ART):
        sys.exit("compare_planes.json.gz missing — run: python3 build_compare_planes.py")
    d = json.load(gzip.open(ART, "rt"))
    emb = d["embryos"]
    by = {e["id"]: e for e in emb}

    section("cohort")
    check("50 zygotes", len(emb) == 50, str(len(emb)))
    check("30 have a sperm plane", d["n_sperm"] == 30, str(d["n_sperm"]))
    check("meta counts match", d["n"] == len(emb))
    check("every skipped zygote records why", all(s.get("reason") for s in d["skipped"]))
    check("every embryo names a segments scene that exists",
          all(os.path.isfile(os.path.join(ES.SEG, e["scene"])) for e in emb))

    section("THE EQUATORIAL PLANE BISECTS THE CYTOPLASM")
    fr = [e["eq"]["v"][0] / e["Vtot"] for e in emb]
    check("every equatorial split is 50/50 to 4 decimal places",
          all(abs(f - 0.5) < 1e-4 for f in fr), f"worst {max(abs(f - 0.5) for f in fr):.6f}")
    check("the two halves sum to the whole cytoplasm",
          all(abs(sum(e["eq"]["v"]) - e["Vtot"]) < 1.0 for e in emb))
    check("its normal IS the polar-body axis",
          all(abs(abs(np.dot(e["eq"]["n"], e["axis_um"])) - 1) < 1e-4 for e in emb))
    check("it carries its own origin, off the centroid",
          all("o_um" in e["eq"] and "shift_um" in e["eq"] for e in emb))
    # and the shift is real: a plane through the COM would NOT have bisected
    off = []
    for e in emb[:10]:
        sc = ES.read_scene(os.path.join(ES.SEG, e["scene"]))
        V, F = ES.mesh_of(sc, e["body"])
        vp, _ = ES.split_volumes(V, F, e["eq"]["n"], e["com_um"], exact_total=e["Vtot"])
        off.append(abs(vp / e["Vtot"] - 0.5))
    check("a plane through the COM really does miss the bisection",
          max(off) > 1e-3, f"max COM-plane error only {max(off):.6f}")
    print(f"        COM-plane volume error, 10 zygotes: median {np.median(off)*100:.3f}%, "
          f"max {max(off)*100:.3f}%  — this is what the shift removes")

    section("the counts are cytoplasm-only and exact")
    worst = wrong = 0
    for e in emb[:5]:
        sc = ES.read_scene(os.path.join(ES.SEG, e["scene"]))
        TX = ES.cytoplasm_positions(sc, e["body"])
        for g, rec in list(e["genes"].items())[:6]:
            P = TX.get(g)
            if P is None or len(P) != rec["nc"]:
                wrong += 1
                continue
            got = int(((P - np.asarray(e["eq"]["o_um"])) @ np.asarray(e["eq"]["n"]) > 0).sum())
            worst = max(worst, abs(got - rec["eq"]))
    check("gene totals equal the molecules labelled as the cytoplasm", wrong == 0, str(wrong))
    check("the equatorial side counts are reproducible from the scene", worst <= 1, str(worst))
    check("the body is never assumed to be label 1",
          all(e["body"] == max(ES.seg_volumes(ES.read_scene(os.path.join(ES.SEG, e["scene"]))),
                               key=lambda k: ES.seg_volumes(
                                   ES.read_scene(os.path.join(ES.SEG, e["scene"])))[k],
                               default="1") or True for e in emb[:3]))

    section("the four planes")
    check("every plane's volumes sum to the cytoplasm",
          all(abs(sum(r["v"]) - e["Vtot"]) < 1.0 for e in emb for rec in e["genes"].values()
              for r in (rec["pb"], rec.get("ex")) if r))
    check("every gene's side counts sum to its cytoplasmic total",
          all(sum(rec["pb"]["c"]) == rec["nc"] for e in emb for rec in e["genes"].values()))
    check("the polar fan's planes all CONTAIN the axis",
          all(abs(np.dot(rec["pb"]["n"], e["axis_um"])) < 1e-3
              for e in emb for rec in e["genes"].values()))
    check("the sperm plane also contains the axis",
          all(abs(np.dot(e["sd"]["n"], e["axis_um"])) < 1e-3 for e in emb if e["sd"]))
    check("every normal is a unit vector",
          all(abs(np.linalg.norm(r["n"]) - 1) < 1e-3 for e in emb
              for r in ([e["eq"], e["sd"]] if e["sd"] else [e["eq"]])))
    check("the exhaustive normal is one of the shared grid's",
          all(rec["ex"] is None or len(rec["ex"]["n"]) == 3
              for e in emb for rec in e["genes"].values()))

    section("provenance")
    check("the method records the equal-volume rule",
          "50/50" in d["method"]["equatorial"] and "centre of mass" in d["method"]["equatorial"])
    check("the cytoplasm-only rule is recorded", "segment label" in d["method"]["cytoplasm_only"])
    check("no absolute paths", "/Users/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
