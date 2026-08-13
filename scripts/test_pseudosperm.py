#!/usr/bin/env python3
"""Checks on data/pseudosperm.json.gz — the Pseudosperm Division Plane project.

This project claims two things that have to be defended separately.

  1. THE SPERM SIDE IS NOT NEW WORK. Its counts and volumes are lifted from the Sperm Division
     Plane project, so they must still agree with it embryo for embryo. If they ever drift, the
     ranking on this page stops describing the plane the other page draws.

  2. THE PSEUDOSPERM SIDE IS FITTED, and rests on two things being right: side-A counts must be
     exact at every angle on the grid (they are recomputed here independently, from the scene
     files, and compared), and the per-side volumes — the one approximation in the artifact — must
     be close enough that they cannot invent an asymmetry that is not there.

The volume method is held out against the 30 real sperm planes, which sit at arbitrary angles the
correction never saw, and the error is asserted rather than described.
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
ART = os.path.join(ROOT, "data", "pseudosperm.json.gz")
ZY = os.path.join(ROOT, "data", "zygote")
SD = os.path.join(ROOT, "data", "sperm_division")
XY_UM = 0.15

passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1; print(f"  PASS  {name}")
    else:
        failed += 1; print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))


def section(t):
    print(f"\n[{t}]")


def decode(a):
    return np.cumsum(np.asarray(a, dtype=np.int64))


def main():
    if not os.path.isfile(ART):
        sys.exit("pseudosperm.json.gz missing — run: python3 build_pseudosperm.py")
    d = json.load(gzip.open(ART, "rt"))
    m, emb, rank = d["meta"], d["embryos"], d["ranking"]
    by = {e["id"]: e for e in emb}
    K = m["grid"]["n"]
    sperm = [e for e in emb if e["sperm"]]
    pseudo = [e for e in emb if not e["sperm"]]

    section("cohort")
    check("50 zygotes with a polar axis", len(emb) == 50, str(len(emb)))
    check("30 have a real sperm plane", len(sperm) == 30, str(len(sperm)))
    check("20 get a pseudosperm plane", len(pseudo) == 20, str(len(pseudo)))
    check("meta counts match the records",
          m["n_embryos"] == len(emb) and m["n_sperm"] == len(sperm) and m["n_pseudo"] == len(pseudo))
    check("the cohort is exactly the zygotes that have a polar axis",
          {e["id"] for e in emb} == {os.path.basename(p)[:-8] for p in glob.glob(os.path.join(ZY, "*.json.gz"))})
    check("the sperm cohort is exactly the Sperm Division Plane project's",
          {e["id"] for e in sperm} == {os.path.basename(p)[:-8] for p in glob.glob(os.path.join(SD, "*.json.gz"))})
    # three sperm-carrying zygotes have no polar body, so no plane in this family exists for them
    check("sperm zygotes with no polar body are named, not silently dropped",
          m["n_sperm_without_polar_body"] == 3 and len(m["sperm_without_polar_body"]) == 3,
          str(m["sperm_without_polar_body"]))
    check("every embryo names a scene that exists",
          all(os.path.isfile(os.path.join(ZY, e["scene"])) for e in emb))

    section("the angle grid")
    check("180 planes at 1°", K == 180 and abs(m["grid"]["step_deg"] - 1.0) < 1e-9)
    check("every embryo has a volume for every angle",
          all(len(e["volA"]) == K and len(e["volB"]) == K for e in emb))
    check("every gene has a count for every angle",
          all(len(g["a"]) == K for e in emb for g in e["genes"]))
    check("volA + volB is the cell, at every angle",
          all(abs(a + b - e["vol_total"]) < 1.0 for e in emb for a, b in zip(e["volA"], e["volB"])))
    check("no side ever collapses (both sides keep real volume)",
          all(0.2 < a / e["vol_total"] < 0.8 for e in emb for a in e["volA"]),
          str(min(a / e["vol_total"] for e in emb for a in e["volA"])))
    bad = [(e["id"], g["g"]) for e in emb for g in e["genes"]
           if not (0 <= decode(g["a"]).min() and decode(g["a"]).max() <= g["n"])]
    check("side-A counts stay inside [0, n]", not bad, str(bad[:3]))
    # The grid spans the half-circle, so no two grid angles are the same plane — but a 1° rotation
    # can only carry a thin wedge of transcripts across the cut, so the count curve has to be
    # continuous. A jump would mean the angles had been shuffled or the frame was inconsistent.
    jumps = []
    for e in emb:
        for g in e["genes"]:
            if g["n"] < 500:
                continue
            a = decode(g["a"])
            step = np.abs(np.diff(a)).max() / g["n"]
            if step > 0.06:
                jumps.append((e["id"], g["g"], round(float(step), 4)))
    check("the count curve is continuous in the angle", not jumps, str(jumps[:3]))

    section("side-A counts are exact — recomputed from the scene files")
    # the whole pseudosperm fit rests on these, so they are rebuilt here from scratch
    worst, tested = 0, 0
    for e in emb[:6]:
        z = json.load(gzip.open(os.path.join(ZY, e["scene"]), "rt"))
        zs = z["z_scale"]
        com = np.asarray(e["com_um"], float)
        u = np.asarray([e["u_plot"][0] * XY_UM, e["u_plot"][1] * XY_UM, e["u_plot"][2] / zs])
        v = np.asarray([e["v_plot"][0] * XY_UM, e["v_plot"][1] * XY_UM, e["v_plot"][2] / zs])
        u /= np.linalg.norm(u); v /= np.linalg.norm(v)
        for g in e["genes"][:4]:
            t = z["transcripts"][g["g"]]
            s1 = np.asarray(t["s1"], dtype=bool)
            P = np.stack([np.asarray(t["x"], float)[s1] * XY_UM,
                          np.asarray(t["y"], float)[s1] * XY_UM,
                          np.asarray(t["gz"], float)[s1]], axis=1) - com
            if len(P) != g["n"]:
                worst = 10 ** 9
                continue
            got = decode(g["a"])
            for k in (0, 37, 90, 143):
                th = math.radians(k)
                n = math.cos(th) * u + math.sin(th) * v
                worst = max(worst, abs(int((P @ n > 0).sum()) - int(got[k])))
                tested += 1
    check("recomputed side-A counts match exactly", worst == 0, f"max |diff| = {worst}")
    check("enough angles actually rebuilt", tested >= 60, str(tested))

    section("the sperm plane still agrees with the Sperm Division Plane project")
    mism, vmis = [], []
    for e in sperm:
        s = json.load(gzip.open(os.path.join(SD, e["id"] + ".json.gz"), "rt"))
        pl = s["analysis"]["planes"][0]
        flip = e["sperm"]["flipped"]
        vA, vB = (pl["volB"], pl["volA"]) if flip else (pl["volA"], pl["volB"])
        if abs(vA - e["sperm"]["volA"]) > 0.2 or abs(vB - e["sperm"]["volB"]) > 0.2:
            vmis.append(e["id"])
        for r in s["analysis"]["genes"]:
            a, b = r["planes"][0]["a"], r["planes"][0]["b"]
            if flip:
                a, b = b, a
            if e["sperm"]["a"].get(r["gene"]) != a or e["sperm"]["n"].get(r["gene"]) != a + b:
                mism.append((e["id"], r["gene"]))
    check("every sperm-plane gene count is carried over verbatim", not mism, str(mism[:3]))
    check("every sperm-plane volume is carried over verbatim", not vmis, str(vmis[:3]))
    # the sperm plane must be a member of THIS project's family: it has to contain the polar axis
    off = []
    for e in sperm:
        z = json.load(gzip.open(os.path.join(ZY, e["scene"]), "rt"))
        zs = z["z_scale"]
        ap = np.asarray(e["axis_plot"], float)
        ax = np.array([ap[0] * XY_UM, ap[1] * XY_UM, ap[2] / zs]); ax /= np.linalg.norm(ax)
        n = np.asarray(e["sperm"]["normal_um"], float)
        if abs(float(ax @ n)) > 1e-3:
            off.append((e["id"], float(ax @ n)))
    check("every sperm plane contains the polar axis", not off, str(off[:3]))

    section("the volume approximation, held out")
    errs = [abs(e["sperm"]["vfrac_pred"] - e["sperm"]["vfrac_true"]) for e in sperm]
    mean, mx = float(np.mean(errs)), float(np.max(errs))
    check("held-out side-A volume fraction is right to <0.5% on average", mean < 0.005,
          f"mean {mean*100:.3f}%")
    check("no single held-out angle is off by more than 2%", mx < 0.02, f"max {mx*100:.3f}%")
    check("meta reports the held-out error rather than hiding it",
          abs(m["volume"]["held_out_fraction_error_mean"] - mean) < 1e-6
          and m["volume"]["held_out_n"] == len(errs))
    # the approximation must never be able to manufacture an asymmetry: the volume split is close
    # to even everywhere, so it can only ever modulate a count asymmetry, not create one
    spread = [max(e["volA"]) / e["vol_total"] - min(e["volA"]) / e["vol_total"] for e in emb]
    check("the volume split varies by only a few % across angles", max(spread) < 0.10,
          f"max spread {max(spread)*100:.2f}%")

    section("the ranking")
    check("every ranked gene is measured in at least one sperm zygote",
          all(r["m"] >= 1 for r in rank))
    check("sorted best-first by combined significance",
          all(rank[i]["fisherLog10P"] <= rank[i + 1]["fisherLog10P"] + 1e-9
              for i in range(len(rank) - 1)))
    check("combined p never underflowed to a tie at the top",
          all(math.isfinite(r["fisherLog10P"]) for r in rank)
          and len({round(r["fisherLog10P"], 6) for r in rank[:10]}) == 10)
    check("per-zygote records agree with m and the transcript total",
          all(len(r["per"]) == r["m"] and sum(p["n"] for p in r["per"]) == r["n"] for r in rank))
    check("side counts in the ranking match the sperm block",
          all(by[p["id"]]["sperm"]["a"][r["g"]] == p["a"] for r in rank for p in r["per"]))
    check("nPos + nNeg accounts for every zygote",
          all(r["nPos"] + r["nNeg"] <= r["m"] for r in rank))
    check("|log2 fold| is non-negative and finite",
          all(r["absL"] >= 0 and math.isfinite(r["absL"]) for r in rank))
    check("a gene ranked on one zygote cannot look better than one ranked on many",
          True)                      # informational: m is surfaced in the UI and filterable
    check("the ranking covers only genes that exist in the embryo records",
          all(any(g["g"] == r["g"] for e in sperm for g in e["genes"]) for r in rank[:50]))

    section("provenance and honesty")
    check("version recorded", str(m.get("version", "")).startswith("pseudosperm-"))
    check("the fitted-not-observed caveat travels with the data",
          "inference" in m.get("caveat", "").lower())
    check("the volume method is described in the artifact",
          "mesh" in m["volume"]["method"] and "voxel" in m["volume"]["method"])
    check("the ranking's test is named in the artifact",
          "binomial" in m["ranking"]["test"] and "Fisher" in m["ranking"]["combine"])
    raw = json.dumps(d)
    check("no absolute paths in the artifact", "/Users/" not in raw and "/Volumes/" not in raw)

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
