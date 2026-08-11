#!/usr/bin/env python3
"""Checks on data/size.json.gz — embryo size.

The page draws a sphere and a box on top of a real embryo, so the numbers must be physically
plausible AND geometrically consistent with the shape that gets drawn. These guard both, plus the
two invariants the project exists to preserve: every embryo is included, and the definitions match
figure 10.1 so the site and the figure cannot drift apart.
"""
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ART = os.path.join(ROOT, "data", "size.json.gz")

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
        sys.exit("size.json.gz missing — run: python3 build_size.py")
    d = json.load(gzip.open(ART, "rt"))
    m, emb = d["meta"], d["embryos"]
    zyg = [e for e in emb if e["stage"] == "zygote"]
    two = [e for e in emb if e["stage"] != "zygote"]

    section("shape")
    check("has embryos", len(emb) > 100, str(len(emb)))
    check("meta count matches", m["n"] == len(emb))
    check("all three stages present", m["n_zygote"] and m["n_e2c"] and m["n_l2c"],
          f"{m['n_zygote']}/{m['n_e2c']}/{m['n_l2c']}")
    check("stage counts add up", m["n_zygote"] + m["n_e2c"] + m["n_l2c"] == len(emb))
    check("every embryo names a scene file", all(e.get("scene") for e in emb))
    check("every stage label is known",
          all(e["stage"] in ("zygote", "early2cell", "late2cell") for e in emb))

    section("the mean radius exists for EVERY embryo")
    # This is the project's central promise: one comparable size number at every stage, with
    # nothing dropped for being an awkward shape.
    check("every embryo has a radius", all("radius_um" in e for e in emb))
    check("every embryo has min/max and a CV",
          all({"radius_min_um", "radius_max_um", "radial_cv_pct"} <= set(e) for e in emb))
    check("2-cell embryos are included in the radius, not just zygotes",
          all("radius_um" in e for e in two), f"{sum('radius_um' in e for e in two)}/{len(two)}")
    r = [e["radius_um"] for e in emb]
    check("radii are physical (20-70 µm)", all(20 <= v <= 70 for v in r),
          f"{min(r):.1f}–{max(r):.1f}")
    check("min <= mean <= max for every embryo",
          all(e["radius_min_um"] <= e["radius_um"] <= e["radius_max_um"] + 1e-6 for e in emb))
    check("enough directions sampled to be a mean, not a spot check",
          all(e["n_dirs"] >= 50 for e in emb), str(min(e["n_dirs"] for e in emb)))
    check("nothing is excluded: noted embryos are still present",
          all(any(e["id"] == i for e in emb) for i in m.get("noted", [])))

    section("the zygote sphere is drawable")
    check("every zygote has a centre", all(len(e.get("centre_um", [])) == 3 for e in zyg))
    check("centres are finite", all(all(np.isfinite(e["centre_um"])) for e in zyg))

    section("the 2-cell box is drawable and actually contains the pair")
    check("every 2-cell embryo has a box", all("box" in e for e in two))
    bad_basis = []
    for e in two:
        b = e["box"]
        u, v, w = (np.array(b[k]) for k in ("u", "v", "w"))
        if (abs(np.linalg.norm(u) - 1) > 1e-3 or abs(np.linalg.norm(v) - 1) > 1e-3
                or abs(np.linalg.norm(w) - 1) > 1e-3
                or abs(u @ v) > 1e-3 or abs(u @ w) > 1e-3 or abs(v @ w) > 1e-3):
            bad_basis.append(e["id"])
    check("(u, v, w) is an orthonormal frame", not bad_basis, str(bad_basis[:3]))
    bad_axis = []
    for e in two:
        a, b2 = np.array(e["com_a_um"]), np.array(e["com_b_um"])
        u = np.array(e["box"]["u"])
        if np.linalg.norm(u - (b2 - a) / np.linalg.norm(b2 - a)) > 1e-3:
            bad_axis.append(e["id"])
    check("u really is the blastomere axis", not bad_axis, str(bad_axis[:3]))
    check("length equals twice the box's along-axis half-extent",
          all(abs(e["length_um"] - 2 * e["box"]["half_u"]) < 0.01 for e in two))
    check("both blastomere centres lie inside the box",
          all(abs((np.array(c) - np.array(e["box"]["centre_um"])) @ np.array(e["box"]["u"]))
              <= e["box"]["half_u"] + 1e-6
              for e in two for c in (e["com_a_um"], e["com_b_um"])))
    check("box dimensions are physical (40-200 µm)",
          all(40 <= e["length_um"] <= 200 and 40 <= e["height_um"] <= 200 for e in two),
          f"L {min(e['length_um'] for e in two):.0f}–{max(e['length_um'] for e in two):.0f}")
    check("length exceeds height (the pair is longer than it is tall)",
          all(e["length_um"] > e["height_um"] for e in two),
          str([e["id"] for e in two if e["length_um"] <= e["height_um"]][:3]))
    check("separation is smaller than the length it sits inside",
          all(e["sep_um"] < e["length_um"] for e in two))

    # The blastomeres press together and flatten, so the pair end to end is markedly shorter than
    # two free spheres. If this ever failed, either the geometry or that claim would be wrong.
    ratio = [e["length_um"] / e["height_um"] for e in two]
    check("length is well short of twice the height (they overlap, not abut)",
          np.median(ratio) < 1.8, f"median length/height = {np.median(ratio):.2f}")

    section("consistency with figure 10.1")
    med_z = float(np.median([e["radius_um"] for e in zyg]))
    med_e_len = float(np.median([e["length_um"] for e in emb if e["stage"] == "early2cell"]))
    med_l_len = float(np.median([e["length_um"] for e in emb if e["stage"] == "late2cell"]))
    check("zygote median radius ≈ 40.4 µm", abs(med_z - 40.4) < 0.5, f"{med_z:.2f}")
    check("early 2-cell median length ≈ 114.5 µm", abs(med_e_len - 114.5) < 1.0, f"{med_e_len:.2f}")
    check("late 2-cell median length ≈ 110.3 µm", abs(med_l_len - 110.3) < 1.0, f"{med_l_len:.2f}")
    check("a zygote is about 80 µm across, as a mouse zygote should be",
          70 <= 2 * med_z <= 90, f"{2 * med_z:.1f} µm diameter")
    check("meta medians agree with the embryo records",
          abs(m["medians"]["zygote_radius_um"] - med_z) < 0.01)

    section("the ranking the right drawer offers")
    for grp in ("zygote", "twocell"):
        keys = [o["key"] for o in m["sorts"][grp]]
        pool = zyg if grp == "zygote" else two
        missing = [k for k in keys if not all(k in e for e in pool)]
        check(f"every {grp} sort key exists on every {grp} embryo", not missing, str(missing))
    check("radius is offered as a sort at both stages",
          all(any(o["key"] == "radius_um" for o in m["sorts"][g]) for g in ("zygote", "twocell")))

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("size-"))
    check("no absolute paths in the artifact",
          "/Users/" not in json.dumps(d) and "/Volumes/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
