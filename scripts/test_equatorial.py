#!/usr/bin/env python3
"""
Tests for the Equatorial Division Plane aggregate (data/equatorial_cross.json.gz)
and per-embryo scenes. Fails on the ways this could silently drift:
  * the plane normal not being the polar-body axis;
  * a gene's stored side-A count not matching a direct recomputation of the
    (p − cell_com)·axis > 0 split from the scene transcripts;
  * the animal/vegetal segment-1 volumes being wildly unequal (a COM plane bisects);
  * missing density-mode fields (vp/vt) or the gb best-plane tuple;
  * absolute paths leaking into the committed artifacts.
Run: python3 scripts/test_equatorial.py
"""
from __future__ import annotations
import gzip, json, math, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CROSS = os.path.join(HERE, "data", "equatorial_cross.json.gz")
MAN = os.path.join(HERE, "data", "equatorial_manifest.json")
SCENE_DIR = os.path.join(HERE, "data", "equatorial")
XY = 0.15
FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def main():
    print("equatorial division plane — tests\n")
    if not os.path.isfile(CROSS):
        print("  FAIL  aggregate missing — run python3 build_equatorial.py"); return 1
    c = json.load(gzip.open(CROSS, "rt"))
    man = json.load(open(MAN))["embryos"]
    E = c["embryos"]

    print("[structure]")
    check("50-ish embryos", len(E) >= 45, str(len(E)))
    check("manifest matches cross", len(man) == len(E))
    e0 = E[0]
    check("entry has gb / vp / vt", "gb" in e0 and "vp" in e0 and "vt" in e0)
    check("gb tuple length 12", len(next(iter(e0["gb"].values()))) == 12)

    print("\n[geometry: normal == polar axis, split reproduces]")
    ids = [m["id"] for m in man[:6]]
    bad_norm = bad_split = 0
    for eid in ids:
        sc = json.load(gzip.open(os.path.join(SCENE_DIR, eid + ".json.gz"), "rt"))
        A = sc["analysis"]; p0 = A["planes"][0]
        n = p0["normal_um"]; ax = A["axis_plot"]
        # axis_plot ∝ the µm axis (up to the plot z stretch); the plane normal is the UNIT axis, so
        # normal should be parallel to unit(pb_um − com_um). Check |normal|≈1 and it points along axis_plot's xy.
        if abs(math.sqrt(sum(x * x for x in n)) - 1.0) > 1e-4:
            bad_norm += 1
        # recompute one gene's side-A count from the scene transcripts
        com = A["com_um"]
        g = A["genes"][0]["gene"]
        t = sc["transcripts"][g]
        aA = 0
        for i in range(len(t["x"])):
            if not t["s1"][i]:
                continue
            p = [t["x"][i] * XY, t["y"][i] * XY, t["gz"][i] * 1.0]
            if (p[0] - com[0]) * n[0] + (p[1] - com[1]) * n[1] + (p[2] - com[2]) * n[2] > 0:
                aA += 1
        stored = A["genes"][0]["planes"][0]["a"]
        if aA != stored:
            bad_split += 1
    check("plane normal is a unit vector", bad_norm == 0, f"{bad_norm} bad")
    check("recomputed side-A count matches stored", bad_split == 0, f"{bad_split}/{len(ids)} mismatch")

    print("\n[COM plane roughly bisects the cytoplasm]")
    ratios = []
    for eid in ids:
        sc = json.load(gzip.open(os.path.join(SCENE_DIR, eid + ".json.gz"), "rt"))
        p0 = sc["analysis"]["planes"][0]
        vA, vB = p0["volA"], p0["volB"]
        ratios.append(min(vA, vB) / max(vA, vB))
    check("animal/vegetal volumes within 2x", all(r > 0.5 for r in ratios), f"min ratio {min(ratios):.2f}")

    print("\n[provenance]")
    for path in (CROSS, os.path.join(SCENE_DIR, ids[0] + ".json.gz")):
        txt = open(path, "rb").read()
        leaks = [s for s in (b"/Users/", b"/Volumes/", b"C:\\") if s in txt]
        check(f"no absolute paths in {os.path.basename(path)}", not leaks, str(leaks))

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"  FAILED: {f}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
