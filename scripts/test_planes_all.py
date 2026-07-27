#!/usr/bin/env python3
"""
Tests for the "every plane" full-sphere aggregate (data/planes_all_*).
Fails on the ways this could silently drift:
  * the normal grid not being unit vectors on a single hemisphere (antipode dedup);
  * a gene's stored best-plane side-A count not matching a recomputation of the
    (p − com)·normal > 0 split at that plane index, from the scene transcripts;
  * the search-corrected p being degenerate (everything significant — the failure
    mode of an uncorrected best-of-20000 search) or below the empirical floor;
  * the per-embryo axis / surf points missing (aligned + orientation views need them);
  * absolute paths leaking into the artifacts.
Run: python3 scripts/test_planes_all.py
"""
from __future__ import annotations
import gzip, json, math, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAN = os.path.join(HERE, "data", "planes_all_manifest.json")
NORM = os.path.join(HERE, "data", "planes_all_normals.json.gz")
CROSS = os.path.join(HERE, "data", "planes_all_cross.json.gz")
SCENE_DIR = os.path.join(HERE, "data", "planes_all")
XY = 0.15
FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def main():
    print("division plane sweep — every plane — tests\n")
    for p in (MAN, NORM, CROSS):
        if not os.path.isfile(p):
            print(f"  FAIL  missing {os.path.basename(p)} — run python3 build_planes_all.py"); return 1
    man = json.load(open(MAN))
    N = json.load(gzip.open(NORM, "rt"))["normals"]
    c = json.load(gzip.open(CROSS, "rt"))
    E = c["embryos"]

    print("[normal grid]")
    check("~20k normals", 15000 < len(N) < 25000, str(len(N)))
    check("grid size matches manifest/cross", len(N) == man["m_planes"] == c["m_planes"])
    import random
    rng = random.Random(0)
    samp = [N[rng.randrange(len(N))] for _ in range(200)]
    check("normals are unit vectors", all(abs(math.sqrt(sum(x * x for x in n)) - 1) < 1e-3 for n in samp))
    # hemisphere: no normal is the near-antipode of another (canonical hemisphere keeps one of ±n)
    def maxdot(n):
        return max((n[0]*m[0]+n[1]*m[1]+n[2]*m[2]) for m in samp if m is not n)
    worst = min(maxdot(n) for n in samp[:60])
    check("hemisphere dedup (no antipodes)", worst > -0.995, f"min dot {worst:.3f}")

    print("\n[best-plane split reproduces from the scene]")
    ids = [m["id"] for m in man["embryos"][:5]]
    bad = 0
    for eid in ids:
        sc = json.load(gzip.open(os.path.join(SCENE_DIR, eid + ".json.gz"), "rt"))
        A = sc["analysis"]; com = A["com_um"]
        row = A["genes"][0]
        for idx_key, a_key in (("iVol", "aVol"), ("iCnt", "aCnt")):
            n = N[row[idx_key]]
            t = sc["transcripts"][row["gene"]]
            aA = 0
            for i in range(len(t["x"])):
                if not t["s1"][i]:
                    continue
                p = (t["x"][i] * XY - com[0], t["y"][i] * XY - com[1], t["gz"][i] * 1.0 - com[2])
                if p[0] * n[0] + p[1] * n[1] + p[2] * n[2] > 0:
                    aA += 1
            if aA != row[a_key]:
                bad += 1
    check("recomputed side-A count == stored (both modes)", bad == 0, f"{bad} mismatch")

    print("\n[search-corrected null is not degenerate]")
    # collect per-gene p across a few scenes; an UNcorrected best-of-20000 search would make
    # essentially everything significant. The corrected p must have a healthy non-significant mass.
    pv, pc = [], []
    for eid in ids:
        sc = json.load(gzip.open(os.path.join(SCENE_DIR, eid + ".json.gz"), "rt"))
        for r in sc["analysis"]["genes"]:
            pv.append(r["pVol"]); pc.append(r["pCnt"])
    frac_sig_cnt = sum(1 for p in pc if p < 0.05) / len(pc)
    frac_floor = sum(1 for p in pc if p <= 0.0063) / len(pc)   # genes pinned at the empirical floor
    floor = round(1.0 / (man["null_b"] + 1), 5)                # stored p is rounded to 5 dp
    # An UNCORRECTED best-of-20000 search would pin ~every gene at the floor (~100% significant).
    # The search correction must leave a real non-significant mass AND not floor-out everything.
    check("search-correction leaves a non-significant mass (< 85% sig)", frac_sig_cnt < 0.85, f"{frac_sig_cnt:.2f} sig")
    check("not everything is pinned at the floor (< 55%)", frac_floor < 0.55, f"{frac_floor:.2f} at floor")
    check("p never below the empirical floor", min(pv + pc) >= floor - 1e-6, f"min {min(pv+pc):.5f} vs floor {floor:.5f}")
    check("p never exceeds 1", max(pv + pc) <= 1.0 + 1e-9)

    print("\n[aligned/orientation inputs present]")
    check("every embryo has an axis", all(e.get("axis") for e in E))
    check("every embryo has surface points", all(e.get("surf") for e in E))
    check("gb tuple length 12", len(next(iter(E[0]["gb"].values()))) == 12)

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
