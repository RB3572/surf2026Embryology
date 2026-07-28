#!/usr/bin/env python3
"""
Tests for the Sperm-Entry-Site Enrichment aggregate (data/sperm_sphere.json.gz).

The build now ships RAW per-SEGMENT counts and voxel volumes; the fold, binomial p
and 95% band are all computed in the browser for whatever segment set is selected.
So the tests guard the invariants that keep that client math valid:
  * radii strictly increasing; enough zygotes; segMeta matches segs;
  * per segment, the sphere voxel count grows with radius and never exceeds the cell;
  * per (gene, segment), in-sphere counts grow with radius and never exceed in-cell;
  * every stored gene clears the min-count floor (Σ nc ≥ minCount);
  * a gene only carries counts in segments the zygote actually has (nc>0 ⇒ vc>0);
  * re-deriving the front-end statistic (all segments on, default radius) gives a
    finite fold, a volume fraction p0 in (0,1], and an occupancy in [0,1];
  * labels are canonical Z-P?-fov?; no absolute paths leak into the artifact.
Run: python3 scripts/test_sperm_sphere.py
"""
from __future__ import annotations
import gzip, json, os, re, sys

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
    R = d["radii"]; E = d["embryos"]; SEGS = d["segs"]; nseg = len(SEGS); nr = len(R)
    minc = d["meta"]["minCount"]

    print("[structure]")
    check("radii are increasing", all(R[i] < R[i + 1] for i in range(nr - 1)), str(R))
    check("has sperm zygotes", len(E) >= 20, str(len(E)))
    check("meta.nZygotes matches", d["meta"].get("nZygotes") == len(E), f"{d['meta'].get('nZygotes')} vs {len(E)}")
    check("meta defaultRadiusIdx valid", 0 <= d["meta"]["defaultRadiusIdx"] < nr)
    check("segMeta matches segs", [m["key"] for m in d["segMeta"]] == SEGS, str(SEGS))

    print("\n[per-zygote sphere volumes]")
    bad_shape = bad_vmono = bad_vsubset = bad_present = 0
    for e in E:
        vc, vs = e["vc"], e["vs"]
        if len(vc) != nseg or len(vs) != nseg or any(len(row) != nr for row in vs):
            bad_shape += 1; continue
        for s in range(nseg):
            if any(vs[s][i] > vs[s][i + 1] for i in range(nr - 1)):
                bad_vmono += 1
            if any(vs[s][i] > vc[s] for i in range(nr)):
                bad_vsubset += 1
        if [SEGS[s] for s in range(nseg) if vc[s] > 0] != list(e["present"]):
            bad_present += 1
    check("vc/vs shapes are [nseg]/[nseg][nr]", bad_shape == 0, f"{bad_shape} bad")
    check("sphere voxels grow with radius", bad_vmono == 0, f"{bad_vmono} bad")
    check("sphere voxels never exceed the cell", bad_vsubset == 0, f"{bad_vsubset} bad")
    check("present == segments with vc>0", bad_present == 0, f"{bad_present} bad")

    print("\n[per-gene counts]")
    bad_gshape = bad_nmono = bad_nsubset = bad_minc = bad_seg = 0
    ngene = 0
    for e in E:
        vc = e["vc"]
        for rec in e["genes"].values():
            ngene += 1
            nc, ns = rec["nc"], rec["ns"]
            if len(nc) != nseg or len(ns) != nseg or any(len(row) != nr for row in ns):
                bad_gshape += 1; continue
            if sum(nc) < minc:
                bad_minc += 1
            for s in range(nseg):
                if any(ns[s][i] > ns[s][i + 1] for i in range(nr - 1)):
                    bad_nmono += 1
                if any(ns[s][i] > nc[s] for i in range(nr)):
                    bad_nsubset += 1
                if nc[s] > 0 and vc[s] == 0:          # a gene can't sit in a segment the zygote lacks
                    bad_seg += 1
    check("nc/ns shapes are [nseg]/[nseg][nr]", bad_gshape == 0, f"{bad_gshape} bad")
    check("in-sphere counts grow with radius", bad_nmono == 0, f"{bad_nmono} bad")
    check("in-sphere counts never exceed in-cell", bad_nsubset == 0, f"{bad_nsubset} bad")
    check("every gene clears the min-count floor", bad_minc == 0, f"{bad_minc}/{ngene}")
    check("gene counts only in present segments", bad_seg == 0, f"{bad_seg} bad")

    print("\n[front-end statistic — all segments on, default radius]")
    ri = d["meta"]["defaultRadiusIdx"]
    bad_p0 = bad_fold = bad_occ = 0
    for e in E:
        vc, vs = e["vc"], e["vs"]
        Vc = sum(vc); Vs = sum(vs[s][ri] for s in range(nseg))
        p0 = Vs / Vc if Vc else 0.0
        if not (0 < p0 <= 1.0001):
            bad_p0 += 1
        tot_ns = tot_nc = 0
        for rec in e["genes"].values():
            nS = sum(rec["ns"][s][ri] for s in range(nseg))
            nC = sum(rec["nc"][s] for s in range(nseg))
            tot_ns += nS; tot_nc += nC
            mu = nC * p0
            fold = (nS / mu) if mu > 0 else 0.0     # = (nS/Vs)/(nC/Vc), the client's fold
            if not (fold >= 0 and fold == fold and fold != float("inf")):
                bad_fold += 1
        occ = (tot_ns / tot_nc) if tot_nc else 0.0
        if not (0 <= occ <= 1.0001):
            bad_occ += 1
    check("volume fraction p0 in (0,1]", bad_p0 == 0, f"{bad_p0} bad")
    check("fold is finite and non-negative", bad_fold == 0, f"{bad_fold} bad")
    check("occupancy fraction in [0,1]", bad_occ == 0, f"{bad_occ} bad")

    print("\n[labels + provenance]")
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
