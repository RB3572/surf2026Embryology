#!/usr/bin/env python3
"""
Enrich the cross-embryo aggregates with per-side SEGMENT-1 VOLUMES at all 18 planes.

The concordance grids + per-side bars gained a DENSITY mode (count ÷ side volume, so a plane that
splits a zygote into very uneven halves is still comparable). Density needs the side-A volume at
each plane; the aggregate only stored counts. This adds, per embryo (gene-independent):

    vp = [volA@plane0, …, volA@plane17]     (segment-1 volume, µm³, on the positive-normal side A)
    vt = volA + volB                         (total segment-1 volume — constant across planes)

so side-B volume at plane k = vt − vp[k]. Read straight from the built scenes
(analysis.planes[k].volA/volB) — no re-analysis. Run from the deploy repo root:
    python3 scripts/add_volumes_to_cross_agg.py
"""
import gzip
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCENES = os.path.join(ROOT, "data", "zygote")
TARGETS = [
    (os.path.join(ROOT, "data", "zygote_cross.json.gz"), "real"),
    (os.path.join(ROOT, "data", "zygote_cross_circ.json.gz"), "circ"),
]


def scene_volumes(scene, which):
    analysis = scene.get("analysis") if which == "real" else (scene.get("circ") or {}).get("analysis")
    planes = (analysis or {}).get("planes")
    if not planes:
        return None, None
    vp = [round(float(p.get("volA", 0)), 1) for p in planes]
    vt = round(float(planes[0].get("volA", 0)) + float(planes[0].get("volB", 0)), 1)
    return vp, vt


def main():
    for path, which in TARGETS:
        if not os.path.exists(path):
            print(f"  -- {os.path.basename(path)} missing, skipped"); continue
        agg = json.load(gzip.open(path, "rt"))
        before = os.path.getsize(path); n_ok = 0
        for emb in agg.get("embryos", []):
            sp = os.path.join(SCENES, emb["id"] + ".json.gz")
            if not os.path.exists(sp):
                print(f"  !! no scene for {emb['id']}"); continue
            vp, vt = scene_volumes(json.load(gzip.open(sp, "rt")), which)
            if vp is None:
                print(f"  !! no {which} analysis for {emb['id']}"); continue
            emb["vp"] = vp; emb["vt"] = vt; n_ok += 1
        with gzip.open(path, "wt") as fh:
            json.dump(agg, fh, separators=(",", ":"))
        print(f"  {os.path.basename(path)} ({which}): +vp/vt on {n_ok} embryos · "
              f"{before/1024:.0f}→{os.path.getsize(path)/1024:.0f} KB")


if __name__ == "__main__":
    main()
