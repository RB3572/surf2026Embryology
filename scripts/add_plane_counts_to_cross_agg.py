#!/usr/bin/env python3
"""
Enrich the cross-embryo aggregates with per-gene side-A counts at ALL 18 planes.

The γ/μ concordance grid needs a negative control: pick a RANDOM plane and randomly
call one side γ for each zygote. The aggregate only stored side-A counts at the four
BEST planes, which is not enough — so this adds, per embryo:

    gp[gene] = [a@plane0, a@plane1, ..., a@plane17]     (side-A transcript counts)

`total` stays in the existing g[gene][0], and the existing g[...] best-plane entries
are left untouched so nothing downstream changes. Values are read straight out of the
already-built per-embryo scenes (analysis.genes[].planes[].a) — no re-analysis.

Run from the deploy repo root:  python3 scripts/add_plane_counts_to_cross_agg.py
"""
import glob
import gzip
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCENES = os.path.join(ROOT, "data", "zygote")

# (aggregate file, which analysis inside the scene to read)
TARGETS = [
    (os.path.join(ROOT, "data", "zygote_cross.json.gz"), "real"),
    (os.path.join(ROOT, "data", "zygote_cross_circ.json.gz"), "circ"),
]


def scene_plane_counts(scene, which):
    """{gene: [side-A count at each of the 18 planes]} for the real or circ analysis."""
    analysis = scene.get("analysis") if which == "real" else (scene.get("circ") or {}).get("analysis")
    if not analysis:
        return None
    out = {}
    for row in analysis.get("genes", []):
        planes = row.get("planes") or []
        out[row["gene"]] = [int(p.get("a", 0)) for p in planes]
    return out


def main():
    for path, which in TARGETS:
        if not os.path.exists(path):
            print(f"  -- {os.path.basename(path)} missing, skipped")
            continue
        agg = json.load(gzip.open(path, "rt"))
        before = os.path.getsize(path)
        n_ok = n_planes = 0
        for emb in agg.get("embryos", []):
            sp = os.path.join(SCENES, emb["id"] + ".json.gz")
            if not os.path.exists(sp):
                print(f"  !! no scene for {emb['id']}")
                continue
            counts = scene_plane_counts(json.load(gzip.open(sp, "rt")), which)
            if counts is None:
                print(f"  !! no {which} analysis for {emb['id']}")
                continue
            # keep only the genes this aggregate already carries, so the two stay in sync
            gp = {g: counts[g] for g in emb.get("g", {}) if g in counts}
            if gp:
                n_planes = max(n_planes, max(len(v) for v in gp.values()))
            emb["gp"] = gp
            n_ok += 1
        agg["n_planes"] = n_planes
        with gzip.open(path, "wt") as fh:
            json.dump(agg, fh, separators=(",", ":"))
        after = os.path.getsize(path)
        print(f"  {os.path.basename(path)}: {n_ok} embryos got gp ({n_planes} planes) "
              f"· {before/1024:.0f} KB → {after/1024:.0f} KB")


if __name__ == "__main__":
    main()
