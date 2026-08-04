#!/usr/bin/env python3
"""Build data/contact.json.gz — genes enriched where the two blastomeres touch.

THE QUESTION
------------
A 2-cell embryo has an interface: the flat face where the two blastomeres press against each
other. Is anything concentrated there? Cell-cell contacts are where junctional and polarity
machinery would be expected to sit, so a gene over-represented at the interface is worth a look.

DEFINING "THE CONTACT REGION"
-----------------------------
The two blastomeres are the two largest segments (the scene ships explicit voxel volumes, so
this is read off rather than assumed). Their mesh centroids give an axis

    û = unit(COM_B − COM_A),      M = midpoint

and every transcript gets a signed coordinate t = (p − M)·û — its position along the axis that
runs from one cell, through the interface, into the other. The contact region is the slab
|t| ≤ D. D is swept so the page can offer a thickness slider rather than baking one number in.

Only transcripts inside the two blastomeres count. Nuclei and the polar body are separate
compartments; including them would let a nucleus that happens to sit near the interface
masquerade as contact enrichment.

WHY ENRICHMENT IS MEASURED AGAINST THE OTHER TRANSCRIPTS, NOT AGAINST VOLUME
---------------------------------------------------------------------------
The obvious test — transcripts per µm³ in the slab vs the rest of the cell — needs the volume of
the slab∩cell intersection, and we ship meshes rather than voxel masks, so that volume would have
to be estimated. Instead each gene is compared to **where transcripts generally are in that same
embryo**: the null is that the gene is distributed like the rest of the transcriptome. That
requires no geometry beyond the axis, and it automatically absorbs cell shape, size, orientation
and detection efficiency.

    f0 = (all transcripts with |t| ≤ D) / (all transcripts)        per embryo, per D
    k  = that gene's transcripts in the slab,  n = that gene's total
    fold = (k/n) / f0        p = two-sided binomial(k; n, f0)

fold > 1 means the gene is over-represented at the interface relative to the transcriptome.

Output: data/contact.json.gz
"""
import glob
import gzip
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SRC = os.path.join(DATA, "segments")
OUT = os.path.join(DATA, "contact.json.gz")

XY_UM = 0.15
RADII = [2, 3, 4, 5, 6, 8, 10, 12, 16, 20]     # slab half-thickness D, µm
MIN_TX = 20                                     # a gene needs this many in the embryo to be scored


def unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / n if n else v


def mesh_centroid(mesh):
    v = np.asarray(mesh["verts"], float).reshape(-1, 3)
    return v.mean(axis=0)


def stage_of(eid):
    low = eid.lower()
    if "e2c" in low or "early2cell" in low:
        return "e2c"
    if "l2c" in low or "late2cell" in low:
        return "l2c"
    return None


def process(path):
    s = json.load(gzip.open(path, "rt"))
    scene_file = os.path.basename(path)
    eid = s["id"].split("__")[-1]
    stage = stage_of(eid) or stage_of(os.path.basename(path))
    if stage is None:
        return None

    segs = s.get("segments") or []
    if len(segs) < 2:
        return None
    # the two blastomeres = the two largest segments by voxel volume
    order = sorted(segs, key=lambda x: -float(x.get("volume") or 0))
    a_lbl, b_lbl = int(order[0]["label"]), int(order[1]["label"])
    va, vb = float(order[0]["volume"]), float(order[1]["volume"])
    # a real 2-cell has two comparable blastomeres; a lopsided pair means the segmentation
    # merged or split something, and the interface would be meaningless
    if vb <= 0 or va / vb > 3.0:
        return None

    meshes = s.get("region_meshes") or {}
    if str(a_lbl) not in meshes or str(b_lbl) not in meshes:
        return None
    zs = s.get("z_scale", 6.667)
    # mesh verts are in PLOT space (x,y px; z scaled) — convert to µm to match transcripts
    def to_um(c):
        return np.array([c[0] * XY_UM, c[1] * XY_UM, c[2] / zs])
    ca, cb = to_um(mesh_centroid(meshes[str(a_lbl)])), to_um(mesh_centroid(meshes[str(b_lbl)]))
    axis = unit(cb - ca)
    if not np.any(axis):
        return None
    mid = (ca + cb) / 2.0
    sep = float(np.linalg.norm(cb - ca))

    # signed position along the interface axis, blastomere transcripts only
    per_gene = {}
    all_t = []
    for g, t in (s.get("transcripts") or {}).items():
        seg = t.get("s")
        if seg is None:
            continue
        seg = np.asarray(seg)
        m = (seg == a_lbl) | (seg == b_lbl)
        if not m.any():
            continue
        P = np.stack([np.asarray(t["x"], float)[m] * XY_UM,
                      np.asarray(t["y"], float)[m] * XY_UM,
                      np.asarray(t["gz"], float)[m]], axis=1)
        tt = (P - mid) @ axis
        per_gene[g] = tt
        all_t.append(tt)
    if not all_t:
        return None
    all_t = np.concatenate(all_t)
    n_all = len(all_t)
    if n_all < 200:
        return None

    f0 = [float((np.abs(all_t) <= D).sum()) / n_all for D in RADII]
    genes = {}
    for g, tt in per_gene.items():
        n = len(tt)
        if n < MIN_TX:
            continue
        genes[g] = {"n": int(n), "k": [int((np.abs(tt) <= D).sum()) for D in RADII]}
    if not genes:
        return None

    return {
        "id": eid, "stage": stage, "scene": scene_file,
        "a": a_lbl, "b": b_lbl, "sep_um": round(sep, 2),
        "n_all": int(n_all), "f0": [round(v, 6) for v in f0],
        "n_genes": len(genes), "genes": genes,
        # everything the page needs to draw the embryo + interface without re-deriving it
        "com_a_um": [round(float(v), 3) for v in ca],
        "com_b_um": [round(float(v), 3) for v in cb],
        "axis_um": [round(float(v), 5) for v in axis],
        "mid_um": [round(float(v), 3) for v in mid],
    }


def main():
    files = sorted(glob.glob(os.path.join(SRC, "*.json.gz")))
    files = [f for f in files if stage_of(os.path.basename(f))]
    if not files:
        sys.exit("no 2-cell scenes in data/segments — nothing to do")
    print(f"contact: scanning {len(files)} two-cell scenes\n")

    out, skipped = [], 0
    for f in files:
        try:
            r = process(f)
        except Exception as e:                    # noqa: BLE001
            print(f"  !! {os.path.basename(f)}: {e}"); skipped += 1; continue
        if r is None:
            skipped += 1; continue
        out.append(r)

    if not out:
        sys.exit("no usable embryos")
    out.sort(key=lambda r: (r["stage"], r["id"]))
    n_e2c = sum(1 for r in out if r["stage"] == "e2c")
    n_l2c = len(out) - n_e2c
    universe = sorted({g for r in out for g in r["genes"]})

    payload = {
        "meta": {
            "version": "contact-1.0.0",
            "n_embryos": len(out), "n_e2c": n_e2c, "n_l2c": n_l2c,
            "skipped": skipped, "radii": RADII, "min_tx": MIN_TX,
            "default_radius": 6, "n_genes": len(universe),
        },
        "genes": universe,
        "embryos": out,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    print(f"  embryos: {len(out)}  ({n_e2c} early-2c, {n_l2c} late-2c)   skipped {skipped}")
    print(f"  gene universe: {len(universe)}")
    print(f"  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
