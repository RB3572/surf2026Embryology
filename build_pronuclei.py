#!/usr/bin/env python3
"""
Build the "Pronuclei Distance vs Transcripts" project data — ZYGOTES only.

A fertilised zygote has two pronuclei (male + female). Their segmentation LABEL
NUMBERS are NOT consistent (sometimes 3 & 4, sometimes 2 & 3, sometimes 4 & 5, …),
so we auto-detect them: segment 1 is always the cytoplasm, and the pronuclei are the
segments that sit INSIDE it. A segment is "inside seg 1" when more of its dilation
shell borders seg 1 than borders background (label 0) — this cleanly excludes the
polar body / perivitelline debris, which border the outside. Among the inside
segments we take the two largest by volume as the pronuclei. For every zygote with
>= 2 such inside segments we compute:
  * the minimum distance between the two pronuclei (the shortest line between their
    voxel clouds), in µm, via a KD-tree; and
  * the total number of detected transcripts in the zygote.

The bottom-drawer scatter plots distance vs total transcripts (with a regression).
The 3-D view shows the cytoplasm + the two pronuclei and the min-distance line.

Outputs: data/pronuclei/<id>.json.gz (slim scene + the two closest points) and
data/pronuclei_manifest.json (the scatter points: id, label, distance, total).
"""
import glob
import gzip
import json
import os

import numpy as np
from embryo_naming import embryo_label
import tifffile
from scipy.spatial import cKDTree
from scipy.ndimage import binary_dilation

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = "/Users/rishib/Desktop/MERFISH/Website2/MerfishAtlasWebsite/public/data/Zygote"
SRC = os.path.join(HERE, "..", "TranscriptomicsData", "JustTifAndCSVData", "Zygote")
OUT_DIR = os.path.join(HERE, "public", "data", "pronuclei")
OUT_MANIFEST = os.path.join(HERE, "public", "data", "pronuclei_manifest.json")
OUT_GENES = os.path.join(HERE, "public", "data", "pronuclei_genes.json.gz")

XY_UM = 0.15
Z_UM = 1.0
DS_XY = 4          # finer than the other builds — the pronuclei are small
DS_Z = 2
CYTO = 1           # segment 1 is always the cytoplasm
MAX_DOTS = 1000    # cap stored transcript dots per gene (subsampled) for file size


def load_sub(label_path):
    """Downsampled label volume (Z, Y, X)."""
    with tifffile.TiffFile(label_path) as t:
        mm = t.series[0].asarray(out="memmap")
        sub = np.asarray(mm[::DS_Z, ::DS_XY, ::DS_XY]); del mm
    return sub


def detect_pronuclei(sub):
    """The two pronuclei = the two largest segments (by volume) that sit INSIDE the
    cytoplasm (seg 1). A segment is inside when its dilation shell borders seg 1 more
    than it borders background (label 0). Returns (labelA, labelB) largest-first, or
    None if fewer than two inside segments exist."""
    labs = [int(v) for v in np.unique(sub) if v >= 1]
    if CYTO not in labs:
        return None
    seg1 = (sub == CYTO)
    bg = (sub == 0)
    inside = []
    for s in labs:
        if s == CYTO:
            continue
        m = (sub == s)
        shell = binary_dilation(m, iterations=2) & ~m
        n1 = int((shell & seg1).sum())
        n0 = int((shell & bg).sum())
        if n1 > n0:                                   # more cytoplasm than background around it
            inside.append((int(m.sum()), s))
    inside.sort(reverse=True)                          # largest volume first
    if len(inside) < 2:
        return None
    return inside[0][1], inside[1][1]


def coords_of(sub, label):
    """Voxel positions (µm) for one segment label."""
    iz, iy, ix = np.nonzero(sub == label)
    return np.stack([ix * DS_XY * XY_UM, iy * DS_XY * XY_UM, iz * DS_Z * Z_UM], axis=1).astype(np.float32)


def round_mesh(m):
    return {"verts": [round(float(v), 1) for v in m["verts"]], "faces": m["faces"]}


def um_to_plot(p, zs):
    return [round(p[0] / XY_UM, 2), round(p[1] / XY_UM, 2), round(p[2] * zs, 2)]


def process(eid):
    scene_p = os.path.join(ATLAS, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(SRC, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    sub = load_sub(lab[0])
    pron = detect_pronuclei(sub)
    if not pron:
        return None                                    # fewer than two pronuclei inside seg 1
    la, lb = pron
    A = coords_of(sub, la)
    B = coords_of(sub, lb)
    tree = cKDTree(B)
    dist, idx = tree.query(A)                           # nearest B for each A
    j = int(np.argmin(dist))
    d_min = float(dist[j])
    pa = A[j]; pb = B[int(idx[j])]

    d = json.load(gzip.open(scene_p, "rt"))
    tx = d.get("transcripts", {})
    total = int(d.get("n_transcripts") or sum(len(t["x"]) for t in tx.values()))
    # per-gene transcript count in this zygote (for the gene↔distance correlation)
    gene_counts = {g: int(v) for g, v in d.get("gene_totals", {}).items() if v}
    if not gene_counts:
        gene_counts = {g: len(t["x"]) for g, t in tx.items() if len(t["x"])}
    # per-gene transcript locations for the 3-D dots (subsampled to MAX_DOTS)
    tx_dots = {}
    for g, t in tx.items():
        n = len(t["x"])
        if n == 0:
            continue
        if n > MAX_DOTS:
            idx = np.linspace(0, n - 1, MAX_DOTS).astype(int)
            tx_dots[g] = {"x": [t["x"][k] for k in idx], "y": [t["y"][k] for k in idx], "gz": [t["gz"][k] for k in idx]}
        else:
            tx_dots[g] = {"x": t["x"], "y": t["y"], "gz": t["gz"]}
    zs = d.get("z_scale", 7.0)
    rm = d.get("region_meshes", {})
    seg_ids = sorted(int(s) for s in np.unique(sub) if s >= 1)
    scene = {
        "id": eid, "z_scale": zs, "extents": d["extents"],
        "mask_labels": seg_ids, "region_defaults": d.get("region_defaults", {}),
        "region_meshes": {str(s): round_mesh(rm[str(s)]) for s in seg_ids if str(s) in rm},
        "pron_labels": [la, lb],                        # [larger, smaller] by volume; auto-detected
        "line_plot": [um_to_plot(pa, zs), um_to_plot(pb, zs)],
        "distance_um": round(d_min, 2), "total_transcripts": total,
        "transcripts": tx_dots,                         # per-gene dot locations (x px, y px, gz frame)
    }
    return scene, d_min, total, gene_counts


def _json_default(o):
    if isinstance(o, np.floating):
        return round(float(o), 6)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not serializable: {type(o)}")


def short_label(eid):
    return embryo_label(eid, "zygote")


def date_short(eid):
    date = eid[:8]
    MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    return f"{MON[int(date[4:6]) - 1]} {int(date[6:8])}" if date.isdigit() else ""


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    ids = sorted(e for e in os.listdir(ATLAS) if not e.startswith("."))
    points, gene_agg = [], []
    for i, eid in enumerate(ids):
        try:
            res = process(eid)
        except Exception as e:            # noqa: BLE001
            print(f"  !! {eid}: {e}"); continue
        if not res:
            print(f"  -- skipped {eid} (needs both pronuclei)"); continue
        scene, dist, total, gene_counts = res
        out = os.path.join(OUT_DIR, eid + ".json.gz")
        with gzip.open(out, "wt") as fh:
            json.dump(scene, fh, separators=(",", ":"), default=_json_default)
        points.append({"id": eid, "label": short_label(eid), "date_short": date_short(eid),
                       "distance": round(dist, 2), "total": total, "pron_labels": scene["pron_labels"],
                       "size_kb": round(os.path.getsize(out) / 1024)})
        gene_agg.append({"id": eid, "distance": round(dist, 2), "total": total, "genes": gene_counts})
        print(f"  [{i+1}/{len(ids)}] {eid}  pronuclei=seg{scene['pron_labels']}  dist={dist:.1f}um  total={total}")
    points.sort(key=lambda p: p["id"]); gene_agg.sort(key=lambda p: p["id"])
    with open(OUT_MANIFEST, "w") as fh:
        json.dump({"embryos": points}, fh, indent=1)
    # per-gene per-zygote counts for the gene↔distance correlation (loaded once by the UI)
    with gzip.open(OUT_GENES, "wt") as fh:
        json.dump({"embryos": gene_agg}, fh, separators=(",", ":"), default=_json_default)
    tot = sum(p["size_kb"] for p in points)
    print(f"\nwrote {len(points)} zygotes with 2 pronuclei  ({tot/1024:.1f} MB)  + gene aggregate "
          f"({os.path.getsize(OUT_GENES)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
