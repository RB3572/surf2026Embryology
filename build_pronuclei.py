#!/usr/bin/env python3
"""
Build the "Pronuclei Distance vs Transcripts" project data — ZYGOTES only.

A fertilised zygote has two pronuclei (male + female), segmented as labels 3 and 4
(label 1 = cytoplasm, 2 = polar body). For every zygote that has BOTH pronuclei we
compute:
  * the minimum distance between the two pronuclei (the shortest line that can be
    drawn between the seg-3 and seg-4 voxel clouds), in µm, via a KD-tree; and
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
import tifffile
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = "/Users/rishib/Desktop/MERFISH/Website2/MerfishAtlasWebsite/public/data/Zygote"
SRC = os.path.join(HERE, "..", "TranscriptomicsData", "JustTifAndCSVData", "Zygote")
OUT_DIR = os.path.join(HERE, "public", "data", "pronuclei")
OUT_MANIFEST = os.path.join(HERE, "public", "data", "pronuclei_manifest.json")

XY_UM = 0.15
Z_UM = 1.0
DS_XY = 4          # finer than the other builds — the pronuclei are small
DS_Z = 2
PRON = (3, 4)      # the two pronuclei labels


def mask_coords(label_path):
    """Downsampled nonzero voxel positions (µm) + their segment labels."""
    with tifffile.TiffFile(label_path) as t:
        mm = t.series[0].asarray(out="memmap")
        sub = np.asarray(mm[::DS_Z, ::DS_XY, ::DS_XY]); del mm
    iz, iy, ix = np.nonzero(sub)
    labels = sub[iz, iy, ix].astype(np.int16)
    pos = np.stack([ix * DS_XY * XY_UM, iy * DS_XY * XY_UM, iz * DS_Z * Z_UM], axis=1)
    return pos.astype(np.float32), labels


def round_mesh(m):
    return {"verts": [round(float(v), 1) for v in m["verts"]], "faces": m["faces"]}


def um_to_plot(p, zs):
    return [round(p[0] / XY_UM, 2), round(p[1] / XY_UM, 2), round(p[2] * zs, 2)]


def process(eid):
    scene_p = os.path.join(ATLAS, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(SRC, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    pos, labels = mask_coords(lab[0])
    if not (PRON[0] in labels and PRON[1] in labels):
        return None                                    # need BOTH pronuclei
    A = pos[labels == PRON[0]]
    B = pos[labels == PRON[1]]
    tree = cKDTree(B)
    dist, idx = tree.query(A)                           # nearest B for each A
    j = int(np.argmin(dist))
    d_min = float(dist[j])
    pa = A[j]; pb = B[int(idx[j])]

    d = json.load(gzip.open(scene_p, "rt"))
    total = int(d.get("n_transcripts") or sum(len(t["x"]) for t in d.get("transcripts", {}).values()))
    zs = d.get("z_scale", 7.0)
    rm = d.get("region_meshes", {})
    seg_ids = sorted(int(s) for s in np.unique(labels) if s >= 1)
    scene = {
        "id": eid, "z_scale": zs, "extents": d["extents"],
        "mask_labels": seg_ids, "region_defaults": d.get("region_defaults", {}),
        "region_meshes": {str(s): round_mesh(rm[str(s)]) for s in seg_ids if str(s) in rm},
        "pron_labels": list(PRON),
        "line_plot": [um_to_plot(pa, zs), um_to_plot(pb, zs)],
        "distance_um": round(d_min, 2), "total_transcripts": total,
    }
    return scene, d_min, total


def _json_default(o):
    if isinstance(o, np.floating):
        return round(float(o), 6)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not serializable: {type(o)}")


def short_label(eid):
    import re
    s = re.sub(r"^\d{8}_", "", eid)
    s = re.sub(r"^zygote_?", "", s, flags=re.I)
    return (s.replace("sample", "s").replace("_", " ").strip()) or eid


def date_short(eid):
    date = eid[:8]
    MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    return f"{MON[int(date[4:6]) - 1]} {int(date[6:8])}" if date.isdigit() else ""


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    ids = sorted(e for e in os.listdir(ATLAS) if not e.startswith("."))
    points = []
    for i, eid in enumerate(ids):
        try:
            res = process(eid)
        except Exception as e:            # noqa: BLE001
            print(f"  !! {eid}: {e}"); continue
        if not res:
            print(f"  -- skipped {eid} (needs both pronuclei)"); continue
        scene, dist, total = res
        out = os.path.join(OUT_DIR, eid + ".json.gz")
        with gzip.open(out, "wt") as fh:
            json.dump(scene, fh, separators=(",", ":"), default=_json_default)
        points.append({"id": eid, "label": short_label(eid), "date_short": date_short(eid),
                       "distance": round(dist, 2), "total": total,
                       "size_kb": round(os.path.getsize(out) / 1024)})
        print(f"  [{i+1}/{len(ids)}] {eid}  dist={dist:.1f}um  total={total}")
    points.sort(key=lambda p: p["id"])
    with open(OUT_MANIFEST, "w") as fh:
        json.dump({"embryos": points, "pron_labels": list(PRON)}, fh, indent=1)
    tot = sum(p["size_kb"] for p in points)
    print(f"\nwrote {len(points)} zygotes with 2 pronuclei  ({tot/1024:.1f} MB)")


if __name__ == "__main__":
    main()
