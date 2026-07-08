#!/usr/bin/env python3
"""
Build the "Segment Gene Enrichment" project data — EVERY embryo in the atlas
(Oocyte, Zygote, Early-2-cell, Late-2-cell).

For each embryo we assign every transcript to a segmentation label (the cell body,
polar body, pronuclei, nuclei, … — the meaning of a label varies by stage, so we
treat them generically as "Segment 1..N"), and rank the genes by how ENRICHED they
are in each segment:

    enrichment(g, s) = (n_{g,s} / V_s) / (n_g / V_tot)
                     = density of gene g in segment s  ÷  its embryo-wide density

where n_{g,s} = # transcripts of g in segment s, V_s = segment volume, n_g = # of g
in ANY segment, V_tot = total segmented volume. Fold > 1 ⇒ concentrated in s.
Genes with < MIN_COUNT transcripts in a segment are dropped (single-molecule noise).

Outputs: data/segments/<stage>__<id>.json.gz (slim scene: rounded region meshes +
per-segment ranked gene lists) and data/segments_manifest.json.
"""
import glob
import gzip
import json
import os
from collections import Counter

import numpy as np
import tifffile

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = "/Users/rishib/Desktop/MERFISH/Website2/MerfishAtlasWebsite/public/data"
SRC = os.path.join(HERE, "..", "TranscriptomicsData", "JustTifAndCSVData")
OUT_DIR = os.path.join(HERE, "public", "data", "segments")
OUT_MANIFEST = os.path.join(HERE, "public", "data", "segments_manifest.json")

STAGES = ["Oocyte", "Zygote", "Early2Cell", "Late2Cell"]
STAGE_LABEL = {"Oocyte": "Oocyte", "Zygote": "Zygote", "Early2Cell": "Early 2-cell",
               "Late2Cell": "Late 2-cell"}
XY_UM = 0.15
Z_UM = 1.0
DS_XY = 6
DS_Z = 2
MIN_COUNT = 3          # drop genes with fewer than this many transcripts in a segment


def mask_and_transcripts(label_path, tx, genes):
    """Downsampled nonzero voxel positions (µm) + segment labels + per-voxel volume,
    and `seg_of`: the segment label under every transcript (full-res lookup)."""
    with tifffile.TiffFile(label_path) as t:
        mm = t.series[0].asarray(out="memmap")            # (Z, Y, X)
        Zn, Yn, Xn = mm.shape
        sub = np.asarray(mm[::DS_Z, ::DS_XY, ::DS_XY])
        lens = [len(tx[g]["x"]) for g in genes]
        seg_of = {}
        if sum(lens):
            gx = np.concatenate([np.asarray(tx[g]["x"], float) for g in genes])
            gy = np.concatenate([np.asarray(tx[g]["y"], float) for g in genes])
            gzf = np.concatenate([np.asarray(tx[g]["gz"], float) for g in genes])
            ix = np.clip(np.round(gx).astype(np.int64), 0, Xn - 1)
            iy = np.clip(np.round(gy).astype(np.int64), 0, Yn - 1)
            iz = np.clip(np.round(gzf).astype(np.int64), 0, Zn - 1)
            labs = np.asarray(mm[iz, iy, ix]).astype(np.int16)
        else:
            labs = np.empty(0, np.int16)
        off = 0
        for g, L in zip(genes, lens):
            seg_of[g] = labs[off:off + L]; off += L
        del mm
    iz, iy, ix = np.nonzero(sub)
    labels = sub[iz, iy, ix].astype(np.int16)
    voxel_vol = (DS_XY * XY_UM) ** 2 * (DS_Z * Z_UM)
    return labels, voxel_vol, seg_of


def round_mesh(m):
    """Round mesh vertices to 0.1 plot-units (sub-pixel) so they gzip well."""
    return {"verts": [round(float(v), 1) for v in m["verts"]], "faces": m["faces"]}


def process(stage, eid):
    scene_p = os.path.join(ATLAS, stage, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(SRC, stage, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    d = json.load(gzip.open(scene_p, "rt"))
    tx = d.get("transcripts", {})
    genes = sorted(tx.keys())
    if not genes:
        return None
    labels, voxvol, seg_of = mask_and_transcripts(lab[0], tx, genes)
    seg_ids = sorted(int(s) for s in np.unique(labels) if s >= 1)
    if not seg_ids:
        return None
    vol = {s: float((labels == s).sum()) * voxvol for s in seg_ids}
    v_tot = sum(vol.values()) or 1.0

    # per-segment gene counts + per-gene total-in-any-segment
    seg_counts = {s: Counter() for s in seg_ids}
    gene_seg_total = Counter()
    for g in genes:
        segs = seg_of[g]
        if len(segs) == 0:
            continue
        inseg = segs[segs >= 1]
        if len(inseg) == 0:
            continue
        gene_seg_total[g] = int(len(inseg))
        vals, cnts = np.unique(inseg, return_counts=True)
        for s, c in zip(vals.tolist(), cnts.tolist()):
            seg_counts[int(s)][g] = int(c)

    ranked = {}
    for s in seg_ids:
        rows = []
        for g, c in seg_counts[s].items():
            if c < MIN_COUNT:
                continue
            ng = gene_seg_total[g]
            enrich = (c / vol[s]) / (ng / v_tot) if ng else 0.0
            rows.append({"gene": g, "enrich": round(enrich, 3), "count": c, "ntot": ng})
        rows.sort(key=lambda r: -r["enrich"])
        ranked[str(s)] = rows

    rm = d.get("region_meshes", {})
    scene = {
        "id": eid, "stage": stage, "z_scale": d.get("z_scale", 7.0),
        "extents": d["extents"], "mask_labels": seg_ids,
        "region_defaults": d.get("region_defaults", {}),
        "region_meshes": {str(s): round_mesh(rm[str(s)]) for s in seg_ids if str(s) in rm},
        "segments": [{"label": s, "volume": round(vol[s], 1)} for s in seg_ids],
        "ranked": ranked,
    }
    return scene


def _json_default(o):
    if isinstance(o, np.floating):
        return round(float(o), 6)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not serializable: {type(o)}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    manifest = []
    for stage in STAGES:
        sdir = os.path.join(ATLAS, stage)
        if not os.path.isdir(sdir):
            continue
        ids = sorted(e for e in os.listdir(sdir) if not e.startswith("."))
        for i, eid in enumerate(ids):
            try:
                scene = process(stage, eid)
            except Exception as e:            # noqa: BLE001
                print(f"  !! {stage}/{eid}: {e}"); continue
            if not scene:
                print(f"  -- skipped {stage}/{eid}"); continue
            key = f"{stage}__{eid}"
            out = os.path.join(OUT_DIR, key + ".json.gz")
            with gzip.open(out, "wt") as fh:
                json.dump(scene, fh, separators=(",", ":"), default=_json_default)
            manifest.append({
                "id": key, "stage": stage, "stage_label": STAGE_LABEL.get(stage, stage),
                "eid": eid, "label": short_label(eid),
                "n_segments": len(scene["mask_labels"]), "segments": scene["mask_labels"],
                "size_kb": round(os.path.getsize(out) / 1024),
            })
            print(f"  [{stage} {i+1}/{len(ids)}] {eid}  segs={scene['mask_labels']}  "
                  f"{manifest[-1]['size_kb']}KB")
    with open(OUT_MANIFEST, "w") as fh:
        json.dump({"embryos": manifest, "stages": STAGES,
                   "stage_labels": STAGE_LABEL, "min_count": MIN_COUNT}, fh, indent=1)
    tot = sum(m["size_kb"] for m in manifest)
    print(f"\nwrote {len(manifest)} embryos across {len(STAGES)} stages  ({tot/1024:.1f} MB)")


def short_label(eid):
    """Compact nav label from an atlas id (drops the date prefix)."""
    import re
    s = re.sub(r"^\d{8}_", "", eid)
    s = re.sub(r"^(zygote|oocyte|l2c|e2c|late2cell|early2cell)_?", "", s, flags=re.I)
    s = s.replace("sample", "s").replace("_", " ")
    return s.strip() or eid


if __name__ == "__main__":
    main()
