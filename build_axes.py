#!/usr/bin/env python3
"""
Build the "Fertilization Geometry & the First Embryonic Axis" project.

Tests the classic pre-patterning question (Piotrowska & Zernicka-Goetz 2001 vs.
Hiiragi & Solter 2004 / Motosugi 2005): does the geometry set up at fertilization —
the sperm entry position (SEP) and the polar body / animal pole — predict the first
cleavage plane and hence the embryonic axis?

Cohorts (the 45 sperm-positive embryos from data/sperm_transcriptomics/embryos.csv):
  * ZYGOTE (26): cleavage not yet occurred → we report the *predictive* geometry
    (animal-vegetal axis, sperm axis, pronuclear axis).
  * TWO-CELL (19): the first cleavage plane is physically realised as the interface
    between the two blastomeres → we run the direct test.

Landmarks (auto-detected, per stage; validated on samples):
  * cell bodies = the largest segment(s): 1 for a zygote (cytoplasm), 2 for a
    two-cell (the two blastomeres).
  * nuclei/pronuclei = segments INSIDE a cell body (dilation shell borders a cell
    body more than background).
  * polar body / animal pole = the largest segment OUTSIDE the cell bodies (shell
    borders background more). The polar body is the accepted animal-pole landmark.
  * sperm = the manually-labelled GFP midpiece coordinate from embryos.csv. The
    midpiece stays near the fusion/entry site (the head decondenses into the
    pronucleus and migrates), so it is a proxy for the sperm entry position.

Axes & angles (physical µm; axes are undirected → we use |cos θ| of the acute angle):
  * AV axis        = COM → polar body.
  * sperm axis     = COM → sperm midpiece.
  * pronuclear axis (zygote) = pronucleus_A → pronucleus_B  (Hiiragi's proposed
    determinant; also the *predicted* cleavage-plane normal).
  * cleavage normal (two-cell) = blastomere_A centroid → blastomere_B centroid; the
    cleavage plane is perpendicular to it, through the interface midpoint.
  * "the plane contains landmark X" ⇔ (COM→X) ⟂ normal ⇔ |cos θ| ≈ 0.
  * SHAPE CONTROL: the embryo's own longest principal axis (PCA of the mask). If the
    cleavage plane merely follows embryo elongation, the AV correlation is a shape
    artefact (Motosugi's critique) — so we store the plane-vs-long-axis and
    landmark-vs-long-axis angles too, for transparency.

Molecular layer: every transcript is assigned a SIDE. Two-cell: 0 = the blastomere
that inherits the sperm midpiece (nearest body centroid), 1 = the other → tests the
P&ZG claim that the sperm-associated blastomere has a distinct identity. Zygote:
0 = animal half, 1 = vegetal half of the AV axis → tests an AV transcriptome gradient.

Outputs: data/axes/<id>.json.gz (meshes + landmarks + axes + per-embryo angles +
transcripts with side), data/axes_manifest.json (all per-embryo angles for the
aggregate figure), data/axes_molecular.json.gz (per-gene per-embryo side counts).
"""
import csv
import glob
import gzip
import json
import os
from collections import Counter

import numpy as np
from embryo_naming import embryo_label
import tifffile
from scipy.ndimage import binary_dilation

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = "/Users/rishib/Desktop/MERFISH/Website2/MerfishAtlasWebsite/public/data"
SRC = os.path.join(HERE, "..", "TranscriptomicsData", "JustTifAndCSVData")
CSV = os.path.join(HERE, "..", "data", "sperm_transcriptomics", "embryos.csv")
OUT_DIR = os.path.join(HERE, "public", "data", "axes")
OUT_MANIFEST = os.path.join(HERE, "public", "data", "axes_manifest.json")
OUT_MOLECULAR = os.path.join(HERE, "public", "data", "axes_molecular.json.gz")

XY_UM = 0.15
Z_UM = 1.0
DS_XY = 6
DS_Z = 2
MAX_DOTS = 1200        # cap stored dots per gene (subsampled) for file size; counts stay exact
STAGE_DIR = {"Zygote": "Zygote", "2-cell (early)": "Early2Cell", "2-cell (late)": "Late2Cell"}


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def load_sub(label_path):
    with tifffile.TiffFile(label_path) as t:
        mm = t.series[0].asarray(out="memmap")
        sub = np.asarray(mm[::DS_Z, ::DS_XY, ::DS_XY]); del mm
    return sub


def seg_pos(sub, L):
    iz, iy, ix = np.nonzero(sub == L)
    return np.stack([ix * DS_XY * XY_UM, iy * DS_XY * XY_UM, iz * DS_Z * Z_UM], axis=1)


def detect_landmarks(sub, n_bodies):
    """n_bodies = 1 (zygote) or 2 (two-cell). Returns cell-body labels, nucleus
    labels (inside), and the polar-body label (largest outside), all by volume."""
    labs = [int(v) for v in np.unique(sub) if v >= 1]
    if not labs:
        return None
    vols = {L: int((sub == L).sum()) for L in labs}
    bodies = sorted(labs, key=lambda L: -vols[L])[:n_bodies]
    cell_mask = np.isin(sub, bodies)
    bg = (sub == 0)
    inside, outside = [], []
    for L in labs:
        if L in bodies:
            continue
        m = (sub == L)
        shell = binary_dilation(m, iterations=2) & ~m
        n_in = int((shell & cell_mask).sum()); n_bg = int((shell & bg).sum())
        (inside if n_in > n_bg else outside).append((vols[L], L))
    inside.sort(reverse=True); outside.sort(reverse=True)
    return {"bodies": bodies, "nuclei": [L for _, L in inside],
            "polar": (outside[0][1] if outside else None)}


def principal_axes(pos, com):
    """PCA of the mask voxels → (long-axis unit vector, aspect ratio = sqrt(λ1/λ3))."""
    d = pos - com
    cov = (d.T @ d) / len(d)
    w, V = np.linalg.eigh(cov)                       # ascending eigenvalues
    long_axis = unit(V[:, -1])
    aspect = float(np.sqrt(w[-1] / max(w[0], 1e-9)))
    return long_axis, aspect


def acos_deg(c):
    return float(np.degrees(np.arccos(np.clip(abs(c), 0, 1))))


def um_to_plot(p, zs):
    return [round(float(p[0]) / XY_UM, 1), round(float(p[1]) / XY_UM, 1), round(float(p[2]) * zs, 1)]


def round_mesh(m):
    return {"verts": [round(float(v), 1) for v in m["verts"]], "faces": m["faces"]}


def process(row):
    eid, stage = row["embryo_id"], row["stage"]
    sdir = STAGE_DIR.get(stage)
    if not sdir:
        return None
    scene_p = os.path.join(ATLAS, sdir, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(SRC, sdir, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    sub = load_sub(lab[0])
    is_zyg = (stage == "Zygote")
    lm = detect_landmarks(sub, 1 if is_zyg else 2)
    if not lm:
        return None

    allpos = np.concatenate([seg_pos(sub, L) for L in [int(v) for v in np.unique(sub) if v >= 1]])
    com = allpos.mean(axis=0)
    long_axis, aspect = principal_axes(allpos, com)

    body_c = [seg_pos(sub, L).mean(axis=0) for L in lm["bodies"]]
    nuc_c = [seg_pos(sub, L).mean(axis=0) for L in lm["nuclei"]]
    pb = seg_pos(sub, lm["polar"]).mean(axis=0) if lm["polar"] is not None else None
    sp = np.array([float(row["sperm_x_um"]), float(row["sperm_y_um"]), float(row["sperm_z_um"])])

    av = unit(pb - com) if pb is not None else None
    sperm_ax = unit(sp - com)

    d = json.load(gzip.open(scene_p, "rt"))
    zs = d.get("z_scale", 7.0)
    rm = d.get("region_meshes", {})
    seg_ids = sorted(int(s) for s in np.unique(sub) if s >= 1)
    tx = d.get("transcripts", {})
    genes = sorted(tx.keys())

    angles = {"aspect_ratio": round(aspect, 3)}
    if is_zyg:
        normal = unit(nuc_c[0] - nuc_c[1]) if len(nuc_c) >= 2 else None   # pronuclear axis
    else:
        normal = unit(body_c[1] - body_c[0])                              # cleavage normal
    if normal is not None:
        if av is not None:
            angles["plane_vs_pb"] = {"deg": acos_deg(av @ normal), "cos": round(abs(float(av @ normal)), 4)}
        angles["plane_vs_sperm"] = {"deg": acos_deg(sperm_ax @ normal), "cos": round(abs(float(sperm_ax @ normal)), 4)}
        angles["plane_vs_shape"] = {"deg": acos_deg(long_axis @ normal), "cos": round(abs(float(long_axis @ normal)), 4)}
    if av is not None:
        angles["pb_vs_shape"] = {"deg": acos_deg(av @ long_axis), "cos": round(abs(float(av @ long_axis)), 4)}
        angles["av_vs_sperm"] = {"deg": acos_deg(av @ sperm_ax), "cos": round(abs(float(av @ sperm_ax)), 4)}

    # ---- molecular side assignment ----
    sperm_body = None
    if not is_zyg and len(body_c) == 2:
        sperm_body = int(np.argmin([np.linalg.norm(sp - c) for c in body_c]))   # blastomere nearest the sperm
    gene_sides = {}
    tx_out = {}
    for g in genes:
        t = tx[g]
        n = len(t["x"])
        if n == 0:
            gene_sides[g] = [0, 0]; tx_out[g] = {"x": [], "y": [], "gz": [], "side": []}; continue
        P = np.stack([np.asarray(t["x"], float) * XY_UM, np.asarray(t["y"], float) * XY_UM,
                      np.asarray(t["gz"], float) * Z_UM], axis=1)
        if not is_zyg and len(body_c) == 2:
            dA = np.linalg.norm(P - body_c[0], axis=1); dB = np.linalg.norm(P - body_c[1], axis=1)
            near = (dB < dA).astype(int)                      # 0 = body0, 1 = body1
            side = (near != sperm_body).astype(int)           # 0 = sperm blastomere, 1 = other
        elif av is not None:
            side = ((P - com) @ av <= 0).astype(int)          # 0 = animal half, 1 = vegetal half
        else:
            side = np.zeros(n, int)
        gene_sides[g] = [int((side == 0).sum()), int((side == 1).sum())]      # exact, from all transcripts
        if n > MAX_DOTS:                                                       # subsample dots for file size
            idx = np.linspace(0, n - 1, MAX_DOTS).astype(int)
            xs = [t["x"][k] for k in idx]; ys = [t["y"][k] for k in idx]
            gzs = [t["gz"][k] for k in idx]; sd = side[idx]
        else:
            xs, ys, gzs, sd = t["x"], t["y"], t["gz"], side
        tx_out[g] = {"x": xs, "y": ys, "gz": gzs, "s": np.asarray(sd, np.uint8).tolist()}

    scene = {
        "id": eid, "stage": "zygote" if is_zyg else "twocell", "stage_label": stage,
        "z_scale": zs, "extents": d["extents"], "mask_labels": seg_ids,
        "region_defaults": d.get("region_defaults", {}),
        "region_meshes": {str(s): round_mesh(rm[str(s)]) for s in seg_ids if str(s) in rm},
        "com_plot": um_to_plot(com, zs),
        "landmarks": {
            "polar_plot": um_to_plot(pb, zs) if pb is not None else None,
            "body_plots": [um_to_plot(c, zs) for c in body_c],
            "nuclei_plots": [um_to_plot(c, zs) for c in nuc_c],
            "sperm_plot": um_to_plot(sp, zs), "sperm_body": sperm_body,
        },
        "long_axis_plot": um_to_plot(com + long_axis * 30, zs),   # a point 30µm along the long axis (for drawing)
        "normal_kind": "pronuclear axis" if is_zyg else "cleavage-plane normal",
        "angles": angles,
        "transcripts": tx_out,
    }
    return scene, angles, gene_sides, sperm_body


def _json_default(o):
    if isinstance(o, np.floating):
        return round(float(o), 6)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not serializable: {type(o)}")


def short_label(eid):
    return embryo_label(eid)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rows = list(csv.DictReader(open(CSV)))
    manifest, molecular = [], []
    for i, row in enumerate(rows):
        try:
            res = process(row)
        except Exception as e:              # noqa: BLE001
            print(f"  !! {row['embryo_id']}: {e}"); continue
        if not res:
            print(f"  -- skipped {row['embryo_id']}"); continue
        scene, angles, gene_sides, sperm_body = res
        out = os.path.join(OUT_DIR, scene["id"] + ".json.gz")
        with gzip.open(out, "wt") as fh:
            json.dump(scene, fh, separators=(",", ":"), default=_json_default)
        manifest.append({
            "id": scene["id"], "stage": scene["stage"], "stage_label": scene["stage_label"],
            "label": short_label(scene["id"]), "angles": angles,
            "has_pb": scene["landmarks"]["polar_plot"] is not None,
            "size_kb": round(os.path.getsize(out) / 1024),
        })
        molecular.append({"id": scene["id"], "stage": scene["stage"], "sides": gene_sides})
        pv = angles.get("plane_vs_pb", {}).get("deg")
        ps = angles.get("plane_vs_sperm", {}).get("deg")
        print(f"  [{i+1}/{len(rows)}] {scene['id']} ({scene['stage']})  "
              f"plane·PB={pv:.0f}°  plane·sperm={ps:.0f}°" if pv else f"  [{i+1}] {scene['id']} ({scene['stage']})")
    with open(OUT_MANIFEST, "w") as fh:
        json.dump({"embryos": manifest, "n_zygote": sum(1 for m in manifest if m["stage"] == "zygote"),
                   "n_twocell": sum(1 for m in manifest if m["stage"] == "twocell")}, fh, indent=1)
    with gzip.open(OUT_MOLECULAR, "wt") as fh:
        json.dump({"embryos": molecular}, fh, separators=(",", ":"), default=_json_default)
    print(f"\nwrote {len(manifest)} embryos  "
          f"({sum(m['size_kb'] for m in manifest)/1024:.1f} MB)")


if __name__ == "__main__":
    main()
