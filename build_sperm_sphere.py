#!/usr/bin/env python3
"""
Build the "Sperm-Entry-Site Enrichment" project data.

For every zygote that has a LABELLED SPERM, test whether each gene's transcripts are
CONCENTRATED (enriched) or DEPLETED inside a sphere of radius r (µm) around the sperm
entry site (SES), across a radius sweep. Everything is in true isotropic µm
(µm = plot × 0.15; transcript z-frame = µm).

This build ships RAW per-SEGMENT counts and voxel volumes so the front-end can
include/exclude segments (cytoplasm, pronuclei, polar body, other) on the fly and
recompute the concentration, fold and a binomial null band for whatever segment set
is selected. Segment of each transcript / voxel comes from the full-resolution label
TIFF; segment identity (pronuclei, polar body) comes from pronuclei_assignments.json.

Per (zygote, gene) we store, per segment s:
    nc[s]        gene transcripts inside segment s (whole cell)          [radius-independent]
    ns[s][r]     those within r µm of the sperm site                     [per radius]
Per zygote, per segment s:
    vc[s]        downsampled voxels of segment s                         [radius-independent]
    vs[s][r]     those within r µm of the sperm                          [per radius]
The front-end then computes, for the selected segment set S:
    concentration_sphere = (Σ_S ns[s][r]) / (Σ_S vs[s][r])
    concentration_cell   = (Σ_S nc[s])    / (Σ_S vc[s])
    fold = concentration_sphere / concentration_cell
    p, 95% band from Binomial(Σ_S nc[s],  Σ_S vs[s][r] / Σ_S vc[s]).

Reuses build_pronuclei (BP) for the atlas/label readers + constants. The 3-D scene is
reused from data/zygote/<id>.json.gz. Output: data/sperm_sphere.json.gz.
Run from the deploy repo root:  python3 build_sperm_sphere.py
"""
import glob
import gzip
import json
import os

import numpy as np
import tifffile

import build_pronuclei as BP
from embryo_naming import embryo_label

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SPERM = os.path.join(DATA, "zygote_sperm.json")
ASSIGN = os.path.join(DATA, "pronuclei_assignments.json")
ZY_MAN = os.path.join(DATA, "zygote_manifest.json")
OUT = os.path.join(DATA, "sperm_sphere.json.gz")

XY_UM, Z_UM = BP.XY_UM, BP.Z_UM           # 0.15, 1.0
DS_XY, DS_Z = BP.DS_XY, BP.DS_Z           # 4, 2
CYTO = BP.CYTO                            # 1
RADII = [5, 8, 10, 12, 15, 18, 20, 24, 28, 32, 36, 40, 45]   # µm sweep (default 20)
MIN_COUNT = 10                           # min in-cell transcripts (all segments) to report a gene

# segment categories, in display order
SEGS = ["cyto", "pron", "polar", "other"]
SEG_META = [
    {"key": "cyto",  "label": "Cytoplasm",  "color": "#64748b"},
    {"key": "pron",  "label": "Pronuclei",  "color": "#a855f7"},
    {"key": "polar", "label": "Polar body", "color": "#f59e0b"},
    {"key": "other", "label": "Other body", "color": "#0d9488"},
]


def transcript_labels_and_um(label_path, tx, genes):
    """Full-res segment label AND µm position for every transcript, per gene."""
    with tifffile.TiffFile(label_path) as t:
        mm = t.series[0].asarray(out="memmap")
        Zn, Yn, Xn = mm.shape
        lens = [len(tx[g]["x"]) for g in genes]
        if sum(lens):
            gx = np.concatenate([np.asarray(tx[g]["x"], float) for g in genes])
            gy = np.concatenate([np.asarray(tx[g]["y"], float) for g in genes])
            gz = np.concatenate([np.asarray(tx[g]["gz"], float) for g in genes])
            ix = np.clip(np.round(gx).astype(np.int64), 0, Xn - 1)
            iy = np.clip(np.round(gy).astype(np.int64), 0, Yn - 1)
            iz = np.clip(np.round(gz).astype(np.int64), 0, Zn - 1)
            labs = np.asarray(mm[iz, iy, ix]).astype(np.int32)
            P = np.stack([gx * XY_UM, gy * XY_UM, gz * Z_UM], axis=1)
        else:
            labs = np.empty(0, np.int32); P = np.empty((0, 3))
        del mm
    out_l, out_p, off = {}, {}, 0
    for g, L in zip(genes, lens):
        out_l[g] = labs[off:off + L]; out_p[g] = P[off:off + L]; off += L
    return out_l, out_p


def seg_category(label_arr, pron_labels, polar_label):
    """Map an int label array to a category string array (cyto/pron/polar/other/'')."""
    cat = np.full(label_arr.shape, "", dtype="<U5")
    cat[label_arr == CYTO] = "cyto"
    for pl in pron_labels:
        cat[label_arr == pl] = "pron"
    if polar_label is not None:
        cat[label_arr == polar_label] = "polar"
    other = (label_arr > 0) & (cat == "")
    cat[other] = "other"
    return cat


def process(eid, sperm_um, pron_labels, polar_label):
    scene_p = os.path.join(BP.ATLAS, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(BP.SRC, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    d = json.load(gzip.open(scene_p, "rt"))
    tx = d.get("transcripts", {})
    genes = [g for g in tx if len(tx[g]["x"])]
    if not genes:
        return None
    labs, pos = transcript_labels_and_um(lab[0], tx, genes)
    sub = BP.load_sub(lab[0])

    # cell voxels (all non-zero) + their segment + µm positions
    iz, iy, ix = np.nonzero(sub > 0)
    if len(iz) == 0:
        return None
    vlab = np.asarray(sub[iz, iy, ix])
    vcat = seg_category(vlab, pron_labels, polar_label)
    Vpos = np.stack([ix * DS_XY * XY_UM, iy * DS_XY * XY_UM, iz * DS_Z * Z_UM], axis=1)
    sp = np.asarray(sperm_um, float)
    dvox = np.linalg.norm(Vpos - sp, axis=1)
    com = Vpos.mean(axis=0)
    r_com = float(np.linalg.norm(sp - com))

    # per-segment cell + sphere VOLUMES (downsampled voxel counts)
    vc = {s: int((vcat == s).sum()) for s in SEGS}
    vs = {s: [int(((vcat == s) & (dvox <= r)).sum()) for r in RADII] for s in SEGS}
    present = [s for s in SEGS if vc[s] > 0]

    gene_out = {}
    for g in genes:
        lg = labs[g]
        gcat = seg_category(lg, pron_labels, polar_label)
        in_cell = gcat != ""
        if int(in_cell.sum()) < MIN_COUNT:
            continue
        Pg = pos[g]
        dsp = np.linalg.norm(Pg - sp, axis=1)
        nc = {s: int((gcat == s).sum()) for s in SEGS}
        ns = {s: [int(((gcat == s) & (dsp <= r)).sum()) for r in RADII] for s in SEGS}
        gene_out[g] = {"nc": [nc[s] for s in SEGS], "ns": [ns[s] for s in SEGS]}

    zs = d.get("z_scale", 7.0)
    return {
        "id": eid, "z_scale": zs,
        "sperm_um": [round(float(x), 2) for x in sp],
        "sperm_plot": [round(float(sp[0] / XY_UM), 2), round(float(sp[1] / XY_UM), 2), round(float(sp[2] * zs), 2)],
        "com_plot": [round(float(com[0] / XY_UM), 2), round(float(com[1] / XY_UM), 2), round(float(com[2] * zs), 2)],
        "r_com_um": round(r_com, 2),
        "present": present,
        "vc": [vc[s] for s in SEGS],
        "vs": [vs[s] for s in SEGS],
        "genes": gene_out,
    }


def main():
    sperm = {e["id"]: e for e in json.load(open(SPERM))["embryos"]}
    assign = {e["id"]: e for e in json.load(open(ASSIGN))["embryos"]}
    zman = {m["id"]: m for m in json.load(open(ZY_MAN))["embryos"]}

    embryos = []
    for eid, se in sorted(sperm.items()):
        sp_plot = se.get("sperm_plot")
        if not sp_plot:
            continue
        sperm_um = np.asarray(sp_plot, float) * XY_UM
        a = assign.get(eid, {})
        pron_labels = [int(p["label"]) for p in (a.get("pron") or []) if p.get("label") is not None]
        polar_label = (a.get("polar") or {}).get("label")
        polar_label = int(polar_label) if polar_label is not None else None
        try:
            r = process(eid, sperm_um, pron_labels, polar_label)
        except Exception as e:                                    # noqa: BLE001
            print(f"  !! {eid}: {e}")
            continue
        if not r:
            print(f"  -- skipped {eid}")
            continue
        r["label"] = (zman.get(eid, {}).get("label")) or embryo_label(eid, "zygote") or eid
        r["date_short"] = zman.get(eid, {}).get("date_short", "")
        embryos.append(r)
        print(f"  {eid}  {r['label']}  {len(r['genes'])} genes  segs={r['present']}  r(sperm→COM)={r['r_com_um']}µm")

    doc = {
        "radii": RADII, "segs": SEGS, "segMeta": SEG_META, "embryos": embryos,
        "meta": {"nZygotes": len(embryos), "minCount": MIN_COUNT,
                 "unit_um_per_plot": XY_UM, "defaultRadiusIdx": RADII.index(20)},
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    ngenes = len({g for e in embryos for g in e["genes"]})
    print(f"\nwrote {len(embryos)} sperm zygotes · {ngenes} genes ({os.path.getsize(OUT)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
