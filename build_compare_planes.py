#!/usr/bin/env python3
"""
Build data/compare_planes.json.gz — the "Compare Division Planes" project.

For every zygote, four candidate dividing planes are compared HEAD-TO-HEAD under ONE
consistent metric:

  1. polar      — best plane CONTAINING the polar-body axis (the 18-plane fan, 10° steps),
                  chosen per gene by the concentration asymmetry below.
  2. exhaustive — best plane from the every-plane search (reused from build_planes_all:
                  each gene's stored density-best normal over the ~20k-orientation grid).
  3. equatorial — the plane whose normal IS the polar-body (animal–vegetal) axis.
  4. sperm      — the plane through the sperm, the cell centre of mass and the polar-body
                  centroid (null for zygotes with no labelled sperm).

The metric is a CONCENTRATION asymmetry on the CELL BODY ONLY (segment 1) — pronuclei,
polar body and any other segment are excluded from BOTH the transcript counts and the
side volumes, exactly as the user asked. Per (zygote, gene, plane) we ship the raw
per-side counts and per-side cell-body volumes; the front-end derives, live:

    concentration_side = count_side / volume_side           (transcripts · µm⁻³)
    asymmetry = |cA − cB| / (cA + cB)   ∈ [0, 1]            (0 = balanced, 1 = one-sided)
    p = two-sided Binomial(n = cA+cB, f = vA/(vA+vB))       (asymmetry vs a uniform null)

Fixed planes (equatorial, sperm) have one normal + one volume split per zygote; the
gene-dependent planes (polar, exhaustive) carry their own normal + volume split per gene.

Reuses build_zygote's readers (mask_and_transcripts, detect_polar_body) and
build_planes_all's geometry (side_counts, in_plane_basis, the normal grid). µm space is
isotropic: x,y × 0.15, z(frame) × 1.0.  Run from the deploy repo root:
    python3 build_compare_planes.py
"""
import glob
import gzip
import json
import os

import numpy as np

from build_zygote import XY_UM, Z_UM, unit, mask_and_transcripts, detect_polar_body
from build_planes_all import ATLAS, SRC, side_counts, in_plane_basis
from embryo_naming import embryo_label

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SPERM = os.path.join(DATA, "zygote_sperm.json")
ZY_MAN = os.path.join(DATA, "zygote_manifest.json")
PA_NORMALS = os.path.join(DATA, "planes_all_normals.json.gz")
PA_DIR = os.path.join(DATA, "planes_all")
OUT = os.path.join(DATA, "compare_planes.json.gz")

N_POLAR = 18            # planes fanning around the polar-body axis (10° steps, 0–170°)
MIN_COUNT = 10          # min cell-body transcripts for a gene to be reported

PLANE_META = [
    {"key": "polar",      "label": "Polar-axis best",   "color": "#a855f7", "perGene": True},
    {"key": "exhaustive", "label": "Exhaustive best",   "color": "#f59e0b", "perGene": True},
    {"key": "equatorial", "label": "Equatorial",        "color": "#0ea5e9", "perGene": False},
    {"key": "sperm",      "label": "Sperm · COM · PB",  "color": "#ff2d95", "perGene": False},
]


def sideA(P, com, N):
    """Side-A count(s) of points P (K×3) for each plane normal in N (M×3), through com."""
    if len(P) == 0:
        return np.zeros(len(N), np.int64)
    return side_counts(P, com, N)


def process(eid, sperm_um, pa_iVol, normals):
    scene_p = os.path.join(ATLAS, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(SRC, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    d = json.load(gzip.open(scene_p, "rt"))
    tx = d.get("transcripts", {})
    genes = [g for g in tx if len(tx[g]["x"])]
    if not genes:
        return None
    pos, labels, voxvol, seg_of, inside = mask_and_transcripts(lab[0], tx, genes)
    if 1 not in labels:
        return None
    zs = d.get("z_scale", 7.0)
    pb_label, pb_det = detect_polar_body(d, labels, zs, inside)
    if pb_label is None:
        return None
    com = pos.mean(axis=0).astype(float)
    pb_com = pos[labels == pb_label].mean(axis=0).astype(float)
    axis = unit(pb_com - com)                                    # polar-body (animal–vegetal) axis
    pos1 = pos[labels == 1].astype(float)                        # CELL BODY only (excl pron/polar/other)
    Vtot = float(len(pos1) * voxvol)
    if Vtot <= 0:
        return None

    # ── the four plane normals ──
    n_eq = axis.copy()
    n_sd = unit(np.cross(sperm_um - com, pb_com - com)) if sperm_um is not None else None
    u, w = in_plane_basis(axis)                                  # span ⊥ axis
    th = np.arange(N_POLAR) * (np.pi / N_POLAR)                  # 0..170°
    N_pol = np.stack([np.cos(t) * u + np.sin(t) * w for t in th])   # (18,3) normals ⊥ axis → planes CONTAIN axis

    # ── per-normal cell-body volume splits (side A), computed ONCE per zygote ──
    volA_eq = float(sideA(pos1, com, n_eq[None])[0] * voxvol)
    volA_sd = float(sideA(pos1, com, n_sd[None])[0] * voxvol) if n_sd is not None else None
    volA_pol = sideA(pos1, com, N_pol).astype(float) * voxvol    # (18,)
    volA_grid = sideA(pos1, com, normals).astype(float) * voxvol  # (M,) for the exhaustive lookup

    def asym(nA, vA):                                            # |cA−cB|/(cA+cB) over a gene's split
        nB, vB = n_in - nA, Vtot - vA
        cA, cB = nA / max(vA, voxvol), nB / max(vB, voxvol)
        s = cA + cB
        return abs(cA - cB) / s if s > 0 else 0.0

    gene_out = {}
    for g in genes:
        seg = seg_of[g]
        in1 = seg == 1                                          # CELL-BODY transcripts only
        n_in = int(in1.sum())
        if n_in < MIN_COUNT:
            continue
        t = tx[g]
        P = np.stack([np.asarray(t["x"], float)[in1] * XY_UM,
                      np.asarray(t["y"], float)[in1] * XY_UM,
                      np.asarray(t["gz"], float)[in1] * Z_UM], axis=1)
        rec = {"nc": n_in}
        rec["eq"] = int(sideA(P, com, n_eq[None])[0])
        rec["sd"] = int(sideA(P, com, n_sd[None])[0]) if n_sd is not None else None
        # polar-axis best: pick the 18-plane orientation with the strongest concentration asymmetry
        aA_pol = sideA(P, com, N_pol).astype(float)
        best_k = int(max(range(N_POLAR), key=lambda k: asym(aA_pol[k], volA_pol[k])))
        rec["pb"] = {"n": [round(float(x), 5) for x in N_pol[best_k]],
                     "v": [round(volA_pol[best_k], 1), round(Vtot - volA_pol[best_k], 1)],
                     "c": [int(aA_pol[best_k]), n_in - int(aA_pol[best_k])], "ang": round(best_k * 180.0 / N_POLAR, 1)}
        # exhaustive best: reuse build_planes_all's density-best normal for this gene
        idx = pa_iVol.get(g)
        if idx is not None:
            nA = int(sideA(P, com, normals[idx][None])[0]); vA = volA_grid[idx]
            rec["ex"] = {"n": [round(float(x), 5) for x in normals[idx]],
                         "v": [round(vA, 1), round(Vtot - vA, 1)], "c": [nA, n_in - nA]}
        else:
            rec["ex"] = None
        gene_out[g] = rec

    com_plot = [round(com[0] / XY_UM, 2), round(com[1] / XY_UM, 2), round(com[2] * zs, 2)]
    ex = d.get("extents") or {}
    L_um = 0.5 * 1.3 * max((ex.get("x", [0, 1])[1] - ex.get("x", [0, 1])[0]) * XY_UM,
                           (ex.get("y", [0, 1])[1] - ex.get("y", [0, 1])[0]) * XY_UM,
                           1.0) if ex else 46.0
    return {
        "id": eid, "z_scale": zs,
        "com_um": [round(float(x), 4) for x in com], "com_plot": com_plot,
        "L_um": round(float(L_um), 1), "Vtot": round(Vtot, 1),
        "axis_um": [round(float(x), 5) for x in axis],
        "eq": {"n": [round(float(x), 5) for x in n_eq], "v": [round(volA_eq, 1), round(Vtot - volA_eq, 1)]},
        "sd": ({"n": [round(float(x), 5) for x in n_sd], "v": [round(volA_sd, 1), round(Vtot - volA_sd, 1)]}
               if n_sd is not None else None),
        "genes": gene_out,
    }


def main():
    sperm = {e["id"]: e for e in json.load(open(SPERM))["embryos"]}
    zman = {m["id"]: m for m in json.load(open(ZY_MAN))["embryos"]}
    normals = np.asarray(json.load(gzip.open(PA_NORMALS, "rt"))["normals"], float)

    embryos = []
    for eid in sorted(os.listdir(PA_DIR)):
        if not eid.endswith(".json.gz"):
            continue
        zid = eid[:-len(".json.gz")]
        # exhaustive density-best normal index per gene, from build_planes_all's per-zygote scene
        pa = json.load(gzip.open(os.path.join(PA_DIR, eid), "rt"))
        pa_iVol = {r["gene"]: r["iVol"] for r in (pa.get("analysis", {}).get("genes") or []) if r.get("iVol") is not None}
        se = sperm.get(zid, {})
        sp_plot = se.get("sperm_plot")
        sperm_um = (np.asarray(sp_plot, float) * XY_UM) if sp_plot else None
        try:
            r = process(zid, sperm_um, pa_iVol, normals)
        except Exception as e:                                  # noqa: BLE001
            print(f"  !! {zid}: {e}")
            continue
        if not r:
            print(f"  -- skipped {zid}")
            continue
        r["label"] = (zman.get(zid, {}).get("label")) or embryo_label(zid, "zygote") or zid
        r["date_short"] = zman.get(zid, {}).get("date_short", "")
        r["has_sperm"] = r["sd"] is not None
        embryos.append(r)
        print(f"  {zid}  {r['label']}  {len(r['genes'])} genes  sperm={r['has_sperm']}")

    doc = {
        "planes": PLANE_META, "unit_um_per_plot": XY_UM, "minCount": MIN_COUNT,
        "n": len(embryos), "n_sperm": sum(1 for e in embryos if e["has_sperm"]),
        "embryos": embryos,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    ngenes = len({g for e in embryos for g in e["genes"]})
    print(f"\nwrote {len(embryos)} zygotes ({doc['n_sperm']} with sperm) · {ngenes} genes ({os.path.getsize(OUT)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
