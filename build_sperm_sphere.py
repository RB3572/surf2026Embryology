#!/usr/bin/env python3
"""
Build the "Sperm-Entry-Site Enrichment" project data.

For every zygote that has a LABELLED SPERM, and every gene, test whether that gene's
transcripts are ENRICHED or DEPLETED inside a sphere of radius r (µm) centred on the
sperm entry site (SES), across a sweep of radii. Everything is done in true isotropic
µm (µm = plot × 0.15; transcript z-frame = µm), and the sphere is intersected with the
cell so the cortical sperm site isn't penalised for the part of the sphere outside the
cell.

Per (zygote, gene, radius r):
    n_cell  = the gene's transcripts inside the cell (cytoplasm + pronuclei; polar body
              and background excluded), via the full-resolution label TIFF.
    n_sph   = those within r µm of the sperm site.
    V_cell  = cell voxels;  V_sph = cell voxels within r of the sperm  (clipped sphere).
    p_null  = V_sph / V_cell                                  (uniform-in-cell chance)
    fold    = (n_sph / V_sph) / (n_cell / V_cell)             (>1 enriched, <1 depleted)
    p_enr   = P(Binomial(n_cell, p_null) >= n_sph)            (upper tail — enrichment)
    p_dep   = P(Binomial(n_cell, p_null) <= n_sph)            (lower tail — depletion)
and a SPATIAL null that asks whether the sperm site is special vs OTHER cortical sites:
    draw B random centres among cell voxels at a similar distance-from-COM as the sperm
    (a cortical shell), recount n_sph for each; report the empirical enrich/deplete p and
    the null fold band (2.5–97.5 pct). This controls for "cortical" and for local density.

Reuses build_pronuclei (imported as BP) for the atlas/label readers, the downsampled
mask, and the constants. Reuses the per-zygote 3-D scenes in data/zygote/<id>.json.gz
for rendering (no meshes re-shipped here). Output:
  data/sperm_sphere.json.gz   {radii, embryos:[{…per-gene sweep…}], byGene, meta}
Run from the deploy repo root:  python3 build_sperm_sphere.py
"""
import glob
import gzip
import json
import os

import numpy as np
import tifffile
from scipy.stats import binom

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
RADII = [5, 7, 9, 11, 13, 15, 18]        # µm sweep
MIN_COUNT = 10                           # min in-cell transcripts for a gene to be reported
FOLD_THRESH = 1.5                        # enriched-list density-fold threshold
DEP_THRESH = 1.0 / FOLD_THRESH           # depleted-list threshold (~0.67)
NULL_B = 400                             # cortical random-centre draws for the spatial null
SHELL_UM = 10.0                          # cortical-shell half-width for random centres (µm)
RNG = np.random.default_rng(20260737)


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


def cell_voxels_um(sub, polar_label):
    """µm positions of the cell voxels = every nonzero label except the polar body."""
    mask = sub > 0
    if polar_label is not None:
        mask &= (sub != polar_label)
    iz, iy, ix = np.nonzero(mask)
    return np.stack([ix * DS_XY * XY_UM, iy * DS_XY * XY_UM, iz * DS_Z * Z_UM], axis=1)


def process(eid, sperm_um, polar_label):
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
    Vc = cell_voxels_um(sub, polar_label)
    if len(Vc) == 0:
        return None
    com = Vc.mean(axis=0)
    sp = np.asarray(sperm_um, float)
    r_com = float(np.linalg.norm(sp - com))
    dcell = np.linalg.norm(Vc - sp, axis=1)                 # cell voxel → sperm distances (µm)

    # cortical random-centre pool for the spatial null: cell voxels at a similar
    # distance-from-COM as the sperm (compare the sperm site to OTHER cortical sites)
    dcom = np.linalg.norm(Vc - com, axis=1)
    shell = np.abs(dcom - r_com) < SHELL_UM
    pool = Vc[shell] if shell.sum() >= 50 else Vc           # fall back to all cell voxels
    ci = RNG.integers(0, len(pool), size=NULL_B)
    C = pool[ci]                                            # (B,3) random cortical centres

    V_cell = len(Vc)
    v_sph = {r: int((dcell <= r).sum()) for r in RADII}     # clipped-to-cell sphere volumes

    # in-cell membership of every transcript: nonzero label AND not the polar body
    def in_cell(lg):
        m = lg > 0
        if polar_label is not None:
            m &= (lg != polar_label)
        return m

    gene_out = {}
    for g in genes:
        m = in_cell(labs[g])
        n_cell = int(m.sum())
        if n_cell < MIN_COUNT:
            continue
        Pg = pos[g][m]                                      # in-cell transcript µm
        dsp = np.linalg.norm(Pg - sp, axis=1)               # → sperm
        dnull = np.linalg.norm(Pg[:, None, :] - C[None, :, :], axis=2)  # (n_cell, B) → random centres
        rec = {"n": n_cell, "nsph": [], "fold": [], "pE": [], "pD": [],
               "pSE": [], "pSD": [], "nlo": [], "nhi": []}
        for r in RADII:
            vsph = v_sph[r]
            p_null = vsph / V_cell if V_cell else 0.0
            nsph = int((dsp <= r).sum())
            fold = (nsph / vsph) / (n_cell / V_cell) if (vsph and n_cell) else 0.0
            p_enr = float(binom.sf(nsph - 1, n_cell, p_null)) if p_null > 0 else 1.0
            p_dep = float(binom.cdf(nsph, n_cell, p_null)) if p_null > 0 else 1.0
            null_counts = (dnull <= r).sum(axis=0)          # (B,) gene counts in random spheres
            ge = int((null_counts >= nsph).sum()); le = int((null_counts <= nsph).sum())
            p_se = (1 + ge) / (1 + NULL_B); p_sd = (1 + le) / (1 + NULL_B)
            null_fold = (null_counts / max(vsph, 1)) / (n_cell / V_cell) if n_cell else null_counts * 0.0
            nlo, nhi = np.percentile(null_fold, [2.5, 97.5])
            rec["nsph"].append(nsph); rec["fold"].append(round(fold, 3))
            rec["pE"].append(round(p_enr, 5)); rec["pD"].append(round(p_dep, 5))
            rec["pSE"].append(round(p_se, 4)); rec["pSD"].append(round(p_sd, 4))
            rec["nlo"].append(round(float(nlo), 3)); rec["nhi"].append(round(float(nhi), 3))
        gene_out[g] = rec

    zs = d.get("z_scale", 7.0)
    return {
        "id": eid, "z_scale": zs,
        "sperm_um": [round(float(x), 2) for x in sp],
        "sperm_plot": [round(float(sp[0] / XY_UM), 2), round(float(sp[1] / XY_UM), 2), round(float(sp[2] * zs), 2)],
        "com_plot": [round(float(com[0] / XY_UM), 2), round(float(com[1] / XY_UM), 2), round(float(com[2] * zs), 2)],
        "r_com_um": round(r_com, 2), "V_cell": V_cell,
        "v_sph": {str(r): v_sph[r] for r in RADII},
        "p_null": {str(r): round(v_sph[r] / V_cell, 5) for r in RADII},
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
        sperm_um = np.asarray(sp_plot, float) * XY_UM             # plot → isotropic µm
        pol = (assign.get(eid, {}).get("polar") or {}).get("label")
        pol = int(pol) if pol is not None else None
        try:
            r = process(eid, sperm_um, pol)
        except Exception as e:                                    # noqa: BLE001
            print(f"  !! {eid}: {e}")
            continue
        if not r:
            print(f"  -- skipped {eid}")
            continue
        r["label"] = (zman.get(eid, {}).get("label")) or embryo_label(eid, "zygote") or eid
        r["date_short"] = zman.get(eid, {}).get("date_short", "")
        embryos.append(r)
        ng = len(r["genes"])
        print(f"  {eid}  {r['label']}  {ng} genes  r(sperm→COM)={r['r_com_um']}µm  Vcell={r['V_cell']}")

    # cross-gene aggregation: per gene, how many zygotes it is enriched / depleted in at each radius
    by_gene = {}
    for emb in embryos:
        for g, rec in emb["genes"].items():
            bg = by_gene.setdefault(g, {"nz": 0, "enr": [0] * len(RADII), "dep": [0] * len(RADII)})
            bg["nz"] += 1
            for ri in range(len(RADII)):
                if rec["fold"][ri] >= FOLD_THRESH and rec["pSE"][ri] <= 0.05:
                    bg["enr"][ri] += 1
                if rec["fold"][ri] <= DEP_THRESH and rec["pSD"][ri] <= 0.05:
                    bg["dep"][ri] += 1

    doc = {
        "radii": RADII, "embryos": embryos, "byGene": by_gene,
        "meta": {"nZygotes": len(embryos), "minCount": MIN_COUNT, "foldThresh": FOLD_THRESH,
                 "depThresh": round(DEP_THRESH, 3), "nullB": NULL_B, "shellUm": SHELL_UM,
                 "unit_um_per_plot": XY_UM, "defaultRadiusIdx": RADII.index(9)},
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(f"\nwrote {len(embryos)} sperm zygotes · {len(by_gene)} genes "
          f"({os.path.getsize(OUT)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
