#!/usr/bin/env python3
"""
Build the "Zygote Division Planes" project data.

For every zygote (all 60 in the atlas, sperm or not) we test whether a single
division plane through the embryo→polar-body axis describes each gene's spatial
distribution.

Geometry (physical µm space; xy = px·0.15, z = frame·1.0):
  * embryo COM     = centroid of ALL segmented voxels (labels 1–5)
  * polar body COM = centroid of segment 2 (small peripheral body)
  * axis           = unit(COM → polar-body COM)
  * 17 planes, each CONTAINS the axis, rotated 10° about it (0°,10°,…,160°).
    Plane k normal nₖ = cosθ·u + sinθ·v (u,v ⊥ axis). Side of point p =
    sign((p−COM)·nₖ).

Per plane we get, per side, the embryo VOLUME (mask voxels on that side) and, per
gene, the transcript COUNT. For each gene×plane we record counts, counts
normalized by volume and by total, their side-differences, a coin-flip NULL
(n fair flips → one representative realization) and — from N_NULL flips — two
permutation p-values (volume- and count-normalized):
    T = |normalized side-difference|;   p = (1 + #{T_null ≥ T_obs}) / (1 + N_null)

Best planes:  per gene → min-p (volume p) and max |a−b|;  per embryo → the plane
minimizing the transcript-weighted mean p and the plane maximizing Σ|a−b|/ΣN.

Outputs (per embryo, gzipped): data/zygote/<id>.json.gz  (slim render scene +
the plane analysis), and data/zygote_manifest.json.
"""
import glob
import gzip
import json
import os

import numpy as np
import tifffile

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = "/Users/rishib/Desktop/MERFISH/Website2/MerfishAtlasWebsite/public/data/Zygote"
SRC = os.path.join(HERE, "..", "TranscriptomicsData", "JustTifAndCSVData", "Zygote")
OUT_DIR = os.path.join(HERE, "public", "data", "zygote")
OUT_MANIFEST = os.path.join(HERE, "public", "data", "zygote_manifest.json")

XY_UM = 0.15          # µm per pixel
Z_UM = 1.0            # µm per frame
DS_XY = 6             # mask downsample (xy) for volume/COM
DS_Z = 2              # mask downsample (z)
N_PLANES = 17
STEP_DEG = 10.0
N_NULL = 10000        # coin-flip trials for the permutation p-value
SLAB_UM = 6.0         # cross-section slab half-thickness (µm)
N_ANG = 120           # angular bins for the cross-section outline
RNG = np.random.default_rng(20260707)


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def label_mask_coords(label_path):
    """Downsampled nonzero voxel positions (µm) + labels, from the label TIFF."""
    with tifffile.TiffFile(label_path) as t:
        arr = t.series[0]                     # (Z, Y, X)
        mm = arr.asarray(out="memmap")
        sub = np.asarray(mm[::DS_Z, ::DS_XY, ::DS_XY])   # (z,y,x) downsampled
        del mm
    iz, iy, ix = np.nonzero(sub)
    labels = sub[iz, iy, ix].astype(np.int16)
    # full-res index → µm position
    pos = np.stack([ix * DS_XY * XY_UM, iy * DS_XY * XY_UM, iz * DS_Z * Z_UM], axis=1)
    voxel_vol = (DS_XY * XY_UM) ** 2 * (DS_Z * Z_UM)      # µm³ per downsampled voxel
    return pos.astype(np.float32), labels, voxel_vol


def plane_normals(u, v):
    th = np.deg2rad(np.arange(N_PLANES) * STEP_DEG)
    normals = np.cos(th)[:, None] * u[None, :] + np.sin(th)[:, None] * v[None, :]   # (K,3)
    m = (-np.sin(th))[:, None] * u[None, :] + np.cos(th)[:, None] * v[None, :]       # in-plane dir
    return th, normals, m


def cross_section_outline(pos, com, a, u, v):
    """Star-shaped outline (max radius per angular bin) of the cell cross-section
    at the COM, perpendicular to the axis a. Returns [(uu,vv)] in µm."""
    d = pos - com
    along = d @ a
    slab = pos[np.abs(along) < SLAB_UM]
    if len(slab) < 20:
        slab = pos[np.abs(along) < SLAB_UM * 3]
    if len(slab) < 20:
        return []
    dd = slab - com
    uu = dd @ u
    vv = dd @ v
    ang = np.arctan2(vv, uu)
    rad = np.hypot(uu, vv)
    bins = ((ang + np.pi) / (2 * np.pi) * N_ANG).astype(int) % N_ANG
    out = []
    for b in range(N_ANG):
        m = bins == b
        if not m.any():
            continue
        r = rad[m].max()
        th = (b + 0.5) / N_ANG * 2 * np.pi - np.pi
        out.append([float(r * np.cos(th)), float(r * np.sin(th))])
    return out


def perm_pvals(a_obs, n, V_A, V_B, null_a):
    """Two permutation p-values (volume- and count-normalized) for one gene×plane.
    null_a: array of coin-flip side-A counts (Binomial(n,0.5), reused across planes)."""
    b_obs = n - a_obs
    # observed |normalized difference|
    T_vol_obs = abs(a_obs / V_A - b_obs / V_B)
    T_cnt_obs = abs((a_obs - b_obs) / n) if n else 0.0
    nb = n - null_a
    T_vol_null = np.abs(null_a / V_A - nb / V_B)
    T_cnt_null = np.abs((null_a - nb) / n) if n else np.zeros_like(null_a, dtype=float)
    p_vol = (1 + int(np.count_nonzero(T_vol_null >= T_vol_obs - 1e-12))) / (1 + len(null_a))
    p_cnt = (1 + int(np.count_nonzero(T_cnt_null >= T_cnt_obs - 1e-12))) / (1 + len(null_a))
    return p_vol, p_cnt


def process(eid):
    scene_p = os.path.join(ATLAS, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(SRC, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    d = json.load(gzip.open(scene_p, "rt"))
    zs = d.get("z_scale", 7.0)

    pos, labels, voxvol = label_mask_coords(lab[0])          # µm
    if 2 not in labels:
        return None
    com = pos.mean(axis=0)                                    # all voxels (µm)
    pb_com = pos[labels == 2].mean(axis=0)                    # segment 2 (µm)
    a = unit(pb_com - com)
    ref = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = unit(np.cross(a, ref)); v = unit(np.cross(a, u))
    th, normals, m = plane_normals(u, v)

    # per-plane volumes (mask voxels per side)
    proj = (pos - com) @ normals.T                           # (Nvox, K)
    volA = (proj > 0).sum(axis=0).astype(float) * voxvol
    volB = (proj <= 0).sum(axis=0).astype(float) * voxvol
    volA = np.maximum(volA, voxvol); volB = np.maximum(volB, voxvol)

    # transcripts (µm) per gene, from the atlas scene
    tx = d.get("transcripts", {})
    genes = sorted(tx.keys(), key=lambda g: -d.get("gene_totals", {}).get(g, 0))
    # side of every transcript for every plane, per gene
    K = N_PLANES
    weighted_p_vol = np.zeros(K); weighted_p_cnt = np.zeros(K)
    diff_sum = np.zeros(K); total_n = 0
    gene_rows = []
    for gi, g in enumerate(genes):
        t = tx[g]
        P = np.stack([np.asarray(t["x"], float) * XY_UM,
                      np.asarray(t["y"], float) * XY_UM,
                      np.asarray(t["gz"], float) * Z_UM], axis=1)
        n = len(P)
        if n == 0:
            continue
        gproj = (P - com) @ normals.T                        # (n, K)
        aK = (gproj > 0).sum(axis=0).astype(int)             # side-A count per plane
        bK = n - aK
        null_a = RNG.binomial(n, 0.5, N_NULL)                # reused across planes
        null_a1 = int(RNG.binomial(n, 0.5))                  # one representative realization
        planes_g = []
        pv = np.empty(K); pc = np.empty(K)
        for k in range(K):
            VA, VB = volA[k], volB[k]
            a_o, b_o = int(aK[k]), int(bK[k])
            nb1 = n - null_a1
            row = {
                "a": a_o, "b": b_o,
                "aV": a_o / VA, "bV": b_o / VB,
                "aC": a_o / n, "bC": b_o / n,
                "dCount": a_o - b_o,
                "dVol": a_o / VA - b_o / VB,
                "dNorm": (a_o - b_o) / n,
                "na": null_a1, "nb": nb1,
                "naV": null_a1 / VA, "nbV": nb1 / VB,
                "ndCount": null_a1 - nb1,
                "ndVol": null_a1 / VA - nb1 / VB,
                "ndNorm": (null_a1 - nb1) / n,
            }
            p_vol, p_cnt = perm_pvals(a_o, n, VA, VB, null_a)
            row["pVol"] = p_vol; row["pCnt"] = p_cnt
            pv[k] = p_vol; pc[k] = p_cnt
            planes_g.append(row)
        # per-gene best planes
        bestP = int(np.argmin(pv))                           # min volume-normalized p
        bestD = int(np.argmax(np.abs(aK - bK)))              # max raw |a−b|
        gene_rows.append({"gene": g, "idx": gi, "total": n, "planes": planes_g,
                          "bestP": bestP, "bestDiff": bestD})
        # per-embryo aggregates
        weighted_p_vol += n * pv; weighted_p_cnt += n * pc
        diff_sum += np.abs(aK - bK)
        total_n += n

    if total_n == 0:
        return None
    weighted_p_vol /= total_n; weighted_p_cnt /= total_n
    diff_sum_norm = diff_sum / total_n
    best_p_plane = int(np.argmin(weighted_p_cnt))            # min weighted mean p (count)
    best_p_plane_vol = int(np.argmin(weighted_p_vol))
    best_diff_plane = int(np.argmax(diff_sum_norm))

    # ---- render geometry (plot space: x_px, y_px, frame·z_scale) ----
    com_plot = [com[0] / XY_UM, com[1] / XY_UM, com[2] * zs]
    pb_plot = [pb_com[0] / XY_UM, pb_com[1] / XY_UM, pb_com[2] * zs]
    # plane render length ≈ embryo radius
    ex = d["extents"]
    L_um = 0.62 * 0.5 * max(ex["x"][1] - ex["x"][0], ex["y"][1] - ex["y"][0],
                            ex["z"][1] - ex["z"][0]) * XY_UM
    planes_geo = []
    for k in range(K):
        # quad spanned by axis a and in-plane m[k], centered at COM (µm→plot)
        planes_geo.append({
            "angle": round(k * STEP_DEG, 1),
            "a_plot": [a[0] / XY_UM, a[1] / XY_UM, a[2] * zs],
            "m_plot": [m[k][0] / XY_UM, m[k][1] / XY_UM, m[k][2] * zs],
            # normal (µm) + COM (µm) let the front-end split transcripts the SAME
            # way this precompute did (side A = (p−COM)·n > 0).
            "normal_um": [round(float(x), 6) for x in normals[k]],
            "L": L_um,
            "volA": round(float(volA[k]), 1), "volB": round(float(volB[k]), 1),
            "wpVol": round(float(weighted_p_vol[k]), 5),
            "wpCnt": round(float(weighted_p_cnt[k]), 5),
            "diffSum": round(float(diff_sum_norm[k]), 5),
        })
    outline = cross_section_outline(pos, com, a, u, v)

    scene = {
        "id": eid, "z_scale": zs, "extents": ex,
        "region_meshes": d["region_meshes"], "region_defaults": d["region_defaults"],
        "mask_labels": d["mask_labels"],
        "genes": genes, "gene_totals": d.get("gene_totals", {}),
        "transcripts": {g: {"x": tx[g]["x"], "y": tx[g]["y"], "gz": tx[g]["gz"]} for g in tx},
        "analysis": {
            "com_plot": [round(c, 2) for c in com_plot],
            "com_um": [round(float(c), 4) for c in com],
            "pb_plot": [round(c, 2) for c in pb_plot],
            "axis_plot": [round(x, 5) for x in [a[0] / XY_UM, a[1] / XY_UM, a[2] * zs]],
            "planes": planes_geo,
            "best_p_plane": best_p_plane, "best_p_plane_vol": best_p_plane_vol,
            "best_diff_plane": best_diff_plane,
            "n_null": N_NULL,
            "cross_section": {
                "u_plot": [u[0] / XY_UM, u[1] / XY_UM, u[2] * zs],
                "v_plot": [v[0] / XY_UM, v[1] / XY_UM, v[2] * zs],
                "outline": [[round(p[0], 2), round(p[1], 2)] for p in outline],
            },
            "genes": gene_rows,
        },
    }
    return scene


def _json_default(o):
    """Serialize numpy scalars/arrays (from the float32 mask math) as JSON."""
    if isinstance(o, np.floating):
        return round(float(o), 7)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not serializable: {type(o)}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    ids = sorted(os.listdir(ATLAS))
    manifest = []
    for i, eid in enumerate(ids):
        try:
            scene = process(eid)
        except Exception as e:              # noqa: BLE001
            print(f"  !! {eid}: {e}")
            continue
        if not scene:
            print(f"  -- skipped {eid}")
            continue
        out = os.path.join(OUT_DIR, eid + ".json.gz")
        with gzip.open(out, "wt") as fh:
            json.dump(scene, fh, separators=(",", ":"), default=_json_default)
        date = eid[:8]
        MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        ds = f"{MON[int(date[4:6]) - 1]} {int(date[6:8])}" if date.isdigit() else ""
        import re
        m2 = re.search(r"_p(\d+)_(.+)$", eid)
        plate = f"p{m2.group(1)}" if m2 else "?"
        fovsub = m2.group(2) if m2 else eid
        manifest.append({
            "id": eid, "label": f"{plate.upper()} · {fovsub}", "date_short": ds,
            "n_genes": len(scene["genes"]),
            "n_transcripts": sum(len(t["x"]) for t in scene["transcripts"].values()),
            "size_kb": round(os.path.getsize(out) / 1024),
            "best_p": scene["analysis"]["best_p_plane"],
            "best_diff": scene["analysis"]["best_diff_plane"],
        })
        print(f"  [{i+1}/{len(ids)}] {eid}  bestP={scene['analysis']['best_p_plane']*10}° "
              f"bestDiff={scene['analysis']['best_diff_plane']*10}°  {manifest[-1]['size_kb']}KB")
    manifest.sort(key=lambda m: m["id"])
    with open(OUT_MANIFEST, "w") as fh:
        json.dump({"embryos": manifest, "n_planes": N_PLANES, "step_deg": STEP_DEG}, fh, indent=1)
    tot = sum(m["size_kb"] for m in manifest)
    print(f"\nwrote {len(manifest)} zygotes  ({tot/1024:.1f} MB)")


if __name__ == "__main__":
    main()
