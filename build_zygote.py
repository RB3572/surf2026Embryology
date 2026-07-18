#!/usr/bin/env python3
"""
Build the "Zygote Division Planes" project data.

For every zygote (all 60 in the atlas, sperm or not) we test whether a single
division plane through the embryo→polar-body axis describes each gene's spatial
distribution.

Geometry (physical µm space; xy = px·0.15, z = frame·1.0):
  * embryo COM     = centroid of ALL segmented voxels
  * polar body COM = centroid of the most peripheral segmented object that lies
                     outside the filled segment-1 cell body
  * axis           = unit(COM → polar-body COM)
  * 18 planes, each CONTAINS the axis, rotated 10° about it (0°,10°,…,170°).
    Plane k normal nₖ = cosθ·u + sinθ·v (u,v ⊥ axis). Side of point p =
    sign((p−COM)·nₖ).

Only SEGMENT 1 (the cell body) counts: transcripts inside the polar body or the
pronuclei are excluded from every count, and the per-side
VOLUME used for normalization is the volume of segment 1 on that side only.

Per plane we get, per side, the segment-1 VOLUME (label-1 voxels on that side) and,
per gene, the segment-1 transcript COUNT. For each gene×plane we record counts,
counts normalized by volume and by total, their side-differences, a coin-flip NULL
(n fair flips → one representative realization) and — from N_NULL flips — two
permutation p-values (volume- and count-normalized):
    T = |normalized side-difference|;   p = 2·(1 + #{T_null ≥ T_obs}) / (1 + N_null)

Best planes (per embryo): FOUR planes = {min transcript-weighted-mean p} and
{max Σ side-difference / ΣN}, each computed under BOTH normalizations (volume /
count) → best_planes = {pVol, pCnt, diffVol, diffCnt}. Per gene we also store the
min-p and max-|diff| plane under each normalization.

Outputs (per embryo, gzipped): data/zygote/<id>.json.gz  (slim render scene +
the plane analysis), and data/zygote_manifest.json.
"""
import glob
import gzip
import json
import os
from collections import Counter

import numpy as np
import tifffile
from scipy.ndimage import binary_fill_holes
from scipy.spatial import ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = "/Users/rishib/Desktop/MERFISH/Website2/MerfishAtlasWebsite/public/data/Zygote"
SRC = os.path.join(HERE, "..", "TranscriptomicsData", "JustTifAndCSVData", "Zygote")
OUT_DIR = os.path.join(HERE, "data", "zygote")
OUT_MANIFEST = os.path.join(HERE, "data", "zygote_manifest.json")

XY_UM = 0.15          # µm per pixel
Z_UM = 1.0            # µm per frame
DS_XY = 6             # mask downsample (xy) for volume/COM
DS_Z = 2              # mask downsample (z)
N_PLANES = 18          # 0°–170° in 10° steps (tiles the half-circle)
STEP_DEG = 10.0
N_NULL = 10000        # coin-flip trials for the permutation p-value
SLAB_UM = 6.0         # cross-section slab half-thickness (µm)
N_ANG = 120           # angular bins for the cross-section outline
RNG = np.random.default_rng(20260707)

# Independent cell-mask and mesh checks determine whether a candidate is external.
# Distances are in physical micrometres.
PB_OUTSIDE_TOL_UM = 0.75
PB_MIN_OUTSIDE_FRACTION = 0.02
PB_MIN_MAX_DISTANCE_UM = 1.5
PB_MAX_INSIDE_FRACTION = 0.5
PB_MIN_RADIAL_RATIO = 0.6
PB_SURFACE_CONE_DEG = 20.0


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def mask_and_transcripts(label_path, tx, genes):
    """From the label TIFF, return (a) downsampled nonzero voxel positions (µm),
    their segment labels, and the per-voxel volume, and (b) `seg_of`: the segment
    label under every transcript (full-res nearest-voxel lookup), per gene.
    Segment 1 is the cell body. Other label identities vary by embryo and must
    be inferred from geometry rather than label number."""
    with tifffile.TiffFile(label_path) as t:
        mm = t.series[0].asarray(out="memmap")            # (Z, Y, X)
        Zn, Yn, Xn = mm.shape
        sub = np.asarray(mm[::DS_Z, ::DS_XY, ::DS_XY])     # (z,y,x) downsampled
        # segment label under each transcript (one full-res fancy-index for all)
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
    cell_envelope = binary_fill_holes(sub == 1)
    inside_fractions = {}
    for label in np.unique(labels):
        if label == 1:
            continue
        segment = sub == label
        inside_fractions[int(label)] = float((segment & cell_envelope).sum() / segment.sum())
    # full-res index → µm position
    pos = np.stack([ix * DS_XY * XY_UM, iy * DS_XY * XY_UM, iz * DS_Z * Z_UM], axis=1)
    voxel_vol = (DS_XY * XY_UM) ** 2 * (DS_Z * Z_UM)      # µm³ per downsampled voxel
    return pos.astype(np.float32), labels, voxel_vol, seg_of, inside_fractions


def mesh_vertices_um(mesh, z_scale):
    """Convert an atlas render mesh from plot coordinates to physical µm."""
    vertices = np.asarray(mesh["verts"], dtype=float).reshape(-1, 3).copy()
    vertices[:, 0] *= XY_UM
    vertices[:, 1] *= XY_UM
    vertices[:, 2] = vertices[:, 2] / z_scale * Z_UM
    return vertices


def detect_polar_body(scene, labels, z_scale, inside_fractions):
    """Return the most peripheral segment demonstrably outside the cell body.

    The filled 3-D segment-1 mask distinguishes enclosed pronuclear cavities from
    background-connected objects. A convex-hull and directional-radius check then
    guards against small segmentation leaks and identifies the most peripheral
    valid object. No label number is assumed.
    Returns ``(label, diagnostics)`` or ``(None, diagnostics)`` when the scene has
    no defensible external segment; callers must fail closed in that case.
    """
    region_meshes = scene.get("region_meshes", {})
    body_mesh = region_meshes.get("1")
    if not body_mesh:
        return None, {"reason": "missing cell-body mesh", "candidates": []}

    body_vertices = mesh_vertices_um(body_mesh, z_scale)
    if len(body_vertices) < 4:
        return None, {"reason": "insufficient cell-body mesh", "candidates": []}
    try:
        hull = ConvexHull(body_vertices, qhull_options="QJ")
    except Exception as exc:  # noqa: BLE001
        return None, {"reason": f"cell-body hull failed: {exc}", "candidates": []}

    plane_normals = hull.equations[:, :3]
    plane_offsets = hull.equations[:, 3]
    body_center = body_vertices.mean(axis=0)
    body_delta = body_vertices - body_center
    body_radius = np.linalg.norm(body_delta, axis=1)
    body_directions = body_delta / np.maximum(body_radius[:, None], 1e-9)
    label_counts = Counter(int(label) for label in labels)
    candidates = []
    for label in sorted(label_counts):
        if label == 1:
            continue
        mesh = region_meshes.get(str(label))
        if not mesh:
            continue
        vertices = mesh_vertices_um(mesh, z_scale)
        # Positive signed distance means a point is beyond at least one hull face.
        signed = np.max(plane_normals @ vertices.T + plane_offsets[:, None], axis=0)
        outside_fraction = float(np.mean(signed > PB_OUTSIDE_TOL_UM))
        max_outside_um = float(np.max(signed))
        hull_external = (
            outside_fraction >= PB_MIN_OUTSIDE_FRACTION
            or max_outside_um >= PB_MIN_MAX_DISTANCE_UM
        )
        candidate_delta = vertices.mean(axis=0) - body_center
        candidate_radius = float(np.linalg.norm(candidate_delta))
        direction = candidate_delta / max(candidate_radius, 1e-9)
        directional = body_radius[
            body_directions @ direction > np.cos(np.deg2rad(PB_SURFACE_CONE_DEG))
        ]
        surface_radius = float(
            np.percentile(directional, 90) if len(directional) >= 10 else np.percentile(body_radius, 90)
        )
        radial_ratio = candidate_radius / max(surface_radius, 1e-9)
        inside_fraction = float(inside_fractions.get(int(label), 1.0))
        is_external = (
            inside_fraction < PB_MAX_INSIDE_FRACTION
            and (hull_external or radial_ratio >= PB_MIN_RADIAL_RATIO)
        )
        candidates.append({
            "label": int(label),
            "voxel_count": int(label_counts[label]),
            "inside_cell_fraction": round(inside_fraction, 4),
            "radial_ratio": round(radial_ratio, 4),
            "outside_fraction": round(outside_fraction, 4),
            "max_outside_um": round(max_outside_um, 3),
            "external": bool(is_external),
        })

    external = [candidate for candidate in candidates if candidate["external"]]
    if not external:
        return None, {"reason": "no segment passes outside-cell geometry checks", "candidates": candidates}
    selected = max(
        external,
        key=lambda candidate: (
            candidate["radial_ratio"],
            candidate["outside_fraction"],
            candidate["max_outside_um"],
            candidate["voxel_count"],
        ),
    )
    return int(selected["label"]), {"reason": "most peripheral external segment", "candidates": candidates}


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
    # two-sided: p = 2 * (1 + #{T_null >= T_obs}) / (1 + N_null)
    p_vol = 2 * (1 + int(np.count_nonzero(T_vol_null >= T_vol_obs - 1e-12))) / (1 + len(null_a))
    p_cnt = 2 * (1 + int(np.count_nonzero(T_cnt_null >= T_cnt_obs - 1e-12))) / (1 + len(null_a))
    return p_vol, p_cnt


CIRC_W = 0.95          # balloon inflation weight: 0.95 = 95% toward the average-radius sphere (smooth, round)
CIRC_NBIN = 22         # (theta) surface-grid resolution for the radial surface function
CIRC_TOL = 0.05        # min sphericity: no direction may bulge more than ±5% off the average-radius
                       # sphere. Any part that still deviates past this after the 95% pull is clamped
                       # (circularized further) so residual lobes/spikes can't survive.
RNG_C = np.random.default_rng(20260715)   # separate RNG so circ nulls don't shift the real p-values


def _surface_grid(pos1, center, nbin=CIRC_NBIN):
    """Max-radius-per-direction surface R_surf(û) of segment 1, as a (theta, phi) grid
    (empties filled with the mean, lightly smoothed). Returns (grid, R_avg)."""
    from scipy.ndimage import uniform_filter
    d = pos1 - center
    r = np.linalg.norm(d, axis=1)
    u = d / np.maximum(r, 1e-9)[:, None]
    theta = np.arccos(np.clip(u[:, 2], -1, 1))
    phi = np.arctan2(u[:, 1], u[:, 0])
    ti = np.clip((theta / np.pi * nbin).astype(int), 0, nbin - 1)
    pj = ((phi + np.pi) / (2 * np.pi) * (2 * nbin)).astype(int) % (2 * nbin)
    grid = np.zeros((nbin, 2 * nbin))
    np.maximum.at(grid, (ti, pj), r)
    filled = grid[grid > 0]
    grid[grid == 0] = float(np.median(filled)) if filled.size else 1.0
    grid = uniform_filter(grid, size=3, mode="wrap")
    # median (not mean) target radius — robust to thin spikes that would inflate the sphere
    return grid, float(np.median(grid))


def balloon(pos1, w=CIRC_W):
    """Blow-up-the-balloon transform for segment 1: inflate toward the average-radius
    sphere by weight w, keeping interior points proportionally inside.
      ρ = r / R_surf(û);  R_new(û) = (1−w)·R_surf(û) + w·R_avg;  r' = ρ·R_new(û).
    Returns a function fn(P µm) -> P' µm (P are absolute µm positions), plus (center, R_avg)."""
    C = pos1.mean(axis=0)
    grid, R_avg = _surface_grid(pos1, C)
    nbin = grid.shape[0]

    def R_surf(uu):
        th = np.arccos(np.clip(uu[:, 2], -1, 1))
        ph = np.arctan2(uu[:, 1], uu[:, 0])
        ti = np.clip((th / np.pi * nbin).astype(int), 0, nbin - 1)
        pj = ((ph + np.pi) / (2 * np.pi) * (2 * nbin)).astype(int) % (2 * nbin)
        return grid[ti, pj]

    def fn(P):
        d = np.asarray(P, float) - C
        r = np.linalg.norm(d, axis=1)
        u = d / np.maximum(r, 1e-9)[:, None]
        Rs = R_surf(u)
        # clamp ρ ≤ 1: a point can't be past the surface (mis-binned spikes would otherwise
        # overshoot outward and stay lumpy) — this pulls protrusions onto the sphere.
        rho = np.minimum(r / np.maximum(Rs, 1e-9), 1.0)
        Rnew = (1 - w) * Rs + w * R_avg
        # enforce a minimum sphericity: clamp any direction that still deviates past ±CIRC_TOL
        # from the average-radius sphere, so a residual lobe/spike gets circularized further.
        Rnew = np.clip(Rnew, R_avg * (1 - CIRC_TOL), R_avg * (1 + CIRC_TOL))
        return C + u * (rho * Rnew)[:, None]

    return fn, C, R_avg


def analyze(pos1, com, pb_com, genes, tx1, voxvol, ex, zs, rng):
    """The full division-plane analysis for one geometry. `pos1` = segment-1 voxels (µm),
    `com`/`pb_com` = cell / polar-body centroids (µm), `tx1` = {gene: seg-1 transcript µm
    positions}. Returns the `analysis` dict (or None if no segment-1 transcripts)."""
    a = unit(pb_com - com)
    ref = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = unit(np.cross(a, ref)); v = unit(np.cross(a, u))
    th, normals, m = plane_normals(u, v)
    K = N_PLANES

    proj = (pos1 - com) @ normals.T
    volA = (proj > 0).sum(axis=0).astype(float) * voxvol
    volB = (proj <= 0).sum(axis=0).astype(float) * voxvol
    volA = np.maximum(volA, voxvol); volB = np.maximum(volB, voxvol)

    weighted_p_vol = np.zeros(K); weighted_p_cnt = np.zeros(K)
    diff_cnt_sum = np.zeros(K); diff_vol_sum = np.zeros(K)
    total_n = 0; gene_rows = []
    for gi, g in enumerate(genes):
        P = tx1[g]
        n = len(P)
        if n == 0:
            continue
        gproj = (P - com) @ normals.T
        aK = (gproj > 0).sum(axis=0).astype(int)
        bK = n - aK
        null_a = rng.binomial(n, 0.5, N_NULL)
        null_a1 = int(rng.binomial(n, 0.5))
        planes_g = []; pv = np.empty(K); pc = np.empty(K)
        for k in range(K):
            VA, VB = volA[k], volB[k]
            a_o, b_o = int(aK[k]), int(bK[k])
            nb1 = n - null_a1
            row = {
                "a": a_o, "b": b_o, "aV": a_o / VA, "bV": b_o / VB, "aC": a_o / n, "bC": b_o / n,
                "dCount": a_o - b_o, "dVol": a_o / VA - b_o / VB, "dNorm": (a_o - b_o) / n,
                "na": null_a1, "nb": nb1, "naV": null_a1 / VA, "nbV": nb1 / VB,
                "ndCount": null_a1 - nb1, "ndVol": null_a1 / VA - nb1 / VB, "ndNorm": (null_a1 - nb1) / n,
            }
            p_vol, p_cnt = perm_pvals(a_o, n, VA, VB, null_a)
            row["pVol"] = p_vol; row["pCnt"] = p_cnt
            pv[k] = p_vol; pc[k] = p_cnt
            planes_g.append(row)
        dVolK = aK / volA - bK / volB
        gene_rows.append({"gene": g, "idx": gi, "total": n, "planes": planes_g,
                          "bestP_vol": int(np.argmin(pv)), "bestP_cnt": int(np.argmin(pc)),
                          "bestDiff_vol": int(np.argmax(np.abs(dVolK))),
                          "bestDiff_cnt": int(np.argmax(np.abs(aK - bK)))})
        weighted_p_vol += n * pv; weighted_p_cnt += n * pc
        diff_cnt_sum += np.abs(aK - bK); diff_vol_sum += n * np.abs(dVolK)
        total_n += n

    if total_n == 0:
        return None
    weighted_p_vol /= total_n; weighted_p_cnt /= total_n
    diff_cnt_norm = diff_cnt_sum / total_n; diff_vol_norm = diff_vol_sum / total_n
    best_planes = {"pVol": int(np.argmin(weighted_p_vol)), "pCnt": int(np.argmin(weighted_p_cnt)),
                   "diffVol": int(np.argmax(diff_vol_norm)), "diffCnt": int(np.argmax(diff_cnt_norm))}

    com_plot = [com[0] / XY_UM, com[1] / XY_UM, com[2] * zs]
    pb_plot = [pb_com[0] / XY_UM, pb_com[1] / XY_UM, pb_com[2] * zs]
    L_um = 0.62 * 0.5 * max(ex["x"][1] - ex["x"][0], ex["y"][1] - ex["y"][0],
                            ex["z"][1] - ex["z"][0]) * XY_UM
    planes_geo = []
    for k in range(K):
        planes_geo.append({
            "angle": round(k * STEP_DEG, 1),
            "a_plot": [a[0] / XY_UM, a[1] / XY_UM, a[2] * zs],
            "m_plot": [m[k][0] / XY_UM, m[k][1] / XY_UM, m[k][2] * zs],
            "normal_um": [round(float(x), 6) for x in normals[k]], "L": L_um,
            "volA": round(float(volA[k]), 1), "volB": round(float(volB[k]), 1),
            "wpVol": round(float(weighted_p_vol[k]), 5), "wpCnt": round(float(weighted_p_cnt[k]), 5),
            "dmVol": round(float(diff_vol_norm[k]), 7), "dmCnt": round(float(diff_cnt_norm[k]), 5),
        })
    outline = cross_section_outline(pos1, com, a, u, v)
    return {
        "com_plot": [round(c, 2) for c in com_plot], "com_um": [round(float(c), 4) for c in com],
        "pb_plot": [round(c, 2) for c in pb_plot],
        "axis_plot": [round(x, 5) for x in [a[0] / XY_UM, a[1] / XY_UM, a[2] * zs]],
        "planes": planes_geo, "best_planes": best_planes, "n_null": N_NULL,
        "cross_section": {"u_plot": [u[0] / XY_UM, u[1] / XY_UM, u[2] * zs],
                          "v_plot": [v[0] / XY_UM, v[1] / XY_UM, v[2] * zs],
                          "outline": [[round(p[0], 2), round(p[1], 2)] for p in outline]},
        "genes": gene_rows,
    }


def _circ_transcripts(tx, seg_of, fn):
    """Circularized display transcripts: seg-1 molecules moved by the balloon fn (µm),
    everything else unchanged. Kept in the same (x_px, y_px, gz_frame) units + s1 flag."""
    out = {}
    for g, t in tx.items():
        x = np.asarray(t["x"], float); y = np.asarray(t["y"], float); gz = np.asarray(t["gz"], float)
        s1 = seg_of[g] == 1
        if s1.any():
            P = np.stack([x[s1] * XY_UM, y[s1] * XY_UM, gz[s1] * Z_UM], axis=1)
            Pc = fn(P)
            x = x.copy(); y = y.copy(); gz = gz.copy()
            x[s1] = Pc[:, 0] / XY_UM; y[s1] = Pc[:, 1] / XY_UM; gz[s1] = Pc[:, 2] / Z_UM
        out[g] = {"x": [round(float(a), 2) for a in x], "y": [round(float(a), 2) for a in y],
                  "gz": [round(float(a), 2) for a in gz], "s1": s1.astype(np.uint8).tolist()}
    return out


def _circ_mesh(mesh, fn, zs):
    """Balloon-transform a region mesh (plot-space verts x_px,y_px,frame·zs)."""
    verts = np.asarray(mesh["verts"], float).reshape(-1, 3)
    P = np.stack([verts[:, 0] * XY_UM, verts[:, 1] * XY_UM, (verts[:, 2] / zs) * Z_UM], axis=1)
    Pc = fn(P)
    vc = np.stack([Pc[:, 0] / XY_UM, Pc[:, 1] / XY_UM, (Pc[:, 2] / Z_UM) * zs], axis=1)
    return {"verts": [round(float(a), 1) for a in vc.ravel()], "faces": mesh["faces"]}


def _circularize_outline(outline, tol=CIRC_TOL):
    """Keep a circularized cross-section within a tolerance of its sphere WITHOUT flattening
    it into a perfect circle: clamp every radius into a ±tol band around the median (fills
    inward notches like the polar-body cleft, caps outward spikes) so no inward spike can
    survive, then apply only a LIGHT 3-bin circular smooth to soften the clamp corners —
    leaving the cell's natural bump texture intact. Angles preserved; returns [[x,y]] µm."""
    if len(outline) < 8:
        return outline
    P = np.asarray(outline, float)
    ang = np.arctan2(P[:, 1], P[:, 0])
    r = np.hypot(P[:, 0], P[:, 1])
    order = np.argsort(ang)
    ang, r = ang[order], r[order]
    Rmed = float(np.median(r))
    r = np.clip(r, Rmed * (1.0 - tol), Rmed * (1.0 + tol))   # within tolerance of the sphere; no inward spike
    k = 1                                                     # light 3-bin smooth (wrap): soften corners, keep texture
    rp = np.concatenate([r[-k:], r, r[:k]])
    r = np.convolve(rp, np.ones(2 * k + 1) / (2 * k + 1), mode="same")[k:-k]
    return [[round(float(rr * np.cos(a)), 2), round(float(rr * np.sin(a)), 2)] for a, rr in zip(ang, r)]


def process(eid):
    scene_p = os.path.join(ATLAS, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(SRC, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    d = json.load(gzip.open(scene_p, "rt"))
    zs = d.get("z_scale", 7.0)

    # transcripts (µm) per gene, from the atlas scene
    tx = d.get("transcripts", {})
    genes = sorted(tx.keys(), key=lambda g: -d.get("gene_totals", {}).get(g, 0))

    pos, labels, voxvol, seg_of, inside_fractions = mask_and_transcripts(lab[0], tx, genes)   # µm
    if 1 not in labels:
        return None
    polar_body_label, polar_body_detection = detect_polar_body(d, labels, zs, inside_fractions)
    if polar_body_label is None:
        print(f"  -- {eid}: {polar_body_detection['reason']}")
        return None
    com = pos.mean(axis=0)                                    # all voxels (µm)
    pb_com = pos[labels == polar_body_label].mean(axis=0)     # external segment (µm)
    a = unit(pb_com - com)
    ref = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = unit(np.cross(a, ref)); v = unit(np.cross(a, u))
    th, normals, m = plane_normals(u, v)

    # Per-plane volumes use SEGMENT 1 only. The geometry-selected polar body and
    # all internal segments are excluded from volume normalization.
    pos1 = pos[labels == 1]
    proj = (pos1 - com) @ normals.T                          # (Nvox1, K)
    volA = (proj > 0).sum(axis=0).astype(float) * voxvol
    volB = (proj <= 0).sum(axis=0).astype(float) * voxvol
    volA = np.maximum(volA, voxvol); volB = np.maximum(volB, voxvol)

    # side of every segment-1 transcript for every plane, per gene
    K = N_PLANES
    weighted_p_vol = np.zeros(K); weighted_p_cnt = np.zeros(K)
    diff_cnt_sum = np.zeros(K)          # Σ|a−b|                (count difference)
    diff_vol_sum = np.zeros(K)          # Σ n·|a/V_A − b/V_B|   (volume difference)
    total_n = 0
    gene_rows = []
    for gi, g in enumerate(genes):
        t = tx[g]
        in1 = seg_of[g] == 1                                 # count only segment-1 transcripts
        P = np.stack([np.asarray(t["x"], float)[in1] * XY_UM,
                      np.asarray(t["y"], float)[in1] * XY_UM,
                      np.asarray(t["gz"], float)[in1] * Z_UM], axis=1)
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
        # per-gene best planes (min p / max |diff|, each by both normalizations)
        dVolK = aK / volA - bK / volB
        gene_rows.append({"gene": g, "idx": gi, "total": n, "planes": planes_g,
                          "bestP_vol": int(np.argmin(pv)), "bestP_cnt": int(np.argmin(pc)),
                          "bestDiff_vol": int(np.argmax(np.abs(dVolK))),
                          "bestDiff_cnt": int(np.argmax(np.abs(aK - bK)))})
        # per-embryo aggregates
        weighted_p_vol += n * pv; weighted_p_cnt += n * pc
        diff_cnt_sum += np.abs(aK - bK)
        diff_vol_sum += n * np.abs(dVolK)
        total_n += n

    if total_n == 0:
        return None
    weighted_p_vol /= total_n; weighted_p_cnt /= total_n
    diff_cnt_norm = diff_cnt_sum / total_n      # count-normalized difference metric
    diff_vol_norm = diff_vol_sum / total_n      # volume-normalized difference metric
    # 4 per-embryo best planes: {min weighted-p, max Σ|diff|} × {volume, count}
    best_planes = {
        "pVol": int(np.argmin(weighted_p_vol)),
        "pCnt": int(np.argmin(weighted_p_cnt)),
        "diffVol": int(np.argmax(diff_vol_norm)),
        "diffCnt": int(np.argmax(diff_cnt_norm)),
    }

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
            "dmVol": round(float(diff_vol_norm[k]), 7),
            "dmCnt": round(float(diff_cnt_norm[k]), 5),
        })
    # Cross-section of the CYTOPLASM ONLY (segment 1). Everything else — polar body,
    # pronuclei, and any extra peripheral bodies — is excluded,
    # so the outline is the embryo cell boundary with no external blobs bulging it out.
    outline = cross_section_outline(pos[labels == 1], com, a, u, v)

    scene = {
        "id": eid, "z_scale": zs, "extents": ex,
        "region_meshes": d["region_meshes"], "region_defaults": d["region_defaults"],
        "mask_labels": d["mask_labels"],
        "genes": genes, "gene_totals": d.get("gene_totals", {}),
        # s1[i] = 1 if transcript i is in segment 1 (counted), else 0 (rendered green)
        "transcripts": {g: {"x": tx[g]["x"], "y": tx[g]["y"], "gz": tx[g]["gz"],
                            "s1": (seg_of[g] == 1).astype(np.uint8).tolist()} for g in tx},
        "analysis": {
            "com_plot": [round(c, 2) for c in com_plot],
            "com_um": [round(float(c), 4) for c in com],
            "pb_plot": [round(c, 2) for c in pb_plot],
            "polar_body_label": int(polar_body_label),
            "polar_body_detection": polar_body_detection,
            "axis_plot": [round(x, 5) for x in [a[0] / XY_UM, a[1] / XY_UM, a[2] * zs]],
            "planes": planes_geo,
            "best_planes": best_planes,
            "n_null": N_NULL,
            "cross_section": {
                "u_plot": [u[0] / XY_UM, u[1] / XY_UM, u[2] * zs],
                "v_plot": [v[0] / XY_UM, v[1] / XY_UM, v[2] * zs],
                "outline": [[round(p[0], 2), round(p[1], 2)] for p in outline],
            },
            "genes": gene_rows,
        },
    }

    # ---- circularized ("blow up the balloon") variant: segment 1 only ----
    # Inflate seg-1 to its average-radius sphere (90%), carrying seg-1 transcripts inside,
    # and recompute the WHOLE analysis on that geometry. A separate RNG keeps the real
    # p-values above byte-identical. Stored under scene["circ"]; the front-end toggles it.
    tx1 = {}
    for g in genes:
        t = tx[g]; in1 = seg_of[g] == 1
        tx1[g] = np.stack([np.asarray(t["x"], float)[in1] * XY_UM,
                           np.asarray(t["y"], float)[in1] * XY_UM,
                           np.asarray(t["gz"], float)[in1] * Z_UM], axis=1)
    fn, _C, _Ravg = balloon(pos1)
    pos1_c = fn(pos1)
    pos_c = pos.copy(); pos_c[labels == 1] = pos1_c
    com_c = pos_c.mean(axis=0)
    tx1_c = {g: (fn(P) if len(P) else P) for g, P in tx1.items()}
    analysis_c = analyze(pos1_c, com_c, pb_com, genes, tx1_c, voxvol, ex, zs, RNG_C)
    if analysis_c is not None:
        # min-sphericity pass: the circularized cross-section can still carry an inward
        # notch (the polar-body cleft measured from the whole-embryo centroid); clamp +
        # smooth it into the sphericity band so no spike survives.
        cs = analysis_c.get("cross_section")
        if cs and cs.get("outline"):
            cs["outline"] = _circularize_outline(cs["outline"])
        m1 = d["region_meshes"].get("1")
        scene["circ"] = {"analysis": analysis_c,
                         "transcripts": _circ_transcripts(tx, seg_of, fn),
                         "mesh1": _circ_mesh(m1, fn, zs) if m1 else None,
                         "R_avg_um": round(float(_Ravg), 2)}
    else:
        scene["circ"] = None

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


BEST_KEYS = ["pVol", "pCnt", "diffVol", "diffCnt"]


def agg_entry(eid, label, A):
    """Slim per-embryo record for the cross-embryo bottom-drawer visuals: the
    cross-section outline, the 4 best-plane indices, and — per gene — its total n
    and side-A count `a` at each of the 4 best planes (b = n − a). `A` = an analysis
    dict (scene['analysis'] for real, scene['circ']['analysis'] for circularized)."""
    bp = A["best_planes"]
    plns = A["planes"]
    best = [bp[k] for k in BEST_KEYS]
    g = {}
    for row in A["genes"]:
        planes = row["planes"]
        g[row["gene"]] = [row["total"]] + [planes[bi]["a"] for bi in best]
    # per-best-plane significance (transcript-weighted-mean permutation p) for colouring the
    # cross-section outlines by how significantly that embryo's transcriptome splits.
    _sigk = {"pVol": "wpVol", "pCnt": "wpCnt", "diffVol": "wpVol", "diffCnt": "wpCnt"}
    sig = {k: round(float(plns[bp[k]][_sigk[k]]), 5) for k in BEST_KEYS}
    return {"id": eid, "label": label, "outline": A["cross_section"]["outline"],
            "best": best, "sig": sig, "g": g}


def write_cross_aggregate(entries, path):
    """One file for the bottom drawer. The dataset spans multiple gene panels, so
    NO gene is in all embryos; the alignment-gene dropdown is therefore the UNION
    of all genes (alphabetical) with each gene's embryo coverage. Selecting a gene
    shows only the embryos that contain it. Default = widest-coverage gene (ties
    broken by summed count)."""
    if not entries:
        return
    cov = Counter()
    sums = Counter()
    for e in entries:
        for gn, arr in e["g"].items():
            cov[gn] += 1
            sums[gn] += arr[0]
    genes_all = sorted(cov.keys(), key=str.lower)
    gene_cov = {gn: cov[gn] for gn in genes_all}
    default_align = max(genes_all, key=lambda gn: (cov[gn], sums[gn])) if genes_all else None
    agg = {"step_deg": STEP_DEG, "best_keys": BEST_KEYS, "n_embryos": len(entries),
           "genes_all": genes_all, "gene_cov": gene_cov,
           "default_align_gene": default_align, "embryos": entries}
    with gzip.open(path, "wt") as fh:
        json.dump(agg, fh, separators=(",", ":"), default=_json_default)
    dcov = cov[default_align] if default_align else 0
    print(f"  cross-aggregate: {len(entries)} embryos, {len(genes_all)} union genes, "
          f"default align gene = {default_align} (cov {dcov}/{len(entries)})  "
          f"({os.path.getsize(path)/1024:.0f} KB)")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    ids = sorted(os.listdir(ATLAS))
    manifest = []
    agg_entries = []
    agg_circ = []
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
        label = f"{plate.upper()} · {fovsub}"
        bp = scene["analysis"]["best_planes"]
        manifest.append({
            "id": eid, "label": label, "date_short": ds,
            "n_genes": len(scene["genes"]),
            "n_transcripts": sum(len(t["x"]) for t in scene["transcripts"].values()),
            "size_kb": round(os.path.getsize(out) / 1024),
            "polar_body_label": scene["analysis"]["polar_body_label"],
            "best_planes": bp,
        })
        agg_entries.append(agg_entry(eid, label, scene["analysis"]))
        if scene.get("circ"):
            agg_circ.append(agg_entry(eid, label, scene["circ"]["analysis"]))
        print(f"  [{i+1}/{len(ids)}] {eid}  pVol={bp['pVol']*10}° pCnt={bp['pCnt']*10}° "
              f"diffVol={bp['diffVol']*10}° diffCnt={bp['diffCnt']*10}°  {manifest[-1]['size_kb']}KB")
    manifest.sort(key=lambda m: m["id"])
    agg_entries.sort(key=lambda e: e["id"])
    with open(OUT_MANIFEST, "w") as fh:
        json.dump({"embryos": manifest, "n_planes": N_PLANES, "step_deg": STEP_DEG}, fh, indent=1)
    write_cross_aggregate(agg_entries, os.path.join(HERE, "data", "zygote_cross.json.gz"))
    if agg_circ:
        agg_circ.sort(key=lambda e: e["id"])
        write_cross_aggregate(agg_circ, os.path.join(HERE, "data", "zygote_cross_circ.json.gz"))
    tot = sum(m["size_kb"] for m in manifest)
    print(f"\nwrote {len(manifest)} zygotes  ({tot/1024:.1f} MB)")


if __name__ == "__main__":
    main()
