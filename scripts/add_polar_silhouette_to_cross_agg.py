#!/usr/bin/env python3
"""
Enrich the Zygote Division Planes cross-embryo aggregates with each embryo's
POLAR-BODY SILHOUETTE, in the SAME (u,v) frame as the stored cell `outline`
(the plane ⊥ the polar-body axis, relative to the cell COM).

The aligned-outlines overlay draws every zygote's cell cross-section; this adds
`pb_outline` so the front-end can overlay the polar-body silhouettes too. Because
the frame is ⊥ the polar-body axis, the silhouettes sit near the centre (the
polar body is on that axis) — that is the honest geometry, and the reason a
silhouette (not a scattered marker) is the right representation.

Surgical, like add_volumes_to_cross_agg.py: reads the built scenes, adds one
field, preserves everything else (vp/vt/gp). Self-validates by recomputing each
cell outline from its mesh and checking it matches the stored one, so a wrong
frame cannot pass silently.

Run from the deploy repo root:  python3 scripts/add_polar_silhouette_to_cross_agg.py
"""
import gzip
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
XY_UM = 0.15
SLAB_UM = 6.0
N_ANG = 120

SETS = [
    (os.path.join(ROOT, "data", "zygote"),
     [os.path.join(ROOT, "data", "zygote_cross.json.gz"),
      os.path.join(ROOT, "data", "zygote_cross_circ.json.gz")]),
]


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n else v


def mesh_um(verts_plot, zs):
    v = np.asarray(verts_plot, float).reshape(-1, 3)
    return np.stack([v[:, 0] * XY_UM, v[:, 1] * XY_UM, v[:, 2] / zs], axis=1)


def basis_um(uv_plot, zs):
    return unit(np.array([uv_plot[0] * XY_UM, uv_plot[1] * XY_UM, uv_plot[2] / zs], float))


def star_outline(uu, vv, n_ang=N_ANG):
    """Max-radius-per-angular-bin outline around the origin (matches build_zygote)."""
    ang = np.arctan2(vv, uu)
    rad = np.hypot(uu, vv)
    bins = ((ang + np.pi) / (2 * np.pi) * n_ang).astype(int) % n_ang
    out = []
    for b in range(n_ang):
        m = bins == b
        if not m.any():
            continue
        r = float(rad[m].max())
        th = (b + 0.5) / n_ang * 2 * np.pi - np.pi
        out.append([r * np.cos(th), r * np.sin(th)])
    return out


def convex_hull(uu, vv):
    """Closed 2D convex-hull polygon of the projected points (the silhouette)."""
    from scipy.spatial import ConvexHull
    pts = np.column_stack([uu, vv])
    if len(pts) < 3:
        return []
    try:
        h = ConvexHull(pts)
    except Exception:                                              # noqa: BLE001
        return []
    poly = [[float(pts[i, 0]), float(pts[i, 1])] for i in h.vertices]
    poly.append(poly[0])
    return poly


def scene_frame(scene):
    A = scene["analysis"]
    zs = scene["z_scale"]
    com = np.array(A["com_um"], float)
    u = basis_um(A["cross_section"]["u_plot"], zs)
    v = basis_um(A["cross_section"]["v_plot"], zs)
    a = unit(np.array([A["axis_plot"][0] * XY_UM, A["axis_plot"][1] * XY_UM,
                       A["axis_plot"][2] / zs], float))
    return com, u, v, a, zs


def cell_outline_from_mesh(scene, com, u, v, a):
    """Recompute the cell cross-section outline from the mesh, to validate the frame."""
    m = scene["region_meshes"].get("1")
    if not m:
        return None
    verts = mesh_um(m["verts"], scene["z_scale"])
    rel = verts - com
    along = rel @ a
    slab = rel[np.abs(along) < SLAB_UM]
    if len(slab) < 20:
        slab = rel[np.abs(along) < SLAB_UM * 3]
    if len(slab) < 20:
        return None
    return np.array(star_outline(slab @ u, slab @ v))


def outline_mismatch(recomputed, stored):
    """Max radial deviation (µm) between two star outlines, matched by angle bin."""
    if recomputed is None or not stored:
        return None
    st = np.array(stored, float)
    ra = np.arctan2(recomputed[:, 1], recomputed[:, 0]); rr = np.hypot(recomputed[:, 0], recomputed[:, 1])
    sa = np.arctan2(st[:, 1], st[:, 0]); sr = np.hypot(st[:, 0], st[:, 1])
    order = np.argsort(sa)
    interp = np.interp(ra, sa[order], sr[order], period=2 * np.pi)
    return float(np.median(np.abs(rr - interp)))


def pb_silhouette(scene, com, u, v):
    A = scene["analysis"]
    pbk = str(A["polar_body_label"])
    m = scene["region_meshes"].get(pbk)
    if not m:
        return None, None
    verts = mesh_um(m["verts"], scene["z_scale"])
    rel = verts - com
    uu, vv = rel @ u, rel @ v
    poly = convex_hull(uu, vv)
    centroid = [float(uu.mean()), float(vv.mean())]
    return poly, centroid


def main():
    scenes_cache = {}
    for scenes_dir, agg_paths in SETS:
        for agg_path in agg_paths:
            if not os.path.isfile(agg_path):
                print(f"  skip (absent): {os.path.relpath(agg_path, ROOT)}")
                continue
            agg = json.load(gzip.open(agg_path, "rt"))
            n_ok = n_pb = 0
            worst = 0.0
            for e in agg["embryos"]:
                sp = os.path.join(scenes_dir, f"{e['id']}.json.gz")
                if not os.path.isfile(sp):
                    continue
                if sp not in scenes_cache:
                    scenes_cache[sp] = json.load(gzip.open(sp, "rt"))
                scene = scenes_cache[sp]
                com, u, v, a, zs = scene_frame(scene)
                # validate the frame against the stored cell outline
                mm = outline_mismatch(cell_outline_from_mesh(scene, com, u, v, a), e.get("outline"))
                if mm is not None:
                    worst = max(worst, mm); n_ok += 1
                poly, centroid = pb_silhouette(scene, com, u, v)
                if poly:
                    e["pb_outline"] = [[round(x, 2), round(y, 2)] for x, y in poly]
                    e["pb_centroid"] = [round(centroid[0], 2), round(centroid[1], 2)]
                    n_pb += 1
            with gzip.open(agg_path, "wt") as fh:
                json.dump(agg, fh, separators=(",", ":"))
            print(f"  {os.path.relpath(agg_path, ROOT)}: "
                  f"pb_outline added to {n_pb}/{len(agg['embryos'])} · "
                  f"frame validated on {n_ok} (worst cell-outline mismatch {worst:.2f} µm)")
    print("done")


if __name__ == "__main__":
    main()
