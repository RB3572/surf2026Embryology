#!/usr/bin/env python3
"""
Build data/stage_expr.json.gz for the Stage Expression Explorer.

For every sample in the three developmental stages Harry's page compares — Zygote,
Early 2-cell, Late 2-cell (Oocyte excluded to match) — this precomputes, from the
already-built segments scenes:

  * a 2-D silhouette of the cell body, projected with the same camera direction the
    site's 3-D viewer uses, so the grid thumbnails read as mini versions of the live
    view;
  * the cell volume (for volume-density normalization) and the sample's total
    transcript count (for CPM);
  * per gene, the raw count in that sample and a downsampled set of 2-D-projected
    transcript dots (in the same silhouette frame) for the grid.

Counts reproduce segments_genes.json.gz exactly (validated against Harry's Pard3
numbers). The page derives CPM = count / total_tx * 1e6 and density = count / volume
on the fly, so no normalization is baked in here.

Run from the deploy repo root:  python3 build_stage_expr.py
"""
from __future__ import annotations

import glob
import gzip
import json
import os

import numpy as np
from scipy.spatial import ConvexHull

from embryo_naming import embryo_label

ROOT = os.path.dirname(os.path.abspath(__file__))
SEG_DIR = os.path.join(ROOT, "data", "segments")
GENES = os.path.join(ROOT, "data", "segments_genes.json.gz")
MANIFEST = os.path.join(ROOT, "data", "segments_manifest.json")
OUT = os.path.join(ROOT, "data", "stage_expr.json.gz")

# The three stages the explorer compares, in developmental order, with the colours
# Harry's page uses so the two read as siblings.
STAGES = ["Zygote", "Early2Cell", "Late2Cell"]
STAGE_LABEL = {"Zygote": "Zygote", "Early2Cell": "Early 2-cell", "Late2Cell": "Late 2-cell"}
STAGE_COLOR = {"Zygote": "#55B3CB", "Early2Cell": "#D56E2E", "Late2Cell": "#574AA8"}

DOT_CAP = 70            # max transcript dots per gene per sample in the grid
HULL_MAX = 54           # max silhouette polygon points
SILH_LABEL = "1"        # cytoplasm / whole-cell body mesh label
SEED = 7


def camera_basis():
    """Screen basis for the viewer's default camera (eye ≈ (1.5, 1.5, 1.15))."""
    eye = np.array([1.5, 1.5, 1.15], float)
    fwd = eye / np.linalg.norm(eye)                 # look toward origin
    up0 = np.array([0.0, 0.0, 1.0])
    right = np.cross(up0, fwd); right /= np.linalg.norm(right)
    up = np.cross(fwd, right)
    return right, up


RIGHT, UP = camera_basis()


def _norm_scene(verts, ex):
    """Aspect-normalize like the viewer: center, then scale each axis by span/maxspan."""
    lo = np.array([ex["x"][0], ex["y"][0], ex["z"][0]], float)
    hi = np.array([ex["x"][1], ex["y"][1], ex["z"][1]], float)
    span = np.maximum(hi - lo, 1e-6)
    m = span.max()
    c = (lo + hi) / 2
    return (verts - c) / m * (span / m)             # aspect-scaled, centered


def project(verts, ex):
    n = _norm_scene(np.asarray(verts, float), ex)
    u = n @ RIGHT
    v = n @ UP
    return np.column_stack([u, v])


def mesh_volume(verts, faces):
    """Signed volume via the divergence theorem (in scene pixel/z units)."""
    v = np.asarray(verts, float).reshape(-1, 3)
    f = np.asarray(faces, int).reshape(-1, 3)
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    return float(abs(np.sum(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0))


def silhouette(verts2d):
    try:
        h = ConvexHull(verts2d)
        poly = verts2d[h.vertices]
    except Exception:
        return None
    if len(poly) > HULL_MAX:
        idx = np.linspace(0, len(poly) - 1, HULL_MAX).astype(int)
        poly = poly[idx]
    return poly


def frame(sample_pts_list):
    """Common [0,1] framing box for a sample from all its point sets."""
    allpts = np.vstack([p for p in sample_pts_list if p is not None and len(p)])
    lo = allpts.min(axis=0); hi = allpts.max(axis=0)
    span = np.maximum(hi - lo, 1e-6)
    s = span.max()
    # center within a square, margin 6%
    cen = (lo + hi) / 2
    def to01(p):
        return np.clip((p - cen) / s * 0.88 + 0.5, 0.0, 1.0)
    return to01


def scene_path(eid):
    p = os.path.join(SEG_DIR, eid + ".json.gz")
    return p if os.path.exists(p) else None


def main():
    rng = np.random.default_rng(SEED)
    genes = json.load(gzip.open(GENES, "rt"))
    embInfo = genes["embInfo"]
    geneMap = genes["genes"]

    # keep only the 3 comparison stages; map old index -> new sample index
    keep = [(i, e) for i, e in enumerate(embInfo) if e["stage"] in
            (STAGE_LABEL[s] for s in STAGES)]
    old2new = {i: n for n, (i, _e) in enumerate(keep)}
    label2stage = {v: k for k, v in STAGE_LABEL.items()}

    samples = []
    tx_by_sample = {}          # new idx -> {gene: projected transcript array}
    n_ok = 0
    for i, e in keep:
        eid = e["id"]
        sp = scene_path(eid)
        hull = None; vol = 0.0; total_tx = 0; to01 = None
        gene_tx_here = {}
        if sp:
            sc = json.load(gzip.open(sp, "rt"))
            ex = sc["extents"]
            mesh = sc["region_meshes"].get(SILH_LABEL)
            if mesh:
                verts = np.asarray(mesh["verts"], float).reshape(-1, 3)
                proj = project(verts, ex)
                vol = mesh_volume(mesh["verts"], mesh["faces"])
                # gather this sample's transcripts (all genes) for a shared frame
                zs = float(sc.get("z_scale", 1.0))
                tx_raw = {}
                for g, t in sc.get("transcripts", {}).items():
                    if not t or "x" not in t:
                        continue
                    # z is stored as gz; the viewer renders it as gz * z_scale, which is
                    # the same space as the mesh verts, so project it the same way.
                    pts = np.column_stack([t["x"], t["y"], np.asarray(t["gz"], float) * zs])
                    total_tx += len(pts)
                    p2 = project(pts, ex)
                    if len(p2) > DOT_CAP:
                        sel = rng.choice(len(p2), DOT_CAP, replace=False)
                        p2 = p2[sel]
                    tx_raw[g] = p2
                to01 = frame([proj] + list(tx_raw.values()) or [proj])
                hpoly = silhouette(proj)
                hull = to01(hpoly) if hpoly is not None else None
                for g, p2 in tx_raw.items():
                    gene_tx_here[g] = np.round(to01(p2) * 1000).astype(int)  # 0..1000 ints
            n_ok += 1
        samples.append({
            "id": eid, "label": embryo_label(eid),
            "stage": label2stage.get(e["stage"], e["stage"]),
            "hull": (np.round(hull * 1000).astype(int).tolist() if hull is not None else None),
            "vol": round(vol, 1), "total_tx": int(total_tx),
        })
        tx_by_sample[old2new[i]] = gene_tx_here

    # per-gene raw counts by NEW sample index (validated to reproduce Harry's numbers)
    gene_counts = {}
    gene_tx = {}
    for g, entries in geneMap.items():
        counts = {}
        dots = {}
        for embIdx, segs in entries:
            if embIdx not in old2new:
                continue
            ni = old2new[embIdx]
            counts[str(ni)] = int(segs[0][3])       # embTotal
            arr = tx_by_sample.get(ni, {}).get(g)
            if arr is not None and len(arr):
                dots[str(ni)] = arr.tolist()
        if counts:
            gene_counts[g] = counts
            if dots:
                gene_tx[g] = dots

    out = {
        "stages": STAGES,
        "stage_labels": STAGE_LABEL,
        "stage_colors": STAGE_COLOR,
        "samples": samples,
        "gene_counts": gene_counts,
        "gene_tx": gene_tx,
        "genes_all": sorted(gene_counts.keys(), key=str.lower),
        "dot_cap": DOT_CAP,
        "note": "counts from segments_genes.json.gz; hulls/dots projected with the viewer camera",
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(out, fh, separators=(",", ":"))
    from collections import Counter
    sc = Counter(s["stage"] for s in samples)
    print(f"samples: {len(samples)} ({dict(sc)}), scenes rendered: {n_ok}")
    print(f"genes: {len(gene_counts)}  ·  {os.path.getsize(OUT)/1024/1024:.2f} MB")


if __name__ == "__main__":
    main()
