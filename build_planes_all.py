#!/usr/bin/env python3
"""
Build the "Division Plane Sweep — every plane" project data.

This is the UNRESTRICTED sibling of build_zygote.py. Instead of testing 18 planes
that all contain the polar-body axis, it searches EVERY plane orientation through
the cell centre of mass — the plane normal swept over the whole (hemi)sphere at
~1° spacing (~20 000 candidate planes) — and reports, for each gene, that gene's
single globally-best dividing plane (its own 3-D normal), under two normalizations
(count and segment-1-volume).

Because a plane and its flipped normal describe the same partition, the search set
is a Fibonacci HEMISPHERE of unit normals (shared across all zygotes, so a plane is
referenced everywhere by its integer index into that one grid).

Honesty under a huge search: the best-of-~20 000-planes split is large even for
random data, so per-gene significance is a SEARCH-CORRECTED empirical p — the null
is the distribution of the best-plane split for n uniformly-random in-cell points
(CSR), evaluated over a coarse plane grid and amortized across genes by transcript
count n (a monotone reference computed at log-spaced n-anchors, interpolated).

Geometry, units, polar-body detection, cross-sections, the balloon variant and the
cross-aggregate writer are all reused verbatim from build_zygote (imported as BZ);
only the plane SET and the best-plane search differ.

Outputs (deploy repo root):
  data/planes_all/<id>.json.gz     per-embryo render scene + per-gene best planes
  data/planes_all_manifest.json    embryo list + the shared normal grid size
  data/planes_all_normals.json.gz  the shared M×3 unit-normal grid (µm)
  data/planes_all_cross.json.gz    cross-embryo aggregate for the bottom drawer
Run from the deploy repo root:  python3 build_planes_all.py
"""
import glob
import gzip
import json
import os
from collections import Counter

import numpy as np

import build_zygote as BZ           # reuse ALL the heavy machinery
from build_zygote import (XY_UM, Z_UM, N_NULL, unit, mask_and_transcripts,
                          detect_polar_body, cross_section_outline, _json_default)
from embryo_naming import embryo_label

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = BZ.ATLAS
SRC = BZ.SRC
OUT_DIR = os.path.join(HERE, "data", "planes_all")
OUT_MANIFEST = os.path.join(HERE, "data", "planes_all_manifest.json")
OUT_NORMALS = os.path.join(HERE, "data", "planes_all_normals.json.gz")
OUT_CROSS = os.path.join(HERE, "data", "planes_all_cross.json.gz")

M_PLANES = 20000       # candidate plane normals over a hemisphere (~1° spacing) — the search grid
M_NULL = 2000          # grid for the search-corrected p (both the observed-max AND the CSR-max
                       # are taken over THIS grid, so they are consistent; ~2.5° spacing already
                       # resolves the smooth split field, so max-over-2000 ≈ max-over-20000)
NULL_B = 160           # CSR draws per n-anchor
NULL_ANCHORS = 14      # log-spaced transcript-count anchors for the amortized null
N_SURF = 320           # cell-surface points shipped per zygote (browser aligned-outlines)
CHUNK = 4096           # plane-block size for the big (points × planes) projections
RNG = np.random.default_rng(20260731)
BEST_KEYS = BZ.BEST_KEYS          # ["pVol","pCnt","diffVol","diffCnt"]


# ─────────────────────────── plane sets ───────────────────────────
def fib_hemisphere(m):
    """~m near-uniform unit normals on a hemisphere (Fibonacci sphere, kept where the
    canonical hemisphere test holds so a normal and its antipode aren't both present —
    n and −n describe the same plane)."""
    n = 2 * m + 8                                  # oversample the full sphere, then halve
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)                 # polar angle, uniform in cos
    gold = np.pi * (1 + 5 ** 0.5)
    theta = gold * i
    xyz = np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)], axis=1)
    # canonical hemisphere: keep z>0, or (z==0 & y>0), or (z==0 & y==0 & x>0)
    z, y, x = xyz[:, 2], xyz[:, 1], xyz[:, 0]
    keep = (z > 1e-9) | ((np.abs(z) <= 1e-9) & (y > 1e-9)) | ((np.abs(z) <= 1e-9) & (np.abs(y) <= 1e-9) & (x > 0))
    H = xyz[keep]
    return (H / np.linalg.norm(H, axis=1, keepdims=True)).astype(np.float64)


def side_counts(P, com, N):
    """Side-A counts of points P against every plane normal in N (through com).
    Chunked over planes to bound memory. Returns int array (len N)."""
    D = (P - com)
    out = np.empty(len(N), np.int64)
    for s in range(0, len(N), CHUNK):
        blk = N[s:s + CHUNK]
        out[s:s + len(blk)] = (D @ blk.T > 0).sum(axis=0)
    return out


def in_plane_basis(nrm):
    """An orthonormal (t, w) spanning the plane with unit normal nrm (t is a stable
    in-plane 'vertical'; w = nrm × t is the viewing direction that makes the plane edge-on)."""
    ref = np.array([0.0, 0.0, 1.0]) if abs(nrm[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    t = unit(np.cross(nrm, ref))
    w = unit(np.cross(nrm, t))
    return t, w


# ─────────────────── search-corrected CSR null (per zygote) ───────────────────
def build_null_reference(pos1, com, voxvol, Nn, volA_n, volB_n):
    """Amortized CSR null of the BEST-plane split over a coarse plane grid Nn.
    For each of NULL_ANCHORS log-spaced n-values, draw NULL_B sets of n random seg-1
    voxels, and record the distribution of max_k|dNorm| and max_k|dVol|. Returns
    (anchors[n], cnt_null[A,B], vol_null[A,B]) — sorted null arrays per anchor for
    interpolated empirical p by transcript count."""
    nv = len(pos1)
    D1 = pos1 - com                                        # (nv,3)
    # precompute side membership of every seg-1 voxel against every null plane, once
    S = np.empty((nv, len(Nn)), bool)
    for s in range(0, len(Nn), CHUNK):
        blk = Nn[s:s + CHUNK]
        S[:, s:s + len(blk)] = (D1 @ blk.T) > 0
    lo = 5
    hi = max(lo + 1, nv)
    anchors = np.unique(np.round(np.geomspace(lo, min(hi, 20000), NULL_ANCHORS)).astype(int))
    cnt_null = np.zeros((len(anchors), NULL_B))
    vol_null = np.zeros((len(anchors), NULL_B))
    for ai, n in enumerate(anchors):
        for b in range(NULL_B):
            idx = RNG.integers(0, nv, size=int(n))
            aK = S[idx].sum(axis=0).astype(float)          # (M_null,) side-A counts
            bK = n - aK
            cnt_null[ai, b] = np.max(np.abs(aK - bK) / n)
            vol_null[ai, b] = np.max(np.abs(aK / volA_n - bK / volB_n))
    cnt_null.sort(axis=1); vol_null.sort(axis=1)
    return anchors, cnt_null, vol_null


def _emp_p(sorted_null_row, obs):
    """Empirical upper-tail p from a sorted null sample: (1+#{null>=obs})/(1+B)."""
    B = len(sorted_null_row)
    ge = B - int(np.searchsorted(sorted_null_row, obs, side="left"))
    return (1 + ge) / (1 + B)


def null_p(anchors, null_rows, n, obs):
    """Interpolated search-corrected p for a gene of size n and observed best-split obs.
    Interpolates the empirical p between the two bracketing n-anchors (log scale)."""
    j = int(np.searchsorted(anchors, n))
    if j <= 0:
        return _emp_p(null_rows[0], obs)
    if j >= len(anchors):
        return _emp_p(null_rows[-1], obs)
    a0, a1 = anchors[j - 1], anchors[j]
    p0, p1 = _emp_p(null_rows[j - 1], obs), _emp_p(null_rows[j], obs)
    w = (np.log(n) - np.log(a0)) / (np.log(a1) - np.log(a0)) if a1 > a0 else 0.0
    return float(min(1.0, max(0.0, p0 + w * (p1 - p0))))


# ─────────────────────────── per-embryo build ───────────────────────────
def process(eid, N, Nn):
    scene_p = os.path.join(ATLAS, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(SRC, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    d = json.load(gzip.open(scene_p, "rt"))
    zs = d.get("z_scale", 7.0)
    tx = d.get("transcripts", {})
    genes = sorted(tx.keys(), key=lambda g: -d.get("gene_totals", {}).get(g, 0))

    pos, labels, voxvol, seg_of, inside_fractions = mask_and_transcripts(lab[0], tx, genes)
    if 1 not in labels:
        return None
    pb_label, pb_det = detect_polar_body(d, labels, zs, inside_fractions)
    if pb_label is None:
        print(f"  -- {eid}: {pb_det['reason']}")
        return None
    com = pos.mean(axis=0)
    pb_com = pos[labels == pb_label].mean(axis=0)
    a = unit(pb_com - com)
    pos1 = pos[labels == 1]

    # per-plane segment-1 volumes over the FULL search grid + the coarse null grid
    volA = side_counts(pos1, com, N).astype(float) * voxvol
    volB = (len(pos1) - side_counts(pos1, com, N)).astype(float) * voxvol
    volA = np.maximum(volA, voxvol); volB = np.maximum(volB, voxvol)
    aNv = side_counts(pos1, com, Nn).astype(float)
    volA_n = np.maximum(aNv * voxvol, voxvol)
    volB_n = np.maximum((len(pos1) - aNv) * voxvol, voxvol)

    anchors, cnt_null, vol_null = build_null_reference(pos1, com, voxvol, Nn, volA_n, volB_n)

    gene_rows = []
    weff_cnt = np.zeros(len(N)); weff_vol = np.zeros(len(N)); total_n = 0
    for gi, g in enumerate(genes):
        t = tx[g]; in1 = seg_of[g] == 1
        P = np.stack([np.asarray(t["x"], float)[in1] * XY_UM,
                      np.asarray(t["y"], float)[in1] * XY_UM,
                      np.asarray(t["gz"], float)[in1] * Z_UM], axis=1)
        n = len(P)
        if n == 0:
            continue
        # fine grid (~1°): the reported best plane geometry + effect size
        aK = side_counts(P, com, N).astype(float)
        bK = n - aK
        dCnt = np.abs(aK - bK) / n                    # |dNorm| per plane
        dVol = np.abs(aK / volA - bK / volB)          # |dVol| per plane
        i_cnt = int(np.argmax(dCnt))
        i_vol = int(np.argmax(dVol))
        # coarse grid (matches the null grid): the SEARCH-CORRECTED p is computed from the
        # observed max over the SAME grid the CSR-max null used, so they are directly comparable.
        aKn = side_counts(P, com, Nn).astype(float)
        obs_cnt = float(np.max(np.abs(aKn - (n - aKn)) / n))
        obs_vol = float(np.max(np.abs(aKn / volA_n - (n - aKn) / volB_n)))
        p_cnt = null_p(anchors, cnt_null, n, obs_cnt)
        p_vol = null_p(anchors, vol_null, n, obs_vol)
        # coin-flip representative for the mini counts chart (one fair-flip realization)
        na = int(RNG.binomial(n, 0.5))
        gene_rows.append({
            "gene": g, "idx": gi, "total": n,
            "iCnt": i_cnt, "iVol": i_vol,
            "aCnt": int(aK[i_cnt]), "aVol": int(aK[i_vol]),
            "volA_cnt": round(float(volA[i_cnt]), 1), "volB_cnt": round(float(volB[i_cnt]), 1),
            "volA_vol": round(float(volA[i_vol]), 1), "volB_vol": round(float(volB[i_vol]), 1),
            "effCnt": round(float(dCnt[i_cnt]), 5), "effVol": round(float(dVol[i_vol]), 7),
            "pCnt": round(float(p_cnt), 5), "pVol": round(float(p_vol), 5),
            "na": na,
        })
        weff_cnt += n * dCnt; weff_vol += n * dVol; total_n += n

    if total_n == 0:
        return None
    weff_cnt /= total_n; weff_vol /= total_n
    emb_i_cnt = int(np.argmax(weff_cnt))              # transcript-weighted embryo-consensus plane
    emb_i_vol = int(np.argmax(weff_vol))
    best_idx = {"pVol": emb_i_vol, "pCnt": emb_i_cnt, "diffVol": emb_i_vol, "diffCnt": emb_i_cnt}

    # ── render geometry (plot space) ──
    com_plot = [com[0] / XY_UM, com[1] / XY_UM, com[2] * zs]
    pb_plot = [pb_com[0] / XY_UM, pb_com[1] / XY_UM, pb_com[2] * zs]
    ex = d["extents"]
    L_um = 0.62 * 0.5 * max(ex["x"][1] - ex["x"][0], ex["y"][1] - ex["y"][0],
                            ex["z"][1] - ex["z"][0]) * XY_UM
    # cell-surface points (µm, relative to COM) for the browser aligned-outlines figure
    m1 = d["region_meshes"].get("1")
    surf = []
    if m1:
        V = np.asarray(m1["verts"], float).reshape(-1, 3)
        Vum = np.stack([V[:, 0] * XY_UM, V[:, 1] * XY_UM, (V[:, 2] / zs) * Z_UM], axis=1) - com
        if len(Vum) > N_SURF:
            sel = np.linspace(0, len(Vum) - 1, N_SURF).astype(int)
            Vum = Vum[sel]
        surf = [[round(float(c), 2) for c in p] for p in Vum]

    # equatorial (⊥ axis) outline for the 3-D cross-section fallback + axis display
    ref = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = unit(np.cross(a, ref)); v = unit(np.cross(a, u))
    outline = cross_section_outline(pos1, com, a, u, v)

    analysis = {
        "com_plot": [round(c, 2) for c in com_plot], "com_um": [round(float(c), 4) for c in com],
        "pb_plot": [round(c, 2) for c in pb_plot], "polar_body_label": int(pb_label),
        "polar_body_detection": pb_det,
        "axis_plot": [round(x, 5) for x in [a[0] / XY_UM, a[1] / XY_UM, a[2] * zs]],
        "axis_um": [round(float(x), 6) for x in a],       # polar-body axis unit (µm) for the aligned/orientation views
        "best_idx": best_idx, "n_null": NULL_B, "L_um": round(float(L_um), 2),
        "cross_section": {"u_plot": [u[0] / XY_UM, u[1] / XY_UM, u[2] * zs],
                          "v_plot": [v[0] / XY_UM, v[1] / XY_UM, v[2] * zs],
                          "outline": [[round(p[0], 2), round(p[1], 2)] for p in outline]},
        "surf_um": surf, "genes": gene_rows,
    }
    scene = {
        "id": eid, "z_scale": zs, "extents": ex,
        "region_meshes": d["region_meshes"], "region_defaults": d["region_defaults"],
        "mask_labels": d["mask_labels"], "genes": genes, "gene_totals": d.get("gene_totals", {}),
        "transcripts": {g: {"x": tx[g]["x"], "y": tx[g]["y"], "gz": tx[g]["gz"],
                            "s1": (seg_of[g] == 1).astype(np.uint8).tolist()} for g in tx},
        "analysis": analysis, "circ": None,
    }
    return scene


def agg_entry(eid, label, A):
    """Slim per-embryo record for the cross-embryo bottom drawer.
    `gb[gene]` mirrors build_zygote's layout so the front-end reuses geneAlign():
      [ iVol, iCnt, iVol, iCnt,  aVol, aCnt, aVol, aCnt,  pVol, pCnt, pVol, pCnt ]
    (only two distinct planes exist per gene — the volume-best and the count-best —
    duplicated into the 4 BEST_KEYS slots). `g[gene]` = [total, a@each best]. `vp`/`vt`
    give per-side seg-1 volumes for the density modes. `surf` = cell-surface µm points
    (relative to COM) so the aligned-outlines figure can be re-projected per gene."""
    best = [A["best_idx"][k] for k in BEST_KEYS]
    g, gb, vp, vt = {}, {}, {}, {}
    for row in A["genes"]:
        gn = row["gene"]
        aV, aC = row["aVol"], row["aCnt"]
        g[gn] = [row["total"], aV, aC, aV, aC]
        gb[gn] = ([row["iVol"], row["iCnt"], row["iVol"], row["iCnt"]]
                  + [aV, aC, aV, aC]
                  + [row["pVol"], row["pCnt"], row["pVol"], row["pCnt"]])
        vp[gn] = [row["volA_vol"], row["volA_cnt"], row["volA_vol"], row["volA_cnt"]]
        vt[gn] = [round(row["volA_vol"] + row["volB_vol"], 1),
                  round(row["volA_cnt"] + row["volB_cnt"], 1),
                  round(row["volA_vol"] + row["volB_vol"], 1),
                  round(row["volA_cnt"] + row["volB_cnt"], 1)]
    sig = {k: 0.0 for k in BEST_KEYS}
    return {"id": eid, "label": label, "outline": A["cross_section"]["outline"],
            "surf": A.get("surf_um", []), "axis": A.get("axis_um"), "best": best, "sig": sig,
            "g": g, "gb": gb, "vp": vp, "vt": vt}


def write_cross(entries, normals, path):
    cov, sums = Counter(), Counter()
    for e in entries:
        for gn, arr in e["g"].items():
            cov[gn] += 1; sums[gn] += arr[0]
    genes_all = sorted(cov.keys(), key=str.lower)
    default_align = max(genes_all, key=lambda gn: (cov[gn], sums[gn])) if genes_all else None
    agg = {"m_planes": len(normals), "best_keys": BEST_KEYS, "n_embryos": len(entries),
           "genes_all": genes_all, "gene_cov": {gn: cov[gn] for gn in genes_all},
           "default_align_gene": default_align,
           "normals": [[round(float(x), 5) for x in nrm] for nrm in normals],
           "embryos": entries}
    with gzip.open(path, "wt") as fh:
        json.dump(agg, fh, separators=(",", ":"), default=_json_default)
    print(f"  cross-aggregate: {len(entries)} embryos, {len(genes_all)} union genes, "
          f"default align = {default_align}  ({os.path.getsize(path)/1024:.0f} KB)")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    N = fib_hemisphere(M_PLANES)
    Nn = fib_hemisphere(M_NULL)
    print(f"plane grid: {len(N)} search normals, {len(Nn)} null normals")
    with gzip.open(OUT_NORMALS, "wt") as fh:
        json.dump({"normals": [[round(float(x), 5) for x in nrm] for nrm in N]}, fh, separators=(",", ":"))

    ids = sorted(os.listdir(ATLAS))
    manifest, agg_entries = [], []
    MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    for i, eid in enumerate(ids):
        try:
            scene = process(eid, N, Nn)
        except Exception as e:                       # noqa: BLE001
            print(f"  !! {eid}: {e}")
            continue
        if not scene:
            print(f"  -- skipped {eid}")
            continue
        out = os.path.join(OUT_DIR, eid + ".json.gz")
        with gzip.open(out, "wt") as fh:
            json.dump(scene, fh, separators=(",", ":"), default=_json_default)
        date = eid[:8]
        ds = f"{MON[int(date[4:6]) - 1]} {int(date[6:8])}" if date.isdigit() else ""
        label = embryo_label(eid, "zygote")
        manifest.append({
            "id": eid, "label": label, "date_short": ds,
            "n_genes": len(scene["genes"]),
            "n_transcripts": sum(len(t["x"]) for t in scene["transcripts"].values()),
            "size_kb": round(os.path.getsize(out) / 1024),
            "polar_body_label": scene["analysis"]["polar_body_label"],
        })
        agg_entries.append(agg_entry(eid, label, scene["analysis"]))
        print(f"  [{i+1}/{len(ids)}] {eid}  {manifest[-1]['n_genes']}g  {manifest[-1]['size_kb']}KB")
    manifest.sort(key=lambda m: m["id"]); agg_entries.sort(key=lambda e: e["id"])
    with open(OUT_MANIFEST, "w") as fh:
        json.dump({"embryos": manifest, "m_planes": len(N), "m_null": len(Nn),
                   "null_b": NULL_B}, fh, indent=1)
    write_cross(agg_entries, N, OUT_CROSS)
    tot = sum(m["size_kb"] for m in manifest)
    print(f"\nwrote {len(manifest)} zygotes  ({tot/1024:.1f} MB)")


if __name__ == "__main__":
    main()
