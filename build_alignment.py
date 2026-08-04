#!/usr/bin/env python3
"""Build data/alignment.json.gz — where the sperm sits once you pick an ANCHOR.

THE QUESTION
------------
A 2-cell embryo gives you one unambiguous axis: the line joining the two blastomere centres.
Everything else — which blastomere is "first", and which way is up — is unlabelled. So a pile of
2-cell embryos cannot be overlaid until you choose a rule that orients them. That rule is the
ANCHOR, and the point of this project is that the answer to "where does the sperm sit?" is not a
property of the embryos alone: it is a property of the embryos AND the anchor you picked.

THE FRAME
---------
Per embryo, from the two largest segments (the scene ships voxel volumes, so this is read off,
and it matches build_contact.py so the two projects agree):

    û = unit(COM_B − COM_A),   M = midpoint

That fixes the horizontal. An anchor supplies the two things û cannot:

  * WHICH BLASTOMERE IS RIGHT (+x).  For a gene anchor, the blastomere holding more of it per
    unit volume. For the polar-body anchor, the blastomere whose centre is nearer the polar body.
    This is the same kind of rule build_alphabeta.py calls an alpha/beta method.
  * WHICH WAY IS UP (+y).  The component of the anchor perpendicular to û: the gene's transcript
    centroid, or the polar body. So the cross-section plane is the plane containing û and the
    anchor, and the anchor always lands in the upper half.

Every embryo carrying the anchor can then be drawn in one shared frame, and each sperm gets an
angle θ about its own blastomere's centre. θ is the number the rankings and the heatmap are about.

⚠️  THE PROBESETS ARE DISJOINT, AND HERE IT REALLY HURTS
--------------------------------------------------------
A gene is only measured in the embryos whose panel contains it, and no gene is in every panel.
Across the 93 two-cell embryos the best-covered gene reaches 48; the median reaches 21. Of the 21
embryos that ALSO have a labelled sperm, the best gene anchor reaches 13 and most reach 7-10.
So a gene anchor is never a statement about all 2-cell embryos — it is a statement about the
handful it can even orient. Every anchor therefore carries its own n, nothing is ranked below
MIN_SPERM, and the page is expected to show n next to every number.

The polar-body anchor is the exception: it needs no panel, so it covers every embryo that has a
polar-body segment.

OUTLINES
--------
For the two blastomeres a radius map R(t, ψ) is precomputed from the mesh — the distance from the
axis to the surface at along-axis position t and azimuth ψ. A cross-section at azimuth φ is then
just two slices of that map (φ for the top, φ+π for the bottom), so the browser can re-cut every
embryo the moment the anchor changes without shipping meshes. Nuclei and the polar body are small
and near-convex, so they ship as ellipsoids (centroid + covariance); the silhouette of an
ellipsoid in any plane is an exact ellipse, which is cheap and orientation-correct.

Units: x,y pixels × 0.15 = µm; mesh z is PLOT space so ÷ z_scale; transcript gz is ALREADY µm.

Output: data/alignment.json.gz
"""
import csv
import glob
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SRC = os.path.join(DATA, "segments")
SPERM_CSV = os.path.join(HERE, "..", "data", "merfish_sperm.csv")
OUT = os.path.join(DATA, "alignment.json.gz")

VERSION = "alignment-1.0.0"
XY_UM = 0.15
T_BINS = 40                 # along-axis bins of the blastomere radius map
PSI_BINS = 36               # azimuth bins (10°)
MIN_TX = 12                 # transcripts a gene needs in an embryo to orient it
MIN_EMB = 5                 # embryos a gene must orient at all to be offered as an anchor
MIN_SPERM = 4               # sperm-carrying embryos an anchor needs before it may be RANKED
POLAR = "__polar__"         # the polar-body anchor's key


def unit(v):
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-12 else v


def verts_um(mesh, zs):
    v = np.asarray(mesh["verts"], float).reshape(-1, 3)
    return np.stack([v[:, 0] * XY_UM, v[:, 1] * XY_UM, v[:, 2] / zs], axis=1)


def radius_map(V, com_a, com_b, e1, e2):
    """R(t, ψ) for one blastomere: max distance from the embryo axis to the surface.

    t is normalised so 0 = COM_A and 1 = COM_B, then extended a little either side, which keeps
    the two blastomeres on a common scale whatever the embryo's absolute size.
    """
    u = unit(com_b - com_a)
    d = V - com_a
    t = d @ u
    L = float(np.linalg.norm(com_b - com_a)) or 1.0
    tn = t / L
    perp = d - np.outer(t, u)
    rho = np.linalg.norm(perp, axis=1)
    psi = np.arctan2(perp @ e2, perp @ e1)                    # −π..π

    lo, hi = -0.75, 1.75                                       # t range covered by the map
    ti = np.clip(((tn - lo) / (hi - lo) * T_BINS).astype(int), 0, T_BINS - 1)
    pi_ = np.clip(((psi + np.pi) / (2 * np.pi) * PSI_BINS).astype(int), 0, PSI_BINS - 1)
    R = np.zeros((T_BINS, PSI_BINS), float)
    np.maximum.at(R, (ti, pi_), rho)
    # A vertex mesh leaves gaps; fill each t-row's empty azimuths from their neighbours so the
    # sliced outline is continuous rather than dropping to the axis.
    for i in range(T_BINS):
        row = R[i]
        if not row.any():
            continue
        idx = np.nonzero(row)[0]
        if len(idx) < PSI_BINS:
            ang = (idx + 0.5) / PSI_BINS * 2 * np.pi
            allang = (np.arange(PSI_BINS) + 0.5) / PSI_BINS * 2 * np.pi
            # circular interpolation
            ext_a = np.concatenate([ang - 2 * np.pi, ang, ang + 2 * np.pi])
            ext_v = np.concatenate([row[idx], row[idx], row[idx]])
            R[i] = np.interp(allang, ext_a, ext_v)
    return R, lo, hi, L


def ellipsoid(V):
    """centroid + covariance — the silhouette of this in any plane is an exact ellipse."""
    c = V.mean(axis=0)
    d = V - c
    C = (d.T @ d) / max(len(d) - 1, 1)
    return c, C


def load_sperm():
    """resolved id -> (segment label, sperm point in PLOT space from the axes scene)."""
    out = {}
    for r in csv.DictReader(open(SPERM_CSV)):
        eid = (r.get("resolved_embryo_id") or "").strip()
        seg = (r.get("segment") or "").strip()
        if not eid or not seg:
            continue
        try:
            out[eid] = int(float(seg))
        except ValueError:
            pass
    return out


def axes_landmark(eid, key, zs):
    """A curated landmark (sperm_plot, polar_plot) from the axes scene, in µm.

    The axes build already identified these by hand for every sperm-labelled embryo, and it is
    far more reliable than guessing from segment sizes — the size heuristic finds a polar body in
    6 of the 21 sperm-carrying 2-cell embryos, the curated landmark in 18.
    """
    p = os.path.join(DATA, "axes", eid + ".json.gz")
    if not os.path.isfile(p):
        return None
    s = json.load(gzip.open(p, "rt"))
    v = (s.get("landmarks") or {}).get(key)
    if not v:
        return None
    z = s.get("z_scale", zs) or zs
    return np.array([v[0] * XY_UM, v[1] * XY_UM, v[2] / z], float)


def build_embryo(path, sperm_seg):
    s = json.load(gzip.open(path, "rt"))
    eid = s["id"]
    zs = s.get("z_scale") or 7.0
    segs = sorted(s.get("segments") or [], key=lambda d: -d["volume"])
    meshes = s.get("region_meshes") or {}
    if len(segs) < 2:
        return None, "fewer than two segments"
    la, lb = segs[0]["label"], segs[1]["label"]
    if str(la) not in meshes or str(lb) not in meshes:
        return None, "blastomere mesh missing"

    VA = verts_um(meshes[str(la)], zs)
    VB = verts_um(meshes[str(lb)], zs)
    ca, cb = VA.mean(axis=0), VB.mean(axis=0)
    u = unit(cb - ca)
    if not np.isfinite(u).all() or np.linalg.norm(cb - ca) < 5:
        return None, "blastomere centres too close"
    mid = (ca + cb) / 2.0

    # a fixed perpendicular basis for this embryo — azimuths are all measured against it
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(u @ tmp) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    e1 = unit(tmp - (tmp @ u) * u)
    e2 = np.cross(u, e1)

    RA, tlo, thi, L = radius_map(VA, ca, cb, e1, e2)
    RB, _, _, _ = radius_map(VB, ca, cb, e1, e2)

    def q(R):                       # quantise to uint8 against this map's own max
        mx = float(R.max()) or 1.0
        return mx, np.clip(np.round(R / mx * 255), 0, 255).astype(np.uint8).ravel().tolist()

    mxa, qa = q(RA)
    mxb, qb = q(RB)

    # everything that is not a blastomere, as an ellipsoid, ranked by size so the page can style
    # the two big ones as nuclei and anything smaller as a polar body
    blobs = []
    for sg in segs[2:]:
        k = str(sg["label"])
        if k not in meshes:
            continue
        V = verts_um(meshes[k], zs)
        if len(V) < 12:
            continue
        c, C = ellipsoid(V)
        blobs.append({"label": int(sg["label"]), "vol": round(float(sg["volume"]), 1),
                      "c": [round(float(v), 2) for v in c],
                      "cov": [round(float(C[i, j]), 3) for i, j in
                              ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))]})
    blobs.sort(key=lambda b: -b["vol"])

    # per-gene: counts in each blastomere, and the azimuth of the gene's centroid
    genes = {}
    tx = s.get("transcripts") or {}
    for g, rec in tx.items():
        seg = np.asarray(rec.get("s") or [], int)
        if seg.size == 0:
            continue
        x = np.asarray(rec["x"], float) * XY_UM
        y = np.asarray(rec["y"], float) * XY_UM
        z = np.asarray(rec["gz"], float)                      # already µm
        P = np.stack([x, y, z], axis=1)
        na = int((seg == la).sum())
        nb = int((seg == lb).sum())
        if na + nb < MIN_TX:
            continue
        inside = P[(seg == la) | (seg == lb)]
        d = inside.mean(axis=0) - mid
        perp = d - (d @ u) * u
        if np.linalg.norm(perp) < 1e-6:
            continue
        genes[g] = [na, nb, round(float(np.degrees(np.arctan2(perp @ e2, perp @ e1))), 1)]

    out = {
        "id": eid,
        "stage": "e2c" if "Early" in os.path.basename(path) else "l2c",
        "vol_a": round(float(segs[0]["volume"]), 1), "vol_b": round(float(segs[1]["volume"]), 1),
        "label_a": int(la), "label_b": int(lb),
        "com_a": [round(float(v), 2) for v in ca], "com_b": [round(float(v), 2) for v in cb],
        "mid": [round(float(v), 2) for v in mid],
        "u": [round(float(v), 5) for v in u],
        "e1": [round(float(v), 5) for v in e1], "e2": [round(float(v), 5) for v in e2],
        "sep": round(float(np.linalg.norm(cb - ca)), 2),
        "map": {"t0": tlo, "t1": thi, "nt": T_BINS, "npsi": PSI_BINS, "L": round(L, 3),
                "a_max": round(mxa, 3), "a": qa, "b_max": round(mxb, 3), "b": qb},
        "blobs": blobs,
        "genes": genes,
    }

    # the polar body: the curated landmark when the axes build has one, else the smallest blob
    # provided it is clearly smaller than the nuclei rather than just the third-largest thing in
    # a noisy segmentation
    pb = axes_landmark(eid, "polar_plot", zs)
    if pb is not None:
        out["polar"] = [round(float(v), 2) for v in pb]
        out["polar_src"] = "landmark"
    elif len(blobs) >= 3 and blobs[-1]["vol"] < 0.5 * blobs[-2]["vol"]:
        out["polar"] = blobs[-1]["c"]
        out["polar_src"] = "segment"

    seg = sperm_seg.get(eid)
    if seg is not None:
        p = axes_landmark(eid, "sperm_plot", zs)
        if p is not None:
            out["sperm"] = [round(float(v), 2) for v in p]
            out["sperm_seg"] = int(seg)
            # which blastomere: trust the labelled segment, fall back to the nearer centre
            if seg == la:
                out["sperm_side"] = "a"
            elif seg == lb:
                out["sperm_side"] = "b"
            else:
                out["sperm_side"] = "a" if np.linalg.norm(p - ca) <= np.linalg.norm(p - cb) else "b"
    return out, None


def main():
    files = sorted(glob.glob(os.path.join(SRC, "Early2Cell__*.json.gz"))) + \
            sorted(glob.glob(os.path.join(SRC, "Late2Cell__*.json.gz")))
    if not files:
        sys.exit("no 2-cell scenes in data/segments — nothing to do")
    sperm_seg = load_sperm()
    print(f"alignment: scanning {len(files)} two-cell scenes")

    embryos, skipped = [], []
    for p in files:
        e, why = build_embryo(p, sperm_seg)
        if e is None:
            skipped.append((os.path.basename(p), why)); continue
        embryos.append(e)
    print(f"  {len(embryos)} usable, {len(skipped)} skipped")
    for n, w in skipped:
        print(f"    skip {n}: {w}")

    # anchor inventory: which embryos each one can orient, and how many of those carry a sperm
    with_sperm = {e["id"] for e in embryos if "sperm" in e}
    counts = {}
    for e in embryos:
        for g in e["genes"]:
            c = counts.setdefault(g, [0, 0])
            c[0] += 1
            if e["id"] in with_sperm:
                c[1] += 1
    anchors = []
    for g, (n, ns) in sorted(counts.items()):
        if n < MIN_EMB:
            continue
        anchors.append({"key": g, "kind": "gene", "n": n, "n_sperm": ns})
    n_polar = sum(1 for e in embryos if "polar" in e)
    n_polar_sp = sum(1 for e in embryos if "polar" in e and "sperm" in e)
    if n_polar:
        anchors.append({"key": POLAR, "kind": "polar", "n": n_polar, "n_sperm": n_polar_sp})
    anchors.sort(key=lambda a: (-a["n_sperm"], -a["n"], a["key"]))

    meta = {
        "version": VERSION,
        "n_embryos": len(embryos),
        "n_e2c": sum(1 for e in embryos if e["stage"] == "e2c"),
        "n_l2c": sum(1 for e in embryos if e["stage"] == "l2c"),
        "n_sperm": len(with_sperm),
        "n_sperm_e2c": sum(1 for e in embryos if e["stage"] == "e2c" and "sperm" in e),
        "n_sperm_l2c": sum(1 for e in embryos if e["stage"] == "l2c" and "sperm" in e),
        "n_anchors": len(anchors),
        "n_polar": n_polar,
        "min_tx": MIN_TX, "min_emb": MIN_EMB, "min_sperm": MIN_SPERM,
        "t_bins": T_BINS, "psi_bins": PSI_BINS,
        "polar_key": POLAR,
        "best_anchor_sperm": max((a["n_sperm"] for a in anchors), default=0),
    }
    doc = {"meta": meta, "anchors": anchors, "embryos": embryos}
    os.makedirs(DATA, exist_ok=True)
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    mb = os.path.getsize(OUT) / 1e6
    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({mb:.2f} MB)")
    print(f"  {len(embryos)} embryos ({meta['n_e2c']} early + {meta['n_l2c']} late), "
          f"{len(with_sperm)} with a sperm")
    print(f"  {len(anchors)} anchors; best gene anchor orients {meta['best_anchor_sperm']} "
          f"sperm-carrying embryos; polar-body anchor covers {n_polar_sp}")


if __name__ == "__main__":
    main()
