#!/usr/bin/env python3
"""The house definitions, in one place.

Ported from `HighResSlideshowExports/_work/scene_export.py` + `_work/halfspace.py`, which are the
specification for how this project's analysis is conducted. Every build that measures a gene across
a plane should import from here rather than re-deriving it, because the audit found five different
conventions drifting apart across the site.

THE FIVE RULES THIS MODULE EXISTS TO ENFORCE

  1. CYTOPLASM ONLY, EXACTLY. Counts come from the per-molecule segment label the segments scenes
     carry (`transcripts[g]["s"]`), so a pronuclear or polar-body molecule is excluded by identity,
     never by a containment test. Volume is the body label's own voxel volume from the scene's
     `segments` block, which by construction excludes every other label.

  2. THE BODY IS FOUND BY VOLUME, NEVER BY LABEL NUMBER. "Labels are not consistent across
     embryos" — scene_export.classify. Every build that hard-codes `label == 1` is one relabelled
     embryo away from measuring a pronucleus.

  3. SIDES ARE ORIENTED BEFORE ANYTHING IS AVERAGED. A plane's two sides are an arbitrary
     geometric convention with no anatomical meaning, so averaging a signed log-ratio over embryos
     whose sides are labelled arbitrarily cancels a real effect to zero. The house convention is
     the FULLER half — whichever holds more cytoplasmic transcripts across the whole panel. It is
     intrinsic to the embryo and needs no landmark.

  4. THE BULK CORRECTION IS A MEDIAN OF RATIOS. One gene can carry 30% of an embryo's cytoplasmic
     transcripts (up to 83%), and abundant genes are less asymmetric than typical, so dividing by
     the total under-corrects everything else. Subtracting the MEDIAN per-gene log ratio — the same
     idea as DESeq's size factors — centres the typical gene at zero instead.

  5. THE NULL IS VOLUME-MATCHED, NOT A FAIR COIN. Half the cell is not half the volume, so a
     Binomial(n, 0.5) null attributes the cell's own geometry to the gene. The null proportion is
     volA/(volA+volB).

and two more that follow from them:

  · +0.5 PSEUDOCOUNT on every count before a log, so a zero-count half is finite.
  · ACROSS EMBRYOS: a one-sample t-test of the per-embryo log2 fold against zero, with floors
    (MIN_TOTAL summed, MIN_EMBRYOS carrying), then Benjamini-Hochberg.

Nothing here reads the atlas or the label TIFFs: everything comes from data/segments/, which ships
in this repo. That is deliberate — the atlas lives on another machine.
"""
import gzip
import json
import math
import os

import numpy as np

PX = 0.15                 # µm per pixel, in ALL THREE axes of the segments scenes' pixel space
EPS = 0.5                 # pseudocount
MIN_TOTAL = 20            # a gene needs this many cytoplasmic transcripts summed over embryos
MIN_EMBRYOS = 5           # ...in at least this many embryos
N_NULL_DRAWS = 200        # count-matched null draws per gene


# ───────────────────────── scenes ─────────────────────────
def read_scene(path):
    return json.load(gzip.open(path, "rt"))


def seg_volumes(sc):
    """{label(str): volume µm³} from the scene's own segments block."""
    return {str(s["label"]): float(s["volume"]) for s in sc.get("segments", [])}


def classify_body(sc):
    """The cell-body label(s), by VOLUME — never by label number.

    A zygote has one body then a big step down to the pronuclei; a 2-cell embryo has two
    comparable ones. 0.4 sits in the empty gap between those cases (real sister blastomeres are
    within ~2x of each other, a pronucleus is under a tenth of the cell)."""
    vol = seg_volumes(sc)
    if not vol:
        return []
    order = sorted(vol, key=lambda k: vol[k], reverse=True)
    two = len(order) > 1 and vol[order[1]] / vol[order[0]] > 0.4
    return order[:2] if two else order[:1]


def mesh_of(sc, label):
    """A segment's mesh in µm."""
    r = sc["region_meshes"][str(label)]
    return (np.asarray(r["verts"], float).reshape(-1, 3) * PX,
            np.asarray(r["faces"], int).reshape(-1, 3))


def mesh_volume(V, F):
    T = V[F]
    return abs(float(np.einsum("ij,ij->i", T[:, 0], np.cross(T[:, 1], T[:, 2])).sum()) / 6.0)


def vol_centroid(V, F):
    """Centroid of the enclosed solid — not the mean vertex, which uneven triangulation biases."""
    T = V[F]
    w = np.einsum("ij,ij->i", T[:, 0], np.cross(T[:, 1], T[:, 2])) / 6.0
    return (((T[:, 0] + T[:, 1] + T[:, 2]) / 4.0) * w[:, None]).sum(0) / w.sum()


PB_CONE_DEG = 25.0        # cone half-angle used to measure the cortex in a candidate's direction
PB_MIN_REACH = 0.90       # a polar body's centroid must reach ~the cortex, not sit inside it


def polar_label(sc, exclude=()):
    """The polar body, or None.

    ⚠️ RETURNING None IS THE POINT. A zygote with no polar body has no animal–vegetal axis, so
    every meridional and equatorial construction is undefined for it. An earlier version of this
    just took the most peripheral non-body segment, which on the nine zygotes that genuinely have
    none picked a PRONUCLEUS (radial ratio 0.35-0.76, i.e. well inside the cell) and handed back a
    confident, meaningless axis.

    So a candidate must actually REACH THE CORTEX: its centroid, projected onto its own direction
    from the body centre, must be at least PB_MIN_REACH of how far the body surface extends in
    that same direction (measured in a cone, so an elongated cell is not judged by its long axis).
    """
    bodies = set(classify_body(sc))
    vol = seg_volumes(sc)
    cand = [k for k in vol if k not in bodies and k not in {str(x) for x in exclude}
            and str(k) in sc.get("region_meshes", {})]
    if not cand:
        return None
    bV, bF = mesh_of(sc, sorted(bodies, key=lambda k: -vol[k])[0])
    com = vol_centroid(bV, bF)
    D = bV - com
    dn = np.linalg.norm(D, axis=1)
    cos_cone = math.cos(math.radians(PB_CONE_DEG))

    best, best_reach = None, -np.inf
    for k in cand:
        c = mesh_of(sc, k)[0].mean(0) - com
        d = float(np.linalg.norm(c))
        if d <= 0:
            continue
        u = c / d
        sel = (D @ u) >= cos_cone * np.maximum(dn, 1e-9)          # body verts in this direction
        surface = float(dn[sel].max()) if sel.any() else float(dn.max())
        reach = d / max(surface, 1e-9)
        if reach > best_reach:
            best, best_reach = k, reach
    return best if best_reach >= PB_MIN_REACH else None


def cytoplasm_positions(sc, body_label):
    """{gene: (n,3) µm positions} for molecules labelled as this body. Exact by construction."""
    zs = sc["z_scale"]
    out = {}
    for g, t in sc["transcripts"].items():
        sel = np.asarray(t["s"], int) == int(body_label)
        if not sel.any():
            continue
        out[g] = np.stack([np.asarray(t["x"], float),
                           np.asarray(t["y"], float),
                           np.asarray(t["gz"], float) * zs], axis=1)[sel] * PX
    return out


# ───────────────────────── half-space geometry ─────────────────────────
def half_volume(V, F, normal, origin):
    """Volume of the closed mesh on the (p−origin)·normal > 0 side.

    Putting the origin ON the cutting plane is what makes this exact: the planar cap needed to
    close the clipped surface then lies in a plane through the origin, so every cap triangle has a
    zero scalar triple product and drops out. No polygon ordering, no shapely."""
    n = np.asarray(normal, float)
    n = n / (np.linalg.norm(n) or 1.0)
    P = np.asarray(V, float) - np.asarray(origin, float)
    T = P[np.asarray(F, int)]
    d = T @ n
    inside = d >= 0
    k = inside.sum(1)
    tris = [T[k == 3]]

    def cut(a, b, da, db):
        return a + (da / (da - db))[:, None] * (b - a)

    for want in (1, 2):
        sel = k == want
        if not sel.any():
            continue
        t3, d3 = T[sel], d[sel]
        odd = (d3 >= 0) if want == 1 else (d3 < 0)
        idx = np.argmax(odd, axis=1)
        r = (np.arange(3)[None, :] + idx[:, None]) % 3
        t3 = np.take_along_axis(t3, r[:, :, None], axis=1)
        d3 = np.take_along_axis(d3, r, axis=1)
        v0, v1, v2 = t3[:, 0], t3[:, 1], t3[:, 2]
        d0, d1, d2 = d3[:, 0], d3[:, 1], d3[:, 2]
        p01, p02 = cut(v0, v1, d0, d1), cut(v0, v2, d0, d2)
        if want == 1:
            tris.append(np.stack([v0, p01, p02], 1))
        else:
            tris.append(np.stack([p01, v1, v2], 1))
            tris.append(np.stack([p01, v2, p02], 1))
    A = np.concatenate([t for t in tris if len(t)], 0) if any(len(t) for t in tris) \
        else np.zeros((0, 3, 3))
    if not len(A):
        return 0.0
    return float(np.einsum("ij,ij->i", A[:, 0], np.cross(A[:, 1], A[:, 2])).sum() / 6.0)


def split_volumes(V, F, normal, origin, exact_total=None):
    """(vPos, vNeg), optionally rescaled so they sum to an EXACT total.

    The split is geometric (from the mesh) and the total is exact (from the scene's voxel volume),
    which is the reference's construction."""
    vp = half_volume(V, F, normal, origin)
    vm = half_volume(V, F, -np.asarray(normal, float), origin)
    if exact_total is not None:
        k = exact_total / max(vp + vm, 1e-9)
        vp, vm = vp * k, vm * k
    return vp, vm


def equal_volume_plane(V, F, axis, com, exact_total=None):
    """The equatorial plane: ⟂ `axis`, SHIFTED along it until cytoplasm volume splits 50/50.

    A plane merely through the centre of mass does NOT bisect an irregular cell — the reference
    measures it landing between 0.490 and 0.509 across the 50 zygotes — so the shift is the whole
    point. Returns (normal, origin)."""
    from scipy.optimize import brentq
    a = np.asarray(axis, float)
    a = a / (np.linalg.norm(a) or 1.0)
    total = exact_total if exact_total is not None else mesh_volume(V, F)
    span = float(np.abs((V - com) @ a).max())

    def excess(t):
        vp, _ = split_volumes(V, F, a, com + t * a, exact_total=total)
        return vp - total / 2.0

    lo, hi = -span, span
    if excess(lo) * excess(hi) > 0:                  # degenerate; fall back to the centroid
        return a, np.asarray(com, float)
    t = brentq(excess, lo, hi, xtol=1e-6)
    return a, com + t * a


# ───────────────────────── the measurement ─────────────────────────
def orient_by_count(counts, volPos, volNeg):
    """Rule 3. `counts` is {gene: (nPos, nNeg)}. Returns ({gene: (nFuller, nEmptier)}, vF, vE,
    flipped) with the FULLER half — more cytoplasmic transcripts overall — first."""
    tp = sum(a for a, _ in counts.values())
    tm = sum(b for _, b in counts.values())
    if tp >= tm:
        return dict(counts), volPos, volNeg, False
    return {g: (b, a) for g, (a, b) in counts.items()}, volNeg, volPos, True


def log_folds(counts, vF, vE, eps=EPS):
    """Rules 4 and 6: bulk-corrected log2 concentration fold, fuller over emptier.

    Returns {gene: lfc}. The bulk term is the MEDIAN of the per-gene raw log ratios in this
    embryo, so what survives is asymmetry relative to the rest of the panel."""
    if vF <= 0 or vE <= 0:
        return {}
    raw = {g: math.log2(((a + eps) / vF) / ((b + eps) / vE))
           for g, (a, b) in counts.items() if a + b > 0}
    if not raw:
        return {}
    bulk = float(np.median(list(raw.values())))
    return {g: v - bulk for g, v in raw.items()}


def side_counts(P, com, normals):
    """Side-positive counts of points P for each unit normal. `normals` is (M,3) or (3,)."""
    N = np.atleast_2d(np.asarray(normals, float))
    if not len(P):
        return np.zeros(len(N), np.int64)
    return ((np.asarray(P, float) - np.asarray(com, float)) @ N.T > 0).sum(axis=0).astype(np.int64)


# ───────────────────────── across embryos ─────────────────────────
def bh(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    if not len(p):
        return p
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    run = 1.0
    for i in range(m - 1, -1, -1):
        run = min(run, p[o[i]] * m / (i + 1))
        q[o[i]] = run
    return q


def rank_genes(per_embryo, min_total=MIN_TOTAL, min_embryos=MIN_EMBRYOS):
    """Rule 7. `per_embryo` is a list of {id, gene, lfc, count}. One-sample t-test of each gene's
    per-embryo log2 fold against zero, floors applied, BH added, sorted by P."""
    from scipy import stats
    by = {}
    for r in per_embryo:
        by.setdefault(r["gene"], []).append(r)
    rows = []
    for g, sub in by.items():
        if len(sub) < min_embryos or sum(r["count"] for r in sub) < min_total:
            continue
        v = np.array([r["lfc"] for r in sub], float)
        t, p = stats.ttest_1samp(v, 0.0)
        rows.append({
            "g": g, "m": len(sub), "n": int(sum(r["count"] for r in sub)),
            "meanCount": round(float(np.mean([r["count"] for r in sub])), 1),
            "lfc": round(float(v.mean()), 5), "fold": round(float(2 ** abs(v.mean())), 4),
            "sd": round(float(v.std(ddof=1)), 5), "p": float(p),
            "per": [{"id": r["id"], "lfc": round(r["lfc"], 4), "count": r["count"]}
                    for r in sorted(sub, key=lambda r: -abs(r["lfc"]))],
        })
    if rows:
        for r, q in zip(rows, bh(np.array([r["p"] for r in rows]))):
            r["q"] = float(q)
    rows.sort(key=lambda r: r["p"])
    for i, r in enumerate(rows, start=1):
        r["rank"] = i
        r["side"] = "fuller" if r["lfc"] > 0 else "emptier"
        r["weight"] = round(-math.log10(max(r["p"], 1e-12)), 4)
    return rows


def count_matched_null(per_embryo, vols, seed=20260813, draws=N_NULL_DRAWS, eps=EPS):
    """Rule 5. Per gene, the mean |log2 fold| a gene of the SAME counts would show with no real
    asymmetry — drawn against the volume split, with the same median bulk correction applied.

    `vols` is {embryo_id: (vF, vE)}. Returns {gene: null_fold}."""
    if not per_embryo:
        return {}
    rng = np.random.default_rng(seed)
    cnt = np.array([r["count"] for r in per_embryo])
    vF = np.array([vols[r["id"]][0] for r in per_embryo], float)
    vE = np.array([vols[r["id"]][1] for r in per_embryo], float)
    A = rng.binomial(cnt[:, None], (vF / (vF + vE))[:, None], size=(len(per_embryo), draws))
    B = cnt[:, None] - A
    L = np.log2(((A + eps) / vF[:, None]) / ((B + eps) / vE[:, None]))
    idx = {}
    for i, r in enumerate(per_embryo):
        idx.setdefault(r["id"], []).append(i)
    for _, ii in idx.items():                                   # the same median bulk correction
        L[ii] -= np.median(L[ii], axis=0)[None, :]
    gidx = {}
    for i, r in enumerate(per_embryo):
        gidx.setdefault(r["gene"], []).append(i)
    return {g: float(np.mean(2.0 ** np.abs(L[ii].mean(0)))) for g, ii in gidx.items()}


# ───────────────────────── repo layout ─────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SEG = os.path.join(DATA, "segments")


def scene_path(stage, eid):
    return os.path.join(SEG, f"{stage}__{eid}.json.gz")


def stage_ids(stage):
    import glob
    pre = f"{stage}__"
    return sorted(os.path.basename(p)[len(pre):-len(".json.gz")]
                  for p in glob.glob(os.path.join(SEG, pre + "*.json.gz")))


def probesets():
    p = os.path.join(DATA, "probesets.json")
    return json.load(open(p)) if os.path.isfile(p) else {}
