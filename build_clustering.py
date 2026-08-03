#!/usr/bin/env python3
"""Build data/clustering.json.gz — genes that occupy the same places inside the zygote.

WHAT THIS ASKS
--------------
Do some genes consistently put their transcripts in the same part of the cell? Not "are they
expressed together" — *where* they sit. Two genes cluster here if, embryo after embryo, they are
enriched in the same shell (cortex ↔ interior) and the same end of the polar-body axis
(animal ↔ vegetal).

WHY A PER-GENE SIGNATURE, AND NOT CO-VARIATION ACROSS EMBRYOS
------------------------------------------------------------
Each embryo is imaged with ONE probeset, and the probesets are disjoint: panel 3 shares zero
genes with panels 0-2. Genes from different probesets are therefore never measured in the same
embryo, so any method built on "these two genes rise and fall together across embryos" can only
ever relate genes inside a single panel. Giving every gene its own signature in a shared,
normalized frame is what lets all 420 genes be compared at once.

THE SIGNATURE
-------------
Per embryo, per gene, over CELL-BODY transcripts only (segment 1 — pronuclei and polar body
excluded, they are separate compartments and would dominate):

  radial  u = |p − COM| / R          → N_RAD equal-VOLUME shells, so a uniform cloud is flat
  axial   v = (p − COM)·â / R        → N_AX equal bins, −1 (vegetal) … +1 (animal, polar body)

R is the 99th percentile of |p − COM| over that embryo's whole cell-body cloud — taken from the
data rather than from a stored radius, so it cannot drift out of step with the segmentation.

Raw fractions would cluster everything into one blob, because every gene inherits the cell's own
shape. So each gene's profile is divided by THAT EMBRYO'S all-transcript profile and logged:

  signature[bin] = log2( shrink(gene_counts, background) / background_frac )

which is enrichment relative to where transcripts generally are in that embryo. Embryo-level
effects (size, orientation, detection efficiency) divide out. Signatures are then averaged over
the embryos where the gene is measured well enough.

CLUSTERING
----------
Correlation distance between mean signatures → Ward hierarchical clustering. Assignments are
emitted for every k in K_RANGE so the page can change k without a rebuild; the default k is the
one with the best silhouette. Two 2-D layouts ship: metric MDS (faithful to the distances) and
t-SNE (separates groups more legibly).

Output: data/clustering.json.gz
"""
import gzip
import json
import os
import sys

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS, TSNE
from sklearn.metrics import silhouette_score

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SCENES = os.path.join(DATA, "planes_all")
MANIFEST = os.path.join(DATA, "planes_all_manifest.json")
OUT = os.path.join(DATA, "clustering.json.gz")

XY_UM = 0.15          # plot px → µm in x/y
N_RAD = 6             # equal-volume radial shells
N_AX = 6              # bins along the polar-body axis
MIN_TX = 30           # cell-body transcripts for an embryo to contribute a profile
MIN_EMB = 5           # embryos a gene needs before it can be clustered
K_RANGE = list(range(4, 13))
ALPHA = 20.0          # pseudo-transcripts of shrinkage (see enrichment())
SEED = 0

# Equal-volume shell edges: u_i = (i/N)^(1/3). A spatially uniform cloud then puts the same
# number of points in every shell, so a flat profile means "no radial preference".
RAD_EDGES = np.array([(i / N_RAD) ** (1.0 / 3.0) for i in range(N_RAD + 1)])
AX_EDGES = np.linspace(-1.0, 1.0, N_AX + 1)


def unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / n if n else v


def counts(P, com, axis, R):
    """(radial, axial) histograms of a point cloud, as raw COUNTS. P is (n,3) in µm."""
    d = P - com
    u = np.linalg.norm(d, axis=1) / R
    v = (d @ axis) / R
    rad = np.histogram(np.clip(u, 0, 1 - 1e-9), bins=RAD_EDGES)[0].astype(float)
    ax = np.histogram(np.clip(v, -1 + 1e-9, 1 - 1e-9), bins=AX_EDGES)[0].astype(float)
    return np.concatenate([rad, ax])


def enrichment(c_gene, f_bg):
    """log2 enrichment of a gene's bin counts over the embryo's own background fractions,
    shrunk toward "no enrichment" in proportion to how little evidence a bin carries.

    A flat pseudocount is NOT good enough here. The axial bins near ±1 are geometrically thin,
    so they legitimately hold few transcripts; an empty one would read as log2(eps/f) ≈ −7 and
    swamp the distance metric with what is really just sampling noise. Adding ALPHA pseudo-
    transcripts distributed like the background (an empirical-Bayes / Dirichlet prior) pulls a
    zero-count bin toward 0 instead, while a well-sampled bin is left essentially untouched.
    """
    # radial and axial are separate compositions over the same transcripts — normalise each half
    out = np.empty_like(c_gene)
    for sl in (slice(0, N_RAD), slice(N_RAD, N_RAD + N_AX)):
        c, f = c_gene[sl], f_bg[sl]
        n = c.sum()
        shrunk = (c + ALPHA * f) / (n + ALPHA)
        out[sl] = np.log2(shrunk / np.maximum(f, 1e-9))
    return out


def load_embryo(eid):
    """Cell-body point clouds per gene, in µm, plus the embryo's frame. None if unusable."""
    p = os.path.join(SCENES, eid + ".json.gz")
    if not os.path.isfile(p):
        return None
    s = json.load(gzip.open(p, "rt"))
    a = s.get("analysis") or {}
    com, axis = a.get("com_um"), a.get("axis_um")
    if com is None or axis is None:
        return None
    com = np.asarray(com, float)
    axis = unit(axis)

    clouds = {}
    allpts = []
    for g, t in (s.get("transcripts") or {}).items():
        x = np.asarray(t["x"], float) * XY_UM
        y = np.asarray(t["y"], float) * XY_UM
        z = np.asarray(t["gz"], float)          # already µm (1 µm z-frames)
        s1 = t.get("s1")
        m = np.asarray(s1, bool) if s1 is not None else np.ones(len(x), bool)
        if not m.any():
            continue
        P = np.stack([x[m], y[m], z[m]], axis=1)
        clouds[g] = P
        allpts.append(P)
    if not allpts:
        return None

    bg = np.concatenate(allpts)
    R = float(np.percentile(np.linalg.norm(bg - com, axis=1), 99))
    if not np.isfinite(R) or R <= 0:
        return None
    return {"com": com, "axis": axis, "R": R, "clouds": clouds, "bg": bg}


def main():
    ids = [e["id"] for e in json.load(open(MANIFEST))["embryos"]]
    print(f"clustering: {len(ids)} zygotes\n")

    # gene -> list of per-embryo signatures; and the embryo showing each gene best
    sigs = {}
    best = {}          # gene -> (n_transcripts, embryo_id)
    n_used = 0
    for eid in ids:
        emb = load_embryo(eid)
        if emb is None:
            print(f"  -- {eid}: unusable"); continue
        n_used += 1
        bg = counts(emb["bg"], emb["com"], emb["axis"], emb["R"])
        # background as fractions, each half (radial / axial) its own composition
        f_bg = np.empty_like(bg)
        for sl in (slice(0, N_RAD), slice(N_RAD, N_RAD + N_AX)):
            f_bg[sl] = bg[sl] / max(bg[sl].sum(), 1.0)
        for g, P in emb["clouds"].items():
            if len(P) < MIN_TX:
                continue
            sigs.setdefault(g, []).append(enrichment(counts(P, emb["com"], emb["axis"], emb["R"]), f_bg))
            if len(P) > best.get(g, (0, None))[0]:
                best[g] = (len(P), eid)

    genes = sorted(g for g, v in sigs.items() if len(v) >= MIN_EMB)
    if len(genes) < 10:
        sys.exit(f"only {len(genes)} genes clear MIN_TX={MIN_TX} / MIN_EMB={MIN_EMB} — nothing to cluster")
    X = np.stack([np.mean(sigs[g], axis=0) for g in genes])
    n_emb = [len(sigs[g]) for g in genes]
    print(f"  embryos used: {n_used}   genes clustered: {len(genes)}"
          f"   (≥{MIN_TX} transcripts in ≥{MIN_EMB} embryos)")
    print(f"  signature: {N_RAD} radial + {N_AX} axial = {X.shape[1]} dims\n")

    # ---- distances + clustering ------------------------------------------------
    D = pdist(X, metric="correlation")
    Z = linkage(D, method="ward")
    Dsq = squareform(D)

    clusters, sil = {}, {}
    for k in K_RANGE:
        lab = fcluster(Z, k, criterion="maxclust") - 1
        clusters[k] = lab
        sil[k] = float(silhouette_score(Dsq, lab, metric="precomputed")) if len(set(lab)) > 1 else -1.0
    best_k = max(K_RANGE, key=lambda k: sil[k])
    print("  silhouette by k: " + "  ".join(f"k={k}:{sil[k]:+.3f}" for k in K_RANGE))
    print(f"  default k = {best_k}\n")

    # ---- 2-D layouts -----------------------------------------------------------
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=SEED,
              normalized_stress="auto", n_init=4).fit_transform(Dsq)
    tsne = TSNE(n_components=2, metric="precomputed", init="random", random_state=SEED,
                perplexity=min(30, max(5, len(genes) // 12))).fit_transform(Dsq)

    def norm(E):
        """Centre and scale a layout to roughly [-1,1] so the page needs no axis logic."""
        E = np.asarray(E, float) - np.asarray(E, float).mean(0)
        s = np.abs(E).max() or 1.0
        return np.round(E / s, 4)

    mds, tsne = norm(mds), norm(tsne)

    # ---- per-cluster summaries -------------------------------------------------
    def cluster_info(lab):
        out = []
        for c in range(int(lab.max()) + 1):
            idx = np.flatnonzero(lab == c)
            centroid = X[idx].mean(0)
            # order members by how central they are — the first is the cluster's exemplar
            order = idx[np.argsort([np.linalg.norm(X[i] - centroid) for i in idx])]
            out.append({
                "n": int(len(idx)),
                "profile": [round(float(v), 4) for v in centroid],
                "members": [genes[i] for i in order],
            })
        return out

    payload = {
        "meta": {
            "version": "clustering-1.0.0",
            "n_embryos": n_used,
            "n_genes": len(genes),
            "n_rad": N_RAD, "n_ax": N_AX,
            "min_tx": MIN_TX, "min_emb": MIN_EMB,
            "k_range": K_RANGE, "default_k": int(best_k),
            "silhouette": {str(k): round(sil[k], 4) for k in K_RANGE},
            "rad_edges": [round(float(v), 4) for v in RAD_EDGES],
            "ax_edges": [round(float(v), 4) for v in AX_EDGES],
        },
        "genes": [
            {
                "gene": g,
                "n_emb": int(n_emb[i]),
                "profile": [round(float(v), 4) for v in X[i]],
                "best_emb": best[g][1], "best_n": int(best[g][0]),
                "mds": [float(mds[i][0]), float(mds[i][1])],
                "tsne": [float(tsne[i][0]), float(tsne[i][1])],
            }
            for i, g in enumerate(genes)
        ],
        # cluster labels + summaries for every k, so the page can change k instantly
        "k": {
            str(k): {
                "labels": [int(v) for v in clusters[k]],
                "clusters": cluster_info(clusters[k]),
            } for k in K_RANGE
        },
    }

    os.makedirs(DATA, exist_ok=True)
    with gzip.open(OUT, "wt") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    kb = os.path.getsize(OUT) / 1024
    sizes = [c["n"] for c in payload["k"][str(best_k)]["clusters"]]
    print(f"  wrote {os.path.relpath(OUT, HERE)} — {len(genes)} genes, k={best_k} "
          f"sizes {sizes} ({kb:.0f} KB)")


if __name__ == "__main__":
    main()
