#!/usr/bin/env python3
"""Build data/halves.json.gz — the two halves of a zygote, four ways of cutting it
(figures 4.14, 4.17, 4.18, 4.19 and 4.15).

THE ONE DECISION EVERYTHING ELSE HANGS OFF: WHICH SIDE IS WHICH.

Cut a zygote in two and you get an A and a B, and nothing yet says which is which. Get that wrong
and every gene's fold averages toward zero across embryos, because half the cells are flipped. Two
obvious landmarks are both unusable here, for reasons worth stating rather than discovering:

  · THE SPERM cannot name a side of the sperm plane, because the plane is drawn THROUGH it — the
    sperm, the cytoplasm centroid and the polar body are the three points that define it. Its
    signed distance to its own plane is 0.00 µm up to rounding, and the sign of that rounding is
    positive for about half the zygotes. The polar body is on the plane for the same reason.
  · THE PRONUCLEI are not independent of the question. The paternal pronucleus descends from the
    sperm and sits near where it entered, so orienting by the maternal→paternal axis orients by
    the sperm at one remove, on a plane the sperm defines.

So the side is named by TOTAL CYTOPLASMIC TRANSCRIPT COUNT: side F is whichever half holds more
transcripts summed over the whole panel. It is intrinsic, needs no landmark, has no free sign, and
is the same quantity on every plane definition — which is what makes the four comparable at all.

  ⚠️ ONE CONSEQUENCE. The bulk correction subtracts each embryo's whole-cytoplasm log ratio, which
  is positive on the fuller half by construction. A gene that merely tracks total transcript
  density therefore reads ≈ 0, and what survives is enrichment BEYOND the bulk. That is the point,
  but it does mean the most abundant genes — the ones that dominate the sum and so fix the
  orientation — sit near zero by construction.

THE FOUR PLANE DEFINITIONS, all through the cytoplasm centroid:

  random       a fixed pseudo-random direction per zygote — the control the others must beat
  equatorial   normal along the COM→polar-body axis (the animal/vegetal cut)
  sperm        through {sperm, COM, polar body} — only the 30 zygotes that have a sperm
  polar18      the best of the meridional family: planes CONTAINING the polar axis, scanned at 1°,
               the one with the largest median per-gene fold asymmetry

  polar18 is chosen to maximise asymmetry, so it is not a fair comparator on its own — it is
  reported next to the random plane precisely so the size of that selection is visible.

WHAT IS COMPUTED

  4.19  half-enrichment volcanoes on the sperm plane and on polar18, per gene: the per-embryo
        bulk-corrected log2(fuller ÷ emptier concentration), a one-sample t-test across embryos,
        floors of ≥20 transcripts and ≥5 embryos. Plus an ALIGNMENT NULL: re-run the whole thing
        with each embryo's side assignment flipped at random, 200 draws, and count how many genes
        get called. That asks whether the calls come from coordinated polarity or from the fact
        that any consistent labelling of halves produces some.
  4.17  per-gene mean fold asymmetry against a COUNT-MATCHED null — a gene with few transcripts
        looks asymmetric by sampling alone, so the null draws the same counts from the same
        volumes and reports the excess over it.
  4.18  the median fold over all testable genes for each of the four plane definitions, with a
        bootstrap interval and a paired test against the random control.
  4.15  do a zygote's two halves sit together? Each half is a vector of per-gene SHARE (a/(a+b)),
        embedded per probeset — the four panels are disjoint gene sets, so one embedding over all
        of them would cluster by panel, not by biology. The pairing statistic is mean within-zygote
        distance ÷ mean all-pairs distance, with P from permuting which A-half meets which B-half.

        ⚠️ AND IT HAS A TRAP IN IT. Under share features the two halves of a zygote are exact
        complements — every feature in B is 1 minus the matching feature in A — so any embedding
        pushes them to OPPOSITE sides and pairing comes out ABOVE 1. That is geometry, not
        biology. The concentration-normalised version of the same halves is the one that answers
        "do halves cluster", and both are built here so the artefact is visible rather than
        reported as a result.

        The reference embeds with UMAP; umap-learn is not installed on this machine, so what is
        shipped is the PCA the reference also ships, and the page says which it is drawing. The
        pairing statistic is computed on the FULL feature vectors, not on the 2-D coordinates, so
        it does not depend on the embedding at all.

Geometry, counts and the meridional family all come from build_pseudosperm.py's own loader, so
this build and the Pseudosperm project cut the same zygotes the same way by construction.

Output: data/halves.json.gz
"""
import collections
import glob
import gzip
import json
import math
import os
import sys

import numpy as np
from scipy import stats

import build_pseudosperm as PS
import embryo_stats as ES

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(DATA, "halves.json.gz")

VERSION = "halves-1.0.0"
MIN_TOTAL = 20            # transcripts of a gene in an embryo for that embryo to count
MIN_EMBRYOS = 5           # embryos a gene needs to be tested at all
CALL_P = 0.05
N_ALIGN_DRAWS = 200       # random side-flip draws behind the alignment null
N_BOOT = 2000             # bootstrap draws for 4.18's intervals
SHARE_PRESENT = 0.90      # a gene must appear in this fraction of a panel's halves to be a feature
RNG = np.random.default_rng(20260813)

PLANES = ["random", "equatorial", "sperm", "polar18"]
PLANE_LABEL = {"random": "Random plane", "equatorial": "Equatorial (polar axis)",
               "sperm": "Sperm plane", "polar18": "Best meridional"}


# ───────────────────────── per-embryo halves ─────────────────────────
def halves_for(rec, err):
    """{plane: {"a": {gene: count}, "n": {gene: total}, "volA", "volB"}} for one zygote.

    Every split is CYTOPLASM-ONLY and every pair of half-volumes is rescaled to the body label's
    own exact voxel volume, exactly as build_pseudosperm does it — the split is geometric, the
    total is exact."""
    if rec is None:
        return None
    K = len(rec["volA"])
    n_of = {}
    a_of = {}
    for row in rec["genes"]:
        n_of[row["g"]] = row["n"]
        a_of[row["g"]] = PS_decode(row["a"])
    out = {}

    # polar18: the meridional plane with the largest median per-gene |log2 fold|
    best_k, best_score = None, -np.inf
    for k in range(K):
        vA, vB = rec["volA"][k], rec["volB"][k]
        if vA <= 0 or vB <= 0:
            continue
        lf = [abs(math.log2(((a_of[g][k] + ES.EPS) / vA) / ((n_of[g] - a_of[g][k] + ES.EPS) / vB)))
              for g in n_of if n_of[g] >= MIN_TOTAL]
        if not lf:
            continue
        s = float(np.median(lf))
        if s > best_score:
            best_k, best_score = k, s
    if best_k is not None:
        out["polar18"] = {"a": {g: int(a_of[g][best_k]) for g in n_of}, "n": dict(n_of),
                          "volA": rec["volA"][best_k], "volB": rec["volB"][best_k],
                          "angle_deg": round(best_k * 180.0 / K, 3)}

    # equatorial: normal along the polar axis. Not in the meridional family (those CONTAIN the
    # axis), so it is cut here directly, on the same mesh and rescaled the same way.
    sc = ES.read_scene(os.path.join(DATA, "segments", rec["scene"]))
    V, F = ES.mesh_of(sc, rec["body"])
    com = np.asarray(rec["com_um"], float)
    axis = np.asarray(rec["axis_um"], float)
    TX = ES.cytoplasm_positions(sc, rec["body"])
    cyto = rec["cyto_vol"]

    def cut(normal, name, extra=None):
        vp = ES.half_volume(V, F, normal, com)
        vm = ES.half_volume(V, F, -normal, com)
        s = cyto / max(vp + vm, 1e-9)
        a = {g: int(((P - com) @ normal > 0).sum()) for g, P in TX.items()}
        r = {"a": a, "n": {g: len(P) for g, P in TX.items()},
             "volA": round(float(vp * s), 1), "volB": round(float(vm * s), 1)}
        if extra:
            r.update(extra)
        out[name] = r

    cut(axis, "equatorial")
    # the random control: a direction that depends only on the embryo id, so a rerun reproduces it
    rng = np.random.default_rng(abs(hash(rec["id"])) % (2 ** 32))
    v = rng.normal(size=3)
    cut(v / np.linalg.norm(v), "random")

    if rec.get("sperm"):
        sp = rec["sperm"]
        out["sperm"] = {"a": {g: int(c) for g, c in sp["a"].items()},
                        "n": {g: int(c) for g, c in sp["n"].items()},
                        "volA": sp["volA"], "volB": sp["volB"],
                        "angle_deg": sp["angle_deg"]}
    return out


def PS_decode(a):
    """build_pseudosperm delta-encodes the per-angle counts; undo it."""
    return np.cumsum(np.asarray(a, np.int64))


# ───────────────────────── the statistics ─────────────────────────
def oriented_lfc(rec):
    """Per-gene bulk-corrected log2(fuller ÷ emptier concentration) for one embryo, one plane.

    The FULLER half is the one with more cytoplasmic transcripts over the whole panel — the only
    orientation rule available that is not the quantity being measured."""
    a, n, vA, vB = rec["a"], rec["n"], rec["volA"], rec["volB"]
    tot_a = sum(a.values())
    tot_b = sum(n[g] - a[g] for g in n)
    flipped = tot_a < tot_b
    if flipped:
        cnt = {g: (n[g] - a[g], a[g]) for g in n}
        vF, vE = vB, vA
    else:
        cnt = {g: (a[g], n[g] - a[g]) for g in n}
        vF, vE = vA, vB
    if vF <= 0 or vE <= 0:
        return None, None, None
    lf = {g: math.log2(((f + ES.EPS) / vF) / ((e + ES.EPS) / vE)) for g, (f, e) in cnt.items()}
    bulk = float(np.median(list(lf.values())))       # median of ratios, not ratio of totals
    return ({g: lf[g] - bulk for g in lf}, cnt,
            {"flipped": bool(flipped), "vF": vF, "vE": vE,
             "totF": max(tot_a, tot_b), "totE": min(tot_a, tot_b), "bulk": round(bulk, 5)})


def test_genes(per):
    """per: gene -> [(embryo, lfc, total)]. One-sample t-test across embryos, BH within the set."""
    rows = []
    for g, sub in per.items():
        keep = [x for x in sub if x[2] >= MIN_TOTAL]
        if len(keep) < MIN_EMBRYOS:
            continue
        y = np.array([x[1] for x in keep], float)
        t, p = stats.ttest_1samp(y, 0.0)
        rows.append({"g": g, "n": len(keep), "total": int(sum(x[2] for x in keep)),
                     "lfc": round(float(y.mean()), 5), "sd": round(float(y.std(ddof=1)), 5),
                     "p": float(p),
                     "per": [{"id": x[0], "lfc": round(x[1], 4), "n": x[2]}
                             for x in sorted(keep, key=lambda x: -abs(x[1]))]})
    if rows:
        for r, q in zip(rows, ES.bh(np.array([r["p"] for r in rows]))):
            r["q"] = float(q)
    rows.sort(key=lambda r: r["p"])
    for i, r in enumerate(rows, start=1):
        r["rank"] = i
        r["called"] = bool(r["p"] < CALL_P)
        r["side"] = "fuller" if r["lfc"] > 0 else "emptier"
    return rows


def alignment_null(by_embryo, n_draws=N_ALIGN_DRAWS):
    """How many genes get called if each embryo's side assignment is flipped at random?

    THIS IS THE QUESTION THE VOLCANO CANNOT ANSWER ON ITS OWN. Any consistent labelling of halves
    produces some calls; what matters is whether the real labelling produces more."""
    embryos = sorted(by_embryo)
    counts = []
    for _ in range(n_draws):
        flip = {e: (1 if RNG.random() < 0.5 else -1) for e in embryos}
        per = collections.defaultdict(list)
        for e in embryos:
            for g, (lfc, tot) in by_embryo[e].items():
                per[g].append((e, lfc * flip[e], tot))
        counts.append(sum(1 for r in test_genes(per) if r["called"]))
    return counts


def count_matched_fold(per_plane_gene, vols):
    """4.17's null: the same counts drawn from the same volumes, so sampling noise alone is
    subtracted rather than argued about.

    ES.count_matched_null wants a FLAT list of {id, gene, lfc, count} — one record per
    embryo x gene — because it applies the same median bulk correction per embryo that the
    observed folds get. Building it any other way would compare a corrected observation against
    an uncorrected null."""
    flat = [{"id": e, "gene": g, "lfc": l, "count": t}
            for g, sub in per_plane_gene.items() for e, l, t in sub]
    return ES.count_matched_null(flat, vols)


# ───────────────────────── 4.15: do the halves pair? ─────────────────────────
def pairing(feats):
    """(within, allpair, ratio, P). feats: {embryo: (vecA, vecB)} over a shared gene set.

    P permutes WHICH A-half is matched to WHICH B-half, so it asks whether a zygote's own two
    halves are closer than an arbitrary pairing — the null the figure actually needs."""
    ids = sorted(feats)
    A = np.array([feats[i][0] for i in ids])
    B = np.array([feats[i][1] for i in ids])
    n = len(ids)
    if n < 4:
        return None
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)     # (i, j): A_i to B_j
    within = float(np.mean(np.diag(D)))
    allpair = float(np.mean(D))
    obs = within / allpair if allpair else None
    # BOTH TAILS, because both are meaningful and only one of them is "clustering". A ratio above
    # 1 is not a null result — it means a zygote's own two halves sit FARTHER apart than an
    # arbitrary pairing, which is exactly what complementary features force.
    lo = hi = 0
    NP = 2000
    for _ in range(NP):
        pm = RNG.permutation(n)
        v = float(np.mean(D[np.arange(n), pm]))
        if v <= within:
            lo += 1
        if v >= within:
            hi += 1
    return {"n": n, "within": round(within, 4), "allpair": round(allpair, 4),
            "ratio": round(obs, 4),
            # P(a random pairing is at least this close) — small means the halves really do
            # cluster. `lo` counts the permutations that got at least as close, so it is the one
            # that goes here; putting `hi` here reads every clustering result as P = 1.
            "p_closer": round((lo + 1) / (NP + 1), 5),
            "p_farther": round((hi + 1) / (NP + 1), 5),     # halves farther → the complementarity
            "p": round(min(1.0, 2 * min((hi + 1) / (NP + 1), (lo + 1) / (NP + 1))), 5)}


def embed(feats):
    """PCA to 2-D on the stacked halves. Deterministic, and the pairing statistic above does not
    use it — it is a picture, not the measurement."""
    ids = sorted(feats)
    X = np.array([v for i in ids for v in feats[i]], float)
    X = X - X.mean(0)
    if X.shape[0] < 3 or X.shape[1] < 2:
        return None
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    Y = U[:, :2] * S[:2]
    ev = (S ** 2) / max(float((S ** 2).sum()), 1e-12)
    return {"ids": ids, "xy": [[round(float(a), 4), round(float(b), 4)] for a, b in Y],
            "explained": [round(float(ev[0]), 4), round(float(ev[1]), 4)]}


def main():
    # the same three inputs build_pseudosperm reads, resolved the same way — the assignments file
    # wins over the manifest for the polar body, because that is the one a human has looked at
    sperm_of = {e["id"]: e["sperm_plot"] for e in json.load(open(PS.SPERM))["embryos"]
                if e.get("sperm_plot")}
    man = {m["id"]: m for m in json.load(open(PS.ZY_MAN))["embryos"]}
    polar_of = {m["id"]: m.get("polar_body_label") for m in man.values()}
    for e in json.load(open(PS.ASSIGN))["embryos"]:
        if e.get("polar"):
            polar_of[e["id"]] = e["polar"]["label"]
    zsrc_of = {}
    for f in glob.glob(os.path.join(DATA, "zygote", "*.json.gz")):
        zsrc_of[os.path.basename(f)[:-len(".json.gz")]] = json.load(gzip.open(f, "rt"))["z_scale"]
    probeset = ES.probesets()

    cache = os.environ.get("HALVES_CACHE")
    stage1 = None
    if cache and os.path.isfile(cache):
        print(f"  reading the stage-1 cache: {cache}")
        stage1 = json.load(gzip.open(cache, "rt"))

    ids = ES.stage_ids("Zygote")
    per_plane = {p: collections.defaultdict(list) for p in PLANES}   # plane -> gene -> [(id,lfc,n)]
    by_embryo = {p: {} for p in PLANES}                              # plane -> id -> {g:(lfc,tot)}
    vols = {p: {} for p in PLANES}
    orient = {p: [] for p in PLANES}
    shares = {}                                                      # id -> {gene: share} sperm/polar18
    emb_meta, skipped = [], []

    cache_out = {}
    for i, eid in enumerate(ids, start=1):
        if stage1 is not None:
            H = stage1["halves"].get(eid)
            if H is None:
                skipped.append({"id": eid, "reason": stage1["skipped"].get(eid, "cached skip")})
                continue
        else:
            try:
                rec, err = PS.load_embryo(eid, sperm_of.get(eid), polar_of.get(eid),
                                          zsrc_of.get(eid))
            except Exception as exc:                              # noqa: BLE001
                skipped.append({"id": eid, "reason": str(exc)}); continue
            if rec is None:
                skipped.append({"id": eid, "reason": err or "no geometry"})
                continue
            H = halves_for(rec, err)
            if H:
                cache_out[eid] = H
        if not H:
            skipped.append({"id": eid, "reason": "no usable plane"})
            continue
        got = []
        for p in PLANES:
            if p not in H:
                continue
            lf, cnt, meta = oriented_lfc(H[p])
            if lf is None:
                continue
            got.append(p)
            for g, v in lf.items():
                per_plane[p][g].append((eid, v, H[p]["n"][g]))
            by_embryo[p][eid] = {g: (lf[g], H[p]["n"][g]) for g in lf}
            vols[p][eid] = (meta["vF"], meta["vE"])
            orient[p].append({"id": eid, "flipped": meta["flipped"],
                              "totF": meta["totF"], "totE": meta["totE"],
                              "vF": round(meta["vF"], 1), "vE": round(meta["vE"], 1),
                              "bulk": meta["bulk"],
                              "frac": round(meta["totF"] / max(meta["totF"] + meta["totE"], 1), 4)})
            if p == ("sperm" if "sperm" in H else "polar18"):
                shares[eid] = {"plane": p,
                               "share": {g: c[0] / (c[0] + c[1])
                                         for g, c in cnt.items() if c[0] + c[1] > 0},
                               "cnt": {g: [c[0], c[1]] for g, c in cnt.items() if c[0] + c[1] > 0}}
        emb_meta.append({"id": eid, "label": man.get(eid, {}).get("label") or eid,
                         "probeset": probeset.get(eid, "?"), "planes": got,
                         "has_sperm": "sperm" in H,
                         "polar18_deg": H.get("polar18", {}).get("angle_deg"),
                         "sperm_deg": H.get("sperm", {}).get("angle_deg")})
        print(f"  [{i}/{len(ids)}] {eid:32s} {' '.join(got)}")
    if cache and stage1 is None:
        with gzip.open(cache, "wt") as fh:
            json.dump({"halves": cache_out,
                       "skipped": {s_["id"]: s_["reason"] for s_ in skipped}}, fh)
        print(f"  wrote the stage-1 cache: {cache}")

    # ---- 4.19 volcanoes + the alignment null ----
    volcano = {}
    for p in ("sperm", "polar18"):
        rows = test_genes(per_plane[p])
        n_called = sum(1 for r in rows if r["called"])
        null = alignment_null(by_embryo[p])
        nz = np.array(null, float)
        volcano[p] = {
            "label": PLANE_LABEL[p], "genes": rows, "n_called": n_called,
            "n_embryos": len(by_embryo[p]), "orientation": orient[p],
            "n_flipped": sum(1 for o in orient[p] if o["flipped"]),
            "null": {"median": float(np.median(nz)), "p95": float(np.percentile(nz, 95)),
                     "max": int(nz.max()), "draws": len(null),
                     "p": float((np.sum(nz >= n_called) + 1) / (len(null) + 1)),
                     "hist": [int(x) for x in np.bincount(np.asarray(null, int),
                                                          minlength=max(int(nz.max()), n_called) + 1)]},
        }
        print(f"  {p}: {len(rows)} genes, {n_called} called, "
              f"alignment null median {np.median(nz):.0f} P={volcano[p]['null']['p']:.4f}")

    # ---- 4.17 per-gene fold vs a count-matched null (on the best meridional plane) ----
    fold_plane = "polar18"
    pe = {g: [(e, l, t) for e, l, t in sub] for g, sub in per_plane[fold_plane].items()}
    per_gene = []
    for r in test_genes(pe):
        vals = [abs(x["lfc"]) for x in r["per"]]
        obs = float(np.mean([2 ** v for v in vals]))                 # fold, not log fold
        per_gene.append({"g": r["g"], "n": r["n"], "total": r["total"],
                         "fold": round(obs, 5), "p": r["p"], "q": r.get("q"),
                         "lfc": r["lfc"]})
    nullfold = count_matched_fold(per_plane[fold_plane], vols[fold_plane])
    for r in per_gene:
        nf = nullfold.get(r["g"])
        r["null"] = round(float(nf), 5) if nf is not None else None
        r["excess"] = round(r["fold"] - r["null"], 5) if r["null"] is not None else None
    per_gene.sort(key=lambda r: -(r["excess"] if r["excess"] is not None else -9))

    # ---- 4.18 median fold by plane definition, bootstrap + paired vs random ----
    fold_by_gene = {}
    for p in PLANES:
        d = {}
        for r in test_genes(per_plane[p]):
            d[r["g"]] = float(np.mean([2 ** abs(x["lfc"]) for x in r["per"]]))
        fold_by_gene[p] = d
    shared = sorted(set(fold_by_gene["random"]) & set(fold_by_gene["equatorial"]) &
                    set(fold_by_gene["polar18"]))
    byplane = []
    for p in PLANES:
        d = fold_by_gene[p]
        gs = [g for g in shared if g in d] if p != "sperm" else sorted(set(shared) & set(d))
        v = np.array([d[g] for g in gs], float)
        if not len(v):
            continue
        boot = np.array([np.median(RNG.choice(v, len(v), replace=True)) for _ in range(N_BOOT)])
        row = {"plane": p, "label": PLANE_LABEL[p], "n_genes": len(gs),
               "median": round(float(np.median(v)), 5),
               "ci_lo": round(float(np.percentile(boot, 2.5)), 5),
               "ci_hi": round(float(np.percentile(boot, 97.5)), 5)}
        if p != "random":
            pair = [(fold_by_gene["random"][g], d[g]) for g in gs if g in fold_by_gene["random"]]
            if len(pair) >= 5:
                a = np.array([x[0] for x in pair]); b = np.array([x[1] for x in pair])
                try:
                    st = stats.wilcoxon(b - a, alternative="two-sided")
                    row["p_vs_random"] = float(st.pvalue)
                except ValueError:
                    row["p_vs_random"] = None
                row["n_paired"] = len(pair)
                row["median_gain"] = round(float(np.median(b - a)), 5)
        byplane.append(row)
    heat_genes = [r["g"] for r in sorted(per_gene, key=lambda r: -(r["excess"] or -9))[:15]]
    heat = [{"g": g, **{p: (round(fold_by_gene[p][g], 4) if g in fold_by_gene[p] else None)
                        for p in PLANES}} for g in heat_genes]

    # ---- 4.15 pairing, per probeset, both normalisations ----
    pair_out = {}
    for norm in ("ratio", "conc"):
        panels = {}
        by_panel = collections.defaultdict(list)
        for eid, s in shares.items():
            by_panel[str(probeset.get(eid, "?"))].append(eid)
        for panel, members in sorted(by_panel.items()):
            if len(members) < 4:
                continue
            counts = collections.Counter()
            for eid in members:
                counts.update(shares[eid]["share"].keys())
            genes = sorted(g for g, c in counts.items() if c >= SHARE_PRESENT * len(members))
            if len(genes) < 5:
                continue
            feats = {}
            for eid in members:
                sh = shares[eid]["share"]
                a = np.array([sh.get(g, 0.5) for g in genes], float)
                if norm == "ratio":
                    feats[eid] = (a, 1.0 - a)          # exact complements — see the docstring
                else:
                    # CONCENTRATION, and deliberately NOT renormalised to sum to 1. Rescaling each
                    # half to a unit sum would put the compositional complementarity straight back
                    # in and there would be no contrast with the share panels at all. log1p of
                    # transcripts per µm³ leaves each half carrying its embryo's own abundance
                    # profile, which is the thing that could make two halves of one zygote cluster.
                    pl = shares[eid]["plane"]
                    vF, vE = vols[pl][eid]
                    c = shares[eid]["cnt"]
                    x = np.array([c.get(g, [0, 0])[0] for g in genes], float) / vF
                    y = np.array([c.get(g, [0, 0])[1] for g in genes], float) / vE
                    feats[eid] = (np.log1p(x * 1e4), np.log1p(y * 1e4))
            st = pairing(feats)
            em = embed(feats)
            if st and em:
                panels[panel] = {"stat": st, "embed": em, "n_genes": len(genes),
                                 "members": sorted(members)}
        pair_out[norm] = panels

    doc = {
        "meta": {
            "version": VERSION,
            "method": "figures 4.14, 4.15, 4.17, 4.18 and 4.19",
            "params": {"MIN_TOTAL": MIN_TOTAL, "MIN_EMBRYOS": MIN_EMBRYOS, "CALL_P": CALL_P,
                       "N_ALIGN_DRAWS": N_ALIGN_DRAWS, "N_BOOT": N_BOOT,
                       "SHARE_PRESENT": SHARE_PRESENT},
            "orientation": "side F is the half holding MORE cytoplasmic transcripts over the whole "
                           "panel. The sperm cannot name a side of the plane drawn THROUGH it, and "
                           "the pronuclei are not independent of the sperm, so total count is the "
                           "only rule left that is intrinsic and has no free sign.",
            "bulk": "each embryo's MEDIAN per-gene log ratio is subtracted, not the ratio of "
                    "totals — a single gene can carry 30% of the panel. A gene that merely tracks "
                    "total density therefore reads about 0, by construction.",
            "planes": {p: PLANE_LABEL[p] for p in PLANES},
            "polar18": "the best of the meridional family — planes CONTAINING the polar axis, "
                       "scanned at 1°, keeping the largest median per-gene fold. It is SELECTED to "
                       "be asymmetric, which is why it is always reported next to the random "
                       "plane.",
            "alignment_null": "the volcano is re-run with each embryo's side assignment flipped at "
                              f"random, {N_ALIGN_DRAWS} times. Any consistent labelling of halves "
                              "produces some calls; this asks whether the real one produces more.",
            "pairing_trap": "under SHARE features a zygote's two halves are exact complements — "
                            "every feature in B is 1 minus the one in A — so any embedding pushes "
                            "them apart and the pairing ratio comes out ABOVE 1. That is geometry. "
                            "The concentration-normalised panels are the ones that answer whether "
                            "halves cluster; both are shipped so the artefact is visible.",
            "embedding": "PCA, not UMAP: umap-learn is not installed on this machine. The pairing "
                         "statistic is computed on the FULL feature vectors, never on the 2-D "
                         "coordinates, so it does not depend on the embedding at all.",
            "n_embryos": len(emb_meta),
            "n_sperm": sum(1 for e in emb_meta if e["has_sperm"]),
            "skipped": skipped,
        },
        "embryos": emb_meta,
        "volcano": volcano,
        "perGene": {"plane": fold_plane, "rows": per_gene},
        "byPlane": byplane,
        "heat": {"planes": PLANES, "rows": heat},
        "pairing": pair_out,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e6:.2f} MB)")
    print(f"  {len(emb_meta)} zygotes ({doc['meta']['n_sperm']} with sperm)")
    for r in byplane:
        extra = (f"  vs random {r['median_gain']:+.3f} P={r.get('p_vs_random'):.2g}"
                 if r.get("p_vs_random") is not None else "")
        print(f"  {r['plane']:11s} median fold {r['median']:.4f} "
              f"[{r['ci_lo']:.4f}, {r['ci_hi']:.4f}] over {r['n_genes']} genes{extra}")
    for norm, panels in pair_out.items():
        print(f"  pairing/{norm}: " + ", ".join(
            f"panel {k} {v['stat']['ratio']:.2f} "
            f"(closer P {v['stat']['p_closer']:.3f}, n {v['stat']['n']})"
            for k, v in sorted(panels.items())))


if __name__ == "__main__":
    main()
