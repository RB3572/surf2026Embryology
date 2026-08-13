#!/usr/bin/env python3
"""Build data/pseudosperm.json.gz — the Pseudosperm Division Plane project.

This is a port of figure 4.21's method (_work/cells/f4_21.py in the slideshow repo), which is the
specification. Where the two could differ, that file wins.

THE IDEA, IN TWO HALVES.

  A. THE REAL SPERM PLANE. In each of the 30 zygotes with a labelled sperm, cut with the plane
     through {sperm, cytoplasm COM, polar-body COM} and ask, for every gene, whether it AGREES or
     OPPOSES the bulk transcript gradient across that cut. Rank the genes by how reproducibly they
     do so across embryos.

  B. THE PSEUDOSPERM PLANE. In each of the 20 zygotes with no sperm, scan every plane containing
     the COM→polar-body axis and keep the one where those same genes reproduce the same
     agreements and oppositions. That is the pseudosperm plane.

ORIENTATION IS BY TOTAL CYTOPLASMIC COUNT, NOT BY THE SPERM. The sperm lies ON its own plane (it
is one of the three points defining it), so it cannot pick a side. Instead the "+" half is
whichever holds MORE cytoplasmic transcripts across the whole panel — the FULLER half. That is
intrinsic to the zygote, needs no landmark, and is the same quantity on a sperm plane and on a
candidate pseudosperm plane, which is what makes the two comparable at all.

THE BULK CORRECTION IS A MEDIAN OF RATIOS, NOT A RATIO OF TOTALS. A single gene can carry 30% of an
embryo's cytoplasmic transcripts, and abundant genes are less asymmetric than typical, so
log(Σa/Σb) under-corrects everything else. Subtracting the MEDIAN per-gene log ratio (the same idea
as DESeq's size factors) centres the typical gene at zero, so what survives is asymmetry relative
to the rest of the panel.

EVERYTHING IS CYTOPLASM-ONLY. Counts come from the per-transcript segment label in the segments
scenes, so only molecules labelled as the cell body are ever counted — pronuclei and polar body are
excluded exactly, not by a containment test. Volume is the body label's own voxel volume, which by
construction excludes every other label; the two half-volumes come from exact half-space clipping
of the body mesh and are then rescaled to sum to that exact total, so the SPLIT is geometric and
the TOTAL is exact.

WHAT THE ARTIFACT SHIPS. Per-angle cytoplasm half-volumes and per-gene half-counts, for all 180
candidate planes of every zygote. The log-fold, the bulk correction, the template and the fit are
all recomputed in the browser from those, so the p-value cutoff that defines the template is a live
control rather than something baked in here.

Output: data/pseudosperm.json.gz
"""
import glob
import gzip
import json
import math
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SEG = os.path.join(DATA, "segments")
SPERM = os.path.join(DATA, "zygote_sperm.json")
ASSIGN = os.path.join(DATA, "pronuclei_assignments.json")
ZY_MAN = os.path.join(DATA, "zygote_manifest.json")
OUT = os.path.join(DATA, "pseudosperm.json.gz")

VERSION = "pseudosperm-2.0.0"
PX = 0.15                # µm per pixel, in all three axes of the segments scenes' pixel space
K = 180                  # candidate meridional planes per zygote (1° steps over the half-circle)
STEP_DEG = 180.0 / K
MIN_TOTAL = 20           # cytoplasmic transcripts of a gene, summed over embryos, to be ranked
MIN_EMBRYOS = 5          # embryos that must carry the gene at all
ALPHA = 0.05             # the default template cutoff on the raw one-sample P
MIN_ABS_LFC = 0.5        # minimum |mean log2 fold| to join the template
N_FOLD_NULL_DRAWS = 200  # count-matched null draws behind each gene
RNG = np.random.default_rng(20260812)


# ───────────────────────── geometry ─────────────────────────
def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def mesh_of(sc, label):
    r = sc["region_meshes"][str(label)]
    return (np.asarray(r["verts"], float).reshape(-1, 3) * PX,
            np.asarray(r["faces"], int).reshape(-1, 3))


def vol_centroid(V, F):
    """Centroid of the enclosed solid — not the mean vertex, which uneven triangulation biases."""
    T = V[F]
    w = np.einsum("ij,ij->i", T[:, 0], np.cross(T[:, 1], T[:, 2])) / 6.0
    return (((T[:, 0] + T[:, 1] + T[:, 2]) / 4.0) * w[:, None]).sum(0) / w.sum()


def half_volume(V, F, n, origin):
    """Volume of the closed mesh on the (p−origin)·n > 0 side.

    Translating the origin ONTO the cutting plane is what makes this exact: the planar cap needed
    to close the clipped surface then lies in a plane through the origin, so every cap triangle has
    a zero scalar triple product and drops out. No polygon ordering, no shapely."""
    P = V - origin
    s = P @ n
    tri = P[F]
    sg = s[F]
    cnt = (sg > 0).sum(axis=1)
    total = 0.0
    whole = tri[cnt == 3]
    if len(whole):
        total += float(np.einsum("ij,ij->i", np.cross(whole[:, 0], whole[:, 1]), whole[:, 2]).sum())
    for c in (1, 2):
        for t in np.where(cnt == c)[0]:
            p = tri[t]
            q = sg[t]
            poly = []
            for i in range(3):
                j = (i + 1) % 3
                if q[i] > 0:
                    poly.append(p[i])
                if (q[i] > 0) != (q[j] > 0):
                    poly.append(p[i] + (p[j] - p[i]) * (q[i] / (q[i] - q[j])))
            for m in range(1, len(poly) - 1):
                total += float(np.dot(np.cross(poly[0], poly[m]), poly[m + 1]))
    return total / 6.0


def classify_body(sc):
    """The cytoplasm label: the largest segment. A zygote has one body and then a big step down to
    the pronuclei, so this needs no threshold — but it is read off volume, never off label number,
    because label numbering is not consistent across embryos."""
    vol = {str(s["label"]): float(s["volume"]) for s in sc.get("segments", [])}
    if not vol:
        return None, None
    body = max(vol, key=lambda k: vol[k])
    return body, vol


def bh(p):
    """Benjamini–Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    run = 1.0
    for i in range(m - 1, -1, -1):
        run = min(run, p[o[i]] * m / (i + 1))
        q[o[i]] = run
    return q


def delta_encode(a):
    return [int(x) for x in np.diff(np.asarray(a, dtype=np.int64), prepend=0)]


# ───────────────────────── per embryo ─────────────────────────
def sperm_um(sperm_plot, zs_scene, zs_source, eid):
    """The sperm landmark in THIS scene's µm.

    ⚠️ zygote_sperm.json stores z as `frame × z_scale`, and the z_scale it used is the ZYGOTE
    scene's (6.667), while the segments scenes this build reads use 7.0. Multiplying that stored
    number by PX would put the sperm ~5% too shallow — enough to tilt the plane and move tens of
    transcripts across it. So the frame index is recovered and re-scaled into this scene's space,
    which is what the reference does by reading `z_frame` from the raw inventory."""
    frame = sperm_plot[2] / zs_source
    if abs(frame - round(frame)) > 0.02:
        raise ValueError(f"{eid}: sperm z {sperm_plot[2]} is not frame×{zs_source} "
                         f"(got {frame:.3f}) — the storage convention has changed")
    return np.array([sperm_plot[0] * PX, sperm_plot[1] * PX, round(frame) * zs_scene * PX])


def load_embryo(eid, sperm_plot, polar_label, zs_source):
    """Everything one zygote contributes: geometry, the 180-plane sweep, and the sperm plane."""
    p = os.path.join(SEG, f"Zygote__{eid}.json.gz")
    if not os.path.isfile(p):
        return None, "no segments scene"
    sc = json.load(gzip.open(p, "rt"))
    zs = sc["z_scale"]
    body, vols = classify_body(sc)
    if body is None:
        return None, "no segments volume block"
    cyto_vol = vols[body]
    Vb, Fb = mesh_of(sc, body)
    com = vol_centroid(Vb, Fb)

    if polar_label is None or str(polar_label) not in sc["region_meshes"]:
        return None, "no polar body: no meridional family exists"
    pb = mesh_of(sc, polar_label)[0].mean(0)
    axis = unit(pb - com)

    # cytoplasm-only transcripts, exact by the per-molecule segment label
    genes, TX = [], {}
    for g, t in sc["transcripts"].items():
        sel = np.asarray(t["s"], int) == int(body)
        if not sel.any():
            continue
        P = np.stack([np.asarray(t["x"], float),
                      np.asarray(t["y"], float),
                      np.asarray(t["gz"], float) * zs], axis=1)[sel] * PX
        genes.append(g)
        TX[g] = P - com
    if not genes:
        return None, "no cytoplasmic transcripts"

    # the meridional family: planes containing the COM→polar-body axis
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(tmp @ axis) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = unit(np.cross(axis, tmp))
    e2 = unit(np.cross(axis, e1))

    thetas = np.deg2rad(np.arange(K) * STEP_DEG)
    normals = np.cos(thetas)[:, None] * e1[None, :] + np.sin(thetas)[:, None] * e2[None, :]

    volP = np.empty(K)
    for k in range(K):
        vp = half_volume(Vb, Fb, normals[k], com)
        vm = half_volume(Vb, Fb, -normals[k], com)
        volP[k] = cyto_vol * vp / max(vp + vm, 1e-9)          # split geometric, total exact
    volM = cyto_vol - volP

    rows = []
    for g in genes:
        proj = TX[g] @ normals.T                              # (n, K)
        rows.append({"g": g, "n": int(len(TX[g])),
                     "a": delta_encode((proj > 0).sum(axis=0).astype(np.int64))})

    # star-shaped outline of the cytoplasm in the (e1, e2) frame, for the cross-section panel
    along = (Vb - com) @ axis
    slab = Vb[np.abs(along) < 6.0]
    if len(slab) < 20:
        slab = Vb[np.abs(along) < 18.0]
    outline = []
    if len(slab) >= 20:
        dd = slab - com
        uu, vv = dd @ e1, dd @ e2
        ang = np.arctan2(vv, uu)
        rad = np.hypot(uu, vv)
        NB = 120
        b = ((ang + np.pi) / (2 * np.pi) * NB).astype(int) % NB
        for j in range(NB):
            m = b == j
            if not m.any():
                continue
            th = (j + 0.5) / NB * 2 * np.pi - np.pi
            r = float(rad[m].max())
            outline.append([round(r * math.cos(th), 2), round(r * math.sin(th), 2)])

    rec = {
        "id": eid, "scene": f"Zygote__{eid}.json.gz",
        "outline": outline,
        "com_um": [round(float(c), 4) for c in com],
        "pb_um": [round(float(c), 4) for c in pb],
        # 12 dp, not 6: the page rebuilds the cut from this frame to colour transcripts by side,
        # and at 6 dp a molecule sitting on the cut can land on the other half of it
        "axis_um": [round(float(c), 12) for c in axis],
        "e1_um": [round(float(c), 12) for c in e1],
        "e2_um": [round(float(c), 12) for c in e2],
        "z_scale": zs, "body": int(body),
        "cyto_vol": round(float(cyto_vol), 1),
        "volA": [round(float(v), 1) for v in volP],
        "volB": [round(float(v), 1) for v in volM],
        "genes": rows,
        "sperm": None,
    }

    if sperm_plot is not None:
        sp = sperm_um(sperm_plot, zs, zs_source, eid)
        n = np.cross(sp - com, pb - com)
        if np.linalg.norm(n) < 1e-9:
            return rec, "sperm collinear with the polar axis"
        n = unit(n)
        vp = half_volume(Vb, Fb, n, com)
        vm = half_volume(Vb, Fb, -n, com)
        s = cyto_vol / max(vp + vm, 1e-9)
        cnt = {g: int((TX[g] @ n > 0).sum()) for g in genes}
        # the sperm plane is a member of the meridional family, so it has an angle in the (e1, e2)
        # frame — that is what the validation panel compares a refit against
        th_s = math.atan2(float(n @ e2), float(n @ e1)) % math.pi
        rec["sperm"] = {
            "angle_deg": round(math.degrees(th_s), 4),
            "axis_dot": round(float(abs(n @ axis)), 8),   # must be ~0: the plane holds the axis
            "normal_um": [round(float(x), 6) for x in n],
            "volA": round(float(vp * s), 1), "volB": round(float(vm * s), 1),
            "a": cnt, "n": {g: int(len(TX[g])) for g in genes},
            "sperm_um": [round(float(x), 3) for x in sp],
            "dist_to_plane_um": round(float((sp - com) @ n), 4),
        }
    return rec, None


# ───────────────────────── the scoring the browser mirrors ─────────────────────────
def oriented(a, n_of, volA, volB):
    """Orient a split so side F is the FULLER half — more cytoplasmic transcripts overall."""
    tot_a = sum(a.values())
    tot_b = sum(n_of[g] - a[g] for g in a)
    if tot_a >= tot_b:
        return {g: (a[g], n_of[g] - a[g]) for g in a}, volA, volB
    return {g: (n_of[g] - a[g], a[g]) for g in a}, volB, volA


def lfcs(cnt, vF, vE):
    """Bulk-corrected log2 fold of concentration, fuller half over emptier half."""
    if vF <= 0 or vE <= 0:
        return {}
    raw = {}
    for g, (a, b) in cnt.items():
        if a + b <= 0:
            continue
        raw[g] = math.log2(((a + 0.5) / vF) / ((b + 0.5) / vE))
    if not raw:
        return {}
    bulk = float(np.median(list(raw.values())))
    return {g: v - bulk for g, v in raw.items()}


def main():
    for p, what in ((SEG, "data/segments"), (SPERM, "zygote_sperm.json"), (ASSIGN, "assignments")):
        if not os.path.exists(p):
            sys.exit(f"missing {what}: {p}")

    sperm_of = {e["id"]: e["sperm_plot"] for e in json.load(open(SPERM))["embryos"]
                if e.get("sperm_plot")}
    man = {m["id"]: m for m in json.load(open(ZY_MAN))["embryos"]}
    # The ASSIGNMENTS file wins where the two disagree, as the reference does: it is the one a
    # human has looked at (and it carries the manual overrides), whereas the zygote manifest's
    # label comes from the automatic peripheral-segment detector alone.
    polar_of = {m["id"]: m.get("polar_body_label") for m in man.values()}
    for e in json.load(open(ASSIGN))["embryos"]:
        if e.get("polar"):
            polar_of[e["id"]] = e["polar"]["label"]

    zsrc_of = {}
    for f in glob.glob(os.path.join(DATA, "zygote", "*.json.gz")):
        zsrc_of[os.path.basename(f)[:-len(".json.gz")]] = \
            json.load(gzip.open(f, "rt"))["z_scale"]

    ids = sorted(os.path.basename(p)[len("Zygote__"):-len(".json.gz")]
                 for p in glob.glob(os.path.join(SEG, "Zygote__*.json.gz")))
    embryos, skipped = [], []
    for i, eid in enumerate(ids, start=1):
        try:
            rec, why = load_embryo(eid, sperm_of.get(eid), polar_of.get(eid),
                                   zsrc_of.get(eid))
        except Exception as exc:                                  # noqa: BLE001
            rec, why = None, str(exc)
        if rec is None:
            skipped.append({"id": eid, "has_sperm": eid in sperm_of, "reason": why})
            print(f"  -- [{i}/{len(ids)}] {eid}: {why}")
            continue
        m = man.get(eid, {})
        rec["label"] = m.get("label") or eid
        rec["date"] = m.get("date_short", "")
        embryos.append(rec)
        tag = "sperm plane" if rec["sperm"] else "pseudosperm"
        print(f"  [{i}/{len(ids)}] {eid:34s} {len(rec['genes']):3d} genes  {tag}")

    # ---- PART A: the real sperm plane, per gene per embryo ----
    per, geom = [], {}
    for e in embryos:
        sp = e["sperm"]
        if not sp:
            continue
        cnt, vF, vE = oriented(sp["a"], sp["n"], sp["volA"], sp["volB"])
        L = lfcs(cnt, vF, vE)
        geom[e["id"]] = (vF, vE)
        for g, v in L.items():
            a, b = cnt[g]
            per.append({"id": e["id"], "g": g, "lfc": v, "n": a + b, "a": a, "b": b})

    by_gene = {}
    for r in per:
        by_gene.setdefault(r["g"], []).append(r)

    ranking = []
    for g, sub in by_gene.items():
        if len(sub) < MIN_EMBRYOS or sum(r["n"] for r in sub) < MIN_TOTAL:
            continue
        v = np.array([r["lfc"] for r in sub], float)
        t, p = stats.ttest_1samp(v, 0.0)
        ranking.append({
            "g": g, "m": len(sub),
            "n": int(sum(r["n"] for r in sub)),
            "meanCount": round(float(np.mean([r["n"] for r in sub])), 1),
            "lfc": round(float(v.mean()), 5),
            "fold": round(float(2 ** abs(v.mean())), 4),
            "sd": round(float(v.std(ddof=1)), 5),
            "p": float(p),
            "per": [{"id": r["id"], "lfc": round(r["lfc"], 4), "n": r["n"], "a": r["a"]}
                    for r in sorted(sub, key=lambda r: -abs(r["lfc"]))],
        })

    # count-matched null: same counts, same half volumes, no real asymmetry
    if ranking and per:
        cnts = np.array([r["n"] for r in per])
        vF_ = np.array([geom[r["id"]][0] for r in per])
        vE_ = np.array([geom[r["id"]][1] for r in per])
        A = RNG.binomial(cnts[:, None], (vF_ / (vF_ + vE_))[:, None],
                         size=(len(per), N_FOLD_NULL_DRAWS))
        B = cnts[:, None] - A
        LN = np.log2(((A + 0.5) / vF_[:, None]) / ((B + 0.5) / vE_[:, None]))
        idx_of = {}
        for i, r in enumerate(per):
            idx_of.setdefault(r["id"], []).append(i)
        for eid, idx in idx_of.items():                      # the same median bulk correction
            LN[idx] -= np.median(LN[idx], axis=0)[None, :]
        gidx = {}
        for i, r in enumerate(per):
            gidx.setdefault(r["g"], []).append(i)
        for row in ranking:
            L = LN[gidx[row["g"]]]
            row["nullFold"] = round(float(np.mean(2.0 ** np.abs(L.mean(0)))), 4)
            row["excess"] = round(row["fold"] - row["nullFold"], 4)

    q = bh(np.array([r["p"] for r in ranking])) if ranking else np.array([])
    for r, qq in zip(ranking, q):
        r["q"] = float(qq)
    ranking.sort(key=lambda r: r["p"])
    for i, r in enumerate(ranking, start=1):
        r["rank"] = i
        r["side"] = "fuller" if r["lfc"] > 0 else "emptier"
        r["weight"] = round(-math.log10(max(r["p"], 1e-12)), 4)

    n_sperm = sum(1 for e in embryos if e["sperm"])
    n_sig = sum(1 for r in ranking if r["p"] < ALPHA)
    n_tmpl = sum(1 for r in ranking if r["p"] < ALPHA and abs(r["lfc"]) >= MIN_ABS_LFC)

    # the smallest p cutoff at which every sperm-free zygote carries at least one template gene —
    # swept, not guessed, so the page can offer the minimum widening that reaches full coverage
    free_sets = [{r["g"] for r in e["genes"]} for e in embryos if not e["sperm"]]
    p_cover = None
    for pt in sorted({round(x, 5) for x in
                      list(np.arange(0.05, 1.001, 0.005)) + [r["p"] for r in ranking]}):
        tg = {r["g"] for r in ranking if r["p"] < pt and abs(r["lfc"]) >= MIN_ABS_LFC}
        if tg and all(gs & tg for gs in free_sets):
            p_cover = float(pt)
            break

    meta = {
        "version": VERSION,
        "method": "figure 4.21 (_work/cells/f4_21.py) — that file is the specification",
        "n_embryos": len(embryos), "n_sperm": n_sperm, "n_pseudo": len(embryos) - n_sperm,
        "n_skipped": len(skipped), "skipped": skipped,
        "grid": {"n": K, "step_deg": STEP_DEG, "span_deg": 180.0},
        "params": {"MIN_TOTAL": MIN_TOTAL, "MIN_EMBRYOS": MIN_EMBRYOS, "ALPHA": ALPHA,
                   "MIN_ABS_LFC": MIN_ABS_LFC, "N_FOLD_NULL_DRAWS": N_FOLD_NULL_DRAWS},
        "n_genes_ranked": len(ranking),
        "n_significant": n_sig, "n_template_default": n_tmpl,
        "p_full_coverage": p_cover,
        "orientation": "the + half is the FULLER one — more cytoplasmic transcripts across the "
                       "whole panel. The sperm lies ON its own plane, so it cannot pick a side.",
        "bulk": "each gene's log2 fold has that embryo's MEDIAN per-gene log ratio subtracted, so "
                "what is left is asymmetry relative to the rest of the panel",
        "cytoplasm_only": "counts use the per-molecule segment label; volume is the body label's "
                          "own voxel volume. Pronuclei and polar body enter neither.",
        "caveat": "A pseudosperm plane is fitted to a gene template. It is an inference about "
                  "those genes' geometry, not an observation of a sperm.",
    }
    doc = {"meta": meta, "embryos": embryos, "ranking": ranking}
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))

    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e6:.2f} MB)")
    print(f"  {len(embryos)} zygotes: {n_sperm} with a sperm plane, {len(embryos)-n_sperm} pseudosperm"
          f"  ({len(skipped)} skipped)")
    print(f"  ranked {len(ranking)} genes: {n_sig} at raw P < {ALPHA}, "
          f"{n_tmpl} also clear |log2 fold| >= {MIN_ABS_LFC} and form the default template")
    print(f"  full template coverage of the sperm-free zygotes first reached at P < {p_cover}")
    if ranking:
        print("  strongest: " + ", ".join(
            f"{r['g']} (P={r['p']:.1e}, {r['lfc']:+.2f} log2, {r['side']}, m={r['m']})"
            for r in ranking[:5]))


if __name__ == "__main__":
    main()
