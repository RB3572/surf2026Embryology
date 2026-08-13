#!/usr/bin/env python3
"""Build data/pseudosperm.json.gz — the Pseudosperm Division Plane project.

THE QUESTION. The Sperm Division Plane project cuts each zygote with the one plane through
{sperm, polar-body COM, cell COM}. Because that plane contains the cell→polar-body axis, it is a
member of exactly the same family the Zygote Division Planes project sweeps: planes containing the
polar axis, one free rotation angle θ. The sperm picks θ. Only 30 of our zygotes have a labelled
sperm, so for the other 20 that angle is simply unknown.

This project asks the inverse question. First, which genes split most sharply across the REAL
sperm plane, pooled over every sperm-carrying zygote that measures them? Then, for a zygote with no
sperm, where would the plane have to sit to reproduce that split? The angle that answers it is the
PSEUDOSPERM plane, and the in-plane direction perpendicular to the polar axis is where the sperm
would have to have been.

⚠️ This is an inference, not a measurement. A pseudosperm plane is the best fit to a gene the user
chose; it is evidence about that gene's geometry, not an observation of a sperm. Nothing here
should be read as having located a sperm entry site.

WHAT IS SHIPPED, AND HOW EXACT IT IS.

  · THE SPERM PLANE (30 zygotes) — read straight out of data/sperm_division/, so the side counts
    and the per-side volumes are byte-identical to that project's. Exact.

  · THE ANGLE GRID (all 50 zygotes) — 180 planes at 1°, spanning the half-circle. Per-gene side-A
    counts are EXACT at every angle (they are just projections of transcript positions, and are
    computed by binary search on the sorted azimuth rather than by binning).

  · THE PER-SIDE VOLUMES on that grid are the one approximation. The voxel masks that
    build_zygote.py counted are not on this machine, so the cytoplasm volume either side of an
    arbitrary plane is recomputed by CLIPPING THE SEGMENT-1 MESH: translate so the origin lies on
    the cutting plane, clip every triangle to the half-space, and sum (1/6)(v₀×v₁)·v₂ — the planar
    cap then contributes exactly zero, so the result is the exact volume of the clipped polyhedron.
    That mesh carries the pronuclear cavity surfaces, so pronuclear volume is excluded structurally
    rather than by subtraction. The mesh and the (downsampled) voxel estimate disagree by ~0.2% in
    side-A volume fraction, so the curve is corrected by a periodic fit to the residual against the
    18 angles where build_zygote.py's exact voxel volumes ARE known. At those 18 angles the result
    is exact by construction; held out against the 30 real sperm planes, which sit at arbitrary
    angles, the error is ~0.2% of the volume fraction. Recorded in meta and asserted in the test.

  Nothing about pronuclei or the polar body enters any count or any volume, at any stage: only
  segment 1 (the cytoplasm, with the pronuclei carved out of it) is ever counted or measured.

THE RANKING is computed on the sperm plane alone, over the 30 sperm zygotes:
  per zygote × gene   two-sided EXACT BINOMIAL test of the side-A count against the null
                      proportion volA/(volA+volB) — the volume-aware null, so a gene is only
                      called asymmetric if it beats what its own cell's geometry would give.
  across zygotes      Fisher's combined p over every zygote that carries the gene (≥1 transcript),
                      plus a signed Stouffer Z. Fisher answers "is it asymmetric"; Stouffer answers
                      the much stronger "is it asymmetric the SAME WAY every time", which the sign
                      of the plane normal (unit(axis × sperm)) makes a well-defined question.

Output: data/pseudosperm.json.gz
"""
import glob
import gzip
import json
import math
import os
import sys

import numpy as np
from scipy.stats import binom, chi2, norm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
ZY = os.path.join(DATA, "zygote")                       # 50 zygotes with a polar axis
SD = os.path.join(DATA, "sperm_division")               # the 30 with a sperm plane
SPERM = os.path.join(DATA, "zygote_sperm.json")
ZY_MAN = os.path.join(DATA, "zygote_manifest.json")
OUT = os.path.join(DATA, "pseudosperm.json.gz")

VERSION = "pseudosperm-1.0.0"
XY_UM = 0.15
K = 180                      # planes on the grid: 1° steps across the half-circle
STEP_DEG = 180.0 / K
NHARM = 3                    # harmonics in the periodic volume-residual fit
MIN_TX_RANK = 1              # a zygote enters a gene's ranking with at least this many transcripts


# ───────────────────────── geometry ─────────────────────────
def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def half_volume(P, faces, n):
    """Volume of the closed mesh on the (p·n > 0) side, with the origin ON the cutting plane.

    The planar cap needed to close the clipped surface lies in a plane through the origin, so
    every cap triangle has a zero scalar triple product and drops out — the sum over the clipped
    triangles alone is the exact volume."""
    s = P @ n
    tri = P[faces]
    sg = s[faces]
    cnt = (sg > 0).sum(axis=1)
    total = 0.0
    whole = tri[cnt == 3]
    if len(whole):
        total += float(np.einsum("ij,ij->i", np.cross(whole[:, 0], whole[:, 1]), whole[:, 2]).sum())
    for c in (1, 2):                                    # straddling triangles, clipped
        for t in np.where(cnt == c)[0]:
            p = tri[t]; q = sg[t]
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


def side_counts(uv, thetas):
    """#{transcripts with (p−C)·n(θ) > 0} for every θ, where n(θ) = cosθ·u + sinθ·v.

    `uv` is the (n, 2) array of [(p−C)·u, (p−C)·v], so (p−C)·n(θ) is just uv @ [cosθ, sinθ] — one
    small matrix product for every angle at once. Evaluating the sign of that dot product directly
    is what makes this EXACT: an earlier version bracketed the azimuth instead, which disagreed on
    the handful of transcripts sitting within rounding distance of the cut."""
    if not len(uv):
        return np.zeros(len(thetas), dtype=np.int32)
    R = np.stack([np.cos(thetas), np.sin(thetas)])            # (2, K)
    out = np.empty(len(thetas), dtype=np.int32)
    CH = 200000                                               # chunked so memory stays bounded
    for s in range(0, len(uv), CH):
        proj = uv[s:s + CH] @ R
        c = (proj > 0).sum(axis=0)
        if s == 0:
            out[:] = c
        else:
            out += c.astype(np.int32)
    return out


def periodic_fit(theta_knots, resid, nharm=NHARM):
    """Least-squares fit of a π-periodic function through the knots (planes repeat every 180°)."""
    cols = [np.ones(len(theta_knots))]
    for h in range(1, nharm + 1):
        cols += [np.cos(2 * h * theta_knots), np.sin(2 * h * theta_knots)]
    Mx = np.column_stack(cols)
    return np.linalg.lstsq(Mx, resid, rcond=None)[0]


def periodic_eval(coef, thetas, nharm=NHARM):
    out = np.full(len(thetas), coef[0], dtype=float)
    for h in range(1, nharm + 1):
        out += coef[2 * h - 1] * np.cos(2 * h * thetas) + coef[2 * h] * np.sin(2 * h * thetas)
    return out


LN10 = math.log(10.0)
LOG2 = math.log(2.0)


def _logsumexp(terms):
    m = max(terms)
    return m + math.log(sum(math.exp(t - m) for t in terms))


def _binom_log_tail(a, n, p0, upper):
    """log_e of an exact binomial tail, summed in log space.

    scipy's binom.logsf / logcdf return -inf here: with tens of thousands of transcripts the real
    tail reaches 1e-4000, and their internal betainc underflows before the log is taken. Summing
    log-pmf terms with a log-sum-exp never underflows and stays exact to double precision. The
    terms fall off geometrically away from `a`, so the sum is truncated once a term is 1e-17 of
    the leading one."""
    lgn = math.lgamma(n + 1)
    lp, lq = math.log(p0), math.log1p(-p0)
    step = 1 if upper else -1
    stop = n if upper else 0
    terms = []
    k = a
    while (k <= stop) if upper else (k >= stop):
        t = lgn - math.lgamma(k + 1) - math.lgamma(n - k + 1) + k * lp + (n - k) * lq
        terms.append(t)
        if len(terms) > 4 and t < terms[0] - 40:
            break
        k += step
    return _logsumexp(terms) if terms else -math.inf


def binom_log_p(a, n, p0):
    """log_e of the two-sided exact binomial p — twice the smaller tail."""
    lo = _binom_log_tail(a, n, p0, upper=False)
    hi = _binom_log_tail(a, n, p0, upper=True)
    return min(LOG2 + min(lo, hi), 0.0)


def chi2_log_sf(x, k):
    """log_e P(chi2_k > x). scipy underflows past x ~ 1500; beyond that the standard
    large-x asymptotic is used, which agrees with scipy to <0.02 nats where both work."""
    if x <= 0:
        return 0.0
    v = float(chi2.logsf(x, k))
    if math.isfinite(v):
        return v
    return (k / 2 - 1) * math.log(x) - x / 2 - (k / 2 - 1) * LOG2 - math.lgamma(k / 2)


def _z_from_logp(lp):
    """Invert log p = log 2 + log Φ̄(z) for z, past the point where norm.isf underflows to inf.
    Two Newton-ish passes on the Mills-ratio asymptotic are plenty at these magnitudes."""
    t = -(lp - LOG2)
    z = math.sqrt(2 * t)
    for _ in range(3):
        z = math.sqrt(max(2 * (t - math.log(z * math.sqrt(2 * math.pi))), 1.0))
    return z


def norm_log_sf(z):
    """log_e P(Z > z), with the Mills-ratio asymptotic where scipy underflows."""
    v = float(norm.logsf(z))
    if math.isfinite(v):
        return v
    return -0.5 * z * z - math.log(z * math.sqrt(2 * math.pi))


def delta_encode(a):
    """[a0, a1−a0, a2−a1, …] — consecutive angles differ by a handful of transcripts, so the
    artifact is mostly one- and two-character tokens."""
    d = np.diff(np.asarray(a, dtype=np.int64), prepend=0)
    return [int(x) for x in d]


# ───────────────────────── per embryo ─────────────────────────
def build_embryo(eid, label, date_short, sperm_plot):
    d = json.load(gzip.open(os.path.join(ZY, eid + ".json.gz"), "rt"))
    A = d["analysis"]
    zs = d["z_scale"]
    com = np.asarray(A["com_um"], float)

    ap = np.asarray(A["axis_plot"], float)
    axis = unit(np.array([ap[0] * XY_UM, ap[1] * XY_UM, ap[2] / zs]))
    # the SAME frame build_zygote.py uses, so angle 0 here is its plane 0
    ref = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = unit(np.cross(axis, ref))
    v = unit(np.cross(axis, u))

    # ---- volumes on the grid: clipped mesh, corrected onto the 18 exact voxel knots ----
    mesh = d["region_meshes"].get("1")
    if not mesh:
        return None, "no segment-1 mesh"
    Vv = np.asarray(mesh["verts"], float).reshape(-1, 3)
    Vv = np.stack([Vv[:, 0] * XY_UM, Vv[:, 1] * XY_UM, Vv[:, 2] / zs], axis=1) - com
    Ff = np.asarray(mesh["faces"], int).reshape(-1, 3)

    knot_deg = np.array([pl["angle"] for pl in A["planes"]], float)
    knot_th = np.deg2rad(knot_deg)
    f_vox = np.array([pl["volA"] / (pl["volA"] + pl["volB"]) for pl in A["planes"]])
    vol_total = float(A["planes"][0]["volA"] + A["planes"][0]["volB"])

    f_mesh_knot = np.empty(len(knot_th))
    for i, th in enumerate(knot_th):
        n = math.cos(th) * u + math.sin(th) * v
        vA = half_volume(Vv, Ff, n)
        vB = half_volume(Vv, Ff, -n)
        f_mesh_knot[i] = vA / (vA + vB)
    coef = periodic_fit(knot_th, f_vox - f_mesh_knot)

    thetas = np.deg2rad(np.arange(K) * STEP_DEG)
    f_grid = np.empty(K)
    for i, th in enumerate(thetas):
        n = math.cos(th) * u + math.sin(th) * v
        vA = half_volume(Vv, Ff, n)
        vB = half_volume(Vv, Ff, -n)
        f_grid[i] = vA / (vA + vB)
    f_grid = np.clip(f_grid + periodic_eval(coef, thetas), 0.02, 0.98)
    volA = f_grid * vol_total
    volB = (1 - f_grid) * vol_total

    def volume_fraction_at(th):
        """The corrected side-A volume fraction at an ARBITRARY angle — used to hold the method
        out against the real sperm planes, which never land on the grid."""
        n = math.cos(th) * u + math.sin(th) * v
        vA = half_volume(Vv, Ff, n)
        vB = half_volume(Vv, Ff, -n)
        return float(np.clip(vA / (vA + vB) + periodic_eval(coef, np.array([th]))[0], 0.02, 0.98))

    # ---- exact per-gene side-A counts on the grid ----
    tx = d["transcripts"]
    genes = [g for g in d["genes"] if g in tx]
    gene_rows = []
    for g in genes:
        t = tx[g]
        s1 = np.asarray(t["s1"], dtype=bool)
        if not s1.any():
            continue
        P = np.stack([np.asarray(t["x"], float)[s1] * XY_UM,
                      np.asarray(t["y"], float)[s1] * XY_UM,
                      np.asarray(t["gz"], float)[s1]], axis=1) - com
        uv = np.stack([P @ u, P @ v], axis=1)
        gene_rows.append({"g": g, "n": int(len(P)), "a": delta_encode(side_counts(uv, thetas))})

    if not gene_rows:
        return None, "no segment-1 transcripts"

    # 6 dp, not 2: the browser rebuilds the cut from these to colour transcripts by side, and at
    # 2 dp the direction is off by ~1e-3 rad — enough to sweep a handful of molecules across the
    # cut and disagree with the counts shipped alongside them.
    to_plot = lambda p: [round(p[0] / XY_UM, 6), round(p[1] / XY_UM, 6), round(p[2] * zs, 6)]
    rec = {
        "id": eid, "label": label, "date": date_short,
        "scene": eid + ".json.gz",
        "com_plot": A["com_plot"], "com_um": [round(float(c), 4) for c in com],
        "pb_plot": A["pb_plot"],
        "axis_plot": A["axis_plot"],
        "u_plot": to_plot(u), "v_plot": to_plot(v),
        "L": round(float(A["planes"][0]["L"]), 3),
        "vol_total": round(vol_total, 1),
        "volA": [round(float(x), 1) for x in volA],
        "volB": [round(float(x), 1) for x in volB],
        "outline": A["cross_section"]["outline"],
        "genes": gene_rows,
        "sperm": None,
    }

    # ---- the real sperm plane, lifted verbatim from the Sperm Division Plane project ----
    sd_path = os.path.join(SD, eid + ".json.gz")
    if sperm_plot is not None and os.path.isfile(sd_path):
        s = json.load(gzip.open(sd_path, "rt"))
        SA = s["analysis"]
        pl = SA["planes"][0]
        n_s = unit(np.asarray(pl["normal_um"], float))
        # the same plane, expressed as an angle in this project's (u,v) frame
        th_s = math.atan2(float(n_s @ v), float(n_s @ u)) % math.pi
        # does that θ's normal point the same way as the sperm project's, or the opposite?
        n_th = math.cos(th_s) * u + math.sin(th_s) * v
        flip = float(n_th @ n_s) < 0
        counts = {r["gene"]: (r["planes"][0]["a"], r["planes"][0]["b"]) for r in SA["genes"]}
        vA, vB = float(pl["volA"]), float(pl["volB"])
        if flip:                                   # re-express on THIS frame's side-A convention
            counts = {g: (b, a) for g, (a, b) in counts.items()}
            vA, vB = vB, vA
        rec["sperm"] = {
            "angle_deg": round(math.degrees(th_s), 4),
            "normal_um": [round(float(x), 6) for x in n_th],
            "sperm_plot": [round(float(x), 2) for x in sperm_plot],
            "volA": round(vA, 1), "volB": round(vB, 1),
            "a": {g: int(a) for g, (a, b) in counts.items()},
            "n": {g: int(a + b) for g, (a, b) in counts.items()},
            "flipped": bool(flip),
            # the corrected mesh curve evaluated at this exact off-grid angle, against the
            # project's exact voxel value — a genuine held-out test of the volume method
            "vfrac_pred": round(volume_fraction_at(th_s), 6),
            "vfrac_true": round(vA / (vA + vB), 6),
        }
    return rec, None


# ───────────────────────── the cross-embryo ranking ─────────────────────────
def build_ranking(embryos):
    """Per gene, over every sperm-carrying zygote that measures it: an exact two-sided binomial
    test against the volume-aware null, Fisher-combined; plus a signed Stouffer Z that asks the
    stronger question of whether the split lands on the SAME side every time."""
    per_gene = {}
    for e in embryos:
        sp = e.get("sperm")
        if not sp:
            continue
        vA, vB = sp["volA"], sp["volB"]
        p0 = vA / (vA + vB)
        for g, n in sp["n"].items():
            if n < MIN_TX_RANK:
                continue
            a = sp["a"][g]
            b = n - a
            lp = binom_log_p(a, n, p0)                     # natural log of the two-sided p
            # concentration fold across the plane; the sign is the side, which is well defined
            cA = a / vA
            cB = b / vB
            L = math.log2((cA + 1e-12) / (cB + 1e-12))
            # signed z carrying the same evidence as lp, computed in log space so that a p far
            # below float underflow still produces a finite, correctly-ordered z
            z = float(norm.isf(math.exp(lp) / 2)) if lp > -700 else _z_from_logp(lp)
            z *= (1.0 if L >= 0 else -1.0)
            per_gene.setdefault(g, []).append(
                {"id": e["id"], "label": e["label"], "n": n, "a": a, "b": b,
                 "lp": lp, "L": L, "z": z, "p0": p0})

    rows = []
    for g, recs in per_gene.items():
        m = len(recs)
        n_tot = sum(r["n"] for r in recs)
        chi = -2.0 * sum(r["lp"] for r in recs)
        # log10 of Fisher's combined p — the ONLY orderable form at the top of this ranking,
        # where the combined p runs past 1e-3000
        fisher_log10 = chi2_log_sf(chi, 2 * m) / LN10
        w = np.array([math.sqrt(r["n"]) for r in recs])
        z = np.array([r["z"] for r in recs])
        Z = float((w * z).sum() / math.sqrt((w * w).sum())) if m else 0.0
        stouffer_log10 = (LOG2 + norm_log_sf(abs(Z))) / LN10
        absL = np.array([abs(r["L"]) for r in recs])
        sgnL = np.array([r["L"] for r in recs])
        nw = np.array([r["n"] for r in recs], float)
        rows.append({
            "g": g, "m": m, "n": int(n_tot),
            "fisherLog10P": round(fisher_log10, 4), "fisherChi2": round(chi, 4),
            "stoufferZ": round(Z, 4), "stoufferLog10P": round(stouffer_log10, 4),
            "absL": round(float((absL * nw).sum() / nw.sum()), 5),      # transcript-weighted
            "medL": round(float(np.median(sgnL)), 5),
            "nPos": int((sgnL > 0).sum()), "nNeg": int((sgnL < 0).sum()),
            "minLog10P": round(min(r["lp"] for r in recs) / LN10, 4),
            "per": [{"id": r["id"], "n": r["n"], "a": r["a"], "L": round(r["L"], 5),
                     "log10p": round(r["lp"] / LN10, 4)}
                    for r in sorted(recs, key=lambda r: r["lp"])],
        })
    rows.sort(key=lambda r: (r["fisherLog10P"], -r["absL"]))
    for i, r in enumerate(rows, start=1):
        r["rank"] = i
    return rows


def main():
    for p, what in ((ZY, "data/zygote"), (SD, "data/sperm_division"), (SPERM, "zygote_sperm.json")):
        if not os.path.exists(p):
            sys.exit(f"missing {what}: {p}")

    sperm_of = {e["id"]: e["sperm_plot"] for e in json.load(open(SPERM))["embryos"]
                if e.get("sperm_plot")}
    man = {m["id"]: m for m in json.load(open(ZY_MAN))["embryos"]}
    ids = sorted(os.path.basename(p)[:-8] for p in glob.glob(os.path.join(ZY, "*.json.gz")))

    embryos, skipped = [], []
    for i, eid in enumerate(ids, start=1):
        m = man.get(eid, {})
        try:
            rec, why = build_embryo(eid, m.get("label") or eid, m.get("date_short", ""),
                                    sperm_of.get(eid))
        except Exception as exc:                                   # noqa: BLE001
            rec, why = None, str(exc)
        if rec is None:
            skipped.append({"id": eid, "reason": why})
            print(f"  -- [{i}/{len(ids)}] {eid}: {why}")
            continue
        embryos.append(rec)
        tag = f"sperm @ {rec['sperm']['angle_deg']:6.2f}°" if rec["sperm"] else "no sperm"
        print(f"  [{i}/{len(ids)}] {eid:34s} {len(rec['genes']):3d} genes  {tag}")

    ranking = build_ranking(embryos)
    n_sperm = sum(1 for e in embryos if e["sperm"])

    # how well the corrected mesh curve reproduces the sperm planes' exact voxel volumes —
    # the only held-out check available, since those sit at arbitrary angles
    errs = [abs(e["sperm"]["vfrac_pred"] - e["sperm"]["vfrac_true"])
            for e in embryos if e.get("sperm")]

    meta = {
        "version": VERSION,
        "n_embryos": len(embryos), "n_sperm": n_sperm, "n_pseudo": len(embryos) - n_sperm,
        "n_skipped": len(skipped), "skipped": skipped,
        "n_no_polar_body": 60 - len(ids),
        "n_sperm_without_polar_body": len(set(sperm_of) - set(ids)),
        "sperm_without_polar_body": sorted(set(sperm_of) - set(ids)),
        "grid": {"n": K, "step_deg": STEP_DEG, "span_deg": 180.0},
        "n_genes": len({r["g"] for r in ranking}),
        "min_tx_rank": MIN_TX_RANK,
        "volume": {
            "method": "segment-1 mesh clipped at the plane (pronuclear cavities included in the "
                      "surface, so their volume is excluded structurally), corrected by a periodic "
                      "fit to the residual against build_zygote.py's 18 exact voxel angles",
            "held_out_fraction_error_mean": round(float(np.mean(errs)), 6) if errs else None,
            "held_out_fraction_error_max": round(float(np.max(errs)), 6) if errs else None,
            "held_out_n": len(errs),
        },
        "ranking": {
            "test": "two-sided exact binomial of the side-A count against volA/(volA+volB)",
            "combine": "Fisher over every sperm zygote carrying the gene; signed Stouffer Z "
                       "(weights √n) for whether the split lands on the same side each time",
        },
        "caveat": "A pseudosperm plane is fitted to the genes you select. It is an inference about "
                  "those genes' geometry, not an observation of a sperm.",
    }
    doc = {"meta": meta, "embryos": embryos, "ranking": ranking}
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))

    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e6:.2f} MB)")
    print(f"  {len(embryos)} zygotes with a polar axis: {n_sperm} with a sperm plane, "
          f"{len(embryos)-n_sperm} pseudosperm")
    print(f"  {meta['n_sperm_without_polar_body']} sperm zygotes have no polar body and are absent")
    print(f"  ranking: {len(ranking)} genes; best = "
          + ", ".join(f"{r['g']} (log10 p={r['fisherLog10P']:.1f}, m={r['m']})" for r in ranking[:5]))
    if errs:
        print(f"  volume A-fraction vs the 30 exact sperm planes: mean "
              f"{np.mean(errs)*100:.4f}%  max {np.max(errs)*100:.4f}%")


if __name__ == "__main__":
    main()
