#!/usr/bin/env python3
"""
Sperm alpha/beta analysis for 2-cell embryos.

Recreates the mentor's "alpha/beta assignments across methods" grid (one row per
assignment method, one column per embryo, each cell = which ORIGINAL blastomere that
method calls alpha: A = segmentation labels 1+3, B = labels 2+4).

Two fully-defined methods:
  * "Sperm entry"            — the blastomere the sperm entered is alpha (seg 1/3 -> A, 2/4 -> B).
  * "Higher total transcript"— the blastomere with more transcripts is alpha (tie -> T).
    (verified against Harry's late grids: 20/20.)

The remaining rows recreate Harry's methods.  His deck names the methods and shows the
results but never the exact parameters (gene-panel membership, the per-gene t-statistic,
the PCA axes), so these are BEST-GUESS reconstructions of the most likely implementation:

  * "Expression-axis PCA"        — project each blastomere's per-gene volume-normalized
                                    log-ratio vector onto PC1 of the stage-wide asymmetry
                                    covariance over ALL genes; alpha = +side.
  * "Decreased/Increased panel"  — genes split into maternal (fraction falls zygote->2-cell)
                                    vs zygotic (fraction rises), from the stage-wise
                                    counts_summary matrices.
  * "... : ratio-sum"            — sum the volume-normalized per-gene log-ratio q_g over the
                                    panel; alpha = sign.
  * "... : PCA"                  — project that same q vector onto the panel's leading
                                    covariance axis (genes weighted by the dominant covarying
                                    asymmetry pattern instead of uniformly).
  * "Exhaustive: unfiltered"     — per-gene Poisson deviation t_g of A's count from its
                                    volume-expected share; exhaustive search over the top-k
                                    |t| cutoff (max |vote|/sqrt(k)); alpha = majority sign.
  * "Exhaustive: mean count>=20" — same, restricted to genes with mean blastomere count >= 20.

Per-blastomere transcript counts come from the already-built segments data
(public/data/segments/<Stage>__<id>.json.gz, which carries each transcript's segment).
Sperm blastomere comes from data/merfish_sperm.csv (the labelled sperm's `segment`).
Gene panels come from data/../TranscriptomicsData/*/*_{zyg,e2c,l2c}_counts_summary.csv.

Output: public/data/alphabeta.json   {stages:{early:[...], late:[...]}, methods:[...], panels:{...}}
"""
import csv
import glob
import gzip
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SEG_DIR = os.path.join(HERE, "public", "data", "segments")
SPERM_CSV = os.path.join(HERE, "..", "data", "merfish_sperm.csv")
TX_DIR = os.path.join(HERE, "..", "TranscriptomicsData")
OUT = os.path.join(HERE, "public", "data", "alphabeta.json")

# blastomere A = segmentation labels {1,3}; blastomere B = {2,4}. label 5 = extra body (ignored).
A_LABELS, B_LABELS = {1, 3}, {2, 4}
STAGES = [("early", "Early2Cell"), ("late", "Late2Cell")]

# methods produced here (order = row order in the grid). Each has an id + a display label.
# order mirrors Harry's grid, with our Sperm-entry row on top.
METHODS = [
    {"id": "sperm", "label": "Sperm entry", "guess": False},
    {"id": "total", "label": "Higher total transcript", "guess": False},
    {"id": "pca_expr", "label": "Expression-axis PCA", "guess": True},
    {"id": "pca_dec", "label": "Decreased panel: PCA", "guess": True},
    {"id": "ratio_dec", "label": "Decreased panel: ratio-sum", "guess": True},
    {"id": "pca_inc", "label": "Increased panel: PCA", "guess": True},
    {"id": "ratio_inc", "label": "Increased panel: ratio-sum", "guess": True},
    {"id": "exh_unf", "label": "Exhaustive: unfiltered", "guess": True},
    {"id": "exh_20", "label": "Exhaustive: mean count ≥ 20", "guess": True},
]

EPS = 1e-9


def blastomere_of_segment(seg):
    if seg in A_LABELS:
        return "A"
    if seg in B_LABELS:
        return "B"
    return None


def sperm_blastomere_by_embryo():
    """embryo_id -> {'blast','seg'} for each sperm-positive 2-cell embryo."""
    out = {}
    if not os.path.exists(SPERM_CSV):
        return out
    for r in csv.DictReader(open(SPERM_CSV, newline="")):
        eid = (r.get("resolved_embryo_id") or "").strip()
        md = r.get("mefish_dir", "")
        if not eid or not ("e2c" in md or "l2c" in md):
            continue
        if not r.get("x_px", "").strip() or r.get("no_sperm", "").strip().lower() == "yes":
            continue
        try:
            seg = int(float(r.get("segment", "")))
        except (TypeError, ValueError):
            continue
        bl = blastomere_of_segment(seg)
        if bl:
            out[eid] = {"blast": bl, "seg": seg}
    return out


# ---------- gene panels (maternal / decreased vs zygotic / increased) ----------
def _pool_fraction(suffixes):
    """gene -> fraction of total transcripts, pooled across all counts_summary of these stages."""
    tot, grand = {}, 0
    for suf in suffixes:
        for f in glob.glob(os.path.join(TX_DIR, "*", f"*_{suf}_counts_summary.csv")):
            with open(f, newline="") as fh:
                rd = csv.reader(fh)
                next(rd, None)
                for row in rd:
                    if not row:
                        continue
                    g = row[0]
                    s = sum(int(x) for x in row[1:] if x.strip().lstrip("-").isdigit())
                    tot[g] = tot.get(g, 0) + s
                    grand += s
    return {g: (c / grand if grand else 0.0) for g, c in tot.items()}, tot


def build_panels():
    """Split genes into 'decreased' (maternal: fraction falls zygote->2-cell) and
    'increased' (zygotic: fraction rises).  Returns (dec:set, inc:set)."""
    fz, _ = _pool_fraction(["zyg"])
    ftc, tc = _pool_fraction(["e2c", "l2c"])
    genes = set(fz) | set(ftc)
    dec, inc = set(), set()
    for g in genes:
        # ignore ultra-low genes that never appear meaningfully in 2-cell
        if tc.get(g, 0) < 1:
            dec.add(g)
            continue
        (inc if ftc.get(g, 0.0) >= fz.get(g, 0.0) else dec).add(g)
    return dec, inc


# ---------- per-embryo blastomere expression ----------
def embryo_expression(seg_file):
    """Return per-gene (a,b) counts, blastomere volumes (VA,VB), mask labels."""
    d = json.load(gzip.open(seg_file, "rt"))
    vol = {s.get("label"): float(s.get("volume", 0) or 0) for s in d.get("segments", [])}
    VA = sum(vol.get(l, 0.0) for l in A_LABELS) or 1.0
    VB = sum(vol.get(l, 0.0) for l in B_LABELS) or 1.0
    genes = {}
    for g, t in d.get("transcripts", {}).items():
        a = b = 0
        for s in t.get("s", []):
            if s in A_LABELS:
                a += 1
            elif s in B_LABELS:
                b += 1
        if a or b:
            genes[g] = (a, b)
    return genes, VA, VB, d.get("mask_labels", [])


def logratio(a, b, VA, VB):
    """Volume-normalized per-gene log2 density ratio; 0 == symmetric split."""
    return math.log2((a + 0.5) / (b + 0.5)) - math.log2(VA / VB)


def sign_call(x, tie=1e-9):
    if x > tie:
        return "A"
    if x < -tie:
        return "B"
    return "T"


def pc1(matrix):
    """Leading right singular vector (gene weights) of a centered embryos x genes matrix."""
    if matrix.shape[0] < 2 or matrix.shape[1] < 1:
        return None
    M = matrix - matrix.mean(axis=0, keepdims=True)
    if not np.any(np.abs(M) > EPS):
        return None
    try:
        _, _, Vt = np.linalg.svd(M, full_matrices=False)
    except np.linalg.LinAlgError:
        return None
    v = Vt[0]
    # deterministic orientation
    if v.sum() < 0:
        v = -v
    return v


def exhaustive_call(genes, VA, VB, min_mean=0.0):
    """Per-gene Poisson deviation of A's count from its volume-expected share; exhaustive
    top-k |t| cutoff search maximizing |vote|/sqrt(k); alpha = majority sign of that subset."""
    ts = []
    for a, b in genes.values():
        n = a + b
        if n == 0 or (a + b) / 2.0 < min_mean:
            continue
        eA = n * VA / (VA + VB)                 # volume-expected count in A under symmetry
        t = (a - eA) / math.sqrt(eA + 1.0)      # Poisson-ish z of A's excess
        ts.append(t)
    if not ts:
        return "NA", 0.0, 0
    ts.sort(key=lambda t: -abs(t))
    best_vote, best_score, best_k = 0, -1.0, 0
    run = 0
    for k, t in enumerate(ts, 1):
        run += 1 if t > 0 else (-1 if t < 0 else 0)
        score = abs(run) / math.sqrt(k)
        if score > best_score:
            best_score, best_vote, best_k = score, run, k
    return sign_call(best_vote), float(best_vote), best_k


def main():
    sperm = sperm_blastomere_by_embryo()
    dec_panel, inc_panel = build_panels()
    print(f"  panels: decreased(maternal)={len(dec_panel)}  increased(zygotic)={len(inc_panel)}")

    stages_out = {}
    n_total = 0
    for key, stage_dir in STAGES:
        # 1) load every valid embryo (needs both blastomeres)
        embryos = []
        for f in sorted(glob.glob(os.path.join(SEG_DIR, f"{stage_dir}__*.json.gz"))):
            eid = os.path.basename(f).split("__", 1)[1][: -len(".json.gz")]
            genes, VA, VB, mask = embryo_expression(f)
            if not ({1, 2, 3, 4} <= set(mask)):
                continue
            q = {g: logratio(a, b, VA, VB) for g, (a, b) in genes.items()}
            embryos.append({"id": eid, "genes": genes, "VA": VA, "VB": VB, "q": q})

        # 2) PCA axes (all genes / decreased / increased) from the stage-wide asymmetry matrix
        universe = sorted({g for e in embryos for g in e["q"]})
        idx = {g: i for i, g in enumerate(universe)}
        R = np.zeros((len(embryos), len(universe)))
        for ei, e in enumerate(embryos):
            for g, v in e["q"].items():
                R[ei, idx[g]] = v
        sub = {
            "pca_expr": np.array([True] * len(universe)),
            "pca_dec": np.array([g in dec_panel for g in universe]),
            "pca_inc": np.array([g in inc_panel for g in universe]),
        }
        axes = {}
        for mid, msk in sub.items():
            cols = np.where(msk)[0]
            axes[mid] = (cols, pc1(R[:, cols])) if len(cols) else (cols, None)

        # 3) per-embryo method calls
        rows = []
        for ei, e in enumerate(embryos):
            genes, VA, VB, q = e["genes"], e["VA"], e["VB"], e["q"]
            nA = sum(a for a, _ in genes.values())
            nB = sum(b for _, b in genes.values())
            calls, scores = {}, {}

            calls["total"] = "A" if nA > nB else ("B" if nB > nA else "T")
            sp = sperm.get(e["id"])
            calls["sperm"] = sp["blast"] if sp else "NA"

            for mid, panel in (("ratio_dec", dec_panel), ("ratio_inc", inc_panel)):
                vals = [q[g] for g in q if g in panel]
                if not vals:
                    calls[mid], scores[mid] = "NA", None
                else:
                    s = float(np.sum(vals))
                    calls[mid], scores[mid] = sign_call(s), round(s, 2)

            for mid in ("pca_expr", "pca_dec", "pca_inc"):
                cols, v = axes[mid]
                if v is None or len(cols) == 0:
                    calls[mid], scores[mid] = "NA", None
                else:
                    s = float(R[ei, cols] @ v)
                    calls[mid], scores[mid] = sign_call(s), round(s, 3)

            for mid, mm in (("exh_unf", 0.0), ("exh_20", 20.0)):
                c, vote, k = exhaustive_call(genes, VA, VB, mm)
                calls[mid] = c
                scores[mid] = None if c == "NA" else f"{int(vote):+d}/{k}"

            rows.append({
                "id": e["id"],
                "calls": calls,
                "raw": {"nA": nA, "nB": nB, "spermSeg": (sp["seg"] if sp else None),
                        "scores": scores},
            })

        rows.sort(key=lambda r: r["id"])
        stages_out[key] = rows
        n_total += len(rows)
        n_sp = sum(1 for r in rows if r["calls"]["sperm"] != "NA")
        print(f"  {key:5s} 2-cell: {len(rows)} embryos ({n_sp} with a labelled sperm)")

    payload = {
        "methods": METHODS,
        "legend": {"A": "alpha = original A (labels 1+3)",
                   "B": "alpha = original B (labels 2+4)",
                   "T": "tie", "NA": "unavailable"},
        "panels": {"decreased": sorted(dec_panel), "increased": sorted(inc_panel)},
        "stages": stages_out,
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(payload, f, separators=(",", ":"))
    print(f"wrote {OUT}  ({n_total} embryos, {len(METHODS)} methods)")


if __name__ == "__main__":
    main()
