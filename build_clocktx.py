#!/usr/bin/env python3
"""Build data/clocktx.json.gz — "Transcriptome vs the Clock".

A port of figures 4.8, 4.11 and 5.4 from HighResSlideshowExports, which is the specification.

THE QUESTION. Every fixed zygote carries an inferred pseudotime τ. Does any gene's presence in the
transcriptome move with it?

WHAT IS CORRELATED, AND WHY IT IS A SHARE. A gene's SHARE is its transcript count divided by that
embryo's total — composition, not abundance. Raw counts differ by an order of magnitude between
probesets and between embryos for reasons that have nothing to do with time, and a share removes
that. Concentration (count ÷ cytoplasm volume) ships alongside because it asks a genuinely
different question: a cell that shrinks raises concentration at constant share.

⚠️ THE REGION SET IS NOT CYTOPLASM-ONLY, AND THAT IS DELIBERATE. Unlike every plane analysis on
this site, a share is about the composition of the embryo's transcriptome, so it counts the
CYTOPLASM AND THE PRONUCLEI. The polar body is excluded — it is a discarded cell, not part of the
zygote's transcriptome — but a variant including it ships too, so the choice is inspectable rather
than buried. The denominator is always summed over exactly the same region set as the numerator: a
share whose numerator and denominator have different physical scope is not a share.

CENTRED WITHIN PROBESET. For each gene, its share is demeaned separately within each probeset
group before correlating. The four panels differ by more than an order of magnitude in total
counts, and without this a gene's apparent trend would just track which panel it was measured on.
For a gene seen in one probeset this is a constant shift and changes neither ρ nor P.

SMALL n IS HANDLED EXACTLY. scipy's Spearman P is an asymptotic approximation that returns a
literal 0.0 for a perfect rank match — impossible at n = 3, where only 3! = 6 orderings exist and
the smallest two-sided P is 1/3. Below n = 6, and whenever scipy returns a non-positive P at
n ≤ 10, the P is computed by exhaustive permutation instead.

THE ENTRY FLOOR is ≥ 7 counts in a zygote for that zygote to qualify, and ≥ 2 qualifying zygotes
for the gene to be testable. 7 is the reference's swept value, not a round number: it is where
MuERV-L first reaches a non-trivial sample without its trend being diluted by noisier low-count
zygotes.

Output: data/clocktx.json.gz
"""
import gzip
import json
import math
import os
import sys
from itertools import permutations

import numpy as np
from scipy import stats

import embryo_stats as ES
# The embryo label is LOOKED UP (data/embryo_ids.json via embryo_naming), never derived and
# never read off a manifest — rebuilding an artifact must not quietly reintroduce a legacy
# name. embryo_label() falls back conspicuously when an embryo is missing from the lookup.
from embryo_naming import embryo_label

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
PT = os.path.join(HERE, "..", "embryo_viewer", "public", "data", "pronuclei_pseudotime.json")
PROBESETS = os.path.join(DATA, "probesets.json")
ZY_MAN = os.path.join(DATA, "zygote_manifest.json")
OUT = os.path.join(DATA, "clocktx.json.gz")

VERSION = "clocktx-1.0.0"
MIN_COUNT = 7             # a zygote qualifies for a gene at this many transcripts
MIN_ZYGOTES = 2           # ...and a gene needs this many qualifying zygotes to be testable
EXACT_BELOW_N = 6         # below this n, the P is exact by permutation
EXACT_MAX_N = 10          # ...and up to this n if scipy returns a non-positive P
EXEMPLARS = ["Zp1", "Zp2", "Pin1", "Cdc42"]     # the reference's four (figure 4.11)


def exact_spearman(x, y):
    """Exhaustive permutation P for Spearman ρ. Exact, and at this n, cheap."""
    from scipy.stats import rankdata
    rx, ry = rankdata(x), rankdata(y)
    rho_obs = float(np.corrcoef(rx, ry)[0, 1])
    tot = hit = 0
    for perm in permutations(ry):
        tot += 1
        if abs(float(np.corrcoef(rx, np.asarray(perm))[0, 1])) >= abs(rho_obs) - 1e-12:
            hit += 1
    return rho_obs, hit / tot


def correlate(rows, tau, probeset, ycol):
    """ρ and P per gene, on the probeset-centred value. `rows` is a list of dicts."""
    by = {}
    for r in rows:
        if r["total"] >= MIN_COUNT:
            by.setdefault(r["gene"], []).append(r)
    out, n_single = [], 0
    for g, sub in by.items():
        n = len(sub)
        rec = {"g": g, "n": n,
               "medCount": float(np.median([r["total"] for r in sub])),
               "nProbesets": len(({probeset.get(r["id"], "?") for r in sub}))}
        if n < MIN_ZYGOTES:
            n_single += 1
            rec.update(rho=None, p=None)
            out.append(rec)
            continue
        y = np.array([r[ycol] for r in sub], float)
        t = np.array([tau[r["id"]] for r in sub], float)
        ps = [probeset.get(r["id"], "?") for r in sub]
        # centre within probeset — the four panels differ by >10x in total counts
        c = y.copy()
        for p in set(ps):
            m = np.array([q == p for q in ps])
            c[m] -= c[m].mean()
        if c.std() == 0 or t.std() == 0:
            rec.update(rho=None, p=None)
            out.append(rec)
            continue
        if n < EXACT_BELOW_N:
            rho, p = exact_spearman(c, t)
        else:
            rho, p = stats.spearmanr(c, t)
            if not (p > 0) and n <= EXACT_MAX_N:      # scipy reports 0 for a perfect rank match
                rho, p = exact_spearman(c, t)
        rec.update(rho=float(rho), p=float(p))
        out.append(rec)
    testable = [r for r in out if r["p"] is not None]
    if testable:
        for r, q in zip(testable, ES.bh(np.array([r["p"] for r in testable]))):
            r["fdr"] = float(q)
    for r in out:
        r.setdefault("fdr", None)
    out.sort(key=lambda r: (r["p"] is None, r["p"] if r["p"] is not None else 1.0))
    return out, n_single


def main():
    if not os.path.isfile(PT):
        sys.exit(f"missing pseudotime: {PT}")
    pt = {e["id"]: e for e in json.load(open(PT))["embryos"] if e.get("tau") is not None}
    probeset = json.load(open(PROBESETS)) if os.path.isfile(PROBESETS) else {}
    man = {m["id"]: m for m in json.load(open(ZY_MAN))["embryos"]}

    # ---- per (embryo, gene) counts over the two region sets ----
    per, per_polar, emb_meta, skipped = [], [], [], []
    for eid in sorted(pt):
        p = ES.scene_path("Zygote", eid)
        if not os.path.isfile(p):
            skipped.append({"id": eid, "reason": "no segments scene"})
            continue
        sc = ES.read_scene(p)
        bodies = ES.classify_body(sc)
        if not bodies:
            skipped.append({"id": eid, "reason": "no body segment"})
            continue
        body = bodies[0]
        polar = ES.polar_label(sc)
        # body + pronuclei; the polar body is a discarded cell, not the zygote's transcriptome
        keep = {str(k) for k in ES.seg_volumes(sc)} - ({str(polar)} if polar else set())
        keep_polar = {str(k) for k in ES.seg_volumes(sc)}
        V, F = ES.mesh_of(sc, body)
        cyto = ES.seg_volumes(sc)[body]
        for g, t in sc["transcripts"].items():
            s = np.asarray(t["s"], int).astype(str)
            per.append({"id": eid, "gene": g, "total": int(np.isin(s, list(keep)).sum())})
            per_polar.append({"id": eid, "gene": g,
                              "total": int(np.isin(s, list(keep_polar)).sum())})
        emb_meta.append({"id": eid, "label": embryo_label(eid),
                         "tau": round(float(pt[eid]["tau"]), 5),
                         "qc": pt[eid].get("qc"),
                         "probeset": probeset.get(eid, "?"),
                         "cyto_vol": round(float(cyto), 1),
                         "polar_excluded": polar is not None})
    if not per:
        sys.exit("no zygote produced counts")

    tau = {e["id"]: e["tau"] for e in emb_meta}
    cyto_of = {e["id"]: e["cyto_vol"] for e in emb_meta}
    keep_ids = set(tau)
    per = [r for r in per if r["id"] in keep_ids]
    per_polar = [r for r in per_polar if r["id"] in keep_ids]

    def assemble(rows):
        tot = {}
        for r in rows:
            tot[r["id"]] = tot.get(r["id"], 0) + r["total"]
        for r in rows:
            r["share"] = r["total"] / max(tot[r["id"]], 1)
            r["conc"] = r["total"] / cyto_of[r["id"]]
        return tot

    tot_main = assemble(per)
    tot_polar = assemble(per_polar)

    variants = {}
    for name, rows in (("main", per), ("withPolar", per_polar)):
        for ycol in ("share", "conc"):
            g, n_single = correlate(rows, tau, probeset, ycol)
            variants[f"{name}.{ycol}"] = {"genes": g, "n_single": n_single}
            sig = [r for r in g if r["p"] is not None and r["p"] < 0.05]
            print(f"  {name}.{ycol}: {len(g)} genes clear the floor, "
                  f"{sum(1 for r in g if r['p'] is not None)} testable, "
                  f"{len(sig)} nominal at P<0.05 "
                  f"({sum(1 for r in sig if r['rho'] < 0)} down, "
                  f"{sum(1 for r in sig if r['rho'] > 0)} up)")

    # per-gene trajectories, for the exemplar panel and for click-through
    traj = {}
    for r in per:
        traj.setdefault(r["gene"], []).append(
            {"id": r["id"], "n": r["total"], "share": round(r["share"], 8),
             "conc": round(r["conc"], 8)})
    for g in traj:
        traj[g].sort(key=lambda r: tau[r["id"]])

    for e in emb_meta:
        e["total_tx"] = int(tot_main[e["id"]])
        e["total_tx_polar"] = int(tot_polar[e["id"]])

    doc = {
        "meta": {
            "version": VERSION,
            "method": "figures 4.8 / 4.11 / 5.4 (_work/cells/f4_8.py) — the specification",
            "n_embryos": len(emb_meta), "n_genes": len(traj),
            "skipped": skipped,
            "params": {"MIN_COUNT": MIN_COUNT, "MIN_ZYGOTES": MIN_ZYGOTES,
                       "EXACT_BELOW_N": EXACT_BELOW_N, "EXACT_MAX_N": EXACT_MAX_N},
            "exemplars": EXEMPLARS,
            "regions": "cytoplasm + pronuclei; the polar body is excluded (a discarded cell, not "
                       "the zygote's transcriptome). The 'with polar body' variant includes it.",
            "denominator": "summed over exactly the same region set as the numerator",
            "centring": "each gene's value is demeaned within probeset before correlating",
            "smallN": "P is exact by exhaustive permutation below n=6, and whenever scipy "
                      "returns a non-positive P at n<=10",
        },
        "embryos": emb_meta,
        "variants": variants,
        "traj": traj,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e6:.2f} MB)")
    print(f"  {len(emb_meta)} zygotes with a tau, {len(traj)} genes")
    top = [r for r in variants["main.share"]["genes"] if r["p"] is not None][:6]
    print("  strongest (share): " + ", ".join(
        f"{r['g']} (rho={r['rho']:+.2f}, P={r['p']:.1e}, n={r['n']})" for r in top))


if __name__ == "__main__":
    main()
