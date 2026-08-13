#!/usr/bin/env python3
"""Build data/stages.json.gz — "Across the Stages".

A port of figures 8.3, 8.6, 9.1 and 9.2 from HighResSlideshowExports, which is the specification.

THE QUESTION. Take one gene and ask how unevenly it sits between the two halves of an embryo, at
each of the three stages in turn. Does polarization appear, persist, or fall away?

THE TWO HALVES ARE NOT THE SAME OBJECT AT EVERY STAGE, and cannot be:
  · ZYGOTE — the best meridional plane, the pipeline's own `best_planes.pVol` pick from the
    18-plane fan about the polar-body axis. Cytoplasm-only on both the counts and the volumes.
  · 2-CELL — the two blastomeres. Plane-independent; nuclei and the polar body carry their own
    segment labels and never enter either half.

THREE RULES THAT DECIDE WHAT A NUMBER MEANS HERE

  1. A GENE MISSING FROM AN EMBRYO'S RECORDS MEANS TWO DIFFERENT THINGS. Not on that probeset (no
     information, must not enter the average) or on it and never detected (a real zero, which
     must). Panel membership is read off the probeset's gene set, never off the records.

  2. THE ZERO RULE. An embryo with no transcripts of a gene has no defined density ratio, so it
     contributes |log2 ratio| = 0 — "no measurable asymmetry" — rather than being dropped. A gene
     never detected at a stage therefore lands on exactly 1.0 rather than vanishing.

  3. A COUNT FLOOR OF 20. Below ~20 transcripts in an embryo the split is counting noise: 11 vs 8
     is not a measurement of asymmetry. Anything under the floor is treated exactly like a gene
     that was never detected, rather than contributing a large ratio built on single digits.

⚠️ THE BULK CORRECTION HERE IS THE RATIO OF TOTALS, NOT THE MEDIAN OF RATIOS. That differs from
figures 4.19/4.21, and it is deliberate: this is what f8_3.py does, and the reference is followed
per figure rather than harmonised by assumption. The two answer slightly different questions and
the site now contains both, each where its own figure puts it.

WHAT IS RECOMPUTED AND WHAT IS IMPORTED. Figure 8.3 is recomputed here in full from data/zygote/
and data/segments/, and validated against the reference's own shipped table. Figures 8.6/9.1/9.2
rest on two things the reference states cannot be re-derived — the curated retained/lost calls
(which survive both count adjustments, and from which Zbed3 was deliberately dropped) and the
k-means trajectory cluster (fitted on stage-to-stage CHANGES, so it cannot be read off the fold
values). Those columns are imported from the reference tables and labelled as imported. Everything
derivable beside them, including the count-matched percentile, is recomputed.

Output: data/stages.json.gz
"""
import collections
import csv
import glob
import gzip
import json
import math
import os
import sys

import numpy as np

import embryo_stats as ES

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
ZY = os.path.join(DATA, "zygote")
PROBESETS = os.path.join(DATA, "probesets.json")
OUT = os.path.join(DATA, "stages.json.gz")

REF = "/Users/rishib/Desktop/EmbyroPlayground/HighResSlideshowExports"
REF_83 = os.path.join(REF, "Index8", "8.3_fold_change_across_stages", "data.csv")
REF_86 = os.path.join(REF, "Index8", "8.6_retained_lost", "data.csv")
REF_91 = os.path.join(REF, "Index9", "9.1_fold_change_heatmap", "data.csv")
REF_92 = os.path.join(REF, "Index9", "9.2_count_matched_percentile_heatmap", "data.csv")

VERSION = "stages-1.0.0"
EPS = 0.5                 # pseudocount
COUNT_FLOOR = 20          # below this in an embryo, the split is counting noise
NULL_REPS = 40            # count-matched null draws per gene
NULL_SEED = 11
ABUNDANCE_WINDOW = 3.0    # "similar abundance" for the count-matched percentile: within 3x
STAGES = ["zygote", "early2cell", "late2cell"]
STAGE_LABEL = {"zygote": "Zygote", "early2cell": "Early 2-cell", "late2cell": "Late 2-cell"}


def lfc(a, b, vA, vB):
    return math.log2(((a + EPS) / vA) / ((b + EPS) / vB))


def panel_sets(probeset):
    """Which genes each probeset COULD have measured, per source. Membership decides whether an
    absent gene is a real zero or simply not on the panel."""
    zyg, seg = collections.defaultdict(set), collections.defaultdict(set)
    for fp in sorted(glob.glob(os.path.join(ZY, "*.json.gz"))):
        d = json.load(gzip.open(fp, "rt"))
        p = probeset.get(d["id"])
        if p:
            zyg[p] |= {g["gene"] for g in d["analysis"]["genes"]}
    for fp in sorted(glob.glob(os.path.join(ES.SEG, "*.json.gz"))):
        d = json.load(gzip.open(fp, "rt"))
        p = probeset.get(d["id"])
        if p:
            seg[p] |= set(d["transcripts"])
    return zyg, seg


def zygote_rows(probeset, panel):
    """Zygote halves from the best meridional plane. The stored a/b are cytoplasm-only and volA/volB
    are the matching cytoplasm half-volumes, so density is cytoplasmic top and bottom."""
    rows = []
    for fp in sorted(glob.glob(os.path.join(ZY, "*.json.gz"))):
        d = json.load(gzip.open(fp, "rt"))
        eid, an = d["id"], d["analysis"]
        bi = an["best_planes"]["pVol"]
        pl = an["planes"][bi]
        vA, vB = pl["volA"], pl["volB"]
        recs = {g["gene"]: g for g in an["genes"]}
        A = sum(g["planes"][bi]["a"] for g in an["genes"])
        B = sum(g["planes"][bi]["b"] for g in an["genes"])
        if A == 0 or B == 0:
            continue
        bulk = math.log2((A / vA) / (B / vB))
        for gn in sorted(panel.get(probeset.get(eid), set()) | set(recs)):
            g = recs.get(gn)
            if g is None:
                rows.append({"embryo": eid, "gene": gn, "total": 0, "vA": vA, "vB": vB,
                             "lfc": 0.0, "zero": 1})
            else:
                P = g["planes"][bi]
                rows.append({"embryo": eid, "gene": gn, "total": g["total"], "vA": vA, "vB": vB,
                             "lfc": lfc(P["a"], P["b"], vA, vB) - bulk, "zero": 0})
    return rows


def twocell_rows(stage, probeset, panel):
    """Blastomere halves — plane-independent. The two bodies are found BY VOLUME."""
    rows, skipped = [], []
    for fp in sorted(glob.glob(os.path.join(ES.SEG, stage + "__*.json.gz"))):
        d = json.load(gzip.open(fp, "rt"))
        eid = d["id"]
        bodies = ES.classify_body(d)
        if len(bodies) != 2:
            skipped.append({"id": eid, "reason": f"{len(bodies)} body segment(s), not 2"})
            continue
        vol = ES.seg_volumes(d)
        a_lbl, b_lbl = bodies
        vA, vB = vol[a_lbl], vol[b_lbl]
        if min(vA, vB) / max(vA, vB) < 0.5:
            skipped.append({"id": eid, "reason": "a 2:1 volume split is a segmentation failure"})
            continue
        cnt = {}
        for g, T in d["transcripts"].items():
            s = np.asarray(T["s"], int).astype(str)
            cnt[g] = (int((s == a_lbl).sum()), int((s == b_lbl).sum()))
        A = sum(v[0] for v in cnt.values())
        B = sum(v[1] for v in cnt.values())
        if A == 0 or B == 0:
            skipped.append({"id": eid, "reason": "a blastomere holds no transcripts"})
            continue
        bulk = math.log2((A / vA) / (B / vB))
        for gn in sorted(panel.get(probeset.get(eid), set()) | set(cnt)):
            a, b = cnt.get(gn, (0, 0))
            if a + b == 0:
                rows.append({"embryo": eid, "gene": gn, "total": 0, "vA": vA, "vB": vB,
                             "lfc": 0.0, "zero": 1})
            else:
                rows.append({"embryo": eid, "gene": gn, "total": a + b, "vA": vA, "vB": vB,
                             "lfc": lfc(a, b, vA, vB) - bulk, "zero": 0})
    return rows, skipped


def per_gene(rows):
    """fold = 2 ** mean|log2 ratio|, with the floor and the zero rule applied."""
    for r in rows:
        if r["total"] < COUNT_FLOOR:
            r["lfc"], r["zero"] = 0.0, 1
    by = collections.defaultdict(list)
    for r in rows:
        by[r["gene"]].append(r)
    out = {}
    for g, sub in by.items():
        out[g] = {
            "n": len(sub),
            "nMeas": sum(1 for r in sub if not r["zero"]),
            "medTotal": float(np.median([r["total"] for r in sub])),
            "fold": float(2 ** np.mean([abs(r["lfc"]) for r in sub])),
        }
    return out


def null_fold(rows, reps=NULL_REPS, seed=NULL_SEED):
    """What the same arithmetic returns on pure counting noise: keep every embryo's total and
    volume ratio, redraw which half each transcript fell in. The null proportion is the VOLUME
    split, not a fair coin — half the cell is not half the volume."""
    rng = np.random.default_rng(seed)
    by = collections.defaultdict(list)
    for r in rows:
        by[r["gene"]].append(r)
    out = {}
    for g, sub in by.items():
        tot = np.array([r["total"] for r in sub], int)
        vA = np.array([r["vA"] for r in sub], float)
        vB = np.array([r["vB"] for r in sub], float)
        ok = tot > 0
        p = np.where(ok, vA / np.where(ok, vA + vB, 1.0), 0.5)
        vals = []
        for _ in range(reps):
            a = rng.binomial(np.maximum(tot, 0), p)
            b = tot - a
            L = np.where(ok, np.log2(((a + EPS) / np.where(ok, vA, 1.0)) /
                                     ((b + EPS) / np.where(ok, vB, 1.0))), 0.0)
            vals.append(float(np.abs(np.where(ok, L, 0.0)).mean()))
        out[g] = float(2 ** np.mean(vals))
    return out


def count_matched_percentile(stats_):
    """Each gene ranked against genes of SIMILAR ABUNDANCE at the same stage.

    The asymmetry-to-noise ratio rises with transcript count, so an unranked axis would mostly
    plot abundance. "Similar" is within 3x of the gene's own median total."""
    genes = [g for g, s in stats_.items() if s["nMeas"] > 0 and s["medTotal"] > 0]
    ratio = {g: stats_[g]["fold"] / max(stats_[g]["nullFold"], 1e-9) for g in genes}
    out = {}
    for g in genes:
        n = stats_[g]["medTotal"]
        peers = [h for h in genes
                 if n / ABUNDANCE_WINDOW <= stats_[h]["medTotal"] <= n * ABUNDANCE_WINDOW]
        if len(peers) < 3:
            peers = genes
        below = sum(1 for h in peers if ratio[h] < ratio[g])
        out[g] = {"pct": below / max(len(peers) - 1, 1), "ratio": ratio[g], "nPeers": len(peers)}
    return out


def read_ref(path):
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        return list(csv.DictReader(fh))


def main():
    probeset = json.load(open(PROBESETS))
    zyg_panel, seg_panel = panel_sets(probeset)

    rows = {"zygote": zygote_rows(probeset, zyg_panel)}
    skipped = []
    for st, key in (("Early2Cell", "early2cell"), ("Late2Cell", "late2cell")):
        rows[key], sk = twocell_rows(st, probeset, seg_panel)
        skipped += [dict(s, stage=key) for s in sk]

    stats_ = {}
    for st in STAGES:
        s = per_gene(rows[st])
        nf = null_fold(rows[st])
        for g in s:
            s[g]["nullFold"] = nf.get(g, 1.0)
            s[g]["excess"] = s[g]["fold"] - s[g]["nullFold"]
        pct = count_matched_percentile(s)
        for g in s:
            s[g].update(pct.get(g, {"pct": None, "ratio": None, "nPeers": 0}))
        stats_[st] = s
        print(f"  {STAGE_LABEL[st]:14s} {len(s)} genes, "
              f"{sum(1 for v in s.values() if v['nMeas'] > 0)} measured somewhere, "
              f"median fold {np.median([v['fold'] for v in s.values()]):.3f}, "
              f"median null {np.median([v['nullFold'] for v in s.values()]):.3f}")

    genes = sorted(set().union(*[set(stats_[st]) for st in STAGES]))

    # ---- imported: what the reference states cannot be re-derived ----
    imported = {"available": {}, "note": ""}
    r91 = read_ref(REF_91)
    r92 = read_ref(REF_92)
    r86 = read_ref(REF_86)
    cluster = {r["gene"]: r["group"] for r in r91} if r91 else {}
    curated = {r["gene"]: r["group"] for r in (r92 or r86 or [])}
    ref_pct = {}
    if r92:
        for r in r92:
            ref_pct[r["gene"]] = {"e": float(r["e_pct"]), "l": float(r["l_pct"]),
                                  "eN": int(float(r["e_n"])), "lN": int(float(r["l_n"])),
                                  "eRatio": float(r["e_ratio"]), "lRatio": float(r["l_ratio"])}
    imported["available"] = {"cluster": bool(cluster), "curated": bool(curated),
                             "refPct": bool(ref_pct)}
    imported["note"] = (
        "The k-means trajectory CLUSTER is fitted on stage-to-stage log2 CHANGES, so it cannot be "
        "read off the three fold values; the retained/lost CALL is curated and survives both count "
        "adjustments (Zbed3 was dropped for failing that check). The reference states both must "
        "travel as columns rather than be re-derived, so they are imported here and marked as "
        "imported. The count-matched percentile beside them IS recomputed.")

    # the derived third fate, by the reference's own symmetric rule
    if curated and ref_pct:
        lost = [g for g, v in curated.items() if v == "lost" and g in ref_pct]
        drop = float(np.median([ref_pct[g]["e"] - ref_pct[g]["l"] for g in lost])) if lost else 0.0
        for g, v in ref_pct.items():
            if curated.get(g, "other") == "other" and (v["l"] - v["e"]) >= drop:
                curated[g] = "gained"
        imported["gained_rule"] = (
            f"DERIVED, not curated: an unclassified gene whose percentile rises by at least as "
            f"much as the median lost gene's falls ({drop:.4f}). Symmetric by construction, so "
            f"the three fates are commensurable.")

    out_genes = []
    for g in genes:
        rec = {"g": g, "cluster": cluster.get(g), "group": curated.get(g, "other")}
        for st in STAGES:
            s = stats_[st].get(g)
            rec[st] = None if s is None else {
                "fold": round(s["fold"], 5), "null": round(s["nullFold"], 5),
                "excess": round(s["excess"], 5), "n": s["n"], "nMeas": s["nMeas"],
                "medTotal": s["medTotal"],
                "pct": None if s["pct"] is None else round(s["pct"], 5),
                "ratio": None if s["ratio"] is None else round(s["ratio"], 5),
            }
        rec["refPct"] = ref_pct.get(g)
        out_genes.append(rec)

    # ---- validation against the reference's own 8.3 table ----
    val = {"available": False}
    r83 = read_ref(REF_83)
    if r83:
        ref = {r["gene"]: r for r in r83}
        d = {st: [] for st in STAGES}
        for rec in out_genes:
            r = ref.get(rec["g"])
            if not r:
                continue
            for st in STAGES:
                key = {"zygote": "fold_zygote", "early2cell": "fold_early2cell",
                       "late2cell": "fold_late2cell"}[st]
                if rec[st] and r.get(key):
                    d[st].append(abs(rec[st]["fold"] - float(r[key])))
        val = {"available": True, "n_shared": len(set(ref) & {g["g"] for g in out_genes}),
               "max_abs_diff": {st: (round(max(d[st]), 6) if d[st] else None) for st in STAGES},
               "median_abs_diff": {st: (round(float(np.median(d[st])), 6) if d[st] else None)
                                   for st in STAGES}}
        print(f"\n  vs the reference's 8.3 table ({val['n_shared']} shared genes):")
        for st in STAGES:
            print(f"    {STAGE_LABEL[st]:14s} median |Δfold| {val['median_abs_diff'][st]}, "
                  f"max {val['max_abs_diff'][st]}")

    doc = {
        "meta": {
            "version": VERSION,
            "method": "figures 8.3 / 8.6 / 9.1 / 9.2 — the specification",
            "stages": STAGES, "stageLabel": STAGE_LABEL,
            "params": {"EPS": EPS, "COUNT_FLOOR": COUNT_FLOOR, "NULL_REPS": NULL_REPS,
                       "ABUNDANCE_WINDOW": ABUNDANCE_WINDOW},
            "halves": {"zygote": "the best meridional plane (best_planes.pVol of the 18-plane fan "
                                 "about the polar-body axis); cytoplasm-only",
                       "twocell": "the two blastomeres, found by volume; nuclei and the polar "
                                  "body carry their own labels and enter neither half"},
            "bulk": "the RATIO OF TOTALS, per figure 8.3 — not the median of ratios that "
                    "figures 4.19/4.21 use. Followed per figure rather than harmonised.",
            "zeroRule": "an embryo with no transcripts of a gene contributes |log2 ratio| = 0, "
                        "rather than being dropped",
            "floorRule": f"under {COUNT_FLOOR} transcripts in an embryo, the split is counting "
                         f"noise and is treated exactly like a non-detection",
            "panelRule": "a gene absent from an embryo's records is a real zero only if it is on "
                         "that probeset; otherwise it carries no information and is not counted",
            "n_genes": len(out_genes), "skipped": skipped,
            "imported": imported, "validation": val,
        },
        "genes": out_genes,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e6:.2f} MB)")
    print(f"  {len(out_genes)} genes; {sum(1 for g in out_genes if g['cluster'])} carry an "
          f"imported trajectory cluster, {sum(1 for g in out_genes if g['group'] != 'other')} a "
          f"retained/lost/gained call")


if __name__ == "__main__":
    main()
