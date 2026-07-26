#!/usr/bin/env python3
"""
Tests for the Stage Expression Explorer aggregate (data/stage_expr.json.gz).

Designed to FAIL on the ways this could silently drift from the source data:
  * counts not reproducing segments_genes.json.gz (the number the plots show);
  * samples landing in the wrong stage, or Oocyte leaking into the 3 compared stages;
  * a sample missing its 3-D silhouette (an empty grid card);
  * the Welch p-values or per-stage means not matching a direct computation;
  * absolute paths leaking into the committed artifact.

Run: python3 scripts/test_stage_expr.py
"""
from __future__ import annotations

import gzip
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AGG = os.path.join(HERE, "data", "stage_expr.json.gz")
SRC = os.path.join(HERE, "data", "segments_genes.json.gz")

FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def mean(v):
    return sum(v) / len(v) if v else None


def variance(v):
    if len(v) < 2:
        return 0.0
    m = mean(v)
    return sum((x - m) ** 2 for x in v) / (len(v) - 1)


def welch_p(a, b):
    """Two-sided Welch p via scipy if present, else the t-survival used on the page."""
    try:
        from scipy import stats
        return float(stats.ttest_ind(a, b, equal_var=False).pvalue)
    except Exception:
        sa, sb = variance(a) / len(a), variance(b) / len(b)
        t = (mean(a) - mean(b)) / math.sqrt(sa + sb)
        df = (sa + sb) ** 2 / (sa ** 2 / (len(a) - 1) + sb ** 2 / (len(b) - 1))
        # regularized incomplete beta (same as the page)
        import mpmath  # noqa
        return float(mpmath.betainc(df / 2, 0.5, 0, df / (df + t * t), regularized=True))


def main():
    print("stage expression explorer — tests\n")
    if not os.path.isfile(AGG):
        print("  FAIL  aggregate missing — run python3 build_stage_expr.py")
        return 1
    d = json.load(gzip.open(AGG, "rt"))
    S = d["samples"]

    # ── stages ──
    print("[stages & samples]")
    check("three compared stages, in order", d["stages"] == ["Zygote", "Early2Cell", "Late2Cell"])
    stages = {s["stage"] for s in S}
    check("no Oocyte leaked into the compared stages", "Oocyte" not in stages, str(stages))
    from collections import Counter
    sc = Counter(s["stage"] for s in S)
    check("expected per-stage sample counts (60/55/38)",
          (sc["Zygote"], sc["Early2Cell"], sc["Late2Cell"]) == (60, 55, 38), str(dict(sc)))
    check("every sample has a 3-D silhouette", all(s["hull"] and len(s["hull"]) >= 3 for s in S),
          f"{sum(1 for s in S if not s['hull'])} without a hull")
    check("every sample has a positive volume and total_tx",
          all(s["vol"] > 0 and s["total_tx"] > 0 for s in S))

    # ── counts reproduce the source ──
    print("\n[counts reproduce segments_genes.json.gz]")
    src = json.load(gzip.open(SRC, "rt"))
    src_emb = src["embInfo"]
    label2stage = {"Zygote": "Zygote", "Early 2-cell": "Early2Cell", "Late 2-cell": "Late2Cell"}
    # rebuild the expected per-sample count for a few genes straight from source
    mism = 0
    for g in ("Pard3", "Nlrp5", "MuERV-L", "Yap1"):
        if g not in d["gene_counts"] or g not in src["genes"]:
            continue
        exp = {}
        for embIdx, segs in src["genes"][g]:
            e = src_emb[embIdx]
            if e["stage"] not in label2stage:
                continue
            exp[e["id"]] = int(segs[0][3])
        got = {S[int(si)]["id"]: c for si, c in d["gene_counts"][g].items()}
        if exp != got:
            mism += 1
    check("per-sample counts match the source for Pard3/Nlrp5/MuERV-L/Yap1", mism == 0,
          f"{mism} genes differ")

    # ── Pard3 reproduces Harry's published numbers ──
    print("\n[Pard3 reproduces the source explorer]")
    by = {"Zygote": [], "Early2Cell": [], "Late2Cell": []}
    for si, c in d["gene_counts"]["Pard3"].items():
        by[S[int(si)]["stage"]].append(c)
    zn, en, ln = len(by["Zygote"]), len(by["Early2Cell"]), len(by["Late2Cell"])
    check("Pard3 sample counts 15 / 15 / 9", (zn, en, ln) == (15, 15, 9), f"{zn}/{en}/{ln}")
    zm, em = mean(by["Zygote"]), mean(by["Early2Cell"])
    check("Pard3 zygote mean ≈ 2291", abs(zm - 2291) < 5, f"{zm:.0f}")
    check("Pard3 early-2c mean ≈ 4162", abs(em - 4162) < 5, f"{em:.0f}")
    try:
        p_ze = welch_p(by["Zygote"], by["Early2Cell"])
        p_el = welch_p(by["Early2Cell"], by["Late2Cell"])
        check("Welch Zyg vs Early ≈ 0.015", abs(p_ze - 0.015) < 0.004, f"{p_ze:.4f}")
        check("Welch Early vs Late highly significant", p_el < 0.001, f"{p_el:.5f}")
    except Exception as e:
        check("Welch p-values computable", False, str(e))

    # ── transcript dots present and in-frame ──
    print("\n[grid transcript dots]")
    tx = d["gene_tx"].get("Pard3", {})
    check("Pard3 has projected dots for detected samples", len(tx) > 20, str(len(tx)))
    allpts = [p for arr in tx.values() for p in arr]
    check("dots are inside the 0..1000 frame",
          all(0 <= x <= 1000 and 0 <= y <= 1000 for x, y in allpts))
    check("dots downsampled to the cap", all(len(arr) <= d["dot_cap"] for arr in tx.values()))

    # ── provenance ──
    print("\n[provenance]")
    txt = open(AGG, "rb").read()
    leaks = [s for s in (b"/Users/", b"/Volumes/", b"C:\\") if s in txt]
    check("no absolute paths in the artifact", not leaks, str(leaks))

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"  FAILED: {f}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
