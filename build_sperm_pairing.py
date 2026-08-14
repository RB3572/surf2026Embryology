#!/usr/bin/env python3
"""Build data/sperm_pairing.json — maternal vs paternal, paired within each zygote (figure 4.16).

THE QUESTION. The two pronuclei of a zygote are not interchangeable: one came from the oocyte and
one from the sperm. Once each is labelled, three geometric quantities can be compared BETWEEN
PARENTS INSIDE THE SAME CELL:

  · com_sperm  — how far each pronucleus's centre of mass sits from the sperm entry point
  · volume     — each pronucleus's enclosed volume
  · com_polar  — how far each sits from the polar body's centre of mass

Pairing within the embryo is what makes this worth doing: embryo size, orientation, segmentation
quality and detection efficiency all divide out, so a Wilcoxon signed-rank test on the within-cell
differences is the whole analysis. n is small; the test is exact.

⚠️ EVERY ONE OF THESE COMPARISONS IS CIRCULAR AGAINST THE CONSENSUS THAT LABELLED THE PRONUCLEI.
The site's maternal/paternal call is a majority vote of four tests — pb_com, pb_shell, sperm and
volume — and each comparison above is exactly one of those tests, re-asked as a measurement. Run
naively, "the sperm is closer to the paternal pronucleus" would be partly a restatement of the
vote that decided which one is paternal.

So each comparison is ALSO run on a LEAVE-ONE-OUT consensus that drops the test(s) it would
restate, and both numbers are shipped side by side:

  com_sperm  drops  sperm
  volume     drops  volume
  com_polar  drops  pb_com AND pb_shell   (both read the same landmark)

A zygote whose remaining tests tie is dropped from the leave-one-out variant rather than guessed
at. Hand-made calls are not votes at all — they come from a person looking at the images — so they
are carried into both variants unchanged and counted separately.

  THE LEAVE-ONE-OUT NUMBER IS THE ONE TO QUOTE. The naive one is shipped so the size of the
  circularity is visible rather than argued about.

Geometry comes from data/pronuclei_assignments.json, in the viewer's isotropic plot space;
µm = plot distance × 0.15. Volumes are the segmentation's own voxel volumes, already µm³.

Output: data/sperm_pairing.json
"""
import json
import math
import os

from scipy import stats

# The embryo label is LOOKED UP (data/embryo_ids.json via embryo_naming), never derived and
# never read off a manifest — rebuilding an artifact must not quietly reintroduce a legacy
# name. embryo_label() falls back conspicuously when an embryo is missing from the lookup.
from embryo_naming import embryo_label

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
ASSIGN = os.path.join(DATA, "pronuclei_assignments.json")
OUT = os.path.join(DATA, "sperm_pairing.json")

VERSION = "sperm-pairing-1.0.0"
UM = 0.15                 # plot-space unit → µm (isotropic), as in build_sperm_pseudotime.py
ALL_TESTS = ("pb_com", "pb_shell", "sperm", "volume")

# key -> (label, which tests it would restate, the units, what "maternal larger" means)
COMPARISONS = {
    "com_sperm": {
        "label": "Centre of mass → sperm entry point",
        "drops": ("sperm",),
        "unit": "µm",
        "needs": "sperm",
        "expect": "The sperm-derived pronucleus should be the nearer one, so maternal − paternal "
                  "should be POSITIVE.",
    },
    "volume": {
        "label": "Enclosed volume of the pronucleus",
        "drops": ("volume",),
        "unit": "µm³",
        "needs": None,
        "expect": "This one has no agreed direction. The pipeline's volume heuristic ASSUMES "
                  "the smaller pronucleus is maternal, which is the weakest of the four tests "
                  "(it disagrees with its own consensus on one zygote in five), and the hand "
                  "calls on this cohort went the other way. The leave-one-out column is the only "
                  "informative one here.",
    },
    "com_polar": {
        "label": "Centre of mass → polar body",
        "drops": ("pb_com", "pb_shell"),
        "unit": "µm",
        "needs": "polar",
        "expect": "The maternal pronucleus should sit nearer the polar body — both came from the "
                  "same meiotic division — so maternal − paternal should be NEGATIVE.",
    },
}


def d3(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def loo_female(rec, drop):
    """The consensus with `drop` removed. None when the rest tie or nothing is left.

    A hand call is not a vote — it stands whatever the tests say, so it survives every
    leave-one-out."""
    if (rec.get("consensus") or {}).get("manual"):
        return rec["consensus"].get("female"), "manual"
    votes = [t["female"] for k, t in (rec.get("tests") or {}).items()
             if k not in drop and isinstance(t, dict) and t.get("female") is not None]
    if not votes:
        return None, "no tests left"
    n0 = sum(1 for v in votes if v == 0)
    n1 = len(votes) - n0
    if n0 == n1:
        return None, "the remaining tests tie"
    return (0 if n0 > n1 else 1), f"{max(n0, n1)}/{len(votes)}"


def measure(rec, key, female):
    """(maternal, paternal) for one comparison, or None if the geometry is not there."""
    pron = rec["pron"]
    mat, pat = pron[female], pron[1 - female]
    if key == "volume":
        if mat.get("volume") is None or pat.get("volume") is None:
            return None
        return float(mat["volume"]), float(pat["volume"])
    if key == "com_sperm":
        sp = rec.get("sperm_plot")
        if not sp:
            return None
        return d3(mat["com_plot"], sp) * UM, d3(pat["com_plot"], sp) * UM
    if key == "com_polar":
        pb = (rec.get("polar") or {}).get("com_plot")
        if not pb:
            return None
        return d3(mat["com_plot"], pb) * UM, d3(pat["com_plot"], pb) * UM
    raise KeyError(key)


def wilcoxon(pairs):
    """Paired signed-rank on maternal − paternal. Exact: n is small."""
    d = [m - p for m, p in pairs]
    nz = [x for x in d if x != 0]
    if len(nz) < 5:
        return {"p": None, "stat": None, "note": f"only {len(nz)} non-zero differences — no test"}
    try:
        st = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox",
                            mode="exact" if len(nz) <= 25 else "approx")
    except TypeError:                                   # newer scipy renamed the argument
        st = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox",
                            method="exact" if len(nz) <= 25 else "approx")
    return {"p": float(st.pvalue), "stat": float(st.statistic),
            "exact": bool(len(nz) <= 25), "n_nonzero": len(nz)}


def summarise(rows):
    d = [r["m"] - r["p"] for r in rows]
    d_sorted = sorted(d)
    n = len(d)
    med = (d_sorted[n // 2] if n % 2 else (d_sorted[n // 2 - 1] + d_sorted[n // 2]) / 2) if n else None
    out = {"n": n, "rows": rows,
           "median_diff": round(med, 4) if med is not None else None,
           "mean_diff": round(sum(d) / n, 4) if n else None,
           "n_maternal_larger": sum(1 for x in d if x > 0),
           "n_paternal_larger": sum(1 for x in d if x < 0)}
    out.update(wilcoxon([(r["m"], r["p"]) for r in rows]) if n else {"p": None})
    return out


def main():
    recs = json.load(open(ASSIGN))["embryos"]
    manual = {r["id"] for r in recs if (r.get("consensus") or {}).get("manual")}

    comparisons = {}
    for key, spec in COMPARISONS.items():
        variants = {}
        for vname, drop in (("all", ()), ("loo", spec["drops"])):
            rows, dropped = [], []
            for rec in recs:
                c = rec.get("consensus") or {}
                if vname == "all":
                    female = None if c.get("split") else c.get("female")
                    why = "manual" if rec["id"] in manual else f"{max(c.get('n0', 0), c.get('n1', 0))}/4"
                else:
                    female, why = loo_female(rec, drop)
                if female is None:
                    dropped.append({"id": rec["id"],
                                    "reason": why if vname == "loo" else "the four tests split"})
                    continue
                mp = measure(rec, key, female)
                if mp is None:
                    dropped.append({"id": rec["id"], "reason": f"no {spec['needs'] or key} geometry"})
                    continue
                rows.append({"id": rec["id"], "label": embryo_label(rec["id"]),
                             "m": round(mp[0], 4), "p": round(mp[1], 4),
                             "call": why, "manual": rec["id"] in manual})
            rows.sort(key=lambda r: r["m"] - r["p"])
            v = summarise(rows)
            v["dropped"] = dropped
            v["n_manual"] = sum(1 for r in rows if r["manual"])
            variants[vname] = v
        comparisons[key] = {
            "label": spec["label"], "unit": spec["unit"], "expect": spec["expect"],
            "drops": list(spec["drops"]), "variants": variants,
            "circular": f"This comparison IS the {' and '.join(spec['drops'])} test"
                        f"{'s' if len(spec['drops']) > 1 else ''}, re-asked as a measurement. "
                        f"Quote the leave-one-out column.",
        }

    doc = {
        "meta": {
            "version": VERSION,
            "method": "figure 4.16 — maternal vs paternal, paired within each zygote",
            "unit_um_per_plot": UM,
            "pairing": "every comparison is WITHIN one cell, so embryo size, orientation and "
                       "detection efficiency divide out; the test is a paired Wilcoxon "
                       "signed-rank, exact because n is small",
            "circularity": "the maternal/paternal call is a majority of four tests — pb_com, "
                           "pb_shell, sperm, volume — and each comparison here is one of them "
                           "re-asked as a measurement. Every comparison is therefore ALSO run on "
                           "a consensus with that test (or, for the polar-body distance, both "
                           "polar-body tests) removed. THE LEAVE-ONE-OUT NUMBER IS THE ONE TO "
                           "QUOTE; the naive one is shipped so the size of the circularity is "
                           "visible rather than argued about.",
            "ties": "a zygote whose remaining tests tie is dropped from the leave-one-out variant "
                    "rather than guessed at",
            "manual": "hand-made calls are not votes — they come from a person looking at the "
                      "images — so they survive every leave-one-out unchanged and are counted "
                      "separately",
            "n_embryos": len(recs), "n_manual": len(manual),
            "tests": list(ALL_TESTS),
        },
        "comparisons": comparisons,
    }
    with open(OUT, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(f"  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1024:.0f} KB)")
    for k, c in comparisons.items():
        for vn in ("all", "loo"):
            v = c["variants"][vn]
            p = "–" if v.get("p") is None else f"{v['p']:.4f}"
            print(f"  {k:10s} {vn:4s}  n={v['n']:2d}  median Δ {v['median_diff']:+9.3f} "
                  f"{c['unit']:4s}  {v['n_maternal_larger']}/{v['n']} maternal larger  P={p}"
                  + (f"  ({v['n_manual']} manual)" if v["n_manual"] else ""))


if __name__ == "__main__":
    main()
