#!/usr/bin/env python3
"""Checks on data/halves.json.gz — the two halves of a zygote, four ways of cutting it.

This artifact has one load-bearing assumption and three traps, and the checks are aimed at them
rather than at arithmetic:

  · THE ORIENTATION. Side F must be the half with MORE cytoplasmic transcripts, in every embryo,
    on every plane. If a single embryo were oriented the other way its genes would enter the
    average with the wrong sign, and nothing downstream would notice.
  · THE ALIGNMENT NULL must be a real null: it has to move when the sides are shuffled, and the
    reported P has to be the tail it claims. A null that always returned the observed count would
    pass every other check on this page.
  · THE RANDOM PLANE must NOT sit at fold 1.0. Counting noise alone produces asymmetry, and if
    the control read 1.0 it would mean the fold is being computed against nothing.
  · THE PAIRING TRAP must be present, not fixed. Under share features the halves are exact
    complements, so the ratio MUST come out above 1 — if it did not, the features would not be
    what the page says they are. The concentration variant is the one allowed to fall below 1.
"""
import gzip
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import embryo_stats as ES                                          # noqa: E402

ART = os.path.join(ROOT, "data", "halves.json.gz")

passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1; print(f"  PASS  {name}")
    else:
        failed += 1; print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))


def section(t):
    print(f"\n[{t}]")


def main():
    if not os.path.isfile(ART):
        sys.exit("halves.json.gz missing — run: python3 build_halves.py")
    d = json.load(gzip.open(ART, "rt"))
    m, P = d["meta"], d["meta"]["params"]

    section("shape")
    check("both volcano planes present", sorted(d["volcano"]) == ["polar18", "sperm"])
    for k, s in d["volcano"].items():
        check(f"{k}: no duplicate genes", len({g['g'] for g in s['genes']}) == len(s["genes"]))
        check(f"{k}: ranks are 1..n in P order",
              [g["rank"] for g in s["genes"]] == list(range(1, len(s["genes"]) + 1))
              and all(s["genes"][i]["p"] <= s["genes"][i + 1]["p"] + 1e-12
                      for i in range(len(s["genes"]) - 1)))
        check(f"{k}: the called count is the real count",
              s["n_called"] == sum(1 for g in s["genes"] if g["p"] < P["CALL_P"]))
        check(f"{k}: one orientation record per embryo",
              len(s["orientation"]) == s["n_embryos"])
    check("every skipped embryo records why", all(x.get("reason") for x in m.get("skipped", [])))
    check("the sperm plane has fewer embryos than the meridional one — not every zygote has one",
          d["volcano"]["sperm"]["n_embryos"] < d["volcano"]["polar18"]["n_embryos"])

    section("side F really is the fuller half")
    for k, s in d["volcano"].items():
        bad = [o["id"] for o in s["orientation"] if o["totF"] < o["totE"]]
        check(f"{k}: every embryo's F half holds at least as many transcripts", not bad,
              str(bad[:4]))
        check(f"{k}: the fraction in the fuller half is at least a half",
              all(o["frac"] >= 0.5 - 1e-9 for o in s["orientation"]))
        # the flip must actually happen sometimes, or the stored side was already the fuller one
        # by luck and the rule is untested
        nf = sum(1 for o in s["orientation"] if o["flipped"])
        check(f"{k}: the orientation rule bites — some embryos had to be flipped", nf > 0, str(nf))
        check(f"{k}: ...but not all of them", nf < len(s["orientation"]), str(nf))
        print(f"        {k}: {nf}/{len(s['orientation'])} flipped")
    check("the orientation rule is stated, with why the sperm cannot do it",
          "more cytoplasmic transcripts" in m["orientation"].lower()
          and "through" in m["orientation"].lower())

    section("the fold is recomputable and it is bulk-centred")
    for k, s in d["volcano"].items():
        worst = 0.0
        for g in s["genes"]:
            worst = max(worst, abs(float(np.mean([x["lfc"] for x in g["per"]])) - g["lfc"]))
        check(f"{k}: the reported mean is the mean of the per-embryo values", worst < 1e-3,
              f"max |Δ| {worst:.2e}")
        check(f"{k}: n is the number of embryos contributing",
              all(g["n"] == len(g["per"]) for g in s["genes"]))
        check(f"{k}: every embryo clears the {P['MIN_TOTAL']}-transcript floor",
              all(x["n"] >= P["MIN_TOTAL"] for g in s["genes"] for x in g["per"]))
        check(f"{k}: every gene clears {P['MIN_EMBRYOS']} embryos",
              all(g["n"] >= P["MIN_EMBRYOS"] for g in s["genes"]))
        # the bulk correction is a MEDIAN of per-gene ratios, so the median gene must sit near 0
        med = float(np.median([g["lfc"] for g in s["genes"]]))
        check(f"{k}: the median gene sits near zero — that is what the bulk correction does",
              abs(med) < 0.15, f"median lfc {med:.4f}")
    check("the bulk rule names the median, not the ratio of totals",
          "MEDIAN" in m["bulk"].upper() and "ratio of totals" in m["bulk"])

    section("the alignment null is a real null")
    for k, s in d["volcano"].items():
        n = s["null"]
        check(f"{k}: the histogram has the declared number of draws",
              sum(n["hist"]) == n["draws"] == P["N_ALIGN_DRAWS"], str(sum(n["hist"])))
        # a null that just echoed the observation would sit exactly on it
        check(f"{k}: the null is not the observation", n["median"] != s["n_called"],
              f"null median {n['median']} vs observed {s['n_called']}")
        check(f"{k}: shuffling the sides destroys most of the signal",
              n["median"] < s["n_called"], f"{n['median']} vs {s['n_called']}")
        # P must be the upper tail of the shipped histogram
        tail = sum(c for v, c in enumerate(n["hist"]) if v >= s["n_called"])
        want = (tail + 1) / (n["draws"] + 1)
        check(f"{k}: the reported P is the histogram's own upper tail",
              abs(want - n["p"]) < 1e-9, f"{want:.5f} vs {n['p']}")
        check(f"{k}: the 95th percentile sits between the median and the max",
              n["median"] <= n["p95"] <= n["max"])
        print(f"        {k}: observed {s['n_called']}, null median {n['median']:.0f}, "
              f"p95 {n['p95']:.1f}, max {n['max']}, P {n['p']:.4f}")
    check("the null's purpose is stated", "flipped at random" in m["alignment_null"])

    section("the random plane is the control, and it is not at 1.0")
    rows = {r["plane"]: r for r in d["byPlane"]}
    check("the random control is present", "random" in rows)
    R = rows["random"]
    # THE POINT OF THE CONTROL: counting noise alone produces a fold. If this read 1.0 the fold
    # would be being computed against nothing and every other bar would be meaningless.
    check("counting noise alone already produces a fold above 1", R["median"] > 1.05,
          f"random median fold {R['median']}")
    check("every plane's interval brackets its own median",
          all(r["ci_lo"] <= r["median"] <= r["ci_hi"] for r in d["byPlane"]))
    check("the selected meridional plane beats the random one",
          rows["polar18"]["median"] > R["median"],
          f"{rows['polar18']['median']} vs {R['median']}")
    check("...and the artifact says that lead is partly the selection",
          "SELECTED" in m["polar18"].upper())
    for r in d["byPlane"]:
        print(f"        {r['plane']:11s} {r['median']:.4f} [{r['ci_lo']:.4f}, {r['ci_hi']:.4f}] "
              f"over {r['n_genes']} genes")

    section("the count-matched null is subtracted, not assumed")
    pg = [r for r in d["perGene"]["rows"] if r["null"] is not None]
    check("most genes have a null", len(pg) > 0.8 * len(d["perGene"]["rows"]))
    check("excess is observed minus null",
          all(abs(r["excess"] - (r["fold"] - r["null"])) < 1e-4 for r in pg))
    check("every null fold is above 1 — sampling alone is asymmetric",
          all(r["null"] > 1 - 1e-9 for r in pg), f"min {min(r['null'] for r in pg):.4f}")
    # and it must depend on count: a low-count gene should have a LARGER null than a high-count one
    lo = [r["null"] for r in pg if r["total"] < np.percentile([x["total"] for x in pg], 25)]
    hi = [r["null"] for r in pg if r["total"] > np.percentile([x["total"] for x in pg], 75)]
    check("the null is count-matched — sparse genes get a bigger one",
          float(np.median(lo)) > float(np.median(hi)),
          f"low-count {np.median(lo):.3f} vs high-count {np.median(hi):.3f}")
    check("some genes fail to beat their own null", any(r["excess"] <= 0 for r in pg),
          str(sum(1 for r in pg if r["excess"] <= 0)))

    section("the pairing trap is present, not fixed")
    pr = d["pairing"]
    check("both normalisations are shipped", sorted(pr) == ["conc", "ratio"])
    for norm, panels in pr.items():
        for k, v in panels.items():
            st = v["stat"]
            check(f"{norm}/panel {k}: the ratio is within/allpair",
                  abs(st["ratio"] - st["within"] / st["allpair"]) < 1e-3)
            check(f"{norm}/panel {k}: both tails are probabilities",
                  0 < st["p_closer"] <= 1 and 0 < st["p_farther"] <= 1)
            # the two tails must disagree — if both were small the permutation test is broken
            check(f"{norm}/panel {k}: exactly one tail is the significant one",
                  (st["p_closer"] < 0.05) != (st["p_farther"] < 0.05),
                  f"closer {st['p_closer']}, farther {st['p_farther']}")
            check(f"{norm}/panel {k}: two points per zygote",
                  len(v["embed"]["xy"]) == 2 * st["n"])
    # THE TRAP: complements must push the halves apart. If a share panel came out below 1 the
    # features are not what the page claims they are.
    ratios = [v["stat"]["ratio"] for v in pr["ratio"].values()]
    concs = [v["stat"]["ratio"] for v in pr["conc"].values()]
    check("every SHARE panel pairs above 1 — that is the complementarity artefact",
          all(r > 1 for r in ratios), str([round(r, 2) for r in ratios]))
    check("every CONCENTRATION panel falls below 1 — the halves really do cluster",
          all(r < 1 for r in concs), str([round(r, 2) for r in concs]))
    check("...and it is the clustering tail that is significant there",
          all(v["stat"]["p_closer"] < 0.05 for v in pr["conc"].values()))
    check("...while the share panels are significant in the other direction",
          all(v["stat"]["p_farther"] < 0.05 for v in pr["ratio"].values()))
    check("the trap is documented", "complement" in m["pairing_trap"].lower()
          and "geometry" in m["pairing_trap"].lower())
    check("the embedding is labelled as PCA, not UMAP",
          "PCA" in m["embedding"] and "UMAP" in m["embedding"])
    check("...and the statistic is stated not to use it",
          "FULL feature vectors" in m["embedding"])
    print("        share " + ", ".join(f"{r:.2f}" for r in ratios) +
          "  |  conc " + ", ".join(f"{r:.2f}" for r in concs))

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("halves-"))
    check("the specifications are named",
          all(x in m["method"] for x in ("4.14", "4.15", "4.17", "4.18", "4.19")))
    check("no absolute paths in the artifact", "/Users/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
