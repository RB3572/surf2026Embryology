#!/usr/bin/env python3
"""Checks on data/clocktx.json.gz — Transcriptome vs the Clock (figures 4.8 / 4.11 / 5.4).

The three things this analysis can get quietly wrong, and so the three things asserted hardest:

  1. A SHARE WHOSE DENOMINATOR HAS DIFFERENT SCOPE FROM ITS NUMERATOR IS NOT A SHARE. The region
     set counted for a gene and the region set summed for the total must be identical, and the
     shares in an embryo must therefore sum to 1.

  2. THE REGION SET IS DELIBERATELY NOT CYTOPLASM-ONLY. Every plane analysis on this site is; this
     one counts the pronuclei too, because it is about composition rather than a spatial split. If
     that ever silently reverted to cytoplasm-only the numbers would still look plausible.

  3. A SMALL-n P MUST NOT BE scipy's ASYMPTOTIC ONE. At n = 3 the smallest achievable two-sided P
     is 1/3, so any P below that at small n is impossible and means the exact path was skipped.
"""
import gzip
import json
import math
import os
import sys
from itertools import permutations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import embryo_stats as ES                                          # noqa: E402

ART = os.path.join(ROOT, "data", "clocktx.json.gz")

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
        sys.exit("clocktx.json.gz missing — run: python3 build_clocktx.py")
    d = json.load(gzip.open(ART, "rt"))
    m, emb, traj = d["meta"], d["embryos"], d["traj"]
    by = {e["id"]: e for e in emb}
    V = d["variants"]

    section("cohort")
    check("every zygote has a τ in [0, 1]", all(0 <= e["tau"] <= 1 for e in emb),
          str([e["tau"] for e in emb if not 0 <= e["tau"] <= 1][:3]))
    check("51 zygotes carry a τ", len(emb) == 51, str(len(emb)))
    check("every zygote has a probeset", all(e.get("probeset") for e in emb))
    check("every zygote has a positive cytoplasm volume", all(e["cyto_vol"] > 0 for e in emb))
    check("meta count matches", m["n_embryos"] == len(emb))
    check("every skipped zygote records why", all(s.get("reason") for s in m.get("skipped", [])))

    section("a share is a share")
    # the shares in an embryo must sum to 1 — that is what "the denominator has the same scope"
    # means, and it is the one arithmetic identity this whole project rests on
    tot = {}
    for g, rows in traj.items():
        for r in rows:
            tot[r["id"]] = tot.get(r["id"], 0.0) + r["share"]
    worst = max(abs(v - 1.0) for v in tot.values())
    check("every embryo's shares sum to 1", worst < 1e-4, f"worst |sum-1| = {worst:.2e}")
    check("no share is negative or above 1",
          all(0 <= r["share"] <= 1 for rows in traj.values() for r in rows))
    # 1e-7, not 1e-9: conc ships at 8 decimal places
    check("concentration is count ÷ that embryo's cytoplasm volume",
          all(abs(r["conc"] - r["n"] / by[r["id"]]["cyto_vol"]) < 1e-7
              for rows in traj.values() for r in rows))
    check("total_tx equals the summed counts",
          all(abs(sum(r["n"] for rows in traj.values() for r in rows if r["id"] == e["id"])
                  - e["total_tx"]) < 0.5 for e in emb[:6]))

    section("the region set includes the pronuclei, and excludes the polar body")
    # recompute one embryo from the scene: cytoplasm-only would be strictly smaller
    e0 = emb[0]
    sc = ES.read_scene(ES.scene_path("Zygote", e0["id"]))
    body = ES.classify_body(sc)[0]
    polar = ES.polar_label(sc)
    cyto_only = 0
    body_plus_pn = 0
    with_polar = 0
    for g, t in sc["transcripts"].items():
        s = np.asarray(t["s"], int).astype(str)
        cyto_only += int((s == str(body)).sum())
        # NOT len(s): some molecules carry segment 0 — outside every segment, i.e. outside the
        # cell — and the build is right to drop them from both numerator and denominator
        with_polar += int(np.isin(s, list(ES.seg_volumes(sc))).sum())
        body_plus_pn += int(np.isin(s, [k for k in ES.seg_volumes(sc)
                                        if k != str(polar)]).sum()) if polar else len(s)
    check("the total is bigger than cytoplasm alone — the pronuclei ARE counted",
          e0["total_tx"] > cyto_only, f"total {e0['total_tx']} vs cytoplasm {cyto_only}")
    check("the total matches cytoplasm + pronuclei exactly",
          e0["total_tx"] == body_plus_pn, f"{e0['total_tx']} vs {body_plus_pn}")
    check("the polar-body variant is bigger still, and counts only in-segment molecules",
          e0["total_tx_polar"] == with_polar and with_polar >= e0["total_tx"],
          f"{e0['total_tx_polar']} vs {with_polar}")
    check("molecules outside every segment are excluded",
          e0["total_tx_polar"] < sum(len(t["s"]) for t in sc["transcripts"].values()),
          "nothing was outside the cell, which would be surprising")
    check("meta says so", "pronuclei" in m["regions"] and "polar body is excluded" in m["regions"])

    section("the four variants")
    check("all four exist", set(V) == {"main.share", "main.conc",
                                       "withPolar.share", "withPolar.conc"}, str(sorted(V)))
    for k, v in V.items():
        g = v["genes"]
        check(f"{k}: sorted by P with untestable genes last",
              all((g[i]["p"] is None) >= (g[i - 1]["p"] is None) for i in range(1, len(g)))
              and all(g[i]["p"] >= g[i - 1]["p"] - 1e-12 for i in range(1, len(g))
                      if g[i]["p"] is not None and g[i - 1]["p"] is not None))
        check(f"{k}: every ρ is in [-1, 1]",
              all(r["rho"] is None or -1.0001 <= r["rho"] <= 1.0001 for r in g))
        check(f"{k}: every P is in (0, 1]",
              all(r["p"] is None or 0 < r["p"] <= 1.0000001 for r in g),
              str([r["g"] for r in g if r["p"] is not None and not 0 < r["p"] <= 1][:3]))
        check(f"{k}: a gene under the zygote floor has no ρ",
              all(r["rho"] is None for r in g if r["n"] < m["params"]["MIN_ZYGOTES"]))
        check(f"{k}: BH q is never below the raw P",
              all(r["fdr"] is None or r["fdr"] >= r["p"] - 1e-12 for r in g if r["p"] is not None))

    section("small samples get an EXACT P")
    # at n the smallest achievable two-sided permutation P is 2/n! (and 1/3 at n=3), so anything
    # below that is scipy's asymptotic formula leaking through
    bad = []
    for k, v in V.items():
        for r in v["genes"]:
            if r["p"] is None or r["n"] >= m["params"]["EXACT_BELOW_N"]:
                continue
            floor = 2.0 / math.factorial(r["n"])
            if r["p"] < floor - 1e-12:
                bad.append((k, r["g"], r["n"], r["p"], floor))
    check("no small-n P is below what permutation can produce", not bad, str(bad[:3]))
    check("no P is exactly zero",
          not [r["g"] for v in V.values() for r in v["genes"] if r["p"] == 0])

    section("the correlation is reproducible")
    # rebuild rho for a few genes straight from the shipped trajectories
    from scipy import stats
    tau = {e["id"]: e["tau"] for e in emb}
    ps = {e["id"]: e["probeset"] for e in emb}
    worst_rho = 0.0
    tested = 0
    for r in V["main.share"]["genes"][:25]:
        if r["p"] is None:
            continue
        rows = [x for x in traj[r["g"]] if x["n"] >= m["params"]["MIN_COUNT"]]
        if len(rows) != r["n"]:
            continue
        y = np.array([x["share"] for x in rows], float)
        for p in {ps[x["id"]] for x in rows}:
            msk = np.array([ps[x["id"]] == p for x in rows])
            y[msk] -= y[msk].mean()
        t = np.array([tau[x["id"]] for x in rows], float)
        rho = float(stats.spearmanr(y, t).statistic)
        worst_rho = max(worst_rho, abs(rho - r["rho"]))
        tested += 1
    check("recomputed ρ matches the artifact", worst_rho < 1e-6, f"max |diff| {worst_rho:.2e}")
    check("enough genes rebuilt", tested >= 10, str(tested))

    section("the exemplars the reference names")
    for g in m["exemplars"]:
        check(f"{g} has a trajectory", g in traj)
    check("all four exemplars are testable",
          all(any(r["g"] == g and r["p"] is not None for r in V["main.share"]["genes"])
              for g in m["exemplars"]))

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("clocktx-"))
    check("the specification is named", "f4_8" in m.get("method", ""))
    check("the centring rule is recorded", "probeset" in m["centring"])
    check("the small-n rule is recorded", "permutation" in m["smallN"])
    check("no absolute paths", "/Users/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
