#!/usr/bin/env python3
"""Checks on data/contacthalves.json.gz — the contact region on the reference's own definition.

The strongest check is at the bottom: the reference ships its own 7.2 table, and this build is
compared against it gene by gene. It is not an exact-match check and should not be — this build
places the splitting plane on the mesh, so folds land within about 0.02 rather than on top of each
other. What must match exactly is the FAMILY: every gene the reference tests is tested here and no
gene it excludes is included, because the per-embryo detection floor was recovered from that table
rather than guessed.

Above it, the things that would quietly manufacture a contact bias in EVERY gene:

  · THE SPLIT IS EQUAL-VOLUME, per blastomere. If it drifted to the centroid plane instead, the
    junction-side half would be systematically the smaller one and every gene would lean.
  · THE EMBRYO IS CENTRED ON ITS OWN BULK SPLIT before averaging, so a marginally misplaced plane
    cannot push a whole embryo's genes the same way.
  · THE PROFILE FALLS OFF AT BOTH ENDS. That is cell shape, and if it did not appear the frame
    would not be the frame the analysis claims.
  · THE MAPS ARE RATIOS TO THE ALL-GENE MAP, so they must straddle zero rather than trace the cell.
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

ART = os.path.join(ROOT, "data", "contacthalves.json.gz")
KEYS = ["early2cell", "late2cell"]

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
        sys.exit("contacthalves.json.gz missing — run: python3 build_contacthalves.py")
    d = json.load(gzip.open(ART, "rt"))
    m, P = d["meta"], d["meta"]["params"]

    section("shape")
    check("both stages present", sorted(d["stages"]) == sorted(KEYS))
    for k in KEYS:
        s = d["stages"][k]
        check(f"{k}: no duplicate genes", len({g['g'] for g in s['genes']}) == len(s["genes"]))
        check(f"{k}: ranks are 1..n in P order",
              [g["rank"] for g in s["genes"]] == list(range(1, len(s["genes"]) + 1))
              and all(s["genes"][i]["p"] <= s["genes"][i + 1]["p"] + 1e-12
                      for i in range(len(s["genes"]) - 1)))
        check(f"{k}: embryo count matches meta", m["n_embryos"][k] == len(s["embryos"]))
    check("every skipped embryo records why", all(x.get("reason") for x in m.get("skipped", [])))

    section("the split is equal-volume, per blastomere")
    # both halves of BOTH blastomeres, summed: if the plane drifted to the centroid this would
    # sit off 0.5 and every gene would inherit the imbalance
    for k in KEYS:
        fr = [e["vIn"] / (e["vIn"] + e["vOut"]) for e in d["stages"][k]["embryos"]]
        check(f"{k}: every embryo splits at 0.5 to 4 decimals",
              all(abs(f - 0.5) < 1e-4 for f in fr), f"worst {max(abs(f-0.5) for f in fr):.2e}")
    check("the artifact says the plane is slid, not centred",
          "equal volume" in m["plane"].lower() and "centroid" in m["plane"].lower())
    check("the artifact says nuclei and the polar body are excluded",
          "nuclei" in m["excluded"] and "polar body" in m["excluded"])

    section("the fold is recomputable, and it is bulk-centred")
    for k in KEYS:
        s = d["stages"][k]
        bulk = {e["id"]: e["bulk"] for e in s["embryos"]}
        worst = worst_mean = 0.0
        for g in s["genes"]:
            for r in g["per"]:
                want = math.log2((r["c"] + ES.EPS) / (r["e"] + ES.EPS)) - bulk[r["id"]]
                worst = max(worst, abs(want - r["lfc"]))
            worst_mean = max(worst_mean,
                             abs(float(np.mean([r["lfc"] for r in g["per"]])) - g["lfc"]))
        check(f"{k}: every per-embryo value is the bulk-centred log2 ratio", worst < 2e-3,
              f"max |Δ| {worst:.2e}")
        check(f"{k}: the reported mean is the mean of those", worst_mean < 1e-3)
        check(f"{k}: n is the number of embryos contributing",
              all(g["n"] == len(g["per"]) for g in s["genes"]))
        check(f"{k}: total is the sum over those embryos",
              all(g["total"] == sum(r["c"] + r["e"] for r in g["per"]) for g in s["genes"]))
    # the centring must actually do something: bulk splits are not all zero
    b = [abs(e["bulk"]) for k in KEYS for e in d["stages"][k]["embryos"]]
    check("the bulk correction bites — embryos are not already balanced",
          float(np.median(b)) > 0.005, f"median |bulk| {np.median(b):.4f}")
    check("the centring rule is recorded", "bulk" in m["centring"])

    section("the floors are enforced")
    for k in KEYS:
        s = d["stages"][k]
        check(f"{k}: every embryo clears the {P['MIN_TX']}-transcript floor",
              all(r["c"] + r["e"] >= P["MIN_TX"] for g in s["genes"] for r in g["per"]))
        check(f"{k}: every gene clears {P['MIN_EMBRYOS']} embryos",
              all(g["n"] >= P["MIN_EMBRYOS"] for g in s["genes"]))
        check(f"{k}: side follows the sign of the fold",
              all(g["side"] == ("contact" if g["lfc"] > 0 else "edge") for g in s["genes"]))

    section("the headline is the count beside its own chance expectation")
    for k in KEYS:
        s = d["stages"][k]
        check(f"{k}: the nominal count is the real count",
              s["n_nominal"] == sum(1 for g in s["genes"] if g["p"] < 0.05), str(s["n_nominal"]))
        check(f"{k}: the expectation is 0.05 x genes tested",
              abs(s["expected"] - 0.05 * len(s["genes"])) < 0.05)
        check(f"{k}: the FDR count is the real count",
              s["n_fdr"] == sum(1 for g in s["genes"] if g.get("q", 1) < 0.05))
        q = np.array([g["q"] for g in s["genes"]])
        check(f"{k}: q is monotone down the ranking", bool(np.all(np.diff(q) >= -1e-12)))
        check(f"{k}: q never falls below its own P",
              all(g["q"] >= g["p"] - 1e-12 for g in s["genes"]))
        print(f"        {k}: {s['n_nominal']} nominal vs {s['expected']} expected, "
              f"{s['n_fdr']} FDR")

    section("the profile is the frame, and it shows cell shape")
    for k in KEYS:
        rows = [p for p in d["profile"] if p["stage"] == k]
        check(f"{k}: one profile per embryo", len(rows) == len(d["stages"][k]["embryos"]))
        check(f"{k}: every profile has {P['NBIN']} bins", all(len(p["f"]) == P["NBIN"] for p in rows))
        check(f"{k}: every profile sums to 1", all(abs(sum(p["f"]) - 1) < 1e-3 for p in rows))
        mu = np.mean([p["f"] for p in rows], axis=0)
        # THE FALL-OFF: the middle bins must hold more than the two end bins. That is the cell
        # shape the volcano then divides out; without it the frame is not what it claims to be.
        check(f"{k}: the profile falls off at both ends (cell shape, not localisation)",
              mu[len(mu) // 2] > mu[0] and mu[len(mu) // 2] > mu[-1],
              f"ends {mu[0]:.4f}/{mu[-1]:.4f} vs middle {mu[len(mu)//2]:.4f}")
        check(f"{k}: and it is not flat", float(mu.max() / mu.min()) > 1.5,
              f"max/min {mu.max()/mu.min():.2f}")
    check("the profile's normalisation is recorded",
          "own axial reach" in m["profile"].lower() and "cell shape" in m["profile"].lower())

    section("the density maps are ratios, not pictures of the cell")
    maps = d["maps"]
    check("one map per direction", sorted(maps) == ["contact", "edge"])
    for k, M in maps.items():
        check(f"{k}: the declared grid", len(M["z"]) == P["NMAP"]
              and all(len(r) == P["NMAP"] for r in M["z"]))
        check(f"{k}: it names the genes it pooled", len(M["genes"]) == P["TOP_MAP"])
        check(f"{k}: enough pooled transcripts to draw", M["enough"], str(M["n"]))
    fin = [x for M in maps.values() for r in M["z"] for x in r if x is not None]
    check("map values straddle zero — they are a ratio to the all-gene map",
          min(fin) < -0.3 and max(fin) > 0.3, f"{min(fin):.2f}..{max(fin):.2f}")
    check("and they are centred near zero", abs(float(np.median(fin))) < 0.6,
          f"median {np.median(fin):.3f}")
    check("low-coverage bins are masked out, which is what gives the stepped outline",
          any(x is None for M in maps.values() for r in M["z"] for x in r))
    # the sets were chosen BY their position on this axis, so they must lean the way they were
    # chosen — this is a self-consistency check on the frame, not evidence of anything
    e2c = d["stages"]["early2cell"]["genes"]
    top = {g["g"] for g in sorted(e2c, key=lambda r: -r["lfc"])[:P["TOP_MAP"]]}
    check("the contact map's gene set is the top of the volcano by fold",
          set(maps["contact"]["genes"]) == top)
    check("the maps are labelled illustrative, not a test",
          "NOT A TEST" in m["maps"].upper())

    section("what is deliberately absent is said so")
    check("the missing GO panel is explained", "7.5" in m["no_go"] and "annotation" in m["no_go"])

    section("VALIDATION against the reference's own 7.2 table")
    v = m.get("validation") or {}
    if not v.get("available"):
        check("the reference table was read", False, "not present on this machine")
    else:
        for k in KEYS:
            s = v["stages"][k]
            # THE FAMILY MUST MATCH EXACTLY. The per-embryo floor was recovered from this table,
            # so any drift here means the floor is wrong and every P is computed over the wrong
            # family — which is what the chance expectation is divided by.
            check(f"{k}: exactly the reference's gene family, no more and no fewer",
                  s["n_shared"] == s["n_mine"] == s["n_ref"],
                  f"shared {s['n_shared']}, mine {s['n_mine']}, ref {s['n_ref']}")
            check(f"{k}: the folds track the reference", s["r"] > 0.95, f"r {s['r']}")
            check(f"{k}: and land close to it", s["median_abs_diff"] < 0.03,
                  f"median |Δ| {s['median_abs_diff']}")
            check(f"{k}: they mostly agree in direction", s["sign_agree"] > 0.85,
                  f"{s['sign_agree']}")
            print(f"        {k}: {s['n_shared']} genes, r = {s['r']}, "
                  f"median |Δ| {s['median_abs_diff']}, sign agree {s['sign_agree']}")

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("contacthalves-"))
    check("the specification is named",
          all(x in m["method"] for x in ("7.1", "7.2", "7.3")))
    check("no absolute paths in the artifact", "/Users/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
