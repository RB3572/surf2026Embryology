#!/usr/bin/env python3
"""Checks on data/animalveg.json.gz — Animal–Vegetal Enrichment (figures 4.3 and 4.4).

The result this page ships is a NULL result, and a null result has a failure mode that a
positive one does not: almost any bug makes it look more interesting than it is. So the checks
below are aimed at the things that would manufacture a signal.

  · THE PLANE IS AN EQUAL-CYTOPLASMIC-VOLUME SPLIT, not the plane through the centre of mass.
    The pronuclei and polar body are a median ~6% of the cell and sit off-centre; charging their
    volume to one side is enough to invent an asymmetry. Both halves must hold the same volume
    to floating-point precision, and the slide that achieves it must be non-zero.
  · ANIMAL IS THE POLAR-BODY SIDE by construction, so the sign of a log2 FC has a fixed meaning.
  · THE FOLD IS A DENSITY RATIO, recomputable from the per-zygote counts and volumes it ships.
  · THE HEADLINE MUST BE HONEST — the significant count has to be compared against 0.05 x the
    number of genes tested, and the BH q has to be present and monotone with P.
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

ART = os.path.join(ROOT, "data", "animalveg.json.gz")

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
        sys.exit("animalveg.json.gz missing — run: python3 build_animalveg.py")
    d = json.load(gzip.open(ART, "rt"))
    m, genes, emb, maps = d["meta"], d["genes"], d["embryos"], d["maps"]
    P = m["params"]

    section("shape")
    check("embryo count matches", m["n_embryos"] == len(emb))
    check("gene count matches", m["n_genes"] == len(genes))
    check("map count matches", m["n_maps"] == len(maps))
    check("no duplicate genes", len({g["g"] for g in genes}) == len(genes))
    check("every skipped embryo records why", all(s.get("reason") for s in m.get("skipped", [])))
    check("ranks are 1..n in P order",
          [g["rank"] for g in genes] == list(range(1, len(genes) + 1))
          and all(genes[i]["p"] <= genes[i + 1]["p"] + 1e-12 for i in range(len(genes) - 1)))

    section("the plane is an equal-CYTOPLASMIC-volume split")
    # the whole point of sliding the plane: both halves must hold the same volume once the
    # pronuclei and polar body are gone. A COM plane would leave a visible imbalance here.
    frac = [e["vAn"] / (e["vAn"] + e["vVeg"]) for e in emb]
    check("every zygote splits at 0.5 to 4 decimals",
          all(abs(f - 0.5) < 1e-4 for f in frac), f"worst {max(abs(f-0.5) for f in frac):.2e}")
    check("the two halves sum to the cytoplasm volume",
          all(abs(e["vAn"] + e["vVeg"] - e["cyto_vol"]) < 0.5 for e in emb))
    shifts = [abs(e["shift_um"]) for e in emb]
    check("the slide is real — the COM plane is NOT the equal-volume plane",
          float(np.median(shifts)) > 0.05, f"median |shift| {np.median(shifts):.3f} µm")
    check("but the slide is a correction, not a relocation", max(shifts) < 8.0,
          f"max {max(shifts):.2f} µm")
    check("the artifact says which plane it used",
          "EQUAL-CYTOPLASMIC-VOLUME" in m["plane"].upper())
    check("the artifact says which side is animal",
          "polar" in m["animal"] and "construction" in m["animal"])

    section("the fold is a density ratio, and it is recomputable")
    vol = {e["id"]: (e["vAn"], e["vVeg"]) for e in emb}
    worst_lfc = worst_mean = 0.0
    for g in genes:
        for r in g["per"]:
            vAn, vVeg = vol[r["id"]]
            want = math.log2(((r["an"] + ES.EPS) / vAn) / ((r["veg"] + ES.EPS) / vVeg))
            worst_lfc = max(worst_lfc, abs(want - r["lfc"]))
        worst_mean = max(worst_mean,
                         abs(float(np.mean([r["lfc"] for r in g["per"]])) - g["lfc"]))
    check("every per-zygote log2 FC is the density ratio it claims to be", worst_lfc < 1e-3,
          f"max |Δ| {worst_lfc:.2e}")
    check("the reported mean is the mean of those", worst_mean < 1e-3, f"max |Δ| {worst_mean:.2e}")
    check("animal + vegetal is the zygote's total for that gene",
          all(r["an"] + r["veg"] == r["n"] for g in genes for r in g["per"]))
    check("n is the number of zygotes contributing", all(g["n"] == len(g["per"]) for g in genes))
    check("total is the sum over those zygotes",
          all(g["total"] == sum(r["n"] for r in g["per"]) for g in genes))

    section("the floors are enforced")
    check(f"every tested gene clears {P['MIN_EMBRYOS']} zygotes",
          all(g["n"] >= P["MIN_EMBRYOS"] for g in genes))
    check(f"every tested gene clears {P['MIN_TOTAL']} transcripts",
          all(g["total"] >= P["MIN_TOTAL"] for g in genes))
    check("the calling rule is P and effect size together",
          all(g["called"] == (g["p"] < P["CALL_P"] and abs(g["lfc"]) >= P["CALL_LFC"])
              for g in genes))
    check("side follows the sign of the fold",
          all(g["side"] == ("animal" if g["lfc"] > 0 else "vegetal") for g in genes))

    section("the headline is honest")
    n_sig = sum(1 for g in genes if g["p"] < P["CALL_P"])
    n_called = sum(1 for g in genes if g["called"])
    check("the significant count is the real count", m["n_significant"] == n_sig, str(n_sig))
    check("the called count is the real count", m["n_called"] == n_called, str(n_called))
    check("the chance expectation is 0.05 x genes tested",
          abs(m["expected_by_chance"] - len(genes) * P["CALL_P"]) < 0.05,
          f"{m['expected_by_chance']} vs {len(genes) * P['CALL_P']:.1f}")
    check("the meta says the P-values are unadjusted", "UNADJUSTED" in m["unadjusted"])
    # THE RESULT ITSELF: this is a null result, and the page is built to say so. If a future
    # rebuild ever clears this, the page's framing is wrong and must be revisited.
    check("this is still a null result — significance does not beat chance",
          n_sig <= m["expected_by_chance"] * 1.5,
          f"{n_sig} significant vs {m['expected_by_chance']} expected")
    q = np.array([g["q"] for g in genes])
    check("BH q is present on every gene", np.isfinite(q).all())
    check("q never falls below its own P", all(g["q"] >= g["p"] - 1e-12 for g in genes))
    check("q is monotone down the ranking", bool(np.all(np.diff(q) >= -1e-12)))
    check("nothing survives correction", float(q.min()) >= 0.05, f"min q {q.min():.3f}")

    section("the density maps are relative to the panel, not pictures of the cell")
    check("every map clears the pooled-count floor",
          all(v["n"] >= P["MIN_MAP_COUNT"] for v in maps.values()))
    check("every map has the declared grid",
          all(len(v["z"]) == P["NA"] and all(len(r) == P["NR"] for r in v["z"])
              for v in maps.values()))
    # log2 vs the all-gene background: a map of RAW density would be positive almost everywhere
    # (it would just trace the cell), so the values must straddle zero
    fin = [x for v in maps.values() for r in v["z"] for x in r if x is not None]
    check("map values straddle zero — they are a ratio to the panel",
          min(fin) < -0.3 and max(fin) > 0.3, f"{min(fin):.2f}..{max(fin):.2f}")
    med = float(np.median(fin))
    check("and they are centred near zero", abs(med) < 1.0, f"median {med:.3f}")
    check("empty bins are shipped as null, not as zero",
          any(x is None for v in maps.values() for r in v["z"] for x in r))
    check("the map's meaning is recorded", "log2" in m["map"] and "all-gene" in m["map"])
    # a map can exist for a gene the volcano did not test (the map floor is pooled counts, the
    # test floor is per-embryo detections) but never the reverse for a well-detected gene
    extra = sorted(set(maps) - {g["g"] for g in genes})
    check("mapped genes not in the volcano are rare and explained by the different floors",
          len(extra) < 0.05 * len(maps), f"{len(extra)}: {extra[:5]}")
    check("most tested genes have a map", len(maps) > 0.8 * len(genes),
          f"{len(maps)}/{len(genes)}")

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("animalveg-"))
    check("the specification is named", "4.3" in m["method"] and "4.4" in m["method"])
    check("the pairing argument is stated", "paired" in m["paired"].lower())
    check("no absolute paths in the artifact", "/Users/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
