#!/usr/bin/env python3
"""Checks on data/contact.json.gz — blastomere contact-region enrichment.

Guards the assumptions the page cannot show you: that k is nested in n, that the expected
fraction rises with slab thickness, that both 2-cell stages are present, and that the result
is not simply "the slab is where all the transcripts are".
"""
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ART = os.path.join(ROOT, "data", "contact.json.gz")

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
        sys.exit("contact.json.gz missing — run: python3 build_contact.py")
    d = json.load(gzip.open(ART, "rt"))
    m, emb, genes = d["meta"], d["embryos"], d["genes"]
    R = m["radii"]

    section("shape")
    check("has embryos", len(emb) > 20, f"{len(emb)}")
    check("meta counts match", m["n_embryos"] == len(emb))
    check("both stages present", m["n_e2c"] > 0 and m["n_l2c"] > 0, f"e2c={m['n_e2c']} l2c={m['n_l2c']}")
    check("stage counts add up", m["n_e2c"] + m["n_l2c"] == len(emb))
    check("gene universe non-empty", len(genes) > 50, f"{len(genes)}")
    check("radii ascending", all(R[i] < R[i + 1] for i in range(len(R) - 1)))
    check("default radius is one of them", m["default_radius"] in R)

    section("per-embryo geometry")
    bad_axis = [e["id"] for e in emb if abs(np.linalg.norm(e["axis_um"]) - 1) > 1e-3]
    check("interface axis is a unit vector", not bad_axis, str(bad_axis[:3]))
    mids = [np.allclose(np.array(e["mid_um"]),
                        (np.array(e["com_a_um"]) + np.array(e["com_b_um"])) / 2, atol=1e-2) for e in emb]
    check("midpoint is between the two blastomere centres", all(mids))
    check("blastomere separation is physical (10-120 um)",
          all(10 <= e["sep_um"] <= 120 for e in emb),
          str(sorted(e["sep_um"] for e in emb)[:2]))
    check("every embryo names two distinct blastomere labels",
          all(e["a"] != e["b"] for e in emb))
    check("scene file recorded for every embryo", all(e.get("scene") for e in emb))

    section("counts are self-consistent")
    bad = []
    for e in emb:
        for g, r in e["genes"].items():
            k = r["k"]
            if len(k) != len(R): bad.append((e["id"], g, "len")); break
            if any(v < 0 or v > r["n"] for v in k): bad.append((e["id"], g, "k>n")); break
            if any(k[i] > k[i + 1] for i in range(len(k) - 1)): bad.append((e["id"], g, "not monotone")); break
    check("k is within [0,n] and rises with slab thickness", not bad, str(bad[:3]))
    check("every gene clears min_tx",
          all(r["n"] >= m["min_tx"] for e in emb for r in e["genes"].values()))

    f0 = np.array([e["f0"] for e in emb])
    check("f0 in (0,1]", ((f0 > 0) & (f0 <= 1)).all())
    check("f0 rises with slab thickness", all((f0[:, i] <= f0[:, i + 1] + 1e-9).all() for i in range(f0.shape[1] - 1)))
    di = R.index(m["default_radius"])
    check("the default slab is a MINORITY of the embryo (a real sub-region)",
          f0[:, di].mean() < 0.25, f"mean f0 = {f0[:, di].mean():.3f}")

    section("there is signal, and it is not one-sided")
    # pooled fold per gene at the default slab — a working analysis should find both
    # enriched and depleted genes, not push everything one way
    acc = {}
    for e in emb:
        f = e["f0"][di]
        for g, r in e["genes"].items():
            a = acc.setdefault(g, [0, 0, 0.0])
            a[0] += r["k"][di]; a[1] += r["n"]; a[2] += r["n"] * f
    folds = [a[0] / a[2] for a in acc.values() if a[1] >= 300 and a[2] > 0]
    folds = np.array(folds)
    check("enough genes to compare", len(folds) > 30, f"{len(folds)}")
    check("some genes enriched (fold > 1.15)", (folds > 1.15).any(), f"max {folds.max():.2f}")
    check("some genes depleted (fold < 0.85)", (folds < 0.85).any(), f"min {folds.min():.2f}")
    check("folds centre near 1 (the null is calibrated)",
          0.85 < float(np.median(folds)) < 1.15, f"median {np.median(folds):.3f}")

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("contact-"))
    check("no absolute paths in the artifact",
          "/Users/" not in json.dumps(d) and "\\Users\\" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
