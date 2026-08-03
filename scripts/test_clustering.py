#!/usr/bin/env python3
"""Checks on data/clustering.json.gz — the spatial gene clustering artifact.

Guards the things that would be invisible on the page but wrong: signatures that secretly
encode how well-measured a gene is rather than where it sits, cluster labels that disagree
with the member lists, and the sparse-bin blow-up that a flat pseudocount would reintroduce.
"""
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ART = os.path.join(ROOT, "data", "clustering.json.gz")

passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        print(f"  PASS  {name}")
    else:
        failed += 1
        print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))


def section(t):
    print(f"\n[{t}]")


def main():
    if not os.path.isfile(ART):
        sys.exit("clustering.json.gz missing — run: python3 build_clustering.py")
    d = json.load(gzip.open(ART, "rt"))
    m, genes = d["meta"], d["genes"]
    names = [g["gene"] for g in genes]
    X = np.array([g["profile"] for g in genes], float)
    dim = m["n_rad"] + m["n_ax"]

    section("shape")
    check("has genes", len(genes) > 50, f"{len(genes)}")
    check("meta n_genes matches", m["n_genes"] == len(genes))
    check("gene names unique", len(set(names)) == len(names))
    check(f"every signature is {dim}-D", X.shape[1] == dim and np.isfinite(X).all())
    check("genes sorted (stable ordering for the page)", names == sorted(names))

    section("signature sanity")
    # A flat pseudocount let empty bins reach log2(eps/f) ~ -7; shrinkage must keep them sane.
    check("no sparse-bin blow-up (|log2| < 4 everywhere)", np.abs(X).max() < 4.0,
          f"max |value| = {np.abs(X).max():.2f}")
    check("signatures are not all identical", X.std(axis=0).max() > 0.02)
    check("some gene is enriched somewhere (positive values exist)", (X > 0.05).any())
    check("some gene is depleted somewhere (negative values exist)", (X < -0.05).any())

    section("no measurement-depth confound")
    # If clusters simply tracked how many embryos a gene appears in, this would be a
    # coverage map wearing a biology costume.
    n_emb = np.array([g["n_emb"] for g in genes], float)
    check(f"every gene clears min_emb={m['min_emb']}", (n_emb >= m["min_emb"]).all())
    amp = np.abs(X).mean(axis=1)
    r = float(np.corrcoef(n_emb, amp)[0, 1])
    check("signature strength not driven by embryo count (|r| < 0.5)", abs(r) < 0.5, f"r = {r:+.3f}")

    section("clusterings")
    for k in m["k_range"]:
        kd = d["k"][str(k)]
        lab = kd["labels"]
        cl = kd["clusters"]
        ok_len = len(lab) == len(genes)
        ok_k = len(set(lab)) == k == len(cl)
        # the member lists must agree with the labels, or the right drawer lies
        by_label = {}
        for i, l in enumerate(lab):
            by_label.setdefault(l, set()).add(names[i])
        ok_mem = all(set(cl[c]["members"]) == by_label.get(c, set()) for c in range(len(cl)))
        ok_n = all(cl[c]["n"] == len(cl[c]["members"]) for c in range(len(cl)))
        ok_prof = all(len(cl[c]["profile"]) == dim for c in range(len(cl)))
        check(f"k={k}: labels/clusters/members consistent",
              ok_len and ok_k and ok_mem and ok_n and ok_prof)

    section("layouts")
    for key in ("mds", "tsne"):
        E = np.array([g[key] for g in genes], float)
        check(f"{key} is 2-D, finite, normalised to ~[-1,1]",
              E.shape == (len(genes), 2) and np.isfinite(E).all() and np.abs(E).max() <= 1.0001,
              f"max |coord| = {np.abs(E).max():.3f}")
        check(f"{key} is not degenerate", E.std(axis=0).min() > 1e-3)

    section("defaults + provenance")
    check("default_k is in k_range", m["default_k"] in m["k_range"])
    best = max(m["k_range"], key=lambda k: m["silhouette"][str(k)])
    check("default_k is the best silhouette", m["default_k"] == best,
          f"default {m['default_k']} vs best {best}")
    check("silhouettes are positive (real structure)",
          all(v > 0 for v in m["silhouette"].values()))
    # Edges ship rounded to 4 dp, so cubing them drifts ~1e-4. 1e-3 absorbs that while still
    # catching a genuinely wrong scheme: equal-WIDTH edges would be off by 0.375.
    check("radial edges are equal-volume",
          np.allclose([e ** 3 for e in m["rad_edges"]],
                      np.linspace(0, 1, m["n_rad"] + 1), atol=1e-3))
    check("version recorded", str(m.get("version", "")).startswith("clustering-"))
    check("no absolute paths in the artifact",
          "/Users/" not in json.dumps(d) and "\\Users\\" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
