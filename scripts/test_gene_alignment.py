#!/usr/bin/env python3
"""
Tests for gene-specific alignment of the cross-embryo aligned-outlines figure.

The figure asks "where does THIS GENE split this zygote best?", but the aggregate's
`best` / `sig` answer a different question — "where does this zygote's WHOLE
transcriptome split best?". Using the latter rotates each zygote by an angle chosen
from other genes and colours it by other genes' significance, and draws zygotes that
never detected the selected gene at all.

Designed to FAIL on the ways that regression could come back:
  * `gb` missing, so the viewer silently falls back to the all-gene plane;
  * `gb` disagreeing with the per-plane counts the same file already ships in `gp`;
  * a gene's stored best plane not actually being that gene's best plane;
  * the p paired with a Δ-mode using the wrong normalization;
  * non-carrier zygotes being counted as alignable.

Run: python3 scripts/test_gene_alignment.py
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")

BEST_KEYS = ["pVol", "pCnt", "diffVol", "diffCnt"]
BEST_FIELD = {"pVol": "bestP_vol", "pCnt": "bestP_cnt",
              "diffVol": "bestDiff_vol", "diffCnt": "bestDiff_cnt"}
SIG_FIELD = {"pVol": "pVol", "pCnt": "pCnt", "diffVol": "pVol", "diffCnt": "pCnt"}

# (aggregate, scene dir, sweeps multiple candidate planes?)
AGGS = [
    ("zygote_cross.json.gz", "zygote", True),
    ("zygote_cross_circ.json.gz", "zygote", True),
    ("sperm_division_cross.json.gz", "sperm_division", False),
    ("sperm_division_cross_circ.json.gz", "sperm_division", False),
]

FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def load(path):
    return json.load(gzip.open(path, "rt"))


def main():
    print("gene-specific alignment — tests\n")
    for fname, scene_dir, sweeps in AGGS:
        path = os.path.join(DATA, fname)
        print(f"[{fname}]")
        if not os.path.isfile(path):
            check(f"{fname} exists", False, "aggregate missing")
            continue
        agg = load(path)
        embs = agg["embryos"]
        keys = agg.get("best_keys", BEST_KEYS)

        # ── every embryo carries gb, or the viewer silently uses the all-gene plane ──
        missing = [e["id"] for e in embs if "gb" not in e]
        check("every embryo has gb", not missing, f"{len(missing)} without gb")

        # ── gb covers exactly the genes the aggregate already reports ──
        mismatch = [e["id"] for e in embs if set(e.get("gb", {})) != set(e.get("g", {}))]
        check("gb covers exactly the genes in g", not mismatch, f"{len(mismatch)} embryos differ")

        # ── shape: 4 planes + 4 counts + 4 p-values ──
        bad_len = [(e["id"], gn) for e in embs for gn, r in e.get("gb", {}).items() if len(r) != 12]
        check("each gb row is 4 planes + 4 counts + 4 p-values", not bad_len, str(bad_len[:2]))

        # ── planes are in range; the single-plane project must always pick plane 0 ──
        npl = agg.get("n_planes") or 1
        oob = [(e["id"], gn, r[:4]) for e in embs for gn, r in e.get("gb", {}).items()
               if any(not (0 <= p < npl) for p in r[:4])]
        check("best planes are within range", not oob, str(oob[:2]))
        if not sweeps:
            nz = [(e["id"], gn) for e in embs for gn, r in e.get("gb", {}).items() if set(r[:4]) != {0}]
            check("single-plane project always aligns to plane 0", not nz, str(nz[:2]))

        # ── the stored side-A count must be that gene's count AT the stored plane ──
        bad_a, checked_a = [], 0
        for e in embs:
            gp = e.get("gp") or {}
            for gn, r in e.get("gb", {}).items():
                per = gp.get(gn)
                if not per:
                    continue
                for i in range(4):
                    checked_a += 1
                    if r[4 + i] != per[r[i]]:
                        bad_a.append((e["id"], gn, keys[i]))
        check(f"side-A count matches gp at the stored plane ({checked_a} checks)",
              not bad_a, str(bad_a[:2]))

        # ── counts never exceed the gene total, p-values are probabilities ──
        bad_n = [(e["id"], gn) for e in embs for gn, r in e.get("gb", {}).items()
                 if any(not (0 <= a <= e["g"][gn][0]) for a in r[4:8])]
        check("side-A counts lie within [0, gene total]", not bad_n, str(bad_n[:2]))
        # The build's two-sided permutation p is 2 x an uncapped one-sided tail
        # (build_zygote.perm_pvals), so a thoroughly non-significant split legitimately
        # scores above 1 and tops out at 2 — the same convention the embryo-level `sig`
        # already uses, and what the plot legend's "n.s. (p >= 1)" refers to. The bound
        # to enforce is [0, 2]; anything outside that is a real corruption.
        bad_p = [(e["id"], gn, r[8:]) for e in embs for gn, r in e.get("gb", {}).items()
                 if any(not (0.0 <= p <= 2.0) for p in r[8:])]
        check("p-values lie in [0, 2] (2x one-sided permutation tail)", not bad_p, str(bad_p[:2]))
        sig_max = max((v for e in embs for v in e.get("sig", {}).values()), default=0)
        check("gb p-values share the same convention as the existing sig",
              sig_max <= 2.0, f"sig max {sig_max}")

        # ── the per-gene plane must genuinely differ from the all-gene plane ──
        # (if it never did, the whole point of gb would be moot and a regression
        #  reverting to `best` would pass every other test here)
        if sweeps:
            diff = sum(1 for e in embs for gn, r in e.get("gb", {}).items()
                       for i in range(4) if r[i] != e["best"][i])
            tot = sum(4 for e in embs for _ in e.get("gb", {}))
            frac = diff / tot if tot else 0
            check("per-gene planes differ from the all-gene plane", frac > 0.5,
                  f"only {frac:.0%} differ — gb may have been rebuilt from `best`")

        # ── gb must agree with the per-embryo scenes it was derived from ──
        scenes = os.path.join(DATA, scene_dir)
        which = "circ" if "circ" in fname else "real"
        bad_scene, n_scene = [], 0
        for e in embs[:6]:                                  # sample: full sweep is slow
            sp = os.path.join(scenes, e["id"] + ".json.gz")
            if not os.path.isfile(sp):
                continue
            sc = load(sp)
            A = sc.get("analysis") if which == "real" else (sc.get("circ") or {}).get("analysis")
            if not A:
                continue
            rows = {r["gene"]: r for r in A.get("genes", [])}
            for gn, r in e.get("gb", {}).items():
                row = rows.get(gn)
                if not row:
                    continue
                n_scene += 1
                for i, k in enumerate(keys):
                    bp = int(row[BEST_FIELD[k]])
                    if r[i] != bp or r[4 + i] != int(row["planes"][bp]["a"]):
                        bad_scene.append((e["id"], gn, k))
                    elif abs(r[8 + i] - float(row["planes"][bp][SIG_FIELD[k]])) > 1e-4:
                        bad_scene.append((e["id"], gn, k + " p"))
        check(f"gb matches the source scenes ({n_scene} embryo×gene pairs)",
              not bad_scene, str(bad_scene[:3]))

        # ── a Δ-mode's p must use its matching normalization, not the other one ──
        # diffVol pairs with pVol and diffCnt with pCnt, exactly as `sig` already does.
        ivol, icnt = keys.index("diffVol"), keys.index("diffCnt")
        ipv, ipc = keys.index("pVol"), keys.index("pCnt")
        crossed = 0
        for e in embs:
            for gn, r in e.get("gb", {}).items():
                # when a Δ-mode lands on the same plane as its p-mode, the p must match
                if r[ivol] == r[ipv] and abs(r[8 + ivol] - r[8 + ipv]) > 1e-9:
                    crossed += 1
                if r[icnt] == r[ipc] and abs(r[8 + icnt] - r[8 + ipc]) > 1e-9:
                    crossed += 1
        check("Δ-modes report their matching normalization's p", crossed == 0,
              f"{crossed} rows disagree")

        # ── carriers: a zygote with no counts for a gene must not be alignable ──
        ghosts = [(e["id"], gn) for e in embs for gn, arr in e.get("g", {}).items() if not arr[0]]
        check("no zero-count genes recorded as carriers", not ghosts, str(ghosts[:2]))
        cov = agg.get("gene_cov") or {}
        if cov:
            gn0 = max(cov, key=lambda k: cov[k])
            real = sum(1 for e in embs if e.get("g", {}).get(gn0, [0])[0])
            check(f"gene_cov matches carrier count for {gn0}", cov[gn0] == real,
                  f"{cov[gn0]} vs {real}")
            check("no gene is present in every embryo (so filtering always matters)",
                  max(cov.values()) < len(embs), f"max cov {max(cov.values())}/{len(embs)}")
        print()

    # ── the viewers must actually consume gb ──
    print("[viewer wiring]")
    for js in ("zygote.js", "sperm-division.js"):
        src = open(os.path.join(HERE, js)).read()
        check(f"{js} reads gb", "e.gb" in src)
        check(f"{js} excludes non-carriers", "if (!ga) { nSkipped++; return; }" in src
              or "if (!gs) { nSkipped++; return; }" in src)
        # p can reach 2, so the colour ramp must clamp rather than run off its scale
        check(f"{js} clamps the significance colour at p >= 1",
              "Math.min(1, (Math.log10(Math.max(p, 1e-12)) + 3) / 3)" in src)

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"  FAILED: {f}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
