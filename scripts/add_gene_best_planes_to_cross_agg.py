#!/usr/bin/env python3
"""
Enrich the cross-embryo aggregates with each gene's OWN best plane and p-value.

The aggregate's `best` / `sig` are TRANSCRIPT-WEIGHTED OVER ALL GENES — they answer
"where does this zygote's whole transcriptome split best?". The aligned-outlines
figure asks a different question: "where does THIS GENE split best?". Using the
all-gene plane there rotates every zygote by an angle chosen from other genes.

The per-embryo scenes already carry the exact per-gene answers (analysis.genes[]
.bestP_vol / bestP_cnt / bestDiff_vol / bestDiff_cnt, and the permutation p-values
in .planes[].pVol / .pCnt), so this only lifts them into the aggregate — the numbers
are the same ones the rest of the site reports, not a client-side approximation.

Adds, per embryo:

    gb[gene] = [ bp,  bp,  bp,  bp,      # best plane index, one per BEST_KEY
                 a,   a,   a,   a,       # side-A count at that plane (flip side)
                 p,   p,   p,   p ]      # permutation p at that plane (outline colour)

ordered by BEST_KEYS = [pVol, pCnt, diffVol, diffCnt]. The p paired with a Δ-mode is
that mode's matching normalization (diffVol→pVol, diffCnt→pCnt), exactly as the
embryo-level `sig` already does. Everything else is left untouched.

Run from the deploy repo root:  python3 scripts/add_gene_best_planes_to_cross_agg.py
"""
import gzip
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")

BEST_KEYS = ["pVol", "pCnt", "diffVol", "diffCnt"]
# per-gene best-plane field in the scene, and which permutation p to report there
BEST_FIELD = {"pVol": "bestP_vol", "pCnt": "bestP_cnt",
              "diffVol": "bestDiff_vol", "diffCnt": "bestDiff_cnt"}
SIG_FIELD = {"pVol": "pVol", "pCnt": "pCnt", "diffVol": "pVol", "diffCnt": "pCnt"}

# (aggregate, scene dir, which analysis). Division Plane Sweep sweeps 18 candidate
# planes; the sperm-plane project has a single geometrically-fixed plane, so its
# per-gene best plane is trivially 0 — but its per-gene p still matters for colour.
TARGETS = [
    (os.path.join(DATA, "zygote_cross.json.gz"), os.path.join(DATA, "zygote"), "real"),
    (os.path.join(DATA, "zygote_cross_circ.json.gz"), os.path.join(DATA, "zygote"), "circ"),
    (os.path.join(DATA, "sperm_division_cross.json.gz"),
     os.path.join(DATA, "sperm_division"), "real"),
    (os.path.join(DATA, "sperm_division_cross_circ.json.gz"),
     os.path.join(DATA, "sperm_division"), "circ"),
]


def scene_gene_best(scene, which):
    """{gene: [4 best planes, 4 side-A counts, 4 p-values]} for the real or circ analysis."""
    analysis = scene.get("analysis") if which == "real" else (scene.get("circ") or {}).get("analysis")
    if not analysis:
        return None
    out = {}
    for row in analysis.get("genes", []):
        planes = row.get("planes") or []
        if not planes:
            continue
        bps, aas, pps = [], [], []
        for key in BEST_KEYS:
            bp = int(row[BEST_FIELD[key]])
            bp = max(0, min(bp, len(planes) - 1))
            bps.append(bp)
            aas.append(int(planes[bp].get("a", 0)))
            pps.append(round(float(planes[bp].get(SIG_FIELD[key], 1.0)), 5))
        out[row["gene"]] = bps + aas + pps
    return out


def main():
    for path, scenes, which in TARGETS:
        if not os.path.exists(path):
            print(f"  -- {os.path.basename(path)} missing, skipped")
            continue
        agg = json.load(gzip.open(path, "rt"))
        before = os.path.getsize(path)
        n_ok = n_moved = n_pairs = 0
        for emb in agg.get("embryos", []):
            sp = os.path.join(scenes, emb["id"] + ".json.gz")
            if not os.path.exists(sp):
                print(f"  !! no scene for {emb['id']}")
                continue
            best = scene_gene_best(json.load(gzip.open(sp, "rt")), which)
            if best is None:
                print(f"  !! no {which} analysis for {emb['id']}")
                continue
            # keep only the genes this aggregate already carries, so the two stay in sync
            gb = {g: best[g] for g in emb.get("g", {}) if g in best}
            # sanity: the count stored at a gene's own best plane must match the
            # scene's per-plane counts the aggregate already ships in gp
            gp = emb.get("gp") or {}
            for g, row in gb.items():
                per = gp.get(g)
                if per:
                    for i in range(4):
                        assert row[4 + i] == per[row[i]], f"{emb['id']}/{g}: a mismatch at plane {row[i]}"
                n_pairs += 1
                n_moved += sum(1 for i in range(4) if row[i] != emb["best"][i])
            emb["gb"] = gb
            n_ok += 1
        with gzip.open(path, "wt") as fh:
            json.dump(agg, fh, separators=(",", ":"))
        after = os.path.getsize(path)
        pct = (100.0 * n_moved / (4 * n_pairs)) if n_pairs else 0.0
        print(f"  {os.path.basename(path)}: {n_ok} embryos got gb over {n_pairs} embryo×gene pairs "
              f"· {pct:.0f}% of planes differ from the all-gene plane "
              f"· {before/1024:.0f} KB → {after/1024:.0f} KB")


if __name__ == "__main__":
    main()
