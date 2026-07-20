#!/usr/bin/env python3
"""
Pronuclei enrichment lists for the right drawer of the pronuclei project.

For every zygote (that has two detected pronuclei) and every gene, we assign each
transcript to a segment via the full-resolution label TIFF (the pronuclei are small,
so the downsampled mask is too coarse for the counts). Then, per (zygote, gene):

  * n_pn   = transcripts inside the two pronuclei (labels stored in the scene)
  * n_cyto = transcripts inside the cytoplasm (segment 1)
  * n_cell = n_pn + n_cyto            (transcripts within the cell)
  * V_pn, V_cell = pronuclei / (cytoplasm+pronuclei) volumes (downsampled voxels)
  * p_null = V_pn / V_cell            (chance a random transcript lands in a pronucleus)
  * fold   = (n_pn/V_pn) / (n_cell/V_cell) = pronuclei density ÷ cell-average density
  * p      = P(Binomial(n_cell, p_null) >= n_pn)   — the exact null the prompt asks for:
             scatter the same number of transcripts at random and count how many land in
             the pronuclei segment.

List 1  "Enriched in pronuclei": n_cell >= MIN_COUNT and fold >= FOLD_THRESH.  Report fold + p.
List 2  "Only in pronuclei":     n_cyto == 0 (every cell-transcript is in a pronucleus),
                                 n_cell >= MIN_COUNT_2.  Report p.

Output: public/data/pronuclei_enrichment.json.gz  {enriched:[...], onlyPn:[...], meta}
Run from embryo_viewer/:  python3 build_pronuclei_enrichment.py
"""
import glob
import gzip
import json
import os

import numpy as np
import tifffile
from scipy.stats import binom

import build_pronuclei as BP

HERE = BP.HERE
DATA = os.path.join(HERE, "public", "data")
OUT = os.path.join(DATA, "pronuclei_enrichment.json.gz")

MIN_COUNT = 10      # list 1: minimum total (in-cell) transcripts to be testable
FOLD_THRESH = 1.5   # list 1: >= 1.5x density fold increase
MIN_COUNT_2 = 5     # list 2: minimum transcripts for an "only in pronuclei" call
MIN_FRAC_2 = 0.75   # list 2: >= 75% of in-cell transcripts inside the pronuclei. The pronuclei
                    # are ~5% of the cell volume, so 75% there is a ~15x enrichment — effectively
                    # "only in the pronuclei". Strict 100%/90% is essentially never met (a few
                    # stray transcripts always land in the cytoplasm).


def transcript_labels(label_path, tx, genes):
    """Full-res segment label under each transcript, per gene (one fancy-index for all)."""
    with tifffile.TiffFile(label_path) as t:
        mm = t.series[0].asarray(out="memmap")
        Zn, Yn, Xn = mm.shape
        lens = [len(tx[g]["x"]) for g in genes]
        if sum(lens):
            gx = np.concatenate([np.asarray(tx[g]["x"], float) for g in genes])
            gy = np.concatenate([np.asarray(tx[g]["y"], float) for g in genes])
            gz = np.concatenate([np.asarray(tx[g]["gz"], float) for g in genes])
            ix = np.clip(np.round(gx).astype(np.int64), 0, Xn - 1)
            iy = np.clip(np.round(gy).astype(np.int64), 0, Yn - 1)
            iz = np.clip(np.round(gz).astype(np.int64), 0, Zn - 1)
            labs = np.asarray(mm[iz, iy, ix]).astype(np.int32)
        else:
            labs = np.empty(0, np.int32)
        del mm
    out, off = {}, 0
    for g, L in zip(genes, lens):
        out[g] = labs[off:off + L]; off += L
    return out


def process(eid, pron_labels):
    scene_p = os.path.join(BP.ATLAS, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(BP.SRC, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    d = json.load(gzip.open(scene_p, "rt"))
    tx = d.get("transcripts", {})
    genes = [g for g in tx if len(tx[g]["x"])]
    if not genes:
        return None
    labs = transcript_labels(lab[0], tx, genes)
    la, lb = int(pron_labels[0]), int(pron_labels[1])
    sub = BP.load_sub(lab[0])
    v_pn = int(((sub == la) | (sub == lb)).sum())
    v_cyto = int((sub == BP.CYTO).sum())
    v_cell = v_pn + v_cyto
    if v_pn == 0 or v_cell == 0:
        return None
    p_null = v_pn / v_cell
    rows = []
    for g in genes:
        lg = labs[g]
        npn = int(((lg == la) | (lg == lb)).sum())
        ncyto = int((lg == BP.CYTO).sum())
        ncell = npn + ncyto
        if ncell == 0:
            continue
        rows.append((g, ncell, npn, ncyto))
    return rows, p_null, v_pn, v_cell


def main():
    genes_agg = json.load(gzip.open(os.path.join(DATA, "pronuclei_genes.json.gz"), "rt"))
    man = json.load(open(os.path.join(DATA, "pronuclei_manifest.json")))
    labels = {m["id"]: m.get("label", m["id"]) for m in (man.get("embryos") or man.get("points") or man)}

    enriched, only_pn, candidates = [], [], []
    n_ok = 0
    for e in genes_agg["embryos"]:
        eid = e["id"]
        sp = os.path.join(DATA, "pronuclei", eid + ".json.gz")
        if not os.path.exists(sp):
            continue
        pron_labels = json.load(gzip.open(sp, "rt")).get("pron_labels")
        if not pron_labels or len(pron_labels) < 2:
            continue
        r = process(eid, pron_labels)
        if not r:
            print(f"  -- skipped {eid}")
            continue
        rows, p_null, v_pn, v_cell = r
        n_ok += 1
        lab = labels.get(eid, eid)
        for g, ncell, npn, ncyto in rows:
            fold = (npn / v_pn) / (ncell / v_cell) if npn else 0.0
            p = float(binom.sf(npn - 1, ncell, p_null)) if npn > 0 else 1.0
            frac = npn / ncell
            if npn == 0 or ncell < MIN_COUNT_2:
                continue
            row = {"id": eid, "label": lab, "gene": g, "n": ncell, "npn": npn,
                   "frac": round(frac, 3), "fold": round(fold, 2), "p": p, "pnFrac": round(p_null, 4)}
            # keep the pronuclei-leaning candidates so the front-end thresholds can change without
            # re-reading every TIFF (all fold>=1.5 hits, plus anything >=50% in the pronuclei)
            if fold >= FOLD_THRESH or frac >= 0.5:
                candidates.append(row)
            if ncell >= MIN_COUNT and fold >= FOLD_THRESH:
                enriched.append(row)
            if frac >= MIN_FRAC_2:
                only_pn.append(row)

    enriched.sort(key=lambda r: (-r["fold"], r["p"]))
    only_pn.sort(key=lambda r: (-r["frac"], r["p"]))
    candidates.sort(key=lambda r: (-r["frac"], r["p"]))
    fr = sorted((r["frac"] for r in candidates), reverse=True)
    print("  frac in pronuclei ≥ : " + " ".join(
        f"{t}:{sum(f >= t for f in fr)}" for t in (0.9, 0.8, 0.75, 0.7, 0.6, 0.5)))
    out = {"enriched": enriched, "onlyPn": only_pn, "candidates": candidates,
           "meta": {"minCount": MIN_COUNT, "foldThresh": FOLD_THRESH, "minCount2": MIN_COUNT_2,
                    "onlyFrac": MIN_FRAC_2, "nZygotes": n_ok,
                    "nEnriched": len(enriched), "nOnlyPn": len(only_pn)}}
    with gzip.open(OUT, "wt") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(f"  {n_ok} zygotes · {len(enriched)} enriched (fold>={FOLD_THRESH}, n>={MIN_COUNT}) · "
          f"{len(only_pn)} only-in-pronuclei (n>={MIN_COUNT_2})")
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
