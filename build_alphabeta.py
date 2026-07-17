#!/usr/bin/env python3
"""
Sperm alpha/beta analysis for 2-cell embryos.

Recreates the mentor's "alpha/beta assignments across methods" grid (one row per
assignment method, one column per embryo, each cell = which ORIGINAL blastomere that
method calls alpha: A = segmentation labels 1+3, B = labels 2+4). This build produces
the two fully-defined methods:

  * "Sperm entry"          — the blastomere the sperm entered is alpha (seg 1/3 -> A, 2/4 -> B).
  * "Higher total transcript" — the blastomere with more transcripts is alpha (tie -> T).

More methods (PCA, ratio-sum, exhaustive, panel) will be added once the definitions are
in hand; the front-end + schema already support extra rows.

Per-blastomere transcript counts come from the already-built segments data
(public/data/segments/<Stage>__<id>.json.gz, which carries each transcript's segment).
Sperm blastomere comes from data/merfish_sperm.csv (the labelled sperm's `segment`).

Output: public/data/alphabeta.json   {stages:{early:[...], late:[...]}, methods:[...]}
"""
import csv
import glob
import gzip
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SEG_DIR = os.path.join(HERE, "public", "data", "segments")
SPERM_CSV = os.path.join(HERE, "..", "data", "merfish_sperm.csv")
OUT = os.path.join(HERE, "public", "data", "alphabeta.json")

# blastomere A = segmentation labels {1,3}; blastomere B = {2,4}. label 5 = extra body (ignored).
A_LABELS, B_LABELS = {1, 3}, {2, 4}
STAGES = [("early", "Early2Cell"), ("late", "Late2Cell")]

# methods produced here (order = row order in the grid). Each has an id + a display label.
METHODS = [
    {"id": "sperm", "label": "Sperm entry"},
    {"id": "total", "label": "Higher total transcript"},
]


def blastomere_of_segment(seg):
    """Which original blastomere a segment label belongs to ('A', 'B', or None)."""
    if seg in A_LABELS:
        return "A"
    if seg in B_LABELS:
        return "B"
    return None


def sperm_blastomere_by_embryo():
    """embryo_id -> 'A'/'B' for each sperm-positive 2-cell embryo (from the labelled segment)."""
    out = {}
    if not os.path.exists(SPERM_CSV):
        return out
    for r in csv.DictReader(open(SPERM_CSV, newline="")):
        eid = (r.get("resolved_embryo_id") or "").strip()
        md = r.get("mefish_dir", "")
        if not eid or not ("e2c" in md or "l2c" in md):
            continue
        if not r.get("x_px", "").strip() or r.get("no_sperm", "").strip().lower() == "yes":
            continue
        try:
            seg = int(float(r.get("segment", "")))
        except (TypeError, ValueError):
            continue
        bl = blastomere_of_segment(seg)
        if bl:
            out[eid] = {"blast": bl, "seg": seg}
    return out


def counts_per_blastomere(seg_file):
    """(nA, nB) transcript counts in blastomere A (1+3) and B (2+4) from a segments json.gz."""
    d = json.load(gzip.open(seg_file, "rt"))
    nA = nB = 0
    for t in d.get("transcripts", {}).values():
        for s in t.get("s", []):
            if s in A_LABELS:
                nA += 1
            elif s in B_LABELS:
                nB += 1
    return nA, nB, d.get("mask_labels", [])


def main():
    sperm = sperm_blastomere_by_embryo()
    stages_out = {}
    n_total = 0
    for key, stage_dir in STAGES:
        rows = []
        for f in sorted(glob.glob(os.path.join(SEG_DIR, f"{stage_dir}__*.json.gz"))):
            eid = os.path.basename(f).split("__", 1)[1][: -len(".json.gz")]
            nA, nB, mask = counts_per_blastomere(f)
            if not ({1, 2, 3, 4} <= set(mask)):
                continue  # need both blastomeres to compare
            # method: higher total transcript
            total_call = "A" if nA > nB else ("B" if nB > nA else "T")
            # method: sperm entry
            sp = sperm.get(eid)
            sperm_call = sp["blast"] if sp else "NA"
            rows.append({
                "id": eid,
                "calls": {"sperm": sperm_call, "total": total_call},
                "raw": {"nA": nA, "nB": nB,
                        "spermSeg": (sp["seg"] if sp else None)},
            })
        rows.sort(key=lambda r: r["id"])
        stages_out[key] = rows
        n_total += len(rows)
        n_sp = sum(1 for r in rows if r["calls"]["sperm"] != "NA")
        print(f"  {key:5s} 2-cell: {len(rows)} embryos ({n_sp} with a labelled sperm)")

    payload = {
        "methods": METHODS,
        "legend": {"A": "alpha = original A (labels 1+3)",
                   "B": "alpha = original B (labels 2+4)",
                   "T": "tie", "NA": "unavailable"},
        "stages": stages_out,
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(payload, f, separators=(",", ":"))
    print(f"wrote {OUT}  ({n_total} embryos, {len(METHODS)} methods)")


if __name__ == "__main__":
    main()
