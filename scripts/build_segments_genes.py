#!/usr/bin/env python3
"""
Cross-embryo gene aggregate for the segment-enrichment project.

The per-embryo scenes (data/segments/<id>.json.gz) each carry ranked[seg] =
[{gene, enrich, count, ntot}, ...] — a gene's fold enrichment and count in each
segment of THAT embryo. This pools them so the bottom drawer can show a selected
gene's enrichment across EVERY embryo it appears in (density fold or count fraction).

Output: data/segments_genes.json.gz
  { embInfo: [{id, label, stage, segs:[labels]}],           # index -> embryo
    genes:   { gene: [ [embIdx, [[seg, enrich, count, ntot], ...]], ... ] } }

Run:  python3 scripts/build_segments_genes.py   (from the deploy repo root)
"""
import glob
import gzip
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCENES = os.path.join(ROOT, "data", "segments")
MAN = os.path.join(ROOT, "data", "segments_manifest.json")
OUT = os.path.join(ROOT, "data", "segments_genes.json.gz")


def main():
    man = json.load(open(MAN))
    emb_info, genes = [], {}
    for e in man["embryos"]:
        sp = os.path.join(SCENES, e["id"] + ".json.gz")
        if not os.path.exists(sp):
            print(f"  -- no scene for {e['id']}")
            continue
        d = json.load(gzip.open(sp, "rt"))
        idx = len(emb_info)
        emb_info.append({"id": e["id"], "label": e.get("label", e["id"]),
                         "stage": e.get("stage_label", e.get("stage", "")),
                         "segs": e.get("segments", d.get("mask_labels", []))})
        per_gene = {}   # gene -> {seg: [enrich, count, ntot]}
        for seg in d.get("mask_labels", []):
            for r in d.get("ranked", {}).get(str(seg), []):
                per_gene.setdefault(r["gene"], {})[int(seg)] = [round(r["enrich"], 3), r["count"], r["ntot"]]
        for g, segvals in per_gene.items():
            rows = sorted([[s] + v for s, v in segvals.items()])
            genes.setdefault(g, []).append([idx, rows])

    out = {"embInfo": emb_info, "genes": genes}
    with gzip.open(OUT, "wt") as fh:
        json.dump(out, fh, separators=(",", ":"))
    ng = len(genes)
    cov = sorted((len(v) for v in genes.values()), reverse=True)
    print(f"  {len(emb_info)} embryos · {ng} genes · widest coverage {cov[0] if cov else 0} embryos · "
          f"{os.path.getsize(OUT)/1024:.0f} KB")


if __name__ == "__main__":
    main()
