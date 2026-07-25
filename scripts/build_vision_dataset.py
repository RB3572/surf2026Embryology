#!/usr/bin/env python3
"""
Build the cached vision-model dataset: DAPI image features + targets, one row per
fixed zygote. Committed CSV is small (features + targets only, NO pixels), so the
model training and the website are self-contained and never need the raw stacks.

Targets (all from the ALREADY-VALIDATED segmentation pipeline, used only as
labels here):
  * distance_sum_um       cell-centre distance sum (the geometry clock's feature)
  * tau                   published pseudotime for the zygote
  * pron_min_distance_um  pronuclei min surface-to-surface gap (pronuclei project)
  * transcript_total      total transcripts (for the distance-vs-transcripts view)

Usage:  python3 scripts/build_vision_dataset.py
Needs the raw DAPI stacks under ../TranscriptomicsData/.../Zygote/<id>/dapi.tif.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from scripts.vision_pseudotime import vision_features as VF  # noqa: E402

Z = os.path.join(HERE, "..", "TranscriptomicsData", "JustTifAndCSVData", "Zygote")
GEOM = os.path.join(HERE, "calibration_data", "fixed_cohort_geometry.csv")
FIXED = os.path.join(HERE, "data", "pronuclei_pseudotime.json")
MANIFEST = os.path.join(HERE, "data", "pronuclei_manifest.json")
OUT = os.path.join(HERE, "calibration_data", "vision_pseudotime", "vision_features.csv")


def main() -> int:
    geom = {r["id"]: r for r in csv.DictReader(open(GEOM)) if r.get("distance_sum_um")}
    fx = {r["id"]: r for r in json.load(open(FIXED))["embryos"] if r.get("features")}
    man = {e["id"]: e for e in json.load(open(MANIFEST))["embryos"]}

    ids = sorted(set(geom) & set(fx) & set(man))
    ids = [i for i in ids if os.path.isfile(os.path.join(Z, i, "dapi.tif"))]
    print(f"{len(ids)} zygotes with DAPI + geometry + tau + manifest")

    rows = []
    for k, eid in enumerate(ids, 1):
        t0 = time.time()
        try:
            feats = VF.extract_from_path(os.path.join(Z, eid, "dapi.tif"))
        except Exception as e:                                     # noqa: BLE001
            print(f"  [{k}/{len(ids)}] {eid}: FEATURE ERROR {e}")
            continue
        g, f, m = geom[eid], fx[eid], man[eid]
        row = {"id": eid, "label": m.get("label", eid), "date_short": m.get("date_short", ""),
               "distance_sum_um": float(g["distance_sum_um"]),
               "nearer_to_center_um": float(g["nearer_to_center_um"]),
               "farther_to_center_um": float(g["farther_to_center_um"]),
               "tau": float(f["tau"]), "qc": f.get("qc", ""),
               "pron_min_distance_um": float(m["distance"]),
               "transcript_total": int(m["total"]), **feats}
        rows.append(row)
        print(f"  [{k}/{len(ids)}] {m.get('label', eid):14s} "
              f"spread={feats['spread_rms']:6.1f} dsum={row['distance_sum_um']:5.1f} "
              f"tau={row['tau']:.3f} mindist={row['pron_min_distance_um']:5.1f} ({time.time()-t0:.1f}s)")

    cols = (["id", "label", "date_short", "distance_sum_um", "nearer_to_center_um",
             "farther_to_center_um", "tau", "qc", "pron_min_distance_um", "transcript_total"]
            + list(VF.FEATURE_NAMES))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nwrote {os.path.relpath(OUT, HERE)} ({len(rows)} rows, {len(cols)} cols, "
          f"{os.path.getsize(OUT)/1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
