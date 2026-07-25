#!/usr/bin/env python3
"""
Run the structure-aware segmentation audit over the mounted dataset and cache the
dimensionless geometry per embryo (task 16 -> feeds 17/18/19).

Audits all zygotes (the pronuclear target) and a sample of 2-cell embryos (used
only as post-NEBD / out-of-distribution reference — never as pronuclear labels).
Writes a committed, redacted JSON; no pixels, no absolute paths.

Usage: python3 scripts/build_pn3d_geometry.py [--all-2cell]
"""
from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from scripts.pn3d import config, segment_audit as SA  # noqa: E402

OUT = os.path.join(config.DATA_DIR, "segmentation_geometry.json")


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-2cell", action="store_true", help="audit every 2-cell too (slow)")
    ap.add_argument("--n-2cell", type=int, default=10, help="how many 2-cell embryos to sample")
    a = ap.parse_args()
    config.ensure_dirs()

    man = json.load(open(os.path.join(config.DERIVED_DIR, "manifest_local.json")))
    zyg = [e for e in man["embryos"] if e["stage"] == "zygote" and e["complete"]]
    two = [e for e in man["embryos"] if e["stage"] in ("e2c", "l2c") and e["complete"]]
    if not a.all_2cell:
        two = two[:a.n_2cell]
    todo = zyg + two
    print(f"auditing {len(zyg)} zygotes + {len(two)} 2-cell (OOD reference)")

    records = []
    for k, e in enumerate(todo, 1):
        lp = e["_paths"]["label_tifs"][0]
        t0 = time.time()
        try:
            au = SA.audit_path(lp)
        except Exception as ex:                                    # noqa: BLE001
            print(f"  [{k}/{len(todo)}] {e['embryo_id']}: ERROR {ex}")
            records.append({"embryo_id": e["embryo_id"], "stage": e["stage"],
                            "status": "error", "reason": str(ex)})
            continue
        rec = {
            "embryo_id": e["embryo_id"], "experiment": e["experiment"],
            "batch_date": e["batch_date"], "stage": e["stage"], "plate": e["plate"],
            "fov": e["fov"], "sub_index": e["sub_index"], "shared_fov": e["shared_fov"],
            "n_segments": e.get("n_segments"), "spacing_grid": e.get("spacing_grid"),
            "status": au["status"], "confidence": au["confidence"], "flags": au["flags"],
            "cell_body_label": au.get("cell_body_label"),
            "pronucleus_labels": au.get("pronucleus_labels"),
            "polar_body_label": au.get("polar_body_label"),
            "geometry": au.get("geometry"),
            "is_pronuclear_target": e["stage"] == "zygote",
        }
        records.append(rec)
        g = au.get("geometry")
        tag = f"rms/R={g['rms_over_R']} nPN={g['n_pronuclei']}" if g else au["status"]
        print(f"  [{k}/{len(todo)}] {e['embryo_id']:40s} {au['status']:10s} "
              f"conf={au['confidence']:.2f} {tag} ({time.time()-t0:.1f}s)")

    zres = [r for r in records if r["stage"] == "zygote" and r.get("geometry")]
    payload = {
        "schema_version": 1, "package_version": "pn3d-0.1.0",
        "generated_at_utc": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc).replace(microsecond=0).isoformat(),
        "audit_params": {"z_step": SA.Z_STEP, "xy_step": SA.XY_STEP,
                         "voxel_um": [round(float(x), 3) for x in SA.voxel_um()],
                         "inside_rule": "fill-containment >= 0.6 OR border-ratio > 3",
                         "semantics_source": "geometric audit — NOT segment/file names"},
        "n_zygotes_audited": sum(r["stage"] == "zygote" for r in records),
        "n_zygotes_resolved": len(zres),
        "n_2cell_ood": sum(r["stage"] in ("e2c", "l2c") for r in records),
        "embryos": records,
    }
    json.dump(payload, open(OUT, "w"), indent=1)
    print(f"\nscored {len(zres)}/{sum(r['stage']=='zygote' for r in records)} zygotes (geometry computed)")
    print(f"wrote {os.path.relpath(OUT, HERE)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
