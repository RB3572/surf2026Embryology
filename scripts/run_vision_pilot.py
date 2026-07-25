#!/usr/bin/env python3
"""
One-command reproducible run of the fixed-image pseudotime PILOT.

  python3 scripts/run_vision_pilot.py            # validate + audit + report
  python3 scripts/run_vision_pilot.py --previews # also write 2.5D preview PNGs (gitignored)

Steps, in order:
  1. validate manifests (hard-fail on schema / out-of-git / supervised-fixed errors)
  2. metadata-only TIFF audit of the intensity/label stacks
  3. (optional) stream the 2.5D projections and write preview PNGs to derived/
  4. hybrid cross-check + example-stack gating, movie inventory, encoder feasibility
  5. write the redacted report JSON the website reads

Nothing is deployed or pushed. Large data is never copied into the repo.
"""
from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from scripts.vision_pseudotime import (config, manifest, projection,  # noqa: E402
                                       report, tiff_audit)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--previews", action="store_true",
                    help="write 2.5D preview PNGs to the gitignored derived dir")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if the manifest has any error")
    a = ap.parse_args()
    config.ensure_dirs()

    print("== 1. manifest validation ==")
    rep = manifest.validate_all()
    manifest._print(rep)
    if a.strict and not rep["ok"]:
        return 1

    print("\n== 2. metadata-only TIFF audit ==")
    for sid, path, kind in tiff_audit.stack_sources():
        print(tiff_audit._fmt(tiff_audit.audit(path, declared_kind=kind)))

    if a.previews:
        print("\n== 3. streaming 2.5D projections (previews) ==")
        ip = config.resolve_source("target_intensity_example")
        lp = config.resolve_source("target_labels_example")
        if os.path.isfile(ip) and os.path.isfile(lp):
            tensor, meta = projection.build_2p5d(ip, lp, out=256, z_step=1)
            base, _ = projection.single_channel_mip_baseline(ip, out=256, z_step=1)
            for c, name in enumerate(meta["channels"]):
                projection.save_png(tensor[..., c],
                                    os.path.join(config.DERIVED_DIR, f"target01_ch{c}_{name}.png"))
            projection.save_png(base, os.path.join(config.DERIVED_DIR, "target01_mip_baseline.png"))
            if "boundary_plane" in meta:
                projection.save_png(projection.crop_resize(meta["boundary_plane"], meta["bbox"], 256),
                                    os.path.join(config.DERIVED_DIR, "target01_boundary.png"))
            print(f"   wrote previews to {config.rel_to_repo(config.DERIVED_DIR)} "
                  f"(bbox {meta['bbox']}, channels {meta['channels']})")
        else:
            print("   matched stacks not present — skipped")

    print("\n== 4-5. hybrid + encoder feasibility + report ==")
    r = report.write_report()
    H = r["hybrid_geometry_clock"]
    E = r["exploratory_encoder"]
    print(f"   wrote {config.rel_to_repo(report.REPORT_PATH)}")
    print(f"   hybrid: reused {H.get('clock_reused')}, cross-checked "
          f"{H.get('validated_zygotes')} zygotes, max|Δτ|={H.get('max_tau_error')}")
    print(f"   example stack: {H.get('example_stack', {}).get('status')} "
          f"({H.get('example_stack', {}).get('blockers')})")
    print(f"   encoder: {E['status']} — synthetic grouped MAE {E['held_out_embryo_grouped_mae']}, "
          f"refuses real data = {E['real_data_refusal'] is not None}")
    print("\nPILOT run complete. Nothing deployed or pushed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
