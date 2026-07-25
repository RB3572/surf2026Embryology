#!/usr/bin/env python3
"""
Build the versioned data manifest for the 3D Pronuclear Pseudotime project.

Writes two files:
  * data/pn3d/manifest.json                 committed, redacted paths (site + provenance)
  * calibration_data/pn3d/derived/manifest_local.json   gitignored, resolved paths (pipeline)

Metadata-only; safe to run repeatedly. Usage: python3 scripts/build_pn3d_manifest.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from scripts.pn3d import config, manifest  # noqa: E402


def main() -> int:
    config.ensure_dirs()
    if not config.data_available():
        print(f"data root not mounted ({config.redact_path(config.data_root())}); cannot inventory")
        return 1
    m = manifest.build(include_paths=True)

    local_p = os.path.join(config.DERIVED_DIR, "manifest_local.json")
    json.dump(m, open(local_p, "w"), indent=1)              # resolved paths (gitignored)

    committed = json.loads(json.dumps(m))
    for e in committed["embryos"]:
        e.pop("_paths", None)
    committed_p = os.path.join(config.DATA_DIR, "manifest.json")
    json.dump(committed, open(committed_p, "w"), indent=1)

    print(f"experiments : {m['n_experiments']}")
    print(f"embryos     : {m['n_embryos']}")
    print(f"by stage    : {m['by_stage']}")
    print(f"by batch    : {m['by_batch']}")
    print(f"spacing seen: {m['spacing_grid_values_seen']}")
    z = [e for e in m["embryos"] if e["stage"] == "zygote"]
    print(f"zygotes     : {len(z)}  (complete: {sum(e['complete'] for e in z)}, "
          f"4-segment: {sum(e.get('n_segments') == 4 for e in z)})")
    print(f"shared-FOV embryos: {sum(e['shared_fov'] for e in m['embryos'])}")
    print(f"\nwrote data/pn3d/manifest.json (redacted) + derived/manifest_local.json (resolved)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
