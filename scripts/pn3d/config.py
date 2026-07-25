"""
Configuration + path handling for the 3D Pronuclear Pseudotime project.

The authoritative dataset lives OUTSIDE the repository on a large mounted volume
(~733 GB). Nothing raw is ever copied into git. `redact_path` turns any absolute
local path into a safe token so committed manifests / artifacts carry no
machine-specific or private path.

Data root resolution:
  1. env  PN3D_DATA_ROOT
  2. gitignored  pn3d_local.json  in the repo root: {"data_root": "..."}
  3. default mount  /Volumes/HW/MEFISH Labels
"""
from __future__ import annotations

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))

DEFAULT_DATA_ROOT = "/Volumes/HW/MEFISH Labels"

PN3D_DIR = os.path.join(REPO_ROOT, "calibration_data", "pn3d")
DERIVED_DIR = os.path.join(PN3D_DIR, "derived")           # gitignored pixel/preview outputs
ARTIFACTS_DIR = os.path.join(PN3D_DIR, "artifacts")       # small committed JSON/CSV
DATA_DIR = os.path.join(REPO_ROOT, "data", "pn3d")        # site-served small JSON

# Physical scale: the MERFISH acquisition xy pixel size (µm). The segmentation
# NRRDs record only the voxel-grid anisotropy (z = 7 × xy); this documented
# constant maps to microns for DISPLAY. All model features are DIMENSIONLESS
# (divided by cell radius), so they do not depend on this value.
XY_UM = 0.15
Z_OVER_XY = 7.0                                            # from NRRD space directions (1,1,7)
Z_UM = XY_UM * Z_OVER_XY

LOCAL_OVERRIDES = os.path.join(REPO_ROOT, "pn3d_local.json")


def data_root() -> str:
    env = os.environ.get("PN3D_DATA_ROOT")
    if env:
        return os.path.expanduser(env)
    if os.path.isfile(LOCAL_OVERRIDES):
        try:
            v = json.load(open(LOCAL_OVERRIDES)).get("data_root")
            if v:
                return os.path.expanduser(v)
        except Exception:                                  # noqa: BLE001
            pass
    return DEFAULT_DATA_ROOT


def data_available() -> bool:
    return os.path.isdir(data_root())


def redact_path(path: str) -> str:
    """Safe token for committed output: never an absolute local path."""
    if not path:
        return ""
    root = data_root()
    if path.startswith(root):
        rel = path[len(root):].lstrip("/\\")
        return f"<data-root>/{rel}" if rel else "<data-root>"
    if path.startswith(REPO_ROOT):
        return os.path.relpath(path, REPO_ROOT)
    base = os.path.basename(path.rstrip("/\\"))
    return f"<local>/{base}"


def ensure_dirs() -> None:
    for d in (DERIVED_DIR, ARTIFACTS_DIR, DATA_DIR):
        os.makedirs(d, exist_ok=True)
