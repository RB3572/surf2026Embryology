"""
Metadata-only TIFF audit (brief item 2).

Reads TIFF *tags* and series geometry WITHOUT loading pixel data, so a 2.5 GB
volume is inspected in milliseconds. Reports shape/dtype/axes, whether the stack
is an intensity or an integer-label stack, and — crucially — whether physical
voxel spacing and channel/class identity are trustworthy from the file itself.

The example stacks were written by InsightToolkit with no ImageJ/OME metadata,
a bogus inch-based XResolution, and (as always for baseline TIFF) no Z spacing.
So a sidecar is required. This module states that requirement instead of
silently inventing a scale.

Sidecar JSON schema (see sidecars/*.example.json):
  {
    "embryo_id": "str",
    "voxel_um":  [x, y, z],          # physical size of one voxel, microns
    "channels":  ["DAPI", ...],      # channel identity, in order (intensity stacks)
    "label_classes": {"1": "cytoplasm", "2": "pronucleus"},  # label stacks
    "acquisition_time": null,        # or ISO string; fixed snapshots have none
    "treatment": "untreated",
    "notes": "..."
  }
"""
from __future__ import annotations

import json
import math
import os

import tifffile

from . import config

# Only embryo_id is strictly required. Physical spacing is often not yet known
# for these lab stacks; the sidecar records that honestly with voxel_um=null and
# voxel_um_verified=false rather than inventing a scale.
REQUIRED_SIDECAR_KEYS = ("embryo_id",)
# ResolutionUnit: 1=none, 2=inch, 3=cm. Microscopy spacing is never inches.
_TRUSTWORTHY_UNITS = {3}                       # only centimetre-based is even plausible


def _rational(v):
    try:
        if isinstance(v, tuple) and len(v) == 2 and v[1]:
            return v[0] / v[1]
    except Exception:                                              # noqa: BLE001
        pass
    return None


def sidecar_path_for(tiff_path: str) -> str:
    base = os.path.splitext(os.path.basename(tiff_path))[0]
    return os.path.join(config.SIDECAR_DIR, base + ".sidecar.json")


def load_sidecar(tiff_path: str, explicit: str | None = None) -> dict | None:
    p = explicit or sidecar_path_for(tiff_path)
    if os.path.isfile(p):
        try:
            return json.load(open(p))
        except Exception as e:                                     # noqa: BLE001
            return {"_error": f"unreadable sidecar {config.rel_to_repo(p)}: {e}"}
    return None


def validate_sidecar(sc: dict) -> list[str]:
    errs = []
    if sc is None:
        return ["no sidecar found"]
    if "_error" in sc:
        return [sc["_error"]]
    for k in REQUIRED_SIDECAR_KEYS:
        if k not in sc:
            errs.append(f"sidecar missing required key {k!r}")
    vox = sc.get("voxel_um")
    if vox is not None:                                            # null is allowed (unknown scale)
        if not (isinstance(vox, (list, tuple)) and len(vox) == 3
                and all(isinstance(v, (int, float)) and v > 0 for v in vox)):
            errs.append("voxel_um, if present, must be three positive numbers [x, y, z]")
    return errs


def spacing_status(sc: dict | None) -> str:
    """'verified' | 'provisional' | 'missing' — never invents a scale."""
    if not sc or "_error" in sc:
        return "missing"
    if sc.get("voxel_um") is None:
        return "missing"
    return "verified" if sc.get("voxel_um_verified") is True else "provisional"


def audit(tiff_path: str, sidecar: str | None = None, declared_kind: str | None = None) -> dict:
    """Inspect a TIFF using tags only. Never reads pixels.

    `declared_kind` is the manifest source_kind ('confocal_stack' /
    'segmentation_stack'). Metadata alone CANNOT distinguish an intensity stack
    from a label stack — both are integer dtypes here — so intensity-vs-label is
    taken from the manifest/sidecar, never guessed from dtype.
    """
    out: dict = {"path_redacted": config.redact_path(tiff_path),
                 "basename": os.path.basename(tiff_path),
                 "exists": os.path.isfile(tiff_path)}
    if not out["exists"]:
        out["error"] = "file not present on this machine"
        return out

    out["size_bytes"] = os.path.getsize(tiff_path)
    with tifffile.TiffFile(tiff_path) as tf:
        s = tf.series[0]
        shape, dtype, axes = tuple(int(x) for x in s.shape), str(s.dtype), s.axes
        out["shape"] = shape
        out["axes"] = axes
        out["dtype"] = dtype
        out["n_pages"] = len(s.pages)
        out["dtype_is_integer"] = dtype in (
            "int8", "int16", "int32", "uint8", "uint16", "uint32")
        # z, y, x for a plain grayscale stack (axes like IYX / ZYX)
        zyx = None
        if len(shape) == 3 and axes in ("IYX", "ZYX", "QYX"):
            zyx = {"z": shape[0], "y": shape[1], "x": shape[2]}
        out["zyx"] = zyx
        itemsize = int(str(s.dtype).replace("uint", "").replace("int", "").replace("float", "")) // 8
        out["est_full_load_bytes"] = int(math.prod(shape)) * max(itemsize, 1)

        pg = s.pages[0]
        tags = {t.name: t.value for t in pg.tags}
        out["software"] = tags.get("Software")
        out["resolution_unit"] = tags.get("ResolutionUnit")
        xr, yr = _rational(tags.get("XResolution")), _rational(tags.get("YResolution"))
        out["xresolution_px_per_unit"] = round(xr, 6) if xr else None
        out["yresolution_px_per_unit"] = round(yr, 6) if yr else None
        out["has_imagej_metadata"] = tf.imagej_metadata is not None
        out["has_ome_metadata"] = getattr(tf, "ome_metadata", None) is not None
        # ImageJ 'spacing' is the only baseline place a z-step could live
        ij_spacing = (tf.imagej_metadata or {}).get("spacing") if tf.imagej_metadata else None
        out["imagej_z_spacing"] = ij_spacing

    # spacing trustworthiness: need real xy scale AND a z step
    unit_ok = out["resolution_unit"] in _TRUSTWORTHY_UNITS
    xy_ok = unit_ok and bool(out["xresolution_px_per_unit"]) and out["xresolution_px_per_unit"] > 100
    z_ok = out["imagej_z_spacing"] is not None or out["has_ome_metadata"]
    out["spacing_trustworthy"] = bool(xy_ok and z_ok)
    out["needs_sidecar"] = not out["spacing_trustworthy"]

    sc = load_sidecar(tiff_path, sidecar)
    out["sidecar_present"] = sc is not None and "_error" not in (sc or {})
    out["sidecar_errors"] = validate_sidecar(sc) if out["needs_sidecar"] else []
    if out["sidecar_present"]:
        out["voxel_um"] = sc.get("voxel_um")
        out["channels"] = sc.get("channels")
        out["label_classes"] = sc.get("label_classes")
        out["embryo_id"] = sc.get("embryo_id")
        out["treatment"] = sc.get("treatment")
        out["acquisition_time"] = sc.get("acquisition_time")
    # intensity vs label comes from the manifest/sidecar, NOT from dtype
    has_label_classes = bool(out.get("label_classes"))
    if declared_kind == "segmentation_stack" or has_label_classes:
        out["stack_role"] = "labels"
    elif declared_kind == "confocal_stack" or out.get("channels"):
        out["stack_role"] = "intensity"
    else:
        out["stack_role"] = "unknown"
    out["spacing_status"] = spacing_status(sc)
    out["voxel_um_verified"] = out["spacing_status"] == "verified"
    # physical-scale geometry (µm distances -> the µm clock) is only sound with a
    # VERIFIED voxel size. A provisional sidecar lets projections/segmentation run
    # in voxel space, clearly labelled, but does not unlock physical-unit claims.
    out["metadata_sufficient"] = bool(
        out["spacing_trustworthy"] or (out["sidecar_present"] and out["voxel_um_verified"]))
    out["projection_possible"] = out["sidecar_present"] and not out["sidecar_errors"]
    return out


def _fmt(a: dict) -> str:
    if not a.get("exists"):
        return f"  {a['basename']}: {a.get('error')}"
    gb = a["est_full_load_bytes"] / 1e9
    lines = [
        f"  {a['basename']}",
        f"    shape={a['shape']} axes={a['axes']} dtype={a['dtype']}  "
        f"(full load ≈ {gb:.2f} GB — never done)",
        f"    stack_role={a.get('stack_role')}  dtype_is_integer={a['dtype_is_integer']}  "
        f"software={a['software']}",
        f"    spacing_trustworthy={a['spacing_trustworthy']}  needs_sidecar={a['needs_sidecar']}  "
        f"sidecar_present={a['sidecar_present']}  metadata_sufficient={a['metadata_sufficient']}",
    ]
    if a.get("sidecar_errors"):
        for e in a["sidecar_errors"]:
            lines.append(f"    sidecar: {e}")
    if a.get("voxel_um"):
        lines.append(f"    voxel_um={a['voxel_um']} channels={a.get('channels')} "
                     f"classes={a.get('label_classes')}")
    return "\n".join(lines)


def stack_sources() -> list[tuple[str, str, str]]:
    """[(source_id, resolved_path, source_kind)] for the intensity/label stacks."""
    import csv
    out = []
    with open(config.SOURCES_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("source_kind") in ("confocal_stack", "segmentation_stack"):
                out.append((r["source_id"],
                            config.resolve_source_path(r["source_id"], r.get("location", "")),
                            r["source_kind"]))
    return out


if __name__ == "__main__":
    import sys
    for sid, path, kind in stack_sources():
        print(_fmt(audit(path, declared_kind=kind)))
    sys.exit(0)
