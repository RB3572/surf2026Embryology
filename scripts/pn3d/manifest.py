"""
Versioned data manifest + inventory for the mounted MEFISH dataset (task 15).

Metadata-only: walks the volume, identifies every embryo (a FOV or a nested
sub-embryo directory carrying image + segmentation), and records provenance,
batch, stage, channels, segmentation, voxel spacing (from the NRRD header),
completeness, and duplicate / multi-instance fields of view. No large file is
loaded. Committed manifest carries only redacted paths.
"""
from __future__ import annotations

import glob
import json
import os
import re
from datetime import datetime, timezone

from . import SCHEMA_VERSION, __version__, config

STAGE_RE = re.compile(r"(zygote|e2c|l2c)")
DATE_RE = re.compile(r"^(\d{8})")


def _nrrd_header(path: str) -> dict:
    """Parse a NRRD text header (spacing, sizes, segment names/extents). Tiny read."""
    out: dict = {}
    try:
        with open(path, "rb") as f:
            buf = b""
            while b"\n\n" not in buf and len(buf) < 65536:
                chunk = f.read(4096)
                if not chunk:
                    break
                buf += chunk
        txt = buf.split(b"\n\n")[0].decode("latin-1")
    except Exception:                                              # noqa: BLE001
        return out
    seg_names, seg_labels = {}, {}
    for line in txt.splitlines():
        if line.startswith("sizes:"):
            out["sizes"] = [int(x) for x in line.split(":", 1)[1].split()]
        elif line.startswith("space directions:"):
            dirs = re.findall(r"\(([^)]*)\)", line)
            vals = []
            for d in dirs:
                comps = [float(c) for c in d.split(",")]
                vals.append(max(abs(c) for c in comps))       # magnitude per axis
            out["spacing_grid"] = vals                        # (x, y, z) grid steps
        elif line.startswith("space:"):
            out["space"] = line.split(":", 1)[1].strip()
        else:
            m = re.match(r"Segment(\d+)_Name:=(.*)", line)
            if m:
                seg_names[int(m.group(1))] = m.group(2).strip()
            m = re.match(r"Segment(\d+)_LabelValue:=(\d+)", line)
            if m:
                seg_labels[int(m.group(1))] = int(m.group(2))
    if seg_names:
        out["n_segments"] = len(seg_names)
        out["segment_names"] = [seg_names.get(i, "") for i in sorted(seg_names)]
        out["segment_names_all_generic"] = all(
            re.fullmatch(r"Segment_\d+", n or "") for n in seg_names.values())
        out["segment_label_values"] = [seg_labels.get(i) for i in sorted(seg_labels)]
    return out


def _is_embryo_dir(d: str) -> bool:
    has_img = any(os.path.isfile(os.path.join(d, n)) for n in
                  ("all.tif", "all.tiff", "dapi.tif", "dapi_1.tif"))
    has_seg = bool(glob.glob(os.path.join(d, "Segmentation-label*.tif")))
    return has_img and has_seg


def _first(d: str, *names) -> str:
    for n in names:
        p = os.path.join(d, n)
        if os.path.isfile(p):
            return p
    # glob fallback
    for n in names:
        g = glob.glob(os.path.join(d, n))
        if g:
            return g[0]
    return ""


def _embryo_record(exp: str, date: str, stage: str, plate: str,
                   fov: str, sub: str, d: str) -> dict:
    dna = _first(d, "dapi.tif", "dapi_1.tif")
    allc = _first(d, "all.tif", "all.tiff", "all_1.tiff", "all_1.tif")
    labels = sorted(glob.glob(os.path.join(d, "Segmentation-label*.tif")))
    seg_nrrd = _first(d, "Segmentation.seg.nrrd")
    scene = _first(d, "*.mrml")
    ctbl = _first(d, "Segmentation-label_ColorTable.ctbl", "Segmentation-label_ColorTable.txt")
    hdr = _nrrd_header(seg_nrrd) if seg_nrrd else {}

    eid = "/".join([exp, fov] + ([sub] if sub else []))
    completeness = {
        "dna": bool(dna), "all": bool(allc), "label_tif": bool(labels),
        "seg_nrrd": bool(seg_nrrd), "scene": bool(scene), "color_table": bool(ctbl),
    }
    return {
        "embryo_id": eid, "experiment": exp, "batch_date": date, "stage": stage,
        "plate": plate, "fov": fov, "sub_index": sub or None,
        "n_label_instances": len(labels),          # >1 = multiple segmented cells in one FOV
        "spacing_grid": hdr.get("spacing_grid"),    # (x,y,z) voxel-grid steps; z≈7×xy
        "sizes_xyz": hdr.get("sizes"),
        "n_segments": hdr.get("n_segments"),
        "segment_names_all_generic": hdr.get("segment_names_all_generic"),
        "space": hdr.get("space"),
        "completeness": completeness,
        "complete": all(completeness[k] for k in ("dna", "all", "label_tif", "seg_nrrd")),
        "paths_redacted": {
            "dna": config.redact_path(dna), "all": config.redact_path(allc),
            "label_tifs": [config.redact_path(p) for p in labels],
            "seg_nrrd": config.redact_path(seg_nrrd), "scene": config.redact_path(scene),
        },
        "_paths": {"dna": dna, "all": allc, "label_tifs": labels,
                   "seg_nrrd": seg_nrrd, "scene": scene},   # resolved; stripped from committed copy
    }


def build(include_paths: bool = False) -> dict:
    root = config.data_root()
    if not os.path.isdir(root):
        raise SystemExit(f"data root not mounted: {config.redact_path(root)}")

    embryos = []
    experiments = []
    for exp in sorted(os.listdir(root)):
        ed = os.path.join(root, exp)
        if not os.path.isdir(ed):
            continue
        date = (DATE_RE.match(exp) or [None, ""])[1] if DATE_RE.match(exp) else ""
        stage_m = STAGE_RE.search(exp)
        stage = stage_m.group(1) if stage_m else "unknown"
        plate_m = re.search(r"_(p\d+)", exp)
        plate = plate_m.group(1) if plate_m else ""
        exp_embryos = 0
        for fov in sorted(os.listdir(ed)):
            fd = os.path.join(ed, fov)
            if not os.path.isdir(fd):
                continue
            if _is_embryo_dir(fd):
                embryos.append(_embryo_record(exp, date, stage, plate, fov, "", fd))
                exp_embryos += 1
            # nested sub-embryo dirs (multiple cells / re-segmentations in a FOV)
            for sub in sorted(os.listdir(fd)):
                sd = os.path.join(fd, sub)
                if os.path.isdir(sd) and _is_embryo_dir(sd):
                    embryos.append(_embryo_record(exp, date, stage, plate, fov, sub, sd))
                    exp_embryos += 1
        experiments.append({"experiment": exp, "batch_date": date, "stage": stage,
                            "plate": plate, "n_embryos": exp_embryos})

    # duplicate / shared-FOV detection: same experiment+fov with >1 embryo record
    from collections import Counter
    fov_counts = Counter((e["experiment"], e["fov"]) for e in embryos)
    for e in embryos:
        e["shared_fov"] = fov_counts[(e["experiment"], e["fov"])] > 1

    by_stage = Counter(e["stage"] for e in embryos)
    by_batch = Counter(e["batch_date"] for e in embryos)
    spacing_seen = Counter(tuple(e["spacing_grid"]) if e["spacing_grid"] else None
                           for e in embryos)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "package_version": __version__,
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "data_root_redacted": config.redact_path(root),
        "n_experiments": len(experiments),
        "n_embryos": len(embryos),
        "by_stage": dict(by_stage),
        "by_batch": dict(by_batch),
        "spacing_grid_values_seen": {str(k): v for k, v in spacing_seen.items()},
        "physical_scale": {"xy_um_assumed": config.XY_UM, "z_over_xy": config.Z_OVER_XY,
                           "z_um_assumed": round(config.Z_UM, 4),
                           "note": "xy µm is a documented display constant from the MERFISH "
                                   "pipeline; all model features are dimensionless (÷ cell radius)."},
        "time_supervision": {
            "in_this_dataset": False,
            "note": "The mounted dataset is fixed MERFISH with no time-lapse / PNF / NEBD "
                    "annotations. True developmental tau comes from the external Scheffler 2021 "
                    "live-imaging trajectories; fixed stacks are the target imaging domain "
                    "(segmentation / adaptation / inference), never time supervision."},
        "experiments": experiments,
        "embryos": embryos,
    }
    if not include_paths:
        for e in manifest["embryos"]:
            e.pop("_paths", None)
    return manifest


def resolved_paths(manifest: dict, embryo_id: str) -> dict | None:
    for e in manifest["embryos"]:
        if e["embryo_id"] == embryo_id:
            return e.get("_paths")
    return None
