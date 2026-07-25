#!/usr/bin/env python3
"""
Generate orthogonal-view previews with segmentation overlays for the website
(task 19). For a CURATED set of embryos (spanning tau + failure cases) render:

  * DAPI MIP in XY, XZ, YZ (the 3D exploration planes);
  * the same three with cell-body / pronucleus-1 / pronucleus-2 / polar-body
    overlays (colours assigned from the GEOMETRIC audit, not names).

Previews are downsampled, embryo-cropped PNGs written to the gitignored derived
dir; a small committed index lists which previews exist so the site can load
them locally and degrade gracefully in production. No raw pixels enter git.

Usage: python3 scripts/build_pn3d_previews.py [--n 14]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import tifffile
from PIL import Image

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from scripts.pn3d import config, segment_audit as SA  # noqa: E402

OUT_DIR = os.path.join(config.DERIVED_DIR, "previews")
INDEX = os.path.join(config.DATA_DIR, "preview_index.json")
Z_STEP, XY_STEP = 2, 6
COLORS = {"cell": (150, 160, 175), "pn_near": (60, 130, 246),
          "pn_far": (220, 60, 60), "polar": (245, 175, 60)}


def _load(path, is_label):
    with tifffile.TiffFile(path) as t:
        mm = t.series[0].asarray(out="memmap")
        v = np.asarray(mm[::Z_STEP, ::XY_STEP, ::XY_STEP])
        del mm
    return v.astype(np.float32 if not is_label else np.int16)


def _norm(a):
    lo, hi = np.percentile(a, [1, 99.5])
    return np.clip((a - lo) / (hi - lo + 1e-6), 0, 1)


def _crop_bbox(mask, pad=6):
    idx = np.nonzero(mask)
    if len(idx[0]) == 0:
        return None
    b = [(int(ax.min()), int(ax.max()) + 1) for ax in idx]
    return b


def _project(vol, axis):
    return vol.max(axis=axis)


def _overlay(gray01, structures, axis, bbox):
    """RGB overlay of grayscale projection + coloured structure projections."""
    g = (gray01 * 255).astype(np.uint8)
    rgb = np.stack([g, g, g], -1).astype(np.float32)
    for mask, col in structures:
        m2 = mask.max(axis=axis) > 0
        for c in range(3):
            rgb[..., c] = np.where(m2, 0.45 * rgb[..., c] + 0.55 * col[c], rgb[..., c])
    z0, z1 = bbox
    return np.clip(rgb, 0, 255).astype(np.uint8)


def _save(arr, path, out_px=320):
    im = Image.fromarray(arr)
    h, w = arr.shape[:2]
    s = out_px / max(h, w)
    im = im.resize((max(1, int(w * s)), max(1, int(h * s))), Image.BILINEAR)
    im.save(path)


def render_embryo(eid, dna_path, label_path, audit):
    dna = _load(dna_path, False)
    L = _load(label_path, True)
    # crop to cell body bbox (+pad) so views are embryo-centric
    cell_label = audit["cell_body_label"]
    filled = (L == cell_label)
    for pl in (audit.get("pronucleus_labels") or []):
        filled = filled | (L == pl)
    bb = _crop_bbox(filled, pad=8)
    if bb:
        (z0, z1), (y0, y1), (x0, x1) = bb
        dna = dna[z0:z1, y0:y1, x0:x1]; L = L[z0:z1, y0:y1, x0:x1]
    structs = []
    pl = audit.get("pronucleus_labels") or []
    if audit["cell_body_label"] is not None:
        structs.append(((L == audit["cell_body_label"]), COLORS["cell"]))
    if len(pl) >= 1:
        structs.append(((L == pl[0]), COLORS["pn_near"]))
    if len(pl) >= 2:
        structs.append(((L == pl[1]), COLORS["pn_far"]))
    if audit.get("polar_body_label"):
        structs.append(((L == audit["polar_body_label"]), COLORS["polar"]))

    planes = {"xy": 0, "xz": 1, "yz": 2}
    made = {}
    safe = eid.replace("/", "__")
    os.makedirs(OUT_DIR, exist_ok=True)
    for name, axis in planes.items():
        gray = _norm(_project(dna, axis))
        p1 = os.path.join(OUT_DIR, f"{safe}_{name}.png")
        _save((gray * 255).astype(np.uint8), p1)
        ov = _overlay(gray, structs, axis, (0, dna.shape[axis]))
        p2 = os.path.join(OUT_DIR, f"{safe}_{name}_seg.png")
        _save(ov, p2)
        made[name] = {"raw": os.path.relpath(p1, HERE), "seg": os.path.relpath(p2, HERE)}
    return {"embryo_id": eid, "safe": safe, "planes": made}


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=14)
    a = ap.parse_args()
    config.ensure_dirs()

    man = json.load(open(os.path.join(config.DERIVED_DIR, "manifest_local.json")))
    paths = {e["embryo_id"]: e["_paths"] for e in man["embryos"]}
    inf = json.load(open(os.path.join(config.DATA_DIR, "inference.json")))["embryos"]
    geo = {e["embryo_id"]: e for e in json.load(
        open(os.path.join(config.DATA_DIR, "segmentation_geometry.json")))["embryos"]}

    resolved = [e for e in inf if e["stage"] == "zygote" and e["inferable"]]
    resolved.sort(key=lambda e: e["pseudotime"]["tau_mean"])
    # curated: spread across tau + a couple OOD + a couple unresolved (failure cases)
    pick = []
    if resolved:
        idx = np.linspace(0, len(resolved) - 1, min(a.n - 4, len(resolved))).astype(int)
        pick += [resolved[i]["embryo_id"] for i in sorted(set(idx))]
    ood = [e["embryo_id"] for e in inf if e["stage"] == "zygote"
           and e["ood_level"] == "out_of_domain" and e["inferable"]][:2]
    unres = [e["embryo_id"] for e in inf if e["stage"] == "zygote" and not e["inferable"]][:2]
    for x in ood + unres:
        if x not in pick:
            pick.append(x)

    index = []
    for k, eid in enumerate(pick, 1):
        if eid not in paths:
            continue
        au = geo.get(eid)
        if not au or au["cell_body_label"] is None:
            continue
        t0 = time.time()
        try:
            rec = render_embryo(eid, paths[eid]["dna"], paths[eid]["label_tifs"][0], au)
            index.append(rec)
            print(f"  [{k}/{len(pick)}] {eid:40s} rendered ({time.time()-t0:.1f}s)")
        except Exception as ex:                                    # noqa: BLE001
            print(f"  [{k}/{len(pick)}] {eid}: ERROR {ex}")
    json.dump({"schema_version": 1, "n": len(index),
               "note": "previews are gitignored derived pixels; absent in production",
               "colors": COLORS, "embryos": index}, open(INDEX, "w"), indent=1)
    print(f"\nwrote {os.path.relpath(INDEX, HERE)} ({len(index)} embryos) + PNGs to derived/previews")
    return 0


if __name__ == "__main__":
    sys.exit(main())
