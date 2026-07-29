#!/usr/bin/env python3
"""
Build data/equatorial_sperm.json — each labelled sperm projected into the equatorial
project's aligned-outline frame, so the "plot sperm locations" toggle can drop a mark
on every carrier zygote's silhouette.

The aligned outline lives in 2-D coordinates (along-axis a, chord1), exactly as
build_equatorial builds it (BZ.cross_section_outline(pos1, cell_com, chord2, a, chord1)
→ points [d·a, d·chord1]). Here we reproduce ONLY that geometry (cell-body COM, the
polar-body axis and chord1) and project the sperm the same way — no per-gene work, so it
reads just the label mask for the ~23 sperm zygotes.

Run from the deploy repo root:  python3 build_equatorial_sperm.py
"""
import glob
import gzip
import json
import os

import numpy as np

from build_zygote import XY_UM, unit, mask_and_transcripts, detect_polar_body
from build_planes_all import ATLAS, SRC
from embryo_naming import embryo_label

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SPERM = os.path.join(DATA, "zygote_sperm.json")
ZY_MAN = os.path.join(DATA, "zygote_manifest.json")
OUT = os.path.join(DATA, "equatorial_sperm.json")


def main():
    sperm = {e["id"]: e for e in json.load(open(SPERM))["embryos"] if e.get("sperm_plot")}
    zman = {m["id"]: m for m in json.load(open(ZY_MAN))["embryos"]}
    out = []
    for zid, se in sorted(sperm.items()):
        scene_p = os.path.join(ATLAS, zid, "scene.json.gz")
        lab = glob.glob(os.path.join(SRC, zid, "*_label.tif"))
        if not (os.path.isfile(scene_p) and lab):
            continue
        d = json.load(gzip.open(scene_p, "rt"))
        try:
            pos, labels, voxvol, seg_of, inside = mask_and_transcripts(lab[0], {}, [])
        except Exception as e:                          # noqa: BLE001
            print(f"  !! {zid}: {e}"); continue
        if 1 not in labels:
            continue
        pb_label, _ = detect_polar_body(d, labels, d.get("z_scale", 7.0), inside)
        if pb_label is None:
            continue
        cell_com = pos[labels == 1].mean(axis=0)
        pb_com = pos[labels == pb_label].mean(axis=0)
        a = unit(pb_com - cell_com)
        ref = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        chord1 = unit(np.cross(a, ref))                 # matches build_equatorial's chord1
        sp_um = np.asarray(se["sperm_plot"], float) * XY_UM
        dsp = sp_um - cell_com
        out.append({
            "id": zid,
            "label": (zman.get(zid, {}).get("label")) or embryo_label(zid, "zygote") or zid,
            "uv": [round(float(dsp @ a), 2), round(float(dsp @ chord1), 2)],   # (along-axis, chord1) µm
        })
        print(f"  {zid}  uv={out[-1]['uv']}")

    with open(OUT, "w") as fh:
        json.dump({"embryos": out}, fh, separators=(",", ":"))
    print(f"\nwrote {len(out)} sperm zygotes → {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)} bytes)")


if __name__ == "__main__":
    main()
