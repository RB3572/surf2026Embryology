#!/usr/bin/env python3
"""
Build data/sperm_pseudotime.json for the Sperm Location vs Pseudotime project.

For every zygote that has a LABELLED SPERM, this joins:
  * the sperm 3-D position and the two pronuclei + polar-body centroids
    (from pronuclei_assignments.json — the same geometry the Maternal/Paternal
    Pronucleus ID project uses),
  * the maternal / paternal identity of each pronucleus (that project's consensus),
  * the calibrated pronuclear pseudotime tau + its 95% interval
    (pronuclei_pseudotime.json, pnpt-3.0.0),

and precomputes the sperm-to-object distance (µm) for each selectable object:
polar body, maternal pronucleus, paternal pronucleus.

Coordinates are in the viewer's isotropic plot space; µm = plot-distance × 0.15
(XY_UM_PER_PIXEL, validated against the pronuclei project's reported distances to
within a few %). Zygotes whose maternal/paternal call is an even split have no
defined maternal/paternal object (those distances are null); the polar-body
distance is still defined.

Run from the deploy repo root:  python3 build_sperm_pseudotime.py
"""
from __future__ import annotations

import json
import math
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
ASSIGN = os.path.join(ROOT, "data", "pronuclei_assignments.json")
PT = os.path.join(ROOT, "data", "pronuclei_pseudotime.json")
OUT = os.path.join(ROOT, "data", "sperm_pseudotime.json")

UM_PER_UNIT = 0.15          # plot-space unit → µm (isotropic; XY_UM_PER_PIXEL)


def d3(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def main():
    assign = {r["id"]: r for r in json.load(open(ASSIGN))["embryos"]}
    pt = {e["id"]: e for e in json.load(open(PT))["embryos"]}

    out = []
    n_polar = n_split = 0
    for eid, r in assign.items():
        sp = r.get("sperm_plot")
        if sp is None:
            continue                        # only embryos with a labelled sperm
        e = pt.get(eid)
        if not e or e.get("tau") is None:
            continue                        # need a pseudotime
        c = r.get("consensus") or {}
        split = bool(c.get("split"))
        n_split += 1 if split else 0
        female = None if split else c.get("female")
        pron = [{"com": p["com_plot"], "vol": p.get("volume"), "label": p.get("label")}
                for p in r["pron"]]
        polar = (r.get("polar") or {}).get("com_plot")
        polar_label = (r.get("polar") or {}).get("label")
        if polar:
            n_polar += 1

        def um(obj):
            return round(d3(sp, obj) * UM_PER_UNIT, 2) if obj else None

        maternal = pron[female]["com"] if female is not None else None
        paternal = pron[1 - female]["com"] if female is not None else None
        out.append({
            "id": eid, "label": r["label"], "date": e.get("date_short", ""),
            "tau": round(float(e["tau"]), 4),
            "tau_lo": round(float(e.get("lo95", e["tau"])), 4),
            "tau_hi": round(float(e.get("hi95", e["tau"])), 4),
            "split": split, "female": female,
            "sperm": sp, "pron": pron, "polar": polar, "polar_label": polar_label,
            "dist_um": {"polar": um(polar), "maternal": um(maternal), "paternal": um(paternal)},
        })

    out.sort(key=lambda x: x["tau"])
    doc = {
        "unit_um_per_plot": UM_PER_UNIT,
        "objects": [
            {"key": "polar", "label": "Polar body", "color": "#8b6fc4"},
            {"key": "maternal", "label": "Maternal pronucleus ♀", "color": "#d6336c"},
            {"key": "paternal", "label": "Paternal pronucleus ♂", "color": "#2563eb"},
        ],
        "sperm_color": "#ff2d95",
        "model_version": json.load(open(PT))["meta"].get("model_version", ""),
        "n": len(out), "n_split": n_split, "n_with_polar": n_polar,
        "embryos": out,
    }
    with open(OUT, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(f"{len(out)} sperm-labelled zygotes  ·  {n_split} split (no M/P)  ·  {n_polar} with a polar body")
    print(f"tau range {out[0]['tau']}–{out[-1]['tau']}  ·  {os.path.getsize(OUT)/1024:.0f} KB")


if __name__ == "__main__":
    main()
