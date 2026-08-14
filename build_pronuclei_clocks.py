#!/usr/bin/env python3
"""
Build data/pronuclei_clocks.json — two extra *geometric* time axes for the
"Transcripts vs Pronuclear Distance" project (pronuclei.html).

Alongside the established time axes (the min pronuclei surface gap, and the
frozen calibrated pseudotime τ), the page's Time-axis menu can order zygotes by
either of two migration distances, both of which GROW as the zygote ages:

  * mat_polar — maternal ♀ pronucleus  →  polar body.
      The maternal pronucleus forms at the animal pole beside the polar body and
      migrates inward toward the cell centre, so this distance increases with time.
  * sperm_pat — sperm entry point  →  paternal ♂ pronucleus.
      The paternal pronucleus forms at the (fixed, cortical) sperm-entry site and
      migrates inward, so this distance increases with time.

Geometry is taken from pronuclei_assignments.json (the same centroids the
Maternal/Paternal Pronucleus-ID project uses); maternal/paternal identity is that
project's consensus. µm = plot-distance × 0.15 (XY_UM_PER_PIXEL, isotropic — the
same scale build_sperm_pseudotime.py validated against the reported distances).

A distance is only defined when the inputs exist:
  * mat_polar needs a polar body AND a non-split maternal/paternal consensus;
  * sperm_pat needs a labelled sperm AND a non-split consensus.
Undefined distances are emitted as null (never 0); a split (no-consensus) zygote
has no maternal/paternal identity, so both are null for it.

Run from the deploy repo root:  python3 build_pronuclei_clocks.py
"""
from __future__ import annotations

import json
import math
import os

from embryo_naming import embryo_label

ROOT = os.path.dirname(os.path.abspath(__file__))
ASSIGN = os.path.join(ROOT, "data", "pronuclei_assignments.json")
PT = os.path.join(ROOT, "data", "pronuclei_pseudotime.json")     # only for the τ sanity print
OUT = os.path.join(ROOT, "data", "pronuclei_clocks.json")

UM_PER_UNIT = 0.15          # plot-space unit → µm (isotropic; XY_UM_PER_PIXEL)


def d3(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def um(a, b):
    return round(d3(a, b) * UM_PER_UNIT, 2) if (a and b) else None


def main():
    assign = json.load(open(ASSIGN))["embryos"]

    rows = []
    n_split = n_mp = n_sp = 0
    for r in assign:
        pron = r.get("pron") or []
        if len(pron) != 2:
            continue                                    # this project is strictly two-pronucleus
        c = r.get("consensus") or {}
        split = bool(c.get("split"))
        female = None if split else c.get("female")
        n_split += 1 if split else 0

        polar = (r.get("polar") or {}).get("com_plot")
        sperm = r.get("sperm_plot")
        maternal = pron[female]["com_plot"] if female is not None else None
        paternal = pron[1 - female]["com_plot"] if female is not None else None

        mat_polar = um(maternal, polar)                 # ♀ pronucleus → polar body
        sperm_pat = um(sperm, paternal)                 # sperm entry → ♂ pronucleus
        if mat_polar is None and sperm_pat is None:
            continue                                    # nothing to add for this zygote
        n_mp += 1 if mat_polar is not None else 0
        n_sp += 1 if sperm_pat is not None else 0
        rows.append({
            "id": r["id"], "label": embryo_label(r["id"]),
            "mat_polar": mat_polar, "sperm_pat": sperm_pat,
        })

    rows.sort(key=lambda x: x["id"])
    doc = {
        "unit_um_per_plot": UM_PER_UNIT,
        # every axis here is a migration distance that GROWS with developmental time
        "clocks": [
            {"key": "mat_polar", "label": "Maternal ♀ pronucleus → polar body",
             "unit": "µm", "larger_is_later": True, "n": n_mp},
            {"key": "sperm_pat", "label": "Sperm entry → paternal ♂ pronucleus",
             "unit": "µm", "larger_is_later": True, "n": n_sp},
        ],
        "n": len(rows), "n_split": n_split,
        "embryos": rows,
    }
    with open(OUT, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))

    print(f"{len(rows)} zygotes  ·  {n_mp} with ♀PN→polar  ·  {n_sp} with sperm→♂PN  ·  "
          f"{n_split} split (no M/P)  ·  {os.path.getsize(OUT)/1024:.1f} KB")

    # developmental-direction sanity: both distances should rise with the calibrated τ
    try:
        pt = {e["id"]: e for e in json.load(open(PT))["embryos"]}

        def corr(key):
            xy = [(pt[r["id"]]["tau"], r[key]) for r in rows
                  if r[key] is not None and pt.get(r["id"], {}).get("tau") is not None]
            n = len(xy)
            if n < 3:
                return None, n
            mx = sum(p[0] for p in xy) / n
            my = sum(p[1] for p in xy) / n
            sxx = sum((p[0] - mx) ** 2 for p in xy)
            syy = sum((p[1] - my) ** 2 for p in xy)
            sxy = sum((p[0] - mx) * (p[1] - my) for p in xy)
            return (sxy / math.sqrt(sxx * syy) if sxx and syy else 0.0), n

        for k in ("mat_polar", "sperm_pat"):
            r_, n_ = corr(k)
            print(f"   {k}: Pearson r vs τ = {r_:+.2f} (n={n_})" if r_ is not None
                  else f"   {k}: too few for a τ check")
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()
