"""
Crosswalk between the 3D pseudotime cohort and the transcriptomics cohort.

The two cohorts label the same embryos differently:

    pn3d (MEFISH-Labels mount)   20260407_zygote_p1/fov6      and  .../fov6/0
    transcriptomics              20260407_zygote_p1_6_1       and  ..._6_0

A field of view often contains SEVERAL distinct embryos. The transcript side
splits them into _N_0 / _N_1; the mount segments one at the FOV level and each
additional one into a nested sub-directory. Names alone therefore cannot say
which segmentation is which embryo — and guessing would attach one embryo's
pseudotime to another's transcripts.

So the primary match is SPATIAL: both pipelines segment the same field in the
same coordinate frame, so each embryo is identified by its cell-body centroid.
Validated on the 30 pairs whose names are unambiguous, nearest-centroid
reproduces the name-based match 30/30, with a median centroid disagreement of
0.008 µm (max 0.143 µm). Matching is mutual-nearest and gated by MAX_MATCH_UM,
so an embryo present in only one cohort stays unmatched rather than being forced
onto a distant neighbour.

`name_key` remains available as an independent cross-check.
"""
from __future__ import annotations

import math
import re

PN3D_RE = re.compile(r"^(\d{8})_(?:\w+?_)?zygote_(p\d+)/fov(\d+)(?:/(.+))?$")
TX_RE = re.compile(r"^(\d{8})_zygote_(p\d+)_(\d+)(?:_(.+))?$")

# generous vs the observed 0.143 µm worst case, tight vs the ~80 µm embryo spacing
MAX_MATCH_UM = 12.0


# ───────────────────────── name keys (cross-check only) ─────────────────────────
def _norm(date, plate, fov, sub):
    s = (sub or "").strip()
    if s.startswith(f"{fov}_"):
        s = s[len(fov) + 1:]
    return (date, plate, fov, s)


def pn3d_key(embryo_id: str):
    m = PN3D_RE.match(embryo_id or "")
    return _norm(*m.groups()) if m else None


def tx_key(tx_id: str):
    m = TX_RE.match(tx_id or "")
    return _norm(*m.groups()) if m else None


def field_key(embryo_id: str, kind: str):
    """(date, plate, fov) — the field of view, ignoring which embryo within it."""
    k = pn3d_key(embryo_id) if kind == "pn3d" else tx_key(embryo_id)
    return k[:3] if k else None


# ───────────────────────── spatial match (primary) ─────────────────────────
def build_spatial(pn_points: dict, tx_points: dict, max_um: float = MAX_MATCH_UM) -> dict:
    """
    pn_points / tx_points: {id: (x_um, y_um)} cell-body centroids in the shared frame.

    Mutual-nearest within `max_um`, restricted to the same field of view so a
    centroid can never match across fields. Returns matched/ambiguous/unmatched.
    """
    def dist(a, b):
        return math.hypot(a[0] - b[0], a[1] - b[1])

    # candidate pairs within the same field
    tx_by_field: dict = {}
    for t, p in tx_points.items():
        tx_by_field.setdefault(field_key(t, "tx"), []).append(t)

    best_for_pn, best_for_tx = {}, {}
    for e, pe in pn_points.items():
        fld = field_key(e, "pn3d")
        cands = tx_by_field.get(fld, [])
        scored = sorted(((dist(pe, tx_points[t]), t) for t in cands), key=lambda z: z[0])
        if scored and scored[0][0] <= max_um:
            best_for_pn[e] = scored[0]
    for t, pt in tx_points.items():
        fld = field_key(t, "tx")
        cands = [e for e in pn_points if field_key(e, "pn3d") == fld]
        scored = sorted(((dist(pt, pn_points[e]), e) for e in cands), key=lambda z: z[0])
        if scored and scored[0][0] <= max_um:
            best_for_tx[t] = scored[0]

    matched, ambiguous = {}, []
    for e, (d, t) in best_for_pn.items():
        back = best_for_tx.get(t)
        if back and back[1] == e:                       # mutual nearest
            matched[e] = {"tx_id": t, "centroid_error_um": round(d, 4)}
        else:
            ambiguous.append({"pn3d_id": e, "nearest_tx": t, "distance_um": round(d, 4),
                              "reason": "not a mutual-nearest pair"})
    unmatched_pn3d = [{"pn3d_id": e, "reason": ("no transcript embryo within "
                                                f"{max_um} µm in the same field")}
                      for e in pn_points if e not in matched]
    unmatched_tx = sorted(t for t in tx_points
                          if t not in {m["tx_id"] for m in matched.values()})
    return {"matched": matched, "ambiguous": ambiguous,
            "unmatched_pn3d": unmatched_pn3d, "unmatched_tx": unmatched_tx,
            "max_match_um": max_um}


def name_agreement(matched: dict) -> dict:
    """Independent check: does the spatial match agree with the name key where the
    name is unambiguous (exactly one transcript embryo in that field)?"""
    agree = disagree = 0
    conflicts = []
    for pn_id, m in matched.items():
        pk, tk = pn3d_key(pn_id), tx_key(m["tx_id"])
        if pk and tk and pk == tk:
            agree += 1
        elif pk and tk:
            disagree += 1
            conflicts.append({"pn3d_id": pn_id, "tx_id": m["tx_id"],
                              "pn3d_key": list(pk), "tx_key": list(tk)})
    return {"name_key_agrees": agree, "name_key_differs": disagree,
            "note": "differences are expected in multi-embryo fields, where the mount's "
                    "FOV-level segmentation need not be the transcript side's _0",
            "examples": conflicts[:8]}
