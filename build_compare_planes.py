#!/usr/bin/env python3
"""Build data/compare_planes.json.gz — the "Compare Division Planes" project.

Four candidate dividing planes, compared head-to-head under ONE consistent metric:

  1. polar      — best plane CONTAINING the polar-body axis (the 18-plane fan, 10° steps),
                  chosen per gene by concentration asymmetry.
  2. exhaustive — best plane from the every-plane search (each gene's stored density-best normal
                  over the ~20k-orientation grid, reused from build_planes_all).
  3. equatorial — ⟂ the polar-body axis and SHIFTED along it until the cytoplasm splits 50/50.
  4. sperm      — through the sperm, the cytoplasm COM and the polar-body centroid.

REBUILT ON THE HOUSE DEFINITIONS (embryo_stats.py). Two things changed, and both move numbers:

  · THE EQUATORIAL PLANE WAS IN THE WRONG PLACE. It used to pass through the centre of mass. A
    plane through the COM does not bisect an irregular cell — on 20251226_zygote_p1_1 it splits
    the cytoplasm 0.4937 / 0.5063, and the reference measures the same error across the cohort
    (0.490 to 0.509). It is now shifted along the axis until the split is exactly 50/50, which is
    what makes "equatorial" mean anything. Every equatorial number in this artifact changed.

  · IT NOW READS data/segments/, NOT THE ATLAS. Counts come from the per-molecule segment label
    and volume from the body label's own voxel volume, so cytoplasm-only is exact rather than
    inferred, and the body is found by volume rather than assumed to be label 1. This also makes
    the build runnable in this repo — the atlas it used to need is on another machine.

The metric is unchanged: a CONCENTRATION asymmetry on the cell body only,
|c_A − c_B| / (c_A + c_B) with c = count ÷ that side's cytoplasm volume.

Output: data/compare_planes.json.gz
"""
import gzip
import json
import math
import os
import sys

import numpy as np

import embryo_stats as ES

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
PLANES_ALL = os.path.join(DATA, "planes_all")
NORMALS = os.path.join(DATA, "planes_all_normals.json.gz")
SPERM = os.path.join(DATA, "zygote_sperm.json")
ASSIGN = os.path.join(DATA, "pronuclei_assignments.json")
ZY_MAN = os.path.join(DATA, "zygote_manifest.json")
OUT = os.path.join(DATA, "compare_planes.json.gz")

N_POLAR = 18              # the meridional fan, 10° steps
MIN_COUNT = 10            # cytoplasmic transcripts a gene needs in an embryo to be compared

PLANE_META = [
    {"key": "polar",      "label": "Polar-axis best",   "color": "#a855f7", "perGene": True},
    {"key": "exhaustive", "label": "Exhaustive best",   "color": "#f59e0b", "perGene": True},
    {"key": "equatorial", "label": "Equatorial",        "color": "#0ea5e9", "perGene": False},
    {"key": "sperm",      "label": "Sperm · COM · PB",  "color": "#ff2d95", "perGene": False},
]


def in_plane_basis(a):
    ref = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(a, ref); u /= np.linalg.norm(u)
    return u, np.cross(a, u)


def process(eid, sperm_plot, zs_source, polar_lbl, pa_iVol, normals):
    p = ES.scene_path("Zygote", eid)
    if not os.path.isfile(p):
        return None, "no segments scene"
    sc = ES.read_scene(p)
    zs = sc["z_scale"]
    bodies = ES.classify_body(sc)
    if len(bodies) != 1:
        return None, f"expected one body, found {len(bodies)}"
    body = bodies[0]
    Vtot = ES.seg_volumes(sc)[body]
    V, F = ES.mesh_of(sc, body)
    com = ES.vol_centroid(V, F)

    if polar_lbl is None or str(polar_lbl) not in sc["region_meshes"]:
        polar_lbl = ES.polar_label(sc)
    if polar_lbl is None:
        return None, "no polar body: no axis"
    pb = ES.mesh_of(sc, polar_lbl)[0].mean(0)
    axis = pb - com
    axis /= np.linalg.norm(axis)

    TX = ES.cytoplasm_positions(sc, body)
    if not TX:
        return None, "no cytoplasmic transcripts"

    # ---- the four plane families ----
    n_eq, o_eq = ES.equal_volume_plane(V, F, axis, com, exact_total=Vtot)
    u, w = in_plane_basis(axis)
    th = np.arange(N_POLAR) * (np.pi / N_POLAR)
    N_pol = np.stack([np.cos(t) * u + np.sin(t) * w for t in th])

    n_sd = None
    if sperm_plot is not None:
        frame = sperm_plot[2] / zs_source
        sp = np.array([sperm_plot[0] * ES.PX, sperm_plot[1] * ES.PX,
                       round(frame) * zs * ES.PX])
        n_sd = np.cross(sp - com, pb - com)
        n_sd = n_sd / np.linalg.norm(n_sd) if np.linalg.norm(n_sd) > 1e-9 else None

    volA_eq, volB_eq = ES.split_volumes(V, F, n_eq, o_eq, exact_total=Vtot)
    volA_pol = np.array([ES.split_volumes(V, F, N_pol[k], com, exact_total=Vtot)[0]
                         for k in range(N_POLAR)])
    volA_sd = ES.split_volumes(V, F, n_sd, com, exact_total=Vtot)[0] if n_sd is not None else None
    grid_cache = {}

    def asym(nA, nTot, vA):
        nB, vB = nTot - nA, Vtot - vA
        cA, cB = nA / max(vA, 1e-9), nB / max(vB, 1e-9)
        s = cA + cB
        return abs(cA - cB) / s if s > 0 else 0.0

    gene_out = {}
    for g, P in TX.items():
        n_in = len(P)
        if n_in < MIN_COUNT:
            continue
        rec = {"nc": n_in}
        # the equatorial plane is the one that is NOT through the COM
        rec["eq"] = int(((P - o_eq) @ n_eq > 0).sum())
        rec["sd"] = int(((P - com) @ n_sd > 0).sum()) if n_sd is not None else None
        aA = ES.side_counts(P, com, N_pol).astype(float)
        best = int(max(range(N_POLAR), key=lambda k: asym(aA[k], n_in, volA_pol[k])))
        rec["pb"] = {"n": [round(float(x), 5) for x in N_pol[best]],
                     "v": [round(float(volA_pol[best]), 1), round(float(Vtot - volA_pol[best]), 1)],
                     "c": [int(aA[best]), n_in - int(aA[best])],
                     "ang": round(best * 180.0 / N_POLAR, 1)}
        idx = pa_iVol.get(g)
        if idx is not None:
            if idx not in grid_cache:
                grid_cache[idx] = ES.split_volumes(V, F, normals[idx], com, exact_total=Vtot)[0]
            vA = grid_cache[idx]
            nA = int(ES.side_counts(P, com, normals[idx])[0])
            rec["ex"] = {"n": [round(float(x), 5) for x in normals[idx]],
                         "v": [round(float(vA), 1), round(float(Vtot - vA), 1)],
                         "c": [nA, n_in - nA]}
        else:
            rec["ex"] = None
        gene_out[g] = rec

    if not gene_out:
        return None, f"no gene reaches {MIN_COUNT} cytoplasmic transcripts"

    ex = sc.get("extents") or {}
    L_um = 0.5 * 1.3 * max((ex.get("x", [0, 1])[1] - ex.get("x", [0, 1])[0]) * ES.PX,
                           (ex.get("y", [0, 1])[1] - ex.get("y", [0, 1])[0]) * ES.PX, 1.0) \
        if ex else 46.0
    return {
        "id": eid, "z_scale": zs,
        # the page renders from THIS scene: the zygote scenes use a different z_scale (6.667 vs
        # 7.0), so mixing the two would misplace every plane by ~5% in z
        "scene": f"Zygote__{eid}.json.gz", "body": int(body),
        "com_um": [round(float(x), 4) for x in com],
        "com_plot": [round(float(x) / ES.PX, 2) for x in com],
        "L_um": round(float(L_um), 1), "Vtot": round(float(Vtot), 1),
        "axis_um": [round(float(x), 5) for x in axis],
        # the equatorial plane carries its own ORIGIN now — it no longer passes through the COM
        "eq": {"n": [round(float(x), 5) for x in n_eq],
               "o_um": [round(float(x), 4) for x in o_eq],
               "o_plot": [round(float(x) / ES.PX, 2) for x in o_eq],
               "shift_um": round(float((o_eq - com) @ n_eq), 4),
               "v": [round(float(volA_eq), 1), round(float(volB_eq), 1)]},
        "sd": ({"n": [round(float(x), 5) for x in n_sd],
                "v": [round(float(volA_sd), 1), round(float(Vtot - volA_sd), 1)]}
               if n_sd is not None else None),
        "genes": gene_out,
    }, None


def main():
    normals = np.asarray(json.load(gzip.open(NORMALS, "rt"))["normals"], float)
    sperm_of = {e["id"]: e["sperm_plot"] for e in json.load(open(SPERM))["embryos"]
                if e.get("sperm_plot")}
    man = {m["id"]: m for m in json.load(open(ZY_MAN))["embryos"]}
    polar_of = {k: v.get("polar_body_label") for k, v in man.items()}
    for e in json.load(open(ASSIGN))["embryos"]:
        if e.get("polar"):
            polar_of[e["id"]] = e["polar"]["label"]
    zsrc = {}
    import glob
    for f in glob.glob(os.path.join(DATA, "zygote", "*.json.gz")):
        zsrc[os.path.basename(f)[:-len(".json.gz")]] = json.load(gzip.open(f, "rt"))["z_scale"]

    out, skipped = [], []
    ids = ES.stage_ids("Zygote")
    for i, eid in enumerate(ids, start=1):
        pa = {}
        pp = os.path.join(PLANES_ALL, eid + ".json.gz")
        if os.path.isfile(pp):
            for r in json.load(gzip.open(pp, "rt"))["analysis"]["genes"]:
                pa[r["gene"]] = r["iVol"]
        try:
            rec, why = process(eid, sperm_of.get(eid), zsrc.get(eid, 6.667),
                               polar_of.get(eid), pa, normals)
        except Exception as exc:                                   # noqa: BLE001
            rec, why = None, str(exc)
        if rec is None:
            skipped.append({"id": eid, "reason": why})
            print(f"  -- [{i}/{len(ids)}] {eid}: {why}")
            continue
        m = man.get(eid, {})
        rec["label"] = m.get("label") or eid
        rec["date_short"] = m.get("date_short", "")
        rec["has_sperm"] = rec["sd"] is not None
        out.append(rec)
        f = rec["eq"]["v"][0] / rec["Vtot"]
        print(f"  [{i}/{len(ids)}] {eid:34s} {len(rec['genes']):3d} genes  "
              f"eq split {f:.6f}  shift {rec['eq']['shift_um']:+.2f} µm"
              f"{'  sperm' if rec['has_sperm'] else ''}")

    doc = {
        "planes": PLANE_META,
        "unit_um_per_plot": ES.PX,
        "minCount": MIN_COUNT,
        "n": len(out), "n_sperm": sum(1 for r in out if r["has_sperm"]),
        "skipped": skipped,
        "method": {
            "equatorial": "perpendicular to the polar-body axis and SHIFTED along it until the "
                          "cytoplasm splits exactly 50/50 — a plane through the centre of mass "
                          "does not bisect an irregular cell",
            "cytoplasm_only": "counts from the per-molecule segment label; volume is the body "
                              "label's own voxel volume, with the half-split from mesh clipping",
            "metric": "|cA-cB|/(cA+cB), c = count / that side's cytoplasm volume",
        },
        "embryos": out,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    fr = [r["eq"]["v"][0] / r["Vtot"] for r in out]
    sh = [abs(r["eq"]["shift_um"]) for r in out]
    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e6:.2f} MB)")
    print(f"  {len(out)} zygotes ({doc['n_sperm']} with a sperm plane), {len(skipped)} skipped")
    print(f"  equatorial split now {min(fr):.6f}–{max(fr):.6f} (was ~0.49–0.51 through the COM)")
    print(f"  the shift it needed: median {np.median(sh):.2f} µm, max {max(sh):.2f} µm")


if __name__ == "__main__":
    main()
