#!/usr/bin/env python3
"""
Build the "Equatorial Division Plane" project data.

Same per-gene analysis as the Division Plane projects, but the single plane per
zygote is the EQUATORIAL plane: the plane that passes through the cell centre of
mass and is PERPENDICULAR to the animal–vegetal (polar-body) axis. It divides the
zygote into an animal half (toward the polar body) and a vegetal half.

    a       = unit(polar-body COM − cell COM)      # animal–vegetal axis
    normal  = a                                     # the equatorial plane's normal IS the axis
    side A  = (p − cell_com)·a > 0                   # animal (polar-body) side

Because the plane's normal is the axis itself, the natural 2-D view is MERIDIONAL
(a slice containing the axis): x = distance from the equatorial plane (= along the
axis, + toward the polar body), y = an equatorial chord. The equatorial plane is
the vertical line x = 0. The cross-section outline is built in that (along-axis,
chord) frame.

Only zygotes with a detected polar body and segment 1 get a plane (≈ the whole
atlas, like the 18-plane project — no sperm needed).

Reuses build_zygote (imported as BZ) for the atlas readers, polar-body detection,
permutation p-values, the balloon variant, and the cross aggregate. Outputs:
  data/equatorial/<id>.json.gz,  data/equatorial_manifest.json,
  data/equatorial_cross.json.gz (+ _circ).
Run from the deploy repo root:  python3 build_equatorial.py
"""
import glob
import gzip
import json
import os

import numpy as np

import build_zygote as BZ
from embryo_naming import embryo_label

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = BZ.ATLAS
SRC = BZ.SRC
OUT_DIR = os.path.join(HERE, "data", "equatorial")
OUT_MANIFEST = os.path.join(HERE, "data", "equatorial_manifest.json")
OUT_CROSS = os.path.join(HERE, "data", "equatorial_cross.json.gz")
OUT_CROSS_CIRC = os.path.join(HERE, "data", "equatorial_cross_circ.json.gz")

XY_UM = BZ.XY_UM
Z_UM = BZ.Z_UM
N_NULL = BZ.N_NULL
RNG = np.random.default_rng(20260733)      # independent reproducible stream
RNG_C = np.random.default_rng(20260734)    # circularized variant


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def analyze_equatorial(pos1, cell_com, a, genes, tx1, voxvol, ex, zs, rng):
    """Division-plane analysis at the SINGLE equatorial plane (normal = a, through
    cell_com). Mirrors build_sperm_division.analyze_single but with n = a and a
    MERIDIONAL section frame. Returns the analysis dict or None."""
    n = a                                             # equatorial plane normal = the axis
    # meridional view basis: chord1 (equatorial chord) + chord2 (⊥ both). The 2-D outline
    # lives in (along-a, chord1); the slab that produces it is thin along chord2.
    ref = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    chord1 = unit(np.cross(a, ref))
    chord2 = unit(np.cross(a, chord1))

    # The plane is placed where it splits the CYTOPLASM INTO EQUAL VOLUMES, not through the
    # centroid. Every voxel of pos1 carries the same volume, so the equal-volume plane is exactly
    # the MEDIAN of the voxel projections onto the axis — no search, no tolerance. The centroid is
    # pulled toward whichever end is fatter, so the two differ whenever the zygote is not
    # symmetric about its equator, which is most of them.
    #
    # pos1 is segment 1 only, so any SEPARATELY SEGMENTED pronucleus or polar body is already
    # excluded from this volume. Where the segmentation did not resolve the pronuclei they sit
    # inside label 1 and cannot be excluded; `pn_resolved` records which is the case per embryo
    # so the page can say so rather than imply an exclusion that did not happen.
    t = pos1 @ n
    t_star = float(np.median(t))
    plane_pt = cell_com + (t_star - float(cell_com @ n)) * n
    com_offset = float(cell_com @ n) - t_star        # + = centroid sits animal of the equator

    proj = t - t_star
    VA = max(float((proj > 0).sum()) * voxvol, voxvol)      # animal (polar-body) half volume
    VB = max(float((proj <= 0).sum()) * voxvol, voxvol)     # vegetal half volume

    total_n = 0
    wp_vol = wp_cnt = 0.0
    diff_cnt_sum = diff_vol_sum = 0.0
    gene_rows = []
    for gi, g in enumerate(genes):
        P = tx1[g]
        nn = len(P)
        if nn == 0:
            continue
        gproj = (P - plane_pt) @ n
        a_o = int((gproj > 0).sum()); b_o = nn - a_o
        null_a = rng.binomial(nn, 0.5, N_NULL)
        null_a1 = int(rng.binomial(nn, 0.5)); nb1 = nn - null_a1
        p_vol, p_cnt = BZ.perm_pvals(a_o, nn, VA, VB, null_a)
        row = {
            "a": a_o, "b": b_o, "aV": a_o / VA, "bV": b_o / VB, "aC": a_o / nn, "bC": b_o / nn,
            "dCount": a_o - b_o, "dVol": a_o / VA - b_o / VB, "dNorm": (a_o - b_o) / nn,
            "na": null_a1, "nb": nb1, "naV": null_a1 / VA, "nbV": nb1 / VB,
            "ndCount": null_a1 - nb1, "ndVol": null_a1 / VA - nb1 / VB, "ndNorm": (null_a1 - nb1) / nn,
            "pVol": p_vol, "pCnt": p_cnt,
        }
        gene_rows.append({"gene": g, "idx": gi, "total": nn, "planes": [row],
                          "bestP_vol": 0, "bestP_cnt": 0, "bestDiff_vol": 0, "bestDiff_cnt": 0})
        wp_vol += nn * p_vol; wp_cnt += nn * p_cnt
        diff_cnt_sum += abs(a_o - b_o); diff_vol_sum += nn * abs(a_o / VA - b_o / VB)
        total_n += nn

    if total_n == 0:
        return None
    wp_vol /= total_n; wp_cnt /= total_n
    dm_cnt = diff_cnt_sum / total_n; dm_vol = diff_vol_sum / total_n

    def to_plot(p):
        return [p[0] / XY_UM, p[1] / XY_UM, p[2] * zs]

    L_um = 0.62 * 0.5 * max(ex["x"][1] - ex["x"][0], ex["y"][1] - ex["y"][0],
                            ex["z"][1] - ex["z"][0]) * XY_UM
    plane_geo = {
        "angle": 0.0,
        "a_plot": to_plot(chord1),          # plane in-plane basis → cross-section "axis" (y)
        "m_plot": to_plot(chord2),          # depth (slab) direction
        "normal_um": [round(float(x), 6) for x in n], "L": L_um,
        "volA": round(VA, 1), "volB": round(VB, 1),
        "comOffset": round(com_offset, 3),          # how far the centroid sits off the equator
        "plane_plot": [round(c, 2) for c in to_plot(plane_pt)],
        "plane_um": [round(float(c), 4) for c in plane_pt],
        "wpVol": round(wp_vol, 5), "wpCnt": round(wp_cnt, 5),
        "dmVol": round(dm_vol, 7), "dmCnt": round(dm_cnt, 5),
    }
    # MERIDIONAL outline: slab thin along chord2, projected onto (along-axis a, chord1).
    outline = BZ.cross_section_outline(pos1, plane_pt, chord2, a, chord1)
    return {
        "com_plot": [round(c, 2) for c in to_plot(cell_com)],
        "com_um": [round(float(c), 4) for c in cell_com],
        "axis_plot": [round(x, 5) for x in to_plot(a)],   # the true polar axis (for the 3-D axis line)
        "planes": [plane_geo], "best_planes": {"pVol": 0, "pCnt": 0, "diffVol": 0, "diffCnt": 0},
        "n_null": N_NULL,
        "cross_section": {"u_plot": to_plot(a), "v_plot": to_plot(chord1),
                          "outline": [[round(p[0], 2), round(p[1], 2)] for p in outline]},
        "genes": gene_rows,
    }


def process(eid):
    scene_p = os.path.join(ATLAS, eid, "scene.json.gz")
    lab = glob.glob(os.path.join(SRC, eid, "*_label.tif"))
    if not (os.path.isfile(scene_p) and lab):
        return None
    d = json.load(gzip.open(scene_p, "rt"))
    zs = d.get("z_scale", 7.0)

    tx = d.get("transcripts", {})
    genes = sorted(tx.keys(), key=lambda g: -d.get("gene_totals", {}).get(g, 0))
    pos, labels, voxvol, seg_of, inside_fractions = BZ.mask_and_transcripts(lab[0], tx, genes)
    if 1 not in labels:
        return None
    pb_label, pb_det = BZ.detect_polar_body(d, labels, zs, inside_fractions)
    if pb_label is None:
        print(f"  -- {eid}: {pb_det['reason']}")
        return None

    pos1 = pos[labels == 1]
    cell_com = pos1.mean(axis=0)                       # segment-1 COM (µm)
    pb_com = pos[labels == pb_label].mean(axis=0)
    a = unit(pb_com - cell_com)

    tx1 = {}
    for g in genes:
        t = tx[g]; in1 = seg_of[g] == 1
        tx1[g] = np.stack([np.asarray(t["x"], float)[in1] * XY_UM,
                           np.asarray(t["y"], float)[in1] * XY_UM,
                           np.asarray(t["gz"], float)[in1] * Z_UM], axis=1)

    analysis = analyze_equatorial(pos1, cell_com, a, genes, tx1, voxvol, d["extents"], zs, RNG)
    if analysis is None:
        return None

    def to_plot(p):
        return [round(p[0] / XY_UM, 2), round(p[1] / XY_UM, 2), round(p[2] * zs, 2)]

    analysis["pb_plot"] = to_plot(pb_com)
    analysis["polar_body_label"] = int(pb_label)
    analysis["polar_body_detection"] = pb_det

    scene = {
        "id": eid, "z_scale": zs, "extents": d["extents"],
        "region_meshes": d["region_meshes"], "region_defaults": d["region_defaults"],
        "mask_labels": d["mask_labels"], "genes": genes, "gene_totals": d.get("gene_totals", {}),
        "transcripts": {g: {"x": tx[g]["x"], "y": tx[g]["y"], "gz": tx[g]["gz"],
                            "s1": (seg_of[g] == 1).astype(np.uint8).tolist()} for g in tx},
        "analysis": analysis,
    }

    # circularized ("blow up the balloon") variant — segment 1 only
    fn, _C, _Ravg = BZ.balloon(pos1)
    pos1_c = fn(pos1)
    cell_com_c = pos1_c.mean(axis=0)
    a_c = unit(pb_com - cell_com_c)
    tx1_c = {g: (fn(P) if len(P) else P) for g, P in tx1.items()}
    analysis_c = analyze_equatorial(pos1_c, cell_com_c, a_c, genes, tx1_c, voxvol,
                                    d["extents"], zs, RNG_C)
    if analysis_c is not None:
        cs = analysis_c.get("cross_section")
        if cs and cs.get("outline"):
            cs["outline"] = BZ._circularize_outline(cs["outline"])
        analysis_c["pb_plot"] = to_plot(pb_com)
        m1 = d["region_meshes"].get("1")
        scene["circ"] = {"analysis": analysis_c,
                         "transcripts": BZ._circ_transcripts(tx, seg_of, fn),
                         "mesh1": BZ._circ_mesh(m1, fn, zs) if m1 else None,
                         "R_avg_um": round(float(_Ravg), 2)}
    else:
        scene["circ"] = None

    return scene


def agg_entry_eq(eid, label, A):
    """BZ.agg_entry plus per-side segment-1 volumes (vp/vt) for the density modes
    (single plane → one-entry arrays)."""
    e = BZ.agg_entry(eid, label, A)
    p0 = A["planes"][0]
    e["vp"] = [round(float(p0["volA"]), 1)]
    e["vt"] = round(float(p0["volA"] + p0["volB"]), 1)
    return e


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    ids = sorted(os.listdir(ATLAS))
    manifest = []; agg_entries = []; agg_circ = []
    MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    for i, eid in enumerate(ids):
        try:
            scene = process(eid)
        except Exception as e:              # noqa: BLE001
            print(f"  !! {eid}: {e}")
            continue
        if not scene:
            print(f"  -- skipped {eid}")
            continue
        out = os.path.join(OUT_DIR, eid + ".json.gz")
        with gzip.open(out, "wt") as fh:
            json.dump(scene, fh, separators=(",", ":"), default=BZ._json_default)
        date = eid[:8]
        ds = f"{MON[int(date[4:6]) - 1]} {int(date[6:8])}" if date.isdigit() else ""
        label = embryo_label(eid, "zygote")
        A = scene["analysis"]
        manifest.append({
            "id": eid, "label": label, "date_short": ds,
            "n_genes": len(scene["genes"]),
            "n_transcripts": sum(len(t["x"]) for t in scene["transcripts"].values()),
            "size_kb": round(os.path.getsize(out) / 1024),
            "polar_body_label": A["polar_body_label"],
            "wpVol": A["planes"][0]["wpVol"], "wpCnt": A["planes"][0]["wpCnt"],
        })
        agg_entries.append(agg_entry_eq(eid, label, A))
        if scene.get("circ"):
            agg_circ.append(agg_entry_eq(eid, label, scene["circ"]["analysis"]))
        print(f"  [{i+1}/{len(ids)}] {eid}  wpVol={A['planes'][0]['wpVol']:.3f} "
              f"wpCnt={A['planes'][0]['wpCnt']:.3f}  {manifest[-1]['size_kb']}KB")
    manifest.sort(key=lambda m: m["id"])
    agg_entries.sort(key=lambda e: e["id"])
    with open(OUT_MANIFEST, "w") as fh:
        json.dump({"embryos": manifest}, fh, indent=1)
    BZ.write_cross_aggregate(agg_entries, OUT_CROSS)
    if agg_circ:
        agg_circ.sort(key=lambda e: e["id"])
        BZ.write_cross_aggregate(agg_circ, OUT_CROSS_CIRC)
    tot = sum(m["size_kb"] for m in manifest)
    print(f"\nwrote {len(manifest)} equatorial-plane zygotes  ({tot/1024:.1f} MB)")


if __name__ == "__main__":
    main()
