#!/usr/bin/env python3
"""Build data/size.json.gz — how big is each embryo, and what shape is the space it occupies.

Two measurements, both from the segmentation meshes, both in microns, and both defined
identically to figure 10.1 (Index10) so the site and the figure can never disagree.

ZYGOTE — MEAN CORTEX RADIUS
--------------------------
Measured DIRECTIONALLY, not from volume. The cytoplasm surface is binned into NT x NP
directions about its own centroid, and the radius in each direction is the FURTHEST vertex
that way; the reported radius is the mean over occupied directions.

Furthest, not nearest, and this matters: the pronuclei are carved out of the cytoplasm label,
so that label's mesh also carries interior surfaces. The nearest vertex in a given direction is
frequently a pronuclear wall rather than the cortex, and a nearest-vertex rule would report the
pronucleus as if it were the cell boundary.

EVERY embryo gets this number, at every stage. For a 2-cell embryo it is measured over both
blastomeres taken as one body — a dumbbell has a mean radius just as a ball does, and having one
comparable size figure across all 151 embryos is worth more than restricting it to the shapes it
describes most naturally. Shape statistics (radial CV, min/max) are reported alongside so the
spread is visible, but they are never used to exclude an embryo.

2-CELL — THE BOX THAT HOLDS THE PAIR
------------------------------------
An oriented bounding box along the axis joining the two blastomere centres:

    LENGTH  extent along that axis                (the pair, end to end)
    HEIGHT  extent across it, i.e. twice the greatest perpendicular distance any blastomere
            vertex reaches from the axis — the minimum height a box must have to contain them

LENGTH is NOT the sum of the two blastomere diameters and should not be expected to match it:
the blastomeres press together and flatten, so the pair end-to-end is markedly shorter than two
free spheres would be. That is a real property of the embryo, not a measurement error.

THE POLAR BODY IS EXCLUDED THROUGHOUT, as are the nuclei. Both are separate segments and never
enter the cytoplasm or blastomere labels, so the exclusion is structural rather than a filter.

SEGMENTATION QUALITY
--------------------
Some cortex surfaces have holes or nicks, which drags the nearest radius in a direction well
below anything a cell can be. Those embryos are marked `noted: "irregular"` for information, but
they are NOT excluded and NOT filtered anywhere — every embryo is measured and every embryo is
plotted. The mark exists so an implausible entry near the top of a ranked list can be recognised
for what it is, not so it can be hidden.

Units: x,y pixels x 0.15 = um; mesh z is PLOT space, so divide by z_scale.

Output: data/size.json.gz
"""
import glob
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SRC = os.path.join(DATA, "segments")
OUT = os.path.join(DATA, "size.json.gz")

VERSION = "size-1.0.0"
XY_UM = 0.15
NT, NP = 18, 36                 # cortex direction bins (polar x azimuth) — matches figure 10.1
STAGES = {"Zygote": "zygote", "Early2Cell": "early2cell", "Late2Cell": "late2cell"}

# Flagging thresholds. Chosen from the observed distribution rather than picked a priori:
# radial CV is tight across the cohort (median 7.9%, max 17.8%) so it is a poor discriminator,
# whereas the min/mean radius ratio has a long low tail — the worst embryos have the cortex
# dipping to 3-8 um in some direction against a ~40 um mean, which is a hole in the label
# surface, not a shape any cell can take. MIN_RATIO catches those; MAX_CV is a backstop.
MIN_RATIO = 0.30                # nearest cortex radius, as a fraction of the mean
MAX_CV = 15.0                   # radial coefficient of variation, %


def unit(v):
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-12 else v


def verts_um(mesh, zs):
    v = np.asarray(mesh["verts"], float).reshape(-1, 3)
    return np.stack([v[:, 0] * XY_UM, v[:, 1] * XY_UM, v[:, 2] / zs], axis=1)


def cortex_radii(V):
    """Furthest-vertex radius per direction bin about the centroid (µm)."""
    c = V.mean(axis=0)
    d = V - c
    r = np.linalg.norm(d, axis=1)
    keep = r > 1e-9
    d, r = d[keep], r[keep]
    u = d / r[:, None]
    ti = np.clip((np.arccos(np.clip(u[:, 2], -1, 1)) / np.pi * NT).astype(int), 0, NT - 1)
    pi_ = np.clip(((np.arctan2(u[:, 1], u[:, 0]) + np.pi) / (2 * np.pi) * NP).astype(int), 0, NP - 1)
    R = np.zeros((NT, NP))
    np.maximum.at(R, (ti, pi_), r)
    return R[R > 0], c


def do_zygote(s, segs, meshes, zs):
    k = str(segs[0]["label"])                          # cytoplasm = largest segment
    if k not in meshes:
        return None, "no cytoplasm mesh"
    V = verts_um(meshes[k], zs)
    r, c = cortex_radii(V)
    if r.size < 50:
        return None, "too few cortex directions"
    mean_r = float(r.mean())
    lo, hi = float(r.min()), float(r.max())
    cv = float(r.std() / mean_r * 100)
    aspect = hi / max(lo, 1e-9)
    ratio = lo / max(mean_r, 1e-9)
    # Shape statistics are REPORTED, never used to drop an embryo: the mean radius is defined for
    # every embryo and every embryo is plotted. `noted` is informational only.
    noted = "irregular" if (ratio < MIN_RATIO or cv > MAX_CV) else ""
    return {
        "radius_um": round(mean_r, 3),
        "radius_min_um": round(lo, 3), "radius_max_um": round(hi, 3),
        "radial_cv_pct": round(cv, 2), "aspect": round(aspect, 2),
        "min_over_mean": round(ratio, 3),
        "n_dirs": int(r.size),
        "centre_um": [round(float(v), 3) for v in c],
        "noted": noted,
    }, None


def do_twocell(s, segs, meshes, zs):
    if len(segs) < 2:
        return None, "fewer than two segments"
    ka, kb = str(segs[0]["label"]), str(segs[1]["label"])
    if ka not in meshes or kb not in meshes:
        return None, "blastomere mesh missing"
    A, B = verts_um(meshes[ka], zs), verts_um(meshes[kb], zs)
    ca, cb = A.mean(axis=0), B.mean(axis=0)
    sep = float(np.linalg.norm(cb - ca))
    if sep < 5:
        return None, "blastomere centres too close"
    u = unit(cb - ca)

    P = np.vstack([A, B])
    t = (P - ca) @ u
    t0, t1 = float(t.min()), float(t.max())
    length = t1 - t0

    # a perpendicular frame, and the box's true half-extents in it
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(u @ tmp) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    v = unit(tmp - (tmp @ u) * u)
    w = np.cross(u, v)
    pv, pw = (P - ca) @ v, (P - ca) @ w
    # HEIGHT per the agreed definition: twice the greatest perpendicular reach of ANY blastomere
    # vertex, i.e. the minimum height a box must have to contain the pair.
    perp = P - ca
    perp = perp - np.outer(perp @ u, u)
    height = float(2 * np.linalg.norm(perp, axis=1).max())
    # the same quantity restricted to the LARGER blastomere, which is what figure 10.1 reports
    pA = A - ca
    pA = pA - np.outer(pA @ u, u)
    height_larger = float(2 * np.linalg.norm(pA, axis=1).max())

    centre = ca + ((t0 + t1) / 2) * u + ((pv.min() + pv.max()) / 2) * v + ((pw.min() + pw.max()) / 2) * w

    # The SAME mean cortex radius the zygotes get, measured over both blastomeres as one body,
    # so a single comparable size number exists for every embryo at every stage. A 2-cell embryo
    # is a dumbbell rather than a ball, so this is a mean radius of a non-spherical shape - it is
    # a size summary, not a claim that the embryo is round.
    r2, c2 = cortex_radii(P)
    return {
        "radius_um": round(float(r2.mean()), 3),
        "radius_min_um": round(float(r2.min()), 3), "radius_max_um": round(float(r2.max()), 3),
        "radial_cv_pct": round(float(r2.std() / r2.mean() * 100), 2),
        "n_dirs": int(r2.size),
        "centre_um": [round(float(x), 3) for x in c2],
        "length_um": round(length, 3),
        "height_um": round(height, 3),
        "height_larger_um": round(height_larger, 3),
        "sep_um": round(sep, 3),
        "com_a_um": [round(float(x), 3) for x in ca],
        "com_b_um": [round(float(x), 3) for x in cb],
        "box": {
            "centre_um": [round(float(x), 3) for x in centre],
            "u": [round(float(x), 5) for x in u],
            "v": [round(float(x), 5) for x in v],
            "w": [round(float(x), 5) for x in w],
            "half_u": round(length / 2, 3),
            "half_v": round(float(pv.max() - pv.min()) / 2, 3),
            "half_w": round(float(pw.max() - pw.min()) / 2, 3),
        },
        "noted": "",
    }, None


def main():
    files = sorted(glob.glob(os.path.join(SRC, "*.json.gz")))
    if not files:
        sys.exit("no scenes in data/segments — nothing to do")
    out, skipped = [], []
    for path in files:
        base = os.path.basename(path)
        stage = STAGES.get(base.split("__")[0])
        if not stage:
            continue
        s = json.load(gzip.open(path, "rt"))
        zs = s.get("z_scale") or 7.0
        meshes = s.get("region_meshes") or {}
        segs = sorted(s.get("segments") or [], key=lambda x: -x["volume"])
        if not segs:
            skipped.append((base, "no segments")); continue
        rec, why = (do_zygote if stage == "zygote" else do_twocell)(s, segs, meshes, zs)
        if rec is None:
            skipped.append((base, why)); continue
        rec.update({"id": s["id"], "stage": stage, "scene": base,
                    "z_scale": zs, "n_segments": len(segs)})
        out.append(rec)

    zyg = [e for e in out if e["stage"] == "zygote"]
    e2c = [e for e in out if e["stage"] == "early2cell"]
    l2c = [e for e in out if e["stage"] == "late2cell"]
    noted = [e["id"] for e in out if e.get("noted")]

    def med(v):
        return round(float(np.median(v)), 2) if v else None

    meta = {
        "version": VERSION,
        "n": len(out), "n_zygote": len(zyg), "n_e2c": len(e2c), "n_l2c": len(l2c),
        "n_noted": len(noted), "noted": noted,
        "skipped": len(skipped),
        "t_bins": NT, "psi_bins": NP,
        "min_ratio": MIN_RATIO, "max_cv_pct": MAX_CV,
        "medians": {
            "zygote_radius_um": med([e["radius_um"] for e in zyg]),
            "e2c_radius_um": med([e["radius_um"] for e in e2c]),
            "l2c_radius_um": med([e["radius_um"] for e in l2c]),
            "e2c_length_um": med([e["length_um"] for e in e2c]),
            "e2c_height_um": med([e["height_um"] for e in e2c]),
            "l2c_length_um": med([e["length_um"] for e in l2c]),
            "l2c_height_um": med([e["height_um"] for e in l2c]),
        },
        "sorts": {
            "zygote": [{"key": "radius_um", "label": "Mean cortex radius"}],
            "twocell": [{"key": "radius_um", "label": "Mean cortex radius"},
                        {"key": "length_um", "label": "Box length"},
                        {"key": "height_um", "label": "Box height"},
                        {"key": "sep_um", "label": "Blastomere separation"}],
        },
    }
    os.makedirs(DATA, exist_ok=True)
    with gzip.open(OUT, "wt") as fh:
        json.dump({"meta": meta, "embryos": out}, fh, separators=(",", ":"))

    print(f"  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e3:.0f} KB)")
    print(f"  {len(out)} embryos: {len(zyg)} zygote, {len(e2c)} early 2-cell, {len(l2c)} late 2-cell")
    print(f"  medians: zygote r {meta['medians']['zygote_radius_um']} µm · "
          f"e2c {meta['medians']['e2c_height_um']}×{meta['medians']['e2c_length_um']} · "
          f"l2c {meta['medians']['l2c_height_um']}×{meta['medians']['l2c_length_um']}")
    print(f"  {len(noted)} embryos noted as irregular (reported only — all {len(out)} are included)")
    for b, why in skipped:
        print(f"    skip {b}: {why}")
    # does the box-containing height differ from figure 10.1's larger-blastomere height?
    diff = [abs(e["height_um"] - e["height_larger_um"]) for e in out if "height_um" in e]
    if diff:
        print(f"  height (pair) vs height (larger blastomere): max difference {max(diff):.2f} µm")


if __name__ == "__main__":
    main()
