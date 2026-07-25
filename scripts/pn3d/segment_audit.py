"""
Structure-aware segmentation audit (task 16).

The segment labels are GENERIC (Segment_1..N). This module establishes each
segment's biological identity GEOMETRICALLY — never from a name — and enforces
the biological constraints of a pronuclear-stage zygote:

  * cell body  = the largest segment that encloses the interior structures;
  * pronucleus = a COMPACT segment INSIDE the cell body; exactly two are required;
  * polar body = a compact segment OUTSIDE the cell body (borders background);
  * the polar body must never be counted as a pronucleus.

Every audit returns a status (resolved / unresolved), per-structure confidence,
and explicit ambiguity flags. Samples that fail the constraints are marked, not
forced. The "inside vs outside the cell" test is the audited dilation-shell test
from the validated pipeline, generalized to full 4-structure classification.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from . import config

# downsample for the audit: robust for centroids/volumes, bounded memory
Z_STEP, XY_STEP = 2, 6


def load_label(path: str, z_step: int = Z_STEP, xy_step: int = XY_STEP) -> np.ndarray:
    import tifffile
    with tifffile.TiffFile(path) as t:
        mm = t.series[0].asarray(out="memmap")
        L = np.asarray(mm[::z_step, ::xy_step, ::xy_step])
        del mm
    return L.astype(np.int16)


def voxel_um(z_step: int = Z_STEP, xy_step: int = XY_STEP):
    return np.array([config.XY_UM * xy_step, config.XY_UM * xy_step,
                     config.Z_UM * z_step], float)          # (x, y, z) µm per downsampled voxel


def _centroid_um(mask, vox):
    iz, iy, ix = np.nonzero(mask)
    return np.array([ix.mean() * vox[0], iy.mean() * vox[1], iz.mean() * vox[2]], float)


def _compactness(mask):
    """Sphericity in [0,1]: 1 = perfect sphere. Robust roundness proxy."""
    v = float(mask.sum())
    if v < 8:
        return 0.0
    surf = float((ndi.binary_dilation(mask) & ~mask).sum()) + 1e-6
    # sphericity = pi^(1/3) (6V)^(2/3) / A
    return float(np.pi ** (1 / 3) * (6 * v) ** (2 / 3) / surf)


def _containment(mask, cell_mask):
    """
    Fraction of a segment's shell enclosed once the cell body and this segment are
    filled together. ~1 for a pronucleus (even a peripheral 'bay' one carved out of
    the cytoplasm); low for a polar body detached across the perivitelline gap.
    Robust where a raw dilation-border ratio fails on boundary-touching pronuclei.
    """
    filled = ndi.binary_fill_holes(cell_mask | mask)
    shell = ndi.binary_dilation(mask, iterations=2) & ~mask
    tot = int(shell.sum())
    if tot == 0:
        return 0.0
    return float((shell & filled).sum()) / tot


def _inside_score(mask, cell_mask, bg_mask):
    """Border ratio (diagnostic) alongside fill containment (decision)."""
    shell = ndi.binary_dilation(mask, iterations=2) & ~mask
    n_cell = int((shell & cell_mask).sum())
    n_bg = int((shell & bg_mask).sum())
    return (n_cell + 1) / (n_bg + 1), n_cell, n_bg


def audit(L: np.ndarray, vox=None) -> dict:
    """Classify segments of one downsampled label volume. Never uses names."""
    if vox is None:
        vox = voxel_um()
    labs = [int(v) for v in np.unique(L) if v >= 1]
    flags, out = [], {"labels_present": labs, "voxel_um": [round(float(x), 3) for x in vox]}
    if not labs:
        return {**out, "status": "unresolved", "reason": "no segments", "flags": ["empty"]}

    vols = {s: int((L == s).sum()) for s in labs}
    # cell body = largest segment (cytoplasm); verified to border background
    cell = max(vols, key=vols.get)
    cell_mask = L == cell
    bg = L == 0
    cell_shell_bg = int((ndi.binary_dilation(cell_mask) & ~cell_mask & bg).sum())
    if cell_shell_bg == 0:
        flags.append("cell_body_does_not_border_background")

    # classify the remaining segments as inside (pronucleus candidate) or outside (polar body)
    inside, outside = [], []
    detail = {}
    for s in labs:
        if s == cell:
            continue
        m = L == s
        score, nc, nb = _inside_score(m, cell_mask, bg)
        contain = _containment(m, cell_mask)
        comp = _compactness(m)
        rec = {"label": s, "volume_vox": vols[s], "inside_score": round(score, 2),
               "containment": round(contain, 3), "compactness": round(comp, 3)}
        detail[s] = rec
        # inside if the fill encloses it (bay-safe) OR the border ratio is clearly interior
        (inside if (contain >= 0.6 or score > 3.0) else outside).append(rec)

    # pronuclei = the two most compact / largest INSIDE segments
    inside.sort(key=lambda r: (r["compactness"] >= 0.35, r["volume_vox"]), reverse=True)
    pron = inside[:2]
    extra_inside = inside[2:]
    if extra_inside:
        flags.append(f"{len(extra_inside)} extra inside segment(s) beyond two pronuclei")

    # polar body = an OUTSIDE compact segment (external to the cell body)
    outside.sort(key=lambda r: r["volume_vox"], reverse=True)
    polar = outside[0] if outside else None
    if polar and polar["compactness"] < 0.2:
        flags.append("external segment is not compact (uncertain polar body)")

    # A pronuclear-stage zygote should show TWO pronuclei. One is still a real,
    # measurable configuration (a late/merged mass, or an incompletely annotated
    # pair), so it is measured and clearly labelled rather than refused — the
    # model returns an answer for every input, with the uncertainty and QC that
    # the input deserves. Zero DNA inside the cell body is genuinely unmeasurable.
    if len(pron) >= 2:
        status = "resolved"
    elif len(pron) == 1:
        status = "single_pronucleus"
        flags.append("only 1 pronucleus inside the cell body — geometry uses the single "
                     "annotated mass, so tau is an extrapolation (flagged, not hidden)")
    else:
        status = "unresolved"
        flags.append("no pronucleus inside the cell body — nothing to measure")
    # a would-be pronucleus that actually sits outside must never be substituted
    if len(inside) < 2 and outside:
        flags.append("a candidate lies OUTSIDE the cell body; not substituted for a pronucleus")

    # ---- geometry (whenever at least one pronucleus was found) ----
    geom = None
    if pron:
        masks = [L == p["label"] for p in pron]
        filled = ndi.binary_fill_holes(cell_mask | np.logical_or.reduce(masks))
        c = _centroid_um(filled, vox)
        cens = [_centroid_um(m, vox) for m in masks]
        dists = [float(np.linalg.norm(x - c)) for x in cens]
        vox_um3 = float(np.prod(vox))
        cell_vol = float(filled.sum()) * vox_um3
        cell_r = float((3 * cell_vol / (4 * np.pi)) ** (1 / 3))
        vols = [float(m.sum()) * vox_um3 for m in masks]
        # GENERALIZED feature — the RMS distance of pronuclear DNA from the cell
        # centre. For two pronuclei this is sqrt((d_near² + d_far²)/2), which is
        # equally computable from the live-imaging cohort, so the clock is trained
        # on the identical quantity. For one it is that mass's distance. Defined
        # for any pronucleus count, which is what lets every embryo get a tau.
        rms = float(np.sqrt(np.mean(np.square(dists))))
        geom = {
            "n_pronuclei": len(pron),
            "cell_center_um": [round(float(x), 2) for x in c],
            "cell_volume_um3": round(cell_vol, 1), "cell_radius_um": round(cell_r, 2),
            "pron_labels": [int(p["label"]) for p in pron],
            "pron_distances_um": [round(d, 3) for d in dists],
            "pron_volumes_um3": [round(v, 1) for v in vols],
            "rms_to_center_um": round(rms, 3),
            "rms_over_R": round(rms / cell_r, 4),          # ← the universal model input
            "pron_vol_frac": round(sum(vols) / (cell_vol + 1e-6), 5),
            "polar_body_present": polar is not None,
            "polar_body_external": polar is not None,      # by construction (outside score)
        }
        if len(pron) >= 2:                                  # two-pronucleus extras
            near, far = sorted(dists)[:2]
            inter = float(np.linalg.norm(cens[0] - cens[1]))
            vA, vB = vols[0], vols[1]
            geom.update({
                "nearer_to_center_um": round(near, 3), "farther_to_center_um": round(far, 3),
                "distance_sum_um": round(near + far, 3), "inter_pn_um": round(inter, 3),
                "near_over_R": round(near / cell_r, 4), "far_over_R": round(far / cell_r, 4),
                "sum_over_R": round((near + far) / cell_r, 4),
                "inter_over_R": round(inter / cell_r, 4),
                "vol_asymmetry": round(abs(vA - vB) / (vA + vB + 1e-6), 4),
            })

    # confidence: compact instances + a clear inside/outside margin + the expected
    # count. A single-pronucleus read is capped, so downstream uncertainty widens
    # for it rather than the sample being dropped.
    conf = 0.0
    if pron:
        comp_ok = float(np.mean([p["compactness"] for p in pron]))
        margin = min(p["inside_score"] for p in pron)
        conf = float(np.clip(0.4 * min(comp_ok / 0.5, 1) + 0.4 * min(margin / 3, 1)
                             + 0.2 * (len(extra_inside) == 0), 0, 1))
        if status == "single_pronucleus":
            conf = round(min(conf, 0.45), 3)
    out.update({
        "status": status, "cell_body_label": cell,
        "pronucleus_labels": [p["label"] for p in pron],
        "polar_body_label": (polar["label"] if polar else None),
        "n_inside": len(inside), "n_outside": len(outside),
        "confidence": round(conf, 3), "flags": flags,
        "segment_detail": detail, "geometry": geom,
    })
    return out


def audit_path(label_path: str, z_step: int = Z_STEP, xy_step: int = XY_STEP) -> dict:
    return audit(load_label(label_path, z_step, xy_step), voxel_um(z_step, xy_step))
