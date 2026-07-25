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

    status = "resolved"
    if len(pron) < 2:
        status = "unresolved"
        flags.append(f"only {len(pron)} pronucleus/pronuclei inside the cell body (need 2)")
    # a would-be pronucleus that actually sits outside must never be substituted
    if len(inside) < 2 and outside:
        flags.append("a candidate lies OUTSIDE the cell body; not substituted for a pronucleus")

    # ---- geometry (only when resolved) ----
    geom = None
    if status == "resolved":
        pA, pB = pron[0]["label"], pron[1]["label"]
        mA, mB = L == pA, L == pB
        filled = ndi.binary_fill_holes(cell_mask | mA | mB)
        c = _centroid_um(filled, vox)
        cA, cB = _centroid_um(mA, vox), _centroid_um(mB, vox)
        dA = float(np.linalg.norm(cA - c)); dB = float(np.linalg.norm(cB - c))
        near, far = (dA, dB) if dA <= dB else (dB, dA)
        vox_um3 = float(np.prod(vox))
        cell_vol = float(filled.sum()) * vox_um3
        cell_r = float((3 * cell_vol / (4 * np.pi)) ** (1 / 3))
        volA = float(mA.sum()) * vox_um3; volB = float(mB.sum()) * vox_um3
        inter = float(np.linalg.norm(cA - cB))
        geom = {
            "cell_center_um": [round(float(x), 2) for x in c],
            "cell_volume_um3": round(cell_vol, 1), "cell_radius_um": round(cell_r, 2),
            "nearer_to_center_um": round(near, 3), "farther_to_center_um": round(far, 3),
            "distance_sum_um": round(near + far, 3),
            "inter_pn_um": round(inter, 3),
            "pron_volume_near_um3": round(min(volA, volB) if dA <= dB else max(volA, volB), 1),
            "pron_volume_far_um3": round(max(volA, volB) if dA <= dB else min(volA, volB), 1),
            "pron_labels": [int(pA), int(pB)],
            # DIMENSIONLESS features (÷ cell radius) — scale/anisotropy-robust
            "near_over_R": round(near / cell_r, 4), "far_over_R": round(far / cell_r, 4),
            "sum_over_R": round((near + far) / cell_r, 4),
            "inter_over_R": round(inter / cell_r, 4),
            "vol_asymmetry": round(abs(volA - volB) / (volA + volB + 1e-6), 4),
            "pron_vol_frac": round((volA + volB) / (cell_vol + 1e-6), 5),
            "polar_body_present": polar is not None,
            "polar_body_external": polar is not None,     # by construction (outside score)
        }

    # confidence: clean two-inside + compact + a clear inside/outside margin
    conf = 0.0
    if status == "resolved":
        comp_ok = np.mean([p["compactness"] for p in pron])
        margin = min(p["inside_score"] for p in pron)
        conf = float(np.clip(0.4 * min(comp_ok / 0.5, 1) + 0.4 * min(margin / 3, 1)
                             + 0.2 * (len(extra_inside) == 0), 0, 1))
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
