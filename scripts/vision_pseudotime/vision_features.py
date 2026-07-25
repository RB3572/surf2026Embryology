"""
DAPI image feature extraction for the fixed-zygote vision model.

The DAPI (DNA) channel images the two pronuclei directly as DNA-dense masses.
As the pronuclei migrate together with developmental time, the DNA cloud goes
from two separated masses (spread out) to one compact mass. This module turns a
zygote's DAPI z-stack into a rotation/reflection-invariant feature vector that
captures that geometry — WITHOUT using the segmentation label file, so it is a
genuine image→geometry predictor (the segmentation is used only to make the
training targets).

Every feature is invariant to the embryo's absolute position and in-plane pose
(moments, eigen-spreads, multi-threshold blob separations, radial DNA profile,
intensity statistics), so the regressor cannot exploit pose or field position.

Loading is bounded: a strided (downsampled) read, never the full ~284 MB volume.
Deterministic: identical input + parameters -> identical features.
"""
from __future__ import annotations

import numpy as np
import tifffile
from scipy import ndimage as ndi

Z_STEP = 2
XY_STEP = 2
N_RINGS = 10
THRESH_PCTS = (95.0, 98.0, 99.0, 99.5)

FEATURE_NAMES = (
    ["spread_rms", "spread_axis0", "spread_axis1", "spread_axis2",
     "aniso_01", "aniso_02", "planarity"]
    + [f"blobsep_p{int(p)}" for p in THRESH_PCTS]
    + [f"blobsep_p{int(p)}_norm" for p in THRESH_PCTS]
    + ["blob_count_p98", "top2_massfrac_p98", "principal_bimodality"]
    + [f"dnavol_p{int(p)}" for p in THRESH_PCTS]
    + [f"radial_{k}" for k in range(N_RINGS)]
    + ["intensity_cv", "intensity_p90_over_p50", "dna_extent_rms"]
)


def load_dapi(path: str, z_step: int = Z_STEP, xy_step: int = XY_STEP) -> np.ndarray:
    with tifffile.TiffFile(path) as t:
        mm = t.series[0].asarray(out="memmap")
        v = np.asarray(mm[::z_step, ::xy_step, ::xy_step]).astype(np.float32)
        del mm
    return v


def _weighted_cov_eig(v: np.ndarray, floor_pct: float = 50.0):
    """Eigenvalues of the intensity-weighted position covariance (spread of DNA)."""
    nz = v[v > 0]
    floor = np.percentile(nz, floor_pct) if nz.size else 0.0
    w = np.clip(v - floor, 0.0, None)
    W = float(w.sum())
    if W <= 0:
        return np.zeros(3), np.zeros(3)
    zz, yy, xx = np.mgrid[0:v.shape[0], 0:v.shape[1], 0:v.shape[2]]
    cz = float((zz * w).sum() / W); cy = float((yy * w).sum() / W); cx = float((xx * w).sum() / W)
    dz, dy, dx = zz - cz, yy - cy, xx - cx
    cov = np.array([
        [(w * dz * dz).sum(), (w * dz * dy).sum(), (w * dz * dx).sum()],
        [(w * dy * dz).sum(), (w * dy * dy).sum(), (w * dy * dx).sum()],
        [(w * dx * dz).sum(), (w * dx * dy).sum(), (w * dx * dx).sum()],
    ]) / W
    ev = np.sort(np.linalg.eigvalsh(cov))[::-1]
    return ev, np.array([cz, cy, cx])


def _blob_sep(v: np.ndarray, pct: float):
    """Separation (voxels) of the two largest DNA components above a percentile."""
    nz = v[v > 0]
    if nz.size == 0:
        return 0.0, 0, 0.0
    thr = np.percentile(nz, pct)
    mask = v > thr
    lab, n = ndi.label(mask)
    if n < 1:
        return 0.0, 0, 0.0
    sizes = np.asarray(ndi.sum(np.ones_like(lab, np.float32), lab, range(1, n + 1)))
    order = np.argsort(sizes)[::-1] + 1
    if n < 2:
        return 0.0, int(n), 1.0
    coms = ndi.center_of_mass(mask, lab, [int(order[0]), int(order[1])])
    (z1, y1, x1), (z2, y2, x2) = coms
    sep = float(np.sqrt((z1 - z2) ** 2 + (y1 - y2) ** 2 + (x1 - x2) ** 2))
    top2_frac = float((sizes[order[0] - 1] + sizes[order[1] - 1]) / max(sizes.sum(), 1))
    return sep, int(n), top2_frac


def _radial_profile(v: np.ndarray, center: np.ndarray, n_rings: int = N_RINGS):
    """Intensity-weighted mean DAPI in concentric shells about the DNA centroid."""
    zz, yy, xx = np.mgrid[0:v.shape[0], 0:v.shape[1], 0:v.shape[2]]
    r = np.sqrt((zz - center[0]) ** 2 + (yy - center[1]) ** 2 + (xx - center[2]) ** 2)
    rmax = float(r.max()) if r.max() > 0 else 1.0
    edges = np.linspace(0, rmax, n_rings + 1)
    prof = []
    for k in range(n_rings):
        m = (r >= edges[k]) & (r < edges[k + 1])
        prof.append(float(v[m].mean()) if m.any() else 0.0)
    prof = np.asarray(prof)
    return prof / (prof.max() + 1e-8)


def _principal_bimodality(v: np.ndarray, center: np.ndarray, axis: np.ndarray):
    """Sarle's bimodality coefficient of the DNA mass projected on the principal axis.
    Two separated pronuclei -> bimodal projection -> high coefficient."""
    nz = v > np.percentile(v[v > 0], 80) if (v > 0).any() else v > 0
    idx = np.nonzero(nz)
    if len(idx[0]) < 8:
        return 0.0
    pts = np.stack([idx[0] - center[0], idx[1] - center[1], idx[2] - center[2]], axis=1)
    w = v[nz]
    proj = pts @ axis
    m = np.average(proj, weights=w)
    var = np.average((proj - m) ** 2, weights=w)
    if var <= 1e-9:
        return 0.0
    sd = np.sqrt(var)
    skew = np.average(((proj - m) / sd) ** 3, weights=w)
    kurt = np.average(((proj - m) / sd) ** 4, weights=w)
    return float((skew ** 2 + 1.0) / max(kurt, 1e-6))


def extract(v: np.ndarray) -> dict:
    """Full feature dict for a downsampled DAPI volume. Deterministic."""
    ev, center = _weighted_cov_eig(v)
    spread_rms = float(np.sqrt(ev.sum()))
    ax = np.sqrt(np.clip(ev, 0, None))
    aniso01 = float(ax[0] / (ax[1] + 1e-6))
    aniso02 = float(ax[0] / (ax[2] + 1e-6))
    planarity = float((ax[1] - ax[2]) / (ax[0] + 1e-6))

    # principal axis for bimodality
    _, evec_center = ev, center
    # recompute eigenvectors for the projection axis
    axis = _principal_axis(v, center)

    seps, norms = [], []
    n98 = top2_98 = 0
    for p in THRESH_PCTS:
        s, n, t2 = _blob_sep(v, p)
        seps.append(s); norms.append(s / (spread_rms + 1e-6))
        if p == 98.0:
            n98, top2_98 = n, t2
    bim = _principal_bimodality(v, center, axis)

    dnavol = []
    nz = v[v > 0]
    for p in THRESH_PCTS:
        thr = np.percentile(nz, p) if nz.size else 0
        dnavol.append(float((v > thr).mean()))

    radial = _radial_profile(v, center)

    cv = float(nz.std() / (nz.mean() + 1e-6)) if nz.size else 0.0
    p90_p50 = float(np.percentile(nz, 90) / (np.percentile(nz, 50) + 1e-6)) if nz.size else 0.0

    vals = ([spread_rms, ax[0], ax[1], ax[2], aniso01, aniso02, planarity]
            + seps + norms + [float(n98), top2_98, bim] + dnavol
            + list(radial) + [cv, p90_p50, spread_rms])
    return dict(zip(FEATURE_NAMES, [float(x) for x in vals]))


def _principal_axis(v: np.ndarray, center: np.ndarray) -> np.ndarray:
    nz = v[v > 0]
    floor = np.percentile(nz, 50) if nz.size else 0.0
    w = np.clip(v - floor, 0.0, None)
    W = float(w.sum())
    if W <= 0:
        return np.array([1.0, 0.0, 0.0])
    zz, yy, xx = np.mgrid[0:v.shape[0], 0:v.shape[1], 0:v.shape[2]]
    dz, dy, dx = zz - center[0], yy - center[1], xx - center[2]
    cov = np.array([
        [(w * dz * dz).sum(), (w * dz * dy).sum(), (w * dz * dx).sum()],
        [(w * dy * dz).sum(), (w * dy * dy).sum(), (w * dy * dx).sum()],
        [(w * dx * dz).sum(), (w * dx * dy).sum(), (w * dx * dx).sum()],
    ]) / W
    evals, evecs = np.linalg.eigh(cov)
    return evecs[:, int(np.argmax(evals))]


def extract_from_path(path: str, z_step: int = Z_STEP, xy_step: int = XY_STEP) -> dict:
    return extract(load_dapi(path, z_step, xy_step))
