"""
Hybrid path: image -> segmentation -> symmetric geometry -> FROZEN clock
(brief item 7).

This is the deployable route for the fixed-image problem. It does NOT train a
new image model; it reuses the already-validated geometry-to-tau clock
(pnpt-3.0.0) unchanged. The only new work is turning a label volume into the
clock's symmetric geometric features, then reporting tau with the clock's own
uncertainty interval and QC status.

Single source of truth: tau prediction and QC come from
`build_pronuclei_pseudotime` (B.predict / B.qc_status). The geometry extraction
mirrors that module's audited convention (cell = fill_holes(cyto | pronucleiA |
pronucleiB); pronuclei = the two largest inside-cytoplasm segments).

Honest gating — the clock speaks microns, so the hybrid REFUSES to emit a tau
unless BOTH hold:
  * voxel spacing is VERIFIED (sidecar voxel_um_verified=true, or trustworthy tags);
  * the segmentation actually contains a cell + two pronuclei with verified
    class semantics.
Otherwise it returns a 'blocked' result naming exactly what is missing, plus the
scale-free geometry it *could* compute (in voxel units) so the plumbing is
visible without a fabricated physical number.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from . import config, tiff_audit

HERE = config.REPO_ROOT
sys.path.insert(0, HERE)
import build_pronuclei as BP                     # noqa: E402  geometry constants + detect_pronuclei
import build_pronuclei_pseudotime as B           # noqa: E402  frozen predict + qc (single source)

FEATURE_COLS = ["nearer_to_center_um", "farther_to_center_um",
                "distance_sum_um", "distance_difference_um"]
CORE = ["nearer_to_center_um", "farther_to_center_um"]
MODEL_P = os.path.join(HERE, "data", "pseudotime_calibration", "model.json")


# ───────────────────────────── segmentation -> geometry ─────────────────────────────
def _centroid(mask: np.ndarray, voxel_um) -> np.ndarray:
    """Centroid of a boolean mask in physical microns given per-axis voxel size."""
    iz, iy, ix = np.nonzero(mask)
    vx, vy, vz = voxel_um
    return np.array([ix.mean() * vx, iy.mean() * vy, iz.mean() * vz], float)


def segment_to_geometry(sub: np.ndarray, voxel_um, cyto_label: int = 1,
                        pron_labels=None) -> dict:
    """
    Symmetric geometric features from a (downsampled) integer label volume.

    voxel_um is the physical size of ONE voxel of `sub` along (x, y, z). For the
    project's standard downsampled frame that is (DS_XY*XY_UM, DS_XY*XY_UM,
    DS_Z*Z_UM); callers pass whatever matches `sub`.

    Returns {ok, features?, reason?, structure...}. Never invents a scale.
    """
    from scipy.ndimage import binary_fill_holes
    labs = [int(v) for v in np.unique(sub) if v >= 1]
    if cyto_label not in labs:
        return {"ok": False, "reason": f"cytoplasm label {cyto_label} absent (labels={labs})",
                "labels_present": labs}
    if pron_labels is None:
        det = BP.detect_pronuclei(sub)
        if not det:
            return {"ok": False, "reason": "fewer than two pronuclei inside the cytoplasm",
                    "labels_present": labs}
        pron_labels = list(det)
    la, lb = pron_labels
    seg1, pA, pB = (sub == cyto_label), (sub == la), (sub == lb)
    if pA.sum() < 4 or pB.sum() < 4:
        return {"ok": False, "reason": "a detected pronucleus has too few voxels for a centroid",
                "labels_present": labs}
    filled = binary_fill_holes(seg1 | pA | pB)
    c = _centroid(filled, voxel_um)
    dA = float(np.linalg.norm(_centroid(pA, voxel_um) - c))
    dB = float(np.linalg.norm(_centroid(pB, voxel_um) - c))
    near, far = (dA, dB) if dA <= dB else (dB, dA)
    feats = {"nearer_to_center_um": round(near, 4), "farther_to_center_um": round(far, 4),
             "distance_sum_um": round(near + far, 4),
             "distance_difference_um": round(abs(far - near), 4)}
    return {"ok": True, "features": feats, "pron_labels": [int(la), int(lb)],
            "labels_present": labs,
            "cell_voxels": int(filled.sum())}


# ───────────────────────────── geometry -> frozen clock ─────────────────────────────
_MODEL = None


def _model():
    global _MODEL
    if _MODEL is None:
        _MODEL = json.load(open(MODEL_P))
    return _MODEL


def geometry_to_tau(features: dict, training_csv: str | None = None) -> dict:
    """
    Apply the frozen clock (reused from build_pronuclei_pseudotime) to symmetric
    features. Returns tau, empirical interval, and QC — the clock's own outputs.
    """
    m = _model()
    raw = B.predict(m["spec"], features)
    tau = float(np.clip(raw, 0.0, 1.0))
    hw = m["halfwidth_95"]
    # QC vs the training feature ranges + Mahalanobis (same as the deployed apply step)
    stats = m["feature_stats"]
    tr_csv = training_csv or B.TRAIN_CSV
    mahal = None
    thr = {"caution": float("inf"), "extreme": float("inf")}
    if os.path.isfile(tr_csv):
        import csv as _csv
        tr = list(_csv.DictReader(open(tr_csv)))
        M = np.array([[float(r[c]) for c in CORE] for r in tr], float)
        mu, inv = M.mean(axis=0), np.linalg.pinv(np.cov(M, rowvar=False))
        d = np.array([features[c] for c in CORE], float) - mu
        mahal = float(np.sqrt(max(0.0, d @ inv @ d)))
        tm = np.sqrt(np.maximum(0.0, np.einsum("ij,jk,ik->i", M - mu, inv, M - mu)))
        thr = {"caution": float(np.quantile(tm, 0.975)), "extreme": float(np.quantile(tm, 0.999))}
    level, reasons = B.qc_status(features, stats, mahal, thr)
    return {
        "tau": round(tau, 5), "tau_raw": round(float(raw), 5),
        "lo95": round(max(0.0, tau - hw), 5), "hi95": round(min(1.0, tau + hw), 5),
        "halfwidth_95": hw, "qc": level, "qc_reasons": reasons,
        "mahalanobis": (round(mahal, 3) if mahal is not None else None),
        "model_version": m["model_version"], "interval_kind": "empirical (not guaranteed coverage)",
    }


def hybrid_predict(features: dict, spacing_status: str = "verified",
                   classes_verified: bool = True) -> dict:
    """
    The full hybrid decision. Emits a tau ONLY with verified µm spacing AND
    verified class semantics; otherwise returns a blocked result naming the gaps.
    """
    blockers = []
    if spacing_status != "verified":
        blockers.append(f"voxel spacing {spacing_status} — the clock needs microns")
    if not classes_verified:
        blockers.append("label class semantics unverified — cannot trust cell/pronucleus assignment")
    if blockers:
        return {"status": "blocked", "blockers": blockers, "features_units": "voxel (uncalibrated)",
                "note": "geometry plumbing ran; physical-unit tau withheld until blockers clear"}
    out = geometry_to_tau(features)
    out["status"] = "predicted"
    return out


# ───────────────────────────── end-to-end from a label stack ─────────────────────────────
def _load_sub(label_path: str) -> np.ndarray:
    """Downsampled label volume (Z,Y,X) via strided memmap — never the full stack."""
    import tifffile
    with tifffile.TiffFile(label_path) as t:
        mm = t.series[0].asarray(out="memmap")
        sub = np.asarray(mm[::BP.DS_Z, ::BP.DS_XY, ::BP.DS_XY])
        del mm
    return sub


def hybrid_from_label_stack(label_path: str, sidecar: str | None = None,
                            cyto_label: int = 1) -> dict:
    """
    Run the hybrid path on a segmentation stack, honestly gated by its sidecar.
    Extracts geometry in the downsampled voxel frame; converts to microns ONLY if
    the sidecar verifies voxel_um; predicts tau ONLY if spacing + classes verified.
    """
    audit = tiff_audit.audit(label_path, sidecar, declared_kind="segmentation_stack")
    res: dict = {"stack": audit["basename"], "spacing_status": audit.get("spacing_status"),
                 "sidecar_present": audit.get("sidecar_present")}
    sc = tiff_audit.load_sidecar(label_path, sidecar) or {}
    classes_verified = bool(sc.get("label_classes_verified"))
    res["classes_verified"] = classes_verified

    sub = _load_sub(label_path)
    # voxel size of the DOWNSAMPLED array: verified sidecar scale x downsample, else voxel units
    if audit.get("spacing_status") == "verified" and sc.get("voxel_um"):
        vx, vy, vz = sc["voxel_um"]
        vox = (vx * BP.DS_XY, vy * BP.DS_XY, vz * BP.DS_Z)
        units = "um"
    else:
        vox = (BP.DS_XY, BP.DS_XY, BP.DS_Z)      # scale-free voxel counts
        units = "voxel"
    geom = segment_to_geometry(sub, vox, cyto_label=cyto_label)
    res["geometry"] = geom
    res["geometry_units"] = units
    if not geom["ok"]:
        res["status"] = "blocked"
        res["blockers"] = [geom["reason"]]
        return res
    res.update(hybrid_predict(geom["features"],
                              spacing_status=audit.get("spacing_status", "missing"),
                              classes_verified=classes_verified))
    return res
