#!/usr/bin/env python3
"""
Apply the frozen pronuclear pseudotime clock to the FIXED MERFISH zygotes.

SELF-CONTAINED BY DESIGN.  Geometry extraction needs the raw label TIFFs, which are (correctly) not
in this repository. So the extraction is CACHED: `calibration_data/fixed_cohort_geometry.csv` holds
the derived per-zygote geometry (distances, volumes, legacy gap) and nothing else — no transcripts,
no gene identities. The apply step therefore runs from the checked-in cache alone, in under a
second, with no dependency on the sibling embryo_viewer checkout, on an absolute iCloud path, or on
the raw microscopy.

  * default            : read the cached geometry, apply the frozen model  (self-contained)
  * --extract          : re-derive geometry from the raw TIFFs and refresh the cache
                         (needs build_pronuclei.py + the TranscriptomicsData tree; ~7 min)

GEOMETRY DEFINITIONS (used by --extract; already baked into the cache otherwise)
  * pronuclei   = the two largest segments INSIDE the cytoplasm, via the project's audited
                  dilation-shell test (which excludes the polar body — it borders background).
  * whole cell  = binary_fill_holes(seg1 | pronucleusA | pronucleusB).  The union is REQUIRED and
                  this was verified rather than assumed: the label volumes are mutually exclusive
                  (segment 1 has the pronuclei carved out) and fill_holes(seg1) ALONE fails to
                  enclose a pronucleus in 10 of 51 embryos — one touching the boundary is a bay,
                  not a hole.
  * features    = both pronuclear-centroid-to-cell-centre distances in physical 3-D um, sorted into
                  nearer/farther (identity-free: fixed zygotes have no reliable male/female call),
                  plus their sum and absolute difference. Exactly the training schema.
  * legacy      = the pre-existing minimum surface-to-surface pronuclear gap, carried through
                  UNCHANGED as a separate, clearly-labelled legacy field. The Scheffler workbook
                  contains no surface-gap measurement, so this calibration says NOTHING about it.

NO transcript, gene, probe or expression quantity is read anywhere in this file.

QC: pass / caution / out-of-domain from per-feature training range and a Mahalanobis check.
Out-of-domain embryos are still given a tau, but downstream analyses exclude them by default.

Usage:  python3 build_pronuclei_pseudotime.py [--extract]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_P = os.path.join(HERE, "data", "pseudotime_calibration", "model.json")
TRAIN_CSV = os.path.join(HERE, "calibration_data", "scheffler2021",
                         "scheffler_2021_control_zygote_trajectories.csv")
GEOM_CSV = os.path.join(HERE, "calibration_data", "fixed_cohort_geometry.csv")
OUT = os.path.join(HERE, "data", "pronuclei_pseudotime.json")

CORE = ["nearer_to_center_um", "farther_to_center_um"]
FEATURE_COLS = ["nearer_to_center_um", "farther_to_center_um",
                "distance_sum_um", "distance_difference_um"]


# ───────────────────────────── optional re-extraction from raw TIFFs ─────────────────────────────
def extract_geometry():
    """Re-derive the geometry cache from the raw label TIFFs. Only runs with --extract."""
    import glob
    from scipy.ndimage import binary_fill_holes
    sys.path.insert(0, HERE)
    try:
        import build_pronuclei as BP
    except ImportError as e:
        raise SystemExit(f"--extract needs build_pronuclei.py importable from {HERE}: {e}")
    man_p = os.path.join(HERE, "data", "pronuclei_manifest.json")
    if not os.path.isfile(man_p):
        raise SystemExit(f"--extract needs {man_p}")
    man = json.load(open(man_p))["embryos"]

    def centroid_um(mask):
        iz, iy, ix = np.nonzero(mask)
        return np.array([ix.mean() * BP.DS_XY * BP.XY_UM, iy.mean() * BP.DS_XY * BP.XY_UM,
                         iz.mean() * BP.DS_Z * BP.Z_UM], float)

    rows, n_bay = [], 0
    for i, e in enumerate(man):
        eid = e["id"]
        row = {"id": eid, "label": e["label"], "date_short": e.get("date_short", ""),
               "legacy_surface_gap_um": e["distance"], "extract_error": ""}
        try:
            lab = glob.glob(os.path.join(BP.SRC, eid, "*_label.tif"))
            if not lab:
                raise RuntimeError("no label TIFF found")
            sub = BP.load_sub(lab[0])
            la, lb = e["pron_labels"]
            seg1, pA, pB = (sub == BP.CYTO), (sub == la), (sub == lb)
            if not seg1.any():
                raise RuntimeError("segment 1 (cytoplasm) absent")
            if pA.sum() < 4 or pB.sum() < 4:
                raise RuntimeError("a detected pronucleus has too few voxels for a centroid")
            filled = binary_fill_holes(seg1 | pA | pB)
            f1 = binary_fill_holes(seg1)
            if min((f1 & pA).sum() / max(pA.sum(), 1), (f1 & pB).sum() / max(pB.sum(), 1)) < 0.5:
                n_bay += 1
            c = centroid_um(filled)
            dA = float(np.linalg.norm(centroid_um(pA) - c))
            dB = float(np.linalg.norm(centroid_um(pB) - c))
            near, far = (dA, dB) if dA <= dB else (dB, dA)
            vox = (BP.DS_XY * BP.XY_UM) ** 2 * BP.DS_Z * BP.Z_UM
            row.update(nearer_to_center_um=round(near, 4), farther_to_center_um=round(far, 4),
                       distance_sum_um=round(near + far, 4),
                       distance_difference_um=round(abs(far - near), 4),
                       cell_volume_um3=round(float(filled.sum()) * vox, 1),
                       cell_radius_equiv_um=round(float((3 * filled.sum() * vox / (4 * np.pi))
                                                        ** (1 / 3)), 2),
                       pron_volume_a_um3=round(float(pA.sum()) * vox, 1),
                       pron_volume_b_um3=round(float(pB.sum()) * vox, 1))
            print(f"  [{i+1}/{len(man)}] {e['label']:<12} near={near:6.2f} far={far:6.2f}")
        except Exception as ex:                                        # noqa: BLE001
            row["extract_error"] = str(ex)
            print(f"  -- {eid}: {ex}")
        rows.append(row)
    cols = ["id", "label", "date_short", "legacy_surface_gap_um", *FEATURE_COLS,
            "cell_volume_um3", "cell_radius_equiv_um", "pron_volume_a_um3", "pron_volume_b_um3",
            "extract_error"]
    with open(GEOM_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    print(f"  refreshed {os.path.relpath(GEOM_CSV, HERE)} ({len(rows)} zygotes; "
          f"{n_bay} would be missed by fill_holes(seg1) alone)")


# ───────────────────────────── frozen-model application ─────────────────────────────
def apply_isotonic(spec, f):
    x = spec["sign"] * f[spec["features"][0]]
    return float(np.interp(x, spec["knots_x"], spec["knots_y"],
                           left=spec["knots_y"][0], right=spec["knots_y"][-1]))


def apply_monospline(spec, f):
    """Evaluate the frozen I-spline: cumulative-reversed B-spline basis with non-negative coefs."""
    from scipy.interpolate import BSpline
    x = np.array([spec["sign"] * f[spec["features"][0]]], float)
    sp = spec["spline"]
    t, k = np.asarray(sp["knots"], float), int(sp["degree"])
    n_basis = len(t) - k - 1
    B = np.empty((1, n_basis))
    for j in range(n_basis):
        c = np.zeros(n_basis); c[j] = 1.0
        B[:, j] = BSpline(t, c, k, extrapolate=True)(x)
    B = np.nan_to_num(B)
    I = np.cumsum(B[:, ::-1], axis=1)[:, ::-1]
    return float(spec["intercept"] + I @ np.asarray(spec["coef"], float))


def apply_linear(spec, f):
    return float(spec["a"] + spec["b"] * f[spec["features"][0]])


def apply_ridge(spec, f):
    p = spec["pipeline"]
    x = np.array([f[k] for k in spec["features"]], float)
    z = (x - np.array(p["scaler_mean"])) / np.array(p["scaler_scale"])
    return float(p["intercept"] + float(np.dot(np.array(p["coef"]), z)))


def predict(spec, f):
    kind = spec.get("kind", "")
    form = spec.get("form") or (spec.get("pipeline") or {}).get("form", "")
    if kind == "IsotonicOnFeature" or "isotonic" in form:
        return apply_isotonic(spec, f)
    if kind == "MonotoneSplineOnFeature" or "Ispline" in form:
        return apply_monospline(spec, f)
    if kind == "LinearOnFeature":
        return apply_linear(spec, f)
    if "coef" in (spec.get("pipeline") or {}):
        return apply_ridge(spec, f)
    raise RuntimeError(f"cannot apply a frozen model of kind {kind!r} — tree ensembles are not "
                       "analytically serialisable; re-select or export ONNX")


def qc_status(f, stats, mahal, thr):
    reasons, level = [], "pass"
    for k in CORE:
        v, s = f[k], stats[k]
        if v < s["min"] or v > s["max"]:
            level = "out-of-domain"
            reasons.append(f"{k}={v:.1f}um outside training range [{s['min']:.1f}, {s['max']:.1f}]")
        elif v < s["p01"] or v > s["p99"]:
            if level == "pass":
                level = "caution"
            reasons.append(f"{k}={v:.1f}um outside training p01-p99 [{s['p01']:.1f}, {s['p99']:.1f}]")
    if mahal is not None:
        if mahal > thr["extreme"]:
            level = "out-of-domain"
            reasons.append(f"Mahalanobis {mahal:.2f} > {thr['extreme']:.2f} (training p99.9)")
        elif mahal > thr["caution"]:
            if level == "pass":
                level = "caution"
            reasons.append(f"Mahalanobis {mahal:.2f} > {thr['caution']:.2f} (training p97.5)")
    return level, reasons


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", action="store_true",
                    help="re-derive the geometry cache from raw TIFFs (needs the raw data)")
    a = ap.parse_args()
    if a.extract:
        extract_geometry()
    for p, what in ((MODEL_P, "run scripts/train_pronuclear_pseudotime.py first"),
                    (GEOM_CSV, "run with --extract to build it from the raw TIFFs"),
                    (TRAIN_CSV, "the checked-in Scheffler derived CSV is missing")):
        if not os.path.isfile(p):
            raise SystemExit(f"missing {os.path.relpath(p, HERE)} — {what}")

    model = json.load(open(MODEL_P))
    spec, stats_, hw = model["spec"], model["feature_stats"], model["halfwidth_95"]
    print(f"applying {model['label']} ({model['model_version']}) to the fixed MERFISH zygotes")

    # training cloud for the multivariate domain check (geometry only)
    tr = list(csv.DictReader(open(TRAIN_CSV)))
    M = np.array([[float(r[c]) for c in CORE] for r in tr], float)
    mu, inv = M.mean(axis=0), np.linalg.pinv(np.cov(M, rowvar=False))

    def mahal(f):
        d = np.array([f[c] for c in CORE], float) - mu
        return float(np.sqrt(max(0.0, d @ inv @ d)))

    tm = np.array([mahal({CORE[0]: r[0], CORE[1]: r[1]}) for r in M])
    thr = {"caution": float(np.quantile(tm, 0.975)), "extreme": float(np.quantile(tm, 0.999))}

    rows, n_ok = [], 0
    for g in csv.DictReader(open(GEOM_CSV)):
        base = {"id": g["id"], "label": g["label"], "date_short": g.get("date_short", ""),
                "legacy_surface_gap_um": (float(g["legacy_surface_gap_um"])
                                          if g.get("legacy_surface_gap_um") else None),
                "model_version": model["model_version"], "data_version": model["data_version"],
                "feature_schema": model["features"]}
        if g.get("extract_error") or not g.get(FEATURE_COLS[0]):
            rows.append({**base, "tau": None, "lo95": None, "hi95": None, "qc": "out-of-domain",
                         "reason": g.get("extract_error") or "no geometry in the cache",
                         "features": None})
            continue
        f = {c: float(g[c]) for c in FEATURE_COLS}
        extra = {c: (float(g[c]) if g.get(c) else None)
                 for c in ("cell_volume_um3", "cell_radius_equiv_um",
                           "pron_volume_a_um3", "pron_volume_b_um3")}
        raw = predict(spec, f)
        tau = float(np.clip(raw, 0.0, 1.0))
        m = mahal(f)
        lvl, why = qc_status(f, stats_, m, thr)
        rows.append({**base, "tau": round(tau, 5), "tau_raw": round(raw, 5),
                     "out_of_range": bool(raw < 0 or raw > 1),
                     "lo95": round(max(0.0, tau - hw), 5), "hi95": round(min(1.0, tau + hw), 5),
                     "qc": lvl, "reason": "; ".join(why), "mahalanobis": round(m, 3),
                     "features": {**f, **extra}})
        n_ok += 1

    # cohort-level domain shift vs the live training distribution
    okr = [r for r in rows if r.get("features")]
    shift = {}
    for c in ("nearer_to_center_um", "farther_to_center_um", "distance_sum_um"):
        fv = np.array([r["features"][c] for r in okr], float)
        tv = np.array([float(r[c]) for r in tr], float)
        shift[c] = {"fixed": {"min": round(float(fv.min()), 2),
                              "median": round(float(np.median(fv)), 2),
                              "max": round(float(fv.max()), 2)},
                    "training": {"min": round(float(tv.min()), 2),
                                 "median": round(float(np.median(tv)), 2),
                                 "max": round(float(tv.max()), 2)},
                    "fixed_median_percentile_in_training":
                        round(float(np.mean(tv <= np.median(fv))) * 100, 1),
                    "n_above_training_max": int((fv > tv.max()).sum()),
                    "n_below_training_min": int((fv < tv.min()).sum())}
    taus = np.array([r["tau"] for r in okr], float)
    floor_n = int((taus <= 0.06).sum())
    shift["consequence"] = {
        "n_distinct_tau": int(len(set(np.round(taus, 4)))),
        "n_at_or_below_tau_0_06": floor_n,
        "note": ("The fixed cohort's cell-centred distances are systematically LARGER than the live "
                 "training distribution (median near the "
                 f"{shift['distance_sum_um']['fixed_median_percentile_in_training']:.0f}th training "
                 "percentile). Two indistinguishable explanations: (a) these fixed zygotes really "
                 "are sampled earlier in pronuclear migration, or (b) a systematic measurement "
                 "offset between Scheffler's live-imaging geometry and this project's fixed 3-D "
                 "segmentation. Separating them needs raw stacks processed through the identical "
                 f"pipeline. Consequence: {floor_n} of {len(okr)} fixed zygotes sit at or below "
                 "tau=0.06, so the early end is compressed and weakly ordered there.")}

    qc = {}
    for r in rows:
        qc[r["qc"]] = qc.get(r["qc"], 0) + 1
    payload = {"meta": {
        "model_version": model["model_version"], "data_version": model["data_version"],
        "model_label": model["label"], "features": model["features"],
        "model_class": model.get("model_class", "empirical geometry-to-time calibration"),
        "trained_at_utc": model["trained_at_utc"], "seed": model["seed"],
        "applied_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "halfwidth_95": hw, "nested_outer_test": model.get("nested_outer_test"),
        "n_total": len(rows), "n_predicted": n_ok, "qc_counts": qc,
        "geometry_source": os.path.relpath(GEOM_CSV, HERE),
        "cell_center_definition": ("centroid of binary_fill_holes(segment1 | pronucleusA | "
                                   "pronucleusB) in physical um; segment 1 has the pronuclei carved "
                                   "out, so the union is required"),
        "qc_thresholds": {"mahalanobis_caution": round(thr["caution"], 3),
                          "mahalanobis_extreme": round(thr["extreme"], 3), "core_features": CORE},
        "downstream_default": ("out-of-domain zygotes are EXCLUDED from downstream regressions by "
                               "default; they can be shown explicitly in the UI"),
        "domain_shift": shift, "transcript_data_used": False,
        "legacy_note": ("legacy_surface_gap_um is the pre-existing minimum surface-to-surface "
                        "pronuclear gap. The Scheffler workbook calibrates CELL-CENTRED DISTANCE "
                        "features and contains no surface-gap measurement, so this calibration "
                        "neither validates nor calibrates the legacy score.")},
        "embryos": rows}
    with open(OUT, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT) / 1024:.0f} KB)")
    print(f"  {n_ok}/{len(rows)} predicted · QC {qc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
