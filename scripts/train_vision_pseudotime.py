#!/usr/bin/env python3
"""
Train the fixed-zygote DAPI vision model and evaluate it with leave-one-
embryo-out cross-validation.

Two deployable predictors, both from DAPI image features only (the segmentation
is used only to make the targets):

  A. distance_sum_um  (the geometry clock's feature)  ->  tau via the FROZEN
     validated clock (pnpt-3.0.0).  Also a direct image->tau model for comparison.
  B. pron_min_distance_um  (the pronuclei project's surface-to-surface gap).

Evaluation is honest: nested LOO. Each of the 51 zygotes is held out once; inside
each fold RidgeCV picks its regularization from the other 50 only. A final model
is then fit on all 51 for inference on new images. Deterministic (no RNG in
ridge; fixed feature order).

Artifacts:
  data/vision_pseudotime_model.json      (A: coeffs, LOO metrics, per-zygote preds)
  data/vision_pronuclei_distance.json    (B: per-zygote predicted gap + transcripts)

Usage:  python3 scripts/train_vision_pseudotime.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import build_pronuclei_pseudotime as CLOCK  # noqa: E402  frozen geometry->tau
from scripts.vision_pseudotime import evaluate as V  # noqa: E402
from scripts.vision_pseudotime import vision_features as VF  # noqa: E402

DATA = os.path.join(HERE, "calibration_data", "vision_pseudotime", "vision_features.csv")
MODEL_A = os.path.join(HERE, "data", "vision_pseudotime_model.json")
MODEL_B = os.path.join(HERE, "data", "vision_pronuclei_distance.json")
MODEL_VERSION = "vpt-vision-1.0.0"
SEED = 20260725


def load():
    rows = list(csv.DictReader(open(DATA)))
    X = np.array([[float(r[f]) for f in VF.FEATURE_NAMES] for r in rows], float)
    return rows, X


def loo_predict(X, y):
    """Nested leave-one-out: RidgeCV picks alpha inside each training fold."""
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    n = len(y)
    pred = np.zeros(n)
    alphas = np.logspace(-2, 4, 25)
    for i in range(n):
        tr = np.arange(n) != i
        model = make_pipeline(StandardScaler(), RidgeCV(alphas=alphas))
        model.fit(X[tr], y[tr])
        pred[i] = float(model.predict(X[i:i + 1])[0])
    return pred


def final_model(X, y):
    """Model refit on ALL rows, returned as plain coefficients for inference."""
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler().fit(X)
    rc = RidgeCV(alphas=np.logspace(-2, 4, 25)).fit(sc.transform(X), y)
    return {"feature_names": list(VF.FEATURE_NAMES),
            "scaler_mean": sc.mean_.tolist(), "scaler_scale": sc.scale_.tolist(),
            "coef": rc.coef_.tolist(), "intercept": float(rc.intercept_),
            "alpha": float(rc.alpha_)}


def metrics(y_true, y_pred, groups):
    return {
        "mae": float(np.mean(np.abs(y_true - y_pred))),
        "rmse": float(np.sqrt(np.mean((y_true - y_pred) ** 2))),
        "spearman": V.spearman(y_true, y_pred),
        "pearson": float(np.corrcoef(y_true, y_pred)[0, 1]),
        "pair_order_accuracy": V.pair_order_accuracy(y_true, y_pred)["strict_accuracy"],
        "n": int(len(y_true)),
    }


def main() -> int:
    rows, X = load()
    ids = [r["id"] for r in rows]
    labels = [r["label"] for r in rows]
    qc = [r.get("qc", "") for r in rows]
    dsum = np.array([float(r["distance_sum_um"]) for r in rows])
    tau = np.array([float(r["tau"]) for r in rows])
    mind = np.array([float(r["pron_min_distance_um"]) for r in rows])
    total = np.array([int(r["transcript_total"]) for r in rows])
    n = len(rows)
    print(f"loaded {n} zygotes, {X.shape[1]} DAPI features")

    # ---- A: image -> distance_sum -> frozen clock -> tau ----
    dsum_hat = loo_predict(X, dsum)
    spec = json.load(open(os.path.join(HERE, "data", "pseudotime_calibration", "model.json")))["spec"]
    # tau via the frozen clock, using predicted distance_sum (near/far unknown from image;
    # the clock only needs distance_sum for the isotonic-on-sum model)
    tau_via_clock = np.array([float(np.clip(CLOCK.predict(spec, {"distance_sum_um": float(d)}), 0, 1))
                              for d in dsum_hat])
    # B-style direct image -> tau (comparison)
    tau_direct = np.clip(loo_predict(X, tau), 0, 1)

    m_dsum = metrics(dsum, dsum_hat, ids)
    m_tau_clock = metrics(tau, tau_via_clock, ids)
    m_tau_direct = metrics(tau, tau_direct, ids)

    # empirical LOO interval half-width for tau (via clock path)
    resid = np.abs(tau - tau_via_clock)
    hw = float(np.quantile(resid, 0.95))
    lo = np.clip(tau_via_clock - hw, 0, 1); hi = np.clip(tau_via_clock + hw, 0, 1)
    cov = V.interval_coverage(tau, lo, hi)

    # baselines
    base_const = float(np.mean(np.abs(tau - np.mean(tau))))
    spread = X[:, VF.FEATURE_NAMES.index("spread_rms")]
    base_spread = metrics(tau, np.clip(np.polyval(np.polyfit(spread, tau, 1), spread), 0, 1), ids)["mae"]

    perzyg = [{"id": ids[i], "label": labels[i], "qc": qc[i],
               "true_tau": round(float(tau[i]), 4), "pred_tau": round(float(tau_via_clock[i]), 4),
               "pred_tau_direct": round(float(tau_direct[i]), 4),
               "true_distance_sum_um": round(float(dsum[i]), 3),
               "pred_distance_sum_um": round(float(dsum_hat[i]), 3),
               "lo95": round(float(lo[i]), 4), "hi95": round(float(hi[i]), 4)} for i in range(n)]

    payload_a = {
        "meta": {"model_version": MODEL_VERSION, "seed": SEED,
                 "trained_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
                 "n_zygotes": n, "n_features": X.shape[1],
                 "input": "DAPI z-stack image features only (no segmentation at inference)",
                 "target_source": "distance_sum & tau from the validated segmentation clock (labels only)",
                 "evaluation": "nested leave-one-embryo-out CV (each zygote held out once)",
                 "clock_reused": "pnpt-3.0.0",
                 "tau_path": "image -> distance_sum (vision) -> frozen clock -> tau",
                 "interval_halfwidth_95_loo": round(hw, 4),
                 "interval_coverage_loo": round(float(cov), 4)},
        "metrics": {"distance_sum_um": m_dsum, "tau_via_clock": m_tau_clock,
                    "tau_direct": m_tau_direct},
        "baselines": {"tau_constant_mae": round(base_const, 4),
                      "tau_single_feature_spread_mae": round(base_spread, 4)},
        "final_model_distance_sum": final_model(X, dsum),
        "final_model_tau_direct": final_model(X, tau),
        "per_zygote": perzyg,
    }
    json.dump(payload_a, open(MODEL_A, "w"), indent=1)

    # ---- B: image -> pronuclei min surface-to-surface distance ----
    mind_hat = loo_predict(X, mind)
    m_mind = metrics(mind, mind_hat, ids)
    perzyg_b = [{"id": ids[i], "label": labels[i], "qc": qc[i],
                 "true_distance": round(float(mind[i]), 3),
                 "pred_distance": round(float(max(0.0, mind_hat[i])), 3),
                 "transcript_total": int(total[i])} for i in range(n)]
    payload_b = {
        "meta": {"model_version": MODEL_VERSION, "seed": SEED,
                 "trained_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
                 "n_zygotes": n, "input": "DAPI image features only",
                 "target": "pronuclei min surface-to-surface distance (pronuclei project metric)",
                 "evaluation": "nested leave-one-embryo-out CV"},
        "metrics": {"pron_min_distance_um": m_mind},
        "final_model": final_model(X, mind),
        "per_zygote": perzyg_b,
    }
    json.dump(payload_b, open(MODEL_B, "w"), indent=1)

    print(f"\n[A] image->distance_sum : MAE {m_dsum['mae']:.2f} µm  ρ {m_dsum['spearman']:.3f}")
    print(f"[A] ->clock-> tau       : MAE {m_tau_clock['mae']:.3f}  ρ {m_tau_clock['spearman']:.3f}  "
          f"pair-order {m_tau_clock['pair_order_accuracy']:.3f}  coverage {cov:.3f}")
    print(f"[A] direct image->tau   : MAE {m_tau_direct['mae']:.3f}  ρ {m_tau_direct['spearman']:.3f}")
    print(f"    baselines: constant {base_const:.3f}  single-feature {base_spread:.3f}")
    print(f"[B] image->min-distance : MAE {m_mind['mae']:.2f} µm  ρ {m_mind['spearman']:.3f}  "
          f"pair-order {m_mind['pair_order_accuracy']:.3f}")
    print(f"\nwrote {os.path.relpath(MODEL_A, HERE)} and {os.path.relpath(MODEL_B, HERE)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
