#!/usr/bin/env python3
"""
Pronuclear pseudotime clock — NESTED embryo-grouped training, validation and uncertainty.

WHAT THIS IS.  An EMPIRICAL GEOMETRY-TO-TIME CALIBRATION: a monotone/regularised map from a fixed
zygote's pronuclear geometry to normalised migration pseudotime.  It is NOT a mechanistic model of
forces or of pronuclear migration; it makes no claim about the physics that moves the pronuclei.

TARGET      tau = (t - pronuclear formation) / (NEBD - pronuclear formation), on [0,1].
COHORT      Scheffler et al. 2021 public source data: 53 untreated live-imaged zygotes, 2057 frames.

WHY v3 (pnpt-3.0.0) — interval semantics and attenuation scenarios changed again:

  A. THE PRIMARY INTERVAL IS NOT CONFORMAL AND NO LONGER CLAIMS TO BE.  v2 applied a frame-count
     conformal correction while its own metadata said the correction was over the EMBRYO count, and
     asserted that a finite exact 95% cluster certificate was "attainable" with 15 calibration
     embryos when the threshold is 19.  Both were wrong.  The primary interval is now named an
     EMPIRICAL DISJOINT-EMBRYO PREDICTION INTERVAL: a plain 95th percentile on calibration embryos,
     evaluated on different test embryos, with bootstrap uncertainty and NO formal guarantee.
  B. A genuinely rigorous CLUSTER-LEVEL conformal interval is reported separately, using ONE
     nonconformity score per calibration embryo (max |residual|), which protects SIMULTANEOUS
     per-embryo coverage.  With 53 embryos the valid quantile is the MAXIMUM score; the resulting
     width is honestly wide and reported as such.
  C. The nested score estimates the INNER-SELECTION PROCEDURE (whose chosen family varies by fold),
     not the locked production family.  A separate, clearly post-selection fixed-family grouped-CV
     sensitivity is now reported for the locked family.

WHY v2 — the validation semantics changed, so the version had to change:

  1. NESTED grouped CV.  v1 chose the model and reported its score on the SAME 5 out-of-fold
     folds, which is selection-optimistic: the winner is partly chosen for fitting those folds.
     Now an OUTER grouped 5-fold loop holds embryos completely untouched; inside each outer-train
     set an INNER grouped 4-fold loop does ALL model/feature/hyperparameter choice.  Only outer-test
     predictions are reported as performance.  The production family is then locked by a
     predeclared aggregation over the inner results and refit on all 53 embryos.

  2. DISJOINT-EMBRYO INTERVAL.  v1 took the 95th percentile of the same residuals whose coverage it
     then reported — circular, establishing nothing.  Now three DISJOINT embryo roles per replicate:
     fit / calibrate / test.  The width comes from calibration embryos, the coverage from different
     test embryos, with bootstrap CIs.  (v3 corrects how this is NAMED — see A above.)

  3. EXPLICIT PAIR COUNTING.  v1 reported (Kendall tau-b + 1)/2 as "pairwise ordering accuracy".
     With a step-function predictor that is wrong — tau-b absorbs ties into a correction rather than
     counting them.  Now concordant / discordant / tied-prediction / tied-truth are counted
     explicitly under a documented convention, pooled AND macro-within-embryo.

HARD CONSTRAINTS (enforced in code, checked by scripts/test_pronuclear_pseudotime.py):
  * No transcript / gene / probe / expression quantity is read anywhere in this file.
  * Deployable models may not use male_*/female_* (fixed zygotes lack a reliable identity call) or
    any *volume* column (published relative volumes are normalised to each pronucleus's own FUTURE
    endpoint volume: unmeasurable in a snapshot and directly leaky).
  * Every preprocessing step, direction choice and isotonic knot is fitted inside its training split.

Usage:  python3 scripts/train_pronuclear_pseudotime.py [--seed 20260723] [--quick]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_CSV = os.path.join(HERE, "calibration_data", "scheffler2021",
                       "scheffler_2021_control_zygote_trajectories.csv")
OUT_DIR = os.path.join(HERE, "data", "pseudotime_calibration")

MODEL_VERSION = "pnpt-3.0.0"
DATA_VERSION = "scheffler2021-v1"
SEED = 20260723
N_OUTER, N_INNER = 5, 4
TARGET, GROUP = "normalized_time_tau", "embryo_id"
EXPECTED_EMBRYOS, EXPECTED_ROWS = 53, 2057

SYMMETRIC_FEATURES = ["nearer_to_center_um", "farther_to_center_um",
                      "distance_sum_um", "distance_difference_um"]
BANNED_SUBSTRINGS = ("volume", "male_", "female_", "transcript", "gene", "probe",
                     "expression", "count", "time_h", "duration", "tau")

# Predeclared selection rule, fixed BEFORE any result was seen. Applied identically inside every
# inner loop and for the final production lock.
SELECTION_RULE = ("lowest inner-CV macro MAE (mean of per-embryo MAE, so each embryo counts once); "
                  "ties within 0.002 broken by higher strict pairwise ordering accuracy, then by "
                  "lower complexity rank")
TIE_TOL = 0.002
COMPLEXITY_RANK = {"linear_sum": 0, "linear_nearer": 0, "isotonic_sum": 1, "isotonic_nearer": 1,
                   "monospline_sum": 2, "ridge_core": 3, "ridge_symmetric": 4, "hgb_symmetric": 5}


def assert_deployable(features, name):
    for f in features:
        for bad in BANNED_SUBSTRINGS:
            if bad in f.lower():
                raise ValueError(f"model '{name}' requests non-deployable feature '{f}' ('{bad}')")


# ═══════════════════════════════════════ data ═══════════════════════════════════════
def load_and_validate(path=SRC_CSV):
    if not os.path.isfile(path):
        raise SystemExit(f"calibration CSV missing: {path}")
    df = pd.read_csv(path)
    problems, notes = [], []
    n_emb, n_row = df[GROUP].nunique(), len(df)
    if n_emb != EXPECTED_EMBRYOS:
        problems.append(f"expected {EXPECTED_EMBRYOS} embryos, found {n_emb}")
    if n_row != EXPECTED_ROWS:
        problems.append(f"expected {EXPECTED_ROWS} rows, found {n_row}")
    required = set(SYMMETRIC_FEATURES) | {GROUP, TARGET, "time_h", "migration_duration_h"}
    missing = sorted(required - set(df.columns))
    if missing:
        problems.append(f"missing columns: {missing}")
    if not problems:
        for f in SYMMETRIC_FEATURES:
            if df[f].isna().any():
                problems.append(f"feature '{f}' has {int(df[f].isna().sum())} missing values")
        bad_mono = [e for e, g in df.groupby(GROUP) if not g["time_h"].is_monotonic_increasing]
        if bad_mono:
            problems.append(f"non-monotonic time_h in {len(bad_mono)} embryos")
        if df[TARGET].min() < 0 or df[TARGET].max() > 1:
            problems.append("tau outside [0,1]")
        if not np.allclose(df["time_h"] / df["migration_duration_h"], df[TARGET], atol=1e-9):
            problems.append("tau != time_h / migration_duration_h")
    leaky = [c for c in df.columns if "volume" in c.lower() or c.startswith(("male_", "female_"))]
    n_leaky_na = int(df[leaky].isna().sum().sum()) if leaky else 0
    if n_leaky_na:
        notes.append(f"{n_leaky_na} missing values in non-deployable columns (never used)")
    if problems:
        raise SystemExit("SOURCE VALIDATION FAILED:\n  - " + "\n  - ".join(problems))
    with open(path, "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    per = df.groupby(GROUP).size()
    return df, {
        "source_csv": os.path.relpath(path, HERE), "sha256": sha, "data_version": DATA_VERSION,
        "n_embryos": int(n_emb), "n_frames": int(n_row),
        "frames_per_embryo": {"min": int(per.min()), "max": int(per.max()),
                              "median": float(per.median())},
        "frame_interval_h": 0.25,
        "migration_duration_h": {"min": round(float(df["migration_duration_h"].min()), 3),
                                 "max": round(float(df["migration_duration_h"].max()), 3),
                                 "median": round(float(df["migration_duration_h"].median()), 3)},
        "tau_range": [round(float(df[TARGET].min()), 4), round(float(df[TARGET].max()), 4)],
        "deployable_features": SYMMETRIC_FEATURES, "excluded_columns": leaky,
        "exclusion_reason": ("male_*/female_* need pronuclear identity, absent in fixed zygotes; "
                             "*volume* columns are normalised to each pronucleus's own future "
                             "endpoint volume, so they leak elapsed time"),
        "validation": {"checks_passed": [
            f"{EXPECTED_EMBRYOS} embryos", f"{EXPECTED_ROWS} rows", "required columns present",
            "no missing values in deployable features", "time_h monotonic within embryo",
            "tau within [0,1]", "tau == time_h / migration_duration_h"], "notes": notes},
    }


# ══════════════════════════════ explicit ordering metrics ══════════════════════════════
def ordering_counts(y_true, y_pred, tol=0.0):
    """Count EVERY unordered pair explicitly. No tau-b shortcut.

    A pair is `comparable` when the truth is not tied (tied-truth pairs carry no ordering
    information and are excluded from the denominator). Among comparable pairs a prediction is
    concordant / discordant / tied.

    CONVENTION (documented, and asserted by a hand-checked test):
      strict_accuracy = concordant / comparable        -- prediction ties score 0, the conservative
                                                          reading, because a tie fails to order.
      tie_rate        = tied_pred / comparable
      inversion_rate  = discordant / comparable
      so strict_accuracy + tie_rate + inversion_rate == 1 exactly.
    `half_credit_accuracy` = (concordant + 0.5*tied_pred)/comparable is also reported for
    comparability with tau-b-style figures, but is NOT the headline number."""
    y_true = np.asarray(y_true, float); y_pred = np.asarray(y_pred, float)
    n = len(y_true)
    if n < 2:
        return {"comparable": 0, "concordant": 0, "discordant": 0, "tied_pred": 0, "tied_true": 0,
                "strict_accuracy": None, "half_credit_accuracy": None,
                "tie_rate": None, "inversion_rate": None}
    c = d = tp = tt = 0
    for i in range(n - 1):
        dt = y_true[i + 1:] - y_true[i]
        dp = y_pred[i + 1:] - y_pred[i]
        true_tie = np.abs(dt) <= tol
        tt += int(true_tie.sum())
        m = ~true_tie
        if not m.any():
            continue
        dtm, dpm = dt[m], dp[m]
        pred_tie = np.abs(dpm) <= tol
        tp += int(pred_tie.sum())
        agree = (np.sign(dtm) == np.sign(dpm)) & ~pred_tie
        c += int(agree.sum())
        d += int((~agree & ~pred_tie).sum())
    comp = c + d + tp
    return {"comparable": comp, "concordant": c, "discordant": d, "tied_pred": tp, "tied_true": tt,
            "strict_accuracy": (c / comp) if comp else None,
            "half_credit_accuracy": ((c + 0.5 * tp) / comp) if comp else None,
            "tie_rate": (tp / comp) if comp else None,
            "inversion_rate": (d / comp) if comp else None}


def rnd(d, k=5):
    return {kk: (round(v, k) if isinstance(v, float) else v) for kk, v in d.items()}


def metric_block(y_true, y_pred, groups):
    """Pooled + macro (per-embryo) metrics. Macro is primary: one vote per embryo, so a long
    trajectory cannot dominate, and the near-monotone within-embryo sequence cannot inflate the
    headline pooled ordering number."""
    y_true = np.asarray(y_true, float); y_pred = np.asarray(y_pred, float)
    groups = np.asarray(groups)
    err = np.abs(y_pred - y_true)
    per = {}
    for e in np.unique(groups):
        m = groups == e
        if m.sum() < 2:
            continue
        rho = spearmanr(y_true[m], y_pred[m]).statistic
        oc = ordering_counts(y_true[m], y_pred[m])
        per[str(e)] = {"n": int(m.sum()), "mae": round(float(err[m].mean()), 5),
                       "rmse": round(float(np.sqrt(((y_pred[m] - y_true[m]) ** 2).mean())), 5),
                       "bias": round(float((y_pred[m] - y_true[m]).mean()), 5),
                       "spearman": None if not np.isfinite(rho) else round(float(rho), 5),
                       "strict_ordering": None if oc["strict_accuracy"] is None
                       else round(oc["strict_accuracy"], 5),
                       "tie_rate": None if oc["tie_rate"] is None else round(oc["tie_rate"], 5)}
    macro_mae = float(np.mean([v["mae"] for v in per.values()])) if per else float("nan")
    ords = [v["strict_ordering"] for v in per.values() if v["strict_ordering"] is not None]
    rhos = [v["spearman"] for v in per.values() if v["spearman"] is not None]
    pooled = ordering_counts(y_true, y_pred)
    r = float(np.corrcoef(y_true, y_pred)[0, 1]) if len(y_true) > 2 else float("nan")
    return {
        "macro_mae": round(macro_mae, 5), "pooled_mae": round(float(err.mean()), 5),
        "median_ae": round(float(np.median(err)), 5),
        "pooled_rmse": round(float(np.sqrt(((y_pred - y_true) ** 2).mean())), 5),
        "pearson_r": round(r, 5), "pearson_r2": round(r * r, 5),
        "pooled_spearman": round(float(spearmanr(y_true, y_pred).statistic), 5),
        "pooled_kendall_tau_b": round(float(kendalltau(y_true, y_pred).statistic), 5),
        "pooled_ordering": rnd(pooled),
        "macro_within_embryo_ordering": round(float(np.mean(ords)), 5) if ords else None,
        "macro_within_embryo_spearman": round(float(np.mean(rhos)), 5) if rhos else None,
        "per_embryo": per,
    }


# ═════════════════════════════════════ candidates ═════════════════════════════════════
class LinearOnFeature:
    def __init__(self, feature): self.feature = feature
    def fit(self, X, y):
        v = np.asarray(X[self.feature], float)
        self.b_, self.a_ = np.polyfit(v, y, 1); return self
    def predict(self, X):
        return self.a_ + self.b_ * np.asarray(X[self.feature], float)


class IsotonicOnFeature:
    """Monotone step calibration of ONE feature. Direction learned on the training split only."""
    def __init__(self, feature): self.feature = feature
    def fit(self, X, y):
        v = np.asarray(X[self.feature], float)
        self.sign_ = -1.0 if np.corrcoef(v, y)[0, 1] < 0 else 1.0
        self.iso_ = IsotonicRegression(out_of_bounds="clip", increasing=True).fit(self.sign_ * v, y)
        return self
    def predict(self, X):
        return self.iso_.predict(self.sign_ * np.asarray(X[self.feature], float))


class MonotoneSplineOnFeature:
    """SMOOTH monotone alternative to isotonic: I-spline basis (SplineTransformer with
    extrapolation='linear') + non-negative-coefficient least squares, so the fit is monotone in the
    learned direction but continuous — no step ties. Added in v2 because isotonic's ties are a real
    liability for ordering a snapshot cohort. It competes under the SAME inner-CV rule as every
    other candidate; it is not privileged."""
    def __init__(self, feature, n_knots=6, degree=3):
        self.feature, self.n_knots, self.degree = feature, n_knots, degree
    def fit(self, X, y):
        v = np.asarray(X[self.feature], float).reshape(-1, 1)
        self.sign_ = -1.0 if np.corrcoef(v.ravel(), y)[0, 1] < 0 else 1.0
        v = self.sign_ * v
        self.st_ = SplineTransformer(n_knots=self.n_knots, degree=self.degree,
                                     extrapolation="linear", include_bias=False).fit(v)
        B = np.cumsum(self.st_.transform(v)[:, ::-1], axis=1)[:, ::-1]   # I-spline (monotone) basis
        from scipy.optimize import nnls
        A = np.hstack([np.ones((len(B), 1)), B])
        lo = y.min()
        coef, _ = nnls(A[:, 1:], y - lo)                                  # non-negative => monotone
        self.coef_, self.intercept_ = coef, lo
        return self
    def predict(self, X):
        v = self.sign_ * np.asarray(X[self.feature], float).reshape(-1, 1)
        B = np.cumsum(self.st_.transform(v)[:, ::-1], axis=1)[:, ::-1]
        return self.intercept_ + B @ self.coef_


def ridge_pipe(features, alpha=1.0):
    return Pipeline([("impute", SimpleImputer(strategy="median")),
                     ("scale", StandardScaler()), ("model", Ridge(alpha=alpha))])


def hgb_pipe(seed):
    return Pipeline([("impute", SimpleImputer(strategy="median")),
                     ("model", HistGradientBoostingRegressor(
                         max_depth=3, max_iter=200, learning_rate=0.06, min_samples_leaf=40,
                         l2_regularization=1.0, early_stopping=False, random_state=seed))])


def candidates(seed):
    core = ["nearer_to_center_um", "farther_to_center_um"]
    return [
        {"key": "linear_sum", "label": "Linear · distance sum", "features": ["distance_sum_um"],
         "complexity": "2 parameters", "note": "transparent one-feature baseline",
         "make": lambda s: LinearOnFeature("distance_sum_um")},
        {"key": "linear_nearer", "label": "Linear · nearer distance",
         "features": ["nearer_to_center_um"], "complexity": "2 parameters",
         "note": "one-feature baseline on the nearer pronucleus",
         "make": lambda s: LinearOnFeature("nearer_to_center_um")},
        {"key": "isotonic_sum", "label": "Isotonic · distance sum", "features": ["distance_sum_um"],
         "complexity": "monotone step function",
         "note": "monotone, assumes direction not constant speed; produces prediction TIES",
         "make": lambda s: IsotonicOnFeature("distance_sum_um")},
        {"key": "isotonic_nearer", "label": "Isotonic · nearer distance",
         "features": ["nearer_to_center_um"], "complexity": "monotone step function",
         "note": "monotone step calibration of the nearer pronucleus",
         "make": lambda s: IsotonicOnFeature("nearer_to_center_um")},
        {"key": "monospline_sum", "label": "Monotone spline · distance sum",
         "features": ["distance_sum_um"], "complexity": "I-spline, 6 knots, non-negative coefs",
         "note": "smooth monotone alternative to isotonic — continuous, no step ties",
         "make": lambda s: MonotoneSplineOnFeature("distance_sum_um")},
        {"key": "ridge_core", "label": "Ridge · near+far", "features": core,
         "complexity": "2 features, L2", "note": "regularised linear on the two sorted distances",
         "make": lambda s: ridge_pipe(core)},
        {"key": "ridge_symmetric", "label": "Ridge · symmetric distances",
         "features": SYMMETRIC_FEATURES, "complexity": "4 features, L2",
         "note": "adds the (redundant) sum and difference",
         "make": lambda s: ridge_pipe(SYMMETRIC_FEATURES)},
        {"key": "hgb_symmetric", "label": "Gradient boosting · symmetric distances",
         "features": SYMMETRIC_FEATURES, "complexity": "depth-3 trees, 200 iters",
         "note": "the one nonlinear candidate", "make": lambda s: hgb_pipe(s)},
    ]


def leaky_candidate():
    """NON-DEPLOYABLE upper bound. Never enters selection; reported separately in the UI."""
    feats = SYMMETRIC_FEATURES + ["male_relative_volume", "female_relative_volume",
                                  "volume_sum", "volume_difference"]
    return {"key": "LEAKY_ridge_volumes", "features": feats,
            "label": "Ridge · distances + endpoint-normalised volumes",
            "complexity": "8 features, L2", "deployable": False,
            "note": ("NON-DEPLOYABLE upper bound: relative volumes are normalised to each "
                     "pronucleus's own future endpoint volume, so they encode elapsed time and "
                     "cannot be measured in a fixed snapshot."),
            "make": lambda s: ridge_pipe(feats)}


def fit_predict(spec, tr_df, te_df, y_tr, seed):
    m = spec["make"](seed)
    m.fit(tr_df[spec["features"]], y_tr)
    return np.asarray(m.predict(te_df[spec["features"]]), float)


def pick(scores):
    """Apply SELECTION_RULE to {key: {macro_mae, strict_ordering}}. Deterministic."""
    best = min(v["macro_mae"] for v in scores.values())
    near = [k for k, v in scores.items() if v["macro_mae"] <= best + TIE_TOL]
    return sorted(near, key=lambda k: (-(scores[k]["strict_ordering"] or 0),
                                       COMPLEXITY_RANK.get(k, 99), k))[0]


# ═══════════════════════════════ grouped folds ═══════════════════════════════
def grouped_folds(df, embryos, n_splits, seed):
    """Deterministic embryo-grouped folds. Returns fold index per ROW."""
    rng = np.random.default_rng(seed)
    shuffled = np.array(embryos)[rng.permutation(len(embryos))]
    order = {e: i for i, e in enumerate(shuffled)}
    d = df.assign(_o=df[GROUP].map(order)).sort_values(["_o"])
    idx = np.full(len(df), -1)
    g = d[GROUP].to_numpy()
    for k, (_, te) in enumerate(GroupKFold(n_splits=n_splits).split(d, d[TARGET], g)):
        idx[d.index.to_numpy()[te]] = k
    assert (idx >= 0).all()
    return idx


def boot_ci(values, groups, stat, seed, B=2000, level=0.95):
    """Embryo-level (cluster) bootstrap CI. Resamples EMBRYOS with replacement, never frames, so
    the CI reflects the 53 independent trajectories rather than 2057 correlated frames."""
    rng = np.random.default_rng(seed)
    embs = np.unique(groups)
    idx_of = {e: np.flatnonzero(groups == e) for e in embs}
    out = []
    for _ in range(B):
        pick_e = rng.choice(embs, size=len(embs), replace=True)
        sel = np.concatenate([idx_of[e] for e in pick_e])
        gg = np.concatenate([np.full(len(idx_of[e]), f"{e}#{i}") for i, e in enumerate(pick_e)])
        try:
            v = stat(sel, gg)
        except Exception:                                              # noqa: BLE001
            continue
        if v is not None and np.isfinite(v):
            out.append(v)
    if not out:
        return None
    a = (1 - level) / 2
    return [round(float(np.quantile(out, a)), 5), round(float(np.quantile(out, 1 - a)), 5)]


# ═════════════════════════════════════ main ═════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--quick", action="store_true", help="fewer bootstrap/conformal replicates")
    args = ap.parse_args()
    seed = args.seed
    B_BOOT = 400 if args.quick else 2000
    N_CONF_REP = 8 if args.quick else 40
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"pronuclear pseudotime · {MODEL_VERSION} · seed {seed}")
    print("  EMPIRICAL geometry-to-time calibration (not a mechanistic migration model)")
    df, summary = load_and_validate()
    df = df.reset_index(drop=True)
    y = df[TARGET].to_numpy(float)
    groups = df[GROUP].to_numpy()
    embryos = np.array(sorted(df[GROUP].unique()))
    print(f"  validated {summary['n_embryos']} embryos / {summary['n_frames']} frames "
          f"(sha256 {summary['sha256'][:16]})")

    specs = candidates(seed)
    for s in specs:
        assert_deployable(s["features"], s["key"])

    # ---------------- NESTED grouped CV ----------------
    outer = grouped_folds(df, embryos, N_OUTER, seed)
    oof = np.full(len(df), np.nan)          # outer-test predictions ONLY (unbiased)
    oof_raw = np.full(len(df), np.nan)
    outer_records, inner_choice = [], []
    for ko in range(N_OUTER):
        te_o = outer == ko
        tr_o = ~te_o
        tr_embs = np.array(sorted(pd.unique(groups[tr_o])))
        sub = df.loc[tr_o].reset_index(drop=True)
        inner = grouped_folds(sub, tr_embs, N_INNER, seed + 100 + ko)
        ysub = sub[TARGET].to_numpy(float)
        gsub = sub[GROUP].to_numpy()

        # inner loop: model selection ONLY, never touching the outer-test embryos
        inner_scores = {}
        for sp in specs:
            p = np.full(len(sub), np.nan)
            for ki in range(N_INNER):
                m_te = inner == ki
                p[m_te] = fit_predict(sp, sub.loc[~m_te], sub.loc[m_te], ysub[~m_te], seed)
            pc = np.clip(p, 0, 1)
            mb = metric_block(ysub, pc, gsub)
            inner_scores[sp["key"]] = {"macro_mae": mb["macro_mae"],
                                       "strict_ordering": mb["pooled_ordering"]["strict_accuracy"]}
        chosen = pick(inner_scores)
        inner_choice.append(chosen)
        sp = next(s for s in specs if s["key"] == chosen)
        pred = fit_predict(sp, df.loc[tr_o], df.loc[te_o], y[tr_o], seed)
        oof_raw[te_o] = pred
        oof[te_o] = np.clip(pred, 0, 1)
        mb = metric_block(y[te_o], oof[te_o], groups[te_o])
        outer_records.append({
            "outer_fold": ko, "chosen_by_inner_cv": chosen,
            "chosen_label": sp["label"],
            "n_train_embryos": int(len(tr_embs)), "n_test_embryos": int(pd.unique(groups[te_o]).size),
            "n_train_frames": int(tr_o.sum()), "n_test_frames": int(te_o.sum()),
            "test_embryos": sorted(pd.unique(groups[te_o]).tolist()),
            "inner_scores": {k: rnd(v) for k, v in inner_scores.items()},
            "outer_test_metrics": {k: v for k, v in mb.items() if k != "per_embryo"},
        })
        print(f"  outer fold {ko + 1}/{N_OUTER}: inner-CV chose {chosen:<16} "
              f"outer-test macroMAE={mb['macro_mae']:.4f}")

    nested = metric_block(y, oof, groups)
    print(f"\n  NESTED outer-test (unbiased): macroMAE={nested['macro_mae']:.4f}  "
          f"r={nested['pearson_r']:.3f}  rho={nested['pooled_spearman']:.3f}  "
          f"strict-order={nested['pooled_ordering']['strict_accuracy']:.3f}  "
          f"ties={nested['pooled_ordering']['tie_rate']:.3f}")

    # embryo-bootstrap CIs on the nested outer-test estimate
    def _mae(sel, gg):
        e = np.abs(oof[sel] - y[sel])
        return float(np.mean([e[gg == u].mean() for u in np.unique(gg)]))
    def _r(sel, gg): return float(np.corrcoef(y[sel], oof[sel])[0, 1])
    def _r2(sel, gg): return float(np.corrcoef(y[sel], oof[sel])[0, 1] ** 2)
    def _rho(sel, gg): return float(spearmanr(y[sel], oof[sel]).statistic)
    def _ord(sel, gg): return ordering_counts(y[sel], oof[sel])["strict_accuracy"]
    def _macro_ord(sel, gg):
        vals = [ordering_counts(y[sel][gg == u], oof[sel][gg == u])["strict_accuracy"]
                for u in np.unique(gg)]
        vals = [v for v in vals if v is not None]
        return float(np.mean(vals)) if vals else None
    ci = {"macro_mae": boot_ci(None, groups, _mae, seed, B_BOOT),
          "pearson_r": boot_ci(None, groups, _r, seed, B_BOOT),
          "pearson_r2": boot_ci(None, groups, _r2, seed, B_BOOT),
          "pooled_spearman": boot_ci(None, groups, _rho, seed, B_BOOT),
          "pooled_strict_ordering": boot_ci(None, groups, _ord, seed, B_BOOT),
          "macro_within_embryo_ordering": boot_ci(None, groups, _macro_ord, seed, B_BOOT)}
    print("  embryo-bootstrap 95% CIs: " + ", ".join(
        f"{k}={v}" for k, v in ci.items() if v))

    # ---------------- production lock: predeclared aggregation of INNER results ----------------
    agg = {}
    for sp in specs:
        maes = [r["inner_scores"][sp["key"]]["macro_mae"] for r in outer_records]
        ords = [r["inner_scores"][sp["key"]]["strict_ordering"] or 0 for r in outer_records]
        agg[sp["key"]] = {"macro_mae": float(np.mean(maes)),
                          "strict_ordering": float(np.mean(ords))}
    prod_key = pick(agg)
    prod = next(s for s in specs if s["key"] == prod_key)
    print(f"\n  PRODUCTION LOCK (mean inner-CV across outer folds, {SELECTION_RULE.split(';')[0]}):"
          f" {prod['label']}")
    print("    inner-CV means: " + "  ".join(
        f"{k}={v['macro_mae']:.4f}" for k, v in sorted(agg.items(), key=lambda x: x[1]['macro_mae'])))

    # honest disclosure when the outer folds disagreed with the locked family
    stability = {k: inner_choice.count(k) for k in set(inner_choice)}

    # ───── SECONDARY: fixed-family grouped CV for the LOCKED production family ─────
    # The nested estimate above scores the SELECTION PROCEDURE, whose chosen family varied across
    # outer folds. It is NOT a direct unbiased score of the locked family. This block runs a plain
    # grouped 5-fold CV holding the family FIXED at the production choice, which answers "how does
    # isotonic-on-distance-sum itself behave out-of-fold?" — but it is POST-SELECTION and therefore
    # DESCRIPTIVE ONLY: the family was chosen using these same 53 embryos, so this number is
    # optimistically biased and must never be quoted as the headline performance.
    ff_pred = np.full(len(df), np.nan)
    for ko in range(N_OUTER):
        te_o = outer == ko
        ff_pred[te_o] = np.clip(fit_predict(prod, df.loc[~te_o], df.loc[te_o], y[~te_o], seed), 0, 1)
    ff_m = metric_block(y, ff_pred, groups)
    def _ff_mae(sel, gg):
        e = np.abs(ff_pred[sel] - y[sel])
        return float(np.mean([e[gg == u].mean() for u in np.unique(gg)]))
    fixed_family = {
        "label": prod["label"], "key": prod_key,
        "status": "POST-SELECTION — descriptive only, NOT an unbiased estimate",
        "why": ("The production family was locked using all 53 embryos (via the inner-CV aggregate), "
                "so a grouped CV of that same family on the same embryos is optimistically biased. "
                "The unbiased development-time number is the nested SELECTION-PROCEDURE estimate."),
        "protocol": f"plain grouped {N_OUTER}-fold CV, family held fixed at the production choice",
        "metrics": {k: v for k, v in ff_m.items() if k != "per_embryo"},
        "bootstrap_ci95_macro_mae": boot_ci(None, groups, _ff_mae, seed, B_BOOT),
        "difference_vs_nested_macro_mae": round(ff_m["macro_mae"] - nested["macro_mae"], 5),
    }
    print(f"  fixed-family sensitivity (POST-SELECTION, {prod['label']}): "
          f"macroMAE={ff_m['macro_mae']:.4f} "
          f"(nested selection-procedure estimate was {nested['macro_mae']:.4f})")

    # ─────────── UNCERTAINTY: an EMPIRICAL interval + a RIGOROUS cluster-level sensitivity ───────
    # HONEST FRAMING (corrected in v3). The primary interval below is NOT a formal conformal
    # certificate and is not described as one. Frames within a trajectory are strongly correlated,
    # so frame-level exchangeability — which any frame-count conformal correction assumes — does not
    # hold. What we actually do is: estimate a width on CALIBRATION embryos, evaluate coverage on
    # DIFFERENT test embryos, and report the measured coverage with bootstrap uncertainty. That is
    # an empirical disjoint-embryo prediction interval, nothing stronger.
    #
    # A genuinely valid cluster-level statement needs ONE nonconformity score per embryo. With n_cal
    # calibration embryos the conformal index is k = ceil((n_cal+1)*(1-alpha)); a finite exact 95%
    # statement requires k <= n_cal, i.e. n_cal >= 19. Moreover k == n_cal for every n_cal in
    # [19, 38], so with 53 embryos the valid quantile IS THE MAXIMUM of the per-embryo scores. That
    # is rigorous but very wide, and it is reported as such rather than hidden.
    ALPHA = 0.05
    N_CAL_MIN = int(np.ceil(ALPHA ** -1)) - 1          # smallest n with ceil((n+1)(1-a)) <= n  -> 19
    rng = np.random.default_rng(seed + 7)

    # --- primary: empirical disjoint-embryo interval (fit 50% / calibrate 30% / test 20%) ---
    n_fit, n_cal = int(0.50 * len(embryos)), int(0.30 * len(embryos))
    n_te = len(embryos) - n_fit - n_cal
    cov_frame, cov_macro, widths, per_phase = [], [], [], {"early": [], "mid": [], "late": []}
    for _ in range(N_CONF_REP):
        perm = rng.permutation(embryos)
        e_fit, e_cal, e_te = perm[:n_fit], perm[n_fit:n_fit + n_cal], perm[n_fit + n_cal:]
        m_fit, m_cal, m_te = (np.isin(groups, e) for e in (e_fit, e_cal, e_te))
        if not (m_fit.any() and m_cal.any() and m_te.any()):
            continue
        p_cal = np.clip(fit_predict(prod, df.loc[m_fit], df.loc[m_cal], y[m_fit], seed), 0, 1)
        p_te = np.clip(fit_predict(prod, df.loc[m_fit], df.loc[m_te], y[m_fit], seed), 0, 1)
        res_cal = np.abs(p_cal - y[m_cal])
        # PLAIN empirical 95th percentile of the pooled calibration-embryo residuals. No conformal
        # finite-sample correction is applied, because none would be valid under clustering; calling
        # it a "frame-count corrected conformal quantile" (as v2 did) was the error being fixed.
        hw = float(np.quantile(res_cal, 1 - ALPHA))
        res_te = np.abs(p_te - y[m_te]); gte = groups[m_te]
        cov_frame.append(float(np.mean(res_te <= hw)))
        cov_macro.append(float(np.mean([np.mean(res_te[gte == e] <= hw) for e in np.unique(gte)])))
        widths.append(2 * hw)
        for name, lo, hi in (("early", 0, 1 / 3), ("mid", 1 / 3, 2 / 3), ("late", 2 / 3, 1.001)):
            mm = (y[m_te] >= lo) & (y[m_te] < hi)
            if mm.any():
                per_phase[name].append(float(np.mean(res_te[mm] <= hw)))

    # --- sensitivity: RIGOROUS cluster-level split conformal, one score per calibration embryo ---
    # Score = max |residual| over the embryo's frames  =>  the interval is SIMULTANEOUS over all
    # frames of a new embryo: P(every frame of an exchangeable new embryo is covered) >= 1 - alpha.
    # Split forced to n_cal >= N_CAL_MIN so the quantile index is attainable.
    n_cal_r = max(N_CAL_MIN, int(0.36 * len(embryos)))
    n_fit_r = int(0.45 * len(embryos))
    n_te_r = len(embryos) - n_fit_r - n_cal_r
    k_idx = int(np.ceil((n_cal_r + 1) * (1 - ALPHA)))
    feasible = k_idx <= n_cal_r
    uses_max = feasible and k_idx == n_cal_r
    cov_simul, cov_frame_r, widths_r = [], [], []
    if feasible and n_te_r >= 3:
        for _ in range(N_CONF_REP):
            perm = rng.permutation(embryos)
            e_fit, e_cal, e_te = (perm[:n_fit_r], perm[n_fit_r:n_fit_r + n_cal_r],
                                  perm[n_fit_r + n_cal_r:])
            m_fit, m_cal, m_te = (np.isin(groups, e) for e in (e_fit, e_cal, e_te))
            if not (m_fit.any() and m_cal.any() and m_te.any()):
                continue
            p_cal = np.clip(fit_predict(prod, df.loc[m_fit], df.loc[m_cal], y[m_fit], seed), 0, 1)
            p_te = np.clip(fit_predict(prod, df.loc[m_fit], df.loc[m_te], y[m_fit], seed), 0, 1)
            rc, gc = np.abs(p_cal - y[m_cal]), groups[m_cal]
            scores = np.sort(np.array([rc[gc == e].max() for e in np.unique(gc)]))   # ONE per embryo
            hw_r = float(scores[k_idx - 1])                                          # k-th smallest
            rt, gt = np.abs(p_te - y[m_te]), groups[m_te]
            cov_simul.append(float(np.mean([rt[gt == e].max() <= hw_r for e in np.unique(gt)])))
            cov_frame_r.append(float(np.mean(rt <= hw_r)))
            widths_r.append(2 * hw_r)

    def _s(a):
        a = np.array(a, float)
        if not len(a):
            return None
        return {"mean": round(float(a.mean()), 4),
                "ci95": [round(float(np.quantile(a, .025)), 4), round(float(np.quantile(a, .975)), 4)],
                "n_replicates": len(a)}

    halfwidth = float(np.mean(widths) / 2)
    hw_rig = float(np.mean(widths_r) / 2) if widths_r else None
    conformal = {
        "interval_name": "empirical disjoint-embryo prediction interval",
        "formal_guarantee": False,
        "method": (
            "Repeated 3-way EMBRYO split (disjoint fit / calibration / test). The half-width is the "
            "PLAIN empirical 95th percentile of pooled absolute residuals on the CALIBRATION "
            "embryos; coverage is then measured on different, untouched TEST embryos and reported "
            "with a bootstrap interval. NO finite-sample conformal correction is applied and NO "
            "formal 95% guarantee is claimed: frames within a trajectory are correlated, so the "
            "frame-level exchangeability such a correction assumes does not hold."),
        "why_not_conformal": (
            "A frame-count conformal correction would assume frames are exchangeable, which is false "
            "here (adjacent frames of one embryo are near-duplicates). A valid CLUSTER-level "
            f"conformal statement needs one score per embryo and n_cal >= {N_CAL_MIN}; this primary "
            f"split uses n_cal = {n_cal}, which is below that threshold. The rigorous version is "
            "reported separately below."),
        "interval_type": (
            "MARGINAL PER FRAME and EMPIRICAL: it targets 95% of frames on average across "
            "replicates. It is not simultaneous over an embryo's frames and carries no coverage "
            "certificate."),
        "applies_to_fixed_cohort_as": (
            "EMPIRICAL TRANSFER UNCERTAINTY. Applying this width to a fixed MERFISH snapshot assumes "
            "the live-imaged residual distribution transfers to a different segmentation pipeline "
            "and a shifted geometry cohort. That assumption is unverified, so the fixed-cohort "
            "interval is an empirical indication of spread, NOT a coverage guarantee."),
        "estimator_note": (
            "The published half-width is the MEAN of the per-replicate widths, applied to the final "
            "model refitted on all 53 embryos. Averaging widths across random splits and reusing "
            "them for a differently-fitted model is a reasonable empirical summary but is not "
            "itself a split-conformal construction."),
        "level": 0.95, "n_replicates": len(widths),
        "quantile_computed_over": "pooled calibration-embryo FRAMES (plain empirical percentile)",
        "split_per_replicate": {"fit_embryos": n_fit, "calibration_embryos": n_cal,
                                "test_embryos": n_te},
        "halfwidth_mean": round(halfwidth, 5),
        "width_mean": round(float(np.mean(widths)), 5),
        "coverage_frame_marginal": _s(cov_frame),
        "coverage_embryo_macro": _s(cov_macro),
        "coverage_by_phase": {k: _s(v) for k, v in per_phase.items() if v},
        "rigorous_cluster_sensitivity": {
            "name": "cluster-conformal construction sensitivity",
            "what_is_reported": (
                "the MEAN half-width over independent per-split cluster-conformal constructions, "
                "plus the coverage measured on each split's held-out embryos"),
            # The distinction the aggregate hides: EACH construction is a valid conformal interval;
            # their AVERAGE is not, and neither is applying that average to a model refitted on all
            # 53 embryos (which saw every calibration embryo).
            "per_split_construction_has_finite_sample_guarantee": bool(feasible),
            "reported_mean_halfwidth_has_formal_guarantee": False,
            "guarantee_scope": (
                f"Each of the {N_CONF_REP} individual constructions is a valid split-conformal "
                "interval: with one nonconformity score per calibration embryo and "
                f"k = {k_idx} <= n_cal = {n_cal_r}, it attains >= 95% SIMULTANEOUS per-embryo "
                "coverage. The MEAN of those widths reported below is a summary statistic, NOT a "
                "conformal interval, and attaching the guarantee to it would be an overclaim."),
            "assumptions": [
                "embryos are exchangeable (one lab, one protocol, untreated controls — plausible "
                "within this cohort, untested beyond it)",
                "the model is fitted independently of the calibration embryos (enforced by the "
                "disjoint fit/calibration/test split)",
                "the guarantee is per-construction and does NOT transfer to the final model "
                "refitted on all 53 embryos, which has seen every calibration embryo",
            ],
            "nonconformity_score": "max |residual| over all frames of a calibration embryo",
            "protects": ("per construction: P(EVERY frame of a new exchangeable embryo lies inside "
                         "the interval) >= 0.95, i.e. SIMULTANEOUS per-embryo coverage — a strictly "
                         "stronger statement than the frame-marginal primary interval"),
            "split_per_replicate": {"fit_embryos": n_fit_r, "calibration_embryos": n_cal_r,
                                    "test_embryos": n_te_r},
            "quantile_index_k": k_idx, "n_calibration_embryos": n_cal_r,
            "uses_maximum_score": bool(uses_max),
            "feasibility_note": (
                f"k = ceil((n_cal+1)*0.95) = {k_idx} and n_cal = {n_cal_r}. "
                + (f"k == n_cal, so the valid quantile IS THE MAXIMUM of the {n_cal_r} per-embryo "
                   "scores. That is exact but deliberately conservative: with 53 embryos, k < n_cal "
                   "would require n_cal >= 39, leaving too few embryos to fit and test. The width "
                   "below is therefore honestly wide."
                   if uses_max else "k < n_cal, so a sub-maximal order statistic is used.")),
            "halfwidth_mean": round(hw_rig, 5) if hw_rig is not None else None,
            "halfwidth_per_split_range": ([round(float(np.min(widths_r) / 2), 5),
                                           round(float(np.max(widths_r) / 2), 5)]
                                          if widths_r else None),
            "width_mean": round(float(np.mean(widths_r)), 5) if widths_r else None,
            "n_constructions": len(widths_r),
            "coverage_simultaneous_per_embryo": _s(cov_simul),
            "coverage_frame_marginal": _s(cov_frame_r),
            "width_ratio_vs_primary": (round(hw_rig / halfwidth, 3)
                                       if (hw_rig and halfwidth) else None),
        },
        "limitations": [
            "The primary interval has NO formal coverage guarantee — it is an empirical "
            "disjoint-embryo interval whose coverage is measured, not certified.",
            "53 embryos is too few for a non-conservative cluster-level conformal quantile: any "
            f"feasible n_cal in [{N_CAL_MIN}, 38] forces the MAXIMUM per-embryo score.",
            "The cluster-conformal guarantee holds for EACH individual split construction, not for "
            "the averaged width reported, and not for the production model refitted on all 53 "
            "embryos — that model has seen every calibration embryo, so no split-conformal "
            "statement applies to it.",
            "The interval is constant width in tau, so genuinely harder phases are under-covered "
            "and easy phases over-covered — see coverage_by_phase.",
            "Estimated on live-imaged embryos; transfer to fixed MERFISH geometry is assumed, not "
            "measured, because no raw 3-D stacks are public.",
        ],
    }
    print(f"  interval (EMPIRICAL, disjoint fit/cal/test, {len(widths)} reps): half-width "
          f"±{halfwidth:.3f} · frame coverage {conformal['coverage_frame_marginal']['mean']:.3f} "
          f"{conformal['coverage_frame_marginal']['ci95']}")
    if hw_rig is not None:
        rc_ = conformal["rigorous_cluster_sensitivity"]
        print(f"  rigorous cluster conformal (n_cal={n_cal_r}, k={k_idx}"
              f"{', = MAX' if uses_max else ''}): half-width ±{hw_rig:.3f} "
              f"({rc_['width_ratio_vs_primary']}x wider) · simultaneous per-embryo coverage "
              f"{rc_['coverage_simultaneous_per_embryo']['mean']:.3f}")

    # ---------------- refit production on ALL embryos (only after selection + evaluation) --------
    final = prod["make"](seed)
    final.fit(df[prod["features"]], y)

    def describe(m, feats):
        d = {"kind": type(m).__name__, "features": feats}
        if isinstance(m, LinearOnFeature):
            d.update(form="tau = a + b*x", a=round(float(m.a_), 8), b=round(float(m.b_), 8))
        elif isinstance(m, IsotonicOnFeature):
            d.update(form="tau = isotonic(sign*x)", sign=float(m.sign_),
                     knots_x=[round(float(v), 6) for v in m.iso_.X_thresholds_],
                     knots_y=[round(float(v), 6) for v in m.iso_.y_thresholds_])
        elif isinstance(m, MonotoneSplineOnFeature):
            d.update(form="tau = intercept + Ispline(sign*x) . coef  (coef >= 0 => monotone)",
                     sign=float(m.sign_), intercept=round(float(m.intercept_), 8),
                     coef=[round(float(v), 8) for v in m.coef_],
                     spline={"n_knots": m.n_knots, "degree": m.degree,
                             "knots": [round(float(v), 6) for v in m.st_.bsplines_[0].t],
                             "extrapolation": "linear"})
        elif isinstance(m, Pipeline):
            st = {}
            if "impute" in m.named_steps:
                st["imputer_median"] = [round(float(v), 8) for v in m.named_steps["impute"].statistics_]
            if "scale" in m.named_steps:
                sc = m.named_steps["scale"]
                st["scaler_mean"] = [round(float(v), 8) for v in sc.mean_]
                st["scaler_scale"] = [round(float(v), 8) for v in sc.scale_]
            mm = m.named_steps["model"]
            if isinstance(mm, Ridge):
                st.update(form="tau = intercept + coef . ((x-mean)/scale)", alpha=float(mm.alpha),
                          coef=[round(float(v), 8) for v in np.ravel(mm.coef_)],
                          intercept=round(float(mm.intercept_), 8))
            else:
                st.update(form="gradient-boosted trees (not analytically serialisable)")
            d["pipeline"] = st
        return d
    spec_json = describe(final, prod["features"])

    # ---------------- leaky upper bound, scored the SAME nested way (never selectable) -----------
    lk = leaky_candidate()
    lk_oof = np.full(len(df), np.nan)
    for ko in range(N_OUTER):
        te_o = outer == ko
        lk_oof[te_o] = np.clip(fit_predict(lk, df.loc[~te_o], df.loc[te_o], y[~te_o], seed), 0, 1)
    lk_m = metric_block(y, lk_oof, groups)

    # ---------------- artifacts ----------------
    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    fold_of = {str(e): int(outer[groups == e][0]) for e in embryos}
    for e in embryos:
        assert len(set(outer[groups == e].tolist())) == 1, f"embryo {e} spans outer folds"

    resid = oof - y
    frames = [{
        "embryo_id": str(groups[i]), "outer_fold": int(outer[i]),
        "time_h": round(float(df["time_h"].iat[i]), 4),
        "tau_true": round(float(y[i]), 5), "tau_pred": round(float(oof[i]), 5),
        "tau_pred_raw": round(float(oof_raw[i]), 5), "residual": round(float(resid[i]), 5),
        "lo95": round(float(max(0.0, oof[i] - halfwidth)), 5),
        "hi95": round(float(min(1.0, oof[i] + halfwidth)), 5),
        "covered": bool(abs(resid[i]) <= halfwidth),
        "out_of_range": bool(oof_raw[i] < 0 or oof_raw[i] > 1),
        "nearer_um": round(float(df["nearer_to_center_um"].iat[i]), 4),
        "farther_um": round(float(df["farther_to_center_um"].iat[i]), 4),
        "migration_duration_h": round(float(df["migration_duration_h"].iat[i]), 3),
    } for i in range(len(df))]

    traj = []
    for e in embryos:
        m = groups == e
        traj.append({"embryo_id": str(e), "outer_fold": fold_of[str(e)], "n": int(m.sum()),
                     "migration_duration_h": round(float(df.loc[m, "migration_duration_h"].iat[0]), 3),
                     "time_h": [round(float(v), 4) for v in df.loc[m, "time_h"]],
                     "tau_true": [round(float(v), 5) for v in y[m]],
                     "tau_pred": [round(float(v), 5) for v in oof[m]],
                     "nearer_um": [round(float(v), 3) for v in df.loc[m, "nearer_to_center_um"]],
                     "farther_um": [round(float(v), 3) for v in df.loc[m, "farther_to_center_um"]]})

    feature_stats = {f: {"min": round(float(df[f].min()), 4), "max": round(float(df[f].max()), 4),
                         "mean": round(float(df[f].mean()), 4), "std": round(float(df[f].std()), 4),
                         "p01": round(float(df[f].quantile(.01)), 4),
                         "p99": round(float(df[f].quantile(.99)), 4)} for f in SYMMETRIC_FEATURES}

    models_out = []
    for sp in specs:
        maes = [r["inner_scores"][sp["key"]]["macro_mae"] for r in outer_records]
        ords = [r["inner_scores"][sp["key"]]["strict_ordering"] for r in outer_records]
        models_out.append({
            "key": sp["key"], "label": sp["label"], "features": sp["features"],
            "n_features": len(sp["features"]), "complexity": sp["complexity"],
            "deployable": True, "note": sp["note"],
            "inner_cv_macro_mae_mean": round(float(np.mean(maes)), 5),
            "inner_cv_macro_mae_per_outer": [round(v, 5) for v in maes],
            "inner_cv_strict_ordering_mean": round(float(np.mean([o or 0 for o in ords])), 5),
            "chosen_in_n_outer_folds": inner_choice.count(sp["key"]),
            "selected_for_production": sp["key"] == prod_key})
    models_out.append({
        "key": lk["key"], "label": lk["label"], "features": lk["features"],
        "n_features": len(lk["features"]), "complexity": lk["complexity"], "deployable": False,
        "note": lk["note"], "selected_for_production": False,
        "nested_outer_test_metrics": {k: v for k, v in lk_m.items() if k != "per_embryo"}})

    artifact = {
        "meta": {
            "model_version": MODEL_VERSION, "data_version": DATA_VERSION, "seed": seed,
            "trained_at_utc": now, "n_outer_folds": N_OUTER, "n_inner_folds": N_INNER,
            "split_strategy": ("NESTED GroupKFold on embryo_id: outer folds are untouched "
                               "evaluation embryos; inner folds inside each outer-train set do all "
                               "model/feature/hyperparameter selection"),
            "target": "tau = (t - pronuclear formation) / (NEBD - pronuclear formation)",
            "selection_rule": SELECTION_RULE,
            "model_class": ("empirical geometry-to-time calibration; NOT a mechanistic model of "
                            "forces or pronuclear migration"),
            "sklearn_version": __import__("sklearn").__version__,
            "changes_from_v2": [
                "the primary interval is renamed an EMPIRICAL DISJOINT-EMBRYO PREDICTION INTERVAL "
                "and no longer claims any conformal guarantee; v2 mixed a frame-count correction "
                "with embryo-count metadata and wrongly said a 15-embryo calibration set could give "
                "a finite exact 95% cluster certificate (the threshold is 19)",
                "a rigorous cluster-level split-conformal sensitivity is added, using one "
                "nonconformity score per calibration embryo and protecting SIMULTANEOUS per-embryo "
                "coverage; with 53 embryos its valid quantile is the maximum score, so it is wide",
                "the nested estimate is explicitly labelled as the performance of the SELECTION "
                "PROCEDURE, with a separate post-selection fixed-family sensitivity for the locked "
                "production family",
            ],
            "changes_from_v1": [
                "nested outer/inner grouped CV replaces single-level OOF used for both selection "
                "and reporting (v1 was selection-optimistic)",
                "group-aware split-conformal with disjoint fit/calibration/test embryos replaces "
                "the circular same-residual quantile-and-coverage",
                "explicit pair counting replaces (Kendall tau-b + 1)/2 as ordering accuracy",
                "added a smooth monotone spline candidate that competes under the same inner-CV rule",
            ]},
        "dataset": summary,
        "folds": {"outer_assignment": fold_of,
                  "outer_sizes": [{"fold": k,
                                   "n_embryos": sum(1 for v in fold_of.values() if v == k),
                                   "n_frames": int((outer == k).sum())} for k in range(N_OUTER)]},
        "nested_evaluation": {
            "estimates": ("the performance of the COMPLETE INNER-SELECTION PROCEDURE, not of the "
                          "locked production family: the family chosen by the inner CV varied "
                          "across outer folds, so this is the expected performance of 'run the "
                          "selection procedure, then predict'"),
            "description": ("Unbiased for that procedure: every prediction comes from a model "
                            "chosen by an inner CV that never saw this embryo, then fitted without "
                            "it. See fixed_family_sensitivity for the locked family (post-selection)."),
            "outer_test_metrics": {k: v for k, v in nested.items() if k != "per_embryo"},
            "bootstrap_ci95_by_embryo": ci,
            "per_outer_fold": outer_records,
            "inner_selection_stability": stability,
            "per_embryo": nested["per_embryo"]},
        "fixed_family_sensitivity": fixed_family,
        "models": models_out,
        "production": {
            "key": prod_key, "label": prod["label"], "features": prod["features"],
            "spec": spec_json, "locked_by": SELECTION_RULE,
            "inner_cv_aggregate": {k: rnd(v) for k, v in agg.items()},
            "feature_stats": feature_stats,
            "note": ("Refit on all 53 embryos AFTER nested evaluation and selection. Its honest "
                     "performance is the nested outer-test estimate above, not a refit score.")},
        "uncertainty": conformal,
        "residual_diagnostics": {
            "overall_mae": round(float(np.abs(resid).mean()), 5),
            "overall_bias": round(float(resid.mean()), 5),
            "by_phase": [{"phase": n, "tau_range": [lo, round(hi, 3)],
                          "n": int(((y >= lo) & (y < hi)).sum()),
                          "mae": round(float(np.abs(resid[(y >= lo) & (y < hi)]).mean()), 5),
                          "bias": round(float(resid[(y >= lo) & (y < hi)].mean()), 5)}
                         for n, lo, hi in (("early", 0, 1 / 3), ("mid", 1 / 3, 2 / 3), ("late", 2 / 3, 1.001))]},
    }

    with open(os.path.join(OUT_DIR, "calibration.json"), "w") as fh:
        json.dump(artifact, fh, indent=1)
    with gzip.open(os.path.join(OUT_DIR, "oof_frames.json.gz"), "wt") as fh:
        json.dump({"model_version": MODEL_VERSION, "halfwidth_95": round(halfwidth, 5),
                   "note": "outer-test (nested) predictions only", "frames": frames}, fh,
                  separators=(",", ":"))
    with gzip.open(os.path.join(OUT_DIR, "trajectories.json.gz"), "wt") as fh:
        json.dump({"model_version": MODEL_VERSION, "embryos": traj}, fh, separators=(",", ":"))
    with open(os.path.join(OUT_DIR, "model.json"), "w") as fh:
        json.dump({"model_version": MODEL_VERSION, "data_version": DATA_VERSION,
                   "trained_at_utc": now, "seed": seed, "selected_key": prod_key,
                   "label": prod["label"], "features": prod["features"], "spec": spec_json,
                   "halfwidth_95": round(halfwidth, 5), "feature_stats": feature_stats,
                   "nested_outer_test": {k: v for k, v in nested.items() if k != "per_embryo"},
                   "model_class": "empirical geometry-to-time calibration"}, fh, indent=1)

    for f in ("calibration.json", "oof_frames.json.gz", "trajectories.json.gz", "model.json"):
        p = os.path.join(OUT_DIR, f)
        print(f"  wrote {os.path.relpath(p, HERE)} ({os.path.getsize(p) / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
