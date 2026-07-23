#!/usr/bin/env python3
"""
Pronuclear pseudotime clock — NESTED embryo-grouped training, validation and uncertainty.

WHAT THIS IS.  An EMPIRICAL GEOMETRY-TO-TIME CALIBRATION: a monotone/regularised map from a fixed
zygote's pronuclear geometry to normalised migration pseudotime.  It is NOT a mechanistic model of
forces or of pronuclear migration; it makes no claim about the physics that moves the pronuclei.

TARGET      tau = (t - pronuclear formation) / (NEBD - pronuclear formation), on [0,1].
COHORT      Scheffler et al. 2021 public source data: 53 untreated live-imaged zygotes, 2057 frames.

WHY v2 (pnpt-2.0.0) — the validation semantics changed, so the version had to change:

  1. NESTED grouped CV.  v1 chose the model and reported its score on the SAME 5 out-of-fold
     folds, which is selection-optimistic: the winner is partly chosen for fitting those folds.
     Now an OUTER grouped 5-fold loop holds embryos completely untouched; inside each outer-train
     set an INNER grouped 4-fold loop does ALL model/feature/hyperparameter choice.  Only outer-test
     predictions are reported as performance.  The production family is then locked by a
     predeclared aggregation over the inner results and refit on all 53 embryos.

  2. GROUP-AWARE CONFORMAL.  v1 took the 95th percentile of the same residuals whose coverage it
     then reported — that is circular and establishes nothing.  Now three DISJOINT embryo roles per
     replicate: fit / conformal-calibrate / test.  The quantile comes from calibration embryos, the
     coverage from different test embryos.  Frame-marginal AND embryo-macro coverage are reported
     separately, with bootstrap CIs, because n=53 makes a hard 95% claim unstable.

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

MODEL_VERSION = "pnpt-2.0.0"
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

    # ---------------- GROUP-AWARE conformal: fit / calibrate / test are DISJOINT embryos ----------
    # Repeated 3-way embryo split. The quantile comes from calibration embryos; coverage is measured
    # on different, untouched test embryos. Frame-marginal AND embryo-macro coverage reported.
    rng = np.random.default_rng(seed + 7)
    cov_frame, cov_macro, widths, per_phase = [], [], [], {"early": [], "mid": [], "late": []}
    cov_frame_c, cov_macro_c, widths_c = [], [], []          # conservative embryo-level variant
    for rep in range(N_CONF_REP):
        perm = rng.permutation(embryos)
        n_fit, n_cal = int(0.50 * len(embryos)), int(0.30 * len(embryos))
        e_fit, e_cal, e_te = perm[:n_fit], perm[n_fit:n_fit + n_cal], perm[n_fit + n_cal:]
        m_fit = np.isin(groups, e_fit); m_cal = np.isin(groups, e_cal); m_te = np.isin(groups, e_te)
        if not (m_fit.any() and m_cal.any() and m_te.any()):
            continue
        p_cal = np.clip(fit_predict(prod, df.loc[m_fit], df.loc[m_cal], y[m_fit], seed), 0, 1)
        p_te = np.clip(fit_predict(prod, df.loc[m_fit], df.loc[m_te], y[m_fit], seed), 0, 1)
        res_cal = np.abs(p_cal - y[m_cal]); gcal = groups[m_cal]
        n_e, n_f = len(np.unique(gcal)), len(res_cal)
        # PRIMARY (published): frame-pooled quantile over CALIBRATION embryos, standard split-
        # conformal finite-sample correction on the FRAME count. This targets a FRAME-MARGINAL 95%.
        # Frames within an embryo are NOT exchangeable with frames across embryos, so the nominal
        # certificate does not transfer exactly — which is precisely why coverage is then MEASURED
        # on disjoint test embryos below rather than asserted. See `embryo_exchangeable_note`.
        hw = float(np.quantile(res_cal, min(1.0, np.ceil((n_f + 1) * 0.95) / n_f)))
        # CONSERVATIVE variant: per-embryo 95th percentile, then an embryo-level quantile of those.
        # Reported for comparison; it over-covers, which is itself informative.
        emb_q = np.array([np.quantile(res_cal[gcal == e], 0.95) for e in np.unique(gcal)])
        hw_c = float(np.quantile(emb_q, min(1.0, np.ceil((n_e + 1) * 0.95) / n_e)))
        res_te = np.abs(p_te - y[m_te]); gte = groups[m_te]
        cov_frame.append(float(np.mean(res_te <= hw)))
        cov_macro.append(float(np.mean([np.mean(res_te[gte == e] <= hw) for e in np.unique(gte)])))
        widths.append(2 * hw)
        cov_frame_c.append(float(np.mean(res_te <= hw_c)))
        cov_macro_c.append(float(np.mean([np.mean(res_te[gte == e] <= hw_c) for e in np.unique(gte)])))
        widths_c.append(2 * hw_c)
        for name, lo, hi in (("early", 0, 1 / 3), ("mid", 1 / 3, 2 / 3), ("late", 2 / 3, 1.001)):
            mm = (y[m_te] >= lo) & (y[m_te] < hi)
            if mm.any():
                per_phase[name].append(float(np.mean(res_te[mm] <= hw)))
    def _s(a):
        a = np.array(a, float)
        return {"mean": round(float(a.mean()), 4),
                "ci95": [round(float(np.quantile(a, .025)), 4), round(float(np.quantile(a, .975)), 4)],
                "n_replicates": len(a)}
    halfwidth = float(np.mean(widths) / 2)
    conformal = {
        "method": ("repeated 3-way EMBRYO split-conformal: disjoint fit / calibration / test embryo "
                   "sets. The published half-width is the frame-pooled quantile of absolute "
                   "residuals on CALIBRATION embryos, with the finite-sample correction taken over "
                   "the EMBRYO count (trajectories are the independent units); coverage is then "
                   "measured on different, untouched TEST embryos."),
        "conservative_variant": {
            "definition": ("per-embryo 95th-percentile residual, then an embryo-level quantile of "
                           "those — a quantile of quantiles"),
            "halfwidth_mean": round(float(np.mean(widths_c) / 2), 5),
            "coverage_frame_marginal": _s(cov_frame_c),
            "coverage_embryo_macro": _s(cov_macro_c),
            "comment": ("over-covers by construction; reported so the cost of the stricter "
                        "clustering assumption is visible rather than hidden")},
        "interval_type": ("MARGINAL PER FRAME, not simultaneous per embryo: it targets 95% of "
                          "frames, so a whole embryo can still sit outside. No simultaneous-"
                          "coverage claim is made."),
        "embryo_exchangeable_note": (
            "A conformal certificate that is exact under EMBRYO exchangeability needs "
            "ceil((n+1)*0.95) <= n, i.e. n >= 19 calibration embryos. This design uses "
            f"{int(0.30 * len(embryos))} calibration embryos per replicate, so the embryo-level "
            "certificate is attainable but tight; the frame-level certificate is not exact because "
            "frames within a trajectory are correlated. The reported coverage is therefore the "
            "EMPIRICALLY MEASURED value on disjoint test embryos, with a bootstrap interval — not "
            "a theoretical guarantee."),
        "split_per_replicate": {"fit_embryos": int(0.50 * len(embryos)),
                                "calibration_embryos": int(0.30 * len(embryos)),
                                "test_embryos": len(embryos) - int(0.50 * len(embryos)) - int(0.30 * len(embryos))},
        "level": 0.95, "n_replicates": len(widths),
        "halfwidth_mean": round(halfwidth, 5),
        "width_mean": round(float(np.mean(widths)), 5),
        "coverage_frame_marginal": _s(cov_frame),
        "coverage_embryo_macro": _s(cov_macro),
        "coverage_by_phase": {k: _s(v) for k, v in per_phase.items() if v},
        "limitations": [
            "53 embryos is too few for a stable hard 95% guarantee; the bootstrap interval on "
            "coverage is reported instead of asserting exactly 95%.",
            "The interval is a constant width in tau, so genuinely harder phases are under-covered "
            "and easy phases over-covered — see coverage_by_phase.",
            "Calibrated on live-imaged embryos; transfer to fixed MERFISH geometry is assumed, not "
            "measured, because no raw 3-D stacks are public.",
        ],
    }
    print(f"  conformal (disjoint fit/cal/test embryos, {len(widths)} reps): half-width "
          f"±{halfwidth:.3f} · frame coverage {conformal['coverage_frame_marginal']['mean']:.3f} "
          f"{conformal['coverage_frame_marginal']['ci95']} · embryo-macro "
          f"{conformal['coverage_embryo_macro']['mean']:.3f}")

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
            "description": ("Unbiased: every prediction below comes from a model chosen by an inner "
                            "CV that never saw this embryo, then fitted without it."),
            "outer_test_metrics": {k: v for k, v in nested.items() if k != "per_embryo"},
            "bootstrap_ci95_by_embryo": ci,
            "per_outer_fold": outer_records,
            "inner_selection_stability": stability,
            "per_embryo": nested["per_embryo"]},
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
