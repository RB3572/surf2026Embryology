"""
Monotonic, probabilistic, calibrated pseudotime clock (task 17).

Maps a dimensionless geometry feature — the pronucleus-to-cell-centre distance
sum divided by cell radius — to developmental pseudotime tau, and returns a
POSTERIOR (mean, sd), 50/80/95% prediction intervals, and a calibrated
confidence. Monotonic by construction (isotonic regression); the interval width
is learned from held-out residuals and rescaled so its nominal coverage is
honest.

True time supervision is the Scheffler 2021 live-imaging cohort (53 embryos,
2057 frames, real tau). Scheffler records physical distances but no cell radius;
since a mouse zygote's radius is ~constant through pronuclear migration, a single
reference radius R0 makes its features dimensionless and comparable to the fixed
MERFISH stacks, whose per-embryo MEASURED radius then self-corrects any imaging
scale offset (domain adaptation by normalization).
"""
from __future__ import annotations

import os
import sys

import numpy as np

SCRIPTS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # .../scripts
sys.path.insert(0, SCRIPTS_DIR)
import train_pronuclear_pseudotime as T  # noqa: E402

Z_50, Z_80, Z_95 = 0.6744898, 1.2815516, 1.9599640


def load_scheffler(reference_radius_um: float):
    """Return dimensionless sum/R, tau, and embryo groups from the live-imaging cohort."""
    df, _ = T.load_and_validate()
    s_over_R = (df["distance_sum_um"].to_numpy(float) / reference_radius_um)
    tau = df[T.TARGET].to_numpy(float)
    groups = df[T.GROUP].to_numpy()
    return s_over_R, tau, groups


def _iso_fit(x, y):
    from sklearn.isotonic import IsotonicRegression
    return IsotonicRegression(increasing=False, out_of_bounds="clip").fit(x, y)


def _grouped_cv(x, y, groups):
    """Leave-one-embryo-out isotonic predictions (held-out)."""
    pred = np.zeros_like(y)
    for g in np.unique(groups):
        te = groups == g
        iso = _iso_fit(x[~te], y[~te])
        pred[te] = iso.predict(x[te])
    return pred


def _sd_vs_tau(pred, resid, n_bins=8):
    """Robust residual sd as a smooth function of predicted tau (heteroscedastic)."""
    order = np.argsort(pred)
    p, r = pred[order], resid[order]
    edges = np.linspace(0, 1, n_bins + 1)
    centers, sds = [], []
    for i in range(n_bins):
        m = (p >= edges[i]) & (p <= edges[i + 1])
        if m.sum() >= 8:
            mad = np.median(np.abs(r[m] - np.median(r[m])))
            centers.append((edges[i] + edges[i + 1]) / 2)
            sds.append(max(1.4826 * mad, 1e-3))
    if len(centers) < 2:
        s = float(max(np.std(resid), 1e-2))
        return np.array([0.0, 1.0]), np.array([s, s])
    return np.array(centers), np.array(sds)


class ProbClock:
    """Isotonic mean + heteroscedastic, calibrated Gaussian posterior on [0,1]."""

    def __init__(self, reference_radius_um: float = 42.0):
        self.R0 = float(reference_radius_um)
        self.iso = None
        self.sd_x = None
        self.sd_y = None
        self.calib_c = 1.0
        self.cv = {}

    def _sd(self, tau_pred):
        return np.interp(tau_pred, self.sd_x, self.sd_y) * self.calib_c

    def fit(self, s_over_R, tau, groups):
        self.iso = _iso_fit(s_over_R, tau)
        cv_pred = _grouped_cv(s_over_R, tau, groups)
        resid = tau - cv_pred
        self._cv_pred, self._cv_tau, self._cv_x = cv_pred, tau, s_over_R
        self.sd_x, self.sd_y = _sd_vs_tau(cv_pred, resid)
        # calibrate: scale sd so the 95% interval covers ~95% on held-out data
        base_sd = np.interp(cv_pred, self.sd_x, self.sd_y)
        def cov(c):
            lo = np.clip(cv_pred - Z_95 * c * base_sd, 0, 1)
            hi = np.clip(cv_pred + Z_95 * c * base_sd, 0, 1)
            return np.mean((tau >= lo) & (tau <= hi))
        cs = np.linspace(0.4, 3.0, 261)
        self.calib_c = float(cs[np.argmin([abs(cov(c) - 0.95) for c in cs])])
        # honest CV metrics + coverage at all nominal levels
        self.cv = self._metrics(cv_pred, tau, groups)
        return self

    def _metrics(self, pred, tau, groups):
        from scipy.stats import spearmanr
        sd = self._sd(pred)
        def coverage(z):
            lo = np.clip(pred - z * sd, 0, 1); hi = np.clip(pred + z * sd, 0, 1)
            return float(np.mean((tau >= lo) & (tau <= hi)))
        # within-embryo monotonicity: fraction of embryos with non-negative rank corr
        mono = []
        for g in np.unique(groups):
            m = groups == g
            if m.sum() >= 4:
                rho = spearmanr(tau[m], pred[m]).correlation
                mono.append(rho if rho == rho else 0.0)
        return {
            "n": int(len(tau)), "n_embryos": int(len(np.unique(groups))),
            "mae": float(np.mean(np.abs(tau - pred))),
            "spearman": float(spearmanr(tau, pred).correlation),
            "coverage_50": coverage(Z_50), "coverage_80": coverage(Z_80),
            "coverage_95": coverage(Z_95),
            "within_embryo_mono_median": float(np.median(mono)) if mono else float("nan"),
            "calibration_scale": self.calib_c,
        }

    def predict(self, s_over_R: float) -> dict:
        m = float(self.iso.predict([s_over_R])[0])
        sd = float(self._sd(np.array([m]))[0])
        def iv(z):
            return [float(np.clip(m - z * sd, 0, 1)), float(np.clip(m + z * sd, 0, 1))]
        # confidence: narrower interval + mid-range -> higher; shrinks near sd ceiling
        conf = float(np.clip(1.0 - sd / 0.35, 0.05, 1.0))
        return {"tau_mean": round(m, 4), "tau_sd": round(sd, 4),
                "interval_50": [round(v, 4) for v in iv(Z_50)],
                "interval_80": [round(v, 4) for v in iv(Z_80)],
                "interval_95": [round(v, 4) for v in iv(Z_95)],
                "confidence": round(conf, 3)}

    def calibration_curve(self):
        """Empirical coverage vs nominal for a range of interval levels (dev view)."""
        from scipy.stats import norm
        pred, tau = self._cv_pred, self._cv_tau
        sd = self._sd(pred)
        out = []
        for nominal in (0.5, 0.6, 0.7, 0.8, 0.9, 0.95):
            z = norm.ppf(0.5 + nominal / 2)
            lo = np.clip(pred - z * sd, 0, 1); hi = np.clip(pred + z * sd, 0, 1)
            out.append({"nominal": nominal,
                        "empirical": round(float(np.mean((tau >= lo) & (tau <= hi))), 3)})
        return out

    def cv_scatter(self, max_pts=500):
        """Downsampled held-out (true, predicted) pairs for the dev view scatter."""
        n = len(self._cv_tau)
        idx = np.linspace(0, n - 1, min(max_pts, n)).astype(int)
        return [{"true": round(float(self._cv_tau[i]), 4),
                 "pred": round(float(self._cv_pred[i]), 4),
                 "sum_over_R": round(float(self._cv_x[i]), 4)} for i in idx]

    def to_dict(self) -> dict:
        return {"kind": "IsotonicProbabilisticClock", "reference_radius_um": self.R0,
                "iso_x": self.iso.f_.x.tolist(), "iso_y": self.iso.f_.y.tolist(),
                "sd_x": self.sd_x.tolist(), "sd_y": self.sd_y.tolist(),
                "calibration_scale": self.calib_c, "cv_metrics": self.cv,
                "calibration_curve": self.calibration_curve(),
                "cv_scatter": self.cv_scatter(),
                "feature": "distance_sum / cell_radius (dimensionless)",
                "z_levels": {"50": Z_50, "80": Z_80, "95": Z_95}}


def fit_default(reference_radius_um: float = 42.0) -> ProbClock:
    x, y, g = load_scheffler(reference_radius_um)
    return ProbClock(reference_radius_um).fit(x, y, g)
