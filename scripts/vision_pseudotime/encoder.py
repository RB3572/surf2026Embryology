"""
Exploratory (feasibility-only) direct image encoder (brief item 8).

THIS IS NOT A PRODUCTION MODEL. It is a harness that shows whether a direct
image -> tau regressor is even learnable, and it is wired to REFUSE to train
unless there are enough independent, timestamped embryos to mean anything.

Hard guards, straight from the brief and README:
  * training needs >= MIN_INDEPENDENT_EMBRYOS distinct biological embryos, each
    with a defensible per-frame timestamp (annotated pronuclear formation + NEBD).
    Movie 1 is ONE embryo, so it can never satisfy this.
  * perturbation / overshoot frames are refused as training data (OOD only).
  * tau labels must exist per frame; frame index alone is rejected.

The encoder features are deliberately rotation/reflection-invariant (radial
intensity rings about the embryo centroid), so a fitted head reads the
concentric separation geometry that carries time — not absolute pose and not a
corner timestamp. The head is a plain, seedless ridge (deterministic).

On the real material this module refuses (1 embryo, no per-frame tau). Its
end-to-end behaviour is demonstrated on synthetic multi-embryo phantoms with
known tau, clearly labelled synthetic, so the plumbing is proven without any
fabricated performance claim on real data.
"""
from __future__ import annotations

import numpy as np

MODEL_STATUS = "exploratory-feasibility"
MIN_INDEPENDENT_EMBRYOS = 3          # < this cannot support a held-out-embryo test
N_RINGS = 12


class InsufficientData(RuntimeError):
    """Raised when the material cannot support even a feasibility fit."""


# ───────────────────────────── rotation-invariant features ─────────────────────────────
def _centroid_xy(gray: np.ndarray):
    h, w = gray.shape
    ys, xs = np.mgrid[0:h, 0:w]
    tot = gray.sum()
    if tot <= 0:
        return (h - 1) / 2.0, (w - 1) / 2.0
    return float((ys * gray).sum() / tot), float((xs * gray).sum() / tot)


def radial_features(image: np.ndarray, n_rings: int = N_RINGS) -> np.ndarray:
    """
    Orientation-free descriptor: mean intensity in concentric rings about the
    embryo centroid, per channel, plus global occupancy. Invariant to rotation
    and reflection, so the head cannot exploit pose or a fixed-corner overlay.
    """
    a = np.asarray(image, np.float32)
    if a.ndim == 2:
        a = a[..., None]
    h, w, c = a.shape
    gray = a.mean(axis=2)
    cy, cx = _centroid_xy(gray)
    ys, xs = np.mgrid[0:h, 0:w]
    r = np.sqrt((ys - cy) ** 2 + (xs - cx) ** 2)
    rmax = r.max() if r.max() > 0 else 1.0
    edges = np.linspace(0, rmax, n_rings + 1)
    feats = []
    for ch in range(c):
        v = a[..., ch]
        for k in range(n_rings):
            m = (r >= edges[k]) & (r < edges[k + 1])
            feats.append(float(v[m].mean()) if m.any() else 0.0)
    feats.append(float((gray > gray.mean()).mean()))          # occupancy-ish
    feats.append(float(np.percentile(gray, 90)))
    return np.asarray(feats, np.float32)


# ───────────────────────────── the feasibility model ─────────────────────────────
class ExploratoryImageEncoder:
    """Handcrafted-feature ridge regressor. Deterministic. Feasibility only."""

    def __init__(self, n_rings: int = N_RINGS, alpha: float = 1.0):
        self.n_rings = n_rings
        self.alpha = alpha
        self.mean_ = None
        self.scale_ = None
        self.coef_ = None
        self.intercept_ = None
        self.status = MODEL_STATUS
        self.fitted_ = False

    def features(self, image) -> np.ndarray:
        return radial_features(image, self.n_rings)

    def fit(self, images, taus, groups, treatments=None):
        """
        Fit ONLY if the data can support a feasibility claim: enough independent
        embryos and real per-frame tau. Refuses perturbation frames.
        """
        taus = np.asarray(taus, float)
        groups = np.asarray(groups)
        n_emb = len(set(groups.tolist()))
        if n_emb < MIN_INDEPENDENT_EMBRYOS:
            raise InsufficientData(
                f"{n_emb} independent embryo(s) < required {MIN_INDEPENDENT_EMBRYOS}; "
                "a direct image clock cannot be trained or validated. Movie 1 is one embryo.")
        if not np.isfinite(taus).all():
            raise InsufficientData("per-frame tau labels missing; frame index alone is not tau.")
        if treatments is not None:
            bad = [t for t in treatments if _is_perturbation(t)]
            if bad:
                raise InsufficientData(
                    "perturbation/overshoot frames present in training set — OOD only.")
        X = np.stack([self.features(im) for im in images])
        self.mean_ = X.mean(0); self.scale_ = X.std(0) + 1e-8
        Z = (X - self.mean_) / self.scale_
        # closed-form ridge (deterministic, no solver randomness)
        d = Z.shape[1]
        A = Z.T @ Z + self.alpha * np.eye(d)
        self.coef_ = np.linalg.solve(A, Z.T @ (taus - taus.mean()))
        self.intercept_ = float(taus.mean())
        self.fitted_ = True
        return self

    def predict(self, image) -> float:
        if not self.fitted_:
            raise RuntimeError("encoder not fitted")
        z = (self.features(image) - self.mean_) / self.scale_
        return float(np.clip(self.intercept_ + z @ self.coef_, 0.0, 1.0))

    def predict_many(self, images) -> np.ndarray:
        return np.array([self.predict(im) for im in images], float)


def _is_perturbation(text: str) -> bool:
    t = (text or "").lower()
    return any(w in t for w in ("nocodazole", "overshoot", "overexpression",
                                "cytochalasin", "s25n", "droplet", "p150"))


# ───────────────────────────── synthetic phantom (labelled synthetic) ─────────────────────────────
def synthetic_embryo_frame(tau: float, seed: int, out: int = 96,
                           overlay: bool = False) -> np.ndarray:
    """
    A deterministic SYNTHETIC 2.5D phantom whose two 'pronuclei' move together as
    tau -> 1 (separation shrinks), mimicking the real time-bearing geometry. For
    proving the harness only — never presented as real data.

    overlay=True stamps a CONSTANT bright corner square (a fake burned-in
    timestamp) that is identical in every frame and therefore uncorrelated with
    tau. A model reading the overlay cannot predict tau; the radial encoder,
    centred on the embryo, ignores it. Used to exercise the anti-overlay
    saliency check.
    """
    r = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:out, 0:out].astype(np.float32)
    cy = cx = out / 2
    body = np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * (out * 0.33) ** 2)))
    sep = (1 - tau) * out * 0.28 + out * 0.05      # separation shrinks with tau
    ang = r.uniform(0, 2 * np.pi)                   # random pose (invariance test)
    dx, dy = np.cos(ang) * sep / 2, np.sin(ang) * sep / 2
    def blob(oy, ox):
        return np.exp(-(((yy - (cy + oy)) ** 2 + (xx - (cx + ox)) ** 2) / (2 * (out * 0.09) ** 2)))
    pron = blob(dy, dx) + blob(-dy, -dx)
    img = np.stack([body * 0.5 + pron, body, np.clip(pron, 0, 1)], axis=-1).astype(np.float32)
    img += r.normal(0, 0.01, img.shape).astype(np.float32)
    img = np.clip(img / img.max(), 0, 1)
    if overlay:
        s = max(4, out // 12)
        img[:s, :s, :] = 1.0                       # constant corner "timestamp"
    return img, (cy, cx)


OVERLAY_BOX_FRAC = 1.0 / 12.0                        # matches the stamped corner size


def overlay_box(out: int):
    s = max(4, int(out * OVERLAY_BOX_FRAC))
    return (0, 0, s, s)                              # y0, x0, y1, x1


def synthetic_dataset(n_embryos: int, frames_per: int, seed: int = 0, overlay: bool = False):
    """Multiple synthetic embryos, each a tau sweep. Returns images, taus, groups."""
    images, taus, groups = [], [], []
    for e in range(n_embryos):
        for f in range(frames_per):
            tau = f / (frames_per - 1)
            img, _ = synthetic_embryo_frame(tau, seed=seed * 1000 + e * 100 + f, overlay=overlay)
            images.append(img); taus.append(tau); groups.append(f"synthetic_emb{e}")
    return images, np.array(taus, float), np.array(groups)
