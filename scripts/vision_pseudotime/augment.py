"""
Deterministic, geometry-preserving augmentation (brief item 5).

Augmentation regularizes a model; it does NOT manufacture new independent
embryos (see the model card). Every transform here is chosen to leave the
time-bearing geometry — the pronucleus-to-cell-centre distances — either
unchanged or globally (isotropically) rescaled, so an augmented view of an
embryo still corresponds to the same tau.

ALLOWED (default): rotation, reflection, translation, mild ISOTROPIC scaling,
intensity/gamma, additive noise, blur, simulated photobleaching, z-slice
dropout (as a channel/projection-window effect), projection-window change.

DISABLED (by default, and refused unless explicitly forced with a biological
justification): anisotropic stretch, shear, elastic deformation, independent
pronuclear motion — each changes the distance geometry that encodes time.

Determinism: augment(img, seed) is a pure function of (img, seed, params).
The same seed reproduces the same output byte-for-byte.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, rotate, shift, zoom

ALLOWED = ("rotate", "reflect", "translate", "scale_isotropic", "intensity", "gamma",
           "noise", "blur", "bleach", "zdropout", "proj_window")
# transforms that alter the distance geometry encoding time — never on by default
DISABLED_GEOMETRY_BREAKING = ("anisotropic_stretch", "shear", "elastic", "pronuclear_motion")


class GeometryBreakingTransform(ValueError):
    """Raised if a caller asks for a transform that would corrupt the tau geometry."""


def default_params() -> dict:
    return {
        "rotate": {"max_deg": 180.0},        # full in-plane rotation is geometry-preserving
        "reflect": {"p": 0.5},
        "translate": {"max_frac": 0.06},
        "scale_isotropic": {"range": [0.92, 1.08]},   # SAME factor on x and y
        "intensity": {"range": [0.85, 1.15]},
        "gamma": {"range": [0.8, 1.25]},
        "noise": {"sigma": 0.02},
        "blur": {"max_sigma": 0.8},
        "bleach": {"range": [0.75, 1.0]},    # global multiplicative dimming
        "zdropout": {"p_channel": 0.15},     # drop the seg/percentile channel occasionally
        "proj_window": {"shift_frac": 0.1},  # small gamma-like reweight of a channel
    }


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(int(seed) & 0xFFFFFFFF)


def augment(img: np.ndarray, seed: int, params: dict | None = None,
            enabled: tuple | None = None) -> tuple[np.ndarray, dict]:
    """
    Apply a deterministic sequence of geometry-preserving transforms.

    img: HxW or HxWxC float array in [0,1]. Returns (aug_img, applied_record).
    A label channel (if present) is transformed with the SAME geometric ops
    (nearest-neighbour) so pixel correspondence — the label — is preserved.
    """
    p = params or default_params()
    enabled = enabled or ALLOWED
    for t in enabled:
        if t in DISABLED_GEOMETRY_BREAKING:
            raise GeometryBreakingTransform(
                f"{t} changes the distance geometry that encodes time; it is disabled by default. "
                "Enable only with a biological justification and recomputed targets.")
        if t not in ALLOWED:
            raise ValueError(f"unknown transform {t!r}")
    r = _rng(seed)
    a = np.asarray(img, np.float32).copy()
    has_c = a.ndim == 3
    applied = {}

    # ---- geometric (isotropic only) ----
    if "rotate" in enabled:
        deg = float(r.uniform(-p["rotate"]["max_deg"], p["rotate"]["max_deg"]))
        a = rotate(a, deg, axes=(0, 1), reshape=False, order=1, mode="reflect")
        applied["rotate_deg"] = round(deg, 3)
    if "reflect" in enabled and r.random() < p["reflect"]["p"]:
        a = a[:, ::-1, ...].copy()
        applied["reflect"] = True
    if "translate" in enabled:
        h, w = a.shape[:2]
        mf = p["translate"]["max_frac"]
        dy = float(r.uniform(-mf, mf) * h); dx = float(r.uniform(-mf, mf) * w)
        sh = (dy, dx, 0) if has_c else (dy, dx)
        a = shift(a, sh, order=1, mode="reflect")
        applied["translate_px"] = [round(dy, 2), round(dx, 2)]
    if "scale_isotropic" in enabled:
        s = float(r.uniform(*p["scale_isotropic"]["range"]))     # ONE factor -> isotropic
        a = _zoom_center(a, s)
        applied["scale_isotropic"] = round(s, 4)

    # ---- photometric (do not touch geometry) ----
    if "intensity" in enabled:
        f = float(r.uniform(*p["intensity"]["range"]))
        a = np.clip(a * f, 0, 1); applied["intensity"] = round(f, 4)
    if "gamma" in enabled:
        g = float(r.uniform(*p["gamma"]["range"]))
        a = np.clip(a, 0, 1) ** g; applied["gamma"] = round(g, 4)
    if "bleach" in enabled:
        b = float(r.uniform(*p["bleach"]["range"]))
        a = np.clip(a * b, 0, 1); applied["bleach"] = round(b, 4)
    if "blur" in enabled:
        sg = float(r.uniform(0, p["blur"]["max_sigma"]))
        if sg > 0:
            sigma = (sg, sg, 0) if has_c else (sg, sg)
            a = gaussian_filter(a, sigma=sigma); applied["blur_sigma"] = round(sg, 4)
    if "noise" in enabled:
        sd = p["noise"]["sigma"]
        a = np.clip(a + r.normal(0, sd, a.shape).astype(np.float32), 0, 1)
        applied["noise_sigma"] = sd
    if "zdropout" in enabled and has_c and a.shape[-1] >= 3:
        if r.random() < p["zdropout"]["p_channel"]:
            a = a.copy(); a[..., 2] = 0.0; applied["zdropout_channel"] = 2
    if "proj_window" in enabled and has_c:
        sfrac = p["proj_window"]["shift_frac"]
        g = float(1.0 + r.uniform(-sfrac, sfrac))
        a = a.copy(); a[..., 0] = np.clip(a[..., 0], 0, 1) ** g
        applied["proj_window_gamma"] = round(g, 4)

    return np.clip(a, 0, 1).astype(np.float32), applied


def _zoom_center(a: np.ndarray, s: float) -> np.ndarray:
    """Isotropic zoom about the centre, cropped/padded back to the original size."""
    h, w = a.shape[:2]
    has_c = a.ndim == 3
    zf = (s, s, 1) if has_c else (s, s)
    z = zoom(a, zf, order=1, mode="reflect")
    zh, zw = z.shape[:2]
    out = np.zeros_like(a)
    y0 = (zh - h) // 2; x0 = (zw - w) // 2
    if s >= 1.0:                                    # crop centre
        out[:] = z[y0:y0 + h, x0:x0 + w, ...]
    else:                                           # pad centre
        oy = (h - zh) // 2; ox = (w - zw) // 2
        out[oy:oy + zh, ox:ox + zw, ...] = z
    return out


def augment_with_labels(img: np.ndarray, label: np.ndarray, seed: int,
                        params: dict | None = None):
    """
    Augment an image and its integer label map with the SAME geometric transform,
    so the label (pixel correspondence) is preserved. Photometric ops apply to the
    image only; the label uses nearest-neighbour geometry with matched seed.
    """
    aug_img, applied = augment(img, seed, params)
    # replay only the geometric parts on the label with nearest-neighbour
    r = _rng(seed)
    lab = np.asarray(label).astype(np.float32)
    p = params or default_params()
    # rotate
    deg = float(r.uniform(-p["rotate"]["max_deg"], p["rotate"]["max_deg"]))
    lab = rotate(lab, deg, axes=(0, 1), reshape=False, order=0, mode="constant")
    # reflect
    if r.random() < p["reflect"]["p"]:
        lab = lab[:, ::-1, ...].copy()
    # translate
    h, w = lab.shape[:2]; mf = p["translate"]["max_frac"]
    dy = float(r.uniform(-mf, mf) * h); dx = float(r.uniform(-mf, mf) * w)
    lab = shift(lab, (dy, dx) if lab.ndim == 2 else (dy, dx, 0), order=0, mode="constant")
    # isotropic scale
    s = float(r.uniform(*p["scale_isotropic"]["range"]))
    lab = np.rint(_zoom_center(lab, s)).astype(label.dtype if hasattr(label, "dtype") else int)
    return aug_img, lab, applied
