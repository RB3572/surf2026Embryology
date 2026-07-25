"""
Evaluation: grouped metrics, trivial baselines, saliency (brief item 9).

Everything here is embryo-GROUPED where grouping matters, so a metric can never
be inflated by adjacent frames of one embryo landing on both sides. Trivial
baselines (constant, mean-brightness, corner-overlay) are provided so a model
must be shown to BEAT them; a model that only matches the corner-overlay baseline
is reading the timestamp, not the embryo. Saliency (occlusion) checks that
prediction mass sits on the embryo, not on overlays/borders.

No metric here fabricates a value; each is computed from arrays the caller
supplies.
"""
from __future__ import annotations

import numpy as np


# ───────────────────────────── pseudotime metrics ─────────────────────────────
def embryo_grouped_mae(y_true, y_pred, groups) -> dict:
    y_true, y_pred = np.asarray(y_true, float), np.asarray(y_pred, float)
    groups = np.asarray(groups)
    per = {}
    for g in np.unique(groups):
        m = groups == g
        per[str(g)] = float(np.mean(np.abs(y_true[m] - y_pred[m])))
    macro = float(np.mean(list(per.values()))) if per else float("nan")
    overall = float(np.mean(np.abs(y_true - y_pred))) if len(y_true) else float("nan")
    return {"macro_mae": macro, "overall_mae": overall, "per_group_mae": per,
            "n_groups": len(per)}


def spearman(y_true, y_pred) -> float:
    from scipy.stats import spearmanr
    if len(y_true) < 3:
        return float("nan")
    r = spearmanr(y_true, y_pred).correlation
    return float(r) if r == r else float("nan")


def pair_order_accuracy(y_true, y_pred, tol: float = 0.0) -> dict:
    """Explicit concordant/discordant pair counting over comparable pairs."""
    yt, yp = np.asarray(y_true, float), np.asarray(y_pred, float)
    n = len(yt)
    conc = disc = tied_pred = comparable = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs(yt[i] - yt[j]) <= tol:
                continue                                # not comparable on truth
            comparable += 1
            dt = yt[i] - yt[j]; dp = yp[i] - yp[j]
            if dp == 0:
                tied_pred += 1
            elif (dt > 0) == (dp > 0):
                conc += 1
            else:
                disc += 1
    acc = conc / comparable if comparable else float("nan")
    return {"comparable": comparable, "concordant": conc, "discordant": disc,
            "tied_pred": tied_pred, "strict_accuracy": acc}


def interval_coverage(y_true, lo, hi) -> float:
    yt, lo, hi = np.asarray(y_true, float), np.asarray(lo, float), np.asarray(hi, float)
    if not len(yt):
        return float("nan")
    return float(np.mean((yt >= lo) & (yt <= hi)))


# ───────────────────────────── segmentation metrics ─────────────────────────────
def dice(a, b) -> float:
    a, b = np.asarray(a).astype(bool), np.asarray(b).astype(bool)
    s = a.sum() + b.sum()
    return float(2 * (a & b).sum() / s) if s else 1.0


def iou(a, b) -> float:
    a, b = np.asarray(a).astype(bool), np.asarray(b).astype(bool)
    u = (a | b).sum()
    return float((a & b).sum() / u) if u else 1.0


def centroid_error(mask_pred, mask_true, voxel_um=(1.0, 1.0, 1.0)) -> float:
    """Euclidean distance between mask centroids in physical units."""
    def cen(m):
        idx = np.nonzero(np.asarray(m).astype(bool))
        if len(idx[0]) == 0:
            return None
        # idx order matches array dims; scale by matching voxel sizes (reversed to x,y,z if 3D)
        coords = [ax.mean() for ax in idx]
        return np.array(coords, float)
    cp, ct = cen(mask_pred), cen(mask_true)
    if cp is None or ct is None:
        return float("nan")
    vox = np.asarray(voxel_um[:len(cp)], float)
    return float(np.linalg.norm((cp - ct) * vox))


def ood_rejection_rate(is_ood, flagged) -> dict:
    """How well an OOD flag catches known-OOD samples (and spares in-domain ones)."""
    is_ood = np.asarray(is_ood).astype(bool)
    flagged = np.asarray(flagged).astype(bool)
    tp = int((is_ood & flagged).sum()); fn = int((is_ood & ~flagged).sum())
    fp = int((~is_ood & flagged).sum()); tn = int((~is_ood & ~flagged).sum())
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    specificity = tn / (tn + fp) if (tn + fp) else float("nan")
    return {"ood_recall": recall, "in_domain_specificity": specificity,
            "tp": tp, "fn": fn, "fp": fp, "tn": tn}


# ───────────────────────────── trivial / artifact baselines ─────────────────────────────
def constant_baseline(y_true_train) -> float:
    return float(np.mean(y_true_train)) if len(y_true_train) else 0.5


def brightness_baseline(images, y_true_train, images_train):
    """Rank test images by mean brightness fit to train tau (monotone). A model
    that only matches this is using global exposure, not embryo structure."""
    bt = np.array([np.asarray(im, float).mean() for im in images_train])
    yt = np.asarray(y_true_train, float)
    a, b = np.polyfit(bt, yt, 1) if len(bt) >= 2 else (0.0, float(np.mean(yt)))
    return np.clip(a * np.array([np.asarray(im, float).mean() for im in images]) + b, 0, 1)


def corner_overlay_baseline(images, y_true_train, images_train, frac=0.18):
    """Predict tau from a CORNER patch only (where timestamps live). If a model
    ties this, it is reading the overlay. Fit corner-brightness -> tau."""
    def corner(im):
        a = np.asarray(im, float)
        h, w = a.shape[:2]
        ch, cw = int(h * frac), int(w * frac)
        return a[:ch, :cw].mean()
    ct = np.array([corner(im) for im in images_train]); yt = np.asarray(y_true_train, float)
    a, b = np.polyfit(ct, yt, 1) if len(ct) >= 2 and ct.std() > 0 else (0.0, float(np.mean(yt)))
    return np.clip(a * np.array([corner(im) for im in images]) + b, 0, 1)


# ───────────────────────────── saliency (occlusion) ─────────────────────────────
def occlusion_saliency(predict_fn, image, patch: int = 24, stride: int = 12,
                       fill: float = 0.0) -> np.ndarray:
    """
    Map |Δprediction| when each patch is occluded. High values = the model
    relies on that region. Deterministic.

    `fill` defaults to 0.0, which is correct for dark-background fluorescence:
    occluding background (already ~0) changes nothing, so saliency concentrates
    on real structure. Filling with the image mean would inject signal into the
    background and spuriously light up empty regions.
    """
    a = np.asarray(image, np.float32)
    h, w = a.shape[:2]
    base = float(predict_fn(a))
    sal = np.zeros((h, w), np.float32)
    cnt = np.zeros((h, w), np.float32)
    for y in range(0, h, stride):
        for x in range(0, w, stride):
            occ = a.copy()
            occ[y:y + patch, x:x + patch, ...] = fill
            d = abs(float(predict_fn(occ)) - base)
            sal[y:y + patch, x:x + patch] += d
            cnt[y:y + patch, x:x + patch] += 1
    return sal / np.maximum(cnt, 1)


def saliency_on_embryo_fraction(saliency, embryo_mask) -> float:
    """Fraction of total saliency mass falling inside the embryo mask.
    Near 1 = attends to the embryo; low = attends to borders/overlays."""
    s = np.asarray(saliency, float)
    m = np.asarray(embryo_mask).astype(bool)
    tot = s.sum()
    return float(s[m].sum() / tot) if tot > 0 else float("nan")


def saliency_contrast(saliency, embryo_mask) -> float:
    """Mean saliency DENSITY inside the embryo / density outside it.
    >1 means the model relies on the embryo more than the background per pixel."""
    s = np.asarray(saliency, float)
    m = np.asarray(embryo_mask).astype(bool)
    inside = s[m].mean() if m.any() else 0.0
    outside = s[~m].mean() if (~m).any() else 0.0
    return float(inside / outside) if outside > 0 else float("inf")


def saliency_box_mass(saliency, box) -> float:
    """Fraction of saliency mass in a box [y0,x0,y1,x1] — e.g. a timestamp corner.
    Near 0 means the model ignores that overlay region."""
    s = np.asarray(saliency, float)
    y0, x0, y1, x1 = box
    tot = s.sum()
    return float(s[y0:y1, x0:x1].sum() / tot) if tot > 0 else float("nan")
