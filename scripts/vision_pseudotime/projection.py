"""
Incremental crop + deterministic 2.5D projections (brief items 2, 3).

Never loads a 2.5 GB volume in full. A streaming pass reads ONE z-page at a time
(bounded memory: a handful of H x W planes) and accumulates:

  * MIP            — running per-pixel maximum            (intensity)
  * robust mean    — running nonzero sum / nonzero count  (intensity)
  * occupancy      — running count of label>0 over z       (labels)
  * boundary       — running max of per-slice label edges  (labels)
  * z-profile      — per-page signal, to find the z-extent

From those planes we crop to the embryo, then a SECOND bounded read loads only
the cropped z-range (downsampled) to compute a true per-pixel percentile
projection on the small ROI. The full-resolution volume is never resident.

Deterministic: identical input + parameters -> byte-identical float32 output.

Rule (README): NO raw sum projection — it saturates signal and preserves
background. `sum` is intentionally absent from PROJECTIONS.
"""
from __future__ import annotations

import numpy as np
import tifffile

PROJECTIONS = ("mip", "robust_mean", "percentile", "occupancy", "boundary")
DEFAULT_OUT = 256                       # square edge of the final representation


# ───────────────────────────── streaming accumulation ─────────────────────────────
def _iter_pages(path: str, z_step: int = 1):
    """Yield (z_index, plane) one page at a time. Bounded memory."""
    with tifffile.TiffFile(path) as tf:
        pages = tf.series[0].pages
        for z in range(0, len(pages), z_step):
            yield z, np.asarray(pages[z].asarray())


def scan(path: str, is_labels: bool, z_step: int = 1, label_id: int | None = None) -> dict:
    """One bounded streaming pass. Returns accumulator planes + z profile."""
    mip = None
    nz_sum = nz_cnt = None
    occ = None
    bnd = None
    zprof = []
    n = 0
    for z, plane in _iter_pages(path, z_step):
        n += 1
        if is_labels:
            m = (plane != 0) if label_id is None else (plane == label_id)
            m = m.astype(np.float32)
            occ = m.copy() if occ is None else occ + m
            # per-slice boundary: label pixels adjacent to a different value
            edge = np.zeros_like(m)
            edge[:-1, :] += (plane[:-1, :] != plane[1:, :])
            edge[1:, :] += (plane[:-1, :] != plane[1:, :])
            edge[:, :-1] += (plane[:, :-1] != plane[:, 1:])
            edge[:, 1:] += (plane[:, :-1] != plane[:, 1:])
            edge = ((edge > 0) & (m > 0)).astype(np.float32)
            bnd = edge if bnd is None else np.maximum(bnd, edge)
            zprof.append(float(m.sum()))
        else:
            p = plane.astype(np.float32)
            mip = p.copy() if mip is None else np.maximum(mip, p)
            nzs = p * (p > 0)
            nz_sum = nzs if nz_sum is None else nz_sum + nzs
            c = (p > 0).astype(np.float32)
            nz_cnt = c if nz_cnt is None else nz_cnt + c
            zprof.append(float(p.sum()))
    return {"n_pages_read": n, "z_step": z_step, "is_labels": is_labels,
            "mip": mip, "nz_sum": nz_sum, "nz_cnt": nz_cnt,
            "occupancy": occ, "boundary": bnd, "z_profile": np.asarray(zprof, np.float32)}


# ───────────────────────────── projections ─────────────────────────────
def robust_normalize(a: np.ndarray, lo: float = 1.0, hi: float = 99.5) -> np.ndarray:
    """Percentile-clip to [0,1]. Deterministic; robust to a few hot pixels."""
    a = np.asarray(a, np.float32)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros_like(a, np.float32)
    plo, phi = np.percentile(finite, [lo, hi])
    if phi <= plo:
        phi = plo + 1.0
    return np.clip((a - plo) / (phi - plo), 0.0, 1.0).astype(np.float32)


def mip_projection(acc: dict) -> np.ndarray:
    return robust_normalize(acc["mip"])


def robust_mean_projection(acc: dict) -> np.ndarray:
    m = acc["nz_sum"] / np.maximum(acc["nz_cnt"], 1.0)
    return robust_normalize(m)


def occupancy_projection(acc: dict) -> np.ndarray:
    o = acc["occupancy"]
    return (o / max(acc["n_pages_read"], 1)).astype(np.float32)


def boundary_projection(acc: dict) -> np.ndarray:
    return np.clip(acc["boundary"], 0.0, 1.0).astype(np.float32)


def bbox_from_plane(plane: np.ndarray, thresh: float | None = None, margin_frac: float = 0.06):
    """Tight embryo bounding box from an occupancy / normalized-intensity plane."""
    a = np.asarray(plane, np.float32)
    if thresh is None:
        pos = a[a > 0]
        thresh = float(np.percentile(pos, 25)) if pos.size else 0.0
    mask = a > thresh
    if not mask.any():
        return (0, a.shape[0], 0, a.shape[1])
    ys, xs = np.nonzero(mask)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    my = int((y1 - y0) * margin_frac); mx = int((x1 - x0) * margin_frac)
    return (max(0, y0 - my), min(a.shape[0], y1 + my),
            max(0, x0 - mx), min(a.shape[1], x1 + mx))


def percentile_projection(path: str, bbox, z_lo: int, z_hi: int,
                          q: float = 90.0, max_edge: int = 256) -> np.ndarray:
    """
    True per-pixel percentile ALONG z, computed on the cropped + downsampled ROI
    only. The ROI is small (<= max_edge^2 x z_span), so this is bounded even
    though a per-pixel percentile needs all z for each pixel.
    """
    y0, y1, x0, x1 = bbox
    ds = max(1, int(max(y1 - y0, x1 - x0) / max_edge))
    planes = []
    with tifffile.TiffFile(path) as tf:
        pages = tf.series[0].pages
        for z in range(max(0, z_lo), min(len(pages), z_hi)):
            p = np.asarray(pages[z].asarray())[y0:y1:ds, x0:x1:ds].astype(np.float32)
            planes.append(p)
    if not planes:
        return np.zeros((1, 1), np.float32)
    vol = np.stack(planes, axis=0)                      # (z_roi, h, w) — small
    proj = np.percentile(vol, q, axis=0).astype(np.float32)
    return robust_normalize(proj)


def z_extent(zprof: np.ndarray, frac: float = 0.02):
    """[z_lo, z_hi) covering pages whose signal exceeds frac of the peak."""
    if zprof.size == 0 or zprof.max() <= 0:
        return 0, int(zprof.size)
    keep = np.nonzero(zprof >= frac * zprof.max())[0]
    return int(keep.min()), int(keep.max()) + 1


# ───────────────────────────── resizing / assembly ─────────────────────────────
def _resize_square(a: np.ndarray, out: int) -> np.ndarray:
    """Deterministic anti-aliased resize to (out, out), values kept in [0,1]."""
    from skimage.transform import resize
    r = resize(np.asarray(a, np.float32), (out, out), order=1, mode="edge",
               anti_aliasing=True, preserve_range=True)
    return np.clip(r, 0.0, 1.0).astype(np.float32)


def crop_resize(plane: np.ndarray, bbox, out: int) -> np.ndarray:
    y0, y1, x0, x1 = bbox
    return _resize_square(np.asarray(plane, np.float32)[y0:y1, x0:x1], out)


def single_channel_mip_baseline(intensity_path: str, out: int = DEFAULT_OUT,
                                z_step: int = 1) -> tuple[np.ndarray, dict]:
    """The README 'single normalized MIP' baseline representation."""
    acc = scan(intensity_path, is_labels=False, z_step=z_step)
    mip = mip_projection(acc)
    bbox = bbox_from_plane(mip)
    img = crop_resize(mip, bbox, out)
    return img, {"bbox": bbox, "n_pages_read": acc["n_pages_read"], "z_step": z_step,
                 "channels": ["mip_norm"]}


def build_2p5d(intensity_path: str, label_path: str | None = None, out: int = DEFAULT_OUT,
               z_step: int = 1, label_id: int | None = None) -> tuple[np.ndarray, dict]:
    """
    Deterministic 2.5D representation, cropped to the embryo:
      ch0 = normalized MIP            (intensity)
      ch1 = robust nonzero-mean       (intensity)
      ch2 = segmentation occupancy    (labels)  OR  intensity percentile if no labels

    Also returns provenance (bbox, z-extent, per-channel meaning). No raw sum.
    """
    ai = scan(intensity_path, is_labels=False, z_step=z_step)
    mip = mip_projection(ai)
    rmean = robust_mean_projection(ai)

    meta: dict = {"z_step": z_step, "n_pages_read_intensity": ai["n_pages_read"]}
    if label_path:
        al = scan(label_path, is_labels=True, z_step=z_step, label_id=label_id)
        occ = occupancy_projection(al)
        bnd = boundary_projection(al)
        # crop to where the segmentation says the embryo is
        bbox = bbox_from_plane(occ, thresh=0.0)
        ch2 = occ
        meta["channels"] = ["mip_norm", "robust_mean_norm", "seg_occupancy"]
        meta["n_pages_read_labels"] = al["n_pages_read"]
        meta["boundary_available"] = True
        meta["boundary_plane"] = bnd            # kept for the website/QC, not a model channel
        meta["z_extent_labels"] = z_extent(al["z_profile"])
    else:
        bbox = bbox_from_plane(mip)
        zlo, zhi = z_extent(ai["z_profile"])
        ch2 = percentile_projection(intensity_path, bbox, zlo, zhi, q=90.0, max_edge=out)
        # percentile ROI is already cropped+resized; bring it to full-plane crop frame
        meta["channels"] = ["mip_norm", "robust_mean_norm", "intensity_p90"]
        meta["z_extent_intensity"] = (zlo, zhi)

    c0 = crop_resize(mip, bbox, out)
    c1 = crop_resize(rmean, bbox, out)
    if ch2.shape == (out, out):
        c2 = np.clip(ch2, 0, 1).astype(np.float32)
    else:
        c2 = crop_resize(ch2, bbox, out)
    tensor = np.stack([c0, c1, c2], axis=-1).astype(np.float32)   # (out, out, 3)
    meta["bbox"] = [int(v) for v in bbox]
    meta["out"] = out
    meta["shape"] = list(tensor.shape)
    meta["no_raw_sum"] = True
    return tensor, meta


def save_png(img: np.ndarray, path: str) -> None:
    """Write a preview PNG (uint8) of a [0,1] plane or HxWx3 tensor. Deterministic."""
    from PIL import Image
    a = np.asarray(img, np.float32)
    a = np.clip(a, 0, 1)
    arr = (a * 255 + 0.5).astype(np.uint8)
    Image.fromarray(arr).save(path)
