"""
Movie frame extraction with overlay stripping (brief item 4).

Rendered supplementary movies carry burned-in overlays (timestamps, scale bars,
panel labels, borders). This module extracts frames and records, per frame:
source movie, panel, frame index, visible time, treatment, channel, and
inclusion status — and strips overlays out of the MODEL pixels while keeping any
parsed/manual timestamp ONLY in the metadata table.

Two hard rules from the brief and README:
  * tau is NEVER inferred from frame index. tau stays null unless a manual
    annotation supplies the pronuclear-formation frame AND the NEBD frame.
  * perturbation / overshoot movies never enter normal-development training;
    their inclusion status is 'ood' or 'excluded'.

Decoding uses the system ffmpeg/ffprobe binaries (no Python video dependency).
The pixel-processing functions are pure numpy and are unit-tested on synthetic
frames, so the tests do not depend on the movie files.
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess

import numpy as np

from . import config

FFMPEG = shutil.which("ffmpeg")
FFPROBE = shutil.which("ffprobe")

# default_role (movies.csv) -> per-frame inclusion status for tau work
ROLE_INCLUSION = {
    "single_embryo_pilot": "pilot_normal_dev",       # 1 embryo — feasibility only, not validation
    "split_control_and_ood": "needs_panel_isolation",  # isolate + stage-verify the control panel first
    "ood_only": "ood",
    "exclude_from_tau_training": "excluded",
}


# ───────────────────────────── decode (system ffmpeg) ─────────────────────────────
def probe(path: str) -> dict:
    if not (FFPROBE and os.path.isfile(path)):
        return {"available": False}
    out = subprocess.run(
        [FFPROBE, "-v", "error", "-select_streams", "v:0", "-show_entries",
         "stream=width,height,nb_frames,avg_frame_rate,pix_fmt", "-of", "json", path],
        capture_output=True, text=True)
    try:
        s = json.loads(out.stdout)["streams"][0]
        return {"available": True, "width": int(s["width"]), "height": int(s["height"]),
                "n_frames": int(s.get("nb_frames") or 0), "pix_fmt": s.get("pix_fmt"),
                "fps": s.get("avg_frame_rate")}
    except Exception:                                              # noqa: BLE001
        return {"available": False}


def iter_frames(path: str, w: int, h: int, max_frames: int | None = None):
    """Yield rgb frames (H,W,3 uint8) one at a time from a single ffmpeg pipe."""
    if not FFMPEG:
        raise RuntimeError("ffmpeg binary not found")
    # stderr silenced: when max_frames closes the pipe early, ffmpeg logs a
    # harmless broken-pipe error we do not want in the output.
    proc = subprocess.Popen(
        [FFMPEG, "-v", "error", "-i", path, "-f", "rawvideo", "-pix_fmt", "rgb24", "-"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    fsize = w * h * 3
    i = 0
    try:
        while True:
            if max_frames is not None and i >= max_frames:
                break
            buf = proc.stdout.read(fsize)
            if len(buf) < fsize:
                break
            yield i, np.frombuffer(buf, np.uint8).reshape(h, w, 3)
            i += 1
    finally:
        proc.stdout.close()
        proc.terminate()
        proc.wait()


# ───────────────────────────── pixel processing (pure, testable) ─────────────────────────────
def detect_uniform_border(frame: np.ndarray, tol: int = 4) -> tuple[int, int, int, int]:
    """Crop rows/cols that are (near-)uniform from each edge — kills solid borders."""
    g = frame.astype(np.int16)
    if g.ndim == 3:
        g = g.mean(axis=2)
    h, w = g.shape
    def uniform(line):
        return (line.max() - line.min()) <= tol
    t = 0
    while t < h - 1 and uniform(g[t, :]):
        t += 1
    b = h
    while b > t + 1 and uniform(g[b - 1, :]):
        b -= 1
    l = 0
    while l < w - 1 and uniform(g[:, l]):
        l += 1
    r = w
    while r > l + 1 and uniform(g[:, r - 1]):
        r -= 1
    return t, b, l, r


def strip_overlays(frame: np.ndarray, border_frac: float = 0.0,
                   overlay_boxes: list | None = None,
                   crop_uniform: bool = True) -> tuple[np.ndarray, dict]:
    """
    Remove overlays from MODEL pixels, deterministically:
      1. crop any solid uniform border,
      2. crop a fractional margin (scale bars / labels usually hug the edges),
      3. blank listed overlay boxes (normalized [y0,x0,y1,x1]) e.g. a timestamp corner.
    Returns (clean_frame, record). Records exactly what was removed. This is NOT
    OCR; residual-overlay risk is flagged so nothing pretends the pixels are pristine.
    """
    h0, w0 = frame.shape[:2]
    t, b, l, r = (0, h0, 0, w0)
    if crop_uniform:
        t, b, l, r = detect_uniform_border(frame)
    if border_frac > 0:
        my, mx = int((b - t) * border_frac), int((r - l) * border_frac)
        t, b, l, r = t + my, b - my, l + mx, r - mx
    clean = frame[t:b, l:r].copy()
    ch, cw = clean.shape[:2]
    boxes_applied = []
    for box in (overlay_boxes or []):
        y0, x0, y1, x1 = box
        yy0, yy1 = int(y0 * ch), int(y1 * ch)
        xx0, xx1 = int(x0 * cw), int(x1 * cw)
        clean[yy0:yy1, xx0:xx1] = 0
        boxes_applied.append([yy0, xx0, yy1, xx1])
    return clean, {
        "orig_hw": [h0, w0], "crop_tblr": [t, b, l, r], "clean_hw": [ch, cw],
        "overlay_boxes_zeroed": boxes_applied, "border_frac": border_frac,
        "method": "deterministic border-crop + corner-mask (not OCR)",
        "residual_overlay_risk": "possible — verify manually before publication use",
    }


def split_panels(frame: np.ndarray, rows: int = 1, cols: int = 1) -> list[np.ndarray]:
    """Split a multi-panel comparison frame into an ordered list of sub-panels."""
    h, w = frame.shape[:2]
    ph, pw = h // rows, w // cols
    out = []
    for ri in range(rows):
        for ci in range(cols):
            out.append(frame[ri * ph:(ri + 1) * ph, ci * pw:(ci + 1) * pw])
    return out


def to_gray01(frame: np.ndarray) -> np.ndarray:
    g = np.asarray(frame, np.float32)
    if g.ndim == 3:
        g = g.mean(axis=2)
    m = g.max()
    return (g / m).astype(np.float32) if m > 0 else g.astype(np.float32)


# ───────────────────────────── tau ONLY from annotation ─────────────────────────────
def tau_from_annotation(frame_index: int, pn_formation_frame, nebd_frame):
    """tau = (f - pn)/(nebd - pn) — ONLY with both annotated frames. Else None."""
    if pn_formation_frame is None or nebd_frame is None:
        return None
    if nebd_frame <= pn_formation_frame:
        return None
    tau = (frame_index - pn_formation_frame) / (nebd_frame - pn_formation_frame)
    return float(min(1.0, max(0.0, tau)))


# ───────────────────────────── orchestration ─────────────────────────────
def _movie_rows() -> list[dict]:
    with open(config.MOVIES_CSV, newline="") as fh:
        return [r for r in csv.DictReader(fh) if r.get("movie_id")]


def load_panel_layouts() -> dict:
    """Optional committed movie_panels.json: {movie_id:{rows,cols,overlay_boxes,...}}."""
    p = os.path.join(config.VISION_DIR, "movie_panels.json")
    return json.load(open(p)) if os.path.isfile(p) else {}


def load_annotations() -> dict:
    """Optional committed movie_annotations.json with manual pn/nebd frames + times."""
    p = os.path.join(config.VISION_DIR, "movie_annotations.json")
    return json.load(open(p)) if os.path.isfile(p) else {}


def extract(movie_id: str, out_dir: str | None = None, max_frames: int | None = None,
            save_frames: bool = False) -> dict:
    """
    Extract one movie into a per-frame metadata table with overlays stripped.
    Frame PNGs (if save_frames) go to the gitignored derived dir; the returned
    metadata carries NO pixels and NO absolute paths.
    """
    rows = {r["movie_id"]: r for r in _movie_rows()}
    if movie_id not in rows:
        raise KeyError(f"{movie_id} not in movies.csv")
    m = rows[movie_id]
    role = m.get("default_role", "")
    inclusion = ROLE_INCLUSION.get(role, "excluded")

    # resolve the movie file via the scheffler_movies source dir
    mv_dir = config.resolve_source("scheffler_movies")
    path = os.path.join(mv_dir, m["file_name"]) if mv_dir else ""
    pr = probe(path)

    layouts = load_panel_layouts().get(movie_id, {})
    rows_n, cols_n = int(layouts.get("rows", 1)), int(layouts.get("cols", 1))
    overlay_boxes = layouts.get("overlay_boxes", [])
    ann = load_annotations().get(movie_id, {})
    pn_f, nebd_f = ann.get("pn_formation_frame"), ann.get("nebd_frame")
    visible_times = ann.get("visible_time_by_frame", {})   # manual, per-frame; never OCR

    frames_meta = []
    derived = None
    if save_frames:
        derived = out_dir or os.path.join(config.DERIVED_DIR, f"movie_{movie_id}")
        os.makedirs(derived, exist_ok=True)

    if pr.get("available"):
        for fi, frame in iter_frames(path, pr["width"], pr["height"], max_frames=max_frames):
            clean, rec = strip_overlays(frame, border_frac=float(layouts.get("border_frac", 0.0)),
                                        overlay_boxes=overlay_boxes)
            panels = split_panels(clean, rows_n, cols_n)
            for pi, panel in enumerate(panels):
                tau = tau_from_annotation(fi, pn_f, nebd_f)
                frames_meta.append({
                    "movie_id": movie_id, "panel": pi, "n_panels": len(panels),
                    "frame": fi, "visible_time": visible_times.get(str(fi)),
                    "treatment": m.get("condition_or_subject", ""),
                    "channel": layouts.get("channel"),
                    "inclusion": inclusion,
                    "tau": tau,
                    "tau_source": ("annotation" if tau is not None else
                                   "unavailable — no annotated pn-formation/NEBD"),
                    "clean_hw": rec["clean_hw"],
                })
                if save_frames and derived:
                    from PIL import Image
                    Image.fromarray((to_gray01(panel) * 255 + 0.5).astype(np.uint8)).save(
                        os.path.join(derived, f"f{fi:04d}_p{pi}.png"))
    return {
        "movie_id": movie_id, "file_name": m["file_name"],
        "default_role": role, "inclusion": inclusion,
        "trainable_normal_dev": inclusion == "pilot_normal_dev",
        "probe": {k: pr.get(k) for k in ("available", "width", "height", "n_frames", "fps")},
        "panel_layout": {"rows": rows_n, "cols": cols_n},
        "n_frame_panel_rows": len(frames_meta),
        "any_tau_available": any(fr["tau"] is not None for fr in frames_meta),
        "frames": frames_meta,
    }


if __name__ == "__main__":
    import sys
    mid = sys.argv[1] if len(sys.argv) > 1 else "SM01"
    r = extract(mid, max_frames=int(sys.argv[2]) if len(sys.argv) > 2 else None)
    print(f"{r['movie_id']} ({r['file_name']}) role={r['default_role']} inclusion={r['inclusion']}")
    print(f"  probe={r['probe']}  panel_rows={r['n_frame_panel_rows']}  "
          f"any_tau={r['any_tau_available']}")
    if r["frames"]:
        print(f"  frame0 meta: {r['frames'][0]}")
