#!/usr/bin/env python3
"""Build data/sperm_locations.json (+ data/gfp/*.jpg) — find embryos by WHERE the sperm is.

A browser over the sperm-positive embryos, sorted by how the sperm sits relative to the
structures that might matter:

  zygote   distance to the cortex · to the polar body · to the maternal pronucleus ·
           to the paternal pronucleus
  2-cell   distance to the JUNCTION between the two blastomeres (the one the PI cares about) ·
           to the polar body · to the cortex, plus which blastomere the sperm ended up in

Geometry comes from data/axes/<id>.json.gz, whose landmarks already carry the manually-labelled
sperm (`sperm_plot`), the polar body, the nuclei and the body centroids. Maternal/paternal
identity for zygotes is joined from pronuclei_assignments.json (consensus.female indexes which
of the two pronuclei is the female one).

THE JUNCTION
------------
A 2-cell axes scene ships TWO body centroids. The junction is taken as the plane through their
midpoint with normal û = unit(B − A), and the sperm's distance to it is |(s − M)·û| — signed, so
the sign also says which blastomere the sperm is in. This is the same construction the contact
project uses for the interface, kept deliberately identical so the two projects agree.

GFP STILLS
----------
The site is static and cannot reach /Volumes/HW, so the GFP frames are rendered here into
data/gfp/ as PNGs: for every embryo whose sperm was actually pointed at in the GFP image, the
z-slice the labeller marked, in BOTH channels, plus a max-Z projection of each channel. The
labelled sperm pixel travels in the JSON so the page can draw the marker as an overlay rather
than burning it into the image.

Units: landmarks are in PLOT space (x,y px; z × z_scale) — everything is converted to µm here
(x,y × 0.15; z ÷ z_scale) so every distance printed on the page is a real micron.
"""
import csv
import glob
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
AXES = os.path.join(DATA, "axes")
MANIFEST = os.path.join(DATA, "axes_manifest.json")
ASSIGN = os.path.join(DATA, "pronuclei_assignments.json")
SPERM_CSV = os.path.join(HERE, "..", "data", "merfish_sperm.csv")
GFP_DIR = os.path.join(DATA, "gfp")
OUT = os.path.join(DATA, "sperm_locations.json")

XY_UM = 0.15
GFP_W = 1600                # rendered width; source frames are 2304². JPEG, not PNG: these are
GFP_JPEG_Q = 88             # single-channel fluorescence, and PNG would cost ~185 MB vs ~35 MB
GFP_LO, GFP_HI = 50.0, 99.7   # the same percentile window the labelling tool uses…
GFP_GAMMA = 2.24              # …and the same gamma, so these look like what was labelled
CH_COLOR = {0: (0, 1, 0), 1: (0.25, 0.45, 1.0)}   # 488 → green, 405 → blue


def um(p, zs):
    """plot → µm."""
    return np.array([p[0] * XY_UM, p[1] * XY_UM, p[2] / zs], float)


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n else v


def mesh_verts_um(mesh, zs):
    v = np.asarray(mesh["verts"], float).reshape(-1, 3)
    return np.stack([v[:, 0] * XY_UM, v[:, 1] * XY_UM, v[:, 2] / zs], axis=1)


# ───────────────────────── geometry ─────────────────────────
def measure(eid, stage, assign):
    p = os.path.join(AXES, eid + ".json.gz")
    if not os.path.isfile(p):
        return None
    s = json.load(gzip.open(p, "rt"))
    lm = s.get("landmarks") or {}
    sp = lm.get("sperm_plot")
    if not sp:
        return None
    zs = s.get("z_scale", 6.667)
    sperm = um(sp, zs)

    out = {"id": eid, "label": s.get("stage_label") or stage, "stage": stage,
           "sperm_plot": [round(float(v), 2) for v in sp],
           "z_scale": zs, "metrics": {}}
    M = out["metrics"]

    # distance to the cortex = nearest point on the cell-body surface. For a 2-cell that is the
    # surface of whichever blastomere the sperm is in, so the two stages stay comparable.
    meshes = s.get("region_meshes") or {}
    bodies = [k for k in meshes if int(k) in (1, 2)] if stage == "twocell" else \
             [k for k in meshes if int(k) == 1]
    best = None
    for k in bodies:
        V = mesh_verts_um(meshes[k], zs)
        d = float(np.min(np.linalg.norm(V - sperm, axis=1)))
        if best is None or d < best:
            best = d
    if best is not None:
        M["cortex"] = round(best, 2)

    polar = lm.get("polar_plot")
    if polar:
        M["polar"] = round(float(np.linalg.norm(um(polar, zs) - sperm)), 2)

    if stage == "twocell":
        bp = lm.get("body_plots") or []
        if len(bp) >= 2:
            a, b = um(bp[0], zs), um(bp[1], zs)
            axis = unit(b - a)
            mid = (a + b) / 2.0
            signed = float((sperm - mid) @ axis)
            M["junction"] = round(abs(signed), 2)
            out["junction"] = {
                "mid_um": [round(float(v), 3) for v in mid],
                "axis_um": [round(float(v), 5) for v in axis],
                "sep_um": round(float(np.linalg.norm(b - a)), 2),
                "side": 1 if signed > 0 else 0,        # which blastomere the sperm sits in
                "signed_um": round(signed, 2),
            }
        nuc = lm.get("nuclei_plots") or []
        if nuc:
            M["nucleus"] = round(min(float(np.linalg.norm(um(n, zs) - sperm)) for n in nuc), 2)
    else:
        a = assign.get(eid)
        if a and a.get("pron") and len(a["pron"]) == 2:
            fem = (a.get("consensus") or {}).get("female")
            if fem in (0, 1):
                mat = um(a["pron"][fem]["com_plot"], zs)
                pat = um(a["pron"][1 - fem]["com_plot"], zs)
                M["maternal"] = round(float(np.linalg.norm(mat - sperm)), 2)
                M["paternal"] = round(float(np.linalg.norm(pat - sperm)), 2)
                out["pron_um"] = {"maternal": [round(float(v), 3) for v in mat],
                                  "paternal": [round(float(v), 3) for v in pat]}
    return out


# ───────────────────────── GFP stills ─────────────────────────
def stretch(a):
    a = a.astype(np.float32)
    lo, hi = np.percentile(a, GFP_LO), np.percentile(a, GFP_HI)
    return np.clip((a - lo) / max(1.0, hi - lo), 0, 1) ** GFP_GAMMA


def save_png(arr01, ch, path):
    from PIL import Image
    r, g, b = CH_COLOR.get(ch, (1, 1, 1))
    H, W = arr01.shape
    rgb = np.zeros((H, W, 3), np.uint8)
    rgb[..., 0] = (arr01 * 255 * r).astype(np.uint8)
    rgb[..., 1] = (arr01 * 255 * g).astype(np.uint8)
    rgb[..., 2] = (arr01 * 255 * b).astype(np.uint8)
    Image.fromarray(rgb).resize((GFP_W, int(H * GFP_W / W))).save(
        path, "JPEG", quality=GFP_JPEG_Q, optimize=True, subsampling=0)


def render_gfp(rows_by_id):
    """One z-slice per channel at the labelled z, plus a max-Z projection per channel."""
    try:
        import tifffile
    except ImportError:
        print("  !! tifffile missing — skipping GFP stills"); return {}
    os.makedirs(GFP_DIR, exist_ok=True)
    out = {}
    for eid, r in sorted(rows_by_id.items()):
        path = (r.get("associated_gfp_path") or "").strip()
        zf = str(r.get("gfp_z_frame") or "").strip()
        if not path or not zf or not os.path.isfile(path):
            continue
        try:
            z = int(float(zf))
            m = tifffile.memmap(path, mode="r")
            H, W = m.shape[-2], m.shape[-1]
            fl = m.reshape(-1, H, W)
            nz = fl.shape[0] // 2
            files = {}
            for ch in (0, 1):
                k = min(fl.shape[0] - 1, max(0, 2 * (z - 1) + ch))
                fn = f"{eid}_ch{ch}.jpg"
                save_png(stretch(np.asarray(fl[k])), ch, os.path.join(GFP_DIR, fn))
                files[f"ch{ch}"] = fn
                # max-Z projection over the whole stack, sampled so a deep stack stays quick
                step = max(1, nz // 60)
                acc = None
                for zz in range(0, nz, step):
                    idx = 2 * zz + ch
                    if idx >= fl.shape[0]:
                        break
                    fr = np.asarray(fl[idx])
                    acc = fr if acc is None else np.maximum(acc, fr)
                if acc is not None:
                    fn = f"{eid}_mip{ch}.jpg"
                    save_png(stretch(acc), ch, os.path.join(GFP_DIR, fn))
                    files[f"mip{ch}"] = fn
            out[eid] = {
                "files": files, "z": z, "nz": int(nz),
                "src_w": int(W), "src_h": int(H), "w": GFP_W,
                "x": int(float(r["gfp_x_px"])) if str(r.get("gfp_x_px") or "").strip() else None,
                "y": int(float(r["gfp_y_px"])) if str(r.get("gfp_y_px") or "").strip() else None,
                "name": os.path.basename(path),
            }
            print(f"    gfp {eid}: z {z}/{nz}  ({len(files)} images)")
        except Exception as e:                       # noqa: BLE001
            print(f"    !! gfp {eid}: {e}")
    return out


def main():
    man = json.load(open(MANIFEST))["embryos"]
    assign = {e["id"]: e for e in json.load(open(ASSIGN))["embryos"]}
    print(f"sperm locations: {len(man)} embryos in the axes cohort\n")

    embryos = []
    for e in man:
        r = measure(e["id"], e.get("stage"), assign)
        if r:
            r["label"] = e.get("label") or r["label"]
            r["stage_label"] = e.get("stage_label") or ""
            embryos.append(r)

    # early vs late 2-cell, off the id (the axes manifest only says "twocell")
    for r in embryos:
        if r["stage"] == "twocell":
            low = r["id"].lower()
            r["sub"] = "e2c" if ("e2c" in low or "early" in low) else \
                       ("l2c" if ("l2c" in low or "late" in low) else "2c")
        else:
            r["sub"] = "zygote"

    n_z = sum(1 for r in embryos if r["stage"] == "zygote")
    n_e = sum(1 for r in embryos if r["sub"] == "e2c")
    n_l = sum(1 for r in embryos if r["sub"] == "l2c")
    print(f"  measured: {len(embryos)}  ({n_z} zygote, {n_e} early-2C, {n_l} late-2C)")
    print(f"  with a junction distance: {sum(1 for r in embryos if 'junction' in r['metrics'])}")
    print(f"  with maternal/paternal:   {sum(1 for r in embryos if 'maternal' in r['metrics'])}\n")

    # Only render stills for embryos that actually appear in the browser. Some labelled embryos
    # have no axes scene (no usable segmentation), and rendering those would ship megabytes of
    # images nothing can ever link to — plus a dangling reference in the artifact.
    known = {e["id"] for e in embryos}
    rows, skipped_gfp = {}, 0
    if os.path.isfile(SPERM_CSV):
        for r in csv.DictReader(open(SPERM_CSV)):
            eid = (r.get("resolved_embryo_id") or "").strip()
            if eid and str(r.get("gfp_z_frame") or "").strip():
                if eid in known:
                    rows[eid] = r
                else:
                    skipped_gfp += 1
    if skipped_gfp:
        print(f"  ({skipped_gfp} GFP-labelled embryos have no axes scene — not rendered)")
    print(f"  rendering GFP stills for {len(rows)} labelled embryos…")
    gfp = render_gfp(rows)

    payload = {
        "meta": {
            "version": "sperm-locations-1.0.0",
            "n": len(embryos), "n_zygote": n_z, "n_e2c": n_e, "n_l2c": n_l,
            "n_gfp": len(gfp), "gfp_dir": "data/gfp",
            "xy_um": XY_UM,
            "sorts": {
                "zygote": [
                    {"key": "cortex", "label": "Distance to the cortex"},
                    {"key": "polar", "label": "Distance to the polar body"},
                    {"key": "maternal", "label": "Distance to the maternal pronucleus"},
                    {"key": "paternal", "label": "Distance to the paternal pronucleus"},
                ],
                "twocell": [
                    {"key": "junction", "label": "Distance to the blastomere junction"},
                    {"key": "polar", "label": "Distance to the polar body"},
                    {"key": "cortex", "label": "Distance to the cortex"},
                    {"key": "nucleus", "label": "Distance to the nearest nucleus"},
                ],
            },
        },
        "gfp": gfp,
        "embryos": embryos,
    }
    with open(OUT, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    total = sum(os.path.getsize(os.path.join(GFP_DIR, f))
                for f in os.listdir(GFP_DIR)) / 1e6 if os.path.isdir(GFP_DIR) else 0
    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1024:.0f} KB)"
          f" + {len(gfp)} embryos of GFP stills ({total:.1f} MB)")


if __name__ == "__main__":
    main()
