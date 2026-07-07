#!/usr/bin/env python3
"""
Build the viewer's per-embryo scene data for the 45 sperm-positive embryos.

Source of the 3-D geometry: the MERFISH atlas already precomputes, for every
embryo, a `scene.json.gz` (region meshes via marching-cubes on the segmentation
label, per-gene transcript point clouds, extents, colors) using the pipeline in
  /Users/rishib/Desktop/MERFISH/Website2/MerfishAtlasWebsite
We reuse that geometry VERBATIM (identical rendering style) but:
  * keep only the 45 embryos in data/sperm_transcriptomics/embryos.csv,
  * slim each scene to just what this minimal viewer needs
    (meshes + selected-gene clouds + sperm), and
  * overwrite the sperm marker with THIS project's authoritative coordinates
    (data/sperm_transcriptomics/embryos.csv), which supersede the atlas's.

Outputs:
  public/data/scenes/<embryo_id>.json.gz   slim per-embryo scene
  public/data/manifest.json                the 45-embryo index for the nav bar
"""
import csv
import gzip
import json
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
EMB_CSV = os.path.join(REPO, "data", "sperm_transcriptomics", "embryos.csv")
ATLAS = "/Users/rishib/Desktop/MERFISH/Website2/MerfishAtlasWebsite/public/data"
OUT_SCENES = os.path.join(HERE, "public", "data", "scenes")
OUT_MANIFEST = os.path.join(HERE, "public", "data", "manifest.json")

# canonical developmental order for grouping the nav bar
STAGE_ORDER = ["Oocyte", "Zygote", "2-cell (early)", "2-cell (late)"]

SPERM_COLOR = "#ff2d95"
SPERM_SIZE = 11

XY_UM = 0.15   # µm per pixel (xy)
Z_UM = 1.0     # µm per z-frame


def _unit(v):
    n = float(np.linalg.norm(v))
    return [round(float(c), 5) for c in (v / n)] if n > 0 else [0.0, 0.0, 0.0]


def compute_analysis(d, sperm, zs):
    """Pre-compute the three per-gene display vectors, in physical µm space.

    Points live in the viewer's plot space (x,y = px, plot-z = frame*z_scale).
    We convert to µm (x*0.15, y*0.15, frame*1.0) so directions are physically
    isotropic, compute unit vectors there, and keep anchors in plot space for
    rendering.

      1) PCA  : first principal-component axis of the gene's point cloud.
      2) g2e  : unit vector from the gene cloud's centre of mass to the embryo
                centre of mass (all segmentation voxels, uniform density).
      3) sperm_to_emb : unit vector from the sperm to the embryo centre of mass
                (gene-independent).
    """
    dc = d.get("data_centroid")
    if not dc:
        return None
    emb_plot = [round(dc["x"], 2), round(dc["y"], 2), round(dc["z"], 2)]
    emb_um = np.array([dc["x"] * XY_UM, dc["y"] * XY_UM, dc["z"] / zs * Z_UM])

    sp_plot = [round(sperm["x"], 2), round(sperm["y"], 2), round(sperm["z"] * zs, 2)]
    sp_um = np.array([sperm["x"] * XY_UM, sperm["y"] * XY_UM, sperm["z"] * Z_UM])

    genes = {}
    for g, t in d.get("transcripts", {}).items():
        x = np.asarray(t["x"], float)
        y = np.asarray(t["y"], float)
        z = np.asarray(t["gz"], float)          # z-frame
        n = len(x)
        if n == 0:
            continue
        com_plot = [round(float(x.mean()), 2), round(float(y.mean()), 2),
                    round(float(z.mean() * zs), 2)]
        com_um = np.array([x.mean() * XY_UM, y.mean() * XY_UM, z.mean() * Z_UM])
        g2e = _unit(emb_um - com_um)
        pca = None      # PC1 unit axis
        evr = None      # explained-variance ratio [PC1, PC2, PC3] (sums to 1)
        if n >= 2:
            P = np.stack([x * XY_UM, y * XY_UM, z * Z_UM], axis=1)
            P = P - P.mean(axis=0)
            try:
                _, S, Vt = np.linalg.svd(P, full_matrices=False)
                v = Vt[0]
                if v[int(np.argmax(np.abs(v)))] < 0:      # deterministic sign
                    v = -v
                pca = [round(float(c), 5) for c in v]
                var = np.asarray(S, float) ** 2
                var = np.concatenate([var, np.zeros(3 - len(var))])[:3]
                tot = float(var.sum())
                if tot > 0:
                    evr = [round(float(e), 5) for e in (var / tot)]
            except np.linalg.LinAlgError:
                pca = None
        genes[g] = {"com": com_plot, "pca": pca, "g2e": g2e, "evr": evr, "n": n}

    ex = d["extents"]
    ranges = [ex["x"][1] - ex["x"][0], ex["y"][1] - ex["y"][0], ex["z"][1] - ex["z"][0]]
    arrow_len = round(0.4 * (sum(ranges) / 3.0), 1)

    return {
        "embryo_com": emb_plot,
        "sperm_plot": sp_plot,
        "sperm_to_emb": _unit(emb_um - sp_um),
        "arrow_len": arrow_len,
        "genes": genes,
    }


def find_atlas_scene(embryo_id):
    for cat in os.listdir(ATLAS):
        p = os.path.join(ATLAS, cat, embryo_id, "scene.json.gz")
        if os.path.isfile(p):
            return p
    return None


def main():
    os.makedirs(OUT_SCENES, exist_ok=True)
    rows = list(csv.DictReader(open(EMB_CSV, newline="")))
    manifest = []
    index = []          # compact cross-embryo vectors for the violin plots
    patched = 0

    for r in rows:
        eid = r["embryo_id"]
        src = find_atlas_scene(eid)
        if not src:
            print(f"  !! no atlas scene for {eid} -- skipped")
            continue
        d = json.load(gzip.open(src, "rt"))

        # ---- authoritative sperm coordinate (this project) ----
        sperm = {
            "x": float(r["sperm_x_px"]),
            "y": float(r["sperm_y_px"]),
            "z": float(r["sperm_z_frame"]),
            "segment": r["sperm_segment"],
        }
        atlas_mk = (d.get("initial", {}).get("markers") or [{}])[0]
        if (abs(atlas_mk.get("x", -1) - sperm["x"]) > 1 or
                abs(atlas_mk.get("y", -1) - sperm["y"]) > 1 or
                abs(atlas_mk.get("z", -1) - sperm["z"]) > 1):
            patched += 1

        # ---- genes sorted by total abundance (default = most abundant) ----
        totals = d.get("gene_totals", {})
        genes = sorted(d.get("genes", []),
                       key=lambda g: (-totals.get(g, 0), g))

        # ---- slim transcript clouds: keep x, y, gz per gene (drop mask) ----
        tx = {}
        for g, t in d.get("transcripts", {}).items():
            tx[g] = {"x": t["x"], "y": t["y"], "gz": t["gz"]}

        title = f"{r['stage']} · {r['plate'].upper()} · FOV " + \
                "".join(ch for ch in r["fov"] if ch.isdigit())

        analysis = compute_analysis(d, sperm, d.get("z_scale", 7.0))
        scene = {
            "id": eid,
            "stage": r["stage"],
            "title": title,
            "merfish_index": r["merfish_index"],
            "z_scale": d.get("z_scale", 7.0),
            "extents": d["extents"],
            "region_meshes": d["region_meshes"],
            "region_defaults": d["region_defaults"],
            "mask_labels": d["mask_labels"],
            "genes": genes,
            "gene_totals": totals,
            "transcripts": tx,
            "sperm": sperm,
            "sperm_color": SPERM_COLOR,
            "sperm_size": SPERM_SIZE,
            "analysis": analysis,
        }
        out = os.path.join(OUT_SCENES, eid + ".json.gz")
        with gzip.open(out, "wt") as fh:
            json.dump(scene, fh, separators=(",", ":"))

        # compact, unique-within-stage tab label, e.g. "P1 · 0_2" / "P2 · 18"
        m2 = re.search(r"_p\d+_(.+)$", eid)
        fovsub = m2.group(1) if m2 else "".join(c for c in r["fov"] if c.isdigit())
        label = f"{r['plate'].upper()} · {fovsub}"
        # imaging date (embryos of the same plate/fov but different runs exist)
        ymd = eid[:8]
        MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        date_short = f"{MON[int(ymd[4:6]) - 1]} {int(ymd[6:8])}" if ymd.isdigit() else ""

        manifest.append({
            "id": eid,
            "stage": r["stage"],
            "plate": r["plate"],
            "fov": r["fov"],
            "label": label,
            "date_short": date_short,
            "title": title,
            "merfish_index": int(r["merfish_index"]),
            "n_transcripts": int(r["n_transcripts"]),
            "n_genes": int(r["n_genes"]),
            "default_gene": genes[0] if genes else None,
            "scene": f"data/scenes/{eid}.json.gz",
            "size_kb": round(os.path.getsize(out) / 1024),
        })

        # cross-embryo index: per-embryo sperm→COM axis + per-gene PC1 / gene→COM
        # unit vectors, so the front-end can compute |cos| dot products across all
        # embryos for the selected gene (violin plots).
        if analysis:
            index.append({
                "id": eid, "stage": r["stage"], "label": label,
                "sperm_to_emb": analysis["sperm_to_emb"],
                "genes": {g: {"pca": gv["pca"], "g2e": gv["g2e"], "n": gv["n"]}
                          for g, gv in analysis["genes"].items()},
            })

    # order the nav bar by developmental stage, then plate, then FOV number
    def sort_key(m):
        st = STAGE_ORDER.index(m["stage"]) if m["stage"] in STAGE_ORDER else 99
        fovn = int("".join(ch for ch in m["fov"] if ch.isdigit()) or 0)
        return (st, m["plate"], fovn, m["id"])
    manifest.sort(key=sort_key)

    with open(OUT_MANIFEST, "w") as fh:
        json.dump({"embryos": manifest, "stage_order": STAGE_ORDER}, fh, indent=1)

    idx_path = os.path.join(HERE, "public", "data", "analysis_index.json.gz")
    with gzip.open(idx_path, "wt") as fh:
        json.dump({"embryos": index}, fh, separators=(",", ":"))
    print(f"wrote analysis_index.json.gz  ({os.path.getsize(idx_path)/1024:.0f} KB)")

    total_kb = sum(m["size_kb"] for m in manifest)
    print(f"wrote {len(manifest)} scenes  ({total_kb/1024:.1f} MB total)")
    print(f"sperm markers patched to this project's coords: {patched}")
    biggest = max(manifest, key=lambda m: m["size_kb"])
    print(f"largest scene: {biggest['id']}  {biggest['size_kb']} KB "
          f"({biggest['n_transcripts']:,} transcripts)")


if __name__ == "__main__":
    main()
