#!/usr/bin/env python3
"""Build data/animalveg.json.gz — the animal–vegetal analysis (figures 4.3 and 4.4).

Two tabs for the Equatorial Division Plane project.

THE PLANE. Not the equatorial plane through the centre of mass, and not even a plane that halves
the whole cell: the EQUAL-CYTOPLASMIC-VOLUME split. Same orientation as the equatorial plane —
normal along the cell→polar-body axis — but slid along that axis until the two halves hold equal
volume AFTER the pronuclei and polar body are removed. The reference measures those at a median 6%
of the cell, and they sit off-centre, so an uncorrected split charges their volume to one side.

  ANIMAL is the polar-body side. The axis points from the cytoplasm centroid toward the polar
  body, so the positive half is animal by construction and needs no separate convention.

4.4 — THE VOLCANO. Per zygote, per gene: log2(animal density / vegetal density), density being
count ÷ that half's cytoplasm volume. Across zygotes: a one-sample t-test of those per-embryo
values against 0. The test is PAIRED BY CONSTRUCTION — both halves come from the same cell — which
is why a 5-embryo floor is enough to call anything. Floors are ≥20 transcripts summed and ≥5
embryos, and the calling rule is P < 0.05 with |mean log2 FC| ≥ 0.25.

  ⚠️ THOSE P-VALUES ARE UNADJUSTED. At ~380 genes, P < 0.05 alone would return ~19 by chance. The
  BH q is computed and shipped beside them so the page can say so; read the called genes as a
  ranked list, not as hits.

4.3 — THE DENSITY MAP. Where a gene sits inside the cell, pooled over every zygote in one
normalised meridional frame: axial position along the polar axis (−1 vegetal, +1 animal) against
radial distance from that axis (0 centre, 1 cortex). What is shipped, and what the page draws, is
log2(this gene's density / the all-gene density) in that frame, so it is a map of where a gene sits
RELATIVE TO THE REST OF THE PANEL rather than a map of where transcripts are in general — which
would just be a picture of the cell.

  Generalised from the reference deliberately: figure 4.3 draws two fixed gene CLUSTERS (C1
  animal-graded, C4 vegetal-graded) whose membership comes from a table that is not in this repo.
  Shipping the per-gene maps instead lets the same picture be drawn for any gene, and the two
  cluster panels can be reproduced by selecting their members if that table ever arrives.

Output: data/animalveg.json.gz
"""
import collections
import glob
import gzip
import json
import math
import os
import sys

import numpy as np
from scipy import stats

import embryo_stats as ES
# The embryo label is LOOKED UP (data/embryo_ids.json via embryo_naming), never derived and
# never read off a manifest — rebuilding an artifact must not quietly reintroduce a legacy
# name. embryo_label() falls back conspicuously when an embryo is missing from the lookup.
from embryo_naming import embryo_label

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
ASSIGN = os.path.join(DATA, "pronuclei_assignments.json")
ZY_MAN = os.path.join(DATA, "zygote_manifest.json")
PROBESETS = os.path.join(DATA, "probesets.json")
OUT = os.path.join(DATA, "animalveg.json.gz")

VERSION = "animalveg-1.0.0"
MIN_TOTAL = 20            # transcripts summed over embryos for a gene to be tested
MIN_EMBRYOS = 5           # ...in at least this many zygotes
CALL_P = 0.05
CALL_LFC = 0.25
NA, NR = 24, 12           # meridional map: axial bins x radial bins
MIN_MAP_COUNT = 40        # a gene needs this many pooled transcripts to get a map worth drawing


def main():
    man = {m["id"]: m for m in json.load(open(ZY_MAN))["embryos"]}
    probeset = json.load(open(PROBESETS)) if os.path.isfile(PROBESETS) else {}
    polar_of = {k: v.get("polar_body_label") for k, v in man.items()}
    for e in json.load(open(ASSIGN))["embryos"]:
        if e.get("polar"):
            polar_of[e["id"]] = e["polar"]["label"]

    per, emb_meta, skipped = [], [], []
    # pooled maps, in the normalised meridional frame
    gmap = collections.defaultdict(lambda: np.zeros((NA, NR)))
    allmap = np.zeros((NA, NR))

    ids = ES.stage_ids("Zygote")
    for i, eid in enumerate(ids, start=1):
        sc = ES.read_scene(ES.scene_path("Zygote", eid))
        bodies = ES.classify_body(sc)
        if len(bodies) != 1:
            skipped.append({"id": eid, "reason": f"{len(bodies)} body segments"})
            continue
        body = bodies[0]
        cyto = ES.seg_volumes(sc)[body]
        V, F = ES.mesh_of(sc, body)
        com = ES.vol_centroid(V, F)
        pl = polar_of.get(eid)
        if pl is None or str(pl) not in sc["region_meshes"]:
            pl = ES.polar_label(sc)
        if pl is None:
            skipped.append({"id": eid, "reason": "no polar body: no animal-vegetal axis"})
            continue
        pb = ES.mesh_of(sc, pl)[0].mean(0)
        axis = pb - com
        axis /= np.linalg.norm(axis)

        # the equal-CYTOPLASMIC-volume plane: same orientation, slid to a true half
        n_eq, o_eq = ES.equal_volume_plane(V, F, axis, com, exact_total=cyto)
        vAn, vVeg = ES.split_volumes(V, F, n_eq, o_eq, exact_total=cyto)

        TX = ES.cytoplasm_positions(sc, body)
        if not TX:
            skipped.append({"id": eid, "reason": "no cytoplasmic transcripts"})
            continue

        # the normalised meridional frame, from the cytoplasm mesh itself
        d = V - com
        a_ext = float(np.abs(d @ axis).max())
        r_ext = float(np.linalg.norm(d - np.outer(d @ axis, axis), axis=1).max())

        for g, P in TX.items():
            q = P - o_eq
            an = int((q @ n_eq > 0).sum())
            veg = len(P) - an
            per.append({"id": eid, "gene": g, "an": an, "veg": veg, "total": len(P),
                        "vAn": vAn, "vVeg": vVeg,
                        "lfc": math.log2(((an + ES.EPS) / vAn) / ((veg + ES.EPS) / vVeg))})
            dd = P - com
            a = np.clip((dd @ axis) / max(a_ext, 1e-9), -0.999999, 0.999999)
            r = np.clip(np.linalg.norm(dd - np.outer(dd @ axis, axis), axis=1) / max(r_ext, 1e-9),
                        0, 0.999999)
            ia = ((a + 1) / 2 * NA).astype(int)
            ir = (r * NR).astype(int)
            np.add.at(gmap[g], (ia, ir), 1)
            np.add.at(allmap, (ia, ir), 1)

        emb_meta.append({"id": eid, "label": embryo_label(eid),
                         "probeset": probeset.get(eid, "?"),
                         "cyto_vol": round(float(cyto), 1),
                         "vAn": round(float(vAn), 1), "vVeg": round(float(vVeg), 1),
                         "shift_um": round(float((o_eq - com) @ n_eq), 4)})
        print(f"  [{i}/{len(ids)}] {eid:34s} {len(TX):3d} genes  "
              f"split {vAn / cyto:.6f}  shift {(o_eq - com) @ n_eq:+.2f} µm")

    # ---- 4.4: the volcano ----
    by = collections.defaultdict(list)
    for r in per:
        by[r["gene"]].append(r)
    rows = []
    for g, sub in by.items():
        if len(sub) < MIN_EMBRYOS or sum(r["total"] for r in sub) < MIN_TOTAL:
            continue
        y = np.array([r["lfc"] for r in sub], float)
        t, p = stats.ttest_1samp(y, 0.0)
        rows.append({"g": g, "n": len(sub), "total": int(sum(r["total"] for r in sub)),
                     "lfc": round(float(y.mean()), 5), "sd": round(float(y.std(ddof=1)), 5),
                     "p": float(p),
                     "per": [{"id": r["id"], "lfc": round(r["lfc"], 4), "n": r["total"],
                              "an": r["an"], "veg": r["veg"]}
                             for r in sorted(sub, key=lambda r: -abs(r["lfc"]))]})
    if rows:
        for r, q in zip(rows, ES.bh(np.array([r["p"] for r in rows]))):
            r["q"] = float(q)
    rows.sort(key=lambda r: r["p"])
    for i, r in enumerate(rows, start=1):
        r["rank"] = i
        r["called"] = bool(r["p"] < CALL_P and abs(r["lfc"]) >= CALL_LFC)
        r["side"] = "animal" if r["lfc"] > 0 else "vegetal"

    # ---- 4.3: the maps, as log2 vs the all-gene background ----
    all_tot = allmap.sum()
    bg = allmap / max(all_tot, 1)
    maps = {}
    for g, M in gmap.items():
        tot = M.sum()
        if tot < MIN_MAP_COUNT:
            continue
        f = M / tot
        with np.errstate(divide="ignore", invalid="ignore"):
            L = np.log2(np.where((f > 0) & (bg > 0), f / np.where(bg > 0, bg, 1), np.nan))
        maps[g] = {"n": int(tot),
                   "z": [[None if not np.isfinite(v) else round(float(v), 3) for v in row]
                         for row in L]}

    n_called = sum(1 for r in rows if r["called"])
    n_sig = sum(1 for r in rows if r["p"] < CALL_P)
    expected = round(len(rows) * CALL_P, 1)

    doc = {
        "meta": {
            "version": VERSION,
            "method": "figures 4.3 and 4.4 — the specification",
            "n_embryos": len(emb_meta), "n_genes": len(rows), "n_maps": len(maps),
            "skipped": skipped,
            "params": {"MIN_TOTAL": MIN_TOTAL, "MIN_EMBRYOS": MIN_EMBRYOS,
                       "CALL_P": CALL_P, "CALL_LFC": CALL_LFC,
                       "NA": NA, "NR": NR, "MIN_MAP_COUNT": MIN_MAP_COUNT},
            "plane": "the EQUAL-CYTOPLASMIC-VOLUME split: the equatorial orientation slid along "
                     "the polar axis until both halves hold equal volume once the pronuclei and "
                     "polar body are removed",
            "animal": "the polar-body side; the axis points cytoplasm-centroid -> polar body, so "
                      "the positive half is animal by construction",
            "paired": "the t-test is paired by construction — both halves come from the same "
                      "cell — which is why a 5-embryo floor can call anything",
            "unadjusted": f"P is UNADJUSTED: {len(rows)} genes at P < {CALL_P} would give about "
                          f"{expected} by chance and {n_sig} are significant, {n_called} of them "
                          f"also clearing |log2 FC| >= {CALL_LFC}. Read them as a ranked list.",
            "map": "log2(this gene's density / the all-gene density) in a normalised meridional "
                   "frame: axial position on the polar axis (-1 vegetal, +1 animal) against "
                   "radial distance from it (0 centre, 1 cortex), pooled over every zygote",
            "n_called": n_called, "n_significant": n_sig, "expected_by_chance": expected,
        },
        "embryos": emb_meta,
        "genes": rows,
        "maps": maps,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e6:.2f} MB)")
    print(f"  {len(emb_meta)} zygotes, {len(rows)} genes tested, {len(maps)} with a density map")
    print(f"  {n_sig} at P < {CALL_P} ({expected} expected by chance), {n_called} also clear "
          f"|log2 FC| >= {CALL_LFC}")
    if rows:
        print("  strongest: " + ", ".join(
            f"{r['g']} ({r['lfc']:+.2f}, P={r['p']:.1e}, {r['side']})" for r in rows[:5]))


if __name__ == "__main__":
    main()
