#!/usr/bin/env python3
"""Build data/contacthalves.json.gz — the contact region, on the reference's own definition
(figures 7.1, 7.2 and 7.3).

WHY THIS EXISTS ALONGSIDE build_contact.py
------------------------------------------
The site already has a contact project. It asks the same question with a different instrument: a
SLAB of half-thickness D straddling the junction, scored against where transcripts generally are
in that embryo. That is a good test and it needs no volume estimate — but it is not what the
reference does, and the two are not interchangeable. A slab is one region shared by both
blastomeres and its thickness is a free parameter; the reference splits EACH BLASTOMERE IN TWO
and lets geometry fix the boundary:

  ⓵ THE AXIS is the junction→edge direction, per blastomere: û = unit(COM_self − COM_other).
  ⓶ THE PLANE is perpendicular to û and slid along it until THIS BLASTOMERE'S TWO HALVES HOLD
    EQUAL VOLUME. Not the centroid plane. Blastomeres are flattened against each other at the
    junction and rounded at the edge, so the equal-volume plane sits off the centroid — and any
    volume imbalance would otherwise appear as a gene-independent contact bias in every gene.
  ⓷ CONTACT is the junction-side half, EDGE the outer half, and there are two of each per embryo.

  Nuclei and the polar body are excluded by segment label, not by geometry: a nucleus sitting near
  the junction would otherwise read as contact enrichment.

7.1 — THE PROFILE. Before asking about any gene, what do ALL transcripts do along that axis?
Position is normalised per blastomere by that cell's OWN axial reach (the 99.5th percentile of
|axial|, NOT the split threshold, so embryo size divides out) and binned into 20. The profile
falls off at both ends because a bin at the junction and a bin at the far edge each sample less
cytoplasm than one in the middle — that is CELL SHAPE, and it is exactly what 7.2 and 7.3 then
divide out.

7.2 — THE VOLCANO. Per embryo per gene: log2(contact ÷ edge), then the EMBRYO'S OWN BULK SPLIT
subtracted, so an embryo whose plane sat slightly off does not push every gene the same way.
Across embryos: a one-sample t-test against 0, BH per stage.

  ⚠️ THE RESULT IS AT CHANCE. The reference prints the hit count beside its own expectation for
  exactly this reason — 11 nominal against 10 expected at early 2-cell. That comparison is
  computed here and shipped, and the page leads with it.

7.3 — THE DENSITY MAPS. The top 15 leaning genes each way, pooled over every early 2-cell embryo
in a frame whose x is the junction→edge axis, divided by the same map built from ALL genes so the
shared two-lobed shape cancels. Illustrative, NOT a test: the sets were selected by their position
on this very axis, so the gradients necessarily run in the selected direction.

7.5 (the GO dot plot) is NOT built. It needs a gene→term annotation source that does not exist
anywhere on this machine, and 12 imported rows whose k/K/n cannot be re-derived would be a picture
of someone else's computation. The page says so rather than pretending.

VALIDATION: the reference ships its own 7.2 table, and every gene's fold and P is compared against
it at the bottom of this build.

Output: data/contacthalves.json.gz
"""
import collections
import csv
import gzip
import json
import math
import os
import sys

import numpy as np
from scipy import stats

import embryo_stats as ES

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(DATA, "contacthalves.json.gz")
REF = ("/Users/rishib/Desktop/EmbyroPlayground/HighResSlideshowExports/Index7/"
       "7.2_contact_edge_volcanoes/data.csv")

VERSION = "contacthalves-1.0.0"
STAGES = [("Early2Cell", "early2cell"), ("Late2Cell", "late2cell")]
MIN_TX = 50               # a gene needs this many transcripts in an embryo to be measured there
                          # (recovered from the reference's own 7.2 table: at 50 every gene's
                          #  n_embryos AND total_counts match it exactly, 201/201 and 126/126)
MIN_EMBRYOS = 5           # ...in at least this many embryos to be tested
NBIN = 20                 # 7.1 axial bins
REACH_PCT = 99.5          # per-blastomere axial reach, as a percentile of |axial|
NMAP = 34                 # 7.3 grid, matching the reference
TOP_MAP = 15              # genes per direction in the density maps
MIN_MAP_TX = 200          # pooled transcripts a map needs to be worth drawing
MIN_BIN_ALL = 20          # a bin needs this many ALL-GENE transcripts to be drawn at all —
                          # below it the ratio is a handful of counts over a handful of counts,
                          # which reads as structure and is not. This is what gives the maps
                          # their stepped outline: the mask IS the pooled cell.


def blastomere_frame(sc, a, b):
    """(û, o, vIn, vOut) for blastomere `a` against its sister `b`.

    û points from the junction toward `a`'s own edge, and o is the equal-volume plane's origin,
    ON the plane — which is what makes ES.split_volumes exact."""
    Va, Fa = ES.mesh_of(sc, a)
    ca = ES.vol_centroid(Va, Fa)
    cb = ES.vol_centroid(*ES.mesh_of(sc, b))
    u = ca - cb
    u /= np.linalg.norm(u)
    vol = ES.seg_volumes(sc)[a]
    n, o = ES.equal_volume_plane(Va, Fa, u, ca, exact_total=vol)
    vOut, vIn = ES.split_volumes(Va, Fa, n, o, exact_total=vol)
    return u, o, ca, vIn, vOut, Va


def main():
    man_p = os.path.join(DATA, "zygote_manifest.json")
    label_of = {}
    if os.path.isfile(man_p):
        label_of = {m["id"]: m.get("label") or m["id"]
                    for m in json.load(open(man_p))["embryos"]}
    probeset = ES.probesets()

    per = collections.defaultdict(list)        # stage -> [{id, gene, cont, edge, lfc}]
    prof = []                                  # 7.1
    emb_meta = collections.defaultdict(list)
    skipped = []
    # 7.3, pooled over early 2-cell only (the reference's own choice)
    gmap = collections.defaultdict(lambda: np.zeros((NMAP, NMAP)))
    allmap = np.zeros((NMAP, NMAP))
    map_pos = []                               # (gene, ia, ir) deferred until the top sets exist

    for stage, key in STAGES:
        for eid in ES.stage_ids(stage):
            sc = ES.read_scene(ES.scene_path(stage, eid))
            bodies = ES.classify_body(sc)
            if len(bodies) != 2:
                skipped.append({"id": eid, "stage": key,
                                "reason": f"{len(bodies)} blastomere(s), need exactly 2"})
                continue
            cont = collections.Counter()
            edge = collections.Counter()
            vIn = vOut = 0.0
            axial_all, prof_hist = [], np.zeros(NBIN)
            n_prof = 0
            for a, b in ((bodies[0], bodies[1]), (bodies[1], bodies[0])):
                u, o, ca, vi, vo, Va = blastomere_frame(sc, a, b)
                vIn += vi; vOut += vo
                TX = ES.cytoplasm_positions(sc, a)
                if not TX:
                    continue
                P = np.concatenate(list(TX.values())) if TX else np.zeros((0, 3))
                # the cell's own axial reach, measured from the SPLIT PLANE outward
                t_all = (P - o) @ u
                reach = float(np.percentile(np.abs(t_all), REACH_PCT)) or 1.0
                for g, Q in TX.items():
                    t = (Q - o) @ u
                    c = int((t <= 0).sum())
                    cont[g] += c
                    edge[g] += len(t) - c
                # 7.1: x runs 0 (contact plane) -> 1 (cell edge), 0.5 = the equal-volume split.
                # The split sits at reach*? in raw units, so normalise the junction-side extent
                # onto [0, 0.5] and the edge-side onto [0.5, 1] against this cell's own reach.
                tin = t_all[t_all <= 0]
                tout = t_all[t_all > 0]
                rin = float(np.percentile(np.abs(tin), REACH_PCT)) if len(tin) else 1.0
                rout = float(np.percentile(tout, REACH_PCT)) if len(tout) else 1.0
                xin = 0.5 - np.clip(np.abs(tin) / (rin or 1.0), 0, 1) * 0.5
                xout = 0.5 + np.clip(tout / (rout or 1.0), 0, 1) * 0.5
                x = np.concatenate([xin, xout])
                prof_hist += np.histogram(np.clip(x, 0, 0.999999),
                                          bins=NBIN, range=(0, 1))[0]
                n_prof += len(x)

                if key == "early2cell":
                    # 7.3: x is the junction->edge axis (0.5 at the split), y the transverse
                    # distance from the axis, both normalised by this cell's own extents
                    d = P - o
                    ax_ = d @ u
                    rad = np.linalg.norm(d - np.outer(ax_, u), axis=1)
                    rr = float(np.percentile(rad, REACH_PCT)) or 1.0
                    xin_ = 0.5 - np.clip(np.abs(ax_[ax_ <= 0]) / (rin or 1.0), 0, 1) * 0.5
                    xout_ = 0.5 + np.clip(ax_[ax_ > 0] / (rout or 1.0), 0, 1) * 0.5
                    xx = np.empty(len(ax_)); xx[ax_ <= 0] = xin_; xx[ax_ > 0] = xout_
                    yy = np.clip(rad / rr, 0, 0.999999)
                    ia = np.clip((xx * NMAP).astype(int), 0, NMAP - 1)
                    ir = np.clip((yy * NMAP).astype(int), 0, NMAP - 1)
                    np.add.at(allmap, (ia, ir), 1)
                    off = 0
                    for g, Q in TX.items():
                        sl = slice(off, off + len(Q)); off += len(Q)
                        map_pos.append((g, ia[sl], ir[sl]))

            if not cont and not edge:
                skipped.append({"id": eid, "stage": key, "reason": "no cytoplasmic transcripts"})
                continue

            # the embryo's OWN bulk split, subtracted from every gene below
            bulk = math.log2((sum(cont.values()) + ES.EPS) / (sum(edge.values()) + ES.EPS))
            for g in set(cont) | set(edge):
                n = cont[g] + edge[g]
                if n < MIN_TX:
                    continue
                per[key].append({
                    "id": eid, "gene": g, "cont": cont[g], "edge": edge[g], "n": n,
                    "lfc": math.log2((cont[g] + ES.EPS) / (edge[g] + ES.EPS)) - bulk})
            if n_prof:
                prof.append({"id": eid, "stage": key, "n": int(n_prof),
                             "f": [round(float(v), 6) for v in prof_hist / n_prof]})
            emb_meta[key].append({
                "id": eid, "label": label_of.get(eid, eid), "probeset": probeset.get(eid, "?"),
                "n_tx": int(sum(cont.values()) + sum(edge.values())),
                "vIn": round(float(vIn), 1), "vOut": round(float(vOut), 1),
                "bulk": round(float(bulk), 4)})
            print(f"  {key:10s} {eid:30s} {sum(cont.values()) + sum(edge.values()):6d} tx  "
                  f"split {vIn / (vIn + vOut):.6f}  bulk {bulk:+.4f}")

    # ---- 7.2 ----
    stages_out = {}
    for _, key in STAGES:
        by = collections.defaultdict(list)
        for r in per[key]:
            by[r["gene"]].append(r)
        rows = []
        for g, sub in by.items():
            if len(sub) < MIN_EMBRYOS:
                continue
            y = np.array([r["lfc"] for r in sub], float)
            t, p = stats.ttest_1samp(y, 0.0)
            rows.append({"g": g, "n": len(sub), "total": int(sum(r["n"] for r in sub)),
                         "lfc": round(float(y.mean()), 5), "sd": round(float(y.std(ddof=1)), 5),
                         "p": float(p),
                         "per": [{"id": r["id"], "lfc": round(r["lfc"], 4),
                                  "c": r["cont"], "e": r["edge"]}
                                 for r in sorted(sub, key=lambda r: -abs(r["lfc"]))]})
        if rows:
            for r, q in zip(rows, ES.bh(np.array([r["p"] for r in rows]))):
                r["q"] = float(q)
        rows.sort(key=lambda r: r["p"])
        for i, r in enumerate(rows, start=1):
            r["rank"] = i
            r["side"] = "contact" if r["lfc"] > 0 else "edge"
        n_nom = sum(1 for r in rows if r["p"] < 0.05)
        stages_out[key] = {
            "genes": rows, "embryos": emb_meta[key],
            "n_nominal": n_nom, "expected": round(0.05 * len(rows), 1),
            "n_fdr": sum(1 for r in rows if r.get("q", 1) < 0.05)}

    # ---- 7.3 ----
    e2c = stages_out["early2cell"]["genes"]
    lean = {"contact": [r["g"] for r in sorted(e2c, key=lambda r: -r["lfc"])[:TOP_MAP]],
            "edge": [r["g"] for r in sorted(e2c, key=lambda r: r["lfc"])[:TOP_MAP]]}
    sets = {k: set(v) for k, v in lean.items()}
    grids = {k: np.zeros((NMAP, NMAP)) for k in lean}
    for g, ia, ir in map_pos:
        for k, S in sets.items():
            if g in S:
                np.add.at(grids[k], (ia, ir), 1)
    bg = allmap / max(allmap.sum(), 1)
    maps = {}
    for k, M in grids.items():
        tot = M.sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            f = M / max(tot, 1)
            L = np.log2(np.where((f > 0) & (bg > 0), f / np.where(bg > 0, bg, 1), np.nan))
        L = np.where(allmap >= MIN_BIN_ALL, L, np.nan)
        maps[k] = {"n": int(tot), "genes": lean[k],
                   "z": [[None if not np.isfinite(v) else round(float(v), 3) for v in row]
                         for row in L],
                   "enough": bool(tot >= MIN_MAP_TX)}

    # ---- validation against the reference's own 7.2 table ----
    val = {"available": False, "path": os.path.basename(REF)}
    if os.path.isfile(REF):
        ref = collections.defaultdict(dict)
        with open(REF) as fh:
            for row in csv.DictReader(fh):
                ref[row["stage"]][row["gene"]] = (float(row["lfc_contact_vs_edge"]),
                                                 float(row["p"]))
        val = {"available": True, "path": os.path.basename(REF), "stages": {}}
        for st, key in STAGES:
            mine = {r["g"]: (r["lfc"], r["p"]) for r in stages_out[key]["genes"]}
            shared = sorted(set(mine) & set(ref.get(st, {})))
            if not shared:
                val["stages"][key] = {"n_shared": 0}
                continue
            dl = np.array([abs(mine[g][0] - ref[st][g][0]) for g in shared])
            a = np.array([mine[g][0] for g in shared])
            b = np.array([ref[st][g][0] for g in shared])
            val["stages"][key] = {
                "n_shared": len(shared), "n_mine": len(mine), "n_ref": len(ref.get(st, {})),
                "median_abs_diff": round(float(np.median(dl)), 5),
                "max_abs_diff": round(float(dl.max()), 5),
                "r": round(float(np.corrcoef(a, b)[0, 1]), 5),
                "sign_agree": round(float((np.sign(a) == np.sign(b)).mean()), 4)}

    doc = {
        "meta": {
            "version": VERSION,
            "method": "figures 7.1, 7.2 and 7.3 — the reference's own contact-region definition",
            "params": {"MIN_TX": MIN_TX, "MIN_EMBRYOS": MIN_EMBRYOS, "NBIN": NBIN,
                       "NMAP": NMAP, "TOP_MAP": TOP_MAP, "REACH_PCT": REACH_PCT,
                       "MIN_BIN_ALL": MIN_BIN_ALL},
            "plane": "each blastomere is split by a plane PERPENDICULAR to its own junction->edge "
                     "axis and SLID until that blastomere's two halves hold equal volume — not "
                     "the centroid plane; blastomeres are flattened at the junction and rounded "
                     "at the edge, so the two differ and the difference would read as a "
                     "gene-independent contact bias",
            "contact": "the junction-side half of each blastomere; two per embryo",
            "excluded": "nuclei and the polar body, by segment label — a nucleus near the "
                        "junction would otherwise read as contact enrichment",
            "centring": "each embryo's own bulk log2(contact/edge) is subtracted from every gene, "
                        "so an embryo whose plane sat slightly off does not push every gene the "
                        "same way",
            "profile": "7.1: position normalised per blastomere by that cell's OWN axial reach "
                       f"({REACH_PCT}th percentile of |axial|), NOT the split threshold, so "
                       "embryo size divides out. The fall-off at both ends is cell shape.",
            "maps": "7.3: pooled density of the top 15 leaning genes each way, divided by the "
                    "all-gene map so the shared two-lobed shape cancels. ILLUSTRATIVE, NOT A "
                    "TEST — the sets were selected by their position on this very axis.",
            "no_go": "figure 7.5's GO dot plot is not built: it needs a gene->term annotation "
                     "source that is not present on this machine, and importing 12 rows whose "
                     "k/K/n cannot be re-derived would be a picture of someone else's "
                     "computation rather than a result of this one",
            "skipped": skipped,
            "validation": val,
            "n_embryos": {k: len(emb_meta[k]) for _, k in STAGES},
        },
        "stages": stages_out,
        "profile": prof,
        "maps": maps,
    }
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(f"\n  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e6:.2f} MB)")
    for _, key in STAGES:
        s = stages_out[key]
        print(f"  {key}: {len(s['embryos'])} embryos, {len(s['genes'])} genes, "
              f"{s['n_nominal']} nominal vs {s['expected']} expected, {s['n_fdr']} FDR")
    if val.get("available"):
        for _, key in STAGES:
            v = val["stages"][key]
            if v.get("n_shared"):
                print(f"  vs reference {key}: {v['n_shared']} shared, r={v['r']}, "
                      f"median |Δlfc| {v['median_abs_diff']}, sign agree {v['sign_agree']}")
    print("  maps: " + ", ".join(f"{k} {v['n']} tx" for k, v in maps.items()))


if __name__ == "__main__":
    main()
