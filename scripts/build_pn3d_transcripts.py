#!/usr/bin/env python3
"""
Join the 3D pseudotime model's tau to per-embryo transcript counts.

Produces data/pn3d_transcripts.json: one row per zygote that appears in BOTH
cohorts (strict 1:1 crosswalk), carrying the calibrated tau posterior + interval
+ QC/OOD from the 3D model and the transcript totals / per-gene counts.

Also computes, per gene, the correlation of its count with tau across embryos —
the pseudotime analogue of the existing gene<->distance ranking.

Everything the page needs is in this one committed file; no raw data, no
absolute paths. The mounted dataset is NOT required to run this.

Usage: python3 scripts/build_pn3d_transcripts.py
"""
from __future__ import annotations

import gzip
import json
import math
import os
import sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from scripts.pn3d import crosswalk  # noqa: E402

INFER = os.path.join(HERE, "data", "pn3d", "inference.json")
MODEL = os.path.join(HERE, "data", "pn3d", "model.json")
GEOM = os.path.join(HERE, "data", "pn3d", "segmentation_geometry.json")
GENES = os.path.join(HERE, "data", "pronuclei_genes.json.gz")
MANIFEST = os.path.join(HERE, "data", "pronuclei_manifest.json")
TX_GEOM = os.path.join(HERE, "calibration_data", "fixed_cohort_geometry.csv")
OUT = os.path.join(HERE, "data", "pn3d_transcripts.json")

MIN_EMBRYOS_PER_GENE = 8          # below this a correlation is not worth ranking
XY_UM = 0.15                      # transcript centroids are stored in plot units


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / math.sqrt(sxx * syy)


def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(ranks(xs), ranks(ys))


def main() -> int:
    import csv
    for p in (INFER, MODEL, GEOM, GENES, MANIFEST, TX_GEOM):
        if not os.path.isfile(p):
            raise SystemExit(f"missing {os.path.relpath(p, HERE)}")

    inf = json.load(open(INFER))
    model = json.load(open(MODEL))
    geo = {e["embryo_id"]: e for e in json.load(open(GEOM))["embryos"]}
    genes_doc = json.load(gzip.open(GENES, "rt"))
    man = {e["id"]: e for e in json.load(open(MANIFEST))["embryos"]}
    tx_rows = {e["id"]: e for e in genes_doc["embryos"]}

    zygotes = [e for e in inf["embryos"] if e["stage"] == "zygote"]
    by_id = {e["embryo_id"]: e for e in zygotes}

    # cell-body centroids (µm) in the shared field frame, both cohorts
    pn_pts = {}
    for e in zygotes:
        g = (geo.get(e["embryo_id"]) or {}).get("geometry")
        if g and g.get("cell_center_um"):
            pn_pts[e["embryo_id"]] = (g["cell_center_um"][0], g["cell_center_um"][1])
    tx_pts = {}
    for r in csv.DictReader(open(TX_GEOM)):
        if r.get("center_plot_x") and r["id"] in tx_rows:
            tx_pts[r["id"]] = (float(r["center_plot_x"]) * XY_UM,
                               float(r["center_plot_y"]) * XY_UM)

    cw = crosswalk.build_spatial(pn_pts, tx_pts)
    agree = crosswalk.name_agreement(cw["matched"])

    rows = []
    for pn_id, m in sorted(cw["matched"].items()):
        tx_id = m["tx_id"]
        e = by_id[pn_id]
        if not e.get("pseudotime"):
            continue                                  # segmentation unresolved -> no tau
        t = tx_rows.get(tx_id, {})
        mrow = man.get(tx_id, {})
        pt = e["pseudotime"]
        rows.append({
            "pn3d_id": pn_id, "tx_id": tx_id,
            "centroid_error_um": m["centroid_error_um"],
            "label": mrow.get("label", tx_id), "date_short": mrow.get("date_short", ""),
            "batch_date": e.get("batch_date"), "experiment": e.get("experiment"),
            "tau": pt["tau_mean"], "tau_sd": pt["tau_sd"],
            "lo50": pt["interval_50"][0], "hi50": pt["interval_50"][1],
            "lo95": pt["interval_95"][0], "hi95": pt["interval_95"][1],
            "confidence": pt["confidence"],
            "qc": e.get("segmentation_status"), "ood": e.get("ood_level"),
            "seg_confidence": e.get("segmentation_confidence"),
            "sum_over_R": (e.get("geometry") or {}).get("sum_over_R"),
            "cell_radius_um": (e.get("geometry") or {}).get("cell_radius_um"),
            "legacy_pron_distance_um": m.get("distance"),
            "total": int(t.get("total") or m.get("total") or 0),
            "genes": t.get("genes", {}),
        })

    # per-gene correlation with tau (across the joined embryos)
    all_genes: dict = {}
    for r in rows:
        for g, c in r["genes"].items():
            all_genes.setdefault(g, 0)
            all_genes[g] += 1
    # The per-embryo gene dicts are SPARSE: a gene absent from an embryo means zero
    # transcripts were detected there, which is a real measurement, not missing data.
    # Both readings are legitimate and give different correlations, so BOTH are
    # computed and the page uses one consistently for the plot AND the ranking:
    #   detected : only embryos where the gene was detected  (n_detected)
    #   all      : every joined embryo, absent counted as 0  (n_all)
    # Counts also scale with sequencing depth (totals span ~4k-500k here), so the
    # depth-corrected "fraction of total" variants are reported alongside.
    gene_stats = []
    for g in sorted(all_genes):
        det = [(r["tau"], float(r["genes"].get(g, 0)), r["total"]) for r in rows if g in r["genes"]]
        alln = [(r["tau"], float(r["genes"].get(g, 0)), r["total"]) for r in rows]
        if len(det) < MIN_EMBRYOS_PER_GENE:
            continue
        def stats(pts):
            t = [p[0] for p in pts]
            c = [p[1] for p in pts]
            fr = [p[1] / p[2] if p[2] else 0.0 for p in pts]
            return (_r(pearson(t, c)), _r(spearman(t, c)),
                    _r(pearson(t, fr)), _r(spearman(t, fr)))
        rc_d, rho_d, rfc_d, rhof_d = stats(det)
        rc_a, rho_a, rfc_a, rhof_a = stats(alln)
        gene_stats.append({
            "gene": g,
            "n_detected": len(det), "n_all": len(alln),
            "total_count": int(sum(p[1] for p in det)),
            # detected-only
            "r_count": rc_d, "rho_count": rho_d, "r_frac": rfc_d, "rho_frac": rhof_d,
            # all joined embryos, absent = 0
            "r_count_all": rc_a, "rho_count_all": rho_a,
            "r_frac_all": rfc_a, "rho_frac_all": rhof_a,
        })
    gene_stats.sort(key=lambda s: -(abs(s["rho_count"]) if s["rho_count"] is not None else -1))

    taus = [r["tau"] for r in rows]
    payload = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "model": {
            "package_version": model["package_version"],
            "clock": model["clock"]["kind"],
            "feature": model["clock"]["feature"],
            "reference_radius_um": model["reference_radius_um"],
            "cv_metrics": model["clock"]["cv_metrics"],
            "time_supervision": "Scheffler 2021 live-imaging (53 embryos, 2057 frames); "
                                "leave-one-embryo-out CV. Fixed stacks have no true time.",
        },
        "join": {
            "method": "spatial: mutual-nearest cell-body centroid within the same field of view "
                      f"(max {cw['max_match_um']} µm); names used only as a cross-check",
            "n_pn3d_zygotes": len(zygotes),
            "n_transcript_embryos": len(tx_rows),
            "n_joined": len(rows),
            "n_ambiguous": len(cw["ambiguous"]),
            "n_unmatched_pn3d": len(cw["unmatched_pn3d"]),
            "n_unmatched_tx": len(cw["unmatched_tx"]),
            "centroid_error_um": {
                "median": round(sorted(r["centroid_error_um"] for r in rows)[len(rows) // 2], 4),
                "max": round(max(r["centroid_error_um"] for r in rows), 4),
            } if rows else None,
            "name_cross_check": agree,
            "ambiguous": cw["ambiguous"],
            "unmatched_pn3d": cw["unmatched_pn3d"][:60],
            "note": "A field of view can hold several distinct embryos. The transcript side splits "
                    "them into _N_0/_N_1; the mount segments one at the FOV level and each extra "
                    "one into a nested sub-directory (segmentation only, sharing the FOV image). "
                    "They are matched by cell-body centroid, so each embryo keeps its own "
                    "pseudotime and its own transcripts.",
        },
        "tau_range": {"min": min(taus), "max": max(taus)} if taus else None,
        "qc_counts": _counts(rows, "ood"),
        "min_embryos_per_gene": MIN_EMBRYOS_PER_GENE,
        "gene_stat_note": ("Gene dicts are sparse: absence means zero detected, which is a real "
                           "measurement. 'detected' statistics use only embryos where the gene was "
                           "seen; 'all' statistics use every joined embryo with absent = 0. The "
                           "page applies one choice to both the plot and the ranking so they "
                           "always agree. Raw counts scale with sequencing depth, so the "
                           "fraction-of-total variants are the depth-corrected reading."),
        "n_genes_ranked": len(gene_stats),
        "gene_stats": gene_stats,
        "embryos": rows,
    }
    json.dump(payload, open(OUT, "w"), separators=(",", ":"))
    kb = os.path.getsize(OUT) / 1024
    print(f"joined {len(rows)} embryos  (ambiguous {len(cw['ambiguous'])}, "
          f"unmatched pn3d {len(cw['unmatched_pn3d'])})")
    print(f"ranked {len(gene_stats)} genes (>= {MIN_EMBRYOS_PER_GENE} embryos)")
    print(f"tau range {min(taus):.3f}-{max(taus):.3f} · QC {payload['qc_counts']}")
    print(f"wrote {os.path.relpath(OUT, HERE)} ({kb:.0f} KB)")
    return 0


def _r(v):
    return None if v is None else round(v, 4)


def _counts(rows, key):
    out: dict = {}
    for r in rows:
        out[r.get(key)] = out.get(r.get(key), 0) + 1
    return out


if __name__ == "__main__":
    sys.exit(main())
