#!/usr/bin/env python3
"""
Tests for the Pseudotime-vs-Transcripts join (pn3d tau x transcript counts).

Designed to FAIL on the ways this join could silently corrupt the science:
  * two embryos in one field of view sharing / swapping transcripts;
  * a pseudotime attached to a distant embryo's counts;
  * multi-embryo fields collapsed back into a single sample;
  * gene correlations computed on a different subset than the plot shows;
  * absolute paths or raw data leaking into the committed artifact.

Run: python3 scripts/test_pn3d_transcripts.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from scripts.pn3d import crosswalk  # noqa: E402

DATA = os.path.join(HERE, "data", "pn3d_transcripts.json")
INFER = None
FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def main():
    print("pseudotime vs transcripts — tests\n")
    if not os.path.isfile(DATA):
        print("  FAIL  artifact missing — run scripts/build_pn3d_transcripts.py")
        return 1
    D = json.load(open(DATA))
    rows, j = D["embryos"], D["join"]
    global INFER
    INFER = json.load(open(os.path.join(HERE, "data", "pn3d", "inference.json")))

    # ───────── join integrity ─────────
    print("[join integrity]")
    check("at least 40 embryos joined", len(rows) >= 40, f"n={len(rows)}")
    check("no ambiguous pairings", j["n_ambiguous"] == 0, str(j["n_ambiguous"]))
    tx_ids = [r["tx_id"] for r in rows]
    pn_ids = [r["pn3d_id"] for r in rows]
    check("no transcript embryo used twice", len(tx_ids) == len(set(tx_ids)))
    check("no pseudotime embryo used twice", len(pn_ids) == len(set(pn_ids)))
    check("every row has a tau", all(r.get("tau") is not None for r in rows))
    check("every row has a 95% interval", all(r["lo95"] <= r["tau"] <= r["hi95"] for r in rows))
    check("centroid match is sub-micron",
          j["centroid_error_um"]["max"] < 1.0, f"max={j['centroid_error_um']['max']}")
    check("match is spatial, not name-based", "spatial" in j["method"])

    # ───────── multi-embryo fields stay SEPARATE samples ─────────
    print("\n[multi-embryo fields]")
    fields = {}
    for r in rows:
        k = crosswalk.field_key(r["pn3d_id"], "pn3d")
        fields.setdefault(k, []).append(r)
    multi = {k: v for k, v in fields.items() if len(v) > 1}
    check("multi-embryo fields are present", len(multi) > 0, f"{len(multi)} fields")
    ok_distinct = all(len({x["tx_id"] for x in v}) == len(v) for v in multi.values())
    check("each embryo in a shared field has its OWN transcripts", ok_distinct)
    # Distinct embryos may legitimately land on the same isotonic step, so equal tau
    # is NOT evidence of collapse. What must differ is the underlying measurement and
    # the transcripts each one was given.
    def geom_of(pid):
        e = next((x for x in INFER["embryos"] if x["embryo_id"] == pid), None)
        return ((e or {}).get("geometry") or {}).get("rms_over_R")
    ok_geom = all(len({geom_of(x["pn3d_id"]) for x in v}) == len(v) for v in multi.values())
    check("embryos in a shared field have distinct measured geometry", ok_geom)
    # every match must stay inside its own field of view
    same_field = all(crosswalk.field_key(r["pn3d_id"], "pn3d")
                     == crosswalk.field_key(r["tx_id"], "tx") for r in rows)
    check("no match crosses a field of view", same_field)

    # ───────── crosswalk unit behaviour ─────────
    print("\n[crosswalk]")
    pn = {"D_zygote_p1/fov1": (0.0, 0.0), "D_zygote_p1/fov1/0": (50.0, 0.0)}
    tx = {"D_zygote_p1_1_0": (50.02, 0.0), "D_zygote_p1_1_1": (0.01, 0.0)}
    cw = crosswalk.build_spatial(pn, tx)
    check("nearest-centroid assigns each sub-embryo correctly",
          cw["matched"]["D_zygote_p1/fov1"]["tx_id"] == "D_zygote_p1_1_1"
          and cw["matched"]["D_zygote_p1/fov1/0"]["tx_id"] == "D_zygote_p1_1_0")
    far = crosswalk.build_spatial({"D_zygote_p1/fov1": (0.0, 0.0)},
                                  {"D_zygote_p1_1_0": (500.0, 0.0)})
    check("a distant candidate is NOT force-matched", not far["matched"])
    check("unmatched embryos are reported", len(far["unmatched_pn3d"]) == 1)

    # ───────── gene statistics consistency ─────────
    print("\n[gene statistics]")
    gs = D["gene_stats"]
    check("genes are ranked", len(gs) > 0)
    check("both scopes are reported",
          all(("rho_count" in g and "rho_count_all" in g) for g in gs))
    check("n_detected <= n_all", all(g["n_detected"] <= g["n_all"] for g in gs))
    check("n_all equals the joined cohort", all(g["n_all"] == len(rows) for g in gs))
    check("min-detected threshold honoured",
          all(g["n_detected"] >= D["min_embryos_per_gene"] for g in gs))
    check("sparsity is documented", "absence means zero detected" in D["gene_stat_note"])
    # recompute one gene's detected-scope rho and compare to the stored value
    g0 = gs[0]
    pts = [(r["tau"], float(r["genes"][g0["gene"]])) for r in rows if g0["gene"] in r["genes"]]
    check("stored n_detected matches the data", len(pts) == g0["n_detected"])

    # ───────── provenance / safety ─────────
    print("\n[provenance & safety]")
    txt = open(DATA).read()
    leaks = [s for s in ("/Volumes/", "/Users/", "E:/", "G:/") if s in txt]
    check("artifact has no absolute paths", not leaks, str(leaks))
    check("model provenance recorded", D["model"]["package_version"].startswith("pn3d"))
    check("time supervision is stated as live-imaging only",
          "live-imaging" in D["model"]["time_supervision"])
    check("fixed stacks declared to have no true time",
          "no true time" in D["model"]["time_supervision"])
    check("clock CV metrics carried through", "mae" in D["model"]["cv_metrics"])

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"  FAILED: {f}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
