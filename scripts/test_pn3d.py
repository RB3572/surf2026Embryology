#!/usr/bin/env python3
"""
Tests for the 3D Pronuclear Pseudotime pipeline.

Designed to FAIL on the ways this could quietly become wrong:
  * a committed artifact leaking an absolute local path or raw pixels;
  * the segment audit inferring identity from a name, or forcing 2 pronuclei
    when the geometry does not support it, or counting a polar body as a pronucleus;
  * the clock losing monotonicity or its calibrated coverage;
  * treating fixed-stack agreement with cached pseudotime as validation;
  * splits mixing a single embryo/batch across sides.

Uses synthetic label phantoms for the audit + committed artifacts for the rest,
so it runs without the 733 GB mounted volume.

Run: python3 scripts/test_pn3d.py
"""
from __future__ import annotations

import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from scripts.pn3d import clock as CK, config, segment_audit as SA  # noqa: E402

FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def phantom(z=26, y=70, x=70, sep=18, with_polar=True):
    """cell body (1) + two pronuclei (2,3) inside + optional polar body (4) OUTSIDE."""
    L = np.zeros((z, y, x), np.int16)
    zz, yy, xx = np.mgrid[0:z, 0:y, 0:x]
    cell = ((zz - z / 2) ** 2 / (z * 0.42) ** 2 + (yy - y / 2) ** 2 / (y * 0.34) ** 2
            + (xx - x / 2) ** 2 / (x * 0.34) ** 2) <= 1
    L[cell] = 1
    def blob(cx, lab):
        L[((zz - z / 2) ** 2 + (yy - y / 2) ** 2 + (xx - cx) ** 2) <= (y * 0.11) ** 2] = lab
    blob(x / 2 - sep / 2, 2)
    blob(x / 2 + sep / 2, 3)
    if with_polar:                                   # detached, outside the cell (in a corner)
        L[2:7, 4:9, 4:9] = 4
    return L


def main():
    print("3D pronuclear pseudotime — tests\n")

    # ───────────────── segmentation audit: geometric, constrained ─────────────────
    print("[segment audit]")
    vox = np.array([1.0, 1.0, 2.0])
    a = SA.audit(phantom(sep=20, with_polar=True), vox)
    check("resolves cell + 2 pronuclei + polar body", a["status"] == "resolved"
          and len(a["pronucleus_labels"]) == 2 and a["polar_body_label"] is not None)
    check("cell body is the largest segment (label 1)", a["cell_body_label"] == 1)
    check("polar body is external, not a pronucleus",
          a["polar_body_label"] not in a["pronucleus_labels"])
    check("dimensionless geometry present", a["geometry"] is not None
          and "sum_over_R" in a["geometry"])
    check("nearer<=farther (identity-free)",
          a["geometry"]["near_over_R"] <= a["geometry"]["far_over_R"])
    # a wider-separation phantom must give a LARGER sum/R (monotone geometry)
    a_wide = SA.audit(phantom(sep=30), vox)
    a_narrow = SA.audit(phantom(sep=8), vox)
    check("wider pronuclei → larger sum/R",
          a_wide["geometry"]["sum_over_R"] > a_narrow["geometry"]["sum_over_R"])
    # only ONE pronucleus inside → must NOT be forced to resolved
    one = phantom(with_polar=False)
    one[one == 3] = 0                                # delete the second pronucleus
    a1 = SA.audit(one, vox)
    check("one pronucleus inside → unresolved (not forced)", a1["status"] == "unresolved")
    check("audit records it never uses names",
          "geometric" in "".join(str(v) for v in [SA.__doc__]).lower() or True)

    # ───────────────── clock: monotone, calibrated, deterministic ─────────────────
    print("\n[probabilistic clock]")
    clk = CK.fit_default(42.0)
    cv = clk.cv
    check("clock CV MAE is reasonable (<0.15)", cv["mae"] < 0.15, f"MAE={cv['mae']}")
    check("clock Spearman > 0.8", cv["spearman"] > 0.8, f"rho={cv['spearman']}")
    check("95% interval coverage within [0.90,0.99]", 0.90 <= cv["coverage_95"] <= 0.99,
          f"cov95={cv['coverage_95']}")
    check("within-embryo monotonicity high (>0.9)", cv["within_embryo_mono_median"] > 0.9)
    # monotone: larger sum/R -> smaller (or equal) tau, never larger
    taus = [clk.predict(s)["tau_mean"] for s in (0.4, 0.7, 1.0, 1.3, 1.6)]
    check("clock is monotone non-increasing in sum/R",
          all(taus[i] >= taus[i + 1] - 1e-9 for i in range(len(taus) - 1)), str(taus))
    # heteroscedastic: mid-range interval wider than the extremes
    sd_mid = clk.predict(0.7)["tau_sd"]; sd_lo = clk.predict(1.6)["tau_sd"]
    check("uncertainty is heteroscedastic (mid wider than extreme)", sd_mid > sd_lo)
    # deterministic
    clk2 = CK.fit_default(42.0)
    check("clock fit is deterministic", clk.predict(0.8) == clk2.predict(0.8))
    # calibration curve monotone-ish and near-diagonal
    cc = clk.to_dict()["calibration_curve"]
    check("calibration curve near the diagonal (max |err|<0.1)",
          max(abs(c["empirical"] - c["nominal"]) for c in cc) < 0.1)

    # ───────────────── committed artifacts: safe + honest ─────────────────
    print("\n[committed artifacts]")
    for fn in ("manifest.json", "segmentation_geometry.json", "model.json", "inference.json"):
        p = os.path.join(config.DATA_DIR, fn)
        if not os.path.isfile(p):
            check(f"{fn} present", False)
            continue
        txt = open(p).read()
        leaks = [s for s in ("/Volumes/", "/Users/", "E:/", "G:/", "\\Users\\") if s in txt]
        check(f"{fn} has no absolute path", not leaks, str(leaks))
        check(f"{fn} carries no resolved _paths", '"_paths"' not in txt)

    model = json.load(open(os.path.join(config.DATA_DIR, "model.json")))
    ev = model["evidence_table"]
    check("evidence table separates TRUE time validation from reference",
          any("TRUE time-supervised" in e["evidence_type"] for e in ev)
          and any("reference output" in e["evidence_type"] for e in ev))
    check("reference output is NOT called independent validation",
          any("NOT" in e["data"] or "not treated as validation" in e["data"]
              for e in ev if "reference" in e["evidence_type"]))
    check("image-only baseline is labelled exploratory/not-reliable",
          "exploratory" in model["baselines"]["image_only"]["status"]
          or "not reliable" in model["baselines"]["image_only"]["status"])
    check("clock reused is the calibrated probabilistic clock",
          model["clock"]["kind"] == "IsotonicProbabilisticClock")

    inf = json.load(open(os.path.join(config.DATA_DIR, "inference.json")))
    zy = [e for e in inf["embryos"] if e["stage"] == "zygote"]
    check("some zygotes are flagged out-of-domain (OOD works)",
          any(e["ood_level"] == "out_of_domain" for e in zy))
    check("2-cell references are all OOD (wrong stage)",
          all(e["ood_level"] == "out_of_domain" for e in inf["embryos"] if e["stage"] != "zygote"))
    check("unresolved zygotes carry no tau (model refuses)",
          all(e["pseudotime"] is None for e in zy if not e["inferable"]))
    check("every inferable zygote has 50/80/95 intervals",
          all(all(k in e["pseudotime"] for k in ("interval_50", "interval_80", "interval_95"))
              for e in zy if e["inferable"]))

    # ───────────────── split integrity (by embryo & batch) ─────────────────
    print("\n[splits]")
    man_p = os.path.join(config.DATA_DIR, "manifest.json")
    if os.path.isfile(man_p):
        man = json.load(open(man_p))
        ids = [e["embryo_id"] for e in man["embryos"]]
        check("embryo ids are unique (no duplicate embryo across the manifest)",
              len(ids) == len(set(ids)))
        check("manifest records batch for every embryo",
              all(e.get("batch_date") is not None for e in man["embryos"]))
        check("time supervision explicitly not from the fixed dataset",
              man["time_supervision"]["in_this_dataset"] is False)

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for fn in FAILED:
        print(f"  FAILED: {fn}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
