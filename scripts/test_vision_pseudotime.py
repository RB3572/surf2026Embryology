#!/usr/bin/env python3
"""
Tests for the fixed-image pseudotime PILOT pipeline (vision_pseudotime).

Designed to FAIL on the specific ways this pipeline could quietly become wrong:
  * an absolute local path or raw pixels leaking into the committed report
  * a TIFF audit loading pixels, or guessing intensity-vs-label from dtype
  * a raw SUM projection sneaking in, or a non-deterministic projection
  * a label augmentation inventing label values, or a geometry-breaking transform
    being enabled by default
  * a frame or augmented derivative of one embryo crossing a split boundary
  * tau being inferred from a frame index with no annotated pn-formation/NEBD
  * an overlay/timestamp region surviving into model pixels
  * the hybrid clock drifting from the deployed calibration
  * the exploratory encoder training on too few embryos / perturbation frames

Synthetic fixtures keep the suite independent of the 2.5 GB stacks; a few extra
checks run only if the real data happens to be present.

Run:  python3 scripts/test_vision_pseudotime.py
"""
from __future__ import annotations

import json
import os
import re
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from scripts.vision_pseudotime import (augment, config, encoder, evaluate,  # noqa: E402
                                       hybrid, manifest, movies, projection, report,
                                       splits, tiff_audit)

FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def write_tiff(path, arr):
    import tifffile
    tifffile.imwrite(path, arr)


# ───────────────────────────── synthetic label phantom ─────────────────────────────
def phantom_labels(z=24, y=64, x=64, sep=18):
    """A 3-label volume: cytoplasm(1) + two pronuclei(2,3) separated by `sep` in x."""
    vol = np.zeros((z, y, x), np.int16)
    zz, yy, xx = np.mgrid[0:z, 0:y, 0:x]
    cell = ((zz - z / 2) ** 2 / (z * 0.4) ** 2 + (yy - y / 2) ** 2 / (y * 0.4) ** 2
            + (xx - x / 2) ** 2 / (x * 0.4) ** 2) <= 1
    vol[cell] = 1
    def blob(cx, lab):
        m = ((zz - z / 2) ** 2 + (yy - y / 2) ** 2 + (xx - cx) ** 2) <= (y * 0.12) ** 2
        vol[m] = lab
    blob(x / 2 - sep / 2, 2)
    blob(x / 2 + sep / 2, 3)
    return vol


def main():
    print("fixed-image pseudotime pilot — tests\n")

    # ───────────────────────────── config / manifest ─────────────────────────────
    print("[config & manifest]")
    check("redact_path hides absolute paths",
          "/Users" not in config.redact_path("/Users/x/all.tif")
          and config.redact_path("/Users/x/all.tif").endswith("all.tif"))
    check("placeholder passes through redaction", config.redact_path("not_yet_available") == "not_yet_available")
    check("is_inside_repo true for repo path", config.is_inside_repo(config.SOURCES_CSV))
    check("is_inside_repo false for /tmp", not config.is_inside_repo("/tmp/x.tif"))
    rep = manifest.validate_all()
    check("committed manifests validate cleanly", rep["ok"],
          str(rep["sources"]["errors"] + rep["movies"]["errors"]))
    fixed_sup = [e for e in rep["sources"]["entries"]
                 if e["source_kind"] in ("confocal_stack", "segmentation_stack")
                 and e["supervised_tau_use"] == "yes"]
    check("no fixed snapshot is a supervised tau source", not fixed_sup)
    # a source that puts a LARGE stack inside the repo must be rejected
    with tempfile.TemporaryDirectory() as td:
        bad_csv = os.path.join(td, "sources.csv")
        inside = os.path.join(config.REPO_ROOT, "README.md")   # a real path inside the repo
        with open(bad_csv, "w") as fh:
            fh.write("source_id,source_kind,location,independent_embryos,time_truth,"
                     "pixels_or_geometry,supervised_tau_use,intended_use,limitations\n")
            fh.write(f"bad_stack,confocal_stack,{inside},1,none,pixels,no,test,inside-repo\n")
        orig = config.SOURCES_CSV
        try:
            config.SOURCES_CSV = bad_csv
            manifest.config.SOURCES_CSV = bad_csv
            r = manifest.validate_sources()
            check("large source inside the repo is rejected",
                  any("INSIDE the repo" in e for e in r["errors"]), str(r["errors"]))
        finally:
            config.SOURCES_CSV = orig
            manifest.config.SOURCES_CSV = orig

    # ───────────────────────────── TIFF audit + sidecar ─────────────────────────────
    print("\n[tiff audit]")
    with tempfile.TemporaryDirectory() as td:
        ip = os.path.join(td, "syn_intensity.tif"); lp = os.path.join(td, "syn_labels.tif")
        write_tiff(ip, (np.random.default_rng(0).random((10, 32, 32)) * 500).astype(np.uint16))
        write_tiff(lp, phantom_labels(10, 32, 32).astype(np.int16))
        ai = tiff_audit.audit(ip, declared_kind="confocal_stack")
        al = tiff_audit.audit(lp, declared_kind="segmentation_stack")
        check("audit reports shape without loading pixels", ai["shape"] == (10, 32, 32))
        check("audit reports est full-load bytes", ai["est_full_load_bytes"] == 10 * 32 * 32 * 2)
        check("intensity vs label from manifest, not dtype",
              ai["stack_role"] == "intensity" and al["stack_role"] == "labels"
              and ai["dtype"] == "uint16")
        check("stack with untrustworthy tags needs a sidecar", ai["needs_sidecar"])
        check("no verified spacing => metadata not sufficient", not ai["metadata_sufficient"])
        # sidecar validation
        check("sidecar missing embryo_id is flagged",
              "embryo_id" in " ".join(tiff_audit.validate_sidecar({"voxel_um": [1, 1, 1]})))
        check("sidecar bad voxel_um is flagged",
              any("voxel_um" in e for e in tiff_audit.validate_sidecar(
                  {"embryo_id": "e", "voxel_um": [1, 1]})))
        check("null voxel_um allowed (unknown scale)",
              tiff_audit.validate_sidecar({"embryo_id": "e", "voxel_um": None}) == [])

    # ───────────────────────────── projection + label handling ─────────────────────────────
    print("\n[projection & labels]")
    check("no raw sum projection is offered", "sum" not in projection.PROJECTIONS)
    with tempfile.TemporaryDirectory() as td:
        ip = os.path.join(td, "i.tif"); lp = os.path.join(td, "l.tif")
        vol = phantom_labels(24, 64, 64)
        inten = (vol > 0).astype(np.uint16) * 300 + (np.random.default_rng(1).random((24, 64, 64))
                                                      * 20).astype(np.uint16)
        write_tiff(ip, inten); write_tiff(lp, vol.astype(np.int16))
        acc = projection.scan(lp, is_labels=True)
        occ = projection.occupancy_projection(acc)
        bnd = projection.boundary_projection(acc)
        check("occupancy in [0,1]", float(occ.min()) >= 0 and float(occ.max()) <= 1)
        check("occupancy nonzero only where labels exist", occ.max() > 0)
        check("boundary is binary-ish and bounded", float(bnd.max()) <= 1)
        t1, m1 = projection.build_2p5d(ip, lp, out=48, z_step=1)
        t2, _ = projection.build_2p5d(ip, lp, out=48, z_step=1)
        check("2.5D tensor shape (out,out,3)", t1.shape == (48, 48, 3))
        check("2.5D values in [0,1]", float(t1.min()) >= 0 and float(t1.max()) <= 1)
        check("2.5D projection is deterministic", np.array_equal(t1, t2))
        check("2.5D declares no-raw-sum", m1["no_raw_sum"] is True)
        check("MIP baseline is single channel", projection.single_channel_mip_baseline(
            ip, out=48)[0].shape == (48, 48))
        # streaming never holds the whole volume: robust mean uses per-plane accumulators
        check("robust-mean projection defined", projection.robust_mean_projection(
            projection.scan(ip, is_labels=False)).shape == (64, 64))

    # ───────────────────────────── segmentation -> geometry ─────────────────────────────
    print("\n[hybrid geometry]")
    vol = phantom_labels(24, 64, 64, sep=20)
    g = hybrid.segment_to_geometry(vol, (0.6, 0.6, 2.0))
    check("phantom geometry extracted", g["ok"], g.get("reason", ""))
    check("symmetric features present",
          all(k in g["features"] for k in hybrid.FEATURE_COLS))
    check("nearer <= farther (sorted, identity-free)",
          g["features"]["nearer_to_center_um"] <= g["features"]["farther_to_center_um"])
    two_label = np.where(vol == 3, 1, vol)          # collapse to cyto + one pronucleus
    g2 = hybrid.segment_to_geometry(two_label, (0.6, 0.6, 2.0))
    check("fewer than two pronuclei is refused", not g2["ok"])
    # blocked hybrid when spacing/classes unverified
    hb = hybrid.hybrid_predict(g["features"], spacing_status="missing", classes_verified=False)
    check("hybrid blocks without verified spacing+classes", hb["status"] == "blocked"
          and len(hb["blockers"]) == 2)
    hp = hybrid.hybrid_predict(g["features"], spacing_status="verified", classes_verified=True)
    check("hybrid predicts with verified inputs", hp["status"] == "predicted"
          and 0 <= hp["tau"] <= 1)

    # hybrid clock reproduces the DEPLOYED calibration on the committed fixed cohort
    fx_p = os.path.join(HERE, "data", "pronuclei_pseudotime.json")
    if os.path.isfile(fx_p):
        fx = json.load(open(fx_p))["embryos"]
        worst = 0.0
        for r in fx:
            if not r.get("features"):
                continue
            f = {k: r["features"][k] for k in hybrid.FEATURE_COLS}
            worst = max(worst, abs(hybrid.geometry_to_tau(f)["tau"] - r["tau"]))
        check("hybrid clock matches deployed tau within rounding", worst <= 1e-4,
              f"max|Δτ|={worst:.2e}")

    # ───────────────────────────── movies / timestamp masking ─────────────────────────────
    print("\n[movies & overlays]")
    frame = np.zeros((60, 60, 3), np.uint8)
    frame[5:55, 5:55] = 80                          # content
    frame[0:6, 0:6] = 255                           # a bright corner "timestamp"
    clean, rec = movies.strip_overlays(frame, overlay_boxes=[[0, 0, 0.2, 0.2]])
    y0, x0, y1, x1 = rec["overlay_boxes_zeroed"][0]
    check("uniform border cropped", rec["clean_hw"][0] < 60)
    check("overlay box zeroed in model pixels", float(clean[y0:y1, x0:x1].max()) == 0)
    check("strip method is not OCR (honest)", "not OCR" in rec["method"])
    check("tau None without annotated frames", movies.tau_from_annotation(5, None, None) is None)
    check("tau None if only one boundary", movies.tau_from_annotation(5, 0, None) is None)
    check("tau computed with both boundaries",
          abs(movies.tau_from_annotation(5, 0, 10) - 0.5) < 1e-9)
    panels = movies.split_panels(np.zeros((40, 80)), rows=1, cols=2)
    check("panel split yields 2 panels", len(panels) == 2 and panels[0].shape == (40, 40))
    check("perturbation role maps to non-training inclusion",
          movies.ROLE_INCLUSION["ood_only"] == "ood"
          and movies.ROLE_INCLUSION["exclude_from_tau_training"] == "excluded")

    # ───────────────────────────── augmentation ─────────────────────────────
    print("\n[augmentation]")
    img = np.random.default_rng(2).random((48, 48, 3)).astype(np.float32)
    a1, ap1 = augment.augment(img, 99)
    a2, _ = augment.augment(img, 99)
    a3, _ = augment.augment(img, 100)
    check("augment deterministic for a seed", np.array_equal(a1, a2))
    check("augment varies with seed", not np.array_equal(a1, a3))
    check("augment output in [0,1]", float(a1.min()) >= 0 and float(a1.max()) <= 1)
    lab = np.zeros((48, 48), int); lab[10:30, 12:34] = 1; lab[18:24, 20:26] = 2
    _, al2, _ = augment.augment_with_labels(img, lab, 5)
    check("label augmentation invents no new labels", set(np.unique(al2)).issubset({0, 1, 2}))
    check("label augmentation deterministic",
          np.array_equal(augment.augment_with_labels(img, lab, 5)[1], al2))
    refused = 0
    for bad in augment.DISABLED_GEOMETRY_BREAKING:
        try:
            augment.augment(img, 1, enabled=(bad,))
        except augment.GeometryBreakingTransform:
            refused += 1
    check("all geometry-breaking transforms refused by default",
          refused == len(augment.DISABLED_GEOMETRY_BREAKING))

    # ───────────────────────────── split leakage (the core item-6 test) ─────────────────────────────
    print("\n[split leakage]")
    samples = [{"embryo_id": f"emb{e}", "frame": f, "panel": 0}
               for e in range(8) for f in range(15)]
    assigned = splits.assign_samples(samples, seed=20260724)
    recs = splits.expand_with_augmentation(assigned, n_aug=5, aug_seed=20260724)
    chk = splits.check_no_leakage(recs)
    check("no group appears in multiple splits", not chk["cross_split_groups"], str(chk["cross_split_groups"]))
    check("no augmented derivative crosses its parent split", not chk["augmented_split_mismatches"])
    check("leakage check overall ok", chk["ok"])
    # every embryo entirely within one split
    from collections import defaultdict
    per = defaultdict(set)
    for r in recs:
        per[r["embryo_id"]].add(r["split"])
    check("each embryo in exactly one split", all(len(v) == 1 for v in per.values()))
    # determinism across recompute (stable hash, not salted hash())
    recs2 = splits.expand_with_augmentation(
        splits.assign_samples(samples, seed=20260724), 5, 20260724)
    check("split assignment is reproducible",
          [r["split"] for r in recs] == [r["split"] for r in recs2])
    check("splits happen BEFORE augmentation (order enforced by API)",
          all("split" in a for a in assigned))

    # ───────────────────────────── evaluation + saliency ─────────────────────────────
    print("\n[evaluation]")
    check("dice identity", evaluate.dice(np.ones((4, 4)), np.ones((4, 4))) == 1.0)
    check("iou disjoint is zero", evaluate.iou(np.array([1, 0]), np.array([0, 1])) == 0.0)
    check("coverage counts inside interval",
          evaluate.interval_coverage([0.5], [0.4], [0.6]) == 1.0)
    po = evaluate.pair_order_accuracy([0.1, 0.2, 0.3], [0.1, 0.2, 0.3])
    check("pair-order perfect on monotone", po["strict_accuracy"] == 1.0)
    ood = evaluate.ood_rejection_rate([1, 1, 0, 0], [1, 0, 0, 0])
    check("ood recall computed", abs(ood["ood_recall"] - 0.5) < 1e-9)
    gm = evaluate.embryo_grouped_mae([0, 1, 0, 1], [0, 1, 1, 0], ["a", "a", "b", "b"])
    check("grouped MAE macro over groups", abs(gm["macro_mae"] - 0.5) < 1e-9 and gm["n_groups"] == 2)

    # ───────────────────────────── exploratory encoder ─────────────────────────────
    print("\n[exploratory encoder]")
    imgs, taus, groups = encoder.synthetic_dataset(5, 8, seed=3)
    try:
        encoder.ExploratoryImageEncoder().fit(imgs[:8], taus[:8], ["one"] * 8)
        one = False
    except encoder.InsufficientData:
        one = True
    check("encoder refuses < min independent embryos", one)
    try:
        encoder.ExploratoryImageEncoder().fit(imgs, [float("nan")] * len(imgs), groups)
        nan_ok = False
    except encoder.InsufficientData:
        nan_ok = True
    check("encoder refuses missing per-frame tau", nan_ok)
    try:
        encoder.ExploratoryImageEncoder().fit(imgs, taus, groups,
                                              treatments=["1 uM nocodazole"] * len(imgs))
        pert = False
    except encoder.InsufficientData:
        pert = True
    check("encoder refuses perturbation frames", pert)
    enc = encoder.ExploratoryImageEncoder().fit(imgs, taus, groups)
    enc_b = encoder.ExploratoryImageEncoder().fit(imgs, taus, groups)
    p1, p2 = enc.predict_many(imgs), enc_b.predict_many(imgs)
    check("encoder fit+predict deterministic", np.allclose(p1, p2) and np.allclose(enc.coef_, enc_b.coef_))
    check("encoder marked feasibility-only", enc.status == "exploratory-feasibility")

    # ───────────────────────────── report redaction ─────────────────────────────
    print("\n[report redaction]")
    rp = report.assemble(now_iso="2026-01-01T00:00:00+00:00")
    js = json.dumps(rp)
    for pat in ("/Users/", "/home/", "iCloud", "Mobile Documents"):
        check(f"report has no absolute path fragment {pat!r}", pat not in js)
    check("report declares pilot status", "PILOT" in rp["meta"]["status"])
    check("report lists data blockers", len(rp["data_blockers"]) >= 3)
    check("report warns fixed snapshots have no tau",
          any("no true tau" in w or "no supervised" in w.lower() for w in rp["warnings"]))

    # ───────────────────────────── optional: real data present ─────────────────────────────
    ipath = config.resolve_source("target_intensity_example")
    lpath = config.resolve_source("target_labels_example")
    if os.path.isfile(ipath):
        print("\n[real data — present]")
        a = tiff_audit.audit(ipath, declared_kind="confocal_stack")
        check("real intensity audit is metadata-only fast", a["shape"] is not None
              and a["est_full_load_bytes"] > 1e9)
        if os.path.isfile(lpath):
            hbk = hybrid.hybrid_from_label_stack(lpath)
            check("real example stack is honestly blocked", hbk["status"] == "blocked")

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"  FAILED: {f}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
