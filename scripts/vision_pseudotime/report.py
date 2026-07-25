"""
Redacted provenance / metrics report for the website research view (brief item 10).

Assembles a SMALL, SAFE JSON: no absolute paths (redacted via config.redact_path),
no pixels, no raw data. The website reads this for provenance, preprocessing,
extracted geometry, predicted tau + uncertainty, held-out metrics, and explicit
pilot / OOD / blocker warnings.

Measured numbers (TIFF shapes, projection stats, hybrid cross-check) are filled
in only when the raw data is present on this machine; otherwise the report marks
them unavailable rather than inventing them.
"""
from __future__ import annotations

import csv
import glob
import json
import os
from datetime import datetime, timezone

import numpy as np

from . import __version__, config, encoder, evaluate, hybrid, manifest, movies, splits, tiff_audit

REPORT_PATH = os.path.join(config.REPO_ROOT, "data", "vision_pseudotime_report.json")


# ───────────────────────────── measured pieces (only if data present) ─────────────────────────────
def audit_section() -> list[dict]:
    out = []
    for sid, path, kind in tiff_audit.stack_sources():
        a = tiff_audit.audit(path, declared_kind=kind)
        keep = ("basename", "path_redacted", "exists", "shape", "axes", "dtype",
                "stack_role", "est_full_load_bytes", "spacing_trustworthy", "needs_sidecar",
                "sidecar_present", "spacing_status", "metadata_sufficient", "voxel_um",
                "channels", "label_classes", "software")
        out.append({"source_id": sid, **{k: a.get(k) for k in keep}})
    return out


def projection_section() -> dict:
    """Run the streaming 2.5D build if the matched stacks are present. No pixels stored."""
    ipath = config.resolve_source("target_intensity_example")
    lpath = config.resolve_source("target_labels_example")
    if not (os.path.isfile(ipath) and os.path.isfile(lpath)):
        return {"available": False, "note": "matched stacks not present on this machine"}
    from . import projection as P
    tensor, meta = P.build_2p5d(ipath, lpath, out=256, z_step=1)
    stats = {meta["channels"][c]: {"mean": round(float(tensor[..., c].mean()), 4),
                                   "p05": round(float(np.percentile(tensor[..., c], 5)), 4),
                                   "p95": round(float(np.percentile(tensor[..., c], 95)), 4)}
             for c in range(tensor.shape[-1])}
    return {"available": True, "channels": meta["channels"], "bbox": meta["bbox"],
            "out": meta["out"], "z_extent_labels": meta.get("z_extent_labels"),
            "no_raw_sum": meta["no_raw_sum"], "n_pages_read": meta.get("n_pages_read_labels"),
            "peak_memory_note": "streamed one z-page at a time; full-volume load avoided",
            "channel_stats": stats,
            "preview_pngs_dir_redacted": config.redact_path(config.DERIVED_DIR),
            "previews_published": False}


def movie_section() -> dict:
    rows = movies._movie_rows()
    inv, counts = [], {"pilot_normal_dev": 0, "ood": 0, "excluded": 0, "needs_panel_isolation": 0}
    for r in rows:
        inc = movies.ROLE_INCLUSION.get(r.get("default_role", ""), "excluded")
        counts[inc] = counts.get(inc, 0) + 1
        inv.append({"movie_id": r["movie_id"], "default_role": r.get("default_role"),
                    "inclusion": inc, "frames": r.get("frames"),
                    "condition": r.get("condition_or_subject", "")[:70]})
    return {"n_movies": len(rows), "inclusion_counts": counts,
            "normal_dev_pilot_movies": [m["movie_id"] for m in inv if m["inclusion"] == "pilot_normal_dev"],
            "tau_from_frame_index": "forbidden without annotated pn-formation + NEBD",
            "inventory": inv}


def hybrid_section(max_zygotes: int = 6) -> dict:
    """Cross-check the hybrid path against the deployed cache on real zygotes."""
    cache_p = os.path.join(config.REPO_ROOT, "calibration_data", "fixed_cohort_geometry.csv")
    fx_p = os.path.join(config.REPO_ROOT, "data", "pronuclei_pseudotime.json")
    out = {"clock_reused": hybrid._model()["model_version"], "validated_zygotes": 0}
    if not (os.path.isfile(cache_p) and os.path.isfile(fx_p)):
        out["note"] = "cache/published tau not present"
    else:
        import build_pronuclei as BP
        cache = {r["id"]: r for r in csv.DictReader(open(cache_p)) if r.get("nearer_to_center_um")}
        fx = {r["id"]: r for r in json.load(open(fx_p))["embryos"]}
        dfeat_max = dtau_max = 0.0
        n = 0
        for eid in cache:
            lab = glob.glob(os.path.join(BP.SRC, eid, "*_label.tif"))
            if not lab:
                continue
            sub = hybrid._load_sub(lab[0])
            vox = (BP.DS_XY * BP.XY_UM, BP.DS_XY * BP.XY_UM, BP.DS_Z * BP.Z_UM)
            g = hybrid.segment_to_geometry(sub, vox)
            if not g["ok"]:
                continue
            fc = {k: float(cache[eid][k]) for k in hybrid.FEATURE_COLS}
            dfeat_max = max(dfeat_max, max(abs(g["features"][k] - fc[k]) for k in hybrid.FEATURE_COLS))
            pred = hybrid.hybrid_predict(g["features"], "verified", True)
            dtau_max = max(dtau_max, abs(pred["tau"] - fx[eid]["tau"]))
            n += 1
            if n >= max_zygotes:
                break
        out.update(validated_zygotes=n, max_feature_error_um=round(dfeat_max, 6),
                   max_tau_error=round(dtau_max, 6),
                   verdict=("hybrid reproduces the deployed pipeline exactly"
                            if dfeat_max < 1e-4 and dtau_max < 1e-4 else "MISMATCH — investigate"))
    # example stack: must be blocked
    lpath = config.resolve_source("target_labels_example")
    if os.path.isfile(lpath):
        r = hybrid.hybrid_from_label_stack(lpath)
        out["example_stack"] = {"status": r["status"], "blockers": r.get("blockers"),
                                "spacing_status": r.get("spacing_status"),
                                "classes_verified": r.get("classes_verified")}
    return out


def encoder_feasibility_demo() -> dict:
    """Deterministic SYNTHETIC demonstration that the direct-encoder harness works,
    plus the real-data refusal. Clearly labelled synthetic — no real claim."""
    imgs, taus, groups = encoder.synthetic_dataset(6, 10, seed=1, overlay=True)
    box = encoder.overlay_box(imgs[0].shape[0])
    stripped = []
    for im in imgs:
        c = im.copy(); y0, x0, y1, x1 = box; c[y0:y1, x0:x1] = 0; stripped.append(c)
    smap = splits.split_groups(set(groups.tolist()), (0.6, 0.0, 0.4), seed=7)
    tr = [i for i, g in enumerate(groups) if smap[g] == "train"]
    te = [i for i, g in enumerate(groups) if smap[g] == "test"]
    enc = encoder.ExploratoryImageEncoder().fit([stripped[i] for i in tr], taus[tr], groups[tr])
    pred = enc.predict_many([stripped[i] for i in te])
    im0 = stripped[te[0]]
    sal = evaluate.occlusion_saliency(lambda a: enc.predict(a), im0, patch=8, stride=6)
    body = im0.mean(axis=2) > 0.1 * im0.mean(axis=2).max()
    o_pred = evaluate.corner_overlay_baseline([stripped[i] for i in te], taus[tr],
                                              [stripped[i] for i in tr])
    refusal = None
    try:
        encoder.ExploratoryImageEncoder().fit([im0], [float("nan")], ["single_real_embryo"])
    except encoder.InsufficientData as e:
        refusal = str(e)
    return {
        "status": encoder.MODEL_STATUS, "data": "SYNTHETIC phantoms (labelled synthetic)",
        "min_independent_embryos_required": encoder.MIN_INDEPENDENT_EMBRYOS,
        "held_out_embryo_grouped_mae": round(evaluate.embryo_grouped_mae(
            taus[te], pred, groups[te])["macro_mae"], 4),
        "held_out_spearman": round(evaluate.spearman(taus[te], pred), 4),
        "held_out_pair_order_accuracy": round(evaluate.pair_order_accuracy(
            taus[te], pred)["strict_accuracy"], 4),
        "overlay_baseline_mae": round(float(np.mean(np.abs(taus[te] - o_pred))), 4),
        "saliency_corner_mass_after_strip": round(evaluate.saliency_box_mass(sal, box), 4),
        "saliency_embryo_contrast": round(evaluate.saliency_contrast(sal, body), 3),
        "real_data_refusal": refusal,
        "interpretation": ("harness learns synthetic tau and ignores the stripped overlay; on the "
                           "REAL material it refuses — 1 embryo, no per-frame tau."),
    }


# ───────────────────────────── assembly ─────────────────────────────
def assemble(now_iso: str | None = None) -> dict:
    man = manifest.validate_all()
    now = now_iso or datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    return {
        "meta": {
            "package_version": __version__,
            "geometry_baseline_model": hybrid._model()["model_version"],
            "generated_at_utc": now,
            "status": "PILOT — data pipeline + feasibility, not a validated image clock",
            "no_absolute_paths": True, "no_raw_pixels_in_report": True,
        },
        "manifest": {"ok": man["ok"],
                     "sources": man["sources"]["entries"],
                     "movies_summary": {"n": man["movies"]["n"],
                                        "errors": man["movies"]["errors"]}},
        "tiff_audits": audit_section(),
        "preprocessing": projection_section(),
        "movies": movie_section(),
        "hybrid_geometry_clock": hybrid_section(),
        "exploratory_encoder": encoder_feasibility_demo(),
        "evaluation_metrics_available": [
            "embryo_grouped_mae", "spearman", "pair_order_accuracy", "interval_coverage",
            "segmentation_dice", "segmentation_iou", "centroid_error", "ood_rejection",
            "trivial_baselines(constant/brightness/corner-overlay)", "occlusion_saliency"],
        "warnings": [
            "PILOT: no independently validated end-to-end image clock exists yet.",
            "Fixed snapshots have no true tau; they support segmentation/OOD, not supervised labels.",
            "Augmented views of one embryo are NOT independent samples.",
            "Public rendered movies cannot be reliably mapped to the 53 numeric trajectories.",
            "The example lab stack is BLOCKED: unverified voxel spacing and unverified label classes.",
        ],
        "data_blockers": [
            "Verified voxel spacing (xy pixel size + z step) for the lab stacks.",
            "Verified label-class semantics (which ID is cell / pronucleus).",
            "Multiple independent untreated RAW time-lapse embryos with pn-formation + NEBD frames.",
            "A held-out embryo/batch test set and frozen preprocessing.",
        ],
        "acquisition_checklist": "calibration_data/vision_pseudotime/ACQUISITION_CHECKLIST.md",
        "model_card": "calibration_data/vision_pseudotime/MODEL_CARD.md",
    }


def write_report(path: str = REPORT_PATH, now_iso: str | None = None) -> dict:
    rep = assemble(now_iso)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(rep, fh, indent=1)
    return rep


if __name__ == "__main__":
    r = write_report()
    print(f"wrote {config.rel_to_repo(REPORT_PATH)} "
          f"({os.path.getsize(REPORT_PATH) / 1024:.0f} KB)")
    print(f"  manifest ok={r['manifest']['ok']}  "
          f"hybrid validated={r['hybrid_geometry_clock'].get('validated_zygotes')} zygotes "
          f"(max Δτ={r['hybrid_geometry_clock'].get('max_tau_error')})")
    print(f"  example stack: {r['hybrid_geometry_clock'].get('example_stack')}")
    print(f"  encoder demo (synthetic): grouped_mae="
          f"{r['exploratory_encoder']['held_out_embryo_grouped_mae']}  refuses_real="
          f"{r['exploratory_encoder']['real_data_refusal'] is not None}")
