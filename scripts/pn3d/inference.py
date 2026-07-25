"""
Apply the calibrated clock to the fixed MERFISH cohort, with OOD detection (task 18).

For every resolved zygote the dimensionless geometry (distance_sum / MEASURED cell
radius) is fed to the probabilistic clock, giving a tau posterior + intervals +
confidence. Because the fixed radius is measured per embryo, an imaging scale
offset relative to the live-imaging cohort is normalized away (domain adaptation).

Out-of-distribution detection marks a stack as uncertain rather than assigning
false confidence when: the segmentation is unresolved / low-confidence; the
dimensionless feature falls outside the live-training support; or the stage is
not a zygote (2-cell = post-NEBD, structurally OOD for a pronuclear clock).

Nothing here is independent TIME validation — the fixed stacks have no true tau.
These are inference outputs, and agreement with any prior cached pseudotime is
NOT treated as validation.
"""
from __future__ import annotations

import numpy as np

from . import clock as C


def scheffler_support(reference_radius_um: float):
    x, _, _ = C.load_scheffler(reference_radius_um)
    return {"p01": float(np.percentile(x, 1)), "p99": float(np.percentile(x, 99)),
            "min": float(x.min()), "max": float(x.max()),
            "mean": float(x.mean()), "std": float(x.std())}


def ood_flags(sum_over_R, support, seg_status, seg_conf, stage) -> dict:
    reasons = []
    level = "in_domain"
    if stage != "zygote":
        level = "out_of_domain"; reasons.append(f"stage {stage} is not a pronuclear zygote")
    if seg_status != "resolved":
        level = "out_of_domain"; reasons.append(f"segmentation {seg_status}")
    elif seg_conf < 0.5:
        level = "caution" if level == "in_domain" else level
        reasons.append(f"low segmentation confidence {seg_conf:.2f}")
    if sum_over_R is not None:
        if sum_over_R < support["min"] or sum_over_R > support["max"]:
            level = "out_of_domain"
            reasons.append(f"sum/R {sum_over_R:.2f} outside live-training range "
                           f"[{support['min']:.2f}, {support['max']:.2f}]")
        elif sum_over_R < support["p01"] or sum_over_R > support["p99"]:
            if level == "in_domain":
                level = "caution"
            reasons.append(f"sum/R {sum_over_R:.2f} outside live p01-p99 "
                           f"[{support['p01']:.2f}, {support['p99']:.2f}]")
        z = abs(sum_over_R - support["mean"]) / (support["std"] + 1e-9)
        reasons_z = z
    else:
        reasons_z = None
    return {"ood_level": level, "ood_reasons": reasons,
            "feature_zscore_vs_live": (round(reasons_z, 2) if reasons_z is not None else None)}


def infer_record(geom_rec: dict, clk: C.ProbClock, support: dict) -> dict:
    stage = geom_rec.get("stage")
    status = geom_rec.get("status")
    conf = geom_rec.get("confidence", 0.0)
    g = geom_rec.get("geometry")
    out = {"embryo_id": geom_rec["embryo_id"], "experiment": geom_rec.get("experiment"),
           "batch_date": geom_rec.get("batch_date"), "stage": stage,
           "segmentation_status": status, "segmentation_confidence": conf,
           "flags": geom_rec.get("flags", [])}
    sum_over_R = g["sum_over_R"] if g else None
    ood = ood_flags(sum_over_R, support, status, conf, stage)
    out.update(ood)
    if g and status == "resolved":
        post = clk.predict(sum_over_R)
        out["geometry"] = {k: g[k] for k in ("sum_over_R", "near_over_R", "far_over_R",
                                             "inter_over_R", "cell_radius_um", "distance_sum_um",
                                             "vol_asymmetry", "polar_body_present")}
        out["pseudotime"] = post
        out["inferable"] = True
    else:
        out["inferable"] = False
        out["pseudotime"] = None
    return out
