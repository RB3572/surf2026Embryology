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


def scheffler_support(reference_radius_um: float, feature: str = "rms"):
    x, _, _ = C.load_scheffler(reference_radius_um, feature)
    return {"p01": float(np.percentile(x, 1)), "p99": float(np.percentile(x, 99)),
            "min": float(x.min()), "max": float(x.max()),
            "mean": float(x.mean()), "std": float(x.std())}


def ood_flags(sum_over_R, support, seg_status, seg_conf, stage) -> dict:
    reasons = []
    level = "in_domain"
    if stage != "zygote":
        level = "out_of_domain"; reasons.append(f"stage {stage} is not a pronuclear zygote")
    if seg_status == "single_pronucleus":
        # scored, but the clock was calibrated on two-pronucleus geometry, so this
        # is an extrapolation: flagged, never silently mixed in with the rest
        level = "out_of_domain"
        reasons.append("single annotated pronucleus — geometry extrapolates beyond the "
                       "two-pronucleus configuration the clock was calibrated on")
    elif seg_status != "resolved":
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


def single_pn_sensitivity(d_over_R: float, clk: C.ProbClock, peer_d_over_R) -> dict:
    """
    For a single-annotated-pronucleus embryo, the honest question is: what if a
    second pronucleus IS there but was not labelled?

    The model input is rms = sqrt((d1² + d2²)/2)/R. Sweeping the missing d2 over
    the distances actually observed in the two-pronucleus cohort gives the range
    of tau that hypothesis spans. That converts an unquantified annotation risk
    into a stated interval, instead of a narrow clock-only interval that would
    look more confident than the sample deserves.
    """
    import numpy as np
    peers = np.asarray(list(peer_d_over_R), float)
    if peers.size == 0:
        return {}
    rms_alt = np.sqrt((d_over_R ** 2 + peers ** 2) / 2.0)
    taus = np.array([clk.predict(float(r))["tau_mean"] for r in rms_alt])
    return {
        "if_second_pronucleus_present": {
            "tau_lo": round(float(np.percentile(taus, 2.5)), 4),
            "tau_hi": round(float(np.percentile(taus, 97.5)), 4),
            "tau_median": round(float(np.median(taus)), 4),
            "n_hypotheses": int(peers.size),
            "note": "range of tau if the unlabelled second pronucleus sat at any distance "
                    "observed in the two-pronucleus cohort; this is annotation uncertainty, "
                    "which the calibrated clock interval does NOT cover",
        }
    }


def infer_record(geom_rec: dict, clk: C.ProbClock, support: dict,
                 peer_d_over_R=None) -> dict:
    stage = geom_rec.get("stage")
    status = geom_rec.get("status")
    conf = geom_rec.get("confidence", 0.0)
    g = geom_rec.get("geometry")
    out = {"embryo_id": geom_rec["embryo_id"], "experiment": geom_rec.get("experiment"),
           "batch_date": geom_rec.get("batch_date"), "stage": stage,
           "segmentation_status": status, "segmentation_confidence": conf,
           "n_pronuclei": (g or {}).get("n_pronuclei"),
           "flags": geom_rec.get("flags", [])}
    # the model's universal input: defined for one OR two pronuclei
    x = g.get("rms_over_R") if g else None
    ood = ood_flags(x, support, status, conf, stage)
    out.update(ood)
    if g and x is not None:
        post = clk.predict(x)
        keep = ("n_pronuclei", "rms_over_R", "rms_to_center_um", "cell_radius_um",
                "pron_vol_frac", "polar_body_present", "sum_over_R", "near_over_R",
                "far_over_R", "inter_over_R", "distance_sum_um", "vol_asymmetry")
        out["geometry"] = {k: g[k] for k in keep if k in g}
        if g.get("n_pronuclei") == 1 and peer_d_over_R is not None:
            d1 = (g.get("pron_distances_um") or [0])[0] / max(g.get("cell_radius_um") or 1, 1e-6)
            post.update(single_pn_sensitivity(d1, clk, peer_d_over_R))
        out["pseudotime"] = post
        out["inferable"] = True
    else:
        out["inferable"] = False
        out["pseudotime"] = None
    return out
