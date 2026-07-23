#!/usr/bin/env python3
"""
What does pseudotime noise do to transcript R-squared?

This is the analysis Matt asked for in the 23 Jul meeting: calibrate the clock against videos with
known real time, quantify its error on held-out embryos, and then use that measured error to say
what transcript-vs-pseudotime R-squared one should EXPECT — so a value like MERVL R^2 ~ 0.1-0.15 can
be interpreted instead of dismissed.

THREE PARTS
  1. CLOCK RELIABILITY   from the frozen clock's UNTOUCHED nested outer-test predictions only.
  2. ATTENUATION CEILING generic simulation: for a grid of LATENT true R^2 values, generate an
                         outcome from TRUE tau, then regress it on the EMPIRICAL predicted tau and
                         record the OBSERVED R^2. Answers "how much can a real trend be attenuated
                         by this clock's measured error?"
  3. FROZEN ILLUSTRATION predeclared genes (MuERV-L + ZGA markers) against calibrated tau on the
                         fixed cohort, placed on the attenuation curve.

STRICT SEPARATION.  Transcript counts are read ONLY in part 3, and only after the clock, its
feature set, its hyperparameters and its uncertainty model are already frozen on disk. Nothing in
parts 1-2 touches transcripts. The clock never saw a transcript at any stage.

HONESTY REQUIREMENTS BUILT IN
  * Pseudotime error is NOT claimed to be the only cause of low transcript R^2 — measurement noise,
    batch/probe-set effects, biological heterogeneity and small n all contribute, and the artifact
    says so explicitly.
  * The live-cohort error is a LOWER BOUND on the fixed-cohort error (domain shift, and fixed
    snapshots are segmented by a different pipeline).
  * Genes are predeclared, not chosen after seeing results; an all-gene scan is reported separately
    with FDR correction and is explicitly labelled discovery, not confirmation.

Usage:  python3 scripts/analyze_pseudotime_noise_ceiling.py [--seed 20260723] [--quick]
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAL_DIR = os.path.join(HERE, "data", "pseudotime_calibration")
FIXED_P = os.path.join(HERE, "data", "pronuclei_pseudotime.json")
GENES_P = os.path.join(HERE, "data", "pronuclei_genes.json.gz")

# Predeclared BEFORE running: the confirmatory gene is MuERV-L (the MERVL element Matt named).
# The others are the standard ZGA markers already used elsewhere in this project. This list is
# fixed in source, so it cannot be tuned to the result.
PREDECLARED = ["MuERV-L", "Zscan4d", "Zfp352", "Obox3"]

# EXPLICIT, STABLE per-scenario seed offsets. Never derive these from hash(): Python randomises
# string hashing per process (PYTHONHASHSEED), so hash-derived seeds silently make the artifact
# non-reproducible across fresh runs even with the same --seed. Any new scenario must be given an
# offset here; _scenario_seed() raises rather than falling back to something unstable.
SCENARIO_SEED_OFFSET = {"validation_uniform": 101, "fixed_cohort_matched": 202}


def _scenario_seed(base, name):
    if name not in SCENARIO_SEED_OFFSET:
        raise KeyError(f"scenario '{name}' has no stable seed offset in SCENARIO_SEED_OFFSET; "
                       "add one rather than deriving a seed from hash()")
    return int(base) + SCENARIO_SEED_OFFSET[name]
CONFIRMATORY = "MuERV-L"
MIN_N = 8                      # below this an embryo-level fit is not interpretable


def r2_of(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    r = float(np.corrcoef(x, y)[0, 1])
    return r * r


def load_frozen():
    cal = json.load(open(os.path.join(CAL_DIR, "calibration.json")))
    frames = json.load(gzip.open(os.path.join(CAL_DIR, "oof_frames.json.gz"), "rt"))
    return cal, frames


# ═══════════════════════════ 1 · clock reliability (outer-test only) ═══════════════════════════
def reliability(cal, frames, seed, B):
    fr = frames["frames"]
    y = np.array([f["tau_true"] for f in fr], float)
    p = np.array([f["tau_pred"] for f in fr], float)
    g = np.array([f["embryo_id"] for f in fr])
    dur = np.array([f["migration_duration_h"] for f in fr], float)
    err = np.abs(p - y)

    r = float(np.corrcoef(y, p)[0, 1])
    rho = float(stats.spearmanr(y, p).statistic)
    # hours conversion is valid ONLY for the live cohort, where each embryo's duration is observed
    err_h = err * dur

    rng = np.random.default_rng(seed)
    embs = np.unique(g)
    idx = {e: np.flatnonzero(g == e) for e in embs}
    boots = {"pearson_r": [], "pearson_r2": [], "spearman": [], "macro_mae": [],
             "median_ae": [], "median_ae_hours": []}
    for _ in range(B):
        pick = rng.choice(embs, size=len(embs), replace=True)
        sel = np.concatenate([idx[e] for e in pick])
        rr = float(np.corrcoef(y[sel], p[sel])[0, 1])
        boots["pearson_r"].append(rr)
        boots["pearson_r2"].append(rr * rr)
        boots["spearman"].append(float(stats.spearmanr(y[sel], p[sel]).statistic))
        boots["macro_mae"].append(float(np.mean([err[idx[e]].mean() for e in pick])))
        boots["median_ae"].append(float(np.median(err[sel])))
        boots["median_ae_hours"].append(float(np.median(err_h[sel])))
    ci = {k: [round(float(np.quantile(v, .025)), 5), round(float(np.quantile(v, .975)), 5)]
          for k, v in boots.items()}

    nested = cal["nested_evaluation"]["outer_test_metrics"]
    return {
        "source": "nested outer-test predictions only (untouched evaluation embryos)",
        "model_version": cal["meta"]["model_version"],
        "n_embryos": int(len(embs)), "n_frames": int(len(fr)),
        "pearson_r": round(r, 5), "pearson_r2": round(r * r, 5),
        "spearman": round(rho, 5),
        "strict_pairwise_ordering": nested["pooled_ordering"]["strict_accuracy"],
        "macro_within_embryo_ordering": nested["macro_within_embryo_ordering"],
        "tie_rate": nested["pooled_ordering"]["tie_rate"],
        "inversion_rate": nested["pooled_ordering"]["inversion_rate"],
        "macro_mae_tau": nested["macro_mae"], "median_ae_tau": nested["median_ae"],
        "rmse_tau": nested["pooled_rmse"],
        "median_ae_hours": round(float(np.median(err_h)), 4),
        "mean_ae_hours": round(float(err_h.mean()), 4),
        "hours_caveat": ("Hours are computed as |tau error| x that embryo's OBSERVED migration "
                         "duration, which exists only for the live-imaging validation cohort. "
                         "Fixed MERFISH zygotes have NO known absolute duration, so their tau "
                         "cannot be converted to hours."),
        "bootstrap_ci95_by_embryo": ci,
    }


# ═══════════════════════════ 2 · attenuation / noise ceiling ═══════════════════════════
def _calibrate_noise(sig, t_true, target_r2, rng):
    """Choose the noise sd that TARGETS `target_r2` for the R^2 against TRUE tau, under the SAME
    model form later used against predicted tau (OLS on a single predictor).

    With standardized signal s and outcome = s + e, e ~ N(0, sd):
        E[R^2(outcome ~ t_true)] ~= rho^2 / (1 + sd^2),   rho = corr(t_true, s)
    so sd^2 = rho^2/target - 1.

    THIS IS A TARGET, NOT AN IDENTITY. The noise is a finite random draw, so the REALIZED R^2 in any
    single replicate scatters around the target (and is very slightly upward-biased, as sample R^2
    always is). The realized value is therefore MEASURED per replicate and reported as
    `achieved_true_r2_median` / `achieved_true_r2_ci95`; the plot's x-axis uses the achieved value,
    not the request.

    For a NONLINEAR trend rho^2 < 1, which caps the achievable R^2 — v2 silently ignored this and
    mislabelled `rho^2 * target` as the latent R^2. Targets above rho^2 are now reported as NOT
    ACHIEVABLE rather than quietly produced."""
    rho = float(np.corrcoef(t_true, sig)[0, 1])
    rho2 = rho * rho
    if not np.isfinite(rho2) or rho2 <= 0:
        return None, rho2
    if target_r2 >= rho2:
        return None, rho2                     # unreachable with this trend shape
    sd = float(np.sqrt(rho2 / target_r2 - 1.0))
    return sig + rng.normal(0, sd, len(sig)), rho2


def _sampler_uniform(idx, embs, n, rng):
    """VALIDATION-UNIFORM: one uniformly random frame per sampled embryo. Because the live videos
    are sampled evenly in time, this yields an approximately UNIFORM true-tau distribution."""
    pick = rng.choice(embs, size=n, replace=True)
    return np.array([rng.choice(idx[e]) for e in pick])


def _make_matched_sampler(p_all, g_all, fixed_taus):
    """FIXED-COHORT-MATCHED: draw validation frames whose PREDICTED tau reproduces the fixed
    cohort's predicted-tau distribution.

    IMPORTANT AND UNAVOIDABLE CAVEAT: only the PREDICTOR distribution can be matched. The fixed
    cohort's TRUE tau distribution is unknown and unknowable (that is the entire reason the clock
    exists), so this scenario assumes the residual structure conditional on predicted tau transfers
    from live to fixed. It is a sensitivity analysis, not a correction."""
    order = np.argsort(p_all)
    ps = p_all[order]

    def sampler(idx, embs, n, rng):
        targets = rng.choice(fixed_taus, size=n, replace=True)
        out = []
        for t in targets:
            j = int(np.searchsorted(ps, t))
            lo, hi = max(0, j - 12), min(len(ps), j + 12)
            cand = order[lo:hi]
            out.append(int(rng.choice(cand)) if len(cand) else int(rng.integers(0, len(p_all))))
        return np.array(out)
    return sampler


def attenuation(frames, seed, B, n_snapshot, latent_grid, fixed_taus=None):
    """For each latent true R^2, simulate an outcome from TRUE tau, then measure the R^2 you would
    observe regressing it on the clock's PREDICTED tau.

    Design choices that keep this realistic rather than flattering:
      * SNAPSHOT design: one frame per embryo per replicate, because the fixed cohort is one
        snapshot per zygote — not a trajectory.
      * The (true, predicted) PAIR is resampled together, so the empirical residual structure
        (heteroscedasticity, step ties, phase-dependent bias) is preserved exactly.
      * Embryo-level resampling with replacement.
      * TWO SAMPLING SCENARIOS (see below) — R^2 attenuation depends strongly on predictor range,
        and the fixed cohort is compressed near early tau, so a single curve would be misleading.
      * The realized true-time R^2 is calibrated per replicate and REPORTED, not assumed.
    """
    fr = frames["frames"]
    y = np.array([f["tau_true"] for f in fr], float)
    p = np.array([f["tau_pred"] for f in fr], float)
    g = np.array([f["embryo_id"] for f in fr])
    embs = np.unique(g)
    idx = {e: np.flatnonzero(g == e) for e in embs}

    scenarios = {"validation_uniform": {
        "sampler": _sampler_uniform,
        "label": "validation-uniform snapshot sampling",
        "description": ("one uniformly random frame per sampled validation embryo; the live videos "
                        "are evenly sampled in time, so true tau is approximately uniform"),
        "relevance": ("characterises the clock over its full calibrated range; NOT matched to the "
                      "fixed MERFISH cohort")}}
    if fixed_taus is not None and len(fixed_taus) >= 5:
        scenarios["fixed_cohort_matched"] = {
            "sampler": _make_matched_sampler(p, g, np.asarray(fixed_taus, float)),
            "label": "fixed-cohort-matched sensitivity",
            "description": ("validation frames drawn so their PREDICTED tau reproduces the fixed "
                            "MERFISH cohort's predicted-tau distribution, out-of-domain zygotes "
                            "excluded"),
            "relevance": ("the scenario relevant to interpreting a fixed-cohort gene result such as "
                          "MuERV-L, because the fixed cohort is compressed near early tau and R^2 "
                          "depends strongly on predictor range"),
            "caveat": ("only the PREDICTED-tau distribution is matched; the fixed cohort's TRUE tau "
                       "distribution is unknown, so this assumes the residual structure conditional "
                       "on predicted tau transfers from live to fixed imaging")}

    shapes = {"linear": lambda t: t,
              "monotone_saturating": lambda t: 1.0 - np.exp(-3.0 * t)}
    out = {}
    for sc_name, sc in scenarios.items():
        rng = np.random.default_rng(_scenario_seed(seed, sc_name))
        out[sc_name] = {k: v for k, v in sc.items() if k != "sampler"}
        curves = {}
        for shape_name, f_shape in shapes.items():
            rows = []
            for r2_true in latent_grid:
                obs, ach, n_unreachable = [], [], 0
                for _ in range(B):
                    sel = sc["sampler"](idx, embs, n_snapshot, rng)
                    t_true, t_pred = y[sel], p[sel]
                    sig = f_shape(t_true)
                    if np.std(sig) == 0 or np.std(t_true) == 0:
                        continue
                    sig = (sig - sig.mean()) / np.std(sig)
                    outcome, rho2 = _calibrate_noise(sig, t_true, r2_true, rng)
                    if outcome is None:
                        n_unreachable += 1
                        continue
                    a = r2_of(t_true, outcome)          # REALIZED true-time R^2, same model form
                    o = r2_of(t_pred, outcome)          # observed, against predicted tau
                    if np.isfinite(a) and np.isfinite(o):
                        ach.append(a); obs.append(o)
                obs, ach = np.array(obs), np.array(ach)
                if len(obs) < max(10, 0.05 * B):
                    rows.append({"latent_true_r2": round(float(r2_true), 4), "achievable": False,
                                 "reason": ("requested R^2 exceeds what this trend shape can produce "
                                            "against true tau (rho^2 cap)"),
                                 "n_unreachable_replicates": n_unreachable})
                    continue
                rows.append({
                    "latent_true_r2": round(float(r2_true), 4), "achievable": True,
                    "achieved_true_r2_median": round(float(np.median(ach)), 5),
                    "achieved_true_r2_ci95": [round(float(np.quantile(ach, .025)), 5),
                                              round(float(np.quantile(ach, .975)), 5)],
                    "observed_r2_median": round(float(np.median(obs)), 5),
                    "observed_r2_ci95": [round(float(np.quantile(obs, .025)), 5),
                                         round(float(np.quantile(obs, .975)), 5)],
                    "observed_r2_q25_q75": [round(float(np.quantile(obs, .25)), 5),
                                            round(float(np.quantile(obs, .75)), 5)],
                    "attenuation_ratio": (round(float(np.median(obs) / np.median(ach)), 4)
                                          if np.median(ach) > 0 else None),
                    "n_replicates_used": int(len(obs))})
            curves[shape_name] = rows
        out[sc_name]["curves"] = curves
    return {
        "description": ("Attenuation: an outcome is generated from TRUE tau at a calibrated latent "
                        "R^2, then regressed on the clock's empirical PREDICTED tau."),
        "design": {"cohort_design": "snapshot — one frame per sampled embryo",
                   "n_snapshot_embryos": n_snapshot, "n_replicates": B,
                   "resampling": "embryo-level with replacement",
                   "error_model": "empirical (true, predicted) pairs, not a Gaussian assumption",
                   "scenario_seed_offsets": dict(SCENARIO_SEED_OFFSET),
                   "seed_note": ("scenario seeds are explicit stable offsets from --seed; they are "
                                 "never derived from hash(), which Python randomises per process"),
                   "latent_r2_definition": (
                       "the noise sd is solved per replicate to TARGET the requested OLS R^2 "
                       "against TRUE tau, using the same single-predictor model form later applied "
                       "to predicted tau. This is a target in expectation, not an identity: a "
                       "finite random draw makes the realized R^2 scatter around it (and sample "
                       "R^2 is slightly upward-biased). The REALIZED value is measured per "
                       "replicate and reported as achieved_true_r2_median / achieved_true_r2_ci95, "
                       "and the curve is plotted against the achieved value, not the request.")},
        "scenarios": out,
        "no_unique_inversion": (
            "The two scenarios give DIFFERENT mappings from observed to latent R^2, so there is no "
            "unique inversion. An observed value must be read against the scenario that matches the "
            "cohort it came from, and even then the interval is wide."),
        "interpretation_guard": (
            "This quantifies ONLY the attenuation caused by pseudotime error. A low observed R^2 "
            "is NOT thereby explained: transcript measurement noise, probe-set and batch effects, "
            "biological heterogeneity between embryos, and small n all reduce R^2 independently. "
            "The curve gives an UPPER envelope for what a given latent trend could look like "
            "through this clock, not a claim that the trend exists."),
    }


def invert_curve(curve_rows, observed_r2):
    """Range of ACHIEVED true-time R^2 whose simulated 95% band contains the observation. Reported
    per scenario — the two scenarios disagree, so there is NO unique inversion."""
    rows = [r for r in curve_rows if r.get("achievable")]
    compat = [r["achieved_true_r2_median"] for r in rows
              if r["observed_r2_ci95"][0] <= observed_r2 <= r["observed_r2_ci95"][1]]
    if not compat:
        return None
    return {"min": round(min(compat), 4), "max": round(max(compat), 4)}


# ═══════════════════════════ 3 · frozen downstream illustration ═══════════════════════════
def downstream(fixed, seed, B, atten_curves):
    if not os.path.isfile(GENES_P):
        return {"status": "not estimable", "reason": "pronuclei_genes.json.gz not present"}
    ga = json.load(gzip.open(GENES_P, "rt"))["embryos"]
    counts = {e["id"]: e.get("genes", {}) for e in ga}
    fx = {r["id"]: r for r in fixed["embryos"]}
    rng = np.random.default_rng(seed)

    def fit_gene(gene, include_ood):
        ids, tau, cnt = [], [], []
        for eid, r in fx.items():
            if r.get("tau") is None:
                continue
            if not include_ood and r.get("qc") == "out-of-domain":
                continue
            c = counts.get(eid, {}).get(gene)
            if c is None:
                continue
            ids.append(eid); tau.append(r["tau"]); cnt.append(float(c))
        tau = np.array(tau, float); cnt = np.array(cnt, float)
        if len(tau) < MIN_N or np.std(tau) == 0:
            return {"status": "not estimable",
                    "reason": (f"only {len(tau)} zygotes carry {gene} with a usable tau "
                               f"(need >= {MIN_N})" if len(tau) < MIN_N
                               else "no variation in tau among the zygotes carrying this gene"),
                    "n": int(len(tau))}
        sl, ic, r, pv, se = stats.linregress(tau, cnt)
        # embryo-level permutation: shuffle tau against counts, keeping n fixed
        perm = np.array([r2_of(rng.permutation(tau), cnt) for _ in range(B)])
        p_perm = float((np.sum(perm >= r * r) + 1) / (B + 1))
        # embryo bootstrap CI on R^2
        bs = []
        for _ in range(B):
            k = rng.integers(0, len(tau), len(tau))
            if np.std(tau[k]) == 0:
                continue
            bs.append(r2_of(tau[k], cnt[k]))
        bs = np.array([v for v in bs if np.isfinite(v)])

        # ---- ROBUSTNESS. Transcript counts are heavy-tailed and a single ZGA burst can carry an
        # OLS R^2 on its own. These diagnostics are mandatory, not optional: without them an
        # outlier-driven R^2 would be reported as if it were a stable trend.
        order = np.argsort(cnt)[::-1]
        loo = []
        for k in range(1, min(4, len(cnt) - 3)):
            m = np.ones(len(cnt), bool); m[order[:k]] = False
            if np.std(tau[m]) == 0:
                continue
            s2, _, r2_, p2, _ = stats.linregress(tau[m], cnt[m])
            loo.append({"dropped_top_n_by_count": k, "n": int(m.sum()),
                        "r2": round(float(r2_ * r2_), 5), "p": round(float(p2), 5)})
        infl = []
        for i in range(len(cnt)):
            m = np.ones(len(cnt), bool); m[i] = False
            if np.std(tau[m]) == 0:
                continue
            infl.append(abs(r2_of(tau[m], cnt[m]) - r * r))
        max_infl = float(np.nanmax(infl)) if infl else None
        rho, rho_p = stats.spearmanr(tau, cnt)
        perm_rho = np.array([abs(stats.spearmanr(rng.permutation(tau), cnt).statistic)
                             for _ in range(B)])
        p_perm_rho = float((np.sum(perm_rho >= abs(rho)) + 1) / (B + 1))
        lg = np.log1p(cnt)
        slg, _, rlg, plg, _ = stats.linregress(tau, lg)
        robust = {
            "why": ("transcript counts are heavy-tailed; MERVL in particular fires as a burst at "
                    "ZGA, so one late embryo can dominate an OLS fit"),
            "max_count": int(cnt.max()), "median_count": float(np.median(cnt)),
            "count_ratio_max_to_next": (round(float(np.sort(cnt)[-1] / max(np.sort(cnt)[-2], 1)), 1)
                                        if len(cnt) > 1 else None),
            "leave_top_out": loo,
            "max_single_point_influence_on_r2": (round(max_infl, 5)
                                                 if max_infl is not None else None),
            "spearman_rho": round(float(rho), 5), "spearman_p": round(float(rho_p), 5),
            "spearman_permutation_p": p_perm_rho,
            "log1p_ols_r2": round(float(rlg * rlg), 5), "log1p_ols_p": round(float(plg), 5),
            "log1p_direction": "increases with tau" if slg > 0 else "decreases with tau",
            "verdict": ("OUTLIER-DRIVEN: dropping the single highest-count embryo changes R^2 by "
                        f"{max_infl:.3f}" if (max_infl is not None and max_infl > 0.15)
                        else "stable to single-point removal"),
        }
        return {
            "status": "estimable", "gene": gene, "n": int(len(tau)),
            "include_out_of_domain": include_ood,
            "model_form": "ordinary least squares: count ~ tau (single predictor, no covariates)",
            "observed_r2": round(float(r * r), 5), "pearson_r": round(float(r), 5),
            "slope_counts_per_unit_tau": round(float(sl), 3), "intercept": round(float(ic), 3),
            "direction": "increases with tau" if sl > 0 else "decreases with tau",
            "ols_p": float(pv),
            "permutation_p_embryo_level": p_perm, "n_permutations": B,
            "r2_bootstrap_ci95": ([round(float(np.quantile(bs, .025)), 5),
                                   round(float(np.quantile(bs, .975)), 5)] if len(bs) else None),
            "robustness": robust,
            "probe_set_coverage": {"n_zygotes_with_gene": int(len(tau)),
                                   "n_zygotes_with_tau": int(sum(1 for r_ in fx.values()
                                                                 if r_.get("tau") is not None)),
                                   "note": ("MERFISH panels are disjoint, so a gene is measured "
                                            "only in the subset of zygotes whose panel contains it")},
        }

    results = {}
    for g in PREDECLARED:
        res = fit_gene(g, include_ood=False)
        if res.get("status") == "estimable":
            res["sensitivity_including_out_of_domain"] = {
                k: v for k, v in fit_gene(g, include_ood=True).items()
                if k in ("n", "observed_r2", "pearson_r", "permutation_p_embryo_level")}
            res["compatible_latent_true_r2_by_scenario"] = {
                sc: invert_curve(v["curves"]["linear"], res["observed_r2"])
                for sc, v in atten_curves["scenarios"].items()}
        results[g] = res

    # ---- separate DISCOVERY scan, FDR-corrected, explicitly not confirmatory ----
    scan = []
    for gene in sorted({g for c in counts.values() for g in c}):
        r_ = fit_gene(gene, include_ood=False)
        if r_.get("status") == "estimable":
            scan.append({"gene": gene, "n": r_["n"], "r2": r_["observed_r2"], "p": r_["ols_p"]})
    scan.sort(key=lambda s: s["p"])
    m = len(scan)
    for i, s in enumerate(scan, 1):                       # Benjamini-Hochberg
        s["q_bh"] = round(float(min(1.0, s["p"] * m / i)), 6)
    for i in range(m - 2, -1, -1):
        scan[i]["q_bh"] = round(float(min(scan[i]["q_bh"], scan[i + 1]["q_bh"])), 6)

    conf = results.get(CONFIRMATORY, {})
    return {
        "confirmatory_gene": CONFIRMATORY,
        "predeclared_genes": PREDECLARED,
        "predeclared_note": ("This list is fixed in the source file and was NOT chosen after seeing "
                             "results. Out-of-domain zygotes are EXCLUDED by default; a sensitivity "
                             "fit including them is reported alongside."),
        "results": results,
        "discovery_scan": {
            "status": "DISCOVERY — not confirmatory; reported separately and FDR-corrected",
            "n_genes_tested": m,
            "n_q_below_0_05": int(sum(1 for s in scan if s["q_bh"] < 0.05)),
            "top20": scan[:20],
            "note": ("A screen over every gene is exploratory. Any hit here would need independent "
                     "confirmation; it does not license reinterpreting the predeclared result."),
        },
        "headline": build_headline(conf, atten_curves),
    }


def build_headline(conf, atten):
    if conf.get("status") != "estimable":
        return {"verdict": "not estimable",
                "text": (f"The predeclared confirmatory gene {CONFIRMATORY} could not be evaluated: "
                         f"{conf.get('reason', 'unknown reason')}.")}
    obs = conf["observed_r2"]
    bysc = conf.get("compatible_latent_true_r2_by_scenario", {})
    sc_name = ("fixed_cohort_matched" if "fixed_cohort_matched" in atten["scenarios"]
               else "validation_uniform")
    comp = bysc.get(sc_name)
    rb = conf.get("robustness", {})
    lin = [r for r in atten["scenarios"][sc_name]["curves"]["linear"] if r.get("achievable")]
    def med(x):
        if not lin:
            return float("nan")
        r = min(lin, key=lambda v: abs(v["achieved_true_r2_median"] - x))
        return r["observed_r2_median"]

    fragile = (rb.get("max_single_point_influence_on_r2") or 0) > 0.15
    loo1 = next((d for d in rb.get("leave_top_out", []) if d["dropped_top_n_by_count"] == 1), None)

    txt = (f"{CONFIRMATORY} against calibrated tau gives OLS R^2 = {obs:.3f} "
           f"(n = {conf['n']}, embryo-level permutation p = "
           f"{conf['permutation_p_embryo_level']:.3f}). ")
    if fragile:
        txt += (f"THIS NUMBER IS NOT ROBUST: one embryo carries {rb.get('max_count')} counts "
                f"({rb.get('count_ratio_max_to_next')}x the next highest), and removing that single "
                f"embryo moves R^2 by {rb['max_single_point_influence_on_r2']:.3f}"
                + (f" (to {loo1['r2']:.3f}, p = {loo1['p']:.3f})" if loo1 else "") + ". "
                f"The outlier-robust statistics are Spearman rho = {rb.get('spearman_rho'):.3f} "
                f"(permutation p = {rb.get('spearman_permutation_p'):.3f}) and log1p-OLS "
                f"R^2 = {rb.get('log1p_ols_r2'):.3f}. Biologically the burst is plausible — MERVL "
                f"fires at ZGA and that embryo sits latest on the clock — but a single point cannot "
                f"carry the claim. ")
    txt += (f"For scale, under the measured clock error and the {sc_name.replace('_', '-')} "
            f"scenario, a latent trend of true R^2 = 0.50 is expected to present as observed R^2 "
            f"~ {med(0.5):.3f}, true R^2 = 0.30 as ~ {med(0.3):.3f}, and true R^2 = 0.10 as "
            f"~ {med(0.1):.3f}. ")
    if comp:
        other = {k: v for k, v in bysc.items() if k != sc_name and v}
        txt += (f"Taken at face value the observed R^2 is compatible with latent true R^2 in "
                f"[{comp['min']:.2f}, {comp['max']:.2f}] under the {sc_name.replace('_', '-')} "
                f"scenario; the simulation cannot narrow it further")
        if other:
            o = list(other.items())[0]
            txt += (f", and the {o[0].replace('_', '-')} scenario gives a DIFFERENT range "
                    f"[{o[1]['min']:.2f}, {o[1]['max']:.2f}] — there is no unique inversion.")
        else:
            txt += "."
    else:
        txt += ("The observed value falls outside the simulated band at every latent R^2 on the "
                "grid, indicating variance the clock-error model does not capture.")

    if conf["permutation_p_embryo_level"] >= 0.05:
        verdict = "not-significant"
        txt += (" The embryo-level permutation test is NOT significant, so this is not "
                "distinguishable from chance in this cohort regardless of attenuation.")
    elif fragile:
        verdict = "significant-but-outlier-driven"
    elif comp:
        verdict = "compatible-but-not-demonstrated"
    else:
        verdict = "outside-simulated-band"

    return {"verdict": verdict, "observed_r2": obs,
            "robust_spearman_rho": rb.get("spearman_rho"),
            "robust_log1p_r2": rb.get("log1p_ols_r2"),
            "outlier_driven": bool(fragile),
            "scenario_used": sc_name,
            "answer_to_the_low_r2_question": (
                "Yes — an observed R^2 of 0.10-0.15 is compatible with a real, moderately strong "
                f"latent trend under this clock's measured error. In the {sc_name.replace('_', '-')} "
                f"scenario the simulation puts a latent true R^2 of 0.15 at an observed median of "
                f"{med(0.15):.3f} and a latent 0.30 at {med(0.3):.3f}, with wide intervals at "
                f"n~{atten['design']['n_snapshot_embryos']}. A low observed R^2 is therefore NOT "
                "evidence against a trend. It is equally NOT evidence for one: clock error is only "
                "one of several variance sources, this analysis bounds only that one, and the two "
                "sampling scenarios give different mappings, so no unique inversion exists."),
            "text": txt}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20260723)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    B = 300 if a.quick else 2000
    B_SIM = 200 if a.quick else 1200

    if not os.path.isfile(os.path.join(CAL_DIR, "calibration.json")):
        raise SystemExit("run scripts/train_pronuclear_pseudotime.py first")
    cal, frames = load_frozen()
    ver = cal["meta"]["model_version"]
    print(f"noise-ceiling analysis · frozen clock {ver}")

    rel = reliability(cal, frames, a.seed, B)
    print(f"  reliability (outer-test): r={rel['pearson_r']:.3f} R2={rel['pearson_r2']:.3f} "
          f"rho={rel['spearman']:.3f} macroMAE={rel['macro_mae_tau']:.4f} "
          f"medianAE={rel['median_ae_tau']:.4f} tau ({rel['median_ae_hours']:.2f} h)")

    n_fixed = 0
    fixed = None
    fixed_taus = None
    if os.path.isfile(FIXED_P):
        fixed = json.load(open(FIXED_P))
        keep = [r for r in fixed["embryos"]
                if r.get("tau") is not None and r.get("qc") != "out-of-domain"]
        n_fixed = len(keep)
        fixed_taus = [r["tau"] for r in keep]          # OOD excluded by default, as downstream
    n_snapshot = max(12, n_fixed or 40)
    grid = [0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95]
    att = attenuation(frames, a.seed + 1, B_SIM, n_snapshot, grid, fixed_taus)
    print(f"  attenuation (n={n_snapshot} snapshot embryos, {B_SIM} reps, "
          f"{len(att['scenarios'])} scenarios):")
    for sc, v in att["scenarios"].items():
        print(f"    [{sc}]")
        for r in v["curves"]["linear"]:
            if r["latent_true_r2"] in (0.1, 0.3, 0.5, 0.8):
                if not r.get("achievable"):
                    print(f"      requested {r['latent_true_r2']:.2f} -> NOT ACHIEVABLE")
                else:
                    print(f"      requested {r['latent_true_r2']:.2f} "
                          f"(achieved {r['achieved_true_r2_median']:.3f}) -> observed median "
                          f"{r['observed_r2_median']:.3f} CI {r['observed_r2_ci95']}")
        sat = v["curves"]["monotone_saturating"]
        n_un = sum(1 for r in sat if not r.get("achievable"))
        if n_un:
            print(f"      (saturating trend: {n_un}/{len(sat)} grid points unreachable — rho^2 cap)")

    down = {"status": "not estimable", "reason": "data/pronuclei_pseudotime.json not present"}
    if fixed:
        if fixed["meta"]["model_version"] != ver:
            down = {"status": "stale", "reason": (
                f"fixed predictions were made with {fixed['meta']['model_version']} but the frozen "
                f"clock is {ver}; re-run build_pronuclei_pseudotime.py")}
            print(f"  !! fixed cohort is stale ({fixed['meta']['model_version']} != {ver})")
        else:
            down = downstream(fixed, a.seed + 2, B, att)
            h = down["headline"]
            print(f"  downstream {CONFIRMATORY}: {h.get('verdict')}")
            print(f"    {h.get('text','')[:200]}")

    payload = {"meta": {
        "model_version": ver, "data_version": cal["meta"]["data_version"],
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "seed": a.seed, "quick": bool(a.quick),
        "purpose": ("quantify how much the measured pseudotime error attenuates transcript-vs-"
                    "pseudotime R^2, so an observed value can be interpreted rather than dismissed"),
        "transcripts_used_in_clock": False,
        "transcripts_used_here": "only in the frozen downstream illustration (part 3)"},
        "clock_reliability": rel, "attenuation": att, "downstream_illustration": down}
    out = os.path.join(CAL_DIR, "noise_ceiling.json")
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"  wrote {os.path.relpath(out, HERE)} ({os.path.getsize(out) / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
