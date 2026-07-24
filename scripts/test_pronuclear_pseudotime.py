#!/usr/bin/env python3
"""
Tests for the pronuclear pseudotime calibration (pnpt-3.x).

Designed to FAIL on the specific ways this analysis could quietly become wrong:
  * an embryo crossing a train / conformal-calibration / test boundary
  * the same embryos both setting an interval and testing its coverage
  * model selection seeing outer-test outcomes
  * transcript / gene / count fields entering any clock-training path
  * reported pairwise accuracy disagreeing with a hand-checked example containing ties
  * this repository being unable to regenerate artifacts without the sibling embryo_viewer checkout
  * an out-of-domain fixed embryo entering downstream regressions by default
  * reliability / noise-ceiling artifacts not matching the frozen model version

Run:  python3 scripts/test_pronuclear_pseudotime.py
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))
import train_pronuclear_pseudotime as T  # noqa: E402

ART = os.path.join(HERE, "data", "pseudotime_calibration")
FIXED = os.path.join(HERE, "data", "pronuclei_pseudotime.json")
FAILED, PASSED = [], []


def check(name, cond, detail=""):
    (PASSED if cond else FAILED).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  — {detail}" if detail and not cond else ""))


def load(p, gz=False):
    return json.load(gzip.open(p, "rt")) if gz else json.load(open(p))


def main():
    print("pronuclear pseudotime calibration — tests\n")

    # ───────────────────────────── source schema ─────────────────────────────
    print("[source validation]")
    df, summary = T.load_and_validate()
    check("53 embryos", summary["n_embryos"] == 53)
    check("2057 frames", summary["n_frames"] == 2057)
    check("tau within [0,1]", df[T.TARGET].between(0, 1).all())
    check("tau == time_h / duration",
          np.allclose(df["time_h"] / df["migration_duration_h"], df[T.TARGET], atol=1e-9))
    check("time monotonic within every embryo",
          all(g["time_h"].is_monotonic_increasing for _, g in df.groupby(T.GROUP)))
    with open(T.SRC_CSV, "rb") as fh:
        check("sha256 matches the file on disk",
              hashlib.sha256(fh.read()).hexdigest() == summary["sha256"])
    tmp = "/tmp/pnpt_bad.csv"
    try:
        df.iloc[:-1].to_csv(tmp, index=False)
        rejected = False
        try:
            T.load_and_validate(tmp)
        except SystemExit:
            rejected = True
        check("validator rejects a wrong row count", rejected)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    art = load(os.path.join(ART, "calibration.json"))
    ver = art["meta"]["model_version"]

    # ───────────────── hand-checked ordering with TIES (the v1 bug) ─────────────────
    print("\n[ordering metrics — hand-checked]")
    # truth 1<2<3<4 ; prediction ties the first two, orders the rest correctly.
    # Pairs: (1,2) pred tie ; (1,3)(1,4)(2,3)(2,4)(3,4) all concordant  -> 5 concordant, 1 tied.
    yt = [1.0, 2.0, 3.0, 4.0]
    yp = [5.0, 5.0, 6.0, 7.0]
    oc = T.ordering_counts(yt, yp)
    check("hand example: 6 comparable pairs", oc["comparable"] == 6, str(oc))
    check("hand example: 5 concordant", oc["concordant"] == 5, str(oc))
    check("hand example: 1 tied prediction", oc["tied_pred"] == 1, str(oc))
    check("hand example: 0 discordant", oc["discordant"] == 0, str(oc))
    check("hand example: strict accuracy == 5/6", abs(oc["strict_accuracy"] - 5 / 6) < 1e-12,
          f"got {oc['strict_accuracy']}")
    check("hand example: half-credit == 5.5/6", abs(oc["half_credit_accuracy"] - 5.5 / 6) < 1e-12)
    check("hand example: strict + tie + inversion == 1",
          abs(oc["strict_accuracy"] + oc["tie_rate"] + oc["inversion_rate"] - 1) < 1e-12)
    # tau-b would report a DIFFERENT number here; the point of v2 is that we no longer use it
    from scipy.stats import kendalltau
    taub_acc = (kendalltau(yt, yp).statistic + 1) / 2
    check("strict accuracy differs from the old (tau_b+1)/2 formula under ties",
          abs(oc["strict_accuracy"] - taub_acc) > 1e-6,
          f"strict={oc['strict_accuracy']:.4f} taub_formula={taub_acc:.4f}")
    # tied TRUTH pairs must be excluded from the denominator
    oc2 = T.ordering_counts([1.0, 1.0, 2.0], [3.0, 4.0, 5.0])
    check("tied-truth pairs excluded from comparable", oc2["comparable"] == 2 and oc2["tied_true"] == 1,
          str(oc2))
    # a perfectly reversed prediction is 0% strict, 100% inversions
    oc3 = T.ordering_counts([1.0, 2.0, 3.0], [3.0, 2.0, 1.0])
    check("reversed prediction => strict 0.0, inversion 1.0",
          oc3["strict_accuracy"] == 0.0 and oc3["inversion_rate"] == 1.0, str(oc3))

    # ───────────────────────── nested CV: no leakage anywhere ─────────────────────────
    print("\n[nested CV leakage]")
    fold_of = art["folds"]["outer_assignment"]
    check("every embryo has exactly one outer fold", len(fold_of) == 53)
    frames = load(os.path.join(ART, "oof_frames.json.gz"), gz=True)["frames"]
    by = {}
    for f in frames:
        by.setdefault(f["embryo_id"], set()).add(f["outer_fold"])
    check("no embryo ID appears in more than one outer fold",
          all(len(v) == 1 for v in by.values()),
          str([k for k, v in by.items() if len(v) > 1][:5]))
    check("frame outer-fold matches the assignment table",
          all(next(iter(v)) == fold_of[k] for k, v in by.items()))
    ok = True
    for rec in art["nested_evaluation"]["per_outer_fold"]:
        te = set(rec["test_embryos"])
        if rec["n_train_embryos"] + rec["n_test_embryos"] != 53:
            ok = False
        # the fold's test embryos must all carry that fold id
        if any(fold_of[e] != rec["outer_fold"] for e in te):
            ok = False
    check("each outer fold partitions all 53 embryos, test set matches its fold id", ok)
    check("inner folds exist and are reported per outer fold",
          all("inner_scores" in r and len(r["inner_scores"]) >= 4
              for r in art["nested_evaluation"]["per_outer_fold"]))
    check("model choice is recorded per outer fold (made by inner CV)",
          all(r.get("chosen_by_inner_cv") for r in art["nested_evaluation"]["per_outer_fold"]))

    # selection must not have seen outer-test outcomes: the inner score keys must be exactly the
    # candidate keys, and no per-fold record may carry an outer-test metric inside inner_scores
    leaked = []
    for r in art["nested_evaluation"]["per_outer_fold"]:
        for k, v in r["inner_scores"].items():
            if set(v.keys()) - {"macro_mae", "strict_ordering"}:
                leaked.append((r["outer_fold"], k, sorted(v.keys())))
    check("inner selection scores carry only inner-CV quantities (no outer-test metrics)",
          not leaked, str(leaked[:3]))

    # the production lock must be derivable from the INNER aggregate alone
    agg = art["production"]["inner_cv_aggregate"]
    recomputed = T.pick({k: {"macro_mae": v["macro_mae"], "strict_ordering": v["strict_ordering"]}
                         for k, v in agg.items()})
    check("production model is reproducible from the inner-CV aggregate alone",
          recomputed == art["production"]["key"],
          f"recomputed={recomputed} stored={art['production']['key']}")

    # ───────────── interval: honest naming, no false conformal claim, code==metadata ─────────────
    print("\n[prediction interval — naming and validity]")
    u = art["uncertainty"]
    sp = u["split_per_replicate"]
    check("fit / calibration / test embryo counts are disjoint and sum to 53",
          sp["fit_embryos"] + sp["calibration_embryos"] + sp["test_embryos"] == 53, str(sp))
    check("all three roles are non-empty", min(sp.values()) > 0, str(sp))
    check("primary interval is named an EMPIRICAL disjoint-embryo interval",
          "empirical" in u["interval_name"].lower() and "disjoint-embryo" in u["interval_name"].lower(),
          u.get("interval_name"))
    check("primary interval explicitly disclaims a formal guarantee",
          u["formal_guarantee"] is False)

    # THE v2 BUG: n_cal below the cluster-conformal threshold while claiming a finite 95% certificate
    ALPHA = 0.05
    N_MIN = int(np.ceil(1 / ALPHA)) - 1                       # 19
    ncal = sp["calibration_embryos"]
    txt_all = json.dumps(u)
    claims_cert = bool(re.search(r"certificate is attainable|guarantee(?:d|s)? 95|"
                                 r"exact 95% (?:cluster|conformal)", txt_all, re.I))
    check(f"primary interval (n_cal={ncal}) makes no exact-95% cluster-conformal claim",
          not (ncal < N_MIN and claims_cert),
          f"n_cal={ncal} < {N_MIN} yet the text claims a finite exact 95% certificate")
    # the word may appear in a negation ("NO conformal correction is applied"); what must not
    # happen is the method DESCRIBING ITSELF as conformal
    m_low = u["method"].lower()
    self_describes = bool(re.search(r"(?<!no )(?<!not )(?:split-)?conformal (?:interval|quantile|"
                                    r"construction|certificate)|is a (?:split-)?conformal", m_low))
    check("primary method does not describe ITSELF as conformal",
          not self_describes, u["method"][:140])
    check("primary method explicitly negates the conformal correction",
          bool(re.search(r"no .{0,40}conformal correction", m_low)))
    check("primary states why a conformal correction would be invalid",
          "exchangeab" in u["why_not_conformal"].lower())

    # code/metadata agreement: v2 used the FRAME count while saying EMBRYO count
    src = open(os.path.join(HERE, "scripts", "train_pronuclear_pseudotime.py")).read()
    prim = re.search(r"hw = float\(np\.quantile\(res_cal,([^)]*)\)\)", src)
    check("primary half-width is a plain percentile (no hidden finite-sample correction)",
          bool(prim) and "ceil" not in prim.group(1), prim.group(1) if prim else "not found")
    says_embryo_count = "correction taken over" in u["method"] and "EMBRYO count" in u["method"]
    check("metadata does not claim an embryo-count correction the code never applies",
          not says_embryo_count)
    check("metadata states exactly what the quantile is computed over",
          "FRAME" in u["quantile_computed_over"].upper()
          or "frame" in u["quantile_computed_over"].lower(), u.get("quantile_computed_over"))
    check("estimator note admits averaging widths is not a conformal construction",
          "not" in u["estimator_note"].lower() and "conformal" in u["estimator_note"].lower())
    check("fixed-cohort application is labelled empirical transfer uncertainty",
          "EMPIRICAL TRANSFER UNCERTAINTY" in u["applies_to_fixed_cohort_as"])

    # ADVERSARIAL: these two guards must FIRE on the exact v2 defects, so prove they do.
    print("\n[guards fire on the v2 defects]")
    def guard_ncal_vs_claim(n_cal_, text_):
        """True = guard fires (a false exact-95% cluster claim with too few calibration embryos)."""
        claims = bool(re.search(r"certificate is attainable|guarantee(?:d|s)? 95|"
                                r"exact 95% (?:cluster|conformal)", text_, re.I))
        return n_cal_ < N_MIN and claims
    v2_text = ("A conformal certificate that is exact under EMBRYO exchangeability needs "
               "ceil((n+1)*0.95) <= n, i.e. n >= 19 calibration embryos. This design uses 15 "
               "calibration embryos per replicate, so the embryo-level certificate is attainable "
               "but tight;")
    check("guard FIRES on the v2 text (15 cal embryos claiming an attainable certificate)",
          guard_ncal_vs_claim(15, v2_text))
    check("guard does NOT fire on the current v3 artifact",
          not guard_ncal_vs_claim(ncal, txt_all))

    def guard_count_mismatch(code_expr, meta_text):
        """True = guard fires (code uses the FRAME count while metadata claims EMBRYO count)."""
        uses_frame = "n_f" in code_expr
        says_embryo = "EMBRYO count" in meta_text
        return uses_frame and says_embryo
    check("guard FIRES on the v2 frame-count-vs-embryo-count mismatch",
          guard_count_mismatch("np.ceil((n_f + 1) * 0.95) / n_f",
                               "the finite-sample correction taken over the EMBRYO count"))
    prim_expr = prim.group(1) if prim else ""
    check("guard does NOT fire on the current v3 code/metadata pair",
          not guard_count_mismatch(prim_expr, u["method"]))

    print("\n[cluster-conformal construction sensitivity]")
    rc = u["rigorous_cluster_sensitivity"]
    ns = rc["nonconformity_score"].lower()
    check("uses ONE nonconformity score per calibration embryo (an embryo-level aggregate)",
          "embryo" in ns and any(k in ns for k in ("max", "maximum", "quantile", "per-embryo")),
          rc["nonconformity_score"])
    check("the score is an aggregate over the embryo's frames, not a per-frame score",
          "frames of" in ns or "over all frames" in ns, rc["nonconformity_score"])
    check(f"n_cal >= {N_MIN} so the conformal quantile index is attainable",
          rc["n_calibration_embryos"] >= N_MIN, str(rc["n_calibration_embryos"]))
    k_expect = int(np.ceil((rc["n_calibration_embryos"] + 1) * (1 - ALPHA)))
    check("quantile index k matches ceil((n_cal+1)*(1-alpha))",
          rc["quantile_index_k"] == k_expect, f"{rc['quantile_index_k']} != {k_expect}")
    check("k <= n_cal (each construction is finite and exact)",
          rc["quantile_index_k"] <= rc["n_calibration_embryos"])
    check("states precisely what each construction protects (simultaneous per-embryo)",
          "SIMULTANEOUS" in rc["protects"].upper())
    check("reports simultaneous per-embryo coverage, not just frame coverage",
          rc.get("coverage_simultaneous_per_embryo") is not None)
    check("discloses whether the maximum score is used",
          isinstance(rc.get("uses_maximum_score"), bool))
    if rc.get("uses_maximum_score"):
        check("explains why the maximum is forced at this sample size",
              "MAXIMUM" in rc["feasibility_note"].upper())
    check("rigorous construction is wider than the empirical interval (honest cost)",
          rc["halfwidth_mean"] >= u["halfwidth_mean"],
          f"rigorous {rc['halfwidth_mean']} vs empirical {u['halfwidth_mean']}")
    check("no quantile-of-quantiles is described as simultaneous",
          "quantile of quantiles" not in json.dumps(rc).lower())

    # ── THE AGGREGATE OVERCLAIM: a guarantee must not attach to an AVERAGED width ──
    print("\n[aggregate width carries no guarantee]")
    check("per-split construction guarantee is stated separately from the reported mean",
          "per_split_construction_has_finite_sample_guarantee" in rc
          and "reported_mean_halfwidth_has_formal_guarantee" in rc,
          str(sorted(k for k in rc if "guarantee" in k)))
    check("the REPORTED MEAN half-width does NOT claim a formal guarantee",
          rc["reported_mean_halfwidth_has_formal_guarantee"] is False)
    check("each individual construction DOES carry the finite-sample guarantee",
          rc["per_split_construction_has_finite_sample_guarantee"] is True)
    check("no bare `formal_guarantee: true` remains on the aggregate block",
          rc.get("formal_guarantee") is not True, str(rc.get("formal_guarantee")))
    check("guarantee scope says the mean is a summary, not a conformal interval",
          "NOT a conformal interval" in rc["guarantee_scope"]
          or "not a conformal interval" in rc["guarantee_scope"].lower())
    check("assumptions list exchangeability and calibration-independent fitting",
          any("exchangeab" in a.lower() for a in rc["assumptions"])
          and any("independent" in a.lower() for a in rc["assumptions"]))
    check("states the guarantee does NOT transfer to the all-embryo production model",
          any("all 53 embryos" in a or "refitted on all" in a for a in rc["assumptions"])
          or "refitted on all 53" in json.dumps(u))
    check("the mean is labelled as a mean over multiple constructions",
          "MEAN" in rc["what_is_reported"].upper() and rc.get("n_constructions", 0) > 1)
    check("per-split width range is reported so the spread behind the mean is visible",
          isinstance(rc.get("halfwidth_per_split_range"), list))
    # the UI must not present the averaged width as a rigorous guaranteed interval
    js = open(os.path.join(HERE, "pseudotime-calibration.js")).read()
    check("UI panel is not titled as a rigorous guaranteed interval",
          "Sensitivity — rigorous cluster-level conformal" not in js)
    check("UI states the guarantee is per construction, not for the displayed number",
          "guarantee is per construction" in js)

    print("\n[interval coverage reporting]")
    check("frame-marginal and embryo-macro coverage reported separately",
          u["coverage_frame_marginal"] is not None and u["coverage_embryo_macro"] is not None)
    check("coverage carries a bootstrap interval",
          isinstance(u["coverage_frame_marginal"].get("ci95"), list))
    check("coverage reported by tau phase", len(u.get("coverage_by_phase", {})) >= 2)
    check("multiple replicates were run", u["n_replicates"] >= 5)

    # ───────── attenuation seeds must be stable across PROCESSES, not just within one ─────────
    print("\n[attenuation seed reproducibility]")
    A_PATH = os.path.join(HERE, "scripts", "analyze_pseudotime_noise_ceiling.py")
    asrc = open(A_PATH).read()
    # AST, not substring: the source deliberately WARNS about hash() in comments and strings, so a
    # text search false-positives. Only an actual call to the builtin `hash` matters.
    hash_calls = [n for n in ast.walk(ast.parse(asrc))
                  if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "hash"]
    check("scenario seeds are NOT derived from hash() (Python randomises it per process)",
          not hash_calls, f"{len(hash_calls)} literal hash() call(s) at lines "
                          f"{[n.lineno for n in hash_calls]}")
    # and prove the guard would fire on the v3 defect it replaces
    bad_src = "rng = np.random.default_rng(seed + abs(hash(sc_name)) % 10000)"
    fires = [n for n in ast.walk(ast.parse(bad_src))
             if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "hash"]
    check("guard FIRES on the previous hash-derived seed expression", len(fires) == 1)
    check("explicit stable scenario seed offsets are declared",
          "SCENARIO_SEED_OFFSET" in asrc)
    # run the seed computation in FRESH processes under different PYTHONHASHSEED values
    prog = ("import sys,json; sys.path.insert(0,%r); "
            "import analyze_pseudotime_noise_ceiling as A; "
            "print(json.dumps({k: A._scenario_seed(20260723, k) "
            "for k in sorted(A.SCENARIO_SEED_OFFSET)}))" % os.path.join(HERE, "scripts"))
    seeds = []
    for hs in ("0", "1", "12345", "random"):
        env = dict(os.environ, PYTHONHASHSEED=hs)
        r = subprocess.run([sys.executable, "-c", prog], capture_output=True, text=True, env=env)
        seeds.append(r.stdout.strip())
    check("scenario seeds identical across PYTHONHASHSEED=0/1/12345/random",
          len(set(seeds)) == 1, str(seeds))
    # an unknown scenario must raise rather than silently fall back to something unstable
    prog2 = ("import sys; sys.path.insert(0,%r); "
             "import analyze_pseudotime_noise_ceiling as A; "
             "\ntry:\n A._scenario_seed(1,'nope'); print('NO_RAISE')\n"
             "except KeyError: print('RAISED')" % os.path.join(HERE, "scripts"))
    r2 = subprocess.run([sys.executable, "-c", prog2], capture_output=True, text=True)
    check("an unregistered scenario raises instead of using an unstable seed",
          "RAISED" in r2.stdout, r2.stdout.strip() + r2.stderr[-120:])

    nc_probe = os.path.join(ART, "noise_ceiling.json")
    if os.path.isfile(nc_probe):
        ncp = load(nc_probe)
        off = ncp["attenuation"]["design"].get("scenario_seed_offsets")
        check("artifact records the scenario seed offsets actually used", bool(off), str(off))
        check("artifact notes seeds are not hash-derived",
              "hash()" in ncp["attenuation"]["design"].get("seed_note", ""))

    # ────────── nested = SELECTION PROCEDURE, distinct from the locked production family ──────────
    print("\n[what the nested estimate estimates]")
    ne = art["nested_evaluation"]
    check("nested block states it estimates the SELECTION PROCEDURE",
          "PROCEDURE" in ne["estimates"].upper())
    stab = ne.get("inner_selection_stability", {})
    check("inner selection stability is recorded", bool(stab))
    if len(stab) > 1:
        check("nested estimate is NOT described as an unbiased score of the locked family",
              "locked production family" in ne["estimates"] or "not of the" in ne["estimates"])
    ff = art.get("fixed_family_sensitivity")
    check("a fixed-family sensitivity is reported for the locked family", ff is not None)
    if ff:
        check("fixed-family result is labelled POST-SELECTION / descriptive",
              "POST-SELECTION" in ff["status"].upper())
        check("fixed-family explains why it is biased",
              "optimistic" in ff["why"].lower() or "biased" in ff["why"].lower())
        check("fixed-family targets the locked production family",
              ff["key"] == art["production"]["key"])
        check("fixed-family reports its own bootstrap CI",
              isinstance(ff.get("bootstrap_ci95_macro_mae"), list))
        check("the two estimates are reported separately, not merged",
              ff["metrics"]["macro_mae"] != ne["outer_test_metrics"]["macro_mae"]
              or ff["difference_vs_nested_macro_mae"] == 0)

    # ───────────────────── feature governance / no transcripts in the clock ─────────────────────
    print("\n[feature governance]")
    prod_feats = art["production"]["features"]
    banned = ("volume", "male_", "female_", "transcript", "gene", "probe", "expression")
    check("production features contain no banned token",
          not any(b in f.lower() for f in prod_feats for b in banned), str(prod_feats))
    check("no leaky model is marked for production",
          all(not m.get("selected_for_production") for m in art["models"] if not m["deployable"]))
    check("the leaky upper bound is present but flagged non-deployable",
          any(not m["deployable"] for m in art["models"]))
    for bad in (["male_relative_volume"], ["nearer_to_center_um", "volume_sum"]):
        fired = False
        try:
            T.assert_deployable(bad, "unit-test")
        except ValueError:
            fired = True
        check(f"assert_deployable rejects {bad}", fired)

    # AST-level: the clock-training path must not evaluate any transcript identifier
    for src_name in ("scripts/train_pronuclear_pseudotime.py", "build_pronuclei_pseudotime.py"):
        p = os.path.join(HERE, src_name)
        tree = ast.parse(open(p).read())
        names = set()
        for n in ast.walk(tree):
            if isinstance(n, ast.Name):
                names.add(n.id)
            elif isinstance(n, ast.Attribute):
                names.add(n.attr)
        hits = sorted(x for x in names
                      if re.search(r"transcript|gene|probe|expression", x, re.I))
        check(f"{src_name} evaluates no transcript/gene identifier", not hits, str(hits))

    # ───────────────────── self-containment (no embryo_viewer dependency) ─────────────────────
    print("\n[repository self-containment]")
    for f in ("scripts/train_pronuclear_pseudotime.py", "scripts/test_pronuclear_pseudotime.py",
              "scripts/analyze_pseudotime_noise_ceiling.py", "build_pronuclei_pseudotime.py",
              "calibration_data/scheffler2021/scheffler_2021_control_zygote_trajectories.csv",
              "calibration_data/scheffler2021/extract_source_data.py",
              "calibration_data/scheffler2021/README.md",
              "calibration_data/fixed_cohort_geometry.csv"):
        check(f"present in this repo: {f}", os.path.isfile(os.path.join(HERE, f)))
    def code_strings(path):
        """Every string CONSTANT actually evaluated, excluding docstrings — so a prose mention of
        the old layout in a comment/docstring does not mask a real path dependency, and does not
        cause a false failure either."""
        tree = ast.parse(open(path).read())
        docs = set()
        for n in ast.walk(tree):
            if isinstance(n, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                d = ast.get_docstring(n, clean=False)
                if d:
                    docs.add(d)
        return [n.value for n in ast.walk(tree)
                if isinstance(n, ast.Constant) and isinstance(n.value, str) and n.value not in docs]

    for f in ("scripts/train_pronuclear_pseudotime.py", "build_pronuclei_pseudotime.py",
              "scripts/analyze_pseudotime_noise_ceiling.py"):
        p = os.path.join(HERE, f)
        strs = code_strings(p)
        check(f"{f}: no embryo_viewer path in executable code",
              not any("embryo_viewer" in s for s in strs),
              str([s for s in strs if "embryo_viewer" in s][:2]))
        check(f"{f}: no iCloud absolute path in executable code",
              not any("Mobile Documents" in s for s in strs))
        check(f"{f}: no '..' parent escape in executable path constants",
              not any(s.strip() == ".." for s in strs))
    # the default (no --extract) apply path must not import build_pronuclei at module scope
    bsrc = open(os.path.join(HERE, "build_pronuclei_pseudotime.py")).read()
    top = [n for n in ast.walk(ast.parse(bsrc))
           if isinstance(n, ast.Import) and any("build_pronuclei" in a.name for a in n.names)]
    top_level = [n for n in ast.parse(bsrc).body if isinstance(n, ast.Import)]
    check("build_pronuclei is imported lazily (only under --extract)",
          not any("build_pronuclei" in a.name for n in top_level for a in n.names))

    # ───────────────────── artifacts + version consistency ─────────────────────
    print("\n[artifact integrity + versioning]")
    mver = re.match(r"pnpt-(\d+)\.", ver)
    check("model version bumped past v1 (interval + attenuation semantics changed)",
          bool(mver) and int(mver.group(1)) >= 3, ver)
    model = load(os.path.join(ART, "model.json"))
    check("model.json version matches calibration.json", model["model_version"] == ver)
    check("model.json features match the production features", model["features"] == prod_feats)
    check("calibration.json documents the v2->v3 semantic change",
          len(art["meta"].get("changes_from_v2", [])) >= 3)
    check("calibration.json still documents the v1->v2 change",
          len(art["meta"].get("changes_from_v1", [])) >= 3)
    check("model is described as an empirical calibration, not mechanistic",
          "empirical" in art["meta"]["model_class"].lower()
          and "not a mechanistic" in art["meta"]["model_class"].lower())
    check("nested outer-test metrics carry embryo-bootstrap CIs",
          all(art["nested_evaluation"]["bootstrap_ci95_by_embryo"].get(k)
              for k in ("macro_mae", "pearson_r", "pooled_spearman")))
    check("outer-test ordering reports explicit pair counts",
          set(("concordant", "discordant", "tied_pred", "comparable")).issubset(
              art["nested_evaluation"]["outer_test_metrics"]["pooled_ordering"].keys()))

    nc_p = os.path.join(ART, "noise_ceiling.json")
    if os.path.isfile(nc_p):
        nc = load(nc_p)
        check("noise-ceiling artifact matches the frozen model version",
              nc["meta"]["model_version"] == ver,
              f"{nc['meta']['model_version']} != {ver}")
        check("noise-ceiling reliability uses outer-test predictions only",
              "outer-test" in nc["clock_reliability"]["source"])
        check("noise-ceiling states transcripts were not used in the clock",
              nc["meta"]["transcripts_used_in_clock"] is False)
        att = nc["attenuation"]
        check("attenuation reports named SCENARIOS, not a single curve",
              "scenarios" in att and len(att["scenarios"]) >= 1)
        check("validation-uniform scenario present", "validation_uniform" in att["scenarios"])
        if "fixed_cohort_matched" in att["scenarios"]:
            fm = att["scenarios"]["fixed_cohort_matched"]
            check("fixed-cohort-matched scenario documents that only the PREDICTED tau is matched",
                  "TRUE tau" in fm.get("caveat", "") or "true tau" in fm.get("caveat", "").lower())
            check("fixed-cohort-matched scenario says OOD zygotes are excluded",
                  "out-of-domain" in fm["description"].lower())
        check("artifact denies a unique observed->latent inversion",
              "no unique" in att["no_unique_inversion"].lower()
              or "NO unique" in att["no_unique_inversion"])
        for sc, v in att["scenarios"].items():
            check(f"{sc}: covers both linear and a nonlinear trend", len(v["curves"]) >= 2)
        check("attenuation carries an interpretation guard against over-claiming",
              "NOT thereby explained" in att["interpretation_guard"])
        d_ = att["design"]["latent_r2_definition"]
        check("latent R^2 definition says TARGET, not an identity",
              "TARGET" in d_.upper() and "not an identity" in d_.lower(), d_[:140])
        check("latent R^2 definition says the realized value is measured and reported",
              "measured" in d_.lower() and "achieved_true_r2_median" in d_)

        # ACHIEVED vs REQUESTED latent R^2 must agree (the v2 nonlinear mislabelling)
        print("\n[latent R^2 calibration]")
        for sc, v in att["scenarios"].items():
            for shape, rows in v["curves"].items():
                ok_rows = [r for r in rows if r.get("achievable")]
                check(f"{sc}/{shape}: achieved true R^2 recorded for every achievable point",
                      all("achieved_true_r2_median" in r for r in ok_rows))
                worst = max((abs(r["achieved_true_r2_median"] - r["latent_true_r2"])
                             for r in ok_rows), default=0.0)
                check(f"{sc}/{shape}: median achieved true R^2 within 0.03 of requested",
                      worst <= 0.03, f"max deviation {worst:.4f}")
                unreach = [r for r in rows if not r.get("achievable")]
                if unreach:
                    check(f"{sc}/{shape}: unreachable grid points carry a reason",
                          all(r.get("reason") for r in unreach))
        d = nc.get("downstream_illustration", {})
        if d.get("results"):
            check("downstream genes are predeclared in source",
                  d["predeclared_genes"] == T_pre())
            check("discovery scan is separated and FDR-corrected",
                  "DISCOVERY" in d["discovery_scan"]["status"]
                  and all("q_bh" in s for s in d["discovery_scan"]["top20"]))
            conf = d["results"].get(d["confirmatory_gene"], {})
            if conf.get("status") == "estimable":
                check("confirmatory result reports robustness diagnostics",
                      "robustness" in conf and "leave_top_out" in conf["robustness"])
                check("confirmatory result reports n, model form and permutation p",
                      all(k in conf for k in ("n", "model_form", "permutation_p_embryo_level")))
                check("headline states whether the result is outlier-driven",
                      "outlier_driven" in d["headline"])
                check("headline answers the low-R^2 question explicitly",
                      "answer_to_the_low_r2_question" in d["headline"])
    else:
        print("  (noise_ceiling.json not built — run scripts/analyze_pseudotime_noise_ceiling.py)")

    # ───────────────────── fixed cohort + downstream default ─────────────────────
    if os.path.isfile(FIXED):
        print("\n[fixed cohort]")
        fx = load(FIXED)
        rows = fx["embryos"]
        check("fixed predictions match the frozen model version",
              all(r.get("model_version") == ver for r in rows),
              f"artifact={ver} fixed={rows[0].get('model_version') if rows else '?'}")
        check("fixed predictions record the feature schema", fx["meta"]["features"] == prod_feats)
        check("every fixed row has a QC status",
              all(r.get("qc") in ("pass", "caution", "out-of-domain") for r in rows))
        okr = [r for r in rows if r.get("tau") is not None]
        check("every published fixed tau within [0,1]", all(0 <= r["tau"] <= 1 for r in okr))
        check("intervals present and ordered", all(r["lo95"] <= r["tau"] <= r["hi95"] for r in okr))
        check("unavailable predictions carry a reason",
              all(r.get("reason") for r in rows if r.get("tau") is None))
        check("fixed rows carry no transcript/gene field",
              not any(re.search(r"transcript|gene|probe|expression", k, re.I)
                      for r in rows for k in r))
        check("domain shift is reported", "domain_shift" in fx["meta"])
        check("downstream default excludes out-of-domain embryos",
              "EXCLUDED" in fx["meta"]["downstream_default"])
        check("legacy note states the calibration does not validate the surface gap",
              "neither validates nor calibrates" in fx["meta"]["legacy_note"])

        # the downstream analysis must actually honour that default
        if os.path.isfile(nc_p):
            nc = load(nc_p)
            d = nc.get("downstream_illustration", {})
            conf = (d.get("results") or {}).get(d.get("confirmatory_gene", ""), {})
            if conf.get("status") == "estimable":
                ood = {r["id"] for r in rows if r.get("qc") == "out-of-domain"}
                n_ood_with_tau = len([r for r in rows
                                      if r.get("qc") == "out-of-domain" and r.get("tau") is not None])
                check("confirmatory fit ran with include_out_of_domain = False",
                      conf.get("include_out_of_domain") is False)
                check("a sensitivity fit including out-of-domain embryos is also reported",
                      "sensitivity_including_out_of_domain" in conf)
                if n_ood_with_tau:
                    s = conf["sensitivity_including_out_of_domain"]
                    check("including out-of-domain embryos changes n (so the default really excludes them)",
                          s["n"] >= conf["n"], f"default n={conf['n']} incl-OOD n={s['n']}")

    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"  FAILED: {f}")
    return 1 if FAILED else 0


def T_pre():
    import analyze_pseudotime_noise_ceiling as A
    return A.PREDECLARED


if __name__ == "__main__":
    sys.exit(main())
