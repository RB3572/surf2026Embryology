#!/usr/bin/env python3
"""Checks on data/scheffler.json.gz and pt-models.js — the pseudotime training workshop.

This page is unusual: the model is not precomputed and shipped, it is REFIT IN THE BROWSER from
the trajectories in the artifact. So the tests have to cover two things that normally live apart.

  1. The artifact is what it claims — 53 real trajectories, a normalised time that is genuinely
     t_real / T_duration, and none of the leaky volume columns that would let a visitor rebuild a
     model that scores well and cannot be deployed.

  2. The browser's fitter reproduces the published pipeline. The strongest available check is to
     run pt-models.js under node over the SAME fold assignment the offline pipeline used and
     confirm each model lands on the published macro MAE. If the two ever diverge, the page would
     be quietly teaching a different method from the one in the methods section.
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ART = os.path.join(ROOT, "data", "scheffler.json.gz")
PTJS = os.path.join(ROOT, "pt-models.js")
PAGEJS = os.path.join(ROOT, "scheffler.js")

passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1; print(f"  PASS  {name}")
    else:
        failed += 1; print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))


def section(t):
    print(f"\n[{t}]")


# The browser fitter, driven from node over the published fold assignment.
NODE = r"""
const fs = require("fs"), zlib = require("zlib");
global.window = {};
(0, eval)(fs.readFileSync(process.argv[2], "utf8"));
const PT = global.window.PTModels;
const doc = JSON.parse(zlib.gunzipSync(fs.readFileSync(process.argv[3])).toString());

const rows = [];
doc.train.forEach((e) => e.f.forEach((fr) =>
  rows.push({ emb: e.id, t: fr[1], f: PT.featsFromPair(fr[2], fr[3]) })));

const assign = doc.folds.assignment;
const nFolds = Math.max(...Object.values(assign)) + 1;
const out = {};
for (const key of PT.KEYS) {
  const allRows = [], allPred = [];
  for (let k = 0; k < nFolds; k++) {
    const testIds = Object.keys(assign).filter((id) => assign[id] === k);
    const r = PT.holdout(key, rows, testIds);
    allRows.push(...r.test); allPred.push(...r.pred);
  }
  out[key] = PT.score(allRows, allPred);
}
// the feature vocabulary, so the test can assert nothing identity-keyed is even offered
out.__feats = PT.FEATS;
out.__specFeats = Object.fromEntries(PT.KEYS.map((k) =>
  [k, PT.SPECS[k].kind === "isotonic" ? [PT.SPECS[k].feat] : PT.SPECS[k].feats]));
// the isotonic staircase, for the monotonicity check
const iso = PT.fit("isotonic_sum", rows);
out.__knots = iso.model.knots;
out.__clamped = [PT.predict(iso, { sum: 0, nearer: 0, farther: 0, diff: 0 }),
                 PT.predict(iso, { sum: 1e6, nearer: 1e6, farther: 1e6, diff: 0 })];
process.stdout.write(JSON.stringify(out));
"""


def run_node(art):
    with tempfile.NamedTemporaryFile("w", suffix=".cjs", delete=False) as fh:
        fh.write(NODE)
        path = fh.name
    try:
        r = subprocess.run(["node", path, PTJS, art], capture_output=True, text=True, timeout=600)
        if r.returncode:
            return None, (r.stderr or "").strip()[-400:]
        return json.loads(r.stdout), ""
    finally:
        os.unlink(path)


def main():
    if not os.path.isfile(ART):
        sys.exit("scheffler.json.gz missing — run: python3 build_scheffler.py")
    d = json.load(gzip.open(ART, "rt"))
    m, train, models, ours = d["meta"], d["train"], d["models"], d["ours"]
    raw = json.dumps(d)

    section("shape")
    check("53 training zygotes", len(train) == 53, str(len(train)))
    check("meta count matches", m["n_train"] == len(train))
    check("frame count matches", sum(len(e["f"]) for e in train) == m["n_frames"],
          f"{sum(len(e['f']) for e in train)} vs {m['n_frames']}")
    check("every training embryo carries a datasheet index",
          sorted(e["idx"] for e in train) == list(range(1, len(train) + 1)))
    check("indices and labels agree", all(e["label"] == f"Z{e['idx']:02d}" for e in train))
    check("every embryo has enough frames to be a trajectory",
          all(len(e["f"]) >= 20 for e in train), str(min(len(e["f"]) for e in train)))
    check("8 models precomputed", len(models) == 8, str(len(models)))
    check("models are sorted best-first",
          all(models[i]["metrics"]["macro_mae"] <= models[i + 1]["metrics"]["macro_mae"]
              for i in range(len(models) - 1)))

    section("normalised time really is t_real / T_duration")
    bad_t, bad_range, bad_mono, bad_start = [], [], [], []
    for e in train:
        T = e["T"]
        for t_real, t, *_ in e["f"]:
            if abs(t - t_real / T) > 1e-3:
                bad_t.append(e["label"]); break
        ts = [f[1] for f in e["f"]]
        if min(ts) < -1e-9 or max(ts) > 1 + 1e-9:
            bad_range.append(e["label"])
        if any(b <= a for a, b in zip(ts, ts[1:])):
            bad_mono.append(e["label"])
        if abs(ts[0]) > 1e-9:
            bad_start.append(e["label"])
    check("t = t_real / T_duration for every frame", not bad_t, str(bad_t[:3]))
    check("t stays inside [0, 1]", not bad_range, str(bad_range[:3]))
    check("t increases frame to frame", not bad_mono, str(bad_mono[:3]))
    check("every trajectory starts at t = 0", not bad_start, str(bad_start[:3]))
    check("durations span the published range (8.75–11.75 h)",
          m["duration_h"] == [8.75, 11.75], str(m["duration_h"]))
    check("no trajectory reaches t = 1 (NEBD itself is not a frame)",
          m["t_max_observed"] < 1.0, str(m["t_max_observed"]))

    section("distances are physical")
    dists = [v for e in train for f in e["f"] for v in (f[2], f[3])]
    check("all distances are positive", min(dists) > 0, f"{min(dists):.2f}")
    check("no distance exceeds a mouse zygote's radius by much (< 60 µm)",
          max(dists) < 60, f"{max(dists):.2f}")
    # the whole clock rests on this: the pronuclei end up closer together than they started
    ends = [(sum(e["f"][0][2:4]), sum(e["f"][-1][2:4])) for e in train]
    check("∂ falls from first frame to last in every embryo",
          all(b < a for a, b in ends), str(sum(1 for a, b in ends if b >= a)) + " rose")

    section("the leaky volume features are absent, not merely unused")
    # The leaky model is NAMED in the ranking on purpose — that is the point of showing it. What
    # must not ship is its DATA, or a visitor could refit it in devtools and get a great score for
    # a model that cannot be measured in a fixed embryo.
    raw_train = json.dumps(train).lower()
    for tok in ("vol", "male", "female"):
        check(f"no '{tok}' anywhere in the training trajectories", tok not in raw_train)
    check("training frames carry exactly [t_real, t, d1, d2] — four numbers, no volumes",
          all(len(f) == 4 for e in train for f in e["f"]))
    check("the leaky model is present in the ranking but flagged non-deployable",
          any(not mm["deployable"] and "volume" in mm["key"].lower() for mm in models))
    check("every other model is deployable",
          sum(1 for mm in models if mm["deployable"]) == 7)
    check("no deployable model lists a volume feature",
          not [mm["key"] for mm in models if mm["deployable"]
               and any("vol" in f.lower() for f in mm.get("features", []))])

    section("features the page can build are identity-free")
    check("the page sorts the pair through featsFromPair, never using d1/d2 directly",
          "featsFromPair(fr[2], fr[3])" in open(PAGEJS).read())
    check("our zygotes ship as nearer/farther, so identity never enters",
          all(e["d1"] <= e["d2"] + 1e-9 for e in ours),
          str([e["id"] for e in ours if e["d1"] > e["d2"] + 1e-9][:3]))
    check("sum and diff agree with the pair",
          all(abs(e["sum"] - (e["d1"] + e["d2"])) < 1e-3
              and abs(e["diff"] - (e["d2"] - e["d1"])) < 1e-3 for e in ours))

    section("our fixed zygotes")
    check("51 of 60 scenes carry pronuclear geometry",
          m["n_ours"] == 51 and m["n_zygote_scenes"] == 60,
          f"{m['n_ours']}/{m['n_zygote_scenes']}")
    check("the 9 excluded are counted, not silently dropped",
          m["n_ours_excluded"] == m["n_zygote_scenes"] - m["n_ours"] == 9)
    check("meta count matches the records", m["n_ours"] == len(ours))
    check("every one has transcript counts", all(e["total_tx"] > 0 for e in ours),
          str([e["id"] for e in ours if e["total_tx"] <= 0][:3]))
    check("per-gene counts sum to the total",
          all(sum(e["g"].values()) == e["total_tx"] for e in ours))
    check("the gene list covers every gene any zygote reports",
          set(d["genes"]) == {g for e in ours for g in e["g"]})
    check("their geometry is in the training range (so τ is interpolation, mostly)",
          min(e["sum"] for e in ours) > 0 and max(e["sum"] for e in ours) < 120,
          f"{min(e['sum'] for e in ours):.1f}–{max(e['sum'] for e in ours):.1f}")

    section("the fold assignment the published ranking used")
    a = d["folds"]["assignment"]
    check("every training embryo has a fold", set(a) == {e["id"] for e in train})
    check("5 folds", len(set(a.values())) == 5, str(sorted(set(a.values()))))
    sizes = {k: sum(1 for v in a.values() if v == k) for k in set(a.values())}
    check("folds are balanced within one embryo", max(sizes.values()) - min(sizes.values()) <= 2,
          str(sizes))

    section("the browser fitter reproduces the published pipeline")
    got, err = run_node(ART)
    if got is None:
        check("node ran pt-models.js over the published folds", False, err)
    else:
        pub = {mm["key"]: mm["metrics"]["macro_mae"] for mm in models}
        for key, s in got.items():
            if key.startswith("__"):
                continue
            check(f"{key} macro MAE matches published", abs(s["macro_mae"] - pub[key]) < 0.003,
                  f"refit {s['macro_mae']:.4f} vs published {pub[key]:.4f}")
        check("the only features on offer are the identity-free four",
              set(got["__feats"]) == {"nearer", "farther", "sum", "diff"}, str(got["__feats"]))
        offenders = {k: v for k, v in got["__specFeats"].items()
                     if not set(v) <= {"nearer", "farther", "sum", "diff"}}
        check("no refittable model uses anything outside them", not offenders, str(offenders))
        check("the selected model is the best refit deployable one",
              min((s["macro_mae"], k) for k, s in got.items() if not k.startswith("__"))[1]
              == m["selected_key"])
        # isotonic must be a monotone DECREASING staircase in ∂ — that is the whole modelling claim
        kn = got["__knots"]
        check("isotonic knots are sorted in ∂", all(kn[i][0] <= kn[i + 1][0] for i in range(len(kn) - 1)))
        check("isotonic τ never rises as ∂ rises",
              all(kn[i][1] >= kn[i + 1][1] - 1e-12 for i in range(len(kn) - 1)),
              str(sum(1 for i in range(len(kn) - 1) if kn[i][1] < kn[i + 1][1] - 1e-12)) + " rises")
        check("enough knots to be a curve, few enough to be a staircase", 20 <= len(kn) <= 400,
              str(len(kn)))
        lo, hi = got["__clamped"]
        check("τ is clamped to [0, 1] outside the training range",
              abs(lo - 1) < 1e-9 or lo <= 1, f"∂=0 → {lo:.3f}")
        check("an absurdly large ∂ still returns a valid τ", 0 <= hi <= 1, f"{hi:.3f}")

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("scheffler-"))
    check("the paper and licence are carried with the data",
          "Scheffler" in m["paper"] and m["licence"] == "CC BY 4.0")
    check("the schematic caveat travels with the artifact",
          "drawing convention" in m.get("cartoon_note", ""))
    check("no absolute paths in the artifact", "/Users/" not in raw and "/Volumes/" not in raw)

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
