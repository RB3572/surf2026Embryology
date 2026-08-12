#!/usr/bin/env python3
"""Build data/scheffler.json.gz — everything the pseudotime-training page needs, in one file.

The page lets you refit the clock yourself: pick which embryos are held out, run the fit in the
browser, and see the knots, the held-out accuracy, and what the resulting model says about our own
fixed zygotes. So this artifact ships three things.

1. THE TRAINING COHORT — the live-imaged trajectories
   53 untreated zygotes, 2,057 frames, from the public source-data workbook of Scheffler et al.
   2021 (Nat Commun 12:841, CC BY 4.0), untreated condition only. Per frame: real elapsed time
   t_real, the two pronucleus-to-cell-centre distances, and the normalised true time

       t = t_real / T_duration        0 at pronuclear formation, 1 at NEBD

   ⚠️ The workbook gives DISTANCES ONLY — no angles, no cell outline, no cell radius. The cartoon
   the page draws is therefore a schematic: the radial distances are real data, where the two
   pronuclei sit around the circle is a drawing convention, and the cell outline is invented.

2. THE PRECOMPUTED MODEL COMPARISON — the eight configurations already evaluated offline, with
   their held-out metrics and the fold assignment they were scored under, so the page's static
   ranking is the published one rather than something recomputed and subtly different.

   Two feature families are marked non-deployable and stay marked here: anything keyed to
   male/female pronuclear identity (a fixed zygote gives no reliable way to tell them apart) and
   anything using the published relative volumes (they are normalised to each pronucleus's own
   FUTURE endpoint volume, so they encode elapsed time and cannot be measured in a snapshot).

3. OUR FIXED ZYGOTES — the cohort a refitted model gets applied to, with their geometry and their
   transcript counts so the page can plot abundance against inferred pseudotime.

   ⚠️ 51 of our 60 zygote scenes carry both pronuclei. The other nine come from label stacks that
   resolve only the cytoplasm, so no pronuclear distance exists and no pseudotime can be assigned
   to them by any model. They are excluded here and the count is recorded, rather than being
   quietly absorbed into a "60".

Output: data/scheffler.json.gz
"""
import csv
import glob
import gzip
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SEG = os.path.join(DATA, "segments")
TRAIN_CSV = os.path.join(HERE, "..", "embryo_viewer", "calibration_data", "scheffler2021",
                         "scheffler_2021_control_zygote_trajectories.csv")
CALIB = os.path.join(HERE, "..", "embryo_viewer", "public", "data", "pseudotime_calibration",
                     "calibration.json")
OURS = os.path.join(HERE, "..", "embryo_viewer", "public", "data", "pronuclei_pseudotime.json")
OUT = os.path.join(DATA, "scheffler.json.gz")

VERSION = "scheffler-1.0.0"
MIN_TX = 1          # a gene needs at least this many transcripts in an embryo to be listed


def load_training():
    rows = list(csv.DictReader(open(TRAIN_CSV)))
    by = defaultdict(list)
    dur = {}
    for r in rows:
        eid = r["embryo_id"]
        dur[eid] = float(r["migration_duration_h"])
        by[eid].append((
            float(r["time_h"]),
            float(r["normalized_time_tau"]),
            float(r["male_to_center_um"]),
            float(r["female_to_center_um"]),
        ))
    out = []
    for i, eid in enumerate(sorted(by), start=1):
        fr = sorted(by[eid])
        out.append({
            "id": eid,
            "idx": i,                                   # the datasheet index, Z01..Z53
            "label": f"Z{i:02d}",
            "T": round(dur[eid], 4),
            # [t_real, t, d1 (male), d2 (female)] — d1+d2 is derived in the browser, and
            # nearer/farther are min/max of the pair, so nothing identity-keyed ships as a feature
            "f": [[round(a, 4), round(b, 6), round(c, 4), round(d, 4)] for a, b, c, d in fr],
        })
    return out, len(rows)


def load_models():
    d = json.load(open(CALIB))
    keep = ("key", "label", "features", "n_features", "complexity", "deployable", "note")
    models = []
    for m in d["models"]:
        rec = {k: m[k] for k in keep if k in m}
        rec["metrics"] = m["metrics"]
        models.append(rec)
    models.sort(key=lambda m: m["metrics"]["macro_mae"])
    return models, d.get("folds", {}), d.get("selected", {}).get("key")


def load_ours():
    ours = json.load(open(OURS))["embryos"]
    gene_of = {}
    for path in sorted(glob.glob(os.path.join(SEG, "Zygote__*.json.gz"))):
        s = json.load(gzip.open(path, "rt"))
        counts = {g: len(rec["x"]) for g, rec in (s.get("transcripts") or {}).items()
                  if len(rec.get("x", [])) >= MIN_TX}
        gene_of[s["id"]] = counts

    out, no_counts = [], 0
    for e in ours:
        f = e.get("features") or {}
        if "distance_sum_um" not in f:
            continue
        counts = gene_of.get(e["id"])
        if counts is None:
            no_counts += 1
            counts = {}
        out.append({
            "id": e["id"], "label": e.get("label") or e["id"],
            "d1": round(f["nearer_to_center_um"], 4),      # identity-free: nearer / farther
            "d2": round(f["farther_to_center_um"], 4),
            "sum": round(f["distance_sum_um"], 4),
            "diff": round(f["distance_difference_um"], 4),
            "tau_ref": e.get("tau"),                       # the shipped model's answer, for reference
            "qc": e.get("qc"),
            "total_tx": int(sum(counts.values())),
            "g": counts,
        })
    return out, no_counts


def main():
    for p, what in ((TRAIN_CSV, "training CSV"), (CALIB, "calibration.json"), (OURS, "our zygotes")):
        if not os.path.isfile(p):
            sys.exit(f"missing {what}: {p}")

    train, n_frames = load_training()
    models, folds, selected = load_models()
    ours, no_counts = load_ours()

    n_zyg_scenes = len(glob.glob(os.path.join(SEG, "Zygote__*.json.gz")))
    genes = sorted({g for e in ours for g in e["g"]})
    durs = [e["T"] for e in train]
    taus = [f[1] for e in train for f in e["f"]]

    meta = {
        "version": VERSION,
        "paper": "Scheffler et al. 2021, Two mechanisms drive pronuclear migration in mouse "
                 "zygotes, Nat Commun 12:841",
        "doi": "10.1038/s41467-021-21020-x",
        "licence": "CC BY 4.0",
        "sheets": ["Figure 1b (pronucleus-to-cell-centre distance)",
                   "Figure S1m (per-embryo migration duration)"],
        "n_train": len(train), "n_frames": n_frames,
        "frames_per_embryo": [min(len(e["f"]) for e in train), max(len(e["f"]) for e in train)],
        "duration_h": [min(durs), max(durs)],
        "t_max_observed": round(max(taus), 4),
        "n_ours": len(ours),
        "n_zygote_scenes": n_zyg_scenes,
        "n_ours_excluded": n_zyg_scenes - len(ours),
        "n_ours_without_counts": no_counts,
        "n_genes": len(genes),
        "selected_key": selected,
        "notation": {
            "t_real": "real hours since pronuclear formation",
            "T_duration": "that embryo's pronuclear formation -> NEBD interval",
            "t": "normalised TRUE time = t_real / T_duration",
            "tau": "PREDICTED normalised time",
            "d": "summed pronucleus-to-cell-centre distance",
        },
        # the cartoon is a schematic; only the radii are data
        "cartoon_note": "The workbook gives distances only. Radial distances are real; the angular "
                        "placement of the two pronuclei and the cell outline are drawing "
                        "conventions, not measurements.",
    }

    doc = {"meta": meta, "train": train, "models": models, "folds": folds,
           "ours": ours, "genes": genes}
    with gzip.open(OUT, "wt") as fh:
        json.dump(doc, fh, separators=(",", ":"))

    print(f"  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1e3:.0f} KB)")
    print(f"  training: {len(train)} zygotes, {n_frames} frames, "
          f"{meta['frames_per_embryo'][0]}-{meta['frames_per_embryo'][1]} per embryo, "
          f"durations {durs and min(durs)}-{max(durs)} h")
    print(f"  models: {len(models)} precomputed ({sum(1 for m in models if m['deployable'])} deployable), "
          f"selected = {selected}")
    print(f"  ours: {len(ours)} zygotes with pronuclear geometry of {n_zyg_scenes} scenes "
          f"({meta['n_ours_excluded']} have no resolvable pronuclei), {len(genes)} genes")
    if no_counts:
        print(f"  !! {no_counts} of ours have geometry but no transcript counts")


if __name__ == "__main__":
    main()
