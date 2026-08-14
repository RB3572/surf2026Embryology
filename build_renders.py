#!/usr/bin/env python3
"""Build data/renders.json — every 3-D render in the deck, recomputed from this repo's own data.

WHAT THIS IS FOR. The deck's render figures are the ones a reader trusts most and can check least:
a picture of one embryo with a number printed under it. This page names the exact embryo and gene
behind each one, redraws it live from the repo's own scenes, and — the point — RECOMPUTES THE
PRINTED NUMBER and says whether it agrees.

  A render that disagrees with its own caption is the failure mode worth catching. It cannot be
  seen by looking at the picture, which is why every panel here carries a verdict rather than an
  image.

THE FAMILIES, and how each number is re-derived:

  1.7   three stages, nine embryos, one gene each — the printed count is that gene's transcripts
        in that embryo. ⚠️ The deck also carries `amplify_to` and `n_synthetic`: sparse panels are
        drawn with DUPLICATED transcripts so the point cloud reads at print size. That is a
        display device, not data, and this build records it rather than reproducing it.
  5.1   MuERV-L, almost entirely inside the pronuclei — counts split by segment label, with the
        maternal/paternal identity taken from the site's own consensus.
  6.1   Ddx20 in one early 2-cell, split by blastomere — counts per body label, and the fold as a
        ratio of DENSITIES, since the two blastomeres are not the same size.
  6.2   four genes x two stages, same rule.
  8.7   three genes x two stages, same rule.
  7.4   Trib3 and Dusp5 split contact half vs cell-edge half — the equal-volume split of each
        blastomere, the same construction the Contact Halves project uses.
  8.2   three genes x two stages, split by a plane whose normal the deck ships — so the plane is
        not re-derived here, only the counts on either side of it and the signed log fold.

WHERE THE NUMBERS COME FROM. Counts are cytoplasm-by-segment-label, never by a containment test.
Volumes are each label's own voxel volume from the scene's `segments` block. Both are exactly what
`embryo_stats` uses everywhere else, so a disagreement here is a real disagreement and not two
conventions passing each other.

Output: data/renders.json
"""
import collections
import csv
import json
import math
import os
import sys

import numpy as np

import embryo_stats as ES
# The embryo label is LOOKED UP (data/embryo_ids.json via embryo_naming), never derived and
# never read off a manifest — rebuilding an artifact must not quietly reintroduce a legacy
# name. embryo_label() falls back conspicuously when an embryo is missing from the lookup.
from embryo_naming import embryo_label

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(DATA, "renders.json")
REF = "/Users/rishib/Desktop/EmbyroPlayground/HighResSlideshowExports"

VERSION = "renders-1.0.0"
TOL_COUNT = 0            # counts must match exactly, or the disagreement is real
TOL_FOLD = 0.02          # folds are printed rounded; 2% is print precision, not slack

STAGE_DIR = {"Zygote": "Zygote", "Early2Cell": "Early2Cell", "Late2Cell": "Late2Cell",
             "Early 2 cell": "Early2Cell", "Late 2 cell": "Late2Cell"}


def scene_for(stage, eid):
    st = STAGE_DIR.get(stage, stage)
    p = os.path.join(DATA, "segments", f"{st}__{eid}.json.gz")
    return (st, p) if os.path.isfile(p) else (st, None)


def counts_by_label(sc, gene):
    """{segment label(str): count} for one gene, from the per-molecule segment label."""
    t = sc["transcripts"].get(gene)
    if t is None:
        return None, None
    s = np.asarray(t["s"], int)
    return collections.Counter(str(int(x)) for x in s), len(s)


def positions(sc, gene):
    t = sc["transcripts"].get(gene)
    if t is None:
        return None, None
    zs = sc["z_scale"]
    P = np.stack([np.asarray(t["x"], float), np.asarray(t["y"], float),
                  np.asarray(t["gz"], float) * zs], axis=1) * ES.PX
    return P, np.asarray(t["s"], int)


def read(path):
    p = os.path.join(REF, path)
    if not os.path.isfile(p):
        return None
    with open(p) as fh:
        return list(csv.DictReader(fh))


def verdict(checks):
    """A panel agrees only if every one of its checks does. Silence is not agreement."""
    if not checks:
        return "not checked"
    return "agrees" if all(c["ok"] for c in checks) else "disagrees"


def chk(name, deck, ours, ok, note=""):
    """A check carries its own SIZE as well as its verdict.

    One molecule sitting exactly on a cutting plane and a genuinely different measurement both
    read as "disagrees"; without the relative size the page cannot tell the reader which it is
    looking at, and every near-miss would be as alarming as a real error."""
    rel = adiff = None
    try:
        d = float(deck)
        adiff = abs(float(ours) - d)
        rel = adiff / max(abs(d), 1e-9)
    except (TypeError, ValueError):
        pass
    return {"name": name, "deck": deck, "ours": ours, "ok": bool(ok), "note": note,
            "rel": None if rel is None else round(rel, 6),
            "diff": None if adiff is None else round(adiff, 6)}


# ───────────────────────── the families ─────────────────────────
def family_1_7(panels):
    out = []
    for r in panels:
        st, p = scene_for(r["stage"], r["embryo"])
        if not p:
            continue
        sc = ES.read_scene(p)
        _, n = counts_by_label(sc, r["gene"])
        deck = int(r["found_n"])
        amp = int(r.get("amplify_to") or 0)
        syn = int(r.get("n_synthetic") or 0)
        checks = [chk("transcripts of this gene", deck, n, n == deck)]
        out.append({
            "fig": "1.7", "panel": r["panel"], "title": f"{r['gene']} · {r['stage']}",
            "stage": st, "embryo": r["embryo"], "label": embryo_label(r["embryo"]), "genes": [r["gene"]],
            "scene": f"{st}__{r['embryo']}.json.gz",
            "readout": [{"k": "transcripts", "v": f"{n:,}" if n is not None else "–"}],
            "checks": checks, "verdict": verdict(checks),
            "caveat": (f"The deck DUPLICATES this panel's transcripts up to {amp:,} for display "
                       f"({syn:,} synthetic points) so the cloud reads at print size. That is a "
                       f"drawing device; the count above is the real one."
                       if amp and syn else ""),
        })
    return out


def family_5_1(panels):
    out = []
    assign = {r["id"]: r for r in json.load(open(os.path.join(DATA,
                                                              "pronuclei_assignments.json")))["embryos"]}
    for r in panels:
        st, p = scene_for("Zygote", r["embryo"])
        if not p:
            continue
        sc = ES.read_scene(p)
        by, n = counts_by_label(sc, r["gene"])
        if by is None:
            continue
        a = assign.get(r["embryo"]) or {}
        pron = [str(x["label"]) for x in (a.get("pron") or [])]
        fem = (a.get("consensus") or {}).get("female")
        n_pron = sum(by.get(k, 0) for k in pron)
        checks = [chk("total transcripts", int(r["n_total"]), n, n == int(r["n_total"]))]
        if r.get("n_pronuclei"):
            checks.append(chk("inside the pronuclei", int(r["n_pronuclei"]), n_pron,
                              n_pron == int(r["n_pronuclei"])))
        ro = [{"k": "total", "v": f"{n:,}"},
              {"k": "inside the pronuclei", "v": f"{n_pron:,} ({100 * n_pron / max(n, 1):.1f}%)"}]
        if fem is not None and len(pron) == 2:
            ro.append({"k": "maternal ♀", "v": f"{by.get(pron[fem], 0):,}"})
            ro.append({"k": "paternal ♂", "v": f"{by.get(pron[1 - fem], 0):,}"})
        out.append({"fig": "5.1", "panel": r["panel"], "title": f"{r['gene']} · zygote",
                    "stage": st, "embryo": r["embryo"], "label": embryo_label(r["embryo"]), "genes": [r["gene"]],
                    "scene": f"{st}__{r['embryo']}.json.gz",
                    "highlight": {"kind": "labels", "labels": pron},
                    "readout": ro, "checks": checks, "verdict": verdict(checks)})
    return out


def family_blastomeres(fig, panels, title_of):
    """6.1, 6.2 and 8.7 all split one gene between the two blastomeres.

    The fold is a ratio of DENSITIES, not of counts: sister blastomeres are routinely 20% apart in
    volume, and a count ratio would read that as expression."""
    out = []
    seen = set()
    # 6.1's top2/top3/top4 are SUB-SAMPLED to the main panel's counts so the four read as one
    # series — a drawing device, like 1.7's amplification. Their shipped n_hi/n_lo are therefore
    # display numbers, not measurements, and checking a real count against them would be wrong.
    main = next((x for x in panels if x.get("panel") == "main"), None)
    for r in panels:
        key = (r.get("stage") or "Early2Cell", r["embryo"], r["gene"])
        if key in seen:                       # 8.7 ships the same panel twice under two scalings
            continue
        seen.add(key)
        st, p = scene_for(r.get("stage") or "Early2Cell", r["embryo"])
        if not p:
            continue
        sc = ES.read_scene(p)
        by, n = counts_by_label(sc, r["gene"])
        if by is None:
            continue
        vol = ES.seg_volumes(sc)
        hi, lo = str(r["hi_region"]), str(r["lo_region"])
        n_hi, n_lo = by.get(hi, 0), by.get(lo, 0)
        v_hi, v_lo = vol.get(hi, 0.0), vol.get(lo, 0.0)
        fold = ((n_hi / v_hi) / (n_lo / v_lo)) if n_lo and v_hi and v_lo else None
        subsampled = bool(main and r is not main and r["gene"] != main["gene"]
                          and r["n_hi"] == main["n_hi"] and r["n_lo"] == main["n_lo"])
        if subsampled:
            checks = []
        else:
            checks = [chk("count in the higher half", int(r["n_hi"]), n_hi, n_hi == int(r["n_hi"])),
                      chk("count in the lower half", int(r["n_lo"]), n_lo, n_lo == int(r["n_lo"]))]
            if r.get("fold") and fold:
                d = float(r["fold"])
                checks.append(chk("density fold", round(d, 4), round(fold, 4),
                                  abs(fold - d) / max(d, 1e-9) < TOL_FOLD))
        out.append({
            "fig": fig, "panel": r["panel"], "title": title_of(r),
            "stage": st, "embryo": r["embryo"], "label": embryo_label(r["embryo"]), "genes": [r["gene"]],
            "scene": f"{st}__{r['embryo']}.json.gz",
            "highlight": {"kind": "blastomeres", "hi": hi, "lo": lo},
            "readout": [
                {"k": "higher blastomere", "v": f"{n_hi:,} in {v_hi:,.0f} µm³"},
                {"k": "lower blastomere", "v": f"{n_lo:,} in {v_lo:,.0f} µm³"},
                {"k": "density fold", "v": f"{fold:.3f}" if fold else "–"},
            ],
            "checks": checks,
            "verdict": ("display sub-sample" if subsampled else verdict(checks)),
            "caveat": (("This panel is SUB-SAMPLED down to the main panel's transcript counts so "
                        "the four read as one series, so its printed counts are a drawing device "
                        "and there is nothing to check them against. The counts above are this "
                        "gene's real ones. ") if subsampled else "") +
                      ("The fold is a ratio of DENSITIES. Sister blastomeres differ in volume by "
                       f"{abs(v_hi - v_lo) / max(v_hi, v_lo, 1) * 100:.0f}% here, so a ratio of "
                       "raw counts would read that difference as expression."),
        })
    return out


def family_7_4(panels):
    """Trib3 and Dusp5, contact half vs cell-edge half — the same equal-volume split of each
    blastomere that the Contact Halves project uses."""
    out = []
    for r in panels:
        st, p = scene_for(r["stage"], r["embryo"])
        if not p:
            continue
        sc = ES.read_scene(p)
        bodies = ES.classify_body(sc)
        if len(bodies) != 2:
            continue
        P, S = positions(sc, r["gene"])
        if P is None:
            continue
        n_cont = n_edge = 0
        for a, b in ((bodies[0], bodies[1]), (bodies[1], bodies[0])):
            Va, Fa = ES.mesh_of(sc, a)
            ca = ES.vol_centroid(Va, Fa)
            cb = ES.vol_centroid(*ES.mesh_of(sc, b))
            u = ca - cb
            u /= np.linalg.norm(u)
            nrm, o = ES.equal_volume_plane(Va, Fa, u, ca, exact_total=ES.seg_volumes(sc)[a])
            sel = S == int(a)
            if not sel.any():
                continue
            t = (P[sel] - o) @ nrm
            n_cont += int((t <= 0).sum())
            n_edge += int((t > 0).sum())
        n_excl = int(len(P) - n_cont - n_edge)
        checks = [chk("contact half", int(r["n_contact"]), n_cont, n_cont == int(r["n_contact"])),
                  chk("cell-edge half", int(r["n_edge"]), n_edge, n_edge == int(r["n_edge"])),
                  chk("excluded (nucleus or polar body)", int(r["n_excluded"]), n_excl,
                      n_excl == int(r["n_excluded"]))]
        out.append({
            "fig": "7.4", "panel": r["panel"], "title": f"{r['gene']} · early 2-cell",
            "stage": st, "embryo": r["embryo"], "label": embryo_label(r["embryo"]), "genes": [r["gene"]],
            "scene": f"{st}__{r['embryo']}.json.gz",
            "highlight": {"kind": "contact"},
            "readout": [{"k": "contact half", "v": f"{n_cont:,}"},
                        {"k": "cell-edge half", "v": f"{n_edge:,}"},
                        {"k": "excluded", "v": f"{n_excl:,}"},
                        {"k": "deck: mean LFC", "v": f"{float(r['mean_lfc']):+.3f} "
                                                     f"over {r['n_embryos']} embryos"}],
            "checks": checks, "verdict": verdict(checks),
            "caveat": ("The deck's own caption prints "
                       f"{r.get('deck_contact')}/{r.get('deck_edge')} while its shipped table says "
                       f"{r.get('ref_contact')}/{r.get('ref_edge')} — a disagreement inside the "
                       "reference itself, noted in its notebook."
                       if r.get("deck_contact") and r.get("deck_contact") != r.get("ref_contact")
                       else ""),
        })
    return out


def family_8_2(panels):
    """The plane's normal is shipped by the deck, so the plane is NOT re-derived here — only the
    counts either side of it and the signed log fold. That keeps the check on the arithmetic
    rather than on a plane definition this build would be re-guessing."""
    out = []
    for r in panels:
        stage = "Zygote" if r["side"] == "zygote" else "Early2Cell"
        st, p = scene_for(stage, r["embryo"])
        if not p:
            continue
        sc = ES.read_scene(p)
        bodies = ES.classify_body(sc)
        P, S = positions(sc, r["gene"])
        if P is None or not bodies:
            continue
        inside = np.isin(S, [int(b) for b in bodies])
        V, F = ES.mesh_of(sc, bodies[0])
        com = ES.vol_centroid(V, F)
        try:
            nrm = np.array([float(r["normal_x"]), float(r["normal_y"]), float(r["normal_z"])])
        except (TypeError, ValueError):
            # the 2-cell panels split by BLASTOMERE, not by a plane, so they ship no normal
            nrm = None
        if nrm is None:
            if len(bodies) != 2:
                continue
            n0 = int((S == int(bodies[0])).sum())
            n1 = int((S == int(bodies[1])).sum())
            # `hi_is_a` is the deck's own statement of which side it called a; guessing by size
            # would silently swap the pair on every panel where the caption says otherwise
            hi_is_a = str(r.get("hi_is_a", "")).lower() in ("true", "1")
            n_a, n_b = (max(n0, n1), min(n0, n1)) if hi_is_a else (min(n0, n1), max(n0, n1))
            checks = [chk("count on side a", int(r["n_a"]), n_a, n_a == int(r["n_a"])),
                      chk("count on side b", int(r["n_b"]), n_b, n_b == int(r["n_b"]))]
            out.append({
                "fig": "8.2", "panel": r["panel"], "title": f"{r['gene']} · {r['side']}",
                "stage": st, "embryo": r["embryo"], "label": embryo_label(r["embryo"]), "genes": [r["gene"]],
                "scene": f"{st}__{r['embryo']}.json.gz",
                "highlight": {"kind": "blastomeres", "hi": str(bodies[0]), "lo": str(bodies[1])},
                "readout": [{"k": "side a", "v": f"{n_a:,}"}, {"k": "side b", "v": f"{n_b:,}"},
                            {"k": "deck: signed LFC", "v": f"{float(r['lfc']):+.3f}"}],
                "checks": checks, "verdict": verdict(checks),
                "caveat": "This 2-cell panel is split by BLASTOMERE, not by a plane — the deck "
                          "ships no normal for it, so the two sides are the two cell bodies.",
            })
            continue
        nrm /= np.linalg.norm(nrm)
        t = (P[inside] - com) @ nrm
        n_pos, n_neg = int((t > 0).sum()), int((t <= 0).sum())
        hi_is_a = str(r.get("hi_is_a", "")).lower() in ("true", "1")
        n_a, n_b = (n_pos, n_neg) if hi_is_a == (n_pos >= n_neg) else (n_neg, n_pos)
        checks = [chk("count on side a", int(r["n_a"]), n_a, n_a == int(r["n_a"])),
                  chk("count on side b", int(r["n_b"]), n_b, n_b == int(r["n_b"]))]
        out.append({
            "fig": "8.2", "panel": r["panel"], "title": f"{r['gene']} · {r['side']}",
            "stage": st, "embryo": r["embryo"], "label": embryo_label(r["embryo"]), "genes": [r["gene"]],
            "scene": f"{st}__{r['embryo']}.json.gz",
            "highlight": {"kind": "plane", "normal": [round(float(x), 6) for x in nrm],
                          "origin": [round(float(x), 4) for x in com]},
            "readout": [{"k": "side a", "v": f"{n_a:,}"}, {"k": "side b", "v": f"{n_b:,}"},
                        {"k": "deck: signed LFC", "v": f"{float(r['lfc']):+.3f}"}],
            "checks": checks, "verdict": verdict(checks),
            "caveat": "The plane's normal is the deck's own; only the counts either side of it "
                      "are recomputed, so this checks the arithmetic and not the plane rule.",
        })
    return out


def main():
    if not os.path.isdir(REF):
        sys.exit(f"the reference export is not on this machine: {REF}")
    fams = []
    p = read("Index1/1.7_three_stage_render_grid/data_panels.csv")
    if p:
        fams += family_1_7(p)
    p = read("Index5/5.1_muervl_render/data_panels.csv")
    if p:
        fams += family_5_1(p)
    p = read("Index6/6.1_split_transcript_render/data_panels.csv")
    if p:
        fams += family_blastomeres("6.1", p, lambda r: f"{r['gene']} · early 2-cell")
    p = read("Index6/6.2_gene_panel_twocell_renders/data_panels.csv")
    if p:
        fams += family_blastomeres("6.2", p, lambda r: f"{r['gene']} · {r['stage']}")
    p = read("Index8/8.7_muervl_ddx20_ltbp1_renders/data_panels.csv")
    if p:
        fams += family_blastomeres("8.7", p, lambda r: f"{r['gene']} · {r['stage']}")
    p = read("Index7/7.4_trib3_dusp5_renders/data_panels.csv")
    if p:
        fams += family_7_4(p)
    p = read("Index8/8.2_pard3_yap1_padi6_renders/data_panels.csv")
    if p:
        fams += family_8_2(p)

    for i, r in enumerate(fams, start=1):
        r["id"] = f"{r['fig']}-{r['panel']}"
        r["idx"] = i
    n_ok = sum(1 for r in fams if r["verdict"] == "agrees")
    n_bad = sum(1 for r in fams if r["verdict"] == "disagrees")
    n_sub = sum(1 for r in fams if r["verdict"] == "display sub-sample")
    # a disagreement of one molecule is a molecule sitting on a cutting plane, not an error, and
    # the summary says so rather than lumping it in with a real mismatch
    # "one or two molecules" needs BOTH tests: 1 out of 66,548 is small in relative terms and 1
    # out of 27 is not, but both are the same single molecule landing on the wrong side of a tie.
    n_tiny = sum(1 for r in fams if r["verdict"] == "disagrees"
                 and all((c["diff"] or 0) <= 2 or (c["rel"] or 0) < 0.01
                         for c in r["checks"] if not c["ok"]))

    doc = {
        "meta": {
            "version": VERSION,
            "method": "every 3-D render in the deck, its embryo and gene named and its printed "
                      "number recomputed from this repo's own scenes",
            "why": "a render is the figure a reader trusts most and can check least: one embryo, "
                   "one picture, one number. A caption that disagrees with its own data cannot be "
                   "seen by looking at the picture, so every panel here carries a verdict.",
            "counts": "cytoplasm-by-segment-label, never by a containment test; volumes are each "
                      "label's own voxel volume from the scene. Exactly what embryo_stats uses, "
                      "so a disagreement here is a real one and not two conventions passing.",
            "folds": "ratios of DENSITIES. Sister blastomeres are routinely 20% apart in volume, "
                     "and a ratio of raw counts would read that difference as expression.",
            "amplify": "the deck duplicates transcripts on sparse panels (`amplify_to`) so the "
                       "cloud reads at print size. That is a drawing device, recorded here and "
                       "not reproduced.",
            "tolerance": f"counts must match EXACTLY; folds within {TOL_FOLD:.0%}, which is the "
                         "precision they are printed at.",
            "n_panels": len(fams), "n_agree": n_ok, "n_disagree": n_bad,
            "n_subsampled": n_sub, "n_near_miss": n_tiny,
            "near_miss": "a disagreement of a molecule or two — or under 1% — is a molecule "
                         "sitting exactly on a "
                         "cutting plane, which is a tie-break convention rather than an error. "
                         "The page separates those from real mismatches instead of reporting one "
                         "count of failures.",
            "figures": sorted({r["fig"] for r in fams}),
        },
        "panels": fams,
    }
    with open(OUT, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(f"  wrote {os.path.relpath(OUT, HERE)} ({os.path.getsize(OUT)/1024:.0f} KB)")
    print(f"  {len(fams)} panels · {n_ok} agree · {n_bad} disagree "
          f"({n_tiny} of them by a molecule or two) · {n_sub} are display sub-samples")
    for r in fams:
        if r["verdict"] == "disagrees":
            bad = [c for c in r["checks"] if not c["ok"]]
            print(f"    {r['id']:34s} {r['embryo']:28s} " +
                  "; ".join(f"{c['name']}: deck {c['deck']} vs ours {c['ours']}" for c in bad))


if __name__ == "__main__":
    main()
