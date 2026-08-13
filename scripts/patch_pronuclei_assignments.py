#!/usr/bin/env python3
"""Apply data/pronuclei_assignments_manual.json to data/pronuclei_assignments.json.

WHY THIS EXISTS. build_pronuclei_assignments.py needs the full-resolution label TIFFs and the
MERFISH atlas scenes, neither of which lives in this repo, so it cannot be re-run here. It also
drops any zygote whose two pronuclei its detector does not accept — which is the whole reason an
override is needed in the first place. This script adds those zygotes back from data the repo DOES
have: the per-label meshes in data/zygote/<id>.json.gz.

The geometry is recomputed rather than typed in. Volumes and centroids come from the closed
segmentation mesh by the divergence theorem, and shell distances from nearest vertex pairs — which
is a finer surface than the downsampled voxel cloud the original build used, so the two agree in
substance but not to the last decimal. Every number the page draws is derived; only the SEX CALL
is hand-made, and it is carried in a `manual` block so the page can say so.

Idempotent: re-running replaces the patched entries rather than duplicating them.
"""
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
ASSIGN = os.path.join(DATA, "pronuclei_assignments.json")
MANUAL = os.path.join(DATA, "pronuclei_assignments_manual.json")
ZY = os.path.join(DATA, "zygote")
SPERM = os.path.join(DATA, "zygote_sperm.json")
XY_UM = 0.15


def mesh_um(scene, label):
    """A segment's mesh vertices and faces in µm."""
    m = (scene.get("region_meshes") or {}).get(str(label))
    if not m:
        return None, None
    zs = scene["z_scale"]
    V = np.asarray(m["verts"], float).reshape(-1, 3)
    V = np.stack([V[:, 0] * XY_UM, V[:, 1] * XY_UM, V[:, 2] / zs], axis=1)
    return V, np.asarray(m["faces"], int).reshape(-1, 3)


def solid_stats(V, F):
    """Enclosed volume (µm³) and solid centroid of a closed mesh, by the divergence theorem —
    not the mean vertex position, which would be biased by uneven triangulation."""
    t = V[F]
    w = np.einsum("ij,ij->i", np.cross(t[:, 0], t[:, 1]), t[:, 2]) / 6.0
    vol = float(w.sum())
    c = (t[:, 0] + t[:, 1] + t[:, 2]) / 4.0
    com = (c * w[:, None]).sum(0) / w.sum()
    return abs(vol), com


def nearest_pair(A, B):
    """Closest vertex pair between two surfaces: (distance, point in A, point in B)."""
    best = (np.inf, None, None)
    for s in range(0, len(A), 4000):
        chunk = A[s:s + 4000]
        d = np.linalg.norm(chunk[:, None, :] - B[None, :, :], axis=2)
        i, j = np.unravel_index(int(np.argmin(d)), d.shape)
        if d[i, j] < best[0]:
            best = (float(d[i, j]), chunk[i], B[j])
    return best


def main():
    for p, what in ((ASSIGN, "pronuclei_assignments.json"), (MANUAL, "the overrides file")):
        if not os.path.isfile(p):
            sys.exit(f"missing {what}: {p}")
    doc = json.load(open(ASSIGN))
    manual = json.load(open(MANUAL))
    sperm_of = {e["id"]: e["sperm_plot"] for e in json.load(open(SPERM))["embryos"]
                if e.get("sperm_plot")}
    by_id = {e["id"]: e for e in doc["embryos"]}

    for ov in manual["overrides"]:
        eid = ov["id"]

        # ---- the common case: the build already produced this embryo ----
        # Keep its geometry and its tests EXACTLY as the build computed them (from the full-res
        # voxel masks) and change only the sex call. Recomputing from meshes here would quietly
        # put one embryo on a different estimator from the other fifty.
        if eid in by_id and not ov.get("include"):
            rec = by_id[eid]
            have = {p["label"] for p in rec["pron"]}
            if have != set(ov["pronuclei"]):
                print(f"  !! {eid}: override names pronuclei {ov['pronuclei']} but the artifact "
                      f"has {sorted(have)} — skipped")
                continue
            female_label = next(int(l) for l, s in ov["assign"].items() if s == "female")
            male_label = next(int(l) for l, s in ov["assign"].items() if s == "male")
            fi = next(i for i, p in enumerate(rec["pron"]) if p["label"] == female_label)
            calls = [t["female"] for t in rec["tests"].values() if t]
            n0, n1 = calls.count(0), calls.count(1)
            agrees = all(c == fi for c in calls) if calls else False
            rec["consensus"] = {"female": fi, "split": False, "n0": n0, "n1": n1, "manual": True}
            rec["manual"] = {
                "female_label": female_label, "male_label": male_label,
                "female_index": fi, "male_index": 1 - fi,
                "by": ov.get("by"), "date": ov.get("date"),
                "reason": ov.get("reason"), "caveat": ov.get("caveat"),
                "auto_agrees": agrees,
            }
            print(f"  updated {eid} (geometry untouched)")
            print(f"    manual: M{female_label} female / M{male_label} male  →  female index {fi}")
            ran = ", ".join(f"{k}=F{v['female']}" for k, v in rec["tests"].items() if v) or "none"
            print(f"    tests that ran: {ran}  "
                  f"({'all agree' if agrees else 'THE HAND CALL OVERRIDES THEM'})")
            continue

        # ---- the build dropped this embryo entirely: rebuild its geometry from the meshes ----
        scene_p = os.path.join(ZY, eid + ".json.gz")
        if not os.path.isfile(scene_p):
            print(f"  !! {eid}: no scene at {scene_p}")
            continue
        sc = json.load(gzip.open(scene_p, "rt"))
        zs = sc["z_scale"]
        to_plot = lambda p: [round(float(p[0]) / XY_UM, 1), round(float(p[1]) / XY_UM, 1),
                             round(float(p[2]) * zs, 1)]

        # ---- geometry, recomputed from the meshes ----
        pn = {}
        for lbl in ov["pronuclei"]:
            V, F = mesh_um(sc, lbl)
            if V is None:
                sys.exit(f"{eid}: label {lbl} has no mesh")
            vol, com = solid_stats(V, F)
            pn[lbl] = {"V": V, "vol": vol, "com": com}
        pbV, pbF = mesh_um(sc, ov["polar_body"])
        pb_vol, pb_com = solid_stats(pbV, pbF) if pbV is not None else (None, None)

        # the schema orders pron[0] = LARGER, pron[1] = SMALLER; keep that so the page's index
        # conventions (and every test's 0/1 answer) stay meaningful
        order = sorted(ov["pronuclei"], key=lambda l: -pn[l]["vol"])
        idx_of = {lbl: i for i, lbl in enumerate(order)}
        female_label = next(int(l) for l, s in ov["assign"].items() if s == "female")
        male_label = next(int(l) for l, s in ov["assign"].items() if s == "male")
        fi = idx_of[female_label]

        tests = {}
        if pbV is not None:
            dA = np.linalg.norm(pn[order[0]]["com"] - pb_com)
            dB = np.linalg.norm(pn[order[1]]["com"] - pb_com)
            f = 0 if dA < dB else 1
            tests["pb_com"] = {"female": f,
                               "line": [to_plot(pn[order[f]]["com"]), to_plot(pb_com)]}
            dA2, pA, qA = nearest_pair(pn[order[0]]["V"], pbV)
            dB2, pB, qB = nearest_pair(pn[order[1]]["V"], pbV)
            f2 = 0 if dA2 < dB2 else 1
            tests["pb_shell"] = {"female": f2,
                                 "line": [to_plot(pA if f2 == 0 else pB),
                                          to_plot(qA if f2 == 0 else qB)]}
        tests["volume"] = {"female": 1}                      # the smaller pronucleus
        sp_plot = sperm_of.get(eid)
        if sp_plot:
            sp = np.array([sp_plot[0] * XY_UM, sp_plot[1] * XY_UM, sp_plot[2] / zs])
            dA3 = np.linalg.norm(pn[order[0]]["V"] - sp, axis=1)
            dB3 = np.linalg.norm(pn[order[1]]["V"] - sp, axis=1)
            male_i = 0 if dA3.min() < dB3.min() else 1
            pt = (pn[order[0]]["V"][int(np.argmin(dA3))] if male_i == 0
                  else pn[order[1]]["V"][int(np.argmin(dB3))])
            tests["sperm"] = {"female": 1 - male_i, "male": male_i,
                              "line": [to_plot(pt), to_plot(sp)]}
        else:
            tests["sperm"] = None

        calls = [t["female"] for t in tests.values() if t]
        n0, n1 = calls.count(0), calls.count(1)

        rec = {
            "id": eid,
            "label": by_id.get(eid, {}).get("label") or _label_of(eid),
            "pron": [{"label": int(l), "com_plot": to_plot(pn[l]["com"]),
                      "volume": round(pn[l]["vol"])} for l in order],
            "polar": None if pbV is None else {"label": int(ov["polar_body"]),
                                               "com_plot": to_plot(pb_com),
                                               "volume": round(pb_vol)},
            "sperm_plot": None if not sp_plot else [round(float(x), 1) for x in sp_plot],
            "tests": tests,
            # the consensus IS the manual call — it is not a vote that happened to land here
            "consensus": {"female": fi, "split": False, "n0": n0, "n1": n1, "manual": True},
            "manual": {
                "female_label": female_label, "male_label": male_label,
                "female_index": fi, "male_index": 1 - fi,
                "by": ov.get("by"), "date": ov.get("date"),
                "reason": ov.get("reason"), "caveat": ov.get("caveat"),
                "auto_agrees": (n1 > n0) == (fi == 1) and n0 != n1,
            },
        }
        # the page renders from data/pronuclei/<id>.json.gz, and an embryo the original build
        # dropped has no such scene — emit the slim one it expects, from the meshes we do have
        pn_scene = os.path.join(DATA, "pronuclei", eid + ".json.gz")
        if not os.path.isfile(pn_scene):
            slim = {"id": eid, "z_scale": zs, "extents": sc["extents"],
                    "mask_labels": sc["mask_labels"],
                    "region_defaults": sc["region_defaults"],
                    "region_meshes": sc["region_meshes"],
                    "pron_labels": [int(l) for l in order]}
            os.makedirs(os.path.dirname(pn_scene), exist_ok=True)
            with gzip.open(pn_scene, "wt") as fh:
                json.dump(slim, fh, separators=(",", ":"))
            print(f"    wrote scene data/pronuclei/{eid}.json.gz "
                  f"({os.path.getsize(pn_scene)/1024:.0f} KB)")

        if eid in by_id:
            doc["embryos"][doc["embryos"].index(by_id[eid])] = rec
            print(f"  updated {eid}")
        else:
            doc["embryos"].append(rec)
            print(f"  added   {eid}")
        agree = "every automatic test agrees" if rec["manual"]["auto_agrees"] else "AUTO TESTS DISAGREE"
        print(f"    pron[0]=M{order[0]} ({round(pn[order[0]]['vol'])} µm³), "
              f"pron[1]=M{order[1]} ({round(pn[order[1]]['vol'])} µm³)")
        print(f"    manual: M{female_label} female / M{male_label} male  →  female index {fi}")
        print(f"    tests: " + ", ".join(f"{k}=F{v['female']}" for k, v in tests.items() if v)
              + f"  ({agree})")

    doc["embryos"].sort(key=lambda e: e["id"])
    with open(ASSIGN, "w") as fh:
        json.dump(doc, fh, indent=1)
    n_manual = sum(1 for e in doc["embryos"] if e.get("manual"))
    print(f"\n  wrote {os.path.relpath(ASSIGN, ROOT)} — {len(doc['embryos'])} zygotes "
          f"({n_manual} with a manual call)")


def _label_of(eid):
    """Fall back to the canonical display label if the embryo was never in the artifact."""
    try:
        sys.path.insert(0, ROOT)
        from embryo_naming import embryo_label
        return embryo_label(eid, "zygote")
    except Exception:                                          # noqa: BLE001
        return eid


if __name__ == "__main__":
    main()
