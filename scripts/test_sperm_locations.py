#!/usr/bin/env python3
"""Checks on data/sperm_locations.json — the sperm-location browser.

Guards the things the page would show wrongly rather than loudly: distances that are not
physical, a junction frame that does not agree with its own signed distance, and GFP entries
pointing at files that were never written.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ART = os.path.join(ROOT, "data", "sperm_locations.json")
passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok: passed += 1; print(f"  PASS  {name}")
    else: failed += 1; print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))


def main():
    if not os.path.isfile(ART):
        sys.exit("sperm_locations.json missing — run: python3 build_sperm_locations.py")
    d = json.load(open(ART))
    m, emb, gfp = d["meta"], d["embryos"], d["gfp"]

    print("\n[shape]")
    check("has embryos", len(emb) > 20, f"{len(emb)}")
    check("meta count matches", m["n"] == len(emb))
    check("stage counts add up", m["n_zygote"] + m["n_e2c"] + m["n_l2c"] <= len(emb))
    check("every embryo has a sperm position", all(len(e.get("sperm_plot") or []) == 3 for e in emb))
    check("ids unique", len({e["id"] for e in emb}) == len(emb))

    print("\n[distances are physical]")
    for key in ("cortex", "polar", "maternal", "paternal", "junction", "nucleus"):
        vals = [e["metrics"][key] for e in emb if key in e["metrics"]]
        if not vals: continue
        check(f"{key}: {len(vals)} values, all finite and 0-200 um",
              all(np.isfinite(v) and 0 <= v <= 200 for v in vals),
              f"range {min(vals):.1f}-{max(vals):.1f}")

    print("\n[the junction frame is self-consistent]")
    bad = []
    for e in emb:
        j = e.get("junction")
        if not j: continue
        if abs(np.linalg.norm(j["axis_um"]) - 1) > 1e-3: bad.append((e["id"], "axis not unit"))
        if abs(abs(j["signed_um"]) - e["metrics"]["junction"]) > 0.02: bad.append((e["id"], "|signed| != junction"))
        if j["side"] != (1 if j["signed_um"] > 0 else 0): bad.append((e["id"], "side disagrees with sign"))
    check("axis unit, |signed| == junction, side matches the sign", not bad, str(bad[:3]))
    njunc = sum(1 for e in emb if "junction" in e["metrics"])
    check("every 2-cell embryo got a junction",
          njunc == sum(1 for e in emb if e["stage"] == "twocell"), f"{njunc}")
    check("only 2-cell embryos have a junction",
          all(e["stage"] == "twocell" for e in emb if "junction" in e["metrics"]))
    check("only zygotes have maternal/paternal",
          all(e["stage"] == "zygote" for e in emb if "maternal" in e["metrics"]))

    print("\n[GFP stills]")
    gdir = os.path.join(ROOT, m["gfp_dir"])
    missing = [f for g in gfp.values() for f in g["files"].values()
               if not os.path.isfile(os.path.join(gdir, f))]
    check("every referenced GFP file exists", not missing, f"{len(missing)} missing")
    check("meta n_gfp matches", m["n_gfp"] == len(gfp))
    check("GFP entries belong to real embryos", all(k in {e['id'] for e in emb} for k in gfp))
    check("labelled z is inside the stack",
          all(1 <= g["z"] <= g["nz"] for g in gfp.values()))
    check("sperm pixel is inside the frame",
          all(g["x"] is None or (0 <= g["x"] < g["src_w"] and 0 <= g["y"] < g["src_h"])
              for g in gfp.values()))
    check("both channels present for every entry",
          all("ch0" in g["files"] and "ch1" in g["files"] for g in gfp.values()))

    # The projection is the whole reason these look like anything: plane 0 of every stack is a
    # saturated calibration frame, so a naive whole-stack max sets the display window ~150x above
    # the sperm and the embryo renders black. Guard the properties that keep that from returning.
    print("\n[max-Z projection window]")
    pj = {k: g["proj"] for k, g in gfp.items() if "proj" in g}
    check("every entry records how it was projected", len(pj) == len(gfp), f"{len(pj)}/{len(gfp)}")
    check("window is ordered and non-empty", all(p["z0"] < p["z1"] for p in pj.values()))
    check("window lies inside the stack",
          all(0 <= p["z0"] and p["z1"] <= gfp[k]["nz"] for k, p in pj.items()))
    check("the labelled sperm plane is inside the window",
          all(p["z0"] <= gfp[k]["z"] - 1 <= p["z1"] for k, p in pj.items()),
          str([k for k, p in pj.items() if not p["z0"] <= gfp[k]["z"] - 1 <= p["z1"]][:3]))
    # Not every stack carries a bad plane, so this is "the detector fires", not "always fires".
    n_drop = sum(1 for p in pj.values() if p["dropped"] > 0)
    check("the artifact-plane detector is actually firing",
          n_drop >= 0.5 * len(pj), f"{n_drop}/{len(pj)} stacks had planes dropped")
    check("but it never eats most of a stack",
          all(p["dropped"] < 0.2 * gfp[k]["nz"] for k, p in pj.items()),
          str(sorted((p["dropped"] / gfp[k]["nz"] for k, p in pj.items()), reverse=True)[:2]))
    check("projection uses fewer planes than the window spans (drops are real)",
          all(p["n"] <= p["z1"] - p["z0"] for p in pj.values()))
    check("the window is a MINORITY of the stack — not a whole-stack max in disguise",
          all((p["z1"] - p["z0"]) < 0.6 * gfp[k]["nz"] for k, p in pj.items()),
          str(sorted(((p["z1"] - p["z0"]) / gfp[k]["nz"] for k, p in pj.items()), reverse=True)[:2]))
    check("enough planes to be a projection, not a slice",
          all(p["n"] >= 15 for p in pj.values()),
          str(sorted(p["n"] for p in pj.values())[:3]))
    check("detection method recorded and known",
          all(p["how"] in ("405", "sperm") for p in pj.values()))
    n405 = sum(1 for p in pj.values() if p["how"] == "405")
    check("most embryos got a real 405-derived extent, not the fallback",
          n405 >= 0.7 * len(pj), f"{n405}/{len(pj)} from 405")

    print("\n[provenance]")
    check("version recorded", str(m.get("version", "")).startswith("sperm-locations-"))
    check("no absolute paths in the artifact", "/Users/" not in json.dumps(d) and "/Volumes/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
