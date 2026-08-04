#!/usr/bin/env python3
"""Checks on data/alignment.json.gz — sperm alignment in the 2-cell stage.

The page's whole claim is that an ANCHOR supplies the two things the blastomere axis cannot, so
these guard the frame (orthonormal, axis really joins the two centres), the radius maps that the
browser slices to draw outlines, and — most importantly — the coverage arithmetic, because with
disjoint probesets an anchor is only ever a statement about the handful of embryos it can orient.
"""
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ART = os.path.join(ROOT, "data", "alignment.json.gz")

passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1; print(f"  PASS  {name}")
    else:
        failed += 1; print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))


def section(t):
    print(f"\n[{t}]")


def main():
    if not os.path.isfile(ART):
        sys.exit("alignment.json.gz missing — run: python3 build_alignment.py")
    d = json.load(gzip.open(ART, "rt"))
    m, emb, anchors = d["meta"], d["embryos"], d["anchors"]

    section("shape")
    check("has embryos", len(emb) > 60, str(len(emb)))
    check("meta count matches", m["n_embryos"] == len(emb))
    check("both 2-cell stages present", m["n_e2c"] > 0 and m["n_l2c"] > 0,
          f"e2c={m['n_e2c']} l2c={m['n_l2c']}")
    check("stage counts add up", m["n_e2c"] + m["n_l2c"] == len(emb))
    check("every embryo is a 2-cell stage", all(e["stage"] in ("e2c", "l2c") for e in emb))
    check("some embryos carry a sperm", m["n_sperm"] > 10, str(m["n_sperm"]))
    check("sperm count matches the embryos", m["n_sperm"] == sum(1 for e in emb if "sperm" in e))
    check("anchor count matches", m["n_anchors"] == len(anchors))

    section("the frame")
    bad_u = [e["id"] for e in emb if abs(np.linalg.norm(e["u"]) - 1) > 1e-3]
    check("axis is a unit vector", not bad_u, str(bad_u[:3]))
    bad_b = []
    for e in emb:
        e1, e2, u = np.array(e["e1"]), np.array(e["e2"]), np.array(e["u"])
        if (abs(np.linalg.norm(e1) - 1) > 1e-3 or abs(np.linalg.norm(e2) - 1) > 1e-3
                or abs(e1 @ e2) > 1e-3 or abs(e1 @ u) > 1e-3 or abs(e2 @ u) > 1e-3):
            bad_b.append(e["id"])
    check("(u, e1, e2) is an orthonormal basis", not bad_b, str(bad_b[:3]))
    bad_axis = []
    for e in emb:
        a, b, u = np.array(e["com_a"]), np.array(e["com_b"]), np.array(e["u"])
        if np.linalg.norm(u - (b - a) / np.linalg.norm(b - a)) > 1e-3:
            bad_axis.append(e["id"])
    check("u really is unit(COM_B − COM_A)", not bad_axis, str(bad_axis[:3]))
    mids = [np.allclose(e["mid"], (np.array(e["com_a"]) + np.array(e["com_b"])) / 2, atol=1e-2)
            for e in emb]
    check("mid is the midpoint of the two centres", all(mids))
    check("blastomere separation is physical (10-120 µm)",
          all(10 <= e["sep"] <= 120 for e in emb),
          str(sorted(e["sep"] for e in emb)[:2]))
    check("the two blastomeres are the two largest segments",
          all(e["vol_a"] >= e["vol_b"] for e in emb))

    section("radius maps (what the browser slices to draw outlines)")
    nt, npsi = m["t_bins"], m["psi_bins"]
    bad_len = [e["id"] for e in emb
               if len(e["map"]["a"]) != nt * npsi or len(e["map"]["b"]) != nt * npsi]
    check("map is nt × npsi for both blastomeres", not bad_len, str(bad_len[:3]))
    bad_q = [e["id"] for e in emb
             if min(e["map"]["a"]) < 0 or max(e["map"]["a"]) > 255
             or min(e["map"]["b"]) < 0 or max(e["map"]["b"]) > 255]
    check("quantised into a byte range", not bad_q, str(bad_q[:3]))
    check("every map reaches its own maximum (the scale is not dead)",
          all(max(e["map"]["a"]) > 200 and max(e["map"]["b"]) > 200 for e in emb))
    radii = [e["map"]["a_max"] for e in emb] + [e["map"]["b_max"] for e in emb]
    check("blastomere radii are physical (10-70 µm)",
          all(10 <= r <= 70 for r in radii), f"{min(radii):.1f}–{max(radii):.1f}")
    check("map L matches the centre separation",
          all(abs(e["map"]["L"] - e["sep"]) < 0.05 for e in emb))

    section("anchors and their coverage")
    ids = {e["id"] for e in emb}
    with_sperm = {e["id"] for e in emb if "sperm" in e}
    bad_cov = []
    for a in anchors:
        if a["kind"] == "polar":
            n = sum(1 for e in emb if "polar" in e)
            ns = sum(1 for e in emb if "polar" in e and "sperm" in e)
        else:
            n = sum(1 for e in emb if a["key"] in e["genes"])
            ns = sum(1 for e in emb if a["key"] in e["genes"] and "sperm" in e)
        if (n, ns) != (a["n"], a["n_sperm"]):
            bad_cov.append((a["key"], (a["n"], a["n_sperm"]), (n, ns)))
    check("every anchor's advertised coverage is its real coverage",
          not bad_cov, str(bad_cov[:3]))
    check("n_sperm never exceeds n", all(a["n_sperm"] <= a["n"] for a in anchors))
    check("n_sperm never exceeds the sperm-carrying embryos",
          all(a["n_sperm"] <= len(with_sperm) for a in anchors))
    check("every anchor clears min_emb", all(a["n"] >= m["min_emb"] for a in anchors))
    check("there is a polar-body anchor", any(a["kind"] == "polar" for a in anchors))

    # This is the constraint the page exists to be honest about: disjoint probesets mean NO gene
    # is measured everywhere. If a gene ever covered every embryo, the panels would have merged
    # and the coverage warnings on the page would be wrong.
    genes = [a for a in anchors if a["kind"] == "gene"]
    check("no gene anchor covers every embryo (probesets really are disjoint)",
          all(a["n"] < len(emb) for a in genes),
          f"max {max((a['n'] for a in genes), default=0)} of {len(emb)}")
    check("the best gene anchor still misses most sperm embryos",
          max((a["n_sperm"] for a in genes), default=0) < len(with_sperm),
          f"best {max((a['n_sperm'] for a in genes), default=0)} of {len(with_sperm)}")

    section("sperm")
    sp = [e for e in emb if "sperm" in e]
    check("every sperm names the blastomere it is in",
          all(e.get("sperm_side") in ("a", "b") for e in sp))
    bad_near = []
    for e in sp:
        p = np.array(e["sperm"])
        da, db = np.linalg.norm(p - np.array(e["com_a"])), np.linalg.norm(p - np.array(e["com_b"]))
        # the labelled segment wins, but a sperm should not be wildly nearer the OTHER centre
        near = "a" if da <= db else "b"
        if near != e["sperm_side"] and abs(da - db) > 0.5 * e["sep"]:
            bad_near.append((e["id"], e["sperm_side"], round(da, 1), round(db, 1)))
    check("the labelled blastomere is not contradicted by the geometry",
          not bad_near, str(bad_near[:3]))
    check("sperm sits within a plausible distance of its centre",
          all(min(np.linalg.norm(np.array(e["sperm"]) - np.array(e["com_a"])),
                  np.linalg.norm(np.array(e["sperm"]) - np.array(e["com_b"]))) < 90 for e in sp))

    section("gene records")
    bad_g = []
    for e in emb:
        for g, r in e["genes"].items():
            if len(r) != 3 or r[0] < 0 or r[1] < 0 or not (-180.1 <= r[2] <= 180.1):
                bad_g.append((e["id"], g, r)); break
    check("every gene record is [nA, nB, azimuth°]", not bad_g, str(bad_g[:3]))
    check("every gene clears min_tx",
          all(r[0] + r[1] >= m["min_tx"] for e in emb for r in e["genes"].values()))

    section("nucleus / polar-body surfaces")
    blobs = [(e["id"], b) for e in emb for b in e.get("blobs", [])]
    check("every embryo ships at least one blob", len(blobs) >= len(emb), str(len(blobs)))
    check("blob maps are nt × np", all(len(b["r"]) == b["nt"] * b["np"] for _, b in blobs))
    check("blob maps are byte-quantised", all(0 <= min(b["r"]) and max(b["r"]) <= 255 for _, b in blobs))
    check("every blob map reaches its own maximum",
          all(max(b["r"]) > 200 for _, b in blobs))
    # These are drawn as real plane-surface intersections, so a radius that is zero anywhere would
    # put a spike in the outline, and an implausible one would put a nucleus outside its cell.
    check("no blob has a zero radius anywhere (slices would spike)",
          all(min(b["r"]) > 0 for _, b in blobs),
          str([i for i, b in blobs if min(b["r"]) == 0][:3]))
    rr = [b["rmax"] for _, b in blobs]
    check("blob radii are physical (2-40 µm)", all(2 <= r <= 40 for r in rr),
          f"{min(rr):.1f}–{max(rr):.1f}")
    check("every blob is named nucleus or polar", all(b["kind"] in ("nucleus", "polar") for _, b in blobs))
    check("at most one polar body per embryo",
          all(sum(1 for b in e.get("blobs", []) if b["kind"] == "polar") <= 1 for e in emb))
    bad_in = []
    for e in emb:
        for b in e.get("blobs", []):
            c = np.array(b["c"])
            if min(np.linalg.norm(c - np.array(e["com_a"])),
                   np.linalg.norm(c - np.array(e["com_b"]))) > 90:
                bad_in.append((e["id"], b["label"]))
    check("blobs sit near the embryo, not off in space", not bad_in, str(bad_in[:3]))

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("alignment-"))
    check("no absolute paths in the artifact",
          "/Users/" not in json.dumps(d) and "/Volumes/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
