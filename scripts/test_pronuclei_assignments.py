#!/usr/bin/env python3
"""Checks on data/pronuclei_assignments.json — maternal / paternal pronucleus calls.

The point of this file is the hand-made calls. They are the one thing here that no script can
regenerate from the raw data on this machine, so they are exactly what a rebuild would silently
drop. Every override in pronuclei_assignments_manual.json must be present, must be marked as a
manual call rather than passing as a consensus, and must still have a scene for the page to draw.
"""
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
ASSIGN = os.path.join(DATA, "pronuclei_assignments.json")
MANUAL = os.path.join(DATA, "pronuclei_assignments_manual.json")

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
    if not os.path.isfile(ASSIGN):
        sys.exit("pronuclei_assignments.json missing")
    doc = json.load(open(ASSIGN))
    emb = doc["embryos"]
    by = {e["id"]: e for e in emb}
    manual = json.load(open(MANUAL))["overrides"]

    section("shape")
    check("every zygote has exactly two pronuclei", all(len(e["pron"]) == 2 for e in emb))
    check("pron[0] is the larger of the pair",
          all(e["pron"][0]["volume"] >= e["pron"][1]["volume"] for e in emb),
          str([e["id"] for e in emb if e["pron"][0]["volume"] < e["pron"][1]["volume"]][:3]))
    check("the two pronuclei are different segments",
          all(e["pron"][0]["label"] != e["pron"][1]["label"] for e in emb))
    check("a pronucleus is never also the polar body",
          all(not e["polar"] or e["polar"]["label"] not in
              {e["pron"][0]["label"], e["pron"][1]["label"]} for e in emb))
    check("every consensus names a real pronucleus index",
          all(e["consensus"]["female"] in (0, 1) for e in emb))
    check("ids are unique", len(by) == len(emb))
    check("every test that ran names a female index",
          all(t["female"] in (0, 1) for e in emb for t in e["tests"].values() if t))

    section("the volume test is what it claims")
    # "the smaller pronucleus is female" — with pron[1] the smaller, that is always index 1
    check("the volume test always calls index 1 female",
          all(e["tests"]["volume"]["female"] == 1 for e in emb))

    section("the sperm test runs exactly where a sperm exists")
    sperm_of = {e["id"] for e in json.load(open(os.path.join(DATA, "zygote_sperm.json")))["embryos"]
                if e.get("sperm_plot")}
    mismatch = [e["id"] for e in emb
                if bool(e["tests"].get("sperm")) != (e["id"] in sperm_of)]
    check("no zygote has a sperm test without a labelled sperm, or vice versa",
          not mismatch, str(mismatch[:4]))
    check("a sperm test always names the male as the other index",
          all(t["male"] == 1 - t["female"] for e in emb
              if (t := e["tests"].get("sperm"))))

    section("hand-made calls survive")
    for ov in manual:
        eid = ov["id"]
        e = by.get(eid)
        check(f"{eid} is present", e is not None)
        if not e:
            continue
        man = e.get("manual")
        check(f"{eid} is marked as a manual call", bool(man) and e["consensus"].get("manual") is True)
        if not man:
            continue
        want_f = next(int(l) for l, s in ov["assign"].items() if s == "female")
        want_m = next(int(l) for l, s in ov["assign"].items() if s == "male")
        check(f"{eid}: M{want_f} is female, M{want_m} is male",
              man["female_label"] == want_f and man["male_label"] == want_m,
              f"got F=M{man['female_label']} M=M{man['male_label']}")
        # the index the page colours by must point at the segment the call names
        check(f"{eid}: the female INDEX points at M{want_f}",
              e["pron"][man["female_index"]]["label"] == want_f)
        check(f"{eid}: the consensus follows the manual call",
              e["consensus"]["female"] == man["female_index"])
        check(f"{eid}: the reason is recorded", bool(man.get("reason")) and bool(man.get("by")))
        check(f"{eid}: both named segments are the pronuclei",
              {p["label"] for p in e["pron"]} == set(ov["pronuclei"]))
        check(f"{eid}: the polar body is the segment the override names",
              e["polar"] and e["polar"]["label"] == ov["polar_body"])
        # the page loads data/pronuclei/<id>.json.gz — an override with no scene is an empty page
        sc = os.path.join(DATA, "pronuclei", eid + ".json.gz")
        check(f"{eid}: a scene exists for the page to draw", os.path.isfile(sc))
        if os.path.isfile(sc):
            s = json.load(gzip.open(sc, "rt"))
            check(f"{eid}: that scene carries both pronuclei meshes",
                  all(str(l) in s["region_meshes"] for l in ov["pronuclei"]))
            check(f"{eid}: and the polar-body mesh",
                  str(ov["polar_body"]) in s["region_meshes"])

    section("the oocyte is not filed as a zygote")
    # 20260425_zygote_p2_3 is an oocyte with one pronucleus; it must not appear here at all
    check("20260425_zygote_p2_3 is absent (it is an oocyte, one pronucleus)",
          "20260425_zygote_p2_3" not in by)
    seg = os.path.join(DATA, "segments")
    check("its segments scene is filed under Oocyte",
          os.path.isfile(os.path.join(seg, "Oocyte__20260425_zygote_p2_3.json.gz"))
          and not os.path.isfile(os.path.join(seg, "Zygote__20260425_zygote_p2_3.json.gz")))
    sman = json.load(open(os.path.join(DATA, "segments_manifest.json")))["embryos"]
    row = next((r for r in sman if r.get("eid") == "20260425_zygote_p2_3"), None)
    check("the segments manifest calls it an Oocyte with an O- label",
          row and row["stage"] == "Oocyte" and row["label"].startswith("O-"),
          json.dumps(row) if row else "missing")

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
