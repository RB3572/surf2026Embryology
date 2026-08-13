#!/usr/bin/env python3
"""Checks on data/renders.json — every deck render, recomputed.

This artifact is itself a checker, so the risk is not that it computes a wrong number: it is that
it QUIETLY STOPS CHECKING. A panel with an empty check list, a verdict that does not follow from
its checks, or a family silently dropped would all leave the page reporting "31 agree" over a
smaller and smaller denominator. So:

  · every verdict must be derivable from that panel's own checks
  · "not checked" must never appear — a panel with nothing to check must say WHY (display
    sub-sample), and only the two families that legitimately have that reason may use it
  · every family in the deck must be present, with the number of panels it actually has
  · the near-miss count must be what the shipped diffs say, so "4 disagree, all by a molecule"
    cannot drift into hiding a real mismatch
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

ART = os.path.join(ROOT, "data", "renders.json")
SEG = os.path.join(ROOT, "data", "segments")

passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1; print(f"  PASS  {name}")
    else:
        failed += 1; print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))


def section(t):
    print(f"\n[{t}]")


def near(c):
    return (c.get("diff") is not None and c["diff"] <= 2) or \
           (c.get("rel") is not None and c["rel"] < 0.01)


def main():
    if not os.path.isfile(ART):
        sys.exit("renders.json missing — run: python3 build_renders.py")
    d = json.load(open(ART))
    m, P = d["meta"], d["panels"]

    section("shape")
    check("panel count matches", m["n_panels"] == len(P))
    check("ids are unique", len({p["id"] for p in P}) == len(P))
    check("every panel names its embryo, gene and scene",
          all(p["embryo"] and p["genes"] and p["scene"] for p in P))
    check("every panel carries a readout", all(p["readout"] for p in P))

    section("every scene it points at is actually here")
    missing = [p["scene"] for p in P if not os.path.isfile(os.path.join(SEG, p["scene"]))]
    check("no panel points at a scene this repo does not have", not missing, str(missing[:3]))

    section("the checker is still checking")
    # a verdict must FOLLOW from the checks; a panel that quietly lost its checks would otherwise
    # keep reporting "agrees" and inflate the headline
    bad = [p["id"] for p in P
           if p["verdict"] == "agrees" and (not p["checks"] or not all(c["ok"] for c in p["checks"]))]
    check("every 'agrees' has checks and they all pass", not bad, str(bad[:4]))
    bad = [p["id"] for p in P
           if p["verdict"] == "disagrees" and all(c["ok"] for c in p["checks"])]
    check("every 'disagrees' has a failing check", not bad, str(bad[:4]))
    check("nothing is silently unchecked",
          not [p["id"] for p in P if p["verdict"] == "not checked"])
    # only a display sub-sample may have no checks, and it must say so
    nochk = [p for p in P if not p["checks"]]
    check("a panel with no checks is a declared display sub-sample",
          all(p["verdict"] == "display sub-sample" for p in nochk),
          str([p["id"] for p in nochk if p["verdict"] != "display sub-sample"]))
    check("...and it explains why in its own caveat",
          all("SUB-SAMPLED" in (p.get("caveat") or "").upper() for p in nochk))
    check("the sub-sampled panels are only from 6.1",
          {p["fig"] for p in nochk} <= {"6.1"}, str({p["fig"] for p in nochk}))

    section("every family is present")
    got = {}
    for p in P:
        got[p["fig"]] = got.get(p["fig"], 0) + 1
    want = {"1.7": 9, "5.1": 1, "6.1": 6, "6.2": 8, "7.4": 2, "8.2": 6}
    for fig, n in want.items():
        check(f"figure {fig}: {n} panels", got.get(fig) == n, f"got {got.get(fig)}")
    check("8.7 is present and deduplicated", got.get("8.7", 0) == 6, str(got.get("8.7")))
    check("the meta lists every figure it built", sorted(m["figures"]) == sorted(got))
    print("        " + ", ".join(f"{k} {v}" for k, v in sorted(got.items())))

    section("the counts are the headline, and the headline is the counts")
    check("n_agree is the real count",
          m["n_agree"] == sum(1 for p in P if p["verdict"] == "agrees"))
    check("n_disagree is the real count",
          m["n_disagree"] == sum(1 for p in P if p["verdict"] == "disagrees"))
    check("n_subsampled is the real count",
          m["n_subsampled"] == sum(1 for p in P if p["verdict"] == "display sub-sample"))
    check("the three verdicts account for every panel",
          m["n_agree"] + m["n_disagree"] + m["n_subsampled"] == len(P))
    tiny = sum(1 for p in P if p["verdict"] == "disagrees"
               and all(near(c) for c in p["checks"] if not c["ok"]))
    check("the near-miss count is what the shipped diffs say", m["n_near_miss"] == tiny,
          f"{m['n_near_miss']} vs {tiny}")
    for p in P:
        if p["verdict"] == "disagrees":
            for c in p["checks"]:
                if not c["ok"]:
                    print(f"        {p['id']:28s} {c['name']}: deck {c['deck']} vs ours "
                          f"{c['ours']} (Δ{c['diff']}, {100 * (c['rel'] or 0):.2f}%)")

    section("every check carries its own size")
    allc = [c for p in P for c in p["checks"]]
    check("every check has both a difference and a relative size",
          all(c["diff"] is not None and c["rel"] is not None for c in allc))
    check("a passing check has zero difference",
          all(c["diff"] == 0 for c in allc if c["ok"] and "fold" not in c["name"]))
    check("the fold checks are the only ones allowed a tolerance",
          all(c["diff"] == 0 or "fold" in c["name"] for c in allc if c["ok"]))

    section("the density rule is applied where it matters")
    blas = [p for p in P if (p.get("highlight") or {}).get("kind") == "blastomeres"]
    check("the blastomere panels exist", len(blas) >= 10, str(len(blas)))
    check("each names both blastomeres",
          all(p["highlight"].get("hi") and p["highlight"].get("lo") for p in blas))
    # only the families that actually PRINT a fold need the density note; 8.2's 2-cell panels are
    # split by blastomere too but print a signed LFC, not a fold
    fold_figs = {"6.1", "6.2", "8.7"}
    check("every panel that prints a fold says it is a ratio of densities",
          all("DENSITIES" in (p.get("caveat") or "") for p in blas if p["fig"] in fold_figs),
          str([p["id"] for p in blas if p["fig"] in fold_figs
               and "DENSITIES" not in (p.get("caveat") or "")]))
    pl = [p for p in P if (p.get("highlight") or {}).get("kind") == "plane"]
    check("the plane panels ship a unit normal and an origin",
          all(abs(sum(x * x for x in p["highlight"]["normal"]) - 1) < 1e-5   # shipped at 6 dp
              and len(p["highlight"]["origin"]) == 3 for p in pl), str(len(pl)))
    check("...and say the plane is the deck's, not re-derived",
          all("deck's own" in (p.get("caveat") or "") for p in pl))

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("renders-"))
    check("the point of the page is stated in the artifact",
          "trusts most" in m["why"] and "check least" in m["why"])
    check("the counting convention is stated", "segment-label" in m["counts"]
          and "containment test" in m["counts"])
    check("the amplification device is recorded", "amplify_to" in m["amplify"])
    check("no absolute paths in the artifact", "/Users/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
