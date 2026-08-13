#!/usr/bin/env python3
"""Checks on data/sperm_pairing.json — maternal vs paternal, paired per zygote (figure 4.16).

The failure mode this artifact has to be defended against is not a wrong number, it is a
CIRCULAR one. Every comparison it ships is one of the four tests that decided which pronucleus is
maternal, re-asked as a measurement, so the leave-one-out variant is not a nicety — it is the
result. These checks make sure it is really leave-one-out:

  · the dropped test must actually be absent from the vote that labelled each zygote
  · a hand call must survive the leave-one-out unchanged, because it is not a vote
  · a zygote whose remaining tests tie must be dropped, not resolved by a coin flip
  · the naive and leave-one-out variants must actually differ somewhere, or the whole exercise
    is decoration

and, separately, that the numbers are recomputable from the rows the artifact itself ships.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

ART = os.path.join(ROOT, "data", "sperm_pairing.json")
ASSIGN = os.path.join(ROOT, "data", "pronuclei_assignments.json")

passed = failed = 0


def check(name, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1; print(f"  PASS  {name}")
    else:
        failed += 1; print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))


def section(t):
    print(f"\n[{t}]")


def d3(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def main():
    if not os.path.isfile(ART):
        sys.exit("sperm_pairing.json missing — run: python3 build_sperm_pairing.py")
    d = json.load(open(ART))
    m, C = d["meta"], d["comparisons"]
    recs = {r["id"]: r for r in json.load(open(ASSIGN))["embryos"]}
    manual = {i for i, r in recs.items() if (r.get("consensus") or {}).get("manual")}
    UM = m["unit_um_per_plot"]

    section("shape")
    check("the three comparisons are present",
          sorted(C) == ["com_polar", "com_sperm", "volume"], str(sorted(C)))
    for k, c in C.items():
        check(f"{k}: both variants", sorted(c["variants"]) == ["all", "loo"])
        check(f"{k}: it names the test(s) it would restate", bool(c["drops"]))
        check(f"{k}: n matches the row count",
              all(v["n"] == len(v["rows"]) for v in c["variants"].values()))
        check(f"{k}: no duplicate zygotes",
              all(len({r['id'] for r in v['rows']}) == v["n"] for v in c["variants"].values()))

    section("the leave-one-out really leaves the test out")
    DROPS = {k: set(c["drops"]) for k, c in C.items()}
    for k, c in C.items():
        bad = []
        for r in c["variants"]["loo"]["rows"]:
            if r["manual"]:
                continue                     # a hand call is not a vote; it is exempt by design
            rec = recs[r["id"]]
            votes = [t["female"] for kk, t in (rec.get("tests") or {}).items()
                     if kk not in DROPS[k] and isinstance(t, dict) and t.get("female") is not None]
            n0 = sum(1 for v in votes if v == 0)
            want = 0 if n0 > len(votes) - n0 else 1
            # the maternal pronucleus is pron[want]; check the shipped measurement is that one
            mat = rec["pron"][want]
            if k == "volume":
                got = abs(float(mat["volume"]) - r["m"]) < 0.5
            elif k == "com_sperm":
                got = abs(d3(mat["com_plot"], rec["sperm_plot"]) * UM - r["m"]) < 1e-3
            else:
                got = abs(d3(mat["com_plot"], rec["polar"]["com_plot"]) * UM - r["m"]) < 1e-3
            if not got or n0 == len(votes) - n0:
                bad.append(r["id"])
        check(f"{k}: every leave-one-out row uses the vote WITHOUT {'/'.join(sorted(DROPS[k]))}",
              not bad, str(bad[:4]))

    section("ties are dropped, not resolved")
    for k, c in C.items():
        loo_ids = {r["id"] for r in c["variants"]["loo"]["rows"]}
        tied = []
        for eid, rec in recs.items():
            if eid in manual:
                continue
            votes = [t["female"] for kk, t in (rec.get("tests") or {}).items()
                     if kk not in DROPS[k] and isinstance(t, dict) and t.get("female") is not None]
            n0 = sum(1 for v in votes if v == 0)
            if votes and n0 == len(votes) - n0 and eid in loo_ids:
                tied.append(eid)
        check(f"{k}: no zygote with a tied remaining vote is used", not tied, str(tied[:4]))
        check(f"{k}: every dropped zygote records why",
              all(x.get("reason") for x in c["variants"]["loo"]["dropped"]))

    section("hand calls survive the leave-one-out")
    for k, c in C.items():
        for vn in ("all", "loo"):
            rows = {r["id"]: r for r in c["variants"][vn]["rows"]}
            flagged = {i for i, r in rows.items() if r["manual"]}
            check(f"{k}/{vn}: every manual zygote present is flagged manual",
                  flagged == (manual & set(rows)), f"{sorted(flagged)} vs {sorted(manual & set(rows))}")
            # and its maternal side must be the HAND CALL, not a recomputed vote
            bad = []
            for i in flagged:
                rec = recs[i]
                mat = rec["pron"][rec["consensus"]["female"]]
                if k == "volume":
                    got = float(mat["volume"])
                elif k == "com_sperm":
                    got = d3(mat["com_plot"], rec["sperm_plot"]) * UM
                else:
                    got = d3(mat["com_plot"], rec["polar"]["com_plot"]) * UM
                if abs(got - rows[i]["m"]) > (0.5 if k == "volume" else 1e-3):
                    bad.append(i)
            check(f"{k}/{vn}: the hand call is what was measured", not bad, str(bad))
        check(f"{k}: the manual count is the real count",
              all(c["variants"][vn]["n_manual"] ==
                  sum(1 for r in c["variants"][vn]["rows"] if r["manual"]) for vn in ("all", "loo")))

    section("the leave-one-out is not decoration")
    # if dropping the test changed nothing anywhere, the circularity claim would be untestable
    moved = []
    for k, c in C.items():
        a, l = c["variants"]["all"], c["variants"]["loo"]
        if a["n"] != l["n"] or abs(a["median_diff"] - l["median_diff"]) > 1e-9:
            moved.append(k)
    check("dropping the test changes at least one comparison", moved, str(moved))
    for k, c in C.items():
        a, l = c["variants"]["all"], c["variants"]["loo"]
        print(f"        {k}: all n={a['n']} med {a['median_diff']:+.2f} P={a['p']:.4g}  |  "
              f"loo n={l['n']} med {l['median_diff']:+.2f} P={l['p']:.4g}")
    # the polar-body comparison is the fully circular one: on the naive consensus it should be
    # near-perfectly one-sided, and the leave-one-out should be visibly less so
    a = C["com_polar"]["variants"]["all"]; l = C["com_polar"]["variants"]["loo"]
    check("the polar-body comparison is near-total on the naive consensus, as circularity predicts",
          a["n_maternal_larger"] <= 1, f"{a['n_maternal_larger']}/{a['n']} maternal farther")
    check("...and the leave-one-out is visibly less one-sided",
          l["n_maternal_larger"] > a["n_maternal_larger"],
          f"{l['n_maternal_larger']}/{l['n']}")

    section("the summary statistics are the rows'")
    for k, c in C.items():
        for vn, v in c["variants"].items():
            diffs = sorted(r["m"] - r["p"] for r in v["rows"])
            n = len(diffs)
            med = diffs[n // 2] if n % 2 else (diffs[n // 2 - 1] + diffs[n // 2]) / 2
            check(f"{k}/{vn}: the median difference is the rows' median",
                  abs(med - v["median_diff"]) < 1e-3, f"{med:.4f} vs {v['median_diff']}")
            check(f"{k}/{vn}: the direction counts are the rows'",
                  v["n_maternal_larger"] == sum(1 for x in diffs if x > 0)
                  and v["n_paternal_larger"] == sum(1 for x in diffs if x < 0))
            check(f"{k}/{vn}: P is a probability", 0 <= v["p"] <= 1, str(v["p"]))

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("sperm-pairing-"))
    check("the specification is named", "4.16" in m["method"])
    check("the circularity is stated in the artifact itself",
          "LEAVE-ONE-OUT" in m["circularity"].upper() and "QUOTE" in m["circularity"].upper())
    check("the tie rule is stated", "tie" in m["ties"])
    check("the manual rule is stated", "not votes" in m["manual"])
    check("no absolute paths in the artifact", "/Users/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
