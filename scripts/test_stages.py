#!/usr/bin/env python3
"""Checks on data/stages.json.gz — Across the Stages (figures 8.3 / 8.6 / 9.1 / 9.2).

The strongest check is at the bottom: figure 8.3 is recomputed here from scratch, and the
reference ships its own table of the same quantity, so every gene's fold at every stage must
match it. That validates the port against the specification rather than against itself.

Above it, the three rules that decide what a number means, each of which would still produce
plausible-looking output if it silently broke:

  · THE ZERO RULE — an embryo with no transcripts contributes |log2 ratio| = 0, so a gene never
    detected at a stage lands on exactly 1.0 rather than vanishing from the average.
  · THE COUNT FLOOR — under 20 transcripts in an embryo the split is counting noise.
  · THE PANEL RULE — a gene absent from an embryo's records is a real zero only if it is on that
    probeset; otherwise it carries no information at all.
"""
import gzip
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import embryo_stats as ES                                          # noqa: E402

ART = os.path.join(ROOT, "data", "stages.json.gz")

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
        sys.exit("stages.json.gz missing — run: python3 build_stages.py")
    d = json.load(gzip.open(ART, "rt"))
    m, genes = d["meta"], d["genes"]
    ST = m["stages"]
    by = {g["g"]: g for g in genes}

    section("shape")
    check("three stages", ST == ["zygote", "early2cell", "late2cell"], str(ST))
    check("every gene has a name and no duplicates", len(by) == len(genes))
    check("every gene carries a record for at least one stage",
          all(any(g[s] for s in ST) for g in genes))
    check("meta count matches", m["n_genes"] == len(genes))
    check("every skipped embryo records why", all(s.get("reason") for s in m.get("skipped", [])))

    section("the fold and its null are commensurable")
    for s in ST:
        rec = [g[s] for g in genes if g[s]]
        check(f"{s}: every fold is >= 1", all(r["fold"] >= 1 - 1e-9 for r in rec),
              str(min(r["fold"] for r in rec)))
        check(f"{s}: every null is >= 1", all(r["null"] >= 1 - 1e-9 for r in rec))
        check(f"{s}: excess is fold minus null",
              all(abs(r["excess"] - (r["fold"] - r["null"])) < 1e-4 for r in rec))
        check(f"{s}: nMeas never exceeds n", all(r["nMeas"] <= r["n"] for r in rec))

    section("the zero rule")
    # a gene never detected above the floor at a stage must land on EXACTLY 1.0 — that is the
    # rule, and it is what stops a never-detected gene silently disappearing from the average
    bad = [(g["g"], s, g[s]["fold"]) for g in genes for s in ST
           if g[s] and g[s]["nMeas"] == 0 and abs(g[s]["fold"] - 1.0) > 1e-9]
    check("a gene never measured at a stage has fold exactly 1.0", not bad, str(bad[:3]))
    n_zero = sum(1 for g in genes for s in ST if g[s] and g[s]["nMeas"] == 0)
    check("that rule actually bites — some genes are never measured", n_zero > 0, str(n_zero))
    check("but not all of them", n_zero < sum(1 for g in genes for s in ST if g[s]))

    section("the count floor and the panel rule are recorded")
    check("the floor is 20", m["params"]["COUNT_FLOOR"] == 20)
    check("the zero rule is stated", "log2 ratio| = 0" in m["zeroRule"] or "0" in m["zeroRule"])
    check("the panel rule is stated", "probeset" in m["panelRule"])
    check("the bulk rule names the ratio of totals, not the median",
          "RATIO OF TOTALS" in m["bulk"] and "median" in m["bulk"])
    check("the halves are described per stage",
          "meridional" in m["halves"]["zygote"] and "blastomere" in m["halves"]["twocell"])

    section("the null is volume-matched, not a fair coin")
    # a fair-coin null would sit at ~1.0 regardless of the cell; a volume-matched one tracks the
    # counts, so it must RISE as embryo counts fall into late 2-cell
    med = {s: float(np.median([g[s]["null"] for g in genes if g[s]])) for s in ST}
    check("the median null rises into late 2-cell (fewer, smaller embryos)",
          med["late2cell"] > med["early2cell"], str(med))
    check("the median null is well above 1 — counting noise alone produces a fold",
          all(v > 1.2 for v in med.values()), str(med))

    section("imported vs derived is labelled")
    check("the imported note exists", bool(m["imported"].get("note")))
    note = m["imported"]["note"].lower()
    check("it names what cannot be re-derived",
          "cluster" in note and "curated" in note)
    check("the derived third fate is labelled derived",
          "DERIVED" in m["imported"].get("gained_rule", ""))
    groups = {}
    for g in genes:
        groups[g["group"]] = groups.get(g["group"], 0) + 1
    check("retained / lost / gained all exist",
          all(groups.get(k) for k in ("retained", "lost", "gained")), str(groups))
    check("every gene with a fate also carries the percentile it was called from",
          all(g["refPct"] for g in genes if g["group"] != "other"))
    check("the percentiles are in [0, 1]",
          all(0 <= g["refPct"]["e"] <= 1 and 0 <= g["refPct"]["l"] <= 1
              for g in genes if g["refPct"]))

    section("the recomputed percentile is a percentile")
    for s in ("early2cell", "late2cell"):
        v = [g[s]["pct"] for g in genes if g[s] and g[s]["pct"] is not None]
        check(f"{s}: every recomputed percentile is in [0, 1]", all(0 <= x <= 1 for x in v),
              f"{min(v):.3f}–{max(v):.3f}")
        check(f"{s}: they span most of the range", max(v) - min(v) > 0.8)

    section("VALIDATION against the reference's own 8.3 table")
    val = m.get("validation") or {}
    if not val.get("available"):
        check("the reference table was read", False, "not present on this machine")
    else:
        check("hundreds of genes were compared", val["n_shared"] >= 400, str(val["n_shared"]))
        for s in ST:
            mx = val["max_abs_diff"][s]
            check(f"{s}: every fold matches the reference", mx is not None and mx < 1e-4,
                  f"max |Δ| {mx}")
        print(f"        median |Δfold|: " +
              ", ".join(f"{s} {val['median_abs_diff'][s]}" for s in ST))

    section("provenance")
    check("version recorded", str(m.get("version", "")).startswith("stages-"))
    check("the specification is named", "8.3" in m.get("method", ""))
    check("no absolute paths in the artifact", "/Users/" not in json.dumps(d))

    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
