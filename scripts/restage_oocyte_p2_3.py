#!/usr/bin/env python3
"""Re-file 20260425_zygote_p2_3 as the OOCYTE it is.

The id says "zygote" because that is what the field of view was called at acquisition, and the
naming layer and data/manifest.json were both corrected long ago (embryo_naming already tests it as
`Oocyte__20260425_zygote_p2_3` → O-P3-fov3). The segments artifact never was, so this embryo was
still sitting in the Zygote group with a Z- label — and every downstream count that globs
`Zygote__*` was counting an oocyte as a zygote, including "60 zygote scenes".

It has ONE pronucleus, not two, which is why it is already absent from every project that needs a
pronuclear pair. Nothing here changes those; this only corrects the stage.

Renames the scene file, fixes the manifest row and the scene's own stage field. Idempotent: if it
has already been re-filed the script says so and changes nothing.

Downstream, re-run afterwards:  build_size.py, build_scheffler.py
"""
import glob
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SEG = os.path.join(ROOT, "data", "segments")
MAN = os.path.join(ROOT, "data", "segments_manifest.json")

EID = "20260425_zygote_p2_3"
OLD_ID, NEW_ID = f"Zygote__{EID}", f"Oocyte__{EID}"
NEW_LABEL = "O-P3-fov3"          # matches data/manifest.json, which was corrected already


def main():
    old_f = os.path.join(SEG, OLD_ID + ".json.gz")
    new_f = os.path.join(SEG, NEW_ID + ".json.gz")
    man = json.load(open(MAN))
    rows = man["embryos"]
    row = next((r for r in rows if r.get("eid") == EID), None)
    if row is None:
        sys.exit(f"{EID} is not in the segments manifest")

    if row["stage"] == "Oocyte" and os.path.isfile(new_f) and not os.path.isfile(old_f):
        print(f"  already re-filed as {NEW_ID} — nothing to do")
        return

    if os.path.isfile(old_f):
        sc = json.load(gzip.open(old_f, "rt"))
        sc["stage"] = "Oocyte"
        with gzip.open(new_f, "wt") as fh:
            json.dump(sc, fh, separators=(",", ":"))
        os.remove(old_f)
        print(f"  scene  {OLD_ID} → {NEW_ID}  (stage Zygote → Oocyte)")

    row["id"] = NEW_ID
    row["stage"] = "Oocyte"
    row["stage_label"] = "Oocyte"
    row["label"] = NEW_LABEL
    row["size_kb"] = round(os.path.getsize(new_f) / 1024)
    rows.sort(key=lambda r: (man["stages"].index(r["stage"]), r["id"]))
    with open(MAN, "w") as fh:
        json.dump(man, fh, indent=1)

    n_z = len(glob.glob(os.path.join(SEG, "Zygote__*.json.gz")))
    n_o = len(glob.glob(os.path.join(SEG, "Oocyte__*.json.gz")))
    print(f"  manifest row → stage Oocyte, label {NEW_LABEL}")
    print(f"  segments now: {n_z} zygote scenes, {n_o} oocyte scenes")
    print("\n  now re-run:  python3 build_size.py && python3 build_scheffler.py")


if __name__ == "__main__":
    main()
