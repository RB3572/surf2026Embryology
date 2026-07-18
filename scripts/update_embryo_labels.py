#!/usr/bin/env python3
"""Apply canonical embryo display names to currently generated site data."""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from embryo_naming import embryo_label  # noqa: E402


def update_entry(entry: dict) -> bool:
    embryo_id = entry.get("eid") or entry.get("id")
    if not embryo_id:
        return False
    label = embryo_label(embryo_id, entry.get("stage") or entry.get("stage_label"))
    changed = False
    if "label" in entry and entry["label"] != label:
        entry["label"] = label
        changed = True
    if "title" in entry and entry["title"] != label:
        entry["title"] = label
        changed = True
    return changed


def update_json(path: Path) -> int:
    payload = json.loads(path.read_text())
    changed = sum(update_entry(entry) for entry in payload.get("embryos", []))
    if changed:
        path.write_text(json.dumps(payload, indent=1) + "\n")
    return changed


def update_gzip(path: Path, *, scene: bool = False) -> int:
    with gzip.open(path, "rt") as handle:
        payload = json.load(handle)
    if scene:
        label = embryo_label(payload["id"], payload.get("stage"))
        changed = int(payload.get("title") != label)
        if changed:
            payload["title"] = label
    else:
        changed = sum(update_entry(entry) for entry in payload.get("embryos", []))
    if changed:
        temporary = path.with_suffix(path.suffix + ".tmp")
        with gzip.open(temporary, "wt") as handle:
            json.dump(payload, handle, separators=(",", ":"))
        temporary.replace(path)
    return changed


def main() -> None:
    changed = 0
    for relative in (
        "data/manifest.json",
        "data/segments_manifest.json",
        "data/axes_manifest.json",
        "data/pronuclei_manifest.json",
        "data/zygote_manifest.json",
    ):
        changed += update_json(ROOT / relative)
    for relative in (
        "data/analysis_index.json.gz",
        "data/zygote_cross.json.gz",
        "data/zygote_cross_circ.json.gz",
    ):
        changed += update_gzip(ROOT / relative)
    for path in (ROOT / "data" / "scenes").glob("*.json.gz"):
        changed += update_gzip(path, scene=True)
    print(f"Updated {changed} generated embryo labels/titles")


if __name__ == "__main__":
    main()
