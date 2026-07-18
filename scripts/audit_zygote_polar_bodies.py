#!/usr/bin/env python3
"""Validate geometry-selected polar bodies in generated zygote analyses."""

from __future__ import annotations

import gzip
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "zygote"
MANIFEST_PATH = ROOT / "data" / "zygote_manifest.json"

# Regression case that exposed the original hard-coded-label defect.
KNOWN_POLAR_BODY_LABELS = {"20260407_zygote_p1_12": 4}


def fail(message: str) -> None:
    raise AssertionError(message)


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text())
    entries = manifest.get("embryos", [])
    expected_ids = {entry["id"] for entry in entries}
    scene_paths = sorted(DATA_DIR.glob("*.json.gz"))
    scene_ids = {path.name[:-8] for path in scene_paths}
    if scene_ids != expected_ids:
        fail(
            f"manifest/scene mismatch: missing={sorted(expected_ids - scene_ids)}, "
            f"stale={sorted(scene_ids - expected_ids)}"
        )

    selected_labels = Counter()
    for path in scene_paths:
        with gzip.open(path, "rt") as handle:
            scene = json.load(handle)
        analysis = scene.get("analysis", {})
        label = analysis.get("polar_body_label")
        diagnostics = analysis.get("polar_body_detection", {})
        candidates = diagnostics.get("candidates", [])
        selected = next((candidate for candidate in candidates if candidate.get("label") == label), None)
        external = [candidate for candidate in candidates if candidate.get("external")]
        if selected is None:
            fail(f"{scene['id']}: selected polar-body label {label} has no diagnostics")
        if not selected.get("external"):
            fail(f"{scene['id']}: selected polar-body label {label} is not external")
        if not external:
            fail(f"{scene['id']}: no external polar-body candidate")
        expected = max(
            external,
            key=lambda candidate: (
                candidate["radial_ratio"],
                candidate["outside_fraction"],
                candidate["max_outside_um"],
                candidate["voxel_count"],
            ),
        )
        if expected["label"] != label:
            fail(f"{scene['id']}: selected label {label}, expected largest external label {expected['label']}")
        if len(analysis.get("pb_plot", [])) != 3:
            fail(f"{scene['id']}: missing polar-body coordinates")
        if any(not 0 <= int(index) < 18 for index in analysis.get("best_planes", {}).values()):
            fail(f"{scene['id']}: best-plane index outside 0..17")
        selected_labels[int(label)] += 1

    by_id = {entry["id"]: entry for entry in entries}
    for embryo_id, expected_label in KNOWN_POLAR_BODY_LABELS.items():
        actual_label = by_id.get(embryo_id, {}).get("polar_body_label")
        if actual_label != expected_label:
            fail(f"{embryo_id}: regression, expected label {expected_label}, found {actual_label}")

    print(f"PASS: {len(scene_paths)} zygotes have externally validated polar bodies")
    print(f"Selected label distribution: {dict(sorted(selected_labels.items()))}")


if __name__ == "__main__":
    main()
