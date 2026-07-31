#!/usr/bin/env python3
"""Regression tests for the four-way embryo probeset display names."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from embryo_naming import embryo_label, probeset_for  # noqa: E402


class EmbryoNamingTests(unittest.TestCase):
    def test_mapping_has_only_the_four_real_probesets(self) -> None:
        mapping = json.loads((ROOT / "data" / "probesets.json").read_text())
        self.assertEqual(len(mapping), 158)
        self.assertEqual(set(mapping.values()), {"1_0", "1_1", "2_0", "2_1"})

    def test_representative_labels(self) -> None:
        cases = {
            "20251226_zygote_p1_0": "Z-P1_1-fov0",
            "20251226_zygote_p1_2_1": "Z-P1_0-fov2_1",
            "20251226_zygote_p1_2_2": "Z-P1_1-fov2_2",
            "20260408_e2c_p2_15_1": "e2c-P2_1-fov15_1",
            "20260420_l2c_p1_3_1": "l2c-P1_1-fov3_1",
            "20260425_oocyte_p2_11": "O-P2_1-fov11",
            "20260425_zygote_p2_1": "Z-P2_0-fov1",
            "20251226_sample1_zygote5": "Z-P1_1-fov5",
        }
        for embryo_id, expected in cases.items():
            with self.subTest(embryo_id=embryo_id):
                self.assertEqual(embryo_label(embryo_id), expected)

    def test_stage_namespace_and_override_are_preserved(self) -> None:
        self.assertEqual(
            embryo_label("Oocyte__20260425_zygote_p2_3", "Oocyte"),
            "O-P2_0-fov3",
        )

    def test_unknown_embryo_is_not_given_a_guessed_probeset(self) -> None:
        embryo_id = "20990101_zygote_p1_0"
        self.assertIsNone(probeset_for(embryo_id))
        self.assertEqual(embryo_label(embryo_id), embryo_id)


if __name__ == "__main__":
    unittest.main()
