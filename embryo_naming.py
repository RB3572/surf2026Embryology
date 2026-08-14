"""Canonical display names for embryo records and figures."""

from __future__ import annotations

import json
import re
from pathlib import Path


_PROBESET_PATH = Path(__file__).resolve().parent / "data" / "probesets.json"
_PROBESETS: dict[str, str] = json.loads(_PROBESET_PATH.read_text())

# The canonical embryo IDs, generated from CompleteEmbryoDataset.xlsx. Both
# this site and the MERFISH atlas LOOK THE ID UP here rather than deriving it,
# so neither can drift from the workbook or from each other -- which is what
# happened with the three incompatible schemes this replaces.
#
#   Z-251226-PL02-PSB1_1-FOV07
#   stage (verified) - date - plate - probeset - FOV
#
# The probeset carries a letter AND its p<n>_<m> name because `p1` in the old
# filenames is the PLATE, and every plate splits across two probesets.
_CANONICAL_PATH = Path(__file__).resolve().parent / "data" / "embryo_ids.json"
try:
    _CANONICAL: dict[str, dict] = json.loads(
        _CANONICAL_PATH.read_text()).get("embryos", {})
except (OSError, ValueError):
    _CANONICAL = {}


def canonical_id(embryo_id: str) -> str | None:
    """The canonical UID for an embryo, or None if it is not in the lookup."""
    return (_CANONICAL.get(raw_embryo_id(embryo_id)) or {}).get("uid")


_STAGE_PREFIX = {
    "o": "O",
    "oocyte": "O",
    "z": "Z",
    "zygote": "Z",
    "e2c": "e2c",
    "early": "e2c",
    "early2cell": "e2c",
    "l2c": "l2c",
    "late": "l2c",
    "late2cell": "l2c",
}


def raw_embryo_id(embryo_id: str) -> str:
    """Strip a stage namespace such as ``Zygote__`` from an embryo ID."""
    return str(embryo_id).split("__", 1)[-1]


def probeset_for(embryo_id: str) -> str | None:
    """Return the displayed probeset assignment (1, 2, 3, or 4)."""
    return _PROBESETS.get(raw_embryo_id(embryo_id))


def embryo_label(embryo_id: str, stage: str | None = None) -> str:
    """The canonical UID, e.g. ``Z-251226-PL02-PSB1_1-FOV07``.

    Falls back to the old ``TYPE-PROBESET-fovN`` form only for an embryo that
    is not in ``data/embryo_ids.json`` -- and, as before, to the raw ID when
    even the probeset is unknown, which is deliberately conspicuous.
    """
    raw_id = raw_embryo_id(embryo_id)
    uid = canonical_id(raw_id)
    if uid:
        return uid
    patterns = (
        r"^\d{8}_(?P<stage>oocyte|zygote|e2c|l2c|early2cell|late2cell)_p(?P<probe>\d+)_(?P<fov>.+)$",
        r"^\d{8}_(?P<stage>l2c)_blastomere_p(?P<probe>\d+)_(?P<fov>.+)$",
        r"^\d{8}_sample(?P<probe>\d+)_(?P<stage>zygote)(?P<fov>\d+(?:_\d+)?)$",
    )
    match = next(
        (candidate for pattern in patterns if (candidate := re.match(pattern, raw_id, flags=re.IGNORECASE))),
        None,
    )
    if match is None:
        return embryo_id

    id_stage = match.group("stage")
    probeset = probeset_for(raw_id)
    if probeset is None:
        # Returning the ID is intentionally conspicuous: a new embryo must be added to
        # data/probesets.json rather than silently displaying its old plate number.
        return embryo_id
    fov = match.group("fov")
    normalized_stage = re.sub(r"[^a-z0-9]", "", (stage or "").lower())
    # Generic two-cell stage labels cannot distinguish early from late; the ID can.
    prefix = _STAGE_PREFIX.get(normalized_stage) or _STAGE_PREFIX[id_stage.lower()]
    return f"{prefix}-P{probeset}-fov{fov}"
