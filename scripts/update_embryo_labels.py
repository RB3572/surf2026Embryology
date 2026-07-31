#!/usr/bin/env python3
"""Apply the authoritative four-way probeset names to all generated site data."""

from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from embryo_naming import probeset_for, raw_embryo_id  # noqa: E402


DISPLAY_NAME = re.compile(r"^(O|Z|e2c|l2c)-P(?:1|2)(?:_[01])?-fov")
ID_KEYS = ("eid", "tx_id", "embryo_id", "id")
JS_MAPPING = ROOT / "embryo-probesets.js"
LABEL_MARKERS = (b"-P1", b"-P2")


def mapped_id(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    candidate = raw_embryo_id(value)
    return candidate if probeset_for(candidate) else None


def replace_display_name(value: str, embryo_id: str) -> str:
    probeset = probeset_for(embryo_id)
    if not probeset:
        return value
    return DISPLAY_NAME.sub(rf"\1-P{probeset}-fov", value, count=1)


def update_node(node: object, inherited_id: str | None = None) -> list[tuple[str, str]]:
    """Return ordered JSON-string replacements while walking the parsed payload."""
    replacements: list[tuple[str, str]] = []
    if isinstance(node, dict):
        local_id = next((candidate for key in ID_KEYS if (candidate := mapped_id(node.get(key)))), inherited_id)
        for key, value in node.items():
            if key in ("label", "title") and local_id:
                old = value
                if isinstance(old, str) and DISPLAY_NAME.match(old):
                    new = replace_display_name(old, local_id)
                    if new != old:
                        node[key] = new
                        replacements.append((old, new))
                continue
            child_id = mapped_id(key) or local_id
            replacements.extend(update_node(value, child_id))
    elif isinstance(node, list):
        for value in node:
            replacements.extend(update_node(value, inherited_id))
    return replacements


def apply_replacements(text: str, replacements: list[tuple[str, str]], path: Path) -> str:
    """Patch only changed JSON string values so unrelated formatting stays byte-for-byte stable."""
    index = 0
    for old, new in replacements:
        token = json.dumps(old)
        found = text.find(token, index)
        if found < 0:
            raise RuntimeError(f"Could not locate {token} in {path}")
        replacement = json.dumps(new)
        text = text[:found] + replacement + text[found + len(token):]
        index = found + len(replacement)
    return text


def update_json(path: Path) -> int:
    raw = path.read_bytes()
    if not any(marker in raw for marker in LABEL_MARKERS):
        return 0
    text = raw.decode()
    payload = json.loads(text)
    replacements = update_node(payload)
    if replacements:
        path.write_text(apply_replacements(text, replacements, path))
    return len(replacements)


def update_gzip(path: Path) -> int:
    raw = gzip.decompress(path.read_bytes())
    if not any(marker in raw for marker in LABEL_MARKERS):
        return 0
    if path.parent.name == "scenes":
        embryo_id = path.name.removesuffix(".json.gz")
        probeset = probeset_for(embryo_id)
        if not probeset:
            raise RuntimeError(f"No probeset assignment for scene {embryo_id}")
        pattern = re.compile(rb'("title"\s*:\s*"(?:O|Z|e2c|l2c)-P)(?:1|2)(?:_[01])?(-fov)')
        updated, count = pattern.subn(
            lambda match: match.group(1) + probeset.encode() + match.group(2), raw, count=1
        )
        if count != 1:
            raise RuntimeError(f"Could not update scene title in {path}")
        if updated == raw:
            return 0
        path.write_bytes(gzip.compress(updated, compresslevel=9, mtime=0))
        return count
    text = raw.decode()
    payload = json.loads(text)
    replacements = update_node(payload)
    if replacements:
        updated = apply_replacements(text, replacements, path).encode()
        path.write_bytes(gzip.compress(updated, compresslevel=9, mtime=0))
    return len(replacements)


def update_csv(path: Path) -> int:
    lines = path.read_text().splitlines(keepends=True)
    changed = 0
    for index, line in enumerate(lines[1:], start=1):
        columns = line.split(",", 2)
        if len(columns) < 2:
            continue
        embryo_id, old = columns[0], columns[1]
        if not mapped_id(embryo_id) or not DISPLAY_NAME.match(old):
            continue
        new = replace_display_name(old, embryo_id)
        if new != old:
            columns[1] = new
            lines[index] = ",".join(columns)
            changed += 1
    if changed:
        path.write_text("".join(lines))
    return changed


def write_js_mapping() -> None:
    probesets = json.loads((ROOT / "data" / "probesets.json").read_text())
    body = json.dumps(probesets, indent=2, sort_keys=True)
    JS_MAPPING.write_text(
        "// Generated from data/probesets.json by scripts/update_embryo_labels.py.\n"
        f"window.EMBRYO_PROBESETS = Object.freeze({body});\n"
    )


def main() -> None:
    changed = 0
    files_changed = 0
    for path in sorted((ROOT / "data").rglob("*.json")):
        count = update_json(path)
        changed += count
        files_changed += bool(count)
    for path in sorted((ROOT / "data").rglob("*.json.gz")):
        count = update_gzip(path)
        changed += count
        files_changed += bool(count)
    for path in sorted((ROOT / "calibration_data").rglob("*.csv")):
        count = update_csv(path)
        changed += count
        files_changed += bool(count)
    write_js_mapping()
    print(f"Updated {changed} embryo labels/titles across {files_changed} data files")


if __name__ == "__main__":
    main()
