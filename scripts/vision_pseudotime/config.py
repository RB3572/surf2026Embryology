"""
Local-data configuration and path handling (brief item 1).

Large TIFFs and movies live OUTSIDE the repository. This module resolves where
they are on a given machine and — critically — provides `redact_path`, so that
nothing written for the website or committed to git ever contains an absolute
local path or points at raw data.

Resolution order for a source's local path:
  1. environment variable  VPT_SRC_<SOURCE_ID>   (upper-case, non-alnum -> _)
  2. a gitignored  vision_local.json  in the repo root: {"source_id": "path"}
  3. the `location` column recorded in sources.csv

None of these are required to import the package; absence is reported, not fatal.
"""
from __future__ import annotations

import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))              # surf2026Embryology/
VISION_DIR = os.path.join(REPO_ROOT, "calibration_data", "vision_pseudotime")

SOURCES_CSV = os.path.join(VISION_DIR, "sources.csv")
MOVIES_CSV = os.path.join(VISION_DIR, "movies.csv")

# generated pixels (projections, thumbnails, extracted frames) — gitignored
DERIVED_DIR = os.path.join(VISION_DIR, "derived")
# tiny JSON/CSV safe to commit (provenance, metrics, redacted metadata)
ARTIFACTS_DIR = os.path.join(VISION_DIR, "artifacts")
# hand-authored sidecar metadata for stacks whose TIFF tags are insufficient
SIDECAR_DIR = os.path.join(VISION_DIR, "sidecars")

LOCAL_OVERRIDES = os.path.join(REPO_ROOT, "vision_local.json")


def _env_key(source_id: str) -> str:
    return "VPT_SRC_" + re.sub(r"[^A-Za-z0-9]", "_", source_id).upper()


def load_overrides() -> dict:
    if os.path.isfile(LOCAL_OVERRIDES):
        try:
            return dict(json.load(open(LOCAL_OVERRIDES)))
        except Exception:                                          # noqa: BLE001
            return {}
    return {}


def _abspath(p: str) -> str:
    """Expand ~ and resolve a relative path against the manifest dir (not CWD)."""
    if not p:
        return ""
    p = os.path.expanduser(p)
    if os.path.isabs(p):
        return p
    return os.path.normpath(os.path.join(VISION_DIR, p))


def resolve_source_path(source_id: str, location_fallback: str = "") -> str:
    """Best-effort absolute path for a source on THIS machine (may not exist).

    Relative locations are resolved against the manifest directory so validation
    does not depend on the current working directory.
    """
    env = os.environ.get(_env_key(source_id))
    if env:
        return _abspath(env)
    ov = load_overrides().get(source_id)
    if ov:
        return _abspath(ov)
    return _abspath(location_fallback or "")


def source_location(source_id: str) -> str:
    """The raw `location` recorded for a source in sources.csv (may be relative)."""
    import csv
    if not os.path.isfile(SOURCES_CSV):
        return ""
    with open(SOURCES_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("source_id") == source_id:
                return r.get("location", "")
    return ""


def resolve_source(source_id: str) -> str:
    """Full resolution: env / vision_local.json / sources.csv location."""
    return resolve_source_path(source_id, source_location(source_id))


def is_inside_repo(path: str) -> bool:
    """True if `path` resolves to somewhere inside the repository tree."""
    if not path:
        return False
    try:
        rp = os.path.realpath(path)
        root = os.path.realpath(REPO_ROOT) + os.sep
        return (rp + os.sep).startswith(root)
    except Exception:                                              # noqa: BLE001
        return False


def redact_path(path: str) -> str:
    """
    Turn an absolute local path into a safe token for committed / production
    output. Never emits a real filesystem path. A path that is a well-known
    non-local placeholder (e.g. 'not_yet_available') is passed through as-is.
    """
    if not path:
        return ""
    if "/" not in path and "\\" not in path:                       # already a token/placeholder
        return path
    base = os.path.basename(path.rstrip("/\\")) or "root"
    # label the neighbourhood without revealing it
    low = path.lower()
    if "icloud" in low or "mobile documents" in low:
        loc = "icloud-vault"
    elif low.startswith(("/users", "/home", os.path.expanduser("~").lower())):
        loc = "user-local"
    elif is_inside_repo(path):
        loc = "repo"
    else:
        loc = "local"
    return f"<{loc}>/{base}"


def rel_to_repo(path: str) -> str:
    """Repo-relative path for committed artifacts (safe: no absolute prefix)."""
    try:
        return os.path.relpath(path, REPO_ROOT)
    except Exception:                                              # noqa: BLE001
        return redact_path(path)


def ensure_dirs() -> None:
    for d in (DERIVED_DIR, ARTIFACTS_DIR, SIDECAR_DIR):
        os.makedirs(d, exist_ok=True)
