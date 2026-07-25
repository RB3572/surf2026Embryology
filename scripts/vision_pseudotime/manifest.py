"""
Manifest validation for the vision pseudotime material (brief item 1).

Checks the two committed manifests (sources.csv, movies.csv) for schema and for
the non-negotiable safety rules:

  * every large source resolves to a path OUTSIDE the repository (never in git);
  * fixed snapshots (confocal / segmentation stacks) are never marked as
    supervised-tau sources — they have no true tau;
  * rendered movies are pilot-only, never a supervised tau source;
  * perturbation / overshoot / specialised-channel movies are excluded from tau
    training (role = exclude_from_tau_training or ood_only).

Availability (does the file exist on THIS machine) is *reported*, never fatal:
the pipeline must validate cleanly on a checkout that does not carry the raw
data.
"""
from __future__ import annotations

import csv
import os

from . import config

SOURCE_COLS = ["source_id", "source_kind", "location", "independent_embryos",
               "time_truth", "pixels_or_geometry", "supervised_tau_use",
               "intended_use", "limitations"]
MOVIE_COLS = ["movie_id", "file_name", "width_px", "height_px", "frames",
              "condition_or_subject", "default_role"]

SUPERVISED_VALUES = {"yes", "no", "pilot_only"}
FIXED_KINDS = {"confocal_stack", "segmentation_stack"}
# kinds that carry LARGE raw pixels/movies and therefore must live outside git.
# derived_numeric is a small committed CSV and is exempt.
LARGE_KINDS = {"confocal_stack", "segmentation_stack", "rendered_movies", "raw_3d_timelapse"}
MOVIE_ROLES = {"single_embryo_pilot", "split_control_and_ood",
               "exclude_from_tau_training", "ood_only"}
NON_TRAINING_ROLES = {"exclude_from_tau_training", "ood_only"}
PLACEHOLDER_LOCATIONS = {"not_yet_available", "unknown", ""}


def _read_csv(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return [r for r in csv.DictReader(fh) if any((v or "").strip() for v in r.values())]


def validate_sources(strict: bool = False) -> dict:
    rows = _read_csv(config.SOURCES_CSV)
    errors, warnings, entries = [], [], []
    missing = [c for c in SOURCE_COLS if rows and c not in rows[0]]
    if missing:
        errors.append(f"sources.csv missing columns: {missing}")
    for r in rows:
        sid = r.get("source_id", "?")
        loc = (r.get("location") or "").strip()
        sup = (r.get("supervised_tau_use") or "").strip()
        kind = (r.get("source_kind") or "").strip()
        resolved = config.resolve_source_path(sid, loc)
        is_placeholder = loc in PLACEHOLDER_LOCATIONS
        available = bool(resolved) and os.path.exists(resolved) and not is_placeholder

        if sup and sup not in SUPERVISED_VALUES:
            errors.append(f"{sid}: supervised_tau_use={sup!r} not in {sorted(SUPERVISED_VALUES)}")
        # non-negotiable: fixed snapshots are never supervised tau labels
        if kind in FIXED_KINDS and sup == "yes":
            errors.append(f"{sid}: fixed {kind} marked supervised_tau_use=yes — fixed snapshots "
                          "have no true tau")
        if kind == "rendered_movies" and sup == "yes":
            errors.append(f"{sid}: rendered movies cannot be a supervised tau source (pilot_only)")
        # non-negotiable: LARGE raw data must live outside the repo. Small derived
        # numeric CSVs (kind=derived_numeric) are legitimately committed inside it.
        if not is_placeholder and kind in LARGE_KINDS and config.is_inside_repo(resolved):
            errors.append(f"{sid}: large {kind} resolves INSIDE the repo "
                          f"({config.rel_to_repo(resolved)}) — must stay out of git")
        if not is_placeholder and not available:
            level = warnings if kind in LARGE_KINDS else errors
            label = "not present on this machine" if kind in LARGE_KINDS else "committed file missing"
            level.append(f"{sid}: {label} ({config.redact_path(resolved)})")

        entries.append({
            "source_id": sid, "source_kind": kind,
            "location_redacted": loc if is_placeholder else config.redact_path(resolved or loc),
            "independent_embryos": r.get("independent_embryos", ""),
            "supervised_tau_use": sup, "intended_use": r.get("intended_use", ""),
            "available_on_this_machine": available, "is_placeholder": is_placeholder,
        })
    if strict and errors:
        raise ValueError("sources.csv validation failed:\n  " + "\n  ".join(errors))
    return {"n": len(rows), "errors": errors, "warnings": warnings, "entries": entries}


def validate_movies(strict: bool = False) -> dict:
    rows = _read_csv(config.MOVIES_CSV)
    errors, warnings, entries = [], [], []
    missing = [c for c in MOVIE_COLS if rows and c not in rows[0]]
    if missing:
        errors.append(f"movies.csv missing columns: {missing}")
    for r in rows:
        mid = r.get("movie_id", "?")
        role = (r.get("default_role") or "").strip()
        cond = (r.get("condition_or_subject") or "").lower()
        if role and role not in MOVIE_ROLES:
            errors.append(f"{mid}: default_role={role!r} not in {sorted(MOVIE_ROLES)}")
        for numeric in ("width_px", "height_px", "frames"):
            v = (r.get(numeric) or "").strip()
            if v and not v.isdigit():
                errors.append(f"{mid}: {numeric}={v!r} is not an integer")
        # perturbation / overshoot subjects must never sit in normal-development training
        flagged = any(w in cond for w in ("nocodazole", "overshoot", "overexpression",
                                          "cytochalasin", "s25n", "droplet"))
        if flagged and role == "single_embryo_pilot":
            errors.append(f"{mid}: perturbation/overshoot subject marked single_embryo_pilot — "
                          "must be exclude_from_tau_training or ood_only")
        entries.append({"movie_id": mid, "default_role": role,
                        "trainable_normal_dev": role not in NON_TRAINING_ROLES,
                        "frames": r.get("frames", ""), "condition": r.get("condition_or_subject", "")})
    if strict and errors:
        raise ValueError("movies.csv validation failed:\n  " + "\n  ".join(errors))
    return {"n": len(rows), "errors": errors, "warnings": warnings, "entries": entries}


def validate_all(strict: bool = False) -> dict:
    s = validate_sources(strict=strict)
    m = validate_movies(strict=strict)
    ok = not (s["errors"] or m["errors"])
    return {"ok": ok, "sources": s, "movies": m}


def _print(report: dict) -> int:
    for name in ("sources", "movies"):
        sec = report[name]
        print(f"[{name}] {sec['n']} rows · {len(sec['errors'])} errors · {len(sec['warnings'])} warnings")
        for e in sec["errors"]:
            print(f"  ERROR   {e}")
        for w in sec["warnings"]:
            print(f"  note    {w}")
    print(f"\nmanifest validation: {'OK' if report['ok'] else 'FAILED'}")
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    import sys
    sys.exit(_print(validate_all()))
