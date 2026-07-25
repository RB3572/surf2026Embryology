"""
Embryo/source-grouped splits BEFORE augmentation (brief item 6).

The split a sample lands in is a PURE FUNCTION of its biological group key
(embryo_id or source_id) — computed before any augmentation exists. Because
every augmented derivative inherits its parent's group key, it is structurally
impossible for a frame or an augmented copy of one embryo to appear in more than
one split. Adjacent frames of one embryo are not independent, so they share a
group and therefore a split.

Assignment is deterministic via a stable hash (hashlib, NOT Python's salted
hash()), so it reproduces across processes and machines.
"""
from __future__ import annotations

import hashlib

SPLITS = ("train", "val", "test")


def group_key(sample: dict) -> str:
    """The biological unit that must never be split. Prefers embryo_id, then source."""
    for k in ("embryo_id", "group", "source_id", "movie_id"):
        v = sample.get(k)
        if v:
            return str(v)
    raise KeyError(f"sample has no biological group key: {sample!r}")


def _unit(key: str, seed: int) -> float:
    """Stable [0,1) hash of (seed, key)."""
    h = hashlib.sha256(f"{seed}:{key}".encode()).hexdigest()
    return int(h[:16], 16) / float(1 << 64)


def assign_split(key: str, ratios=(0.6, 0.2, 0.2), seed: int = 0) -> str:
    """Deterministically map a group key to a split by its stable hash."""
    total = float(sum(ratios))
    r = _unit(key, seed) * total
    c = 0.0
    for name, w in zip(SPLITS, ratios):
        c += w
        if r < c:
            return name
    return SPLITS[-1]


def split_groups(group_keys, ratios=(0.6, 0.2, 0.2), seed: int = 0) -> dict:
    """{group_key: split} for a set of biological groups."""
    return {k: assign_split(k, ratios, seed) for k in sorted(set(group_keys))}


def assign_samples(samples: list[dict], ratios=(0.6, 0.2, 0.2), seed: int = 0) -> list[dict]:
    """Attach a 'split' to each sample from its group key (pre-augmentation)."""
    out = []
    for s in samples:
        gk = group_key(s)
        out.append({**s, "group_key": gk, "split": assign_split(gk, ratios, seed)})
    return out


def expand_with_augmentation(assigned: list[dict], n_aug: int, aug_seed: int = 0) -> list[dict]:
    """
    Produce augmented sample records AFTER the split is fixed. Each derivative
    inherits its parent's group_key and split; only a per-copy augmentation seed
    changes. This ordering is the whole point: augmentation cannot move a sample.
    """
    out = []
    for s in assigned:
        base = {**s, "aug_index": 0, "is_augmented": False, "aug_seed": None}
        out.append(base)
        for j in range(1, n_aug + 1):
            seed = int(_unit(f"{s['group_key']}:{s.get('frame', 0)}:{s.get('panel', 0)}:{j}",
                             aug_seed) * (1 << 32))
            out.append({**s, "aug_index": j, "is_augmented": True, "aug_seed": seed})
    return out


def check_no_leakage(records: list[dict]) -> dict:
    """
    Verify the non-negotiable: no group in >1 split, and every augmented
    derivative shares its parent's split. Returns {ok, violations, ...}.
    """
    group_to_splits: dict[str, set] = {}
    parent_split: dict = {}
    aug_violations = []
    for r in records:
        gk = r.get("group_key") or group_key(r)
        group_to_splits.setdefault(gk, set()).add(r["split"])
        if not r.get("is_augmented"):
            parent_split[(gk, r.get("frame", 0), r.get("panel", 0))] = r["split"]
    for r in records:
        if r.get("is_augmented"):
            key = (r["group_key"], r.get("frame", 0), r.get("panel", 0))
            ps = parent_split.get(key)
            if ps is not None and ps != r["split"]:
                aug_violations.append({**{k: r[k] for k in ("group_key", "frame", "panel", "split")},
                                       "parent_split": ps})
    cross = {g: sorted(s) for g, s in group_to_splits.items() if len(s) > 1}
    return {"ok": not cross and not aug_violations,
            "groups": len(group_to_splits),
            "cross_split_groups": cross,
            "augmented_split_mismatches": aug_violations,
            "per_split_group_counts": _counts(group_to_splits)}


def _counts(group_to_splits: dict) -> dict:
    c = {s: 0 for s in SPLITS}
    for splits in group_to_splits.values():
        for s in splits:
            c[s] = c.get(s, 0) + 1
    return c
