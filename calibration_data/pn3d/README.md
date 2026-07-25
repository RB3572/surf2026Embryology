# 3D Pronuclear Pseudotime Calibration

Estimate a calibrated developmental pseudotime **τ ∈ [0,1]** (0 = pronuclear
formation, 1 = NEBD) for a **static 3D zygote stack**, as a distribution with
uncertainty and structural QC — never an unjustified point estimate.

```
image → audited 3D segmentation → dimensionless geometry → calibrated probabilistic clock → τ ± interval (+ QC / OOD)
```

## Data (mounted, out of git)

Authoritative source: `/Volumes/HW/MEFISH Labels` (~733 GB, never copied into the
repo). Inventoried into a **versioned manifest** (`data/pn3d/manifest.json`,
redacted paths):

- **8 experiments / acquisition batches**, **121 embryos** — 54 zygote, 42 early-
  2-cell, 25 late-2-cell.
- Per embryo: `dapi.tif` (DNA), `all.tif` (all channels), a 3D Slicer
  `Segmentation-label.tif` + `Segmentation.seg.nrrd`, generic segment names.
- Voxel grid anisotropy z = 7 × xy (from the NRRD); physical xy is a documented
  display constant — **all model features are dimensionless** (÷ cell radius).
- **No time-lapse, no PNF/NEBD annotations** anywhere in this dataset. It is the
  *target imaging domain*. True developmental time comes from the external
  **Scheffler 2021 live-imaging** cohort (53 embryos, 2057 frames, real τ).

## Pipeline

1. **Manifest** (`scripts/build_pn3d_manifest.py`) — inventory, provenance,
   spacing, completeness, duplicates/shared-FOV, batch.
2. **Segmentation audit** (`scripts/pn3d/segment_audit.py`) — classify the generic
   segments **geometrically** (never by name) into cell body / pronucleus 1 /
   pronucleus 2 / polar body. Enforce: exactly two compact pronuclei *inside* the
   cell body; the polar body *external* and never substituted for a pronucleus.
   Fill-containment distinguishes a boundary-touching pronucleus from a detached
   polar body. Unresolvable stacks are **marked, not forced** (40/53 zygotes
   resolve; the rest are honestly unresolved).
3. **Dimensionless geometry** — Σ(pronucleus→cell-centre)/radius, inter-pronuclear
   /radius, volume asymmetry, etc.
4. **Probabilistic clock** (`scripts/pn3d/clock.py`) — monotone isotonic mean +
   heteroscedastic, **calibrated** interval, fit on live-imaging τ.
5. **Inference + OOD** (`scripts/pn3d/inference.py`) — per-embryo τ posterior;
   flag stacks outside live support / low-confidence / wrong stage.

## Validated results (leave-one-embryo-out CV on live-imaging τ)

The **only** independent time validation:

| metric | value |
|---|---|
| MAE | **0.090** |
| Spearman ρ | **0.898** |
| within-embryo monotonicity (median ρ) | **0.998** |
| interval coverage 50 / 80 / 95 % | **0.46 / 0.79 / 0.95** (calibrated) |

Calibration is near-diagonal at every level (0.5→0.46 … 0.95→0.95).

## Evidence table (what is and isn't established)

| claim | evidence type | result |
|---|---|---|
| geometry→τ clock generalizes | **TRUE time-supervised** (live CV) | MAE 0.090, 95% coverage 0.95 |
| segmentation resolves the 4 structures | structural | 40/53 resolved with constraints; rest unresolved, not forced |
| fixed domain handled | domain adaptation (partial) | scale-invariant features + OOD gating; residual gap is stage-sampling |
| image-only τ from raw pixels | **exploratory (not reliable)** | ρ≈0.3, batch-confounded — needs a learned 3D segmenter |
| per-embryo τ for fixed MERFISH | **reference output** | NOT independent validation; fixed stacks have no true time |

## Reproduce

```bash
python3 scripts/build_pn3d_manifest.py      # inventory (metadata only, ~1 s)
python3 scripts/build_pn3d_geometry.py      # segmentation audit + geometry (~7 min)
python3 scripts/build_pn3d_model.py         # clock + inference + evidence
python3 scripts/build_pn3d_previews.py      # orthogonal-view previews (gitignored)
python3 scripts/test_pn3d.py                # 35 tests
```

Website: **`pronuclear-pseudotime.html`** — an Explorer (per-embryo 3D XY/XZ/YZ
with segmentation overlays, geometry, τ posterior, QC/OOD, τ 0→1 animation,
JSON/CSV download) and a Model-development view (provenance, clock validation +
calibration, baselines, domain shift, evidence table, failure cases).

## Limitations & next data

- **No independent time labels** on the fixed cohort — τ is inferred, validated
  only on live imaging; agreement with prior cached pseudotime is **not**
  validation.
- The segmentation layer **uses the provided Slicer masks** (audited + QC'd). A
  learned volumetric segmenter (needs GPU/torch, absent here) is future work;
  classical voxel-RF segmentation was tried and failed (Dice ≈ 0.08).
- **Direct image→τ** from raw pixels does not work with available tools
  (ρ≈0.3, batch-confounded) — a genuine negative result, kept as the image-only
  baseline.
- To lift these: multiple **independent raw untreated time-lapse** embryos with
  frame times + PNF/NEBD frames + verified voxel spacing, and per-stack live cell
  radius to remove the reference-radius assumption.
