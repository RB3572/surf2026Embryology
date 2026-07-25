# Model card — fixed-image pronuclear pseudotime (PILOT)

**Package:** `scripts/vision_pseudotime/` (`vpt-0.1.0-pilot`)
**Status:** pilot + data pipeline. **No end-to-end image clock is validated.**
The only deployable predictor here is the pre-existing geometry clock, reused
unchanged through a hybrid path.

This card separates what is **validated** from what is **pilot-only**, and lists
the data required before any publication-grade image-clock claim.

---

## 1. Validated / deployable

### 1a. Geometry-to-τ clock (baseline) — `pnpt-3.0.0`
Unchanged from the existing calibration; this pilot does not modify it.

- Input: symmetric pronucleus-to-cell-centre distances (µm) — identity-free.
- Trained on 53 untreated live-imaged embryos (Scheffler 2021), 2,057 frames.
- **Nested, embryo-grouped CV (untouched outer test):**
  macro MAE **0.091** (95% CI 0.081–0.102), Spearman ρ **0.899**,
  strict pair-order **86.4%**, empirical interval half-width **±0.242 τ**
  (empirical, not a coverage guarantee on the fixed cohort).
- Selected model: isotonic regression on the distance sum.

### 1b. Hybrid path: image → segmentation → geometry → frozen clock
The deployable route for a *segmented* fixed embryo. It reuses `pnpt-3.0.0`
verbatim (`build_pronuclei_pseudotime.predict` / `qc_status`) and only adds
geometry extraction from a label volume.

- **Verified** to reproduce the deployed pipeline exactly on real MERFISH
  zygotes: max feature error **0.0000 µm**, max **|Δτ| = 0.0** vs the published
  `data/pronuclei_pseudotime.json` (cross-check in `report.hybrid_section`,
  re-checked in the test suite over the whole fixed cohort within rounding).
- Emits τ **only** with verified micron spacing **and** verified label-class
  semantics; otherwise it returns a `blocked` result naming the gaps.

**What 1a/1b support:** giving a *segmented* fixed zygote a calibrated τ with the
clock's own uncertainty and QC. Nothing here is a direct pixels→τ model.

---

## 2. Pilot-only (NOT validated)

### 2a. Metadata audit + streaming 2.5D preprocessing
- Audits a 2.66 GB stack from tags alone in ~0.3 s (no pixel load) and builds a
  256² 3-channel 2.5D representation (normalized MIP, robust nonzero mean,
  segmentation occupancy) by streaming one z-page at a time — full pass ~13 s,
  peak RSS ~0.47 GB, never the 2.66 GB full load. No raw sum projection.
- **On the example lab stack this is pilot-only:** voxel spacing is unverified
  and label-class semantics are unverified, so projections are valid as pixels
  but carry no physical scale, and the hybrid clock is correctly **blocked**
  (that stack also has only 2 label IDs, not the cell + 2-pronuclei structure the
  clock needs).

### 2b. Movie extraction
- Extracts frames with overlays stripped (deterministic border-crop + corner
  mask, **not** OCR) and records source/panel/frame/treatment/channel/inclusion.
- τ is **never** inferred from frame index; it stays null without annotated
  pronuclear-formation and NEBD frames. Perturbation/overshoot movies are OOD or
  excluded. Movie 1 is one embryo → cannot validate anything.

### 2c. Exploratory direct image encoder — **feasibility only**
- A deterministic handcrafted-feature (rotation-invariant radial) ridge model.
- **Demonstrated only on labelled-synthetic phantoms** (held-out-embryo grouped
  MAE ≈ 0.025, Spearman ≈ 0.996 on synthetic), where it beats constant /
  brightness / corner-overlay baselines and — after overlay stripping — shows
  near-zero saliency on the timestamp corner.
- **On the real material it refuses to train** (one embryo, no per-frame τ). No
  performance number on real data is claimed, because none can be.

---

## 3. Non-negotiables honoured

- Augmentation regularizes; it does not create independent embryos. Augmented
  views inherit their parent embryo's split (leakage-tested).
- Fixed DAPI/all-stain snapshots have no true τ → segmentation / preprocessing /
  OOD only, never supervised labels.
- Public rendered movies are not mapped to the 53 numeric trajectories.
- Deterministic seeds throughout; model/data versions recorded in artifacts.
- Large TIFFs/movies stay outside git; the committed report carries no absolute
  paths and no raw pixels.

---

## 4. Data required before publication-grade image-clock claims

See `ACQUISITION_CHECKLIST.md`. In short:

1. **Multiple independent untreated RAW time-lapse embryos** (~8–10+), each with
   real frame times, verified xy pixel size + z spacing, treatment, embryo ID,
   and annotated pronuclear-formation + NEBD frames.
2. **Verified voxel spacing** and **verified label-class semantics** for every
   stack/segmentation (sidecars).
3. A **held-out embryo/batch** test with **frozen preprocessing**.
4. Masks/annotations where available (for supervised segmentation + Dice/IoU /
   centroid error).

Until then, the deployable answer is the geometry clock (1a) applied through the
hybrid path (1b); the direct image encoder remains a labelled feasibility study.
