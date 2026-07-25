# Vision pseudotime dataset plan

This directory catalogs the material currently available for a fixed-image
pronuclear pseudotime model. Large microscopy and movie files stay outside Git;
the manifests record provenance and intended use.

## Scientific target

- Input: one fixed zygote snapshot, represented as a confocal z-stack or a
  deterministic 2D/2.5D projection.
- Output: normalized pronuclear-stage pseudotime `tau`, where `tau = 0` is
  pronuclear formation and `tau = 1` is nuclear-envelope breakdown (NEBD), plus
  an uncertainty interval and an out-of-distribution warning.
- The fixed MERFISH snapshots have no true collection time. They are deployment
  data, not supervised labels.

## Available evidence

1. `../scheffler2021/scheffler_2021_control_zygote_trajectories.csv`
   contains 2,057 time-labelled rows from 53 untreated embryos. It has
   pronucleus-to-cell-centre distances and `tau`, but no pixels, masks, 3D
   coordinates, cell diameter, or frame-to-image mapping.
2. The 20 public Scheffler supplementary movies are rendered MJPEG presentation
   movies. `movies.csv` records their metadata and biological conditions. Movie
   1 is the only single-panel untreated full pronuclear-migration sequence. Most
   other movies are comparisons, perturbations, specialised channels, or
   overshoot examples.
3. `/Users/rishib/Desktop/all.tif` and
   `/Users/rishib/Desktop/Segmentation-label_1.tif` are matched 251-slice,
   2303 x 2304 stacks. The first is `uint16` intensity and the second is `int16`
   labels with IDs 1 and 2. They do not contain reliable embedded voxel spacing,
   channel identity, embryo ID, or acquisition time. A sidecar record is required.
4. Raw Scheffler microscopy stacks are not publicly deposited; the paper states
   that they are available from the corresponding author on request.

## What can and cannot be trained now

The current material supports a **pilot**, not an independently validated
end-to-end image clock. Augmented views of one movie remain one embryo and may
not be counted as independent samples. The deployable baseline remains the
existing geometry-to-`tau` model trained on 53 embryos.

The recommended model is therefore staged:

1. Reproduce and freeze the existing geometry baseline.
2. Build target-domain segmentation and projection preprocessing from the
   laboratory stacks and masks.
3. Train a hybrid image -> cell/pronucleus masks -> symmetric geometry -> `tau`
   pipeline. This can use the 53-embryo numeric clock even before raw public
   pixels are obtained.
4. Build an exploratory direct image encoder only from frames with defensible
   biological timestamps. Label it as a feasibility model until multiple raw,
   independent untreated time-lapse embryos are available.
5. Retrain and evaluate the image model after obtaining the raw time-lapse
   stacks and acquisition metadata.

## Input representation

Do not sum all z-slices: it saturates signal and preserves background artifacts.
Crop to the embryo first. Prefer a 2.5D tensor with separately normalized
channels:

1. intensity maximum projection;
2. robust mean or percentile projection;
3. segmentation occupancy/depth or boundary projection.

A single normalized maximum-intensity projection is a baseline only. Preserve
the full stack and voxel metadata so a 3D model remains possible.

## Augmentation rules

Allowed defaults: rotations, flips, translations that preserve the embryo,
mild isotropic scaling, intensity/gamma variation, noise, blur, simulated
photobleaching, z-slice dropout, and small projection-window changes.

Disabled by default: anisotropic stretch, shear, elastic warping, independent
pronuclear movement, or any transform that changes the geometry encoding time.
Such transforms require a biological justification and recomputed targets.

Split by embryo/source **before** augmentation. Every frame and every augmented
copy from one embryo must remain in one split. Adjacent frames are not
independent observations.

## Leakage controls

- Crop timestamps, scale bars, panel labels, and borders before model input.
- Parse timestamp text only into a separate label table; never leave it in the
  training pixels.
- Do not derive `tau` as `frame_index / frame_count` unless pronuclear formation
  and NEBD have been independently annotated.
- Do not use future-normalized pronuclear volumes from the source workbook.
- Keep perturbation movies out of normal-development training; use them for
  out-of-distribution tests unless a specific control panel is isolated.
- Compare against trivial baselines and inspect saliency to verify that the
  embryo, not presentation artifacts, drives predictions.

## Validation

Report embryo-grouped MAE, Spearman ordering, strict pair-order accuracy,
interval coverage, and OOD performance. Also report segmentation IoU/Dice,
pronuclear-centroid error, and geometry error. A split of frames from a single
movie is useful for debugging but is not an external test.

The first production image model requires multiple independent untreated raw
time-lapse embryos, a held-out embryo/batch test, and frozen preprocessing.

