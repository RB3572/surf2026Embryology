# Implementation brief: fixed-image pronuclear pseudotime

Work in `/Users/rishib/Desktop/SpermLabeling/surf2026Embryology`.

Build a reproducible **pilot and data pipeline**, not an overstated production
vision model. Read `calibration_data/vision_pseudotime/README.md`, `sources.csv`,
and `movies.csv` first, then inspect the existing numeric calibration code and
tests before changing anything.

## Required work

1. Add a local-data configuration and manifest validator. Large TIFFs and movies
   must remain outside Git. Never copy them into the repository or frontend.
2. Add a metadata-only TIFF audit and incremental projection tool that never
   loads a 2.5 GB volume in full. Support intensity stacks and integer label
   stacks. Require sidecar voxel/channel/class metadata when TIFF tags are
   insufficient.
3. Implement embryo cropping and a deterministic 2.5D representation:
   normalized maximum-intensity projection, robust nonzero mean/percentile
   projection, and segmentation occupancy/boundary projection. Also emit a
   single-channel MIP baseline. Do not use a raw sum projection.
4. Add movie extraction that records source movie, panel, frame, visible time,
   treatment, channel, and inclusion status. Remove timestamps, scale bars,
   labels, and borders from model pixels. Store parsed/manual timestamps only in
   metadata. Do not infer `tau` from frame index without annotated pronuclear
   formation and NEBD.
5. Build augmentation code with deterministic seeds. Default transforms may use
   rotation, reflection, translation, mild isotropic scaling, intensity/gamma,
   noise, blur, simulated bleaching, z-dropout, and projection-window changes.
   Do not enable anisotropic stretch, shear, elastic deformation, or independent
   pronuclear motion because they alter the time-bearing geometry.
6. Split by biological embryo/source before augmentation. Add tests proving that
   no frame or augmented derivative from one embryo occurs in multiple splits.
7. Preserve the current geometry-to-`tau` clock as the primary baseline. Add a
   hybrid path: image -> segmentation -> symmetric geometric features -> frozen
   clock. Report uncertainty and QC failures.
8. Add a clearly labelled exploratory direct image encoder only for frames with
   defensible timestamps. Movie 1 is one embryo and cannot provide independent
   validation. Use any isolated control panels only after condition and stage
   boundaries are verified. Keep all perturbation/overshoot examples out of
   normal-development training and use them as OOD cases.
9. Evaluate embryo-grouped MAE, Spearman ordering, pair-order accuracy, interval
   coverage, segmentation Dice/IoU, centroid error, and OOD rejection. Compare
   with the frozen geometry baseline and trivial artifact baselines. Add saliency
   checks to verify predictions attend to embryo structures rather than overlays.
10. Add a website research view that shows: source/provenance, preprocessing,
    original stack slices and projections, segmentation, extracted geometry,
    predicted `tau` with uncertainty, held-out metrics, and explicit pilot/OOD
    warnings. Do not expose absolute local paths or raw data in production.

## Non-negotiable interpretation

- Augmentation regularizes a model; it does not create new independent embryos.
- The fixed DAPI/all-stain snapshots have no true `tau`, so they may support
  segmentation, preprocessing, self-supervised target-domain adaptation, and
  OOD testing, but not supervised pseudotime labels.
- The public rendered movies cannot be matched reliably to all 53 numeric
  trajectories. Do not invent mappings.
- The final end-to-end image model remains blocked on multiple independent raw
  untreated time-lapse stacks. Make the pipeline ready for those data and write
  a short acquisition checklist for the authors/lab: raw channels, frame times,
  z spacing, xy pixel size, treatment, embryo ID, pronuclear-formation frame,
  NEBD frame, and masks/annotations if available.

## Verification and deliverables

- Keep existing behavior and generated calibration artifacts intact.
- Add focused tests for TIFF metadata/projection, label handling, split leakage,
  timestamp masking, augmentation label preservation, and deterministic output.
- Run the existing pseudotime test suite plus new tests.
- Produce a short model card that separates validated results from pilot-only
  results and lists the exact data needed before publication-grade claims.
- Do not deploy, push, or fabricate performance numbers. Finish with a summary of
  files changed, commands run, verified results, and remaining data blockers.
