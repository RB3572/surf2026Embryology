# Acquisition checklist — data needed to unblock a fixed-image pseudotime clock

For the corresponding author / imaging lab. The pipeline in
`scripts/vision_pseudotime/` is ready; it is blocked only on data of the
following kind. Please provide, for **multiple independent untreated embryos**
(not one embryo re-imaged, and not perturbation conditions):

## Per embryo / per time-lapse

| Field | Why it is needed | Example |
|---|---|---|
| **Raw channels** (not a rendered movie) | rendered supplementary movies are 2D, overlaid, and cannot be mapped to the numeric trajectories | 16-bit TIFF/OME-TIFF per channel |
| **Channel identity** | to know which channel is DNA vs cell-surface vs other | `["H2B-mCherry", "MyrGFP"]` |
| **Frame times** (real timestamps) | tau must come from real time, never from frame index | seconds or hours since acquisition start, per frame |
| **z spacing (µm)** | 3D geometry and physical distances | e.g. `1.0` |
| **xy pixel size (µm)** | 3D geometry and physical distances | e.g. `0.15` |
| **Treatment / condition** | untreated only for the normal-development clock; perturbations are OOD | `untreated` / `1 uM nocodazole` |
| **Embryo ID / batch** | grouped splits — every frame of an embryo stays in one split | stable unique string |
| **Pronuclear-formation frame** | defines tau = 0 | integer frame index |
| **NEBD frame** | defines tau = 1 | integer frame index |
| **Masks / annotations** (if available) | supervised segmentation, centroid error, Dice/IoU | label TIFF with class semantics |

## Minimums before publication-grade claims

- **≥ ~8–10 independent untreated embryos** with the fields above (more is
  better; the clock needs a held-out-embryo test, so a handful cannot validate).
- **Verified voxel spacing** (xy + z) for every stack — a sidecar or trustworthy
  OME metadata. Without it, geometry stays in uncalibrated voxel units.
- **Verified label-class semantics** for any segmentation (which integer ID is
  cell / pronucleus). IDs alone are not enough.
- A **held-out embryo / batch** never seen in training or preprocessing tuning,
  with **frozen preprocessing**.

## Sidecar format (drop next to each stack)

`calibration_data/vision_pseudotime/sidecars/<stackname>.sidecar.json`:

```json
{
  "embryo_id": "untreated_emb_03",
  "voxel_um": [0.15, 0.15, 1.0],
  "voxel_um_verified": true,
  "channels": ["H2B-mCherry", "MyrGFP"],
  "label_classes": {"1": "cytoplasm", "2": "pronucleus"},
  "label_classes_verified": true,
  "acquisition_time": "2026-01-01T09:00:00",
  "treatment": "untreated"
}
```

Set `voxel_um_verified` / `label_classes_verified` to `true` only when confirmed;
the pipeline treats unverified values as provisional and withholds physical-unit
claims.
