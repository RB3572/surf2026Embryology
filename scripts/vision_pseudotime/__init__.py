"""
Fixed-image pronuclear pseudotime — PILOT pipeline and data plumbing.

This package is a *pilot and data pipeline*, deliberately NOT an end-to-end
production vision clock. It exists to (a) make the current material reproducible
and honestly catalogued, and (b) be ready for the raw untreated time-lapse
stacks that the final image model is blocked on.

Module map (one module per item of the implementation brief):

  config       local-data configuration + path redaction          (brief item 1)
  manifest     sources.csv / movies.csv validation, out-of-git     (brief item 1)
  tiff_audit   metadata-only TIFF audit + sidecar requirement      (brief item 2)
  projection   incremental crop + deterministic 2.5D projections   (brief items 2,3)
  movies       movie frame extraction with overlay masking         (brief item 4)
  augment      deterministic, geometry-preserving augmentation     (brief item 5)
  splits       embryo/source-grouped splits BEFORE augmentation    (brief item 6)
  hybrid       image -> segmentation -> geometry -> frozen clock    (brief item 7)
  encoder      exploratory (feasibility-only) direct image encoder (brief item 8)
  evaluate     grouped metrics, trivial baselines, saliency        (brief item 9)
  report       redacted provenance/metrics JSON for the website    (brief item 10)

Nothing here fabricates performance numbers, invents movie->trajectory
mappings, or copies large data into the repository.
"""

__version__ = "vpt-0.1.0-pilot"

# the deployable, already-validated clock this pilot is built around
GEOMETRY_BASELINE_MODEL_VERSION = "pnpt-3.0.0"
