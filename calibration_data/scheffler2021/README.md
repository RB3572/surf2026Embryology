# Scheffler et al. 2021 — pronuclear-migration calibration cohort (derived data)

Training/validation source for the **Pronuclear Pseudotime Calibration** project.

## Provenance

| | |
|---|---|
| Paper | Scheffler et al. (2021), *Two mechanisms drive pronuclear migration in mouse zygotes*, Nat Commun 12:841 |
| DOI | https://doi.org/10.1038/s41467-021-21020-x |
| Upstream file | Springer Nature `Source Data.xlsx` (public supplementary source-data workbook) |
| Sheets used | `Figure 1b` (male/female pronucleus→cell-centre distance), `Figure S1b`/`S1c` (relative volumes), `Figure S1m` (per-embryo migration duration) |
| Licence | CC BY 4.0 (Nature Communications open access) |

Only the **small derived table** is committed here. The 4 MB source workbook, the
supplementary movie and the PDF are deliberately NOT committed — regenerate the CSV from
the workbook with `extract_source_data.py` (also copied here verbatim for reproducibility).

## Committed files

| File | SHA-256 (first 16) | Description |
|---|---|---|
| `scheffler_2021_control_zygote_trajectories.csv` | `af781a219b3b2e66` | 2,057 frame rows × 53 untreated zygotes |
| `extract_source_data.py` | — | Workbook → CSV extractor (expects `Source Data.xlsx` beside it) |

## Cohort

- **53 untreated (control) zygotes**, 2,057 frames, 35–41 frames per embryo.
- Sampled every **0.25 h** from pronuclear formation.
- Migration duration (pronuclear formation → NEBD) **8.75–11.75 h**.
- `normalized_time_tau = time_h / migration_duration_h`, observed range **0.0 – 0.9756**.
- 1 missing `male_relative_volume` (and hence `volume_sum`/`volume_difference`); all
  distance columns are complete.

## Column semantics

| Column | Deployable? | Notes |
|---|---|---|
| `embryo_id` | group key | never a feature |
| `time_h`, `migration_duration_h` | label only | ground truth |
| `normalized_time_tau` | **target** | `time_h / migration_duration_h` |
| `male_to_center_um`, `female_to_center_um` | ✗ | require male/female identity, which fixed MERFISH zygotes do not reliably have |
| `nearer_to_center_um`, `farther_to_center_um` | ✓ | identity-free (sorted), the deployable core |
| `distance_sum_um`, `distance_difference_um` | ✓ | identity-free |
| `male_relative_volume`, `female_relative_volume`, `volume_sum`, `volume_difference` | ✗ **LEAKY** | normalized to each pronucleus's **own future endpoint volume**; encodes the answer |

## Why the volume columns are excluded from every deployable model

The published relative volumes are each divided by that pronucleus's volume at the *end of
its own trajectory*. A fixed snapshot has no future, so the quantity is unmeasurable at
application time, and within the training table it directly encodes elapsed time. Including
it drops MAE to ≈0.05 — an upper bound reported in the model comparison strictly as a
**non-deployable** reference, never selectable as the production model.

## Known limitations of this calibration source

1. Provides pronucleus→**cell-centre** distances only. It does **not** contain the minimum
   surface-to-surface pronuclear gap used by the legacy site score, so this dataset cannot
   validate that legacy metric.
2. 2D-rendered supplementary movie only; no raw 3D stacks are public, so segmentation
   processing cannot be matched between the live and fixed pipelines.
3. No cell-boundary/cortex distances, polar-body geometry, or absolute volumes.
4. Single laboratory and imaging protocol → no independent batch for external validation.

Raw microscopy is available from the corresponding author on request; obtaining it is the
main path to lifting limitations 1–3.
