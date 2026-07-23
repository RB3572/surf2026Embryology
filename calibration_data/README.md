# Pronuclear pseudotime calibration — reproducing every artifact

This directory holds the **derived, checked-in inputs**. Everything the website shows can be
regenerated from them with no raw microscopy, no sibling checkout, and no absolute paths.

## One command sequence (run from the repository root)

```bash
python3 scripts/train_pronuclear_pseudotime.py      # nested CV, conformal, frozen model
python3 build_pronuclei_pseudotime.py               # apply the frozen clock to fixed zygotes
python3 scripts/analyze_pseudotime_noise_ceiling.py # reliability, attenuation, MERVL illustration
python3 scripts/test_pronuclear_pseudotime.py       # 92 checks — must print "0 failed"
```

Total runtime ≈ 3–4 min (add `--quick` to the first and third for a ~40 s smoke run; `--quick`
reduces bootstrap/simulation replicates only and must **not** be used for reported numbers).

Outputs:

| Path | Written by |
|---|---|
| `data/pseudotime_calibration/calibration.json` | train |
| `data/pseudotime_calibration/oof_frames.json.gz` | train (nested **outer-test** predictions) |
| `data/pseudotime_calibration/trajectories.json.gz` | train |
| `data/pseudotime_calibration/model.json` | train (frozen model card) |
| `data/pseudotime_calibration/noise_ceiling.json` | analyze |
| `data/pronuclei_pseudotime.json` | build_pronuclei_pseudotime |

## Checked-in inputs

### `scheffler2021/`
The training/validation cohort — 53 untreated live-imaged zygotes, 2,057 frames. See
`scheffler2021/README.md` for full provenance, licence and column semantics. Only the small derived
CSV is committed; the 4 MB Springer Nature workbook, the supplementary movie and the PDF are
deliberately **not** in the repository. Regenerate the CSV with `scheffler2021/extract_source_data.py`
if you place `Source Data.xlsx` beside it.

### `fixed_cohort_geometry.csv`
Per-zygote **derived geometry** for the 51 fixed MERFISH zygotes: the two pronuclear-centroid →
cell-centre distances (sorted nearer/farther), their sum and difference, cell and pronuclear
volumes, and the legacy minimum surface gap. **Geometry only — no transcripts, no gene identities.**

This cache exists so the *apply* step is self-contained. Extraction needs the raw label TIFFs
(correctly absent from this repo); caching the derived features means re-applying a new model
version takes 0.2 s instead of ~7 min and needs no raw data at all.

To re-derive it when the raw data *is* available:

```bash
python3 build_pronuclei_pseudotime.py --extract
```

which needs `build_pronuclei.py`, `data/pronuclei_manifest.json` and the
`TranscriptomicsData/JustTifAndCSVData/Zygote` tree.

## What is deliberately NOT here

- Raw microscopy (TIFF stacks) — not public, and not needed for any reported number.
- The Scheffler workbook / movie / PDF — regenerable, and large.
- Any transcript quantity used by the clock — the clock never sees one. Transcript counts are read
  only by `analyze_pseudotime_noise_ceiling.py`, and only *after* the clock is frozen on disk.
