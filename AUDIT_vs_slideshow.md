# Audit: site projects vs `HighResSlideshowExports`

`HighResSlideshowExports` is the truth for how the analysis should be conducted. This file records
where the site's pipelines differ from it, and which of its analyses have no project yet.

Audited 2026-08-13, against `_work/cells/f*.py` (38 figure modules) and `_work/scene_export.py`
(the shared "house definitions"). **No code was changed.**

---

## Part 1 — Pipeline differences in projects that exist

### 1.1 Differences that change numbers

| # | What | Reference | Site | Affects |
|---|---|---|---|---|
| **A** | **Equatorial plane placement** | `scene_export.equal_volume_plane`: ⟂ the PB←COM axis, then **shifted along it (brentq on the mesh) until cytoplasm volume splits exactly 50/50**. Explicitly: "a plane merely through the centre of mass does NOT bisect volume … lands between 0.490 and 0.509 across the 50 zygotes" | `build_equatorial.py` **does** shift — but by the **median of voxel projections**, not a mesh search. `build_compare_planes.py` does **NOT** shift at all: `sideA(pos1, com, n_eq)` puts it **through the COM** | `equatorial-planes` (estimator differs), `compare-planes` (**wrong plane**) |
| **B** | **The null** | Volume-matched: `rng.binomial(cnt, vp/(vp+vm))` — the null proportion is the volume split (4.17, 4.19, 4.21, 9.2), 200 draws | `build_zygote.py` / `build_sperm_division.py`: `rng.binomial(n, 0.5)` — a **fair coin**, ignoring that the two halves have different volumes | `zygote-planes`, `sperm-division`, `equatorial-planes`, `planes-all` |
| **C** | **Cross-embryo statistics** | One-sample **t-test** of the per-embryo log₂ fold against 0, with `MIN_TOTAL = 20` (summed) and `MIN_EMBRYOS = 5`, then **BH** adjustment (4.19, 4.21) | Only `build_pseudosperm.py` does this. Every other project reports **per-embryo permutation p-values** and never combines across embryos, and nothing else applies BH | `zygote-planes`, `sperm-division`, `equatorial-planes`, `compare-planes`, `contact` |
| **D** | **Bulk correction** | Each gene's log ratio has that embryo's **median per-gene log ratio** subtracted (4.19, 4.21, 8.3). Rationale in `f4_21`: one gene can carry 30% of an embryo's cytoplasm (up to 83%), and abundant genes are less asymmetric, so a total-based correction under-corrects everything else | Only `build_pseudosperm.py`. No other project corrects for bulk at all | all half/plane projects |
| **E** | **Pseudocount** | `+0.5` on every count before the log (`EPS = 0.5` in 8.3; `(a+0.5)/vp` in 4.21) | Absent everywhere except `build_pseudosperm.py` | any project taking a log of a count |
| **F** | **Side orientation** | A plane's stored `a`/`b` are an arbitrary geometric convention with **no anatomical meaning** — 4.15 establishes half A is uncorrelated with which side holds more. Signed quantities are oriented before averaging: **by total cytoplasmic count** (4.21) or **by the UMAP pairing eigenvector** (4.19) | No site project orients before aggregating except `pseudosperm` (count) and `alignment` (its anchor). `compare-planes` and the `*_cross` aggregates use \|absolute\| values, which sidesteps but also discards direction | `compare-planes`, cross-embryo drawers |
| **G** | **Contact normalisation** | 7.1–7.3 use the house density: count ÷ **cytoplasm-only volume** on that side | `build_contact.py` deliberately uses a **transcriptome-fraction null** (`fold = (k/n) ÷ f0`) because it ships meshes not voxel masks. Documented in the build, but it is a different quantity from the reference's | `contact` |
| **H** | **Body/segment identification** | `scene_export.classify` — by **geometry and volume, never by label number** ("labels are not consistent across embryos") | `build_zygote.py`, `build_planes_all.py`, `build_equatorial.py` and others hard-code **segment 1 = cytoplasm**. `build_contact.py`, `build_alignment.py`, `build_alphabeta.py` and the new `build_pseudosperm.py` do read it off volume | `zygote-planes`, `planes-all`, `equatorial-planes`, `sperm-division`, `compare-planes` |

### 1.2 Differences that are estimator-level, not method-level

| # | What | Note |
|---|---|---|
| I | **Volume estimator** | Reference integrates the **mesh** (and notes the body mesh already excludes nuclei, because the cavity shells integrate out). Site mostly counts **downsampled voxels** (`DS_XY = 6`, `DS_Z = 2`). Same quantity, ~0.2–1% apart |
| J | **Exhaustive grid** | Reference "20,004-normal orientation grid"; site `M_PLANES = 20000` Fibonacci hemisphere. Effectively the same |
| K | **Meridional fan** | Reference and site both use 18 planes at 10° for the fan; `pseudosperm` uses 180 at 1° (`N_AZIMUTH = 180`, matching 4.21) |

### 1.3 Where the site and the reference already agree

- `pseudosperm` — a verified port of 4.21 (reproduced its ranking exactly).
- `scheffler` + `pronuclei` + `extpt` — the isotonic clock and its validation match 4.5 / 10.5.
- `size` — matches 10.1–10.4 (same definitions, re-exported together).
- `alignment` — the same construction as 4.6 (average 2-cell outline, sperm projected to the
  cortex), generalised: the reference draws one fixed panel, the site makes the anchor a control.
- `sperm-sphere`, `sperm-map`, `clustering`, `stage-expression`, `pronuclei-assignments` — no
  reference figure covers them; they are site-only and nothing conflicts.

---

## Part 2 — Reference analyses with no project

38 figures. ~13 are 3-D renders (galleries), the rest are charts. Grouped into **six** proposed
projects rather than one per figure.

### P1 · "Transcriptome vs the clock" — 4.8, 4.11, 5.4
Per-gene **share** (a gene's count ÷ that embryo's total — composition, not abundance),
**centred within probeset** because the four panels differ by an order of magnitude in total
counts, correlated against τ. 4.8 is the Spearman volcano, 4.11 the four exemplar trajectories
(Zp1, Zp2, Pin1, Cdc42), 5.4 the total-count decline.
*Shape:* volcano + click-through per-gene trajectory, τ on x. Reuses the existing τ.
*Note:* `MIN_COUNT = 7`, `MIN_ZYGOTES_FOR_TEST = 2`; the probeset centring is essential and is a
constraint the site has hit before.

### P2 · "Halves, paired" — 4.15, 4.19, and the honest version of 4.14/4.17/4.18
The half-UMAP pairing (each zygote contributes two half-points; the **side axis is the leading
eigenvector of Σ d·dᵀ**, which is what makes the sign meaningful), then the half-enrichment
volcanoes under two plane definitions, then fold-asymmetry per gene against the count-matched
null and summarised by plane definition.
*Shape:* this is the natural home for differences **B, C, D, E, F** above — it is the reference's
own answer to "how do you aggregate a signed half-difference across embryos". Could **absorb or
replace `compare-planes`**, which currently answers the same question with a weaker method and
the wrong equatorial plane.

### P3 · "Animal–vegetal" — 4.3, 4.4
Density disks (48×48 meridional map, log₂ vs the all-gene map) and the volcano on the
**equal-cytoplasmic-volume** split. The site has the plane (`equatorial-planes`) but neither the
2-D map nor the volcano.
*Shape:* add two drawer tabs to `equatorial-planes` rather than a new project — provided A is
fixed first, since the volcano is defined on the equal-volume split.

### P4 · "Across the stages" — 8.3, 8.6, 9.1, 9.2
Between-half fold change for every gene at zygote / early / late 2-cell (8.3), the retained-vs-lost
percentile plot (8.6), and the two heatmaps (9.1 drops/peaks/rises, 9.2 count-matched percentiles).
*Shape:* one project — a gene × stage matrix with the trajectory classes as the organising idea.
Complements `stage-expression` (which shows **levels**, not **between-half fold**) and would sit
beside it.
*Note:* 8.3 explicitly is **not** built from `fold_per_stage_all_genes.csv` (r = 0.56 against the
right computation) — worth heeding.

### P5 · "Contact, completed" — 7.3, 7.5 (+ 7.1's profile)
The site's `contact` project covers 7.1/7.2 in spirit. Missing: the pooled **spatial density maps**
of the two leaning gene sets (34×34, log₂ vs all-gene) and the **GO dot plot**.
*Shape:* two more tabs on `contact`. The GO panel needs an ontology source the repo does not have
— that is the blocker, not the plotting.

### P6 · "Render gallery" — 1.7, 4.2, 4.9, 4.10, 4.12, 5.1, 5.2, 5.5, 5.6, 6.1, 6.2, 7.4, 8.2, 8.7
Thirteen 3-D renders plus 4.9's 18-panel edge-on grid. Every one is a **hand-picked embryo + gene**
chosen to illustrate a definition, and the site already renders all of these scenes.
*Shape:* one "Figure Gallery" project — pick a figure, see the exact embryo/gene the deck uses,
with the number it prints (fold, LFC, τ, counts) recomputed live so the page and the deck can be
checked against each other. Cheap, and it makes the deck reproducible from the site.
*Note:* several carry deliberate departures that must be preserved as annotations, e.g. 4.2's
polar-body **size augmentation**, 8.7's Ltbp1 panel (2–3 transcripts — flagged weak in the source),
7.4's Trib3 count disagreement with the deck (257/244 printed vs 262/239 computed).

### Not proposed
- **4.7** (observed vs null, sperm near the contact line) — its search-corrected P = 0.016 was
  never written to a table and is transcribed from a README in the reference itself. Would need
  recomputing from scratch; flag before building.
- **4.13, 4.16, 6.3** — largely covered by `sperm-map` / `sperm-pseudotime` / `clustering` +
  `alphabeta` respectively. 4.16 (sperm→maternal vs paternal, paired) is the closest to a real
  gap and could be one tab on `sperm-pseudotime`.

---

## Suggested order

1. **A** — `compare-planes` is using a plane that is not the equatorial plane. Smallest fix,
   clearest wrongness.
2. **B/C/D/E** — adopt the reference's null, t-test, BH, bulk correction and pseudocount as a
   shared helper, then apply project by project. `build_pseudosperm.py` already implements all of
   them and can be the source.
3. **P1** and **P4** — the two largest genuinely-missing analyses.
4. **P2** — the biggest single piece of work, and the one that would let `compare-planes` be
   retired rather than patched.
