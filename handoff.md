# Handoff

Written 2026-08-14, at the end of the session that ported the paper figures onto this site.
It is the things a next agent would otherwise have to rediscover — the conventions, the traps,
and the reasoning behind decisions that look arbitrary from the outside.

`README.md` is the map of the repo. `AUTH.md` is access control. `AUDIT_vs_slideshow.md` is how
this site's methods compare to the paper figures. This file is the rest.

---

## 1. What this site is for, and what it is not

`/Users/rishib/Desktop/EmbyroPlayground/HighResSlideshowExports/` produces the **fixed** figures
for the paper. This repo produces the **interactive** ones: the statistics are recomputed in the
browser so a P-cutoff or a gene selection is a live control rather than something baked in.

**The reference is the authority on method.** When the two disagree about how an analysis should
be done, the reference wins and this site changes. That was the finding of the August 2026 audit
and it is why `embryo_stats.py` exists.

Where the reference cannot be reproduced, the page **says so on its face** rather than
approximating quietly. Two live examples: figure 7.5's GO dot plot is absent because no gene→term
annotation source exists on this machine, and 4.15's UMAP is drawn as PCA because `umap-learn` is
not installed. Both are stated in the page's own help text.

---

## 2. The house rules

These are not style preferences. Each one exists because getting it wrong produced a plausible,
wrong number that no figure would have revealed.

### Statistics — `embryo_stats.py`

Import it. Do not re-derive any of this.

* **Cytoplasm-only, by the per-molecule segment label `s`** — never by a containment test. A
  containment test is both slower and wrong at the boundary.
* **The body is found by volume, never by label number.** Label numbering is not stable across
  embryos. `classify_body()` returns two bodies when the second is >40% of the first, which is
  the empty gap between "sister blastomeres" (within ~2×) and "a pronucleus" (under a tenth).
* **Half-volumes come from exact mesh clipping with the origin ON the plane**, so the closing cap
  has zero scalar triple product and drops out. No shapely, no polygon ordering.
* **The bulk correction is a median of per-gene ratios, not `log(Σa/Σb)`.** One gene can be 30% of
  an embryo's cytoplasmic transcripts (up to 83%), so a ratio of totals is set by whichever gene
  is abundant and under-corrects everything else.
* **The null is volume-matched, not a fair coin.** Counting noise alone produces a fold; see §3.
* **`polar_label()` returns `None` rather than inventing an axis.** An earlier version took the
  most peripheral non-body segment, which on the nine zygotes that genuinely have no polar body
  picked a *pronucleus* and handed back a confident, meaningless axis.

### Naming — `embryo_naming.py` / `embryo-uids.js`

**The embryo label is looked up, never derived, and never read off a manifest or a cache.**
`data/embryo_ids.json` is generated from `CompleteEmbryoDataset.xlsx`; both this site and the
MERFISH atlas read it, so neither can drift.

This bit people twice in one day, so it is worth stating plainly: an artifact whose labels were
rewritten in place will have them **silently reverted by the next rebuild** unless its builder
calls `embryo_label()`. `build_pronuclei_pseudotime.py` did exactly that, because it read the
label from a geometry cache CSV that is only rewritten under `--extract`. Both its extract path
and its read path look the ID up now.

`p1` in a legacy filename is the **plate**, not the probeset, and every plate splits across two
probesets — so a filename token can never stand in for the real probeset. Use `probesets.json`.

### Rendering — the constants in `viewer-core.js`

| constant | value | what it is |
|---|---|---|
| `BODY_OPACITY` | 0.13 | the cell body / segmentation meshes |
| `DOT_SIZE` | 2.5 | a transcript dot, on the slider's 0.5 step grid |
| `DOT_OPACITY` | 0.9 | the **focus** transcript layer |

Prefer `V.bodyTraces(scene)` over building meshes yourself; it gets the scene's own per-segment
colours, the house opacity and the dark-render toggle for free. If a page must build its own —
several do, to recolour pronuclei or blastomeres — import the constants rather than picking a
number.

**Context layers stay dimmer on purpose**, and that contrast carries information: "not counted"
(pronuclei/polar body) on the plane pages at 0.5–0.7, the all-gene background on Clustering at
1.6/0.30, "every other gene" on Render Check at 1.3/0.18, the observed cloud behind the
simulation on Diffusion at 0.5. Do not sweep those up to the house value.

---

## 3. Results that will mislead you if you don't know them

* **A random plane already gives a median fold asymmetry of 1.263.** Counting noise is asymmetric.
  Nothing can be read against 1.0. On that footing neither the equatorial nor the sperm plane
  beats random (P = 0.34, 0.53); only the best meridional plane does, and it is *selected* to.
* **The animal–vegetal analysis is a null result** — 16 genes at P < 0.05 against 18.9 expected by
  chance, nothing surviving BH. The page leads with that.
* **A half-volcano needs an alignment null.** Any consistent labelling of halves produces some
  calls. Re-run with each embryo's sides flipped at random and compare.
* **Gene *share* makes a zygote's two halves exact complements**, so any embedding pushes them
  apart and a pairing statistic comes out above 1. That is geometry. Under concentration (per µm³,
  *not* renormalised to a unit sum, which puts the complementarity straight back) the same halves
  pair at 0.19–0.47.
* **Every pronucleus-vs-pronucleus comparison is circular** against the four-test consensus that
  labelled them. On the full consensus the maternal pronucleus is nearer the polar body in **41 of
  41** zygotes — that is the labelling, not a finding. Re-run with the relevant test removed; the
  signal survives but shrinks honestly (→ 30/38, P = 0.005). Hand calls are not votes and survive
  every leave-one-out.
* **Some numbers printed under the deck's renders are drawing devices.** Figure 1.7 duplicates
  transcripts on sparse panels (`amplify_to`) so the cloud reads at print size; 6.1's top2–top4
  are sub-sampled to the main panel's counts so the four read as one series. Checking a real count
  against either manufactures a failure.
* **Which side is which is not free.** The sperm lies *on* the plane it defines — it is one of the
  three points that draw it — so it cannot name a side, and the pronuclei descend from the sperm.
  Total cytoplasmic count is the only intrinsic rule with no free sign.

---

## 4. Traps in the data

* **Two coordinate conventions.** `data/segments/` scenes are one isotropic pixel space where
  µm = pixel × 0.15 on all three axes (`z_scale` 7.0); `data/zygote/` scenes use `z_scale` 6.667.
  Mixing them misplaces geometry by ~5% in z. To move a z coordinate between them, recover the
  **frame index** by dividing, then rescale — do not scale the µm value directly.
* **`20260425_zygote_p2_3` is an OOCYTE**, despite its name and its `category` in the sperm CSV.
  Uncorrected, the zygote counts read 60 and 34 instead of 59 and 33.
* **Ten zygotes have no polar body**, so every polar-axis analysis is genuinely undefined for them.
  The exclusion is correct, not a truncation.
* **Plane 0 of every GFP stack is a saturated calibration frame** (p99.9 ~17–25k against 150–300
  in a real plane). Max-projecting the whole stack sets the display window ~150× above the sperm
  and the embryo renders black.
* **`data/` is ~450 MB committed to the repo** and is the site's own copy of the scenes. Builds
  read these, not the atlas — which means the two can drift when a scene is regenerated under
  `localdata/`.

---

## 5. How to add a project

Six files, and missing any one of them fails quietly:

1. `<key>.html` / `.js` / `.css`
2. `build_<key>.py` → `data/<key>.json[.gz]`
3. `scripts/test_<key>.py`
4. a `help.js` entry — **without it the `?` modal opens empty**
5. registration in `lib/projects.mjs`
6. a card, in `index.html` (public) or `admin.card.html` (admin-only)

**The two registration tiers are not interchangeable.** Public → the `PROJECTS` array *and* a card
in a landing group. Admin-only → the `ADMIN_ONLY_PAGES` set *and* a card in `admin.card.html`;
middleware 404s the page for everyone else, and it is kept out of `PROJECTS` so the per-user
access matrix does not offer a toggle for a page only an admin can open. Promoting a project means
moving it in **both** places at once.

### Write the test against the rules, not the arithmetic

The per-artifact tests are the most valuable thing in the repo and they are not what a test
usually is. They do not re-add the build's sums. They assert the things that would still produce
plausible output if they silently broke:

* a plane really splits at 0.500000
* a null actually **moves** when its input is shuffled (a null that echoed the observation would
  pass every other check)
* a leave-one-out really left the test out
* a gene never detected lands on exactly 1.0 rather than vanishing from an average
* a checker has not quietly **stopped checking** — a panel with no checks must say why

Several real bugs were found this way and none of them were visible in the figures.

---

## 6. State as of this handoff

27 public projects in four landing groups, 6 admin-only. 30 build scripts, 37 pages, 28 artifact
tests, 4 Node test files (`npm test`, 58 tests, fail-closed auth logic).

The six projects ported from the paper figures — `clocktx`, `stages`, `animalveg`,
`contacthalves`, `halves` — are now public under **4 · Volcano plots**. `renders` (Render Check)
stays admin-only: it is a check on the deck rather than a result.

### Known and open

* **An intermittent `npm test` failure.** One test in 58 failed twice in ~50 runs, both times
  while a CPU-saturating Python build was running, and has not reproduced in ~40 runs since. Not
  identified. The two files with any time dependence are `tests/session.test.mjs` and
  `tests/middleware.test.mjs`, but their margins (3600s/60s/10s) look far too wide to flake, so
  don't assume that is the cause.
* **Twelve builders still take an embryo label from a manifest** rather than looking it up:
  `pronuclei_clocks`, `pronuclei_enrichment`, `sperm_pseudotime`, `stage_expr`, `contact`,
  `pseudosperm`, `sperm_locations`, `alignment` (plus `alphabeta`, `compare_planes`, `scheffler`,
  `size`, which label *series* rather than embryos and are probably fine). Any of the first eight
  can revert the naming migration on its next rebuild. Not fixed because they are not mine to
  assume — ask first.
* **`compare-planes` is superseded by `halves`**, which covers the same four plane definitions
  with the reference's statistics. Retiring it would be a decision, not a tidy-up.
* **`contact` (a slab) and `contacthalves` (equal-volume halves)** answer the same question with
  different instruments. Both are kept deliberately.
* **GFP coverage is the most valuable gap in the dataset** — 56 of 157 embryos have a path, 42
  have coordinates — and it is the vision model's input.

### Rebuild costs

Most builds are seconds. `build_halves.py` is ~5 minutes because it loads every zygote's
meridional family; it honours a `HALVES_CACHE=<path>` env var that caches the per-embryo stage so
iterating on the statistics does not re-pay for the geometry. `build_animalveg.py` and
`build_clocktx.py` are ~2 minutes each.

---

## 7. Working notes

* **Verify in the browser, don't assume.** Several bugs this session were invisible in the data
  and obvious on screen — a plot rendering at 240px because the drawer was still animating open, a
  bottom drawer sitting underneath the ranking panel and clipping a bar.
* **A Plotly panel drawn while a drawer is opening gets the wrong width.** Kick it with
  `Plotly.Plots.resize` on the next frame and again at ~160 ms and ~420 ms.
* **`git add -A` swept an uncommitted file of the user's into one of my commits.** Stage
  deliberately.
* **Check whether something is a regression before "fixing" it.** Twice this session a value
  looked wrong and turned out to date from the project's original commit — worth stating out loud,
  because the fix is the same but the explanation is not.
