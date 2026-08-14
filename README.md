# SURF 2026 Embryology — the analysis site

Deployed to **incrementum.rishib.com** (Vercel, redeploys on push). 33 interactive analyses of the
MERFISH mouse zygote and early-2-cell dataset: 3-D scenes, division planes, pseudotime clocks,
gene enrichment and volcano plots, each recomputing its statistics in the browser so thresholds and
gene selections are live controls rather than baked-in choices.

Companion docs: **[handoff.md](handoff.md)** (the conventions and traps — read this before
changing an analysis), **[AUTH.md](AUTH.md)** (access control),
**[AUDIT_vs_slideshow.md](AUDIT_vs_slideshow.md)** (how these methods compare to the paper figures).

There is **no build step** for the site itself. `vercel.json` pins the project to *No Framework*
and serves the repo root. `.json.gz` artifacts are gunzipped in-browser via `DecompressionStream`.

---

## Anatomy of a project

Every project is the same six files, and missing any one of them fails quietly:

| File | Role |
|---|---|
| `<key>.html` / `.js` / `.css` | the page |
| `build_<key>.py` | produces the artifact; run from the repo root |
| `data/<key>.json[.gz]` | the artifact the page loads |
| `scripts/test_<key>.py` | asserts the artifact's **rules**, not its arithmetic |
| `help.js` entry | the `?` modal — without it the modal opens empty |
| `lib/projects.mjs` + a card | registration, see below |

### Registration has two tiers, and they are not interchangeable

* **Public** — add to the `PROJECTS` array in `lib/projects.mjs` *and* a card in one of
  `index.html`'s landing groups. Appears for every lab member; can be restricted per user from the
  admin console.
* **Admin-only** — add to the `ADMIN_ONLY_PAGES` set *and* a card in `admin.card.html`. Middleware
  404s the page for everyone else. Deliberately kept out of `PROJECTS` so the per-user access
  matrix does not offer a toggle for a page only an admin can open.

Promoting a project means moving it in **both** places at once.

---

## Shared code

| File | What |
|---|---|
| **`embryo_stats.py`** | the house statistics. Cytoplasm-only counting by segment label, body-by-volume, exact half-space clipping, `equal_volume_plane()`, median-of-ratios bulk correction, the volume-matched null, BH. **Import it rather than re-deriving any of it.** |
| **`embryo_naming.py`** / `embryo-uids.js` | the canonical embryo ID, **looked up** from `data/embryo_ids.json`, never derived. Both this site and the MERFISH atlas read it, so neither can drift. |
| **`viewer-core.js`** | shared 3-D and layout: `bodyTraces`, `sceneLayout`, `plotConfig`, `loadGz`, `buildTabs`, `wireWindow`, the drawer chrome, the dark-render toggle, and the rendering constants `BODY_OPACITY` / `DOT_SIZE` / `DOT_OPACITY`. |
| `lib/`, `api/`, `middleware.js` | Lab Logger SSO, per-user project access (Neon), analytics. See `AUTH.md`. |

---

## The projects

**1 · Pseudotime** — `pronuclei-assignments`, `scheffler`, `pronuclei`, `extpt`

**2 · Sperm, division & body axes** — `sperm-map`, `size`, `alignment`, `alphabeta`, `sperm-pca`,
`sperm-pseudotime`, `zygote-planes`, `pseudosperm`, `sperm-division`, `planes-all`,
`equatorial-planes`, `compare-planes`, `sperm-sphere`

**3 · Gene expression & spatial pattern** — `contact`, `clustering`, `stage-expression`,
`segments`, `diffusion`

**4 · Volcano plots** — `clocktx`, `stages`, `animalveg`, `contacthalves`, `halves`

**Admin only** — `renders` (Render Check), `pronuclear-pseudotime`, `pseudotime-calibration`,
`vision-pseudotime`, `pn3d-transcripts`, `axes`, plus the usage-analytics console.

---

## Data

`data/` is ~450 MB committed to the repo — the site's own copy of the scenes, so it builds and
serves without the external volume mounted.

| Directory | Scenes | What |
|---|---|---|
| `segments/` | 157 | **the population.** Every embryo, all four stages. Per-transcript compartment label `s` and each label's exact voxel volume. |
| `zygote/` | 50 | zygotes with a polar-body axis: 18 meridional planes, `best_planes` |
| `equatorial/` | 50 | the same 50, equatorial split |
| `sperm_division/` | 30 | pre-built sperm-plane scenes |
| `planes_all/` `pronuclei/` `axes/` `extpt/` `gene_diffusion/` `scenes/` | 50, 52, 54, 93, 51, 46 | per-analysis products |

⚠️ **Only `segments/` is the full population.** Treating any other directory as "all the embryos"
silently drops real data.

⚠️ Two coordinate conventions: `segments/` scenes are one isotropic pixel space (µm = pixel × 0.15
on all three axes, `z_scale` 7.0); `zygote/` scenes use `z_scale` 6.667. Mixing them misplaces
geometry by ~5% in z.

---

## Tests

```bash
npm test
```

58 Node tests over the auth and access logic — session sign/verify, the `?next=` open-redirect
guard, both branches of the lab-membership check, and the middleware enforcement matrix. These
**fail closed** by design.

```bash
python3 scripts/test_<project>.py
```

28 per-artifact checkers, one per analysis. They do not re-add the build's sums; they assert the
things that would still produce plausible output if they silently broke — that a plane really
splits at 0.5, that a null actually moves when its input is shuffled, that a leave-one-out really
left the test out. See `handoff.md` §5.

---

## Regenerating

Every `build_*.py` runs from the repo root and reads `data/segments/` etc. Most take seconds.
`build_halves.py` is ~5 minutes and honours `HALVES_CACHE=<path>` to cache its per-embryo stage;
`build_animalveg.py` and `build_clocktx.py` are ~2 minutes each.

`build_renders.py` additionally reads the reference figure export at
`~/Desktop/EmbyroPlayground/HighResSlideshowExports/` and will exit if it is not present.
