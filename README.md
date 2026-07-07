# SURF 2026 Embryology — Sperm · Embryo 3D Viewer

Deployed (Vercel) static site for visualizing each sperm-positive embryo in 3-D:
the embryo body (segmentation meshes), a selected gene's transcript point cloud,
and the sperm location — plus PCA / alignment vectors and a cross-embryo gene
ranking. One nav-bar tab per embryo (45 total).

## Structure (pure static site, served from the repo root)

| Path | What |
|---|---|
| `index.html`, `app.js`, `style.css` | the app |
| `plotly.min.js` | 3-D rendering engine |
| `data/manifest.json` | the 45-embryo index for the nav bar |
| `data/analysis_index.json.gz` | per-embryo / per-gene unit vectors (violins + ranking) |
| `data/scenes/<id>.json.gz` | per-embryo geometry (meshes + transcript clouds + sperm + analysis) |

There is **no build step**. `vercel.json` pins the project to *No Framework* and
serves the repo root directly. The `.json.gz` scenes are gunzipped in-browser via
`DecompressionStream` (the loader also handles the case where the host transparently
decodes them).

## Regenerating the data

The `data/` files are produced by [`build_viewer_data.py`](build_viewer_data.py)
from the source dataset in the companion **SpermLabeling** project: it reads the
MERFISH atlas geometry and this project's sperm dataset, then writes the slim
per-embryo scenes, the manifest, and the cross-embryo analysis index. Run it there
(it uses absolute paths to that project), then sync the refreshed `public/` output
into this repo and push — Vercel redeploys on push.
