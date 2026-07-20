/* Sperm · Embryo 3D Viewer
 *
 * Minimal front-end: a top nav-bar of the 46 sperm-positive embryos; selecting
 * one renders its 3-D scene — the embryo body (segmentation meshes), the point
 * cloud for the selected gene, and the sperm location. Nothing else.
 *
 * The 3-D rendering recipe (mesh3d body + scatter3d cloud/sperm, z_scale,
 * pinned-extent scene) mirrors the MERFISH atlas viewer so the visual style
 * matches. Per-embryo scene data is fetched from data/scenes/<id>.json.gz and
 * gunzipped in-browser via DecompressionStream.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  let vcExtras = null;   // dot-size + atlas-link row (VCore.addWindowExtras)
  const tabsEl     = $("#tabs");
  const controlsEl = $("#controls");
  const controlsHeader = $("#controls-header");
  const geneSelect = $("#gene-select");
  const vecTogglesEl = $("#vec-toggles");
  const readoutEl  = $("#readout");
  const plotHost   = $("#plot-host");
  const placeholder= $("#placeholder");
  const loadingEl  = $("#loading");
  const loadingTxt = $("#loading-text");
  const countEl    = $("#embryo-count");
  const pcaPlot    = $("#pca-plot");
  const pcaGeneEl  = $("#pca-gene");
  const drawer     = $("#drawer");
  const drawerHandle = $("#drawer-handle");
  const drawerGeneEl = $("#drawer-gene");
  const violinPca  = $("#violin-pca");
  const violinG2e  = $("#violin-g2e");
  const violinPcaSub = $("#violin-pca-sub");
  const violinG2eSub = $("#violin-g2e-sub");
  const violinWrap = $("#violin-wrap");
  const violinResize = $("#violin-resize");
  const drawerBody = $("#drawer-body");
  const drawerResize = $("#drawer-resize");
  const rdrawer    = $("#rdrawer");
  const rdrawerHandle = $("#rdrawer-handle");
  const rdrawerResize = $("#rdrawer-resize");
  const rtabsEl    = $("#rtabs");
  const rankListEl = $("#rank-list");
  const rankFilterEl = $("#rank-filter");
  const predictEnable = $("#predict-enable");
  const predictConfig = $("#predict-config");
  const predictTopN   = $("#predict-topn");
  const predictLoo    = $("#predict-loo");
  const predictReadout = $("#predict-readout");
  const nullShowEl = $("#null-show");
  const nullRegenEl = $("#null-regen");

  const VIOLIN_W_KEY = "sperm_viewer_violin_w";
  const DRAWER_H_KEY = "sperm_viewer_drawer_h";
  const RDRAWER_W_KEY = "sperm_viewer_rdrawer_w";
  const PREDICT_KEY = "sperm_viewer_predict";
  const NULL_DIRS_KEY = "sperm_viewer_null_dirs";
  const NULL_SHOW_KEY = "sperm_viewer_null_show";
  const RANK_MIN_N = 3;   // min embryos for a gene to be ranked (meaningful σ)

  function loadNullDirs() {
    try { return JSON.parse(localStorage.getItem(NULL_DIRS_KEY) || "{}") || {}; }
    catch (_) { return {}; }
  }

  function loadPredictCfg() {
    try {
      const v = JSON.parse(localStorage.getItem(PREDICT_KEY) || "{}");
      return { enabled: !!v.enabled, topN: v.topN || "10", loo: v.loo !== false };
    } catch (_) { return { enabled: false, topN: "10", loo: true }; }
  }
  function savePredictCfg() {
    try { localStorage.setItem(PREDICT_KEY, JSON.stringify(state.predict)); } catch (_) {}
  }

  // Stage accent colors (match the nav bar), used to color violin points.
  const STAGE_COLOR = {
    "Zygote": "#6366f1", "2-cell (early)": "#0ea5e9",
    "2-cell (late)": "#14b8a6", "Oocyte": "#f59e0b",
  };
  function hexA(hex, a) {
    const h = hex.replace("#", "");
    return `rgba(${parseInt(h.slice(0, 2), 16)},${parseInt(h.slice(2, 4), 16)},` +
           `${parseInt(h.slice(4, 6), 16)},${a})`;
  }

  const STAGE_VAR = {                       // stage -> CSS accent variable
    "Oocyte": "--oocyte", "Zygote": "--zygote",
    "2-cell (early)": "--e2c", "2-cell (late)": "--l2c",
  };
  // Gene point-cloud styling (single selected gene).
  const GENE_COLOR = "#ff7f00";
  const GENE_SIZE  = 0.5;

  // Sperm-condition colors: green = real, red = randomized null.
  const REAL_COLOR = "#16a34a";   // real sperm (diamond) + predicted (circle)
  const NULL_COLOR = "#dc2626";   // null sperm (diamond)

  // The three pre-computed analysis vectors (see build_viewer_data.py).
  const VECTORS = [
    { key: "pca",       label: "PC1 axis",    color: "#7c3aed",
      desc: "First principal-component axis of this gene's point cloud" },
    { key: "gene2emb",  label: "Gene → COM",  color: "#2563eb",
      desc: "Gene cloud centre of mass → embryo centre of mass" },
    { key: "sperm2emb", label: "Sperm → COM", color: "#16a34a",
      desc: "Sperm → embryo centre of mass" },
  ];
  const VEC_STORE_KEY = "sperm_viewer_vectors";

  function loadVectorToggles() {
    try {
      const v = JSON.parse(localStorage.getItem(VEC_STORE_KEY) || "{}");
      return { pca: !!v.pca, gene2emb: !!v.gene2emb, sperm2emb: !!v.sperm2emb };
    } catch (_) { return { pca: false, gene2emb: false, sperm2emb: false }; }
  }

  const state = {
    manifest: [],
    byId: new Map(),
    currentId: null,
    dotSize: 1.5,         // transcript dot marker size (floating-window control)
    scene: null,          // decoded scene for the current embryo
    sceneCache: new Map(),
    // Vector-toggle visibility is GLOBAL: it persists as the user switches
    // embryos via the nav bar (and across reloads, via localStorage).
    vectors: loadVectorToggles(),
    // The user's chosen gene sticks across embryo switches (when that embryo's
    // panel contains it); null until the user picks one explicitly.
    userGene: null,
    index: null,          // cross-embryo analysis vectors (for the violins)
    drawerOpen: false,
    rankAxis: "pca",      // "pca" | "g2e"  (active tab's axis)
    rankCond: "real",     // "real" | "null"  (active tab's condition)
    rankFilter: 0,        // "> N transcripts in every embryo" (0 = all)
    rankings: null,       // { real:{pca,g2e}, null:{pca,g2e} }
    geneToEmbryos: null,  // Map<gene, embryoId[]>  for jump-to-gene
    predict: loadPredictCfg(),   // { enabled, topN, loo } — persisted
    prediction: null,     // computed result for the current embryo
    nullDirs: loadNullDirs(),    // embryoId -> random µm unit direction
    nullPoint: null,      // null-sperm shell hit for the current embryo (plot)
    showNull: localStorage.getItem(NULL_SHOW_KEY) !== "0",   // marker visibility
  };

  // ---------- fetch + gunzip a scene (mirrors the atlas loader) ----------
  async function loadScene(url) {
    const resp = await fetch(url);
    if (!resp.ok) throw new Error("HTTP " + resp.status + " for " + url);
    let buf = new Uint8Array(await resp.arrayBuffer());
    if (buf.length >= 2 && buf[0] === 0x1f && buf[1] === 0x8b) {
      if (typeof DecompressionStream === "undefined")
        throw new Error("This browser lacks DecompressionStream; please update it.");
      const ds = new DecompressionStream("gzip");
      const stream = new Response(buf).body.pipeThrough(ds);
      buf = new Uint8Array(await new Response(stream).arrayBuffer());
    }
    return JSON.parse(new TextDecoder().decode(buf));
  }

  // ---------- nav bar ----------
  function buildTabs(manifest, stageOrder) {
    tabsEl.innerHTML = "";
    const groups = new Map();
    for (const e of manifest) {
      if (!groups.has(e.stage)) groups.set(e.stage, []);
      groups.get(e.stage).push(e);
    }
    const order = stageOrder.filter((s) => groups.has(s));
    for (const stage of order) {
      const g = document.createElement("div");
      g.className = "stage-group";
      g.style.setProperty("--stage", `var(${STAGE_VAR[stage] || "--accent"})`);
      const tag = document.createElement("span");
      tag.className = "stage-tag";
      tag.textContent = stage.replace("2-cell ", "2c ").replace(/[()]/g, "");
      g.appendChild(tag);
      for (const e of groups.get(stage)) {
        const btn = document.createElement("button");
        btn.className = "tab";
        btn.dataset.id = e.id;
        btn.style.setProperty("--stage", `var(${STAGE_VAR[stage] || "--accent"})`);
        btn.title = `${e.title}  ·  ${e.n_transcripts.toLocaleString()} transcripts` +
                    `  ·  record ${e.merfish_index}`;
        btn.innerHTML =
          `<span class="tab-label">${e.label}</span>` +
          `<span class="tab-date">${e.date_short || ""}</span>`;
        btn.addEventListener("click", () => selectEmbryo(e.id));
        g.appendChild(btn);
      }
      tabsEl.appendChild(g);
    }
  }

  function markActiveTab(id) {
    tabsEl.querySelectorAll(".tab").forEach((b) =>
      b.classList.toggle("active", b.dataset.id === id));
    const active = tabsEl.querySelector(".tab.active");
    if (active) active.scrollIntoView({ inline: "center", block: "nearest", behavior: "smooth" });
  }

  // ---------- select + render an embryo ----------
  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    markActiveTab(id);
    const meta = state.byId.get(id);
    showLoading(`Loading ${meta.label}…`);
    try {
      let scene = state.sceneCache.get(id);
      if (!scene) {
        scene = await loadScene(meta.scene);
        state.sceneCache.set(id, scene);
      }
      if (state.currentId !== id) return;   // user switched away mid-load
      state.scene = scene;
      if (vcExtras) vcExtras.setAtlas(id);
      populateGenes(scene);
      controlsEl.hidden = false;
      placeholder.hidden = true;
      drawer.hidden = false;
      computePrediction();          // depends on the new embryo; render() draws it
      computeNullPoint();           // null-sperm shell hit for this embryo
      render();
      updateAnalysisPlots();
      updatePredictReadout();
    } catch (err) {
      showError(err.message || String(err));
    } finally {
      hideLoading();
    }
  }

  function populateGenes(scene) {
    const totals = scene.gene_totals || {};
    const opts = scene.genes.map((g) =>
      `<option value="${g}">${g}  (${(totals[g] || 0).toLocaleString()})</option>`);
    // The selected gene ALWAYS persists across embryos. If this embryo's panel
    // doesn't contain it, keep it as a marked, selectable option so it is never
    // silently reset — the 3-D cloud/PCA just show "not in this embryo".
    const sel = state.userGene || scene.genes[0];
    if (sel && !scene.genes.includes(sel))
      opts.unshift(`<option value="${sel}">${sel}  (not in this embryo)</option>`);
    geneSelect.innerHTML = opts.join("");
    geneSelect.value = sel;
  }

  geneSelect.addEventListener("change", () => {
    state.userGene = geneSelect.value;
    render();
    updateAnalysisPlots();
  });

  // ---------- analysis-vector toggle toolbar (persistent, global) ----------
  function buildVectorToolbar() {
    vecTogglesEl.innerHTML = "";
    VECTORS.forEach((v) => {
      const btn = document.createElement("button");
      btn.type = "button";
      btn.className = "vec-btn";
      btn.dataset.vec = v.key;
      btn.style.setProperty("--vc", v.color);
      btn.title = v.desc;
      btn.setAttribute("aria-pressed", String(!!state.vectors[v.key]));
      btn.classList.toggle("on", !!state.vectors[v.key]);
      btn.innerHTML =
        `<span class="vec-arrow" aria-hidden="true"></span>` +
        `<span class="vec-name">${v.label}</span>`;
      btn.addEventListener("click", () => {
        state.vectors[v.key] = !state.vectors[v.key];
        btn.classList.toggle("on", state.vectors[v.key]);
        btn.setAttribute("aria-pressed", String(state.vectors[v.key]));
        try { localStorage.setItem(VEC_STORE_KEY, JSON.stringify(state.vectors)); } catch (_) {}
        render();
      });
      vecTogglesEl.appendChild(btn);
    });
  }

  // A 3-D arrow = a line shaft (scatter3d) + a cone arrowhead. The direction is
  // a physical (µm) unit vector; we map it into the anisotropic plot space
  // (xy: px = µm/0.15; z: frame*z_scale) so it points correctly on screen, then
  // draw a fixed visual length L from `anchor` (plot space).
  function arrowTraces(anchor, physUnit, L, zs, color, name) {
    if (!physUnit || !anchor) return [];
    let dx = physUnit[0] / 0.15, dy = physUnit[1] / 0.15, dz = physUnit[2] * zs;
    const n = Math.hypot(dx, dy, dz);
    if (!n) return [];
    const ux = dx / n, uy = dy / n, uz = dz / n;         // unit in plot space
    const head = [anchor[0] + ux * L, anchor[1] + uy * L, anchor[2] + uz * L];
    const shaft = {
      type: "scatter3d", mode: "lines",
      x: [anchor[0], head[0]], y: [anchor[1], head[1]], z: [anchor[2], head[2]],
      line: { color, width: 7 },
      name, legendgroup: name, showlegend: true, legendrank: 40000,
      hovertemplate: `${name}<extra></extra>`,
    };
    const cone = {
      type: "cone", x: [head[0]], y: [head[1]], z: [head[2]],
      u: [ux], v: [uy], w: [uz], anchor: "tip",
      sizemode: "absolute", sizeref: L * 0.33,
      colorscale: [[0, color], [1, color]], showscale: false,
      legendgroup: name, showlegend: false, hoverinfo: "skip",
      lighting: { ambient: 0.95, diffuse: 0.4, specular: 0.1 },
    };
    return [shaft, cone];
  }

  // ---------- build the Plotly figure (body + gene cloud + sperm) ----------
  function render() {
    const s = state.scene;
    if (!s) return;
    const zs = s.z_scale || 7.0;
    const traces = [];

    // 1) Embryo body — one translucent mesh per segmentation region.
    for (const lbl of [...s.mask_labels].sort((a, b) => a - b)) {
      const mesh = s.region_meshes[String(lbl)];
      if (!mesh) continue;
      const def = (s.region_defaults || {})[String(lbl)] || { color: "#cccccc", opacity: 0.15 };
      const v = mesh.verts, f = mesh.faces;
      const nV = v.length / 3, nF = f.length / 3;
      const x = new Array(nV), y = new Array(nV), z = new Array(nV);
      for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
      const ii = new Array(nF), jj = new Array(nF), kk = new Array(nF);
      for (let i = 0; i < nF; i++) { ii[i] = f[i * 3]; jj[i] = f[i * 3 + 1]; kk[i] = f[i * 3 + 2]; }
      traces.push({
        type: "mesh3d", x, y, z, i: ii, j: jj, k: kk,
        color: def.color, opacity: def.opacity,
        name: `body M${lbl}`, showlegend: true, flatshading: false,
        hoverinfo: "skip",
        lighting: { ambient: 0.6, diffuse: 0.6, specular: 0.2, roughness: 0.9 },
        legendrank: lbl,
      });
    }

    // 2) Point cloud for the selected gene.
    const gene = geneSelect.value;
    const t = s.transcripts[gene];
    if (t && t.x.length) {
      const n = t.x.length;
      const gz = new Array(n);
      for (let i = 0; i < n; i++) gz[i] = t.gz[i] * zs;
      traces.push({
        type: "scatter3d", mode: "markers",
        name: `${gene}  (n=${n.toLocaleString()})`,
        x: t.x, y: t.y, z: gz,
        marker: { size: state.dotSize, color: GENE_COLOR, opacity: 0.85, line: { width: 0 } },
        hovertemplate: `<b>${gene}</b><br>x=%{x:.0f}, y=%{y:.0f}<extra></extra>`,
        legendrank: 20000,
      });
    }

    // 3) Sperm location — real (green diamond).
    const sp = s.sperm;
    traces.push({
      type: "scatter3d", mode: "markers", name: "Sperm (real)",
      x: [sp.x], y: [sp.y], z: [sp.z * zs],
      marker: {
        size: s.sperm_size || 11, color: REAL_COLOR,
        symbol: "diamond", opacity: 1,
        line: { width: 1.5, color: "#ffffff" },
      },
      hovertemplate: `<b>Sperm (real)</b><br>x=%{x:.0f}, y=%{y:.0f}, ` +
        `z=${sp.z} frame<extra></extra>`,
      legendrank: 30000,
    });

    // 3b) Null sperm — random COM direction ∩ shell (red diamond, toggleable).
    if (state.showNull && state.nullPoint) {
      const P = state.nullPoint;
      traces.push({
        type: "scatter3d", mode: "markers", name: "Sperm (null)",
        x: [P[0]], y: [P[1]], z: [P[2]],
        marker: { size: s.sperm_size || 11, color: NULL_COLOR, symbol: "diamond",
                  opacity: 1, line: { width: 1.5, color: "#ffffff" } },
        hovertemplate: "<b>Null sperm</b> (random direction)<extra></extra>",
        legendrank: 30500,
      });
    }

    // 4) Analysis vectors (pre-computed) — drawn as 3-D arrows if toggled on.
    const A = s.analysis;
    if (A) {
      const gv = A.genes[gene];
      const byKey = Object.fromEntries(VECTORS.map((v) => [v.key, v]));
      // embryo COM marker — the target the COM-directed arrows point at.
      if (state.vectors.gene2emb || state.vectors.sperm2emb) {
        const c = A.embryo_com;
        traces.push({
          type: "scatter3d", mode: "markers", name: "Embryo COM",
          x: [c[0]], y: [c[1]], z: [c[2]],
          marker: { size: 5.5, color: "#334155", opacity: 0.95,
                    line: { width: 1, color: "#ffffff" } },
          hovertemplate: "Embryo COM<extra></extra>", legendrank: 39000,
        });
      }
      if (state.vectors.pca && gv && gv.pca)
        traces.push(...arrowTraces(gv.com, gv.pca, A.arrow_len, zs, byKey.pca.color, byKey.pca.label));
      if (state.vectors.gene2emb && gv)
        traces.push(...arrowTraces(gv.com, gv.g2e, A.arrow_len, zs, byKey.gene2emb.color, byKey.gene2emb.label));
      if (state.vectors.sperm2emb)
        traces.push(...arrowTraces(A.sperm_plot, A.sperm_to_emb, A.arrow_len, zs, byKey.sperm2emb.color, byKey.sperm2emb.label));
    }

    // 5) Predicted sperm location (from the ranked-gene axis→sperm mapping).
    const pr = state.prediction;
    if (state.predict.enabled && pr && pr.point) {
      const P = pr.point, C = pr.com;
      traces.push({
        type: "scatter3d", mode: "lines", name: "Predicted axis",
        x: [C[0], P[0]], y: [C[1], P[1]], z: [C[2], P[2]],
        line: { color: REAL_COLOR, width: 5, dash: "dot" },
        hoverinfo: "skip", showlegend: false, legendrank: 50000,
      });
      traces.push({
        type: "scatter3d", mode: "markers", name: "Predicted sperm",
        x: [P[0]], y: [P[1]], z: [P[2]],
        marker: { size: 15, color: REAL_COLOR, symbol: "circle-open",
                  line: { width: 3, color: REAL_COLOR }, opacity: 1 },
        hovertemplate: "<b>Predicted sperm</b>" +
          (pr.errDeg != null ? `<br>error ${pr.errDeg.toFixed(0)}° · ${pr.used} genes` : "") +
          "<extra></extra>",
        legendrank: 50001,
      });
    }

    Plotly.react(plotHost, traces, sceneLayout(s), {
      responsive: true, displaylogo: false,
      modeBarButtonsToRemove: ["tableRotation", "resetCameraLastSave3d", "hoverClosest3d"],
      toImageButtonOptions: { format: "png", scale: 2 },
    });
    updateReadout(s, gene, t ? t.x.length : 0);
  }

  // Scene box locked to the data extents (+5% pad); axes hidden; physical
  // aspect from extents. uirevision keyed to the embryo id keeps the camera
  // while changing genes, but resets it when switching embryos.
  function sceneLayout(s) {
    const ex = s.extents;
    const pad = (r) => { const p = (r[1] - r[0]) * 0.05 || 1; return [r[0] - p, r[1] + p]; };
    const rx = pad(ex.x), ry = pad(ex.y), rz = pad(ex.z);
    const sx = rx[1] - rx[0], sy = ry[1] - ry[0], sz = rz[1] - rz[0];
    const m = Math.max(sx, sy, sz) || 1;
    return {
      template: "plotly_white",
      paper_bgcolor: "rgba(0,0,0,0)",
      margin: { l: 0, r: 0, t: 0, b: 0 },
      autosize: true,
      showlegend: true,
      legend: {
        itemsizing: "constant", font: { size: 12 },
        bgcolor: "rgba(255,255,255,0.82)", bordercolor: "#e7e9ef", borderwidth: 1,
        x: 0.99, xanchor: "right", y: 0.98, yanchor: "top",
      },
      scene: {
        xaxis: { visible: false, range: rx, autorange: false },
        yaxis: { visible: false, range: ry, autorange: false },
        zaxis: { visible: false, range: rz, autorange: false },
        aspectmode: "manual",
        aspectratio: { x: sx / m, y: sy / m, z: sz / m },
        uirevision: s.id,
        camera: { eye: { x: 1.5, y: 1.5, z: 1.15 } },
      },
    };
  }

  function updateReadout(s, gene, n) {
    const inPanel = s.genes.includes(gene);
    readoutEl.innerHTML =
      `<span><b>${s.title}</b></span>` +
      (inPanel
        ? `<span>${gene}: <b>${n.toLocaleString()}</b> transcripts</span>`
        : `<span>${gene}: <b>—</b> <span class="muted-note">not in this embryo</span></span>`) +
      `<span class="sperm-chip">sperm · segment ${s.sperm.segment || "?"}</span>`;
  }

  // =====================================================================
  // Analysis plots
  // =====================================================================
  // Safe (re)render into a small plot div. Only clears innerHTML when the div
  // is NOT already a Plotly graph — manually blanking a live Plotly div corrupts
  // its internal state and blanks subsequent Plotly.react() updates.
  function plotInto(div, traces, layout) {
    if (!div.classList.contains("js-plotly-plot")) div.innerHTML = "";
    Plotly.react(div, traces, layout, { responsive: true, displayModeBar: false });
  }
  function plotEmpty(div, msg, pad) {
    Plotly.purge(div);
    div.innerHTML = `<div class="violin-empty"${pad ? ` style="padding:${pad}"` : ""}>${msg}</div>`;
  }

  // Refresh the gene-dependent plots: PCA scree (current embryo+gene) and,
  // if the drawer is open, the two cross-embryo violins (gene only).
  function updateAnalysisPlots() {
    const gene = geneSelect.value;
    if (drawerGeneEl) drawerGeneEl.textContent = gene ? "· " + gene : "";
    renderPcaPlot();
    if (state.drawerOpen) renderViolins();
    updateRankingHighlight();
  }

  // ---- PCA scree bar (PC1/PC2/PC3 variance explained) for embryo+gene ----
  function renderPcaPlot() {
    const s = state.scene;
    const gene = geneSelect.value;
    pcaGeneEl.textContent = gene ? "· " + gene : "";
    const inPanel = s && s.genes.includes(gene);
    const gv = s && s.analysis && s.analysis.genes[gene];
    const evr = gv && gv.evr;
    if (!evr) {
      plotEmpty(pcaPlot, inPanel ? "Fewer than 2 transcripts — no PCA"
                                 : "Gene not in this embryo", "38px 4px");
      return;
    }
    const pct = evr.map((e) => +(e * 100).toFixed(1));
    const trace = {
      type: "bar", x: ["PC1", "PC2", "PC3"], y: pct,
      marker: { color: ["#7c3aed", "#c4b5fd", "#e4dcff"], line: { width: 0 } },
      text: pct.map((p) => p + "%"), textposition: "outside",
      textfont: { size: 10, color: "#3c4453" }, cliponaxis: false,
      hovertemplate: "%{x}: %{y:.1f}%<extra></extra>",
    };
    const layout = {
      margin: { l: 30, r: 6, t: 12, b: 20 }, height: 116,
      yaxis: { range: [0, 105], ticksuffix: "%", tickfont: { size: 9 },
               gridcolor: "#eef1f5", fixedrange: true, dtick: 25 },
      xaxis: { tickfont: { size: 11 }, fixedrange: true },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      font: { color: "#1a2233" }, showlegend: false, bargap: 0.35,
    };
    plotInto(pcaPlot, [trace], layout);
  }

  // ---- cross-embryo |cos θ| series for the selected gene ----
  // dirOf(embryo) returns the sperm-axis unit vector for that embryo (real =
  // sperm_to_emb, null = the embryo's random null direction).
  function dotAbs(a, b) { return Math.abs(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]); }
  const realDirOf = (e) => e.sperm_to_emb;
  const nullDirOf = (e) => state.nullDirs[e.id];
  function violinSeries(gene, key, dirOf) {
    const out = { vals: [], labels: [], colors: [], cur: null, curLabel: null };
    if (!state.index) return out;
    for (const e of state.index.embryos) {
      const gv = e.genes[gene];
      if (!gv) continue;
      const v = gv[key];
      if (!v) continue;                          // PC1 may be null (n<2)
      const dir = dirOf(e);
      if (!dir) continue;
      const d = dotAbs(v, dir);
      out.vals.push(+d.toFixed(4));
      out.labels.push(`${e.label} · ${e.stage}`);
      out.colors.push(STAGE_COLOR[e.stage] || "#888");
      if (e.id === state.currentId) { out.cur = d; out.curLabel = e.label; }
    }
    return out;
  }

  // Two violins per plot: real (green) vs randomized null (red).
  function renderOneViolin(div, real, nul) {
    if (!real.vals.length && !nul.vals.length) {
      plotEmpty(div, "This gene is not present in any embryo panel.");
      return;
    }
    const mk = (data, color, name) => ({
      type: "violin", y: data.vals, x: data.vals.map(() => name), name,
      points: "all", pointpos: 0, jitter: 0.35, scalemode: "width", width: 0.7,
      box: { visible: true, width: 0.12 }, meanline: { visible: true },
      span: [0, 1], spanmode: "manual",
      line: { color, width: 1.5 }, fillcolor: hexA(color, 0.14),
      marker: { size: 6, color: data.colors, line: { width: 0.5, color: "#fff" }, opacity: 0.9 },
      customdata: data.labels,
      hovertemplate: `<b>${name}</b><br>%{customdata}<br>|cos θ| = %{y:.3f}<extra></extra>`,
    });
    const traces = [];
    if (real.vals.length) traces.push(mk(real, REAL_COLOR, "Real"));
    if (nul.vals.length) traces.push(mk(nul, NULL_COLOR, "Null"));
    const shapes = [], anns = [];
    const addCur = (data, color, name) => {
      if (data.cur == null) return;
      shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1, y0: data.cur, y1: data.cur,
                    line: { color, width: 1.5, dash: "dot" } });
      anns.push({ xref: "paper", x: 0.99, xanchor: "right", y: data.cur, yanchor: "bottom",
                  text: `${name} ${data.cur.toFixed(2)}`, showarrow: false,
                  font: { size: 9.5, color } });
    };
    addCur(real, REAL_COLOR, "this embryo · real");
    addCur(nul, NULL_COLOR, "null");
    const layout = {
      margin: { l: 40, r: 12, t: 6, b: 18 }, height: 200, violinmode: "group",
      yaxis: { range: [-0.02, 1.04], tickfont: { size: 9 }, gridcolor: "#eef1f5",
               fixedrange: true, dtick: 0.25, zeroline: false,
               title: { text: "|cos θ|", font: { size: 10 } } },
      xaxis: { type: "category", fixedrange: true, tickfont: { size: 11 } },
      shapes, annotations: anns,
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      font: { color: "#1a2233" }, showlegend: false,
    };
    plotInto(div, traces, layout);
  }

  function renderViolins() {
    const gene = geneSelect.value;
    if (!gene) return;
    const aR = violinSeries(gene, "pca", realDirOf), aN = violinSeries(gene, "pca", nullDirOf);
    const bR = violinSeries(gene, "g2e", realDirOf), bN = violinSeries(gene, "g2e", nullDirOf);
    violinPcaSub.textContent = `${gene} · n = ${aR.vals.length} · real vs null`;
    violinG2eSub.textContent = `${gene} · n = ${bR.vals.length} · real vs null`;
    renderOneViolin(violinPca, aR, aN);
    renderOneViolin(violinG2e, bR, bN);
  }

  function setDrawerOpen(open) {
    state.drawerOpen = open;
    drawer.dataset.open = open ? "true" : "false";
    drawerHandle.setAttribute("aria-expanded", String(open));
    if (open) {
      renderViolins();
      requestAnimationFrame(() => {
        try { Plotly.Plots.resize(violinPca); Plotly.Plots.resize(violinG2e); } catch (_) {}
      });
    }
  }

  // ---- resizable violin width (symmetric about the centered group) ----
  function resizeViolinPlots() {
    try { Plotly.Plots.resize(violinPca); Plotly.Plots.resize(violinG2e); } catch (_) {}
  }
  function applyViolinWidth(w) {
    const maxW = (drawerBody.clientWidth || window.innerWidth) - 32;
    w = Math.max(300, Math.min(w, maxW));
    violinWrap.style.width = w + "px";
    return w;
  }
  function loadViolinWidth() {
    const v = parseFloat(localStorage.getItem(VIOLIN_W_KEY));
    if (isFinite(v)) applyViolinWidth(v);
  }
  (function wireViolinResize() {
    let dragging = false, startX = 0, startW = 0, raf = 0;
    const liveResize = () => {
      if (raf) return;
      raf = requestAnimationFrame(() => { raf = 0; resizeViolinPlots(); });
    };
    violinResize.addEventListener("pointerdown", (e) => {
      dragging = true; startX = e.clientX;
      startW = violinWrap.getBoundingClientRect().width;
      violinResize.classList.add("dragging");
      violinResize.setPointerCapture(e.pointerId);
      e.preventDefault();
    });
    violinResize.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      applyViolinWidth(startW + 2 * (e.clientX - startX));   // right edge → ±both sides
      liveResize();
    });
    const end = (e) => {
      if (!dragging) return;
      dragging = false;
      violinResize.classList.remove("dragging");
      try { violinResize.releasePointerCapture(e.pointerId); } catch (_) {}
      try {
        localStorage.setItem(VIOLIN_W_KEY,
          String(Math.round(violinWrap.getBoundingClientRect().width)));
      } catch (_) {}
      resizeViolinPlots();
    };
    violinResize.addEventListener("pointerup", end);
    violinResize.addEventListener("pointercancel", end);
  })();

  // ---- generic pointer-drag resizer (applies synchronously per move) ----
  function wireDragResize(handle, { onStart, onDrag, onEnd }) {
    let dragging = false, sx = 0, sy = 0;
    handle.addEventListener("pointerdown", (e) => {
      dragging = true; sx = e.clientX; sy = e.clientY;
      handle.classList.add("dragging");
      try { handle.setPointerCapture(e.pointerId); } catch (_) {}
      e.preventDefault();
      if (onStart) onStart();
    });
    handle.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      onDrag(e.clientX - sx, e.clientY - sy);        // setting a CSS var is cheap
    });
    const end = (e) => {
      if (!dragging) return;
      dragging = false; handle.classList.remove("dragging");
      try { handle.releasePointerCapture(e.pointerId); } catch (_) {}
      if (onEnd) onEnd();
    };
    handle.addEventListener("pointerup", end);
    handle.addEventListener("pointercancel", end);
  }

  // ---- bottom drawer: variable height (grows upward from a pinned bottom) ----
  function applyDrawerHeight(h) {
    const maxH = window.innerHeight - 100;      // stay clear of the topbar + handle
    h = Math.max(140, Math.min(h, maxH));
    drawer.style.setProperty("--drawer-h", h + "px");
    return h;
  }
  (function () {
    let startH = 0;
    wireDragResize(drawerResize, {
      onStart: () => { startH = drawerBody.getBoundingClientRect().height; },
      onDrag: (dx, dy) => { applyDrawerHeight(startH - dy); },   // drag up → taller
      onEnd: () => {
        try { localStorage.setItem(DRAWER_H_KEY,
          String(Math.round(drawerBody.getBoundingClientRect().height))); } catch (_) {}
        resizeViolinPlots();
      },
    });
  })();

  // ---- right drawer: variable width (grows leftward from a pinned right) ----
  function applyRdrawerWidth(w) {
    const maxW = Math.min(window.innerWidth * 0.7, 640);
    w = Math.max(240, Math.min(w, maxW));
    rdrawer.style.setProperty("--rdrawer-w", w + "px");
    return w;
  }
  (function () {
    let startW = 0;
    wireDragResize(rdrawerResize, {
      onStart: () => { startW = rdrawer.getBoundingClientRect().width; },
      onDrag: (dx) => { applyRdrawerWidth(startW - dx); },        // drag left → wider
      onEnd: () => {
        try { localStorage.setItem(RDRAWER_W_KEY,
          String(Math.round(rdrawer.getBoundingClientRect().width))); } catch (_) {}
      },
    });
  })();

  function loadDrawerSizes() {
    const h = parseFloat(localStorage.getItem(DRAWER_H_KEY));
    if (isFinite(h)) applyDrawerHeight(h);
    const w = parseFloat(localStorage.getItem(RDRAWER_W_KEY));
    if (isFinite(w)) applyRdrawerWidth(w);
  }

  // =====================================================================
  // Gene ranking (right drawer)
  // =====================================================================
  function mean(a) { let s = 0; for (const x of a) s += x; return a.length ? s / a.length : 0; }
  function stdev(a) {
    if (a.length < 2) return 0;
    const m = mean(a);
    let s = 0; for (const x of a) s += (x - m) * (x - m);
    return Math.sqrt(s / a.length);
  }

  // For every gene, gather |cos θ| across all embryos (per dot product), then
  // rank ascending by σ (lowest variation = most consistent = top). Genes must
  // appear in ≥ RANK_MIN_N embryos to be ranked. Computed once (index is static).
  // dirOf(embryo) supplies the sperm axis (real or null), so the same ranking
  // logic serves both conditions.
  function computeRankings(dirOf) {
    const acc = new Map();          // gene -> { pca:[], g2e:[], pcaN:[], g2eN:[] }
    const g2emb = new Map();        // gene -> [embryoId,...]
    for (const e of state.index.embryos) {
      const dir = dirOf(e);
      if (!dir) continue;
      for (const g in e.genes) {
        const gv = e.genes[g];
        if (!acc.has(g)) acc.set(g, { pca: [], g2e: [], pcaN: [], g2eN: [] });
        const a = acc.get(g);
        const cnt = gv.n || 0;      // transcripts of this gene in this embryo
        if (gv.pca) { a.pca.push(dotAbs(gv.pca, dir)); a.pcaN.push(cnt); }
        if (gv.g2e) { a.g2e.push(dotAbs(gv.g2e, dir)); a.g2eN.push(cnt); }
        if (!g2emb.has(g)) g2emb.set(g, []);
        g2emb.get(g).push(e.id);
      }
    }
    state.geneToEmbryos = g2emb;    // gene presence — same for real / null
    const build = (key) => {
      const out = [];
      for (const [g, a] of acc) {
        const vals = a[key];
        if (vals.length >= RANK_MIN_N)
          out.push({ gene: g, n: vals.length, std: stdev(vals), mean: mean(vals),
                     // smallest transcript count across the embryos it appears in —
                     // used by the "> N in every embryo" filter.
                     minN: Math.min(...a[key + "N"]) });
      }
      out.sort((x, y) => x.std - y.std || x.gene.localeCompare(y.gene));
      return out;
    };
    return { pca: build("pca"), g2e: build("g2e") };
  }
  function computeAllRankings() {
    state.rankings = { real: computeRankings(realDirOf), null: computeRankings(nullDirOf) };
  }

  const RANK_COLOR = { real: REAL_COLOR, null: NULL_COLOR };
  // Ranked rows for a condition/axis (with the transcript filter applied).
  function rankedRows(cond, axis) {
    let rows = (state.rankings && state.rankings[cond] && state.rankings[cond][axis]) || [];
    const thr = state.rankFilter || 0;
    if (thr > 0) rows = rows.filter((r) => r.minN > thr);   // > N in EVERY embryo
    return rows;
  }
  // The ranked rows currently shown (active tab). Prediction always uses REAL.
  function currentRankedRows() { return rankedRows(state.rankCond, state.rankAxis); }

  function renderRankings() {
    if (!state.rankings) return;
    const rows = currentRankedRows();
    const cur = geneSelect.value;
    const thr = state.rankFilter || 0;
    if (!rows.length) {
      rankListEl.innerHTML = thr > 0
        ? `<div class="rank-empty">No genes have &gt; ${thr} transcripts in every embryo.</div>`
        : `<div class="rank-empty">No genes appear in ≥ ${RANK_MIN_N} embryos.</div>`;
      return;
    }
    const maxStd = Math.max(...rows.map((r) => r.std)) || 1;
    const color = RANK_COLOR[state.rankCond];
    rankListEl.innerHTML = rows.map((r, i) => {
      const pct = Math.round((r.std / maxStd) * 100);
      return `<div class="rank-row${r.gene === cur ? " current" : ""}" data-gene="${r.gene}"` +
        ` title="σ = ${r.std.toFixed(3)}  ·  mean |cosθ| = ${r.mean.toFixed(2)}  ·  n = ${r.n} embryos">` +
        `<span class="rank-num">${i + 1}</span>` +
        `<span class="rank-gene">${r.gene}</span>` +
        `<span class="rank-std">σ ${r.std.toFixed(2)}</span>` +
        `<span class="rank-bar"><span style="width:${pct}%;background:${color}"></span></span>` +
        `</div>`;
    }).join("");
  }

  function updateRankingHighlight() {
    const cur = geneSelect.value;
    rankListEl.querySelectorAll(".rank-row").forEach((row) =>
      row.classList.toggle("current", row.dataset.gene === cur));
  }

  // Select a gene from the ranking: if the current embryo has it, just switch;
  // otherwise jump to an embryo whose panel contains it so it's visible.
  async function selectGene(gene) {
    state.userGene = gene;
    const s = state.scene;
    if (s && s.genes.includes(gene)) {
      geneSelect.value = gene;
      render();
      updateAnalysisPlots();
    } else {
      const ids = state.geneToEmbryos && state.geneToEmbryos.get(gene);
      const target = ids && state.manifest.find((m) => ids.includes(m.id));
      if (target) await selectEmbryo(target.id);
    }
    updateRankingHighlight();
  }

  // =====================================================================
  // Predicted sperm location
  // =====================================================================
  // Assume the sperm lies on the embryo cortex, so predicting it reduces to a
  // direction from the embryo centre of mass. For each ranked gene we map its
  // (PCA or gene→COM) axis to a sperm-direction estimate via the pre-computed
  // average rotation (Kabsch); the estimates are combined as an inverse-variance
  // weighted sum (lower σ → higher weight), normalized, and intersected with the
  // blastomere shell.
  //
  // --- small 3-D vector / 3×3 matrix helpers (matrices are row-major [9]) ---
  const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const vdot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const vcross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const vnorm = (a) => Math.hypot(a[0], a[1], a[2]);
  const vunit = (a) => { const n = vnorm(a); return n > 0 ? [a[0] / n, a[1] / n, a[2] / n] : [0, 0, 0]; };
  const mvec = (M, v) => [M[0] * v[0] + M[1] * v[1] + M[2] * v[2],
                          M[3] * v[0] + M[4] * v[1] + M[5] * v[2],
                          M[6] * v[0] + M[7] * v[1] + M[8] * v[2]];
  function mmul(A, B) {
    const C = new Array(9);
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
      let s = 0; for (let k = 0; k < 3; k++) s += A[i * 3 + k] * B[k * 3 + j];
      C[i * 3 + j] = s;
    }
    return C;
  }
  const mT = (A) => [A[0], A[3], A[6], A[1], A[4], A[7], A[2], A[5], A[8]];
  const mdet = (A) => A[0] * (A[4] * A[8] - A[5] * A[7]) - A[1] * (A[3] * A[8] - A[5] * A[6]) +
                      A[2] * (A[3] * A[7] - A[4] * A[6]);

  // Symmetric 3×3 eigendecomposition (Jacobi). Returns eigenvalues (descending)
  // and eigenvectors as columns.
  function eigSym3(Ain) {
    let A = [[Ain[0], Ain[1], Ain[2]], [Ain[3], Ain[4], Ain[5]], [Ain[6], Ain[7], Ain[8]]];
    let V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    const mm = (P, Q) => { const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) { let s = 0;
        for (let k = 0; k < 3; k++) s += P[i][k] * Q[k][j]; C[i][j] = s; } return C; };
    const tr = (P) => [[P[0][0], P[1][0], P[2][0]], [P[0][1], P[1][1], P[2][1]], [P[0][2], P[1][2], P[2][2]]];
    for (let it = 0; it < 100; it++) {
      let p = 0, q = 1, mx = Math.abs(A[0][1]);
      if (Math.abs(A[0][2]) > mx) { mx = Math.abs(A[0][2]); p = 0; q = 2; }
      if (Math.abs(A[1][2]) > mx) { mx = Math.abs(A[1][2]); p = 1; q = 2; }
      if (mx < 1e-14) break;
      const phi = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]);
      const c = Math.cos(phi), s = Math.sin(phi);
      const J = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
      J[p][p] = c; J[q][q] = c; J[p][q] = s; J[q][p] = -s;
      A = mm(mm(tr(J), A), J);
      V = mm(V, J);
    }
    const vals = [A[0][0], A[1][1], A[2][2]];
    const cols = [[V[0][0], V[1][0], V[2][0]], [V[0][1], V[1][1], V[2][1]], [V[0][2], V[1][2], V[2][2]]];
    const idx = [0, 1, 2].sort((i, j) => vals[j] - vals[i]);
    return { values: idx.map((i) => vals[i]), vecs: idx.map((i) => cols[i]) };
  }

  // Kabsch rotation from B = Σ gᵢ sᵢᵀ (row-major [9]): R minimizing Σ‖R g − s‖².
  // SVD B = U Σ Vᵀ (V from eig of BᵀB, U = B V Σ⁻¹); R = V diag(1,1,d) Uᵀ.
  function kabschFromB(B) {
    const M = mmul(mT(B), B);                 // BᵀB (symmetric PSD)
    const { values, vecs } = eigSym3(M);
    const sig = values.map((x) => Math.sqrt(Math.max(x, 0)));
    const Vc = vecs;                          // right singular vectors (columns)
    const Uc = [null, null, null];
    for (let k = 0; k < 3; k++) if (sig[k] > 1e-9) {
      const bv = mvec(B, Vc[k]); Uc[k] = [bv[0] / sig[k], bv[1] / sig[k], bv[2] / sig[k]];
    }
    for (let k = 0; k < 3; k++) if (!Uc[k]) {   // orthonormal completion for σ≈0
      const known = Uc.filter(Boolean);
      if (known.length === 2) Uc[k] = vunit(vcross(known[0], known[1]));
      else if (known.length === 1) {
        const a = known[0], t = Math.abs(a[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
        Uc[k] = vunit(vcross(a, t));
      } else Uc[k] = [k === 0 ? 1 : 0, k === 1 ? 1 : 0, k === 2 ? 1 : 0];
    }
    const U = [Uc[0][0], Uc[1][0], Uc[2][0], Uc[0][1], Uc[1][1], Uc[2][1], Uc[0][2], Uc[1][2], Uc[2][2]];
    const V = [Vc[0][0], Vc[1][0], Vc[2][0], Vc[0][1], Vc[1][1], Vc[2][1], Vc[0][2], Vc[1][2], Vc[2][2]];
    const d = mdet(mmul(V, mT(U))) >= 0 ? 1 : -1;
    const Vd = [V[0], V[1], V[2] * d, V[3], V[4], V[5] * d, V[6], V[7], V[8] * d];  // V·diag(1,1,d)
    return mmul(Vd, mT(U));
  }

  // Möller–Trumbore ray/triangle; returns forward t (>ε) or null.
  function rayTri(o, dir, a, b, c) {
    const e1 = vsub(b, a), e2 = vsub(c, a), pv = vcross(dir, e2), det = vdot(e1, pv);
    if (Math.abs(det) < 1e-9) return null;
    const inv = 1 / det, tv = vsub(o, a);
    const u = vdot(tv, pv) * inv; if (u < -1e-6 || u > 1 + 1e-6) return null;
    const qv = vcross(tv, e1), v = vdot(dir, qv) * inv; if (v < -1e-6 || u + v > 1 + 1e-6) return null;
    const t = vdot(e2, qv) * inv; return t > 1e-4 ? t : null;
  }

  // Intersect the ray (origin=COM, dir) with the blastomere shell (segments 1 & 2);
  // return the FARTHEST forward hit (the outer cortex). Falls back to the extent
  // of seg-1/2 vertices projected on the ray if no triangle is hit.
  function rayShellHit(o, dir, meshes) {
    let maxT = -Infinity;
    for (const lbl of ["1", "2"]) {
      const m = meshes && meshes[lbl]; if (!m) continue;
      const v = m.verts, f = m.faces;
      for (let fi = 0; fi < f.length; fi += 3) {
        const i0 = f[fi] * 3, i1 = f[fi + 1] * 3, i2 = f[fi + 2] * 3;
        const t = rayTri(o, dir, [v[i0], v[i0 + 1], v[i0 + 2]],
          [v[i1], v[i1 + 1], v[i1 + 2]], [v[i2], v[i2 + 1], v[i2 + 2]]);
        if (t != null && t > maxT) maxT = t;
      }
    }
    if (!isFinite(maxT) || maxT <= 0) {          // fallback: max projection on ray
      maxT = -Infinity;
      for (const lbl of ["1", "2"]) {
        const m = meshes && meshes[lbl]; if (!m) continue;
        const v = m.verts;
        for (let i = 0; i < v.length; i += 3) {
          const proj = (v[i] - o[0]) * dir[0] + (v[i + 1] - o[1]) * dir[1] + (v[i + 2] - o[2]) * dir[2];
          if (proj > maxT) maxT = proj;
        }
      }
      if (!isFinite(maxT) || maxT <= 0) return null;
    }
    return [o[0] + dir[0] * maxT, o[1] + dir[1] * maxT, o[2] + dir[2] * maxT];
  }

  const PRED_EPS = 0.05;   // inverse-variance floor: w = 1/(σ² + ε²)

  function computePrediction() {
    state.prediction = null;
    const s = state.scene, idx = state.index;
    if (!s || !s.analysis || !idx || !idx.mappings) return;
    const axis = state.rankAxis;                // "pca" | "g2e" (from active tab)
    const rows = rankedRows("real", axis);      // prediction is always vs REAL sperm
    const wantAll = state.predict.topN === "all";
    const N = wantAll ? Infinity : parseInt(state.predict.topN, 10);
    const ste = s.analysis.sperm_to_emb;
    const sE = ste ? [-ste[0], -ste[1], -ste[2]] : null;   // current COM→sperm (µm)

    // Walk the ranking (most consistent first) and take the top-N genes that are
    // actually PRESENT in this embryo (gene panels differ across runs, so the
    // globally top-ranked gene often isn't in a given embryo).
    let acc = [0, 0, 0], totW = 0, used = 0;
    for (const r of rows) {
      if (used >= N) break;
      const gv = s.analysis.genes[r.gene];
      const gE = gv && gv[axis];                // current embryo's gene axis (µm)
      const m = idx.mappings[r.gene] && idx.mappings[r.gene][axis];
      if (!gE || !m) continue;
      let R;
      if (state.predict.loo && sE) {
        const B = m.B, Bloo = [
          B[0] - gE[0] * sE[0], B[1] - gE[0] * sE[1], B[2] - gE[0] * sE[2],
          B[3] - gE[1] * sE[0], B[4] - gE[1] * sE[1], B[5] - gE[1] * sE[2],
          B[6] - gE[2] * sE[0], B[7] - gE[2] * sE[1], B[8] - gE[2] * sE[2],
        ];
        R = kabschFromB(Bloo);
      } else {
        R = m.R;
      }
      const sHat = vunit(mvec(R, gE));          // predicted COM→sperm from this gene
      if (vnorm(sHat) === 0) continue;
      const w = 1 / (r.std * r.std + PRED_EPS * PRED_EPS);
      acc = [acc[0] + w * sHat[0], acc[1] + w * sHat[1], acc[2] + w * sHat[2]];
      totW += w; used++;
    }
    if (used === 0 || totW === 0) { state.prediction = { used: 0 }; return; }
    const resLen = vnorm(acc);
    if (resLen < 1e-9) { state.prediction = { used, weak: true }; return; }
    const dirUm = [acc[0] / resLen, acc[1] / resLen, acc[2] / resLen];    // µm unit
    const zs = s.z_scale || 7.0;
    const dirPlot = vunit([dirUm[0] / 0.15, dirUm[1] / 0.15, dirUm[2] * zs]);
    const com = s.analysis.embryo_com;
    const point = rayShellHit(com, dirPlot, s.region_meshes);
    let errDeg = null;
    if (sE) { const c = Math.max(-1, Math.min(1, vdot(dirUm, vunit(sE)))); errDeg = Math.acos(c) * 180 / Math.PI; }
    state.prediction = { used, point, com, confidence: resLen / totW, errDeg, axis };
  }

  function updatePrediction() {
    computePrediction();
    updatePredictReadout();
    if (state.scene) render();
  }

  const AXIS_LABEL = { pca: "PC1", g2e: "Gene→COM" };
  function updatePredictReadout() {
    if (!predictReadout) return;
    if (!state.predict.enabled) { predictReadout.innerHTML = ""; return; }
    const pr = state.prediction;
    const axisTxt = AXIS_LABEL[state.rankAxis] || state.rankAxis;
    if (!pr || !pr.used) {
      predictReadout.innerHTML =
        `<span class="predict-warn">No ranked genes for this axis are present in this ` +
        `embryo — can't predict.</span>`;
      return;
    }
    const conf = Math.round((pr.confidence || 0) * 100);
    const err = pr.errDeg != null ? `${pr.errDeg.toFixed(0)}°` : "—";
    predictReadout.innerHTML =
      `<div><b>${pr.used}</b> genes · axis <b>${axisTxt}</b>` +
      (state.predict.loo ? " · LOO" : "") + `</div>` +
      `<div>angular error <b>${err}</b> · agreement <b>${conf}%</b></div>`;
  }

  // =====================================================================
  // Null sperm condition (per-embryo random direction from the COM)
  // =====================================================================
  // A uniform random unit direction on the sphere (µm space).
  function randUnitDir() {
    const z = 2 * Math.random() - 1;
    const th = 2 * Math.PI * Math.random();
    const r = Math.sqrt(Math.max(0, 1 - z * z));
    return [r * Math.cos(th), r * Math.sin(th), z];
  }
  // Draw a fresh, independent null direction for every embryo.
  function generateNullDirs() {
    const d = {};
    for (const e of (state.index ? state.index.embryos : [])) d[e.id] = randUnitDir();
    state.nullDirs = d;
    try { localStorage.setItem(NULL_DIRS_KEY, JSON.stringify(d)); } catch (_) {}
  }
  // Shell intersection of the current embryo's null direction (plot space).
  function computeNullPoint() {
    state.nullPoint = null;
    const s = state.scene;
    if (!s || !s.analysis) return;
    const nd = state.nullDirs[s.id];
    if (!nd) return;
    const zs = s.z_scale || 7.0;
    const dirPlot = vunit([nd[0] / 0.15, nd[1] / 0.15, nd[2] * zs]);
    state.nullPoint = rayShellHit(s.analysis.embryo_com, dirPlot, s.region_meshes);
  }
  // Re-draw all null directions and refresh everything that depends on them.
  function regenerateNull() {
    generateNullDirs();
    if (state.rankings) state.rankings.null = computeRankings(nullDirOf);
    computeNullPoint();
    if (state.rankCond === "null") renderRankings();
    if (state.drawerOpen) renderViolins();
    if (state.scene) render();
  }

  // ---------- ui helpers ----------
  function showLoading(txt) { loadingTxt.textContent = txt; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) {
    placeholder.hidden = false;
    placeholder.innerHTML =
      `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div>` +
      `<div class="ph-sub">${msg}</div></div>`;
  }

  drawerHandle.addEventListener("click", () => setDrawerOpen(!state.drawerOpen));

  // Right drawer (gene ranking): open/close, tab switch, row → select gene.
  rdrawerHandle.addEventListener("click", () => {
    const open = rdrawer.dataset.open !== "true";
    rdrawer.dataset.open = open ? "true" : "false";
    rdrawerHandle.setAttribute("aria-expanded", String(open));
  });
  rtabsEl.querySelectorAll(".rtab").forEach((b) =>
    b.addEventListener("click", () => {
      state.rankAxis = b.dataset.axis;
      state.rankCond = b.dataset.cond;
      rtabsEl.querySelectorAll(".rtab").forEach((x) => x.classList.toggle("active", x === b));
      renderRankings();
      updatePrediction();          // prediction axis follows the active tab (always real)
    }));
  rankListEl.addEventListener("click", (e) => {
    const row = e.target.closest(".rank-row");
    if (row) selectGene(row.dataset.gene);
  });
  rankFilterEl.addEventListener("change", () => {
    state.rankFilter = parseInt(rankFilterEl.value, 10) || 0;
    renderRankings();
    updatePrediction();            // the ranked set changed
  });

  // ---------- prediction controls ----------
  function syncPredictControls() {
    predictEnable.checked = state.predict.enabled;
    predictTopN.value = state.predict.topN;
    predictLoo.checked = state.predict.loo;
    predictConfig.hidden = !state.predict.enabled;
  }
  predictEnable.addEventListener("change", () => {
    state.predict.enabled = predictEnable.checked;
    predictConfig.hidden = !state.predict.enabled;
    savePredictCfg();
    updatePrediction();
  });
  predictTopN.addEventListener("change", () => {
    state.predict.topN = predictTopN.value; savePredictCfg(); updatePrediction();
  });
  predictLoo.addEventListener("change", () => {
    state.predict.loo = predictLoo.checked; savePredictCfg(); updatePrediction();
  });

  // ---------- null-sperm controls (in the floating window) ----------
  function syncNullControls() { if (nullShowEl) nullShowEl.checked = state.showNull; }
  if (nullShowEl) nullShowEl.addEventListener("change", () => {
    state.showNull = nullShowEl.checked;
    try { localStorage.setItem(NULL_SHOW_KEY, state.showNull ? "1" : "0"); } catch (_) {}
    if (state.scene) render();
  });
  if (nullRegenEl) nullRegenEl.addEventListener("click", () => regenerateNull());

  // ---------- floating control window: drag (header) + resize (corners) ----------
  const CONTROLS_BOX_KEY = "sperm_viewer_controls_box";
  const resizePcaPlot = () => { try { Plotly.Plots.resize(pcaPlot); } catch (_) {} };
  function stageBox() { return controlsEl.parentElement.getBoundingClientRect(); }
  function saveControlsBox() {
    try {
      localStorage.setItem(CONTROLS_BOX_KEY, JSON.stringify({
        left: parseFloat(controlsEl.style.left), top: parseFloat(controlsEl.style.top),
        width: parseFloat(controlsEl.style.width), height: parseFloat(controlsEl.style.height),
      }));
    } catch (_) {}
  }
  function loadControlsBox() {
    let b; try { b = JSON.parse(localStorage.getItem(CONTROLS_BOX_KEY) || "null"); } catch (_) {}
    if (!b) return;
    if (isFinite(b.width)) controlsEl.style.width = b.width + "px";
    if (isFinite(b.height)) controlsEl.style.height = b.height + "px";
    if (isFinite(b.left)) { controlsEl.style.left = b.left + "px"; controlsEl.style.right = "auto"; }
    if (isFinite(b.top)) controlsEl.style.top = b.top + "px";
  }
  // Drag by the header.
  (function () {
    let start = null;
    controlsHeader.addEventListener("pointerdown", (e) => {
      const r = controlsEl.getBoundingClientRect(), st = stageBox();
      start = { x: e.clientX, y: e.clientY, left: r.left - st.left, top: r.top - st.top };
      controlsHeader.setPointerCapture(e.pointerId); e.preventDefault();
    });
    controlsHeader.addEventListener("pointermove", (e) => {
      if (!start) return;
      const st = stageBox();
      let left = start.left + (e.clientX - start.x), top = start.top + (e.clientY - start.y);
      left = Math.max(0, Math.min(left, st.width - controlsEl.offsetWidth));
      top = Math.max(0, Math.min(top, st.height - 30));      // keep header reachable
      controlsEl.style.left = left + "px"; controlsEl.style.top = top + "px"; controlsEl.style.right = "auto";
    });
    const end = (e) => { if (start) { start = null; try { controlsHeader.releasePointerCapture(e.pointerId); } catch (_) {} saveControlsBox(); } };
    controlsHeader.addEventListener("pointerup", end);
    controlsHeader.addEventListener("pointercancel", end);
  })();
  // Resize by any of the four corners (opposite edges stay pinned).
  const MINW = 240, MINH = 130;
  controlsEl.querySelectorAll(".rz").forEach((h) => {
    const cfg = { nw: [1, 1], ne: [0, 1], sw: [1, 0], se: [0, 0] }[h.dataset.corner];
    const [leftEdge, topEdge] = cfg;
    let start = null;
    h.addEventListener("pointerdown", (e) => {
      const r = controlsEl.getBoundingClientRect(), st = stageBox();
      const left = r.left - st.left, top = r.top - st.top;
      start = { x: e.clientX, y: e.clientY, w: r.width, h: r.height, left, top,
                right: left + r.width, bottom: top + r.height };
      h.setPointerCapture(e.pointerId); e.preventDefault(); e.stopPropagation();
    });
    h.addEventListener("pointermove", (e) => {
      if (!start) return;
      const dx = e.clientX - start.x, dy = e.clientY - start.y;
      let w, hgt, left = start.left, top = start.top;
      if (leftEdge) { w = Math.max(MINW, start.w - dx); left = start.right - w; }
      else { w = Math.max(MINW, start.w + dx); }
      if (topEdge) { hgt = Math.max(MINH, start.h - dy); top = start.bottom - hgt; }
      else { hgt = Math.max(MINH, start.h + dy); }
      controlsEl.style.width = w + "px"; controlsEl.style.height = hgt + "px";
      controlsEl.style.left = left + "px"; controlsEl.style.top = top + "px"; controlsEl.style.right = "auto";
      resizePcaPlot();
    });
    const end = (e) => { if (start) { start = null; try { h.releasePointerCapture(e.pointerId); } catch (_) {} saveControlsBox(); resizePcaPlot(); } };
    h.addEventListener("pointerup", end);
    h.addEventListener("pointercancel", end);
  });

  window.addEventListener("resize", () => { if (state.drawerOpen) resizeViolinPlots(); });

  // ---------- boot ----------
  (async function init() {
    try {
      const m = await (await fetch("data/manifest.json")).json();
      state.manifest = m.embryos;
      state.manifest.forEach((e) => state.byId.set(e.id, e));
      countEl.textContent = `${m.embryos.length} embryos`;
      buildVectorToolbar();
      buildTabs(m.embryos, m.stage_order || ["Zygote", "2-cell (early)", "2-cell (late)"]);
      loadViolinWidth();
      loadDrawerSizes();
      loadControlsBox();
      if (V && V.addWindowExtras) vcExtras = V.addWindowExtras($("#controls-body"), { defaultSize: state.dotSize, onDotSize: (s) => { state.dotSize = s; render(); } });
      // Cross-embryo vectors for the violin plots + gene ranking (load once).
      state.index = await loadScene("data/analysis_index.json.gz").catch(() => null);
      if (state.index) {
        // ensure every embryo has a persisted null direction
        if (state.index.embryos.some((e) => !state.nullDirs[e.id])) generateNullDirs();
        computeAllRankings();       // real + null
        renderRankings();
        rdrawer.hidden = false;
      }
      syncPredictControls();
      syncNullControls();
    } catch (err) {
      showError("Failed to load manifest: " + (err.message || err));
    }
  })();
})();
