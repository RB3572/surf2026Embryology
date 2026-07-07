/* Sperm · Embryo 3D Viewer
 *
 * Minimal front-end: a top nav-bar of the 45 sperm-positive embryos; selecting
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
  const tabsEl     = $("#tabs");
  const controlsEl = $("#controls");
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

  const VIOLIN_W_KEY = "sperm_viewer_violin_w";
  const DRAWER_H_KEY = "sperm_viewer_drawer_h";
  const RDRAWER_W_KEY = "sperm_viewer_rdrawer_w";
  const RANK_MIN_N = 3;   // min embryos for a gene to be ranked (meaningful σ)

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
  const GENE_SIZE  = 2.6;

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
    rankTab: "pca",       // which ranking tab is active
    rankFilter: 0,        // "> N transcripts in every embryo" (0 = all)
    rankings: null,       // { pca: [...], g2e: [...] }  (computed once)
    geneToEmbryos: null,  // Map<gene, embryoId[]>  for jump-to-gene
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
      populateGenes(scene);
      controlsEl.hidden = false;
      placeholder.hidden = true;
      drawer.hidden = false;
      render();
      updateAnalysisPlots();
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
        marker: { size: GENE_SIZE, color: GENE_COLOR, opacity: 0.85, line: { width: 0 } },
        hovertemplate: `<b>${gene}</b><br>x=%{x:.0f}, y=%{y:.0f}<extra></extra>`,
        legendrank: 20000,
      });
    }

    // 3) Sperm location.
    const sp = s.sperm;
    traces.push({
      type: "scatter3d", mode: "markers", name: "Sperm",
      x: [sp.x], y: [sp.y], z: [sp.z * zs],
      marker: {
        size: s.sperm_size || 11, color: s.sperm_color || "#ff2d95",
        symbol: "diamond", opacity: 1,
        line: { width: 1.5, color: "#ffffff" },
      },
      hovertemplate: `<b>Sperm</b><br>x=%{x:.0f}, y=%{y:.0f}, ` +
        `z=${sp.z} frame<extra></extra>`,
      legendrank: 30000,
    });

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
  function dotAbs(a, b) { return Math.abs(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]); }
  function violinSeries(gene, key) {
    const out = { vals: [], labels: [], colors: [], cur: null, curLabel: null };
    if (!state.index) return out;
    for (const e of state.index.embryos) {
      const gv = e.genes[gene];
      if (!gv) continue;
      const v = gv[key];
      if (!v) continue;                          // PC1 may be null (n<2)
      const d = dotAbs(v, e.sperm_to_emb);
      out.vals.push(+d.toFixed(4));
      out.labels.push(`${e.label} · ${e.stage}`);
      out.colors.push(STAGE_COLOR[e.stage] || "#888");
      if (e.id === state.currentId) { out.cur = d; out.curLabel = e.label; }
    }
    return out;
  }

  function renderOneViolin(div, data, color) {
    if (!data.vals.length) {
      plotEmpty(div, "This gene is not present in any embryo panel.");
      return;
    }
    const violin = {
      type: "violin", y: data.vals, name: "", orientation: "v",
      points: "all", pointpos: 0, jitter: 0.35, scalemode: "width", width: 0.7,
      box: { visible: true, width: 0.12 }, meanline: { visible: true },
      span: [0, 1], spanmode: "manual",
      line: { color, width: 1.5 }, fillcolor: hexA(color, 0.14),
      marker: { size: 6, color: data.colors, line: { width: 0.5, color: "#fff" }, opacity: 0.9 },
      customdata: data.labels,
      hovertemplate: "%{customdata}<br>|cos θ| = %{y:.3f}<extra></extra>",
    };
    const shapes = [], anns = [];
    if (data.cur != null) {
      shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1,
                    y0: data.cur, y1: data.cur,
                    line: { color: "#0b0d13", width: 1.5, dash: "dot" } });
      anns.push({ xref: "paper", x: 0.99, xanchor: "right", y: data.cur, yanchor: "bottom",
                  text: `this embryo (${data.curLabel}) = ${data.cur.toFixed(2)}`,
                  showarrow: false, font: { size: 10, color: "#0b0d13" } });
    }
    const layout = {
      margin: { l: 40, r: 12, t: 6, b: 14 }, height: 200,
      yaxis: { range: [-0.02, 1.04], tickfont: { size: 9 }, gridcolor: "#eef1f5",
               fixedrange: true, dtick: 0.25, zeroline: false,
               title: { text: "|cos θ|", font: { size: 10 } } },
      xaxis: { visible: false, fixedrange: true },
      shapes, annotations: anns,
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      font: { color: "#1a2233" }, showlegend: false,
    };
    plotInto(div, [violin], layout);
  }

  function renderViolins() {
    const gene = geneSelect.value;
    if (!gene) return;
    const a = violinSeries(gene, "pca");
    const b = violinSeries(gene, "g2e");
    violinPcaSub.textContent = `${gene} · n = ${a.vals.length} embryos`;
    violinG2eSub.textContent = `${gene} · n = ${b.vals.length} embryos`;
    renderOneViolin(violinPca, a, "#7c3aed");    // matches the PC1 arrow color
    renderOneViolin(violinG2e, b, "#2563eb");    // matches the Gene→COM arrow color
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
  function computeRankings() {
    const acc = new Map();          // gene -> { pca:[], g2e:[], pcaN:[], g2eN:[] }
    const g2emb = new Map();        // gene -> [embryoId,...]
    for (const e of state.index.embryos) {
      for (const g in e.genes) {
        const gv = e.genes[g];
        if (!acc.has(g)) acc.set(g, { pca: [], g2e: [], pcaN: [], g2eN: [] });
        const a = acc.get(g);
        const cnt = gv.n || 0;      // transcripts of this gene in this embryo
        if (gv.pca) { a.pca.push(dotAbs(gv.pca, e.sperm_to_emb)); a.pcaN.push(cnt); }
        if (gv.g2e) { a.g2e.push(dotAbs(gv.g2e, e.sperm_to_emb)); a.g2eN.push(cnt); }
        if (!g2emb.has(g)) g2emb.set(g, []);
        g2emb.get(g).push(e.id);
      }
    }
    state.geneToEmbryos = g2emb;
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

  const RANK_COLOR = { pca: "#7c3aed", g2e: "#2563eb" };
  function renderRankings() {
    if (!state.rankings) return;
    let rows = state.rankings[state.rankTab] || [];
    const thr = state.rankFilter || 0;
    if (thr > 0) rows = rows.filter((r) => r.minN > thr);   // > N in EVERY embryo
    const cur = geneSelect.value;
    if (!rows.length) {
      rankListEl.innerHTML = thr > 0
        ? `<div class="rank-empty">No genes have &gt; ${thr} transcripts in every embryo.</div>`
        : `<div class="rank-empty">No genes appear in ≥ ${RANK_MIN_N} embryos.</div>`;
      return;
    }
    const maxStd = Math.max(...rows.map((r) => r.std)) || 1;
    const color = RANK_COLOR[state.rankTab];
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
      state.rankTab = b.dataset.rtab;
      rtabsEl.querySelectorAll(".rtab").forEach((x) => x.classList.toggle("active", x === b));
      renderRankings();
    }));
  rankListEl.addEventListener("click", (e) => {
    const row = e.target.closest(".rank-row");
    if (row) selectGene(row.dataset.gene);
  });
  rankFilterEl.addEventListener("change", () => {
    state.rankFilter = parseInt(rankFilterEl.value, 10) || 0;
    renderRankings();
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
      // Cross-embryo vectors for the violin plots + gene ranking (load once).
      state.index = await loadScene("data/analysis_index.json.gz").catch(() => null);
      if (state.index) {
        state.rankings = computeRankings();
        renderRankings();
        rdrawer.hidden = false;
      }
    } catch (err) {
      showError("Failed to load manifest: " + (err.message || err));
    }
  })();
})();
