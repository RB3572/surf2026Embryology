/* Zygote Division Planes — analysis model.
 * Built on viewer-core.js (VCore): validated zygotes, the polar-body axis, 18 candidate
 * division planes, per-gene transcript split (blue/red) across the selected plane,
 * a counts+null chart, best-plane tables, and a cross-section plot. All statistics
 * are pre-computed (build_zygote.py); the UI only reads them.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const AXIS_C = "#111827", PLANE_C = "#f97316";
  const BLUE = "#0099a8", RED = "#f05a47", GREEN = "#16a34a";

  // Plotly.react corrupts if the div was cleared with innerHTML="" while it was
  // a live plot. Only clear non-Plotly content (e.g. an empty-state message).
  function plotInto(div, traces, layout, cfg) {
    if (!div.classList.contains("js-plotly-plot")) div.innerHTML = "";
    Plotly.react(div, traces, layout, cfg || { responsive: true, displayModeBar: false });
  }

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const geneSelect = $("#gene-select"), planeSelect = $("#plane-select");
  const axisShow = $("#axis-show"), planeShow = $("#plane-show"), allShow = $("#all-show");
  const circShow = $("#circ-show");
  const chartEl = $("#chart"), chartSub = $("#chart-sub"), chartReadout = $("#chart-readout");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), drawerBody = $("#drawer-body");
  const drawerEmb = $("#drawer-emb");
  const pcolorMode = $("#pcolor-mode");
  const xsPlane = $("#xs-plane"), xsNote = $("#xs-note");
  const xsCaption = $("#xs-caption"), xsBarSub = $("#xs-bar-sub");
  const xsOutlines = $("#xs-outlines"), xsBars = $("#xs-bars");
  const xsBody = $("#xs-body"), xsPb = $("#xs-pb"), xsPronuclei = $("#xs-pronuclei");
  const xsCrossLegend = $("#xs-cross-legend"), xsCrossDownload = $("#xs-cross-download");
  const xsCrossHighColor = $("#xs-cross-high-color"), xsCrossLowColor = $("#xs-cross-low-color");
  const xsCrossFillColor = $("#xs-cross-fill-color");
  const xsCrossBodyColor = $("#xs-cross-body-color"), xsCrossPbColor = $("#xs-cross-pb-color");
  const xsCrossPnColor = $("#xs-cross-pn-color"), xsCrossScale = $("#xs-cross-scale"), xsCrossFormat = $("#xs-cross-format");
  const xsBarLog = $("#xs-bar-log"), xsBarLegend = $("#xs-bar-legend"), xsBarNull = $("#xs-bar-null");
  const xsBarInterval = $("#xs-bar-interval"), xsBarGrid = $("#xs-bar-grid"), xsBarDownload = $("#xs-bar-download");
  const xsBarHighColor = $("#xs-bar-high-color"), xsBarLowColor = $("#xs-bar-low-color");
  const xsBarNullColor = $("#xs-bar-null-color"), xsBarScale = $("#xs-bar-scale"), xsBarFormat = $("#xs-bar-format");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle");
  const rtabsEl = $("#rtabs"), bestListEl = $("#best-list");

  const state = {
    manifest: [], currentId: null, scene: null, userGene: null, planeIdx: 0,
    drawerOpen: false, bestTab: "pVol", crossMode: "vol",
    crossKey: "pVol", agg: null, pronucleusVisible: {},
    circ: false, aggCirc: null, _aggCircP: null,
  };
  // circularize accessors: when the "blow up the balloon" toggle is on, every read of
  // the current embryo's analysis / transcripts, and the cross-embryo aggregate, uses
  // the precomputed circularized (spherical, seg-1-only) version instead of the real one.
  const curA = () => (state.circ && state.scene && state.scene.circ) ? state.scene.circ.analysis : (state.scene && state.scene.analysis);
  const curTX = () => (state.circ && state.scene && state.scene.circ) ? state.scene.circ.transcripts : (state.scene && state.scene.transcripts);
  const curAGG = () => (state.circ && state.aggCirc) ? state.aggCirc : state.agg;
  const BESTKEY_LABEL = { pVol: "min p · vol", pCnt: "min p · count", diffVol: "max Δ · vol", diffCnt: "max Δ · count" };

  // ---------- boot ----------
  (async function init() {
    try {
      const m = await (await fetch("data/zygote_manifest.json")).json();
      state.manifest = m.embryos;
      state.nPlanes = m.n_planes; state.step = m.step_deg;
      countEl.textContent = `${m.embryos.length} zygotes`;
      V.buildTabs(tabsEl, m.embryos, selectEmbryo, (e) => ({
        label: e.label, sub: e.date_short,
        title: `${e.label} · ${e.n_transcripts.toLocaleString()} transcripts`,
      }));
      // plane selector 0..16
      planeSelect.innerHTML = "";
      for (let k = 0; k < m.n_planes; k++)
        planeSelect.innerHTML += `<option value="${k}">${k * m.step_deg}°</option>`;
      V.wireWindow(controlsEl, $("#controls-header"),
        [...controlsEl.querySelectorAll(".rz")], "zygote_controls_box")
        .setResizeCb(() => { try { Plotly.Plots.resize(chartEl); } catch (_) {} });
      wireDrawer(); wireRdrawer(); wireControls();
    } catch (err) { showError("Failed to load manifest: " + (err.message || err)); }
  })();

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    V.markActiveTab(tabsEl, id);
    showLoading(`Loading ${state.byLabel(id)}…`);
    try {
      const scene = await V.loadGz(`data/zygote/${id}.json.gz`);
      if (state.currentId !== id) return;
      state.scene = scene;
      populateGenes(scene);
      controlsEl.hidden = false; placeholder.hidden = true;
      drawer.hidden = false; rdrawer.hidden = false;
      render(); renderChart(); renderBestList();
      if (state.drawerOpen && state.agg) {
        renderCrossAgg();
      }

    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }
  state.byLabel = (id) => (state.manifest.find((m) => m.id === id) || {}).label || id;

  function populateGenes(scene) {
    const tot = scene.gene_totals || {};
    const sorted = [...scene.genes].sort((a, b) =>
      a.localeCompare(b, undefined, { sensitivity: "base", numeric: true }));
    geneSelect.innerHTML = sorted
      .map((g) => `<option value="${g}">${g}  (${(tot[g] || 0).toLocaleString()})</option>`).join("");
    geneSelect.value = (state.userGene && scene.genes.includes(state.userGene)) ? state.userGene : scene.genes[0];
  }

  // ---------- geometry helpers ----------
  const gene = () => geneSelect.value;
  function planeGeo(k) { return curA().planes[k]; }
  // Split a gene's transcripts into three clouds. Only segment-1 transcripts
  // (t.s1[i]) are counted; those get side A (blue) / side B (red) by plane side:
  // (pos_um − com_um)·normal_um  (matches precompute). Non-segment-1 transcripts
  // are not counted and render green.
  function splitCloud(scene, g, k) {
    const t = curTX()[g]; const A = curA();
    const com = A.com_um, nrm = planeGeo(k).normal_um, zs = scene.z_scale;
    const s1 = t && t.s1;
    const bx = [], by = [], bz = [], rx = [], ry = [], rz = [], gx = [], gy = [], gz = [];
    if (t) for (let i = 0; i < t.x.length; i++) {
      if (s1 && !s1[i]) { gx.push(t.x[i]); gy.push(t.y[i]); gz.push(t.gz[i] * zs); continue; }
      const s = (t.x[i] * XY - com[0]) * nrm[0] + (t.y[i] * XY - com[1]) * nrm[1] + (t.gz[i] - com[2]) * nrm[2];
      if (s > 0) { bx.push(t.x[i]); by.push(t.y[i]); bz.push(t.gz[i] * zs); }
      else { rx.push(t.x[i]); ry.push(t.y[i]); rz.push(t.gz[i] * zs); }
    }
    return { bx, by, bz, rx, ry, rz, gx, gy, gz };
  }
  // plane k as an orange quad (physical square in µm → plot space)
  const PLANE_SCALE = 2;                              // rendered plane size (×) vs precomputed L
  function planeQuad(k, color, op, name, rank, showLegend = true) {
    const A = curA(), p = planeGeo(k), zs = state.scene.z_scale;
    const com = A.com_um, L = p.L * PLANE_SCALE;
    const aUm = [p.a_plot[0] * XY, p.a_plot[1] * XY, p.a_plot[2] / zs];
    const mUm = [p.m_plot[0] * XY, p.m_plot[1] * XY, p.m_plot[2] / zs];
    const cor = [[-1, -1], [-1, 1], [1, 1], [1, -1]].map(([sa, sm]) => {
      const um = [com[0] + sa * L * aUm[0] + sm * L * mUm[0],
                  com[1] + sa * L * aUm[1] + sm * L * mUm[1],
                  com[2] + sa * L * aUm[2] + sm * L * mUm[2]];
      return [um[0] / XY, um[1] / XY, um[2] * zs];       // → plot
    });
    return {
      type: "mesh3d", x: cor.map((c) => c[0]), y: cor.map((c) => c[1]), z: cor.map((c) => c[2]),
      i: [0, 0], j: [1, 2], k: [2, 3], color, opacity: op, name, showlegend: showLegend,
      hoverinfo: "name", flatshading: true, legendrank: rank,
    };
  }

  // ---------- 3-D render ----------
  function render() {
    const s = state.scene; if (!s) return;
    const zs = s.z_scale, A = curA(), k = state.planeIdx, g = gene();
    // when circularized, render the balloon-inflated segment-1 mesh in place of the real one
    const circOn = state.circ && s.circ && s.circ.mesh1;
    const bodyScene = circOn
      ? Object.assign({}, s, { region_meshes: Object.assign({}, s.region_meshes, { "1": s.circ.mesh1 }) })
      : s;
    const traces = V.bodyTraces(bodyScene);

    const sp = splitCloud(s, g, k);
    traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side A`,
      x: sp.bx, y: sp.by, z: sp.bz, marker: { size: 0.5, color: BLUE, opacity: 0.85, line: { width: 0 } },
      hovertemplate: `${g} · side A (counted)<extra></extra>`, legendrank: 20000 });
    traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side B`,
      x: sp.rx, y: sp.ry, z: sp.rz, marker: { size: 0.5, color: RED, opacity: 0.85, line: { width: 0 } },
      hovertemplate: `${g} · side B (counted)<extra></extra>`, legendrank: 20001 });
    traces.push({ type: "scatter3d", mode: "markers", name: `${g} · not counted`,
      x: sp.gx, y: sp.gy, z: sp.gz, marker: { size: 0.5, color: GREEN, opacity: 0.7, line: { width: 0 } },
      hovertemplate: `${g} · not counted (not segment 1)<extra></extra>`, legendrank: 20002 });

    if (axisShow.checked) {
      const c = A.com_plot, ax = A.axis_plot;
      const an = Math.hypot(ax[0], ax[1], ax[2]) || 1;
      const ex = s.extents, R = 0.62 * Math.max(ex.x[1] - ex.x[0], ex.y[1] - ex.y[0], ex.z[1] - ex.z[0]);
      const u = [ax[0] / an, ax[1] / an, ax[2] / an];
      traces.push({ type: "scatter3d", mode: "lines", name: "Polar-body axis",
        x: [c[0] - R * u[0], c[0] + R * u[0]],
        y: [c[1] - R * u[1], c[1] + R * u[1]],
        z: [c[2] - R * u[2], c[2] + R * u[2]],
        line: { color: AXIS_C, width: 6 }, hovertemplate: "Polar-body axis<extra></extra>", legendrank: 40000 });
      traces.push({ type: "scatter3d", mode: "markers", name: "Polar body",
        x: [A.pb_plot[0]], y: [A.pb_plot[1]], z: [A.pb_plot[2]],
        marker: { size: 7, color: AXIS_C, symbol: "circle", line: { width: 1, color: "#fff" } },
        hovertemplate: "Polar body<extra></extra>", legendrank: 40001 });
    }
    if (planeShow.checked)
      traces.push(planeQuad(k, PLANE_C, 0.28, `Plane ${k * state.step}°`, 41000));
    if (allShow.checked) {
      // all planes, colored by weighted p (same viridis code as the bottom drawer:
      // low p → dark/blue, high p → yellow). Uses the cross-section's vol/count mode.
      const planes = A.planes;
      const pOf = (kk) => state.crossMode === "vol" ? planes[kk].wpVol : planes[kk].wpCnt;
      const ps = planes.map((_, kk) => pOf(kk));
      const pmin = Math.min(...ps), pmax = Math.max(...ps), span = pmax - pmin || 1;
      let rank = 42000;
      for (let kk = 0; kk < planes.length; kk++) {
        if (planeShow.checked && kk === k) continue;   // selected plane already drawn in orange
        const t = (pmax - pOf(kk)) / span;             // 1 = lowest p = most intense
        traces.push(planeQuad(kk, viridis(1 - t), 0.16,
          `Plane ${kk * state.step}° · p=${fmtP(pOf(kk))}`, rank++, false));
      }
    }
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  // ---------- counts chart (floating window) ----------
  function geneRow(g, k) {
    const rows = curA().genes;
    const r = rows.find((x) => x.gene === g);
    return r ? r.planes[k] : null;
  }
  function renderChart() {
    const s = state.scene; if (!s) return;
    const g = gene(), k = state.planeIdx, row = geneRow(g, k);
    chartSub.textContent = `· ${g} · plane ${k * state.step}°`;
    if (!row) { Plotly.purge(chartEl); chartEl.classList.remove("js-plotly-plot"); chartEl.innerHTML = '<div class="chart-readout">Gene not in this embryo.</div>'; chartReadout.innerHTML = ""; return; }
    const traces = [
      { type: "bar", name: "Real", x: ["Side A", "Side B"], y: [row.a, row.b],
        marker: { color: [BLUE, RED] }, hovertemplate: "%{x}: %{y}<extra>real</extra>" },
      { type: "bar", name: "Null", x: ["Side A", "Side B"], y: [row.na, row.nb],
        marker: { color: "#9ca3af" }, opacity: 0.7, hovertemplate: "%{x}: %{y}<extra>null</extra>" },
    ];
    plotInto(chartEl, traces, {
      barmode: "group", margin: { l: 34, r: 6, t: 6, b: 20 }, height: 150,
      yaxis: { tickfont: { size: 9 }, gridcolor: "#eef1f5", fixedrange: true, title: { text: "count", font: { size: 9 } } },
      xaxis: { tickfont: { size: 11 }, fixedrange: true }, bargap: 0.3, bargroupgap: 0.08,
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", font: { color: "#1a2233" },
      legend: { orientation: "h", font: { size: 10 }, y: 1.15, x: 1, xanchor: "right" },
    });
    const sigV = row.pVol <= 0.05, sigC = row.pCnt <= 0.05;
    chartReadout.innerHTML =
      `<div>n = <b>${row.a + row.b}</b> · Δcount = <b>${row.dCount}</b> ` +
      `(real) vs <b>${row.ndCount}</b> (null)</div>` +
      `<div>Δ / total = <b>${row.dNorm.toFixed(3)}</b> · Δ / volume = <b>${row.dVol.toExponential(2)}</b></div>` +
      `<div>p(vol) = <span class="${sigV ? "sig" : ""}">${fmtP(row.pVol)}</span> · ` +
      `p(count) = <span class="${sigC ? "sig" : ""}">${fmtP(row.pCnt)}</span></div>`;
  }
  const fmtP = (p) => (p == null || !isFinite(p)) ? "n/a" : p < 0.001 ? p.toExponential(1) : p.toFixed(3);

  // ---------- best-plane list (right drawer) ----------
  // Four best planes: {min weighted p, max Σ|diff|} × {volume, count}.
  const BEST_META = {
    pVol:    { plane: "pVol",    isVol: true,  which: "minimizes the transcript-weighted mean p (volume)" },
    pCnt:    { plane: "pCnt",    isVol: false, which: "minimizes the transcript-weighted mean p (count)" },
    diffVol: { plane: "diffVol", isVol: true,  which: "maximizes Σ|side density-difference| / total (volume)" },
    diffCnt: { plane: "diffCnt", isVol: false, which: "maximizes Σ|side count-difference| / total (count)" },
  };
  function renderBestList() {
    const s = state.scene; if (!s) return;
    const A = curA(), meta = BEST_META[state.bestTab];
    const k = A.best_planes[meta.plane];
    const isVol = meta.isVol;
    // real vs null side-difference under the tab's normalization; p under the same
    const real = (r) => isVol ? r.dVol : r.dNorm;
    const nul = (r) => isVol ? r.ndVol : r.ndNorm;
    const pv = (r) => isVol ? r.pVol : r.pCnt;
    const rows = A.genes.map((r) => ({ gene: r.gene, ...r.planes[k], total: r.total }))
      .sort((a, b) => pv(a) - pv(b) || Math.abs(real(b)) - Math.abs(real(a)));
    const cur = gene();
    const fmtD = (x) => isVol ? x.toExponential(1) : x.toFixed(3);
    let html = `<div class="best-plane-note">Plane <b>${k * state.step}°</b> — ${meta.which}.</div>`;
    html += `<div class="best-head"><span></span><span>gene</span><span>Δ real</span><span>Δ null</span>` +
            `<span>p(${isVol ? "vol" : "cnt"})</span></div>`;
    html += rows.map((r, i) =>
      `<div class="best-row${r.gene === cur ? " current" : ""}" data-gene="${r.gene}"` +
      ` title="n=${r.total} · Δcount ${r.dCount} · p(vol) ${fmtP(r.pVol)} · p(count) ${fmtP(r.pCnt)}">` +
      `<span class="best-num">${i + 1}</span>` +
      `<span class="best-gene">${r.gene}</span>` +
      `<span class="best-real">${fmtD(real(r))}</span>` +
      `<span class="best-null">${fmtD(nul(r))}</span>` +
      `<span class="best-p${pv(r) <= 0.05 ? " sig" : ""}">${fmtP(pv(r))}</span></div>`).join("");
    bestListEl.innerHTML = html;
  }

  // ---------- cross-section plot (bottom drawer) ----------
  // Viridis colormap (low p → bright/yellow end = more intense).
  const VIRIDIS = [[68, 1, 84], [71, 44, 122], [59, 81, 139], [44, 113, 142], [33, 144, 141],
                   [39, 173, 129], [92, 200, 99], [170, 220, 50], [253, 231, 37]];
  function viridis(t) {
    t = Math.max(0, Math.min(1, t));
    const x = t * (VIRIDIS.length - 1), i = Math.min(VIRIDIS.length - 2, Math.floor(x)), f = x - i;
    const a = VIRIDIS[i], b = VIRIDIS[i + 1];
    return `rgb(${Math.round(a[0] + (b[0] - a[0]) * f)},${Math.round(a[1] + (b[1] - a[1]) * f)},` +
           `${Math.round(a[2] + (b[2] - a[2]) * f)})`;
  }
  // ---------- cross-embryo bottom drawer (all zygotes) ----------
  // Loads the aggregate (every embryo's aligned cross-section outline + per-gene
  // side counts at the 4 best planes) once, lazily.
  function ensureAgg() {
    if (state.agg) return Promise.resolve(state.agg);
    if (state._aggP) return state._aggP;
    if (!state._aggCircP) state._aggCircP = V.loadGz("data/zygote_cross_circ.json.gz")
      .then((a) => { state.aggCirc = a; }).catch(() => {});
    state._aggP = V.loadGz("data/zygote_cross.json.gz").then((agg) => {
      state.agg = agg;
      return agg;
    });
    return state._aggP;
  }
  const bestKeyIndex = () => curAGG().best_keys.indexOf(state.crossKey);
  const XS_CFG = {
    responsive: true, scrollZoom: true, displaylogo: false, displayModeBar: "hover",
    modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "toggleSpikelines",
      "hoverClosestCartesian", "hoverCompareCartesian", "zoom2d"],
  };
  const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const unit3 = (v) => { const d = Math.hypot(v[0], v[1], v[2]) || 1; return v.map((x) => x / d); };
  const hexAlpha = (hex, alpha) => {
    const h = hex.replace("#", "");
    const n = parseInt(h.length === 3 ? h.split("").map((c) => c + c).join("") : h, 16);
    return `rgba(${n >> 16},${(n >> 8) & 255},${n & 255},${alpha})`;
  };
  function sectionBasis() {
    const p = planeGeo(state.planeIdx), s = state.scene;
    return {
      com: curA().com_um,
      normal: unit3(p.normal_um),
      axis: unit3([p.a_plot[0] * XY, p.a_plot[1] * XY, p.a_plot[2] / s.z_scale]),
      depth: unit3([p.m_plot[0] * XY, p.m_plot[1] * XY, p.m_plot[2] / s.z_scale]),
    };
  }
  function projectedMeshSection(mesh, basis) {
    if (!mesh || !mesh.verts || !mesh.faces) return { x: [], y: [] };
    const v = mesh.verts, f = mesh.faces, zs = state.scene.z_scale, pts = new Array(v.length / 3);
    for (let i = 0; i < pts.length; i++) {
      const rel = [v[i * 3] * XY - basis.com[0], v[i * 3 + 1] * XY - basis.com[1], v[i * 3 + 2] / zs - basis.com[2]];
      pts[i] = { x: dot3(rel, basis.normal), y: dot3(rel, basis.axis), d: dot3(rel, basis.depth) };
    }
    const segments = [], eps = 1e-7;
    for (let i = 0; i < f.length; i += 3) {
      const tri = [pts[f[i]], pts[f[i + 1]], pts[f[i + 2]]], hits = [];
      for (let e = 0; e < 3; e++) {
        const a = tri[e], b = tri[(e + 1) % 3];
        if (Math.abs(a.d) < eps) hits.push([a.x, a.y]);
        if (a.d * b.d < 0) {
          const t = a.d / (a.d - b.d);
          hits.push([a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)]);
        }
      }
      const unique = hits.filter((p, j) => !hits.slice(0, j).some((q) => Math.hypot(p[0] - q[0], p[1] - q[1]) < 1e-5));
      if (unique.length < 2) continue;
      let a = unique[0], b = unique[1], far = -1;
      for (let q = 0; q < unique.length; q++) for (let r = q + 1; r < unique.length; r++) {
        const d = Math.hypot(unique[q][0] - unique[r][0], unique[q][1] - unique[r][1]);
        if (d > far) { far = d; a = unique[q]; b = unique[r]; }
      }
      segments.push([a, b]);
    }
    // Adjacent mesh triangles meet at the same interpolated points. Quantize those
    // endpoints and walk the resulting graph to recover closed cross-section loops.
    const nodes = new Map(), keyOf = (p) => `${p[0].toFixed(3)},${p[1].toFixed(3)}`;
    const nodeFor = (p) => {
      const key = keyOf(p);
      if (!nodes.has(key)) nodes.set(key, { key, x: p[0], y: p[1], links: new Set() });
      return nodes.get(key);
    };
    for (const [pa, pb] of segments) {
      const a = nodeFor(pa), b = nodeFor(pb); if (a.key === b.key) continue;
      a.links.add(b.key); b.links.add(a.key);
    }
    const edgeKey = (a, b) => a < b ? `${a}|${b}` : `${b}|${a}`, used = new Set(), paths = [];
    for (const start of nodes.values()) for (const first of start.links) {
      if (used.has(edgeKey(start.key, first))) continue;
      const path = [[start.x, start.y]]; let current = start.key, next = first;
      for (let guard = 0; guard < nodes.size + 2; guard++) {
        used.add(edgeKey(current, next));
        const n = nodes.get(next); path.push([n.x, n.y]);
        if (next === start.key) break;
        const options = [...n.links].filter((k) => !used.has(edgeKey(next, k)));
        if (!options.length) break;
        current = next; next = options[0];
      }
      if (path.length < 4) continue;
      const d = Math.hypot(path[0][0] - path[path.length - 1][0], path[0][1] - path[path.length - 1][1]);
      if (d > 0.01) path.push(path[0]);
      let area = 0;
      for (let i = 0; i < path.length - 1; i++) area += path[i][0] * path[i + 1][1] - path[i + 1][0] * path[i][1];
      if (Math.abs(area) > 0.02) paths.push(path);
    }
    paths.sort((a, b) => b.length - a.length);
    const xs = [], ys = [];
    for (const path of paths) {
      for (const p of path) { xs.push(p[0]); ys.push(p[1]); }
      xs.push(null); ys.push(null);
    }
    return { x: xs, y: ys };
  }
  function pronucleusLabels() {
    const A = curA(), pb = Number(A.polar_body_label);
    const candidates = (A.polar_body_detection && A.polar_body_detection.candidates) || [];
    const internal = candidates.filter((c) => !c.external && Number(c.label) !== pb).map((c) => Number(c.label));
    return internal.length ? internal : state.scene.mask_labels.filter((x) => x !== 1 && x !== pb);
  }
  function syncPronucleusControls(labels) {
    const key = `${state.currentId}:${labels.join(",")}`;
    if (xsPronuclei.dataset.key === key) return;
    xsPronuclei.dataset.key = key;
    xsPronuclei.innerHTML = labels.map((label, i) =>
      `<label class="xs-tg"><input type="checkbox" data-pronucleus="${label}" checked><span>pronucleus ${i + 1}</span></label>`).join("");
    state.pronucleusVisible = Object.fromEntries(labels.map((label) => [label, true]));
    xsPronuclei.querySelectorAll("input").forEach((input) => input.addEventListener("change", () => {
      state.pronucleusVisible[input.dataset.pronucleus] = input.checked; renderCurrentCrossSection();
    }));
  }
  function renderCurrentCrossSection() {
    if (!state.scene) return;
    const s = state.scene, A = curA(), g = gene(), basis = sectionBasis();
    const pb = Number(A.polar_body_label), pns = pronucleusLabels(); syncPronucleusControls(pns);
    const traces = [], limits = [];
    const addOutline = (label, name, color, width, opacity, visible) => {
      if (!visible) return;
      const mesh = (label === 1 && state.circ && s.circ && s.circ.mesh1) ? s.circ.mesh1 : s.region_meshes[String(label)];
      const sec = projectedMeshSection(mesh, basis); if (!sec.x.length) return;
      sec.x.forEach((x, i) => { if (x != null) limits.push(Math.abs(x), Math.abs(sec.y[i])); });
      traces.push({ type: "scatter", mode: "lines", x: sec.x, y: sec.y, name,
        fill: "toself", fillcolor: hexAlpha(xsCrossFillColor.value, opacity),
        line: { color, width }, hoverinfo: "name", showlegend: xsCrossLegend.checked });
    };
    addOutline(1, "Cell outline", xsCrossBodyColor.value, 1.8, 0.14, xsBody.checked);
    addOutline(pb, "Polar body", xsCrossPbColor.value, 1.5, 0.34, xsPb.checked);
    pns.forEach((label, i) => addOutline(label, `Pronucleus ${i + 1}`, xsCrossPnColor.value, 1.35, 0.30, state.pronucleusVisible[label] !== false));
    const t = curTX()[g], sideA = { x: [], y: [] }, sideB = { x: [], y: [] };
    if (t) for (let i = 0; i < t.x.length; i++) {
      if (t.s1 && !t.s1[i]) continue;
      const rel = [t.x[i] * XY - basis.com[0], t.y[i] * XY - basis.com[1], t.gz[i] - basis.com[2]];
      const x = dot3(rel, basis.normal), y = dot3(rel, basis.axis), dst = x > 0 ? sideA : sideB;
      dst.x.push(x); dst.y.push(y); limits.push(Math.abs(x), Math.abs(y));
    }
    traces.push(
      { type: "scattergl", mode: "markers", x: sideA.x, y: sideA.y, name: "Side A transcripts",
        marker: { color: xsCrossHighColor.value, size: 3, opacity: 0.58 }, hoverinfo: "skip", showlegend: xsCrossLegend.checked },
      { type: "scattergl", mode: "markers", x: sideB.x, y: sideB.y, name: "Side B transcripts",
        marker: { color: xsCrossLowColor.value, size: 3, opacity: 0.58 }, hoverinfo: "skip", showlegend: xsCrossLegend.checked }
    );
    const pbUm = [A.pb_plot[0] * XY, A.pb_plot[1] * XY, A.pb_plot[2] / s.z_scale];
    const pbRel = [pbUm[0] - basis.com[0], pbUm[1] - basis.com[1], pbUm[2] - basis.com[2]];
    const pbCenter = [dot3(pbRel, basis.normal), dot3(pbRel, basis.axis)];
    limits.push(Math.abs(pbCenter[0]), Math.abs(pbCenter[1]));
    traces.push(
      { type: "scatter", mode: "markers", x: [0], y: [0], name: "Cell-body center of mass",
        marker: { symbol: "circle", size: 8, color: "#111827", line: { color: "#ffffff", width: 1.2 } },
        hovertemplate: "Cell-body center of mass<extra></extra>", showlegend: xsCrossLegend.checked },
      { type: "scatter", mode: "markers", x: [pbCenter[0]], y: [pbCenter[1]], name: "Polar-body center of mass",
        marker: { symbol: "diamond", size: 9, color: "#1e3a5f", line: { color: "#ffffff", width: 1.2 } },
        hovertemplate: "Polar-body center of mass<extra></extra>", showlegend: xsCrossLegend.checked }
    );
    const lim = Math.max(20, ...(limits.length ? limits : [20])) * 1.08;
    traces.push({ type: "scatter", mode: "lines", x: [0, 0], y: [-lim, lim], name: "Division plane",
      line: { color: "#111827", width: 1.5, dash: "dash" }, hoverinfo: "skip", showlegend: xsCrossLegend.checked });
    xsCaption.textContent = `· ${s.id} · ${g} · ${state.planeIdx * state.step}°`;
    plotInto(xsOutlines, traces, {
      dragmode: "pan", margin: { l: 42, r: 10, t: 8, b: 38 }, height: xsOutlines.clientHeight || 330,
      xaxis: { title: { text: "Distance from division plane (µm)", font: { size: 10 } }, range: [-lim, lim],
        scaleanchor: "y", scaleratio: 1, gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
      yaxis: { title: { text: "Polar-body axis (µm)", font: { size: 10 } }, range: [-lim, lim],
        gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
      legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", uirevision: `${s.id}-${state.planeIdx}`,
    }, XS_CFG);
  }
  function binomial95(n) {
    if (!n) return [0, 0];
    const mode = Math.floor(n / 2), weights = new Float64Array(n + 1); weights[mode] = 1;
    for (let k = mode - 1; k >= 0; k--) weights[k] = weights[k + 1] * (k + 1) / (n - k);
    for (let k = mode + 1; k <= n; k++) weights[k] = weights[k - 1] * (n - k + 1) / k;
    let sum = 0; for (const w of weights) sum += w;
    let cdf = 0, lo = 0, hi = n;
    for (let k = 0; k <= n; k++) { cdf += weights[k] / sum; if (cdf >= 0.025) { lo = k; break; } }
    cdf = 0; for (let k = 0; k <= n; k++) { cdf += weights[k] / sum; if (cdf >= 0.975) { hi = k; break; } }
    return [lo, hi];
  }
  function renderBars() {
    const agg = curAGG(); if (!agg) return;
    const ki = bestKeyIndex(), g = gene();
    const rows = agg.embryos.map((emb) => {
      const row = emb.g[g]; if (!row) return null;
      const a = row[1 + ki], b = row[0] - a, total = a + b; if (!total) return null;
      const [nullLow, nullHigh] = binomial95(total);
      return { label: emb.label, high: Math.max(a, b), low: Math.min(a, b), total,
        nullMean: total / 2, nullLow, nullHigh };
    }).filter(Boolean);
    if (!rows.length) {
      Plotly.purge(xsBars); xsBars.classList.remove("js-plotly-plot");   // so plotInto re-inits cleanly next time
      xsBars.innerHTML = `<div class="xs-empty"><div><b>${g}</b> is not detected in any retained zygote.</div></div>`;
      xsBarSub.textContent = `· ${g}`;
      return;
    }
    const x = rows.map((_, i) => i), baseline = xsBarLog.checked ? 1 : 0, shapes = [], traces = [];
    rows.forEach((r, i) => {
      if (xsBarLog.checked) {
        // On a log axis both observed halves share the same baseline and sit side by side.
        shapes.push(
          { type: "rect", xref: "x", yref: "y", layer: "above", x0: i - 0.24, x1: i - 0.03,
            y0: baseline, y1: r.high, fillcolor: xsBarHighColor.value, line: { color: "rgba(15,23,42,0.48)", width: 0.65 } },
          { type: "rect", xref: "x", yref: "y", layer: "above", x0: i + 0.03, x1: i + 0.24,
            y0: baseline, y1: r.low, fillcolor: xsBarLowColor.value, line: { color: "rgba(15,23,42,0.48)", width: 0.65 } }
        );
      } else {
        // Linear mode preserves the cumulative stacked comparison.
        shapes.push(
          { type: "rect", xref: "x", yref: "y", layer: "above", x0: i - 0.23, x1: i + 0.23,
            y0: baseline, y1: r.high, fillcolor: xsBarHighColor.value, line: { color: "rgba(15,23,42,0.48)", width: 0.65 } },
          { type: "rect", xref: "x", yref: "y", layer: "above", x0: i - 0.23, x1: i + 0.23,
            y0: r.high, y1: r.total, fillcolor: xsBarLowColor.value, line: { color: "rgba(15,23,42,0.48)", width: 0.65 } }
        );
      }
      if (xsBarNull.checked) shapes.push({ type: "rect", xref: "x", yref: "y", layer: "above",
        x0: i - 0.34, x1: i + 0.34, y0: baseline, y1: r.nullMean,
        fillcolor: hexAlpha(xsBarNullColor.value, 0.30), line: { color: "#5f666d", width: 0.7 } });
    });
    if (xsBarLegend.checked) {
      traces.push(
        { type: "scatter", mode: "markers", name: "Higher-count half", legendrank: 10, x: [null], y: [null],
          marker: { symbol: "square", size: 11, color: xsBarHighColor.value }, hoverinfo: "skip" },
        { type: "scatter", mode: "markers", name: "Lower-count half", legendrank: 20, x: [null], y: [null],
          marker: { symbol: "square", size: 11, color: xsBarLowColor.value }, hoverinfo: "skip" }
      );
      if (xsBarNull.checked) traces.push({ type: "scatter", mode: "markers", name: "Null mean", legendrank: 30,
        x: [null], y: [null], marker: { symbol: "square", size: 11, color: hexAlpha(xsBarNullColor.value, 0.6) }, hoverinfo: "skip" });
    }
    if (xsBarInterval.checked) traces.push({ type: "scatter", mode: "markers", name: "95% null interval", x,
      legendrank: 40,
      y: rows.map((r) => r.nullMean), marker: { size: 1, color: "rgba(0,0,0,0)" },
      error_y: { type: "data", symmetric: false, array: rows.map((r) => r.nullHigh - r.nullMean),
        arrayminus: rows.map((r) => r.nullMean - r.nullLow), color: "#111827", thickness: 1.5, width: 4 },
      hovertemplate: "%{x}<br>95% null interval<extra></extra>" });
    traces.push({ type: "scatter", mode: "markers", showlegend: false, x,
      y: rows.map((r) => xsBarLog.checked ? Math.max(r.high, r.low) : r.total),
      customdata: rows.map((r) => [r.label, r.high, r.low, r.total, r.nullMean, r.nullLow, r.nullHigh]),
      marker: { size: 24, color: "rgba(0,0,0,0)" },
      hovertemplate: "%{customdata[0]}<br>higher half %{customdata[1]}<br>lower half %{customdata[2]}" +
        "<br>total %{customdata[3]}<br>null %{customdata[4]:.1f} (%{customdata[5]}–%{customdata[6]})<extra></extra>" });
    xsBarSub.textContent = `· ${g} · ${rows.length} zygotes`;
    const maxY = Math.max(...rows.map((r) => xsBarLog.checked
      ? Math.max(r.high, r.low, r.nullHigh)
      : Math.max(r.total, r.nullHigh))) * 1.18;
    plotInto(xsBars, traces, {
      shapes, margin: { l: 56, r: 10, t: xsBarLegend.checked ? 48 : 8, b: 92 }, height: xsBars.clientHeight || 330,
      showlegend: xsBarLegend.checked,
      legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      xaxis: { fixedrange: false, tickmode: "array", tickvals: x, ticktext: rows.map((r) => r.label),
        range: [-0.55, rows.length - 0.45], tickangle: -48, tickfont: { size: 8 }, automargin: true },
      yaxis: { title: { text: `${g} transcript count`, font: { size: 10 } }, tickfont: { size: 9 },
               type: xsBarLog.checked ? "log" : "linear", range: xsBarLog.checked ? [0, Math.log10(maxY)] : [0, maxY],
               dtick: xsBarLog.checked ? 1 : undefined,
               gridcolor: xsBarGrid.checked ? "#e2e5e8" : "rgba(0,0,0,0)", fixedrange: false },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, XS_CFG);
  }
  function renderCrossAgg() {
    return ensureAgg().then(() => {
      const cov = curAGG().gene_cov[gene()] || 0, tot = curAGG().n_embryos;
      drawerEmb.textContent = `· ${gene()}`;
      xsNote.innerHTML = `Counts use each zygote's <b>${BESTKEY_LABEL[state.crossKey]}</b> plane; <b>${cov}/${tot}</b> retained zygotes contain <b>${gene()}</b>.`;
      renderCurrentCrossSection(); renderBars();
      requestAnimationFrame(() => { try { Plotly.Plots.resize(xsOutlines); Plotly.Plots.resize(xsBars); } catch (_) {} });
    });
  }

  // ---------- wiring ----------
  function wireControls() {
    geneSelect.addEventListener("change", () => {
      state.userGene = geneSelect.value; render(); renderChart(); highlightBest();
      if (state.drawerOpen && state.agg) renderCrossAgg();
    });
    planeSelect.addEventListener("change", () => {
      state.planeIdx = parseInt(planeSelect.value, 10) || 0; render(); renderChart();
      if (state.drawerOpen) renderCurrentCrossSection();
    });
    [axisShow, planeShow, allShow].forEach((c) => c.addEventListener("change", () => render()));
    // circularize ("blow up the balloon"): switch every plot + analysis to the precomputed
    // spherical (seg-1) version and re-render the 3-D view, chart, best-planes, and drawer.
    circShow.addEventListener("change", () => {
      state.circ = circShow.checked;
      render(); renderChart(); renderBestList();
      if (state.drawerOpen && state.agg) renderCrossAgg();
    });
    // p-value colormap normalization (drives the "All planes" 3-D coloring)
    pcolorMode.addEventListener("change", () => { state.crossMode = pcolorMode.value; if (allShow.checked) render(); });
  }
  function highlightBest() {
    const cur = gene();
    bestListEl.querySelectorAll(".best-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === cur));
  }
  // Make a drawer's main handle pull it open and resize it by dragging, while a
  // plain tap still toggles it. cfg.computeSize/clampSize/applySize map the pointer
  // to the drawer's size CSS var; cfg.setOpen(open) shows/hides it.
  function wireHandleDrag(drawerEl, handleEl, cfg) {
    let start = null, moved = false;
    handleEl.addEventListener("pointerdown", (e) => {
      if (e.button && e.button !== 0) return;
      start = { x: e.clientX, y: e.clientY }; moved = false;
      try { handleEl.setPointerCapture(e.pointerId); } catch (_) {}
    });
    handleEl.addEventListener("pointermove", (e) => {
      if (!start) return;
      if (!moved && Math.hypot(e.clientX - start.x, e.clientY - start.y) < 5) return;
      if (!moved) { moved = true; drawerEl.classList.add("dragging"); if (drawerEl.dataset.open !== "true") cfg.setOpen(true); }
      cfg.applySize(cfg.clampSize(cfg.computeSize(e)));
      e.preventDefault();
    });
    const up = (e) => {
      if (!start) return;
      try { handleEl.releasePointerCapture(e.pointerId); } catch (_) {}
      if (moved) { drawerEl.classList.remove("dragging"); cfg.afterDrag && cfg.afterDrag(); }
      else { cfg.setOpen(drawerEl.dataset.open !== "true"); }   // tap = toggle
      start = null; moved = false;
    };
    handleEl.addEventListener("pointerup", up);
    handleEl.addEventListener("pointercancel", up);
  }
  function wireRdrawer() {
    const setOpen = (open) => {
      rdrawer.dataset.open = open ? "true" : "false";
      rdrawerHandle.setAttribute("aria-expanded", String(open));
    };
    wireHandleDrag(rdrawer, rdrawerHandle, {
      computeSize: (e) => window.innerWidth - e.clientX,
      clampSize: (px) => Math.max(240, Math.min(window.innerWidth - 80, px)),
      applySize: (px) => rdrawer.style.setProperty("--rdrawer-w", px + "px"),
      setOpen,
    });
    // left-edge grabber: fine width adjustment while open
    let sw = 0;
    const rrz = $("#rdrawer-resize");
    rrz.addEventListener("pointerdown", (e) => { sw = rdrawer.getBoundingClientRect().width; rrz._d = { x: e.clientX }; rrz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rrz.addEventListener("pointermove", (e) => { if (!rrz._d) return;
      const w = Math.max(240, Math.min(window.innerWidth - 80, sw - (e.clientX - rrz._d.x)));
      rdrawer.style.setProperty("--rdrawer-w", w + "px"); });
    const rend = (e) => { if (rrz._d) { rrz._d = null; try { rrz.releasePointerCapture(e.pointerId); } catch (_) {} } };
    rrz.addEventListener("pointerup", rend); rrz.addEventListener("pointercancel", rend);
    rtabsEl.querySelectorAll(".rtab").forEach((b) => b.addEventListener("click", () => {
      state.bestTab = b.dataset.best;
      rtabsEl.querySelectorAll(".rtab").forEach((x) => x.classList.toggle("active", x === b));
      renderBestList();
    }));
    bestListEl.addEventListener("click", (e) => {
      const row = e.target.closest(".best-row"); if (!row) return;
      const g = row.dataset.gene;
      if (state.scene.genes.includes(g)) {
        state.userGene = g; geneSelect.value = g; render(); renderChart(); highlightBest();
        if (state.drawerOpen && state.agg) renderCrossAgg();
      }
    });
  }
  const resizeXs = () => { try { Plotly.Plots.resize(xsOutlines); Plotly.Plots.resize(xsBars); } catch (_) {} };
  function wireDrawer() {
    const setOpen = (open) => {
      state.drawerOpen = open;
      drawer.dataset.open = open ? "true" : "false";
      drawerHandle.setAttribute("aria-expanded", String(open));
      if (open) renderCrossAgg();
    };
    wireHandleDrag(drawer, drawerHandle, {
      computeSize: (e) => window.innerHeight - e.clientY - 40,   // body height (handle = 40px)
      clampSize: (px) => Math.max(160, Math.min(window.innerHeight - 100, px)),
      applySize: (px) => drawer.style.setProperty("--drawer-h", px + "px"),
      setOpen,
      afterDrag: resizeXs,
    });
    // top-edge grabber: fine height adjustment while open
    let sh = 0;
    const rz = $("#drawer-resize");
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { x: e.clientX, y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return;
      const h = Math.max(160, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y)));
      drawer.style.setProperty("--drawer-h", h + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} resizeXs(); } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
    // Cross-embryo criterion and publication-plot settings.
    xsPlane.addEventListener("change", () => { state.crossKey = xsPlane.value; if (state.agg) renderCrossAgg(); });
    [xsBody, xsPb, xsCrossLegend, xsCrossHighColor, xsCrossLowColor, xsCrossFillColor, xsCrossBodyColor, xsCrossPbColor, xsCrossPnColor]
      .forEach((el) => el.addEventListener("change", renderCurrentCrossSection));
    [xsBarLog, xsBarLegend, xsBarNull, xsBarInterval, xsBarGrid, xsBarHighColor, xsBarLowColor, xsBarNullColor]
      .forEach((el) => el.addEventListener("change", renderBars));
    xsCrossDownload.addEventListener("click", () => downloadPlot(xsOutlines, "cross_section", xsCrossFormat, xsCrossScale, 1800, 1400));
    xsBarDownload.addEventListener("click", () => downloadPlot(xsBars, "side_counts", xsBarFormat, xsBarScale, 2000, 1250));
  }

  function downloadPlot(div, suffix, formatEl, scaleEl, width, height) {
    if (!state.scene || !div.classList.contains("js-plotly-plot")) return;
    const clean = (x) => String(x || "plot").replace(/[^a-z0-9_-]+/gi, "_").replace(/^_+|_+$/g, "");
    Plotly.downloadImage(div, {
      format: formatEl.value, width, height, scale: Number(scaleEl.value) || 2,
      filename: `${clean(state.scene.id)}_${clean(gene())}_${state.planeIdx * state.step}deg_${suffix}`,
    });
  }

  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false;
    placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
