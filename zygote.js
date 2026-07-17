/* Zygote Division Planes — analysis model.
 * Built on viewer-core.js (VCore): 60 zygotes, the polar-body axis, 18 candidate
 * division planes, per-gene transcript split (blue/red) across the selected plane,
 * a counts+null chart, best-plane tables, and a cross-section plot. All statistics
 * are pre-computed (build_zygote.py); the UI only reads them.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const AXIS_C = "#111827", PLANE_C = "#f97316";
  const BLUE = "#2563eb", RED = "#dc2626", GREEN = "#16a34a";

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
  const xsPlane = $("#xs-plane"), xsAlign = $("#xs-align"), xsNote = $("#xs-note");
  const xsCaption = $("#xs-caption"), xsBarSub = $("#xs-bar-sub");
  const xsOutlines = $("#xs-outlines"), xsBars = $("#xs-bars");
  const xsMean = $("#xs-mean"), xsOnly = $("#xs-only");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle");
  const rtabsEl = $("#rtabs"), bestListEl = $("#best-list");

  const state = {
    manifest: [], currentId: null, scene: null, userGene: null, planeIdx: 0,
    drawerOpen: false, bestTab: "pVol", crossMode: "vol",
    crossKey: "pVol", alignGene: null, agg: null,
    xsShowMean: true, xsOnlyCurrent: false,
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
        // keep the alignment gene in the current embryo's panel so its zygotes can
        // contain the display gene (disjoint panels ⇒ otherwise the bars go empty)
        if (!scene.genes.includes(state.alignGene)) { state.alignGene = pickDefaultAlign(); xsAlign.value = state.alignGene; }
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
      x: sp.bx, y: sp.by, z: sp.bz, marker: { size: 2.6, color: BLUE, opacity: 0.85, line: { width: 0 } },
      hovertemplate: `${g} · side A (counted)<extra></extra>`, legendrank: 20000 });
    traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side B`,
      x: sp.rx, y: sp.ry, z: sp.rz, marker: { size: 2.6, color: RED, opacity: 0.85, line: { width: 0 } },
      hovertemplate: `${g} · side B (counted)<extra></extra>`, legendrank: 20001 });
    traces.push({ type: "scatter3d", mode: "markers", name: `${g} · not counted`,
      x: sp.gx, y: sp.gy, z: sp.gz, marker: { size: 2.6, color: GREEN, opacity: 0.7, line: { width: 0 } },
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
    // lazily load the circularized aggregate too, so the bottom drawer can follow the balloon toggle
    if (!state._aggCircP) state._aggCircP = V.loadGz("data/zygote_cross_circ.json.gz")
      .then((a) => { state.aggCirc = a; }).catch(() => {});
    state._aggP = V.loadGz("data/zygote_cross.json.gz").then((agg) => {
      state.agg = agg;
      // No gene spans all embryos (multiple panels), so the dropdown is the union
      // of genes, alphabetical, each labeled with how many zygotes contain it.
      xsAlign.innerHTML = agg.genes_all
        .map((g) => `<option value="${g}">${g} (${agg.gene_cov[g]}/${agg.n_embryos})</option>`).join("");
      state.alignGene = pickDefaultAlign();
      xsAlign.value = state.alignGene;
      return agg;
    });
    return state._aggP;
  }
  // Default alignment gene = the current embryo's most-widely-covered gene OTHER than
  // the display gene, so the aligned set overlaps whatever display gene is selected
  // (panels are disjoint, so the global widest gene can share no embryo with the
  // current display gene) and the orientation stays independent of the display gene.
  function pickDefaultAlign() {
    const agg = curAGG(), sc = state.scene, skip = gene();
    if (sc && sc.genes) {
      let best = null, bestCov = -1;
      for (const g of sc.genes) {
        if (g === skip) continue;
        const c = agg.gene_cov[g] || 0;
        if (c > bestCov) { bestCov = c; best = g; }
      }
      if (best) return best;
    }
    return agg.default_align_gene;
  }
  const bestKeyIndex = () => curAGG().best_keys.indexOf(state.crossKey);
  // 60 distinct-ish colors via golden-angle hue spread (same index → same embryo
  // in the outline plot and the bar plot).
  const embColor = (i) => `hsl(${((i * 137.508) % 360).toFixed(1)}, 62%, 52%)`;
  // Orientation for one embryo at the chosen best plane: rotate by −θ so the plane
  // is vertical (then +x = plane's side A); flip 180° when the alignment gene's
  // higher-count side is on the left, so it always ends up on the right (+x).
  function embOrient(emb, ki) {
    const theta = emb.best[ki] * curAGG().step_deg * Math.PI / 180;
    const ga = emb.g[state.alignGene];
    let flip = false;
    if (ga) { const a = ga[1 + ki]; flip = (ga[0] - a) > a; }   // sideB (=n−a) > sideA
    return { c: Math.cos(theta), s: Math.sin(theta), flip };
  }
  // interactive config for the cross-section: scroll to zoom, drag to pan, double-click resets.
  // Minimal modebar (only on hover) keeps the interface clean and professional.
  const XS_CFG = {
    responsive: true, scrollZoom: true, displaylogo: false, displayModeBar: "hover",
    modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "toggleSpikelines",
      "hoverClosestCartesian", "hoverCompareCartesian", "zoom2d"],
  };
  // Mean aligned outline: average each embryo's max-radius-per-angle, lightly smoothed —
  // the canonical cross-section the crowd is scattered around.
  function meanOutline(oriented) {
    if (oriented.length < 4) return null;
    const M = 120, sum = new Array(M).fill(0), cnt = new Array(M).fill(0);
    for (const o of oriented) {
      const rb = new Array(M).fill(0);
      for (let k = 0; k < o.xs.length; k++) {
        const r = Math.hypot(o.xs[k], o.ys[k]);
        let bi = Math.floor((Math.atan2(o.ys[k], o.xs[k]) + Math.PI) / (2 * Math.PI) * M) % M;
        if (bi < 0) bi += M;
        if (r > rb[bi]) rb[bi] = r;
      }
      for (let b = 0; b < M; b++) if (rb[b] > 0) { sum[b] += rb[b]; cnt[b]++; }
    }
    const need = oriented.length * 0.4, raw = [];
    for (let b = 0; b < M; b++) raw.push(cnt[b] >= need ? sum[b] / cnt[b] : null);
    const xs = [], ys = [];
    for (let b = 0; b < M; b++) {                      // 3-bin circular smoothing over present bins
      let s = 0, n = 0;
      for (let d = -1; d <= 1; d++) { const v = raw[(b + d + M) % M]; if (v != null) { s += v; n++; } }
      if (!n) continue;
      const r = s / n, th = (b + 0.5) / M * 2 * Math.PI - Math.PI;
      xs.push(r * Math.cos(th)); ys.push(r * Math.sin(th));
    }
    if (xs.length < 12) return null;
    xs.push(xs[0]); ys.push(ys[0]);
    return { xs, ys };
  }
  // significance → viridis position: log10(p) mapped 0.001 → 0 (dark = most significant)
  // through 1 → 1 (yellow = not significant). Matches the "All planes (p-value)" colouring.
  const sigT = (p) => (p == null || !isFinite(p)) ? 1 : Math.max(0, Math.min(1, (Math.log10(Math.max(p, 1e-6)) + 3) / 3));
  const sigOf = (emb) => (emb.sig ? emb.sig[state.crossKey] : null);
  function renderOutlines() {
    const agg = curAGG(); if (!agg) return;
    const ki = bestKeyIndex();
    const oriented = []; let R = 0;
    agg.embryos.forEach((emb) => {
      if (!emb.outline.length || !emb.g[state.alignGene]) return;   // only embryos with the align gene
      const { c, s, flip } = embOrient(emb, ki);
      const xs = [], ys = [];
      for (const p of emb.outline) {
        let x = p[0] * c + p[1] * s, y = -p[0] * s + p[1] * c;   // rotate by −θ
        if (flip) { x = -x; y = -y; }
        xs.push(x); ys.push(y);
        const r = Math.hypot(x, y); if (r > R) R = r;
      }
      xs.push(xs[0]); ys.push(ys[0]);                            // close the loop
      oriented.push({ id: emb.id, label: emb.label, xs, ys, sig: sigOf(emb), isCurrent: emb.id === state.currentId });
    });
    const lim = (R * 1.08) || 1, traces = [], only = state.xsOnlyCurrent;
    // crowd: each embryo coloured by how significantly its transcriptome splits at its best
    // plane (viridis: low p → dark/purple = significant, high p → yellow = not). Log-scaled.
    if (!only) for (const o of oriented) {
      if (o.isCurrent) continue;
      traces.push({ type: "scatter", mode: "lines", x: o.xs, y: o.ys, opacity: 0.6,
        line: { color: viridis(sigT(o.sig)), width: 1.4, shape: "spline", smoothing: 1.0 }, showlegend: false,
        hovertemplate: `${o.label} · p=${fmtP(o.sig)}<extra></extra>` });
    }
    // mean outline (toggleable) — the typical aligned cross-section
    if (!only && state.xsShowMean) {
      const mean = meanOutline(oriented);
      if (mean) traces.push({ type: "scatter", mode: "lines", x: mean.xs, y: mean.ys,
        line: { color: "#0f172a", width: 2.4, shape: "spline", smoothing: 1.0 }, showlegend: false,
        hovertemplate: `mean outline · ${oriented.length} zygotes<extra></extra>` });
    }
    // the currently-viewed embryo — its OWN significance colour, cased in white so it stands out
    // (and is the only thing shown in "only this embryo" mode)
    const cur = oriented.find((o) => o.isCurrent);
    if (cur) {
      traces.push({ type: "scatter", mode: "lines", x: cur.xs, y: cur.ys,
        line: { color: "#ffffff", width: 5.5, shape: "spline", smoothing: 1.0 }, hoverinfo: "skip", showlegend: false });
      traces.push({ type: "scatter", mode: "lines", x: cur.xs, y: cur.ys,
        line: { color: viridis(sigT(cur.sig)), width: 3, shape: "spline", smoothing: 1.0 }, showlegend: false,
        hovertemplate: `${cur.label} · current · p=${fmtP(cur.sig)}<extra></extra>` });
    }
    // division-plane guide (vertical) + centre dot
    traces.push({ type: "scatter", mode: "lines", x: [0, 0], y: [-lim, lim],
      line: { color: "#94a3b8", width: 1.3, dash: "dash" }, hoverinfo: "skip", showlegend: false });
    traces.push({ type: "scatter", mode: "markers", x: [0], y: [0],
      marker: { color: "#94a3b8", size: 4 }, hoverinfo: "skip", showlegend: false });
    // 20 µm scale bar, lower-left
    const sb = 20, bx = -lim * 0.9, by = -lim * 0.9;
    traces.push({ type: "scatter", mode: "lines", x: [bx, bx + sb], y: [by, by],
      line: { color: "#475569", width: 3 }, hoverinfo: "skip", showlegend: false });
    plotInto(xsOutlines, traces, {
      dragmode: "pan", margin: { l: 8, r: 8, t: 8, b: 8 }, height: xsOutlines.clientHeight || 340,
      xaxis: { range: [-lim, lim], visible: false, scaleanchor: "y", scaleratio: 1, constrain: "domain" },
      yaxis: { range: [-lim, lim], visible: false, constrain: "domain" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      annotations: [
        { x: lim, y: 0, xanchor: "right", yanchor: "bottom", showarrow: false,
          text: `higher ${state.alignGene} →`, font: { size: 11, color: "#64748b" } },
        { x: bx + sb / 2, y: by, xanchor: "center", yanchor: "top", yshift: -4, showarrow: false,
          text: "20 µm", font: { size: 10, color: "#475569" } },
      ],
    }, XS_CFG);
  }
  function renderBars() {
    const agg = curAGG(); if (!agg) return;
    const ki = bestKeyIndex(), g = gene();
    const bars = [], connX = [], connY = [];
    let leftCum = 0, rightCum = 0, nEmb = 0;
    agg.embryos.forEach((emb, i) => {
      if (!emb.g[state.alignGene]) return;                      // not orientable (no align gene)
      const row = emb.g[g]; if (!row) return;                   // display gene absent here
      const a = row[1 + ki], b = row[0] - a;
      const { flip } = embOrient(emb, ki);
      const right = flip ? b : a, left = flip ? a : b;
      if (left + right === 0) return;
      nEmb++;
      bars.push({ type: "bar", x: [0, 1], y: [left, right], width: 0.5,
        marker: { color: embColor(i), line: { width: 0 } }, showlegend: false,
        hovertemplate: `${emb.label}<br>left ${left} · right ${right}<extra></extra>` });
      leftCum += left; rightCum += right;
      connX.push(0.25, 0.75, null); connY.push(leftCum, rightCum, null);
    });
    if (!nEmb) {                     // display gene shares no zygote with the align gene
      Plotly.purge(xsBars); xsBars.classList.remove("js-plotly-plot");   // so plotInto re-inits cleanly next time
      xsBars.innerHTML = `<div class="xs-empty"><div><b>${g}</b> is not detected in any of the ` +
        `${agg.gene_cov[state.alignGene] || 0} zygotes aligned by <b>${state.alignGene}</b> ` +
        `(different MERFISH panels).<br>Pick a gene shared with those zygotes, or change the alignment gene above.</div></div>`;
      xsBarSub.textContent = `· ${g}`;
      return;
    }
    xsBarSub.textContent = `· ${g} · ${nEmb} zygotes`;
    bars.push({ type: "scatter", mode: "lines", x: connX, y: connY,
      line: { color: "rgba(100,116,139,0.35)", width: 1 }, hoverinfo: "skip", showlegend: false });
    plotInto(xsBars, bars, {
      barmode: "stack", margin: { l: 46, r: 8, t: 6, b: 24 }, height: xsBars.clientHeight || 190,
      xaxis: { tickvals: [0, 1], ticktext: [`left (lower ${state.alignGene})`, `right (higher ${state.alignGene})`],
               range: [-0.5, 1.5], fixedrange: true, tickfont: { size: 10 } },
      yaxis: { title: { text: `${g} count`, font: { size: 10 } }, tickfont: { size: 9 },
               gridcolor: "#eef1f5", fixedrange: true, rangemode: "tozero" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    });
  }
  function renderCrossAgg() {
    return ensureAgg().then(() => {
      const cov = curAGG().gene_cov[state.alignGene] || 0, tot = curAGG().n_embryos;
      xsCaption.textContent = `· best plane ${BESTKEY_LABEL[state.crossKey]}`;
      xsNote.innerHTML = `No gene spans all zygotes (multiple panels). Showing the <b>${cov}/${tot}</b> ` +
        `zygotes that contain <b>${state.alignGene}</b>, each flipped so its higher-${state.alignGene} side is on the right.`;
      renderOutlines(); renderBars();
      requestAnimationFrame(() => { try { Plotly.Plots.resize(xsOutlines); Plotly.Plots.resize(xsBars); } catch (_) {} });
    });
  }

  // ---------- wiring ----------
  function wireControls() {
    geneSelect.addEventListener("change", () => {
      state.userGene = geneSelect.value; render(); renderChart(); highlightBest();
      if (state.drawerOpen && state.agg) renderBars();          // bottom bar plot = selected gene
    });
    planeSelect.addEventListener("change", () => { state.planeIdx = parseInt(planeSelect.value, 10) || 0; render(); renderChart(); });
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
      if (state.scene.genes.includes(g)) { state.userGene = g; geneSelect.value = g; render(); renderChart(); highlightBest(); }
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
    // cross-embryo controls: best-plane alignment + orientation gene
    xsPlane.addEventListener("change", () => { state.crossKey = xsPlane.value; if (state.agg) renderCrossAgg(); });
    xsAlign.addEventListener("change", () => { state.alignGene = xsAlign.value; if (state.agg) renderCrossAgg(); });
    xsMean.addEventListener("change", () => { state.xsShowMean = xsMean.checked; if (state.agg) { renderOutlines(); resizeXs(); } });
    xsOnly.addEventListener("change", () => { state.xsOnlyCurrent = xsOnly.checked; if (state.agg) { renderOutlines(); resizeXs(); } });
  }

  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false;
    placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
