/* Zygote Division Planes — analysis model.
 * Built on viewer-core.js (VCore): validated zygotes, the polar-body axis, 18 candidate
 * division planes, per-gene transcript split (blue/red) across the selected plane,
 * a counts+null chart, best-plane tables, and a cross-section plot. All statistics
 * are pre-computed (build_zygote.py); the UI only reads them.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const el = (tag, cls) => { const e = document.createElement(tag); if (cls) e.className = cls; return e; };
  const V = window.VCore;
  const XY = 0.15;
  const AXIS_C = "#111827", PLANE_C = "#f97316";
  const BLUE = "#2166ac", RED = "#b2182b", GREEN = "#16a34a";

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
  const xsCaption = $("#xs-caption"), xsBarTitle = $("#xs-bar-title"), xsBarSub = $("#xs-bar-sub");
  const xsOutlines = $("#xs-outlines"), xsBars = $("#xs-bars");
  const xsGm = $("#xs-gm"), xsGmSub = $("#xs-gm-sub"), xsGmNote = $("#xs-gm-note"), xsGmDownload = $("#xs-gm-download");
  const xsGmAnchor = $("#xs-gm-anchor"), xsGmReroll = $("#xs-gm-reroll"), xsGmStats = $("#xs-gm-stats");
  const xsTabsEl = $("#xs-tabs"), xsPanels = $("#xs-panels");
  const xsAlign = $("#xs-align"), xsAlignSub = $("#xs-align-sub"), xsAlignDownload = $("#xs-align-download");
  const xsAlignCells = $("#xs-align-cells"), xsAlignMean = $("#xs-align-mean"), xsAlignOnly = $("#xs-align-only");
  const xsAlignPlane = $("#xs-align-plane"), xsAlignLegend = $("#xs-align-legend");
  const xsSp = $("#xs-sp"), xsSpSub = $("#xs-sp-sub"), xsSpNote = $("#xs-sp-note"), xsSpDownload = $("#xs-sp-download");
  const xsSpDensity = $("#xs-sp-density"), xsGmDensity = $("#xs-gm-density"), xsGmSperm = $("#xs-gm-sperm");
  const xsBody = $("#xs-body"), xsPb = $("#xs-pb"), xsPronuclei = $("#xs-pronuclei");
  const xsCrossLegend = $("#xs-cross-legend"), xsCrossDownload = $("#xs-cross-download");
  const xsCrossHighColor = $("#xs-cross-high-color"), xsCrossLowColor = $("#xs-cross-low-color");
  const xsCrossFillColor = $("#xs-cross-fill-color");
  const xsCrossBodyColor = $("#xs-cross-body-color"), xsCrossPbColor = $("#xs-cross-pb-color");
  const xsCrossPnColor = $("#xs-cross-pn-color"), xsCrossScale = $("#xs-cross-scale"), xsCrossFormat = $("#xs-cross-format");
  const xsCrossDotSize = $("#xs-cross-dot-size"), xsCrossDotSizeValue = $("#xs-cross-dot-size-value");
  const xsBarPercent = $("#xs-bar-percent"), xsBarLog = $("#xs-bar-log"), xsBarAdjacent = $("#xs-bar-adjacent");
  const xsBarDensity = $("#xs-bar-density");
  const xsBarLegend = $("#xs-bar-legend"), xsBarNull = $("#xs-bar-null");
  const xsBarInterval = $("#xs-bar-interval"), xsBarGrid = $("#xs-bar-grid"), xsBarDownload = $("#xs-bar-download");
  const xsBarHighColor = $("#xs-bar-high-color"), xsBarLowColor = $("#xs-bar-low-color");
  const xsBarNullColor = $("#xs-bar-null-color"), xsBarScale = $("#xs-bar-scale"), xsBarFormat = $("#xs-bar-format");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle");
  const rtabsEl = $("#rtabs"), bestListEl = $("#best-list");
  const rconcordControls = $("#rconcord-controls"), rconcordRank = $("#rconcord-rank");

  const state = {
    manifest: [], currentId: null, scene: null, userGene: null, planeIdx: 0,
    drawerOpen: false, bestTab: "pVol", crossMode: "vol", xsTab: "sperm", spermData: null,
    crossKey: "pVol", agg: null, pronucleusVisible: {}, dotSize: 1.5,
    gmAnchor: "gene", gmDraw: null, gmDrawKey: null,   // γ/μ grid: real anchor vs random control
    gmDensity: true, spDensity: true,                  // concordance by density (count ÷ side volume)
    concordRank: "frac", _concord: null,               // right-drawer concordance tab
    circ: false, aggCirc: null, _aggCircP: null,
  };
  // circularize accessors: when the "blow up the balloon" toggle is on, every read of
  // the current embryo's analysis / transcripts, and the cross-embryo aggregate, uses
  // the precomputed circularized (spherical, seg-1-only) version instead of the real one.
  const curA = () => (state.circ && state.scene && state.scene.circ) ? state.scene.circ.analysis : (state.scene && state.scene.analysis);
  const curTX = () => (state.circ && state.scene && state.scene.circ) ? state.scene.circ.transcripts : (state.scene && state.scene.transcripts);
  const curAGG = () => (state.circ && state.aggCirc) ? state.aggCirc : state.agg;
  const BESTKEY_LABEL = { pVol: "min p · vol", pCnt: "min p · count", diffVol: "max Δ · vol", diffCnt: "max Δ · count" };
  let vcExtras = null;   // dot-size + atlas-link row (VCore.addWindowExtras)

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
      vcExtras = V.addWindowExtras($("#controls-body"), { defaultSize: state.dotSize, onDotSize: (s) => { state.dotSize = s; render(); } });
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
      if (vcExtras) vcExtras.setAtlas(id);
      controlsEl.hidden = false; placeholder.hidden = true;
      drawer.hidden = false; rdrawer.hidden = false;
      render(); renderChart(); renderBestList();
      if (state.drawerOpen) renderCrossAgg();

    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }
  state.byLabel = (id) => (state.manifest.find((m) => m.id === id) || {}).label || id;
  // Zygote holding the most transcripts of a gene — the most informative one to jump to when
  // a gene is picked that the current zygote doesn't carry. Only offers loadable embryos.
  function mostAbundantEmbryoFor(g) {
    const agg = curAGG(); if (!agg || !agg.embryos) return null;
    const loadable = new Set(state.manifest.map((m) => m.id));
    let bestId = null, bestN = -1;
    for (const e of agg.embryos) {
      const r = e.g[g];
      if (r && r[0] > bestN && loadable.has(e.id)) { bestN = r[0]; bestId = e.id; }
    }
    return bestId;
  }

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
      x: sp.bx, y: sp.by, z: sp.bz, marker: { size: state.dotSize, color: BLUE, opacity: 0.85, line: { width: 0 } },
      hovertemplate: `${g} · side A (counted)<extra></extra>`, legendrank: 20000 });
    traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side B`,
      x: sp.rx, y: sp.ry, z: sp.rz, marker: { size: state.dotSize, color: RED, opacity: 0.85, line: { width: 0 } },
      hovertemplate: `${g} · side B (counted)<extra></extra>`, legendrank: 20001 });
    traces.push({ type: "scatter3d", mode: "markers", name: `${g} · not counted`,
      x: sp.gx, y: sp.gy, z: sp.gz, marker: { size: state.dotSize, color: GREEN, opacity: 0.7, line: { width: 0 } },
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
  // ---------- concordance ranking (right drawer, "Concordance" tab) ----------
  // Every gene taken as the anchor in turn: what share of the panel stays γ on the half that
  // gene defines, and is it beyond the random-plane/random-side null? Doing this per gene the
  // naive way is O(genes² × planes), so the per-plane totals are computed ONCE over all genes
  // and each anchor's own row is just subtracted back out.
  const CONCORD_MIN_COV = 5, CONCORD_DRAWS = 1000;
  function gmAllPlaneStats(embryos) {
    return embryos.map((e) => {
      const gp = e.gp || {}, g = e.g;
      let nP = 0; for (const k in gp) { nP = gp[k].length; break; }
      const cAll = new Int32Array(nP), nAll = new Int32Array(nP);
      for (const name in gp) {
        const total = g[name] && g[name][0]; if (!total) continue;
        const arr = gp[name];
        for (let p = 0; p < nP; p++) {
          const twice = arr[p] * 2;
          if (twice > total) { cAll[p]++; nAll[p]++; } else if (twice < total) nAll[p]++;
        }
      }
      return { cAll, nAll, nP };
    });
  }
  function concordanceRanking() {
    const agg = curAGG(); if (!agg || !agg.embryos) return [];
    const ki = bestKeyIndex();
    const key = `${state.circ ? "circ" : "real"}|${state.crossKey}|${agg.n_embryos}`;
    if (state._concord && state._concord.key === key) return state._concord.rows;
    const embryos = agg.embryos, all = gmAllPlaneStats(embryos), rows = [];
    for (const name of (agg.genes_all || [])) {
      const idx = [];
      embryos.forEach((e, i) => { if (e.gp && e.gp[name] && e.g[name] && e.g[name][0] > 0) idx.push(i); });
      if (idx.length < CONCORD_MIN_COV) continue;
      const stats = idx.map((i) => {
        const base = all[i], e = embryos[i], arr = e.gp[name], total = e.g[name][0];
        const cA = new Int32Array(base.nP), nD = new Int32Array(base.nP);
        for (let p = 0; p < base.nP; p++) {          // drop the anchor's own row from the totals
          const twice = arr[p] * 2;
          cA[p] = base.cAll[p] - (twice > total ? 1 : 0);
          nD[p] = base.nAll[p] - (twice !== total ? 1 : 0);
        }
        return { cA, nD, nP: base.nP };
      });
      const realPlanes = idx.map((i) => embryos[i].best[ki]);
      const realSides = idx.map((i) => { const r = embryos[i].g[name]; return r[1 + ki] * 2 >= r[0]; });
      const st = gmNullTestFrom(stats, realPlanes, realSides, CONCORD_DRAWS);
      if (st) rows.push({ gene: name, n: idx.length, frac: st.sObs, p: st.p, pSide: st.pSide });
    }
    state._concord = { key, rows };
    return rows;
  }
  function renderConcordList() {
    const rows = concordanceRanking().slice();
    const byP = state.concordRank === "p";
    rows.sort(byP ? (a, b) => (a.p - b.p) || (b.frac - a.frac) || a.gene.localeCompare(b.gene)
                  : (a, b) => (b.frac - a.frac) || (a.p - b.p) || a.gene.localeCompare(b.gene));
    const cur = gene();
    let html = `<div class="best-plane-note">Each gene used as the <b>anchor</b>: the share of the panel ` +
      `that stays <b>γ</b> on the half it defines, with the random plane+side p ` +
      `(${CONCORD_DRAWS} draws). Genes in <b>≥${CONCORD_MIN_COV}</b> zygotes, split at each zygote's ` +
      `<b>${BESTKEY_LABEL[state.crossKey]}</b> plane. Ranked by <b>${byP ? "p value" : "concordance"}</b>.</div>`;
    html += `<div class="best-head"><span></span><span>gene</span><span>n</span><span>γ frac</span><span>p</span></div>`;
    html += rows.map((r, i) =>
      `<div class="best-row${r.gene === cur ? " current" : ""}" data-gene="${r.gene}"` +
      ` title="${r.gene} · anchor in ${r.n} zygotes · ${(r.frac * 100).toFixed(1)}% γ · ` +
      `p ${fmtP(r.p)} (side-flip only: ${fmtP(r.pSide)})">` +
      `<span class="best-num">${i + 1}</span>` +
      `<span class="best-gene">${r.gene}</span>` +
      `<span class="best-null">${r.n}</span>` +
      `<span class="best-real">${(r.frac * 100).toFixed(1)}%</span>` +
      `<span class="best-p${r.p <= 0.05 ? " sig" : ""}">${fmtP(r.p)}</span></div>`).join("");
    bestListEl.innerHTML = rows.length ? html
      : `<div class="best-plane-note">No gene reaches ${CONCORD_MIN_COV} zygotes in this aggregate.</div>`;
  }
  function renderBestList() {
    if (rconcordControls) rconcordControls.hidden = state.bestTab !== "concord";
    if (state.bestTab === "concord") {
      bestListEl.innerHTML = `<div class="best-plane-note">Scoring every gene as an anchor…</div>`;
      ensureAgg().then(() => { if (state.bestTab === "concord") renderConcordList(); });
      return;
    }
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
    if (!state._spermP) state._spermP = fetch("data/zygote_sperm.json")   // sperm concordance grid
      .then((r) => r.json()).then((sp) => { state.spermData = sp; if (state.drawerOpen) renderSperm(); }).catch(() => {});
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
    if (!state.scene || !xsOutlines.offsetParent) return;   // skip when its tab is hidden
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
        marker: { color: xsCrossHighColor.value, size: Number(xsCrossDotSize.value), opacity: 0.58 }, hoverinfo: "skip", showlegend: xsCrossLegend.checked },
      { type: "scattergl", mode: "markers", x: sideB.x, y: sideB.y, name: "Side B transcripts",
        marker: { color: xsCrossLowColor.value, size: Number(xsCrossDotSize.value), opacity: 0.58 }, hoverinfo: "skip", showlegend: xsCrossLegend.checked }
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
      dragmode: "pan", margin: { l: 42, r: 10, t: 8, b: 38 }, autosize: true,
      xaxis: { title: { text: "Distance from division plane (µm)", font: { size: 10 } }, range: [-lim, lim],
        scaleanchor: "y", scaleratio: 1, gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
      yaxis: { title: { text: "Polar-body axis (µm)", font: { size: 10 } }, range: [-lim, lim],
        gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
      legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", uirevision: `${s.id}-${state.planeIdx}`,
    }, XS_CFG);
  }
  // ---------- cross-embryo ALIGNED outlines, coloured by significance (p-value) ----------
  // Every zygote's cell-body cross-section (aggregate `outline`, in u,v ⊥ the polar-body axis),
  // rotated so its best plane (current mode) is vertical and flipped so the selected gene's higher
  // side is on +x, coloured by viridis(sigT(p)) — dark = significant split, yellow = n.s.
  const sigT = (p) => Math.max(0, Math.min(1, (Math.log10(Math.max(p, 1e-12)) + 3) / 3));
  function meanOutline(aligned) {
    const K = 120, sum = new Float64Array(K), cnt = new Int32Array(K);
    aligned.forEach((o) => {
      const rBin = new Float64Array(K), has = new Uint8Array(K);
      o.pts.forEach(([x, y]) => {
        const b = Math.min(K - 1, Math.max(0, ((Math.atan2(y, x) + Math.PI) / (2 * Math.PI) * K) | 0)), r = Math.hypot(x, y);
        if (r > rBin[b]) { rBin[b] = r; has[b] = 1; }
      });
      for (let b = 0; b < K; b++) if (has[b]) { sum[b] += rBin[b]; cnt[b]++; }
    });
    const out = [];
    for (let b = 0; b < K; b++) if (cnt[b]) { const r = sum[b] / cnt[b], a = (b + 0.5) / K * 2 * Math.PI - Math.PI; out.push([r * Math.cos(a), r * Math.sin(a)]); }
    if (out.length < 3) return null; out.push(out[0]); return out;
  }
  function renderAlignedOutlines() {
    if (!xsAlign || !xsAlign.offsetParent) return;               // skip when its tab is hidden
    const agg = curAGG(); if (!agg) return;
    const ki = bestKeyIndex(), g = gene(), step = state.step || 10, key = state.crossKey;
    const onlyCur = xsAlignOnly && xsAlignOnly.checked;
    const showCells = (xsAlignCells ? xsAlignCells.checked : true) && !onlyCur;
    const showMean = (xsAlignMean ? xsAlignMean.checked : true) && !onlyCur;
    const showLegend = xsAlignLegend ? xsAlignLegend.checked : true;
    const showPlane = xsAlignPlane ? xsAlignPlane.checked : true;
    const aligned = [];
    agg.embryos.forEach((e) => {
      if (!e.outline || !e.outline.length) return;
      const th = e.best[ki] * step * Math.PI / 180, c = Math.cos(th), s = Math.sin(th);
      const gRow = e.g[g], flip = gRow ? gRow[1 + ki] * 2 < gRow[0] : false;   // higher gene side → +x
      const pts = e.outline.map(([u, v]) => { let nx = u * c + v * s; const ny = -u * s + v * c; if (flip) nx = -nx; return [nx, ny]; });
      pts.push(pts[0]);
      aligned.push({ pts, p: (e.sig && e.sig[key] != null) ? e.sig[key] : 1, id: e.id, label: e.label, isCur: e.id === state.currentId });
    });
    if (!aligned.length) {
      Plotly.purge(xsAlign); xsAlign.classList.remove("js-plotly-plot");
      xsAlign.innerHTML = `<div class="xs-empty"><div>No aligned outlines for this selection.</div></div>`;
      if (xsAlignSub) xsAlignSub.textContent = `· ${g}`; return;
    }
    const traces = []; let lim = 20;
    const bump = (o) => o.pts.forEach(([x, y]) => { const m = Math.max(Math.abs(x), Math.abs(y)); if (m > lim) lim = m; });
    if (showCells) aligned.forEach((o) => { if (o.isCur) return; bump(o);
      traces.push({ type: "scatter", mode: "lines", x: o.pts.map((q) => q[0]), y: o.pts.map((q) => q[1]),
        line: { color: viridis(sigT(o.p)), width: 1, shape: "spline" }, opacity: 0.5, hoverinfo: "skip", showlegend: false }); });
    if (showMean) { const mo = meanOutline(aligned); if (mo) traces.push({ type: "scatter", mode: "lines",
      x: mo.map((q) => q[0]), y: mo.map((q) => q[1]), name: "mean outline",
      line: { color: "#0f172a", width: 2.6, shape: "spline" }, hoverinfo: "skip", showlegend: showLegend }); }
    const cur = aligned.find((o) => o.isCur);
    if (cur) { bump(cur);
      traces.push({ type: "scatter", mode: "lines", x: cur.pts.map((q) => q[0]), y: cur.pts.map((q) => q[1]),
        line: { color: "#fff", width: 5, shape: "spline" }, hoverinfo: "skip", showlegend: false });
      traces.push({ type: "scatter", mode: "lines", x: cur.pts.map((q) => q[0]), y: cur.pts.map((q) => q[1]),
        name: `${cur.label} (this embryo)`, line: { color: viridis(sigT(cur.p)), width: 2.6, shape: "spline" }, hoverinfo: "skip", showlegend: showLegend }); }
    lim *= 1.08;
    if (showPlane) traces.push({ type: "scatter", mode: "lines", x: [0, 0], y: [-lim, lim], name: "Division plane",
      line: { color: "#111827", width: 1.4, dash: "dash" }, hoverinfo: "skip", showlegend: showLegend });
    if (xsAlignSub) xsAlignSub.textContent = `· ${g} · ${aligned.length} zygotes · ${BESTKEY_LABEL[key]} plane`;
    plotInto(xsAlign, traces, {
      dragmode: "pan", margin: { l: 42, r: 10, t: 8, b: 38 }, autosize: true, showlegend: showLegend,
      xaxis: { title: { text: "Distance from division plane (µm)", font: { size: 10 } }, range: [-lim, lim],
        scaleanchor: "y", scaleratio: 1, gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
      yaxis: { range: [-lim, lim], gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
      legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
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
    if (!xsBars.offsetParent) return;                       // skip when its tab is hidden
    const agg = curAGG(); if (!agg) return;
    const ki = bestKeyIndex(), g = gene();
    // DENSITY mode: each side's count ÷ that side's volume (adjacent bars); needs per-side volumes.
    const density = !!(xsBarDensity && xsBarDensity.checked && agg.embryos[0] && agg.embryos[0].vp);
    const percentMode = !density && xsBarPercent.checked;
    const DSCALE = 1e4;                                      // show density as ×10⁻⁴ per µm³
    const rows = agg.embryos.map((emb) => {
      const row = emb.g[g]; if (!row) return null;
      const a = row[1 + ki], b = row[0] - a, totalCount = a + b; if (!totalCount) return null;
      const highCount = Math.max(a, b), lowCount = Math.min(a, b);
      if (density) {
        const p = emb.best[ki], vA = emb.vp[p], vB = Math.max(emb.vt - emb.vp[p], 1);
        const dA = a / vA * DSCALE, dB = b / vB * DSCALE, hi = Math.max(dA, dB), lo = Math.min(dA, dB);
        return { label: emb.label, high: hi, low: lo, total: hi + lo, highCount, lowCount, totalCount,
          nullMean: 0, nullLow: 0, nullHigh: 0 };
      }
      const scale = percentMode ? 100 / totalCount : 1;
      const [nullLowCount, nullHighCount] = binomial95(totalCount);
      return { label: emb.label, high: highCount * scale, low: lowCount * scale,
        total: totalCount * scale, highCount, lowCount, totalCount,
        nullMean: (totalCount / 2) * scale, nullLow: nullLowCount * scale,
        nullHigh: nullHighCount * scale };
    }).filter(Boolean);
    if (!rows.length) {
      Plotly.purge(xsBars); xsBars.classList.remove("js-plotly-plot");   // so plotInto re-inits cleanly next time
      xsBars.innerHTML = `<div class="xs-empty"><div><b>${g}</b> is not detected in any retained zygote.</div></div>`;
      xsBarSub.textContent = `· ${g}`;
      return;
    }
    const useLog = !percentMode && xsBarLog.checked;
    // Bar centres are offset by BAR_X0 (and the axis starts at 0) so the first bar sits clear
    // of the vertical axis instead of straddling it. Widest element is the null bar (±0.34),
    // so BAR_X0 = 0.6 leaves the same 0.26 gap at both ends of the axis.
    const BAR_X0 = 0.6;
    const x = rows.map((_, i) => i + BAR_X0), baseline = useLog ? 1 : 0, shapes = [], traces = [];
    const stacked = !density && !xsBarAdjacent.checked;     // density is always adjacent side-by-side
    rows.forEach((r, i) => {
      const cx = i + BAR_X0;
      if (!stacked) {
        shapes.push(
          { type: "rect", xref: "x", yref: "y", layer: "above", x0: cx - 0.24, x1: cx - 0.03,
            y0: baseline, y1: r.high, fillcolor: xsBarHighColor.value, line: { color: "rgba(15,23,42,0.48)", width: 0.65 } },
          { type: "rect", xref: "x", yref: "y", layer: "above", x0: cx + 0.03, x1: cx + 0.24,
            y0: baseline, y1: r.low, fillcolor: xsBarLowColor.value, line: { color: "rgba(15,23,42,0.48)", width: 0.65 } }
        );
      } else {
        shapes.push(
          { type: "rect", xref: "x", yref: "y", layer: "above", x0: cx - 0.23, x1: cx + 0.23,
            y0: baseline, y1: r.high, fillcolor: xsBarHighColor.value, line: { color: "rgba(15,23,42,0.48)", width: 0.65 } },
          { type: "rect", xref: "x", yref: "y", layer: "above", x0: cx - 0.23, x1: cx + 0.23,
            y0: r.high, y1: r.total, fillcolor: xsBarLowColor.value, line: { color: "rgba(15,23,42,0.48)", width: 0.65 } }
        );
      }
      if (!density && xsBarNull.checked) shapes.push({ type: "rect", xref: "x", yref: "y", layer: "above",
        x0: cx - 0.34, x1: cx + 0.34, y0: baseline, y1: r.nullMean,
        fillcolor: hexAlpha(xsBarNullColor.value, 0.30), line: { color: "#5f666d", width: 0.55 } });
    });
    if (xsBarLegend.checked) {
      traces.push(
        { type: "scatter", mode: "markers", name: density ? "Higher-density half" : "Higher-count half", legendrank: 10, x: [null], y: [null],
          marker: { symbol: "square", size: 11, color: xsBarHighColor.value }, hoverinfo: "skip" },
        { type: "scatter", mode: "markers", name: density ? "Lower-density half" : "Lower-count half", legendrank: 20, x: [null], y: [null],
          marker: { symbol: "square", size: 11, color: xsBarLowColor.value }, hoverinfo: "skip" }
      );
      if (!density && xsBarNull.checked) traces.push({ type: "scatter", mode: "markers", name: percentMode ? "50% null expectation" : "Null mean", legendrank: 30,
        x: [null], y: [null], marker: { symbol: "square", size: 11, color: hexAlpha(xsBarNullColor.value, 0.6) }, hoverinfo: "skip" });
    }
    if (!density && xsBarInterval.checked) traces.push({ type: "scatter", mode: "markers", name: "95% null interval", x,
      legendrank: 40,
      y: rows.map((r) => r.nullMean), marker: { size: 1, color: "rgba(0,0,0,0)" },
      error_y: { type: "data", symmetric: false, array: rows.map((r) => r.nullHigh - r.nullMean),
        arrayminus: rows.map((r) => r.nullMean - r.nullLow), color: "#111827", thickness: 0.8, width: 2.5 },
      hovertemplate: "%{x}<br>95% null interval<extra></extra>" });
    traces.push({ type: "scatter", mode: "markers", showlegend: false, x,
      y: rows.map((r) => stacked ? r.total : Math.max(r.high, r.low)),
      customdata: rows.map((r) => [r.label, r.high, r.low, r.total, r.nullMean, r.nullLow, r.nullHigh,
        r.highCount, r.lowCount, r.totalCount]),
      marker: { size: 24, color: "rgba(0,0,0,0)" },
      hovertemplate: density
        ? "%{customdata[0]}<br>higher-density half %{customdata[1]:.1f}<br>lower-density half %{customdata[2]:.1f}" +
          " ×10⁻⁴/µm³<br>counts %{customdata[7]} / %{customdata[8]}<extra></extra>"
        : percentMode
          ? "%{customdata[0]}<br>higher half %{customdata[1]:.1f}% (%{customdata[7]} transcripts)" +
            "<br>lower half %{customdata[2]:.1f}% (%{customdata[8]} transcripts)" +
            "<br>95% null interval %{customdata[5]:.1f}%–%{customdata[6]:.1f}%<extra></extra>"
          : "%{customdata[0]}<br>higher half %{customdata[1]}<br>lower half %{customdata[2]}" +
            "<br>total %{customdata[3]}<br>null %{customdata[4]:.1f} (%{customdata[5]}–%{customdata[6]})<extra></extra>" });
    xsBarTitle.textContent = density ? "Per-side density" : percentMode ? "Per-side percentage" : "Per-side counts";
    xsBarSub.textContent = `· ${g} · ${rows.length} zygotes${density ? " · count ÷ side volume" : ""}`;
    const maxY = percentMode ? 100 : Math.max(...rows.map((r) => stacked
      ? Math.max(r.total, r.nullHigh)
      : Math.max(r.high, r.low, r.nullHigh))) * 1.18;
    plotInto(xsBars, traces, {
      shapes, margin: { l: 56, r: 10, t: xsBarLegend.checked ? 48 : 8, b: 92 }, autosize: true,
      showlegend: xsBarLegend.checked,
      legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      xaxis: { fixedrange: false, tickmode: "array", tickvals: x, ticktext: rows.map((r) => r.label),
        range: [0, rows.length - 1 + 2 * BAR_X0], tickangle: -48, tickfont: { size: 8 }, automargin: true },
      yaxis: { title: { text: density ? `${g} density (×10⁻⁴ per µm³)` : percentMode ? `${g} share of zygote transcripts (%)` : `${g} transcript count`, font: { size: 10 } }, tickfont: { size: 9 },
               type: useLog ? "log" : "linear", range: useLog ? [0, Math.log10(maxY)] : [0, maxY],
               dtick: useLog ? 1 : undefined, ticksuffix: percentMode ? "%" : "",
               gridcolor: xsBarGrid.checked ? "#e2e5e8" : "rgba(0,0,0,0)", fixedrange: false },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, XS_CFG);
  }
  // γ / μ concordance grid: the anchor gene (floating-window selection) splits each zygote
  // by its best plane into a γ half (more anchor transcripts) and a μ half (fewer). Then, for
  // EVERY gene, each cell asks: on that same γ half, is this gene also enriched (γ) or does it
  // flip to the μ side (μ)? Rows = genes, columns = zygotes containing the anchor gene.
  // NEGATIVE CONTROL (Harry's null test): instead of an anchor gene, give each zygote a RANDOM
  // plane and call one side γ at random. Concordance is then scored against that null to get a
  // p-value. Row ORDER always comes from the real anchor — re-sorting the control by its own
  // γ-fraction would manufacture a gradient out of pure noise and make the null look structured.
  // γ if the gene sits higher on the anchor's γ half — by raw count, or by DENSITY (count ÷ side
  // volume) when `density` + the per-side volumes (volA/volB) are supplied.
  function gmCall(a, total, gammaSideA, volA, volB, density) {
    if (a == null || !total) return null;           // gene absent in this zygote
    const b = total - a;
    const vA = density ? a / volA : a, vB = density ? b / volB : b;
    const gv = gammaSideA ? vA : vB, ov = gammaSideA ? vB : vA;
    return gv > ov ? "G" : (gv < ov ? "M" : "T");
  }
  function gmNPlanes(cols) {
    for (const e of cols) { const gp = e.gp; if (gp) { for (const k in gp) return gp[k].length; } }
    return 0;                                       // aggregate predates the per-plane counts
  }
  function gmColsKey(cols) { return cols.map((e) => e.id).join("|"); }
  function gmEnsureDraw(cols) {
    const key = gmColsKey(cols), nP = gmNPlanes(cols);
    if (!state.gmDraw || state.gmDrawKey !== key) {
      state.gmDraw = cols.map(() => ({ p: nP ? (Math.random() * nP) | 0 : 0, sideA: Math.random() < 0.5 }));
      state.gmDrawKey = key;
    }
    return state.gmDraw;
  }
  // per zygote, per plane: how many genes sit on side A and how many are decided (non-tie).
  // The anchor's own row is excluded — it is γ by construction and would inflate concordance.
  function gmPlaneStats(cols, anchor, density) {
    return cols.map((e) => {
      const gp = e.gp || {}, g = e.g, vp = e.vp, vt = e.vt, canD = !!(density && vp && vt);
      let nP = 0; for (const k in gp) { nP = gp[k].length; break; }
      const cA = new Int32Array(nP), nD = new Int32Array(nP);
      for (const name in gp) {
        if (name === anchor) continue;
        const total = g[name] && g[name][0]; if (!total) continue;
        const arr = gp[name];
        for (let p = 0; p < nP; p++) {
          const a = arr[p], b = total - a; let onA;
          if (canD) { const dA = a / vp[p], dB = b / Math.max(vt - vp[p], 1); onA = dA > dB ? 1 : (dA < dB ? -1 : 0); }
          else { const twice = a * 2; onA = twice > total ? 1 : (twice < total ? -1 : 0); }
          if (onA > 0) { cA[p]++; nD[p]++; } else if (onA < 0) nD[p]++;
        }
      }
      return { cA, nD, nP };
    });
  }
  function gmGammaFrac(stats, planes, sides) {
    let gam = 0, tot = 0;
    for (let i = 0; i < stats.length; i++) {
      const s = stats[i], p = planes[i];
      if (!s.nP || p >= s.nP) continue;
      const a = s.cA[p], n = s.nD[p];
      gam += sides[i] ? a : n - a; tot += n;
    }
    return tot ? gam / tot : 0.5;
  }
  // Sign-flip / random-plane permutation test. Cells within a zygote are NOT independent (one
  // γ side is chosen per zygote), so the null must randomise per zygote — a naive per-cell
  // binomial would wildly overstate significance.
  function gmNullTest(cols, anchor, realPlanes, realSides, N, density) {
    return gmNullTestFrom(gmPlaneStats(cols, anchor, density), realPlanes, realSides, N);
  }
  function gmNullTestFrom(stats, realPlanes, realSides, N) {
    const nP = stats.length ? stats[0].nP : 0;
    if (!nP) return null;
    const sObs = gmGammaFrac(stats, realPlanes, realSides), dev = Math.abs(sObs - 0.5);
    const planes = new Array(stats.length), sides = new Array(stats.length);
    let ge = 0, geSide = 0, sum = 0, sum2 = 0;
    for (let it = 0; it < N; it++) {
      for (let i = 0; i < stats.length; i++) { planes[i] = (Math.random() * nP) | 0; sides[i] = Math.random() < 0.5; }
      const s = gmGammaFrac(stats, planes, sides);
      sum += s; sum2 += s * s;
      if (Math.abs(s - 0.5) >= dev) ge++;
      // conservative variant: keep each zygote's real best plane, flip only the γ side
      for (let i = 0; i < stats.length; i++) sides[i] = Math.random() < 0.5;
      if (Math.abs(gmGammaFrac(stats, realPlanes, sides) - 0.5) >= dev) geSide++;
    }
    const mean = sum / N, sd = Math.sqrt(Math.max(0, sum2 / N - mean * mean));
    return { sObs, nullMean: mean, nullSd: sd, p: (1 + ge) / (1 + N), pSide: (1 + geSide) / (1 + N), N };
  }
  function gmModel() {
    const agg = curAGG(); if (!agg) return null;
    const A = gene(), ki = bestKeyIndex();
    const cols = agg.embryos.filter((e) => e.g[A] && e.g[A][0] > 0);
    const isNull = state.gmAnchor === "null";
    if (!cols.length) return { anchor: A, ki, isNull, cols: [], rows: [], stats: null, hasGp: false };
    const hasGp = gmNPlanes(cols) > 0;
    const density = !!(state.gmDensity && cols[0].vp && cols[0].vt);   // needs per-side volumes
    const volAt = (e, p) => (e.vp ? e.vp[p] : 1);
    const volBt = (e, p) => (e.vp ? Math.max(e.vt - e.vp[p], 1) : 1);
    // real anchor: each zygote's best plane + the side the anchor is enriched on (count or density)
    const realPlanes = cols.map((e) => e.best[ki]);
    const realSides = cols.map((e, i) => {
      const r = e.g[A], a = r[1 + ki], b = r[0] - a, p = realPlanes[i];
      const vA = density ? a / volAt(e, p) : a, vB = density ? b / volBt(e, p) : b;
      return vA >= vB;
    });
    const useNull = isNull && hasGp;
    const draw = useNull ? gmEnsureDraw(cols) : null;
    const planes = useNull ? draw.map((d) => d.p) : realPlanes;
    const sides = useNull ? draw.map((d) => d.sideA) : realSides;
    const aOf = (e, i, g) => {
      if (useNull) { const arr = e.gp && e.gp[g]; return arr ? arr[planes[i]] : null; }
      return e.g[g] ? e.g[g][1 + ki] : null;
    };
    const totalOf = (e, g) => (e.g[g] ? e.g[g][0] : 0);

    const geneSet = new Set(); cols.forEach((e) => Object.keys(e.g).forEach((g) => geneSet.add(g)));
    // fixed row order, always from the REAL anchor (see note above)
    const realFrac = {};
    for (const g of geneSet) {
      let nG = 0, n = 0;
      cols.forEach((e, i) => {
        const p = realPlanes[i];
        const c = gmCall(e.g[g] ? e.g[g][1 + ki] : null, totalOf(e, g), realSides[i], volAt(e, p), volBt(e, p), density);
        if (c === "G") { nG++; n++; } else if (c === "M") n++;
      });
      realFrac[g] = n ? nG / n : 0;
    }
    const rowNames = [A, ...[...geneSet].filter((g) => g !== A)
      .sort((x, y) => (realFrac[y] - realFrac[x]) || x.localeCompare(y))];
    const rows = rowNames.map((g) => {
      const cells = cols.map((e, i) => gmCall(aOf(e, i, g), totalOf(e, g), sides[i], volAt(e, planes[i]), volBt(e, planes[i]), density));
      const present = cells.filter((c) => c !== null);
      const nG = present.filter((c) => c === "G").length;
      return { gene: g, cells, cov: present.length,
               gammaFrac: present.length ? nG / present.length : 0, realFrac: realFrac[g] };
    });
    return { anchor: A, ki, isNull, hasGp, density, cols, rows, realSides,
             stats: gmNullTest(cols, A, realPlanes, realSides, 2000, density) };
  }
  // NOTE: injected via innerHTML, so the "<" must be escaped.
  function gmFmtP(p) { return p < 0.001 ? "&lt; 0.001" : "= " + p.toFixed(3); }
  function renderGammaMu() {
    if (!xsGm || !xsGm.offsetParent) return;                // skip when its tab is hidden
    const m = gmModel();
    if (!m) return;
    const A = m.anchor, showingNull = m.isNull && m.hasGp;
    if (xsGmReroll) xsGmReroll.disabled = !showingNull;
    if (!m.cols.length) {
      xsGm.dataset.anchor = A;
      xsGm.innerHTML = `<div class="xs-empty"><div><b>${A}</b> is not detected in any retained zygote.</div></div>`;
      if (xsGmSub) xsGmSub.textContent = `· ${A}`;
      if (xsGmNote) xsGmNote.textContent = "";
      if (xsGmStats) xsGmStats.textContent = "";
      return;
    }
    if (xsGmNote) xsGmNote.innerHTML = m.isNull && !m.hasGp
      ? `<b>Negative control unavailable:</b> this aggregate predates the per-plane counts — rerun ` +
        `<code>scripts/add_plane_counts_to_cross_agg.py</code> to enable it.`
      : showingNull
        ? `<b>Negative control</b> — no anchor gene: every zygote uses a <b>random plane</b> and a ` +
          `<b>randomly chosen γ side</b>. Rows keep the real anchor's order, so real structure would ` +
          `still line up here; noise will not. Hit <b>Re-roll</b> for another draw.`
        : `Each zygote is split at its <b>${BESTKEY_LABEL[state.crossKey]}</b> plane; the half with the higher ` +
          `<b>${A}</b> ${m.density ? "<b>density</b> (count ÷ side volume)" : "count"} is <b>γ</b>. Every gene is then ` +
          `<b>γ</b> if it is also higher on that half, <b>μ</b> if it flips.`;
    if (xsGmStats) {
      const st = m.stats;
      xsGmStats.innerHTML = st
        ? `<b>${A}</b>: observed <b>${(st.sObs * 100).toFixed(1)}% γ</b> vs null ` +
          `<b>${(st.nullMean * 100).toFixed(1)}% ± ${(st.nullSd * 100).toFixed(1)}</b> · ` +
          `<b>p ${gmFmtP(st.p)}</b> <span class="xs-gm-p2">(${st.N} random plane+side draws; ` +
          `side-flip only: p ${gmFmtP(st.pSide)})</span>`
        : "";
    }
    const frag = document.createDocumentFragment();
    const corner = el("div", "xs-gm-corner"); corner.textContent = "gene \\ zygote"; frag.appendChild(corner);
    m.cols.forEach((e) => { const c = el("div", "xs-gm-collabel"); const s = document.createElement("span"); s.textContent = e.label; c.appendChild(s); frag.appendChild(c); });
    m.rows.forEach((rd) => {
      const isA = rd.gene === A && !showingNull;    // no anchor exists in the control
      const lab = el("div", "xs-gm-rowlabel" + (isA ? " anchor" : ""));
      lab.textContent = rd.gene;
      lab.title = isA ? `${rd.gene} — anchor (defines the γ half)` :
        `${rd.gene} — γ in ${Math.round(rd.gammaFrac * 100)}% of ${rd.cov} zygote(s)` +
        (showingNull ? ` (control draw; real anchor: ${Math.round(rd.realFrac * 100)}%)` : "");
      frag.appendChild(lab);
      rd.cells.forEach((c, i) => {
        const kind = c === "G" ? "g-G" : c === "M" ? "g-M" : c === "T" ? "g-T" : "g-NA";
        const box = el("div", "xs-gm-cell " + kind);
        box.textContent = c === "G" ? "γ" : c === "M" ? "μ" : c === "T" ? "·" : "";
        box.title = `${rd.gene} · ${m.cols[i].label}: ` +
          (c === "G" ? "γ — shares the anchor's γ half" : c === "M" ? "μ — flipped to the μ half" :
           c === "T" ? "tie (50/50)" : "gene not present");
        frag.appendChild(box);
      });
    });
    // optional bottom SPERM row — gene-relative: is the sperm on the anchor's γ half (Ω) or μ half (Σ)?
    // Depends on the anchor gene (via realSides), so it updates whenever the gene changes.
    let spN = 0, spOm = 0;
    if (xsGmSperm && xsGmSperm.checked && !showingNull && state.spermData && m.realSides) {
      const spById = {}; state.spermData.embryos.forEach((r) => (spById[r.id] = r));
      const lab = el("div", "xs-gm-rowlabel anchor"); lab.textContent = "Sperm side";
      lab.title = "Ω = sperm on the anchor gene's γ (higher) half · Σ = its μ (lower) half";
      frag.appendChild(lab);
      m.cols.forEach((e, i) => {
        const sp = spById[e.id], sm = sp && sp.modes[state.crossKey];
        let kind = "g-NA", glyph = "", title = `${e.label}: no located sperm`;
        if (sm) {
          const spermSideA = sm.spermSide > 0, onGamma = spermSideA === m.realSides[i];
          kind = onGamma ? "g-O" : "g-S"; glyph = onGamma ? "Ω" : "Σ"; spN++; if (onGamma) spOm++;
          title = `${e.label}: sperm on the anchor's ${onGamma ? "γ (higher)" : "μ (lower)"} half of ${A}`;
        }
        const box = el("div", "xs-gm-cell " + kind); box.textContent = glyph; box.title = title;
        frag.appendChild(box);
      });
    }
    xsGm.style.setProperty("--gm-cols", m.cols.length);
    xsGm.dataset.anchor = A;
    xsGm.innerHTML = ""; xsGm.appendChild(frag);
    if (xsGmSub) xsGmSub.textContent = `· ${A} · ${m.rows.length} genes × ${m.cols.length} zygotes` +
      (spN ? ` · sperm Ω ${spOm}/${spN} on the γ half` : "");
  }
  function downloadGammaMuCSV() {
    const m = gmModel(); if (!m || !m.cols.length) return;
    const glyph = (c) => (c === "G" ? "gamma" : c === "M" ? "mu" : c === "T" ? "tie" : "");
    const showingNull = m.isNull && m.hasGp, st = m.stats;
    const lines = [];
    lines.push(`# anchor=${showingNull ? "RANDOM negative control (random plane + random side)" : m.anchor}`);
    lines.push(`# plane=${BESTKEY_LABEL[state.crossKey]}${showingNull ? " (control uses random planes)" : ""}`);
    if (st) lines.push(`# observed_gamma_fraction=${st.sObs.toFixed(4)} null_mean=${st.nullMean.toFixed(4)} ` +
      `null_sd=${st.nullSd.toFixed(4)} p_random_plane_and_side=${st.p.toFixed(4)} p_side_flip_only=${st.pSide.toFixed(4)} draws=${st.N}`);
    lines.push(["gene", ...m.cols.map((e) => e.label.replace(/,/g, ""))].join(","));
    m.rows.forEach((rd) => lines.push([rd.gene, ...rd.cells.map(glyph)].join(",")));
    const blob = new Blob([lines.join("\n")], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = `gamma_mu_${showingNull ? "NULL_control" : m.anchor}.csv`;
    document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(a.href);
  }
  function renderCrossAgg() {
    return ensureAgg().then(() => {
      const cov = curAGG().gene_cov[gene()] || 0, tot = curAGG().n_embryos;
      drawerEmb.textContent = `· ${gene()}`;
      xsNote.innerHTML = `Counts use each zygote's <b>${BESTKEY_LABEL[state.crossKey]}</b> plane; <b>${cov}/${tot}</b> retained zygotes contain <b>${gene()}</b>.`;
      renderCurrentCrossSection(); renderAlignedOutlines(); renderBars(); renderGammaMu(); renderSperm();   // each self-gates on its tab
      requestAnimationFrame(() => { try { Plotly.Plots.resize(xsOutlines); Plotly.Plots.resize(xsBars); } catch (_) {} });
    });
  }
  // ---------- bottom-drawer tabs (one graph per tab) ----------
  const XS_RENDER = { cross: renderCurrentCrossSection, align: renderAlignedOutlines, bars: renderBars, gm: renderGammaMu, sperm: renderSperm };
  function switchXsTab(which) {
    if (!XS_RENDER[which]) which = "cross";
    state.xsTab = which;
    xsTabsEl.querySelectorAll(".xs-gtab").forEach((t) => {
      const on = t.dataset.tab === which; t.classList.toggle("active", on); t.setAttribute("aria-selected", String(on));
    });
    xsPanels.querySelectorAll(".xs-panel").forEach((p) => (p.hidden = p.dataset.tab !== which));
    if (state.agg || which === "cross") XS_RENDER[which]();
    requestAnimationFrame(() => { try { Plotly.Plots.resize(xsOutlines); Plotly.Plots.resize(xsBars); } catch (_) {} });
  }

  // ---------- sperm concordance grid ----------
  // For each zygote with a located sperm, mark whether the sperm sits on the side of its current
  // best plane that holds MORE transcripts (Ω omega, green) or FEWER (Σ sigma, purple). Data from
  // build_zygote_sperm.py (data/zygote_sperm.json); the mark follows the "Compare each zygote at" plane.
  function binomTwoSided(k, n) {
    if (!n) return NaN;
    const lg = (m) => { let s = 0; for (let i = 2; i <= m; i++) s += Math.log(i); return s; };
    const pmf = (i) => Math.exp(lg(n) - lg(i) - lg(n - i) - n * Math.LN2);
    const kk = Math.min(k, n - k); let tail = 0;
    for (let i = 0; i <= kk; i++) tail += pmf(i);
    return Math.min(1, 2 * tail);
  }
  function renderSperm() {
    if (!xsSp || !xsSp.offsetParent) return;
    const sp = state.spermData;
    if (!sp) { xsSp.innerHTML = `<div class="xs-gm-empty">Loading sperm data…</div>`; return; }
    const mode = state.crossKey;
    const density = xsSpDensity ? xsSpDensity.checked : state.spDensity;
    const cross = curAGG(); const known = cross ? new Set(cross.embryos.map((e) => e.id)) : null;
    // only zygotes that have sperm AND are in the current (real/circularized) cross set
    const rows = sp.embryos.filter((r) => r.modes && r.modes[mode] && (!known || known.has(r.id)));
    if (!rows.length) { xsSp.innerHTML = `<div class="xs-gm-empty">No zygotes with a located sperm.</div>`; xsSpSub.textContent = ""; return; }
    // Ω (omega) = sperm on the side with GREATER transcript density (or count); Σ (sigma) = lesser
    const vSperm = (m) => (density ? m.cntSperm / (m.volSperm || 1) : m.cntSperm);
    const vOther = (m) => (density ? m.cntOther / (m.volOther || 1) : m.cntOther);
    const isOmega = (m) => vSperm(m) > vOther(m);
    const om = rows.filter((r) => isOmega(r.modes[mode])).length, sg = rows.length - om;
    const p = binomTwoSided(Math.min(om, sg), rows.length);
    const unit = density ? "density" : "transcripts";
    xsSpSub.textContent = `· ${BESTKEY_LABEL[mode]} plane · by ${density ? "density" : "count"} · total transcripts (not gene-specific) · ${rows.length} zygotes`;
    xsSpNote.innerHTML =
      `Sperm on the <b style="color:#16a34a">greater-${unit}</b> side (Ω) in <b>${om}</b> of ${rows.length}; ` +
      `on the <b style="color:#7c3aed">lesser</b> side (Σ) in <b>${sg}</b>. ` +
      `Binomial vs 50/50: <b>p ${p < 0.001 ? "< 0.001" : "= " + p.toFixed(3)}</b>.`;
    // grid: a row-label column + one column per zygote; top row = Ω/Σ marks, below = zygote labels
    xsSp.style.gridTemplateColumns = `minmax(74px,auto) repeat(${rows.length}, minmax(20px, 1fr))`;
    const cell = (r) => {
      const m = r.modes[mode], om1 = isOmega(m); const cls = om1 ? "g-O" : "g-S", gl = om1 ? "Ω" : "Σ";
      const share = m.cntSperm + m.cntOther ? (100 * m.cntSperm / (m.cntSperm + m.cntOther)).toFixed(0) : "–";
      const dS = (m.cntSperm / (m.volSperm || 1) * 1e4).toFixed(2), dO = (m.cntOther / (m.volOther || 1) * 1e4).toFixed(2);
      return `<div class="xs-gm-cell ${cls}" title="${r.label}: sperm on the ${om1 ? "greater" : "lesser"}-${unit} side — ${share}% of transcripts on the sperm's side · density sperm ${dS} vs other ${dO} (×10⁻⁴/µm³, plane ${m.plane})">${gl}</div>`;
    };
    xsSp.innerHTML =
      `<div class="xs-gm-corner"></div>` +
      rows.map((r) => `<div class="xs-gm-collabel"><span>${r.label}</span></div>`).join("") +
      `<div class="xs-gm-rowlabel anchor">Sperm side</div>` +
      rows.map(cell).join("");
  }

  // ---------- wiring ----------
  function wireControls() {
    geneSelect.addEventListener("change", () => {
      state.userGene = geneSelect.value; render(); renderChart(); highlightBest();
      if (state.drawerOpen) renderCrossAgg();
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
      if (state.drawerOpen) renderCrossAgg();
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
    if (rconcordRank) rconcordRank.addEventListener("change", () => {
      state.concordRank = rconcordRank.value;
      if (state.bestTab === "concord") renderConcordList();
    });
    bestListEl.addEventListener("click", (e) => {
      const row = e.target.closest(".best-row"); if (!row) return;
      const g = row.dataset.gene;
      state.userGene = g;                       // remembered, so it survives an embryo switch
      if (state.scene && state.scene.genes.includes(g)) {
        geneSelect.value = g; render(); renderChart(); highlightBest();
        if (state.drawerOpen) renderCrossAgg();
        return;
      }
      // The gene dropdown only holds THIS zygote's genes, so a cross-embryo pick (the
      // Concordance tab especially) can't just be selected — jump to the zygote where the
      // gene is most abundant; populateGenes() then picks it up from state.userGene.
      ensureAgg().then(() => {
        const id = mostAbundantEmbryoFor(g);
        if (id && id !== state.currentId) selectEmbryo(id);
      });
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
    const syncBarModeControls = () => {
      const dens = !!(xsBarDensity && xsBarDensity.checked);   // density → percent/adjacent don't apply
      xsBarPercent.disabled = dens; xsBarAdjacent.disabled = dens;
      xsBarLog.disabled = xsBarPercent.checked && !dens;
      if (xsBarLog.disabled) xsBarLog.checked = false;
      renderBars();
    };
    syncBarModeControls();
    xsPlane.addEventListener("change", () => { state.crossKey = xsPlane.value;
      if (state.agg) renderCrossAgg();
      if (state.bestTab === "concord") renderBestList();   // ranking depends on the plane choice
    });
    [xsBody, xsPb, xsCrossLegend, xsCrossHighColor, xsCrossLowColor, xsCrossFillColor, xsCrossBodyColor, xsCrossPbColor, xsCrossPnColor]
      .forEach((el) => el.addEventListener("change", renderCurrentCrossSection));
    xsCrossDotSize.addEventListener("input", () => {
      xsCrossDotSizeValue.value = xsCrossDotSize.value;
      renderCurrentCrossSection();
    });
    [xsAlignCells, xsAlignMean, xsAlignOnly, xsAlignPlane, xsAlignLegend]
      .forEach((el) => el && el.addEventListener("change", renderAlignedOutlines));
    if (xsAlignDownload) xsAlignDownload.addEventListener("click", () => {
      if (xsAlign.classList.contains("js-plotly-plot"))
        Plotly.downloadImage(xsAlign, { format: "png", scale: 4, width: 1800, height: 1400,
          filename: `aligned_outlines_${(gene() || "gene").replace(/[^a-z0-9]+/gi, "_")}` });
    });
    xsBarPercent.addEventListener("change", syncBarModeControls);
    if (xsBarDensity) xsBarDensity.addEventListener("change", syncBarModeControls);
    [xsBarLog, xsBarAdjacent, xsBarLegend, xsBarNull, xsBarInterval, xsBarGrid, xsBarHighColor, xsBarLowColor, xsBarNullColor]
      .forEach((el) => el.addEventListener("change", renderBars));
    xsCrossDownload.addEventListener("click", () => downloadPlot(xsOutlines, "cross_section", xsCrossFormat, xsCrossScale, 1800, 1400));
    xsBarDownload.addEventListener("click", () => downloadPlot(xsBars,
      xsBarPercent.checked ? "side_percentages" : "side_counts", xsBarFormat, xsBarScale, 2000, 1250));
    if (xsGmDownload) xsGmDownload.addEventListener("click", downloadGammaMuCSV);
    if (xsGmAnchor) xsGmAnchor.addEventListener("change", () => {
      state.gmAnchor = xsGmAnchor.value; renderGammaMu();
    });
    if (xsGmReroll) xsGmReroll.addEventListener("click", () => {
      state.gmDrawKey = null;                       // force a fresh random draw
      renderGammaMu();
    });
    // bottom-drawer tabs (one graph per tab) + per-graph corner resize + sperm CSV
    xsTabsEl.addEventListener("click", (e) => { const t = e.target.closest(".xs-gtab"); if (t) switchXsTab(t.dataset.tab); });
    xsPanels.querySelectorAll(".xs-resizable").forEach((box) => {
      const plot = box.querySelector(".xs-plot"); if (!plot) return;
      let raf = 0;
      new ResizeObserver(() => {
        if (!box.offsetParent) return;
        cancelAnimationFrame(raf); raf = requestAnimationFrame(() => { try { Plotly.Plots.resize(plot); } catch (_) {} });
      }).observe(box);
    });
    if (xsSpDownload) xsSpDownload.addEventListener("click", downloadSpermCSV);
    if (xsSpDensity) xsSpDensity.addEventListener("change", () => { state.spDensity = xsSpDensity.checked; renderSperm(); });
    if (xsGmDensity) xsGmDensity.addEventListener("change", () => { state.gmDensity = xsGmDensity.checked; renderGammaMu(); });
    if (xsGmSperm) xsGmSperm.addEventListener("change", renderGammaMu);
  }
  function downloadSpermCSV() {
    const sp = state.spermData; if (!sp) return;
    const mode = state.crossKey;
    const cross = curAGG(), known = cross ? new Set(cross.embryos.map((e) => e.id)) : null;
    const rows = sp.embryos.filter((r) => r.modes && r.modes[mode] && (!known || known.has(r.id)));
    const lines = [`# sperm_concordance plane=${mode} (Ω=sperm on greater-transcript side, Σ=lesser)`,
      "zygote,mark,sperm_side_count,other_side_count,plane"];
    rows.forEach((r) => { const m = r.modes[mode]; lines.push([r.label.replace(/,/g, ""), m.mark, m.cntSperm, m.cntOther, m.plane].join(",")); });
    const blob = new Blob([lines.join("\n")], { type: "text/csv" });
    const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = `sperm_concordance_${mode}.csv`;
    document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(a.href);
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
