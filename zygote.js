/* Zygote Division Planes — analysis model.
 * Built on viewer-core.js (VCore): 60 zygotes, the polar-body axis, 17 candidate
 * division planes, per-gene transcript split (blue/red) across the selected plane,
 * a counts+null chart, best-plane tables, and a cross-section plot. All statistics
 * are pre-computed (build_zygote.py); the UI only reads them.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const AXIS_C = "#111827", PLANE_C = "#f97316", BEST_C = "#16a34a";
  const BLUE = "#2563eb", RED = "#dc2626";

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const geneSelect = $("#gene-select"), planeSelect = $("#plane-select");
  const axisShow = $("#axis-show"), planeShow = $("#plane-show"), bestShow = $("#best-show");
  const chartEl = $("#chart"), chartSub = $("#chart-sub"), chartReadout = $("#chart-readout");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), drawerBody = $("#drawer-body");
  const drawerEmb = $("#drawer-emb"), crossPlot = $("#cross-plot"), crossPmode = $("#cross-pmode");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle");
  const rtabsEl = $("#rtabs"), bestListEl = $("#best-list");

  const state = {
    manifest: [], currentId: null, scene: null, userGene: null, planeIdx: 0,
    drawerOpen: false, bestTab: "pVol", crossMode: "vol",
  };

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
      drawerEmb.textContent = "· " + state.byLabel(id);
      render(); renderChart(); renderBestList(); if (state.drawerOpen) renderCross();
    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }
  state.byLabel = (id) => (state.manifest.find((m) => m.id === id) || {}).label || id;

  function populateGenes(scene) {
    const tot = scene.gene_totals || {};
    geneSelect.innerHTML = scene.genes
      .map((g) => `<option value="${g}">${g}  (${(tot[g] || 0).toLocaleString()})</option>`).join("");
    geneSelect.value = (state.userGene && scene.genes.includes(state.userGene)) ? state.userGene : scene.genes[0];
  }

  // ---------- geometry helpers ----------
  const gene = () => geneSelect.value;
  function planeGeo(k) { return state.scene.analysis.planes[k]; }
  // side of a transcript for plane k: (pos_um − com_um)·normal_um  (matches precompute)
  function splitCloud(scene, g, k) {
    const t = scene.transcripts[g]; const A = scene.analysis;
    const com = A.com_um, nrm = planeGeo(k).normal_um, zs = scene.z_scale;
    const bx = [], by = [], bz = [], rx = [], ry = [], rz = [];
    if (t) for (let i = 0; i < t.x.length; i++) {
      const s = (t.x[i] * XY - com[0]) * nrm[0] + (t.y[i] * XY - com[1]) * nrm[1] + (t.gz[i] - com[2]) * nrm[2];
      if (s > 0) { bx.push(t.x[i]); by.push(t.y[i]); bz.push(t.gz[i] * zs); }
      else { rx.push(t.x[i]); ry.push(t.y[i]); rz.push(t.gz[i] * zs); }
    }
    return { bx, by, bz, rx, ry, rz };
  }
  // plane k as an orange quad (physical square in µm → plot space)
  const PLANE_SCALE = 2;                              // rendered plane size (×) vs precomputed L
  function planeQuad(k, color, op, name, rank) {
    const A = state.scene.analysis, p = planeGeo(k), zs = state.scene.z_scale;
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
      i: [0, 0], j: [1, 2], k: [2, 3], color, opacity: op, name, showlegend: true,
      hoverinfo: "name", flatshading: true, legendrank: rank,
    };
  }

  // ---------- 3-D render ----------
  function render() {
    const s = state.scene; if (!s) return;
    const zs = s.z_scale, A = s.analysis, k = state.planeIdx, g = gene();
    const traces = V.bodyTraces(s);

    const sp = splitCloud(s, g, k);
    traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side A`,
      x: sp.bx, y: sp.by, z: sp.bz, marker: { size: 2.6, color: BLUE, opacity: 0.85, line: { width: 0 } },
      hovertemplate: `${g} · side A<extra></extra>`, legendrank: 20000 });
    traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side B`,
      x: sp.rx, y: sp.ry, z: sp.rz, marker: { size: 2.6, color: RED, opacity: 0.85, line: { width: 0 } },
      hovertemplate: `${g} · side B<extra></extra>`, legendrank: 20001 });

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
    if (bestShow.checked) {
      const bp = A.best_planes;
      const items = [["pVol", "min p·vol"], ["pCnt", "min p·cnt"], ["diffVol", "max Δ·vol"], ["diffCnt", "max Δ·cnt"]];
      const seen = new Set(); let rank = 42000;
      for (const [key] of items) {
        const kk = bp[key];
        if (seen.has(kk)) continue;                 // dedup coincident best planes
        seen.add(kk);
        const which = items.filter(([k2]) => bp[k2] === kk).map(([, l]) => l).join(", ");
        traces.push(planeQuad(kk, BEST_C, 0.22, `Best ${kk * state.step}° (${which})`, rank++));
      }
    }
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  // ---------- counts chart (floating window) ----------
  function geneRow(g, k) {
    const rows = state.scene.analysis.genes;
    const r = rows.find((x) => x.gene === g);
    return r ? r.planes[k] : null;
  }
  function renderChart() {
    const s = state.scene; if (!s) return;
    const g = gene(), k = state.planeIdx, row = geneRow(g, k);
    chartSub.textContent = `· ${g} · plane ${k * state.step}°`;
    if (!row) { Plotly.purge(chartEl); chartEl.innerHTML = '<div class="chart-readout">Gene not in this embryo.</div>'; chartReadout.innerHTML = ""; return; }
    chartEl.innerHTML = "";
    const traces = [
      { type: "bar", name: "Real", x: ["Side A", "Side B"], y: [row.a, row.b],
        marker: { color: [BLUE, RED] }, hovertemplate: "%{x}: %{y}<extra>real</extra>" },
      { type: "bar", name: "Null", x: ["Side A", "Side B"], y: [row.na, row.nb],
        marker: { color: "#9ca3af" }, opacity: 0.7, hovertemplate: "%{x}: %{y}<extra>null</extra>" },
    ];
    Plotly.react(chartEl, traces, {
      barmode: "group", margin: { l: 34, r: 6, t: 6, b: 20 }, height: 150,
      yaxis: { tickfont: { size: 9 }, gridcolor: "#eef1f5", fixedrange: true, title: { text: "count", font: { size: 9 } } },
      xaxis: { tickfont: { size: 11 }, fixedrange: true }, bargap: 0.3, bargroupgap: 0.08,
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", font: { color: "#1a2233" },
      legend: { orientation: "h", font: { size: 10 }, y: 1.15, x: 1, xanchor: "right" },
    }, { responsive: true, displayModeBar: false });
    const sigV = row.pVol <= 0.05, sigC = row.pCnt <= 0.05;
    chartReadout.innerHTML =
      `<div>n = <b>${row.a + row.b}</b> · Δcount = <b>${row.dCount}</b> ` +
      `(real) vs <b>${row.ndCount}</b> (null)</div>` +
      `<div>Δ / total = <b>${row.dNorm.toFixed(3)}</b> · Δ / volume = <b>${row.dVol.toExponential(2)}</b></div>` +
      `<div>p(vol) = <span class="${sigV ? "sig" : ""}">${fmtP(row.pVol)}</span> · ` +
      `p(count) = <span class="${sigC ? "sig" : ""}">${fmtP(row.pCnt)}</span></div>`;
  }
  const fmtP = (p) => p < 0.001 ? p.toExponential(1) : p.toFixed(3);

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
    const A = s.analysis, meta = BEST_META[state.bestTab];
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
  function renderCross() {
    const s = state.scene; if (!s) return;
    const A = s.analysis, cs = A.cross_section, planes = A.planes;
    const outline = cs.outline;
    if (!outline.length) { Plotly.purge(crossPlot); crossPlot.innerHTML = '<div class="chart-readout" style="padding:20px;text-align:center">No cross-section.</div>'; return; }
    crossPlot.innerHTML = "";
    let R = 0; for (const p of outline) R = Math.max(R, Math.hypot(p[0], p[1]));
    // normalize p across the planes so the colormap spans the full range (low p = intense)
    const pOf = (k) => state.crossMode === "vol" ? planes[k].wpVol : planes[k].wpCnt;
    const ps = planes.map((_, k) => pOf(k));
    const pmin = Math.min(...ps), pmax = Math.max(...ps), span = pmax - pmin || 1;
    const traces = [];
    for (let k = 0; k < planes.length; k++) {
      const th = k * state.step * Math.PI / 180;
      const dir = [-Math.sin(th), Math.cos(th)];   // plane ∩ cross-section direction (u,v)
      const p = pOf(k), t = (pmax - p) / span;     // 1 = lowest p = most intense
      const sel = k === state.planeIdx;
      traces.push({ type: "scatter", mode: "lines", showlegend: false,
        x: [-R * 1.05 * dir[0], R * 1.05 * dir[0]], y: [-R * 1.05 * dir[1], R * 1.05 * dir[1]],
        line: { color: sel ? PLANE_C : viridis(t), width: sel ? 4 : (2 + 2.8 * t) },
        hovertemplate: `plane ${k * state.step}°<br>weighted p = ${fmtP(p)}<extra></extra>` });
    }
    // outline (bold black, closed)
    traces.push({ type: "scatter", mode: "lines", showlegend: false,
      x: outline.map((p) => p[0]).concat([outline[0][0]]),
      y: outline.map((p) => p[1]).concat([outline[0][1]]),
      line: { color: "#0b0d13", width: 3 }, hoverinfo: "skip" });
    const lim = R * 1.12;
    Plotly.react(crossPlot, traces, {
      margin: { l: 10, r: 10, t: 8, b: 10 }, height: crossPlot.clientHeight || 280,
      xaxis: { range: [-lim, lim], visible: false, fixedrange: true, scaleanchor: "y", scaleratio: 1 },
      yaxis: { range: [-lim, lim], visible: false, fixedrange: true },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, { responsive: true, displayModeBar: false });
  }

  // ---------- wiring ----------
  function wireControls() {
    geneSelect.addEventListener("change", () => { state.userGene = geneSelect.value; render(); renderChart(); highlightBest(); });
    planeSelect.addEventListener("change", () => { state.planeIdx = parseInt(planeSelect.value, 10) || 0; render(); renderChart(); if (state.drawerOpen) renderCross(); });
    [axisShow, planeShow, bestShow].forEach((c) => c.addEventListener("change", () => render()));
  }
  function highlightBest() {
    const cur = gene();
    bestListEl.querySelectorAll(".best-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === cur));
  }
  function wireRdrawer() {
    rdrawerHandle.addEventListener("click", () => {
      const open = rdrawer.dataset.open !== "true";
      rdrawer.dataset.open = open ? "true" : "false";
      rdrawerHandle.setAttribute("aria-expanded", String(open));
    });
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
  function wireDrawer() {
    drawerHandle.addEventListener("click", () => {
      state.drawerOpen = drawer.dataset.open !== "true";
      drawer.dataset.open = state.drawerOpen ? "true" : "false";
      drawerHandle.setAttribute("aria-expanded", String(state.drawerOpen));
      if (state.drawerOpen) { renderCross(); requestAnimationFrame(() => { try { Plotly.Plots.resize(crossPlot); } catch (_) {} }); }
    });
    // height resize (reuse the shell's drawer-resize)
    let sh = 0;
    const rz = $("#drawer-resize");
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { x: e.clientX, y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return;
      const h = Math.max(160, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y)));
      drawer.style.setProperty("--drawer-h", h + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} try { Plotly.Plots.resize(crossPlot); } catch (_) {} } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
    crossPmode.addEventListener("change", () => { state.crossMode = crossPmode.value; renderCross(); });
  }

  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false;
    placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
