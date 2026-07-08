/* Pronuclei Distance vs Transcripts — zygotes with two pronuclei (segments 3 & 4).
 * Built on viewer-core.js (VCore). The 3-D view shows the two pronuclei and the
 * shortest line between them; the bottom drawer scatters that distance against the
 * zygote's total transcript count, with a least-squares regression. Precompute in
 * build_pronuclei.py; the UI only reads it.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const LINE_C = "#111827", CUR_C = "#0891b2";

  function plotInto(div, traces, layout, cfg) {
    if (!div.classList.contains("js-plotly-plot")) div.innerHTML = "";
    Plotly.react(div, traces, layout, cfg || { responsive: true, displayModeBar: false });
  }

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const pnReadout = $("#pn-readout"), pnFit = $("#pn-fit");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), drawerBody = $("#drawer-body");
  const scatterPlot = $("#scatter-plot");

  const state = { points: [], currentId: null, scene: null, fit: null, drawerOpen: false };

  (async function init() {
    try {
      const m = await (await fetch("data/pronuclei_manifest.json")).json();
      state.points = m.embryos;
      countEl.textContent = `${m.embryos.length} zygotes · pronuclei auto-detected inside the cytoplasm`;
      state.fit = linreg(state.points.map((p) => p.total), state.points.map((p) => p.distance));
      V.buildTabs(tabsEl, m.embryos, selectEmbryo, (e) => ({
        label: e.label, sub: e.date_short,
        title: `${e.label} · dist ${e.distance} µm · ${e.total.toLocaleString()} transcripts`,
      }));
      wireDrawer();
    } catch (err) { showError("Failed to load manifest: " + (err.message || err)); }
  })();

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    V.markActiveTab(tabsEl, id);
    const meta = state.points.find((p) => p.id === id) || {};
    showLoading(`Loading ${meta.label || id}…`);
    try {
      const scene = await V.loadGz(`data/pronuclei/${id}.json.gz`);
      if (state.currentId !== id) return;
      state.scene = scene;
      controlsEl.hidden = false; placeholder.hidden = true; drawer.hidden = false;
      render(); renderReadout(meta);
      if (!state.drawerOpen) openDrawer(true);        // reveal the scatter on first pick
      else renderScatter();
    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }

  // ---------- 3-D ----------
  function segMesh(scene, lbl, color, opacity, name) {
    const mesh = scene.region_meshes[String(lbl)]; if (!mesh) return null;
    const v = mesh.verts, f = mesh.faces, nV = v.length / 3, nF = f.length / 3;
    const x = new Array(nV), y = new Array(nV), z = new Array(nV);
    for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
    const ii = new Array(nF), jj = new Array(nF), kk = new Array(nF);
    for (let i = 0; i < nF; i++) { ii[i] = f[i * 3]; jj[i] = f[i * 3 + 1]; kk[i] = f[i * 3 + 2]; }
    return { type: "mesh3d", x, y, z, i: ii, j: jj, k: kk, color, opacity,
      name, showlegend: true, flatshading: false, hoverinfo: "name",
      lighting: { ambient: 0.7, diffuse: 0.55, specular: 0.12, roughness: 0.9 }, legendrank: lbl };
  }
  function render() {
    const s = state.scene; if (!s) return;
    const [la, lb] = s.pron_labels;                    // auto-detected pronuclei labels
    const pcolor = { [la]: "#2563eb", [lb]: "#dc2626" };
    const traces = [];
    for (const lbl of s.mask_labels) {
      const pron = lbl === la || lbl === lb;
      const name = pron ? `Pronucleus (seg ${lbl})` : `Segment ${lbl}`;
      const t = segMesh(s, lbl, pron ? pcolor[lbl] : "#9aa3b2", pron ? 0.5 : 0.08, name);
      if (t) traces.push(t);
    }
    const [a, b] = s.line_plot;
    traces.push({ type: "scatter3d", mode: "lines+markers", name: `min dist ${s.distance_um} µm`,
      x: [a[0], b[0]], y: [a[1], b[1]], z: [a[2], b[2]],
      line: { color: LINE_C, width: 7 }, marker: { size: 4, color: LINE_C },
      hovertemplate: `min distance ${s.distance_um} µm<extra></extra>`, legendrank: 100 });
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }
  function renderReadout(meta) {
    const s = state.scene, pred = state.fit.a + state.fit.b * s.total_transcripts;
    pnReadout.innerHTML =
      `<div class="pn-big"><span>${s.distance_um}</span> µm <span class="pn-lbl">pronuclei distance</span></div>` +
      `<div class="pn-big"><span>${s.total_transcripts.toLocaleString()}</span> <span class="pn-lbl">total transcripts</span></div>` +
      `<div class="pn-resid">pronuclei auto-detected as segments <b>${s.pron_labels[0]}</b> &amp; <b>${s.pron_labels[1]}</b></div>` +
      `<div class="pn-resid">fit predicts ${pred.toFixed(1)} µm here (residual ${(s.distance_um - pred >= 0 ? "+" : "") + (s.distance_um - pred).toFixed(1)} µm)</div>`;
  }

  // ---------- bottom-drawer scatter ----------
  function linreg(xs, ys) {
    const n = xs.length; let sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
    for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; syy += ys[i] * ys[i]; }
    const den = n * sxx - sx * sx;
    const b = den ? (n * sxy - sx * sy) / den : 0;
    const a = (sy - b * sx) / n;
    const rden = Math.sqrt(den * (n * syy - sy * sy));
    const r = rden ? (n * sxy - sx * sy) / rden : 0;
    return { a, b, r, r2: r * r, n };
  }
  function renderScatter() {
    const pts = state.points, f = state.fit;
    const xs = pts.map((p) => p.total), ys = pts.map((p) => p.distance);
    const cur = state.currentId;
    const others = pts.filter((p) => p.id !== cur), curP = pts.find((p) => p.id === cur);
    const xmin = Math.min(...xs), xmax = Math.max(...xs);
    const traces = [
      { type: "scatter", mode: "markers", name: "zygotes",
        x: others.map((p) => p.total), y: others.map((p) => p.distance),
        marker: { size: 8, color: "#94a3b8", opacity: 0.75, line: { width: 0 } },
        text: others.map((p) => p.label), hovertemplate: "%{text}<br>%{x:,} transcripts<br>%{y} µm<extra></extra>" },
      { type: "scatter", mode: "lines", name: "least-squares fit",
        x: [xmin, xmax], y: [f.a + f.b * xmin, f.a + f.b * xmax],
        line: { color: "#0891b2", width: 2, dash: "solid" }, hoverinfo: "skip" },
    ];
    if (curP) traces.push({ type: "scatter", mode: "markers", name: curP.label,
      x: [curP.total], y: [curP.distance],
      marker: { size: 15, color: CUR_C, line: { width: 2, color: "#fff" }, symbol: "circle" },
      hovertemplate: `${curP.label}<br>%{x:,} transcripts<br>%{y} µm<extra></extra>` });
    plotInto(scatterPlot, traces, {
      margin: { l: 52, r: 12, t: 8, b: 44 }, height: scatterPlot.clientHeight || 300,
      xaxis: { title: { text: "total transcripts", font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", zeroline: false },
      yaxis: { title: { text: "min pronuclei distance (µm)", font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", rangemode: "tozero" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      legend: { orientation: "h", font: { size: 10 }, y: 1.12, x: 1, xanchor: "right" },
    });
    const slopePer100k = (f.b * 1e5).toFixed(1);
    pnFit.innerHTML = `<b>${f.n}</b> zygotes · Pearson r = <b>${f.r.toFixed(3)}</b> ` +
      `(r² = ${f.r2.toFixed(3)}) · slope ${slopePer100k} µm per 100k transcripts`;
  }

  // ---------- drawer ----------
  function openDrawer(open) {
    state.drawerOpen = open;
    drawer.dataset.open = open ? "true" : "false";
    drawerHandle.setAttribute("aria-expanded", String(open));
    if (open) { renderScatter(); requestAnimationFrame(() => { try { Plotly.Plots.resize(scatterPlot); } catch (_) {} }); }
  }
  function wireDrawer() {
    wireHandleDrag(drawer, drawerHandle, {
      computeSize: (e) => window.innerHeight - e.clientY - 40,
      clampSize: (px) => Math.max(180, Math.min(window.innerHeight - 100, px)),
      applySize: (px) => drawer.style.setProperty("--drawer-h", px + "px"),
      setOpen: openDrawer,
      afterDrag: () => { try { Plotly.Plots.resize(scatterPlot); } catch (_) {} },
    });
    let sh = 0; const rz = $("#drawer-resize");
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return;
      const h = Math.max(180, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y)));
      drawer.style.setProperty("--drawer-h", h + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} try { Plotly.Plots.resize(scatterPlot); } catch (_) {} } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }
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
      cfg.applySize(cfg.clampSize(cfg.computeSize(e))); e.preventDefault();
    });
    const up = (e) => {
      if (!start) return;
      try { handleEl.releasePointerCapture(e.pointerId); } catch (_) {}
      if (moved) { drawerEl.classList.remove("dragging"); cfg.afterDrag && cfg.afterDrag(); }
      else cfg.setOpen(drawerEl.dataset.open !== "true");
      start = null; moved = false;
    };
    handleEl.addEventListener("pointerup", up); handleEl.addEventListener("pointercancel", up);
  }

  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false;
    placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
