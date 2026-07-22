/* Gene Diffusion Rates — a passive-diffusion null model. mRNA is born at the nucleus (the two
 * pronuclei) and spreads by Brownian motion until it matches the gene's observed distribution;
 * the time that takes (set by the mRNA's size → diffusion coefficient) is the read-out. The 3-D
 * scene animates one precomputed diffusion trajectory (data/gene_diffusion.json.gz), stopping at
 * the frame that matches the observed spread. Right drawer = per-gene times (longest first) for
 * the embryo; bottom drawer = a cross-embryo scatter of mean diffusion time vs a chosen property. */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const XY = 0.15, ACCENT = "#0284c7";
  const CELL_C = "#c9d3df", PRON_C = "#f59e0b", PART_C = "#0284c7", OBS_C = "#94a3b8";
  const FPS = 15;   // deliberate playback so the fine-grained diffusion is watchable

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const geneSelect = $("#gene-select"), metricSelect = $("#metric-select");
  const dfReadout = $("#df-readout"), dfNote = $("#df-note");
  const playBtn = $("#df-play"), stopBtn = $("#df-stop"), resetBtn = $("#df-reset");
  const progFill = $("#df-progress-fill"), icSeg = $("#df-ic");
  const dfTimeplot = $("#df-timeplot"), dfTimeplotSub = $("#df-timeplot-sub");
  const obsToggle = $("#df-observed"), cellToggle = $("#df-cell");
  const cellColorInp = $("#df-cell-color"), cellOpacityInp = $("#df-cell-opacity"), cellOpacityVal = $("#df-cell-opacity-val");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), dfScatter = $("#df-scatter");
  const dfYaxis = $("#df-yaxis"), dfScatterNote = $("#df-scatter-note"), dfScatterSub = $("#df-scatter-sub");
  const dfModel = $("#df-model"), dfFitStats = $("#df-fit-stats");
  const drawerBody = $("#drawer-body"), dfDrawerTabs = $("#df-drawer-tabs");
  const dfAgeSrc = $("#df-age-src"), dfAgeNote = $("#df-age-note"), dfAgeLegend = $("#df-age-legend"), dfAgeGrid = $("#df-age-grid");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle"), dfList = $("#df-list"), dfEmbAvg = $("#df-embavg");

  const state = {
    data: null, genesD: null, meta: null, points: [], byId: {}, currentId: null,
    scene: null, rec: null, ic: "com", metric: "ks", gene: null,
    playing: false, frame: 0, stopFrame: 0, sampleIdx: null, raf: 0, lastT: 0,
    dotSize: 2.5, drawerOpen: false, drawerTab: "scatter", cellColor: "#c9d3df", cellOpacity: 0.16,
  };
  let vcExtras = null;

  const icRec = () => state.rec && state.rec.ic[state.ic];
  const geneRec = (g) => { const r = icRec(); return r && r.genes[g]; };
  const D_of = (g) => (state.genesD[g] || {}).D || null;
  // real time (seconds) for a gene to reach its observed distribution: t = τ_match / (2·D)
  const timeSec = (g, ic, metric) => {
    const r = state.rec && state.rec.ic[ic]; const gr = r && r.genes[g]; const D = D_of(g);
    if (!gr || !D) return null;
    return gr.tau[metric] / (2 * D);
  };
  function fmtTime(sec) {
    if (sec == null || !isFinite(sec)) return { v: "—", u: "" };
    const min = sec / 60;
    if (min < 90) return { v: min.toFixed(min < 10 ? 1 : 0), u: "min" };
    return { v: (min / 60).toFixed(1), u: "hours" };
  }

  (async function init() {
    try {
      const [man, data] = await Promise.all([
        (await fetch("data/pronuclei_manifest.json")).json(),
        V.loadGz("data/gene_diffusion.json.gz"),
      ]);
      state.data = data.embryos; state.genesD = data.genes; state.meta = data.meta;
      state.points = man.embryos.filter((p) => data.embryos[p.id]);
      state.points.forEach((p) => (state.byId[p.id] = p));
      countEl.textContent = `${state.points.length} zygotes · passive-diffusion null model`;
      V.buildTabs(tabsEl, state.points, selectEmbryo, (e) => ({
        label: e.label, sub: e.date_short,
        title: `${e.label} · mean diffusion time ${fmtEmbAvg(e.id)}`,
      }));
      wireControls(); wireDrawer(); wireRdrawer();
      vcExtras = V.addWindowExtras($("#controls-body"), { defaultSize: state.dotSize, sizeLabel: "Particle size",
        onDotSize: (s) => { state.dotSize = s; if (state.scene) drawFrame(); } });
    } catch (err) { showError("Failed to load: " + (err.message || err)); }
  })();

  function fmtEmbAvg(id) {
    const a = embAvgSec(id, state.ic, state.metric); const f = fmtTime(a); return a == null ? "—" : `${f.v} ${f.u}`;
  }
  function embAvgSec(id, ic, metric) {
    const r = state.data[id] && state.data[id].ic[ic]; if (!r) return null;
    let s = 0, n = 0;
    for (const g in r.genes) { const D = (state.genesD[g] || {}).D; if (!D) continue; s += r.genes[g].tau[metric] / (2 * D); n++; }
    return n ? s / n : null;
  }

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    stopAnim();
    state.currentId = id; V.markActiveTab(tabsEl, id);
    showLoading(`Loading ${(state.byId[id] || {}).label || id}…`);
    try {
      state.rec = state.data[id];
      const jobs = [V.loadGz(`data/pronuclei/${id}.json.gz`)];
      if (!state.rec.ic.com.traj) jobs.push(V.loadGz(`data/gene_diffusion/${id}.json.gz`));   // lazy animation trajectories
      const [scene, traj] = await Promise.all(jobs);
      if (state.currentId !== id) return;
      if (traj) for (const ic in traj.ic) if (state.rec.ic[ic]) state.rec.ic[ic].traj = traj.ic[ic].traj;
      state.scene = scene;
      populateGenes();
      if (vcExtras) vcExtras.setAtlas(id);
      controlsEl.hidden = false; placeholder.hidden = true; drawer.hidden = false; rdrawer.hidden = false;
      state.frame = 0; pickSample(); computeStop();
      drawFrame(); renderReadout(); renderList(); renderTimePlot();
      renderActiveDrawer();
    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }

  function populateGenes() {
    const r = icRec(); const genes = r ? Object.keys(r.genes) : [];
    genes.sort((a, b) => (timeSec(b, state.ic, state.metric) || 0) - (timeSec(a, state.ic, state.metric) || 0));
    geneSelect.innerHTML = genes.map((g) => {
      const D = state.genesD[g] || {}; const t = fmtTime(timeSec(g, state.ic, state.metric));
      return `<option value="${g}">${g}  (${t.v} ${t.u})</option>`;
    }).join("");
    if (!genes.includes(state.gene)) state.gene = genes[0];
    geneSelect.value = state.gene;
  }

  // stable particle subsample sized to the gene's transcript count (≤ stored trajectory width)
  function pickSample() {
    const r = icRec(); if (!r) { state.sampleIdx = []; return; }
    const P = r.traj[0].length; const gr = geneRec(state.gene);
    const k = Math.max(1, Math.min(gr ? gr.n : P, P));
    state.sampleIdx = Array.from({ length: k }, (_, i) => Math.floor(i * P / k));
  }
  function computeStop() {
    const r = icRec(); const gr = geneRec(state.gene); if (!r || !gr) { state.stopFrame = 0; return; }
    const tauM = gr.tau[state.metric]; const tf = r.tau_frames;
    let f = tf.findIndex((t) => t >= tauM); if (f < 0) f = tf.length - 1;
    state.stopFrame = f;
  }

  // ---------- 3-D ----------
  function segMesh(lbl, color, opacity, name, rank) {
    const mesh = state.scene.region_meshes[String(lbl)]; if (!mesh) return null;
    const v = mesh.verts, f = mesh.faces, nV = v.length / 3, nF = f.length / 3;
    const x = new Array(nV), y = new Array(nV), z = new Array(nV);
    for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
    const ii = new Array(nF), jj = new Array(nF), kk = new Array(nF);
    for (let i = 0; i < nF; i++) { ii[i] = f[i * 3]; jj[i] = f[i * 3 + 1]; kk[i] = f[i * 3 + 2]; }
    return { type: "mesh3d", x, y, z, i: ii, j: jj, k: kk, color, opacity, name, showlegend: true,
      hoverinfo: "name", flatshading: false, lighting: { ambient: 0.75, diffuse: 0.5, specular: 0.1, roughness: 0.9 }, legendrank: rank };
  }
  function staticTraces() {
    const s = state.scene, r = state.rec, traces = [];
    if (cellToggle.checked) {
      const t1 = segMesh(1, state.cellColor, state.cellOpacity, "Cell (cytoplasm)", 10); if (t1) traces.push(t1);
      const pl = (s.pron_labels || []);
      pl.forEach((lbl, i) => { const t = segMesh(lbl, PRON_C, 0.28, `Nucleus · pronucleus ${i + 1}`, 20 + i); if (t) traces.push(t); });
    }
    if (obsToggle.checked && s.transcripts && s.transcripts[state.gene]) {
      const t = s.transcripts[state.gene], zs = s.z_scale;
      traces.push({ type: "scatter3d", mode: "markers", name: `${state.gene} · observed`,
        x: t.x, y: t.y, z: t.gz.map((g) => g * zs),
        marker: { size: Math.max(1.5, state.dotSize * 0.7), color: OBS_C, opacity: 0.5 },
        hovertemplate: `${state.gene} observed<extra></extra>`, legendrank: 30 });
    }
    return traces;
  }
  function particleTrace() {
    const r = icRec(); if (!r) return null;
    const pts = r.traj[state.frame], idx = state.sampleIdx;
    const x = idx.map((i) => pts[i][0]), y = idx.map((i) => pts[i][1]), z = idx.map((i) => pts[i][2]);
    return { type: "scatter3d", mode: "markers", name: `${state.gene} · simulated (${idx.length})`,
      x, y, z, marker: { size: state.dotSize, color: PART_C, opacity: 0.9, line: { width: 0 } },
      hovertemplate: `${state.gene} diffusing<extra></extra>`, legendrank: 5 };
  }
  function drawFrame() {
    const s = state.scene; if (!s) return;
    const traces = staticTraces(); const p = particleTrace(); if (p) traces.push(p);
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
    const r = icRec(); progFill.style.width = r ? `${100 * state.frame / (r.traj.length - 1)}%` : "0%";
  }
  // fast per-frame marker update during playback (restyle only the particle trace)
  function updateParticles() {
    const r = icRec(); if (!r) return;
    const pts = r.traj[state.frame], idx = state.sampleIdx, pi = plotHost.data.length - 1;
    Plotly.restyle(plotHost, {
      x: [idx.map((i) => pts[i][0])], y: [idx.map((i) => pts[i][1])], z: [idx.map((i) => pts[i][2])],
    }, [pi]);
    progFill.style.width = `${100 * state.frame / (r.traj.length - 1)}%`;
  }

  // ---------- animation ----------
  function play() {
    if (state.playing || !icRec()) return;
    if (state.frame >= state.stopFrame) state.frame = 0;
    state.playing = true; playBtn.classList.add("playing"); playBtn.textContent = "❚❚ Playing";
    drawFrame(); updateTimeMarker(); state.lastT = 0; state.raf = requestAnimationFrame(tick);
  }
  function tick(ts) {
    if (!state.playing) return;
    if (ts - state.lastT >= 1000 / FPS) {
      state.lastT = ts;
      if (state.frame >= state.stopFrame) { finishAnim(); return; }
      state.frame++; updateParticles(); updateTimeMarker();
    }
    state.raf = requestAnimationFrame(tick);
  }
  function stopAnim() {
    state.playing = false; if (state.raf) cancelAnimationFrame(state.raf); state.raf = 0;
    playBtn.classList.remove("playing"); playBtn.textContent = "▶ Start";
  }
  function finishAnim() { stopAnim(); renderReadout(true); }

  // ---------- read-out ----------
  function renderReadout(done) {
    const g = state.gene, gr = geneRec(g), D = D_of(g); if (!gr) { dfReadout.innerHTML = ""; return; }
    const t = fmtTime(timeSec(g, state.ic, state.metric));
    const gd = state.genesD[g] || {}; const mrna = gd.mrna_nt;
    dfReadout.innerHTML =
      `<div class="df-big"><b>${t.v}</b> ${t.u} <span class="df-muted">· diffusion time${done ? " ✓ matched" : ""}</span></div>` +
      `<div class="df-sub"><b>${g}</b> · ${gr.n.toLocaleString()} transcripts · observed mean reach <b>${gr.obs_mean} µm</b></div>` +
      `<div class="df-sub df-muted">mRNA ${mrna ? mrna.toLocaleString() + " nt" : "length n/a (panel median)"} · D = ${D ? D.toFixed(4) : "—"} µm²/s</div>`;
  }

  // ---------- distance-to-nucleus over time (floating-window mini plot) ----------
  // The cloud's distance-to-nucleus distribution (10–90% band + mean) as it evolves, on the
  // gene's real-time axis (τ/2D). Target = observed mean; green dot = the match; the dark
  // vertical line tracks where the 3-D animation currently is.
  function comOriginUm() {
    const zs = state.scene.z_scale, op = (state.rec.ic.com || {}).origin_plot;
    return op ? [op[0] * XY, op[1] * XY, op[2] / zs] : [0, 0, 0];
  }
  function simBands(ic) {
    const r = state.rec && state.rec.ic[ic]; if (!r || !r.traj) return null;
    if (r._bands) return r._bands;
    const zs = state.scene.z_scale, o = comOriginUm(), p10 = [], p50 = [], p90 = [];
    for (let f = 0; f < r.traj.length; f++) {
      const d = r.traj[f].map((p) => Math.hypot(p[0] * XY - o[0], p[1] * XY - o[1], p[2] / zs - o[2])).sort((a, b) => a - b);
      const q = (t) => d[Math.min(d.length - 1, Math.floor(t * d.length))];
      p10.push(q(0.1)); p90.push(q(0.9)); p50.push(d.reduce((a, b) => a + b, 0) / d.length);
    }
    r._bands = { p10, p50, p90 }; return r._bands;
  }
  function renderTimePlot() {
    if (!state.rec || !state.scene || !dfTimeplot) return;
    const b = simBands(state.ic), gr = geneRec(state.gene), D = D_of(state.gene);
    if (!b || !gr || !D) { try { Plotly.purge(dfTimeplot); } catch (_) {} return; }
    const tf = state.rec.ic[state.ic].tau_frames;
    const secMax = tf[tf.length - 1] / (2 * D);
    const u = secMax / 60 < 90 ? { div: 60, lab: "minutes" } : { div: 3600, lab: "hours" };
    dfTimeplot._unit = u.div;
    const x = tf.map((t) => t / (2 * D) / u.div), stopF = state.stopFrame, obs = gr.obs_mean;
    const ymax = Math.max(Math.max(...b.p90), obs) * 1.05, C = "#0284c7";
    dfTimeplotSub.textContent = `· ${state.gene}`;
    const traces = [
      { x, y: b.p10, mode: "lines", line: { width: 0 }, hoverinfo: "skip", showlegend: false },
      { x, y: b.p90, mode: "lines", line: { width: 0 }, fill: "tonexty", fillcolor: "rgba(2,132,199,.13)", hoverinfo: "skip", showlegend: false },
      { x, y: b.p50, mode: "lines", line: { color: C, width: 2 }, hovertemplate: "%{x:.1f} · mean %{y:.0f} µm<extra></extra>", showlegend: false },
      { x: [0, x[x.length - 1]], y: [obs, obs], mode: "lines", line: { color: "#64748b", width: 1.4, dash: "dash" }, hoverinfo: "skip", showlegend: false },
      { x: [x[stopF]], y: [b.p50[stopF]], mode: "markers", marker: { size: 8, color: "#16a34a", line: { width: 1.5, color: "#fff" } }, hovertemplate: `matched · %{x:.1f} ${u.lab}<extra></extra>`, showlegend: false },
      { x: [x[state.frame], x[state.frame]], y: [0, ymax], mode: "lines", line: { color: "#0f172a", width: 1.4 }, hoverinfo: "skip", showlegend: false },
    ];
    Plotly.react(dfTimeplot, traces, {
      margin: { l: 36, r: 8, t: 4, b: 28 }, height: 150,
      xaxis: { title: { text: u.lab, font: { size: 9 } }, tickfont: { size: 8.5 }, gridcolor: "#eef1f5", zeroline: false, fixedrange: true },
      yaxis: { title: { text: "µm", font: { size: 9 } }, tickfont: { size: 8.5 }, gridcolor: "#eef1f5", zeroline: false, rangemode: "tozero", fixedrange: true },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", showlegend: false, font: { color: "#1a2233" },
    }, { displayModeBar: false, responsive: true });
  }
  function updateTimeMarker() {
    const r = state.rec && state.rec.ic[state.ic];
    if (!r || !dfTimeplot || !dfTimeplot.data || !dfTimeplot._unit) return;
    const D = D_of(state.gene); if (!D) return;
    const t = r.tau_frames[state.frame] / (2 * D) / dfTimeplot._unit;
    try { Plotly.restyle(dfTimeplot, { x: [[t, t]] }, [dfTimeplot.data.length - 1]); } catch (_) {}
  }

  // ---------- right drawer: per-gene times ----------
  function renderList() {
    const r = icRec(); if (!r) { dfList.innerHTML = ""; return; }
    const rows = Object.keys(r.genes).map((g) => ({ g, sec: timeSec(g, state.ic, state.metric), n: r.genes[g].n }))
      .filter((x) => x.sec != null).sort((a, b) => b.sec - a.sec);
    const avg = embAvgSec(state.currentId, state.ic, state.metric); const fa = fmtTime(avg);
    dfEmbAvg.innerHTML = `<div class="df-embavg-l">embryo average</div><div class="df-embavg-t">${fa.v} <span style="font-size:14px">${fa.u}</span></div>`;
    const max = rows.length ? rows[0].sec : 1;
    dfList.innerHTML = rows.map((x) => {
      const f = fmtTime(x.sec);
      return `<div class="df-row${x.g === state.gene ? " current" : ""}" data-gene="${x.g}">` +
        `<span class="df-g">${x.g}<small>${x.n.toLocaleString()} tx</small></span>` +
        `<span class="df-t">${f.v} ${f.u}</span>` +
        `<span class="df-bar" style="width:${Math.max(4, 100 * x.sec / max)}%"></span></div>`;
    }).join("");
  }

  // ---------- bottom drawer: cross-embryo scatter ----------
  const Y_LABEL = { min_dist: "Minimum pronuclei distance (µm)", total_tx: "Total transcript number", gene_count: "Selected gene's transcript count" };
  function yVal(id) {
    const e = state.data[id]; const which = dfYaxis.value;
    if (which === "min_dist") return e.min_dist;
    if (which === "total_tx") return e.total_tx;
    return (e.gene_counts || {})[state.gene];   // selected-gene count
  }
  function renderScatter() {
    if (!dfScatter.offsetParent && dfScatter.dataset.drawn) return;
    const ids = Object.keys(state.data);
    const xs = [], ys = [], txt = [], cols = [], sz = [];
    ids.forEach((id) => {
      const x = embAvgSec(id, state.ic, state.metric); const y = yVal(id);
      if (x == null || y == null) return;
      xs.push(x / 60); ys.push(y); txt.push(state.data[id].label);
      cols.push(id === state.currentId ? ACCENT : "#b8c4d4"); sz.push(id === state.currentId ? 13 : 8);
    });
    const which = dfYaxis.value;
    dfScatterSub.textContent = `· ${which === "gene_count" ? state.gene + " count" : Y_LABEL[which]}`;
    dfScatterNote.textContent = `${xs.length} zygotes · x = mean gene diffusion time (min) · ${state.ic === "com" ? "from nucleus centre" : "from nucleus surface"} · match by ${state.metric}`;
    const traces = [{ type: "scatter", mode: "markers", name: "zygotes", x: xs, y: ys, text: txt,
      marker: { size: sz, color: cols, line: { width: 1, color: "#fff" } },
      hovertemplate: "%{text}<br>mean time %{x:.0f} min<br>" + Y_LABEL[which] + " %{y}<extra></extra>" }];
    // fit the chosen regression model (same family as the Pronuclei project) and draw the curve
    const R = window.Regressions;
    if (R && xs.length >= 3) {
      const fit = R.fitModel(dfModel.value, xs, ys);
      const lo = Math.min(...xs), hi = Math.max(...xs), cx = [], cy = [];
      for (let i = 0; i < 90; i++) { const x = lo + (hi - lo) * i / 89, p = fit.predict(x); if (isFinite(p)) { cx.push(x); cy.push(p); } }
      traces.push({ type: "scatter", mode: "lines", name: fit.label, x: cx, y: cy,
        line: { color: ACCENT, width: 2.4, shape: "spline" }, hoverinfo: "skip" });
      dfFitStats.innerHTML = R.statsHtml(xs, ys, fit, { xName: "mean diffusion time", yName: Y_LABEL[which] });
    } else if (dfFitStats) dfFitStats.textContent = "";
    const layout = { margin: { l: 60, r: 12, t: 8, b: 44 }, paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      font: { color: "#1a2233", size: 12 }, hovermode: "closest", showlegend: false,
      xaxis: { title: "mean diffusion time (min)", gridcolor: "#eef1f5", zeroline: false },
      yaxis: { title: Y_LABEL[which], gridcolor: "#eef1f5", zeroline: false } };
    Plotly.react(dfScatter, traces, layout, { responsive: true, displaylogo: false, displayModeBar: "hover",
      modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "toggleSpikelines"] });
    dfScatter.dataset.drawn = "1";
  }

  // ---------- bottom drawer: gene-time vs embryo-age grid ----------
  // Each gene gets a tile coloured by ratio r = (its diffusion time) / (the embryo's age). The age is
  // either the embryo's mean gene diffusion time, or its pronuclei-distance pseudotime (smaller distance
  // = older) rescaled onto the cohort's mean-diffusion-time range so the two live on the same clock.
  // Per the request, intensity grows as the ratio drops below 1 (diffusion time well under the age →
  // flagged as a candidate for active/transport-assisted spreading); genes older than the age fade to grey.
  function ageSecFor(id) {
    const e = state.data[id]; if (!e) return null;
    if (dfAgeSrc.value === "pron") {
      const ids = Object.keys(state.data);
      const dists = ids.map((k) => state.data[k].min_dist).filter((d) => d != null);
      if (e.min_dist == null || dists.length < 2) return null;
      const dmin = Math.min(...dists), dmax = Math.max(...dists);
      const pt = dmax > dmin ? (dmax - e.min_dist) / (dmax - dmin) : 0.5;   // 1 = oldest (smallest pronuclei distance)
      const avgs = ids.map((k) => embAvgSec(k, state.ic, state.metric)).filter((a) => a != null);
      if (!avgs.length) return null;
      const aMin = Math.min(...avgs), aMax = Math.max(...avgs);
      return aMin + pt * (aMax - aMin);
    }
    return embAvgSec(id, state.ic, state.metric);
  }
  const lerp = (a, b, t) => Math.round(a + (b - a) * t);
  const mix = (c1, c2, t) => `rgb(${lerp(c1[0], c2[0], t)},${lerp(c1[1], c2[1], t)},${lerp(c1[2], c2[2], t)})`;
  function ageColor(r) {
    if (r <= 1) { const t = Math.max(0, Math.min(1, 1 - r)); return mix([255, 255, 255], [136, 19, 55], t); }  // white → crimson
    const t = Math.max(0, Math.min(1, (r - 1) / 2)); return mix([255, 255, 255], [203, 213, 225], t);          // white → slate grey
  }
  const ageText = (r) => (r <= 1 && 1 - r > 0.52) ? "#fff" : "#0f172a";
  function renderAgeGrid() {
    if (!state.rec || !dfAgeGrid) return;
    const r = icRec(); if (!r) { dfAgeGrid.innerHTML = ""; return; }
    const src = dfAgeSrc.value, age = ageSecFor(state.currentId);
    if (age == null) {
      dfAgeLegend.innerHTML = ""; dfAgeNote.textContent = "";
      dfAgeGrid.innerHTML = `<div class="df-age-empty">No ${src === "pron" ? "pronuclei-distance value" : "mean diffusion time"} available for this embryo.</div>`;
      return;
    }
    const fa = fmtTime(age);
    dfAgeNote.textContent = `age = ${fa.v} ${fa.u} · ` + (src === "pron"
      ? "pronuclei-distance pseudotime (smaller distance = older), rescaled onto the cohort's mean-diffusion-time range"
      : "this embryo's mean gene diffusion time");
    dfAgeLegend.innerHTML =
      `<span class="df-lg-bar"><i>≪ age</i><i>= age</i><i>&gt; age</i></span>` +
      `<span class="df-lg-cap">colour intensity = how far a gene's diffusion time sits <b>below</b> the embryo's age — deepest = passive diffusion alone over-explains its spread (transport-assisted candidate); grey = older than the age</span>`;
    const rows = Object.keys(r.genes).map((g) => {
      const sec = timeSec(g, state.ic, state.metric);
      return sec == null ? null : { g, sec, n: r.genes[g].n, ratio: sec / age };
    }).filter(Boolean).sort((a, b) => a.ratio - b.ratio);
    dfAgeGrid.innerHTML = rows.map((x) => {
      const f = fmtTime(x.sec), bg = ageColor(x.ratio), fg = ageText(x.ratio);
      const mark = x.ratio < 0.97 ? "↓" : x.ratio > 1.03 ? "↑" : "=";
      const rel = x.ratio < 0.97 ? "less than age" : x.ratio > 1.03 ? "older than age" : "≈ age";
      return `<button type="button" class="df-tile${x.g === state.gene ? " current" : ""}" data-gene="${x.g}" ` +
        `style="background:${bg};color:${fg}" title="${x.g} · ${f.v} ${f.u} · ${x.ratio.toFixed(2)}× age (${rel}) · ${x.n.toLocaleString()} tx">` +
        `<span class="df-tile-top"><span class="df-tile-g">${x.g}</span><span class="df-tile-m">${mark}</span></span>` +
        `<span class="df-tile-t">${f.v} ${f.u}</span>` +
        `<span class="df-tile-r">${x.ratio.toFixed(2)}× age</span></button>`;
    }).join("");
  }
  function renderActiveDrawer() {
    if (!state.drawerOpen) return;
    if (state.drawerTab === "age") renderAgeGrid();
    else { renderScatter(); setTimeout(() => { try { Plotly.Plots.resize(dfScatter); } catch (_) {} }, 40); }
  }
  function setDrawerTab(tab) {
    state.drawerTab = tab;
    dfDrawerTabs.querySelectorAll(".df-dtab").forEach((b) => {
      const on = b.dataset.tab === tab; b.classList.toggle("active", on); b.setAttribute("aria-selected", String(on));
    });
    drawerBody.querySelectorAll(".df-dpanel").forEach((p) => { p.hidden = p.dataset.tab !== tab; });
    renderActiveDrawer();
  }

  // ---------- wiring ----------
  function onGeneChange(g) {
    if (!g || g === state.gene) return; state.gene = g; geneSelect.value = g;
    stopAnim(); state.frame = 0; pickSample(); computeStop();
    drawFrame(); renderReadout(); renderList(); renderTimePlot();
    if (state.drawerOpen) { if (state.drawerTab === "age") renderAgeGrid(); else if (dfYaxis.value === "gene_count") renderScatter(); }
  }
  function onMetricOrIc() {
    stopAnim(); state.frame = 0; populateGenes(); pickSample(); computeStop();
    drawFrame(); renderReadout(); renderList(); renderTimePlot(); renderActiveDrawer();
    V.buildTabs(tabsEl, state.points, selectEmbryo, (e) => ({ label: e.label, sub: e.date_short,
      title: `${e.label} · mean diffusion time ${fmtEmbAvg(e.id)}` }));
    V.markActiveTab(tabsEl, state.currentId);
  }
  function wireControls() {
    geneSelect.addEventListener("change", () => onGeneChange(geneSelect.value));
    metricSelect.addEventListener("change", () => { state.metric = metricSelect.value; onMetricOrIc(); });
    icSeg.addEventListener("click", (e) => { const b = e.target.closest(".df-seg-btn"); if (!b) return;
      icSeg.querySelectorAll(".df-seg-btn").forEach((x) => x.classList.toggle("active", x === b));
      state.ic = b.dataset.ic; onMetricOrIc(); });
    playBtn.addEventListener("click", () => (state.playing ? stopAnim() : play()));
    stopBtn.addEventListener("click", stopAnim);
    resetBtn.addEventListener("click", () => { stopAnim(); state.frame = 0; drawFrame(); renderReadout(); updateTimeMarker(); });
    [obsToggle, cellToggle].forEach((c) => c.addEventListener("change", () => { if (state.scene) drawFrame(); }));
    cellColorInp.addEventListener("input", () => { state.cellColor = cellColorInp.value; if (state.scene && cellToggle.checked) drawFrame(); });
    cellOpacityInp.addEventListener("input", () => {
      state.cellOpacity = +cellOpacityInp.value; cellOpacityVal.value = (+cellOpacityInp.value).toFixed(2);
      if (state.scene && cellToggle.checked) drawFrame();
    });
    const corners = [...controlsEl.querySelectorAll(".rz")];
    try { V.wireWindow(controlsEl, $("#controls-header"), corners, "df.win"); } catch (_) {}
  }
  function wireDrawer() {
    dfYaxis.addEventListener("change", renderScatter);
    dfModel.addEventListener("change", renderScatter);
    dfAgeSrc.addEventListener("change", renderAgeGrid);
    dfDrawerTabs.addEventListener("click", (e) => { const b = e.target.closest(".df-dtab"); if (b) setDrawerTab(b.dataset.tab); });
    dfAgeGrid.addEventListener("click", (e) => { const t = e.target.closest(".df-tile"); if (t) onGeneChange(t.dataset.gene); });
    let start = null, moved = false;
    drawerHandle.addEventListener("pointerdown", (e) => { if (e.button && e.button !== 0) return; start = { x: e.clientX, y: e.clientY }; moved = false; try { drawerHandle.setPointerCapture(e.pointerId); } catch (_) {} });
    drawerHandle.addEventListener("pointermove", (e) => { if (!start) return;
      if (!moved && Math.hypot(e.clientX - start.x, e.clientY - start.y) < 5) return;
      if (!moved) { moved = true; drawer.classList.add("dragging"); if (drawer.dataset.open !== "true") openDrawer(true); }
      drawer.style.setProperty("--drawer-h", Math.max(220, Math.min(innerHeight - 100, innerHeight - e.clientY - 40)) + "px"); e.preventDefault(); });
    const up = (e) => { if (!start) return; try { drawerHandle.releasePointerCapture(e.pointerId); } catch (_) {}
      if (moved) drawer.classList.remove("dragging"); else openDrawer(drawer.dataset.open !== "true"); start = null; moved = false; };
    drawerHandle.addEventListener("pointerup", up); drawerHandle.addEventListener("pointercancel", up);
    const rz = $("#drawer-resize"); let sh = 0;
    rz.addEventListener("pointerdown", (e) => { sh = innerHeight - drawer.getBoundingClientRect().top; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return; drawer.style.setProperty("--drawer-h", Math.max(220, Math.min(innerHeight - 100, sh - (e.clientY - rz._d.y))) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }
  function openDrawer(open) {
    state.drawerOpen = open; drawer.dataset.open = open ? "true" : "false"; drawerHandle.setAttribute("aria-expanded", String(open));
    if (open) renderActiveDrawer();
  }
  function wireRdrawer() {
    dfList.addEventListener("click", (e) => { const r = e.target.closest(".df-row"); if (r) onGeneChange(r.dataset.gene); });
    let open = false;
    const set = (o) => { open = o; rdrawer.dataset.open = o ? "true" : "false"; rdrawerHandle.setAttribute("aria-expanded", String(o)); };
    rdrawerHandle.addEventListener("click", () => set(!open));
    const rz = $("#rdrawer-resize");
    rz.addEventListener("pointerdown", (e) => { rz._d = { x: e.clientX }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return; rdrawer.style.setProperty("--rdrawer-w", Math.max(240, Math.min(560, innerWidth - e.clientX)) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }

  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false;
    placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
