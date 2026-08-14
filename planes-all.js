/* Division Plane Sweep — every plane. The unrestricted sibling of the 18-plane project:
 * for each gene we searched EVERY plane orientation through the cell centre (~20,000 normals at
 * ~1°) and stored that gene's single globally-best dividing plane (its own 3-D normal), under two
 * normalizations (density / count), with a search-corrected significance. The 3-D view shows the
 * selected gene's best plane and its transcript split. The bottom drawer overlays every carrier
 * zygote's cell outline rotated into that gene's own best-plane frame, the per-side split, and the
 * distribution of best-plane orientations relative to the polar-body axis. Precomputed in
 * build_planes_all.py; the JS reads the shared normal grid and the per-gene best planes.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const AXIS_C = "#111827", PLANE_C = "#f59e0b", SPERM_C = V.SPERM_COLOR;
  const BLUE = "#2166ac", RED = "#b2182b", GREEN = "#16a34a";
  const VIRIDIS = [[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],[31,158,137],[53,183,121],[110,206,88],[181,222,43],[253,231,37]];

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const geneSelect = $("#gene-select"), modeSelect = $("#mode-select");
  const planeShow = $("#plane-show"), axisShow = $("#axis-show"), spermShow = $("#sperm-show"), dotsShow = $("#dots-show");
  const allPlanesShow = $("#all-planes-show");
  const readoutEl = $("#pa-readout");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), drawerBody = $("#drawer-body"), drawerGene = $("#drawer-gene");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle"), bestList = $("#best-list");
  const xsTabs = $("#xs-tabs"), xsPanels = $("#xs-panels");

  const state = { manifest: [], byId: {}, normals: null, currentId: null, scene: null,
    gene: null, mode: "vol", drawerOpen: false, tab: "align", sort: "p", dotSize: V.DOT_SIZE,
    agg: null, _aggP: null, spermData: null, _spermP: null, vcExtras: null, _sceneCache: {} };

  const vec = (a) => [a[0], a[1], a[2]];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const norm = (a) => Math.hypot(a[0], a[1], a[2]);
  const unit = (a) => { const n = norm(a) || 1; return [a[0] / n, a[1] / n, a[2] / n]; };
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const scal = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const lerp3 = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
  function viridis(t) { t = Math.max(0, Math.min(1, t)); const x = t * 8, i = Math.min(7, Math.floor(x)); const c = lerp3(VIRIDIS[i], VIRIDIS[i + 1], x - i); return `rgb(${c.map(Math.round).join(",")})`; }
  const sigT = (p) => Math.max(0, Math.min(1, (Math.log10(Math.max(p, 1e-4)) + 3) / 3));   // p≤.001 dark → n.s. yellow
  const fmtP = (p) => (p == null || !isFinite(p) ? "–" : p < 1e-4 ? p.toExponential(1) : p < 0.1 ? p.toPrecision(2) : p.toFixed(2));
  const sigOk = (p) => p != null && p <= 0.05;

  function plotInto(div, traces, layout, cfg) {
    if (!div.classList.contains("js-plotly-plot")) div.innerHTML = "";
    Plotly.react(div, traces, layout, cfg || { responsive: true, displayModeBar: false });
  }
  const shown = (el) => !!(el && el.offsetParent);

  (async function init() {
    try {
      const [m, nrm] = await Promise.all([
        (await fetch("data/planes_all_manifest.json")).json(),
        V.loadGz("data/planes_all_normals.json.gz"),
      ]);
      state.manifest = m.embryos; state.normals = nrm.normals;
      m.embryos.forEach((e) => (state.byId[e.id] = e));
      countEl.textContent = `${m.embryos.length} zygotes · ${(m.m_planes / 1000).toFixed(0)}k candidate planes · ~1°`;
      V.buildTabs(tabsEl, m.embryos, selectEmbryo, (e) => ({ label: e.label, sub: e.date_short, title: `${e.label} · ${e.n_genes} genes` }));
      wireControls(); wireDrawer(); wireRdrawer();
    } catch (err) { showError("Failed to load: " + (err.message || err)); }
  })();

  const gene = () => geneSelect.value;
  const ki = () => (state.mode === "vol" ? 0 : 1);
  function geneRow(g) { const a = state.scene && state.scene.analysis; return a && a.genes.find((r) => r.gene === g); }
  // per-gene best-plane fields for the current mode
  function best(g) {
    const r = geneRow(g); if (!r) return null;
    const vol = state.mode === "vol";
    return { idx: vol ? r.iVol : r.iCnt, a: vol ? r.aVol : r.aCnt, total: r.total,
      eff: vol ? r.effVol : r.effCnt, p: vol ? r.pVol : r.pCnt,
      volA: vol ? r.volA_vol : r.volA_cnt, volB: vol ? r.volB_vol : r.volB_cnt };
  }
  const normalOf = (idx) => state.normals[idx];

  function populateGenes(scene) {
    const genes = scene.genes.slice().sort((a, b) => a.localeCompare(b, undefined, { sensitivity: "base", numeric: true }));
    const tot = scene.gene_totals || {};
    geneSelect.innerHTML = genes.map((g) => `<option value="${g}">${g} (${(tot[g] || 0).toLocaleString()})</option>`).join("");
    if (state.gene && genes.includes(state.gene)) geneSelect.value = state.gene;
    else { state.gene = scene.genes[0]; geneSelect.value = state.gene; }
  }

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id; V.markActiveTab(tabsEl, id);
    const e = state.byId[id]; showLoading(`Loading ${e.label}…`);
    try {
      let sc = state._sceneCache[id];
      if (!sc) { sc = await V.loadGz(`data/planes_all/${id}.json.gz`); state._sceneCache[id] = sc; }
      if (state.currentId !== id) return;
      state.scene = sc; populateGenes(sc);
      if (!state.vcExtras) state.vcExtras = V.addWindowExtras($("#controls-body"), { defaultSize: state.dotSize, onDotSize: (s) => { state.dotSize = s; render(); } });
      state.vcExtras.setAtlas(id);
      controlsEl.hidden = false; placeholder.hidden = true; drawer.hidden = false; rdrawer.hidden = false;
      render(); renderReadout(); renderRanks();
      if (!state.drawerOpen) openDrawer(true); else renderActive();
    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }

  function ensureAgg() { if (!state._aggP) state._aggP = V.loadGz("data/planes_all_cross.json.gz").then((a) => (state.agg = a)); return state._aggP; }
  function ensureSperm() { if (!state._spermP) state._spermP = fetch("data/zygote_sperm.json").then((r) => r.json()).then((d) => (state.spermData = d)).catch(() => null); return state._spermP; }

  // ───────────────────────── 3-D ─────────────────────────
  const toPlot = (pUm, zs) => [pUm[0] / XY, pUm[1] / XY, pUm[2] * zs];
  function planeQuad(comUm, n, L, zs, color, op, name, rank) {
    const ref = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const t = unit(cross(n, ref)), w = unit(cross(n, t));
    const corners = [[1, 1], [1, -1], [-1, -1], [-1, 1]].map(([s1, s2]) => {
      const pUm = [comUm[0] + L * (s1 * t[0] + s2 * w[0]), comUm[1] + L * (s1 * t[1] + s2 * w[1]), comUm[2] + L * (s1 * t[2] + s2 * w[2])];
      return toPlot(pUm, zs);
    });
    return { type: "mesh3d", x: corners.map((c) => c[0]), y: corners.map((c) => c[1]), z: corners.map((c) => c[2]),
      i: [0, 0], j: [1, 2], k: [2, 3], color, opacity: op, name, showlegend: true, legendrank: rank, hoverinfo: "name", flatshading: true };
  }
  function render() {
    const s = state.scene; if (!s) return;
    const A = s.analysis, g = gene(), zs = s.z_scale, b = best(g);
    const traces = V.bodyTraces(s);
    if (dotsShow.checked && s.transcripts[g] && b) {
      const n = normalOf(b.idx), com = A.com_um, t = s.transcripts[g];
      const ax = [], ay = [], az = [], bx = [], by = [], bz = [], gx = [], gy = [], gz = [];
      for (let k = 0; k < t.x.length; k++) {
        const zp = t.gz[k] * zs;
        if (!t.s1[k]) { gx.push(t.x[k]); gy.push(t.y[k]); gz.push(zp); continue; }
        const side = (t.x[k] * XY - com[0]) * n[0] + (t.y[k] * XY - com[1]) * n[1] + (t.gz[k] * 1.0 - com[2]) * n[2];
        if (side > 0) { ax.push(t.x[k]); ay.push(t.y[k]); az.push(zp); } else { bx.push(t.x[k]); by.push(t.y[k]); bz.push(zp); }
      }
      traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side A`, x: ax, y: ay, z: az, marker: { size: state.dotSize, color: BLUE, opacity: V.DOT_OPACITY, line: { width: 0 } }, hovertemplate: `${g} · side A<extra></extra>`, legendrank: 20000 });
      traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side B`, x: bx, y: by, z: bz, marker: { size: state.dotSize, color: RED, opacity: V.DOT_OPACITY, line: { width: 0 } }, hovertemplate: `${g} · side B<extra></extra>`, legendrank: 20001 });
      if (gx.length) traces.push({ type: "scatter3d", mode: "markers", name: `${g} · not counted`, x: gx, y: gy, z: gz, marker: { size: state.dotSize, color: GREEN, opacity: 0.6, line: { width: 0 } }, hovertemplate: `${g} · not counted<extra></extra>`, legendrank: 20002 });
    }
    if (axisShow.checked) {
      const c = A.com_plot, ap = A.axis_plot, an = norm(ap) || 1, ex = s.extents;
      const R = 0.62 * Math.max(ex.x[1] - ex.x[0], ex.y[1] - ex.y[0], ex.z[1] - ex.z[0]);
      const u = [ap[0] / an, ap[1] / an, ap[2] / an];
      traces.push({ type: "scatter3d", mode: "lines", name: "Polar-body axis", x: [c[0] - R * u[0], c[0] + R * u[0]], y: [c[1] - R * u[1], c[1] + R * u[1]], z: [c[2] - R * u[2], c[2] + R * u[2]], line: { color: AXIS_C, width: 6 }, hovertemplate: "Polar-body axis<extra></extra>", legendrank: 40000 });
      traces.push({ type: "scatter3d", mode: "markers", name: "Polar body", x: [A.pb_plot[0]], y: [A.pb_plot[1]], z: [A.pb_plot[2]], marker: { size: 7, color: AXIS_C, line: { width: 1, color: "#fff" } }, hovertemplate: "Polar body<extra></extra>", legendrank: 40001 });
    }
    // every gene's best plane for this zygote, as ONE mesh coloured per-plane by significance
    if (allPlanesShow.checked && A.genes) {
      const com = A.com_um, L = A.L_um * 1.5, vol = state.mode === "vol";
      const X = [], Y = [], Z = [], I = [], J = [], K = [], FC = []; let base = 0;
      A.genes.forEach((r) => {
        const idx = vol ? r.iVol : r.iCnt, p = vol ? r.pVol : r.pCnt;
        if (idx == null) return;
        const n = normalOf(idx), ref = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
        const t = unit(cross(n, ref)), w = unit(cross(n, t)), col = viridis(sigT(p));
        [[1, 1], [1, -1], [-1, -1], [-1, 1]].forEach(([s1, s2]) => {
          const pp = toPlot([com[0] + L * (s1 * t[0] + s2 * w[0]), com[1] + L * (s1 * t[1] + s2 * w[1]), com[2] + L * (s1 * t[2] + s2 * w[2])], zs);
          X.push(pp[0]); Y.push(pp[1]); Z.push(pp[2]);
        });
        I.push(base, base); J.push(base + 1, base + 2); K.push(base + 2, base + 3); FC.push(col, col); base += 4;
      });
      if (X.length) traces.push({ type: "mesh3d", x: X, y: Y, z: Z, i: I, j: J, k: K, facecolor: FC,
        opacity: 0.13, name: "all gene planes", showlegend: true, legendrank: 30000, hoverinfo: "skip", flatshading: true });
    }
    if (planeShow.checked && b) traces.push(planeQuad(A.com_um, normalOf(b.idx), A.L_um * 1.6, zs, PLANE_C, 0.3, `${g} best plane`, 41000));
    if (spermShow.checked && state.spermData) {
      const se = (state.spermData.embryos || []).find((x) => x.id === state.currentId);
      if (se && se.sperm_plot) traces.push({ type: "scatter3d", mode: "markers", name: "sperm", x: [se.sperm_plot[0]], y: [se.sperm_plot[1]], z: [se.sperm_plot[2]], marker: { size: 8, color: SPERM_C, symbol: "diamond", line: { width: 1, color: "#fff" } }, hovertemplate: "sperm<extra></extra>", legendrank: 40002 });
    }
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  function angleToAxis(nrm, axisUm) {
    if (!axisUm) return null;
    const c = Math.abs(dot(nrm, axisUm)) / (norm(nrm) * norm(axisUm) || 1);
    return Math.acos(Math.max(0, Math.min(1, c))) * 180 / Math.PI;
  }
  function renderReadout() {
    const g = gene(), b = best(g), A = state.scene.analysis;
    if (!b) { readoutEl.innerHTML = `<div class="pa-r-na"><b>${g}</b> — not in this zygote's panel</div>`; return; }
    const enr = b.a * 2 >= b.total, hi = Math.max(b.a, b.total - b.a), frac = hi / b.total;
    const ang = angleToAxis(normalOf(b.idx), A.axis_um);
    const geom = ang == null ? "" : ang > 60 ? "≈ meridional (contains the axis)" : ang < 30 ? "≈ equatorial (⟂ the axis)" : "oblique";
    readoutEl.innerHTML =
      `<div class="pa-r-line"><b>${g}</b> · best ${state.mode === "vol" ? "density" : "count"} plane of ${(state.normals.length / 1000).toFixed(0)}k</div>` +
      `<div class="pa-big"><span>${(frac * 100).toFixed(1)}%</span> <span class="pa-lbl">on the higher side (${b.total.toLocaleString()} transcripts)</span></div>` +
      `<div class="pa-r-line">split |Δ| = <b>${state.mode === "vol" ? b.eff.toExponential(2) : (b.eff).toFixed(3)}</b> · search-corrected <span class="pa-pval">p = <b class="${sigOk(b.p) ? "pa-sig" : ""}">${fmtP(b.p)}</b></span></div>` +
      (ang == null ? "" : `<div class="pa-r-line">plane vs polar axis: <b>${ang.toFixed(0)}°</b> · ${geom}</div>`);
  }

  // ───────────────────────── drawer ─────────────────────────
  const RENDER = { align: renderAlign, bars: renderBars, orient: renderOrient };
  function renderActive() { (RENDER[state.tab] || renderAlign)(); }
  function switchTab(which) {
    if (!RENDER[which]) which = "align";
    state.tab = which;
    xsTabs.querySelectorAll(".xs-gtab").forEach((t) => { const on = t.dataset.tab === which; t.classList.toggle("active", on); t.setAttribute("aria-selected", String(on)); });
    xsPanels.querySelectorAll(".xs-panel").forEach((p) => (p.hidden = p.dataset.tab !== which));
    renderActive(); requestAnimationFrame(resizeAll);
  }
  const resizeAll = () => ["align-plot", "bars-plot", "orient-plot"].forEach((id) => { try { Plotly.Plots.resize($("#" + id)); } catch (_) {} });
  function baseLayout(xt, yt) {
    return { margin: { l: 50, r: 12, t: 8, b: 40 }, autosize: true,
      xaxis: { title: { text: xt, font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", zeroline: false },
      yaxis: { title: { text: yt, font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", zeroline: false },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", legend: { orientation: "h", font: { size: 10 }, y: 1.14, x: 1, xanchor: "right" } };
  }

  // silhouette (max radius per angular bin) of 2-D points centred at origin
  function silhouette(pts, nbins) {
    nbins = nbins || 96; const rad = new Array(nbins).fill(0);
    for (const [x, y] of pts) { const a = Math.atan2(y, x); const bi = ((a + Math.PI) / (2 * Math.PI) * nbins) | 0; const r = Math.hypot(x, y); const j = ((bi % nbins) + nbins) % nbins; if (r > rad[j]) rad[j] = r; }
    // fill empty bins from nearest non-empty neighbours
    for (let b = 0; b < nbins; b++) if (rad[b] === 0) {
      let l = 0, rr = 0;
      for (let k = 1; k < nbins; k++) { const v = rad[(b - k + nbins) % nbins]; if (v > 0) { l = v; break; } }
      for (let k = 1; k < nbins; k++) { const v = rad[(b + k) % nbins]; if (v > 0) { rr = v; break; } }
      rad[b] = (l + rr) / 2;
    }
    // light circular smoothing (5-tap) to tame the max-radius spikes from a sparse point set
    const sm = rad.map((_, b) => { let s = 0, w = 0; for (let k = -2; k <= 2; k++) { const wt = 3 - Math.abs(k); s += wt * rad[(b + k + nbins) % nbins]; w += wt; } return s / w; });
    const out = [];
    for (let b = 0; b < nbins; b++) { const th = (b + 0.5) / nbins * 2 * Math.PI - Math.PI; out.push([sm[b] * Math.cos(th), sm[b] * Math.sin(th)]); }
    out.push(out[0]); return out;
  }
  // project a zygote's surface points into its gene-best-plane frame; flip so the higher side is +x
  function alignedOutline(e, gRow) {
    const surf = e.surf, axis = e.axis; if (!surf || !surf.length) return null;
    const idx = gRow[ki()], aAt = gRow[4 + ki()], total = (e.g[gene()] || [0])[0];
    const n = normalOf(idx);
    let t = axis ? sub(axis, scal(n, dot(axis, n))) : [0, 0, 0];
    if (norm(t) < 1e-6) t = cross(n, Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0]);
    t = unit(t);
    const flip = aAt * 2 < total ? -1 : 1;
    const pts = surf.map((s) => [flip * dot(s, n), dot(s, t)]);
    return { outline: silhouette(pts), p: gRow[8 + ki()], flip };
  }
  function renderAlign() {
    const div = $("#align-plot"); if (!shown(div)) return;
    ensureAgg().then(() => {
      if (!shown(div)) return;
      const g = gene(), agg = state.agg;
      const carriers = agg.embryos.filter((e) => e.gb && e.gb[g]);
      $("#align-sub").textContent = `· ${g} · ${carriers.length} of ${agg.n_embryos} zygotes carry this gene`;
      if (!carriers.length) { Plotly.purge(div); div.classList.remove("js-plotly-plot"); div.innerHTML = `<div class="pa-empty"><b>${g}</b> is not detected in enough zygotes.</div>`; return; }
      const traces = [], showCells = $("#align-cells").checked, showMean = $("#align-mean").checked, showPlane = $("#align-plane").checked;
      const aligned = [];
      let lim = 20;
      for (const e of carriers) { const o = alignedOutline(e, e.gb[g]); if (!o) continue; aligned.push({ e, ...o }); for (const p of o.outline) lim = Math.max(lim, Math.abs(p[0]), Math.abs(p[1])); }
      if (showCells) for (const o of aligned) { if (o.e.id === state.currentId) continue;
        traces.push({ type: "scatter", mode: "lines", x: o.outline.map((p) => p[0]), y: o.outline.map((p) => p[1]), line: { color: viridis(sigT(o.p)), width: 1.1, shape: "spline" }, opacity: 0.5, hoverinfo: "skip", showlegend: false }); }
      if (showMean && aligned.length > 1) {
        const nb = 120, acc = new Array(nb).fill(0), cnt = new Array(nb).fill(0);
        for (const o of aligned) o.outline.slice(0, nb).forEach((p, b) => { acc[b] += Math.hypot(p[0], p[1]); cnt[b]++; });
        const mo = []; for (let b = 0; b < nb; b++) { const r = cnt[b] ? acc[b] / cnt[b] : 0; const th = (b + 0.5) / nb * 2 * Math.PI - Math.PI; mo.push([r * Math.cos(th), r * Math.sin(th)]); } mo.push(mo[0]);
        traces.push({ type: "scatter", mode: "lines", name: "mean", x: mo.map((p) => p[0]), y: mo.map((p) => p[1]), line: { color: "#0f172a", width: 2.4, shape: "spline" }, hoverinfo: "skip" });
      }
      const cur = aligned.find((o) => o.e.id === state.currentId);
      if (cur) { traces.push({ type: "scatter", mode: "lines", x: cur.outline.map((p) => p[0]), y: cur.outline.map((p) => p[1]), line: { color: "#fff", width: 5 }, hoverinfo: "skip", showlegend: false });
        traces.push({ type: "scatter", mode: "lines", name: "this zygote", x: cur.outline.map((p) => p[0]), y: cur.outline.map((p) => p[1]), line: { color: viridis(sigT(cur.p)), width: 2.6, shape: "spline" }, hoverinfo: "skip" }); }
      lim *= 1.08;
      if (showPlane) traces.push({ type: "scatter", mode: "lines", name: "best plane", x: [0, 0], y: [-lim, lim], line: { color: "#111827", width: 1.4, dash: "dash" }, hoverinfo: "skip" });
      const lay = baseLayout("distance from best plane (µm) · higher side →", "polar-body axis (µm)");
      lay.xaxis.range = [-lim, lim]; lay.yaxis.range = [-lim, lim]; lay.xaxis.scaleanchor = "y"; lay.xaxis.scaleratio = 1;
      plotInto(div, traces, lay);
    });
  }
  function renderBars() {
    const div = $("#bars-plot"); if (!shown(div)) return;
    ensureAgg().then(() => {
      if (!shown(div)) return;
      const g = gene(), agg = state.agg;
      const rows = agg.embryos.filter((e) => e.gb && e.gb[g]).map((e) => { const gb = e.gb[g], aAt = gb[4 + ki()], total = (e.g[g] || [0])[0]; return { e, frac: total ? Math.max(aAt, total - aAt) / total : 0.5, total, p: gb[8 + ki()] }; });
      $("#bars-sub").textContent = `· ${g} · ${rows.length} zygotes`;
      if (!rows.length) { Plotly.purge(div); div.classList.remove("js-plotly-plot"); div.innerHTML = `<div class="pa-empty">No carrier zygotes.</div>`; return; }
      rows.sort((a, b) => b.frac - a.frac);
      const colr = rows.map((r) => (sigOk(r.p) ? viridis(sigT(r.p)) : "#cbd5e1"));
      const tr = { type: "bar", x: rows.map((r) => r.e.label), y: rows.map((r) => r.frac * 100), marker: { color: colr },
        customdata: rows.map((r) => [r.total, r.p]), hovertemplate: "%{x}: %{y:.1f}% on higher side · n=%{customdata[0]} · p=%{customdata[1]:.2g}<extra></extra>" };
      const lay = baseLayout("", "% of transcripts on the higher-density side");
      lay.xaxis.tickangle = -48; lay.xaxis.tickfont = { size: 9 }; lay.yaxis.range = [40, 100]; lay.showlegend = false;
      lay.shapes = [{ type: "line", x0: -0.5, x1: rows.length - 0.5, y0: 50, y1: 50, line: { color: "#94a3b8", width: 1, dash: "dash" } }];
      plotInto(div, [tr], lay);
    });
  }
  function renderOrient() {
    const div = $("#orient-plot"); if (!shown(div)) return;
    ensureAgg().then(() => {
      if (!shown(div)) return;
      const g = gene(), agg = state.agg;
      const angs = [];
      for (const e of agg.embryos) { if (!e.gb || !e.gb[g] || !e.axis) continue; const a = angleToAxis(normalOf(e.gb[g][ki()]), e.axis); if (a != null) angs.push(a); }
      $("#orient-sub").textContent = `· ${g} · ${angs.length} zygotes`;
      if (!angs.length) { Plotly.purge(div); div.classList.remove("js-plotly-plot"); div.innerHTML = `<div class="pa-empty">No carrier zygotes with an axis.</div>`; return; }
      const tr = { type: "histogram", x: angs, xbins: { start: 0, end: 90, size: 10 }, marker: { color: "#f59e0b", line: { color: "#b45309", width: 1 } }, hovertemplate: "%{x}°: %{y} zygotes<extra></extra>" };
      const lay = baseLayout("angle between best plane & polar-body axis (°)", "zygotes");
      lay.xaxis.range = [0, 90]; lay.xaxis.tickvals = [0, 30, 45, 60, 90]; lay.showlegend = false;
      lay.annotations = [{ x: 0, y: 1.02, yref: "paper", xref: "x", text: "equatorial", showarrow: false, font: { size: 9, color: "#64748b" }, xanchor: "left" },
        { x: 90, y: 1.02, yref: "paper", xref: "x", text: "meridional", showarrow: false, font: { size: 9, color: "#64748b" }, xanchor: "right" }];
      plotInto(div, [tr], lay);
    });
  }

  // ───────────────────────── right drawer ─────────────────────────
  function renderRanks() {
    const A = state.scene && state.scene.analysis; if (!A) return;
    const vol = state.mode === "vol";
    const rows = A.genes.map((r) => ({ gene: r.gene, total: r.total, eff: vol ? r.effVol : r.effCnt, p: vol ? r.pVol : r.pCnt }));
    rows.sort((a, b) => (state.sort === "p" ? (a.p - b.p || b.eff - a.eff) : (b.eff - a.eff || a.p - b.p)));
    const top = rows.slice(0, 80), curG = gene();
    let html = `<div class="best-head"><span></span><span>gene</span><span>${state.mode === "vol" ? "|Δ|·V" : "|Δ|"}</span><span>p</span></div>`;
    html += top.map((r, i) => {
      const eff = vol ? r.eff.toExponential(1) : r.eff.toFixed(2);
      return `<div class="best-row${r.gene === curG ? " current" : ""}" data-gene="${r.gene}" title="${r.gene}: n=${r.total}, split |Δ|=${eff}, search-corrected p=${fmtP(r.p)}">` +
        `<span class="best-num">${i + 1}</span><span class="best-gene">${r.gene}</span><span class="best-real">${eff}</span>` +
        `<span class="best-p${sigOk(r.p) ? " sig" : ""}">${fmtP(r.p)}</span></div>`;
    }).join("");
    bestList.innerHTML = html;
  }
  function highlightRank() { const g = gene(); bestList.querySelectorAll(".best-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === g)); }

  // ───────────────────────── wiring ─────────────────────────
  function onGene() { state.gene = gene(); drawerGene.textContent = `· ${gene()}`; render(); renderReadout(); renderActive(); highlightRank(); }
  function wireControls() {
    geneSelect.addEventListener("change", onGene);
    modeSelect.addEventListener("change", () => { state.mode = modeSelect.value; render(); renderReadout(); renderActive(); renderRanks(); });
    [planeShow, axisShow, dotsShow, allPlanesShow].forEach((c) => c.addEventListener("change", render));
    spermShow.addEventListener("change", () => { if (spermShow.checked) ensureSperm().then(render); else render(); });
  }
  function openDrawer(open) { state.drawerOpen = open; drawer.dataset.open = open ? "true" : "false"; drawerHandle.setAttribute("aria-expanded", String(open)); drawerGene.textContent = `· ${gene()}`; if (open) { renderActive(); requestAnimationFrame(resizeAll); } }
  function wireDrawer() {
    xsTabs.addEventListener("click", (e) => { const t = e.target.closest(".xs-gtab"); if (t) switchTab(t.dataset.tab); });
    ["align-cells", "align-mean", "align-plane"].forEach((id) => $("#" + id).addEventListener("change", renderAlign));
    $("#align-download").addEventListener("click", () => dl($("#align-plot"), "aligned"));
    $("#bars-download").addEventListener("click", () => dl($("#bars-plot"), "per-side"));
    $("#orient-download").addEventListener("click", () => dl($("#orient-plot"), "orientations"));
    wireHandleDrag(drawer, drawerHandle, { computeSize: (e) => window.innerHeight - e.clientY - 40, clampSize: (px) => Math.max(200, Math.min(window.innerHeight - 100, px)), applySize: (px) => drawer.style.setProperty("--drawer-h", px + "px"), setOpen: openDrawer, afterDrag: resizeAll });
    const rz = $("#drawer-resize"); let sh = 0;
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return; drawer.style.setProperty("--drawer-h", Math.max(200, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y))) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} resizeAll(); } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
    xsPanels.querySelectorAll(".xs-resizable").forEach((box) => { const plot = box.querySelector(".xs-plot"); new ResizeObserver(() => { if (box.offsetParent) { try { Plotly.Plots.resize(plot); } catch (_) {} } }).observe(box); });
  }
  function dl(div, name) { try { Plotly.downloadImage(div, { format: "png", scale: 4, width: 1600, height: 1100, filename: `planes-all-${state.currentId}-${gene()}-${name}` }); } catch (_) {} }
  function wireRdrawer() {
    $("#rtabs").addEventListener("click", (e) => { const t = e.target.closest(".rtab"); if (!t) return; state.sort = t.dataset.sort; $("#rtabs").querySelectorAll(".rtab").forEach((b) => b.classList.toggle("active", b === t)); renderRanks(); });
    bestList.addEventListener("click", (e) => { const row = e.target.closest(".best-row"); if (row && row.dataset.gene) { geneSelect.value = row.dataset.gene; onGene(); } });
    wireHandleDrag(rdrawer, rdrawerHandle, { computeSize: (e) => window.innerWidth - e.clientX, clampSize: (px) => Math.max(240, Math.min(window.innerWidth - 80, px)), applySize: (px) => rdrawer.style.setProperty("--rdrawer-w", px + "px"), setOpen: (o) => { rdrawer.dataset.open = o ? "true" : "false"; rdrawerHandle.setAttribute("aria-expanded", String(o)); } });
    const rrz = $("#rdrawer-resize"); let sw = 0;
    rrz.addEventListener("pointerdown", (e) => { sw = rdrawer.getBoundingClientRect().width; rrz._d = { x: e.clientX }; rrz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rrz.addEventListener("pointermove", (e) => { if (!rrz._d) return; rdrawer.style.setProperty("--rdrawer-w", Math.max(240, Math.min(window.innerWidth - 80, sw - (e.clientX - rrz._d.x))) + "px"); });
    const end = (e) => { if (rrz._d) { rrz._d = null; try { rrz.releasePointerCapture(e.pointerId); } catch (_) {} } };
    rrz.addEventListener("pointerup", end); rrz.addEventListener("pointercancel", end);
  }
  function wireHandleDrag(el, handle, cfg) {
    let start = null, moved = false;
    handle.addEventListener("pointerdown", (e) => { if (e.button && e.button !== 0) return; start = { x: e.clientX, y: e.clientY }; moved = false; try { handle.setPointerCapture(e.pointerId); } catch (_) {} });
    handle.addEventListener("pointermove", (e) => { if (!start) return; if (!moved && Math.hypot(e.clientX - start.x, e.clientY - start.y) < 5) return; if (!moved) { moved = true; el.classList.add("dragging"); if (el.dataset.open !== "true") cfg.setOpen(true); } cfg.applySize(cfg.clampSize(cfg.computeSize(e))); e.preventDefault(); });
    const up = (e) => { if (!start) return; try { handle.releasePointerCapture(e.pointerId); } catch (_) {} if (moved) { el.classList.remove("dragging"); cfg.afterDrag && cfg.afterDrag(); } else cfg.setOpen(el.dataset.open !== "true"); start = null; moved = false; };
    handle.addEventListener("pointerup", up); handle.addEventListener("pointercancel", up);
  }
  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false; placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
