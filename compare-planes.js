/* Compare Division Planes.
 *
 * Four candidate dividing planes per zygote — polar-axis best, exhaustive best, equatorial,
 * and the sperm·COM·polar-body plane — compared under ONE metric: the concentration asymmetry
 * of a gene's transcripts on the CELL BODY (pronuclei + polar body excluded) across the plane.
 * Data: data/compare_planes.json.gz (build_compare_planes.py) ships per-side counts + volumes;
 * everything (asymmetry, p, angles, alignment) is derived here. The 3-D meshes + transcript
 * clouds are reused from data/planes_all/<id>.json.gz.
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const XY = 0.15, Z_UM = 1.0;
  const BLUE = "#2563eb", RED = "#dc2626", GREEN = "#10b981", AXIS_C = "#6d28d9";

  const state = {
    data: null, byId: {}, planes: [], colorOf: {}, labelOf: {},
    currentId: null, scene: null, sceneCache: {},
    gene: null, planesOn: { polar: true, exhaustive: true, equatorial: true, sperm: true },
    splitBy: "equatorial", dotsOn: true, dotSize: 1.5,
    tab: "arrange", pairA: "equatorial", pairB: "polar",
    rankMode: "agree", rankPlanes: { polar: false, exhaustive: true, equatorial: true, sperm: false },
    minCount: 10, drawerOpen: false, rdrawerOpen: false, vcExtras: null, _rankKey: null, _rankRows: null,
  };

  // ───────── binomial statistics (asymmetry vs a uniform split) ─────────
  function gammaln(x) { const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, t = x + 5.5; t -= (x + 0.5) * Math.log(t); let s = 1.000000000190015; for (let j = 0; j < 6; j++) s += c[j] / ++y; return -t + Math.log(2.5066282746310005 * s / x); }
  function betacf(a, b, x) { const FP = 1e-300; let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FP) d = FP; d = 1 / d; let h = d; for (let m = 1; m <= 200; m++) { const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2)); d = 1 + aa * d; if (Math.abs(d) < FP) d = FP; c = 1 + aa / c; if (Math.abs(c) < FP) c = FP; d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)); d = 1 + aa * d; if (Math.abs(d) < FP) d = FP; c = 1 + aa / c; if (Math.abs(c) < FP) c = FP; d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < 3e-12) break; } return h; }
  function betai(a, b, x) { if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(gammaln(a + b) - gammaln(a) - gammaln(b) + a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b; }
  const binomSF = (k, n, p) => (k <= 0 ? 1 : k > n ? 0 : betai(k, n - k + 1, p));          // P(X≥k)
  const binomCDF = (k, n, p) => (k >= n ? 1 : k < 0 ? 0 : 1 - betai(k + 1, n - k, p));       // P(X≤k)
  function binomTwoSided(obs, n, f) { if (!n) return 1; const lo = binomCDF(obs, n, f), hi = binomSF(obs, n, f); return Math.min(1, 2 * Math.min(lo, hi)); }
  const fmtP = (p) => (p == null || !isFinite(p)) ? "–" : p < 1e-4 ? p.toExponential(1) : p < 0.1 ? p.toPrecision(2) : p.toFixed(2);
  const sig = (p) => p != null && p <= 0.05;

  // ───────── plane records (unify fixed + per-gene planes) ─────────
  function planeRec(e, g, key) {
    const rec = e.genes[g]; if (!rec) return null;
    if (key === "polar") { const b = rec.pb; return b ? { n: b.n, vA: b.v[0], vB: b.v[1], nA: b.c[0], nB: b.c[1] } : null; }
    if (key === "exhaustive") { const b = rec.ex; return b ? { n: b.n, vA: b.v[0], vB: b.v[1], nA: b.c[0], nB: b.c[1] } : null; }
    if (key === "equatorial") { if (!e.eq) return null; const nA = rec.eq; return { n: e.eq.n, vA: e.eq.v[0], vB: e.eq.v[1], nA, nB: rec.nc - nA }; }
    if (key === "sperm") { if (!e.sd || rec.sd == null) return null; const nA = rec.sd; return { n: e.sd.n, vA: e.sd.v[0], vB: e.sd.v[1], nA, nB: rec.nc - nA }; }
    return null;
  }
  const conc = (r) => [r.nA / Math.max(r.vA, 1e-9), r.nB / Math.max(r.vB, 1e-9)];
  function asym(r) { const [a, b] = conc(r); const s = a + b; return s > 0 ? Math.abs(a - b) / s : 0; }
  function pOf(r) { const n = r.nA + r.nB, f = r.vA / Math.max(r.vA + r.vB, 1e-9); return binomTwoSided(r.nA, n, f); }
  const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const norm3 = (a) => Math.hypot(a[0], a[1], a[2]) || 1;
  function angleDeg(n1, n2) { const d = Math.abs(dot3(n1, n2)) / (norm3(n1) * norm3(n2)); return Math.acos(Math.max(0, Math.min(1, d))) * 180 / Math.PI; } // 0=same plane, 90=orthogonal

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await V.loadGz("data/compare_planes.json.gz"); }
    catch (e) { $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div><div class="ph-sub">Run <code>python3 build_compare_planes.py</code>.</div></div>`; return; }
    const d = state.data;
    d.embryos.forEach((e) => (state.byId[e.id] = e));
    state.planes = d.planes;
    d.planes.forEach((p) => { state.colorOf[p.key] = p.color; state.labelOf[p.key] = p.label; });
    $("#embryo-count").textContent = `${d.n} zygotes · ${d.n_sperm} with sperm · 4 candidate planes`;
    fillPairSelects(); buildRankPlanes();
    V.buildTabs($("#tabs"), d.embryos, selectEmbryo, (e) => ({ label: e.label, sub: e.date_short, title: `${e.label} · ${Object.keys(e.genes).length} genes${e.has_sperm ? "" : " · no sperm"}` }));
    wire();
    selectEmbryo(d.embryos[0].id);
  })();

  function fillGenes() {
    const e = cur(); const genes = Object.keys(e.genes).sort((a, b) => a.localeCompare(b, undefined, { sensitivity: "base", numeric: true }));
    const sel = $("#gene-select");
    sel.innerHTML = genes.map((g) => `<option value="${g}">${g} (${e.genes[g].nc})</option>`).join("");
    if (!state.gene || !(state.gene in e.genes)) state.gene = genes[0];
    sel.value = state.gene;
  }
  const cur = () => state.byId[state.currentId];
  const gene = () => state.gene;

  async function selectEmbryo(id) {
    if (id === state.currentId && state.scene) return;
    state.currentId = id; V.markActiveTab($("#tabs"), id);
    const e = state.byId[id]; $("#loading").hidden = false; $("#loading-text").textContent = `Loading ${e.label}…`;
    try {
      let sc = state.sceneCache[id];
      if (!sc) { sc = await V.loadGz(`data/planes_all/${id}.json.gz`); state.sceneCache[id] = sc; }
      if (state.currentId !== id) return;
      state.scene = sc; fillGenes();
      if (!state.vcExtras) state.vcExtras = V.addWindowExtras($("#controls-body"), { defaultSize: state.dotSize, onDotSize: (s) => { state.dotSize = s; render3D(); } });
      state.vcExtras.setAtlas && state.vcExtras.setAtlas(id);
      $("#controls").hidden = false; $("#placeholder").hidden = true; $("#drawer").hidden = false; $("#rdrawer").hidden = false;
      render3D(); renderReadout(); renderActive();
      if (!state.drawerOpen) openDrawer(true);
    } catch (err) { $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Scene missing</div><div class="ph-sub">${err.message || err}</div></div>`; }
    finally { $("#loading").hidden = true; }
  }

  // ───────── 3-D ─────────
  const toPlot = (pUm, zs) => [pUm[0] / XY, pUm[1] / XY, pUm[2] * zs];
  function unitv(n) { const m = norm3(n); return [n[0] / m, n[1] / m, n[2] / m]; }
  function planeQuad(comUm, n, L, zs, color, op, name, rank) {
    const nn = unitv(n), ref = Math.abs(nn[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const cr = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    const t = unitv(cr(nn, ref)), w = unitv(cr(nn, t));
    const C = [[1, 1], [1, -1], [-1, -1], [-1, 1]].map(([s1, s2]) => toPlot([comUm[0] + L * (s1 * t[0] + s2 * w[0]), comUm[1] + L * (s1 * t[1] + s2 * w[1]), comUm[2] + L * (s1 * t[2] + s2 * w[2])], zs));
    return { type: "mesh3d", x: C.map((c) => c[0]), y: C.map((c) => c[1]), z: C.map((c) => c[2]), i: [0, 0], j: [1, 2], k: [2, 3], color, opacity: op, name, showlegend: true, legendrank: rank, hoverinfo: "name", flatshading: true };
  }
  function render3D() {
    const s = state.scene, e = cur(), g = gene(); if (!s || !e) return;
    const zs = s.z_scale, com = e.com_um, traces = V.bodyTraces(s);
    // gene dots, split by the chosen plane
    const splitRec = planeRec(e, g, state.splitBy);
    if (state.dotsOn && s.transcripts[g]) {
      const t = s.transcripts[g], n = splitRec ? unitv(splitRec.n) : null;
      const ax = [], ay = [], az = [], bx = [], by = [], bz = [], gx = [], gy = [], gz = [];
      for (let k = 0; k < t.x.length; k++) {
        const zp = t.gz[k] * zs;
        if (t.s1 && !t.s1[k]) { gx.push(t.x[k]); gy.push(t.y[k]); gz.push(zp); continue; }
        if (!n) { bx.push(t.x[k]); by.push(t.y[k]); bz.push(zp); continue; }
        const side = (t.x[k] * XY - com[0]) * n[0] + (t.y[k] * XY - com[1]) * n[1] + (t.gz[k] * Z_UM - com[2]) * n[2];
        if (side > 0) { ax.push(t.x[k]); ay.push(t.y[k]); az.push(zp); } else { bx.push(t.x[k]); by.push(t.y[k]); bz.push(zp); }
      }
      if (ax.length) traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side A`, x: ax, y: ay, z: az, marker: { size: state.dotSize, color: BLUE, opacity: 0.85, line: { width: 0 } }, hovertemplate: `${g} · A<extra></extra>`, legendrank: 20000 });
      if (bx.length) traces.push({ type: "scatter3d", mode: "markers", name: `${g} · side B`, x: bx, y: by, z: bz, marker: { size: state.dotSize, color: RED, opacity: 0.85, line: { width: 0 } }, hovertemplate: `${g} · B<extra></extra>`, legendrank: 20001 });
      if (gx.length) traces.push({ type: "scatter3d", mode: "markers", name: `${g} · not counted`, x: gx, y: gy, z: gz, marker: { size: state.dotSize, color: GREEN, opacity: 0.5, line: { width: 0 } }, hovertemplate: `${g} · pron/polar/other<extra></extra>`, legendrank: 20002 });
    }
    // planes
    state.planes.forEach((p, i) => {
      if (!state.planesOn[p.key]) return;
      const r = planeRec(e, g, p.key); if (!r) return;
      traces.push(planeQuad(com, r.n, e.L_um * 1.4, zs, p.color, 0.28, `${p.label} plane`, 41000 + i));
    });
    // polar axis line
    const a = e.axis_um, R = e.L_um; const cp = toPlot(com, zs);
    const p0 = toPlot([com[0] - R * a[0], com[1] - R * a[1], com[2] - R * a[2]], zs), p1 = toPlot([com[0] + R * a[0], com[1] + R * a[1], com[2] + R * a[2]], zs);
    traces.push({ type: "scatter3d", mode: "lines", name: "Polar-body axis", x: [p0[0], p1[0]], y: [p0[1], p1[1]], z: [p0[2], p1[2]], line: { color: AXIS_C, width: 5 }, hovertemplate: "polar-body axis<extra></extra>", legendrank: 40000 });
    Plotly.react("plot-host", traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  function renderReadout() {
    const e = cur(), g = gene(); if (!e) return;
    const rows = state.planes.map((p) => { const r = planeRec(e, g, p.key); return { p, r, a: r ? asym(r) : null, pv: r ? pOf(r) : null }; });
    $("#cp-readout").innerHTML = `<div class="cp-r-head"><b>${g}</b> · ${e.genes[g].nc} cell-body transcripts</div>` +
      rows.map((o) => o.r
        ? `<div class="cp-r-line"><span class="cp-sw" style="background:${o.p.color}"></span>${o.p.label}: asym <b>${o.a.toFixed(3)}</b> · <span class="${sig(o.pv) ? "cp-sig" : ""}">p ${fmtP(o.pv)}</span></div>`
        : `<div class="cp-r-line cp-r-na"><span class="cp-sw" style="background:${o.p.color}"></span>${o.p.label}: —</div>`).join("");
  }

  // ───────── bottom-drawer tabs ─────────
  const RENDER = { arrange: renderArrange, asym: renderAsym, pair: renderPair, heat: renderHeat, orient: renderOrient };
  function renderActive() { $("#drawer-gene").textContent = `· ${gene()}`; (RENDER[state.tab] || renderArrange)(); }
  const shown = (el) => !!(el && el.offsetParent);
  function plotInto(div, traces, layout) { Plotly.react(div, traces, layout, { responsive: true, displaylogo: false, displayModeBar: false }); }
  /** Prepare a panel for plotInto().
   *
   *  Plotly.react() PATCHES the existing graph DOM rather than rebuilding it. Blanking
   *  innerHTML while gd._fullLayout is still attached therefore makes the next react() a silent
   *  no-op — the first render works and every re-render comes back empty. So: leave a live graph
   *  alone (react updates it in place, which is also why it animates), and only clear the div
   *  when it holds a placeholder. The "no data" branches call Plotly.purge() first, which drops
   *  _fullLayout, so their placeholder is correctly cleared here on the way back. */
  function resetPlot(div) {
    if (div._fullLayout && div.classList.contains("js-plotly-plot")) return;
    div.innerHTML = "";
  }
  const baseLayout = (xt, yt) => ({ margin: { l: 54, r: 12, t: 10, b: 46 }, showlegend: false, paper_bgcolor: "transparent", plot_bgcolor: "#fcfdfe",
    xaxis: { title: { text: xt, font: { size: 11 } }, gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
    yaxis: { title: { text: yt, font: { size: 11 } }, gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } }, font: { size: 11, color: "#334155" } });
  // embryos that carry the current gene, with all four plane records resolved
  function geneRows(g) {
    return state.data.embryos.map((e) => {
      if (!(g in e.genes)) return null;
      const recs = {}; state.planes.forEach((p) => (recs[p.key] = planeRec(e, g, p.key)));
      return { e, recs };
    }).filter(Boolean);
  }
  const meanPairAngle = (recs, keys) => { const ns = keys.map((k) => recs[k]).filter(Boolean).map((r) => r.n); let s = 0, c = 0; for (let i = 0; i < ns.length; i++) for (let j = i + 1; j < ns.length; j++) { s += angleDeg(ns[i], ns[j]); c++; } return c ? s / c : null; };

  function renderArrange() {
    const div = $("#cp-arrange"); if (!shown(div)) return;
    const g = gene(), rows = geneRows(g);
    $("#cp-arrange-sub").textContent = `· ${g} · ${rows.length} zygotes`;
    if (rows.length < 2) { Plotly.purge(div); div.innerHTML = `<div class="cp-empty">Need ≥2 zygotes carrying ${g}.</div>`; $("#cp-arrange-note").textContent = ""; return; }
    resetPlot(div);
    // per embryo: each plane's angle to the polar-body axis (0–90°), one series per plane
    const order = rows.map((o) => ({ o, m: meanPairAngle(o.recs, state.planes.map((p) => p.key)) || 0 })).sort((a, b) => a.m - b.m).map((x) => x.o);
    const labels = order.map((o) => o.e.label);
    const traces = state.planes.map((p) => ({
      type: "scatter", mode: "markers", name: p.label,
      x: labels, y: order.map((o) => { const r = o.recs[p.key]; return r ? angleDeg(r.n, o.e.axis_um) : null; }),
      marker: { size: 9, color: p.color, line: { color: "#fff", width: 1 } },
      hovertemplate: `${p.label} · %{x} · %{y:.0f}° from polar axis<extra></extra>`,
    }));
    const lay = baseLayout("", "angle to polar-body axis (°)"); lay.xaxis.tickangle = -45; lay.xaxis.tickfont = { size: 8 }; lay.yaxis.range = [-3, 93];
    // Legend sits ABOVE the axes, so the top margin has to make room for it — baseLayout's
    // default t:10 leaves none and the keys land on top of the data.
    lay.showlegend = true; lay.margin.t = 46;
    lay.legend = { orientation: "h", y: 1.03, x: 0, yanchor: "bottom", font: { size: 9 } };
    plotInto(div, traces, lay);
    // summary: mean pairwise angle among all 4 planes across embryos
    const mps = rows.map((o) => meanPairAngle(o.recs, state.planes.map((p) => p.key))).filter((v) => v != null);
    const mean = mps.reduce((a, b) => a + b, 0) / (mps.length || 1);
    $("#cp-arrange-note").innerHTML = `Each dot is one pathway's plane in one zygote, placed by its angle to the polar-body axis (0° = the plane contains the axis). Embryos sorted by how tightly the four pathways agree. <b>Mean pairwise angle between the four planes ≈ ${mean.toFixed(0)}°</b> (0° = identical orientation, 90° = orthogonal).`;
  }

  function renderAsym() {
    const div = $("#cp-asym"); if (!shown(div)) return;
    const e = cur(), g = gene();
    $("#cp-asym-sub").textContent = `· ${g} · ${e.label}`;
    const rows = state.planes.map((p) => ({ p, r: planeRec(e, g, p.key) })).filter((o) => o.r);
    if (!rows.length) { Plotly.purge(div); div.innerHTML = `<div class="cp-empty">${g} has no plane data in ${e.label}.</div>`; return; }
    resetPlot(div);
    const tr = { type: "bar", x: rows.map((o) => o.p.label), y: rows.map((o) => asym(o.r)), marker: { color: rows.map((o) => o.p.color) },
      text: rows.map((o) => { const pv = pOf(o.r); return `p ${fmtP(pv)}${sig(pv) ? " *" : ""}`; }), textposition: "outside", textfont: { size: 10 },
      customdata: rows.map((o) => { const [a, b] = conc(o.r); return [o.r.nA, o.r.nB, pOf(o.r)]; }),
      hovertemplate: "%{x}<br>asymmetry %{y:.3f}<br>%{customdata[0]}/%{customdata[1]} · p %{customdata[2]:.2g}<extra></extra>" };
    const lay = baseLayout("", "concentration asymmetry"); lay.yaxis.range = [0, Math.max(0.05, Math.max(...rows.map((o) => asym(o.r))) * 1.25)];
    plotInto(div, [tr], lay);
  }

  function fillPairSelects() {
    const opts = state.planes.map((p) => `<option value="${p.key}">${p.label}</option>`).join("");
    ["#pair-a", "#pair-b"].forEach((s, i) => { const el = $(s); if (el) { el.innerHTML = opts; el.value = i ? state.pairB : state.pairA; } });
  }
  function renderPair() {
    const div = $("#cp-pair"); if (!shown(div)) return;
    const g = gene(), A = state.pairA, B = state.pairB;
    const rows = geneRows(g).filter((o) => o.recs[A] && o.recs[B]);
    $("#cp-pair-sub").textContent = `· ${g} · ${state.labelOf[A]} vs ${state.labelOf[B]} · ${rows.length} zygotes`;
    if (!rows.length) { Plotly.purge(div); div.innerHTML = `<div class="cp-empty">No zygote has both planes for ${g}.</div>`; $("#cp-pair-note").textContent = ""; return; }
    resetPlot(div);
    const data = rows.map((o) => ({ label: o.e.label, ang: angleDeg(o.recs[A].n, o.recs[B].n), aA: asym(o.recs[A]), aB: asym(o.recs[B]), cur: o.e.id === state.currentId }))
      .sort((a, b) => a.ang - b.ang);
    const tr = { type: "bar", x: data.map((d) => d.label), y: data.map((d) => d.ang),
      marker: { color: data.map((d) => d.cur ? "#0f172a" : "#64748b") },
      hovertemplate: "%{x}<br>angle %{y:.0f}° · alignment %{customdata:.2f}<extra></extra>", customdata: data.map((d) => 1 - d.ang / 90) };
    const lay = baseLayout("", `angle between the two planes (°)`); lay.xaxis.tickangle = -45; lay.xaxis.tickfont = { size: 8 }; lay.yaxis.range = [0, 93];
    lay.shapes = [{ type: "line", x0: -0.5, x1: data.length - 0.5, y0: 45, y1: 45, line: { color: "#cbd5e1", width: 1, dash: "dot" } }];
    plotInto(div, [tr], lay);
    const mAng = data.reduce((s, d) => s + d.ang, 0) / data.length, mAlign = 1 - mAng / 90;
    const mAsymA = data.reduce((s, d) => s + d.aA, 0) / data.length, mAsymB = data.reduce((s, d) => s + d.aB, 0) / data.length;
    $("#cp-pair-note").innerHTML = `Mean angle <b>${mAng.toFixed(0)}°</b> → alignment score <b>${mAlign.toFixed(2)}</b> (1 = identical orientation, 0 = orthogonal). ` +
      `Average concentration asymmetry created: <b style="color:${state.colorOf[A]}">${state.labelOf[A]} ${mAsymA.toFixed(3)}</b> · <b style="color:${state.colorOf[B]}">${state.labelOf[B]} ${mAsymB.toFixed(3)}</b>.`;
  }

  function renderHeat() {
    const div = $("#cp-heat"); if (!shown(div)) return;
    const g = gene(), rows = geneRows(g);
    $("#cp-heat-sub").textContent = `· ${g} · ${rows.length} zygotes`;
    if (!rows.length) { Plotly.purge(div); div.innerHTML = `<div class="cp-empty">No zygote carries ${g}.</div>`; return; }
    resetPlot(div);
    const order = rows.map((o) => ({ o, m: meanPairAngle(o.recs, state.planes.map((p) => p.key)) || 0 })).sort((a, b) => a.m - b.m).map((x) => x.o);
    const z = state.planes.map((p) => order.map((o) => { const r = o.recs[p.key]; if (!r) return null; return -Math.log10(Math.max(pOf(r), 1e-6)); }));
    const tr = { type: "heatmap", z, x: order.map((o) => o.e.label), y: state.planes.map((p) => p.label),
      colorscale: [[0, "#f8fafc"], [0.25, "#c7d2fe"], [0.5, "#818cf8"], [0.75, "#6d28d9"], [1, "#3b0764"]], zmin: 0, zmax: 4,
      colorbar: { title: { text: "−log₁₀ p", side: "right", font: { size: 9 } }, thickness: 10, len: 0.8, tickfont: { size: 8 } },
      hovertemplate: "%{y} · %{x}<br>−log₁₀ p %{z:.1f}<extra></extra>", xgap: 1, ygap: 1 };
    const lay = baseLayout("", ""); lay.margin.l = 110; lay.margin.b = 70; lay.xaxis.tickangle = -45; lay.xaxis.tickfont = { size: 8 };
    plotInto(div, [tr], lay);
  }

  function renderOrient() {
    const div = $("#cp-orient"); if (!shown(div)) return;
    const g = gene(), rows = geneRows(g);
    $("#cp-orient-sub").textContent = `· ${g} · ${rows.length} zygotes`;
    if (rows.length < 2) { Plotly.purge(div); div.innerHTML = `<div class="cp-empty">Need ≥2 zygotes carrying ${g}.</div>`; return; }
    resetPlot(div);
    // box/strip per plane of angle-to-polar-axis
    const box = state.planes.map((p) => ({
      type: "box", name: p.label, boxpoints: "all", jitter: 0.5, pointpos: 0, marker: { size: 6, color: p.color }, line: { color: p.color }, fillcolor: p.color + "22",
      y: rows.map((o) => { const r = o.recs[p.key]; return r ? angleDeg(r.n, o.e.axis_um) : null; }).filter((v) => v != null),
      x: rows.map((o) => (o.recs[p.key] ? p.label : null)).filter(Boolean),
      hovertemplate: `${p.label} · %{y:.0f}°<extra></extra>`,
    }));
    const lay = baseLayout("", "angle to polar-body axis (°)"); lay.yaxis.range = [-3, 93]; lay.boxmode = "group";
    plotInto(div, box, lay);
  }
  function rowsAxis() { return [0, 0, 1]; } // unused fallback

  // ───────── right drawer: gene ranking by cross-plane agreement ─────────
  function buildRankPlanes() {
    const host = $("#cp-rank-planes"); if (!host) return;
    host.innerHTML = state.planes.map((p) => `<label class="zt"><input type="checkbox" data-plane="${p.key}" ${state.rankPlanes[p.key] ? "checked" : ""}><span class="zt-sw" style="background:${p.color}"></span>${p.label}</label>`).join("");
  }
  function rankGenes() {
    const keys = state.planes.map((p) => p.key).filter((k) => state.rankPlanes[k]);
    const cacheKey = `${keys.join(",")}|${state.minCount}`;
    if (state._rankKey === cacheKey) return state._rankRows;
    const per = {};   // gene -> {sum, n, tot}
    state.data.embryos.forEach((e) => {
      for (const g in e.genes) {
        const recs = {}; let ok = 0; keys.forEach((k) => { const r = planeRec(e, g, k); if (r) { recs[k] = r; ok++; } });
        if (ok < 2) continue;
        const m = meanPairAngle(recs, keys); if (m == null) continue;
        const a = per[g] || (per[g] = { sum: 0, n: 0, tot: 0 });
        a.sum += m; a.n++; a.tot += e.genes[g].nc;
      }
    });
    const rows = Object.entries(per).map(([g, a]) => ({ gene: g, ang: a.sum / a.n, n: a.n, tot: a.tot })).filter((r) => r.tot >= state.minCount);
    state._rankKey = cacheKey; state._rankRows = rows; return rows;
  }
  function renderRank() {
    const host = $("#cp-rank-list"); if (!host) return;
    const keys = state.planes.map((p) => p.key).filter((k) => state.rankPlanes[k]);
    if (keys.length < 2) { host.innerHTML = `<div class="best-plane-note">Pick at least two plane types above.</div>`; return; }
    const rows = rankGenes().slice().sort((a, b) => state.rankMode === "agree" ? a.ang - b.ang : b.ang - a.ang);
    const curG = gene();
    let html = `<div class="best-plane-note">Ranked by <b>${state.rankMode === "agree" ? "agreement (smallest angle)" : "disagreement (largest angle)"}</b> among ${keys.length} planes, mean over each gene's zygotes (≥ ${state.minCount} transcripts).</div>`;
    html += `<div class="best-head cp-best-head"><span></span><span>gene</span><span title="mean pairwise angle between the chosen planes (°)">angle</span><span title="zygotes averaged">zyg</span></div>`;
    html += rows.slice(0, 200).map((r, i) =>
      `<div class="best-row${r.gene === curG ? " current" : ""}" data-gene="${r.gene}" title="${r.gene} · mean ${r.ang.toFixed(1)}° over ${r.n} zygotes · ${r.tot} transcripts">` +
      `<span class="best-num">${i + 1}</span><span class="best-gene">${r.gene}</span>` +
      `<span class="best-real">${r.ang.toFixed(1)}°</span><span class="best-p">${r.n}</span></div>`).join("");
    host.innerHTML = rows.length ? html : `<div class="best-plane-note">No gene has ≥2 of the chosen planes in ≥1 zygote at this threshold.</div>`;
  }

  // ───────── wiring ─────────
  function setGene(g) { state.gene = g; $("#gene-select").value = g; render3D(); renderReadout(); renderActive(); renderRank(); }
  function wire() {
    $("#gene-select").addEventListener("change", (e) => setGene(e.target.value));
    $("#plane-toggles").addEventListener("change", (e) => { const c = e.target.closest("input[data-plane]"); if (!c) return; state.planesOn[c.dataset.plane] = c.checked; render3D(); });
    $("#dots-show").addEventListener("change", (e) => { state.dotsOn = e.target.checked; render3D(); });
    $("#split-by").addEventListener("change", (e) => { state.splitBy = e.target.value; render3D(); });
    $("#cp-tabs").addEventListener("click", (e) => { const t = e.target.closest(".xs-gtab"); if (t) switchTab(t.dataset.tab); });
    $("#pair-a").addEventListener("change", (e) => { state.pairA = e.target.value; renderPair(); });
    $("#pair-b").addEventListener("change", (e) => { state.pairB = e.target.value; renderPair(); });
    $("#cp-rank-planes").addEventListener("change", (e) => { const c = e.target.closest("input[data-plane]"); if (!c) return; state.rankPlanes[c.dataset.plane] = c.checked; state._rankKey = null; renderRank(); });
    $("#cp-rank-mode").addEventListener("click", (e) => { const b = e.target.closest("[data-mode]"); if (!b) return; state.rankMode = b.dataset.mode; $("#cp-rank-mode").querySelectorAll("button").forEach((x) => x.classList.toggle("active", x === b)); renderRank(); });
    $("#cp-rank-min").addEventListener("input", (e) => { state.minCount = Math.max(1, parseInt(e.target.value, 10) || 1); state._rankKey = null; renderRank(); });
    $("#cp-rank-list").addEventListener("click", (e) => { const row = e.target.closest(".best-row"); if (row && row.dataset.gene && (row.dataset.gene in cur().genes)) setGene(row.dataset.gene); });
    wireDrawer(); wireRdrawer();
    renderRank();
  }
  function switchTab(which) {
    state.tab = which;
    $("#cp-tabs").querySelectorAll(".xs-gtab").forEach((t) => { const on = t.dataset.tab === which; t.classList.toggle("active", on); t.setAttribute("aria-selected", String(on)); });
    $("#cp-panels").querySelectorAll(".xs-panel").forEach((p) => { p.hidden = p.dataset.tab !== which; });
    renderActive(); requestAnimationFrame(() => { try { Plotly.Plots.resize($("#cp-" + (which === "asym" ? "asym" : which === "pair" ? "pair" : which === "heat" ? "heat" : which === "orient" ? "orient" : "arrange"))); } catch (_) {} });
  }
  function openDrawer(open) { state.drawerOpen = open; $("#drawer").dataset.open = open ? "true" : "false"; $("#drawer-handle").setAttribute("aria-expanded", String(open)); if (open) renderActive(); }
  function wireDrawer() {
    $("#drawer-handle").addEventListener("click", () => openDrawer($("#drawer").dataset.open !== "true"));
    const rz = $("#drawer-resize"); let sh = 0;
    rz.addEventListener("pointerdown", (e) => { sh = $("#drawer-body").getBoundingClientRect().height; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return; $("#drawer").style.setProperty("--drawer-h", Math.max(200, Math.min(window.innerHeight - 120, sh + (rz._d.y - e.clientY))) + "px"); });
    const end = (e) => { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} }; rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }
  function wireRdrawer() {
    const rd = $("#rdrawer"), h = $("#rdrawer-handle");
    h.addEventListener("click", () => { const o = rd.dataset.open !== "true"; rd.dataset.open = o ? "true" : "false"; h.setAttribute("aria-expanded", String(o)); if (o) renderRank(); });
    const rrz = $("#rdrawer-resize");
    rrz.addEventListener("pointerdown", (e) => { rrz._d = { x: e.clientX, w: rd.getBoundingClientRect().width }; rrz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rrz.addEventListener("pointermove", (e) => { if (!rrz._d) return; rd.style.setProperty("--rdrawer-w", Math.max(260, Math.min(window.innerWidth - 80, rrz._d.w - (e.clientX - rrz._d.x))) + "px"); });
    const end = (e) => { rrz._d = null; try { rrz.releasePointerCapture(e.pointerId); } catch (_) {} }; rrz.addEventListener("pointerup", end); rrz.addEventListener("pointercancel", end);
  }
})();
