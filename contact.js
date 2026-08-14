/* Blastomere Contact Enrichment.
 *
 * Is anything concentrated at the face where the two blastomeres press together? Data:
 * data/contact.json.gz (build_contact.py) ships, per embryo and per slab thickness D, the
 * fraction of ALL that embryo's blastomere transcripts within |t| <= D of the interface (f0)
 * and each gene's count there (k) out of n. Every statistic below is derived from those.
 *
 * The null is deliberately "distributed like the rest of the transcriptome in this same
 * embryo", not "uniform in space" — that is what makes a fold of 1.3 mean something, since
 * cell shape and detection efficiency divide out.
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const CT = "#7c3aed", GENE = "#2563eb", E2C = "#7c3aed", L2C = "#f97316";

  const state = {
    data: null, byId: {}, ri: 5, gene: null, currentId: null, scene: null, sceneCache: {},
    tab: "one", slabOn: true, dotsOn: true, axisOn: false, dotSize: V.DOT_SIZE,
    genesStage: "both", genesMin: 300, rankStage: "both", rankDir: "enr", rankMin: 300,
    drawerOpen: false, vcExtras: null,
  };

  // ───────── binomial statistics (same helpers the other projects use) ─────────
  function gammaln(x) { const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, t = x + 5.5; t -= (x + 0.5) * Math.log(t); let s = 1.000000000190015;
    for (let j = 0; j < 6; j++) s += c[j] / ++y;
    return -t + Math.log(2.5066282746310005 * s / x); }
  function betacf(a, b, x) { const FP = 1e-300; let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FP) d = FP; d = 1 / d; let h = d;
    for (let m = 1; m <= 200; m++) { const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2)); d = 1 + aa * d; if (Math.abs(d) < FP) d = FP;
      c = 1 + aa / c; if (Math.abs(c) < FP) c = FP; d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)); d = 1 + aa * d; if (Math.abs(d) < FP) d = FP;
      c = 1 + aa / c; if (Math.abs(c) < FP) c = FP; d = 1 / d; const del = d * c; h *= del;
      if (Math.abs(del - 1) < 3e-12) break; } return h; }
  function betai(a, b, x) { if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(gammaln(a + b) - gammaln(a) - gammaln(b) + a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b; }
  const binomSF = (k, n, p) => (k <= 0 ? 1 : k > n ? 0 : betai(k, n - k + 1, p));
  const binomCDF = (k, n, p) => (k >= n ? 1 : k < 0 ? 0 : 1 - betai(k + 1, n - k, p));
  function binomTwo(k, n, p) { if (!n) return 1;
    return Math.min(1, 2 * Math.min(binomCDF(k, n, p), binomSF(k, n, p))); }
  const fmtP = (p) => (p == null || !isFinite(p)) ? "–" : p < 1e-4 ? p.toExponential(1) : p < 0.1 ? p.toPrecision(2) : p.toFixed(2);
  const sig = (p) => p != null && p <= 0.05;
  /** 95% range of k/n under the null, as a fold. */
  function nullBand(n, p) {
    const sd = Math.sqrt(Math.max(p * (1 - p) / Math.max(n, 1), 0));
    return [Math.max(0, (p - 1.96 * sd) / p), (p + 1.96 * sd) / p];
  }

  const D = () => state.data.meta.radii[state.ri];
  const cur = () => state.byId[state.currentId];
  const gene = () => state.gene;
  /** {k,n,f0,fold,p} for one gene in one embryo at the current slab, or null. */
  function rec(e, g) {
    const r = e && e.genes[g];
    if (!r) return null;
    const f0 = e.f0[state.ri], k = r.k[state.ri], n = r.n;
    return { k, n, f0, fold: f0 > 0 ? (k / n) / f0 : null, p: binomTwo(k, n, f0) };
  }

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await V.loadGz("data/contact.json.gz"); }
    catch (e) {
      $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_contact.py</code>.</div></div>`;
      return;
    }
    const d = state.data, m = d.meta;
    d.embryos.forEach((e) => (state.byId[e.id] = e));
    state.ri = m.radii.indexOf(m.default_radius);
    if (state.ri < 0) state.ri = Math.floor(m.radii.length / 2);
    $("#radius-range").max = m.radii.length - 1;
    $("#radius-range").value = state.ri;

    $("#embryo-count").textContent =
      `${m.n_embryos} two-cell embryos · ${m.n_e2c} early + ${m.n_l2c} late · ${m.n_genes} genes`;
    $("#gene-select").innerHTML = d.genes.map((g) => `<option value="${g}">${g}</option>`).join("");
    state.gene = d.genes.includes("Gpsm1") ? "Gpsm1" : d.genes[0];
    $("#gene-select").value = state.gene;

    V.buildTabs($("#tabs"), d.embryos, selectEmbryo, (e) => ({
      label: V.embryoLabel ? V.embryoLabel(e.id, e.stage) : e.id,
      sub: e.stage === "e2c" ? "early 2C" : "late 2C",
      title: `${e.id} · ${e.n_all.toLocaleString()} blastomere transcripts · ${e.n_genes} genes`,
    }));
    wire();
    selectEmbryo(d.embryos[0].id);
  })();

  async function selectEmbryo(id) {
    if (id === state.currentId && state.scene) return;
    state.currentId = id; V.markActiveTab($("#tabs"), id);
    const e = state.byId[id];
    $("#loading").hidden = false; $("#loading-text").textContent = `Loading ${id}…`;
    try {
      let sc = state.sceneCache[id];
      if (!sc) { sc = await V.loadGz(`data/segments/${e.scene}`); state.sceneCache[id] = sc; }
      if (state.currentId !== id) return;
      state.scene = sc;
      // keep the gene if this embryo has it, else fall back to its best-measured one
      if (!(state.gene in e.genes)) {
        const best = Object.keys(e.genes).sort((a, b) => e.genes[b].n - e.genes[a].n)[0];
        if (best) { state.gene = best; $("#gene-select").value = best; }
      }
      if (!state.vcExtras) state.vcExtras = V.addWindowExtras($("#controls-body"),
        { defaultSize: state.dotSize, onDotSize: (s) => { state.dotSize = s; render3D(); } });
      state.vcExtras.setAtlas && state.vcExtras.setAtlas(id);
      $("#controls").hidden = false; $("#placeholder").hidden = true;
      $("#drawer").hidden = false; $("#rdrawer").hidden = false;
      render3D(); renderReadout(); renderActive(); renderRank();
      if (!state.drawerOpen) openDrawer(true);
    } catch (err) {
      $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Scene missing</div>` +
        `<div class="ph-sub">${err.message || err}</div></div>`;
    } finally { $("#loading").hidden = true; }
  }

  // ───────── 3-D ─────────
  const toPlot = (pUm, zs) => [pUm[0] / XY, pUm[1] / XY, pUm[2] * zs];
  function unitv(n) { const m = Math.hypot(n[0], n[1], n[2]) || 1; return [n[0] / m, n[1] / m, n[2] / m]; }

  /** A square quad of the plane through `c` with normal `n`, in plot space. */
  function quad(cUm, nUm, half, zs, color, op, name, rank) {
    const nn = unitv(nUm), ref = Math.abs(nn[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const cr = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    const t = unitv(cr(nn, ref)), w = unitv(cr(nn, t));
    const P = [], sgn = [[-1, -1], [1, -1], [1, 1], [-1, 1]];
    for (const [a, b] of sgn) {
      P.push(toPlot([cUm[0] + half * (a * t[0] + b * w[0]), cUm[1] + half * (a * t[1] + b * w[1]),
                     cUm[2] + half * (a * t[2] + b * w[2])], zs));
    }
    return { type: "mesh3d", x: P.map((p) => p[0]), y: P.map((p) => p[1]), z: P.map((p) => p[2]),
      i: [0, 0], j: [1, 2], k: [2, 3], color, opacity: op, name, showlegend: true,
      hoverinfo: "skip", flatshading: true, legendrank: rank };
  }

  function render3D() {
    const s = state.scene, e = cur(), g = gene(); if (!s || !e) return;
    const zs = s.z_scale, traces = V.bodyTraces(s);
    const mid = e.mid_um, ax = e.axis_um, d = D();
    const half = e.sep_um * 0.85;

    if (state.dotsOn && s.transcripts[g]) {
      const t = s.transcripts[g], seg = t.s || null;
      const inx = [], iny = [], inz = [], ox = [], oy = [], oz = [];
      for (let i = 0; i < t.x.length; i++) {
        if (seg && seg[i] !== e.a && seg[i] !== e.b) continue;      // blastomeres only
        const pUm = [t.x[i] * XY, t.y[i] * XY, t.gz[i]];
        const tt = (pUm[0] - mid[0]) * ax[0] + (pUm[1] - mid[1]) * ax[1] + (pUm[2] - mid[2]) * ax[2];
        const P = [t.x[i], t.y[i], t.gz[i] * zs];
        if (Math.abs(tt) <= d) { inx.push(P[0]); iny.push(P[1]); inz.push(P[2]); }
        else { ox.push(P[0]); oy.push(P[1]); oz.push(P[2]); }
      }
      if (ox.length) traces.push({ type: "scatter3d", mode: "markers", name: `${g} · elsewhere`,
        x: ox, y: oy, z: oz, marker: { size: state.dotSize, color: "#94a3b8", opacity: .5, line: { width: 0 } },
        hovertemplate: `${g} · outside the slab<extra></extra>`, legendrank: 20001 });
      if (inx.length) traces.push({ type: "scatter3d", mode: "markers", name: `${g} · at the contact`,
        x: inx, y: iny, z: inz, marker: { size: state.dotSize + 0.8, color: GENE, opacity: V.DOT_OPACITY, line: { width: 0 } },
        hovertemplate: `${g} · in the contact slab<extra></extra>`, legendrank: 20000 });
    }

    if (state.slabOn) {
      // the two faces of the slab, so its thickness is visible rather than implied
      traces.push(quad([mid[0] + ax[0] * d, mid[1] + ax[1] * d, mid[2] + ax[2] * d], ax, half, zs,
        CT, 0.16, `Contact slab ±${d} µm`, 41000));
      traces.push(quad([mid[0] - ax[0] * d, mid[1] - ax[1] * d, mid[2] - ax[2] * d], ax, half, zs,
        CT, 0.16, "slab (far face)", 41001));
    }
    if (state.axisOn) {
      const a0 = e.com_a_um, b0 = e.com_b_um;
      const p0 = toPlot(a0, zs), p1 = toPlot(b0, zs);
      traces.push({ type: "scatter3d", mode: "lines+markers", name: "Interface axis",
        x: [p0[0], p1[0]], y: [p0[1], p1[1]], z: [p0[2], p1[2]],
        line: { color: "#111", width: 5 }, marker: { size: 5, color: "#111" },
        hovertemplate: "blastomere centre<extra></extra>", legendrank: 40000 });
    }
    Plotly.react($("#plot-host"), traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  function renderReadout() {
    const e = cur(), g = gene(); if (!e) return;
    const r = rec(e, g);
    const cls = r && r.fold != null ? (r.fold > 1 ? "up" : "dn") : "";
    $("#ct-readout").innerHTML =
      `<div class="ct-r-head"><b>${g}</b> · ${e.stage === "e2c" ? "early" : "late"} 2-cell</div>` +
      (r
        ? `<div class="ct-r-line"><span class="k">in the slab</span><span class="v">${r.k} / ${r.n}</span></div>` +
          `<div class="ct-r-line"><span class="k">expected</span><span class="v">${(r.f0 * 100).toFixed(1)}%</span></div>` +
          `<div class="ct-r-line"><span class="k">fold</span><span class="v ${cls}">${r.fold == null ? "–" : r.fold.toFixed(2)}×</span></div>` +
          `<div class="ct-r-line"><span class="k">p</span><span class="v">${fmtP(r.p)}</span></div>`
        : `<div class="ct-r-line"><span class="k">not measured in this embryo</span></div>`) +
      `<div class="ct-r-line"><span class="k">blastomere gap</span><span class="v">${e.sep_um} µm</span></div>`;
  }

  // ───────── drawer panels ─────────
  const shown = (el) => !!(el && el.offsetParent);
  function resetPlot(div) {
    if (div._fullLayout && div.classList.contains("js-plotly-plot")) return;
    div.innerHTML = "";
  }
  const baseLayout = (xt, yt) => ({
    margin: { l: 58, r: 12, t: 10, b: 46 }, showlegend: false,
    paper_bgcolor: "transparent", plot_bgcolor: "#fcfdfe",
    xaxis: { title: { text: xt, font: { size: 11 } }, gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
    yaxis: { title: { text: yt, font: { size: 11 } }, gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
    font: { size: 11, color: "#334155" },
  });

  /** Panel 1 — the axis profile for this embryo. */
  function renderOne() {
    const div = $("#ct-one"); if (!shown(div)) return;
    const s = state.scene, e = cur(), g = gene(); if (!s || !e) return;
    $("#ct-one-sub").textContent = `· ${g} · ${e.id.replace(/^\d{8}_/, "")}`;
    resetPlot(div);
    const mid = e.mid_um, ax = e.axis_um;
    const proj = (t, i) => (t.x[i] * XY - mid[0]) * ax[0] + (t.y[i] * XY - mid[1]) * ax[1] + (t.gz[i] - mid[2]) * ax[2];
    const geneT = [], allT = [];
    for (const [gn, t] of Object.entries(s.transcripts)) {
      const seg = t.s || null;
      for (let i = 0; i < t.x.length; i++) {
        if (seg && seg[i] !== e.a && seg[i] !== e.b) continue;
        const v = proj(t, i);
        allT.push(v);
        if (gn === g) geneT.push(v);
      }
    }
    const lim = Math.max(1, e.sep_um * 0.9);
    const bins = { start: -lim, end: lim, size: (2 * lim) / 40 };
    const traces = [
      { type: "histogram", name: "all transcripts", x: allT, xbins: bins, histnorm: "probability",
        marker: { color: "#cbd5e1" }, hovertemplate: "all · %{x:.1f} µm: %{y:.3f}<extra></extra>" },
      { type: "histogram", name: g, x: geneT, xbins: bins, histnorm: "probability",
        marker: { color: GENE, opacity: .75 }, hovertemplate: `${g} · %{x:.1f} µm: %{y:.3f}<extra></extra>` },
    ];
    const lay = baseLayout("position along the interface axis (µm) — 0 = the contact", "fraction of transcripts");
    lay.barmode = "overlay"; lay.showlegend = true; lay.margin.t = 34;
    lay.legend = { orientation: "h", y: 1.03, x: 0, yanchor: "bottom", font: { size: 9 } };
    lay.shapes = [{ type: "rect", x0: -D(), x1: D(), yref: "paper", y0: 0, y1: 1,
      fillcolor: CT, opacity: .10, line: { width: 0 } },
      { type: "line", x0: 0, x1: 0, yref: "paper", y0: 0, y1: 1,
        line: { color: CT, width: 1.5, dash: "dot" } }];
    Plotly.react(div, traces, lay, { responsive: true, displaylogo: false, displayModeBar: false });

    const r = rec(e, g);
    $("#ct-one-stat").innerHTML = r
      ? `<b>${r.k}</b> of <b>${r.n}</b> ${g} transcripts lie within ±${D()} µm of the interface — ` +
        `<b>${((r.k / r.n) * 100).toFixed(1)}%</b> against an expected <b>${(r.f0 * 100).toFixed(1)}%</b>, ` +
        `fold <b>${r.fold == null ? "–" : r.fold.toFixed(2)}×</b>, ` +
        `<span class="${sig(r.p) ? "sig" : ""}">p ${fmtP(r.p)}</span>.`
      : `${g} has fewer than ${state.data.meta.min_tx} blastomere transcripts here.`;
  }

  /** Panels 2-4 — this gene across a set of embryos. */
  function acrossPanel(divId, subId, stages, statId) {
    const div = $("#" + divId); if (!shown(div)) return;
    const g = gene();
    const rows = state.data.embryos
      .filter((e) => stages.includes(e.stage) && e.genes[g])
      .map((e) => ({ e, r: rec(e, g) }))
      .filter((o) => o.r && o.r.fold != null)
      .sort((a, b) => b.r.fold - a.r.fold);
    $("#" + subId).textContent = `· ${g} · ±${D()} µm · ${rows.length} embryos`;
    resetPlot(div);
    if (!rows.length) {
      Plotly.purge(div);
      div.innerHTML = `<div class="ct-empty">${g} is not measured in any of these embryos.</div>`;
      if (statId) $("#" + statId).textContent = "";
      return;
    }
    const bands = rows.map((o) => nullBand(o.r.n, o.r.f0));
    const traces = [{
      type: "bar", x: rows.map((o) => V.embryoLabel ? V.embryoLabel(o.e.id, o.e.stage) : o.e.id),
      y: rows.map((o) => o.r.fold),
      marker: { color: rows.map((o) => (stages.length > 1
        ? (o.e.stage === "e2c" ? E2C : L2C)
        : (sig(o.r.p) ? (o.r.fold > 1 ? "#16a34a" : "#dc2626") : "#cbd5e1"))) },
      error_y: { type: "data", symmetric: false,
        array: bands.map((b, i) => b[1] - rows[i].r.fold),
        arrayminus: bands.map((b, i) => rows[i].r.fold - b[0]),
        color: "#94a3b8", thickness: 1, width: 2 },
      customdata: rows.map((o) => [o.r.k, o.r.n, fmtP(o.r.p), o.e.stage]),
      hovertemplate: "<b>%{x}</b> (%{customdata[3]})<br>fold %{y:.2f}×<br>" +
        "%{customdata[0]}/%{customdata[1]} in the slab · p %{customdata[2]}<extra></extra>",
    }];
    const lay = baseLayout("", "fold enrichment at the contact");
    lay.xaxis.tickangle = -45; lay.xaxis.tickfont = { size: 8 };
    lay.shapes = [{ type: "line", x0: -0.5, x1: rows.length - 0.5, y0: 1, y1: 1,
      line: { color: "#94a3b8", width: 1, dash: "dot" } }];
    Plotly.react(div, traces, lay, { responsive: true, displaylogo: false, displayModeBar: false });

    if (statId) {
      const k = rows.reduce((a, o) => a + o.r.k, 0);
      const exp = rows.reduce((a, o) => a + o.r.n * o.r.f0, 0);
      const n = rows.reduce((a, o) => a + o.r.n, 0);
      const fold = exp > 0 ? k / exp : null;
      const p = binomTwo(k, n, exp / n);
      const nE = rows.filter((o) => o.e.stage === "e2c").length;
      $("#" + statId).innerHTML =
        `Pooled over ${rows.length} embryos (${nE} early, ${rows.length - nE} late): ` +
        `<b>${k}</b> of <b>${n}</b> against an expected <b>${exp.toFixed(0)}</b> — ` +
        `fold <b>${fold == null ? "–" : fold.toFixed(2)}×</b>, ` +
        `<span class="${sig(p) ? "sig" : ""}">p ${fmtP(p)}</span>.`;
    }
  }

  /** Pooled fold per gene over a stage selection. */
  function poolGenes(stage, min) {
    const acc = {};
    for (const e of state.data.embryos) {
      if (stage !== "both" && e.stage !== stage) continue;
      const f0 = e.f0[state.ri];
      for (const [g, r] of Object.entries(e.genes)) {
        const a = acc[g] || (acc[g] = { k: 0, n: 0, exp: 0, emb: 0 });
        a.k += r.k[state.ri]; a.n += r.n; a.exp += r.n * f0; a.emb++;
      }
    }
    return Object.entries(acc)
      .filter(([, a]) => a.n >= min && a.exp > 0)
      .map(([g, a]) => ({ gene: g, fold: a.k / a.exp, k: a.k, n: a.n,
                          p: binomTwo(a.k, a.n, a.exp / a.n), emb: a.emb }));
  }

  /** Panel 5 — every gene compared. */
  function renderGenes() {
    const div = $("#ct-genes"); if (!shown(div)) return;
    const rows = poolGenes(state.genesStage, state.genesMin).sort((a, b) => b.fold - a.fold);
    const label = { both: "both 2-cell stages", e2c: "early 2-cell", l2c: "late 2-cell" }[state.genesStage];
    $("#ct-genes-sub").textContent = `· ${label} · ±${D()} µm · ${rows.length} genes`;
    resetPlot(div);
    if (!rows.length) {
      Plotly.purge(div);
      div.innerHTML = `<div class="ct-empty">No gene reaches ${state.genesMin} transcripts here.</div>`;
      return;
    }
    const top = rows.slice(0, 30).concat(rows.slice(-10)).filter((v, i, a) => a.indexOf(v) === i);
    top.sort((a, b) => a.fold - b.fold);
    const traces = [{
      type: "bar", orientation: "h", x: top.map((r) => r.fold), y: top.map((r) => r.gene),
      marker: { color: top.map((r) => (sig(r.p) ? (r.fold > 1 ? "#16a34a" : "#dc2626") : "#cbd5e1")) },
      customdata: top.map((r) => [r.k, r.n, fmtP(r.p), r.emb]),
      hovertemplate: "<b>%{y}</b><br>fold %{x:.2f}×<br>%{customdata[0]}/%{customdata[1]} · " +
        "%{customdata[3]} embryos · p %{customdata[2]}<extra></extra>",
    }];
    const lay = baseLayout("fold enrichment at the contact", "");
    lay.margin.l = 92; lay.yaxis.tickfont = { size: 9 };
    lay.shapes = [{ type: "line", x0: 1, x1: 1, yref: "paper", y0: 0, y1: 1,
      line: { color: "#94a3b8", width: 1, dash: "dot" } }];
    Plotly.react(div, traces, lay, { responsive: true, displaylogo: false, displayModeBar: false });
    div.on("plotly_click", (ev) => {
      const y = ev.points && ev.points[0] && ev.points[0].y;
      if (y) setGene(y);
    });
  }

  // ───────── right drawer ─────────
  function renderRank() {
    const el = $("#ct-rank-list");
    let rows = poolGenes(state.rankStage, state.rankMin);
    rows.sort((a, b) => (state.rankDir === "enr" ? b.fold - a.fold : a.fold - b.fold));
    rows = rows.slice(0, 120);
    if (!rows.length) { el.innerHTML = `<div class="ct-empty">Nothing clears ${state.rankMin} transcripts.</div>`; return; }
    el.innerHTML = `<div class="ct-head"><span></span><span>gene</span><span>fold</span><span>p</span></div>` +
      rows.map((r, i) =>
        `<div class="ct-row${r.gene === state.gene ? " current" : ""}" data-gene="${r.gene}" ` +
        `title="${r.gene} · ${r.k}/${r.n} in the slab across ${r.emb} embryos">` +
        `<span class="n">${i + 1}</span><span class="g">${r.gene}</span>` +
        `<span class="f">${r.fold.toFixed(2)}×</span>` +
        `<span class="p${sig(r.p) ? " sig" : ""}">${fmtP(r.p)}</span></div>`).join("");
    el.querySelectorAll(".ct-row").forEach((n) => n.addEventListener("click", () => setGene(n.dataset.gene)));
  }

  // ───────── plumbing ─────────
  const RENDER = {
    one: renderOne,
    e2c: () => acrossPanel("ct-e2c", "ct-e2c-sub", ["e2c"]),
    l2c: () => acrossPanel("ct-l2c", "ct-l2c-sub", ["l2c"]),
    all: () => acrossPanel("ct-all", "ct-all-sub", ["e2c", "l2c"], "ct-all-stat"),
    genes: renderGenes,
  };
  function renderActive() { $("#drawer-gene").textContent = `· ${gene()}`; (RENDER[state.tab] || renderOne)(); }

  function setGene(g) {
    if (!g) return;
    state.gene = g; $("#gene-select").value = g;
    render3D(); renderReadout(); renderActive(); renderRank();
  }

  function openDrawer(open) {
    state.drawerOpen = open;
    $("#drawer").dataset.open = open ? "true" : "false";
    $("#drawer-handle").setAttribute("aria-expanded", String(open));
    if (open) renderActive();
  }

  function wire() {
    $("#gene-select").addEventListener("change", (e) => setGene(e.target.value));
    $("#radius-range").addEventListener("input", (e) => {
      state.ri = +e.target.value;
      $("#radius-val").textContent = `± ${D()} µm`;
      render3D(); renderReadout(); renderActive(); renderRank();
    });
    $("#slab-show").addEventListener("change", (e) => { state.slabOn = e.target.checked; render3D(); });
    $("#dots-show").addEventListener("change", (e) => { state.dotsOn = e.target.checked; render3D(); });
    $("#axis-show").addEventListener("change", (e) => { state.axisOn = e.target.checked; render3D(); });

    $("#genes-stage").addEventListener("change", (e) => { state.genesStage = e.target.value; renderGenes(); });
    $("#genes-min").addEventListener("change", (e) => {
      state.genesMin = Math.max(20, +e.target.value || 300); renderGenes();
    });

    $("#drawer-handle").addEventListener("click", () => openDrawer($("#drawer").dataset.open !== "true"));
    $("#ct-tabs").addEventListener("click", (e) => {
      const b = e.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#ct-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#ct-panels").querySelectorAll(".xs-panel").forEach((p) => (p.hidden = p.dataset.tab !== state.tab));
      renderActive();
    });
    const rz = $("#drawer-resize"); let sh = 0;
    rz.addEventListener("pointerdown", (e) => {
      sh = $("#drawer-body").getBoundingClientRect().height; rz._d = { y: e.clientY };
      rz.setPointerCapture(e.pointerId); e.preventDefault();
    });
    rz.addEventListener("pointermove", (e) => {
      if (!rz._d) return;
      $("#drawer").style.setProperty("--drawer-h",
        Math.max(200, Math.min(window.innerHeight - 120, sh + (rz._d.y - e.clientY))) + "px");
    });
    const end = (e) => { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);

    const rd = $("#rdrawer"), rh = $("#rdrawer-handle");
    rh.addEventListener("click", () => {
      const o = rd.dataset.open !== "true";
      rd.dataset.open = o ? "true" : "false";
      rh.setAttribute("aria-expanded", String(o));
      if (o) renderRank();
    });
    $("#ct-rank-stage").addEventListener("click", (e) => {
      const b = e.target.closest("button"); if (!b) return;
      state.rankStage = b.dataset.stage;
      [...e.currentTarget.children].forEach((x) => x.classList.toggle("active", x === b));
      renderRank();
    });
    $("#ct-rank-dir").addEventListener("click", (e) => {
      const b = e.target.closest("button"); if (!b) return;
      state.rankDir = b.dataset.dir;
      [...e.currentTarget.children].forEach((x) => x.classList.toggle("active", x === b));
      renderRank();
    });
    $("#rank-min").addEventListener("change", (e) => {
      state.rankMin = Math.max(20, +e.target.value || 300); renderRank();
    });

    window.addEventListener("resize", () => { try { Plotly.Plots.resize($("#plot-host")); } catch (_) {} });
  }
})();
