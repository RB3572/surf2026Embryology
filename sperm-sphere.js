/* Sperm-Entry-Site Enrichment — for zygotes with a labelled sperm, test whether each
 * gene's transcripts are CONCENTRATED (enriched) or DEPLETED within a sphere (radius set
 * in the floating window) around the sperm entry site. The build ships raw per-SEGMENT
 * transcript counts + voxel volumes (cytoplasm / pronuclei / polar body / other); this
 * script sums whichever segments are toggled on and computes the concentration, fold and a
 * binomial null (p + 95% band) live, so every chart updates as you include/exclude segments.
 * 3-D: the cell, the sperm dot, a sphere at the site, and the gene's molecules in/out of it.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const SPERM_C = "#db2777", SPHERE_C = "#f472b6";
  const IN_C = "#db2777", OUT_C = "#c7ccd4";
  const ENR_C = "#dc2626", DEP_C = "#2563eb", REF_C = "#94a3b8", NULL_C = "rgba(148,163,184,0.28)";

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const geneSelect = $("#gene-select"), radiusRange = $("#radius-range"), radiusVal = $("#radius-val");
  const sphereShow = $("#sphere-show"), spermShow = $("#sperm-show"), dotsShow = $("#dots-show");
  const readoutEl = $("#ss-readout");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), drawerBody = $("#drawer-body");
  const drawerGene = $("#drawer-gene");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle"), bestList = $("#best-list");
  const rankDirEl = $("#rank-dir");
  const xsTabs = $("#xs-tabs"), xsPanels = $("#xs-panels");
  const segBox = $("#seg-controls"), segBoxR = $("#seg-controls-r"), minCountEl = $("#min-count"), minCountREl = $("#min-count-r");

  const state = {
    data: null, byId: {}, radii: [], meta: null, segs: [], segMeta: [], segsOn: {},
    currentId: null, scene: null, gene: null, ri: 2, minCount: 10,
    drawerOpen: false, tab: "curve", rankDir: "enr", dotSize: 2, vcExtras: null, _sceneCache: {},
  };
  const cur = () => state.byId[state.currentId];
  const geneRec = (g) => { const e = cur(); return e && e.genes[g]; };

  function plotInto(div, traces, layout, cfg) {
    if (!div.classList.contains("js-plotly-plot")) div.innerHTML = "";
    Plotly.react(div, traces, layout, cfg || { responsive: true, displayModeBar: false });
  }
  const shown = (el) => !!(el && el.offsetParent);
  const fmtP = (p) => (p == null || !isFinite(p) ? "–" : p < 1e-4 ? p.toExponential(1) : p < 0.1 ? p.toPrecision(2) : p.toFixed(2));
  const sig = (p) => (p != null && p <= 0.05);

  // ── binomial null (regularized incomplete beta; matches scipy) ──
  function gammaln(x) {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015; for (let j = 0; j < 6; j++) { y++; ser += c[j] / y; }
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }
  function betacf(a, b, x) {
    const MAXIT = 200, EPS = 3e-12, FPMIN = 1e-300;
    const qab = a + b, qap = a + 1, qam = a - 1; let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d;
      const del = d * c; h *= del; if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }
  function betai(a, b, x) {
    if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(gammaln(a + b) - gammaln(a) - gammaln(b) + a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
  }
  const binomSF = (k, n, p) => (k <= 0 ? 1 : k > n ? 0 : betai(k, n - k + 1, p));   // P(X >= k)
  const binomCDF = (k, n, p) => (k >= n ? 1 : k < 0 ? 0 : 1 - betai(k + 1, n - k, p)); // P(X <= k)

  (async function init() {
    try {
      const d = await V.loadGz("data/sperm_sphere.json.gz");
      state.data = d; state.radii = d.radii; state.meta = d.meta;
      state.segs = d.segs; state.segMeta = d.segMeta;
      state.minCount = d.meta.minCount || 10;
      state.ri = d.meta.defaultRadiusIdx ?? 2;
      d.embryos.forEach((e) => (state.byId[e.id] = e));
      // segments present in ANY zygote → the toggle set (default all on = "any transcripts")
      const un = new Set(); d.embryos.forEach((e) => (e.present || []).forEach((s) => un.add(s)));
      state.presentUnion = d.segs.filter((s) => un.has(s));
      state.presentUnion.forEach((s) => (state.segsOn[s] = true));
      countEl.textContent = `${d.embryos.length} sperm-positive zygotes · concentration in a sphere at the entry site`;
      radiusRange.max = String(d.radii.length - 1); radiusRange.value = String(state.ri);
      radiusVal.textContent = `${d.radii[state.ri]} µm`;
      if (minCountEl) minCountEl.value = String(state.minCount);
      if (minCountREl) minCountREl.value = String(state.minCount);
      buildSegControls();
      V.buildTabs(tabsEl, d.embryos, selectEmbryo, (e) => ({
        label: e.label, sub: e.date_short,
        title: `${e.label} · ${Object.keys(e.genes).length} genes · sperm ${e.r_com_um}µm from COM`,
      }));
      populateGenes();
      wireControls(); wireDrawer(); wireRdrawer();
    } catch (err) { showError("Failed to load: " + (err.message || err)); }
  })();

  // segment include/exclude checkboxes (shared: applies to sphere AND cell, all charts).
  // rendered into BOTH the bottom-drawer row and the right drawer, kept in sync.
  function renderSegInto(box) {
    if (!box) return;
    const meta = Object.fromEntries(state.segMeta.map((m) => [m.key, m]));
    box.innerHTML = `<span class="seg-lbl">Count transcripts in:</span>` +
      state.presentUnion.map((s) => {
        const m = meta[s] || { label: s, color: "#94a3b8" };
        return `<label class="seg-tg"><input type="checkbox" data-seg="${s}" ${state.segsOn[s] ? "checked" : ""}>` +
          `<span class="seg-sw" style="background:${m.color}"></span>${m.label}</label>`;
      }).join("");
  }
  function syncSegBoxes() {
    [segBox, segBoxR].forEach((box) => { if (!box) return;
      box.querySelectorAll("input[data-seg]").forEach((c) => { c.checked = !!state.segsOn[c.dataset.seg]; }); });
  }
  function onSegChange(ev) {
    const c = ev.target.closest("input[data-seg]"); if (!c) return;
    state.segsOn[c.dataset.seg] = c.checked;
    // never allow zero segments
    if (!state.presentUnion.some((s) => state.segsOn[s])) state.segsOn[c.dataset.seg] = true;
    syncSegBoxes();
    render(); renderReadout(); renderActive(); renderRanks();
  }
  function buildSegControls() {
    renderSegInto(segBox); renderSegInto(segBoxR);
    [segBox, segBoxR].forEach((box) => { if (box && !box._ssBound) { box._ssBound = true; box.addEventListener("change", onSegChange); } });
  }
  const segMask = () => state.segs.map((s) => (state.segsOn[s] ? 1 : 0));
  const sumSeg = (arr, m) => arr.reduce((a, v, i) => a + (m[i] ? v : 0), 0);
  const sumSegR = (a2, m, ri) => a2.reduce((a, row, i) => a + (m[i] ? row[ri] : 0), 0);

  // the concentration statistic for one (zygote, gene, radius) over the selected segments
  function stat(e, g, ri) {
    const rec = e.genes[g]; if (!rec) return null;
    const m = segMask();
    const nC = sumSeg(rec.nc, m), nS = sumSegR(rec.ns, m, ri);
    const vC = sumSeg(e.vc, m), vS = sumSegR(e.vs, m, ri);
    if (!vC || !vS || !nC) return { nC, nS, vC, vS, p0: 0, fold: 0, pEnr: 1, pDep: 1, foldLo: 0, foldHi: 0 };
    const p0 = vS / vC, mu = nC * p0, sd = Math.sqrt(nC * p0 * (1 - p0));
    const klo = Math.max(0, mu - 1.96 * sd), khi = Math.min(nC, mu + 1.96 * sd);
    return {
      nC, nS, vC, vS, p0, fold: nS / mu,                    // fold = (nS/vS)/(nC/vC) = nS/(nC·p0)
      pEnr: binomSF(nS, nC, p0), pDep: binomCDF(nS, nC, p0),
      foldLo: mu > 0 ? klo / mu : 0, foldHi: mu > 0 ? khi / mu : 0,   // 95% band on the fold
    };
  }
  // whole-transcriptome occupancy (all genes) over the selected segments
  function occStat(e, ri) {
    const m = segMask();
    let nS = 0, nC = 0;
    for (const g in e.genes) { const r = e.genes[g]; nS += sumSegR(r.ns, m, ri); nC += sumSeg(r.nc, m); }
    const vC = sumSeg(e.vc, m), vS = sumSegR(e.vs, m, ri), p0 = vC ? vS / vC : 0;
    const se = nC ? Math.sqrt(p0 * (1 - p0) / nC) : 0;
    return { obs: nC ? nS / nC : 0, exp: p0, lo: Math.max(0, p0 - 1.96 * se), hi: Math.min(1, p0 + 1.96 * se) };
  }

  function populateGenes() {
    const cov = {};
    state.data.embryos.forEach((e) => Object.keys(e.genes).forEach((g) => (cov[g] = (cov[g] || 0) + 1)));
    const genes = Object.keys(cov).sort((a, b) => a.localeCompare(b, undefined, { sensitivity: "base", numeric: true }));
    geneSelect.innerHTML = genes.map((g) => `<option value="${g}">${g} (${cov[g]})</option>`).join("");
    state.gene = Object.keys(cov).sort((a, b) => cov[b] - cov[a])[0] || genes[0];
    geneSelect.value = state.gene;
  }
  const gene = () => geneSelect.value;
  const radiusUm = () => state.radii[state.ri];

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id; V.markActiveTab(tabsEl, id);
    const e = state.byId[id]; showLoading(`Loading ${e.label}…`);
    try {
      let sc = state._sceneCache[id];
      if (!sc) { sc = await V.loadGz(`data/zygote/${id}.json.gz`); state._sceneCache[id] = sc; }
      if (state.currentId !== id) return;
      state.scene = sc;
      if (!state.vcExtras)
        state.vcExtras = V.addWindowExtras($("#controls-body"), { defaultSize: state.dotSize, onDotSize: (s) => { state.dotSize = s; render(); } });
      state.vcExtras.setAtlas(id);
      controlsEl.hidden = false; placeholder.hidden = true; drawer.hidden = false; rdrawer.hidden = false;
      render(); renderReadout();
      if (!state.drawerOpen) openDrawer(true); else renderActive();
      renderRanks();
    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }

  // ───────────────────────── 3-D ─────────────────────────
  function sphereMesh(center, rUm, zs, color, op) {
    const nu = 24, nv = 16, rx = rUm / XY, rz = rUm * zs;
    const x = [], y = [], z = [];
    for (let iv = 0; iv <= nv; iv++) { const th = Math.PI * iv / nv;
      for (let iu = 0; iu < nu; iu++) { const ph = 2 * Math.PI * iu / nu;
        x.push(center[0] + rx * Math.sin(th) * Math.cos(ph)); y.push(center[1] + rx * Math.sin(th) * Math.sin(ph)); z.push(center[2] + rz * Math.cos(th)); } }
    const idx = (a, b) => a * nu + (b % nu); const I = [], J = [], K = [];
    for (let iv = 0; iv < nv; iv++) for (let iu = 0; iu < nu; iu++) {
      const a = idx(iv, iu), b = idx(iv, iu + 1), c = idx(iv + 1, iu), dd = idx(iv + 1, iu + 1); I.push(a, b); J.push(b, dd); K.push(c, c); }
    return { type: "mesh3d", x, y, z, i: I, j: J, k: K, color, opacity: op, hoverinfo: "skip",
      flatshading: false, name: "sperm sphere", showlegend: true, legendrank: 30000,
      lighting: { ambient: 0.8, diffuse: 0.4, specular: 0.05 } };
  }
  function render() {
    const s = state.scene, e = cur(); if (!s || !e) return;
    const zs = s.z_scale, g = gene(), r = radiusUm(), sp = e.sperm_plot, spUm = e.sperm_um;
    const traces = V.bodyTraces(s);
    if (dotsShow.checked && s.transcripts[g]) {
      const t = s.transcripts[g], ix = [], iy = [], iz = [], ox = [], oy = [], oz = [];
      for (let k = 0; k < t.x.length; k++) {                 // ANY transcript of the gene (all segments)
        const px = t.x[k], py = t.y[k], pz = t.gz[k];
        const dx = px * XY - spUm[0], dy = py * XY - spUm[1], dz = pz * 1.0 - spUm[2], zplot = pz * zs;
        if (dx * dx + dy * dy + dz * dz <= r * r) { ix.push(px); iy.push(py); iz.push(zplot); }
        else { ox.push(px); oy.push(py); oz.push(zplot); }
      }
      traces.push({ type: "scatter3d", mode: "markers", name: `${g} · outside`, x: ox, y: oy, z: oz,
        marker: { size: state.dotSize, color: OUT_C, opacity: 0.5, line: { width: 0 } }, hovertemplate: `${g} · outside sphere<extra></extra>`, legendrank: 20001 });
      traces.push({ type: "scatter3d", mode: "markers", name: `${g} · within ${r}µm`, x: ix, y: iy, z: iz,
        marker: { size: state.dotSize + 1.5, color: IN_C, opacity: 0.95, line: { width: 0 } }, hovertemplate: `${g} · within ${r}µm of the sperm<extra></extra>`, legendrank: 20000 });
    }
    if (sphereShow.checked) traces.push(sphereMesh(sp, r, zs, SPHERE_C, 0.16));
    if (spermShow.checked) traces.push({ type: "scatter3d", mode: "markers", name: "sperm entry site",
      x: [sp[0]], y: [sp[1]], z: [sp[2]], marker: { size: 8, color: SPERM_C, symbol: "diamond", line: { width: 1.4, color: "#fff" } },
      hovertemplate: "sperm entry site<extra></extra>", legendrank: 30001 });
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  function renderReadout() {
    const e = cur(), g = gene(), r = radiusUm(), ri = state.ri, st = stat(e, g, ri);
    const occ = occStat(e, ri);
    let html = `<div class="ss-r-line"><b>${g}</b> · sphere r = <b>${r} µm</b></div>`;
    if (!st || !st.nC) {
      html += `<div class="ss-r-na">not detected in the selected segments of this zygote</div>`;
    } else {
      const enr = st.fold >= 1, pB = enr ? st.pEnr : st.pDep;
      const cls = st.fold >= 1.5 ? "ss-enr" : st.fold <= 0.67 ? "ss-dep" : "";
      html += `<div class="ss-big ${cls}"><span>${st.fold.toFixed(2)}×</span> <span class="ss-lbl">concentration fold ${enr ? "(enriched)" : "(depleted)"}</span></div>`;
      html += `<div class="ss-r-line">${st.nS} of ${st.nC} in the sphere · expected ${(st.nC * st.p0).toFixed(1)}</div>`;
      html += `<div class="ss-r-line">binomial p ${enr ? "enr" : "dep"} = <b class="${sig(pB) ? "ss-sig" : ""}">${fmtP(pB)}</b></div>`;
    }
    html += `<div class="ss-r-occ">whole cell: <b class="${occ.obs < occ.exp ? "ss-dep" : ""}">${(occ.obs * 100).toFixed(2)}%</b> of transcripts in the sphere ` +
            `vs ${(occ.exp * 100).toFixed(2)}% expected</div>`;
    readoutEl.innerHTML = html;
  }

  // ───────────────────────── bottom drawer ─────────────────────────
  const RENDER = { curve: renderCurve, rank: renderRank, across: renderAcross, occ: renderOcc };
  function renderActive() { (RENDER[state.tab] || renderCurve)(); }
  function switchTab(which) {
    if (!RENDER[which]) which = "curve";
    state.tab = which;
    xsTabs.querySelectorAll(".xs-gtab").forEach((t) => { const on = t.dataset.tab === which; t.classList.toggle("active", on); t.setAttribute("aria-selected", String(on)); });
    xsPanels.querySelectorAll(".xs-panel").forEach((p) => (p.hidden = p.dataset.tab !== which));
    renderActive(); requestAnimationFrame(resizeAll);
  }
  const resizeAll = () => ["curve-plot", "rank-plot", "across-plot", "occ-plot"].forEach((id) => { try { Plotly.Plots.resize($("#" + id)); } catch (_) {} });
  function baseLayout(xt, yt) {
    return { margin: { l: 54, r: 12, t: 8, b: 40 }, autosize: true,
      xaxis: { title: { text: xt, font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", zeroline: false },
      yaxis: { title: { text: yt, font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", zeroline: false },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", legend: { orientation: "h", font: { size: 10 }, y: 1.14, x: 1, xanchor: "right" } };
  }

  // fold vs radius for the current (zygote, gene) with the binomial 95% band (grey only)
  function renderCurve() {
    const div = $("#curve-plot"); if (!shown(div)) return;
    const e = cur(), g = gene(), R = state.radii;
    $("#curve-sub").textContent = `· ${g} · ${e.label}`;
    const sts = R.map((r, i) => stat(e, g, i));
    if (!sts[0] || !sts[0].nC) { Plotly.purge(div); div.classList.remove("js-plotly-plot");
      div.innerHTML = `<div class="ss-empty"><b>${g}</b> not detected in the selected segments of ${e.label}.</div>`; return; }
    const lo = sts.map((s) => s.foldLo), hi = sts.map((s) => s.foldHi), fold = sts.map((s) => s.fold);
    const band = { type: "scatter", mode: "lines", x: R.concat([...R].reverse()), y: hi.concat([...lo].reverse()),
      fill: "toself", fillcolor: NULL_C, line: { width: 0, color: NULL_C }, hoverinfo: "skip", name: "95% null (binomial)" };
    const line = { type: "scatter", mode: "lines+markers", x: R, y: fold, name: "sperm site",
      line: { color: SPERM_C, width: 2.6 }, marker: { size: 7, color: SPERM_C },
      customdata: sts.map((s) => [s.nS, s.nC]), hovertemplate: "r=%{x}µm · fold %{y:.2f} · %{customdata[0]}/%{customdata[1]} in sphere<extra></extra>" };
    const ref = { type: "scatter", mode: "lines", x: [R[0], R[R.length - 1]], y: [1, 1], line: { color: REF_C, width: 1, dash: "dash" }, hoverinfo: "skip", name: "fold = 1 (no change)" };
    const lay = baseLayout("sphere radius (µm)", "concentration fold (sphere ÷ cell)");
    lay.shapes = [{ type: "line", x0: radiusUm(), x1: radiusUm(), y0: 0, y1: 1, yref: "paper", line: { color: "#0f172a", width: 1, dash: "dot" } }];
    plotInto(div, [band, ref, line], lay);
  }

  // genes ranked at the selected radius, filtered by min transcript count, click to view
  function renderRank() {
    const div = $("#rank-plot"); if (!shown(div)) return;
    const e = cur(), ri = state.ri, dir = state.rankDir, mc = state.minCount;
    $("#rank-r").textContent = `${radiusUm()} µm`;
    const rows = Object.keys(e.genes).map((g) => ({ g, s: stat(e, g, ri) })).filter((o) => o.s && o.s.nC >= mc);
    const filt = dir === "enr" ? rows.filter((o) => o.s.fold >= 1).sort((a, b) => b.s.fold - a.s.fold)
      : rows.filter((o) => o.s.fold < 1).sort((a, b) => a.s.fold - b.s.fold);
    const top = filt.slice(0, 20).reverse();
    $("#rank-sub").textContent = `· ${e.label} · ${filt.length} ${dir === "enr" ? "enriched" : "depleted"} (≥ ${mc} transcripts)`;
    if (!top.length) { Plotly.purge(div); div.classList.remove("js-plotly-plot"); div.innerHTML = `<div class="ss-empty">No ${dir === "enr" ? "enriched" : "depleted"} genes with ≥ ${mc} transcripts at this radius.</div>`; return; }
    const colr = top.map((o) => (sig(dir === "enr" ? o.s.pEnr : o.s.pDep) ? (dir === "enr" ? ENR_C : DEP_C) : "#cbd5e1"));
    const tr = { type: "bar", orientation: "h", x: top.map((o) => o.s.fold), y: top.map((o) => o.g),
      marker: { color: colr }, customdata: top.map((o) => [o.s.nS, o.s.nC, dir === "enr" ? o.s.pEnr : o.s.pDep]),
      hovertemplate: "%{y}: fold %{x:.2f} · %{customdata[0]}/%{customdata[1]} · p %{customdata[2]:.2g}  (click to view)<extra></extra>" };
    const lay = baseLayout("concentration fold (sphere ÷ cell)", "");
    lay.margin.l = 96; lay.showlegend = false;
    lay.shapes = [{ type: "line", x0: 1, x1: 1, y0: -0.5, y1: top.length - 0.5, line: { color: REF_C, width: 1, dash: "dash" } }];
    plotInto(div, [tr], lay);
    if (!div._ssBound) { div._ssBound = true; div.on("plotly_click", (ev) => { const p = ev.points && ev.points[0]; if (p && p.y) { geneSelect.value = p.y; onGene(); } }); }
  }

  // this gene across every sperm zygote (legend + 95% null error bars)
  function renderAcross() {
    const div = $("#across-plot"); if (!shown(div)) return;
    const g = gene(), ri = state.ri;
    const rows = state.data.embryos.map((e) => ({ e, s: stat(e, g, ri) })).filter((o) => o.s && o.s.nC);
    $("#across-sub").textContent = `· ${g} · ${rows.length} of ${state.data.embryos.length} zygotes · r=${radiusUm()}µm`;
    if (!rows.length) { Plotly.purge(div); div.classList.remove("js-plotly-plot"); div.innerHTML = `<div class="ss-empty"><b>${g}</b> is not in the selected segments of any sperm zygote.</div>`; return; }
    rows.sort((a, b) => b.s.fold - a.s.fold);
    const grp = { enrSig: [], enrNs: [], depSig: [], depNs: [] };
    rows.forEach((o) => { const f = o.s.fold, s = f >= 1 ? (sig(o.s.pEnr) ? "enrSig" : "enrNs") : (sig(o.s.pDep) ? "depSig" : "depNs"); grp[s].push(o); });
    const cats = [
      ["enrSig", ENR_C, "enriched (p ≤ 0.05)"], ["enrNs", "#f9a8c9", "enriched (n.s.)"],
      ["depSig", DEP_C, "depleted (p ≤ 0.05)"], ["depNs", "#bcd0f0", "depleted (n.s.)"],
    ];
    const traces = cats.filter(([k]) => grp[k].length).map(([k, col, nm]) => ({
      type: "bar", name: nm, x: grp[k].map((o) => o.e.label), y: grp[k].map((o) => o.s.fold),
      marker: { color: col },
      error_y: { type: "data", symmetric: false, array: grp[k].map((o) => o.s.foldHi - o.s.fold), arrayminus: grp[k].map((o) => o.s.fold - o.s.foldLo), color: "rgba(71,85,105,0.55)", thickness: 1, width: 2 },
      customdata: grp[k].map((o) => [o.s.nS, o.s.nC]),
      hovertemplate: "%{x}: fold %{y:.2f} · %{customdata[0]}/%{customdata[1]}<extra>" + nm + "</extra>",
    }));
    const lay = baseLayout("", "concentration fold (sphere ÷ cell)");
    lay.xaxis.tickangle = -48; lay.xaxis.tickfont = { size: 9 }; lay.barmode = "group"; lay.showlegend = true;
    lay.legend = { orientation: "h", font: { size: 9 }, y: 1.13, x: 0, xanchor: "left" };
    lay.shapes = [{ type: "line", x0: -0.5, x1: rows.length - 0.5, y0: 1, y1: 1, line: { color: REF_C, width: 1, dash: "dash" } }];
    plotInto(div, traces, lay);
  }

  // whole-transcriptome occupancy vs radius, with a 95% shaded null band
  function renderOcc() {
    const div = $("#occ-plot"); if (!shown(div)) return;
    const e = cur(), R = state.radii;
    $("#occ-sub").textContent = `· ${e.label}`;
    const os = R.map((r, i) => occStat(e, i));
    const obs = os.map((o) => o.obs * 100), exp = os.map((o) => o.exp * 100);
    const lo = os.map((o) => o.lo * 100), hi = os.map((o) => o.hi * 100);
    const band = { type: "scatter", mode: "lines", x: R.concat([...R].reverse()), y: hi.concat([...lo].reverse()),
      fill: "toself", fillcolor: NULL_C, line: { width: 0, color: NULL_C }, hoverinfo: "skip", name: "95% null (binomial)" };
    const tExp = { type: "scatter", mode: "lines", x: R, y: exp, name: "expected (uniform)", line: { color: REF_C, width: 1.6, dash: "dash" } };
    const tObs = { type: "scatter", mode: "lines+markers", x: R, y: obs, name: "observed", line: { color: SPERM_C, width: 2.6 }, marker: { size: 7 } };
    plotInto(div, [band, tExp, tObs], baseLayout("sphere radius (µm)", "% of transcripts in the sphere"));
  }

  // ───────────────────────── right drawer: cross-embryo consistency ─────────────────────────
  function renderRanks() {
    if (!state.data) return;
    const ri = state.ri, dir = state.rankDir === "dep" ? "dep" : "enr", mc = state.minCount, curG = gene();
    // per gene: # zygotes it is significantly enriched/depleted in (selected segments, this radius)
    const agg = {};
    state.data.embryos.forEach((e) => {
      for (const g in e.genes) {
        const s = stat(e, g, ri); if (!s || s.nC < mc) continue;
        const a = agg[g] || (agg[g] = { nz: 0, enr: 0, dep: 0 });
        a.nz++;
        if (s.fold >= 1.5 && sig(s.pEnr)) a.enr++;
        if (s.fold <= 0.67 && sig(s.pDep)) a.dep++;
      }
    });
    const rows = Object.entries(agg).map(([g, a]) => ({ g, k: a[dir], nz: a.nz })).filter((r) => r.k > 0)
      .sort((a, b) => b.k - a.k || b.nz - a.nz).slice(0, 80);
    const word = dir === "enr" ? "enriched" : "depleted";
    let html = `<div class="best-head best-head-ss"><span></span><span>gene</span>` +
      `<span title="zygotes where this gene is significantly ${word} near the sperm (fold ${dir === "enr" ? "≥ 1.5" : "≤ 0.67"}, p ≤ 0.05)"># ${word}</span>` +
      `<span title="zygotes that contain the gene (≥ ${mc} transcripts)">of total</span></div>`;
    html += rows.map((r, i) =>
      `<div class="best-row${r.g === curG ? " current" : ""}" data-gene="${r.g}" title="${r.g}: significantly ${word} near the sperm in ${r.k} of ${r.nz} zygotes at ${radiusUm()}µm (≥ ${mc} transcripts each)">` +
      `<span class="best-num">${i + 1}</span><span class="best-gene">${r.g}</span>` +
      `<span class="best-real" style="color:${dir === "enr" ? ENR_C : DEP_C}">${r.k}</span><span class="best-p">${r.nz}</span></div>`).join("");
    bestList.innerHTML = html + (rows.length ? "" : `<div class="ss-empty">No gene is significantly ${word} near the sperm at ${radiusUm()}µm (≥ ${mc} transcripts) in any zygote.</div>`);
  }

  // ───────────────────────── wiring ─────────────────────────
  function onGene() { state.gene = gene(); render(); renderReadout(); renderActive(); highlightRank(); }
  function highlightRank() { const g = gene(); bestList.querySelectorAll(".best-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === g)); }
  function onRadius() { state.ri = parseInt(radiusRange.value, 10) || 0; radiusVal.textContent = `${radiusUm()} µm`; render(); renderReadout(); renderActive(); renderRanks(); }
  function onMinCount(v) {
    const n = Math.max(1, parseInt(v, 10) || 1); state.minCount = n;
    if (minCountEl && minCountEl.value !== String(n)) minCountEl.value = String(n);
    if (minCountREl && minCountREl.value !== String(n)) minCountREl.value = String(n);
    renderRank(); renderRanks();
  }
  function wireControls() {
    geneSelect.addEventListener("change", onGene);
    radiusRange.addEventListener("input", onRadius);
    [sphereShow, spermShow, dotsShow].forEach((c) => c.addEventListener("change", () => render()));
    rankDirEl.addEventListener("change", () => { state.rankDir = rankDirEl.value; renderRank(); });
    if (minCountEl) minCountEl.addEventListener("change", () => onMinCount(minCountEl.value));
    if (minCountREl) minCountREl.addEventListener("change", () => onMinCount(minCountREl.value));
  }

  function openDrawer(open) {
    state.drawerOpen = open; drawer.dataset.open = open ? "true" : "false";
    drawerHandle.setAttribute("aria-expanded", String(open)); drawerGene.textContent = `· ${gene()}`;
    if (open) { renderActive(); requestAnimationFrame(resizeAll); }
  }
  function wireDrawer() {
    xsTabs.addEventListener("click", (e) => { const t = e.target.closest(".xs-gtab"); if (t) switchTab(t.dataset.tab); });
    wireHandleDrag(drawer, drawerHandle, {
      computeSize: (e) => window.innerHeight - e.clientY - 40,
      clampSize: (px) => Math.max(200, Math.min(window.innerHeight - 100, px)),
      applySize: (px) => drawer.style.setProperty("--drawer-h", px + "px"), setOpen: openDrawer, afterDrag: resizeAll,
    });
    const rz = $("#drawer-resize"); let sh = 0;
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return; drawer.style.setProperty("--drawer-h", Math.max(200, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y))) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} resizeAll(); } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
    xsPanels.querySelectorAll(".xs-resizable").forEach((box) => { const plot = box.querySelector(".xs-plot"); new ResizeObserver(() => { if (box.offsetParent) { try { Plotly.Plots.resize(plot); } catch (_) {} } }).observe(box); });
    $("#curve-download").addEventListener("click", () => dl($("#curve-plot"), "fold-vs-radius"));
    $("#across-download").addEventListener("click", () => dl($("#across-plot"), "across-zygotes"));
    $("#occ-download").addEventListener("click", () => dl($("#occ-plot"), "occupancy"));
  }
  function dl(div, name) { try { Plotly.downloadImage(div, { format: "png", scale: 4, width: 1400, height: 900, filename: `sperm-sphere-${state.currentId}-${name}` }); } catch (_) {} }

  function wireRdrawer() {
    $("#rtabs").addEventListener("click", (e) => { const t = e.target.closest(".rtab"); if (!t) return;
      state.rankDir = t.dataset.dir; $("#rtabs").querySelectorAll(".rtab").forEach((b) => b.classList.toggle("active", b === t));
      if (rankDirEl) rankDirEl.value = state.rankDir; renderRanks(); if (state.tab === "rank") renderRank(); });
    bestList.addEventListener("click", (e) => { const row = e.target.closest(".best-row"); if (row && row.dataset.gene) { geneSelect.value = row.dataset.gene; onGene(); } });
    wireHandleDrag(rdrawer, rdrawerHandle, {
      computeSize: (e) => window.innerWidth - e.clientX, clampSize: (px) => Math.max(240, Math.min(window.innerWidth - 80, px)),
      applySize: (px) => rdrawer.style.setProperty("--rdrawer-w", px + "px"),
      setOpen: (o) => { rdrawer.dataset.open = o ? "true" : "false"; rdrawerHandle.setAttribute("aria-expanded", String(o)); },
    });
    const rrz = $("#rdrawer-resize"); let sw = 0;
    rrz.addEventListener("pointerdown", (e) => { sw = rdrawer.getBoundingClientRect().width; rrz._d = { x: e.clientX }; rrz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rrz.addEventListener("pointermove", (e) => { if (!rrz._d) return; rdrawer.style.setProperty("--rdrawer-w", Math.max(240, Math.min(window.innerWidth - 80, sw - (e.clientX - rrz._d.x))) + "px"); });
    const end = (e) => { if (rrz._d) { rrz._d = null; try { rrz.releasePointerCapture(e.pointerId); } catch (_) {} } };
    rrz.addEventListener("pointerup", end); rrz.addEventListener("pointercancel", end);
  }
  function wireHandleDrag(el, handle, cfg) {
    let start = null, moved = false;
    handle.addEventListener("pointerdown", (e) => { if (e.button && e.button !== 0) return; start = { x: e.clientX, y: e.clientY }; moved = false; try { handle.setPointerCapture(e.pointerId); } catch (_) {} });
    handle.addEventListener("pointermove", (e) => { if (!start) return;
      if (!moved && Math.hypot(e.clientX - start.x, e.clientY - start.y) < 5) return;
      if (!moved) { moved = true; el.classList.add("dragging"); if (el.dataset.open !== "true") cfg.setOpen(true); }
      cfg.applySize(cfg.clampSize(cfg.computeSize(e))); e.preventDefault(); });
    const up = (e) => { if (!start) return; try { handle.releasePointerCapture(e.pointerId); } catch (_) {}
      if (moved) { el.classList.remove("dragging"); cfg.afterDrag && cfg.afterDrag(); } else cfg.setOpen(el.dataset.open !== "true");
      start = null; moved = false; };
    handle.addEventListener("pointerup", up); handle.addEventListener("pointercancel", up);
  }

  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false; placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
