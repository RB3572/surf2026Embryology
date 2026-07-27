/* Sperm-Entry-Site Enrichment — for zygotes with a labelled sperm, test whether each
 * gene's transcripts are enriched or depleted within a sphere (radius set in the floating
 * window) around the sperm entry site, clipped to the cell. 3-D: the cell, the sperm dot,
 * a sphere at the site, and the selected gene's molecules coloured in/out of the sphere.
 * Drawer: fold vs radius (with a cortical-random-site null band), a per-gene ranking, this
 * gene across every sperm zygote, and the whole-transcriptome occupancy. All statistics are
 * precomputed (build_sperm_sphere.py); the 3-D scene reuses data/zygote/<id>.json.gz.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const SPERM_C = "#db2777", SPHERE_C = "#f472b6";
  const IN_C = "#db2777", OUT_C = "#c7ccd4", BODY_LB = "#cbd5e1";
  const ENR_C = "#dc2626", DEP_C = "#2563eb", REF_C = "#94a3b8";

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

  const state = {
    data: null, byId: {}, radii: [], meta: null,
    currentId: null, scene: null, gene: null, ri: 2,
    drawerOpen: false, tab: "curve", rankDir: "enr", dotSize: 2, vcExtras: null,
    _sceneCache: {},
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

  (async function init() {
    try {
      const d = await V.loadGz("data/sperm_sphere.json.gz");
      state.data = d; state.radii = d.radii; state.meta = d.meta;
      state.ri = d.meta.defaultRadiusIdx ?? 2;
      d.embryos.forEach((e) => (state.byId[e.id] = e));
      countEl.textContent = `${d.embryos.length} sperm-positive zygotes · sphere around the entry site`;
      radiusRange.max = String(d.radii.length - 1);
      radiusRange.value = String(state.ri);
      radiusVal.textContent = `${d.radii[state.ri]} µm`;
      V.buildTabs(tabsEl, d.embryos, selectEmbryo, (e) => ({
        label: e.label, sub: e.date_short,
        title: `${e.label} · ${Object.keys(e.genes).length} genes · sperm ${e.r_com_um}µm from COM`,
      }));
      populateGenes();
      wireControls(); wireDrawer(); wireRdrawer();
    } catch (err) { showError("Failed to load: " + (err.message || err)); }
  })();

  // union of genes across all zygotes, sorted; default = the most widely-covered gene
  function populateGenes() {
    const cov = {};
    state.data.embryos.forEach((e) => Object.keys(e.genes).forEach((g) => (cov[g] = (cov[g] || 0) + 1)));
    const genes = Object.keys(cov).sort((a, b) => a.localeCompare(b, undefined, { sensitivity: "base", numeric: true }));
    geneSelect.innerHTML = genes.map((g) => `<option value="${g}">${g} (${cov[g]})</option>`).join("");
    state.gene = Object.keys(cov).sort((a, b) => cov[b] - cov[a])[0] || genes[0];
    geneSelect.value = state.gene;
  }
  // if the selected gene isn't in the current zygote, keep it (readout shows "not detected")
  const gene = () => geneSelect.value;
  const radiusUm = () => state.radii[state.ri];

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    V.markActiveTab(tabsEl, id);
    const e = state.byId[id];
    showLoading(`Loading ${e.label}…`);
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
    const nu = 24, nv = 16, rx = rUm / XY, rz = rUm * zs;   // µm→plot: xy ÷0.15, z ×z_scale
    const x = [], y = [], z = [];
    for (let iv = 0; iv <= nv; iv++) {
      const th = Math.PI * iv / nv;
      for (let iu = 0; iu < nu; iu++) {
        const ph = 2 * Math.PI * iu / nu;
        x.push(center[0] + rx * Math.sin(th) * Math.cos(ph));
        y.push(center[1] + rx * Math.sin(th) * Math.sin(ph));
        z.push(center[2] + rz * Math.cos(th));
      }
    }
    const idx = (a, b) => a * nu + (b % nu);
    const I = [], J = [], K = [];
    for (let iv = 0; iv < nv; iv++) for (let iu = 0; iu < nu; iu++) {
      const a = idx(iv, iu), b = idx(iv, iu + 1), c = idx(iv + 1, iu), dd = idx(iv + 1, iu + 1);
      I.push(a, b); J.push(b, dd); K.push(c, c);
    }
    return { type: "mesh3d", x, y, z, i: I, j: J, k: K, color, opacity: op, hoverinfo: "skip",
      flatshading: false, name: "sperm sphere", showlegend: true, legendrank: 30000,
      lighting: { ambient: 0.8, diffuse: 0.4, specular: 0.05 } };
  }

  function render() {
    const s = state.scene, e = cur(); if (!s || !e) return;
    const zs = s.z_scale, g = gene(), r = radiusUm();
    const sp = e.sperm_plot, spUm = e.sperm_um;
    const traces = V.bodyTraces(s);
    // gene transcripts, split by in/out of the sperm sphere (µm distance)
    if (dotsShow.checked && s.transcripts[g]) {
      const t = s.transcripts[g];
      const ix = [], iy = [], iz = [], ox = [], oy = [], oz = [];
      for (let k = 0; k < t.x.length; k++) {
        if (!(t.s1[k])) continue;                        // in-cell (segment 1) only for clarity
        const px = t.x[k], py = t.y[k], pz = t.gz[k];
        const dx = px * XY - spUm[0], dy = py * XY - spUm[1], dz = pz * 1.0 - spUm[2];
        const zplot = pz * zs;
        if (dx * dx + dy * dy + dz * dz <= r * r) { ix.push(px); iy.push(py); iz.push(zplot); }
        else { ox.push(px); oy.push(py); oz.push(zplot); }
      }
      traces.push({ type: "scatter3d", mode: "markers", name: `${g} · outside`, x: ox, y: oy, z: oz,
        marker: { size: state.dotSize, color: OUT_C, opacity: 0.5, line: { width: 0 } },
        hovertemplate: `${g} · outside sphere<extra></extra>`, legendrank: 20001 });
      traces.push({ type: "scatter3d", mode: "markers", name: `${g} · within ${r}µm`, x: ix, y: iy, z: iz,
        marker: { size: state.dotSize + 1.5, color: IN_C, opacity: 0.95, line: { width: 0 } },
        hovertemplate: `${g} · within ${r}µm of the sperm<extra></extra>`, legendrank: 20000 });
    }
    if (sphereShow.checked) traces.push(sphereMesh(sp, r, zs, SPHERE_C, 0.16));
    if (spermShow.checked) traces.push({ type: "scatter3d", mode: "markers", name: "sperm entry site",
      x: [sp[0]], y: [sp[1]], z: [sp[2]],
      marker: { size: 8, color: SPERM_C, symbol: "diamond", line: { width: 1.4, color: "#fff" } },
      hovertemplate: "sperm entry site<extra></extra>", legendrank: 30001 });
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  function renderReadout() {
    const e = cur(), g = gene(), r = radiusUm(), ri = state.ri;
    const rec = geneRec(g);
    const occFrac = occupancy(e, ri);
    let html = `<div class="ss-r-line"><b>${g}</b> · sphere r = <b>${r} µm</b></div>`;
    if (!rec) {
      html += `<div class="ss-r-na">not detected (or &lt; ${state.meta.minCount} transcripts) in this zygote</div>`;
    } else {
      const fold = rec.fold[ri], nsph = rec.nsph[ri], enr = fold >= 1;
      const pB = enr ? rec.pE[ri] : rec.pD[ri], pS = enr ? rec.pSE[ri] : rec.pSD[ri];
      const cls = fold >= state.meta.foldThresh ? "ss-enr" : fold <= state.meta.depThresh ? "ss-dep" : "";
      html += `<div class="ss-big ${cls}"><span>${fold.toFixed(2)}×</span> <span class="ss-lbl">density fold ${enr ? "(enriched)" : "(depleted)"}</span></div>`;
      html += `<div class="ss-r-line">${nsph} of ${rec.n} in the sphere · expected ${(rec.n * e.p_null[String(r)]).toFixed(1)}</div>`;
      html += `<div class="ss-r-line">binomial p ${enr ? "enr" : "dep"} = <b class="${sig(pB) ? "ss-sig" : ""}">${fmtP(pB)}</b> · ` +
              `cortical-site p = <b class="${sig(pS) ? "ss-sig" : ""}">${fmtP(pS)}</b></div>`;
    }
    html += `<div class="ss-r-occ">whole cell: <b class="${occFrac.dep ? "ss-dep" : ""}">${(occFrac.obs * 100).toFixed(2)}%</b> of transcripts in the sphere ` +
            `vs ${(occFrac.exp * 100).toFixed(2)}% expected</div>`;
    readoutEl.innerHTML = html;
  }
  function occupancy(e, ri) {
    let sph = 0, cell = 0;
    for (const g in e.genes) { sph += e.genes[g].nsph[ri]; cell += e.genes[g].n; }
    const exp = e.p_null[String(state.radii[ri])];
    return { obs: cell ? sph / cell : 0, exp, dep: cell && (sph / cell) < exp };
  }

  // ───────────────────────── bottom drawer ─────────────────────────
  const RENDER = { curve: renderCurve, rank: renderRank, across: renderAcross, occ: renderOcc };
  function renderActive() { (RENDER[state.tab] || renderCurve)(); }
  function switchTab(which) {
    if (!RENDER[which]) which = "curve";
    state.tab = which;
    xsTabs.querySelectorAll(".xs-gtab").forEach((t) => { const on = t.dataset.tab === which; t.classList.toggle("active", on); t.setAttribute("aria-selected", String(on)); });
    xsPanels.querySelectorAll(".xs-panel").forEach((p) => (p.hidden = p.dataset.tab !== which));
    renderActive();
    requestAnimationFrame(() => resizeAll());
  }
  const resizeAll = () => ["curve-plot", "rank-plot", "across-plot", "occ-plot"].forEach((id) => { try { Plotly.Plots.resize($("#" + id)); } catch (_) {} });

  function baseLayout(xt, yt) {
    return { margin: { l: 52, r: 12, t: 8, b: 40 }, autosize: true,
      xaxis: { title: { text: xt, font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", zeroline: false },
      yaxis: { title: { text: yt, font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", zeroline: false },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      legend: { orientation: "h", font: { size: 10 }, y: 1.14, x: 1, xanchor: "right" } };
  }

  // fold vs radius for the current (zygote, gene) with the cortical-null band
  function renderCurve() {
    const div = $("#curve-plot"); if (!shown(div)) return;
    const e = cur(), g = gene(), rec = geneRec(g), R = state.radii;
    $("#curve-sub").textContent = `· ${g} · ${e.label}`;
    if (!rec) { Plotly.purge(div); div.classList.remove("js-plotly-plot");
      div.innerHTML = `<div class="ss-empty"><b>${g}</b> not detected (or &lt; ${state.meta.minCount} transcripts) in ${e.label}.</div>`; return; }
    const band = { type: "scatter", x: R.concat([...R].reverse()), y: rec.nhi.concat([...rec.nlo].reverse()),
      fill: "toself", fillcolor: "rgba(148,163,184,0.22)", line: { width: 0 }, hoverinfo: "skip", name: "cortical null 95%" };
    const line = { type: "scatter", mode: "lines+markers", x: R, y: rec.fold, name: "sperm site",
      line: { color: SPERM_C, width: 2.6 }, marker: { size: 7, color: rec.fold.map((f, i) => (f >= 1 ? (sig(rec.pSE[i]) ? ENR_C : SPERM_C) : (sig(rec.pSD[i]) ? DEP_C : SPERM_C))) },
      customdata: R.map((r, i) => [rec.nsph[i], rec.n, rec.pSE[i], rec.pSD[i]]),
      hovertemplate: "r=%{x}µm · fold %{y:.2f} · %{customdata[0]}/%{customdata[1]} · cortical pE %{customdata[2]:.2f} pD %{customdata[3]:.2f}<extra></extra>" };
    const ref = { type: "scatter", mode: "lines", x: [R[0], R[R.length - 1]], y: [1, 1], line: { color: REF_C, width: 1, dash: "dash" }, hoverinfo: "skip", name: "fold = 1" };
    const lay = baseLayout("sphere radius (µm)", "density fold (sphere ÷ cell)");
    // mark the currently-selected radius
    lay.shapes = [{ type: "line", x0: radiusUm(), x1: radiusUm(), y0: 0, y1: 1, yref: "paper", line: { color: "#0f172a", width: 1, dash: "dot" } }];
    plotInto(div, [band, ref, line], lay);
  }

  // genes ranked at the selected radius (enriched or depleted) for the current zygote
  function renderRank() {
    const div = $("#rank-plot"); if (!shown(div)) return;
    const e = cur(), ri = state.ri, dir = state.rankDir;
    $("#rank-r").textContent = `${radiusUm()} µm`;
    const rows = Object.entries(e.genes).map(([g, rec]) => ({ g, fold: rec.fold[ri], n: rec.n, nsph: rec.nsph[ri], pS: dir === "enr" ? rec.pSE[ri] : rec.pSD[ri] }));
    const filt = dir === "enr"
      ? rows.filter((r) => r.fold >= 1).sort((a, b) => b.fold - a.fold)
      : rows.filter((r) => r.fold < 1).sort((a, b) => a.fold - b.fold);
    const top = filt.slice(0, 20).reverse();
    $("#rank-sub").textContent = `· ${e.label} · ${filt.length} ${dir === "enr" ? "enriched" : "depleted"}`;
    if (!top.length) { Plotly.purge(div); div.classList.remove("js-plotly-plot"); div.innerHTML = `<div class="ss-empty">No ${dir === "enr" ? "enriched" : "depleted"} genes at this radius.</div>`; return; }
    const colr = top.map((r) => (sig(r.pS) ? (dir === "enr" ? ENR_C : DEP_C) : "#cbd5e1"));
    const tr = { type: "bar", orientation: "h", x: top.map((r) => r.fold), y: top.map((r) => r.g),
      marker: { color: colr }, customdata: top.map((r) => [r.nsph, r.n, r.pS]),
      hovertemplate: "%{y}: fold %{x:.2f} · %{customdata[0]}/%{customdata[1]} · cortical p %{customdata[2]:.2f}<extra></extra>" };
    const lay = baseLayout("density fold (sphere ÷ cell)", "");
    lay.margin.l = 92; lay.shapes = [{ type: "line", x0: 1, x1: 1, y0: -0.5, y1: top.length - 0.5, line: { color: REF_C, width: 1, dash: "dash" } }];
    lay.showlegend = false;
    plotInto(div, [tr], lay);
    div.removeAllListeners && div.removeAllListeners("plotly_click");
    if (!div._ssBound) { div._ssBound = true; div.on("plotly_click", (ev) => { const p = ev.points && ev.points[0]; if (p && p.y) { geneSelect.value = p.y; onGene(); } }); }
  }

  // this gene across every sperm zygote, at the selected radius
  function renderAcross() {
    const div = $("#across-plot"); if (!shown(div)) return;
    const g = gene(), ri = state.ri;
    const rows = state.data.embryos.map((e) => ({ e, rec: e.genes[g] })).filter((o) => o.rec);
    $("#across-sub").textContent = `· ${g} · ${rows.length} of ${state.data.embryos.length} zygotes · r=${radiusUm()}µm`;
    if (!rows.length) { Plotly.purge(div); div.classList.remove("js-plotly-plot"); div.innerHTML = `<div class="ss-empty"><b>${g}</b> is not in ≥ ${state.meta.minCount} transcripts in any sperm zygote.</div>`; return; }
    rows.sort((a, b) => b.rec.fold[ri] - a.rec.fold[ri]);
    const colr = rows.map((o) => { const f = o.rec.fold[ri]; return f >= 1 ? (sig(o.rec.pSE[ri]) ? ENR_C : "#f9a8c9") : (sig(o.rec.pSD[ri]) ? DEP_C : "#bcd0f0"); });
    const tr = { type: "bar", x: rows.map((o) => o.e.label), y: rows.map((o) => o.rec.fold[ri]),
      marker: { color: colr }, customdata: rows.map((o) => [o.rec.nsph[ri], o.rec.n]),
      hovertemplate: "%{x}: fold %{y:.2f} · %{customdata[0]}/%{customdata[1]}<extra></extra>" };
    const lay = baseLayout("", "density fold (sphere ÷ cell)");
    lay.xaxis.tickangle = -48; lay.xaxis.tickfont = { size: 9 };
    lay.shapes = [{ type: "line", x0: -0.5, x1: rows.length - 0.5, y0: 1, y1: 1, line: { color: REF_C, width: 1, dash: "dash" } }];
    lay.showlegend = false;
    plotInto(div, [tr], lay);
  }

  // whole-transcriptome occupancy vs radius (observed vs expected)
  function renderOcc() {
    const div = $("#occ-plot"); if (!shown(div)) return;
    const e = cur(), R = state.radii;
    $("#occ-sub").textContent = `· ${e.label}`;
    const obs = R.map((r, i) => occupancy(e, i).obs * 100);
    const exp = R.map((r) => e.p_null[String(r)] * 100);
    const tObs = { type: "scatter", mode: "lines+markers", x: R, y: obs, name: "observed", line: { color: SPERM_C, width: 2.6 }, marker: { size: 7 } };
    const tExp = { type: "scatter", mode: "lines", x: R, y: exp, name: "expected (uniform)", line: { color: REF_C, width: 1.6, dash: "dash" } };
    plotInto(div, [tExp, tObs], baseLayout("sphere radius (µm)", "% of in-cell transcripts in the sphere"));
  }

  // ───────────────────────── right drawer: cross-embryo consistency ─────────────────────────
  function renderRanks() {
    if (!state.data) return;
    const ri = state.ri, dir = state.rankDir === "dep" ? "dep" : "enr";
    const rows = Object.entries(state.data.byGene)
      .map(([g, bg]) => ({ g, k: bg[dir][ri], nz: bg.nz }))
      .filter((r) => r.k > 0)
      .sort((a, b) => b.k - a.k || b.nz - a.nz).slice(0, 60);
    const curG = gene();
    let html = `<div class="best-head"><span></span><span>gene</span><span>${dir === "enr" ? "enr" : "dep"}</span><span>n</span></div>`;
    html += rows.map((r, i) =>
      `<div class="best-row${r.g === curG ? " current" : ""}" data-gene="${r.g}" title="${r.g}: significantly ${dir === "enr" ? "enriched" : "depleted"} near the sperm in ${r.k} of ${r.nz} zygotes at ${radiusUm()}µm">` +
      `<span class="best-num">${i + 1}</span><span class="best-gene">${r.g}</span>` +
      `<span class="best-real" style="color:${dir === "enr" ? ENR_C : DEP_C}">${r.k}</span><span class="best-p">${r.nz}</span></div>`).join("");
    bestList.innerHTML = html || `<div class="ss-empty">No gene is significantly ${dir === "enr" ? "enriched" : "depleted"} near the sperm at ${radiusUm()}µm in any zygote.</div>`;
  }

  // ───────────────────────── wiring ─────────────────────────
  function onGene() { state.gene = gene(); render(); renderReadout(); renderActive(); highlightRank(); }
  function highlightRank() { const g = gene(); bestList.querySelectorAll(".best-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === g)); }
  function onRadius() {
    state.ri = parseInt(radiusRange.value, 10) || 0;
    radiusVal.textContent = `${radiusUm()} µm`;
    render(); renderReadout(); renderActive(); renderRanks();
  }
  function wireControls() {
    geneSelect.addEventListener("change", onGene);
    radiusRange.addEventListener("input", onRadius);
    [sphereShow, spermShow, dotsShow].forEach((c) => c.addEventListener("change", () => render()));
    rankDirEl.addEventListener("change", () => { state.rankDir = rankDirEl.value; renderRank(); });
  }

  function openDrawer(open) {
    state.drawerOpen = open; drawer.dataset.open = open ? "true" : "false";
    drawerHandle.setAttribute("aria-expanded", String(open));
    drawerGene.textContent = `· ${gene()}`;
    if (open) { renderActive(); requestAnimationFrame(resizeAll); }
  }
  function wireDrawer() {
    xsTabs.addEventListener("click", (e) => { const t = e.target.closest(".xs-gtab"); if (t) switchTab(t.dataset.tab); });
    wireHandleDrag(drawer, drawerHandle, {
      computeSize: (e) => window.innerHeight - e.clientY - 40,
      clampSize: (px) => Math.max(200, Math.min(window.innerHeight - 100, px)),
      applySize: (px) => drawer.style.setProperty("--drawer-h", px + "px"),
      setOpen: openDrawer, afterDrag: resizeAll,
    });
    const rz = $("#drawer-resize"); let sh = 0;
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return; drawer.style.setProperty("--drawer-h", Math.max(200, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y))) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} resizeAll(); } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
    // per-panel resize observers
    xsPanels.querySelectorAll(".xs-resizable").forEach((box) => {
      const plot = box.querySelector(".xs-plot");
      new ResizeObserver(() => { if (box.offsetParent) { try { Plotly.Plots.resize(plot); } catch (_) {} } }).observe(box);
    });
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
      computeSize: (e) => window.innerWidth - e.clientX,
      clampSize: (px) => Math.max(240, Math.min(window.innerWidth - 80, px)),
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
