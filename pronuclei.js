/* Pronuclei Distance vs Transcripts — zygotes with two auto-detected pronuclei.
 * Built on viewer-core.js (VCore). The 3-D view shows the two pronuclei and the
 * shortest line between them. The bottom drawer scatters pronuclei distance against
 * (top) the SELECTED GENE's count and (below) the total transcript count, each with
 * a least-squares fit. The right drawer ranks genes by the Pearson correlation of
 * their count with the distance (most positive / most negative). Precompute in
 * build_pronuclei.py; the UI only reads it.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const LINE_C = "#111827", CUR_C = "#0891b2";
  const MIN_ZYG = 5;                    // min zygotes containing a gene to rank/correlate it

  function plotInto(div, traces, layout, cfg) {
    if (!div.classList.contains("js-plotly-plot")) div.innerHTML = "";
    Plotly.react(div, traces, layout, cfg || { responsive: true, displayModeBar: false });
  }

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const pnReadout = $("#pn-readout"), pnFit = $("#pn-fit");
  const geneSelect = $("#gene-select"), geneScatter = $("#gene-scatter"), geneFit = $("#gene-fit");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), drawerBody = $("#drawer-body");
  const scatterPlot = $("#scatter-plot");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle");
  const rankNEl = $("#rank-n"), rankPosEl = $("#rank-pos"), rankNegEl = $("#rank-neg");

  const state = { points: [], byId: {}, genesAgg: null, geneCorr: [], userGene: null, rankN: 10,
                  currentId: null, scene: null, fit: null, drawerOpen: false };

  (async function init() {
    try {
      const [m, ga] = await Promise.all([
        (await fetch("data/pronuclei_manifest.json")).json(),
        V.loadGz("data/pronuclei_genes.json.gz"),
      ]);
      state.points = m.embryos; state.genesAgg = ga;
      state.points.forEach((p) => (state.byId[p.id] = p));
      countEl.textContent = `${m.embryos.length} zygotes · pronuclei auto-detected inside the cytoplasm`;
      state.fit = linreg(state.points.map((p) => p.total), state.points.map((p) => p.distance));
      computeGeneCorr();
      populateGenes();
      V.buildTabs(tabsEl, m.embryos, selectEmbryo, (e) => ({
        label: e.label, sub: e.date_short,
        title: `${e.label} · dist ${e.distance} µm · ${e.total.toLocaleString()} transcripts`,
      }));
      wireDrawer(); wireRdrawer(); renderRanks();
      geneSelect.addEventListener("change", () => selectGene(geneSelect.value));
      rankNEl.addEventListener("change", () => { state.rankN = parseInt(rankNEl.value, 10) || 10; renderRanks(); });
    } catch (err) { showError("Failed to load: " + (err.message || err)); }
  })();

  // ---------- gene ↔ distance correlations ----------
  function geneSeries(g) {
    const xs = [], ys = [], ids = [], labels = [];
    for (const e of state.genesAgg.embryos) {
      const c = e.genes[g]; if (c == null) continue;
      xs.push(c); ys.push(e.distance); ids.push(e.id); labels.push((state.byId[e.id] || {}).label || e.id);
    }
    return { xs, ys, ids, labels };
  }
  function computeGeneCorr() {
    const genes = new Set();
    state.genesAgg.embryos.forEach((e) => Object.keys(e.genes).forEach((g) => genes.add(g)));
    const out = [];
    for (const g of genes) {
      const s = geneSeries(g); if (s.xs.length < MIN_ZYG) continue;
      const f = linreg(s.xs, s.ys); if (!isFinite(f.r)) continue;
      out.push({ gene: g, r: f.r, n: s.xs.length });
    }
    state.geneCorr = out;
  }
  function populateGenes() {
    const genes = [...new Set(state.genesAgg.embryos.flatMap((e) => Object.keys(e.genes)))]
      .sort((a, b) => a.localeCompare(b, undefined, { sensitivity: "base", numeric: true }));
    geneSelect.innerHTML = genes.map((g) => `<option value="${g}">${g}</option>`).join("");
    const top = [...state.geneCorr].sort((a, b) => b.r - a.r)[0];   // default = strongest positive correlation
    state.userGene = (top ? top.gene : genes[0]) || "";
    geneSelect.value = state.userGene;
  }
  const gene = () => geneSelect.value;
  function selectGene(g) {
    state.userGene = g; geneSelect.value = g;
    renderGeneScatter(); highlightRank();
    if (state.scene) renderReadout(state.byId[state.currentId] || {});
  }

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    V.markActiveTab(tabsEl, id);
    const meta = state.byId[id] || {};
    showLoading(`Loading ${meta.label || id}…`);
    try {
      const scene = await V.loadGz(`data/pronuclei/${id}.json.gz`);
      if (state.currentId !== id) return;
      state.scene = scene;
      controlsEl.hidden = false; placeholder.hidden = true; drawer.hidden = false; rdrawer.hidden = false;
      render(); renderReadout(meta);
      if (!state.drawerOpen) openDrawer(true);
      else { renderScatter(); renderGeneScatter(); }
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
    const [la, lb] = s.pron_labels;
    const pcolor = { [la]: "#2563eb", [lb]: "#dc2626" };
    const traces = [];
    for (const lbl of s.mask_labels) {
      const pron = lbl === la || lbl === lb;
      const t = segMesh(s, lbl, pron ? pcolor[lbl] : "#9aa3b2", pron ? 0.5 : 0.08,
        pron ? `Pronucleus (seg ${lbl})` : `Segment ${lbl}`);
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
    const gc = (state.genesAgg.embryos.find((e) => e.id === s.id) || { genes: {} }).genes[gene()];
    pnReadout.innerHTML =
      `<div class="pn-big"><span>${s.distance_um}</span> µm <span class="pn-lbl">pronuclei distance</span></div>` +
      `<div class="pn-big"><span>${s.total_transcripts.toLocaleString()}</span> <span class="pn-lbl">total transcripts</span></div>` +
      `<div class="pn-resid"><b>${gene()}</b> here: ${gc != null ? gc.toLocaleString() + " transcripts" : "not in this zygote's panel"}</div>` +
      `<div class="pn-resid">pronuclei auto-detected as segments <b>${s.pron_labels[0]}</b> &amp; <b>${s.pron_labels[1]}</b></div>`;
  }

  // ---------- scatters ----------
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
  function scatter(div, xs, ys, ids, labels, fit, xTitle, unit, curId) {
    const others = xs.map((x, i) => ({ x, y: ys[i], lab: labels[i], id: ids[i] })).filter((o) => o.id !== curId);
    const cur = xs.map((x, i) => ({ x, y: ys[i], lab: labels[i], id: ids[i] })).find((o) => o.id === curId);
    const xmin = Math.min(...xs), xmax = Math.max(...xs);
    const traces = [
      { type: "scatter", mode: "markers", name: "zygotes",
        x: others.map((o) => o.x), y: others.map((o) => o.y),
        marker: { size: 8, color: "#94a3b8", opacity: 0.75, line: { width: 0 } },
        text: others.map((o) => o.lab), hovertemplate: `%{text}<br>%{x:,} ${unit}<br>%{y} µm<extra></extra>` },
      { type: "scatter", mode: "lines", name: "least-squares fit",
        x: [xmin, xmax], y: [fit.a + fit.b * xmin, fit.a + fit.b * xmax],
        line: { color: "#0891b2", width: 2 }, hoverinfo: "skip" },
    ];
    if (cur) traces.push({ type: "scatter", mode: "markers", name: cur.lab,
      x: [cur.x], y: [cur.y], marker: { size: 15, color: CUR_C, line: { width: 2, color: "#fff" } },
      hovertemplate: `${cur.lab}<br>%{x:,} ${unit}<br>%{y} µm<extra></extra>` });
    plotInto(div, traces, {
      margin: { l: 52, r: 12, t: 6, b: 40 }, height: div.clientHeight || 220,
      xaxis: { title: { text: xTitle, font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", zeroline: false },
      yaxis: { title: { text: "min pronuclei distance (µm)", font: { size: 10 } }, tickfont: { size: 9 }, gridcolor: "#eef1f5", rangemode: "tozero" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      legend: { orientation: "h", font: { size: 10 }, y: 1.16, x: 1, xanchor: "right" },
    });
  }
  function renderScatter() {
    const pts = state.points;
    scatter(scatterPlot, pts.map((p) => p.total), pts.map((p) => p.distance), pts.map((p) => p.id),
      pts.map((p) => p.label), state.fit, "total transcripts", "transcripts", state.currentId);
    const f = state.fit;
    pnFit.innerHTML = `· <b>${f.n}</b> zygotes · r = <b>${f.r.toFixed(3)}</b> (r²=${f.r2.toFixed(3)}) · slope ${(f.b * 1e5).toFixed(1)} µm / 100k`;
  }
  function renderGeneScatter() {
    const g = gene(), s = geneSeries(g);
    geneFit.innerHTML = `· <b>${g}</b>`;
    if (s.xs.length < 2) {
      Plotly.purge(geneScatter); geneScatter.classList.remove("js-plotly-plot");
      geneScatter.innerHTML = `<div class="pn-empty">Only ${s.xs.length} zygote${s.xs.length === 1 ? "" : "s"} contain${s.xs.length === 1 ? "s" : ""} <b>${g}</b> — too few to correlate.</div>`;
      return;
    }
    const f = linreg(s.xs, s.ys);
    scatter(geneScatter, s.xs, s.ys, s.ids, s.labels, f, `${g} transcript count`, g, state.currentId);
    geneFit.innerHTML = `· <b>${g}</b> · r = <b style="color:${f.r >= 0 ? "#dc2626" : "#2563eb"}">${f.r >= 0 ? "+" : ""}${f.r.toFixed(3)}</b> (n=${s.xs.length} zygotes)`;
  }

  // ---------- right drawer: correlation ranking ----------
  function rankRows(rows) {
    const cur = gene();
    let html = `<div class="best-head"><span></span><span>gene</span><span>r</span><span>n</span></div>`;
    html += rows.map((r, i) =>
      `<div class="best-row pn-row${r.gene === cur ? " current" : ""}" data-gene="${r.gene}" ` +
      `title="Pearson r = ${r.r.toFixed(3)} of ${r.gene} count vs pronuclei distance across ${r.n} zygotes">` +
      `<span class="best-num">${i + 1}</span><span class="best-gene">${r.gene}</span>` +
      `<span class="best-real" style="color:${r.r >= 0 ? "#dc2626" : "#2563eb"}">${r.r >= 0 ? "+" : ""}${r.r.toFixed(2)}</span>` +
      `<span class="best-p">${r.n}</span></div>`).join("");
    return html || `<div class="pn-empty">No genes in ≥ ${MIN_ZYG} zygotes.</div>`;
  }
  function renderRanks() {
    const n = state.rankN;
    rankPosEl.innerHTML = rankRows([...state.geneCorr].sort((a, b) => b.r - a.r).slice(0, n));
    rankNegEl.innerHTML = rankRows([...state.geneCorr].sort((a, b) => a.r - b.r).slice(0, n));
  }
  function highlightRank() {
    const cur = gene();
    rdrawer.querySelectorAll(".pn-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === cur));
  }

  // ---------- drawers ----------
  const resizeScatters = () => { try { Plotly.Plots.resize(scatterPlot); Plotly.Plots.resize(geneScatter); } catch (_) {} };
  function openDrawer(open) {
    state.drawerOpen = open;
    drawer.dataset.open = open ? "true" : "false";
    drawerHandle.setAttribute("aria-expanded", String(open));
    if (open) { renderScatter(); renderGeneScatter(); requestAnimationFrame(resizeScatters); }
  }
  function wireDrawer() {
    wireHandleDrag(drawer, drawerHandle, {
      computeSize: (e) => window.innerHeight - e.clientY - 40,
      clampSize: (px) => Math.max(240, Math.min(window.innerHeight - 100, px)),
      applySize: (px) => drawer.style.setProperty("--drawer-h", px + "px"),
      setOpen: openDrawer, afterDrag: resizeScatters,
    });
    let sh = 0; const rz = $("#drawer-resize");
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return;
      drawer.style.setProperty("--drawer-h", Math.max(240, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y))) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} resizeScatters(); } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }
  function wireRdrawer() {
    wireHandleDrag(rdrawer, rdrawerHandle, {
      computeSize: (e) => window.innerWidth - e.clientX,
      clampSize: (px) => Math.max(260, Math.min(window.innerWidth - 80, px)),
      applySize: (px) => rdrawer.style.setProperty("--rdrawer-w", px + "px"),
      setOpen: (o) => { rdrawer.dataset.open = o ? "true" : "false"; rdrawerHandle.setAttribute("aria-expanded", String(o)); },
    });
    let sw = 0; const rrz = $("#rdrawer-resize");
    rrz.addEventListener("pointerdown", (e) => { sw = rdrawer.getBoundingClientRect().width; rrz._d = { x: e.clientX }; rrz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rrz.addEventListener("pointermove", (e) => { if (!rrz._d) return;
      rdrawer.style.setProperty("--rdrawer-w", Math.max(260, Math.min(window.innerWidth - 80, sw - (e.clientX - rrz._d.x))) + "px"); });
    const end = (e) => { if (rrz._d) { rrz._d = null; try { rrz.releasePointerCapture(e.pointerId); } catch (_) {} } };
    rrz.addEventListener("pointerup", end); rrz.addEventListener("pointercancel", end);
    rdrawer.addEventListener("click", (e) => {
      const row = e.target.closest(".pn-row"); if (row) selectGene(row.dataset.gene);
    });
  }
  function wireHandleDrag(el, handle, cfg) {
    let start = null, moved = false;
    handle.addEventListener("pointerdown", (e) => { if (e.button && e.button !== 0) return; start = { x: e.clientX, y: e.clientY }; moved = false; try { handle.setPointerCapture(e.pointerId); } catch (_) {} });
    handle.addEventListener("pointermove", (e) => {
      if (!start) return;
      if (!moved && Math.hypot(e.clientX - start.x, e.clientY - start.y) < 5) return;
      if (!moved) { moved = true; el.classList.add("dragging"); if (el.dataset.open !== "true") cfg.setOpen(true); }
      cfg.applySize(cfg.clampSize(cfg.computeSize(e))); e.preventDefault();
    });
    const up = (e) => { if (!start) return; try { handle.releasePointerCapture(e.pointerId); } catch (_) {}
      if (moved) { el.classList.remove("dragging"); cfg.afterDrag && cfg.afterDrag(); } else cfg.setOpen(el.dataset.open !== "true");
      start = null; moved = false; };
    handle.addEventListener("pointerup", up); handle.addEventListener("pointercancel", up);
  }

  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false;
    placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
