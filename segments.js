/* Segment Gene Enrichment — every embryo (oocyte → late 2-cell).
 * Built on viewer-core.js (VCore). For the selected embryo we render its
 * segmentation meshes and, per segment, the genes ranked by density fold-change
 * (precompute in build_segments.py). Selecting a gene colours the segments by its
 * enrichment and charts it across segments. The UI reads only the precompute.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;

  function plotInto(div, traces, layout, cfg) {
    if (!div.classList.contains("js-plotly-plot")) div.innerHTML = "";
    Plotly.react(div, traces, layout, cfg || { responsive: true, displayModeBar: false });
  }

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const geneSelect = $("#gene-select"), chartEl = $("#chart"), chartSub = $("#chart-sub");
  const chartReadout = $("#chart-readout"), segLegend = $("#seg-legend");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle");
  const rtabsEl = $("#rtabs"), segListEl = $("#seg-list"), rdrawerDesc = $("#rdrawer-desc");
  const minCountEl = $("#min-count"), metricModeEl = $("#metric-mode");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), drawerBody = $("#drawer-body");
  const crossGeneEl = $("#cross-gene"), crossMetricEl = $("#cross-metric"), crossNote = $("#cross-note"), crossChart = $("#cross-chart");
  const geneChips = $("#gene-chips"), segDotToggles = $("#seg-dot-toggles");

  const state = { manifest: [], currentId: null, scene: null, userGene: null,
                  segTab: null, geneMap: null, cutoff: 0, metric: "density", dotSize: 1.5,
                  segGenes: null, crossMetric: "density", drawerOpen: false,
                  genes: [], segShow: {}, pnEnr: null };
  const GENE_PALETTE = ["#e11d48", "#2563eb", "#16a34a", "#f97316", "#7c3aed", "#0891b2", "#db2777", "#ca8a04"];
  let vcExtras = null;   // dot-size + atlas-link row (VCore.addWindowExtras)

  const segColor = (lbl) => {
    const d = (state.scene.region_defaults || {})[String(lbl)];
    return (d && d.color) || "#94a3b8";
  };
  const segName = (lbl) => `Segment ${lbl}`;

  (async function init() {
    try {
      const m = await (await fetch("data/segments_manifest.json")).json();
      state.manifest = m.embryos;
      V.loadGz("data/segments_genes.json.gz").then((g) => { state.segGenes = g; if (state.drawerOpen) renderCrossGene(); }).catch(() => {});
      V.loadGz("data/pronuclei_enrichment.json.gz").then((e) => { state.pnEnr = e; if (state.segTab === "enr" || state.segTab === "only") renderSegList(); }).catch(() => {});
      countEl.textContent = `${m.embryos.length} embryos · ${m.stages.length} stages`;
      V.buildTabs(tabsEl, m.embryos, selectEmbryo, (e) => ({
        label: e.label, sub: e.stage_label,
        title: `${e.stage_label} · ${e.eid} · segments ${e.segments.join(",")}`,
      }));
      V.wireWindow(controlsEl, $("#controls-header"),
        [...controlsEl.querySelectorAll(".rz")], "segments_controls_box")
        .setResizeCb(() => { try { Plotly.Plots.resize(chartEl); } catch (_) {} });
      vcExtras = V.addWindowExtras($("#controls-body"), { defaultSize: state.dotSize, onDotSize: (s) => { state.dotSize = s; if (state.scene) render(); } });
      wireRdrawer(); wireCrossDrawer();
      geneSelect.addEventListener("change", () => {
        if (!geneSelect.value) return;
        addGene(geneSelect.value); geneSelect.value = "";
        render(); renderChart(); highlightGene(); if (state.drawerOpen) renderCrossGene();
      });
      geneChips.addEventListener("click", (e) => {
        const x = e.target.closest(".seg-gx");
        if (x) { removeGene(x.dataset.gene); render(); renderChart(); highlightGene(); if (state.drawerOpen) renderCrossGene(); return; }
        const chip = e.target.closest(".seg-gchip");
        if (chip) { state.userGene = chip.dataset.gene; renderChips(); highlightGene(); renderChart(); if (state.drawerOpen) renderCrossGene(); }
      });
      minCountEl.addEventListener("input", () => {
        state.cutoff = Math.max(0, parseInt(minCountEl.value, 10) || 0);
        if (state.scene) renderSegList();
      });
      metricModeEl.addEventListener("change", () => {
        state.metric = metricModeEl.value;
        if (state.scene) { renderSegList(); renderChart(); }
      });
    } catch (err) { showError("Failed to load manifest: " + (err.message || err)); }
  })();

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    V.markActiveTab(tabsEl, id);
    const meta = state.manifest.find((m) => m.id === id) || {};
    showLoading(`Loading ${meta.label || id}…`);
    try {
      const scene = await V.loadGz(`data/segments/${id}.json.gz`);
      if (state.currentId !== id) return;
      state.scene = scene;
      if (vcExtras) vcExtras.setAtlas(id);
      buildGeneMap(scene);
      populateGenes(scene);
      if (state.segTab !== "enr" && state.segTab !== "only" && !scene.mask_labels.includes(state.segTab)) state.segTab = scene.mask_labels[0];
      controlsEl.hidden = false; placeholder.hidden = true; rdrawer.hidden = false; drawer.hidden = false;
      renderLegend(); buildSegDotToggles(); buildSegTabs(); render(); renderChart(); renderSegList();
    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }

  // gene -> { seg -> {enrich, count} }, from the per-segment ranked lists
  function buildGeneMap(scene) {
    const map = {};
    for (const s of scene.mask_labels) {
      for (const r of (scene.ranked[String(s)] || [])) {
        (map[r.gene] = map[r.gene] || {})[s] = { enrich: r.enrich, count: r.count, ntot: r.ntot };
      }
    }
    state.geneMap = map;
  }
  function populateGenes(scene) {
    const genes = Object.keys(state.geneMap).sort((a, b) =>
      a.localeCompare(b, undefined, { sensitivity: "base", numeric: true }));
    geneSelect.innerHTML = `<option value="">＋ add a gene…</option>` +
      genes.map((g) => `<option value="${g}">${g}</option>`).join("");
    geneSelect.value = "";
    // keep only genes still present in this embryo; always plot at least one
    state.genes = state.genes.filter((x) => state.geneMap[x.gene]);
    if (!state.genes.length && genes[0]) addGene(genes[0]);
    if (!state.userGene || !state.geneMap[state.userGene]) state.userGene = state.genes.length ? state.genes[0].gene : "";
    renderChips();
  }

  // ---------- 3-D ----------
  function segMesh(scene, lbl, color, opacity) {
    const mesh = scene.region_meshes[String(lbl)]; if (!mesh) return null;
    const v = mesh.verts, f = mesh.faces, nV = v.length / 3, nF = f.length / 3;
    const x = new Array(nV), y = new Array(nV), z = new Array(nV);
    for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
    const ii = new Array(nF), jj = new Array(nF), kk = new Array(nF);
    for (let i = 0; i < nF; i++) { ii[i] = f[i * 3]; jj[i] = f[i * 3 + 1]; kk[i] = f[i * 3 + 2]; }
    return { type: "mesh3d", x, y, z, i: ii, j: jj, k: kk, color, opacity,
      name: segName(lbl), showlegend: true, flatshading: false, hoverinfo: "name",
      lighting: { ambient: 0.68, diffuse: 0.6, specular: 0.12, roughness: 0.9 }, legendrank: lbl };
  }
  // ---------- multi-gene set + per-segment dot toggles ----------
  function addGene(g) {
    if (!g || state.genes.some((x) => x.gene === g)) return;
    const used = new Set(state.genes.map((x) => x.color));
    const color = GENE_PALETTE.find((c) => !used.has(c)) || GENE_PALETTE[state.genes.length % GENE_PALETTE.length];
    state.genes.push({ gene: g, color });
    state.userGene = g; renderChips();
  }
  function removeGene(g) {
    state.genes = state.genes.filter((x) => x.gene !== g);
    if (state.userGene === g) state.userGene = state.genes.length ? state.genes[state.genes.length - 1].gene : "";
    renderChips();
  }
  function renderChips() {
    geneChips.innerHTML = state.genes.map((x) =>
      `<span class="seg-gchip${x.gene === gene() ? " current" : ""}" data-gene="${x.gene}" title="click to focus, × to remove">` +
      `<span class="seg-gdot" style="background:${x.color}"></span>${x.gene}` +
      `<button class="seg-gx" data-gene="${x.gene}" aria-label="remove ${x.gene}">×</button></span>`).join("") ||
      `<span class="seg-gchips-empty">pick a gene above to plot it</span>`;
  }
  const gene = () => state.userGene;                          // the focused gene (right-drawer highlight)
  function buildSegDotToggles() {
    const s = state.scene; if (!s) return;
    segDotToggles.innerHTML = s.mask_labels.map((l) =>
      `<label class="seg-dtog" title="show/hide the plotted genes' dots that fall in ${segName(l)}">` +
      `<input type="checkbox" data-seg="${l}" ${state.segShow[String(l)] === false ? "" : "checked"}>` +
      `<span class="seg-dot" style="background:${segColor(l)}"></span>${l}</label>`).join("");
    segDotToggles.querySelectorAll("input").forEach((c) => c.addEventListener("change", () => {
      state.segShow[c.dataset.seg] = c.checked; render();
    }));
  }
  function render() {
    const s = state.scene; if (!s) return;
    const traces = [];
    for (const lbl of s.mask_labels) { const t = segMesh(s, lbl, segColor(lbl), 0.1); if (t) traces.push(t); }
    const zs = s.z_scale;
    // each plotted gene's transcript dots, coloured by gene, filtered to the shown segments
    for (const { gene: g, color } of state.genes) {
      const tx = s.transcripts ? s.transcripts[g] : null;
      if (!tx || !tx.x.length) continue;
      const xi = [], yi = [], zi = [];
      for (let i = 0; i < tx.x.length; i++) {
        if (state.segShow[String(tx.s[i])] === false) continue;
        xi.push(tx.x[i]); yi.push(tx.y[i]); zi.push(tx.gz[i] * zs);
      }
      if (!xi.length) continue;
      traces.push({ type: "scatter3d", mode: "markers", name: `${g} · ${xi.length}`,
        x: xi, y: yi, z: zi, marker: { size: state.dotSize, color, opacity: 0.85, line: { width: 0 } },
        hovertemplate: `${g}<extra></extra>`, legendrank: 20000 });
    }
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  // ---------- floating-window chart: gene enrichment across segments ----------
  // Layered bar across segments: a wide translucent grey bar for ALL plotted genes combined,
  // with skinny per-gene coloured bars sitting inside it.
  function renderChart() {
    const s = state.scene; if (!s) return;
    const genes = state.genes.filter((x) => state.geneMap[x.gene]);
    chartSub.textContent = genes.length ? `· ${genes.map((x) => x.gene).join(", ")}` : "";
    const frac = isFrac(), labels = s.mask_labels, xcat = labels.map(segName);
    const vtot = s.segments.reduce((a, q) => a + (q.volume || 0), 0);
    const volOf = (l) => (s.segments.find((q) => q.label === l) || {}).volume || 0;
    const gv = (g, l) => { const r = (state.geneMap[g] || {})[l]; return r ? metricVal(r) : 0; };
    const combined = (l) => {                                  // all plotted genes pooled as one
      let cs = 0, ct = 0;
      for (const { gene: g } of genes) { const r = (state.geneMap[g] || {})[l]; if (r) { cs += r.count; ct += r.ntot; } }
      if (!ct) return 0;
      if (frac) return cs / ct;
      const vol = volOf(l); return vol ? (cs / vol) / (ct / vtot) : 0;
    };
    const traces = [];
    traces.push({ type: "bar", name: "all genes", x: xcat, y: labels.map(combined),
      marker: { color: "rgba(100,110,125,0.28)", line: { width: 0 } }, width: 0.74, offset: -0.37,
      hovertemplate: frac ? "%{x}: %{y:.0%} combined<extra></extra>" : "%{x}: %{y:.2f}× combined<extra></extra>" });
    const n = Math.max(1, genes.length), bw = 0.66 / n;
    genes.forEach(({ gene: g, color }, gi) => {
      traces.push({ type: "bar", name: g, x: xcat, y: labels.map((l) => gv(g, l)),
        marker: { color }, width: bw, offset: -0.33 + gi * bw,
        hovertemplate: `${g} · %{x}: ${frac ? "%{y:.0%}" : "%{y:.2f}×"}<extra></extra>` });
    });
    if (!frac) traces.push({ type: "scatter", mode: "lines", x: [xcat[0], xcat[xcat.length - 1]], y: [1, 1],
      line: { color: "#94a3b8", width: 1, dash: "dot" }, hoverinfo: "skip", showlegend: false });
    plotInto(chartEl, traces, {
      barmode: "overlay", margin: { l: 38, r: 6, t: genes.length > 1 ? 22 : 6, b: 34 }, height: 158,
      yaxis: { tickfont: { size: 9 }, gridcolor: "#eef1f5", fixedrange: true,
        title: { text: frac ? "fraction" : "fold", font: { size: 9 } }, rangemode: "tozero",
        tickformat: frac ? ".0%" : undefined },
      xaxis: { tickfont: { size: 9 }, fixedrange: true, tickangle: -30, type: "category" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      showlegend: genes.length > 1, legend: { orientation: "h", font: { size: 8 }, y: 1.02, x: 0, yanchor: "bottom" },
    });
    const gc = gene();
    const gm = state.geneMap[gc] || {};
    const segs = Object.entries(gm).sort((a, b) => metricVal(b[1]) - metricVal(a[1]));
    if (segs.length) {
      const [top, ti] = segs[0];
      chartReadout.innerHTML = `<div><b>${gc}</b> — ${ti.ntot} transcripts · ${genes.length} gene${genes.length === 1 ? "" : "s"} plotted</div>` +
        `<div>${frac ? "most concentrated" : "most enriched"} in <b>${segName(top)}</b>: <b>${metricStr(ti)}</b> (${ti.count} here)</div>`;
    } else chartReadout.innerHTML = `<div>${genes.length} gene${genes.length === 1 ? "" : "s"} plotted.</div>`;
  }

  function renderLegend() {
    const s = state.scene;
    segLegend.innerHTML = s.mask_labels.map((l) => {
      const vol = (s.segments.find((q) => q.label === l) || {}).volume || 0;
      return `<span class="seg-chip" title="volume ${Math.round(vol).toLocaleString()} µm³">` +
             `<span class="seg-dot" style="background:${segColor(l)}"></span>${segName(l)}</span>`;
    }).join("");
  }

  // ---------- right drawer: ranked genes per segment ----------
  function buildSegTabs() {
    const s = state.scene;
    let html = s.mask_labels.map((l) =>
      `<button class="rtab${l === state.segTab ? " active" : ""}" data-seg="${l}">Seg ${l}</button>`).join("");
    // cross-embryo pronuclei-enrichment lists (all zygotes)
    html += `<button class="rtab${state.segTab === "enr" ? " active" : ""}" data-seg="enr" title="Genes enriched in the pronuclei, across all zygotes">Enriched in PN</button>`;
    html += `<button class="rtab${state.segTab === "only" ? " active" : ""}" data-seg="only" title="Genes only in the pronuclei, across all zygotes">Only in PN</button>`;
    rtabsEl.innerHTML = html;
    rtabsEl.querySelectorAll(".rtab").forEach((b) => b.addEventListener("click", () => {
      const v = b.dataset.seg;
      state.segTab = (v === "enr" || v === "only") ? v : parseInt(v, 10);
      rtabsEl.querySelectorAll(".rtab").forEach((x) => x.classList.toggle("active", x === b));
      renderSegList();
    }));
  }
  const isFrac = () => state.metric === "fraction";
  const metricVal = (r) => isFrac() ? r.count / r.ntot : r.enrich;               // fraction of gene in seg, or density fold
  const metricStr = (r) => isFrac() ? `${(100 * r.count / r.ntot).toFixed(0)}%` : `${r.enrich}×`;
  function updateDesc() {
    rdrawerDesc.innerHTML = isFrac()
      ? "Fraction of the gene's transcripts that fall in this segment (segment count ÷ the gene's total across segments) — no volume correction."
      : "Density fold-change vs the embryo-wide average (a gene entirely in one segment tops out at V<sub>total</sub>/V<sub>segment</sub>).";
  }
  function renderSegList() {
    const s = state.scene; if (!s) return;
    if (state.segTab === "enr" || state.segTab === "only") { renderPnList(state.segTab); return; }
    segListEl.classList.remove("seg-pnlist");
    updateDesc();
    const cut = state.cutoff || 0, frac = isFrac();
    // cutoff = minimum total transcripts of the gene in this embryo (ntot); re-sort by the chosen metric
    const rows = (s.ranked[String(state.segTab)] || []).filter((r) => r.ntot > cut)
      .sort((a, b) => metricVal(b) - metricVal(a));
    const cur = gene();
    const vol = (s.segments.find((q) => q.label === state.segTab) || {}).volume || 0;
    const note = cut > 0 ? `> ${cut} transcripts in embryo` : "≥3 transcripts in segment";
    let html = `<div class="best-plane-note">${segName(state.segTab)} · volume ` +
      `<b>${Math.round(vol).toLocaleString()}</b> µm³ · ${rows.length} genes (${note})</div>`;
    html += `<div class="best-head"><span></span><span>gene</span><span>${frac ? "% here" : "fold"}</span><span>count</span></div>`;
    html += rows.map((r, i) =>
      `<div class="best-row seg-row${r.gene === cur ? " current" : ""}" data-gene="${r.gene}" ` +
      `title="${r.gene}: ${r.count} of ${r.ntot} transcripts here${frac ? ` · ${r.enrich}× density` : ` · ${(100 * r.count / r.ntot).toFixed(0)}% of the gene`}">` +
      `<span class="best-num">${i + 1}</span>` +
      `<span class="best-gene">${r.gene}</span>` +
      `<span class="best-real">${metricStr(r)}</span>` +
      `<span class="best-p">${r.count}</span></div>`).join("");
    segListEl.innerHTML = html || `<div class="best-plane-note">No genes above the threshold.</div>`;
  }
  function highlightGene() {
    const cur = gene();
    segListEl.querySelectorAll(".seg-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === cur));
    renderChips();
  }
  // Cross-embryo pronuclei enrichment (from build_pronuclei_enrichment.py), shown in the right drawer.
  function renderPnList(which) {
    const enr = state.pnEnr;
    rdrawerDesc.innerHTML = which === "enr"
      ? "Across all zygotes: genes whose transcripts are <b>≥1.5× denser</b> in the two pronuclei than the cell average (≥10 in-cell transcripts). Fold = pronuclei density ÷ cell-average; p = the random-scatter null. Click to open that zygote &amp; gene."
      : "Genes with <b>≥75%</b> of their in-cell transcripts inside the pronuclei (≥5) — effectively <i>only</i> there (pronuclei ≈5% of the cell volume). Click to open that zygote &amp; gene.";
    segListEl.classList.add("seg-pnlist");
    if (!enr) { segListEl.innerHTML = `<div class="best-plane-note">Loading pronuclei enrichment…</div>`; return; }
    const all = (which === "enr" ? enr.enriched : enr.onlyPn) || [];
    const cur = gene();
    const fmtP = (p) => (p == null ? "—" : p < 1e-4 ? p.toExponential(1) : p.toFixed(4));
    let html = `<div class="best-plane-note"><b>${all.length}</b> gene·zygote hits · ${which === "enr" ? "fold ≥ 1.5, n ≥ 10" : "≥ 75% in pronuclei"}</div>`;
    html += `<div class="best-head"><span></span><span>gene · zygote</span><span>${which === "enr" ? "fold" : "in PN"}</span><span>p</span></div>`;
    html += all.slice(0, 250).map((r, i) =>
      `<div class="best-row seg-row seg-pnrow${r.gene === cur ? " current" : ""}" data-gene="${r.gene}" data-zid="${r.id}" ` +
      `title="${r.gene} · ${r.label} · ${r.npn}/${r.n} in pronuclei${r.fold != null ? " · fold " + r.fold + "×" : ""}${r.frac != null ? " · " + Math.round(r.frac * 100) + "%" : ""} · p ${fmtP(r.p)}">` +
      `<span class="best-num">${i + 1}</span>` +
      `<span class="best-gene"><b>${r.gene}</b> <span class="seg-pnzyg">${r.label}</span></span>` +
      `<span class="best-real">${which === "enr" ? r.fold.toFixed(2) + "×" : Math.round(r.frac * 100) + "%"}</span>` +
      `<span class="best-p${r.p <= 0.05 ? " sig" : ""}">${fmtP(r.p)}</span></div>`).join("");
    segListEl.innerHTML = html || `<div class="best-plane-note">None found.</div>`;
  }
  // ---------- bottom drawer: the selected gene's enrichment across every embryo it appears in ----------
  const STAGE_RANK = (s) => /oocyte/i.test(s) ? 0 : /zygote/i.test(s) ? 1 : /early/i.test(s) ? 2 : /late/i.test(s) ? 3 : 4;
  const SEG_PALETTE = ["#8dd3c7", "#ffd92f", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5"];
  function renderCrossGene() {
    if (!crossChart) return;
    const g = gene(), sg = state.segGenes;
    if (crossGeneEl) crossGeneEl.textContent = g ? `· ${g}` : "";
    if (!sg) { crossChart.innerHTML = `<div class="seg-cross-empty">Loading cross-embryo data…</div>`; return; }
    if (!g || !sg.genes[g]) {
      crossChart.innerHTML = `<div class="seg-cross-empty"><b>${g || "—"}</b> is not present in any embryo.</div>`;
      if (crossNote) crossNote.textContent = ""; return;
    }
    const entries = sg.genes[g].map(([idx, rows]) => ({ e: sg.embInfo[idx], rows }))
      .sort((a, b) => STAGE_RANK(a.e.stage) - STAGE_RANK(b.e.stage) || a.e.label.localeCompare(b.e.label));
    const frac = state.crossMetric === "fraction";
    const xlabels = entries.map((x) => x.e.label);
    const segNums = [...new Set(entries.flatMap((x) => x.rows.map((r) => r[0])))].sort((a, b) => a - b);
    const traces = segNums.map((seg) => ({
      type: "bar", name: `Segment ${seg}`, x: xlabels,
      y: entries.map((x) => { const r = x.rows.find((rr) => rr[0] === seg); if (!r) return null;
        return frac ? (r[3] ? r[2] / r[3] : 0) : r[1]; }),                        // [seg, enrich, count, ntot]
      marker: { color: SEG_PALETTE[(seg - 1) % SEG_PALETTE.length] },
      hovertemplate: `%{x} · seg ${seg}<br>${frac ? "%{y:.0%} of the gene" : "%{y}× enrichment"}<extra></extra>`,
    }));
    const shapes = frac ? [] : [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: 1, y1: 1,
      line: { color: "#94a3b8", width: 1, dash: "dot" } }];
    if (crossNote) crossNote.innerHTML = `<b>${g}</b> · ${entries.length} embryos · bars grouped by segment ` +
      `— segment meanings vary by stage, hover for detail`;
    plotInto(crossChart, traces, {
      barmode: "group", bargap: 0.25, bargroupgap: 0.04, shapes,
      margin: { l: 46, r: 8, t: 6, b: 104 }, height: crossChart.clientHeight || 300,
      showlegend: true, legend: { orientation: "h", y: 1.06, x: 0, font: { size: 9 } },
      xaxis: { type: "category", tickangle: -55, tickfont: { size: 7 }, automargin: true, fixedrange: false },
      yaxis: { title: { text: frac ? "fraction of the gene" : "fold enrichment", font: { size: 10 } },
        tickfont: { size: 9 }, gridcolor: "#eef1f5", rangemode: "tozero", tickformat: frac ? ".0%" : undefined, fixedrange: false },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, { responsive: true, displayModeBar: false, scrollZoom: true });
  }
  function setDrawerOpen(open) {
    state.drawerOpen = open;
    drawer.dataset.open = open ? "true" : "false";
    drawerHandle.setAttribute("aria-expanded", String(open));
    if (open) { renderCrossGene(); requestAnimationFrame(() => { try { Plotly.Plots.resize(crossChart); } catch (_) {} }); }
  }
  function wireCrossDrawer() {
    drawerHandle.addEventListener("click", () => setDrawerOpen(!state.drawerOpen));
    let sh = 0; const rz = $("#drawer-resize");
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return;
      drawer.style.setProperty("--drawer-h", Math.max(220, Math.min(window.innerHeight - 120, sh + (rz._d.y - e.clientY))) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {}
      requestAnimationFrame(() => { try { Plotly.Plots.resize(crossChart); } catch (_) {} }); } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
    if (crossMetricEl) crossMetricEl.addEventListener("change", () => { state.crossMetric = crossMetricEl.value; renderCrossGene(); });
  }
  function wireRdrawer() {
    // pull-open + resize via the handle (tap toggles); grabber fine-tunes width
    wireHandleDrag(rdrawer, rdrawerHandle, {
      computeSize: (e) => window.innerWidth - e.clientX,
      clampSize: (px) => Math.max(240, Math.min(window.innerWidth - 80, px)),
      applySize: (px) => rdrawer.style.setProperty("--rdrawer-w", px + "px"),
      setOpen: (open) => { rdrawer.dataset.open = open ? "true" : "false"; rdrawerHandle.setAttribute("aria-expanded", String(open)); },
    });
    let sw = 0; const rrz = $("#rdrawer-resize");
    rrz.addEventListener("pointerdown", (e) => { sw = rdrawer.getBoundingClientRect().width; rrz._d = { x: e.clientX }; rrz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rrz.addEventListener("pointermove", (e) => { if (!rrz._d) return;
      const w = Math.max(240, Math.min(window.innerWidth - 80, sw - (e.clientX - rrz._d.x)));
      rdrawer.style.setProperty("--rdrawer-w", w + "px"); });
    const rend = (e) => { if (rrz._d) { rrz._d = null; try { rrz.releasePointerCapture(e.pointerId); } catch (_) {} } };
    rrz.addEventListener("pointerup", rend); rrz.addEventListener("pointercancel", rend);
    segListEl.addEventListener("click", (e) => {
      const row = e.target.closest(".seg-row"); if (!row) return;
      const g = row.dataset.gene, zid = row.dataset.zid;
      if (zid) {                                          // a pronuclei-enrichment hit → open that zygote + gene
        const segId = "Zygote__" + zid;
        if (segId !== state.currentId && state.manifest.some((m) => m.id === segId)) {
          selectEmbryo(segId).then(() => { state.genes = []; addGene(g); render(); renderChart(); highlightGene(); });
        } else { addGene(g); render(); renderChart(); highlightGene(); }
        return;
      }
      addGene(g); render(); renderChart(); highlightGene(); if (state.drawerOpen) renderCrossGene();
    });
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
      if (moved) drawerEl.classList.remove("dragging");
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
