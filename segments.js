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
  // Viridis colormap for enrichment (low = dark/blue, high = yellow).
  const VIRIDIS = [[68, 1, 84], [59, 81, 139], [33, 144, 141], [92, 200, 99], [253, 231, 37]];
  function viridis(t) {
    t = Math.max(0, Math.min(1, t));
    const x = t * (VIRIDIS.length - 1), i = Math.min(VIRIDIS.length - 2, Math.floor(x)), f = x - i;
    const a = VIRIDIS[i], b = VIRIDIS[i + 1];
    return `rgb(${Math.round(a[0] + (b[0] - a[0]) * f)},${Math.round(a[1] + (b[1] - a[1]) * f)},${Math.round(a[2] + (b[2] - a[2]) * f)})`;
  }

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const geneSelect = $("#gene-select"), chartEl = $("#chart"), chartSub = $("#chart-sub");
  const chartReadout = $("#chart-readout"), segLegend = $("#seg-legend");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle");
  const rtabsEl = $("#rtabs"), segListEl = $("#seg-list");

  const state = { manifest: [], currentId: null, scene: null, userGene: null,
                  segTab: null, geneMap: null };

  const segColor = (lbl) => {
    const d = (state.scene.region_defaults || {})[String(lbl)];
    return (d && d.color) || "#94a3b8";
  };
  const segName = (lbl) => `Segment ${lbl}`;

  (async function init() {
    try {
      const m = await (await fetch("data/segments_manifest.json")).json();
      state.manifest = m.embryos;
      countEl.textContent = `${m.embryos.length} embryos · ${m.stages.length} stages`;
      V.buildTabs(tabsEl, m.embryos, selectEmbryo, (e) => ({
        label: e.label, sub: e.stage_label,
        title: `${e.stage_label} · ${e.eid} · segments ${e.segments.join(",")}`,
      }));
      V.wireWindow(controlsEl, $("#controls-header"),
        [...controlsEl.querySelectorAll(".rz")], "segments_controls_box")
        .setResizeCb(() => { try { Plotly.Plots.resize(chartEl); } catch (_) {} });
      wireRdrawer();
      geneSelect.addEventListener("change", () => { state.userGene = geneSelect.value; render(); renderChart(); highlightGene(); });
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
      buildGeneMap(scene);
      populateGenes(scene);
      if (!scene.mask_labels.includes(state.segTab)) state.segTab = scene.mask_labels[0];
      controlsEl.hidden = false; placeholder.hidden = true; rdrawer.hidden = false;
      renderLegend(); buildSegTabs(); render(); renderChart(); renderSegList();
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
    geneSelect.innerHTML = genes.map((g) => `<option value="${g}">${g}</option>`).join("");
    geneSelect.value = (state.userGene && state.geneMap[state.userGene]) ? state.userGene : (genes[0] || "");
    state.userGene = geneSelect.value;
  }
  const gene = () => geneSelect.value;

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
  function render() {
    const s = state.scene; if (!s) return;
    const g = gene(), gm = state.geneMap[g] || {};
    const enrichVals = s.mask_labels.map((l) => (gm[l] ? gm[l].enrich : null)).filter((x) => x != null);
    const emax = enrichVals.length ? Math.max(...enrichVals) : 1;
    const emin = enrichVals.length ? Math.min(...enrichVals) : 0;
    const span = (emax - emin) || 1;
    const traces = [];
    for (const lbl of s.mask_labels) {
      const info = gm[lbl];
      let color, opacity;
      if (g && info) {                       // colour by this gene's enrichment
        color = viridis((info.enrich - emin) / span);
        opacity = 0.35 + 0.5 * ((info.enrich - emin) / span);
      } else if (g) {                        // gene not in this segment (< min count)
        color = "#c7ccd6"; opacity = 0.12;
      } else {
        color = segColor(lbl); opacity = 0.28;
      }
      const t = segMesh(s, lbl, color, opacity);
      if (t) { t.name = `${segName(lbl)}${info ? ` · ${info.enrich}×` : ""}`; traces.push(t); }
    }
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  // ---------- floating-window chart: gene enrichment across segments ----------
  function renderChart() {
    const s = state.scene; if (!s) return;
    const g = gene(); chartSub.textContent = g ? `· ${g}` : "";
    const gm = state.geneMap[g] || {};
    const labels = s.mask_labels;
    const y = labels.map((l) => (gm[l] ? gm[l].enrich : 0));
    const colors = labels.map((l) => (gm[l] ? segColor(l) : "#d5d9e2"));
    plotInto(chartEl, [
      { type: "bar", x: labels.map(segName), y, marker: { color: colors },
        hovertemplate: "%{x}: %{y}× enrichment<extra></extra>" },
      { type: "scatter", mode: "lines", x: [segName(labels[0]), segName(labels[labels.length - 1])],
        y: [1, 1], line: { color: "#94a3b8", width: 1, dash: "dot" }, hoverinfo: "skip" },
    ], {
      margin: { l: 34, r: 6, t: 6, b: 34 }, height: 150,
      yaxis: { tickfont: { size: 9 }, gridcolor: "#eef1f5", fixedrange: true, title: { text: "fold", font: { size: 9 } }, rangemode: "tozero" },
      xaxis: { tickfont: { size: 9 }, fixedrange: true, tickangle: -30 }, bargap: 0.45,
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", showlegend: false,
    });
    const segs = Object.entries(gm).sort((a, b) => b[1].enrich - a[1].enrich);
    if (segs.length) {
      const [top, ti] = segs[0];
      chartReadout.innerHTML = `<div><b>${g}</b> — ${ti.ntot} transcripts in segments</div>` +
        `<div>most enriched in <b>${segName(top)}</b>: <b>${ti.enrich}×</b> (${ti.count} here)</div>`;
    } else chartReadout.innerHTML = `<div><b>${g}</b> — below the count threshold in every segment.</div>`;
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
    rtabsEl.innerHTML = s.mask_labels.map((l) =>
      `<button class="rtab${l === state.segTab ? " active" : ""}" data-seg="${l}">Seg ${l}</button>`).join("");
    rtabsEl.querySelectorAll(".rtab").forEach((b) => b.addEventListener("click", () => {
      state.segTab = parseInt(b.dataset.seg, 10);
      rtabsEl.querySelectorAll(".rtab").forEach((x) => x.classList.toggle("active", x === b));
      renderSegList();
    }));
  }
  function renderSegList() {
    const s = state.scene; if (!s) return;
    const rows = s.ranked[String(state.segTab)] || [];
    const cur = gene();
    const vol = (s.segments.find((q) => q.label === state.segTab) || {}).volume || 0;
    let html = `<div class="best-plane-note">${segName(state.segTab)} · volume ` +
      `<b>${Math.round(vol).toLocaleString()}</b> µm³ · ${rows.length} genes (≥3 transcripts)</div>`;
    html += `<div class="best-head"><span></span><span>gene</span><span>fold</span><span>count</span></div>`;
    html += rows.map((r, i) =>
      `<div class="best-row seg-row${r.gene === cur ? " current" : ""}" data-gene="${r.gene}" ` +
      `title="${r.gene}: ${r.count} of ${r.ntot} transcripts here">` +
      `<span class="best-num">${i + 1}</span>` +
      `<span class="best-gene">${r.gene}</span>` +
      `<span class="best-real">${r.enrich}×</span>` +
      `<span class="best-p">${r.count}</span></div>`).join("");
    segListEl.innerHTML = html || `<div class="best-plane-note">No genes above the threshold.</div>`;
  }
  function highlightGene() {
    const cur = gene();
    segListEl.querySelectorAll(".seg-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === cur));
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
      const g = row.dataset.gene;
      if (state.geneMap[g]) { state.userGene = g; geneSelect.value = g; render(); renderChart(); highlightGene(); }
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
