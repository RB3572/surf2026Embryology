/* Stage Expression Explorer.
 *
 * Modelled on Harry Wang's Stage Expression Explorer (same stages, axis labels, stats),
 * rebuilt in this site's style: the MAIN view is a grid of 3-D sample views (cell
 * silhouette + located transcripts) for every sample that detected the gene, and the
 * BOTTOM DRAWER carries the per-stage plots — mean bars, box, or violin — with Welch
 * pairwise p-values and a per-stage summary.
 *
 * Data: data/stage_expr.json.gz (built by build_stage_expr.py). Raw counts reproduce
 * segments_genes.json.gz exactly; CPM and volume-density are derived on the fly.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;

  const AGG_URL = "data/stage_expr.json.gz";
  const STAGE_ORDER = ["Zygote", "Early2Cell", "Late2Cell"];

  const els = {
    count: $("#samp-count"), form: $("#gene-form"), search: $("#gene-search"),
    datalist: $("#gene-options"), legend: $("#stage-legend"),
    gridHead: $("#grid-head"), gridGene: $("#grid-gene"), gridSub: $("#grid-sub"),
    grid: $("#stage-grid"), placeholder: $("#placeholder"), loading: $("#loading"),
    drawer: $("#drawer"), drawerHandle: $("#drawer-handle"), drawerGene: $("#drawer-gene"),
    statStrip: $("#stat-strip"), chart: $("#expr-chart"), summary: $("#summary-band"),
    exportCsv: $("#export-csv"), logToggle: $("#log-toggle"),
    focus: $("#focus-modal"), focusTitle: $("#focus-title"), focusPlot: $("#focus-plot"),
    focusClose: $("#focus-close"), focusBackdrop: $("#focus-backdrop"), focusLoading: $("#focus-loading"),
  };

  const state = {
    agg: null, gene: null, mode: "bar", norm: "raw", log: true,
    sceneCache: {},
  };

  // ---------- boot ----------
  (async function init() {
    try {
      state.agg = await V.loadGz(AGG_URL);
    } catch (e) {
      els.placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_stage_expr.py</code> to generate ` +
        `data/stage_expr.json.gz.</div></div>`;
      return;
    }
    const a = state.agg;
    els.count.textContent = `${a.samples.length} samples · ${a.stages.length} stages · ${a.genes_all.length} genes`;
    els.datalist.innerHTML = a.genes_all.map((g) => `<option value="${g}"></option>`).join("");
    els.legend.innerHTML = STAGE_ORDER.map((s) =>
      `<span><i style="background:${a.stage_colors[s]}"></i>${a.stage_labels[s]}</span>`).join("");
    wireControls();
    const start = a.gene_counts["Pard3"] ? "Pard3" : a.genes_all[0];
    els.search.value = start;
    selectGene(start);
  })();

  function wireControls() {
    els.form.addEventListener("submit", (e) => {
      e.preventDefault();
      const g = (els.search.value || "").trim();
      if (state.agg.gene_counts[g]) selectGene(g);
      else flashInvalid();
    });
    els.search.addEventListener("change", () => {
      const g = (els.search.value || "").trim();
      if (state.agg.gene_counts[g]) selectGene(g);
    });
    document.querySelectorAll("[data-mode]").forEach((b) =>
      b.addEventListener("click", () => setSeg("mode", b.dataset.mode, "[data-mode]", b)));
    document.querySelectorAll("[data-norm]").forEach((b) =>
      b.addEventListener("click", () => setSeg("norm", b.dataset.norm, "[data-norm]", b)));
    els.logToggle.addEventListener("change", () => { state.log = els.logToggle.checked; renderDrawer(); });
    els.exportCsv.addEventListener("click", exportCsv);
    // focus modal
    els.focusClose.addEventListener("click", closeFocus);
    els.focusBackdrop.addEventListener("click", closeFocus);
    document.addEventListener("keydown", (e) => { if (e.key === "Escape") closeFocus(); });
    wireDrawer();
  }
  function setSeg(key, val, sel, btn) {
    state[key] = val;
    document.querySelectorAll(sel).forEach((b) => b.classList.toggle("active", b === btn));
    if (key === "norm") renderDrawer();       // norm only affects the plots
    else renderDrawer();
  }
  function flashInvalid() {
    els.search.style.borderColor = "#dc2626";
    setTimeout(() => { els.search.style.borderColor = ""; }, 900);
  }

  // ---------- per-sample value + grouping ----------
  const S = () => state.agg.samples;
  function rawCount(si) { const m = state.agg.gene_counts[state.gene]; return m ? (m[String(si)] ?? null) : null; }
  function value(si) {
    const c = rawCount(si); if (c == null) return null;
    const s = S()[si];
    if (state.norm === "cpm") return s.total_tx ? c / s.total_tx * 1e6 : 0;
    if (state.norm === "density") return s.vol ? c / s.vol * 1e9 : 0;    // per 1e9 voxels, readable
    return c;
  }
  const NORM_UNIT = { raw: "Transcript count", cpm: "Transcripts per million (CPM)",
                      density: "Transcripts per 10⁹ voxels" };

  function byStage() {
    const m = state.agg.gene_counts[state.gene] || {};
    const out = { Zygote: [], Early2Cell: [], Late2Cell: [] };
    Object.keys(m).forEach((si) => {
      const i = +si, st = S()[i].stage;
      if (out[st]) out[st].push(i);
    });
    for (const st of STAGE_ORDER)
      out[st].sort((a, b) => V.cmpEmbryo(S()[a], S()[b]));
    return out;
  }

  // ---------- select gene ----------
  function selectGene(g) {
    state.gene = g;
    els.search.value = g;
    els.placeholder.hidden = true;
    els.gridHead.hidden = false;
    els.gridGene.textContent = g;
    if (els.drawer.hidden) { els.drawer.hidden = false; }
    renderGrid();
    if (state.drawerOpen) renderDrawer(); else openDrawer(true);
  }

  // ---------- grid of 3-D thumbnails ----------
  function renderGrid() {
    const groups = byStage();
    const a = state.agg;
    const total = STAGE_ORDER.reduce((n, s) => n + groups[s].length, 0);
    els.gridSub.textContent = `· ${total} samples`;
    const frag = document.createDocumentFragment();
    for (const st of STAGE_ORDER) {
      const ids = groups[st];
      if (!ids.length) continue;
      const block = el("div", "stage-block");
      const head = el("div", "stage-block-head");
      head.innerHTML = `<h3 style="color:${a.stage_colors[st]}">${a.stage_labels[st]}</h3>` +
        `<span class="sb-n">${ids.length} sample${ids.length === 1 ? "" : "s"}</span>` +
        `<span class="sb-bar" style="background:${a.stage_colors[st]}"></span>`;
      block.appendChild(head);
      const row = el("div", "thumb-row");
      ids.forEach((si) => row.appendChild(makeThumb(si, st)));
      block.appendChild(row);
      frag.appendChild(block);
    }
    els.grid.innerHTML = "";
    els.grid.appendChild(frag);
  }

  function makeThumb(si, st) {
    const a = state.agg, s = S()[si], color = a.stage_colors[st];
    const card = el("div", "thumb");
    card.style.borderTopColor = color;
    const cv = document.createElement("canvas");
    cv.width = 256; cv.height = 256;
    card.appendChild(cv);
    const rot = el("div", "t-rotate"); rot.textContent = "rotate ⤢"; card.appendChild(rot);
    const meta = el("div", "t-meta");
    const c = rawCount(si);
    meta.innerHTML = `<div class="t-label" title="${s.id}">${s.label}</div>` +
      `<div class="t-count">${c == null ? "—" : c.toLocaleString()} transcripts</div>`;
    card.appendChild(meta);
    drawThumb(cv, s, (a.gene_tx[state.gene] || {})[String(si)], color);
    card.addEventListener("click", () => openFocus(si));
    return card;
  }

  // silhouette (hull) + located transcript dots, both in the aggregate's 0..1000 frame
  function drawThumb(cv, s, dots, color) {
    const ctx = cv.getContext("2d");
    const W = cv.width, H = cv.height, pad = 14;
    const map = ([x, y]) => [pad + x / 1000 * (W - 2 * pad), pad + (1 - y / 1000) * (H - 2 * pad)];
    ctx.clearRect(0, 0, W, H);
    if (s.hull && s.hull.length > 2) {
      ctx.beginPath();
      s.hull.forEach((p, i) => { const q = map(p); i ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]); });
      ctx.closePath();
      ctx.fillStyle = hexA(color, 0.09); ctx.fill();
      ctx.lineWidth = 1.6; ctx.strokeStyle = hexA(color, 0.55); ctx.stroke();
    }
    if (dots && dots.length) {
      ctx.fillStyle = hexA(color, 0.9);
      for (const p of dots) { const q = map(p); ctx.beginPath(); ctx.arc(q[0], q[1], 2.1, 0, 7); ctx.fill(); }
    } else {
      ctx.fillStyle = "#9aa3b2"; ctx.font = "12px -apple-system, sans-serif"; ctx.textAlign = "center";
      ctx.fillText("not detected", W / 2, H / 2);
    }
  }

  // ---------- focused interactive 3-D ----------
  async function openFocus(si) {
    const s = S()[si];
    els.focus.hidden = false;
    els.focusTitle.innerHTML = `${s.label} <span class="ft-sub">${state.agg.stage_labels[s.stage]} · ${state.gene}</span>`;
    els.focusLoading.hidden = false;
    Plotly.purge(els.focusPlot);
    let scene = state.sceneCache[s.id];
    try {
      if (!scene) { scene = await V.loadGz(`data/segments/${s.id}.json.gz`); state.sceneCache[s.id] = scene; }
    } catch (e) { els.focusLoading.innerHTML = "<span>Could not load scene.</span>"; return; }
    if (els.focus.hidden) return;                // closed while loading
    const traces = V.bodyTraces(scene);
    const tx = (scene.transcripts || {})[state.gene];
    if (tx && tx.x) {
      const zs = scene.z_scale || 1;
      traces.push({
        type: "scatter3d", mode: "markers", name: `${state.gene} (${tx.x.length.toLocaleString()})`,
        x: tx.x, y: tx.y, z: tx.gz.map((v) => v * zs),
        marker: { size: 2.4, color: state.agg.stage_colors[s.stage], opacity: V.DOT_OPACITY, line: { width: 0 } },
        hovertemplate: `${state.gene}<extra></extra>`,
      });
    }
    els.focusLoading.hidden = true;
    Plotly.react(els.focusPlot, traces, V.sceneLayout(scene.extents, s.id + state.gene), V.plotConfig);
  }
  function closeFocus() { if (els.focus.hidden) return; els.focus.hidden = true; Plotly.purge(els.focusPlot); }

  // ---------- bottom drawer: per-stage plots ----------
  function openDrawer(open) {
    state.drawerOpen = open;
    els.drawer.dataset.open = open ? "true" : "false";
    els.drawerHandle.setAttribute("aria-expanded", String(open));
    if (open) renderDrawer();
  }
  function wireDrawer() {
    els.drawerHandle.addEventListener("click", () => openDrawer(els.drawer.dataset.open !== "true"));
    // drag-to-resize height
    const rz = $("#drawer-resize");
    let sh = 0, dy = 0, dragging = false;
    rz.addEventListener("pointerdown", (e) => { dragging = true; sh = els.drawer.getBoundingClientRect().height; dy = e.clientY; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!dragging) return;
      const h = Math.max(220, Math.min(window.innerHeight - 90, sh + (dy - e.clientY)));
      els.drawer.style.setProperty("--drawer-h", h + "px"); });
    const end = () => { dragging = false; };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }

  function renderDrawer() {
    if (els.drawer.hidden || !state.gene) return;
    els.drawerGene.textContent = state.gene;
    const groups = byStage();
    const data = {};                            // stage -> [values]
    const detail = {};                          // stage -> [{label,id,value,raw}]
    for (const st of STAGE_ORDER) {
      const arr = [], det = [];
      for (const si of groups[st]) {
        const v = value(si); if (v == null) continue;
        arr.push(v); det.push({ label: S()[si].label, id: S()[si].id, value: v, raw: rawCount(si) });
      }
      data[st] = arr; detail[st] = det;
    }
    renderStats(data);
    drawChart(data, detail);
    renderSummary(data);
    state._csv = detail;
  }

  // ---- statistics strip (Welch pairwise) ----
  function renderStats(data) {
    const a = state.agg;
    const pairs = [["Zygote", "Early2Cell"], ["Zygote", "Late2Cell"], ["Early2Cell", "Late2Cell"]];
    let html = `<span class="st-method">Welch two-sided p-value</span>`;
    for (const [x, y] of pairs) {
      const p = welchP(data[x], data[y]);
      const txt = p == null ? "—" : (p < 0.0001 ? "<0.0001" : p.toFixed(p < 0.01 ? 4 : 3));
      const sig = p != null && p < 0.05 ? " st-sig" : "";
      html += `<span class="st-pair">${a.stage_labels[x]} vs ${a.stage_labels[y]}` +
        `<strong class="${sig.trim()}">${txt}</strong></span>`;
    }
    els.statStrip.innerHTML = html;
  }

  // ---- summary band ----
  function renderSummary(data) {
    const a = state.agg;
    els.summary.innerHTML = STAGE_ORDER.map((st) => {
      const v = data[st], n = v.length, det = v.filter((x) => x > 0).length;
      const m = mean(v), sem = stderr(v), med = median(v);
      const f = (x) => x == null ? "—" : fmt(x);
      return `<article><div class="s-title"><i style="background:${a.stage_colors[st]}"></i>` +
        `<h4>${a.stage_labels[st]}</h4></div><dl>` +
        `<div><dt>Samples</dt><dd>${n}</dd></div>` +
        `<div><dt>Detected</dt><dd>${det}</dd></div>` +
        `<div><dt>Mean ± SEM</dt><dd>${f(m)} ± ${f(sem)}</dd></div>` +
        `<div><dt>Median</dt><dd>${f(med)}</dd></div></dl></article>`;
    }).join("");
  }

  // ---- the SVG chart (bar / box / violin), Harry's structure ----
  function drawChart(data, detail) {
    const a = state.agg;
    const Wv = 1120, Hv = 560, L = 92, R = 1092, T = 34, B = 490;
    const tf = state.log ? ((x) => Math.log10(x + 1)) : ((x) => x);
    const allV = STAGE_ORDER.flatMap((s) => data[s]);
    const maxRaw = allV.length ? Math.max(...allV) : 1;
    const ymaxT = tf(maxRaw) * 1.02 || 1;
    const yToPx = (t) => B - (t / (ymaxT || 1)) * (B - T);
    const bandW = (R - L) / 3, barW = 132;
    const cx = (i) => L + bandW * i + bandW / 2;

    let svg = `<svg class="expr-chart" viewBox="0 0 ${Wv} ${Hv}" role="img" ` +
      `aria-label="${state.mode} plot of ${state.gene} across stages"><rect x="${L}" y="${T}" ` +
      `width="${R - L}" height="${B - T}" class="plot-bg"/>`;

    // y grid + ticks
    const ticks = axisTicks(maxRaw, tf, ymaxT);
    for (const tk of ticks) {
      const y = yToPx(tf(tk));
      svg += `<line x1="${L}" x2="${R}" y1="${y}" y2="${y}" class="grid-line"/>` +
        `<text x="${L - 12}" y="${y + 4}" text-anchor="end" class="axis-tick">${fmtTick(tk)}</text>`;
    }

    STAGE_ORDER.forEach((st, i) => {
      const v = data[st], color = a.stage_colors[st], x = cx(i);
      if (!v.length) return;
      if (state.mode === "bar") svg += barGlyph(v, x, barW, color, tf, yToPx);
      else if (state.mode === "box") svg += boxGlyph(v, x, barW, color, tf, yToPx);
      else svg += violinGlyph(v, x, barW, color, tf, yToPx, ymaxT);
      svg += dotsGlyph(detail[st], x, barW, color, tf, yToPx);
    });

    // axes + stage labels
    svg += `<line x1="${L}" x2="${R}" y1="${B}" y2="${B}" class="axis-line"/>` +
      `<line x1="${L}" x2="${L}" y1="${T}" y2="${B}" class="axis-line"/>`;
    STAGE_ORDER.forEach((st, i) => {
      const x = cx(i);
      svg += `<rect x="${x - 8}" y="512" width="16" height="5" rx="2" fill="${a.stage_colors[st]}"/>` +
        `<text x="${x}" y="540" text-anchor="middle" class="stage-label" fill="${a.stage_colors[st]}">${a.stage_labels[st]}</text>`;
    });
    const ytitle = NORM_UNIT[state.norm] + (state.log ? " (log10(value + 1))" : "");
    svg += `<text x="26" y="${(T + B) / 2}" transform="rotate(-90 26 ${(T + B) / 2})" text-anchor="middle" class="axis-title">${ytitle}</text>`;
    svg += `</svg>`;
    els.chart.innerHTML = svg;
    // dot tooltips → title elements already inline; clicking a dot opens its 3-D
    els.chart.querySelectorAll(".sample-dot").forEach((d) =>
      d.addEventListener("click", () => { const i = +d.dataset.si; if (i >= 0) openFocus(i); }));
  }

  function barGlyph(v, x, w, color, tf, yToPx) {
    const m = mean(v), sem = stderr(v);
    const yTop = yToPx(tf(m)), y0 = yToPx(0);
    let s = `<rect x="${x - w / 2}" y="${yTop}" width="${w}" height="${Math.max(0, y0 - yTop)}" rx="2" ` +
      `fill="${color}" fill-opacity="0.7"><title>Mean ${fmt(m)}; SEM ${fmt(sem)}</title></rect>`;
    if (sem > 0) {
      const yhi = yToPx(tf(m + sem)), ylo = yToPx(tf(Math.max(0, m - sem)));
      s += `<line x1="${x}" x2="${x}" y1="${yhi}" y2="${ylo}" class="err-bar"/>` +
        `<line x1="${x - 22}" x2="${x + 22}" y1="${yhi}" y2="${yhi}" class="err-bar"/>` +
        `<line x1="${x - 22}" x2="${x + 22}" y1="${ylo}" y2="${ylo}" class="err-bar"/>`;
    }
    return s;
  }

  function boxGlyph(v, x, w, color, tf, yToPx) {
    const q = quantiles(v);
    const bw = w * 0.7, x0 = x - bw / 2, x1 = x + bw / 2;
    const yq1 = yToPx(tf(q.q1)), yq3 = yToPx(tf(q.q3)), ymed = yToPx(tf(q.med));
    const ywlo = yToPx(tf(q.wlo)), ywhi = yToPx(tf(q.whi));
    return `<rect x="${x0}" y="${yq3}" width="${bw}" height="${Math.max(1, yq1 - yq3)}" ` +
      `fill="${color}" fill-opacity="0.55" class="box-line"/>` +
      `<line x1="${x0}" x2="${x1}" y1="${ymed}" y2="${ymed}" class="box-line"/>` +
      `<line x1="${x}" x2="${x}" y1="${yq3}" y2="${ywhi}" class="box-line"/>` +
      `<line x1="${x}" x2="${x}" y1="${yq1}" y2="${ywlo}" class="box-line"/>` +
      `<line x1="${x - bw / 4}" x2="${x + bw / 4}" y1="${ywhi}" y2="${ywhi}" class="box-line"/>` +
      `<line x1="${x - bw / 4}" x2="${x + bw / 4}" y1="${ywlo}" y2="${ywlo}" class="box-line"/>`;
  }

  function violinGlyph(v, x, w, color, tf, yToPx, ymaxT) {
    if (v.length < 2) return boxGlyph(v, x, w, color, tf, yToPx);
    const tv = v.map(tf);
    const lo = 0, hi = ymaxT;
    const N = 48, bw = silverman(tv);
    const half = w * 0.46;
    const dens = [];
    let dmax = 0;
    for (let i = 0; i <= N; i++) {
      const y = lo + (hi - lo) * i / N;
      let d = 0; for (const t of tv) d += gauss((y - t) / bw);
      d /= (tv.length * bw); dens.push([y, d]); dmax = Math.max(dmax, d);
    }
    if (dmax === 0) return "";
    const pts = dens.map(([y, d]) => [x - half * d / dmax, yToPx(y)])
      .concat(dens.slice().reverse().map(([y, d]) => [x + half * d / dmax, yToPx(y)]));
    const path = pts.map((p, i) => `${i ? "L" : "M"}${p[0].toFixed(1)} ${p[1].toFixed(1)}`).join(" ") + "Z";
    const q = quantiles(v);
    return `<path d="${path}" fill="${color}" fill-opacity="0.35" stroke="${color}" stroke-width="1.2"/>` +
      `<line x1="${x}" x2="${x}" y1="${yToPx(tf(q.q1))}" y2="${yToPx(tf(q.q3))}" class="box-line"/>` +
      `<circle cx="${x}" cy="${yToPx(tf(q.med))}" r="3" fill="#1f2937"/>`;
  }

  function dotsGlyph(det, x, w, color, tf, yToPx) {
    if (!det.length) return "";
    let rng = 987 + det.length * 13;
    const rand = () => { rng = (rng * 1103515245 + 12345) & 0x7fffffff; return rng / 0x7fffffff; };
    const jw = w * 0.42;
    return det.map((d) => {
      const cxx = x + (rand() - 0.5) * jw, cyy = yToPx(tf(d.value));
      return `<circle cx="${cxx.toFixed(1)}" cy="${cyy.toFixed(1)}" r="5" fill="${color}" ` +
        `stroke="#fff" stroke-width="1.5" class="sample-dot" data-si="${sampleIdx(d.id)}" tabindex="0">` +
        `<title>${d.id}: ${fmt(d.raw)} transcripts${state.norm !== "raw" ? ` (${fmt(d.value)} ${state.norm})` : ""}</title></circle>`;
    }).join("");
  }

  // ---------- CSV ----------
  function exportCsv() {
    if (!state._csv) return;
    const rows = [["sample", "stage", "raw_count", `value_${state.norm}`]];
    for (const st of STAGE_ORDER)
      for (const d of state._csv[st]) rows.push([d.id, state.agg.stage_labels[st], d.raw, d.value]);
    const blob = new Blob([rows.map((r) => r.join(",")).join("\n")], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = `${state.gene}_stage_expression.csv`;
    document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(a.href);
  }

  // ---------- helpers ----------
  const sampleIdxCache = {};
  function sampleIdx(id) {
    if (sampleIdxCache[id] == null) sampleIdxCache[id] = S().findIndex((s) => s.id === id);
    return sampleIdxCache[id];
  }
  function el(tag, cls) { const e = document.createElement(tag); if (cls) e.className = cls; return e; }
  function hexA(hex, a) {
    const n = parseInt(hex.slice(1), 16);
    return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
  }
  const mean = (v) => v.length ? v.reduce((s, x) => s + x, 0) / v.length : null;
  function variance(v) { if (v.length < 2) return 0; const m = mean(v); return v.reduce((s, x) => s + (x - m) ** 2, 0) / (v.length - 1); }
  const stderr = (v) => v.length < 2 ? 0 : Math.sqrt(variance(v) / v.length);
  function median(v) { return quantile(v, 0.5); }
  function quantile(v, p) {
    if (!v.length) return null;
    const s = [...v].sort((a, b) => a - b), idx = (s.length - 1) * p, lo = Math.floor(idx), hi = Math.ceil(idx);
    return s[lo] + (s[hi] - s[lo]) * (idx - lo);
  }
  function quantiles(v) {
    const q1 = quantile(v, 0.25), q3 = quantile(v, 0.75), med = quantile(v, 0.5);
    const iqr = q3 - q1, s = [...v].sort((a, b) => a - b);
    const wlo = Math.max(s[0], q1 - 1.5 * iqr), whi = Math.min(s[s.length - 1], q3 + 1.5 * iqr);
    return { q1, q3, med, wlo, whi };
  }
  const gauss = (z) => Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);
  function silverman(v) {
    const n = v.length, sd = Math.sqrt(variance(v));
    const q1 = quantile(v, 0.25), q3 = quantile(v, 0.75), iqr = (q3 - q1) / 1.349;
    const sigma = Math.min(sd || iqr || 1, iqr || sd || 1) || 1;
    return Math.max(1e-3, 0.9 * sigma * Math.pow(n, -0.2));
  }

  // Welch's t-test, two-sided p via the t-distribution (regularized incomplete beta).
  function welchP(a, b) {
    if (a.length < 2 || b.length < 2) return null;
    const ma = mean(a), mb = mean(b), va = variance(a), vb = variance(b);
    const sa = va / a.length, sb = vb / b.length;
    const se = Math.sqrt(sa + sb); if (se === 0) return 1;
    const t = (ma - mb) / se;
    const df = (sa + sb) ** 2 / (sa ** 2 / (a.length - 1) + sb ** 2 / (b.length - 1));
    return studentTwoSided(t, df);
  }
  function studentTwoSided(t, df) {
    const x = df / (df + t * t);
    return clamp01(betai(df / 2, 0.5, x));
  }
  const clamp01 = (x) => Math.max(0, Math.min(1, x));
  function betai(a, b, x) {
    if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(gammaln(a + b) - gammaln(a) - gammaln(b) + a * Math.log(x) + b * Math.log(1 - x));
    if (x < (a + 1) / (a + b + 2)) return bt * betacf(a, b, x) / a;
    return 1 - bt * betacf(b, a, 1 - x) / b;
  }
  function betacf(a, b, x) {
    const MAXIT = 200, EPS = 3e-12, FPMIN = 1e-300;
    let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d;
      const del = d * c; h *= del; if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }
  function gammaln(x) {
    const cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
      -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += cof[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }

  // formatting
  function fmt(x) {
    if (x == null) return "—";
    if (x === 0) return "0";
    const abs = Math.abs(x);
    if (abs >= 1000) return (x / 1000).toFixed(x >= 10000 ? 0 : 1).replace(/\.0$/, "") + "k";
    if (abs >= 100) return x.toFixed(0);
    if (abs >= 10) return x.toFixed(1);
    if (abs >= 1) return x.toFixed(2);
    return x.toPrecision(2);
  }
  function fmtTick(x) {
    if (x === 0) return "0";
    if (x >= 1000) return (x / 1000).toFixed(x >= 10000 ? 0 : 1).replace(/\.0$/, "") + "k";
    if (x >= 10) return x.toFixed(x >= 100 ? 0 : 1);
    return x.toFixed(2);
  }
  function axisTicks(maxRaw, tf, ymaxT) {
    if (state.log) {
      // even steps in log space, labelled with the original value (Harry's style)
      const n = 5, out = [];
      for (let i = 0; i <= n; i++) { const t = ymaxT * i / n; out.push(Math.pow(10, t) - 1); }
      return out.map((x) => x < 0 ? 0 : x);
    }
    const n = 5, out = [];
    for (let i = 0; i <= n; i++) out.push(maxRaw * i / n);
    return out;
  }
})();
