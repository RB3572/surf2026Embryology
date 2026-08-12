/* Embryo Size.
 *
 * How big is each embryo, and what shape is the space it occupies. Two measurements, both from
 * the segmentation meshes and both defined identically to figure 10.1 so the site and the figure
 * cannot disagree:
 *
 *   MEAN CORTEX RADIUS   for EVERY embryo at every stage — distance from the centre of mass out
 *                        to the cortex, averaged over ~650 directions. Drawn as a sphere.
 *   THE BOX              for 2-cell embryos — the oriented box that just contains both
 *                        blastomeres, along the axis joining their centres. Drawn as a wireframe
 *                        with its length and height labelled.
 *
 * Polar body and nuclei are excluded throughout: they are separate segments and never enter the
 * cytoplasm or blastomere labels, so the exclusion is structural rather than a filter.
 *
 * Every embryo is plotted. Shape irregularity is reported but never used to drop one.
 *
 * Data: data/size.json.gz (build_size.py).
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const C_SPHERE = "#b45309", C_BOX = "#0d9488", C_LABEL = "#0f172a";
  // The figure 10.1–10.4 palette, used verbatim so a panel here and the exported figure read as
  // the same picture: green zygote, blue early 2-cell, orange late 2-cell.
  const STAGE_C = { zygote: "#2f8f6b", early2cell: "#3f76b5", late2cell: "#e5734f" };
  const STAGE_LABEL = { zygote: "zygote", early2cell: "early 2-cell", late2cell: "late 2-cell" };
  // and the figure's line convention: solid = height / radius, dashed = length / diameter
  const LS_H = "solid", LS_L = "dash";

  const state = {
    data: null, byId: {}, stage: "zygote", sort: "radius_um", dir: "desc",
    currentId: null, scene: null, tab: "dist",
    shape: true, labels: true, mesh: true,
  };

  const meta = () => state.data.meta;
  const cur = () => state.byId[state.currentId];
  const isTwo = (e) => e && e.stage !== "zygote";

  // µm -> the PLOT space bodyTraces/sceneLayout use (x,y in px; z scaled)
  const toPlot = (p, zs) => [p[0] / XY, p[1] / XY, p[2] * zs];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const mul = (a, k) => [a[0] * k, a[1] * k, a[2] * k];

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await V.loadGz("data/size.json.gz"); }
    catch (err) {
      $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_size.py</code>.</div></div>`;
      return;
    }
    const m = meta();
    state.data.embryos.forEach((e) => (state.byId[e.id] = e));
    $("#embryo-count").textContent =
      `${m.n} embryos · ${m.n_zygote} zygote · ${m.n_e2c} early + ${m.n_l2c} late 2-cell · ` +
      `median radius ${m.medians.zygote_radius_um} µm (zygote)`;
    $("#sw-shape").style.background = C_SPHERE;
    fillSorts();
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    openDrawer(true);
    renderList();
    const first = ranked()[0];
    if (first) selectEmbryo(first.id);
  })();

  const sortsFor = () => meta().sorts[state.stage === "zygote" ? "zygote" : "twocell"];

  function fillSorts(reset) {
    const opts = sortsFor();
    $("#sort-sel").innerHTML = opts.map((o) => `<option value="${o.key}">${o.label}</option>`).join("");
    if (reset || !opts.some((o) => o.key === state.sort)) state.sort = opts[0].key;
    $("#sort-sel").value = state.sort;
  }

  /** The embryos in the current stage group, ordered. Nothing is ever filtered out. */
  function ranked() {
    const want = state.stage;
    const rows = state.data.embryos.filter((e) => {
      if (want === "twocell") return e.stage !== "zygote";
      return e.stage === want;
    }).filter((e) => e[state.sort] != null)
      .map((e) => ({ id: e.id, e, v: e[state.sort] }));
    rows.sort((a, b) => (state.dir === "asc" ? a.v - b.v : b.v - a.v));
    return rows;
  }

  function renderList() {
    const rows = ranked(), el = $("#sz-list");
    const lab = (sortsFor().find((o) => o.key === state.sort) || {}).label || state.sort;
    if (!rows.length) { el.innerHTML = `<div class="sz-empty">Nothing in this group.</div>`; return; }
    el.innerHTML =
      `<div class="sz-head-row"><span></span><span>embryo</span><span>µm</span><span></span></div>` +
      rows.map((r, i) => {
        const st = r.e.stage;
        return `<div class="sz-row${r.id === state.currentId ? " current" : ""}" data-id="${r.id}"
            title="${r.id} · ${STAGE_LABEL[st]}${r.e.noted ? " · irregular cortex" : ""}">
          <span class="n">${i + 1}</span>
          <span class="e">${V.embryoLabel ? V.embryoLabel(r.id) : r.id}</span>
          <span class="d">${r.v.toFixed(1)}</span>
          <span class="g" style="background:${STAGE_C[st]}" title="${STAGE_LABEL[st]}"></span></div>`;
      }).join("");
    el.querySelectorAll(".sz-row").forEach((row) =>
      row.addEventListener("click", () => selectEmbryo(row.dataset.id)));
    $("#sz-dist-sub").textContent = `· sorted by ${lab.toLowerCase()}`;
  }

  // ───────── 3-D ─────────
  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    const e = state.byId[id];
    $("#loading").hidden = false; $("#loading-text").textContent = `Loading ${id}…`;
    try { state.scene = await V.loadGz(`data/segments/${e.scene}`); }
    catch (err) { state.scene = null; }
    $("#loading").hidden = true;
    $("#controls").hidden = false; $("#placeholder").hidden = true;
    $("#drawer-emb").textContent = V.embryoLabel ? V.embryoLabel(id) : id;
    render3D(); renderReadout(); renderList(); renderActive();
  }

  /** A lat/long wireframe sphere as ONE trace, using null breaks between rings. */
  function sphereTrace(centre, r, zs, color) {
    const X = [], Y = [], Z = [];
    const push = (p) => { const q = toPlot(p, zs); X.push(q[0]); Y.push(q[1]); Z.push(q[2]); };
    const brk = () => { X.push(null); Y.push(null); Z.push(null); };
    const NLAT = 7, NLON = 12, NSEG = 48;
    for (let i = 1; i <= NLAT; i++) {                 // latitude rings
      const phi = Math.PI * i / (NLAT + 1), rr = r * Math.sin(phi), zz = r * Math.cos(phi);
      for (let k = 0; k <= NSEG; k++) {
        const t = 2 * Math.PI * k / NSEG;
        push([centre[0] + rr * Math.cos(t), centre[1] + rr * Math.sin(t), centre[2] + zz]);
      }
      brk();
    }
    for (let j = 0; j < NLON; j++) {                  // longitude rings
      const th = Math.PI * j / NLON;
      for (let k = 0; k <= NSEG; k++) {
        const t = 2 * Math.PI * k / NSEG;
        push([centre[0] + r * Math.sin(t) * Math.cos(th),
              centre[1] + r * Math.sin(t) * Math.sin(th),
              centre[2] + r * Math.cos(t)]);
      }
      brk();
    }
    return { type: "scatter3d", mode: "lines", name: "mean-radius sphere", x: X, y: Y, z: Z,
             line: { color, width: 1.6 }, opacity: 0.75, hoverinfo: "skip", legendrank: 10 };
  }

  /** The 12 edges of the oriented box as one trace, plus its corner set. */
  function boxCorners(b) {
    const c = b.centre_um, u = b.u, v = b.v, w = b.w;
    const out = [];
    for (const su of [-1, 1]) for (const sv of [-1, 1]) for (const sw of [-1, 1]) {
      out.push(add(add(add(c, mul(u, su * b.half_u)), mul(v, sv * b.half_v)), mul(w, sw * b.half_w)));
    }
    return out;                                   // index bits: u<<2 | v<<1 | w
  }

  function boxTrace(b, zs, color) {
    const P = boxCorners(b);
    const idx = (u, v, w) => (u << 2) | (v << 1) | w;
    const edges = [];
    for (const v of [0, 1]) for (const w of [0, 1]) edges.push([idx(0, v, w), idx(1, v, w)]);
    for (const u of [0, 1]) for (const w of [0, 1]) edges.push([idx(u, 0, w), idx(u, 1, w)]);
    for (const u of [0, 1]) for (const v of [0, 1]) edges.push([idx(u, v, 0), idx(u, v, 1)]);
    const X = [], Y = [], Z = [];
    edges.forEach(([a, b2]) => {
      [P[a], P[b2]].forEach((p) => { const q = toPlot(p, zs); X.push(q[0]); Y.push(q[1]); Z.push(q[2]); });
      X.push(null); Y.push(null); Z.push(null);
    });
    return { type: "scatter3d", mode: "lines", name: "bounding box", x: X, y: Y, z: Z,
             line: { color, width: 3 }, hoverinfo: "skip", legendrank: 10 };
  }

  function render3D() {
    const e = cur(), sc = state.scene;
    if (!e || !sc) return;
    const zs = sc.z_scale || 7;
    const traces = state.mesh ? V.bodyTraces(sc) : [];

    if (state.shape) {
      if (isTwo(e)) {
        traces.push(boxTrace(e.box, zs, C_BOX));
      } else {
        traces.push(sphereTrace(e.centre_um, e.radius_um, zs, C_SPHERE));
      }
    }
    if (state.labels) {
      const T = { type: "scatter3d", mode: "text", x: [], y: [], z: [], text: [],
                  textfont: { size: 13, color: C_LABEL }, hoverinfo: "skip", showlegend: false };
      const put = (p, t) => { const q = toPlot(p, zs); T.x.push(q[0]); T.y.push(q[1]); T.z.push(q[2]); T.text.push(t); };
      if (isTwo(e)) {
        const b = e.box;
        // length label at the middle of a bottom-front length edge; height at a vertical edge
        const cL = add(add(b.centre_um, mul(b.v, -b.half_v)), mul(b.w, -b.half_w));
        put(add(cL, mul(b.u, b.half_u * 0.15)), `length ${e.length_um.toFixed(1)} µm`);
        const cH = add(add(b.centre_um, mul(b.u, b.half_u)), mul(b.w, -b.half_w));
        put(cH, `height ${e.height_um.toFixed(1)} µm`);
      } else {
        put(add(e.centre_um, [0, 0, e.radius_um]), `r̄ ${e.radius_um.toFixed(1)} µm`);
      }
      traces.push(T);
    }
    // The scene extents are sized to the MESH. A bounding box reaches past the blastomeres by
    // construction, and a mean-radius sphere can too, so the drawn shape gets clipped unless the
    // view is widened to contain it.
    const ex = JSON.parse(JSON.stringify(sc.extents));
    const grow = (pts) => pts.forEach((p) => {
      const q = toPlot(p, zs);
      ["x", "y", "z"].forEach((k, i) => {
        ex[k][0] = Math.min(ex[k][0], q[i]);
        ex[k][1] = Math.max(ex[k][1], q[i]);
      });
    });
    if (state.shape) {
      if (isTwo(e)) grow(boxCorners(e.box));
      else grow([[-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1]]
                  .map((d) => add(e.centre_um, mul(d, e.radius_um))));
    }
    const host = $("#plot-host");
    host.innerHTML = "";
    Plotly.newPlot(host, traces, V.sceneLayout(ex, sc.id), V.plotConfig);
  }

  function renderReadout() {
    const e = cur(); if (!e) return;
    const L = [`<div class="sz-r-head">${V.embryoLabel ? V.embryoLabel(e.id) : e.id}</div>`];
    const line = (k, v, cls = "") => `<div class="sz-r-line ${cls}"><span class="k">${k}</span><span class="v">${v}</span></div>`;
    L.push(line("stage", STAGE_LABEL[e.stage]));
    L.push(line("mean cortex radius", `${e.radius_um.toFixed(1)} µm`, "is-key"));
    L.push(line("radius range", `${e.radius_min_um.toFixed(1)} – ${e.radius_max_um.toFixed(1)} µm`));
    L.push(line("radial CV", `${e.radial_cv_pct.toFixed(1)} %`));
    if (isTwo(e)) {
      L.push(`<div class="sz-r-sep">bounding box</div>`);
      L.push(line("length", `${e.length_um.toFixed(1)} µm`, "is-key"));
      L.push(line("height", `${e.height_um.toFixed(1)} µm`, "is-key"));
      L.push(line("blastomere separation", `${e.sep_um.toFixed(1)} µm`));
    }
    const rows = ranked(), pos = rows.findIndex((r) => r.id === e.id);
    if (pos >= 0) L.push(line("rank in this group", `${pos + 1} of ${rows.length}`));
    if (e.noted) {
      L.push(`<div class="sz-warn">The cortex dips to <b>${e.radius_min_um.toFixed(1)} µm</b> in some
        direction against a ${e.radius_um.toFixed(1)} µm mean. No cell is that shape, so this is very
        likely a hole in the label surface — the radius is still plotted, but treat it with care.</div>`);
    }
    $("#sz-readout").innerHTML = L.join("");
  }

  // ───────── bottom drawer ─────────
  const baseLayout = (xt, xr) => ({
    margin: { l: 52, r: 10, t: 26, b: 40 }, showlegend: false,
    xaxis: { title: { text: xt, font: { size: 10 } }, range: xr, zeroline: false,
      gridcolor: "#eef1f5", tickfont: { size: 9 } },
    yaxis: { title: { text: "density (scaled to peak)", font: { size: 10 } }, range: [-0.075, 1.08],
      zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
    paper_bgcolor: "transparent", plot_bgcolor: "transparent",
  });
  // the figure's axis windows, so a curve occupies the same span of page here as in the export
  const XLIM_R = [10, 80], XLIM_BOX = [25, 240], XLIM_ALL = [10, 240];

  // ---- density curves, matching figure 10.1-10.4 exactly ----------------------------------
  // A Gaussian KDE with Scott's rule — scipy's gaussian_kde default, so a curve drawn here and
  // one drawn by make_index10.py are the same curve rather than two smoothings that merely
  // resemble each other. h = sd * n^(-1/5) in one dimension.
  function kde(v, xs) {
    const n = v.length, mean = v.reduce((s, x) => s + x, 0) / n;
    const sd = Math.sqrt(v.reduce((s, x) => s + (x - mean) ** 2, 0) / (n - 1));
    const h = sd * Math.pow(n, -1 / 5) || 1e-6;
    const k = 1 / (n * h * Math.sqrt(2 * Math.PI));
    return xs.map((x) => k * v.reduce((s, xi) => s + Math.exp(-0.5 * ((x - xi) / h) ** 2), 0));
  }
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));

  /** One labelled distribution: the peak-normalised curve, its mean ± 1 SD band, a mean rule, and
   *  a rug of the actual values so the smooth curve is never mistaken for the data. */
  function densityTraces(vRaw, color, label, xlim,
                         { ls = LS_H, fill = 0.28, rug = true, mean = true } = {}) {
    const v = vRaw.filter((x) => Number.isFinite(x));
    if (v.length < 2) return { traces: [], mu: NaN, sd: NaN };
    const xs = linspace(xlim[0], xlim[1], 600);
    const y0 = kde(v, xs);
    const peak = Math.max(...y0) || 1;
    const y = y0.map((q) => q / peak);                 // peak = 1, so widths stay comparable
    const n = v.length, mu = v.reduce((s, x) => s + x, 0) / n;
    const sd = Math.sqrt(v.reduce((s, x) => s + (x - mu) ** 2, 0) / (n - 1));
    const rgb = hexRgb(color);
    const band = xs.map((x, i) => (x >= mu - sd && x <= mu + sd ? y[i] : null));
    const traces = [
      { type: "scatter", mode: "lines", x: xs, y: band, fill: "tozeroy",
        fillcolor: `rgba(${rgb},${fill})`, line: { width: 0 }, hoverinfo: "skip",
        showlegend: false, connectgaps: false },
      { type: "scatter", mode: "lines", x: xs, y, name: `${label}  (n = ${n})`,
        line: { color, width: 2, dash: ls },
        hovertemplate: `${label}<br>%{x:.1f} µm<extra></extra>` },
    ];
    if (mean) {
      traces.push({ type: "scatter", mode: "lines", x: [mu, mu], y: [0, kde(v, [mu])[0] / peak],
        line: { color, width: 1.2, dash: "dash" }, hoverinfo: "skip", showlegend: false });
    }
    if (rug) {
      traces.push({ type: "scatter", mode: "markers", x: v, y: v.map(() => -0.035),
        marker: { symbol: "line-ns-open", size: 7, color, opacity: 0.55,
                  line: { color, width: 1 } },
        hovertemplate: `%{x:.1f} µm<extra></extra>`, showlegend: false });
    }
    return { traces, mu, sd };
  }

  const hexRgb = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16)).join(",");

  const CFG = { displaylogo: false, responsive: true,
                modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                toImageButtonOptions: { format: "png", scale: 4 } };

  function plotInto(el, traces, layout) {
    el.innerHTML = ""; el.classList.remove("js-plotly-plot");
    Plotly.newPlot(el, traces, layout, CFG);
  }

  const vals = (stage, key) => state.data.embryos.filter((e) => e.stage === stage && e[key] != null)
                                                 .map((e) => e[key]);

  function marker(x, layout, color) {
    layout.shapes = [{ type: "line", x0: x, x1: x, yref: "paper", y0: 0, y1: 1,
                       line: { color: color || "#111827", width: 2, dash: "dot" } }];
    return layout;
  }

  // inside the axes at upper right, as in the exported figures — and it keeps the panel title,
  // which sits above the plot, from colliding with it
  const legendIn = { x: 0.99, y: 0.99, xanchor: "right", yanchor: "top", font: { size: 9 },
                     bgcolor: "rgba(255,255,255,0.72)", borderwidth: 0 };

  function renderDist() {
    const host = $("#sz-rad");
    if (!host || !host.offsetParent) return;
    const e = cur();

    // panel 1 — mean cortex radius, EVERY embryo, all three stages (figure 10.1)
    const t1 = [];
    ["zygote", "early2cell", "late2cell"].forEach((st) => {
      t1.push(...densityTraces(vals(st, "radius_um"), STAGE_C[st], STAGE_LABEL[st], XLIM_R,
                               { fill: 0.22 }).traces);
    });
    const l1 = baseLayout("mean radius, COM → cortex (µm)", XLIM_R);
    l1.showlegend = true; l1.legend = legendIn; l1.margin.t = 22;
    if (e) marker(e.radius_um, l1);
    plotInto(host, t1, l1);

    // panels 2 and 3 — the 2-cell box, per stage (figures 10.2, 10.3)
    // solid = height, dashed = length, both in that stage's colour — the figure's convention
    [["#sz-e2c", "early2cell", "Early 2-cell"], ["#sz-l2c", "late2cell", "Late 2-cell"]].forEach(
      ([sel, st, ttl]) => {
        const el = $(sel); if (!el) return;
        const col = STAGE_C[st];
        const H = vals(st, "height_um"), L = vals(st, "length_um");
        const t = [
          ...densityTraces(H, col, "box height", XLIM_BOX, { ls: LS_H, fill: 0.28 }).traces,
          ...densityTraces(L, col, "box length", XLIM_BOX, { ls: LS_L, fill: 0.14 }).traces,
        ];
        const lay = baseLayout("box dimension (µm)", XLIM_BOX);
        lay.showlegend = true; lay.legend = legendIn; lay.margin.t = 22;
        lay.title = { text: `${ttl} · n = ${L.length}`, x: 0, xanchor: "left",
                      font: { size: 10, color: "#64748b" } };
        if (e && e.stage === st) {
          lay.shapes = [
            { type: "line", x0: e.height_um, x1: e.height_um, yref: "paper", y0: 0, y1: 1,
              line: { color: C_LABEL, width: 1.5, dash: "dot" } },
            { type: "line", x0: e.length_um, x1: e.length_um, yref: "paper", y0: 0, y1: 1,
              line: { color: C_LABEL, width: 1.5, dash: "dot" } },
          ];
        }
        plotInto(el, t, lay);
      });
  }

  /** Figure 10.4: every distribution on one micron scale, each curve peak-normalised. */
  function renderAll() {
    const host = $("#sz-all");
    if (!host || !host.offsetParent) return;
    const e = cur();
    const zr = vals("zygote", "radius_um");
    const series = [
      ["zygote radius", zr, STAGE_C.zygote, LS_H],
      ["zygote diameter (2r)", zr.map((v) => 2 * v), STAGE_C.zygote, LS_L],
      ["early 2-cell height", vals("early2cell", "height_um"), STAGE_C.early2cell, LS_H],
      ["early 2-cell length", vals("early2cell", "length_um"), STAGE_C.early2cell, LS_L],
      ["late 2-cell height", vals("late2cell", "height_um"), STAGE_C.late2cell, LS_H],
      ["late 2-cell length", vals("late2cell", "length_um"), STAGE_C.late2cell, LS_L],
    ];
    const traces = [];
    series.forEach(([name, v, color, ls]) => {
      // no rug and no mean rule here, exactly as figure 10.4 draws it: six overlapping rugs would
      // be noise, and a dashed vertical among dashed curves would read as another distribution
      traces.push(...densityTraces(v, color, name, XLIM_ALL,
                                   { ls, fill: 0.16, rug: false, mean: false }).traces);
    });
    const lay = baseLayout("µm", XLIM_ALL);
    lay.showlegend = true; lay.margin.t = 10; lay.margin.l = 56;
    lay.legend = legendIn;
    lay.yaxis.range = [-0.03, 1.08];
    if (e) lay.shapes = [{ type: "line", x0: e.radius_um, x1: e.radius_um, yref: "paper",
                           y0: 0, y1: 1, line: { color: C_LABEL, width: 1.5, dash: "dot" } }];
    plotInto(host, traces, lay);
    $("#sz-all-sub").textContent =
      `· each curve scaled to its own peak · shaded = mean ± 1 SD · solid = height / radius, ` +
      `dashed = length / diameter` + (e ? ` · dotted line = this embryo's radius` : "");
  }

  const RENDER = { dist: renderDist, all: renderAll };
  function renderActive() { safely(state.tab, () => (RENDER[state.tab] || renderDist)()); }

  function safely(what, fn) {
    try { fn(); } catch (err) {
      console.error(`[size] ${what} failed`, err);
      const host = what === "all" ? $("#sz-all") : $("#sz-rad");
      if (host) host.innerHTML = `<div class="sz-empty">This panel hit an error: ${err.message}</div>`;
    }
  }

  // ───────── chrome ─────────
  function openDrawer(open) {
    $("#drawer").dataset.open = open ? "true" : "false";
    $("#drawer-handle").setAttribute("aria-expanded", String(open));
    if (open) setTimeout(renderActive, 30);
  }

  function wire() {
    $("#stage-sel").addEventListener("change", (ev) => {
      state.stage = ev.target.value; fillSorts(true); renderList();
      const first = ranked()[0];
      if (first && !ranked().some((r) => r.id === state.currentId)) selectEmbryo(first.id);
    });
    $("#sort-sel").addEventListener("change", (ev) => { state.sort = ev.target.value; renderList(); renderReadout(); });
    $("#sz-dir").addEventListener("click", (ev) => {
      const b = ev.target.closest("button[data-dir]"); if (!b) return;
      state.dir = b.dataset.dir;
      $("#sz-dir").querySelectorAll("button").forEach((x) => x.classList.toggle("active", x === b));
      renderList(); renderReadout();
    });
    [["tg-shape", "shape"], ["tg-labels", "labels"], ["tg-mesh", "mesh"]].forEach(([id, k]) =>
      $("#" + id).addEventListener("change", (ev) => { state[k] = ev.target.checked; render3D(); }));

    $("#drawer-handle").addEventListener("click", () => openDrawer($("#drawer").dataset.open !== "true"));
    $("#sz-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#sz-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b; x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#sz-panels").querySelectorAll(".xs-panel").forEach((p) => (p.hidden = p.dataset.tab !== state.tab));
      setTimeout(renderActive, 20);
    });
    const rh = $("#rdrawer-handle");
    rh.addEventListener("click", () => {
      const r = $("#rdrawer"), open = r.dataset.open !== "true";
      r.dataset.open = String(open); rh.setAttribute("aria-expanded", String(open));
    });
    V.wireWindow($("#controls"), $("#controls-header"), [...$("#controls").querySelectorAll(".rz")], "size_controls_box");
    window.addEventListener("resize", () => {
      try { Plotly.Plots.resize($("#plot-host")); } catch (_) {}
      renderActive();
    });
    window.addEventListener("vcore:dark", () => render3D());
  }
})();
