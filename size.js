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
  const STAGE_C = { zygote: "#334155", early2cell: "#c98a3f", late2cell: "#4a7db5" };
  const STAGE_LABEL = { zygote: "zygote", early2cell: "early 2-cell", late2cell: "late 2-cell" };

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
  const baseLayout = (xt) => ({
    margin: { l: 46, r: 10, t: 26, b: 40 }, showlegend: false, bargap: 0.04,
    xaxis: { title: { text: xt, font: { size: 10 } }, zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
    yaxis: { title: { text: "embryos", font: { size: 10 } }, zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
    paper_bgcolor: "transparent", plot_bgcolor: "transparent",
  });

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

  function renderDist() {
    const host = $("#sz-rad");
    if (!host || !host.offsetParent) return;
    const e = cur();

    // panel 1 — mean cortex radius, EVERY embryo, all three stages
    const t1 = ["zygote", "early2cell", "late2cell"].map((st) => ({
      type: "histogram", x: vals(st, "radius_um"), name: STAGE_LABEL[st],
      marker: { color: STAGE_C[st], line: { color: "#fff", width: .5 } },
      opacity: 0.72, nbinsx: 18,
      hovertemplate: `${STAGE_LABEL[st]}<br>%{x:.1f} µm · %{y} embryos<extra></extra>`,
    }));
    const l1 = baseLayout("mean cortex radius (µm)");
    l1.barmode = "overlay"; l1.showlegend = true;
    l1.legend = { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } };
    l1.margin.t = 34;
    if (e) marker(e.radius_um, l1);
    plotInto(host, t1, l1);

    // panels 2 and 3 — the 2-cell box, per stage
    [["#sz-e2c", "early2cell", "Early 2-cell"], ["#sz-l2c", "late2cell", "Late 2-cell"]].forEach(
      ([sel, st, ttl]) => {
        const el = $(sel); if (!el) return;
        const H = vals(st, "height_um"), L = vals(st, "length_um");
        const t = [
          { type: "histogram", x: H, name: "height", marker: { color: "#e5734f", line: { color: "#fff", width: .5 } }, opacity: .75, nbinsx: 16 },
          { type: "histogram", x: L, name: "length", marker: { color: "#3f76b5", line: { color: "#fff", width: .5 } }, opacity: .75, nbinsx: 16 },
        ];
        const lay = baseLayout("box dimension (µm)");
        lay.barmode = "overlay"; lay.showlegend = true; lay.margin.t = 34;
        lay.legend = { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } };
        lay.title = { text: `${ttl} · n = ${L.length}`, x: 0, xanchor: "left", font: { size: 10, color: "#64748b" } };
        if (e && e.stage === st) {
          lay.shapes = [
            { type: "line", x0: e.height_um, x1: e.height_um, yref: "paper", y0: 0, y1: 1, line: { color: "#e5734f", width: 2, dash: "dot" } },
            { type: "line", x0: e.length_um, x1: e.length_um, yref: "paper", y0: 0, y1: 1, line: { color: "#3f76b5", width: 2, dash: "dot" } },
          ];
        }
        plotInto(el, t, lay);
      });
  }

  function renderAll() {
    const host = $("#sz-all");
    if (!host || !host.offsetParent) return;
    const e = cur();
    const series = [
      ["zygote radius", vals("zygote", "radius_um"), STAGE_C.zygote],
      ["zygote diameter (2r)", vals("zygote", "radius_um").map((v) => 2 * v), "#94a3b8"],
      ["early 2-cell radius", vals("early2cell", "radius_um"), STAGE_C.early2cell],
      ["early 2-cell height", vals("early2cell", "height_um"), "#e5734f"],
      ["early 2-cell length", vals("early2cell", "length_um"), "#3f76b5"],
      ["late 2-cell radius", vals("late2cell", "radius_um"), STAGE_C.late2cell],
      ["late 2-cell height", vals("late2cell", "height_um"), "#e5734f"],
      ["late 2-cell length", vals("late2cell", "length_um"), "#3f76b5"],
    ];
    const traces = [];
    series.forEach(([name, v, color], i) => {
      const y = series.length - 1 - i;
      traces.push({ type: "scatter", mode: "markers", name, x: v,
        y: v.map(() => y + (Math.random() - 0.5) * 0.28),
        marker: { size: 6, color, opacity: 0.6, line: { width: 0 } },
        hovertemplate: `${name}<br>%{x:.1f} µm<extra></extra>` });
      const med = v.slice().sort((a, b) => a - b)[Math.floor(v.length / 2)];
      traces.push({ type: "scatter", mode: "lines", x: [med, med], y: [y - 0.34, y + 0.34],
        line: { color: "#111827", width: 2 }, hoverinfo: "skip", showlegend: false });
    });
    const lay = {
      margin: { l: 150, r: 14, t: 10, b: 42 }, showlegend: false,
      xaxis: { title: { text: "µm", font: { size: 10 } }, zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { tickmode: "array", tickvals: series.map((_, i) => series.length - 1 - i),
               ticktext: series.map((s) => s[0]), zeroline: false, gridcolor: "#f6f8fb", tickfont: { size: 10 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    };
    if (e) lay.shapes = [{ type: "line", x0: e.radius_um, x1: e.radius_um, yref: "paper",
                           y0: 0, y1: 1, line: { color: "#111827", width: 1.5, dash: "dot" } }];
    plotInto(host, traces, lay);
    $("#sz-all-sub").textContent = `· ${meta().n} embryos` + (e ? ` · dotted line = this embryo's radius` : "");
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
