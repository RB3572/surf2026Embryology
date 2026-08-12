/* Scheffler Pseudotime Model Training.
 *
 * The calibration cohort, and the fitting, made touchable: watch a real live-imaged zygote's
 * pronuclei migrate, then hold embryos out, refit any of the candidate models in the browser, and
 * see what the model you just built says about our own fixed zygotes.
 *
 * ⚠️ THE CARTOON IS A SCHEMATIC. The published workbook gives DISTANCES ONLY — no angles, no cell
 * outline, no cell radius. The radial distances drawn are real data; where the two pronuclei sit
 * around the circle, and the circle itself, are drawing conventions. The page says so on screen
 * rather than letting the picture imply measurements that were never made.
 *
 * Data: data/scheffler.json.gz (build_scheffler.py). Fitting: pt-models.js.
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore, PT = window.PTModels;
  const C_A = "#b58900", C_B = "#2aa198", C_CELL = "#cbd5e1", C_TRUE = "#111827", C_PRED = "#4f46e5";

  const state = {
    data: null, byId: {}, curId: null, frame: 0, playing: null,
    model: "isotonic_sum", test: new Set(), fits: {}, runs: {}, rtab: "pub", yaxis: "__total__",
  };
  const meta = () => state.data.meta;
  const cur = () => state.byId[state.curId];

  /** Every training frame as a flat row: {emb, t, f:{nearer,farther,sum,diff}}. */
  let ROWS = [];

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await V.loadGz("data/scheffler.json.gz"); }
    catch (e) {
      $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_scheffler.py</code>.</div></div>`;
      return;
    }
    const m = meta();
    state.data.train.forEach((e) => (state.byId[e.id] = e));
    ROWS = [];
    state.data.train.forEach((e) => e.f.forEach((fr) =>
      ROWS.push({ emb: e.id, t: fr[1], f: PT.featsFromPair(fr[2], fr[3]) })));

    $("#embryo-count").textContent =
      `${m.n_train} live-imaged zygotes · ${m.n_frames} frames · ${m.duration_h[0]}–${m.duration_h[1]} h ` +
      `· applied to ${m.n_ours} of our fixed zygotes`;

    V.buildTabs($("#tabs"), state.data.train.map((e) => ({ id: e.id, label: e.label, idx: e.idx })),
      selectEmbryo, (e) => ({ label: e.label, sub: `#${e.idx}`,
        title: `${e.id} · datasheet index ${e.idx}` }));

    buildModelTabs();
    buildChips();
    pickFold();
    buildYAxis();
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    renderRank();
    selectEmbryo(state.data.train[0].id);
  })();

  // ───────── the cartoon ─────────
  function selectEmbryo(id) {
    stopPlay();
    state.curId = id; state.frame = 0;
    V.markActiveTab($("#tabs"), id);
    const e = cur();
    $("#controls").hidden = false; $("#placeholder").hidden = true;
    const sl = $("#t-slider");
    sl.max = String(e.f.length - 1); sl.value = "0";
    renderCartoon(); renderTrack();
  }

  /** Distances are data; the angular placement is not. The two pronuclei are drawn on a fixed
   *  vertical diameter, opposed, because the workbook records no angle between them. */
  function renderCartoon() {
    const e = cur(); if (!e) return;
    const fr = e.f[state.frame];
    const [tReal, t, d1, d2] = fr;
    const near = Math.min(d1, d2), far = Math.max(d1, d2);
    // a schematic cell: big enough to contain the largest distance this embryo ever shows
    const R = Math.max(...e.f.map((q) => Math.max(q[2], q[3]))) * 1.22;

    const ring = (cx, cy, r, n = 120) => {
      const x = [], y = [];
      for (let i = 0; i <= n; i++) {
        const a = 2 * Math.PI * i / n;
        x.push(cx + r * Math.cos(a)); y.push(cy + r * Math.sin(a));
      }
      return { x, y };
    };
    const cell = ring(0, 0, R);
    const traces = [
      // dashed on purpose: a solid outline would read as a measured cell boundary, and no cell
      // radius exists anywhere in the source workbook
      { type: "scatter", mode: "lines", x: cell.x, y: cell.y, name: "cell (schematic)",
        line: { color: C_CELL, width: 2, dash: "dot" }, fill: "toself",
        fillcolor: "rgba(203,213,225,0.15)", hoverinfo: "skip", showlegend: false },
      // the two distances, drawn as spokes from the centre
      { type: "scatter", mode: "lines", x: [0, 0], y: [0, far], line: { color: C_A, width: 2, dash: "dot" },
        hoverinfo: "skip", showlegend: false },
      { type: "scatter", mode: "lines", x: [0, 0], y: [0, -near], line: { color: C_B, width: 2, dash: "dot" },
        hoverinfo: "skip", showlegend: false },
    ];
    const pn = (cy, r, color, label, dist) => {
      const c = ring(0, cy, r);
      traces.push({ type: "scatter", mode: "lines", x: c.x, y: c.y, name: label,
        line: { color, width: 2 }, fill: "toself", fillcolor: color, opacity: 0.55,
        hovertemplate: `${label}<br>${dist.toFixed(1)} µm from centre<extra></extra>` });
    };
    pn(far, R * 0.20, C_A, "farther pronucleus", far);
    pn(-near, R * 0.20, C_B, "nearer pronucleus", near);
    traces.push({ type: "scatter", mode: "markers", x: [0], y: [0], name: "cell centre",
      marker: { symbol: "cross", size: 9, color: "#111827", line: { width: 1 } },
      hoverinfo: "skip", showlegend: false });

    const lim = R * 1.12;
    Plotly.react($("#sf-cartoon"), traces, {
      margin: { l: 8, r: 8, t: 8, b: 8 }, showlegend: false, dragmode: false,
      xaxis: { range: [-lim, lim], visible: false, scaleanchor: "y", scaleratio: 1 },
      yaxis: { range: [-lim, lim], visible: false },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      annotations: [{ x: 0, y: -lim * 0.93, xref: "x", yref: "y", showarrow: false,
        text: `∂ = ${(d1 + d2).toFixed(1)} µm`, font: { size: 13, color: "#111827" } }],
    }, { displaylogo: false, responsive: true, staticPlot: false, displayModeBar: false });

    $("#t-lab").textContent =
      `frame ${state.frame + 1} / ${e.f.length} · t_real ${tReal.toFixed(2)} h of ${e.T} h · t = ${t.toFixed(3)}`;
    const L = [];
    L.push(`<div class="sf-r-head">${e.label} <span class="sf-idx">datasheet #${e.idx}</span></div>`);
    const row = (k, v, cls = "") => `<div class="sf-r-line ${cls}"><span class="k">${k}</span><span class="v">${v}</span></div>`;
    L.push(row("T_duration", `${e.T} h`));
    L.push(row("t_real", `${tReal.toFixed(2)} h`));
    L.push(row("t (normalised)", t.toFixed(3), "is-key"));
    L.push(row("nearer ∂₁", `${near.toFixed(1)} µm`));
    L.push(row("farther ∂₂", `${far.toFixed(1)} µm`));
    L.push(row("∂ = ∂₁ + ∂₂", `${(d1 + d2).toFixed(1)} µm`, "is-key"));
    L.push(`<div class="sf-warn">The circle and the pronuclei's positions around it are a
      <b>drawing convention</b> — the source data records distances from the centre and nothing
      else. Only the two radii are measurements.</div>`);
    $("#sf-readout").innerHTML = L.join("");
  }

  /** The whole trajectory as a line, with the current frame marked. */
  function renderTrack() {
    const e = cur(); if (!e) return;
    const t = e.f.map((q) => q[1]);
    const sum = e.f.map((q) => q[2] + q[3]);
    const near = e.f.map((q) => Math.min(q[2], q[3]));
    const far = e.f.map((q) => Math.max(q[2], q[3]));
    const here = state.frame;
    Plotly.react($("#sf-track"), [
      { type: "scatter", mode: "lines", x: t, y: sum, name: "∂ = ∂₁+∂₂",
        line: { color: C_TRUE, width: 2.4 }, hovertemplate: "t %{x:.3f}<br>∂ %{y:.1f} µm<extra></extra>" },
      { type: "scatter", mode: "lines", x: t, y: far, name: "farther",
        line: { color: C_A, width: 1.4 }, hoverinfo: "skip" },
      { type: "scatter", mode: "lines", x: t, y: near, name: "nearer",
        line: { color: C_B, width: 1.4 }, hoverinfo: "skip" },
      { type: "scatter", mode: "markers", x: [t[here]], y: [sum[here]], name: "this frame",
        marker: { size: 11, color: C_TRUE, line: { color: "#fff", width: 2 } }, hoverinfo: "skip" },
    ], {
      margin: { l: 52, r: 12, t: 24, b: 40 }, showlegend: true,
      legend: { orientation: "h", x: 0, y: 1.04, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      xaxis: { title: { text: "t — normalised true time", font: { size: 10 } }, range: [-0.02, 1.02],
        zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { title: { text: "distance to centre (µm)", font: { size: 10 } }, rangemode: "tozero",
        zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, { displaylogo: false, responsive: true, displayModeBar: false });
  }

  function setFrame(i) {
    const e = cur(); if (!e) return;
    state.frame = Math.max(0, Math.min(e.f.length - 1, i));
    $("#t-slider").value = String(state.frame);
    renderCartoon(); renderTrack();
  }
  function stopPlay() {
    if (state.playing) { clearInterval(state.playing); state.playing = null; }
    const b = $("#t-play"); if (b) b.textContent = "▶ play";
  }
  function togglePlay() {
    if (state.playing) return stopPlay();
    $("#t-play").textContent = "❚❚ pause";
    state.playing = setInterval(() => {
      const e = cur(); if (!e) return stopPlay();
      setFrame(state.frame >= e.f.length - 1 ? 0 : state.frame + 1);
    }, 90);
  }

  // ───────── cohort selection ─────────
  function buildChips() {
    $("#sf-chips").innerHTML = state.data.train.map((e) =>
      `<button type="button" class="sf-chip" data-id="${e.id}" title="${e.id}">${e.label}</button>`).join("");
    $("#sf-chips").addEventListener("click", (ev) => {
      const b = ev.target.closest(".sf-chip"); if (!b) return;
      const id = b.dataset.id;
      if (state.test.has(id)) state.test.delete(id); else state.test.add(id);
      syncChips();
    });
  }
  function syncChips() {
    $("#sf-chips").querySelectorAll(".sf-chip").forEach((b) =>
      b.classList.toggle("on", state.test.has(b.dataset.id)));
    const n = state.test.size, tot = state.data.train.length;
    $("#sf-cohort-count").textContent = `${n} held out · ${tot - n} used to fit`;
    $("#sf-run").disabled = n === 0 || n === tot;
    $("#sf-run-all").disabled = n === 0 || n === tot;
  }
  function pickFold() {
    const a = (state.data.folds || {}).assignment || {};
    state.test = new Set(Object.keys(a).filter((k) => a[k] === 0));
    if (!state.test.size) pickRandom(10);
    syncChips();
  }
  function pickRandom(k) {
    const ids = state.data.train.map((e) => e.id);
    for (let i = ids.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1)); [ids[i], ids[j]] = [ids[j], ids[i]];
    }
    state.test = new Set(ids.slice(0, k));
    syncChips();
  }

  // ───────── model tabs ─────────
  function buildModelTabs() {
    const order = ["isotonic_sum", "hgb_symmetric", "ridge_symmetric", "ridge_symmetric_core",
                   "linear_sum", "isotonic_nearer", "linear_nearer"];
    $("#sf-tabs").innerHTML = order.map((k) =>
      `<button type="button" class="xs-gtab${k === state.model ? " active" : ""}" data-model="${k}"
         role="tab">${PT.SPECS[k].label}</button>`).join("");
    $("#sf-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.model = b.dataset.model;
      $("#sf-tabs").querySelectorAll(".xs-gtab").forEach((x) =>
        x.classList.toggle("active", x === b));
      renderModelHead(); renderFit(); renderTest(); renderOurs();
    });
  }
  function renderModelHead() {
    const s = PT.SPECS[state.model];
    $("#drawer-model").textContent = s.label;
    $("#sf-model-title").textContent = s.label;
    const feats = s.kind === "isotonic" ? [s.feat] : s.feats;
    $("#sf-model-note").innerHTML =
      `<span class="sf-pill">${s.complexity}</span>` +
      feats.map((f) => `<span class="sf-pill feat">${f}</span>`).join("") +
      `<span class="sf-why">${s.note}</span>`;
  }

  function buildYAxis() {
    // gene options, ordered by how many of our zygotes detected them
    const cnt = {};
    state.data.ours.forEach((o) => Object.keys(o.g).forEach((g) => (cnt[g] = (cnt[g] || 0) + 1)));
    const genes = Object.keys(cnt).sort((a, b) => cnt[b] - cnt[a] || a.localeCompare(b));
    $("#sf-yaxis").innerHTML = `<option value="__total__">Total transcripts</option>` +
      genes.map((g) => `<option value="${g}">${g} · in ${cnt[g]} zygotes</option>`).join("");
  }

  // ───────── fitting ─────────
  function runOne(key) {
    const r = PT.holdout(key, ROWS, [...state.test]);
    if (r) { state.fits[key] = r; state.runs[key] = r.metrics; }
    return r;
  }
  function run(all) {
    const keys = all ? PT.KEYS : [state.model];
    $("#sf-status").textContent = "fitting…";
    setTimeout(() => {
      const t0 = Date.now();
      keys.forEach(runOne);
      $("#sf-status").textContent =
        `fitted ${keys.length} model${keys.length > 1 ? "s" : ""} in ${Date.now() - t0} ms`;
      renderFit(); renderTest(); renderOurs();
      if (state.rtab === "mine") renderRank();
    }, 10);
  }

  const CFG = { displaylogo: false, responsive: true,
                modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                toImageButtonOptions: { format: "png", scale: 4 } };
  const empty = (el, msg) => { Plotly.purge(el); el.innerHTML = `<div class="sf-empty">${msg}</div>`; };

  function renderFit() {
    const el = $("#sf-fit"), r = state.fits[state.model];
    if (!r) return empty(el, "Choose a testing cohort and press <b>Fit this model</b>.");
    el.innerHTML = "";
    const s = PT.SPECS[state.model];
    $("#sf-fit-sub").textContent = `· fitted on ${r.n_train_embryos} embryos, ${r.n_train_frames} frames`;

    if (s.kind === "isotonic") {
      const k = r.fitted.model.knots;
      Plotly.newPlot(el, [{
        type: "scatter", mode: "lines", x: k.map((p) => p[0]), y: k.map((p) => p[1]),
        line: { color: C_PRED, width: 2, shape: "hv" }, name: `${k.length} knots`,
        hovertemplate: "∂ %{x:.1f} µm → τ %{y:.3f}<extra></extra>",
      }], {
        margin: { l: 50, r: 12, t: 10, b: 42 }, showlegend: false,
        xaxis: { title: { text: `${s.feat} (µm)`, font: { size: 10 } }, gridcolor: "#eef1f5", tickfont: { size: 9 } },
        yaxis: { title: { text: "τ", font: { size: 10 } }, range: [-0.02, 1.02], gridcolor: "#eef1f5", tickfont: { size: 9 } },
        paper_bgcolor: "transparent", plot_bgcolor: "transparent",
        annotations: [{ x: 0.98, y: 0.96, xref: "paper", yref: "paper", showarrow: false,
          text: `${k.length} knots · monotone decreasing`, font: { size: 10, color: "#64748b" },
          xanchor: "right" }],
      }, CFG);
    } else if (s.kind === "linear") {
      const w = r.fitted.model.w;
      const names = ["intercept"].concat(s.feats);
      Plotly.newPlot(el, [{
        type: "bar", x: names, y: w, marker: { color: C_PRED },
        hovertemplate: "%{x}<br>%{y:.5f}<extra></extra>",
      }], {
        margin: { l: 60, r: 12, t: 10, b: 48 }, showlegend: false,
        xaxis: { tickfont: { size: 9 } },
        yaxis: { title: { text: "coefficient", font: { size: 10 } }, gridcolor: "#eef1f5",
                 zeroline: true, zerolinecolor: "#94a3b8", tickfont: { size: 9 } },
        paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      }, CFG);
    } else {
      // boosted trees have no coefficients; show how often each feature is split on
      const used = {}; s.feats.forEach((f) => (used[f] = 0));
      const walk = (n) => { if (n.leaf !== undefined) return; used[s.feats[n.f]]++; walk(n.L); walk(n.R); };
      r.fitted.model.trees.forEach(walk);
      Plotly.newPlot(el, [{
        type: "bar", x: s.feats, y: s.feats.map((f) => used[f]), marker: { color: C_PRED },
        hovertemplate: "%{x}<br>%{y} splits<extra></extra>",
      }], {
        margin: { l: 56, r: 12, t: 10, b: 48 }, showlegend: false,
        xaxis: { tickfont: { size: 9 } },
        yaxis: { title: { text: "times split on", font: { size: 10 } }, gridcolor: "#eef1f5", tickfont: { size: 9 } },
        paper_bgcolor: "transparent", plot_bgcolor: "transparent",
        annotations: [{ x: 0.98, y: 0.96, xref: "paper", yref: "paper", showarrow: false,
          text: `${r.fitted.model.trees.length} trees · no closed form`, xanchor: "right",
          font: { size: 10, color: "#64748b" } }],
      }, CFG);
    }
  }

  function renderTest() {
    const el = $("#sf-test"), r = state.fits[state.model];
    if (!r) return empty(el, "Nothing fitted yet.");
    el.innerHTML = "";
    const m = r.metrics;
    $("#sf-test-sub").textContent =
      `· ${m.n_embryos} held-out embryos, ${m.n_frames} frames`;
    Plotly.newPlot(el, [
      { type: "scatter", mode: "markers", x: r.test.map((q) => q.t), y: r.pred,
        marker: { size: 5, color: C_PRED, opacity: 0.45, line: { width: 0 } },
        text: r.test.map((q) => state.byId[q.emb].label),
        hovertemplate: "%{text}<br>true t %{x:.3f} → τ %{y:.3f}<extra></extra>", name: "held-out frames" },
      { type: "scatter", mode: "lines", x: [0, 1], y: [0, 1], line: { color: "#111827", width: 1.2, dash: "dash" },
        hoverinfo: "skip", name: "identity" },
    ], {
      margin: { l: 50, r: 12, t: 10, b: 42 }, showlegend: false,
      xaxis: { title: { text: "true t", font: { size: 10 } }, range: [-0.02, 1.02], gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { title: { text: "predicted τ", font: { size: 10 } }, range: [-0.02, 1.02], gridcolor: "#eef1f5", tickfont: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      annotations: [{ x: 0.02, y: 0.98, xref: "paper", yref: "paper", showarrow: false, align: "left",
        xanchor: "left", yanchor: "top", font: { size: 10, color: "#111827" },
        text: `MAE ${m.macro_mae.toFixed(4)} · Spearman ${m.pooled_spearman.toFixed(3)}<br>` +
              `ordering ${(m.pairwise_ordering_accuracy * 100).toFixed(1)}% pooled, ` +
              `${(m.within_embryo_ordering_accuracy * 100).toFixed(1)}% within embryo` }],
    }, CFG);
  }

  function renderOurs() {
    const el = $("#sf-ours"), r = state.fits[state.model];
    if (!r) return empty(el, "Fit a model to see what it infers for our zygotes.");
    el.innerHTML = "";
    const y = state.yaxis;
    const pts = state.data.ours.map((o) => {
      const f = { nearer: o.d1, farther: o.d2, sum: o.sum, diff: o.diff };
      return { o, tau: PT.predict(r.fitted, f), v: y === "__total__" ? o.total_tx : (o.g[y] || 0) };
    }).filter((p) => y === "__total__" || p.o.g[y] !== undefined);
    const lab = y === "__total__" ? "total transcripts" : `${y} transcripts`;
    $("#sf-ours-sub").textContent =
      `· ${pts.length} zygotes · τ ${Math.min(...pts.map((p) => p.tau)).toFixed(2)}–` +
      `${Math.max(...pts.map((p) => p.tau)).toFixed(2)} · ${PT.SPECS[state.model].label}`;
    const tau = pts.map((p) => p.tau), vv = pts.map((p) => p.v);
    const rho = PT.spearman(tau, vv);
    Plotly.newPlot(el, [{
      type: "scatter", mode: "markers", x: tau, y: vv,
      marker: { size: 8, color: C_PRED, opacity: 0.7, line: { color: "#fff", width: 0.6 } },
      text: pts.map((p) => p.o.label),
      hovertemplate: `%{text}<br>τ %{x:.3f}<br>${lab} %{y}<extra></extra>`,
    }], {
      margin: { l: 62, r: 12, t: 10, b: 44 }, showlegend: false,
      xaxis: { title: { text: "τ — inferred pseudotime", font: { size: 10 } }, range: [-0.03, 1.03],
        gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { title: { text: lab, font: { size: 10 } }, rangemode: "tozero", gridcolor: "#eef1f5", tickfont: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      annotations: [{ x: 0.98, y: 0.98, xref: "paper", yref: "paper", showarrow: false, xanchor: "right",
        yanchor: "top", font: { size: 10, color: "#64748b" },
        text: isNaN(rho) ? "" : `Spearman ρ = ${rho.toFixed(2)}` }],
    }, CFG);
  }

  // ───────── right drawer: rankings ─────────
  function renderRank() {
    const el = $("#sf-rank");
    const rows = state.rtab === "pub"
      ? state.data.models.map((m) => ({ key: m.key, label: m.label, dep: m.deployable, m: m.metrics }))
      : Object.keys(state.runs).map((k) => ({ key: k, label: PT.SPECS[k].label, dep: true, m: state.runs[k] }));
    if (!rows.length) {
      el.innerHTML = `<div class="sf-empty-list">Nothing fitted yet. Choose a cohort and press
        <b>Fit every model on this cohort</b> to build your own ranking.</div>`;
      return;
    }
    rows.sort((a, b) => a.m.macro_mae - b.m.macro_mae);
    const head = state.rtab === "pub"
      ? `<div class="sf-rnote">The published comparison: nested five-fold cross-validation grouped by
           embryo, scored offline. Fixed, so it does not move when you fit your own.
           A <span class="sf-leak">leaky</span> row scores well and <b>is not a clock</b> — its
           features are pronuclear volumes normalised to each pronucleus's own future endpoint, so
           they encode elapsed time and cannot be measured in a fixed embryo. It is shown as an
           upper bound, has no tab, and its data are not in this page.</div>`
      : `<div class="sf-rnote">Your runs, on the cohort you defined (${state.test.size} held out).
           Directly comparable to each other, but not to the published table above, which used all
           five folds rather than one split.</div>`;
    el.innerHTML = head +
      `<div class="sf-rank-head"><span></span><span>model</span><span>MAE</span><span>order</span></div>` +
      rows.map((r, i) => `<div class="sf-rank-row${r.key === state.model ? " current" : ""}${r.dep ? "" : " leaky"}"
          data-key="${r.key}" title="${r.label}${r.dep ? "" : " — NON-DEPLOYABLE"}">
        <span class="n">${i + 1}</span>
        <!-- the badge leads, so a clipped label can never hide that this row is not deployable -->
        <span class="e">${r.dep ? "" : '<span class="sf-leak">leaky</span> '}${r.label}</span>
        <span class="d">${r.m.macro_mae.toFixed(4)}</span>
        <span class="o">${(r.m.pairwise_ordering_accuracy * 100).toFixed(0)}%</span></div>`).join("");
    el.querySelectorAll(".sf-rank-row").forEach((row) => row.addEventListener("click", () => {
      const k = row.dataset.key;
      if (!PT.SPECS[k]) return;                        // the leaky model has no refittable tab
      state.model = k;
      $("#sf-tabs").querySelectorAll(".xs-gtab").forEach((x) =>
        x.classList.toggle("active", x.dataset.model === k));
      renderModelHead(); renderFit(); renderTest(); renderOurs(); renderRank();
    }));
  }

  // ───────── chrome ─────────
  function openDrawer(open) {
    $("#drawer").dataset.open = open ? "true" : "false";
    $("#drawer-handle").setAttribute("aria-expanded", String(open));
    if (open) setTimeout(() => { renderFit(); renderTest(); renderOurs(); }, 30);
  }

  /** The drag-to-resize handle both drawers carry, in one place. */
  function dragResize(el, move, start) {
    if (!el) return;
    el.addEventListener("pointerdown", (e) => {
      el._d = Object.assign({ x: e.clientX, y: e.clientY }, start());
      el.setPointerCapture(e.pointerId); e.preventDefault();
      el.classList.add("dragging");
    });
    el.addEventListener("pointermove", (e) => { if (el._d) move(el._d, e); });
    const end = (e) => {
      el._d = null; el.classList.remove("dragging");
      try { el.releasePointerCapture(e.pointerId); } catch (_) {}
    };
    el.addEventListener("pointerup", end);
    el.addEventListener("pointercancel", end);
  }

  function wire() {
    $("#t-slider").addEventListener("input", (e) => { stopPlay(); setFrame(+e.target.value); });
    $("#t-prev").addEventListener("click", () => { stopPlay(); setFrame(state.frame - 1); });
    $("#t-next").addEventListener("click", () => { stopPlay(); setFrame(state.frame + 1); });
    $("#t-play").addEventListener("click", togglePlay);

    $(".sf-cohort-btns").addEventListener("click", (ev) => {
      const b = ev.target.closest("[data-pick]"); if (!b) return;
      const p = b.dataset.pick;
      if (p === "10") pickRandom(10);
      else if (p === "fold") pickFold();
      else if (p === "none") { state.test = new Set(); syncChips(); }
      else if (p === "invert") {
        const all = state.data.train.map((e) => e.id);
        state.test = new Set(all.filter((id) => !state.test.has(id)));
        syncChips();
      }
    });
    $("#sf-run").addEventListener("click", () => run(false));
    $("#sf-run-all").addEventListener("click", () => run(true));
    $("#sf-yaxis").addEventListener("change", (e) => { state.yaxis = e.target.value; renderOurs(); });

    $("#drawer-handle").addEventListener("click", () => openDrawer($("#drawer").dataset.open !== "true"));
    $("#sf-rtabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".sf-rtab"); if (!b) return;
      state.rtab = b.dataset.rtab;
      $("#sf-rtabs").querySelectorAll(".sf-rtab").forEach((x) => x.classList.toggle("active", x === b));
      renderRank();
    });
    const rh = $("#rdrawer-handle");
    rh.addEventListener("click", () => {
      const r = $("#rdrawer"), open = r.dataset.open !== "true";
      r.dataset.open = String(open); rh.setAttribute("aria-expanded", String(open));
    });
    dragResize($("#drawer-resize"), (d, e) => {
      $("#drawer").style.setProperty("--drawer-h",
        Math.max(200, Math.min(window.innerHeight - 120, d.v + (d.y - e.clientY))) + "px");
    }, () => ({ v: $("#drawer-body").getBoundingClientRect().height }));
    dragResize($("#rdrawer-resize"), (d, e) => {
      $("#rdrawer").style.setProperty("--rdrawer-w",
        Math.max(260, Math.min(window.innerWidth - 80, d.v - (e.clientX - d.x))) + "px");
    }, () => ({ v: $("#rdrawer").getBoundingClientRect().width }));
    V.wireWindow($("#controls"), $("#controls-header"), [...$("#controls").querySelectorAll(".rz")],
      "scheffler_controls_box");
    window.addEventListener("resize", () => {
      ["#sf-cartoon", "#sf-track", "#sf-fit", "#sf-test", "#sf-ours"].forEach((s) => {
        try { Plotly.Plots.resize($(s)); } catch (_) {}
      });
    });
    renderModelHead();
  }
})();
