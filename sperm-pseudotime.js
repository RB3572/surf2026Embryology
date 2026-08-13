/* Sperm Location vs Pseudotime.
 *
 * For every zygote with a labelled sperm, plots the distance from the sperm to a
 * selectable object — polar body, maternal pronucleus, or paternal pronucleus —
 * against the calibrated pronuclear pseudotime τ. The maternal/paternal identity of
 * each pronucleus is the consensus from the Maternal/Paternal Pronucleus ID project.
 *
 * Main view: the current zygote in 3-D (sperm + the two pronuclei + polar body), with
 * the distance to the selected object drawn. Bottom drawer: that distance vs τ across
 * every sperm-labelled zygote, this one highlighted.
 *
 * Data: data/sperm_pseudotime.json (build_sperm_pseudotime.py) + per-embryo meshes
 * from data/pronuclei/<id>.json.gz. Distances are µm (plot units × 0.15).
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;

  const F = "#d6336c", M = "#2563eb", PB = "#8b6fc4", SP = "#ff2d95", SPLIT = "#a855f7";
  const els = {
    count: $("#embryo-count"), tabs: $("#tabs"), controls: $("#controls"),
    plotHost: $("#plot-host"), placeholder: $("#placeholder"),
    loading: $("#loading"), loadingText: $("#loading-text"),
    objSelect: $("#obj-select"), showAll: $("#show-all"), showBody: $("#show-body"),
    clockSelect: $("#clock-select"),
    readout: $("#sp-readout"), legend: $("#sp-legend"),
    drawer: $("#drawer"), drawerHandle: $("#drawer-handle"), drawerObj: $("#drawer-obj"),
    stat: $("#sp-stat"), scatter: $("#sp-scatter"), note: $("#sp-note"), exportCsv: $("#export-csv"),
  };
  const OBJ_LABEL = { polar: "polar body", maternal: "maternal pronucleus ♀", paternal: "paternal pronucleus ♂" };
  const OBJ_COLOR = { polar: PB, maternal: F, paternal: M };

  const state = { doc: null, byId: {}, obj: "polar", currentId: null, scene: null,
                  drawerOpen: false, sceneCache: {}, clock: "tau",
                  // figure 4.16: the leave-one-out consensus is the default on purpose — see
                  // renderPair. `pair` stays undefined until first fetched, null if it 404s.
                  drawerTab: "clock", pair: undefined, pairKey: null, pairVar: "loo" };
  const activeClock = () => (state.doc.clocks || []).find((c) => c.key === state.clock) || (state.doc.clocks || [])[0];

  // ---------- boot ----------
  (async function init() {
    try { state.doc = await load(); }
    catch (e) { els.placeholder.innerHTML =
      `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
      `<div class="ph-sub">Run <code>python3 build_sperm_pseudotime.py</code>.</div></div>`; return; }
    const d = state.doc;
    d.embryos.forEach((e) => (state.byId[e.id] = e));
    els.count.textContent = `${d.n} sperm-labelled zygotes · pseudotime ${d.model_version}`;
    renderLegend();
    V.buildTabs(els.tabs, d.embryos, selectEmbryo, (e) => ({
      label: e.label, sub: e.date,
      title: `${e.label} · τ ${e.tau}${e.split ? " · split (no M/P)" : ""}`,
    }));
    wire();
    selectEmbryo(d.embryos[0].id);
  })();
  async function load() {
    const r = await fetch("data/sperm_pseudotime.json"); return r.json();
  }

  function wire() {
    els.objSelect.addEventListener("click", (e) => {
      const b = e.target.closest("[data-obj]"); if (!b) return;
      state.obj = b.dataset.obj;
      els.objSelect.querySelectorAll("button").forEach((x) => x.classList.toggle("active", x === b));
      render3D(); renderReadout(); renderScatter();
    });
    els.showAll.addEventListener("change", render3D);
    els.showBody.addEventListener("change", render3D);
    // clock toggle: pronuclei-to-COM τ (default) vs interpronuclei distance → re-plots the scatter x-axis
    if (els.clockSelect) {
      els.clockSelect.addEventListener("click", (e) => {
        const b = e.target.closest("[data-clock]"); if (!b) return;
        state.clock = b.dataset.clock;
        els.clockSelect.querySelectorAll("button").forEach((x) => x.classList.toggle("active", x === b));
        renderReadout(); renderScatter();
      });
    }
    els.exportCsv.addEventListener("click", exportCsv);
    wireDrawer();
  }

  // ---------- select embryo ----------
  async function selectEmbryo(id) {
    state.currentId = id; V.markActiveTab(els.tabs, id);
    const e = state.byId[id];
    els.loading.hidden = false; els.loadingText.textContent = `Loading ${e.label}…`;
    let scene = state.sceneCache[id];
    try {
      if (!scene) { scene = await V.loadGz(`data/pronuclei/${id}.json.gz`); state.sceneCache[id] = scene; }
    } catch (err) { els.loading.hidden = true; return; }
    if (state.currentId !== id) return;
    state.scene = scene;
    els.controls.hidden = false; els.placeholder.hidden = true; els.drawer.hidden = false;
    els.loading.hidden = true;
    render3D(); renderReadout();
    if (!state.drawerOpen) openDrawer(true); else renderDrawerPanel();
  }

  // ---------- 3-D ----------
  function segMesh(scene, lbl, color, opacity, name, rank) {
    const mesh = scene.region_meshes[String(lbl)]; if (!mesh) return null;
    const v = mesh.verts, f = mesh.faces, nV = v.length / 3, nF = f.length / 3;
    const x = new Array(nV), y = new Array(nV), z = new Array(nV);
    for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
    const ii = new Array(nF), jj = new Array(nF), kk = new Array(nF);
    for (let i = 0; i < nF; i++) { ii[i] = f[i * 3]; jj[i] = f[i * 3 + 1]; kk[i] = f[i * 3 + 2]; }
    return { type: "mesh3d", x, y, z, i: ii, j: jj, k: kk, color, opacity, name, showlegend: true,
      flatshading: false, hoverinfo: "name", legendrank: rank || lbl,
      lighting: { ambient: 0.7, diffuse: 0.55, specular: 0.12, roughness: 0.9 } };
  }
  // object centroid + colour for the current embryo
  function objCentroid(e, key) {
    if (key === "polar") return e.polar;
    if (e.female == null) return null;                 // split → no M/P
    return key === "maternal" ? e.pron[e.female].com : e.pron[1 - e.female].com;
  }
  function render3D() {
    const s = state.scene, e = state.byId[state.currentId]; if (!s || !e) return;
    const showAll = els.showAll.checked, showBody = els.showBody.checked;
    const matLbl = e.female == null ? null : e.pron[e.female].label;
    const patLbl = e.female == null ? null : e.pron[1 - e.female].label;
    const splitLbls = e.split ? [e.pron[0].label, e.pron[1].label] : [];
    const traces = [];
    for (const lbl of s.mask_labels) {
      let color, op, name, rank, show;
      if (splitLbls.includes(lbl)) {
        color = SPLIT; op = 0.5; name = "Pronucleus (unassigned)"; rank = 10;
        show = showAll || state.obj !== "polar";
      } else if (lbl === matLbl) {
        color = F; op = 0.55; name = "Maternal ♀"; rank = 11; show = showAll || state.obj === "maternal";
      } else if (lbl === patLbl) {
        color = M; op = 0.55; name = "Paternal ♂"; rank = 12; show = showAll || state.obj === "paternal";
      } else if (lbl === e.polar_label) {
        color = PB; op = 0.4; name = "Polar body"; rank = 13; show = showAll || state.obj === "polar";
      } else {
        color = "#9aa3b2"; op = 0.06; name = `Segment ${lbl}`; rank = lbl; show = showBody;   // cell body
      }
      if (!show) continue;
      const t = segMesh(s, lbl, color, op, name, rank); if (t) traces.push(t);
    }
    // sperm
    const sp = e.sperm;
    traces.push({ type: "scatter3d", mode: "markers", name: "sperm",
      x: [sp[0]], y: [sp[1]], z: [sp[2]],
      marker: { size: 7, color: SP, symbol: "diamond", line: { color: "#fff", width: 1 } },
      hovertemplate: "sperm<extra></extra>", legendrank: 5 });
    // distance line to the selected object
    const c = objCentroid(e, state.obj);
    if (c) {
      const um = e.dist_um[state.obj];
      traces.push({ type: "scatter3d", mode: "lines+markers", name: `${um} µm → ${OBJ_LABEL[state.obj]}`,
        x: [sp[0], c[0]], y: [sp[1], c[1]], z: [sp[2], c[2]],
        line: { color: OBJ_COLOR[state.obj], width: 6, dash: "dash" },
        marker: { size: 3, color: OBJ_COLOR[state.obj] },
        hovertemplate: `${um} µm<extra></extra>`, legendrank: 6 });
    }
    Plotly.react(els.plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  function renderReadout() {
    const e = state.byId[state.currentId]; if (!e) return;
    const um = e.dist_um[state.obj];
    const split = e.split && state.obj !== "polar";
    const noPolar = state.obj === "polar" && e.polar == null;
    let big;
    if (split) big = `<span class="sp-split">no M/P consensus</span>`;
    else if (noPolar) big = `<span class="sp-split">no polar body</span>`;
    else big = `${um} <small>µm</small>`;
    els.readout.innerHTML =
      `<div class="sp-big">${big}</div>` +
      `<div class="sp-sub">sperm → <b>${OBJ_LABEL[state.obj]}</b> centre of mass<br>` +
      `pronuclei→COM clock τ = <b>${e.tau}</b> <small>(0 = PN formation → 1 = NEBD)</small>` +
      (e.gap != null ? `<br>interpronuclei distance = <b>${e.gap}</b> µm` : "") +
      (e.split ? `<br><span class="sp-split">pronuclei are an even split — maternal/paternal undefined</span>` : "") +
      `</div>`;
  }

  function renderLegend() {
    els.legend.innerHTML =
      `<div class="lg"><span class="di" style="background:${SP}"></span>sperm</div>` +
      `<div class="lg"><i style="background:${F}"></i>maternal pronucleus ♀</div>` +
      `<div class="lg"><i style="background:${M}"></i>paternal pronucleus ♂</div>` +
      `<div class="lg"><i style="background:${PB}"></i>polar body</div>`;
  }

  // ---------- drawer scatter ----------
  function openDrawer(open) {
    state.drawerOpen = open;
    els.drawer.dataset.open = open ? "true" : "false";
    els.drawerHandle.setAttribute("aria-expanded", String(open));
    if (open) renderDrawerPanel();
  }
  function wireDrawer() {
    els.drawerHandle.addEventListener("click", () => openDrawer(els.drawer.dataset.open !== "true"));
    $("#sp-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.drawerTab = b.dataset.tab;
      $("#sp-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#sp-panels").querySelectorAll(".xs-panel").forEach((p) =>
        (p.hidden = p.dataset.tab !== state.drawerTab));
      renderDrawerPanel();
    });
    $("#sp-pairctl").addEventListener("click", (ev) => {
      const b = ev.target.closest("[data-cmp]"); if (!b) return;
      state.pairKey = b.dataset.cmp; renderPair();
    });
    $("#sp-varctl").addEventListener("click", (ev) => {
      const b = ev.target.closest("[data-var]"); if (!b) return;
      state.pairVar = b.dataset.var; renderPair();
    });
    const rz = $("#drawer-resize"); let sh = 0, dy = 0, drag = false;
    rz.addEventListener("pointerdown", (ev) => { drag = true; sh = els.drawer.getBoundingClientRect().height; dy = ev.clientY; rz.setPointerCapture(ev.pointerId); ev.preventDefault(); });
    rz.addEventListener("pointermove", (ev) => { if (!drag) return;
      const h = Math.max(220, Math.min(window.innerHeight - 90, sh + (dy - ev.clientY)));
      els.drawer.style.setProperty("--drawer-h", h + "px"); requestAnimationFrame(() => { try { Plotly.Plots.resize(state.drawerTab === "pair" ? $("#sp-pair") : els.scatter); } catch (_) {} }); });
    const end = () => { drag = false; }; rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }

  function renderScatter() {
    if (els.drawer.hidden || !state.doc) return;
    els.drawerObj.textContent = OBJ_LABEL[state.obj];
    const col = OBJ_COLOR[state.obj];
    const clk = activeClock();
    const hint = $("#drawer-hint");
    if (hint) hint.textContent = `one point per sperm-labelled zygote · x = ${clk.key === "tau" ? "calibrated τ" : "interpronuclei distance"}, y = sperm-to-object distance`;
    const pts = state.doc.embryos
      .map((e) => ({ e, y: e.dist_um[state.obj], x: e[clk.field] }))
      .filter((p) => p.y != null && p.x != null);
    const xs = pts.map((p) => p.x), ys = pts.map((p) => p.y);
    const cur = state.currentId;
    const isCur = pts.map((p) => p.e.id === cur);

    const markers = {
      type: "scatter", mode: "markers", name: OBJ_LABEL[state.obj],
      x: xs, y: ys,
      error_x: clk.hasInterval ? {
        type: "data", symmetric: false,
        array: pts.map((p) => Math.max(0, p.e.tau_hi - p.e.tau)),
        arrayminus: pts.map((p) => Math.max(0, p.e.tau - p.e.tau_lo)),
        color: "rgba(120,130,150,0.35)", thickness: 1, width: 0,
      } : undefined,
      marker: {
        size: isCur.map((c) => (c ? 15 : 9)),
        color: col, opacity: 0.85,
        line: { color: isCur.map((c) => (c ? "#111827" : "#fff")), width: isCur.map((c) => (c ? 2.5 : 1)) },
      },
      text: pts.map((p) => `${p.e.label}  ·  ${clk.key === "tau" ? "τ " + p.e.tau : p.e.gap + " µm gap"}  ·  ${p.y} µm`),
      hovertemplate: "%{text}<extra></extra>", customdata: pts.map((p) => p.e.id),
    };
    const traces = [markers];
    // OLS trend
    const fit = ols(xs, ys);
    if (fit) {
      const gx = [Math.min(...xs), Math.max(...xs)];
      traces.push({ type: "scatter", mode: "lines", name: "linear fit", hoverinfo: "skip",
        x: gx, y: gx.map((x) => fit.b0 + fit.b1 * x),
        line: { color: "#111827", width: 2, dash: "dot" } });
    }
    const layout = {
      margin: { l: 54, r: 12, t: 8, b: 42 }, showlegend: false,
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "#fcfdfe",
      xaxis: { title: { text: clk.axis, font: { size: 12 } },
        gridcolor: "#eef1f5", zeroline: false,
        range: clk.key === "tau" ? [-0.03, Math.max(...xs) * 1.05 + 0.02]
                                  : [Math.min(...xs) * 0.9, Math.max(...xs) * 1.05] },
      yaxis: { title: { text: `sperm → ${OBJ_LABEL[state.obj]} distance (µm)`, font: { size: 12 } },
        gridcolor: "#eef1f5", zeroline: false, rangemode: "tozero" },
    };
    Plotly.react(els.scatter, traces, layout, { responsive: true, displaylogo: false,
      modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"] });
    els.scatter.on("plotly_click", (ev) => {
      const id = ev.points && ev.points[0] && ev.points[0].customdata;
      if (id && state.byId[id]) selectEmbryo(id);
    });
    // stats
    const r = fit ? fit.r : null;
    const dropped = state.doc.embryos.length - pts.length;
    els.stat.innerHTML = fit
      ? `n = <b>${pts.length}</b> · Pearson r = <b>${r.toFixed(2)}</b> · ` +
        `slope <b>${fit.b1.toFixed(clk.key === "tau" ? 0 : 2)}</b> µm/${clk.key === "tau" ? "τ" : "µm"} · p ${fmtP(fit.p)}` +
        (Math.abs(r) > 0 && fit.p < 0.05 ? ` <span class="sig">significant</span>` : "")
      : `n = <b>${pts.length}</b>`;
    els.note.innerHTML =
      `Distance is measured in 3-D from the sperm to the ${OBJ_LABEL[state.obj]} centre of mass (COM). ` +
      (clk.key === "tau"
        ? `τ is the calibrated pronuclei-to-COM clock (${state.doc.model_version}); horizontal bars are its 95% interval. `
        : `The x-axis is the interpronuclei distance — the raw surface gap between the two pronuclei (larger = earlier). `) +
      (dropped > 0
        ? `<b>${dropped}</b> zygote${dropped === 1 ? "" : "s"} omitted ` +
          (state.obj === "polar" ? "(no polar body labelled)." : "(pronuclei are an even split, so maternal/paternal is undefined).")
        : "All sperm-labelled zygotes shown.");
  }

  // ---------- maternal vs paternal, paired within the zygote (figure 4.16) ----------
  //
  // The two pronuclei of one cell are compared against each other, so embryo size, orientation
  // and segmentation quality all divide out and the whole analysis is a paired Wilcoxon.
  //
  // ⚠️ EACH COMPARISON IS ONE OF THE FOUR TESTS THAT DECIDED WHICH PRONUCLEUS IS WHICH, re-asked
  // as a measurement. The build ships each one twice — once on the full consensus, once on a
  // consensus with that test removed — and this panel DEFAULTS TO THE LEAVE-ONE-OUT column,
  // because that is the only one that is not partly a restatement of its own labelling.
  async function loadPairing() {
    if (state.pair !== undefined) return state.pair;
    try { state.pair = await (await fetch("data/sperm_pairing.json")).json(); }
    catch (_) { state.pair = null; }
    return state.pair;
  }

  async function renderPair() {
    const host = $("#sp-pair");
    const doc = await loadPairing();
    if (!doc) {
      host.innerHTML = `<div class="sp-empty">Run <code>python3 build_sperm_pairing.py</code>.</div>`;
      return;
    }
    const keys = Object.keys(doc.comparisons);
    if (!state.pairKey || !doc.comparisons[state.pairKey]) state.pairKey = keys[0];
    const ctl = $("#sp-pairctl");
    if (!ctl.dataset.built) {
      ctl.innerHTML = keys.map((k) =>
        `<button type="button" data-cmp="${k}"${k === state.pairKey ? ' class="on"' : ""}>${doc.comparisons[k].label}</button>`).join("");
      ctl.dataset.built = "1";
    }
    ctl.querySelectorAll("[data-cmp]").forEach((b) =>
      b.classList.toggle("on", b.dataset.cmp === state.pairKey));
    $("#sp-varctl").querySelectorAll("[data-var]").forEach((b) =>
      b.classList.toggle("on", b.dataset.var === state.pairVar));

    const C = doc.comparisons[state.pairKey];
    const v = C.variants[state.pairVar];
    const rows = v.rows;
    const naive = C.variants.all, loo = C.variants.loo;

    $("#sp-circ").innerHTML =
      `<b>${C.circular}</b> On the full consensus it gives ` +
      `${naive.n_maternal_larger}/${naive.n} and P ${fmtP(naive.p)}; with ` +
      `${C.drops.join(" and ")} removed, ${loo.n_maternal_larger}/${loo.n} and P ${fmtP(loo.p)}. ` +
      `${C.expect}`;

    // one line per zygote, maternal on the left and paternal on the right — the paired plot,
    // because the pairing IS the analysis and a pair of box plots would hide it
    const traces = [];
    rows.forEach((r, i) => traces.push({
      type: "scatter", mode: "lines+markers", x: [0, 1], y: [r.m, r.p],
      line: { color: r.m > r.p ? "rgba(214,51,108,0.45)" : "rgba(37,99,235,0.45)", width: 1.4 },
      marker: { size: r.manual ? 11 : 7, symbol: r.manual ? "diamond" : "circle",
                color: [F, M], line: { color: "#fff", width: 1 } },
      hovertemplate: `${r.label}<br>maternal %{y:.1f}<extra></extra>`,
      text: [r.label, r.label], showlegend: false,
      customdata: [r.id, r.id],
    }));
    const med = (a) => { const s = a.slice().sort((x, y) => x - y); const n = s.length;
      return n ? (n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2) : 0; };
    const mm = med(rows.map((r) => r.m)), mp = med(rows.map((r) => r.p));
    traces.push({ type: "scatter", mode: "lines+markers", x: [0, 1], y: [mm, mp],
      line: { color: "#111827", width: 3 }, marker: { size: 13, color: "#111827" },
      hovertemplate: "median %{y:.1f}<extra></extra>", showlegend: false });

    Plotly.react(host, traces, {
      margin: { l: 66, r: 12, t: 10, b: 40 }, showlegend: false,
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "#fcfdfe",
      xaxis: { range: [-0.35, 1.35], tickvals: [0, 1],
        ticktext: ["maternal ♀", "paternal ♂"], tickfont: { size: 12 },
        gridcolor: "#f4f6f9", zeroline: false },
      yaxis: { title: { text: `${C.label} (${C.unit})`, font: { size: 12 } },
        gridcolor: "#eef1f5", zeroline: false, rangemode: "tozero" },
    }, { responsive: true, displaylogo: false,
         modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"] });
    host.on("plotly_click", (ev) => {
      const id = ev.points && ev.points[0] && ev.points[0].customdata;
      if (id && state.byId[id]) selectEmbryo(id);
    });

    // "larger" is right for a volume and wrong for a distance, where the reader wants "farther"
    const more = C.unit === "µm³" ? "larger" : "farther";
    const dir = v.median_diff > 0 ? `maternal is ${more}` : `paternal is ${more}`;
    $("#sp-pairnote").innerHTML =
      `<b>n = ${v.n}</b> zygotes · median maternal − paternal ` +
      `<b>${v.median_diff > 0 ? "+" : ""}${v.median_diff.toFixed(C.unit === "µm³" ? 0 : 2)} ${C.unit}</b> ` +
      `(${dir} in ${Math.max(v.n_maternal_larger, v.n_paternal_larger)} of ${v.n}) · ` +
      `paired Wilcoxon signed-rank ${v.exact ? "(exact)" : "(normal approximation)"} P ${fmtP(v.p)}. ` +
      `Each thin line is one zygote; the heavy line is the median pair. ` +
      (v.n_manual ? `Diamonds are the ${v.n_manual} hand-made call${v.n_manual === 1 ? "" : "s"} — ` +
        `those are not votes, so they survive every leave-one-out unchanged. ` : "") +
      (v.dropped.length ? `${v.dropped.length} zygote${v.dropped.length === 1 ? "" : "s"} omitted ` +
        `(${[...new Set(v.dropped.map((d) => d.reason))].join("; ")}).` : "");
  }

  function renderDrawerPanel() {
    const title = document.querySelector(".drawer-title"), hint = $("#drawer-hint");
    if (state.drawerTab === "pair") {
      if (title) title.innerHTML = 'Maternal vs paternal <span class="drawer-gene">paired</span>';
      if (hint) hint.textContent = "one line per zygote · the two pronuclei of the same cell, " +
        "compared against each other";
      renderPair();
    } else {
      if (title) title.innerHTML =
        'Sperm distance vs pseudotime <span class="drawer-gene" id="drawer-obj"></span>';
      els.drawerObj = $("#drawer-obj");
      renderScatter();
    }
  }

  // ---------- export + stats ----------
  function exportCsv() {
    const rows = [["zygote", "tau", "tau_lo", "tau_hi", "interpronuclei_gap_um", "sperm_to_polar_um", "sperm_to_maternal_um", "sperm_to_paternal_um", "split"]];
    state.doc.embryos.forEach((e) => rows.push([e.id, e.tau, e.tau_lo, e.tau_hi, e.gap ?? "",
      e.dist_um.polar ?? "", e.dist_um.maternal ?? "", e.dist_um.paternal ?? "", e.split]));
    const blob = new Blob([rows.map((r) => r.join(",")).join("\n")], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = "sperm_pseudotime.csv";
    document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(a.href);
  }
  function ols(xs, ys) {
    const n = xs.length; if (n < 3) return null;
    const mx = xs.reduce((s, x) => s + x, 0) / n, my = ys.reduce((s, y) => s + y, 0) / n;
    let sxx = 0, sxy = 0, syy = 0;
    for (let i = 0; i < n; i++) { const dx = xs[i] - mx, dy = ys[i] - my; sxx += dx * dx; sxy += dx * dy; syy += dy * dy; }
    if (sxx === 0 || syy === 0) return null;
    const b1 = sxy / sxx, b0 = my - b1 * mx, r = sxy / Math.sqrt(sxx * syy);
    const t = r * Math.sqrt((n - 2) / (1 - r * r));
    return { b0, b1, r, p: studentTwoSided(Math.abs(t), n - 2) };
  }
  function studentTwoSided(t, df) { return clamp01(betai(df / 2, 0.5, df / (df + t * t))); }
  const clamp01 = (x) => Math.max(0, Math.min(1, x));
  function betai(a, b, x) {
    if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(gln(a + b) - gln(a) - gln(b) + a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
  }
  function betacf(a, b, x) {
    const FPMIN = 1e-300; let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
    for (let m = 1; m <= 200; m++) { const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d;
      const del = d * c; h *= del; if (Math.abs(del - 1) < 3e-12) break; }
    return h;
  }
  function gln(x) {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp); let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += c[j] / ++y; return -tmp + Math.log(2.5066282746310005 * ser / x);
  }
  function fmtP(p) { return p < 0.0001 ? "< 0.0001" : p.toFixed(p < 0.01 ? 4 : 3); }
})();
