/* Across the Stages — figures 8.3 / 8.6 / 9.1 / 9.2.
 *
 * How unevenly does a gene sit between the two halves of an embryo, at each of the three stages?
 * The halves are the best meridional plane in the zygote and the two blastomeres at 2-cell, so the
 * object being split changes — which is why the null matters as much as the fold.
 *
 * Data: data/stages.json.gz (build_stages.py).
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const CCOL = { "drops after the zygote": "#0f766e", "peaks at early 2-cell": "#c2540b",
                 "rises into late 2-cell": "#5b3fa6" };
  const FATE = { retained: "#dc2626", lost: "#2563eb", gained: "#0f766e", other: "#cbd5e1" };
  const C_NULL = "#94a3b8", C_SEL = "#111827";

  const state = { data: null, gene: null, tab: "heat", filter: "measured", showNull: true, find: "" };
  const meta = () => state.data.meta;
  const ST = () => meta().stages;
  const LAB = () => meta().stageLabel;

  const measured = (g) => ST().some((s) => g[s] && g[s].nMeas > 0);
  function shown() {
    const G = state.data.genes;
    if (state.filter === "all") return G;
    if (state.filter === "clustered") return G.filter((g) => g.cluster);
    if (state.filter === "curated") return G.filter((g) => g.refPct);
    return G.filter(measured);
  }
  const recOf = (g) => state.data.genes.find((r) => r.g === g) || null;
  /** How far a gene's fold sits above what counting noise alone would give, worst stage first. */
  const excess = (g) => Math.max(...ST().map((s) => (g[s] ? g[s].excess : -Infinity)));

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await V.loadGz("data/stages.json.gz"); }
    catch (err) {
      $("#plot").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_stages.py</code>.</div></div>`;
      return;
    }
    const m = meta();
    $("#gene-count").textContent =
      `${m.n_genes} genes · zygote = best meridional plane · 2-cell = the two blastomeres` +
      (m.validation && m.validation.available
        ? ` · matches the reference to ${m.validation.max_abs_diff.zygote}` : "");
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    const s = shown().slice().sort((a, b) => excess(b) - excess(a))[0];
    state.gene = s && s.g;
    refresh();
  })();

  // ───────── the stage plot (8.3) ─────────
  function renderMain() {
    const G = shown();
    const x = ST().map((s) => LAB()[s]);
    const traces = [];
    G.forEach((g) => {
      const y = ST().map((s) => (g[s] ? g[s].fold : null));
      if (y.every((v) => v == null)) return;
      // the selected gene must be visible even when it carries no cluster or fate,
      // otherwise clicking an unclassified line appears to do nothing
      const col = g.g === state.gene ? C_SEL
        : g.cluster ? CCOL[g.cluster] : (g.group !== "other" ? FATE[g.group] : "#dde3ea");
      traces.push({
        type: "scatter", mode: "lines+markers", x, y, name: g.g,
        line: { color: col, width: g.g === state.gene ? 3 : 1.1 },
        marker: { size: g.g === state.gene ? 8 : 4, color: col },
        opacity: g.g === state.gene ? 1 : (g.cluster || g.group !== "other" ? 0.75 : 0.30),
        hovertemplate: `<b>${g.g}</b><br>%{x}<br>fold %{y:.3f}<extra></extra>`,
        showlegend: false,
      });
    });
    if (state.showNull) {
      // the median null at each stage: below this line, a fold is arithmetic, not asymmetry
      const med = ST().map((s) => {
        const v = G.filter((g) => g[s]).map((g) => g[s].null).sort((a, b) => a - b);
        return v.length ? v[Math.floor(v.length / 2)] : null;
      });
      traces.push({ type: "scatter", mode: "lines+markers", x, y: med, name: "median null",
        line: { color: C_NULL, width: 2.5, dash: "dash" }, marker: { size: 8, color: C_NULL },
        hovertemplate: "median count-matched null<br>%{x}: %{y:.3f}<extra></extra>",
        showlegend: true });
    }
    const lay = {
      margin: { l: 66, r: 18, t: 20, b: 46 }, showlegend: true,
      legend: { x: 0.99, y: 0.99, xanchor: "right", yanchor: "top", font: { size: 10 },
                bgcolor: "rgba(255,255,255,0.75)" },
      xaxis: { type: "category", gridcolor: "#eef1f5", tickfont: { size: 11 } },
      yaxis: { title: { text: "between-half fold change", font: { size: 11 } },
        gridcolor: "#eef1f5", tickfont: { size: 10 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [{ type: "line", x0: -0.5, x1: 2.5, y0: 1, y1: 1,
                 line: { color: "#111827", width: 1, dash: "dot" } }],
      annotations: [{ x: -0.45, y: 1, xanchor: "left", yanchor: "bottom", showarrow: false,
        text: "1.0 = the two halves carry equal density", font: { size: 9, color: "#64748b" } }],
    };
    Plotly.react($("#plot"), traces, lay,
      { displaylogo: false, responsive: true,
        modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
        toImageButtonOptions: { format: "png", scale: 4 } });
    $("#plot").on("plotly_click", (ev) => {
      const p = ev.points && ev.points[0];
      if (p && p.data && p.data.name && p.data.name !== "median null") {
        state.gene = p.data.name; refresh();
      }
    });
  }

  // ───────── readout + list ─────────
  function renderReadout() {
    const g = state.gene && recOf(state.gene);
    const el = $("#sg-readout");
    if (!g) { el.innerHTML = `<div class="sg-hint">Click a line to select a gene.</div>`; return; }
    const L = [`<div class="sg-gene">${g.g}</div>`];
    if (g.cluster) L.push(`<div class="sg-tag" style="background:${CCOL[g.cluster]}">${g.cluster}</div>`);
    if (g.group !== "other") L.push(`<div class="sg-tag" style="background:${FATE[g.group]}">${g.group}</div>`);
    L.push(`<table class="sg-mini"><tr><th></th><th>fold</th><th>null</th><th>n</th><th>meas</th></tr>` +
      ST().map((s) => {
        const r = g[s];
        if (!r) return `<tr><td>${LAB()[s]}</td><td colspan="4">not on any panel</td></tr>`;
        const strong = r.excess > 0.05;
        return `<tr><td>${LAB()[s]}</td>` +
          `<td class="${strong ? "hi" : ""}">${r.fold.toFixed(3)}</td>` +
          `<td class="nul">${r.null.toFixed(3)}</td>` +
          `<td>${r.n}</td><td class="${r.nMeas ? "" : "zero"}">${r.nMeas}</td></tr>`;
      }).join("") + `</table>`);
    const anyBelow = ST().some((s) => g[s] && g[s].excess <= 0);
    if (anyBelow) {
      L.push(`<div class="sg-warn">At least one stage's fold sits <b>at or below its own
        count-matched null</b> — at that stage this is arithmetic, not asymmetry.</div>`);
    }
    const noMeas = ST().filter((s) => g[s] && g[s].nMeas === 0);
    if (noMeas.length) {
      L.push(`<div class="sg-warn">Never detected above the floor at
        ${noMeas.map((s) => LAB()[s]).join(", ")} — the fold there is exactly 1.0 by the zero
        rule, not a measurement.</div>`);
    }
    el.innerHTML = L.join("");
  }

  function renderRank() {
    const rows = shown().filter((g) =>
      !state.find || g.g.toLowerCase().includes(state.find.toLowerCase()))
      .slice().sort((a, b) => excess(b) - excess(a));
    $("#sg-rank-desc").innerHTML =
      `${shown().length} shown of ${state.data.genes.length}. Sorted by how far the best stage's
       fold rises <b>above its own null</b>, which is the only fair way to rank across stages.`;
    $("#sg-rank").innerHTML =
      `<div class="sg-rank-head"><span>#</span><span>gene</span><span>zyg</span><span>e2c</span><span>l2c</span></div>` +
      rows.slice(0, 400).map((g, i) => `<div class="sg-rank-row${g.g === state.gene ? " on" : ""}"
          data-g="${g.g}" title="${g.g}${g.cluster ? " · " + g.cluster : ""}${g.group !== "other" ? " · " + g.group : ""}">
        <span class="n">${i + 1}</span>
        <span class="e">${g.g}</span>` +
        ST().map((s) => `<span class="f${g[s] && g[s].excess > 0.05 ? " hi" : ""}">${
          g[s] ? g[s].fold.toFixed(2) : "–"}</span>`).join("") + `</div>`).join("");
    $("#sg-rank").querySelectorAll(".sg-rank-row").forEach((el) =>
      el.addEventListener("click", () => { state.gene = el.dataset.g; refresh(); }));
  }

  // ───────── drawer panels ─────────
  const CFG = { displaylogo: false, responsive: true,
                modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                toImageButtonOptions: { format: "png", scale: 4 } };
  const empty = (el, msg) => { Plotly.purge(el); el.innerHTML = `<div class="sg-empty">${msg}</div>`; };

  function renderHeat() {
    const el = $("#sg-heat");
    const G = state.data.genes.filter((g) => g.cluster)
      .sort((a, b) => (a.cluster.localeCompare(b.cluster)) ||
                      ((b.zygote ? b.zygote.fold : 0) - (a.zygote ? a.zygote.fold : 0)));
    if (!G.length) return empty(el, "No gene carries an imported trajectory cluster.");
    el.innerHTML = "";
    const z = G.map((g) => ST().map((s) => (g[s] ? g[s].fold : null)));
    Plotly.newPlot(el, [{
      type: "heatmap", z, x: ST().map((s) => LAB()[s]), y: G.map((g) => g.g),
      colorscale: "RdBu", reversescale: true, zmid: 1,
      colorbar: { title: { text: "fold", font: { size: 10 } }, thickness: 12 },
      hovertemplate: "%{y}<br>%{x}<br>fold %{z:.3f}<extra></extra>",
    }], {
      margin: { l: 92, r: 14, t: 14, b: 44 },
      xaxis: { side: "top", tickfont: { size: 10 } },
      yaxis: { autorange: "reversed", tickfont: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, CFG);
    const byC = {};
    G.forEach((g) => (byC[g.cluster] = (byC[g.cluster] || 0) + 1));
    $("#sg-heat-sub").textContent = "· " + Object.entries(byC)
      .map(([k, v]) => `${v} ${k}`).join(" · ") + " · grouping imported";
  }

  function renderPct() {
    const el = $("#sg-pct");
    const G = state.data.genes.filter((g) => g.refPct);
    if (!G.length) return empty(el, "The reference's percentile table is not present.");
    el.innerHTML = "";
    const byFate = {};
    G.forEach((g) => (byFate[g.group] = byFate[g.group] || []).push(g));
    const traces = Object.entries(byFate).map(([f, arr]) => ({
      type: "scatter", mode: "markers", name: f,
      x: arr.map((g) => g.refPct.e), y: arr.map((g) => g.refPct.l),
      marker: { size: f === "other" ? 7 : 11, color: FATE[f] || "#94a3b8",
                opacity: f === "other" ? 0.55 : 0.95, line: { color: "#fff", width: 1 } },
      text: arr.map((g) => `${g.g} · ${g.refPct.eN}→${g.refPct.lN} embryos`),
      hovertemplate: "%{text}<br>early %{x:.2f} → late %{y:.2f}<extra></extra>",
    }));
    traces.push({ type: "scatter", mode: "lines", x: [0, 1], y: [0, 1], name: "no change",
      line: { color: "#111827", width: 1.2, dash: "dash" }, hoverinfo: "skip" });
    Plotly.newPlot(el, traces, {
      margin: { l: 62, r: 14, t: 14, b: 46 }, showlegend: true,
      legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      xaxis: { title: { text: "early 2-cell — count-matched percentile", font: { size: 10 } },
        range: [-0.03, 1.03], gridcolor: "#eef1f5", tickfont: { size: 9 }, scaleanchor: "y", scaleratio: 1 },
      yaxis: { title: { text: "late 2-cell", font: { size: 10 } }, range: [-0.03, 1.03],
        gridcolor: "#eef1f5", tickfont: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, CFG);
    $("#sg-pct-sub").textContent =
      `· ${G.length} genes measurable at both 2-cell stages · above the line = polarization gained`;
  }

  function renderFate() {
    const el = $("#sg-fate");
    const G = state.data.genes.filter((g) => g.refPct);
    if (!G.length) return empty(el, "The reference's percentile table is not present.");
    el.innerHTML = "";
    const traces = [];
    G.filter((g) => g.group === "other").forEach((g) => traces.push({
      type: "scatter", mode: "lines", x: ["early 2-cell", "late 2-cell"],
      y: [g.refPct.e, g.refPct.l], line: { color: "#dde3ea", width: 1 },
      hovertemplate: `${g.g}<extra></extra>`, showlegend: false }));
    ["retained", "lost", "gained"].forEach((f) => {
      const arr = G.filter((g) => g.group === f);
      arr.forEach((g, i) => traces.push({
        type: "scatter", mode: "lines+markers", x: ["early 2-cell", "late 2-cell"],
        y: [g.refPct.e, g.refPct.l],
        line: { color: FATE[f], width: 2.2 }, marker: { size: 7, color: FATE[f] },
        name: f, legendgroup: f, showlegend: i === 0,
        hovertemplate: `<b>${g.g}</b> · ${f}<br>%{x} %{y:.2f}<extra></extra>` }));
    });
    Plotly.newPlot(el, traces, {
      margin: { l: 62, r: 14, t: 20, b: 44 }, showlegend: true,
      legend: { orientation: "h", x: 0, y: 1.03, xanchor: "left", yanchor: "bottom", font: { size: 10 } },
      xaxis: { type: "category", gridcolor: "#eef1f5", tickfont: { size: 11 } },
      yaxis: { title: { text: "count-matched percentile", font: { size: 10 } }, range: [-0.03, 1.03],
        gridcolor: "#eef1f5", tickfont: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, CFG);
    const c = {};
    G.forEach((g) => (c[g.group] = (c[g.group] || 0) + 1));
    $("#sg-fate-sub").textContent =
      `· ${c.retained || 0} retained · ${c.lost || 0} lost · ${c.gained || 0} gained (derived) · ` +
      `${c.other || 0} unclassified in grey`;
  }

  function renderTable() {
    const rows = shown().slice().sort((a, b) => excess(b) - excess(a));
    $("#sg-table").innerHTML =
      `<table class="sg-tab"><thead><tr><th>#</th><th>gene</th><th>group</th>` +
      ST().map((s) => `<th colspan="4">${LAB()[s]}</th>`).join("") + `</tr><tr><th></th><th></th><th></th>` +
      ST().map(() => `<th>fold</th><th>null</th><th>n</th><th>meas</th>`).join("") +
      `</tr></thead><tbody>` +
      rows.map((g, i) => `<tr class="${g.g === state.gene ? "on" : ""}" data-g="${g.g}">
        <td>${i + 1}</td><td class="g">${g.g}</td>
        <td class="grp">${g.cluster || (g.group !== "other" ? g.group : "")}</td>` +
        ST().map((s) => {
          const r = g[s];
          if (!r) return `<td colspan="4" class="na">–</td>`;
          return `<td class="${r.excess > 0.05 ? "hi" : ""}">${r.fold.toFixed(3)}</td>` +
                 `<td class="nul">${r.null.toFixed(3)}</td><td>${r.n}</td>` +
                 `<td class="${r.nMeas ? "" : "zero"}">${r.nMeas}</td>`;
        }).join("") + `</tr>`).join("") + `</tbody></table>`;
    $("#sg-table").querySelectorAll("tr[data-g]").forEach((tr) =>
      tr.addEventListener("click", () => { state.gene = tr.dataset.g; refresh(); }));
    $("#sg-table-sub").textContent = `· ${rows.length} genes`;
  }

  const RENDER = { heat: renderHeat, pct: renderPct, fate: renderFate, table: renderTable };
  function renderPanel() {
    const fn = RENDER[state.tab];
    if (!fn) return;
    try { fn(); } catch (err) { console.error("[stages]", state.tab, err); }
  }
  function refresh() {
    $("#drawer-gene").textContent = state.gene || "";
    renderMain(); renderReadout(); renderRank(); renderPanel();
  }

  // ───────── chrome ─────────
  function dragResize(el, move, start) {
    if (!el) return;
    el.addEventListener("pointerdown", (ev) => {
      el._d = Object.assign({ x: ev.clientX, y: ev.clientY }, start());
      el.setPointerCapture(ev.pointerId); ev.preventDefault(); el.classList.add("dragging");
    });
    el.addEventListener("pointermove", (ev) => { if (el._d) move(el._d, ev); });
    const end = (ev) => {
      el._d = null; el.classList.remove("dragging");
      try { el.releasePointerCapture(ev.pointerId); } catch (_) {}
    };
    el.addEventListener("pointerup", end); el.addEventListener("pointercancel", end);
  }

  function wire() {
    $("#sg-filter").addEventListener("change", (e) => { state.filter = e.target.value; refresh(); });
    $("#t-null").addEventListener("change", (e) => { state.showNull = e.target.checked; renderMain(); });
    $("#sg-find").addEventListener("input", (e) => { state.find = e.target.value; renderRank(); });
    $("#sg-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#sg-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#sg-panels").querySelectorAll(".xs-panel").forEach((p) =>
        (p.hidden = p.dataset.tab !== state.tab));
      renderPanel();
    });
    const openDrawer = (open) => {
      $("#drawer").dataset.open = open ? "true" : "false";
      $("#drawer-handle").setAttribute("aria-expanded", String(open));
      if (open) setTimeout(renderPanel, 30);
    };
    $("#drawer-handle").addEventListener("click", () => openDrawer($("#drawer").dataset.open !== "true"));
    const rh = $("#rdrawer-handle");
    rh.addEventListener("click", () => {
      const r = $("#rdrawer"), open = r.dataset.open !== "true";
      r.dataset.open = String(open); rh.setAttribute("aria-expanded", String(open));
    });
    dragResize($("#drawer-resize"), (d, ev) => {
      $("#drawer").style.setProperty("--drawer-h",
        Math.max(200, Math.min(window.innerHeight - 120, d.v + (d.y - ev.clientY))) + "px");
    }, () => ({ v: $("#drawer-body").getBoundingClientRect().height }));
    dragResize($("#rdrawer-resize"), (d, ev) => {
      $("#rdrawer").style.setProperty("--rdrawer-w",
        Math.max(260, Math.min(window.innerWidth - 80, d.v - (ev.clientX - d.x))) + "px");
    }, () => ({ v: $("#rdrawer").getBoundingClientRect().width }));
    V.wireWindow($("#controls"), $("#controls-header"),
                 [...$("#controls").querySelectorAll(".rz")], "stages_controls_box");
    window.addEventListener("resize", () => {
      ["#plot", "#sg-heat", "#sg-pct", "#sg-fate"].forEach((s) => {
        try { Plotly.Plots.resize($(s)); } catch (_) {}
      });
    });
  }
})();
