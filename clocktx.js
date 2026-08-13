/* Transcriptome vs the Clock — figures 4.8 / 4.11 / 5.4.
 *
 * Does any gene's presence in the transcriptome move with pseudotime? Every point in the volcano
 * is a gene: Spearman rho of its probeset-centred value against tau, on the zygotes carrying at
 * least the entry floor of its transcripts.
 *
 * Two axes of choice, both live: SHARE (fraction of the mix) vs CONCENTRATION (per µm³ of
 * cytoplasm) — a shrinking cell raises the second at constant first — and whether the polar body
 * is counted. build_clocktx.py computes all four combinations; the page just switches between them.
 *
 * Data: data/clocktx.json.gz (build_clocktx.py).
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const C_UP = "#b45309", C_DOWN = "#0f766e", C_NS = "#cbd5e1", C_SEL = "#dc2626";

  const state = { data: null, ycol: "share", region: "main", gene: null, tab: "one",
                  fdr: false, labels: true, find: "" };
  const meta = () => state.data.meta;
  const variant = () => state.data.variants[`${state.region}.${state.ycol}`];
  const tauOf = {};
  const embOf = {};

  const YLAB = () => (state.ycol === "share" ? "share of the transcriptome"
                                             : "transcripts per µm³ of cytoplasm");

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await V.loadGz("data/clocktx.json.gz"); }
    catch (err) {
      $("#plot").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_clocktx.py</code>.</div></div>`;
      return;
    }
    const m = meta();
    state.data.embryos.forEach((e) => { tauOf[e.id] = e.tau; embOf[e.id] = e; });
    $("#embryo-count").textContent =
      `${m.n_embryos} zygotes with a τ · ${m.n_genes} genes · floor ${m.params.MIN_COUNT} transcripts`;
    $("#ct-floor").textContent = m.params.MIN_COUNT;
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    const first = variant().genes.find((g) => g.p !== null);
    state.gene = (m.exemplars || []).find((g) => variant().genes.some((r) => r.g === g)) ||
                 (first && first.g);
    refresh();
  })();

  const testable = () => variant().genes.filter((r) => r.p !== null);
  const recOf = (g) => variant().genes.find((r) => r.g === g) || null;

  // ───────── the volcano ─────────
  function renderVolcano() {
    const rows = testable();
    const sig = (r) => r.p < 0.05;
    const x = rows.map((r) => r.rho);
    const y = rows.map((r) => -Math.log10(Math.max(r.p, 1e-12)));
    const col = rows.map((r) => (!sig(r) ? C_NS : r.rho > 0 ? C_UP : C_DOWN));
    const size = rows.map((r) => 5 + 5 * Math.min(1, Math.log10(Math.max(r.medCount, 1)) / 3));
    const traces = [{
      type: "scattergl", mode: "markers", x, y,
      marker: { size, color: col, opacity: 0.82, line: { color: "#fff", width: 0.5 } },
      text: rows.map((r) => r.g),
      customdata: rows.map((r) => [r.n, r.medCount, r.p, r.fdr == null ? NaN : r.fdr]),
      hovertemplate: "<b>%{text}</b><br>ρ %{x:.3f}<br>P %{customdata[2]:.2g} · " +
                     "q %{customdata[3]:.2g}<br>%{customdata[0]} zygotes · median " +
                     "%{customdata[1]:.0f} transcripts<extra></extra>",
      name: "genes",
    }];
    const sel = state.gene && rows.find((r) => r.g === state.gene);
    if (sel) {
      traces.push({
        type: "scattergl", mode: "markers", x: [sel.rho],
        y: [-Math.log10(Math.max(sel.p, 1e-12))],
        marker: { size: 15, color: "rgba(0,0,0,0)", line: { color: C_SEL, width: 2.5 } },
        hoverinfo: "skip", name: sel.g, showlegend: false });
    }
    const lay = {
      margin: { l: 62, r: 18, t: 16, b: 48 }, showlegend: false,
      xaxis: { title: { text: "Spearman ρ  (probeset-centred value vs τ)", font: { size: 11 } },
        range: [-1.08, 1.08], zeroline: true, zerolinecolor: "#94a3b8", gridcolor: "#eef1f5",
        tickfont: { size: 10 } },
      yaxis: { title: { text: "−log₁₀ P", font: { size: 11 } }, rangemode: "tozero",
        gridcolor: "#eef1f5", tickfont: { size: 10 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [{ type: "line", x0: -1.08, x1: 1.08, y0: -Math.log10(0.05), y1: -Math.log10(0.05),
                 line: { color: "#94a3b8", width: 1, dash: "dot" } }],
      annotations: [
        { x: -1.05, y: -Math.log10(0.05), xanchor: "left", yanchor: "bottom", showarrow: false,
          text: "P = 0.05", font: { size: 9, color: "#94a3b8" } },
        { x: -1.0, y: 1.0, xref: "paper", yref: "paper", xanchor: "left", yanchor: "top",
          showarrow: false, text: "", font: { size: 10 } },
      ],
    };
    // the FDR line only exists if anything survives it
    if (state.fdr) {
      const q = rows.filter((r) => r.fdr != null && r.fdr < 0.05).map((r) => r.p);
      if (q.length) {
        const cut = -Math.log10(Math.max(...q));
        lay.shapes.push({ type: "line", x0: -1.08, x1: 1.08, y0: cut, y1: cut,
                          line: { color: C_UP, width: 1.4, dash: "dash" } });
        lay.annotations.push({ x: 1.05, y: cut, xanchor: "right", yanchor: "bottom",
          showarrow: false, text: `FDR 0.05 (${q.length} genes)`,
          font: { size: 9, color: C_UP } });
      } else {
        lay.annotations.push({ x: 0.5, y: 1.0, xref: "paper", yref: "paper", yanchor: "top",
          showarrow: false, text: "nothing survives FDR < 0.05",
          font: { size: 10, color: C_UP } });
      }
    }
    if (state.labels) {
      rows.slice(0, 12).forEach((r) => lay.annotations.push({
        x: r.rho, y: -Math.log10(Math.max(r.p, 1e-12)), text: r.g, showarrow: false,
        yshift: 9, font: { size: 9, color: "#334155" } }));
    }
    Plotly.react($("#plot"), traces, lay,
                 { displaylogo: false, responsive: true,
                   modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                   toImageButtonOptions: { format: "png", scale: 4 } });
    $("#plot").removeAllListeners && $("#plot").removeAllListeners("plotly_click");
    $("#plot").on("plotly_click", (ev) => {
      const p = ev.points && ev.points[0];
      if (p && p.text) { state.gene = p.text; refresh(); }
    });
  }

  // ───────── readout + ranking ─────────
  function renderReadout() {
    const r = state.gene && recOf(state.gene);
    const el = $("#ct-readout");
    if (!r) { el.innerHTML = `<div class="ct-hint">Click a point in the volcano.</div>`; return; }
    const row = (k, v, cls = "") =>
      `<div class="ct-line ${cls}"><span class="k">${k}</span><span class="v">${v}</span></div>`;
    const L = [`<div class="ct-gene">${r.g}</div>`];
    if (r.p === null) {
      L.push(`<div class="ct-warn">Only ${r.n} zygote${r.n === 1 ? "" : "s"} clears the floor —
        a correlation needs at least ${meta().params.MIN_ZYGOTES}, so this gene has no ρ.</div>`);
    } else {
      L.push(row("Spearman ρ", (r.rho >= 0 ? "+" : "") + r.rho.toFixed(3), "is-key"));
      L.push(row("P", r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4),
                 r.p < 0.05 ? "is-key" : ""));
      L.push(row("q (BH)", r.fdr == null ? "–" : (r.fdr < 1e-3 ? r.fdr.toExponential(1) : r.fdr.toFixed(3))));
      L.push(row("direction", r.rho > 0 ? "rises with τ" : "falls with τ"));
    }
    L.push(row("zygotes over the floor", `${r.n}`));
    L.push(row("median transcripts", r.medCount.toFixed(0)));
    L.push(row("probesets", `${r.nProbesets}`));
    if (r.n < 6 && r.p !== null) {
      L.push(`<div class="ct-warn">At n = ${r.n} the P is computed by <b>exhaustive
        permutation</b>, not scipy's asymptotic formula — which would report an impossible
        P = 0 for a perfect rank match.</div>`);
    }
    el.innerHTML = L.join("");
  }

  function renderRank() {
    const rows = variant().genes.filter((r) =>
      !state.find || r.g.toLowerCase().includes(state.find.toLowerCase()));
    const m = meta();
    $("#ct-rank-desc").innerHTML =
      `${testable().length} testable of ${variant().genes.length} clearing the floor.
       <b>${testable().filter((r) => r.p < 0.05).length}</b> at P&nbsp;&lt;&nbsp;0.05.`;
    $("#ct-rank").innerHTML =
      `<div class="ct-rank-head"><span>#</span><span>gene</span><span>ρ</span><span>P</span><span>n</span></div>` +
      rows.slice(0, 400).map((r, i) => `<div class="ct-rank-row${r.g === state.gene ? " on" : ""}"
          data-g="${r.g}" title="${r.g} · median ${r.medCount.toFixed(0)} transcripts · ${r.nProbesets} probeset(s)">
        <span class="n">${r.p === null ? "–" : i + 1}</span>
        <span class="e">${r.g}</span>
        <span class="r" style="color:${r.rho == null ? "#94a3b8" : r.rho > 0 ? C_UP : C_DOWN}">${
          r.rho == null ? "–" : (r.rho >= 0 ? "+" : "") + r.rho.toFixed(2)}</span>
        <span class="p">${r.p === null ? "–" : r.p < 0.001 ? r.p.toExponential(0) : r.p.toFixed(3)}</span>
        <span class="m">${r.n}</span></div>`).join("");
    $("#ct-rank").querySelectorAll(".ct-rank-row").forEach((el) =>
      el.addEventListener("click", () => { state.gene = el.dataset.g; refresh(); }));
  }

  // ───────── drawer panels ─────────
  const CFG = { displaylogo: false, responsive: true,
                modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                toImageButtonOptions: { format: "png", scale: 4 } };
  const empty = (el, msg) => { Plotly.purge(el); el.innerHTML = `<div class="ct-empty">${msg}</div>`; };
  const PS_COLOR = { "1": "#2563eb", "2": "#0f766e", "3": "#b45309", "4": "#db2777", "?": "#94a3b8" };

  function traceFor(g, showLegend) {
    const pts = (state.data.traj[g] || []).filter((p) => p.n >= meta().params.MIN_COUNT);
    const byPs = {};
    pts.forEach((p) => {
      const ps = (embOf[p.id] || {}).probeset || "?";
      (byPs[ps] = byPs[ps] || []).push(p);
    });
    return Object.entries(byPs).map(([ps, arr]) => ({
      type: "scatter", mode: "markers", name: `probeset ${ps}`,
      x: arr.map((p) => tauOf[p.id]), y: arr.map((p) => p[state.ycol]),
      marker: { size: 9, color: PS_COLOR[ps] || "#94a3b8", opacity: 0.85,
                line: { color: "#fff", width: 0.8 } },
      text: arr.map((p) => `${(embOf[p.id] || {}).label || p.id} · ${p.n} transcripts`),
      hovertemplate: "%{text}<br>τ %{x:.3f}<br>%{y:.3g}<extra></extra>",
      showlegend: !!showLegend, legendgroup: ps,
    }));
  }

  function renderOne() {
    const el = $("#ct-one");
    if (!state.gene) return empty(el, "Click a point in the volcano.");
    const r = recOf(state.gene);
    el.innerHTML = "";
    Plotly.newPlot(el, traceFor(state.gene, true), {
      margin: { l: 68, r: 16, t: 14, b: 46 }, showlegend: true,
      legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      xaxis: { title: { text: "τ — inferred pseudotime", font: { size: 10 } }, range: [-0.03, 1.03],
        gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { title: { text: YLAB(), font: { size: 10 } }, rangemode: "tozero",
        gridcolor: "#eef1f5", tickfont: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, CFG);
    $("#ct-one-sub").textContent = r
      ? `· ${state.gene} · ${r.n} zygotes` +
        (r.p === null ? " · not testable" :
         ` · ρ ${r.rho >= 0 ? "+" : ""}${r.rho.toFixed(2)}, P ${r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(3)}`)
      : "";
  }

  function renderFour() {
    const el = $("#ct-four");
    const genes = (meta().exemplars || []).filter((g) => state.data.traj[g]);
    if (!genes.length) return empty(el, "None of the reference's exemplar genes is in this dataset.");
    el.innerHTML = "";
    const traces = [];
    genes.forEach((g, i) => {
      traceFor(g, i === 0).forEach((t) => {
        t.xaxis = `x${i + 1}`; t.yaxis = `y${i + 1}`;
        traces.push(t);
      });
    });
    const lay = { margin: { l: 60, r: 16, t: 30, b: 44 },
      grid: { rows: 2, columns: 2, pattern: "independent", roworder: "top to bottom" },
      showlegend: true,
      legend: { orientation: "h", x: 0, y: 1.06, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", annotations: [] };
    genes.forEach((g, i) => {
      const r = recOf(g);
      lay[`xaxis${i + 1}`] = { title: { text: i >= 2 ? "τ" : "", font: { size: 10 } },
        range: [-0.03, 1.03], gridcolor: "#eef1f5", tickfont: { size: 9 } };
      lay[`yaxis${i + 1}`] = { rangemode: "tozero", gridcolor: "#eef1f5", tickfont: { size: 9 } };
      lay.annotations.push({ xref: `x${i + 1} domain`, yref: `y${i + 1} domain`, x: 0, y: 1.0,
        xanchor: "left", yanchor: "bottom", showarrow: false,
        text: `<b>${g}</b>` + (r && r.p !== null
          ? `  ρ ${r.rho >= 0 ? "+" : ""}${r.rho.toFixed(2)} · P ${r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(3)}`
          : "  (not testable)"),
        font: { size: 10, color: "#334155" } });
    });
    Plotly.newPlot(el, traces, lay, CFG);
    $("#ct-four-sub").textContent = `· ${YLAB()}`;
  }

  function renderTotal() {
    const el = $("#ct-total");
    el.innerHTML = "";
    const E = state.data.embryos;
    const key = state.region === "withPolar" ? "total_tx_polar" : "total_tx";
    const byPs = {};
    E.forEach((e) => { (byPs[e.probeset || "?"] = byPs[e.probeset || "?"] || []).push(e); });
    const traces = Object.entries(byPs).map(([ps, arr]) => ({
      type: "scatter", mode: "markers", name: `probeset ${ps}`,
      x: arr.map((e) => e.tau), y: arr.map((e) => e[key]),
      marker: { size: 10, color: PS_COLOR[ps] || "#94a3b8", opacity: 0.85,
                line: { color: "#fff", width: 0.8 } },
      text: arr.map((e) => e.label),
      hovertemplate: "%{text}<br>τ %{x:.3f}<br>%{y:,} transcripts<extra></extra>",
    }));
    // one trend line per probeset, because the panels differ by an order of magnitude
    Object.entries(byPs).forEach(([ps, arr]) => {
      if (arr.length < 3) return;
      const xs = arr.map((e) => e.tau), ys = arr.map((e) => e[key]);
      const mx = xs.reduce((a, b) => a + b, 0) / xs.length;
      const my = ys.reduce((a, b) => a + b, 0) / ys.length;
      let num = 0, den = 0;
      xs.forEach((x, i) => { num += (x - mx) * (ys[i] - my); den += (x - mx) ** 2; });
      if (den <= 0) return;
      const b = num / den, a = my - b * mx;
      traces.push({ type: "scatter", mode: "lines", x: [0, 1], y: [a, a + b],
        line: { color: PS_COLOR[ps] || "#94a3b8", width: 1.4, dash: "dash" },
        hoverinfo: "skip", showlegend: false });
    });
    Plotly.newPlot(el, traces, {
      margin: { l: 74, r: 16, t: 14, b: 46 }, showlegend: true,
      legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      xaxis: { title: { text: "τ — inferred pseudotime", font: { size: 10 } }, range: [-0.03, 1.03],
        gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { title: { text: "total detected transcripts", font: { size: 10 } },
        rangemode: "tozero", gridcolor: "#eef1f5", tickfont: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, CFG);
    $("#ct-total-sub").textContent =
      `· ${E.length} zygotes · dashed line = least squares within each probeset`;
  }

  function renderTable() {
    const rows = testable();
    $("#ct-table").innerHTML =
      `<table class="ct-tab"><thead><tr><th>#</th><th>gene</th><th>ρ</th><th>P</th><th>q</th>
         <th>n</th><th>median count</th><th>probesets</th></tr></thead><tbody>` +
      rows.map((r, i) => `<tr class="${r.g === state.gene ? "on" : ""}" data-g="${r.g}">
        <td>${i + 1}</td><td class="g">${r.g}</td>
        <td style="color:${r.rho > 0 ? C_UP : C_DOWN}">${(r.rho >= 0 ? "+" : "") + r.rho.toFixed(3)}</td>
        <td>${r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4)}</td>
        <td>${r.fdr == null ? "–" : r.fdr < 1e-3 ? r.fdr.toExponential(1) : r.fdr.toFixed(3)}</td>
        <td>${r.n}</td><td>${r.medCount.toFixed(0)}</td><td>${r.nProbesets}</td></tr>`).join("") +
      `</tbody></table>`;
    $("#ct-table").querySelectorAll("tr[data-g]").forEach((tr) =>
      tr.addEventListener("click", () => { state.gene = tr.dataset.g; refresh(); }));
    $("#ct-table-sub").textContent = `· ${rows.length} testable genes`;
  }

  const RENDER = { one: renderOne, four: renderFour, total: renderTotal, table: renderTable };
  function renderPanel() {
    const fn = RENDER[state.tab];
    if (!fn) return;
    try { fn(); } catch (err) { console.error("[clocktx]", state.tab, err); }
  }
  function refresh() {
    $("#drawer-gene").textContent = state.gene || "";
    renderVolcano(); renderReadout(); renderRank(); renderPanel();
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
    $("#ct-ycol").addEventListener("change", (e) => { state.ycol = e.target.value; refresh(); });
    $("#ct-region").addEventListener("change", (e) => { state.region = e.target.value; refresh(); });
    $("#t-fdr").addEventListener("change", (e) => { state.fdr = e.target.checked; renderVolcano(); });
    $("#t-labels").addEventListener("change", (e) => { state.labels = e.target.checked; renderVolcano(); });
    $("#ct-find").addEventListener("input", (e) => { state.find = e.target.value; renderRank(); });
    $("#ct-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#ct-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#ct-panels").querySelectorAll(".xs-panel").forEach((p) =>
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
                 [...$("#controls").querySelectorAll(".rz")], "clocktx_controls_box");
    window.addEventListener("resize", () => {
      ["#plot", "#ct-one", "#ct-four", "#ct-total"].forEach((s) => {
        try { Plotly.Plots.resize($(s)); } catch (_) {}
      });
    });
  }
})();
