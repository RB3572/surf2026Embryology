/* Animal–Vegetal Enrichment — figures 4.3 and 4.4.
 *
 * The plane is the EQUAL-CYTOPLASMIC-VOLUME split: the equatorial orientation slid along the polar
 * axis until both halves hold equal volume once the pronuclei and polar body are removed. Animal
 * is the polar-body side, by construction.
 *
 * ⚠️ THE HEADLINE IS A NEGATIVE RESULT, and the page says so before anything else: the number of
 * genes reaching P < 0.05 is about the number chance would give, and nothing survives BH. The
 * volcano is a ranked list, not a set of hits.
 *
 * Data: data/animalveg.json.gz (build_animalveg.py).
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const C_AN = "#b45309", C_VEG = "#0284c7", C_NS = "#cbd5e1", C_SEL = "#dc2626";

  const state = { data: null, gene: null, tab: "map", find: "" };
  const meta = () => state.data.meta;
  const recOf = (g) => state.data.genes.find((r) => r.g === g) || null;

  (async function init() {
    try { state.data = await V.loadGz("data/animalveg.json.gz"); }
    catch (err) {
      $("#plot").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_animalveg.py</code>.</div></div>`;
      return;
    }
    const m = meta();
    $("#gene-count").textContent =
      `${m.n_embryos} zygotes · ${m.n_genes} genes · equal-cytoplasmic-volume split`;
    // the negative result leads, before any gene is named
    $("#av-chance").innerHTML =
      `<b>Read this as a ranked list, not as hits.</b> ${m.n_significant} genes reach
       P&nbsp;&lt;&nbsp;${m.params.CALL_P} and chance alone would give about
       ${m.expected_by_chance}. ${allQ() ? "Nothing survives BH correction." :
       "Only what survives BH below is worth quoting."}`;
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    state.gene = state.data.genes.length ? state.data.genes[0].g : null;
    refresh();
  })();

  const allQ = () => state.data.genes.every((r) => (r.q == null ? true : r.q >= 0.05));

  // ───────── the volcano ─────────
  function renderVolcano() {
    const rows = state.data.genes;
    const col = rows.map((r) => (!r.called ? C_NS : r.lfc > 0 ? C_AN : C_VEG));
    const traces = [{
      type: "scattergl", mode: "markers",
      x: rows.map((r) => r.lfc), y: rows.map((r) => -Math.log10(Math.max(r.p, 1e-12))),
      marker: { size: rows.map((r) => 5 + 5 * Math.min(1, Math.log10(Math.max(r.total, 1)) / 3.5)),
                color: col, opacity: 0.85, line: { color: "#fff", width: 0.5 } },
      text: rows.map((r) => r.g),
      customdata: rows.map((r) => [r.n, r.total, r.p, r.q == null ? NaN : r.q]),
      hovertemplate: "<b>%{text}</b><br>mean log₂ FC %{x:.3f}<br>P %{customdata[2]:.2g} · " +
                     "q %{customdata[3]:.2g}<br>%{customdata[0]} zygotes · " +
                     "%{customdata[1]} transcripts<extra></extra>",
    }];
    const sel = state.gene && rows.find((r) => r.g === state.gene);
    if (sel) traces.push({ type: "scattergl", mode: "markers", x: [sel.lfc],
      y: [-Math.log10(Math.max(sel.p, 1e-12))],
      marker: { size: 15, color: "rgba(0,0,0,0)", line: { color: C_SEL, width: 2.5 } },
      hoverinfo: "skip", showlegend: false });
    const m = meta();
    const lim = Math.max(1.6, ...rows.map((r) => Math.abs(r.lfc))) * 1.08;
    Plotly.react($("#plot"), traces, {
      margin: { l: 62, r: 18, t: 16, b: 48 }, showlegend: false,
      xaxis: { title: { text: "← vegetal        mean log₂(animal ÷ vegetal density)        animal →",
        font: { size: 11 } }, range: [-lim, lim], zeroline: true, zerolinecolor: "#94a3b8",
        gridcolor: "#eef1f5", tickfont: { size: 10 } },
      yaxis: { title: { text: "−log₁₀ P", font: { size: 11 } }, rangemode: "tozero",
        gridcolor: "#eef1f5", tickfont: { size: 10 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [
        { type: "line", x0: -lim, x1: lim, y0: -Math.log10(m.params.CALL_P),
          y1: -Math.log10(m.params.CALL_P), line: { color: "#94a3b8", width: 1, dash: "dot" } },
        { type: "line", x0: m.params.CALL_LFC, x1: m.params.CALL_LFC, y0: 0, y1: 1, yref: "paper",
          line: { color: "#cbd5e1", width: 1, dash: "dot" } },
        { type: "line", x0: -m.params.CALL_LFC, x1: -m.params.CALL_LFC, y0: 0, y1: 1, yref: "paper",
          line: { color: "#cbd5e1", width: 1, dash: "dot" } },
      ],
      annotations: [{ x: -lim * 0.98, y: -Math.log10(m.params.CALL_P), xanchor: "left",
        yanchor: "bottom", showarrow: false, text: `P = ${m.params.CALL_P} (unadjusted)`,
        font: { size: 9, color: "#94a3b8" } }].concat(
        rows.filter((r) => r.called).slice(0, 14).map((r) => ({
          x: r.lfc, y: -Math.log10(Math.max(r.p, 1e-12)), text: r.g, showarrow: false,
          yshift: 9, font: { size: 9, color: "#334155" } }))),
    }, { displaylogo: false, responsive: true,
         modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
         toImageButtonOptions: { format: "png", scale: 4 } });
    $("#plot").on("plotly_click", (ev) => {
      const p = ev.points && ev.points[0];
      if (p && p.text) { state.gene = p.text; refresh(); }
    });
  }

  function renderReadout() {
    const r = state.gene && recOf(state.gene);
    const el = $("#av-readout");
    if (!r) { el.innerHTML = `<div class="av-hint">Click a point to select a gene.</div>`; return; }
    const row = (k, v, cls = "") =>
      `<div class="av-line ${cls}"><span class="k">${k}</span><span class="v">${v}</span></div>`;
    const L = [`<div class="av-gene">${r.g}</div>`,
      `<div class="av-side" style="background:${r.lfc > 0 ? C_AN : C_VEG}">${r.side}</div>`];
    L.push(row("mean log₂ FC", (r.lfc >= 0 ? "+" : "") + r.lfc.toFixed(3), "is-key"));
    L.push(row("across zygotes SD", r.sd.toFixed(3)));
    L.push(row("P (unadjusted)", r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4),
               r.p < 0.05 ? "is-key" : ""));
    L.push(row("q (BH)", r.q == null ? "–" : r.q.toFixed(3)));
    L.push(row("zygotes", `${r.n}`));
    L.push(row("transcripts", r.total.toLocaleString()));
    if (r.q != null && r.q >= 0.05) {
      L.push(`<div class="av-warn">Does <b>not</b> survive multiple-testing correction
        (q = ${r.q.toFixed(2)}). At ${meta().n_genes} genes tested, a P this size is not
        surprising on its own.</div>`);
    }
    el.innerHTML = L.join("");
  }

  function renderRank() {
    const rows = state.data.genes.filter((r) =>
      !state.find || r.g.toLowerCase().includes(state.find.toLowerCase()));
    $("#av-rank-desc").innerHTML =
      `${state.data.genes.length} genes tested. <b>${meta().n_called}</b> clear both calling rules
       — against about ${meta().expected_by_chance} expected by chance.`;
    $("#av-rank").innerHTML =
      `<div class="av-rank-head"><span>#</span><span>gene</span><span>log₂FC</span><span>P</span><span>n</span></div>` +
      rows.slice(0, 400).map((r) => `<div class="av-rank-row${r.g === state.gene ? " on" : ""}"
          data-g="${r.g}" title="${r.g} · ${r.total.toLocaleString()} transcripts · q ${r.q == null ? "–" : r.q.toFixed(3)}">
        <span class="n">${r.rank}</span><span class="e">${r.g}</span>
        <span class="f" style="color:${r.lfc > 0 ? C_AN : C_VEG}">${(r.lfc >= 0 ? "+" : "") + r.lfc.toFixed(2)}</span>
        <span class="p">${r.p < 0.001 ? r.p.toExponential(0) : r.p.toFixed(3)}</span>
        <span class="m">${r.n}</span></div>`).join("");
    $("#av-rank").querySelectorAll(".av-rank-row").forEach((el) =>
      el.addEventListener("click", () => { state.gene = el.dataset.g; refresh(); }));
  }

  // ───────── drawer panels ─────────
  const CFG = { displaylogo: false, responsive: true,
                modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                toImageButtonOptions: { format: "png", scale: 4 } };
  const empty = (el, msg) => { Plotly.purge(el); el.innerHTML = `<div class="av-empty">${msg}</div>`; };

  function renderMap() {
    const el = $("#av-map");
    const M = state.gene && state.data.maps[state.gene];
    if (!M) return empty(el,
      `<b>${state.gene || "This gene"}</b> has fewer than ${meta().params.MIN_MAP_COUNT} pooled ` +
      `transcripts — too few for a map that would mean anything.`);
    el.innerHTML = "";
    const NA = meta().params.NA, NR = meta().params.NR;
    const y = Array.from({ length: NA }, (_, i) => -1 + (2 * (i + 0.5)) / NA);
    const x = Array.from({ length: NR }, (_, j) => (j + 0.5) / NR);
    const lim = Math.max(0.5, ...M.z.flat().filter((v) => v != null).map(Math.abs));
    Plotly.newPlot(el, [{
      type: "heatmap", z: M.z, x, y, colorscale: "RdBu", reversescale: true,
      zmid: 0, zmin: -lim, zmax: lim,
      colorbar: { title: { text: "log₂ vs panel", font: { size: 9 } }, thickness: 11 },
      hovertemplate: "axial %{y:.2f}<br>radial %{x:.2f}<br>log₂ %{z:.2f}<extra></extra>",
    }], {
      margin: { l: 56, r: 10, t: 22, b: 46 },
      xaxis: { title: { text: "centre → cortex", font: { size: 10 } }, tickfont: { size: 9 },
        range: [0, 1] },
      yaxis: { title: { text: "vegetal → animal", font: { size: 10 } }, tickfont: { size: 9 },
        range: [-1, 1] },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [{ type: "line", x0: 0, x1: 1, y0: 0, y1: 0,
                 line: { color: "#111827", width: 1, dash: "dash" } }],
    }, CFG);
    $("#av-map-sub").textContent =
      `· ${M.n.toLocaleString()} transcripts pooled · dashed line = the equal-volume split`;
  }

  function renderPer() {
    const el = $("#av-per");
    const r = state.gene && recOf(state.gene);
    if (!r) return empty(el, "Click a gene.");
    el.innerHTML = "";
    const per = r.per.slice().sort((a, b) => a.lfc - b.lfc);
    Plotly.newPlot(el, [{
      type: "scatter", mode: "markers", x: per.map((p) => p.lfc),
      y: per.map((_, i) => i), name: "zygotes",
      marker: { size: 9, color: per.map((p) => (p.lfc > 0 ? C_AN : C_VEG)),
                line: { color: "#fff", width: 0.8 } },
      text: per.map((p) => `${p.id} · ${p.an} animal / ${p.veg} vegetal`),
      hovertemplate: "%{text}<br>log₂ FC %{x:.3f}<extra></extra>",
    }, {
      type: "scatter", mode: "lines", x: [r.lfc, r.lfc], y: [-1, per.length],
      line: { color: C_SEL, width: 2, dash: "dash" }, name: "mean", hoverinfo: "skip",
    }], {
      margin: { l: 40, r: 14, t: 14, b: 46 }, showlegend: false,
      xaxis: { title: { text: "log₂(animal ÷ vegetal density)", font: { size: 10 } },
        zeroline: true, zerolinecolor: "#111827", gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { showticklabels: false, gridcolor: "#f6f8fb", range: [-1, per.length] },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, CFG);
    $("#av-per-sub").textContent =
      `· ${r.g} · ${r.n} zygotes · mean ${r.lfc >= 0 ? "+" : ""}${r.lfc.toFixed(3)} ` +
      `· SD ${r.sd.toFixed(3)} · P ${r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4)}`;
  }

  function renderTable() {
    const rows = state.data.genes;
    $("#av-table").innerHTML =
      `<table class="av-tab"><thead><tr><th>#</th><th>gene</th><th>log₂FC</th><th>SD</th>
        <th>P</th><th>q</th><th>zygotes</th><th>transcripts</th><th>called</th></tr></thead><tbody>` +
      rows.map((r) => `<tr class="${r.g === state.gene ? "on" : ""}" data-g="${r.g}">
        <td>${r.rank}</td><td class="g">${r.g}</td>
        <td style="color:${r.lfc > 0 ? C_AN : C_VEG}">${(r.lfc >= 0 ? "+" : "") + r.lfc.toFixed(3)}</td>
        <td>${r.sd.toFixed(3)}</td>
        <td>${r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4)}</td>
        <td>${r.q == null ? "–" : r.q.toFixed(3)}</td>
        <td>${r.n}</td><td>${r.total.toLocaleString()}</td>
        <td>${r.called ? "✓" : ""}</td></tr>`).join("") + `</tbody></table>`;
    $("#av-table").querySelectorAll("tr[data-g]").forEach((tr) =>
      tr.addEventListener("click", () => { state.gene = tr.dataset.g; refresh(); }));
    $("#av-table-sub").textContent = `· ${rows.length} genes · ${meta().n_called} called`;
  }

  const RENDER = { map: renderMap, per: renderPer, table: renderTable };
  function renderPanel() {
    const fn = RENDER[state.tab];
    if (!fn) return;
    try { fn(); } catch (err) { console.error("[animalveg]", state.tab, err); }
  }
  function refresh() {
    $("#drawer-gene").textContent = state.gene || "";
    renderVolcano(); renderReadout(); renderRank(); renderPanel();
  }

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
    $("#av-find").addEventListener("input", (e) => { state.find = e.target.value; renderRank(); });
    $("#av-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#av-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#av-panels").querySelectorAll(".xs-panel").forEach((p) =>
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
                 [...$("#controls").querySelectorAll(".rz")], "animalveg_controls_box");
    window.addEventListener("resize", () => {
      ["#plot", "#av-map", "#av-per"].forEach((s) => {
        try { Plotly.Plots.resize($(s)); } catch (_) {}
      });
    });
  }
})();
