/* Contact Halves — the blastomere contact region on the reference's own definition
 * (figures 7.1, 7.2 and 7.3).
 *
 * Each blastomere is split by a plane perpendicular to its own junction→edge axis, slid until
 * THAT BLASTOMERE'S two halves hold equal volume. Contact is the junction-side half — two per
 * embryo. Nuclei and the polar body are excluded by segment label.
 *
 * ⚠️ THE HEADLINE IS THE HIT COUNT BESIDE ITS OWN CHANCE EXPECTATION, which is the whole point of
 * the reference's panel: 11 nominal against 10 expected at early 2-cell. A page headed "nothing is
 * enriched" would pre-empt the reader; 11-vs-10 lets them weigh it.
 *
 * Data: data/contacthalves.json.gz (build_contacthalves.py).
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const C_CON = "#7c3aed", C_EDGE = "#ea580c", C_NS = "#cbd5e1", C_SEL = "#dc2626";

  const state = { data: null, stage: "early2cell", gene: null, tab: "profile", find: "" };
  const meta = () => state.data.meta;
  const S = () => state.data.stages[state.stage];
  const recOf = (g) => S().genes.find((r) => r.g === g) || null;

  (async function init() {
    try { state.data = await V.loadGz("data/contacthalves.json.gz"); }
    catch (err) {
      $("#plot").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_contacthalves.py</code>.</div></div>`;
      return;
    }
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    renderValidation();
    selectStage(state.stage);
  })();

  function selectStage(key) {
    state.stage = key;
    $("#ch-stage").querySelectorAll("button").forEach((b) =>
      b.classList.toggle("on", b.dataset.stage === key));
    const s = S();
    $("#gene-count").textContent =
      `${s.embryos.length} embryos · ${s.genes.length} genes · equal-volume half of each blastomere`;
    // the count beside its own chance expectation, exactly as the reference frames it
    const excess = s.n_nominal - s.expected;
    $("#ch-chance").innerHTML =
      `<b>${s.n_nominal} nominal hits, ${s.expected} expected by chance.</b> ` +
      (s.n_fdr ? `${s.n_fdr} survive${s.n_fdr === 1 ? "s" : ""} FDR.`
               : "Nothing survives FDR.") +
      (excess <= 1.5 ? " That is the result — read the volcano as a ranked list." : "");
    state.gene = s.genes.length ? s.genes[0].g : null;
    refresh();
  }

  // ───────── the volcano (7.2) ─────────
  function renderVolcano() {
    const rows = S().genes;
    const sig = (r) => r.p < 0.05;
    const col = rows.map((r) => (!sig(r) ? C_NS : r.lfc > 0 ? C_CON : C_EDGE));
    const traces = [{
      type: "scattergl", mode: "markers",
      x: rows.map((r) => r.lfc), y: rows.map((r) => -Math.log10(Math.max(r.p, 1e-12))),
      marker: { size: rows.map((r) => 5 + 5 * Math.min(1, Math.log10(Math.max(r.total, 1)) / 5.5)),
                color: col, opacity: 0.85, line: { color: "#fff", width: 0.5 } },
      text: rows.map((r) => r.g),
      customdata: rows.map((r) => [r.n, r.total, r.p, r.q == null ? NaN : r.q]),
      hovertemplate: "<b>%{text}</b><br>mean log₂(contact ÷ edge) %{x:.3f}<br>" +
                     "P %{customdata[2]:.2g} · q %{customdata[3]:.2g}<br>" +
                     "%{customdata[0]} embryos · %{customdata[1]} transcripts<extra></extra>",
    }];
    const sel = state.gene && rows.find((r) => r.g === state.gene);
    if (sel) traces.push({ type: "scattergl", mode: "markers", x: [sel.lfc],
      y: [-Math.log10(Math.max(sel.p, 1e-12))],
      marker: { size: 15, color: "rgba(0,0,0,0)", line: { color: C_SEL, width: 2.5 } },
      hoverinfo: "skip", showlegend: false });
    const lim = Math.max(0.35, ...rows.map((r) => Math.abs(r.lfc))) * 1.12;
    Plotly.react($("#plot"), traces, {
      margin: { l: 62, r: 18, t: 16, b: 48 }, showlegend: false,
      xaxis: { title: { text: "← cell edge        mean log₂(contact ÷ edge)        contact →",
        font: { size: 11 } }, range: [-lim, lim], zeroline: true, zerolinecolor: "#94a3b8",
        gridcolor: "#eef1f5", tickfont: { size: 10 } },
      yaxis: { title: { text: "−log₁₀ P", font: { size: 11 } }, rangemode: "tozero",
        gridcolor: "#eef1f5", tickfont: { size: 10 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [{ type: "line", x0: -lim, x1: lim, y0: -Math.log10(0.05), y1: -Math.log10(0.05),
                 line: { color: "#94a3b8", width: 1, dash: "dot" } }],
      annotations: [{ x: -lim * 0.98, y: -Math.log10(0.05), xanchor: "left", yanchor: "bottom",
        showarrow: false, text: "P = 0.05 (nominal)", font: { size: 9, color: "#94a3b8" } }].concat(
        rows.filter(sig).slice(0, 6).map((r) => ({
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
    const el = $("#ch-readout");
    if (!r) { el.innerHTML = `<div class="ch-hint">Click a point to select a gene.</div>`; return; }
    const row = (k, v, cls = "") =>
      `<div class="ch-line ${cls}"><span class="k">${k}</span><span class="v">${v}</span></div>`;
    const L = [`<div class="ch-gene">${r.g}</div>`,
      `<div class="ch-side" style="background:${r.lfc > 0 ? C_CON : C_EDGE}">${r.side}</div>`];
    L.push(row("mean log₂(contact ÷ edge)", (r.lfc >= 0 ? "+" : "") + r.lfc.toFixed(3), "is-key"));
    L.push(row("across embryos SD", r.sd.toFixed(3)));
    L.push(row("P (nominal)", r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4),
               r.p < 0.05 ? "is-key" : ""));
    L.push(row("q (BH, this stage)", r.q == null ? "–" : r.q.toFixed(3)));
    L.push(row("embryos", `${r.n}`));
    L.push(row("transcripts", r.total.toLocaleString()));
    if (r.q != null && r.q >= 0.05) {
      L.push(`<div class="ch-warn">Does <b>not</b> survive FDR (q = ${r.q.toFixed(2)}). This
        stage returns ${S().n_nominal} nominal hits against ${S().expected} expected by chance, so
        a P this size is what the family produces on its own.</div>`);
    }
    el.innerHTML = L.join("");
  }

  function renderRank() {
    const all = S().genes;
    const rows = all.filter((r) =>
      !state.find || r.g.toLowerCase().includes(state.find.toLowerCase()));
    $("#ch-rank-desc").innerHTML =
      `${all.length} genes tested at this stage. <b>${S().n_nominal}</b> reach P&nbsp;&lt;&nbsp;0.05
       — against about ${S().expected} expected by chance.`;
    $("#ch-rank").innerHTML =
      `<div class="ch-rank-head"><span>#</span><span>gene</span><span>log₂FC</span><span>P</span><span>n</span></div>` +
      rows.slice(0, 400).map((r) => `<div class="ch-rank-row${r.g === state.gene ? " on" : ""}"
          data-g="${r.g}" title="${r.g} · ${r.total.toLocaleString()} transcripts · q ${r.q == null ? "–" : r.q.toFixed(3)}">
        <span class="n">${r.rank}</span><span class="e">${r.g}</span>
        <span class="f" style="color:${r.lfc > 0 ? C_CON : C_EDGE}">${(r.lfc >= 0 ? "+" : "") + r.lfc.toFixed(2)}</span>
        <span class="p">${r.p < 0.001 ? r.p.toExponential(0) : r.p.toFixed(3)}</span>
        <span class="m">${r.n}</span></div>`).join("");
    $("#ch-rank").querySelectorAll(".ch-rank-row").forEach((el) =>
      el.addEventListener("click", () => { state.gene = el.dataset.g; refresh(); }));
  }

  function renderValidation() {
    const v = meta().validation;
    const el = $("#ch-val");
    if (!v || !v.available) {
      el.innerHTML = `<b>Not validated on this machine.</b> The reference's own 7.2 table was not
        readable here, so the folds below stand on this build alone.`;
      return;
    }
    const parts = Object.entries(v.stages).map(([k, s]) =>
      `${k === "early2cell" ? "early" : "late"} 2-cell: ${s.n_shared} genes, r = ${s.r},
       median |Δ| ${s.median_abs_diff}`);
    el.innerHTML = `<b>Checked against the reference's own table.</b> Every gene it tests is tested
      here and no more — ${parts.join("; ")}. The residual is the plane itself: this build places
      it on the mesh, and small placement differences move a fold by about 0.02.
      <br><br><b>The GO dot plot (7.5) is not here.</b> It needs a gene→term annotation source that
      does not exist on this machine, and 12 imported rows whose counts cannot be re-derived would
      be a picture of someone else's computation.`;
  }

  // ───────── drawer panels ─────────
  const CFG = { displaylogo: false, responsive: true,
                modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                toImageButtonOptions: { format: "png", scale: 4 } };
  const empty = (el, msg) => { Plotly.purge(el); el.innerHTML = `<div class="ch-empty">${msg}</div>`; };

  function renderProfile() {
    const el = $("#ch-profile");
    const rows = state.data.profile.filter((p) => p.stage === state.stage);
    if (!rows.length) return empty(el, "No profile for this stage.");
    el.innerHTML = "";
    const NB = meta().params.NBIN;
    const x = Array.from({ length: NB }, (_, i) => (i + 0.5) / NB);
    const M = rows.map((r) => r.f);
    const mu = x.map((_, i) => M.reduce((a, f) => a + f[i], 0) / M.length);
    const sd = x.map((_, i) =>
      Math.sqrt(M.reduce((a, f) => a + (f[i] - mu[i]) ** 2, 0) / M.length));
    const band = (sign) => mu.map((m, i) => m + sign * sd[i]);
    Plotly.newPlot(el, [
      { type: "scatter", mode: "lines", x, y: band(+1), line: { width: 0 },
        hoverinfo: "skip", showlegend: false },
      { type: "scatter", mode: "lines", x, y: band(-1), line: { width: 0 }, fill: "tonexty",
        fillcolor: "rgba(126,136,150,0.24)", hoverinfo: "skip", name: "± s.d." },
      { type: "scatter", mode: "lines", x, y: mu, line: { color: "#334155", width: 2 },
        name: "mean", hovertemplate: "x %{x:.3f}<br>fraction %{y:.4f}<extra></extra>" },
      { type: "scatter", mode: "lines", x: [0, 1], y: [1 / NB, 1 / NB],
        line: { color: "#94a3b8", width: 1, dash: "dot" }, name: "flat (no axial structure)" },
    ], {
      margin: { l: 62, r: 14, t: 14, b: 52 }, showlegend: false,
      xaxis: { title: { text: "junction  →  equal-volume split  →  cell edge", font: { size: 10 } },
        range: [0, 1], tickvals: [0, 0.5, 1],
        ticktext: ["contact plane", "equal-volume split", "cell edge"],
        tickfont: { size: 9 }, gridcolor: "#f2f5f9" },
      yaxis: { title: { text: "fraction of transcripts", font: { size: 10 } },
        tickfont: { size: 9 }, gridcolor: "#eef1f5", rangemode: "tozero" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [{ type: "line", x0: 0.5, x1: 0.5, y0: 0, y1: 1, yref: "paper",
                 line: { color: "#6b7280", width: 1, dash: "dash" } }],
    }, CFG);
    $("#ch-profile-sub").textContent =
      `· ${rows.length} embryos · mean ± s.d. · a flat line at ${(1 / NB).toFixed(2)} would mean ` +
      `no axial structure at all`;
  }

  function renderMaps() {
    const el = $("#ch-maps");
    const maps = state.data.maps;
    if (state.stage !== "early2cell") return empty(el,
      "The reference pools its density maps over the early 2-cell embryos only, so that is what " +
      "this build ships. Switch to Early 2-cell.");
    el.innerHTML = "";
    const N = meta().params.NMAP;
    const ax = Array.from({ length: N }, (_, i) => (i + 0.5) / N);
    const all = ["contact", "edge"].flatMap((k) => maps[k].z.flat()).filter((v) => v != null);
    const s = all.map(Math.abs).sort((a, b) => a - b);
    const lim = Math.max(0.4, s[Math.floor(s.length * 0.97)] || 1);
    const traces = ["contact", "edge"].map((k, i) => ({
      type: "heatmap", z: maps[k].z[0].map((_, c) => maps[k].z.map((r) => r[c])),
      x: ax, y: ax, xaxis: i ? "x2" : "x", yaxis: i ? "y2" : "y",
      colorscale: "RdBu", reversescale: true, zmid: 0, zmin: -lim, zmax: lim,
      zsmooth: "best", showscale: i === 1,
      colorbar: { title: { text: "log₂ vs all genes", font: { size: 9 } }, thickness: 11,
                  len: 0.9 },
      hovertemplate: "junction→edge %{x:.2f}<br>out from axis %{y:.2f}<br>log₂ %{z:.2f}<extra></extra>",
    }));
    const line = (xr) => ({ type: "line", xref: xr, yref: xr === "x" ? "y" : "y2",
      x0: 0.5, x1: 0.5, y0: 0, y1: 1, line: { color: "#111827", width: 1, dash: "dash" } });
    Plotly.newPlot(el, traces, {
      margin: { l: 54, r: 14, t: 34, b: 46 }, showlegend: false,
      xaxis: { domain: [0, 0.44], anchor: "y",
        title: { text: "junction → edge", font: { size: 10 } }, tickfont: { size: 9 } },
      xaxis2: { domain: [0.52, 0.92], anchor: "y2",
        title: { text: "junction → edge", font: { size: 10 } }, tickfont: { size: 9 } },
      yaxis: { domain: [0, 1], anchor: "x",
        title: { text: "out from the axis", font: { size: 10 } }, tickfont: { size: 9 } },
      yaxis2: { domain: [0, 1], anchor: "x2", tickfont: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [line("x"), line("x2")],
      annotations: [
        { x: 0.22, y: 1.06, xref: "paper", yref: "paper", showarrow: false,
          text: `Contact-leaning · top ${meta().params.TOP_MAP}`,
          font: { size: 11, color: C_CON } },
        { x: 0.735, y: 1.06, xref: "paper", yref: "paper", showarrow: false,
          text: `Edge-leaning · top ${meta().params.TOP_MAP}`, font: { size: 11, color: C_EDGE } },
      ],
    }, CFG);
    $("#ch-maps-sub").textContent =
      `· ${maps.contact.n.toLocaleString()} + ${maps.edge.n.toLocaleString()} transcripts pooled ` +
      `· dashed line = the equal-volume split`;
  }

  function renderPer() {
    const el = $("#ch-per");
    const r = state.gene && recOf(state.gene);
    if (!r) return empty(el, "Click a gene.");
    el.innerHTML = "";
    const per = r.per.slice().sort((a, b) => a.lfc - b.lfc);
    Plotly.newPlot(el, [{
      type: "scatter", mode: "markers", x: per.map((p) => p.lfc), y: per.map((_, i) => i),
      marker: { size: 9, color: per.map((p) => (p.lfc > 0 ? C_CON : C_EDGE)),
                line: { color: "#fff", width: 0.8 } },
      text: per.map((p) => `${p.id} · ${p.c} contact / ${p.e} edge`),
      hovertemplate: "%{text}<br>log₂ FC %{x:.3f}<extra></extra>",
    }, {
      type: "scatter", mode: "lines", x: [r.lfc, r.lfc], y: [-1, per.length],
      line: { color: C_SEL, width: 2, dash: "dash" }, hoverinfo: "skip",
    }], {
      margin: { l: 40, r: 14, t: 14, b: 46 }, showlegend: false,
      xaxis: { title: { text: "log₂(contact ÷ edge), bulk-centred", font: { size: 10 } },
        zeroline: true, zerolinecolor: "#111827", gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { showticklabels: false, gridcolor: "#f6f8fb", range: [-1, per.length] },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, CFG);
    $("#ch-per-sub").textContent =
      `· ${r.g} · ${r.n} embryos · mean ${r.lfc >= 0 ? "+" : ""}${r.lfc.toFixed(3)} ` +
      `· SD ${r.sd.toFixed(3)} · P ${r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4)}`;
  }

  function renderTable() {
    const rows = S().genes;
    $("#ch-table").innerHTML =
      `<table class="ch-tab"><thead><tr><th>#</th><th>gene</th><th>log₂FC</th><th>SD</th>
        <th>P</th><th>q</th><th>embryos</th><th>transcripts</th><th>side</th></tr></thead><tbody>` +
      rows.map((r) => `<tr class="${r.g === state.gene ? "on" : ""}" data-g="${r.g}">
        <td>${r.rank}</td><td class="g">${r.g}</td>
        <td style="color:${r.lfc > 0 ? C_CON : C_EDGE}">${(r.lfc >= 0 ? "+" : "") + r.lfc.toFixed(3)}</td>
        <td>${r.sd.toFixed(3)}</td>
        <td>${r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4)}</td>
        <td>${r.q == null ? "–" : r.q.toFixed(3)}</td>
        <td>${r.n}</td><td>${r.total.toLocaleString()}</td>
        <td>${r.p < 0.05 ? r.side : ""}</td></tr>`).join("") + `</tbody></table>`;
    $("#ch-table").querySelectorAll("tr[data-g]").forEach((tr) =>
      tr.addEventListener("click", () => { state.gene = tr.dataset.g; refresh(); }));
    $("#ch-table-note").textContent = S().n_fdr
      ? `${S().n_fdr} gene${S().n_fdr === 1 ? "" : "s"} survive${S().n_fdr === 1 ? "s" : ""} it.`
      : "Nothing survives it at this stage.";
    $("#ch-table-sub").textContent =
      `· ${rows.length} genes · ${S().n_nominal} nominal · ${S().expected} expected by chance`;
  }

  const RENDER = { profile: renderProfile, maps: renderMaps, per: renderPer, table: renderTable };
  function renderPanel() {
    const fn = RENDER[state.tab];
    if (!fn) return;
    try { fn(); } catch (err) { console.error("[contacthalves]", state.tab, err); }
    // the panel is drawn while the drawer is still animating open, so the container can still be
    // narrower than its final width; one frame later it is not, and Plotly needs telling
    const kick = () => {
      const el = $("#ch-" + state.tab);
      if (el && el.querySelector(".main-svg")) { try { Plotly.Plots.resize(el); } catch (_) {} }
    };
    requestAnimationFrame(kick);
    setTimeout(kick, 160);
    setTimeout(kick, 420);
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
    $("#ch-stage").addEventListener("click", (ev) => {
      const b = ev.target.closest("button[data-stage]");
      if (b && b.dataset.stage !== state.stage) selectStage(b.dataset.stage);
    });
    $("#ch-find").addEventListener("input", (e) => { state.find = e.target.value; renderRank(); });
    $("#ch-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#ch-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#ch-panels").querySelectorAll(".xs-panel").forEach((p) =>
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
                 [...$("#controls").querySelectorAll(".rz")], "contacthalves_controls_box");
    window.addEventListener("resize", () => {
      ["#plot", "#ch-profile", "#ch-maps", "#ch-per"].forEach((s) => {
        try { Plotly.Plots.resize($(s)); } catch (_) {}
      });
    });
  }
})();
