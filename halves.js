/* Halves — the two halves of a zygote, four ways of cutting it
 * (figures 4.14, 4.15, 4.17, 4.18 and 4.19).
 *
 * The one decision everything hangs off is WHICH SIDE IS WHICH, and neither obvious landmark works:
 * the sperm lies ON the plane it defines, and the pronuclei are not independent of the sperm. So
 * the fuller half — more cytoplasmic transcripts over the whole panel — names the side. It is
 * intrinsic, has no free sign, and is the same quantity on every plane definition.
 *
 * ⚠️ THE PAGE LEADS WITH THE ALIGNMENT NULL, not with the volcano. Any consistent labelling of
 * halves produces some calls; the question is whether the real one produces more.
 *
 * Data: data/halves.json.gz (build_halves.py).
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const C_F = "#059669", C_E = "#b45309", C_NS = "#cbd5e1", C_SEL = "#dc2626";

  const state = { data: null, plane: "polar18", gene: null, tab: "null", find: "", norm: "conc" };
  const meta = () => state.data.meta;
  const S = () => state.data.volcano[state.plane];
  const recOf = (g) => S().genes.find((r) => r.g === g) || null;

  (async function init() {
    try { state.data = await V.loadGz("data/halves.json.gz"); }
    catch (err) {
      $("#plot").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_halves.py</code>.</div></div>`;
      return;
    }
    $("#hv-plane").innerHTML = Object.keys(state.data.volcano).map((k) =>
      `<button type="button" data-plane="${k}"${k === state.plane ? ' class="on"' : ""}>${state.data.volcano[k].label}</button>`).join("");
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    renderSidebarNote();
    selectPlane(state.plane);
  })();

  function selectPlane(k) {
    state.plane = k;
    $("#hv-plane").querySelectorAll("button").forEach((b) =>
      b.classList.toggle("on", b.dataset.plane === k));
    const s = S();
    $("#gene-count").textContent =
      `${s.n_embryos} zygotes · ${s.genes.length} genes · fuller half vs emptier half`;
    // the alignment null, not the call count, is the headline
    const n = s.null;
    $("#hv-chance").innerHTML =
      `<b>${s.n_called} genes called; flipping the halves at random gives ${Math.round(n.median)}.</b> ` +
      (n.p < 0.05
        ? `Beyond the null (P = ${n.p.toFixed(3)}).`
        : `That is inside the null (P = ${n.p.toFixed(3)}) — read this as a ranked list.`);
    state.gene = s.genes.length ? s.genes[0].g : null;
    refresh();
  }

  // ───────── the volcano (4.19) ─────────
  function renderVolcano() {
    const rows = S().genes;
    const col = rows.map((r) => (!r.called ? C_NS : r.lfc > 0 ? C_F : C_E));
    const traces = [{
      type: "scattergl", mode: "markers",
      x: rows.map((r) => r.lfc), y: rows.map((r) => -Math.log10(Math.max(r.p, 1e-12))),
      marker: { size: rows.map((r) => 5 + 5 * Math.min(1, Math.log10(Math.max(r.total, 1)) / 5.5)),
                color: col, opacity: 0.85, line: { color: "#fff", width: 0.5 } },
      text: rows.map((r) => r.g),
      customdata: rows.map((r) => [r.n, r.total, r.p, r.q == null ? NaN : r.q]),
      hovertemplate: "<b>%{text}</b><br>mean log₂(fuller ÷ emptier) %{x:.3f}<br>" +
                     "P %{customdata[2]:.2g} · q %{customdata[3]:.2g}<br>" +
                     "%{customdata[0]} zygotes · %{customdata[1]} transcripts<extra></extra>",
    }];
    const sel = state.gene && rows.find((r) => r.g === state.gene);
    if (sel) traces.push({ type: "scattergl", mode: "markers", x: [sel.lfc],
      y: [-Math.log10(Math.max(sel.p, 1e-12))],
      marker: { size: 15, color: "rgba(0,0,0,0)", line: { color: C_SEL, width: 2.5 } },
      hoverinfo: "skip", showlegend: false });
    const lim = Math.max(0.4, ...rows.map((r) => Math.abs(r.lfc))) * 1.12;
    Plotly.react($("#plot"), traces, {
      margin: { l: 62, r: 18, t: 16, b: 48 }, showlegend: false,
      xaxis: { title: { text: "← emptier half        mean log₂(fuller ÷ emptier concentration)        fuller half →",
        font: { size: 11 } }, range: [-lim, lim], zeroline: true, zerolinecolor: "#94a3b8",
        gridcolor: "#eef1f5", tickfont: { size: 10 } },
      yaxis: { title: { text: "−log₁₀ P", font: { size: 11 } }, rangemode: "tozero",
        gridcolor: "#eef1f5", tickfont: { size: 10 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [{ type: "line", x0: -lim, x1: lim, y0: -Math.log10(0.05), y1: -Math.log10(0.05),
                 line: { color: "#94a3b8", width: 1, dash: "dot" } }],
      annotations: [{ x: -lim * 0.98, y: -Math.log10(0.05), xanchor: "left", yanchor: "bottom",
        showarrow: false, text: "P = 0.05 (nominal)", font: { size: 9, color: "#94a3b8" } },
        { x: 0, y: 1, xref: "x", yref: "paper", yanchor: "top", yshift: -4, showarrow: false,
          text: "the bulk correction puts a gene that merely tracks total density here",
          font: { size: 9, color: "#b0b8c4" } }].concat(
        rows.filter((r) => r.called).slice(0, 12).map((r) => ({
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
    const el = $("#hv-readout");
    if (!r) { el.innerHTML = `<div class="hv-hint">Click a point to select a gene.</div>`; return; }
    const row = (k, v, cls = "") =>
      `<div class="hv-line ${cls}"><span class="k">${k}</span><span class="v">${v}</span></div>`;
    const L = [`<div class="hv-gene">${r.g}</div>`,
      `<div class="hv-side" style="background:${r.lfc > 0 ? C_F : C_E}">${r.side} half</div>`];
    L.push(row("mean log₂ FC", (r.lfc >= 0 ? "+" : "") + r.lfc.toFixed(3), "is-key"));
    L.push(row("across zygotes SD", r.sd.toFixed(3)));
    L.push(row("P (nominal)", r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4),
               r.p < 0.05 ? "is-key" : ""));
    L.push(row("q (BH)", r.q == null ? "–" : r.q.toFixed(3)));
    L.push(row("zygotes", `${r.n}`));
    L.push(row("transcripts", r.total.toLocaleString()));
    const pg = (state.data.perGene.rows || []).find((x) => x.g === r.g);
    if (pg && pg.excess != null) {
      L.push(row("fold vs count-matched null", `${pg.fold.toFixed(3)} vs ${pg.null.toFixed(3)}`));
      L.push(row("excess over the null", (pg.excess >= 0 ? "+" : "") + pg.excess.toFixed(3),
                 pg.excess > 0 ? "" : "is-warn"));
    }
    if (S().null.p >= 0.05) {
      L.push(`<div class="hv-warn">This plane's whole call set is inside the alignment null
        (P = ${S().null.p.toFixed(3)}), so a single gene's P is not evidence on its own.</div>`);
    }
    el.innerHTML = L.join("");
  }

  function renderRank() {
    const all = S().genes;
    const rows = all.filter((r) =>
      !state.find || r.g.toLowerCase().includes(state.find.toLowerCase()));
    $("#hv-rank-desc").innerHTML =
      `${all.length} genes on this plane. <b>${S().n_called}</b> reach P&nbsp;&lt;&nbsp;0.05 —
       against ${Math.round(S().null.median)} from randomly flipped halves.`;
    $("#hv-rank").innerHTML =
      `<div class="hv-rank-head"><span>#</span><span>gene</span><span>log₂FC</span><span>P</span><span>n</span></div>` +
      rows.slice(0, 400).map((r) => `<div class="hv-rank-row${r.g === state.gene ? " on" : ""}"
          data-g="${r.g}" title="${r.g} · ${r.total.toLocaleString()} transcripts · q ${r.q == null ? "–" : r.q.toFixed(3)}">
        <span class="n">${r.rank}</span><span class="e">${r.g}</span>
        <span class="f" style="color:${r.lfc > 0 ? C_F : C_E}">${(r.lfc >= 0 ? "+" : "") + r.lfc.toFixed(2)}</span>
        <span class="p">${r.p < 0.001 ? r.p.toExponential(0) : r.p.toFixed(3)}</span>
        <span class="m">${r.n}</span></div>`).join("");
    $("#hv-rank").querySelectorAll(".hv-rank-row").forEach((el) =>
      el.addEventListener("click", () => { state.gene = el.dataset.g; refresh(); }));
  }

  function renderSidebarNote() {
    $("#hv-val").innerHTML =
      `<b>Which side is which is not free.</b> The sperm lies <i>on</i> the plane it defines — it is
       one of the three points that draw it — so it cannot name a side, and the pronuclei descend
       from the sperm, so they are not independent either. The fuller half (more cytoplasmic
       transcripts over the whole panel) is the only rule left that is intrinsic and has no free
       sign.<br><br>
       <b>A consequence worth keeping in view.</b> The bulk correction subtracts each embryo's
       median per-gene log ratio, which is positive on the fuller half by construction — so a gene
       that merely tracks total transcript density lands near zero, and the abundant genes that fix
       the orientation sit there too. What survives is enrichment <i>beyond</i> the bulk.`;
  }

  // ───────── drawer panels ─────────
  const CFG = { displaylogo: false, responsive: true,
                modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                toImageButtonOptions: { format: "png", scale: 4 } };
  const empty = (el, msg) => { Plotly.purge(el); el.innerHTML = `<div class="hv-empty">${msg}</div>`; };

  function renderNull() {
    const el = $("#hv-null");
    const s = S(), n = s.null;
    el.innerHTML = "";
    const x = n.hist.map((_, i) => i);
    Plotly.newPlot(el, [
      { type: "bar", x, y: n.hist, marker: { color: "#cbd5e1" }, name: "flipped at random",
        hovertemplate: "%{y} of " + n.draws + " draws called %{x} genes<extra></extra>" },
      { type: "scatter", mode: "lines", x: [s.n_called, s.n_called],
        y: [0, Math.max(...n.hist) * 1.05], line: { color: C_SEL, width: 3 },
        hovertemplate: `the real halves: ${s.n_called} called<extra></extra>`, name: "observed" },
    ], {
      margin: { l: 52, r: 14, t: 14, b: 46 }, showlegend: false, bargap: 0.05,
      xaxis: { title: { text: "genes called at P < 0.05", font: { size: 10 } },
        tickfont: { size: 9 }, gridcolor: "#f2f5f9" },
      yaxis: { title: { text: `draws (of ${n.draws})`, font: { size: 10 } },
        tickfont: { size: 9 }, gridcolor: "#eef1f5" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      annotations: [{ x: s.n_called, y: Math.max(...n.hist) * 1.05, yanchor: "bottom",
        showarrow: false, text: `observed ${s.n_called}`, font: { size: 10, color: C_SEL } }],
    }, CFG);
    $("#hv-null-sub").textContent =
      `· ${n.draws} draws · null median ${Math.round(n.median)}, 95th ${n.p95.toFixed(0)}, ` +
      `max ${n.max} · P = ${n.p.toFixed(4)}`;
  }

  function renderPlanePanel() {
    const el = $("#hv-planeplot");
    const rows = state.data.byPlane;
    el.innerHTML = "";
    Plotly.newPlot(el, [{
      type: "bar", orientation: "h",
      y: rows.map((r) => r.label), x: rows.map((r) => r.median),
      error_x: { type: "data", symmetric: false,
        array: rows.map((r) => r.ci_hi - r.median), arrayminus: rows.map((r) => r.median - r.ci_lo),
        color: "#64748b", thickness: 1.4, width: 5 },
      marker: { color: rows.map((r) => (r.plane === "random" ? "#94a3b8" : C_F)) },
      text: rows.map((r) => r.median.toFixed(3)), textposition: "outside",
      hovertemplate: "%{y}<br>median fold %{x:.4f}<extra></extra>",
    }], {
      margin: { l: 150, r: 60, t: 12, b: 44 }, showlegend: false,
      xaxis: { title: { text: "median fold asymmetry over the shared genes", font: { size: 10 } },
        tickfont: { size: 9 }, gridcolor: "#eef1f5",
        range: [1, Math.max(...rows.map((r) => r.ci_hi)) * 1.02] },
      yaxis: { tickfont: { size: 10 }, automargin: true },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [{ type: "line", x0: 1, x1: 1, y0: 0, y1: 1, yref: "paper",
                 line: { color: "#111827", width: 1, dash: "dot" } }],
    }, CFG);
    const R = rows.find((r) => r.plane === "random");
    $("#hv-plane-sub").textContent =
      `· ${rows[0].n_genes} genes · bootstrap 95% intervals · the random plane already sits at ` +
      `${R ? R.median.toFixed(3) : "?"}, which is counting noise, not asymmetry`;

    const H = state.data.heat;
    const best = (row) => {
      let bk = null, bv = -Infinity;
      H.planes.forEach((p) => { if (row[p] != null && row[p] > bv) { bv = row[p]; bk = p; } });
      return bk;
    };
    $("#hv-heat").innerHTML =
      `<table><thead><tr><th>gene</th>` +
      H.planes.map((p) => `<th>${state.data.meta.planes[p]}</th>`).join("") +
      `</tr></thead><tbody>` +
      H.rows.map((row) => {
        const b = best(row);
        return `<tr><td>${row.g}</td>` + H.planes.map((p) =>
          `<td class="${p === b ? "best" : ""}">${row[p] == null ? "–" : row[p].toFixed(3)}</td>`).join("") +
          `</tr>`;
      }).join("") + `</tbody></table>`;
  }

  function renderGenePanel() {
    const el = $("#hv-gene");
    const rows = (state.data.perGene.rows || []).filter((r) => r.excess != null).slice(0, 30);
    if (!rows.length) return empty(el, "No gene clears the floors on this plane.");
    el.innerHTML = "";
    Plotly.newPlot(el, [
      { type: "bar", x: rows.map((r) => r.g), y: rows.map((r) => r.null),
        marker: { color: "#e2e8f0" }, name: "count-matched null",
        hovertemplate: "%{x}<br>null fold %{y:.3f}<extra></extra>" },
      { type: "bar", x: rows.map((r) => r.g), y: rows.map((r) => r.excess),
        marker: { color: rows.map((r) => (r.g === state.gene ? C_SEL : C_F)) },
        name: "excess over the null",
        customdata: rows.map((r) => [r.fold, r.total, r.n]),
        hovertemplate: "%{x}<br>observed fold %{customdata[0]:.3f}<br>excess %{y:.3f}<br>" +
                       "%{customdata[1]} transcripts in %{customdata[2]} zygotes<extra></extra>" },
    ], {
      margin: { l: 52, r: 14, t: 12, b: 78 }, showlegend: false, barmode: "stack",
      xaxis: { tickfont: { size: 9 }, tickangle: -60 },
      yaxis: { title: { text: "fold asymmetry", font: { size: 10 } }, tickfont: { size: 9 },
        gridcolor: "#eef1f5", rangemode: "tozero" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      shapes: [{ type: "line", x0: -0.5, x1: rows.length - 0.5, y0: 1, y1: 1,
                 line: { color: "#111827", width: 1, dash: "dot" } }],
    }, CFG);
    el.on("plotly_click", (ev) => {
      const p = ev.points && ev.points[0];
      if (p && p.x) { state.gene = p.x; refresh(); }
    });
    $("#hv-gene-sub").textContent =
      `· ${state.data.meta.planes[state.data.perGene.plane]} · top ${rows.length} by excess · ` +
      `grey is the null the observed fold has to beat`;
  }

  function renderPair() {
    const el = $("#hv-pair");
    const panels = state.data.pairing[state.norm] || {};
    const keys = Object.keys(panels).sort();
    $("#hv-pairctl").querySelectorAll("[data-norm]").forEach((b) =>
      b.classList.toggle("on", b.dataset.norm === state.norm));
    $("#hv-embed-note").innerHTML =
      ` The coordinates are <b>PCA</b>, not UMAP — umap-learn is not installed here — but the
        pairing statistic is computed on the <b>full feature vectors</b>, never on these two
        coordinates, so it does not depend on the picture.`;
    if (!keys.length) return empty(el, "No probeset panel has enough zygotes.");
    el.innerHTML = "";
    const traces = [];
    const nx = keys.length;
    keys.forEach((k, i) => {
      const P = panels[k], em = P.embed;
      const ax = i ? `x${i + 1}` : "x", ay = i ? `y${i + 1}` : "y";
      em.ids.forEach((id, j) => {
        const a = em.xy[2 * j], b = em.xy[2 * j + 1];
        traces.push({ type: "scatter", mode: "lines+markers", xaxis: ax, yaxis: ay,
          x: [a[0], b[0]], y: [a[1], b[1]],
          line: { color: "rgba(100,116,139,0.35)", width: 1 },
          marker: { size: 7, color: [C_F, C_E], line: { color: "#fff", width: 0.8 } },
          text: [id, id], hovertemplate: "%{text}<extra></extra>", showlegend: false });
      });
    });
    const layout = { margin: { l: 40, r: 14, t: 52, b: 36 }, showlegend: false,
      paper_bgcolor: "transparent", plot_bgcolor: "transparent", annotations: [] };
    keys.forEach((k, i) => {
      const lo = i / nx + 0.02, hi = (i + 1) / nx - 0.02;
      layout[i ? `xaxis${i + 1}` : "xaxis"] = { domain: [lo, hi], anchor: i ? `y${i + 1}` : "y",
        showticklabels: false, gridcolor: "#f4f6f9", zeroline: false };
      layout[i ? `yaxis${i + 1}` : "yaxis"] = { domain: [0, 1], anchor: i ? `x${i + 1}` : "x",
        showticklabels: false, gridcolor: "#f4f6f9", zeroline: false };
      const st = panels[k].stat;
      // quote the tail that matches what actually happened: below 1 the halves cluster, above 1
      // they are pushed apart — reporting the clustering tail on a ratio of 1.3 would read as
      // "P = 1, nothing here" when the real finding is the artefact itself
      const below = st.ratio < 1;
      const pv = below ? st.p_closer : st.p_farther;
      layout.annotations.push({ x: (lo + hi) / 2, y: 1.10, xref: "paper", yref: "paper",
        showarrow: false, font: { size: 10, color: "#334155" },
        text: `panel ${k} · n ${st.n} · <b>${st.ratio.toFixed(2)}</b>` +
              `<br>${below ? "closer" : "farther"}, P ${pv.toFixed(3)}` });
    });
    Plotly.newPlot(el, traces, layout, CFG);
    const rs = keys.map((k) => panels[k].stat.ratio);
    $("#hv-pair-sub").textContent =
      `· ${keys.length} probeset panels · pairing ${Math.min(...rs).toFixed(2)}–${Math.max(...rs).toFixed(2)} ` +
      (state.norm === "ratio"
        ? "· above 1 is the complementarity artefact, not a result"
        : "· below 1: a zygote's two halves really do sit together, because both carry that "
          + "embryo's own abundance profile");
  }

  function renderPer() {
    const el = $("#hv-per");
    const r = state.gene && recOf(state.gene);
    if (!r) return empty(el, "Click a gene.");
    el.innerHTML = "";
    const per = r.per.slice().sort((a, b) => a.lfc - b.lfc);
    Plotly.newPlot(el, [{
      type: "scatter", mode: "markers", x: per.map((p) => p.lfc), y: per.map((_, i) => i),
      marker: { size: 9, color: per.map((p) => (p.lfc > 0 ? C_F : C_E)),
                line: { color: "#fff", width: 0.8 } },
      text: per.map((p) => `${p.id} · ${p.n} transcripts`),
      hovertemplate: "%{text}<br>log₂ FC %{x:.3f}<extra></extra>",
    }, {
      type: "scatter", mode: "lines", x: [r.lfc, r.lfc], y: [-1, per.length],
      line: { color: C_SEL, width: 2, dash: "dash" }, hoverinfo: "skip",
    }], {
      margin: { l: 40, r: 14, t: 14, b: 46 }, showlegend: false,
      xaxis: { title: { text: "log₂(fuller ÷ emptier concentration), bulk-centred", font: { size: 10 } },
        zeroline: true, zerolinecolor: "#111827", gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { showticklabels: false, gridcolor: "#f6f8fb", range: [-1, per.length] },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    }, CFG);
    $("#hv-per-sub").textContent =
      `· ${r.g} · ${r.n} zygotes · mean ${r.lfc >= 0 ? "+" : ""}${r.lfc.toFixed(3)} ` +
      `· SD ${r.sd.toFixed(3)} · P ${r.p < 1e-3 ? r.p.toExponential(1) : r.p.toFixed(4)}`;

    const O = S().orientation.slice().sort((a, b) => b.frac - a.frac);
    $("#hv-orient").innerHTML =
      `<table><thead><tr><th>zygote</th><th>fuller half</th><th>emptier half</th>
        <th>fraction in the fuller half</th><th>bulk log₂</th><th>stored side flipped?</th></tr></thead><tbody>` +
      O.map((o) => `<tr><td>${o.id}</td><td>${o.totF.toLocaleString()}</td>
        <td>${o.totE.toLocaleString()}</td><td>${(o.frac * 100).toFixed(1)}%</td>
        <td>${o.bulk >= 0 ? "+" : ""}${o.bulk.toFixed(3)}</td>
        <td>${o.flipped ? "yes" : ""}</td></tr>`).join("") + `</tbody></table>`;
  }

  const RENDER = { null: renderNull, plane: renderPlanePanel, gene: renderGenePanel,
                   pair: renderPair, per: renderPer };
  function renderPanel() {
    const fn = RENDER[state.tab];
    if (!fn) return;
    try { fn(); } catch (err) { console.error("[halves]", state.tab, err); }
    // the panel is drawn while the drawer is still animating open, so the container can still be
    // narrower than its final width; a couple of frames later it is not
    const kick = () => {
      const el = $("#hv-" + (state.tab === "plane" ? "planeplot" : state.tab));
      if (el && el.querySelector(".main-svg")) { try { Plotly.Plots.resize(el); } catch (_) {} }
    };
    requestAnimationFrame(kick); setTimeout(kick, 160); setTimeout(kick, 420);
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
    $("#hv-plane").addEventListener("click", (ev) => {
      const b = ev.target.closest("[data-plane]");
      if (b && b.dataset.plane !== state.plane) selectPlane(b.dataset.plane);
    });
    $("#hv-pairctl").addEventListener("click", (ev) => {
      const b = ev.target.closest("[data-norm]"); if (!b) return;
      state.norm = b.dataset.norm; renderPair();
    });
    $("#hv-find").addEventListener("input", (e) => { state.find = e.target.value; renderRank(); });
    $("#hv-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#hv-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#hv-panels").querySelectorAll(".xs-panel").forEach((p) =>
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
                 [...$("#controls").querySelectorAll(".rz")], "halves_controls_box");
    window.addEventListener("resize", () => {
      ["#plot", "#hv-null", "#hv-planeplot", "#hv-gene", "#hv-pair", "#hv-per"].forEach((s) => {
        try { Plotly.Plots.resize($(s)); } catch (_) {}
      });
    });
  }
})();
