/* Spatial Gene Clustering.
 *
 * One dot per gene, placed by WHERE its transcripts sit inside the zygote rather than by how
 * much of it there is. Data: data/clustering.json.gz (build_clustering.py) ships each gene's
 * mean 12-D signature (6 equal-volume radial shells + 6 axial bins, as log2 enrichment over
 * each embryo's own transcript distribution), cluster labels for every k, and two 2-D layouts.
 *
 * Colour carries exactly one message. Every cluster is a grey; the selected gene's cluster
 * darkens so you can see who it travels with; the selected gene itself is the only vermilion
 * dot on the page. Distances for the "nearest genes" panel are computed here from the
 * signatures — the same correlation distance the clustering used — so no matrix need ship.
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const SEL = "#dc2626";                       // the selected gene — the one vibrant colour
  const SEL_DIM = "#f3b0b0";

  const state = {
    data: null, genes: [], byGene: {}, k: 8, layout: "mds",
    gene: null, labels: [], clusters: [], labels_: null, tab: "sig",
    drawerOpen: false, rdrawerOpen: false, showLabels: false,
  };

  /* Cluster greys. Spread across a light band so no cluster shouts; the SELECTED cluster is
   * drawn near-black instead. Cluster order is arbitrary, so walking the ramp straight through
   * would imply a ranking the data does not have — instead adjacent cluster indices are pushed
   * to opposite halves of the ramp. The permutation must be a BIJECTION: a modulo-based
   * interleave collides, and two clusters sharing a grey makes them indistinguishable. */
  function GREY(i, n) {
    n = Math.max(n, 1);
    // even indices fill the light half in order, odd indices fill the dark half — every
    // cluster lands on its own slot, and neighbours are ~half the ramp apart
    const slot = i % 2 === 0 ? i / 2 : Math.ceil(n / 2) + (i - 1) / 2;
    const t = n <= 1 ? 0.5 : Math.min(slot, n - 1) / (n - 1);
    const l = Math.round(214 - t * 74);        // #d6d6d3 → #8c8c88
    return `rgb(${l},${l},${l - 3})`;
  }
  const DARK = "#2b2b28";                      // the selected gene's cluster
  const DARK_SOFT = "#5a5a55";

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await V.loadGz("data/clustering.json.gz"); }
    catch (e) {
      $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_clustering.py</code>.</div></div>`;
      return;
    }
    const d = state.data, m = d.meta;
    state.genes = d.genes.map((g) => g.gene);
    d.genes.forEach((g, i) => (state.byGene[g.gene] = Object.assign({ i }, g)));
    state.k = m.default_k;

    $("#embryo-count").textContent =
      `${m.n_genes} genes · ${m.n_embryos} zygotes · clustered by where they sit`;
    const ksel = $("#k-range");
    ksel.min = m.k_range[0]; ksel.max = m.k_range[m.k_range.length - 1]; ksel.value = state.k;

    $("#gene-select").innerHTML = state.genes
      .map((g) => `<option value="${g}">${g}</option>`).join("");

    applyK(state.k);
    setGene(bestSeedGene());
    wire();
    renderMethod();
    $("#controls").hidden = false; $("#placeholder").hidden = true;
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    if (!state.drawerOpen) openDrawer(true);
  })();

  /** Open on a gene with a strong, well-measured signature so the first view shows the idea. */
  function bestSeedGene() {
    let best = state.genes[0], score = -1;
    for (const g of state.genes) {
      const r = state.byGene[g];
      const amp = Math.max(...r.profile.map(Math.abs));
      const s = amp * Math.log10(1 + r.n_emb);
      if (s > score) { score = s; best = g; }
    }
    return best;
  }

  function applyK(k) {
    state.k = k;
    const kd = state.data.k[String(k)];
    state.labels = kd.labels;
    state.clusters = kd.clusters;
    $("#k-val").textContent = k;
    const sil = state.data.meta.silhouette[String(k)];
    const best = state.data.meta.default_k;
    $("#k-hint").innerHTML =
      `silhouette ${sil >= 0 ? "+" : ""}${sil.toFixed(3)}` +
      (k === best ? " · best of the range" : ` · best is k=${best}`) +
      `<br>sizes ${state.clusters.map((c) => c.n).join(" · ")}`;
  }

  const clusterOf = (g) => state.labels[state.byGene[g].i];

  function setGene(g) {
    if (!(g in state.byGene)) return;
    state.gene = g;
    $("#gene-select").value = g;
    renderMap(); renderReadout(); renderMembers(); renderActive();
    $("#drawer-gene").textContent = `· ${g}`;
  }

  // ───────── the map ─────────
  function renderMap() {
    const d = state.data, key = state.layout;
    const sel = state.gene, selC = clusterOf(sel);
    const n = state.clusters.length;

    // One trace per cluster keeps the legend meaningful and lets the selected cluster sit on top.
    const traces = [];
    for (let c = 0; c < n; c++) {
      const idx = state.labels.map((l, i) => (l === c ? i : -1)).filter((i) => i >= 0);
      if (!idx.length) continue;
      const isSel = c === selC;
      const pts = idx.map((i) => d.genes[i]);
      traces.push({
        type: "scattergl", mode: state.showLabels && isSel ? "markers+text" : "markers",
        name: `Cluster ${c + 1}${isSel ? " ·  selected" : ""} (${idx.length})`,
        x: pts.map((p) => p[key][0]), y: pts.map((p) => p[key][1]),
        text: pts.map((p) => p.gene),
        textposition: "top center",
        textfont: { size: 9, color: DARK_SOFT },
        marker: {
          size: isSel ? 9 : 7,
          color: isSel ? DARK : GREY(c, n),
          line: { width: isSel ? 1 : 0.5, color: isSel ? "#fff" : "rgba(255,255,255,.75)" },
          opacity: isSel ? 1 : 0.9,
        },
        customdata: pts.map((p) => [p.gene, c + 1, p.n_emb]),
        hovertemplate: "<b>%{customdata[0]}</b><br>cluster %{customdata[1]} · %{customdata[2]} embryos<extra></extra>",
        legendrank: isSel ? 0 : 10 + c,
      });
    }
    // the selected gene, last so it draws above everything
    const p = state.byGene[sel];
    traces.push({
      type: "scattergl", mode: "markers+text", name: sel,
      x: [p[key][0]], y: [p[key][1]], text: [sel], textposition: "top center",
      textfont: { size: 12, color: SEL, family: "system-ui, sans-serif" },
      marker: { size: 15, color: SEL, line: { width: 2, color: "#fff" }, symbol: "circle" },
      hovertemplate: `<b>${sel}</b> · selected<extra></extra>`, legendrank: -1,
    });

    const lay = {
      margin: { l: 10, r: 10, t: 10, b: 10 },
      paper_bgcolor: "transparent", plot_bgcolor: "#fcfdfe",
      xaxis: { visible: false, zeroline: false, scaleanchor: "y", scaleratio: 1 },
      yaxis: { visible: false, zeroline: false },
      showlegend: true,
      legend: { x: 1, xanchor: "right", y: 1, font: { size: 10 }, bgcolor: "rgba(255,255,255,.82)",
                bordercolor: "#e7e9ef", borderwidth: 1 },
      hoverlabel: { bgcolor: "#fff", bordercolor: "#e7e9ef", font: { size: 11 } },
      font: { size: 11, color: "#334155" },
      uirevision: key,      // keep pan/zoom when only the selection changes
    };
    Plotly.react($("#plot-host"), traces, lay, { responsive: true, displaylogo: false, displayModeBar: false });

    const host = $("#plot-host");
    if (!host._clWired) {
      host.on("plotly_click", (ev) => {
        const pt = ev.points && ev.points[0];
        const g = pt && pt.customdata ? pt.customdata[0] : (pt && pt.text);
        if (g && g in state.byGene) setGene(g);
      });
      host._clWired = true;
    }
  }

  function renderReadout() {
    const g = state.gene, r = state.byGene[g], c = clusterOf(g), info = state.clusters[c];
    const NR = state.data.meta.n_rad;
    const rad = r.profile.slice(0, NR), ax = r.profile.slice(NR);
    const arg = (a) => a.indexOf(Math.max(...a));
    const radWord = ["innermost core", "inner", "mid-inner", "mid-outer", "outer", "cortex"][arg(rad)] || "—";
    const axWord = ["far vegetal", "vegetal", "mid-vegetal", "mid-animal", "animal", "far animal"][arg(ax)] || "—";
    $("#cl-readout").innerHTML =
      `<div class="cl-r-gene">${g}</div>` +
      `<div class="cl-r-line"><span class="cl-sw" style="background:${DARK}"></span>` +
        `<span class="k">cluster</span><span class="v">${c + 1} of ${state.clusters.length}</span></div>` +
      `<div class="cl-r-line"><span class="k">cluster size</span><span class="v">${info.n} genes</span></div>` +
      `<div class="cl-r-line"><span class="k">measured in</span><span class="v">${r.n_emb} zygotes</span></div>` +
      `<div class="cl-r-line"><span class="k">most enriched</span><span class="v">${radWord}</span></div>` +
      `<div class="cl-r-line"><span class="k">along the axis</span><span class="v">${axWord}</span></div>`;
  }

  // ───────── right drawer: who else is in this cluster ─────────
  function renderMembers() {
    const g = state.gene, c = clusterOf(g), info = state.clusters[c];
    $("#cl-mem-title").textContent = `Cluster ${c + 1} · ${info.n} genes`;
    $("#cl-mem-desc").textContent =
      "Ordered by how central each gene is to the cluster — the first is the cluster's most typical " +
      "member. Click any gene to select it.";
    const cen = info.profile;
    const rows = info.members.map((m) => ({ g: m, d: dist(state.byGene[m].profile, cen) }));
    $("#cl-members").innerHTML = rows.map((r, i) =>
      `<div class="cl-mem${r.g === g ? " current" : ""}" data-gene="${r.g}" title="${r.g} · distance to cluster centre ${r.d.toFixed(3)}">` +
      `<span class="n">${i + 1}</span><span class="g" data-gene="${r.g}">${r.g}</span>` +
      `<span class="d">${r.d.toFixed(2)}</span></div>`).join("");
    $("#cl-members").querySelectorAll(".cl-mem").forEach((el) =>
      el.addEventListener("click", () => setGene(el.dataset.gene)));
  }

  /** Correlation distance — the measure the clustering itself used. */
  function dist(a, b) {
    const n = a.length;
    let ma = 0, mb = 0;
    for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
    ma /= n; mb /= n;
    let num = 0, da = 0, db = 0;
    for (let i = 0; i < n; i++) {
      const x = a[i] - ma, y = b[i] - mb;
      num += x * y; da += x * x; db += y * y;
    }
    const den = Math.sqrt(da * db);
    return den > 0 ? 1 - num / den : 1;
  }

  // ───────── drawer panels ─────────
  const baseLayout = (yt) => ({
    margin: { l: 58, r: 12, t: 10, b: 42 }, showlegend: false,
    paper_bgcolor: "transparent", plot_bgcolor: "#fcfdfe",
    xaxis: { gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
    yaxis: { title: { text: yt, font: { size: 11 } }, gridcolor: "#eef1f5",
             zerolinecolor: "#cbd5e1", tickfont: { size: 9 } },
    font: { size: 11, color: "#334155" },
  });
  const binLabels = () => {
    const m = state.data.meta;
    const rad = ["core", "", "", "", "", "cortex"].map((s, i) => s || `r${i + 1}`);
    const ax = ["vegetal", "", "", "", "", "animal"].map((s, i) => s || `a${i + 1}`);
    return { rad, ax, NR: m.n_rad };
  };

  function renderSig() {
    const div = $("#cl-sig"); if (!shown(div)) return;
    const g = state.gene, r = state.byGene[g], c = clusterOf(g), info = state.clusters[c];
    const { rad, ax, NR } = binLabels();
    const labels = rad.concat(ax);
    $("#cl-sig-sub").textContent = `· ${g} · cluster ${c + 1} · ${r.n_emb} zygotes`;
    resetPlot(div);
    const traces = [
      { type: "bar", name: `cluster ${c + 1} mean`, x: labels, y: info.profile,
        marker: { color: DARK, opacity: .35 },
        hovertemplate: "cluster mean · %{x}: %{y:.3f}<extra></extra>" },
      { type: "scatter", mode: "lines+markers", name: g, x: labels, y: r.profile,
        line: { color: SEL, width: 2.5 }, marker: { size: 7, color: SEL },
        hovertemplate: `${g} · %{x}: %{y:.3f}<extra></extra>` },
    ];
    const lay = baseLayout("log2 enrichment vs the embryo's own transcripts");
    lay.showlegend = true; lay.margin.t = 34;
    lay.legend = { orientation: "h", y: 1.03, x: 0, yanchor: "bottom", font: { size: 9 } };
    lay.shapes = [{ type: "line", x0: -0.5, x1: labels.length - 0.5, y0: 0, y1: 0,
                    line: { color: "#94a3b8", width: 1, dash: "dot" } },
                  { type: "line", x0: NR - 0.5, x1: NR - 0.5, yref: "paper", y0: 0, y1: 1,
                    line: { color: "#cbd5e1", width: 1 } }];
    lay.annotations = [
      { x: (NR - 1) / 2, y: 1, yref: "paper", text: "RADIAL", showarrow: false,
        font: { size: 9, color: "#94a3b8" }, yanchor: "bottom" },
      { x: NR + (NR - 1) / 2, y: 1, yref: "paper", text: "AXIAL", showarrow: false,
        font: { size: 9, color: "#94a3b8" }, yanchor: "bottom" },
    ];
    Plotly.react(div, traces, lay, { responsive: true, displaylogo: false, displayModeBar: false });
  }

  function renderProfiles() {
    const div = $("#cl-prof"); if (!shown(div)) return;
    const { rad, ax, NR } = binLabels();
    const labels = rad.concat(ax);
    const selC = clusterOf(state.gene);
    $("#cl-prof-sub").textContent = `· k = ${state.k}`;
    resetPlot(div);
    const z = state.clusters.map((c) => c.profile);
    const yl = state.clusters.map((c, i) => `C${i + 1} (${c.n})${i === selC ? "  ◀" : ""}`);
    const traces = [{
      type: "heatmap", z, x: labels, y: yl,
      colorscale: [[0, "#3b4a6b"], [0.5, "#f7f7f5"], [1, "#8c3b2f"]],
      zmid: 0, colorbar: { title: { text: "log2", font: { size: 9 } }, thickness: 10,
                           tickfont: { size: 9 }, len: .9 },
      hovertemplate: "%{y} · %{x}: %{z:.3f}<extra></extra>",
    }];
    const lay = baseLayout("");
    lay.margin.l = 76;
    lay.shapes = [{ type: "line", x0: NR - 0.5, x1: NR - 0.5, yref: "paper", y0: 0, y1: 1,
                    line: { color: "#94a3b8", width: 1 } }];
    Plotly.react(div, traces, lay, { responsive: true, displaylogo: false, displayModeBar: false });
  }

  function renderNear() {
    const div = $("#cl-near"); if (!shown(div)) return;
    const g = state.gene, p = state.byGene[g].profile, selC = clusterOf(g);
    const rows = state.genes.filter((x) => x !== g)
      .map((x) => ({ g: x, d: dist(state.byGene[x].profile, p), c: clusterOf(x) }))
      .sort((a, b) => a.d - b.d).slice(0, 24).reverse();
    $("#cl-near-sub").textContent = `· nearest 24 to ${g}`;
    resetPlot(div);
    const traces = [{
      type: "bar", orientation: "h",
      x: rows.map((r) => r.d), y: rows.map((r) => r.g),
      marker: {
        color: rows.map((r) => (r.c === selC ? DARK : GREY(r.c, state.clusters.length))),
        line: { width: 0.5, color: "rgba(255,255,255,.7)" },
      },
      customdata: rows.map((r) => r.c + 1),
      hovertemplate: "<b>%{y}</b> · cluster %{customdata}<br>distance %{x:.3f}<extra></extra>",
    }];
    const lay = baseLayout("");
    lay.margin.l = 92; lay.xaxis.title = { text: "correlation distance (0 = identical)", font: { size: 11 } };
    lay.yaxis.tickfont = { size: 9 };
    Plotly.react(div, traces, lay, { responsive: true, displaylogo: false, displayModeBar: false });
    div.on("plotly_click", (ev) => {
      const y = ev.points && ev.points[0] && ev.points[0].y;
      if (y && y in state.byGene) setGene(y);
    });
  }

  function renderMethod() {
    const m = state.data.meta;
    $("#cl-method").innerHTML = `
      <h4>The question</h4>
      <p>Not which genes are expressed together — <b>which genes sit in the same place</b>. Two genes
      land in the same cluster when, embryo after embryo, they are enriched in the same shell of the
      cell and the same end of the polar-body axis.</p>

      <h4>Each gene's signature</h4>
      <p>For every zygote and every gene, over <b>cell-body transcripts only</b> (pronuclei and polar
      body are separate compartments and would dominate):
      <b>${m.n_rad} radial shells</b> of equal <i>volume</i> — so a uniformly scattered gene gives a flat
      line — and <b>${m.n_ax} bins</b> along the polar-body axis. Each gene's profile is divided by
      <b>that embryo's own all-transcript profile</b> and logged, which is what removes cell shape,
      size, orientation and detection efficiency. Signatures are then averaged over the embryos
      carrying the gene.</p>

      <h4>Why not "genes that co-vary across embryos"</h4>
      <div class="warn">Each embryo is imaged with one probeset, and the probesets are <b>disjoint</b> —
      one panel shares zero genes with the others. Genes from different panels are therefore never
      measured in the same embryo, so a co-variation method could only ever relate genes inside a
      single panel. Giving every gene its own signature in a shared frame is what lets all
      ${m.n_genes} be compared at once.</p></div>

      <h4>Sparse bins are shrunk, not trusted</h4>
      <p>The axial bins nearest the poles are geometrically thin and legitimately hold few transcripts.
      A plain pseudocount would turn an empty one into a huge negative spike and let sampling noise
      drive the clustering, so bins are shrunk toward "no enrichment" in proportion to how little
      evidence they carry (<code>${"α"} = 20</code> pseudo-transcripts distributed like the background).</p>

      <h4>Clustering</h4>
      <p>Correlation distance between mean signatures, then Ward hierarchical clustering. A gene needs
      <b>≥${m.min_tx} cell-body transcripts in ≥${m.min_emb} zygotes</b> to be included, which is what
      reduces the ${420} probed genes to <b>${m.n_genes}</b>. Labels ship for every k in
      ${m.k_range[0]}–${m.k_range[m.k_range.length - 1]}; k=${m.default_k} scores best by silhouette.</p>

      <h4>Reading it honestly</h4>
      <p>The best silhouette here is <b>${m.silhouette[String(m.default_k)].toFixed(3)}</b>. That is real
      structure, not a razor-sharp partition — spatial preference in the zygote is a gradient, and the
      clusters are cuts through it. Treat cluster membership as "these genes lean the same way",
      and use the <b>Nearest genes</b> panel rather than the cut when a specific pair matters.</p>`;
  }

  // ───────── plumbing ─────────
  const shown = (el) => !!(el && el.offsetParent);
  function resetPlot(div) {
    if (div._fullLayout && div.classList.contains("js-plotly-plot")) return;
    div.innerHTML = "";
  }
  const RENDER = { sig: renderSig, profiles: renderProfiles, near: renderNear, method: () => {} };
  function renderActive() { (RENDER[state.tab] || renderSig)(); }

  function openDrawer(open) {
    state.drawerOpen = open;
    $("#drawer").dataset.open = open ? "true" : "false";
    $("#drawer-handle").setAttribute("aria-expanded", String(open));
    if (open) renderActive();
  }

  function wire() {
    $("#gene-select").addEventListener("change", (e) => setGene(e.target.value));
    $("#k-range").addEventListener("input", (e) => {
      applyK(parseInt(e.target.value, 10));
      renderMap(); renderReadout(); renderMembers(); renderActive();
    });
    $("#layout-sel").addEventListener("change", (e) => { state.layout = e.target.value; renderMap(); });
    $("#show-labels").addEventListener("change", (e) => { state.showLabels = e.target.checked; renderMap(); });

    // bottom drawer
    $("#drawer-handle").addEventListener("click", () => openDrawer($("#drawer").dataset.open !== "true"));
    $("#cl-tabs").addEventListener("click", (e) => {
      const b = e.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#cl-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#cl-panels").querySelectorAll(".xs-panel").forEach((p) =>
        (p.hidden = p.dataset.tab !== state.tab));
      renderActive();
    });
    const rz = $("#drawer-resize"); let sh = 0;
    rz.addEventListener("pointerdown", (e) => {
      sh = $("#drawer-body").getBoundingClientRect().height; rz._d = { y: e.clientY };
      rz.setPointerCapture(e.pointerId); e.preventDefault();
    });
    rz.addEventListener("pointermove", (e) => {
      if (!rz._d) return;
      $("#drawer").style.setProperty("--drawer-h",
        Math.max(200, Math.min(window.innerHeight - 120, sh + (rz._d.y - e.clientY))) + "px");
    });
    const end = (e) => { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);

    // right drawer
    const rd = $("#rdrawer"), rh = $("#rdrawer-handle");
    rh.addEventListener("click", () => {
      const o = rd.dataset.open !== "true";
      rd.dataset.open = o ? "true" : "false";
      rh.setAttribute("aria-expanded", String(o));
      if (o) renderMembers();
    });

    // header search
    const box = $("#gene-search"), sug = $("#gene-suggest");
    const close = () => { sug.hidden = true; sug.innerHTML = ""; };
    box.addEventListener("input", () => {
      const q = box.value.trim().toLowerCase();
      if (!q) return close();
      const hits = state.genes.filter((g) => g.toLowerCase().includes(q)).slice(0, 12);
      if (!hits.length) return close();
      const n = state.clusters.length;
      sug.innerHTML = hits.map((g) => {
        const c = clusterOf(g);
        return `<div class="cl-sg" data-gene="${g}">` +
          `<span class="cdot" style="background:${GREY(c, n)}"></span>${g}` +
          `<span class="cn">C${c + 1}</span></div>`;
      }).join("");
      sug.hidden = false;
      sug.querySelectorAll(".cl-sg").forEach((el) => el.addEventListener("mousedown", (ev) => {
        ev.preventDefault(); setGene(el.dataset.gene); box.value = ""; close();
      }));
    });
    box.addEventListener("blur", () => setTimeout(close, 120));
    box.addEventListener("keydown", (e) => {
      if (e.key === "Escape") { box.value = ""; close(); box.blur(); }
      if (e.key === "Enter") {
        const first = sug.querySelector(".cl-sg");
        if (first) { setGene(first.dataset.gene); box.value = ""; close(); box.blur(); }
      }
    });

    window.addEventListener("resize", () => { try { Plotly.Plots.resize($("#plot-host")); } catch (_) {} });
  }
})();
