/* Fertilization Geometry & the First Embryonic Axis.
 * Does the sperm entry position (GFP midpiece) and the polar body (animal pole)
 * predict the first cleavage plane? Geometry + a molecular-asymmetry layer, tested
 * against an exact random-orientation null, with an embryo-shape control. All
 * numbers are precomputed (build_axes.py); the front-end reads them and pools the
 * per-embryo angles into the aggregate test.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const C = { av: "#7c3aed", sperm: "#059669", plane: "#f97316", shape: "#64748b",
              com: "#111827", nuc: "#2563eb", body: "#f97316", side0: "#e11d48", side1: "#2563eb" };

  function plotInto(div, traces, layout, cfg) {
    if (!div.classList.contains("js-plotly-plot")) div.innerHTML = "";
    Plotly.react(div, traces, layout, cfg || { responsive: true, displayModeBar: false });
  }
  // Φ(z): standard-normal CDF via an Abramowitz–Stegun erf approximation.
  function normCdf(z) {
    const t = 1 / (1 + 0.2316419 * Math.abs(z));
    const d = 0.3989422804 * Math.exp(-z * z / 2);
    let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
    return z > 0 ? 1 - p : p;
  }
  // two-sided binomial sign test: k of n successes vs p=0.5
  function binomTwoSided(k, n) {
    if (n === 0) return 1;
    const logC = (n, k) => { let s = 0; for (let i = 0; i < k; i++) s += Math.log(n - i) - Math.log(i + 1); return s; };
    let cum = 0; for (let i = 0; i <= n; i++) cum += Math.exp(logC(n, i) - n * Math.log(2));  // ~1
    const pmf = (i) => Math.exp(logC(n, i) - n * Math.log(2));
    let lo = 0; for (let i = 0; i <= Math.min(k, n - k); i++) lo += pmf(i);
    return Math.min(1, 2 * lo);
  }

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const togglesEl = $("#ax-toggles"), geneSelect = $("#gene-select"), readoutEl = $("#ax-readout");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle"), rdrawerDesc = $("#rdrawer-desc"), geneListEl = $("#gene-list");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), drawerBody = $("#drawer-body");
  const verdictEl = $("#ax-verdict"), histsEl = $("#ax-hists");

  const state = { manifest: null, mol: null, currentId: null, scene: null, userGene: null, drawerOpen: false,
                  show: { av: true, sperm: true, plane: true, shape: false, dots: false } };

  const TOGGLES = [
    { key: "av", label: "Animal–Vegetal axis", desc: "COM → polar body — the polar body is the accepted animal-pole landmark." },
    { key: "sperm", label: "Sperm axis", desc: "COM → GFP midpiece — a proxy for the sperm entry position (stays near the fusion site)." },
    { key: "plane", label: "Cleavage plane / pronuclear axis", desc: "two-cell: the interface between the two blastomeres (the realized first cleavage plane); zygote: the pronuclear axis (predicted normal)." },
    { key: "shape", label: "Embryo long axis (shape control)", desc: "the embryo's own elongation axis — if the plane just follows this, the correlation is a shape artefact, not pre-patterning." },
    { key: "dots", label: "Gene transcript dots", desc: "the selected gene's transcripts, coloured by which side of the plane / blastomere they fall in." },
  ];

  (async function init() {
    try {
      const [m, mol] = await Promise.all([
        (await fetch("data/axes_manifest.json")).json(),
        V.loadGz("data/axes_molecular.json.gz"),
      ]);
      state.manifest = m; state.mol = mol;
      countEl.textContent = `${m.n_zygote} zygotes + ${m.n_twocell} two-cell · sperm-positive`;
      V.buildTabs(tabsEl, m.embryos, selectEmbryo, (e) => ({
        label: e.label, sub: e.stage === "zygote" ? "zygote" : "2-cell",
        title: `${e.stage} · ${e.id}`,
      }));
      buildToggles();
      V.wireWindow(controlsEl, $("#controls-header"), [...controlsEl.querySelectorAll(".rz")], "axes_controls_box");
      wireRdrawer(); wireDrawer();
      geneSelect.addEventListener("change", () => { state.userGene = geneSelect.value; render(); renderReadout(); });
    } catch (err) { showError("Failed to load: " + (err.message || err)); }
  })();

  function buildToggles() {
    togglesEl.innerHTML = TOGGLES.map((t) =>
      `<label class="ax-toggle"><input type="checkbox" data-key="${t.key}" ${state.show[t.key] ? "checked" : ""}/>` +
      `<span class="ax-t-main"><span class="ax-t-dot" style="background:${C[t.key] || "#999"}"></span>${t.label}</span>` +
      `<span class="ax-t-desc">${t.desc}</span></label>`).join("");
    togglesEl.querySelectorAll("input").forEach((c) => c.addEventListener("change", () => {
      state.show[c.dataset.key] = c.checked; render();
    }));
  }

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    V.markActiveTab(tabsEl, id);
    showLoading(`Loading ${id}…`);
    try {
      const scene = await V.loadGz(`data/axes/${id}.json.gz`);
      if (state.currentId !== id) return;
      state.scene = scene;
      populateGenes(scene);
      controlsEl.hidden = false; placeholder.hidden = true; rdrawer.hidden = false; drawer.hidden = false;
      render(); renderReadout(); renderGeneList(); renderAggregate();
    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }

  function populateGenes(scene) {
    const genes = Object.keys(scene.transcripts).sort((a, b) => a.localeCompare(b, undefined, { sensitivity: "base", numeric: true }));
    geneSelect.innerHTML = genes.map((g) => `<option value="${g}">${g}</option>`).join("");
    geneSelect.value = (state.userGene && scene.transcripts[state.userGene]) ? state.userGene : genes[0];
    state.userGene = geneSelect.value;
  }
  const gene = () => geneSelect.value;

  // ---------- 3-D ----------
  function segMesh(scene, lbl, color, opacity) {
    const mesh = scene.region_meshes[String(lbl)]; if (!mesh) return null;
    const v = mesh.verts, f = mesh.faces, nV = v.length / 3, nF = f.length / 3;
    const x = new Array(nV), y = new Array(nV), z = new Array(nV);
    for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
    const ii = new Array(nF), jj = new Array(nF), kk = new Array(nF);
    for (let i = 0; i < nF; i++) { ii[i] = f[i * 3]; jj[i] = f[i * 3 + 1]; kk[i] = f[i * 3 + 2]; }
    return { type: "mesh3d", x, y, z, i: ii, j: jj, k: kk, color, opacity, name: `seg ${lbl}`,
      showlegend: false, flatshading: false, hoverinfo: "skip",
      lighting: { ambient: 0.75, diffuse: 0.5, specular: 0.1, roughness: 0.9 } };
  }
  const sub3 = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const add3 = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const mul3 = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const norm3 = (a) => Math.hypot(a[0], a[1], a[2]);
  const unit3 = (a) => { const n = norm3(a) || 1; return [a[0] / n, a[1] / n, a[2] / n]; };
  const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  function line(p0, p1, color, width, name, dash) {
    return { type: "scatter3d", mode: "lines", x: [p0[0], p1[0]], y: [p0[1], p1[1]], z: [p0[2], p1[2]],
      line: { color, width, dash: dash || "solid" }, name, hovertemplate: `${name}<extra></extra>`, showlegend: true };
  }
  function marker(p, color, symbol, size, name) {
    return { type: "scatter3d", mode: "markers", x: [p[0]], y: [p[1]], z: [p[2]],
      marker: { color, symbol: symbol || "circle", size: size || 6, line: { width: 1, color: "#fff" } },
      name, hovertemplate: `${name}<extra></extra>`, showlegend: true };
  }
  function planeDisc(center, normal, R, color) {
    const n = unit3(normal);
    let ref = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const u = unit3(cross3(n, ref)), v = unit3(cross3(n, u));
    const N = 28, cx = [], cy = [], cz = [], I = [], J = [], K = [];
    cx.push(center[0]); cy.push(center[1]); cz.push(center[2]);              // 0 = center
    for (let a = 0; a < N; a++) {
      const t = a / N * 2 * Math.PI;
      const p = add3(center, add3(mul3(u, R * Math.cos(t)), mul3(v, R * Math.sin(t))));
      cx.push(p[0]); cy.push(p[1]); cz.push(p[2]);
    }
    for (let a = 0; a < N; a++) { I.push(0); J.push(1 + a); K.push(1 + (a + 1) % N); }
    return { type: "mesh3d", x: cx, y: cy, z: cz, i: I, j: J, k: K, color, opacity: 0.16,
      name: "cleavage plane", showlegend: true, hoverinfo: "name", flatshading: true };
  }
  function render() {
    const s = state.scene; if (!s) return;
    const ex = s.extents, R = 0.55 * Math.max(ex.x[1] - ex.x[0], ex.y[1] - ex.y[0], ex.z[1] - ex.z[0]);
    const com = s.com_plot, L = s.landmarks, isZyg = s.stage === "zygote";
    const traces = [];
    for (const lbl of s.mask_labels) { const t = segMesh(s, lbl, "#9aa3b2", 0.07); if (t) traces.push(t); }
    traces.push(marker(com, C.com, "circle", 5, "Embryo COM"));

    if (state.show.av && L.polar_plot) {
      const dir = unit3(sub3(L.polar_plot, com));
      traces.push(line(add3(com, mul3(dir, -R)), add3(com, mul3(dir, R)), C.av, 6, "Animal–Vegetal axis"));
      traces.push(marker(L.polar_plot, C.av, "diamond", 7, "Polar body (animal pole)"));
    }
    if (state.show.sperm) {
      const dir = unit3(sub3(L.sperm_plot, com));
      traces.push(line(com, add3(com, mul3(dir, R * 1.15)), C.sperm, 6, "Sperm axis"));
      traces.push(marker(L.sperm_plot, C.sperm, "diamond", 8, "Sperm midpiece"));
    }
    if (state.show.plane) {
      if (isZyg && L.nuclei_plots.length >= 2) {
        traces.push(line(L.nuclei_plots[0], L.nuclei_plots[1], C.plane, 6, "Pronuclear axis"));
        L.nuclei_plots.forEach((p, i) => traces.push(marker(p, C.nuc, "circle", 6, `Pronucleus ${i + 1}`)));
      } else if (!isZyg && L.body_plots.length === 2) {
        const mid = mul3(add3(L.body_plots[0], L.body_plots[1]), 0.5);
        const nrm = sub3(L.body_plots[1], L.body_plots[0]);
        traces.push(planeDisc(mid, nrm, R * 0.9, C.plane));
        traces.push(line(L.body_plots[0], L.body_plots[1], C.plane, 4, "Cleavage normal", "dot"));
        L.body_plots.forEach((p, i) => traces.push(marker(p, i === L.sperm_body ? C.sperm : C.body, "circle", 7,
          i === L.sperm_body ? "Blastomere (sperm side)" : "Blastomere")));
      }
    }
    if (state.show.shape) {
      const dir = unit3(sub3(s.long_axis_plot, com));
      traces.push(line(add3(com, mul3(dir, -R)), add3(com, mul3(dir, R)), C.shape, 4, "Embryo long axis", "dash"));
    }
    if (state.show.dots) {
      const tx = s.transcripts[gene()];
      if (tx && tx.x.length) {
        const zs = s.z_scale;
        traces.push({ type: "scatter3d", mode: "markers", name: `${gene()} · side 0`,
          x: tx.x.filter((_, i) => tx.s[i] === 0), y: tx.y.filter((_, i) => tx.s[i] === 0),
          z: tx.gz.filter((_, i) => tx.s[i] === 0).map((z) => z * zs),
          marker: { size: 0.5, color: C.side0, opacity: 0.8 }, hovertemplate: `${gene()}<extra></extra>` });
        traces.push({ type: "scatter3d", mode: "markers", name: `${gene()} · side 1`,
          x: tx.x.filter((_, i) => tx.s[i] === 1), y: tx.y.filter((_, i) => tx.s[i] === 1),
          z: tx.gz.filter((_, i) => tx.s[i] === 1).map((z) => z * zs),
          marker: { size: 0.5, color: C.side1, opacity: 0.8 }, hovertemplate: `${gene()}<extra></extra>` });
      }
    }
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  // ---------- per-embryo readout ----------
  const interp = (a) => a == null ? "—" : (a.cos < 0.5
    ? `${a.deg.toFixed(0)}° — near-perpendicular → landmark lies IN the plane`
    : `${a.deg.toFixed(0)}° — landmark toward a pole (plane ⟂ to it)`);
  function renderReadout() {
    const s = state.scene; if (!s) return;
    const a = s.angles, isZyg = s.stage === "zygote";
    let html = `<div class="ax-r-title">${s.stage_label} · ${s.id}</div>`;
    html += `<div class="ax-r-row"><span>${isZyg ? "Pronuclear axis" : "Cleavage plane"} vs polar body</span><b>${interp(a.plane_vs_pb)}</b></div>`;
    html += `<div class="ax-r-row"><span>${isZyg ? "Pronuclear axis" : "Cleavage plane"} vs sperm</span><b>${interp(a.plane_vs_sperm)}</b></div>`;
    html += `<div class="ax-r-row ax-r-shape"><span>vs embryo long axis (shape)</span><b>${a.plane_vs_shape ? a.plane_vs_shape.deg.toFixed(0) + "°" : "—"}${a.plane_vs_shape && a.plane_vs_shape.cos > 0.5 ? " — plane follows shape" : ""}</b></div>`;
    html += `<div class="ax-r-row"><span>embryo elongation (aspect)</span><b>${a.aspect_ratio}×</b></div>`;
    // gene asymmetry in this embryo
    const g = gene(), row = state.mol.embryos.find((e) => e.id === s.id);
    const sd = row && row.sides[g];
    if (sd) {
      const tot = sd[0] + sd[1], f0 = tot ? (100 * sd[0] / tot).toFixed(0) : 0;
      const s0lbl = isZyg ? "animal half" : "sperm-side blastomere";
      html += `<div class="ax-r-gene"><b>${g}</b>: ${f0}% on the ${s0lbl} (${sd[0]} vs ${sd[1]})</div>`;
    }
    readoutEl.innerHTML = html;
  }

  // ---------- bottom drawer: aggregate test ----------
  function statFor(cohort, metric) {
    const cs = state.manifest.embryos.filter((e) => e.stage === cohort && e.angles[metric])
      .map((e) => e.angles[metric].cos);
    const n = cs.length; if (!n) return null;
    const mean = cs.reduce((a, b) => a + b, 0) / n;
    const z = (mean - 0.5) / Math.sqrt((1 / 12) / n);            // null: |cos| ~ U(0,1), var = 1/12
    const p = 2 * (1 - normCdf(Math.abs(z)));
    const meanDeg = cs.map((c) => Math.acos(Math.min(1, c)) * 180 / Math.PI).reduce((a, b) => a + b, 0) / n;
    return { cs, n, mean, p, meanDeg };
  }
  function renderAggregate() {
    const s = state.scene; if (!s) return;
    const cohort = s.stage;
    const label = cohort === "zygote" ? "zygotes (pronuclear axis)" : "two-cell embryos (realized cleavage plane)";
    const metrics = [
      { key: "plane_vs_pb", title: "vs polar body (animal–vegetal axis)" },
      { key: "plane_vs_sperm", title: "vs sperm entry (midpiece)" },
      { key: "plane_vs_shape", title: "vs embryo long axis — SHAPE CONTROL" },
    ];
    const stats = metrics.map((m) => ({ ...m, s: statFor(cohort, m.key) }));
    const pb = stats[0].s, sp = stats[1].s, sh = stats[2].s;
    const planeName = cohort === "zygote" ? "The pronuclear axis" : "The first cleavage plane";
    const stat = (st) => `mean angle ${st.meanDeg.toFixed(0)}°, mean|cos|=${st.mean.toFixed(2)}, p=${fmtP(st.p)}, n=${st.n}`;
    const pbLine = !pb ? "" : (pb.mean < 0.5
      ? `<b style="color:${C.av}">passes near the polar body</b> (the animal–vegetal axis lies in the plane) — ${stat(pb)}.`
      : `<b style="color:${C.av}">does not pass near the polar body</b> (it sits toward a pole) — ${stat(pb)}.`);
    const spLine = !sp ? "" : (sp.mean < 0.5
      ? `<b style="color:${C.sperm}">passes near the sperm entry point</b> — ${stat(sp)}.`
      : `<b style="color:${C.sperm}">does NOT pass through the sperm</b> (the midpiece sits toward one blastomere) — ${stat(sp)}.`);
    const shLine = !sh ? "" : (sh.mean > 0.5
      ? `lies almost <b>along the embryo's own long axis</b> (${stat(sh)}) — so it strongly follows embryo elongation, and the animal–vegetal alignment above is largely a consequence of shape, not independent pre-patterning.`
      : `is <b>not simply following embryo shape</b> (${stat(sh)}) — so the alignments above are not merely a shape artefact.`);
    verdictEl.innerHTML =
      `<div class="ax-v-title">Across ${pb ? pb.n : "?"} ${label}:</div>` +
      `<div class="ax-v-line">${planeName} ${pbLine}</div>` +
      `<div class="ax-v-line">It ${spLine}</div>` +
      `<div class="ax-v-line ax-v-shape"><b style="color:${C.shape}">Shape control:</b> the plane ${shLine}</div>` +
      `<div class="ax-v-null">Null: for a randomly oriented plane, |cos θ| between the plane normal and a landmark direction is Uniform(0,1) (mean 0.5); mean &lt; 0.5 ⇒ the landmark lies in the plane, &gt; 0.5 ⇒ toward a pole. Tested by z-score against the exact null variance 1/12.</div>`;
    // histograms of |cos| with the flat null
    plotInto(histsEl, stats.filter((m) => m.s).map((m, idx) => ({
      type: "histogram", x: m.s.cs, xbins: { start: 0, end: 1, size: 0.1 },
      marker: { color: [C.av, C.sperm, C.shape][idx] }, opacity: 0.8, name: m.title, xaxis: `x${idx + 1}`, yaxis: `y${idx + 1}`,
    })), histLayout(stats));
  }
  function histLayout(stats) {
    const lay = { showlegend: false, height: histsEl.clientHeight || 220, margin: { l: 30, r: 8, t: 24, b: 30 },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      grid: { rows: 1, columns: 3, pattern: "independent" }, annotations: [] };
    stats.forEach((m, i) => {
      const ax = i === 0 ? "" : i + 1;
      lay[`xaxis${ax}`] = { title: { text: "|cos θ|", font: { size: 9 } }, range: [0, 1], tickfont: { size: 8 }, gridcolor: "#eef1f5" };
      lay[`yaxis${ax}`] = { tickfont: { size: 8 }, gridcolor: "#eef1f5" };
      lay.annotations.push({ text: m.title + (m.s ? `<br>mean|cos|=${m.s.mean.toFixed(2)} · p=${fmtP(m.s.p)}` : ""),
        showarrow: false, x: 0.5, xref: `x${ax} domain`, y: 1.18, yref: `y${ax} domain`, font: { size: 9, color: "#475569" } });
    });
    return lay;
  }
  const fmtP = (p) => p < 0.001 ? p.toExponential(1) : p.toFixed(3);

  // ---------- right drawer: molecular asymmetry ranking ----------
  function renderGeneList() {
    const s = state.scene; if (!s) return;
    const cohort = s.stage, isZyg = cohort === "zygote";
    const s0 = isZyg ? "animal half" : "sperm-side blastomere", s1 = isZyg ? "vegetal half" : "other blastomere";
    rdrawerDesc.innerHTML = `Per gene, the fraction of transcripts on the <b>${s0}</b> vs <b>${s1}</b>, pooled across ` +
      `${cohort === "zygote" ? "zygotes" : "two-cell embryos"}. Ranked by how consistently a gene is asymmetric ` +
      `(sign test across embryos). Tests whether the transcriptome is pre-patterned along this axis.`;
    // aggregate per gene
    const embs = state.mol.embryos.filter((e) => e.stage === cohort);
    const agg = {};
    for (const e of embs) for (const [g, sd] of Object.entries(e.sides)) {
      const tot = sd[0] + sd[1]; if (tot < 10) continue;                 // ignore too-sparse
      (agg[g] = agg[g] || []).push(sd[0] / tot);
    }
    const rows = Object.entries(agg).filter(([, fr]) => fr.length >= 3).map(([g, fr]) => {
      const mean = fr.reduce((a, b) => a + b, 0) / fr.length;
      const k = fr.filter((x) => x > 0.5).length;
      return { gene: g, mean, n: fr.length, p: binomTwoSided(k, fr.length), bias: Math.abs(mean - 0.5) };
    }).sort((a, b) => a.p - b.p || b.bias - a.bias);
    const cur = gene();
    let html = `<div class="best-head"><span></span><span>gene</span><span>% s0</span><span>n</span><span>p</span></div>`;
    html += rows.slice(0, 60).map((r, i) =>
      `<div class="best-row seg-row${r.gene === cur ? " current" : ""}" data-gene="${r.gene}" ` +
      `title="${r.gene}: ${(100 * r.mean).toFixed(0)}% on ${s0} across ${r.n} embryos">` +
      `<span class="best-num">${i + 1}</span><span class="best-gene">${r.gene}</span>` +
      `<span class="best-real">${(100 * r.mean).toFixed(0)}%</span>` +
      `<span class="best-p">${r.n}</span>` +
      `<span class="best-p${r.p <= 0.05 ? " sig" : ""}">${fmtP(r.p)}</span></div>`).join("");
    geneListEl.innerHTML = html || `<div class="best-plane-note">No genes with enough coverage.</div>`;
  }
  function highlightGene() {
    const cur = gene();
    geneListEl.querySelectorAll(".seg-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === cur));
  }

  // ---------- drawers ----------
  function wireRdrawer() {
    wireHandleDrag(rdrawer, rdrawerHandle, {
      computeSize: (e) => window.innerWidth - e.clientX,
      clampSize: (px) => Math.max(260, Math.min(window.innerWidth - 80, px)),
      applySize: (px) => rdrawer.style.setProperty("--rdrawer-w", px + "px"),
      setOpen: (o) => { rdrawer.dataset.open = o ? "true" : "false"; rdrawerHandle.setAttribute("aria-expanded", String(o)); },
    });
    let sw = 0; const rrz = $("#rdrawer-resize");
    rrz.addEventListener("pointerdown", (e) => { sw = rdrawer.getBoundingClientRect().width; rrz._d = { x: e.clientX }; rrz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rrz.addEventListener("pointermove", (e) => { if (!rrz._d) return;
      rdrawer.style.setProperty("--rdrawer-w", Math.max(260, Math.min(window.innerWidth - 80, sw - (e.clientX - rrz._d.x))) + "px"); });
    const end = (e) => { if (rrz._d) { rrz._d = null; try { rrz.releasePointerCapture(e.pointerId); } catch (_) {} } };
    rrz.addEventListener("pointerup", end); rrz.addEventListener("pointercancel", end);
    geneListEl.addEventListener("click", (e) => {
      const row = e.target.closest(".seg-row"); if (!row) return;
      const g = row.dataset.gene;
      if (state.scene.transcripts[g]) { state.userGene = g; geneSelect.value = g; render(); renderReadout(); highlightGene(); }
    });
  }
  function wireDrawer() {
    const setOpen = (o) => {
      state.drawerOpen = o; drawer.dataset.open = o ? "true" : "false"; drawerHandle.setAttribute("aria-expanded", String(o));
      if (o) { renderAggregate(); requestAnimationFrame(() => { try { Plotly.Plots.resize(histsEl); } catch (_) {} }); }
    };
    wireHandleDrag(drawer, drawerHandle, {
      computeSize: (e) => window.innerHeight - e.clientY - 40,
      clampSize: (px) => Math.max(200, Math.min(window.innerHeight - 100, px)),
      applySize: (px) => drawer.style.setProperty("--drawer-h", px + "px"),
      setOpen, afterDrag: () => { try { Plotly.Plots.resize(histsEl); } catch (_) {} },
    });
    let sh = 0; const rz = $("#drawer-resize");
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return;
      drawer.style.setProperty("--drawer-h", Math.max(200, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y))) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} try { Plotly.Plots.resize(histsEl); } catch (_) {} } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }
  function wireHandleDrag(el, handle, cfg) {
    let start = null, moved = false;
    handle.addEventListener("pointerdown", (e) => { if (e.button && e.button !== 0) return; start = { x: e.clientX, y: e.clientY }; moved = false; try { handle.setPointerCapture(e.pointerId); } catch (_) {} });
    handle.addEventListener("pointermove", (e) => {
      if (!start) return;
      if (!moved && Math.hypot(e.clientX - start.x, e.clientY - start.y) < 5) return;
      if (!moved) { moved = true; el.classList.add("dragging"); if (el.dataset.open !== "true") cfg.setOpen(true); }
      cfg.applySize(cfg.clampSize(cfg.computeSize(e))); e.preventDefault();
    });
    const up = (e) => { if (!start) return; try { handle.releasePointerCapture(e.pointerId); } catch (_) {}
      if (moved) { el.classList.remove("dragging"); cfg.afterDrag && cfg.afterDrag(); } else cfg.setOpen(el.dataset.open !== "true");
      start = null; moved = false; };
    handle.addEventListener("pointerup", up); handle.addEventListener("pointercancel", up);
  }

  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false;
    placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
