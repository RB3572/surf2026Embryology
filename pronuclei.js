/* Pronuclei Distance vs Transcripts — zygotes with two auto-detected pronuclei.
 * Built on viewer-core.js (VCore). The 3-D view shows the two pronuclei and the
 * shortest line between them. The bottom drawer plots transcript count (y) against
 * pronuclei distance (x, a proxy for developmental time — the pronuclei migrate
 * together as the zygote ages) for (top) the SELECTED GENE and (below) all transcripts,
 * fit by a USER-SELECTABLE regression model (linear, quadratic, exponential, log, power,
 * logistic, Poisson/negative-binomial/binomial GLMs, or LOESS) — fitted client-side; the
 * R² is reported on each model's natural scale. The right drawer ranks genes by the
 * Pearson correlation of their count with the distance. Precompute in build_pronuclei.py.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const LINE_C = "#111827", CUR_C = "#0891b2", DOT_C = "#f59e0b";
  const MIN_ZYG = 5;                    // min zygotes containing a gene to rank/correlate it

  function plotInto(div, traces, layout, cfg) {
    if (!div.classList.contains("js-plotly-plot")) div.innerHTML = "";
    Plotly.react(div, traces, layout, cfg || { responsive: true, displayModeBar: false });
  }

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const pnReadout = $("#pn-readout"), pnFit = $("#pn-fit");
  const regTypeEl = $("#reg-type"), regNoteEl = $("#reg-note");
  const geneSelect = $("#gene-select"), geneScatter = $("#gene-scatter"), geneFit = $("#gene-fit");
  const dotsShow = $("#dots-show");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle"), drawerBody = $("#drawer-body");
  const scatterPlot = $("#scatter-plot");
  const rdrawer = $("#rdrawer"), rdrawerHandle = $("#rdrawer-handle");
  const rankNEl = $("#rank-n"), rankPosEl = $("#rank-pos"), rankNegEl = $("#rank-neg");
  const ptToggle = $("#pseudotime-toggle");

  const state = { points: [], byId: {}, genesAgg: null, geneCorr: [], userGene: null, rankN: 10,
                  currentId: null, scene: null, fit: null, drawerOpen: false, showDots: false,
                  regType: "linear", dotSize: 1.5, pseudotime: false };
  let vcExtras = null;   // dot-size + atlas-link row (VCore.addWindowExtras)

  (async function init() {
    try {
      const [m, ga] = await Promise.all([
        (await fetch("data/pronuclei_manifest.json")).json(),
        V.loadGz("data/pronuclei_genes.json.gz"),
      ]);
      state.points = m.embryos; state.genesAgg = ga;
      state.points.forEach((p) => (state.byId[p.id] = p));
      countEl.textContent = `${m.embryos.length} zygotes · pronuclei auto-detected inside the cytoplasm`;
      computeGeneCorr();
      populateGenes();
      V.buildTabs(tabsEl, m.embryos, selectEmbryo, (e) => ({
        label: e.label, sub: e.date_short,
        title: `${e.label} · dist ${e.distance} µm · ${e.total.toLocaleString()} transcripts`,
      }));
      wireDrawer(); wireRdrawer(); renderRanks();
      ptToggle.addEventListener("change", () => { state.pseudotime = ptToggle.checked; renderScatter(); renderGeneScatter(); });
      vcExtras = V.addWindowExtras($("#controls-body"), { defaultSize: state.dotSize, onDotSize: (s) => { state.dotSize = s; if (state.scene) render(); } });
      geneSelect.addEventListener("change", () => selectGene(geneSelect.value));
      rankNEl.addEventListener("change", () => { state.rankN = parseInt(rankNEl.value, 10) || 10; renderRanks(); });
      dotsShow.addEventListener("change", () => { state.showDots = dotsShow.checked; ensureDotGene(); if (state.scene) render(); });
      regTypeEl.addEventListener("change", () => { state.regType = regTypeEl.value; updRegNote(); renderScatter(); renderGeneScatter(); });
      updRegNote();
    } catch (err) { showError("Failed to load: " + (err.message || err)); }
  })();

  // ---------- gene ↔ distance correlations ----------
  function geneSeries(g) {
    const xs = [], ys = [], ids = [], labels = [];
    for (const e of state.genesAgg.embryos) {
      const c = e.genes[g]; if (c == null) continue;
      xs.push(c); ys.push(e.distance); ids.push(e.id); labels.push((state.byId[e.id] || {}).label || e.id);
    }
    return { xs, ys, ids, labels };
  }
  function computeGeneCorr() {
    const genes = new Set();
    state.genesAgg.embryos.forEach((e) => Object.keys(e.genes).forEach((g) => genes.add(g)));
    const out = [];
    for (const g of genes) {
      const s = geneSeries(g); if (s.xs.length < MIN_ZYG) continue;
      const f = linreg(s.xs, s.ys); if (!isFinite(f.r)) continue;
      out.push({ gene: g, r: f.r, n: s.xs.length });
    }
    state.geneCorr = out;
  }
  function populateGenes() {
    const genes = [...new Set(state.genesAgg.embryos.flatMap((e) => Object.keys(e.genes)))]
      .sort((a, b) => a.localeCompare(b, undefined, { sensitivity: "base", numeric: true }));
    geneSelect.innerHTML = genes.map((g) => `<option value="${g}">${g}</option>`).join("");
    const top = [...state.geneCorr].sort((a, b) => b.r - a.r)[0];   // default = strongest positive correlation
    state.userGene = (top ? top.gene : genes[0]) || "";
    geneSelect.value = state.userGene;
  }
  const gene = () => geneSelect.value;
  function selectGene(g) {
    state.userGene = g; geneSelect.value = g;
    renderGeneScatter(); highlightRank();
    if (state.scene) { renderReadout(state.byId[state.currentId] || {}); if (state.showDots) render(); }
  }
  // The gene dropdown is the UNION across all zygotes (disjoint MERFISH panels), so the
  // selected gene often has NO molecules in the current embryo — which silently showed no
  // dots. These keep "show dots" pointed at a gene that actually exists in this zygote.
  const geneHasDots = (g) => { const s = state.scene; return !!(s && s.transcripts && s.transcripts[g] && s.transcripts[g].x && s.transcripts[g].x.length); };
  function topEmbryoGene() {
    const s = state.scene; if (!s || !s.transcripts) return null;
    let best = null, bn = -1;
    for (const g in s.transcripts) { const t = s.transcripts[g], n = t && t.x ? t.x.length : 0; if (n > bn) { bn = n; best = g; } }
    return best;
  }
  function ensureDotGene() {   // when dots are on, snap the selection to a gene present in this zygote
    if (!state.showDots || !state.scene) return;
    if (geneHasDots(gene())) return;
    const g = topEmbryoGene();
    if (g && g !== gene()) { state.userGene = g; geneSelect.value = g; renderGeneScatter(); highlightRank(); }
  }

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    V.markActiveTab(tabsEl, id);
    const meta = state.byId[id] || {};
    showLoading(`Loading ${meta.label || id}…`);
    try {
      const scene = await V.loadGz(`data/pronuclei/${id}.json.gz`);
      if (state.currentId !== id) return;
      state.scene = scene;
      if (vcExtras) vcExtras.setAtlas(id);
      controlsEl.hidden = false; placeholder.hidden = true; drawer.hidden = false; rdrawer.hidden = false;
      ensureDotGene(); render(); renderReadout(meta);
      if (!state.drawerOpen) openDrawer(true);
      else { renderScatter(); renderGeneScatter(); }
    } catch (err) { showError(err.message || String(err)); }
    finally { hideLoading(); }
  }

  // ---------- 3-D ----------
  function segMesh(scene, lbl, color, opacity, name) {
    const mesh = scene.region_meshes[String(lbl)]; if (!mesh) return null;
    const v = mesh.verts, f = mesh.faces, nV = v.length / 3, nF = f.length / 3;
    const x = new Array(nV), y = new Array(nV), z = new Array(nV);
    for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
    const ii = new Array(nF), jj = new Array(nF), kk = new Array(nF);
    for (let i = 0; i < nF; i++) { ii[i] = f[i * 3]; jj[i] = f[i * 3 + 1]; kk[i] = f[i * 3 + 2]; }
    return { type: "mesh3d", x, y, z, i: ii, j: jj, k: kk, color, opacity,
      name, showlegend: true, flatshading: false, hoverinfo: "name",
      lighting: { ambient: 0.7, diffuse: 0.55, specular: 0.12, roughness: 0.9 }, legendrank: lbl };
  }
  function render() {
    const s = state.scene; if (!s) return;
    const [la, lb] = s.pron_labels;
    const pcolor = { [la]: "#2563eb", [lb]: "#dc2626" };
    const traces = [];
    for (const lbl of s.mask_labels) {
      const pron = lbl === la || lbl === lb;
      const t = segMesh(s, lbl, pron ? pcolor[lbl] : "#9aa3b2", pron ? 0.5 : 0.08,
        pron ? `Pronucleus (seg ${lbl})` : `Segment ${lbl}`);
      if (t) traces.push(t);
    }
    const [a, b] = s.line_plot;
    traces.push({ type: "scatter3d", mode: "lines+markers", name: `min dist ${s.distance_um} µm`,
      x: [a[0], b[0]], y: [a[1], b[1]], z: [a[2], b[2]],
      line: { color: LINE_C, width: 7 }, marker: { size: 4, color: LINE_C },
      hovertemplate: `min distance ${s.distance_um} µm<extra></extra>`, legendrank: 100 });
    if (state.showDots) {
      const sel = gene(), g = geneHasDots(sel) ? sel : topEmbryoGene();
      const tx = g ? s.transcripts[g] : null;
      if (tx && tx.x.length) {
        const zs = s.z_scale;
        const nm = g === sel ? `${g} · ${tx.x.length} dots`
                             : `${g} · ${tx.x.length} dots · (${sel} not in this zygote)`;
        traces.push({ type: "scatter3d", mode: "markers", name: nm,
          x: tx.x, y: tx.y, z: tx.gz.map((z) => z * zs),
          marker: { size: state.dotSize, color: DOT_C, opacity: 0.85, line: { width: 0 } },
          hovertemplate: `${g}<extra></extra>`, legendrank: 200 });
      }
    }
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }
  function renderReadout(meta) {
    const s = state.scene;
    const gc = (state.genesAgg.embryos.find((e) => e.id === s.id) || { genes: {} }).genes[gene()];
    pnReadout.innerHTML =
      `<div class="pn-big"><span>${s.distance_um}</span> µm <span class="pn-lbl">pronuclei distance</span></div>` +
      `<div class="pn-big"><span>${s.total_transcripts.toLocaleString()}</span> <span class="pn-lbl">total transcripts</span></div>` +
      `<div class="pn-resid"><b>${gene()}</b> here: ${gc != null ? gc.toLocaleString() + " transcripts" : "not in this zygote's panel"}</div>` +
      `<div class="pn-resid">pronuclei auto-detected as segments <b>${s.pron_labels[0]}</b> &amp; <b>${s.pron_labels[1]}</b></div>`;
  }

  // ---------- scatters ----------
  function linreg(xs, ys) {
    const n = xs.length; let sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
    for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; syy += ys[i] * ys[i]; }
    const den = n * sxx - sx * sx;
    const b = den ? (n * sxy - sx * sy) / den : 0;
    const a = (sy - b * sx) / n;
    const rden = Math.sqrt(den * (n * syy - sy * sy));
    const r = rden ? (n * sxy - sx * sy) / rden : 0;
    return { a, b, r, r2: r * r, n };
  }
  // ---------- regression models: transcript count (y) vs pronuclei distance (x) ----------
  // The two pronuclei migrate together as the zygote ages, so distance is a proxy for
  // developmental time (smaller = later) and transcript count is the response. Each model
  // returns a predictor y(x) + a fit statistic on its NATURAL scale (labelled in the UI).
  const avg = (a) => a.reduce((s, v) => s + v, 0) / (a.length || 1);
  const clmp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
  const sci = (v) => { const a = Math.abs(v); if (!isFinite(v)) return "–"; if (a === 0) return "0";
    return (a >= 1e4 || a < 1e-3) ? v.toExponential(2) : String(+v.toPrecision(3)); };
  function wls(xs, ys, ws) {                      // weighted least squares: y = b0 + b1·x
    let sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
    for (let i = 0; i < xs.length; i++) { const w = ws ? ws[i] : 1;
      sw += w; swx += w * xs[i]; swy += w * ys[i]; swxx += w * xs[i] * xs[i]; swxy += w * xs[i] * ys[i]; }
    const den = sw * swxx - swx * swx, b1 = den ? (sw * swxy - swx * swy) / den : 0;
    return [(swy - b1 * swx) / sw, b1];
  }
  function r2on(xs, ys, pred) {                   // R² of pred() vs ys over xs
    const yb = avg(ys); let sr = 0, st = 0;
    for (let i = 0; i < ys.length; i++) { const e = ys[i] - pred(xs[i]); sr += e * e; st += (ys[i] - yb) ** 2; }
    return st > 0 ? 1 - sr / st : 0;
  }
  const r2lin = (xs, ys) => { const [a, b] = wls(xs, ys); return r2on(xs, ys, (x) => a + b * x); };
  function solve3(A, d) {                         // 3×3 Gaussian elimination
    A = A.map((r, i) => r.concat(d[i]));
    for (let c = 0; c < 3; c++) { let p = c;
      for (let r = c + 1; r < 3; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
      [A[c], A[p]] = [A[p], A[c]]; if (Math.abs(A[c][c]) < 1e-12) A[c][c] = 1e-12;
      for (let r = 0; r < 3; r++) { if (r === c) continue; const f = A[r][c] / A[c][c]; for (let k = c; k < 4; k++) A[r][k] -= f * A[c][k]; } }
    return [A[0][3] / A[0][0], A[1][3] / A[1][1], A[2][3] / A[2][2]];
  }
  function quadFit(xs, ys) {                       // y = c0 + c1·x + c2·x²
    const S = [0, 0, 0, 0, 0], T = [0, 0, 0];
    for (let i = 0; i < xs.length; i++) { let xp = 1; for (let k = 0; k < 5; k++) { S[k] += xp; if (k < 3) T[k] += xp * ys[i]; xp *= xs[i]; } }
    return solve3([[S[0], S[1], S[2]], [S[1], S[2], S[3]], [S[2], S[3], S[4]]], T);
  }
  function irls(xs, ys, fam, init, iter = 60) {   // GLM via iteratively-reweighted least squares
    let [b0, b1] = init;
    for (let k = 0; k < iter; k++) {
      const z = [], w = [];
      for (let i = 0; i < xs.length; i++) { const e = b0 + b1 * xs[i], mu = fam.li(e), g = fam.me(e, mu), V = fam.v(mu);
        w.push(g * g / Math.max(V, 1e-9)); z.push(e + (ys[i] - mu) / (g || 1e-9)); }
      const [n0, n1] = wls(xs, z, w);
      if (!isFinite(n0) || !isFinite(n1)) break;
      if (Math.abs(n0 - b0) + Math.abs(n1 - b1) < 1e-10) { b0 = n0; b1 = n1; break; }
      b0 = n0; b1 = n1;
    }
    return [b0, b1];
  }
  const POIS = { li: (e) => Math.exp(Math.min(e, 30)), me: (e, m) => m, v: (m) => Math.max(m, 1e-9) };
  const NBfam = (t) => ({ li: (e) => Math.exp(Math.min(e, 30)), me: (e, m) => m, v: (m) => m + m * m / t });
  const BINfam = (N) => ({ li: (e) => N / (1 + Math.exp(-clmp(e, -30, 30))), me: (e, m) => { const p = m / N; return N * p * (1 - p); }, v: (m) => { const p = m / N; return Math.max(N * p * (1 - p), 1e-9); } });
  const logInit = (xs, ys) => wls(xs, ys.map((y) => Math.log(Math.max(y, 1))));
  function devR2(xs, ys, mu) {                     // Poisson-deviance pseudo-R²
    const yb = avg(ys), d = (y, m) => 2 * ((y > 0 ? y * Math.log(y / m) : 0) - (y - m));
    let dr = 0, dn = 0; for (let i = 0; i < ys.length; i++) { dr += d(ys[i], mu[i]); dn += d(ys[i], yb); }
    return dn > 0 ? 1 - dr / dn : 0;
  }
  function binDevR2(ys, mu, N) {                    // binomial-deviance pseudo-R²
    const yb = avg(ys);
    const d = (y, m) => 2 * ((y > 0 ? y * Math.log(y / m) : 0) + (N - y > 0 ? (N - y) * Math.log((N - y) / (N - m)) : 0));
    let dr = 0, dn = 0; for (let i = 0; i < ys.length; i++) { dr += d(ys[i], mu[i]); dn += d(ys[i], yb); }
    return dn > 0 ? 1 - dr / dn : 0;
  }
  function loessPredictor(xs, ys, span) {          // local linear regression (tricube weights)
    const n = xs.length, k = Math.max(3, Math.round(span * n));
    return (x0) => {
      const sorted = xs.map((x) => Math.abs(x - x0)).sort((a, b) => a - b);
      const h = sorted[Math.min(n - 1, k - 1)] || 1e-9;
      const wx = [], wy = [], ww = [];
      for (let i = 0; i < n; i++) { const u = Math.abs(xs[i] - x0) / h; if (u >= 1) continue; wx.push(xs[i]); wy.push(ys[i]); ww.push((1 - u ** 3) ** 3); }
      if (wx.length < 2) return avg(ys);
      const [a, b] = wls(wx, wy, ww); return a + b * x0;
    };
  }
  const MODELS = {
    linear: { label: "Linear", scale: "raw", bio: "constant rate of change with developmental time — the simplest baseline.",
      fit(xs, ys) { const [a, b] = wls(xs, ys), p = (x) => a + b * x;
        return { predict: p, r2: r2on(xs, ys, p), params: `y = ${sci(a)} ${b >= 0 ? "+" : "−"} ${sci(Math.abs(b))}·x` }; } },
    quadratic: { label: "Quadratic", scale: "raw", bio: "one peak or trough — captures the maternal-to-zygotic hand-off (maternal store depletes, then zygotic transcription rises).",
      fit(xs, ys) { const [c0, c1, c2] = quadFit(xs, ys), p = (x) => c0 + c1 * x + c2 * x * x;
        return { predict: p, r2: r2on(xs, ys, p), params: `y = ${sci(c0)} + ${sci(c1)}·x + ${sci(c2)}·x²` }; } },
    exp: { label: "Exponential", scale: "log", bio: "constant fractional change per unit time — first-order kinetics of maternal-mRNA decay (or exponential accumulation).",
      fit(xs, ys) { const px = [], py = []; for (let i = 0; i < xs.length; i++) if (ys[i] > 0) { px.push(xs[i]); py.push(Math.log(ys[i])); }
        const [la, b] = wls(px, py), a = Math.exp(la);
        return { predict: (x) => a * Math.exp(b * x), r2: r2lin(px, py), params: `y = ${sci(a)}·e^(${b.toFixed(4)}·x)` }; } },
    log: { label: "Logarithmic", scale: "raw", bio: "fast early change that levels off — diminishing returns.",
      fit(xs, ys) { const lx = xs.map((x) => Math.log(Math.max(x, 1e-6))); const [a, b] = wls(lx, ys), p = (x) => a + b * Math.log(Math.max(x, 1e-6));
        return { predict: p, r2: r2on(xs, ys, p), params: `y = ${sci(a)} ${b >= 0 ? "+" : "−"} ${sci(Math.abs(b))}·ln(x)` }; } },
    power: { label: "Power law", scale: "log-log", bio: "scale-free (allometric) — a fixed % change in count per % change in the time proxy.",
      fit(xs, ys) { const px = [], py = []; for (let i = 0; i < xs.length; i++) if (xs[i] > 0 && ys[i] > 0) { px.push(Math.log(xs[i])); py.push(Math.log(ys[i])); }
        const [la, b] = wls(px, py), a = Math.exp(la);
        return { predict: (x) => a * Math.pow(Math.max(x, 1e-6), b), r2: r2lin(px, py), params: `y = ${sci(a)}·x^${b.toFixed(3)}` }; } },
    logistic: { label: "Logistic (sigmoid)", scale: "logit", bio: "a saturating switch — the shape of zygotic genome activation (transcription turns on, then plateaus).",
      fit(xs, ys) { const L = 1.02 * Math.max(...ys); const px = [], py = []; for (let i = 0; i < xs.length; i++) { const q = clmp(ys[i] / L, 0.001, 0.999); px.push(xs[i]); py.push(Math.log(q / (1 - q))); }
        const [a, b] = wls(px, py);
        return { predict: (x) => L / (1 + Math.exp(-clmp(a + b * x, -30, 30))), r2: r2lin(px, py), params: `L=${sci(L)}, k=${b.toFixed(3)}, x₀=${b ? (-a / b).toFixed(1) : "–"} µm` }; } },
    poisson: { label: "Poisson (GLM)", scale: "deviance", bio: "the canonical model for count data — a log-linear rate with variance equal to the mean.",
      fit(xs, ys) { const [b0, b1] = irls(xs, ys, POIS, logInit(xs, ys)); const mu = xs.map((x) => Math.exp(b0 + b1 * x));
        return { predict: (x) => Math.exp(b0 + b1 * x), r2: devR2(xs, ys, mu), params: `log μ = ${b0.toFixed(2)} ${b1 >= 0 ? "+" : "−"} ${Math.abs(b1).toFixed(4)}·x` }; } },
    negbin: { label: "Negative binomial (GLM)", scale: "deviance", bio: "count model for OVER-dispersed data (variance ≫ mean) — the standard for transcript counts (as in DESeq2 / edgeR).",
      fit(xs, ys) { const [p0, p1] = irls(xs, ys, POIS, logInit(xs, ys)); const m0 = xs.map((x) => Math.exp(p0 + p1 * x));
        let nu = 0, de = 0; for (let i = 0; i < ys.length; i++) { nu += m0[i] * m0[i]; de += Math.max((ys[i] - m0[i]) ** 2 - m0[i], 0); } const th = de > 0 ? nu / de : 1e6;
        const [b0, b1] = irls(xs, ys, NBfam(th), [p0, p1]); const mu = xs.map((x) => Math.exp(b0 + b1 * x));
        return { predict: (x) => Math.exp(b0 + b1 * x), r2: devR2(xs, ys, mu), params: `log μ = ${b0.toFixed(2)} ${b1 >= 0 ? "+" : "−"} ${Math.abs(b1).toFixed(4)}·x · θ=${sci(th)}` }; } },
    binomial: { label: "Binomial (GLM, logit)", scale: "deviance", bio: "models the count as a fraction of a ceiling via a logit link — a saturating S-curve on the count scale.",
      fit(xs, ys) { const N = Math.round(1.02 * Math.max(...ys)); const init = wls(xs, ys.map((y) => Math.log((y + 0.5) / (N - y + 0.5))));
        const [b0, b1] = irls(xs, ys, BINfam(N), init); const pred = (x) => N / (1 + Math.exp(-clmp(b0 + b1 * x, -30, 30)));
        return { predict: pred, r2: binDevR2(ys, xs.map(pred), N), params: `N=${sci(N)}, logit p = ${b0.toFixed(2)} ${b1 >= 0 ? "+" : "−"} ${Math.abs(b1).toFixed(4)}·x` }; } },
    loess: { label: "LOESS (local smoother)", scale: "raw", bio: "non-parametric — lets the data show its own trend with no assumed functional form.",
      fit(xs, ys) { const p = loessPredictor(xs, ys, 0.6); return { predict: p, r2: r2on(xs, ys, p), params: "local linear · span 0.6" }; } },
  };
  const SCALE_LABEL = { raw: "R²", log: "R² (log)", "log-log": "R² (log-log)", logit: "R² (logit)", deviance: "pseudo-R²" };
  function fitModel(type, xs, ys) {
    const m = MODELS[type] || MODELS.linear;
    let res; try { res = m.fit(xs.slice(), ys.slice()); } catch (_) { res = MODELS.linear.fit(xs.slice(), ys.slice()); }
    return { ...res, type, label: m.label, bio: m.bio, scale: m.scale, n: xs.length };
  }

  // ---------- scatters (x = pronuclei distance µm, y = transcript count) ----------
  function scatter(div, xs, ys, ids, labels, model, xTitle, yTitle, yUnit, curId, xUnit) {
    xUnit = xUnit == null ? "µm" : xUnit;
    const pts = xs.map((x, i) => ({ x, y: ys[i], lab: labels[i], id: ids[i] }));
    const others = pts.filter((o) => o.id !== curId), cur = pts.find((o) => o.id === curId);
    const xmin = Math.min(...xs), xmax = Math.max(...xs);
    const cx = [], cy = [], N = 90;                // smooth model curve
    for (let i = 0; i < N; i++) { const x = xmin + (xmax - xmin) * i / (N - 1), y = model.predict(x);
      if (isFinite(y)) { cx.push(x); cy.push(y); } }
    const hover = `%{text}<br>%{x:.1f}${xUnit ? " " + xUnit : ""} · %{y:,} ${yUnit}<extra></extra>`;
    const traces = [
      { type: "scatter", mode: "markers", name: "zygotes", x: others.map((o) => o.x), y: others.map((o) => o.y),
        marker: { size: 8, color: "#94a3b8", opacity: 0.72, line: { width: 0 } },
        text: others.map((o) => o.lab), hovertemplate: hover },
      { type: "scatter", mode: "lines", name: model.label, x: cx, y: cy,
        line: { color: "#0891b2", width: 2.4, shape: "spline" }, hoverinfo: "skip" },
    ];
    if (cur) traces.push({ type: "scatter", mode: "markers", name: cur.lab, x: [cur.x], y: [cur.y],
      marker: { size: 15, color: CUR_C, line: { width: 2, color: "#fff" } },
      text: [cur.lab], hovertemplate: hover });
    plotInto(div, traces, {
      margin: { l: 64, r: 12, t: 6, b: 40 }, height: div.clientHeight || 220,
      xaxis: { title: { text: xTitle, font: { size: 11 } }, tickfont: { size: 10 }, gridcolor: "#eef1f5", zeroline: false },
      yaxis: { title: { text: yTitle, font: { size: 10 } }, tickfont: { size: 9 }, gridcolor: "#eef1f5", rangemode: "tozero" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      legend: { orientation: "h", font: { size: 10 }, y: 1.16, x: 1, xanchor: "right" },
    });
  }
  function updRegNote() { const m = MODELS[state.regType] || MODELS.linear; if (regNoteEl) regNoteEl.textContent = m.bio; }
  // pseudotime = max(distance among the plotted embryos) − distance, so larger = later in
  // development (pronuclei approach over time). maxD is per-chart (the embryos being plotted).
  const ptx = (dists, maxD) => (state.pseudotime ? dists.map((d) => maxD - d) : dists);
  const PT_X_TITLE = "pseudotime  ·  max distance − pronuclei distance (larger = later)";
  function renderScatter() {
    const pts = state.points, dists = pts.map((p) => p.distance), Y = pts.map((p) => p.total);
    const X = ptx(dists, Math.max(...dists));
    state.fit = fitModel(state.regType, X, Y);
    scatter(scatterPlot, X, Y, pts.map((p) => p.id), pts.map((p) => p.label), state.fit,
      state.pseudotime ? PT_X_TITLE : "min pronuclei distance (µm)  ·  smaller = later in development",
      "total transcripts", "transcripts", state.currentId, state.pseudotime ? "" : "µm");
    const f = state.fit;
    pnFit.innerHTML = `· <b>${f.n}</b> zygotes · <b>${f.label}</b> · ${SCALE_LABEL[f.scale]} = <b>${f.r2.toFixed(3)}</b> · <span class="pn-params">${f.params}</span>`;
  }
  function renderGeneScatter() {
    const g = gene(), s = geneSeries(g);
    if (s.xs.length < 2) {
      Plotly.purge(geneScatter); geneScatter.classList.remove("js-plotly-plot");
      geneScatter.innerHTML = `<div class="pn-empty">Only ${s.xs.length} zygote${s.xs.length === 1 ? "" : "s"} contain${s.xs.length === 1 ? "s" : ""} <b>${g}</b> — too few to fit.</div>`;
      geneFit.innerHTML = `· <b>${g}</b>`;
      return;
    }
    const X = ptx(s.ys, Math.max(...s.ys)), Y = s.xs;   // x = distance (or pseudotime), y = this gene's count
    const model = fitModel(state.regType, X, Y);
    scatter(geneScatter, X, Y, s.ids, s.labels, model,
      state.pseudotime ? "pseudotime" : "min pronuclei distance (µm)", `${g} transcript count`, g, state.currentId,
      state.pseudotime ? "" : "µm");
    geneFit.innerHTML = `· <b>${g}</b> · ${model.label} · ${SCALE_LABEL[model.scale]} = <b>${model.r2.toFixed(3)}</b> · n=${s.xs.length}`;
  }

  // ---------- right drawer: correlation ranking ----------
  function rankRows(rows) {
    const cur = gene();
    let html = `<div class="best-head"><span></span><span>gene</span><span>r</span><span>n</span></div>`;
    html += rows.map((r, i) =>
      `<div class="best-row pn-row${r.gene === cur ? " current" : ""}" data-gene="${r.gene}" ` +
      `title="Pearson r = ${r.r.toFixed(3)} of ${r.gene} count vs pronuclei distance across ${r.n} zygotes">` +
      `<span class="best-num">${i + 1}</span><span class="best-gene">${r.gene}</span>` +
      `<span class="best-real" style="color:${r.r >= 0 ? "#dc2626" : "#2563eb"}">${r.r >= 0 ? "+" : ""}${r.r.toFixed(2)}</span>` +
      `<span class="best-p">${r.n}</span></div>`).join("");
    return html || `<div class="pn-empty">No genes in ≥ ${MIN_ZYG} zygotes.</div>`;
  }
  function renderRanks() {
    const n = state.rankN;
    rankPosEl.innerHTML = rankRows([...state.geneCorr].sort((a, b) => b.r - a.r).slice(0, n));
    rankNegEl.innerHTML = rankRows([...state.geneCorr].sort((a, b) => a.r - b.r).slice(0, n));
  }
  function highlightRank() {
    const cur = gene();
    rdrawer.querySelectorAll(".pn-row").forEach((r) => r.classList.toggle("current", r.dataset.gene === cur));
  }

  // ---------- drawers ----------
  const resizeScatters = () => { try { Plotly.Plots.resize(scatterPlot); Plotly.Plots.resize(geneScatter); } catch (_) {} };
  function openDrawer(open) {
    state.drawerOpen = open;
    drawer.dataset.open = open ? "true" : "false";
    drawerHandle.setAttribute("aria-expanded", String(open));
    if (open) { renderScatter(); renderGeneScatter(); requestAnimationFrame(resizeScatters); }
  }
  function wireDrawer() {
    wireHandleDrag(drawer, drawerHandle, {
      computeSize: (e) => window.innerHeight - e.clientY - 40,
      clampSize: (px) => Math.max(240, Math.min(window.innerHeight - 100, px)),
      applySize: (px) => drawer.style.setProperty("--drawer-h", px + "px"),
      setOpen: openDrawer, afterDrag: resizeScatters,
    });
    let sh = 0; const rz = $("#drawer-resize");
    rz.addEventListener("pointerdown", (e) => { sh = drawerBody.getBoundingClientRect().height; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return;
      drawer.style.setProperty("--drawer-h", Math.max(240, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y))) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} resizeScatters(); } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }
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
    rdrawer.addEventListener("click", (e) => {
      const row = e.target.closest(".pn-row"); if (row) selectGene(row.dataset.gene);
    });
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
