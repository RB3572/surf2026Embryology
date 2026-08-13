/* Pseudosperm Division Plane.
 *
 * The Sperm Division Plane project cuts a zygote with the plane through {sperm, polar-body COM,
 * cell COM}. That plane contains the polar axis, so it is one member of the family the Zygote
 * Division Planes project sweeps — the sperm just picks the rotation angle. 30 zygotes have a
 * sperm and therefore an angle. 20 do not.
 *
 * This page runs that backwards. The right drawer ranks genes by how sharply they split across the
 * REAL sperm plane, pooled over every sperm zygote that measures them. Select some, and for a
 * zygote with no sperm the page finds the angle that best reproduces that split — the PSEUDOSPERM
 * plane — and draws where the sperm would have had to be.
 *
 * ⚠️ A pseudosperm plane is FITTED, not observed. The Validation tab exists to keep that honest:
 * it runs the same fit on the zygotes that do have a sperm, with the sperm hidden, and shows how
 * far off it lands. Chance is 45°, not 0°.
 *
 * Data: data/pseudosperm.json.gz (build_pseudosperm.py) + data/zygote/<id>.json.gz for the scene.
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const C_PLANE = "#7c3aed", C_AXIS = "#0f172a", C_SPERM = "#16a34a", C_PSEUDO = "#0891b2";
  const C_A = "#b45309", C_B = "#0d9488";

  const state = {
    data: null, byId: {}, rank: [], curId: null, scene: null,
    sel: new Set(), metric: "vol", sort: "fisher", minM: 5, tab: "cross",
    show: { plane: true, axis: true, sperm: true, tx: false, mesh: true },
  };
  const meta = () => state.data.meta;
  const cur = () => state.byId[state.curId];
  const K = () => meta().grid.n;

  // ───────── log-space statistics, mirroring build_pseudosperm.py ─────────
  const LG_C = [676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
                12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
                1.5056327351493116e-7];
  function lgamma(z) {
    if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
    z -= 1;
    let x = 0.99999999999980993;
    for (let i = 0; i < 8; i++) x += LG_C[i] / (z + i + 1);
    const t = z + 7.5;
    return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
  }
  /** log of an exact binomial tail, summed in log space — the real p reaches 1e-3000 here, so a
   *  plain probability would underflow to 0 and tie every strong gene together. */
  function binomLogTail(a, n, p0, upper) {
    if (p0 <= 0 || p0 >= 1) return 0;
    const lgn = lgamma(n + 1), lp = Math.log(p0), lq = Math.log1p(-p0);
    const step = upper ? 1 : -1, stop = upper ? n : 0;
    const terms = [];
    for (let k = a; upper ? k <= stop : k >= stop; k += step) {
      const t = lgn - lgamma(k + 1) - lgamma(n - k + 1) + k * lp + (n - k) * lq;
      terms.push(t);
      if (terms.length > 4 && t < terms[0] - 40) break;
    }
    if (!terms.length) return -Infinity;
    const m = Math.max(...terms);
    return m + Math.log(terms.reduce((s, t) => s + Math.exp(t - m), 0));
  }
  const binomLogP = (a, n, p0) =>
    Math.min(Math.LN2 + Math.min(binomLogTail(a, n, p0, false), binomLogTail(a, n, p0, true)), 0);
  const log10p = (a, n, p0) => (n > 0 ? binomLogP(a, n, p0) / Math.LN10 : 0);
  function fmtP(l10) {
    if (!isFinite(l10)) return "0";
    if (l10 > -3) return Math.pow(10, l10).toFixed(3);
    const e = Math.floor(l10), m = Math.pow(10, l10 - e);
    return `${m.toFixed(1)}e${e}`;
  }

  // ───────── per-embryo angle sweep ─────────
  /** Side-A counts per angle for one gene, decoded from the delta encoding. */
  function counts(g) {
    if (g._a) return g._a;
    const out = new Int32Array(g.a.length);
    let s = 0;
    for (let i = 0; i < g.a.length; i++) { s += g.a[i]; out[i] = s; }
    g._a = out;
    return out;
  }
  /** Pooled side-A counts over the selected genes, at every angle. */
  function pooled(e, genes) {
    const k = K(), A = new Float64Array(k);
    let n = 0;
    for (const g of genes) {
      const c = counts(g);
      for (let i = 0; i < k; i++) A[i] += c[i];
      n += g.n;
    }
    return { A, n };
  }
  /** The asymmetry objective at every angle, for the current measure. */
  function sweep(e, genes) {
    const { A, n } = pooled(e, genes);
    const k = K(), out = new Float64Array(k);
    if (!n) return { curve: out, n: 0, best: 0 };
    for (let i = 0; i < k; i++) {
      const a = A[i], b = n - a;
      out[i] = state.metric === "vol"
        ? Math.abs(a / e.volA[i] - b / e.volB[i]) * (e.vol_total / 2)   // per half-cell volume
        : Math.abs(a - b) / n;
    }
    let best = 0;
    for (let i = 1; i < k; i++) if (out[i] > out[best]) best = i;
    return { curve: out, n, best, A };
  }
  const genesOf = (e, names) => names.map((g) => e.genes.find((r) => r.g === g)).filter(Boolean);
  const selGenes = (e) => genesOf(e, [...state.sel]);

  /** The plane this embryo is currently showing: the real sperm plane if it has one, otherwise
   *  the fitted pseudosperm plane. Returns null when nothing is selected and there is no sperm. */
  function activePlane(e) {
    if (e.sperm) {
      return { kind: "sperm", angle: e.sperm.angle_deg, k: Math.round(e.sperm.angle_deg) % K(),
               volA: e.sperm.volA, volB: e.sperm.volB, exact: true };
    }
    const gs = selGenes(e);
    if (!gs.length) return null;
    const s = sweep(e, gs);
    return { kind: "pseudo", angle: s.best * meta().grid.step_deg, k: s.best,
             volA: e.volA[s.best], volB: e.volB[s.best], exact: false, sweepRes: s };
  }

  /** Per-gene numbers at a given angle index (or at the exact sperm plane). */
  function geneStats(e, pl) {
    const out = [];
    for (const name of state.sel) {
      const g = e.genes.find((r) => r.g === name);
      if (!g) continue;
      let a, n;
      if (pl.kind === "sperm" && e.sperm.n[name] !== undefined) {
        a = e.sperm.a[name]; n = e.sperm.n[name];              // exact, from the sperm project
      } else {
        a = counts(g)[pl.k]; n = g.n;
      }
      const b = n - a, p0 = pl.volA / (pl.volA + pl.volB);
      const cA = a / pl.volA, cB = b / pl.volB;
      out.push({ g: name, a, b, n, cA, cB,
                 l2: Math.log2((cA + 1e-12) / (cB + 1e-12)),
                 log10p: log10p(a, n, p0) });
    }
    return out;
  }

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await V.loadGz("data/pseudosperm.json.gz"); }
    catch (err) {
      $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_pseudosperm.py</code>.</div></div>`;
      return;
    }
    const m = meta();
    state.data.embryos.forEach((e) => (state.byId[e.id] = e));
    state.rank = state.data.ranking;
    $("#embryo-count").textContent =
      `${m.n_embryos} zygotes · ${m.n_sperm} with a sperm plane · ${m.n_pseudo} pseudosperm · ` +
      `${m.n_genes} genes ranked`;
    $("#ps-rank-desc").innerHTML =
      `Each gene tested across the <b>real</b> sperm plane in every one of the ${m.n_sperm} zygotes
       that has one and measures it, then combined. Tick genes to fit the pseudosperm plane.`;

    V.buildTabs($("#tabs"), state.data.embryos, selectEmbryo, (e) => ({
      label: e.label, sub: e.date, cls: e.sperm ? "ps-has-sperm" : "",
      title: `${e.id}${e.sperm ? " · sperm plane @ " + e.sperm.angle_deg.toFixed(1) + "°" : " · no sperm — pseudosperm plane"}`,
    }));
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    selectTop(5);                       // seed the selection BEFORE the list is drawn
    renderRank();
    selectEmbryo(state.data.embryos[0].id);
  })();

  // ───────── the ranking ─────────
  function rankRows() {
    const rows = state.rank.filter((r) => r.m >= state.minM);
    const key = {
      fisher: (r) => r.fisherLog10P,
      stouffer: (r) => -Math.abs(r.stoufferZ),
      fold: (r) => -r.absL,
      m: (r) => -r.m,
    }[state.sort];
    return rows.slice().sort((a, b) => key(a) - key(b));
  }
  function renderRank() {
    const rows = rankRows(), el = $("#ps-rank");
    if (!rows.length) {
      el.innerHTML = `<div class="ps-empty">No gene is measured in ${state.minM}+ sperm zygotes.
        The probesets are disjoint, so no gene is in every embryo — lower the threshold.</div>`;
      return;
    }
    el.innerHTML =
      `<div class="ps-rank-head"><span></span><span>gene</span><span>zygotes</span>
         <span>log₁₀p</span><span>|log₂ fold|</span><span>Z</span></div>` +
      rows.map((r, i) => `<div class="ps-rank-row${state.sel.has(r.g) ? " on" : ""}" data-g="${r.g}"
          title="${r.g} · ${r.n.toLocaleString()} transcripts in ${r.m} zygotes · ${r.nPos} on side A, ${r.nNeg} on side B">
        <span class="c"><input type="checkbox" ${state.sel.has(r.g) ? "checked" : ""} tabindex="-1" /></span>
        <span class="e">${r.g}</span>
        <span class="m">${r.m}</span>
        <span class="p">${r.fisherLog10P.toFixed(0)}</span>
        <span class="f">${r.absL.toFixed(2)}</span>
        <span class="z">${r.stoufferZ.toFixed(1)}</span></div>`).join("");
    el.querySelectorAll(".ps-rank-row").forEach((row) => row.addEventListener("click", () => {
      const g = row.dataset.g;
      if (state.sel.has(g)) state.sel.delete(g); else state.sel.add(g);
      syncSel(); renderRank(); refresh();
    }));
  }
  function syncSel() {
    const n = state.sel.size;
    $("#ps-sel-count").textContent = n ? `${n} gene${n > 1 ? "s" : ""} selected` : "nothing selected";
  }
  function selectTop(k) {
    state.sel = new Set(rankRows().slice(0, k).map((r) => r.g));
    syncSel();
  }

  // ───────── 3-D ─────────
  async function selectEmbryo(id) {
    if (id === state.curId) return;
    state.curId = id;
    V.markActiveTab($("#tabs"), id);
    const e = cur();
    $("#loading").hidden = false; $("#loading-text").textContent = `Loading ${e.label}…`;
    try { state.scene = await V.loadGz(`data/zygote/${e.scene}`); }
    catch (err) { state.scene = null; }
    $("#loading").hidden = true;
    $("#controls").hidden = false; $("#placeholder").hidden = true;
    $("#drawer-emb").textContent = e.label;
    refresh();
  }

  const addv = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const mulv = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
  const norm3 = (p) => { const n = Math.hypot(p[0], p[1], p[2]) || 1; return [p[0] / n, p[1] / n, p[2] / n]; };
  /* Plot space scales x,y and z DIFFERENTLY (x,y = µm/0.15, z = µm·z_scale), so a direction
   * expressed there is not a unit vector and cannot be scaled by a length. Every construction
   * below is therefore done in µm and converted per corner at the end, exactly as zygote.js does. */
  const umOf = (p, zs) => norm3([p[0] * XY, p[1] * XY, p[2] / zs]);
  const toPlot = (p, zs) => [p[0] / XY, p[1] / XY, p[2] * zs];

  const PLANE_SCALE = 2;                    // rendered plane size (×) vs the precomputed L

  /** The plane at angle index k. Its normal is n(θ) = cosθ·u + sinθ·v, so the plane itself is
   *  spanned by the polar axis and m(θ) = −sinθ·u + cosθ·v — the direction a sperm defining it
   *  would have to lie along. All in µm; the quad corners are converted to plot space. */
  function planeGeom(e, k) {
    const zs = state.scene ? state.scene.z_scale : 1;
    const th = k * meta().grid.step_deg * Math.PI / 180;
    const u = umOf(e.u_plot, zs), v = umOf(e.v_plot, zs), a = umOf(e.axis_plot, zs);
    const C = e.com_um;
    const m = norm3(addv(mulv(u, -Math.sin(th)), mulv(v, Math.cos(th))));
    const nrm = norm3(addv(mulv(u, Math.cos(th)), mulv(v, Math.sin(th))));
    const L = e.L * PLANE_SCALE;
    const quad = [[-1, -1], [-1, 1], [1, 1], [1, -1]].map(([sa, sm]) =>
      toPlot(addv(C, addv(mulv(a, sa * L), mulv(m, sm * L))), zs));
    return { quad, m, nrm, a, C, L, zs };
  }

  function render3D() {
    const e = cur(), host = $("#plot");
    if (!e || !state.scene) return;
    const sc = state.scene, zs = sc.z_scale;
    const pl = activePlane(e);
    const traces = [];
    if (state.show.mesh) traces.push(...V.bodyTraces(sc));

    if (pl && state.show.plane) {
      const G = planeGeom(e, pl.k);
      const q = G.quad;
      traces.push({
        type: "mesh3d", x: q.map((p) => p[0]), y: q.map((p) => p[1]), z: q.map((p) => p[2]),
        i: [0, 0], j: [1, 2], k: [2, 3], color: pl.kind === "sperm" ? C_PLANE : C_PSEUDO,
        opacity: 0.3, name: pl.kind === "sperm" ? "sperm plane" : "pseudosperm plane",
        showlegend: true, hoverinfo: "skip", flatshading: true,
      });
    }
    if (state.show.axis) {
      const a = umOf(e.axis_plot, zs), C = e.com_um, L = e.L * PLANE_SCALE;
      const p0 = toPlot(addv(C, mulv(a, -L)), zs), p1 = toPlot(addv(C, mulv(a, L)), zs);
      traces.push({ type: "scatter3d", mode: "lines", x: [p0[0], p1[0]], y: [p0[1], p1[1]],
        z: [p0[2], p1[2]], line: { color: C_AXIS, width: 5, dash: "dot" },
        name: "polar-body axis", hoverinfo: "skip" });
      traces.push({ type: "scatter3d", mode: "markers", x: [e.pb_plot[0]], y: [e.pb_plot[1]],
        z: [e.pb_plot[2]], marker: { size: 6, color: C_AXIS, symbol: "circle" },
        name: "polar body", hovertemplate: "polar body<extra></extra>" });
    }
    if (state.show.sperm && pl) {
      if (e.sperm) {
        const s = e.sperm.sperm_plot;
        traces.push({ type: "scatter3d", mode: "markers", x: [s[0]], y: [s[1]], z: [s[2]],
          marker: { size: 7, color: C_SPERM, symbol: "diamond",
                    line: { color: "#fff", width: 1 } },
          name: "sperm", hovertemplate: "sperm<extra></extra>" });
      } else {
        // the direction a sperm would have to lie along for this to BE its plane. Drawn both
        // ways: the plane is a plane, so it fixes the line, never which end of it.
        const G = planeGeom(e, pl.k);
        const t1 = toPlot(addv(G.C, mulv(G.m, e.L)), zs);
        const t2 = toPlot(addv(G.C, mulv(G.m, -e.L)), zs);
        const c0 = toPlot(G.C, zs);
        traces.push({ type: "scatter3d", mode: "lines", x: [t2[0], c0[0], t1[0]],
          y: [t2[1], c0[1], t1[1]], z: [t2[2], c0[2], t1[2]],
          line: { color: C_PSEUDO, width: 5 }, name: "pseudosperm direction",
          hovertemplate: "a sperm would have to lie along this line<extra></extra>" });
        traces.push({ type: "scatter3d", mode: "markers", x: [t1[0], t2[0]], y: [t1[1], t2[1]],
          z: [t1[2], t2[2]], marker: { size: 7, color: C_PSEUDO, symbol: "diamond-open",
            line: { width: 2 } }, name: "either end", showlegend: false,
          hovertemplate: "either end is equally consistent<extra></extra>" });
      }
    }
    if (state.show.tx && pl) {
      const G = planeGeom(e, pl.k);
      for (const name of state.sel) {
        const t = sc.transcripts[name];
        if (!t) continue;
        const xs = [[], []], ys = [[], []], zs2 = [[], []];
        for (let i = 0; i < t.x.length; i++) {
          if (!t.s1[i]) continue;                         // cytoplasm only, never pronuclei/PB
          // the normal is a µm unit vector, so the offset has to be in µm too
          const d = [t.x[i] * XY - e.com_um[0], t.y[i] * XY - e.com_um[1], t.gz[i] - e.com_um[2]];
          const s = d[0] * G.nrm[0] + d[1] * G.nrm[1] + d[2] * G.nrm[2] > 0 ? 0 : 1;
          xs[s].push(t.x[i]); ys[s].push(t.y[i]); zs2[s].push(t.gz[i] * zs);
        }
        [0, 1].forEach((s) => xs[s].length && traces.push({
          type: "scatter3d", mode: "markers", x: xs[s], y: ys[s], z: zs2[s],
          marker: { size: 1.7, color: s ? C_B : C_A, opacity: 0.6 },
          name: `${name} · side ${s ? "B" : "A"}`, hoverinfo: "skip" }));
      }
    }
    Plotly.react(host, traces, V.sceneLayout(sc.extents, sc.id), V.plotConfig);
  }

  // ───────── readout ─────────
  function renderReadout() {
    const e = cur(); if (!e) return;
    const pl = activePlane(e);
    const head = $("#ps-plane-head");
    if (!pl) {
      head.innerHTML = `<div class="ps-badge none">no plane yet</div>
        <div class="ps-headnote">This zygote has no sperm. Tick one or more genes in the ranking
        and the pseudosperm plane is fitted to them.</div>`;
      $("#ps-readout").innerHTML = "";
      return;
    }
    head.innerHTML = pl.kind === "sperm"
      ? `<div class="ps-badge sperm">real sperm plane</div>
         <div class="ps-headnote">Through the sperm, the polar-body COM and the cell COM — a
         measurement, taken verbatim from the Sperm Division Plane project.</div>`
      : `<div class="ps-badge pseudo">pseudosperm plane · fitted</div>
         <div class="ps-headnote">No sperm here. This is the angle that maximises the selected
         genes' asymmetry — an <b>inference about those genes</b>, not a located sperm.</div>`;

    const L = [];
    const row = (k, v, cls = "") => `<div class="ps-line ${cls}"><span class="k">${k}</span><span class="v">${v}</span></div>`;
    L.push(row("plane angle", `${pl.angle.toFixed(1)}°`, "is-key"));
    if (pl.kind === "pseudo") {
      const s = pl.sweepRes;
      const mean = s.curve.reduce((x, y) => x + y, 0) / s.curve.length;
      const peak = s.curve[pl.k];
      L.push(row("peak / mean asymmetry", `${(peak / (mean || 1)).toFixed(2)}×`,
                 peak / (mean || 1) < 1.3 ? "is-weak" : ""));
      L.push(row("fitted on", `${state.sel.size} gene${state.sel.size > 1 ? "s" : ""}, ` +
                 `${s.n.toLocaleString()} transcripts`));
    } else if (state.sel.size) {
      L.push(row("sperm angle is a measurement", "not fitted"));
    }
    L.push(row("side-A cytoplasm", `${(pl.volA / 1000).toFixed(1)}k`));
    L.push(row("side-B cytoplasm", `${(pl.volB / 1000).toFixed(1)}k`));
    L.push(row("volume split", `${(100 * pl.volA / (pl.volA + pl.volB)).toFixed(2)} %`,
               pl.exact ? "" : "is-approx"));
    const gs = geneStats(e, pl);
    // the probesets are disjoint, so a selected gene is often simply not measured in this embryo
    if (state.sel.size && gs.length < state.sel.size) {
      L.push(row("selected genes measured here", `${gs.length} of ${state.sel.size}`, "is-weak"));
    }
    if (gs.length) {
      L.push(`<div class="ps-sep">selected genes at this plane</div>`);
      gs.sort((a, b) => a.log10p - b.log10p).forEach((g) =>
        L.push(row(g.g, `${g.l2 >= 0 ? "+" : ""}${g.l2.toFixed(2)} log₂ · p ${fmtP(g.log10p)}`)));
    }
    if (!pl.exact) {
      L.push(`<div class="ps-warn">Per-side volume at an off-grid angle is recovered by clipping
        the cytoplasm mesh, not by counting the original voxels — worth about
        <b>${(meta().volume.held_out_fraction_error_mean * 100).toFixed(2)}%</b> of the volume
        split on average. The sperm planes use the exact voxel volumes.</div>`);
    }
    $("#ps-readout").innerHTML = L.join("");
  }

  // ───────── drawer panels ─────────
  const CFG = { displaylogo: false, responsive: true,
                modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                toImageButtonOptions: { format: "png", scale: 4 } };
  const empty = (el, msg) => { Plotly.purge(el); el.innerHTML = `<div class="ps-empty">${msg}</div>`; };
  const baseLay = (xt, yt) => ({
    margin: { l: 58, r: 12, t: 12, b: 44 }, showlegend: true,
    legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
    xaxis: { title: { text: xt, font: { size: 10 } }, zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
    yaxis: { title: { text: yt, font: { size: 10 } }, zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
    paper_bgcolor: "transparent", plot_bgcolor: "transparent",
  });

  function renderCross() {
    const el = $("#ps-cross"), e = cur();
    if (!e) return;
    const pl = activePlane(e);
    if (!pl) return empty(el, "Tick a gene in the ranking to fit this zygote's plane.");
    el.innerHTML = "";
    const out = e.outline || [];
    // the outline lives in the (u,v) frame; the cut at angle θ has normal cosθ·u + sinθ·v, so in
    // that 2-D frame it is the line through the origin along (−sinθ, cosθ)
    const th = pl.angle * Math.PI / 180;
    const nx = Math.cos(th), ny = Math.sin(th);
    const R = Math.max(...out.map((p) => Math.hypot(p[0], p[1])), 1);
    const traces = [{
      type: "scatter", mode: "lines", x: out.concat([out[0] || [0, 0]]).map((p) => p[0]),
      y: out.concat([out[0] || [0, 0]]).map((p) => p[1]), name: "cytoplasm",
      line: { color: "#94a3b8", width: 2 }, fill: "toself", fillcolor: "rgba(148,163,184,0.10)",
      hoverinfo: "skip",
    }, {
      type: "scatter", mode: "lines", x: [-ny * R * 1.15, ny * R * 1.15],
      y: [nx * R * 1.15, -nx * R * 1.15],
      line: { color: pl.kind === "sperm" ? C_PLANE : C_PSEUDO, width: 2.5 },
      name: pl.kind === "sperm" ? "sperm plane" : "pseudosperm plane", hoverinfo: "skip",
    }];
    if (state.scene) {
      const sc = state.scene, zs = sc.z_scale;
      const C = e.com_um, u = e.u_plot, v = e.v_plot;
      const uu = [u[0] * XY, u[1] * XY, u[2] / zs], vv = [v[0] * XY, v[1] * XY, v[2] / zs];
      for (const name of state.sel) {
        const t = sc.transcripts[name];
        if (!t) continue;
        const X = [[], []], Y = [[], []];
        for (let i = 0; i < t.x.length; i++) {
          if (!t.s1[i]) continue;
          const d = [t.x[i] * XY - C[0], t.y[i] * XY - C[1], t.gz[i] - C[2]];
          const a = d[0] * uu[0] + d[1] * uu[1] + d[2] * uu[2];
          const b = d[0] * vv[0] + d[1] * vv[1] + d[2] * vv[2];
          const s = a * nx + b * ny > 0 ? 0 : 1;
          X[s].push(a); Y[s].push(b);
        }
        [0, 1].forEach((s) => X[s].length && traces.push({
          type: "scatter", mode: "markers", x: X[s], y: Y[s],
          marker: { size: 3, color: s ? C_B : C_A, opacity: 0.45, line: { width: 0 } },
          name: `${name} · ${s ? "B" : "A"}`, hoverinfo: "skip" }));
      }
    }
    const lay = baseLay("µm", "µm");
    lay.xaxis.scaleanchor = "y"; lay.xaxis.scaleratio = 1;
    lay.xaxis.range = [-R * 1.2, R * 1.2]; lay.yaxis.range = [-R * 1.2, R * 1.2];
    Plotly.newPlot(el, traces, lay, CFG);
    $("#ps-cross-sub").textContent =
      `· ${pl.kind === "sperm" ? "real sperm plane" : "fitted"} at ${pl.angle.toFixed(1)}°` +
      ` · projected down the polar axis`;
  }

  function renderSweep() {
    const el = $("#ps-sweep"), e = cur();
    if (!e) return;
    const gs = selGenes(e);
    if (!gs.length) return empty(el, "Tick a gene in the ranking to sweep it.");
    el.innerHTML = "";
    const s = sweep(e, gs);
    const step = meta().grid.step_deg;
    const x = Array.from({ length: K() }, (_, i) => i * step);
    const y = Array.from(s.curve);
    const traces = [{
      type: "scatter", mode: "lines", x, y, name: "asymmetry",
      line: { color: C_PSEUDO, width: 2 },
      hovertemplate: "%{x:.0f}°<br>%{y:.4g}<extra></extra>",
    }, {
      type: "scatter", mode: "markers", x: [s.best * step], y: [y[s.best]], name: "best angle",
      marker: { size: 11, color: C_PSEUDO, line: { color: "#fff", width: 2 } },
      hovertemplate: `peak %{x:.0f}°<extra></extra>`,
    }];
    const lay = baseLay("plane angle about the polar axis (°)",
                        state.metric === "vol" ? "|Δ concentration| (per half-cell)" : "|count fraction difference|");
    lay.xaxis.range = [0, 180]; lay.xaxis.dtick = 30;
    if (e.sperm) {
      lay.shapes = [{ type: "line", x0: e.sperm.angle_deg, x1: e.sperm.angle_deg, yref: "paper",
                      y0: 0, y1: 1, line: { color: C_PLANE, width: 2, dash: "dash" } }];
      lay.annotations = [{ x: e.sperm.angle_deg, y: 1.0, yref: "paper", yanchor: "bottom",
                           text: "real sperm plane", showarrow: false,
                           font: { size: 10, color: C_PLANE } }];
    }
    Plotly.newPlot(el, traces, lay, CFG);
    const mean = y.reduce((a, b) => a + b, 0) / y.length;
    const ratio = y[s.best] / (mean || 1);
    const off = e.sperm
      ? ` · ${angDiff(s.best * step, e.sperm.angle_deg).toFixed(0)}° from the real sperm plane`
      : "";
    $("#ps-sweep-sub").textContent =
      `· peak ${(s.best * step).toFixed(0)}° · ${ratio.toFixed(2)}× the mean` +
      (ratio < 1.3 ? " (weak — the angle is barely determined)" : "") + off;
  }
  /** Planes are 180°-periodic, so the largest possible separation is 90°. */
  const angDiff = (a, b) => { let d = Math.abs(a - b) % 180; return d > 90 ? 180 - d : d; };

  function renderBars() {
    const el = $("#ps-bars"), e = cur();
    if (!e) return;
    const pl = activePlane(e);
    if (!pl || !state.sel.size) return empty(el, "Tick a gene in the ranking.");
    el.innerHTML = "";
    const gs = geneStats(e, pl).sort((a, b) => a.log10p - b.log10p);
    const lab = gs.map((g) => g.g);
    const traces = [
      { type: "bar", x: lab, y: gs.map((g) => g.cA * 1000), name: "side A",
        marker: { color: C_A },
        customdata: gs.map((g) => [g.a, g.log10p]),
        hovertemplate: "%{x} · side A<br>%{customdata[0]} transcripts<br>%{y:.3g} per 1000 vol<extra></extra>" },
      { type: "bar", x: lab, y: gs.map((g) => g.cB * 1000), name: "side B",
        marker: { color: C_B },
        customdata: gs.map((g) => [g.b, g.log10p]),
        hovertemplate: "%{x} · side B<br>%{customdata[0]} transcripts<br>%{y:.3g} per 1000 vol<extra></extra>" },
    ];
    const lay = baseLay("", "concentration (transcripts per 1000 volume units)");
    lay.barmode = "group";
    lay.margin.b = 60;
    lay.annotations = gs.map((g, i) => ({
      x: i, y: Math.max(g.cA, g.cB) * 1000, yanchor: "bottom", showarrow: false,
      text: `p ${fmtP(g.log10p)}`, font: { size: 9, color: g.log10p < -2 ? "#b91c1c" : "#94a3b8" },
    }));
    Plotly.newPlot(el, traces, lay, CFG);
    $("#ps-bars-sub").textContent =
      `· ${pl.kind === "sperm" ? "real sperm plane" : "fitted plane"} at ${pl.angle.toFixed(1)}°` +
      ` · exact binomial against a ${(100 * pl.volA / (pl.volA + pl.volB)).toFixed(1)}% volume split`;
  }

  /** The honest test: refit every SPERM zygote from the selected genes alone, with its sperm
   *  ignored, and see how far the fit lands from the plane the sperm actually defines. */
  function renderValid() {
    const el = $("#ps-valid"), el2 = $("#ps-valid2");
    if (!state.sel.size) {
      empty(el, "Tick a gene in the ranking to test it.");
      empty(el2, "");
      $("#ps-valid-sub").textContent = "";
      return;
    }
    el.innerHTML = ""; el2.innerHTML = "";
    const step = meta().grid.step_deg;
    const rows = [];
    for (const e of state.data.embryos) {
      if (!e.sperm) continue;
      const gs = selGenes(e);
      if (!gs.length) continue;
      const s = sweep(e, gs);
      const fit = s.best * step;
      rows.push({ id: e.id, label: e.label, fit, real: e.sperm.angle_deg,
                  err: angDiff(fit, e.sperm.angle_deg), n: s.n });
    }
    if (!rows.length) {
      empty(el, "None of the sperm zygotes measures the selected genes — the probesets are disjoint.");
      empty(el2, "");
      return;
    }
    rows.sort((a, b) => a.err - b.err);
    Plotly.newPlot(el, [{
      type: "scatter", mode: "markers", x: rows.map((r) => r.real), y: rows.map((r) => r.fit),
      marker: { size: 9, color: rows.map((r) => r.err), colorscale: "Viridis", reversescale: true,
                cmin: 0, cmax: 90, line: { color: "#fff", width: 1 },
                colorbar: { title: { text: "°off", font: { size: 9 } }, thickness: 10, len: 0.8 } },
      text: rows.map((r) => `${r.label} · ${r.err.toFixed(0)}° off`),
      hovertemplate: "%{text}<br>real %{x:.0f}° → fitted %{y:.0f}°<extra></extra>", name: "zygotes",
    }, {
      type: "scatter", mode: "lines", x: [0, 180], y: [0, 180], name: "exact agreement",
      line: { color: "#111827", width: 1.2, dash: "dash" }, hoverinfo: "skip",
    }], Object.assign(baseLay("real sperm-plane angle (°)", "fitted pseudosperm angle (°)"),
                      { xaxis: { title: { text: "real sperm-plane angle (°)", font: { size: 10 } }, range: [0, 180], dtick: 45, gridcolor: "#eef1f5", tickfont: { size: 9 } },
                        yaxis: { title: { text: "fitted angle (°)", font: { size: 10 } }, range: [0, 180], dtick: 45, gridcolor: "#eef1f5", tickfont: { size: 9 } } }), CFG);

    const errs = rows.map((r) => r.err);
    const med = errs.slice().sort((a, b) => a - b)[Math.floor(errs.length / 2)];
    const mean = errs.reduce((a, b) => a + b, 0) / errs.length;
    // under the null the two angles are independent and uniform, so the separation is uniform
    // on [0,90] with mean 45 — that, not zero, is what the fit has to beat
    const bins = Array.from({ length: 9 }, (_, i) => i * 10);
    const hist = new Array(9).fill(0);
    errs.forEach((v) => hist[Math.min(8, Math.floor(v / 10))]++);
    Plotly.newPlot(el2, [
      { type: "bar", x: bins.map((b) => b + 5), y: hist, name: "observed",
        marker: { color: C_PSEUDO }, width: 9,
        hovertemplate: "%{y} zygotes at %{x:.0f}±5°<extra></extra>" },
      { type: "scatter", mode: "lines", x: bins.map((b) => b + 5),
        y: new Array(9).fill(errs.length / 9), name: "chance",
        line: { color: "#111827", width: 1.5, dash: "dash" }, hoverinfo: "skip" },
    ], Object.assign(baseLay("angle between fitted and real plane (°)", "zygotes"),
                     { bargap: 0.08 }), CFG);
    const better = mean < 45;
    $("#ps-valid-sub").textContent =
      `· ${rows.length} sperm zygotes measure these genes · median ${med.toFixed(0)}°, ` +
      `mean ${mean.toFixed(0)}° off · chance is 45° ` +
      (better ? `(better than chance by ${(45 - mean).toFixed(0)}°)` : "(no better than chance)");
  }

  const RENDER = { cross: renderCross, sweep: renderSweep, bars: renderBars, valid: renderValid };
  function renderPanel() {
    const fn = RENDER[state.tab];
    if (!fn) return;
    try { fn(); } catch (err) { console.error("[pseudosperm]", state.tab, err); }
  }
  function refresh() { render3D(); renderReadout(); renderPanel(); }

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
    $("#ps-metric").addEventListener("change", (e) => { state.metric = e.target.value; refresh(); });
    $("#ps-sort").addEventListener("change", (e) => { state.sort = e.target.value; renderRank(); });
    $("#ps-minm").addEventListener("input", (e) => {
      state.minM = +e.target.value; $("#ps-minm-val").textContent = e.target.value; renderRank();
    });
    $("#ps-top5").addEventListener("click", () => { selectTop(5); renderRank(); refresh(); });
    $("#ps-clear").addEventListener("click", () => { state.sel.clear(); syncSel(); renderRank(); refresh(); });
    [["t-plane", "plane"], ["t-axis", "axis"], ["t-sperm", "sperm"], ["t-tx", "tx"],
     ["t-mesh", "mesh"]].forEach(([id, key]) =>
      $("#" + id).addEventListener("change", (e) => { state.show[key] = e.target.checked; render3D(); }));

    $("#ps-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#ps-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#ps-panels").querySelectorAll(".xs-panel").forEach((p) =>
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
        Math.max(280, Math.min(window.innerWidth - 80, d.v - (ev.clientX - d.x))) + "px");
    }, () => ({ v: $("#rdrawer").getBoundingClientRect().width }));
    V.wireWindow($("#controls"), $("#controls-header"),
                 [...$("#controls").querySelectorAll(".rz")], "pseudosperm_controls_box");
    window.addEventListener("resize", () => {
      ["#ps-cross", "#ps-sweep", "#ps-bars", "#ps-valid", "#ps-valid2"].forEach((s) => {
        try { Plotly.Plots.resize($(s)); } catch (_) {}
      });
    });
  }
})();
