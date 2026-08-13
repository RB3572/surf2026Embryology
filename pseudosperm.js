/* Pseudosperm Division Plane.
 *
 * A port of figure 4.21's method; build_pseudosperm.py documents it in full and reproduces that
 * figure's ranking exactly. Two halves:
 *
 *   A. In the 30 zygotes with a sperm, cut with the plane through {sperm, cytoplasm COM,
 *      polar-body COM} and ask, per gene, whether it AGREES or OPPOSES the bulk transcript
 *      gradient across that cut. Rank genes by how reproducibly they do so.
 *
 *   B. In the 20 zygotes with no sperm, scan every plane containing the COM→polar-body axis and
 *      keep the one where those same genes reproduce the same agreements and oppositions.
 *
 * ORIENTATION IS BY TOTAL CYTOPLASMIC COUNT. The sperm lies ON its own plane, so it cannot pick a
 * side; the "+" half is whichever holds more cytoplasmic transcripts overall — the FULLER half.
 * That is intrinsic to the zygote and means the same thing on a sperm plane and a candidate one.
 *
 * The template is chosen by a P CUTOFF, not by picking genes: everything below the cutoff (and
 * clearing |log2 fold| >= 0.5) is in, weighted by -log10(P).
 *
 * Everything is cytoplasm-only — pronuclei and polar body enter no count and no volume.
 *
 * ⚠️ A pseudosperm plane is FITTED. The Validation tab refits the sperm zygotes with the sperm
 * hidden and shows how far off it lands; chance is 45°, not 0°.
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const PX = 0.15;
  const C_PLANE = "#7c3aed", C_AXIS = "#0f172a", C_SPERM = "#16a34a", C_PSEUDO = "#0891b2";
  const C_AGREE = "#0d9488", C_OPPOSE = "#dc2626";
  const C_F = "#b45309", C_E = "#0d9488";

  const state = {
    data: null, byId: {}, rank: [], curId: null, scene: null,
    pCut: 0.05, tab: "cross",
    show: { plane: true, axis: true, sperm: true, tx: false, mesh: true },
  };
  const meta = () => state.data.meta;
  const cur = () => state.byId[state.curId];
  const K = () => meta().grid.n;
  const STEP = () => meta().grid.step_deg;
  const MIN_ABS = () => meta().params.MIN_ABS_LFC;

  // ───────── the method, mirroring build_pseudosperm.py ─────────
  /** Side-A counts per angle for one gene, from the delta encoding. */
  function counts(g) {
    if (g._a) return g._a;
    const out = new Int32Array(g.a.length);
    let s = 0;
    for (let i = 0; i < g.a.length; i++) { s += g.a[i]; out[i] = s; }
    g._a = out;
    return out;
  }

  /** Bulk-corrected log2 fold per gene: concentration on the fuller half over the emptier one,
   *  minus this embryo's MEDIAN per-gene log ratio. The median, not the ratio of totals — one
   *  abundant gene can carry a third of the panel and would otherwise set the correction. */
  function lfcsOf(cnt, vF, vE) {
    const raw = new Map();
    for (const [g, ab] of cnt) {
      const a = ab[0], b = ab[1];
      if (a + b <= 0) continue;
      raw.set(g, Math.log2(((a + 0.5) / vF) / ((b + 0.5) / vE)));
    }
    if (!raw.size) return raw;
    const v = [...raw.values()].sort((x, y) => x - y);
    const n = v.length;
    const bulk = n % 2 ? v[(n - 1) / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
    for (const [g, x] of raw) raw.set(g, x - bulk);
    return raw;
  }

  /** One candidate plane: orient by total count, then the per-gene folds. */
  function planeAt(e, k) {
    let totA = 0, totB = 0;
    const pair = new Map();
    for (const g of e.genes) {
      const a = counts(g)[k], b = g.n - a;
      totA += a; totB += b;
      pair.set(g.g, [a, b]);
    }
    const flipped = totA < totB;
    let vF = e.volA[k], vE = e.volB[k];
    if (flipped) {
      vF = e.volB[k]; vE = e.volA[k];
      for (const [g, ab] of pair) pair.set(g, [ab[1], ab[0]]);
    }
    return { cnt: pair, vF, vE, flipped, lfc: lfcsOf(pair, vF, vE),
             totF: Math.max(totA, totB), totE: Math.min(totA, totB) };
  }

  /** The real sperm plane, from the exact counts the build carried over. */
  function spermPlane(e) {
    const s = e.sperm;
    let totA = 0, totB = 0;
    const pair = new Map();
    for (const g in s.a) {
      const a = s.a[g], b = s.n[g] - a;
      totA += a; totB += b;
      pair.set(g, [a, b]);
    }
    const flipped = totA < totB;
    let vF = s.volA, vE = s.volB;
    if (flipped) {
      vF = s.volB; vE = s.volA;
      for (const [g, ab] of pair) pair.set(g, [ab[1], ab[0]]);
    }
    return { cnt: pair, vF, vE, flipped, lfc: lfcsOf(pair, vF, vE),
             totF: Math.max(totA, totB), totE: Math.min(totA, totB) };
  }

  /** The template: every ranked gene under the P cutoff that also clears |log2 fold|. */
  function template() {
    return state.rank.filter((r) => r.p < state.pCut && Math.abs(r.lfc) >= MIN_ABS())
      .map((r) => ({ g: r.g, sign: r.lfc > 0 ? 1 : -1, weight: r.weight, lfc: r.lfc,
                     side: r.side, p: r.p }));
  }

  /** Weighted mean of sign × log2 fold over the template genes this embryo measures. A gene
   *  reproduces its own sperm-plane direction in proportion to how significant that direction
   *  was, so a gene at P = 0.001 pulls three times as hard as one at P = 0.1. */
  function scoreOf(pl, tmpl) {
    let num = 0, den = 0, used = 0;
    for (const t of tmpl) {
      const v = pl.lfc.get(t.g);
      if (v === undefined) continue;
      num += t.weight * t.sign * v;
      den += t.weight;
      used++;
    }
    return used ? { score: num / den, used } : null;
  }

  /** Scan every meridional plane and keep the best-scoring one. */
  function fit(e, tmpl) {
    const k = K(), curve = new Float64Array(k).fill(NaN);
    let best = -Infinity, bestK = -1, used = 0;
    for (let i = 0; i < k; i++) {
      const s = scoreOf(planeAt(e, i), tmpl);
      if (!s) continue;
      curve[i] = s.score; used = s.used;
      if (s.score > best) { best = s.score; bestK = i; }
    }
    return bestK < 0 ? null : { k: bestK, score: best, curve, used };
  }

  /** The plane this zygote is showing: the real sperm plane if it has one, else the fit. */
  function activePlane(e) {
    const tmpl = template();
    if (e.sperm) {
      const pl = spermPlane(e);
      return { kind: "sperm", pl, angle: e.sperm.angle_deg, k: Math.round(e.sperm.angle_deg / STEP()) % K(),
               tmpl, score: scoreOf(pl, tmpl) };
    }
    if (!tmpl.length) return null;
    const f = fit(e, tmpl);
    if (!f) return null;
    return { kind: "pseudo", pl: planeAt(e, f.k), angle: f.k * STEP(), k: f.k,
             tmpl, fitres: f, score: { score: f.score, used: f.used } };
  }

  const angDiff = (a, b) => { const d = Math.abs(a - b) % 180; return d > 90 ? 180 - d : d; };

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
    state.pCut = m.params.ALPHA;
    $("#embryo-count").textContent =
      `${m.n_embryos} zygotes · ${m.n_sperm} sperm plane · ${m.n_pseudo} pseudosperm · ` +
      `${m.n_genes_ranked} genes ranked`;
    $("#ps-rank-desc").innerHTML =
      `Every gene's asymmetry across the <b>real</b> sperm plane, in the ${m.n_sperm} zygotes that
       have one: a one-sample t-test of its per-embryo log₂ fold against zero. <b>Set a P cutoff</b>
       and everything below it that also clears |log₂ fold| ≥ ${MIN_ABS()} becomes the template.`;

    V.buildTabs($("#tabs"), state.data.embryos, selectEmbryo, (e) => ({
      label: e.label, sub: e.date, cls: e.sperm ? "ps-has-sperm" : "",
      title: `${e.id}${e.sperm ? " · sperm plane @ " + e.sperm.angle_deg.toFixed(1) + "°"
                              : " · no sperm — pseudosperm plane"}`,
    }));
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    syncCut();
    renderRank();
    selectEmbryo(state.data.embryos[0].id);
  })();

  // ───────── the ranking ─────────
  function syncCut() {
    const tmpl = template();
    $("#ps-pcut-val").textContent = state.pCut.toFixed(3);
    $("#ps-tmpl-count").innerHTML =
      `<b>${tmpl.length}</b> gene${tmpl.length === 1 ? "" : "s"} in the template` +
      (tmpl.length ? ` · ${tmpl.filter((t) => t.sign > 0).length} fuller, ` +
                     `${tmpl.filter((t) => t.sign < 0).length} emptier` : "");
  }
  function renderRank() {
    const el = $("#ps-rank"), inT = new Set(template().map((t) => t.g));
    el.innerHTML =
      `<div class="ps-rank-head"><span>#</span><span>gene</span><span>log₂FC</span>
         <span>P</span><span>side</span><span>n</span></div>` +
      state.rank.map((r) => `<div class="ps-rank-row${inT.has(r.g) ? " on" : ""}" data-g="${r.g}"
          title="${r.g} · ${r.n.toLocaleString()} cytoplasmic transcripts in ${r.m} sperm zygotes · fold ${r.fold.toFixed(2)}× vs null ${(r.nullFold || 0).toFixed(2)}×">
        <span class="n">${r.rank}</span>
        <span class="e">${r.g}</span>
        <span class="f" style="color:${r.lfc > 0 ? C_F : C_E}">${r.lfc >= 0 ? "+" : ""}${r.lfc.toFixed(2)}</span>
        <span class="p">${r.p < 0.001 ? r.p.toExponential(0) : r.p.toFixed(3)}</span>
        <span class="s">${r.side === "fuller" ? "F" : "E"}</span>
        <span class="m">${r.m}</span></div>`).join("");
  }

  // ───────── 3-D ─────────
  async function selectEmbryo(id) {
    if (id === state.curId) return;
    state.curId = id;
    V.markActiveTab($("#tabs"), id);
    const e = cur();
    $("#loading").hidden = false; $("#loading-text").textContent = `Loading ${e.label}…`;
    try { state.scene = await V.loadGz(`data/segments/${e.scene}`); }
    catch (err) { state.scene = null; }
    $("#loading").hidden = true;
    $("#controls").hidden = false; $("#placeholder").hidden = true;
    $("#drawer-emb").textContent = e.label;
    refresh();
  }

  const addv = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const mulv = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
  // the segments scenes are one isotropic pixel space: µm = pixel × 0.15 on all three axes
  const toPlot = (p) => [p[0] / PX, p[1] / PX, p[2] / PX];
  /** Half-extent for the drawn plane and axis: the radius of a sphere of the same cytoplasm
   *  volume, with a little overhang. (Cube-rooting the volume directly gives a box edge, which is
   *  ~1.6x too big and pushes the quad well outside the cell.) */
  const halfExtent = (e) => Math.cbrt((3 * e.cyto_vol) / (4 * Math.PI)) * 1.25;

  function planeGeom(e, k) {
    const th = k * STEP() * Math.PI / 180;
    const e1 = e.e1_um, e2 = e.e2_um, a = e.axis_um, C = e.com_um;
    const m = addv(mulv(e1, -Math.sin(th)), mulv(e2, Math.cos(th)));
    const nrm = addv(mulv(e1, Math.cos(th)), mulv(e2, Math.sin(th)));
    const L = halfExtent(e);
    const quad = [[-1, -1], [-1, 1], [1, 1], [1, -1]].map(([sa, sm]) =>
      toPlot(addv(C, addv(mulv(a, sa * L), mulv(m, sm * L)))));
    return { quad, m, nrm, a, C, L };
  }

  function render3D() {
    const e = cur(), host = $("#plot");
    if (!e || !state.scene) return;
    const sc = state.scene;
    const A = activePlane(e);
    const traces = [];
    if (state.show.mesh) traces.push(...V.bodyTraces(sc));

    if (A && state.show.plane) {
      const q = planeGeom(e, A.kind === "sperm" ? A.angle / STEP() : A.k).quad;
      traces.push({ type: "mesh3d", x: q.map((p) => p[0]), y: q.map((p) => p[1]),
        z: q.map((p) => p[2]), i: [0, 0], j: [1, 2], k: [2, 3],
        color: A.kind === "sperm" ? C_PLANE : C_PSEUDO, opacity: 0.3,
        name: A.kind === "sperm" ? "sperm plane" : "pseudosperm plane",
        showlegend: true, hoverinfo: "skip", flatshading: true });
    }
    if (state.show.axis) {
      const L = halfExtent(e);
      const p0 = toPlot(addv(e.com_um, mulv(e.axis_um, -L)));
      const p1 = toPlot(addv(e.com_um, mulv(e.axis_um, L)));
      traces.push({ type: "scatter3d", mode: "lines", x: [p0[0], p1[0]], y: [p0[1], p1[1]],
        z: [p0[2], p1[2]], line: { color: C_AXIS, width: 5, dash: "dot" },
        name: "COM → polar-body axis", hoverinfo: "skip" });
      const pb = toPlot(e.pb_um);
      traces.push({ type: "scatter3d", mode: "markers", x: [pb[0]], y: [pb[1]], z: [pb[2]],
        marker: { size: 6, color: C_AXIS }, name: "polar body",
        hovertemplate: "polar body<extra></extra>" });
    }
    if (state.show.sperm && A) {
      if (e.sperm) {
        const s = toPlot(e.sperm.sperm_um);
        traces.push({ type: "scatter3d", mode: "markers", x: [s[0]], y: [s[1]], z: [s[2]],
          marker: { size: 7, color: C_SPERM, symbol: "diamond", line: { color: "#fff", width: 1 } },
          name: "sperm", hovertemplate: "sperm<extra></extra>" });
      } else {
        // the line a sperm would have to lie along. A plane fixes the line, never which end.
        const G = planeGeom(e, A.k);
        const t1 = toPlot(addv(G.C, mulv(G.m, G.L)));
        const t2 = toPlot(addv(G.C, mulv(G.m, -G.L)));
        const c0 = toPlot(G.C);
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
    if (state.show.tx && A) {
      const G = planeGeom(e, A.kind === "sperm" ? A.angle / STEP() : A.k);
      const nrm = e.sperm && A.kind === "sperm" ? e.sperm.normal_um : G.nrm;
      for (const t of A.tmpl.slice(0, 6)) {
        const rec = sc.transcripts[t.g];
        if (!rec) continue;
        const X = [[], []], Y = [[], []], Z = [[], []];
        for (let i = 0; i < rec.x.length; i++) {
          if (rec.s[i] !== e.body) continue;                 // cytoplasm only
          const p = [rec.x[i] * PX, rec.y[i] * PX, rec.gz[i] * e.z_scale * PX];
          const d = [p[0] - e.com_um[0], p[1] - e.com_um[1], p[2] - e.com_um[2]];
          const s = d[0] * nrm[0] + d[1] * nrm[1] + d[2] * nrm[2] > 0 ? 0 : 1;
          X[s].push(rec.x[i]); Y[s].push(rec.y[i]); Z[s].push(rec.gz[i] * e.z_scale);
        }
        [0, 1].forEach((s) => X[s].length && traces.push({
          type: "scatter3d", mode: "markers", x: X[s], y: Y[s], z: Z[s],
          marker: { size: 1.8, color: s ? C_E : C_F, opacity: 0.6 },
          name: `${t.g} · ${s ? "B" : "A"}`, hoverinfo: "skip" }));
      }
    }
    Plotly.react(host, traces, V.sceneLayout(sc.extents, sc.id), V.plotConfig);
  }

  // ───────── readout ─────────
  function renderReadout() {
    const e = cur(); if (!e) return;
    const A = activePlane(e);
    const head = $("#ps-plane-head");
    if (!A) {
      head.innerHTML = `<div class="ps-badge none">no plane</div>
        <div class="ps-headnote">The template is empty at this P cutoff, or none of its genes is on
        this embryo's probeset. Raise the cutoff.</div>`;
      $("#ps-readout").innerHTML = "";
      return;
    }
    head.innerHTML = A.kind === "sperm"
      ? `<div class="ps-badge sperm">real sperm plane</div>
         <div class="ps-headnote">Through the sperm, the cytoplasm COM and the polar-body COM. A
         measurement — the template is not used to place it.</div>`
      : `<div class="ps-badge pseudo">pseudosperm plane · fitted</div>
         <div class="ps-headnote">The meridional plane where the template genes best reproduce the
         agreements and oppositions they showed across the real sperm planes. An
         <b>inference</b>, not a located sperm.</div>`;

    const L = [];
    const row = (k, v, cls = "") => `<div class="ps-line ${cls}"><span class="k">${k}</span><span class="v">${v}</span></div>`;
    L.push(row("plane angle", `${A.angle.toFixed(1)}°`, "is-key"));
    if (A.score) {
      L.push(row("template score", A.score.score.toFixed(4), "is-key"));
      L.push(row("template genes here", `${A.score.used} of ${A.tmpl.length}`,
                 A.score.used < 3 ? "is-weak" : ""));
    }
    if (A.kind === "pseudo" && A.fitres) {
      const c = [...A.fitres.curve].filter((x) => !isNaN(x));
      const spread = Math.max(...c) - Math.min(...c);
      L.push(row("peak above trough", spread.toFixed(4), spread < 0.05 ? "is-weak" : ""));
    }
    if (e.sperm) {
      L.push(row("sperm sits on its plane", `${Math.abs(e.sperm.dist_to_plane_um).toFixed(3)} µm`));
    }
    L.push(`<div class="ps-sep">cytoplasm, this plane</div>`);
    L.push(row("fuller half", `${A.pl.totF.toLocaleString()} tx · ${(A.pl.vF / 1000).toFixed(1)}k µm³`));
    L.push(row("emptier half", `${A.pl.totE.toLocaleString()} tx · ${(A.pl.vE / 1000).toFixed(1)}k µm³`));
    L.push(row("whole cytoplasm", `${(e.cyto_vol / 1000).toFixed(1)}k µm³`));
    L.push(`<div class="ps-warn">Pronuclei and the polar body are excluded from every count and
      every volume — counts by the per-molecule segment label, volume by the cytoplasm label's own
      voxel volume.</div>`);
    $("#ps-readout").innerHTML = L.join("");
  }

  // ───────── drawer panels ─────────
  const CFG = { displaylogo: false, responsive: true,
                modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
                toImageButtonOptions: { format: "png", scale: 4 } };
  const empty = (el, msg) => { Plotly.purge(el); el.innerHTML = `<div class="ps-empty">${msg}</div>`; };
  const baseLay = (xt, yt) => ({
    margin: { l: 60, r: 12, t: 12, b: 46 }, showlegend: true,
    legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
    xaxis: { title: { text: xt, font: { size: 10 } }, zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
    yaxis: { title: { text: yt, font: { size: 10 } }, zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
    paper_bgcolor: "transparent", plot_bgcolor: "transparent",
  });

  function renderCross() {
    const el = $("#ps-cross"), e = cur(); if (!e) return;
    const A = activePlane(e);
    if (!A) return empty(el, "No plane at this P cutoff.");
    el.innerHTML = "";
    const th = A.angle * Math.PI / 180;
    const nx = Math.cos(th), ny = Math.sin(th);
    const out = e.outline || [];
    const R = Math.max(...out.map((p) => Math.hypot(p[0], p[1])), 1);
    const ring = out.concat([out[0] || [0, 0]]);
    const traces = [
      { type: "scatter", mode: "lines", x: ring.map((p) => p[0]), y: ring.map((p) => p[1]),
        name: "cytoplasm", line: { color: "#94a3b8", width: 2 }, fill: "toself",
        fillcolor: "rgba(148,163,184,0.10)", hoverinfo: "skip" },
      { type: "scatter", mode: "lines", x: [-ny * R * 1.15, ny * R * 1.15],
        y: [nx * R * 1.15, -nx * R * 1.15],
        line: { color: A.kind === "sperm" ? C_PLANE : C_PSEUDO, width: 2.5 },
        name: A.kind === "sperm" ? "sperm plane" : "pseudosperm plane", hoverinfo: "skip" },
    ];
    if (state.scene) {
      const sc = state.scene, e1 = e.e1_um, e2 = e.e2_um, C = e.com_um;
      for (const t of A.tmpl.slice(0, 6)) {
        const rec = sc.transcripts[t.g];
        if (!rec) continue;
        const X = [[], []], Y = [[], []];
        for (let i = 0; i < rec.x.length; i++) {
          if (rec.s[i] !== e.body) continue;
          const d = [rec.x[i] * PX - C[0], rec.y[i] * PX - C[1],
                     rec.gz[i] * e.z_scale * PX - C[2]];
          const a = d[0] * e1[0] + d[1] * e1[1] + d[2] * e1[2];
          const b = d[0] * e2[0] + d[1] * e2[1] + d[2] * e2[2];
          const s = a * nx + b * ny > 0 ? 0 : 1;
          X[s].push(a); Y[s].push(b);
        }
        [0, 1].forEach((s) => X[s].length && traces.push({
          type: "scatter", mode: "markers", x: X[s], y: Y[s],
          marker: { size: 3, color: s ? C_E : C_F, opacity: 0.45, line: { width: 0 } },
          name: `${t.g} · ${s ? "B" : "A"}`, hoverinfo: "skip" }));
      }
    }
    const lay = baseLay("µm", "µm");
    lay.xaxis.scaleanchor = "y"; lay.xaxis.scaleratio = 1;
    lay.xaxis.range = [-R * 1.2, R * 1.2]; lay.yaxis.range = [-R * 1.2, R * 1.2];
    Plotly.newPlot(el, traces, lay, CFG);
    $("#ps-cross-sub").textContent =
      `· ${A.kind === "sperm" ? "real sperm plane" : "fitted"} at ${A.angle.toFixed(1)}°` +
      ` · looking down the polar axis · up to 6 template genes shown`;
  }

  function renderSweep() {
    const el = $("#ps-sweep"), e = cur(); if (!e) return;
    const tmpl = template();
    if (!tmpl.length) return empty(el, "The template is empty at this P cutoff.");
    const f = fit(e, tmpl);
    if (!f) return empty(el, "None of the template genes is on this embryo's probeset.");
    el.innerHTML = "";
    const x = Array.from({ length: K() }, (_, i) => i * STEP());
    const y = Array.from(f.curve);
    const traces = [
      { type: "scatter", mode: "lines", x, y, name: "template score",
        line: { color: C_PSEUDO, width: 2 },
        hovertemplate: "%{x:.0f}°<br>score %{y:.4f}<extra></extra>" },
      { type: "scatter", mode: "markers", x: [f.k * STEP()], y: [y[f.k]], name: "best plane",
        marker: { size: 11, color: C_PSEUDO, line: { color: "#fff", width: 2 } },
        hovertemplate: "peak %{x:.0f}°<extra></extra>" },
    ];
    const lay = baseLay("plane angle about the COM → polar-body axis (°)",
                        "weighted mean of sign × log₂ fold");
    lay.xaxis.range = [0, 180]; lay.xaxis.dtick = 30;
    lay.shapes = [{ type: "line", x0: 0, x1: 180, y0: 0, y1: 0,
                    line: { color: "#94a3b8", width: 1, dash: "dot" } }];
    if (e.sperm) {
      lay.shapes.push({ type: "line", x0: e.sperm.angle_deg, x1: e.sperm.angle_deg, yref: "paper",
                        y0: 0, y1: 1, line: { color: C_PLANE, width: 2, dash: "dash" } });
      lay.annotations = [{ x: e.sperm.angle_deg, y: 1.0, yref: "paper", yanchor: "bottom",
                           text: "real sperm plane", showarrow: false,
                           font: { size: 10, color: C_PLANE } }];
    }
    Plotly.newPlot(el, traces, lay, CFG);
    const c = y.filter((v) => !isNaN(v));
    const spread = Math.max(...c) - Math.min(...c);
    $("#ps-sweep-sub").textContent =
      `· peak ${(f.k * STEP()).toFixed(0)}° · ${f.used} template gene${f.used === 1 ? "" : "s"} here` +
      ` · peak−trough ${spread.toFixed(3)}` + (spread < 0.05 ? " (flat — barely determined)" : "") +
      (e.sperm ? ` · ${angDiff(f.k * STEP(), e.sperm.angle_deg).toFixed(0)}° from the real plane` : "");
  }

  /** Which template genes agree with their sperm-plane direction at this plane, and which oppose. */
  function renderAgree() {
    const el = $("#ps-agree"), e = cur(); if (!e) return;
    const A = activePlane(e);
    if (!A) return empty(el, "No plane at this P cutoff.");
    el.innerHTML = "";
    const rows = [];
    for (const t of A.tmpl) {
      const v = A.pl.lfc.get(t.g);
      if (v === undefined) continue;
      rows.push({ g: t.g, v, sign: t.sign, side: t.side, w: t.weight, agree: t.sign * v > 0,
                  contrib: t.sign * v });
    }
    if (!rows.length) return empty(el, "None of the template genes is on this embryo's probeset.");
    rows.sort((a, b) => b.contrib - a.contrib);
    Plotly.newPlot(el, [{
      type: "bar", orientation: "h",
      y: rows.map((r) => `${r.g} (${r.side === "fuller" ? "F" : "E"})`),
      x: rows.map((r) => r.contrib),
      marker: { color: rows.map((r) => (r.agree ? C_AGREE : C_OPPOSE)) },
      customdata: rows.map((r) => [r.v, r.sign > 0 ? "fuller" : "emptier", r.w]),
      hovertemplate: "%{y}<br>log₂ fold here %{customdata[0]:.2f}<br>" +
                     "wants the %{customdata[1]} half<br>weight %{customdata[2]:.2f}<extra></extra>",
      name: "sign × log₂ fold",
    }], Object.assign(baseLay("sign × log₂ fold  (positive = agrees with its sperm-plane direction)", ""),
      { showlegend: false, margin: { l: 120, r: 14, t: 12, b: 48 },
        shapes: [{ type: "line", x0: 0, x1: 0, yref: "paper", y0: 0, y1: 1,
                   line: { color: "#111827", width: 1.2 } }] }), CFG);
    const nA = rows.filter((r) => r.agree).length;
    $("#ps-agree-sub").textContent =
      `· ${nA} of ${rows.length} template genes agree here · ` +
      `score ${A.score ? A.score.score.toFixed(4) : "—"}` +
      (A.kind === "sperm" ? " · at the REAL sperm plane, not a fit" : "");
  }

  /** The honest test: refit each SPERM zygote from the template alone, sperm ignored. */
  function renderValid() {
    const el = $("#ps-valid"), el2 = $("#ps-valid2");
    const tmpl = template();
    if (!tmpl.length) { empty(el, "The template is empty."); empty(el2, ""); return; }
    el.innerHTML = ""; el2.innerHTML = "";
    const rows = [];
    for (const e of state.data.embryos) {
      if (!e.sperm) continue;
      const f = fit(e, tmpl);
      if (!f) continue;
      rows.push({ label: e.label, fit: f.k * STEP(), real: e.sperm.angle_deg,
                  err: angDiff(f.k * STEP(), e.sperm.angle_deg) });
    }
    if (!rows.length) {
      empty(el, "No sperm zygote measures a template gene — the probesets are disjoint.");
      empty(el2, "");
      $("#ps-valid-sub").textContent = "";
      return;
    }
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
    }], Object.assign(baseLay("real sperm-plane angle (°)", "fitted angle (°)"), {
      xaxis: { title: { text: "real sperm-plane angle (°)", font: { size: 10 } }, range: [0, 180],
               dtick: 45, gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { title: { text: "fitted angle (°)", font: { size: 10 } }, range: [0, 180],
               dtick: 45, gridcolor: "#eef1f5", tickfont: { size: 9 } } }), CFG);

    const errs = rows.map((r) => r.err);
    const med = errs.slice().sort((a, b) => a - b)[Math.floor(errs.length / 2)];
    const mean = errs.reduce((a, b) => a + b, 0) / errs.length;
    const hist = new Array(9).fill(0);
    errs.forEach((v) => hist[Math.min(8, Math.floor(v / 10))]++);
    const bins = Array.from({ length: 9 }, (_, i) => i * 10 + 5);
    Plotly.newPlot(el2, [
      { type: "bar", x: bins, y: hist, name: "observed", marker: { color: C_PSEUDO }, width: 9,
        hovertemplate: "%{y} zygotes at %{x:.0f}±5°<extra></extra>" },
      { type: "scatter", mode: "lines", x: bins, y: new Array(9).fill(errs.length / 9),
        name: "chance", line: { color: "#111827", width: 1.5, dash: "dash" }, hoverinfo: "skip" },
    ], Object.assign(baseLay("angle between fitted and real plane (°)", "zygotes"),
                     { bargap: 0.08 }), CFG);
    $("#ps-valid-sub").textContent =
      `· ${rows.length} sperm zygotes carry a template gene · median ${med.toFixed(0)}°, ` +
      `mean ${mean.toFixed(0)}° off · chance is 45° ` +
      (mean < 45 ? `(better by ${(45 - mean).toFixed(0)}°)` : "(no better than chance)");
  }

  const RENDER = { cross: renderCross, sweep: renderSweep, agree: renderAgree, valid: renderValid };
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
    // the P cutoff is the one control: it defines the template, and everything follows
    const slider = $("#ps-pcut");
    const CUTS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2,
                  meta().p_full_coverage || 0.255, 0.3, 0.4, 0.5, 0.75, 1.0]
      .filter((v, i, a) => v && a.indexOf(v) === i).sort((a, b) => a - b);
    slider.max = String(CUTS.length - 1);
    slider.value = String(Math.max(0, CUTS.indexOf(meta().params.ALPHA)));
    slider.addEventListener("input", (ev) => {
      state.pCut = CUTS[+ev.target.value];
      syncCut(); renderRank(); refresh();
    });
    $("#ps-cover").addEventListener("click", () => {
      const t = meta().p_full_coverage;
      if (!t) return;
      slider.value = String(CUTS.indexOf(t));
      state.pCut = t; syncCut(); renderRank(); refresh();
    });
    [["t-plane", "plane"], ["t-axis", "axis"], ["t-sperm", "sperm"], ["t-tx", "tx"],
     ["t-mesh", "mesh"]].forEach(([id, key]) =>
      $("#" + id).addEventListener("change", (ev) => {
        state.show[key] = ev.target.checked; render3D();
      }));

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
      ["#ps-cross", "#ps-sweep", "#ps-agree", "#ps-valid", "#ps-valid2"].forEach((s) => {
        try { Plotly.Plots.resize($(s)); } catch (_) {}
      });
    });
  }
})();
