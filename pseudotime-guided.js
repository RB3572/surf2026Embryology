/* Pronuclear Pseudotime Calibration — GUIDED VISUAL NARRATIVE.
 *
 * A seven-step explanation of the frozen clock, sitting in front of the existing technical views
 * (which remain untouched under "Advanced analysis"). Everything here is READ FROM THE ARTIFACTS:
 * no metric, coefficient or embryo id is hard-coded, and nothing is retrained or recomputed.
 *
 * Honesty rules baked into this file:
 *   · The τ schematic in step 1 is labelled a SCHEMATIC. Its pronuclear separation is driven by the
 *     cohort-median MEASURED distance sum, so even the cartoon is data-derived — but the dataset
 *     has no 3-D positions or angles, so it is never presented as observed microscopy.
 *   · The 3-D panels draw the ACTUAL measured centroids and cell centre from the geometry cache,
 *     so the drawn segments and the printed d_near / d_far cannot disagree.
 *   · τ colour never touches microscopy: meshes stay grey/blue/red, and τ lives in borders,
 *     timelines, dots and labels only.
 *   · The ±interval is described as an empirical interval measured on held-out embryos of the same
 *     published live cohort — never as a coverage guarantee for our fixed MERFISH zygotes.
 *   · caution / out-of-domain embryos are always shown; filters default to everything on.
 */
window.PTGuided = (() => {
  const $ = (s) => document.querySelector(s);
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const fmt = (v, d = 3) => (v == null || !isFinite(v) ? "—" : Number(v).toFixed(d));
  const pct = (v, d = 1) => (v == null || !isFinite(v) ? "—" : (100 * v).toFixed(d) + "%");

  // perceptually-uniform, colour-blind-safe τ scale (viridis), used ONLY for τ — never for tissue
  const VIRIDIS = [[0, "#440154"], [0.25, "#3b528b"], [0.5, "#21918c"], [0.75, "#5ec962"], [1, "#fde725"]];
  const viridisAt = (t) => {
    t = Math.max(0, Math.min(1, t));
    const hex2 = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
    for (let i = 1; i < VIRIDIS.length; i++) {
      if (t <= VIRIDIS[i][0]) {
        const [a, b] = [VIRIDIS[i - 1], VIRIDIS[i]];
        const f = (t - a[0]) / (b[0] - a[0]);
        const ca = hex2(a[1]), cb = hex2(b[1]);
        return `rgb(${Math.round(ca[0] + (cb[0] - ca[0]) * f)},${Math.round(ca[1] + (cb[1] - ca[1]) * f)},${Math.round(ca[2] + (cb[2] - ca[2]) * f)})`;
      }
    }
    return VIRIDIS[VIRIDIS.length - 1][1];
  };
  const QC_C = { pass: "#15803d", caution: "#b45309", "out-of-domain": "#b42318" };
  const INK = "#1a2233", GRID = "#eef1f5", MUTED = "#64748b";

  const S = { cal: null, frames: null, traj: null, fixed: null, geom: null, V: null,
              sel: null, scene: null, sceneId: null, anim: 0, replay: 0, drawn: {} };

  const LAY = (o) => Object.assign({
    margin: { l: 56, r: 14, t: 10, b: 44 }, paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    font: { color: INK, size: 11.5 }, hovermode: "closest", showlegend: false,
    xaxis: { gridcolor: GRID, zeroline: false }, yaxis: { gridcolor: GRID, zeroline: false },
  }, o || {});
  const CFG = { responsive: true, displaylogo: false, displayModeBar: false };
  const plot = (id, tr, lay) => { const el = $("#" + id); if (el) Plotly.react(el, tr, LAY(lay), CFG); };

  // ── frozen isotonic model, evaluated exactly as build_pronuclei_pseudotime.py does ──
  let SPEC = null;
  function tauOf(sumUm) {
    if (!SPEC) return null;
    const x = SPEC.sign * sumUm, kx = SPEC.knots_x, ky = SPEC.knots_y;
    if (x <= kx[0]) return ky[0];
    if (x >= kx[kx.length - 1]) return ky[ky.length - 1];
    let lo = 0, hi = kx.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (kx[m] <= x) lo = m; else hi = m; }
    const f = (x - kx[lo]) / (kx[hi] - kx[lo] || 1);
    return ky[lo] + f * (ky[hi] - ky[lo]);
  }

  const okFixed = () => S.fixed.embryos.filter((e) => e.tau != null && e.features);
  const selRec = () => S.fixed.embryos.find((e) => e.id === S.sel) || null;

  // ═══════════════════════════════════ init ═══════════════════════════════════
  function init(d) {
    S.cal = d.cal; S.frames = d.frames; S.traj = d.traj; S.fixed = d.fixed; S.V = d.V;
    S.geom = d.geom || null;
    SPEC = S.cal.production.spec;
    if (!SPEC || !SPEC.knots_x) { SPEC = null; }
    buildRail();
    const rows = okFixed().slice().sort((a, b) => a.tau - b.tau);
    $("#gd-embryo").innerHTML = rows.map((r) =>
      `<option value="${esc(r.id)}">${esc(r.label)} · τ ${fmt(r.tau, 2)} · ${esc(r.qc)}</option>`).join("");
    S.sel = rows.length ? rows[Math.floor(rows.length / 2)].id : null;
    $("#gd-embryo").value = S.sel;
    $("#gd-embryo").addEventListener("change", () => selectEmbryo($("#gd-embryo").value));
    $("#gd-scrub").addEventListener("input", drawTraining);
    $("#gd-project").addEventListener("click", projectAnim);
    $("#gd-replay").addEventListener("click", toggleReplay);
    $("#gd-replay-sel").addEventListener("change", () => { stopReplay(); drawReplay(0); });
    ["gd-qc-pass", "gd-qc-caution", "gd-qc-ood"].forEach((i) =>
      $("#" + i).addEventListener("change", drawFixed));
    wireExports();

    drawTimeline(); drawUnknown(); drawTraining(); drawCurve();
    drawFolds(); drawKpis(); drawHeldout(); drawResid();
    fillReplaySel(); drawReplay(0);
    drawFixed(); drawShift(); drawExamples(); drawLimits(); drawCaptions();
    selectEmbryo(S.sel);

    let rs = 0;
    window.addEventListener("resize", () => { clearTimeout(rs); rs = setTimeout(resizeAll, 140); });
  }
  function resizeAll() {
    document.querySelectorAll("#gd-main .gd-fig").forEach((p) => {
      try { if (p.offsetParent && p.data) Plotly.Plots.resize(p); } catch (_) {}
    });
  }

  // .gd-main scrolls internally, so jumps are computed against it rather than the window
  function jumpTo(id) {
    const el = document.getElementById(id), sc = $("#gd-main");
    if (!el || !sc) return;
    sc.scrollTo({ top: sc.scrollTop + el.getBoundingClientRect().top - sc.getBoundingClientRect().top - 8,
                  behavior: "smooth" });
  }

  function buildRail() {
    const steps = [["gd-1", "τ scale"], ["gd-2", "Snapshot"], ["gd-3", "Measurement"],
                   ["gd-4", "Training data"], ["gd-5", "Model"], ["gd-6", "Held-out test"],
                   ["gd-7", "Our zygotes"]];
    $("#gd-rail").innerHTML = steps.map(([id, l], i) =>
      `<button class="gd-rb" data-t="${id}"><span>${i + 1}</span>${esc(l)}</button>`).join("");
    $("#gd-rail").addEventListener("click", (e) => {
      const b = e.target.closest(".gd-rb"); if (!b) return;
      jumpTo(b.dataset.t);
    });
    const obs = new IntersectionObserver((es) => {
      es.forEach((en) => {
        if (!en.isIntersecting) return;
        document.querySelectorAll(".gd-rb").forEach((b) =>
          b.classList.toggle("active", b.dataset.t === en.target.id));
      });
    }, { rootMargin: "-45% 0px -50% 0px" });
    document.querySelectorAll(".gd-step[id]").forEach((s) => obs.observe(s));
  }

  // ═════════════════════ step 1 · τ timeline (SVG schematic) ═════════════════════
  // Separation is driven by the cohort-median MEASURED distance sum at each τ, so the cartoon
  // tracks real data. It is still explicitly a schematic: the source has no positions or angles.
  function medianSumAt(t) {
    const v = [];
    S.traj.forEach((e) => {
      let bi = 0, bd = 9e9;
      for (let i = 0; i < e.tau_true.length; i++) {
        const d = Math.abs(e.tau_true[i] - t); if (d < bd) { bd = d; bi = i; }
      }
      v.push(e.nearer_um[bi] + e.farther_um[bi]);
    });
    v.sort((a, b) => a - b);
    return v.length ? v[v.length >> 1] : 0;
  }
  function drawTimeline() {
    const host = $("#fig-timeline");
    const W = 1120, H = 300;
    const sums = [];
    for (let i = 0; i <= 100; i++) sums.push(medianSumAt(i / 100));
    const sMax = Math.max(...sums), sMin = Math.min(...sums);
    const stops = [];
    for (let i = 0; i <= 10; i++) stops.push(`<stop offset="${i * 10}%" stop-color="${viridisAt(i / 10)}"/>`);
    const x0 = 90, x1 = W - 60, ybar = 232;
    const marks = [[0, "τ = 0", "pronuclear formation"], [1 / 3, "0.33", "early"],
                   [2 / 3, "0.67", "late"], [1, "τ = 1", "NEBD"]];
    host.innerHTML =
      `<svg class="gd-svg" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" role="img"
            aria-label="Normalized zygote timeline from pronuclear formation to nuclear envelope breakdown">
        <defs><linearGradient id="gdTau" x1="0" x2="1">${stops.join("")}</linearGradient></defs>
        <text x="${x0}" y="26" class="s-t">Normalized pronuclear-migration time τ</text>
        <text x="${x0}" y="46" class="s-m">every zygote is rescaled onto the same 0→1 axis, so embryos with different absolute durations are comparable</text>
        <g id="gd-cells"></g>
        <rect x="${x0}" y="${ybar}" width="${x1 - x0}" height="13" rx="6.5" fill="url(#gdTau)"/>
        ${marks.map(([t, a, b]) => {
          const x = x0 + (x1 - x0) * t;
          return `<line x1="${x}" y1="${ybar - 7}" x2="${x}" y2="${ybar + 20}" stroke="#334155" stroke-width="1.2"/>
                  <text x="${x}" y="${ybar + 36}" class="s-l">${a}</text>
                  <text x="${x}" y="${ybar + 50}" class="s-s">${b}</text>`;
        }).join("")}
        <circle id="gd-play-dot" cx="${x0}" cy="${ybar + 6.5}" r="7" fill="#fff" stroke="#0f172a" stroke-width="2"/>
        <text x="${x1}" y="26" style="text-anchor:end" class="s-b">SCHEMATIC</text>
        <text x="${x1}" y="44" style="text-anchor:end" class="s-s">illustration, not microscopy</text>
        <text x="${x1}" y="59" style="text-anchor:end" class="s-s">separation follows the cohort-median measured distance sum</text>
      </svg>`;
    // five schematic cells across the timeline
    const g = host.querySelector("#gd-cells");
    const N = 5, R = 42;
    let html = "";
    for (let i = 0; i < N; i++) {
      const t = i / (N - 1), cx = x0 + (x1 - x0) * t, cy = 140;
      const sep = (medianSumAt(t) - sMin) / ((sMax - sMin) || 1);      // 1 = far apart, 0 = together
      const off = 6 + sep * 22;
      html += `<g>
        <circle cx="${cx}" cy="${cy}" r="${R}" fill="#f1f5f9" stroke="#cbd5e1" stroke-width="1.3"/>
        <circle cx="${cx - off}" cy="${cy}" r="13" fill="#dbeafe" stroke="#2563eb" stroke-width="1.6"/>
        <circle cx="${cx + off}" cy="${cy}" r="13" fill="#fee2e2" stroke="#dc2626" stroke-width="1.6"/>
        <circle cx="${cx}" cy="${cy}" r="2.4" fill="#0f172a"/>
        <circle cx="${cx}" cy="${cy + R + 16}" r="5" fill="${viridisAt(t)}"/>
        <text x="${cx}" y="${cy - R - 10}" class="s-l">τ ${t.toFixed(2)}</text>
      </g>`;
    }
    g.innerHTML = html;
    $("#cap-timeline").innerHTML =
      `τ = 0 at <b>pronuclear formation</b>, τ = 1 at <b>nuclear-envelope breakdown (NEBD)</b>. ` +
      `Across the ${S.cal.dataset.n_embryos} live-imaged training zygotes the migration takes ` +
      `<b>${S.cal.dataset.migration_duration_h.min}–${S.cal.dataset.migration_duration_h.max} h</b> ` +
      `(median ${S.cal.dataset.migration_duration_h.median} h), which is exactly why time is ` +
      `normalized rather than measured in hours. The cells above are a <b>schematic</b>: the two ` +
      `pronuclei are drawn moving together because the measured cohort-median distance sum falls ` +
      `from ${sMax.toFixed(1)} µm to ${sMin.toFixed(1)} µm, but the source data contain no 3-D ` +
      `positions or orientations, so this is an illustration, not microscopy.`;
  }

  // ═════════════════════ step 2 · the snapshot question ═════════════════════
  function drawUnknown() {
    const r = selRec();
    const band = r ? [Math.max(0, r.lo95), Math.min(1, r.hi95)] : [0, 1];
    const tr = [
      { type: "scatter", mode: "lines", x: [0, 1], y: [0, 0],
        line: { color: "#e2e8f0", width: 16 }, hoverinfo: "skip" },
    ];
    if (r) {
      tr.push(
        { type: "scatter", mode: "lines", x: band, y: [0, 0],
          line: { color: "rgba(33,145,140,0.30)", width: 16 }, hoverinfo: "skip" },
        { type: "scatter", mode: "markers", x: [r.tau], y: [0],
          marker: { size: 15, color: viridisAt(r.tau), line: { color: "#fff", width: 2 } },
          hovertemplate: `estimated τ ${fmt(r.tau)}<extra></extra>` });
    }
    plot("fig-unknown", tr, {
      margin: { l: 12, r: 12, t: 6, b: 34 }, height: 96,
      xaxis: { range: [-0.02, 1.02], gridcolor: GRID, zeroline: false,
               tickvals: [0, 0.25, 0.5, 0.75, 1], title: { text: "τ", font: { size: 10 } } },
      yaxis: { visible: false, range: [-1, 1] } });
  }

  // ═════════════════════ step 3 · inputs → output ═════════════════════
  function drawInputs() {
    const r = selRec(); if (!r) return;
    const f = r.features;
    $("#gd-eq").innerHTML =
      `<div class="gd-eq-row"><span class="gd-eq-k">d<sub>near</sub></span>` +
      `<span class="gd-eq-v">${fmt(f.nearer_to_center_um, 2)} µm</span></div>` +
      `<div class="gd-eq-row"><span class="gd-eq-k">d<sub>far</sub></span>` +
      `<span class="gd-eq-v">${fmt(f.farther_to_center_um, 2)} µm</span></div>` +
      `<div class="gd-eq-row gd-eq-sum"><span class="gd-eq-k">S = d<sub>near</sub> + d<sub>far</sub></span>` +
      `<span class="gd-eq-v">${fmt(f.distance_sum_um, 2)} µm</span></div>`;
    const cls = r.qc === "pass" ? "ok" : r.qc === "caution" ? "warn" : "bad";
    $("#gd-flow").innerHTML =
      `<div class="gd-flow-b">geometry<span>S = ${fmt(f.distance_sum_um, 1)} µm</span></div>` +
      `<div class="gd-flow-a">→</div>` +
      `<div class="gd-flow-b">frozen model<span>${esc(S.cal.production.label)}</span></div>` +
      `<div class="gd-flow-a">→</div>` +
      `<div class="gd-flow-b gd-flow-out" style="border-color:${viridisAt(r.tau)}">` +
      `τ ≈ ${fmt(r.tau)}<span>${fmt(r.lo95, 2)} – ${fmt(r.hi95, 2)} · <b class="gd-qc-${cls}">${esc(r.qc)}</b></span></div>`;
  }

  // 3-D: reuse the pronuclei scene + the MEASURED centroids from the geometry cache
  async function loadScene(id) {
    if (S.sceneId === id && S.scene) return S.scene;
    try {
      S.scene = await S.V.loadGz(`data/pronuclei/${id}.json.gz`); S.sceneId = id;
    } catch (_) { S.scene = null; S.sceneId = null; }
    return S.scene;
  }
  function segMesh(sc, lbl, color, opacity, name) {
    const m = sc.region_meshes[String(lbl)]; if (!m) return null;
    const v = m.verts, f = m.faces, nV = v.length / 3, nF = f.length / 3;
    const x = new Array(nV), y = new Array(nV), z = new Array(nV);
    for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
    const a = new Array(nF), b = new Array(nF), c = new Array(nF);
    for (let i = 0; i < nF; i++) { a[i] = f[i * 3]; b[i] = f[i * 3 + 1]; c[i] = f[i * 3 + 2]; }
    return { type: "mesh3d", x, y, z, i: a, j: b, k: c, color, opacity, name,
      showlegend: false, hoverinfo: "name", flatshading: false,
      lighting: { ambient: 0.72, diffuse: 0.5, specular: 0.1, roughness: 0.9 } };
  }
  // geometry cache row -> plot-space points (null when the cache predates the coordinate columns)
  function geomPts(id) {
    const g = S.geom && S.geom[id];
    if (!g || g.center_plot_x === "" || g.center_plot_x == null) return null;
    const P = (t) => [+g[t + "_plot_x"], +g[t + "_plot_y"], +g[t + "_plot_z"]];
    const c = P("center"), n = P("near"), f = P("far");
    if ([...c, ...n, ...f].some((v) => !isFinite(v))) return null;
    return { c, n, f };
  }
  async function draw3D(divId, withMeasure) {
    const r = selRec(); if (!r) return;
    const sc = await loadScene(r.id);
    const el = $("#" + divId); if (!el) return;
    if (!sc) { el.innerHTML = `<div class="gd-empty">3-D scene unavailable for ${esc(r.label)}.</div>`; return; }
    const [la, lb] = sc.pron_labels;
    const traces = [];
    for (const lbl of sc.mask_labels) {
      const pron = lbl === la || lbl === lb;
      // microscopy keeps its own neutral palette — τ colour never tints tissue
      const t = segMesh(sc, lbl, pron ? (lbl === la ? "#2563eb" : "#dc2626") : "#9aa3b2",
                        pron ? 0.42 : 0.07, pron ? `pronucleus ${lbl}` : `segment ${lbl}`);
      if (t) traces.push(t);
    }
    const pts = withMeasure ? geomPts(r.id) : null;
    if (pts) {
      const seg = (p, q, col, nm, dv) => ({ type: "scatter3d", mode: "lines", name: nm,
        x: [p[0], q[0]], y: [p[1], q[1]], z: [p[2], q[2]],
        line: { color: col, width: 5 }, showlegend: false,
        hovertemplate: `${nm} ${fmt(dv, 2)} µm<extra></extra>` });
      traces.push(
        seg(pts.c, pts.n, "#0f172a", "d_near", r.features.nearer_to_center_um),
        seg(pts.c, pts.f, "#64748b", "d_far", r.features.farther_to_center_um),
        { type: "scatter3d", mode: "markers", name: "cell centre",
          x: [pts.c[0]], y: [pts.c[1]], z: [pts.c[2]], showlegend: false,
          marker: { size: 6, color: "#0f172a", symbol: "x" },
          hovertemplate: "cell centre<extra></extra>" },
        { type: "scatter3d", mode: "markers", name: "pronuclear centroids",
          x: [pts.n[0], pts.f[0]], y: [pts.n[1], pts.f[1]], z: [pts.n[2], pts.f[2]],
          showlegend: false, marker: { size: 5, color: ["#2563eb", "#dc2626"],
            line: { color: "#fff", width: 1 } },
          hovertemplate: "pronuclear centroid<extra></extra>" });
    }
    Plotly.react(el, traces, S.V.sceneLayout(sc.extents, "gd" + divId), S.V.plotConfig);
    if (withMeasure) {
      $("#cap-inputs").innerHTML = pts
        ? `Black × marks the cell centre — the centroid of the filled cell (cytoplasm ∪ both ` +
          `pronuclei). The two thin segments are the <b>measured</b> centroid-to-centre distances ` +
          `that produced the numbers on the left; they are read from the geometry cache, not ` +
          `re-derived here, so the drawing and the numbers cannot disagree. Only their <b>sum</b> ` +
          `enters the model — the split into near/far is what makes the feature identity-free, ` +
          `which matters because fixed zygotes carry no reliable male/female pronuclear call.`
        : `<b>Distance overlay unavailable.</b> The geometry cache for this embryo predates the ` +
          `centroid-coordinate columns, so the measured segments cannot be drawn without ` +
          `inventing positions. The distances on the left are still the real measured values. ` +
          `Re-run <code>build_pronuclei_pseudotime.py --extract</code> to populate them.`;
    }
  }

  function selectEmbryo(id) {
    if (!id) return;
    S.sel = id;
    const r = selRec(); if (!r) return;
    if ($("#gd-embryo").value !== id) $("#gd-embryo").value = id;
    $("#gd-q-label").textContent = r.label;
    $("#gd-emb-note").innerHTML =
      `${esc(r.label)} · legacy surface gap ${fmt(r.legacy_surface_gap_um, 2)} µm · ` +
      `QC <b class="gd-qc-${r.qc === "pass" ? "ok" : r.qc === "caution" ? "warn" : "bad"}">${esc(r.qc)}</b>`;
    drawUnknown(); drawInputs(); drawCurve(); drawFixed();
    draw3D("fig-snapshot3d", false); draw3D("fig-inputs3d", true);
    $("#cap-snapshot").innerHTML =
      `This zygote was chemically fixed, so there is <b>no timestamp</b> and no ground-truth τ — ` +
      `not for this embryo, and not for any of the ${S.fixed.meta.n_total} in our cohort. The band ` +
      `is the model's estimate once the geometry is measured; before that, the honest answer is ` +
      `"somewhere on the axis". Everything after this step is about how wide that band has to be.`;
  }

  // ═════════════════════ step 4 · training trajectories ═════════════════════
  function drawTraining() {
    const t = +$("#gd-scrub").value;
    const traces = [];
    S.traj.forEach((e) => {
      traces.push({ type: "scatter", mode: "lines", x: e.tau_true,
        y: e.nearer_um.map((v, i) => v + e.farther_um[i]),
        line: { color: "rgba(100,116,139,0.30)", width: 1 }, hoverinfo: "skip" });
    });
    // cohort median
    const gx = [], gy = [];
    for (let i = 0; i <= 100; i++) { gx.push(i / 100); gy.push(medianSumAt(i / 100)); }
    traces.push({ type: "scatter", mode: "lines", x: gx, y: gy,
      line: { color: "#0f172a", width: 3 }, name: "cohort median",
      hovertemplate: "τ %{x:.2f}<br>median S %{y:.1f} µm<extra></extra>" });
    // scrubber: the measured frame nearest this τ in every trajectory
    const hx = [], hy = [], hc = [];
    S.traj.forEach((e) => {
      let bi = 0, bd = 9e9;
      for (let i = 0; i < e.tau_true.length; i++) {
        const d = Math.abs(e.tau_true[i] - t); if (d < bd) { bd = d; bi = i; }
      }
      hx.push(e.tau_true[bi]); hy.push(e.nearer_um[bi] + e.farther_um[bi]); hc.push(viridisAt(e.tau_true[bi]));
    });
    traces.push({ type: "scatter", mode: "markers", x: hx, y: hy,
      marker: { size: 7, color: hc, line: { color: "#0f172a", width: 0.8 } },
      hovertemplate: "τ %{x:.2f}<br>S %{y:.1f} µm<extra></extra>" });
    plot("fig-training", traces, {
      margin: { l: 62, r: 14, t: 8, b: 46 },
      xaxis: { title: { text: "true τ (known from the live imaging)" }, range: [-0.02, 1.02], gridcolor: GRID, zeroline: false },
      yaxis: { title: { text: "distance sum S = d_near + d_far (µm)" }, gridcolor: GRID, zeroline: false } });
    const sums = hy.slice().sort((a, b) => a - b);
    $("#gd-scrub-note").innerHTML =
      `at τ ≈ <b>${t.toFixed(2)}</b> the ${S.traj.length} embryos span ` +
      `<b>${sums[0].toFixed(1)}–${sums[sums.length - 1].toFixed(1)} µm</b> ` +
      `(median ${sums[sums.length >> 1].toFixed(1)} µm) — that spread is the irreducible noise`;
  }

  // ═════════════════════ step 5 · learned curve ═════════════════════
  function drawCurve(overlay) {
    const fr = S.frames;
    const xs = fr.map((f) => f.nearer_um + f.farther_um);
    const traces = [{ type: "scatter", mode: "markers", x: xs, y: fr.map((f) => f.tau_true),
      marker: { size: 3.4, color: fr.map((f) => viridisAt(f.tau_true)), opacity: 0.30 },
      hovertemplate: "S %{x:.1f} µm<br>true τ %{y:.2f}<extra></extra>" }];
    if (SPEC) {
      const lo = Math.min(...xs), hi = Math.max(...xs), cx = [], cy = [];
      for (let i = 0; i < 260; i++) { const s = lo + (hi - lo) * i / 259; cx.push(s); cy.push(tauOf(s)); }
      traces.push({ type: "scatter", mode: "lines", x: cx, y: cy,
        line: { color: "#0f172a", width: 3, shape: "hv" },
        hovertemplate: "S %{x:.1f} µm → τ %{y:.3f}<extra></extra>" });
    }
    const r = selRec();
    if (r && overlay) {
      const s = r.features.distance_sum_um, tt = r.tau;
      const p = overlay.p;                                   // 0→1 projection progress
      if (p > 0) traces.push({ type: "scatter", mode: "lines", x: [s, s], y: [0, tt * Math.min(1, p * 2)],
        line: { color: viridisAt(tt), width: 2, dash: "dot" }, hoverinfo: "skip" });
      if (p > 0.5) traces.push({ type: "scatter", mode: "lines",
        x: [s, s - (s - Math.min(...xs)) * Math.min(1, (p - 0.5) * 2)], y: [tt, tt],
        line: { color: viridisAt(tt), width: 2, dash: "dot" }, hoverinfo: "skip" });
      traces.push({ type: "scatter", mode: "markers", x: [s], y: [tt],
        marker: { size: 13, color: viridisAt(tt), symbol: "diamond", line: { color: "#fff", width: 2 } },
        hovertemplate: `${esc(r.label)}<br>S ${fmt(s, 1)} µm → τ ${fmt(tt)}<extra></extra>` });
    } else if (r) {
      traces.push({ type: "scatter", mode: "markers", x: [r.features.distance_sum_um], y: [r.tau],
        marker: { size: 12, color: viridisAt(r.tau), symbol: "diamond", line: { color: "#fff", width: 2 } },
        hovertemplate: `${esc(r.label)}<br>S ${fmt(r.features.distance_sum_um, 1)} µm → τ ${fmt(r.tau)}<extra></extra>` });
    }
    plot("fig-curve", traces, {
      margin: { l: 62, r: 14, t: 8, b: 46 },
      xaxis: { title: { text: "distance sum S (µm) — larger = earlier" }, gridcolor: GRID, zeroline: false },
      yaxis: { title: { text: "τ" }, range: [-0.03, 1.03], gridcolor: GRID, zeroline: false } });
  }
  function projectAnim() {
    const r = selRec(); if (!r) return;
    cancelAnimationFrame(S.anim);
    const t0 = performance.now(), D = 1100;
    const step = () => {
      const p = Math.min(1, (performance.now() - t0) / D);
      drawCurve({ p });
      if (p < 1) S.anim = requestAnimationFrame(step);
    };
    S.anim = requestAnimationFrame(step);
    $("#gd-proj-note").innerHTML =
      `${esc(r.label)}: S = <b>${fmt(r.features.distance_sum_um, 1)} µm</b> → up to the curve → ` +
      `across to <b>τ ≈ ${fmt(r.tau)}</b>`;
  }

  // ═════════════════════ step 6 · held-out validation ═════════════════════
  function drawFolds() {
    const n = S.cal.meta.n_outer_folds, ni = S.cal.meta.n_inner_folds;
    const W = 1120, x0 = 176, x1 = W - 30, rowH = 22, gap = 6, top = 52;
    const H = top + n * (rowH + gap) + 54;
    let rows = "";
    for (let k = 0; k < n; k++) {
      const y = top + k * (rowH + gap);
      let cells = "";
      for (let j = 0; j < n; j++) {
        const w = (x1 - x0) / n, x = x0 + j * w, held = j === k;
        cells += `<rect x="${x + 2}" y="${y}" width="${w - 4}" height="${rowH}" rx="3"
          fill="${held ? "#0f172a" : "#e2e8f0"}"/>` +
          `<text x="${x + w / 2}" y="${y + 15}" class="${held ? "s-w" : "s-d"}"
             style="text-anchor:middle">${held ? "held out" : "train + inner CV"}</text>`;
      }
      // CSS text-anchor beats the presentation attribute, so end-anchoring is set inline
      rows += `<text x="${x0 - 12}" y="${y + 15}" class="s-l" style="text-anchor:end">outer fold ${k + 1}</text>${cells}`;
    }
    $("#fig-folds").innerHTML =
      `<svg class="gd-svg" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" role="img"
            aria-label="Nested embryo-grouped cross-validation schematic">
        <text x="20" y="24" class="s-t">Nested, embryo-grouped cross-validation</text>
        <text x="20" y="42" class="s-m" style="text-anchor:start">Every frame of an embryo stays on one side of the split — whole embryos are held out, never individual frames.</text>
        ${rows}
        <text x="20" y="${H - 26}" class="s-m" style="text-anchor:start">The ${ni} inner folds run inside the training embryos only and choose the model family.</text>
        <text x="20" y="${H - 10}" class="s-m" style="text-anchor:start">The chosen model then predicts the untouched outer-fold embryos exactly once — those predictions are what every number below is computed from.</text>
      </svg>`;
  }
  function drawKpis() {
    const m = S.cal.nested_evaluation.outer_test_metrics;
    const ci = S.cal.nested_evaluation.bootstrap_ci95_by_embryo || {};
    const u = S.cal.uncertainty;
    const k = (n, t, s) => `<div class="gd-kpi"><div class="gd-kpi-n">${n}</div>` +
      `<div class="gd-kpi-t">${t}</div><div class="gd-kpi-s">${s}</div></div>`;
    $("#gd-kpis").innerHTML =
      k(fmt(m.macro_mae), "average error in τ",
        `typical miss ≈ ${(m.macro_mae * 100).toFixed(0)}% of the whole window` +
        (ci.macro_mae ? ` · 95% CI ${fmt(ci.macro_mae[0])}–${fmt(ci.macro_mae[1])}` : "")) +
      k(pct(m.pooled_ordering.strict_accuracy), "ordering accuracy",
        "of embryo pairs put in the correct order" +
        (ci.pooled_strict_ordering ? ` · CI ${pct(ci.pooled_strict_ordering[0])}–${pct(ci.pooled_strict_ordering[1])}` : "")) +
      k(fmt(m.pooled_spearman), "rank correlation ρ",
        "monotone agreement with true τ" +
        (ci.pooled_spearman ? ` · CI ${fmt(ci.pooled_spearman[0])}–${fmt(ci.pooled_spearman[1])}` : "")) +
      k("±" + fmt(u.halfwidth_mean), "empirical interval",
        `measured on held-out embryos of the same live cohort — <b>not</b> a coverage guarantee for our fixed zygotes`);
  }
  function drawHeldout() {
    const fr = S.frames, hw = S.cal.uncertainty.halfwidth_mean;
    plot("fig-heldout", [
      { type: "scatter", mode: "lines", x: [0, 1], y: [hw, 1 + hw], line: { width: 0 }, hoverinfo: "skip" },
      { type: "scatter", mode: "lines", x: [0, 1], y: [-hw, 1 - hw], line: { width: 0 },
        fill: "tonexty", fillcolor: "rgba(33,145,140,0.12)", hoverinfo: "skip" },
      { type: "scatter", mode: "lines", x: [0, 1], y: [0, 1],
        line: { color: MUTED, width: 1.4, dash: "dash" }, hoverinfo: "skip" },
      { type: "scatter", mode: "markers", x: fr.map((f) => f.tau_true), y: fr.map((f) => f.tau_pred),
        marker: { size: 3.4, color: fr.map((f) => viridisAt(f.tau_true)), opacity: 0.5 },
        text: fr.map((f) => f.embryo_id),
        hovertemplate: "%{text}<br>true τ %{x:.3f}<br>predicted τ %{y:.3f}<extra></extra>" },
    ], { margin: { l: 56, r: 14, t: 24, b: 44 },
         title: { text: "predicted vs true τ (held-out embryos)", font: { size: 11.5 }, x: 0, xanchor: "left" },
         xaxis: { title: { text: "true τ" }, range: [-0.03, 1.03], gridcolor: GRID, zeroline: false },
         yaxis: { title: { text: "predicted τ" }, range: [-0.03, 1.03], gridcolor: GRID, zeroline: false } });
  }
  function drawResid() {
    const fr = S.frames, hw = S.cal.uncertainty.halfwidth_mean;
    plot("fig-resid", [
      { type: "histogram", x: fr.map((f) => f.residual), nbinsx: 44,
        marker: { color: "#21918c", opacity: 0.8, line: { color: "#fff", width: 0.5 } },
        hovertemplate: "residual %{x}<br>%{y} frames<extra></extra>" },
    ], { margin: { l: 56, r: 14, t: 24, b: 44 },
         title: { text: "error distribution (predicted − true τ)", font: { size: 11.5 }, x: 0, xanchor: "left" },
         xaxis: { title: { text: "error in τ" }, gridcolor: GRID, zeroline: false },
         yaxis: { title: { text: "frames" }, gridcolor: GRID, zeroline: false },
         shapes: [-hw, 0, hw].map((v) => ({ type: "line", x0: v, x1: v, yref: "paper", y0: 0, y1: 1,
           line: { color: v === 0 ? MUTED : "#21918c", width: 1.3, dash: v === 0 ? "dash" : "dot" } })) });
  }
  function fillReplaySel() {
    const t = S.traj.slice().sort((a, b) => a.embryo_id.localeCompare(b.embryo_id));
    $("#gd-replay-sel").innerHTML = t.map((e) =>
      `<option value="${esc(e.embryo_id)}">${esc(e.embryo_id)} · outer fold ${e.outer_fold + 1}</option>`).join("");
  }
  function curTraj() {
    const id = $("#gd-replay-sel").value;
    return S.traj.find((e) => e.embryo_id === id) || S.traj[0];
  }
  function drawReplay(k) {
    const e = curTraj(); if (!e) return;
    const hw = S.cal.uncertainty.halfwidth_mean;
    const n = e.tau_true.length, i = Math.max(0, Math.min(n - 1, k | 0));
    plot("fig-replay", [
      { type: "scatter", mode: "lines", x: e.time_h, y: e.tau_pred.map((v) => Math.min(1, v + hw)),
        line: { width: 0 }, hoverinfo: "skip" },
      { type: "scatter", mode: "lines", x: e.time_h, y: e.tau_pred.map((v) => Math.max(0, v - hw)),
        line: { width: 0 }, fill: "tonexty", fillcolor: "rgba(33,145,140,0.12)", hoverinfo: "skip" },
      { type: "scatter", mode: "lines", x: e.time_h, y: e.tau_true,
        line: { color: "#334155", width: 2 }, name: "true τ",
        hovertemplate: "t %{x:.2f} h<br>true τ %{y:.3f}<extra></extra>" },
      { type: "scatter", mode: "lines", x: e.time_h, y: e.tau_pred,
        line: { color: "#21918c", width: 2 }, name: "predicted τ",
        hovertemplate: "t %{x:.2f} h<br>predicted τ %{y:.3f}<extra></extra>" },
      { type: "scatter", mode: "markers", x: [e.time_h[i], e.time_h[i]], y: [e.tau_true[i], e.tau_pred[i]],
        marker: { size: [10, 10], color: ["#334155", "#21918c"], line: { color: "#fff", width: 1.6 } },
        hovertemplate: "%{y:.3f}<extra></extra>" },
    ], { margin: { l: 56, r: 14, t: 24, b: 44 }, showlegend: true,
         legend: { orientation: "h", x: 0, y: 1.04, yanchor: "bottom", font: { size: 10 } },
         title: { text: "", font: { size: 11 } },
         xaxis: { title: { text: "hours since pronuclear formation" }, gridcolor: GRID, zeroline: false },
         yaxis: { title: { text: "τ" }, range: [-0.05, 1.05], gridcolor: GRID, zeroline: false } });
    const err = Math.abs(e.tau_pred[i] - e.tau_true[i]);
    $("#gd-replay-note").innerHTML =
      `frame ${i + 1}/${n} · true τ <b>${fmt(e.tau_true[i])}</b> · predicted <b>${fmt(e.tau_pred[i])}</b> · ` +
      `error <b>${fmt(err)}</b> — this embryo was <b>never seen</b> by the model that predicted it`;
  }
  function toggleReplay() {
    if (S.replay) { stopReplay(); return; }
    const e = curTraj(); let k = 0;
    $("#gd-replay").textContent = "❚❚ Pause";
    S.replay = setInterval(() => {
      drawReplay(k); k = (k + 1) % e.tau_true.length;
    }, 90);
  }
  function stopReplay() {
    if (S.replay) clearInterval(S.replay);
    S.replay = 0; const b = $("#gd-replay"); if (b) b.textContent = "▶ Play";
  }

  // ═════════════════════ step 7 · our fixed cohort ═════════════════════
  function fixedRows() {
    const on = { pass: $("#gd-qc-pass").checked, caution: $("#gd-qc-caution").checked,
                 "out-of-domain": $("#gd-qc-ood").checked };
    return okFixed().filter((r) => on[r.qc]).sort((a, b) => a.tau - b.tau);
  }
  function drawFixed() {
    const rows = fixedRows();
    const el = $("#fig-fixed");
    const tr = [{ type: "scatter", mode: "markers", x: rows.map((r) => r.tau),
      y: rows.map((r) => r.label),
      error_x: { type: "data", symmetric: false, array: rows.map((r) => r.hi95 - r.tau),
                 arrayminus: rows.map((r) => r.tau - r.lo95), color: "rgba(15,23,42,.28)",
                 thickness: 1, width: 0 },
      marker: { size: 9, color: rows.map((r) => viridisAt(r.tau)),
                line: { color: rows.map((r) => QC_C[r.qc]), width: 2 } },
      customdata: rows.map((r) => r.id),
      text: rows.map((r) => `${r.label}<br>τ ${fmt(r.tau)} [${fmt(r.lo95, 2)}, ${fmt(r.hi95, 2)}]` +
        `<br>S ${fmt(r.features.distance_sum_um, 1)} µm · QC ${r.qc}` +
        (r.reason ? `<br>${r.reason.replace(/; /g, "<br>")}` : "")),
      hovertemplate: "%{text}<extra></extra>" }];
    const H = Math.max(430, 15 * rows.length + 70);
    el.style.height = H + "px";
    plot("fig-fixed", tr, { margin: { l: 104, r: 24, t: 8, b: 46 }, height: H,
      xaxis: { title: { text: "estimated τ · bars = empirical interval (not a guarantee)" },
               range: [-0.04, 1.04], gridcolor: GRID, zeroline: false },
      yaxis: { automargin: true, tickfont: { size: 9 }, type: "category" } });
    if (!el._gdBound) {
      el._gdBound = true;
      el.on("plotly_click", (ev) => {
        const id = ev && ev.points && ev.points[0] && ev.points[0].customdata;
        if (id) { selectEmbryo(id); jumpTo("gd-2"); }
      });
    }
    const qc = okFixed().reduce((m, r) => (m[r.qc] = (m[r.qc] || 0) + 1, m), {});
    $("#cap-fixed").innerHTML =
      `All ${S.fixed.meta.n_predicted} of ${S.fixed.meta.n_total} retained zygotes with two ` +
      `detected pronuclei — <b>${qc.pass || 0} pass, ${qc.caution || 0} caution, ` +
      `${qc["out-of-domain"] || 0} out-of-domain</b>. None are hidden. These are <b>approximate ` +
      `developmental estimates, not timestamps</b>: no fixed embryo has a known true τ, so nothing ` +
      `here can be checked against ground truth the way step 6 could.`;
  }
  function drawShift() {
    const ds = S.fixed.meta.domain_shift || {};
    const live = S.frames.map((f) => f.nearer_um + f.farther_um);
    const fx = okFixed().map((r) => r.features.distance_sum_um);
    plot("fig-shift", [
      { type: "histogram", x: live, histnorm: "probability density", nbinsx: 40,
        name: "live training frames", marker: { color: "rgba(100,116,139,0.55)" },
        hovertemplate: "live S %{x:.0f} µm<extra></extra>" },
      { type: "histogram", x: fx, histnorm: "probability density", nbinsx: 26,
        name: "our fixed zygotes", marker: { color: "rgba(180,35,24,0.55)" },
        hovertemplate: "fixed S %{x:.0f} µm<extra></extra>" },
    ], { barmode: "overlay", margin: { l: 56, r: 14, t: 26, b: 46 }, showlegend: true,
         legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", font: { size: 10 } },
         xaxis: { title: { text: "distance sum S (µm)" }, gridcolor: GRID, zeroline: false },
         yaxis: { title: { text: "density" }, gridcolor: GRID, zeroline: false } });
    const c = ds.consequence || {};
    const p = (ds.distance_sum_um || {}).fixed_median_percentile_in_training;
    $("#gd-shift-warn").innerHTML =
      `<b>Known limitation — the two cohorts do not overlap well.</b> Our fixed zygotes sit at ` +
      `larger distance sums than the live training frames` +
      (p != null ? ` (fixed median ≈ the <b>${p}th percentile</b> of the training distribution)` : "") +
      `, so the model is being asked to work near the edge of what it saw. ` +
      (c.n_at_or_below_tau_0_06 != null
        ? `<b>${c.n_at_or_below_tau_0_06} of ${okFixed().length}</b> fixed zygotes land at or below ` +
          `τ = 0.06, where the fitted curve flattens — their ordering relative to one another is ` +
          `weakly determined. ` : "") +
      `Two explanations cannot be told apart with the data available: these embryos really were ` +
      `fixed early, or there is a systematic offset between live-imaging geometry and our fixed ` +
      `segmentation. Separating them needs raw stacks put through the identical pipeline.`;
  }
  function drawExamples() {
    // representative early / middle / late estimates, chosen dynamically from the QC-pass cohort
    const pass = okFixed().filter((r) => r.qc === "pass").sort((a, b) => a.tau - b.tau);
    if (!pass.length) { $("#gd-examples").innerHTML = ""; return; }
    const pick = (q) => pass[Math.min(pass.length - 1, Math.round(q * (pass.length - 1)))];
    const set = [["earliest estimate", pick(0)], ["middle of the range", pick(0.5)],
                 ["latest estimate", pick(1)]];
    $("#gd-examples").innerHTML = set.map(([k, r]) =>
      `<button class="gd-ex" data-id="${esc(r.id)}" style="border-left-color:${viridisAt(r.tau)}">
        <span class="gd-ex-k">${esc(k)}</span>
        <b>${esc(r.label)}</b>
        <span class="gd-ex-v">τ ${fmt(r.tau)} <span>[${fmt(r.lo95, 2)}, ${fmt(r.hi95, 2)}]</span></span>
        <span class="gd-ex-s">S = ${fmt(r.features.distance_sum_um, 1)} µm</span>
      </button>`).join("");
    $("#gd-examples").onclick = (e) => {
      const b = e.target.closest(".gd-ex"); if (!b) return;
      selectEmbryo(b.dataset.id); jumpTo("gd-2");
    };
  }

  // ═════════════════════ limitations ═════════════════════
  function drawLimits() {
    const u = S.cal.uncertainty, m = S.cal.nested_evaluation.outer_test_metrics;
    const items = [
      ["It is a calibration, not a mechanism",
       "The model is a fitted monotone map from measured distances to normalized time. It says " +
       "nothing about the forces that move the pronuclei and does not simulate their motion."],
      ["The interval is empirical, not guaranteed",
       `±${fmt(u.halfwidth_mean)} τ was measured on held-out embryos from the same published live ` +
       "cohort. It is not a formal coverage guarantee, and applying it to our fixed zygotes " +
       "additionally assumes the live→fixed transfer holds — which is unverified."],
      ["Ordering is stronger than timing",
       `${pct(m.pooled_ordering.strict_accuracy)} of pairs are ordered correctly while the typical ` +
       `error is ${fmt(m.macro_mae)} τ. Treat the output as an approximate ordering with a wide ` +
       "band, not as a clock reading."],
      ["Our cohort sits at the edge of the training range",
       "See step 7: the fixed distance distribution is shifted relative to the live one, and a " +
       "block of embryos is compressed at the early end."],
      ["It does not validate the legacy surface-gap score",
       "The Scheffler source contains cell-centred distances only and no minimum surface-to-surface " +
       "gap, so nothing here calibrates that older metric. It remains available and unchanged."],
      ["No transcript data was involved",
       "Gene and transcript measurements were used nowhere in training, feature choice, model " +
       "selection or uncertainty calibration — so any downstream transcript analysis is not circular."],
    ];
    $("#gd-lims").innerHTML = items.map(([t, b]) =>
      `<div class="gd-lim"><b>${esc(t)}</b><span>${esc(b)}</span></div>`).join("");
  }

  function drawCaptions() {
    const d = S.cal.dataset, meta = S.cal.meta;
    $("#gd-4-sub").innerHTML =
      `${d.n_embryos} untreated live-imaged zygotes · ${d.n_frames.toLocaleString()} frames · ` +
      `Scheffler et al. 2021, sampled every ${d.frame_interval_h} h.`;
    $("#cap-training").innerHTML =
      `Each thin line is one embryo's measured distance sum against its <b>known</b> τ; the dark ` +
      `line is the cohort median. The relationship is clear but far from tight — at any τ the ` +
      `embryos differ by tens of microns, and that biological spread, not measurement error, is ` +
      `what limits the clock. This is measured trajectory data: only one rendered supplementary ` +
      `movie exists for this cohort, so there is no 53-video grid to show.`;
    $("#cap-curve").innerHTML =
      `The frozen model is a <b>${esc(S.cal.production.label)}</b> — a monotone step function fitted ` +
      `to the training embryos only, then locked. It was chosen by nested cross-validation from ` +
      `${S.cal.models.filter((x) => x.deployable).length} candidate families; the comparison is in ` +
      `Advanced analysis. Monotone means it only ever assumes the <i>direction</i> (larger distance ` +
      `sum = earlier), never a constant speed.`;
    $("#cap-heldout").innerHTML =
      `Because whole embryos are held out, these predictions are for zygotes the fitted model had ` +
      `never seen — the honest estimate of how it behaves on a new embryo. These held-out folds are ` +
      `<b>part of the same published cohort</b>, not an independent external dataset.`;
  }

  // ═════════════════════ figure export (SVG preferred, PNG fallback) ═════════════════════
  const FIGS = {
    "01_tau_timeline": { svg: "#fig-timeline" },
    "02_snapshot_question": { plotly: "fig-unknown" },
    "03_model_inputs_output": { plotly: "fig-inputs3d" },
    "04_training_trajectories": { plotly: "fig-training" },
    "05_isotonic_model_curve": { plotly: "fig-curve" },
    "06_heldout_validation": { plotly: "fig-heldout" },
    "07_fixed_cohort_application": { plotly: "fig-fixed" },
  };
  const SLIDE = { w: 1600, h: 900 };                 // 16:9, consistent across every figure
  function wireExports() {
    document.querySelectorAll(".gd-dl").forEach((d) => {
      const key = d.dataset.fig;
      d.innerHTML = `<button class="gd-dlb" data-k="${key}" data-f="svg" title="Download this figure as SVG for slides">SVG</button>` +
                    `<button class="gd-dlb" data-k="${key}" data-f="png" title="Download this figure as a high-resolution PNG">PNG</button>`;
      d.addEventListener("click", (e) => {
        const b = e.target.closest(".gd-dlb"); if (!b) return;
        exportFig(b.dataset.k, b.dataset.f, b);
      });
    });
  }
  function saveBlob(blob, name) {
    const url = URL.createObjectURL(blob), a = document.createElement("a");
    a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1500);
  }
  function exportFig(key, fmtKind, btn) {
    const spec = FIGS[key]; if (!spec) return;
    btn.classList.add("busy");
    const done = () => setTimeout(() => btn.classList.remove("busy"), 400);
    // inline SVG figures serialize directly — no navigation or controls are ever included
    if (spec.svg) {
      const src = document.querySelector(spec.svg + " svg");
      if (!src) { done(); return; }
      const c = src.cloneNode(true);
      c.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      c.setAttribute("width", SLIDE.w); c.setAttribute("height", SLIDE.h);
      const st = document.createElementNS("http://www.w3.org/2000/svg", "style");
      st.textContent = SVG_CSS; c.insertBefore(st, c.firstChild);
      const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      bg.setAttribute("width", "100%"); bg.setAttribute("height", "100%"); bg.setAttribute("fill", "#ffffff");
      c.insertBefore(bg, c.firstChild.nextSibling);
      const txt = new XMLSerializer().serializeToString(c);
      if (fmtKind === "svg") { saveBlob(new Blob([txt], { type: "image/svg+xml" }), key + ".svg"); done(); return; }
      const img = new Image();
      img.onload = () => {
        const cv = document.createElement("canvas");
        cv.width = SLIDE.w * 2; cv.height = SLIDE.h * 2;
        const g = cv.getContext("2d");
        g.fillStyle = "#fff"; g.fillRect(0, 0, cv.width, cv.height);
        g.drawImage(img, 0, 0, cv.width, cv.height);
        cv.toBlob((b) => { saveBlob(b, key + ".png"); done(); });
      };
      img.onerror = done;
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(txt);
      return;
    }
    const gd = document.getElementById(spec.plotly);
    if (!gd || !gd.data) { done(); return; }
    // export on white with the page chrome excluded — Plotly renders only the figure itself
    const layout = Object.assign({}, gd.layout,
      { paper_bgcolor: "#ffffff", plot_bgcolor: "#ffffff", width: SLIDE.w, height: SLIDE.h,
        font: Object.assign({}, gd.layout.font, { size: 18 }) });
    Plotly.toImage({ data: gd.data, layout, config: { displayModeBar: false } },
      { format: fmtKind === "svg" ? "svg" : "png", width: SLIDE.w, height: SLIDE.h,
        scale: fmtKind === "svg" ? 1 : 2 })
      .then((url) => {
        if (fmtKind === "svg") {
          const txt = decodeURIComponent(url.replace(/^data:image\/svg\+xml,/, ""));
          saveBlob(new Blob([txt], { type: "image/svg+xml" }), key + ".svg");
        } else {
          const a = document.createElement("a"); a.href = url; a.download = key + ".png";
          document.body.appendChild(a); a.click(); a.remove();
        }
      }).catch(() => {}).then(done);
  }
  const SVG_CSS =
    ".s-t{font:600 22px system-ui,sans-serif;fill:#0f172a}" +
    ".s-m{font:400 14px system-ui,sans-serif;fill:#64748b}" +
    ".s-l{font:600 13px system-ui,sans-serif;fill:#334155;text-anchor:middle}" +
    ".s-s{font:400 12px system-ui,sans-serif;fill:#64748b;text-anchor:middle}" +
    ".s-b{font:700 13px system-ui,sans-serif;fill:#b45309;letter-spacing:.08em}" +
    ".s-w{font:600 11px system-ui,sans-serif;fill:#fff}" +
    ".s-d{font:600 11px system-ui,sans-serif;fill:#64748b}";

  return { init, resizeAll, stopReplay, selectEmbryo };
})();
