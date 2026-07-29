/* Pronuclei-to-COM Clock (structure-clock τ) — viewer.
 *
 * Organised like the other zygote viewers: a scrolling embryo bar at the top, the
 * selected embryo's 3-D structure in the main window, and a bottom drawer whose
 * FIRST panel spells out exactly what the pseudotime model is fed — the DAPI
 * planes it reads, the structures it finds, the dimensionless numbers it derives,
 * and the tau posterior it returns — before any cross-embryo analysis.
 *
 * Every embryo the model scored appears in the bar, including single-pronucleus
 * and out-of-domain ones; they are marked, never hidden.
 */
(() => {
  const V = window.VCore;
  const $ = (s, r = document) => r.querySelector(s);
  const el = (t, c, h) => { const e = document.createElement(t); if (c) e.className = c; if (h != null) e.innerHTML = h; return e; };
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>]/g, (m) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[m]));
  const f = (v, d = 3) => (v == null || Number.isNaN(v)) ? "—" : (+v).toFixed(d);

  const COL = { cell: "#96a0af", pn1: "#3c82f6", pn2: "#dc3c3c", polar: "#f5af3c" };
  const GRID = "#eef1f5", MUTED = "#64748b";
  const VIRIDIS = [[0, "#440154"], [0.25, "#3b528b"], [0.5, "#21918c"], [0.75, "#5ec962"], [1, "#fde725"]];
  const viridisAt = (t) => {
    t = Math.max(0, Math.min(1, t));
    const hx = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
    for (let i = 1; i < VIRIDIS.length; i++) if (t <= VIRIDIS[i][0]) {
      const a = VIRIDIS[i - 1], b = VIRIDIS[i], q = (t - a[0]) / (b[0] - a[0]);
      const ca = hx(a[1]), cb = hx(b[1]);
      return `rgb(${Math.round(ca[0] + (cb[0] - ca[0]) * q)},${Math.round(ca[1] + (cb[1] - ca[1]) * q)},${Math.round(ca[2] + (cb[2] - ca[2]) * q)})`;
    }
    return VIRIDIS[VIRIDIS.length - 1][1];
  };

  const state = {
    tx: null, inf: null, prev: null, model: null,
    id: null, gene: "__total__", norm: "count", fit: "linear",
    scope: "detected", showOOD: true, tab: "input", drawerOpen: false,
  };
  const statKey = (base) => base + (state.norm === "frac" ? "_frac" : "_count") + (state.scope === "all" ? "_all" : "");

  // ───────────────────────── load ─────────────────────────
  Promise.all([
    fetch("data/pn3d_transcripts.json").then((r) => r.json()),
    fetch("data/pn3d/inference.json").then((r) => r.json()),
    fetch("data/pn3d/preview_index.json").then((r) => r.ok ? r.json() : { embryos: [] }).catch(() => ({ embryos: [] })),
    fetch("data/pn3d/model.json").then((r) => r.json()),
  ]).then(([tx, inf, prev, model]) => {
    state.tx = tx; state.inf = inf; state.prev = prev; state.model = model;
    $("#embryo-count").textContent =
      `${tx.join.n_joined} with transcripts · ${scoredZygotes().length} scored by the clock`;
    buildTabs();
    fillGenes();
    wireControls();
    wireDrawer();
    wirePanels();
    renderJoin();
    renderGeneTable();
    // open on the embryo the model is most confident about, so the input panel is meaningful
    const first = joined().slice().sort((a, b) => a.tau - b.tau)[Math.floor(joined().length / 2)];
    if (first) selectEmbryo(pn3dIdOf(first));
  }).catch((e) => {
    $("#placeholder").innerHTML =
      `<div class="ph-inner"><div class="ph-title">Could not load data</div>
       <div class="ph-sub">Run <code>python3 scripts/build_pn3d_transcripts.py</code>. ${esc(e.message)}</div></div>`;
  });

  const joined = () => state.tx.embryos;
  const scoredZygotes = () => state.inf.embryos.filter((e) => e.stage === "zygote" && e.inferable);
  const pn3dIdOf = (row) => row.pn3d_id;
  const txRowOf = (pnId) => joined().find((r) => r.pn3d_id === pnId);
  const infOf = (pnId) => state.inf.embryos.find((e) => e.embryo_id === pnId);
  const prevOf = (pnId) => (state.prev.embryos || []).find((e) => e.embryo_id === pnId);
  const nPN = (r) => (r && r.geometry && r.geometry.n_pronuclei) || 0;

  // ───────────────────────── top embryo bar ─────────────────────────
  function buildTabs() {
    // every scored zygote, ordered by tau; ones with transcripts are marked
    const list = scoredZygotes().slice().sort((a, b) => a.pseudotime.tau_mean - b.pseudotime.tau_mean);
    const tabsEl = $("#tabs");
    tabsEl.innerHTML = "";
    list.forEach((e) => {
      const row = txRowOf(e.embryo_id);
      const b = document.createElement("button");
      b.className = "tab"; b.dataset.id = e.embryo_id;
      const label = row ? row.label : e.embryo_id.split("/").slice(-2).join("/");
      const one = nPN(e) === 1;
      b.title = `${e.embryo_id} · τ≈${f(e.pseudotime.tau_mean, 2)}` +
        (row ? ` · ${row.total.toLocaleString()} transcripts` : " · no transcript data") +
        (one ? " · single pronucleus (extrapolated)" : "");
      b.innerHTML = `<span class="tab-label">${esc(label)}</span>` +
        `<span class="tab-date">τ ${f(e.pseudotime.tau_mean, 2)}${one ? " ⚠" : ""}</span>` +
        `<span class="tab-year">${row ? (row.total / 1000).toFixed(0) + "k tx" : "no tx"}</span>`;
      b.addEventListener("click", () => selectEmbryo(e.embryo_id));
      tabsEl.appendChild(b);
    });
  }

  // ───────────────────────── selection ─────────────────────────
  function selectEmbryo(pnId) {
    state.id = pnId;
    V.markActiveTab($("#tabs"), pnId);
    $("#placeholder").hidden = true;
    $("#controls").hidden = false;
    $("#drawer").hidden = false;
    const inf = infOf(pnId), row = txRowOf(pnId);
    $("#drawer-emb").textContent = "· " + (row ? row.label : pnId);
    renderMain(inf);
    renderSelected(inf, row);
    renderInputPanel(inf, row);
    renderScatter();
  }

  // main window: the embryo's structure, big
  function renderMain(inf) {
    const host = $("#plot-host");
    const p = prevOf(inf.embryo_id);
    host.innerHTML = "";
    const wrap = el("div", "tx-main-wrap");
    wrap.style.cssText = "padding:18px 22px;overflow:auto;height:100%";
    wrap.appendChild(el("div", "mini-title",
      `${esc(inf.embryo_id)} <span class="mini-sub">3-D structure · the model's view of this embryo</span>`));
    const planes = el("div", "pn-planes");
    [["xy", "XY (axial)"], ["xz", "XZ (coronal)"], ["yz", "YZ (sagittal)"]].forEach(([k, lab]) => {
      const d = el("div", "pn-plane");
      if (p && p.planes && p.planes[k]) {
        d.innerHTML = `<img src="${esc(p.planes[k].seg)}" alt="${esc(lab)}" loading="lazy">
                       <div class="pn-cap">${esc(lab)}</div>`;
      } else { d.className = "pn-plane pn-missing"; d.textContent = `${lab} — no preview`; }
      planes.appendChild(d);
    });
    wrap.appendChild(planes);
    wrap.appendChild(el("div", "tx-legend", legendHTML()));
    wrap.appendChild(el("div", "tx-note",
      `Open the drawer below to see <b>exactly what is fed into the pseudotime model</b> for this ` +
      `embryo — the images, the structures found in them, the dimensionless numbers derived, and the ` +
      `τ that comes back.`));
    host.appendChild(wrap);
  }
  const legendHTML = () =>
    `<span><span class="pn-sw" style="background:${COL.cell}"></span>cell body</span>
     <span><span class="pn-sw" style="background:${COL.pn1}"></span>pronucleus 1</span>
     <span><span class="pn-sw" style="background:${COL.pn2}"></span>pronucleus 2</span>
     <span><span class="pn-sw" style="background:${COL.polar}"></span>polar body</span>`;

  function renderSelected(inf, row) {
    const p = inf.pseudotime;
    const g = state.gene === "__total__" ? null : state.gene;
    const cnt = row ? (g ? (row.genes[g] || 0) : row.total) : null;
    $("#tx-selected").innerHTML =
      `<span class="tx-big">τ ${f(p.tau_mean, 3)}</span>
       <span class="tx-sub">95% ${f(p.interval_95[0], 2)} – ${f(p.interval_95[1], 2)} · ` +
      `${nPN(inf)} pronucle${nPN(inf) === 1 ? "us" : "i"} · ${esc(inf.ood_level)}</span>` +
      (row ? `<div style="margin-top:6px">${g ? esc(g) : "total"}: <b>${cnt.toLocaleString()}</b> transcripts</div>`
           : `<div style="margin-top:6px;color:#b45309">no transcript data for this embryo</div>`);
  }

  // ───────────────────────── drawer · MODEL INPUT ─────────────────────────
  function renderInputPanel(inf, row) {
    const p = prevOf(inf.embryo_id), g = inf.geometry || {}, pt = inf.pseudotime;

    $("#tx-pipe").innerHTML =
      `<span class="tx-stepbox on">DAPI z-stack</span><span class="tx-arrow">→</span>` +
      `<span class="tx-stepbox on">3-D segmentation</span><span class="tx-arrow">→</span>` +
      `<span class="tx-stepbox on">dimensionless geometry</span><span class="tx-arrow">→</span>` +
      `<span class="tx-stepbox on">calibrated clock</span><span class="tx-arrow">→</span>` +
      `<span class="tx-stepbox on">τ + interval</span>`;

    const strip = (hostId, which) => {
      const host = $(hostId); host.innerHTML = "";
      [["xy", "XY"], ["xz", "XZ"], ["yz", "YZ"]].forEach(([k, lab]) => {
        const d = el("div", "pn-plane");
        if (p && p.planes && p.planes[k]) {
          d.innerHTML = `<img src="${esc(p.planes[k][which])}" alt="${lab}" loading="lazy">
                         <div class="pn-cap">${lab}</div>`;
        } else { d.className = "pn-plane pn-missing"; d.textContent = `${lab} — none`; }
        host.appendChild(d);
      });
    };
    strip("#tx-planes-raw", "raw");
    strip("#tx-planes-seg", "seg");
    $("#tx-legend").innerHTML = legendHTML();

    const d = g.pron_distances_um || [];
    $("#tx-geom").innerHTML = [
      ["pronuclei found", nPN(inf) || "—"],
      ["cell radius R", g.cell_radius_um != null ? f(g.cell_radius_um, 1) + " µm" : "—"],
      ...d.map((x, i) => [`pronucleus ${i + 1} → centre`, f(x, 1) + " µm"]),
      ["rms → centre", g.rms_to_center_um != null ? f(g.rms_to_center_um, 1) + " µm" : "—"],
      ["<b>rms ÷ R — model input</b>", `<b>${f(g.rms_over_R, 3)}</b>`],
      ["pronuclear volume fraction", f(g.pron_vol_frac, 4)],
      ["polar body", g.polar_body_present ? "found (external)" : "not found"],
    ].map(([k, v]) => `<div class="k">${k}</div><div class="v">${v}</div>`).join("");

    $("#tx-tau").innerHTML = [
      ["τ (posterior mean)", f(pt.tau_mean, 3)],
      ["sd", f(pt.tau_sd, 3)],
      ["50% interval", `${f(pt.interval_50[0], 2)} – ${f(pt.interval_50[1], 2)}`],
      ["80% interval", `${f(pt.interval_80[0], 2)} – ${f(pt.interval_80[1], 2)}`],
      ["95% interval", `${f(pt.interval_95[0], 2)} – ${f(pt.interval_95[1], 2)}`],
      ["confidence", f(pt.confidence, 2)],
      ["domain", esc(inf.ood_level)],
    ].map(([k, v]) => `<div class="k">${k}</div><div class="v">${v}</div>`).join("");

    // little tau bar
    const W = 420, H = 54, x = (t) => 16 + t * (W - 32);
    const bar = (iv, y, col, h) => `<rect x="${x(iv[0])}" y="${y}" width="${Math.max(1, x(iv[1]) - x(iv[0]))}" height="${h}" rx="3" fill="${col}"/>`;
    $("#tx-tauplot").innerHTML = `<svg viewBox="0 0 ${W} ${H}" style="max-width:100%;margin-top:8px">
      <line x1="16" y1="${H - 14}" x2="${W - 16}" y2="${H - 14}" stroke="#e5e9f0"/>
      ${[0, 0.5, 1].map((t) => `<text x="${x(t)}" y="${H - 3}" font-size="9" fill="${MUTED}" text-anchor="middle">${t}</text>`).join("")}
      ${bar(pt.interval_95, 8, "rgba(124,58,237,.16)", 26)}
      ${bar(pt.interval_80, 13, "rgba(124,58,237,.28)", 16)}
      ${bar(pt.interval_50, 17, "rgba(124,58,237,.45)", 8)}
      <line x1="${x(pt.tau_mean)}" y1="4" x2="${x(pt.tau_mean)}" y2="${H - 14}" stroke="#7c3aed" stroke-width="2"/>
    </svg>`;

    let note = `The model never sees transcripts. It reads the <b>DAPI stack</b>, finds the cell body and ` +
      `pronuclei, reduces them to <b>one dimensionless number</b> (rms distance ÷ cell radius) and maps ` +
      `that through a clock calibrated on live-imaged embryos. Transcript counts are joined only ` +
      `afterwards, so the τ axis on the next tab is independent of the expression data plotted against it.`;
    if (pt.if_second_pronucleus_present) {
      const s = pt.if_second_pronucleus_present;
      note += `<div class="tx-warn"><b>Single annotated pronucleus.</b> This embryo is still scored — a ` +
        `computer-vision model should return an answer for atypical input — but the clock was calibrated ` +
        `on two-pronucleus geometry, so it is flagged out-of-domain. If a second pronucleus were present ` +
        `but unlabelled, τ would fall between <b>${f(s.tau_lo, 2)}</b> and <b>${f(s.tau_hi, 2)}</b>. The ` +
        `interval above covers clock noise only, not that annotation uncertainty.</div>`;
    }
    $("#tx-input-note").innerHTML = note;
  }

  // ───────────────────────── drawer · SCATTER ─────────────────────────
  function rows() {
    let rs = joined().filter((r) => state.showOOD || r.ood !== "out_of_domain");
    if (state.scope === "detected" && state.gene !== "__total__") rs = rs.filter((r) => state.gene in r.genes);
    return rs;
  }
  const yval = (r) => {
    const c = state.gene === "__total__" ? r.total : (r.genes[state.gene] || 0);
    return state.norm === "frac" ? (r.total ? c / r.total : 0) : c;
  };
  const yLabel = () => {
    const g = state.gene === "__total__" ? "total transcripts" : state.gene + " transcripts";
    return state.norm === "frac" ? g.replace("transcripts", "fraction of total") : g;
  };

  function renderScatter() {
    const rs = rows();
    const xs = rs.map((r) => r.tau), ys = rs.map(yval);
    const traces = [{
      type: "scatter", mode: "markers", name: "zygote", x: xs, y: ys,
      error_x: {
        type: "data", symmetric: false,
        array: rs.map((r) => r.hi95 - r.tau), arrayminus: rs.map((r) => r.tau - r.lo95),
        color: "rgba(15,23,42,.20)", thickness: 1.3, width: 0,
      },
      marker: {
        size: rs.map((r) => r.pn3d_id === state.id ? 16 : 11),
        color: rs.map((r) => viridisAt(r.tau)),
        line: { color: rs.map((r) => r.pn3d_id === state.id ? "#0f172a" : "#fff"),
                width: rs.map((r) => r.pn3d_id === state.id ? 3 : 1.4) },
        symbol: rs.map((r) => r.ood === "out_of_domain" ? "diamond" : "circle"),
      },
      customdata: rs.map((r) => r.pn3d_id),
      text: rs.map((r) => `${r.label}<br>τ ${f(r.tau, 3)} [${f(r.lo95, 2)}, ${f(r.hi95, 2)}]<br>` +
        `${state.gene === "__total__" ? "total" : state.gene} ` +
        `${(state.gene === "__total__" ? r.total : (r.genes[state.gene] || 0)).toLocaleString()}<br>${r.ood}`),
      hovertemplate: "%{text}<extra></extra>",
    }];
    let stat = `n = <b>${xs.length}</b> embryos`;
    if (state.fit !== "none" && xs.length >= 3) {
      const ty = ys.map((y) => state.fit === "log" ? Math.log(Math.max(y, 1e-9)) : y);
      const tx = state.fit === "power" ? xs.map((v) => Math.log(Math.max(v, 1e-6))) : xs;
      const L = linreg(tx, ty);
      if (L) {
        const lo = Math.min(...xs), hi = Math.max(...xs), fx = [], fy = [];
        for (let i = 0; i < 60; i++) {
          const x = lo + (hi - lo) * i / 59;
          const xe = state.fit === "power" ? Math.log(Math.max(x, 1e-6)) : x;
          let y = L.a + L.b * xe;
          if (state.fit === "log" || state.fit === "power") y = Math.exp(y);
          fx.push(x); fy.push(y);
        }
        traces.push({ type: "scatter", mode: "lines", x: fx, y: fy, name: "fit",
          line: { color: "#0f172a", width: 2 }, hoverinfo: "skip" });
        stat = `fit <b>${state.fit}</b> · r = <b>${f(L.r, 3)}</b> · R² = <b>${f(L.r * L.r, 3)}</b> · ` +
               `Spearman ρ = <b>${f(spearman(xs, ys), 3)}</b> · n = <b>${xs.length}</b>`;
      }
    }
    Plotly.react("tx-scatter", traces, {
      margin: { l: 70, r: 16, t: 10, b: 46 }, showlegend: false,
      xaxis: { title: { text: "pseudotime τ  (0 = pronuclear formation → 1 = NEBD)" }, range: [-0.05, 1.05], gridcolor: GRID, zeroline: false },
      yaxis: { title: { text: yLabel() }, gridcolor: GRID, zeroline: false },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      font: { family: "system-ui,-apple-system,sans-serif", size: 11.5, color: "#334155" },
    }, { responsive: true, displaylogo: false, displayModeBar: false });
    $("#tx-stat").innerHTML = stat +
      ` &nbsp;<span style="color:${MUTED}">— diamonds are out-of-domain; the outlined point is the ` +
      `selected embryo. Bars are the clock's 95% interval, which is why single-embryo ordering is uncertain.</span>`;
    $("#tx-scatter-sub").textContent =
      `· ${state.gene === "__total__" ? "all transcripts" : state.gene} · ${rows().length} embryos`;
    const gd = document.getElementById("tx-scatter");
    if (gd && !gd._txBound) {
      gd._txBound = true;
      gd.on("plotly_click", (ev) => {
        const id = ev && ev.points && ev.points[0] && ev.points[0].customdata;
        if (id) selectEmbryo(id);
      });
    }
  }

  function linreg(xs, ys) {
    const n = xs.length; if (n < 3) return null;
    const mx = xs.reduce((a, b) => a + b, 0) / n, my = ys.reduce((a, b) => a + b, 0) / n;
    let sxx = 0, syy = 0, sxy = 0;
    for (let i = 0; i < n; i++) { const dx = xs[i] - mx, dy = ys[i] - my; sxx += dx * dx; syy += dy * dy; sxy += dx * dy; }
    if (sxx <= 0 || syy <= 0) return null;
    return { a: my - (sxy / sxx) * mx, b: sxy / sxx, r: sxy / Math.sqrt(sxx * syy) };
  }
  function spearman(xs, ys) {
    const rk = (v) => {
      const o = v.map((x, i) => [x, i]).sort((p, q) => p[0] - q[0]), r = new Array(v.length);
      let i = 0;
      while (i < o.length) {
        let j = i; while (j + 1 < o.length && o[j + 1][0] === o[i][0]) j++;
        const avg = (i + j) / 2 + 1;
        for (let k = i; k <= j; k++) r[o[k][1]] = avg;
        i = j + 1;
      }
      return r;
    };
    const L = linreg(rk(xs), rk(ys)); return L ? L.r : null;
  }

  // ───────────────────────── drawer · GENES + JOIN ─────────────────────────
  function fillGenes() {
    const sel = $("#gene-select");
    sel.innerHTML = "";
    sel.appendChild(new Option("all transcripts (total)", "__total__"));
    state.tx.gene_stats.forEach((g) => sel.appendChild(new Option(`${g.gene} (in ${g.n_detected})`, g.gene)));
  }
  function renderGeneTable(filter) {
    const q = (filter || "").toLowerCase();
    const nKey = state.scope === "all" ? "n_all" : "n_detected";
    const rhoK = statKey("rho"), rK = statKey("r");
    const gs = state.tx.gene_stats.filter((g) => !q || g.gene.toLowerCase().includes(q))
      .slice().sort((a, b) => Math.abs(b[rhoK] ?? 0) - Math.abs(a[rhoK] ?? 0));
    const t = el("table", "pn-table");
    t.innerHTML = "<thead><tr><th>gene</th><th class='num'>n</th><th class='num'>ρ</th><th class='num'>r</th><th class='num'>detected</th></tr></thead>";
    const tb = el("tbody");
    gs.slice(0, 400).forEach((g) => {
      const tr = el("tr");
      tr.innerHTML = `<td>${esc(g.gene)}</td><td class="num">${g[nKey]}</td>
        <td class="num">${f(g[rhoK], 3)}</td><td class="num">${f(g[rK], 3)}</td>
        <td class="num">${g.n_detected}</td>`;
      tr.onclick = () => { state.gene = g.gene; $("#gene-select").value = g.gene; renderScatter(); renderSelectedSafe(); switchTab("scatter"); };
      tb.appendChild(tr);
    });
    t.appendChild(tb);
    const host = $("#tx-genetbl"); host.innerHTML = ""; host.appendChild(t);
    $("#tx-genes-sub").textContent = `· ${gs.length} genes · ${state.scope === "all" ? "all embryos" : "detected only"}`;
  }
  function renderSelectedSafe() {
    if (!state.id) return;
    renderSelected(infOf(state.id), txRowOf(state.id));
  }
  function renderJoin() {
    const j = state.tx.join, m = state.tx.model;
    $("#tx-join").innerHTML = [
      ["embryos joined", j.n_joined],
      ["clock-scored zygotes", j.n_pn3d_zygotes],
      ["transcript embryos", j.n_transcript_embryos],
      ["ambiguous pairings", j.n_ambiguous],
      ["centroid error (median)", j.centroid_error_um ? f(j.centroid_error_um.median, 3) + " µm" : "—"],
      ["centroid error (max)", j.centroid_error_um ? f(j.centroid_error_um.max, 3) + " µm" : "—"],
      ["clock", m.clock],
      ["clock CV MAE", f(m.cv_metrics.mae, 3)],
    ].map(([k, v]) => `<div class="k">${esc(k)}</div><div class="v">${esc(String(v))}</div>`).join("");
    $("#tx-join-note").innerHTML =
      `${esc(j.method)}<br><br>${esc(j.note)}<br><br><b>Time supervision:</b> ${esc(m.time_supervision)}`;
  }

  // ───────────────────────── wiring ─────────────────────────
  function wireControls() {
    $("#gene-select").onchange = (e) => { state.gene = e.target.value; renderScatter(); renderSelectedSafe(); };
    $("#norm-sel").onchange = (e) => { state.norm = e.target.value; renderScatter(); renderGeneTable($("#tx-gene-filter").value); };
    $("#fit-sel").onchange = (e) => { state.fit = e.target.value; renderScatter(); };
    $("#scope-sel").onchange = (e) => { state.scope = e.target.value; renderScatter(); renderGeneTable($("#tx-gene-filter").value); };
    $("#show-ood").onchange = (e) => { state.showOOD = e.target.checked; renderScatter(); };
    $("#tx-gene-filter").oninput = (e) => renderGeneTable(e.target.value);
    $("#tx-scatter-dl").onclick = () => {
      const rs = rows();
      const head = ["embryo", "pn3d_id", "tau", "lo95", "hi95", "n_pronuclei", "ood",
                    state.gene === "__total__" ? "total" : state.gene, "total"];
      const lines = [head.join(",")].concat(rs.map((r) => [
        r.label, r.pn3d_id, r.tau, r.lo95, r.hi95,
        (infOf(r.pn3d_id) || {}).n_pronuclei ?? "", r.ood,
        state.gene === "__total__" ? r.total : (r.genes[state.gene] || 0), r.total,
      ].join(",")));
      const b = new Blob([lines.join("\n")], { type: "text/csv" });
      const a = document.createElement("a");
      a.href = URL.createObjectURL(b);
      a.download = `transcripts_vs_tau_${state.gene}.csv`;
      a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1500);
    };
    if (V.wireWindow) { try { V.wireWindow($("#controls"), $("#controls-header")); } catch (_) { /* optional */ } }
  }

  function switchTab(which) {
    state.tab = which;
    $("#tx-tabs").querySelectorAll(".xs-gtab").forEach((t) => {
      const on = t.dataset.tab === which;
      t.classList.toggle("active", on); t.setAttribute("aria-selected", String(on));
    });
    $("#tx-panels").querySelectorAll(".xs-panel").forEach((p) => { p.hidden = p.dataset.tab !== which; });
    if (which === "scatter") requestAnimationFrame(() => { try { Plotly.Plots.resize($("#tx-scatter")); } catch (_) {} });
  }
  function wirePanels() {
    $("#tx-tabs").querySelectorAll(".xs-gtab").forEach((t) => t.onclick = () => switchTab(t.dataset.tab));
  }

  function wireDrawer() {
    const drawer = $("#drawer"), handle = $("#drawer-handle");
    const setOpen = (open) => {
      state.drawerOpen = open;
      drawer.dataset.open = open ? "true" : "false";
      handle.setAttribute("aria-expanded", String(open));
      if (open) requestAnimationFrame(() => { try { Plotly.Plots.resize($("#tx-scatter")); } catch (_) {} });
    };
    handle.addEventListener("click", () => setOpen(drawer.dataset.open !== "true"));
    let sh = 0; const rz = $("#drawer-resize");
    rz.addEventListener("pointerdown", (e) => {
      sh = $("#drawer-body").getBoundingClientRect().height;
      rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault();
    });
    rz.addEventListener("pointermove", (e) => {
      if (!rz._d) return;
      const h = Math.max(180, Math.min(window.innerHeight - 140, sh + (rz._d.y - e.clientY)));
      drawer.style.setProperty("--drawer-h", h + "px");
    });
    rz.addEventListener("pointerup", (e) => { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} });
    setOpen(true);           // the model-input panel is the point of the page — start open
  }
})();
