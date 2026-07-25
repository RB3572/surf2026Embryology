/* Pronuclear Structure Clock — GUIDED VISUAL NARRATIVE.
 *
 * Mirrors the Pronuclear Distance Clock's guided view (same CSS, same figure and
 * export idiom) but tells this model's story: a fixed 3-D stack → audited
 * segmentation → dimensionless geometry → a monotone probabilistic clock → tau
 * with calibrated intervals and QC.
 *
 * Every number is read from the frozen artifacts (data/pn3d/*.json). Nothing is
 * hard-coded and nothing is hidden: single-pronucleus and out-of-domain embryos
 * appear alongside the clean ones, marked.
 */
window.PN3DGuided = (() => {
  const $ = (s, r = document) => r.querySelector(s);
  const el = (t, c, h) => { const e = document.createElement(t); if (c) e.className = c; if (h != null) e.innerHTML = h; return e; };
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>]/g, (m) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[m]));
  const fmt = (v, d = 3) => (v == null || Number.isNaN(v)) ? "—" : (+v).toFixed(d);

  const VIRIDIS = [[0, "#440154"], [0.25, "#3b528b"], [0.5, "#21918c"], [0.75, "#5ec962"], [1, "#fde725"]];
  const viridisAt = (t) => {
    t = Math.max(0, Math.min(1, t));
    const hx = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
    for (let i = 1; i < VIRIDIS.length; i++) {
      if (t <= VIRIDIS[i][0]) {
        const a = VIRIDIS[i - 1], b = VIRIDIS[i], f = (t - a[0]) / (b[0] - a[0]);
        const ca = hx(a[1]), cb = hx(b[1]);
        return `rgb(${Math.round(ca[0] + (cb[0] - ca[0]) * f)},${Math.round(ca[1] + (cb[1] - ca[1]) * f)},${Math.round(ca[2] + (cb[2] - ca[2]) * f)})`;
      }
    }
    return VIRIDIS[VIRIDIS.length - 1][1];
  };
  const GRID = "#eef1f5", INK = "#0f172a", MUTED = "#64748b";
  const COL = { cell: "#96a0af", pn1: "#3c82f6", pn2: "#dc3c3c", polar: "#f5af3c" };
  const S = { model: null, inf: null, geo: null, prev: null, man: null, sel: null };

  // ───────────────────────── boot ─────────────────────────
  function init(data) {
    Object.assign(S, data);
    buildRail();
    buildViewSwitch();
    fillEmbryoSelect();
    drawTimeline();
    drawClock();
    drawValidation();
    drawCohort();
    drawLimits();
    selectEmbryo(defaultEmbryo());
    wireExports();
    $("#gd-seg").addEventListener("change", () => selectEmbryo(S.sel));
    $("#gd-embryo").addEventListener("change", (e) => selectEmbryo(e.target.value));
    ["gd-qc-two", "gd-qc-one"].forEach((id) => $("#" + id).addEventListener("change", drawCohort));
    window.addEventListener("resize", resizeAll);
  }

  const zygotes = () => S.inf.embryos.filter((e) => e.stage === "zygote");
  const scored = () => zygotes().filter((e) => e.inferable);
  const recOf = (id) => S.inf.embryos.find((e) => e.embryo_id === id);
  const prevOf = (id) => (S.prev.embryos || []).find((e) => e.embryo_id === id);
  const nPN = (r) => (r && r.geometry && r.geometry.n_pronuclei) || 0;

  function defaultEmbryo() {
    const withPrev = scored().filter((e) => prevOf(e.embryo_id) && nPN(e) === 2)
      .sort((a, b) => a.pseudotime.tau_mean - b.pseudotime.tau_mean);
    const pool = withPrev.length ? withPrev : scored();
    return (pool[Math.floor(pool.length / 2)] || zygotes()[0]).embryo_id;
  }

  function buildRail() {
    const steps = [["1", "τ scale"], ["2", "The stack"], ["3", "Structures"], ["4", "Geometry"],
                   ["5", "The clock"], ["6", "Validation"], ["7", "Our zygotes"], ["!", "Limits"]];
    const rail = $("#gd-rail");
    rail.innerHTML = "";
    steps.forEach(([n, label], i) => {
      const b = el("button", "gd-rb", `<b>${n}</b>${esc(label)}`);
      b.onclick = () => jumpTo("gd-" + (i + 1));
      rail.appendChild(b);
    });
    const main = $("#gd-main");
    main.addEventListener("scroll", () => {
      const bs = [...rail.querySelectorAll(".gd-rb")];
      let active = 0;
      steps.forEach((_, i) => {
        const s = $("#gd-" + (i + 1));
        if (s && s.getBoundingClientRect().top < window.innerHeight * 0.42) active = i;
      });
      bs.forEach((b, i) => b.classList.toggle("active", i === active));
    }, { passive: true });
  }
  function jumpTo(id) {
    const sc = $("#gd-main"), t = $("#" + id);
    if (sc && t) sc.scrollTop = t.getBoundingClientRect().top - sc.getBoundingClientRect().top + sc.scrollTop - 8;
  }

  function buildViewSwitch() {
    document.querySelectorAll(".pt-vb").forEach((b) => {
      b.onclick = () => {
        const v = b.dataset.view;
        document.querySelectorAll(".pt-vb").forEach((x) => x.classList.toggle("active", x === b));
        $("#gd-main").hidden = v !== "guided";
        $("#pt-main").hidden = v !== "advanced";
        $("#pt-nav").hidden = v !== "advanced";
        if (v === "guided") { $("#gd-main").scrollTop = 0; resizeAll(); }
        else if (window.PN3DAdvanced) window.PN3DAdvanced.show();
      };
    });
  }

  function fillEmbryoSelect() {
    const sel = $("#gd-embryo");
    const two = scored().filter((e) => nPN(e) === 2).sort((a, b) => a.pseudotime.tau_mean - b.pseudotime.tau_mean);
    const one = scored().filter((e) => nPN(e) === 1).sort((a, b) => a.pseudotime.tau_mean - b.pseudotime.tau_mean);
    const grp = (label, list) => {
      if (!list.length) return;
      const g = el("optgroup"); g.label = label;
      list.forEach((e) => {
        const o = el("option"); o.value = e.embryo_id;
        o.textContent = `${e.embryo_id}  ·  τ≈${fmt(e.pseudotime.tau_mean, 2)}` +
          (prevOf(e.embryo_id) ? "" : "  (no preview)");
        g.appendChild(o);
      });
      sel.appendChild(g);
    };
    sel.innerHTML = "";
    grp(`two pronuclei (${two.length})`, two);
    grp(`single pronucleus — extrapolated (${one.length})`, one);
  }

  // ───────────────────────── 1 · τ timeline ─────────────────────────
  function drawTimeline() {
    const W = 1120, H = 330, host = $("#fig-timeline");
    const traj = S.model.canonical_trajectory || [];
    const rmsAt = (t) => {
      if (!traj.length) return 0.8;
      let best = traj[0];
      for (const p of traj) if (Math.abs(p.tau - t) < Math.abs(best.tau - t)) best = p;
      return best.rms_over_R;
    };
    const stops = [];
    for (let i = 0; i <= 10; i++) stops.push(`<stop offset="${i * 10}%" stop-color="${viridisAt(i / 10)}"/>`);
    const x0 = 90, x1 = W - 60, ybar = 250;
    const marks = [[0, "τ = 0", "pronuclear formation"], [1 / 3, "0.33", "early"],
                   [2 / 3, "0.67", "late"], [1, "τ = 1", "NEBD"]];
    const cells = [0, 0.25, 0.5, 0.75, 1].map((t, i) => {
      const cx = 150 + i * ((W - 260) / 4), cy = 142, R = 52;
      const sep = Math.min(R * 0.86, rmsAt(t) * R * 0.92);
      return `<g>
        <circle cx="${cx}" cy="${cy}" r="${R}" fill="#eef2f7" stroke="#cbd5e1" stroke-width="1.5"/>
        <circle cx="${cx - sep}" cy="${cy}" r="20" fill="none" stroke="${COL.pn1}" stroke-width="2.4"/>
        <circle cx="${cx + sep}" cy="${cy}" r="20" fill="none" stroke="${COL.pn2}" stroke-width="2.4"/>
        <circle cx="${cx}" cy="${cy}" r="2.6" fill="${INK}"/>
        <text class="s-l" x="${cx}" y="${cy - R - 14}">τ ${t.toFixed(2)}</text>
        <circle cx="${cx}" cy="${cy + R + 20}" r="6" fill="${viridisAt(t)}"/></g>`;
    }).join("");
    host.innerHTML = `<svg class="gd-svg" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" role="img">
      <defs><linearGradient id="pn3dtau" x1="0" x2="1">${stops.join("")}</linearGradient></defs>
      <text class="s-t" x="60" y="34">Normalized pronuclear-migration time τ</text>
      <text class="s-m" x="60" y="56">one 0→1 scale for every zygote, so embryos of different absolute duration compare directly</text>
      <text class="s-b" x="${W - 60}" y="34" style="text-anchor:end">SCHEMATIC</text>
      <text class="s-m" x="${W - 60}" y="56" style="text-anchor:end">separation follows the cohort-median measured geometry</text>
      ${cells}
      <rect x="${x0}" y="${ybar}" width="${x1 - x0}" height="14" rx="7" fill="url(#pn3dtau)"/>
      ${marks.map(([t, a, b]) => {
        const x = x0 + t * (x1 - x0);
        return `<line x1="${x}" y1="${ybar - 8}" x2="${x}" y2="${ybar + 22}" stroke="#94a3b8"/>
                <text class="s-l" x="${x}" y="${ybar + 42}">${a}</text>
                <text class="s-s" x="${x}" y="${ybar + 58}">${b}</text>`;
      }).join("")}
    </svg>`;
    $("#cap-timeline").innerHTML =
      `τ is <b>normalized</b> time: 0 at pronuclear formation, 1 at nuclear-envelope breakdown. ` +
      `Zygotes differ in how many hours that takes, so absolute time is not comparable across embryos — ` +
      `τ is. The circles are a schematic drawn from the cohort-median measured geometry, not microscopy.`;
  }

  // ───────────────────────── 2 · the stack + 3 · structures ─────────────────────────
  function selectEmbryo(id) {
    if (!id) return;
    S.sel = id;
    const r = recOf(id); if (!r) return;
    if ($("#gd-embryo").value !== id) $("#gd-embryo").value = id;
    const g = r.geometry || {};
    const showSeg = $("#gd-seg").checked;
    const p = prevOf(id);

    $("#gd-emb-note").innerHTML =
      `${esc(id)} · <b>${nPN(r)}</b> pronucle${nPN(r) === 1 ? "us" : "i"} · ` +
      `segmentation ${qcPill(r)} · domain ${oodPill(r)}`;

    // --- planes ---
    const planes = el("div", "pn-planes");
    [["xy", "XY (axial)"], ["xz", "XZ (coronal)"], ["yz", "YZ (sagittal)"]].forEach(([k, lab]) => {
      const d = el("div", "pn-plane");
      if (p && p.planes && p.planes[k]) {
        const src = showSeg ? p.planes[k].seg : p.planes[k].raw;
        d.innerHTML = `<img src="${esc(src)}" alt="${esc(lab)}" loading="lazy"
            onerror="this.closest('.pn-plane').classList.add('pn-missing');this.closest('.pn-plane').textContent='preview not published'">
          <div class="pn-cap">${esc(lab)}</div>`;
      } else { d.className = "pn-plane pn-missing"; d.textContent = `${lab} — no preview`; }
      planes.appendChild(d);
    });
    const host = $("#fig-planes"); host.innerHTML = ""; host.appendChild(planes);
    host.appendChild(el("div", "pn-legend",
      `<span><span class="pn-sw" style="background:${COL.cell}"></span>cell body</span>
       <span><span class="pn-sw" style="background:${COL.pn1}"></span>pronucleus 1</span>
       <span><span class="pn-sw" style="background:${COL.pn2}"></span>pronucleus 2</span>
       <span><span class="pn-sw" style="background:${COL.polar}"></span>polar body</span>`));
    $("#cap-planes").innerHTML =
      `DAPI intensity projections (90th percentile through depth) in the three orthogonal planes — the whole of what the model ` +
      `sees for this embryo. There is no earlier or later frame to compare against, which is exactly ` +
      `why a single stack is developmentally ambiguous and why the geometry has to carry the answer.`;

    drawStructures(r);
    drawGeometry(r);
    drawClock();                     // re-mark this embryo on the curve
  }

  function qcPill(r) {
    const s = r.segmentation_status;
    if (s === "resolved") return `<span class="pn-pill pn-pill-ok">two pronuclei</span>`;
    if (s === "single_pronucleus") return `<span class="pn-pill pn-pill-warn">single pronucleus</span>`;
    return `<span class="pn-pill pn-pill-bad">${esc(s)}</span>`;
  }
  function oodPill(r) {
    const l = r.ood_level;
    if (l === "in_domain") return `<span class="pn-pill pn-pill-ok">in-domain</span>`;
    if (l === "caution") return `<span class="pn-pill pn-pill-warn">caution</span>`;
    return `<span class="pn-pill pn-pill-bad">out-of-domain</span>`;
  }

  function drawStructures(r) {
    const g = r.geometry || {};
    const p = prevOf(r.embryo_id);
    const host = $("#fig-struct"); host.innerHTML = "";
    if (p && p.planes && p.planes.xy) {
      host.innerHTML = `<img src="${esc(p.planes.xy.seg)}" alt="segmented embryo"
        style="width:100%;border-radius:10px;background:#eef2f7" loading="lazy">`;
    } else host.innerHTML = `<div class="gd-empty">No published preview for this embryo.</div>`;

    const has2 = nPN(r) === 2;
    const checks = [
      [has2 ? "ok" : "warn", has2 ? "Exactly two pronuclei inside the cell body"
        : "Only one pronucleus inside the cell body — measured anyway, and flagged"],
      [g.polar_body_present ? "ok" : "warn", g.polar_body_present
        ? "A polar body was found OUTSIDE the cell body" : "No external polar body found"],
      ["ok", "Polar body never substituted for a pronucleus (it fails the containment test)"],
      ["ok", "Identity assigned from geometry — the label file only says “Segment_1…4”"],
    ];
    $("#gd-constraints").innerHTML =
      `<div class="pn-checks">` + checks.map(([k, t]) =>
        `<div class="pn-check"><span class="pn-mark pn-${k}">${k === "ok" ? "✓" : "!"}</span><span>${esc(t)}</span></div>`
      ).join("") + `</div>`;
    $("#cap-struct").innerHTML =
      `The cell body is the largest segment; a pronucleus is a compact segment <b>enclosed</b> by it; ` +
      `the polar body sits outside, across the perivitelline gap. Containment is tested by filling the ` +
      `cell together with the candidate — which keeps a pronucleus that touches the cortex, and still ` +
      `rejects a detached polar body.`;
  }

  // ───────────────────────── 4 · geometry ─────────────────────────
  function drawGeometry(r) {
    const g = r.geometry || {};
    const R = g.cell_radius_um, d = g.pron_distances_um || [];
    const host = $("#fig-geom");
    const W = 560, H = 260, cx = W / 2, cy = H / 2 - 8, RR = 92;
    const scale = R ? RR / R : 1;
    const blobs = d.map((dist, i) => {
      const ang = d.length > 1 ? (i === 0 ? Math.PI : 0) : Math.PI * 0.75;
      const x = cx + Math.cos(ang) * dist * scale, y = cy + Math.sin(ang) * dist * scale;
      const c = i === 0 ? COL.pn1 : COL.pn2;
      return `<line x1="${cx}" y1="${cy}" x2="${x}" y2="${y}" stroke="${c}" stroke-width="2" stroke-dasharray="4 3"/>
              <circle cx="${x}" cy="${y}" r="17" fill="${c}" fill-opacity=".25" stroke="${c}" stroke-width="2"/>
              <text class="s-s" x="${x}" y="${y - 22}">d${i + 1} = ${fmt(dist, 1)} µm</text>`;
    }).join("");
    host.innerHTML = `<svg class="gd-svg" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" role="img">
      <circle cx="${cx}" cy="${cy}" r="${RR}" fill="#eef2f7" stroke="#cbd5e1" stroke-width="1.6"/>
      <line x1="${cx}" y1="${cy}" x2="${cx + RR}" y2="${cy}" stroke="#94a3b8" stroke-width="1.2"/>
      <text class="s-s" x="${cx + RR / 2}" y="${cy - 8}">R = ${fmt(R, 1)} µm</text>
      ${blobs}
      <circle cx="${cx}" cy="${cy}" r="3.4" fill="${INK}"/>
      <text class="s-s" x="${cx}" y="${cy + 20}">cell centre</text>
    </svg>`;
    const rms = g.rms_over_R;
    $("#gd-eq").innerHTML =
      `<div class="gd-eq-row"><span class="gd-eq-k">cell radius R</span><span class="gd-eq-v">${fmt(R, 1)} µm</span></div>` +
      d.map((x, i) => `<div class="gd-eq-row"><span class="gd-eq-k">pronucleus ${i + 1} → centre</span><span class="gd-eq-v">${fmt(x, 1)} µm</span></div>`).join("") +
      `<div class="gd-eq-row gd-eq-sum"><span class="gd-eq-k">rms ÷ R <b>(model input)</b></span><span class="gd-eq-v">${fmt(rms, 3)}</span></div>`;
    $("#cap-geom").innerHTML =
      `The model input is the <b>root-mean-square</b> pronucleus-to-centre distance divided by the ` +
      `embryo's own cell radius. Dividing by R makes it dimensionless, so it cannot be thrown off by ` +
      `microns-per-pixel or by a larger embryo. Using the rms (rather than a sum of exactly two ` +
      `distances) is what lets a one-pronucleus embryo be scored on the same scale — and on ` +
      `two-pronucleus embryos it is the identical quantity the clock was trained on.`;
  }

  // ───────────────────────── 5 · the clock ─────────────────────────
  function drawClock() {
    const c = S.model.clock;
    const xs = c.iso_x, ys = c.iso_y;
    const traces = [{
      type: "scatter", mode: "lines", x: xs, y: ys, name: "clock",
      line: { color: INK, width: 3, shape: "hv" },
      hovertemplate: "rms/R %{x:.2f} → τ %{y:.3f}<extra></extra>",
    }];
    const r = S.sel && recOf(S.sel);
    if (r && r.geometry && r.pseudotime) {
      const x = r.geometry.rms_over_R, t = r.pseudotime.tau_mean;
      traces.push({
        type: "scatter", mode: "lines", x: [x, x], y: [0, t],
        line: { color: viridisAt(t), width: 2, dash: "dot" }, hoverinfo: "skip", showlegend: false,
      }, {
        type: "scatter", mode: "lines", x: [xs[0], x], y: [t, t],
        line: { color: viridisAt(t), width: 2, dash: "dot" }, hoverinfo: "skip", showlegend: false,
      }, {
        type: "scatter", mode: "markers", x: [x], y: [t], name: "this embryo",
        marker: { size: 13, color: viridisAt(t), symbol: "diamond", line: { color: "#fff", width: 2 } },
        hovertemplate: `${esc(r.embryo_id)}<br>rms/R %{x:.3f} → τ %{y:.3f}<extra></extra>`,
      });
    }
    plot("fig-curve", traces, {
      margin: { l: 62, r: 14, t: 8, b: 46 }, showlegend: false,
      xaxis: { title: { text: "rms pronucleus→centre distance ÷ cell radius  (larger = earlier)" }, gridcolor: GRID, zeroline: false },
      yaxis: { title: { text: "τ" }, range: [-0.03, 1.03], gridcolor: GRID, zeroline: false },
    });
    $("#cap-curve").innerHTML =
      `A <b>monotone</b> (isotonic) curve: more separation always means earlier, never the reverse. ` +
      `Nothing else about the shape is assumed. It is fit only to the live-imaged cohort, where the ` +
      `true time is known. To read a fixed embryo: go up from its measured rms/R, then across to τ.`;
  }

  // ───────────────────────── 6 · validation ─────────────────────────
  function drawValidation() {
    const cv = S.model.clock.cv_metrics;
    $("#gd-kpis").innerHTML = [
      [fmt(cv.mae, 3), "average error in τ", `${cv.n_embryos} embryos held out one at a time`],
      [fmt(cv.spearman, 3), "rank correlation ρ", "ordering agreement with true τ"],
      [fmt(cv.within_embryo_mono_median, 3), "within-embryo monotonicity", "median ρ inside a trajectory"],
      [`${(cv.coverage_95 * 100).toFixed(0)}%`, "95% interval coverage", "stated 95% — measured on held-out frames"],
    ].map(([n, t, s]) => `<div class="gd-kpi"><div class="gd-kpi-n">${n}</div>
        <div class="gd-kpi-t">${esc(t)}</div><div class="gd-kpi-s">${esc(s)}</div></div>`).join("");

    const sc = S.model.clock.cv_scatter || [];
    plot("fig-heldout", [
      { type: "scatter", mode: "lines", x: [0, 1], y: [0, 1], line: { color: MUTED, width: 1.4, dash: "dash" }, hoverinfo: "skip" },
      { type: "scatter", mode: "markers", x: sc.map((p) => p.true), y: sc.map((p) => p.pred),
        marker: { size: 4, color: sc.map((p) => viridisAt(p.true)), opacity: 0.55 },
        hovertemplate: "true τ %{x:.3f}<br>predicted τ %{y:.3f}<extra></extra>" },
    ], {
      margin: { l: 54, r: 12, t: 24, b: 44 }, showlegend: false,
      title: { text: "held-out predicted vs true τ", font: { size: 11.5 }, x: 0, xanchor: "left" },
      xaxis: { title: { text: "true τ" }, range: [-0.03, 1.03], gridcolor: GRID, zeroline: false },
      yaxis: { title: { text: "predicted τ" }, range: [-0.03, 1.03], gridcolor: GRID, zeroline: false },
    });

    const cc = S.model.clock.calibration_curve || [];
    plot("fig-calib", [
      { type: "scatter", mode: "lines", x: [0.4, 1], y: [0.4, 1], line: { color: MUTED, width: 1.4, dash: "dash" }, hoverinfo: "skip" },
      { type: "scatter", mode: "lines+markers", x: cc.map((c) => c.nominal), y: cc.map((c) => c.empirical),
        line: { color: "#4338ca", width: 2.4 }, marker: { size: 7, color: "#4338ca" },
        hovertemplate: "stated %{x:.0%}<br>actual %{y:.1%}<extra></extra>" },
    ], {
      margin: { l: 54, r: 12, t: 24, b: 44 }, showlegend: false,
      title: { text: "calibration — stated vs actual coverage", font: { size: 11.5 }, x: 0, xanchor: "left" },
      xaxis: { title: { text: "stated coverage" }, gridcolor: GRID, zeroline: false },
      yaxis: { title: { text: "actual" }, gridcolor: GRID, zeroline: false },
    });
    $("#cap-valid").innerHTML =
      `Each embryo is held out in turn, so no embryo is ever tested by a clock that saw it. The right-hand ` +
      `plot is the honesty check: a stated 80% interval should contain the truth about 80% of the time, and ` +
      `it lands on the diagonal. <b>This is the only cohort with real time</b>; the fixed embryos in step 7 ` +
      `have no ground truth to be scored against.`;
  }

  // ───────────────────────── 7 · every zygote ─────────────────────────
  function drawCohort() {
    const showTwo = $("#gd-qc-two").checked, showOne = $("#gd-qc-one").checked;
    const rows = scored().filter((e) => (nPN(e) === 2 ? showTwo : showOne))
      .sort((a, b) => a.pseudotime.tau_mean - b.pseudotime.tau_mean);
    const two = scored().filter((e) => nPN(e) === 2).length;
    const one = scored().filter((e) => nPN(e) === 1).length;
    $("#gd-cohort-note").innerHTML =
      `<b>${scored().length}</b> of ${zygotes().length} zygotes scored — ${two} with two pronuclei, ` +
      `${one} with a single pronucleus (extrapolated, shown hollow)`;

    const H = Math.max(420, 15 * rows.length + 70);
    const host = $("#fig-cohort"); host.style.height = H + "px";
    plot("fig-cohort", [{
      type: "scatter", mode: "markers", x: rows.map((r) => r.pseudotime.tau_mean),
      y: rows.map((r) => r.embryo_id),
      error_x: {
        type: "data", symmetric: false,
        array: rows.map((r) => r.pseudotime.interval_95[1] - r.pseudotime.tau_mean),
        arrayminus: rows.map((r) => r.pseudotime.tau_mean - r.pseudotime.interval_95[0]),
        color: "rgba(15,23,42,.26)", thickness: 1.2, width: 0,
      },
      marker: {
        size: 9,
        color: rows.map((r) => nPN(r) === 2 ? viridisAt(r.pseudotime.tau_mean) : "rgba(255,255,255,0.9)"),
        line: { color: rows.map((r) => nPN(r) === 2 ? "#ffffff" : "#b45309"), width: 2 },
      },
      text: rows.map((r) => `${r.embryo_id}<br>τ ${fmt(r.pseudotime.tau_mean)} ` +
        `[${fmt(r.pseudotime.interval_95[0], 2)}, ${fmt(r.pseudotime.interval_95[1], 2)}]` +
        `<br>${nPN(r)} pronucle${nPN(r) === 1 ? "us" : "i"} · ${r.ood_level}`),
      hovertemplate: "%{text}<extra></extra>",
      customdata: rows.map((r) => r.embryo_id),
    }], {
      margin: { l: 210, r: 24, t: 8, b: 46 }, height: H, showlegend: false,
      xaxis: { title: { text: "estimated τ · bars = calibrated 95% interval" }, range: [-0.04, 1.04], gridcolor: GRID, zeroline: false },
      yaxis: { automargin: true, type: "category", tickfont: { size: 8.5 } },
    });
    const gd = document.getElementById("fig-cohort");
    if (gd && !gd._pnBound) {
      gd._pnBound = true;
      gd.on("plotly_click", (ev) => {
        const id = ev && ev.points && ev.points[0] && ev.points[0].customdata;
        if (id) { selectEmbryo(id); jumpTo("gd-2"); }
      });
    }
    $("#cap-cohort").innerHTML =
      `Every zygote the model can see, ordered by estimated τ — click any row to load it above. ` +
      `Filled markers have two pronuclei; <b>hollow amber</b> markers had a single annotated pronucleus, ` +
      `so their τ extrapolates past the configuration the clock was calibrated on and they are flagged ` +
      `out-of-domain. None of them are hidden, and none of them have a known true τ to check against.`;
  }

  // ───────────────────────── 8 · limits ─────────────────────────
  function drawLimits() {
    const c = S.model.counts;
    const lims = [
      ["No fixed embryo has a measured time.", "τ here is a model estimate. The clock is validated " +
        "only on the live-imaged cohort; agreement with any earlier pseudotime is not independent proof."],
      ["Single-pronucleus embryos are extrapolations.", `${c.zygotes_single_pronucleus || 0} zygotes ` +
        "had one annotated pronucleus. They are scored so nothing is dropped, but the clock was " +
        "calibrated on two-pronucleus geometry — the page also shows what τ would be if the missing " +
        "pronucleus were present."],
      ["The segmentation is given, not learned.", "Structures come from the project's 3-D Slicer " +
        "labels, audited geometrically. A learned volumetric segmenter is future work."],
      ["Physical scale is a display convention.", "Model features are dimensionless (÷ cell radius), " +
        "so microns affect only what is printed, never the estimate."],
      ["Ordering is stronger than any single value.", "Intervals overlap between neighbouring embryos; " +
        "widely separated embryos can be compared, adjacent ones should not be."],
    ];
    $("#gd-lims").innerHTML = lims.map(([t, b]) =>
      `<div class="gd-lim"><b>${esc(t)}</b><span>${esc(b)}</span></div>`).join("");
  }

  // ───────────────────────── plotting + export ─────────────────────────
  function plot(id, traces, layout) {
    const gd = document.getElementById(id); if (!gd) return;
    Plotly.react(gd, traces, Object.assign({
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
      font: { family: "system-ui,-apple-system,sans-serif", size: 11.5, color: "#334155" },
    }, layout), { responsive: true, displaylogo: false, displayModeBar: false });
  }
  function resizeAll() {
    ["fig-curve", "fig-heldout", "fig-calib", "fig-cohort"].forEach((id) => {
      const gd = document.getElementById(id);
      if (gd && gd.data) { try { Plotly.Plots.resize(gd); } catch (_) { /* not yet drawn */ } }
    });
  }

  const FIGS = {
    "01_tau_timeline": { svg: "#fig-timeline" },
    "02_stack": { dom: "#fig-planes" },
    "03_structures": { dom: "#fig-struct" },
    "04_geometry": { svg: "#fig-geom" },
    "05_clock_curve": { plotly: "fig-curve" },
    "06_validation": { plotly: "fig-heldout" },
    "07_cohort": { plotly: "fig-cohort" },
  };
  const SLIDE = { w: 1600, h: 900 };
  const SVG_CSS =
    ".s-t{font:600 22px system-ui,sans-serif;fill:#0f172a}" +
    ".s-m{font:400 14px system-ui,sans-serif;fill:#64748b}" +
    ".s-l{font:600 13px system-ui,sans-serif;fill:#334155;text-anchor:middle}" +
    ".s-s{font:400 12px system-ui,sans-serif;fill:#64748b;text-anchor:middle}" +
    ".s-b{font:700 13px system-ui,sans-serif;fill:#b45309;letter-spacing:.08em}";

  function wireExports() {
    document.querySelectorAll(".gd-dl").forEach((d) => {
      const key = d.dataset.fig;
      if (!FIGS[key]) { d.remove(); return; }
      d.innerHTML = `<button class="gd-dlb" data-k="${key}" data-f="svg" title="Download as SVG">SVG</button>` +
                    `<button class="gd-dlb" data-k="${key}" data-f="png" title="Download as PNG">PNG</button>`;
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
  function exportFig(key, kind, btn) {
    const spec = FIGS[key]; if (!spec) return;
    btn.classList.add("busy");
    const done = () => setTimeout(() => btn.classList.remove("busy"), 400);
    if (spec.svg) {
      const src = document.querySelector(spec.svg + " svg"); if (!src) return done();
      const c = src.cloneNode(true);
      c.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      c.setAttribute("width", SLIDE.w); c.setAttribute("height", SLIDE.h);
      const st = document.createElementNS("http://www.w3.org/2000/svg", "style");
      st.textContent = SVG_CSS; c.insertBefore(st, c.firstChild);
      const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      bg.setAttribute("width", "100%"); bg.setAttribute("height", "100%"); bg.setAttribute("fill", "#ffffff");
      c.insertBefore(bg, st.nextSibling);
      const txt = new XMLSerializer().serializeToString(c);
      if (kind === "svg") { saveBlob(new Blob([txt], { type: "image/svg+xml" }), key + ".svg"); return done(); }
      const img = new Image();
      img.onload = () => {
        const cv = document.createElement("canvas");
        cv.width = SLIDE.w * 2; cv.height = SLIDE.h * 2;
        const g = cv.getContext("2d"); g.fillStyle = "#fff"; g.fillRect(0, 0, cv.width, cv.height);
        g.drawImage(img, 0, 0, cv.width, cv.height);
        cv.toBlob((b) => { saveBlob(b, key + ".png"); done(); });
      };
      img.onerror = done;
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(txt);
      return;
    }
    if (spec.dom) {            // image strip — hand the user the underlying PNGs
      const imgs = document.querySelectorAll(spec.dom + " img");
      imgs.forEach((im, i) => {
        const a = document.createElement("a");
        a.href = im.getAttribute("src"); a.download = `${key}_${i}.png`;
        document.body.appendChild(a); a.click(); a.remove();
      });
      return done();
    }
    const gd = document.getElementById(spec.plotly);
    if (!gd || !gd.data) return done();
    const layout = Object.assign({}, gd.layout, {
      paper_bgcolor: "#ffffff", plot_bgcolor: "#ffffff", width: SLIDE.w, height: SLIDE.h,
      font: Object.assign({}, gd.layout.font, { size: 18 }),
    });
    Plotly.toImage({ data: gd.data, layout, config: { displayModeBar: false } },
      { format: kind === "svg" ? "svg" : "png", width: SLIDE.w, height: SLIDE.h, scale: kind === "svg" ? 1 : 2 })
      .then((url) => {
        if (kind === "svg") {
          saveBlob(new Blob([decodeURIComponent(url.replace(/^data:image\/svg\+xml,/, ""))],
            { type: "image/svg+xml" }), key + ".svg");
        } else {
          const a = document.createElement("a"); a.href = url; a.download = key + ".png";
          document.body.appendChild(a); a.click(); a.remove();
        }
      }).catch(() => {}).then(done);
  }

  return { init, selectEmbryo, resizeAll, jumpTo, viridisAt, qcPill, oodPill, plot };
})();
