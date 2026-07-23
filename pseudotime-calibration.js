/* Pronuclear Pseudotime Calibration — replay of the offline training/validation artifacts
 * (scripts/train_pronuclear_pseudotime.py) plus the frozen model's application to the fixed
 * MERFISH cohort (build_pronuclei_pseudotime.py).
 *
 * Everything shown here is PRECOMPUTED. The "training replay" steps through grouped folds that
 * were fitted offline; no model is fitted in the browser. Three cohorts/quantities are kept
 * visually and structurally separate throughout:
 *   1. training/validation  — Scheffler 2021 live-imaged control zygotes (have ground-truth tau)
 *   2. application          — this project's fixed MERFISH zygotes (no ground truth)
 *   3. legacy               — the pre-existing minimum surface-gap score (a different quantity) */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const C = { train: "#4f46e5", apply: "#0891b2", legacy: "#94a3b8", ok: "#15803d",
              warn: "#b45309", bad: "#b42318", ink: "#1a2233", grid: "#eef1f5" };
  const FOLD_C = ["#4f46e5", "#0891b2", "#0d9488", "#b45309", "#be185d"];
  const SECTIONS = [
    ["summary", "Summary"], ["replay", "Nested CV replay"], ["oof", "Predicted vs true"],
    ["traj", "Trajectories"], ["models", "Model comparison"], ["diag", "Diagnostics"],
    ["noise", "R² attenuation"], ["fixed", "Fixed cohort"], ["limits", "Limitations"],
  ];
  const ORD = (m) => (m.pooled_ordering || {}).strict_accuracy;   // explicit strict pair accuracy

  const state = { cal: null, frames: null, traj: null, fixed: null, noise: null,
                  sec: "summary", rpModel: null, rpFold: 0, playing: false, timer: 0, drawn: {} };

  const fmt = (v, d = 3) => (v == null || !isFinite(v) ? "—" : Number(v).toFixed(d));
  const pct = (v, d = 1) => (v == null || !isFinite(v) ? "—" : (100 * v).toFixed(d) + "%");
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

  const LAY = (over) => Object.assign({
    margin: { l: 58, r: 14, t: 10, b: 46 }, paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    font: { color: C.ink, size: 11.5 }, hovermode: "closest", showlegend: false,
    xaxis: { gridcolor: C.grid, zeroline: false }, yaxis: { gridcolor: C.grid, zeroline: false },
  }, over || {});
  const CFG = { responsive: true, displaylogo: false, displayModeBar: false };

  function plot(id, traces, layout) {
    const el = $("#" + id); if (!el) return;
    Plotly.react(el, traces, LAY(layout), CFG);
  }
  function downloadText(text, filename, mime) {
    const blob = new Blob([text], { type: mime || "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = filename;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }

  // ═════════════════════════════════════ init ═════════════════════════════════════
  (async function init() {
    try {
      const [cal, frames, traj, fixed, noise] = await Promise.all([
        (await fetch("data/pseudotime_calibration/calibration.json")).json(),
        V.loadGz("data/pseudotime_calibration/oof_frames.json.gz"),
        V.loadGz("data/pseudotime_calibration/trajectories.json.gz"),
        fetch("data/pronuclei_pseudotime.json").then((r) => (r.ok ? r.json() : null)).catch(() => null),
        fetch("data/pseudotime_calibration/noise_ceiling.json").then((r) => (r.ok ? r.json() : null)).catch(() => null),
      ]);
      state.cal = cal; state.frames = frames.frames; state.traj = traj.embryos; state.fixed = fixed;
      // a stale artifact must never be shown as if it matched the frozen clock
      state.noise = (noise && noise.meta.model_version === cal.meta.model_version) ? noise : null;
      state.noiseStale = !!(noise && noise.meta.model_version !== cal.meta.model_version);
      state.rpModel = cal.production.key;
      buildNav();
      $("#pt-subtitle").textContent =
        `${cal.dataset.n_embryos} live-imaged zygotes · ${cal.dataset.n_frames.toLocaleString()} frames · model ${cal.meta.model_version}`;
      renderSummary(); renderLimits();
      show("summary");
    } catch (e) {
      $("#pt-main").innerHTML =
        `<div class="pt-empty"><b>Failed to load calibration artifacts.</b><br>${esc(e.message || e)}
         <br><br>Run <code>scripts/train_pronuclear_pseudotime.py</code> then
         <code>build_pronuclei_pseudotime.py</code>.</div>`;
    }
  })();

  function buildNav() {
    const nav = $("#pt-nav");
    nav.innerHTML = SECTIONS.map(([k, l]) =>
      `<button class="pt-navb" data-sec="${k}">${l}</button>`).join("");
    nav.addEventListener("click", (e) => {
      const b = e.target.closest(".pt-navb"); if (b) show(b.dataset.sec);
    });
  }
  function show(sec) {
    state.sec = sec; stopPlay();
    document.querySelectorAll(".pt-navb").forEach((b) => b.classList.toggle("active", b.dataset.sec === sec));
    document.querySelectorAll(".pt-sec").forEach((s) => { s.hidden = s.dataset.sec !== sec; });
    const R = { replay: renderReplay, oof: renderOof, traj: renderTraj, models: renderModels,
                diag: renderDiag, noise: renderNoise, fixed: renderFixed };
    if (R[sec] && !state.drawn[sec]) { R[sec](); state.drawn[sec] = true; }
    resizeIn(sec);
    window.scrollTo({ top: 0 });
  }
  // Plotly cannot measure a display:none container, so a plot drawn at one width keeps that width
  // when its section is hidden during a viewport change. Resize on show (after layout settles) and
  // on window resize, for the visible section only.
  function resizeIn(sec) {
    const el = $("#sec-" + sec); if (!el) return;
    requestAnimationFrame(() => el.querySelectorAll(".pt-plot").forEach((p) => {
      try { if (p.offsetParent) Plotly.Plots.resize(p); } catch (_) {}
    }));
  }
  let rsTimer = 0;
  window.addEventListener("resize", () => {
    clearTimeout(rsTimer);
    rsTimer = setTimeout(() => resizeIn(state.sec), 120);
  });

  // ══════════════════════════════════ 1 · summary ══════════════════════════════════
  function renderSummary() {
    const c = state.cal, d = c.dataset, s = c.production;
    const m = c.nested_evaluation.outer_test_metrics, u = c.uncertainty;
    const ci = c.nested_evaluation.bootstrap_ci95_by_embryo;

    $("#sum-train-title").textContent = "Scheffler et al. 2021 — live-imaged control zygotes";
    $("#sum-train-desc").innerHTML =
      `${d.n_embryos} untreated zygotes · ${d.n_frames.toLocaleString()} frames · every ${d.frame_interval_h} h · ` +
      `migration ${d.migration_duration_h.min}–${d.migration_duration_h.max} h. ` +
      `<b>Has ground-truth τ</b> from real timestamps, so it is the only cohort that can validate the clock.`;
    const fx = state.fixed;
    $("#sum-apply-title").textContent = "This project's fixed MERFISH zygotes";
    $("#sum-apply-desc").innerHTML = fx
      ? `${fx.meta.n_predicted} of ${fx.meta.n_total} zygotes with two detected pronuclei. ` +
        `<b>No ground-truth time</b> — τ here is a model output carrying the held-out uncertainty, not a measurement.`
      : `Not built yet — run <code>build_pronuclei_pseudotime.py</code>.`;

    const kpi = (n, k, cls, sub) => `<div class="pt-kpi ${cls || ""}"><div class="pt-kpi-n">${n}</div>` +
      `<div class="pt-kpi-k">${k}</div>${sub ? `<div class="pt-kpi-ci">${sub}</div>` : ""}</div>`;
    const ciTxt = (a) => (a ? `95% CI ${fmt(a[0])} – ${fmt(a[1])}` : "");
    const cov = u.coverage_frame_marginal;
    $("#sum-kpis").innerHTML =
      kpi(fmt(m.macro_mae), "macro MAE<br><b>nested outer-test</b>", m.macro_mae <= 0.10 ? "pt-kpi-ok" : "pt-kpi-warn", ciTxt(ci.macro_mae)) +
      kpi(fmt(m.median_ae), "median absolute error", m.median_ae <= 0.10 ? "pt-kpi-ok" : "") +
      kpi(fmt(m.pooled_spearman), "Spearman ρ", m.pooled_spearman >= 0.8 ? "pt-kpi-ok" : "pt-kpi-warn", ciTxt(ci.pooled_spearman)) +
      kpi(pct(ORD(m)), "strict pairwise ordering<br><span class=\"pt-kpi-nb\">explicit pair counting</span>", ORD(m) >= 0.8 ? "pt-kpi-ok" : "pt-kpi-warn", ciTxt(ci.pooled_strict_ordering)) +
      kpi(pct(cov.mean), "95% coverage<br><span class=\"pt-kpi-nb\">measured on separate embryos</span>", Math.abs(cov.mean - 0.95) <= 0.05 ? "pt-kpi-ok" : "pt-kpi-warn", `95% CI ${pct(cov.ci95[0])} – ${pct(cov.ci95[1])}`) +
      kpi("±" + fmt(u.halfwidth_mean), "95% interval half-width", u.halfwidth_mean > 0.2 ? "pt-kpi-warn" : "");

    $("#sum-warn").innerHTML =
      `<b>Endpoint-normalized volumes are excluded from every deployable model.</b> The published
       relative pronuclear volumes are each divided by that pronucleus's volume at the <i>end of its
       own trajectory</i>, so they encode elapsed time and cannot be measured in a fixed snapshot.
       Including them would cut MAE to ≈0.05, but that number is a
       <b>non-deployable upper bound</b> only — it is reported in Model comparison, visually
       separated, and is structurally barred from selection. Male/female-specific distances are also
       excluded, because fixed zygotes have no reliable pronuclear identity call.`;

    const row = (k, v) => `<tr><td class="pt-k">${k}</td><td class="pt-wrap">${v}</td></tr>`;
    $("#sum-dataset").innerHTML =
      row("Source", `Scheffler et al. 2021, Nat Commun 12:841 — public source data`) +
      row("Data version", `<span class="pt-mono">${esc(d.data_version)}</span>`) +
      row("SHA-256", `<span class="pt-mono">${esc(d.sha256.slice(0, 24))}…</span>`) +
      row("Embryos / frames", `${d.n_embryos} / ${d.n_frames.toLocaleString()}`) +
      row("Frames per embryo", `${d.frames_per_embryo.min}–${d.frames_per_embryo.max}`) +
      row("τ observed range", `${d.tau_range[0]} – ${d.tau_range[1]}`) +
      row("Split", `${esc(c.meta.split_strategy)} — <b>${c.meta.n_outer_folds} outer × ` +
        `${c.meta.n_inner_folds} inner</b>, every embryo entirely within one group at every stage`) +
      row("Validation", d.validation.checks_passed.map((x) => `✓ ${esc(x)}`).join("<br>"));

    $("#sum-model").innerHTML =
      row("Production model", `<b>${esc(s.label)}</b>`) +
      row("Model class", esc(c.meta.model_class)) +
      row("Model version", `<span class="pt-mono">${esc(c.meta.model_version)}</span>`) +
      row("Features", s.features.map((f) => `<span class="pt-mono">${esc(f)}</span>`).join("<br>")) +
      row("Selection rule", esc(c.meta.selection_rule)) +
      row("Uncertainty", esc(u.method)) +
      row("Interval type", esc(u.interval_type)) +
      row("Seed", `<span class="pt-mono">${c.meta.seed}</span>`) +
      row("Trained (UTC)", `<span class="pt-mono">${esc(c.meta.trained_at_utc)}</span>`) +
      row("scikit-learn", `<span class="pt-mono">${esc(c.meta.sklearn_version)}</span>`);

    const deployable = new Set(s.features);
    const all = [
      ["nearer_to_center_um", "nearer pronucleus → cell centre"],
      ["farther_to_center_um", "farther pronucleus → cell centre"],
      ["distance_sum_um", "sum of the two distances"],
      ["distance_difference_um", "|difference| of the two distances"],
      ["male_to_center_um", "male pronucleus → cell centre"],
      ["female_to_center_um", "female pronucleus → cell centre"],
      ["male_relative_volume", "male volume ÷ its own endpoint volume"],
      ["female_relative_volume", "female volume ÷ its own endpoint volume"],
      ["volume_sum", "sum of relative volumes"],
      ["volume_difference", "|difference| of relative volumes"],
    ];
    const why = (f) =>
      f.startsWith("male_") || f.startsWith("female_")
        ? (f.includes("volume") ? "leaks the future endpoint AND needs identity" : "needs pronuclear identity — unavailable in fixed zygotes")
        : f.includes("volume") ? "normalized to a future endpoint — leaks elapsed time" : "identity-free geometry available in a fixed snapshot";
    $("#sum-features").innerHTML =
      `<tr><th>Feature</th><th>Deployable</th><th>In selected model</th><th>Reason</th></tr>` +
      all.map(([f, desc]) => {
        const ok = !/volume|male_|female_/.test(f);
        return `<tr><td><span class="pt-mono">${esc(f)}</span><br><span class="pt-note">${esc(desc)}</span></td>` +
          `<td><span class="pt-pill ${ok ? "pt-pill-yes" : "pt-pill-no"}">${ok ? "yes" : "no"}</span></td>` +
          `<td>${deployable.has(f) ? '<span class="pt-pill pt-pill-sel">selected</span>' : "—"}</td>` +
          `<td class="pt-wrap pt-note">${esc(why(f))}</td></tr>`;
      }).join("");
  }

  // ══════════════════════════ 2 · nested CV replay ══════════════════════════
  // Honest replay of the NESTED design: the model selector shows which candidate the INNER CV
  // picked for each outer fold; the outer-test embryos were never seen by that selection.
  function renderReplay() {
    const c = state.cal;
    $("#rp-fold").innerHTML = c.nested_evaluation.per_outer_fold.map((f) =>
      `<option value="${f.outer_fold}">Outer fold ${f.outer_fold + 1} of ${c.meta.n_outer_folds} · ` +
      `${f.n_test_embryos} embryos untouched · inner CV chose ${esc(f.chosen_by_inner_cv)}</option>`).join("");
    const sel = $("#rp-model");
    sel.innerHTML = `<option value="inner">Inner-CV choice for this fold (what actually happened)</option>` +
      c.models.filter((m) => m.deployable).map((m) =>
        `<option value="${esc(m.key)}">Inner-CV score of: ${esc(m.label)}` +
        `${m.selected_for_production ? "  ★ production" : ""}</option>`).join("");
    state.rpModel = "inner";
    sel.addEventListener("change", () => { state.rpModel = sel.value; drawReplay(); });
    $("#rp-fold").addEventListener("change", () => { state.rpFold = +$("#rp-fold").value; drawReplay(); });
    $("#rp-play").addEventListener("click", togglePlay);
    $("#rp-step").addEventListener("click", stepFold);
    drawReplay();
  }
  function stepFold() {
    state.rpFold = (state.rpFold + 1) % state.cal.meta.n_outer_folds;
    $("#rp-fold").value = state.rpFold; drawReplay();
  }
  function togglePlay() {
    if (state.playing) { stopPlay(); return; }
    state.playing = true; $("#rp-play").textContent = "❚❚ Pause";
    state.timer = setInterval(stepFold, 2400);
  }
  function stopPlay() {
    if (state.timer) clearInterval(state.timer);
    state.timer = 0; state.playing = false;
    const b = $("#rp-play"); if (b) b.textContent = "▶ Play";
  }

  function drawReplay() {
    const c = state.cal, k = state.rpFold;
    const pf = c.nested_evaluation.per_outer_fold[k];
    const assign = c.folds.outer_assignment;

    $("#rp-note").innerHTML =
      `<b>Inner CV</b> (${c.meta.n_inner_folds}-fold, grouped) ran inside the ` +
      `${pf.n_train_embryos} training embryos and chose <b>${esc(pf.chosen_by_inner_cv)}</b>. ` +
      `That model was then fitted on those ${pf.n_train_embryos} embryos and scored once on the ` +
      `<b>${pf.n_test_embryos} untouched</b> outer-test embryos ` +
      `(${pf.n_test_frames.toLocaleString()} frames). ` +
      `<span class="pt-note">Replay of precomputed offline artifacts — nothing is being fitted in ` +
      `your browser.</span>`;

    $("#rp-strip").innerHTML = Object.keys(assign).sort().map((e) => {
      const h = assign[e] === k;
      return `<span class="pt-fs ${h ? "pt-fs-held" : "pt-fs-train"}" title="${esc(e)} · outer fold ` +
        `${assign[e] + 1}${h ? " · OUTER TEST (untouched)" : " · outer training (inner CV ran here)"}"></span>`;
    }).join("");

    const fr = state.frames.filter((f) => f.outer_fold === k);
    $("#rp-fit-sub").textContent = `· outer fold ${k + 1} · ${pf.n_test_embryos} untouched embryos`;
    plot("rp-fit", [
      { type: "scatter", mode: "lines", x: [0, 1], y: [0, 1],
        line: { color: "#94a3b8", width: 1.4, dash: "dash" }, hoverinfo: "skip", showlegend: false },
      { type: "scatter", mode: "markers", x: fr.map((f) => f.tau_true), y: fr.map((f) => f.tau_pred),
        marker: { size: 5, color: FOLD_C[k % 5], opacity: 0.62, line: { width: 0 } },
        text: fr.map((f) => f.embryo_id),
        hovertemplate: "%{text}<br>true τ %{x:.3f}<br>predicted τ %{y:.3f}<extra></extra>" },
    ], { xaxis: { title: "true τ (untouched outer-test embryos)", range: [-0.03, 1.03], gridcolor: C.grid, zeroline: false },
         yaxis: { title: "predicted τ", range: [-0.03, 1.03], gridcolor: C.grid, zeroline: false } });

    // second panel: either the outer-test residuals, or the inner-CV scores of one candidate
    if (state.rpModel === "inner") {
      plot("rp-res", [
        { type: "scatter", mode: "lines", x: [0, 1], y: [0, 0], line: { color: "#94a3b8", width: 1.4 }, hoverinfo: "skip" },
        { type: "scatter", mode: "markers", x: fr.map((f) => f.tau_true), y: fr.map((f) => f.residual),
          marker: { size: 5, color: FOLD_C[k % 5], opacity: 0.62 }, text: fr.map((f) => f.embryo_id),
          hovertemplate: "%{text}<br>true τ %{x:.3f}<br>residual %{y:+.3f}<extra></extra>" },
      ], { xaxis: { title: "true τ", range: [-0.03, 1.03], gridcolor: C.grid, zeroline: false },
           yaxis: { title: "outer-test residual (predicted − true)", gridcolor: C.grid, zeroline: false } });
    } else {
      const keys = Object.keys(pf.inner_scores);
      const vals = keys.map((kk) => pf.inner_scores[kk].macro_mae);
      plot("rp-res", [{ type: "bar", orientation: "h", x: vals.slice().reverse(),
        y: keys.slice().reverse(),
        marker: { color: keys.slice().reverse().map((kk) =>
          kk === pf.chosen_by_inner_cv ? C.train : (kk === state.rpModel ? "#0d9488" : "#c7cddb")) },
        text: vals.slice().reverse().map((v) => v.toFixed(4)), textposition: "outside",
        textfont: { size: 9.5 },
        hovertemplate: "%{y}<br>inner-CV macro MAE %{x:.4f}<extra></extra>" }],
        { margin: { l: 130, r: 46, t: 8, b: 40 },
          xaxis: { title: "inner-CV macro MAE (selection only)", gridcolor: C.grid, zeroline: false,
                   range: [0, Math.max(...vals) * 1.25] },
          yaxis: { automargin: true, tickfont: { size: 9.5 } } });
    }

    const rows = c.nested_evaluation.per_outer_fold.map((f) => {
      const m = f.outer_test_metrics;
      return `<tr class="${f.outer_fold === k ? "pt-row-sel" : ""}"><td>${f.outer_fold + 1}</td>` +
        `<td class="pt-mono" style="font-size:11px">${esc(f.chosen_by_inner_cv)}</td>` +
        `<td>${f.n_train_embryos}</td><td>${f.n_test_embryos}</td>` +
        `<td>${f.n_test_frames.toLocaleString()}</td><td>${fmt(m.macro_mae)}</td>` +
        `<td>${fmt(m.pooled_rmse)}</td><td>${fmt(m.pooled_spearman)}</td>` +
        `<td>${pct(ORD(m))}</td></tr>`;
    }).join("");
    const nm = c.nested_evaluation.outer_test_metrics;
    $("#rp-metrics").innerHTML =
      `<tr><th>Outer fold</th><th>Inner CV chose</th><th>Train emb.</th><th>Untouched emb.</th>` +
      `<th>Test frames</th><th>macro MAE</th><th>RMSE</th><th>Spearman ρ</th><th>Strict ordering</th></tr>` +
      rows +
      `<tr class="pt-row-sel" style="border-top:2px solid var(--line)"><td><b>All outer-test</b></td>` +
      `<td class="pt-note">pooled</td><td>—</td><td>53</td>` +
      `<td>${c.dataset.n_frames.toLocaleString()}</td><td><b>${fmt(nm.macro_mae)}</b></td>` +
      `<td>${fmt(nm.pooled_rmse)}</td><td>${fmt(nm.pooled_spearman)}</td>` +
      `<td>${pct(ORD(nm))}</td></tr>`;
  }

  // ══════════════════════════════ 3 · predicted vs true ══════════════════════════════
  function renderOof() {
    ["oof-color", "oof-band", "oof-ident"].forEach((id) =>
      $("#" + id).addEventListener("change", drawOof));
    $("#oof-csv").addEventListener("click", () => {
      const h = "embryo_id,outer_fold,time_h,tau_true,tau_pred,residual,lo95,hi95,covered,out_of_range,nearer_um,farther_um";
      const lines = state.frames.map((f) =>
        [f.embryo_id, f.fold, f.time_h, f.tau_true, f.tau_pred, f.residual, f.lo95, f.hi95,
         f.covered, f.out_of_range, f.nearer_um, f.farther_um].join(","));
      downloadText([h].concat(lines).join("\n"),
        `pseudotime_oof_frames_${state.cal.meta.model_version}.csv`);
    });
    drawOof();
  }
  function drawOof() {
    const fr = state.frames, u = state.cal.uncertainty, hw = u.halfwidth_mean;
    const mode = $("#oof-color").value;
    const traces = [];
    if ($("#oof-band").checked) {
      traces.push(
        { type: "scatter", mode: "lines", x: [0, 1], y: [hw, 1 + hw], line: { width: 0 },
          hoverinfo: "skip", showlegend: false },
        { type: "scatter", mode: "lines", x: [0, 1], y: [-hw, 1 - hw], line: { width: 0 },
          fill: "tonexty", fillcolor: "rgba(79,70,229,0.10)", hoverinfo: "skip",
          showlegend: false, name: "95% interval" });
    }
    if ($("#oof-ident").checked) {
      traces.push({ type: "scatter", mode: "lines", x: [0, 1], y: [0, 1],
        line: { color: "#64748b", width: 1.5, dash: "dash" }, hoverinfo: "skip" });
    }
    const hover = "%{text}<extra></extra>";
    const txt = fr.map((f) => `${f.embryo_id} · outer fold ${f.outer_fold + 1}<br>true τ ${f.tau_true.toFixed(3)}` +
      `<br>predicted τ ${f.tau_pred.toFixed(3)}<br>residual ${f.residual >= 0 ? "+" : ""}${f.residual.toFixed(3)}` +
      `<br>${f.covered ? "within" : "OUTSIDE"} 95% interval${f.out_of_range ? "<br>clipped to [0,1]" : ""}`);
    let marker;
    if (mode === "fold") {
      marker = { size: 5, opacity: 0.6, color: fr.map((f) => FOLD_C[f.outer_fold % 5]) };
    } else if (mode === "residual") {
      marker = { size: 5, opacity: 0.75, color: fr.map((f) => f.residual), colorscale: "RdBu",
        cmin: -0.4, cmax: 0.4, reversescale: true,
        colorbar: { title: { text: "residual", font: { size: 10 } }, thickness: 10, len: 0.7 } };
    } else if (mode === "embryo") {
      const ids = [...new Set(fr.map((f) => f.embryo_id))].sort();
      const idx = new Map(ids.map((e, i) => [e, i]));
      marker = { size: 5, opacity: 0.6, color: fr.map((f) => idx.get(f.embryo_id)),
        colorscale: "Turbo", cmin: 0, cmax: ids.length - 1 };
    } else {
      marker = { size: 5, opacity: 0.55, color: C.train };
    }
    traces.push({ type: "scatter", mode: "markers", x: fr.map((f) => f.tau_true),
      y: fr.map((f) => f.tau_pred), marker, text: txt, hovertemplate: hover });
    plot("oof-plot", traces, {
      xaxis: { title: "true τ (ground truth from live imaging)", range: [-0.03, 1.03], gridcolor: C.grid, zeroline: false },
      yaxis: { title: "out-of-fold predicted τ", range: [-0.03, 1.03], gridcolor: C.grid, zeroline: false } });
    const m = state.cal.nested_evaluation.outer_test_metrics;
    $("#oof-note").innerHTML =
      `${fr.length.toLocaleString()} frames from ${state.cal.dataset.n_embryos} embryos, each predicted by a model ` +
      `that never saw that embryo. macro MAE <b>${fmt(m.macro_mae)}</b> · Spearman <b>${fmt(m.pooled_spearman)}</b> · ` +
      `strict ordering <b>${pct(ORD(m))}</b> · inversions <b>${pct(m.pooled_ordering.inversion_rate)}</b> · ` +
      `independent 95% coverage <b>${pct(u.coverage_frame_marginal.mean)}</b> (measured on embryos that did not set the interval).`;
  }

  // ═════════════════════════════════ 4 · trajectories ═════════════════════════════════
  function renderTraj() {
    const sel = $("#tj-embryo");
    sel.innerHTML = state.traj.map((t, i) =>
      `<option value="${i}">${esc(t.embryo_id)} · fold ${t.outer_fold + 1} · ${t.migration_duration_h} h</option>`).join("");
    ["tj-embryo", "tj-dist", "tj-marks"].forEach((id) => $("#" + id).addEventListener("change", drawTraj));
    drawTraj();
  }
  function drawTraj() {
    const t = state.traj[+$("#tj-embryo").value];
    const showD = $("#tj-dist").checked, showM = $("#tj-marks").checked;
    const hw = state.cal.uncertainty.halfwidth;
    const resid = t.tau_pred.map((p, i) => p - t.tau_true[i]);
    const traces = [
      { type: "scatter", mode: "lines", x: t.time_h, y: t.tau_pred.map((p) => Math.min(1, p + hw)),
        line: { width: 0 }, hoverinfo: "skip", showlegend: false },
      { type: "scatter", mode: "lines", name: "95% interval", x: t.time_h,
        y: t.tau_pred.map((p) => Math.max(0, p - hw)),
        line: { width: 0 }, fill: "tonexty", fillcolor: "rgba(79,70,229,0.10)", hoverinfo: "skip" },
      { type: "scatter", mode: "lines+markers", name: "true τ", x: t.time_h, y: t.tau_true,
        line: { color: "#334155", width: 2 }, marker: { size: 4 },
        hovertemplate: "t = %{x:.2f} h<br>true τ %{y:.3f}<extra></extra>" },
      { type: "scatter", mode: "lines+markers", name: "predicted τ (out-of-fold)", x: t.time_h, y: t.tau_pred,
        line: { color: C.train, width: 2 }, marker: { size: 4 },
        hovertemplate: "t = %{x:.2f} h<br>predicted τ %{y:.3f}<extra></extra>" },
    ];
    if (showM) {
      const big = t.time_h.map((_, i) => i).filter((i) => Math.abs(resid[i]) > hw);
      if (big.length) traces.push({ type: "scatter", mode: "markers", name: "|residual| > 95% half-width",
        x: big.map((i) => t.time_h[i]), y: big.map((i) => t.tau_pred[i]),
        marker: { size: 9, symbol: "circle-open", color: C.bad, line: { width: 1.8 } },
        hovertemplate: "t = %{x:.2f} h<br>outside the 95% interval<extra></extra>" });
      const inv = [];
      for (let i = 1; i < t.tau_pred.length; i++) if (t.tau_pred[i] < t.tau_pred[i - 1] - 1e-9) inv.push(i);
      if (inv.length) traces.push({ type: "scatter", mode: "markers", name: "ordering inversion",
        x: inv.map((i) => t.time_h[i]), y: inv.map((i) => t.tau_pred[i]),
        marker: { size: 7, symbol: "x", color: C.warn },
        hovertemplate: "t = %{x:.2f} h<br>predicted τ decreased vs the previous frame<extra></extra>" });
    }
    const layout = {
      showlegend: true,
      legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", font: { size: 10 } },
      margin: { l: 58, r: showD ? 58 : 14, t: 54, b: 46 },
      xaxis: { title: "hours since pronuclear formation", gridcolor: C.grid, zeroline: false },
      yaxis: { title: "τ", range: [-0.05, 1.05], gridcolor: C.grid, zeroline: false },
    };
    if (showD) {
      traces.push(
        { type: "scatter", mode: "lines", name: "nearer → centre (µm)", x: t.time_h, y: t.nearer_um,
          yaxis: "y2", line: { color: "#0d9488", width: 1.4, dash: "dot" },
          hovertemplate: "t = %{x:.2f} h<br>nearer %{y:.1f} µm<extra></extra>" },
        { type: "scatter", mode: "lines", name: "farther → centre (µm)", x: t.time_h, y: t.farther_um,
          yaxis: "y2", line: { color: "#b45309", width: 1.4, dash: "dot" },
          hovertemplate: "t = %{x:.2f} h<br>farther %{y:.1f} µm<extra></extra>" });
      layout.yaxis2 = { title: "distance to cell centre (µm)", overlaying: "y", side: "right",
        showgrid: false, zeroline: false, rangemode: "tozero" };
    }
    plot("tj-plot", traces, layout);

    const mae = resid.reduce((s, r) => s + Math.abs(r), 0) / resid.length;
    const cov = resid.filter((r) => Math.abs(r) <= hw).length / resid.length;
    let inv = 0;
    for (let i = 1; i < t.tau_pred.length; i++) if (t.tau_pred[i] < t.tau_pred[i - 1] - 1e-9) inv++;
    const per = (state.cal.nested_evaluation.per_embryo || {})[t.embryo_id] || {};
    const kpi = (n, k) => `<div class="pt-kpi"><div class="pt-kpi-n">${n}</div><div class="pt-kpi-k">${k}</div></div>`;
    $("#tj-kpis").innerHTML =
      kpi(t.n, "frames") + kpi(t.migration_duration_h + " h", "migration duration") +
      kpi(fmt(mae), "MAE for this embryo") +
      kpi(per.spearman == null ? "—" : fmt(per.spearman), "within-embryo Spearman") +
      kpi(pct(cov), "within its 95% interval") +
      kpi(inv, "predicted-τ inversions<br>between consecutive frames");
  }

  // ═══════════════════════════════ 5 · model comparison ═══════════════════════════════
  // Ranked by INNER-CV macro MAE — the quantity selection actually used. The outer-test number is
  // NOT shown per candidate, because only the inner-CV choice per fold was ever evaluated there;
  // quoting an outer-test score for every candidate would reintroduce exactly the selection
  // optimism this version removes.
  function renderModels() {
    const c = state.cal;
    const dep = c.models.filter((m) => m.deployable)
      .slice().sort((a, b) => a.inner_cv_macro_mae_mean - b.inner_cv_macro_mae_mean);
    const leaky = c.models.filter((m) => !m.deployable);

    plot("cmp-plot", [{
      type: "bar", orientation: "h",
      x: dep.map((m) => m.inner_cv_macro_mae_mean).reverse(),
      y: dep.map((m) => m.label).reverse(),
      marker: { color: dep.map((m) => (m.selected_for_production ? C.train : "#c7cddb")).reverse(),
                line: { color: "rgba(15,23,42,.25)", width: 0.7 } },
      text: dep.map((m) => fmt(m.inner_cv_macro_mae_mean, 4) +
        (m.selected_for_production ? "  ★" : "")).reverse(),
      textposition: "outside", textfont: { size: 10.5 },
      hovertemplate: "%{y}<br>inner-CV macro MAE %{x:.4f}<extra></extra>",
    }], {
      margin: { l: 250, r: 60, t: 8, b: 44 }, showlegend: false,
      xaxis: { title: "mean inner-CV macro MAE across outer folds (selection criterion)",
               gridcolor: C.grid, zeroline: false,
               range: [0, Math.max(...dep.map((m) => m.inner_cv_macro_mae_mean)) * 1.22] },
      yaxis: { automargin: true, tickfont: { size: 10.5 } },
    });

    const nOuter = c.meta.n_outer_folds;
    $("#cmp-table").innerHTML =
      `<tr><th>Model</th><th>Features</th><th>Complexity</th><th>Inner-CV macro MAE</th>` +
      `<th>Inner-CV strict ordering</th><th>Chosen in</th><th>Status</th></tr>` +
      dep.map((m) =>
        `<tr class="${m.selected_for_production ? "pt-row-sel" : ""}">` +
        `<td><b>${esc(m.label)}</b><br><span class="pt-note">${esc(m.note)}</span></td>` +
        `<td>${m.n_features}</td><td class="pt-wrap pt-note">${esc(m.complexity)}</td>` +
        `<td><b>${fmt(m.inner_cv_macro_mae_mean, 4)}</b></td>` +
        `<td>${pct(m.inner_cv_strict_ordering_mean)}</td>` +
        `<td>${m.chosen_in_n_outer_folds} / ${nOuter} outer folds</td>` +
        `<td>${m.selected_for_production ? '<span class="pt-pill pt-pill-sel">PRODUCTION</span>'
             : '<span class="pt-pill pt-pill-yes">deployable</span>'}</td></tr>`).join("");

    const stab = c.nested_evaluation.inner_selection_stability || {};
    const stabTxt = Object.entries(stab).sort((a, b) => b[1] - a[1])
      .map(([k, v]) => `<b>${esc(k)}</b> in ${v}/${nOuter}`).join(", ");
    $("#cmp-note").innerHTML =
      `Ranked by the criterion selection actually used: <b>mean inner-CV macro MAE</b>. ` +
      `The unbiased performance number is the single nested outer-test estimate ` +
      `(<b>${fmt(c.nested_evaluation.outer_test_metrics.macro_mae)}</b>), not any per-candidate ` +
      `score here. <br><b>Selection stability:</b> the inner CV did not always pick the same ` +
      `family — ${stabTxt}. That instability is real and is reported rather than hidden; it means ` +
      `the top few candidates are close enough that the choice is not strongly determined.`;

    $("#cmp-leaky").innerHTML =
      `<tr><th>Model</th><th>Features</th><th>Nested outer-test macro MAE</th>` +
      `<th>Spearman ρ</th><th>Strict ordering</th><th>Status</th></tr>` +
      leaky.map((m) => {
        const mm = m.nested_outer_test_metrics || {};
        return `<tr><td><b>${esc(m.label)}</b><br><span class="pt-note">${esc(m.note)}</span></td>` +
          `<td>${m.n_features}</td><td><b>${fmt(mm.macro_mae)}</b></td>` +
          `<td>${fmt(mm.pooled_spearman)}</td><td>${pct(ORD(mm))}</td>` +
          `<td><span class="pt-pill pt-pill-no">non-deployable</span></td></tr>`;
      }).join("");
    $("#cmp-leaky-note").innerHTML =
      `These use pronuclear volumes normalised to each pronucleus's <i>own endpoint volume</i>, which
       (a) cannot be measured in a fixed snapshot and (b) directly encode elapsed time. They are
       excluded from selection in code (<span class="pt-mono">deployable: false</span>), never appear
       in the production selector, and must not be quoted as achievable accuracy. They are shown only
       as an upper bound on what the geometry could do if the future were known.`;
  }

  // ══════════════════════════════════ 6 · diagnostics ══════════════════════════════════
  function renderDiag() {
    const fr = state.frames, u = state.cal.uncertainty, hw = u.halfwidth_mean;
    plot("dg-scatter", [
      { type: "scatter", mode: "lines", x: [0, 1], y: [0, 0], line: { color: "#64748b", width: 1.4 }, hoverinfo: "skip" },
      { type: "scatter", mode: "lines", x: [0, 1], y: [hw, hw], line: { color: C.train, width: 1.2, dash: "dash" }, hoverinfo: "skip" },
      { type: "scatter", mode: "lines", x: [0, 1], y: [-hw, -hw], line: { color: C.train, width: 1.2, dash: "dash" }, hoverinfo: "skip" },
      { type: "scatter", mode: "markers", x: fr.map((f) => f.tau_true), y: fr.map((f) => f.residual),
        marker: { size: 4.5, opacity: 0.5, color: fr.map((f) => (f.covered ? C.train : C.bad)) },
        text: fr.map((f) => f.embryo_id),
        hovertemplate: "%{text}<br>true τ %{x:.3f}<br>residual %{y:+.3f}<extra></extra>" },
    ], { xaxis: { title: "true τ", range: [-0.03, 1.03], gridcolor: C.grid, zeroline: false },
         yaxis: { title: "residual (predicted − true)", gridcolor: C.grid, zeroline: false } });

    plot("dg-hist", [{ type: "histogram", x: fr.map((f) => f.residual), nbinsx: 46,
      marker: { color: C.train, opacity: 0.82, line: { color: "#fff", width: 0.5 } },
      hovertemplate: "residual %{x}<br>%{y} frames<extra></extra>" }],
      { xaxis: { title: "residual (predicted − true τ)", gridcolor: C.grid, zeroline: false },
        yaxis: { title: "frames", gridcolor: C.grid, zeroline: false },
        shapes: [{ type: "line", x0: 0, x1: 0, yref: "paper", y0: 0, y1: 1,
                   line: { color: "#64748b", width: 1.4, dash: "dash" } }] });

    const ph = Object.entries(u.coverage_by_phase || {}).map(([k, v]) => ({ phase: k, ...v }));
    const rd = state.cal.residual_diagnostics.by_phase || [];
    plot("dg-cov", [
      { type: "bar", x: ph.map((p) => p.phase), y: ph.map((p) => p.mean),
        error_y: { type: "data", symmetric: false,
                   array: ph.map((p) => p.ci95[1] - p.mean),
                   arrayminus: ph.map((p) => p.mean - p.ci95[0]),
                   color: "#334155", thickness: 1.1, width: 5 },
        marker: { color: ph.map((p) => (Math.abs(p.mean - 0.95) <= 0.05 ? C.ok : C.warn)) },
        text: ph.map((p) => pct(p.mean)), textposition: "outside", textfont: { size: 10.5 },
        hovertemplate: "%{x}<br>coverage %{y:.3f}<extra></extra>" },
    ], { margin: { l: 52, r: 14, t: 8, b: 40 },
         xaxis: { title: "migration phase (by true τ)", gridcolor: C.grid },
         yaxis: { title: "95% coverage", range: [0, 1.12], gridcolor: C.grid, zeroline: false },
         shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: 0.95, y1: 0.95,
                    line: { color: "#64748b", width: 1.4, dash: "dash" } }] });

    $("#dg-cov-table").innerHTML =
      `<tr><th>Phase</th><th>Frames</th><th>Outer-test MAE</th><th>Bias</th>` +
      `<th>95% coverage (independent embryos)</th><th>Coverage 95% CI</th></tr>` +
      ph.map((p) => {
        const r = rd.find((x) => x.phase === p.phase) || {};
        return `<tr><td>${esc(p.phase)}</td><td>${(r.n || 0).toLocaleString()}</td>` +
          `<td>${fmt(r.mae)}</td><td>${r.bias >= 0 ? "+" : ""}${fmt(r.bias)}</td>` +
          `<td>${pct(p.mean)}</td><td class="pt-note">${pct(p.ci95[0])} – ${pct(p.ci95[1])}</td></tr>`;
      }).join("") +
      `<tr class="pt-row-sel"><td><b>Overall</b></td><td>${fr.length.toLocaleString()}</td>` +
      `<td>${fmt(state.cal.residual_diagnostics.overall_mae)}</td>` +
      `<td>${state.cal.residual_diagnostics.overall_bias >= 0 ? "+" : ""}${fmt(state.cal.residual_diagnostics.overall_bias)}</td>` +
      `<td>${pct(u.coverage_frame_marginal.mean)}</td>` +
      `<td class="pt-note">${pct(u.coverage_frame_marginal.ci95[0])} – ${pct(u.coverage_frame_marginal.ci95[1])}</td></tr>`;

    const cv = u.conservative_variant || {};
    $("#dg-unc").innerHTML =
      `<p class="pt-lede" style="margin-bottom:10px"><b>Method.</b> ${esc(u.method)}</p>` +
      `<div class="pt-two" style="margin-bottom:12px">
         <div class="pt-block"><h2>Published interval</h2>
           <div class="pt-kpis pt-kpis-sm">
             <div class="pt-kpi"><div class="pt-kpi-n">±${fmt(u.halfwidth_mean)}</div>
               <div class="pt-kpi-k">half-width in τ</div></div>
             <div class="pt-kpi"><div class="pt-kpi-n">${pct(u.coverage_frame_marginal.mean)}</div>
               <div class="pt-kpi-k">frame-marginal coverage</div>
               <div class="pt-kpi-ci">95% CI ${pct(u.coverage_frame_marginal.ci95[0])} – ${pct(u.coverage_frame_marginal.ci95[1])}</div></div>
             <div class="pt-kpi"><div class="pt-kpi-n">${pct(u.coverage_embryo_macro.mean)}</div>
               <div class="pt-kpi-k">embryo-macro coverage</div>
               <div class="pt-kpi-ci">95% CI ${pct(u.coverage_embryo_macro.ci95[0])} – ${pct(u.coverage_embryo_macro.ci95[1])}</div></div>
           </div></div>
         <div class="pt-block"><h2>Conservative variant</h2>
           <p class="pt-note">${esc(cv.definition || "")}</p>
           <div class="pt-kpis pt-kpis-sm">
             <div class="pt-kpi"><div class="pt-kpi-n">±${fmt(cv.halfwidth_mean)}</div>
               <div class="pt-kpi-k">half-width in τ</div></div>
             <div class="pt-kpi"><div class="pt-kpi-n">${pct((cv.coverage_frame_marginal || {}).mean)}</div>
               <div class="pt-kpi-k">frame-marginal coverage</div></div>
           </div>
           <p class="pt-note">${esc(cv.comment || "")}</p></div>
       </div>` +
      `<div class="pt-lim"><b>Interval type</b><span>${esc(u.interval_type)}</span></div>` +
      `<div class="pt-lim"><b>Small-sample honesty</b><span>${esc(u.embryo_exchangeable_note)}</span></div>` +
      u.limitations.map((l) => `<div class="pt-lim"><span>${esc(l)}</span></div>`).join("");
  }

  // ═════════════════════════════════ 7 · fixed cohort ═════════════════════════════════
  function renderFixed() {
    if (!state.fixed) {
      $("#sec-fixed").innerHTML =
        `<h1>Fixed MERFISH cohort</h1><div class="pt-empty">
         <b>Not built.</b><br>Run <code>python3 build_pronuclei_pseudotime.py</code>.</div>`;
      return;
    }
    ["fx-order", "fx-ood"].forEach((id) => $("#" + id).addEventListener("change", drawFixed));
    $("#fx-csv").addEventListener("click", () => {
      const h = "embryo_id,label,tau,lo95,hi95,qc,reason,legacy_surface_gap_um,nearer_to_center_um," +
        "farther_to_center_um,distance_sum_um,distance_difference_um,mahalanobis,model_version,data_version,feature_schema";
      const lines = state.fixed.embryos.map((r) => {
        const f = r.features || {};
        return [r.id, r.label, r.tau ?? "", r.lo95 ?? "", r.hi95 ?? "", r.qc,
          `"${(r.reason || "").replace(/"/g, "'")}"`, r.legacy_surface_gap_um,
          f.nearer_to_center_um ?? "", f.farther_to_center_um ?? "", f.distance_sum_um ?? "",
          f.distance_difference_um ?? "", r.mahalanobis ?? "", r.model_version, r.data_version,
          `"${(r.feature_schema || []).join(" ")}"`].join(",");
      });
      downloadText([h].concat(lines).join("\n"),
        `fixed_cohort_pseudotime_${state.fixed.meta.model_version}.csv`);
    });
    $("#fx-json").addEventListener("click", () =>
      downloadText(JSON.stringify(state.fixed, null, 1),
        `fixed_cohort_pseudotime_${state.fixed.meta.model_version}.json`, "application/json"));
    drawFixed();
  }

  const QC_C = { pass: C.ok, caution: C.warn, "out-of-domain": C.bad };

  function fixedRows() {
    let rows = state.fixed.embryos.slice();
    // Out-of-domain zygotes are EXCLUDED BY DEFAULT from the displayed cohort and from anything
    // downstream; the toggle reveals them explicitly rather than hiding them silently.
    if (!$("#fx-ood").checked) rows = rows.filter((r) => r.qc !== "out-of-domain");
    const by = $("#fx-order").value;
    rows.sort((a, b) => {
      if (by === "label") return String(a.label).localeCompare(String(b.label), undefined, { numeric: true });
      if (by === "legacy") return (a.legacy_surface_gap_um ?? 1e9) - (b.legacy_surface_gap_um ?? 1e9);
      return (a.tau ?? 1e9) - (b.tau ?? 1e9);
    });
    return rows;
  }

  function drawFixed() {
    const rows = fixedRows(), meta = state.fixed.meta;
    const ok = rows.filter((r) => r.tau != null);
    const nOod = state.fixed.embryos.filter((r) => r.qc === "out-of-domain").length;
    const y = ok.map((r) => r.label);
    plot("fx-plot", [{
      type: "scatter", mode: "markers", x: ok.map((r) => r.tau), y,
      error_x: { type: "data", symmetric: false,
        array: ok.map((r) => r.hi95 - r.tau), arrayminus: ok.map((r) => r.tau - r.lo95),
        color: "rgba(15,23,42,.32)", thickness: 1, width: 0 },
      marker: { size: 8, color: ok.map((r) => QC_C[r.qc] || C.apply),
                line: { color: "#fff", width: 1 } },
      text: ok.map((r) => `${r.label}<br>τ ${r.tau.toFixed(3)} [${r.lo95.toFixed(2)}, ${r.hi95.toFixed(2)}]` +
        `<br>QC: ${r.qc}${r.reason ? "<br>" + r.reason.replace(/; /g, "<br>") : ""}` +
        `<br>legacy gap ${r.legacy_surface_gap_um} µm`),
      hovertemplate: "%{text}<extra></extra>",
    }], {
      margin: { l: 108, r: 20, t: 8, b: 46 },
      height: Math.max(420, 15 * ok.length + 70),
      xaxis: { title: "calibrated τ (0 = pronuclear formation, 1 = NEBD)", range: [-0.04, 1.04],
               gridcolor: C.grid, zeroline: false },
      yaxis: { automargin: true, tickfont: { size: 9.5 }, type: "category" },
    });
    $("#fx-plot").style.height = Math.max(420, 15 * ok.length + 70) + "px";

    // τ vs legacy score — two different quantities, never merged
    const both = rows.filter((r) => r.tau != null && r.legacy_surface_gap_um != null);
    plot("fx-cmp", [{
      type: "scatter", mode: "markers", x: both.map((r) => r.legacy_surface_gap_um),
      y: both.map((r) => r.tau),
      marker: { size: 8, color: both.map((r) => QC_C[r.qc] || C.apply), line: { color: "#fff", width: 1 } },
      text: both.map((r) => `${r.label}<br>legacy gap ${r.legacy_surface_gap_um} µm<br>calibrated τ ${r.tau.toFixed(3)}<br>QC ${r.qc}`),
      hovertemplate: "%{text}<extra></extra>",
    }], {
      xaxis: { title: "legacy minimum surface-to-surface pronuclear gap (µm)", gridcolor: C.grid, zeroline: false },
      yaxis: { title: "calibrated τ", range: [-0.04, 1.04], gridcolor: C.grid, zeroline: false } });

    const ds = meta.domain_shift || {};
    $("#fx-legacy-note").innerHTML =
      `<b>These are different quantities and are never merged.</b> The legacy score is the minimum
       surface-to-surface gap between the two pronuclei; calibrated τ is a cell-centred distance clock.
       The Scheffler workbook contains <b>no surface-gap measurement</b>, so this calibration neither
       validates nor calibrates the legacy score — the panel is a comparison, not a conversion.` +
      (ds.consequence ? `<br><br><b>Domain shift.</b> ${esc(ds.consequence.note)}` : "");

    const qc = rows.reduce((m, r) => (m[r.qc] = (m[r.qc] || 0) + 1, m), {});
    $("#fx-table").innerHTML =
      `<tr><th>Embryo</th><th>τ</th><th>95% interval</th><th>QC</th><th>nearer µm</th><th>farther µm</th>` +
      `<th>Σ µm</th><th>Mahalanobis</th><th>Legacy gap µm</th><th>Model</th><th>Notes</th></tr>` +
      rows.map((r) => {
        const f = r.features || {};
        const cls = r.qc === "pass" ? "pt-qc-pass" : r.qc === "caution" ? "pt-qc-caution" : "pt-qc-out";
        return `<tr data-id="${esc(r.id)}" title="Open ${esc(r.label)} in the Pronuclei viewer">` +
          `<td><b>${esc(r.label)}</b></td>` +
          `<td>${r.tau == null ? '<span class="pt-note">not available</span>' : fmt(r.tau)}</td>` +
          `<td>${r.tau == null ? "—" : `${fmt(r.lo95, 2)} – ${fmt(r.hi95, 2)}`}</td>` +
          `<td><span class="pt-pill ${cls}">${esc(r.qc)}</span></td>` +
          `<td>${f.nearer_to_center_um == null ? "—" : fmt(f.nearer_to_center_um, 1)}</td>` +
          `<td>${f.farther_to_center_um == null ? "—" : fmt(f.farther_to_center_um, 1)}</td>` +
          `<td>${f.distance_sum_um == null ? "—" : fmt(f.distance_sum_um, 1)}</td>` +
          `<td>${r.mahalanobis == null ? "—" : fmt(r.mahalanobis, 2)}</td>` +
          `<td>${r.legacy_surface_gap_um ?? "—"}</td>` +
          `<td class="pt-mono" style="font-size:10.5px">${esc(r.model_version)}</td>` +
          `<td class="pt-wrap pt-note">${esc(r.reason || "")}</td></tr>`;
      }).join("");
    $("#fx-table").onclick = (e) => {
      const tr = e.target.closest("tr[data-id]");
      if (tr) window.location.href = `pronuclei.html?embryo=${encodeURIComponent(tr.dataset.id)}`;
    };

    const kpiWrap = $("#sec-fixed").querySelector(".pt-kpis") || (() => {
      const d = document.createElement("div"); d.className = "pt-kpis pt-kpis-sm";
      $("#fx-plot").parentNode.insertBefore(d, $("#fx-plot")); return d;
    })();
    const kpi = (n, k, cls) => `<div class="pt-kpi ${cls || ""}"><div class="pt-kpi-n">${n}</div><div class="pt-kpi-k">${k}</div></div>`;
    kpiWrap.innerHTML =
      kpi(meta.n_predicted + " / " + meta.n_total, "zygotes with a τ estimate") +
      kpi(qc.pass || 0, "QC pass", "pt-kpi-ok") +
      kpi(qc.caution || 0, "QC caution", (qc.caution ? "pt-kpi-warn" : "")) +
      kpi(qc["out-of-domain"] || 0, "out-of-domain", (qc["out-of-domain"] ? "pt-kpi-warn" : "")) +
      kpi("±" + fmt(meta.halfwidth_95), "95% half-width (from live validation)");
    const ex = $("#fx-exclude");
    if (ex) ex.innerHTML = $("#fx-ood").checked
      ? `<b>Showing all ${state.fixed.embryos.length} zygotes, including ${nOod} out-of-domain.</b> ` +
        `Out-of-domain geometry lies outside the training range, so those τ values are ` +
        `extrapolations and are excluded from every downstream regression regardless of this toggle.`
      : `<b>${nOod} out-of-domain ${nOod === 1 ? "zygote is" : "zygotes are"} excluded by default.</b> ` +
        `Their geometry lies outside the training range, so their τ would be an extrapolation. ` +
        `Tick the toggle above to display them.`;
  }

  // ═════════════════════════ 8 · R² attenuation / noise ceiling ═════════════════════════
  // The analysis Matt asked for: use the measured clock error to say what transcript-vs-pseudotime
  // R² should be EXPECTED, so an observed value can be interpreted instead of dismissed.
  function renderNoise() {
    const nc = state.noise;
    if (!nc) {
      $("#sec-noise").innerHTML =
        `<h1>What pseudotime noise does to transcript R²</h1><div class="pt-empty">
         <b>Not built.</b><br>Run <code>python3 scripts/analyze_pseudotime_noise_ceiling.py</code>.</div>`;
      return;
    }
    const rel = nc.clock_reliability, att = nc.attenuation, dn = nc.downstream_illustration || {};

    $("#nc-guard").innerHTML =
      `<b>What this does and does not show.</b> ${esc(att.interpretation_guard)}`;

    const ci = rel.bootstrap_ci95_by_embryo || {};
    const kpi = (n, k, sub) => `<div class="pt-kpi"><div class="pt-kpi-n">${n}</div>` +
      `<div class="pt-kpi-k">${k}</div>${sub ? `<div class="pt-kpi-ci">${sub}</div>` : ""}</div>`;
    const cit = (a, f) => (a ? `95% CI ${(f || fmt)(a[0])} – ${(f || fmt)(a[1])}` : "");
    $("#nc-rel").innerHTML =
      kpi(fmt(rel.pearson_r2), "R² of predicted vs true τ", cit(ci.pearson_r2)) +
      kpi(fmt(rel.pearson_r), "Pearson r", cit(ci.pearson_r)) +
      kpi(fmt(rel.spearman), "Spearman ρ", cit(ci.spearman)) +
      kpi(fmt(rel.macro_mae_tau), "macro MAE (τ)", cit(ci.macro_mae)) +
      kpi(fmt(rel.median_ae_tau), "median |error| (τ)", cit(ci.median_ae)) +
      kpi(fmt(rel.median_ae_hours, 2) + " h", "median |error| (hours)", cit(ci.median_ae_hours, (v) => fmt(v, 2)));
    $("#nc-hours").innerHTML =
      `<b>Hours caveat.</b> ${esc(rel.hours_caveat)} Ordering on the same untouched predictions: ` +
      `strict pairwise ${pct(rel.strict_pairwise_ordering)}, within-embryo ` +
      `${pct(rel.macro_within_embryo_ordering)}, inversions ${pct(rel.inversion_rate)}, ` +
      `prediction ties ${pct(rel.tie_rate)}.`;

    const shapes = Object.keys(att.curves);
    $("#nc-shape").innerHTML = shapes.map((k) =>
      `<option value="${esc(k)}">${esc(k.replace(/_/g, " "))}</option>`).join("");
    ["nc-shape", "nc-band", "nc-ident", "nc-marker"].forEach((id) =>
      $("#" + id).addEventListener("change", drawNoise));
    drawNoise();

    // ---- predeclared downstream illustration ----
    const h = dn.headline || {};
    const vcls = h.verdict === "not-significant" ? "pt-v-bad"
      : h.verdict === "significant-but-outlier-driven" ? "pt-v-warn" : "pt-v-ok";
    $("#nc-verdict").innerHTML = dn.results
      ? `<div class="pt-verdict-box ${vcls}"><div class="pt-verdict-k">` +
        `${esc(dn.confirmatory_gene)} · verdict: ${esc(h.verdict || "—")}</div>` +
        `<div class="pt-verdict-t">${esc(h.text || "")}</div>` +
        (h.answer_to_the_low_r2_question
          ? `<div class="pt-verdict-a"><b>Does this rescue an R² of 0.1–0.15?</b> ` +
            `${esc(h.answer_to_the_low_r2_question)}</div>` : "") +
        `</div>`
      : `<div class="pt-empty">${esc(dn.reason || dn.status || "not estimable")}</div>`;

    if (dn.results) {
      $("#nc-genes").innerHTML =
        `<tr><th>Gene</th><th>n</th><th>Observed R²</th><th>R² 95% CI</th><th>Spearman ρ</th>` +
        `<th>log1p R²</th><th>Permutation p</th><th>Direction</th><th>Robustness</th></tr>` +
        Object.entries(dn.results).map(([g, r]) => {
          if (r.status !== "estimable") {
            return `<tr><td><b>${esc(g)}</b></td><td>${r.n ?? "—"}</td>` +
              `<td colspan="7" class="pt-wrap pt-note">not estimable — ${esc(r.reason || "")}</td></tr>`;
          }
          const rb = r.robustness || {};
          const bad = rb.verdict && rb.verdict.startsWith("OUTLIER");
          return `<tr class="${g === dn.confirmatory_gene ? "pt-row-sel" : ""}">` +
            `<td><b>${esc(g)}</b>${g === dn.confirmatory_gene ? ' <span class="pt-pill pt-pill-sel">confirmatory</span>' : ""}</td>` +
            `<td>${r.n}</td><td><b>${fmt(r.observed_r2)}</b></td>` +
            `<td class="pt-note">${r.r2_bootstrap_ci95 ? fmt(r.r2_bootstrap_ci95[0]) + " – " + fmt(r.r2_bootstrap_ci95[1]) : "—"}</td>` +
            `<td>${fmt(rb.spearman_rho)}</td><td>${fmt(rb.log1p_ols_r2)}</td>` +
            `<td>${fmt(r.permutation_p_embryo_level)}</td>` +
            `<td class="pt-note">${esc(r.direction || "")}</td>` +
            `<td class="pt-wrap"><span class="pt-pill ${bad ? "pt-pill-no" : "pt-pill-yes"}">` +
            `${bad ? "outlier-driven" : "stable"}</span></td></tr>`;
        }).join("");
      const sc = dn.discovery_scan || {};
      $("#nc-scan").innerHTML =
        `<b>Genes are predeclared.</b> ${esc(dn.predeclared_note || "")} ` +
        `A separate all-gene <b>discovery</b> scan (${sc.n_genes_tested || 0} genes, ` +
        `Benjamini–Hochberg) found <b>${sc.n_q_below_0_05 ?? 0}</b> with q &lt; 0.05. ${esc(sc.note || "")}`;
    }
  }

  function drawNoise() {
    const att = state.noise.attenuation;
    const dn = state.noise.downstream_illustration || {};
    const rows = att.curves[$("#nc-shape").value] || [];
    const x = rows.map((r) => r.latent_true_r2);
    const traces = [];
    if ($("#nc-band").checked) {
      traces.push(
        { type: "scatter", mode: "lines", x, y: rows.map((r) => r.observed_r2_ci95[1]),
          line: { width: 0 }, hoverinfo: "skip", showlegend: false },
        { type: "scatter", mode: "lines", name: "95% simulation band", x,
          y: rows.map((r) => r.observed_r2_ci95[0]), line: { width: 0 }, fill: "tonexty",
          fillcolor: "rgba(79,70,229,0.13)", hoverinfo: "skip" });
    }
    if ($("#nc-ident").checked) {
      traces.push({ type: "scatter", mode: "lines", name: "no attenuation (y = x)", x: [0, 1], y: [0, 1],
        line: { color: "#94a3b8", width: 1.5, dash: "dash" }, hoverinfo: "skip" });
    }
    traces.push({ type: "scatter", mode: "lines+markers", name: "median observed R²", x,
      y: rows.map((r) => r.observed_r2_median),
      line: { color: C.train, width: 2.4 }, marker: { size: 6 },
      hovertemplate: "latent true R² %{x:.2f}<br>observed median %{y:.3f}<extra></extra>" });

    const conf = (dn.results || {})[dn.confirmatory_gene];
    if ($("#nc-marker").checked && conf && conf.status === "estimable") {
      const o = conf.observed_r2;
      traces.push({ type: "scatter", mode: "lines", x: [0, 1], y: [o, o],
        line: { color: C.bad, width: 1.6, dash: "dot" },
        name: `${dn.confirmatory_gene} observed R² = ${fmt(o)}`,
        hovertemplate: `${esc(dn.confirmatory_gene)} observed R² ${fmt(o)}<extra></extra>` });
      const comp = conf.compatible_latent_true_r2;
      if (comp) {
        traces.push({ type: "scatter", mode: "markers",
          x: [comp.min, comp.max], y: [o, o],
          marker: { size: 11, color: C.bad, symbol: "diamond", line: { color: "#fff", width: 1.4 } },
          name: "compatible latent range",
          hovertemplate: "compatible latent true R² %{x:.2f}<extra></extra>" });
      }
    }
    plot("nc-plot", traces, {
      showlegend: true,
      legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", font: { size: 10 } },
      margin: { l: 60, r: 16, t: 54, b: 48 },
      xaxis: { title: "latent TRUE R² (biological trend vs real τ)", range: [0, 1], gridcolor: C.grid, zeroline: false },
      yaxis: { title: "observed R² (same trend vs PREDICTED τ)", range: [0, 1], gridcolor: C.grid, zeroline: false },
    });
    const d = att.design;
    $("#nc-note").innerHTML =
      `Simulated as a <b>snapshot cohort</b>: ${d.n_snapshot_embryos} embryos, one frame each, ` +
      `${d.n_replicates.toLocaleString()} replicates, ${esc(d.resampling)}. The error model is ` +
      `<b>${esc(d.error_model)}</b> — real (true, predicted) pairs are resampled together, so the ` +
      `step ties, heteroscedasticity and phase-dependent bias of the actual clock are preserved ` +
      `rather than replaced by a Gaussian approximation.`;
  }

  // ═══════════════════════════════════ limitations ═══════════════════════════════════
  function renderLimits() {
    const c = state.cal, u = c.uncertainty;
    const ds = (state.fixed && state.fixed.meta.domain_shift) || null;
    const nc = state.noise;
    const conf = nc && (nc.downstream_illustration.results || {})[nc.downstream_illustration.confirmatory_gene];
    const items = [
      ["This is an empirical calibration, not a mechanistic model",
       "The clock is a monotone map from pronuclear geometry to normalised time, fitted to data. It " +
       "makes no claim about the forces that move the pronuclei and should not be described as a " +
       "model of migration."],
      ["The calibration source cannot validate the legacy surface-gap score",
       "The public Scheffler workbook provides pronucleus→cell-centre distances only. It contains no " +
       "minimum surface-to-surface gap, so this work calibrates CELL-CENTRED DISTANCE features and " +
       "says nothing about the legacy score, which remains the site default and is unchanged."],
      ["Live→fixed transfer is assumed, not measured",
       "Training uses live-imaged embryos measured by Scheffler's pipeline; application uses fixed " +
       "MERFISH zygotes segmented by this project's pipeline. No raw 3-D stacks are public, so the " +
       "two geometries cannot be processed identically and the transfer error is unquantified. The " +
       "live-cohort error is therefore a LOWER BOUND on the fixed-cohort error."],
      ["Prediction intervals are wide, constant-width, and frame-marginal",
       `The 95% half-width is ±${fmt(u.halfwidth_mean)} in τ — about ±${(u.halfwidth_mean * 100).toFixed(0)}% ` +
       "of the whole pronuclear-formation→NEBD interval. It is marginal per frame, not simultaneous " +
       "per embryo, and constant in τ, so hard phases are under-covered. With 53 embryos a hard 95% " +
       "certificate is not attainable; the measured coverage and its bootstrap interval are reported " +
       "instead of asserting the nominal level."],
      ["Model selection is not strongly determined",
       "The inner CV did not choose the same family in every outer fold — the top candidates differ " +
       "by less than the fold-to-fold spread. The production lock is a documented tie-break, not a " +
       "clear winner, and a different seed could plausibly lock a different family."],
      ["One laboratory, one protocol, 53 trajectories",
       "There is no independent batch, imaging modality or laboratory to hold out, so generalisation " +
       "beyond this cohort is untested. This is a pilot-scale benchmark."],
      ["τ is normalised, not absolute time",
       "Migration duration varies 8.75–11.75 h. Errors are converted to hours only for the live " +
       "cohort, where each embryo's duration is observed; a fixed snapshot has no known duration, so " +
       "its τ cannot be expressed in hours."],
      ["Fixed zygotes have no pronuclear identity",
       "Only identity-free sorted near/far features are deployable. An identity-aware model might do " +
       "better but could not be applied to this cohort."],
    ];
    if (ds && ds.consequence) {
      items.splice(3, 0, ["Measured domain shift between the two cohorts", ds.consequence.note]);
    }
    // limitations that specifically block the low-R² argument from being publication-ready
    const blockers = [
      ["Attenuation bounds ONE variance source only",
       "The simulation quantifies how much clock error alone can shrink an R². It does not account " +
       "for transcript counting noise, probe-set and batch effects, or biological heterogeneity. A " +
       "low observed R² is therefore NOT explained by this analysis — only shown to be compatible " +
       "with a real trend."],
      ["The attenuation curve uses live-cohort error on a fixed cohort",
       "Residuals come from live-imaged embryos. If fixed-image segmentation is noisier — likely, " +
       "and unmeasured — the true attenuation is stronger than simulated, so the curve is optimistic."],
      ["MERFISH panels are disjoint, so n is small and gene-dependent",
       "Each gene is measured only in the subset of zygotes whose panel contains it. The " +
       "confirmatory gene is available in roughly 20 zygotes, which is too few to separate " +
       "attenuation from sampling noise with any precision."],
    ];
    if (conf && conf.status === "estimable" && (conf.robustness || {}).verdict &&
        conf.robustness.verdict.startsWith("OUTLIER")) {
      blockers.unshift(["The confirmatory result is outlier-driven",
        `${nc.downstream_illustration.confirmatory_gene} reaches R² = ${fmt(conf.observed_r2)} only ` +
        `because one embryo carries ${conf.robustness.max_count} counts ` +
        `(${conf.robustness.count_ratio_max_to_next}× the next highest). Dropping that single ` +
        `embryo moves R² by ${fmt(conf.robustness.max_single_point_influence_on_r2)}. The burst is ` +
        "biologically plausible for a ZGA element at late τ, but a publication claim cannot rest " +
        "on one point; more zygotes carrying this gene are required."]);
    }
    $("#lim-body").innerHTML =
      `<h2 style="margin-bottom:10px">General limitations</h2>` +
      items.map(([t, b]) => `<div class="pt-lim"><b>${esc(t)}</b><span>${esc(b)}</span></div>`).join("") +
      `<h2 style="margin:20px 0 10px">Specifically blocking a publication-ready low-R² argument</h2>` +
      blockers.map(([t, b]) => `<div class="pt-lim" style="border-left-color:#b42318"><b>${esc(t)}</b><span>${esc(b)}</span></div>`).join("");
  }
})();
