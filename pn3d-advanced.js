/* Pronuclear Structure Clock — ADVANCED ANALYSIS view + data loading.
 *
 * Loads the frozen artifacts, hands them to the guided narrative, and renders the
 * two advanced sections: a per-embryo explorer (every embryo, with its input
 * planes, geometry, posterior and QC) and the model-development view (provenance,
 * validation, baselines, domain shift, evidence table, failure cases).
 */
window.PN3DAdvanced = (() => {
  const $ = (s, r = document) => r.querySelector(s);
  const el = (t, c, h) => { const e = document.createElement(t); if (c) e.className = c; if (h != null) e.innerHTML = h; return e; };
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>]/g, (m) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[m]));
  const f = (v, d = 3) => (v == null || Number.isNaN(v)) ? "—" : (+v).toFixed(d);
  let M, INF, GEO, PREV, MAN, built = false, curId = null;

  const zygotes = () => INF.embryos.filter((e) => e.stage === "zygote");
  const recOf = (id) => INF.embryos.find((e) => e.embryo_id === id);
  const prevOf = (id) => (PREV.embryos || []).find((e) => e.embryo_id === id);
  const nPN = (r) => (r && r.geometry && r.geometry.n_pronuclei) || 0;

  Promise.all([
    fetch("data/pn3d/model.json").then((r) => r.json()),
    fetch("data/pn3d/inference.json").then((r) => r.json()),
    fetch("data/pn3d/segmentation_geometry.json").then((r) => r.json()),
    fetch("data/pn3d/preview_index.json").then((r) => r.ok ? r.json() : { embryos: [] }).catch(() => ({ embryos: [] })),
    fetch("data/pn3d/manifest.json").then((r) => r.ok ? r.json() : null).catch(() => null),
  ]).then(([m, inf, geo, prev, man]) => {
    M = m; INF = inf; GEO = geo; PREV = prev; MAN = man;
    const c = m.counts;
    $("#pt-subtitle").textContent =
      `${m.package_version} · ${c.zygotes_scored ?? c.zygotes_resolved}/${c.zygotes_audited} zygotes scored · clock MAE ${f(m.clock.cv_metrics.mae)}`;
    buildNav();
    window.PN3DGuided.init({ model: m, inf, geo, prev, man });
  }).catch((e) => {
    document.querySelector(".gd-intro").insertAdjacentHTML("beforeend",
      `<p class="gd-warn">Could not load the model artifacts — run <code>python3 scripts/build_pn3d_model.py</code>. ${esc(e.message)}</p>`);
  });

  function buildNav() {
    const nav = $("#pt-nav"); nav.innerHTML = "";
    [["explorer", "Per-embryo explorer"], ["model", "Model development"]].forEach(([k, label], i) => {
      const b = el("button", "pt-navb" + (i === 0 ? " active" : ""), esc(label));
      b.onclick = () => {
        nav.querySelectorAll(".pt-navb").forEach((x) => x.classList.toggle("active", x === b));
        const t = document.querySelector(`.pt-sec[data-sec="${k}"]`);
        if (t) t.scrollIntoView({ block: "start", behavior: "smooth" });
      };
      nav.appendChild(b);
    });
  }

  function show() { if (!built) { buildExplorer(); buildDev(); built = true; } }

  // ───────────────────────── per-embryo explorer ─────────────────────────
  function buildExplorer() {
    const host = $("#adv-explorer"); host.innerHTML = "";
    const rows = zygotes().slice().sort((a, b) =>
      (a.inferable ? a.pseudotime.tau_mean : 9e9) - (b.inferable ? b.pseudotime.tau_mean : 9e9));

    const wrap = el("div", "pt-two");
    const left = el("div");
    left.appendChild(el("div", "pt-sub", `All ${rows.length} zygotes · click one`));
    const sc = el("div", "pn-scroll");
    const t = el("table", "pn-table");
    t.innerHTML = "<thead><tr><th>embryo</th><th>PN</th><th class='num'>τ</th><th class='num'>rms/R</th><th>domain</th></tr></thead>";
    const tb = el("tbody");
    rows.forEach((r) => {
      const tr = el("tr");
      tr.innerHTML = `<td>${esc(r.embryo_id)}</td>
        <td>${nPN(r) || "—"}</td>
        <td class="num">${r.inferable ? f(r.pseudotime.tau_mean, 3) : "—"}</td>
        <td class="num">${r.geometry ? f(r.geometry.rms_over_R, 3) : "—"}</td>
        <td>${window.PN3DGuided.oodPill(r)}</td>`;
      tr.onclick = () => detail(r.embryo_id);
      tb.appendChild(tr);
    });
    t.appendChild(tb); sc.appendChild(t); left.appendChild(sc);
    wrap.appendChild(left);
    const right = el("div"); right.id = "adv-detail"; wrap.appendChild(right);
    host.appendChild(wrap);
    detail(rows[Math.floor(rows.length / 2)].embryo_id);
  }

  function detail(id) {
    curId = id;
    const r = recOf(id), host = $("#adv-detail"); if (!r || !host) return;
    host.innerHTML = "";
    const g = r.geometry || {}, p = r.pseudotime;
    host.appendChild(el("div", "pt-sub", esc(id)));

    const pl = el("div", "pn-planes");
    const pv = prevOf(id);
    [["xy", "XY"], ["xz", "XZ"], ["yz", "YZ"]].forEach(([k, lab]) => {
      const d = el("div", "pn-plane");
      if (pv && pv.planes && pv.planes[k]) {
        d.innerHTML = `<img src="${esc(pv.planes[k].seg)}" alt="${lab}" loading="lazy"><div class="pn-cap">${lab}</div>`;
      } else { d.className = "pn-plane pn-missing"; d.textContent = `${lab} — no preview`; }
      pl.appendChild(d);
    });
    host.appendChild(pl);

    const kv = [
      ["pronuclei found", nPN(r) || "—"],
      ["cell radius", g.cell_radius_um != null ? f(g.cell_radius_um, 1) + " µm" : "—"],
      ["rms → centre", g.rms_to_center_um != null ? f(g.rms_to_center_um, 1) + " µm" : "—"],
      ["rms ÷ R (model input)", f(g.rms_over_R, 3)],
      ["pronuclear volume fraction", f(g.pron_vol_frac, 4)],
      ["polar body", g.polar_body_present ? "found (external)" : "not found"],
    ];
    if (p) kv.push(["τ", `${f(p.tau_mean, 3)} ± ${f(p.tau_sd, 3)}`],
                   ["50% interval", `${f(p.interval_50[0], 2)} – ${f(p.interval_50[1], 2)}`],
                   ["95% interval", `${f(p.interval_95[0], 2)} – ${f(p.interval_95[1], 2)}`],
                   ["confidence", f(p.confidence, 2)]);
    host.appendChild(el("div", "gd-eq", kv.map(([k, v]) =>
      `<div class="gd-eq-row"><span class="gd-eq-k">${esc(k)}</span><span class="gd-eq-v">${esc(String(v))}</span></div>`).join("")));

    if (p && p.if_second_pronucleus_present) {
      const s = p.if_second_pronucleus_present;
      host.appendChild(el("p", "gd-warn",
        `<b>Single annotated pronucleus.</b> If a second pronucleus were present but unlabelled, τ would ` +
        `lie between <b>${f(s.tau_lo, 2)}</b> and <b>${f(s.tau_hi, 2)}</b> (median ${f(s.tau_median, 2)}), ` +
        `sweeping the missing pronucleus over every distance seen in the two-pronucleus cohort. ` +
        `The calibrated interval above does <i>not</i> cover that annotation uncertainty.`));
    }
    const fl = (r.flags || []).concat(r.ood_reasons || []);
    if (fl.length) host.appendChild(el("ul", "gd-lims",
      fl.map((x) => `<div class="gd-lim"><span>${esc(x)}</span></div>`).join("")));
  }

  // ───────────────────────── model development ─────────────────────────
  function buildDev() {
    const host = $("#adv-dev"); host.innerHTML = "";
    const c = M.counts, cv = M.clock.cv_metrics;

    host.appendChild(el("div", "pt-sub", "Data provenance & coverage"));
    host.appendChild(el("div", "gd-eq", [
      ["experiments / batches", MAN ? MAN.n_experiments : "—"],
      ["embryos inventoried", MAN ? MAN.n_embryos : "—"],
      ["zygotes audited", c.zygotes_audited],
      ["zygotes scored", c.zygotes_scored ?? c.zygotes_resolved],
      ["… with two pronuclei", c.zygotes_resolved],
      ["… single pronucleus (extrapolated)", c.zygotes_single_pronucleus ?? 0],
      ["flagged out-of-domain", c.zygotes_ood],
      ["2-cell OOD reference", c.two_cell_reference],
      ["time supervision", "live imaging (Scheffler 2021)"],
    ].map(([k, v]) => `<div class="gd-eq-row"><span class="gd-eq-k">${esc(k)}</span><span class="gd-eq-v">${esc(String(v))}</span></div>`).join("")));

    host.appendChild(el("div", "pt-sub", "Clock validation — the only independent time test"));
    host.appendChild(el("div", "gd-eq", Object.entries({
      "leave-one-embryo-out MAE": f(cv.mae, 4),
      "Spearman ρ": f(cv.spearman, 4),
      "within-embryo monotonicity": f(cv.within_embryo_mono_median, 4),
      "coverage 50 / 80 / 95%": `${f(cv.coverage_50, 2)} / ${f(cv.coverage_80, 2)} / ${f(cv.coverage_95, 2)}`,
      "calibration scale": f(cv.calibration_scale, 3),
      "frames / embryos": `${cv.n} / ${cv.n_embryos}`,
      "model input": M.clock.input_feature || M.clock.feature,
    }).map(([k, v]) => `<div class="gd-eq-row"><span class="gd-eq-k">${esc(k)}</span><span class="gd-eq-v">${esc(String(v))}</span></div>`).join("")));

    host.appendChild(el("div", "pt-sub", "Baselines"));
    const bt = el("table", "pn-table");
    bt.innerHTML = "<thead><tr><th>model</th><th>description</th><th>result</th></tr></thead>";
    const bb = el("tbody");
    Object.entries(M.baselines).forEach(([k, v]) => {
      const res = v.scheffler_cv_mae != null ? `CV MAE ${f(v.scheffler_cv_mae, 3)}`
        : `${v.status || ""}${v.spearman_approx ? ` (ρ≈${v.spearman_approx})` : ""}`;
      bb.appendChild(el("tr", null, `<td><b>${esc(k)}</b></td><td>${esc(v.description || "")}</td><td>${esc(res)}</td>`));
    });
    bt.appendChild(bb); host.appendChild(bt);

    host.appendChild(el("div", "pt-sub", "Domain: fixed MERFISH vs live imaging"));
    const ds = M.domain_shift;
    host.appendChild(el("div", "gd-eq", [
      ["fixed median (physical)", ds.physical_distance_sum.fixed_median + " µm"],
      ["live median (physical)", ds.physical_distance_sum.live_median + " µm"],
      ["KS distance (physical)", ds.physical_distance_sum.ks_statistic],
      ["KS distance (dimensionless)", ds.dimensionless_rms_over_R
        ? ds.dimensionless_rms_over_R.ks_statistic : "—"],
      ["fixed stacks above live support", ds.n_fixed_above_live_support],
    ].map(([k, v]) => `<div class="gd-eq-row"><span class="gd-eq-k">${esc(k)}</span><span class="gd-eq-v">${esc(String(v))}</span></div>`).join("")));
    host.appendChild(el("p", "gd-cap", esc(ds.verdict)));

    host.appendChild(el("div", "pt-sub", "Evidence table"));
    const et = el("table", "pn-table");
    et.innerHTML = "<thead><tr><th>claim</th><th>evidence type</th><th>result</th></tr></thead>";
    const eb = el("tbody");
    M.evidence_table.forEach((e) => eb.appendChild(el("tr", null,
      `<td>${esc(e.claim)}</td><td>${esc(e.evidence_type)}</td><td>${esc(e.result)}</td>`)));
    et.appendChild(eb); host.appendChild(et);

    host.appendChild(el("div", "pt-sub", "Failure cases & flagged embryos"));
    const ft = el("table", "pn-table");
    ft.innerHTML = "<thead><tr><th>embryo</th><th>status</th><th>why</th></tr></thead>";
    const fb = el("tbody");
    zygotes().filter((e) => e.ood_level !== "in_domain").slice(0, 40).forEach((e) => {
      fb.appendChild(el("tr", null,
        `<td>${esc(e.embryo_id)}</td><td>${window.PN3DGuided.oodPill(e)}</td>
         <td>${esc((e.flags || []).concat(e.ood_reasons || []).slice(0, 2).join("; "))}</td>`));
    });
    ft.appendChild(fb); host.appendChild(ft);
  }

  return { show };
})();
