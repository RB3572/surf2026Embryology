/* Render Check — every 3-D render in the deck, redrawn here and its printed number recomputed.
 *
 * A render is the figure a reader trusts most and can check least: one embryo, one picture, one
 * number. A caption that disagrees with its own data is invisible in the picture, so every panel
 * here carries a verdict.
 *
 * The scenes are the site's own (data/segments/<Stage>__<id>.json.gz) and the transcripts are
 * coloured by whatever rule that panel's figure used — blastomere, plane side, or segment label.
 *
 * Data: data/renders.json (build_renders.py) + the scene files, loaded on demand.
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  // The scene's own `extents` — and therefore V.sceneLayout's axis ranges — are in the viewer's
  // isotropic PLOT space, so everything drawn here stays in plot units and only the build's µm
  // geometry is converted in. µm = plot x 0.15 on all three axes.
  const PX = 0.15;
  const HI = "#b01b12", LO = "#14539e", OTHER = "#9aa3b2", GONE = "#cbd5e1";
  const PLAIN = "#c2410c";      // a panel with no split rule still needs a visible cloud

  const state = { doc: null, byId: {}, cur: null, scene: null, cache: {} };

  (async function init() {
    try { state.doc = await (await fetch("data/renders.json")).json(); }
    catch (e) {
      $("#placeholder").innerHTML =
        `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_renders.py</code>.</div></div>`;
      return;
    }
    const m = state.doc.meta;
    state.doc.panels.forEach((p) => (state.byId[p.id] = p));
    $("#panel-count").textContent =
      `${m.n_panels} renders · ${m.n_agree} agree · ${m.n_disagree} disagree`;
    V.buildTabs($("#tabs"), state.doc.panels, select, (p) => ({
      label: `${p.fig} ${p.title}`, sub: p.label || p.embryo,
      title: `${p.fig} · ${p.title} · ${p.label || p.embryo} · ${p.verdict}`,
      cls: p.verdict === "disagrees" ? "tab-warn"
         : p.verdict === "display sub-sample" ? "tab-note" : "",
    }));
    wire();
    renderSummary();
    select(state.doc.panels[0].id);
  })();

  async function select(id) {
    state.cur = id;
    V.markActiveTab($("#tabs"), id);
    const p = state.byId[id];
    $("#loading").hidden = false;
    $("#loading-text").textContent = `Loading ${p.label || p.embryo}…`;
    let sc = state.cache[p.scene];
    try {
      if (!sc) { sc = await V.loadGz(`data/segments/${p.scene}`); state.cache[p.scene] = sc; }
    } catch (_) { $("#loading").hidden = true; return; }
    if (state.cur !== id) return;
    state.scene = sc;
    $("#controls").hidden = false; $("#placeholder").hidden = true; $("#drawer").hidden = false;
    $("#loading").hidden = true;
    $("#drawer-panel").textContent = `${p.fig} · ${p.title}`;
    render3D(); renderReadout();
    // the drawer stays shut: the render is the thing to look at, and the table is one click away
    if ($("#drawer").dataset.open === "true") renderSummary();
  }

  // ───────── the scene ─────────
  function txPositions(sc, gene) {
    const t = sc.transcripts[gene];
    if (!t) return null;
    const zs = sc.z_scale, n = t.x.length;
    const x = new Array(n), y = new Array(n), z = new Array(n);
    for (let i = 0; i < n; i++) {
      x[i] = t.x[i]; y[i] = t.y[i]; z[i] = t.gz[i] * zs;
    }
    return { x, y, z, s: t.s, n };
  }

  // the colour rule IS the figure — each family split its transcripts a different way, and using
  // one generic colouring would quietly hide what each panel was about
  function sideOf(p, T, i) {
    const h = p.highlight || {};
    if (h.kind === "blastomeres") {
      const s = String(T.s[i]);
      return s === String(h.hi) ? 0 : s === String(h.lo) ? 1 : 2;
    }
    if (h.kind === "labels") return (h.labels || []).includes(String(T.s[i])) ? 0 : 2;
    if (h.kind === "plane") {
      // the build's origin is µm; the points here are plot units, so bring the origin across
      const o = h.origin;
      const d = (T.x[i] * PX - o[0]) * h.normal[0] + (T.y[i] * PX - o[1]) * h.normal[1] +
                (T.z[i] * PX - o[2]) * h.normal[2];
      return d > 0 ? 0 : 1;
    }
    return 2;
  }

  function render3D() {
    const sc = state.scene, p = state.byId[state.cur];
    if (!sc || !p) return;
    const showBody = $("#show-body").checked, showTx = $("#show-tx").checked;
    const showOther = $("#show-other").checked;
    const h = p.highlight || {};
    const traces = [];
    if (showBody) {
      // THE BODIES COME FROM V.bodyTraces, not from a local mesh builder. Rolling one here is
      // what made this page look unlike every other one: it drew the cytoplasm flat grey at
      // opacity 0.06, against the house 0.13 and the scene's own per-segment colour, and it
      // could not follow the dark-render toggle at all. Only the highlighted labels are
      // recoloured on top, so the panel still says which half is which.
      for (const t of V.bodyTraces(sc)) {
        const lbl = String(t.name).replace(/^body M/, "");
        const isHi = h.kind === "blastomeres" && lbl === String(h.hi);
        const isLo = h.kind === "blastomeres" && lbl === String(h.lo);
        const isMark = h.kind === "labels" && (h.labels || []).includes(lbl);
        if (isHi || isLo) {
          t.color = isHi ? HI : LO;
          t.opacity = Math.max(t.opacity, 0.15);
          t.name = isHi ? "higher blastomere" : "lower blastomere";
        } else if (isMark) {
          t.color = "#7c3aed";
          t.opacity = Math.max(t.opacity, 0.34);
          t.name = `region ${lbl} (highlighted)`;
        }
        t.hoverinfo = "name";
        traces.push(t);
      }
    }
    if (showTx) {
      const gene = p.genes[0];
      const T = txPositions(sc, gene);
      if (T) {
        const groups = [[], [], []];
        for (let i = 0; i < T.n; i++) groups[sideOf(p, T, i)].push(i);
        const NAME = h.kind === "blastomeres" ? ["higher blastomere", "lower blastomere", "elsewhere"]
                   : h.kind === "labels" ? ["inside the marked region", "—", "outside it"]
                   : h.kind === "plane" ? ["side a", "side b", "outside the cell"]
                   : [gene, "—", "—"];
        // With no split rule (1.7 draws one gene, whole embryo) every point lands in group 2,
        // and painting those "elsewhere" grey left the cloud almost invisible. Grey only means
        // "outside the region this panel is about", which needs a region to be about.
        const plain = !h.kind;
        [[0, plain ? PLAIN : HI], [1, LO], [2, plain ? PLAIN : GONE]].forEach(([g, col]) => {
          if (!groups[g].length) return;
          traces.push({ type: "scatter3d", mode: "markers",
            name: plain ? `${gene} · ${groups[g].length}` : `${gene} · ${NAME[g]}`,
            x: groups[g].map((i) => T.x[i]), y: groups[g].map((i) => T.y[i]),
            z: groups[g].map((i) => T.z[i]),
            marker: { size: 2.6, color: col,
                      opacity: (g === 2 && !plain) ? 0.35 : 0.85, line: { width: 0 } },
            hovertemplate: `${gene}<extra></extra>`, legendrank: 10 + g });
        });
      }
      if (showOther) {
        const x = [], y = [], z = [];
        const zs = sc.z_scale;
        for (const g of Object.keys(sc.transcripts)) {
          if (g === gene) continue;
          const t = sc.transcripts[g];
          for (let i = 0; i < t.x.length; i += 3) {          // every third, for the browser's sake
            x.push(t.x[i]); y.push(t.y[i]); z.push(t.gz[i] * zs);
          }
        }
        traces.push({ type: "scatter3d", mode: "markers", name: "every other gene (1 in 3)",
          x, y, z, marker: { size: 1.3, color: OTHER, opacity: 0.18 },
          hoverinfo: "skip", legendrank: 60 });
      }
    }
    // the cut, drawn as a quad, when the panel is defined by a plane
    if (h.kind === "plane") {
      const n = h.normal, o = h.origin.map((v) => v / PX);     // µm -> plot units
      let a = Math.abs(n[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
      const e1 = norm(cross(n, a)), e2 = norm(cross(n, e1)), L = 55 / PX;
      const c = [[-1, -1], [1, -1], [1, 1], [-1, 1]].map(([u, v]) =>
        [0, 1, 2].map((k) => o[k] + (e1[k] * u + e2[k] * v) * L));
      traces.push({ type: "mesh3d", x: c.map((q) => q[0]), y: c.map((q) => q[1]),
        z: c.map((q) => q[2]), i: [0, 0], j: [1, 2], k: [2, 3], color: "#111827",
        opacity: 0.14, name: "the deck's plane", hoverinfo: "name", legendrank: 5 });
    }
    Plotly.react($("#plot-host"), traces, V.sceneLayout(sc.extents, sc.id),
      V.plotConfig);
  }
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                           a[0] * b[1] - a[1] * b[0]];
  const norm = (v) => { const n = Math.hypot(...v) || 1; return v.map((x) => x / n); };

  // ───────── the verdict ─────────
  function badge(v) {
    const cls = v === "agrees" ? "ok" : v === "disagrees" ? "bad" : "note";
    return `<span class="rn-badge ${cls}">${v}</span>`;
  }
  const nearMiss = (c) => (c.diff != null && c.diff <= 2) || (c.rel != null && c.rel < 0.01);

  function renderReadout() {
    const p = state.byId[state.cur];
    const L = [`<div class="rn-fig">${p.fig}</div><div class="rn-title">${p.title}</div>`,
      `<div class="rn-emb">${p.label || p.embryo}</div>`, badge(p.verdict)];
    p.readout.forEach((r) =>
      L.push(`<div class="rn-line"><span class="k">${r.k}</span><span class="v">${r.v}</span></div>`));
    if (p.checks.length) {
      L.push(`<div class="rn-checks">` + p.checks.map((c) =>
        `<div class="rn-check ${c.ok ? "ok" : nearMiss(c) ? "near" : "bad"}">
           <span class="n">${c.name}</span>
           <span class="d">${c.ok ? "✓" : `deck ${c.deck} · ours ${c.ours}`}</span></div>`).join("") +
        `</div>`);
      const bad = p.checks.filter((c) => !c.ok);
      if (bad.length && bad.every(nearMiss)) {
        L.push(`<div class="rn-note"><b>Off by ${bad.map((c) => c.diff).join(" and ")}
          molecule${bad[0].diff === 1 ? "" : "s"}.</b> A transcript sitting exactly on a cutting
          plane goes one way here and the other way there; that is a tie-break convention, not a
          disagreement about the measurement.</div>`);
      }
    }
    if (p.caveat) L.push(`<div class="rn-caveat">${p.caveat}</div>`);
    $("#rn-readout").innerHTML = L.join("");
  }

  function renderSummary() {
    const m = state.doc.meta;
    $("#rn-summary").innerHTML =
      `<b>${m.n_agree} of ${m.n_panels} renders reproduce exactly.</b>
       ${m.n_disagree} disagree — ${m.n_near_miss} of those by a molecule or two, which is a
       transcript landing on a cutting plane rather than a different measurement.
       ${m.n_subsampled} panels print counts that are a <i>drawing device</i> (sub-sampled so a
       series reads consistently), so there is nothing to check them against and the page says so
       instead of scoring them.
       <br><br><b>Counts must match exactly</b> — cytoplasm by segment label, never a containment
       test — and folds within ${(0.02 * 100).toFixed(0)}%, which is the precision they are
       printed at. Folds are ratios of <b>densities</b>: sister blastomeres are routinely 20%
       apart in volume, and a ratio of raw counts would read that difference as expression.`;
    const rows = state.doc.panels;
    $("#rn-table").innerHTML =
      `<table class="rn-tab"><thead><tr><th>fig</th><th>panel</th><th>embryo</th><th>gene</th>
        <th>verdict</th><th>what was checked</th></tr></thead><tbody>` +
      rows.map((p) => `<tr data-id="${p.id}" class="${p.id === state.cur ? "on" : ""}">
        <td>${p.fig}</td><td class="g">${p.title}</td><td class="e">${p.label || p.embryo}</td>
        <td>${p.genes.join(", ")}</td><td>${badge(p.verdict)}</td>
        <td class="c">${p.checks.length
          ? p.checks.map((c) => `${c.name}: ${c.ok ? "✓" : `<b>${c.deck} vs ${c.ours}</b>`}`).join(" · ")
          : "nothing to check — the printed numbers are a drawing device"}</td></tr>`).join("") +
      `</tbody></table>`;
    $("#rn-table").querySelectorAll("tr[data-id]").forEach((tr) =>
      tr.addEventListener("click", () => select(tr.dataset.id)));
  }

  function openDrawer(open) {
    $("#drawer").dataset.open = open ? "true" : "false";
    $("#drawer-handle").setAttribute("aria-expanded", String(open));
    if (open) renderSummary();
  }

  function wire() {
    ["#show-body", "#show-tx", "#show-other"].forEach((s) =>
      $(s).addEventListener("change", render3D));
    $("#drawer-handle").addEventListener("click", () =>
      openDrawer($("#drawer").dataset.open !== "true"));
    const rz = $("#drawer-resize");
    let sh = 0, dy = 0, drag = false;
    rz.addEventListener("pointerdown", (ev) => {
      drag = true; sh = $("#drawer").getBoundingClientRect().height; dy = ev.clientY;
      rz.setPointerCapture(ev.pointerId); ev.preventDefault();
    });
    rz.addEventListener("pointermove", (ev) => {
      if (!drag) return;
      $("#drawer").style.setProperty("--drawer-h",
        Math.max(200, Math.min(window.innerHeight - 90, sh + (dy - ev.clientY))) + "px");
    });
    const end = () => { drag = false; };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
    V.wireWindow($("#controls"), $("#controls-header"),
                 [...$("#controls").querySelectorAll(".rz")], "renders_controls_box");
  }
})();
