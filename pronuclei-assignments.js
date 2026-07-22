/* Pronuclei Assignments — assign male/female to each zygote's two pronuclei by several tests and
 * take the consensus. Built on viewer-core.js (VCore). The 3-D view shows the two pronuclei (coloured
 * by the consensus call), the polar body, the sperm, and the per-test distance lines. The bottom
 * drawer is a concordance grid: rows = tests, columns = the two pronuclei, ♀ pink / ♂ blue per cell,
 * consensus at the bottom (half-and-half on an even split). Data from build_pronuclei_assignments.py;
 * the 3-D meshes are reused from the Pronuclei project's scenes. */
(() => {
  const $ = (s) => document.querySelector(s);
  const V = window.VCore;
  const F = "#ec4899", M = "#3b82f6", PB_C = "#9aa3b2", SP_C = "#f59e0b", LINE_C = "#111827";
  const TESTS = [
    { key: "pb_com", label: "Closer to polar body · centre of mass", line: "#0ea5e9" },
    { key: "pb_shell", label: "Closer to polar body · nearest shells", line: "#8b5cf6" },
    { key: "volume", label: "Smaller pronucleus = female" },
    { key: "sperm", label: "Closer to sperm = male", line: SP_C },
  ];

  const tabsEl = $("#tabs"), countEl = $("#embryo-count");
  const controlsEl = $("#controls"), plotHost = $("#plot-host");
  const placeholder = $("#placeholder"), loadingEl = $("#loading"), loadingTxt = $("#loading-text");
  const paReadout = $("#pa-readout"), paEmb = $("#pa-emb"), paGrid = $("#pa-grid"), paGridNote = $("#pa-grid-note");
  const paColor = $("#pa-color"), paPb = $("#pa-pb"), paSperm = $("#pa-sperm"), paLines = $("#pa-lines");
  const drawer = $("#drawer"), drawerHandle = $("#drawer-handle");

  const state = { points: [], byId: {}, assign: {}, currentId: null, scene: null, rec: null, drawerOpen: false, dotSize: 1.5 };
  let vcExtras = null;

  (async function init() {
    try {
      const [m, asg] = await Promise.all([
        (await fetch("data/pronuclei_manifest.json")).json(),
        (await fetch("data/pronuclei_assignments.json")).json(),
      ]);
      state.points = m.embryos.filter((p) => asg.embryos.some((a) => a.id === p.id));
      state.points.forEach((p) => (state.byId[p.id] = p));
      asg.embryos.forEach((a) => (state.assign[a.id] = a));
      countEl.textContent = `${state.points.length} zygotes · male/female pronuclei by consensus`;
      V.buildTabs(tabsEl, state.points, selectEmbryo, (e) => {
        const a = state.assign[e.id] || {};
        return { label: e.label, sub: e.date_short,
          title: `${e.label} · consensus female = pronucleus ${(a.consensus || {}).female === 0 ? "1" : "2"}` };
      });
      wireControls(); wireDrawer();
      vcExtras = V.addWindowExtras($("#controls-body"), { defaultSize: state.dotSize, onDotSize: (s) => { state.dotSize = s; if (state.scene) render(); } });
    } catch (err) { showError("Failed to load: " + (err.message || err)); }
  })();

  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    V.markActiveTab(tabsEl, id);
    const meta = state.byId[id] || {};
    showLoading(`Loading ${meta.label || id}…`);
    try {
      const scene = await V.loadGz(`data/pronuclei/${id}.json.gz`);
      if (state.currentId !== id) return;
      state.scene = scene; state.rec = state.assign[id] || null;
      if (vcExtras) vcExtras.setAtlas(id);
      controlsEl.hidden = false; placeholder.hidden = true; drawer.hidden = false;
      render(); renderReadout(); renderGrid();
      if (!state.drawerOpen) openDrawer(true);
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
    return { type: "mesh3d", x, y, z, i: ii, j: jj, k: kk, color, opacity, name, showlegend: true,
      flatshading: false, hoverinfo: "name", lighting: { ambient: 0.7, diffuse: 0.55, specular: 0.12, roughness: 0.9 }, legendrank: lbl };
  }
  const pronColor = (idx) => { const c = (state.rec || {}).consensus || {}; if (!paColor.checked) return idx === 0 ? "#2563eb" : "#dc2626";
    if (c.split) return "#a855f7"; return c.female === idx ? F : M; };
  function render() {
    const s = state.scene, r = state.rec; if (!s || !r) return;
    const pl = r.pron.map((p) => p.label);
    const traces = [];
    for (const lbl of s.mask_labels) {
      const pi = pl.indexOf(lbl);
      const isPb = r.polar && lbl === r.polar.label;
      let color = "#9aa3b2", op = 0.07, name = `Segment ${lbl}`;
      if (pi >= 0) { color = pronColor(pi); op = 0.55; name = `Pronucleus ${pi + 1} · ${paColor.checked ? (((r.consensus || {}).split) ? "split" : ((r.consensus || {}).female === pi ? "♀ female" : "♂ male")) : "seg " + lbl}`; }
      else if (isPb && paPb.checked) { color = PB_C; op = 0.32; name = `Polar body`; }
      else if (isPb) continue;
      const t = segMesh(s, lbl, color, op, name); if (t) traces.push(t);
    }
    // assignment distance lines
    if (paLines.checked) TESTS.forEach((t) => {
      const rt = r.tests[t.key]; if (!rt || !rt.line) return;
      const [a, b] = rt.line;
      traces.push({ type: "scatter3d", mode: "lines+markers", name: t.label,
        x: [a[0], b[0]], y: [a[1], b[1]], z: [a[2], b[2]],
        line: { color: t.line || LINE_C, width: 5 }, marker: { size: 3, color: t.line || LINE_C },
        hovertemplate: `${t.label}<extra></extra>`, legendrank: 100 });
    });
    // sperm marker
    if (paSperm.checked && r.sperm_plot) {
      const sp = r.sperm_plot;
      traces.push({ type: "scatter3d", mode: "markers", name: "sperm",
        x: [sp[0]], y: [sp[1]], z: [sp[2]], marker: { size: 6, color: SP_C, symbol: "diamond", line: { color: "#fff", width: 1 } },
        hovertemplate: "sperm<extra></extra>", legendrank: 200 });
    }
    Plotly.react(plotHost, traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }
  function renderReadout() {
    const r = state.rec; if (!r) return;
    const c = r.consensus || {};
    const fem = c.split ? "even split" : `pronucleus ${c.female + 1}`;
    paReadout.innerHTML =
      `<div class="pn-big"><span style="color:${c.split ? "#a855f7" : F}">${fem}</span> <span class="pn-lbl">consensus female</span></div>` +
      `<div class="pn-resid">pronucleus 1 = seg <b>${r.pron[0].label}</b> (${r.pron[0].volume.toLocaleString()} µm³, larger) · pronucleus 2 = seg <b>${r.pron[1].label}</b> (${r.pron[1].volume.toLocaleString()} µm³)</div>` +
      `<div class="pn-resid">polar body: ${r.polar ? `seg <b>${r.polar.label}</b>` : "not detected"} · sperm: ${r.sperm_plot ? "located" : "not available"}</div>`;
  }

  // ---------- concordance grid ----------
  const cell = (isFemale) => isFemale == null
    ? `<div class="pa-cell pa-na" title="test unavailable"></div>`
    : `<div class="pa-cell ${isFemale ? "pa-f" : "pa-m"}">${isFemale ? "♀" : "♂"}</div>`;
  function renderGrid() {
    const r = state.rec; if (!r) { paGrid.innerHTML = ""; return; }
    paEmb.textContent = `· ${r.label}`;
    const cols = [`Pronucleus 1<span class="pa-sub">seg ${r.pron[0].label}</span>`, `Pronucleus 2<span class="pa-sub">seg ${r.pron[1].label}</span>`];
    let html = `<div class="pa-corner">test \\ pronucleus</div>` + cols.map((c) => `<div class="pa-collabel">${c}</div>`).join("");
    TESTS.forEach((t) => {
      const rt = r.tests[t.key];
      html += `<div class="pa-rowlabel" title="${t.label}">${t.label}</div>`;
      if (!rt) html += cell(null) + cell(null);
      else { const fem = rt.female; html += cell(fem === 0) + cell(fem === 1); }
    });
    // consensus row
    const c = r.consensus || {};
    html += `<div class="pa-rowlabel pa-consensus">Consensus</div>`;
    if (c.split) html += `<div class="pa-cell pa-split" title="even split">♀♂</div><div class="pa-cell pa-split" title="even split">♀♂</div>`;
    else html += cell(c.female === 0) + cell(c.female === 1);
    paGrid.innerHTML = html;
    const votes = TESTS.map((t) => r.tests[t.key]).filter(Boolean).length;
    paGridNote.innerHTML = c.split
      ? `<b>Even split</b> — ${c.n0} test(s) call pronucleus 1 female, ${c.n1} call pronucleus 2. No consensus.`
      : `<b>Consensus: pronucleus ${c.female + 1} is female</b> — ${Math.max(c.n0, c.n1)} of ${votes} tests agree.`;
  }

  // ---------- drawer ----------
  function openDrawer(open) {
    state.drawerOpen = open;
    drawer.dataset.open = open ? "true" : "false";
    drawerHandle.setAttribute("aria-expanded", String(open));
  }
  function wireControls() {
    [paColor, paPb, paSperm, paLines].forEach((c) => c.addEventListener("change", () => { if (state.scene) render(); }));
    const corners = [...controlsEl.querySelectorAll(".rz")];
    try { V.wireWindow(controlsEl, $("#controls-header"), corners, "pa.win"); } catch (_) {}
  }
  function wireDrawer() {
    let start = null, moved = false;
    drawerHandle.addEventListener("pointerdown", (e) => { if (e.button && e.button !== 0) return; start = { x: e.clientX, y: e.clientY }; moved = false; try { drawerHandle.setPointerCapture(e.pointerId); } catch (_) {} });
    drawerHandle.addEventListener("pointermove", (e) => { if (!start) return;
      if (!moved && Math.hypot(e.clientX - start.x, e.clientY - start.y) < 5) return;
      if (!moved) { moved = true; drawer.classList.add("dragging"); if (drawer.dataset.open !== "true") openDrawer(true); }
      drawer.style.setProperty("--drawer-h", Math.max(200, Math.min(window.innerHeight - 100, window.innerHeight - e.clientY - 40)) + "px"); e.preventDefault(); });
    const up = (e) => { if (!start) return; try { drawerHandle.releasePointerCapture(e.pointerId); } catch (_) {}
      if (moved) drawer.classList.remove("dragging"); else openDrawer(drawer.dataset.open !== "true"); start = null; moved = false; };
    drawerHandle.addEventListener("pointerup", up); drawerHandle.addEventListener("pointercancel", up);
    const rz = $("#drawer-resize"); let sh = 0;
    rz.addEventListener("pointerdown", (e) => { sh = window.innerHeight - drawer.getBoundingClientRect().top; rz._d = { y: e.clientY }; rz.setPointerCapture(e.pointerId); e.preventDefault(); });
    rz.addEventListener("pointermove", (e) => { if (!rz._d) return; drawer.style.setProperty("--drawer-h", Math.max(200, Math.min(window.innerHeight - 100, sh - (e.clientY - rz._d.y))) + "px"); });
    const end = (e) => { if (rz._d) { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} } };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);
  }

  function showLoading(t) { loadingTxt.textContent = t; loadingEl.hidden = false; }
  function hideLoading() { loadingEl.hidden = true; }
  function showError(msg) { placeholder.hidden = false;
    placeholder.innerHTML = `<div class="ph-inner"><div class="ph-title" style="color:#c0392b">Error</div><div class="ph-sub">${msg}</div></div>`; }
})();
