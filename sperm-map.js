/* Sperm Location Browser.
 *
 * Find an embryo by WHERE its sperm ended up. The right drawer is the entry point: pick a stage,
 * pick a measure, and the embryos sort by it — for 2-cell that includes distance to the junction
 * between the two blastomeres, which is the question this was built for.
 *
 * Data: data/sperm_locations.json (build_sperm_locations.py) carries, per embryo, the labelled
 * sperm position and its distance to the cortex / polar body / pronuclei / junction, plus the
 * junction frame itself. The GFP stills are pre-rendered into data/gfp/ because the site is
 * static and cannot reach the acquisition drive.
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const SPERM = V.SPERM_COLOR, MARK = "#7c3aed", PLANE = "#0ea5e9";

  // The GFP frames ship GREYSCALE — colour, gain and black point are applied here, so the two
  // channels can be retuned without rebuilding a single image. 488 is green and 405 blue by
  // default, which is how the labelling tool draws them.
  const GFP_DEF = { 0: { color: "#22c55e", gain: 1, black: 0 },
                    1: { color: "#3b82f6", gain: 1, black: 0 } };
  const GFP_KEY = "surf_sperm_map_gfp";
  const SWATCHES = ["#22c55e", "#22d3ee", "#3b82f6", "#ff2fd0", "#f59e0b", "#ef4444", "#e8ecf5"];

  const state = {
    data: null, byId: {}, stage: "zygote", sub: "all", sort: "cortex", dir: "near",
    currentId: null, scene: null, cache: {}, tab: "gfp",
    spermOn: true, marksOn: true, planeOn: true, gfpMark: true, drawerOpen: false,
    gfpCh: loadGfpPrefs(),
  };

  function loadGfpPrefs() {
    const d = { 0: { ...GFP_DEF[0] }, 1: { ...GFP_DEF[1] } };
    try {
      const s = JSON.parse(localStorage.getItem(GFP_KEY) || "{}");
      [0, 1].forEach((c) => {
        if (!s[c]) return;
        if (typeof s[c].color === "string") d[c].color = s[c].color;
        if (isFinite(s[c].gain)) d[c].gain = Math.min(6, Math.max(0.2, +s[c].gain));
        if (isFinite(s[c].black)) d[c].black = Math.min(0.9, Math.max(0, +s[c].black));
      });
    } catch (_) {}
    return d;
  }
  const saveGfpPrefs = () => {
    try { localStorage.setItem(GFP_KEY, JSON.stringify(state.gfpCh)); } catch (_) {}
  };

  const hexRgb = (h) => {
    const m = /^#([0-9a-f]{6})$/i.exec(String(h).trim());
    if (!m) return [255, 255, 255];
    return [0, 2, 4].map((i) => parseInt(m[1].slice(i, i + 2), 16));
  };

  /* The marker ring has to stay legible whatever tint the channel is wearing — a green ring on a
   * green image is no ring at all. Pick whichever candidate sits furthest from the tint in RGB. */
  const RINGS = ["#ff2d55", "#ffffff", "#00e5ff"];
  function ringFor(tint) {
    const t = hexRgb(tint);
    const dist = (c) => { const p = hexRgb(c);
      return Math.hypot(p[0] - t[0], p[1] - t[1], p[2] - t[2]); };
    return dist(RINGS[0]) >= 150 ? RINGS[0]
         : RINGS.slice(1).reduce((a, b) => (dist(b) > dist(a) ? b : a));
  }

  /* Colour + window/level as an SVG filter, so moving a slider restyles the images already on
   * screen instead of re-fetching them. feComponentTransfer does out = gain·(in − black) on the
   * grey value; feColorMatrix then maps that grey onto the channel's colour. */
  const NS = "http://www.w3.org/2000/svg";
  const gfpFilters = {};
  function ensureFilters() {
    if (gfpFilters[0]) return;
    const svg = document.createElementNS(NS, "svg");
    svg.setAttribute("class", "sm-defs");
    svg.setAttribute("aria-hidden", "true");
    const defs = document.createElementNS(NS, "defs");
    svg.appendChild(defs);
    [0, 1].forEach((ch) => {
      const f = document.createElementNS(NS, "filter");
      f.setAttribute("id", "smf" + ch);
      f.setAttribute("color-interpolation-filters", "sRGB");
      const ct = document.createElementNS(NS, "feComponentTransfer");
      const funcs = ["feFuncR", "feFuncG", "feFuncB"].map((t) => {
        const fn = document.createElementNS(NS, t);
        fn.setAttribute("type", "linear");
        ct.appendChild(fn);
        return fn;
      });
      const mat = document.createElementNS(NS, "feColorMatrix");
      mat.setAttribute("type", "matrix");
      f.appendChild(ct); f.appendChild(mat); defs.appendChild(f);
      gfpFilters[ch] = { funcs, mat };
    });
    document.body.appendChild(svg);
  }

  function applyChannel(ch) {
    ensureFilters();
    const c = state.gfpCh[ch], f = gfpFilters[ch];
    f.funcs.forEach((fn) => {
      fn.setAttribute("slope", c.gain.toFixed(4));
      fn.setAttribute("intercept", (-c.black * c.gain).toFixed(4));
    });
    const [r, g, b] = hexRgb(c.color).map((v) => v / 255);
    f.mat.setAttribute("values",
      `${r} 0 0 0 0  ${g} 0 0 0 0  ${b} 0 0 0 0  0 0 0 1 0`);
    // the ring is drawn in SVG on top of the image, outside the filter, so it is repainted here
    document.querySelectorAll(`.sm-gfp-cell[data-ch="${ch}"] .sm-ring`).forEach((el) => {
      el.setAttribute("stroke", ringFor(c.color));
    });
  }

  const meta = () => state.data.meta;
  const cur = () => state.byId[state.currentId];
  const toPlot = (pUm, zs) => [pUm[0] / XY, pUm[1] / XY, pUm[2] * zs];

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await (await fetch("data/sperm_locations.json")).json(); }
    catch (e) {
      $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_sperm_locations.py</code>.</div></div>`;
      return;
    }
    const m = meta();
    state.data.embryos.forEach((e) => (state.byId[e.id] = e));
    $("#embryo-count").textContent =
      `${m.n} sperm-labelled embryos · ${m.n_zygote} zygote · ${m.n_e2c} early + ${m.n_l2c} late 2-cell` +
      ` · ${m.n_gfp} with GFP frames`;
    // the legend swatch has to track the shared sperm colour, not a hardcoded hex
    $("#sw-sperm").style.background = SPERM;
    fillSorts();
    wire();
    $("#rdrawer").hidden = false; $("#drawer").hidden = false;
    renderList();
    // open on the most extreme example of the default sort — the point of the project
    const first = ranked()[0];
    if (first) select(first.id);
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
  })();

  /** `reset` forces the stage's headline measure. Changing stage must reset even when the old
   *  key still exists in the new list — "cortex" is in both, so without this, switching to
   *  2-cell would silently keep sorting by cortex and bury the junction, which is the whole
   *  reason the 2-cell view exists. */
  function fillSorts(reset) {
    const opts = meta().sorts[state.stage] || [];
    $("#sort-sel").innerHTML = opts.map((o) => `<option value="${o.key}">${o.label}</option>`).join("");
    if (reset || !opts.some((o) => o.key === state.sort)) state.sort = opts[0] ? opts[0].key : "cortex";
    $("#sort-sel").value = state.sort;
    $("#sub-wrap").hidden = state.stage !== "twocell";
  }

  /** The embryos in the current group that actually have the current measure, in order. */
  function ranked() {
    const rows = state.data.embryos.filter((e) => {
      if (e.stage !== state.stage) return false;
      if (state.stage === "twocell" && state.sub !== "all" && e.sub !== state.sub) return false;
      return e.metrics[state.sort] != null;
    }).map((e) => ({ id: e.id, e, v: e.metrics[state.sort] }));
    rows.sort((a, b) => (state.dir === "near" ? a.v - b.v : b.v - a.v));
    return rows;
  }

  function renderList() {
    const rows = ranked(), el = $("#sm-list");
    const unit = "µm";
    if (!rows.length) {
      el.innerHTML = `<div class="sm-empty">No embryo in this group has that measurement.</div>`;
      return;
    }
    el.innerHTML =
      `<div class="sm-head"><span></span><span>embryo</span><span>${unit}</span><span></span></div>` +
      rows.map((r, i) => {
        const hasGfp = !!state.data.gfp[r.id];
        const sub = r.e.stage === "twocell" ? (r.e.sub === "e2c" ? " · early" : r.e.sub === "l2c" ? " · late" : "") : "";
        return `<div class="sm-row${r.id === state.currentId ? " current" : ""}" data-id="${r.id}" ` +
          `title="${r.id}${sub}${hasGfp ? " · has GFP frames" : ""}">` +
          `<span class="n">${i + 1}</span><span class="e">${r.e.label || r.id}${sub}</span>` +
          `<span class="d">${r.v.toFixed(1)}</span>` +
          `<span class="g">${hasGfp ? "●" : ""}</span></div>`;
      }).join("");
    el.querySelectorAll(".sm-row").forEach((n) => n.addEventListener("click", () => select(n.dataset.id)));
  }

  async function select(id) {
    if (!state.byId[id]) return;
    state.currentId = id;
    const e = cur();
    $("#loading").hidden = false; $("#loading-text").textContent = `Loading ${e.label || id}…`;
    try {
      let sc = state.cache[id];
      if (!sc) { sc = await V.loadGz(`data/axes/${id}.json.gz`); state.cache[id] = sc; }
      if (state.currentId !== id) return;
      state.scene = sc;
      $("#controls").hidden = false; $("#placeholder").hidden = true;
      $("#drawer-emb").textContent = `· ${e.label || id}`;
      render3D(); renderReadout(); renderList(); renderActive();
      if (!state.drawerOpen) openDrawer(true);
    } catch (err) {
      $("#placeholder").hidden = false;
      $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Scene missing</div>` +
        `<div class="ph-sub">${err.message || err}</div></div>`;
    } finally { $("#loading").hidden = true; }
  }

  // ───────── 3-D ─────────
  function planeQuad(cUm, nUm, half, zs) {
    const n = nUm, m = Math.hypot(n[0], n[1], n[2]) || 1;
    const nn = [n[0] / m, n[1] / m, n[2] / m];
    const ref = Math.abs(nn[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const cr = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    const u = (v) => { const q = Math.hypot(v[0], v[1], v[2]) || 1; return [v[0] / q, v[1] / q, v[2] / q]; };
    const t = u(cr(nn, ref)), w = u(cr(nn, t));
    const P = [[-1, -1], [1, -1], [1, 1], [-1, 1]].map(([a, b]) =>
      toPlot([cUm[0] + half * (a * t[0] + b * w[0]), cUm[1] + half * (a * t[1] + b * w[1]),
              cUm[2] + half * (a * t[2] + b * w[2])], zs));
    return { type: "mesh3d", x: P.map((p) => p[0]), y: P.map((p) => p[1]), z: P.map((p) => p[2]),
      i: [0, 0], j: [1, 2], k: [2, 3], color: PLANE, opacity: 0.22, name: "Junction plane",
      showlegend: true, hoverinfo: "skip", flatshading: true, legendrank: 41000 };
  }

  function render3D() {
    const s = state.scene, e = cur(); if (!s || !e) return;
    const zs = s.z_scale, traces = V.bodyTraces(s);
    const lm = s.landmarks || {};

    if (state.planeOn && e.junction) {
      traces.push(planeQuad(e.junction.mid_um, e.junction.axis_um, e.junction.sep_um * 0.8, zs));
    }
    if (state.marksOn) {
      const put = (p, name, color, sym, rank) => {
        if (!p) return;
        traces.push({ type: "scatter3d", mode: "markers", name,
          x: [p[0]], y: [p[1]], z: [p[2]],
          marker: { size: 7, color, symbol: sym, line: { width: 1, color: "#fff" } },
          hovertemplate: `${name}<extra></extra>`, legendrank: rank });
      };
      put(lm.polar_plot, "Polar body", MARK, "circle", 40001);
      (lm.nuclei_plots || []).forEach((n, i) => put(n, `Nucleus ${i + 1}`, "#0891b2", "circle", 40010 + i));
      if (e.pron_um) {
        put(toPlot(e.pron_um.maternal, zs), "Maternal pronucleus", "#db2777", "diamond", 40020);
        put(toPlot(e.pron_um.paternal, zs), "Paternal pronucleus", "#2563eb", "diamond", 40021);
      }
    }
    if (state.spermOn && e.sperm_plot) {
      traces.push({ type: "scatter3d", mode: "markers", name: "Sperm",
        x: [e.sperm_plot[0]], y: [e.sperm_plot[1]], z: [e.sperm_plot[2]],
        marker: { size: 11, color: SPERM, symbol: "diamond", line: { width: 2, color: "#fff" } },
        hovertemplate: "labelled sperm<extra></extra>", legendrank: 40000 });
    }
    Plotly.react($("#plot-host"), traces, V.sceneLayout(s.extents, s.id), V.plotConfig);
  }

  const LBL = { cortex: "to the cortex", polar: "to the polar body", maternal: "to the maternal pronucleus",
                paternal: "to the paternal pronucleus", junction: "to the junction", nucleus: "to the nearest nucleus" };

  function renderReadout() {
    const e = cur(); if (!e) return;
    const M = e.metrics;
    const rows = Object.keys(LBL).filter((k) => M[k] != null).map((k) =>
      `<div class="sm-r-line${k === state.sort ? " is-sort" : ""}">` +
      `<span class="k">${LBL[k]}</span><span class="v">${M[k].toFixed(1)} µm</span></div>`).join("");
    const side = e.junction ? `<div class="sm-r-line"><span class="k">blastomere</span>` +
      `<span class="v">${e.junction.side === 0 ? "A" : "B"}</span></div>` : "";
    const g = state.data.gfp[e.id];
    $("#sm-readout").innerHTML =
      `<div class="sm-r-head">${e.label || e.id}</div>` + rows + side +
      (g ? `<div class="sm-r-line"><span class="k">GFP z</span><span class="v">${g.z} / ${g.nz}</span></div>`
         : `<div class="sm-r-line"><span class="k">GFP</span><span class="v">not pointed at</span></div>`);
  }

  // ───────── bottom drawer: the GFP evidence ─────────
  function renderGfp() {
    const wrap = $("#sm-gfp"), e = cur(); if (!e) return;
    const g = state.data.gfp[e.id];
    $("#sm-gfp-sub").textContent = g ? `· ${g.name} · z ${g.z} of ${g.nz}` : "";
    if (!g) {
      wrap.innerHTML = `<div class="sm-gfp-empty">The sperm for this embryo was located on the ` +
        `MERFISH image but never pointed at in the GFP stack, so there is no frame to show.</div>`;
      $("#sm-gfp-note").textContent = "";
      return;
    }
    const dir = meta().gfp_dir;
    const pj = g.proj || null;
    const CH = [
      { ch: 0, name: "488 · GFP", slice: "ch0", mip: "mip0" },
      { ch: 1, name: "405", slice: "ch1", mip: "mip1" },
    ];
    const projCap = pj ? `max-Z · planes ${pj.z0}–${pj.z1}` : "max-Z projection";
    // the marker only means anything on the slice the sperm was called on, not on a projection
    // through the stack — but the x,y is still where it is, so it is drawn on both.
    const sx = g.x != null ? (g.x / g.src_w) * 100 : null;
    const sy = g.y != null ? (g.y / g.src_h) * 100 : null;
    const ring = (ch) => (state.gfpMark && sx != null)
      ? `<svg viewBox="0 0 100 100" preserveAspectRatio="none">
           <circle class="sm-ring" cx="${sx}" cy="${sy}" r="2.2" fill="none"
                   stroke="${ringFor(state.gfpCh[ch].color)}" stroke-width="0.7"
                   vector-effect="non-scaling-stroke"/>
           <line class="sm-ring" x1="${sx - 5}" y1="${sy}" x2="${sx - 3}" y2="${sy}"
                 stroke="${ringFor(state.gfpCh[ch].color)}" stroke-width="0.6" vector-effect="non-scaling-stroke"/>
           <line class="sm-ring" x1="${sx + 3}" y1="${sy}" x2="${sx + 5}" y2="${sy}"
                 stroke="${ringFor(state.gfpCh[ch].color)}" stroke-width="0.6" vector-effect="non-scaling-stroke"/>
         </svg>` : "";

    const tile = (ch, key, cap) => {
      if (!g.files[key]) return "";
      const href = `${dir}/${g.files[key]}`;
      return `<div class="sm-gfp-cell" data-ch="${ch}" data-key="${key}">
        <div class="sm-gfp-cap">${cap}
          <button type="button" class="dl sm-dl1" data-key="${key}" data-ch="${ch}"
                  title="Save this frame at ${g.w}px, exactly as it looks here">↓ full size</button>
        </div>
        <div class="sm-gfp-wrap"><img src="${href}" alt="${cap}" loading="lazy"
             style="filter:url(#smf${ch})">${ring(ch)}</div>
      </div>`;
    };

    const controls = (ch) => {
      const c = state.gfpCh[ch];
      const sw = SWATCHES.map((s) => `<button type="button" class="sm-sw${s.toLowerCase() === c.color.toLowerCase() ? " on" : ""}" ` +
        `data-ch="${ch}" data-color="${s}" style="background:${s}" title="${s}" aria-label="Colour ${s}"></button>`).join("");
      return `<div class="sm-ch-ctl" data-ch="${ch}">
        <span class="sm-cl">colour</span><span class="sm-sws">${sw}</span>
        <label class="sm-sl"><span class="sm-cl">brightness</span>
          <input type="range" class="sm-gain" data-ch="${ch}" min="0.2" max="6" step="0.05" value="${c.gain}">
          <output class="sm-out">${(+c.gain).toFixed(2)}×</output></label>
        <label class="sm-sl"><span class="sm-cl">black point</span>
          <input type="range" class="sm-black" data-ch="${ch}" min="0" max="0.9" step="0.01" value="${c.black}">
          <output class="sm-out">${Math.round(c.black * 100)}%</output></label>
        <button type="button" class="sm-rst" data-ch="${ch}" title="Back to the default colour and brightness">reset</button>
      </div>`;
    };

    wrap.innerHTML = CH.map((c) => `
      <section class="sm-ch" data-ch="${c.ch}">
        <div class="sm-ch-head">${c.name}</div>
        <div class="sm-gfp-grid">${tile(c.ch, c.slice, `z-slice · z ${g.z}`)}${tile(c.ch, c.mip, projCap)}</div>
        ${controls(c.ch)}
      </section>`).join("");
    [0, 1].forEach(applyChannel);

    const dropped = pj && pj.dropped ? pj.dropped : 0;
    $("#sm-gfp-note").innerHTML =
      `Rendered at <b>${g.w}px</b> from the ${g.src_w}px source, in greyscale — the colour and ` +
      `brightness above are applied in your browser, so retuning them costs nothing and the ` +
      `download matches what you see. ` +
      (pj
        ? `The projection is a max over <b>${pj.n} planes (${pj.z0}–${pj.z1} of ${g.nz})</b>, the ` +
          `embryo's own z extent` + (pj.how === "405" ? " as marked by the 405 channel" : "") +
          (dropped ? `, with ${dropped} acquisition-artifact plane${dropped > 1 ? "s" : ""} dropped` : "") +
          `. Projecting the whole stack instead lets a single saturated calibration frame set the ` +
          `contrast window and the sperm disappears. `
        : "") +
      `The ring is the pixel that was clicked.`;
  }

  // ───────── image export ─────────
  /* Archives are built by zip.js (STORE method, shared so the test suite validates the same
   * bytes the browser ships). */
  const makeZip = (files) => window.SurfZip.zipBlob(files);

  const loadImg = (src) => new Promise((res, rej) => {
    const i = new Image();
    i.onload = () => res(i); i.onerror = () => rej(new Error("load " + src));
    i.src = src;
  });

  const saveBlob = (blob, name) => {
    const u = URL.createObjectURL(blob), a = document.createElement("a");
    a.href = u; a.download = name; document.body.appendChild(a); a.click();
    a.remove(); setTimeout(() => URL.revokeObjectURL(u), 4000);
  };

  /* Bake the channel's colour + window/level into pixels. The SVG filter does this for the
   * screen, but a downloaded file has to carry the look with it, so the same arithmetic is
   * repeated here rather than hoping a canvas honours filter:url(). */
  async function renderExport(g, key, ch, withRing) {
    const img = await loadImg(`${meta().gfp_dir}/${g.files[key]}`);
    const cv = document.createElement("canvas");
    cv.width = img.naturalWidth; cv.height = img.naturalHeight;
    const cx = cv.getContext("2d");
    cx.drawImage(img, 0, 0);
    const d = cx.getImageData(0, 0, cv.width, cv.height), p = d.data;
    const c = state.gfpCh[ch], [tr, tg, tb] = hexRgb(c.color);
    for (let i = 0; i < p.length; i += 4) {
      let v = (p[i] / 255 - c.black) * c.gain;
      v = v < 0 ? 0 : v > 1 ? 1 : v;
      p[i] = v * tr; p[i + 1] = v * tg; p[i + 2] = v * tb;
    }
    cx.putImageData(d, 0, 0);
    if (withRing && g.x != null) {
      const s = cv.width / g.src_w, x = g.x * s, y = g.y * s, r = Math.max(10, cv.width * 0.022);
      cx.strokeStyle = ringFor(c.color);
      cx.lineWidth = Math.max(1.5, cv.width / 900);
      cx.beginPath(); cx.arc(x, y, r, 0, Math.PI * 2); cx.stroke();
      cx.beginPath();
      cx.moveTo(x - r * 2.3, y); cx.lineTo(x - r * 1.35, y);
      cx.moveTo(x + r * 1.35, y); cx.lineTo(x + r * 2.3, y);
      cx.stroke();
    }
    return cv;
  }

  const toBytes = (cv, type, q) => new Promise((res) =>
    cv.toBlob((b) => b.arrayBuffer().then((a) => res(new Uint8Array(a))), type, q));

  const KEY_CAP = { ch0: "488_zslice", ch1: "405_zslice", mip0: "488_maxZ", mip1: "405_maxZ" };

  async function embryoFiles(e, { png = true, prefix = "" } = {}) {
    const g = state.data.gfp[e.id];
    if (!g) return [];
    const out = [];
    for (const [key, ch] of [["ch0", 0], ["ch1", 1], ["mip0", 0], ["mip1", 1]]) {
      if (!g.files[key]) continue;
      const cv = await renderExport(g, key, ch, state.gfpMark);
      out.push({ name: `${prefix}${e.id}_${KEY_CAP[key]}.${png ? "png" : "jpg"}`,
                 data: await toBytes(cv, png ? "image/png" : "image/jpeg", 0.92) });
    }
    return out;
  }

  function readmeFor(list) {
    // width comes off an actual entry rather than meta, which does not carry it
    const w = (state.data.gfp[list[0] && list[0].id] || {}).w || 1600;
    const L = [
      "Sperm Location Browser — GFP frames", "",
      `embryos in this archive: ${list.length}`,
      `source frames: ${w}px, rendered from the acquisition stacks`,
      "",
      "Each embryo contributes up to four frames:",
      "  <id>_488_zslice  the 488 (GFP) plane the sperm was called on",
      "  <id>_405_zslice  the 405 plane at the same z",
      "  <id>_488_maxZ    max projection over the embryo's z extent, 488",
      "  <id>_405_maxZ    the same for 405",
      "",
      "The projection covers only the embryo's own z extent, with acquisition-artifact planes",
      "excluded — plane 0 of these stacks saturates and would otherwise set the contrast window.",
      "Colour and brightness are as set in the browser at export time.",
      "", "embryo\tstage\t" + (meta().sorts[state.stage] || []).map((s) => s.key).join("\t"),
    ];
    list.forEach((e) => L.push(e.id + "\t" + (e.label || e.stage) + "\t" +
      (meta().sorts[state.stage] || []).map((s) => (e.metrics[s.key] != null ? e.metrics[s.key] : "")).join("\t")));
    return new TextEncoder().encode(L.join("\n"));
  }

  async function zipCurrent(btn) {
    const e = cur(); if (!e) return;
    const old = btn.textContent;
    btn.disabled = true; btn.textContent = "packing…";
    try {
      const files = await embryoFiles(e, { png: true });
      if (!files.length) { btn.textContent = "no frames"; return; }
      files.push({ name: "README.txt", data: readmeFor([e]) });
      saveBlob(makeZip(files), `sperm_gfp_${e.id}.zip`);
    } catch (err) {
      btn.textContent = "failed"; console.error(err); return;
    } finally {
      setTimeout(() => { btn.disabled = false; btn.textContent = old; }, 1200);
    }
  }

  async function zipAll(btn) {
    const list = ranked().map((r) => r.e).filter((e) => state.data.gfp[e.id]);
    if (!list.length) return;
    const old = btn.textContent;
    btn.disabled = true;
    const files = [];
    try {
      for (let i = 0; i < list.length; i++) {
        btn.textContent = `packing ${i + 1}/${list.length}…`;
        // JPEG for the bulk archive: the same frames as PNG run to a few hundred MB, which the
        // browser has to hold in memory all at once before it can hand over a single blob.
        files.push(...await embryoFiles(list[i], { png: false }));
      }
      files.push({ name: "README.txt", data: readmeFor(list) });
      saveBlob(makeZip(files), `sperm_gfp_${state.stage}_${list.length}embryos.zip`);
    } catch (err) {
      btn.textContent = "failed"; console.error(err); return;
    } finally {
      setTimeout(() => { btn.disabled = false; btn.textContent = old; }, 1500);
    }
  }

  async function saveOne(btn) {
    const e = cur(), g = e && state.data.gfp[e.id];
    if (!g) return;
    const key = btn.dataset.key, ch = +btn.dataset.ch;
    const cv = await renderExport(g, key, ch, state.gfpMark);
    cv.toBlob((b) => saveBlob(b, `${e.id}_${KEY_CAP[key]}.png`), "image/png");
  }

  // ───────── bottom drawer: the distribution ─────────
  const baseLayout = (xt, yt) => ({
    margin: { l: 52, r: 12, t: 28, b: 44 }, showlegend: false,
    paper_bgcolor: "transparent", plot_bgcolor: "#fcfdfe",
    xaxis: { title: { text: xt, font: { size: 11 } }, gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
    yaxis: { title: { text: yt, font: { size: 11 } }, gridcolor: "#eef1f5", zeroline: false, tickfont: { size: 9 } },
    font: { size: 11, color: "#334155" },
  });

  function renderDist() {
    const hist = $("#sm-hist"), cross = $("#sm-cross");
    if (!hist.offsetParent) return;
    const rows = ranked(), e = cur();
    const label = (meta().sorts[state.stage].find((o) => o.key === state.sort) || {}).label || state.sort;
    $("#sm-dist-sub").textContent =
      `· ${rows.length} embryos · ${state.stage === "zygote" ? "zygote" :
        (state.sub === "all" ? "early + late 2-cell" : state.sub === "e2c" ? "early 2-cell" : "late 2-cell")}`;

    // ---- left: the distribution of the sort measure ----
    const vals = rows.map((r) => r.v);
    const mine = e && e.metrics[state.sort];
    const t1 = [{ type: "histogram", x: vals, marker: { color: "#cbd5e1" },
      nbinsx: Math.max(8, Math.min(22, Math.ceil(rows.length / 2))),
      hovertemplate: "%{x:.1f} µm: %{y}<extra></extra>" }];
    const l1 = baseLayout(label + " (µm)", "embryos");
    l1.title = { text: "Distribution", font: { size: 11.5, color: "#64748b" }, x: 0, xanchor: "left" };
    if (mine != null) l1.shapes = [{ type: "line", x0: mine, x1: mine, yref: "paper", y0: 0, y1: 1,
      line: { color: SPERM, width: 2 } }];
    if (mine != null) l1.annotations = [{ x: mine, y: 1, yref: "paper", text: e.label || e.id,
      showarrow: false, font: { size: 9, color: SPERM }, yanchor: "bottom", xanchor: "left" }];
    Plotly.react(hist, t1, l1, { responsive: true, displaylogo: false, displayModeBar: false });

    // ---- right: every sperm in one shared cross-section ----
    // 2-cell: the junction frame — x = signed distance along the junction normal (0 = interface),
    //         y = how far off that axis the sperm sits. Zygotes have no junction, so the frame is
    //         radial instead: x = distance from the centre normalised by the cell radius.
    let t2, l2;
    if (state.stage === "twocell") {
      const pts = rows.filter((r) => r.e.junction);
      const x = pts.map((r) => r.e.junction.signed_um);
      const y = pts.map((r) => Math.max(0, Math.sqrt(Math.max(0,
        Math.pow(r.e.metrics.cortex || 0, 2)))));   // offset proxy: distance to the cortex
      t2 = [
        { type: "scatter", mode: "markers", x, y,
          marker: { size: 9, color: pts.map((r) => (r.e.sub === "e2c" ? "#7c3aed" : "#f97316")),
                    line: { width: 1, color: "#fff" } },
          text: pts.map((r) => r.e.label || r.e.id),
          hovertemplate: "<b>%{text}</b><br>%{x:.1f} µm from the junction<br>%{y:.1f} µm from the cortex<extra></extra>" },
      ];
      if (e && e.junction) t2.push({ type: "scatter", mode: "markers",
        x: [e.junction.signed_um], y: [e.metrics.cortex || 0],
        marker: { size: 15, color: SPERM, line: { width: 2, color: "#fff" } },
        hovertemplate: `${e.label || e.id} · selected<extra></extra>` });
      l2 = baseLayout("signed distance from the junction (µm) — 0 = the interface", "distance to the cortex (µm)");
      l2.shapes = [{ type: "line", x0: 0, x1: 0, yref: "paper", y0: 0, y1: 1,
        line: { color: PLANE, width: 2, dash: "dot" } }];
      l2.title = { text: "Junction frame · purple early, orange late", font: { size: 11.5, color: "#64748b" }, x: 0, xanchor: "left" };
    } else {
      const pts = rows.filter((r) => r.e.metrics.cortex != null && r.e.metrics.polar != null);
      t2 = [{ type: "scatter", mode: "markers",
        x: pts.map((r) => r.e.metrics.cortex), y: pts.map((r) => r.e.metrics.polar),
        marker: { size: 9, color: "#94a3b8", line: { width: 1, color: "#fff" } },
        text: pts.map((r) => r.e.label || r.e.id),
        hovertemplate: "<b>%{text}</b><br>%{x:.1f} µm from the cortex<br>%{y:.1f} µm from the polar body<extra></extra>" }];
      if (e && e.metrics.cortex != null && e.metrics.polar != null) t2.push({ type: "scatter", mode: "markers",
        x: [e.metrics.cortex], y: [e.metrics.polar],
        marker: { size: 15, color: SPERM, line: { width: 2, color: "#fff" } },
        hovertemplate: `${e.label || e.id} · selected<extra></extra>` });
      l2 = baseLayout("distance to the cortex (µm)", "distance to the polar body (µm)");
      l2.title = { text: "Where every sperm sits", font: { size: 11.5, color: "#64748b" }, x: 0, xanchor: "left" };
    }
    Plotly.react(cross, t2, l2, { responsive: true, displaylogo: false, displayModeBar: false });
  }

  // ───────── plumbing ─────────
  const RENDER = { gfp: renderGfp, dist: renderDist };
  function renderActive() { (RENDER[state.tab] || renderGfp)(); }
  function openDrawer(open) {
    state.drawerOpen = open;
    $("#drawer").dataset.open = open ? "true" : "false";
    $("#drawer-handle").setAttribute("aria-expanded", String(open));
    if (open) renderActive();
  }
  function regroup(reset) { fillSorts(reset); renderList(); renderActive();
    const rows = ranked();
    if (rows.length && !rows.some((r) => r.id === state.currentId)) select(rows[0].id);
  }

  function wire() {
    $("#stage-sel").addEventListener("change", (e) => { state.stage = e.target.value; regroup(true); });
    $("#sub-sel").addEventListener("change", (e) => { state.sub = e.target.value; regroup(); });
    $("#sort-sel").addEventListener("change", (e) => { state.sort = e.target.value; renderList(); renderReadout(); renderActive(); });
    $("#sm-dir").addEventListener("click", (e) => {
      const b = e.target.closest("button"); if (!b) return;
      state.dir = b.dataset.dir;
      [...e.currentTarget.children].forEach((x) => x.classList.toggle("active", x === b));
      renderList();
    });
    $("#sperm-show").addEventListener("change", (e) => { state.spermOn = e.target.checked; render3D(); });
    $("#marks-show").addEventListener("change", (e) => { state.marksOn = e.target.checked; render3D(); });
    $("#plane-show").addEventListener("change", (e) => { state.planeOn = e.target.checked; render3D(); });
    $("#gfp-mark").addEventListener("change", (e) => { state.gfpMark = e.target.checked; renderGfp(); });

    // The channel controls live inside the re-rendered grid, so they are handled by delegation
    // on the container rather than rebound on every embryo change.
    const gfp = $("#sm-gfp");
    gfp.addEventListener("click", (ev) => {
      const sw = ev.target.closest(".sm-sw");
      if (sw) {
        const ch = +sw.dataset.ch;
        state.gfpCh[ch].color = sw.dataset.color;
        saveGfpPrefs(); applyChannel(ch);
        gfp.querySelectorAll(`.sm-sw[data-ch="${ch}"]`)
           .forEach((b) => b.classList.toggle("on", b === sw));
        return;
      }
      const rst = ev.target.closest(".sm-rst");
      if (rst) {
        const ch = +rst.dataset.ch;
        state.gfpCh[ch] = { ...GFP_DEF[ch] };
        saveGfpPrefs(); renderGfp();
        return;
      }
      const dl = ev.target.closest(".sm-dl1");
      if (dl) saveOne(dl);
    });
    gfp.addEventListener("input", (ev) => {
      const t = ev.target;
      if (!t.classList.contains("sm-gain") && !t.classList.contains("sm-black")) return;
      const ch = +t.dataset.ch, gain = t.classList.contains("sm-gain");
      state.gfpCh[ch][gain ? "gain" : "black"] = +t.value;
      const out = t.parentElement.querySelector(".sm-out");
      if (out) out.textContent = gain ? (+t.value).toFixed(2) + "×" : Math.round(t.value * 100) + "%";
      applyChannel(ch);
    });
    gfp.addEventListener("change", (ev) => {
      if (ev.target.classList.contains("sm-gain") || ev.target.classList.contains("sm-black")) saveGfpPrefs();
    });
    $("#sm-zip1").addEventListener("click", (ev) => zipCurrent(ev.currentTarget));
    $("#sm-zipall").addEventListener("click", (ev) => zipAll(ev.currentTarget));

    $("#drawer-handle").addEventListener("click", () => openDrawer($("#drawer").dataset.open !== "true"));
    $("#sm-tabs").addEventListener("click", (e) => {
      const b = e.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#sm-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#sm-panels").querySelectorAll(".xs-panel").forEach((p) => (p.hidden = p.dataset.tab !== state.tab));
      renderActive();
    });
    const rz = $("#drawer-resize"); let sh = 0;
    rz.addEventListener("pointerdown", (e) => {
      sh = $("#drawer-body").getBoundingClientRect().height; rz._d = { y: e.clientY };
      rz.setPointerCapture(e.pointerId); e.preventDefault();
    });
    rz.addEventListener("pointermove", (e) => {
      if (!rz._d) return;
      $("#drawer").style.setProperty("--drawer-h",
        Math.max(200, Math.min(window.innerHeight - 120, sh + (rz._d.y - e.clientY))) + "px");
    });
    const end = (e) => { rz._d = null; try { rz.releasePointerCapture(e.pointerId); } catch (_) {} };
    rz.addEventListener("pointerup", end); rz.addEventListener("pointercancel", end);

    const rd = $("#rdrawer"), rh = $("#rdrawer-handle");
    rh.addEventListener("click", () => {
      const o = rd.dataset.open !== "true";
      rd.dataset.open = o ? "true" : "false";
      rh.setAttribute("aria-expanded", String(o));
      if (o) renderList();
    });
    window.addEventListener("resize", () => {
      for (const id of ["#plot-host", "#sm-hist", "#sm-cross"]) { try { Plotly.Plots.resize($(id)); } catch (_) {} }
    });
  }
})();
