/* viewer-core.js — shared foundation for the analysis-model viewers.
 *
 * The reusable "blank project setup": fetch/gunzip a per-embryo scene, build the
 * embryo nav-bar, wire a draggable + corner-resizable floating window, and draw
 * the base 3-D scene (translucent body meshes + a selected gene's point cloud).
 * A project layers its own analysis on top via the returned hooks.
 */
window.VCore = (function () {
  "use strict";

  // Fetch a .json.gz scene and gunzip it in-browser (handles host auto-decode).
  async function loadGz(url) {
    const resp = await fetch(url);
    if (!resp.ok) throw new Error("HTTP " + resp.status + " for " + url);
    let buf = new Uint8Array(await resp.arrayBuffer());
    if (buf.length >= 2 && buf[0] === 0x1f && buf[1] === 0x8b) {
      const ds = new DecompressionStream("gzip");
      buf = new Uint8Array(await new Response(new Response(buf).body.pipeThrough(ds)).arrayBuffer());
    }
    return JSON.parse(new TextDecoder().decode(buf));
  }

  // Embryo nav tabs. `spec(m)` -> { label, sub, title }.
  // the acquisition year is encoded in the embryo id (YYYYMMDD, after any "Type__" prefix)
  function idYear(id) {
    const m = String(id || "").split("__").pop().match(/^(\d{4})\d{4}/);
    return m ? m[1] : "";
  }
  // canonical display name (TYPE-PROBESET-fovN), matching embryo_naming.py
  const _STAGE_PREFIX = { o: "O", oocyte: "O", z: "Z", zygote: "Z", e2c: "e2c", early: "e2c",
    early2cell: "e2c", l2c: "l2c", late: "l2c", late2cell: "l2c" };
  function embryoLabel(id, stage) {
    const raw = String(id == null ? "" : id).split("__").pop();
    const probe = (window.EMBRYO_PROBESETS || {})[raw];
    if (!probe) return id;
    let m, idStage, fov;
    if ((m = raw.match(/^\d{8}_(oocyte|zygote|e2c|l2c|early2cell|late2cell)_p\d+_(.+)$/i))) { idStage = m[1]; fov = m[2]; }
    else if ((m = raw.match(/^\d{8}_(l2c)_blastomere_p\d+_(.+)$/i))) { idStage = m[1]; fov = m[2]; }
    else if ((m = raw.match(/^\d{8}_sample\d+_(zygote)(\d+(?:_\d+)?)$/i))) { idStage = m[1]; fov = m[2]; }
    else return id;
    const norm = (stage || "").toLowerCase().replace(/[^a-z0-9]/g, "");
    const prefix = _STAGE_PREFIX[norm] || _STAGE_PREFIX[idStage.toLowerCase()];
    return `${prefix}-P${probe}-fov${fov}`;
  }
  function buildTabs(tabsEl, embryos, onSelect, spec) {
    tabsEl.innerHTML = "";
    embryos.forEach((m) => {
      const s = spec(m);
      const b = document.createElement("button");
      // spec.cls lets a project mark tabs that are selectable but degraded — e.g. an embryo
      // with no labelled sperm in a sperm-defined analysis.
      b.className = "tab" + (s.cls ? " " + s.cls : ""); b.dataset.id = m.id; b.title = s.title || s.label;
      const yr = idYear(m.id);
      b.innerHTML = `<span class="tab-label">${s.label}</span>` +
                    `<span class="tab-date">${s.sub || ""}</span>` +
                    (yr ? `<span class="tab-year">${yr}</span>` : "");
      b.addEventListener("click", () => onSelect(m.id));
      tabsEl.appendChild(b);
    });
  }
  function markActiveTab(tabsEl, id) {
    tabsEl.querySelectorAll(".tab").forEach((b) => b.classList.toggle("active", b.dataset.id === id));
    const a = tabsEl.querySelector(".tab.active");
    if (a) a.scrollIntoView({ inline: "center", block: "nearest", behavior: "smooth" });
  }

  // ---- dark "glass render" mode -------------------------------------------------------------
  // A rendering mode for the 3-D scenes only — the page chrome stays light. Cells become
  // translucent white glass envelopes and the smaller nested structures (pronuclei, polar body)
  // glow, on a near-black scene. Ported from the MERFISH atlas viewer, with one change: that
  // viewer classifies regions by voxel count, which we do not ship, so bounding-box volume of the
  // mesh is used instead — same ">=40% of the largest" rule, derived rather than hardcoded, so it
  // works for any embryo whatever its label numbering.
  const DARK_BG = "#04050a";
  const DARK_NUCLEUS_COLORS = ["#22d3ee", "#ff2fd0", "#facc15", "#4ade80",
                               "#fb923c", "#a78bfa", "#f87171", "#38bdf8"];
  const DARK_DOT = 0.5;            // dark render reads best with tiny bright dots
  const DARK_KEY = "surf_dark";
  let darkOn = false;
  try { darkOn = localStorage.getItem(DARK_KEY) === "1"; } catch (_) {}
  const dotSliders = new Set();    // every Dot size control on the page, so all stay in step
  const isDark = () => darkOn;

  function bboxVol(verts) {
    let xn = Infinity, yn = Infinity, zn = Infinity, xx = -Infinity, yx = -Infinity, zx = -Infinity;
    for (let i = 0; i < verts.length; i += 3) {
      const a = verts[i], b = verts[i + 1], c = verts[i + 2];
      if (a < xn) xn = a; if (a > xx) xx = a;
      if (b < yn) yn = b; if (b > yx) yx = b;
      if (c < zn) zn = c; if (c > zx) zx = c;
    }
    if (!isFinite(xn)) return 0;
    return Math.max(xx - xn, 0) * Math.max(yx - yn, 0) * Math.max(zx - zn, 0);
  }

  /** label -> {type:'cell'|'nucleus', color} for a scene, by relative mesh size. */
  function classifyDark(scene) {
    const out = new Map();
    const labels = [...(scene.mask_labels || [])].sort((a, b) => a - b);
    const vol = {};
    let maxV = 0;
    for (const l of labels) {
      const m = (scene.region_meshes || {})[String(l)];
      vol[l] = m && m.verts ? bboxVol(m.verts) : 0;
      if (vol[l] > maxV) maxV = vol[l];
    }
    let ni = 0;
    for (const l of labels) {
      if (maxV <= 0 || vol[l] >= 0.4 * maxV) out.set(l, { type: "cell", color: "#ffffff" });
      else out.set(l, { type: "nucleus", color: DARK_NUCLEUS_COLORS[ni++ % DARK_NUCLEUS_COLORS.length] });
    }
    return out;
  }

  const DARK_CELL = { opacity: 0.10,
    lighting: { ambient: 0.12, diffuse: 0.5, specular: 1.0, roughness: 0.06, fresnel: 2.0 } };
  const DARK_NUC = { opacity: 0.34,
    lighting: { ambient: 0.28, diffuse: 0.65, specular: 1.0, roughness: 0.18 } };
  const DARK_LIGHTPOS = { x: 100, y: 200, z: 300 };

  /* Restyle the 3-D plots already on screen, so toggling is instant and no project needs to
   * listen for anything. bodyTraces/sceneLayout below are dark-aware too, so a later re-render
   * by the project stays consistent. Original light values are stashed on the trace the first
   * time we touch it, which is what lets the toggle go back. */
  function applyDarkToPlot(gd) {
    if (!gd || !gd._fullLayout || !gd.data) return;
    if (!gd._fullLayout.scene) return;            // 3-D scenes only; drawer charts stay light
    // ONLY the segmentation meshes get the glass treatment. Projects also draw analysis meshes
    // — division planes, spheres — and a tilted flat plane has a full-cube bounding box, so
    // sizing would rank it above the cytoplasm and hand the plane the "cell" style while the
    // actual embryo turned into a nucleus. bodyTraces names anatomy "body M<label>"; anything
    // else is an overlay and keeps its own colour, which is right: it is data, not anatomy.
    const ANATOMY = /^body M\d+$/;
    const meshIdx = [], meshVol = [];
    gd.data.forEach((t, i) => {
      if (t.type !== "mesh3d" || !ANATOMY.test(String(t.name || ""))) return;
      meshIdx.push(i);
      const v = [];
      for (let k = 0; k < t.x.length; k++) v.push(t.x[k], t.y[k], t.z[k]);
      meshVol.push(bboxVol(v));
    });
    const maxV = Math.max(0, ...meshVol);
    let ni = 0;
    meshIdx.forEach((i, n) => {
      const t = gd.data[i];
      if (t._vcLight === undefined) {
        t._vcLight = { color: t.color, opacity: t.opacity, lighting: t.lighting,
                       lightposition: t.lightposition };
      }
      if (darkOn) {
        const cell = maxV <= 0 || meshVol[n] >= 0.4 * maxV;
        const st = cell ? DARK_CELL : DARK_NUC;
        t.color = cell ? "#ffffff" : DARK_NUCLEUS_COLORS[ni++ % DARK_NUCLEUS_COLORS.length];
        t.opacity = st.opacity; t.lighting = st.lighting; t.lightposition = DARK_LIGHTPOS;
      } else {
        Object.assign(t, t._vcLight);
      }
    });
    gd.data.forEach((t) => {
      if (t.type !== "scatter3d" || !t.marker) return;
      if (t.marker._vcSize === undefined) {
        t.marker._vcSize = t.marker.size; t.marker._vcOp = t.marker.opacity;
      }
      t.marker.size = darkOn ? DARK_DOT : t.marker._vcSize;
      t.marker.opacity = darkOn ? 1.0 : t.marker._vcOp;
    });
    try {
      Plotly.update(gd, {}, {
        "scene.bgcolor": darkOn ? DARK_BG : "rgba(0,0,0,0)",
        paper_bgcolor: darkOn ? DARK_BG : "rgba(0,0,0,0)",
        template: darkOn ? "plotly_dark" : "plotly_white",
        "legend.bgcolor": darkOn ? "rgba(12,14,20,0.82)" : "rgba(255,255,255,0.82)",
        "legend.bordercolor": darkOn ? "#2a3040" : "#e7e9ef",
        "legend.font.color": darkOn ? "#e6ecff" : "#334155",
      });
    } catch (_) {}
  }

  function applyDark(on, { silent = false } = {}) {
    darkOn = !!on;
    try { localStorage.setItem(DARK_KEY, darkOn ? "1" : "0"); } catch (_) {}
    document.body.classList.toggle("vc-dark", darkOn);
    document.querySelectorAll(".vc-darkbtn").forEach((b) => {
      b.classList.toggle("on", darkOn);
      b.setAttribute("aria-pressed", String(darkOn));
      b.title = darkOn ? "Switch back to the light render" : "Dark glass render of the 3-D scene";
    });
    dotSliders.forEach((s) => s.sync());
    document.querySelectorAll(".js-plotly-plot").forEach(applyDarkToPlot);
    if (!silent) window.dispatchEvent(new CustomEvent("vcore:dark", { detail: { dark: darkOn } }));
  }

  // Navbar toggle, injected into whatever .topbar the page has — every project gets it without
  // touching 21 HTML files.
  (function darkToggle() {
    if (typeof document === "undefined" || window.__vcDark) return;
    window.__vcDark = true;
    const css = document.createElement("style");
    css.textContent =
      ".vc-darkbtn{margin-left:10px;flex:none;display:inline-flex;align-items:center;gap:6px;" +
      "font:600 12px/1 var(--sans,system-ui);color:var(--ink-2,#3c4453);background:var(--surface,#f3f3f1);" +
      "border:1px solid var(--line,#e7e9ef);border-radius:999px;padding:6px 11px;cursor:pointer;" +
      "transition:background .15s,color .15s,border-color .15s}" +
      ".vc-darkbtn:hover{border-color:#cfd6e0}" +
      ".vc-darkbtn.on{background:#0b0d13;border-color:#0b0d13;color:#e6ecff}" +
      ".vc-darkbtn .vc-dd{width:9px;height:9px;border-radius:50%;background:currentColor;opacity:.75}" +
      "body.vc-dark .plot-host,body.vc-dark .stage{background:" + DARK_BG + "}";
    (document.head || document.documentElement).appendChild(css);

    const mount = () => {
      const bar = document.querySelector(".topbar");
      if (!bar || bar.querySelector(".vc-darkbtn")) return;
      const b = document.createElement("button");
      b.type = "button"; b.className = "vc-darkbtn"; b.setAttribute("aria-pressed", "false");
      b.innerHTML = '<span class="vc-dd"></span>Dark render';
      b.addEventListener("click", () => applyDark(!darkOn));
      b.hidden = true;                          // revealed once a 3-D scene actually exists
      bar.appendChild(b);
      applyDark(darkOn, { silent: true });      // reflect the persisted preference
      // Some pages are pure 2-D (grids, charts). Offering a "dark render" there would be a
      // control that does nothing, so the button only appears once a 3-D scene is on the page.
      const has3d = () => [...document.querySelectorAll(".js-plotly-plot")]
        .some((gd) => gd._fullLayout && gd._fullLayout.scene);
      const reveal = () => {
        if (!has3d()) return false;
        b.hidden = false;
        applyDark(darkOn, { silent: true });    // style the scene that just appeared
        return true;
      };
      if (!reveal()) {
        const obs = new MutationObserver(() => { if (reveal()) obs.disconnect(); });
        obs.observe(document.body, { childList: true, subtree: true });
        setTimeout(() => obs.disconnect(), 20000);
      }
    };
    if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", mount);
    else mount();
  })();

  // Plotly 3-D scene layout locked to the embryo extents (axes hidden).
  function sceneLayout(ex, uirev) {
    const pad = (r) => { const p = (r[1] - r[0]) * 0.05 || 1; return [r[0] - p, r[1] + p]; };
    const rx = pad(ex.x), ry = pad(ex.y), rz = pad(ex.z);
    const sx = rx[1] - rx[0], sy = ry[1] - ry[0], sz = rz[1] - rz[0];
    const m = Math.max(sx, sy, sz) || 1;
    const dk = darkOn;
    return {
      template: dk ? "plotly_dark" : "plotly_white",
      paper_bgcolor: dk ? DARK_BG : "rgba(0,0,0,0)",
      margin: { l: 0, r: 0, t: 0, b: 0 }, autosize: true, showlegend: true,
      legend: { itemsizing: "constant",
                font: { size: 12, color: dk ? "#e6ecff" : undefined },
                bgcolor: dk ? "rgba(12,14,20,0.82)" : "rgba(255,255,255,0.82)",
                bordercolor: dk ? "#2a3040" : "#e7e9ef", borderwidth: 1,
                x: 0.99, xanchor: "right", y: 0.98, yanchor: "top" },
      scene: {
        bgcolor: dk ? DARK_BG : undefined,
        xaxis: { visible: false, range: rx, autorange: false },
        yaxis: { visible: false, range: ry, autorange: false },
        zaxis: { visible: false, range: rz, autorange: false },
        aspectmode: "manual", aspectratio: { x: sx / m, y: sy / m, z: sz / m },
        uirevision: uirev, camera: { eye: { x: 1.5, y: 1.5, z: 1.15 } },
      },
    };
  }
  const plotConfig = {
    responsive: true, displaylogo: false,
    modeBarButtonsToRemove: ["tableRotation", "resetCameraLastSave3d", "hoverClosest3d"],
    toImageButtonOptions: { format: "png", scale: 4 },
  };

  // Translucent segmentation-body meshes.
  function bodyTraces(scene) {
    const out = [];
    let dkMap = null;                    // classified lazily, only when the dark render is on
    for (const lbl of [...scene.mask_labels].sort((a, b) => a - b)) {
      const mesh = scene.region_meshes[String(lbl)];
      if (!mesh) continue;
      const def = (scene.region_defaults || {})[String(lbl)] || { color: "#cccccc", opacity: 0.13 };
      const v = mesh.verts, f = mesh.faces, nV = v.length / 3, nF = f.length / 3;
      const x = new Array(nV), y = new Array(nV), z = new Array(nV);
      for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
      const ii = new Array(nF), jj = new Array(nF), kk = new Array(nF);
      for (let i = 0; i < nF; i++) { ii[i] = f[i * 3]; jj[i] = f[i * 3 + 1]; kk[i] = f[i * 3 + 2]; }
      const dc = darkOn ? (dkMap || (dkMap = classifyDark(scene))).get(lbl) : null;
      const st = dc ? (dc.type === "cell" ? DARK_CELL : DARK_NUC) : null;
      out.push({ type: "mesh3d", x, y, z, i: ii, j: jj, k: kk,
        color: dc ? dc.color : def.color,
        opacity: st ? st.opacity : Math.min(def.opacity, 0.13),
        name: `body M${lbl}`, showlegend: true,
        flatshading: false, hoverinfo: "skip",
        lighting: st ? st.lighting
                     : { ambient: 0.65, diffuse: 0.6, specular: 0.15, roughness: 0.9 },
        lightposition: st ? DARK_LIGHTPOS : undefined,
        legendrank: lbl });
    }
    return out;
  }

  // Draggable (header) + 4-corner-resizable floating window, position/size persisted.
  function wireWindow(win, header, corners, storeKey) {
    const stageBox = () => win.parentElement.getBoundingClientRect();
    const save = () => { try { localStorage.setItem(storeKey, JSON.stringify({
      left: parseFloat(win.style.left), top: parseFloat(win.style.top),
      width: parseFloat(win.style.width), height: parseFloat(win.style.height) })); } catch (_) {} };
    const load = () => { let b; try { b = JSON.parse(localStorage.getItem(storeKey) || "null"); } catch (_) {}
      if (!b) return;
      if (isFinite(b.width)) win.style.width = b.width + "px";
      if (isFinite(b.height)) win.style.height = b.height + "px";
      if (isFinite(b.left)) { win.style.left = b.left + "px"; win.style.right = "auto"; }
      if (isFinite(b.top)) win.style.top = b.top + "px"; };
    let onResize = null;
    (function () {   // drag
      let st = null;
      header.addEventListener("pointerdown", (e) => {
        const r = win.getBoundingClientRect(), s = stageBox();
        st = { x: e.clientX, y: e.clientY, left: r.left - s.left, top: r.top - s.top };
        header.setPointerCapture(e.pointerId); e.preventDefault();
      });
      header.addEventListener("pointermove", (e) => { if (!st) return;
        const s = stageBox();
        let l = st.left + (e.clientX - st.x), t = st.top + (e.clientY - st.y);
        l = Math.max(0, Math.min(l, s.width - win.offsetWidth)); t = Math.max(0, Math.min(t, s.height - 30));
        win.style.left = l + "px"; win.style.top = t + "px"; win.style.right = "auto"; });
      const end = (e) => { if (st) { st = null; try { header.releasePointerCapture(e.pointerId); } catch (_) {} save(); } };
      header.addEventListener("pointerup", end); header.addEventListener("pointercancel", end);
    })();
    const MINW = 260, MINH = 140;
    corners.forEach((h) => {
      const cfg = { nw: [1, 1], ne: [0, 1], sw: [1, 0], se: [0, 0] }[h.dataset.corner];
      let st = null;
      h.addEventListener("pointerdown", (e) => {
        const r = win.getBoundingClientRect(), s = stageBox();
        const left = r.left - s.left, top = r.top - s.top;
        st = { x: e.clientX, y: e.clientY, w: r.width, h: r.height, left, top, right: left + r.width, bottom: top + r.height };
        h.setPointerCapture(e.pointerId); e.preventDefault(); e.stopPropagation();
      });
      h.addEventListener("pointermove", (e) => { if (!st) return;
        const dx = e.clientX - st.x, dy = e.clientY - st.y;
        let w, hh, left = st.left, top = st.top;
        if (cfg[0]) { w = Math.max(MINW, st.w - dx); left = st.right - w; } else { w = Math.max(MINW, st.w + dx); }
        if (cfg[1]) { hh = Math.max(MINH, st.h - dy); top = st.bottom - hh; } else { hh = Math.max(MINH, st.h + dy); }
        win.style.width = w + "px"; win.style.height = hh + "px";
        win.style.left = left + "px"; win.style.top = top + "px"; win.style.right = "auto";
        if (onResize) onResize(); });
      const end = (e) => { if (st) { st = null; try { h.releasePointerCapture(e.pointerId); } catch (_) {} save(); if (onResize) onResize(); } };
      h.addEventListener("pointerup", end); h.addEventListener("pointercancel", end);
    });
    load();
    return { save, load, setResizeCb: (cb) => { onResize = cb; } };
  }

  // Physical µm ↔ plot-space (x_px, y_px, frame·z_scale) helpers.
  const XY = 0.15;
  const umToPlot = (p, zs) => [p[0] / XY, p[1] / XY, p[2] * zs];
  const plotToUm = (p, zs) => [p[0] * XY, p[1] * XY, p[2] / zs];

  // Deep-link to the partner MERFISH atlas (merfish.rishib.com) for one embryo id, e.g.
  // "20260425_zygote_p2_2" -> https://merfish.rishib.com/?embryo=Zygote%2F20260425_zygote_p2_2.html
  // Robust minimum distance between the two pronuclei, computed from the DISPLAYED mesh vertices
  // (marching-cubes surfaces) rather than the raw label voxels — so the drawn line always TOUCHES
  // both pronuclei. Returns { line:[a,b] in plot space, distUm } or null. Vert counts are small
  // (≤~1250 per pronucleus) so the exact closest-vertex-pair brute force is fast.
  function pronMinDist(scene) {
    const pl = scene && scene.pron_labels;
    if (!pl || pl.length < 2 || !scene.region_meshes) return null;
    const ma = scene.region_meshes[String(pl[0])], mb = scene.region_meshes[String(pl[1])];
    if (!ma || !mb || !ma.verts || !mb.verts || !ma.verts.length || !mb.verts.length) return null;
    const va = ma.verts, vb = mb.verts;
    let best = Infinity, ai = 0, bi = 0;
    for (let i = 0; i < va.length; i += 3) {
      const ax = va[i], ay = va[i + 1], az = va[i + 2];
      for (let j = 0; j < vb.length; j += 3) {
        const dx = ax - vb[j], dy = ay - vb[j + 1], dz = az - vb[j + 2];
        const d = dx * dx + dy * dy + dz * dz;
        if (d < best) { best = d; ai = i; bi = j; }
      }
    }
    const a = [va[ai], va[ai + 1], va[ai + 2]], b = [vb[bi], vb[bi + 1], vb[bi + 2]];
    const zs = scene.z_scale || 7;
    // plot → µm: x,y × XY (0.15 µm/px); z is stored as frame·z_scale and 1 frame = 1 µm ⇒ µm = plot / z_scale
    const um = (p) => [p[0] * XY, p[1] * XY, p[2] / zs];
    const ua = um(a), ub = um(b);
    const distUm = Math.round(Math.hypot(ua[0] - ub[0], ua[1] - ub[1], ua[2] - ub[2]) * 100) / 100;
    return { line: [a, b], distUm };
  }

  function atlasLink(id) {
    const stage = /oocyte/i.test(id) ? "Oocyte"
      : /zygote/i.test(id) ? "Zygote"
        : /e2c/i.test(id) ? "Early2Cell"
          : /l2c/i.test(id) ? "Late2Cell"
            : "Zygote";
    return "https://merfish.rishib.com/?embryo=" + encodeURIComponent(stage + "/" + id + ".html");
  }

  // Append the shared "Dot size + Atlas link" row to a floating control window's body.
  // opts.onDotSize(size) fires on change; returns { size(), setAtlas(id) }.
  function addWindowExtras(body, opts) {
    opts = opts || {};
    const size = opts.defaultSize == null ? 1.5 : opts.defaultSize;
    const row = document.createElement("div");
    row.className = "controls-row vc-extras";
    row.innerHTML =
      '<label class="ctl vc-dotsize"><span class="ctl-label">Dot size</span>' +
      '<input type="range" min="0.5" max="6" step="0.5" value="' + size + '"><output>' + size + '</output></label>' +
      '<a class="vc-atlas" href="https://merfish.rishib.com" target="_blank" rel="noopener"' +
      ' title="Open this embryo in the MERFISH atlas">Atlas ↗</a>';
    body.appendChild(row);
    const range = row.querySelector("input"), out = row.querySelector("output"), link = row.querySelector("a");
    // Light and dark keep independent dot sizes — the glass render reads best with tiny dots, and
    // it would be tedious to re-set the slider every time you toggle.
    const sizes = { light: size, dark: DARK_DOT };
    const put = (v) => { range.value = v; out.textContent = v; };
    range.addEventListener("input", () => {
      out.textContent = range.value;
      sizes[darkOn ? "dark" : "light"] = +range.value;
      if (opts.onDotSize) opts.onDotSize(+range.value);
    });
    const entry = { sync: () => { put(sizes[darkOn ? "dark" : "light"]);
                                  if (opts.onDotSize) opts.onDotSize(+range.value); } };
    dotSliders.add(entry);
    put(sizes[darkOn ? "dark" : "light"]);      // open in the right mode's size
    return { size: () => +range.value, setAtlas: (id) => { if (id) link.href = atlasLink(id); } };
  }

  // ---- global: a high-resolution PNG download button on EVERY Plotly figure ----
  // Every analysis page loads viewer-core.js, so this one block gives all figures — the 3-D scenes
  // and every 2-D drawer chart alike — a one-click high-res (4×) export, with no per-page wiring.
  // A MutationObserver catches plots as they are created/re-created (Plotly adds `.js-plotly-plot`).
  (function figureDownloads() {
    if (typeof document === "undefined" || window.__vcFigDl) return;
    window.__vcFigDl = true;
    const SCALE = 4;                                       // 4× the on-screen pixels → crisp/print res
    const ICON = '<svg viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor"' +
      ' stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">' +
      '<path d="M12 3v11"/><path d="M7.5 10 12 14.5 16.5 10"/><path d="M5 20h14"/></svg>';
    const style = document.createElement("style");
    style.textContent =
      ".vc-figdl{position:absolute;top:8px;left:8px;z-index:6;display:inline-flex;align-items:center;" +
      "justify-content:center;width:27px;height:27px;padding:0;border-radius:7px;border:1px solid rgba(20,25,35,.14);" +
      "background:rgba(255,255,255,.8);color:#334155;cursor:pointer;opacity:.5;line-height:0;" +
      "transition:opacity .12s,background .12s,box-shadow .12s;-webkit-backdrop-filter:blur(3px);backdrop-filter:blur(3px)}" +
      ".vc-figdl:hover{opacity:1;background:#fff;box-shadow:0 2px 8px rgba(20,25,35,.16)}" +
      ".vc-figdl:active{transform:translateY(1px)}" +
      ".js-plotly-plot:hover>.vc-figdl{opacity:.92}" +
      ".vc-figdl.busy{opacity:.7;cursor:progress}" +
      "@media print{.vc-figdl{display:none}}";
    (document.head || document.documentElement).appendChild(style);

    const slug = (s) => (s || "figure").toLowerCase().replace(/[^a-z0-9]+/g, "-")
      .replace(/^-+|-+$/g, "").slice(0, 60) || "figure";
    const filenameFor = (gd) =>
      slug((document.title || "figure").split("·")[0]) + "-" + slug(gd.id || "figure");

    function download(gd, btn) {
      if (!window.Plotly) return;
      btn.classList.add("busy");
      const fl = gd._fullLayout || {};
      const w = Math.max(fl.width || gd.clientWidth || 900, 320);
      const h = Math.max(fl.height || gd.clientHeight || 600, 240);
      const done = () => setTimeout(() => btn.classList.remove("busy"), 400);
      const name = filenameFor(gd);
      // The pages draw on a TRANSPARENT background so plots sit on the page colour. A transparent
      // PNG then renders black in most slide/doc tools, so the EXPORT is forced onto white. Done on
      // a cloned figure via toImage, so the on-screen plot never flickers or changes.
      const layout = Object.assign({}, gd.layout,
        { paper_bgcolor: "#ffffff", plot_bgcolor: "#ffffff" });
      const legacy = () => Plotly.downloadImage(gd,
        { format: "png", scale: SCALE, width: w, height: h, filename: name });
      if (!Plotly.toImage) { Promise.resolve(legacy()).catch(() => {}).then(done); return; }
      Plotly.toImage({ data: gd.data, layout, config: { displayModeBar: false } },
                     { format: "png", scale: SCALE, width: w, height: h })
        .then((url) => {
          const a = document.createElement("a");
          a.href = url; a.download = name + ".png";
          document.body.appendChild(a); a.click(); a.remove();
        })
        .catch(() => legacy())
        .catch(() => {})
        .then(done);
    }
    function addBtn(gd) {
      if (!gd || !gd.classList || !gd.classList.contains("js-plotly-plot")) return;
      if (gd.querySelector(":scope > .vc-figdl")) return;   // already has one (survives re-render/purge)
      if (getComputedStyle(gd).position === "static") gd.style.position = "relative";
      const btn = document.createElement("button");
      btn.type = "button"; btn.className = "vc-figdl";
      btn.title = "Download this figure — high-resolution PNG";
      btn.setAttribute("aria-label", "Download this figure as a high-resolution PNG");
      btn.innerHTML = ICON;
      btn.addEventListener("click", (e) => { e.preventDefault(); e.stopPropagation(); download(gd, btn); });
      gd.appendChild(btn);
    }
    const scan = (root) => {
      if (!root || root.nodeType !== 1) return;
      if (root.classList && root.classList.contains("js-plotly-plot")) addBtn(root);
      if (root.querySelectorAll) root.querySelectorAll(".js-plotly-plot").forEach(addBtn);
    };
    function start() {
      scan(document.body);
      new MutationObserver((muts) => {
        for (const m of muts) {
          if (m.type === "attributes") {
            if (m.target.classList && m.target.classList.contains("js-plotly-plot")) addBtn(m.target);
          } else if (m.addedNodes) {
            m.addedNodes.forEach(scan);
            if (m.target && m.target.classList && m.target.classList.contains("js-plotly-plot")) addBtn(m.target);
          }
        }
      }).observe(document.body, { childList: true, subtree: true, attributes: true, attributeFilter: ["class"] });
    }
    if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
    else start();
  })();

  return { isDark, applyDark, classifyDark, DARK_BG,
           loadGz, buildTabs, markActiveTab, sceneLayout, plotConfig, bodyTraces,
           wireWindow, XY, umToPlot, plotToUm, atlasLink, addWindowExtras, pronMinDist,
           embryoLabel, idYear };
})();
