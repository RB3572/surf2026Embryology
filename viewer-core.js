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
  function buildTabs(tabsEl, embryos, onSelect, spec) {
    tabsEl.innerHTML = "";
    embryos.forEach((m) => {
      const s = spec(m);
      const b = document.createElement("button");
      b.className = "tab"; b.dataset.id = m.id; b.title = s.title || s.label;
      b.innerHTML = `<span class="tab-label">${s.label}</span>` +
                    `<span class="tab-date">${s.sub || ""}</span>`;
      b.addEventListener("click", () => onSelect(m.id));
      tabsEl.appendChild(b);
    });
  }
  function markActiveTab(tabsEl, id) {
    tabsEl.querySelectorAll(".tab").forEach((b) => b.classList.toggle("active", b.dataset.id === id));
    const a = tabsEl.querySelector(".tab.active");
    if (a) a.scrollIntoView({ inline: "center", block: "nearest", behavior: "smooth" });
  }

  // Plotly 3-D scene layout locked to the embryo extents (axes hidden).
  function sceneLayout(ex, uirev) {
    const pad = (r) => { const p = (r[1] - r[0]) * 0.05 || 1; return [r[0] - p, r[1] + p]; };
    const rx = pad(ex.x), ry = pad(ex.y), rz = pad(ex.z);
    const sx = rx[1] - rx[0], sy = ry[1] - ry[0], sz = rz[1] - rz[0];
    const m = Math.max(sx, sy, sz) || 1;
    return {
      template: "plotly_white", paper_bgcolor: "rgba(0,0,0,0)",
      margin: { l: 0, r: 0, t: 0, b: 0 }, autosize: true, showlegend: true,
      legend: { itemsizing: "constant", font: { size: 12 }, bgcolor: "rgba(255,255,255,0.82)",
                bordercolor: "#e7e9ef", borderwidth: 1, x: 0.99, xanchor: "right", y: 0.98, yanchor: "top" },
      scene: {
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
    for (const lbl of [...scene.mask_labels].sort((a, b) => a - b)) {
      const mesh = scene.region_meshes[String(lbl)];
      if (!mesh) continue;
      const def = (scene.region_defaults || {})[String(lbl)] || { color: "#cccccc", opacity: 0.13 };
      const v = mesh.verts, f = mesh.faces, nV = v.length / 3, nF = f.length / 3;
      const x = new Array(nV), y = new Array(nV), z = new Array(nV);
      for (let i = 0; i < nV; i++) { x[i] = v[i * 3]; y[i] = v[i * 3 + 1]; z[i] = v[i * 3 + 2]; }
      const ii = new Array(nF), jj = new Array(nF), kk = new Array(nF);
      for (let i = 0; i < nF; i++) { ii[i] = f[i * 3]; jj[i] = f[i * 3 + 1]; kk[i] = f[i * 3 + 2]; }
      out.push({ type: "mesh3d", x, y, z, i: ii, j: jj, k: kk, color: def.color,
        opacity: Math.min(def.opacity, 0.13), name: `body M${lbl}`, showlegend: true,
        flatshading: false, hoverinfo: "skip",
        lighting: { ambient: 0.65, diffuse: 0.6, specular: 0.15, roughness: 0.9 }, legendrank: lbl });
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
    range.addEventListener("input", () => { out.textContent = range.value; if (opts.onDotSize) opts.onDotSize(+range.value); });
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
      if (!window.Plotly || !Plotly.downloadImage) return;
      btn.classList.add("busy");
      const fl = gd._fullLayout || {};
      const w = Math.max(fl.width || gd.clientWidth || 900, 320);
      const h = Math.max(fl.height || gd.clientHeight || 600, 240);
      Promise.resolve(Plotly.downloadImage(gd, { format: "png", scale: SCALE, width: w, height: h, filename: filenameFor(gd) }))
        .catch(() => {})
        .then(() => setTimeout(() => btn.classList.remove("busy"), 400));
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

  return { loadGz, buildTabs, markActiveTab, sceneLayout, plotConfig, bodyTraces,
           wireWindow, XY, umToPlot, plotToUm, atlasLink, addWindowExtras, pronMinDist };
})();
