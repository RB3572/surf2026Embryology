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
    toImageButtonOptions: { format: "png", scale: 2 },
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

  return { loadGz, buildTabs, markActiveTab, sceneLayout, plotConfig, bodyTraces,
           wireWindow, XY, umToPlot, plotToUm, atlasLink, addWindowExtras };
})();
