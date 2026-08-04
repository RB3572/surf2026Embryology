/* Sperm Alignment (2-cell).
 *
 * A 2-cell embryo gives you one unambiguous axis — the line joining the two blastomere centres.
 * Which blastomere is "first" and which way is up are not in the data, so a pile of 2-cell
 * embryos cannot be overlaid until you pick a rule that supplies them. That rule is the ANCHOR.
 * Everything on this page is downstream of that choice, which is the point: "where does the sperm
 * sit" is a fact about the embryos AND the anchor, not about the embryos alone.
 *
 * Data: data/alignment.json.gz (build_alignment.py) — per embryo the axis frame, a radius map
 * R(t, ψ) for each blastomere (so a cross-section at any azimuth is two slices of a table rather
 * than a mesh cut), ellipsoids for the nuclei and polar body, the sperm, and per gene the counts
 * in each blastomere plus the azimuth of its centroid.
 */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const V = window.VCore;
  const XY = 0.15;
  const SPERM = V.SPERM_COLOR;
  const C_A = "#c98a3f", C_B = "#4a7db5";        // right (alpha) and left (beta) blastomeres
  const C_NUC = "#6b7280", C_PB = "#16a34a";
  const NTHETA = 72;                              // resampling of the mean outline

  const state = {
    data: null, byId: {}, anchor: null, stage: "all", currentId: null,
    scene: null, rank: "side", minSperm: 6, tab: "cross", agreeN: 30,
    outlines: true, mean: true, blobs: true, project: true,
  };

  const meta = () => state.data.meta;
  const POLAR = () => meta().polar_key;
  const cur = () => state.byId[state.currentId];

  // ───────── vector helpers ─────────
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const mul = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];

  /** The frame this anchor puts an embryo in, or null if it cannot orient it.
   *  s = +1 keeps the axis as built (alpha is blastomere b), −1 mirrors it so alpha lands right.
   *  yhat = the anchor's own direction, perpendicular to the axis. */
  function frameFor(e, anchor) {
    let phiDeg, alphaIsB;
    if (anchor === POLAR()) {
      if (!e.polar) return null;
      const d = sub(e.polar, e.mid), along = dot(d, e.u);
      const perp = sub(d, mul(e.u, along));
      const py = dot(perp, e.e2), px = dot(perp, e.e1);
      if (Math.hypot(px, py) < 1e-6) return null;
      phiDeg = Math.atan2(py, px) * 180 / Math.PI;
      alphaIsB = along > 0;                       // the blastomere the polar body sits toward
    } else {
      const g = e.genes[anchor];
      if (!g) return null;
      const [na, nb, phi] = g;
      phiDeg = phi;
      alphaIsB = (nb / Math.max(e.vol_b, 1)) > (na / Math.max(e.vol_a, 1));
    }
    const s = alphaIsB ? 1 : -1;
    const r = phiDeg * Math.PI / 180;
    const yhat = add(mul(e.e1, Math.cos(r)), mul(e.e2, Math.sin(r)));
    return { s, yhat, alphaIsB, phiDeg };
  }

  /** A 3-D point in µm → the anchor frame's 2-D (x = along the axis, y = toward the anchor). */
  const toXY = (e, f, p) => {
    const d = sub(p, e.mid);
    return [dot(d, e.u) * f.s, dot(d, f.yhat)];
  };

  /** Slice a blastomere's radius map at the anchor azimuth → a closed outline in frame coords. */
  function sliceOutline(e, f, which) {
    const m = e.map, n = m.nt, np_ = m.npsi;
    const raw = which === "a" ? m.a : m.b, mx = which === "a" ? m.a_max : m.b_max;
    const bin = (deg) => {
      let d = ((deg + 180) % 360 + 360) % 360;
      return Math.min(np_ - 1, Math.floor(d / 360 * np_));
    };
    const top = bin(f.phiDeg), bot = bin(f.phiDeg + 180);
    const pts = [], back = [];
    for (let i = 0; i < n; i++) {
      const tn = m.t0 + (i + 0.5) / n * (m.t1 - m.t0);
      const x = (tn - 0.5) * m.L * f.s;
      const rt = raw[i * np_ + top] / 255 * mx, rb = raw[i * np_ + bot] / 255 * mx;
      if (rt > 0) pts.push([x, rt]);
      if (rb > 0) back.push([x, -rb]);
    }
    if (pts.length < 3 || back.length < 3) return null;
    back.reverse();
    return pts.concat(back);
  }

  /** r(θ) about a centre, resampled on a fixed θ grid — what makes outlines averageable. */
  function polarProfile(pts, cx) {
    const acc = new Array(NTHETA).fill(0), cnt = new Array(NTHETA).fill(0);
    pts.forEach(([x, y]) => {
      const dx = x - cx, dy = y;
      const th = Math.atan2(dy, dx), r = Math.hypot(dx, dy);
      const k = Math.min(NTHETA - 1, Math.floor(((th + Math.PI) / (2 * Math.PI)) * NTHETA));
      if (r > acc[k]) { acc[k] = r; }
      cnt[k]++;
    });
    // fill gaps by circular interpolation so the mean is defined at every θ
    const have = [];
    for (let k = 0; k < NTHETA; k++) if (acc[k] > 0) have.push(k);
    if (have.length < 6) return null;
    for (let k = 0; k < NTHETA; k++) {
      if (acc[k] > 0) continue;
      let lo = null, hi = null;
      for (let d = 1; d <= NTHETA; d++) {
        const a = (k - d + NTHETA) % NTHETA, b = (k + d) % NTHETA;
        if (lo === null && acc[a] > 0) lo = { k: a, d };
        if (hi === null && acc[b] > 0) hi = { k: b, d };
        if (lo && hi) break;
      }
      acc[k] = (acc[lo.k] * hi.d + acc[hi.k] * lo.d) / (lo.d + hi.d);
    }
    return acc;
  }

  const thetaOf = (k) => -Math.PI + (k + 0.5) / NTHETA * 2 * Math.PI;
  const ringPts = (prof, cx) => {
    const out = prof.map((r, k) => [cx + r * Math.cos(thetaOf(k)), r * Math.sin(thetaOf(k))]);
    out.push(out[0]);
    return out;
  };

  /** Everything the cross-section and the statistics need, for one anchor + stage filter. */
  function assemble(anchor, stage) {
    const rows = [];
    state.data.embryos.forEach((e) => {
      if (stage !== "all" && e.stage !== stage) return;
      const f = frameFor(e, anchor);
      if (!f) return;
      const oa = sliceOutline(e, f, "a"), ob = sliceOutline(e, f, "b");
      if (!oa || !ob) return;
      // centres in frame coords: a sits at −L/2·s, b at +L/2·s
      const xa = -e.map.L / 2 * f.s, xb = e.map.L / 2 * f.s;
      const pa = polarProfile(oa, xa), pb = polarProfile(ob, xb);
      if (!pa || !pb) return;
      const row = { e, f, oa, ob, xa, xb, pa, pb };
      if (e.sperm) {
        const p = toXY(e, f, e.sperm);
        const side = e.sperm_side;                       // 'a' or 'b' from the labelled segment
        const cx = side === "a" ? xa : xb;
        row.sperm = { p, side, cx, theta: Math.atan2(p[1], p[0] - cx),
                      isAlpha: (side === "b") === f.alphaIsB };
      }
      rows.push(row);
    });
    return rows;
  }

  const meanProfile = (rows, key) => {
    if (!rows.length) return null;
    const acc = new Array(NTHETA).fill(0);
    rows.forEach((r) => r[key].forEach((v, k) => (acc[k] += v)));
    return acc.map((v) => v / rows.length);
  };

  /** Circular concentration of the sperm angles, and the alpha/beta split. */
  function stats(rows) {
    const sp = rows.filter((r) => r.sperm);
    const n = sp.length;
    if (!n) return { n: 0 };
    let sx = 0, sy = 0, nAlpha = 0;
    sp.forEach((r) => { sx += Math.cos(r.sperm.theta); sy += Math.sin(r.sperm.theta);
                        if (r.sperm.isAlpha) nAlpha++; });
    const R = Math.hypot(sx / n, sy / n);
    return { n, R, spread: 1 - R, mean: Math.atan2(sy, sx), nAlpha,
             fracAlpha: nAlpha / n, imbalance: Math.abs(nAlpha / n - 0.5) * 2 };
  }

  // ───────── boot ─────────
  (async function init() {
    try { state.data = await V.loadGz("data/alignment.json.gz"); }
    catch (e) {
      $("#placeholder").innerHTML = `<div class="ph-inner"><div class="ph-title">Could not load data</div>` +
        `<div class="ph-sub">Run <code>python3 build_alignment.py</code>.</div></div>`;
      return;
    }
    const m = meta();
    state.data.embryos.forEach((e) => (state.byId[e.id] = e));
    $("#embryo-count").textContent =
      `${m.n_embryos} two-cell embryos · ${m.n_e2c} early + ${m.n_l2c} late · ` +
      `${m.n_sperm} with a labelled sperm · ${m.n_anchors} anchors`;
    state.anchor = POLAR();
    fillAnchors();
    V.buildTabs($("#tabs"), state.data.embryos, selectEmbryo, (e) => ({
      label: V.embryoLabel ? V.embryoLabel(e.id) : e.id,
      sub: e.stage === "e2c" ? "early" : "late",
      title: `${e.id} · ${e.stage === "e2c" ? "early" : "late"} 2-cell` +
             (e.sperm ? " · sperm located" : " · no sperm"),
      cls: e.sperm ? "" : "is-muted",
    }));
    wire();
    $("#drawer").hidden = false; $("#rdrawer").hidden = false;
    openDrawer(true);
    $("#rdrawer").dataset.open = "true";
    $("#rdrawer-handle").setAttribute("aria-expanded", "true");
    renderAll();
    const first = state.data.embryos.find((e) => e.sperm);
    if (first) selectEmbryo(first.id);
  })();

  function fillAnchors() {
    const sel = $("#anchor-select");
    const list = state.data.anchors;
    const opt = (a) => {
      const nm = a.kind === "polar" ? "Polar body" : a.key;
      return `<option value="${a.key}">${nm} · ${a.n_sperm} sperm / ${a.n} embryos</option>`;
    };
    sel.innerHTML = list.map(opt).join("");
    sel.value = state.anchor;
  }

  const anchorName = (k) => (k === POLAR() ? "Polar body" : k);

  // ───────── 3-D: the selected embryo and its sperm ─────────
  async function selectEmbryo(id) {
    if (id === state.currentId) return;
    state.currentId = id;
    V.markActiveTab($("#tabs"), id);
    const e = state.byId[id];
    $("#loading").hidden = false;
    $("#loading-text").textContent = `Loading ${id}…`;
    try {
      const file = (e.stage === "e2c" ? "Early2Cell__" : "Late2Cell__") + id + ".json.gz";
      state.scene = await V.loadGz(`data/segments/${file}`);
    } catch (err) { state.scene = null; }
    $("#loading").hidden = true;
    $("#controls").hidden = false; $("#placeholder").hidden = true;
    render3D();
    renderReadout();
    renderCross();
  }

  function render3D() {
    const e = cur(), sc = state.scene;
    if (!e || !sc) return;
    const traces = V.bodyTraces(sc);
    if (e.sperm) {
      const zs = sc.z_scale || 7;
      traces.push({
        type: "scatter3d", mode: "markers", name: "Sperm",
        x: [e.sperm[0] / XY], y: [e.sperm[1] / XY], z: [e.sperm[2] * zs],
        marker: { size: 11, color: SPERM, symbol: "diamond", line: { width: 2, color: "#fff" } },
        hovertemplate: "sperm<extra></extra>", legendrank: 1,
      });
    }
    const host = $("#plot-host");
    host.innerHTML = "";
    Plotly.newPlot(host, traces, V.sceneLayout(sc.extents, sc.id), V.plotConfig);
  }

  function renderReadout() {
    const e = cur(); if (!e) return;
    const f = frameFor(e, state.anchor);
    const rowsAll = assemble(state.anchor, state.stage);
    const st = stats(rowsAll);
    const L = [];
    L.push(`<div class="al-r-head">${e.id}</div>`);
    L.push(`<div class="al-r-line"><span class="k">stage</span><span class="v">${e.stage === "e2c" ? "early" : "late"} 2-cell</span></div>`);
    L.push(`<div class="al-r-line"><span class="k">blastomere separation</span><span class="v">${e.sep.toFixed(1)} µm</span></div>`);
    if (!f) {
      L.push(`<div class="al-warn">This embryo does not carry <b>${anchorName(state.anchor)}</b>,
        so the anchor cannot orient it — it is left out of the cross-section.</div>`);
    } else {
      L.push(`<div class="al-r-line"><span class="k">right blastomere (α)</span><span class="v">${f.alphaIsB ? "B" : "A"}</span></div>`);
      if (e.sperm) {
        const row = rowsAll.find((r) => r.e.id === e.id);
        L.push(`<div class="al-r-line is-sperm"><span class="k">sperm is in</span><span class="v">${row && row.sperm.isAlpha ? "α (right)" : "β (left)"}</span></div>`);
        if (row) L.push(`<div class="al-r-line"><span class="k">angle on the outline</span><span class="v">${(row.sperm.theta * 180 / Math.PI).toFixed(0)}°</span></div>`);
      } else {
        L.push(`<div class="al-r-line"><span class="k">sperm</span><span class="v">not located</span></div>`);
      }
    }
    if (st.n) {
      L.push(`<div class="al-r-sep">this anchor, over ${st.n} sperm</div>`);
      L.push(`<div class="al-r-line"><span class="k">on the α side</span><span class="v">${st.nAlpha}/${st.n}</span></div>`);
      L.push(`<div class="al-r-line"><span class="k">position spread (1−R)</span><span class="v">${st.spread.toFixed(2)}</span></div>`);
    }
    $("#al-readout").innerHTML = L.join("");
  }

  // ───────── bottom drawer: the aligned cross-section ─────────
  function renderCross() {
    const host = $("#al-cross");
    if (!host || !host.offsetParent) return;
    const rows = assemble(state.anchor, state.stage);
    const sub = $("#al-cross-sub");
    if (!rows.length) {
      Plotly.purge(host); host.innerHTML =
        `<div class="al-empty">No 2-cell embryo carries <b>${anchorName(state.anchor)}</b> in this stage.</div>`;
      if (sub) sub.textContent = "";
      return;
    }
    const nSperm = rows.filter((r) => r.sperm).length;
    const total = state.data.embryos.filter((e) => state.stage === "all" || e.stage === state.stage).length;
    if (sub) sub.textContent = `· ${anchorName(state.anchor)} · orients ${rows.length} of ${total} embryos · ${nSperm} carry a sperm`;

    const traces = [];
    if (state.outlines) rows.forEach((r) => {
      [["oa", r.e.id + " A", C_A], ["ob", r.e.id + " B", C_B]].forEach(([k, nm, col]) => {
        const p = r[k].concat([r[k][0]]);
        traces.push({ type: "scatter", mode: "lines", x: p.map((q) => q[0]), y: p.map((q) => q[1]),
          line: { color: col, width: 1 }, opacity: 0.28, hoverinfo: "skip", showlegend: false });
      });
    });
    // the average outline of each blastomere, and the sperm snapped onto it
    const ma = meanProfile(rows, "pa"), mb = meanProfile(rows, "pb");
    const xa = rows[0].xa < 0 ? -avg(rows.map((r) => Math.abs(r.xa))) : avg(rows.map((r) => Math.abs(r.xa)));
    const xb = -xa;
    if (state.mean) {
      [[ma, xa, "β (left)", C_B], [mb, xb, "α (right)", C_A]].forEach(([prof, cx, nm, col]) => {
        const p = ringPts(prof, cx);
        traces.push({ type: "scatter", mode: "lines", name: `average outline · ${nm}`,
          x: p.map((q) => q[0]), y: p.map((q) => q[1]),
          line: { color: col, width: 4.5, shape: "spline" }, hoverinfo: "skip" });
      });
    }
    if (state.blobs) blobTraces(rows).forEach((t) => traces.push(t));

    const sp = rows.filter((r) => r.sperm);
    if (sp.length) {
      const X = [], Y = [], T = [];
      sp.forEach((r) => {
        const isA = r.sperm.side === "a";
        const prof = isA ? ma : mb, cx = isA ? xa : xb;
        let p = r.sperm.p;
        if (state.project && prof) {
          const th = Math.atan2(p[1], p[0] - (isA ? r.xa : r.xb));
          const k = Math.min(NTHETA - 1, Math.floor(((th + Math.PI) / (2 * Math.PI)) * NTHETA));
          p = [cx + prof[k] * Math.cos(th), prof[k] * Math.sin(th)];
        }
        X.push(p[0]); Y.push(p[1]); T.push(r.e.id);
      });
      traces.push({ type: "scatter", mode: "markers", name: `sperm (${sp.length})`,
        x: X, y: Y, text: T,
        marker: { symbol: "diamond", size: 12, color: SPERM, line: { color: "#fff", width: 1.6 } },
        hovertemplate: "sperm · %{text}<extra></extra>" });
    }
    let lim = 10;
    traces.forEach((t) => (t.x || []).forEach((v, i) => {
      lim = Math.max(lim, Math.abs(v), Math.abs(t.y[i]));
    }));
    lim *= 1.06;
    plotInto(host, traces, {
      dragmode: "pan", margin: { l: 44, r: 10, t: 8, b: 40 }, autosize: true, showlegend: true,
      xaxis: { title: { text: "along the blastomere axis (µm) · + = α, the anchor's blastomere", font: { size: 10 } },
        range: [-lim, lim], scaleanchor: "y", scaleratio: 1, zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
      yaxis: { title: { text: "toward the anchor (µm)", font: { size: 10 } },
        range: [-lim, lim], zeroline: false, gridcolor: "#eef1f5", tickfont: { size: 9 } },
      legend: { orientation: "h", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", font: { size: 9 } },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    });
  }

  const avg = (a) => a.reduce((x, y) => x + y, 0) / (a.length || 1);

  /** Nuclei and polar body as projected ellipses — exact for an ellipsoid in any plane. */
  function blobTraces(rows) {
    const out = [];
    let namedNuc = false, namedPb = false;
    rows.forEach((r) => {
      const e = r.e, f = r.f;
      (e.blobs || []).forEach((b, i) => {
        const isPb = e.blobs.length >= 3 && i === e.blobs.length - 1;
        const c = toXY(e, f, b.c);
        // the 2-D covariance of the projection is P Σ Pᵀ with P the frame's two in-plane axes
        const xh = mul(e.u, f.s), yh = f.yhat;
        const S = [[b.cov[0], b.cov[1], b.cov[2]], [b.cov[1], b.cov[3], b.cov[4]], [b.cov[2], b.cov[4], b.cov[5]]];
        const mv = (M, v) => [dot(M[0], v), dot(M[1], v), dot(M[2], v)];
        const sxx = dot(xh, mv(S, xh)), syy = dot(yh, mv(S, yh)), sxy = dot(xh, mv(S, yh));
        const pts = [];
        for (let k = 0; k <= 40; k++) {
          const t = k / 40 * 2 * Math.PI;
          // 2σ contour of the projected covariance
          const ct = Math.cos(t), st = Math.sin(t);
          const tr = sxx + syy, det = Math.max(sxx * syy - sxy * sxy, 1e-9);
          const l1 = tr / 2 + Math.sqrt(Math.max(tr * tr / 4 - det, 0));
          const l2 = tr / 2 - Math.sqrt(Math.max(tr * tr / 4 - det, 0));
          const ang = 0.5 * Math.atan2(2 * sxy, sxx - syy);
          const a1 = 2 * Math.sqrt(Math.max(l1, 0)), a2 = 2 * Math.sqrt(Math.max(l2, 0));
          pts.push([c[0] + a1 * ct * Math.cos(ang) - a2 * st * Math.sin(ang),
                    c[1] + a1 * ct * Math.sin(ang) + a2 * st * Math.cos(ang)]);
        }
        const showLegend = isPb ? !namedPb : !namedNuc;
        if (isPb) namedPb = true; else namedNuc = true;
        out.push({ type: "scatter", mode: "lines", x: pts.map((q) => q[0]), y: pts.map((q) => q[1]),
          name: isPb ? "polar body" : "nuclei",
          line: { color: isPb ? C_PB : C_NUC, width: isPb ? 1.6 : 1.1, dash: isPb ? "solid" : "dash" },
          opacity: isPb ? 0.7 : 0.5, hoverinfo: "skip", showlegend: showLegend });
      });
    });
    return out;
  }

  // ───────── bottom drawer: which anchors agree ─────────
  /** Per-anchor sperm angles, keyed by embryo — the raw material for both the ranking and the map. */
  function anglesFor(anchor) {
    const out = {};
    state.data.embryos.forEach((e) => {
      if (!e.sperm) return;
      if (state.stage !== "all" && e.stage !== state.stage) return;
      const f = frameFor(e, anchor);
      if (!f) return;
      const p = toXY(e, f, e.sperm);
      const cx = (e.sperm_side === "a" ? -1 : 1) * e.map.L / 2 * f.s;
      out[e.id] = Math.atan2(p[1], p[0] - cx);
    });
    return out;
  }

  function renderAgree() {
    const host = $("#al-agree");
    if (!host || !host.offsetParent) return;
    const N = state.agreeN;
    const ranked = rankedAnchors().slice(0, N);
    if (ranked.length < 2) {
      Plotly.purge(host);
      host.innerHTML = `<div class="al-empty">Not enough anchors clear the minimum of
        ${state.minSperm} sperm to compare.</div>`;
      return;
    }
    const ang = ranked.map((a) => anglesFor(a.key));
    const n = ranked.length;
    const Z = [], TXT = [];
    for (let i = 0; i < n; i++) {
      const row = [], trow = [];
      for (let j = 0; j < n; j++) {
        const A = ang[i], B = ang[j];
        let s = 0, k = 0;
        for (const id in A) if (id in B) { s += Math.cos(A[id] - B[id]); k++; }
        if (k < 3) { row.push(null); trow.push(`${ranked[i].key} vs ${ranked[j].key}<br>only ${k} shared embryos`); }
        else { row.push(s / k); trow.push(`${ranked[i].key} vs ${ranked[j].key}<br>mean cos Δθ = ${(s / k).toFixed(2)} over ${k} embryos`); }
      }
      Z.push(row); TXT.push(trow);
    }
    const labels = ranked.map((a) => anchorName(a.key));
    plotInto(host, [{
      type: "heatmap", z: Z, x: labels, y: labels, text: TXT, hoverinfo: "text",
      colorscale: [[0, "#f7f9fc"], [0.5, "#9dc3e6"], [1, "#0b3d66"]], zmin: -1, zmax: 1,
      colorbar: { title: { text: "mean cos Δθ", font: { size: 9 } }, thickness: 10, tickfont: { size: 8 } },
      hoverongaps: false,
    }], {
      margin: { l: 92, r: 10, t: 8, b: 92 }, autosize: true,
      xaxis: { tickfont: { size: 8 }, tickangle: -90, automargin: false },
      yaxis: { tickfont: { size: 8 }, automargin: false, autorange: "reversed" },
      paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    });
    const sub = $("#al-agree-sub");
    if (sub) sub.textContent = `· ${n} anchors with ≥ ${state.minSperm} sperm`;
  }

  // ───────── right drawer: the ranking ─────────
  function rankedAnchors() {
    const out = [];
    state.data.anchors.forEach((a) => {
      const rows = assemble(a.key, state.stage);
      const st = stats(rows);
      if (!st.n || st.n < state.minSperm) return;
      out.push({ key: a.key, kind: a.kind, n: rows.length, ...st });
    });
    const cmp = {
      side: (x, y) => y.imbalance - x.imbalance || y.n - x.n,
      tight: (x, y) => x.spread - y.spread || y.n - x.n,
      spread: (x, y) => y.spread - x.spread || y.n - x.n,
    }[state.rank];
    out.sort(cmp);
    return out;
  }

  function renderList() {
    const rows = rankedAnchors(), el = $("#al-list");
    if (!rows.length) {
      el.innerHTML = `<div class="al-empty-list">No anchor orients ${state.minSperm} or more
        sperm-carrying embryos in this stage. Lower the minimum.</div>`;
      return;
    }
    const val = (r) => state.rank === "side"
      ? `${r.nAlpha}/${r.n}`
      : r.spread.toFixed(2);
    const head = state.rank === "side" ? "α SIDE" : "SPREAD";
    el.innerHTML =
      `<div class="al-head-row"><span></span><span>anchor</span><span>${head}</span><span>n</span></div>` +
      rows.map((r, i) => `<div class="al-row${r.key === state.anchor ? " current" : ""}" data-key="${r.key}"
          title="${anchorName(r.key)} · orients ${r.n} embryos, ${r.n} with sperm">
        <span class="n">${i + 1}</span>
        <span class="e">${anchorName(r.key)}${r.kind === "polar" ? ' <span class="al-pb">landmark</span>' : ""}</span>
        <span class="d">${val(r)}</span>
        <span class="g">${r.n}</span></div>`).join("");
    el.querySelectorAll(".al-row").forEach((row) => row.addEventListener("click", () => {
      state.anchor = row.dataset.key;
      $("#anchor-select").value = state.anchor;
      renderAll();
    }));
  }

  /* One failing panel must not take the page down with it. Before this, a throw inside the
   * cross-section aborted init before any embryo was selected, so the 3-D scene and the controls
   * never appeared either — the symptom was three broken things, the cause was one. */
  function safely(what, fn) {
    try { fn(); } catch (err) {
      console.error(`[alignment] ${what} failed`, err);
      const host = what === "cross" ? $("#al-cross") : what === "agree" ? $("#al-agree") : null;
      if (host) host.innerHTML = `<div class="al-empty">This panel hit an error: ${err.message}</div>`;
    }
  }

  function renderAll() {
    $("#drawer-anchor").textContent = anchorName(state.anchor);
    safely("readout", renderReadout);
    safely("list", renderList);
    safely(state.tab, () => (state.tab === "cross" ? renderCross() : renderAgree()));
  }

  // ───────── chrome ─────────
  function plotInto(el, traces, layout) {
    el.innerHTML = "";
    el.classList.remove("js-plotly-plot");
    Plotly.newPlot(el, traces, layout,
      { displaylogo: false, responsive: true, scrollZoom: true,
        modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
        toImageButtonOptions: { format: "png", scale: 4 } });
  }

  function openDrawer(open) {
    const d = $("#drawer");
    d.dataset.open = open ? "true" : "false";
    $("#drawer-handle").setAttribute("aria-expanded", String(open));
    if (open) setTimeout(() => (state.tab === "cross" ? renderCross() : renderAgree()), 30);
  }

  function wire() {
    $("#anchor-select").addEventListener("change", (ev) => { state.anchor = ev.target.value; renderAll(); });
    $("#stage-select").addEventListener("change", (ev) => { state.stage = ev.target.value; renderAll(); });
    $("#rank-select").addEventListener("change", (ev) => { state.rank = ev.target.value; renderList(); if (state.tab === "agree") renderAgree(); });
    $("#min-sperm").addEventListener("input", (ev) => {
      state.minSperm = +ev.target.value; $("#min-sperm-val").textContent = ev.target.value;
      renderList(); if (state.tab === "agree") renderAgree();
    });
    $("#agree-n").addEventListener("change", (ev) => { state.agreeN = +ev.target.value; renderAgree(); });
    [["tg-outlines", "outlines"], ["tg-mean", "mean"], ["tg-blobs", "blobs"], ["tg-project", "project"]]
      .forEach(([id, k]) => $("#" + id).addEventListener("change", (ev) => {
        state[k] = ev.target.checked; renderCross();
      }));
    $("#drawer-handle").addEventListener("click", () => openDrawer($("#drawer").dataset.open !== "true"));
    $("#al-tabs").addEventListener("click", (ev) => {
      const b = ev.target.closest(".xs-gtab"); if (!b) return;
      state.tab = b.dataset.tab;
      $("#al-tabs").querySelectorAll(".xs-gtab").forEach((x) => {
        const on = x === b;
        x.classList.toggle("active", on); x.setAttribute("aria-selected", String(on));
      });
      $("#al-panels").querySelectorAll(".xs-panel").forEach((p) => (p.hidden = p.dataset.tab !== state.tab));
      setTimeout(() => (state.tab === "cross" ? renderCross() : renderAgree()), 20);
    });
    const rh = $("#rdrawer-handle");
    rh.addEventListener("click", () => {
      const r = $("#rdrawer"), open = r.dataset.open !== "true";
      r.dataset.open = String(open); rh.setAttribute("aria-expanded", String(open));
    });
    V.wireWindow($("#controls"), $("#controls-header"), [...$("#controls").querySelectorAll(".rz")],
      "alignment_controls_box");
    window.addEventListener("resize", () => {
      try { Plotly.Plots.resize($("#plot-host")); } catch (_) {}
      if (state.tab === "cross") renderCross(); else renderAgree();
    });
    window.addEventListener("vcore:dark", () => render3D());
  }
})();
