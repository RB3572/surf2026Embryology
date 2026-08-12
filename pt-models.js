/* pt-models.js — the pseudotime model family, refittable in the browser.
 *
 * Every model here maps pronuclear geometry to a normalised time. The training rows carry the
 * true t; a fitted model predicts tau.
 *
 * FEATURES ARE IDENTITY-FREE BY CONSTRUCTION. A fixed zygote gives no reliable way to tell the
 * maternal from the paternal pronucleus, so the two distances are only ever used sorted
 * (nearer / farther) or combined (sum / difference). Nothing keyed to male/female is offered.
 *
 * THE LEAKY MODEL CANNOT BE REFIT HERE, and that is deliberate rather than an omission. Its
 * features are the published relative volumes, which are normalised to each pronucleus's own
 * FUTURE endpoint volume — unmeasurable in a snapshot, and inside the training table they encode
 * elapsed time directly. Those columns are simply not in the shipped artifact, so the page can
 * show the number it reaches in the static ranking while making it impossible to reproduce as a
 * deployable model.
 */
window.PTModels = (function () {
  "use strict";

  // ---- features ---------------------------------------------------------------------------
  // A training frame is [t_real, t, d1, d2]; d1/d2 are the raw pair and are immediately sorted.
  function featsFromPair(a, b) {
    const nearer = Math.min(a, b), farther = Math.max(a, b);
    return { nearer, farther, sum: a + b, diff: Math.abs(a - b) };
  }
  const FEATS = ["nearer", "farther", "sum", "diff"];
  const vec = (f, names) => names.map((n) => f[n]);

  // ---- linear algebra: solve (XᵀX + λI) w = Xᵀy by Gaussian elimination with partial pivoting
  function solve(A, b) {
    const n = b.length, M = A.map((r, i) => r.concat([b[i]]));
    for (let c = 0; c < n; c++) {
      let p = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      if (Math.abs(M[p][c]) < 1e-12) return null;
      [M[c], M[p]] = [M[p], M[c]];
      for (let r = 0; r < n; r++) {
        if (r === c) continue;
        const k = M[r][c] / M[c][c];
        for (let j = c; j <= n; j++) M[r][j] -= k * M[c][j];
      }
    }
    return M.map((r, i) => r[n] / r[i]);      // r[i] is the diagonal after elimination
  }

  function fitLinear(X, y, lambda) {
    const p = X[0].length, n = X.length;
    // design matrix with an intercept column; the intercept is never penalised
    const A = Array.from({ length: p + 1 }, () => new Array(p + 1).fill(0));
    const bb = new Array(p + 1).fill(0);
    for (let i = 0; i < n; i++) {
      const xi = [1].concat(X[i]);
      for (let a = 0; a <= p; a++) {
        bb[a] += xi[a] * y[i];
        for (let b2 = 0; b2 <= p; b2++) A[a][b2] += xi[a] * xi[b2];
      }
    }
    for (let a = 1; a <= p; a++) A[a][a] += (lambda || 0);
    const w = solve(A, bb);
    return w ? { w } : null;
  }
  const predLinear = (m, x) => {
    let s = m.w[0];
    for (let i = 0; i < x.length; i++) s += m.w[i + 1] * x[i];
    return s;
  };

  // ---- isotonic: pool adjacent violators ---------------------------------------------------
  /* Sort by x, then walk left to right merging any block whose mean would go DOWN into the block
   * before it. What survives is a staircase: where the data agree the fit follows them, where they
   * contradict each other it flattens and refuses to commit. `sign` is -1 for a decreasing
   * relationship (bigger distance = earlier), which is folded in by fitting on sign*x. */
  function fitIsotonic(x, y, sign) {
    const idx = x.map((v, i) => i).sort((a, b) => (sign * x[a]) - (sign * x[b]));
    const vals = [], wts = [], xs = [];
    for (const i of idx) {
      let v = y[i], w = 1, xv = sign * x[i];
      // merge back while the new block would break monotonicity
      while (vals.length && vals[vals.length - 1] > v) {
        const pv = vals.pop(), pw = wts.pop(); xs.pop();
        v = (v * w + pv * pw) / (w + pw); w += pw;
      }
      vals.push(v); wts.push(w); xs.push(xv);
    }
    // xs holds the LAST x of each block; that is the block's right edge in sign*x
    const knots = xs.map((xv, i) => [sign * xv, vals[i]]);
    knots.sort((a, b) => a[0] - b[0]);
    return { knots, sign };
  }
  function predIsotonic(m, xv) {
    const k = m.knots;
    if (!k.length) return 0;
    if (m.sign < 0) {
      // decreasing in x: find the last knot whose x is <= xv
      if (xv <= k[0][0]) return k[0][1];
      let lo = 0, hi = k.length - 1;
      while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (k[mid][0] <= xv) lo = mid; else hi = mid - 1; }
      return k[lo][1];
    }
    if (xv <= k[0][0]) return k[0][1];
    let lo = 0, hi = k.length - 1;
    while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (k[mid][0] <= xv) lo = mid; else hi = mid - 1; }
    return k[lo][1];
  }

  // ---- gradient-boosted regression trees ---------------------------------------------------
  function buildTree(X, g, rows, depth, minLeaf) {
    const node = {};
    const mean = rows.reduce((s, i) => s + g[i], 0) / rows.length;
    if (depth === 0 || rows.length < 2 * minLeaf) { node.leaf = mean; return node; }
    let best = null;
    for (let f = 0; f < X[0].length; f++) {
      const ord = rows.slice().sort((a, b) => X[a][f] - X[b][f]);
      let sl = 0, nl = 0;
      const tot = ord.reduce((s, i) => s + g[i], 0), N = ord.length;
      for (let i = 0; i < N - 1; i++) {
        sl += g[ord[i]]; nl++;
        if (X[ord[i]][f] === X[ord[i + 1]][f]) continue;
        if (nl < minLeaf || N - nl < minLeaf) continue;
        // variance reduction for squared loss == this gain
        const gain = (sl * sl) / nl + ((tot - sl) * (tot - sl)) / (N - nl) - (tot * tot) / N;
        if (!best || gain > best.gain) {
          best = { gain, f, thr: (X[ord[i]][f] + X[ord[i + 1]][f]) / 2 };
        }
      }
    }
    if (!best || best.gain <= 1e-12) { node.leaf = mean; return node; }
    const L = [], R = [];
    for (const i of rows) (X[i][best.f] <= best.thr ? L : R).push(i);
    if (!L.length || !R.length) { node.leaf = mean; return node; }
    node.f = best.f; node.thr = best.thr;
    node.L = buildTree(X, g, L, depth - 1, minLeaf);
    node.R = buildTree(X, g, R, depth - 1, minLeaf);
    return node;
  }
  const evalTree = (n, x) => (n.leaf !== undefined ? n.leaf : evalTree(x[n.f] <= n.thr ? n.L : n.R, x));

  function fitHGB(X, y, { iters = 120, lr = 0.1, depth = 3, minLeaf = 20 } = {}) {
    const n = X.length, base = y.reduce((s, v) => s + v, 0) / n;
    const pred = new Array(n).fill(base), trees = [];
    const rows = X.map((_, i) => i);
    for (let it = 0; it < iters; it++) {
      const g = new Array(n);
      for (let i = 0; i < n; i++) g[i] = y[i] - pred[i];          // squared-loss residual
      const t = buildTree(X, g, rows, depth, minLeaf);
      trees.push(t);
      for (let i = 0; i < n; i++) pred[i] += lr * evalTree(t, X[i]);
    }
    return { base, lr, trees };
  }
  const predHGB = (m, x) => m.trees.reduce((s, t) => s + m.lr * evalTree(t, x), m.base);

  // ---- the registry -------------------------------------------------------------------------
  const SPECS = {
    isotonic_sum: { label: "Isotonic · distance sum", kind: "isotonic", feat: "sum", sign: -1,
      complexity: "monotone step function",
      note: "assumes the direction of travel, not a constant speed" },
    isotonic_nearer: { label: "Isotonic · nearer distance", kind: "isotonic", feat: "nearer", sign: -1,
      complexity: "monotone step function",
      note: "the same fit on the nearer pronucleus alone" },
    linear_sum: { label: "Linear · distance sum", kind: "linear", feats: ["sum"], lambda: 0,
      complexity: "2 parameters", note: "one-feature straight-line baseline" },
    linear_nearer: { label: "Linear · nearer distance", kind: "linear", feats: ["nearer"], lambda: 0,
      complexity: "2 parameters", note: "transparent one-feature baseline" },
    ridge_symmetric: { label: "Ridge · symmetric distances", kind: "linear",
      feats: ["nearer", "farther", "sum", "diff"], lambda: 1.0, complexity: "4 features, L2",
      note: "sum and diff are exact linear combinations of nearer and farther, so they add no information" },
    ridge_symmetric_core: { label: "Ridge · near + far only", kind: "linear",
      feats: ["nearer", "farther"], lambda: 1.0, complexity: "2 features, L2",
      note: "drops the redundant pair; should match the 4-feature ridge almost exactly" },
    hgb_symmetric: { label: "Gradient boosting · symmetric distances", kind: "hgb",
      feats: ["nearer", "farther", "sum", "diff"], complexity: "depth-3 trees",
      note: "the one nonlinear candidate; justified only if it clearly beats the linear models" },
  };
  const KEYS = Object.keys(SPECS);

  function fit(key, rows) {
    const s = SPECS[key];
    const y = rows.map((r) => r.t);
    if (s.kind === "isotonic") {
      return { key, spec: s, model: fitIsotonic(rows.map((r) => r.f[s.feat]), y, s.sign) };
    }
    const X = rows.map((r) => vec(r.f, s.feats));
    if (s.kind === "linear") return { key, spec: s, model: fitLinear(X, y, s.lambda) };
    return { key, spec: s, model: fitHGB(X, y) };
  }

  function predict(fitted, f) {
    const s = fitted.spec;
    if (!fitted.model) return 0;
    let v;
    if (s.kind === "isotonic") v = predIsotonic(fitted.model, f[s.feat]);
    else if (s.kind === "linear") v = predLinear(fitted.model, vec(f, s.feats));
    else v = predHGB(fitted.model, vec(f, s.feats));
    return Math.max(0, Math.min(1, v));            // tau lives on [0,1]
  }

  // ---- metrics --------------------------------------------------------------------------------
  function rank(v) {
    const o = v.map((_, i) => i).sort((a, b) => v[a] - v[b]);
    const r = new Array(v.length);
    o.forEach((i, p) => (r[i] = p));
    return r;
  }
  function spearman(a, b) {
    const ra = rank(a), rb = rank(b), n = a.length;
    if (n < 3) return NaN;
    const ma = (n - 1) / 2, mb = ma;
    let num = 0, da = 0, db = 0;
    for (let i = 0; i < n; i++) {
      num += (ra[i] - ma) * (rb[i] - mb); da += (ra[i] - ma) ** 2; db += (rb[i] - mb) ** 2;
    }
    return num / Math.sqrt(da * db);
  }
  /** Fraction of pairs whose predicted order matches the true order. Sampled above `cap` pairs so
   *  a big cohort cannot lock the UI; the sample is deterministic. */
  function ordering(t, p, cap = 400000) {
    const n = t.length;
    let ok = 0, tot = 0;
    const total = (n * (n - 1)) / 2;
    const step = total > cap ? Math.ceil(total / cap) : 1;
    let c = 0;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (c++ % step) continue;
        if (t[i] === t[j]) continue;
        tot++;
        if ((t[i] < t[j]) === (p[i] < p[j])) ok++;
        else if (p[i] === p[j]) ok += 0.5;
      }
    }
    return tot ? ok / tot : NaN;
  }

  /** Score predictions on held-out rows. `rows` carry {emb, t}; `pred` is parallel. */
  function score(rows, pred) {
    const n = rows.length;
    const ae = rows.map((r, i) => Math.abs(pred[i] - r.t));
    const byEmb = new Map();
    rows.forEach((r, i) => {
      if (!byEmb.has(r.emb)) byEmb.set(r.emb, []);
      byEmb.get(r.emb).push(i);
    });
    let macro = 0;
    byEmb.forEach((ix) => { macro += ix.reduce((s, i) => s + ae[i], 0) / ix.length; });
    macro /= byEmb.size;
    const sorted = ae.slice().sort((a, b) => a - b);
    const t = rows.map((r) => r.t);
    let within = 0, wn = 0;
    byEmb.forEach((ix) => {
      if (ix.length < 2) return;
      const o = ordering(ix.map((i) => t[i]), ix.map((i) => pred[i]));
      if (!isNaN(o)) { within += o; wn++; }
    });
    return {
      n_frames: n, n_embryos: byEmb.size,
      macro_mae: macro,
      pooled_mae: ae.reduce((s, v) => s + v, 0) / n,
      median_ae: sorted[Math.floor(n / 2)],
      pooled_rmse: Math.sqrt(rows.reduce((s, r, i) => s + (pred[i] - r.t) ** 2, 0) / n),
      pooled_spearman: spearman(t, pred),
      pairwise_ordering_accuracy: ordering(t, pred),
      within_embryo_ordering_accuracy: wn ? within / wn : NaN,
    };
  }

  /** Fit on the embryos NOT in `testIds`, score on those that are. */
  function holdout(key, allRows, testIds) {
    const test = new Set(testIds);
    const tr = allRows.filter((r) => !test.has(r.emb));
    const te = allRows.filter((r) => test.has(r.emb));
    if (!tr.length || !te.length) return null;
    const fitted = fit(key, tr);
    const pred = te.map((r) => predict(fitted, r.f));
    return { fitted, test: te, pred, metrics: score(te, pred),
             n_train_embryos: new Set(tr.map((r) => r.emb)).size,
             n_train_frames: tr.length };
  }

  return { SPECS, KEYS, FEATS, featsFromPair, fit, predict, score, holdout, spearman, ordering };
})();
