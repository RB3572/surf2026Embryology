/* Shared regression-model + correlation-stats toolkit (window.Regressions).
 * Extracted from the Pronuclei / Extended-Pseudotime scatter code so other projects
 * (e.g. Gene Diffusion Rates) can offer the SAME family of fits. Each model returns a
 * predictor y(x) + a fit statistic on its natural scale; statsHtml renders the shared
 * "model · equation · R² · Pearson r · p" line. Pure functions — no project state. */
window.Regressions = (() => {
  const avg = (a) => a.reduce((s, v) => s + v, 0) / (a.length || 1);
  const clmp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
  const sci = (v) => { const a = Math.abs(v); if (!isFinite(v)) return "–"; if (a === 0) return "0";
    return (a >= 1e4 || a < 1e-3) ? v.toExponential(2) : String(+v.toPrecision(3)); };

  // ---- correlation stats (Pearson r + exact t-test p + label-permutation null) ----
  function linreg(xs, ys) {
    const n = xs.length; let sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
    for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; syy += ys[i] * ys[i]; }
    const den = n * sxx - sx * sx;
    const b = den ? (n * sxy - sx * sy) / den : 0;
    const a = (sy - b * sx) / n;
    const rden = Math.sqrt(den * (n * syy - sy * sy));
    const r = rden ? (n * sxy - sx * sy) / rden : 0;
    return { a, b, r, r2: r * r, n };
  }
  function gammaln(x) {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
               -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015; for (let j = 0; j < 6; j++) { y++; ser += c[j] / y; }
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }
  function betacf(a, b, x) {
    const MAXIT = 200, EPS = 3e-12, FPMIN = 1e-300;
    const qab = a + b, qap = a + 1, qam = a - 1; let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d;
      const del = d * c; h *= del; if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }
  function betai(a, b, x) {
    if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(gammaln(a + b) - gammaln(a) - gammaln(b) + a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
  }
  function pearsonP(r, n) {                        // exact two-sided p, H0: ρ = 0 (Student-t, df = n−2)
    if (n < 3) return NaN;
    if (Math.abs(r) >= 1) return 0;
    const df = n - 2, t2 = r * r * df / (1 - r * r);
    return betai(0.5 * df, 0.5, df / (df + t2));
  }
  const pearsonR = (xs, ys) => linreg(xs, ys).r;
  function mulberry32(a) { return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
  function permP(xs, ys, robs, B) {               // label-permutation null (deterministic seed)
    const n = xs.length; if (n < 3) return NaN;
    const a = Math.abs(robs), y = ys.slice(), rnd = mulberry32(0x5eed);
    let count = 0;
    for (let b = 0; b < B; b++) {
      for (let i = n - 1; i > 0; i--) { const j = (rnd() * (i + 1)) | 0; const t = y[i]; y[i] = y[j]; y[j] = t; }
      if (Math.abs(pearsonR(xs, y)) >= a - 1e-12) count++;
    }
    return (count + 1) / (B + 1);
  }
  function fmtP(p) {
    if (p == null || !isFinite(p)) return "–";
    if (p < 1e-4) return p.toExponential(1);
    if (p < 0.1) return p.toPrecision(2);
    return p.toFixed(2);
  }
  const pStars = (p) => (!isFinite(p) ? "" : p < 1e-3 ? " ***" : p < 0.01 ? " **" : p < 0.05 ? " *" : " ns");

  // ---- regression models ----
  function wls(xs, ys, ws) {                      // weighted least squares: y = b0 + b1·x
    let sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
    for (let i = 0; i < xs.length; i++) { const w = ws ? ws[i] : 1;
      sw += w; swx += w * xs[i]; swy += w * ys[i]; swxx += w * xs[i] * xs[i]; swxy += w * xs[i] * ys[i]; }
    const den = sw * swxx - swx * swx, b1 = den ? (sw * swxy - swx * swy) / den : 0;
    return [(swy - b1 * swx) / sw, b1];
  }
  function r2on(xs, ys, pred) {
    const yb = avg(ys); let sr = 0, st = 0;
    for (let i = 0; i < ys.length; i++) { const e = ys[i] - pred(xs[i]); sr += e * e; st += (ys[i] - yb) ** 2; }
    return st > 0 ? 1 - sr / st : 0;
  }
  const r2lin = (xs, ys) => { const [a, b] = wls(xs, ys); return r2on(xs, ys, (x) => a + b * x); };
  function solve3(A, d) {
    A = A.map((r, i) => r.concat(d[i]));
    for (let c = 0; c < 3; c++) { let p = c;
      for (let r = c + 1; r < 3; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
      [A[c], A[p]] = [A[p], A[c]]; if (Math.abs(A[c][c]) < 1e-12) A[c][c] = 1e-12;
      for (let r = 0; r < 3; r++) { if (r === c) continue; const f = A[r][c] / A[c][c]; for (let k = c; k < 4; k++) A[r][k] -= f * A[c][k]; } }
    return [A[0][3] / A[0][0], A[1][3] / A[1][1], A[2][3] / A[2][2]];
  }
  function quadFit(xs, ys) {
    const S = [0, 0, 0, 0, 0], T = [0, 0, 0];
    for (let i = 0; i < xs.length; i++) { let xp = 1; for (let k = 0; k < 5; k++) { S[k] += xp; if (k < 3) T[k] += xp * ys[i]; xp *= xs[i]; } }
    return solve3([[S[0], S[1], S[2]], [S[1], S[2], S[3]], [S[2], S[3], S[4]]], T);
  }
  function irls(xs, ys, fam, init, iter = 60) {   // GLM via iteratively-reweighted least squares
    let [b0, b1] = init;
    for (let k = 0; k < iter; k++) {
      const z = [], w = [];
      for (let i = 0; i < xs.length; i++) { const e = b0 + b1 * xs[i], mu = fam.li(e), g = fam.me(e, mu), V = fam.v(mu);
        w.push(g * g / Math.max(V, 1e-9)); z.push(e + (ys[i] - mu) / (g || 1e-9)); }
      const [n0, n1] = wls(xs, z, w);
      if (!isFinite(n0) || !isFinite(n1)) break;
      if (Math.abs(n0 - b0) + Math.abs(n1 - b1) < 1e-10) { b0 = n0; b1 = n1; break; }
      b0 = n0; b1 = n1;
    }
    return [b0, b1];
  }
  const POIS = { li: (e) => Math.exp(Math.min(e, 30)), me: (e, m) => m, v: (m) => Math.max(m, 1e-9) };
  const NBfam = (t) => ({ li: (e) => Math.exp(Math.min(e, 30)), me: (e, m) => m, v: (m) => m + m * m / t });
  const BINfam = (N) => ({ li: (e) => N / (1 + Math.exp(-clmp(e, -30, 30))), me: (e, m) => { const p = m / N; return N * p * (1 - p); }, v: (m) => { const p = m / N; return Math.max(N * p * (1 - p), 1e-9); } });
  const logInit = (xs, ys) => wls(xs, ys.map((y) => Math.log(Math.max(y, 1))));
  function devR2(xs, ys, mu) {
    const yb = avg(ys), d = (y, m) => 2 * ((y > 0 ? y * Math.log(y / m) : 0) - (y - m));
    let dr = 0, dn = 0; for (let i = 0; i < ys.length; i++) { dr += d(ys[i], mu[i]); dn += d(ys[i], yb); }
    return dn > 0 ? 1 - dr / dn : 0;
  }
  function binDevR2(ys, mu, N) {
    const yb = avg(ys);
    const d = (y, m) => 2 * ((y > 0 ? y * Math.log(y / m) : 0) + (N - y > 0 ? (N - y) * Math.log((N - y) / (N - m)) : 0));
    let dr = 0, dn = 0; for (let i = 0; i < ys.length; i++) { dr += d(ys[i], mu[i]); dn += d(ys[i], yb); }
    return dn > 0 ? 1 - dr / dn : 0;
  }
  function loessPredictor(xs, ys, span) {
    const n = xs.length, k = Math.max(3, Math.round(span * n));
    return (x0) => {
      const sorted = xs.map((x) => Math.abs(x - x0)).sort((a, b) => a - b);
      const h = sorted[Math.min(n - 1, k - 1)] || 1e-9;
      const wx = [], wy = [], ww = [];
      for (let i = 0; i < n; i++) { const u = Math.abs(xs[i] - x0) / h; if (u >= 1) continue; wx.push(xs[i]); wy.push(ys[i]); ww.push((1 - u ** 3) ** 3); }
      if (wx.length < 2) return avg(ys);
      const [a, b] = wls(wx, wy, ww); return a + b * x0;
    };
  }
  const MODELS = {
    linear: { label: "Linear", scale: "raw", bio: "a straight line — constant rate of change; the simplest baseline.",
      fit(xs, ys) { const [a, b] = wls(xs, ys), p = (x) => a + b * x;
        return { predict: p, r2: r2on(xs, ys, p), params: `y = ${sci(a)} ${b >= 0 ? "+" : "−"} ${sci(Math.abs(b))}·x` }; } },
    quadratic: { label: "Quadratic", scale: "raw", bio: "a parabola — one peak or trough.",
      fit(xs, ys) { const [c0, c1, c2] = quadFit(xs, ys), p = (x) => c0 + c1 * x + c2 * x * x;
        return { predict: p, r2: r2on(xs, ys, p), params: `y = ${sci(c0)} + ${sci(c1)}·x + ${sci(c2)}·x²` }; } },
    exp: { label: "Exponential", scale: "log", bio: "constant fractional change per unit x — first-order kinetics.",
      fit(xs, ys) { const px = [], py = []; for (let i = 0; i < xs.length; i++) if (ys[i] > 0) { px.push(xs[i]); py.push(Math.log(ys[i])); }
        const [la, b] = wls(px, py), a = Math.exp(la);
        return { predict: (x) => a * Math.exp(b * x), r2: r2lin(px, py), params: `y = ${sci(a)}·e^(${b.toFixed(4)}·x)` }; } },
    log: { label: "Logarithmic", scale: "raw", bio: "fast early change that levels off — diminishing returns.",
      fit(xs, ys) { const lx = xs.map((x) => Math.log(Math.max(x, 1e-6))); const [a, b] = wls(lx, ys), p = (x) => a + b * Math.log(Math.max(x, 1e-6));
        return { predict: p, r2: r2on(xs, ys, p), params: `y = ${sci(a)} ${b >= 0 ? "+" : "−"} ${sci(Math.abs(b))}·ln(x)` }; } },
    power: { label: "Power law", scale: "log-log", bio: "scale-free (allometric) — a fixed % change in y per % change in x.",
      fit(xs, ys) { const px = [], py = []; for (let i = 0; i < xs.length; i++) if (xs[i] > 0 && ys[i] > 0) { px.push(Math.log(xs[i])); py.push(Math.log(ys[i])); }
        const [la, b] = wls(px, py), a = Math.exp(la);
        return { predict: (x) => a * Math.pow(Math.max(x, 1e-6), b), r2: r2lin(px, py), params: `y = ${sci(a)}·x^${b.toFixed(3)}` }; } },
    logistic: { label: "Logistic (sigmoid)", scale: "logit", bio: "a saturating switch — rises then plateaus.",
      fit(xs, ys) { const L = 1.02 * Math.max(...ys); const px = [], py = []; for (let i = 0; i < xs.length; i++) { const q = clmp(ys[i] / L, 0.001, 0.999); px.push(xs[i]); py.push(Math.log(q / (1 - q))); }
        const [a, b] = wls(px, py);
        return { predict: (x) => L / (1 + Math.exp(-clmp(a + b * x, -30, 30))), r2: r2lin(px, py), params: `L=${sci(L)}, k=${b.toFixed(3)}, x₀=${b ? (-a / b).toFixed(1) : "–"}` }; } },
    poisson: { label: "Poisson (GLM)", scale: "deviance", bio: "the canonical model for count data — a log-linear rate with variance equal to the mean.",
      fit(xs, ys) { const [b0, b1] = irls(xs, ys, POIS, logInit(xs, ys)); const mu = xs.map((x) => Math.exp(b0 + b1 * x));
        return { predict: (x) => Math.exp(b0 + b1 * x), r2: devR2(xs, ys, mu), params: `log μ = ${b0.toFixed(2)} ${b1 >= 0 ? "+" : "−"} ${Math.abs(b1).toFixed(4)}·x` }; } },
    negbin: { label: "Negative binomial (GLM)", scale: "deviance", bio: "count model for OVER-dispersed data (variance ≫ mean) — the standard for transcript counts.",
      fit(xs, ys) { const [p0, p1] = irls(xs, ys, POIS, logInit(xs, ys)); const m0 = xs.map((x) => Math.exp(p0 + p1 * x));
        let nu = 0, de = 0; for (let i = 0; i < ys.length; i++) { nu += m0[i] * m0[i]; de += Math.max((ys[i] - m0[i]) ** 2 - m0[i], 0); } const th = de > 0 ? nu / de : 1e6;
        const [b0, b1] = irls(xs, ys, NBfam(th), [p0, p1]); const mu = xs.map((x) => Math.exp(b0 + b1 * x));
        return { predict: (x) => Math.exp(b0 + b1 * x), r2: devR2(xs, ys, mu), params: `log μ = ${b0.toFixed(2)} ${b1 >= 0 ? "+" : "−"} ${Math.abs(b1).toFixed(4)}·x · θ=${sci(th)}` }; } },
    binomial: { label: "Binomial (GLM, logit)", scale: "deviance", bio: "models y as a fraction of a ceiling via a logit link — a saturating S-curve.",
      fit(xs, ys) { const N = Math.round(1.02 * Math.max(...ys)); const init = wls(xs, ys.map((y) => Math.log((y + 0.5) / (N - y + 0.5))));
        const [b0, b1] = irls(xs, ys, BINfam(N), init); const pred = (x) => N / (1 + Math.exp(-clmp(b0 + b1 * x, -30, 30)));
        return { predict: pred, r2: binDevR2(ys, xs.map(pred), N), params: `N=${sci(N)}, logit p = ${b0.toFixed(2)} ${b1 >= 0 ? "+" : "−"} ${Math.abs(b1).toFixed(4)}·x` }; } },
    loess: { label: "LOESS (local smoother)", scale: "raw", bio: "non-parametric — lets the data show its own trend with no assumed functional form.",
      fit(xs, ys) { const p = loessPredictor(xs, ys, 0.6); return { predict: p, r2: r2on(xs, ys, p), params: "local linear · span 0.6" }; } },
  };
  const SCALE_LABEL = { raw: "R²", log: "R² (log)", "log-log": "R² (log-log)", logit: "R² (logit)", deviance: "pseudo-R²" };
  function fitModel(type, xs, ys) {
    const m = MODELS[type] || MODELS.linear;
    let res; try { res = m.fit(xs.slice(), ys.slice()); } catch (_) { res = MODELS.linear.fit(xs.slice(), ys.slice()); }
    return { ...res, type, label: m.label, bio: m.bio, scale: m.scale, n: xs.length };
  }
  // shared "model · equation · R² · Pearson r · p" line; opts.xName / opts.yName name the variables in the p tooltip
  function statsHtml(xs, ys, fit, opts) {
    opts = opts || {};
    const xName = opts.xName || "x", yName = opts.yName || "y";
    const n = xs.length, r = pearsonR(xs, ys);
    const pA = pearsonP(r, n), pP = permP(xs, ys, r, 2000);
    const eq = fit.params ? ` · <span class="pn-params">${fit.params}</span>` : "";
    const rTxt = isFinite(r) ? (r >= 0 ? "+" : "") + r.toFixed(3) : "–";
    const tip = `p tests H0: ${xName} & ${yName} are uncorrelated (two-sided). ` +
      `Pearson t-test p = ${fmtP(pA)}, df ${Math.max(n - 2, 0)}. Label-permutation null p ≈ ${fmtP(pP)} (2000 shuffles). ` +
      `R² is the ${fit.label} fit on its ${SCALE_LABEL[fit.scale]} scale.`;
    return `<b>${fit.label}</b>${eq} · ${SCALE_LABEL[fit.scale]} = <b>${fit.r2.toFixed(3)}</b>` +
      ` · Pearson r = <b>${rTxt}</b> · <span class="pn-pval" title="${tip}">p = <b>${fmtP(pA)}</b>${pStars(pA)}</span>`;
  }

  return { MODELS, SCALE_LABEL, fitModel, statsHtml, pearsonR, pearsonP, permP, fmtP, pStars, linreg };
})();
