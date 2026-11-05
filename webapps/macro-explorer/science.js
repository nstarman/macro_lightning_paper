// Pure science for the macro explorer: no DOM, so `figpress-apps test` can check it against Python.
// Grids are row-major (n_sigma, n_mass): element [i, j] is sigma index i, mass index j.

/**
 * Projected straight-lightning f_DM (Starkman+2021 recast; src/scripts/lightning.py) on one planet:
 * f = n_obs / N_exp, N_exp = eff * (rho/m) * pi R^2 f_TS * T * <v; sigma>, only for sigma/m <= max.
 * `p` holds the shared constants, `planet` its radius `r` [cm] and thunderstorm fraction `f_ts`.
 * Unconstrained cells are Infinity.
 */
export function lightningGrid(log10m, log10s, vbar, p, planet, years, eff = 1) {
  const nm = log10m.length, ns = log10s.length, out = new Float32Array(nm * ns);
  const exposure = eff * p.rho_g * Math.PI * planet.r ** 2 * planet.f_ts * years * p.year_s; // g s / cm
  for (let i = 0; i < ns; i++) {
    const s = 10 ** log10s[i];
    for (let j = 0; j < nm; j++) {
      const mg = 10 ** log10m[j] * p.gev2g;
      const n = exposure * vbar[i] / mg;
      out[i * nm + j] = n > 0 && s / mg <= p.sig_over_m_max ? p.n_obs / n : Infinity;
    }
  }
  return out;
}

/** Cell-wise minimum over the given grids → {f, who}: who = index into `grids` (-1 if none ≤ 1). */
export function combine(grids, n) {
  const f = new Float32Array(n).fill(Infinity), who = new Int16Array(n).fill(-1);
  grids.forEach((g, k) => {
    for (let c = 0; c < n; c++) if (g[c] < f[c]) { f[c] = g[c]; who[c] = g[c] <= 1 ? k : who[c]; }
  });
  for (let c = 0; c < n; c++) if (!(f[c] <= 1)) who[c] = -1;
  return { f, who };
}

/** RGBA pixels (row 0 = largest sigma, for drawing) from f; transparent above 10^vmax. */
export function colorize(f, nm, ns, lut, [vmin, vmax]) {
  const px = new Uint8ClampedArray(4 * nm * ns), nl = lut.length / 4;
  for (let i = 0; i < ns; i++) {
    for (let j = 0; j < nm; j++) {
      const v = Math.log10(f[i * nm + j]);
      if (!(v <= vmax)) continue;
      const t = Math.min(1, Math.max(0, (v - vmin) / (vmax - vmin)));
      const k = 4 * Math.round(t * (nl - 1)), o = 4 * ((ns - 1 - i) * nm + j);
      px[o] = lut[k]; px[o + 1] = lut[k + 1]; px[o + 2] = lut[k + 2]; px[o + 3] = 255;
    }
  }
  return px;
}

/** Alpha mask (row 0 = largest sigma) of cells with f ≤ 1: 255 inside, 0 outside. */
export function insideMask(f, nm, ns) {
  const px = new Uint8ClampedArray(4 * nm * ns);
  for (let i = 0; i < ns; i++) for (let j = 0; j < nm; j++) {
    if (f[i * nm + j] <= 1) px[4 * ((ns - 1 - i) * nm + j) + 3] = 255;
  }
  return px;
}

// Marching-squares edge pairs per case (bit 1 = corner (i,j), 2 = (i,j+1), 4 = (i+1,j+1), 8 = (i+1,j)).
// Edges: 0 bottom (i), 1 right (j+1), 2 top (i+1), 3 left (j). Saddles 5 and 10 give two segments.
const CASES = [[], [3, 0], [0, 1], [3, 1], [1, 2], [3, 0, 1, 2], [0, 2], [3, 2],
  [2, 3], [0, 2], [0, 1, 2, 3], [1, 2], [1, 3], [0, 1], [3, 0], []];

/**
 * The f = level contour as line segments, flat [j0, i0, j1, i1, ...] in fractional grid indices.
 * Works in log10 f; unconstrained cells (≥ 99 or Infinity) count as log10 f = 2.
 */
export function contourSegments(f, nm, ns, level = 1) {
  const L = Math.log10(level), out = [];
  const v = (i, j) => { const x = f[i * nm + j]; return x > 0 && x < 99 ? Math.log10(x) : 2; };
  const t = (a, b) => (a === b ? 0.5 : (L - a) / (b - a));
  for (let i = 0; i < ns - 1; i++) {
    for (let j = 0; j < nm - 1; j++) {
      const a = v(i, j), b = v(i, j + 1), c = v(i + 1, j + 1), d = v(i + 1, j);
      const k = (a <= L) | ((b <= L) << 1) | ((c <= L) << 2) | ((d <= L) << 3);
      const edges = CASES[k];
      if (!edges.length) continue;
      const pt = [
        [j + t(a, b), i], [j + 1, i + t(b, c)], [j + t(d, c), i + 1], [j, i + t(a, d)],
      ];
      for (let e = 0; e < edges.length; e += 2) out.push(...pt[edges[e]], ...pt[edges[e + 1]]);
    }
  }
  return Float32Array.from(out);
}

/** log10 geometric cross-section [cm^2] of a constant-density sphere: sigma = pi (3m / 4 pi rho)^(2/3). */
export function densitySigma(log10mGeV, rho, gev2g) {
  return Math.log10(Math.PI) + (2 / 3) * Math.log10((3 * 10 ** log10mGeV * gev2g) / (4 * Math.PI * rho));
}

/** Nearest grid index of x in an evenly spaced increasing array. */
export function nearest(arr, x) {
  const k = Math.round(((x - arr[0]) / (arr[arr.length - 1] - arr[0])) * (arr.length - 1));
  return Math.min(arr.length - 1, Math.max(0, k));
}

const sci = (x) => x.toExponential(1).replace("e+", "e");

/**
 * Hover text for (log10 m, log10 sigma): the limiting existing constraint, then any projection
 * (`projected` = [[label, f], ...]) that reaches f ≤ 1 there.
 */
export function readout(lm, ls, f, who, labels, gev2g, projected = []) {
  const where = `m = ${sci(10 ** lm)} GeV (${sci(10 ** lm * gev2g)} g), σ = ${sci(10 ** ls)} cm²`;
  const parts = [who < 0 ? "unconstrained" : `f_DM ≤ ${sci(f)} (${labels[who]})`];
  for (const [lab, fp] of projected) if (fp <= 1) parts.push(`${lab}: ${sci(fp)}`);
  return `${where}: ${parts.join(" · ")}`;
}
