// Run after export.py: uv run figpress-apps test webapps/macro_explorer
// Checks the JS lightning recast against scripts/lightning.py (exported as check_lightning_<planet>).
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { loadData } from "figpress-apps/data.js";
import { lightningGrid, combine, colorize, contourSegments, densitySigma, insideMask } from "app/science.js";

const fileFetch = async (url) => {
  const buf = await readFile(new URL(url, `file://${process.cwd()}/`));
  return { ok: true, json: async () => JSON.parse(buf), arrayBuffer: async () => buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength) };
};
// APP_DATA: where the dataset was exported (syw-site: docs/.syw-site/webapps/macro-explorer/data/).
const data = await loadData("macro-constraints", { fetch: fileFetch, base: process.env.APP_DATA ?? "build/data/" });
const A = data.arrays, P = data.meta.lightning;

for (const planet of Object.keys(P.planets)) {
  test(`${planet} lightning grid matches Python at the default exposure`, () => {
    const f = lightningGrid(A.log10m.data, A.log10s.data, A[`vbar_${planet}`].data, P, P.planets[planet], P.years_default);
    const ref = A[`check_lightning_${planet}`].data;
    let compared = 0;
    for (let c = 0; c < f.length; c++) {
      const py = ref[c];
      if (py >= 99) { assert.ok(f[c] >= 99, `cell ${c}: js ${f[c]} should be unconstrained`); continue; } // Python caps at 99
      assert.ok(Math.abs(f[c] / py - 1) < 1e-4, `cell ${c}: js ${f[c]} py ${py}`);
      compared++;
    }
    assert.ok(compared > 100, `only ${compared} constrained cells`);
  });
}

test("lightning scales as 1/exposure", () => {
  const g = (yr) => lightningGrid(A.log10m.data, A.log10s.data, A.vbar_earth.data, P, P.planets.earth, yr);
  const f1 = g(2), f2 = g(20), c = f1.findIndex((x) => x < 1);
  assert.ok(Math.abs(f1[c] / f2[c] - 10) < 1e-3);
});

test("combine takes the minimum and records the limiting grid", () => {
  const { f, who } = combine([new Float32Array([0.5, 3, 99]), new Float32Array([0.1, 2, 99])], 3);
  assert.deepEqual(Array.from(f), [Float32Array.of(0.1)[0], 2, 99]);
  assert.deepEqual(Array.from(who), [1, -1, -1]);
});

test("colorize and insideMask leave unconstrained cells transparent and flip rows", () => {
  const lut = new Uint8Array(256 * 4).fill(7);
  const px = colorize(new Float32Array([99, 1e-3]), 1, 2, lut, [-15, 0.04]);
  assert.equal(px[3], 255); // row 0 drawn = largest sigma = input row 1 (constrained)
  assert.equal(px[7], 0);
  const m = insideMask(new Float32Array([99, 1e-3]), 1, 2);
  assert.deepEqual([m[3], m[7]], [255, 0]);
});

test("contour of a square region closes and sits halfway between cells", () => {
  // 4x4 grid, inside (f = 0.01) on the central 2x2: log f = -2 inside, +2 outside -> crossings at t = 0.5
  const n = 4, f = new Float32Array(n * n).fill(99);
  for (const i of [1, 2]) for (const j of [1, 2]) f[i * n + j] = 0.01;
  const seg = contourSegments(f, n, n);
  assert.equal(seg.length % 4, 0);
  const xs = [], ys = [];
  for (let k = 0; k < seg.length; k += 2) { xs.push(seg[k]); ys.push(seg[k + 1]); }
  assert.deepEqual([Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)], [0.5, 2.5, 0.5, 2.5]);
  // closed: every endpoint appears an even number of times
  const count = new Map();
  for (let k = 0; k < seg.length; k += 2) { const key = `${seg[k]},${seg[k + 1]}`; count.set(key, (count.get(key) ?? 0) + 1); }
  for (const c of count.values()) assert.equal(c % 2, 0);
});

test("nuclear-density line matches macrolimits' densities()", () => {
  // macrolimits: sigma = pi ((3/4pi) m / (rho/1.8e-24))^(2/3) with m in GeV
  const rho = 2.8e14, m = 1e30;
  const ref = Math.log10(Math.PI * ((3 / (4 * Math.PI)) * (m / (rho / 1.8e-24))) ** (2 / 3));
  assert.ok(Math.abs(densitySigma(30, rho, 1.8e-24) - ref) < 1e-9);
});
