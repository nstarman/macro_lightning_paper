import { test } from "node:test";
import assert from "node:assert/strict";
import { poisson, limit, rng, arrivals } from "app/science.js";

test("Poisson probabilities sum to one and start at exp(-lam)", () => {
  const p = poisson(2.5, 40);
  assert.ok(Math.abs(p.reduce((a, b) => a + b, 0) - 1) < 1e-12);
  assert.equal(p[0], Math.exp(-2.5));
  assert.ok(Math.abs(p[3] - (Math.exp(-2.5) * 2.5 ** 3) / 6) < 1e-15);
});

test("seeing none excludes lam > -ln(1 - CL): the paper's 3 at 95%, 2.3 at 90%", () => {
  assert.ok(Math.abs(limit(0.95) - 2.9957) < 1e-4);
  assert.ok(Math.abs(limit(0.9) - 2.3026) < 1e-4);
  assert.ok(Math.abs(poisson(limit(0.95), 0)[0] - 0.05) < 1e-12);
});

test("arrivals are in the window, sorted, and average lam", () => {
  const r = rng(1);
  let total = 0;
  for (let k = 0; k < 4000; k++) {
    const a = arrivals(3, r);
    assert.ok(a.every((t, i) => t >= 0 && t < 1 && (i === 0 || t > a[i - 1])));
    total += a.length;
  }
  assert.ok(Math.abs(total / 4000 - 3) < 0.1);
});

test("the generator is reproducible from its seed", () => {
  assert.deepEqual(arrivals(5, rng(7)), arrivals(5, rng(7)));
});
