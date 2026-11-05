// The counting behind the paper's limit: macro-induced strikes arrive at random, so the number seen in
// an exposure is Poisson with mean lam = n_ml * Delta t. Seeing none has probability exp(-lam), so
// seeing none rules out every lam with exp(-lam) < 1 - CL, that is lam > -ln(1 - CL).

/** P(n) for a Poisson mean `lam`, for n = 0 .. nMax. */
export function poisson(lam, nMax) {
  const p = [Math.exp(-lam)];
  for (let n = 1; n <= nMax; n++) p.push((p[n - 1] * lam) / n);
  return p;
}

/** The mean above which seeing no strikes is excluded at confidence `cl` (e.g. 0.95 -> 3.0). */
export const limit = (cl) => -Math.log(1 - cl);

/** A small seeded generator (mulberry32), so a "new exposure" is reproducible in tests. */
export function rng(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Arrival times in [0, 1) of a Poisson process with `lam` expected events over the window. */
export function arrivals(lam, random) {
  const out = [];
  for (let t = -Math.log(1 - random()) / lam; t < 1; t += -Math.log(1 - random()) / lam) out.push(t);
  return out;
}
