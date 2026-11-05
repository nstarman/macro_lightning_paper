// Seeing no straight lightning: the talk's "Placing Constraints" slides as an app. Top, one exposure:
// strikes arriving at random, as many as a Poisson process with the chosen mean happens to give.
// Bottom, the chance of each count, with "none" in red. Seeing none excludes every mean whose P(0) is
// below 1 - CL; the app says whether the chosen mean is excluded.
import { loadData } from "figpress-apps/data.js";
import { createStage, cssVar } from "figpress-apps/stage.js";
import { poisson, limit, rng, arrivals } from "app/science.js";

const el_ = (tag, props = {}) => Object.assign(document.createElement(tag), props);
const RED = "crimson", BAR = "#8e5bb3", BOLT = "#f2b705", BOLT_EDGE = "#8a6500", GREEN = "#1b8a2e";
const pct = (p) => (p >= 0.995 ? ">99%" : p < 0.005 ? "<1%" : `${Math.round(100 * p)}%`);

export default {
  async render({ model, el, signal }) {
    const source = model.get("manifest") ? model : (model.get("data_source") ?? "strike-count");
    const data = await loadData(source, { base: model.get("data_base") ?? "data/" });
    const meta = data.meta, nMax = data.arrays.n.data.length - 1;

    if (model.get("lam") == null) model.set("lam", meta.lam_default);
    if (model.get("cl") == null) model.set("cl", meta.levels[0]);
    if (model.get("seed") == null) model.set("seed", 18); // a first exposure with its strikes spread out
    model.save_changes();

    // ---- layout: controls, the drawing, and a fixed-height sentence that never moves the layout
    el.classList.add("fk-app");
    const toolbar = el_("div", { className: "fk-toolbar" });
    const lam = el_("input", { type: "range", min: -1, max: 1, step: 0.01, value: Math.log10(model.get("lam")) });
    const lamOut = el_("span", { className: "fk-readout" });
    lamOut.style.cssText = "display:inline-block; min-width:3.5em; font-variant-numeric:tabular-nums";
    const cl = el_("select");
    for (const c of meta.levels) cl.append(el_("option", { value: c, textContent: `${Math.round(100 * c)}%`, selected: c === model.get("cl") }));
    const again = el_("button", { textContent: "New exposure" });
    toolbar.append("Expected strikes", lam, lamOut, " confidence", cl, again);
    const panel = el_("div", { className: "fk-panel" });
    const plot = el_("div", { className: "fk-plot" });
    plot.style.cssText = "aspect-ratio: 16 / 9; width: 100%;";
    const say = el_("p", { className: "fk-readout" });
    say.style.cssText = "height:2.8em; line-height:1.4em; margin:6px 0 0; overflow:hidden; white-space:normal; font-variant-numeric:tabular-nums";
    panel.append(plot, say);
    el.replaceChildren(toolbar, panel);

    lam.oninput = () => { model.set("lam", 10 ** +lam.value); model.save_changes(); };
    cl.onchange = () => { model.set("cl", +cl.value); model.save_changes(); };
    again.onclick = () => { model.set("seed", model.get("seed") + 1); model.save_changes(); };

    // ---- drawing
    function bolt(ctx, x, y0, y1, color) { // a small zigzag bolt from y0 down to y1
      const h = y1 - y0, w = h * 0.22;
      ctx.beginPath(); ctx.moveTo(x + w, y0); ctx.lineTo(x - w * 0.6, y0 + h * 0.55); ctx.lineTo(x + w * 0.5, y0 + h * 0.5);
      ctx.lineTo(x - w, y1); ctx.lineTo(x + w * 0.8, y0 + h * 0.42); ctx.lineTo(x - w * 0.3, y0 + h * 0.47); ctx.closePath();
      ctx.fillStyle = color; ctx.fill(); ctx.strokeStyle = BOLT_EDGE; ctx.lineWidth = 1; ctx.stroke();
    }

    function draw(st) {
      const ctx = st.ctx, W = st.width, H = st.height;
      const fg = cssVar(el, "--fk-fg", "#222"), muted = cssVar(el, "--fk-muted", "#666"), font = cssVar(el, "--fk-font", "sans-serif");
      const m = model.get("lam"), p = poisson(m, nMax);
      const seen = arrivals(m, rng(model.get("seed")));
      ctx.fillStyle = cssVar(el, "--fk-bg", "#fff"); ctx.fillRect(0, 0, W, H);

      // top: one exposure, strikes where they happened
      const tx0 = 40, tx1 = W - 40, ty = H * 0.3;
      ctx.strokeStyle = fg; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(tx0, ty); ctx.lineTo(tx1, ty);
      for (const x of [tx0, tx1]) { ctx.moveTo(x, ty - 7); ctx.lineTo(x, ty + 7); }
      ctx.stroke();
      for (const t of seen) bolt(ctx, tx0 + t * (tx1 - tx0), ty - H * 0.17, ty - 2, BOLT);
      ctx.font = `13px ${font}`; ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillStyle = muted;
      ctx.fillText("one exposure Δt", (tx0 + tx1) / 2, ty + 10);
      ctx.textAlign = "right"; ctx.fillStyle = seen.length ? fg : RED;
      ctx.fillText(`seen: ${seen.length}`, tx1, ty + 10);

      // bottom: P(n), with n = 0 in red
      const bx0 = 40, bx1 = W - 40, by1 = H - 26, by0 = H * 0.5, n = nMax + 1, slot = (bx1 - bx0) / n;
      const top = Math.max(...p);
      for (let k = 0; k < n; k++) {
        const h = (p[k] / top) * (by1 - by0);
        ctx.fillStyle = k === 0 ? RED : BAR; ctx.fillRect(bx0 + k * slot + slot * 0.15, by1 - h, slot * 0.7, h);
      }
      ctx.strokeStyle = fg; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(bx0, by1); ctx.lineTo(bx1, by1); ctx.stroke();
      ctx.fillStyle = fg; ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.font = `12px ${font}`;
      for (let k = 0; k < n; k++) ctx.fillText(String(k), bx0 + (k + 0.5) * slot, by1 + 4);
      ctx.textAlign = "left"; ctx.textBaseline = "bottom"; ctx.fillStyle = RED; ctx.font = `13px ${font}`;
      ctx.fillText(`P(0) = e^−λ = ${pct(p[0])}`, bx0 + slot * 0.15, by0 - 4);
      ctx.textAlign = "right"; ctx.fillStyle = muted;
      ctx.fillText("chance of seeing n strikes", bx1, by0 - 4);
    }
    const stage = createStage(plot, { draw });

    function update() {
      const m = model.get("lam"), c = model.get("cl"), lim = limit(c), p0 = Math.exp(-m);
      lamOut.textContent = `λ = ${m < 1 ? m.toFixed(2) : m.toFixed(1)}`;
      const level = `${Math.round(100 * c)}%`, ruled = m > lim;
      say.innerHTML = `Seeing none has probability ${pct(p0)}. At ${level} confidence, seeing none rules out more than `
        + `${lim.toFixed(1)} expected strikes, so this rate is <b style="color:${ruled ? GREEN : RED}">${ruled ? "excluded" : "not excluded"}</b>.`;
      stage.requestDraw();
    }
    for (const key of ["lam", "cl", "seed"]) model.on(`change:${key}`, update);
    update();
    signal?.addEventListener("abort", () => {
      for (const key of ["lam", "cl", "seed"]) model.off(`change:${key}`, update);
      stage.dispose();
    });
  },
};
