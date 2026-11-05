// Macro constraint explorer, styled like the slides' constraint plane: the Picker+2026 constraints as
// an f_DM heat map with outlined, labelled regions; ice / nuclear / chosen-density lines; the
// black-hole region; and projected straight lightning on Earth and Jupiter (hatched).
// Toggle constraints, vary the lightning exposure, hover for f_DM. A region's name shows only while it
// is hovered, on the plot or in the list; where regions overlap the smallest wins, so each one that
// is not wholly covered by smaller ones can be reached.
import { loadData } from "figpress-apps/data.js";
import { createStage, cssVar } from "figpress-apps/stage.js";
import {
  lightningGrid, combine, colorize, insideMask, contourSegments, densitySigma, nearest, readout,
} from "app/science.js";

const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
const pow10 = (k) => "10" + String(k).split("").map((c) => SUP[c]).join("");
const el_ = (tag, props = {}) => Object.assign(document.createElement(tag), props);
const HOT = "crimson";
const PLANETS = [ // projected straight lightning: id, label, colour, hatch direction
  { id: "earth", label: "Earth lightning (proj.)", color: "crimson", dir: 1, side: "below-left" },
  { id: "jupiter", label: "Jupiter lightning (proj.)", color: "#0e9bb0", dir: -1, side: "right" },
];
const REFS0 = { ice: true, nuclear: true, bh: true }; // reference layers shown by default
const REF_LINES = [ // fixed density lines, as in macrolimits' densities()
  { label: "Ice", rho: 0.92, color: "darkorange", dash: [10, 6] },
  { label: "Nuclear", rho: 2.8e14, color: "firebrick", dash: [12, 4, 3, 4] },
];

export default {
  async render({ model, el, signal }) {
    const source = model.get("manifest") ? model : (model.get("data_source") ?? "macro-constraints");
    const data = await loadData(source, { base: model.get("data_base") ?? "data/" });
    const A = data.arrays, meta = data.meta, P = meta.lightning;
    const lm = A.log10m.data, ls = A.log10s.data, nm = lm.length, ns = ls.length;
    const keys = meta.constraints.map((c) => c.key);
    const labels = meta.constraints.map((c) => c.label);

    if (!model.get("enabled")) model.set("enabled", Object.fromEntries(keys.map((k) => [k, true])));
    const L0 = { earth: true, jupiter: true, years: P.years_default, eff: 1 };
    if (!model.get("lightning") || !("jupiter" in model.get("lightning"))) model.set("lightning", L0);
    if (!model.get("refs")) model.set("refs", REFS0);
    model.save_changes();
    let hot = null; // constraint key or planet id under the pointer in the list
    let under = null; // constraint key or planet id under the pointer on the plot

    // ---- layout
    el.classList.add("fk-app");
    const toolbar = el_("div", { className: "fk-toolbar" });
    const years = el_("input", { type: "range", min: -1, max: 3, step: 0.05, value: Math.log10(model.get("lightning").years) });
    const eff = el_("input", { type: "range", min: Math.log10(1 / 8), max: 0, step: 0.01, value: Math.log10(model.get("lightning").eff) });
    const yearsOut = el_("span", { className: "fk-readout" }), effOut = el_("span", { className: "fk-readout" });
    // fixed widths and digits that don't change width, so dragging a slider doesn't shift the toolbar
    for (const o of [yearsOut, effOut]) o.style.cssText = "display:inline-block; min-width:4.5em; font-variant-numeric:tabular-nums";
    const reset = el_("button", { textContent: "Reset" });
    toolbar.append("Lightning exposure", years, yearsOut, " efficiency", eff, effOut, reset);

    const split = el_("div", { className: "fk-split" });
    const panel = el_("div", { className: "fk-panel" });
    const head = el_("div", { className: "fk-panel-head" });
    const info = el_("span", { className: "fk-readout", textContent: "hover the plot" });
    // the hover readout gets a fixed two-line box and wraps inside it, so its changing length never moves the plot
    info.style.cssText = "display:block; height:2.6em; line-height:1.3em; overflow:hidden; white-space:normal; overflow-wrap:anywhere; font-variant-numeric:tabular-nums";
    head.append(el_("span", { textContent: meta.title }), info);
    head.style.display = "block"; // title on one row, the readout always on its own row below
    const plot = el_("div", { className: "fk-plot" });
    plot.style.cssText = "aspect-ratio: 4 / 3; width: 100%;";
    panel.append(head, plot);

    const list = el_("div", { className: "fk-panel" });
    // legend swatches, drawn as on the plot: a dashed line, a filled box or a hatched box
    const swatch = (kind, color, style) => {
      const c = el_("canvas", { width: 56, height: 24 }), x = c.getContext("2d");
      c.style.cssText = "width:28px;height:12px;vertical-align:middle;margin:0 4px 0 6px";
      x.scale(2, 2);
      if (kind === "line") { x.strokeStyle = color; x.lineWidth = 2; x.setLineDash(style); x.beginPath(); x.moveTo(1, 6); x.lineTo(27, 6); x.stroke(); }
      if (kind === "box") { x.fillStyle = color; x.fillRect(1, 1, 26, 10); x.fillStyle = "#777"; for (let i = 3; i < 27; i += 5) for (let j = 3; j < 11; j += 4) x.fillRect(i, j, 1, 1); }
      if (kind === "hatch") { x.fillStyle = hatchPattern(x, color, style); x.fillRect(1, 1, 26, 10); x.strokeStyle = color; x.lineWidth = 1.2; x.strokeRect(1, 1, 26, 10); }
      return c;
    };
    const row = (input, text, key, sw) => {
      const r = el_("label", { style: "display:block; cursor:pointer; break-inside:avoid; white-space:nowrap" });
      if (input) r.append(input); else r.style.paddingLeft = "20px";
      if (sw) r.append(sw);
      r.append(" " + text);
      r.onpointerenter = () => { hot = key; stage.requestDraw(); };
      r.onpointerleave = () => { hot = null; stage.requestDraw(); };
      return r;
    };
    // The constraints, the long group, flow into as many columns as fit; Projected and Reference stack
    // in a column beside them when there is room, and drop below when there is not.
    list.style.cssText = "display:flex; flex-wrap:wrap; gap:10px 24px; align-content:start";
    const side = el_("div", { style: "flex:1 1 12em; min-width:0; display:flex; flex-direction:column; gap:10px" });
    const group = (title, parent) => {
      const g = el_("div", { style: parent === list ? "flex:3 1 20em; min-width:0" : "" });
      g.append(el_("div", { className: "fk-panel-head", textContent: title }));
      parent.append(g);
      return g;
    };
    const boxes = {};
    const gCons = group("Constraints (Picker+2026)", list);
    list.append(side);
    const cols = el_("div", { style: "columns:9em; column-gap:20px" });
    for (const c of meta.constraints) {
      const cb = el_("input", { type: "checkbox", checked: model.get("enabled")[c.key] });
      cb.onchange = () => { model.set("enabled", { ...model.get("enabled"), [c.key]: cb.checked }); model.save_changes(); };
      boxes[c.key] = cb;
      cols.append(row(cb, c.label, c.key));
    }
    const allOn = el_("button", { textContent: "All" }), allOff = el_("button", { textContent: "None" });
    allOn.onclick = () => setAll(true); allOff.onclick = () => setAll(false);
    const buttons = el_("div", { style: "margin-top:6px" });
    buttons.append(allOn, " ", allOff);
    gCons.append(cols, buttons);
    const gProj = group("Projected (Starkman+2021)", side);
    const pboxes = {};
    for (const p of PLANETS) {
      const cb = el_("input", { type: "checkbox", checked: model.get("lightning")[p.id] });
      cb.onchange = () => { setLightning({ [p.id]: cb.checked }); model.save_changes(); };
      pboxes[p.id] = cb;
      gProj.append(row(cb, p.label, p.id, swatch("hatch", p.color, p.dir)));
    }
    const gRefs = group("Reference", side);
    const rboxes = {};
    const refRow = (id, text, sw) => {
      const cb = el_("input", { type: "checkbox", checked: model.get("refs")[id] });
      cb.onchange = () => { model.set("refs", { ...model.get("refs"), [id]: cb.checked }); model.save_changes(); };
      rboxes[id] = cb;
      return row(cb, text, null, sw);
    };
    for (const r of REF_LINES) gRefs.append(refRow(r.label.toLowerCase(), `${r.label} density`, swatch("line", r.color, r.dash)));
    gRefs.append(refRow("bh", "Black holes", swatch("box", "#111")));
    split.append(panel, list);
    el.replaceChildren(toolbar, split);

    function setAll(v) {
      for (const k of keys) boxes[k].checked = v;
      model.set("enabled", Object.fromEntries(keys.map((k) => [k, v])));
      model.save_changes();
    }
    const setLightning = (patch) => { model.set("lightning", { ...model.get("lightning"), ...patch }); };
    years.oninput = () => setLightning({ years: 10 ** +years.value });
    eff.oninput = () => setLightning({ eff: 10 ** +eff.value });
    for (const s of [years, eff]) s.onchange = () => model.save_changes();
    reset.onclick = () => {
      years.value = Math.log10(P.years_default); eff.value = 0;
      for (const p of PLANETS) pboxes[p.id].checked = true;
      model.set("lightning", L0);
      for (const id in rboxes) rboxes[id].checked = true;
      model.set("refs", REFS0);
      setAll(true);
    };

    // ---- science: static outlines once; heat map and projections when their inputs change
    const outline = Object.fromEntries(keys.map((k) => [k, contourSegments(A[`f_${k}`].data, nm, ns)]));
    const areaOf = (f) => { let n = 0; for (let c = 0; c < f.length; c++) if (f[c] <= 1) n++; return n; };
    const area = Object.fromEntries(keys.map((k) => [k, areaOf(A[`f_${k}`].data)])); // planets: in recompute
    const canvasOf = () => { const c = el_("canvas", { width: nm, height: ns }); return [c, c.getContext("2d")]; };
    const [image, ictx] = canvasOf();
    const proj = Object.fromEntries(PLANETS.map((p) => {
      const [mask, mctx] = canvasOf();
      return [p.id, { f: null, key: "", mask, mctx, seg: null, at: null }];
    }));
    let combined = null;
    function recompute() {
      const L = model.get("lightning");
      for (const p of PLANETS) {
        const pr = proj[p.id], k = `${L.years}|${L.eff}`;
        if (!L[p.id] || pr.key === k) continue;
        pr.f = lightningGrid(lm, ls, A[`vbar_${p.id}`].data, P, P.planets[p.id], L.years, L.eff);
        pr.key = k;
        pr.seg = contourSegments(pr.f, nm, ns);
        pr.mctx.putImageData(new ImageData(insideMask(pr.f, nm, ns), nm, ns), 0, 0);
        let j0 = Infinity, j1 = -1, i0 = Infinity, i1 = -1; // bounding box of the constrained cells
        for (let i = 0; i < ns; i++) for (let j = 0; j < nm; j++) if (pr.f[i * nm + j] <= 1) {
          j0 = Math.min(j0, j); j1 = Math.max(j1, j); i0 = Math.min(i0, i); i1 = Math.max(i1, i);
        }
        pr.at = j1 >= 0 ? { j0, j1, i0, i1 } : null;
        area[p.id] = areaOf(pr.f);
      }
      const en = model.get("enabled");
      const idx = keys.map((k, i) => (en[k] ? i : -1)).filter((i) => i >= 0);
      const c = combine(idx.map((i) => A[`f_${keys[i]}`].data), nm * ns);
      c.who = c.who.map((w) => (w < 0 ? -1 : idx[w]));
      combined = c;
      ictx.putImageData(new ImageData(colorize(c.f, nm, ns, A.viridis.data, meta.f_range), nm, ns), 0, 0);
      yearsOut.textContent = `${L.years.toPrecision(2)} yr`;
      effOut.textContent = `${L.eff.toFixed(2)}`;
      stage.requestDraw();
    }

    // ---- drawing
    const pad = { l: 64, r: 88, t: 34, b: 44 };
    const bar = el_("canvas", { width: 1, height: 256 }); // colour bar, top = f_DM 1
    const lut = A.viridis.data, bpx = new Uint8ClampedArray(4 * 256);
    for (let k = 0; k < 256; k++) bpx.set(lut.subarray(4 * (255 - k), 4 * (255 - k) + 4), 4 * k);
    bar.getContext("2d").putImageData(new ImageData(bpx, 1, 256), 0, 0);
    const frame = (st) => ({ x0: pad.l, y0: pad.t, w: st.width - pad.l - pad.r, h: st.height - pad.t - pad.b });
    const dm = (lm[nm - 1] - lm[0]) / (nm - 1), ds = (ls[ns - 1] - ls[0]) / (ns - 1);
    const toX = (f, v) => f.x0 + ((v - lm[0]) / (lm[nm - 1] - lm[0])) * f.w;
    const toY = (f, v) => f.y0 + f.h - ((v - ls[0]) / (ls[ns - 1] - ls[0])) * f.h;
    const gx = (f, j) => toX(f, lm[0] + j * dm), gy = (f, i) => toY(f, ls[0] + i * ds); // grid index → px

    function strokeSegments(ctx, f, seg, color, width, alpha = 1) {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.globalAlpha = alpha; ctx.beginPath();
      for (let k = 0; k < seg.length; k += 4) {
        ctx.moveTo(gx(f, seg[k]), gy(f, seg[k + 1])); ctx.lineTo(gx(f, seg[k + 2]), gy(f, seg[k + 3]));
      }
      ctx.stroke(); ctx.globalAlpha = 1;
    }
    function hatchPattern(ctx, color, dir) {
      const t = el_("canvas", { width: 10, height: 10 }), c = t.getContext("2d");
      c.strokeStyle = color; c.lineWidth = 1.5; c.beginPath();
      for (const o of [-10, 0, 10]) {
        if (dir > 0) { c.moveTo(o, 10); c.lineTo(o + 10, 0); } else { c.moveTo(o, 0); c.lineTo(o + 10, 10); }
      }
      c.stroke();
      return ctx.createPattern(t, "repeat");
    }
    const off = el_("canvas"); // scratch for masked hatching, in CSS px
    function drawHatched(ctx, f, mask, color, dir) {
      off.width = Math.ceil(f.w); off.height = Math.ceil(f.h);
      const o = off.getContext("2d");
      o.imageSmoothingEnabled = false;
      o.drawImage(mask, 0, 0, f.w, f.h);
      o.globalCompositeOperation = "source-in";
      o.fillStyle = hatchPattern(o, color, dir); o.fillRect(0, 0, f.w, f.h);
      o.globalCompositeOperation = "source-over";
      ctx.drawImage(off, f.x0, f.y0, f.w, f.h);
    }
    // Labels stay inside the plot frame `f`: unrotated ones are shifted back in, all are clipped to it.
    function label(ctx, text, x, y, rot, color, bold, f) {
      ctx.save();
      ctx.beginPath(); ctx.rect(f.x0, f.y0, f.w, f.h); ctx.clip();
      ctx.font = `${bold ? "600 " : ""}13px ${cssVar(el, "--fk-font", "sans-serif")}`;
      if (!rot) {
        const w = ctx.measureText(text).width;
        x = Math.min(Math.max(x, f.x0 + 3), f.x0 + f.w - w - 3);
        y = Math.min(Math.max(y, f.y0 + 16), f.y0 + f.h - 3);
      }
      ctx.translate(x, y); ctx.rotate(-rot * Math.PI / 180);
      ctx.textAlign = "left"; ctx.textBaseline = "bottom";
      ctx.lineWidth = 3; ctx.strokeStyle = cssVar(el, "--fk-bg", "#fff"); ctx.strokeText(text, 0, 0);
      ctx.fillStyle = color; ctx.fillText(text, 0, 0); ctx.restore();
    }
    // Math set on a canvas, centred on the origin: [text, "sub" | "it"] runs, subscripts smaller and lowered.
    function mathText(ctx, runs, size, family) {
      const fontOf = (k) => (k === "sub" ? `${Math.round(size * 0.7)}px ${family}` : `${k === "it" ? "italic " : ""}${size}px ${family}`);
      const w = runs.map(([t, k]) => { ctx.font = fontOf(k); return ctx.measureText(t).width; });
      let x = -w.reduce((a, b) => a + b, 0) / 2;
      ctx.textAlign = "left"; ctx.textBaseline = "alphabetic";
      runs.forEach(([t, k], i) => { ctx.font = fontOf(k); ctx.fillText(t, x, k === "sub" ? size * 0.25 + 4 : 4); x += w[i]; });
    }
    function densityLine(ctx, f, rhoV, color, dash, width = 2) {
      ctx.strokeStyle = color; ctx.setLineDash(dash); ctx.lineWidth = width; ctx.beginPath();
      ctx.moveTo(toX(f, lm[0]), toY(f, densitySigma(lm[0], rhoV, P.gev2g)));
      ctx.lineTo(toX(f, lm[nm - 1]), toY(f, densitySigma(lm[nm - 1], rhoV, P.gev2g)));
      ctx.stroke(); ctx.setLineDash([]);
    }

    function draw(st) {
      const ctx = st.ctx, f = frame(st), fg = cssVar(el, "--fk-fg", "#222"), muted = cssVar(el, "--fk-muted", "#666");
      const font = cssVar(el, "--fk-font", "sans-serif"), en = model.get("enabled"), L = model.get("lightning");
      ctx.fillStyle = cssVar(el, "--fk-bg", "#fff");
      ctx.fillRect(0, 0, st.width, st.height);
      ctx.save(); ctx.beginPath(); ctx.rect(f.x0, f.y0, f.w, f.h); ctx.clip();
      ctx.imageSmoothingEnabled = false;
      ctx.drawImage(image, f.x0, f.y0, f.w, f.h);
      const lit = hot ?? under; // the region named and outlined in red
      for (const k of keys) if (en[k] && k !== lit) strokeSegments(ctx, f, outline[k], fg, 1, 0.45);
      // black holes: sigma below the Schwarzschild area (as in macrolimits), black with grey dots
      const refs = model.get("refs");
      if (refs.bh) {
      const bh = (m) => Math.log10(Math.PI * (100 * 2 * 6.67e-11 * 10 ** m / 5.62e23 * 1e-3 / 9e16) ** 2);
      ctx.beginPath(); ctx.moveTo(toX(f, lm[nm - 1]), toY(f, ls[0]));
      for (let j = 0; j < nm; j += 2) ctx.lineTo(toX(f, lm[j]), toY(f, Math.max(ls[0] - 1, Math.min(ls[ns - 1] + 1, bh(lm[j])))));
      ctx.lineTo(toX(f, lm[nm - 1]), toY(f, Math.min(ls[ns - 1] + 1, bh(lm[nm - 1]))));
      ctx.closePath(); ctx.fillStyle = "#111"; ctx.fill();
      const dots = el_("canvas", { width: 8, height: 8 }), dc = dots.getContext("2d");
      dc.fillStyle = "#777"; dc.fillRect(3, 3, 1.5, 1.5);
      ctx.fillStyle = ctx.createPattern(dots, "repeat"); ctx.fill();
      }
      // density lines
      for (const r of REF_LINES) if (refs[r.label.toLowerCase()]) densityLine(ctx, f, r.rho, r.color, r.dash);
      // projected straight lightning: hatched, outlined
      for (const p of PLANETS) {
        const pr = proj[p.id];
        if (!L[p.id] || !pr.f) continue;
        drawHatched(ctx, f, pr.mask, p.color, p.dir);
        strokeSegments(ctx, f, pr.seg, p.color, lit === p.id ? 3.5 : 2.2);
      }
      if (lit && outline[lit] && en[lit]) strokeSegments(ctx, f, outline[lit], HOT, 3);
      ctx.restore();
      // the lit region's name only (constraints at the slide positions; projections beside their region)
      for (const c of meta.constraints) {
        if (c.key !== lit || !en[c.key]) continue;
        const [am, as, rot] = c.at;
        label(ctx, c.label, toX(f, am), toY(f, as), rot, HOT, true, f);
      }
      for (const p of PLANETS) {
        const pr = proj[p.id];
        if (p.id !== lit || !L[p.id] || !pr.at) continue;
        const b = pr.at, text = p.label.replace(" (proj.)", "");
        if (p.side === "right") label(ctx, text, gx(f, b.j1) + 6, gy(f, (b.i0 + b.i1) / 2) + 6, 0, p.color, true, f);
        else label(ctx, text, gx(f, b.j0) - 20, gy(f, b.i0) + 20, 0, p.color, true, f);
      }
      // colour bar: log10 f_DM from f_range[1] (top) to f_range[0] (bottom), label rotated beside it
      const [vmin, vmax] = meta.f_range, bx = f.x0 + f.w + 14;
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(bar, bx, f.y0, 12, f.h);
      ctx.fillStyle = fg; ctx.textAlign = "left"; ctx.textBaseline = "middle"; ctx.font = `11px ${font}`;
      for (let k = 0; k >= vmin; k -= 3) ctx.fillText(String(k).replace("-", "−"), bx + 16, f.y0 + ((vmax - k) / (vmax - vmin)) * f.h);
      ctx.save(); ctx.translate(bx + 48, f.y0 + f.h / 2); ctx.rotate(-Math.PI / 2);
      mathText(ctx, [["log"], ["10", "sub"], ["\u2009"], ["f", "it"], ["DM", "sub"]], 13, font); ctx.restore();
      // axes: GeV below, grams above
      ctx.strokeStyle = fg; ctx.lineWidth = 1; ctx.strokeRect(f.x0, f.y0, f.w, f.h);
      ctx.font = `12px ${font}`; ctx.textAlign = "center"; ctx.textBaseline = "top";
      for (let k = 0; k <= 50; k += 10) ctx.fillText(pow10(k), toX(f, k), f.y0 + f.h + 4);
      ctx.fillText("dark matter mass [GeV]", f.x0 + f.w / 2, f.y0 + f.h + 24);
      ctx.textBaseline = "bottom"; ctx.fillStyle = muted; ctx.strokeStyle = muted;
      const g0 = Math.log10(P.gev2g);
      for (let k = -20; k <= 25; k += 5) {
        const x = toX(f, k - g0);
        if (x < f.x0 || x > f.x0 + f.w) continue;
        ctx.beginPath(); ctx.moveTo(x, f.y0); ctx.lineTo(x, f.y0 - 4); ctx.stroke();
        ctx.fillText(pow10(k), x, f.y0 - 5);
      }
      ctx.textAlign = "left"; ctx.fillText("g", f.x0 + f.w + 4, f.y0 - 5);
      ctx.fillStyle = fg; ctx.textAlign = "right"; ctx.textBaseline = "middle";
      for (let k = -40; k <= 20; k += 10) ctx.fillText(pow10(k), f.x0 - 4, toY(f, k));
      ctx.save(); ctx.translate(14, f.y0 + f.h / 2); ctx.rotate(-Math.PI / 2);
      ctx.textAlign = "center"; ctx.fillStyle = muted; ctx.fillText("geometric cross section [cm²]", 0, 0); ctx.restore();
    }
    const stage = createStage(plot, { draw });

    // The region under cell c: of those enabled and constraining there (f ≤ 1), the smallest.
    function regionAt(c) {
      const en = model.get("enabled"), L = model.get("lightning");
      let best = null;
      for (const k of keys) if (en[k] && A[`f_${k}`].data[c] <= 1 && (best === null || area[k] < area[best])) best = k;
      for (const p of PLANETS) {
        const pf = proj[p.id].f;
        if (L[p.id] && pf && pf[c] <= 1 && (best === null || area[p.id] < area[best])) best = p.id;
      }
      return best;
    }
    const setUnder = (k) => { if (k !== under) { under = k; stage.requestDraw(); } };
    const hover = (e) => {
      const b = plot.getBoundingClientRect(), f = frame(stage);
      const x = e.clientX - b.left, y = e.clientY - b.top;
      if (!combined || x < f.x0 || x > f.x0 + f.w || y < f.y0 || y > f.y0 + f.h) { setUnder(null); return; }
      const vm = lm[0] + ((x - f.x0) / f.w) * (lm[nm - 1] - lm[0]);
      const vs = ls[0] + ((f.y0 + f.h - y) / f.h) * (ls[ns - 1] - ls[0]);
      const c = nearest(ls, vs) * nm + nearest(lm, vm), L = model.get("lightning");
      const projected = PLANETS.filter((p) => L[p.id] && proj[p.id].f).map((p) => [p.label, proj[p.id].f[c]]);
      info.textContent = readout(vm, vs, combined.f[c], combined.who[c], labels, P.gev2g, projected);
      setUnder(regionAt(c));
    };
    const leave = () => setUnder(null);
    plot.addEventListener("pointermove", hover);
    plot.addEventListener("pointerleave", leave);

    for (const key of ["enabled", "lightning"]) model.on(`change:${key}`, recompute);
    const redraw = () => stage.requestDraw();
    model.on("change:refs", redraw);
    recompute();
    signal?.addEventListener("abort", () => {
      for (const key of ["enabled", "lightning"]) model.off(`change:${key}`, recompute);
      model.off("change:refs", redraw);
      plot.removeEventListener("pointermove", hover);
      plot.removeEventListener("pointerleave", leave);
      stage.dispose();
    });
  },
};
