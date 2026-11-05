// The paper's LaTeX as HTML: math through KaTeX (an error stops the build), siunitx's \qty, \unit,
// \num and their ranges, and the text commands papers use. Not a TeX engine: extend it when a paper
// needs more.
import katex from 'katex';
import { facts } from './facts.js';

const UNITS = {
  kpc: 'kpc', pc: 'pc', parsec: 'pc', Mpc: 'Mpc', km: 'km', meter: 'm', metre: 'm', gram: 'g',
  second: 's', Gyr: 'Gyr', Myr: 'Myr', year: 'yr', degree: '^{\\circ}', deg: '^{\\circ}',
  arcminute: '^{\\prime}', arcsecond: '^{\\prime\\prime}', Msun: 'M_{\\odot}', percent: '\\%',
};
const PREFIXES = { giga: 'G', mega: 'M', kilo: 'k', centi: 'c', milli: 'm', micro: '\\mu' };
const rm = (u) => (/^[A-Za-z]+$/.test(u) ? `\\mathrm{${u}}` : u);

// \qty{1.2e-10}{\kilo\meter\per\second\squared} → 1.2\times10^{-10}\,\mathrm{km}\,\mathrm{s}^{-2}
function units(spec) {
  const parts = [];
  let per = false, prefix = '', power = 1;
  for (const [, name] of spec.matchAll(/\\([A-Za-z]+)/g)) {
    if (name === 'unit') continue;  // a paper macro such as \kpc → \unit{\kilo\parsec}, inside \qty
    else if (name === 'per') per = true;
    else if (name === 'square' || name === 'cubic') power = name === 'square' ? 2 : 3;
    else if (name === 'squared' || name === 'cubed') parts.at(-1).exp *= name === 'squared' ? 2 : 3;
    else if (name in PREFIXES) prefix += PREFIXES[name];
    else {
      const u = UNITS[name] ?? name;
      parts.push({ body: /^[A-Za-z]+$/.test(prefix + u) ? rm(prefix + u) : rm(prefix) + rm(u), exp: per ? -power : power });
      per = false; prefix = ''; power = 1;
    }
  }
  return parts.map(({ body, exp }) => (exp === 1 ? body : `{${body}}^{${exp}}`)).join('\\,');
}
// 1.2e-10 → 1.2\times10^{-10}; e12 → 10^{12}
const num = (n) => n.trim().replace(/^(.*?)e([+-]?\d+)$/, (_, m, e) => `${m ? `${m}\\times` : ''}10^{${Number(e)}}`);
const U = '((?:[^{}]|\\{[^}]*\\})*)';  // a unit spec: one level of braces
const SI = [
  [`\\\\qtyrange\\{([^}]*)\\}\\{([^}]*)\\}\\{${U}\\}`, (_, a, b, u) => `${num(a)}\\text{--}${num(b)}\\,${units(u)}`],
  [`\\\\qty\\{([^}]*)\\}\\{${U}\\}`, (_, v, u) => `${num(v)}\\,${units(u)}`],
  [`\\\\unit\\{${U}\\}`, (_, u) => units(u)],
  ['\\\\numrange\\{([^}]*)\\}\\{([^}]*)\\}', (_, a, b) => `${num(a)}\\text{--}${num(b)}`],
  ['\\\\num\\{([^}]*)\\}', (_, v) => num(v)],
];
export const siunitx = (s) => SI.reduce((s, [re, f]) => s.replace(new RegExp(re, 'g'), f), s);

// The paper's own \newcommand / \DeclareMathOperator definitions (facts.json, from ms.tex and its
// own packages), with their siunitx; \ensuremath (in math, a no-op) and \textsl are LaTeX's, not KaTeX's.
const macros = {
  '\\ensuremath': '#1', '\\textsl': '\\textit{#1}',
  ...Object.fromEntries(Object.entries(facts.macros ?? {}).map(([k, v]) => [k, siunitx(v)])),
};
// Those without arguments, for text: \survey → \textsl{Survey}.
const TEXT_MACROS = Object.entries(facts.macros ?? {}).filter(([, v]) => !v.includes('#'));

const esc = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
export const math = (s) => katex.renderToString(siunitx(s), { throwOnError: true, output: 'html', macros });

// A textual citation from a `Surname:2023` or `Surname+:2023` (et al.) key: "Surname et al. (2023)".
// Other keys go: the site has no bibliography.
const citet = (keys) => keys.split(',').map((k) => k.trim().match(/^([A-Za-z-]+)(\+?):(\d{4})/))
  .filter(Boolean).map(([, n, plus, y]) => `${n}${plus ? ' et al.' : ''} (${y})`).join('; ');

function text(s) {
  return esc(s)
    .replace(/\\(?:citet|citetalias)\*?(?:\[[^\]]*\])*\{([^}]*)\}/g, (_, k) => citet(k))
    .replace(/[\s~]*\\(?:cite[pt]?|citealp|citeauthor|citepalias)\*?(?:\[[^\]]*\])*\{[^}]*\}/g, '')
    .replace(/\\(?:emph|textit|textsl)\{([^}]*)\}/g, '<em>$1</em>')
    .replace(/\\textbf\{([^}]*)\}/g, '<strong>$1</strong>')
    .replace(/\\url\{([^}]*)\}/g, (_, u) => { u = u.replace(/"/g, '&quot;').replace(/-/g, '&#45;'); return `<a href="${u}">${u}</a>`; })
    .replace(/\\([%&_#$])/g, '$1')
    .replace(/---/g, '—').replace(/--/g, '–').replace(/~/g, '&nbsp;')
    .replace(/``/g, '“').replace(/''/g, '”')
    .replace(/\\[A-Za-z]+\*?\{([^}]*)\}/g, '$1')
    .replace(/\\[A-Za-z]+\*?/g, '')
    .replace(/[{}]/g, '');
}

const MATH_IN_TEXT = new RegExp([...SI.map(([re]) => re), `\\\\ensuremath\\{${U}\\}`].join('|'), 'g');

// Running text: the paper's macros expanded, siunitx and \ensuremath rendered as math, the rest as text.
function prose(s) {
  for (const [name, body] of TEXT_MACROS) s = s.replace(new RegExp(`\\${name}(?![A-Za-z])(?:\\{\\})?`, 'g'), () => body);
  const rendered = [];
  s = s.replace(MATH_IN_TEXT, (m) => `\u0000${rendered.push(math(m)) - 1}\u0000`);
  return text(s).replace(/\u0000(\d+)\u0000/g, (_, i) => rendered[i]);
}

/** Text with inline $math$; with refs (a figure's label → its number), \ref{label} becomes the
 * number and \autoref/\cref "Figure N". Other references go: \ref with the ~ before it
 * ("(Appendix~\ref{app:x})" reads "(Appendix)"), \autoref/\cref to "the paper" ("(\cref{app:x})"
 * reads "(see the paper)"). */
export function tex(s, refs = {}) {
  s = s.replace(/(~?)\\(auto|c|C)?ref\{([^}]*)\}/g, (_, tie, kind, l) => (
    l in refs ? tie + (kind ? 'Figure~' : '') + refs[l] : kind ? `${tie}the paper` : ''));
  s = s.replace(/(^\s*|[.!?]\s+)the paper/g, '$1The paper');
  s = s.replace(/\(the paper\)/g, '(see the paper)');
  return s.split(/(\$\$[^$]*\$\$)/).map((chunk) => (chunk.startsWith('$$') && chunk.endsWith('$$') && chunk.length > 3
    ? katex.renderToString(siunitx(chunk.slice(2, -2)), { throwOnError: true, displayMode: true, macros })
    : chunk.split(/(?<!\\)\$/).map((part, i) => (i % 2 ? math(part) : prose(part))).join(''))).join('');
}

/** A caption's short title, for a gallery: its leading \emph{…} title if it has one, else its first
 * sentence, either cut at the first colon, comma, semicolon or parenthesis outside math and braces
 * ("Mock observation of the GD-1 stream: the stream's…" → "Mock observation of the GD-1 stream"). */
export function short(caption) {
  let s = caption.trim();
  const lead = s.match(/^\\(?:emph|textit|textbf)\{/);
  if (lead) {
    let depth = 1, i = lead[0].length;
    for (; i < s.length && depth; i++) depth += s[i] === '{' ? 1 : s[i] === '}' ? -1 : 0;
    s = s.slice(lead[0].length, i - 1);
  }
  let depth = 0, inMath = false;
  for (let i = 0; i < s.length; i++) {
    const c = s[i];
    if (c === '\\') { i++; continue; } // a command or an escaped character (\, \$ …)
    if (c === '$') inMath = !inMath;
    else if (c === '{') depth++;
    else if (c === '}') depth--;
    else if (!inMath && depth === 0 && (':;,('.includes(c) || (c === '.' && /\s/.test(s[i + 1] ?? ' ')))) return s.slice(0, i).trim();
  }
  return s;
}

/** A \variable value: math when it has TeX in it, text otherwise. */
export const value = (s) => (s.includes('$') ? tex(s) : /[\\^_]/.test(s) ? math(s) : prose(s));
