// The paper's facts, written by `syw-site build` (Snakemake rule syw_site_facts). Every page reads
// them from here; a name that is not in them stops the build (see lookup below).
import fs from 'node:fs';
import path from 'node:path';

const file = path.join(process.cwd(), 'public', 'generated', 'facts.json');
if (!fs.existsSync(file)) throw new Error('syw-site: no public/generated/facts.json — run `syw-site build` in the repository root');

export const facts = JSON.parse(fs.readFileSync(file, 'utf8'));

export function variable(name) {
  if (!(name in facts.variables)) throw new Error(`syw-site: no variable "${name}" in facts.json`);
  return facts.variables[name];
}

export function figure(label) {
  const fig = facts.figures.find((f) => f.label === label || f.slug === label);
  if (!fig) throw new Error(`syw-site: no figure "${label}" in facts.json`);
  return fig;
}

export function webapp(name) {
  if (!(name in (facts.webapps ?? {}))) throw new Error(`syw-site: no webapp "${name}" in docs/syw-site.yml`);
  return name;
}

/** A URL under the site's base, for files in public/. */
export const url = (p) => `${import.meta.env.BASE_URL.replace(/\/$/, '')}/${p.replace(/^\//, '')}`;
