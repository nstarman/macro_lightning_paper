// The paper site's tabs, in order; a tab with nothing in it is left out.
import { facts } from './facts.js';
import { config, resolve, root } from './store.js';

const infoPages = import.meta.glob('../pages/info/*.{md,mdx,astro}');

/** The webapps on the Explore tab, in syw-site.yml's order. */
export const exploreApps = () => Object.entries(facts.webapps ?? {}).filter(([, w]) => w.explore).map(([name, w]) => ({ name, ...w }));

/** What funded the paper (its `funding` on nstarkman.space), each award with its recipients' names. */
export function funding() {
  if (!root) return [];
  const paper = resolve(root);
  const name = (family) => {
    const a = paper.authors.find((x) => x.family === family);
    return a ? `${a.given[0]}. ${a.family}` : family;
  };
  return (paper.funding ?? []).map((f) => ({ ...f, who: (f.recipients ?? []).map(name).join(', ') }));
}

/** Co-authors on the collaborator map: those the snapshot places. */
export const hasMap = () => config.collaborators.people.length > 0;

/** The author's own Info pages (src/pages/info/*), by URL. */
export const infoLinks = () => Object.keys(infoPages).map((f) => f.replace('../pages/', '').replace(/\.(mdx?|astro)$/, '/')).filter((p) => p !== 'info/index/');

export function tabs() {
  return [
    { key: 'overview', label: 'Overview', href: '' },
    exploreApps().length && { key: 'explore', label: 'Explore', href: 'explore/' },
    { key: 'figures', label: 'Figures', href: 'figures/' },
    { key: 'info', label: 'Info', href: 'info/' }, // always: it has the showyourwork button
    { key: 'cite', label: 'Cite', href: 'cite/' },
  ].filter(Boolean);
}
