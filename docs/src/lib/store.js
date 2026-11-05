// The website's data for this paper: card.json, the snapshot `syw-site card` keeps of
// what nstarkman.space knows about it (its record, what its refs reach, the people and
// places on its map). @nstarkman/research-kit's components read it as the
// website reads its own data, so the header is the website's entry for this paper.
import fs from 'node:fs';
import path from 'node:path';
import { createStore } from '@nstarkman/research-kit/store.js';

const docs = process.cwd();
// facts.json carries the settings' website_id (src/lib/facts.js says how to make it).
const factsFile = path.join(docs, 'public', 'generated', 'facts.json');
const settings = fs.existsSync(factsFile) ? JSON.parse(fs.readFileSync(factsFile, 'utf8')) : {};
const file = path.join(docs, 'card.json');
const EMPTY = { items: [], config: { collaborators: { people: [] }, acknowledgements: { papers: [] }, journals: {}, stars: { stars: {} }, places: { places: {} } } };

let snapshot = EMPTY;
if (settings.website_id) {
  if (!fs.existsSync(file)) throw new Error('syw-site: no card.json — run `syw-site card` in the repository root');
  snapshot = JSON.parse(fs.readFileSync(file, 'utf8'));
}

// Images the card draws at the site's root (/highlights/…) are served from this
// site's own copies (public/card/), so they show locally and published alike.
const base = import.meta.env.BASE_URL.replace(/^\/|\/$/g, '');
const own = (p) => (p ? `${base}/card/${path.basename(p)}` : p);
// A link to this site itself (the paper's `homepage`, nstarkman.space/<repo>/) is dropped: here
// it would only reload the page.
const self = (l) => l.rel === 'homepage' && new URL(l.url).pathname.replace(/\/$/, '') === `/${base}`;
for (const item of snapshot.items) {
  if (item.highlight?.image) item.highlight = { ...item.highlight, image: own(item.highlight.image) };
  if (item.photo) item.photo = own(item.photo);
  if (item.links?.some(self)) item.links = item.links.filter((l) => !self(l));
}

/** The paper's id on nstarkman.space, or null before it has a record there. */
export const root = settings.website_id ? snapshot.root : null;
/** The card's preset (cards/<id>.json on the website, else paperSite). */
export const card = snapshot.card ?? { preset: 'paperSite' };

export const store = createStore(snapshot);
export const {
  items, lists, config, byType, featured, resolve, titleOf, splitTitle, HIGHLIGHT_TOPICS,
  SOFTWARE_GROUPS, displayName, authorPosition, venueUrl, softwarePapers, softwareUsed,
  affiliationAt, authors, positionStops, money, venueLine, year,
  dateLabel, relKey, REL_ICON, primaryLink, links, assists,
  assistItem, groupedByDefault, institutionGroups, groupTitle, groupDate,
} = store;
