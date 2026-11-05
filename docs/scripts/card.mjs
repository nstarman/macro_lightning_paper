// syw-site card: what nstarkman.space knows about this paper, kept here so the site
// builds offline. Writes card.json (the snapshot: research-kit's collect()), the
// images it draws (public/card/) and the author's portrait (portrait.jpg, stamped on
// the figures), all from one commit of nstarman/nstarman.github.io.
//   node scripts/card.mjs <website_id> [ref]
import fs from 'node:fs';
import path from 'node:path';
import { collect } from '@nstarkman/research-kit/collect.js';
import { fromGitHub } from '@nstarkman/research-kit/readers.js';

const REPO = 'nstarman/nstarman.github.io';
const [id, ref = 'main'] = process.argv.slice(2);
if (!id) throw new Error('usage: node scripts/card.mjs <website_id> [ref]');

const snapshot = await collect(id, fromGitHub({ repo: REPO, ref }));
const sha = snapshot.source.sha;
const raw = async (p) => {
  const r = await fetch(`https://raw.githubusercontent.com/${REPO}/${sha}/${p}`);
  if (!r.ok) throw new Error(`${p}: ${r.status}`);
  return Buffer.from(await r.arrayBuffer());
};

fs.mkdirSync('public/card', { recursive: true });
const images = snapshot.items.flatMap((i) => [i.highlight?.image, i.photo]).filter(Boolean);
for (const image of images) fs.writeFileSync(path.join('public/card', path.basename(image)), await raw(`public/${image}`));
fs.writeFileSync('portrait.jpg', await raw('cv/assets/portrait.jpg'));
fs.writeFileSync('card.json', `${JSON.stringify(snapshot, null, 2)}\n`);
console.log(`card.json: ${id} at ${REPO}@${sha.slice(0, 7)}, ${snapshot.items.length} items, ${images.length} images`);
