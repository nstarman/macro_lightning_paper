// The docs site of nstarman/macro_lightning_paper, published at nstarkman.space/macro_lightning_paper/ (made by syw-site).
import fs from 'node:fs';
import { fileURLToPath } from 'node:url';
import { defineConfig } from 'astro/config';
import mdx from '@astrojs/mdx';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';

// $math$ in pages gets the paper's own macros, as components' math does (src/lib/tex.js); an
// unknown command stops the build. facts.json is checked with a clear message in src/lib/facts.js.
const facts = 'public/generated/facts.json';
const macros = fs.existsSync(facts) ? JSON.parse(fs.readFileSync(facts, 'utf8')).macros ?? {} : {};

export default defineConfig({
  site: 'https://nstarkman.space',
  base: '/macro_lightning_paper',
  build: { format: 'directory' },
  // The paper's card, entry and collaborator map are the website's own components
  // (@nstarkman/research-kit), reading the paper's snapshot of the website's data
  // (card.json) through src/lib/store.js. Vite must load those packages itself, in dev
  // and in the build (Astro prerenders in its own environment), so the alias applies.
  vite: {
    resolve: { alias: { '@nstarkman/research-store': fileURLToPath(new URL('./src/lib/store.js', import.meta.url)) } },
    ssr: { noExternal: [/^@nstarkman\//] },
    environments: { prerender: { resolve: { noExternal: [/^@nstarkman\//] } } },
  },
  integrations: [mdx()],
  markdown: { remarkPlugins: [remarkMath], rehypePlugins: [[rehypeKatex, { macros, throwOnError: true }]] },
});
