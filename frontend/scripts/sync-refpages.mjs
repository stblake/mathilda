/*
 * Copy the generated reference pages into public/ so the app can fetch one on
 * demand, and build the name -> path index it needs to do so.
 *
 * The pages are produced by `python3 site/generate.py`, which mines
 * docs/spec/builtins/<slug>.md and RE-VERIFIES every example against the built
 * ./Mathilda binary. This script does not author anything; it only mirrors that
 * output. Re-run it after regenerating the site, or a symbol added since the
 * last sync will 404 in the app.
 *
 * Why an index rather than a flat <Name>.md directory: symbol names are NOT
 * unique across categories. Det, Factor, Inverse, LinearSolve, PolyGamma and a
 * few others are documented under both their own category and `flint`, so a
 * flat copy would have one silently overwrite the other. The index resolves a
 * name to exactly one page, preferring the category that builtins.json assigns
 * to the symbol -- i.e. the same page the published site links to.
 */

import { readFileSync, writeFileSync, mkdirSync, copyFileSync, readdirSync, statSync, rmSync } from 'node:fs';
import { join, dirname, basename } from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = dirname(fileURLToPath(import.meta.url));
const REPO = join(HERE, '..', '..');
const SRC = join(REPO, 'site', 'docs', 'documentation');
const DEST = join(HERE, '..', 'public', 'refpages');
const BUILTINS = join(REPO, 'site', 'docs', 'assets', 'builtins.json');

/* Preferred category per symbol, straight from the site's own manifest. */
const preferred = new Map();
try {
  for (const b of JSON.parse(readFileSync(BUILTINS, 'utf8'))) {
    if (b.name && b.slug) preferred.set(b.name, b.slug);
  }
} catch {
  console.warn('sync-refpages: no builtins.json; falling back to first-seen category');
}

function walk(dir, ext = '.md') {
  const out = [];
  for (const entry of readdirSync(dir)) {
    const p = join(dir, entry);
    if (statSync(p).isDirectory()) out.push(...walk(p, ext));
    else if (entry.endsWith(ext) && entry !== 'index.md') out.push(p);
  }
  return out;
}

const pages = walk(SRC);

/* name -> slug/Name.md, resolving the cross-category duplicates. */
const index = {};
const chosenSlug = {};
for (const p of pages) {
  const name = basename(p, '.md');
  const slug = basename(dirname(p));
  const want = preferred.get(name);
  const better =
    !(name in index) ||                       /* first one wins by default */
    (want && slug === want) ||                /* the manifest's category always wins */
    (!want && slug < chosenSlug[name]);       /* else deterministic, not filesystem order */
  if (better) {
    index[name] = `${slug}/${name}.md`;
    chosenSlug[name] = slug;
  }
}

/* Rebuild from scratch so a page deleted upstream does not linger and keep
   answering fetches with stale content. */
rmSync(DEST, { recursive: true, force: true });
let copied = 0;
for (const [name, rel] of Object.entries(index)) {
  const dst = join(DEST, rel);
  mkdirSync(dirname(dst), { recursive: true });
  copyFileSync(join(SRC, rel), dst);
  copied++;
}
writeFileSync(join(DEST, 'index.json'), JSON.stringify(index));

/* Saved figures, so a graphics example shows its plot before anything is run.
   Same <slug>/<Name> layout as the pages, so a page's figures are found by
   swapping the extension. */
const PLOTS_SRC = join(REPO, 'site', 'docs', 'assets', 'plots');
let figures = 0;
try {
  for (const p of walk(PLOTS_SRC, '.json')) {
    const rel = p.slice(PLOTS_SRC.length + 1);
    const dst = join(DEST, 'plots', rel);
    mkdirSync(dirname(dst), { recursive: true });
    copyFileSync(p, dst);
    figures++;
  }
} catch { /* no figures generated yet */ }

const dupes = pages.length - copied;
console.log(`sync-refpages: ${copied} pages` +
            (figures ? `, ${figures} figure set(s)` : '') +
            ` -> public/refpages` +
            (dupes > 0 ? ` (${dupes} cross-category duplicate(s) resolved)` : ''));
