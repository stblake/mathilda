// check-search.mjs — behavioural checks for notebook search.
//
// Compiles src/lib/search.ts with the project's own tsc and imports the RESULT,
// so these run the shipped functions rather than a paraphrase of them. Covers the
// pure half only: matching and index stepping. The DOM half (revealing a match in
// CodeMirror or in a prose cell) needs a real pointer in the WKWebView.
//
//     npm run check:search

import { execFileSync } from 'child_process';
import { mkdtempSync, rmSync } from 'fs';
import { join } from 'path';
import { pathToFileURL } from 'url';

/* Emitted INSIDE the project, not in the system temp dir: search.ts imports
   svelte/store, and node resolves bare specifiers by walking up from the importing
   file -- from /tmp there is no node_modules to find. */
const out = mkdtempSync(join('node_modules', '.mathilda-search-'));
try {
  execFileSync('node_modules/.bin/tsc', [
    'src/lib/search.ts',
    '--target', 'es2020', '--module', 'esnext', '--moduleResolution', 'bundler',
    '--skipLibCheck', '--ignoreConfig', '--outDir', out,
  ], { stdio: 'inherit' });
} catch {
  console.error('tsc failed to compile src/lib/search.ts');
  process.exit(1);
}

const { findMatches, stepIndex } = await import(pathToFileURL(join(out, 'search.js')).href);

let fail = 0;
const ok = (name, cond, got) => {
  console.log(`${cond ? 'PASS' : 'FAIL'}  ${name}${cond ? '' : `   got: ${JSON.stringify(got)}`}`);
  if (!cond) fail++;
};
const eq = (name, got, want) =>
  ok(`${name}`, JSON.stringify(got) === JSON.stringify(want), got);

const cells = [
  { id: 'c1', source: 'Sin[x] + sin[y]' },
  { id: 'c2', source: 'Table[i, {i, 4}]' },
  { id: 'c3', source: 'aaaa' },
];

eq('empty query finds nothing', findMatches(cells, '', false), []);

// NON-OVERLAPPING is the property, and "aa" in "aaaa" is the case that separates
// it from the naive from = at + 1 loop, which would give three matches.
eq('"aa" in "aaaa" gives 2 non-overlapping matches',
   findMatches([{ id: 'c3', source: 'aaaa' }], 'aa', false),
   [{ cellId: 'c3', start: 0, end: 2 }, { cellId: 'c3', start: 2, end: 4 }]);

// Case folding, both ways, on a source that contains both spellings.
eq('case-insensitive finds both spellings',
   findMatches([cells[0]], 'sin', false).map(m => m.start), [0, 9]);
eq('case-sensitive finds only the exact one',
   findMatches([cells[0]], 'sin', true).map(m => m.start), [9]);
eq('case-sensitive finds only the other exact one',
   findMatches([cells[0]], 'Sin', true).map(m => m.start), [0]);

// Offsets must be usable directly as source offsets -- that is what the CodeMirror
// dispatch relies on.
const m0 = findMatches([cells[0]], 'sin', true)[0];
ok('offsets slice out the match', cells[0].source.slice(m0.start, m0.end) === 'sin',
   cells[0].source.slice(m0.start, m0.end));

// Document order: cells in notebook order, offsets ascending within a cell.
const multi = findMatches(cells, 'a', false);
ok('matches come in document order',
   multi.every((m, i) => i === 0 || cellRank(m) > cellRank(multi[i - 1]) ||
                         (m.cellId === multi[i - 1].cellId && m.start > multi[i - 1].start)),
   multi.map(m => `${m.cellId}:${m.start}`));
function cellRank(m) { return cells.findIndex(c => c.id === m.cellId); }

ok('every match carries the id of the cell it is in',
   multi.every(m => cells.find(c => c.id === m.cellId)?.source
                         .slice(m.start, m.end).toLowerCase() === 'a'),
   multi);

// A single space is a real search; only the empty string is rejected.
ok('a space is a legitimate query',
   findMatches([cells[1]], ' ', false).length === 2,
   findMatches([cells[1]], ' ', false).length);

ok('a query longer than the source finds nothing',
   findMatches([{ id: 'x', source: 'ab' }], 'abc', false).length === 0, 'found some');

// stepIndex: the negative-modulo trap. In JavaScript -1 % 3 is -1, so the naive
// form sends Shift+Enter at the first match to a negative index.
eq('stepping back from the first match wraps to the last', stepIndex(0, -1, 3), 2);
eq('stepping forward from the last wraps to the first', stepIndex(2, 1, 3), 0);
eq('stepping forward in the middle', stepIndex(0, 1, 3), 1);
eq('no matches stays at 0', stepIndex(0, 1, 0), 0);
eq('a single match stays put going forward', stepIndex(0, 1, 1), 0);
eq('a single match stays put going back', stepIndex(0, -1, 1), 0);
// Walking a full cycle in each direction must return to the start exactly.
let i = 0; for (let k = 0; k < 7; k++) i = stepIndex(i, 1, 7);
eq('a full forward cycle returns to the start', i, 0);
let j = 0; for (let k = 0; k < 7; k++) j = stepIndex(j, -1, 7);
eq('a full backward cycle returns to the start', j, 0);

rmSync(out, { recursive: true, force: true });
console.log(fail === 0 ? '\nall search checks passed' : `\n${fail} FAILED`);
process.exit(fail === 0 ? 0 : 1);
