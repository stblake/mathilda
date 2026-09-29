// check-notebook.mjs — behavioural checks for three notebook fixes.
//
//   1. Cell > Convert to ... : convertCell really retypes the cell in the store
//      (the menu used to update only the toolbar's record of the active cell),
//      and the native menu items route through it.
//   2. Stacking: the focused-mode overlays (properties panel, find bar) sit
//      ABOVE the focused notebook view on the one z-index scale in app.css.
//   3. Find bar: the current-match mark lives in CodeMirror state -- set, moved,
//      mapped through edits, clamped, cleared -- so the find field can keep focus.
//
// Like check:search, the TypeScript is compiled with the project's own tsc and
// the RESULT is imported, so these run the shipped functions, not a paraphrase.
//
//     npm run check:notebook

import { execFileSync } from 'child_process';
import { mkdtempSync, rmSync, readFileSync } from 'fs';
import { join } from 'path';
import { pathToFileURL } from 'url';

/* Emitted inside the project so bare imports (svelte/store, @codemirror/*)
   resolve from node_modules. */
const out = mkdtempSync(join('node_modules', '.mathilda-notebook-'));
try {
  execFileSync('node_modules/.bin/tsc', [
    'src/lib/notebook.ts', 'src/lib/cellCommands.ts', 'src/lib/searchHighlight.ts',
    '--target', 'es2020', '--module', 'esnext', '--moduleResolution', 'bundler',
    '--skipLibCheck', '--ignoreConfig', '--outDir', out,
  ], { stdio: 'inherit' });
} catch {
  console.error('tsc failed to compile the notebook modules');
  process.exit(1);
}

const imp = (f) => import(pathToFileURL(join(out, f)).href);
const { createNotebook } = await imp('notebook.js');
const { convertCell } = await imp('cellCommands.js');
const { searchMarkExtension, setSearchMark } = await imp('searchHighlight.js');
const { EditorState } = await import('@codemirror/state');
const { EditorView } = await import('@codemirror/view');

let fail = 0;
const ok = (name, cond, got) => {
  console.log(`${cond ? 'PASS' : 'FAIL'}  ${name}${cond ? '' : `   got: ${JSON.stringify(got)}`}`);
  if (!cond) fail++;
};

/* ---- 1. Convert to ---------------------------------------------------------- */

const store = createNotebook();
const id = store.allCells()[0].id;
store.updateSource(id, 'Integrate[x, x]');
store.appendOutput(id, { kind: 'expr', text: 'x^2/2' });

ok('convertCell code -> text changes the cell in the store',
   convertCell(store, id, 'text') === true && store.allCells()[0].type === 'text',
   store.allCells()[0].type);
ok('the source survives the conversion', store.allCells()[0].source === 'Integrate[x, x]',
   store.allCells()[0].source);
ok('text -> section', convertCell(store, id, 'section') && store.allCells()[0].type === 'section',
   store.allCells()[0].type);
ok('converting to the current type is a no-op that reports false',
   convertCell(store, id, 'section') === false, 'true');
ok('section -> code restores the kept output',
   convertCell(store, id, 'code') && store.allCells()[0].output.length === 1,
   store.allCells()[0].output);
ok('an unknown cell id reports false', convertCell(store, 'nope', 'text') === false, 'true');

const refStore = createNotebook();
refStore.setCellSourceAndType('Sin', 'ref');
const refId = refStore.allCells()[0].id;
ok('a reference-page cell is not convertible',
   convertCell(refStore, refId, 'code') === false && refStore.allCells()[0].type === 'ref',
   refStore.allCells()[0].type);

/* The native menu's three items must go through the store conversion, not only
   through retypeActiveCell (which relabels the toolbar and nothing else). */
const menuSrc = readFileSync('src/lib/menuCommands.ts', 'utf8');
for (const [mid, t] of [['cell.toInput', 'code'], ['cell.toText', 'text'], ['cell.toSection', 'section']]) {
  const re = new RegExp(`case '${mid.replace('.', '\\.')}':\\s*convertActiveCell\\('${t}'\\)`);
  ok(`menu ${mid} calls convertActiveCell('${t}')`, re.test(menuSrc), 'not wired');
}
ok('convertActiveCell converts in the store before relabelling the toolbar',
   /convertCell\(act\.store, cell\.cellId, type\)\)\s*retypeActiveCell\(type\)/.test(menuSrc),
   'not found');
const toolbarSrc = readFileSync('src/lib/Toolbar.svelte', 'utf8');
ok('the toolbar cell-style control uses the same convertCell',
   /convertCell\(pane\.store, activeCellObj\.id, type\)/.test(toolbarSrc), 'not found');

/* ---- 2. Stacking ------------------------------------------------------------ */

const css = readFileSync('src/app.css', 'utf8');
const z = (name) => {
  const m = css.match(new RegExp(`--z-${name}:\\s*(-?\\d+)`));
  return m ? Number(m[1]) : NaN;
};
ok('properties panel sits above the focused notebook view',
   z('overlay') > z('focused-view'), { overlay: z('overlay'), view: z('focused-view') });
ok('find bar sits above the properties panel', z('overlay-top') > z('overlay'),
   { top: z('overlay-top'), overlay: z('overlay') });
ok('overlays sit above the status dock', z('overlay') > z('status-dock'),
   { overlay: z('overlay'), dock: z('status-dock') });
ok('menus sit above everything else', z('menu') > Math.max(z('appbar'), z('banner'), z('overlay-top')),
   { menu: z('menu') });

const usesVar = (file, v) => new RegExp(`z-index:\\s*var\\(--z-${v}\\)`).test(readFileSync(file, 'utf8'));
ok('PropertiesPanel takes its z-index from the scale', usesVar('src/lib/PropertiesPanel.svelte', 'overlay'), 'raw number');
ok('SearchBar takes its z-index from the scale', usesVar('src/lib/SearchBar.svelte', 'overlay-top'), 'raw number');
ok('the focused view takes its z-index from the scale', usesVar('src/lib/Canvas.svelte', 'focused-view'), 'raw number');

/* ---- 3. The find bar's current-match mark ---------------------------------- */

const marks = (state) => {
  const found = [];
  for (const d of state.facet(EditorView.decorations)) {
    const set = typeof d === 'function' ? null : d;
    if (!set) continue;
    for (const c = set.iter(); c.value; c.next()) found.push([c.from, c.to]);
  }
  return found;
};

let st = EditorState.create({ doc: 'Sin[x] + Cos[x]', extensions: [searchMarkExtension] });
ok('no mark initially', marks(st).length === 0, marks(st));
st = st.update({ effects: setSearchMark.of({ start: 9, end: 12 }) }).state;
ok('setSearchMark paints the range', JSON.stringify(marks(st)) === '[[9,12]]', marks(st));
st = st.update({ effects: setSearchMark.of({ start: 0, end: 3 }) }).state;
ok('a new mark replaces the old one (one current match)', JSON.stringify(marks(st)) === '[[0,3]]', marks(st));
st = st.update({ changes: { from: 0, insert: '2 ' } }).state;
ok('the mark maps through an edit before it', JSON.stringify(marks(st)) === '[[2,5]]', marks(st));
st = st.update({ effects: setSearchMark.of({ start: 10, end: 999 }) }).state;
ok('a range past the end is clamped to the document',
   JSON.stringify(marks(st)) === `[[10,${st.doc.length}]]`, marks(st));
st = st.update({ effects: setSearchMark.of(null) }).state;
ok('null clears the mark', marks(st).length === 0, marks(st));
st = st.update({ effects: setSearchMark.of({ start: 4, end: 4 }) }).state;
ok('an empty range paints nothing', marks(st).length === 0, marks(st));

rmSync(out, { recursive: true, force: true });
console.log(fail === 0 ? '\nall notebook checks passed' : `\n${fail} FAILED`);
process.exit(fail === 0 ? 0 : 1);
