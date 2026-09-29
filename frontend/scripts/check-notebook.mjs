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
const { createNotebook, CELL_STYLES, isHeading, headingTag } = await imp('notebook.js');
const { convertCell, appendCellIfLast, previousInputSource } = await imp('cellCommands.js');
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
for (const [mid, t] of [['cell-toInput', 'code'], ['cell-toText', 'text'], ['cell-toSection', 'section']]) {
  const re = new RegExp(`case '${mid}':\\s*convertActiveCell\\('${t}'\\)`);
  ok(`menu ${mid} calls convertActiveCell('${t}')`, re.test(menuSrc), 'not wired');
}
ok('convertActiveCell converts in the store before relabelling the toolbar',
   /convertCell\(act\.store, cell\.cellId, type\)\)\s*retypeActiveCell\(type\)/.test(menuSrc),
   'not found');
const toolbarSrc = readFileSync('src/lib/Toolbar.svelte', 'utf8');
ok('the toolbar cell-style control uses the same convertCell',
   /convertCell\(pane\.store, activeCellObj\.id, type\)/.test(toolbarSrc), 'not found');

/* ---- 1b. Evaluating the last cell leaves a fresh one below -------------------
   The runaway cases are the ones worth pinning: an empty cell must not append
   (or holding Shift+Enter fills the notebook), and a cell that is not last must
   not, however many cells sit beside it in its own row. */

const ap = createNotebook();
const apFirst = ap.allCells()[0].id;
ok('an empty last cell appends nothing', appendCellIfLast(ap, apFirst) === null, 'appended');

ap.updateSource(apFirst, 'Integrate[x, x]');
const appended = appendCellIfLast(ap, apFirst);
ok('a non-empty last code cell appends one code cell', typeof appended === 'string', appended);
ok('the appended cell is code and empty',
   ap.allCells().length === 2 && ap.allCells()[1].type === 'code' && ap.allCells()[1].source === '',
   ap.allCells().map((c) => [c.type, c.source]));
ok('the appended cell is the one whose id came back', ap.allCells()[1].id === appended,
   ap.allCells()[1].id);

ok('the same cell does not append twice (it is no longer last)',
   appendCellIfLast(ap, apFirst) === null && ap.allCells().length === 2, ap.allCells().length);
ok('the new empty cell appends nothing either — no runaway',
   appendCellIfLast(ap, appended) === null && ap.allCells().length === 2, ap.allCells().length);

/* A sibling BESIDE the cell is not below it: a side-by-side row is still the last row. */
const sbs = createNotebook();
const sbsId = sbs.allCells()[0].id;
sbs.updateSource(sbsId, 'x + 1');
sbs.insertCellInRow(sbs.getRows()[0].id, 1, 'code', 'y + 1');
ok('a cell with a side-by-side sibling still counts as last',
   typeof appendCellIfLast(sbs, sbsId) === 'string' && sbs.getRows().length === 2,
   sbs.getRows().length);

/* Only code evaluates, so only code continues. */
const tx = createNotebook();
const txId = tx.allCells()[0].id;
tx.updateSource(txId, 'Some prose.');
convertCell(tx, txId, 'text');
ok('a text cell appends nothing', appendCellIfLast(tx, txId) === null, 'appended');
ok('an unknown cell id appends nothing', appendCellIfLast(tx, 'nope') === null, 'appended');

/* The two single-cell run paths call it; runAll/runRange must NOT, or "Evaluate
   Notebook" would end by appending a stray empty cell every time. */
const cardSrc = readFileSync('src/lib/NotebookCard.svelte', 'utf8');
const body = (name) => {
  const i = cardSrc.indexOf(`function ${name}(`);
  return i < 0 ? '' : cardSrc.slice(i, cardSrc.indexOf('\n  }', i));
};
ok('handleRun continues below', /continueBelowIfLast\(/.test(body('handleRun')), 'not wired');
ok('runCellById continues below', /continueBelowIfLast\(/.test(body('runCellById')), 'not wired');
ok('runAll does NOT append', !/continueBelowIfLast\(/.test(body('runAll')), 'appends');
ok('runRange does NOT append', !/continueBelowIfLast\(/.test(body('runRange')), 'appends');

/* Mod-Enter inserts a cell itself, so it must insert BEFORE it runs — otherwise
   run's own append fires too and the last cell yields two empty cells. */
const shellSrc = readFileSync('src/lib/CellShell.svelte', 'utf8');
/* Up to `return true`, not to the first `}` — the body contains `{ rowId }`. */
const modEnter = shellSrc.match(/key: 'Mod-Enter',\s*run\(\)\s*\{([\s\S]*?)return true;/);
ok('Mod-Enter dispatches addBelow before run',
   !!modEnter && modEnter[1].indexOf('addBelow') < modEnter[1].indexOf("'run'"),
   modEnter?.[1]?.trim());

/* ---- 1c. Output wraps rather than scrolling right --------------------------- */

const outSrc = readFileSync('src/lib/Output.svelte', 'utf8');
const rule = (sel) => {
  const i = outSrc.indexOf(`  .${sel} {`);
  return i < 0 ? '' : outSrc.slice(i, outSrc.indexOf('}', i));
};
for (const sel of ['out-error', 'out-expected', 'out-stream', 'out-message']) {
  ok(`.${sel} wraps instead of scrolling`, /white-space:\s*pre-wrap/.test(rule(sel)), rule(sel));
}
ok('.out-error no longer scrolls horizontally', !/overflow-x:\s*auto/.test(rule('out-error')),
   rule('out-error'));
ok('a typeset expression too wide for the card falls back to wrapping code',
   /renderOutput\(item\.text, item\.latex, wideExpr\[idx\]\)/.test(outSrc), 'not wired');
ok('the width fallback is measured, not guessed from string length',
   /scrollWidth > w \+ 1/.test(outSrc), 'not measured');
ok('kernel text is escaped before reaching {@html}',
   /replace\(\/</.test(outSrc) && /\$\{esc\(/.test(outSrc), 'unescaped');

/* ---- 1d. Cell styles ---------------------------------------------------------
   The lossy failure is the front end offering a style that the Rust `.mathilda`
   writer does not know: normalise_type runs on SERIALIZE too, so the heading is
   written out as `code` and destroyed by a save-and-reload rather than merely
   misread. That is what this diff is for. */

const styleIds = CELL_STYLES.map((s) => s.id);

const rustSrc = readFileSync('src-tauri/src/notebook_format.rs', 'utf8');
const knownSeg = rustSrc.slice(rustSrc.indexOf('const KNOWN_TYPES'), rustSrc.indexOf('];', rustSrc.indexOf('const KNOWN_TYPES')));
const rustTypes = [...knownSeg.matchAll(/"([a-z]+)"/g)].map((m) => m[1]);

/* `ref` is generated, never chosen, so it is in the Rust list and not in CELL_STYLES. */
const missingInRust = styleIds.filter((id) => !rustTypes.includes(id));
ok('every cell style survives a .mathilda save (present in Rust KNOWN_TYPES)',
   missingInRust.length === 0, missingInRust);
const strayInRust = rustTypes.filter((t) => t !== 'ref' && !styleIds.includes(t));
ok('Rust knows no type the front end cannot produce', strayInRust.length === 0, strayInRust);
ok('ref is known to Rust but is not an offered style',
   rustTypes.includes('ref') && !styleIds.includes('ref'), { rustTypes, styleIds });

for (const want of ['title', 'subtitle', 'chapter', 'section', 'subsection', 'subsubsection']) {
  ok(`${want} is an offered style and a heading`,
     styleIds.includes(want) && isHeading(want) === true, styleIds);
}
ok('code and text are offered and are not headings',
   isHeading('code') === false && isHeading('text') === false, 'heading');
ok('ref is not a heading', isHeading('ref') === false, 'heading');
ok('every heading style names a real heading element',
   CELL_STYLES.filter((s) => s.tag).every((s) => /^h[1-6]$/.test(s.tag)),
   CELL_STYLES.map((s) => [s.id, s.tag]));
ok('headingTag is null for a non-heading', headingTag('code') === null, headingTag('code'));
ok('styles are listed in outline order',
   styleIds.join() === 'code,text,title,subtitle,chapter,section,subsection,subsubsection',
   styleIds.join());
ok('every style has a label and a description',
   CELL_STYLES.every((s) => s.label && s.desc), CELL_STYLES.map((s) => s.id));

/* The store must accept each of them, and convertCell must actually retype. */
const styleStore = createNotebook();
const styleId = styleStore.allCells()[0].id;
for (const s of styleIds.filter((i) => i !== 'code')) {
  ok(`convertCell retypes to ${s}`,
     convertCell(styleStore, styleId, s) === true && styleStore.allCells()[0].type === s,
     styleStore.allCells()[0].type);
}

/* Each style needs its own CSS rule, or two levels render identically — which is
   exactly what a tag-keyed ladder did to Title and Section, both <h1>. */
const shellCss = readFileSync('src/lib/CellShell.svelte', 'utf8');
for (const s of CELL_STYLES.filter((x) => x.tag)) {
  ok(`.heading-${s.id} has its own rule`,
     new RegExp(`\\.heading-${s.id}\\s*\\{`).test(shellCss), 'no rule');
}
ok('the heading ladder is no longer keyed on the tag',
   !/h1\.heading-cell\s*\{/.test(shellCss) && !/h2\.heading-cell\s*\{/.test(shellCss),
   'still tag-keyed');
ok('one heading branch renders every level',
   /<svelte:element\s+this=\{headingTag\(cell\.type\)\}/.test(shellCss), 'not collapsed');

/* Every style reachable from the native menu, and the menu ids legal (see
   tools/check_menu_ids.py, which owns the grammar). */
for (const s of styleIds) {
  const cap = s.charAt(0).toUpperCase() + s.slice(1);
  const mid = s === 'code' ? 'cell-toInput' : `cell-to${cap}`;
  ok(`menu ${mid} exists for style ${s}`, menuSrc.includes(`'${mid}'`), 'missing');
}

/* ---- 1e. Cmd+L — copy input from above -------------------------------------- */

const ci = createNotebook();
const ciFirst = ci.allCells()[0].id;
ok('no input above the first cell', previousInputSource(ci, ciFirst) === null,
   previousInputSource(ci, ciFirst));

ci.updateSource(ciFirst, 'D[Sin[x], x]');
const ciSecond = ci.addRow('code');
ok('the cell above is found', previousInputSource(ci, ciSecond) === 'D[Sin[x], x]',
   previousInputSource(ci, ciSecond));

/* An empty cell in between must be skipped — the cell that evaluating the last
   cell appends is empty, and it would otherwise always be "the input above". */
const ciThird = ci.addRow('code');
ok('an empty cell in between is skipped', previousInputSource(ci, ciThird) === 'D[Sin[x], x]',
   previousInputSource(ci, ciThird));

/* So must a heading or prose cell: "input" means a code cell. */
const ciHead = ci.addRow('section', 'A heading');
const ciAfter = ci.addRow('code');
ok('a heading is not input', previousInputSource(ci, ciAfter) === 'D[Sin[x], x]',
   previousInputSource(ci, ciAfter));
ok('a heading cell itself has the code above it',
   previousInputSource(ci, ciHead) === 'D[Sin[x], x]', previousInputSource(ci, ciHead));
ok('an unknown cell id has no input above', previousInputSource(ci, 'nope') === null, 'found');

/* Multi-line input comes back whole. */
const ml = createNotebook();
const mlFirst = ml.allCells()[0].id;
ml.updateSource(mlFirst, 'a = 1\nb = 2\na + b');
const mlNext = ml.addRow('code');
ok('multi-line input is copied whole',
   previousInputSource(ml, mlNext) === 'a = 1\nb = 2\na + b', previousInputSource(ml, mlNext));

ok('Mod-l is bound in the cell keymap', /key: 'Mod-l'/.test(shellSrc), 'unbound');
ok('Mod-l is bound before defaultKeymap so it wins',
   shellSrc.indexOf("key: 'Mod-l'") < shellSrc.indexOf('...defaultKeymap'), 'after defaultKeymap');
ok('Edit > Copy Input from Above is wired', menuSrc.includes("'edit-copyInputAbove'"), 'missing');

/* ---- 1f. The canvas opens with one empty notebook ---------------------------- */

const canvasSrc = readFileSync('src/lib/canvas.ts', 'utf8');
ok('canvasState is seeded with exactly one notebook',
   /notebooks:\s*\[_nb1\]\s*as CanvasNotebook\[\]/.test(canvasSrc), 'more than one');
ok('that notebook has no title of its own, so it reads as File > New',
   /const _nb1 = makeCard\(''/.test(canvasSrc), 'titled');
/* Comment-blanked before the test, like tools/check_message_routing.py: the comment that replaced
   the call NAMES it, so a raw grep matches the explanation and reports the call it documents. */
const decomment = (s) => s.replace(/\/\*[\s\S]*?\*\//g, '').replace(/(^|[^:])\/\/.*$/gm, '$1');
const canvasComp = decomment(readFileSync('src/lib/Canvas.svelte', 'utf8'));
ok('no tour is loaded on mount', !/loadStartupContent\(\)/.test(canvasComp), 'still loads');
ok('both tours survive as opt-in calls',
   /export function loadStartupContent/.test(canvasSrc) &&
   /export function loadDemoContent/.test(canvasSrc), 'a tour was deleted');
ok('a tour builds its own cards rather than assuming nine exist',
   /function ensureStarterCards/.test(canvasSrc) &&
   /ensureStarterCards\(\)/.test(canvasSrc.slice(canvasSrc.indexOf('loadDemoContent'))),
   'still indexes a nine-card array');

/* ---- 1g. In[n] comes from the KERNEL, so `%n` addresses what the label says ---
   `%`, `%%`, `%3`, In[3] and Out[3] all resolve against the kernel's $Line
   (src/repl.c; behaviour gated by tools/check_pipe_protocol.py). One kernel
   serves every notebook on the canvas, so a per-notebook counter would label a
   cell In[2] while `%2` addressed a line typed in another pane. The local stamp
   stays as the placeholder that appears the instant Shift+Enter is pressed. */

const ex = createNotebook();
const exId = ex.allCells()[0].id;
ok('stampExec still numbers locally for the instant label', ex.stampExec(exId) === 1, 'no stamp');
ex.setExec(exId, 7);
ok('setExec replaces the label with the kernel line', ex.allCells()[0].execIdx === 7,
   ex.allCells()[0].execIdx);
ok('setExec on an unknown id changes nothing',
   (ex.setExec('nope', 99), ex.allCells()[0].execIdx === 7), ex.allCells()[0].execIdx);

const ipcSrc = readFileSync('src/lib/ipc.ts', 'utf8');
ok('the protocol declares the line message',
   /type:\s*'line';\s*line:\s*number/.test(ipcSrc), 'undeclared');
const runCellBody = body('runCell');
ok('runCell reconciles the label from the first line message',
   /msg\.type === 'line'/.test(runCellBody) && /setExec\(cellId, msg\.line\)/.test(runCellBody),
   'not wired');
ok('only the FIRST line labels the cell — a multi-statement cell takes several',
   /kernelLine === null/.test(runCellBody), 'relabels on every statement');
ok('the timing record uses the kernel line too',
   /label: `In\[\$\{kernelLine \?\? execIdx\}\]`/.test(cardSrc), 'local counter');
ok('a line message renders no output of its own',
   /case 'line':\s*return null;/.test(cardSrc), 'shown as output');

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
