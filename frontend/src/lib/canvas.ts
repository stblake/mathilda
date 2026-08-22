// canvas.ts — infinite canvas state: pan, zoom, and named notebooks

import { writable, derived, get } from 'svelte/store';

// Re-export the notebook factory so each card gets its own store instance.
export { createNotebook } from './notebook';
import { createNotebook } from './notebook';
import { fetchRefpageMarkdown, splitRefpage, buildToc, fetchRefpageFigures } from './refpages';
import { clearActiveCell } from './active';

export type NotebookStore = ReturnType<typeof createNotebook>;

export type CanvasNotebook = {
  id:        string;
  title:     string;
  x:         number;
  y:         number;
  width:     number;   // card width in canvas pixels
  height:    number | null; // card height override; null = auto
  collapsed: boolean;
  store:     NotebookStore;
  /** Set when the card was opened to answer a query (e.g. clicking a name in a
   *  `?pat*` result). NotebookCard runs its first cell once on mount and clears
   *  the flag, so the answer is already there when the card appears. */
  autoRun?:  boolean;
  /** A generated reference page rather than a user notebook. Its headings are
   *  documentation structure, not editable prose, so clicking one folds the
   *  section instead of placing a caret. */
  refpage?:  boolean;
};

let _nextId = 1;
const PALETTE = ['#4a90e2','#e67e22','#27ae60','#8e44ad','#e74c3c','#16a085'];

function makeCard(title: string, x: number, y: number): CanvasNotebook {
  // CRITICAL: always capture and increment _nextId so every card gets a unique id.
  // The old pattern `title || \`Notebook ${_nextId++}\`` never incremented when
  // title was non-empty, causing every subsequent makeCard('', ...) to get id="nb-1"
  // and Svelte's keyed {#each} to collapse them into one rendered card.
  const n = _nextId++;
  return {
    id: `nb-${n}`,
    title: title || `Notebook ${n}`,
    x, y,
    width: 640,
    height: null,
    collapsed: false,
    store: createNotebook(),
  };
}

// 8 starter notebooks in 3 clusters (loadStartupContent fills them)
// Cluster 1 — Calculus (left)
const _nb1 = makeCard('Derivatives',        50,  50);
const _nb2 = makeCard('Integration',       760,  50);
const _nb3 = makeCard('Function Plots',     50, 720);
// Cluster 2 — Algebra & Numbers (right)
const _nb4 = makeCard('Polynomial Algebra',1620,  50);
const _nb5 = makeCard('Number Theory',     1620, 720);
const _nb6 = makeCard('Linear Algebra',    2290,  50);
// Cluster 3 — Special Topics (bottom center)
const _nb7 = makeCard('Special Functions',  750,1650);
const _nb8 = makeCard('Applied Math',      1430,1650);
const _nb9 = makeCard('Associations',      2110,1650);

/* ---------------------------------------------------------------------------
 * Pane registries
 *
 * A card's controls live in the card (its run loop, its layout flag, its rename
 * state) but in focused mode the buttons belong in the top bar rather than in a
 * second bar directly underneath it. So each card publishes its controls here
 * and the toolbar reads them.
 *
 * This used to be a single `focusedActions` writable holding one record. That
 * only worked because exactly one notebook could be focused at a time; with
 * panes, N cards would fight over one slot. Now it is keyed by notebook id and
 * each card owns exactly its own key.
 *
 * Split into TWO maps on purpose. `panes` holds the stable method table,
 * registered once on mount. `paneFlags` holds the handful of values that change
 * constantly -- `allSectionsCollapsed` derives from the notebook's rows, so it
 * is invalidated on every keystroke. Publishing both together meant a keystroke
 * replaced the whole record and re-ran every toolbar button's subscriber; with
 * seven glyph buttons that was invisible, with labelled groups and a
 * validate-on-read cell lookup it would not be.
 * ------------------------------------------------------------------------- */

/** Stable per-card methods. Registered once on mount; never rebuilt. */
export interface PaneActions {
  notebookId: string;
  store: NotebookStore;
  runAll: () => Promise<void>;
  runCell: (cellId: string) => Promise<void>;
  /** Inclusive index range over store.allCells(); skips non-code and blank. */
  runRange: (fromIdx: number, toIdx: number) => Promise<void>;
  /** Move the caret to a cell. Wraps the card's private focus registry, which
   *  the toolbar has no prop path to. */
  focusCell: (cellId: string) => void;
  toggleLayout: () => void;
  rename: () => void;
  toggleAllSections: () => void;
  toggleCollapse: () => void;
  close: () => void;
}

/** The volatile half: values a keystroke can change. */
export interface PaneFlags {
  horizontal: boolean;
  allSectionsCollapsed: boolean;
  hasSections: boolean;
  collapsed: boolean;
}

export const panes     = writable<Map<string, PaneActions>>(new Map());
export const paneFlags = writable<Map<string, PaneFlags>>(new Map());

/* Which card instance currently owns each key. Plain Map, not a store -- it is
   bookkeeping, and nothing renders from it.

   This exists because a card can be destroyed and re-created for the same
   notebook id when the view switches between canvas and focused mode, and
   Svelte does not guarantee the old block is torn down before the new one
   mounts. Without an ownership token, a destroy running after the replacement
   card mounted would delete the entry that card had just registered, leaving a
   focused notebook whose toolbar buttons do nothing. */
const paneTokens = new Map<string, object>();

/** Publish a card's stable methods. `token` identifies the card instance. */
export function registerPane(id: string, token: object, actions: PaneActions) {
  paneTokens.set(id, token);
  panes.update(m => new Map(m).set(id, actions));
}

/** Publish a card's volatile flags. Ignored if a newer card owns the key. */
export function publishPaneFlags(id: string, token: object, flags: PaneFlags) {
  if (paneTokens.get(id) !== token) return;
  paneFlags.update(m => new Map(m).set(id, flags));
}

/** Withdraw a card's entries — unless a newer card has already taken the key. */
export function retractPane(id: string, token: object) {
  if (paneTokens.get(id) !== token) return;
  paneTokens.delete(id);
  panes.update(m     => { const n = new Map(m); n.delete(id); return n; });
  paneFlags.update(m => { const n = new Map(m); n.delete(id); return n; });
}

// ---------------------------------------------------------------------------
// Focus state

export type FocusLayout = 'h' | 'v' | 'grid';

/** Tiling more than four notebooks in one window is not usable at any realistic
 *  window size, and the grid layout has no meaning beyond four. */
export const MAX_PANES = 4;

/** The focus half of canvasState, as one literal.
 *
 *  Shared by the store's initial value AND loadLibraryData's full-object set().
 *  That set() is the most dangerous line in this file: it rebuilds the whole
 *  state from a literal, so a field omitted there silently becomes `undefined`
 *  and the next render throws on `focusedIds.length`. One literal, two callers,
 *  no way to forget a field. */
function initialFocus() {
  return {
    /** Notebooks filling the window. Empty = canvas mode. Order = pane order.
     *  Unique: a notebook may appear in at most one pane (see addPaneToFocus). */
    focusedIds:      [] as string[],
    focusedLayout:   'h' as FocusLayout,
    /** The pane the toolbar drives. Always a member of focusedIds, or null when
     *  focusedIds is empty. Deliberately separate from `activeId` below, which
     *  is canvas z-order: conflating them would make the toolbar's target change
     *  as a side effect of opening a reference page. */
    focusedActiveId: null as string | null,
    /** h/v: percent per pane along the flow axis. Length === focusedIds.length. */
    focusedSizes:    [] as number[],
    /** grid: the two divider positions, as percentages. Kept separate from
     *  focusedSizes because a flat array cannot express a 2x2 -- a grid needs two
     *  independent dividers -- and because keeping both means switching
     *  h -> grid -> h round-trips without losing either. */
    focusedGrid:     { x: 50, y: 50 },
  };
}

export const canvasState = writable({
  notebooks:    [_nb1,_nb2,_nb3,_nb4,_nb5,_nb6,_nb7,_nb8,_nb9] as CanvasNotebook[],
  panX:         0,
  panY:         0,
  zoom:         1.0,
  ...initialFocus(),
  /* The card drawn on top. Lives in the store rather than in Canvas.svelte
     because openRefpage has to raise the card it creates -- a documentation
     page that opens behind the notebook you asked from is invisible. */
  activeId:     null as string | null,
  selectedIds:  [] as string[],   // IDs of notebooks selected by rubber-band
});

/** The active pane's methods — what the toolbar's buttons call. */
export const activeActions = derived(
  [panes, canvasState],
  ([$panes, $s]) => ($s.focusedActiveId ? $panes.get($s.focusedActiveId) ?? null : null),
);

/** The active pane's volatile flags — what the toolbar's icons reflect. */
export const activeFlags = derived(
  [paneFlags, canvasState],
  ([$flags, $s]) => ($s.focusedActiveId ? $flags.get($s.focusedActiveId) ?? null : null),
);

/** True while any notebook fills the window. */
export const isFocused = derived(canvasState, s => s.focusedIds.length > 0);

/** Pre-fill the starter notebooks with rich example cells. Call from onMount. */
/* What the canvas opens with.
 *
 * The image subsystem is what is being built and read right now, so the canvas opens on its
 * documentation: nine reference pages, each a real notebook whose examples are live cells that
 * can be edited and re-run. `loadDemoContent` below is the previous tour (calculus, plots,
 * linear algebra) and is one call away -- swap the call in Canvas.svelte's mount. */
const STARTUP_DOCS = [
  'Image', 'Image3D', 'ImageConvolve',
  'GaussianFilter', 'EdgeDetect', 'Binarize',
  'Dilation', 'ImagePad', 'Import',
];

export function loadStartupContent() {
  const s = get(canvasState);
  /* Retitle the starter cards in place and mark them as reference pages. Mutating the existing
     cards rather than creating new ones keeps their laid-out positions, which is the whole
     reason the 3x3 arrangement reads as a grid. */
  STARTUP_DOCS.forEach((name, i) => {
    const nb = s.notebooks[i];
    if (!nb) return;
    nb.title = name;
    nb.refpage = true;
    nb.store.setCellSourceAndType(`Loading the reference page for \`${name}\`...`, 'ref');
    void fillRefpage(nb.store, name);
  });
  /* One update so every retitled card renders once, not nine times. */
  canvasState.update(st => ({ ...st, notebooks: [...st.notebooks] }));
}

export function loadDemoContent() {
  const s = get(canvasState);

  const [nb1,nb2,nb3,nb4,nb5,nb6,nb7,nb8,nb9] = s.notebooks;

  // Cluster 1 — Calculus
  if (nb1) nb1.store.load([
    { cells: [{ type: 'section', source: 'Derivatives' }] },
    { cells: [{ type: 'text',    source: 'Symbolic differentiation. Press Shift+Enter to evaluate.' }] },
    { cells: [{ type: 'code',    source: 'D[Sin[x]^2, x]' }] },
    { cells: [{ type: 'code',    source: 'D[x^x, x]' }] },
    { cells: [{ type: 'code',    source: 'D[Exp[x] Sin[x], x]' }] },
    { cells: [{ type: 'code',    source: 'D[Log[x^2 + 1], x]' }] },
    { cells: [{ type: 'section', source: 'Higher Order & Chain Rule' }] },
    { cells: [{ type: 'code',    source: 'D[Sin[x], {x, 4}]' }] },
    { cells: [{ type: 'code',    source: 'D[Sin[Exp[x]], x]' }] },
    { cells: [{ type: 'code',    source: 'D[Sqrt[1 + x^2], x]' }] },
  ]);

  if (nb2) nb2.store.load([
    { cells: [{ type: 'section', source: 'Definite Integration' }] },
    { cells: [{ type: 'text',    source: 'Compute exact areas and antiderivatives.' }] },
    { cells: [{ type: 'code',    source: 'Integrate[x^2, {x, 0, 1}]' }] },
    { cells: [{ type: 'code',    source: 'Integrate[Sin[x], {x, 0, Pi}]' }] },
    { cells: [{ type: 'code',    source: 'Integrate[Exp[-x^2], {x, -Infinity, Infinity}]' }] },
    { cells: [{ type: 'code',    source: 'Integrate[Log[x], x]' }] },
    { cells: [{ type: 'section', source: 'Series & Limits' }] },
    { cells: [{ type: 'code',    source: 'Series[Sin[x], {x, 0, 7}]' }] },
    { cells: [{ type: 'code',    source: 'Limit[Sin[x]/x, x -> 0]' }] },
    { cells: [{ type: 'code',    source: 'Limit[(1 + 1/n)^n, n -> Infinity]' }] },
  ]);

  if (nb3) nb3.store.load([
    { cells: [{ type: 'section', source: 'Function Plots' }] },
    { cells: [{ type: 'text',    source: 'Adaptive 2D plots. Multiple functions can be overlaid.' }] },
    { cells: [{ type: 'code',    source: 'Plot[Sin[x], {x, 0, 2 Pi}]' }] },
    { cells: [{ type: 'code',    source: 'Plot[{Sin[x], Cos[x]}, {x, 0, 2 Pi}]' }] },
    { cells: [{ type: 'section', source: 'Filled Plots' }] },
    { cells: [{ type: 'code',    source: 'Plot[Sin[x] + Sin[5 x], {x, 0, 4 Pi}, Filling -> Axis]' }] },
    { cells: [{ type: 'code',    source: 'Plot[Exp[-x^2/2]/Sqrt[2 Pi], {x, -4, 4}, Filling -> Axis]' }] },
    { cells: [{ type: 'section', source: '3D Plots' }] },
    { cells: [{ type: 'code',    source: 'Plot3D[Sin[x] Cos[y], {x, 0, 3}, {y, 0, 3}]' }] },
  ]);

  // Cluster 2 — Algebra & Numbers
  if (nb4) nb4.store.load([
    { cells: [{ type: 'section', source: 'Factoring & Expanding' }] },
    { cells: [{ type: 'text',    source: 'Polynomial manipulation.' }] },
    { cells: [{ type: 'code',    source: 'Factor[x^4 - 1]' }] },
    { cells: [{ type: 'code',    source: 'Factor[x^6 - y^6]' }] },
    { cells: [{ type: 'code',    source: 'Expand[(x + y)^5]' }] },
    { cells: [{ type: 'section', source: 'Solving Equations' }] },
    { cells: [{ type: 'code',    source: 'Solve[x^2 - 5 x + 6 == 0, x]' }] },
    { cells: [{ type: 'code',    source: 'Solve[{x + y == 5, x - y == 1}, {x, y}]' }] },
    { cells: [{ type: 'section', source: 'Simplification' }] },
    { cells: [{ type: 'code',    source: 'Simplify[(x^2 - 1)/(x - 1)]' }] },
    { cells: [{ type: 'code',    source: 'Together[1/x + 1/(x+1) + 1/(x+2)]' }] },
  ]);

  if (nb5) nb5.store.load([
    { cells: [{ type: 'section', source: 'Primes & Factoring' }] },
    { cells: [{ type: 'text',    source: 'Integer structure and divisibility.' }] },
    { cells: [{ type: 'code',    source: 'Select[Range[50], PrimeQ]' }] },
    { cells: [{ type: 'code',    source: 'FactorInteger[720720]' }] },
    { cells: [{ type: 'code',    source: 'NextPrime[1000]' }] },
    { cells: [{ type: 'section', source: 'Arithmetic Functions' }] },
    { cells: [{ type: 'code',    source: 'EulerPhi[100]' }] },
    { cells: [{ type: 'code',    source: 'Divisors[360]' }] },
    { cells: [{ type: 'code',    source: 'GCD[144, 89]' }] },
    { cells: [{ type: 'section', source: 'Digit Functions' }] },
    { cells: [{ type: 'code',    source: 'DigitSum[123456789]' }] },
    { cells: [{ type: 'code',    source: 'Table[DigitSum[2^n], {n, 1, 15}]' }] },
  ]);

  if (nb6) nb6.store.load([
    { cells: [{ type: 'section', source: 'Matrix Operations' }] },
    { cells: [{ type: 'text',    source: 'Matrices encode linear transformations.' }] },
    { cells: [{ type: 'code',    source: 'Det[{{1,2,3},{4,5,6},{7,8,10}}]' }] },
    { cells: [{ type: 'code',    source: 'Inverse[{{2,1},{1,3}}]' }] },
    { cells: [{ type: 'code',    source: 'MatrixPower[{{1,1},{1,0}}, 10]' }] },
    { cells: [{ type: 'section', source: 'Eigenvalues' }] },
    { cells: [{ type: 'code',    source: 'Eigenvalues[{{2,1},{1,2}}]' }] },
    { cells: [{ type: 'code',    source: 'Eigenvectors[{{3,1},{1,3}}]' }] },
    { cells: [{ type: 'section', source: 'Systems of Equations' }] },
    { cells: [{ type: 'code',    source: 'LinearSolve[{{1,2},{3,4}},{5,6}]' }] },
    { cells: [{ type: 'code',    source: 'NullSpace[{{1,2,3},{4,5,6},{7,8,9}}]' }] },
  ]);

  // Cluster 3 — Special Topics
  if (nb7) nb7.store.load([
    { cells: [{ type: 'section', source: 'Gamma & Zeta' }] },
    { cells: [{ type: 'text',    source: 'Higher transcendental functions.' }] },
    { cells: [{ type: 'code',    source: 'Gamma[5]' }] },
    { cells: [{ type: 'code',    source: 'Gamma[1/2]' }] },
    { cells: [{ type: 'code',    source: 'Zeta[2]' }] },
    { cells: [{ type: 'code',    source: 'Zeta[4]' }] },
    { cells: [{ type: 'code',    source: 'N[Pi, 50]' }] },
    { cells: [{ type: 'section', source: 'Combinatorics' }] },
    { cells: [{ type: 'code',    source: 'Table[Fibonacci[n], {n, 1, 15}]' }] },
    { cells: [{ type: 'code',    source: 'Table[Binomial[n, 2], {n, 1, 10}]' }] },
    { cells: [{ type: 'code',    source: 'LucasL[10]' }] },
  ]);

  if (nb8) nb8.store.load([
    { cells: [{ type: 'section', source: 'Sums & Series' }] },
    { cells: [{ type: 'text',    source: 'Symbolic summation and closed forms.' }] },
    { cells: [{ type: 'code',    source: 'Sum[k^2, {k, 1, n}]' }] },
    { cells: [{ type: 'code',    source: 'Sum[1/k^2, {k, 1, Infinity}]' }] },
    { cells: [{ type: 'code',    source: 'Series[1/(1-x), {x, 0, 8}]' }] },
    { cells: [{ type: 'section', source: 'Modular Arithmetic' }] },
    { cells: [{ type: 'code',    source: 'PowerMod[2, 100, 13]' }] },
    { cells: [{ type: 'code',    source: 'Table[Mod[2^n, 7], {n, 0, 12}]' }] },
    { cells: [{ type: 'section', source: 'Rational Functions' }] },
    { cells: [{ type: 'code',    source: 'Apart[1/(x^2 - 1)]' }] },
    { cells: [{ type: 'code',    source: 'Apart[(x^2+1)/((x-1)(x+2))]' }] },
  ]);

  if (nb9) nb9.store.load([
    { cells: [{ type: 'section', source: 'Associations' }] },
    { cells: [{ type: 'text',    source: 'Key-value data: <|key -> value, ...|>. Keys are unique and ordered.' }] },
    { cells: [{ type: 'code',    source: 'data = <|"apples" -> 3, "pears" -> 5, "plums" -> 2|>' }] },
    { cells: [{ type: 'code',    source: 'data[["pears"]]' }] },
    { cells: [{ type: 'code',    source: 'data["plums"]' }] },
    { cells: [{ type: 'code',    source: 'nested = <|"r1" -> <|"x" -> 1, "y" -> 2|>|>; nested["r1", "y"]' }] },
    { cells: [{ type: 'code',    source: 'Keys[data]' }] },
    { cells: [{ type: 'code',    source: 'Values[data]' }] },
    { cells: [{ type: 'code',    source: 'First[data]' }] },
    { cells: [{ type: 'code',    source: 'Lookup[data, "figs", 0]' }] },
    { cells: [{ type: 'code',    source: 'KeyFreeQ[data, "figs"]' }] },
    { cells: [{ type: 'section', source: 'Aggregation' }] },
    { cells: [{ type: 'text',    source: 'Counts, GroupBy and Merge are hash-backed — O(n) over large lists.' }] },
    { cells: [{ type: 'code',    source: 'Counts[{1, 2, 2, 3, 3, 3, 1}]' }] },
    { cells: [{ type: 'code',    source: 'GroupBy[Range[10], EvenQ]' }] },
    { cells: [{ type: 'code',    source: 'Merge[{<|"a" -> 1|>, <|"a" -> 2, "b" -> 3|>}, Total]' }] },
    { cells: [{ type: 'code',    source: 'PositionIndex[{a, b, a, c, a, b}]' }] },
    { cells: [{ type: 'code',    source: 'Position[<|"a" -> 1, "b" -> 2, "c" -> 1|>, 1]' }] },
    { cells: [{ type: 'code',    source: 'MapAt[-# &, <|"a" -> 1, "b" -> 9|>, First[Position[<|"a" -> 1, "b" -> 9|>, 9]]]' }] },
    { cells: [{ type: 'section', source: 'Functional threading' }] },
    { cells: [{ type: 'text',    source: 'Map and Select thread over values, keeping keys (Wolfram style).' }] },
    { cells: [{ type: 'code',    source: 'Map[#^2 &, <|"x" -> 3, "y" -> 4|>]' }] },
    { cells: [{ type: 'code',    source: 'Select[<|"a" -> 1, "b" -> 2, "c" -> 3|>, # > 1 &]' }] },
    { cells: [{ type: 'code',    source: 'KeySortBy[<|"bbb" -> 1, "a" -> 2, "cc" -> 3|>, StringLength]' }] },
    { cells: [{ type: 'code',    source: 'SortBy[{{1, 3}, {1, 1}, {2, 0}, {1, 2}}, {First, Last}]' }] },
    { cells: [{ type: 'section', source: 'Aggregation & mutation' }] },
    { cells: [{ type: 'text',    source: 'Total/Min/Max/Mean reduce over values; Part assignment updates in place.' }] },
    { cells: [{ type: 'code',    source: 'Total[<|"a" -> 3, "b" -> 1, "c" -> 2|>]' }] },
    { cells: [{ type: 'code',    source: 'Fold[Plus, 0, <|"a" -> 1, "b" -> 2, "c" -> 3|>]' }] },
    { cells: [{ type: 'code',    source: 'Table[v^2, {v, <|"a" -> 2, "b" -> 3|>}]' }] },
    { cells: [{ type: 'code',    source: 'Sort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]' }] },
    { cells: [{ type: 'code',    source: 'inv = <|"gold" -> 3|>; inv[["silver"]] = 10; inv[["gold"]] = inv[["gold"]] + 1; inv' }] },
    { cells: [{ type: 'code',    source: 'Cases[<|"a" -> 1, "b" -> 2, "c" -> 3|>, x_ /; x > 1]' }] },
    { cells: [{ type: 'code',    source: 'DeleteCases[<|"a" -> 1, "b" -> 2, "c" -> 3|>, x_ /; x > 1]' }] },
    { cells: [{ type: 'code',    source: 'AllTrue[<|"a" -> 2, "b" -> 4, "c" -> 6|>, EvenQ]' }] },
    { cells: [{ type: 'code',    source: 'SortBy[<|"a" -> 3, "b" -> 1, "c" -> 2|>, Identity]' }] },
    { cells: [{ type: 'code',    source: 'MaximalBy[<|"a" -> 1, "b" -> 3, "c" -> 3|>, Identity]' }] },
    { cells: [{ type: 'code',    source: 'TakeLargest[<|"a" -> 3, "b" -> 9, "c" -> 1, "d" -> 6|>, 2]' }] },
    { cells: [{ type: 'code',    source: 'GroupBy[Range[10], EvenQ, Total]' }] },
    { cells: [{ type: 'code',    source: 'ReverseSort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]' }] },
    { cells: [{ type: 'code',    source: 'SelectFirst[<|"a" -> 1, "b" -> 4, "c" -> 6|>, EvenQ]' }] },
    { cells: [{ type: 'code',    source: 'DeleteMissing[Lookup[<|"a" -> 1, "b" -> 2|>, {"a", "z", "b"}]]' }] },
    { cells: [{ type: 'text',    source: 'Pipelines compose: group by key, sum amounts, rank descending.' }] },
    { cells: [{ type: 'code',    source: 'ReverseSort[GroupBy[{{"x", 1}, {"y", 2}, {"x", 3}}, First -> Last, Total]]' }] },
    { cells: [{ type: 'section', source: 'Windowed & statistical' }] },
    { cells: [{ type: 'text',    source: 'MinMax/Median/Variance/Tally reduce over values; Accumulate/Ratios keep the keys aligned (leading key drops for Ratios).' }] },
    { cells: [{ type: 'code',    source: 'MinMax[<|"a" -> 3, "b" -> 1, "c" -> 9|>]' }] },
    { cells: [{ type: 'code',    source: 'Median[<|"a" -> 1, "b" -> 3, "c" -> 5|>]' }] },
    { cells: [{ type: 'code',    source: 'Tally[<|"a" -> 1, "b" -> 1, "c" -> 2|>]' }] },
    { cells: [{ type: 'code',    source: 'Accumulate[<|"jan" -> 10, "feb" -> 20, "mar" -> 5|>]' }] },
    { cells: [{ type: 'code',    source: 'Ratios[<|"q1" -> 10, "q2" -> 20, "q3" -> 15|>]' }] },
    { cells: [{ type: 'code',    source: 'Plus @@ <|"a" -> 1, "b" -> 2, "c" -> 3|>' }] },
    { cells: [{ type: 'section', source: 'Records & alignment' }] },
    { cells: [{ type: 'text',    source: 'Lookup pulls a field across records; KeyUnion aligns key sets; Catenate merges a list of associations.' }] },
    { cells: [{ type: 'code',    source: 'Lookup[{<|"name" -> "Ada", "age" -> 36|>, <|"name" -> "Alan", "age" -> 41|>}, "age"]' }] },
    { cells: [{ type: 'code',    source: 'GroupBy[{<|"team" -> "A", "pts" -> 3|>, <|"team" -> "B", "pts" -> 5|>, <|"team" -> "A", "pts" -> 4|>}, Key["team"]]' }] },
    { cells: [{ type: 'code',    source: 'SortBy[{<|"n" -> "b", "age" -> 41|>, <|"n" -> "a", "age" -> 36|>}, Key["age"]]' }] },
    { cells: [{ type: 'code',    source: 'KeyUnion[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>}]' }] },
    { cells: [{ type: 'code',    source: 'Catenate[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>}]' }] },
    { cells: [{ type: 'code',    source: 'TakeWhile[<|"a" -> 1, "b" -> 2, "c" -> 5, "d" -> 1|>, # < 3 &]' }] },
    { cells: [{ type: 'code',    source: 'MapIndexed[{#2, #1} &, <|"x" -> 3, "y" -> 4|>]' }] },
    { cells: [{ type: 'text',    source: 'The Unicode arrow → (\\[Rule]) parses as ->, so pasted Wolfram code just works.' }] },
    { cells: [{ type: 'code',    source: '<|"gold" → 5, "silver" → 12|>' }] },
    { cells: [{ type: 'section', source: 'Pattern matching' }] },
    { cells: [{ type: 'text',    source: 'Destructure and filter associations with KeyValuePattern.' }] },
    { cells: [{ type: 'code',    source: 'Cases[{<|"t" -> 1|>, <|"t" -> 2|>, <|"x" -> 3|>}, KeyValuePattern[{"t" -> _}]]' }] },
    { cells: [{ type: 'code',    source: 'area[KeyValuePattern[{"w" -> w_, "h" -> h_}]] := w h; area[<|"w" -> 3, "h" -> 4|>]' }] },
  ]);
}

export function addNotebook(title?: string) {
  canvasState.update(s => {
    const n  = s.notebooks.length;
    const x  = 80 + (n % 4) * 680;
    const y  = 60 + Math.floor(n / 4) * 500;
    const nb = makeCard(title ?? '', x, y);
    return { ...s, notebooks: [...s.notebooks, nb] };
  });
}

/** Add a notebook at specific world coordinates — single atomic update. */
export function addNotebookAt(worldX: number, worldY: number, title?: string) {
  const nb = makeCard(title ?? '', worldX, worldY);
  canvasState.update(s => ({ ...s, notebooks: [...s.notebooks, nb] }));
  return nb.id;
}

/* Open a card immediately to the right of `fromId` holding `query`, and mark it
 * to evaluate on mount. Used by the name grid in a `?pat*` result: clicking a
 * symbol should surface its documentation as its own card, not append a cell to
 * the bottom of whatever notebook happened to run the search -- in a long
 * notebook that lands off-screen and reads as nothing having happened. */
/* Open a symbol's reference page as its own card to the right of `fromId`.
 *
 * The page becomes a real notebook rather than one block of rendered text: prose
 * sections are Markdown cells, and every In[n]:= / Out[n]= example in the page
 * becomes an actual code cell, pre-filled with the input and seeded with the
 * output the generator verified against the built binary. So the page reads as
 * complete on arrival, and any example can be edited and re-run in place.
 *
 * The card appears immediately with a placeholder and is filled in when the
 * fetch lands, rather than the click doing nothing until the network settles. */
/* Where the canvas stage sits on screen, published by Canvas.svelte. Needed to
   turn a viewport point (where a symbol was clicked) into a canvas coordinate. */
let stageOrigin = { left: 0, top: 0 };
export function setStageOrigin(left: number, top: number) {
  stageOrigin = { left, top };
}

/** Fill a card with a symbol's reference page: prose becomes Markdown rows, headings become
 *  section rows the notebook can fold, and each recorded example becomes a real code cell with
 *  its recorded result as output.
 *
 *  Extracted from openRefpage so the startup content can build reference cards too, rather than
 *  having a second, subtly different copy of the same page-to-notebook mapping. */
export async function fillRefpage(store: NotebookStore | undefined, name: string) {
  const target: NotebookStore | undefined = store;
  if (!target) return;
  let segments;
  let figures: Record<string, object> = {};
  try {
    const md = await fetchRefpageMarkdown(name);
    segments = splitRefpage(md);
    /* Place the contents after the definition, not above it: the first thing
       a reader wants is what the function does, and an index of a page they
       have not started reading yet is noise at the top. Inserted before the
       heading that follows Description, or at the front if there is no
       Description section. */
    const toc = buildToc(md);
    if (toc) {
      const tocSeg = { kind: 'md' as const, text: toc };
      const descAt = segments.findIndex(
        seg => seg.kind === 'heading' && /^description$/i.test(seg.text));
      let at = 0;
      if (descAt >= 0) {
        const next = segments.findIndex(
          (seg, i) => i > descAt && seg.kind === 'heading');
        at = next >= 0 ? next : segments.length;
      }
      segments = [...segments.slice(0, at), tocSeg, ...segments.slice(at)];
    }
    figures = await fetchRefpageFigures(name);
  } catch (e) {
    const why = e instanceof Error ? e.message : String(e);
    target.setCellSourceAndType(
      `No reference page for \`${name}\`.\n\n` +
      'Pages are generated by `python3 site/generate.py` and mirrored by ' +
      '`npm run sync:refpages`; a symbol added since the last run will not ' +
      `have one yet.\n\n*(${why})*`, 'ref');
    return;
  }
  if (!segments.length) return;

  /* The card already holds one cell, so the first segment rewrites it and the
     rest are appended -- otherwise the page would open under a stray empty
     cell. */
  let first = true;
  for (const seg of segments) {
    if (seg.kind === 'heading') {
      /* 'section' for H2, 'subsection' for H3 -- the notebook folds the rows
         under each, so every part of the page collapses independently. */
      const type = seg.level === 2 ? 'section' : 'subsection';
      if (first) target.setCellSourceAndType(seg.text, type);
      else target.addRow(type, seg.text);
    } else if (seg.kind === 'md') {
      if (first) target.setCellSourceAndType(seg.text, 'ref');
      else target.addRow('ref', seg.text);
    } else {
      let id: string;
      if (first) {
        target.setCellSourceAndType(seg.input, 'code');
        id = get(target)[0].cells[0].id;
      } else {
        id = target.addRow('code', seg.input);
      }
      /* 'expected', not 'expr': this is the recorded result, and an expr with
         no kernel LaTeX gets typeset by KaTeX into something that looks like
         mathematics but is really Mathilda syntax. Running the cell replaces
         it with the live output. */
      const fig = figures[seg.input];
      if (seg.image) {
        /* The recorded picture IS the output of an image-valued example, so it goes where the
           output goes. Running the cell clears it and draws the live image in its place. */
        target.appendOutput(id, {
          kind: 'html',
          html: `<img class="ref-shot" alt="recorded result" src="${seg.image}">`,
        });
      } else if (fig) {
        /* The saved plot, so a graphics example arrives drawn rather than as
           the word -Graphics-. Running the cell redraws it live. */
        target.appendOutput(id, { kind: 'plot', data: fig });
      } else if (seg.output) {
        target.appendOutput(id, { kind: 'expected', text: seg.output });
      }
    }
    first = false;
  }
}

export function openRefpage(fromId: string, name: string,
                            at?: { x: number; y: number } | null) {
  let newId = '';
  /* Assigned inside the synchronous canvasState.update below; TypeScript cannot
     see that a store callback runs immediately, so the type is widened here. */
  let store: NotebookStore | undefined;

  canvasState.update(s => {
    const from = s.notebooks.find(nb => nb.id === fromId);
    /* Open beside the SYMBOL when we know where it was, not beside the whole
       card: on a wide notebook the far edge can be most of a screen away from
       what the reader just pointed at. Falls back to the card's right edge. */
    let x = from ? from.x + from.width + 40 : 0;
    let y = from ? from.y : 0;
    if (at) {
      const wx = (at.x - stageOrigin.left - s.panX) / s.zoom;
      const wy = (at.y - stageOrigin.top - s.panY) / s.zoom;
      x = wx + 90;                 /* clear of the pointer, not under it */
      y = wy - 60;                 /* the symbol sits near the card's top */
    }
    const nb = makeCard(name, x, y);
    nb.refpage = true;
    nb.store.setCellSourceAndType(`Loading the reference page for \`${name}\`...`, 'ref');
    store = nb.store;
    newId = nb.id;
    /* On top: it is what the reader just asked for. */
    const next = { ...s, notebooks: [...s.notebooks, nb], activeId: nb.id };

    return attachAsPane(next, nb.id, fromId, true);
  });

  void fillRefpage(store, name);

  return newId;
}

export function addQueryNotebook(fromId: string, query: string, title?: string) {
  let newId = '';
  canvasState.update(s => {
    const from = s.notebooks.find(nb => nb.id === fromId);
    const x = from ? from.x + from.width + 40 : 0;
    const y = from ? from.y : 0;
    const nb = makeCard(title ?? query, x, y);
    nb.autoRun = true;
    nb.store.addRow('code', query);
    newId = nb.id;
    /* Same reasoning as openRefpage: on the canvas the position is enough, in
       focused mode it would be invisible. Not treated as a docs pane -- a query
       result is a working notebook the reader may well want to keep. */
    return attachAsPane({ ...s, notebooks: [...s.notebooks, nb] }, nb.id, fromId, false);
  });
  return newId;
}

export function removeNotebook(id: string) {
  canvasState.update(s => {
    const remaining = s.notebooks.filter(nb => nb.id !== id);
    /* Closing the last notebook is a no-op -- an empty canvas has nothing to
       act on and no way back. */
    if (remaining.length === 0) return s;
    /* Closing a notebook that is currently on screen has to drop its pane.
       Leaving a focus id pointing at a card that no longer exists rendered an
       empty window with a toolbar whose buttons acted on nothing.
       normalizeFocus handles all of it: dropping the dead id, re-equalizing the
       survivors, promoting a new active pane if this was the active one, and
       falling back to canvas mode when the last pane goes. */
    return normalizeFocus({ ...s, notebooks: remaining });
  });
}

export function setNotebookPos(id: string, x: number, y: number) {
  canvasState.update(s => ({
    ...s,
    notebooks: s.notebooks.map(nb => nb.id === id ? { ...nb, x, y } : nb),
  }));
}

export function setNotebookWidth(id: string, width: number) {
  const clamped = Math.max(320, Math.min(1600, width));
  canvasState.update(s => ({
    ...s,
    notebooks: s.notebooks.map(nb => nb.id === id ? { ...nb, width: clamped } : nb),
  }));
}

export function setNotebookHeight(id: string, height: number | null) {
  const clamped = height === null ? null : Math.max(100, Math.min(4000, height));
  canvasState.update(s => ({
    ...s,
    notebooks: s.notebooks.map(nb => nb.id === id ? { ...nb, height: clamped } : nb),
  }));
}

export function toggleCollapse(id: string) {
  canvasState.update(s => ({
    ...s,
    notebooks: s.notebooks.map(nb =>
      nb.id === id ? { ...nb, collapsed: !nb.collapsed } : nb
    ),
  }));
}

export function renameNotebook(id: string, title: string) {
  canvasState.update(s => ({
    ...s,
    notebooks: s.notebooks.map(nb =>
      nb.id === id ? { ...nb, title } : nb
    ),
  }));
}

export function setPan(panX: number, panY: number) {
  canvasState.update(s => ({ ...s, panX, panY }));
}

/* ---------------------------------------------------------------------------
 * Focus mutation
 *
 * Every change to the focus fields returns through normalizeFocus, so the
 * invariants hold no matter which entry point was used. Without one place to do
 * this, removeNotebook and loadLibraryData each had to remember the whole set.
 * ------------------------------------------------------------------------- */

type FocusState = ReturnType<typeof initialFocus>;

function normalizeFocus<T extends FocusState & { notebooks: CanvasNotebook[] }>(s: T): T {
  /* Drop ids whose notebook is gone, and any duplicate: a notebook rendered in
     two panes would mean two NotebookCards over one store, and `selectedCells`
     and `kernelStatus` are module-global singletons -- so the two panes would
     share cell selection, both register focus callbacks for the same cell ids,
     and both handle one Cmd+click, opening two reference pages. */
  const live = s.focusedIds.filter(id => s.notebooks.some(nb => nb.id === id));
  const focusedIds = [...new Set(live)].slice(0, MAX_PANES);

  /* Re-equalize when the pane count changed; otherwise keep what the user
     dragged. A new pane the user just asked for should be visible at a usable
     width, so equalizing beats splitting one pane's share. */
  const focusedSizes = focusedIds.length === s.focusedSizes.length
    ? s.focusedSizes
    : Array(focusedIds.length).fill(100 / Math.max(1, focusedIds.length));

  /* The toolbar must always point at a real pane. */
  const focusedActiveId = s.focusedActiveId && focusedIds.includes(s.focusedActiveId)
    ? s.focusedActiveId
    : (focusedIds[0] ?? null);

  /* A 2x2 with two panes is just a side-by-side split with dead space. */
  const focusedLayout = s.focusedLayout === 'grid' && focusedIds.length < 3
    ? 'h' as FocusLayout
    : s.focusedLayout;

  return { ...s, focusedIds, focusedSizes, focusedActiveId, focusedLayout };
}

/** Put a just-created card on screen when the window is in focused mode.
 *
 *  On the canvas a new card's x/y is the whole story. In focused mode that
 *  position is invisible, which is exactly why a Cmd+click on a symbol used to
 *  appear to do nothing: the reference page opened perfectly, behind a
 *  full-window notebook.
 *
 *  The new pane goes on the RIGHT and the caret stays where it was. You asked
 *  what a symbol means while writing, so the next Shift+Enter has to still
 *  belong to your own notebook and not to a read-only reference page.
 *
 *  `reuseDocsPane` makes repeated lookups replace the pane already showing
 *  documentation instead of adding one per lookup -- a docs panel is something
 *  you glance at repeatedly, not something you accumulate. */
function attachAsPane<T extends FocusState & { notebooks: CanvasNotebook[] }>(
  s: T, newId: string, fromId: string, reuseDocsPane: boolean,
): T {
  if (s.focusedIds.length === 0) return s;

  const reusable = reuseDocsPane
    ? s.focusedIds.find(id => s.notebooks.find(nb => nb.id === id)?.refpage)
    : undefined;

  let focusedIds: string[];
  if (reusable) {
    focusedIds = s.focusedIds.map(id => (id === reusable ? newId : id));
  } else if (s.focusedIds.length < MAX_PANES) {
    focusedIds = [...s.focusedIds, newId];
  } else {
    /* Full, with nothing to reuse: take the pane the request came from. Its
       notebook stays on the canvas, and that pane is where the reader's
       attention already is. */
    focusedIds = s.focusedIds.map(id => (id === fromId ? newId : id));
  }

  return normalizeFocus({
    ...s,
    focusedIds,
    /* Side by side: "beside what I am reading" is the whole point. Leave a 2x2
       alone, since it is already showing everything at once. */
    focusedLayout: focusedIds.length >= 3 ? s.focusedLayout : 'h',
    focusedActiveId: focusedIds.includes(s.focusedActiveId ?? '')
      ? s.focusedActiveId
      : focusedIds[0],
  });
}

/** Enter focused mode on one notebook, or return to the canvas with null.
 *
 *  Kept with its original signature deliberately. Three callers should not have
 *  to know panes exist: the card's full-screen button, the toolbar's
 *  back-to-canvas control, and the pinch-out gesture. */
export function setFocused(id: string | null) {
  canvasState.update(s => normalizeFocus({
    ...s,
    focusedIds:      id ? [id] : [],
    focusedLayout:   'h',
    focusedActiveId: id,
    focusedSizes:    id ? [100] : [],
  }));
  /* Leaving focused mode is the one place a stale active cell is genuinely
     wrong -- everywhere else the toolbar's validate-on-read handles it. */
  if (id === null) clearActiveCell();
}

/** Add a notebook as another pane. No-op if already shown, or at the cap. */
export function addPaneToFocus(id: string) {
  canvasState.update(s => {
    if (s.focusedIds.includes(id) || s.focusedIds.length >= MAX_PANES) return s;
    return normalizeFocus({ ...s, focusedIds: [...s.focusedIds, id], focusedActiveId: id });
  });
}

/** Open several notebooks tiled, in one step.
 *
 *  This is the route in from the canvas: rubber-band a few cards, then open them
 *  together. Going through addPaneToFocus one at a time would work but would
 *  re-equalise sizes and re-point the active pane on every call. */
export function openPanes(ids: string[], layout: FocusLayout = 'h') {
  canvasState.update(s => {
    const live = ids.filter(id => s.notebooks.some(nb => nb.id === id));
    const kept = [...new Set(live)].slice(0, MAX_PANES);
    if (kept.length === 0) return s;
    return normalizeFocus({
      ...s,
      focusedIds:      kept,
      /* A 2x2 needs three panes; asking for one with two is a side-by-side. */
      focusedLayout:   kept.length >= 3 ? layout : (layout === 'grid' ? 'h' : layout),
      focusedActiveId: kept[0],
      focusedSizes:    Array(kept.length).fill(100 / kept.length),
      activeId:        kept[0],
    });
  });
}

/** Split the window with whichever notebook is next in canvas order and not
 *  already on screen. The one-click path: no picker, and the pane header names
 *  what arrived so a wrong guess is obvious and fixable with + / x. */
export function splitWithNext(layout: FocusLayout = 'h'): string | null {
  let picked: string | null = null;
  canvasState.update(s => {
    if (s.focusedIds.length >= MAX_PANES) return s;
    const next = s.notebooks.find(nb => !s.focusedIds.includes(nb.id));
    if (!next) return s;
    picked = next.id;
    return normalizeFocus({
      ...s,
      focusedIds:      [...s.focusedIds, next.id],
      focusedLayout:   layout,
      /* Keep the caret where it was: the user asked for a second view, not to
         start typing in it. */
      focusedActiveId: s.focusedActiveId,
    });
  });
  return picked;
}

/** Remove a pane, keeping its notebook on the canvas. Emptying returns to canvas.
 *
 *  Distinct from the toolbar's close button, which calls removeNotebook and
 *  deletes the notebook outright. */
export function removePane(id: string) {
  canvasState.update(s => {
    const focusedIds = s.focusedIds.filter(x => x !== id);
    /* Rescale survivors proportionally so the ratios the user dragged survive. */
    const kept  = s.focusedIds.map((x, i) => ({ x, size: s.focusedSizes[i] ?? 0 }))
                              .filter(e => e.x !== id);
    const total = kept.reduce((a, e) => a + e.size, 0);
    const focusedSizes = total > 0
      ? kept.map(e => (e.size * 100) / total)
      : Array(focusedIds.length).fill(100 / Math.max(1, focusedIds.length));
    return normalizeFocus({ ...s, focusedIds, focusedSizes });
  });
  if (get(canvasState).focusedIds.length === 0) clearActiveCell();
}

/** Point the toolbar at a pane. Also raises that card on the canvas, so
 *  returning from focused mode leaves the last-used notebook on top. */
export function setFocusedActive(id: string) {
  canvasState.update(s =>
    s.focusedActiveId === id ? s : normalizeFocus({ ...s, focusedActiveId: id, activeId: id }));
}

export function setFocusLayout(focusedLayout: FocusLayout) {
  canvasState.update(s => normalizeFocus({ ...s, focusedLayout }));
}

export function setFocusedSizes(focusedSizes: number[]) {
  canvasState.update(s => ({ ...s, focusedSizes }));
}

export function setFocusedGrid(focusedGrid: { x: number; y: number }) {
  canvasState.update(s => ({ ...s, focusedGrid }));
}

export function setZoom(zoom: number, cx: number, cy: number) {
  canvasState.update(s => {
    const clamped = Math.max(0.08, Math.min(3, zoom));
    const factor  = clamped / s.zoom;
    return {
      ...s,
      zoom: clamped,
      panX: cx - factor * (cx - s.panX),
      panY: cy - factor * (cy - s.panY),
    };
  });
}

// ---------------------------------------------------------------------------
// Library serialization / deserialization

/** Serialize the entire canvas to a .lb JSON string. Outputs are not saved (ephemeral). */
export function serializeLibrary(title: string): string {
  const state = get(canvasState);
  const notebooks = state.notebooks.map(nb => ({
    id: nb.id,
    title: nb.title,
    x: nb.x,
    y: nb.y,
    width: nb.width,
    collapsed: nb.collapsed,
    rows: nb.store.serialize(),
  }));
  return JSON.stringify({
    version: '1',
    type: 'mathilda-library',
    title,
    notebooks,
  }, null, 2);
}

export type LibraryData = {
  version: string;
  type: string;
  title: string;
  notebooks: Array<{
    id: string;
    title: string;
    x: number;
    y: number;
    width: number;
    collapsed: boolean;
    rows: Array<{ cells: Array<{ type: string; source: string }> }>;
  }>;
};

/** Deserialize a .lb JSON string and replace the current canvas state. */
export function loadLibraryData(json: string): string {
  const data: LibraryData = JSON.parse(json);
  const notebooks: CanvasNotebook[] = data.notebooks.map(nb => {
    const store = createNotebook();
    store.load(nb.rows);
    return {
      id: nb.id,
      title: nb.title,
      x: nb.x,
      y: nb.y,
      width: nb.width ?? 640,
      height: null,
      collapsed: nb.collapsed ?? false,
      store,
    };
  });
  // Reset nextId past the highest id in the file to avoid collisions
  canvasState.set({
    notebooks:   notebooks.length > 0 ? notebooks : [makeCard('Notebook 1', 60, 60)],
    panX:        0,
    panY:        0,
    zoom:        1.0,
    /* Spread rather than listed field by field: this is a full-object set(), so
       any focus field omitted here would silently become undefined and the next
       render would throw on focusedIds.length. */
    ...initialFocus(),
    activeId:    null,
    selectedIds: [],
  });
  /* The stores this load just replaced are gone, so any remembered cell is too. */
  clearActiveCell();
  return data.title ?? 'Untitled Library';
}
