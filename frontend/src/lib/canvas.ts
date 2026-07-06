// canvas.ts — infinite canvas state: pan, zoom, and named notebooks

import { writable, get } from 'svelte/store';

// Re-export the notebook factory so each card gets its own store instance.
export { createNotebook } from './notebook';
import { createNotebook } from './notebook';

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

export const canvasState = writable({
  notebooks: [_nb1,_nb2,_nb3,_nb4,_nb5,_nb6,_nb7,_nb8,_nb9] as CanvasNotebook[],
  panX: 0,
  panY: 0,
  zoom:  1.0,
  focusedId: null as string | null,
});

/** Pre-fill the starter notebooks with rich example cells. Call from onMount. */
export function loadStartupContent() {
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

export function removeNotebook(id: string) {
  canvasState.update(s => {
    const remaining = s.notebooks.filter(nb => nb.id !== id);
    return { ...s, notebooks: remaining.length > 0 ? remaining : s.notebooks };
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

export function setFocused(id: string | null) {
  canvasState.update(s => ({ ...s, focusedId: id }));
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
    notebooks: notebooks.length > 0 ? notebooks : [makeCard('Notebook 1', 60, 60)],
    panX: 0,
    panY: 0,
    zoom: 1.0,
  });
  return data.title ?? 'Untitled Library';
}
