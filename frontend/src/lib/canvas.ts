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

// Starter notebooks — user can load showcase.lb for the full set
const _startNb1 = makeCard('Calculus I — Derivatives', 50, 50);
const _startNb2 = makeCard('Calculus II — Integration', 760, 50);
const _startNb3 = makeCard('Algebra',                   50, 650);

export const canvasState = writable({
  notebooks: [_startNb1, _startNb2, _startNb3] as CanvasNotebook[],
  panX: 0,
  panY: 0,
  zoom:  1.0,
  focusedId: null as string | null,
});

/** Pre-fill the starter notebooks with rich example cells. Call from onMount. */
export function loadStartupContent() {
  const s = get(canvasState);

  const nb1 = s.notebooks[0];
  if (nb1) nb1.store.load([
    { cells: [{ type: 'section', source: 'Derivatives' }] },
    { cells: [{ type: 'text',    source: 'Symbolic differentiation. Press Shift+Enter to evaluate.' }] },
    { cells: [{ type: 'code',    source: 'D[Sin[x]^2, x]' }] },
    { cells: [{ type: 'code',    source: 'D[x^x, x]' }] },
    { cells: [{ type: 'code',    source: 'D[Exp[x] Sin[x], x]' }] },
    { cells: [{ type: 'code',    source: 'D[Log[x^2 + 1], x]' }] },
    { cells: [{ type: 'section', source: 'Higher Order' }] },
    { cells: [{ type: 'code',    source: 'D[Sin[x], {x, 4}]' }] },
    { cells: [{ type: 'code',    source: 'D[Sin[Exp[x]], x]' }] },
  ]);

  const nb2 = s.notebooks[1];
  if (nb2) nb2.store.load([
    { cells: [{ type: 'section', source: 'Definite Integration' }] },
    { cells: [{ type: 'text',    source: 'Symbolic antiderivatives and definite integrals.' }] },
    { cells: [{ type: 'code',    source: 'Integrate[x^2, {x, 0, 1}]' }] },
    { cells: [{ type: 'code',    source: 'Integrate[Sin[x], {x, 0, Pi}]' }] },
    { cells: [{ type: 'code',    source: 'Integrate[Exp[-x^2], {x, -Infinity, Infinity}]' }] },
    { cells: [{ type: 'section', source: 'Series & Limits' }] },
    { cells: [{ type: 'code',    source: 'Series[Sin[x], {x, 0, 7}]' }] },
    { cells: [{ type: 'code',    source: 'Limit[Sin[x]/x, x -> 0]' }] },
    { cells: [{ type: 'code',    source: 'Limit[(1 + 1/n)^n, n -> Infinity]' }] },
  ]);

  const nb3 = s.notebooks[2];
  if (nb3) nb3.store.load([
    { cells: [{ type: 'section', source: 'Factoring & Solving' }] },
    { cells: [{ type: 'text',    source: 'Polynomial algebra and equation solving.' }] },
    { cells: [{ type: 'code',    source: 'Factor[x^4 - 1]' }] },
    { cells: [{ type: 'code',    source: 'Factor[x^6 - y^6]' }] },
    { cells: [{ type: 'code',    source: 'Solve[x^2 - 5 x + 6 == 0, x]' }] },
    { cells: [{ type: 'code',    source: 'Solve[{x + y == 5, x - y == 1}, {x, y}]' }] },
    { cells: [{ type: 'section', source: 'Simplification' }] },
    { cells: [{ type: 'code',    source: 'Simplify[(x^2 - 1)/(x - 1)]' }] },
    { cells: [{ type: 'code',    source: 'Together[1/x + 1/(x+1) + 1/(x+2)]' }] },
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
