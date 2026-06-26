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
  collapsed: boolean;
  store:     NotebookStore;
};

let _nextId = 1;
const PALETTE = ['#4a90e2','#e67e22','#27ae60','#8e44ad','#e74c3c','#16a085'];

function makeCard(title: string, x: number, y: number): CanvasNotebook {
  const id = `nb-${_nextId}`;
  return {
    id,
    title: title || `Notebook ${_nextId++}`,
    x, y,
    width: 640,
    collapsed: false,
    store: createNotebook(),
  };
}

// ---------------------------------------------------------------------------

export const canvasState = writable({
  notebooks: [makeCard('Notebook 1', 60, 60)] as CanvasNotebook[],
  panX: 0,
  panY: 0,
  zoom:  1.0,
});

export function addNotebook(title?: string) {
  canvasState.update(s => {
    // Offset new cards so they don't stack exactly.
    const n   = s.notebooks.length;
    const x   = 80  + (n % 4) * 680;
    const y   = 60  + Math.floor(n / 4) * 500;
    const nb  = makeCard(title ?? '', x, y);
    return { ...s, notebooks: [...s.notebooks, nb] };
  });
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
