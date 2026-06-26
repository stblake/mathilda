// notebook.ts — reactive notebook state

import { writable, get } from 'svelte/store';

export type CellStatus = 'idle' | 'running' | 'done' | 'error';

export type OutputItem =
  | { kind: 'expr';   text: string }
  | { kind: 'error';  text: string }
  | { kind: 'stream'; text: string }
  | { kind: 'plot';   data: object }
  | { kind: 'html';   html: string };

export type Cell = {
  id: string;
  type: 'code' | 'prose';
  source: string;
  status: CellStatus;
  output: OutputItem[];
  execIdx?: number;   // In[n] index, set when cell runs
};

let _nextCellId = 1;
let _execCounter = 0;

function newId(): string { return `cell-${_nextCellId++}`; }

function makeCell(type: Cell['type'] = 'code', source = ''): Cell {
  return { id: newId(), type, source, status: 'idle', output: [] };
}

// ---------------------------------------------------------------------------
// Selection store — separate so CodeCell can subscribe independently

export const selectedCells = writable<Set<string>>(new Set());
export let lastSelectedId: string | null = null;

// ---------------------------------------------------------------------------
// Notebook store

function createNotebook() {
  const { subscribe, update, set } = writable<Cell[]>([makeCell()]);

  return {
    subscribe,

    addCell(afterId?: string, type: Cell['type'] = 'code') {
      update(cells => {
        const cell = makeCell(type);
        if (!afterId) return [...cells, cell];
        const idx = cells.findIndex(c => c.id === afterId);
        const pos = idx >= 0 ? idx + 1 : cells.length;
        return [...cells.slice(0, pos), cell, ...cells.slice(pos)];
      });
    },

    removeCell(id: string) {
      update(cells => cells.length <= 1 ? cells : cells.filter(c => c.id !== id));
      selectedCells.update(s => { s.delete(id); return new Set(s); });
    },

    removeCells(ids: Set<string>) {
      update(cells => {
        const result = cells.filter(c => !ids.has(c.id));
        return result.length === 0 ? [makeCell()] : result;
      });
      selectedCells.update(s => { ids.forEach(id => s.delete(id)); return new Set(s); });
    },

    updateSource(id: string, source: string) {
      update(cells => cells.map(c => c.id === id ? { ...c, source } : c));
    },

    setStatus(id: string, status: CellStatus) {
      update(cells => cells.map(c => c.id === id ? { ...c, status } : c));
    },

    clearOutput(id: string) {
      update(cells => cells.map(c => c.id === id ? { ...c, output: [] } : c));
    },

    appendOutput(id: string, item: OutputItem) {
      update(cells => cells.map(c =>
        c.id === id ? { ...c, output: [...c.output, item] } : c
      ));
    },

    appendStream(id: string, text: string) {
      update(cells => cells.map(c => {
        if (c.id !== id) return c;
        const out = [...c.output];
        if (out.length > 0 && out[out.length - 1].kind === 'stream') {
          out[out.length - 1] = { kind: 'stream', text: (out[out.length - 1] as any).text + text };
        } else {
          out.push({ kind: 'stream', text });
        }
        return { ...c, output: out };
      }));
    },

    /** Stamp execIdx and return the index (In[n] number). */
    stampExec(id: string): number {
      _execCounter++;
      const n = _execCounter;
      update(cells => cells.map(c => c.id === id ? { ...c, execIdx: n } : c));
      return n;
    },

    resetExecCounter() {
      _execCounter = 0;
      update(cells => cells.map(c => ({ ...c, execIdx: undefined })));
    },

    serialize(): Array<{ type: string; source: string }> {
      return get({ subscribe }).map(c => ({ type: c.type, source: c.source }));
    },

    load(cells: Array<{ type: string; source: string }>) {
      _execCounter = 0;
      set(cells.map(c => makeCell(c.type as Cell['type'], c.source)));
      selectedCells.set(new Set());
    },

    reset() {
      _execCounter = 0;
      set([makeCell()]);
      selectedCells.set(new Set());
    },

    getCells(): Cell[] { return get({ subscribe }); },
  };
}

export const notebook = createNotebook();

// ---------------------------------------------------------------------------
// Selection helpers

export function selectOnly(id: string) {
  lastSelectedId = id;
  selectedCells.set(new Set([id]));
}

export function toggleSelect(id: string) {
  selectedCells.update(s => {
    const n = new Set(s);
    if (n.has(id)) n.delete(id); else n.add(id);
    lastSelectedId = id;
    return n;
  });
}

export function rangeSelect(toId: string) {
  const cells = notebook.getCells();
  const fromId = lastSelectedId;
  if (!fromId) { selectOnly(toId); return; }
  const a = cells.findIndex(c => c.id === fromId);
  const b = cells.findIndex(c => c.id === toId);
  if (a < 0 || b < 0) { selectOnly(toId); return; }
  const lo = Math.min(a, b), hi = Math.max(a, b);
  lastSelectedId = toId;
  selectedCells.set(new Set(cells.slice(lo, hi + 1).map(c => c.id)));
}

export function clearSelection() {
  lastSelectedId = null;
  selectedCells.set(new Set());
}

// ---------------------------------------------------------------------------
// Kernel status store

export type KernelStatusValue = 'starting' | 'ready' | 'busy' | 'restarting' | 'dead';
export const kernelStatus = writable<KernelStatusValue>('starting');
