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
};

let _nextId = 1;
function newId(): string { return `cell-${_nextId++}`; }

function makeCell(type: Cell['type'] = 'code', source = ''): Cell {
  return { id: newId(), type, source, status: 'idle', output: [] };
}

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
      update(cells => {
        if (cells.length <= 1) return cells; // keep at least one cell
        return cells.filter(c => c.id !== id);
      });
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

    /** Append to the last stream item or create a new one. */
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

    serialize(): Array<{ type: string; source: string }> {
      return get({ subscribe }).map(c => ({ type: c.type, source: c.source }));
    },

    load(cells: Array<{ type: string; source: string }>) {
      set(cells.map(c => makeCell(c.type as Cell['type'], c.source)));
    },

    reset() {
      set([makeCell()]);
    },
  };
}

export const notebook = createNotebook();

// ---------------------------------------------------------------------------
// Kernel status store

export type KernelStatusValue = 'starting' | 'ready' | 'busy' | 'restarting' | 'dead';
export const kernelStatus = writable<KernelStatusValue>('starting');
