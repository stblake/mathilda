// notebook.ts — 2D notebook: list of rows, each row holds 1+ cells side-by-side

import { writable, get } from 'svelte/store';

export type CellType = 'code' | 'text' | 'section' | 'subsection';
export type CellStatus = 'idle' | 'running' | 'done' | 'error';

export type OutputItem =
  | { kind: 'expr';   text: string }
  | { kind: 'error';  text: string }
  | { kind: 'stream'; text: string }
  | { kind: 'plot';   data: object }
  | { kind: 'html';   html: string };

export type Cell = {
  id: string;
  type: CellType;
  source: string;
  status: CellStatus;
  output: OutputItem[];
  execIdx?: number;
};

export type NotebookRow = {
  id: string;
  cells: Cell[];
};

let _nextCellId = 1;
let _nextRowId  = 1;
let _execCounter = 0;

function newCellId(): string { return `c${_nextCellId++}`; }
function newRowId():  string { return `r${_nextRowId++}`; }

function makeCell(type: CellType = 'code', source = ''): Cell {
  return { id: newCellId(), type, source, status: 'idle', output: [] };
}

function makeRow(type: CellType = 'code', source = ''): NotebookRow {
  return { id: newRowId(), cells: [makeCell(type, source)] };
}

// ---------------------------------------------------------------------------
// Selection

export const selectedCells = writable<Set<string>>(new Set());
export let lastSelectedId: string | null = null;

export function selectOnly(id: string) {
  lastSelectedId = id;
  selectedCells.set(new Set([id]));
}
export function toggleSelect(id: string) {
  selectedCells.update(s => {
    const n = new Set(s);
    n.has(id) ? n.delete(id) : n.add(id);
    lastSelectedId = id;
    return n;
  });
}
export function clearSelection() {
  lastSelectedId = null;
  selectedCells.set(new Set());
}

export function rangeSelect(toId: string) {
  const cells = notebook.allCells();
  const fromId = lastSelectedId;
  if (!fromId) { selectOnly(toId); return; }
  const a = cells.findIndex(c => c.id === fromId);
  const b = cells.findIndex(c => c.id === toId);
  if (a < 0 || b < 0) { selectOnly(toId); return; }
  const lo = Math.min(a, b), hi = Math.max(a, b);
  lastSelectedId = toId;
  selectedCells.set(new Set(cells.slice(lo, hi + 1).map(c => c.id)));
}

// ---------------------------------------------------------------------------
// Notebook store (holds NotebookRow[])

function createNotebook() {
  const { subscribe, update, set } = writable<NotebookRow[]>([makeRow()]);

  // --- helpers ---

  function findCell(cells: NotebookRow[], cellId: string): { row: NotebookRow; rowIdx: number; cellIdx: number } | null {
    for (let ri = 0; ri < cells.length; ri++) {
      const ci = cells[ri].cells.findIndex(c => c.id === cellId);
      if (ci >= 0) return { row: cells[ri], rowIdx: ri, cellIdx: ci };
    }
    return null;
  }

  return {
    subscribe,

    // --- row operations ---

    /** Insert a new row at absolute row index. Returns new cell id. */
    insertRowAt(rowIdx: number, type: CellType = 'code', source = ''): string {
      const row = makeRow(type, source);
      update(rows => [...rows.slice(0, rowIdx), row, ...rows.slice(rowIdx)]);
      return row.cells[0].id;
    },

    /** Append a row at the end. Returns new cell id. */
    addRow(type: CellType = 'code', source = ''): string {
      const row = makeRow(type, source);
      update(rows => [...rows, row]);
      return row.cells[0].id;
    },

    // --- cell-within-row operations ---

    /** Insert a cell at position cellIdx inside the row identified by rowId. Returns new cell id. */
    insertCellInRow(rowId: string, cellIdx: number, type: CellType = 'code', source = ''): string {
      const cell = makeCell(type, source);
      update(rows => rows.map(row => {
        if (row.id !== rowId) return row;
        const cells = [...row.cells.slice(0, cellIdx), cell, ...row.cells.slice(cellIdx)];
        return { ...row, cells };
      }));
      return cell.id;
    },

    // --- removal ---

    removeCell(cellId: string) {
      update(rows => {
        const next = rows.map(row => ({
          ...row,
          cells: row.cells.filter(c => c.id !== cellId),
        })).filter(row => row.cells.length > 0);
        return next.length === 0 ? [makeRow()] : next;
      });
      selectedCells.update(s => { s.delete(cellId); return new Set(s); });
    },

    removeCells(ids: Set<string>) {
      update(rows => {
        const next = rows.map(row => ({
          ...row,
          cells: row.cells.filter(c => !ids.has(c.id)),
        })).filter(row => row.cells.length > 0);
        return next.length === 0 ? [makeRow()] : next;
      });
      selectedCells.update(s => { ids.forEach(id => s.delete(id)); return new Set(s); });
    },

    // --- mutation ---

    updateSource(id: string, source: string) {
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => c.id === id ? { ...c, source } : c),
      })));
    },

    setCellType(id: string, type: CellType) {
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => c.id === id ? { ...c, type, output: [], execIdx: undefined } : c),
      })));
    },

    setStatus(id: string, status: CellStatus) {
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => c.id === id ? { ...c, status } : c),
      })));
    },

    clearOutput(id: string) {
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => c.id === id ? { ...c, output: [] } : c),
      })));
    },

    appendOutput(id: string, item: OutputItem) {
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => c.id === id ? { ...c, output: [...c.output, item] } : c),
      })));
    },

    appendStream(id: string, text: string) {
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => {
          if (c.id !== id) return c;
          const out = [...c.output];
          if (out.length > 0 && out[out.length - 1].kind === 'stream') {
            out[out.length - 1] = { kind: 'stream', text: (out[out.length - 1] as any).text + text };
          } else {
            out.push({ kind: 'stream', text });
          }
          return { ...c, output: out };
        }),
      })));
    },

    stampExec(id: string): number {
      _execCounter++;
      const n = _execCounter;
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => c.id === id ? { ...c, execIdx: n } : c),
      })));
      return n;
    },

    resetExecCounter() {
      _execCounter = 0;
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => ({ ...c, execIdx: undefined })),
      })));
    },

    // --- serialization ---

    serialize() {
      return get({ subscribe }).map(row => ({
        cells: row.cells.map(c => ({ type: c.type, source: c.source })),
      }));
    },

    load(data: Array<{ cells: Array<{ type: string; source: string }> }>) {
      _execCounter = 0;
      const rows: NotebookRow[] = data.map(rowData => ({
        id: newRowId(),
        cells: rowData.cells.map(cd => makeCell(cd.type as CellType, cd.source)),
      }));
      set(rows.length > 0 ? rows : [makeRow()]);
      selectedCells.set(new Set());
    },

    // Legacy single-cell serialization (for .mathilda files without row structure)
    serializeLegacy() {
      return get({ subscribe }).flatMap(row =>
        row.cells.map(c => ({ type: c.type, source: c.source }))
      );
    },

    loadLegacy(cells: Array<{ type: string; source: string }>) {
      _execCounter = 0;
      set(cells.map(cd => ({ id: newRowId(), cells: [makeCell(cd.type as CellType, cd.source)] })));
      selectedCells.set(new Set());
    },

    allCells(): Cell[] {
      return get({ subscribe }).flatMap(row => row.cells);
    },

    findCell(cellId: string) {
      return findCell(get({ subscribe }), cellId);
    },

    getRows(): NotebookRow[] { return get({ subscribe }); },
  };
}

export const notebook = createNotebook();

// ---------------------------------------------------------------------------

export type KernelStatusValue = 'starting' | 'ready' | 'busy' | 'restarting' | 'dead';
export const kernelStatus = writable<KernelStatusValue>('starting');
