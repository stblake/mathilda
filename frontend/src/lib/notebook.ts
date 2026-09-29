// notebook.ts — 2D notebook: list of rows, each row holds 1+ cells side-by-side

import { writable, get } from 'svelte/store';

/* 'ref' renders a symbol's generated reference page read-only; its `source` is
   the symbol name, not code, and it never goes to the kernel. */
export type CellType = 'code' | 'text' | 'title' | 'subtitle' | 'chapter'
                     | 'section' | 'subsection' | 'subsubsection' | 'ref';

export interface CellStyle {
  id: CellType;
  label: string;
  /** The heading element to render, or null for a style that is not a heading. Doubles as the
   *  isHeading test, so there is no second list to keep in step. */
  tag: 'h1' | 'h2' | 'h3' | null;
  /** Badge glyph in the cell's type picker. */
  icon: string;
  desc: string;
}

/** Every style a reader can put a cell into, in OUTLINE ORDER.
 *
 *  ONE list, because four surfaces present these and each of them fails differently and quietly
 *  when they disagree: the cell's own type picker, the toolbar's Cell Style control, the native
 *  Cell > Convert to items, and -- the one that loses work -- the `.mathilda` file format, whose
 *  Rust `KNOWN_TYPES` rewrites any type it does not recognise to `code` on SAVE as well as on load.
 *  A style added to the UI and not to that list is a cell that silently becomes code the first time
 *  the notebook is saved and reopened. `npm run check:notebook` diffs the two.
 *
 *  `ref` is deliberately absent: a reference-page cell is generated documentation, not a style
 *  anyone chooses, and every control that offers these already refuses to retype one.
 *
 *  ON THE TAGS. They give the DOM a real heading element; they are NOT the visual ladder, which is
 *  carried by CSS keyed on the style id (`.heading-title`, ...). So the tags need not be monotonic
 *  and are not: `section` stays h1 and `subsection` stays h2, exactly as before this list existed,
 *  because every generated reference page is built from those two and re-tagging them would silently
 *  restyle all of them. */
export const CELL_STYLES: CellStyle[] = [
  { id: 'code',          label: 'Code',          tag: null, icon: '▶',   desc: 'Evaluate Mathilda expressions' },
  { id: 'text',          label: 'Text',          tag: null, icon: 'T',   desc: 'Prose / markdown' },
  { id: 'title',         label: 'Title',         tag: 'h1', icon: 'T',   desc: "The notebook's title" },
  { id: 'subtitle',      label: 'Subtitle',      tag: 'h2', icon: 't',   desc: 'Sits under the title' },
  { id: 'chapter',       label: 'Chapter',       tag: 'h2', icon: 'C',   desc: 'Top-level division' },
  { id: 'section',       label: 'Section',       tag: 'h1', icon: '#',   desc: 'Heading' },
  { id: 'subsection',    label: 'Subsection',    tag: 'h2', icon: '##',  desc: 'Inside a section' },
  { id: 'subsubsection', label: 'Subsubsection', tag: 'h3', icon: '###', desc: 'Inside a subsection' },
];

/** True for the styles rendered as a one-line editable heading rather than as code or prose. */
export function isHeading(type: CellType): boolean {
  return CELL_STYLES.some(s => s.id === type && s.tag !== null);
}

/** The heading element for a style, or null when it is not a heading. */
export function headingTag(type: CellType): 'h1' | 'h2' | 'h3' | null {
  return CELL_STYLES.find(s => s.id === type)?.tag ?? null;
}
export type CellStatus = 'idle' | 'running' | 'done' | 'error';

export type OutputItem =
  | { kind: 'expr';   text: string; latex?: string }  // latex = StandardForm LaTeX from kernel
  | { kind: 'usage';  text: string; symbol?: string }  // ?sym help text: preformatted, never math
  | { kind: 'expected'; text: string }                // reference-page example: the
                                                      // recorded, verified result,
                                                      // shown until the cell is run
  | { kind: 'names';  names: string[] }               // ?pat* symbol search: laid out as a grid
  | { kind: 'error';  text: string }
  | { kind: 'stream'; text: string }                 // Print output
  | { kind: 'message'; text: string }                // kernel warning, Head::tag: ...
  | { kind: 'plot';   data: object }
  /* A raster result. `data` is base64 RGBA, w*h*4 bytes, ready for putImageData -- so the
     browser does no per-pixel work. A volume sends ONE slice (the middle) and carries `depth`
     and `slice` so a scrubber can ask for others later. */
  /** One face of a volume: its own pixel size plus base64 RGBA. */
  | { kind: 'image';  w: number; h: number; channels: number; data: string;
      faces?: Record<string, { w: number; h: number; data: string }>;
                      depth?: number; slice?: number }
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

export function createNotebook() {
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

    /** Retype the notebook's FIRST cell in place and give it a source.
     *
     * createNotebook() seeds one empty code row, so a card built for a single
     * purpose (a reference page) would otherwise carry a stray empty cell above
     * its content. Rewrites that row instead of appending after it. */
    setCellSourceAndType(source: string, type: CellType) {
      update(rows => {
        if (rows.length === 0 || rows[0].cells.length === 0) return rows;
        const first = rows[0];
        return [
          { ...first, cells: [{ ...first.cells[0], source, type }, ...first.cells.slice(1)] },
          ...rows.slice(1),
        ];
      });
    },

    /** Change a cell's type, keeping its output.
     *
     *  This used to clear `output` and `execIdx`, which made retyping a
     *  destructive act: a code cell with a result on screen lost it, silently,
     *  and switching back gave you an empty cell. That was tolerable while the
     *  only way to retype was a 12px badge in the gutter; it is not now that the
     *  toolbar's cell-style control puts it one click away.
     *
     *  Keeping the output costs nothing. Only code cells render an output area,
     *  so a retained result is invisible on a text or heading cell, and
     *  serialize() persists only {type, source} so nothing extra reaches disk.
     *  Switching code -> text -> code now restores the result instead of
     *  discarding it. */
    setCellType(id: string, type: CellType) {
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => c.id === id ? { ...c, type } : c),
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

    /** Replace a cell's In[n] with the line the KERNEL used.
     *
     *  stampExec's counter is a placeholder that appears the instant you press
     *  Shift+Enter; this is the real number. It has to come from the kernel
     *  because `%3` / `In[3]` / `Out[3]` resolve against the kernel's $Line, and
     *  one kernel serves every notebook on the canvas -- a local counter would
     *  label a cell In[2] while `%2` addressed a line typed in another pane. */
    setExec(id: string, n: number) {
      update(rows => rows.map(row => ({
        ...row,
        cells: row.cells.map(c => c.id === id ? { ...c, execIdx: n } : c),
      })));
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
