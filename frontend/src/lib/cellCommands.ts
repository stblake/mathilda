// cellCommands.ts — structural edits on cells, and text edits inside one.
//
// Plain functions of (store, cellId) or (view), not methods on the pane's action
// table: they need nothing from the card's private closures, so keeping them here
// makes them readable in isolation and testable without mounting a component.
//
// Every one of these is destructive in some degree, so each states what it does
// when it cannot do the obvious thing rather than throwing or silently doing
// half the job.

import { indentMore, indentLess, toggleComment } from '@codemirror/commands';
import type { EditorView } from '@codemirror/view';
import type { createNotebook, CellType } from './notebook';

type Store = ReturnType<typeof createNotebook>;

/** Split a cell in two at `at`, keeping the head and moving the tail into a new
 *  cell of the same type directly below. Returns the new cell's id.
 *
 *  `at` of null (or a caret that is not live) splits at the end, which produces
 *  an empty cell below -- the same thing pressing Enter at the end of a cell
 *  would give you, and better than refusing. */
export function splitCell(store: Store, cellId: string, at: number | null): string | null {
  const found = store.findCell(cellId);
  if (!found) return null;
  const cell = found.row.cells[found.cellIdx];
  if (!cell) return null;

  const pos  = at == null ? cell.source.length : Math.max(0, Math.min(at, cell.source.length));
  const head = cell.source.slice(0, pos);
  const tail = cell.source.slice(pos);

  store.updateSource(cellId, head);
  return store.insertRowAt(found.rowIdx + 1, cell.type, tail);
}

/** Join a cell with the first cell of the row below it, newline-separated.
 *
 *  Returns false when there is no row below, which is the caller's cue to keep
 *  the control disabled rather than appear to work. */
export function mergeCellDown(store: Store, cellId: string): boolean {
  const found = store.findCell(cellId);
  if (!found) return false;
  const rows = store.getRows();
  const below = rows[found.rowIdx + 1];
  const target = below?.cells[0];
  if (!target) return false;

  const cell = found.row.cells[found.cellIdx];
  if (!cell) return false;

  /* Joined with a newline rather than concatenated: two code cells run as two
     statements, and gluing them into one line would change what they mean. */
  store.updateSource(cellId, `${cell.source}\n${target.source}`);
  store.removeCell(target.id);
  return true;
}

/** Copy a cell into a new row directly below. Outputs are deliberately NOT
 *  copied: a duplicated cell has not been evaluated, and showing the original's
 *  result under it would be a lie about what the kernel has computed. */
export function duplicateCell(store: Store, cellId: string): string | null {
  const found = store.findCell(cellId);
  if (!found) return null;
  const cell = found.row.cells[found.cellIdx];
  if (!cell) return null;
  return store.insertRowAt(found.rowIdx + 1, cell.type, cell.source);
}

/** Delete a cell. The store reseeds an empty notebook rather than leaving one
 *  with nothing in it, so deleting the last cell is safe. */
export function deleteCell(store: Store, cellId: string) {
  store.removeCell(cellId);
}

/** Convert a cell to another style, keeping its source and its output.
 *
 *  The ONE implementation behind every "convert" control -- the toolbar's
 *  cell-style combo, the native Cell > Convert to ... items, and anything later.
 *  Two routes used to exist and only one of them touched the store: the menu
 *  called retypeActiveCell, which updates the toolbar's record of the active
 *  cell and nothing else, so Convert to Text changed the combo's label while the
 *  cell stayed code.
 *
 *  Returns true when the cell's type actually changed. A reference-page cell is
 *  generated documentation, not the reader's to retype, so it is refused (as the
 *  toolbar's locked combo already did); an unknown id or a no-op conversion also
 *  returns false. store.setCellType keeps the output, so code -> text -> code
 *  restores the result rather than discarding it. */
export function convertCell(store: Store, cellId: string, type: CellType): boolean {
  const found = store.findCell(cellId);
  const cell = found?.row.cells[found.cellIdx];
  if (!cell || cell.type === type || cell.type === 'ref' || type === 'ref') return false;
  store.setCellType(cellId, type);
  return true;
}

// ---------------------------------------------------------------------------
// Text edits inside a code cell. These act on a live EditorView, so they need
// the caret to still be in the editor -- which is why every toolbar button
// suppresses pointerdown's default and never lets the editor blur.

export function indentCode(view: EditorView) { indentMore(view); view.focus(); }
export function outdentCode(view: EditorView) { indentLess(view); view.focus(); }

/** Comment or uncomment the selection.
 *
 *  Works only because CellShell declares Mathilda's (* ... *) block comment as
 *  language data: every comment command in @codemirror/commands reads
 *  commentTokens from language data and silently does nothing without it. */
export function commentCode(view: EditorView) { toggleComment(view); view.focus(); }

/** Duplicate the line the caret is on, below itself. */
export function duplicateLine(view: EditorView) {
  const { state } = view;
  const line = state.doc.lineAt(state.selection.main.head);
  view.dispatch({
    changes: { from: line.to, insert: `\n${line.text}` },
    /* Caret to the same column of the copy, so a repeated press stacks copies
       rather than editing the original. */
    selection: { anchor: line.to + 1 + (state.selection.main.head - line.from) },
  });
  view.focus();
}
