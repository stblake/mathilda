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

/** Leave a fresh input cell below the last row, so evaluating the bottom cell of a notebook does not
 *  leave the caret with nowhere to go. Returns the new cell's id, or null when nothing should be
 *  added -- which is the caller's cue to leave the caret where it is.
 *
 *  ROWS, not cells: a row can hold several cells side by side, and a sibling beside the evaluated
 *  cell is not BELOW it, so the test is "is this the last row" rather than "is this the last cell".
 *
 *  Refused for anything but a non-empty code cell, which is deliberately the SAME condition runCell
 *  uses to decide whether to evaluate at all. Otherwise Shift+Enter on the empty cell this just
 *  created would append another, and holding the key would fill the notebook. */
export function appendCellIfLast(store: Store, cellId: string): string | null {
  const found = store.findCell(cellId);
  if (!found) return null;
  const cell = found.row.cells[found.cellIdx];
  if (!cell || cell.type !== 'code' || !cell.source.trim()) return null;
  if (found.rowIdx !== store.getRows().length - 1) return null;
  return store.addRow('code');
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

/** The source of the nearest code cell ABOVE `cellId` that has any, or null when there is none.
 *
 *  Document order across rows, so a side-by-side row is searched right-to-left before moving up, and
 *  empty cells are skipped -- "the input above" means the last thing you actually typed, and the
 *  freshly appended empty cell that evaluating the last cell leaves behind would otherwise always be
 *  the answer. Split out from copyInputFromAbove so the lookup is testable without an editor. */
export function previousInputSource(store: Store, cellId: string): string | null {
  const cells = store.allCells();
  const i = cells.findIndex(c => c.id === cellId);
  if (i <= 0) return null;
  for (let k = i - 1; k >= 0; k--) {
    if (cells[k].type === 'code' && cells[k].source.trim()) return cells[k].source;
  }
  return null;
}

/** Mathematica's Edit > Copy Input from Above (Cmd+L): put the previous input at the caret.
 *
 *  Inserted at the insertion point, replacing the selection, rather than replacing the cell: on the
 *  empty cell this is nearly always used from the two are identical, and on a cell you have started
 *  typing in, replacing the lot would destroy work. Returns false when there is no input above, so
 *  the caller can leave the keystroke alone instead of clearing the cell. */
export function copyInputFromAbove(store: Store, cellId: string, view: EditorView): boolean {
  const src = previousInputSource(store, cellId);
  if (src === null) return false;
  const sel = view.state.selection.main;
  view.dispatch({
    changes: { from: sel.from, to: sel.to, insert: src },
    selection: { anchor: sel.from + src.length },
  });
  view.focus();
  return true;
}

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
