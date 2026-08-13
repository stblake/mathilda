// active.ts — which cell the caret is in, and how to reach its editor.
//
// Nothing tracked this before. `selectedCells` in notebook.ts is a different
// notion: it is populated only by gutter/header clicks, and clicking a cell BODY
// calls clearSelection(), so entering edit mode actively EMPTIES it. A toolbar
// with a cell-style combo and a context-sensitive Text/Code group needs the
// opposite -- the single cell the caret is in.
//
// This is a leaf module on purpose. It could have gone in notebook.ts, but:
//   * createNotebook() has no notebook id, and multi-pane needs one; and
//   * notebook.ts's selectedCells is a global singleton shared by all nine store
//     instances -- a latent multi-pane bug. Don't sit next to it.
// It imports only a TYPE from notebook.ts, so it can be imported from anywhere
// (CellShell, NotebookCard, Toolbar, canvas.ts) without a cycle.

import { writable, get } from 'svelte/store';
import type { EditorView } from '@codemirror/view';
import type { CellType } from './notebook';

/** Plain data only. Every reactive block in the toolbar reads this, so it must
 *  be cheap to compare and must not churn. */
export type ActiveCell = {
  /** Which notebook owns the cell. Required, not decorative: in split mode two
   *  notebooks are on screen at once, so the toolbar has to know that the cell
   *  it remembers belongs to the pane it is currently driving. */
  notebookId: string;
  cellId: string;
  cellType: CellType;
  /** Does this cell's editor hold DOM focus RIGHT NOW?
   *
   *  Used for EXACTLY ONE thing: commands that need a live caret offset, such as
   *  split-at-caret, which fall back to end-of-document when it is false. Never
   *  for enabling or disabling a control -- see markBlurred(). */
  focused: boolean;
};

export const activeCell = writable<ActiveCell | null>(null);

/** An imperative handle onto one cell's editor.
 *
 *  Deliberately NOT in a store: it holds a live EditorView, and a store would
 *  re-run every subscriber each time any cell mounts or unmounts. Commands fetch
 *  the handle on demand instead. */
export type CellHandle = {
  /** Code cells only. */
  view: EditorView | null;
  /** The contenteditable of a text/section/subsection cell. */
  el: HTMLElement | null;
  focus: () => void;
};

const handles = new Map<string, CellHandle>();

export function registerHandle(cellId: string, h: CellHandle) { handles.set(cellId, h); }
export function unregisterHandle(cellId: string) { handles.delete(cellId); }
export function getHandle(cellId: string): CellHandle | null { return handles.get(cellId) ?? null; }

/** The handle for whatever the toolbar currently considers active. */
export function activeHandle(): CellHandle | null {
  const a = get(activeCell);
  return a ? getHandle(a.cellId) : null;
}

export function setActiveCell(notebookId: string, cellId: string, cellType: CellType) {
  activeCell.set({ notebookId, cellId, cellType, focused: true });
}

/** Blur: keep the record, flip the flag. NEVER clears.
 *
 *  This is the crux of making a toolbar work at all. Clicking a toolbar button
 *  blurs the editor before the click handler runs, so a store that emptied on
 *  blur would leave every button acting on nothing. Two mechanisms guard this,
 *  and both are needed:
 *
 *   1. Toolbar buttons use on:pointerdown|preventDefault, which suppresses the
 *      focus transfer entirely -- so the editor never blurs and the live text
 *      selection survives, which is what the formatting commands actually need.
 *      That is the primary fix.
 *   2. This function, for the controls that legitimately must take focus (a menu
 *      needs it for arrow-key navigation).
 *
 *  Nothing ever clears a stale record. That is intentional: the toolbar
 *  re-resolves the cell from its pane's store on read, which self-heals across
 *  cell deletion, notebook close, and a library load replacing every store --
 *  with no teardown hooks to forget. The one genuine clear is leaving focused
 *  mode (see setFocused in canvas.ts). */
export function markBlurred(cellId: string) {
  activeCell.update(a => (a && a.cellId === cellId ? { ...a, focused: false } : a));
}

/** The cell-style combo retyped the active cell under us. Belt-and-braces: the
 *  toolbar's validate-on-read already reads the type off the store. */
export function retypeActiveCell(cellType: CellType) {
  activeCell.update(a => (a ? { ...a, cellType } : a));
}

/** Forget the active cell entirely. Called when leaving focused mode, so
 *  returning to the canvas does not leave a phantom cell for the toolbar. */
export function clearActiveCell() {
  activeCell.set(null);
}
