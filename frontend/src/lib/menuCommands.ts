/* menuCommands.ts -- one implementation behind every menu command.
 *
 * The menu bar itself is NATIVE (built in src-tauri/src/lib.rs, so on macOS it lives in the system
 * bar at the top of the screen and the window's own top strip is free for other things). A native
 * item cannot touch the notebook store directly: it fires in Rust, which emits `menu:<id>`, which
 * the webview receives. So the commands live here, keyed by exactly the ids the Rust menu uses, and
 * the listener in App.svelte is a loop rather than a growing list of one-line handlers.
 *
 * WHY A TABLE AND NOT A SWITCH IN App.svelte. The ids are a contract with the Rust side. Keeping
 * them in one exported list means the wiring can be CHECKED -- MENU_IDS is what App.svelte
 * subscribes to, and an id in the Rust menu with no entry here shows up as a menu item that does
 * nothing, which is the failure this arrangement makes findable instead of silent.
 *
 * A few commands need things only App.svelte has (the file dialogs, which are its own state), so
 * they arrive as hooks rather than being reached for through a module-level singleton.
 */

import { get } from 'svelte/store';
import { canvasState, activeActions, addNotebook, openRefpage } from './canvas';
import { activeCell, activeHandle, retypeActiveCell } from './active';
import { splitCell, mergeCellDown, duplicateCell, deleteCell, convertCell,
         indentCode, outdentCode, commentCode, duplicateLine,
         copyInputFromAbove } from './cellCommands';
import { restart, abortEvaluation } from './kernelActions';
import { darkMode } from './theme';
import type { CellType } from './notebook';

/** What App.svelte lends the dispatcher: the library-level file operations it owns. */
export interface MenuHooks {
  openFile: () => void;
  saveFile: () => void;
  saveFileAs: () => void;
  /** Open entry `i` of the recent-files list (the `i` in the `recent-<i>` menu id). */
  openRecent: (i: number) => void;
  clearRecent: () => void;
}

/* Every id the native menu can emit. App.svelte subscribes to `menu:<id>` for each.
 *
 * HYPHENS, NOT DOTS. An id becomes the Tauri event name `menu:<id>`, and Tauri accepts only
 * `[A-Za-z0-9-/:_]` there -- a `.` makes `listen()` THROW. Every id here was once dotted
 * (`file.new`), which killed the whole menu bar: App.svelte subscribed in one loop, so the first
 * throw aborted it and the legal ids after it were never reached either. `make check-menu-ids`
 * now validates the grammar as well as the two-sided contract, because the old check compared the
 * Rust and TypeScript sides against each other and both agreed on a name Tauri rejects. */
export const MENU_IDS = [
  'file-new', 'open', 'save', 'save-as', 'file-close', 'file-print',
  /* File > Open Recent. A FIXED pool of ten slot ids, spelled out here and in lib.rs rather
     than generated, because check_menu_ids.py reads both sides as source text — a
     `format!("recent-{i}")` would be invisible to it and these would read as dead wiring.
     The digit is an INDEX into the list recentFiles() returns, so the run must stay
     contiguous from 0; the gate checks that too, since a gap would silently open the wrong
     file rather than fail. */
  'recent-0', 'recent-1', 'recent-2', 'recent-3', 'recent-4',
  'recent-5', 'recent-6', 'recent-7', 'recent-8', 'recent-9',
  'recent-clear',
  'edit-comment', 'edit-indent', 'edit-outdent', 'edit-dupLine', 'edit-findDoc',
  'edit-copyInputAbove',
  'insert-code', 'insert-text', 'insert-section',
  'cell-toInput', 'cell-toText',
  'cell-toTitle', 'cell-toSubtitle', 'cell-toChapter',
  'cell-toSection', 'cell-toSubsection', 'cell-toSubsubsection',
  'cell-divide', 'cell-merge', 'cell-duplicate', 'cell-delete',
  'cell-clearOutput', 'cell-clearAllOutput',
  'eval-cell', 'run-all', 'interrupt', 'restart',
  'gfx-plot', 'gfx-image', 'gfx-image3d', 'gfx-graphics',
  'toggle-dark',
] as const;

/* The notebook a command acts on: the active pane, or the first one on the canvas so that a
   documentation command still has somewhere to open a page beside. */
function anchorNotebookId(): string | null {
  const s = get(canvasState);
  return s.focusedIds[0] ?? s.notebooks[0]?.id ?? null;
}

function docFor(name: string) {
  const from = anchorNotebookId();
  if (from) openRefpage(from, name);
}

/* A new cell goes after the row holding the caret, and at the end when there is no caret. Rows
   rather than cells because a row can hold several cells side by side. */
function insertCell(type: 'code' | 'text' | 'section') {
  const act = get(activeActions);
  if (!act) return;
  const cell = get(activeCell);
  const rows = act.store.getRows();
  const ri = cell ? rows.findIndex(r => r.cells.some(x => x.id === cell.cellId)) : -1;
  const newId = ri >= 0 ? act.store.insertRowAt(ri + 1, type) : act.store.addRow(type);
  if (newId) act.focusCell(newId);
}

/* Convert the cell the caret is in. The record in `activeCell` only says WHICH cell; the
   conversion itself is convertCell on the pane's store, exactly what the toolbar's cell-style
   control calls. Calling retypeActiveCell alone -- as these items once did -- relabelled the
   toolbar and left the cell untouched. The notebook-id check is the toolbar's own: in split mode
   a remembered cell from another pane is not this pane's to convert. */
export function convertActiveCell(type: CellType) {
  const act = get(activeActions);
  const cell = get(activeCell);
  if (!act || !cell || cell.notebookId !== act.notebookId) return;
  if (convertCell(act.store, cell.cellId, type)) retypeActiveCell(type);
}

export function runMenuCommand(id: string, hooks: MenuHooks) {
  const act = get(activeActions);
  const cell = get(activeCell);
  const view = activeHandle()?.view ?? null;

  switch (id) {
    /* ---- File ---- */
    case 'file-new':    addNotebook(); break;
    case 'open':        hooks.openFile(); break;
    case 'save':        hooks.saveFile(); break;
    case 'save-as':     hooks.saveFileAs(); break;
    case 'file-close':  act?.close(); break;
    case 'file-print':  window.print(); break;

    /* File > Open Recent. Ten literal cases rather than a prefix test on `id`, because the
       gate reads `case '<id>':` out of this file as text: a `if (id.startsWith('recent-'))`
       would leave all ten looking unhandled and the check would report ten dead menu items.
       The digit is the index into the list — see MENU_IDS. */
    case 'recent-0': hooks.openRecent(0); break;
    case 'recent-1': hooks.openRecent(1); break;
    case 'recent-2': hooks.openRecent(2); break;
    case 'recent-3': hooks.openRecent(3); break;
    case 'recent-4': hooks.openRecent(4); break;
    case 'recent-5': hooks.openRecent(5); break;
    case 'recent-6': hooks.openRecent(6); break;
    case 'recent-7': hooks.openRecent(7); break;
    case 'recent-8': hooks.openRecent(8); break;
    case 'recent-9': hooks.openRecent(9); break;
    case 'recent-clear': hooks.clearRecent(); break;

    /* ---- Edit ----
       Undo, redo, cut, copy, paste and select-all are PREDEFINED native items with no id: they
       route through the macOS responder chain into the focused editor, which CodeMirror already
       binds. Handling them here as well would be dead code pretending to be wiring -- the events
       never arrive, because the items do not emit any. Only the commands with no native
       equivalent live below. */
    case 'edit-comment':  if (view) commentCode(view); break;
    case 'edit-indent':   if (view) indentCode(view); break;
    case 'edit-outdent':  if (view) outdentCode(view); break;
    case 'edit-dupLine':  if (view) duplicateLine(view); break;
    /* Needs the store as well as the view — it reads the cell above, not the selection. */
    case 'edit-copyInputAbove':
      if (act && cell && view) copyInputFromAbove(act.store, cell.cellId, view);
      break;
    case 'edit-findDoc': {
      /* The selected word, or the symbol the caret sits in, opened as its own reference page. */
      const sel = view
        ? view.state.sliceDoc(view.state.selection.main.from, view.state.selection.main.to).trim()
        : '';
      const word = (sel || 'Image').replace(/[^A-Za-z0-9$]/g, '');
      docFor(word || 'Image');
      break;
    }

    /* ---- Insert ---- */
    case 'insert-code':    insertCell('code'); break;
    case 'insert-text':    insertCell('text'); break;
    case 'insert-section': insertCell('section'); break;

    /* ---- Cell ---- */
    case 'cell-toInput':          convertActiveCell('code'); break;
    case 'cell-toText':           convertActiveCell('text'); break;
    case 'cell-toTitle':          convertActiveCell('title'); break;
    case 'cell-toSubtitle':       convertActiveCell('subtitle'); break;
    case 'cell-toChapter':        convertActiveCell('chapter'); break;
    case 'cell-toSection':        convertActiveCell('section'); break;
    case 'cell-toSubsection':     convertActiveCell('subsection'); break;
    case 'cell-toSubsubsection':  convertActiveCell('subsubsection'); break;
    case 'cell-divide': {
      if (!act || !cell) break;
      /* The caret offset exists only while the editor holds focus; splitting at the end beats
         refusing when a native menu has just taken it. */
      const h = activeHandle();
      const at = cell.focused && h?.view ? h.view.state.selection.main.head : null;
      const nid = splitCell(act.store, cell.cellId, at);
      if (nid) act.focusCell(nid);
      break;
    }
    case 'cell-merge':     if (act && cell) mergeCellDown(act.store, cell.cellId); break;
    case 'cell-duplicate': if (act && cell) duplicateCell(act.store, cell.cellId); break;
    case 'cell-delete':    if (act && cell) deleteCell(act.store, cell.cellId); break;
    case 'cell-clearOutput': if (act && cell) act.store.clearOutput(cell.cellId); break;
    case 'cell-clearAllOutput':
      if (act) for (const c of act.store.allCells()) act.store.clearOutput(c.id);
      break;

    /* ---- Evaluation ---- */
    case 'eval-cell':   if (act && cell) act.runCell(cell.cellId); break;
    case 'run-all':     act?.runAll(); break;
    case 'interrupt':   abortEvaluation(); break;
    case 'restart':     restart(); break;

    /* ---- Graphics: documentation entry points. Mathilda's own pages, never an external site. */
    case 'gfx-plot':     docFor('Plot'); break;
    case 'gfx-image':    docFor('Image'); break;
    case 'gfx-image3d':  docFor('Image3D'); break;
    case 'gfx-graphics': docFor('Graphics'); break;

    /* ---- View ---- */
    case 'toggle-dark': darkMode.update(v => !v); break;

    default:
      /* An id the Rust menu emits with no case here: a menu item that does nothing. Say so once
         rather than failing silently. */
      console.warn(`menu command "${id}" has no handler`);
  }
}
