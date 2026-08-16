<!--
  MenuBar.svelte — the application menu bar: File, Edit, Insert, Cell, Evaluation, Graphics.

  IN-APP, NOT NATIVE. A native macOS menu bar would mean the menu definitions live in Rust
  (`tauri::menu`) and every item round-trips an event into the webview to reach the store the
  command actually operates on. That split is worth paying for once the commands are settled;
  while they are still moving, one definition in one language that can call the store directly
  is the cheaper half, and it is also the half that works identically on Linux and Windows.

  WHY IT RENDERS ITS OWN <Menu>. Menu.svelte uses position:fixed, which does not escape an
  ancestor that creates a containing block via backdrop-filter -- and .nb-card has one. This
  component sits in the app bar, outside every card, so it is a legal render site (the same
  reason Toolbar.svelte renders the toolbar's menus rather than the cards doing it).

  DISABLED IS HONEST. Save, Open and the rest of the file layer are shown greyed rather than
  hidden, because a notebook file format does not exist yet: hiding them would misrepresent
  the menu as complete, and wiring them to something that half-works would be worse than
  either. Every enabled item calls a command that exists.
-->
<script lang="ts">
  import Menu from './Menu.svelte';
  import type { MenuItem } from './Menu.svelte';
  import { canvasState, activeActions, addNotebook, openRefpage } from './canvas';
  import { activeCell, activeHandle, retypeActiveCell } from './active';
  import { splitCell, mergeCellDown, duplicateCell, deleteCell,
           indentCode, outdentCode, commentCode, duplicateLine } from './cellCommands';
  import { restart, abortEvaluation } from './kernelActions';
  import { undo, redo } from '@codemirror/commands';

  /* Which top-level title is open, by name. null when the bar is at rest. */
  let openTitle: string | null = null;
  const anchors: Record<string, HTMLButtonElement> = {};

  $: focused = $canvasState.focusedIds.length > 0;
  $: act = $activeActions;
  $: cell = $activeCell;
  /* A code cell is the precondition for most of Edit and all of Cell: the caret has to be
     somewhere for "divide this" to mean anything. */
  $: hasCell = !!cell;
  $: hasCode = !!cell && cell.cellType === 'code';

  const SEP: MenuItem = { kind: 'sep' };
  const item = (id: string, label: string, hint?: string, disabled = false): MenuItem =>
    ({ kind: 'item', id, label, hint, disabled });

  /* The item sets. Built reactively so that "Divide Cell" is greyed when there is no caret
     rather than failing silently when chosen. */
  $: fileItems = [
    item('file.new', 'New Notebook', '⌘N'),
    item('file.open', 'Open…', '⌘O', true),
    item('file.openRecent', 'Open Recent', '', true),
    SEP,
    item('file.close', 'Close', '⌘W', !focused),
    item('file.save', 'Save', '⌘S', true),
    item('file.saveAs', 'Save As…', '⇧⌘S', true),
    SEP,
    item('file.print', 'Print…', '⌘P'),
  ] as MenuItem[];

  $: editItems = [
    item('edit.undo', 'Undo', '⌘Z', !hasCode),
    item('edit.redo', 'Redo', '⇧⌘Z', !hasCode),
    SEP,
    item('edit.cut', 'Cut', '⌘X', !hasCell),
    item('edit.copy', 'Copy', '⌘C', !hasCell),
    item('edit.selectAll', 'Select All', '⌘A', !hasCode),
    SEP,
    item('edit.comment', 'Un/Comment Selection', '⌘/', !hasCode),
    item('edit.indent', 'Indent Selected Lines', '', !hasCode),
    item('edit.outdent', 'Outdent Selected Lines', '', !hasCode),
    item('edit.dupLine', 'Duplicate Line', '⇧⌘D', !hasCode),
    SEP,
    item('edit.find', 'Find Symbol Documentation…', '', !hasCode),
  ] as MenuItem[];

  $: insertItems = [
    item('insert.code', 'Input Cell', '', !focused),
    item('insert.text', 'Text Cell', '', !focused),
    item('insert.section', 'Section Cell', '', !focused),
    SEP,
    item('insert.image', 'Image from File…', '', true),
  ] as MenuItem[];

  /* Convert To is flattened rather than nested: Menu.svelte has no submenu level, and three
     flat items read better than a submenu holding three. */
  $: cellItems = [
    item('cell.toInput', 'Convert to Input', '', !hasCell || cell?.cellType === 'code'),
    item('cell.toText', 'Convert to Text', '', !hasCell || cell?.cellType === 'text'),
    item('cell.toSection', 'Convert to Section', '', !hasCell || cell?.cellType === 'section'),
    SEP,
    item('cell.divide', 'Divide Cell', '⇧⌘D', !hasCell),
    item('cell.merge', 'Merge Cells', '⇧⌘M', !hasCell),
    item('cell.duplicate', 'Duplicate Cell', '', !hasCell),
    item('cell.delete', 'Delete Cell', '', !hasCell),
    SEP,
    item('cell.clearOutput', 'Delete Output', '', !hasCell),
    item('cell.clearAllOutput', 'Delete All Output', '', !focused),
  ] as MenuItem[];

  $: evalItems = [
    item('eval.cell', 'Evaluate Cell', '⇧↵', !hasCell),
    item('eval.notebook', 'Evaluate Notebook', '', !focused),
    SEP,
    item('eval.abort', 'Abort Evaluation', '', !focused),
    item('eval.restart', 'Restart Kernel', ''),
  ] as MenuItem[];

  /* Graphics is a documentation entry point rather than a set of editing commands: the
     renderer has no interactive object model to act on, so pretending otherwise with greyed
     drawing tools would be noise. */
  $: graphicsItems = [
    item('gfx.plot', 'Plot Documentation'),
    item('gfx.image', 'Image Documentation'),
    item('gfx.image3d', 'Image3D Documentation'),
    SEP,
    item('gfx.graphics', 'Graphics Documentation'),
  ] as MenuItem[];

  const TITLES = ['File', 'Edit', 'Insert', 'Cell', 'Evaluation', 'Graphics'];

  $: itemsFor = (t: string | null): MenuItem[] =>
    t === 'File'       ? fileItems :
    t === 'Edit'       ? editItems :
    t === 'Insert'     ? insertItems :
    t === 'Cell'       ? cellItems :
    t === 'Evaluation' ? evalItems :
    t === 'Graphics'   ? graphicsItems : [];

  function toggle(t: string) { openTitle = openTitle === t ? null : t; }
  /* Once one menu is open, sliding along the bar switches menus without a second click --
     the behaviour every desktop menu bar has. */
  function hover(t: string) { if (openTitle && openTitle !== t) openTitle = t; }

  /* The notebook the menu acts on: the active pane in focused mode, and in canvas mode
     whichever card was last active. */
  function firstNotebookId(): string | null {
    return $canvasState.focusedIds[0] ?? null;
  }

  function docFor(name: string) {
    /* openRefpage needs a notebook to open FROM -- it places the reference card next to it. On
       the canvas with nothing focused there is no such anchor, so the first existing notebook
       stands in; a canvas with no notebooks at all has nowhere to put a page and does nothing. */
    const from = firstNotebookId() ?? $canvasState.notebooks[0]?.id ?? null;
    if (from) openRefpage(from, name);
  }

  function onSelect(e: CustomEvent<{ id: string }>) {
    const id = e.detail.id;
    openTitle = null;
    const view = activeHandle()?.view ?? null;
    const c = cell;

    switch (id) {
      /* ---- File ---- */
      case 'file.new':    addNotebook(); break;
      case 'file.close':  act?.close(); break;
      case 'file.print':  window.print(); break;

      /* ---- Edit ---- */
      case 'edit.undo':   if (view) { undo(view); view.focus(); } break;
      case 'edit.redo':   if (view) { redo(view); view.focus(); } break;
      case 'edit.cut':    document.execCommand('cut'); break;
      case 'edit.copy':   document.execCommand('copy'); break;
      case 'edit.selectAll':
        if (view) {
          view.dispatch({ selection: { anchor: 0, head: view.state.doc.length } });
          view.focus();
        }
        break;
      case 'edit.comment':  if (view) commentCode(view); break;
      case 'edit.indent':   if (view) indentCode(view); break;
      case 'edit.outdent':  if (view) outdentCode(view); break;
      case 'edit.dupLine':  if (view) duplicateLine(view); break;
      case 'edit.find': {
        /* The word under the caret is the query -- the same thing the toolbar's docs button
           resolves, reached from the keyboard-first path. */
        const sel = view ? view.state.sliceDoc(view.state.selection.main.from,
                                               view.state.selection.main.to).trim() : '';
        const word = sel || 'Image';
        docFor(word.replace(/[^A-Za-z0-9$]/g, '') || 'Image');
        break;
      }

      /* ---- Insert ---- */
      case 'insert.code':    insertCell('code'); break;
      case 'insert.text':    insertCell('text'); break;
      case 'insert.section': insertCell('section'); break;

      /* ---- Cell ---- */
      case 'cell.toInput':   retypeActiveCell('code'); break;
      case 'cell.toText':    retypeActiveCell('text'); break;
      case 'cell.toSection': retypeActiveCell('section'); break;
      case 'cell.divide': {
        if (!act || !c) break;
        /* The caret offset exists only while the editor holds focus, which the title button's
           suppressed pointerdown preserves. Splitting at the end beats refusing when it does not. */
        const h = activeHandle();
        const at = c.focused && h?.view ? h.view.state.selection.main.head : null;
        const nid = splitCell(act.store, c.cellId, at);
        if (nid) act.focusCell(nid);
        break;
      }
      case 'cell.merge':     if (act && c) mergeCellDown(act.store, c.cellId); break;
      case 'cell.duplicate': if (act && c) duplicateCell(act.store, c.cellId); break;
      case 'cell.delete':    if (act && c) deleteCell(act.store, c.cellId); break;
      case 'cell.clearOutput': if (act && c) act.store.clearOutput(c.cellId); break;
      case 'cell.clearAllOutput':
        if (act) for (const cc of act.store.allCells()) act.store.clearOutput(cc.id);
        break;

      /* ---- Evaluation ---- */
      case 'eval.cell':     if (act && c) act.runCell(c.cellId); break;
      case 'eval.notebook': act?.runAll(); break;
      case 'eval.abort':    abortEvaluation(); break;
      case 'eval.restart':  restart(); break;

      /* ---- Graphics ---- */
      case 'gfx.plot':     docFor('Plot'); break;
      case 'gfx.image':    docFor('Image'); break;
      case 'gfx.image3d':  docFor('Image3D'); break;
      case 'gfx.graphics': docFor('Graphics'); break;
    }
  }

  /* A new cell goes after the row holding the caret, and at the end when there is no caret.
     Rows rather than cells because a row can hold several cells side by side, and "insert after
     the caret" means after the whole row it is in. */
  function insertCell(type: 'code' | 'text' | 'section') {
    if (!act) return;
    const rows = act.store.getRows();
    const ri = cell ? rows.findIndex(r => r.cells.some(x => x.id === cell!.cellId)) : -1;
    const newId = ri >= 0 ? act.store.insertRowAt(ri + 1, type) : act.store.addRow(type);
    if (newId) act.focusCell(newId);
  }
</script>

<div class="menubar" role="menubar">
  {#each TITLES as t}
    <button
      class="mb-title"
      class:open={openTitle === t}
      role="menuitem"
      aria-haspopup="menu"
      aria-expanded={openTitle === t}
      tabindex="-1"
      bind:this={anchors[t]}
      on:pointerdown|preventDefault
      on:click={() => toggle(t)}
      on:pointerenter={() => hover(t)}
    >{t}</button>
  {/each}
</div>

<Menu
  open={openTitle !== null}
  items={itemsFor(openTitle)}
  anchor={openTitle ? anchors[openTitle] : null}
  align="start"
  minWidth={210}
  on:select={onSelect}
  on:close={() => (openTitle = null)}
/>

<style>
  .menubar {
    display: flex;
    align-items: center;
    gap: 1px;
    height: var(--menubar-h, 28px);
    flex: 0 0 auto;
  }
  .mb-title {
    height: 22px;
    padding: 0 9px;
    font: inherit;
    font-size: 12px;
    color: var(--text);
    background: none;
    border: 0;
    border-radius: 4px;
    cursor: default;
    white-space: nowrap;
  }
  .mb-title:hover { background: var(--surface-2); }
  .mb-title.open  { background: var(--surface-3); color: var(--accent); }
</style>
