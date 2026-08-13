<!--
  NotebookCard.svelte
  One draggable glass-dark notebook card on the infinite canvas.

  Props:
    nb          — CanvasNotebook from canvas.ts
    currentZoom — the current animated zoom level (for drag delta correction)

  Responsibilities:
    • Drag title bar → update nb.x/nb.y via canvasState (delta / currentZoom)
    • Glass-dark aesthetic with per-notebook colour accent
    • Collapse toggle (⊞) and close (✕) buttons in title bar
    • Double-click title to rename inline
    • Expanded: full notebook cell UI (rows of cells via CellShell)
    • Collapsed: compact card with "N cells" and last output preview
    • Mount animation: scale(0.94)/opacity(0) → scale(1)/opacity(1), 180ms
    • All cell execution logic self-contained
    • kernelStatus store from notebook.ts for status display
-->
<script lang="ts">
  import { onMount, onDestroy, tick, createEventDispatcher } from 'svelte';
  import { setFocused, setNotebookWidth, setNotebookHeight, addQueryNotebook, openRefpage,
           registerPane, publishPaneFlags, retractPane } from './canvas';
  const dispatch = createEventDispatcher();
  import { get } from 'svelte/store';
  import CellShell from './CellShell.svelte';
  import type { CanvasNotebook } from './canvas';
  import {
    canvasState,
    removeNotebook,
    toggleCollapse,
    renameNotebook,
  } from './canvas';
  import {
    kernelStatus,
    selectedCells,
    clearSelection,
  } from './notebook';
  import type { OutputItem, CellType } from './notebook';
  import {
    evaluateCell,
  } from './ipc';
  import type { OutputMessage } from './ipc';

  export let nb: CanvasNotebook;
  export let currentZoom: number = 1.0;
  export let focused: boolean = false;  // true when rendered full-screen
  export let isSelected: boolean = false;

  // Per-notebook accent colour from the deep-space palette
  const PALETTE = ['#89b4fa','#a6e3a1','#f38ba8','#fab387','#cba6f7','#94e2d5'];
  $: accentColor = PALETTE[parseInt(nb.id.replace('nb-', ''), 10) % PALETTE.length] ?? '#89b4fa';

  // ---- Store subscription ----
  // nb.store is a Svelte store (has .subscribe). Assign it to a local `let`
  // so Svelte 4's `$nbStore` auto-subscription works in the template.
  let nbStore = nb.store;
  $: nbStore = nb.store;

  /* Is this a reference page? Derived from the content as well as the flag:
     a card opened before the flag existed, or restored from a saved library
     (which does not persist it), still has 'ref' cells and should behave as
     documentation -- non-editable headings, no cell-type control. */
  $: isRefpage = !!nb.refpage || $nbStore.some(r => r.cells[0]?.type === 'ref');

  // ---------------------------------------------------------------------------
  // Mount animation

  let mounted = false;
  onMount(() => {
    requestAnimationFrame(() => { mounted = true; });
    /* A card opened to answer a query arrives with its cell already filled;
     * evaluate it once so the answer is there when the card appears. */
    if (nb.autoRun) {
      nb.autoRun = false;
      tick().then(() => {
        const rows = get(nb.store);
        const first = rows.find(r => r.cells[0]?.source.trim());
        if (first) runCell(first.cells[0].id, first.cells[0].source);
      });
    }
  });

  // ---------------------------------------------------------------------------
  // Title bar drag → move card on canvas

  let dragging   = false;
  let dragStartX = 0;
  let dragStartY = 0;
  let dragNbX0   = 0;
  let dragNbY0   = 0;

  let cardEl: HTMLElement;

  let dragMoved = false;

  // ---------------------------------------------------------------------------
  // Right-edge resize handle

  let resizing = false;
  let resizeStartX = 0;
  let resizeStartW = 0;

  function onResizePointerDown(e: PointerEvent) {
    if (e.pointerType === "touch" || e.pointerType === "pen") return;  // finger: pan the canvas, do not resize
    if (e.button !== 0) return;
    e.stopPropagation();
    resizing = true;
    resizeStartX = e.clientX;
    resizeStartW = nb.width;
    (e.currentTarget as HTMLElement).setPointerCapture(e.pointerId);
  }

  function onResizePointerMove(e: PointerEvent) {
    if (!resizing) return;
    e.stopPropagation();
    const dx = (e.clientX - resizeStartX) / (currentZoom || 1);
    setNotebookWidth(nb.id, resizeStartW + dx);
  }

  function onResizePointerUp(_e: PointerEvent) { resizing = false; }

  // ---------------------------------------------------------------------------
  // Bottom / corner resize handles

  let resizingBottom  = false;
  let resizingCorner  = false;
  let resizeStartY    = 0;
  let resizeStartH    = 0;

  const TITLE_BAR_H = 36;

  function currentHeight(): number {
    return nb.height ?? (cardEl?.offsetHeight ?? 400);
  }

  function onBottomResizeDown(e: PointerEvent) {
    if (e.pointerType === "touch" || e.pointerType === "pen") return;  // finger: pan the canvas, do not resize
    if (e.button !== 0) return;
    e.stopPropagation();
    resizingBottom = true;
    resizeStartY   = e.clientY;
    resizeStartH   = currentHeight();
    resizeStartX   = e.clientX;
    resizeStartW   = nb.width;
    (e.currentTarget as HTMLElement).setPointerCapture(e.pointerId);
  }

  function onCornerResizeDown(e: PointerEvent) {
    if (e.pointerType === "touch" || e.pointerType === "pen") return;  // finger: pan the canvas, do not resize
    if (e.button !== 0) return;
    e.stopPropagation();
    resizingCorner = true;
    resizeStartY   = e.clientY;
    resizeStartH   = currentHeight();
    resizeStartX   = e.clientX;
    resizeStartW   = nb.width;
    (e.currentTarget as HTMLElement).setPointerCapture(e.pointerId);
  }

  function onBottomResizeMove(e: PointerEvent) {
    if (!resizingBottom && !resizingCorner) return;
    e.stopPropagation();
    const dy = (e.clientY - resizeStartY) / (currentZoom || 1);
    setNotebookHeight(nb.id, resizeStartH + dy);
    if (resizingCorner) {
      const dx = (e.clientX - resizeStartX) / (currentZoom || 1);
      setNotebookWidth(nb.id, resizeStartW + dx);
    }
  }

  function onBottomResizeUp(_e: PointerEvent) {
    resizingBottom = false;
    resizingCorner = false;
  }

  function onTitlePointerDown(e: PointerEvent) {
    // Finger: no card dragging — a touch on the titlebar should pan the canvas
    // (the event bubbles to Canvas). A mouse still drags (even on a touchscreen
    // laptop). Titlebar buttons keep working via click.
    if (e.pointerType === "touch" || e.pointerType === "pen") return;
    if (e.button !== 0) return;  // left-click drag only
    if (focused) return;         // no drag in full-screen mode
    if ((e.target as HTMLElement).closest('button, input')) return;
    e.stopPropagation();
    dragging   = true;
    dragMoved  = false;
    dragStartX = e.clientX;
    dragStartY = e.clientY;
    dragNbX0   = nb.x;
    dragNbY0   = nb.y;
    cardEl.setPointerCapture(e.pointerId);
  }

  function onTitlePointerMove(e: PointerEvent) {
    if (!dragging) return;
    e.stopPropagation();
    const dxScreen = e.clientX - dragStartX;
    const dyScreen = e.clientY - dragStartY;
    if (Math.abs(dxScreen) > 4 || Math.abs(dyScreen) > 4) dragMoved = true;
    const effectiveZoom = currentZoom || 1;
    const dxWorld = dxScreen / effectiveZoom;
    const dyWorld = dyScreen / effectiveZoom;
    if (isSelected) {
      // Move all selected notebooks by the same world-space delta
      dispatch('groupMove', { dx: dxWorld, dy: dyWorld, originX: dragNbX0, originY: dragNbY0, id: nb.id });
    } else {
      canvasState.update(s => ({
        ...s,
        notebooks: s.notebooks.map(n => n.id === nb.id ? { ...n, x: dragNbX0 + dxWorld, y: dragNbY0 + dyWorld } : n),
      }));
    }
  }

  function onTitlePointerUp(_e: PointerEvent) {
    if (dragging && isSelected) dispatch('groupMoveEnd', {});
    dragging = false;
  }

  // ---------------------------------------------------------------------------
  // Inline rename

  let renaming = false;
  let renameInput: HTMLInputElement;
  let renameValue = '';

  function startRename() {
    renaming    = true;
    renameValue = nb.title;
    tick().then(() => {
      if (!renameInput) return;
      renameInput.focus();
      // Select all text in the contenteditable
      const range = document.createRange();
      range.selectNodeContents(renameInput);
      const sel = window.getSelection();
      if (sel) { sel.removeAllRanges(); sel.addRange(range); }
    });
  }

  function commitRename() {
    const text = renameInput?.innerText?.trim() ?? renameValue.trim();
    if (text) renameNotebook(nb.id, text);
    renaming = false;
  }

  function onRenameKeydown(e: KeyboardEvent) {
    if (e.key === 'Enter') { e.preventDefault(); commitRename(); }
    else if (e.key === 'Escape') { renaming = false; }
  }

  // ---------------------------------------------------------------------------
  // Collapse toggle

  function onToggleCollapse() {
    toggleCollapse(nb.id);
  }

  // ---------------------------------------------------------------------------
  // Cell focus registry

  const cellFocusFns: Record<string, () => void> = {};

  function handleRegister(e: CustomEvent<{ id: string; fn: () => void }>) {
    cellFocusFns[e.detail.id] = e.detail.fn;
  }

  // ---------------------------------------------------------------------------
  // Insertion point — blinking horizontal cursor between cells
  // 0 = before first row, i = between row[i-1] and row[i]

  let insertionIdx: number | null = null;

  function handleFocusPrev(e: CustomEvent<{ id: string }>) {
    const rows = nb.store.getRows();
    const ri = rows.findIndex((r: any) => r.cells.some((c: any) => c.id === e.detail.id));
    if (ri >= 0) { insertionIdx = ri; tick().then(() => cardEl?.focus()); }
  }

  function handleFocusNext(e: CustomEvent<{ id: string }>) {
    const rows = nb.store.getRows();
    const ri = rows.findIndex((r: any) => r.cells.some((c: any) => c.id === e.detail.id));
    if (ri >= 0) { insertionIdx = ri + 1; tick().then(() => cardEl?.focus()); }
  }

  // ---------------------------------------------------------------------------
  // Wheel: scroll card when an editor inside is focused; pan canvas otherwise.

  function onCardBodyWheel(e: WheelEvent) {
    // All wheel events (pan and pinch) pass through to the canvas.
    // The card body's native overflow-y:auto scroll still works when focused
    // because the browser handles it before the event reaches the canvas handler.
  }

  // ---------------------------------------------------------------------------
  // Section collapse — track which section rows are collapsed

  // Layout toggle: vertical (default) vs horizontal (input|output side-by-side)
  let horizontal = false;

  let collapsedSections = new Set<string>();

  /* Every heading row, and whether any is currently open. One button rather
     than two: it collapses everything while anything is open, and expands
     everything once all are closed. */
  $: sectionRowIds = $nbStore
    .filter(r => r.cells[0]?.type === 'section' || r.cells[0]?.type === 'subsection')
    .map(r => r.id);
  $: allSectionsCollapsed =
    sectionRowIds.length > 0 && sectionRowIds.every(id => collapsedSections.has(id));

  function toggleAllSections() {
    collapsedSections = allSectionsCollapsed ? new Set() : new Set(sectionRowIds);
  }

  function toggleSection(rowId: string) {
    collapsedSections = collapsedSections.has(rowId)
      ? new Set([...collapsedSections].filter(id => id !== rowId))
      : new Set([...collapsedSections, rowId]);
  }

  // Reactive set of hidden row IDs — depends on BOTH $nbStore AND collapsedSections.
  // Plain function (not isHidden()) so the $: makes Svelte track collapsedSections,
  // ensuring {#if !hiddenRows.has(row.id)} re-evaluates immediately on toggle.
  /* How deeply each row sits under a heading, so content reads as belonging to
     its section rather than sitting flush with it. 0 for a section heading,
     1 under a section (and for a subsection heading), 2 under a subsection. */
  $: rowDepth = (() => {
    const depth = new Map<string, number>();
    const rows = $nbStore;
    let cur = 0;
    for (const row of rows) {
      const t = row.cells[0]?.type;
      if (t === 'section') { depth.set(row.id, 0); cur = 1; }
      else if (t === 'subsection') { depth.set(row.id, 1); cur = 2; }
      else depth.set(row.id, cur);
    }
    return depth;
  })();

  $: hiddenRows = (() => {
    const hidden = new Set<string>();
    const rows = $nbStore;
    for (let ri = 0; ri < rows.length; ri++) {
      const currentType = rows[ri]?.cells[0]?.type ?? 'code';
      let hide = false;
      /* Only the NEAREST subsection above a row is its parent. Any subsection
       * before that one is a SIBLING, and its fold state says nothing about
       * this row -- but the scan used to test every subsection it walked past,
       * so collapsing "Basic examples" also emptied "Options" and
       * "Applications" below it. Once the parent has been seen and is open, keep
       * walking only to find the enclosing section. */
      let parentSubsectionSeen = false;
      for (let i = ri - 1; i >= 0; i--) {
        const cell = rows[i]?.cells[0];
        if (!cell) continue;
        if (cell.type === 'section') {
          if (currentType !== 'section') hide = collapsedSections.has(rows[i].id);
          break; // stop at first section boundary
        }
        if (cell.type === 'subsection') {
          // A section is never nested inside a subsection.
          if (currentType === 'section') break;
          // A subsection's parent is the SECTION, not the subsection before it.
          // Breaking here meant a subsection heading that follows a sibling
          // never asked whether its section was collapsed, so it stayed visible
          // with all its content hidden -- a heading that could not be opened.
          if (currentType === 'subsection') continue;
          if (!parentSubsectionSeen) {
            parentSubsectionSeen = true;
            if (collapsedSections.has(rows[i].id)) { hide = true; break; }
          }
          // A sibling, or an open parent: keep scanning for the section.
        }
      }
      if (hide) hidden.add(rows[ri].id);
    }
    return hidden;
  })();

  // ---------------------------------------------------------------------------
  // 4-directional row/cell add

  async function addRowAbove(e: CustomEvent<{ rowId: string }>) {
    const rowList = nb.store.getRows();
    const ri = rowList.findIndex((r: any) => r.id === e.detail.rowId);
    if (ri < 0) return;
    const id = nb.store.insertRowAt(ri);
    await tick(); cellFocusFns[id]?.();
  }

  async function addRowBelow(e: CustomEvent<{ rowId: string }>) {
    const rowList = nb.store.getRows();
    const ri = rowList.findIndex((r: any) => r.id === e.detail.rowId);
    const id = nb.store.insertRowAt(ri + 1);
    await tick(); cellFocusFns[id]?.();
  }

  /* `mathilda-refpage` is a custom DOM event bubbled up from the name grid in
     Output.svelte (and from refpages.ts's document-level Cmd+click and F1
     handlers), so the grid needs no prop threaded down through CellShell.
     Attached as an action rather than `on:mathilda-refpage` because Svelte's
     element typings only know standard DOM events. */
  function nameQueryListener(node: HTMLElement) {
    const handler = (e: Event) => onNameQuery(e);
    node.addEventListener('mathilda-refpage', handler);
    return { destroy() { node.removeEventListener('mathilda-refpage', handler); } };
  }

  /* A name in a `?pat*` result was clicked: open its documentation as a new
     card to the right. Appending a cell to this notebook instead put the answer
     at the very bottom, which in a long notebook is off-screen and reads as the
     click having done nothing. */
  async function onNameQuery(e: Event) {
    const detail = (e as CustomEvent<{ name: string; at?: { x: number; y: number } | null }>).detail;
    const name = detail?.name;
    if (!name) return;
    /* The full generated reference page, not `?name`: the terse docstring is
       what the ? query already gives, and a click asking for documentation
       should land on the examples and notes. */
    openRefpage(nb.id, name, detail?.at);
  }

  async function addCellLeft(e: CustomEvent<{ rowId: string; cellIdx: number }>) {
    const id = nb.store.insertCellInRow(e.detail.rowId, e.detail.cellIdx);
    await tick(); cellFocusFns[id]?.();
  }

  async function addCellRight(e: CustomEvent<{ rowId: string; cellIdx: number }>) {
    const id = nb.store.insertCellInRow(e.detail.rowId, e.detail.cellIdx);
    await tick(); cellFocusFns[id]?.();
  }

  function addRow(type: CellType = 'code') {
    const id = nb.store.addRow(type);
    tick().then(() => {
      // CodeMirror needs an extra frame; contenteditable needs focus() directly
      const fn = cellFocusFns[id];
      if (fn) { fn(); } else {
        setTimeout(() => cellFocusFns[id]?.(), 60);
      }
    });
  }

  // ---------------------------------------------------------------------------
  // Cell execution

  async function runAll() {
    nb.store.resetExecCounter();
    const cells = nb.store.allCells();
    for (const cell of cells) {
      if (cell.type === 'code' && cell.source.trim()) {
        await runCell(cell.id, cell.source);
      }
    }
  }

  /** Run an inclusive index range over allCells(). Backs the toolbar's
   *  "evaluate from here" and "evaluate above".
   *
   *  Deliberately does NOT reset the exec counter the way runAll does: a partial
   *  run should not renumber the cells it did not touch. */
  async function runRange(fromIdx: number, toIdx: number) {
    const cells = nb.store.allCells();
    const lo = Math.max(0, fromIdx);
    const hi = Math.min(cells.length - 1, toIdx);
    for (let i = lo; i <= hi; i++) {
      const cell = cells[i];
      if (cell?.type === 'code' && cell.source.trim()) {
        await runCell(cell.id, cell.source);
      }
    }
  }

  /** Run one cell by id — the toolbar has an id, not a source string. */
  async function runCellById(cellId: string) {
    const cell = nb.store.allCells().find((c: any) => c.id === cellId);
    if (cell) await runCell(cell.id, cell.source);
  }

  async function runCell(cellId: string, source: string) {
    if (!source.trim()) return;
    nb.store.stampExec(cellId);
    nb.store.clearOutput(cellId);
    nb.store.setStatus(cellId, 'running');
    kernelStatus.set('busy');
    try {
      await evaluateCell(source, (msg: OutputMessage) => {
        const item = msgToOutputItem(msg);
        if (!item) return;
        if (msg.type === 'stream') nb.store.appendStream(cellId, (msg as any).text ?? '');
        else nb.store.appendOutput(cellId, item);
      });
      nb.store.setStatus(cellId, 'done');
    } catch (err) {
      nb.store.appendOutput(cellId, { kind: 'error', text: String(err) });
      nb.store.setStatus(cellId, 'error');
      kernelStatus.set('dead');
    } finally {
      if (get(kernelStatus) !== 'dead') kernelStatus.set('ready');
    }
  }

  function msgToOutputItem(msg: OutputMessage): OutputItem | null {
    switch (msg.type) {
      case 'expr':   return { kind: 'expr', text: msg.payload, latex: (msg as any).latex };
      case 'usage':  return { kind: 'usage',  text: msg.payload };
      case 'names':  return { kind: 'names',  names: (msg as any).payload ?? [] };
      case 'error':  return { kind: 'error',  text: msg.message };
      case 'stream': return { kind: 'stream', text: (msg as any).text ?? '' };
      case 'plot':   return { kind: 'plot',   data: msg.payload };
      case 'html':   return { kind: 'html',   html: (msg as any).payload ?? '' };
      default:       return null;
    }
  }

  function handleRun(e: CustomEvent<{ id: string }>) {
    const cell = nb.store.allCells().find((c: any) => c.id === e.detail.id);
    if (cell) runCell(cell.id, cell.source);
  }

  function handleChange(e: CustomEvent<{ id: string; source: string }>) {
    nb.store.updateSource(e.detail.id, e.detail.source);
  }

  // ---------------------------------------------------------------------------
  // Collapsed preview: first 60 chars of last output expression

  $: lastOutputPreview = (() => {
    const allRows = $nbStore;
    if (!allRows || !allRows.length) return '';
    const cells = allRows.flatMap((r: any) => r.cells ?? []);
    const last = [...cells].reverse().find((c: any) => c.output?.length > 0);
    if (!last) return '';
    const item = last.output.find((o: any) => o.kind === 'expr' || o.kind === 'stream');
    const text = item?.text ?? '';
    return text.length > 60 ? text.slice(0, 60) + '…' : text;
  })();

  $: totalCells = (() => {
    const allRows = $nbStore;
    if (!allRows || !allRows.length) return 0;
    return allRows.reduce((acc: number, r: any) => acc + (r.cells?.length ?? 0), 0);
  })();

  // ---------------------------------------------------------------------------
  // Keyboard (scoped to card — stop propagation so canvas doesn't eat keys)

  async function onCardKeydown(e: KeyboardEvent) {
    // Let Cmd+0 and Cmd+N pass through to Canvas for fit-all / add notebook
    // Let ALL Cmd/Ctrl combos pass through to App.svelte's window handler
    // (Cmd+0 fit-all, Cmd+N new notebook, Cmd+=/- ui scale, Cmd+S/O save/open)
    if (e.metaKey || e.ctrlKey) return;
    e.stopPropagation();

    // If the event came from a CHILD element (e.g. CodeMirror editor), the
    // editor's keymap already handled it and fired handleFocusPrev/Next which
    // set insertionIdx. Don't immediately process navigation on the SAME event —
    // just let the cursor appear. The user must press again (while the card
    // itself has focus) to navigate further.
    if (insertionIdx === null || e.target !== cardEl) return;

    if (e.key === 'Escape') { e.preventDefault(); insertionIdx = null; return; }
    if (e.key === 'ArrowUp') {
      e.preventDefault();
      const rows = nb.store.getRows();
      // Arrow up from insertion point → enter the cell ABOVE the cursor
      // insertionIdx N = gap between row[N-1] and row[N].
      // "Above" = row[N-1], last cell in that row.
      if (insertionIdx > 0) {
        const targetRow = rows[insertionIdx - 1];
        if (targetRow) {
          insertionIdx = null;
          const lastCell = targetRow.cells[targetRow.cells.length - 1];
          if (lastCell) cellFocusFns[lastCell.id]?.();
        } else {
          insertionIdx = Math.max(0, insertionIdx - 1);
        }
      } else {
        // Already at top — dismiss cursor
        insertionIdx = null;
      }
      return;
    }
    if (e.key === 'ArrowDown') {
      e.preventDefault();
      const rows = nb.store.getRows();
      // Arrow down from insertion point → enter the cell BELOW the cursor
      // insertionIdx N = gap between row[N-1] and row[N].
      // "Below" = row[N], first cell in that row.
      if (insertionIdx < rows.length) {
        const targetRow = rows[insertionIdx];
        if (targetRow) {
          insertionIdx = null;
          const firstCell = targetRow.cells[0];
          if (firstCell) cellFocusFns[firstCell.id]?.();
        } else {
          insertionIdx = Math.min(rows.length, insertionIdx + 1);
        }
      } else {
        insertionIdx = null;
      }
      return;
    }
    if (e.key === 'Enter') {
      e.preventDefault();
      const idx = insertionIdx; insertionIdx = null;
      const id = nb.store.insertRowAt(idx);
      await tick(); cellFocusFns[id]?.();
      return;
    }
    // Printable char → create cell and type into it
    if (e.key.length === 1 && !e.ctrlKey && !e.metaKey) {
      e.preventDefault();
      const idx = insertionIdx; insertionIdx = null;
      const id = nb.store.insertRowAt(idx, 'code', e.key);
      await tick(); cellFocusFns[id]?.();
    }
  }

  /* ---------------------------------------------------------------------------
     Hand this card's controls to the toolbar.

     Registered for every card, not only the focused one: the registry is keyed
     by notebook id and the toolbar reads only the ACTIVE pane's entry, so a
     canvas card's presence costs nothing and means commands that target a
     notebook by id have somewhere to go.

     The stable methods go up ONCE, in onMount. The volatile flags go up in a
     reactive block. Publishing both together meant every keystroke replaced the
     whole record -- `allSectionsCollapsed` derives from the notebook's rows --
     and re-ran every toolbar subscriber.

     `paneToken` identifies this card instance. Switching between canvas and
     focused mode destroys one card and creates another for the SAME notebook id,
     and Svelte does not promise the teardown runs first; without the token, a
     late destroy would delete the entry the replacement card had just
     registered. */
  const paneToken = {};

  onMount(() => {
    registerPane(nb.id, paneToken, {
      notebookId: nb.id,
      store: nb.store,
      runAll,
      runCell: runCellById,
      runRange,
      focusCell: (cellId: string) => { cellFocusFns[cellId]?.(); },
      toggleLayout: () => { horizontal = !horizontal; },
      rename: startRename,
      toggleAllSections,
      toggleCollapse: onToggleCollapse,
      close: () => removeNotebook(nb.id),
    });
  });

  /* Every identifier this block reads must stay written out here. Svelte 4
     derives a reactive block's dependencies from the identifiers it
     syntactically mentions, so hoisting this object into a helper function would
     stop the block re-running and freeze the toolbar's icon states. */
  $: publishPaneFlags(nb.id, paneToken, {
    horizontal,
    allSectionsCollapsed,
    hasSections: sectionRowIds.length > 0,
    collapsed: !!nb.collapsed,
  });

  onDestroy(() => retractPane(nb.id, paneToken));

</script>

<!-- svelte-ignore a11y-no-static-element-interactions -->
<div
  class="nb-card"
  use:nameQueryListener
  class:mounted
  class:collapsed={nb.collapsed}
  class:focused-card={focused}
  class:selected={isSelected}
  style="--accent: {accentColor}; --accent-glow: {accentColor}1a;"
  bind:this={cardEl}
  tabindex="-1"
  on:keydown={onCardKeydown}
  on:pointermove={(e) => { onTitlePointerMove(e); onBottomResizeMove(e); }}
  on:pointerup={(e) => { onTitlePointerUp(e); onBottomResizeUp(e); }}
  on:pointercancel={(e) => { onTitlePointerUp(e); onBottomResizeUp(e); }}
>
  <!-- Title bar — only pointerdown here; move/up handled by cardEl after setPointerCapture -->
  <!-- svelte-ignore a11y-no-static-element-interactions -->
  <!-- In focused mode this bar is not rendered at all: its four controls and
       its title move up into the app bar, which would otherwise sit directly
       on top of it showing the same notebook name twice. -->
  {#if !focused}
  <div
    class="card-titlebar"
    on:pointerdown={onTitlePointerDown}
  >
    {#if renaming}
      <!-- In-place editable title — looks identical to the plain title -->
      <!-- svelte-ignore a11y-click-events-have-key-events -->
      <span
        class="card-title renaming"
        contenteditable="true"
        bind:this={renameInput}
        on:blur={commitRename}
        on:keydown={onRenameKeydown}
        on:pointerdown|stopPropagation
        on:click|stopPropagation
      >{renameValue}</span>
    {:else}
      <!-- svelte-ignore a11y-click-events-have-key-events -->
      <!-- pointer-events restored so dblclick rename works -->
      <span
        class="card-title"
        style="pointer-events:auto;"
        on:dblclick|stopPropagation={startRename}
        title="Double-click to rename"
      >{nb.title}</span>
    {/if}

    <div class="titlebar-actions">
      <button class="tb-btn tb-run-all" title="Run all cells" on:click|stopPropagation={runAll}>▶▶</button>
      <!-- Layout toggle always visible — input/output side-by-side or stacked -->
      <button
        class="tb-btn"
        title={horizontal ? 'Vertical layout' : 'Horizontal layout'}
        on:click|stopPropagation={() => { horizontal = !horizontal; }}
      >{horizontal ? '↕' : '⇄'}</button>
      {#if sectionRowIds.length > 0}
        <button
          class="tb-btn"
          title={allSectionsCollapsed ? 'Expand all sections' : 'Collapse all sections'}
          on:click|stopPropagation={toggleAllSections}
        >{allSectionsCollapsed ? '⌄' : '⌃'}</button>
      {/if}
      {#if !focused}
        <button class="tb-btn" title="Rename" on:click|stopPropagation={startRename}>✎</button>
        <button class="tb-btn tb-focus" title="Full screen" on:click|stopPropagation={() => setFocused(nb.id)}>⤢</button>
        <button class="tb-btn" title="Collapse / expand" on:click|stopPropagation={onToggleCollapse}>
          {nb.collapsed ? '⊟' : '⊞'}
        </button>
        <button class="tb-btn tb-close" title="Close" on:click|stopPropagation={() => removeNotebook(nb.id)}>✕</button>
      {/if}
    </div>
  </div>
  {/if}

  <!-- Right-edge resize handle -->
  {#if !focused}
    <!-- svelte-ignore a11y-no-static-element-interactions -->
    <div
      class="resize-handle"
      on:pointerdown={onResizePointerDown}
      on:pointermove={onResizePointerMove}
      on:pointerup={onResizePointerUp}
      on:pointercancel={onResizePointerUp}
    ></div>
    <!-- Bottom resize handle -->
    <!-- svelte-ignore a11y-no-static-element-interactions -->
    <div class="resize-handle-bottom" on:pointerdown={onBottomResizeDown}></div>
    <!-- Corner resize handle -->
    <!-- svelte-ignore a11y-no-static-element-interactions -->
    <div class="resize-handle-corner" on:pointerdown={onCornerResizeDown}></div>
  {/if}

  <!-- Collapse wrapper -->
  <div class="collapse-wrapper">
    {#if nb.collapsed}
      <!-- Collapsed summary -->
      <div class="collapsed-body">
        <span class="cell-count">{totalCells} cell{totalCells !== 1 ? 's' : ''}</span>
        {#if lastOutputPreview}
          <span class="last-preview">{lastOutputPreview}</span>
        {/if}
      </div>
    {:else}
      <!-- Full notebook UI -->
      <!-- svelte-ignore a11y-no-static-element-interactions -->
      <div
        class="card-body"
        style={nb.height != null ? `max-height: none; height: ${nb.height - TITLE_BAR_H}px; overflow-y: auto;` : ''}

      >
        <!-- Insertion point before first row -->
        {#if insertionIdx === 0}
          <div class="insertion-cursor active"></div>
        {/if}

        {#each $nbStore as row, ri (row.id)}
          {#if !hiddenRows.has(row.id)}
            {@const firstCell = row.cells[0]}
            {@const isSection = firstCell?.type === 'section' || firstCell?.type === 'subsection'}

            {#if isSection}
              <!-- Section/Subsection row with collapse toggle -->
              <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
              <div class="section-row" class:subsection={firstCell.type === 'subsection'}
                   style="--row-depth:{rowDepth.get(row.id) ?? 0}">
                <button
                  class="section-collapse-btn"
                  on:click|stopPropagation={() => toggleSection(row.id)}
                  title={collapsedSections.has(row.id) ? 'Expand' : 'Collapse'}
                >{collapsedSections.has(row.id) ? '▶' : '▼'}</button>
                <div class="cell-col" style="flex:1;min-width:0">
                  <CellShell
                    cell={firstCell}
                    store={nb.store}
                    rowId={row.id}
                    cellIdx={0}
                    notebookId={nb.id}
                    headingReadonly={isRefpage}
                    on:headingClick={(e) => toggleSection(e.detail.rowId)}
                    on:run={handleRun}
                    on:change={handleChange}
                    on:addAbove={addRowAbove}
                    on:addBelow={addRowBelow}
                    on:addLeft={addCellLeft}
                    on:addRight={addCellRight}
                    on:register={handleRegister}
                    on:focusPrev={handleFocusPrev}
                    on:focusNext={handleFocusNext}
                    {horizontal}
                  />
                </div>
              </div>
            {:else}
              <!-- Indent by section depth so a subsection's content reads as
                   belonging to it rather than sitting flush with the heading. -->
              <div class="nb-row" style="--row-depth:{rowDepth.get(row.id) ?? 0}">
                {#each row.cells as cell, ci (cell.id)}
                  <div class="cell-col" style="flex: 1 1 0; min-width: 0;">
                    <CellShell
                      {cell}
                      store={nb.store}
                      rowId={row.id}
                      cellIdx={ci}
                      notebookId={nb.id}
                      headingReadonly={isRefpage}
                      on:run={handleRun}
                      on:change={handleChange}
                      on:addAbove={addRowAbove}
                      on:addBelow={addRowBelow}
                      on:addLeft={addCellLeft}
                      on:addRight={addCellRight}
                      on:register={handleRegister}
                      on:focusPrev={handleFocusPrev}
                      on:focusNext={handleFocusNext}
                    {horizontal}
                    />
                  </div>
                  {#if ci < row.cells.length - 1}
                    <div class="col-divider"></div>
                  {/if}
                {/each}
              </div>
            {/if}

            <!-- Insertion point after this row. INSIDE the visibility guard:
                 rendered per row regardless, a collapsed section still emitted
                 one 3px strip for every row it was hiding, so folding
                 "Examples (26)" left ~90px of dead space behind it while a
                 section with no content left none. -->
            {#if insertionIdx === ri + 1}
              <div class="insertion-cursor active"></div>
            {:else if ri < $nbStore.length - 1}
              <div class="insertion-cursor"></div>
            {/if}
          {/if}
        {/each}

        <!-- Add-row bar -->
        <div class="add-row-bar">
          <button on:click|stopPropagation={() => addRow('code')}>＋ Code</button>
          <button on:click|stopPropagation={() => addRow('text')}>＋ Text</button>
          <button on:click|stopPropagation={() => addRow('section')}>＋ Section</button>
        </div>
      </div>
    {/if}
  </div>
</div>

<style>
  /* ---- Card shell ---- */
  .nb-card {
    position: relative;
    border-radius: 12px;
    background: var(--card-bg, rgba(12,15,28,0.85));
    border: 1px solid var(--card-border, rgba(255,255,255,0.08));
    box-shadow:
      0 0 0 1px rgba(137,180,250,0.10),
      0 24px 64px rgba(0,0,0,0.65),
      inset 0 1px 0 rgba(255,255,255,0.05);
    backdrop-filter: blur(20px);
    -webkit-backdrop-filter: blur(20px);
    overflow: hidden;
    /* Mount animation */
    opacity: 0;
    transform: scale(0.94);
    transition:
      opacity 180ms ease,
      transform 180ms ease;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    color: var(--text, #cdd6f4);
  }

  .nb-card.mounted {
    opacity: 1;
    transform: scale(1);
  }

  /* Multi-selected via rubber-band */
  .nb-card.selected {
    outline: 2px solid var(--accent, #89b4fa);
    outline-offset: 2px;
    box-shadow: 0 0 0 4px rgba(137,180,250,0.18), 0 8px 32px rgba(0,0,0,0.45);
  }

  /* Full-screen focused mode — edge to edge, no card chrome */
  .nb-card.focused-card {
    border-radius: 0;
    border: none;
    box-shadow: none;
    backdrop-filter: none;
    -webkit-backdrop-filter: none;
    background: var(--bg, #050810);
    /* No min-height here: Canvas.svelte's .focused-view is a flex column and
       grows this card to fill it. A 100vh here overflowed the inset view by the
       bar height; a percentage min-height could not resolve at all. See the note
       on .focused-view in Canvas.svelte. */
    width: 100%;
    display: flex;
    flex-direction: column;
  }

  /* Title centered absolutely in focused mode */
  .nb-card.focused-card .card-titlebar {
    position: sticky;
    top: 0;
    z-index: 10;
    background: rgba(5, 8, 16, 0.92);
    backdrop-filter: blur(12px);
    -webkit-backdrop-filter: blur(12px);
    cursor: default;
    /* Extend the bar up under the device status bar / notch so its buttons sit
       BELOW it (they otherwise collide with the clock/battery on mobile). The
       bar's dark background fills the inset strip for a native look; env() is 0
       on desktop, leaving the 36px bar unchanged. box-sizing is border-box, so
       height must be auto to grow by the padding rather than shrink content. */
    height: auto;
    padding-top: env(safe-area-inset-top, 0px);
  }
  /* Title centered in the 36px button row, below the status-bar inset */
  .nb-card.focused-card .card-title {
    position: absolute;
    left: 50%;
    top: calc(env(safe-area-inset-top, 0px) + 18px);
    transform: translate(-50%, -50%);
    pointer-events: none;
    max-width: 60%;
  }

  /* Card body fills remaining height in focused mode */
  .nb-card.focused-card .collapse-wrapper {
    flex: 1;
    display: flex;
    flex-direction: column;
  }
  .nb-card.focused-card .card-body {
    flex: 1;
  }

  /* ---- Title bar ---- */
  .card-titlebar {
    position: relative;    /* needed for absolute-centered title */
    display: flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0 10px 0 0;
    height: 36px;
    background: var(--gutter-bg, rgba(255,255,255,0.03));
    border-bottom: 1px solid var(--border, rgba(255,255,255,0.06));
    cursor: grab;
    user-select: none;
    -webkit-user-select: none;
  }

  .card-titlebar:active { cursor: grabbing; }

  .titlebar-accent {
    width: 3px;
    align-self: stretch;
    border-radius: 12px 0 0 0;
    flex-shrink: 0;
  }

  /* Title absolutely centered in the titlebar — bold, ignores button widths */
  .card-title {
    position: absolute;
    left: 50%;
    transform: translateX(-50%);
    font-size: 0.84rem;
    font-weight: 700;
    color: var(--text, #cdd6f4);
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    max-width: 60%;
    pointer-events: none;   /* let drag pass through to the titlebar */
    cursor: grab;
  }

  /* Renaming state: look identical to the title, just editable */
  .card-title.renaming {
    outline: none;
    border-bottom: 1.5px solid var(--accent, #89b4fa);
    cursor: text;
    min-width: 40px;
  }

  .titlebar-actions {
    display: flex;
    align-items: center;
    gap: 4px;
    flex-shrink: 0;
  }

  .tb-btn {
    background: none;
    border: none;
    color: #585b70;
    cursor: pointer;
    font-size: 0.9rem;
    padding: 2px 5px;
    border-radius: 4px;
    line-height: 1;
    transition: color 0.12s, background 0.12s;
  }
  .tb-btn:hover { color: var(--text, #cdd6f4); background: rgba(128,128,128,0.12); }
  .tb-run-all {
    font-size: 0.72rem;
    padding: 2px 8px;
    background: rgba(166,227,161,0.12);
    border: 1px solid rgba(166,227,161,0.3);
    color: #a6e3a1;
    border-radius: 5px;
  }
  .tb-run-all:hover { background: rgba(166,227,161,0.22) !important; }
  .tb-close:hover { color: #f38ba8; }
  .tb-focus:hover { color: var(--accent, #89b4fa); }
  /* ---- Collapse wrapper ---- */
  .collapse-wrapper {
    overflow: hidden;
  }

  /* ---- Collapsed body ---- */
  .collapsed-body {
    padding: 8px 12px;
    display: flex;
    align-items: center;
    gap: 0.75rem;
  }

  .cell-count {
    font-size: 0.75rem;
    color: #585b70;
    font-variant-numeric: tabular-nums;
  }

  .last-preview {
    font-size: 0.75rem;
    color: #585b70;
    font-family: 'SF Mono', 'Fira Code', monospace;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    max-width: 340px;
  }

  /* ---- Expanded card body ---- */
  .card-body {
    max-height: 70vh;
    overflow-y: auto;
    scrollbar-width: thin;
    scrollbar-color: rgba(137,180,250,0.2) transparent;
    /* The canvas sets touch-action:none to own pan/zoom; re-enable vertical
       finger-scroll inside the card body. */
    touch-action: pan-y;
  }

  /* In focused (full-screen) mode, remove max-height — let .focused-view scroll */
  .focused-card .card-body {
    max-height: none;
    overflow-y: visible;
  }

  /* Also remove overflow:hidden clip from the focused card so nothing gets cut */
  .focused-card {
    overflow: visible !important;
  }

  .card-body::-webkit-scrollbar { width: 5px; }
  .card-body::-webkit-scrollbar-track { background: transparent; }
  .card-body::-webkit-scrollbar-thumb { background: rgba(137,180,250,0.2); border-radius: 3px; }

  /* ---- Right-edge resize handle ---- */
  .resize-handle {
    position: absolute;
    top: 0;
    right: -4px;
    width: 8px;
    height: 100%;
    cursor: ew-resize;
    z-index: 10;
  }
  .resize-handle:hover { background: rgba(137,180,250,0.25); border-radius: 4px; }

  /* ---- Bottom resize handle ---- */
  .resize-handle-bottom {
    position: absolute;
    bottom: -4px;
    left: 0;
    right: 8px;
    height: 8px;
    cursor: ns-resize;
    z-index: 10;
  }
  .resize-handle-bottom:hover { background: rgba(137,180,250,0.25); border-radius: 4px; }

  /* ---- Corner resize handle ---- */
  .resize-handle-corner {
    position: absolute;
    bottom: -4px;
    right: -4px;
    width: 12px;
    height: 12px;
    cursor: nwse-resize;
    z-index: 11;
  }
  .resize-handle-corner:hover { background: rgba(137,180,250,0.4); border-radius: 2px; }

  /* ---- Insertion cursor between rows ---- */
  .insertion-cursor {
    height: 3px;
    margin: 0;
    border-radius: 2px;
    transition: background 0.1s;
    cursor: text;
  }
  .insertion-cursor.active {
    background: var(--accent, #89b4fa);
    animation: ins-blink 1s ease-in-out infinite;
  }
  @keyframes ins-blink { 0%,100%{opacity:1} 50%{opacity:0.35} }

  /* ---- Section rows with collapse toggle ---- */
  .section-row {
    display: flex;
    align-items: flex-start;
    border-bottom: 1px solid rgba(255,255,255,0.05);
  }
  /* The subsection's own fold arrow is smaller and dimmer, matching its
     heading rather than the section above it. */
  .section-row.subsection .section-collapse-btn {
    font-size: 0.62rem;
    opacity: 0.65;
  }

  .section-collapse-btn {
    flex-shrink: 0;
    width: 22px;
    padding: 8px 4px 0;
    background: none;
    border: none;
    color: var(--text-muted);
    font-size: 0.6rem;
    cursor: pointer;
    transition: color 0.1s;
    line-height: 1;
  }
  .section-collapse-btn:hover { color: var(--accent); }

  /* ---- Notebook rows ---- */
  /* One step of indentation per level of section nesting. Applied as padding so
     the row's own background and dividers still span the full card width. */
  .nb-row, .section-row { padding-left: calc(var(--row-depth, 0) * 1.15rem); }

  .nb-row {
    display: flex;
    flex-direction: row;
    align-items: stretch;
    border-bottom: 1px solid rgba(255,255,255,0.05);
  }

  .cell-col { min-width: 0; }

  .col-divider {
    width: 1px;
    background: rgba(255,255,255,0.05);
    flex-shrink: 0;
  }

  /* ---- Add-row bar ---- */
  .add-row-bar {
    padding: 0.5rem 0.5rem;
    display: flex;
    gap: 0.35rem;
  }

  .add-row-bar button {
    background: none;
    border: 1px dashed rgba(255,255,255,0.12);
    color: #585b70;
    border-radius: 4px;
    padding: 0.2rem 0.65rem;
    font-size: 0.75rem;
    cursor: pointer;
    transition: all 0.12s;
  }

  .add-row-bar button:hover {
    border-color: var(--accent, #89b4fa);
    color: var(--accent, #89b4fa);
  }
</style>
