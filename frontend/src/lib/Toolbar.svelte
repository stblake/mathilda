<!--
  Toolbar.svelte — the focused-mode toolbar.

  Replaces the app bar's contents while a notebook fills the window: labelled,
  vertically-ruled groups instead of a row of unlabelled glyphs.

  Two rules the whole file depends on:

  1. Every button carries on:pointerdown|preventDefault and tabindex="-1".
     Suppressing pointerdown's default stops the focus transfer, so the editor
     never blurs and its live text selection survives. That is what makes a
     toolbar work on a text editor at all -- without it, every command would run
     against a cell that had just lost the caret. It is the primary mechanism;
     active.ts's sticky record is the backstop for controls that must take focus.

  2. Exactly one <Menu> is rendered, at the end, driven by `openMenu`. One
     backdrop, one keyboard handler, one measure path -- rather than one of each
     per dropdown. It must live here and not inside a card: .nb-card has
     backdrop-filter, which creates a containing block that would trap a
     position:fixed popover inside the card's bounds.
-->
<script lang="ts">
  import { onDestroy } from 'svelte';
  import Icon from './Icon.svelte';
  import ToolbarGroup from './ToolbarGroup.svelte';
  import Menu from './Menu.svelte';
  import type { MenuItem } from './Menu.svelte';
  import { canvasState, setFocused, activeActions, activeFlags, openRefpage } from './canvas';
  import { activeCell, retypeActiveCell } from './active';
  import { darkMode } from './theme';
  import { symbolAtSelection } from './refpages';
  import type { Cell, CellType, NotebookRow } from './notebook';

  type MenuId = 'style' | 'overflow' | null;
  let openMenu: MenuId = null;

  /* One anchor per trigger, so the menu can measure against the button that
     opened it. */
  let styleAnchor: HTMLElement;
  let overflowAnchor: HTMLElement;

  function toggleMenu(which: Exclude<MenuId, null>) {
    openMenu = openMenu === which ? null : which;
  }

  // ---------------------------------------------------------------------------
  // The active cell, re-resolved from its notebook's store on every read.
  //
  // This is what makes a stale active-cell record harmless. Nothing cleans up
  // after a deleted cell, a closed notebook, or a library load that replaced
  // every store -- instead the id is looked up fresh, and a miss simply reads as
  // "no active cell". No teardown hooks to forget.
  //
  // The pane's store is subscribed manually because it changes with the active
  // pane, and Svelte's `$store` auto-subscription only works on a fixed
  // identifier. Without this subscription the toolbar would not notice rows
  // being added or removed.

  let rows: NotebookRow[] = [];
  let unsubRows: (() => void) | null = null;

  $: {
    unsubRows?.();
    unsubRows = null;
    const store = $activeActions?.store;
    if (store) unsubRows = store.subscribe((r: NotebookRow[]) => { rows = r; });
  }
  onDestroy(() => unsubRows?.());

  $: activeCellObj = ($activeCell && $activeActions && $activeCell.notebookId === $activeActions.notebookId)
    ? (rows.flatMap(r => r.cells).find((c: Cell) => c.id === $activeCell!.cellId) ?? null)
    : null;

  // ---------------------------------------------------------------------------
  // Cell Style

  const STYLES: { id: CellType; label: string }[] = [
    { id: 'code',       label: 'Code' },
    { id: 'text',       label: 'Text' },
    { id: 'section',    label: 'Section' },
    { id: 'subsection', label: 'Subsection' },
  ];

  /* A reference page's cells are generated documentation structure, not the
     reader's to retype. Shown, so the control still says what the cell IS. */
  $: styleLabel = activeCellObj
    ? (activeCellObj.type === 'ref'
        ? 'Reference'
        : STYLES.find(s => s.id === activeCellObj!.type)?.label ?? 'Cell')
    : 'Insert Cell…';

  $: styleLocked = activeCellObj?.type === 'ref';

  $: styleItems = STYLES.map(s => ({
    kind: 'item' as const,
    id: s.id,
    label: activeCellObj ? s.label : `${s.label} cell`,
    checked: activeCellObj?.type === s.id,
    disabled: styleLocked,
  })) as MenuItem[];

  function onStyleSelect(id: string) {
    const type = id as CellType;
    const pane = $activeActions;
    if (!pane) return;
    if (activeCellObj) {
      if (activeCellObj.type === 'ref') return;
      pane.store.setCellType(activeCellObj.id, type);
      retypeActiveCell(type);
    } else {
      /* Nothing is active, so the control reads "Insert Cell…" and this INSERTS
         rather than retypes -- which is what the label promises. */
      const newId = pane.store.addRow(type);
      pane.focusCell(newId);
    }
  }

  // ---------------------------------------------------------------------------
  // Docs
  //
  // Whether a documented symbol is under the caret is not reactive state -- it
  // depends on the DOM selection. Recomputed on selectionchange (and whenever
  // the active cell changes) so the control can honestly show a disabled state
  // instead of looking live and doing nothing.

  let docsTarget: string | null = null;

  function refreshDocsTarget() {
    try { docsTarget = symbolAtSelection()?.name ?? null; }
    catch { docsTarget = null; }
  }
  $: { void $activeCell; refreshDocsTarget(); }

  function openDocs() {
    const hit = symbolAtSelection();
    const pane = $activeActions;
    if (!hit || !pane) return;
    openRefpage(pane.notebookId, hit.name, hit.at);
  }

  // ---------------------------------------------------------------------------
  // Overflow — the controls with no Wolfram group of their own.

  $: overflowItems = [
    { kind: 'item', id: 'layout', label: $activeFlags?.horizontal ? 'Stack input and output' : 'Input and output side by side', icon: $activeFlags?.horizontal ? 'layoutV' : 'layoutH' },
    ...($activeFlags?.hasSections
      ? [{ kind: 'item' as const, id: 'sections', label: $activeFlags?.allSectionsCollapsed ? 'Expand all sections' : 'Collapse all sections', icon: $activeFlags?.allSectionsCollapsed ? 'caret' : 'caretUp' }]
      : []),
    { kind: 'sep' },
    { kind: 'item', id: 'rename', label: 'Rename notebook…', icon: 'rename' },
    { kind: 'item', id: 'collapse', label: $activeFlags?.collapsed ? 'Expand notebook' : 'Collapse notebook', icon: 'collapse' },
    { kind: 'sep' },
    { kind: 'item', id: 'close', label: 'Close notebook', icon: 'trash' },
  ] as MenuItem[];

  function onOverflowSelect(id: string) {
    const pane = $activeActions;
    if (!pane) return;
    switch (id) {
      case 'layout':   pane.toggleLayout(); break;
      case 'sections': pane.toggleAllSections(); break;
      case 'rename':   pane.rename(); break;
      case 'collapse': pane.toggleCollapse(); break;
      case 'close':    pane.close(); break;
    }
  }

  // ---------------------------------------------------------------------------

  $: activeTitle = $canvasState.focusedActiveId
    ? ($canvasState.notebooks.find(n => n.id === $canvasState.focusedActiveId)?.title ?? '')
    : '';

  $: menuItems = openMenu === 'style' ? styleItems : openMenu === 'overflow' ? overflowItems : [];
  $: menuAnchor = openMenu === 'style' ? styleAnchor : openMenu === 'overflow' ? overflowAnchor : null;

  function onMenuSelect(e: CustomEvent<{ id: string }>) {
    if (openMenu === 'style') onStyleSelect(e.detail.id);
    else if (openMenu === 'overflow') onOverflowSelect(e.detail.id);
  }
</script>

<svelte:document on:selectionchange={refreshDocsTarget} />

<!-- Leaving focused mode. NOT in a group and never in the overflow: this is the
     only way back to the canvas that works with a plain mouse -- the alternative
     is a pinch-out gesture, which needs a trackpad or a touchscreen. Putting it
     anywhere collapsible would make focused mode a trap. -->
<button
  class="tb-btn tb-back"
  title="Back to canvas (or pinch out)"
  tabindex="-1"
  on:pointerdown|preventDefault
  on:click={() => setFocused(null)}
><Icon name="back" size={17} /></button>

<ToolbarGroup label="Cell Style">
  <!-- svelte-ignore a11y-no-static-element-interactions -->
  <button
    class="tb-combo"
    class:muted={!activeCellObj}
    title={activeCellObj ? 'Change the active cell’s type' : 'Insert a new cell'}
    aria-haspopup="menu"
    aria-expanded={openMenu === 'style'}
    tabindex="-1"
    bind:this={styleAnchor}
    on:pointerdown|preventDefault
    on:click={() => toggleMenu('style')}
  >
    <span class="tb-combo-label">{styleLabel}</span>
    <Icon name="caret" size={13} />
  </button>
</ToolbarGroup>

<ToolbarGroup label="Notebook">
  <button
    class="tb-btn"
    title={$darkMode ? 'Switch to light theme' : 'Switch to dark theme'}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={() => darkMode.update(v => !v)}
  >{$darkMode ? '◑' : '☀'}</button>

  <button
    class="tb-btn"
    title={docsTarget ? `Documentation for ${docsTarget}` : 'Documentation — put the caret on a symbol first'}
    disabled={!docsTarget}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={openDocs}
  ><Icon name="docs" /></button>
</ToolbarGroup>

<span class="tb-spacer"></span>

<!-- Names what the toolbar's buttons will act on. With one pane that is the
     notebook filling the window; with several it is the active one. -->
{#if activeTitle}<span class="tb-title" title={activeTitle}>{activeTitle}</span>{/if}

<button
  class="tb-btn tb-overflow"
  title="More notebook actions"
  aria-haspopup="menu"
  aria-expanded={openMenu === 'overflow'}
  tabindex="-1"
  bind:this={overflowAnchor}
  on:pointerdown|preventDefault
  on:click={() => toggleMenu('overflow')}
>…</button>

<Menu
  open={openMenu !== null}
  items={menuItems}
  anchor={menuAnchor}
  align={openMenu === 'overflow' ? 'end' : 'start'}
  minWidth={openMenu === 'style' ? 170 : 0}
  on:select={onMenuSelect}
  on:close={() => (openMenu = null)}
/>

<style>
  /* ---- Buttons ---- */
  .tb-btn {
    display: flex;
    align-items: center;
    justify-content: center;
    min-width: var(--tb-btn-sz, 24px);
    height: var(--tb-btn-sz, 24px);
    padding: 0 4px;
    background: none;
    border: 1px solid transparent;
    border-radius: 5px;
    color: var(--text);
    font-size: 0.82rem;
    line-height: 1;
    cursor: pointer;
    transition: background 0.1s, color 0.1s;
  }
  .tb-btn:hover:not(:disabled) { background: var(--surface-2); color: var(--text-h); }
  .tb-btn:active:not(:disabled) { background: var(--surface-3); }
  .tb-btn:disabled { opacity: 0.32; cursor: default; }

  /* The exit control sits outside every group and gets its own vertical rule,
     so it reads as chrome rather than as one of the notebook's commands. */
  .tb-back {
    align-self: center;
    margin-right: 4px;
    padding: 0 6px;
    border-inline-end: 1px solid var(--tb-rule);
    border-radius: 5px 0 0 5px;
    color: var(--accent);
  }
  .tb-back:hover { background: var(--surface-2); }

  .tb-overflow {
    align-self: center;
    font-size: 1rem;
    letter-spacing: 0.06em;
    margin-left: 2px;
  }

  /* ---- Cell-style combo ---- */
  .tb-combo {
    display: flex;
    align-items: center;
    gap: 4px;
    height: var(--tb-btn-sz, 24px);
    min-width: 122px;
    padding: 0 5px 0 8px;
    background: var(--surface-2);
    border: 1px solid var(--tb-rule);
    border-radius: 5px;
    color: var(--text-h);
    font: 500 0.76rem/1 var(--sans);
    cursor: pointer;
    text-align: left;
  }
  .tb-combo:hover { border-color: var(--accent); }
  /* No active cell: the control still works (it inserts) but should not claim
     to be describing something. */
  .tb-combo.muted { color: var(--tb-caption); }
  .tb-combo-label { flex: 1; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }

  /* ---- Right side ---- */
  .tb-spacer { flex: 1; min-width: 8px; }

  .tb-title {
    align-self: center;
    font: 500 0.78rem/1 var(--sans);
    color: var(--text);
    max-width: 22vw;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
    padding-right: 4px;
    user-select: none;
    -webkit-user-select: none;
  }
</style>
