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
  import { canvasState, setFocused, activeActions, activeFlags, openRefpage, addQueryNotebook,
           addPaneToFocus, removePane, setFocusLayout, splitWithNext, MAX_PANES } from './canvas';
  import type { FocusLayout } from './canvas';
  import { activeCell, retypeActiveCell, activeHandle } from './active';
  import { splitCell, mergeCellDown, duplicateCell, deleteCell,
           indentCode, outdentCode, commentCode, duplicateLine } from './cellCommands';
  import { darkMode } from './theme';
  import { symbolAtSelection } from './refpages';
  import { kernelStatus } from './notebook';
  import { restart, abortEvaluation } from './kernelActions';
  import { showStatusBar, resetSessionStats } from './status';
  import { propertiesOpen } from './properties';
  import type { Cell, CellType, NotebookRow } from './notebook';

  type MenuId = 'eval' | 'kernel' | 'docs' | 'style' | 'addpane' | 'overflow' | null;
  let openMenu: MenuId = null;

  /* One anchor per trigger, so the menu can measure against the button that
     opened it. */
  let evalAnchor: HTMLElement;
  let kernelAnchor: HTMLElement;
  let docsAnchor: HTMLElement;
  let styleAnchor: HTMLElement;
  let addPaneAnchor: HTMLElement;
  let overflowAnchor: HTMLElement;

  function toggleMenu(which: Exclude<MenuId, null>) {
    /* Read the caret's symbol at open time. Must happen BEFORE the menu takes
       focus, while the editor's selection is still the document selection. */
    if (which === 'docs' && openMenu !== 'docs') refreshDocsTarget();
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

  $: cells = rows.flatMap(r => r.cells) as Cell[];

  $: activeCellObj = ($activeCell && $activeActions && $activeCell.notebookId === $activeActions.notebookId)
    ? (cells.find((c: Cell) => c.id === $activeCell!.cellId) ?? null)
    : null;

  /* Position in the flattened cell order, which is what runRange indexes. */
  $: activeIdx = activeCellObj ? cells.findIndex(c => c.id === activeCellObj!.id) : -1;

  // ---------------------------------------------------------------------------
  // Evaluation

  $: canRunCell = activeCellObj?.type === 'code' && activeCellObj.source.trim().length > 0;

  /* The bare Run button does the obvious thing for where the caret is: run the
     cell you are in, or the whole notebook if you are not in one. */
  function runPrimary() {
    const pane = $activeActions;
    if (!pane) return;
    if (canRunCell && activeCellObj) pane.runCell(activeCellObj.id);
    else pane.runAll();
  }

  $: evalItems = [
    { kind: 'item', id: 'cell', label: 'Evaluate cell', icon: 'run', hint: '⇧↵', disabled: !canRunCell },
    { kind: 'item', id: 'notebook', label: 'Evaluate notebook', icon: 'runAll' },
    { kind: 'sep' },
    { kind: 'item', id: 'from', label: 'Evaluate from here down', disabled: activeIdx < 0 },
    { kind: 'item', id: 'above', label: 'Evaluate above', disabled: activeIdx <= 0 },
  ] as MenuItem[];

  function onEvalSelect(id: string) {
    const pane = $activeActions;
    if (!pane) return;
    switch (id) {
      case 'cell':     if (activeCellObj) pane.runCell(activeCellObj.id); break;
      case 'notebook': pane.runAll(); break;
      case 'from':     pane.runRange(activeIdx, cells.length - 1); break;
      case 'above':    pane.runRange(0, activeIdx - 1); break;
    }
  }

  const KERNEL_SHORT: Record<string, string> = {
    starting: 'Starting', ready: 'Ready', busy: 'Running',
    restarting: 'Restarting', dead: 'Not running',
  };

  /* "Abort" is spelled out because it is not free: interrupt_kernel kills the
     process without respawning it, so regaining control costs the session's
     definitions. A control whose consequence surprises you is worse than one
     that reads long. */
  $: kernelItems = [
    { kind: 'item', id: 'abort', label: 'Abort evaluation (restarts kernel)', icon: 'abort', hint: '⌘.', disabled: $kernelStatus !== 'busy' },
    { kind: 'item', id: 'restart', label: 'Restart kernel', icon: 'restart' },
  ] as MenuItem[];

  function onKernelSelect(id: string) {
    if (id === 'abort') { abortEvaluation(); resetSessionStats(); }
    else if (id === 'restart') { restart(); resetSessionStats(); }
  }

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
  // Layout — how many notebooks are on screen, and how they are arranged.

  $: paneCount = $canvasState.focusedIds.length;
  $: paneLayout = $canvasState.focusedLayout;
  $: atPaneCap = paneCount >= MAX_PANES;

  /* Notebooks not already on screen. A notebook may only appear in one pane:
     two cards over one store would share cell selection and open two reference
     pages from a single Cmd+click, since those stores are module-global. */
  $: addablePanes = $canvasState.notebooks.filter(n => !$canvasState.focusedIds.includes(n.id));

  $: addPaneItems = (addablePanes.length === 0
    ? [{ kind: 'item' as const, id: '', label: 'Every notebook is already open', disabled: true }]
    : addablePanes.map(n => ({ kind: 'item' as const, id: n.id, label: n.title }))) as MenuItem[];

  /** Add a notebook beside the current one, arranged the given way. */
  function addPane(id: string, layout?: FocusLayout) {
    if (layout && paneCount >= 1) setFocusLayout(layout);
    addPaneToFocus(id);
  }

  /* "Split" is the verb people reach for, so the two arrangement buttons double
     as the way IN to a split: with one notebook open they bring the next one in
     rather than doing nothing visible.
     One click, no picker. Which notebook arrives is deterministic (canvas order)
     rather than chosen, and the pane header names it -- so a wrong guess is
     immediately visible and one click from fixed with + or x. Making the user
     answer a menu before seeing any split at all was the friction. */
  function chooseLayout(layout: FocusLayout) {
    if (paneCount === 1 && addablePanes.length > 0) splitWithNext(layout);
    else setFocusLayout(layout);
  }

  function closeActivePane() {
    const id = $canvasState.focusedActiveId;
    if (id) removePane(id);
  }

  // ---------------------------------------------------------------------------
  // Cells — structural edits on the active cell.

  /* Merge needs a row below to merge with; the others only need a cell. */
  $: canMerge = activeIdx >= 0 && rows.length > 1
                && rows.findIndex(r => r.cells.some(c => c.id === activeCellObj?.id)) < rows.length - 1;

  function doSplit() {
    const pane = $activeActions;
    if (!pane || !activeCellObj) return;
    /* The caret offset only exists while the editor holds focus. Toolbar buttons
       suppress pointerdown's default precisely so it still does -- but if
       something took focus anyway, splitting at the end beats refusing. */
    const h = activeHandle();
    const at = $activeCell?.focused
      ? (h?.view ? h.view.state.selection.main.head
                 : (h?.el ? readContentEditableCaret(h.el) : null))
      : null;
    const newId = splitCell(pane.store, activeCellObj.id, at);
    if (newId) pane.focusCell(newId);
  }

  /** Caret offset inside a contenteditable, as a plain character index. */
  function readContentEditableCaret(el: HTMLElement): number | null {
    const sel = window.getSelection();
    if (!sel || !sel.rangeCount || !el.contains(sel.anchorNode)) return null;
    const r = sel.getRangeAt(0).cloneRange();
    r.selectNodeContents(el);
    r.setEnd(sel.getRangeAt(0).endContainer, sel.getRangeAt(0).endOffset);
    return r.toString().length;
  }

  function doMerge() {
    const pane = $activeActions;
    if (pane && activeCellObj) mergeCellDown(pane.store, activeCellObj.id);
  }

  function doDuplicate() {
    const pane = $activeActions;
    if (!pane || !activeCellObj) return;
    const newId = duplicateCell(pane.store, activeCellObj.id);
    if (newId) pane.focusCell(newId);
  }

  function doDelete() {
    const pane = $activeActions;
    if (pane && activeCellObj) deleteCell(pane.store, activeCellObj.id);
  }

  // ---------------------------------------------------------------------------
  // Code — the context-sensitive group, when the active cell is code.
  //
  // The Text half of this group needs text cells to render Markdown before B/I/U
  // mean anything: until then a bold button would insert literal ** and leave the
  // asterisks on screen. So it is not rendered yet rather than shipped inert, and
  // the group is hidden -- not disabled -- when no code cell is active, because a
  // row of greyed-out controls under a caption reads as broken rather than as
  // inapplicable.

  $: codeView = activeCellObj?.type === 'code' ? (activeHandle()?.view ?? null) : null;
  $: showCodeGroup = activeCellObj?.type === 'code';

  function withView(fn: (v: import('@codemirror/view').EditorView) => void) {
    const v = activeHandle()?.view;
    if (v) fn(v);
  }

  // ---------------------------------------------------------------------------
  // Docs
  //
  // Whether a documented symbol is under the caret is not reactive state -- it
  // depends on the DOM selection. Recomputed on selectionchange (and whenever
  // the active cell changes) so the control can honestly show a disabled state
  // instead of looking live and doing nothing.

  let docsTarget: string | null = null;

  /* Computed when the docs menu opens, NOT on every selection change.
     It used to run on document selectionchange, which fires on every caret move
     and every character typed -- and symbolAtSelection does a
     getBoundingClientRect plus up to four caretRangeFromPoint hit-tests, each
     forcing synchronous layout. Five forced layouts per keystroke to keep a
     tooltip and one menu label up to date is not a trade worth making, and the
     answer is only ever read at the moment the menu is opened. */
  function refreshDocsTarget() {
    try { docsTarget = symbolAtSelection()?.name ?? null; }
    catch { docsTarget = null; }
  }

  /* A menu rather than a bare button.
   *
   * As a button this was disabled whenever no documented symbol sat under the
   * caret -- which is most of the time -- and a control that is usually greyed
   * out and never says why reads as broken rather than as conditional. As a menu
   * it always has at least one thing it can actually do. */
  $: docsItems = [
    {
      kind: 'item',
      id: 'symbol',
      label: docsTarget ? `Documentation for ${docsTarget}` : 'Documentation for the symbol at the caret',
      icon: 'docs',
      hint: 'F1',
      disabled: !docsTarget,
    },
    { kind: 'sep' },
    { kind: 'item', id: 'browse', label: 'Browse all symbols…', icon: 'search' },
  ] as MenuItem[];

  function onDocsSelect(id: string) {
    const pane = $activeActions;
    if (!pane) return;
    if (id === 'symbol') {
      const hit = symbolAtSelection();
      if (hit) openRefpage(pane.notebookId, hit.name, hit.at);
    } else if (id === 'browse') {
      /* `?*` is the kernel's own wildcard symbol lookup, and its result already
         renders as a grid of names where clicking one opens that symbol's
         reference page. So this reuses a working path rather than inventing a
         browser. */
      addQueryNotebook(pane.notebookId, '?*', 'All Symbols');
    }
  }

  // ---------------------------------------------------------------------------
  // Overflow — the controls with no Wolfram group of their own.

  $: overflowItems = [
    { kind: 'item', id: 'layout', label: $activeFlags?.horizontal ? 'Stack input and output' : 'Input and output side by side', icon: $activeFlags?.horizontal ? 'layoutV' : 'layoutH' },
    ...($activeFlags?.hasSections
      ? [{ kind: 'item' as const, id: 'sections', label: $activeFlags?.allSectionsCollapsed ? 'Expand all sections' : 'Collapse all sections', icon: $activeFlags?.allSectionsCollapsed ? 'caret' : 'caretUp' }]
      : []),
    { kind: 'sep' },
    { kind: 'item', id: 'statusbar', label: 'Status bar', checked: $showStatusBar },
    { kind: 'sep' },
    { kind: 'item', id: 'rename', label: 'Rename notebook…', icon: 'rename' },
    { kind: 'item', id: 'collapse', label: $activeFlags?.collapsed ? 'Expand notebook' : 'Collapse notebook', icon: 'collapse' },
    { kind: 'sep' },
    { kind: 'item', id: 'close', label: 'Close notebook', icon: 'trash' },
  ] as MenuItem[];

  function onOverflowSelect(id: string) {
    /* The status bar is a view preference, not a notebook command, so it does
       not need a pane and must be handled before the guard below. */
    if (id === 'statusbar') { showStatusBar.update(v => !v); return; }
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

  /* One <Menu> instance, driven by these three tables. Adding a dropdown means
     adding a row to each, not another component with its own backdrop and
     keyboard handling. */
  $: menuItems =
    openMenu === 'eval'     ? evalItems :
    openMenu === 'kernel'   ? kernelItems :
    openMenu === 'docs'     ? docsItems :
    openMenu === 'style'    ? styleItems :
    openMenu === 'addpane'  ? addPaneItems :
    openMenu === 'overflow' ? overflowItems : [];

  $: menuAnchor =
    openMenu === 'eval'     ? evalAnchor :
    openMenu === 'kernel'   ? kernelAnchor :
    openMenu === 'docs'     ? docsAnchor :
    openMenu === 'style'    ? styleAnchor :
    openMenu === 'addpane'  ? addPaneAnchor :
    openMenu === 'overflow' ? overflowAnchor : null;

  function onMenuSelect(e: CustomEvent<{ id: string }>) {
    switch (openMenu) {
      case 'eval':     onEvalSelect(e.detail.id); break;
      case 'kernel':   onKernelSelect(e.detail.id); break;
      case 'docs':     onDocsSelect(e.detail.id); break;
      case 'style':    onStyleSelect(e.detail.id); break;
      case 'addpane':  if (e.detail.id) addPane(e.detail.id); break;
      case 'overflow': onOverflowSelect(e.detail.id); break;
    }
  }
</script>

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

<ToolbarGroup label="Evaluation">
  <button
    class="tb-btn tb-run"
    title={canRunCell ? 'Evaluate the active cell (Shift+Enter)' : 'Evaluate the whole notebook'}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={runPrimary}
  ><Icon name={canRunCell ? 'run' : 'runAll'} /></button>

  <button
    class="tb-caret"
    title="Evaluation options"
    aria-haspopup="menu"
    aria-expanded={openMenu === 'eval'}
    tabindex="-1"
    bind:this={evalAnchor}
    on:pointerdown|preventDefault
    on:click={() => toggleMenu('eval')}
  ><Icon name="caret" size={12} /></button>

  <span class="tb-mini-rule"></span>

  <!-- Kernel state, and the two things you can do about it. The dot carries the
       status so the text can stay short enough for a toolbar. -->
  <button
    class="tb-kernel"
    data-status={$kernelStatus}
    title={`Local kernel — ${KERNEL_SHORT[$kernelStatus] ?? $kernelStatus}`}
    aria-haspopup="menu"
    aria-expanded={openMenu === 'kernel'}
    tabindex="-1"
    bind:this={kernelAnchor}
    on:pointerdown|preventDefault
    on:click={() => toggleMenu('kernel')}
  >
    <span class="dot"></span>
    <span class="tb-kernel-label">{KERNEL_SHORT[$kernelStatus] ?? $kernelStatus}</span>
    <Icon name="caret" size={12} />
  </button>
</ToolbarGroup>

<!-- How many notebooks share the window, and how they sit. With one notebook the
     two arrangement buttons are also the way IN to a split: they set the
     arrangement and then ask which notebook to bring in, rather than being inert
     until a second pane exists. -->
<ToolbarGroup label="Layout">
  <button
    class="tb-btn"
    class:on={paneCount > 1 && paneLayout === 'h'}
    title={paneCount > 1 ? 'Arrange panes side by side' : 'Open a second notebook side by side'}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={() => chooseLayout('h')}
  ><Icon name="layoutH" /></button>

  <button
    class="tb-btn"
    class:on={paneCount > 1 && paneLayout === 'v'}
    title={paneCount > 1 ? 'Stack panes one above the other' : 'Open a second notebook below'}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={() => chooseLayout('v')}
  ><Icon name="layoutV" /></button>

  <button
    class="tb-btn"
    class:on={paneLayout === 'grid'}
    title={paneCount >= 3 ? 'Arrange panes in a 2×2 grid' : 'A 2×2 grid needs at least three notebooks open'}
    disabled={paneCount < 3}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={() => setFocusLayout('grid')}
  ><Icon name="layoutGrid" /></button>

  <span class="tb-mini-rule"></span>

  <button
    class="tb-btn"
    title={atPaneCap ? `At most ${MAX_PANES} notebooks can share the window` : 'Add a notebook to the window'}
    disabled={atPaneCap || addablePanes.length === 0}
    aria-haspopup="menu"
    aria-expanded={openMenu === 'addpane'}
    tabindex="-1"
    bind:this={addPaneAnchor}
    on:pointerdown|preventDefault
    on:click={() => toggleMenu('addpane')}
  ><Icon name="plus" /></button>

  <button
    class="tb-btn"
    title="Remove the active pane (its notebook stays on the canvas)"
    disabled={paneCount < 2}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={closeActivePane}
  ><Icon name="close" size={14} /></button>
</ToolbarGroup>

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

<ToolbarGroup label="Cells">
  <button class="tb-btn" title="Split the cell at the caret" disabled={!activeCellObj}
          tabindex="-1" on:pointerdown|preventDefault on:click={doSplit}
  ><Icon name="split" /></button>

  <button class="tb-btn" title="Merge this cell with the one below" disabled={!canMerge}
          tabindex="-1" on:pointerdown|preventDefault on:click={doMerge}
  ><Icon name="merge" /></button>

  <button class="tb-btn" title="Duplicate the cell" disabled={!activeCellObj}
          tabindex="-1" on:pointerdown|preventDefault on:click={doDuplicate}
  ><Icon name="duplicate" /></button>

  <button class="tb-btn tb-danger" title="Delete the cell" disabled={!activeCellObj}
          tabindex="-1" on:pointerdown|preventDefault on:click={doDelete}
  ><Icon name="trash" /></button>
</ToolbarGroup>

<!-- Context-sensitive: only for code cells, and hidden rather than disabled
     otherwise. The Text half arrives with Markdown text cells. -->
<ToolbarGroup label="Code" visible={showCodeGroup}>
  <button class="tb-btn" title="Indent" disabled={!codeView}
          tabindex="-1" on:pointerdown|preventDefault on:click={() => withView(indentCode)}
  ><Icon name="indent" /></button>

  <button class="tb-btn" title="Outdent" disabled={!codeView}
          tabindex="-1" on:pointerdown|preventDefault on:click={() => withView(outdentCode)}
  ><Icon name="outdent" /></button>

  <button class="tb-btn" title="Duplicate this line" disabled={!codeView}
          tabindex="-1" on:pointerdown|preventDefault on:click={() => withView(duplicateLine)}
  ><Icon name="duplicate" size={14} /></button>

  <button class="tb-btn tb-mono" title="Comment or uncomment the selection" disabled={!codeView}
          tabindex="-1" on:pointerdown|preventDefault on:click={() => withView(commentCode)}
  >(*=*)</button>
</ToolbarGroup>

<!-- One button, not two: Mathematica's sidebar group also carries a chat panel,
     and Mathilda has no chat. A second button that opened nothing would be worse
     than the asymmetry. -->
<ToolbarGroup label="Sidebar">
  <button
    class="tb-btn"
    class:active={$propertiesOpen}
    title={$propertiesOpen ? 'Hide properties' : 'Show properties'}
    aria-pressed={$propertiesOpen}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={() => propertiesOpen.update(v => !v)}
  ><Icon name="sidebar" /></button>
</ToolbarGroup>

<ToolbarGroup label="Notebook">
  <button
    class="tb-btn"
    title={$darkMode ? 'Switch to light theme' : 'Switch to dark theme'}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={() => darkMode.update(v => !v)}
  >{$darkMode ? '◑' : '☀'}</button>

  <!-- A static title: what the caret is on is only resolved when the menu opens,
       and the menu's own first item names it. -->
  <button
    class="tb-btn"
    title="Documentation"
    aria-haspopup="menu"
    aria-expanded={openMenu === 'docs'}
    tabindex="-1"
    bind:this={docsAnchor}
    on:pointerdown|preventDefault
    on:click={() => toggleMenu('docs')}
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
  .tb-btn.tb-danger:hover:not(:disabled) { color: var(--err); background: var(--surface-2); }
  /* Typographic marks stay text, so they should look like the code they insert. */
  .tb-btn.tb-mono { font-family: var(--mono); font-size: 0.66rem; letter-spacing: -0.02em; }

  /* A toggle that is currently the case — the arrangement buttons. */
  .tb-btn.on {
    background: var(--surface-3);
    border-color: var(--accent);
    color: var(--accent);
  }

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

  /* ---- Evaluation ---- */
  .tb-run { color: var(--ok); }
  .tb-run:hover:not(:disabled) { background: color-mix(in srgb, var(--ok) 14%, transparent); }

  /* A caret paired with a button: narrower, and visually attached to its left
     neighbour so the two read as one control with a dropdown. */
  .tb-caret {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 15px;
    height: var(--tb-btn-sz, 24px);
    padding: 0;
    margin-left: -3px;
    background: none;
    border: none;
    border-radius: 4px;
    color: var(--tb-caption);
    cursor: pointer;
  }
  .tb-caret:hover { background: var(--surface-2); color: var(--text-h); }

  /* Separates the run controls from the kernel chip inside one group. */
  .tb-mini-rule {
    width: 1px;
    height: 14px;
    margin: 0 4px;
    background: var(--tb-rule);
    flex-shrink: 0;
  }

  .tb-kernel {
    display: flex;
    align-items: center;
    gap: 5px;
    height: var(--tb-btn-sz, 24px);
    padding: 0 4px 0 7px;
    background: var(--surface-2);
    border: 1px solid var(--tb-rule);
    border-radius: 11px;
    color: var(--text);
    font: 500 0.72rem/1 var(--sans);
    cursor: pointer;
    white-space: nowrap;
  }
  .tb-kernel:hover { border-color: var(--accent); color: var(--text-h); }
  .tb-kernel-label { font-variant-numeric: tabular-nums; }

  .dot {
    width: 7px;
    height: 7px;
    border-radius: 50%;
    flex-shrink: 0;
    background: var(--text-muted);
  }
  [data-status='ready']      .dot { background: var(--ok); }
  [data-status='busy']       .dot { background: var(--warn); animation: blink 1s ease-in-out infinite; }
  [data-status='starting']   .dot { background: var(--accent); animation: blink 1s ease-in-out infinite; }
  [data-status='restarting'] .dot { background: var(--warn); animation: blink 1s ease-in-out infinite; }
  [data-status='dead']       .dot { background: var(--err); }
  [data-status='dead'] .tb-kernel-label { color: var(--err); }

  @keyframes blink { 0%, 100% { opacity: 1; } 50% { opacity: 0.3; } }

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
