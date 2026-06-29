<!--
  CellShell.svelte
  Wraps any cell type. Provides:
    • Cell type picker (dropdown on the type badge)
    • 4-directional add buttons (↑ ↓ ← →) that appear on hover
    • Selection state (left accent bar)
    • Run / output area for code cells

  The `store` prop is the per-notebook store (createNotebook() instance).
  All mutation calls go through store.xxx() instead of the global notebook singleton.
-->
<script lang="ts">
  import { onMount, onDestroy, createEventDispatcher } from 'svelte';
  import { EditorView, keymap } from '@codemirror/view';
  import { EditorState, EditorSelection } from '@codemirror/state';
  import { defaultKeymap, history, historyKeymap } from '@codemirror/commands';
  import { syntaxHighlighting, defaultHighlightStyle } from '@codemirror/language';
  import Output from './Output.svelte';
  import type { Cell, CellType, OutputItem } from './notebook';
  import { selectedCells, selectOnly, toggleSelect, rangeSelect, clearSelection } from './notebook';

  export let cell: Cell;
  export let rowId: string;
  export let cellIdx: number;
  /** Per-notebook store instance — use store.xxx() for all mutations. */
  export let store: any;

  const dispatch = createEventDispatcher<{
    run:       { id: string };
    change:    { id: string; source: string };
    addAbove:  { rowId: string };
    addBelow:  { rowId: string };
    addLeft:   { rowId: string; cellIdx: number };
    addRight:  { rowId: string; cellIdx: number };
    focusPrev: { id: string };
    focusNext: { id: string };
    register:  { id: string; fn: () => void };
  }>();

  $: selected = $selectedCells.has(cell.id);

  // ---- CodeMirror editor ----
  let editorContainer: HTMLElement;
  let view: EditorView;

  onMount(() => {
    if (cell.type !== 'code') return;
    initEditor();
  });

  function initEditor() {
    if (view) { view.destroy(); view = undefined as any; }
    view = new EditorView({
      state: EditorState.create({
        doc: cell.source,
        extensions: [
          history(),
          syntaxHighlighting(defaultHighlightStyle),
          keymap.of([
            { key: 'Shift-Enter', run() { dispatch('run', { id: cell.id }); return true; } },
            { key: 'Mod-Enter',   run() { dispatch('run', { id: cell.id }); dispatch('addBelow', { rowId }); return true; } },
            { key: 'ArrowUp',     run(v) {
              const sel = v.state.selection.main;
              if (sel.head <= v.state.doc.lineAt(0).to) {
                dispatch('focusPrev', { id: cell.id }); return true;
              }
              return false;
            }},
            { key: 'ArrowDown',   run(v) {
              const sel = v.state.selection.main;
              const last = v.state.doc.lineAt(v.state.doc.length);
              if (sel.head >= last.from) {
                dispatch('focusNext', { id: cell.id }); return true;
              }
              return false;
            }},
            ...defaultKeymap,
            ...historyKeymap,
          ]),
          EditorView.updateListener.of(update => {
            if (update.docChanged)
              dispatch('change', { id: cell.id, source: view.state.doc.toString() });
          }),
          EditorView.theme({
            '&': { fontSize: '14px', fontFamily: "'SF Mono','Fira Code','Cascadia Code',monospace", background: 'transparent' },
            '.cm-scroller': { overflow: 'auto', minHeight: '1.8em' },
            '.cm-content': { padding: '6px 8px', textAlign: 'left', caretColor: '#89b4fa' },
            '.cm-focused': { outline: 'none' },
            '.cm-line': { lineHeight: '1.6', textAlign: 'left' },
            '.cm-cursor, .cm-dropCursor': { borderLeftColor: '#89b4fa !important', borderLeftWidth: '2px !important' },
            '.cm-selectionBackground': { background: 'rgba(137,180,250,0.2) !important' },
          }),
        ],
      }),
      parent: editorContainer,
    });

    dispatch('register', {
      id: cell.id,
      fn: () => {
        view.focus();
        const end = view.state.doc.length;
        view.dispatch({ selection: EditorSelection.cursor(end) });
      },
    });
  }

  onDestroy(() => view?.destroy());

  $: if (view && view.state.doc.toString() !== cell.source) {
    view.dispatch({ changes: { from: 0, to: view.state.doc.length, insert: cell.source } });
  }

  // ---- Type picker ----
  let showTypePicker = false;
  const TYPES: { id: CellType; label: string; icon: string; desc: string }[] = [
    { id: 'code',       icon: '▶',  label: 'Code',       desc: 'Evaluate Mathilda expressions' },
    { id: 'text',       icon: 'T',  label: 'Text',       desc: 'Prose / markdown' },
    { id: 'section',    icon: '#',  label: 'Section',    desc: 'H1 heading' },
    { id: 'subsection', icon: '##', label: 'Subsection', desc: 'H2 heading' },
  ];

  function setType(t: CellType) {
    showTypePicker = false;
    // Use the store prop instead of the global notebook singleton
    store.setCellType(cell.id, t);
    if (t === 'code') {
      setTimeout(initEditor, 10);
    }
  }

  // ---- Selection ----
  function onBodyClick(_e: MouseEvent) { clearSelection(); }
  function onHeaderClick(e: MouseEvent) {
    e.stopPropagation();
    if (e.shiftKey) rangeSelect(cell.id);
    else if (e.metaKey || e.ctrlKey) toggleSelect(cell.id);
    else selectOnly(cell.id);
  }

  // ---- Inline text editing (non-code) ----
  function onTextInput(e: Event) {
    dispatch('change', { id: cell.id, source: (e.target as HTMLElement).innerText });
  }

  // Focus ref + contenteditable state management
  let proseEl: HTMLElement;
  let _lastCellId = '';

  // Only update the DOM when the CELL IDENTITY changes (loading a new file),
  // not on every store update caused by the user typing.
  // If we update on every cell.source change, Svelte overwrites the user's
  // typed content, causing the duplication bug.
  $: if (proseEl && cell.id !== _lastCellId) {
    proseEl.innerText = cell.source;
    _lastCellId = cell.id;
    // Register focus fn once the element exists
    dispatch('register', {
      id: cell.id,
      fn: () => { proseEl?.focus(); },
    });
  }
</script>

<!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
<div
  class="cell-shell"
  class:selected
  class:running={cell.status === 'running'}
  class:type-code={cell.type === 'code'}
  class:type-text={cell.type === 'text'}
  class:type-section={cell.type === 'section'}
  class:type-subsection={cell.type === 'subsection'}
  on:click={onBodyClick}
>

  <!-- Directional add buttons removed for cleaner UI -->

  <!-- Cell header: type badge + run button + exec index -->
  <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
  <div class="cell-header" on:click={onHeaderClick}>

    <!-- Type badge / picker trigger -->
    <!-- Compact type switcher: click icon to open a horizontal pill strip -->
    <div class="type-badge-wrap">
      <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
      <button class="type-badge" on:click|stopPropagation={() => showTypePicker = !showTypePicker} title="Change cell type">
        {TYPES.find(t => t.id === cell.type)?.icon ?? '▶'}
      </button>

      {#if showTypePicker}
        <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
        <div class="type-picker-row" on:click|stopPropagation>
          {#each TYPES as t}
            <button
              class="type-pill"
              class:active={cell.type === t.id}
              on:click={() => setType(t.id)}
              title={t.label}
            >{t.icon} {t.label}</button>
          {/each}
        </div>
      {/if}
    </div>

    {#if cell.type === 'code'}
      {#if cell.execIdx != null}
        <span class="exec-label">In[{cell.execIdx}]</span>
      {/if}
      <button
        class="run-btn"
        title="Run (Shift+Enter)"
        disabled={cell.status === 'running'}
        on:click|stopPropagation={() => dispatch('run', { id: cell.id })}
      >
        {#if cell.status === 'running'}<span class="spinner">●</span>{:else}▶{/if}
      </button>
    {/if}
  </div>

  <!-- Cell content -->
  <div class="cell-content">
    {#if cell.type === 'code'}
      <div bind:this={editorContainer}></div>
      {#if cell.output.length > 0}
        {#if cell.execIdx != null}
          <div class="out-label">Out[{cell.execIdx}]=</div>
        {/if}
        <Output items={cell.output} />
      {/if}

    {:else if cell.type === 'text'}
      <!-- svelte-ignore a11y-click-events-have-key-events -->
      <!-- Content set via JS ($: proseEl update) to avoid contenteditable doubling -->
      <div
        class="prose-cell"
        contenteditable="true"
        bind:this={proseEl}
        on:input={onTextInput}
        on:click|stopPropagation
      ></div>

    {:else if cell.type === 'section'}
      <!-- svelte-ignore a11y-click-events-have-key-events -->
      <h1
        class="heading-cell"
        contenteditable="true"
        bind:this={proseEl}
        on:input={onTextInput}
        on:click|stopPropagation
      ></h1>

    {:else if cell.type === 'subsection'}
      <!-- svelte-ignore a11y-click-events-have-key-events -->
      <h2
        class="heading-cell"
        contenteditable="true"
        bind:this={proseEl}
        on:input={onTextInput}
        on:click|stopPropagation
      ></h2>
    {/if}
  </div>
</div>

<!-- Close type picker when clicking outside -->
{#if showTypePicker}
  <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
  <div style="position:fixed;inset:0;z-index:199;" on:click={() => showTypePicker = false}></div>
{/if}

<style>
  /* ---- Shell ---- */
  .cell-shell {
    position: relative;
    border-top: 1px solid transparent;
    border-bottom: 1px solid transparent;
    transition: border-color 0.12s, background 0.12s;
    background: var(--cell-bg, transparent);
  }
  /* Focus indicator: subtle left glow, not a full-width green stripe */
  .cell-shell:focus-within  { box-shadow: inset 2px 0 0 var(--accent, #89b4fa); }
  .cell-shell.selected      { border-left: 3px solid var(--accent, #89b4fa); }
  .cell-shell.running       { background: rgba(243,156,18,0.04); }

  /* Non-code cell backgrounds */
  .type-section    { background: transparent; }
  .type-subsection { background: transparent; }
  .type-text       { background: transparent; }

  /* add buttons removed */

  /* ---- Cell header ---- */
  .cell-header {
    display: flex;
    align-items: center;
    gap: 0.4rem;
    padding: 3px 6px 0;
    cursor: pointer;
    user-select: none;
  }

  .exec-label {
    font-size: 0.67rem;
    color: var(--text-muted, #585b70);
    font-family: 'SF Mono', monospace;
  }

  .run-btn {
    background: none;
    border: none;
    cursor: pointer;
    font-size: 0.75rem;
    color: var(--accent, #89b4fa);
    padding: 1px 3px;
    border-radius: 3px;
    line-height: 1;
    transition: background 0.1s;
  }
  .run-btn:hover:not(:disabled) { background: var(--accent-glow, rgba(137,180,250,0.12)); }
  .run-btn:disabled { opacity: 0.3; cursor: default; }

  .spinner { animation: pulse 0.9s ease-in-out infinite; display: inline-block; color: #f39c12; }
  @keyframes pulse { 0%,100%{opacity:1} 50%{opacity:0.3} }

  /* ---- Type badge / picker ---- */
  .type-badge-wrap { position: relative; }

  .type-badge {
    background: none;
    border: 1px solid var(--border, rgba(255,255,255,0.08));
    border-radius: 4px;
    font-size: 0.65rem;
    font-family: 'SF Mono', monospace;
    color: var(--text-muted, #585b70);
    padding: 1px 5px;
    cursor: pointer;
    line-height: 1.4;
    transition: border-color 0.1s, color 0.1s;
  }
  .type-badge:hover { border-color: var(--accent, #89b4fa); color: var(--accent, #89b4fa); }

  /* Compact horizontal type picker — appears as a pill row below the badge */
  .type-picker-row {
    position: absolute;
    top: calc(100% + 3px);
    left: 0;
    display: flex;
    gap: 3px;
    background: rgba(12, 15, 28, 0.96);
    border: 1px solid rgba(255,255,255,0.1);
    border-radius: 7px;
    padding: 4px;
    z-index: 200;
    box-shadow: 0 6px 20px rgba(0,0,0,0.5);
    white-space: nowrap;
  }

  .type-pill {
    background: none;
    border: 1px solid rgba(255,255,255,0.08);
    border-radius: 5px;
    color: #6c7086;
    font-size: 0.72rem;
    font-family: 'SF Mono', monospace;
    padding: 2px 8px;
    cursor: pointer;
    transition: all 0.1s;
    white-space: nowrap;
  }
  .type-pill:hover  { border-color: var(--accent, #89b4fa); color: var(--accent, #89b4fa); }
  .type-pill.active { border-color: var(--accent, #89b4fa); color: var(--accent, #89b4fa); background: rgba(137,180,250,0.1); }

  /* ---- Cell content ---- */
  .cell-content { padding: 0; }

  /* CodeMirror text colour for dark canvas */
  :global(.cm-editor .cm-content) { color: #cdd6f4; }
  :global(.cm-editor .cm-line)    { color: #cdd6f4; }

  /* prose / heading cells */
  .prose-cell {
    padding: 6px 8px;
    font-size: 0.95rem;
    color: #cdd6f4;
    min-height: 2em;
    outline: none;
    line-height: 1.6;
    text-align: left;
    white-space: pre-wrap;
  }
  .heading-cell {
    padding: 6px 8px;
    margin: 0;
    font-weight: 700;
    color: #cdd6f4;
    outline: none;
    text-align: left;
    border: none;
    background: transparent;
    width: 100%;
  }
  h1.heading-cell { font-size: 1.5rem; border-bottom: 1px solid rgba(255,255,255,0.08); padding-bottom: 0.3rem; }
  h2.heading-cell { font-size: 1.15rem; }

  .out-label {
    font-size: 0.67rem;
    color: #585b70;
    font-family: 'SF Mono', monospace;
    padding: 0 8px;
    margin-top: 2px;
    text-align: left;
    display: block;
  }
</style>
