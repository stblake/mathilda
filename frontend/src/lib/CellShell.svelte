<!--
  CellShell.svelte
  Wraps any cell type. Provides:
    • Cell type picker (dropdown on the type badge)
    • 4-directional add buttons (↑ ↓ ← →) that appear on hover
    • Selection state (left accent bar)
    • Run / output area for code cells
-->
<script lang="ts">
  import { onMount, onDestroy, createEventDispatcher } from 'svelte';
  import { EditorView, keymap } from '@codemirror/view';
  import { EditorState, EditorSelection } from '@codemirror/state';
  import { defaultKeymap, history, historyKeymap } from '@codemirror/commands';
  import { syntaxHighlighting, defaultHighlightStyle } from '@codemirror/language';
  import Output from './Output.svelte';
  import type { Cell, CellType, OutputItem } from './notebook';
  import { selectedCells, selectOnly, toggleSelect, rangeSelect, clearSelection, notebook } from './notebook';

  export let cell: Cell;
  export let rowId: string;
  export let cellIdx: number;

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
            '.cm-content': { padding: '6px 8px', textAlign: 'left' },
            '.cm-focused': { outline: 'none' },
            '.cm-line': { lineHeight: '1.6', textAlign: 'left' },
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
    notebook.setCellType(cell.id, t);
    if (t === 'code') {
      // Re-init editor next tick.
      setTimeout(initEditor, 10);
    }
  }

  // ---- Selection ----
  function onBodyClick(e: MouseEvent) { clearSelection(); }
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

  <!-- Directional add buttons (visible on hover) -->
  <button class="add-btn add-top"    on:click|stopPropagation={() => dispatch('addAbove', { rowId })}>＋</button>
  <button class="add-btn add-bottom" on:click|stopPropagation={() => dispatch('addBelow', { rowId })}>＋</button>
  <button class="add-btn add-left"   on:click|stopPropagation={() => dispatch('addLeft',  { rowId, cellIdx })}>＋</button>
  <button class="add-btn add-right"  on:click|stopPropagation={() => dispatch('addRight', { rowId, cellIdx: cellIdx + 1 })}>＋</button>

  <!-- Cell header: type badge + run button + exec index -->
  <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
  <div class="cell-header" on:click={onHeaderClick}>

    <!-- Type badge / picker trigger -->
    <div class="type-badge-wrap">
      <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
      <button class="type-badge" on:click|stopPropagation={() => showTypePicker = !showTypePicker} title="Change cell type">
        {TYPES.find(t => t.id === cell.type)?.icon ?? '▶'}
      </button>

      {#if showTypePicker}
        <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
        <div class="type-picker" on:click|stopPropagation>
          {#each TYPES as t}
            <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
            <div
              class="type-option"
              class:active={cell.type === t.id}
              on:click={() => setType(t.id)}
            >
              <span class="type-opt-icon">{t.icon}</span>
              <span class="type-opt-label">{t.label}</span>
              <span class="type-opt-desc">{t.desc}</span>
            </div>
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
      <div
        class="prose-cell"
        contenteditable="true"
        on:input={onTextInput}
        on:click|stopPropagation
      >{cell.source}</div>

    {:else if cell.type === 'section'}
      <!-- svelte-ignore a11y-click-events-have-key-events -->
      <h1
        class="heading-cell"
        contenteditable="true"
        on:input={onTextInput}
        on:click|stopPropagation
      >{cell.source}</h1>

    {:else if cell.type === 'subsection'}
      <!-- svelte-ignore a11y-click-events-have-key-events -->
      <h2
        class="heading-cell"
        contenteditable="true"
        on:input={onTextInput}
        on:click|stopPropagation
      >{cell.source}</h2>
    {/if}
  </div>
</div>

<!-- Close type picker on outside click -->
{#if showTypePicker}
  <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
  <div class="picker-backdrop" on:click={() => showTypePicker = false}></div>
{/if}

<style>
  /* ---- Shell ---- */
  .cell-shell {
    position: relative;
    border-top: 1px solid transparent;
    border-bottom: 1px solid transparent;
    transition: border-color 0.12s, background 0.12s;
    background: var(--cell-bg);
  }
  .cell-shell:focus-within  { border-top-color: var(--accent); }
  .cell-shell.selected      { border-left: 3px solid var(--accent); }
  .cell-shell.running       { background: rgba(243,156,18,0.04); }

  /* Non-code cell backgrounds */
  .type-section    { background: transparent; }
  .type-subsection { background: transparent; }
  .type-text       { background: transparent; }

  /* ---- Directional add buttons ---- */
  .add-btn {
    position: absolute;
    display: none;
    align-items: center;
    justify-content: center;
    width: 20px;
    height: 20px;
    background: var(--accent);
    color: #fff;
    border: none;
    border-radius: 50%;
    font-size: 0.75rem;
    cursor: pointer;
    z-index: 20;
    opacity: 0.8;
    line-height: 1;
    padding: 0;
  }
  .cell-shell:hover .add-btn { display: flex; }
  .add-btn:hover { opacity: 1; transform: scale(1.1); }

  .add-top    { top: -10px;  left: 50%; transform: translateX(-50%); }
  .add-bottom { bottom: -10px; left: 50%; transform: translateX(-50%); }
  .add-left   { left: -10px; top: 50%; transform: translateY(-50%); }
  .add-right  { right: -10px; top: 50%; transform: translateY(-50%); }
  .add-top:hover    { transform: translateX(-50%) scale(1.1); }
  .add-bottom:hover { transform: translateX(-50%) scale(1.1); }
  .add-left:hover   { transform: translateY(-50%) scale(1.1); }
  .add-right:hover  { transform: translateY(-50%) scale(1.1); }

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
    color: var(--text-muted);
    font-family: 'SF Mono', monospace;
  }

  .run-btn {
    background: none;
    border: none;
    cursor: pointer;
    font-size: 0.75rem;
    color: var(--accent);
    padding: 1px 3px;
    border-radius: 3px;
    line-height: 1;
    transition: background 0.1s;
  }
  .run-btn:hover:not(:disabled) { background: var(--accent-glow); }
  .run-btn:disabled { opacity: 0.3; cursor: default; }

  .spinner { animation: pulse 0.9s ease-in-out infinite; display: inline-block; color: #f39c12; }
  @keyframes pulse { 0%,100%{opacity:1} 50%{opacity:0.3} }

  /* ---- Type badge / picker ---- */
  .type-badge-wrap { position: relative; }

  .type-badge {
    background: none;
    border: 1px solid var(--border);
    border-radius: 4px;
    font-size: 0.65rem;
    font-family: 'SF Mono', monospace;
    color: var(--text-muted);
    padding: 1px 5px;
    cursor: pointer;
    line-height: 1.4;
    transition: border-color 0.1s, color 0.1s;
  }
  .type-badge:hover { border-color: var(--accent); color: var(--accent); }

  .type-picker {
    position: absolute;
    top: calc(100% + 4px);
    left: 0;
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 8px;
    box-shadow: 0 8px 24px rgba(0,0,0,0.2);
    z-index: 100;
    min-width: 260px;
    overflow: hidden;
    padding: 4px 0;
  }

  .type-option {
    display: flex;
    align-items: center;
    gap: 0.6rem;
    padding: 0.45rem 0.9rem;
    cursor: pointer;
    transition: background 0.1s;
  }
  .type-option:hover  { background: var(--gutter-hover); }
  .type-option.active { background: var(--accent-glow); }

  .type-opt-icon {
    font-size: 0.8rem;
    font-family: 'SF Mono', monospace;
    color: var(--accent);
    width: 24px;
    text-align: center;
  }
  .type-opt-label {
    font-size: 0.88rem;
    font-weight: 500;
    color: var(--text);
    min-width: 80px;
  }
  .type-opt-desc {
    font-size: 0.75rem;
    color: var(--text-muted);
  }

  .picker-backdrop {
    position: fixed;
    inset: 0;
    z-index: 99;
  }

  /* ---- Cell content ---- */
  .cell-content { padding: 0; }

  /* prose / heading cells */
  .prose-cell {
    padding: 6px 8px;
    font-size: 0.95rem;
    color: var(--text);
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
    color: var(--text);
    outline: none;
    text-align: left;
    border: none;
    background: transparent;
    width: 100%;
  }
  h1.heading-cell { font-size: 1.5rem; border-bottom: 1px solid var(--border); padding-bottom: 0.3rem; }
  h2.heading-cell { font-size: 1.15rem; }

  .out-label {
    font-size: 0.67rem;
    color: var(--text-muted);
    font-family: 'SF Mono', monospace;
    padding: 0 8px;
    margin-top: 2px;
  }
</style>
