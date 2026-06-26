<script lang="ts">
  import { onMount, onDestroy, createEventDispatcher } from 'svelte';
  import { EditorView, keymap } from '@codemirror/view';
  import { EditorState } from '@codemirror/state';
  import { defaultKeymap, history, historyKeymap } from '@codemirror/commands';
  import { syntaxHighlighting, defaultHighlightStyle } from '@codemirror/language';
  import Output from './Output.svelte';
  import type { OutputItem, CellStatus } from './notebook';
  import { selectedCells, selectOnly, toggleSelect, rangeSelect } from './notebook';

  export let cellId: string;
  export let source: string;
  export let status: CellStatus;
  export let output: OutputItem[];
  export let execIdx: number | undefined = undefined;  // In[n] label

  const dispatch = createEventDispatcher<{
    run:      { id: string };
    change:   { id: string; source: string };
    addBelow: { id: string };
  }>();

  let editorContainer: HTMLElement;
  let view: EditorView;
  $: selected = $selectedCells.has(cellId);

  onMount(() => {
    view = new EditorView({
      state: EditorState.create({
        doc: source,
        extensions: [
          history(),
          syntaxHighlighting(defaultHighlightStyle),
          keymap.of([
            { key: 'Shift-Enter', run() { dispatch('run', { id: cellId }); return true; } },
            { key: 'Mod-Enter',   run() { dispatch('run', { id: cellId }); dispatch('addBelow', { id: cellId }); return true; } },
            { key: 'Escape',      run() { editorContainer?.closest('.code-cell')?.blur(); return false; } },
            ...defaultKeymap,
            ...historyKeymap,
          ]),
          EditorView.updateListener.of(update => {
            if (update.docChanged)
              dispatch('change', { id: cellId, source: view.state.doc.toString() });
          }),
          EditorView.theme({
            '&': {
              fontSize: '14px',
              fontFamily: "'SF Mono', 'Fira Code', 'Cascadia Code', monospace",
              background: 'transparent',
            },
            '.cm-scroller': { overflow: 'auto', minHeight: '1.8em' },
            '.cm-content': { padding: '6px 0', textAlign: 'left' },
            '.cm-focused': { outline: 'none' },
            '.cm-line': { lineHeight: '1.6', textAlign: 'left' },
            '.cm-cursor': { borderLeftColor: 'var(--accent)' },
          }),
        ],
      }),
      parent: editorContainer,
    });
  });

  onDestroy(() => view?.destroy());

  $: if (view && view.state.doc.toString() !== source) {
    view.dispatch({ changes: { from: 0, to: view.state.doc.length, insert: source } });
  }

  function onGutterClick(e: MouseEvent) {
    // Don't interfere with the run button (it stops propagation itself).
    if (e.shiftKey) rangeSelect(cellId);
    else if (e.metaKey || e.ctrlKey) toggleSelect(cellId);
    else selectOnly(cellId);
  }

  function onCellClick() {
    // Clicking the body deselects all cells.
    import('./notebook').then(m => m.clearSelection());
  }
</script>

<!-- svelte-ignore a11y-click-events-have-key-events -->
<!-- svelte-ignore a11y-no-static-element-interactions -->
<div
  class="code-cell"
  class:running={status === 'running'}
  class:selected
  on:click={onCellClick}
>
  <!-- Left gutter — clickable for selection -->
  <!-- svelte-ignore a11y-click-events-have-key-events -->
  <div
    class="cell-gutter"
    on:click|stopPropagation={onGutterClick}
    title="Click to select • Shift+click to range select"
  >
    {#if execIdx != null}
      <span class="cell-label">In[{execIdx}]</span>
    {/if}
    <button
      class="run-btn"
      title="Run (Shift+Enter)"
      disabled={status === 'running'}
      on:click|stopPropagation={() => dispatch('run', { id: cellId })}
    >
      {#if status === 'running'}
        <span class="spinner">●</span>
      {:else}
        ▶
      {/if}
    </button>
  </div>

  <!-- Cell body -->
  <div class="cell-body">
    <div class="editor-wrap" bind:this={editorContainer}></div>
    {#if output.length > 0}
      <div class="output-label">
        {#if execIdx != null}<span class="out-tag">Out[{execIdx}]=</span>{/if}
      </div>
      <Output items={output} />
    {/if}
  </div>
</div>

<style>
  .code-cell {
    display: flex;
    flex-direction: row;
    border: 1px solid var(--border);
    border-radius: 6px;
    margin-bottom: 0.5rem;
    background: var(--cell-bg);
    transition: border-color 0.15s, box-shadow 0.15s;
    position: relative;
  }
  .code-cell:focus-within {
    border-color: var(--accent);
    box-shadow: 0 0 0 2px var(--accent-glow);
  }
  .code-cell.running  { border-color: #f39c12; }
  .code-cell.selected {
    border-color: var(--accent);
    box-shadow: -3px 0 0 0 var(--accent);
  }

  .cell-gutter {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: flex-start;
    padding: 0.4rem 0.3rem;
    gap: 0.2rem;
    background: var(--gutter-bg);
    border-radius: 5px 0 0 5px;
    border-right: 1px solid var(--border);
    min-width: 52px;
    cursor: pointer;
    user-select: none;
  }
  .cell-gutter:hover { background: var(--gutter-hover); }

  .cell-label {
    font-size: 0.68rem;
    color: var(--text-muted);
    font-family: 'SF Mono', monospace;
    white-space: nowrap;
  }

  .run-btn {
    background: none;
    border: none;
    cursor: pointer;
    font-size: 0.8rem;
    color: var(--accent);
    padding: 2px 4px;
    border-radius: 3px;
    transition: background 0.1s;
    line-height: 1;
  }
  .run-btn:hover:not(:disabled) { background: var(--accent-glow); }
  .run-btn:disabled { opacity: 0.3; cursor: default; }

  .spinner {
    display: inline-block;
    animation: pulse 0.9s ease-in-out infinite;
    color: #f39c12;
  }
  @keyframes pulse {
    0%, 100% { opacity: 1; }
    50%       { opacity: 0.3; }
  }

  .cell-body { flex: 1; min-width: 0; }
  .editor-wrap { padding: 4px 6px 4px 4px; }

  .output-label { padding: 0 0.5rem; }
  .out-tag {
    font-size: 0.68rem;
    color: var(--text-muted);
    font-family: 'SF Mono', monospace;
  }
</style>
