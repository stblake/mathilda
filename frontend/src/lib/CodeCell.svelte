<script lang="ts">
  import { onMount, onDestroy, createEventDispatcher } from 'svelte';
  import { EditorView, keymap, lineNumbers } from '@codemirror/view';
  import { EditorState } from '@codemirror/state';
  import { defaultKeymap, history, historyKeymap } from '@codemirror/commands';
  import { syntaxHighlighting, defaultHighlightStyle } from '@codemirror/language';
  import Output from './Output.svelte';
  import type { OutputItem, CellStatus } from './notebook';

  export let cellId: string;
  export let source: string;
  export let status: CellStatus;
  export let output: OutputItem[];

  const dispatch = createEventDispatcher<{
    run: { id: string };
    change: { id: string; source: string };
    addBelow: { id: string };
    remove: { id: string };
  }>();

  let editorContainer: HTMLElement;
  let view: EditorView;

  onMount(() => {
    view = new EditorView({
      state: EditorState.create({
        doc: source,
        extensions: [
          lineNumbers(),
          history(),
          syntaxHighlighting(defaultHighlightStyle),
          keymap.of([
            // Shift+Enter: run cell
            {
              key: 'Shift-Enter',
              run() {
                dispatch('run', { id: cellId });
                return true;
              },
            },
            // Ctrl+Enter: run and insert new cell below
            {
              key: 'Mod-Enter',
              run() {
                dispatch('run', { id: cellId });
                dispatch('addBelow', { id: cellId });
                return true;
              },
            },
            ...defaultKeymap,
            ...historyKeymap,
          ]),
          EditorView.updateListener.of(update => {
            if (update.docChanged) {
              dispatch('change', { id: cellId, source: view.state.doc.toString() });
            }
          }),
          EditorView.theme({
            '&': {
              fontSize: '14px',
              fontFamily: "'SF Mono', 'Fira Code', 'Cascadia Code', monospace",
              backgroundColor: '#fafafa',
              borderRadius: '4px',
            },
            '.cm-scroller': { overflow: 'auto', minHeight: '2em' },
            '.cm-content': { padding: '8px 4px' },
            '.cm-focused': { outline: 'none' },
            '.cm-line': { lineHeight: '1.6' },
          }),
        ],
      }),
      parent: editorContainer,
    });
  });

  onDestroy(() => view?.destroy());

  // Keep the editor in sync if source is changed externally.
  $: if (view && view.state.doc.toString() !== source) {
    view.dispatch({
      changes: { from: 0, to: view.state.doc.length, insert: source },
    });
  }
</script>

<div class="code-cell" class:running={status === 'running'}>
  <div class="cell-gutter">
    <button
      class="run-btn"
      title="Run (Shift+Enter)"
      disabled={status === 'running'}
      on:click={() => dispatch('run', { id: cellId })}
    >
      {#if status === 'running'}
        <span class="spinner">⏳</span>
      {:else}
        ▶
      {/if}
    </button>
    <button
      class="remove-btn"
      title="Delete cell"
      on:click={() => dispatch('remove', { id: cellId })}
    >✕</button>
  </div>

  <div class="cell-body">
    <div class="editor-wrap" bind:this={editorContainer}></div>
    {#if output.length > 0}
      <Output items={output} />
    {/if}
  </div>
</div>

<style>
  .code-cell {
    display: flex;
    flex-direction: row;
    border: 1px solid #e0e0e0;
    border-radius: 6px;
    margin-bottom: 0.75rem;
    background: #fff;
    transition: border-color 0.15s;
  }
  .code-cell:focus-within {
    border-color: #4a90e2;
    box-shadow: 0 0 0 2px rgba(74,144,226,0.15);
  }
  .code-cell.running {
    border-color: #f39c12;
  }
  .cell-gutter {
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 0.5rem 0.4rem;
    gap: 0.25rem;
    background: #f5f5f5;
    border-radius: 5px 0 0 5px;
    border-right: 1px solid #e0e0e0;
    min-width: 36px;
  }
  .run-btn {
    background: none;
    border: none;
    cursor: pointer;
    font-size: 0.85rem;
    color: #4a90e2;
    padding: 2px 4px;
    border-radius: 3px;
    transition: background 0.1s;
  }
  .run-btn:hover:not(:disabled) { background: #e8f0fb; }
  .run-btn:disabled { opacity: 0.4; cursor: default; }
  .remove-btn {
    background: none;
    border: none;
    cursor: pointer;
    font-size: 0.7rem;
    color: #aaa;
    padding: 2px 4px;
    border-radius: 3px;
    transition: color 0.1s;
  }
  .remove-btn:hover { color: #c0392b; }
  .cell-body {
    flex: 1;
    min-width: 0;
  }
  .editor-wrap {
    padding: 4px;
  }
  .spinner { animation: spin 1s linear infinite; display: inline-block; }
  @keyframes spin { to { transform: rotate(360deg); } }
</style>
