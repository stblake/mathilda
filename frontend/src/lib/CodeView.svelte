<!--
  CodeView.svelte — a read-only CodeMirror view for a converted OUTPUT form.
  Reuses the input editor's machinery so output selection behaves identically:
    • math=true  (InputForm / FullForm): Mathilda highlighting + bottom-up
      structural (precedence-aware) click selection, plus Alt-Up/Down.
    • math=false (TeXForm / MathML): plain, ordinary text selection.
  Read-only (EditorState.readOnly) keeps the content selectable and copyable while
  rejecting edits; it is NOT editable.of(false), which would stop CodeMirror from
  rendering the structural selection.
-->
<script lang="ts">
  import { onDestroy } from 'svelte';
  import { EditorView, drawSelection } from '@codemirror/view';
  import { EditorState } from '@codemirror/state';
  import { mathildaHighlightExtension } from './mathildaLang';
  import { mathildaStructuralSelection } from './structuralSelection';

  export let doc = '';
  /** true for a Mathilda expression (InputForm/FullForm) — highlight + structural
   *  selection; false for raw markup (TeXForm/MathML) — plain text selection. */
  export let math = false;

  let el: HTMLElement;
  let view: EditorView | undefined;

  function build() {
    if (view) { view.destroy(); view = undefined; }
    if (!el) return;
    const extensions = [
      EditorState.readOnly.of(true),
      EditorView.lineWrapping,
      // Draw the selection ourselves so it renders in a read-only view regardless
      // of native-focus quirks (needed for the structural selection to be visible).
      drawSelection(),
      ...(math ? [...mathildaHighlightExtension, ...mathildaStructuralSelection] : []),
      EditorView.theme({
        '&': { background: 'transparent', fontSize: '0.95em' },
        '.cm-content': {
          padding: '2px 0', textAlign: 'left',
          fontFamily: "'SF Mono','Fira Code','Cascadia Code',monospace",
        },
        '.cm-focused': { outline: 'none' },
        '.cm-line': { padding: '0', lineHeight: '1.6' },
        '.cm-scroller': { overflow: 'visible', fontFamily: 'inherit' },
      }),
    ];
    view = new EditorView({ state: EditorState.create({ doc, extensions }), parent: el });
  }

  // (Re)build when the container mounts or the text/mode changes.
  $: if (el && (!view || view.state.doc.toString() !== doc)) build();

  onDestroy(() => view?.destroy());
</script>

<div class="code-view" bind:this={el}></div>

<style>
  .code-view { text-align: left; }
  .code-view :global(.cm-editor) { background: transparent; }
  .code-view :global(.cm-content) { color: var(--out-text, #cdd6f4); }
  /* Guarantee the output is selectable even if an ancestor sets user-select:none. */
  .code-view :global(.cm-content),
  .code-view :global(.cm-line) { user-select: text; -webkit-user-select: text; }
  .code-view :global(.cm-selectionBackground),
  .code-view :global(.cm-editor.cm-focused .cm-selectionBackground) {
    background: rgba(137,180,250,0.25) !important;
  }
</style>
