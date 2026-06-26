<script lang="ts">
  import type { OutputItem } from './notebook';
  import katex from 'katex';
  import 'katex/dist/katex.min.css';
  import { onMount } from 'svelte';

  export let items: OutputItem[] = [];

  function renderKatex(text: string): string {
    try {
      // Try block-mode first if it looks like a full expression.
      return katex.renderToString(text, {
        throwOnError: false,
        displayMode: false,
      });
    } catch {
      return `<code>${text}</code>`;
    }
  }

  function mountPlot(node: HTMLElement, data: object) {
    // Dynamically import Plotly to keep initial bundle lean.
    import('plotly.js-dist-min').then((Plotly: any) => {
      const spec = data as any;
      Plotly.react(node, spec.data ?? [spec], spec.layout ?? {}, {
        responsive: true,
        displayModeBar: true,
      });
    });
  }
</script>

<div class="output">
  {#each items as item (item)}
    {#if item.kind === 'expr'}
      <div class="out-expr">
        <!-- KaTeX render of the expression string -->
        {@html renderKatex(item.text)}
      </div>
    {:else if item.kind === 'error'}
      <div class="out-error">{item.text}</div>
    {:else if item.kind === 'stream'}
      <pre class="out-stream">{item.text}</pre>
    {:else if item.kind === 'plot'}
      <div class="out-plot" use:mountPlot={item.data}></div>
    {:else if item.kind === 'html'}
      <div class="out-html">{@html item.html}</div>
    {/if}
  {/each}
</div>

<style>
  .output {
    padding: 0.5rem 1rem 0.5rem 1rem;
    min-height: 1px;
  }
  .out-expr {
    font-size: 1.1em;
    padding: 0.25rem 0;
    color: #1a1a2e;
  }
  .out-error {
    color: #c0392b;
    font-family: monospace;
    font-size: 0.9em;
    background: #fdf0f0;
    border-left: 3px solid #c0392b;
    padding: 0.4rem 0.8rem;
    border-radius: 3px;
  }
  .out-stream {
    color: #555;
    font-size: 0.85em;
    margin: 0;
    white-space: pre-wrap;
    font-family: 'SF Mono', 'Fira Code', monospace;
  }
  .out-plot {
    width: 100%;
    min-height: 300px;
  }
  .out-html {
    font-size: 0.95em;
  }
</style>
