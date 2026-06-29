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
    import('plotly.js-dist-min').then((Plotly: any) => {
      const spec = data as any;
      const dark = document.documentElement.classList.contains('dark');

      // Apply dark theme overrides when in dark mode
      const layoutOverride = dark ? {
        plot_bgcolor:  '#1e1e2e',
        paper_bgcolor: '#1e1e2e',
        font: { color: '#cdd6f4' },
        xaxis: {
          ...(spec.layout?.xaxis ?? {}),
          gridcolor:    '#313244',
          zerolinecolor:'#6c7086',
          tickfont:     { color: '#cdd6f4' },
        },
        yaxis: {
          ...(spec.layout?.yaxis ?? {}),
          gridcolor:    '#313244',
          zerolinecolor:'#6c7086',
          tickfont:     { color: '#cdd6f4' },
        },
      } : {
        plot_bgcolor:  '#ffffff',
        paper_bgcolor: '#f5f5fa',
      };

      const layout = { ...(spec.layout ?? {}), ...layoutOverride };
      Plotly.react(node, spec.data ?? [spec], layout, {
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
    padding: 0.3rem 0.75rem 0.5rem;
    min-height: 1px;
    text-align: left;
  }
  .out-expr {
    font-size: 1.05em;
    padding: 0.25rem 0;
    color: var(--out-text, #222);
    text-align: left;
  }
  .out-error {
    color: #e74c3c;
    font-family: 'SF Mono', monospace;
    font-size: 0.88em;
    background: rgba(231,76,60,0.08);
    border-left: 3px solid #e74c3c;
    padding: 0.4rem 0.8rem;
    border-radius: 3px;
    text-align: left;
  }
  .out-stream {
    color: var(--text-muted);
    font-size: 0.84em;
    margin: 0;
    white-space: pre-wrap;
    font-family: 'SF Mono', 'Fira Code', monospace;
    text-align: left;
  }
  .out-plot {
    width: 100%;
    min-height: 320px;
  }
  .out-html {
    font-size: 0.95em;
    text-align: left;
    color: var(--text);
  }
</style>
