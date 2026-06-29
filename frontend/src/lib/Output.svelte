<script lang="ts">
  import type { OutputItem } from './notebook';
  import katex from 'katex';
  import 'katex/dist/katex.min.css';

  export let items: OutputItem[] = [];

  // Max height before output is collapsed with a "Show more" toggle
  const MAX_HEIGHT = 180; // px

  // Per-item expanded/overflow state
  let expanded: Record<number, boolean> = {};
  let overflows: Record<number, boolean> = {};

  // Svelte action: measures actual scrollHeight vs offsetHeight.
  // Triggers reactivity only when overflow state changes.
  function measureOverflow(node: HTMLElement, idx: number) {
    function check() {
      const does = node.scrollHeight > node.offsetHeight + 4;
      if (overflows[idx] !== does) {
        overflows[idx] = does;
        overflows = { ...overflows };
      }
    }
    // First check after mount, then again after async content (KaTeX, Plotly)
    requestAnimationFrame(check);
    const t = setTimeout(check, 150);
    return { destroy() { clearTimeout(t); } };
  }

  function renderKatex(text: string): string {
    try {
      return katex.renderToString(text, { throwOnError: false, displayMode: false });
    } catch {
      return `<code>${text}</code>`;
    }
  }

  function mountPlot(node: HTMLElement, data: object) {
    import('plotly.js-dist-min').then((Plotly: any) => {
      const spec = data as any;
      const dark = !document.documentElement.classList.contains('light');
      const layoutOverride = dark ? {
        plot_bgcolor:  '#181825', paper_bgcolor: '#181825',
        font: { color: '#cdd6f4' },
        xaxis: { ...(spec.layout?.xaxis ?? {}), gridcolor: '#313244', zerolinecolor: '#585b70', tickfont: { color: '#cdd6f4' } },
        yaxis: { ...(spec.layout?.yaxis ?? {}), gridcolor: '#313244', zerolinecolor: '#585b70', tickfont: { color: '#cdd6f4' } },
      } : {
        plot_bgcolor:  '#ffffff', paper_bgcolor: '#f5f5fa',
        font: { color: '#1c1c2e' },
        xaxis: { ...(spec.layout?.xaxis ?? {}), gridcolor: '#d8d9e8', zerolinecolor: '#9999bb', tickfont: { color: '#1c1c2e' } },
        yaxis: { ...(spec.layout?.yaxis ?? {}), gridcolor: '#d8d9e8', zerolinecolor: '#9999bb', tickfont: { color: '#1c1c2e' } },
      };
      Plotly.react(node, spec.data ?? [spec], { ...(spec.layout ?? {}), ...layoutOverride }, {
        responsive: true, displayModeBar: true,
      });
    });
  }

  // Measure height of a rendered output element to decide if it needs collapse
  function checkOverflow(node: HTMLElement, idx: number) {
    requestAnimationFrame(() => {
      if (node.scrollHeight > MAX_HEIGHT + 20) {
        // tall enough to warrant collapsing by default — nothing to do,
        // the CSS max-height handles it; the button appears via CSS
      }
    });
    return {};
  }
</script>

<div class="output">
  {#each items as item, idx (idx)}
    <div class="out-item" class:expanded={expanded[idx]}>
      {#if item.kind === 'expr'}
        <div class="out-collapsible" use:measureOverflow={idx}>
          <div class="out-expr">{@html renderKatex(item.text)}</div>
        </div>
      {:else if item.kind === 'error'}
        <div class="out-error">{item.text}</div>
      {:else if item.kind === 'stream'}
        <div class="out-collapsible" use:measureOverflow={idx}>
          <pre class="out-stream">{item.text}</pre>
        </div>
      {:else if item.kind === 'plot'}
        <div class="out-plot" use:mountPlot={item.data}></div>
      {:else if item.kind === 'html'}
        <div class="out-html">{@html item.html}</div>
      {/if}

      <!-- Toggle only when content actually overflows the cap -->
      {#if overflows[idx]}
        <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
        <div class="out-toggle" on:click={() => expanded[idx] = !expanded[idx]}>
          {expanded[idx] ? '▲ collapse' : '▼ show more'}
        </div>
      {/if}
    </div>
  {/each}
</div>

<style>
  .output {
    padding: 0.3rem 0.75rem 0.5rem;
    min-height: 1px;
    text-align: left;
    min-width: 0;    /* prevent output from pushing cell-content wider */
    overflow: hidden; /* clip anything that escapes a collapsible */
  }

  .out-item {
    position: relative;
    margin-bottom: 0.2rem;
    min-width: 0;
  }

  /* Collapsible wrapper: clips vertically, scrolls horizontally */
  .out-collapsible {
    width: 100%;       /* don't expand beyond parent */
    min-width: 0;      /* flex child must have this for overflow-x to work */
    max-height: 180px;
    overflow-x: auto;
    overflow-y: hidden;
    -webkit-mask-image: linear-gradient(to bottom, black 55%, transparent 100%);
    mask-image:         linear-gradient(to bottom, black 55%, transparent 100%);
  }

  .expanded .out-collapsible {
    max-height: none;
    overflow-y: visible;
    -webkit-mask-image: none;
    mask-image: none;
  }

  /* Show-more button: hidden by default, shown only when content overflows */
  .out-toggle {
    display: none;
    font-size: 0.68rem;
    color: var(--accent, #89b4fa);
    cursor: pointer;
    padding: 2px 0 0;
    user-select: none;
    transition: opacity 0.1s;
    text-align: left;
  }
  .out-toggle:hover { opacity: 0.7; }

  /* Toggle visibility is controlled by JS (overflows[] reactive dict) */

  /* Expression output — overflow handled by parent .out-collapsible */
  .out-expr {
    font-size: 1.05em;
    padding: 0.25rem 0;
    color: var(--out-text, #222);
    text-align: left;
  }

  /* Error output */
  .out-error {
    color: #e74c3c;
    font-family: 'SF Mono', monospace;
    font-size: 0.88em;
    background: rgba(231,76,60,0.08);
    border-left: 3px solid #e74c3c;
    padding: 0.4rem 0.8rem;
    border-radius: 3px;
    text-align: left;
    overflow-x: auto;
  }

  /* Stream (print) output */
  .out-stream {
    color: var(--text-muted);
    font-size: 0.84em;
    margin: 0;
    white-space: pre-wrap;
    word-break: break-all;
    font-family: 'SF Mono', 'Fira Code', monospace;
    text-align: left;
    overflow-x: hidden;
  }

  /* Plot output */
  .out-plot {
    width: 100%;
    min-height: 320px;
  }

  /* HTML output */
  .out-html {
    font-size: 0.95em;
    text-align: left;
    color: var(--text);
    overflow-x: auto;
  }
</style>
