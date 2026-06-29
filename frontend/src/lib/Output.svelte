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

  // Heuristic: expressions that are long lists of numbers/symbols don't
  // benefit from KaTeX (no fractions/superscripts) and KaTeX can't wrap them.
  // Render as code with word-break so they don't overflow the card.
  // Heuristic: if output has >4 commas or is long, it's a list/sequence.
  // KaTeX can't wrap math spans so we use plain code with word-break.
  function isListOutput(text: string): boolean {
    const commas = (text.match(/,/g) ?? []).length;
    return commas > 4 || text.length > 200;
  }

  function renderOutput(text: string, latex?: string): string {
    // Prefer LaTeX from the kernel (StandardForm serialiser) when available
    if (latex && latex.length > 0) {
      try {
        return katex.renderToString(latex, { throwOnError: false, displayMode: false });
      } catch {
        /* fall through to text rendering */
      }
    }
    // Long lists: render as wrapping code, not KaTeX (which can't wrap)
    if (isListOutput(text)) {
      const wrapped = text.replace(/,\s+/g, ', ');
      return `<code class="out-code-wrap">${wrapped}</code>`;
    }
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
          <div class="out-expr">{@html renderOutput(item.text, item.latex)}</div>
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

      <!-- Always show toggle for collapsible output so user can expand/collapse -->
      {#if item.kind !== 'plot' && item.kind !== 'error'}
        <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
        <div class="out-toggle" class:hidden={!overflows[idx] && !expanded[idx]} on:click={() => expanded[idx] = !expanded[idx]}>
          {expanded[idx] ? '▲ collapse' : '▼ show all'}
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
    font-size: 0.68rem;
    color: var(--accent, #89b4fa);
    cursor: pointer;
    padding: 2px 0 0;
    user-select: none;
    transition: opacity 0.1s;
    text-align: left;
  }
  .out-toggle:hover { opacity: 0.7; }
  /* Hide when content fits and not yet expanded */
  .out-toggle.hidden { display: none; }

  /* Expression output — overflow handled by parent .out-collapsible */
  .out-expr {
    font-size: 1.05em;
    padding: 0.25rem 0;
    color: var(--out-text, #222);
    text-align: left;
  }

  /* Long list/sequence outputs rendered as wrapping code */
  :global(.out-code-wrap) {
    display: block;
    font-family: 'SF Mono', 'Fira Code', monospace;
    font-size: 0.95em;
    color: var(--out-text, #cdd6f4);
    background: transparent;  /* override browser default <code> background */
    white-space: normal;
    word-break: break-word;
    overflow-wrap: break-word;
    line-height: 1.7;
    padding: 0.15rem 0;
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
