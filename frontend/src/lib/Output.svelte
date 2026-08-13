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
      // While expanded the clamp is lifted, so scrollHeight == offsetHeight and
      // the node reports "no overflow" -- which would hide the collapse toggle
      // and strand the user in the expanded state. Measure against the clamp
      // instead, so the answer means the same thing in both states.
      const does = expanded[idx]
        ? node.scrollHeight > MAX_HEIGHT + 4
        : node.scrollHeight > node.offsetHeight + 4;
      if (overflows[idx] !== does) {
        overflows[idx] = does;
        overflows = { ...overflows };
      }
    }
    // First check after mount, then again after async content (KaTeX, Plotly).
    requestAnimationFrame(check);
    const t = setTimeout(check, 150);

    // Re-measure on resize. Without this the overflow flag is decided once, at
    // the width the card happened to have on mount: widening the card can make
    // text that needed the clamp fit, and narrowing it can make text that fit
    // need the clamp, but neither re-ran the check -- so the toggle went stale
    // and stopped matching what was on screen.
    let ro: ResizeObserver | undefined;
    if (typeof ResizeObserver !== 'undefined') {
      ro = new ResizeObserver(() => check());
      ro.observe(node);
    }
    return {
      update() { check(); },
      destroy() { clearTimeout(t); ro?.disconnect(); }
    };
  }

  /* A usage message is structured, not a blob: alternating signature lines
     (flush left) and their descriptions (indented, hard-wrapped at ~70 columns
     for the terminal REPL). Rendering it as one <pre> loses all of that -- the
     hard wraps rewrap again at the card edge, a wrapped continuation drops back
     to column 0 so the indent stops meaning anything, and prose set in a
     monospace column is simply harder to read than it needs to be.

     Parse it back into blocks instead: each flush-left line is a signature,
     each run of indented lines is one description paragraph. The renderer then
     sets signatures in mono and descriptions as ordinary text, which is how
     Mathematica displays ?sym too. */
  type UsageBlock = { kind: 'sig' | 'body'; text: string };

  function parseUsage(text: string): UsageBlock[] {
    const out: UsageBlock[] = [];
    for (const line of text.split('\n')) {
      if (!line.trim()) continue;
      const indented = /^\s+\S/.test(line);
      const prev = out.length ? out[out.length - 1] : null;
      if (indented && prev && prev.kind === 'body') {
        prev.text += ' ' + line.trim();          // continuation of the paragraph
      } else if (indented) {
        out.push({ kind: 'body', text: line.trim() });
      } else {
        out.push({ kind: 'sig', text: line.trim() });
      }
    }
    return out;
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
    // Long lists: always use wrapping code regardless of latex field.
    // KaTeX renders math spans without line-breaking, so even \{1,2,...\}
    // produces a single wide unbreakable line.
    if (isListOutput(text)) {
      const wrapped = text.replace(/,\s+/g, ', ');
      return `<code class="out-code-wrap">${wrapped}</code>`;
    }
    // Short expressions: prefer LaTeX from the kernel (StandardForm)
    if (latex && latex.length > 0) {
      try {
        return katex.renderToString(latex, { throwOnError: false, displayMode: false });
      } catch { /* fall through */ }
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
      const sceneAxisDark  = { gridcolor: '#313244', zerolinecolor: '#585b70', tickfont: { color: '#cdd6f4' }, backgroundcolor: '#181825', showbackground: true };
      const sceneAxisLight = { gridcolor: '#d8d9e8', zerolinecolor: '#9999bb', tickfont: { color: '#1c1c2e' }, backgroundcolor: '#f0f0f8', showbackground: true };
      const layoutOverride = dark ? {
        plot_bgcolor:  '#181825', paper_bgcolor: '#181825',
        font: { color: '#cdd6f4' },
        xaxis: { ...(spec.layout?.xaxis ?? {}), gridcolor: '#313244', zerolinecolor: '#585b70', tickfont: { color: '#cdd6f4' } },
        yaxis: { ...(spec.layout?.yaxis ?? {}), gridcolor: '#313244', zerolinecolor: '#585b70', tickfont: { color: '#cdd6f4' } },
        scene: { ...(spec.layout?.scene ?? {}), bgcolor: '#181825', xaxis: sceneAxisDark, yaxis: sceneAxisDark, zaxis: sceneAxisDark },
      } : {
        plot_bgcolor:  '#ffffff', paper_bgcolor: '#f5f5fa',
        font: { color: '#1c1c2e' },
        xaxis: { ...(spec.layout?.xaxis ?? {}), gridcolor: '#d8d9e8', zerolinecolor: '#9999bb', tickfont: { color: '#1c1c2e' } },
        yaxis: { ...(spec.layout?.yaxis ?? {}), gridcolor: '#d8d9e8', zerolinecolor: '#9999bb', tickfont: { color: '#1c1c2e' } },
        scene: { ...(spec.layout?.scene ?? {}), bgcolor: '#f5f5fa', xaxis: sceneAxisLight, yaxis: sceneAxisLight, zaxis: sceneAxisLight },
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
    <div class="out-item" class:expanded={expanded[idx]} class:overflowing={overflows[idx]}>
      {#if item.kind === 'expr'}
        <div class="out-collapsible" use:measureOverflow={idx}>
          <div class="out-expr">{@html renderOutput(item.text, item.latex)}</div>
        </div>
      {:else if item.kind === 'expected'}
        <!-- A reference-page example that has not been run yet. Shown as plain
             text, never typeset: the recorded value is Mathilda syntax, and
             handing it to KaTeX (which is what an `expr` with no kernel LaTeX
             does) renders Derivative[1][g][x] as italic mathematics. -->
        <div class="out-collapsible" use:measureOverflow={idx}>
          <pre class="out-expected">{item.text}</pre>
        </div>
      {:else if item.kind === 'usage'}
        <!-- A usage message is documentation, not an expression: render it
             verbatim. It must not reach KaTeX, which cannot typeset it and
             would fall back to showing the InputForm escapes. -->
        <div class="out-collapsible" use:measureOverflow={idx}>
          <div class="out-usage">
            {#each parseUsage(item.text) as blk, bi (bi)}
              {#if blk.kind === 'sig'}
                <div class="usage-sig">{blk.text}</div>
              {:else}
                <p class="usage-body">{blk.text}</p>
              {/if}
            {/each}
          </div>
        </div>
      {:else if item.kind === 'names'}
        <!-- `?pat*` is a symbol search. A grid reads far better than one long
             braced line, and each name is a discrete thing to scan for. -->
        <div class="out-collapsible" use:measureOverflow={idx}>
          {#if item.names.length === 0}
            <div class="names-empty">No symbol matches that pattern.</div>
          {:else}
            <div class="names-grid">
              {#each item.names as nm (nm)}
                <button
                  class="name-chip"
                  title={`Open the reference page for ${nm}`}
                  on:click={(e) => e.currentTarget.dispatchEvent(
                    new CustomEvent('mathilda-refpage',
                      { detail: { name: nm }, bubbles: true }))}
                >{nm}</button>
              {/each}
            </div>
          {/if}
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
    width: 100%;
    min-width: 0;
    max-height: 180px;
    overflow-x: auto;
    overflow-y: hidden;
    /* Fade applied only when content actually overflows (via .overflowing class) */
  }

  /* Only fade when content genuinely overflows the cap */
  .overflowing .out-collapsible {
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

  /* Usage message from `?sym`, rendered as structure rather than a blob:
     signatures in mono so they read as code, descriptions as ordinary text so
     a paragraph reads as a paragraph. The body is indented as a block, so a
     wrapped line stays aligned under the first instead of falling back to
     column 0 the way a <pre> would. */
  .out-usage {
    text-align: left;
    margin: 0.1rem 0 0.2rem;
  }

  .usage-sig {
    font-family: 'SF Mono', 'Fira Code', monospace;
    font-size: 0.84em;
    color: var(--text);
    margin-top: 0.55rem;
    word-break: break-word;
  }
  .usage-sig:first-child { margin-top: 0; }

  .usage-body {
    margin: 0.15rem 0 0 1.6em;
    font-size: 0.88em;
    line-height: 1.5;
    color: var(--text-muted);
    max-width: 78ch;          /* prose stops being readable much past this */
    word-break: break-word;
  }

  /* `?pat*` symbol search: an auto-fitting grid, so a wide card shows more
     columns instead of one tall column. minmax keeps a long name like
     NeighborhoodContraction from being clipped while still packing short
     ones. */
  .names-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(15rem, 1fr));
    gap: 0.15rem 0.9rem;
    text-align: left;
  }

  /* A <button>, not a <code>: the global `code` rule paints a themed
     background and padding, which turned each name into a filled box whose
     text stopped contrasting with it.

     Colour is `inherit`, deliberately, not a theme token. app.css defines its
     dark palette only under @media (prefers-color-scheme: dark) and has no
     .light override, but the app also has its own light/dark toggle that sets
     a class on <html>. Toggle the app to light while the OS is dark and every
     var() still holds its DARK value on a light surface -- which is why
     --accent read as washed out here and --text-h read as white-on-white.
     Inheriting from the surrounding output text sidesteps that entirely: these
     names are exactly as legible as the text beside them, in every combination
     of OS and app theme. The underline, not colour, carries the affordance. */
  .out-expected {
    font: 0.82rem/1.6 var(--mono);
    color: var(--text-dim, var(--text));
    margin: 0;
    padding: 0;
    white-space: pre-wrap;
    word-break: break-word;
    opacity: 0.85;
  }

  .names-grid .name-chip {
    font-family: 'SF Mono', 'Fira Code', monospace;
    font-size: 0.84em;
    color: inherit;
    background: none;
    border: none;
    padding: 0.1rem 0;
    margin: 0;
    text-align: left;
    cursor: pointer;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    text-decoration: underline;
    text-decoration-color: color-mix(in srgb, currentColor 40%, transparent);
    text-underline-offset: 2px;
  }
  .names-grid .name-chip:hover,
  .names-grid .name-chip:focus-visible {
    text-decoration-color: currentColor;
    text-decoration-thickness: 2px;
  }

  .names-empty {
    color: var(--text-muted);
    font-size: 0.86em;
    text-align: left;
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
