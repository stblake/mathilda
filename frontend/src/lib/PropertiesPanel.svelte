<!--
  PropertiesPanel.svelte — the sidebar that slides in from the left of the
  focused view.

  It answers questions about the notebook you are looking at that nothing else in
  the window can: how big it actually is, whether the kernel is alive, and the two
  display preferences that are settings rather than actions. The toolbar is for
  verbs; a preference that persists for the rest of the session does not belong
  in a strip of buttons, which is why the In[n] label toggle lives here.

  WHAT IT DELIBERATELY DOES NOT CLAIM. A notebook on this canvas has a title and
  no file: `saveNotebook` takes a path from a dialog at the moment of saving, and
  nothing writes that path back into the model. So the Location row says the
  notebook is unsaved rather than inventing a path, and Size is measured from the
  content -- cells, and the characters of source they hold -- which is a real
  number about a real object. A "42 KB" for a file that does not exist would be a
  fabrication, and the panel's whole value is that its numbers can be trusted.
-->
<script lang="ts">
  import { onDestroy } from 'svelte';
  import Icon from './Icon.svelte';
  import { canvasState, activeActions, setFocusLayout } from './canvas';
  import type { FocusLayout } from './canvas';
  import { kernelStatus } from './notebook';
  import type { Cell, NotebookRow } from './notebook';
  import { restart, abortEvaluation } from './kernelActions';
  import { darkMode } from './theme';
  import { propertiesOpen, showExecLabels, uiScale, UI_SCALE_STEPS } from './properties';

  const KERNEL_LABEL: Record<string, string> = {
    starting:   'Starting…',
    ready:      'Ready',
    busy:       'Running',
    restarting: 'Restarting…',
    dead:       'Not running',
  };

  /* The active pane's rows, subscribed by hand: the store identity changes when
     the active pane changes, and `$store` auto-subscription only works on a fixed
     identifier. Same reason StatusBar does it this way. */
  let rows: NotebookRow[] = [];
  let unsub: (() => void) | null = null;
  $: {
    unsub?.();
    unsub = null;
    const store = $activeActions?.store;
    if (store) unsub = store.subscribe((r: NotebookRow[]) => { rows = r; });
    else rows = [];
  }
  onDestroy(() => unsub?.());

  $: cells = rows.flatMap(r => r.cells) as Cell[];
  $: codeCells  = cells.filter(c => c.type === 'code').length;
  $: textCells  = cells.length - codeCells;
  $: evaluated  = cells.filter(c => c.execIdx != null).length;
  /* Characters of source, which is what the notebook actually weighs. Outputs are
     excluded on purpose: they are recomputable and would make the number jump
     around as cells run, which reads as a bug rather than a measurement. */
  $: sourceChars = cells.reduce((n, c) => n + (c.source?.length ?? 0), 0);

  $: title = $activeActions
    ? ($canvasState.notebooks.find(nb => nb.id === $activeActions!.notebookId)?.title ?? '—')
    : '—';
  $: isRefpage = $activeActions
    ? ($canvasState.notebooks.find(nb => nb.id === $activeActions!.notebookId)?.refpage ?? false)
    : false;

  $: paneCount = $canvasState.focusedIds.length;
  $: layout = $canvasState.focusedLayout as FocusLayout;

  function fmtChars(n: number): string {
    if (n < 1000) return `${n} characters`;
    return `${(n / 1000).toFixed(1)}k characters`;
  }
</script>

<!-- aria-hidden while closed so the panel's controls leave the tab order with
     it; the transform alone would keep them focusable off-screen. -->
<aside
  class="props"
  class:open={$propertiesOpen}
  aria-hidden={!$propertiesOpen}
  aria-label="Notebook properties"
>
  <header class="props-head">
    <span class="props-title">Properties</span>
    <button
      class="props-close"
      title="Close properties"
      tabindex={$propertiesOpen ? 0 : -1}
      on:pointerdown|preventDefault
      on:click={() => propertiesOpen.set(false)}
    ><Icon name="close" /></button>
  </header>

  <section class="props-sec">
    <h3>Notebook</h3>
    <div class="row"><span class="k">Name</span><span class="v" title={title}>{title}</span></div>
    <div class="row">
      <span class="k">Location</span>
      <!-- Honest: a canvas notebook has no file until a save dialog gives it one,
           and that path is never written back into the model. -->
      <span class="v muted">{isRefpage ? 'Generated reference page' : 'Not saved to a file'}</span>
    </div>
    <div class="row"><span class="k">Size</span><span class="v">{fmtChars(sourceChars)}</span></div>
  </section>

  <section class="props-sec">
    <h3>Cells</h3>
    <div class="row"><span class="k">Total</span><span class="v">{cells.length}</span></div>
    <div class="row"><span class="k">Code</span><span class="v">{codeCells}</span></div>
    <div class="row"><span class="k">Text</span><span class="v">{textCells}</span></div>
    <div class="row">
      <span class="k">Evaluated</span>
      <span class="v">{evaluated} of {codeCells}</span>
    </div>
  </section>

  <section class="props-sec">
    <h3>Kernel</h3>
    <div class="row">
      <span class="k">Status</span>
      <span class="v" class:bad={$kernelStatus === 'dead'}>
        {KERNEL_LABEL[$kernelStatus] ?? $kernelStatus}
      </span>
    </div>
    <div class="btn-row">
      <button
        class="props-btn"
        tabindex={$propertiesOpen ? 0 : -1}
        title="Restart the kernel, discarding every definition in the session"
        on:pointerdown|preventDefault
        on:click={() => restart()}
      ><Icon name="restart" /> Restart</button>
      <!-- Disabled unless there is something to abort: abortEvaluation restarts
           the kernel, so offering it against an idle kernel invites throwing away
           a session's definitions for no reason. -->
      <button
        class="props-btn"
        tabindex={$propertiesOpen ? 0 : -1}
        disabled={$kernelStatus !== 'busy' && $kernelStatus !== 'restarting'}
        title="Abort the running evaluation"
        on:pointerdown|preventDefault
        on:click={() => abortEvaluation()}
      ><Icon name="abort" /> Abort</button>
    </div>
  </section>

  <section class="props-sec">
    <h3>Display</h3>
    <label class="check">
      <input
        type="checkbox"
        tabindex={$propertiesOpen ? 0 : -1}
        checked={$darkMode}
        on:change={() => darkMode.update(v => !v)}
      />
      Dark theme
    </label>
    <label class="check">
      <input
        type="checkbox"
        tabindex={$propertiesOpen ? 0 : -1}
        checked={$showExecLabels}
        on:change={() => showExecLabels.update(v => !v)}
      />
      Show In[n] labels
    </label>

    <!-- Scale drives the SAME store Cmd+= / Cmd+- / Cmd+0 have always driven; the
         panel is a second way in, not a second setting. A step reads as selected
         only on an exact match, because the keyboard moves in 0.1 and would
         otherwise light up the nearest button while sitting between two. -->
    <div class="row scale-row">
      <span class="k">Scale</span>
      <span class="scale-steps">
        {#each UI_SCALE_STEPS as step}
          <button
            class="props-btn scale-btn"
            class:on={$uiScale === step}
            tabindex={$propertiesOpen ? 0 : -1}
            title={`Set the interface scale to ${Math.round(step * 100)}%`}
            on:pointerdown|preventDefault
            on:click={() => uiScale.set(step)}
          >{Math.round(step * 100)}%</button>
        {/each}
      </span>
    </div>
  </section>

  <!-- Layout is only a question when there is more than one pane; with one pane
       the buttons would be inert controls, so the section is absent instead. -->
  {#if paneCount > 1}
    <section class="props-sec">
      <h3>Layout</h3>
      <div class="btn-row">
        <button
          class="props-btn"
          class:on={layout === 'h'}
          tabindex={$propertiesOpen ? 0 : -1}
          title="Side by side"
          on:pointerdown|preventDefault
          on:click={() => setFocusLayout('h')}
        ><Icon name="layoutH" /> Columns</button>
        <button
          class="props-btn"
          class:on={layout === 'v'}
          tabindex={$propertiesOpen ? 0 : -1}
          title="Stacked"
          on:pointerdown|preventDefault
          on:click={() => setFocusLayout('v')}
        ><Icon name="layoutV" /> Rows</button>
      </div>
      <div class="row"><span class="k">Panes</span><span class="v">{paneCount}</span></div>
    </section>
  {/if}
</aside>

<style>
  /* Slides in from the LEFT, under the toolbar and above the panes. Transformed
     rather than width-animated: animating width reflows the pane grid on every
     frame, and the panes contain editors. */
  .props {
    position: absolute;
    top: var(--toolbar-h, 46px);      /* the FOCUSED strip is 46px; --appbar-h is the canvas's 34 */
    left: 0;
    bottom: 0;
    width: 232px;
    z-index: 40;
    display: flex;
    flex-direction: column;
    gap: 2px;
    overflow-y: auto;
    padding: 10px 12px 16px;
    box-sizing: border-box;
    background: var(--surface);
    border-right: 1px solid var(--border);
    transform: translateX(-100%);
    transition: transform 140ms ease;
    font-size: 12px;
  }
  .props.open { transform: translateX(0); }
  @media (prefers-reduced-motion: reduce) {
    .props { transition: none; }
  }

  .props-head {
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin-bottom: 6px;
  }
  .props-title {
    font-size: 11px;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: var(--text-dim);
  }
  .props-close {
    display: flex;
    align-items: center;
    background: none;
    border: none;
    padding: 2px;
    cursor: pointer;
    color: var(--text-dim);
  }
  .props-close:hover { color: var(--text); }

  .props-sec { padding: 6px 0; border-top: 1px solid var(--border); }
  .props-sec h3 {
    margin: 0 0 5px;
    font-size: 11px;
    font-weight: 600;
    color: var(--text-dim);
  }

  .row {
    display: flex;
    justify-content: space-between;
    gap: 8px;
    padding: 2px 0;
    line-height: 1.5;
  }
  .k { color: var(--text-dim); flex: 0 0 auto; }
  /* The value can be a long title; it truncates rather than widening the panel. */
  .v {
    color: var(--text);
    text-align: right;
    min-width: 0;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
    font-variant-numeric: tabular-nums;
  }
  .v.muted { color: var(--text-muted); }
  .v.bad { color: var(--err); }

  .btn-row { display: flex; gap: 6px; margin: 4px 0; }
  .props-btn {
    display: inline-flex;
    align-items: center;
    gap: 4px;
    flex: 1 1 0;
    justify-content: center;
    padding: 4px 6px;
    font: inherit;
    font-size: 11px;
    color: var(--text);
    background: transparent;
    border: 1px solid var(--border);
    border-radius: 4px;
    cursor: pointer;
  }
  .props-btn:hover:not(:disabled) { border-color: var(--accent); }
  .props-btn:disabled { opacity: 0.45; cursor: default; }
  .props-btn.on { border-color: var(--accent); color: var(--accent); }

  .scale-row { align-items: center; }
  .scale-steps { display: flex; gap: 4px; }
  .scale-btn { flex: 0 0 auto; padding: 3px 5px; font-size: 10px; }

  .check {
    display: flex;
    align-items: center;
    gap: 6px;
    padding: 3px 0;
    cursor: pointer;
    color: var(--text);
  }
  .check input { margin: 0; }
</style>
