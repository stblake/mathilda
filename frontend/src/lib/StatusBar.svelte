<!--
  StatusBar.svelte — the optional strip along the bottom of the focused view.

  Reports what the window otherwise cannot tell you: whether the kernel is alive
  and busy, how long the last evaluation took, and what the session has cost so
  far. Evaluation timing was not measured anywhere before this bar existed; see
  status.ts.

  Everything here is read-only and derived. The one exception is the kernel
  segment, which is clickable because "the kernel is dead" is the one status a
  user needs to act on immediately, and making them hunt through a menu for
  Restart at that moment is the wrong answer.
-->
<script lang="ts">
  import { kernelStatus } from './notebook';
  import { lastOp, evalCount, evalTotalMs, formatMs, kernelMemory, formatBytes } from './status';
  import { activeCell } from './active';
  import { activeActions } from './canvas';
  import { restart } from './kernelActions';
  import type { Cell, NotebookRow } from './notebook';
  import { onDestroy } from 'svelte';

  const KERNEL_LABEL: Record<string, string> = {
    starting:   'Starting…',
    ready:      'Ready',
    busy:       'Running',
    restarting: 'Restarting…',
    dead:       'Not running',
  };

  /* Cell counts come from the active pane's store, which changes with the pane,
     so it is subscribed by hand -- $store auto-subscription only works on a
     fixed identifier. */
  let rows: NotebookRow[] = [];
  let unsub: (() => void) | null = null;
  $: {
    unsub?.();
    unsub = null;
    const store = $activeActions?.store;
    if (store) unsub = store.subscribe((r: NotebookRow[]) => { rows = r; });
  }
  onDestroy(() => unsub?.());

  $: cells = rows.flatMap(r => r.cells) as Cell[];
  $: codeCells = cells.filter(c => c.type === 'code').length;
  $: evaluated = cells.filter(c => c.execIdx != null).length;

  /* The active cell, resolved from the store rather than trusted from the
     record -- same reason as the toolbar: a deleted cell must read as absent. */
  $: activeInfo = ($activeCell && $activeActions && $activeCell.notebookId === $activeActions.notebookId)
    ? cells.find(c => c.id === $activeCell!.cellId) ?? null
    : null;

  /* One line, and short: this sits in a 22px strip. */
  $: opSource = $lastOp
    ? (() => {
        const first = $lastOp.source.split('\n')[0].trim();
        return first.length > 42 ? first.slice(0, 42) + '…' : first;
      })()
    : '';
</script>

<div class="status-bar" role="status" aria-live="polite">
  <!-- Kernel. Clickable only when acting on it is the obvious next step. -->
  <button
    class="seg seg-kernel"
    data-status={$kernelStatus}
    class:actionable={$kernelStatus === 'dead'}
    disabled={$kernelStatus !== 'dead'}
    title={$kernelStatus === 'dead' ? 'Kernel is not running — click to restart' : `Local kernel — ${KERNEL_LABEL[$kernelStatus] ?? $kernelStatus}`}
    tabindex="-1"
    on:pointerdown|preventDefault
    on:click={() => $kernelStatus === 'dead' && restart()}
  >
    <span class="dot"></span>
    <span class="seg-label">Local</span>
    <span class="seg-value">{KERNEL_LABEL[$kernelStatus] ?? $kernelStatus}</span>
  </button>

  <span class="rule"></span>

  <!-- Last operation -->
  {#if $lastOp}
    <span class="seg" class:err={!$lastOp.ok} title={$lastOp.source}>
      <span class="seg-label">{$lastOp.label}</span>
      <span class="seg-value strong">{formatMs($lastOp.ms)}</span>
      {#if opSource}<span class="seg-src">{opSource}</span>{/if}
    </span>
  {:else}
    <span class="seg muted"><span class="seg-value">No evaluations yet</span></span>
  {/if}

  <span class="spacer"></span>

  <!-- Session totals -->
  {#if $evalCount > 0}
    <!-- The kernel's resident memory, as of its last evaluation. Absent until the first one:
         the number can only change because something was evaluated, so there is nothing honest
         to show before that. -->
    {#if $kernelMemory != null}
      <span class="seg" title="Resident memory of the kernel process, after its last evaluation">
        <span class="seg-label">Memory</span>
        <span class="seg-value">{formatBytes($kernelMemory)}</span>
      </span>
      <span class="rule"></span>
    {/if}
    <span class="seg" title="Evaluations completed this session, and their total kernel time">
      <span class="seg-label">Session</span>
      <span class="seg-value">{$evalCount} eval{$evalCount === 1 ? '' : 's'} · {formatMs($evalTotalMs)}</span>
    </span>
    <span class="rule"></span>
  {/if}

  <!-- The active cell, so the toolbar's target is always visible -->
  {#if activeInfo}
    <span class="seg" title="The cell the toolbar acts on">
      <span class="seg-label">Cell</span>
      <span class="seg-value">
        {activeInfo.execIdx != null ? `In[${activeInfo.execIdx}] ` : ''}{activeInfo.type}
      </span>
    </span>
    <span class="rule"></span>
  {/if}

  <!-- Notebook shape -->
  <span class="seg" title="Cells in this notebook, code cells, and how many have been evaluated">
    <span class="seg-value">
      {cells.length} cell{cells.length === 1 ? '' : 's'} · {codeCells} code · {evaluated} run
    </span>
  </span>
</div>

<style>
  .status-bar {
    display: flex;
    align-items: center;
    gap: 8px;
    height: var(--statusbar-h, 22px);
    flex-shrink: 0;
    padding: 0 10px;
    background: var(--gutter-bg);
    border-top: 1px solid var(--border);
    /* px, not rem: a fixed-height strip must not grow with Cmd+= zoom. */
    font: 500 11px/1 var(--sans);
    color: var(--text);
    user-select: none;
    -webkit-user-select: none;
    overflow: hidden;
    white-space: nowrap;
  }

  .seg {
    display: flex;
    align-items: center;
    gap: 5px;
    min-width: 0;
    background: none;
    border: none;
    padding: 0;
    font: inherit;
    color: inherit;
  }
  .seg.muted { opacity: 0.5; }
  .seg.err .seg-value.strong { color: var(--err); }

  .seg-label {
    font-size: 9px;
    font-weight: 600;
    letter-spacing: 0.06em;
    text-transform: uppercase;
    color: var(--tb-caption);
  }
  .seg-value { font-variant-numeric: tabular-nums; }
  .seg-value.strong { color: var(--text-h); font-weight: 600; }

  .seg-src {
    font-family: var(--mono);
    font-size: 10px;
    opacity: 0.55;
    overflow: hidden;
    text-overflow: ellipsis;
    min-width: 0;
  }

  .rule { width: 1px; height: 11px; background: var(--tb-rule); flex-shrink: 0; }
  .spacer { flex: 1; min-width: 6px; }

  /* ---- Kernel segment ---- */
  .seg-kernel { cursor: default; }
  .seg-kernel.actionable { cursor: pointer; }
  .seg-kernel.actionable:hover .seg-value { text-decoration: underline; }

  .dot {
    width: 7px;
    height: 7px;
    border-radius: 50%;
    flex-shrink: 0;
    background: var(--text-muted);
  }
  [data-status='ready']      .dot { background: var(--ok); }
  [data-status='busy']       .dot { background: var(--warn); animation: blink 1s ease-in-out infinite; }
  [data-status='starting']   .dot { background: var(--accent); animation: blink 1s ease-in-out infinite; }
  [data-status='restarting'] .dot { background: var(--warn); animation: blink 1s ease-in-out infinite; }
  [data-status='dead']       .dot { background: var(--err); }

  [data-status='dead'] .seg-value { color: var(--err); }

  @keyframes blink { 0%, 100% { opacity: 1; } 50% { opacity: 0.3; } }
</style>
