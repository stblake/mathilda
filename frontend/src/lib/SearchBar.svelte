<!--
  SearchBar.svelte — find across every cell of the focused notebook.

  Sits under the toolbar, opens on Cmd+F, and reports "n of m" over the WHOLE
  notebook rather than the focused cell. See search.ts for why this is not
  @codemirror/search.

  Navigation drives whichever editor owns the match: a code cell gets a real
  CodeMirror selection and scroll, a prose cell is opened for editing and its
  range selected. What it does NOT do is highlight every match at once -- that
  needs a CodeMirror decoration extension per cell, which is its own change; the
  count tells you how many there are and Enter walks them.
-->
<script lang="ts">
  import { tick } from 'svelte';
  import { onDestroy } from 'svelte';
  import Icon from './Icon.svelte';
  import { activeActions } from './canvas';
  import { getHandle } from './active';
  import type { Cell, NotebookRow } from './notebook';
  import { searchOpen, searchQuery, searchCaseSensitive, searchIndex,
           findMatches, stepIndex } from './search';
  import type { SearchMatch } from './search';

  let inputEl: HTMLInputElement | undefined;

  /* The active pane's rows, subscribed by hand: the store identity changes with
     the pane, and `$store` only auto-subscribes a fixed identifier. */
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
  $: matches = findMatches(cells, $searchQuery, $searchCaseSensitive);
  /* Clamped rather than reset: editing the query usually shortens the list, and
     jumping back to the first match every keystroke would fight the typist. */
  $: current = matches.length ? Math.min($searchIndex, matches.length - 1) : 0;

  /* Focus the field when the bar opens. Not on every render -- that would steal
     focus back from the notebook after every jump. */
  let wasOpen = false;
  $: if ($searchOpen !== wasOpen) {
    wasOpen = $searchOpen;
    if ($searchOpen) void openFocus();
  }
  async function openFocus() {
    await tick();
    inputEl?.focus();
    inputEl?.select();
  }

  function close() {
    searchOpen.set(false);
  }

  function go(delta: number) {
    if (!matches.length) return;
    const next = stepIndex(current, delta, matches.length);
    searchIndex.set(next);
    void reveal(matches[next]);
  }

  /** Put the match on screen and select it in whichever editor owns it. */
  async function reveal(m: SearchMatch) {
    const h = getHandle(m.cellId);
    if (!h) return;
    if (h.view) {
      /* A code cell: CodeMirror does the selection and the scrolling, and the
         offsets are already source offsets, which is what it wants. */
      h.view.dispatch({
        selection: { anchor: m.start, head: m.end },
        scrollIntoView: true,
      });
      h.view.focus();
      return;
    }
    /* A prose cell. focus() opens the editor if the cell was showing rendered
       Markdown, and that is asynchronous, so the range is selected after the
       flush. Selecting is only attempted when the element holds exactly one text
       node, which is what the editor paints into it; anything else means the DOM
       is not what this assumes and leaving the caret alone beats guessing. */
    h.focus();
    await tick();
    const el = h.el;
    if (!el) return;
    el.scrollIntoView({ block: 'nearest' });
    const node = el.firstChild;
    if (!node || node.nodeType !== Node.TEXT_NODE || el.childNodes.length !== 1) return;
    const len = node.textContent?.length ?? 0;
    if (m.end > len) return;
    const range = document.createRange();
    range.setStart(node, m.start);
    range.setEnd(node, m.end);
    const sel = window.getSelection();
    sel?.removeAllRanges();
    sel?.addRange(range);
  }

  function onKeydown(e: KeyboardEvent) {
    if (e.key === 'Escape') { e.preventDefault(); close(); return; }
    if (e.key === 'Enter')  { e.preventDefault(); go(e.shiftKey ? -1 : 1); return; }
  }

  /* Typing a new query re-runs the search but does not jump: the count updates as
     you type, and Enter is what moves. Jumping per keystroke scrolls the notebook
     out from under someone who is still typing. */
</script>

{#if $searchOpen}
  <div class="search-bar" role="search">
    <Icon name="search" />
    <input
      class="search-input"
      type="text"
      placeholder="Find in notebook"
      bind:this={inputEl}
      bind:value={$searchQuery}
      on:keydown={onKeydown}
    />

    <span class="search-count" class:none={$searchQuery && !matches.length}>
      {#if !$searchQuery}
        &nbsp;
      {:else if matches.length}
        {current + 1} of {matches.length}
      {:else}
        No matches
      {/if}
    </span>

    <button
      class="search-btn"
      class:on={$searchCaseSensitive}
      title="Match case"
      aria-pressed={$searchCaseSensitive}
      on:pointerdown|preventDefault
      on:click={() => searchCaseSensitive.update(v => !v)}
    >Aa</button>

    <button class="search-btn" title="Previous match (Shift+Enter)" disabled={!matches.length}
            on:pointerdown|preventDefault on:click={() => go(-1)}
    ><Icon name="caretUp" /></button>

    <button class="search-btn" title="Next match (Enter)" disabled={!matches.length}
            on:pointerdown|preventDefault on:click={() => go(1)}
    ><Icon name="caret" /></button>

    <button class="search-btn" title="Close (Escape)"
            on:pointerdown|preventDefault on:click={close}
    ><Icon name="close" /></button>
  </div>
{/if}

<style>
  /* Under the toolbar, spanning the focused view. Absolute rather than in the
     flow: appearing must not resize the pane grid, which would relayout every
     editor in it. */
  .search-bar {
    position: absolute;
    top: var(--toolbar-h, 46px);
    right: 12px;
    z-index: 45;
    display: flex;
    align-items: center;
    gap: 6px;
    padding: 5px 8px;
    background: var(--menu-bg);
    border: 1px solid var(--menu-border);
    border-top: none;
    border-radius: 0 0 6px 6px;
    box-shadow: var(--menu-shadow);
    font-size: 12px;
  }

  .search-input {
    width: 200px;
    padding: 3px 5px;
    font: inherit;
    color: var(--text);
    background: var(--cell-bg);
    border: 1px solid var(--border);
    border-radius: 3px;
  }
  .search-input:focus { outline: none; border-color: var(--accent); }

  .search-count {
    min-width: 74px;
    text-align: right;
    color: var(--text-muted);
    font-variant-numeric: tabular-nums;
  }
  .search-count.none { color: var(--err); }

  .search-btn {
    display: flex;
    align-items: center;
    padding: 2px 4px;
    font: inherit;
    font-size: 11px;
    color: var(--text-dim);
    background: none;
    border: 1px solid transparent;
    border-radius: 3px;
    cursor: pointer;
  }
  .search-btn:hover:not(:disabled) { color: var(--text); border-color: var(--border); }
  .search-btn:disabled { opacity: 0.4; cursor: default; }
  .search-btn.on { color: var(--accent); border-color: var(--accent); }
</style>
