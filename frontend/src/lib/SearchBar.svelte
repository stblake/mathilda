<!--
  SearchBar.svelte — find across every cell of the focused notebook.

  Sits under the toolbar, opens on Cmd+F, and reports "n of m" over the WHOLE
  notebook rather than the focused cell. See search.ts for why this is not
  @codemirror/search.

  Navigation MARKS the current match in whichever cell owns it and scrolls it
  into view, and the find field keeps focus throughout, so Enter / Shift+Enter
  keep walking. It used to select the match and focus the cell's editor; the
  second Enter then went to the cell and replaced the match with a newline.
  Closing the bar (Escape, or the close button) is what finally puts the caret
  on the current match, as every editor's find does. It does not highlight every
  match at once; the count says how many there are and Enter walks them.
-->
<script lang="ts">
  import { tick } from 'svelte';
  import { onDestroy } from 'svelte';
  import Icon from './Icon.svelte';
  import { activeActions } from './canvas';
  import { getHandle } from './active';
  import { rangeForOffsets } from './searchHighlight';
  import type { Cell, NotebookRow } from './notebook';
  import { searchOpen, searchQuery, searchCaseSensitive, searchIndex,
           findMatches, nextIndex } from './search';
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

  /* The match currently painted, if any. Tracked so the old mark is cleared
     BEFORE the next one is painted (a prose mark is one document-wide highlight,
     so clearing after would erase the new one), and so closing knows where to put
     the caret. */
  let marked: SearchMatch | null = null;

  function clearMark() {
    if (marked) getHandle(marked.cellId)?.mark?.(null);
    marked = null;
  }

  /* A new query or case setting makes the painted match stale; the count updates
     as you type but nothing jumps until Enter. */
  $: $searchQuery, $searchCaseSensitive, clearMark();

  /* Leaving focused mode unmounts the bar; do not leave a mark behind. */
  onDestroy(clearMark);

  /** Close the bar. `placeCaret` puts the caret on the current match -- the
   *  handoff from finding to editing -- and is what Escape and the close button
   *  do. */
  function close(placeCaret = true) {
    const m = marked;
    clearMark();
    searchOpen.set(false);
    if (placeCaret && m) void selectInCell(m);
  }

  /* Closed from outside (Cmd+F toggling, a view switch): just drop the mark. */
  $: if (!$searchOpen && marked) clearMark();

  function go(delta: number) {
    if (!matches.length) return;
    /* The first Enter lands on the match the count already shows; later ones
       step from it. Without this the first press skipped match 1. */
    const next = nextIndex(current, delta, matches.length, marked !== null);
    searchIndex.set(next);
    reveal(matches[next]);
  }

  /** Paint the match and scroll it into view WITHOUT moving focus out of the
   *  find field. */
  function reveal(m: SearchMatch) {
    clearMark();
    const h = getHandle(m.cellId);
    if (!h?.mark) return;
    h.mark({ start: m.start, end: m.end });
    marked = m;
    inputEl?.focus();
  }

  /** Put the caret (a real selection) on the match in whichever editor owns it.
   *  Only on close: this moves focus into the cell. */
  async function selectInCell(m: SearchMatch) {
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
       flush. rangeForOffsets maps source offsets through the text nodes and <br>s
       the editor paints, and declines rather than guessing when they do not fit. */
    h.focus();
    await tick();
    const el = h.el;
    if (!el) return;
    el.scrollIntoView({ block: 'nearest' });
    const range = rangeForOffsets(el, m.start, m.end);
    if (!range) return;
    const sel = window.getSelection();
    sel?.removeAllRanges();
    sel?.addRange(range);
  }

  function onKeydown(e: KeyboardEvent) {
    if (e.key === 'Escape') { e.preventDefault(); close(true); return; }
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
            on:pointerdown|preventDefault on:click={() => close(true)}
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
    z-index: var(--z-overlay-top);
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
