<!--
  Menu.svelte — the one popover primitive.

  There is no Dropdown.svelte wrapper on purpose. The evaluation caret, the cell
  style combo and the … overflow have three visually distinct triggers with no
  shared behaviour, so a wrapper would be pure props passthrough. The trigger
  stays a plain <button bind:this={anchor}> in the caller.

  Toolbar.svelte renders ONE instance of this and swaps `items`, so there is one
  backdrop, one keyboard handler and one measure path rather than four.

  position:fixed, so an ancestor's overflow cannot clip it. Note that fixed does
  NOT escape an ancestor that creates a containing block via transform/filter/
  backdrop-filter -- and .nb-card has backdrop-filter: blur(20px). That is why
  this must be rendered from the toolbar and never from inside a card.

  Keyboard focus goes on the CONTAINER, with aria-activedescendant, rather than
  calling focus() per item -- that avoids the roving-tabindex "who actually has
  focus" mess in a four-item menu.
-->
<!-- A module script, not the instance script: `export` in an instance script
     declares a PROP, so a type declared there would not be importable by the
     callers that need to annotate their item arrays. -->
<script context="module" lang="ts">
  export type MenuItem =
    | { kind: 'sep' }
    | {
        kind: 'item';
        id: string;
        label: string;
        icon?: string;
        /** Right-aligned shortcut text, e.g. "⇧↵". */
        hint?: string;
        disabled?: boolean;
        /** Shows a tick — used by the cell-style combo for the current type. */
        checked?: boolean;
      };
</script>

<script lang="ts">
  import { createEventDispatcher, tick } from 'svelte';
  import Icon from './Icon.svelte';
  import { isTouchDevice } from './platform';

  export let open: boolean;
  export let items: MenuItem[] = [];
  export let anchor: HTMLElement | null = null;
  /** Which edge of the anchor the menu lines up with. */
  export let align: 'start' | 'end' = 'start';
  export let minWidth: number = 0;

  const dispatch = createEventDispatcher<{ select: { id: string }; close: void }>();

  let menuEl: HTMLElement;
  let left = 0;
  let top = 0;
  let activeIdx = -1;

  const GAP = 4;
  const EDGE = 8;

  function firstEnabled(from: number, dir: 1 | -1): number {
    const n = items.length;
    if (n === 0) return -1;
    for (let step = 0; step < n; step++) {
      const i = ((from + dir * step) % n + n) % n;
      const it = items[i];
      if (it.kind === 'item' && !it.disabled) return i;
    }
    return -1;
  }

  /** Measure against the anchor and clamp to the viewport, so a right-hand
   *  group's menu (the … overflow, in particular) does not run off-screen. */
  function place() {
    if (!anchor || !menuEl) return;
    const a = anchor.getBoundingClientRect();
    const m = menuEl.getBoundingClientRect();
    const want = align === 'end' ? a.right - m.width : a.left;
    left = Math.max(EDGE, Math.min(want, window.innerWidth - m.width - EDGE));
    top = Math.min(a.bottom + GAP, window.innerHeight - m.height - EDGE);
  }

  /* Open: measure, focus the container, and pre-arm the first enabled item so a
     keyboard user can press Enter immediately. */
  $: if (open) {
    tick().then(() => {
      place();
      menuEl?.focus();
      if (activeIdx < 0) activeIdx = firstEnabled(0, 1);
    });
  } else {
    activeIdx = -1;
  }

  function choose(i: number) {
    const it = items[i];
    if (!it || it.kind !== 'item' || it.disabled) return;
    dispatch('select', { id: it.id });
    dispatch('close');
  }

  function onKeydown(e: KeyboardEvent) {
    switch (e.key) {
      case 'Escape':
        e.preventDefault();
        e.stopPropagation();
        dispatch('close');
        anchor?.focus();           /* return focus where it came from */
        break;
      case 'ArrowDown':
        e.preventDefault();
        activeIdx = firstEnabled(activeIdx < 0 ? 0 : activeIdx + 1, 1);
        break;
      case 'ArrowUp':
        e.preventDefault();
        activeIdx = firstEnabled(activeIdx < 0 ? items.length - 1 : activeIdx - 1, -1);
        break;
      case 'Home':
        e.preventDefault();
        activeIdx = firstEnabled(0, 1);
        break;
      case 'End':
        e.preventDefault();
        activeIdx = firstEnabled(items.length - 1, -1);
        break;
      case 'Enter':
      case ' ':
        e.preventDefault();
        choose(activeIdx);
        break;
      case 'Tab':
        /* Don't trap focus in a four-item menu. */
        dispatch('close');
        break;
    }
  }
</script>

<svelte:window on:resize={() => open && place()} />

{#if open}
  <!-- Click-outside catcher. `click`, not `pointerdown`: pointerdown here would
       swallow the first click on whatever the user was actually aiming at. -->
  <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
  <div class="menu-backdrop" on:click={() => dispatch('close')}></div>

  <div
    class="menu"
    class:touch={isTouchDevice}
    role="menu"
    tabindex="-1"
    aria-activedescendant={activeIdx >= 0 ? `menu-item-${activeIdx}` : undefined}
    bind:this={menuEl}
    on:keydown={onKeydown}
    style="left:{left}px; top:{top}px; {minWidth ? `min-width:${minWidth}px;` : ''}"
  >
    {#each items as it, i}
      {#if it.kind === 'sep'}
        <div class="menu-sep" role="separator"></div>
      {:else}
        <!-- svelte-ignore a11y-click-events-have-key-events -->
        <div
          id="menu-item-{i}"
          class="menu-item"
          class:active={i === activeIdx}
          class:disabled={it.disabled}
          role="menuitem"
          aria-disabled={it.disabled ? 'true' : undefined}
          on:click={() => choose(i)}
          on:pointerenter={() => { if (!it.disabled) activeIdx = i; }}
        >
          <span class="menu-tick">{it.checked ? '✓' : ''}</span>
          {#if it.icon}<span class="menu-icon"><Icon name={it.icon} size={15} /></span>{/if}
          <span class="menu-label">{it.label}</span>
          {#if it.hint}<span class="menu-hint">{it.hint}</span>{/if}
        </div>
      {/if}
    {/each}
  </div>
{/if}

<style>
  /* Above .app-bar (200) and .kernel-banner (300). */
  .menu-backdrop { position: fixed; inset: 0; z-index: 399; }

  .menu {
    position: fixed;
    z-index: 400;
    min-width: 180px;
    padding: 4px;
    border-radius: 8px;
    background: var(--menu-bg);
    border: 1px solid var(--menu-border);
    box-shadow: var(--menu-shadow);
    outline: none;
    font: 500 0.78rem/1 var(--sans);
    color: var(--text);
    max-height: calc(100vh - 80px);
    overflow-y: auto;
  }

  .menu-item {
    display: flex;
    align-items: center;
    gap: 7px;
    height: 26px;
    padding: 0 8px 0 4px;
    border-radius: 5px;
    cursor: pointer;
    white-space: nowrap;
  }
  .menu-item.active:not(.disabled) { background: var(--surface-2); color: var(--text-h); }
  .menu-item.disabled { opacity: 0.38; cursor: default; }

  /* Fixed-width so labels line up whether or not an item is ticked. */
  .menu-tick {
    width: 12px;
    flex-shrink: 0;
    text-align: center;
    color: var(--accent);
    font-size: 0.72rem;
  }
  .menu-icon  { display: flex; opacity: 0.8; }
  .menu-label { flex: 1; }
  .menu-hint  { opacity: 0.5; font-family: var(--mono); font-size: 0.7rem; padding-left: 12px; }

  .menu-sep { height: 1px; margin: 4px 6px; background: var(--menu-border); }

  /* Finger targets, and a wider box so labels aren't cramped. */
  .menu.touch { min-width: 220px; font-size: 0.9rem; }
  .menu.touch .menu-item { height: 44px; padding: 0 12px 0 6px; }
</style>
