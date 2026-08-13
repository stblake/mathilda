<!--
  ToolbarGroup.svelte — one labelled, vertically-ruled group of toolbar controls.

  Deliberately dumb: a caption, a row, and a rule. The rule is a border on the
  group itself rather than separate divider elements, so groups can be added and
  removed without anyone maintaining a parallel list of separators.

  `visible` exists because the context-sensitive Text/Code group has nothing to
  show when no cell is active, and a row of greyed-out letters under a "Text"
  caption reads as broken rather than as inapplicable.
-->
<script lang="ts">
  export let label: string;
  export let visible: boolean = true;
</script>

{#if visible}
  <div class="tb-group" role="group" aria-label={label}>
    <span class="tb-caption">{label}</span>
    <div class="tb-row"><slot /></div>
  </div>
{/if}

<style>
  .tb-group {
    display: flex;
    flex-direction: column;
    justify-content: center;
    gap: 2px;
    /* px, not rem: the bar is a fixed height and Cmd+= scales the root font
       size, so rem padding here would overflow the bar at 2x. */
    padding: 0 9px;
    border-inline-end: 1px solid var(--tb-rule);
    flex-shrink: 0;
  }
  /* The last group before the flexible spacer should not draw a trailing rule. */
  .tb-group:last-of-type { border-inline-end: none; }

  .tb-caption {
    font-size: var(--tb-caption-sz, 9px);
    line-height: 1;
    font-weight: 600;
    letter-spacing: 0.06em;
    text-transform: uppercase;
    color: var(--tb-caption);
    white-space: nowrap;
    user-select: none;
    -webkit-user-select: none;
  }

  .tb-row {
    display: flex;
    align-items: center;
    gap: 3px;
    height: var(--tb-btn-sz, 24px);
  }
</style>
