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
    flex-shrink: 0;
    /* A rule at the START of every group but the first, rather than at the end of
       every group but the last.
       Visually identical -- exactly one rule between any two groups -- but immune
       to what follows the groups in the DOM. The end-of-all-but-last form used
       :last-of-type, and Menu.svelte renders its backdrop and panel as sibling
       divs of the groups: with a menu open the last div was the menu, no group
       matched, and the trailing group grew a rule on every dropdown open. The
       first div child is always the first group, whatever is appended later. */
    border-inline-start: 1px solid var(--tb-rule);
  }
  .tb-group:first-of-type { border-inline-start: none; }

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
