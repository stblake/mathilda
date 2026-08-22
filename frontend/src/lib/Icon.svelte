<!--
  Icon.svelte — a 24x24 inline SVG from a local path map.

  Inline SVG rather than Unicode glyphs, for two reasons that only bite at
  toolbar scale. Unicode has no acceptable glyph for split-cell, merge, or
  comment-toggle. And glyph metrics vary wildly across the four WebViews this app
  ships to (macOS/iOS/Android/Windows), which is exactly why the old app bar's
  `▶▶ ⇄ ⌃ ✎ ⤡ ⊞ ✕` row sits at inconsistent heights -- tolerable in a 34px bar of
  seven buttons, not next to a caption baseline and a vertical rule.

  Inline rather than a <symbol> sprite: a sprite needs a document-level <svg>
  host injected before first paint, and <use> shadow content does not inherit
  Svelte's scoped CSS. That is real machinery for a map whose path data is under
  2KB.

  stroke="currentColor" is the whole payoff: light/dark, hover, pressed and
  disabled all fall out of the button's own `color` with no per-icon work.

  Typographic marks are deliberately NOT here -- B, I, U, f[], (*=*), … and the
  ◑/☀ theme glyph stay as text, because they should read as letters. Wolfram
  mixes the same way.
-->
<script lang="ts">
  export let name: string;
  export let size: number = 16;

  /* Stroke-only outlines on a 24x24 grid. Add entries as controls land. */
  const PATHS: Record<string, string> = {
    /* navigation */
    back:       'M15 6l-6 6 6 6',
    caret:      'M6 9.5l6 6 6-6',
    caretUp:    'M6 14.5l6-6 6 6',
    close:      'M6 6l12 12M18 6L6 18',
    plus:       'M12 5v14M5 12h14',

    /* evaluation */
    restart:    'M20 12a8 8 0 11-2.34-5.66M20 4v4h-4',
    abort:      'M7 7h10v10H7z',

    /* panels */
    sidebar:    'M3 4h18v16H3zM9 4v16',
    gear:       'M12 15.2a3.2 3.2 0 100-6.4 3.2 3.2 0 000 6.4M19.6 12c0 .5-.05 1-.14 1.46l1.83 1.4-1.7 2.95-2.16-.83c-.74.6-1.6 1.05-2.53 1.31L14.5 21.6h-3.4l-.4-2.31a7.7 7.7 0 01-2.53-1.31l-2.16.83-1.7-2.95 1.83-1.4a7.9 7.9 0 010-2.92l-1.83-1.4 1.7-2.95 2.16.83A7.7 7.7 0 0110.7 6.7l.4-2.3h3.4l.4 2.3c.93.26 1.79.71 2.53 1.31l2.16-.83 1.7 2.95-1.83 1.4c.09.46.14.95.14 1.47z',

    /* cell structure */
    split:      'M4 3.5h16v6H4zM4 14.5h16v6H4zM9 12h6M13 10l2 2-2 2',
    merge:      'M4 3.5h16v5H4zM4 15.5h16v5H4zM12 10v4M10 12.2l2 1.8 2-1.8',
    duplicate:  'M9 9h11v11H9zM15 5H4v11',
    trash:      'M4 7h16M9 7V4h6v3M6 7l1 13h10l1-13M10 11v6M14 11v6',

    /* code editing */
    indent:     'M4 6h16M10 12h10M4 18h16M4 10l2.5 2L4 14',
    outdent:    'M4 6h16M10 12h10M4 18h16M6.5 10L4 12l2.5 2',

    /* layout */
    layoutH:    'M3 4h18v16H3zM12 4v16',
    layoutV:    'M3 4h18v16H3zM3 12h18',
    layoutGrid: 'M3 4h18v16H3zM12 4v16M3 12h18',

    /* notebook */
    search:     'M10.5 17a6.5 6.5 0 100-13 6.5 6.5 0 000 13zM15.4 15.4L20 20',
    docs:       'M12 8.5v.01M12 11.5v4.5M12 21a9 9 0 100-18 9 9 0 000 18z',
    link:       'M10 13.5a4 4 0 005.66 0l2.5-2.5a4 4 0 10-5.66-5.66L11.6 6.7M14 10.5a4 4 0 00-5.66 0l-2.5 2.5a4 4 0 105.66 5.66l.9-.9',
    rename:     'M4 20h4L19 9a2.83 2.83 0 10-4-4L4 16v4z',
    collapse:   'M4 5h16v14H4zM4 10h16',

    /* insert */
    sigma:      'M17 5H7l5 7-5 7h10',
    integral:   'M9 19c0-8 0-14 3-14M13.5 5c1.5 0 2 1 2 2',
    tex:        'M4 7V5h7v2M7.5 5v14M13 10h7M16.5 10v9',
  };

  /* Solid glyphs: an outlined play triangle reads as a "next" chevron. */
  const FILLED = new Set(['run', 'runAll']);
  const FILL_PATHS: Record<string, string> = {
    run:    'M8 5.2l11 6.8-11 6.8z',
    runAll: 'M4 5.2l7.5 6.8L4 18.8zM12.5 5.2L20 12l-7.5 6.8z',
  };

  $: filled = FILLED.has(name);
  $: d = filled ? (FILL_PATHS[name] ?? '') : (PATHS[name] ?? '');
</script>

<svg
  width={size}
  height={size}
  viewBox="0 0 24 24"
  fill={filled ? 'currentColor' : 'none'}
  stroke={filled ? 'none' : 'currentColor'}
  stroke-width="1.6"
  stroke-linecap="round"
  stroke-linejoin="round"
  aria-hidden="true"
  focusable="false"
><path {d} /></svg>

<style>
  svg { display: block; flex-shrink: 0; }
</style>
