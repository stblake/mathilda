<!--
  App.svelte — thin shell.
  Renders the full-viewport Canvas. Handles Cmd+S / Cmd+O and
  kernel status / dark mode toggles in a floating corner overlay.
-->
<script lang="ts">
  import { onMount, onDestroy } from 'svelte';
  import { writable } from 'svelte/store';
  import { open, save } from '@tauri-apps/plugin-dialog';
  import { listen } from '@tauri-apps/api/event';
  import { getCurrentWindow } from '@tauri-apps/api/window';
  import Canvas from './lib/Canvas.svelte';
  import Toolbar from './lib/Toolbar.svelte';
  import { MENU_IDS, runMenuCommand } from './lib/menuCommands';
  import { kernelStatus } from './lib/notebook';
  import { darkMode } from './lib/theme';
  import { kernelMemory } from './lib/status';
  import { pingKernel, saveLibrary, loadLibrary, setWindowTitle as setTitleCmd } from './lib/ipc';
  import { restart, abortEvaluation } from './lib/kernelActions';
  import { serializeLibrary, loadLibraryData, canvasState, activeActions, activeFlags, setFocused } from './lib/canvas';
  /* Imported for its side effect: installs the document-level Cmd+click
     handler that opens a symbol's reference page. Importing it here rather
     than relying on a cell to pull it in means the gesture works from the
     moment the app loads. */
  import './lib/refpages';

  // ---------------------------------------------------------------------------
  // Dark mode. The store lives in lib/theme.ts so the toolbar and the
  // properties panel can reach it; the DOM write stays here.

  /* Set BOTH classes, not just .light: app.css keys its palette off
     :root:not(.light) inside the dark media query and :root.dark outside it,
     so an explicit choice has to be stated in whichever direction it differs
     from the OS. Toggling only .light left OS-light + app-dark unstyled. */
  $: if (typeof document !== 'undefined') {
    document.documentElement.classList.toggle('light', !$darkMode);
    document.documentElement.classList.toggle('dark', $darkMode);
  }

  /* Title of the ACTIVE pane, for the app bar. With one pane that is the
     notebook filling the window, as before; with several it names the one the
     bar's own buttons will act on, which is the cheapest and strongest signal
     available that the toolbar has a target. */
  $: focusedTitle = $canvasState.focusedActiveId
    ? ($canvasState.notebooks.find(n => n.id === $canvasState.focusedActiveId)?.title ?? '')
    : '';

  // ---------------------------------------------------------------------------
  // Kernel init

  onMount(async () => {
    kernelStatus.set('starting');
    await new Promise(r => setTimeout(r, 1200));
    try { await pingKernel(); kernelStatus.set('ready'); }
    catch  { kernelStatus.set('dead'); }
  });

  // ---------------------------------------------------------------------------
  // Menu event listeners

  let unlisten: (() => void)[] = [];
  let libraryTitle = 'Untitled Library';
  let libraryPath: string | null = null;

  onMount(async () => {
    try {
      /* Every id the native menu can emit, subscribed from one list. Previously each item needed
         its own one-line listener here, so an item added in Rust silently did nothing until
         someone remembered this file; now the two sides share MENU_IDS and the dispatcher warns
         about an id it has no case for. */
      /* The kernel reports its resident memory with every `done`; the status bar shows the
         latest. Its own event rather than part of a cell's output stream, since a memory reading
         is not output. */
      unlisten.push(await listen<number>('kernel-memory',
                                         (e) => kernelMemory.set(e.payload)));

      const hooks = { openFile, saveFile, saveFileAs };
      for (const id of MENU_IDS) {
        unlisten.push(await listen(`menu:${id}`, () => runMenuCommand(id, hooks)));
      }
    } catch (e) { console.warn('Menu listen error:', e); }
  });

  onDestroy(() => unlisten.forEach(u => u()));

  // ---------------------------------------------------------------------------
  // File I/O — library-level (whole canvas)

  async function openFile() {
    const sel = await open({
      filters: [{ name: 'Mathilda Library', extensions: ['lb', 'mathilda'] }],
    });
    if (!sel) return;
    const path = typeof sel === 'string' ? sel : (sel as string[])[0];
    try {
      // Use loadLibrary (returns raw JSON string) not loadNotebook (parses as cells)
      const json = await loadLibrary(path);
      const title = loadLibraryData(json);
      libraryTitle = title;
      libraryPath  = path;
      const filename = path.split('/').pop()?.replace(/\.lb$/i, '') ?? title;
      setWindowTitle(filename);
    } catch (e) { console.error('Open failed:', e); }
  }

  async function saveFile() {
    if (libraryPath) doSave(libraryPath); else saveFileAs();
  }

  async function saveFileAs() {
    const path = await save({
      defaultPath: (libraryTitle || 'library') + '.lb',
      filters: [{ name: 'Mathilda Library', extensions: ['lb'] }],
    });
    if (!path) return;
    libraryPath = path;
    doSave(path);
  }

  async function doSave(path: string) {
    try {
      const json = serializeLibrary(libraryTitle);
      await saveLibrary(path, json);
      const filename = path.split('/').pop()?.replace(/\.lb$/i, '') ?? 'Library';
      libraryTitle = filename;
      setWindowTitle(filename);
    } catch (e) { console.error('Save failed:', e); }
  }

  function setWindowTitle(name: string) {
    const title = `Mathilda — ${name}`;
    document.title = title;
    // Use Rust command — most reliable way to set native macOS title bar
    setTitleCmd(title).catch(() => {});
  }

  // ---------------------------------------------------------------------------
  // Kernel restart

  /* restart() and abortEvaluation() live in lib/kernelActions.ts so the
     toolbar's kernel menu drives exactly the same paths as the native Kernel
     menu. Two implementations of "abort" is how one of them ends up wrong. */

  // ---------------------------------------------------------------------------
  // Global UI scale (Cmd+= zoom in, Cmd+- zoom out, Cmd+0 reset)
  let uiScale = 1.0;
  $: document.documentElement.style.fontSize = `${uiScale * 16}px`;

  function onKeydown(e: KeyboardEvent) {
    const mod = e.metaKey || e.ctrlKey;
    if (!mod) return;
    if (e.key === 's' || e.key === 'S') { e.preventDefault(); saveFile(); return; }
    if (e.key === 'o' || e.key === 'O') { e.preventDefault(); openFile(); return; }
    // Cmd+= / Cmd++ → scale up; Cmd+- → scale down; Cmd+0 → reset
    if (e.key === '=' || e.key === '+') {
      e.preventDefault(); uiScale = Math.min(2.0, +(uiScale + 0.1).toFixed(1));
    } else if (e.key === '-' || e.key === '_') {
      e.preventDefault(); uiScale = Math.max(0.5, +(uiScale - 0.1).toFixed(1));
    }
  }
</script>

<svelte:window on:keydown={onKeydown} />

<!-- Full-viewport canvas -->
<Canvas />

<!-- App bar. Two quite different things share this strip:

     On the canvas it is a 34px name-and-theme strip.

     In focused mode it becomes the notebook toolbar — labelled, ruled groups at
     46px. Two heights rather than one reactive variable, because making
     --appbar-h itself change would put a JS style write in the middle of the
     {#if} branch swap between .canvas-stage and .focused-view, and would drag
     .canvas-stage into a change it has no stake in. -->
<div class="app-bar" class:toolbar-mode={$canvasState.focusedIds.length > 0}>
  {#if $canvasState.focusedIds.length}
    <!-- One row again. The menus are NATIVE now -- on macOS they live in the system bar at the
         top of the screen -- so this strip is free for controls that belong to the window. -->
    <Toolbar />
  {:else}
    <span class="app-bar-name">Mathilda</span>
    <span class="dark-toggle-spacer"></span>
    <button
      class="dark-toggle"
      title="Toggle dark mode"
      on:click={() => darkMode.update(v => !v)}
    >
      {$darkMode ? '◑' : '☀'}
    </button>
  {/if}
</div>

<!-- Kernel dead banner -->
{#if $kernelStatus === 'dead'}
  <div class="kernel-banner">
    Kernel not running.
    <button on:click={restart}>Restart</button>
  </div>
{/if}

<style>
  /* ---- Dark mode (default — :root always applies) ---- */
  :global(:root) {
    --bg:          #050810;
    --surface:     rgba(8,10,22,0.96);
    --cell-bg:     rgba(12,15,28,0.85);
    --border:      rgba(255,255,255,0.06);
    --text:        #cdd6f4;
    --text-muted:  #45475a;
    --accent:      #89b4fa;
    --accent-glow: rgba(137,180,250,0.10);
    --out-text:    #cdd6f4;
    --gutter-bg:   rgba(255,255,255,0.015);
    --gutter-hover:rgba(255,255,255,0.03);
    --card-bg:     rgba(12,15,28,0.85);
    --card-border: rgba(255,255,255,0.08);

    /* Toolbar / menu surfaces.
       --surface-2 was referenced at .app-bar .tb-btn:hover with NO fallback and
       defined nowhere in the tree, so the app bar's hover state was a literal
       no-op. RefPage.svelte referenced it in five more places. Defining it here
       rather than in app.css because App.svelte's :global(:root) owns the
       surface palette and wins by load order anyway. */
    --surface-2:   rgba(255,255,255,0.055);   /* hover / raised fill */
    --surface-3:   rgba(255,255,255,0.10);    /* pressed / active toggle */
    --tb-rule:     rgba(255,255,255,0.09);    /* group divider, subtler than --border */
    --tb-caption:  #6c7086;                   /* group caption, dimmer than --text-dim */
    --menu-bg:     rgba(18,21,34,0.98);
    --menu-border: rgba(255,255,255,0.10);
    --menu-shadow: 0 10px 32px rgba(0,0,0,0.55);
    /* --ok was only ever used as a var() fallback; --err was hardcoded. */
    --ok:          #4ade80;
    --warn:        #fab387;
    --err:         #f38ba8;
  }
  :global(body) { background: #050810; }

  /* ---- Light mode (html.light class applied when darkMode = false) ---- */
  :global(html.light) {
    --bg:          #e8e9f0;  /* lighter canvas — less contrast with white cards */
    --surface:     #f8f8fc;
    --cell-bg:     #f8f8fc;
    --border:      rgba(0,0,0,0.06);  /* much softer cell dividers */
    --text:        #1c1c2e;
    --text-muted:  #666688;
    --accent:      #3b82f6;
    --accent-glow: rgba(59,130,246,0.15);
    --out-text:    #1c1c2e;
    --gutter-bg:   #eeeef5;
    --gutter-hover:#e4e5f0;
    --card-bg:     #f8f8fc;
    --card-border: rgba(0,0,0,0.08);

    --surface-2:   rgba(0,0,0,0.045);
    --surface-3:   rgba(0,0,0,0.085);
    --tb-rule:     rgba(0,0,0,0.10);
    --tb-caption:  #8a8a9e;
    --menu-bg:     #ffffff;
    --menu-border: rgba(0,0,0,0.12);
    --menu-shadow: 0 10px 32px rgba(0,0,0,0.18);
    --ok:          #16a34a;
    --warn:        #d97706;
    --err:         #dc2626;
  }
  :global(html.light body) { background: #1a1b2e; }

  :global(*, *::before, *::after) { box-sizing: border-box; }
  :global(body) {
    margin: 0;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    color: #cdd6f4;
    overflow: hidden;
  }

  /* ---- Corner overlay ---- */
  .app-bar {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    height: var(--appbar-h, 34px);
    display: flex;
    align-items: center;
    gap: 0.75rem;
    padding: 0 0.75rem;
    background: var(--bg);
    border-bottom: 1px solid var(--border);
    z-index: 200;
    /* The toolbar is a row of fixed-height groups; nothing here may wrap. */
    overflow: hidden;
    /* Nothing here should swallow a drag meant for the window chrome. */
    user-select: none;
    -webkit-user-select: none;
  }

  /* Focused mode: taller, and its own padding in px rather than rem. The bar is
     a fixed height full of fixed-px content, so rem padding would push the
     groups out of it once Cmd+= scales the root font size. */
  .app-bar.toolbar-mode {
    height: var(--toolbar-h, 46px);
    gap: 0;
    padding: 0 8px;
    align-items: stretch;
  }

  .app-bar-name {
    font: 600 0.78rem/1 var(--sans);
    color: var(--text-h);
    letter-spacing: 0.01em;
  }

  /* The centred title and the seven glyph buttons that used to live here now
     belong to Toolbar.svelte, which owns its own styles. A centred overlay title
     cannot coexist with a full-width row of groups. */

  .dark-toggle-spacer { flex: 1; }

  .dark-toggle {
    background: rgba(128,128,128,0.1);
    border: 1px solid rgba(128,128,128,0.2);
    color: var(--text-muted, #585b70);
    cursor: pointer;
    font-size: 0.9rem;
    padding: 4px 8px;
    border-radius: 6px;
    line-height: 1;
    transition: color 0.1s, background 0.1s;
    min-width: 32px;
    text-align: center;
  }
  .dark-toggle:hover { color: var(--text, #cdd6f4); background: rgba(128,128,128,0.2); }

  /* ---- Kernel dead banner ---- */
  .kernel-banner {
    position: fixed;
    bottom: 1.2rem;
    left: 50%;
    transform: translateX(-50%);
    background: #c0392b;
    color: white;
    padding: 0.45rem 1rem;
    border-radius: 6px;
    font-size: 0.84rem;
    display: flex;
    align-items: center;
    gap: 0.6rem;
    box-shadow: 0 4px 12px rgba(0,0,0,0.5);
    z-index: 300;
  }
  .kernel-banner button {
    background: rgba(255,255,255,0.2);
    border: 1px solid rgba(255,255,255,0.35);
    color: white;
    border-radius: 3px;
    padding: 0.15rem 0.6rem;
    cursor: pointer;
    font-size: 0.8rem;
  }
</style>
