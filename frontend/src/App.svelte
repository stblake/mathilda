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
  import KernelStatus from './lib/KernelStatus.svelte';
  import { kernelStatus } from './lib/notebook';
  import { pingKernel, restartKernel, saveLibrary, loadLibrary, setWindowTitle as setTitleCmd } from './lib/ipc';
  import { serializeLibrary, loadLibraryData, canvasState, focusedActions, setFocused } from './lib/canvas';
  /* Imported for its side effect: installs the document-level Cmd+click
     handler that opens a symbol's reference page. Importing it here rather
     than relying on a cell to pull it in means the gesture works from the
     moment the app loads. */
  import './lib/refpages';

  // ---------------------------------------------------------------------------
  // Dark mode (default to dark — canvas is always dark)

  const darkMode = writable(true);
  /* Set BOTH classes, not just .light: app.css keys its palette off
     :root:not(.light) inside the dark media query and :root.dark outside it,
     so an explicit choice has to be stated in whichever direction it differs
     from the OS. Toggling only .light left OS-light + app-dark unstyled. */
  $: if (typeof document !== 'undefined') {
    document.documentElement.classList.toggle('light', !$darkMode);
    document.documentElement.classList.toggle('dark', $darkMode);
  }

  /* Title of the notebook currently filling the window, for the app bar. */
  $: focusedTitle = $canvasState.focusedId
    ? ($canvasState.notebooks.find(n => n.id === $canvasState.focusedId)?.title ?? '')
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
      unlisten.push(await listen('menu:open',        () => openFile()));
      unlisten.push(await listen('menu:save',        () => saveFile()));
      unlisten.push(await listen('menu:save-as',     () => saveFileAs()));
      unlisten.push(await listen('menu:restart',     () => restart()));
      unlisten.push(await listen('menu:interrupt',   () => {}));
      unlisten.push(await listen('menu:toggle-dark', () => darkMode.update(v => !v)));
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

  async function restart() {
    kernelStatus.set('restarting');
    try { await restartKernel(); kernelStatus.set('ready'); }
    catch { kernelStatus.set('dead'); }
  }

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

<!-- App bar: owns the app name, the focused notebook's name, and the theme
     toggle. The toggle used to float over the canvas at top-right, which put
     it on top of a full-screen card's own toolbar. -->
<div class="app-bar">
  {#if $canvasState.focusedId}
    <!-- The full-screen card's controls. Same icons in the same order as a
         card's own toolbar on the canvas, so the row does not reshuffle when a
         notebook is zoomed in; only the full-screen icon flips to its inverse,
         because that is the one action whose meaning reverses. -->
    <button class="tb-btn tb-run-all" title="Run all cells"
            on:click={() => $focusedActions?.runAll()}>▶▶</button>
    <button class="tb-btn"
            title={$focusedActions?.horizontal ? 'Vertical layout' : 'Horizontal layout'}
            on:click={() => $focusedActions?.toggleLayout()}
    >{$focusedActions?.horizontal ? '↕' : '⇄'}</button>
    {#if $focusedActions?.hasSections}
      <button class="tb-btn"
              title={$focusedActions?.allSectionsCollapsed ? 'Expand all sections' : 'Collapse all sections'}
              on:click={() => $focusedActions?.toggleAllSections()}
      >{$focusedActions?.allSectionsCollapsed ? '⌄' : '⌃'}</button>
    {/if}
    <button class="tb-btn" title="Rename"
            on:click={() => $focusedActions?.rename()}>✎</button>
    <button class="tb-btn tb-focus" title="Back to canvas (pinch out)"
            on:click={() => setFocused(null)}>⤡</button>
    <button class="tb-btn" title="Collapse / expand"
            on:click={() => $focusedActions?.toggleCollapse()}
    >{$focusedActions?.collapsed ? '⊟' : '⊞'}</button>
    <button class="tb-btn tb-close" title="Close"
            on:click={() => $focusedActions?.close()}>✕</button>
  {:else}
    <span class="app-bar-name">Mathilda</span>
  {/if}
  {#if focusedTitle}<span class="app-bar-focus">{focusedTitle}</span>{/if}
  <span class="dark-toggle-spacer"></span>
  <button
    class="dark-toggle"
    title="Toggle dark mode"
    on:click={() => darkMode.update(v => !v)}
  >
    {$darkMode ? '◑' : '☀'}
  </button>
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
    /* Nothing here should swallow a drag meant for the window chrome. */
    user-select: none;
    -webkit-user-select: none;
  }

  .app-bar-name {
    font: 600 0.78rem/1 var(--sans);
    color: var(--text-h);
    letter-spacing: 0.01em;
  }

  /* Centred independently of the flanking items, so the notebook name sits in
     the middle of the WINDOW rather than the middle of the leftover space. */
  .app-bar-focus {
    position: absolute;
    left: 50%;
    transform: translateX(-50%);
    font: 500 0.78rem/1 var(--sans);
    color: var(--text);
    max-width: 50vw;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }

  /* Same visual language as the card toolbar these came from. */
  .app-bar .tb-btn {
    background: none;
    border: 1px solid transparent;
    border-radius: 5px;
    color: var(--text);
    font-size: 0.8rem;
    line-height: 1;
    padding: 0.18rem 0.36rem;
    cursor: pointer;
  }
  .app-bar .tb-btn:hover { background: var(--surface-2); color: var(--text-h); }
  .app-bar .tb-run-all { color: var(--ok, #4ade80); }
  /* Match the card toolbar's accents rather than position-based ones. */
  .app-bar .tb-focus { color: var(--accent, #89b4fa); }
  .app-bar .tb-close:hover { color: #f38ba8; }

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
