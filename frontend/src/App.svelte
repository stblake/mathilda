<!--
  App.svelte — thin shell.
  Renders the full-viewport Canvas. Handles Cmd+S / Cmd+O and
  kernel status / dark mode toggles in a floating corner overlay.
-->
<script lang="ts">
  import { onMount, onDestroy } from 'svelte';
  import { writable, get } from 'svelte/store';
  import { open, save } from '@tauri-apps/plugin-dialog';
  import { listen } from '@tauri-apps/api/event';
  import { invoke } from '@tauri-apps/api/core';
  import Canvas from './lib/Canvas.svelte';
  import Toolbar from './lib/Toolbar.svelte';
  import { MENU_IDS, runMenuCommand } from './lib/menuCommands';
  import PropertiesPanel from './lib/PropertiesPanel.svelte';
  import SearchBar from './lib/SearchBar.svelte';
  import { searchOpen } from './lib/search';
  import { uiScale } from './lib/properties';
  import { kernelStatus, dirty, markClean } from './lib/notebook';
  import { darkMode } from './lib/theme';
  import { kernelMemory } from './lib/status';
  import { pingKernel, saveLibrary, loadLibrary, loadNotebook, saveNotebook,
           setWindowTitle as setTitleCmd,
           recentFiles, pushRecentFile, forgetRecentFile, clearRecentFiles } from './lib/ipc';
  import { restart, abortEvaluation } from './lib/kernelActions';
  import { serializeLibrary, loadLibraryData, canvasState, activeActions, activeFlags, setFocused,
           openNotebookCells, currentNotebook } from './lib/canvas';
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

      const hooks = { openFile, saveFile, saveFileAs, openRecent, clearRecent };
      /* One try PER ID, not one around the loop. `listen` throws on a name Tauri's event grammar
         rejects, and with a single try that first throw aborted the whole loop -- so one bad id
         (`file.new`, back when they were dotted) left the ENTIRE menu bar unsubscribed, including
         every legal id after it. Isolated, a bad id costs exactly its own menu item, and the
         warning names it instead of naming whichever happened to be first. */
      for (const id of MENU_IDS) {
        try {
          unlisten.push(await listen(`menu:${id}`, () => runMenuCommand(id, hooks)));
        } catch (e) { console.warn(`menu id "${id}" could not be subscribed:`, e); }
      }
    } catch (e) { console.warn('Menu listen error:', e); }
  });

  onDestroy(() => unlisten.forEach(u => u()));

  // --- Save-on-close prompt --------------------------------------------------
  /* A 3-way modal (Save / Don't Save / Cancel), shown when the app is closed
     with unsaved changes. The native dialog plugin is 2-button only, so this is
     an in-app modal; `savePrompt` holds the pending promise's resolver.

     EVERY close path is gated in Rust (src-tauri/src/lib.rs) and announced as the
     `quit-requested` event: the window's CloseRequested (red button, File > Close
     Window) AND the app's ExitRequested (Cmd+Q, dock > Quit). macOS routes Cmd+Q
     through terminate:, which never fires the window's CloseRequested — so
     listening to that event alone (as this once did) silently missed every Quit.
     We answer the gate by invoking `confirm_and_quit` / `cancel_quit`. */
  let savePrompt: ((v: 'save' | 'dont' | 'cancel') => void) | null = null;
  function promptSaveClose(): Promise<'save' | 'dont' | 'cancel'> {
    return new Promise((resolve) => { savePrompt = resolve; });
  }
  function answerSavePrompt(v: 'save' | 'dont' | 'cancel') {
    const r = savePrompt; savePrompt = null; r?.(v);
  }

  async function handleQuitRequest() {
    // Re-entrancy: a modal is already up (a second close gesture, or the Rust
    // ExitRequested handler re-emitting). The gate stays held; ignore duplicates
    // rather than stacking prompts.
    if (savePrompt) return;
    if (!get(dirty)) { await invoke('confirm_and_quit'); return; }      // nothing to lose
    const choice = await promptSaveClose();
    if (choice === 'cancel') { await invoke('cancel_quit'); return; }   // re-arm the gate, stay open
    if (choice === 'save') {
      const ok = await saveFile();
      if (!ok) { await invoke('cancel_quit'); return; }                 // save failed/cancelled -> stay open
    }
    await invoke('confirm_and_quit');                                   // Save ok, or Don't Save -> quit
  }

  onMount(async () => {
    try {
      unlisten.push(await listen('quit-requested', handleQuitRequest));
    } catch (e) { console.warn('quit handler:', e); }
  });

  // ---------------------------------------------------------------------------
  // File I/O — library-level (whole canvas)

  /* Two formats, told apart by extension:
       .lb   a LIBRARY -- the whole canvas as JSON; opening one replaces it.
       .mnb  ONE notebook ("Mathilda notebook") in the plain-text stanza format
             (src-tauri's notebook_format.rs); opening one adds it to the canvas.
     The dialog offered the notebook extension for a long time while every file
     went through the JSON library loader, so a notebook file could never open.

     `.mathilda` is what a notebook was written as before v0.235. It is still
     opened -- a file on disk outlives the decision to rename the format -- but
     never written, so Save As produces only `.mnb` from here on. */
  const isNotebookFile = (p: string) => /\.(mnb|mathilda)$/i.test(p);
  const baseName = (p: string) => p.split(/[\\/]/).pop() ?? p;

  /* Open one file by path, whichever of the two formats it is. Split out of openFile so the
     dialog and File > Open Recent share ONE implementation of "what opening means" -- a
     second copy of the .mnb/.lb branch is exactly how the two would drift.

     Returns whether it opened, which is what lets openRecent drop an entry that no longer
     works instead of leaving it in the menu to fail again. */
  async function openPath(path: string): Promise<boolean> {
    if (isNotebookFile(path)) {
      try {
        const cells = await loadNotebook(path);
        openNotebookCells(baseName(path).replace(/\.(mnb|mathilda)$/i, ''), cells);
      } catch (e) { console.error('Open failed:', e); return false; }
      pushRecentFile(path).catch(() => {});
      return true;
    }
    try {
      // Use loadLibrary (returns raw JSON string) not loadNotebook (parses as cells)
      const json = await loadLibrary(path);
      const title = loadLibraryData(json);
      libraryTitle = title;
      libraryPath  = path;
      const filename = path.split('/').pop()?.replace(/\.lb$/i, '') ?? title;
      setWindowTitle(filename);
      markClean(); // a freshly opened library matches disk
    } catch (e) { console.error('Open failed:', e); return false; }
    pushRecentFile(path).catch(() => {});
    return true;
  }

  async function openFile() {
    const sel = await open({
      filters: [
        { name: 'Mathilda Library or Notebook', extensions: ['lb', 'mnb', 'mathilda'] },
        { name: 'Mathilda Library', extensions: ['lb'] },
        { name: 'Mathilda Notebook', extensions: ['mnb', 'mathilda'] },
      ],
    });
    if (!sel) return;
    await openPath(typeof sel === 'string' ? sel : (sel as string[])[0]);
  }

  /* File > Open Recent > the i'th item. The native menu knows only the index it was built
     with, so the path comes back from the kernel-side list -- whose order IS the menu's, by
     construction in recent.rs.

     Re-read rather than cache: this window is not the only thing that can change the list
     (a second window, a later save), and a stale cache here would open the wrong file. */
  async function openRecent(i: number) {
    let list: string[];
    try { list = await recentFiles(); }
    catch (e) { console.error('Open Recent failed:', e); return; }
    const path = list[i];
    if (!path) return;
    if (!await openPath(path)) {
      // It was recorded once and will not open now -- moved, deleted, or unreadable. Drop it
      // rather than leave a menu entry that only ever fails.
      forgetRecentFile(path).catch(() => {});
    }
  }

  async function clearRecent() {
    try { await clearRecentFiles(); }
    catch (e) { console.error('Clear Menu failed:', e); }
  }

  async function saveFile(): Promise<boolean> {
    return libraryPath ? await doSave(libraryPath) : await saveFileAs();
  }

  /* Save As offers both formats. Choosing .mnb EXPORTS the current notebook
     (the active pane, else the top card) as sources only; it does not become the
     library's path, so a later Cmd+S still saves the whole canvas as .lb rather
     than overwriting the notebook file with something else. */
  async function saveFileAs(): Promise<boolean> {
    const path = await save({
      defaultPath: (libraryTitle || 'library') + '.lb',
      filters: [
        { name: 'Mathilda Library', extensions: ['lb'] },
        { name: 'Mathilda Notebook (current notebook, no outputs)', extensions: ['mnb'] },
      ],
    });
    if (!path) return false;
    if (isNotebookFile(path)) {
      const nb = currentNotebook();
      if (!nb) { console.error('Save failed: no notebook to save'); return false; }
      try {
        await saveNotebook(path, nb.store.serializeLegacy());
        pushRecentFile(path).catch(() => {});
        markClean();
        return true;
      }
      catch (e) { console.error('Save failed:', e); return false; }
    }
    libraryPath = path;
    return await doSave(path);
  }

  async function doSave(path: string): Promise<boolean> {
    try {
      const json = serializeLibrary(libraryTitle);
      await saveLibrary(path, json);
      const filename = path.split('/').pop()?.replace(/\.lb$/i, '') ?? 'Library';
      libraryTitle = filename;
      setWindowTitle(filename);
      pushRecentFile(path).catch(() => {});
      markClean();
      return true;
    } catch (e) { console.error('Save failed:', e); return false; }
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
  /* Global UI scale (Cmd+= zoom in, Cmd+- zoom out, Cmd+0 reset).
     A store rather than a local, so the properties panel can offer the SAME value
     the keyboard drives. Two independent scales would drift apart. */
  $: document.documentElement.style.fontSize = `${$uiScale * 16}px`;

  function onKeydown(e: KeyboardEvent) {
    /* The save-on-close modal owns the keyboard while it is up: Enter = Save,
       Escape = Cancel. */
    if (savePrompt) {
      if (e.key === 'Escape') { e.preventDefault(); answerSavePrompt('cancel'); }
      else if (e.key === 'Enter') { e.preventDefault(); answerSavePrompt('save'); }
      return;
    }
    const mod = e.metaKey || e.ctrlKey;
    if (!mod) return;
    if (e.key === 's' || e.key === 'S') { e.preventDefault(); saveFile(); return; }
    if (e.key === 'o' || e.key === 'O') { e.preventDefault(); openFile(); return; }
    /* Cmd+F is free: @codemirror/search is not installed, so no editor claims it.
       Focused mode only -- on the canvas there is no notebook to search. */
    if (e.key === 'f' || e.key === 'F') {
      if ($canvasState.focusedIds.length) { e.preventDefault(); searchOpen.set(true); }
      return;
    }
    // Cmd+= / Cmd++ → scale up; Cmd+- → scale down; Cmd+0 → reset
    if (e.key === '=' || e.key === '+') {
      e.preventDefault(); uiScale.update(v => Math.min(2.0, +(v + 0.1).toFixed(1)));
    } else if (e.key === '-' || e.key === '_') {
      e.preventDefault(); uiScale.update(v => Math.max(0.5, +(v - 0.1).toFixed(1)));
    } else if (e.key === '0') {
      /* The comment above has promised this since it was written. */
      e.preventDefault(); uiScale.set(1.0);
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

<!-- Properties sidebar. Focused mode only: it reports on the notebook in the
     active pane, and on the canvas there is no active pane to report on. -->
{#if $canvasState.focusedIds.length}
  <PropertiesPanel />
  <SearchBar />
{/if}

<!-- Kernel dead banner -->
{#if $kernelStatus === 'dead'}
  <div class="kernel-banner">
    Kernel not running.
    <button on:click={restart}>Restart</button>
  </div>
{/if}

<!-- Save-before-closing prompt (unsaved changes on window close) -->
{#if savePrompt}
  <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
  <div class="save-modal-backdrop" on:click={() => answerSavePrompt('cancel')}>
    <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
    <div class="save-modal" role="dialog" aria-modal="true" tabindex="-1" on:click|stopPropagation>
      <div class="save-modal-title">Save changes before closing?</div>
      <div class="save-modal-msg">Your notebook has unsaved changes. If you don’t save, they will be lost.</div>
      <div class="save-modal-buttons">
        <button class="btn-dont" on:click={() => answerSavePrompt('dont')}>Don’t Save</button>
        <span class="save-modal-spacer"></span>
        <button class="btn-cancel" on:click={() => answerSavePrompt('cancel')}>Cancel</button>
        <button class="btn-save" on:click={() => answerSavePrompt('save')}>Save</button>
      </div>
    </div>
  </div>
{/if}

<style>
  /* ---- Save-before-closing modal ---- */
  .save-modal-backdrop {
    position: fixed; inset: 0;
    background: rgba(0, 0, 0, 0.45);
    display: flex; align-items: center; justify-content: center;
    z-index: var(--z-menu, 400);
  }
  .save-modal {
    background: var(--menu-bg, #1a1b2e);
    color: var(--text, #cdd6f4);
    border: 1px solid var(--menu-border, rgba(255, 255, 255, 0.12));
    border-radius: 10px;
    box-shadow: var(--menu-shadow, 0 20px 50px rgba(0, 0, 0, 0.5));
    padding: 20px 22px 16px;
    width: min(420px, 90vw);
  }
  .save-modal-title { font-size: 1.02rem; font-weight: 650; margin-bottom: 8px; }
  .save-modal-msg {
    font-size: 0.86rem; color: var(--text-muted, #9399b2);
    line-height: 1.5; margin-bottom: 18px;
  }
  .save-modal-buttons { display: flex; align-items: center; gap: 8px; }
  .save-modal-spacer { flex: 1; }
  .save-modal-buttons button {
    font-family: inherit; font-size: 0.85rem;
    padding: 6px 14px; border-radius: 6px; cursor: pointer;
    border: 1px solid var(--menu-border, rgba(255, 255, 255, 0.14));
    background: var(--surface-2, rgba(255, 255, 255, 0.06));
    color: var(--text, #cdd6f4);
  }
  .save-modal-buttons button:hover { background: var(--surface-3, rgba(255, 255, 255, 0.12)); }
  .save-modal-buttons .btn-save {
    background: var(--accent, #89b4fa); border-color: var(--accent, #89b4fa);
    color: #0b1020; font-weight: 600;
  }
  .save-modal-buttons .btn-save:hover { filter: brightness(1.07); }
  .save-modal-buttons .btn-dont { color: var(--err, #f38ba8); }

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

    /* Code-cell syntax highlighting (CodeMirror). Dark palette (Catppuccin
       Mocha family, matching --accent). mathildaLang.ts references these. */
    --cm-comment:  #6c7086;
    --cm-string:   #a6e3a1;
    --cm-number:   #fab387;
    --cm-builtin:  #89b4fa;
    --cm-symbol:   #cdd6f4;
    --cm-pattern:  #f38ba8;
    --cm-slot:     #f9e2af;
    --cm-out:      #f9e2af;
    --cm-operator: #89dceb;
    --cm-bracket:  #9399b2;
    --cm-error:    #f38ba8;
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

    /* Code-cell syntax highlighting (CodeMirror). Light palette (Catppuccin
       Latte family), overriding the dark defaults above when html.light is set. */
    --cm-comment:  #8c8fa1;
    --cm-string:   #40a02b;
    --cm-number:   #fe640b;
    --cm-builtin:  #1e66f5;
    --cm-symbol:   #4c4f69;
    --cm-pattern:  #d20f39;
    --cm-slot:     #df8e1d;
    --cm-out:      #df8e1d;
    --cm-operator: #209fb5;
    --cm-bracket:  #7c7f93;
    --cm-error:    #d20f39;
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
    z-index: var(--z-appbar);
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
    z-index: var(--z-banner);
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
