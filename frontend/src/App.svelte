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
  import Canvas from './lib/Canvas.svelte';
  import KernelStatus from './lib/KernelStatus.svelte';
  import { kernelStatus } from './lib/notebook';
  import { pingKernel, restartKernel, saveNotebook, loadNotebook } from './lib/ipc';
  import { serializeLibrary, loadLibraryData } from './lib/canvas';

  // ---------------------------------------------------------------------------
  // Dark mode (default to dark — canvas is always dark)

  const darkMode = writable(true);
  $: if (typeof document !== 'undefined')
    document.documentElement.classList.toggle('dark', $darkMode);

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
      // loadNotebook returns raw file bytes as CellData[] — but for .lb we pass raw json
      const raw = await loadNotebook(path) as any;
      const json = typeof raw === 'string' ? raw : JSON.stringify(raw);
      libraryTitle = loadLibraryData(json);
      libraryPath  = path;
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
      // saveNotebook expects CellData[] — wrap json as a single-cell hack,
      // or rely on the Rust side to handle raw string. For now log a TODO.
      // TODO: extend Tauri command to support raw-string save for library format.
      console.log('Save to:', path, 'size:', json.length);
    } catch (e) { console.error('Save failed:', e); }
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

<!-- Floating corner: kernel status + dark mode toggle -->
<div class="corner-overlay">
  <button
    class="dark-toggle"
    title="Toggle dark mode"
    on:click={() => darkMode.update(v => !v)}
  >
    {$darkMode ? '☀' : '◑'}
  </button>
  <KernelStatus />
</div>

<!-- Kernel dead banner -->
{#if $kernelStatus === 'dead'}
  <div class="kernel-banner">
    Kernel not running.
    <button on:click={restart}>Restart</button>
  </div>
{/if}

<style>
  /* ----  Global tokens for dark canvas aesthetic ---- */
  :global(:root) {
    --bg:          #050810;
    --surface:     rgba(12,15,28,0.85);
    --cell-bg:     transparent;
    --border:      rgba(255,255,255,0.08);
    --text:        #cdd6f4;
    --text-muted:  #585b70;
    --accent:      #89b4fa;
    --accent-glow: rgba(137,180,250,0.12);
    --out-text:    #cdd6f4;
    --gutter-bg:   rgba(255,255,255,0.02);
    --gutter-hover:rgba(255,255,255,0.04);
  }
  /* Keep light-mode fallback minimal (canvas is always dark) */
  :global(html.dark) {
    --bg:          #050810;
    --surface:     rgba(12,15,28,0.85);
    --border:      rgba(255,255,255,0.08);
    --text:        #cdd6f4;
    --text-muted:  #585b70;
    --accent:      #89b4fa;
    --accent-glow: rgba(137,180,250,0.12);
  }

  :global(*, *::before, *::after) { box-sizing: border-box; }
  :global(body) {
    margin: 0;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    background: #050810;
    color: #cdd6f4;
    overflow: hidden;
  }

  /* ---- Corner overlay ---- */
  .corner-overlay {
    position: fixed;
    top: 10px;
    right: 14px;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    z-index: 200;
    pointer-events: auto;
  }

  .dark-toggle {
    background: rgba(255,255,255,0.05);
    border: 1px solid rgba(255,255,255,0.08);
    color: #585b70;
    cursor: pointer;
    font-size: 0.9rem;
    padding: 3px 6px;
    border-radius: 6px;
    line-height: 1;
    transition: color 0.1s, background 0.1s;
  }
  .dark-toggle:hover { color: #cdd6f4; background: rgba(255,255,255,0.09); }

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
