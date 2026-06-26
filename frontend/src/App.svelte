<script lang="ts">
  import { onMount, onDestroy } from 'svelte';
  import { open, save } from '@tauri-apps/plugin-dialog';
  import { listen } from '@tauri-apps/api/event';
  import CodeCell from './lib/CodeCell.svelte';
  import KernelStatus from './lib/KernelStatus.svelte';
  import {
    notebook, kernelStatus, selectedCells, clearSelection
  } from './lib/notebook';
  import type { OutputItem } from './lib/notebook';
  import {
    evaluateCell, restartKernel, interruptKernel,
    pingKernel, saveNotebook, loadNotebook
  } from './lib/ipc';
  import type { OutputMessage } from './lib/ipc';
  import { writable } from 'svelte/store';

  // ---------------------------------------------------------------------------
  // Dark / light mode

  const darkMode = writable(
    window.matchMedia?.('(prefers-color-scheme: dark)').matches ?? false
  );
  $: document.documentElement.classList.toggle('dark', $darkMode);

  // ---------------------------------------------------------------------------
  // Kernel init

  let unlisten: (() => void)[] = [];

  onMount(async () => {
    kernelStatus.set('starting');
    await new Promise(r => setTimeout(r, 1200));
    try {
      await pingKernel();
      kernelStatus.set('ready');
    } catch {
      kernelStatus.set('dead');
    }

    // Listen for native menu events emitted from Rust.
    unlisten.push(await listen('menu:open',       () => openFile()));
    unlisten.push(await listen('menu:save',       () => saveFile()));
    unlisten.push(await listen('menu:save-as',    () => saveFileAs()));
    unlisten.push(await listen('menu:add-cell',   () => addCellBelow()));
    unlisten.push(await listen('menu:run-all',    () => runAll()));
    unlisten.push(await listen('menu:restart',    () => restart()));
    unlisten.push(await listen('menu:interrupt',  () => interrupt()));
    unlisten.push(await listen('menu:toggle-dark',() => darkMode.update(v => !v)));
    unlisten.push(await listen('menu:delete-cells', () => deleteSelected()));
    unlisten.push(await listen('menu:copy-cells',   () => copySelected()));
    unlisten.push(await listen('menu:paste-cells',  () => pasteCells()));
  });

  onDestroy(() => unlisten.forEach(u => u()));

  // ---------------------------------------------------------------------------
  // Clipboard for cell copy/paste
  let copiedSources: string[] = [];

  function copySelected() {
    const sel = $selectedCells;
    if (sel.size === 0) return;
    copiedSources = $notebook
      .filter(c => sel.has(c.id))
      .map(c => c.source);
  }

  function pasteCells() {
    if (copiedSources.length === 0) return;
    const sel = $selectedCells;
    const lastId = sel.size > 0
      ? [...$notebook].reverse().find(c => sel.has(c.id))?.id
      : undefined;
    for (const src of copiedSources) {
      notebook.addCell(lastId, 'code');
      // set source on the just-added cell
      const cells = notebook.getCells();
      const newCell = cells[cells.length - 1];
      notebook.updateSource(newCell.id, src);
    }
  }

  function deleteSelected() {
    const sel = $selectedCells;
    if (sel.size === 0) return;
    notebook.removeCells(sel);
  }

  function addCellBelow() {
    const sel = $selectedCells;
    const lastId = sel.size > 0
      ? [...$notebook].reverse().find(c => sel.has(c.id))?.id
      : undefined;
    notebook.addCell(lastId, 'code');
  }

  // ---------------------------------------------------------------------------
  // Global keyboard handler for selected-cell operations.
  function onKeydown(e: KeyboardEvent) {
    const target = e.target as HTMLElement;
    const inEditor = target.closest('.cm-editor') != null;

    if (!inEditor) {
      if ((e.key === 'Delete' || e.key === 'Backspace') && $selectedCells.size > 0) {
        e.preventDefault();
        deleteSelected();
      }
      if (e.key === 'Escape') { e.preventDefault(); clearSelection(); }
    }

    if ((e.metaKey || e.ctrlKey) && e.key === 'c' && $selectedCells.size > 0 && !inEditor) {
      e.preventDefault(); copySelected();
    }
    if ((e.metaKey || e.ctrlKey) && e.key === 'v' && copiedSources.length > 0 && !inEditor) {
      e.preventDefault(); pasteCells();
    }
  }

  // ---------------------------------------------------------------------------
  // Cell execution

  async function runCell(cellId: string, source: string) {
    if (!source.trim()) return;
    const n = notebook.stampExec(cellId);
    notebook.clearOutput(cellId);
    notebook.setStatus(cellId, 'running');
    kernelStatus.set('busy');
    try {
      await evaluateCell(source, (msg: OutputMessage) => {
        const item = msgToOutputItem(msg);
        if (!item) return;
        if (msg.type === 'stream') notebook.appendStream(cellId, (msg as any).text ?? '');
        else notebook.appendOutput(cellId, item);
      });
      notebook.setStatus(cellId, 'done');
    } catch (e) {
      notebook.appendOutput(cellId, { kind: 'error', text: String(e) });
      notebook.setStatus(cellId, 'error');
      kernelStatus.set('dead');
    } finally {
      if ($kernelStatus !== 'dead') kernelStatus.set('ready');
    }
  }

  function msgToOutputItem(msg: OutputMessage): OutputItem | null {
    switch (msg.type) {
      case 'expr':   return { kind: 'expr',  text: msg.payload };
      case 'error':  return { kind: 'error', text: msg.message };
      case 'stream': return { kind: 'stream', text: (msg as any).text ?? '' };
      case 'plot':   return { kind: 'plot',  data: msg.payload };
      case 'html':   return { kind: 'html',  html: (msg as any).payload ?? '' };
      default:       return null;
    }
  }

  // ---------------------------------------------------------------------------
  // Run All

  async function runAll() {
    kernelStatus.set('restarting');
    notebook.resetExecCounter();
    try { await restartKernel(); kernelStatus.set('ready'); }
    catch { kernelStatus.set('dead'); return; }
    for (const cell of $notebook) {
      if (cell.type === 'code') {
        await runCell(cell.id, cell.source);
        if ($kernelStatus === 'dead') break;
      }
    }
  }

  async function interrupt() {
    try { await interruptKernel(); } catch {}
    kernelStatus.set('dead');
  }

  async function restart() {
    kernelStatus.set('restarting');
    notebook.resetExecCounter();
    for (const cell of $notebook) { notebook.clearOutput(cell.id); notebook.setStatus(cell.id, 'idle'); }
    try { await restartKernel(); kernelStatus.set('ready'); }
    catch { kernelStatus.set('dead'); }
  }

  // ---------------------------------------------------------------------------
  // File I/O

  let notebookPath: string | null = null;

  async function openFile() {
    const selected = await open({ filters: [{ name: 'Mathilda Notebook', extensions: ['mathilda', 'm'] }] });
    if (!selected) return;
    const path = typeof selected === 'string' ? selected : (selected as string[])[0];
    try { notebook.load(await loadNotebook(path)); notebookPath = path; }
    catch (e) { console.error('Open failed:', e); }
  }

  async function saveFile() {
    if (notebookPath) { await doSave(notebookPath); } else { await saveFileAs(); }
  }

  async function saveFileAs() {
    const path = await save({ defaultPath: 'notebook.mathilda', filters: [{ name: 'Mathilda Notebook', extensions: ['mathilda'] }] });
    if (!path) return;
    await doSave(path);
    notebookPath = path;
  }

  async function doSave(path: string) {
    try { await saveNotebook(path, notebook.serialize()); }
    catch (e) { console.error('Save failed:', e); }
  }

  // ---------------------------------------------------------------------------
  // Cell event handlers

  function handleRun(e: CustomEvent<{ id: string }>) {
    const cell = $notebook.find(c => c.id === e.detail.id);
    if (cell) runCell(cell.id, cell.source);
  }

  function handleChange(e: CustomEvent<{ id: string; source: string }>) {
    notebook.updateSource(e.detail.id, e.detail.source);
  }

  function handleAddBelow(e: CustomEvent<{ id: string }>) {
    notebook.addCell(e.detail.id, 'code');
  }
</script>

<svelte:window on:keydown={onKeydown} />

<div class="app">

  <!-- Minimal status bar — everything else is in the native menu -->
  <div class="statusbar">
    <span class="logo">Mathilda</span>
    <div class="statusbar-right">
      <button class="icon-btn" title="Toggle dark mode" on:click={() => darkMode.update(v => !v)}>
        {$darkMode ? '☀' : '◑'}
      </button>
      <KernelStatus />
    </div>
  </div>

  <!-- Notebook cells -->
  <main class="notebook" on:click={() => clearSelection()}>
    {#each $notebook as cell (cell.id)}
      {#if cell.type === 'code'}
        <CodeCell
          cellId={cell.id}
          source={cell.source}
          status={cell.status}
          output={cell.output}
          execIdx={cell.execIdx}
          on:run={handleRun}
          on:change={handleChange}
          on:addBelow={handleAddBelow}
        />
      {/if}
    {/each}

    <div class="add-cell-row">
      <button on:click|stopPropagation={addCellBelow}>＋ Add Cell</button>
    </div>
  </main>

  <!-- Dead kernel banner -->
  {#if $kernelStatus === 'dead'}
    <div class="kernel-banner">
      Kernel not running.
      <button on:click={restart}>Restart</button>
    </div>
  {/if}
</div>

<style>
  /* ---- CSS custom properties (light defaults) ---- */
  :global(:root) {
    --bg:          #f7f7f9;
    --surface:     #ffffff;
    --cell-bg:     #ffffff;
    --gutter-bg:   #f2f2f4;
    --gutter-hover:#e8e8ec;
    --border:      #e0e0e0;
    --text:        #1a1a1a;
    --text-muted:  #888;
    --accent:      #4a90e2;
    --accent-glow: rgba(74,144,226,0.15);
    --statusbar-bg:#ffffff;
    --out-text:    #222;
  }
  :global(.dark:root), :global(html.dark) {
    --bg:          #1e1e2e;
    --surface:     #2a2a3c;
    --cell-bg:     #252535;
    --gutter-bg:   #1e1e2e;
    --gutter-hover:#2e2e44;
    --border:      #3a3a52;
    --text:        #e0e0f0;
    --text-muted:  #888;
    --accent:      #7aa2f7;
    --accent-glow: rgba(122,162,247,0.15);
    --statusbar-bg:#16161e;
    --out-text:    #cdd6f4;
  }

  :global(*, *::before, *::after) { box-sizing: border-box; }
  :global(body) {
    margin: 0;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    background: var(--bg);
    color: var(--text);
  }

  .app {
    display: flex;
    flex-direction: column;
    height: 100vh;
    overflow: hidden;
  }

  /* Minimal status bar */
  .statusbar {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 0.25rem 0.75rem;
    background: var(--statusbar-bg);
    border-bottom: 1px solid var(--border);
    flex-shrink: 0;
    z-index: 10;
    height: 32px;
  }
  .logo {
    font-weight: 700;
    font-size: 0.9rem;
    letter-spacing: -0.3px;
    color: var(--text);
    opacity: 0.7;
  }
  .statusbar-right {
    display: flex;
    align-items: center;
    gap: 0.5rem;
  }
  .icon-btn {
    background: none;
    border: none;
    cursor: pointer;
    font-size: 1rem;
    color: var(--text-muted);
    padding: 2px 4px;
    border-radius: 3px;
    transition: color 0.1s;
  }
  .icon-btn:hover { color: var(--text); }

  /* Notebook scroll area — edge-to-edge, narrow side padding */
  .notebook {
    flex: 1;
    overflow-y: auto;
    padding: 0.75rem 1rem;
    max-width: 1100px;
    width: 100%;
    margin: 0 auto;
  }

  .add-cell-row {
    text-align: left;
    padding: 0.5rem 0 1.5rem;
  }
  .add-cell-row button {
    background: none;
    border: 1px dashed var(--border);
    color: var(--text-muted);
    border-radius: 5px;
    padding: 0.3rem 1rem;
    cursor: pointer;
    font-size: 0.82rem;
    transition: all 0.15s;
  }
  .add-cell-row button:hover { border-color: var(--accent); color: var(--accent); }

  /* Dead kernel banner */
  .kernel-banner {
    position: fixed;
    bottom: 1rem;
    left: 50%;
    transform: translateX(-50%);
    background: #c0392b;
    color: white;
    padding: 0.5rem 1rem;
    border-radius: 6px;
    font-size: 0.85rem;
    display: flex;
    align-items: center;
    gap: 0.7rem;
    box-shadow: 0 4px 12px rgba(0,0,0,0.25);
    z-index: 100;
  }
  .kernel-banner button {
    background: rgba(255,255,255,0.2);
    border: 1px solid rgba(255,255,255,0.4);
    color: white;
    border-radius: 4px;
    padding: 0.15rem 0.6rem;
    cursor: pointer;
    font-size: 0.8rem;
  }
</style>
