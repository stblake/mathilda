<script lang="ts">
  import { onMount, onDestroy, tick } from 'svelte';
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
  // Global keyboard handler for selected-cell operations and insertion point.
  function onKeydown(e: KeyboardEvent) {
    const target = e.target as HTMLElement;
    const inEditor = target.closest('.cm-editor') != null;

    // --- Insertion point mode ---
    if (insertionIdx !== null) {
      if (e.key === 'Escape') {
        e.preventDefault(); insertionIdx = null; return;
      }
      if (e.key === 'ArrowUp') {
        e.preventDefault();
        if (insertionIdx > 0) insertionIdx--;
        else { insertionIdx = null; cellFocusFns[$notebook[0]?.id]?.(); }
        return;
      }
      if (e.key === 'ArrowDown') {
        e.preventDefault();
        const len = $notebook.length;
        if (insertionIdx < len) insertionIdx++;
        else { insertionIdx = null; cellFocusFns[$notebook[len - 1]?.id]?.(); }
        return;
      }
      if (e.key === 'Enter') {
        e.preventDefault(); createCellAtInsertion(''); return;
      }
      // Any printable character → spawn a new cell and type into it.
      if (e.key.length === 1 && !e.ctrlKey && !e.metaKey) {
        e.preventDefault(); createCellAtInsertion(e.key); return;
      }
    }

    // --- Normal mode ---
    if (!inEditor) {
      if ((e.key === 'Delete' || e.key === 'Backspace') && $selectedCells.size > 0) {
        e.preventDefault(); deleteSelected();
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

  // Cell id → focus-editor function, populated via 'register' events.
  const cellFocusFns: Record<string, () => void> = {};

  function handleRegister(e: CustomEvent<{ id: string; fn: () => void }>) {
    cellFocusFns[e.detail.id] = e.detail.fn;
  }

  // ---------------------------------------------------------------------------
  // Insertion point — index of the gap where a new cell would be inserted.
  // null = inactive; 0 = before first cell; N = between cell[N-1] and cell[N].

  let insertionIdx: number | null = null;

  function handleFocusPrev(e: CustomEvent<{ id: string }>) {
    const cells = $notebook;
    const idx = cells.findIndex(c => c.id === e.detail.id);
    // Show insertion cursor in the gap above this cell.
    if (idx >= 0) {
      insertionIdx = idx;
      clearSelection();
    }
  }

  function handleFocusNext(e: CustomEvent<{ id: string }>) {
    const cells = $notebook;
    const idx = cells.findIndex(c => c.id === e.detail.id);
    // Show insertion cursor in the gap below this cell.
    if (idx >= 0 && idx + 1 <= cells.length) {
      insertionIdx = idx + 1;
      clearSelection();
    }
  }

  async function createCellAtInsertion(initialChar: string) {
    const idx = insertionIdx!;
    insertionIdx = null;
    const id = notebook.insertCellAt(idx, 'code', initialChar);
    await tick();
    // After tick the new CodeCell's onMount has run and registered its fn.
    const fn = cellFocusFns[id];
    if (fn) fn();
  }
</script>

<svelte:window on:keydown={onKeydown} />

<div class="app">

  <!-- Notebook — full bleed, no status bar -->
  <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
  <main class="notebook" on:click={() => { clearSelection(); insertionIdx = null; }}>

    <!-- Insertion point before first cell -->
    {#if insertionIdx === 0}
      <div class="insertion-cursor active"></div>
    {/if}

    {#each $notebook as cell, i (cell.id)}
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
          on:focusPrev={handleFocusPrev}
          on:focusNext={handleFocusNext}
          on:register={handleRegister}
        />
      {/if}

      <!-- Insertion point after this cell -->
      {#if insertionIdx === i + 1}
        <div class="insertion-cursor active"></div>
      {:else if i < $notebook.length - 1}
        <div class="insertion-cursor"></div>
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

  /* Full-bleed notebook — no status bar, no side margins */
  .notebook {
    flex: 1;
    overflow-y: auto;
    padding: 0.5rem 0.5rem;
    width: 100%;
  }

  /* Insertion cursor — thin animated line between cells */
  .insertion-cursor {
    height: 4px;
    margin: 1px 0;
    border-radius: 2px;
    transition: background 0.15s;
  }
  .insertion-cursor.active {
    background: var(--accent);
    animation: blink 1s ease-in-out infinite;
  }
  @keyframes blink {
    0%, 100% { opacity: 1; }
    50%       { opacity: 0.4; }
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
