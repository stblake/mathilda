<script lang="ts">
  import { onMount, onDestroy, tick } from 'svelte';
  import { open, save } from '@tauri-apps/plugin-dialog';
  import { listen } from '@tauri-apps/api/event';
  import CellShell from './lib/CellShell.svelte';
  import KernelStatus from './lib/KernelStatus.svelte';
  import { notebook, kernelStatus, selectedCells, clearSelection } from './lib/notebook';
  import type { OutputItem, CellType } from './lib/notebook';
  import {
    evaluateCell, restartKernel, interruptKernel,
    pingKernel, saveNotebook, loadNotebook,
  } from './lib/ipc';
  import type { OutputMessage } from './lib/ipc';
  import { writable } from 'svelte/store';

  // ---------------------------------------------------------------------------
  // Dark mode

  const darkMode = writable(
    typeof window !== 'undefined' && window.matchMedia?.('(prefers-color-scheme: dark)').matches
  );
  $: if (typeof document !== 'undefined')
      document.documentElement.classList.toggle('dark', $darkMode);

  // ---------------------------------------------------------------------------
  // Kernel init + menu listeners

  let unlisten: (() => void)[] = [];

  onMount(async () => {
    kernelStatus.set('starting');
    await new Promise(r => setTimeout(r, 1200));
    try { await pingKernel(); kernelStatus.set('ready'); }
    catch  { kernelStatus.set('dead'); }

    try {
      unlisten.push(await listen('menu:open',        () => openFile()));
      unlisten.push(await listen('menu:save',        () => saveFile()));
      unlisten.push(await listen('menu:save-as',     () => saveFileAs()));
      unlisten.push(await listen('menu:add-cell',    () => addRow()));
      unlisten.push(await listen('menu:run-all',     () => runAll()));
      unlisten.push(await listen('menu:restart',     () => restart()));
      unlisten.push(await listen('menu:interrupt',   () => interrupt()));
      unlisten.push(await listen('menu:toggle-dark', () => darkMode.update(v => !v)));
      unlisten.push(await listen('menu:delete-cells',() => deleteSelected()));
    } catch (e) { console.warn('Menu listen error:', e); }
  });

  onDestroy(() => unlisten.forEach(u => u()));

  // ---------------------------------------------------------------------------
  // Cell focus registry

  const cellFocusFns: Record<string, () => void> = {};

  function handleRegister(e: CustomEvent<{ id: string; fn: () => void }>) {
    cellFocusFns[e.detail.id] = e.detail.fn;
  }

  // ---------------------------------------------------------------------------
  // Insertion point (row-level gaps)

  let insertionRowIdx: number | null = null;

  function handleFocusPrev(e: CustomEvent<{ id: string }>) {
    const rows = notebook.getRows();
    for (let ri = 0; ri < rows.length; ri++) {
      if (rows[ri].cells.some(c => c.id === e.detail.id)) {
        insertionRowIdx = ri;
        clearSelection();
        return;
      }
    }
  }

  function handleFocusNext(e: CustomEvent<{ id: string }>) {
    const rows = notebook.getRows();
    for (let ri = 0; ri < rows.length; ri++) {
      if (rows[ri].cells.some(c => c.id === e.detail.id)) {
        if (ri + 1 <= rows.length) { insertionRowIdx = ri + 1; clearSelection(); }
        return;
      }
    }
  }

  async function createRowAtInsertion(initialChar: string, type: CellType = 'code') {
    const idx = insertionRowIdx!;
    insertionRowIdx = null;
    const id = notebook.insertRowAt(idx, type, initialChar);
    await tick();
    cellFocusFns[id]?.();
  }

  // ---------------------------------------------------------------------------
  // 4-directional add

  async function addRowAbove(e: CustomEvent<{ rowId: string }>) {
    const rows = notebook.getRows();
    const ri = rows.findIndex(r => r.id === e.detail.rowId);
    if (ri < 0) return;
    const id = notebook.insertRowAt(ri);
    await tick(); cellFocusFns[id]?.();
  }

  async function addRowBelow(e: CustomEvent<{ rowId: string }>) {
    const rows = notebook.getRows();
    const ri = rows.findIndex(r => r.id === e.detail.rowId);
    const id = notebook.insertRowAt(ri + 1);
    await tick(); cellFocusFns[id]?.();
  }

  async function addCellLeft(e: CustomEvent<{ rowId: string; cellIdx: number }>) {
    const id = notebook.insertCellInRow(e.detail.rowId, e.detail.cellIdx);
    await tick(); cellFocusFns[id]?.();
  }

  async function addCellRight(e: CustomEvent<{ rowId: string; cellIdx: number }>) {
    const id = notebook.insertCellInRow(e.detail.rowId, e.detail.cellIdx);
    await tick(); cellFocusFns[id]?.();
  }

  function addRow(type: CellType = 'code') {
    const id = notebook.addRow(type);
    tick().then(() => cellFocusFns[id]?.());
  }

  // ---------------------------------------------------------------------------
  // Clipboard

  let copiedSources: Array<{ type: CellType; source: string }> = [];

  function deleteSelected() {
    const sel = $selectedCells;
    if (sel.size === 0) return;
    notebook.removeCells(sel);
  }

  // ---------------------------------------------------------------------------
  // Keyboard

  function onKeydown(e: KeyboardEvent) {
    const target = e.target as HTMLElement;
    const inEditor = target.closest('.cm-editor') != null || target.isContentEditable;

    if (insertionRowIdx !== null) {
      if (e.key === 'Escape') { e.preventDefault(); insertionRowIdx = null; return; }
      if (e.key === 'ArrowUp') {
        e.preventDefault();
        if (insertionRowIdx > 0) insertionRowIdx--;
        else { insertionRowIdx = null; const rows = notebook.getRows(); cellFocusFns[rows[0]?.cells[0]?.id]?.(); }
        return;
      }
      if (e.key === 'ArrowDown') {
        e.preventDefault();
        const len = notebook.getRows().length;
        if (insertionRowIdx < len) insertionRowIdx++;
        else { insertionRowIdx = null; const rows = notebook.getRows(); const last = rows[rows.length-1]; cellFocusFns[last?.cells[last?.cells.length-1]?.id]?.(); }
        return;
      }
      if (e.key === 'Enter') { e.preventDefault(); createRowAtInsertion(''); return; }
      if (e.key.length === 1 && !e.ctrlKey && !e.metaKey) { e.preventDefault(); createRowAtInsertion(e.key); return; }
    }

    if (!inEditor) {
      if ((e.key === 'Delete' || e.key === 'Backspace') && $selectedCells.size > 0) {
        e.preventDefault(); deleteSelected();
      }
      if (e.key === 'Escape') { e.preventDefault(); clearSelection(); }
    }
  }

  // ---------------------------------------------------------------------------
  // Cell execution

  async function runCell(cellId: string, source: string) {
    if (!source.trim()) return;
    notebook.stampExec(cellId);
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

  async function runAll() {
    kernelStatus.set('restarting');
    notebook.resetExecCounter();
    try { await restartKernel(); kernelStatus.set('ready'); }
    catch { kernelStatus.set('dead'); return; }
    for (const cell of notebook.allCells()) {
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
    for (const cell of notebook.allCells()) {
      notebook.clearOutput(cell.id);
      notebook.setStatus(cell.id, 'idle');
    }
    try { await restartKernel(); kernelStatus.set('ready'); }
    catch { kernelStatus.set('dead'); }
  }

  // ---------------------------------------------------------------------------
  // File I/O

  let notebookPath: string | null = null;

  async function openFile() {
    const sel = await open({ filters: [{ name: 'Mathilda Notebook', extensions: ['mathilda', 'm'] }] });
    if (!sel) return;
    const path = typeof sel === 'string' ? sel : (sel as string[])[0];
    try {
      const raw = await loadNotebook(path);
      // raw is Array<{type, source}> (legacy) or could be extended
      notebook.loadLegacy(raw as any);
      notebookPath = path;
    } catch (e) { console.error('Open failed:', e); }
  }

  async function saveFile() {
    if (notebookPath) doSave(notebookPath); else saveFileAs();
  }
  async function saveFileAs() {
    const path = await save({ defaultPath: 'notebook.mathilda', filters: [{ name: 'Mathilda Notebook', extensions: ['mathilda'] }] });
    if (!path) return;
    await doSave(path); notebookPath = path;
  }
  async function doSave(path: string) {
    try { await saveNotebook(path, notebook.serializeLegacy() as any); }
    catch (e) { console.error('Save failed:', e); }
  }

  // ---------------------------------------------------------------------------
  // Event handlers from CellShell

  function handleRun(e: CustomEvent<{ id: string }>) {
    const cell = notebook.allCells().find(c => c.id === e.detail.id);
    if (cell) runCell(cell.id, cell.source);
  }

  function handleChange(e: CustomEvent<{ id: string; source: string }>) {
    notebook.updateSource(e.detail.id, e.detail.source);
  }

  function handleAddBelow(e: CustomEvent<{ rowId: string }>) {
    addRowBelow(e);
  }
</script>

<svelte:window on:keydown={onKeydown} />

<!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
<div class="app" on:click={() => { clearSelection(); insertionRowIdx = null; }}>

  <!-- Kernel status corner -->
  <div class="corner-status">
    <button class="dark-toggle" on:click|stopPropagation={() => darkMode.update(v => !v)} title="Toggle dark mode">
      {$darkMode ? '☀' : '◑'}
    </button>
    <KernelStatus />
  </div>

  <!-- Notebook: list of rows -->
  <main class="notebook">

    {#if insertionRowIdx === 0}
      <div class="row-cursor active"></div>
    {/if}

    {#each $notebook as row, ri (row.id)}
      <div class="nb-row">
        {#each row.cells as cell, ci (cell.id)}
          <div class="cell-col" style="flex: 1 1 0; min-width: 0">
            <CellShell
              {cell}
              rowId={row.id}
              cellIdx={ci}
              on:run={handleRun}
              on:change={handleChange}
              on:addAbove={addRowAbove}
              on:addBelow={addRowBelow}
              on:addLeft={addCellLeft}
              on:addRight={addCellRight}
              on:focusPrev={handleFocusPrev}
              on:focusNext={handleFocusNext}
              on:register={handleRegister}
            />
          </div>
          {#if ci < row.cells.length - 1}
            <div class="col-divider"></div>
          {/if}
        {/each}
      </div>

      {#if insertionRowIdx === ri + 1}
        <div class="row-cursor active"></div>
      {:else if ri < $notebook.length - 1}
        <div class="row-cursor"></div>
      {/if}
    {/each}

    <div class="add-row-bar">
      <button on:click|stopPropagation={() => addRow('code')}>＋ Code</button>
      <button on:click|stopPropagation={() => addRow('text')}>＋ Text</button>
      <button on:click|stopPropagation={() => addRow('section')}>＋ Section</button>
    </div>
  </main>

  {#if $kernelStatus === 'dead'}
    <div class="kernel-banner">
      Kernel not running.
      <button on:click={restart}>Restart</button>
    </div>
  {/if}
</div>

<style>
  :global(:root) {
    --bg:          #f7f7f9;
    --surface:     #ffffff;
    --cell-bg:     #ffffff;
    --gutter-bg:   #f2f2f4;
    --gutter-hover:#e8e8ec;
    --border:      #e4e4e4;
    --text:        #1a1a1a;
    --text-muted:  #999;
    --accent:      #4a90e2;
    --accent-glow: rgba(74,144,226,0.15);
    --out-text:    #222;
  }
  :global(html.dark) {
    --bg:          #1e1e2e;
    --surface:     #2a2a3c;
    --cell-bg:     #1e1e2e;
    --gutter-bg:   #16161e;
    --gutter-hover:#2e2e44;
    --border:      #2e2e44;
    --text:        #cdd6f4;
    --text-muted:  #6c7086;
    --accent:      #89b4fa;
    --accent-glow: rgba(137,180,250,0.12);
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
    position: relative;
  }

  /* Corner status — floats top-right, doesn't take layout space */
  .corner-status {
    position: fixed;
    top: 6px;
    right: 12px;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    z-index: 50;
  }
  .dark-toggle {
    background: none;
    border: none;
    cursor: pointer;
    font-size: 0.95rem;
    color: var(--text-muted);
    padding: 2px 4px;
    border-radius: 4px;
    transition: color 0.1s;
  }
  .dark-toggle:hover { color: var(--text); }

  /* Full-bleed notebook */
  .notebook {
    flex: 1;
    overflow-y: auto;
    padding: 0;
    width: 100%;
  }

  /* Rows */
  .nb-row {
    display: flex;
    flex-direction: row;
    align-items: stretch;
    border-bottom: 1px solid var(--border);
  }

  .col-divider {
    width: 1px;
    background: var(--border);
    flex-shrink: 0;
  }

  /* Row insertion cursor */
  .row-cursor {
    height: 3px;
    margin: 0;
    transition: background 0.12s;
  }
  .row-cursor.active {
    background: var(--accent);
    animation: blink 1s ease-in-out infinite;
  }
  @keyframes blink { 0%,100%{opacity:1} 50%{opacity:0.35} }

  /* Bottom add bar */
  .add-row-bar {
    padding: 0.6rem 0.5rem;
    display: flex;
    gap: 0.4rem;
  }
  .add-row-bar button {
    background: none;
    border: 1px dashed var(--border);
    color: var(--text-muted);
    border-radius: 4px;
    padding: 0.25rem 0.75rem;
    font-size: 0.8rem;
    cursor: pointer;
    transition: all 0.12s;
  }
  .add-row-bar button:hover { border-color: var(--accent); color: var(--accent); }

  /* Kernel banner */
  .kernel-banner {
    position: fixed;
    bottom: 1rem;
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
    box-shadow: 0 4px 12px rgba(0,0,0,0.25);
    z-index: 100;
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
