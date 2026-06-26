<script lang="ts">
  import { onMount } from 'svelte';
  import { open, save } from '@tauri-apps/plugin-dialog';
  import CodeCell from './lib/CodeCell.svelte';
  import KernelStatus from './lib/KernelStatus.svelte';
  import { notebook, kernelStatus } from './lib/notebook';
  import type { OutputItem } from './lib/notebook';
  import {
    evaluateCell,
    restartKernel,
    interruptKernel,
    pingKernel,
    saveNotebook,
    loadNotebook,
  } from './lib/ipc';
  import type { OutputMessage } from './lib/ipc';

  // ---------------------------------------------------------------------------
  // Kernel init

  onMount(async () => {
    kernelStatus.set('starting');
    // Give the async spawn task in Rust a moment to complete.
    await new Promise(r => setTimeout(r, 1500));
    try {
      await pingKernel();
      kernelStatus.set('ready');
    } catch (e) {
      console.error('Kernel not ready:', e);
      kernelStatus.set('dead');
    }
  });

  // ---------------------------------------------------------------------------
  // Cell execution

  async function runCell(cellId: string, source: string) {
    if (!source.trim()) return;
    notebook.clearOutput(cellId);
    notebook.setStatus(cellId, 'running');
    kernelStatus.set('busy');
    try {
      await evaluateCell(source, (msg: OutputMessage) => {
        const item = msgToOutputItem(msg);
        if (!item) return;
        if (msg.type === 'stream') {
          notebook.appendStream(cellId, (msg as any).text ?? '');
        } else {
          notebook.appendOutput(cellId, item);
        }
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
    try {
      await restartKernel();
      kernelStatus.set('ready');
    } catch (e) {
      kernelStatus.set('dead');
      return;
    }
    for (const cell of $notebook) {
      if (cell.type === 'code') {
        await runCell(cell.id, cell.source);
        if ($kernelStatus === 'dead') break;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Interrupt / Restart

  async function interrupt() {
    try {
      await interruptKernel();
    } catch (e) {
      console.error('Interrupt failed:', e);
    }
    kernelStatus.set('dead');
  }

  async function restart() {
    kernelStatus.set('restarting');
    for (const cell of $notebook) {
      notebook.clearOutput(cell.id);
      notebook.setStatus(cell.id, 'idle');
    }
    try {
      await restartKernel();
      kernelStatus.set('ready');
    } catch (e) {
      kernelStatus.set('dead');
    }
  }

  // ---------------------------------------------------------------------------
  // File I/O

  let notebookPath: string | null = null;

  async function openFile() {
    const selected = await open({
      filters: [{ name: 'Mathilda Notebook', extensions: ['mathilda', 'm'] }],
    });
    if (!selected) return;
    const path = typeof selected === 'string' ? selected : (selected as string[])[0];
    try {
      const cells = await loadNotebook(path);
      notebook.load(cells);
      notebookPath = path;
    } catch (e) {
      console.error('Open failed:', e);
    }
  }

  async function saveFile() {
    const path =
      notebookPath ??
      (await save({
        defaultPath: 'notebook.mathilda',
        filters: [{ name: 'Mathilda Notebook', extensions: ['mathilda'] }],
      }));
    if (!path) return;
    try {
      await saveNotebook(path, notebook.serialize());
      notebookPath = path;
    } catch (e) {
      console.error('Save failed:', e);
    }
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

  function handleRemove(e: CustomEvent<{ id: string }>) {
    notebook.removeCell(e.detail.id);
  }
</script>

<div class="app">

  <!-- Toolbar -->
  <header class="toolbar">
    <div class="toolbar-left">
      <span class="logo">Mathilda</span>
    </div>
    <div class="toolbar-center">
      <button on:click={() => notebook.addCell(undefined, 'code')} title="New code cell">＋ Cell</button>
      <button on:click={runAll} title="Restart kernel and run all cells top-to-bottom">▶▶ Run All</button>
      <button on:click={interrupt} title="Kill running computation" class="btn-danger">■ Interrupt</button>
      <button on:click={restart} title="Restart kernel">↺ Restart</button>
      <span class="sep">|</span>
      <button on:click={openFile}>📂 Open</button>
      <button on:click={saveFile}>💾 Save</button>
    </div>
    <div class="toolbar-right">
      <KernelStatus />
    </div>
  </header>

  <!-- Notebook cells -->
  <main class="notebook">
    {#each $notebook as cell (cell.id)}
      {#if cell.type === 'code'}
        <CodeCell
          cellId={cell.id}
          source={cell.source}
          status={cell.status}
          output={cell.output}
          on:run={handleRun}
          on:change={handleChange}
          on:addBelow={handleAddBelow}
          on:remove={handleRemove}
        />
      {/if}
    {/each}

    <div class="add-cell-row">
      <button on:click={() => notebook.addCell(undefined, 'code')}>＋ Add Cell</button>
    </div>
  </main>

  <!-- Dead kernel banner -->
  {#if $kernelStatus === 'dead'}
    <div class="kernel-banner">
      Kernel not running.
      <button on:click={restart}>Restart Kernel</button>
    </div>
  {/if}
</div>

<style>
  :global(*, *::before, *::after) { box-sizing: border-box; }
  :global(body) {
    margin: 0;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    background: #f7f7f9;
    color: #222;
  }

  .app {
    display: flex;
    flex-direction: column;
    height: 100vh;
    overflow: hidden;
  }

  .toolbar {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 0.4rem 1rem;
    background: #fff;
    border-bottom: 1px solid #e0e0e0;
    box-shadow: 0 1px 3px rgba(0,0,0,0.06);
    flex-shrink: 0;
    z-index: 10;
    gap: 0.5rem;
  }

  .logo {
    font-weight: 700;
    font-size: 1.1rem;
    letter-spacing: -0.5px;
    color: #1a1a2e;
  }

  .toolbar-center {
    display: flex;
    align-items: center;
    gap: 0.3rem;
    flex-wrap: wrap;
  }

  .toolbar button {
    background: #f0f0f0;
    border: 1px solid #ddd;
    border-radius: 5px;
    padding: 0.28rem 0.65rem;
    font-size: 0.8rem;
    cursor: pointer;
    transition: background 0.1s;
    white-space: nowrap;
    color: #333;
  }
  .toolbar button:hover { background: #e4e4e4; }
  .btn-danger { color: #c0392b !important; }
  .btn-danger:hover { background: #fdf0f0 !important; }
  .sep { color: #ddd; user-select: none; }

  .notebook {
    flex: 1;
    overflow-y: auto;
    padding: 1.5rem 2rem;
    max-width: 900px;
    width: 100%;
    margin: 0 auto;
  }

  .add-cell-row {
    text-align: center;
    padding: 1rem 0 2rem;
  }
  .add-cell-row button {
    background: none;
    border: 1px dashed #bbb;
    color: #888;
    border-radius: 5px;
    padding: 0.4rem 1.4rem;
    cursor: pointer;
    font-size: 0.85rem;
    transition: all 0.15s;
  }
  .add-cell-row button:hover { border-color: #4a90e2; color: #4a90e2; }

  .kernel-banner {
    position: fixed;
    bottom: 1.2rem;
    left: 50%;
    transform: translateX(-50%);
    background: #c0392b;
    color: white;
    padding: 0.6rem 1.2rem;
    border-radius: 8px;
    font-size: 0.88rem;
    display: flex;
    align-items: center;
    gap: 0.8rem;
    box-shadow: 0 4px 12px rgba(0,0,0,0.2);
    z-index: 100;
  }
  .kernel-banner button {
    background: rgba(255,255,255,0.2);
    border: 1px solid rgba(255,255,255,0.4);
    color: white;
    border-radius: 4px;
    padding: 0.2rem 0.7rem;
    cursor: pointer;
    font-size: 0.82rem;
  }
</style>
