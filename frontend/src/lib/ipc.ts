// ipc.ts — typed wrappers around Tauri invoke/Channel

import { invoke } from '@tauri-apps/api/core';
import { Channel } from '@tauri-apps/api/core';
import type { CellType } from './notebook';

export type OutputMessage =
  | { id: number; type: 'expr';   payload: string; latex?: string }
  /* Which kernel line this statement took, sent before its output.
     NOT decorative and not something the front end can count for itself: `%`,
     `%%`, `%3`, `In[3]` and `Out[3]` all resolve against the KERNEL's $Line, and
     one kernel serves every notebook on the canvas, so a per-notebook counter
     would label a cell In[2] while `%2` addressed a line from another notebook.
     Arrives even for a statement that produces no result (`x = 5;`), which is
     why it is its own message rather than a field on the result. */
  | { id: number; type: 'line';   line: number }
  /* `symbol` is the name `?name` asked about, so the notebook can offer that symbol's
     documentation page. */
  | { id: number; type: 'usage';  payload: string; symbol?: string }
  | { id: number; type: 'names';  payload: string[] }
  | { id: number; type: 'error';  message: string }
  /* Print output of the statement being evaluated, sent before its result. */
  | { id: number; type: 'stream'; text: string }
  /* A kernel message (warning) such as `Power::infy: ...`, sent before the result. */
  | { id: number; type: 'message'; text: string }
  | { id: number; type: 'plot';   payload: object }
  /* A raster result: base64 RGBA plus its shape. A volume also carries `depth` and the 1-based
     `slice` it sent, which is the middle one. */
  | { id: number; type: 'image';  payload: { w: number; h: number; channels: number;
                                             data: string; depth?: number; slice?: number ; faces?: Record<string, { w: number; h: number; data: string }>} }
  | { id: number; type: 'html';   payload: string };

/** One cell as a `.mnb` file stores it: its type and source, no output.
 *  The Rust side normalises unknown types to 'code' (notebook_format.rs). */
export type CellData = {
  type: CellType;
  source: string;
};

/** Undo the typographic substitutions a WebView (or a paste from a word
 *  processor) may have applied, so the kernel always receives ASCII source:
 *  curly quotes → straight, en/em/minus dashes → hyphen, `→` → `->`. These
 *  characters have no distinct meaning in Mathilda input, and leaving them in
 *  would turn `->` into `–>` (a parse error) or break string delimiters. */
export function normalizeInput(expr: string): string {
  return expr
    .replace(/[‘’‛]/g, "'")     // ‘ ’ ‛ → '
    .replace(/[“”‟]/g, '"')     // “ ” ‟ → "
    .replace(/→/g, '->')                  // → (rightwards arrow) → ->
    .replace(/[–—−]/g, '-');    // – — − → -
}

/** Evaluate a Mathilda expression, calling `onMessage` for each
 *  streamed output until the kernel signals "done". */
export async function evaluateCell(
  expr: string,
  onMessage: (msg: OutputMessage) => void
): Promise<void> {
  const channel = new Channel<OutputMessage>();
  channel.onmessage = onMessage;
  await invoke<void>('evaluate_cell', { expr: normalizeInput(expr), channel });
}

export async function restartKernel(): Promise<void> {
  await invoke<void>('restart_kernel');
}

export async function interruptKernel(): Promise<void> {
  await invoke<void>('interrupt_kernel');
}

export async function pingKernel(): Promise<void> {
  await invoke<void>('ping_kernel');
}

export async function saveNotebook(path: string, cells: CellData[]): Promise<void> {
  await invoke<void>('save_notebook', { path, cells });
}

export async function loadNotebook(path: string): Promise<CellData[]> {
  return await invoke<CellData[]>('load_notebook', { path });
}

export async function saveLibrary(path: string, json: string): Promise<void> {
  await invoke<void>('save_library', { path, json });
}

export async function loadLibrary(path: string): Promise<string> {
  return await invoke<string>('load_library', { path });
}

export async function setWindowTitle(title: string): Promise<void> {
  await invoke<void>("set_window_title", { title });
}

/* File > Open Recent. The LIST lives in Rust (src-tauri/src/recent.rs), which also owns the
   native submenu and the JSON store under the app config dir — the menu is native, so only
   Rust can draw it, and it has to be drawn before the webview has booted.

   What stays here is the acting: Rust emits `menu:recent-<i>` like any other menu item, and
   App.svelte indexes `recentFiles()` with that i. So the order of this array IS the order of
   the menu items, and nothing else may sort it. */

export async function recentFiles(): Promise<string[]> {
  return await invoke<string[]>('recent_files');
}

export async function pushRecentFile(path: string): Promise<void> {
  await invoke<void>('push_recent_file', { path });
}

export async function forgetRecentFile(path: string): Promise<void> {
  await invoke<void>('forget_recent_file', { path });
}

export async function clearRecentFiles(): Promise<void> {
  await invoke<void>('clear_recent_files');
}

/** Open a URL in the user's real browser, not the app's webview.
 *
 * Reference pages link out to GitHub source and to the published site; letting
 * those navigate the webview would replace the notebook with a web page and
 * leave no way back. Imported lazily so the shell plugin is only pulled in if a
 * link is actually clicked. */
export async function openUrl(url: string): Promise<void> {
  const { open } = await import('@tauri-apps/plugin-shell');
  await open(url);
}
