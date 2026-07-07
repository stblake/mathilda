// ipc.ts — typed wrappers around Tauri invoke/Channel

import { invoke } from '@tauri-apps/api/core';
import { Channel } from '@tauri-apps/api/core';

export type OutputMessage =
  | { id: number; type: 'expr';   payload: string }
  | { id: number; type: 'error';  message: string }
  | { id: number; type: 'stream'; text: string }
  | { id: number; type: 'plot';   payload: object }
  | { id: number; type: 'html';   payload: string };

export type CellData = {
  type: 'code' | 'prose';
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
