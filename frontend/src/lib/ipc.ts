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

/** Evaluate a Mathilda expression, calling `onMessage` for each
 *  streamed output until the kernel signals "done". */
export async function evaluateCell(
  expr: string,
  onMessage: (msg: OutputMessage) => void
): Promise<void> {
  const channel = new Channel<OutputMessage>();
  channel.onmessage = onMessage;
  await invoke<void>('evaluate_cell', { expr, channel });
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
