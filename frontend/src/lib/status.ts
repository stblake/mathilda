// status.ts — what the status bar reports.
//
// Evaluation timing did not exist anywhere before this: a cell recorded its
// exec index and its output, but nothing measured how long the kernel took. The
// bar is the first consumer, so the measurement lives here rather than being
// buried in the card that happens to run the loop.

import { writable } from 'svelte/store';

export type LastOp = {
  /** How the operation is named in the bar, e.g. "In[3]". */
  label: string;
  /** Wall-clock milliseconds from request to last output. */
  ms: number;
  ok: boolean;
  /** First line of the expression, for context. Trimmed by the bar. */
  source: string;
};

/** The most recently completed evaluation, or null before the first one. */
export const lastOp = writable<LastOp | null>(null);

/** Evaluations completed this session. Reset when the kernel restarts. */
export const evalCount = writable(0);

/** Cumulative kernel time this session, in milliseconds. */
export const evalTotalMs = writable(0);

/** Is the status bar shown? Optional by request — some people want the window
 *  to be nothing but notebook. */
export const showStatusBar = writable(true);

/* ---- Kernel memory -------------------------------------------------------
 *
 * The kernel's resident bytes, from MemoryInUse[]. Off by default: it costs a
 * kernel round-trip every few seconds, which nobody should pay for without
 * asking, and most sessions never care.
 *
 * WHY THIS POLLS ONLY WHEN THE KERNEL IS IDLE, which is the whole design.
 * The C kernel is single-threaded and serializes every evaluation behind one
 * mutex. A poll issued while a long computation is running would sit in that
 * queue, and then the user's NEXT cell would sit behind the poll — so a status
 * bar, a decoration, would be adding latency to real work. That trade is never
 * worth it, so the poll is skipped whenever the kernel is not `ready`.
 *
 * The cost of that choice is real and worth naming: a big computation is
 * exactly when memory is most interesting, and it is exactly when this cannot
 * look. So the last known figure stays on screen and is marked STALE rather
 * than being blanked — a number with a caveat beats no number, and blanking it
 * would read as "memory dropped to nothing" at the very moment it spiked.
 *
 * Reading it without the kernel's help would need the sidecar's RSS read from
 * Rust, which works while busy but takes a different implementation per
 * platform and per target (a child PID on desktop, self on the in-process
 * mobile build). That is the better answer if the staleness ever becomes
 * annoying; it is not worth it for a first cut.
 */

/** Kernel resident bytes, or null before the first successful sample. */
export const memoryBytes = writable<number | null>(null);

/** True when the displayed figure predates the current kernel activity — set
 *  when a sample is skipped because the kernel is busy, or when one fails. */
export const memoryStale = writable(false);

/** Is the memory segment shown? Off by default; polling is not free. */
export const showMemory = writable(false);

/** Bytes as something that fits a 22px strip.
 *
 * Binary units, because that is what MemoryInUse counts and what Activity
 * Monitor and top show, so the three agree. One decimal below 100 and none
 * above keeps the width steady as the value moves, which matters in a strip
 * where a jittering number is more distracting than an imprecise one. */
export function formatBytes(b: number): string {
  if (!isFinite(b) || b < 0) return '—';
  if (b < 1024) return `${Math.round(b)} B`;
  const units = ['KB', 'MB', 'GB', 'TB'];
  let v = b / 1024;
  let i = 0;
  while (v >= 1024 && i < units.length - 1) { v /= 1024; i++; }
  return `${v >= 100 ? Math.round(v) : v.toFixed(1)} ${units[i]}`;
}

/** A kernel restart invalidates the figure along with the session totals: the
 *  process that used those bytes is gone. */
export function resetMemory() {
  memoryBytes.set(null);
  memoryStale.set(false);
}

export function recordOp(op: LastOp) {
  lastOp.set(op);
  if (op.ok) {
    evalCount.update(n => n + 1);
    evalTotalMs.update(t => t + op.ms);
  }
}

/** Called on kernel restart/abort: the session's totals no longer describe the
 *  kernel that is now running. */
export function resetSessionStats() {
  evalCount.set(0);
  evalTotalMs.set(0);
  lastOp.set(null);
  resetMemory();
}

/** Milliseconds as something short enough for a status bar. */
export function formatMs(ms: number): string {
  if (ms < 1) return '<1 ms';
  if (ms < 1000) return `${Math.round(ms)} ms`;
  if (ms < 60_000) return `${(ms / 1000).toFixed(2)} s`;
  const m = Math.floor(ms / 60_000);
  const s = Math.round((ms % 60_000) / 1000);
  return `${m}m ${s}s`;
}
