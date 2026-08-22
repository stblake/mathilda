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
/** The kernel's resident memory in bytes, as of its last `done`. Null until the first
 *  evaluation: the number can only change because something was evaluated, so there is nothing
 *  honest to show before one has been. */
export const kernelMemory = writable<number | null>(null);

export const showStatusBar = writable(true);

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

/** Bytes for a status bar: three significant figures and a binary unit, so a growing kernel reads
 *  as a growing number rather than as a wall of digits. */
export function formatBytes(b: number): string {
  if (b < 1024) return `${b} B`;
  const kb = b / 1024;
  if (kb < 1024) return `${kb.toFixed(0)} KB`;
  const mb = kb / 1024;
  if (mb < 1024) return `${mb < 10 ? mb.toFixed(1) : mb.toFixed(0)} MB`;
  return `${(mb / 1024).toFixed(2)} GB`;
}
