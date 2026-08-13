// kernelActions.ts — kernel lifecycle, with the status transitions that go with it.
//
// ipc.ts is a thin invoke() wrapper and should stay that way; these two calls
// also own the kernelStatus store, which is a different concern. They live here
// rather than as functions inside App.svelte because the toolbar's kernel menu
// has to drive exactly the same paths as the native Kernel menu -- two
// implementations of "abort" is how one of them ends up wrong.

import { get } from 'svelte/store';
import { kernelStatus } from './notebook';
import { restartKernel, interruptKernel } from './ipc';

export async function restart() {
  kernelStatus.set('restarting');
  try { await restartKernel(); kernelStatus.set('ready'); }
  catch { kernelStatus.set('dead'); }
}

/** Abort a running evaluation.
 *
 *  `interrupt_kernel` does NOT send SIGINT despite what its comment says -- it
 *  calls kill() on the child process and does not respawn it
 *  (src-tauri/src/kernel.rs). Abort on its own therefore leaves a dead kernel
 *  and the next cell fails with no explanation. The restart is not optional.
 *
 *  Which means every control offering this must say so in its label: the user is
 *  trading the session's definitions for regaining control. A cooperative abort
 *  -- a flag checked inside the C evaluator, so the computation stops without
 *  losing state -- is the real fix, and it is not a frontend change. Until then
 *  this is the honest behaviour, plainly labelled. */
export async function abortEvaluation() {
  /* Nothing to abort, and nothing to lose: don't restart a healthy idle kernel
     just because the user pressed Cmd+. out of habit. */
  const status = get(kernelStatus);
  if (status !== 'busy' && status !== 'restarting') return;

  kernelStatus.set('restarting');
  try { await interruptKernel(); } catch { /* already gone; restart regardless */ }
  try { await restartKernel(); kernelStatus.set('ready'); }
  catch { kernelStatus.set('dead'); }
}
