/*
 * memacct.h -- kernel memory accounting, MemoryConstrained, and the constraint
 * frames MemoryConstrained and TimeConstrained share.
 *
 * The kernel keeps one NET count of the bytes it holds (g_mth_mem_live):
 * charged when an Expr node, a String payload, an NDArray buffer or a GMP/MPFR
 * limb block is allocated, released when it is freed.  MemoryConstrained[expr,
 * b] stacks a frame whose limit is (live at entry + b); the effective limit
 * g_mth_mem_limit is the smallest active one, and INT64_MAX when no frame is
 * active, so every hot-path test is a single compare.
 *
 * Net, not gross: memory a computation allocates and frees again does not count
 * against its budget (the reference behaviour: MemoryConstrained[Do[Range[10^4],
 * {1000}], 10^6] is Null).  What is counted is an approximation of the kernel's
 * footprint -- the node struct, string bytes, array buffers and big-number
 * limbs -- not of malloc's bookkeeping; FLINT's own allocations and worker-thread
 * allocations are not counted.
 *
 * How a constraint stops the evaluation.  When the live count passes a frame's
 * limit (detected by the evaluator's per-step poll, or up front by
 * mth_mem_request before a large buffer is allocated) that frame is TRIPPED.  A
 * trip becomes the ordinary Abort[] sentinel -- so CheckAbort catches it and
 * every enclosing construct unwinds through its normal cleanup -- once it is
 * ACTIONABLE: when no AbortProtect opened inside the tripped frame's body is
 * still running (an AbortProtect *enclosing* the MemoryConstrained does not
 * defer it).  The tripped frame's MemoryConstrained consumes the abort and
 * returns failexpr (or $Aborted); a plain Abort[] or a user interrupt is not
 * consumed and propagates.  All of this runs on the main thread only.
 *
 * TimeConstrained (core.c) pushes a TIME frame on the same stack: its trip is
 * recorded when the deadline passes and acted on exactly like a memory trip --
 * the cooperative path that unwinds through ordinary cleanup.  Its siglongjmp
 * is only the backstop for C code that never returns to the evaluator.  One
 * stack means one answer to "whose abort is this?" and one AbortProtect rule
 * for both constraints, nested in any order.
 */
#ifndef MATHILDA_MEMACCT_H
#define MATHILDA_MEMACCT_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "expr.h"

/* Net bytes held (may go negative within a scope: memory allocated before a
 * frame opened and freed inside it), and the effective limit. */
extern int64_t g_mth_mem_live;
extern int64_t g_mth_mem_limit;
/* Index of the outermost tripped frame, or -1. */
extern int     g_cf_tripped;
/* Nonzero while any MemoryConstrained frame is open.  Charging happens only
 * then: a budget is a DELTA from the live count at entry, so what is charged
 * outside every frame can never matter, and an idle kernel pays one predicted
 * branch per allocation instead of a read-modify-write of a global. */
extern int     g_mth_mem_active;

#define MTH_MEM_CHARGE(n)  do { if (g_mth_mem_active) g_mth_mem_live += (int64_t)(n); } while (0)
/* Charge a variable-size block (string, array buffer, limbs) and note an
 * overrun at once: such a block may be freed again before the next poll, and
 * the reference counts what was requested, not what survived to be looked at.
 * The per-node charge stays a bare add (nodes are caught by the poll). */
#define MTH_MEM_CHARGE_BLOCK(n) do { if (g_mth_mem_active) { g_mth_mem_live += (int64_t)(n); \
        if (g_mth_mem_live > g_mth_mem_limit) mc_note_exceeded(); } } while (0)
#define MTH_MEM_RELEASE(n) do { if (g_mth_mem_active) g_mth_mem_live -= (int64_t)(n); } while (0)

/* True when a constraint is tripped or over budget: the evaluator's poll calls
 * cf_poll() only then. */
#define MTH_CF_POLL_DUE() \
    (g_cf_tripped >= 0 || (g_mth_mem_active && g_mth_mem_live > g_mth_mem_limit))

/* Ask before allocating `bytes` for a large buffer.  True: go ahead.  False:
 * the allocation would break an active MemoryConstrained whose abort can act
 * now -- the frame is tripped and the caller must decline (free its locals and
 * return NULL / its failure value) WITHOUT falling back to another path that
 * allocates the same memory.  Inside an AbortProtect the answer is true (the
 * protected code runs to completion) and the trip is only recorded. */
bool mth_mem_request(size_t bytes);

/* mth_mem_request for an array of `elem_bytes` elements whose shape is `spec`:
 * a non-negative machine integer n or a List of them.  Any other spec answers
 * true (the producer decides what it means).  Free while no constraint is on. */
bool mth_mem_request_shape(const Expr* spec, size_t elem_bytes);

/* Record any frame now over its limit as tripped (sticky until consumed). */
void mc_note_exceeded(void);
/* A trip exists and is not deferred by an AbortProtect inside its frame. */
bool cf_abort_actionable(void);
/* The evaluator's poll: note an overrun, then report whether the evaluation
 * must turn into the Abort[] sentinel now. */
bool cf_poll(void);
/* Forget the trip (CheckAbort consumed the abort).  A count still over the
 * limit trips again at the next poll. */
void cf_clear_trip(void);
/* Re-instate a trip cf_clear_trip forgot (CheckAbort with PropagateAborts ->
 * True hands the SAME abort on), if its frame is still open. */
void cf_restore_trip(int frame);

/* Frame stack, for the builtins and for compiled code (OP_MC_ENTER/LEAVE).
 * The push functions return the new frame's index (its depth before the push),
 * or -1 when the stack is full (the caller then runs unconstrained). */
int  mc_frame_push(int64_t budget);
/* A TimeConstrained frame expiring at `deadline` (mth_mono_seconds). */
int  cf_time_push(double deadline);
/* Trip the outermost unspent TIME frame whose deadline is <= now. */
void cf_time_note(double now);
/* Earliest deadline of the unspent TIME frames; false when there is none. */
bool cf_time_deadline(double* out);
/* AbortProtect depth at entry of the outermost unspent TIME frame (-1: none):
 * the asynchronous backstop never unwinds code protected inside any of them. */
int  cf_time_protect_base(void);
/* The frame is a TIME frame. */
bool cf_is_time(int frame);
/* A caught-and-absorbed timeout (CheckAbort, PropagateAborts -> False): the
 * frame keeps running and never fires again, as in the reference. */
void cf_time_spend(int frame);
/* Pop frames until `depth` remain; a trip of a popped frame is forgotten. */
void cf_unwind_to(int depth);
int  cf_depth(void);

void  memacct_init(void);
Expr* builtin_memory_constrained(Expr* res);

#endif /* MATHILDA_MEMACCT_H */
