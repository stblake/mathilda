#ifndef CORE_H
#define CORE_H

#include "expr.h"

// Core programming routines
Expr* builtin_length(Expr* res);

/* LeafCount-style recursive node count.  `heads = true` includes
 * function heads in the count (matching the Mathilda LeafCount[expr,
 * True] semantics); `heads = false` counts only argument leaves. */
int64_t leaf_count_internal(Expr* e, bool heads);
Expr* builtin_dimensions(Expr* res);
Expr* builtin_clear(Expr* res);
Expr* builtin_unset(Expr* res);
Expr* builtin_clear_all(Expr* res);
Expr* builtin_remove(Expr* res);
Expr* builtin_protect(Expr* res);
Expr* builtin_unprotect(Expr* res);
Expr* builtin_append(Expr* res);
Expr* builtin_prepend(Expr* res);
Expr* builtin_append_to(Expr* res);
Expr* builtin_prepend_to(Expr* res);
Expr* builtin_own_values(Expr* res);
Expr* builtin_down_values(Expr* res);
Expr* builtin_out(Expr* res);
Expr* builtin_compoundexpression(Expr* res);
Expr* builtin_atomq(Expr* res);
Expr* builtin_identity(Expr* res);
Expr* builtin_composition(Expr* res);
Expr* builtin_compose_list(Expr* res);
Expr* builtin_numberq(Expr* res);
Expr* builtin_numericq(Expr* res);
Expr* builtin_positive(Expr* res);
Expr* builtin_negative(Expr* res);
Expr* builtin_nonnegative(Expr* res);
Expr* builtin_nonpositive(Expr* res);
Expr* builtin_integerq(Expr* res);
Expr* builtin_evenq(Expr* res);
Expr* builtin_oddq(Expr* res);
Expr* builtin_mod(Expr* res);
Expr* builtin_quotient(Expr* res);
Expr* builtin_quotientremainder(Expr* res);
Expr* builtin_level(Expr* res);
Expr* builtin_depth(Expr* res);
Expr* builtin_leafcount(Expr* res);
Expr* builtin_bytecount(Expr* res);
Expr* builtin_information(Expr* res);
Expr* builtin_evaluate(Expr* res);
Expr* builtin_chop(Expr* res);
Expr* builtin_clip(Expr* res);
Expr* builtin_releasehold(Expr* res);
Expr* builtin_tostring(Expr* res);
Expr* builtin_toexpression(Expr* res);
Expr* builtin_symbol(Expr* res);
Expr* builtin_increment(Expr* res);
Expr* builtin_decrement(Expr* res);
Expr* builtin_preincrement(Expr* res);
Expr* builtin_predecrement(Expr* res);
Expr* builtin_addto(Expr* res);
Expr* builtin_subtractfrom(Expr* res);
Expr* builtin_time_constrained(Expr* res);

/* Cooperative wall-clock deadline check used by TimeConstrained.
 *
 * Called by evaluate()'s fixed-point loop once per rewrite iteration.
 * No-op outside any TimeConstrained call. When the innermost
 * TimeConstrained's wall-clock deadline has passed, sets the timeout
 * flag and siglongjmp's back to the matching sigsetjmp -- identical to
 * the SIGPROF path. This is the portable backstop for hosts whose
 * ITIMER_PROF is unreliable (notably WSL 1, which under-counts CPU time
 * and delivers SIGPROF many times late). On real Linux / macOS the
 * SIGPROF normally fires first and this check is a cheap no-op.
 *
 * Granularity: limited to between rewrite steps -- a single long-running
 * C builtin (e.g. FactorInteger on a large composite) cannot be aborted
 * cooperatively. On those hosts, SIGPROF is the only mechanism, with
 * its usual portability caveats. */
void tc_check_deadline(void);

// Initialize core builtins
void core_init(void);

#endif // CORE_H
