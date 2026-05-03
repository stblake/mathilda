#ifndef EVAL_H
#define EVAL_H

#include "expr.h"
#include "attr.h"
#include <stdint.h>

/* Single-step evaluation (one level of rewrite).
 *
 * `changed` is an out-parameter: on return, *changed is true iff a real
 * rewrite fired during this step (head re-evaluation produced a different
 * subtree, an arg was reduced, Sequence/Flat/Orderless reshaped the args,
 * a DownValue / built-in / special primitive matched, etc.). When it is
 * false, the returned Expr* is structurally identical to the input and
 * the outer fixed-point loop can stop without doing an O(tree) expr_eq
 * compare. May be NULL for callers that do not care.
 *
 * False-positives (`*changed=true` but the result happens to match
 * the input) are correctness-safe — the outer loop just runs one extra
 * iteration. False-negatives are a bug: the loop would terminate before
 * reaching a fixed point. Sub-step instrumentation in evaluate_step is
 * therefore conservative — when in doubt, mark the step as a change. */
Expr* evaluate_step(Expr* e, bool* changed);

// Full infinite evaluation loop (stops at fixed point)
Expr* evaluate(Expr* e);

// Recursion-depth limit for nested evaluate() calls. Guards C-stack overflow
// when expressions trigger deeply recursive sub-evaluation (Listable, Flat,
// user DownValues that re-evaluate, etc.). Default 1024 (matches Mathematica).
int  eval_get_recursion_limit(void);
void eval_set_recursion_limit(int n);
int  eval_get_recursion_depth(void);

// Seeds the user-visible $RecursionLimit symbol with its default value
// and installs the docstring. Call after symtab_init().
void eval_init(void);

/* M3 phase-3: global evaluation clock.
 *
 * Monotonic counter incremented whenever the symbol table is mutated in a
 * way that could change the meaning of a previously-evaluated expression
 * (Set/SetDelayed adding an OwnValue or DownValue, Clear / ClearAll
 * removing them, SetAttributes / ClearAttributes adjusting attributes).
 * Pure builtin invocations and docstring writes do NOT bump the clock.
 *
 * The evaluator uses `eval_clock_get()` at the entry of `evaluate()` to
 * decide whether a previously-stamped result is still valid: if
 * `e->last_evaluated_at == eval_clock_get()`, the expression has been
 * fully reduced under the current symbol-table state and `evaluate()`
 * returns an inc-ref'd view immediately, skipping the entire fixed-point
 * loop and the body of `evaluate_step`.
 *
 * This is the conservative single-counter variant from §3.3 of
 * EVAL_IMPROVEMENTS_PLAN.md: any user definition change invalidates
 * every cached evaluation in one shot, which is correct in all cases
 * but coarser than a per-symbol dependency graph would be. The clock
 * starts at 1 so that fresh nodes (whose `last_evaluated_at` is 0)
 * never accidentally hit the cache on first evaluation. */
uint64_t eval_clock_get(void);
void     eval_clock_bump(void);

// Helper to evaluate and free the input expression
static inline Expr* eval_and_free(Expr* e) {
    if (!e) return NULL;
    Expr* res = evaluate(e);
    expr_free(e);
    return res;
}

// Attributes built-in
Expr* builtin_attributes(Expr* res);

/* Flatten nested same-head sub-expressions in place. Returns true iff
 * the args list was actually rewritten; callers in the §3.4 fixed-point
 * detector use the return value to decide whether to count this as a
 * change in the current evaluation step. */
bool eval_flatten_args(Expr* e, const char* head_name);
int eval_compare_expr_ptrs(const void* a, const void* b);
static inline void eval_sort_args(Expr* e) { if (e->type == EXPR_FUNCTION && e->data.function.arg_count > 0) qsort(e->data.function.args, e->data.function.arg_count, sizeof(Expr*), eval_compare_expr_ptrs); }
#endif // EVAL_H
