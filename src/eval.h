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

/* Trace[expr] support (beads-planning-b3k). Evaluate `held_expr` to a fixed
 * point, recording each top-level rewrite, and return a freshly-owned flat
 * List of the intermediate expressions (Trace v1 semantics: {e0, ..., eN}, or
 * {} when no top-level rewrite fires). `held_expr` is BORROWED (not consumed).
 * Reentrant: safe to call from inside a running evaluation, including nested
 * within another trace. */
Expr* eval_collect_trace(Expr* held_expr);

/* True iff `e` is an in-flight Throw sentinel: Throw[v], Throw[v,tag] or
 * Throw[v,tag,f] (head SYM_Throw, arity 1-3). Catch/Throw propagate by
 * returning this expression up through the normal evaluation return paths;
 * evaluate_step's argument loop short-circuits on it so a Throw fired deep
 * inside f[...]/Plus[...]/etc. reaches the nearest enclosing Catch. A held
 * Throw (e.g. Hold[Throw[1]]) is never evaluated and so never becomes a
 * sentinel, which is why reusing the plain Throw[...] node is unambiguous. */
bool eval_is_inflight_throw(const Expr* e);
/* In-flight Abort[] sentinel; and "a Throw or an Abort" -- what every site
 * that passes a Throw upward must pass upward too. */
bool eval_is_inflight_abort(const Expr* e);
bool eval_is_inflight_unwind(const Expr* e);
/* True (once) if an abort reached the top level since the last call. */
bool eval_take_top_abort(void);
/* Evaluate an assignment LHS's arguments as Set does (patterns and Hold
 * positions held, compound heads f[a][b] too).  Returns a new node with
 * *freed = true (caller frees), or `lhs` itself with *freed = false. */
Expr* eval_assignment_target(Expr* lhs, bool* freed);
/* lhs :> body /; t1 /; t2 as the pattern Condition[Condition[lhs, t1], t2]
 * with replacement body (the guards then run during matching, so the matcher
 * can backtrack into other assignments). Returns lhs itself (borrowed) when
 * rhs carries no Condition, else a new pattern the caller frees; *body is set
 * to the replacement to use (borrowed from rhs). */
Expr* eval_hoist_rhs_conditions(Expr* lhs, Expr* rhs, Expr** body);
/* A Catch / CheckAbort that consumes an in-flight sentinel calls this so the
 * enclosing frame does not re-raise it (see g_pending_unwind in eval.c). */
void eval_unwind_consumed(void);
/* Evaluate a user callback (a predicate, comparator, key function, pattern
 * test or guard) for a builtin that only inspects its value: a Throw/Abort it
 * produces is parked so the builtin's frame propagates it even if the builtin
 * drops the value.  eval_park_unwind parks an already-evaluated result. */
Expr* eval_callback(Expr* call);
void eval_park_unwind(const Expr* r);

/* True iff `e` is an in-flight Goto sentinel: Goto[tag] (head SYM_Goto, arity
 * 1). Like Throw, a Goto propagates by being returned up the evaluation stack;
 * evaluate_step's argument loop short-circuits on it so a Goto fired inside a
 * nested call (e.g. an If branch) bubbles up to the enclosing
 * CompoundExpression, which scans its statements for the matching Label[tag].
 * Unlike Throw it is consumed at a CompoundExpression, not a Catch. */
bool eval_is_inflight_goto(const Expr* e);

/* True iff `e` is an uncaught Break[]/Continue[] marker: head SYM_Break or
 * SYM_Continue with zero arguments. Unlike Throw/Goto, these do NOT short-
 * circuit the argument-evaluation loop -- they propagate only through
 * CompoundExpression, If, and the Do/For/While loop boundaries that inspect
 * their evaluated result's head. This predicate is used at the top level of
 * evaluate() to detect a marker that escaped every enclosing loop, so it can
 * be reported (Break::nofwd / Continue::nofwd) and rendered inert. */
bool eval_is_inflight_break_continue(const Expr* e);

/* Finish an in-flight sentinel that reached a top level, as evaluate() does
 * at depth 0: an uncaught Throw is reported (Throw::nocatch) and becomes
 * Hold[Throw[...]] (or f[v, tag]), Abort[] becomes $Aborted (and is noted for
 * eval_take_top_abort), Goto / Break / Continue are reported and made inert.
 * Anything else is returned unchanged. Takes ownership of `e`. For a nested
 * top level -- a Dialog input line -- which evaluate() runs at depth > 0. */
Expr* eval_report_uncaught(Expr* e);

// Recursion-depth limit for nested evaluate() calls. Guards C-stack overflow
// when expressions trigger deeply recursive sub-evaluation (Listable, Flat,
// user DownValues that re-evaluate, etc.). Default 1024 (matches Mathematica).
int  eval_get_recursion_limit(void);
void eval_set_recursion_limit(int n);
int  eval_get_recursion_depth(void);

/* Force the internal recursion-depth counter back to `n`, clearing the
 * sticky overflow flag.  Intended ONLY for siglongjmp-based unwinds
 * (TimeConstrained, pmint timeout) where evaluate()'s normal
 * "decrement on return" did not get a chance to run.  Calling this
 * during normal evaluation will corrupt the recursion guard. */
void eval_reset_recursion_depth(int n);

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
 * plans/EVAL_IMPROVEMENTS_PLAN.md: any user definition change invalidates
 * every cached evaluation in one shot, which is correct in all cases
 * but coarser than a per-symbol dependency graph would be. The clock
 * starts at 1 so that fresh nodes (whose `last_evaluated_at` is 0)
 * never accidentally hit the cache on first evaluation. */
uint64_t eval_clock_get(void);
void     eval_clock_bump(void);
/* Bump the eval clock for an OwnValue binding of `def`: INERT fixed points
 * (eval.c) stay valid unless one of them relied on `def` being unbound. */
struct SymbolDef;
void     eval_clock_bump_binding(struct SymbolDef* def);

/* Id of the current OUTERMOST evaluate() call: bumped once per user command /
 * script statement, and stable across all the nested evaluations and
 * fixed-point iterations reducing that command. Lets a builtin scope
 * per-command state (see the Integrate fail-memo) without an explicit reset. */
uint64_t eval_toplevel_id(void);

/* Ground fixed-point epoch (see the long comment in eval.c). The clock above
 * is coarse -- every mutation invalidates every cached fixed point, so an
 * iterator rebinding re-canonicalises loop-invariant values O(size) per step.
 * `g_last_rule_change_clock` is a finer epoch advanced ONLY by mutations that
 * change how a head evaluates (DownValue add, attribute/Protect change, Clear/
 * Remove). A node built solely from literals under the six pure structural
 * constructors (List/Association/Rule/RuleDelayed/Complex/Rational) carries a
 * GROUND flag in the top bit of `last_evaluated_at` and re-validates as a fixed
 * point whenever its stamp >= this epoch -- surviving ordinary OwnValue churn.
 *   - eval_rule_epoch_get():  read the epoch.
 *   - eval_rule_epoch_mark(): record the CURRENT clock as a rule change (caller
 *                             has already bumped the clock, e.g. via add_rule).
 *   - eval_rule_epoch_bump(): bump the clock AND mark it (attribute sites). */
uint64_t eval_rule_epoch_get(void);
void     eval_rule_epoch_mark(void);
void     eval_rule_epoch_bump(void);

/* Mask-aware views of `last_evaluated_at` for tests/tools (the raw field now
 * carries the GROUND flag in its top bit; read it through these, not directly).
 *   - eval_node_stamp():     the fixed-point stamp, flag masked off.
 *   - eval_node_is_ground(): true iff the GROUND flag is set. */
uint64_t eval_node_stamp(const Expr* e);
/* Would the evaluator reuse `e` as a fixed point right now (stamped under the
 * live clock, or a still-valid GROUND/INERT node)?  Atoms and symbols: true. */
bool     eval_node_reusable(const Expr* e);
/* Stamp `e` as a fixed point under the live clock, with its GROUND/INERT flag.
 * For a builtin that has PROVED its fresh result is one, so the evaluator does
 * not run the builtin a second time only to confirm nothing changes. */
void     eval_mark_fixed_point(Expr* e);
bool     eval_node_is_ground(const Expr* e);

/* For a builtin returning a freshly built tree of Lists over literal atoms
 * (Integer/Real/BigInt/String): stamp it as an evaluated GROUND fixed point so
 * the evaluator does not re-walk it. A no-op (returns false) when List carries
 * rules; sub-trees that are not pure literal Lists are left unstamped. Borrows
 * `e`; it must not yet be shared with an expression that was evaluated. */
bool     eval_stamp_list_literal(Expr* e);

// Helper to evaluate and free the input expression
static inline Expr* eval_and_free(Expr* e) {
    if (!e) return NULL;
    Expr* res = evaluate(e);
    expr_free(e);
    return res;
}

// Attributes built-in
Expr* builtin_attributes(Expr* res);

/* ------------------------------------------------------------------------
 * Return[] flow-control classifier (Withoff §3.1, EVAL_IMPROVEMENTS_PLAN
 * §4.3).
 *
 * Mathematica's Return is a marker that propagates outward through the
 * current evaluation tree until it hits a "scope boundary" -- the body
 * of a Function, a Module/Block/With, or a Do/For/While loop. The
 * boundary strips the marker and yields the wrapped value as its own
 * result. Return takes effect as soon as it is evaluated, even inside
 * other functions (CompoundExpression and other Hold-free heads pass
 * it through).
 *
 * Three argument forms:
 *
 *   Return[]        -- yield Null from the enclosing boundary.
 *   Return[expr]    -- yield expr from the enclosing boundary.
 *   Return[expr, h] -- yield expr from the *nearest enclosing* boundary
 *                      whose head is the symbol h. Boundaries with a
 *                      different head must propagate the Return outward
 *                      unchanged so that h can be reached.
 *
 * `eval_classify_return` is the single decision point used by every
 * scope-boundary builtin. It inspects `e` once and reports whether the
 * caller should:
 *
 *   EVAL_RETURN_NONE       -- `e` is not a Return marker. Treat as a
 *                              normal value (no behavioural change).
 *   EVAL_RETURN_CONSUME    -- `e` is a Return targeting *this*
 *                              boundary. *out_value is set to a freshly
 *                              owned Expr (Null for Return[], an
 *                              expr_copy of the payload otherwise) that
 *                              the caller should yield as its own
 *                              result. The caller still owns `e` and
 *                              must free it.
 *   EVAL_RETURN_PROPAGATE  -- `e` is Return[v, h] but `h` does not match
 *                              `boundary_head`. The caller should hand
 *                              `e` upward as its own result without
 *                              freeing it (the next boundary outward
 *                              gets to classify it).
 *
 * `boundary_head` is the interned canonical pointer of the enclosing
 * construct's head (e.g. SYM_Module, SYM_Function, SYM_Do). It may be
 * NULL, in which case 2-arg Return is always treated as PROPAGATE.
 *
 * The classifier never frees `e` and never bumps the eval clock. */
typedef enum {
    EVAL_RETURN_NONE,
    EVAL_RETURN_CONSUME,
    EVAL_RETURN_PROPAGATE
} EvalReturnAction;

EvalReturnAction eval_classify_return(Expr* e,
                                      const char* boundary_head,
                                      Expr** out_value);

/* Flatten nested same-head sub-expressions in place. Returns true iff
 * the args list was actually rewritten; callers in the §3.4 fixed-point
 * detector use the return value to decide whether to count this as a
 * change in the current evaluation step. */
bool eval_flatten_args(Expr* e, const char* head_name);
int eval_compare_expr_ptrs(const void* a, const void* b);
/* Canonical Orderless sort with per-sort symbol-set memoization (sort.c). */
void expr_orderless_sort(Expr** args, size_t n);
#if UP_VALUES
/* True when an upvalue of `tag` with LHS `lhs` is exactly what `lhs ^= rhs`
 * would create -- tag is the only operand tag -- so Definition may render it in
 * UpSet form rather than as `tag /: lhs = rhs`. */
bool upvalue_renders_as_upset(const char* tag, Expr* lhs);
#endif
static inline void eval_sort_args(Expr* e) { if (e->type == EXPR_FUNCTION && e->data.function.arg_count > 0) expr_orderless_sort(e->data.function.args, e->data.function.arg_count); }
#endif // EVAL_H
