/* Trace[expr] — user-facing builtin (beads-planning-b3k, Phase 3).
 *
 * A thin wrapper over the evaluator-side collector eval_collect_trace()
 * (src/eval.c). Trace[expr] returns a flat List of the successive top-level
 * expressions produced while evaluating expr to a fixed point:
 *   - >=1 top-level rewrite -> {e0, e1, ..., eN}
 *   - no top-level rewrite (inert atom / normal form) -> {}
 *
 * All the depth/clock/ownership subtlety lives in eval_collect_trace; this file
 * only handles arity dispatch, registration, attributes, and the docstring.
 * The 2-argument form Trace[expr, form] (form filtering) is deferred to
 * beads-planning-h4u: builtin_trace returns NULL for arities != 1 so those
 * calls stay unevaluated rather than silently misbehaving.
 */
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

/* Trace requires HoldAll: the argument must reach the builtin unevaluated so
 * the collector can observe its rewrite sequence from the start. Without Hold,
 * the evaluator would fully reduce expr before builtin_trace ever runs.
 *
 * Ownership: the evaluator owns `res` and frees it after a non-NULL return.
 * eval_collect_trace merely BORROWS args[0] (it inc-refs each recorded step),
 * so the returned List is independent of `res` — no NULL-out-before-free dance
 * is needed here. Returning NULL leaves the call unevaluated and the evaluator
 * retains `res`.
 *
 * HoldForm wrapping: builtin_trace's return value flows back through the
 * evaluator's fixed-point loop, which re-evaluates it. `List` is not held, so a
 * bare intermediate like `1 + 1` would reduce a second time (turning the trace
 * into {2, 2}). Wrapping each step in HoldForm makes the list inert on that
 * pass while printing transparently (HoldForm[1+1] displays as `1 + 1`),
 * yielding the Mathematica-style {1 + 1, 2}. The raw collector output is left
 * unwrapped so eval_collect_trace stays a clean data primitive for later
 * phases (TraceDepth, Trace[expr, form]). */
Expr* builtin_trace(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;   /* 2-arg form (h4u) and malformed calls stay unevaluated */

    Expr* raw = eval_collect_trace(res->data.function.args[0]);
    if (!raw) return NULL;

    for (size_t i = 0; i < raw->data.function.arg_count; i++) {
        Expr* step = raw->data.function.args[i];       /* ownership moves into HoldForm */
        Expr* wrap[1] = { step };
        raw->data.function.args[i] =
            expr_new_function(expr_new_symbol(SYM_HoldForm), wrap, 1);
    }
    return raw;
}

void trace_init(void) {
    symtab_add_builtin("Trace", builtin_trace);
    symtab_get_def("Trace")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Trace",
        "Trace[expr]\n\tGenerates a list of the successive top-level expressions\n"
        "\tproduced while evaluating expr. Returns {} when expr needs no\n"
        "\trewriting. Sub-evaluations of arguments are not recorded (flat form).");
}
