/* Trace[expr] / Trace[expr, form] — user-facing builtin
 * (beads-planning-b3k Phase 3; the two-arg form is beads-planning-h4u).
 *
 * A thin wrapper over the evaluator-side collector eval_collect_trace()
 * (src/eval.c). Trace[expr] returns a flat List of the successive top-level
 * expressions produced while evaluating expr to a fixed point:
 *   - >=1 top-level rewrite -> {e0, e1, ..., eN}
 *   - no top-level rewrite (inert atom / normal form) -> {}
 *
 * Trace[expr, form] additionally filters that flat list to the steps whose
 * expression matches the pattern `form`, using the same structural matcher as
 * MatchQ (src/match.c). Because Trace is HoldAll, `form` reaches the builtin
 * unevaluated, so pattern literals like `_Integer` or `f[_]` work directly.
 * This is the flat, top-level analogue of Mathematica's Trace[expr, form]:
 * only the top-level rewrite sequence (never argument sub-evaluations) is
 * produced, then filtered -- consistent with the v1 flat semantics of the
 * one-arg form. A form that matches nothing yields {}.
 *
 * All the depth/clock/ownership subtlety lives in eval_collect_trace; this file
 * only handles arity dispatch, form filtering, HoldForm wrapping, registration,
 * attributes, and the docstring. Arities other than 1 or 2 return NULL so the
 * call stays unevaluated rather than silently misbehaving.
 */
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "match.h"

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
/* Filter `raw` (a flat List owned by the caller) in place, keeping only the
 * steps that structurally match `form` (borrowed). Dropped steps are freed;
 * kept steps are compacted to the front and arg_count is shrunk. The slot array
 * itself is untouched (still owned by `raw`), so no reallocation is needed and
 * expr_free(raw) later reclaims it. */
static void trace_filter_by_form(Expr* raw, Expr* form) {
    size_t keep = 0;
    for (size_t i = 0; i < raw->data.function.arg_count; i++) {
        Expr* step = raw->data.function.args[i];
        MatchEnv* env = env_new();
        bool matched = match(step, form, env);
        env_free(env);
        if (matched) {
            raw->data.function.args[keep++] = step;   /* compact toward front */
        } else {
            expr_free(step);                          /* drop */
        }
    }
    raw->data.function.arg_count = keep;
}

/* Wrap every element of `raw` (a flat List, owned) in HoldForm, in place, so
 * the returned list stays inert under the evaluator's fixed-point re-pass while
 * still printing transparently. Ownership of each step moves into its HoldForm
 * wrapper; the wrapper takes the slot. */
static void trace_wrap_holdform(Expr* raw) {
    for (size_t i = 0; i < raw->data.function.arg_count; i++) {
        Expr* step = raw->data.function.args[i];
        Expr* wrap[1] = { step };
        raw->data.function.args[i] =
            expr_new_function(expr_new_symbol(SYM_HoldForm), wrap, 1);
    }
}

Expr* builtin_trace(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;  /* other arities stay unevaluated */

    Expr* raw = eval_collect_trace(res->data.function.args[0]);
    if (!raw) return NULL;

    /* Trace[expr, form]: drop steps that don't match the (held) pattern form.
     * Filtering runs on the raw, pre-HoldForm steps so `form` matches the plain
     * expression a user wrote (e.g. _Integer against 5, not against
     * HoldForm[5]). The collector has already restored the outer trace state,
     * so match()'s internal evaluation (for /; and ?test) is safe here. */
    if (argc == 2)
        trace_filter_by_form(raw, res->data.function.args[1]);

    trace_wrap_holdform(raw);
    return raw;
}

void trace_init(void) {
    symtab_add_builtin("Trace", builtin_trace);
    symtab_get_def("Trace")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Trace",
        "Trace[expr]\n\tGenerates a list of the successive top-level expressions\n"
        "\tproduced while evaluating expr. Returns {} when expr needs no\n"
        "\trewriting. Sub-evaluations of arguments are not recorded (flat form).\n"
        "Trace[expr, form]\n\tIncludes only the steps whose expression matches the\n"
        "\tpattern form (e.g. Trace[expr, _Integer]).");
}
