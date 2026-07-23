/* Trace[expr] / Trace[expr, form] — user-facing builtin
 * (beads-planning-b3k Phase 3; the two-arg form is beads-planning-h4u).
 *
 * A thin wrapper over the evaluator-side collector eval_collect_trace()
 * (src/eval.c). Trace[expr] returns a NESTED List that mirrors the structure of
 * expr's evaluation: each argument sub-evaluation that takes a step appears as a
 * sublist, and the reassembled intermediate form appears as a step (matching
 * Mathematica). E.g. Trace[2^3 + 4^2 + 1] -> {{2^3,8},{4^2,16},8+16+1,25}.
 *   - >=1 rewrite -> the nested list of forms
 *   - no rewrite (inert atom / normal form) -> {}
 *
 * Trace[expr, form] filters that nested trace to the step leaves whose
 * expression matches the pattern `form` (same structural matcher as MatchQ,
 * src/match.c), flattening the nesting away into a plain List of matches.
 * Because Trace is HoldAll, `form` reaches the builtin unevaluated, so pattern
 * literals like `_Integer` or `f[_]` work directly. A form that matches nothing
 * yields {}.
 *
 * All the nesting/clock/ownership subtlety lives in eval_collect_trace; this
 * file only handles arity dispatch, form filtering, HoldForm wrapping,
 * registration, attributes, and the docstring. Arities other than 1 or 2 return
 * NULL so the call stays unevaluated rather than silently misbehaving.
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
/* True for a List[...] node -- the collector's nesting marker (each argument
 * sub-evaluation that took a step is spliced in as a List). Step forms are any
 * other expression. */
static bool trace_is_list(Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* Recursively collect every step leaf of the nested trace `node` that
 * structurally matches `form` (borrowed) into the growable buffer (*buf,*n,*cap).
 * List nodes are descended into (structure, not steps); matching leaves are
 * copied in. Mathematica's Trace[expr, form] returns the matched forms; we
 * flatten the nesting away, consistent with the earlier flat semantics the
 * tests assert. */
static void trace_collect_matches(Expr* node, Expr* form,
                                  Expr*** buf, size_t* n, size_t* cap) {
    if (trace_is_list(node)) {
        for (size_t i = 0; i < node->data.function.arg_count; i++)
            trace_collect_matches(node->data.function.args[i], form, buf, n, cap);
        return;
    }
    MatchEnv* env = env_new();
    bool matched = match(node, form, env);
    env_free(env);
    if (!matched) return;
    if (*n == *cap) {
        size_t newcap = *cap ? *cap * 2 : 8;
        Expr** grown = realloc(*buf, newcap * sizeof(Expr*));
        if (!grown) return;                  /* OOM: drop this match */
        *buf = grown;
        *cap = newcap;
    }
    (*buf)[(*n)++] = expr_copy(node);         /* keep a copy; caller frees raw */
}

/* Wrap the step leaves of the nested trace `e` (owned) in HoldForm so the
 * returned structure stays inert under the evaluator's fixed-point re-pass while
 * still printing transparently, and returns the transformed tree (owning). List
 * nodes are recursed into and preserved (nesting); every non-List leaf is
 * wrapped. A leaf that is itself a value-List is recursed too -- wrapping its
 * elements prints identically to wrapping the list whole, so no step/nesting
 * disambiguation is needed. */
static Expr* trace_holdform_tree(Expr* e) {
    if (trace_is_list(e)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            e->data.function.args[i] = trace_holdform_tree(e->data.function.args[i]);
        return e;
    }
    Expr* wrap[1] = { e };
    return expr_new_function(expr_new_symbol(SYM_HoldForm), wrap, 1);
}

Expr* builtin_trace(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;  /* other arities stay unevaluated */

    Expr* raw = eval_collect_trace(res->data.function.args[0]);
    if (!raw) return NULL;

    /* Trace[expr, form]: keep only the step leaves matching the (held) pattern
     * form, flattening the nesting into a plain List of matches. Filtering runs
     * on the raw, pre-HoldForm steps so `form` matches the plain expression a
     * user wrote (e.g. _Integer against 5, not against HoldForm[5]). The
     * collector has already restored the outer trace state, so match()'s
     * internal evaluation (for /; and ?test) is safe here. */
    if (argc == 2) {
        Expr** buf = NULL; size_t n = 0, cap = 0;
        trace_collect_matches(raw, res->data.function.args[1], &buf, &n, &cap);
        expr_free(raw);
        raw = expr_new_function(expr_new_symbol(SYM_List), buf, n);
        free(buf);
    }

    return trace_holdform_tree(raw);
}

void trace_init(void) {
    symtab_add_builtin("Trace", builtin_trace);
    symtab_get_def("Trace")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Trace",
        "Trace[expr]\n\tGenerates a nested list of the expressions produced while\n"
        "\tevaluating expr. Each argument sub-evaluation that takes a step appears\n"
        "\tas a sublist, mirroring the structure of the evaluation. Returns {} when\n"
        "\texpr needs no rewriting.\n"
        "Trace[expr, form]\n\tIncludes only the steps whose expression matches the\n"
        "\tpattern form (e.g. Trace[expr, _Integer]).");
}
