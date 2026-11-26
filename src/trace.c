/* Trace[expr] / Trace[expr, form] — user-facing builtin
 * (beads-planning-b3k Phase 3; the two-arg form is beads-planning-h4u).
 *
 * A thin wrapper over the evaluator-side collector eval_collect_trace()
 * (src/eval.c). Trace[expr] returns a NESTED List that mirrors the structure of
 * expr's evaluation: each argument sub-evaluation that takes a step appears as a
 * sublist, and the reassembled intermediate form appears as a step (as in the
 * reference). E.g. Trace[2^3 + 4^2 + 1] -> {{2^3,8},{4^2,16},8+16+1,25}.
 *   - >=1 rewrite -> the nested list of forms
 *   - no rewrite (inert atom / normal form) -> {}
 *
 * Trace[expr, form] filters that nested trace to the steps whose expression
 * matches the pattern `form` (same structural matcher as MatchQ, src/match.c),
 * keeping the nesting and dropping the sublists left empty; Trace[expr, s]
 * keeps the evaluations that use rules for the symbol s, and
 * Trace[expr, form -> rhs] gives rhs for each matching step. Because Trace is
 * HoldAll, `form` reaches the builtin unevaluated, so pattern literals like
 * `_Integer` or `f[_]` work directly. A form that matches nothing yields {}.
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
#include <stdlib.h>

/* The collector runs in marked mode (eval_collect_trace_marked): every step
 * is HoldForm[e], and every List node is nesting -- one sublist per argument
 * or sub-evaluation that took a step. A step whose value is itself a List
 * ({1, 2, 3} from Range[3]) is therefore a HoldForm, never mistaken for
 * nesting. */
static bool trace_is_nest(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* The expression a step records (the body of its HoldForm). */
static Expr* trace_step_body(Expr* e) {
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 1 &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_HoldForm)
        return e->data.function.args[0];
    return e;
}

/* The filter of Trace[expr, form]:
 *   TF_PATTERN  form is a pattern: keep the steps that match it;
 *   TF_RULE     form is lhs -> rhs (or :>): keep rhs, with lhs's bindings, for
 *               each step that matches lhs;
 *   TF_SYMBOL   form is a symbol s: keep the evaluations that use rules for
 *               s -- each step whose head is s (or that is s itself) and the
 *               step it was rewritten to.
 * The nesting of the trace is kept, and a sublist left empty is dropped, so
 * Trace[fac[3], fac[_]] is {fac[3], {fac[2], {fac[1]}}}. */
typedef enum { TF_PATTERN, TF_RULE, TF_SYMBOL } TraceFilterKind;

typedef struct {
    TraceFilterKind kind;
    Expr* pat;   /* the pattern (lhs for a rule); borrowed */
    Expr* rhs;   /* TF_RULE only; borrowed */
    const char* sym;  /* TF_SYMBOL only */
} TraceFilter;

static bool trace_uses_symbol(const Expr* e, const char* sym) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == sym;
    if (e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    while (h->type == EXPR_FUNCTION) h = h->data.function.head;   /* s[a][b] */
    return h->type == EXPR_SYMBOL && h->data.symbol.name == sym;
}

/* The kept form of one step, or NULL. `prev_kept` (TF_SYMBOL) says that the
 * step before it in the same chain used a rule of the symbol. */
static Expr* trace_filter_step(Expr* step, const TraceFilter* tf, bool* uses) {
    Expr* body = trace_step_body(step);
    *uses = false;
    if (tf->kind == TF_SYMBOL) {
        *uses = trace_uses_symbol(body, tf->sym);
        return *uses ? expr_copy(step) : NULL;
    }
    MatchEnv* env = env_new();
    Expr* out = NULL;
    if (match(body, tf->pat, env)) {
        if (tf->kind == TF_PATTERN) {
            out = expr_copy(step);
        } else {
            Expr* inst = replace_bindings(tf->rhs, env);
            out = evaluate(inst);
            expr_free(inst);
        }
    }
    env_free(env);
    return out;
}

/* Filter the nested trace `node` (a nesting List, borrowed); returns the
 * filtered List, or NULL when nothing in it was kept. */
static Expr* trace_filter_tree(Expr* node, const TraceFilter* tf) {
    size_t n = node->data.function.arg_count;
    Expr** kept = n ? malloc(n * sizeof(Expr*)) : NULL;
    size_t k = 0;
    bool prev_uses = false;
    for (size_t i = 0; i < n; i++) {
        Expr* e = node->data.function.args[i];
        if (trace_is_nest(e)) {
            Expr* sub = trace_filter_tree(e, tf);
            if (sub) kept[k++] = sub;
            continue;     /* a sub-evaluation does not break the chain */
        }
        bool uses = false;
        Expr* f = trace_filter_step(e, tf, &uses);
        if (!f && tf->kind == TF_SYMBOL && prev_uses)
            f = expr_copy(e);   /* what a rule of the symbol rewrote it to */
        prev_uses = uses;
        if (f) kept[k++] = f;
    }
    if (k == 0) { free(kept); return NULL; }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), kept, k);
    free(kept);
    return out;
}

/* Trace requires HoldAll: the argument must reach the builtin unevaluated so
 * the collector can observe its rewrite sequence from the start. The steps
 * come back wrapped in HoldForm, which keeps the list inert on the
 * evaluator's re-pass (a bare 1 + 1 would reduce again) while printing
 * transparently: Trace[1 + 1] is {1 + 1, 2}.
 *
 * Ownership: the evaluator owns `res`; eval_collect_trace_marked BORROWS
 * args[0], so the returned List is independent of it. */
Expr* builtin_trace(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;  /* other arities stay unevaluated */

    Expr* raw = eval_collect_trace_marked(res->data.function.args[0]);
    if (!raw || argc == 1) return raw;

    /* Trace[expr, form]: the collector has already restored the outer trace
     * state, so match()'s internal evaluation (for /; and ?test) and a rule's
     * right-hand side are not themselves traced. */
    Expr* form = res->data.function.args[1];
    TraceFilter tf = { TF_PATTERN, form, NULL, NULL };
    if (form->type == EXPR_FUNCTION && form->data.function.arg_count == 2 &&
        form->data.function.head->type == EXPR_SYMBOL &&
        (form->data.function.head->data.symbol.name == SYM_Rule ||
         form->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
        tf.kind = TF_RULE;
        tf.pat = form->data.function.args[0];
        tf.rhs = form->data.function.args[1];
    } else if (form->type == EXPR_SYMBOL) {
        tf.kind = TF_SYMBOL;
        tf.sym = form->data.symbol.name;
    }
    Expr* out = trace_filter_tree(raw, &tf);
    expr_free(raw);
    return out ? out : expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
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
        "\tpattern form (e.g. Trace[expr, _Integer]), keeping the nesting.\n"
        "Trace[expr, s]\n\tIncludes the evaluations that use rules for the symbol s.\n"
        "Trace[expr, form -> rhs]\n\tGives rhs, with form's bindings, for each step\n"
        "\tthat matches form.");
}
