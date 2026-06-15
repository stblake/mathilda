/*
 * iter.c -- Iteration built-ins for Mathilda.
 *
 * This module implements Mathematica-style imperative iteration primitives
 * that sit outside the rule/pattern-based evaluation model:
 *
 *   Do[expr, spec...]            -- counted / iterator-driven loop
 *   For[start, test, incr, body] -- C-style three-part loop
 *   While[test, body]            -- pre-test loop
 *
 * All three have the HoldAll attribute so that their argument expressions
 * (iterator specs, test expressions, bodies) are re-evaluated on every
 * iteration rather than once at call time.
 *
 * Control-flow semantics (shared by all loops in this file):
 *
 *   Return[val]    -- exits the loop; the value val is returned as the
 *                     result of the loop expression.
 *   Break[]        -- exits the loop; the loop yields Null.
 *   Continue[]     -- aborts the current iteration and proceeds with the
 *                     next one (skipping any remaining body for Do/While,
 *                     or proceeding to the increment step for For).
 *   Throw[...],
 *   Abort[],
 *   Quit[]         -- propagate unchanged to the enclosing handler.
 *
 * These head names are detected by inspecting the expression returned
 * from evaluate() on the body; CompoundExpression already bubbles them
 * up unchanged so that `Body1; Body2; Break[]; Body3` short-circuits on
 * Break[] as expected.
 *
 * Memory-ownership contract (see SPEC.md section 4):
 *   Each builtin_* function receives ownership of its input `res` from
 *   the evaluator but MUST NOT free it -- the evaluator is the actual
 *   owner of the argument tree. Every Expr* allocated internally here
 *   (including copies made when restoring iterator variable bindings)
 *   is freed on every exit path, including error returns and control-
 *   flow short-circuits.
 */

#include "iter.h"
#include "symtab.h"
#include "eval.h"
#include "core.h"
#include "arithmetic.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

/*
 * ============================================================================
 *  Shared iterator-spec helpers (see iter.h for the contract)
 * ============================================================================
 */

void iter_spec_free(IterSpec* s) {
    if (!s) return;
    if (s->var)  expr_free(s->var);
    if (s->imin) expr_free(s->imin);
    if (s->imax) expr_free(s->imax);
    if (s->di)   expr_free(s->di);
    if (s->list) expr_free(s->list);
    memset(s, 0, sizeof(*s));
}

static bool is_list_expr(const Expr* e) {
    return e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List;
}

bool iter_spec_parse(Expr* spec, IterSpec* out) {
    memset(out, 0, sizeof(*out));
    out->kind = ITER_KIND_COUNT;

    if (is_list_expr(spec)) {
        size_t len = spec->data.function.arg_count;
        if (len == 1) {
            /* {n} -- repeat count only */
            out->imax = evaluate(spec->data.function.args[0]);
            return true;
        }
        if (len < 2) return false;             /* {} -- malformed */

        /* {i, ...} -- first arg must be the iterator symbol. */
        Expr* v = spec->data.function.args[0];
        if (v->type != EXPR_SYMBOL) return false;
        out->var = expr_copy(v);

        if (len == 2) {
            /* {i, bound}: list iteration if bound is a List, else {i,1,imax}. */
            Expr* bound = evaluate(spec->data.function.args[1]);
            if (is_list_expr(bound)) {
                out->kind = ITER_KIND_LIST;
                out->list = bound;
            } else {
                out->kind = ITER_KIND_RANGE;
                out->imin = expr_new_integer(1);
                out->imax = bound;
                out->di   = expr_new_integer(1);
            }
            return true;
        }
        if (len == 3) {
            out->kind = ITER_KIND_RANGE;
            out->imin = evaluate(spec->data.function.args[1]);
            out->imax = evaluate(spec->data.function.args[2]);
            out->di   = expr_new_integer(1);
            return true;
        }
        if (len == 4) {
            out->kind = ITER_KIND_RANGE;
            out->imin = evaluate(spec->data.function.args[1]);
            out->imax = evaluate(spec->data.function.args[2]);
            out->di   = evaluate(spec->data.function.args[3]);
            return true;
        }
        /* len >= 5 -- malformed; free the var copy already taken. */
        iter_spec_free(out);
        return false;
    }

    /* Bare count: f[..., n] */
    out->imax = evaluate(spec);
    return true;
}

bool iter_spec_resolve_numeric(const IterSpec* s, bool allow_inf,
                               double* min_val, double* max_val,
                               double* di_val, bool* is_real, bool* is_inf) {
    *min_val = 0; *max_val = 0; *di_val = 0;
    *is_real = false; *is_inf = false;
    int64_t n, d;

    if (s->kind == ITER_KIND_COUNT) {
        if (allow_inf && s->imax->type == EXPR_SYMBOL
                      && s->imax->data.symbol == SYM_Infinity) {
            *is_inf = true;
            return true;
        }
        if (s->imax->type != EXPR_INTEGER) return false;
        *max_val = (double)s->imax->data.integer;
        return true;
    }
    if (s->kind != ITER_KIND_RANGE) return false;

    if (allow_inf && s->imax->type == EXPR_SYMBOL
                  && s->imax->data.symbol == SYM_Infinity) {
        *is_inf = true;
    }
    if (s->imin->type == EXPR_REAL || s->imax->type == EXPR_REAL
                                   || s->di->type   == EXPR_REAL) {
        *is_real = true;
    }

    if (s->imin->type == EXPR_INTEGER)      *min_val = (double)s->imin->data.integer;
    else if (s->imin->type == EXPR_REAL)    *min_val = s->imin->data.real;
    else if (is_rational(s->imin, &n, &d))  *min_val = (double)n / d;
    else return false;

    if (!*is_inf) {
        if (s->imax->type == EXPR_INTEGER)      *max_val = (double)s->imax->data.integer;
        else if (s->imax->type == EXPR_REAL)    *max_val = s->imax->data.real;
        else if (is_rational(s->imax, &n, &d))  *max_val = (double)n / d;
        else return false;
    }

    if (s->di->type == EXPR_INTEGER)      *di_val = (double)s->di->data.integer;
    else if (s->di->type == EXPR_REAL)    *di_val = s->di->data.real;
    else if (is_rational(s->di, &n, &d))  *di_val = (double)n / d;
    else return false;

    if (*di_val == 0) return false;        /* zero step never terminates */
    return true;
}

Rule* iter_spec_shadow(Expr* var) {
    if (!var || var->type != EXPR_SYMBOL) return NULL;
    SymbolDef* def = symtab_get_def(var->data.symbol);
    Rule* old = def->own_values;
    def->own_values = NULL;
    return old;
}

void iter_spec_restore(Expr* var, Rule* saved_own) {
    if (!var || var->type != EXPR_SYMBOL) return;
    SymbolDef* def = symtab_get_def(var->data.symbol);
    Rule* r = def->own_values;
    while (r) {
        Rule* next = r->next;
        expr_free(r->pattern);
        expr_free(r->replacement);
        free(r);
        r = next;
    }
    def->own_values = saved_own;
}

/*
 * is_flow_control_head
 *
 * Returns the head-name string if `e` is a function of the form
 * Return[...], Break[], Continue[], Throw[...], Abort[], or Quit[] --
 * i.e. an expression that must abort the current iteration. Returns
 * NULL otherwise. The returned pointer aliases memory owned by `e`
 * and must not be freed by the caller.
 */
static const char* is_flow_control_head(Expr* e) {
    if (e->type != EXPR_FUNCTION) return NULL;
    if (e->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = e->data.function.head->data.symbol;
    if (h == SYM_Return ||
        h == SYM_Break ||
        h == SYM_Continue ||
        h == SYM_Throw ||
        h == SYM_Abort ||
        h == SYM_Quit) {
        return h;
    }
    return NULL;
}

/*
 * IterFlowAction / iter_flow_classify
 *
 * Single decision point for "what should this Do/For/While step do
 * with the value just produced by evaluate(body)?". Wraps:
 *
 *   - eval_classify_return  -- handles Return[], Return[v], Return[v,h]
 *                              including the 2-arg targeted form.
 *   - is_flow_control_head  -- handles Break, Continue, Throw, Abort, Quit.
 *
 * The Return classifier runs first so that Return[v,h] for h !=
 * boundary_head propagates correctly past the loop.
 *
 *   ITER_FLOW_NONE        -- the value is a normal expression. Caller
 *                              should free eval_expr and continue
 *                              iterating.
 *   ITER_FLOW_BREAK       -- Break[]. Caller frees eval_expr and breaks
 *                              out of its loop; the loop yields Null.
 *   ITER_FLOW_CONTINUE    -- Continue[]. Caller frees eval_expr and
 *                              jumps to the next iteration (Do
 *                              arithmetic-progression must still
 *                              advance its counter before looping).
 *   ITER_FLOW_RETURN_VAL  -- Return targeted at this loop. *value_out
 *                              is set to the caller-owned payload
 *                              (Null for Return[], expr_copy of args[0]
 *                              otherwise). Caller stores it as
 *                              returned_val, frees eval_expr, breaks.
 *   ITER_FLOW_PROPAGATE   -- Throw/Abort/Quit, OR Return[v,h] whose
 *                              target h does not match boundary_head.
 *                              Caller hands eval_expr up unchanged
 *                              (assigns to returned_val and breaks)
 *                              WITHOUT freeing it.
 */
typedef enum {
    ITER_FLOW_NONE,
    ITER_FLOW_BREAK,
    ITER_FLOW_CONTINUE,
    ITER_FLOW_RETURN_VAL,
    ITER_FLOW_PROPAGATE
} IterFlowAction;

static IterFlowAction iter_flow_classify(Expr* eval_expr,
                                         const char* boundary_head,
                                         Expr** value_out) {
    if (value_out) *value_out = NULL;

    Expr* rv = NULL;
    EvalReturnAction ra = eval_classify_return(eval_expr, boundary_head, &rv);
    if (ra == EVAL_RETURN_CONSUME) {
        if (value_out) *value_out = rv;
        else expr_free(rv);
        return ITER_FLOW_RETURN_VAL;
    }
    if (ra == EVAL_RETURN_PROPAGATE) return ITER_FLOW_PROPAGATE;

    const char* h = is_flow_control_head(eval_expr);
    if (!h) return ITER_FLOW_NONE;
    if (h == SYM_Break)    return ITER_FLOW_BREAK;
    if (h == SYM_Continue) return ITER_FLOW_CONTINUE;
    /* Throw / Abort / Quit propagate unchanged. Return cannot reach
     * here: eval_classify_return above already dispatched it. */
    return ITER_FLOW_PROPAGATE;
}

/*
 * ============================================================================
 *  Do
 * ============================================================================
 *
 * Forms supported:
 *   Do[expr, n]                         -- evaluate expr n times
 *   Do[expr, {imax}]                    -- same, written with a 1-tuple
 *   Do[expr, {i, imax}]                 -- i = 1, 2, ..., imax
 *   Do[expr, {i, imin, imax}]           -- i = imin, imin+1, ..., imax
 *   Do[expr, {i, imin, imax, di}]       -- i = imin, imin+di, ..., <= imax
 *   Do[expr, {i, list}]                 -- i ranges over list elements
 *   Do[expr, spec1, spec2, ...]         -- nested loops (innermost last)
 *
 * Multi-spec Do is rewritten as nested two-spec Do invocations and re-
 * evaluated; that keeps the arithmetic below a single-iterator problem.
 */
Expr* builtin_do(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;

    /*
     * Multi-spec Do: rewrite Do[expr, s1, s2, ..., sk] as
     *     Do[Do[expr, sk], s1, s2, ..., s_{k-1}]
     * and let the evaluator recurse. We construct new expression trees
     * and free them after evaluation; the returned value is owned by us.
     */
    if (res->data.function.arg_count > 2) {
        Expr** inner_args = malloc(sizeof(Expr*) * 2);
        inner_args[0] = expr_copy(res->data.function.args[0]);
        inner_args[1] = expr_copy(res->data.function.args[res->data.function.arg_count - 1]);
        Expr* inner_do = expr_new_function(expr_new_symbol(SYM_Do), inner_args, 2);
        free(inner_args);

        size_t outer_argc = res->data.function.arg_count - 1;
        Expr** outer_args = malloc(sizeof(Expr*) * outer_argc);
        outer_args[0] = inner_do;
        for (size_t i = 1; i < outer_argc; i++) {
            outer_args[i] = expr_copy(res->data.function.args[i]);
        }
        Expr* outer_do = expr_new_function(expr_new_symbol(SYM_Do), outer_args, outer_argc);
        free(outer_args);

        Expr* eval_outer = evaluate(outer_do);
        expr_free(outer_do);
        return eval_outer;
    }

    Expr* expr = res->data.function.args[0];   /* body, held */
    Expr* spec = res->data.function.args[1];   /* iterator spec, held */

    /* ---- Parse the iterator spec (shared helper) ---- */
    IterSpec s;
    if (!iter_spec_parse(spec, &s)) return NULL;

    int is_n_times   = (s.kind == ITER_KIND_COUNT);
    int is_list_iter = (s.kind == ITER_KIND_LIST);
    double min_val = 0, max_val = 0, di_val = 0;
    bool is_real = false;
    bool is_inf  = false;

    /* ---- Validate and convert iterator bounds (Infinity allowed) ---- */
    if (!is_list_iter) {
        if (!iter_spec_resolve_numeric(&s, /*allow_inf=*/true,
                                       &min_val, &max_val, &di_val,
                                       &is_real, &is_inf)) {
            iter_spec_free(&s);
            return NULL;
        }
    }

    /* Convenience aliases into the owned IterSpec (freed via iter_spec_free). */
    Expr* var_sym = s.var;
    Expr* imin_e  = s.imin;
    Expr* imax_e  = s.imax;
    Expr* di_e    = s.di;
    Expr* list_e  = s.list;

    /*
     * Shadow any existing OwnValue for the iterator: we temporarily clear
     * it so that evaluation of the body sees successive iterator values,
     * and restore the original binding after the loop exits. This mimics
     * Mathematica's localisation of Do's iterator.
     */
    Rule* old_own = iter_spec_shadow(var_sym);

    Expr* returned_val = NULL;

    if (is_n_times) {
        /* Pure repeat-N-times form. */
        int64_t n = is_inf ? 0 : imax_e->data.integer;
        for (int64_t i = 0; is_inf || i < n; i++) {
            Expr* eval_expr = evaluate(expr);
            Expr* rv = NULL;
            IterFlowAction f = iter_flow_classify(eval_expr, SYM_Do, &rv);
            if (f == ITER_FLOW_RETURN_VAL) { returned_val = rv;        expr_free(eval_expr); break; }
            if (f == ITER_FLOW_PROPAGATE)  { returned_val = eval_expr;                       break; }
            if (f == ITER_FLOW_BREAK)      {                            expr_free(eval_expr); break; }
            if (f == ITER_FLOW_CONTINUE)   {                            expr_free(eval_expr); continue; }
            expr_free(eval_expr);
        }
    } else if (is_list_iter) {
        /* {i, {a, b, c, ...}} -- iterate over explicit list elements. */
        for (size_t i = 0; i < list_e->data.function.arg_count; i++) {
            symtab_add_own_value(var_sym->data.symbol, var_sym, list_e->data.function.args[i]);
            Expr* eval_expr = evaluate(expr);
            Expr* rv = NULL;
            IterFlowAction f = iter_flow_classify(eval_expr, SYM_Do, &rv);
            if (f == ITER_FLOW_RETURN_VAL) { returned_val = rv;        expr_free(eval_expr); break; }
            if (f == ITER_FLOW_PROPAGATE)  { returned_val = eval_expr;                       break; }
            if (f == ITER_FLOW_BREAK)      {                            expr_free(eval_expr); break; }
            if (f == ITER_FLOW_CONTINUE)   {                            expr_free(eval_expr); continue; }
            expr_free(eval_expr);
        }
    } else {
        /*
         * {i, imin, imax, di} -- arithmetic progression.
         *
         * We maintain two synchronised iterators:
         *   val    -- a double used for the terminating comparison, and
         *   curr_e -- the exact Expr* value bound to the iterator var.
         * For exact (integer/rational) arithmetic, curr_e is updated via
         * Plus[curr_e, di_e]; val is then refreshed from curr_e. For real
         * iteration we just add di_val to val directly.
         */
        double val = min_val;
        Expr* curr_e = expr_copy(imin_e);
        while (is_inf || (di_val > 0 && val <= max_val + 1e-14) || (di_val < 0 && val >= max_val - 1e-14)) {
            Expr* i_val = is_real ? expr_new_real(val) : expr_copy(curr_e);
            symtab_add_own_value(var_sym->data.symbol, var_sym, i_val);

            Expr* eval_expr = evaluate(expr);
            expr_free(i_val);

            Expr* rv = NULL;
            IterFlowAction f = iter_flow_classify(eval_expr, SYM_Do, &rv);
            if (f == ITER_FLOW_RETURN_VAL) { returned_val = rv;        expr_free(eval_expr); break; }
            if (f == ITER_FLOW_PROPAGATE)  { returned_val = eval_expr;                       break; }
            if (f == ITER_FLOW_BREAK)      {                            expr_free(eval_expr); break; }
            if (f == ITER_FLOW_CONTINUE)   {
                /* Still advance the iterator before re-testing the loop. */
                expr_free(eval_expr);
                Expr* next_args[2] = { expr_copy(curr_e), expr_copy(di_e) };
                Expr* next_expr = expr_new_function(expr_new_symbol(SYM_Plus), next_args, 2);
                Expr* next_e = evaluate(next_expr);
                expr_free(next_expr);
                expr_free(curr_e);
                curr_e = next_e;
                if (!is_real) {
                    int64_t n, d;
                    if (curr_e->type == EXPR_INTEGER) val = (double)curr_e->data.integer;
                    else if (curr_e->type == EXPR_REAL) val = curr_e->data.real;
                    else if (is_rational(curr_e, &n, &d)) val = (double)n / d;
                } else {
                    val += di_val;
                }
                continue;
            }
            expr_free(eval_expr);

            /* Normal end-of-iteration: step curr_e and val forward. */
            Expr* next_args[2] = { expr_copy(curr_e), expr_copy(di_e) };
            Expr* next_expr = expr_new_function(expr_new_symbol(SYM_Plus), next_args, 2);
            Expr* next_e = evaluate(next_expr);
            expr_free(next_expr);
            expr_free(curr_e);
            curr_e = next_e;

            if (!is_real) {
                int64_t n, d;
                if (curr_e->type == EXPR_INTEGER) val = (double)curr_e->data.integer;
                else if (curr_e->type == EXPR_REAL) val = curr_e->data.real;
                else if (is_rational(curr_e, &n, &d)) val = (double)n / d;
            } else {
                val += di_val;
            }
        }
        if (curr_e) expr_free(curr_e);
    }

    /*
     * Restore the iterator variable's original OwnValue binding (frees the
     * per-iteration OwnValue chain we created), then free the spec.
     */
    iter_spec_restore(var_sym, old_own);
    iter_spec_free(&s);

    if (returned_val) return returned_val;
    return expr_new_symbol(SYM_Null);
}

/*
 * ============================================================================
 *  For
 * ============================================================================
 *
 *   For[start, test, incr]        -- body defaulted to Null
 *   For[start, test, incr, body]
 *
 * Semantics: evaluate `start` once, then loop {evaluate `test`; if not True
 * exit; evaluate `body`; evaluate `incr`}. Break/Continue/Return are honoured
 * inside both `body` and `incr` (Continue in body proceeds to incr; Continue
 * in incr proceeds to the next test). Returns the Return value if one is
 * issued, else Null.
 */
Expr* builtin_for(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 3 || argc > 4) return NULL;

    Expr* start = res->data.function.args[0];
    Expr* test  = res->data.function.args[1];
    Expr* incr  = res->data.function.args[2];
    Expr* body  = (argc == 4) ? res->data.function.args[3] : NULL;

    /* Initialisation: evaluate `start` for its side-effects; discard result. */
    Expr* eval_start = evaluate(start);
    expr_free(eval_start);

    Expr* returned_val = NULL;

    while (true) {
        /* --- Test --- */
        Expr* eval_test = evaluate(test);
        bool condition_met = (eval_test->type == EXPR_SYMBOL && eval_test->data.symbol == SYM_True);
        expr_free(eval_test);
        if (!condition_met) break;

        /* --- Body --- */
        if (body) {
            Expr* eval_body = evaluate(body);
            Expr* rv = NULL;
            IterFlowAction f = iter_flow_classify(eval_body, SYM_For, &rv);
            if (f == ITER_FLOW_RETURN_VAL) { returned_val = rv;        expr_free(eval_body); break; }
            if (f == ITER_FLOW_PROPAGATE)  { returned_val = eval_body;                       break; }
            if (f == ITER_FLOW_BREAK)      {                            expr_free(eval_body); break; }
            /* ITER_FLOW_CONTINUE: fall through to the increment step.
             * ITER_FLOW_NONE:     normal value, discard. */
            expr_free(eval_body);
        }

        /* --- Increment --- */
        Expr* eval_incr = evaluate(incr);
        {
            Expr* rv = NULL;
            IterFlowAction f = iter_flow_classify(eval_incr, SYM_For, &rv);
            if (f == ITER_FLOW_RETURN_VAL) { returned_val = rv;        expr_free(eval_incr); break; }
            if (f == ITER_FLOW_PROPAGATE)  { returned_val = eval_incr;                       break; }
            if (f == ITER_FLOW_BREAK)      {                            expr_free(eval_incr); break; }
            /* ITER_FLOW_CONTINUE here: a Continue[] inside the increment
             * proceeds to the next test (same as a normal increment). */
            expr_free(eval_incr);
        }
    }

    if (returned_val) return returned_val;
    return expr_new_symbol(SYM_Null);
}

/*
 * ============================================================================
 *  While
 * ============================================================================
 *
 *   While[test]          -- repeatedly evaluate `test`; body is implicitly Null.
 *                           Useful when `test` has side-effects.
 *   While[test, body]    -- repeatedly evaluate `test`, then `body`, stopping
 *                           as soon as `test` fails to give True.
 *
 * Behaviour:
 *   - HoldAll: both `test` and `body` are re-evaluated every iteration.
 *   - If `test` does not evaluate to the symbol True, the loop exits and
 *     `body` is not evaluated this iteration.
 *   - Break[]    inside `body` exits the loop, returning Null.
 *   - Continue[] inside `body` skips the remainder of `body` and jumps
 *                back to re-evaluating `test`.
 *   - Return[v]  inside `body` exits and causes While to yield v.
 *   - Throw/Abort/Quit propagate.
 *   - If `test` itself returns a flow-control head, it propagates in the
 *     same way (so e.g. Throw inside test still escapes cleanly).
 *   - The loop always returns Null in the absence of an explicit Return.
 */
Expr* builtin_while(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* test = res->data.function.args[0];
    Expr* body = (argc == 2) ? res->data.function.args[1] : NULL;

    Expr* returned_val = NULL;

    while (true) {
        /* --- Evaluate the test --- */
        Expr* eval_test = evaluate(test);

        /* Honour flow-control from within the test expression itself. */
        {
            Expr* rv = NULL;
            IterFlowAction f = iter_flow_classify(eval_test, SYM_While, &rv);
            if (f == ITER_FLOW_RETURN_VAL) { returned_val = rv;        expr_free(eval_test); break; }
            if (f == ITER_FLOW_PROPAGATE)  { returned_val = eval_test;                       break; }
            if (f == ITER_FLOW_BREAK)      {                            expr_free(eval_test); break; }
            if (f == ITER_FLOW_CONTINUE)   { /* restart the loop */     expr_free(eval_test); continue; }
        }

        /* Exit the loop unless the test evaluated to the symbol True. */
        bool condition_met = (eval_test->type == EXPR_SYMBOL && eval_test->data.symbol == SYM_True);
        expr_free(eval_test);
        if (!condition_met) break;

        /* --- Evaluate the body (if any) --- */
        if (body) {
            Expr* eval_body = evaluate(body);
            Expr* rv = NULL;
            IterFlowAction f = iter_flow_classify(eval_body, SYM_While, &rv);
            if (f == ITER_FLOW_RETURN_VAL) { returned_val = rv;        expr_free(eval_body); break; }
            if (f == ITER_FLOW_PROPAGATE)  { returned_val = eval_body;                       break; }
            if (f == ITER_FLOW_BREAK)      {                            expr_free(eval_body); break; }
            if (f == ITER_FLOW_CONTINUE)   {
                /* Skip rest of this iteration; re-evaluate test next. */
                expr_free(eval_body);
                continue;
            }
            /* ITER_FLOW_NONE: discard body value and loop. */
            expr_free(eval_body);
        }
    }

    if (returned_val) return returned_val;
    return expr_new_symbol(SYM_Null);
}

/*
 * Register the iteration builtins and set their attributes. All three are
 * HoldAll (so that test, body, and iterator-spec expressions are not pre-
 * evaluated by the standard evaluator) and Protected to prevent user
 * redefinition.
 */
void iter_init(void) {
    symtab_add_builtin("Do", builtin_do);
    symtab_get_def("Do")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;

    symtab_add_builtin("For", builtin_for);
    symtab_get_def("For")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;

    symtab_add_builtin("While", builtin_while);
    symtab_get_def("While")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}
