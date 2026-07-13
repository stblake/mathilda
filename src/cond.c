#include "cond.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "match.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

Expr* builtin_if(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 4) return NULL;

    Expr* cond = res->data.function.args[0];

    if (cond->type == EXPR_SYMBOL) {
        if (cond->data.symbol.name == SYM_True) {
            return expr_copy(res->data.function.args[1]);
        } else if (cond->data.symbol.name == SYM_False) {
            if (argc >= 3) {
                return expr_copy(res->data.function.args[2]);
            } else {
                return expr_new_symbol(SYM_Null);
            }
        }
    }

    if (argc == 4) {
        return expr_copy(res->data.function.args[3]);
    }

    return NULL;
}

/*
 * builtin_which -- Which[t1, v1, t2, v2, ...]
 *
 * Which has the HoldAll attribute, so its arguments arrive unevaluated.
 * We evaluate the test expressions one pair at a time:
 *
 *   - if t_i evaluates to True   -> return v_i (held; the outer evaluator
 *                                   will continue evaluating it),
 *   - if t_i evaluates to False  -> drop this pair and continue,
 *   - otherwise                  -> return Which[t_i_eval, v_i, ...remaining...]
 *                                   with the inconclusive test in its
 *                                   evaluated form and the remaining
 *                                   arguments left unevaluated.
 *
 * If every test evaluates to False, Which returns Null. Which[] (no
 * arguments) likewise returns Null. An odd number of arguments is a
 * usage error: we leave the call unevaluated by returning NULL so the
 * user can see what they wrote.
 */
Expr* builtin_which(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) return expr_new_symbol(SYM_Null);
    if (argc % 2 != 0) return NULL;

    Expr** args = res->data.function.args;

    for (size_t i = 0; i < argc; i += 2) {
        Expr* test_eval = evaluate(args[i]);
        if (!test_eval) return NULL;

        bool is_true  = (test_eval->type == EXPR_SYMBOL &&
                         test_eval->data.symbol.name == SYM_True);
        bool is_false = (test_eval->type == EXPR_SYMBOL &&
                         test_eval->data.symbol.name == SYM_False);

        if (is_true) {
            expr_free(test_eval);
            return expr_copy(args[i + 1]);
        }

        if (is_false) {
            expr_free(test_eval);
            continue;
        }

        /* Inconclusive: return Which[t_i_eval, v_i, ...remaining...]. */
        size_t remaining = argc - i;
        Expr** new_args = malloc(sizeof(Expr*) * remaining);
        if (!new_args) {
            expr_free(test_eval);
            return NULL;
        }
        new_args[0] = test_eval;
        for (size_t j = 1; j < remaining; j++) {
            new_args[j] = expr_copy(args[i + j]);
        }
        Expr* head = expr_copy(res->data.function.head);
        Expr* out = expr_new_function(head, new_args, remaining);
        free(new_args);
        return out;
    }

    return expr_new_symbol(SYM_Null);
}

/*
 * builtin_switch -- Switch[expr, form_1, value_1, form_2, value_2, ...]
 *
 * Switch carries the HoldRest attribute, so:
 *   - arg 0 (the discriminant expr) is already evaluated by the outer
 *     evaluator,
 *   - args 1.. (the form/value pairs) arrive unevaluated.
 *
 * Per Mathematica semantics:
 *   - Each form_i is evaluated immediately before its match is tried
 *     (and only then; values for later branches are never touched).
 *   - The first form_i that pattern-matches expr causes the corresponding
 *     value_i to be returned (in its still-unevaluated form). The outer
 *     evaluation loop then evaluates value_i, which is exactly what makes
 *     "only the chosen value is evaluated" hold: every other value_i is
 *     simply never returned.
 *   - The pattern `_` (Blank[]) matches everything, so a trailing `_, def`
 *     pair acts as a default branch. This is just the natural consequence
 *     of normal pattern matching, not a special case.
 *   - If no form matches, the Switch call is returned unevaluated.
 *
 * Switch does NOT substitute pattern bindings (e.g. x_) from the matching
 * form into value_i; the form is purely a discriminator. This matches
 * Mathematica: Switch[{1,2}, {x_,y_}, x+y, _, 0] returns the literal
 * x+y, not 3. Bindings live and die inside the local MatchEnv.
 *
 * Argument-count rules: there must be at least one (form, value) pair
 * (so argc >= 3) and the number of arguments after expr must be even
 * (so argc must be odd). Any other arity is a usage error and we leave
 * the call unevaluated.
 */
Expr* builtin_switch(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc < 3 || (argc % 2) == 0) return NULL;

    Expr** args = res->data.function.args;
    Expr* expr = args[0]; /* already evaluated; HoldRest holds only args 1.. */

    for (size_t i = 1; i + 1 < argc; i += 2) {
        Expr* form_eval = evaluate(args[i]);
        if (!form_eval) return NULL;

        MatchEnv* env = env_new();
        bool matched = match(expr, form_eval, env);
        env_free(env);
        expr_free(form_eval);

        if (matched) {
            return expr_copy(args[i + 1]);
        }
    }

    /* No form matched: leave the Switch unevaluated. */
    return NULL;
}

Expr* builtin_trueq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;

    Expr* cond = res->data.function.args[0];
    bool is_true = (cond->type == EXPR_SYMBOL && cond->data.symbol.name == SYM_True);

    return expr_new_symbol(is_true ? "True" : "False");
}

/* ---------- Piecewise ---------------------------------------------------- *
 *
 *   Piecewise[{{v1, c1}, {v2, c2}, ...}]            -> Piecewise[..., 0]
 *   Piecewise[{{v1, c1}, {v2, c2}, ...}, default]
 *
 * HoldAll. Conditions are evaluated in turn. With surviving clauses:
 *   - {v, False} clauses are dropped;
 *   - at the first {v, True}, all later clauses (and the default) are
 *     dropped, and the True clause becomes the final case;
 *   - otherwise the (evaluated) condition is kept.
 * Consecutive surviving clauses with structurally equal values are merged
 * by combining their conditions with Or.
 *
 * Result:
 *   - 0 surviving clauses  -> default (auto-filled to 0 if omitted).
 *   - exactly 1 clause whose condition is True -> v (evaluated by caller).
 *   - otherwise            -> symbolic Piecewise.
 */

static bool piecewise_is_pair(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List
        && e->data.function.arg_count == 2;
}

static Expr* piecewise_make_pair(Expr* val, Expr* cond) {
    Expr** pa = malloc(2 * sizeof(Expr*));
    if (!pa) { expr_free(val); expr_free(cond); return NULL; }
    pa[0] = val;
    pa[1] = cond;
    Expr* pair = expr_new_function(expr_new_symbol(SYM_List), pa, 2);
    free(pa);
    return pair;
}

Expr* builtin_piecewise(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* clauses_arg = res->data.function.args[0];
    Expr* default_arg = (argc == 2) ? res->data.function.args[1] : NULL;

    /* First argument must be a list of {value, condition} pairs. */
    if (clauses_arg->type != EXPR_FUNCTION
        || clauses_arg->data.function.head->type != EXPR_SYMBOL
        || clauses_arg->data.function.head->data.symbol.name != SYM_List) {
        return NULL;
    }
    size_t n_in = clauses_arg->data.function.arg_count;
    Expr** in_pairs = clauses_arg->data.function.args;
    for (size_t i = 0; i < n_in; i++) {
        if (!piecewise_is_pair(in_pairs[i])) return NULL;
    }

    /* Walk clauses, collecting survivors. */
    Expr** out_vals = NULL;
    Expr** out_conds = NULL;
    size_t out_count = 0;
    size_t out_cap = 0;

    for (size_t i = 0; i < n_in; i++) {
        Expr* v = in_pairs[i]->data.function.args[0];
        Expr* c = in_pairs[i]->data.function.args[1];
        Expr* c_eval = evaluate(c);
        if (!c_eval) {
            for (size_t j = 0; j < out_count; j++) {
                expr_free(out_vals[j]);
                expr_free(out_conds[j]);
            }
            free(out_vals); free(out_conds);
            return NULL;
        }

        bool is_false = (c_eval->type == EXPR_SYMBOL && c_eval->data.symbol.name == SYM_False);
        bool is_true  = (c_eval->type == EXPR_SYMBOL && c_eval->data.symbol.name == SYM_True);

        if (is_false) {
            expr_free(c_eval);
            continue;
        }

        if (out_count == out_cap) {
            out_cap = out_cap ? out_cap * 2 : 8;
            out_vals  = realloc(out_vals,  out_cap * sizeof(Expr*));
            out_conds = realloc(out_conds, out_cap * sizeof(Expr*));
        }
        out_vals[out_count]  = expr_copy(v);
        out_conds[out_count] = c_eval;
        out_count++;

        if (is_true) break;
    }

    /* Merge consecutive clauses with structurally equal values:
     *   {v, c1}, {v, c2}, ...  ->  {v, Or[c1, c2, ...]}
     * Evaluate the merged Or so True/False simplifications kick in. */
    if (out_count > 1) {
        size_t w = 0;
        for (size_t r = 0; r < out_count; ) {
            size_t k = r + 1;
            while (k < out_count && expr_eq(out_vals[r], out_vals[k])) {
                k++;
            }
            if (k - r == 1) {
                out_vals[w]  = out_vals[r];
                out_conds[w] = out_conds[r];
                w++;
            } else {
                size_t span = k - r;
                Expr** or_args = malloc(span * sizeof(Expr*));
                for (size_t j = 0; j < span; j++) or_args[j] = out_conds[r + j];
                Expr* or_expr = expr_new_function(expr_new_symbol(SYM_Or), or_args, span);
                free(or_args);
                /* eval_and_free consumes or_expr (which now owns the
                 * per-clause conditions) and returns the reduced form. */
                Expr* or_eval = eval_and_free(or_expr);

                for (size_t j = 1; j < span; j++) expr_free(out_vals[r + j]);
                out_vals[w]  = out_vals[r];
                out_conds[w] = or_eval;
                w++;
            }
            r = k;
        }
        out_count = w;
    }

    /* No surviving clauses: result is the default. */
    if (out_count == 0) {
        free(out_vals); free(out_conds);
        if (default_arg) return expr_copy(default_arg);
        return expr_new_integer(0);
    }

    /* A single clause whose surviving condition is True dictates the value.
     * This covers both a literal first-True hit and a merged Or that
     * simplified to True. The returned value is held; the outer evaluator
     * will reduce it. */
    bool last_is_true = (out_conds[out_count - 1]->type == EXPR_SYMBOL
                         && out_conds[out_count - 1]->data.symbol.name == SYM_True);
    if (out_count == 1 && last_is_true) {
        Expr* v = out_vals[0];
        expr_free(out_conds[0]);
        free(out_vals); free(out_conds);
        return v;
    }

    /* Build symbolic Piecewise[{{v, c}, ...}, default?]. The default is
     * dropped iff the final clause is unconditionally True (it would be
     * unreachable). */
    Expr** pair_nodes = malloc(out_count * sizeof(Expr*));
    for (size_t i = 0; i < out_count; i++) {
        pair_nodes[i] = piecewise_make_pair(out_vals[i], out_conds[i]);
    }
    Expr* clauses_list = expr_new_function(expr_new_symbol(SYM_List), pair_nodes, out_count);
    free(pair_nodes);
    free(out_vals); free(out_conds);

    Expr** pw_args;
    size_t pw_argc;
    if (last_is_true) {
        pw_args = malloc(sizeof(Expr*));
        pw_args[0] = clauses_list;
        pw_argc = 1;
    } else {
        pw_args = malloc(2 * sizeof(Expr*));
        pw_args[0] = clauses_list;
        pw_args[1] = default_arg ? expr_copy(default_arg) : expr_new_integer(0);
        pw_argc = 2;
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Piecewise), pw_args, pw_argc);
    free(pw_args);

    /* If we produced the exact same expression as the input, signal "no
     * change" so the outer evaluator can fix-point efficiently. */
    if (expr_eq(result, res)) {
        expr_free(result);
        return NULL;
    }
    return result;
}

void cond_init(void) {
    symtab_add_builtin("If", builtin_if);
    symtab_get_def("If")->attributes |= ATTR_HOLDREST | ATTR_PROTECTED;

    symtab_add_builtin("Which", builtin_which);
    symtab_get_def("Which")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;

    symtab_add_builtin("Switch", builtin_switch);
    symtab_get_def("Switch")->attributes |= ATTR_HOLDREST | ATTR_PROTECTED;

    symtab_add_builtin("TrueQ", builtin_trueq);
    symtab_get_def("TrueQ")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("Piecewise", builtin_piecewise);
    symtab_get_def("Piecewise")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}