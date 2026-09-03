/*
 * reduce_eq.c
 *
 * Complete univariate equation engine for `Reduce` (REDUCE_PLAN.md, Phase 1).
 * See reduce_eq.h for the algorithm.  Everything numeric is done by reusing the
 * existing Expr-level builtins (Exponent, Coefficient, Expand, Solve) rather
 * than reimplementing polynomial arithmetic.
 */
#include "reduce_eq.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>

/* ------------------------------------------------------------------ *
 *  Expr-level polynomial helpers (via evaluate)                       *
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* Degree of `poly` in `var` via Exponent[poly, var]; -1 if not a plain int. */
static int poly_degree(const Expr* poly, const Expr* var) {
    Expr* e = eval_and_free(expr_new_function(expr_new_symbol(SYM_Exponent),
        (Expr*[]){ expr_copy((Expr*)poly), expr_copy((Expr*)var) }, 2));
    int d = (e->type == EXPR_INTEGER) ? (int)e->data.integer : -1;
    expr_free(e);
    return d;
}

/* Coefficient[poly, var, deg]. */
static Expr* poly_coeff(const Expr* poly, const Expr* var, int deg) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Coefficient),
        (Expr*[]){ expr_copy((Expr*)poly), expr_copy((Expr*)var),
                   expr_new_integer(deg) }, 3));
}

/* Expand[poly - lc*var^deg] -- `poly` with its leading term removed. */
static Expr* poly_drop_leading(const Expr* poly, const Expr* var,
                               const Expr* lc, int deg) {
    Expr* pw = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy((Expr*)var), expr_new_integer(deg) }, 2);
    Expr* term = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)lc), pw }, 3);
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy((Expr*)poly), term }, 2);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
        (Expr*[]){ sum }, 1));
}

static bool is_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* A nonzero numeric literal: the leading coefficient is definitely nonzero, so
 * no `lc != 0 / lc == 0` case split is needed. */
static bool is_nonzero_numeric(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer != 0;
        case EXPR_REAL:    return e->data.real != 0.0;
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint) != 0;
#ifdef USE_MPFR
        case EXPR_MPFR:    return mpfr_sgn(e->data.mpfr) != 0;
#endif
        case EXPR_FUNCTION:
            /* A Rational is always nonzero (a zero would be Integer 0). */
            return is_head(e, SYM_Rational);
        default: return false;
    }
}

/* Collect the distinct parameter symbols of `e` -- leaves that are NOT one of the
 * reduce vars (a condition is a relation among the free parameters).  Owned. */
static void collect_params(const Expr* e, Expr** vars, int nv, Expr*** out, int* n) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < nv; i++) if (expr_eq(e, vars[i])) return;    /* a reduce var */
        for (int i = 0; i < *n; i++) if (expr_eq((*out)[i], e)) return;  /* dup */
        *out = realloc(*out, (size_t)(*n + 1) * sizeof(Expr*));
        (*out)[(*n)++] = expr_copy((Expr*)e);
        return;
    }
    if (e->type == EXPR_FUNCTION)
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_params(e->data.function.args[i], vars, nv, out, n);
}

/* When `poly` is linear in a single parameter `a` with a nonzero numeric
 * coefficient, return the solved display  a REL (-const/coeff)  (owned); else
 * NULL.  Mirrors reduce_sys.c's sys_norm_condition so a leading-coefficient
 * condition prints minimally: `-1 + a != 0` becomes `a != 1`, `-1 + 2 a == 0`
 * becomes `a == 1/2`.  A nonlinear or multi-parameter condition is left verbatim. */
static Expr* cond_solved_display(const char* rel_head, const Expr* poly,
                                 Expr** vars, int nv) {
    Expr** params = NULL; int nparams = 0;
    collect_params(poly, vars, nv, &params, &nparams);
    Expr* result = NULL;
    if (nparams == 1 && poly_degree(poly, params[0]) == 1) {
        Expr* a = params[0];
        Expr* coeff = poly_coeff(poly, a, 1);
        if (is_nonzero_numeric(coeff)) {
            Expr* c0 = poly_coeff(poly, a, 0);              /* constant term */
            Expr* val = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_new_integer(-1), c0,
                           expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_copy(coeff), expr_new_integer(-1) }, 2) }, 3));
            result = expr_new_function(expr_new_symbol(rel_head),
                (Expr*[]){ expr_copy(a), val }, 2);
        }
        expr_free(coeff);
    }
    for (int i = 0; i < nparams; i++) expr_free(params[i]);
    free(params);
    return result;
}

/* Build a bare relational atom  poly REL 0  (rel is Equal or Unequal).  When the
 * condition is single-parameter-linear, attach a minimal solved-form display. */
static RForm* condition_form(const char* rel_head, const Expr* poly,
                             Expr** vars, int nv, bool* ok) {
    Expr* rel = expr_new_function(expr_new_symbol(rel_head),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    RAtom a = reduce_atom_canonicalize(rel, vars, nv, ok);
    expr_free(rel);
    if (*ok) {
        Expr* disp = cond_solved_display(rel_head, poly, vars, nv);
        if (disp) { if (a.display) expr_free(a.display); a.display = disp; }
    }
    return rform_single_atom(a);
}

/* ------------------------------------------------------------------ *
 *  Generic solutions via Solve                                        *
 * ------------------------------------------------------------------ */

/* Solve[poly == 0, var] over Complexes, reformatted as an OR of solved atoms
 * (x == r_k).  Sets *ok=false on any unexpected shape (declined, conditional,
 * multi-rule, wrong variable). */
static RForm* solve_generic(const Expr* poly, const Expr* var,
                            Expr** vars, int nv, bool* ok,
                            const ReduceOpts* opts) {
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    Expr* base[2] = { eqn, expr_copy((Expr*)var) };
    Expr* call = reduce_opts_build_solve(base, 2, opts);
    Expr* sols = eval_and_free(call);

    if (!is_head(sols, SYM_List)) { *ok = false; expr_free(sols); return rform_false(); }

    RForm* out = rform_false();
    size_t nsol = sols->data.function.arg_count;
    for (size_t i = 0; i < nsol; i++) {
        Expr* row = sols->data.function.args[i];
        if (!is_head(row, SYM_List)) { *ok = false; break; }
        if (row->data.function.arg_count == 0) {   /* {{}} -> all values */
            rform_free(out);
            out = rform_true();
            break;
        }
        if (row->data.function.arg_count != 1) { *ok = false; break; }
        Expr* rule = row->data.function.args[0];
        if (!(is_head(rule, SYM_Rule) || is_head(rule, SYM_RuleDelayed))
            || rule->data.function.arg_count != 2) { *ok = false; break; }
        Expr* lhs = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        if (!expr_eq(lhs, var)) { *ok = false; break; }
        if (is_head(val, SYM_ConditionalExpression)) { *ok = false; break; }
        out = rform_or(out, rform_single_atom(ratom_solved(var, val, vars, nv)));
    }
    expr_free(sols);
    if (!*ok) { rform_free(out); return rform_false(); }
    return out;
}

/* ------------------------------------------------------------------ *
 *  The leading-coefficient-vanishing recursion                        *
 * ------------------------------------------------------------------ */

static RForm* solve_case(const Expr* poly, const Expr* var,
                         Expr** vars, int nv, bool* ok,
                         const ReduceOpts* opts) {
    if (!*ok) return rform_false();

    /* An identically-zero polynomial is `0 == 0`, true for every value of the
     * remaining variable.  This is the base case of a leading-coefficient split
     * whose lower-order part vanished: `Reduce[a x == 0, x]` drops the `a x`
     * term in its `a == 0` branch, leaving 0 -> True, so the branch is `a == 0`
     * (not a decline).  Guard it here because Exponent[0, x] is -Infinity, which
     * poly_degree cannot classify. */
    if (is_zero(poly)) return rform_true();

    int deg = poly_degree(poly, var);
    if (deg < 0) { *ok = false; return rform_false(); }

    if (deg == 0) {
        /* No x left: the equation is the parameter condition poly == 0. */
        if (is_zero(poly)) return rform_true();
        return condition_form(SYM_Equal, poly, vars, nv, ok);
    }

    Expr* lc = poly_coeff(poly, var, deg);
    RForm* gen = solve_generic(poly, var, vars, nv, ok, opts);
    if (!*ok) { expr_free(lc); rform_free(gen); return rform_false(); }

    if (is_nonzero_numeric(lc)) {          /* terminal: no case split needed */
        expr_free(lc);
        return gen;
    }

    /* Branch A: lc != 0  &&  (generic solutions of degree deg). */
    RForm* A = rform_and(condition_form(SYM_Unequal, lc, vars, nv, ok), gen);

    /* Branch B: lc == 0  &&  solve(poly with leading term dropped). */
    Expr* reduced = poly_drop_leading(poly, var, lc, deg);
    RForm* B = rform_and(condition_form(SYM_Equal, lc, vars, nv, ok),
                         solve_case(reduced, var, vars, nv, ok, opts));
    expr_free(reduced);
    expr_free(lc);

    return rform_or(A, B);
}

RForm* reduce_eq_univariate(const Expr* poly, const Expr* var,
                            Expr** vars, int nv, bool* ok,
                            const ReduceOpts* opts) {
    *ok = true;
    return solve_case(poly, var, vars, nv, ok, opts);
}

/* ------------------------------------------------------------------ *
 *  Transcendental single equation -> logical formula                  *
 * ------------------------------------------------------------------ */

Expr* reduce_eq_transcendental(const Expr* poly, const Expr* var,
                               const Expr* dom, const ReduceOpts* opts) {
    /* Re-enter Solve[poly == 0, var (, dom)]; base args are adopted by the
     * Solve call.  An explicit `dom` (Reals/Complexes) is forwarded so the
     * inverse-function solver produces domain-appropriate solutions and
     * conditions; a NULL dom leaves Solve at its Complexes default. */
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    Expr* base[3];
    int nb = 0;
    base[nb++] = eqn;
    base[nb++] = expr_copy((Expr*)var);
    if (dom) base[nb++] = expr_copy((Expr*)dom);
    Expr* sols = eval_and_free(reduce_opts_build_solve(base, nb, opts));

    /* Solve declined (unevaluated / non-list): let the caller fall back. */
    if (!is_head(sols, SYM_List)) { expr_free(sols); return NULL; }

    size_t nsol = sols->data.function.arg_count;
    Expr** terms = (Expr**)malloc((nsol ? nsol : 1) * sizeof(Expr*));
    size_t nt = 0;
    bool ok = true;
    Expr* all_values = NULL;   /* set to True when a {{}} row means "any var" */

    for (size_t i = 0; i < nsol; i++) {
        Expr* row = sols->data.function.args[i];
        if (!is_head(row, SYM_List)) { ok = false; break; }
        if (row->data.function.arg_count == 0) {           /* {{}} -> True */
            all_values = expr_new_symbol(SYM_True);
            break;
        }
        if (row->data.function.arg_count != 1) { ok = false; break; }
        Expr* rule = row->data.function.args[0];
        if (!(is_head(rule, SYM_Rule) || is_head(rule, SYM_RuleDelayed))
            || rule->data.function.arg_count != 2) { ok = false; break; }
        Expr* lhs = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        if (!expr_eq(lhs, var)) { ok = false; break; }

        Expr* term;
        if (is_head(val, SYM_ConditionalExpression)
            && val->data.function.arg_count == 2) {
            /* ConditionalExpression[value, cond] -> cond && (var == value). */
            Expr* value = val->data.function.args[0];
            Expr* cond  = val->data.function.args[1];
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ expr_copy((Expr*)var), expr_copy(value) }, 2);
            term = expr_new_function(expr_new_symbol(SYM_And),
                (Expr*[]){ expr_copy(cond), eq }, 2);
        } else {
            term = expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ expr_copy((Expr*)var), expr_copy(val) }, 2);
        }
        terms[nt++] = term;
    }
    expr_free(sols);

    if (!ok) {
        for (size_t i = 0; i < nt; i++) expr_free(terms[i]);
        free(terms);
        if (all_values) expr_free(all_values);
        return NULL;
    }
    if (all_values) {
        for (size_t i = 0; i < nt; i++) expr_free(terms[i]);
        free(terms);
        return all_values;
    }
    if (nt == 0) { free(terms); return expr_new_symbol(SYM_False); }  /* {} */

    Expr* out;
    if (nt == 1) { out = terms[0]; free(terms); }
    else { out = expr_new_function(expr_new_symbol(SYM_Or), terms, nt); free(terms); }
    /* Let the evaluator flatten / tidy the And/Or combination. */
    return eval_and_free(out);
}
