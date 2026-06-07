/*
 * sum_hypergeometric.c -- Sum`Hypergeometric: infinite hypergeometric series.
 *
 * For a summand t(k) whose term ratio t(k+1)/t(k) is a rational function of k,
 * the infinite sum from a concrete lower bound is a generalized hypergeometric
 * function:
 *
 *     Sum_{k=imin}^{Infinity} t(k) = t(imin) * pFq({a}, {b}, z),
 *
 * where the upper/lower parameters and the argument z are read off by factoring
 * the (reindexed) term ratio into monic linear factors over Q:
 *
 *     t(imin+m+1)/t(imin+m) = z * prod_i (m + a_i) / ( prod_j (m + b_j) (m+1) ).
 *
 * The (m+1) is the canonical factorial factor.  Numerator roots give the upper
 * parameters, the remaining denominator roots the lower parameters, and the
 * ratio of leading coefficients is z.  The result is handed to HypergeometricPFQ,
 * which reduces the common cases to elementary closed forms (e.g. 0F0 -> E^z,
 * the geometric 1F0(1) -> 1/(1-z)) and otherwise returns the pFq itself.
 *
 *   Sum`Hypergeometric[f, k, imin, Infinity]  -> t(imin) pFq({a},{b},z)
 *
 * Only the infinite case with a concrete integer lower bound is handled; every
 * other shape (and a non-rational term ratio, an irrational/complex root, or a
 * divergent p > q+1 series) falls through unevaluated so the cascade can try the
 * next stage or leave Sum[...] held.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "poly.h"
#include "arithmetic.h"   /* is_infinity_sym */
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#define SUM_HG_MAX_DEGREE 256   /* cap on the linear-factor deflation loop */

/* e is a polynomial in var? */
static bool hg_poly_in(Expr* e, Expr* var) {
    Expr* vars[1] = { var };
    return is_polynomial(e, vars, 1);
}

/* degree of e in var (expanding first). */
static int hg_pdeg(Expr* e, Expr* var) {
    Expr* ea[1] = { expr_copy(e) };
    Expr* ex = sum_eval("Expand", ea, 1);
    int d = get_degree_poly(ex, var);
    expr_free(ex);
    return d;
}

/* Pull the first root r out of a Solve result List[List[Rule[var, r]], ...].
 * Returns an owned copy of r, or NULL if the shape is anything else (e.g. an
 * unsolved Solve[...] or an empty solution set). */
static Expr* hg_first_root(Expr* sol) {
    if (!sol || sol->type != EXPR_FUNCTION) return NULL;
    if (sol->data.function.head->type != EXPR_SYMBOL
        || sol->data.function.head->data.symbol != SYM_List
        || sol->data.function.arg_count < 1) return NULL;
    Expr* s0 = sol->data.function.args[0];
    if (s0->type != EXPR_FUNCTION
        || s0->data.function.head->type != EXPR_SYMBOL
        || s0->data.function.head->data.symbol != SYM_List
        || s0->data.function.arg_count < 1) return NULL;
    Expr* rule = s0->data.function.args[0];
    if (rule->type != EXPR_FUNCTION
        || rule->data.function.head->type != EXPR_SYMBOL
        || rule->data.function.head->data.symbol != SYM_Rule
        || rule->data.function.arg_count != 2) return NULL;
    return expr_copy(rule->data.function.args[1]);
}

/* Factor polynomial P (in var) as lc * prod_i (var + a_i) with each a_i free
 * of var.  On success returns true and fills *params (owned array, *np entries,
 * each owned) and *lc (owned, the remaining constant).  Returns false if P does
 * not split into rational/parametric linear factors (e.g. an irreducible
 * quadratic). */
static bool hg_linear_params(Expr* P, Expr* var,
                             Expr*** params, size_t* np, Expr** lc) {
    size_t cap = 4, n = 0;
    Expr** ps = malloc(sizeof(Expr*) * cap);
    Expr* work = expr_copy(P);

    for (int guard = 0; guard < SUM_HG_MAX_DEGREE; guard++) {
        int d = hg_pdeg(work, var);
        if (d <= 0) { *params = ps; *np = n; *lc = work; return true; }

        Expr* eq = expr_new_function(expr_new_symbol("Equal"),
                       (Expr*[]){ expr_copy(work), sum_int(0) }, 2);
        Expr* sol = sum_eval("Solve", (Expr*[]){ eq, expr_copy(var) }, 2);
        Expr* r = hg_first_root(sol);
        expr_free(sol);
        if (!r || !sum_free_of(r, var)) {
            if (r) expr_free(r);
            goto fail;
        }

        /* param = -r */
        Expr* param = sum_eval("Times", (Expr*[]){ sum_int(-1), expr_copy(r) }, 2);
        if (n == cap) { cap *= 2; ps = realloc(ps, sizeof(Expr*) * cap); }
        ps[n++] = param;

        /* deflate: work = PolynomialQuotient[work, var - r, var] */
        Expr* lin = expr_new_function(expr_new_symbol("Plus"),
                        (Expr*[]){ expr_copy(var),
                                   expr_new_function(expr_new_symbol("Times"),
                                       (Expr*[]){ sum_int(-1), r }, 2) }, 2);
        Expr* q = sum_eval("PolynomialQuotient",
                      (Expr*[]){ expr_copy(work), lin, expr_copy(var) }, 3);
        expr_free(work);
        work = q;
        if (hg_pdeg(work, var) >= d) goto fail;   /* degree must drop */
    }

fail:
    for (size_t i = 0; i < n; i++) expr_free(ps[i]);
    free(ps);
    expr_free(work);
    return false;
}

/* List[elems...] from an owned array (adopts the elements; copies pointers). */
static Expr* hg_make_list(Expr** elems, size_t n) {
    Expr* l = expr_new_function(expr_new_symbol("List"), elems, n);
    return l;
}

Expr* builtin_sum_hypergeometric(Expr* res);

Expr* builtin_sum_hypergeometric(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;

    /* Only infinite sums with a concrete integer lower bound. */
    if (!is_infinity_sym(imax)) return NULL;
    if (imin->type != EXPR_INTEGER) return NULL;

    /* Term ratio t(k+1)/t(k), reduced to a rational function num/den in var. */
    Expr* nv = expr_new_function(expr_new_symbol("Plus"),
                   (Expr*[]){ expr_copy(var), sum_int(1) }, 2);
    Expr* tshift = sum_subst(f, var, nv);
    expr_free(nv);
    Expr* ratio_raw = expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){ tshift,
                                     expr_new_function(expr_new_symbol("Power"),
                                         (Expr*[]){ expr_copy(f), sum_int(-1) }, 2) }, 2);
    Expr* simp = sum_eval("Simplify", (Expr*[]){ ratio_raw }, 1);
    Expr* ratio = sum_eval("Together", (Expr*[]){ simp }, 1);
    Expr* num = sum_eval("Numerator", (Expr*[]){ expr_copy(ratio) }, 1);
    Expr* den = sum_eval("Denominator", (Expr*[]){ ratio }, 1);

    if (!hg_poly_in(num, var) || !hg_poly_in(den, var)) {
        expr_free(num); expr_free(den);
        return NULL;
    }

    /* Reindex to m >= 0: substitute var -> var + imin (no-op when imin == 0). */
    Expr* shift = expr_new_function(expr_new_symbol("Plus"),
                      (Expr*[]){ expr_copy(var), expr_copy(imin) }, 2);
    Expr* num_m = sum_subst(num, var, shift);
    Expr* den_m = sum_subst(den, var, shift);
    expr_free(shift);
    expr_free(num); expr_free(den);

    Expr **A = NULL, **B = NULL, *lcN = NULL, *lcD = NULL;
    size_t p = 0, q = 0;
    if (!hg_linear_params(num_m, var, &A, &p, &lcN)) { expr_free(num_m); expr_free(den_m); return NULL; }
    if (!hg_linear_params(den_m, var, &B, &q, &lcD)) {
        for (size_t i = 0; i < p; i++) expr_free(A[i]);
        free(A); expr_free(lcN); expr_free(num_m); expr_free(den_m);
        return NULL;
    }
    expr_free(num_m); expr_free(den_m);

    /* Upper parameters start as the numerator roots. */
    size_t up_cap = p + 1, up_n = 0;
    Expr** up = malloc(sizeof(Expr*) * (up_cap ? up_cap : 1));
    for (size_t i = 0; i < p; i++) up[up_n++] = A[i];
    free(A);

    /* Lower parameters: the denominator roots minus one canonical (m+1)
     * factorial factor (param == 1).  If none is present, introduce the
     * factorial by adding an upper parameter 1 (the (m+1) in the denominator). */
    size_t lo_cap = q ? q : 1, lo_n = 0;
    Expr** lo = malloc(sizeof(Expr*) * lo_cap);
    Expr* one = sum_int(1);
    bool factorial_removed = false;
    for (size_t j = 0; j < q; j++) {
        if (!factorial_removed && expr_eq(B[j], one)) {
            factorial_removed = true;
            expr_free(B[j]);
        } else {
            lo[lo_n++] = B[j];
        }
    }
    free(B);
    expr_free(one);
    if (!factorial_removed) up[up_n++] = sum_int(1);

    /* z = lcN / lcD. */
    Expr* z = sum_eval("Times",
                  (Expr*[]){ lcN,
                             expr_new_function(expr_new_symbol("Power"),
                                 (Expr*[]){ lcD, sum_int(-1) }, 2) }, 2);

    /* Convergence gate.  p <= q is entire (sums for all z).  p == q+1 converges
     * only for |z| < 1 -- outside that, the elementary closed form is the
     * analytic continuation, NOT the value of the (divergent) series, so we must
     * NOT return it (e.g. Sum[2^k, {k,0,Infinity}] must stay Infinity/unevaluated,
     * not 1/(1-2) = -1).  p > q+1 always diverges.  When z is symbolic the gate
     * passes and the closed form is returned as the conventional CAS answer. */
    bool diverges = false;
    if (up_n > lo_n + 1) {
        diverges = true;
    } else if (up_n == lo_n + 1) {
        Expr* absz = sum_eval("Abs", (Expr*[]){ expr_copy(z) }, 1);
        Expr* nz = sum_eval("N", (Expr*[]){ absz }, 1);
        if (nz->type == EXPR_REAL && nz->data.real >= 1.0 - 1e-12) diverges = true;
        expr_free(nz);
    }
    if (diverges) {
        for (size_t i = 0; i < up_n; i++) expr_free(up[i]);
        for (size_t j = 0; j < lo_n; j++) expr_free(lo[j]);
        free(up); free(lo); expr_free(z);
        return NULL;
    }

    /* prefactor t(imin). */
    Expr* prefactor = sum_subst(f, var, imin);

    Expr* upList = hg_make_list(up, up_n); free(up);
    Expr* loList = hg_make_list(lo, lo_n); free(lo);
    Expr* pfq = expr_new_function(expr_new_symbol("HypergeometricPFQ"),
                    (Expr*[]){ upList, loList, z }, 3);
    Expr* result = sum_eval("Times", (Expr*[]){ prefactor, pfq }, 2);
    return result;
}

void sum_hypergeometric_init(void) {
    symtab_add_builtin("Sum`Hypergeometric", builtin_sum_hypergeometric);
}
