/*
 * dsolve_lower_reduce.c — DSolve`LowerDerivativeReduction.
 *
 * Reduces an ODE that is missing the dependent variable AND every derivative
 * below some order m >= 1 — i.e. only y^(m), y^(m+1), ..., y^(n) appear (n >= 3)
 * — by substituting p = y^(m).  The equation becomes an order-(n-m) ODE in p,
 * solved by recursing into the scalar engine; y is then recovered by integrating
 * p back m times (each integration introduces one constant, so the m constants
 * build the expected degree-(m-1) polynomial automatically).
 *
 * This is the general order >= 3 counterpart of ReductionOfOrder (which handles
 * the order-2, m=1 missing-y case) and AutonomousReduction (order-2 missing-x).
 * Gated to n >= 3 so it never competes with those at order 2.  Example:
 *   7 y' y''' - 11 y''^2 == 0   (only y', y'', y''' appear, m = 1)
 *     -> 7 p p'' - 11 p'^2 == 0  (autonomous, missing x) -> p = (4x+7C1)^(-7/4)
 *     -> y = Integrate[p, x] + C[3].
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>
#include <stdint.h>

/* True iff `body` satisfies the reduced order-`k` ODE red(p) == 0 identically
 * (p^(j) -> D[body, {x, j}]; structural zero test, PossibleZeroQ backstop).
 * Everything borrowed. */
static bool lr_satisfies(const Expr* red, const char* pfun, const char* xvar, int k,
                         const Expr* body) {
    Expr* t = expr_copy((Expr*)red);
    Expr* dj = expr_copy((Expr*)body);
    for (int j = 0; j <= k; j++) {
        t = ds_subst(t, ds_make_funcapp(pfun, j, xvar), expr_copy(dj));
        if (j < k) dj = ds_d(dj, expr_new_symbol(xvar));
    }
    expr_free(dj);
    t = eval_and_free(ds_call1(SYM_Together, t));
    bool ok = ds_is_zero(t);
    if (!ok) {
        Expr* pz = eval_and_free(ds_call1("PossibleZeroQ", expr_copy(t)));
        ok = (pz->type == EXPR_SYMBOL && pz->data.symbol.name == SYM_True);
        expr_free(pz);
    }
    expr_free(t);
    return ok;
}

/* Append the distinct Log[..] and square-root Power[.., +-1/2] subexpressions of
 * e to *out (borrowed pointers). */
static void lr_keys(const Expr* e, Expr*** out, size_t* n) {
    if (!e || e->type != EXPR_FUNCTION) return;
    bool key = head_is((Expr*)e, SYM_Log);
    if (!key && head_is((Expr*)e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* q = e->data.function.args[1];
        key = head_is((Expr*)q, SYM_Rational) && q->data.function.arg_count == 2
              && q->data.function.args[1]->type == EXPR_INTEGER
              && q->data.function.args[1]->data.integer == 2
              && q->data.function.args[0]->type == EXPR_INTEGER
              && q->data.function.args[0]->data.integer == 1;
    }
    if (key) {
        for (size_t i = 0; i < *n; i++) if (expr_eq((*out)[i], (Expr*)e)) return;
        *out = realloc(*out, (*n + 1) * sizeof(Expr*));
        (*out)[(*n)++] = (Expr*)e;
        if (head_is((Expr*)e, SYM_Log)) return;
    }
    lr_keys(e->data.function.head, out, n);
    for (size_t i = 0; i < e->data.function.arg_count; i++) lr_keys(e->data.function.args[i], out, n);
}

/* Collect[Expand[y], {Log terms, Sqrt terms}, Together]: the Cauchy sum leaves
 * x (..Sqrt[Q].. + ..Log[..]..) products undistributed.  An identity rewrite (no
 * branch choice).  Consumes y. */
static Expr* lr_collect(Expr* y) {
    if (ds_has_inactive_integrate(y)) return y;
    Expr** keys = NULL; size_t nk = 0;
    lr_keys(y, &keys, &nk);
    if (nk == 0) { free(keys); return y; }
    Expr** kc = malloc(nk * sizeof(Expr*));
    for (size_t i = 0; i < nk; i++) kc[i] = expr_copy(keys[i]);
    free(keys);
    Expr* klist = expr_new_function(expr_new_symbol(SYM_List), kc, nk);
    free(kc);
    Expr* ex = eval_and_free(ds_call1(SYM_Expand, expr_copy(y)));
    Expr* c = eval_and_free(expr_new_function(expr_new_symbol(intern_symbol("Collect")),
                  (Expr*[]){ ex, klist, expr_new_symbol(SYM_Together) }, 3));
    if (!c) return y;
    expr_free(y);
    return c;
}

/* y with y^(m) == p, via Cauchy's repeated-integration formula
 *   y = Sum_{k=0}^{m-1} (-1)^k x^(m-1-k) / (k! (m-1-k)!) Integrate[x^k p, x]
 * plus the polynomial C[base+1] x^(m-1) + ... + C[base+m].  ONE quadrature per k
 * instead of m nested ones: an m-fold antiderivative of p re-integrates the
 * previous antiderivative's transcendental parts (Log[L + Sqrt[Q]] after one pass
 * over Sqrt[Q]), which the integrator cannot close, while every x^k p stays in
 * p's own class (x^k Sqrt[Q] closes in QuadraticRadicals).  Each F_k' = x^k p
 * contributes x^(m-1) p Sum_k (-1)^k Binomial[m-1, k] == 0 to y', so y^(m) == p
 * by induction on m.  Each quadrature is bounded and muted; one that does not
 * close is held as Inactive[Integrate] when `allow_inert` (a general solution),
 * else the branch declines (NULL).  p borrowed. */
static Expr* lr_quadrature(const Expr* p, int m, int base, const char* xvar, bool allow_inert) {
    Expr* x = expr_new_symbol(xvar);
    Expr** terms = malloc((size_t)(2 * m) * sizeof(Expr*));
    size_t nt = 0;
    bool ok = true;
    int64_t fact[21]; fact[0] = 1;                 /* m <= 20 (gated by the caller) */
    for (int i = 1; i < 21; i++) fact[i] = fact[i - 1] * i;
    for (int k = 0; k < m && ok; k++) {
        Expr* f = eval_and_free(ds_call2(SYM_Times,
                      ds_call2(SYM_Power, expr_copy(x), expr_new_integer(k)), expr_copy((Expr*)p)));
        Expr* F = ds_integrate_or_inert(f, xvar, 5);     /* consumes f */
        if (!F || (!allow_inert && ds_has_inactive_integrate(F))) {
            if (F) expr_free(F);
            ok = false; break;
        }
        Expr* coef = eval_and_free(ds_call2(SYM_Times,
                         expr_new_integer((k % 2) ? -1 : 1),
                         ds_call2(SYM_Power, expr_new_integer(fact[k] * fact[m - 1 - k]),
                                  expr_new_integer(-1))));
        terms[nt++] = ds_call2(SYM_Times, coef,
                          ds_call2(SYM_Times, ds_call2(SYM_Power, expr_copy(x),
                                                       expr_new_integer(m - 1 - k)), F));
    }
    for (int i = 1; i <= m && ok; i++)
        terms[nt++] = ds_call2(SYM_Times, ds_const(base + i),
                          ds_call2(SYM_Power, expr_copy(x), expr_new_integer(m - i)));
    expr_free(x);
    if (!ok) {
        for (size_t i = 0; i < nt; i++) expr_free(terms[i]);
        free(terms);
        return NULL;
    }
    Expr* y = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, nt));
    free(terms);
    if (m >= 2) y = lr_collect(y);
    return y;
}

Expr** dsolve_lower_reduce_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    int n = P->max_order[0];
    if (n < 3) return NULL;   /* order-2 missing-var is owned by reduce_order/autonomous */
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];
    Expr* resid = P->eq_residuals[0];

    /* m = lowest derivative order actually present.  Detect each order by
     * substituting its funcapp for a marker (a bare-symbol scan would false-
     * positive on the y inside every Derivative[k][y][x]). */
    const char* mark = intern_symbol("DSolve`lrmark");
    int m = -1;
    for (int k = 0; k <= n && m < 0; k++) {
        Expr* t = ds_subst(expr_copy(resid), ds_make_funcapp(yname, k, xvar), expr_new_symbol(mark));
        bool present = ds_contains(t, mark);
        expr_free(t);
        if (present) m = k;
    }
    if (m < 1 || m > 20) return NULL;   /* y[x] itself is present (or nothing) -> not our case */

    /* substitute y^(m+j) -> p^(j) (j = n-m .. 0) to form the order-(n-m) p-ODE */
    const char* pfun = intern_symbol("DSolve`lrp");
    Expr* red = expr_copy(resid);
    for (int j = n - m; j >= 0; j--)
        red = ds_subst(red, ds_make_funcapp(yname, m + j, xvar), ds_make_funcapp(pfun, j, xvar));
    Expr* reduced = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ expr_copy(red), expr_new_integer(0) }, 2);

    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ reduced, ds_make_funcapp(pfun, 0, xvar), expr_new_symbol(xvar) }, 3);
    Expr* r = eval_and_free(call);
    /* EVERY branch of the reduced solution (a q-ODE solved by inversion has two:
     * taking only the first lost half the family). */
    size_t nb = 0;
    Expr** pb = dsolve_extract_applied_bodies(r, pfun, &nb);
    expr_free(r);
    if (!pb || nb == 0) { free(pb); expr_free(red); return NULL; }

    /* Integrate each p back m times; the reduced ODE used C[1..(n-m)], so the m
     * new integration constants are C[(n-m)+1 .. n]. */
    bool allow_inert = (P->ncond == 0) && !getenv("MATHILDA_NO_LR_INERT");
    Expr** out = malloc(nb * sizeof(Expr*));
    size_t nout = 0;
    for (size_t i = 0; i < nb; i++) {
        Expr* p = pb[i];
        /* Radical cleanup: the reduced solve squares out its inversion and leaves
         * Sqrt[s^2 R] with the square inside (the curvature ODE's p-solution is
         * -Sqrt[-u (1 + u)]/u with 1 + u a perfect square), hiding the
         * quadrature's class.  s Sqrt[R] is only +-Sqrt[s^2 R] pointwise, so the
         * rewrite is kept iff it still solves the reduced ODE. */
        Expr* pc = ds_sqrt_extract_squares(p);
        if (!expr_eq(pc, p) && lr_satisfies(red, pfun, xvar, n - m, pc)) {
            expr_free(p); p = pc;
        } else expr_free(pc);
        Expr* y = lr_quadrature(p, m, n - m, xvar, allow_inert);
        expr_free(p);
        if (y) out[nout++] = y;
    }
    free(pb);
    expr_free(red);
    if (nout == 0) { free(out); return NULL; }
    *nbranch = nout;
    return out;
}

static Expr* builtin_dsolve_lower_reduce(Expr* res) {
    return dsolve_method_builtin(res, dsolve_lower_reduce_try);
}

void dsolve_lower_reduce_init(void) {
    symtab_add_builtin("DSolve`LowerDerivativeReduction", builtin_dsolve_lower_reduce);
    symtab_get_def("DSolve`LowerDerivativeReduction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LowerDerivativeReduction",
        "DSolve`LowerDerivativeReduction[eqn, y, x] solves an order n >= 3 ODE in "
        "which y and every derivative below some order m >= 1 are absent: it sets "
        "p = y^(m), solves the resulting order-(n-m) ODE (every branch), and "
        "integrates p back m times by Cauchy's repeated-integration formula (one "
        "quadrature of x^k p per k; held inert in a general solution when it does "
        "not close).");
}
