/*
 * dsolve_euler.c — DSolve`EulerCauchy.
 *
 * Solves the equidimensional (Cauchy-Euler) linear ODE, including the shifted
 * form centred at any constant b,
 *     a_n (x-b)^n y^(n) + ... + a_1 (x-b) y' + a_0 y == g(x)   (a_k constant)
 * by the trial y = (x-b)^r, which gives the indicial polynomial
 *     P(r) = Σ a_k r(r-1)...(r-k+1).
 * A real root r of multiplicity m contributes (x-b)^r (Log(x-b))^j (j=0..m-1);
 * a complex pair p +- q i contributes (x-b)^p Cos[q Log(x-b)](Log(x-b))^j and
 * the Sin form.
 *
 * The centre b is recovered from the leading coefficient: c_n(x) = a_n (x-b)^n
 * gives n c_n / c_n' = (x-b), hence b = x - n c_n / c_n' (0 for the classic
 * equation, the pole for a shifted one).  The equation is Cauchy-Euler about b
 * iff every a_k = c_k / (x-b)^k is free of x (and a_n != 0).
 *
 * The inhomogeneous solution splits by the nature of the roots, because the two
 * particular-solution techniques have complementary blind spots:
 *
 *   - All roots real: variation of parameters in x, with the (x-b)^r basis.
 *     The x-domain integrals keep special-function forcing intact — e.g.
 *     ∫ e^x / x dx -> ExpIntegralEi[x] for  x^2 y'' - 2x y' + 2y == x^2 e^x —
 *     and with no complex roots there are no trig-of-Log products to stall the
 *     final Simplify.
 *
 *   - Any complex root: reduce to a CONSTANT-COEFFICIENT ODE by (x-b) = e^t
 *     (t = Log(x-b)).  Then (x-b)^k y^(k) = θ(θ-1)...(θ-k+1) Y with θ = d/dt, so
 *     the operator is exactly P(θ); writing P(r) = Σ b_j r^j the transformed
 *     equation is  Σ b_j Y^(j)(t) == g(b + e^t), which the constant-coefficient
 *     engine solves with a hang-free particular solution (undetermined
 *     coefficients for exp-poly-trig forcing).  Mapping t -> Log(x-b) recovers
 *     y — (x-b)^r, (Log(x-b))^j and Cos/Sin fall out of ordinary evaluation
 *     (E^(c Log(x-b)) auto-reduces to (x-b)^c).  Doing the particular solution
 *     in t rather than x-domain variation of parameters avoids the products of
 *     trig-of-Log terms with irrational frequencies that make the x-domain
 *     Simplify blow up — e.g. the complex-root inhomogeneous x^2 y'' + y == x^2.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* Log[w]^j (j >= 1); w borrowed. */
static Expr* logpow_w(const Expr* w, int j) {
    Expr* lg = ds_call1("Log", expr_copy((Expr*)w));
    if (j == 1) return lg;
    return expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ lg, expr_new_integer(j) }, 2);
}
/* w^r (Log w)^j ; r, w borrowed. */
static Expr* euler_real(const Expr* w, const Expr* r, int j) {
    Expr* wr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_copy((Expr*)w), expr_copy((Expr*)r) }, 2));
    if (j == 0) return wr;
    return eval_and_free(ds_call2(SYM_Times, wr, logpow_w(w, j)));
}
/* w^p Cos/Sin[q Log w] (Log w)^j ; p,q,w borrowed. */
static Expr* euler_trig(const Expr* w, const Expr* p, const Expr* q, int j, const char* which) {
    Expr* wp = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_copy((Expr*)w), expr_copy((Expr*)p) }, 2));
    Expr* qlogw = ds_call2(SYM_Times, expr_copy((Expr*)q), ds_call1("Log", expr_copy((Expr*)w)));
    Expr* trig = ds_call1(which, qlogw);
    Expr* prod = eval_and_free(ds_call2(SYM_Times, wp, trig));
    if (j == 0) return prod;
    return eval_and_free(ds_call2(SYM_Times, prod, logpow_w(w, j)));
}

/* Homogeneous fundamental set (n functions) in the base w = (x-b) from the
 * indicial roots R (R.total == n).  Returns bc; basis[] owned. */
static size_t euler_basis(const Expr* w, DSolveRoots* R, int n, Expr** basis) {
    bool* used = calloc(R->ndist, sizeof(bool));
    size_t bc = 0;
    for (size_t i = 0; i < R->ndist && bc <= (size_t)n; i++) {
        if (used[i]) continue;
        if (R->isreal[i]) {
            for (int j = 0; j < R->mult[i]; j++) basis[bc++] = euler_real(w, R->roots[i], j);
            used[i] = true;
        } else {
            Expr* conj = eval_and_free(ds_call1("Conjugate", expr_copy(R->roots[i])));
            long cc = -1;
            for (size_t k = 0; k < R->ndist && cc < 0; k++) {
                if (k == i || used[k]) continue;
                Expr* diff = eval_and_free(ds_call2(SYM_Subtract, expr_copy(R->roots[k]), expr_copy(conj)));
                if (ds_is_zero(diff)) cc = (long)k;
                expr_free(diff);
            }
            expr_free(conj);
            if (cc >= 0) {
                used[i] = used[(size_t)cc] = true;
                Expr* re = eval_and_free(ds_call1("Re", expr_copy(R->roots[i])));
                for (int j = 0; j < R->mult[i]; j++) {
                    basis[bc++] = euler_trig(w, re, R->im[i], j, "Cos");
                    basis[bc++] = euler_trig(w, re, R->im[i], j, "Sin");
                }
                expr_free(re);
            } else {
                for (int j = 0; j < R->mult[i]; j++) basis[bc++] = euler_real(w, R->roots[i], j);
                used[i] = true;
            }
        }
    }
    free(used);
    return bc;
}

/* Real-root path: homogeneous basis in w plus x-domain variation of parameters.
 * Returns the general solution (owned) or NULL to decline. */
static Expr* euler_real_solve(const Expr* w, DSolveRoots* R, int n,
                              const Expr* g, const Expr* leadcoef, const char* xvar) {
    Expr** basis = malloc((size_t)n * sizeof(Expr*));
    size_t bc = euler_basis(w, R, n, basis);
    Expr* general = NULL;
    if (bc == (size_t)n) {
        Expr** hterms = malloc((size_t)n * sizeof(Expr*));
        for (int k = 0; k < n; k++)
            hterms[k] = ds_call2(SYM_Times, ds_const(k + 1), expr_copy(basis[k]));
        Expr* Hgen = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                         hterms, (size_t)n));
        free(hterms);
        if (ds_is_zero(g)) {
            general = Hgen;
        } else {
            Expr* yp = dsolve_variation_of_parameters(basis, (size_t)n, g, leadcoef, xvar);
            if (yp) general = eval_and_free(ds_call2(SYM_Plus, Hgen, yp));
            else expr_free(Hgen);
        }
    }
    for (size_t i = 0; i < bc; i++) expr_free(basis[i]);
    free(basis);
    return general;
}

/* Complex-root / fallback path: reduce to a constant-coefficient ODE in
 * t = Log(x-b), solve it, and map back.  `indicial` and `w` borrowed.  Returns
 * the general solution (owned) or NULL to decline. */
static Expr* euler_transform_solve(const Expr* indicial, const char* rv, int n,
                                   const Expr* g, const Expr* bshift, const Expr* w,
                                   const char* xvar) {
    const char* tvar = intern_symbol("DSolve`et");
    const char* Yfun = intern_symbol("DSolve`ecY");
    const char* coef = intern_symbol("Coefficient");
    Expr** lhs = malloc((size_t)(n + 1) * sizeof(Expr*));
    for (int j = 0; j <= n; j++) {
        Expr* bj = eval_and_free(expr_new_function(expr_new_symbol(coef),
                       (Expr*[]){ expr_copy((Expr*)indicial), expr_new_symbol(rv),
                                  expr_new_integer(j) }, 3));
        lhs[j] = ds_call2(SYM_Times, bj, ds_make_funcapp(Yfun, j, tvar));
    }
    Expr* LHS = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                    lhs, (size_t)(n + 1)));
    free(lhs);
    /* RHS = g(x -> b + e^t) */
    Expr* xoft = eval_and_free(ds_call2(SYM_Plus, expr_copy((Expr*)bshift),
                     ds_call1("Exp", expr_new_symbol(tvar))));
    Expr* RHS = ds_subst(expr_copy((Expr*)g), expr_new_symbol(xvar), xoft);
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
                    (Expr*[]){ LHS, RHS }, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, ds_make_funcapp(Yfun, 0, tvar),
                                expr_new_symbol(tvar) }, 3);
    Expr* r = eval_and_free(call);

    size_t nb = 0;
    Expr** bodies = dsolve_extract_applied_bodies(r, Yfun, &nb);
    expr_free(r);
    Expr* general = NULL;
    if (bodies) {
        if (nb >= 1 && bodies[0]) {
            /* undo (x-b) = e^t:  t -> Log(x-b)  (E^(c Log w) auto-reduces to w^c). */
            general = ds_subst(bodies[0], expr_new_symbol(tvar), ds_call1("Log", expr_copy((Expr*)w)));
            for (size_t i = 1; i < nb; i++) expr_free(bodies[i]);
        } else {
            for (size_t i = 0; i < nb; i++) expr_free(bodies[i]);
        }
        free(bodies);
    }
    return general;
}

Expr** dsolve_euler_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] < 1) return NULL;
    const char* xvar = P->ind_names[0];

    Expr** c; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return NULL;
    dsolve_linear_normalize(c, &g, n, P->ind_names[0]);

    /* Centre b of the equidimensional operator: c_n(x) = a_n (x-b)^n gives
     * n c_n / c_n' = (x-b), so b = x - n c_n / c_n'.  b = 0 for the classic
     * equation; the pole for a shifted one.  b must be free of x. */
    Expr* dcn = ds_d(expr_copy(c[n]), expr_new_symbol(xvar));                    /* c_n'(x) */
    Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ dcn, expr_new_integer(-1) }, 2));                 /* 1 / c_n' */
    Expr* ncn = eval_and_free(ds_call2(SYM_Times, expr_new_integer(n), expr_copy(c[n]))); /* n c_n */
    Expr* frac = eval_and_free(ds_call2(SYM_Times, ncn, inv));                   /* n c_n / c_n' */
    Expr* bshift = eval_and_free(ds_call2(SYM_Subtract, expr_new_symbol(xvar), frac)); /* b */
    Expr* w = eval_and_free(ds_call2(SYM_Subtract, expr_new_symbol(xvar), expr_copy(bshift))); /* x - b */

    /* a_k = c_k / w^k must be free of x, and a_n != 0. */
    Expr** a = malloc((size_t)(n + 1) * sizeof(Expr*));
    bool ok = ds_free_of(bshift, xvar);
    for (int k = 0; k <= n; k++) {
        a[k] = eval_and_free(ds_call2(SYM_Times, expr_copy(c[k]),
                   eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(w), expr_new_integer(-k) }, 2))));
        if (ok && !ds_free_of(a[k], xvar)) ok = false;
    }
    if (ok && ds_is_zero(a[n])) ok = false;

    Expr** out = NULL;
    if (ok) {
        /* indicial polynomial P(r) = Σ a_k (r)_k  (falling factorial (r)_k). */
        const char* rv = intern_symbol("DSolve`r");
        Expr* ff = expr_new_integer(1);               /* (r)_0 = 1 */
        Expr** terms = malloc((size_t)(n + 1) * sizeof(Expr*));
        for (int k = 0; k <= n; k++) {
            terms[k] = eval_and_free(ds_call2(SYM_Times, expr_copy(a[k]), expr_copy(ff)));
            ff = eval_and_free(ds_call2(SYM_Times, ff,
                     ds_call2(SYM_Subtract, expr_new_symbol(rv), expr_new_integer(k))));
        }
        Expr* indicial = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                             terms, (size_t)(n + 1)));
        free(terms);
        expr_free(ff);

        DSolveRoots R;
        bool analyzed = dsolve_analyze_roots(indicial, rv, n, &R);
        bool all_real = analyzed && R.total == n;
        if (all_real)
            for (size_t i = 0; i < R.ndist; i++) if (!R.isreal[i]) { all_real = false; break; }

        Expr* general = all_real
            ? euler_real_solve(w, &R, n, g, c[n], xvar)
            : euler_transform_solve(indicial, rv, n, g, bshift, w, xvar);

        if (analyzed) dsolve_roots_free(&R);
        expr_free(indicial);

        if (general) { out = malloc(sizeof(Expr*)); out[0] = general; *nbranch = 1; }
    }

    for (int k = 0; k <= n; k++) { expr_free(a[k]); expr_free(c[k]); }
    free(a); free(c);
    expr_free(g); expr_free(bshift); expr_free(w);
    return out;
}

static Expr* builtin_dsolve_euler(Expr* res) {
    return dsolve_method_builtin(res, dsolve_euler_try);
}

void dsolve_euler_init(void) {
    symtab_add_builtin("DSolve`EulerCauchy", builtin_dsolve_euler);
    symtab_get_def("DSolve`EulerCauchy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`EulerCauchy",
        "DSolve`EulerCauchy[eqn, y, x] solves the equidimensional equation "
        "a_n (x-b)^n y^(n) + ... + a_0 y == g(x) (any constant centre b) via the "
        "indicial polynomial (trial y = (x-b)^r): complex roots give Cos/Sin of "
        "Log[x-b], repeated roots give Log[x-b] powers.  A forcing g(x) is handled "
        "by variation of parameters (real roots; keeps ExpIntegralEi and other "
        "special functions) or, for complex roots, the constant-coefficient "
        "reduction (x-b) = e^t.");
}
