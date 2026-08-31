/*
 * dsolve_euler.c — DSolve`EulerCauchy.
 *
 * Solves the equidimensional (Cauchy-Euler) linear ODE
 *     a_n x^n y^(n) + ... + a_1 x y' + a_0 y == g(x)   (a_k constant)
 * by the trial y = x^r, which gives the indicial polynomial
 *     P(r) = Σ a_k r(r-1)...(r-k+1) = 0.
 * A real root r of multiplicity m contributes x^r (Log x)^j (j=0..m-1); a
 * complex pair a +- b i contributes x^a Cos[b Log x](Log x)^j and the Sin form.
 * Inhomogeneous g(x) is handled by variation of parameters with leading
 * coefficient a_n x^n.
 *
 * The coefficient of y^(k) is read as c_k(x) (dsolve_linear_coeffs); the equation
 * is Cauchy-Euler iff a_k = c_k / x^k is free of x.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* Log[x]^j (j >= 1) */
static Expr* logpow(const char* xvar, int j) {
    Expr* lg = ds_call1("Log", expr_new_symbol(xvar));
    if (j == 1) return lg;
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ lg, expr_new_integer(j) }, 2);
}
/* x^r (Log x)^j ; r borrowed */
static Expr* euler_real(const char* xvar, const Expr* r, int j) {
    Expr* xr = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_new_symbol(xvar), expr_copy((Expr*)r) }, 2);
    xr = eval_and_free(xr);
    if (j == 0) return xr;
    return eval_and_free(ds_call2(SYM_Times, xr, logpow(xvar, j)));
}
/* x^a Cos/Sin[b Log x] (Log x)^j ; a,b borrowed */
static Expr* euler_trig(const char* xvar, const Expr* a, const Expr* b, int j, const char* which) {
    Expr* xa = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_new_symbol(xvar), expr_copy((Expr*)a) }, 2));
    Expr* blogx = ds_call2(SYM_Times, expr_copy((Expr*)b), ds_call1("Log", expr_new_symbol(xvar)));
    Expr* trig = ds_call1(which, blogx);
    Expr* prod = eval_and_free(ds_call2(SYM_Times, xa, trig));
    if (j == 0) return prod;
    return eval_and_free(ds_call2(SYM_Times, prod, logpow(xvar, j)));
}

Expr** dsolve_euler_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] < 1) return NULL;
    const char* xvar = P->ind_names[0];

    Expr** c; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return NULL;

    /* a_k = c_k / x^k must be free of x */
    Expr** a = malloc((size_t)(n + 1) * sizeof(Expr*));
    bool ok = true;
    for (int k = 0; k <= n; k++) {
        a[k] = eval_and_free(ds_call2(SYM_Times, expr_copy(c[k]),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_symbol(xvar), expr_new_integer(-k) }, 2)));
        if (ok && !ds_free_of(a[k], xvar)) ok = false;
    }
    if (ok && ds_is_zero(a[n])) ok = false;

    /* indicial polynomial P(r) = Σ a_k r(r-1)...(r-k+1) */
    Expr* indicial = NULL;
    const char* rv = intern_symbol("DSolve`r");
    if (ok) {
        Expr* ff = expr_new_integer(1);           /* falling factorial (r)_0 = 1 */
        Expr** terms = malloc((size_t)(n + 1) * sizeof(Expr*));
        for (int k = 0; k <= n; k++) {
            terms[k] = eval_and_free(ds_call2(SYM_Times, expr_copy(a[k]), expr_copy(ff)));
            ff = eval_and_free(ds_call2(SYM_Times, ff,
                     ds_call2(SYM_Subtract, expr_new_symbol(rv), expr_new_integer(k))));
        }
        indicial = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(n + 1)));
        free(terms);
        expr_free(ff);
    }

    Expr** out = NULL;
    if (indicial) {
        DSolveRoots R;
        if (dsolve_analyze_roots(indicial, rv, n, &R) && R.total == n) {
            bool* used = calloc(R.ndist, sizeof(bool));
            Expr** basis = malloc((size_t)n * sizeof(Expr*));
            size_t bc = 0;
            for (size_t i = 0; i < R.ndist; i++) {
                if (used[i]) continue;
                if (R.isreal[i]) {
                    for (int j = 0; j < R.mult[i]; j++) basis[bc++] = euler_real(xvar, R.roots[i], j);
                    used[i] = true;
                } else {
                    Expr* conj = eval_and_free(ds_call1("Conjugate", expr_copy(R.roots[i])));
                    long cc = -1;
                    for (size_t k = 0; k < R.ndist && cc < 0; k++) {
                        if (k == i || used[k]) continue;
                        Expr* diff = eval_and_free(ds_call2(SYM_Subtract, expr_copy(R.roots[k]), expr_copy(conj)));
                        if (ds_is_zero(diff)) cc = (long)k;
                        expr_free(diff);
                    }
                    expr_free(conj);
                    if (cc >= 0) {
                        used[i] = used[(size_t)cc] = true;
                        Expr* re = eval_and_free(ds_call1("Re", expr_copy(R.roots[i])));
                        for (int j = 0; j < R.mult[i]; j++) {
                            basis[bc++] = euler_trig(xvar, re, R.im[i], j, "Cos");
                            basis[bc++] = euler_trig(xvar, re, R.im[i], j, "Sin");
                        }
                        expr_free(re);
                    } else {
                        for (int j = 0; j < R.mult[i]; j++) basis[bc++] = euler_real(xvar, R.roots[i], j);
                        used[i] = true;
                    }
                }
            }
            if (bc == (size_t)n) {
                Expr** hterms = malloc((size_t)n * sizeof(Expr*));
                for (int k = 0; k < n; k++)
                    hterms[k] = ds_call2(SYM_Times, ds_const(k + 1), expr_copy(basis[k]));
                Expr* Hgen = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), hterms, (size_t)n));
                free(hterms);

                Expr* general = NULL;
                if (ds_is_zero(g)) {
                    general = Hgen;
                } else {
                    /* leading coefficient of y^(n) is c[n] = a_n x^n */
                    Expr* yp = dsolve_variation_of_parameters(basis, (size_t)n, g, c[n], xvar);
                    if (yp) general = eval_and_free(ds_call2(SYM_Plus, Hgen, yp));
                    else expr_free(Hgen);
                }
                if (general) { out = malloc(sizeof(Expr*)); out[0] = general; *nbranch = 1; }
            }
            for (size_t i = 0; i < bc; i++) expr_free(basis[i]);
            free(basis); free(used);
        }
        dsolve_roots_free(&R);
        expr_free(indicial);
    }

    for (int k = 0; k <= n; k++) { expr_free(a[k]); expr_free(c[k]); }
    free(a); free(c);
    expr_free(g);
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
        "a_n x^n y^(n) + ... + a_0 y == g(x) via the indicial polynomial (trial "
        "y = x^r); complex roots give Cos/Sin of Log[x], repeated roots give Log[x] powers.");
}
