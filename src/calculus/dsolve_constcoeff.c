/*
 * dsolve_constcoeff.c — DSolve`LinearConstantCoefficients.
 *
 * Solves the linear constant-coefficient ODE
 *     a_n y^(n) + ... + a_1 y' + a_0 y == g(x)   (a_k constant)
 * by the characteristic polynomial  char(L) = sum a_k L^k:
 *   - each root r of multiplicity m contributes  x^j E^(r x),  j = 0..m-1;
 *   - a complex-conjugate pair a +- b i contributes the real forms
 *     x^j E^(a x) Cos[b x] and x^j E^(a x) Sin[b x].
 * The inhomogeneous part is handled by variation of parameters over the
 * fundamental set, using the Wronskian (Cramer's rule + Integrate).
 *
 * Coefficients are read off the residual R by replacing y^(k)[x] with plain
 * symbols D_k: a_k = R_{D_k} must be free of x (constant) and of every D_j
 * (linear); g = -(R with all D_k -> 0).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>
#include <stdio.h>

/* x^j (j >= 0) */
static Expr* xpow(const char* xvar, int j) {
    return expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_new_symbol(xvar), expr_new_integer(j) }, 2);
}
/* Exp[coef * x]; coef consumed */
static Expr* exp_lin(Expr* coef, const char* xvar) {
    return eval_and_free(ds_call1("Exp", ds_call2(SYM_Times, coef, expr_new_symbol(xvar))));
}
/* x^j E^(r x); r borrowed */
static Expr* basis_exp(const char* xvar, int j, const Expr* r) {
    Expr* e = exp_lin(expr_copy((Expr*)r), xvar);
    if (j == 0) return e;
    return eval_and_free(ds_call2(SYM_Times, xpow(xvar, j), e));
}
/* x^j E^(a x) Cos/Sin[b x]; a,b borrowed */
static Expr* basis_trig(const char* xvar, int j, const Expr* a, const Expr* b, const char* trig) {
    Expr* e = exp_lin(expr_copy((Expr*)a), xvar);
    Expr* t = ds_call1(trig, ds_call2(SYM_Times, expr_copy((Expr*)b), expr_new_symbol(xvar)));
    Expr* prod = eval_and_free(ds_call2(SYM_Times, e, t));
    if (j == 0) return prod;
    return eval_and_free(ds_call2(SYM_Times, xpow(xvar, j), prod));
}

/* Build the fundamental set for charpoly; returns malloc'd array of *count basis
 * functions (in xvar), or NULL if the roots do not account for degree n. */
static Expr** homog_basis(const Expr* charpoly, const char* lam, const char* xvar,
                          int n, size_t* count) {
    *count = 0;
    DSolveRoots R;
    if (!dsolve_analyze_roots(charpoly, lam, n, &R)) return NULL;
    Expr** basis = NULL;
    if (R.total == n) {
        bool* used = calloc(R.ndist, sizeof(bool));
        basis = malloc((size_t)n * sizeof(Expr*));
        size_t bc = 0;
        for (size_t i = 0; i < R.ndist; i++) {
            if (used[i]) continue;
            if (R.isreal[i]) {
                for (int j = 0; j < R.mult[i]; j++) basis[bc++] = basis_exp(xvar, j, R.roots[i]);
                used[i] = true;
            } else {
                Expr* conj = eval_and_free(ds_call1("Conjugate", expr_copy(R.roots[i])));
                long c = -1;
                for (size_t k = 0; k < R.ndist && c < 0; k++) {
                    if (k == i || used[k]) continue;
                    Expr* diff = eval_and_free(ds_call2(SYM_Subtract, expr_copy(R.roots[k]), expr_copy(conj)));
                    if (ds_is_zero(diff)) c = (long)k;
                    expr_free(diff);
                }
                expr_free(conj);
                if (c >= 0) {
                    used[i] = used[(size_t)c] = true;
                    Expr* a = eval_and_free(ds_call1("Re", expr_copy(R.roots[i])));
                    for (int j = 0; j < R.mult[i]; j++) {
                        basis[bc++] = basis_trig(xvar, j, a, R.im[i], "Cos");
                        basis[bc++] = basis_trig(xvar, j, a, R.im[i], "Sin");
                    }
                    expr_free(a);
                } else {
                    for (int j = 0; j < R.mult[i]; j++) basis[bc++] = basis_exp(xvar, j, R.roots[i]);
                    used[i] = true;
                }
            }
        }
        *count = bc;
        free(used);
    }
    dsolve_roots_free(&R);
    return basis;
}

Expr** dsolve_constcoeff_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    const char* xvar = P->ind_names[0];

    Expr** a; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &a, &g, &n)) return NULL;

    /* constant coefficients: every a_k free of x, leading a_n nonzero */
    bool ok = true;
    for (int k = 0; k <= n && ok; k++) if (!ds_free_of(a[k], xvar)) ok = false;
    if (ok && ds_is_zero(a[n])) ok = false;

    Expr* charpoly = NULL;
    if (ok) {
        const char* lam = intern_symbol("DSolve`lam");
        Expr** ch = malloc((size_t)(n + 1) * sizeof(Expr*));
        for (int k = 0; k <= n; k++)
            ch[k] = ds_call2(SYM_Times, expr_copy(a[k]),
                        expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_new_symbol(lam), expr_new_integer(k) }, 2));
        charpoly = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), ch, (size_t)(n + 1)));
        free(ch);
    }

    Expr** out = NULL;
    if (charpoly) {
        const char* lam = intern_symbol("DSolve`lam");
        size_t nb = 0;
        Expr** basis = homog_basis(charpoly, lam, xvar, n, &nb);
        if (basis && nb == (size_t)n) {
            Expr** hterms = malloc((size_t)n * sizeof(Expr*));
            for (int k = 0; k < n; k++)
                hterms[k] = ds_call2(SYM_Times, ds_const(k + 1), expr_copy(basis[k]));
            Expr* Hgen = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), hterms, (size_t)n));
            free(hterms);

            Expr* general = NULL;
            if (ds_is_zero(g)) {
                general = Hgen;
            } else {
                Expr* yp = dsolve_variation_of_parameters(basis, (size_t)n, g, a[n], xvar);
                if (yp) general = eval_and_free(ds_call2(SYM_Plus, Hgen, yp));
                else expr_free(Hgen);
            }
            if (general) { out = malloc(sizeof(Expr*)); out[0] = general; *nbranch = 1; }
        }
        if (basis) { for (size_t i = 0; i < nb; i++) expr_free(basis[i]); free(basis); }
        expr_free(charpoly);
    }

    for (int k = 0; k <= n; k++) expr_free(a[k]);
    free(a);
    expr_free(g);
    return out;
}

static Expr* builtin_dsolve_constcoeff(Expr* res) {
    return dsolve_method_builtin(res, dsolve_constcoeff_try);
}

void dsolve_constcoeff_init(void) {
    symtab_add_builtin("DSolve`LinearConstantCoefficients", builtin_dsolve_constcoeff);
    symtab_get_def("DSolve`LinearConstantCoefficients")->attributes
        |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LinearConstantCoefficients",
        "DSolve`LinearConstantCoefficients[eqn, y, x] solves a_n y^(n) + ... + a_0 y "
        "== g(x) with constant a_k via the characteristic polynomial (complex roots "
        "give Cos/Sin terms) plus variation of parameters for the forcing g(x).");
}
