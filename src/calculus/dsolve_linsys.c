/*
 * dsolve_linsys.c — DSolve`LinearFirstOrderSystem.
 *
 * Solves a first-order linear system with constant coefficients
 *     Y' == A Y + b(x),   A an n x n constant matrix,
 * by the eigen-decomposition of A.  Each real eigenvalue lambda (eigenvector v)
 * contributes C e^{lambda x} v; each complex-conjugate pair alpha +- beta i
 * (eigenvector p +- q i) contributes the two real solutions
 *     e^{alpha x}(cos(beta x) p - sin(beta x) q),
 *     e^{alpha x}(sin(beta x) p + cos(beta x) q).
 * A constant forcing b adds the particular solution -A^{-1} b.  Defective
 * (non-diagonalizable) matrices and non-constant forcing are left for later.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>
#include <stdio.h>

static Expr* ev(const char* head, Expr* a) { return eval_and_free(ds_call1(head, a)); }
static Expr* mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus, a, b)); }
static Expr* sub(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Subtract, a, b)); }

Expr** dsolve_linsys_solve(DSolveProblem* P) {
    size_t n = P->nfun;
    const char* xvar = P->ind_names[0];
    if (P->neq != n) return NULL;
    for (size_t i = 0; i < n; i++) if (P->max_order[i] != 1) return NULL;

    /* algebraic residuals: y_j'[x] -> Dsym[j], y_j[x] -> Ysym[j] */
    const char** Dn = malloc(n * sizeof(char*));
    const char** Yn = malloc(n * sizeof(char*));
    for (size_t j = 0; j < n; j++) {
        char b1[40], b2[40];
        snprintf(b1, sizeof(b1), "DSolve`sysD%zu", j); Dn[j] = intern_symbol(b1);
        snprintf(b2, sizeof(b2), "DSolve`sysY%zu", j); Yn[j] = intern_symbol(b2);
    }
    Expr** ralg = malloc(n * sizeof(Expr*));
    for (size_t e = 0; e < n; e++) {
        Expr* r = expr_copy(P->eq_residuals[e]);
        for (size_t j = 0; j < n; j++)
            r = ds_subst(r, ds_make_funcapp(P->fun_names[j], 1, xvar), expr_new_symbol(Dn[j]));
        for (size_t j = 0; j < n; j++)
            r = ds_subst(r, ds_make_funcapp(P->fun_names[j], 0, xvar), expr_new_symbol(Yn[j]));
        ralg[e] = r;
    }

    /* solve each equation for its (unique) leading derivative -> RHS[k] */
    Expr** RHS = calloc(n, sizeof(Expr*));
    bool ok = true;
    for (size_t e = 0; e < n && ok; e++) {
        long lead = -1;
        Expr* coeff = NULL;
        for (size_t j = 0; j < n; j++) {
            Expr* d = ds_d(expr_copy(ralg[e]), expr_new_symbol(Dn[j]));
            if (!ds_is_zero(d)) {
                if (lead >= 0) { ok = false; }            /* two derivatives in one eqn */
                lead = (long)j; if (coeff) expr_free(coeff); coeff = d;
            } else expr_free(d);
        }
        if (!ok || lead < 0 || !ds_free_of(coeff, xvar)) { if (coeff) expr_free(coeff); ok = false; break; }
        for (size_t j = 0; j < n; j++) if (!ds_free_of(coeff, Yn[j])) ok = false;
        if (RHS[lead]) ok = false;                        /* two eqns for same function */
        if (!ok) { expr_free(coeff); break; }
        /* RHS_lead = -(ralg|_{all D=0}) / coeff */
        Expr* r0 = expr_copy(ralg[e]);
        for (size_t j = 0; j < n; j++) r0 = ds_subst(r0, expr_new_symbol(Dn[j]), expr_new_integer(0));
        RHS[lead] = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                        ds_call2(SYM_Times, r0,
                            expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ coeff, expr_new_integer(-1) }, 2))));
    }
    for (size_t e = 0; e < n; e++) expr_free(ralg[e]);
    free(ralg);
    for (size_t k = 0; k < n && ok; k++) if (!RHS[k]) ok = false;

    /* A[i][j] = dRHS_i/dY_j (constant), b[i] = RHS_i|_{Y=0} */
    Expr* Amat = NULL; Expr* bvec = NULL; bool b_zero = true;
    if (ok) {
        Expr** rows = malloc(n * sizeof(Expr*));
        Expr** bs = malloc(n * sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr** cols = malloc(n * sizeof(Expr*));
            for (size_t j = 0; j < n; j++) {
                cols[j] = ds_d(expr_copy(RHS[i]), expr_new_symbol(Yn[j]));
                if (!ds_free_of(cols[j], xvar)) ok = false;
            }
            rows[i] = expr_new_function(expr_new_symbol(SYM_List), cols, n);
            free(cols);
            Expr* bi = expr_copy(RHS[i]);
            for (size_t j = 0; j < n; j++) bi = ds_subst(bi, expr_new_symbol(Yn[j]), expr_new_integer(0));
            if (!ds_is_zero(bi)) b_zero = false;
            bs[i] = bi;
        }
        Amat = expr_new_function(expr_new_symbol(SYM_List), rows, n); free(rows);
        bvec = expr_new_function(expr_new_symbol(SYM_List), bs, n); free(bs);
    }
    for (size_t k = 0; k < n; k++) if (RHS[k]) expr_free(RHS[k]);
    free(RHS); free(Dn);

    /* forcing must be constant (or zero) */
    if (ok && !b_zero && !ds_free_of(bvec, xvar)) ok = false;

    Expr** Y = NULL;
    if (ok) {
        Expr* lam = ds_delist(eval_and_free(ds_call1("Eigenvalues", expr_copy(Amat))));
        Expr* vec = ds_delist(eval_and_free(ds_call1("Eigenvectors", expr_copy(Amat))));
        if (head_is(lam, SYM_List) && head_is(vec, SYM_List)
            && lam->data.function.arg_count == n && vec->data.function.arg_count == n) {
            Y = malloc(n * sizeof(Expr*));
            for (size_t i = 0; i < n; i++) Y[i] = expr_new_integer(0);
            bool* used = calloc(n, sizeof(bool));
            int cnum = 0;
            bool built = true;
            for (size_t k = 0; k < n && built; k++) {
                if (used[k]) continue;
                Expr* lk = lam->data.function.args[k];
                Expr* vk = vec->data.function.args[k];
                if (!head_is(vk, SYM_List) || vk->data.function.arg_count != n) { built = false; break; }
                /* reject a zero eigenvector (defective) */
                bool zero_vec = true;
                for (size_t i = 0; i < n; i++) if (!ds_is_zero(vk->data.function.args[i])) zero_vec = false;
                if (zero_vec) { built = false; break; }

                Expr* imk = ev("Im", expr_copy(lk));
                bool real = ds_is_zero(imk); expr_free(imk);
                if (real) {
                    cnum++;
                    Expr* elx = ev("Exp", mul(expr_copy(lk), expr_new_symbol(xvar)));
                    for (size_t i = 0; i < n; i++)
                        Y[i] = add(Y[i], mul(ds_const(cnum),
                                    mul(expr_copy(elx), expr_copy(vk->data.function.args[i]))));
                    expr_free(elx);
                    used[k] = true;
                } else {
                    Expr* conj = ev("Conjugate", expr_copy(lk));
                    long cc = -1;
                    for (size_t m = 0; m < n && cc < 0; m++) {
                        if (m == k || used[m]) continue;
                        Expr* diff = sub(expr_copy(lam->data.function.args[m]), expr_copy(conj));
                        if (ds_is_zero(diff)) cc = (long)m;
                        expr_free(diff);
                    }
                    expr_free(conj);
                    if (cc < 0) { built = false; break; }
                    used[k] = used[(size_t)cc] = true;
                    Expr* al = ev("Re", expr_copy(lk));
                    Expr* be = ev("Im", expr_copy(lk));
                    Expr* eax = ev("Exp", mul(expr_copy(al), expr_new_symbol(xvar)));
                    Expr* cosbx = ev("Cos", mul(expr_copy(be), expr_new_symbol(xvar)));
                    Expr* sinbx = ev("Sin", mul(expr_copy(be), expr_new_symbol(xvar)));
                    int ca = ++cnum, cb = ++cnum;
                    for (size_t i = 0; i < n; i++) {
                        Expr* wi = vk->data.function.args[i];
                        Expr* p = ev("Re", expr_copy(wi));
                        Expr* q = ev("Im", expr_copy(wi));
                        Expr* s1 = mul(expr_copy(eax),
                                       sub(mul(expr_copy(cosbx), expr_copy(p)), mul(expr_copy(sinbx), expr_copy(q))));
                        Expr* s2 = mul(expr_copy(eax),
                                       add(mul(expr_copy(sinbx), expr_copy(p)), mul(expr_copy(cosbx), expr_copy(q))));
                        Y[i] = add(Y[i], add(mul(ds_const(ca), s1), mul(ds_const(cb), s2)));
                        expr_free(p); expr_free(q);
                    }
                    expr_free(al); expr_free(be); expr_free(eax); expr_free(cosbx); expr_free(sinbx);
                }
            }
            free(used);
            if (!built || cnum != (int)n) { for (size_t i = 0; i < n; i++) expr_free(Y[i]); free(Y); Y = NULL; }
        }
        expr_free(lam); expr_free(vec);

        /* particular solution for constant forcing: Y_p = -A^{-1} b */
        if (Y && !b_zero) {
            Expr* yp = ds_delist(eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                           ds_call2(SYM_Dot, ev("Inverse", expr_copy(Amat)), expr_copy(bvec)))));
            if (head_is(yp, SYM_List) && yp->data.function.arg_count == n)
                for (size_t i = 0; i < n; i++) Y[i] = add(Y[i], expr_copy(yp->data.function.args[i]));
            expr_free(yp);
        }
    }

    if (Amat) expr_free(Amat);
    if (bvec) expr_free(bvec);
    free(Yn);
    return Y;   /* NULL on decline */
}

void dsolve_linsys_init(void) {
    /* dispatched directly for nfun>1; no backtick builtin */
}
