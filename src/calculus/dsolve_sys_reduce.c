/*
 * dsolve_sys_reduce.c — DSolve`SystemReduce.
 *
 * Higher-order coupled linear system of ODEs -> first-order by STATE
 * AUGMENTATION, then reuse the constant-coefficient fundamental-matrix engine
 * (dsolve_linsys_assemble: Jordan -> e^{A t} -> variation of parameters).
 *
 * For a system of n functions u_1..u_n with per-function orders m_1..m_n, the
 * augmented state is
 *     Y = ( u_1, u_1', ..., u_1^(m_1-1),  u_2, ..., u_n^(m_n-1) ),   |Y| = N = Σ m_j.
 * Each state advances as  s_{j,k}' = s_{j,k+1}  (k < m_j-1); the top rows
 * s_{j,m_j-1}' = u_j^(m_j) come from solving the n ORIGINAL equations for the n
 * top derivatives u_j^(m_j) (the "leading matrix" L, which must be invertible):
 *     L · (u_j^(m_j))_j == -r(states, x).
 * That yields Y' == A Y + b(x); when A is constant this is the exact input the
 * LinearFirstOrderSystem engine already solves (defective/complex spectra,
 * VoP forcing), so this file only builds A, b and calls the shared assembler,
 * then reads back the u_j = s_{j,0} components.
 *
 * Scope (first cut): the leading matrix L invertible and A constant, i.e. the
 * constant-coefficient higher-order linear systems (Goode/Boyce §-style, the
 * bulk of the 12000.org §2.1.3 second-order-system set).  A variable A (variable-
 * coefficient higher-order system) or a singular L (a differential-algebraic
 * leading form) is declined here — the operator-determinant elimination fallback
 * (P1b) and the variable-coefficient extension are separate.
 *
 * Verification is the substrate's: dsolve_run_system back-substitutes each u_j
 * into the original higher-order equations, so a wrong A/b can only be dropped.
 */
#include "dsolve_common.h"
#include "dsolve_linsys.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../eval.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* head[a,b] built then evaluated; a,b consumed, result owned. */
static Expr* srr_eval2(const char* head, Expr* a, Expr* b) {
    return eval_and_free(ds_call2(head, a, b));
}

Expr** dsolve_sys_reduce_solve(DSolveProblem* P) {
    size_t n = P->nfun;
    if (P->is_pde || n < 2 || P->neq != n) return NULL;
    const char* xvar = P->ind_names[0];

    /* Gate: this method is for higher-order systems; pure first-order ones are
     * the existing LinearFirstOrderSystem's job. */
    bool higher = false;
    for (size_t j = 0; j < n; j++) if (P->max_order[j] >= 2) higher = true;
    if (!higher) return NULL;

    int* m = malloc(n * sizeof(int));
    size_t* off = malloc(n * sizeof(size_t));   /* start index of function j in Y */
    size_t N = 0;
    for (size_t j = 0; j < n; j++) { m[j] = P->max_order[j]; off[j] = N; N += (size_t)m[j]; }

    /* Interned symbols: state S[j][k] (k<m_j) and top D[j]. */
    const char** Ssym = malloc(N * sizeof(char*));         /* Y-indexed state syms  */
    const char** Dsym = malloc(n * sizeof(char*));         /* per-function top syms */
    for (size_t j = 0; j < n; j++) {
        char buf[48];
        for (int k = 0; k < m[j]; k++) {
            snprintf(buf, sizeof(buf), "DSolve`srS_%zu_%d", j, k);
            Ssym[off[j] + (size_t)k] = intern_symbol(buf);
        }
        snprintf(buf, sizeof(buf), "DSolve`srD_%zu", j);
        Dsym[j] = intern_symbol(buf);
    }

    /* Algebraic residuals: Derivative[k][u_j][x] -> (k<m_j ? S[j][k] : D[j]).
     * Substitute high order first so a lower-order funcapp never shadows a
     * higher-order Derivative subtree. */
    Expr** ralg = malloc(n * sizeof(Expr*));
    for (size_t e = 0; e < n; e++) {
        Expr* r = expr_copy(P->eq_residuals[e]);
        for (size_t j = 0; j < n; j++) {
            for (int k = m[j]; k >= 0; k--) {
                const char* to = (k < m[j]) ? Ssym[off[j] + (size_t)k] : Dsym[j];
                r = ds_subst(r, ds_make_funcapp(P->fun_names[j], k, xvar),
                             expr_new_symbol(to));
            }
        }
        ralg[e] = r;
    }

    bool ok = true;

    /* Leading matrix L[e][j] = dR_e/dD_j, and rhs = -(R_e | all D -> 0).
     * calloc so a partial fill on an early `ok = false` leaves NULL slots the
     * cleanup below skips rather than freeing garbage. */
    Expr** Lrows = calloc(n, sizeof(Expr*));
    Expr** rhs   = calloc(n, sizeof(Expr*));
    for (size_t e = 0; e < n && ok; e++) {
        Expr** cols = malloc(n * sizeof(Expr*));
        for (size_t j = 0; j < n; j++)
            cols[j] = ds_d(expr_copy(ralg[e]), expr_new_symbol(Dsym[j]));
        Lrows[e] = expr_new_function(expr_new_symbol(SYM_List), cols, n);
        free(cols);
        /* L must be free of every D (equations linear in the top derivatives). */
        for (size_t j = 0; j < n; j++) if (!ds_free_of(Lrows[e], Dsym[j])) ok = false;
        Expr* r0 = expr_copy(ralg[e]);
        for (size_t j = 0; j < n; j++) r0 = ds_subst(r0, expr_new_symbol(Dsym[j]), expr_new_integer(0));
        rhs[e] = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), r0));
    }
    for (size_t e = 0; e < n; e++) expr_free(ralg[e]);
    free(ralg);

    Expr** Dexpr = NULL;   /* solved top derivatives, length n */
    if (ok) {
        Expr* Lmat = expr_new_function(expr_new_symbol(SYM_List), Lrows, n);   /* owns Lrows elems */
        Expr* rvec = expr_new_function(expr_new_symbol(SYM_List), rhs,   n);   /* owns rhs elems   */
        free(Lrows); Lrows = NULL;   /* free the arrays; elements now owned by Lmat/rvec */
        free(rhs);   rhs   = NULL;
        Expr* sol = srr_eval2("LinearSolve", Lmat, rvec);
        /* sol must be a List of n entries, none still mentioning a D symbol or
         * LinearSolve (singular / unevaluated -> decline). */
        if (sol && sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol.name == SYM_List
            && sol->data.function.arg_count == n
            && !ds_has_head(sol, "LinearSolve")) {
            Dexpr = malloc(n * sizeof(Expr*));
            for (size_t j = 0; j < n; j++) Dexpr[j] = expr_copy(sol->data.function.args[j]);
        } else {
            ok = false;
        }
        if (sol) expr_free(sol);
    } else {
        if (Lrows) { for (size_t e = 0; e < n; e++) if (Lrows[e]) expr_free(Lrows[e]); free(Lrows); }
        if (rhs)   { for (size_t e = 0; e < n; e++) if (rhs[e])   expr_free(rhs[e]);   free(rhs); }
    }

    Expr** bodies = NULL;
    if (ok && Dexpr) {
        /* Companion rows: Y'_row.  Build A (N x N) and b (N) from
         *   row(s_{j,k}) = s_{j,k+1}     (k < m_j-1)
         *   row(s_{j,m_j-1}) = Dexpr[j]. */
        Expr** rowrhs = malloc(N * sizeof(Expr*));
        for (size_t j = 0; j < n; j++) {
            for (int k = 0; k < m[j]; k++) {
                size_t idx = off[j] + (size_t)k;
                rowrhs[idx] = (k < m[j] - 1)
                    ? expr_new_symbol(Ssym[off[j] + (size_t)k + 1])
                    : expr_copy(Dexpr[j]);
            }
        }
        /* A[row][col] = d(rowrhs)/dS_col ; b[row] = rowrhs | all S -> 0.
         * Require A constant (free of x) and free of every state (linear). */
        Expr** Arows = calloc(N, sizeof(Expr*));   /* calloc: partial fill on !aok */
        Expr** bvec  = calloc(N, sizeof(Expr*));
        bool b_zero = true;
        bool aok = true;
        for (size_t rr = 0; rr < N && aok; rr++) {
            Expr** cols = malloc(N * sizeof(Expr*));
            for (size_t cc = 0; cc < N; cc++)
                cols[cc] = ds_d(expr_copy(rowrhs[rr]), expr_new_symbol(Ssym[cc]));
            Arows[rr] = expr_new_function(expr_new_symbol(SYM_List), cols, N);
            free(cols);
            if (!ds_free_of(Arows[rr], xvar)) aok = false;              /* variable A */
            for (size_t cc = 0; cc < N && aok; cc++)
                if (!ds_free_of(Arows[rr], Ssym[cc])) aok = false;      /* nonlinear  */
            Expr* bi = expr_copy(rowrhs[rr]);
            for (size_t cc = 0; cc < N; cc++) bi = ds_subst(bi, expr_new_symbol(Ssym[cc]), expr_new_integer(0));
            bi = ds_simplify(bi);
            if (!ds_is_zero(bi)) b_zero = false;
            bvec[rr] = bi;
        }
        for (size_t rr = 0; rr < N; rr++) expr_free(rowrhs[rr]);
        free(rowrhs);

        if (aok) {
            Expr* Amat = expr_new_function(expr_new_symbol(SYM_List), Arows, N);
            Expr* bmat = expr_new_function(expr_new_symbol(SYM_List), bvec,  N);
            /* elements now owned by Amat/bmat; the arrays are freed at the end */
            Expr* tsym = expr_new_symbol(xvar);   /* assemble BORROWS t */
            Expr** aug = dsolve_linsys_assemble(Amat, tsym, xvar, bmat, b_zero, N);
            expr_free(tsym);
            expr_free(Amat); expr_free(bmat);
            if (aug) {
                /* u_j = s_{j,0} = aug[off[j]]. */
                bodies = malloc(n * sizeof(Expr*));
                for (size_t j = 0; j < n; j++) bodies[j] = expr_copy(aug[off[j]]);
                for (size_t i = 0; i < N; i++) expr_free(aug[i]);
                free(aug);
            }
        } else {
            for (size_t rr = 0; rr < N; rr++) { if (Arows && Arows[rr]) expr_free(Arows[rr]);
                                                if (bvec  && bvec[rr])  expr_free(bvec[rr]); }
        }
        if (Arows) free(Arows);
        if (bvec)  free(bvec);
    }

    if (Dexpr) { for (size_t j = 0; j < n; j++) expr_free(Dexpr[j]); free(Dexpr); }
    free(m); free(off); free(Ssym); free(Dsym);
    return bodies;
}

static Expr* builtin_dsolve_sys_reduce(Expr* res) {
    return dsolve_method_builtin_system(res, dsolve_sys_reduce_solve);
}

void dsolve_sys_reduce_init(void) {
    symtab_add_builtin("DSolve`SystemReduce", builtin_dsolve_sys_reduce);
    symtab_get_def("DSolve`SystemReduce")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SystemReduce",
        "DSolve`SystemReduce[eqns, {u1,...}, x] solves a higher-order coupled "
        "linear ODE system by state augmentation to a first-order system "
        "Y' == A Y + b, reusing the constant-coefficient fundamental-matrix "
        "engine; declines a variable A or a singular leading matrix.");
}
