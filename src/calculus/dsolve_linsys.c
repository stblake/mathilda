/*
 * dsolve_linsys.c — DSolve`LinearFirstOrderSystem + shared fundamental-matrix
 *                   machinery (declared in dsolve_linsys.h).
 *
 * Solves a first-order linear system with constant coefficients
 *     Y' == A Y + b(x),   A an n x n constant matrix,
 * for ANY constant A (diagonalizable OR defective) via the fundamental matrix
 * Phi(x) = e^{Ax}, assembled from the Jordan decomposition A = S J S^{-1}:
 *
 *     e^{Jx} = e^{Dx} e^{Nx},   J = D + N  (D diagonal eigenvalues,
 *                                           N strictly-upper nilpotent),
 *     e^{Dx} = DiagonalMatrix[Exp[lambda_i x]],
 *     e^{Nx} = sum_{m<n} N^m x^m / m!      (finite: N is nilpotent, N^n == 0),
 *     Phi    = S e^{Jx} S^{-1}.
 *
 * D and N commute (each Jordan block is scalar on its diagonal), so the split is
 * exact.  The homogeneous solution is Y = Phi . {C[1],...,C[n]}; a forcing b(x)
 * is added by variation of parameters, Y = Phi . (C + Integrate[Phi^{-1} b, x]),
 * which subsumes the older -A^{-1} b particular and stays valid when A is
 * singular.  Each body is realified/tidied by
 * Simplify[ComplexExpand[.]] //. (Cosh[a]+Sinh[a] -> E^a): complex eigenvalues
 * collapse to e^{alpha x} Cos/Sin[beta x], and the Cosh/Sinh form Simplify
 * introduces for a repeated real eigenvalue is folded back to E^{lambda x}.
 *
 * Symbolic MatrixExp is currently inert, so Phi is built from Jordan directly.
 *
 * The extraction (Y' == A Y + b), the matrix exponential e^{Mt}, the realifier,
 * and the "Phi . (C + VoP)" assembler are exposed via dsolve_linsys.h so the
 * variable-coefficient sibling (dsolve_linsys_varcoeff.c) reuses them: for a
 * scalar-factor system Y' == f(x) B Y the change of variable t = Integrate[f, x]
 * turns it into dY/dt == B Y, so Phi = e^{B Integrate[f,x]} is exactly this same
 * assembler with `t` = the antiderivative instead of the bare symbol x.
 */
#include "dsolve_linsys.h"
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../internal.h"
#include <stdlib.h>
#include <stdio.h>

static Expr* ev(const char* head, Expr* a) { return eval_and_free(ds_call1(head, a)); }
static Expr* mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus, a, b)); }

/* The rewrite rule  Cosh[a_] + Sinh[a_] :> E^a  (RuleDelayed), used to fold the
 * hyperbolic form Simplify prefers for a real repeated eigenvalue back to a
 * plain exponential without disturbing genuine Cos/Sin from complex spectra. */
static Expr* cosh_sinh_rule(void) {
    Expr* blank1 = expr_new_function(expr_new_symbol("Blank"), NULL, 0);
    Expr* blank2 = expr_new_function(expr_new_symbol("Blank"), NULL, 0);
    Expr* pa1 = expr_new_function(expr_new_symbol("Pattern"),
                    (Expr*[]){ expr_new_symbol("a"), blank1 }, 2);
    Expr* pa2 = expr_new_function(expr_new_symbol("Pattern"),
                    (Expr*[]){ expr_new_symbol("a"), blank2 }, 2);
    Expr* lhs = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
                    expr_new_function(expr_new_symbol("Cosh"), (Expr*[]){ pa1 }, 1),
                    expr_new_function(expr_new_symbol("Sinh"), (Expr*[]){ pa2 }, 1) }, 2);
    Expr* rhs = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_new_symbol("E"), expr_new_symbol("a") }, 2);
    return expr_new_function(expr_new_symbol(SYM_RuleDelayed), (Expr*[]){ lhs, rhs }, 2);
}

/* True iff `e` contains the imaginary unit (a Complex[...] node) anywhere. */
static bool expr_contains_complex(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is((Expr*)e, SYM_Complex)) return true;
    if (expr_contains_complex(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (expr_contains_complex(e->data.function.args[i])) return true;
    return false;
}

/* Realify + tidy a solution body; consumes `body`, returns owned.
 * ComplexExpand is applied ONLY when the body actually carries the imaginary unit
 * (a complex spectrum needing e^{a+ib x} -> e^{ax}Cos/Sin[bx]).  On an already-real
 * body it is skipped: ComplexExpand assumes real-but-possibly-negative variables and
 * would gratuitously split Log[x] -> Log[Abs[x]] + I Arg[x] (which arises in the
 * variable-coefficient forcing integral, e.g. Integrate[x^{-1}, x] = Log[x]). */
Expr* dsolve_linsys_tidy(Expr* body) {
    Expr* ce   = expr_contains_complex(body)
                     ? eval_and_free(ds_call1("ComplexExpand", body)) : body;
    Expr* si   = eval_and_free(ds_call1("Simplify", ce));
    Expr* rule = cosh_sinh_rule();
    return eval_and_free(internal_replace_repeated((Expr*[]){ si, rule }, 2));
}

/* e^{M t} for M = S J S^{-1}, given the Jordan factors S (change of basis) and
 * J (Jordan form, upper-triangular with the eigenvalues on the diagonal).
 * S, J, t are borrowed; the returned n x n matrix (List of Lists) is owned. */
Expr* dsolve_linsys_matexp(Expr* S, Expr* J, Expr* t, size_t n) {
    /* e^{Dt} (diagonal Exp[J_ii t]) and N = J with the diagonal zeroed */
    Expr** drows = malloc(n * sizeof(Expr*));
    Expr** nrows = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        Expr** dcol = malloc(n * sizeof(Expr*));
        Expr** ncol = malloc(n * sizeof(Expr*));
        Expr* Ji = J->data.function.args[i];
        for (size_t j = 0; j < n; j++) {
            Expr* Jij = Ji->data.function.args[j];
            if (i == j) {
                dcol[j] = ev("Exp", mul(expr_copy(Jij), expr_copy(t)));
                ncol[j] = expr_new_integer(0);
            } else {
                dcol[j] = expr_new_integer(0);
                ncol[j] = expr_copy(Jij);
            }
        }
        drows[i] = expr_new_function(expr_new_symbol(SYM_List), dcol, n); free(dcol);
        nrows[i] = expr_new_function(expr_new_symbol(SYM_List), ncol, n); free(ncol);
    }
    Expr* eD   = expr_new_function(expr_new_symbol(SYM_List), drows, n); free(drows);
    Expr* Nmat = expr_new_function(expr_new_symbol(SYM_List), nrows, n); free(nrows);

    /* e^{Nt} = sum_{m=0}^{n-1} MatrixPower[N,m] t^m / m!  (finite, N nilpotent) */
    Expr* eN = NULL;
    long fact = 1;                                /* fact == m! at loop head */
    for (size_t m = 0; m < n; m++) {
        if (m >= 2) fact *= (long)m;
        Expr* Npow = eval_and_free(ds_call2("MatrixPower",
                         expr_copy(Nmat), expr_new_integer((long)m)));
        Expr* scal = mul(expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_copy(t), expr_new_integer((long)m) }, 2),
                         expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_new_integer(fact), expr_new_integer(-1) }, 2));
        Expr* term = mul(scal, Npow);
        eN = eN ? add(eN, term) : term;
    }
    expr_free(Nmat);

    Expr* eJt  = eval_and_free(ds_call2(SYM_Dot, eD, eN));
    Expr* SeJt = eval_and_free(ds_call2(SYM_Dot, expr_copy(S), eJt));
    return eval_and_free(ds_call2(SYM_Dot, SeJt,
                             eval_and_free(ds_call1("Inverse", expr_copy(S)))));
}

/* Extract Y' == A Y + b from the parsed square system.  See dsolve_linsys.h. */
bool dsolve_linsys_extract_Ab(DSolveProblem* P, Expr** Aout, Expr** bout, bool* b_zero_out) {
    size_t n = P->nfun;
    const char* xvar = P->ind_names[0];
    if (P->neq != n) return false;
    for (size_t i = 0; i < n; i++) if (P->max_order[i] != 1) return false;

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

    /* A[i][j] = dRHS_i/dY_j, b[i] = RHS_i|_{Y=0}.  NB: no x-dependence guard here —
     * A/b may depend on x; the constant-A decision is the caller's. */
    Expr* Amat = NULL; Expr* bvec = NULL; bool b_zero = true;
    if (ok) {
        Expr** rows = malloc(n * sizeof(Expr*));
        Expr** bs = malloc(n * sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr** cols = malloc(n * sizeof(Expr*));
            for (size_t j = 0; j < n; j++)
                cols[j] = ds_d(expr_copy(RHS[i]), expr_new_symbol(Yn[j]));
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
    free(RHS); free(Dn); free(Yn);

    if (!ok) { if (Amat) expr_free(Amat); if (bvec) expr_free(bvec); return false; }
    *Aout = Amat; *bout = bvec; *b_zero_out = b_zero;
    return true;
}

/* Y = tidy( e^{M t} . (C + Integrate[e^{-M t} b, x]) ).  See dsolve_linsys.h. */
Expr** dsolve_linsys_assemble(Expr* M, Expr* t, const char* xvar,
                              Expr* b, bool b_zero, size_t n) {
    Expr** Y = NULL;
    Expr* jd = ds_delist(eval_and_free(ds_call1("JordanDecomposition", expr_copy(M))));
    if (head_is(jd, SYM_List) && jd->data.function.arg_count == 2) {
        Expr* S = jd->data.function.args[0];
        Expr* J = jd->data.function.args[1];
        bool shape = head_is(S, SYM_List) && head_is(J, SYM_List)
            && S->data.function.arg_count == n && J->data.function.arg_count == n;
        for (size_t i = 0; shape && i < n; i++)
            if (!head_is(J->data.function.args[i], SYM_List)
                || J->data.function.args[i]->data.function.arg_count != n) shape = false;
        if (shape) {
            Expr* Phi = dsolve_linsys_matexp(S, J, t, n);

            /* rhs = {C[1..n]} (+ variation of parameters if forced) */
            Expr** rc = malloc(n * sizeof(Expr*));
            for (size_t i = 0; i < n; i++) rc[i] = ds_const((int)i + 1);
            bool force_ok = true;
            if (!b_zero) {
                Expr* negt   = mul(expr_new_integer(-1), expr_copy(t));
                Expr* PhiInv = dsolve_linsys_matexp(S, J, negt, n);
                expr_free(negt);
                Expr* integ  = ds_delist(eval_and_free(
                                   ds_call2(SYM_Dot, PhiInv, expr_copy(b))));
                if (head_is(integ, SYM_List) && integ->data.function.arg_count == n) {
                    for (size_t i = 0; i < n && force_ok; i++) {
                        Expr* anti = ds_integrate(expr_copy(integ->data.function.args[i]),
                                                  expr_new_symbol(xvar));
                        if (ds_has_head(anti, SYM_Integrate)) { expr_free(anti); force_ok = false; }
                        else rc[i] = add(rc[i], anti);
                    }
                } else force_ok = false;
                expr_free(integ);
            }

            if (force_ok) {
                Expr* rhs = expr_new_function(expr_new_symbol(SYM_List), rc, n);
                free(rc);
                Expr* Yvec = ds_delist(eval_and_free(ds_call2(SYM_Dot, expr_copy(Phi), rhs)));
                if (head_is(Yvec, SYM_List) && Yvec->data.function.arg_count == n) {
                    Y = malloc(n * sizeof(Expr*));
                    for (size_t i = 0; i < n; i++)
                        Y[i] = dsolve_linsys_tidy(expr_copy(Yvec->data.function.args[i]));
                }
                expr_free(Yvec);
            } else {
                for (size_t i = 0; i < n; i++) expr_free(rc[i]);
                free(rc);
            }
            expr_free(Phi);
        }
    }
    expr_free(jd);
    return Y;
}

Expr** dsolve_linsys_solve(DSolveProblem* P) {
    Expr* A = NULL; Expr* b = NULL; bool b_zero = true;
    if (!dsolve_linsys_extract_Ab(P, &A, &b, &b_zero)) return NULL;
    size_t n = P->nfun;
    const char* xvar = P->ind_names[0];

    /* This method claims only CONSTANT A with constant (or zero) forcing; a
     * variable-coefficient A(x) is DSolve`LinearSystemVarCoeff's job. */
    if (!ds_free_of(A, xvar) || (!b_zero && !ds_free_of(b, xvar))) {
        expr_free(A); expr_free(b); return NULL;
    }

    Expr* t = expr_new_symbol(xvar);
    Expr** Y = dsolve_linsys_assemble(A, t, xvar, b, b_zero, n);
    expr_free(t);
    expr_free(A); expr_free(b);
    return Y;   /* NULL on decline */
}

static Expr* builtin_dsolve_linsys(Expr* res) {
    return dsolve_method_builtin_system(res, dsolve_linsys_solve);
}

void dsolve_linsys_init(void) {
    symtab_add_builtin("DSolve`LinearFirstOrderSystem", builtin_dsolve_linsys);
    symtab_get_def("DSolve`LinearFirstOrderSystem")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LinearFirstOrderSystem",
        "DSolve`LinearFirstOrderSystem[{eqns}, {y1, y2, ...}, x] solves a "
        "constant-coefficient linear system Y' == A Y + b(x) for ANY constant matrix A "
        "(diagonalizable OR defective). The fundamental matrix Phi = e^{Ax} = S.e^{Jx}.S^{-1} "
        "is built from JordanDecomposition (diagonalizable -> C e^{lambda x} v; defective -> "
        "x^k e^{lambda x} generalized-eigenvector terms; complex pairs -> real "
        "e^{a x}Cos/Sin[b x] via ComplexExpand). Forcing b(x) by variation of parameters "
        "Phi.(C + Integral[Phi^{-1} b]), which stays valid for singular A. The general "
        "backstop for irreducibly-coupled constant systems; tried after DecoupleSystem "
        "and TriangularSystem.");
}
