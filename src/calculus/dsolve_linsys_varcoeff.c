/*
 * dsolve_linsys_varcoeff.c — DSolve`LinearSystemVarCoeff.
 *
 * Solves a genuinely-coupled, non-triangular, VARIABLE-coefficient first-order
 * linear system  Y' == A(x) Y + b(x)  for the tractable "scalar-factor" class
 *     A(x) == f(x) B,   B a constant matrix.
 * For this class the change of variable  t = Integrate[f, x]  turns the system
 * into  dY/dt == B Y  (chain rule: Y' = f dY/dt = f B Y), so the fundamental
 * matrix is exactly the constant-coefficient one evaluated at t:
 *     Phi(x) == e^{B t} == e^{B Integrate[f,x]}.
 * We therefore reuse the shared machinery from dsolve_linsys.c — extract A/b,
 * detect the scalar factor f, form B = A/f and t = Integrate[f,x], and hand
 * (B, t) to dsolve_linsys_assemble, which builds Phi = e^{Bt}, applies variation
 * of parameters for the forcing, and realifies.  Every branch is back-substitution
 * verified by dsolve_run_system, so a mis-detected class can only decline.
 *
 * Cascade slot: after DecoupleSystem / TriangularSystem (which claim the
 * decoupled / DAG-triangular variable-coefficient systems) and after
 * LinearFirstOrderSystem (constant A).  This method declines a constant A.
 *
 * Honest remaining gaps (documented, future work): the wider commutative-
 * antiderivative class where A(x) commutes with Integrate[A,x] but is NOT a
 * scalar multiple of a constant matrix, and the genuinely non-commuting case
 * (Floquet/Magnus).  Both need a symbolic MatrixExp of a variable matrix, which
 * is not available; this method deliberately covers only the scalar-factor class.
 */
#include "dsolve_linsys.h"
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

/* True iff every scalar entry of the n x n matrix `M` (a List of List) is free
 * of the variable `var`.  Robust per-entry test (D[entry,var] == 0), rather than
 * one whole-matrix zero test whose List result the zero-test would not decide. */
static bool matrix_free_of(const Expr* M, size_t n, const char* var) {
    if (!head_is((Expr*)M, SYM_List) || M->data.function.arg_count != n) return false;
    for (size_t i = 0; i < n; i++) {
        Expr* row = M->data.function.args[i];
        if (!head_is(row, SYM_List) || row->data.function.arg_count != n) return false;
        for (size_t j = 0; j < n; j++)
            if (!ds_free_of(row->data.function.args[j], var)) return false;
    }
    return true;
}

Expr** dsolve_linsys_varcoeff_solve(DSolveProblem* P) {
    Expr* A = NULL; Expr* b = NULL; bool b_zero = true;
    if (!dsolve_linsys_extract_Ab(P, &A, &b, &b_zero)) return NULL;
    size_t n = P->nfun;
    const char* xvar = P->ind_names[0];

    /* Constant A is LinearFirstOrderSystem's job (it ran first). */
    if (ds_free_of(A, xvar)) { expr_free(A); expr_free(b); return NULL; }

    /* Scalar factor f = the first nonzero entry of A (any works: for A = f B,
     * A_ij / A_pq = B_ij / B_pq is x-free, and the constant B_pq cancels against
     * t = Integrate[f,x] inside e^{B t}). */
    Expr* f = NULL;
    for (size_t i = 0; i < n && !f; i++) {
        Expr* row = A->data.function.args[i];
        for (size_t j = 0; j < n; j++) {
            Expr* aij = row->data.function.args[j];
            if (!ds_is_zero(aij)) { f = expr_copy(aij); break; }
        }
    }
    if (!f) { expr_free(A); expr_free(b); return NULL; }   /* A == 0 (handled above anyway) */

    /* B = Simplify[A / f]; require B constant (scalar-factor class), else decline. */
    Expr* finv = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy(f), expr_new_integer(-1) }, 2);
    Expr* B = ds_delist(ds_simplify(eval_and_free(ds_call2(SYM_Times, finv, expr_copy(A)))));
    if (!matrix_free_of(B, n, xvar)) {
        expr_free(A); expr_free(b); expr_free(f); expr_free(B); return NULL;
    }

    /* t = Integrate[f, x]; decline if not elementary. */
    Expr* t = ds_integrate(expr_copy(f), expr_new_symbol(xvar));
    if (ds_has_head(t, SYM_Integrate)) {
        expr_free(A); expr_free(b); expr_free(f); expr_free(B); expr_free(t); return NULL;
    }

    /* Phi = e^{B t}, then Y = tidy(Phi . (C + Integrate[Phi^{-1} b, x])). */
    Expr** Y = dsolve_linsys_assemble(B, t, xvar, b, b_zero, n);

    expr_free(A); expr_free(b); expr_free(f); expr_free(B); expr_free(t);
    return Y;   /* NULL on decline (Jordan / non-elementary VoP) */
}

static Expr* builtin_dsolve_linsys_varcoeff(Expr* res) {
    return dsolve_method_builtin_system(res, dsolve_linsys_varcoeff_solve);
}

void dsolve_linsys_varcoeff_init(void) {
    symtab_add_builtin("DSolve`LinearSystemVarCoeff", builtin_dsolve_linsys_varcoeff);
    symtab_get_def("DSolve`LinearSystemVarCoeff")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LinearSystemVarCoeff",
        "DSolve`LinearSystemVarCoeff[{eqns}, {y1, y2, ...}, x] solves a genuinely-coupled "
        "variable-coefficient linear system Y' == A(x) Y + b(x) of the scalar-factor class "
        "A(x) == f(x) B (B a constant matrix). The substitution t = Integrate[f, x] reduces "
        "it to the constant-coefficient system dY/dt == B Y, so Phi == e^{B t} is built from "
        "JordanDecomposition[B] (as in LinearFirstOrderSystem), forcing added by variation of "
        "parameters. Tried after DecoupleSystem / TriangularSystem / LinearFirstOrderSystem; "
        "declines a constant A (that is LinearFirstOrderSystem's job) and any A(x) not of the "
        "scalar-factor form. The wider commutative-antiderivative class and the non-commuting "
        "(Floquet/Magnus) case are future work.");
}
