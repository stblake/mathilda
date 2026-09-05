/*
 * dsolve_linsys_commutative.c — DSolve`LinearSystemCommutative.
 *
 * Variable-coefficient 2x2 linear systems  Y' == A(t) Y + b(t)  of the
 * COMMUTATIVE class  A(t) == a(t) I + b(t) K0  with K0 a CONSTANT matrix — the
 * case [A, ∫A] == 0, where the fundamental matrix is the ordinary exponential
 *     Phi = Exp[∫A] = Exp[∫a] · Exp[(∫b) K0].
 * For a constant traceless 2x2 K0 with mu^2 = -det K0,
 *     Exp[s K0] = Cosh[mu s] I + (Sinh[mu s]/mu) K0        (mu != 0)
 *               = I + s K0                                  (mu == 0, nilpotent)
 * and complex mu (det K0 > 0) realifies to the rotation Cos/Sin form.  This
 * subsumes the scalar-factor class (a proportional to b) and adds the genuine
 * rotation/decay families the scalar-factor reduction misses:
 *     {x' = -x + t y, y' = t x - y}                     (K0 symmetric, cosh/sinh)
 *     {x' = x Cos t - y Sin t, y' = x Sin t + y Cos t}  (K0 rotation, cos/sin)
 *     {x' = x/t + y, y' = -x + y/t}                      (a = 1/t, rotation)
 *     {(t^2+1) x' = -t x + y, (t^2+1) y' = -x - t y}
 *
 * Symbolic Eigenvectors of a variable matrix come back in non-constant Sign[t]
 * scaled forms, so this decomposition (trace for a, one scalar factor for the
 * K0 remainder) is the robust route rather than a symbolic eigendecomposition.
 *
 * Runs after LinearSystemVarCoeff (which claims the scalar-factor A == f B first)
 * and declines constant A / genuinely non-commutative A (K0 not constant).
 * Forcing by variation of parameters; verification is the substrate's.
 */
#include "dsolve_common.h"
#include "dsolve_linsys.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

static Expr* cm_mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* cm_add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus, a, b)); }
static Expr* cm_get(const Expr* M, int i, int j) {   /* copy of the (i,j) entry */
    return expr_copy(M->data.function.args[i]->data.function.args[j]);
}
static Expr* cm_row2(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ a, b }, 2);
}

Expr** dsolve_linsys_commutative_solve(DSolveProblem* P) {
    if (P->nfun != 2 || P->neq != 2) return NULL;
    const char* xvar = P->ind_names[0];

    Expr* A = NULL; Expr* b = NULL; bool b_zero = true;
    if (!dsolve_linsys_extract_Ab(P, &A, &b, &b_zero)) return NULL;

    /* constant A is LinearFirstOrderSystem's job */
    if (ds_free_of(A, xvar)) { expr_free(A); expr_free(b); return NULL; }

    bool ok = true;
    Expr* a   = ds_simplify(cm_mul(cm_add(cm_get(A,0,0), cm_get(A,1,1)),
                                   ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1))));
    /* A0 = A - a I */
    Expr* A0[2][2];
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
        Expr* e = cm_get(A, i, j);
        if (i == j) e = cm_add(e, cm_mul(expr_new_integer(-1), expr_copy(a)));
        A0[i][j] = ds_simplify(e);
    }
    /* first nonzero A0 entry -> b(t); K0 = A0/b, must be CONSTANT */
    Expr* bcoef = NULL;
    for (int i = 0; i < 2 && !bcoef; i++) for (int j = 0; j < 2 && !bcoef; j++)
        if (!ds_is_zero(A0[i][j])) bcoef = expr_copy(A0[i][j]);
    Expr* K0[2][2]; for (int i=0;i<2;i++) for (int j=0;j<2;j++) K0[i][j]=NULL;
    if (!bcoef) ok = false;
    if (ok) {
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
            K0[i][j] = ds_simplify(cm_mul(expr_copy(A0[i][j]),
                           ds_call2(SYM_Power, expr_copy(bcoef), expr_new_integer(-1))));
            if (!ds_free_of(K0[i][j], xvar)) ok = false;
        }
    }
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) expr_free(A0[i][j]);

    /* s = ∫b, inta = ∫a, ea = Exp[inta]; both integrals must be elementary */
    Expr* s = NULL; Expr* ea = NULL;
    if (ok) {
        s = ds_integrate(expr_copy(bcoef), expr_new_symbol(xvar));
        if (ds_has_head(s, SYM_Integrate)) ok = false;
    }
    if (ok) {
        Expr* inta = ds_integrate(expr_copy(a), expr_new_symbol(xvar));
        if (ds_has_head(inta, SYM_Integrate)) { expr_free(inta); ok = false; }
        else ea = eval_and_free(ds_call1("Exp", inta));
    }

    Expr** bodies = NULL;
    if (ok) {
        /* mu^2 = -det K0 */
        Expr* detK0 = ds_simplify(cm_add(
            cm_mul(expr_copy(K0[0][0]), expr_copy(K0[1][1])),
            cm_mul(expr_new_integer(-1), cm_mul(expr_copy(K0[0][1]), expr_copy(K0[1][0])))));
        Expr* mu2 = ds_simplify(cm_mul(expr_new_integer(-1), detK0));

        /* ExpK = I + s K0 (nilpotent) or Cosh[mu s] I + Sinh[mu s]/mu K0 */
        Expr* diagT; Expr* offT;   /* diagonal-add and K0-scale factors */
        if (ds_is_zero(mu2)) {
            diagT = expr_new_integer(1);
            offT  = expr_copy(s);
            expr_free(mu2);
        } else {
            Expr* mu = eval_and_free(ds_call1("Sqrt", mu2));
            Expr* ms = cm_mul(expr_copy(mu), expr_copy(s));
            diagT = eval_and_free(ds_call1("Cosh", expr_copy(ms)));
            offT  = cm_mul(eval_and_free(ds_call1("Sinh", ms)),
                           ds_call2(SYM_Power, mu, expr_new_integer(-1)));   /* Sinh[mu s]/mu */
        }
        /* Phi[i][j] = ea * ( (i==j) diagT + offT K0[i][j] ) */
        Expr* Prow[2];
        for (int i = 0; i < 2; i++) {
            Expr* col[2];
            for (int j = 0; j < 2; j++) {
                Expr* term = cm_mul(expr_copy(offT), expr_copy(K0[i][j]));
                if (i == j) term = cm_add(term, expr_copy(diagT));
                col[j] = cm_mul(expr_copy(ea), term);
            }
            Prow[i] = cm_row2(col[0], col[1]);
        }
        Expr* Phi = cm_row2(Prow[0], Prow[1]);
        expr_free(diagT); expr_free(offT);

        /* rhs = {C[1],C[2]} (+ VoP Integrate[Phi^-1 b]) */
        Expr* c1 = ds_const(1); Expr* c2 = ds_const(2);
        Expr* rhs = cm_row2(c1, c2);
        bool force_ok = true;
        if (!b_zero) {
            Expr* PhiInv = eval_and_free(ds_call1("Inverse", expr_copy(Phi)));
            Expr* prod   = ds_delist(eval_and_free(ds_call2(SYM_Dot, PhiInv, expr_copy(b))));
            if (prod->type == EXPR_FUNCTION && prod->data.function.arg_count == 2) {
                Expr* ivec[2] = { NULL, NULL };
                for (int i = 0; i < 2 && force_ok; i++) {
                    Expr* ig = eval_and_free(ds_call1("Expand",
                                   expr_copy(prod->data.function.args[i])));
                    Expr* anti = ds_integrate(ig, expr_new_symbol(xvar));
                    if (ds_has_head(anti, SYM_Integrate)) { expr_free(anti); force_ok = false; }
                    else ivec[i] = anti;
                }
                if (force_ok) {
                    Expr* newrhs = cm_row2(
                        cm_add(expr_copy(rhs->data.function.args[0]), ivec[0]),
                        cm_add(expr_copy(rhs->data.function.args[1]), ivec[1]));
                    expr_free(rhs); rhs = newrhs;
                }
            } else force_ok = false;
            expr_free(prod);
        }

        if (force_ok) {
            Expr* Y = ds_delist(eval_and_free(ds_call2(SYM_Dot, expr_copy(Phi), rhs)));
            if (Y->type == EXPR_FUNCTION && Y->data.function.arg_count == 2) {
                bodies = malloc(2 * sizeof(Expr*));
                bodies[0] = dsolve_linsys_tidy(expr_copy(Y->data.function.args[0]));
                bodies[1] = dsolve_linsys_tidy(expr_copy(Y->data.function.args[1]));
            }
            expr_free(Y);
        } else expr_free(rhs);
        expr_free(Phi);
    }

    expr_free(a); if (s) expr_free(s); if (ea) expr_free(ea); if (bcoef) expr_free(bcoef);
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) if (K0[i][j]) expr_free(K0[i][j]);
    expr_free(A); expr_free(b);
    return bodies;
}

static Expr* builtin_dsolve_linsys_commutative(Expr* res) {
    return dsolve_method_builtin_system(res, dsolve_linsys_commutative_solve);
}

void dsolve_linsys_commutative_init(void) {
    symtab_add_builtin("DSolve`LinearSystemCommutative", builtin_dsolve_linsys_commutative);
    symtab_get_def("DSolve`LinearSystemCommutative")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LinearSystemCommutative",
        "DSolve`LinearSystemCommutative[eqns, {x,y}, t] solves a 2x2 variable-"
        "coefficient linear system Y' == A(t) Y + b(t) of the commutative class "
        "A == a(t) I + b(t) K0 (K0 constant) via Phi == Exp[∫a] Exp[(∫b) K0].");
}
