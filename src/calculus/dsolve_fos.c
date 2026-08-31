/*
 * dsolve_fos.c — DSolve`FirstOrderSubstitution.
 *
 * Solves the first-order ODE  y'[x] == F(a x + b y + c)  whose right-hand side
 * depends on x and y only through a single linear combination.  Writing the RHS
 * as F(x, Y) (Y a plain symbol standing for y[x]), the level curves of F are the
 * straight lines along which the combination is constant, so F is a function of
 * that combination iff the ratio
 *     r = F_x / F_Y
 * is constant (free of x and Y).  With the substitution v = y + r x we have
 *     v' = y' + r = F + r = H(v) + r,   H(W) := F(x, W - r x)  (free of x),
 * an autonomous first-order equation solved by separation:
 *     Integrate[1/(r + H(W)), W] == x + C[1],
 * solved for W = v, and finally y = v - r x.
 *
 * This is the solved-for-y' companion to Homogeneous (a different substitution):
 * it fires exactly when the RHS mixes x and y in a way no earlier first-order
 * method (linear / separable / homogeneous / exact) recognises, e.g.
 * y' == (x + y)^2, y' == (2 y - x)^2, y' == (x + y + 1)^2.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

Expr** dsolve_fos_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* xvar = P->ind_names[0];
    const char* yname = P->fun_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;

    const char* Yname = intern_symbol("DSolve`fosY");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yname)); /* consumes F */

    /* r = F_x / F_Y, required constant (free of x and Y).  A zero derivative
     * means F is free of x (separable) or of y (linear/quadrature) — declined so
     * the simpler method owns it. */
    Expr* Fx  = ds_d(expr_copy(FY), expr_new_symbol(xvar));
    Expr* FYd = ds_d(expr_copy(FY), expr_new_symbol(Yname));
    if (ds_is_zero(Fx) || ds_is_zero(FYd)) {
        expr_free(Fx); expr_free(FYd); expr_free(FY); return NULL;
    }
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
        Fx,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ FYd, expr_new_integer(-1) }, 2)
    }, 2));
    if (!ds_free_of(r, xvar) || !ds_free_of(r, Yname)) { expr_free(r); expr_free(FY); return NULL; }
    /* Collapse the ratio to its closed constant form: the substitution below
     * replaces Y, so r must be *textually* free of Y (a proved-constant ratio
     * that Simplify cannot reduce is declined rather than corrupted). */
    r = ds_simplify(r);
    if (ds_contains(r, Yname) || ds_contains(r, xvar)) { expr_free(r); expr_free(FY); return NULL; }

    /* H(W) = F(x, W - r x); a genuine linear combination iff free of x. */
    const char* Wname = intern_symbol("DSolve`fosW");
    Expr* Wminus = eval_and_free(ds_call2(SYM_Subtract, expr_new_symbol(Wname),
                       ds_call2(SYM_Times, expr_copy(r), expr_new_symbol(xvar))));
    Expr* H = ds_subst(FY, expr_new_symbol(Yname), Wminus); /* consumes FY, Wminus */
    if (!ds_free_of(H, xvar)) { expr_free(H); expr_free(r); return NULL; }

    /* autonomous v' == r + H(v): separate  Integrate[1/(r + H), W] == x + C[1]. */
    Expr* denom = eval_and_free(ds_call2(SYM_Plus, expr_copy(r), H)); /* consumes H */
    if (ds_is_zero(denom)) { expr_free(denom); expr_free(r); return NULL; }
    Expr* invD = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ denom, expr_new_integer(-1) }, 2));
    Expr* lhsInt = ds_integrate(invD, expr_new_symbol(Wname));
    if (ds_has_head(lhsInt, SYM_Integrate)) { expr_free(lhsInt); expr_free(r); return NULL; }

    Expr* rhs = eval_and_free(ds_call2(SYM_Plus, expr_new_symbol(xvar), ds_const(1)));
    Expr* eq  = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ lhsInt, rhs }, 2);
    Expr* solres = ds_solve(eq, expr_new_symbol(Wname));
    size_t nw = 0;
    Expr** wsol = dsolve_extract_solutions(solres, Wname, &nw);
    if (solres) expr_free(solres);
    if (!wsol) { expr_free(r); return NULL; }

    /* y = v - r x for each solved branch. */
    Expr** out = malloc(nw * sizeof(Expr*));
    for (size_t i = 0; i < nw; i++)
        out[i] = eval_and_free(ds_call2(SYM_Subtract, wsol[i],
                     ds_call2(SYM_Times, expr_copy(r), expr_new_symbol(xvar))));
    free(wsol);
    expr_free(r);
    *nbranch = nw;
    return out;
}

static Expr* builtin_dsolve_fos(Expr* res) {
    return dsolve_method_builtin(res, dsolve_fos_try);
}

void dsolve_fos_init(void) {
    symtab_add_builtin("DSolve`FirstOrderSubstitution", builtin_dsolve_fos);
    symtab_get_def("DSolve`FirstOrderSubstitution")->attributes
        |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`FirstOrderSubstitution",
        "DSolve`FirstOrderSubstitution[eqn, y, x] solves y'[x] == F(a x + b y + c) "
        "by the substitution v = y + (F_x/F_y) x, reducing it to an autonomous "
        "separable equation in v.");
}
