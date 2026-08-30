/*
 * dsolve_separable.c — DSolve`Separable.
 *
 * Solves the first-order separable ODE  y'[x] == g(x) h(y).  The RHS F(x, Y)
 * (Y a plain symbol standing for y[x]) is separated by evaluating it at sample
 * points: for X0, Y0 with F(X0,Y0) != 0,
 *
 *     g(x) = F(x, Y0),   h(Y) = F(X0, Y) / F(X0, Y0),
 *
 * which reconstructs F exactly when F is genuinely separable (verified by
 * PossibleZeroQ[F - g h]).  Then Integrate[1/h, Y] == Integrate[g, x] + C[1] is
 * solved for Y.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

Expr** dsolve_separable_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* xvar = P->ind_names[0];
    const char* yname = P->fun_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;

    const char* Yname = intern_symbol("DSolve`Y");
    Expr* yx = ds_make_funcapp(yname, 0, xvar);
    Expr* FY = ds_subst(F, yx, expr_new_symbol(Yname));   /* consumes F, yx */

    /* Search for a separation via sample points. */
    static const int samples[] = { 2, 3, 5, 7, -2, -3 };
    size_t nsamp = sizeof(samples) / sizeof(samples[0]);
    Expr* g = NULL; Expr* h = NULL;
    for (size_t ai = 0; ai < nsamp && !g; ai++) {
        for (size_t bi = 0; bi < nsamp && !g; bi++) {
            int X0 = samples[ai], Y0 = samples[bi];
            Expr* denom = ds_subst(
                ds_subst(expr_copy(FY), expr_new_symbol(xvar), expr_new_integer(X0)),
                expr_new_symbol(Yname), expr_new_integer(Y0));
            if (!ds_is_nonzero(denom)) { expr_free(denom); continue; }
            Expr* gcand = ds_subst(expr_copy(FY), expr_new_symbol(Yname), expr_new_integer(Y0));
            Expr* hx0 = ds_subst(expr_copy(FY), expr_new_symbol(xvar), expr_new_integer(X0));
            Expr* hcand = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                hx0,
                expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(denom), expr_new_integer(-1) }, 2)
            }, 2));
            expr_free(denom);
            /* verify F == g h */
            Expr* prod = ds_call2(SYM_Times, expr_copy(gcand), expr_copy(hcand));
            Expr* check = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY), prod));
            if (ds_is_zero(check)) { g = gcand; h = hcand; }
            else { expr_free(gcand); expr_free(hcand); }
            expr_free(check);
        }
    }
    expr_free(FY);
    if (!g) return NULL;

    /* Integrate[1/h, Y] == Integrate[g, x] + C[1] */
    Expr* invh = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ h, expr_new_integer(-1) }, 2));
    Expr* lhsInt = ds_integrate(invh, expr_new_symbol(Yname));
    Expr* rhsInt = ds_integrate(g, expr_new_symbol(xvar));
    if (ds_has_head(lhsInt, SYM_Integrate) || ds_has_head(rhsInt, SYM_Integrate)) {
        expr_free(lhsInt); expr_free(rhsInt); return NULL;
    }
    Expr* rhs = eval_and_free(ds_call2(SYM_Plus, rhsInt, ds_const(1)));
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ lhsInt, rhs }, 2);

    Expr* solres = ds_solve(eq, expr_new_symbol(Yname));
    size_t nb = 0;
    Expr** bodies = dsolve_extract_solutions(solres, Yname, &nb);
    if (solres) expr_free(solres);
    if (!bodies) return NULL;
    *nbranch = nb;
    return bodies;
}

static Expr* builtin_dsolve_separable(Expr* res) {
    return dsolve_method_builtin(res, dsolve_separable_try);
}

void dsolve_separable_init(void) {
    symtab_add_builtin("DSolve`Separable", builtin_dsolve_separable);
    symtab_get_def("DSolve`Separable")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("DSolve`Separable",
        "DSolve`Separable[eqn, y, x] solves y'[x] == g(x) h(y) by separating "
        "variables: Integrate[1/h, y] == Integrate[g, x] + C[1], solved for y.");
}
