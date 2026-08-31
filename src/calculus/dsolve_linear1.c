/*
 * dsolve_linear1.c — DSolve`LinearFirstOrder.
 *
 * Solves the first-order linear ODE  y'[x] + p(x) y[x] == q(x)  by the
 * integrating factor mu = Exp[Integrate[p, x]]:
 *
 *     y[x] = (Integrate[mu q, x] + C[1]) / mu.
 *
 * The equation is recognised by solving for y' (dsolve_solve_top_derivative),
 * replacing y[x] with a plain symbol Y, and confirming the RHS is affine in Y:
 * F(x, Y) == q(x) - p(x) Y with p, q free of Y.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

Expr** dsolve_linear1_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;

    const char* Yname = intern_symbol("DSolve`Y");
    Expr* yx = ds_make_funcapp(yname, 0, xvar);
    Expr* FY = ds_subst(F, yx, expr_new_symbol(Yname));   /* consumes F, yx */

    /* p = -D[FY, Y];  q = FY /. Y -> 0;  both must be free of Y */
    Expr* p = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                                     ds_d(expr_copy(FY), expr_new_symbol(Yname))));
    Expr* q = ds_subst(expr_copy(FY), expr_new_symbol(Yname), expr_new_integer(0));
    bool ok = !ds_contains(p, Yname) && !ds_contains(q, Yname);
    if (ok) {
        /* verify FY == q - p Y */
        Expr* pY = ds_call2(SYM_Times, expr_copy(p), expr_new_symbol(Yname));
        Expr* qmpY = eval_and_free(ds_call2(SYM_Subtract, expr_copy(q), pY));
        Expr* check = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY), qmpY));
        ok = ds_is_zero(check);
        expr_free(check);
    }
    expr_free(FY);
    if (!ok) { expr_free(p); expr_free(q); return NULL; }

    /* y' + p y == q  ->  integrating-factor solve (p, q consumed) */
    Expr* body = dsolve_linear_factor_solve(p, q, xvar);
    if (!body) return NULL;

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_linear1(Expr* res) {
    return dsolve_method_builtin(res, dsolve_linear1_try);
}

void dsolve_linear1_init(void) {
    symtab_add_builtin("DSolve`LinearFirstOrder", builtin_dsolve_linear1);
    symtab_get_def("DSolve`LinearFirstOrder")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LinearFirstOrder",
        "DSolve`LinearFirstOrder[eqn, y, x] solves y'[x] + p(x) y[x] == q(x) via "
        "the integrating factor Exp[Integrate[p, x]].");
}
