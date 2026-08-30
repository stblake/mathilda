/*
 * dsolve_homogeneous.c — DSolve`Homogeneous.
 *
 * Solves the homogeneous first-order ODE  y'[x] == F(y/x)  (F of degree 0) by
 * the substitution y = v x.  Then y' = v + x v' = F, so x v' = F(v) - v, which
 * is separable:  Integrate[1/(F(v) - v), v] == Log[x] + C[1].  Solving for v and
 * multiplying by x recovers y.
 *
 * Degree-0 homogeneity is detected by substituting Y -> v x into F(x, Y) and
 * requiring the result to be free of x.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

Expr** dsolve_homogeneous_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;
    const char* Yn = intern_symbol("DSolve`Y");
    const char* vn = intern_symbol("DSolve`v");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    /* Fvx = FY /. Y -> v x ; must be free of x for degree-0 homogeneity */
    Expr* vx = ds_call2(SYM_Times, expr_new_symbol(vn), expr_new_symbol(xvar));
    Expr* Fvx = ds_subst(expr_copy(FY), expr_new_symbol(Yn), vx);
    expr_free(FY);
    if (!ds_free_of(Fvx, xvar)) { expr_free(Fvx); return NULL; }

    /* denom = F(v) - v */
    Expr* denom = eval_and_free(ds_call2(SYM_Subtract, Fvx, expr_new_symbol(vn)));
    if (ds_is_zero(denom)) { expr_free(denom); return NULL; }         /* y = C x, handled elsewhere */

    /* Integrate[1/(F(v)-v), v] == Log[x] + C[1] */
    Expr* intV = ds_integrate(
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ denom, expr_new_integer(-1) }, 2),
        expr_new_symbol(vn));
    if (ds_has_head(intV, SYM_Integrate)) { expr_free(intV); return NULL; }

    Expr* rhs = eval_and_free(ds_call2(SYM_Plus,
                    ds_call1("Log", expr_new_symbol(xvar)), ds_const(1)));
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ intV, rhs }, 2);
    Expr* solres = ds_solve(eq, expr_new_symbol(vn));
    size_t nv = 0;
    Expr** vs = dsolve_extract_solutions(solres, vn, &nv);
    if (solres) expr_free(solres);
    if (!vs) return NULL;

    /* y = x v */
    for (size_t i = 0; i < nv; i++)
        vs[i] = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(xvar), vs[i]));
    *nbranch = nv;
    return vs;
}

static Expr* builtin_dsolve_homogeneous(Expr* res) {
    return dsolve_method_builtin(res, dsolve_homogeneous_try);
}

void dsolve_homogeneous_init(void) {
    symtab_add_builtin("DSolve`Homogeneous", builtin_dsolve_homogeneous);
    symtab_get_def("DSolve`Homogeneous")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("DSolve`Homogeneous",
        "DSolve`Homogeneous[eqn, y, x] solves y'[x] == F(y/x) via the substitution "
        "y = v x, which separates variables in v and x.");
}
