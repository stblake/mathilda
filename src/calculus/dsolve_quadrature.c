/*
 * dsolve_quadrature.c — DSolve`Quadrature.
 *
 * Solves y^(n)[x] == f(x), where the right-hand side is free of the dependent
 * function: integrate f n times and add the general n-constant polynomial
 * C[1] + C[2] x + ... + C[n] x^(n-1).  This is the cleanest first stage and
 * also the base case for the linear/separable order-1 methods (which decline
 * here because their RHS depends on y).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

Expr** dsolve_quadrature_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];
    int n = P->max_order[0];
    if (n < 1) return NULL;

    Expr* F = dsolve_solve_top_derivative(P, n);
    if (!F) return NULL;
    /* pure quadrature requires the RHS to be free of the dependent function */
    if (ds_contains(F, yname)) { expr_free(F); return NULL; }

    /* integrate n times */
    Expr* body = F;
    for (int i = 0; i < n; i++) body = ds_integrate(body, expr_new_symbol(xvar));
    if (ds_has_head(body, SYM_Integrate)) { expr_free(body); return NULL; }

    /* + C[1] + C[2] x + ... + C[n] x^(n-1) */
    Expr** terms = malloc((size_t)(n + 1) * sizeof(Expr*));
    terms[0] = body;
    for (int k = 0; k < n; k++) {
        if (k == 0) terms[k + 1] = ds_const(1);
        else terms[k + 1] = ds_call2(SYM_Times, ds_const(k + 1),
                 ds_call2(SYM_Power, expr_new_symbol(xvar), expr_new_integer(k)));
    }
    Expr* full = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(n + 1)));
    free(terms);

    Expr** out = malloc(sizeof(Expr*));
    out[0] = full;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_quadrature(Expr* res) {
    return dsolve_method_builtin(res, dsolve_quadrature_try);
}

void dsolve_quadrature_init(void) {
    symtab_add_builtin("DSolve`Quadrature", builtin_dsolve_quadrature);
    symtab_get_def("DSolve`Quadrature")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Quadrature",
        "DSolve`Quadrature[eqn, y, x] solves y^(n)[x] == f(x) with f free of y "
        "by integrating n times and adding the general constant polynomial.");
}
