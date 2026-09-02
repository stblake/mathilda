/*
 * dsolve_sepreduced.c — DSolve`SeparableReduced.
 *
 * Solves y' == F(x, y) of the form  x y'/y == G(x^n y)  — i.e. r := x F/y is a
 * function of the single combination w = x^n y — by the substitution w = x^n y.
 * Then  w' = (w/x)(n + G(w)),  which is SEPARABLE:
 *     Integrate[1/(w (n + G(w))), w] == Log[x] + C[1].
 * The power n is recovered as n = x r_x / (y r_y) (constant iff the form holds).
 * The first integral G(x, y) == C[1] with w -> x^n y is returned through
 * dsolve_run_implicit (implicit differentiation verifies it).  Mirrors SymPy's
 * `separable_reduced`.  Runs among the first-order substitution reductions.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* x^p (p an expression), owned; base/p consumed via copies by caller convention. */
static Expr* xpow_expr(const char* xvar, Expr* p) {
    return expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_new_symbol(xvar), p }, 2);
}

Expr** dsolve_sepreduced_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`srY");
    const char* wn = intern_symbol("DSolve`srW");

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;

    /* r = x F / y, with y[x] -> Y */
    Expr* r = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(xvar),
                  ds_call2(SYM_Times, F,
                      expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ ds_make_funcapp(yname, 0, xvar), expr_new_integer(-1) }, 2))));
    Expr* rY = ds_simplify(ds_subst(r, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn)));

    /* n = x r_x / (Y r_Y); require r genuinely depends on Y and n is constant */
    Expr* rx = ds_d(expr_copy(rY), expr_new_symbol(xvar));
    Expr* ry = ds_d(expr_copy(rY), expr_new_symbol(Yn));
    if (ds_is_zero(ry)) { expr_free(rY); expr_free(rx); expr_free(ry); return NULL; }
    Expr* n = ds_simplify(eval_and_free(ds_call2(SYM_Times,
                  ds_call2(SYM_Times, expr_new_symbol(xvar), rx),
                  expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ ds_call2(SYM_Times, expr_new_symbol(Yn), ry),
                                 expr_new_integer(-1) }, 2))));
    if (!ds_free_of(n, xvar) || !ds_free_of(n, Yn)) { expr_free(rY); expr_free(n); return NULL; }

    /* G(w) = rY with Y -> w x^{-n} (must be free of x) */
    Expr* wxn = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(wn),
                    xpow_expr(xvar, eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(n))))));
    Expr* Gw = ds_simplify(ds_subst(rY, expr_new_symbol(Yn), wxn));   /* consumes rY */
    if (!ds_free_of(Gw, xvar)) { expr_free(n); expr_free(Gw); return NULL; }

    /* Integrate[1/(w (n + G(w))), w] */
    Expr* denom = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(wn),
                      ds_call2(SYM_Plus, expr_copy(n), Gw)));
    if (ds_is_zero(denom)) { expr_free(n); expr_free(denom); return NULL; }
    Expr* intW = ds_integrate(
                     expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ denom, expr_new_integer(-1) }, 2),
                     expr_new_symbol(wn));
    if (ds_has_head(intW, SYM_Integrate)) { expr_free(n); expr_free(intW); return NULL; }

    /* G_expr = (intW /. w -> x^n y[x]) - Log[x] ;  solution is G_expr == C[1] */
    Expr* wsub = eval_and_free(ds_call2(SYM_Times, xpow_expr(xvar, n),   /* consumes n */
                     ds_make_funcapp(yname, 0, xvar)));
    Expr* G = ds_subst(intW, expr_new_symbol(wn), wsub);
    G = eval_and_free(ds_call2(SYM_Subtract, G, ds_call1("Log", expr_new_symbol(xvar))));

    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_sepreduced(Expr* res) {
    return dsolve_method_builtin_implicit(res, dsolve_sepreduced_try);
}

void dsolve_sepreduced_init(void) {
    symtab_add_builtin("DSolve`SeparableReduced", builtin_dsolve_sepreduced);
    symtab_get_def("DSolve`SeparableReduced")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SeparableReduced",
        "DSolve`SeparableReduced[eqn, y, x] solves y' == F where x F/y is a function "
        "of w = x^n y, via the substitution w = x^n y giving a separable equation; "
        "returns the implicit first integral.");
}
