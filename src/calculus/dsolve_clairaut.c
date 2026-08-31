/*
 * dsolve_clairaut.c — DSolve`Clairaut.
 *
 * Solves the Clairaut equation  y == x y' + f(y')  whose general solution is the
 * one-parameter family of lines  y = C[1] x + f(C[1]).  With IncludeSingularSolutions
 * the singular envelope is added, obtained by eliminating p from
 * { y = x p + f(p), 0 = x + f'(p) }.
 *
 * Detection works on the algebraic residual R(x, Y, p) (p = y'): R must be linear
 * in Y, and with Yexpr the solution of R == 0 for Y we require d/dx Yexpr == p and
 * f(p) := Yexpr - x p free of x and Y.  Clairaut equations are nonlinear in y', so
 * the linear-in-y' methods decline before this one is reached.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

Expr** dsolve_clairaut_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* xvar = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`Y");
    const char* Pn = intern_symbol("DSolve`p");

    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;

    /* R linear in Y: dRdY free of Y and non-zero */
    Expr* dRdY = ds_d(expr_copy(R), expr_new_symbol(Yn));
    if (!ds_free_of(dRdY, Yn) || ds_is_zero(dRdY)) { expr_free(dRdY); expr_free(R); return NULL; }

    /* Yexpr = -(R|Y=0) / dRdY */
    Expr* R0 = ds_subst(expr_copy(R), expr_new_symbol(Yn), expr_new_integer(0));
    expr_free(R);
    Expr* Yexpr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
        expr_new_integer(-1), R0,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ dRdY, expr_new_integer(-1) }, 2)
    }, 3));

    /* require d/dx Yexpr == p */
    Expr* chk = eval_and_free(ds_call2(SYM_Subtract,
                    ds_d(expr_copy(Yexpr), expr_new_symbol(xvar)), expr_new_symbol(Pn)));
    bool ok = ds_is_zero(chk);
    expr_free(chk);
    if (!ok) { expr_free(Yexpr); return NULL; }
    Yexpr = eval_and_free(ds_call1("Expand", Yexpr));   /* distribute a leading -1 */

    /* f(p) = Yexpr - x p, must be free of x and Y */
    Expr* fp = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Yexpr),
                    ds_call2(SYM_Times, expr_new_symbol(xvar), expr_new_symbol(Pn))));
    if (!ds_free_of(fp, xvar) || !ds_free_of(fp, Yn)) { expr_free(fp); expr_free(Yexpr); return NULL; }

    /* general solution: y = Yexpr /. p -> C[1] */
    size_t nb = 0, cap = 4;
    Expr** bodies = malloc(cap * sizeof(Expr*));
    bodies[nb++] = ds_subst(expr_copy(Yexpr), expr_new_symbol(Pn), ds_const(1));

    /* singular envelope: solve x + f'(p) == 0 for p, substitute into Yexpr */
    if (P->include_singular) {
        Expr* seq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){
            eval_and_free(ds_call2(SYM_Plus, expr_new_symbol(xvar), ds_d(expr_copy(fp), expr_new_symbol(Pn)))),
            expr_new_integer(0)
        }, 2);
        Expr* ssol = ds_solve(seq, expr_new_symbol(Pn));
        size_t np = 0;
        Expr** ps = dsolve_extract_solutions(ssol, Pn, &np);
        if (ssol) expr_free(ssol);
        for (size_t i = 0; i < np; i++) {
            if (nb >= cap) { cap *= 2; bodies = realloc(bodies, cap * sizeof(Expr*)); }
            bodies[nb++] = ds_subst(expr_copy(Yexpr), expr_new_symbol(Pn), ps[i]);
        }
        free(ps);
    }

    expr_free(fp); expr_free(Yexpr);
    *nbranch = nb;
    return bodies;
}

static Expr* builtin_dsolve_clairaut(Expr* res) {
    return dsolve_method_builtin(res, dsolve_clairaut_try);
}

void dsolve_clairaut_init(void) {
    symtab_add_builtin("DSolve`Clairaut", builtin_dsolve_clairaut);
    symtab_get_def("DSolve`Clairaut")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Clairaut",
        "DSolve`Clairaut[eqn, y, x] solves y == x y' + f(y'); the general solution "
        "is y = C[1] x + f(C[1]).  With IncludeSingularSolutions -> True the singular "
        "envelope is also returned.");
}
