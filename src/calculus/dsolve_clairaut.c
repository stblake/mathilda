/*
 * dsolve_clairaut.c — DSolve`Clairaut.
 *
 * Solves the Clairaut equation  y == x y' + f(y')  whose general solution is the
 * one-parameter family of lines  y = C[1] x + f(C[1]).  With IncludeSingularSolutions
 * the singular envelope is added, obtained by eliminating p from
 * { y = x p + f(p), 0 = x + f'(p) }.
 *
 * Detection works on the algebraic residual R(x, Y, p) (p = y'): solve R == 0 for Y
 * and, for each root Yexpr, require d/dx Yexpr == p and f(p) := Yexpr - x p free of
 * x and Y.  Clairaut equations are nonlinear in y', so the linear-in-y' methods
 * decline before this one is reached.
 *
 * R linear in Y is the common spelling and is solved directly.  R merely POLYNOMIAL in
 * Y is also Clairaut — `(y - x y')^2 == 1 + y'^2` is the textbook example, whose roots
 * y = x p ± Sqrt[1 + p^2] are two Clairaut equations (§2.2.34-3331).  That family used
 * to be invisible: this method declined on the nonlinearity, and `SolvableForY`, which
 * does isolate the roots, throws them away on purpose, because for a Clairaut equation
 * the differentiation method's denominator p - dG/dx vanishes identically -- its comment
 * says the family is "owned earlier", which it was not.  So each root is now run through
 * the same per-root test.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* Test one candidate root Yexpr of R == 0 for Y, and on success append the line family
 * (and, under IncludeSingularSolutions, the envelope) to `bodies`.  `Yexpr` consumed. */
static void clairaut_emit(DSolveProblem* P, Expr* Yexpr, const char* xvar,
                          const char* Yn, const char* Pn,
                          Expr*** bodies, size_t* nb, size_t* cap) {
    /* require d/dx Yexpr == p */
    Expr* chk = eval_and_free(ds_call2(SYM_Subtract,
                    ds_d(expr_copy(Yexpr), expr_new_symbol(xvar)), expr_new_symbol(Pn)));
    bool ok = ds_is_zero(chk);
    expr_free(chk);
    if (!ok) { expr_free(Yexpr); return; }
    Yexpr = eval_and_free(ds_call1("Expand", Yexpr));   /* distribute a leading -1 */

    /* f(p) = Yexpr - x p, must be free of x and Y */
    Expr* fp = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Yexpr),
                    ds_call2(SYM_Times, expr_new_symbol(xvar), expr_new_symbol(Pn))));
    if (!ds_free_of(fp, xvar) || !ds_free_of(fp, Yn)) { expr_free(fp); expr_free(Yexpr); return; }

    /* general solution: y = Yexpr /. p -> C[1] */
    if (*nb >= *cap) { *cap *= 2; *bodies = realloc(*bodies, *cap * sizeof(Expr*)); }
    (*bodies)[(*nb)++] = ds_subst(expr_copy(Yexpr), expr_new_symbol(Pn), ds_const(1));

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
            if (*nb >= *cap) { *cap *= 2; *bodies = realloc(*bodies, *cap * sizeof(Expr*)); }
            (*bodies)[(*nb)++] = ds_subst(expr_copy(Yexpr), expr_new_symbol(Pn), ps[i]);
        }
        free(ps);
    }
    expr_free(fp); expr_free(Yexpr);
}

Expr** dsolve_clairaut_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* xvar = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`Y");
    const char* Pn = intern_symbol("DSolve`p");

    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;

    size_t nb = 0, cap = 4;
    Expr** bodies = malloc(cap * sizeof(Expr*));

    Expr* dRdY = ds_d(expr_copy(R), expr_new_symbol(Yn));
    if (ds_is_zero(dRdY)) { expr_free(dRdY); expr_free(R); free(bodies); return NULL; }
    if (ds_free_of(dRdY, Yn)) {
        /* R linear in Y: Yexpr = -(R|Y=0) / dRdY, solved directly (the common spelling;
         * this path is byte-identical to what it always was). */
        Expr* R0 = ds_subst(expr_copy(R), expr_new_symbol(Yn), expr_new_integer(0));
        Expr* Yexpr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
            expr_new_integer(-1), R0,
            expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ dRdY, expr_new_integer(-1) }, 2)
        }, 3));
        clairaut_emit(P, Yexpr, xvar, Yn, Pn, &bodies, &nb, &cap);
    } else {
        /* R polynomial in Y: solve for Y and test each root.  Gated on PolynomialQ so a
         * transcendental residual never reaches Solve here (the inverse-function branches
         * would be spurious and the solve slow). */
        expr_free(dRdY);
        Expr* pq = eval_and_free(ds_call2(intern_symbol("PolynomialQ"), expr_copy(R),
                                          expr_new_symbol(Yn)));
        bool poly = (pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True);
        expr_free(pq);
        if (poly) {
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                           (Expr*[]){ expr_copy(R), expr_new_integer(0) }, 2);
            Expr* sol = ds_solve(eq, expr_new_symbol(Yn));
            size_t nr = 0;
            Expr** roots = dsolve_extract_solutions(sol, Yn, &nr);
            if (sol) expr_free(sol);
            for (size_t i = 0; i < nr; i++)
                clairaut_emit(P, roots[i], xvar, Yn, Pn, &bodies, &nb, &cap);
            free(roots);
        }
    }
    expr_free(R);

    if (nb == 0) { free(bodies); return NULL; }
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
