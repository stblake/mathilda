/*
 * dsolve_almostlinear.c — DSolve`AlmostLinear.
 *
 * An equation A(x,y) y' + B(x,y) == 0 is "almost linear" when A separates as
 * f(x) g(y) and the substitution u = l(y) = Integrate[g, y] linearises it:
 * u' = g y' so f u' = -B, i.e. u' = H(x,y) with H = -B/f, and the equation is
 * linear in u exactly when H = Q(x) - P(x) u, i.e.
 *     P(x) = -H_y / g   and   Q(x) = H + P l(y)
 * are both free of y.  Then u' + P u == Q is a LinearFirstOrder equation
 * (integrating factor), and l(y[x]) == U(x) is solved for y.  Mirrors SymPy's
 * `almost_linear`.  Runs among the first-order substitution reductions.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* e /. sym -> val (sym a symbol name, val consumed); e consumed. */
static Expr* subst_sym(Expr* e, const char* sym, Expr* val) {
    return ds_subst(e, expr_new_symbol(sym), val);
}

Expr** dsolve_almostlinear_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`alY");
    const char* pn = intern_symbol("DSolve`alP");
    Expr* R = P->eq_residuals[0];

    /* residual with y'[x] -> plain symbol p and y[x] -> plain symbol Y, so the
     * coefficient of y' is read with Coefficient (NOT D w.r.t. the funcapp y'[x],
     * whose derivative path leaks a node in compute_deriv). */
    Expr* Rp = ds_subst(expr_copy(R), ds_make_funcapp(yname, 1, xvar), expr_new_symbol(pn));
    Rp = ds_subst(Rp, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    /* AY = coeff of p ; BY = Rp with p -> 0 ; require R linear in y' (Rp == AY p + BY) */
    Expr* AY = eval_and_free(ds_call2("Coefficient", expr_copy(Rp), expr_new_symbol(pn)));
    Expr* BY = ds_subst(expr_copy(Rp), expr_new_symbol(pn), expr_new_integer(0));
    Expr* recon = eval_and_free(ds_call2(SYM_Subtract, Rp,
                      ds_call2(SYM_Plus,
                          ds_call2(SYM_Times, expr_copy(AY), expr_new_symbol(pn)),
                          expr_copy(BY))));
    bool linear = ds_is_zero(recon);
    expr_free(recon);
    if (!linear || ds_is_zero(AY)) { expr_free(AY); expr_free(BY); return NULL; }

    /* separate AY = f(x) g(Y):  f = AY|_{Y=1},  g = AY/f (must be free of x) */
    Expr* f = subst_sym(expr_copy(AY), Yn, expr_new_integer(1));
    if (ds_is_zero(f)) { expr_free(AY); expr_free(BY); expr_free(f); return NULL; }
    Expr* g = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_copy(AY),
                  expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(f), expr_new_integer(-1) }, 2))));
    expr_free(AY);
    if (!ds_free_of(g, xvar)) { expr_free(BY); expr_free(f); expr_free(g); return NULL; }

    /* l(Y) = Integrate[g, Y] */
    Expr* l = ds_integrate(expr_copy(g), expr_new_symbol(Yn));
    if (ds_has_head(l, SYM_Integrate)) { expr_free(BY); expr_free(f); expr_free(g); expr_free(l); return NULL; }

    /* H = -B/f ;  P = -H_Y/g ;  Q = H + P l  (P, Q must be free of Y) */
    Expr* finv = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy(f), expr_new_integer(-1) }, 2);
    Expr* H = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                  ds_call2(SYM_Times, expr_copy(BY), finv)));
    expr_free(BY); expr_free(f);

    Expr* Hy = ds_d(expr_copy(H), expr_new_symbol(Yn));
    Expr* Pc = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                   ds_call2(SYM_Times, Hy,
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ expr_copy(g), expr_new_integer(-1) }, 2)))));
    expr_free(g);
    if (!ds_free_of(Pc, Yn)) { expr_free(H); expr_free(l); expr_free(Pc); return NULL; }

    Expr* Qc = ds_simplify(eval_and_free(ds_call2(SYM_Plus, expr_copy(H),
                   ds_call2(SYM_Times, expr_copy(Pc), expr_copy(l)))));
    expr_free(H);
    if (!ds_free_of(Qc, Yn)) { expr_free(l); expr_free(Pc); expr_free(Qc); return NULL; }

    /* clear the (now inert) Y so the coefficients are pure functions of x */
    Pc = subst_sym(Pc, Yn, expr_new_integer(1));
    Qc = subst_sym(Qc, Yn, expr_new_integer(1));

    /* solve u' + Pc u == Qc (integrating factor) */
    Expr* U = dsolve_linear_factor_solve(Pc, Qc, xvar);   /* consumes Pc, Qc */
    if (!U) { expr_free(l); return NULL; }

    /* l(Y) == U -> Solve for Y (explicit y(x) bodies) */
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ l, U }, 2);
    Expr* sol = ds_solve(eq, expr_new_symbol(Yn));
    size_t nb = 0;
    Expr** bodies = dsolve_extract_solutions(sol, Yn, &nb);
    if (sol) expr_free(sol);
    if (!bodies) return NULL;

    *nbranch = nb;
    return bodies;
}

static Expr* builtin_dsolve_almostlinear(Expr* res) {
    return dsolve_method_builtin(res, dsolve_almostlinear_try);
}

void dsolve_almostlinear_init(void) {
    symtab_add_builtin("DSolve`AlmostLinear", builtin_dsolve_almostlinear);
    symtab_get_def("DSolve`AlmostLinear")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`AlmostLinear",
        "DSolve`AlmostLinear[eqn, y, x] solves f(x)g(y) y' + k(x)l(y) + m(x) == 0 by the "
        "substitution u = Integrate[g, y], which turns it into a linear first-order ODE "
        "u' + P(x) u == Q(x), then solves l(y) == U(x) for y.");
}
