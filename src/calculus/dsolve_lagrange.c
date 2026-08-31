/*
 * dsolve_lagrange.c — DSolve`Lagrange (Lagrange / d'Alembert equation).
 *
 * Solves  y == x phi(y') + psi(y')  (with phi(p) != p, else it is Clairaut).
 * Writing p = y' and treating x as a function of p, differentiating w.r.t. x
 * gives the LINEAR first-order ODE
 *
 *     dx/dt - [phi'(t)/(t - phi(t))] x = psi'(t)/(t - phi(t)),   t = p,
 *
 * whose solution X(t, C) is obtained by the integrating-factor helper.  The
 * general solution is then PARAMETRIC:
 *
 *     x = X(t, C),   y = X(t, C) phi(t) + psi(t),
 *
 * with the slope t = y' as the parameter.  The branch is returned through the
 * wrapper DSolve`Param[X, Y, t]; the substrate (dsolve_run_parametric) verifies
 * it by substituting y' = D[Y,t]/D[X,t] into the residual and assembles
 * {{ x -> Function[{t}, X], y -> Function[{t}, Y] }}.  Recognition mirrors
 * dsolve_clairaut.c: work on the algebraic residual R(x, Y, p), require R linear
 * in Y, solve for Yexpr, then phi = d/dx Yexpr (free of x) and
 * psi = Yexpr - x phi (free of x and Y).
 *
 * Singular-line solutions (roots of phi(p) = p) and IVP constant-fitting are not
 * yet handled (future work); an IVP is declined by the substrate.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* base^-1; base consumed. */
static Expr* powneg1(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

/* A bare parameter symbol (Function args must be bare symbols), collision-safe:
 * the first of a small pool not equal to the independent variable and not
 * occurring in phi or psi. */
static const char* pick_param(const Expr* phi, const Expr* psi, const char* xvar) {
    static const char* pool[] = { "t", "s", "u", "w", "r", "q" };
    for (size_t i = 0; i < sizeof(pool) / sizeof(pool[0]); i++) {
        const char* c = intern_symbol(pool[i]);
        if (c == xvar) continue;
        if (ds_contains(phi, c) || ds_contains(psi, c)) continue;
        return c;
    }
    return intern_symbol("t");
}

Expr** dsolve_lagrange_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    if (P->ncond > 0) return NULL;              /* parametric IVP-fitting is future */
    const char* xvar = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`Y");
    const char* Pn = intern_symbol("DSolve`p");

    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;

    /* R linear in Y */
    Expr* dRdY = ds_d(expr_copy(R), expr_new_symbol(Yn));
    if (!ds_free_of(dRdY, Yn) || ds_is_zero(dRdY)) { expr_free(dRdY); expr_free(R); return NULL; }

    /* Yexpr = -(R|Y=0) / dRdY */
    Expr* R0 = ds_subst(expr_copy(R), expr_new_symbol(Yn), expr_new_integer(0));
    expr_free(R);
    Expr* Yexpr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
        expr_new_integer(-1), R0, powneg1(dRdY) }, 3));
    Yexpr = eval_and_free(ds_call1("Expand", Yexpr));

    /* phi(p) = d/dx Yexpr, must be free of x (Yexpr affine in x) */
    Expr* phi = ds_d(expr_copy(Yexpr), expr_new_symbol(xvar));
    if (!ds_free_of(phi, xvar)) { expr_free(phi); expr_free(Yexpr); return NULL; }

    /* psi(p) = Yexpr - x phi, free of x and Y */
    Expr* psi = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Yexpr),
                    ds_call2(SYM_Times, expr_new_symbol(xvar), expr_copy(phi))));
    expr_free(Yexpr);
    if (!ds_free_of(psi, xvar) || !ds_free_of(psi, Yn)) {
        expr_free(phi); expr_free(psi); return NULL;
    }

    /* not Clairaut: phi(p) != p (would also make t - phi zero) */
    Expr* pdiff = eval_and_free(ds_call2(SYM_Subtract, expr_copy(phi), expr_new_symbol(Pn)));
    bool is_clairaut = ds_is_zero(pdiff);
    expr_free(pdiff);
    if (is_clairaut) { expr_free(phi); expr_free(psi); return NULL; }

    /* not linear in y': a d'Alembert equation is genuinely nonlinear in y'.
     * phi constant (phi'==0) AND psi affine in p (psi''==0) makes the equation
     * y == (const) x + a y' + b, which is linear — owned by LinearFirstOrder.
     * Decline so the pinned method does not hand back a degenerate parametric
     * form for an equation with a clean explicit solution. */
    Expr* dphi  = ds_d(expr_copy(phi), expr_new_symbol(Pn));
    Expr* d2psi = ds_d(ds_d(expr_copy(psi), expr_new_symbol(Pn)), expr_new_symbol(Pn));
    bool is_linear = ds_is_zero(dphi) && ds_is_zero(d2psi);
    expr_free(dphi); expr_free(d2psi);
    if (is_linear) { expr_free(phi); expr_free(psi); return NULL; }

    /* choose the parameter symbol and re-express phi, psi in it (p -> t) */
    const char* tname = pick_param(phi, psi, xvar);
    Expr* phit = ds_subst(phi, expr_new_symbol(Pn), expr_new_symbol(tname));   /* consumes phi */
    Expr* psit = ds_subst(psi, expr_new_symbol(Pn), expr_new_symbol(tname));   /* consumes psi */

    /* linear ODE for x(t): dx/dt + Pcoef x = Qcoef,
     *   Pcoef = -phi'(t)/(t - phi),  Qcoef = psi'(t)/(t - phi) */
    Expr* tmphi = eval_and_free(ds_call2(SYM_Subtract, expr_new_symbol(tname), expr_copy(phit)));
    if (ds_is_zero(tmphi)) { expr_free(tmphi); expr_free(phit); expr_free(psit); return NULL; }
    Expr* phip = ds_d(expr_copy(phit), expr_new_symbol(tname));
    Expr* psip = ds_d(expr_copy(psit), expr_new_symbol(tname));
    Expr* Pcoef = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                     ds_call2(SYM_Times, phip, powneg1(expr_copy(tmphi)))));
    Expr* Qcoef = eval_and_free(ds_call2(SYM_Times, psip, powneg1(expr_copy(tmphi))));
    expr_free(tmphi);

    Expr* X = dsolve_linear_factor_solve(Pcoef, Qcoef, tname);   /* consumes Pcoef, Qcoef */
    if (!X) { expr_free(phit); expr_free(psit); return NULL; }

    /* Y = X phi + psi */
    Expr* Y = ds_simplify(eval_and_free(ds_call2(SYM_Plus,
                  ds_call2(SYM_Times, expr_copy(X), expr_copy(phit)),
                  expr_copy(psit))));
    expr_free(phit); expr_free(psit);

    Expr* wrap = expr_new_function(expr_new_symbol("DSolve`Param"),
                     (Expr*[]){ X, Y, expr_new_symbol(tname) }, 3);   /* consumes X, Y */
    Expr** out = malloc(sizeof(Expr*));
    out[0] = wrap;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_lagrange(Expr* res) {
    return dsolve_method_builtin_parametric(res, dsolve_lagrange_try);
}

void dsolve_lagrange_init(void) {
    symtab_add_builtin("DSolve`Lagrange", builtin_dsolve_lagrange);
    symtab_get_def("DSolve`Lagrange")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Lagrange",
        "DSolve`Lagrange[eqn, y, x] solves the Lagrange (d'Alembert) equation "
        "y == x f(y') + g(y') (with f(y') != y'). The general solution is "
        "parametric: {{x -> Function[{t}, X(t)], y -> Function[{t}, Y(t)]}} with "
        "parameter t = y', where X(t) solves the associated linear ODE.");
}
