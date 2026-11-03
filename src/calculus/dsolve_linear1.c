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
#include "../core.h"   /* tc_register_unwind_counter */
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

    /* Distributional / step forcing on an IVP.  The plain integrating-factor form
     * below integrates mu*q INDEFINITELY, which stays inert for a DiracDelta / step
     * q; the downstream IC-fit then substitutes x -> x0 into that inert Integrate
     * (a capture on its bound variable), yielding a malformed Integrate[.,x0] and a
     * solution WRONG past the impulse / breakpoint (e.g. y'+2y==2 UnitStep[t-1],
     * y(0)=1 was correct only for t<1).  Route an IVP two ways:
     *   - UnitStep / Piecewise / Heaviside (no impulse): decline so the proven,
     *     per-interval-verified DSolve`PiecewiseForcing (interval continuation)
     *     owns it — it is already order-agnostic.
     *   - DiracDelta impulse (no step): use the CAUSAL definite convolution, whose
     *     definite DiracDelta integral sifts to a clean HeavisideTheta.
     * A mixed step+impulse forcing declines (neither path closes it) rather than
     * ship the broken inert form.  A no-IC general-solution request keeps the
     * indefinite form unchanged (inert but not wrong — zero regression). */
    bool has_dirac = ds_contains(q, intern_symbol("DiracDelta"));
    bool has_step  = ds_contains(q, SYM_UnitStep) || ds_contains(q, SYM_Piecewise)
                     || ds_contains(q, SYM_HeavisideTheta);
    Expr* body;
    if (P->ncond >= 1 && (has_dirac || has_step)) {
        if (has_step && !has_dirac) {
            /* Pure step forcing.  With NUMERIC breakpoints, decline so the proven,
             * per-interval-verified DSolve`PiecewiseForcing (interval continuation)
             * returns a clean Piecewise.  A SYMBOLIC breakpoint — UnitStep[t-a] — it
             * cannot resolve (pw_numval is NaN), so route THAT to the causal form.
             * Probe: substitute a non-special real for x; a numeric-breakpoint
             * forcing then collapses to a pure number (NumericQ True), whereas a
             * symbolic breakpoint keeps the parameter inside the step (False). */
            Expr* nq = eval_and_free(ds_call1("NumericQ",
                           ds_subst(expr_copy(q), expr_new_symbol(xvar),
                                    expr_new_real(0.3728810))));
            bool numeric_bp = (nq && nq->type == EXPR_SYMBOL && nq->data.symbol.name == SYM_True);
            expr_free(nq);
            if (numeric_bp) { expr_free(p); expr_free(q); return NULL; }  /* -> PiecewiseForcing */
            /* Symbolic breakpoint.  When Integrate closes mu*q (a step with a
             * linear argument integrates to the CONTINUOUS antiderivative
             * (G(x) - G(a)) UnitStep[x - a]), the plain integrating-factor form
             * is a closed form valid across the breakpoint, and the IC-fit
             * substitutes into no inert integral -- the capture hazard above
             * cannot arise.  Only a non-closing integral needs the causal
             * convolution, whose bound dummy must otherwise stay inside it. */
            body = dsolve_linear_factor_solve(expr_copy(p), expr_copy(q), xvar);
            if (body && !ds_has_active_integrate(body)) {
                expr_free(p); expr_free(q);
            } else {
                if (body) expr_free(body);
                body = dsolve_linear_factor_solve_causal(p, q, xvar);
            }
        } else if (has_dirac && !has_step) {
            /* Impulse: Integrate sifts g(x) DiracDelta[x - a] indefinitely to
             * g(a) HeavisideTheta[x - a], a closed form valid across a; only a
             * non-closing integral needs the causal convolution. */
            body = dsolve_linear_factor_solve(expr_copy(p), expr_copy(q), xvar);
            if (body && !ds_has_active_integrate(body)) {
                expr_free(p); expr_free(q);
            } else {
                if (body) expr_free(body);
                body = dsolve_linear_factor_solve_causal(p, q, xvar);
            }
        } else {
            /* Mixed step + impulse: neither path closes it cleanly — decline
             * (honest UNEVAL) rather than ship the broken inert indefinite form. */
            expr_free(p); expr_free(q); return NULL;
        }
    } else {
        /* y' + p y == q  ->  integrating-factor solve (p, q consumed).  A general
         * solution keeps a non-elementary integrating factor as an inert integral
         * (y' == f[x] y + g[x]); an IVP needs a closed form to fit at x0. */
        body = (P->ncond == 0) ? dsolve_linear_factor_solve_inert(p, q, xvar)
                               : dsolve_linear_factor_solve(p, q, xvar);
    }
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
    tc_register_unwind_counter(&g_linfac_inert_active);   /* restored on a timeout unwind */
    symtab_add_builtin("DSolve`LinearFirstOrder", builtin_dsolve_linear1);
    symtab_get_def("DSolve`LinearFirstOrder")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LinearFirstOrder",
        "DSolve`LinearFirstOrder[eqn, y, x] solves y'[x] + p(x) y[x] == q(x) via "
        "the integrating factor Exp[Integrate[p, x]]; for a general solution a "
        "non-elementary Integrate[p, x] (e.g. y' == f[x] y + g[x]) is kept as "
        "Inactive[Integrate].");
}
