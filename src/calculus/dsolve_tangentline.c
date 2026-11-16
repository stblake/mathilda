/*
 * dsolve_tangentline.c -- DSolve`TangentLine: the first-order IVP whose
 * solution curve is a straight line.
 *
 * For y' == F(x, y) with one point condition y(x0) == y0, the solution curve
 * has slope m = F(x0, y0) at the initial point.  If the tangent line
 *
 *     L(x) = y0 + m (x - x0)
 *
 * satisfies the ODE IDENTICALLY in x (F(x, L(x)) - m == 0, decided by the zero
 * test), it is a solution through the point, and when F_y(x0, y0) is finite F
 * is locally Lipschitz in y there, so by Picard-Lindelof it is THE solution.
 * This is a decision for "is the solution of this IVP a line", not a guess: the
 * candidate is unique, and both the ODE and the uniqueness condition are
 * checked.
 *
 * Why the cascade needs it.  The line is typically an invariant curve that the
 * general-solution methods DIVIDE OUT, so their one-parameter family cannot
 * reach the initial point -- the fitted constant is infinite, and that branch
 * is (correctly) dropped:
 *
 *   y' == x Log[y],      y(1) == 1   separable; Integrate[1/Log[y], y] =
 *                                    LogIntegral[y] is infinite at y = 1.
 *                                    Answer y == 1 (the equilibrium Log[y] == 0).
 *   y' == Sqrt[x - y],   y(2) == 1   u = x - y gives u' == 1 - Sqrt[u], whose
 *                                    first integral has Log[1 - Sqrt[u]];
 *                                    u(2) == 1 is its equilibrium.  Answer
 *                                    y == x - 1, which is not constant in y.
 *   y' == Log[1 + y^2],  y(0) == 0   the y-integral is non-elementary and
 *                                    divergent at 0.  Answer y == 0.
 *
 * Every equilibrium y == y0 of an equation with F(x, y0) == 0 is the m == 0
 * case.  The same candidate serves the substrate's singular-IVP fallback
 * (dsolve_common.c), where a method's family provably cannot meet the
 * condition: y' == 1 + (y - x)^2 with y(0) == 0 has the Riccati family
 * x + 1/(-x - C[1]), which no finite C[1] fits, and the line y == x.  The stage runs LAST in the Automatic cascade, so it can only turn an
 * unevaluated IVP into an answer; it can never displace a method's answer.
 * MATHILDA_NO_TANGENTLINE disables it for A/B.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../zero_test.h"
#include <stdlib.h>

/* True when e (evaluated) carries no infinity / Indeterminate. */
static bool tl_finite(const Expr* e) {
    return e && !ds_contains(e, SYM_Indeterminate) && !ds_contains(e, SYM_ComplexInfinity)
             && !ds_contains(e, SYM_DirectedInfinity) && !ds_contains(e, SYM_Infinity);
}

/* F(x, y) at the point (x0, y0), evaluated.  F is in x and y[x]; consumed. */
static Expr* tl_at_point(Expr* F, const char* yname, const char* xv,
                         const Expr* x0, const Expr* y0) {
    F = ds_subst(F, ds_make_funcapp(yname, 0, xv), expr_copy((Expr*)y0));
    return ds_subst(F, expr_new_symbol(xv), expr_copy((Expr*)x0));
}

Expr* dsolve_tangent_line_solution(DSolveProblem* P) {
    if (getenv("MATHILDA_NO_TANGENTLINE")) return NULL;
    if (P->is_pde || P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    if (P->ncond != 1 || P->conds[0].fi != 0 || P->conds[0].order != 0) return NULL;
    const char* yname = P->fun_names[0];
    const char* xv    = P->ind_names[0];
    const Expr* x0 = P->conds[0].point;
    const Expr* y0 = P->conds[0].value;
    if (!ds_free_of(x0, xv) || !ds_free_of(y0, xv)) return NULL;

    Expr* F = dsolve_solve_top_derivative(P, 1);          /* y' == F(x, y[x]) */
    if (!F) return NULL;

    /* Uniqueness: F_y finite at the point (local Lipschitz condition). */
    const char* Yn = intern_symbol("DSolve`tlY");
    Expr* FY = ds_subst(expr_copy(F), ds_make_funcapp(yname, 0, xv), expr_new_symbol(Yn));
    Expr* Fy = ds_d(FY, expr_new_symbol(Yn));             /* consumes FY */
    Fy = ds_subst(Fy, expr_new_symbol(Yn), ds_make_funcapp(yname, 0, xv));
    Expr* fy0 = tl_at_point(Fy, yname, xv, x0, y0);
    bool lipschitz = tl_finite(fy0);
    expr_free(fy0);
    if (!lipschitz) { expr_free(F); return NULL; }

    /* The slope at the point; the candidate line through it. */
    Expr* m = tl_at_point(expr_copy(F), yname, xv, x0, y0);
    if (!tl_finite(m) || !ds_free_of(m, xv)) { expr_free(m); expr_free(F); return NULL; }
    Expr* line = eval_and_free(ds_call2(SYM_Plus, expr_copy((Expr*)y0),
                     ds_call2(SYM_Times, expr_copy(m),
                         ds_call2(SYM_Plus, expr_new_symbol(xv),
                             ds_call2(SYM_Times, expr_new_integer(-1), expr_copy((Expr*)x0))))));

    /* Strict check: F(x, L(x)) - m must be PROVED zero.  The substrate's own
     * verifier keeps a branch it cannot decide, which is right for a method's
     * derived family but not for a candidate that was only guessed. */
    Expr* r = ds_subst(F, ds_make_funcapp(yname, 0, xv), expr_copy(line));   /* consumes F */
    r = eval_and_free(ds_call2(SYM_Plus, r, ds_call2(SYM_Times, expr_new_integer(-1), m)));
    bool ok = zero_test_decide(r) == ZERO_TEST_TRUE;
    expr_free(r);
    if (!ok) { expr_free(line); return NULL; }
    return line;
}

Expr** dsolve_tangentline_try(DSolveProblem* P, size_t* nbranch);
Expr** dsolve_tangentline_try(DSolveProblem* P, size_t* nbranch) {
    Expr* line = dsolve_tangent_line_solution(P);
    if (!line) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = line;
    *nbranch = 1;
    return out;
}
