/*
 * dsolve_chini.c — DSolve`Chini.
 *
 * Solves the Chini equation  y' = f(x) y^n + g(x) y + h(x)  (n != 0, 1; n = 2 is
 * Riccati) in the sub-class reducible to an autonomous equation.  The scaling
 * y = phi(x) u with phi = f^(-1/(n-1)) makes the y^n coefficient 1:
 *
 *     u' = u^n + B u + C,   B = g + f'/((n-1) f),   C = h f^(1/(n-1)),
 *
 * and when B, C are free of x this is autonomous (separable).  Since the
 * integrand 1/(u^n + B u + C) is a rational function of u, its antiderivative
 * Xi(u) is elementary, giving the implicit first integral
 *
 *     Xi(u) - x == C[1],   u = y f^(1/(n-1)),
 *
 * returned through dsolve_run_implicit as {{ G(x,y[x]) == C[1] }} with
 * G = Xi(y f^(1/(n-1))) - x.  Back-substituting y' = -G_x/G_y reproduces the ODE
 * exactly, so the implicit verifier accepts it.  Abel equations reuse the shared
 * dsolve_chini_first_integral() after removing their y^2 term (see dsolve_abel.c).
 */
#include "dsolve_common.h"
#include "integrate.h"
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

/* Reduction (a) — reduce to autonomous.  The scaling y = f^(-1/(n-1)) u makes
 * the y^n coefficient 1, giving u' = u^n + B u + C with B = g + f'/((n-1) f),
 * C = h f^(1/(n-1)); when B, C are free of x this is autonomous (separable) and
 * its first integral is Xi(u) - x == C[1], u = yvar f^(1/(n-1)).  Returns G or
 * NULL if B/C depend on x or the quadrature is non-elementary.  Args borrowed. */
static Expr* chini_reduce_autonomous(const Expr* f, const Expr* g, const Expr* h,
                                     const Expr* nexp, const Expr* nm1,
                                     const Expr* yvar_expr, const char* xvar) {
    /* fpow = f^(1/(n-1)) = 1/phi;  u_expr = yvar * fpow */
    Expr* fpow = eval_and_free(ds_call2(SYM_Power, expr_copy((Expr*)f), powneg1(expr_copy((Expr*)nm1))));
    Expr* u_expr = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)yvar_expr), expr_copy(fpow)));

    /* B = g + f'/((n-1) f) */
    Expr* fp = ds_d(expr_copy((Expr*)f), expr_new_symbol(xvar));
    Expr* Bc = ds_simplify(eval_and_free(ds_call2(SYM_Plus, expr_copy((Expr*)g),
                   ds_call2(SYM_Times, fp,
                       powneg1(eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)nm1), expr_copy((Expr*)f))))))));
    if (!ds_free_of(Bc, xvar)) { expr_free(Bc); expr_free(fpow); expr_free(u_expr); return NULL; }

    /* C = h * fpow */
    Expr* Cc = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)h), expr_copy(fpow))));
    expr_free(fpow);
    if (!ds_free_of(Cc, xvar)) { expr_free(Cc); expr_free(Bc); expr_free(u_expr); return NULL; }

    /* denom = u^n + B u + C  (in dummy symbol u) */
    const char* uSym = intern_symbol("DSolve`chiU");
    Expr* denom = eval_and_free(ds_call2(SYM_Plus,
        eval_and_free(ds_call2(SYM_Power, expr_new_symbol(uSym), expr_copy((Expr*)nexp))),
        ds_call2(SYM_Plus,
            ds_call2(SYM_Times, Bc, expr_new_symbol(uSym)),      /* consumes Bc */
            Cc)));                                               /* consumes Cc */

    /* Xi = Integrate[1/denom, u] (speculative: a non-elementary result declines) */
    g_integrate_quiet++;
    Expr* Xi = ds_integrate(powneg1(denom), expr_new_symbol(uSym));
    g_integrate_quiet--;
    if (ds_has_head(Xi, SYM_Integrate)) { expr_free(Xi); expr_free(u_expr); return NULL; }

    /* G = (Xi /. u -> u_expr) - x.  NOT ds_simplify'd: Simplify chokes (hangs) on
     * a composed Log/ArcTan with radical arguments (e.g. Sqrt[x] (y+1)), and the
     * relation is correct unsimplified — the implicit verifier differentiates it. */
    Expr* Xisub = ds_subst(Xi, expr_new_symbol(uSym), u_expr);   /* consumes Xi, u_expr */
    return eval_and_free(ds_call2(SYM_Subtract, Xisub, expr_new_symbol(xvar)));
}

/* Reduction (b) — remove the linear term.  The scaling y = mu w with
 * mu = Exp[Integrate[g, x]] eliminates the g y term, giving w' = P w^n + R with
 * P = f mu^(n-1), R = h/mu.  When R/P is constant (= c) the equation is
 * SEPARABLE: Integrate[1/(w^n + c), w] == Integrate[P, x] + C[1], so the first
 * integral is Xi(w) - Integrate[P, x] == C[1] with w = yvar/mu.  This covers the
 * Chini sub-class where (a) fails, e.g. y' = x E^(2x) y^3 - y - x E^-x (there
 * B = 1/(2x) is not constant, but c = -1 is).  Returns G or NULL.  Args
 * borrowed. */
static Expr* chini_reduce_linremoval(const Expr* f, const Expr* g, const Expr* h,
                                     const Expr* nexp, const Expr* nm1,
                                     const Expr* yvar_expr, const char* xvar) {
    /* Ig = Integrate[g, x];  mu = Exp[Ig], muinv = Exp[-Ig] */
    g_integrate_quiet++;
    Expr* Ig = ds_integrate(expr_copy((Expr*)g), expr_new_symbol(xvar));
    g_integrate_quiet--;
    if (ds_has_head(Ig, SYM_Integrate)) { expr_free(Ig); return NULL; }   /* Integrate[g] non-elementary */

    Expr* mu    = eval_and_free(ds_call1("Exp", expr_copy(Ig)));
    Expr* muinv = eval_and_free(ds_call1("Exp",
                      eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), Ig))));   /* consumes Ig */

    /* P = f mu^(n-1),  R = h muinv,  c = R/P (must be free of x for separability) */
    Expr* Pc = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)f),
                   eval_and_free(ds_call2(SYM_Power, expr_copy(mu), expr_copy((Expr*)nm1))))));
    Expr* Rc = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)h), expr_copy(muinv))));
    Expr* c  = ds_simplify(eval_and_free(ds_call2(SYM_Times, Rc, powneg1(expr_copy(Pc)))));  /* consumes Rc */
    if (!ds_free_of(c, xvar)) { expr_free(c); expr_free(Pc); expr_free(mu); expr_free(muinv); return NULL; }

    /* w_expr = yvar * muinv  (w = y/mu) */
    Expr* w_expr = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)yvar_expr), expr_copy(muinv)));
    expr_free(mu); expr_free(muinv);

    /* Xi = Integrate[1/(u^n + c), u] */
    const char* uSym = intern_symbol("DSolve`chiU");
    Expr* denom = eval_and_free(ds_call2(SYM_Plus,
        eval_and_free(ds_call2(SYM_Power, expr_new_symbol(uSym), expr_copy((Expr*)nexp))),
        c));                                                     /* consumes c */
    g_integrate_quiet++;
    Expr* Xi = ds_integrate(powneg1(denom), expr_new_symbol(uSym));
    g_integrate_quiet--;
    if (ds_has_head(Xi, SYM_Integrate)) { expr_free(Xi); expr_free(w_expr); expr_free(Pc); return NULL; }

    /* Rhs = Integrate[P, x] */
    g_integrate_quiet++;
    Expr* Rhs = ds_integrate(Pc, expr_new_symbol(xvar));         /* consumes Pc */
    g_integrate_quiet--;
    if (ds_has_head(Rhs, SYM_Integrate)) { expr_free(Rhs); expr_free(Xi); expr_free(w_expr); return NULL; }

    /* G = (Xi /. u -> w_expr) - Rhs  (unsimplified; the implicit verifier differentiates it) */
    Expr* Xisub = ds_subst(Xi, expr_new_symbol(uSym), w_expr);   /* consumes Xi, w_expr */
    return eval_and_free(ds_call2(SYM_Subtract, Xisub, Rhs));    /* consumes Xisub, Rhs */
}

/* Shared: the implicit first integral G of the reduced Chini equation with
 * coefficients f, g, h, exponent nexp, where the dependent expression is
 * yvar_expr (y[x] for Chini; y[x]+f2/(3 f3) for Abel).  Returns G (meaning
 * G == C[1]) or NULL when neither the autonomous reduction (a) nor the
 * linear-term-removal reduction (b) applies.  f,g,h,nexp,yvar_expr are
 * borrowed. */
Expr* dsolve_chini_first_integral(const Expr* f, const Expr* g, const Expr* h,
                                  const Expr* nexp, const Expr* yvar_expr,
                                  const char* xvar) {
    Expr* nm1 = eval_and_free(ds_call2(SYM_Subtract, expr_copy((Expr*)nexp), expr_new_integer(1)));
    if (ds_is_zero(nm1)) { expr_free(nm1); return NULL; }              /* n == 1 */

    Expr* G = chini_reduce_autonomous(f, g, h, nexp, nm1, yvar_expr, xvar);
    if (!G) G = chini_reduce_linremoval(f, g, h, nexp, nm1, yvar_expr, xvar);
    expr_free(nm1);
    return G;
}

Expr** dsolve_chini_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;
    const char* Yn = intern_symbol("DSolve`Y");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    Expr* h = ds_subst(expr_copy(FY), expr_new_symbol(Yn), expr_new_integer(0));   /* FY|_{Y=0} */
    Expr* res = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY), expr_copy(h)));
    if (ds_is_zero(res)) { expr_free(res); expr_free(h); expr_free(FY); return NULL; }

    /* Q = res - Y res_Y  == (1-n) f Y^n  for res = g Y + f Y^n */
    Expr* Q = eval_and_free(ds_call2(SYM_Subtract, expr_copy(res),
                  ds_call2(SYM_Times, expr_new_symbol(Yn), ds_d(expr_copy(res), expr_new_symbol(Yn)))));
    if (ds_is_zero(Q)) { expr_free(Q); expr_free(res); expr_free(h); expr_free(FY); return NULL; }  /* linear */

    /* n = Y Q_Y / Q, a constant.  Reduce (Cancel) BEFORE the free-of guard (see
     * the Bernoulli note): the unexpanded Q makes the raw ratio textually depend
     * on x/Y, sending ds_free_of to a numeric probe that hits 0^(-1) and misreads
     * a constant exponent as x-dependent, so the classifier wrongly declines.
     * Cancel — not Simplify — because this runs on the DECLINE path of every
     * equation reaching Chini: full Simplify hangs on a symbolic radical ratio
     * (e.g. y' == Sqrt[y^4 + K] from AutonomousReduction), while Cancel is a cheap
     * rational-GCD that still reduces the genuine constant exponent. */
    Expr* nexp = eval_and_free(ds_call1("Cancel", eval_and_free(ds_call2(SYM_Times,
                     eval_and_free(ds_call2(SYM_Times, expr_new_symbol(Yn),
                                            ds_d(expr_copy(Q), expr_new_symbol(Yn)))),
                     powneg1(expr_copy(Q))))));
    if (!ds_free_of(nexp, Yn) || !ds_free_of(nexp, xvar)) {
        expr_free(nexp); expr_free(Q); expr_free(res); expr_free(h); expr_free(FY); return NULL;
    }
    /* n == 2 is Riccati (claimed earlier); n == 0/1 cannot arise here */
    { Expr* d2 = eval_and_free(ds_call2(SYM_Subtract, expr_copy(nexp), expr_new_integer(2)));
      bool is2 = ds_is_zero(d2); expr_free(d2);
      if (is2) { expr_free(nexp); expr_free(Q); expr_free(res); expr_free(h); expr_free(FY); return NULL; } }

    /* f = Q / ((1-n) Y^n),  g = (res - f Y^n)/Y */
    Expr* omn   = eval_and_free(ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(nexp)));
    Expr* Ynpow = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(Yn), expr_copy(nexp)));
    Expr* f = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_copy(Q),
                  powneg1(eval_and_free(ds_call2(SYM_Times, expr_copy(omn), expr_copy(Ynpow)))))));
    Expr* g = ds_simplify(eval_and_free(ds_call2(SYM_Times,
                  eval_and_free(ds_call2(SYM_Subtract, expr_copy(res),
                                         ds_call2(SYM_Times, expr_copy(f), expr_copy(Ynpow)))),
                  powneg1(expr_new_symbol(Yn)))));
    expr_free(Q); expr_free(omn); expr_free(res);

    bool ok = ds_free_of(f, Yn) && ds_free_of(g, Yn) && ds_free_of(h, Yn);
    if (ok) {
        Expr* recon = eval_and_free(ds_call2(SYM_Plus, expr_copy(h),
                          ds_call2(SYM_Plus,
                              ds_call2(SYM_Times, expr_copy(g), expr_new_symbol(Yn)),
                              ds_call2(SYM_Times, expr_copy(f), expr_copy(Ynpow)))));
        Expr* chk = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY), recon));
        ok = ds_is_zero(chk);
        expr_free(chk);
    }
    expr_free(FY); expr_free(Ynpow);
    if (!ok) { expr_free(f); expr_free(g); expr_free(h); expr_free(nexp); return NULL; }

    Expr* yfun = ds_make_funcapp(yname, 0, xvar);
    Expr* G = dsolve_chini_first_integral(f, g, h, nexp, yfun, xvar);
    expr_free(f); expr_free(g); expr_free(h); expr_free(nexp); expr_free(yfun);
    if (!G) return NULL;

    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_chini(Expr* res) {
    return dsolve_method_builtin_implicit(res, dsolve_chini_try);
}

void dsolve_chini_init(void) {
    symtab_add_builtin("DSolve`Chini", builtin_dsolve_chini);
    symtab_get_def("DSolve`Chini")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Chini",
        "DSolve`Chini[eqn, y, x] solves the Chini equation y' == f(x) y^n + g(x) y "
        "+ h(x) (n != 0,1,2) in the sub-class reducible to an autonomous equation by "
        "y = f^(-1/(n-1)) u; returns the implicit first integral {{G(x,y[x]) == C[1]}}.");
}
