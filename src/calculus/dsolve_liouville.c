/*
 * dsolve_liouville.c — DSolve`Liouville.
 *
 * Solves the Liouville equation  y'' + g(y) (y')^2 + h(x) y' == 0  by Liouville's
 * transformation.  Dividing by y' gives  d/dx[ln y' + G(y)] == -h(x)  with
 * G = Integrate[g, y], so  ln y' + G(y) == -H(x) + c1  (H = Integrate[h, x]),
 * hence  Exp[G(y)] dy == e^{c1} Exp[-H(x)] dx  and, integrating both sides,
 *     Integrate[Exp[G], y] == C[1] Integrate[Exp[-H], x] + C[2].
 * Solving that relation for y gives the two-constant general solution.
 *
 * Distinct from AutonomousReduction (which needs h == 0, no explicit x) and
 * ReductionOfOrder (which needs g == 0, no y).  Runs in the 2nd-order nonlinear
 * group, after AutonomousReduction.  Mirrors SymPy's `Liouville`.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* Coefficient[e, p, k]; e borrowed, result owned. */
static Expr* coeff_p(const Expr* e, const char* p, int k) {
    return eval_and_free(expr_new_function(expr_new_symbol("Coefficient"),
        (Expr*[]){ expr_copy((Expr*)e), expr_new_symbol(p), expr_new_integer(k) }, 3));
}

Expr** dsolve_liouville_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 2) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`lvY");
    const char* pn = intern_symbol("DSolve`lvP");

    Expr* F = dsolve_solve_top_derivative(P, 2);          /* y'' == F(x, y, y') */
    if (!F) return NULL;
    /* F with y'[x] -> p, y[x] -> Y */
    Expr* Fp = ds_subst(F, ds_make_funcapp(yname, 1, xvar), expr_new_symbol(pn));
    Fp = ds_subst(Fp, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    /* require F == -g(Y) p^2 - h(x) p (degree-2 in p, no constant term) */
    Expr* c2 = coeff_p(Fp, pn, 2);
    Expr* c1 = coeff_p(Fp, pn, 1);
    Expr* c0 = coeff_p(Fp, pn, 0);
    Expr* recon = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Fp),
                      ds_call2(SYM_Plus,
                          ds_call2(SYM_Plus,
                              ds_call2(SYM_Times, expr_copy(c2),
                                  expr_new_function(expr_new_symbol(SYM_Power),
                                      (Expr*[]){ expr_new_symbol(pn), expr_new_integer(2) }, 2)),
                              ds_call2(SYM_Times, expr_copy(c1), expr_new_symbol(pn))),
                          expr_copy(c0))));
    bool okform = ds_is_zero(recon) && ds_is_zero(c0);
    expr_free(recon); expr_free(Fp); expr_free(c0);
    /* g = -c2 (free of x), h = -c1 (free of Y) */
    Expr* g = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), c2));
    Expr* h = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), c1));
    if (!okform || !ds_free_of(g, xvar) || !ds_free_of(h, Yn) || ds_is_zero(g)) {
        expr_free(g); expr_free(h); return NULL;
    }

    /* EG = Integrate[Exp[Integrate[g, Y]], Y] ; EH = Integrate[Exp[-Integrate[h,x]], x] */
    Expr* G  = ds_integrate(g, expr_new_symbol(Yn));
    if (ds_has_head(G, SYM_Integrate)) { expr_free(G); expr_free(h); return NULL; }
    Expr* EG = ds_integrate(ds_call1("Exp", G), expr_new_symbol(Yn));
    if (ds_has_head(EG, SYM_Integrate)) { expr_free(EG); expr_free(h); return NULL; }

    Expr* H  = ds_integrate(h, expr_new_symbol(xvar));
    if (ds_has_head(H, SYM_Integrate)) { expr_free(H); expr_free(EG); return NULL; }
    Expr* mH = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), H));
    Expr* EH = ds_integrate(ds_call1("Exp", mH), expr_new_symbol(xvar));
    if (ds_has_head(EH, SYM_Integrate)) { expr_free(EH); expr_free(EG); return NULL; }

    /* EG(Y) == C[1] EH(x) + C[2] ; solve for Y */
    Expr* rhs = eval_and_free(ds_call2(SYM_Plus,
                    ds_call2(SYM_Times, ds_const(1), EH), ds_const(2)));
    Expr* eq  = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ EG, rhs }, 2);
    Expr* sol = ds_solve(eq, expr_new_symbol(Yn));
    size_t nb = 0;
    Expr** bodies = dsolve_extract_solutions(sol, Yn, &nb);
    if (sol) expr_free(sol);
    if (!bodies) return NULL;

    *nbranch = nb;
    return bodies;
}

static Expr* builtin_dsolve_liouville(Expr* res) {
    return dsolve_method_builtin(res, dsolve_liouville_try);
}

void dsolve_liouville_init(void) {
    symtab_add_builtin("DSolve`Liouville", builtin_dsolve_liouville);
    symtab_get_def("DSolve`Liouville")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Liouville",
        "DSolve`Liouville[eqn, y, x] solves y'' + g(y)(y')^2 + h(x)y' == 0 by Liouville's "
        "transformation: Integrate[Exp[Integrate[g,y]], y] == C[1] Integrate[Exp[-Integrate[h,x]], x] "
        "+ C[2], solved for y.");
}
