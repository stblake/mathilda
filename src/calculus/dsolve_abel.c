/*
 * dsolve_abel.c — DSolve`Abel (Abel equation of the first kind).
 *
 * Solves  y' = f3(x) y^3 + f2(x) y^2 + f1(x) y + f0(x)  (f3 != 0, f2 != 0) by
 * removing the quadratic term with the shift  z = y + f2/(3 f3), which gives a
 * Chini equation of exponent 3:
 *
 *     z' = f3 z^3 + P z + Q,
 *     P = 3 f3 beta^2 + 2 f2 beta + f1,
 *     Q = f3 beta^3 + f2 beta^2 + f1 beta + f0 - beta',   beta = -f2/(3 f3),
 *
 * then hands (f3, P, Q, n=3, z = y[x] - beta) to the shared
 * dsolve_chini_first_integral() and returns its implicit first integral through
 * dsolve_run_implicit.  The case f2 == 0 is already a Chini equation and is left
 * to DSolve`Chini (which runs first).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

extern Expr* dsolve_chini_first_integral(const Expr* f, const Expr* g, const Expr* h,
                                         const Expr* nexp, const Expr* yvar_expr,
                                         const char* xvar);

static Expr* powneg1(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ base, expr_new_integer(-1) }, 2);
}
/* head[a,b] evaluated; a,b consumed. */
static Expr* ev2(const char* head, Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b }, 2));
}

Expr** dsolve_abel_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;
    const char* Yn = intern_symbol("DSolve`Y");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    /* f0 = FY|0, f1 = FY_Y|0, f2 = FY_YY|0 / 2, f3 = FY_YYY|0 / 6 */
    Expr* dF  = ds_d(expr_copy(FY), expr_new_symbol(Yn));
    Expr* d2F = ds_d(expr_copy(dF), expr_new_symbol(Yn));
    Expr* d3F = ds_d(expr_copy(d2F), expr_new_symbol(Yn));
    Expr* f0 = ds_simplify(ds_subst(expr_copy(FY), expr_new_symbol(Yn), expr_new_integer(0)));
    Expr* f1 = ds_simplify(ds_subst(dF,  expr_new_symbol(Yn), expr_new_integer(0)));
    Expr* f2 = ds_simplify(ev2(SYM_Times, ds_subst(d2F, expr_new_symbol(Yn), expr_new_integer(0)),
                               powneg1(expr_new_integer(2))));
    Expr* f3 = ds_simplify(ev2(SYM_Times, ds_subst(d3F, expr_new_symbol(Yn), expr_new_integer(0)),
                               powneg1(expr_new_integer(6))));

    bool ok = ds_free_of(f0,Yn) && ds_free_of(f1,Yn) && ds_free_of(f2,Yn) && ds_free_of(f3,Yn)
              && !ds_is_zero(f3) && !ds_is_zero(f2);
    if (ok) {
        Expr* recon = eval_and_free(ds_call2(SYM_Plus, expr_copy(f0),
            ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_copy(f1), expr_new_symbol(Yn)),
            ds_call2(SYM_Plus,
                ds_call2(SYM_Times, expr_copy(f2), ev2(SYM_Power, expr_new_symbol(Yn), expr_new_integer(2))),
                ds_call2(SYM_Times, expr_copy(f3), ev2(SYM_Power, expr_new_symbol(Yn), expr_new_integer(3)))))));
        Expr* chk = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY), recon));
        ok = ds_is_zero(chk);
        expr_free(chk);
    }
    expr_free(FY);
    if (!ok) { expr_free(f0); expr_free(f1); expr_free(f2); expr_free(f3); return NULL; }

    /* beta = -f2/(3 f3);  reduced Chini-3 coefficients P, Q */
    Expr* beta  = ds_simplify(ev2(SYM_Times, ev2(SYM_Times, expr_new_integer(-1), expr_copy(f2)),
                                  powneg1(ev2(SYM_Times, expr_new_integer(3), expr_copy(f3)))));
    Expr* betap = ds_d(expr_copy(beta), expr_new_symbol(xvar));
    Expr* beta2 = ev2(SYM_Power, expr_copy(beta), expr_new_integer(2));
    Expr* beta3 = ev2(SYM_Power, expr_copy(beta), expr_new_integer(3));

    /* P = 3 f3 beta^2 + 2 f2 beta + f1 */
    Expr* Pc = ds_simplify(eval_and_free(ds_call2(SYM_Plus,
        ds_call2(SYM_Times, ev2(SYM_Times, expr_new_integer(3), expr_copy(f3)), expr_copy(beta2)),
        ds_call2(SYM_Plus,
            ds_call2(SYM_Times, ev2(SYM_Times, expr_new_integer(2), expr_copy(f2)), expr_copy(beta)),
            expr_copy(f1)))));
    /* Q = f3 beta^3 + f2 beta^2 + f1 beta + f0 - beta' */
    Expr* Qc = ds_simplify(eval_and_free(ds_call2(SYM_Plus,
        ds_call2(SYM_Times, expr_copy(f3), expr_copy(beta3)),
        ds_call2(SYM_Plus,
            ds_call2(SYM_Times, expr_copy(f2), expr_copy(beta2)),
            ds_call2(SYM_Plus,
                ds_call2(SYM_Times, expr_copy(f1), expr_copy(beta)),
                ds_call2(SYM_Subtract, expr_copy(f0), betap))))));   /* betap consumed */
    expr_free(beta2); expr_free(beta3);

    /* z = y[x] - beta  (= y + f2/(3 f3)) */
    Expr* zexpr = eval_and_free(ds_call2(SYM_Subtract, ds_make_funcapp(yname, 0, xvar), expr_copy(beta)));
    Expr* three = expr_new_integer(3);
    Expr* G = dsolve_chini_first_integral(f3, Pc, Qc, three, zexpr, xvar);
    expr_free(three); expr_free(zexpr); expr_free(beta);
    expr_free(Pc); expr_free(Qc); expr_free(f0); expr_free(f1); expr_free(f2); expr_free(f3);
    if (!G) return NULL;

    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_abel(Expr* res) {
    return dsolve_method_builtin_implicit(res, dsolve_abel_try);
}

void dsolve_abel_init(void) {
    symtab_add_builtin("DSolve`Abel", builtin_dsolve_abel);
    symtab_get_def("DSolve`Abel")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Abel",
        "DSolve`Abel[eqn, y, x] solves the first-kind Abel equation "
        "y' == f3 y^3 + f2 y^2 + f1 y + f0 by removing the y^2 term (z = y + f2/(3 f3)) "
        "to a Chini equation of exponent 3; returns the implicit first integral when "
        "that reduces to an autonomous form.");
}
