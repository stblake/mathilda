/*
 * dsolve_riccati.c — DSolve`Riccati.
 *
 * Solves the Riccati equation  y'[x] == q0(x) + q1(x) y + q2(x) y^2  (q2 != 0)
 * by the classical linearising substitution  y = -u'/(q2 u), which turns it into
 * the homogeneous second-order linear ODE
 *
 *     u'' - (q1 + q2'/q2) u' + (q0 q2) u == 0.
 *
 * That equation is handed back to the scalar DSolve cascade (const-coeff, Euler,
 * Airy/Bessel special forms, Kovacic, or — last — a Frobenius/power series), and
 * the recovered u is mapped back through y = -u'/(q2 u).  The general Riccati
 * solution depends only on the ratio C[1]:C[2] of the two linear constants, so
 * C[2] is collapsed to 1, leaving the single arbitrary constant C[1] (this omits
 * only the lone C[1]->Infinity particular sheet u = u1, the standard
 * one-parameter parametrisation, matching Mathematica).  When the linearised ODE
 * has no closed form the recursion declines and Riccati declines — no
 * implicit/wrong answer is shipped, and the substrate's back-substitution verify
 * is the authoritative backstop (permissive: an undecidable special-function
 * residual is kept, exactly as for Bessel/Kovacic).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

/* base^-1; base consumed. */
static Expr* powneg1(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

/* Extract the RHS of {{u[x] -> expr}} (applied form). */
static Expr* extract_applied(Expr* r, const char* ufun) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL
                && lhs->data.function.head->data.symbol.name == ufun)
                return expr_copy(rule->data.function.args[1]);
        }
    }
    return NULL;
}

Expr** dsolve_riccati_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);        /* y' == F(x, y) */
    if (!F) return NULL;

    const char* Yn = intern_symbol("DSolve`Y");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    /* quadratic coefficients:  q0 = FY|_{Y=0},  q1 = FY_Y|_{Y=0},  q2 = FY_YY/2 */
    Expr* dFY = ds_d(expr_copy(FY), expr_new_symbol(Yn));
    Expr* q0  = ds_subst(expr_copy(FY),  expr_new_symbol(Yn), expr_new_integer(0));
    Expr* q1  = ds_subst(expr_copy(dFY), expr_new_symbol(Yn), expr_new_integer(0));
    Expr* q2  = eval_and_free(ds_call2(SYM_Times,
                    ds_d(dFY, expr_new_symbol(Yn)),         /* dFY consumed here */
                    powneg1(expr_new_integer(2))));         /* * 1/2 */

    /* genuine Riccati: coefficients free of Y, q2 != 0, and FY reconstructs
     * exactly as q0 + q1 Y + q2 Y^2 (rejects a cubic-or-higher / rational-in-Y
     * nonlinearity that no first-order quadratic form can carry). */
    bool ok = ds_free_of(q0, Yn) && ds_free_of(q1, Yn) && ds_free_of(q2, Yn)
              && !ds_is_zero(q2);
    if (ok) {
        Expr* Y2 = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(Yn), expr_new_integer(2)));
        Expr* recon = eval_and_free(ds_call2(SYM_Plus,
            expr_copy(q0),
            ds_call2(SYM_Plus,
                ds_call2(SYM_Times, expr_copy(q1), expr_new_symbol(Yn)),
                ds_call2(SYM_Times, expr_copy(q2), Y2))));
        Expr* chk = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY), recon));
        ok = ds_is_zero(chk);
        expr_free(chk);
    }
    expr_free(FY);
    if (!ok) { expr_free(q0); expr_free(q1); expr_free(q2); return NULL; }

    /* build the u-ODE:  u'' - (q1 + q2'/q2) u' + (q0 q2) u == 0 */
    const char* ufun = intern_symbol("DSolve`ricU");
    Expr* q2p = ds_d(expr_copy(q2), expr_new_symbol(xvar));                       /* q2' */
    Expr* b1  = eval_and_free(ds_call2(SYM_Plus, expr_copy(q1),
                    ds_call2(SYM_Times, q2p, powneg1(expr_copy(q2)))));           /* q1 + q2'/q2 */
    Expr* b0  = eval_and_free(ds_call2(SYM_Times, expr_copy(q0), expr_copy(q2))); /* q0 q2 */

    Expr* t2 = ds_make_funcapp(ufun, 2, xvar);
    Expr* t1 = ds_call2(SYM_Times,
                   ds_call2(SYM_Times, expr_new_integer(-1), b1),
                   ds_make_funcapp(ufun, 1, xvar));
    Expr* t0 = ds_call2(SYM_Times, b0, ds_make_funcapp(ufun, 0, xvar));
    Expr* lhs = ds_call2(SYM_Plus, t2, ds_call2(SYM_Plus, t1, t0));
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
                    (Expr*[]){ lhs, expr_new_integer(0) }, 2);

    /* recurse into the scalar cascade for u (applied form -> {{u[x] -> body}}) */
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, ds_make_funcapp(ufun, 0, xvar), expr_new_symbol(xvar) }, 3);
    Expr* r = eval_and_free(call);
    Expr* ubody = extract_applied(r, ufun);
    expr_free(r);
    if (!ubody) { expr_free(q0); expr_free(q1); expr_free(q2); return NULL; }

    /* collapse the redundant constant (Riccati has one parameter): C[2] -> 1 */
    ubody = ds_subst(ubody, ds_const(2), expr_new_integer(1));

    /* y = -u'/(q2 u) */
    Expr* up   = ds_d(expr_copy(ubody), expr_new_symbol(xvar));
    Expr* num  = ds_call2(SYM_Times, expr_new_integer(-1), up);
    Expr* den  = eval_and_free(ds_call2(SYM_Times, expr_copy(q2), ubody));   /* consumes ubody */
    Expr* body = ds_simplify(eval_and_free(ds_call2(SYM_Times, num, powneg1(den))));

    expr_free(q0); expr_free(q1); expr_free(q2);

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_riccati(Expr* res) {
    return dsolve_method_builtin(res, dsolve_riccati_try);
}

void dsolve_riccati_init(void) {
    symtab_add_builtin("DSolve`Riccati", builtin_dsolve_riccati);
    symtab_get_def("DSolve`Riccati")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Riccati",
        "DSolve`Riccati[eqn, y, x] solves the Riccati equation "
        "y'[x] == q0(x) + q1(x) y + q2(x) y^2 (q2 != 0) by the substitution "
        "y = -u'/(q2 u), which linearises it to a second-order linear ODE that "
        "the scalar cascade solves; the recovered u is mapped back to y.");
}
