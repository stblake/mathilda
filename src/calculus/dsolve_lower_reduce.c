/*
 * dsolve_lower_reduce.c — DSolve`LowerDerivativeReduction.
 *
 * Reduces an ODE that is missing the dependent variable AND every derivative
 * below some order m >= 1 — i.e. only y^(m), y^(m+1), ..., y^(n) appear (n >= 3)
 * — by substituting p = y^(m).  The equation becomes an order-(n-m) ODE in p,
 * solved by recursing into the scalar engine; y is then recovered by integrating
 * p back m times (each integration introduces one constant, so the m constants
 * build the expected degree-(m-1) polynomial automatically).
 *
 * This is the general order >= 3 counterpart of ReductionOfOrder (which handles
 * the order-2, m=1 missing-y case) and AutonomousReduction (order-2 missing-x).
 * Gated to n >= 3 so it never competes with those at order 2.  Example:
 *   7 y' y''' - 11 y''^2 == 0   (only y', y'', y''' appear, m = 1)
 *     -> 7 p p'' - 11 p'^2 == 0  (autonomous, missing x) -> p = (4x+7C1)^(-7/4)
 *     -> y = Integrate[p, x] + C[3].
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

/* Extract the RHS of {{p[x] -> expr}} (applied form). */
static Expr* lr_extract_applied(Expr* r, const char* pfun) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL
                && lhs->data.function.head->data.symbol.name == pfun)
                return expr_copy(rule->data.function.args[1]);
        }
    }
    return NULL;
}

Expr** dsolve_lower_reduce_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    int n = P->max_order[0];
    if (n < 3) return NULL;   /* order-2 missing-var is owned by reduce_order/autonomous */
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];
    Expr* resid = P->eq_residuals[0];

    /* m = lowest derivative order actually present.  Detect each order by
     * substituting its funcapp for a marker (a bare-symbol scan would false-
     * positive on the y inside every Derivative[k][y][x]). */
    const char* mark = intern_symbol("DSolve`lrmark");
    int m = -1;
    for (int k = 0; k <= n && m < 0; k++) {
        Expr* t = ds_subst(expr_copy(resid), ds_make_funcapp(yname, k, xvar), expr_new_symbol(mark));
        bool present = ds_contains(t, mark);
        expr_free(t);
        if (present) m = k;
    }
    if (m < 1) return NULL;   /* y[x] itself is present (or nothing) -> not our case */

    /* substitute y^(m+j) -> p^(j) (j = n-m .. 0) to form the order-(n-m) p-ODE */
    const char* pfun = intern_symbol("DSolve`lrp");
    Expr* red = expr_copy(resid);
    for (int j = n - m; j >= 0; j--)
        red = ds_subst(red, ds_make_funcapp(yname, m + j, xvar), ds_make_funcapp(pfun, j, xvar));
    Expr* reduced = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ red, expr_new_integer(0) }, 2);

    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ reduced, ds_make_funcapp(pfun, 0, xvar), expr_new_symbol(xvar) }, 3);
    Expr* r = eval_and_free(call);
    Expr* body = lr_extract_applied(r, pfun);   /* body = y^(m), carries C[1..(n-m)] */
    expr_free(r);
    if (!body) return NULL;

    /* integrate p back m times; the reduced ODE used C[1..(n-m)], so the m new
     * integration constants are C[(n-m)+1 .. n].  Each antiderivative is guarded
     * (D[Integrate[.]] == integrand) exactly as ReductionOfOrder does. */
    int base = n - m;
    for (int i = 1; i <= m; i++) {
        Expr* yint = ds_integrate(expr_copy(body), expr_new_symbol(xvar));
        if (ds_has_head(yint, SYM_Integrate)) { expr_free(yint); expr_free(body); return NULL; }
        Expr* chk = eval_and_free(ds_call2(SYM_Subtract,
                        ds_d(expr_copy(yint), expr_new_symbol(xvar)), expr_copy(body)));
        bool intok = ds_is_zero(chk);
        if (!intok) {
            Expr* pz = eval_and_free(ds_call1("PossibleZeroQ", expr_copy(chk)));
            intok = (pz->type == EXPR_SYMBOL && pz->data.symbol.name == SYM_True);
            expr_free(pz);
        }
        expr_free(chk); expr_free(body);
        if (!intok) { expr_free(yint); return NULL; }
        body = eval_and_free(ds_call2(SYM_Plus, yint, ds_const(base + i)));
    }

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_lower_reduce(Expr* res) {
    return dsolve_method_builtin(res, dsolve_lower_reduce_try);
}

void dsolve_lower_reduce_init(void) {
    symtab_add_builtin("DSolve`LowerDerivativeReduction", builtin_dsolve_lower_reduce);
    symtab_get_def("DSolve`LowerDerivativeReduction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LowerDerivativeReduction",
        "DSolve`LowerDerivativeReduction[eqn, y, x] solves an order n >= 3 ODE in "
        "which y and every derivative below some order m >= 1 are absent: it sets "
        "p = y^(m), solves the resulting order-(n-m) ODE, and integrates p back m "
        "times.");
}
