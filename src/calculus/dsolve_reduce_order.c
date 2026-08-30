/*
 * dsolve_reduce_order.c — DSolve`ReductionOfOrder.
 *
 * Reduces a second-order ODE that is missing the dependent variable y — i.e.
 * y''[x] == F(x, y'[x]) with F free of y — to a first-order equation in p = y':
 *     p'[x] == F(x, p),
 * solved by recursing into the scalar engine.  The second solution then comes
 * from y = Integrate[p, x] + C[2].
 *
 * (The autonomous case y'' == F(y, y') missing x, and the energy integral
 * y'' == f(y), are left for a later pass — they typically yield elliptic /
 * implicit forms.)
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
static Expr* extract_applied(Expr* r, const char* pfun) {
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

Expr** dsolve_reduce_order_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 2) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 2);   /* y'' == F(x, y, y') */
    if (!F) return NULL;
    /* missing y: F must be free of the order-0 funcapp y[x].  (ds_contains on the
     * bare symbol y would false-positive on the y inside y'[x] = Derivative[1][y][x],
     * so detect the funcapp by substituting it for a marker.) */
    const char* ychk = intern_symbol("DSolve`ropY");
    Expr* Fchk = ds_subst(expr_copy(F), ds_make_funcapp(yname, 0, xvar), expr_new_symbol(ychk));
    bool has_y = ds_contains(Fchk, ychk);
    expr_free(Fchk);
    if (has_y) { expr_free(F); return NULL; }

    /* replace y'[x] with p[x], forming the reduced equation p'[x] == F(x, p) */
    const char* pfun = intern_symbol("DSolve`rop");
    Expr* Fp = ds_subst(F, ds_make_funcapp(yname, 1, xvar), ds_make_funcapp(pfun, 0, xvar));
    Expr* reduced = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ ds_make_funcapp(pfun, 1, xvar), Fp }, 2);

    /* solve the first-order equation for p (applied form) */
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ reduced, ds_make_funcapp(pfun, 0, xvar), expr_new_symbol(xvar) }, 3);
    Expr* r = eval_and_free(call);
    Expr* pbody = extract_applied(r, pfun);
    expr_free(r);
    if (!pbody) return NULL;

    /* y = Integrate[p, x] + C[2] (p already carries C[1]) */
    Expr* yint = ds_integrate(expr_copy(pbody), expr_new_symbol(xvar));
    if (ds_has_head(yint, SYM_Integrate)) { expr_free(yint); expr_free(pbody); return NULL; }
    /* Guard against a wrong antiderivative (the Integrate engine returns 0 for
     * some rational integrands with symbolic parameters, which would silently
     * yield the degenerate y = const): require D[yint, x] == p. */
    Expr* chk = eval_and_free(ds_call2(SYM_Subtract,
                    ds_d(expr_copy(yint), expr_new_symbol(xvar)), expr_copy(pbody)));
    bool intok = ds_is_zero(chk);
    expr_free(chk); expr_free(pbody);
    if (!intok) { expr_free(yint); return NULL; }
    Expr* body = eval_and_free(ds_call2(SYM_Plus, yint, ds_const(2)));

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_reduce_order(Expr* res) {
    return dsolve_method_builtin(res, dsolve_reduce_order_try);
}

void dsolve_reduce_order_init(void) {
    symtab_add_builtin("DSolve`ReductionOfOrder", builtin_dsolve_reduce_order);
    symtab_get_def("DSolve`ReductionOfOrder")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("DSolve`ReductionOfOrder",
        "DSolve`ReductionOfOrder[eqn, y, x] reduces a second-order ODE missing y "
        "(y'' == F(x, y')) to first order in p = y', solves it, and integrates.");
}
