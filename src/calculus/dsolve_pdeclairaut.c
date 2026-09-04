/*
 * dsolve_pdeclairaut.c — DSolve`PDEClairaut.
 *
 * The first-order PDE Clairaut form
 *
 *     u == v1 u_{v1} + v2 u_{v2} + f(u_{v1}, u_{v2})
 *
 * (u linear in v1, v2 with the derivatives as the coefficients, plus a function
 * of the derivatives alone).  Its complete integral is the two-parameter plane
 * family obtained by replacing each derivative with an arbitrary constant:
 *
 *     u == C[1] v1 + C[2] v2 + f(C[1], C[2]).
 *
 * With IncludeSingularSolutions -> True the singular integral (the envelope of
 * that family) is added: eliminate C[1], C[2] from
 * { v1 + f_{C1} == 0, v2 + f_{C2} == 0 } and substitute back.
 *
 * Detection works on the algebraic residual with the derivative/funcapp terms
 * replaced by plain symbols (sUx = u_{v1}, sUy = u_{v2}, sU = u): solving the
 * residual for sU must give Uexpr with d(Uexpr)/d v1 == sUx, d(Uexpr)/d v2 ==
 * sUy, and f := Uexpr - v1 sUx - v2 sUy free of v1, v2 (a function of the
 * derivatives alone).  This shape is nonlinear in the derivatives, so the linear
 * first-order PDE method declines before this one is reached.
 *
 * Verification and assembly reuse the ordinary explicit PDE path
 * (dsolve_run_pde_implicit's PDEBranches wrapper): each branch back-substitutes
 * to a zero residual — the arbitrary constants C[1], C[2] survive the
 * concrete-test-function rewrite untouched and zero_test decides the residual.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

/* Derivative[o1,o2][u][v1,v2]. */
static Expr* pdec_deriv(const char* u, int o1, int o2, const char* v1, const char* v2) {
    Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                  (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
    Expr* du = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    return expr_new_function(du, (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
}

Expr** dsolve_pdeclairaut_solve(DSolveProblem* P) {
    if (P->nfun != 1 || P->nind != 2 || P->neq != 1) return NULL;
    const char* uname = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];
    const char* sUx = intern_symbol("DSolve`pdeUx");
    const char* sUy = intern_symbol("DSolve`pdeUy");
    const char* sU  = intern_symbol("DSolve`pdeU");

    /* algebraic residual: u_{v1} -> sUx, u_{v2} -> sUy, u -> sU */
    Expr* R = expr_copy(P->eq_residuals[0]);
    R = ds_subst(R, pdec_deriv(uname, 1, 0, v1, v2), expr_new_symbol(sUx));
    R = ds_subst(R, pdec_deriv(uname, 0, 1, v1, v2), expr_new_symbol(sUy));
    R = ds_subst(R, expr_new_function(expr_new_symbol(uname),
                     (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2),
                 expr_new_symbol(sU));

    /* R linear in u with a nonzero constant coefficient: cU = dR/dsU */
    Expr* cU = ds_d(expr_copy(R), expr_new_symbol(sU));
    if (!ds_free_of(cU, sU) || !ds_free_of(cU, sUx) || !ds_free_of(cU, sUy)
        || !ds_free_of(cU, v1) || !ds_free_of(cU, v2) || ds_is_zero(cU)) {
        expr_free(cU); expr_free(R); return NULL;
    }
    /* Uexpr = -(R|sU=0) / cU  (the equation solved for u) */
    Expr* R0 = ds_subst(expr_copy(R), expr_new_symbol(sU), expr_new_integer(0));
    expr_free(R);
    Expr* Uexpr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
        expr_new_integer(-1), R0,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ cU, expr_new_integer(-1) }, 2)
    }, 3));

    /* Clairaut structure: d(Uexpr)/d v1 == sUx and d(Uexpr)/d v2 == sUy. */
    Expr* c1 = eval_and_free(ds_call2(SYM_Subtract,
                   ds_d(expr_copy(Uexpr), expr_new_symbol(v1)), expr_new_symbol(sUx)));
    Expr* c2 = eval_and_free(ds_call2(SYM_Subtract,
                   ds_d(expr_copy(Uexpr), expr_new_symbol(v2)), expr_new_symbol(sUy)));
    bool ok = ds_is_zero(c1) && ds_is_zero(c2);
    expr_free(c1); expr_free(c2);
    if (!ok) { expr_free(Uexpr); return NULL; }

    /* f(sUx,sUy) = Uexpr - v1 sUx - v2 sUy, must be free of v1, v2, sU */
    Expr* fp = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Uexpr),
                   ds_call2(SYM_Plus,
                       ds_call2(SYM_Times, expr_new_symbol(v1), expr_new_symbol(sUx)),
                       ds_call2(SYM_Times, expr_new_symbol(v2), expr_new_symbol(sUy)))));
    if (!ds_free_of(fp, v1) || !ds_free_of(fp, v2) || !ds_free_of(fp, sU)) {
        expr_free(fp); expr_free(Uexpr); return NULL;
    }

    /* branches: the complete integral first, then any singular envelopes */
    size_t nb = 0, cap = 4;
    Expr** branches = malloc(cap * sizeof(Expr*));

    /* complete integral: Uexpr /. { sUx -> C[1], sUy -> C[2] } */
    Expr* complete = ds_subst(expr_copy(Uexpr), expr_new_symbol(sUx), ds_const(1));
    complete = ds_subst(complete, expr_new_symbol(sUy), ds_const(2));
    branches[nb++] = eval_and_free(ds_call1("Expand", complete));

    /* singular integral: eliminate the constants from the envelope conditions
     * { v1 + f_a == 0, v2 + f_b == 0 } (a = sUx, b = sUy), substituting back.
     * Solved over fresh plain symbols ka, kb so Solve never sees indexed vars. */
    if (P->include_singular) {
        const char* ka = intern_symbol("DSolve`pdecA");
        const char* kb = intern_symbol("DSolve`pdecB");
        Expr* fk = ds_subst(expr_copy(fp), expr_new_symbol(sUx), expr_new_symbol(ka));
        fk = ds_subst(fk, expr_new_symbol(sUy), expr_new_symbol(kb));
        Expr* e1 = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){
            eval_and_free(ds_call2(SYM_Plus, expr_new_symbol(v1),
                ds_d(expr_copy(fk), expr_new_symbol(ka)))), expr_new_integer(0) }, 2);
        Expr* e2 = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){
            eval_and_free(ds_call2(SYM_Plus, expr_new_symbol(v2),
                ds_d(expr_copy(fk), expr_new_symbol(kb)))), expr_new_integer(0) }, 2);
        expr_free(fk);
        Expr* eqs  = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ e1, e2 }, 2);
        Expr* vars = expr_new_function(expr_new_symbol(SYM_List),
                         (Expr*[]){ expr_new_symbol(ka), expr_new_symbol(kb) }, 2);
        Expr* ssol = ds_solve(eqs, vars);
        /* body template in ka, kb: ka v1 + kb v2 + f(ka,kb) */
        Expr* utmpl = ds_subst(expr_copy(Uexpr), expr_new_symbol(sUx), expr_new_symbol(ka));
        utmpl = ds_subst(utmpl, expr_new_symbol(sUy), expr_new_symbol(kb));
        if (ssol && head_is(ssol, SYM_List)) {
            for (size_t i = 0; i < ssol->data.function.arg_count; i++) {
                Expr* br = ssol->data.function.args[i];   /* {ka->.., kb->..} */
                if (!head_is(br, SYM_List)) continue;
                Expr* body = expr_copy(utmpl);
                for (size_t j = 0; j < br->data.function.arg_count; j++) {
                    Expr* rule = br->data.function.args[j];
                    if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2
                        && rule->data.function.args[0]->type == EXPR_SYMBOL) {
                        const char* lhs = rule->data.function.args[0]->data.symbol.name;
                        body = ds_subst(body, expr_new_symbol(lhs),
                                        expr_copy(rule->data.function.args[1]));
                    }
                }
                if (nb >= cap) { cap *= 2; branches = realloc(branches, cap * sizeof(Expr*)); }
                branches[nb++] = eval_and_free(ds_call1("Simplify", body));
            }
        }
        expr_free(utmpl);
        if (ssol) expr_free(ssol);
    }

    expr_free(fp); expr_free(Uexpr);

    /* wrap the branches in DSolve`PDEBranches[...] for the extended PDE runner */
    Expr* wrap = expr_new_function(expr_new_symbol(intern_symbol("DSolve`PDEBranches")),
                                   branches, nb);
    free(branches);
    Expr** out = malloc(sizeof(Expr*));
    out[0] = wrap;
    return out;
}

static Expr* builtin_dsolve_pdeclairaut(Expr* res) {
    return dsolve_method_builtin_pde_implicit(res, dsolve_pdeclairaut_solve);
}

void dsolve_pdeclairaut_init(void) {
    symtab_add_builtin("DSolve`PDEClairaut", builtin_dsolve_pdeclairaut);
    symtab_get_def("DSolve`PDEClairaut")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`PDEClairaut",
        "DSolve`PDEClairaut[eqn, u, {v1, v2}] solves the first-order PDE Clairaut "
        "form u == v1 u_{v1} + v2 u_{v2} + f(u_{v1}, u_{v2}); the complete integral "
        "is u == C[1] v1 + C[2] v2 + f(C[1], C[2]).  With IncludeSingularSolutions "
        "-> True the singular integral (the envelope of that family) is also "
        "returned.  Nonlinear in the derivatives, so the linear first-order PDE "
        "method declines to it.");
}
