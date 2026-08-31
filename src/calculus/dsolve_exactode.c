/*
 * dsolve_exactode.c — DSolve`ExactODE.
 *
 * Solves a linear ODE of order n >= 2 whose left side is a total derivative,
 *
 *     L[y] = Sum_{k=0}^{n} a_k(x) y^(k) == g(x),     L[y] == d/dx( M[y] ),
 *
 * where M is an order-(n-1) linear operator.  Such an equation is EXACT and
 * integrates once to the first integral
 *
 *     M[y] = Sum_{j=0}^{n-1} b_j(x) y^(j) == Integrate[g, x] + C[n],
 *
 * a lower-order ODE handed back to the scalar cascade (recurse into DSolve).
 * The first-integral coefficients follow from matching
 * Sum a_k y^(k) = d/dx( Sum b_j y^(j) ):
 *
 *     b_{n-1} = a_n,   b_{k-1} = a_k - b_k'   (k = n-1 .. 1),
 *
 * and the exactness condition is the leftover a_0 == b_0'  (equivalently
 * Sum (-1)^k a_k^(k) == 0).  The reduced equation, being order n-1, contributes
 * C[1..n-1]; the integration constant is C[n] — contiguous, no renumbering
 * (mirroring dsolve_reduce_order.c's ds_const(2) after a first-order C[1]).
 * Iterated exactness needs no special handling: the recursion re-enters ExactODE
 * on the reduced equation, so a doubly-exact 3rd-order equation reduces twice.
 *
 * Scope of this first cut: linear, order >= 2, genuinely exact.  First-order
 * exact M + N y' == 0 is DSolve`Exact.  Integrating-factor (adjoint) exactness
 * and nonlinear total-derivative detection are future work.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

/* Extract the RHS of {{y[x] -> expr}} (applied form) from a DSolve result. */
static Expr* extract_applied(Expr* r, const char* yfun) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL
                && lhs->data.function.head->data.symbol.name == yfun)
                return expr_copy(rule->data.function.args[1]);
        }
    }
    return NULL;
}

Expr** dsolve_exactode_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr** c; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return NULL;   /* nonlinear -> decline */
    if (n < 2) {                                             /* first-order is DSolve`Exact */
        for (int k = 0; k <= n; k++) expr_free(c[k]);
        free(c); expr_free(g);
        return NULL;
    }

    /* First-integral coefficients b[0..n-1]:  b[n-1]=a_n, b[k-1]=a_k - b[k]'. */
    Expr** b = malloc((size_t)n * sizeof(Expr*));
    b[n - 1] = expr_copy(c[n]);
    for (int k = n - 1; k >= 1; k--) {
        Expr* db = ds_d(expr_copy(b[k]), expr_new_symbol(xvar));
        b[k - 1] = eval_and_free(ds_call2(SYM_Subtract, expr_copy(c[k]), db));
    }

    /* Exactness: a_0 - b_0' == 0 (Together first to help zero_test on rationals). */
    Expr* db0 = ds_d(expr_copy(b[0]), expr_new_symbol(xvar));
    Expr* chk = eval_and_free(ds_call2(SYM_Subtract, expr_copy(c[0]), db0));
    Expr* chk_t = eval_and_free(ds_call1("Together", expr_copy(chk)));
    bool exact = ds_is_zero(chk_t);
    expr_free(chk); expr_free(chk_t);

    Expr* body = NULL;
    if (exact) {
        /* Reduced forcing: Integrate[g, x] + C[n].  Decline on a non-elementary
         * antiderivative (an unevaluated Integrate would leak into the sub-solve). */
        Expr* Gint = ds_integrate(expr_copy(g), expr_new_symbol(xvar));
        if (!ds_has_head(Gint, SYM_Integrate)) {
            Expr* rhs = eval_and_free(ds_call2(SYM_Plus, Gint, ds_const(n)));

            /* M[y] = Sum_{j=0}^{n-1} b[j] * y^(j). */
            Expr** terms = malloc((size_t)n * sizeof(Expr*));
            for (int j = 0; j < n; j++)
                terms[j] = ds_call2(SYM_Times, expr_copy(b[j]), ds_make_funcapp(yname, j, xvar));
            Expr* M = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)n));
            free(terms);

            Expr* reduced = expr_new_function(expr_new_symbol(SYM_Equal),
                                (Expr*[]){ M, rhs }, 2);

            /* Recurse into the scalar cascade on the order-(n-1) equation. */
            Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                             (Expr*[]){ reduced, ds_make_funcapp(yname, 0, xvar),
                                        expr_new_symbol(xvar) }, 3);
            Expr* r = eval_and_free(call);
            body = extract_applied(r, yname);
            expr_free(r);
        } else {
            expr_free(Gint);
        }
    }

    for (int k = 0; k <= n; k++) expr_free(c[k]);
    free(c); expr_free(g);
    for (int j = 0; j < n; j++) expr_free(b[j]);
    free(b);

    if (!body) return NULL;
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_exact_ode(Expr* res) {
    return dsolve_method_builtin(res, dsolve_exactode_try);
}

void dsolve_exactode_init(void) {
    symtab_add_builtin("DSolve`ExactODE", builtin_dsolve_exact_ode);
    symtab_get_def("DSolve`ExactODE")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ExactODE",
        "DSolve`ExactODE[eqn, y, x] solves a linear ODE of order >= 2 whose left "
        "side is a total derivative d/dx(M[y]) (an exact equation): it integrates "
        "once to the first integral M[y] == Integrate[g, x] + C[n] and recurses "
        "into the scalar solver.");
}
