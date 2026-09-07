/*
 * dsolve_nth_algebraic.c — DSolve`NthAlgebraic.
 *
 * Handles an ODE that is *algebraic* (polynomial of degree >= 2) in its highest
 * derivative y^(n): Solve[eqn, y^(n)] gives several branches y^(n) == g_k, each
 * recursed back into the scalar cascade (a branch free of y hits Quadrature) and
 * the branch bodies unioned.  Also the degenerate no-derivative case (an equation
 * purely algebraic in y[x]): Solve[eqn, y[x]] gives the constant/algebraic
 * solutions directly.  Mirrors SymPy's `nth_algebraic`.
 *
 * Runs at the FRONT of the cascade (with Factorable), so a top-derivative-in-a-power
 * form is split before the linear specialists try to match it.  Declines the normal
 * case (degree 1 in the top derivative with derivatives present) so those methods
 * claim it.  The substrate back-substitution verifies every unioned branch.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* Recurse DSolve[subeqn, y[x], x] (applied form) and append every explicit body
 * to *acc (grown by realloc).  `subeqn` is consumed. */
static void recurse_collect(Expr* subeqn, const char* yname, const char* xvar,
                            Expr*** acc, size_t* nacc) {
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ subeqn, ds_make_funcapp(yname, 0, xvar),
                                expr_new_symbol(xvar) }, 3);
    Expr* r = eval_and_free(call);
    size_t nb = 0;
    Expr** bodies = dsolve_extract_applied_bodies(r, yname, &nb);
    expr_free(r);
    if (bodies) {
        *acc = realloc(*acc, (*nacc + nb) * sizeof(Expr*));
        for (size_t i = 0; i < nb; i++) (*acc)[(*nacc)++] = bodies[i];
        free(bodies);
    }
}

Expr** dsolve_nth_algebraic_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    int n = P->max_order[0];
    Expr* R = P->eq_residuals[0];

    Expr** acc = NULL; size_t nacc = 0;

    if (n < 1) {
        /* degenerate: purely algebraic in y[x] — Solve[R==0, y[x]] */
        const char* Yn = intern_symbol("DSolve`naY");
        Expr* Rsub = ds_subst(expr_copy(R), ds_make_funcapp(yname, 0, xvar),
                              expr_new_symbol(Yn));
        Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                       (Expr*[]){ Rsub, expr_new_integer(0) }, 2);
        Expr* sol = ds_solve(eq, expr_new_symbol(Yn));
        size_t ng = 0;
        Expr** gs = dsolve_extract_solutions(sol, Yn, &ng);   /* roots are the bodies */
        if (sol) expr_free(sol);
        if (gs) { acc = gs; nacc = ng; }
    } else {
        /* substitute the top derivative y^(n)[x] -> plain symbol, measure degree */
        const char* Dn = intern_symbol("DSolve`naD");
        Expr* topLit = ds_make_funcapp(yname, n, xvar);
        Expr* Rsub = ds_subst(expr_copy(R), expr_copy(topLit), expr_new_symbol(Dn));
        Expr* degE = eval_and_free(ds_call2("Exponent", expr_copy(Rsub),
                                            expr_new_symbol(Dn)));
        long deg = (degE->type == EXPR_INTEGER) ? degE->data.integer : -1;
        expr_free(degE);

        /* only solve when genuinely polynomial in the top derivative (guards
         * Solve against a radical/transcendental form the specialists own) */
        Expr* pq = eval_and_free(ds_call2("PolynomialQ", expr_copy(Rsub),
                                          expr_new_symbol(Dn)));
        bool ispoly = (pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True);
        expr_free(pq);

        if (deg >= 2 && ispoly) {
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                           (Expr*[]){ Rsub, expr_new_integer(0) }, 2);  /* consumes Rsub */
            Expr* sol = ds_solve(eq, expr_new_symbol(Dn));
            size_t ng = 0;
            Expr** gs = dsolve_extract_solutions(sol, Dn, &ng);
            if (sol) expr_free(sol);
            /* A branch given as an implicit Root[...] object means Solve found no
             * explicit closed form for the top derivative (e.g. a general cubic in
             * y' with symbolic coefficients).  Recursing DSolve on y^(n) == Root[...]
             * is futile and can HANG (the implicit-derivative ODE has no method that
             * closes it, and some rewrite loops), so decline the whole method and let
             * a later specialist (Lagrange for a Lagrange/d'Alembert form, ...) take
             * the original equation. */
            bool has_root = false;
            const char* rootn = intern_symbol("Root");
            if (gs) for (size_t i = 0; i < ng; i++)
                if (ds_has_head(gs[i], rootn)) { has_root = true; break; }
            if (gs && !has_root) {
                for (size_t i = 0; i < ng; i++) {
                    /* y^(n)[x] == g_i */
                    Expr* subeqn = expr_new_function(expr_new_symbol(SYM_Equal),
                                       (Expr*[]){ expr_copy(topLit), gs[i] }, 2);
                    recurse_collect(subeqn, yname, xvar, &acc, &nacc);
                }
            } else if (gs) {
                for (size_t i = 0; i < ng; i++) expr_free(gs[i]);   /* declined */
            }
            if (gs) free(gs);
        } else if (!ispoly) {
            /* The top derivative appears RATIONALLY (in a denominator) rather than
             * polynomially, e.g. A/y' == B.  Clear it: Numerator[Together[A/Dn - B]]
             * = A - B Dn is polynomial in Dn.  When that cleared numerator genuinely
             * carries the top derivative (degree >= 1 in Dn), recurse DSolve on the
             * cleared ODE (top derivative restored) and let the full cascade own it --
             * the linear normalizer handles a y-dependent leading coefficient (y^3 x'
             * = 1 - 4 y^2 x) and Lagrange/NthAlgebraic-poly a higher-degree cleared
             * form.  This closes the 12000.org "x = x(y)" spelling (A/x'[y] == B).
             * The cleared form has the derivative polynomially, so re-entry can never
             * clear again (no infinite recursion), and a normal polynomial-in-y' ODE
             * (ispoly true) never reaches here, so nothing is hijacked. */
            Expr* numDn = eval_and_free(ds_call1("Numerator",
                              eval_and_free(ds_call1("Together", expr_copy(Rsub)))));
            Expr* npq = eval_and_free(ds_call2("PolynomialQ", expr_copy(numDn),
                                               expr_new_symbol(Dn)));
            bool npoly = (npq->type == EXPR_SYMBOL && npq->data.symbol.name == SYM_True);
            expr_free(npq);
            /* The cleared numerator must still carry the top derivative (else clearing
             * produced no ODE).  Tested with FreeQ, NOT Exponent: Exponent[expr, Dn]
             * mis-returns 0 whenever expr contains any function application (Sin[y],
             * x[y], ...), which every "x = x(y)" residual does. */
            Expr* fq = eval_and_free(ds_call2("FreeQ", expr_copy(numDn),
                                              expr_new_symbol(Dn)));
            bool has_deriv = (fq->type == EXPR_SYMBOL && fq->data.symbol.name == SYM_False);
            expr_free(fq);
            if (npoly && has_deriv) {
                /* restore the top derivative (Dn -> y^(n)[x]) and recurse on cleared==0 */
                Expr* clearedR = ds_subst(numDn, expr_new_symbol(Dn), expr_copy(topLit));
                Expr* subeqn = expr_new_function(expr_new_symbol(SYM_Equal),
                                   (Expr*[]){ clearedR, expr_new_integer(0) }, 2);
                recurse_collect(subeqn, yname, xvar, &acc, &nacc);
            } else {
                expr_free(numDn);
            }
            expr_free(Rsub);
        } else {
            expr_free(Rsub);
        }
        expr_free(topLit);
    }

    if (nacc == 0) { if (acc) free(acc); return NULL; }
    *nbranch = nacc;
    return acc;
}

static Expr* builtin_dsolve_nth_algebraic(Expr* res) {
    return dsolve_method_builtin(res, dsolve_nth_algebraic_try);
}

void dsolve_nth_algebraic_init(void) {
    symtab_add_builtin("DSolve`NthAlgebraic", builtin_dsolve_nth_algebraic);
    symtab_get_def("DSolve`NthAlgebraic")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`NthAlgebraic",
        "DSolve`NthAlgebraic[eqn, y, x] solves an ODE algebraic (degree >= 2) in its "
        "highest derivative by solving for that derivative and recursing on each "
        "branch; also the degenerate purely-algebraic case.");
}
