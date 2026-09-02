/*
 * dsolve_factorable.c — DSolve`Factorable.
 *
 * Factors the ODE residual as a polynomial in the dependent function and its
 * derivatives.  When it splits into two or more genuine factors F1 * F2 * ... == 0,
 * each factor Fi == 0 is a lower-complexity ODE recursed back into the scalar
 * cascade, and the branch bodies are unioned (a branch that solves any factor
 * satisfies the product).  Mirrors SymPy's `factorable`.
 *
 * The derivatives y[x], y'[x], ... are first replaced by plain symbols d0, d1, ...
 * so FactorList sees an ordinary polynomial: Factor/FactorList over raw Derivative
 * funcapps both HANGS on a non-polynomial argument (e.g. a Sqrt[y[x]] coefficient)
 * AND misfactors even polynomial funcapp forms.  A PolynomialQ gate on that
 * substituted form is therefore the anti-hang guard — a non-polynomial residual
 * (a radical/transcendental in the funcapps) declines cheaply before FactorList.
 *
 * Runs at the FRONT of the cascade (before Quadrature), so a product form is split
 * before the specialists match the whole thing; an irreducible residual (a single
 * dependent-function factor) declines.  Every branch is back-substitution-verified.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>
#include <stdio.h>

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

Expr** dsolve_factorable_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    int n = P->max_order[0];
    if (n < 1) return NULL;
    Expr* R = P->eq_residuals[0];

    /* interned plain symbols d0..dn for the funcapps y^(k)[x] */
    const char** dsym = malloc((size_t)(n + 1) * sizeof(char*));
    for (int k = 0; k <= n; k++) {
        char nm[32];
        snprintf(nm, sizeof(nm), "DSolve`fD%d", k);
        dsym[k] = intern_symbol(nm);
    }

    /* Rsub = R with y^(k)[x] -> d_k */
    Expr* Rsub = expr_copy(R);
    for (int k = 0; k <= n; k++)
        Rsub = ds_subst(Rsub, ds_make_funcapp(yname, k, xvar), expr_new_symbol(dsym[k]));

    /* PolynomialQ[Rsub, {d0..dn}] — anti-hang / anti-misfactor gate */
    Expr** sl = malloc((size_t)(n + 1) * sizeof(Expr*));
    for (int k = 0; k <= n; k++) sl[k] = expr_new_symbol(dsym[k]);
    Expr* symlist = expr_new_function(expr_new_symbol(SYM_List), sl, (size_t)(n + 1));
    free(sl);
    Expr* pq = eval_and_free(ds_call2("PolynomialQ", expr_copy(Rsub), expr_copy(symlist)));
    bool ispoly = (pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True);
    expr_free(pq);
    if (!ispoly) { expr_free(Rsub); expr_free(symlist); free(dsym); return NULL; }

    Expr* fl = eval_and_free(ds_call1("FactorList", Rsub));   /* consumes Rsub */
    expr_free(symlist);
    if (!head_is(fl, SYM_List)) { expr_free(fl); free(dsym); return NULL; }

    /* keep the genuine DIFFERENTIAL factors: a factor must contain a derivative
     * d_k with k >= 1.  A factor in d_0 only (the bare function, e.g. y[x]) is an
     * algebraic constraint giving the trivial y == 0 branch, not a differential
     * factor — and, worse, factoring it out injects a spurious first branch that
     * breaks the recursive callers (AutonomousReduction / ReductionOfOrder) which
     * take the first sub-solution.  Substitute the symbols back to funcapps. */
    size_t m = fl->data.function.arg_count;
    Expr** facs = malloc((m ? m : 1) * sizeof(Expr*));
    size_t nf = 0;
    for (size_t i = 0; i < m; i++) {
        Expr* entry = fl->data.function.args[i];
        if (!head_is(entry, SYM_List) || entry->data.function.arg_count < 1) continue;
        Expr* factor = entry->data.function.args[0];
        bool has_deriv = false;
        for (int k = 1; k <= n && !has_deriv; k++)
            if (ds_contains(factor, dsym[k])) has_deriv = true;
        if (!has_deriv) continue;
        Expr* f = expr_copy(factor);
        for (int k = 0; k <= n; k++)
            f = ds_subst(f, expr_new_symbol(dsym[k]), ds_make_funcapp(yname, k, xvar));
        facs[nf++] = f;
    }
    expr_free(fl);
    free(dsym);

    /* genuinely factored only when >= 2 dependent-function factors */
    if (nf < 2) { for (size_t i = 0; i < nf; i++) expr_free(facs[i]); free(facs); return NULL; }

    Expr** acc = NULL; size_t nacc = 0;
    for (size_t i = 0; i < nf; i++) {
        Expr* subeqn = expr_new_function(expr_new_symbol(SYM_Equal),
                           (Expr*[]){ facs[i], expr_new_integer(0) }, 2);  /* consumes facs[i] */
        recurse_collect(subeqn, yname, xvar, &acc, &nacc);
    }
    free(facs);

    if (nacc == 0) { if (acc) free(acc); return NULL; }
    *nbranch = nacc;
    return acc;
}

static Expr* builtin_dsolve_factorable(Expr* res) {
    return dsolve_method_builtin(res, dsolve_factorable_try);
}

void dsolve_factorable_init(void) {
    symtab_add_builtin("DSolve`Factorable", builtin_dsolve_factorable);
    symtab_get_def("DSolve`Factorable")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Factorable",
        "DSolve`Factorable[eqn, y, x] factors the equation as a product in the "
        "derivatives (F1 F2 ... == 0), solves each factor Fi == 0 by recursing the "
        "cascade, and unions the solution branches.");
}
