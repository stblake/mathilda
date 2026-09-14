/*
 * dsolve_ratsol2.c — DSolve`PolynomialSolution.
 *
 * A homogeneous second-order linear ODE with rational coefficients,
 *     y''[x] + P(x) y'[x] + Q(x) y[x] == 0,   P, Q rational in x,
 * whose FULL fundamental set is polynomial (both independent solutions are
 * polynomials in x) is solved directly by an undetermined-coefficient search:
 * substitute the degree-bounded ansatz  y = a_0 + a_1 x + ... + a_D x^D, clear the
 * common denominator, and solve the linear system Coefficient[., x, j] == 0 for
 * the a_k.  The null space of that system is exactly the space of polynomial
 * solutions; a 2-dimensional null space gives the clean fundamental set
 * {y1, y2} and the general solution C[1] y1 + C[2] y2.
 *
 * Why this runs BEFORE Kovacic.  For e.g.  (x^2+1) y'' - 2x y' + 2 y == 0  the
 * fundamental set is the clean {x, x^2-1}, but Kovacic returns it wrapped in a
 * complex-radical basis  ( x Sqrt[x-I] Sqrt[x+I] / Sqrt[x^2+1], ... )  that does
 * NOT reduce to the rational form, so a downstream VariationOfParameters for a
 * forcing term cannot close its integrals.  This is precisely the residue that
 * left Nasser Abbasi 2.2.26-2592 and 2.2.25-2410 unsolved: a genuinely
 * elementary answer blocked only by the shape of the homogeneous basis.
 * Producing the clean polynomial basis here lets the nonhomogeneous VoP backstop
 * (dsolve_nonhomog_vop, which recurses DSolve on the homogeneous part) close.
 *
 * Homogeneous only: an inhomogeneous equation declines here (dsolve_second_order_PQ
 * requires zero forcing) and is handled by dsolve_nonhomog_vop, which recurses
 * into this method for the basis.  The method never recurses DSolve itself, so no
 * new recursion path is introduced.  A returned basis is numerically self-verified
 * (ds_branch_num_ok) before it is accepted.
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

/* Degree cap of the polynomial ansatz.  A polynomial solution of degree > CAP is
 * out of scope here (such an ODE is typically an orthogonal-polynomial /
 * special-function family already claimed by SpecialFunctionForm earlier in the
 * cascade); a too-small cap can only make this method DECLINE (fall through to
 * Kovacic), never return a wrong answer, since Solve on the exact coefficient
 * equations yields only genuine polynomial solutions.  Kept modest so the common
 * DECLINE path (no polynomial solution) stays cheap: this method runs on every
 * 2nd-order rational-coefficient equation reaching it, and an oversized ansatz adds
 * measurable per-equation cost to cases another method (Frobenius/Kovacic) then owns. */
#define RATSOL2_DEG_CAP 6

/* C[k] integration constant. */
static Expr* Ck(int k) {
    return expr_new_function(expr_new_symbol(intern_symbol("C")),
                             (Expr*[]){ expr_new_integer(k) }, 1);
}

/* Free the coefficient vector from dsolve_linear_coeffs (length n+1) and the array. */
static void free_coeffs(Expr** c, int n) {
    for (int k = 0; k <= n; k++) expr_free(c[k]);
    free(c);
}

Expr** dsolve_ratsol2_try(DSolveProblem* P, size_t* nbranch) {
    /* Second-order linear ODE (homogeneous OR forced): c[2] y'' + c[1] y' + c[0] y == g. */
    Expr** c = NULL; Expr* g = NULL; int n = 0;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return NULL;    /* nonlinear -> decline */
    if (n != 2) { free_coeffs(c, n); expr_free(g); return NULL; }
    const char* xvar = P->ind_names[0];

    /* Rational-coefficient gate: the ansatz + Coefficient machinery is a
     * polynomial/rational-function transform in x.  A transcendental coefficient
     * (Sin[x], E^(x^2), ...) is not this class and would make Together/Coefficient
     * churn; decline fast and let Kovacic / ChangeOfVariable / Frobenius own it. */
    if (!ds_is_rational_in(c[0], xvar) || !ds_is_rational_in(c[1], xvar)
        || !ds_is_rational_in(c[2], xvar)) {
        free_coeffs(c, n); expr_free(g); return NULL;
    }

    const int dg = RATSOL2_DEG_CAP;

    /* Ansatz Y = sum_{k=0}^{dg} a_k x^k, a_k = DSolve`a<k> (fresh internal symbols). */
    const char** avar = malloc((size_t)(dg + 1) * sizeof(char*));
    Expr** aterms = malloc((size_t)(dg + 1) * sizeof(Expr*));
    for (int k = 0; k <= dg; k++) {
        char buf[32]; snprintf(buf, sizeof(buf), "DSolve`a%d", k);
        avar[k] = intern_symbol(buf);
        aterms[k] = ds_call2(SYM_Times, expr_new_symbol(avar[k]),
                             ds_call2(SYM_Power, expr_new_symbol(xvar), expr_new_integer(k)));
    }
    Expr* Y = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                              aterms, (size_t)(dg + 1)));
    free(aterms);

    /* Homogeneous operator L0 = c2 Y'' + c1 Y' + c0 Y ; clear the denominator
     * (a function of x only, free of the a_k) so the coefficient equations are
     * polynomial in x. */
    Expr* Ypp = ds_d(ds_d(expr_copy(Y), expr_new_symbol(xvar)), expr_new_symbol(xvar));
    Expr* Yp  = ds_d(expr_copy(Y), expr_new_symbol(xvar));
    Expr* L = eval_and_free(ds_call2(SYM_Plus,
                 eval_and_free(ds_call2(SYM_Plus,
                     ds_call2(SYM_Times, expr_copy(c[2]), Ypp),
                     ds_call2(SYM_Times, expr_copy(c[1]), Yp))),
                 ds_call2(SYM_Times, expr_copy(c[0]), expr_copy(Y))));
    /* Keep only the leading coefficient (VoP) and forcing g past this point. */
    Expr* leadcoef = expr_copy(c[2]);
    free_coeffs(c, n);
    Expr* Lnum = eval_and_free(ds_call1("Expand",
                     eval_and_free(ds_call1("Numerator",
                         eval_and_free(ds_call1("Together", L))))));

    /* Coefficient equations: CoefficientList[Lnum, x] returns every x-power coefficient
     * (each linear in the a_k) in ONE pass — far cheaper on the common DECLINE path than
     * a Coefficient[...,x,j] loop, and Exponent is unusable here anyway (it returns 0 on
     * an expression carrying a_k funcapps).  Drop the trivially-zero rows. */
    Expr* clist = eval_and_free(expr_new_function(expr_new_symbol(intern_symbol("CoefficientList")),
                      (Expr*[]){ Lnum, expr_new_symbol(xvar) }, 2));   /* consumes Lnum */
    if (!head_is(clist, SYM_List)) {
        expr_free(clist); free(avar); expr_free(Y); expr_free(leadcoef); expr_free(g); return NULL;
    }
    size_t ncoef = clist->data.function.arg_count;
    Expr** eqs = malloc((ncoef ? ncoef : 1) * sizeof(Expr*));
    size_t neq = 0;
    for (size_t j = 0; j < ncoef; j++) {
        Expr* cj = clist->data.function.args[j];
        if (ds_is_zero(cj)) continue;
        eqs[neq++] = eval_and_free(ds_call2(SYM_Equal, expr_copy(cj), expr_new_integer(0)));
    }
    expr_free(clist);
    if (neq == 0) { free(eqs); free(avar); expr_free(Y); expr_free(leadcoef); expr_free(g); return NULL; }

    Expr* eqlist = eval_and_free(expr_new_function(expr_new_symbol(SYM_List), eqs, neq));
    free(eqs);
    Expr** vs = malloc((size_t)(dg + 1) * sizeof(Expr*));
    for (int k = 0; k <= dg; k++) vs[k] = expr_new_symbol(avar[k]);
    Expr* varlist = eval_and_free(expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)(dg + 1)));
    free(vs);

    Expr* sol = eval_and_free(ds_call2("Solve", eqlist, varlist));   /* consumes both */
    if (!head_is(sol, SYM_List) || sol->data.function.arg_count < 1) {
        expr_free(sol); free(avar); expr_free(Y); expr_free(leadcoef); expr_free(g); return NULL;
    }
    /* gensol = Y /. (first solution rule set); free a_k parametrise the null space. */
    Expr* rules = expr_copy(sol->data.function.args[0]);
    expr_free(sol);
    Expr* gensol = eval_and_free(ds_call2("ReplaceAll", expr_copy(Y), rules));
    expr_free(Y);

    /* Basis: y_k = D[gensol, a_k].  gensol is linear in the surviving free a_k, so
     * each nonzero derivative is one independent polynomial solution and the
     * constrained a_k give 0. */
    Expr* basis[2] = { NULL, NULL };
    size_t dim = 0;
    for (int k = 0; k <= dg; k++) {
        Expr* bk = eval_and_free(ds_call1("Expand", ds_d(expr_copy(gensol), expr_new_symbol(avar[k]))));
        if (ds_is_zero(bk)) { expr_free(bk); continue; }
        if (dim < 2) basis[dim] = bk; else expr_free(bk);
        dim++;
    }
    expr_free(gensol);
    free(avar);

    /* Need a full, genuinely 2-dimensional polynomial fundamental set.  (dim > 2 is
     * impossible for a 2nd-order ODE; guard anyway.) */
    if (dim != 2 || !basis[0] || !basis[1]) {
        if (basis[0]) expr_free(basis[0]);
        if (basis[1]) expr_free(basis[1]);
        expr_free(leadcoef); expr_free(g); return NULL;
    }
    /* Independence: Wronskian y1 y2' - y2 y1' not identically zero. */
    Expr* w = eval_and_free(ds_call2(SYM_Subtract,
                 ds_call2(SYM_Times, expr_copy(basis[0]), ds_d(expr_copy(basis[1]), expr_new_symbol(xvar))),
                 ds_call2(SYM_Times, expr_copy(basis[1]), ds_d(expr_copy(basis[0]), expr_new_symbol(xvar)))));
    bool indep = !ds_is_zero(w);
    expr_free(w);
    if (!indep) { expr_free(basis[0]); expr_free(basis[1]); expr_free(leadcoef); expr_free(g); return NULL; }

    /* Homogeneous general solution C[1] y1 + C[2] y2. */
    Expr* hom = eval_and_free(ds_call2(SYM_Plus,
                    ds_call2(SYM_Times, Ck(1), expr_copy(basis[0])),
                    ds_call2(SYM_Times, Ck(2), expr_copy(basis[1]))));

    Expr* body;
    if (ds_is_zero(g)) {
        expr_free(basis[0]); expr_free(basis[1]); expr_free(leadcoef); expr_free(g);
        body = hom;
    } else {
        /* Nonhomogeneous: add the variation-of-parameters particular over the clean
         * polynomial basis (leadcoef = coefficient of y'').  This is exactly the
         * closure that a complex-radical Kovacic basis blocks (2.2.26-2592,
         * 2.2.25-2410).  Decline if the VoP integral is not elementary. */
        Expr* yp = dsolve_variation_of_parameters(basis, 2, g, leadcoef, xvar);
        expr_free(basis[0]); expr_free(basis[1]); expr_free(leadcoef); expr_free(g);
        if (!yp || ds_has_head(yp, SYM_Integrate)) {
            if (yp) expr_free(yp);
            expr_free(hom); return NULL;
        }
        body = eval_and_free(ds_call2(SYM_Plus, hom, yp));
        if (ds_has_head(body, SYM_Integrate)) { expr_free(body); return NULL; }
    }

    /* Numerically self-verify before accepting (defence: never return a solution
     * that does not satisfy the ODE). */
    if (!ds_branch_num_ok(P, body)) { expr_free(body); return NULL; }

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_ratsol2(Expr* res) {
    return dsolve_method_builtin(res, dsolve_ratsol2_try);
}

void dsolve_ratsol2_init(void) {
    symtab_add_builtin("DSolve`PolynomialSolution", builtin_dsolve_ratsol2);
    symtab_get_def("DSolve`PolynomialSolution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`PolynomialSolution",
        "DSolve`PolynomialSolution[eqn, y, x] solves a second-order linear ODE "
        "c2(x)y''+c1(x)y'+c0(x)y==g(x) with rational coefficients whose homogeneous "
        "fundamental set is polynomial, via a degree-bounded undetermined-coefficient "
        "search plus (when forced) variation of parameters over that clean basis.");
}
