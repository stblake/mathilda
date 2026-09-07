/*
 * dsolve_polyshift.c — DSolve`PolynomialShiftSubstitution.
 *
 * The x-dependent-shift generalisation of FirstOrderSubstitution.  Solves a
 * first-order ODE  y' = f(x, y)  in which y appears only inside a fractional
 * power of a term LINEAR in y with an x-dependent shift,
 *
 *     y' = R(x) + g(x) (phi(x) + c y)^p,     p a non-integer rational, c const,
 *
 * where R(x) = -phi'(x)/c so that the substitution  u = phi(x) + c y  removes
 * the y-dependence entirely:
 *
 *     u' = phi' + c y' = phi' + c R + c g u^p = c g(x) u^p     (separable).
 *
 * Its first integral is  u^(1-p)/(1-p) - Integrate[c g, x] == C[1], and with
 * u = phi + c y that is an implicit relation in (x, y) returned through
 * dsolve_run_implicit (verified by the implicit-function rule y' = -G_x/G_y).
 * The implicit form is BRANCH-SAFE: the radical's sign is never committed, so
 * the answer verifies where an explicit y = ... on one branch would not.
 *
 * These are Maple's [F(x), G(x)]-symmetry cases (invariant u = phi + c y).  The
 * heuristic Lie method (dsolve_lie.c `abaco2_similar`) targets the same class but
 * its Q = omega_y/omega_yy, T = Q_x/Q_y construction differentiates the radical
 * omega, blows past the node budget, and aborts; this deterministic substitution
 * sidesteps the symmetry machinery.  Runs after FirstOrderSubstitution and before
 * the LieSymmetry backstop (see dsolve.c).  Detection is confined to the pure-
 * power reduced form u' = k(x) u^p; a general separable u' = k(x) h(u) is future
 * work (e.g. 2.1.2-365, whose reduced h is 1 + I Sqrt[u]).
 *
 * Solves 2.1.2-402/-371/-372/-376/-378/-403/-424 (previously abaco2_similar aborts).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "integrate.h"
#include <stdlib.h>

/* Half = the Rational 1/2 (the Sqrt exponent). */
static Expr* ps_half(void) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ expr_new_integer(2), expr_new_integer(-1) }, 2));
}

/* Find, in `e`, a fractional-power atom base^p (p a non-integer rational, or a
 * Sqrt) whose `base` is LINEAR in Y with a constant (free of x and Y) nonzero
 * coefficient and actually contains Y.  On success returns base (owned) and sets
 * *pexp to the exponent (owned); else NULL.  First match in a pre-order walk. */
static Expr* ps_find_atom(const Expr* e, const char* Y, const char* xv, Expr** pexp) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const Expr* h = e->data.function.head;
    const Expr* base = NULL; Expr* p = NULL;
    if (h && h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        if (hn == SYM_Power && e->data.function.arg_count == 2) {
            const Expr* q = e->data.function.args[1];
            if (q->type == EXPR_FUNCTION && q->data.function.head &&
                q->data.function.head->type == EXPR_SYMBOL &&
                q->data.function.head->data.symbol.name == SYM_Rational) {
                base = e->data.function.args[0]; p = expr_copy((Expr*)q);   /* non-integer */
            }
        } else if (hn == SYM_Sqrt && e->data.function.arg_count == 1) {
            base = e->data.function.args[0]; p = ps_half();
        }
    }
    if (base && ds_contains(base, Y)) {
        Expr* c = ds_simplify(ds_d(expr_copy((Expr*)base), expr_new_symbol(Y)));  /* dbase/dY */
        bool lin = ds_free_of(c, Y) && ds_free_of(c, xv) && !ds_is_zero(c);
        expr_free(c);
        if (lin) { *pexp = p; return expr_copy((Expr*)base); }
    }
    if (p) expr_free(p);
    /* recurse */
    if (h) { Expr* r = ps_find_atom(h, Y, xv, pexp); if (r) return r; }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* r = ps_find_atom(e->data.function.args[i], Y, xv, pexp);
        if (r) return r;
    }
    return NULL;
}

Expr** dsolve_polyshift_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;

    const char* Y = intern_symbol("DSolve`psY");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Y)); /* consumes F */

    /* locate the fractional-power-of-(linear in Y) atom */
    Expr* pexp = NULL;
    Expr* base = ps_find_atom(FY, Y, xvar, &pexp);
    if (!base) { expr_free(FY); return NULL; }

    /* c = dbase/dY (const), phi = base - c Y (free of Y).  Plain evaluation, not
     * ds_simplify — phi may carry an x-radical on which Simplify can loop, and the
     * free-of / is-zero checks below are robust to an unsimplified form. */
    Expr* c   = ds_d(expr_copy(base), expr_new_symbol(Y));
    Expr* phi = eval_and_free(ds_call2(SYM_Subtract, expr_copy(base),
                    ds_call2(SYM_Times, expr_copy(c), expr_new_symbol(Y))));
    if (!ds_free_of(phi, Y)) { expr_free(FY); expr_free(base); expr_free(pexp); expr_free(c); expr_free(phi); return NULL; }

    /* u = base;  substitute Y -> (u - phi)/c  and form  G(x,u) = phi' + c F. */
    const char* U = intern_symbol("DSolve`psU");
    Expr* Yval = eval_and_free(ds_call2(SYM_Times,
                     ds_call2(SYM_Subtract, expr_new_symbol(U), expr_copy(phi)),
                     eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ expr_copy(c), expr_new_integer(-1) }, 2))));   /* (u-phi)/c */
    Expr* Fsub = ds_subst(expr_copy(FY), expr_new_symbol(Y), Yval);   /* consumes Yval */
    Expr* phip = ds_d(expr_copy(phi), expr_new_symbol(xvar));         /* phi' */
    /* Plain evaluation (NOT ds_simplify): the evaluator combines the same-base
     * powers and the additive cancellation (phi' + c R -> 0) that expose whether
     * the reduced form is the pure power k(x) u^p, WITHOUT Simplify's radical
     * rationalisation, which loops on shapes like (x + Sqrt[u]) u^(-1/2). */
    Expr* G = eval_and_free(ds_call2(SYM_Plus, phip,
                  ds_call2(SYM_Times, expr_copy(c), Fsub)));           /* consumes phip, Fsub */
    expr_free(FY);

    /* pure-power reduced form:  k(x) = G u^(-p) must be free of u (and of Y). */
    Expr* umP = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_new_symbol(U),
                        eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(pexp))) }, 2)); /* u^(-p) */
    Expr* k = eval_and_free(ds_call2(SYM_Times, expr_copy(G), umP));  /* consumes umP; eval, not Simplify */
    expr_free(G);
    bool ok = ds_free_of(k, U) && ds_free_of(k, Y) && !ds_is_zero(k);
    if (!ok) { expr_free(base); expr_free(pexp); expr_free(c); expr_free(phi); expr_free(k); return NULL; }

    /* 1 - p (nonzero for a radical); Xi = u^(1-p)/(1-p), u = phi + c y[x]. */
    Expr* omp = eval_and_free(ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(pexp)));
    if (ds_is_zero(omp)) { expr_free(base); expr_free(pexp); expr_free(c); expr_free(phi); expr_free(k); expr_free(omp); return NULL; }

    /* Kint = Integrate[k, x]  (must be elementary) */
    g_integrate_quiet++;
    Expr* Kint = ds_integrate(expr_copy(k), expr_new_symbol(xvar));
    g_integrate_quiet--;
    expr_free(k);
    if (ds_has_head(Kint, SYM_Integrate)) {
        expr_free(base); expr_free(pexp); expr_free(c); expr_free(phi); expr_free(omp); expr_free(Kint);
        return NULL;
    }

    /* base_y = base with Y -> y[x] ; Xi = base_y^(1-p) / (1-p) */
    Expr* base_y = ds_subst(expr_copy(base), expr_new_symbol(Y), ds_make_funcapp(yname, 0, xvar));
    Expr* Xi = eval_and_free(ds_call2(SYM_Times,
                   eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ base_y, expr_copy(omp) }, 2)),
                   eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(omp), expr_new_integer(-1) }, 2))));   /* base_y^(1-p)/(1-p) */

    /* G_impl = Xi - Kint;  dsolve_run_implicit forms G_impl == C[1].  Not
     * ds_simplify'd (Simplify chokes on radical-composed relations; the implicit
     * verifier differentiates it — cf. dsolve_chini.c). */
    Expr* Gimpl = eval_and_free(ds_call2(SYM_Subtract, Xi, Kint));   /* consumes Xi, Kint */

    expr_free(base); expr_free(pexp); expr_free(c); expr_free(phi); expr_free(omp);

    Expr** out = malloc(sizeof(Expr*));
    out[0] = Gimpl;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_polyshift(Expr* res) {
    return dsolve_method_builtin_implicit(res, dsolve_polyshift_try);
}

void dsolve_polyshift_init(void) {
    symtab_add_builtin("DSolve`PolynomialShiftSubstitution", builtin_dsolve_polyshift);
    symtab_get_def("DSolve`PolynomialShiftSubstitution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`PolynomialShiftSubstitution",
        "DSolve`PolynomialShiftSubstitution[eqn, y, x] solves y' == R(x) + "
        "g(x) (phi(x) + c y)^p (p a non-integer rational, c constant, "
        "R = -phi'/c) by the substitution u = phi(x) + c y, which reduces it to the "
        "separable u' == c g(x) u^p; the first integral u^(1-p)/(1-p) - "
        "Integrate[c g, x] == C[1] is returned implicitly (branch-safe). The "
        "x-dependent-shift generalisation of FirstOrderSubstitution.");
}
