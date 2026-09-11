/*
 * dsolve_abel_air.c — DSolve`AbelAIR (Abel-2nd-kind / rational-in-y scaling reduction).
 *
 * A down payment on the full Abel Invariant Rational method (DSOLVE_PLAN.md M13,
 * previously deferred).  It closes the Abel-2nd-kind and rational-in-y first-order
 * ODEs whose right-hand side, after Together, is
 *
 *     y' == omega(x, y) = N(x, y) / D(x, y),
 *
 * with N a polynomial of degree <= 2 in y and D = c(x) * (P1(x) y + P0(x))^k
 * (k in {1, 2}) — a SINGLE linear y-factor in the denominator.  These forms are
 * exactly the ones the polynomial-RHS specialists (Riccati / Chini / Abel-1st) reject
 * (they require D free of y), which is why they otherwise fall to the heuristic Lie
 * backstop and time out.  The key structural fact is that the scaling substitution
 *
 *     u = y / s(x),   s(x) = -P0/P1  (the y-root of the denominator),
 *
 * renders the equation SEPARABLE in (x, u):  u' == A(x) B(u).  The first integral
 *     G(x, y) = Integrate[1/B, u]|_{u -> y/s}  -  Integrate[A, x]   (relation G == C[1])
 * is then returned through dsolve_run_implicit, whose implicit-function-rule verify
 * (y' == -G_x/G_y) confirms it; the separation itself is verified exactly
 * (omega_u - A B == 0) before integrating, so a non-separable shape declines rather
 * than shipping a wrong integral.
 *
 * Closes e.g. §2.2.17-1604 (u=e^x y -> u'=2x/(1+u)), -1607 (u=e^{-2x} y -> u'=x/(1-u)),
 * -1606 ((e^x+y)^2 denominator), -1676.  The full Abel differential-invariant
 * hierarchy (constant / rational invariant solvable-class table) remains deferred.
 *
 * Runs after Chini/Abel-1st and before the substitution reductions and the Lie
 * backstop (dsolve.c).  No recursion, so no wall-clock/memo bounding is needed.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "integrate.h"          /* g_integrate_quiet */
#include <stdlib.h>
#include <math.h>

/* head[a]/head[a,b]/head[a,b,c] evaluated; args consumed. */
static Expr* c1(const char* h, Expr* a) { return eval_and_free(ds_call1(h, a)); }
static Expr* c2(const char* h, Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2));
}
static Expr* c3(const char* h, Expr* a, Expr* b, Expr* c) {
    return eval_and_free(expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b, c }, 3));
}
static Expr* powneg1(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ base, expr_new_integer(-1) }, 2);
}
/* Integer value of an EXPR_INTEGER, else -1 (used to read Exponent[...] results). */
static long as_int(const Expr* e) { return (e && e->type == EXPR_INTEGER) ? (long)e->data.integer : -1; }

/* True iff PolynomialQ[e, {var}] returns True. */
static bool poly_in(const Expr* e, const char* var) {
    Expr* q = c2("PolynomialQ", expr_copy((Expr*)e),
                 expr_new_function(expr_new_symbol(SYM_List),
                                   (Expr*[]){ expr_new_symbol(var) }, 1));
    bool r = (q && q->type == EXPR_SYMBOL && q->data.symbol.name == intern_symbol("True"));
    expr_free(q);
    return r;
}

/* Numeric pre-filter for a separability residual chk = omega_u - A*B (in x and U):
 * false only when it is clearly non-zero at a generic real sample, so the expensive
 * symbolic ds_is_zero is skipped on a genuinely non-separable shape (mirror of
 * dsolve_separable.c's sep_chk_maybe_zero). */
static bool air_split_maybe_zero(const Expr* chk, const char* xv, const char* Uv) {
    static const double xs[] = { 1.3, 0.7 }, us[] = { 0.9, 1.7 };
    for (int i = 0; i < 2; i++) {
        Expr* e = expr_copy((Expr*)chk);
        e = ds_subst(e, expr_new_symbol(xv), expr_new_real(xs[i]));
        e = ds_subst(e, expr_new_symbol(Uv), expr_new_real(us[i]));
        e = c1("Abs", c1("N", e));
        double m = (e && e->type == EXPR_REAL) ? e->data.real
                 : (e && e->type == EXPR_INTEGER) ? (double)e->data.integer : NAN;
        expr_free(e);
        if (!isnan(m) && isfinite(m) && m > 1e-6) return false;
    }
    return true;
}

Expr** dsolve_abel_air_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`airY");
    const char* Un = intern_symbol("DSolve`airU");

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;
    Expr* omega = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    /* Together -> N / D; both must be polynomial in Y. */
    Expr* T = c1("Together", expr_copy(omega));
    Expr* N = c1("Numerator", expr_copy(T));
    Expr* D = c1("Denominator", T);
    if (!poly_in(N, Yn) || !poly_in(D, Yn)) {
        expr_free(omega); expr_free(N); expr_free(D); return NULL;
    }
    /* D must carry y as a single linear factor to power k in {1,2}; a y-free
     * denominator is Riccati/Chini/Abel-1st territory (already claimed upstream).
     * The numerator degree is unconstrained — the linear term of the ODE (y'-a y = ...)
     * inflates deg_N to 3 once combined over the denominator (§2.2.17-1606), yet the
     * scaling u=y/s still separates it; the exact A*B split below is the real gate. */
    Expr* degDe = c2("Exponent", expr_copy(D), expr_new_symbol(Yn));
    long degD = as_int(degDe);
    expr_free(degDe);
    if (degD < 1 || degD > 2) { expr_free(omega); expr_free(N); expr_free(D); return NULL; }

    /* D == c(x) (P1 Y + P0)^k with a single linear factor: the root has multiplicity k,
     * so s = -Coeff[D,Y,k-1] / (k Coeff[D,Y,k]); verify D == Coeff[D,Y,k] (Y - s)^k. */
    Expr* ck  = c3("Coefficient", expr_copy(D), expr_new_symbol(Yn), expr_new_integer(degD));
    Expr* ck1 = c3("Coefficient", expr_copy(D), expr_new_symbol(Yn), expr_new_integer(degD - 1));
    if (ds_is_zero(ck)) { expr_free(omega); expr_free(N); expr_free(D); expr_free(ck); expr_free(ck1); return NULL; }
    Expr* s = ds_simplify(c2(SYM_Times,
                  c2(SYM_Times, expr_new_integer(-1), ck1),                 /* -Coeff[k-1] */
                  powneg1(c2(SYM_Times, expr_new_integer(degD), expr_copy(ck)))));  /* /(k Coeff[k]) */
    /* structural check: D - Coeff[k] (Y - s)^k == 0 */
    Expr* model = c2(SYM_Times, expr_copy(ck),
                     c2(SYM_Power, c2(SYM_Subtract, expr_new_symbol(Yn), expr_copy(s)),
                        expr_new_integer(degD)));
    Expr* dchk = c2(SYM_Subtract, expr_copy(D), model);
    bool single_linear = ds_is_zero(dchk);
    expr_free(dchk); expr_free(ck); expr_free(D); expr_free(N);
    if (!single_linear || !ds_free_of(s, Yn) || ds_is_zero(s)) {
        expr_free(omega); expr_free(s); return NULL;
    }

    /* u = y/s : substitute Y = s U, form u' = omega_u = (omega(x, sU) - s' U)/s. */
    Expr* sU = c2(SYM_Times, expr_copy(s), expr_new_symbol(Un));
    Expr* omega_sU = ds_subst(omega, expr_new_symbol(Yn), sU);         /* consumes omega */
    Expr* sp = ds_d(expr_copy(s), expr_new_symbol(xvar));
    Expr* omega_u = ds_simplify(c2(SYM_Times,
                        c2(SYM_Subtract, omega_sU, c2(SYM_Times, sp, expr_new_symbol(Un))),
                        powneg1(expr_copy(s))));                       /* /s */

    /* Separate omega_u == A(x) B(U) by sampling (A = omega_u|_{U=U0},
     * B = omega_u|_{x=x0} / omega_u(x0,U0)); accept only on an exact split. */
    static const int samples[] = { 2, 3, 5, -2, -3 };
    size_t nsamp = sizeof(samples) / sizeof(samples[0]);
    Expr* A = NULL; Expr* B = NULL;
    for (size_t ai = 0; ai < nsamp && !A; ai++) {
        for (size_t bi = 0; bi < nsamp && !A; bi++) {
            Expr* denom = ds_subst(ds_subst(expr_copy(omega_u),
                              expr_new_symbol(xvar), expr_new_integer(samples[ai])),
                              expr_new_symbol(Un), expr_new_integer(samples[bi]));
            denom = eval_and_free(denom);
            if (ds_is_zero(denom)) { expr_free(denom); continue; }
            Expr* Acand = ds_subst(expr_copy(omega_u), expr_new_symbol(Un), expr_new_integer(samples[bi]));
            Expr* Bx0   = ds_subst(expr_copy(omega_u), expr_new_symbol(xvar), expr_new_integer(samples[ai]));
            Expr* Bcand = c2(SYM_Times, Bx0, powneg1(expr_copy(denom)));
            expr_free(denom);
            Expr* prod  = c2(SYM_Times, expr_copy(Acand), expr_copy(Bcand));
            Expr* chk   = c2(SYM_Subtract, expr_copy(omega_u), prod);
            if (air_split_maybe_zero(chk, xvar, Un) && ds_is_zero(chk)) { A = Acand; B = Bcand; }
            else { expr_free(Acand); expr_free(Bcand); }
            expr_free(chk);
        }
    }
    expr_free(omega_u);
    if (!A) { expr_free(s); return NULL; }

    /* Integrate[1/B, U] and Integrate[A, x]; both must be elementary. */
    g_integrate_quiet++;
    Expr* lhsInt = ds_integrate(powneg1(B), expr_new_symbol(Un));      /* consumes B */
    Expr* rhsInt = ds_integrate(A, expr_new_symbol(xvar));             /* consumes A */
    g_integrate_quiet--;
    if (ds_has_head(lhsInt, SYM_Integrate) || ds_has_head(rhsInt, SYM_Integrate)) {
        expr_free(s); expr_free(lhsInt); expr_free(rhsInt); return NULL;
    }

    /* G(x,y) = (Integrate[1/B,U] with U -> y[x]/s) - Integrate[A,x];  relation G == C[1]. */
    Expr* ybys = c2(SYM_Times, ds_make_funcapp(yname, 0, xvar), powneg1(s));  /* consumes s */
    Expr* G = ds_subst(lhsInt, expr_new_symbol(Un), ybys);
    G = c2(SYM_Subtract, G, rhsInt);

    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_abel_air(Expr* res) {
    return dsolve_method_builtin_implicit(res, dsolve_abel_air_try);
}

void dsolve_abel_air_init(void) {
    symtab_add_builtin("DSolve`AbelAIR", builtin_dsolve_abel_air);
    symtab_get_def("DSolve`AbelAIR")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`AbelAIR",
        "DSolve`AbelAIR[eqn, y, x] solves an Abel-2nd-kind / rational-in-y first-order "
        "ODE y' == N(x,y)/D(x,y), with D a single linear y-factor c(x)(P1 y+P0)^k, by "
        "the scaling substitution u = y/s (s = -P0/P1) that makes it separable; returns "
        "the implicit first integral. A partial implementation of the Abel Invariant "
        "Rational method (the full rational-invariant solvable-class table is deferred).");
}
