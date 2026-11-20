/* integrate_risch_rde.c — Bronstein Chapter 6 Risch differential equation.
 *
 * Extracted from integrate_risch_transcendental.c to modularise the RDE stack.
 * See integrate_risch_rde.h for the layering (base field C(x) + recursive tower
 * lift) and the shared-helper contract.  The small build-and-evaluate / poly
 * primitives (rt_eval*, rt_degree, rt_is_zero, rt_tower_deriv, ...) remain
 * defined in integrate_risch_transcendental.c and are declared in the header.
 */

#include "integrate_risch_rde.h"

#include "expr.h"
#include "eval.h"
#include "print.h"
#include "symtab.h"
#include "sym_intern.h"
#include "flint_bridge.h"
#include "risch_structure.h"
#include "risch_field_integrate.h"   /* g_rt_decision: the auxiliary-integral scope */

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ==================================================================== */
/* Bronstein Chapter 6 — the Risch differential equation  D y + f y = g. */
/*                                                                       */
/* Base field C(x) with D = d/dx (Stage A).  Replaces the undetermined-  */
/* coefficient ansatz with the rational one-step-reduction algorithm     */
/* (WeakNormalizer -> RdeNormalDenominator -> RdeBoundDegreeBase ->       */
/* SPDE -> PolyRischDENoCancel1/2), all in polynomial-gcd time.  This is  */
/* the fix for the Davenport-progressive-reduction blow-up of the former */
/* SolveAlways ansatz (Bronstein, J. Symbolic Comput. 9 (1990) 49-60,    */
/* and Symbolic Integration I, 2nd ed., Chapter 6).  In the base field    */
/* every irreducible is normal for d/dx, so there is no special part and  */
/* SplitFactor is the identity — the algorithms simplify accordingly.     */
/* Correct by construction: SPDE's polynomial identity certifies          */
/* (q E^(i u))' = p E^(i u); returns NULL for a genuinely non-elementary  */
/* term.                                                                  */
/* ==================================================================== */

/* ---- Algebraic-constant mode (C completeness campaign, step 7) ----------
 * Over Q(alpha) (numeric radicals or Root objects in the coefficients), Together
 * and PolynomialRemainder leave a zero coefficient unreduced (Sqrt[2] Sqrt[3] -
 * Sqrt[6], (3 + 3 Sqrt[2])/(1 + Sqrt[2]) - 3), so degrees, divisibility and the
 * non-cancellation loop's leading terms were wrong, and the base solver could
 * only mark its "no" uncertain.  While g_rde_alg is set (rde_base switches it on
 * when every coefficient is a number), every polynomial in the base variable is
 * kept with RootReduced coefficients -- canonical in Q(alpha) -- and division,
 * gcd and the extended gcd are exact Euclidean loops over them. */
static int   g_rde_alg   = 0;      /* tc_register_unwind_counter (integrate_risch_transcendental_init) */
static Expr* g_rde_alg_x = NULL;   /* borrowed: the base variable while the mode is on */
int* rde_alg_mode_ptr(void) { return &g_rde_alg; }
static Expr* rde_algc(Expr* e);                          /* consumes e */
static bool  rde_alg_frac(Expr* e, Expr** n, Expr** d);   /* e borrowed */
static Expr* rde_alg_divmod(Expr* a, Expr* b, Expr** rem);
static Expr* rde_alg_gcd(Expr* a, Expr* b);
static Expr* rde_alg_xgcd_s(Expr* b, Expr* a);
static bool  rde_has_alg_const(const Expr* e);

/* ---- Parametric mode (C completeness campaign, step 3a) -----------------
 * The parametric Risch DE  D y + f y = g0 + Sum c_i g_i  asks for y AND
 * constants c_i ([SI] ch. 7, ParamRischDE).  rde_tower_param runs the ordinary
 * solver with the c_i as named symbols (g_rde_par_unk) in the right-hand side.
 * The operator side f never contains them, so every decision on f is unchanged;
 * generic degrees and valuations of a c-linear right-hand side are already the
 * max / min over its components, and its Together denominator is their lcm.
 * What changes is every EXIT on the right-hand side ("c != 0 with n < 0", a
 * remainder that must vanish, a coefficient above the degree bound, a
 * non-polynomial reduced equation): while the tested quantity contains an
 * unknown, the offending part is dropped and the solve continues.  The
 * identity gate then records the residual D y + f y - g, linear in the c_i, as
 * a constraint instead of testing it, and rde_tower_param solves the
 * constraints.
 *   Completeness: if (y*, c*) solves the equation, every dropped part vanishes
 *   at c*, so the run computes what the ordinary solver computes on g(c*),
 *   and y(c*) is a solution: the constraints are consistent.
 *   Soundness: a solution of the constraints makes y(c) a solution (it is the
 *   residual), and an inconsistent system proves no c works.
 * The one nonlinear site, the integration in K behind the log-derivative
 * branch (rde_cancel_logderiv, and D H = c), is bypassed: the cancellation case
 * runs the linear top-down pass of rde_polyrischde_cancel_par, which adds each
 * homogeneous term lambda z with lambda a FRESH unknown (rde_par_fresh), solved
 * with the others.  A site that cannot be converted (the tangent top, an
 * undecided log-derivative test) marks the scope instead. */
static int     g_rde_par     = 0;      /* tc_register_unwind_counter (integrate_risch_transcendental_init) */
static Expr*   g_rde_par_unk = NULL;   /* owned List: the caller's unknowns, then the fresh lambdas */
static Expr**  g_rde_par_cons = NULL;  /* owned constraints (each must vanish) */
static size_t  g_rde_par_ncons = 0, g_rde_par_cap = 0, g_rde_par_nfresh = 0;
int* rde_par_mode_ptr(void) { return &g_rde_par; }

/* Does e contain one of the unknowns of the running parametric solve? */
static bool rde_par_in(const Expr* e) {
    if (!g_rde_par || !g_rde_par_unk || !e) return false;
    for (size_t i = 0; i < g_rde_par_unk->data.function.arg_count; i++)
        if (!rt_free_of_x((Expr*)e, g_rde_par_unk->data.function.args[i])) return true;
    return false;
}
static void rde_par_require(Expr* e) {              /* consumes e */
    if (g_rde_par_ncons == g_rde_par_cap) {
        g_rde_par_cap = g_rde_par_cap ? 2 * g_rde_par_cap : 8;
        g_rde_par_cons = realloc(g_rde_par_cons, g_rde_par_cap * sizeof(Expr*));
    }
    g_rde_par_cons[g_rde_par_ncons++] = e;
}
/* A fresh unknown constant (a homogeneous coefficient), appended to the
 * unknowns of the running solve.  Returns a new reference to the symbol. */
static Expr* rde_par_fresh(void) {
    char nm[32]; snprintf(nm, sizeof(nm), "rmPl%zu", g_rde_par_nfresh++);
    size_t n = g_rde_par_unk->data.function.arg_count;
    Expr** a = malloc((n + 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(g_rde_par_unk->data.function.args[i]);
    a[n] = expr_new_symbol(nm);
    Expr* nl = expr_new_function(expr_new_symbol("List"), a, n + 1);
    free(a);
    expr_free(g_rde_par_unk);
    g_rde_par_unk = nl;
    return expr_new_symbol(nm);
}

/* Numerator / Denominator of Together[e] (e not consumed). */
static Expr* rde_num(Expr* e) {
    if (g_rde_alg) {
        Expr *n, *d;
        if (rde_alg_frac(e, &n, &d)) { expr_free(d); return n; }
    }
    Expr* r = rt_eval1("Numerator", rt_eval1("Together", expr_copy(e)));
    return g_rde_alg ? rde_algc(r) : r;
}
static Expr* rde_den(Expr* e) {
    if (g_rde_alg) {
        Expr *n, *d;
        if (rde_alg_frac(e, &n, &d)) { expr_free(n); return d; }
    }
    Expr* r = rt_eval1("Denominator", rt_eval1("Together", expr_copy(e)));
    return g_rde_alg ? rde_algc(r) : r;
}
/* Univariate polynomial primitives in x (args not consumed). */
static Expr* rde_gcd(Expr* a, Expr* b, Expr* x) {
    (void)x;   /* PolynomialGCD infers the variable; a third arg is a 3rd poly */
    /* Mathematically gcd(p, 0) = p (up to a unit), but Mathilda's PolynomialGCD
     * returns 1 for a zero operand.  The SPDE (Thm 6.4.1) relies on the correct
     * value in its b = 0 sub-problem: with the wrong gcd = 1 the divisibility
     * test g | c passes spuriously and the subsequent exact PolynomialQuotient
     * silently truncates the nonzero remainder, fabricating a bogus q = 0
     * "solution" to an unsolvable RDE (RISCH_AUDIT_A4.md Finding 1 — e.g.
     * Risch`RischDE[-2/x, x, x] wrongly returned 0).  Guard the zero-operand
     * case explicitly so the base RDE stack stays authoritative. */
    if (rt_is_zero(a)) return expr_copy(b);
    if (rt_is_zero(b)) return expr_copy(a);
    if (g_rde_alg) return rde_alg_gcd(a, b);
    return rt_eval2("PolynomialGCD", expr_copy(a), expr_copy(b));
}
static Expr* rde_quot(Expr* a, Expr* b, Expr* x) {
    if (g_rde_alg && expr_eq(x, g_rde_alg_x)) {
        Expr* r = NULL;
        Expr* q = rde_alg_divmod(a, b, &r);
        if (r) expr_free(r);
        if (q) return q;
    }
    return rt_eval3("PolynomialQuotient", expr_copy(a), expr_copy(b), expr_copy(x));
}
static Expr* rde_rem(Expr* a, Expr* b, Expr* x) {
    if (g_rde_alg && expr_eq(x, g_rde_alg_x)) {
        Expr* r = NULL;
        Expr* q = rde_alg_divmod(a, b, &r);
        if (q) expr_free(q);
        if (r) return r;
    }
    return rt_eval3("PolynomialRemainder", expr_copy(a), expr_copy(b), expr_copy(x));
}
/* Expanded / cancelled arithmetic (args not consumed). */
static Expr* rde_mul(Expr* a, Expr* b) {
    Expr* r = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(a), expr_copy(b) }, 2));
    return g_rde_alg ? rde_algc(r) : r;
}
static Expr* rde_add(Expr* a, Expr* b) {
    Expr* r = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(a), expr_copy(b) }, 2));
    return g_rde_alg ? rde_algc(r) : r;
}
static Expr* rde_sub(Expr* a, Expr* b) {
    Expr* r = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(b) }, 2) }, 2));
    return g_rde_alg ? rde_algc(r) : r;
}
/* Field quotient a/b, FULLY reduced (args not consumed).  Uses Together, not
 * Cancel: `Cancel` does not reduce a ratio whose operands are themselves nested
 * rationals — e.g. Cancel[(1+x)/(1+1/x)] stays put, while Together gives x.  The
 * coefficient field k(x) of a tower monomial routinely produces such forms
 * (lc(b) = (1+x)/x printed as 1+1/x), and an unreduced quotient makes leading
 * terms fail to cancel in the SPDE / PolyRischDE loops, losing genuine solutions
 * (RISCH_AUDIT_A4.md — the nested-exponential Laurent case). */
static Expr* rde_divq(Expr* a, Expr* b) {
    Expr* r = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(b), expr_new_integer(-1) }, 2) }, 2));
    return g_rde_alg ? rde_algc(r) : r;
}
/* Canonical form of a tower field element (arg not consumed):
 * Expand[Cancel[Together[e]]] — Together first combines nested rationals into a
 * single reduced fraction (so lc(b) = (1+x)/x is not left as 1+1/x), Cancel strips
 * common factors, and Expand distributes a monomial across a sum (so a genuine
 * polynomial-in-tau is not hidden as t^2 (…-1/(t^2 …))).  This is the ONE
 * normalizer the Bronstein field-arithmetic boxes (RdeNormalDenominator,
 * RdeSpecialDenomExp, SPDE, PolyRischDE) must use so that PolynomialQ / degree /
 * leading-coefficient tests see the true polynomial and leading terms cancel;
 * plain Cancel or Expand alone each miss one of the two reductions
 * (RISCH_AUDIT_A4.md — the rational-coefficient and nested-exponential Laurent
 * regressions traced to exactly these two gaps). */
static Expr* rde_canon(Expr* e) {
    Expr* r = rt_eval1("Expand", rt_eval1("Cancel", rt_eval1("Together", expr_copy(e))));
    return g_rde_alg ? rde_algc(r) : r;
}
/* Leading coefficient of poly p in x (0 if p == 0). */
static Expr* rde_lc(Expr* p, Expr* x) {
    long d = rt_degree(p, x);
    return d < 0 ? expr_new_integer(0) : rt_coeff(p, x, d);
}
/* Parametric mode: cc minus its top term in v (the term a bound says must
 * vanish; its vanishing is implied by the final residual constraint).  cc is
 * consumed. */
static Expr* rde_par_strip_top(Expr* cc, Expr* v) {
    long d = rt_degree(cc, v);
    if (d < 0) return cc;
    Expr* top = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1),
        rt_coeff(cc, v, d), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(v), expr_new_integer(d) }, 2) }, 3);
    Expr* sum = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ cc, top }, 2);
    Expr* r = rde_canon(sum);
    expr_free(sum);
    return r;
}
/* Parametric mode: the polynomial part in v of e (a reduced equation's
 * right-hand side whose non-polynomial part must vanish).  e borrowed. */
static Expr* rde_par_polypart(Expr* e, Expr* v) {
    Expr* tog = rt_eval1("Together", expr_copy(e));
    Expr* num = rt_eval1("Numerator", expr_copy(tog));
    Expr* den = rt_eval1("Denominator", tog);
    Expr* q = rt_eval3("PolynomialQuotient", num, den, expr_copy(v));
    Expr* r = q ? rde_canon(q) : NULL;
    if (q) expr_free(q);
    return r;
}
/* d/dx (e not consumed). */
static Expr* rde_dx(Expr* e, Expr* x) {
    return rt_eval2("D", expr_copy(e), expr_copy(x));
}

/* SPDE result tuple: any solution q of a Dq + b q = c with deg(q) <= n is
 * q = alpha*H + beta, where H in k[t], deg(H) <= m and D H + b' H = c'. */
void rde_spde_free(RdeSpde* s) {
    if (s->b) expr_free(s->b);
    if (s->c) expr_free(s->c);
    if (s->alpha) expr_free(s->alpha);
    if (s->beta) expr_free(s->beta);
    s->b = s->c = s->alpha = s->beta = NULL;
}

/* SPDE(a, b, c, D=d/dx, n): Rothstein's degree-reducing recursion (Bronstein
 * Thm 6.4.1 / algorithm SPDE, p.203).  Returns 1 and fills *out on success, or
 * 0 for "no solution" (the equation has no solution of degree <= n in k[x]).
 * Each level reduces the degree bound by deg(a) via one ExtendedEuclidean, so
 * it terminates in O(n / deg(a)) steps — no dense linear solve. */
/* ---- Certification scope (see integrate_risch_rde.h) -------------------- */
static RdeCertState g_rde_cert = { 0, 0, NULL };

RdeCertState rde_cert_begin(void) {
    RdeCertState saved = g_rde_cert;
    g_rde_cert.active = 1; g_rde_cert.uncertain = 0; g_rde_cert.why = NULL;
    return saved;
}
void rde_cert_end(RdeCertState saved) { g_rde_cert = saved; }
bool rde_cert_active(void) { return g_rde_cert.active != 0; }
bool rde_cert_uncertain(void) { return g_rde_cert.uncertain != 0; }
const char* rde_cert_reason(void) { return g_rde_cert.why; }
/* MATHILDA_RDE_CERT_DEBUG=1 traces every uncertain path to stderr: the
 * developer's view of which step of the campaign's ratchet blocked a verdict. */
static int rde_cert_debug(void) {
    static int cached = -1;
    if (cached < 0) { const char* e = getenv("MATHILDA_RDE_CERT_DEBUG"); cached = (e && *e && *e != '0'); }
    return cached;
}
void rde_cert_mark(const char* why) {
    if (rde_cert_debug()) fprintf(stderr, "[rde-cert] %s%s\n", why, g_rde_cert.active ? "" : " (no scope)");
    if (!g_rde_cert.active) return;
    if (!g_rde_cert.uncertain) g_rde_cert.why = why;
    g_rde_cert.uncertain = 1;
}

static int rde_spde(Expr* a, Expr* b, Expr* c, Expr* x, long n, RdeSpde* out) {
    out->b = out->c = out->alpha = out->beta = NULL; out->m = 0;
    if (n < 0) {
        /* parametric: q = 0 is the only candidate; its residual is the constraint */
        if (rt_is_zero(c) || rde_par_in(c)) {
            out->b = expr_new_integer(0); out->c = expr_new_integer(0);
            out->m = 0; out->alpha = expr_new_integer(0); out->beta = expr_new_integer(0);
            return 1;
        }
        return 0;                                   /* c != 0, n < 0: no solution */
    }
    Expr* g = rde_gcd(a, b, x);
    Expr* rem = rde_rem(c, g, x);
    bool gdivc = rt_is_zero(rem) || rde_par_in(rem);   /* parametric: dropped (quotient below) */
    expr_free(rem);
    if (!gdivc) { expr_free(g); return 0; }         /* g does not divide c */
    Expr* a1 = rde_quot(a, g, x);
    Expr* b1 = rde_quot(b, g, x);
    Expr* c1 = rde_quot(c, g, x);
    expr_free(g);
    long da = rt_degree(a1, x);
    if (da == 0) {                                  /* a in k*: base of recursion */
        out->b = rde_divq(b1, a1);
        out->c = rde_divq(c1, a1);
        out->m = n;              /* Bronstein SPDE: return (b/a, c/a, n, 1, 0) */
        out->alpha = expr_new_integer(1);
        out->beta = expr_new_integer(0);
        expr_free(a1); expr_free(b1); expr_free(c1);
        return 1;
    }
    /* (r, z) with b1 r + a1 z = c1, deg(r) < deg(a1); gcd(a1,b1)=1 so the
     * ExtendedEuclidean cofactor s of b1 gives r = (s c1) mod a1. */
    Expr* eg = NULL;
    Expr* s = NULL;
    Expr* s_alg = NULL;                              /* owned, algebraic mode */
    if (g_rde_alg && expr_eq(x, g_rde_alg_x)) s = s_alg = rde_alg_xgcd_s(b1, a1);
    else {
        eg = rt_eval3("PolynomialExtendedGCD", expr_copy(b1), expr_copy(a1), expr_copy(x));
        if (eg && eg->type == EXPR_FUNCTION && eg->data.function.arg_count == 2) {
            Expr* cof = eg->data.function.args[1];
            if (cof && cof->type == EXPR_FUNCTION && cof->data.function.arg_count == 2)
                s = cof->data.function.args[0];
        }
    }
    if (!s) {
        rde_cert_mark("SPDE: extended gcd unavailable");
        if (eg) expr_free(eg);
        expr_free(a1); expr_free(b1); expr_free(c1);
        return 0;
    }
    Expr* sc = rde_mul(s, c1);
    Expr* r = rde_rem(sc, a1, x);
    expr_free(sc);
    if (eg) expr_free(eg);
    if (s_alg) expr_free(s_alg);
    /* z = (c1 - b1 r) / a1 (exact). */
    Expr* b1r = rde_mul(b1, r);
    Expr* c1_b1r = rde_sub(c1, b1r);
    expr_free(b1r);
    Expr* z = rde_quot(c1_b1r, a1, x);
    expr_free(c1_b1r);
    /* recurse SPDE(a1, b1 + Da1, z - Dr, n - da). */
    Expr* da1 = rde_dx(a1, x);
    Expr* b2 = rde_add(b1, da1);
    expr_free(da1);
    Expr* dr = rde_dx(r, x);
    Expr* c2 = rde_sub(z, dr);
    expr_free(dr); expr_free(z);
    RdeSpde sub;
    int rc = rde_spde(a1, b2, c2, x, n - da, &sub);
    expr_free(b1); expr_free(c1); expr_free(b2); expr_free(c2);
    if (!rc) { expr_free(a1); expr_free(r); return 0; }
    /* q = a1*H + (a1*sub.alpha) H_sub + ...  ->  alpha = a1*sub.alpha,
     *    beta = a1*sub.beta + r. */
    out->b = sub.b; out->c = sub.c; out->m = sub.m;   /* adopt inner eqn */
    out->alpha = rde_mul(a1, sub.alpha);
    Expr* a1beta = rde_mul(a1, sub.beta);
    out->beta = rde_add(a1beta, r);
    expr_free(a1beta);
    expr_free(sub.alpha); expr_free(sub.beta);
    expr_free(a1); expr_free(r);
    return 1;
}

/* PolyRischDENoCancel1(b, c, D=d/dx, n) with b != 0 (Bronstein p.208).  The
 * leading terms of Dq and bq never cancel, so the solution is built top-down,
 * one monomial per pass, the degree bound strictly decreasing.  Returns q
 * (owned) or NULL for "no solution". */
/* True when e carries an irrational ALGEBRAIC constant: a numeric radical
 * (Power[number, Rational]) or a Root object, free of every symbol. */
static bool rde_has_alg_const(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == intern_symbol("Root")) return true;
        if (h == intern_symbol("Power") && e->data.function.arg_count == 2) {
            const Expr* b = e->data.function.args[0];
            const Expr* x = e->data.function.args[1];
            bool numb = b->type == EXPR_INTEGER || b->type == EXPR_BIGINT
                || (b->type == EXPR_FUNCTION && b->data.function.head->type == EXPR_SYMBOL
                    && b->data.function.head->data.symbol.name == intern_symbol("Rational"));
            bool ratx = x->type == EXPR_FUNCTION && x->data.function.head->type == EXPR_SYMBOL
                && x->data.function.head->data.symbol.name == intern_symbol("Rational");
            if (numb && ratx) return true;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rde_has_alg_const(e->data.function.args[i])) return true;
    return false;
}

/* A numeric algebraic quotient brought to RootReduce's canonical form (owned
 * in, owned out); anything else is returned unchanged.  Together does not
 * rationalise a numeric denominator: (2 + 2 Sqrt[2])/(2 + Sqrt[2]) stays put,
 * the no-cancellation loop's leading terms then fail to cancel, and it
 * reported "no solution" for q' + (2 + Sqrt[2]) q = p with q = Sqrt[2] u^2/(1+u)
 * (found by the F9 false-certificate hunt, 2026-11-17). */
static Expr* rde_alg_reduce(Expr* e) {
    if (!e || !rde_has_alg_const(e)) return e;
    Expr* nq = rt_eval1("NumericQ", expr_copy(e));
    bool num = nq && rt_is_true(nq);
    if (nq) expr_free(nq);
    if (!num) return e;
    Expr* r = rt_eval1("RootReduce", expr_copy(e));
    if (!r) return e;
    expr_free(e);
    return r;
}

/* ---- the algebraic mode's exact arithmetic (see g_rde_alg above) ---- */
/* p, a polynomial in the base variable with numeric coefficients: each
 * coefficient RootReduced.  Consumes p. */
static Expr* rde_alg_poly(Expr* p) {
    Expr* X = g_rde_alg_x;
    Expr* cl = rt_eval2("CoefficientList", p, expr_copy(X));   /* adopts p */
    if (!cl || !rt_head_is(cl, "List")) return cl;
    size_t n = cl->data.function.arg_count;
    Expr** terms = malloc((n ? n : 1) * sizeof(Expr*));
    size_t nt = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* c = rde_alg_reduce(rt_eval1("Together", expr_copy(cl->data.function.args[i])));
        if (!c) continue;
        if (rt_is_zero(c)) { expr_free(c); continue; }
        terms[nt++] = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ c,
            expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_copy(X), expr_new_integer((long)i) }, 2) }, 2);
    }
    expr_free(cl);
    Expr* sum = nt == 0 ? expr_new_integer(0) : expr_new_function(expr_new_symbol("Plus"), terms, nt);
    free(terms);
    return rt_eval_own(sum);
}
/* e as a reduced fraction n/d over Q(alpha)[x], both canonical (owned), for e
 * a rational function of x with numeric algebraic coefficients.  Together bails
 * on an algebraic constant in a denominator (it returns the sum uncombined --
 * the open core bug of exact_poly_div, see the Together memory note), so the
 * fractions are combined here: each term's Numerator/Denominator, summed over
 * the product denominator, the common factor cancelled by the exact gcd. */
static bool rde_alg_frac(Expr* e, Expr** n, Expr** d) {
    *n = *d = NULL;
    Expr* X = g_rde_alg_x;
    if (rt_head_is(e, "Plus")) {
        Expr* N = expr_new_integer(0);
        Expr* D = expr_new_integer(1);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Expr *ni, *di;
            if (!rde_alg_frac(e->data.function.args[i], &ni, &di)) { expr_free(N); expr_free(D); return false; }
            Expr* N2 = rde_alg_poly(rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                expr_new_function(expr_new_symbol("Times"), (Expr*[]){ N, expr_copy(di) }, 2),
                expr_new_function(expr_new_symbol("Times"), (Expr*[]){ ni, expr_copy(D) }, 2) }, 2)));
            Expr* D2 = rde_alg_poly(rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"), (Expr*[]){ D, di }, 2)));
            N = N2; D = D2;
            if (!N || !D) { if (N) expr_free(N); if (D) expr_free(D); return false; }
        }
        *n = N; *d = D;
    } else {
        Expr* nn = rt_eval1("Numerator", expr_copy(e));
        Expr* dd = rt_eval1("Denominator", expr_copy(e));
        if (!nn || !dd || !rt_is_poly(nn, X) || !rt_is_poly(dd, X)) {
            if (nn) expr_free(nn);
            if (dd) expr_free(dd);
            return false;
        }
        *n = rde_alg_poly(rt_eval1("Expand", nn));
        *d = rde_alg_poly(rt_eval1("Expand", dd));
        if (!*n || !*d) { if (*n) expr_free(*n); if (*d) expr_free(*d); *n = *d = NULL; return false; }
    }
    /* cancel the common factor; a constant denominator moves into n */
    if (!rt_is_zero(*n) && rt_degree(*d, X) >= 1) {
        Expr* G = rde_alg_gcd(*n, *d);
        if (G && rt_degree(G, X) >= 1) {
            Expr* r1 = NULL; Expr* r2 = NULL;
            Expr* nq = rde_alg_divmod(*n, G, &r1);
            Expr* dq = rde_alg_divmod(*d, G, &r2);
            if (nq && dq && r1 && r2 && rt_is_zero(r1) && rt_is_zero(r2)) {
                expr_free(*n); expr_free(*d); *n = nq; *d = dq;
            } else { if (nq) expr_free(nq); if (dq) expr_free(dq); }
            if (r1) expr_free(r1);
            if (r2) expr_free(r2);
        }
        if (G) expr_free(G);
    }
    long dd = rt_degree(*d, X);
    Expr* lc = dd >= 0 ? rt_coeff(*d, X, dd) : NULL;
    if (lc && !rt_is_zero(lc)) {                      /* monic denominator */
        Expr* inv = rde_alg_reduce(rt_eval1("Together", expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ lc, expr_new_integer(-1) }, 2)));
        *n = rde_alg_poly(rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"), (Expr*[]){ *n, expr_copy(inv) }, 2)));
        *d = rde_alg_poly(rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"), (Expr*[]){ *d, inv }, 2)));
    } else if (lc) expr_free(lc);
    return *n && *d;
}
/* Canonical form over Q(alpha)(x): a polynomial coefficient-wise; a rational
 * function as its reduced canonical fraction.  Consumes e. */
static Expr* rde_algc(Expr* e) {
    if (!e || !g_rde_alg || !rde_has_alg_const(e)) return e;
    if (rt_is_poly(e, g_rde_alg_x)) return rde_alg_poly(e);
    Expr *n, *d;
    if (!rde_alg_frac(e, &n, &d)) return e;           /* not rational in x: leave it */
    expr_free(e);
    if (rt_degree(d, g_rde_alg_x) == 0 && expr_eq(d, expr_new_integer(1))) { expr_free(d); return n; }
    return rt_eval_own(expr_new_function(expr_new_symbol("Times"), (Expr*[]){ n,
        expr_new_function(expr_new_symbol("Power"), (Expr*[]){ d, expr_new_integer(-1) }, 2) }, 2));
}
/* Exact division of polynomials over Q(alpha): quotient (owned), *rem owned. */
static Expr* rde_alg_divmod(Expr* a, Expr* b, Expr** rem) {
    *rem = NULL;
    Expr* X = g_rde_alg_x;
    Expr* A = rde_alg_poly(rt_eval1("Expand", expr_copy(a)));
    Expr* B = rde_alg_poly(rt_eval1("Expand", expr_copy(b)));
    long db = B ? rt_degree(B, X) : -1;
    if (!A || db < 0 || rt_is_zero(B)) { if (A) expr_free(A); if (B) expr_free(B); return NULL; }
    Expr* lcb = rt_coeff(B, X, db);
    Expr* Q = expr_new_integer(0);
    long guard = rt_degree(A, X) + 2;
    while (!rt_is_zero(A) && guard-- > 0) {
        long da = rt_degree(A, X);
        if (da < db) break;
        Expr* lca = rt_coeff(A, X, da);
        Expr* c = rde_alg_reduce(rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ lca, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(lcb), expr_new_integer(-1) }, 2) }, 2)));
        Expr* term = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ c,
            expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_copy(X), expr_new_integer(da - db) }, 2) }, 2);
        Expr* Qn = rde_alg_poly(rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ Q, expr_copy(term) }, 2)));
        Expr* An = rde_alg_poly(rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ A,
            expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), term, expr_copy(B) }, 3) }, 2)));
        Q = Qn; A = An;
        if (!Q || !A) break;
    }
    expr_free(lcb); expr_free(B);
    if (!Q || !A) { if (Q) expr_free(Q); if (A) expr_free(A); return NULL; }
    *rem = A;
    return Q;
}
/* Monic gcd over Q(alpha) by the Euclidean algorithm on canonical operands. */
static Expr* rde_alg_gcd(Expr* a, Expr* b) {
    Expr* X = g_rde_alg_x;
    Expr* u = rde_alg_poly(rt_eval1("Expand", expr_copy(a)));
    Expr* v = rde_alg_poly(rt_eval1("Expand", expr_copy(b)));
    while (u && v && !rt_is_zero(v)) {
        Expr* r = NULL;
        Expr* q = rde_alg_divmod(u, v, &r);
        if (q) expr_free(q);
        if (!r) { expr_free(u); expr_free(v); return NULL; }
        expr_free(u); u = v; v = r;
    }
    if (v) expr_free(v);
    if (!u) return NULL;
    long d = rt_degree(u, X);
    if (d >= 0) {
        Expr* lc = rt_coeff(u, X, d);
        if (!rt_is_zero(lc)) {
            Expr* inv = rde_alg_reduce(rt_eval1("Together", expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ lc, expr_new_integer(-1) }, 2)));
            u = rde_alg_poly(rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"), (Expr*[]){ u, inv }, 2)));
        } else expr_free(lc);
    }
    return u;
}
/* s with s b == 1 (mod a) over Q(alpha), for coprime a, b (the SPDE cofactor);
 * NULL when gcd(a, b) is not a constant. */
static Expr* rde_alg_xgcd_s(Expr* b, Expr* a) {
    Expr* X = g_rde_alg_x;
    Expr* r0 = rde_alg_poly(rt_eval1("Expand", expr_copy(b)));
    Expr* r1 = rde_alg_poly(rt_eval1("Expand", expr_copy(a)));
    Expr* s0 = expr_new_integer(1);
    Expr* s1 = expr_new_integer(0);
    while (r0 && r1 && !rt_is_zero(r1)) {
        Expr* r2 = NULL;
        Expr* q = rde_alg_divmod(r0, r1, &r2);
        if (!q || !r2) { if (q) expr_free(q); if (r2) expr_free(r2); expr_free(r0); expr_free(r1); expr_free(s0); expr_free(s1); return NULL; }
        Expr* s2 = rde_alg_poly(rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ s0,
            expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), q, expr_copy(s1) }, 3) }, 2)));
        expr_free(r0); r0 = r1; r1 = r2; s0 = s1; s1 = s2;
        if (!s1) break;
    }
    Expr* res = NULL;
    if (r0 && rt_degree(r0, X) == 0 && !rt_is_zero(r0) && s0) {
        Expr* inv = rde_alg_reduce(rt_eval1("Together", expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(r0), expr_new_integer(-1) }, 2)));
        res = rde_alg_poly(rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_copy(s0), inv }, 2)));
    }
    if (r0) expr_free(r0);
    if (r1) expr_free(r1);
    if (s0) expr_free(s0);
    if (s1) expr_free(s1);
    return res;
}
/* An upper bound on [Q(alpha):Q] for the algebraic atoms of e: the product of
 * the radical indices q (b^(p/q)) and the Root degrees, over distinct atoms.
 * RootReduce works in the compositum, and its cost grows fast with the degree
 * (Sqrt[2], Sqrt[3], 2^(1/3) is degree 12: thousands of qqbar factorizations
 * per solve); above the cap the mode stays off and the marks stand. */
static void rc_collect_alg_atoms(Expr* e, Expr*** out, size_t* na, size_t* cap);
static long rde_alg_degree_bound(Expr* f, Expr* g) {
    Expr** atoms = NULL; size_t na = 0, cap = 0;
    rc_collect_alg_atoms(f, &atoms, &na, &cap);
    rc_collect_alg_atoms(g, &atoms, &na, &cap);
    /* powers of one base b^(p/q) generate Q(b^(1/lcm q)): group radicals by base */
    long deg = 1;
    bool* used = calloc(na ? na : 1, sizeof(bool));
    for (size_t i = 0; i < na && deg <= 1000; i++) {
        if (used[i]) continue;
        Expr* a = atoms[i];
        long d = 1;
        if (rt_head_is(a, "Root")) {
            Expr* p = rt_eval_own(expr_new_function(expr_copy(a->data.function.args[0]),
                (Expr*[]){ expr_new_symbol("Risch`Private`rdeDegZ") }, 1));
            Expr* z = expr_new_symbol("Risch`Private`rdeDegZ");
            d = p ? rt_degree(p, z) : 1000;
            if (p) expr_free(p);
            expr_free(z);
        } else {
            Expr* base = a->data.function.args[0];
            long l = 1;
            for (size_t j = i; j < na; j++) {
                Expr* b = atoms[j];
                if (used[j] || rt_head_is(b, "Root") || !expr_eq(b->data.function.args[0], base)) continue;
                used[j] = true;
                Expr* q = b->data.function.args[1]->data.function.args[1];
                long qi = q->type == EXPR_INTEGER ? (long)q->data.integer : 1000;
                long g = l, h = qi;
                while (h) { long t = g % h; g = h; h = t; }
                l = l / (g ? g : 1) * qi;
                if (l > 1000) { l = 1000; break; }
            }
            d = l;
        }
        if (d < 1) d = 1;
        deg *= d;
    }
    free(used);
    free(atoms);
    return deg;
}

/* Are the only symbols of e the variable x (numeric algebraic coefficients)? */
static bool rde_only_var(Expr* e, Expr* x) {
    Expr* v = rt_eval1("Variables", rt_eval1("Together", expr_copy(e)));
    bool ok = v && rt_head_is(v, "List");
    for (size_t i = 0; ok && i < v->data.function.arg_count; i++)
        if (!expr_eq(v->data.function.args[i], x)) ok = false;
    if (v) expr_free(v);
    return ok;
}

static Expr* rde_polyrischde_nocancel1(Expr* b, Expr* c, Expr* x, long n) {
    long db = rt_degree(b, x);
    Expr* lcb = rde_lc(b, x);
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        long dc = rt_degree(cc, x);
        long m = dc - db;
        if (n < 0 || m < 0 || m > n) {
            /* parametric: the top term must vanish; drop it (rde_par_strip_top) */
            if (rde_par_in(cc)) { cc = rde_par_strip_top(cc, x); continue; }
            fail = true; break;
        }
        Expr* lcc = rde_lc(cc, x);
        Expr* coeff = rde_alg_reduce(rde_divq(lcc, lcb));   /* lc(c)/lc(b) */
        expr_free(lcc);
        Expr* p = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ coeff, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(x), expr_new_integer(m) }, 2) }, 2));  /* adopts coeff */
        Expr* nq = rde_add(q, p);
        expr_free(q); q = nq;
        /* c <- c - Dp - b p. */
        Expr* dp = rde_dx(p, x);
        Expr* bp = rde_mul(b, p);
        Expr* t1 = rde_sub(cc, dp);
        Expr* t2 = rde_sub(t1, bp);
        expr_free(cc); cc = t2;
        expr_free(dp); expr_free(bp); expr_free(t1); expr_free(p);
        n = m - 1;
    }
    expr_free(cc); expr_free(lcb);
    if (fail) { expr_free(q); return NULL; }
    return q;
}

/* PolyRischDENoCancel2(b, c, D=d/dx, n) with b == 0 (Bronstein p.209, D=d/dt
 * case): the equation is Dq = c, i.e. antidifferentiation of c.  Returns the
 * polynomial antiderivative q with deg(q) <= n (owned), or NULL if c has no
 * polynomial antiderivative of bounded degree.  (lambda(x) = lc(Dx) = 1,
 * delta(x) = 0 for the base derivation.) */
static Expr* rde_polyrischde_integrate(Expr* c, Expr* x, long n) {
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        long dc = rt_degree(cc, x);
        long m = dc + 1;                            /* deg(c) - delta + 1 */
        if (n < 0 || m < 0 || m > n) {
            if (rde_par_in(cc)) { cc = rde_par_strip_top(cc, x); continue; }
            fail = true; break;
        }
        Expr* lcc = rde_lc(cc, x);
        Expr* coeff = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ lcc, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_integer(m), expr_new_integer(-1) }, 2) }, 2));  /* adopts lcc */
        Expr* p = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ coeff, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(x), expr_new_integer(m) }, 2) }, 2));
        Expr* nq = rde_add(q, p);
        expr_free(q); q = nq;
        Expr* dp = rde_dx(p, x);
        Expr* t1 = rde_sub(cc, dp);
        expr_free(cc); cc = t1;
        expr_free(dp); expr_free(p);
        n = m - 1;
    }
    expr_free(cc);
    if (fail) { expr_free(q); return NULL; }
    return q;
}

/* WeakNormalizer (Bronstein §6.1, p.183), one implementation for the base field
 * and the tower: defined after the k(x)(tau) helpers below (rc_weak_normalizer). */
static Expr* rde_weak_normalizer(Expr* f, Expr* x);

/* RischDE in the base field C(x): solve D y + f y = g for y in C(x).  Returns y
 * (owned) or NULL for "no solution" (the term is non-elementary in this field).
 * Pipeline: weak-normalize f, then RdeNormalDenominator reduces to the
 * polynomial equation a Dq + b q = c (base field has no special denominator),
 * RdeBoundDegreeBase bounds deg(q), SPDE + PolyRischDENoCancel solve it, and
 * y = q / (h * w). */
Expr* rde_base(Expr* f, Expr* g, Expr* x) {
    if (rt_is_zero(g)) return expr_new_integer(0);
    /* algebraic constants and nothing but x: the exact algebraic mode */
    if (!g_rde_alg && x->type == EXPR_SYMBOL && (rde_has_alg_const(f) || rde_has_alg_const(g))
        && rde_only_var(f, x) && rde_only_var(g, x) && rde_alg_degree_bound(f, g) <= 8) {
        g_rde_alg = 1; g_rde_alg_x = x;
        Expr* y = rde_base(f, g, x);
        g_rde_alg = 0; g_rde_alg_x = NULL;
        /* the mode keeps denominators monic (u^2/(Sqrt[2] (1/Sqrt[2] + u)));
         * Cancel restores the natural form (u^2/(1 + Sqrt[2] u)) */
        if (y && !rt_free_of_head(y, "Root")) {        /* RootReduce's Root forms back to radicals */
            Expr* rr = rt_eval1("ToRadicals", expr_copy(y));
            if (rr && rt_free_of_head(rr, "Root")) { expr_free(y); y = rr; } else if (rr) expr_free(rr);
        }
        if (y) { Expr* c = rt_eval1("Cancel", expr_copy(y)); if (c) { expr_free(y); y = c; } }
        return y;
    }

#ifdef USE_FLINT
    /* Native fast path: for the exponential tower the coefficient f is a
     * polynomial over Q(x) and g is a rational function over Q(x); the whole
     * base-field RDE (weak normalization is then a no-op) runs in fmpq_poly,
     * converting f and g straight to FLINT with no evaluator Together/Cancel.
     * This is the dominant, high-degree case (In16/In17/`poly·e^x`) — it
     * collapses from seconds of Expr rational arithmetic to milliseconds. The
     * verdict is authoritative: 1 -> y, 0 -> genuinely no solution (decline),
     * -1 -> out of scope (f rational / not univariate over Q) -> fall through
     * to the Expr path below. */
    /* (over Q only: with an algebraic constant the conversion does not see the
     * constant field, and its "no" was a false certificate --
     * Risch`RischDECertified[2/(u - Sqrt[2]), g, u] for y = Sqrt[3] u^2 said
     * "NoSolution" up to v0.6605; campaign step 7) */
    if (x->type == EXPR_SYMBOL && !rde_has_alg_const(f) && !rde_has_alg_const(g)
        && !rde_par_in(g)) {                         /* parametric: the Expr path converts its exits */
        Expr* y = NULL;
        int nr = flint_rde_base_solve_fg(f, g, x->data.symbol.name, &y);
        if (nr >= 0) return y;              /* y is NULL on nr == 0 (decline) */
    }
#endif

    /* 1. Weak normalization: f <- f - Dw/w, g <- w g, y = z/w. */
    Expr* w = rde_weak_normalizer(f, x);
    Expr* fbar;
    if (rt_is_zero(w)) { rde_cert_mark("WeakNormalizer: zero"); expr_free(w); return NULL; }
    if (w->type == EXPR_INTEGER && w->data.integer == 1) {
        fbar = expr_copy(f);
    } else {
        Expr* dw = rde_dx(w, x);
        Expr* dww = rde_divq(dw, w);
        fbar = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ expr_copy(f), expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), dww }, 2) }, 2));   /* adopts dww */
        if (g_rde_alg) fbar = rde_algc(fbar);
        expr_free(dw);
    }
    Expr* gbar = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(w), expr_copy(g) }, 2));
    if (g_rde_alg) gbar = rde_algc(gbar);

    /* 2. RdeNormalDenominator(fbar, gbar): base field, ds = es = 1. */
    Expr* dn = rde_den(fbar);
    Expr* en = rde_den(gbar);
    Expr* p = rde_gcd(dn, en, x);
    Expr* den_e = rde_dx(en, x);
    Expr* den_num = rde_gcd(en, den_e, x);          /* gcd(en, en') */
    Expr* den_p = rde_dx(p, x);
    Expr* den_den = rde_gcd(p, den_p, x);           /* gcd(p, p')  */
    expr_free(den_e); expr_free(den_p);
    Expr* h = rde_quot(den_num, den_den, x);
    expr_free(den_num); expr_free(den_den);
    /* Guard: en | dn h^2  (else no solution). */
    Expr* h2 = rde_mul(h, h);
    Expr* dnh2 = rde_mul(dn, h2);
    Expr* rem = rde_rem(dnh2, en, x);
    bool ok = rt_is_zero(rem);
    /* parametric: en is the lcm of the components' denominators, so h bounds
     * the denominator of y for every value of the c_i (a smaller en at c* gives
     * a divisor of h); a pole of g that h does not clear must vanish, and the
     * reduced right-hand side below keeps only its polynomial part */
    bool par_guard = !ok && rde_par_in(gbar);
    if (par_guard) ok = true;
    expr_free(rem); expr_free(h2);
    Expr* result = NULL;
    if (ok) {
        /* a = dn h,  b = dn h fbar - dn Dh,  c = dn h^2 gbar.  These are all
         * polynomials by construction; compute them as such (Numerator and
         * exact PolynomialQuotient) rather than via Cancel of a rational, which
         * can fail to reduce a factored/expanded mismatch (e.g. leaving
         * (x+1)(x+2)^2/(x+2)^2 uncancelled). */
        Expr* num_f = rde_num(fbar);                /* dn fbar = num_f (dn=den_f) */
        Expr* num_g = rde_num(gbar);
        Expr* aa = rde_mul(dn, h);                  /* dn h */
        Expr* dh = rde_dx(h, x);
        Expr* dndh = rde_mul(dn, dh);               /* dn Dh */
        Expr* hnf = rde_mul(h, num_f);              /* dn h fbar = h num_f */
        Expr* bb = rde_sub(hnf, dndh);
        Expr* dnh2_en = rde_quot(dnh2, en, x);      /* dn h^2 / en (exact) */
        Expr* cc = rde_mul(dnh2_en, num_g);         /* (dn h^2 / en) num_g */
        expr_free(num_f); expr_free(num_g); expr_free(dh); expr_free(dndh);
        expr_free(hnf); expr_free(dnh2_en);
        if (rde_par_in(gbar) && (par_guard || !rt_is_poly(cc, x))) {
            Expr* full = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(dnh2), expr_copy(gbar) }, 2);
            Expr* pp = rde_par_polypart(full, x);
            expr_free(full);
            if (pp) { expr_free(cc); cc = pp; }
        }
        if (!(rt_is_poly(aa, x) && rt_is_poly(bb, x) && rt_is_poly(cc, x)))
            rde_cert_mark("base: reduced equation not polynomial");
        if (rt_is_poly(aa, x) && rt_is_poly(bb, x) && rt_is_poly(cc, x)) {
            /* 3. Degree bound (RdeBoundDegreeBase, p.199). */
            long da = rt_degree(aa, x), db = rt_degree(bb, x), dc = rt_degree(cc, x);
            long mx = (db > da - 1) ? db : (da - 1);
            long n = dc - mx; if (n < 0) n = 0;
            if (db == da - 1 && da >= 1) {
                Expr* lcb = rde_lc(bb, x); Expr* lca = rde_lc(aa, x);
                Expr* mm = rde_divq(lcb, lca);
                Expr* mmn = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), mm }, 2));   /* -lc(b)/lc(a), adopts mm */
                if (mmn && mmn->type == EXPR_INTEGER) {
                    long mv = (long)mmn->data.integer;
                    long cand = dc - db; if (cand < 0) cand = 0;
                    if (mv > n) n = mv;
                    if (cand > n) n = cand;
                }
                if (mmn) expr_free(mmn);
                expr_free(lcb); expr_free(lca);
            }
            /* 4. SPDE reduce, 5. solve bounded (Expr fallback path — the
             * exponential common case is handled natively in fmpq_poly at the
             * top of rde_base; see the flint_rde_base_solve_fg dispatch). */
            RdeSpde sp;
            if (rde_spde(aa, bb, cc, x, n, &sp)) {
                Expr* H;
                if (rt_is_zero(sp.b))
                    H = rde_polyrischde_integrate(sp.c, x, sp.m);
                else
                    H = rde_polyrischde_nocancel1(sp.b, sp.c, x, sp.m);
                /* Over algebraic constants the polynomial solve's zero tests
                 * are only as good as the normalisation (rde_alg_reduce): a
                 * failure there is not a theorem.  Rational coefficients (with
                 * parameters) keep the authoritative "no". */
                if (H) {
                    Expr* aH = rde_mul(sp.alpha, H);
                    Expr* q = rde_add(aH, sp.beta);
                    expr_free(aH); expr_free(H);
                    /* y = q / (h w). */
                    Expr* hw = rde_mul(h, w);
                    result = rde_divq(q, hw);
                    expr_free(hw); expr_free(q);
                }
                rde_spde_free(&sp);
            }
        }
        expr_free(aa); expr_free(bb); expr_free(cc);
    }
    expr_free(dnh2);
    expr_free(dn); expr_free(en); expr_free(p); expr_free(h);
    expr_free(fbar); expr_free(gbar); expr_free(w);

    /* Correct-by-construction gate (Expr fallback path).  Unlike rde_tower
     * (which gates every result with the exact tower identity at rde_tower's
     * tail), the base Expr reconstruction above trusted its PolynomialQuotients
     * to be exact — the escape route for RISCH_AUDIT_A4.md Finding 1.  Verify
     * D[y] + f y - g == 0 so any residual algebra slip declines cleanly instead
     * of shipping a wrong closed form.  (The FLINT fast path returned earlier;
     * it is separately audited authoritative — F5.) */
    if (result) {
        Expr* dy = rde_dx(result, x);
        Expr* fy = rde_mul(f, result);
        Expr* s1 = rde_add(dy, fy);
        Expr* resid = rde_sub(s1, g);
        Expr* tog = rt_eval1("Together", resid);
        Expr* num = tog ? rt_eval1("Numerator", tog) : NULL;
        Expr* exn = num ? rt_eval1("Expand", num) : NULL;
        if (exn && g_rde_alg) exn = rde_algc(exn);    /* canonical zero over Q(alpha) */
        bool good = exn && rt_is_zero(exn);
        /* parametric: the residual, linear in the c_i, is the constraint (a
         * residual the unknowns cancelled out of is one too: dropped parts that
         * no c removes, so the system is then inconsistent) */
        if (!good && exn && (rde_par_in(exn) || rde_par_in(g))) { rde_par_require(expr_copy(exn)); good = true; }
        expr_free(dy); expr_free(fy); expr_free(s1);
        if (exn) expr_free(exn);
        if (!good) { rde_cert_mark("base: identity gate"); expr_free(result); result = NULL; }
    }
    /* Outside the exact algebraic mode (a parameter beside the constants, or a
     * compositum above the degree cap), no zero test of this path is canonical
     * over Q(alpha): a "no" from any of its exits -- the normal-denominator
     * guard, SPDE, the polynomial solve -- is not a theorem (the guard said "no"
     * for y = 1 + Sqrt[3]/(u - 2^(1/3)), f = 1/(u - Sqrt[2])). */
    if (!result && !g_rde_alg && (rde_has_alg_const(f) || rde_has_alg_const(g)))
        rde_cert_mark("base: algebraic constants outside the exact mode");
    return result;
}

/* Solve the Risch differential equation  q' + i u' q = p  for q in the base
 * field C(x), via the Bronstein one-step-reduction algorithm above.  The
 * coefficient is f = i u'; returns q (owned) or NULL when the term is
 * non-elementary in this field (e.g. E^(-x^2), E^(1/x)).  Delegates to rde_base,
 * the COMPLETE base-field solver (weak normalization + normal-denominator
 * reduction + SPDE + polynomial non-cancellation solve) — so every NULL it
 * returns is an authoritative "no rational solution", never a bounded-ansatz
 * decline.  u and p may be rational in x (not just polynomial). */
Expr* rt_solve_rde(Expr* p, long i, Expr* u, Expr* x) {
    if (i == 0) return NULL;
    Expr* up = rt_eval2("D", expr_copy(u), expr_copy(x));   /* u' */
    if (!up) return NULL;
    Expr* f = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(i), up }, 2));          /* f = i u', adopts up */
    Expr* q = rde_base(f, p, x);
    expr_free(f);
    return q;
}

/* Bronstein Chapter 6 over the TOWER — recursive Risch differential      */
/* equation  D_tower[y] + f y = g  for y in K_m = C(x)(t_0..t_m)  (Gap 1). */
/*                                                                        */
/* The base-field stack (rde_base and friends, above) hardcodes D = d/dx  */
/* over C(x).  This block lifts the SAME Bronstein algorithm boxes to an   */
/* arbitrary tower field by parameterising the derivation and the         */
/* polynomial variable through an RdeCtx: the current top monomial        */
/* tau = t_m (kind RT_LOG / RT_EXP) with the lower tower vars             */
/* {x, t_0..t_{m-1}} = k treated as a transcendental coefficient field.    */
/* Increment 1a implements the PRIMITIVE (RT_LOG) top, non-cancellation    */
/* branch; the exponential top (RdeSpecialDenomExp) and the cancellation   */
/* branches land in later increments (until then rde_tower returns NULL    */
/* and rt_field_rde falls back to the bounded ansatz, so nothing           */
/* regresses).  Every q returned is certified by the exact tower-variable  */
/* identity D_tower[q] + f q == g before it ships, so a bug in the new      */
/* reduction can only ever DECLINE, never emit a wrong closed form.        */
/*                                                                          */
/* R1 (field-correct polynomial algebra over k(x)[tau]): the               */
/* variable-explicit PolynomialQuotient / PolynomialRemainder /            */
/* PolynomialExtendedGCD do TRUE field division over k(x) (verified: they   */
/* introduce no lower-variable content); but PolynomialGCD is content-      */
/* INCLUDING (gcd(x t, x) = x, not 1), so rc_gcd runs a monic Euclidean     */
/* remainder loop in tau instead.                                          */
/* ==================================================================== */


static long  rc_deg(Expr* e, RdeCtx* C)          { return rt_degree(e, C->tau); }
static Expr* rc_coeff(Expr* e, long i, RdeCtx* C){ return rt_coeff(e, C->tau, i); }
static Expr* rc_lc(Expr* e, RdeCtx* C) {
    long d = rc_deg(e, C);
    return d < 0 ? expr_new_integer(0) : rc_coeff(e, d, C);
}
/* Derivation D (d/dx at base; the full tower derivation for m >= 0 — correct on
 * any element of K_m since d/dt_j vanishes for j > m). */
static Expr* rc_D(Expr* e, RdeCtx* C) {
    if (C->m < 0) return rde_dx(e, C->x);
    /* canonical: the raw tower derivative of a monic (1 + t^2 + x t^2)/(1 + x)
     * keeps a t^2 coefficient 1/(1 + x) - (1 + x)/(1 + x)^2 that is 0 but not
     * written as 0, so rc_deg overstated the degree and rc_gcd divided by it */
    Expr* d = rt_tower_deriv(e, C->T, C->x);
    if (!d) return d;
    Expr* r = rde_canon(d); expr_free(d); return r;
}
/* Quotient / remainder in tau over k(x), CANONICAL (rde_canon).  At m >= 0 the
 * coefficients are rational functions of the lower field, and PolynomialRemainder
 * leaves them uncombined: a remainder that is 0 in k(x)[tau] then fails
 * rt_is_zero, so the en | dn h^2 guard of RdeNormalDenominator reported a false,
 * CERTIFIED "no solution", and the Euclidean gcd divided by a zero leading
 * coefficient (1/0).  Shape: Risch`RischDECertified[f, D[1/x] + f/x, x] with
 * f = D[G]/(2 G), G = 1 + (1 + x) Log[x]^2.  At the base field (m < 0) the
 * coefficients are rational numbers and the plain operation is exact. */
static Expr* rc_quot(Expr* a, Expr* b, RdeCtx* C) {
    Expr* q = rde_quot(a, b, C->tau);
    if (C->m < 0 || !q) return q;
    Expr* r = rde_canon(q); expr_free(q); return r;
}
static Expr* rc_rem(Expr* a, Expr* b, RdeCtx* C) {
    Expr* q = rde_rem(a, b, C->tau);
    if (C->m < 0 || !q) return q;
    Expr* r = rde_canon(q); expr_free(q); return r;
}
static Expr* rc_xgcd(Expr* a, Expr* b, RdeCtx* C) {
    return rt_eval3("PolynomialExtendedGCD", expr_copy(a), expr_copy(b), expr_copy(C->tau));
}
/* gcd in tau over the coefficient field k(x).  At base (tau = x over Q) the
 * multivariate PolynomialGCD IS the field gcd; for m >= 0 it is content-including,
 * so compute the monic Euclidean gcd in tau explicitly (R1). */
static Expr* rc_gcd(Expr* a, Expr* b, RdeCtx* C) {
    if (C->m < 0) return rde_gcd(a, b, C->x);
    /* canonical operands (rde_canon, not Cancel): a zero leading coefficient left
     * unreduced makes the monic step below divide by 0 (see rc_rem) */
    Expr* u = rde_canon(a);
    Expr* v = rde_canon(b);
    while (!rt_is_zero(v)) {
        Expr* r = rc_rem(u, v, C);
        expr_free(u); u = v; v = r;
    }
    expr_free(v);
    long d = rc_deg(u, C);
    if (d >= 0) {                                   /* make monic in tau */
        Expr* lc = rc_coeff(u, d, C);
        if (!rt_is_zero(lc)) { Expr* mo = rde_divq(u, lc); expr_free(u); u = mo; }
        expr_free(lc);
    }
    return u;
}
static bool rc_ispoly(Expr* e, RdeCtx* C) {
    Expr* q = rt_eval2("PolynomialQ", expr_copy(e), expr_copy(C->tau));
    bool r = q && q->type == EXPR_SYMBOL && q->data.symbol.name == intern_symbol("True");
    if (q) expr_free(q);
    return r;
}
/* Monic-in-tau denominator of e as an element of k(x)(tau): making the
 * denominator monic strips all pure-k(x) content into the numerator, so the
 * normal-denominator reduction sees only genuine tau-poles (a spurious x-content
 * factor would otherwise corrupt the gcd(den, D den) squarefree test). */
static Expr* rc_den(Expr* e, RdeCtx* C) {
    Expr* d = rde_den(e);
    if (C->m < 0) return d;
    long dd = rt_degree(d, C->tau);
    if (dd <= 0) { expr_free(d); return expr_new_integer(1); }  /* tau-free: no pole */
    Expr* lc = rt_coeff(d, C->tau, dd);
    Expr* mon = rde_divq(d, lc);
    expr_free(d); expr_free(lc);
    return mon;
}

/* Lowest power of tau with a nonzero coefficient in the polynomial p (the tau-adic
 * order at tau = 0); 0 when p == 0 or tau does not divide p. */
static long rc_ord_tau(Expr* p, RdeCtx* C) {
    if (rt_is_zero(p)) return 0;
    Expr* cl = rt_eval2("CoefficientList", expr_copy(p), expr_copy(C->tau));
    long o = 0;
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol.name == intern_symbol("List")) {
        for (size_t i = 0; i < cl->data.function.arg_count; i++)
            if (!rt_is_zero(cl->data.function.args[i])) { o = (long)i; break; }
    }
    if (cl) expr_free(cl);
    return o;
}
/* tau-adic valuation nu_tau(e) of a field element e in k(x)(tau): order of tau in
 * the numerator minus the order in the denominator (negative == a tau-pole). */
static long rc_val_tau(Expr* e, RdeCtx* C) {
    Expr* num = rde_num(e);
    Expr* den = rde_den(e);
    long v = rc_ord_tau(num, C) - rc_ord_tau(den, C);
    expr_free(num); expr_free(den);
    return v;
}
/* tau^k (k >= 0) as a polynomial. */
static Expr* rc_tau_pow(long k, RdeCtx* C) {
    return expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_copy(C->tau), expr_new_integer(k) }, 2);
}
/* Normal part of a (monic-in-tau) denominator (SplitFactor, Bronstein §3.5): for
 * an exponential monomial the only special irreducible is tau, so strip tau^v; for
 * a primitive (RT_LOG) / base monomial every squarefree factor is normal. */
static Expr* rc_normal_part(Expr* den, RdeCtx* C) {
    if (C->m >= 0 && C->T->kind[C->m] == RT_TAN && C->T->tsg[C->m] == 1) {
        /* hypertangent: the special irreducible is tau^2 + 1 (strip its power) */
        Expr* p = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ rc_tau_pow(2, C), expr_new_integer(1) }, 2);
        Expr* d = expr_copy(den);
        while (d && rc_deg(d, C) >= 2) {
            Expr* r = rc_rem(d, p, C);
            bool div = r && rt_is_zero(r);
            if (r) expr_free(r);
            if (!div) break;
            Expr* q = rc_quot(d, p, C);
            expr_free(d); d = q;
        }
        expr_free(p);
        return d;
    }
    if (C->m < 0 || C->T->kind[C->m] != RT_EXP) return expr_copy(den);
    long v = rc_ord_tau(den, C);
    if (v <= 0) return expr_copy(den);
    Expr* tv = rc_tau_pow(v, C);
    Expr* q = rc_quot(den, tv, C);
    expr_free(tv);
    return q;
}
/* RdeNormalDenominator over k(x)(tau) (Bronstein §6.1, Cor 6.1.1), lifted to the
 * tau-field ops and to an arbitrary monomial via the SplitFactor normal/special
 * split: reduce D y + f y = g to a Dq + b q = c with q = y h in k<tau>, a in
 * k[tau] special-free.  For a primitive top b, c are already polynomials; for an
 * exponential top they may still carry tau-poles (the special part), cleared next
 * by RdeSpecialDenomExp.  Returns 1 and fills *a,*b,*c,*h (owned); 0 for the
 * authoritative "no solution" (en does not divide dn h^2, Cor 6.1.1(ii)). */
static int rde_normal_denominator_field(Expr* f, Expr* g, RdeCtx* C,
                                        Expr** a_out, Expr** b_out,
                                        Expr** c_out, Expr** h_out) {
    Expr* den_f = rc_den(f, C);
    Expr* den_g = rc_den(g, C);
    Expr* dn = rc_normal_part(den_f, C);            /* SplitFactor normal part */
    Expr* en = rc_normal_part(den_g, C);
    expr_free(den_f); expr_free(den_g);
    Expr* p  = rc_gcd(dn, en, C);
    Expr* den_e   = rc_D(en, C);
    Expr* den_num = rc_gcd(en, den_e, C);           /* gcd(en, D en) */
    Expr* den_p   = rc_D(p, C);
    Expr* den_den = rc_gcd(p, den_p, C);            /* gcd(p, D p)   */
    expr_free(den_e); expr_free(den_p);
    Expr* h = rc_quot(den_num, den_den, C);
    expr_free(den_num); expr_free(den_den); expr_free(p);
    /* guard: en | dn h^2 (else no solution). */
    Expr* h2 = rde_mul(h, h);
    Expr* dnh2 = rde_mul(dn, h2);
    Expr* rem = rc_rem(dnh2, en, C);
    bool ok = rt_is_zero(rem) || rde_par_in(g);     /* parametric: see rde_base's guard */
    expr_free(rem); expr_free(h2); expr_free(dnh2);
    if (!ok) { expr_free(dn); expr_free(en); expr_free(h); return 0; }
    /* a = dn h ; b = dn h f - dn Dh ; c = dn h^2 g.  Normalize b, c with Together
     * (a single reduced fraction, tau-poles retained in the denominator): the
     * product dn h^2 g is a genuine POLYNOMIAL in tau over k(x) whenever h clears
     * g's normal tau-poles, but `Cancel` does NOT distribute a monomial into a sum
     * with a cancelling denominator (e.g. it leaves t^2 (2/x - 1/(t^2 x)) as-is
     * instead of (2 t^2 - 1)/x), so the downstream PolynomialQ gate then rejects a
     * genuine polynomial and the whole SPDE solve is skipped (RISCH_AUDIT_A4.md —
     * the rational-lower-field-coefficient tower RDE regression).  Together
     * distributes and reduces, giving the canonical polynomial-over-k(x) form; it
     * keeps (never clears) any special tau-pole an exponential top leaves for
     * RdeSpecialDenomExp. */
    Expr* aa = rde_mul(dn, h);
    Expr* dh = rc_D(h, C);
    Expr* dhf = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(dn), expr_copy(h), expr_copy(f) }, 3));
    Expr* dndh = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(dn), dh }, 2));         /* adopts dh */
    Expr* bbraw = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ dhf, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), dndh }, 2) }, 2);   /* adopts dhf, dndh */
    Expr* bb = rde_canon(bbraw); expr_free(bbraw);
    Expr* ccraw = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(dn), expr_copy(h), expr_copy(h), expr_copy(g) }, 4);
    Expr* cc = rde_canon(ccraw); expr_free(ccraw);
    expr_free(dn); expr_free(en);
    *a_out = aa; *b_out = bb; *c_out = cc; *h_out = h;
    return 1;
}

static int rc_param_log_deriv(Expr* f, Expr* eta, RdeCtx* C, long* m_out, Expr** z_out);
static bool rc_is_tower_rational(Expr* e, RdeCtx* C);

/* The field integrator on an AUXILIARY integrand F (an alpha of LimitedIntegrate,
 * a c/w of the log-derivative branch, the c of D H = c), in its own decision
 * scope.  g_rt_decision is write-once and read by the outer decision: a
 * certificate the integrator raises about F is a fact about F, not about the
 * outer integrand, and must not reach the outer verdict (review step 0.1: the
 * outer RDE could then give up uncertain and still read NONELEMENTARY).  F's own
 * verdict is returned in *f_nonelem for the caller's three-way answer. */
static Expr* rc_aux_integrate(Expr* F, RdeCtx* C, Expr** rem, bool* f_nonelem) {
    RtDecision saved = g_rt_decision;
    g_rt_decision = RT_DEC_UNKNOWN;
    Expr* R = rt_field_integrate(F, C->T, C->m, C->x, rem);
    if (f_nonelem) *f_nonelem = (g_rt_decision == RT_DEC_NONELEMENTARY);
    g_rt_decision = saved;
    return R;
}

/* Is R (tower-variable form, possibly with Log[...] terms) an exact
 * antiderivative of F: D_tower[R] - F == 0? */
static bool rc_is_antiderivative(Expr* R, Expr* F, RdeCtx* C) {
    Expr* dR = rt_tower_deriv(R, C->T, C->x);
    if (!dR) return false;
    Expr* d = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ dR,
        expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
    Expr* t = rt_eval1("Together", d);
    Expr* nm = t ? rt_eval1("Numerator", t) : NULL;
    Expr* ex = nm ? rt_eval1("Expand", nm) : NULL;
    bool z = ex && rt_is_zero(ex);
    if (ex) expr_free(ex);
    return z;
}

/* Every ArcTan[q] (q depending on x or the tower) as logarithms,
 * ArcTan[q] = (I/2) (Log[1 - I q] - Log[1 + I q]) + const -- an identity of
 * differential algebra over C(i).  For the integrate-in-K question only: the
 * logarithms are then decided by the structure theorem like any other. */
static Expr* rc_arctan_to_logs(Expr* R) {
    if (rt_free_of_head(R, "ArcTan") && rt_free_of_head(R, "ArcTanh")) return R;
    Expr* pq = expr_new_symbol("Risch`Private`rcAtq");
    Expr* lhs = expr_new_function(expr_new_symbol("ArcTan"), (Expr*[]){
        expr_new_function(expr_new_symbol("Pattern"), (Expr*[]){ expr_copy(pq),
            expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2) }, 1);
    Expr* iq = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        expr_new_function(expr_new_symbol("Complex"), (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2),
        expr_copy(pq) }, 2);
    Expr* rhs = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        expr_new_function(expr_new_symbol("Complex"), (Expr*[]){
            expr_new_integer(0), expr_new_function(expr_new_symbol("Rational"),
                (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) }, 2),
        expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
            expr_new_function(expr_new_symbol("Log"), (Expr*[]){
                expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_new_integer(1),
                    expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(iq) }, 2) }, 2) }, 1),
            expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1),
                expr_new_function(expr_new_symbol("Log"), (Expr*[]){
                    expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_new_integer(1), iq }, 2) }, 1) }, 2) }, 2) }, 2);
    Expr* rule = expr_new_function(expr_new_symbol("RuleDelayed"), (Expr*[]){ lhs, rhs }, 2);
    /* ArcTanh[q] = (1/2) (Log[1 + q] - Log[1 - q]): the lower integrator writes the
     * same conjugate log pair either way, depending on which route answered
     * (16 ArcTanh[1 - 2 I x] = -16 I ArcTan[I + 2 x] + const) */
    Expr* lhs2 = expr_new_function(expr_new_symbol("ArcTanh"), (Expr*[]){
        expr_new_function(expr_new_symbol("Pattern"), (Expr*[]){ expr_copy(pq),
            expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2) }, 1);
    Expr* rhs2 = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        expr_new_function(expr_new_symbol("Rational"), (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2),
        expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
            expr_new_function(expr_new_symbol("Log"), (Expr*[]){
                expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_new_integer(1), expr_copy(pq) }, 2) }, 1),
            expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1),
                expr_new_function(expr_new_symbol("Log"), (Expr*[]){
                    expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_new_integer(1),
                        expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(pq) }, 2) }, 2) }, 1) }, 2) }, 2) }, 2);
    expr_free(pq);
    Expr* rule2 = expr_new_function(expr_new_symbol("RuleDelayed"), (Expr*[]){ lhs2, rhs2 }, 2);
    Expr* rl = expr_new_function(expr_new_symbol("List"), (Expr*[]){ rule, rule2 }, 2);
    Expr* out = rt_eval2("ReplaceAll", R, rl);   /* adopts R */
    return out ? rt_eval1("Expand", out) : NULL;
}

/* Integrate F (an element of K = K_m, C's field) inside K, three ways (campaign
 * step 1, v0.6605):
 *   1  *Q owned, an antiderivative rational in x, t_0..t_m;
 *   0  PROVED: F has no antiderivative in K (decide mode certified F
 *      non-elementary, or its antiderivative carries new logarithms the
 *      structure theorem proves independent of the tower -- Kolchin-Ostrowski
 *      then puts it outside K + C);
 *  -1  undecided.
 * Pass 1 runs the integrator as it integrates (it finds the most answers; a
 * new logarithm that reduces into the tower is folded).  Pass 2, only when
 * pass 1 found nothing, runs it in decide mode for the proof.  Each pass runs
 * in its own certification scope; an undecided answer propagates the reason. */
static int rc_integrate_in_K(Expr* F, RdeCtx* C, Expr** Q) {
    *Q = NULL;
    if (rt_is_zero(F)) { *Q = expr_new_integer(0); return 1; }
    RdeCertState saved = rde_cert_begin();
    Expr* rem = NULL;
    Expr* R = rc_aux_integrate(F, C, &rem, NULL);
    bool unc1 = rde_cert_uncertain();
    const char* why1 = rde_cert_reason();
    rde_cert_end(saved);
    bool norem = !rem || rt_is_zero(rem);
    if (rem) expr_free(rem);
    if (R && norem) {
        if (rc_is_tower_rational(R, C)) { *Q = R; return 1; }
        /* an EXACT antiderivative whose leftover new logarithms are independent
         * of the tower proves "none in K" whichever mode produced it */
        bool exact = rc_is_antiderivative(R, F, C);
        int st = -1;
        R = rc_arctan_to_logs(R);
        Expr* Rk = R ? rt_reduce_new_logs_ext(R, C->T, C->m, C->x, &st) : NULL;   /* consumes R */
        if (Rk && st == 1 && rc_is_tower_rational(Rk, C)) { *Q = Rk; return 1; }
        if (Rk) expr_free(Rk);
        if (exact && st == 0 && !C->T->dependent) return 0;
    } else if (R) expr_free(R);
    /* pass 2: the decision */
    int save_mode = g_rt_decide_mode;
    g_rt_decide_mode = 1;
    saved = rde_cert_begin();
    rem = NULL;
    bool ne = false;
    Expr* R2 = rc_aux_integrate(F, C, &rem, &ne);
    bool unc2 = rde_cert_uncertain();
    const char* why2 = rde_cert_reason();
    rde_cert_end(saved);
    g_rt_decide_mode = save_mode;
    int verdict = -1;
    bool hasrem = rem && !rt_is_zero(rem);
    if (rem) expr_free(rem);
    if (!unc2) {
        if (!R2 && ne) verdict = 0;                   /* F not elementary at all */
        else if (R2 && hasrem) verdict = 0;           /* a certified remainder */
        else if (R2 && rc_is_antiderivative(R2, F, C)) {
            if (rc_is_tower_rational(R2, C)) { *Q = R2; return 1; }
            int st = -1;
            R2 = rc_arctan_to_logs(R2);
            Expr* Rk = R2 ? rt_reduce_new_logs_ext(R2, C->T, C->m, C->x, &st) : NULL;   /* consumes R2 */
            R2 = NULL;
            if (Rk && st == 1 && rc_is_tower_rational(Rk, C)) { *Q = Rk; return 1; }
            if (st == 0) verdict = 0;                 /* independent new logarithms */
            if (Rk) expr_free(Rk);
        }
    }
    if (R2) expr_free(R2);
    if (verdict < 0 && (unc2 || unc1)) rde_cert_mark(unc2 ? why2 : why1);
    return verdict;
}

/* LimitedIntegrate against the one primitive tau (Bronstein §7.2, the case the
 * degree bound of a primitive top needs): is alpha = Dz + s Dtau with z in k and
 * a constant s?  s is unique (two values would put tau in k + C).  Any
 * antiderivative of alpha that lies in K = k(tau) is z + s tau + const, so the
 * question is "does alpha have an antiderivative in K", answered by the field
 * integrator.  Returns 1 with *s owned, -1 when undecided (see
 * rc_integrate_in_K: a new logarithm in the answer is not a proof). */
static int rc_limited_integrate_prim(Expr* alpha, RdeCtx* C, Expr** s_out) {
    *s_out = NULL;
    Expr* Q = NULL;
    int ik = rc_integrate_in_K(alpha, C, &Q);
    if (ik != 1) return ik;                          /* 0: proved no solution */
    if (rc_deg(Q, C) > 1) { expr_free(Q); return -1; }   /* cannot happen; be safe */
    Expr* sv = rc_deg(Q, C) == 1 ? rc_coeff(Q, 1, C) : expr_new_integer(0);
    expr_free(Q);
    *s_out = rde_canon(sv); expr_free(sv);
    return 1;
}

/* RdeSpecialDenomExp (Bronstein §6.2, p.190): given a Dq + b q = c with a in k[tau]
 * special-free and b, c in k<tau> (tau hyperexponential, eta = Dtau/tau in k),
 * clear the special (tau-power) part of the denominator of any solution q.  Returns
 * (a tau^N, (b + n a eta) tau^N, c tau^{N-n}, tau^{-n}); r = q tau^{-n} in k[tau]
 * solves the returned polynomial equation.  1b-core omits the nb == 0 cancellation
 * refinement (the parametric logarithmic-derivative test), so a tau=0 pole arising
 * purely from cancellation is not lowered — sound (the exact identity gate in
 * rde_tower makes any miss a decline). */
static void rde_special_denom_exp(Expr* a, Expr* b, Expr* c, RdeCtx* C,
                                  Expr** a_out, Expr** b_out, Expr** c_out,
                                  Expr** hs_out) {
    Expr* eta = C->T->Dcoef[C->m];                  /* Dtau/tau in k */
    long nb = rc_val_tau(b, C);
    long nc = rc_val_tau(c, C);
    long minb = nb < 0 ? nb : 0;
    long n = nc - minb; if (n > 0) n = 0;            /* n = min(0, nc - min(0,nb)) */
    /* ParametricLogarithmicDerivative refinement (Bronstein §6.2, nb == 0 case,
     * p.190-191): when b is a tau-UNIT (nb == 0) the solution can still carry a
     * special tau = 0 pole produced by leading-coefficient CANCELLATION, which the
     * valuation formula above cannot see.  Writing q = tau^m u (u a tau-unit), the
     * t^m term of a Dq + b q is (m eta a + b) tau^m u, whose value at tau = 0 is
     * (m eta a(0) + b(0)) u(0); a genuine tau^m factor needs the resonance residue
     *   alpha = -b(0) / (a(0) eta)
     * to be a NONPOSITIVE integer, whereupon alpha + m eta = 0 is the trivial
     * logarithmic derivative Dz/z (z const) that Bronstein's ParametricLogarithmic-
     * Derivative certifies.  Lower n to that resonance.  (A spurious lowering only
     * widens the search; the exact-identity gate in rde_tower rejects any wrong q,
     * so this can never ship a wrong result.)  This closes the nested-exponential
     * Laurent solutions t_L/t_0 etc. (RISCH_AUDIT_A4.md — Bronstein Example 6.2.x). */
    if (nb == 0) {
        /* ParametricLogarithmicDerivative (Bronstein §6.2, p.190-191): with b a
         * tau-unit, a solution can still carry a tau^m pole produced by
         * cancellation, exactly when -b(0)/a(0) = m eta + Dz/z (m in Z, z in k*);
         * then n = min(n, m).  Decided by rc_param_log_deriv. */
        Expr* a0 = rt_eval2("ReplaceAll", expr_copy(a), expr_new_function(
            expr_new_symbol("Rule"), (Expr*[]){ expr_copy(C->tau), expr_new_integer(0) }, 2));
        if (a0 && !rt_is_zero(a0)) {
            Expr* b0 = rt_eval2("ReplaceAll", expr_copy(b), expr_new_function(
                expr_new_symbol("Rule"), (Expr*[]){ expr_copy(C->tau), expr_new_integer(0) }, 2));
            Expr* alpha = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), b0,
                    expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(a0), expr_new_integer(-1) }, 2) }, 3));  /* adopts b0 */
            long mm = 0; Expr* z = NULL;
            int pl = alpha ? rc_param_log_deriv(alpha, eta, C, &mm, &z) : -1;
            if (pl == 1 && mm < n) n = mm;            /* n = min(n, m) */
            if (pl < 0) rde_cert_mark("RdeSpecialDenomExp: parametric log-derivative undecided");
            if (z) expr_free(z);
            if (alpha) expr_free(alpha);
        } else rde_cert_mark("RdeSpecialDenomExp: a(0) = 0");
        if (a0) expr_free(a0);
    }
    long N = 0;
    if (-nb > N) N = -nb;
    if (n - nc > N) N = n - nc;                      /* N = max(0, -nb, n - nc) */
    Expr* tauN = rc_tau_pow(N, C);
    *a_out = rde_mul(a, tauN);                       /* a tau^N */
    /* b + n a eta, then * tau^N. */
    Expr* naeta = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(n), expr_copy(a), expr_copy(eta) }, 3));
    Expr* bsum = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(b), naeta }, 2);        /* adopts naeta */
    Expr* bt = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ bsum, expr_copy(tauN) }, 2);      /* adopts bsum */
    *b_out = rde_canon(bt); expr_free(bt);
    Expr* tauNn = rc_tau_pow(N - n, C);
    Expr* ct = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(c), tauNn }, 2);        /* adopts tauNn */
    *c_out = rde_canon(ct); expr_free(ct);
    *hs_out = rc_tau_pow(-n, C);                     /* tau^{-n}, -n >= 0 */
    expr_free(tauN);
}

/* ---- Hypertangent top (Bronstein §6.2-§6.6, delta(t) = 2), review F8 ----
 * tau = t_m with D tau = eta (tau^2 + 1), eta = T->Dcoef[m] in k (sigma = +1;
 * the hyperbolic sigma = -1 top, whose special part t^2 - 1 splits, is still a
 * marked decline).  The special irreducible is p = tau^2 + 1, Dp/p = 2 eta tau;
 * its root sqrt(-1) is adjoined as the constant I: evaluation at tau = I
 * commutes with D because Dtau vanishes there, so a condition "mod p" becomes
 * an equation over k(i). */

static bool rc_is_tan_top(RdeCtx* C) {
    return C->m >= 0 && C->T->kind[C->m] == RT_TAN && C->T->tsg[C->m] == 1;
}
/* p = tau^2 + 1 */
static Expr* rc_tan_p(RdeCtx* C) {
    return expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ rc_tau_pow(2, C), expr_new_integer(1) }, 2);
}
/* multiplicity of p in the polynomial e (0 when e == 0) */
static long rc_ord_p(Expr* e, Expr* p, RdeCtx* C) {
    if (rt_is_zero(e)) return 0;
    Expr* u = expr_copy(e);
    long k = 0;
    while (rc_deg(u, C) >= 2) {
        Expr* r = rc_rem(u, p, C);
        bool div = r && rt_is_zero(r);
        if (r) expr_free(r);
        if (!div) break;
        Expr* q = rc_quot(u, p, C);
        expr_free(u); u = q; k++;
        if (!u) break;
    }
    if (u) expr_free(u);
    return k;
}
/* p-adic valuation of a field element of k(tau) */
static long rc_val_p(Expr* e, Expr* p, RdeCtx* C) {
    Expr* num = rde_num(e);
    Expr* den = rde_den(e);
    long v = rc_ord_p(num, p, C) - rc_ord_p(den, p, C);
    expr_free(num); expr_free(den);
    return v;
}
/* e at tau = I, canonical */
static Expr* rc_at_i(Expr* e, RdeCtx* C) {
    Expr* r = rt_eval2("ReplaceAll", expr_copy(e), expr_new_function(
        expr_new_symbol("Rule"), (Expr*[]){ expr_copy(C->tau), expr_new_symbol("I") }, 2));
    if (!r) return NULL;
    Expr* c = rde_canon(r); expr_free(r);
    return c;
}
/* Complex conjugate of an element of k(i) (every tower variable is real):
 * Complex[a, b] -> Complex[a, -b], canonical.  Owned. */
static Expr* rc_conj(Expr* e) {
    Expr* pa = expr_new_symbol("Risch`Private`rdeCja");
    Expr* pb = expr_new_symbol("Risch`Private`rdeCjb");
    Expr* lhs = expr_new_function(expr_new_symbol("Complex"), (Expr*[]){
        expr_new_function(expr_new_symbol("Pattern"),
            (Expr*[]){ expr_copy(pa), expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2),
        expr_new_function(expr_new_symbol("Pattern"),
            (Expr*[]){ expr_copy(pb), expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2) }, 2);
    Expr* rhs = expr_new_function(expr_new_symbol("Complex"), (Expr*[]){ pa,
        expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), pb }, 2) }, 2);
    Expr* rule = expr_new_function(expr_new_symbol("RuleDelayed"), (Expr*[]){ lhs, rhs }, 2);
    Expr* r = rt_eval2("ReplaceAll", expr_copy(e), rule);
    if (!r) return NULL;
    Expr* c = rde_canon(r); expr_free(r);
    return c;
}
/* Real and imaginary parts of s in k(i): 1 with both owned and free of I, 0
 * otherwise. */
static int rc_re_im(Expr* s, Expr** re, Expr** im) {
    *re = *im = NULL;
    Expr* sb = rc_conj(s);
    if (!sb) return 0;
    Expr* sum = rde_add(s, sb);
    Expr* dif = rde_sub(s, sb);
    expr_free(sb);
    Expr* half = expr_new_function(expr_new_symbol("Rational"),
        (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2);
    Expr* r = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ sum, expr_copy(half) }, 2));                        /* adopts sum */
    Expr* i2 = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ half,
        expr_new_function(expr_new_symbol("Complex"),
            (Expr*[]){ expr_new_integer(0), expr_new_integer(-1) }, 2) }, 2);  /* -I/2 */
    Expr* m = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ dif, i2 }, 2));                                     /* (s - sb)/(2I) */
    bool ok = r && m && rt_free_of_head(r, "Complex") && rt_free_of_head(m, "Complex");
    if (!ok) { if (r) expr_free(r); if (m) expr_free(m); return 0; }
    *re = rde_canon(r); *im = rde_canon(m);
    expr_free(r); expr_free(m);
    return 1;
}

/* RdeSpecialDenomTan (Bronstein §6.2, p.192): as rde_special_denom_exp, for
 * the special irreducible p = tau^2 + 1.  With q = p^n r, a Dq + b q =
 * p^n (a Dr + (b + 2 n eta tau a) r), so r solves the returned equation
 * (a p^N, (b + 2 n eta tau a) p^N, c p^(N - n)) and h_s = p^(-n).  The bound
 * n = min(0, nu_p(c) - min(0, nu_p(b))) is refined when b is a p-unit: a pole
 * of order -m can then come from cancellation, exactly when
 * -b(i)/a(i) = m (2 eta i) + Dv/v for some v in k(i)* (evaluate the order-m
 * term at tau = i, where Dtau = 0), decided by rc_param_log_deriv over k(i). */
static void rde_special_denom_tan(Expr* a, Expr* b, Expr* c, RdeCtx* C,
                                  Expr** a_out, Expr** b_out, Expr** c_out,
                                  Expr** hs_out) {
    Expr* eta = C->T->Dcoef[C->m];
    Expr* p = rc_tan_p(C);
    long nb = rc_val_p(b, p, C);
    long nc = rc_val_p(c, p, C);
    long minb = nb < 0 ? nb : 0;
    long n = nc - minb; if (n > 0) n = 0;
    if (nb == 0) {
        Expr* ai = rc_at_i(a, C);
        Expr* bi = rc_at_i(b, C);
        if (ai && bi && !rt_is_zero(ai)) {
            Expr* f = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), expr_copy(bi),
                    expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(ai), expr_new_integer(-1) }, 2) }, 3));
            Expr* e2 = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(2), expr_new_symbol("I"), expr_copy(eta) }, 3));
            long mm = 0; Expr* z = NULL;
            int pl = (f && e2) ? rc_param_log_deriv(f, e2, C, &mm, &z) : -1;
            if (pl == 1 && mm < n) n = mm;
            if (pl < 0) rde_cert_mark("RdeSpecialDenomTan: parametric log-derivative undecided");
            if (z) expr_free(z);
            if (f) expr_free(f);
            if (e2) expr_free(e2);
        } else rde_cert_mark("RdeSpecialDenomTan: a(i) = 0");
        if (ai) expr_free(ai);
        if (bi) expr_free(bi);
    }
    long N = 0;
    if (-nb > N) N = -nb;
    if (n - nc > N) N = n - nc;
    Expr* pN = expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_copy(p), expr_new_integer(N) }, 2);
    *a_out = rde_mul(a, pN);
    Expr* nterm = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        expr_new_integer(2 * n), expr_copy(eta), expr_copy(C->tau), expr_copy(a) }, 4);
    Expr* bt = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(b), nterm }, 2),
        expr_copy(pN) }, 2);
    *b_out = rde_canon(bt); expr_free(bt);
    Expr* ct = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_copy(c),
        expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(p), expr_new_integer(N - n) }, 2) }, 2);
    *c_out = rde_canon(ct); expr_free(ct);
    *hs_out = rt_eval1("Expand", expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_copy(p), expr_new_integer(-n) }, 2));
    expr_free(pN); expr_free(p);
}

/* RdeBoundDegreeNonLinear (Bronstein §6.3, p.200) for delta(tau) = 2,
 * lambda = lc(Dtau) = eta: a nonzero q of degree n gives a Dq of degree
 * da + n + 1 with leading coefficient n eta lc(a) lc(q), so deg(a Dq + b q) =
 * max(db, da + 1) + n unless db = da + 1 and n = -lc(b)/(eta lc(a)). */
static long rde_bound_degree_tan(Expr* a, Expr* b, Expr* c, RdeCtx* C) {
    long da = rc_deg(a, C), db = rc_deg(b, C), dc = rc_deg(c, C);
    long mx = db > da + 1 ? db : da + 1;
    long n = dc - mx; if (n < 0) n = 0;
    if (db == da + 1) {
        Expr* lca = rc_lc(a, C);
        Expr* lcb = rc_lc(b, C);
        Expr* den = rde_mul(lca, C->T->Dcoef[C->m]);
        Expr* q = rde_divq(lcb, den);
        Expr* m = q ? rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), q }, 2)) : NULL;
        if (m && m->type == EXPR_INTEGER && (long)m->data.integer > n) n = (long)m->data.integer;
        if (m && m->type == EXPR_BIGINT) rde_cert_mark("RdeBoundDegreeTan: bound out of range");
        if (m) expr_free(m);
        expr_free(lca); expr_free(lcb); expr_free(den);
    }
    return n;
}

/* c - D p - b p, canonical */
static Expr* rc_rde_residual(Expr* c, Expr* p, Expr* b, RdeCtx* C) {
    Expr* dp = rc_D(p, C);
    Expr* bp = rde_mul(b, p);
    Expr* t1 = rde_sub(c, dp);
    Expr* t2 = rde_sub(t1, bp);
    Expr* r = rde_canon(t2);
    expr_free(dp); expr_free(bp); expr_free(t1); expr_free(t2);
    return r;
}
/* lc-term coef * tau^m */
static Expr* rc_term(Expr* coef, long m, RdeCtx* C) {
    return rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ coef, rc_tau_pow(m, C) }, 2));        /* adopts coef */
}
static Expr* rde_cancel_tan(Expr* b0, Expr* c, RdeCtx* C, long n);

/* PolyRischDENoCancel2 (Bronstein §6.5, p.208), delta = 2, b in k*: a q of
 * degree d >= 1 gives D q + b q of degree d + 1 with leading coefficient
 * d eta lc(q), so every coefficient down to degree 1 is forced; the degree-0
 * remainder is the Risch DE D z + b z = c over k (recursion).  Each step is
 * forced, so a failure is a proof (the lower-field recursion marks its own
 * uncertainty). */
static Expr* rde_nocancel2_tan(Expr* b, Expr* c, RdeCtx* C, long n) {
    Expr* eta = C->T->Dcoef[C->m];
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        long dc = rc_deg(cc, C);
        if (dc == 1 || n < 0) { fail = true; break; }
        if (dc == 0) {
            RdeCtx Clo = { C->T, C->m - 1, C->x, (C->m - 1 >= 0) ? C->T->t[C->m - 1] : C->x };
            Expr* z = rde_tower(b, cc, &Clo);
            if (!z) { fail = true; break; }
            Expr* nq = rde_add(q, z); expr_free(q); expr_free(z); q = nq;
            break;
        }
        long m = dc - 1;
        if (m > n) { fail = true; break; }
        Expr* lcc = rc_lc(cc, C);
        Expr* me = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(m), expr_copy(eta) }, 2));
        Expr* co = rde_divq(lcc, me);
        expr_free(lcc); expr_free(me);
        Expr* p = rc_term(co, m, C);
        Expr* nq = rde_add(q, p); expr_free(q); q = nq;
        Expr* r = rc_rde_residual(cc, p, b, C);
        expr_free(cc); expr_free(p); cc = r;
        n = m - 1;
    }
    expr_free(cc);
    if (fail) { expr_free(q); return NULL; }
    return q;
}

/* PolyRischDENoCancel3 (Bronstein §6.5, p.209), delta = 2, deg b = 1: the
 * degree-(m+1) coefficient of D q + b q for q of degree m is
 * (m eta + lc(b)) lc(q), forced unless u = m eta + lc(b) = 0, i.e. at the one
 * resonant degree M = -lc(b)/eta; there the remainder is the cancellation
 * equation D h + (b0 - M eta tau) h = c, deg h <= M (PolyRischDECancelTan). */
static Expr* rde_nocancel3_tan(Expr* b, Expr* c, RdeCtx* C, long n) {
    Expr* eta = C->T->Dcoef[C->m];
    Expr* b1 = rc_coeff(b, 1, C);
    long M = -1;
    {
        Expr* r = rde_divq(b1, eta);
        Expr* mr = r ? rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), r }, 2)) : NULL;
        if (mr && mr->type == EXPR_INTEGER && mr->data.integer >= 0) M = (long)mr->data.integer;
        if (mr && mr->type == EXPR_BIGINT) rde_cert_mark("PolyRischDENoCancel3: resonance out of range");
        if (mr) expr_free(mr);
    }
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        /* a resonance above the bound cannot occur; below the top of c every
         * degree above M is forced to vanish, so the solution starts at M */
        long m = rc_deg(cc, C) - 1;
        if (M >= 0 && M <= n && M > m) m = M;
        if (n < 0 || m < 0 || m > n) { fail = true; break; }
        if (m == M) {
            Expr* b0 = rc_coeff(b, 0, C);
            Expr* h = rde_cancel_tan(b0, cc, C, m);
            expr_free(b0);
            if (!h) { fail = true; break; }
            Expr* nq = rde_add(q, h); expr_free(q); expr_free(h); q = nq;
            expr_free(cc); cc = expr_new_integer(0);
            break;
        }
        Expr* u = rt_eval1("Together", expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(m), expr_copy(eta) }, 2), expr_copy(b1) }, 2));
        Expr* cm = rc_coeff(cc, m + 1, C);
        Expr* co = rde_divq(cm, u);
        expr_free(cm); expr_free(u);
        Expr* p = rc_term(co, m, C);
        Expr* nq = rde_add(q, p); expr_free(q); q = nq;
        Expr* r = rc_rde_residual(cc, p, b, C);
        expr_free(cc); expr_free(p); cc = r;
        n = m - 1;
    }
    expr_free(cc); expr_free(b1);
    if (fail) { expr_free(q); return NULL; }
    return q;
}

/* PolyRischDECancelTan (Bronstein §6.6, p.215): D q + (b0 - n eta tau) q = c,
 * deg q <= n, b0 in k.  At tau = i the equation is the Risch DE
 * D s + (b0 - n eta i) s = c(i) over k(i) for s = q(i) (Bronstein's
 * CoupledDESystem, by the trivial route of risch_coupled.h); then
 * q0 = Im(s) tau + Re(s) agrees with q mod p, p divides c - D q0 - b q0, and
 * q = q0 + p h with D h + (b0 - (n - 2) eta tau) h = (c - D q0 - b q0)/p. */
static Expr* rde_cancel_tan(Expr* b0, Expr* c, RdeCtx* C, long n) {
    if (rt_is_zero(c)) return expr_new_integer(0);
    if (n < 0) return NULL;
    RdeCtx Clo = { C->T, C->m - 1, C->x, (C->m - 1 >= 0) ? C->T->t[C->m - 1] : C->x };
    if (n == 0) {
        if (rc_deg(c, C) > 0) return NULL;
        return rde_tower(b0, c, &Clo);
    }
    Expr* eta = C->T->Dcoef[C->m];
    Expr* bi = rt_eval1("Together", expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(b0), expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-n), expr_copy(eta), expr_new_symbol("I") }, 3) }, 2));
    Expr* ci = rc_at_i(c, C);
    Expr* s = (bi && ci) ? rde_tower(bi, ci, &Clo) : NULL;
    /* s is unique unless the homogeneous equation has a solution in k(i)* */
    int hom = -1;
    if (bi) {
        Expr* mb = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(bi) }, 2));
        long mm = 0; Expr* z = NULL;
        hom = mb ? rc_param_log_deriv(mb, NULL, C, &mm, &z) : -1;
        if (z) expr_free(z);
        if (mb) expr_free(mb);
    }
    if (bi) expr_free(bi);
    if (ci) expr_free(ci);
    Expr* q = NULL;
    if (s) {
        Expr *re, *im;
        if (rc_re_im(s, &re, &im)) {
            Expr* imt = rde_mul(im, C->tau);
            Expr* q0 = rde_add(imt, re);                /* Im(s) tau + Re(s) */
            expr_free(imt); expr_free(re); expr_free(im);
            Expr* bt = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(b0), expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-n), expr_copy(eta), expr_copy(C->tau) }, 3) }, 2));
            Expr* r = rc_rde_residual(c, q0, bt, C);
            expr_free(bt);
            Expr* p = rc_tan_p(C);
            Expr* rr = rc_rem(r, p, C);
            if (rr && rt_is_zero(rr)) {
                Expr* r2 = rc_quot(r, p, C);
                Expr* h = r2 ? rde_cancel_tan(b0, r2, C, n - 2) : NULL;
                if (h) {
                    Expr* ph = rde_mul(p, h);
                    q = rde_add(q0, ph);
                    expr_free(ph); expr_free(h);
                }
                if (r2) expr_free(r2);
            } else rde_cert_mark("PolyRischDECancelTan: residual not divisible by tau^2 + 1");
            if (rr) expr_free(rr);
            expr_free(p); expr_free(r); expr_free(q0);
        } else rde_cert_mark("PolyRischDECancelTan: q(i) not split into real parts");
        expr_free(s);
    }
    if (!q && hom != 0)
        rde_cert_mark("PolyRischDECancelTan: homogeneous solution over k(i) not excluded");
    return q;
}

/* RdeBoundDegreeExp (Bronstein §6.3, p.200): upper bound on deg_tau(q) for
 * a Dq + b q = c with tau hyperexponential (deriv-preserving, delta(tau)=1). */
static long rde_bound_degree_exp(Expr* a, Expr* b, Expr* c, RdeCtx* C) {
    long da = rc_deg(a, C), db = rc_deg(b, C), dc = rc_deg(c, C);
    long mx = db > da ? db : da;
    long n = dc - mx; if (n < 0) n = 0;
    /* da == db cancellation refinement (Bronstein §6.3, Lemma 6.3.4(ii)): when
     * deg(a) == deg(b) the leading terms of a Dq and b q can CANCEL, letting deg(q)
     * exceed the naive bound.  Cancellation occurs at the resonance m where
     * -lc(b)/lc(a) = m eta + Dz/z (parametric logarithmic derivative); in the pure
     * case (z const) m = -lc(b)/(lc(a) eta).  Raise the bound to max(n, m).  A
     * too-large bound only widens the SPDE search — the exact-identity gate in
     * rde_tower rejects any wrong q, so this can never ship a wrong result. */
    if (da == db && da >= 0) {
        /* deg a == deg b: the leading terms of a Dq and b q cancel at degree m
         * exactly when -lc(b)/lc(a) = m eta + Dz/z (m in Z, z in k*), Bronstein
         * Lemma 6.3.4(ii); then n = max(n, m).  Decided by rc_param_log_deriv. */
        Expr* eta = C->T->Dcoef[C->m];
        Expr* lca = rc_lc(a, C);
        Expr* lcb = rc_lc(b, C);
        if (lca && !rt_is_zero(lca)) {
            Expr* alpha = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), expr_copy(lcb),
                    expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(lca), expr_new_integer(-1) }, 2) }, 3));
            long mm = 0; Expr* z = NULL;
            int pl = alpha ? rc_param_log_deriv(alpha, eta, C, &mm, &z) : -1;
            if (pl == 1 && mm > n) n = mm;            /* n = max(n, m) */
            if (pl < 0) rde_cert_mark("RdeBoundDegreeExp: parametric log-derivative undecided");
            if (z) expr_free(z);
            if (alpha) expr_free(alpha);
        } else rde_cert_mark("RdeBoundDegreeExp: lc(a) = 0");
        if (lca) expr_free(lca);
        if (lcb) expr_free(lcb);
    }
    return n;
}

/* RdeBoundDegree for a primitive (RT_LOG) / base monomial (Bronstein §6.3,
 * RdeBoundDegreeBase/Prim): an upper bound on deg_tau(q) for a Dq + b q = c,
 * Dtau = eta in k.  Derivation (leading tau-coefficients, q = Sum q_i tau^i):
 *   db > da:      deg(lhs) = db + n exactly          -> n <= dc - db;
 *   db = da, no cancellation: likewise                -> n <= dc - db;
 *   db < da - 1:  deg(lhs) >= da + n - 1             -> n <= dc - da + 1;
 *   db = da - 1:  the tau^(da+n-1) coefficient with q_n constant is
 *                 lc(a)(Dq_{n-1} + n q_n eta) + lc(b) q_n, which vanishes iff
 *                 alpha = -lc(b)/lc(a) = Dz + n eta (z in k): n is the unique
 *                 constant s of LimitedIntegrate(alpha, eta)  -> n <= max(.., s);
 *   db = da:      lc(a) Dq_n + lc(b) q_n vanishes iff alpha = Dz/z, q_n = c z;
 *                 then deg(lhs) <= da + n - 1, which may equal dc, so
 *                 n = dc - da + 1 needs nothing more; a larger n needs the next
 *                 coefficient to vanish too: -beta = D(u/c) + n eta with
 *                 beta = coeff_{da-1}(a Dz + b z)/(z lc(a)), so
 *                 n <= max(dc - da + 1, s), s = LimitedIntegrate(-beta, eta).
 * At the base field (m < 0, tau = x, eta = 1, k = Q) z is a constant: the
 * db = da - 1 case is alpha = s, an integer test, and db = da cannot cancel. */
static long rde_bound_degree_prim(Expr* a, Expr* b, Expr* c, RdeCtx* C) {
    long da = rc_deg(a, C), db = rc_deg(b, C), dc = rc_deg(c, C);
    long n;
    if (db >= da) n = dc - db;      /* no cancellation: deg(lhs) = db + n */
    else n = dc - da + 1;
    if (n < 0) n = 0;
    if (C->m < 0) {
        /* base field: the historical integer test is exact (eta = 1, z constant) */
        if (db == da - 1 && da >= 1) {
            Expr* lcb = rc_lc(b, C); Expr* lca = rc_lc(a, C);
            Expr* mm = rde_divq(lcb, lca);
            Expr* mmn = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), mm }, 2));
            if (mmn && mmn->type == EXPR_INTEGER && (long)mmn->data.integer > n)
                n = (long)mmn->data.integer;
            if (mmn) expr_free(mmn);
            expr_free(lcb); expr_free(lca);
        }
        return n;
    }
    Expr* lca = rc_lc(a, C);
    Expr* lcb = rc_lc(b, C);
    if (db == da - 1 && da >= 1) {
        Expr* alpha = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(lcb),
                expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(lca), expr_new_integer(-1) }, 2) }, 3));
        Expr* sv = NULL;
        int li = alpha ? rc_limited_integrate_prim(alpha, C, &sv) : -1;
        if (li == 1) {
            if (sv->type == EXPR_INTEGER && (long)sv->data.integer > n) n = (long)sv->data.integer;
            if (sv->type == EXPR_BIGINT) rde_cert_mark("RdeBoundDegreePrim: limited integration out of range");
        } else if (li < 0) rde_cert_mark("RdeBoundDegreePrim: limited integration undecided");
        if (sv) expr_free(sv);
        if (alpha) expr_free(alpha);
    } else if (db == da) {
        Expr* alpha = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(lcb),
                expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(lca), expr_new_integer(-1) }, 2) }, 3));
        long m0 = 0; Expr* z = NULL;
        int pl = alpha ? rc_param_log_deriv(alpha, NULL, C, &m0, &z) : -1;
        if (pl < 0) rde_cert_mark("RdeBoundDegreePrim: log-derivative undecided");
        if (pl != 0 && dc - da + 1 > n) n = dc - da + 1;   /* top coefficient may cancel */
        if (pl == 1) {
            /* beta = coeff_{da-1}(a Dz + b z) / (z lc(a)) */
            Expr* Dz = rc_D(z, C);
            Expr* t1 = rde_mul(a, Dz);
            Expr* t2 = rde_mul(b, z);
            Expr* sum = rde_add(t1, t2);
            Expr* cd = rc_coeff(sum, da - 1, C);
            Expr* den = rde_mul(z, lca);
            Expr* beta = rde_divq(cd, den);
            Expr* mbeta = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), beta }, 2));   /* adopts beta */
            Expr* sv = NULL;
            int li = mbeta ? rc_limited_integrate_prim(mbeta, C, &sv) : -1;
            if (li == 1) {
                if (sv->type == EXPR_INTEGER && (long)sv->data.integer > n) n = (long)sv->data.integer;
                if (sv->type == EXPR_BIGINT) rde_cert_mark("RdeBoundDegreePrim: limited integration out of range");
            } else if (li < 0) rde_cert_mark("RdeBoundDegreePrim: limited integration undecided");
            if (sv) expr_free(sv);
            if (mbeta) expr_free(mbeta);
            expr_free(Dz); expr_free(t1); expr_free(t2); expr_free(sum);
            expr_free(cd); expr_free(den);
        }
        if (z) expr_free(z);
        if (alpha) expr_free(alpha);
    }
    expr_free(lca); expr_free(lcb);
    return n;
}

/* SPDE over k(x)[tau] (Bronstein §6.4 Thm 6.4.1, p.203): the field port of
 * rde_spde, degree-reducing recursion via one ExtendedEuclidean per level.
 * Returns 1 + fills *out, or 0 for "no solution of degree <= n". */
int rde_spde_field(Expr* a, Expr* b, Expr* c, RdeCtx* C, long n, RdeSpde* out) {
    out->b = out->c = out->alpha = out->beta = NULL; out->m = 0;
    if (n < 0) {
        if (rt_is_zero(c) || rde_par_in(c)) {       /* parametric: q = 0, see rde_spde */
            out->b = expr_new_integer(0); out->c = expr_new_integer(0);
            out->m = 0; out->alpha = expr_new_integer(0); out->beta = expr_new_integer(0);
            return 1;
        }
        return 0;
    }
    Expr* g = rc_gcd(a, b, C);
    Expr* rem = rc_rem(c, g, C);
    bool gdivc = rt_is_zero(rem) || rde_par_in(rem);
    expr_free(rem);
    if (!gdivc) { expr_free(g); return 0; }
    Expr* a1 = rc_quot(a, g, C);
    Expr* b1 = rc_quot(b, g, C);
    Expr* c1 = rc_quot(c, g, C);
    expr_free(g);
    long da = rc_deg(a1, C);
    if (da == 0) {                                  /* a in k*: base of recursion */
        out->b = rde_divq(b1, a1);
        out->c = rde_divq(c1, a1);
        out->m = n;
        out->alpha = expr_new_integer(1);
        out->beta = expr_new_integer(0);
        expr_free(a1); expr_free(b1); expr_free(c1);
        return 1;
    }
    Expr* eg = rc_xgcd(b1, a1, C);
    Expr* s = NULL;
    if (eg && eg->type == EXPR_FUNCTION && eg->data.function.arg_count == 2) {
        Expr* cof = eg->data.function.args[1];
        if (cof && cof->type == EXPR_FUNCTION && cof->data.function.arg_count == 2)
            s = cof->data.function.args[0];
    }
    if (!s) {
        rde_cert_mark("SPDE(field): extended gcd unavailable");
        if (eg) expr_free(eg);
        expr_free(a1); expr_free(b1); expr_free(c1);
        return 0;
    }
    Expr* sc = rde_mul(s, c1);
    Expr* r = rc_rem(sc, a1, C);
    expr_free(sc);
    if (eg) expr_free(eg);
    Expr* b1r = rde_mul(b1, r);
    Expr* c1_b1r = rde_sub(c1, b1r);
    expr_free(b1r);
    Expr* z = rc_quot(c1_b1r, a1, C);
    expr_free(c1_b1r);
    Expr* da1 = rc_D(a1, C);
    Expr* b2 = rde_add(b1, da1);
    expr_free(da1);
    Expr* dr = rc_D(r, C);
    Expr* c2 = rde_sub(z, dr);
    expr_free(dr); expr_free(z);
    RdeSpde sub;
    int rc = rde_spde_field(a1, b2, c2, C, n - da, &sub);
    expr_free(b1); expr_free(c1); expr_free(b2); expr_free(c2);
    if (!rc) { expr_free(a1); expr_free(r); return 0; }
    out->b = sub.b; out->c = sub.c; out->m = sub.m;
    out->alpha = rde_mul(a1, sub.alpha);
    Expr* a1beta = rde_mul(a1, sub.beta);
    out->beta = rde_add(a1beta, r);
    expr_free(a1beta);
    expr_free(sub.alpha); expr_free(sub.beta);
    expr_free(a1); expr_free(r);
    return 1;
}

/* PolyRischDENoCancel over k(x)[tau], b != 0, deg_tau(b) >= 1 (the
 * non-cancellation branch, Bronstein §6.5 p.208): the leading terms of Dq and bq
 * never cancel, so q is built top-down one tau-monomial per pass.  Returns q
 * (owned) or NULL for "no solution".  Coefficient arithmetic lc(c)/lc(b) is exact
 * field division in k(x); no sub-Risch-DE is needed (that is the cancellation
 * branch, a later increment). */
Expr* rde_polyrischde_nocancel_field(Expr* b, Expr* c, RdeCtx* C, long n) {
    long db = rc_deg(b, C);
    Expr* lcb = rc_lc(b, C);
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        long dc = rc_deg(cc, C);
        long m = dc - db;
        if (n < 0 || m < 0 || m > n) {
            if (rde_par_in(cc)) { cc = rde_par_strip_top(cc, C->tau); continue; }
            fail = true; break;
        }
        Expr* lcc = rc_lc(cc, C);
        Expr* coeff = rde_divq(lcc, lcb);           /* lc(c)/lc(b) in k(x) */
        expr_free(lcc);
        Expr* p = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ coeff, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(C->tau), expr_new_integer(m) }, 2) }, 2));  /* adopts coeff */
        Expr* nq = rde_add(q, p);
        expr_free(q); q = nq;
        Expr* dp = rc_D(p, C);
        Expr* bp = rde_mul(b, p);
        Expr* t1 = rde_sub(cc, dp);
        Expr* t2 = rde_sub(t1, bp);
        /* Canonicalize: Expand alone leaves a cancelled leading coefficient as
         * a sum of k(x)-fractions that is zero only after Together
         * (x^2 t^2/(1+x)^2 + 2x t^2/(1+x)^2 + t^2/(1+x)^2 - t^2), so the
         * CoefficientList degree stays put and the loop reported a false "no
         * solution" (y = E^x + 1/(1+x) for f = E^x). */
        expr_free(cc); cc = rde_canon(t2); expr_free(t2);
        expr_free(dp); expr_free(bp); expr_free(t1); expr_free(p);
        n = m - 1;
    }
    expr_free(cc); expr_free(lcb);
    if (fail) { expr_free(q); return NULL; }
    return q;
}

/* PolyRischDECancel{Prim,Exp} (Bronstein §6.6, p.212/213): the CANCELLATION case
 * b in k* (deg_tau(b) == 0) of the polynomial Risch DE Dq + b q = c, c in k[tau],
 * deg(q) <= n.  Builds q top-down: at each degree m the coefficient s solves the
 * lower-field Risch DE  D[s] + (b + m*eta) s = lc_tau(c)  (eta = Dtau/tau for an
 * exponential tau, absent for a primitive tau) via a RECURSIVE rde_tower over
 * K_{m-1} — the genuine one-level descent of the tower recursion.  Returns q
 * (owned) or NULL.  The b = Dz/z [+ m eta] antidifferentiation branch (a
 * higher-degree solution when b is a logarithmic derivative) needs LimitedIntegrate
 * (Gap 2) and is omitted: sound, since the main loop is exact by construction (each
 * pass cancels the top tau-term) and rde_tower's identity gate declines any miss. */
/* Is e a rational function of the whole tower {x, t_0..t_m} (no NEW Log/Exp of a
 * tower variable)?  Used to accept an antiderivative only when it is a genuine
 * element of K_m (a new logarithm would mean the Risch DE has no rational
 * solution). */
static bool rc_is_tower_rational(Expr* e, RdeCtx* C) {
    size_t nv = (size_t)C->m + 2;                   /* x, t_0..t_m */
    Expr** vv = malloc(nv * sizeof(Expr*));
    vv[0] = expr_copy(C->x);
    for (long j = 0; j <= C->m; j++) vv[j + 1] = expr_copy(C->T->t[j]);
    Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, nv);
    free(vv);
    Expr* tog = rt_eval1("Together", expr_copy(e));
    Expr* num = tog ? rt_eval1("Numerator", expr_copy(tog)) : NULL;
    Expr* den = tog ? rt_eval1("Denominator", expr_copy(tog)) : NULL;
    if (tog) expr_free(tog);
    Expr* pqn = num ? rt_eval2("PolynomialQ", num, expr_copy(vlist)) : NULL;
    Expr* pqd = den ? rt_eval2("PolynomialQ", den, expr_copy(vlist)) : NULL;
    bool ok = pqn && pqd
        && pqn->type == EXPR_SYMBOL && pqn->data.symbol.name == intern_symbol("True")
        && pqd->type == EXPR_SYMBOL && pqd->data.symbol.name == intern_symbol("True");
    if (pqn) expr_free(pqn);
    if (pqd) expr_free(pqd);
    expr_free(vlist);
    return ok;
}

/* The log-derivative branch of PolyRischDECancel{Prim,Exp} (Bronstein §6.6):
 * -b = Dw/w with w = z tau^mm (z in k*, mm in Z; mm = 0 for a primitive tau).
 * Then D q + b q = w D(q/w), so q = w Int(c/w) + lambda w.  The integral is
 * taken in K = k(tau); q must be a polynomial of degree <= n, and lambda (a
 * constant) is fixed from a coefficient q must not have.  Returns 1 with *q
 * owned, 0 (proved: no polynomial solution of degree <= n), -1 undecided. */
static int rde_cancel_logderiv(Expr* c, Expr* z, long mm, RdeCtx* C, long n, Expr** q_out) {
    *q_out = NULL;
    Expr* w = mm == 0 ? expr_copy(z)
            : rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){ expr_copy(z), expr_new_function(expr_new_symbol("Power"),
                      (Expr*[]){ expr_copy(C->tau), expr_new_integer(mm) }, 2) }, 2));
    Expr* F = rde_divq(c, w);
    Expr* Q = NULL;
    int ok = rc_integrate_in_K(F, C, &Q);
    expr_free(F);
    if (ok != 1) { expr_free(w); return ok; }      /* 0: proved, no antiderivative in K */
    Expr* q = rde_mul(w, Q);
    expr_free(Q);
    {   Expr* qc = rde_canon(q); expr_free(q); q = qc; }
    /* the homogeneous solution lambda w may repair a tau-pole (mm < 0) or the
     * degree (deg w = deg q > n): fix lambda from that one coefficient */
    long dq = rc_ispoly(q, C) ? rc_deg(q, C) : -1;
    if (dq < 0 || dq > n) {
        long at = (dq < 0) ? mm : dq;               /* the offending tau power */
        /* the coefficient of tau^at: a polynomial coefficient, or (a tau-pole of
           order -mm, mm < 0) the Laurent coefficient (q tau^-mm) at tau = 0 */
        Expr* qn;
        if (rc_ispoly(q, C)) qn = rc_coeff(q, at, C);
        else {
            Expr* sh = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(q), expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(C->tau), expr_new_integer(-at) }, 2) }, 2));
            qn = sh ? rt_eval2("ReplaceAll", sh, expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_copy(C->tau), expr_new_integer(0) }, 2)) : NULL;
        }
        Expr* lam = qn ? rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), qn, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(z), expr_new_integer(-1) }, 2) }, 3)) : NULL;
        /* lambda must be a constant */
        bool constq = lam != NULL;
        size_t nv = (size_t)C->m + 2;
        for (size_t j = 0; constq && j < nv; j++) {
            Expr* v = j == 0 ? C->x : C->T->t[j - 1];
            if (!rt_free_of_x(lam, v)) constq = false;
        }
        if (constq) {
            Expr* lw = rde_mul(lam, w);
            Expr* q2 = rde_add(q, lw);
            expr_free(lw); expr_free(q); q = q2;
        }
        if (lam) expr_free(lam);
    }
    expr_free(w);
    if (!rc_ispoly(q, C) || rc_deg(q, C) > n) { expr_free(q); return 0; }
    *q_out = q;
    return 1;
}

/* PolyRischDECancel in parametric mode (campaign step 3a): c is linear in the
 * unknowns, b in k.  The pass runs every degree m = n .. 0 and solves the
 * coefficient equation  D s + (b + m eta) s = coeff_m(cc)  by a recursive
 * (parametric) rde_tower over K_{m-1}; nothing on it is nonlinear.
 *   Exponential top: D(s tau^m) = (Ds + m eta s) tau^m, so the coefficient
 *   equations decouple and a homogeneous solution at a resonant m never
 *   matters (take it 0): the pass is complete and its failures are proofs.
 *   Primitive top: D(s tau^m) = Ds tau^m + m eta s tau^(m-1) couples the
 *   degrees.  When -b = Dz/z (rc_param_log_deriv; b = 0 gives z = 1) each s is
 *   only determined up to lambda z, and lambda moves the next coefficient
 *   ([SI] §5.8's LimitedIntegrate), so lambda_m z is added with lambda_m a fresh
 *   unknown and carried down.  -b proved not a log-derivative: s is unique.
 *   Undecided: marked.
 * A part of cc at degree m that s cannot remove is left for the residual
 * constraint.  Returns q (owned) or NULL. */
static Expr* rde_polyrischde_cancel_par(Expr* b, Expr* c, RdeCtx* C, long n) {
    RtKind kind = C->T->kind[C->m];
    Expr* eta = (kind == RT_EXP) ? C->T->Dcoef[C->m] : NULL;
    Expr* z = NULL;                                  /* homogeneous solution (primitive top) */
    if (kind != RT_EXP) {
        Expr* mb = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(b) }, 2));
        long mm = 0;
        int br = mb ? rc_param_log_deriv(mb, NULL, C, &mm, &z) : -1;
        if (mb) expr_free(mb);
        if (br < 0 || (br == 1 && mm != 0)) {
            if (z) expr_free(z);
            rde_cert_mark("parametric: PolyRischDECancel log-derivative test undecided");
            return NULL;
        }
        if (br == 0 && z) { expr_free(z); z = NULL; }
    }
    long dc = rc_deg(c, C);
    if (dc > n) n = dc;          /* a too-small bound only drops terms to the residual */
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    RdeCtx Clo = { C->T, C->m - 1, C->x, (C->m - 1 >= 0) ? C->T->t[C->m - 1] : C->x };
    for (long m = n; m >= 0; m--) {
        Expr* cm = (rc_deg(cc, C) >= m) ? rc_coeff(cc, m, C) : expr_new_integer(0);
        Expr* coef = (kind == RT_EXP)
            ? rt_eval1("Cancel", expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(b), expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(m), expr_copy(eta) }, 2) }, 2))
            : expr_copy(b);
        Expr* s = rt_is_zero(cm) ? expr_new_integer(0) : rde_tower(coef, cm, &Clo);
        expr_free(coef); expr_free(cm);
        if (!s) { fail = true; break; }
        if (z) {                                     /* + lambda_m z */
            Expr* sz = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ s,
                expr_new_function(expr_new_symbol("Times"), (Expr*[]){ rde_par_fresh(), expr_copy(z) }, 2) }, 2);
            s = rde_canon(sz); expr_free(sz);
        }
        if (rt_is_zero(s)) { expr_free(s); continue; }
        Expr* stm = m == 0 ? expr_copy(s)
            : rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){ expr_copy(s), rc_tau_pow(m, C) }, 2));
        expr_free(s);
        q = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ q, expr_copy(stm) }, 2));
        Expr* Dstm = rc_D(stm, C);
        Expr* bstm = rde_mul(b, stm);
        Expr* sum = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ cc, expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), bstm }, 2),
              expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), Dstm }, 2) }, 3);
        cc = rde_canon(sum);
        expr_free(sum); expr_free(stm);
    }
    expr_free(cc);
    if (z) expr_free(z);
    if (fail) { expr_free(q); return NULL; }
    return q;
}

static Expr* rde_polyrischde_cancel(Expr* b, Expr* c, RdeCtx* C, long n) {
    RtKind kind = C->T->kind[C->m];
    Expr* eta = (kind == RT_EXP) ? C->T->Dcoef[C->m] : NULL;
    if (rt_is_zero(c)) return expr_new_integer(0);
    if (rde_par_in(c)) return rde_polyrischde_cancel_par(b, c, C, n);
    /* Is the homogeneous equation solvable in K, -b = Dw/w, w = z tau^mm?
     * Yes: the log-derivative branch.  Proved no: every coefficient equation
     * of the top-down pass below has at most one solution, so a failed pass is
     * a proof.  Undecided: a failed pass proves nothing. */
    int branch;
    {
        Expr* mb = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(b) }, 2));
        long mm = 0; Expr* z = NULL;
        branch = mb ? rc_param_log_deriv(mb, eta, C, &mm, &z) : -1;
        if (mb) expr_free(mb);
        if (branch == 1) {
            Expr* q = NULL;
            int r = rde_cancel_logderiv(c, z, mm, C, n, &q);
            expr_free(z);
            if (r == 1) return q;
            if (r == 0) return NULL;                 /* proved: no solution */
            rde_cert_mark("PolyRischDECancel: log-derivative branch integral undecided");
            return NULL;
        }
        if (z) expr_free(z);
    }
    Expr* cc = expr_copy(c);
    while (n < rc_deg(cc, C) && rde_par_in(cc)) cc = rde_par_strip_top(cc, C->tau);   /* parametric */
    if (n < rc_deg(cc, C)) {
        expr_free(cc);
        if (branch < 0) rde_cert_mark("PolyRischDECancel: log-derivative test undecided");
        return NULL;
    }
    Expr* q = expr_new_integer(0);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        long m = rc_deg(cc, C);
        if (n < m) {
            if (rde_par_in(cc)) { cc = rde_par_strip_top(cc, C->tau); continue; }
            fail = true; break;
        }
        Expr* lcc = rc_lc(cc, C);                    /* leading coeff, in k */
        Expr* coef;                                  /* b + m*eta (exp) or b (prim) */
        if (kind == RT_EXP)
            coef = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(b), expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(m), expr_copy(eta) }, 2) }, 2));
        else
            coef = expr_copy(b);
        RdeCtx Clo = { C->T, C->m - 1, C->x, (C->m - 1 >= 0) ? C->T->t[C->m - 1] : C->x };
        Expr* s = rde_tower(coef, lcc, &Clo);        /* recursion into the lower field */
        expr_free(coef); expr_free(lcc);
        if (!s) { fail = true; break; }
        Expr* stm;                                   /* s tau^m */
        if (m == 0) stm = expr_copy(s);
        else {
            Expr* tm = rc_tau_pow(m, C);
            stm = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(s), tm }, 2));  /* adopts tm */
        }
        expr_free(s);
        q = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ q, expr_copy(stm) }, 2));      /* q += stm; adopts old q */
        /* cc -= b*stm + D(stm). */
        Expr* Dstm = rc_D(stm, C);
        Expr* bstm = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(b), expr_copy(stm) }, 2));
        {   /* cc -= b stm + D stm, canonicalized (Cancel works term by term on a
             * sum, leaving a cancelled leading coefficient unreduced). */
            Expr* sum = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ cc, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), bstm }, 2),
                  expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), Dstm }, 2) }, 3);  /* adopts cc, bstm, Dstm */
            cc = rde_canon(sum);
            expr_free(sum);
        }
        expr_free(stm);
        n = m - 1;
    }
    expr_free(cc);
    if (fail) {
        /* With the log-derivative branch ruled out (branch == 0) the pass was
         * forced at every step, so its failure is a proof; otherwise not. */
        if (branch < 0) rde_cert_mark("PolyRischDECancel: log-derivative test undecided");
        expr_free(q); return NULL;
    }
    return q;
}

/* Positive integer roots of R(z), a polynomial in z whose coefficients lie in
 * Q(k): an integer root is a common root of every coefficient (as a polynomial
 * in the other symbols), so it is a root of their gcd G(z) over Q[z], and the
 * integer roots of G are read off its linear factors over Q (rational-root
 * exact; no Solve).  Fills roots[0..*n-1] (distinct, >= 1) and returns 1, or 0
 * when undecided (R not a rational function, a factorisation that did not
 * evaluate, a root beyond the long range, more than `cap` roots). */
static int rc_pos_int_roots(Expr* R, Expr* z, long* roots, size_t cap, size_t* n);
/* Collect the distinct numeric radicals Power[number, Rational] and Root objects
 * of e (borrowed pointers). */
static void rc_collect_alg_atoms(Expr* e, Expr*** out, size_t* na, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    bool atom = false;
    if (rt_head_is(e, "Root")) atom = true;
    else if (rt_head_is(e, "Power") && e->data.function.arg_count == 2 && rde_has_alg_const(e)) {
        Expr* b = e->data.function.args[0];
        atom = (b->type == EXPR_INTEGER || b->type == EXPR_BIGINT || rt_head_is(b, "Rational"))
               && rt_head_is(e->data.function.args[1], "Rational");
    }
    if (atom) {
        for (size_t i = 0; i < *na; i++) if (expr_eq((*out)[i], e)) return;
        if (*na == *cap) { *cap = *cap ? 2 * *cap : 4; *out = realloc(*out, *cap * sizeof(Expr*)); }
        (*out)[(*na)++] = e;
        return;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rc_collect_alg_atoms(e->data.function.args[i], out, na, cap);
}
/* Positive integer roots of R over Q(alpha), exactly (campaign step 7): every
 * algebraic atom becomes a symbol s_i with its defining relation (s^q - b^p for
 * b^(p/q), f(s) for Root[f, k]); eliminating them by resultants leaves N in
 * Q[z] whose roots include every root of R (the true values are a common
 * root).  The integer roots of N are the candidates; each is kept only when
 * RootReduce[R(n)] == 0. */
static int rc_pos_int_roots_alg(Expr* R, Expr* z, long* roots, size_t cap, size_t* n) {
    *n = 0;
    Expr* t = rt_eval1("Together", expr_copy(R));
    Expr* num = t ? rt_eval1("Expand", rt_eval1("Numerator", t)) : NULL;
    if (!num) return 0;
    Expr** atoms = NULL; size_t na = 0, acap = 0;
    rc_collect_alg_atoms(num, &atoms, &na, &acap);
    Expr** rels = malloc((na ? na : 1) * sizeof(Expr*));
    Expr** syms = malloc((na ? na : 1) * sizeof(Expr*));
    Expr** rules = malloc((na ? na : 1) * sizeof(Expr*));
    bool ok = true;
    size_t nd = 0;                                    /* entries initialised */
    for (size_t i = 0; i < na && ok; i++, nd++) {
        char nm[48]; snprintf(nm, sizeof nm, "Risch`Private`rcAlg%zu", i);
        syms[i] = expr_new_symbol(nm);
        Expr* a = atoms[i];
        if (rt_head_is(a, "Root")) {
            rels[i] = rt_eval_own(expr_new_function(expr_copy(a->data.function.args[0]),
                (Expr*[]){ expr_copy(syms[i]) }, 1));    /* f[s] */
        } else {
            Expr* ex = a->data.function.args[1];          /* Rational[p, q] */
            Expr* p = ex->data.function.args[0], *q = ex->data.function.args[1];
            rels[i] = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_copy(syms[i]), expr_copy(q) }, 2),
                expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1),
                    expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_copy(a->data.function.args[0]), expr_copy(p) }, 2) }, 2) }, 2));
        }
        if (!rels[i]) ok = false;
        rules[i] = expr_new_function(expr_new_symbol("Rule"), (Expr*[]){ expr_copy(a), expr_copy(syms[i]) }, 2);
    }
    Expr* N = NULL;
    if (ok) {
        Expr* rl = expr_new_function(expr_new_symbol("List"), rules, na);   /* adopts rules[] */
        N = rt_eval2("ReplaceAll", expr_copy(num), rl);
        for (size_t i = 0; N && i < na; i++)
            N = rt_eval3("Resultant", N, expr_copy(rels[i]), expr_copy(syms[i]));
        if (N) N = rt_eval1("Expand", N);
    } else for (size_t i = 0; i < nd; i++) expr_free(rules[i]);
    for (size_t i = 0; i < nd; i++) { if (rels[i]) expr_free(rels[i]); expr_free(syms[i]); }
    free(rels); free(syms); free(rules); free(atoms);
    int res = 0;
    if (N && !rt_is_zero(N) && !rde_has_alg_const(N)) {
        long cand[64]; size_t nc = 0;
        if (rc_pos_int_roots(N, z, cand, 64, &nc)) {
            res = 1;
            for (size_t i = 0; i < nc; i++) {
                Expr* v = rt_eval2("ReplaceAll", expr_copy(num), expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(z), expr_new_integer(cand[i]) }, 2));
                Expr* vr = v ? rt_eval1("RootReduce", v) : NULL;
                bool zero = vr && rt_is_zero(vr);
                if (vr) expr_free(vr);
                if (!vr) { res = 0; break; }
                if (zero) { if (*n >= cap) { res = 0; break; } roots[(*n)++] = cand[i]; }
            }
        }
    }
    if (N) expr_free(N);
    expr_free(num);
    return res;
}
static int rc_pos_int_roots(Expr* R, Expr* z, long* roots, size_t cap, size_t* n) {
    *n = 0;
    if (rde_has_alg_const(R)) return rc_pos_int_roots_alg(R, z, roots, cap, n);
    Expr* tog = rt_eval1("Together", expr_copy(R));
    Expr* num = tog ? rt_eval1("Numerator", tog) : NULL;
    num = num ? rt_eval1("Expand", num) : NULL;
    if (!num) return 0;
    if (rt_is_zero(num)) { expr_free(num); return 0; }   /* R == 0: every z, undecided */
    Expr* vars = rt_eval1("Variables", expr_copy(num));
    Expr* zl = expr_new_function(expr_new_symbol("List"), (Expr*[]){ expr_copy(z) }, 1);
    Expr* rest = vars ? rt_eval2("Complement", vars, zl) : NULL;
    if (!vars) expr_free(zl);
    Expr* cs = NULL;
    if (rest && rest->type == EXPR_FUNCTION && rest->data.function.arg_count == 0) {
        expr_free(rest);
        cs = expr_new_function(expr_new_symbol("List"), (Expr*[]){ num }, 1);
    } else if (rest) {
        Expr* cr = rt_eval2("CoefficientRules", num, rest);
        cs = cr ? rt_eval1("Values", cr) : NULL;
    } else expr_free(num);
    if (!cs) return 0;
    Expr* G = rt_eval2("Apply", expr_new_symbol("PolynomialGCD"), cs);
    if (!G) return 0;
    if (rt_is_zero(G)) { expr_free(G); return 0; }
    Expr* fl = rt_eval1("FactorList", G);
    int ok = fl && rt_head_is(fl, "List");
    for (size_t i = 0; ok && i < fl->data.function.arg_count; i++) {
        Expr* pr = fl->data.function.args[i];               /* {factor, multiplicity} */
        if (!rt_head_is(pr, "List") || pr->data.function.arg_count != 2) { ok = 0; break; }
        Expr* fac = pr->data.function.args[0];
        if (rt_degree(fac, z) != 1) continue;               /* constant or irreducible deg >= 2 */
        Expr* a1 = rt_coeff(fac, z, 1);
        Expr* a0 = rt_coeff(fac, z, 0);
        Expr* root = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), a0, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ a1, expr_new_integer(-1) }, 2) }, 3));   /* adopts a0, a1 */
        if (!root) { ok = 0; break; }
        if (root->type == EXPR_BIGINT) {
            /* a positive bigint residue is a real (if absurd) residue: undecided */
            Expr* pos = rt_eval1("Positive", expr_copy(root));
            bool p = pos && pos->type == EXPR_SYMBOL && pos->data.symbol.name == intern_symbol("True");
            if (pos) expr_free(pos);
            if (p) ok = 0;
        } else if (root->type == EXPR_INTEGER && root->data.integer >= 1) {
            if (*n >= cap) ok = 0; else roots[(*n)++] = (long)root->data.integer;
        }
        expr_free(root);
    }
    if (fl) expr_free(fl);
    return ok;
}

/* ParametricLogarithmicDerivative (Bronstein §7.3), decided exactly through
 * the structure of logarithmic derivatives in the tower.  Question: is
 *     f = m eta + Dz/z,   m in Z,  z in k*,
 * with f, eta in k = K_{C->m - 1} (the field below tau; eta == NULL asks for
 * m = 0, the primitive case)?
 *
 * Why it is a finite linear problem: write z = c Prod p^e_p Prod_{j in E} t_j^s_j
 * over the irreducibles p of Q[x, t_0..t_{m-1}] other than the exponential
 * monomials t_j (which are special: Dt_j/t_j = Deta_j has no pole).  Every such
 * p is normal at its own top variable, so Dp/p has a simple pole along p with
 * residue 1 that no other term of Dz/z can cancel; f = m eta + Dz/z then puts p
 * in den(f) or den(eta).  So z is built from the irreducible factors of
 * den(f) den(eta) and the exponential monomials, and f must lie in the Q-span of
 *     G = {eta} u {Dp/p : p | den(f) den(eta), p not an exp monomial}
 *           u {Deta_j : t_j exponential, j < C->m}.
 * G is Q-linearly independent: a relation would make Prod p^e Prod t_j^s (times
 * a power of the monomial theta = exp(eta) above k) a constant, contradicting
 * unique factorisation in Q[x, t] or the transcendence of theta.  So the span
 * coefficients are unique, and f = m eta + Dz/z holds iff they exist and are
 * all integers.  Returns 1 (and m, z owned), 0 (proved: no such m, z), -1
 * (undecided: a factorisation or the span decision did not complete). */
static int rc_param_log_deriv(Expr* f, Expr* eta, RdeCtx* C, long* m_out, Expr** z_out) {
    *m_out = 0; *z_out = NULL;
    if (C->m < 0) return -1;                         /* k must be a field below tau */
    long lo = C->m - 1;                              /* top index of k */
    if (rt_is_zero(f)) { *z_out = expr_new_integer(1); return 1; }
    size_t nv = (size_t)(lo + 2);
    Expr** vv = malloc(nv * sizeof(Expr*));
    vv[0] = expr_copy(C->x);
    for (long j = 0; j <= lo; j++) vv[j + 1] = expr_copy(C->T->t[j]);
    Expr* vars = expr_new_function(expr_new_symbol("List"), vv, nv);
    free(vv);
    /* candidate irreducibles: factors of den(f) den(eta) over Q in vars */
    Expr* df = rde_den(f);
    Expr* de = eta ? rde_den(eta) : expr_new_integer(1);
    Expr* prod = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ df, de }, 2);
    /* Over k(i) (f or eta carries I) the irreducibles are those of Q(i)[vars]:
     * -2 I/(x^2 + 1) = D log((x + I)/(x - I)) is a logarithmic derivative whose
     * Q-factors x^2 + 1 alone do not span it (review F8). */
    bool gauss = !rt_free_of_head(f, "Complex") || (eta && !rt_free_of_head(eta, "Complex"));
    Expr* fl = gauss
        ? rt_eval2("FactorList", prod, expr_new_function(expr_new_symbol("Rule"),
              (Expr*[]){ expr_new_symbol("Extension"), expr_new_symbol("I") }, 2))
        : rt_eval1("FactorList", prod);
    if (!fl || !rt_head_is(fl, "List")) { if (fl) expr_free(fl); expr_free(vars); return -1; }
    size_t cap = fl->data.function.arg_count + 3 * (size_t)(lo + 2) + 1;
    Expr** gens = malloc(cap * sizeof(Expr*));
    Expr** base = malloc(cap * sizeof(Expr*));       /* the element whose log-derivative gens[i] is */
    int*   kind = malloc(cap * sizeof(int));          /* 0 eta, 1 factor p, 2 exp monomial */
    size_t ng = 0;
    if (eta) { gens[ng] = expr_copy(eta); base[ng] = NULL; kind[ng] = 0; ng++; }
    int rc = 1;
    for (size_t i = 0; i < fl->data.function.arg_count && rc == 1; i++) {
        Expr* pr = fl->data.function.args[i];
        if (!rt_head_is(pr, "List") || pr->data.function.arg_count != 2) { rc = -1; break; }
        Expr* p = pr->data.function.args[0];
        bool dep = false;
        for (size_t j = 0; j < nv && !dep; j++)
            if (!rt_free_of_x(p, vars->data.function.args[j])) dep = true;
        if (!dep) continue;                           /* a constant factor */
        bool expmono = false;
        for (long j = 0; j <= lo && !expmono; j++)
            if (C->T->kind[j] == RT_EXP && expr_eq(p, C->T->t[j])) expmono = true;
        if (expmono) continue;                        /* covered by Deta_j below */
        /* a special factor of a lower tangent monomial (a divisor of t_j^2 + sigma
         * in t_j alone) has a pole-free log-derivative: covered below */
        bool tanspecial = false;
        for (long j = 0; j <= lo && !tanspecial; j++) {
            if (C->T->kind[j] != RT_TAN) continue;
            bool only = true;
            for (size_t k = 0; k < nv && only; k++)
                if (!expr_eq(vars->data.function.args[k], C->T->t[j])
                    && !rt_free_of_x(p, vars->data.function.args[k])) only = false;
            if (!only) continue;
            Expr* sp = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(C->T->t[j]), expr_new_integer(2) }, 2),
                expr_new_integer(C->T->tsg[j]) }, 2);
            Expr* r = rt_eval3("PolynomialRemainder", sp, expr_copy(p), expr_copy(C->T->t[j]));
            if (r && rt_is_zero(r)) tanspecial = true;
            if (r) expr_free(r);
        }
        if (tanspecial) continue;
        RdeCtx Ck = { C->T, lo, C->x, lo >= 0 ? C->T->t[lo] : C->x };
        Expr* Dp = rc_D(p, &Ck);
        if (!Dp) { rc = -1; break; }
        gens[ng] = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ Dp, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(p), expr_new_integer(-1) }, 2) }, 2));
        base[ng] = expr_copy(p); kind[ng] = 1; ng++;
    }
    for (long j = 0; j <= lo && rc == 1; j++) {
        if (C->T->kind[j] != RT_EXP) continue;
        gens[ng] = expr_copy(C->T->Dcoef[j]);          /* Dt_j / t_j */
        base[ng] = expr_copy(C->T->t[j]); kind[ng] = 2; ng++;
    }
    /* Lower tangent monomials, D t_j = eta_j (t_j^2 + sigma): their special
     * factors are pole-free in Dz/z, so they enter z without showing in any
     * denominator.  Over Q (sigma = +1) the factor t_j^2 + 1, D log = 2 eta_j t_j;
     * over Q(i), or for sigma = -1, the linear t_j - r (r^2 = -sigma),
     * D log = eta_j (t_j + r). */
    for (long j = 0; j <= lo && rc == 1; j++) {
        if (C->T->kind[j] != RT_TAN) continue;
        Expr* tj = C->T->t[j];
        Expr* ej = C->T->Dcoef[j];
        if (C->T->tsg[j] == 1 && !gauss) {
            gens[ng] = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(2), expr_copy(ej), expr_copy(tj) }, 3));
            base[ng] = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(tj), expr_new_integer(2) }, 2),
                expr_new_integer(1) }, 2);
            kind[ng] = 1; ng++;
            continue;
        }
        for (int sgn = -1; sgn <= 1; sgn += 2) {
            Expr* r = C->T->tsg[j] == 1
                ? expr_new_function(expr_new_symbol("Complex"),
                      (Expr*[]){ expr_new_integer(0), expr_new_integer(sgn) }, 2)
                : expr_new_integer(sgn);
            gens[ng] = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(ej), expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(tj), expr_copy(r) }, 2) }, 2));
            base[ng] = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(tj), expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), r }, 2) }, 2));      /* adopts r */
            kind[ng] = 1; ng++;
        }
    }
    if (rc == 1) {
        bool failed = false;
        Expr* co = risch_rational_span(f, gens, ng, vars, &failed);
        if (!co) rc = failed ? -1 : 0;
        else {
            /* integrality, and z = Prod base_i^coeff_i */
            Expr** facs = malloc((ng + 1) * sizeof(Expr*));
            size_t nf = 0;
            for (size_t i = 0; i < ng && rc == 1; i++) {
                Expr* ci = co->data.function.args[i];
                if (ci->type != EXPR_INTEGER) {
                    /* rational non-integer, or a bigint: no z in k*, or absurd */
                    rc = (ci->type == EXPR_BIGINT) ? -1 : 0;
                    break;
                }
                if (kind[i] == 0) { *m_out = (long)ci->data.integer; continue; }
                if (ci->data.integer == 0) continue;
                facs[nf++] = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(base[i]), expr_copy(ci) }, 2);
            }
            if (rc == 1)
                *z_out = nf == 0 ? expr_new_integer(1)
                       : rt_eval1("Together", expr_new_function(expr_new_symbol("Times"), facs, nf));
            else for (size_t i = 0; i < nf; i++) expr_free(facs[i]);
            free(facs);
            expr_free(co);
        }
    }
    for (size_t i = 0; i < ng; i++) { expr_free(gens[i]); if (base[i]) expr_free(base[i]); }
    free(gens); free(base); free(kind);
    expr_free(fl); expr_free(vars);
    return rc;
}

/* WeakNormalizer (Bronstein §6.1, p.183) over k(x)(tau), m >= 0, or over the
 * base field C(x), m < 0 (where every squarefree factor is normal): returns q in
 * k[tau] such that f - Dq/q is weakly normalized -- no normal simple pole with a
 * positive-integer residue -- so that RdeNormalDenominator's bound is exact
 * (Thm 6.1.2) and a "no solution" below is a proof.  With f = A/(E d1), d1 the
 * product of the normal simple poles, the residues are the roots of
 *     R(z) = Res_tau(d1, A - z E Dd1)
 * and q = Prod_n gcd(A - n E Dd1, d1)^n over the positive integer roots n.
 * Returns q (owned; 1 when f is already weakly normalized), or NULL when the
 * residues could not be decided. */
static Expr* rc_weak_normalizer(Expr* f, RdeCtx* C) {
    Expr* den = rc_den(f, C);
    Expr* dn = rc_normal_part(den, C);
    Expr* q = expr_new_integer(1);
    if (rc_deg(dn, C) >= 1) {
        Expr* ddn = rc_D(dn, C);
        Expr* g = rc_gcd(dn, ddn, C);
        Expr* dstar = rc_quot(dn, g, C);
        Expr* gds = rc_gcd(dstar, g, C);
        Expr* d1 = rc_quot(dstar, gds, C);
        expr_free(ddn); expr_free(g); expr_free(dstar); expr_free(gds);
        if (rc_deg(d1, C) >= 1) {
            Expr* A = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(f), expr_copy(den) }, 2));
            Expr* E = rc_quot(den, d1, C);
            Expr* Dd1 = rc_D(d1, C);
            Expr* z = expr_new_symbol("Risch`Private`rdeWNz");
            Expr* P = rt_eval1("Together", expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(A), expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(z), expr_copy(E), expr_copy(Dd1) }, 4) }, 2));
            Expr* Pn = P ? rt_eval1("Numerator", P) : NULL;
            Expr* d1t = rt_eval1("Together", expr_copy(d1));
            Expr* d1n = d1t ? rt_eval1("Numerator", d1t) : NULL;
            Expr* R = (Pn && d1n) ? rt_eval3("Resultant", d1n, Pn, expr_copy(C->tau)) : NULL;
            long roots[64]; size_t nr = 0;
            int ok = R && rc_pos_int_roots(R, z, roots, 64, &nr);
            if (R) expr_free(R);
            if (!ok) { expr_free(q); q = NULL; }
            for (size_t i = 0; ok && i < nr; i++) {
                Expr* An = rt_eval1("Together", expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(A), expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-roots[i]), expr_copy(E), expr_copy(Dd1) }, 3) }, 2));
                Expr* Ann = An ? rt_eval1("Numerator", An) : NULL;
                Expr* gn = Ann ? rc_gcd(Ann, d1, C) : NULL;
                if (Ann) expr_free(Ann);
                if (!gn) { expr_free(q); q = NULL; break; }
                if (rc_deg(gn, C) >= 1) {
                    Expr* gp = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ gn, expr_new_integer(roots[i]) }, 2);   /* adopts gn */
                    Expr* nq = rde_mul(q, gp);
                    expr_free(q); expr_free(gp); q = nq;
                } else expr_free(gn);
            }
            expr_free(z); expr_free(A); expr_free(E); expr_free(Dd1);
        }
        expr_free(d1);
    }
    expr_free(dn); expr_free(den);
    if (q && !(q->type == EXPR_INTEGER)) { Expr* c = rde_canon(q); expr_free(q); q = c; }
    return q;
}

/* The base-field entry (rde_base): the same normalizer with D = d/dx. */
static Expr* rde_weak_normalizer(Expr* f, Expr* x) {
    RdeCtx Cb = { NULL, -1, x, x };
    Expr* q = rc_weak_normalizer(f, &Cb);
    if (!q) { rde_cert_mark("WeakNormalizer: residues not decided"); return expr_new_integer(1); }
    return q;
}

/* Recursive Risch DE driver: solve D_tower[y] + f y = g for y in K_m (Bronstein
 * §5.2 dispatch + Chapter 6).  m < 0 is the base field C(x) (delegates to the
 * complete, audited rde_base).  For m >= 0 with a PRIMITIVE (RT_LOG) top this
 * runs RdeNormalDenominator -> RdeBoundDegree -> SPDE -> PolyRischDENoCancel (1a);
 * an EXPONENTIAL top, or a branch needing the b=0 / cancellation solvers, returns
 * NULL so the caller falls back to the bounded ansatz until later increments land.
 * Returns y (owned), or NULL.  Correct by construction: the returned y is checked
 * by the exact tower identity D_tower[y] + f y == g, so a reduction slip declines. */
/* A tangent top whose special polynomial splits -- Tanh (t^2 - 1), or Tan over
 * a field containing I (t^2 + 1 = (t - I)(t + I)) -- is an exponential in
 * disguise (Bronstein, Symbolic Integration I, p. 166).  With alpha = 1 (Tanh)
 * or I (Tan), theta = (alpha + t)/(alpha - t) has
 *     D theta / theta = 2 alpha Dt / (alpha^2 - t^2) = a,  a in k,
 * so theta is an exponential monomial over k (E^(2 v) for Tanh[v], E^(-2 I v)
 * for Tan[v]) and k(t) = k(theta), t = alpha (theta - 1)/(theta + 1).  The
 * Risch DE D y + f y = g over k(t) is the one over k(theta) after that
 * substitution: an isomorphism of differential fields over k, so the
 * exponential stack's answer maps back exactly, and its authoritative "no
 * solution" is one here.  The tower is copied with level m replaced (levels
 * 0..m only; the lower ones are shared, borrowed).  f, g borrowed; y owned. */
static Expr* rde_tower_mobius(Expr* f, Expr* g, RdeCtx* C, Expr* alpha) {
    RtTower* T = C->T;
    long m = C->m;
    Expr* t = T->t[m];
    Expr* a2 = expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_copy(alpha), expr_new_integer(2) }, 2);
    Expr* t2 = expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_copy(t), expr_new_integer(2) }, 2);
    Expr* den = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ a2,
        expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), t2 }, 2) }, 2);
    Expr* a = rt_eval1("Together", rt_eval_own(expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        expr_new_integer(2), expr_copy(alpha), rt_dt_i(T, (size_t)m),
        expr_new_function(expr_new_symbol("Power"), (Expr*[]){ den, expr_new_integer(-1) }, 2) }, 4)));
    if (!a || !rt_free_of_x(a, t)) {
        if (a) expr_free(a);
        rde_cert_mark("tower: split tangent top, Dtheta/theta not in k");
        return NULL;
    }
    /* the copy: levels 0..m, level m the exponential theta */
    size_t n = (size_t)m + 1;
    RtTower T2 = *T;
    T2.n = n;
    T2.kind = malloc(n * sizeof(RtKind));
    T2.t = malloc(n * sizeof(Expr*));
    T2.Dcoef = malloc(n * sizeof(Expr*));
    T2.tsg = malloc(n * sizeof(long));
    T2.kernel = malloc(n * sizeof(Expr*));
    T2.arg = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        T2.kind[i] = T->kind[i]; T2.t[i] = T->t[i]; T2.Dcoef[i] = T->Dcoef[i];
        T2.tsg[i] = T->tsg[i]; T2.kernel[i] = T->kernel[i]; T2.arg[i] = T->arg[i];
    }
    char nm[48];
    snprintf(nm, sizeof nm, "Risch`Private`rmMob%ld", m);
    Expr* th = expr_new_symbol(nm);
    /* theta's exponent, 2 v / alpha: 2 v for Tanh[v], -2 I v for Tan[v]
     * (D of it is a; recorded for completeness, the RDE reads only Dcoef) */
    Expr* argm = rt_eval_own(expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        expr_new_integer(2), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(alpha), expr_new_integer(-1) }, 2), expr_copy(T->arg[m]) }, 3));
    T2.kind[m] = RT_EXP; T2.tsg[m] = 0; T2.t[m] = th; T2.Dcoef[m] = a; T2.arg[m] = argm;
    T2.kernel[m] = expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_new_symbol("E"), expr_copy(argm) }, 2);
    /* t = alpha (theta - 1)/(theta + 1) */
    Expr* tin = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_copy(alpha),
        expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(th), expr_new_integer(-1) }, 2),
        expr_new_function(expr_new_symbol("Power"), (Expr*[]){
            expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(th), expr_new_integer(1) }, 2),
            expr_new_integer(-1) }, 2) }, 3);
    Expr* r1 = expr_new_function(expr_new_symbol("Rule"), (Expr*[]){ expr_copy(t), tin }, 2);
    Expr* f2 = rt_eval1("Together", rt_eval2("ReplaceAll", expr_copy(f), expr_copy(r1)));
    Expr* g2 = rt_eval1("Together", rt_eval2("ReplaceAll", expr_copy(g), r1));
    Expr* y = NULL;
    if (f2 && g2) {
        RdeCtx C2 = { &T2, m, C->x, th };
        Expr* y2 = rde_tower(f2, g2, &C2);
        if (y2) {
            /* theta = (alpha + t)/(alpha - t) */
            Expr* tout = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(alpha), expr_copy(t) }, 2),
                expr_new_function(expr_new_symbol("Power"), (Expr*[]){
                    expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ expr_copy(alpha),
                        expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(t) }, 2) }, 2),
                    expr_new_integer(-1) }, 2) }, 2);
            y = rt_eval1("Together", rt_eval2("ReplaceAll", y2,
                expr_new_function(expr_new_symbol("Rule"), (Expr*[]){ expr_copy(th), tout }, 2)));
        }
    }
    if (f2) expr_free(f2);
    if (g2) expr_free(g2);
    expr_free(T2.kernel[m]); expr_free(argm); expr_free(a); expr_free(th);
    free(T2.kind); free(T2.t); free(T2.Dcoef); free(T2.tsg); free(T2.kernel); free(T2.arg);
    return y;
}

Expr* rde_tower(Expr* f, Expr* g, RdeCtx* C) {
    if (rt_is_zero(g)) return expr_new_integer(0);
    if (C->m < 0) return rde_base(f, g, C->x);
    RtKind kind = C->T->kind[C->m];
    /* A hypertangent top (D tau = eta (tau^2 + 1)) runs the same stack with the
     * delta = 2 boxes (review F8).  The hyperbolic sigma = -1 top, whose special
     * part tau^2 - 1 splits, and a Liouvillian primitive are not yet covered. */
    bool tan = rc_is_tan_top(C);
    /* a Liouvillian primitive top (Erf, Ei, li; D tau = eta in k) runs the
     * primitive stack: no special part, RdeBoundDegreePrim, PolyRischDECancel
     * with eta absent -- every box reads Dcoef, not the kernel (campaign step
     * 3b).  Its LimitedIntegrate goes through the field integrator's RT_PRIM
     * branch. */
    /* a split special part (Tanh, or Tan over k(i)): the Moebius change to the
     * exponential theta (rde_tower_mobius) */
    if (kind == RT_TAN && C->T->tsg[C->m] == -1) {
        Expr* one = expr_new_integer(1);
        Expr* y = rde_tower_mobius(f, g, C, one);
        expr_free(one);
        return y;
    }
    if (kind != RT_LOG && kind != RT_EXP && kind != RT_PRIM && !tan) {
        /* unreachable: every RtKind is handled (a Tanh top just above) */
        rde_cert_mark("tower: a top of no known kind");
        return NULL;
    }
    /* Over algebraic constants the tower stack's arithmetic (Together, Cancel,
     * PolynomialRemainder over k(x)[tau]) is not canonical: Together returns an
     * uncombined sum when an algebraic constant sits in a denominator, and a zero
     * coefficient can stay unreduced.  Only the base field has the exact mode
     * (rde_base's algebraic mode, step 7); here a solution is still found and
     * certified by the identity gate, but a "no" proves nothing.  Hidden until
     * v0.6606 by the base RDE's own marks: D[(1 + Sqrt[3]) E^((1 + Sqrt[2]) u)
     * u^(2 Sqrt[2])/(1 + Sqrt[2] u)] was then certified non-elementary. */
    if (rde_has_alg_const(f) || rde_has_alg_const(g))
        rde_cert_mark("tower: algebraic constants (no exact tower arithmetic)");
    /* The hypertangent boxes assume sqrt(-1) is not in k, so that t^2 + 1 is
     * irreducible (Bronstein 5.10, 6.2).  Over k(i) it splits into t - i and
     * t + i, two special factors with independent orders, and a one-sided pole
     * (1/(Tan[x] - I)) is invisible to the combined t^2 + 1 valuation.  Solve
     * anyway (the identity gate certifies any solution) but prove nothing. */
    if (tan && (!rt_free_of_head(f, "Complex") || !rt_free_of_head(g, "Complex"))) {
        Expr* ii = expr_new_function(expr_new_symbol("Complex"),
            (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2);
        Expr* y = rde_tower_mobius(f, g, C, ii);
        expr_free(ii);
        return y;
    }
    /* the parametric conversion covers the log and exponential stacks; the
     * tangent boxes (NoCancel2/3, CancelTan's re/im split) are not converted */
    if (tan && rde_par_in(g)) { rde_cert_mark("parametric: hypertangent top"); return NULL; }
    /* Weak normalization (Bronstein §6.1): RdeNormalDenominator's bound, and so
     * every "no solution" below, assumes f weakly normalized in tau.  With
     * q = WeakNormalizer(f), y solves D y + f y = g iff z = q y solves
     * D z + (f - Dq/q) z = q g. */
    {
        Expr* wq = rc_weak_normalizer(f, C);
        if (!wq) rde_cert_mark("tower: WeakNormalizer residues not decided");
        else if (!(wq->type == EXPR_INTEGER && wq->data.integer == 1)) {
            Expr* Dq = rc_D(wq, C);
            Expr* fexpr = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(f), expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), Dq, expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(wq), expr_new_integer(-1) }, 2) }, 3) }, 2);  /* adopts Dq */
            Expr* fbar = rde_canon(fexpr);
            expr_free(fexpr);
            Expr* gbar = rde_mul(wq, g);
            Expr* z = rde_tower(fbar, gbar, C);
            expr_free(fbar); expr_free(gbar);
            Expr* y = z ? rde_divq(z, wq) : NULL;
            if (z) expr_free(z);
            expr_free(wq);
            return y;
        } else expr_free(wq);
    }

    Expr *a, *b, *c, *h;
    if (!rde_normal_denominator_field(f, g, C, &a, &b, &c, &h)) return NULL;
    /* Special part of the denominator: primitive (RT_LOG) tops have none
     * (h_s = 1); an exponential top clears the tau-power poles (§6.2). */
    Expr* h_s = expr_new_integer(1);
    if (kind == RT_EXP || tan) {
        Expr *na, *nb, *nc, *nhs;
        if (tan) rde_special_denom_tan(a, b, c, C, &na, &nb, &nc, &nhs);
        else rde_special_denom_exp(a, b, c, C, &na, &nb, &nc, &nhs);
        expr_free(a); expr_free(b); expr_free(c); expr_free(h_s);
        a = na; b = nb; c = nc; h_s = nhs;
    }
    Expr* result = NULL;
    if (rc_ispoly(a, C) && rc_ispoly(b, C) && !rc_ispoly(c, C) && rde_par_in(c)) {
        Expr* pp = rde_par_polypart(c, C->tau);      /* parametric: the rest must vanish */
        if (pp) { expr_free(c); c = pp; }
    }
    if (!(rc_ispoly(a, C) && rc_ispoly(b, C) && rc_ispoly(c, C)))
        rde_cert_mark("tower: reduced equation not polynomial");
    if (rc_ispoly(a, C) && rc_ispoly(b, C) && rc_ispoly(c, C)) {
        long n = tan ? rde_bound_degree_tan(a, b, c, C)
               : (kind == RT_EXP) ? rde_bound_degree_exp(a, b, c, C)
                                  : rde_bound_degree_prim(a, b, c, C);
        RdeSpde sp;
        if (rde_spde_field(a, b, c, C, n, &sp)) {
            Expr* H = NULL;
            if (rt_is_zero(sp.b)) {
                /* D H = sp.c (Gap 2 antidifferentiation).  sp.c == 0 -> H = 0
                 * (q = beta).  Otherwise integrate sp.c within the field K_m via
                 * the integrator (Bronstein IntegratePrimitivePolynomial §5.8 /
                 * IntegrateHyperexponentialPolynomial §5.9 — the mutual recursion):
                 * a rational-in-tower antiderivative IS the polynomial solution h; a
                 * result carrying a NEW logarithm means the Risch DE has no rational
                 * solution -> decline. */
                if (rt_is_zero(sp.c)) H = expr_new_integer(0);
                else if (rde_par_in(sp.c)) {
                    /* parametric: the linear cancellation pass (b = 0) */
                    Expr* zero = expr_new_integer(0);
                    H = rde_polyrischde_cancel(zero, sp.c, C, sp.m);
                    expr_free(zero);
                } else {
                    /* D H = c, deg H <= m: the log-derivative branch with w = 1
                     * (three-way integral in K, the constant of integration
                     * fixed by the degree; campaign step 1) */
                    Expr* one = expr_new_integer(1);
                    int r = rde_cancel_logderiv(sp.c, one, 0, C, sp.m, &H);
                    expr_free(one);
                    if (r < 0) rde_cert_mark("tower: antidifferentiation D H = c");
                }
            } else if (tan) {
                /* delta = 2 (Bronstein §6.5): deg b >= 2 never cancels; deg b = 0
                 * is NoCancel2, deg b = 1 NoCancel3 (with PolyRischDECancelTan at
                 * the resonant degree). */
                long dbt = rc_deg(sp.b, C);
                H = dbt >= 2 ? rde_polyrischde_nocancel_field(sp.b, sp.c, C, sp.m)
                  : dbt == 0 ? rde_nocancel2_tan(sp.b, sp.c, C, sp.m)
                             : rde_nocancel3_tan(sp.b, sp.c, C, sp.m);
            } else if (rc_deg(sp.b, C) >= 1) {
                /* Non-cancellation (deg_tau(b) >= 1): the leading terms of D[q] and
                 * b q never cancel, for a primitive OR exponential top (delta<=1). */
                H = rde_polyrischde_nocancel_field(sp.b, sp.c, C, sp.m);
            } else {
                /* deg_tau(sp.b) == 0, b in k*: the CANCELLATION case
                 * (PolyRischDECancel{Prim,Exp}, §6.6) — recurse into the lower field. */
                H = rde_polyrischde_cancel(sp.b, sp.c, C, sp.m);
            }
            if (H) {
                Expr* aH = rde_mul(sp.alpha, H);
                Expr* q = rde_add(aH, sp.beta);      /* q_bar = alpha H + beta */
                expr_free(aH); expr_free(H);
                /* y = q_bar / (h_s h) */
                Expr* hsh = rde_mul(h_s, h);
                result = rde_divq(q, hsh);
                expr_free(hsh); expr_free(q);
            }
            rde_spde_free(&sp);
        }
    }
    expr_free(a); expr_free(b); expr_free(c); expr_free(h); expr_free(h_s);

    /* Correct-by-construction gate: D_tower[y] + f y - g must be exactly 0. */
    if (result) {
        Expr* Dy = rc_D(result, C);
        Expr* fy = rde_mul(f, result);
        Expr* mg = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(g) }, 2));
        Expr* resid = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ Dy, fy, mg }, 3);
        Expr* tog = rt_eval1("Together", resid);
        Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;
        /* Expand the numerator: Together clears denominators but does NOT expand
         * products, so a genuinely-zero residual with an unexpanded factor (e.g.
         * D_tower of a polynomial-coefficient exp term, t*(x-1)) is only reduced to
         * 0 after Expand. */
        Expr* rexp = rnum ? rt_eval1("Expand", rnum) : NULL;
        bool good = rexp && rt_is_zero(rexp);
        /* parametric: the residual, linear in the c_i, is the constraint (a
         * residual the unknowns cancelled out of is one too: dropped parts that
         * no c removes, so the system is then inconsistent) */
        if (!good && rexp && (rde_par_in(rexp) || rde_par_in(g))) { rde_par_require(expr_copy(rexp)); good = true; }
        if (rexp) expr_free(rexp);
        if (!good) { rde_cert_mark("tower: identity gate"); expr_free(result); result = NULL; }
    }
    return result;
}

/* ---- The parametric Risch DE (campaign step 3a) ---------------------------
 * Does D y + f y - g vanish exactly in K_m (the tower identity, outside
 * parametric mode)?  y, f, g borrowed. */
static bool rde_identity_holds(Expr* y, Expr* f, Expr* g, RdeCtx* C) {
    Expr* Dy = rc_D(y, C);
    if (!Dy) return false;
    Expr* fy = rde_mul(f, y);
    Expr* resid = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ Dy, fy,
        expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(g) }, 2) }, 3);
    Expr* tog = rt_eval1("Together", resid);
    Expr* num = tog ? rt_eval1("Numerator", tog) : NULL;
    Expr* ex = num ? rt_eval1("Expand", num) : NULL;
    bool ok = ex && rt_is_zero(ex);
    if (ex) expr_free(ex);
    return ok;
}

/* e /. rules /. (every unknown of `unk` -> 0); e borrowed, result owned. */
static Expr* rde_par_subst(Expr* e, Expr* rules, Expr* unk) {
    Expr* r = rules ? rt_eval2("ReplaceAll", expr_copy(e), expr_copy(rules)) : expr_copy(e);
    size_t nu = unk->data.function.arg_count;
    Expr** zr = malloc((nu ? nu : 1) * sizeof(Expr*));
    for (size_t j = 0; j < nu; j++)
        zr[j] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(unk->data.function.args[j]), expr_new_integer(0) }, 2);
    Expr* zl = expr_new_function(expr_new_symbol("List"), zr, nu);
    free(zr);
    Expr* out = rt_eval1("Together", rt_eval2("ReplaceAll", r, zl));
    return out;
}

int rde_tower_param(Expr* f, Expr* g, RdeCtx* C, Expr* unknowns,
                    Expr** y_out, Expr** sol_out) {
    *y_out = NULL; *sol_out = NULL;
    /* a stale state left by a TimeConstrained unwind (the counter was reset,
     * the arrays were not) */
    if (!g_rde_par && (g_rde_par_cons || g_rde_par_unk)) {
        for (size_t i = 0; i < g_rde_par_ncons; i++) expr_free(g_rde_par_cons[i]);
        free(g_rde_par_cons);
        if (g_rde_par_unk) expr_free(g_rde_par_unk);
        g_rde_par_cons = NULL; g_rde_par_unk = NULL;
        g_rde_par_ncons = g_rde_par_cap = 0;
    }
    /* nesting: save the outer solve's state */
    int    s_par = g_rde_par;
    Expr*  s_unk = g_rde_par_unk;
    Expr** s_cons = g_rde_par_cons;
    size_t s_n = g_rde_par_ncons, s_cap = g_rde_par_cap;
    g_rde_par = 1;
    g_rde_par_unk = expr_copy(unknowns);
    g_rde_par_cons = NULL; g_rde_par_ncons = g_rde_par_cap = 0;

    RdeCertState saved = rde_cert_begin();
    Expr* y = rde_tower(f, g, C);
    bool unc = rde_cert_uncertain();
    const char* why = rde_cert_reason();
    rde_cert_end(saved);

    Expr*  unk = g_rde_par_unk;                      /* the caller's + fresh unknowns */
    Expr** cons = g_rde_par_cons;
    size_t ncons = g_rde_par_ncons;
    g_rde_par = s_par; g_rde_par_unk = s_unk;
    g_rde_par_cons = s_cons; g_rde_par_ncons = s_n; g_rde_par_cap = s_cap;

    int status = -1;
    Expr* rules = NULL;                              /* List of Rule, or NULL: none */
    if (y) {
        if (ncons == 0) {
            rules = expr_new_function(expr_new_symbol("List"), NULL, 0);
        } else {
            /* the constraints are polynomials in x and the tower variables,
             * linear in the unknowns; one equation sum_k cons_k w^k = 0 in an
             * auxiliary variable w keeps them apart in CoefficientList */
            Expr* w = expr_new_symbol("rmPw");
            Expr** tm = malloc(ncons * sizeof(Expr*));
            for (size_t k = 0; k < ncons; k++)
                tm[k] = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_copy(cons[k]),
                    expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(w), expr_new_integer((int64_t)k) }, 2) }, 2);
            Expr* lhs = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"), tm, ncons));
            free(tm);
            size_t nv = (C->T ? C->T->n : 0) + 2;
            Expr** vv = malloc(nv * sizeof(Expr*));
            vv[0] = expr_copy(C->x);
            for (size_t i = 0; C->T && i < C->T->n; i++) vv[i + 1] = expr_copy(C->T->t[i]);
            vv[nv - 1] = w;
            Expr* vars = expr_new_function(expr_new_symbol("List"), vv, nv);
            free(vv);
            Expr* sol = rt_solve_ansatz(expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ lhs, expr_new_integer(0) }, 2), vars, unk);
            if (sol && rt_head_is(sol, "List")) {
                if (sol->data.function.arg_count == 0) status = 0;       /* inconsistent */
                else {
                    rules = expr_copy(sol->data.function.args[0]);
                }
            }
            if (sol) expr_free(sol);
        }
        if (rules) {
            Expr* y0 = rde_par_subst(y, rules, unk);
            Expr* g0 = rde_par_subst(g, rules, unk);
            if (y0 && g0 && rde_identity_holds(y0, f, g0, C)) {
                /* the caller's unknowns, solved (a free one is 0) */
                size_t m = unknowns->data.function.arg_count;
                Expr** cv = malloc((m ? m : 1) * sizeof(Expr*));
                for (size_t j = 0; j < m; j++)
                    cv[j] = rde_par_subst(unknowns->data.function.args[j], rules, unk);
                *sol_out = expr_new_function(expr_new_symbol("List"), cv, m);
                free(cv);
                *y_out = y0; y0 = NULL;
                status = 1;
            } else if (!why) why = "parametric: solution failed the identity gate";
            if (y0) expr_free(y0);
            if (g0) expr_free(g0);
            expr_free(rules);
        } else if (status < 0 && !why) why = "parametric: constraint system undecided";
        expr_free(y);
    } else status = 0;                               /* a c-free exit: no c works */
    for (size_t i = 0; i < ncons; i++) expr_free(cons[i]);
    free(cons);
    expr_free(unk);
    if (status == 0 && unc) status = -1;             /* a "no" behind a mark proves nothing */
    if (status < 0) rde_cert_mark(why ? why : "parametric: undecided");
    return status;
}
