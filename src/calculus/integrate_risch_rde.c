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

#include <stdbool.h>
#include <stddef.h>
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

/* Numerator / Denominator of Together[e] (e not consumed). */
static Expr* rde_num(Expr* e) {
    return rt_eval1("Numerator", rt_eval1("Together", expr_copy(e)));
}
static Expr* rde_den(Expr* e) {
    return rt_eval1("Denominator", rt_eval1("Together", expr_copy(e)));
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
    return rt_eval2("PolynomialGCD", expr_copy(a), expr_copy(b));
}
static Expr* rde_quot(Expr* a, Expr* b, Expr* x) {
    return rt_eval3("PolynomialQuotient", expr_copy(a), expr_copy(b), expr_copy(x));
}
static Expr* rde_rem(Expr* a, Expr* b, Expr* x) {
    return rt_eval3("PolynomialRemainder", expr_copy(a), expr_copy(b), expr_copy(x));
}
/* Expanded / cancelled arithmetic (args not consumed). */
static Expr* rde_mul(Expr* a, Expr* b) {
    return rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(a), expr_copy(b) }, 2));
}
static Expr* rde_add(Expr* a, Expr* b) {
    return rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(a), expr_copy(b) }, 2));
}
static Expr* rde_sub(Expr* a, Expr* b) {
    return rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(b) }, 2) }, 2));
}
/* Field quotient a/b, FULLY reduced (args not consumed).  Uses Together, not
 * Cancel: `Cancel` does not reduce a ratio whose operands are themselves nested
 * rationals — e.g. Cancel[(1+x)/(1+1/x)] stays put, while Together gives x.  The
 * coefficient field k(x) of a tower monomial routinely produces such forms
 * (lc(b) = (1+x)/x printed as 1+1/x), and an unreduced quotient makes leading
 * terms fail to cancel in the SPDE / PolyRischDE loops, losing genuine solutions
 * (RISCH_AUDIT_A4.md — the nested-exponential Laurent case). */
static Expr* rde_divq(Expr* a, Expr* b) {
    return rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(b), expr_new_integer(-1) }, 2) }, 2));
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
    return rt_eval1("Expand", rt_eval1("Cancel", rt_eval1("Together", expr_copy(e))));
}
/* Leading coefficient of poly p in x (0 if p == 0). */
static Expr* rde_lc(Expr* p, Expr* x) {
    long d = rt_degree(p, x);
    return d < 0 ? expr_new_integer(0) : rt_coeff(p, x, d);
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
void rde_cert_mark(const char* why) {
    if (!g_rde_cert.active) return;
    if (!g_rde_cert.uncertain) g_rde_cert.why = why;
    g_rde_cert.uncertain = 1;
}

static int rde_spde(Expr* a, Expr* b, Expr* c, Expr* x, long n, RdeSpde* out) {
    out->b = out->c = out->alpha = out->beta = NULL; out->m = 0;
    if (n < 0) {
        if (rt_is_zero(c)) {
            out->b = expr_new_integer(0); out->c = expr_new_integer(0);
            out->m = 0; out->alpha = expr_new_integer(0); out->beta = expr_new_integer(0);
            return 1;
        }
        return 0;                                   /* c != 0, n < 0: no solution */
    }
    Expr* g = rde_gcd(a, b, x);
    Expr* rem = rde_rem(c, g, x);
    bool gdivc = rt_is_zero(rem);
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
    Expr* eg = rt_eval3("PolynomialExtendedGCD", expr_copy(b1), expr_copy(a1), expr_copy(x));
    Expr* s = NULL;
    if (eg && eg->type == EXPR_FUNCTION && eg->data.function.arg_count == 2) {
        Expr* cof = eg->data.function.args[1];
        if (cof && cof->type == EXPR_FUNCTION && cof->data.function.arg_count == 2)
            s = cof->data.function.args[0];
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
static Expr* rde_polyrischde_nocancel1(Expr* b, Expr* c, Expr* x, long n) {
    long db = rt_degree(b, x);
    Expr* lcb = rde_lc(b, x);
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        long dc = rt_degree(cc, x);
        long m = dc - db;
        if (n < 0 || m < 0 || m > n) { fail = true; break; }
        Expr* lcc = rde_lc(cc, x);
        Expr* coeff = rde_divq(lcc, lcb);           /* lc(c)/lc(b) */
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
        if (n < 0 || m < 0 || m > n) { fail = true; break; }
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
    if (x->type == EXPR_SYMBOL) {
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
        expr_free(dw);
    }
    Expr* gbar = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(w), expr_copy(g) }, 2));

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
        bool good = exn && rt_is_zero(exn);
        expr_free(dy); expr_free(fy); expr_free(s1);
        if (exn) expr_free(exn);
        if (!good) { rde_cert_mark("base: identity gate"); expr_free(result); result = NULL; }
    }
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
    bool ok = rt_is_zero(rem);
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

/* Integrate F (an element of K = k(tau), C's field) inside K: 1 with *Q owned
 * when the integrator returns an antiderivative rational in x, t_0..t_m (no new
 * logarithm, no remainder); -1 otherwise (not a proof: a new logarithm may
 * still combine into K, and the integrator may decline).  The call runs in its
 * own certification scope; an uncertain step inside it is propagated. */
static int rc_integrate_in_K(Expr* F, RdeCtx* C, Expr** Q) {
    *Q = NULL;
    if (rt_is_zero(F)) { *Q = expr_new_integer(0); return 1; }
    RdeCertState saved = rde_cert_begin();
    Expr* rem = NULL;
    Expr* R = rt_field_integrate(F, C->T, C->m, C->x, &rem);
    bool unc = rde_cert_uncertain();
    const char* why = rde_cert_reason();
    rde_cert_end(saved);
    if (unc) rde_cert_mark(why);
    bool ok = R && (!rem || rt_is_zero(rem)) && rc_is_tower_rational(R, C);
    if (rem) expr_free(rem);
    if (!ok) { if (R) expr_free(R); return -1; }
    *Q = R;
    return 1;
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
    if (rc_integrate_in_K(alpha, C, &Q) != 1) return -1;
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
        } else rde_cert_mark("RdeBoundDegreePrim: limited integration undecided");
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
            } else rde_cert_mark("RdeBoundDegreePrim: limited integration undecided");
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
        if (rt_is_zero(c)) {
            out->b = expr_new_integer(0); out->c = expr_new_integer(0);
            out->m = 0; out->alpha = expr_new_integer(0); out->beta = expr_new_integer(0);
            return 1;
        }
        return 0;
    }
    Expr* g = rc_gcd(a, b, C);
    Expr* rem = rc_rem(c, g, C);
    bool gdivc = rt_is_zero(rem);
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
        if (n < 0 || m < 0 || m > n) { fail = true; break; }
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
    if (ok != 1) { expr_free(w); return -1; }
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

static Expr* rde_polyrischde_cancel(Expr* b, Expr* c, RdeCtx* C, long n) {
    RtKind kind = C->T->kind[C->m];
    Expr* eta = (kind == RT_EXP) ? C->T->Dcoef[C->m] : NULL;
    if (rt_is_zero(c)) return expr_new_integer(0);
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
    if (n < rc_deg(c, C)) {
        if (branch < 0) rde_cert_mark("PolyRischDECancel: log-derivative test undecided");
        return NULL;
    }
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        long m = rc_deg(cc, C);
        if (n < m) { fail = true; break; }
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
static int rc_pos_int_roots(Expr* R, Expr* z, long* roots, size_t cap, size_t* n) {
    *n = 0;
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
    Expr* fl = rt_eval1("FactorList", prod);
    if (!fl || !rt_head_is(fl, "List")) { if (fl) expr_free(fl); expr_free(vars); return -1; }
    size_t cap = fl->data.function.arg_count + (size_t)(lo + 2) + 1;
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
Expr* rde_tower(Expr* f, Expr* g, RdeCtx* C) {
    if (rt_is_zero(g)) return expr_new_integer(0);
    if (C->m < 0) return rde_base(f, g, C->x);
    RtKind kind = C->T->kind[C->m];
    /* Tangent top: pre-empted — every reachable hypertangent RDE obligation is
     * discharged by the §5.10 driver's OWN coupled DE system (risch_coupled_desystem,
     * reachable over towers and verified), so this branch is never the deciding path
     * for a constructible elementary integrand.  See RISCH_BRONSTEIN_GAPS.md / the
     * RT_TAN piece-(c) determination.  Kept as an authoritative decline (Gap 3). */
    if (kind != RT_LOG && kind != RT_EXP) { rde_cert_mark("tower: tangent top"); return NULL; }
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
    if (kind == RT_EXP) {
        Expr *na, *nb, *nc, *nhs;
        rde_special_denom_exp(a, b, c, C, &na, &nb, &nc, &nhs);
        expr_free(a); expr_free(b); expr_free(c); expr_free(h_s);
        a = na; b = nb; c = nc; h_s = nhs;
    }
    Expr* result = NULL;
    if (!(rc_ispoly(a, C) && rc_ispoly(b, C) && rc_ispoly(c, C)))
        rde_cert_mark("tower: reduced equation not polynomial");
    if (rc_ispoly(a, C) && rc_ispoly(b, C) && rc_ispoly(c, C)) {
        long n = (kind == RT_EXP) ? rde_bound_degree_exp(a, b, c, C)
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
                else {
                    Expr* rem = NULL;
                    Expr* Hc = rt_field_integrate(sp.c, C->T, C->m, C->x, &rem);
                    bool ok = Hc && (!rem || rt_is_zero(rem)) && rc_is_tower_rational(Hc, C);
                    if (rem) expr_free(rem);
                    if (ok) H = Hc;
                    else if (Hc) expr_free(Hc);
                    if (!H) rde_cert_mark("tower: antidifferentiation D H = c");
                }
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
        if (rexp) expr_free(rexp);
        if (!good) { rde_cert_mark("tower: identity gate"); expr_free(result); result = NULL; }
    }
    return result;
}
