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
/* Field quotient a/b, cancelled (args not consumed). */
static Expr* rde_divq(Expr* a, Expr* b) {
    return rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(b), expr_new_integer(-1) }, 2) }, 2));
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
    if (!s) { if (eg) expr_free(eg); expr_free(a1); expr_free(b1); expr_free(c1); return 0; }
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

/* WeakNormalizer(f, D=d/dx) in the base field (Bronstein p.183).  Returns w in
 * k[x] such that f - Dw/w is weakly normalized (no normal irreducible with a
 * positive-integer residue), so no spurious logarithm survives.  In the base
 * field SplitFactor is the identity (every squarefree factor is normal).  The
 * common case (no positive-integer residue) returns w = 1. */
static Expr* rde_weak_normalizer(Expr* f, Expr* x) {
    Expr* den = rde_den(f);
    Expr* num = rde_num(f);
    Expr* dden = rde_dx(den, x);
    Expr* g = rde_gcd(den, dden, x);                /* gcd(dn, dn') */
    Expr* dstar = rde_quot(den, g, x);              /* d* = dn / g  */
    Expr* gds = rde_gcd(dstar, g, x);
    Expr* d1 = rde_quot(dstar, gds, x);             /* d1 = d* / gcd(d*, g) */
    expr_free(dden); expr_free(g); expr_free(dstar); expr_free(gds);
    Expr* w = expr_new_integer(1);
    long dd1 = rt_degree(d1, x);
    if (dd1 >= 1) {
        /* a with (den/d1) a  ==  num  (mod d1), deg(a) < deg(d1). */
        Expr* A = rde_quot(den, d1, x);             /* den / d1 */
        Expr* eg = rt_eval3("PolynomialExtendedGCD", expr_copy(A), expr_copy(d1), expr_copy(x));
        Expr* s = NULL; Expr* gg = NULL;
        if (eg && eg->type == EXPR_FUNCTION && eg->data.function.arg_count == 2) {
            gg = eg->data.function.args[0];
            Expr* cof = eg->data.function.args[1];
            if (cof && cof->type == EXPR_FUNCTION && cof->data.function.arg_count == 2)
                s = cof->data.function.args[0];
        }
        /* Only the common gg == 1 case yields the clean residue formula; else
         * leave w = 1 (weak normalization is a no-op that at worst declines). */
        if (s && gg && gg->type == EXPR_INTEGER && gg->data.integer == 1) {
            Expr* sa = rde_mul(s, num);
            Expr* a = rde_rem(sa, d1, x);           /* a = s*num mod d1 */
            expr_free(sa);
            /* r = resultant_x(a - z Dd1, d1), roots z = positive integers. */
            Expr* dd1poly = rde_dx(d1, x);
            Expr* zsym = expr_new_symbol("rmWNz");
            Expr* azd = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), zsym, dd1poly }, 3) }, 2));  /* adopts */
            Expr* res = rt_eval3("Resultant", azd, expr_copy(d1), expr_copy(x));
            /* Positive integer roots of res in rmWNz. */
            Expr* roots = rt_eval2("Solve", expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ res, expr_new_integer(0) }, 2), expr_new_symbol("rmWNz"));
            if (roots && roots->type == EXPR_FUNCTION
                && roots->data.function.head->type == EXPR_SYMBOL
                && roots->data.function.head->data.symbol.name == intern_symbol("List")) {
                for (size_t ri = 0; ri < roots->data.function.arg_count; ri++) {
                    Expr* sol = roots->data.function.args[ri];   /* {rmWNz -> val} */
                    if (!sol || sol->type != EXPR_FUNCTION
                        || sol->data.function.arg_count < 1) continue;
                    Expr* rule = sol->data.function.args[0];
                    if (!rt_head_is(rule, "Rule") || rule->data.function.arg_count != 2) continue;
                    Expr* val = rule->data.function.args[1];
                    long nv;
                    if (val && val->type == EXPR_INTEGER && (nv = (long)val->data.integer) >= 1) {
                        /* w *= gcd(a - nv*Dd1, d1)^nv. */
                        Expr* nd = rde_dx(d1, x);
                        Expr* an = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Times"),
                                (Expr*[]){ expr_new_integer(-nv), nd }, 2) }, 2));  /* adopts nd */
                        Expr* gpow = rde_gcd(an, d1, x);
                        expr_free(an);
                        for (long e = 0; e < nv; e++) {
                            Expr* nw = rde_mul(w, gpow);
                            expr_free(w); w = nw;
                        }
                        expr_free(gpow);
                    }
                }
            }
            if (roots) expr_free(roots);
            expr_free(a);
        }
        if (eg) expr_free(eg);
        expr_free(A);
    }
    expr_free(d1); expr_free(den); expr_free(num);
    return w;
}

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
    if (rt_is_zero(w)) { expr_free(w); return NULL; }
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
    return (C->m < 0) ? rde_dx(e, C->x) : rt_tower_deriv(e, C->T, C->x);
}
static Expr* rc_quot(Expr* a, Expr* b, RdeCtx* C) { return rde_quot(a, b, C->tau); }
static Expr* rc_rem(Expr* a, Expr* b, RdeCtx* C)  { return rde_rem(a, b, C->tau); }
static Expr* rc_xgcd(Expr* a, Expr* b, RdeCtx* C) {
    return rt_eval3("PolynomialExtendedGCD", expr_copy(a), expr_copy(b), expr_copy(C->tau));
}
/* gcd in tau over the coefficient field k(x).  At base (tau = x over Q) the
 * multivariate PolynomialGCD IS the field gcd; for m >= 0 it is content-including,
 * so compute the monic Euclidean gcd in tau explicitly (R1). */
static Expr* rc_gcd(Expr* a, Expr* b, RdeCtx* C) {
    if (C->m < 0) return rde_gcd(a, b, C->x);
    Expr* u = rt_eval1("Cancel", expr_copy(a));
    Expr* v = rt_eval1("Cancel", expr_copy(b));
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
static Expr* rc_num(Expr* e, RdeCtx* C) {
    if (C->m < 0) return rde_num(e);
    Expr* d = rde_den(e);
    long dd = rt_degree(d, C->tau);
    if (dd <= 0) { expr_free(d); return rt_eval1("Cancel", expr_copy(e)); }  /* den == 1 */
    Expr* lc = rt_coeff(d, C->tau, dd);
    Expr* nm = rde_num(e);
    Expr* num = rde_divq(nm, lc);                   /* nm/lc matches den = d/lc */
    expr_free(d); expr_free(lc); expr_free(nm);
    return num;
}

/* RdeNormalDenominator over k(x)(tau) for a PRIMITIVE (RT_LOG) or base monomial:
 * reduce D y + f y = g to a Dq + b q = c with q = y h in k(x)[tau].  Primitive
 * monomials have no special denominator (Bronstein §6.2 p.188), so this mirrors
 * the base RdeNormalDenominator inlined in rde_base, lifted to the tau-field ops.
 * Returns 1 and fills *a,*b,*c,*h (all owned) on success; 0 for the authoritative
 * "no solution" (en does not divide dn h^2, Cor 6.1.1(ii)). */
static int rde_normal_denominator_field(Expr* f, Expr* g, RdeCtx* C,
                                        Expr** a_out, Expr** b_out,
                                        Expr** c_out, Expr** h_out) {
    Expr* dn = rc_den(f, C);
    Expr* en = rc_den(g, C);
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
    expr_free(rem); expr_free(h2);
    if (!ok) { expr_free(dn); expr_free(en); expr_free(h); expr_free(dnh2); return 0; }
    /* a = dn h ; b = dn h f - dn Dh ; c = dn h^2 g (= (dn h^2 / en)(en g)). */
    Expr* num_f = rc_num(f, C);                     /* dn f */
    Expr* num_g = rc_num(g, C);                     /* en g */
    Expr* aa = rde_mul(dn, h);
    Expr* dh = rc_D(h, C);
    Expr* dndh = rde_mul(dn, dh);
    Expr* hnf = rde_mul(h, num_f);                  /* dn h f = h (dn f) */
    Expr* bb = rde_sub(hnf, dndh);
    Expr* dnh2_en = rc_quot(dnh2, en, C);           /* dn h^2 / en (exact) */
    Expr* cc = rde_mul(dnh2_en, num_g);             /* (dn h^2/en)(en g) = dn h^2 g */
    expr_free(num_f); expr_free(num_g); expr_free(dh); expr_free(dndh);
    expr_free(hnf); expr_free(dnh2_en); expr_free(dnh2);
    expr_free(dn); expr_free(en);
    *a_out = aa; *b_out = bb; *c_out = cc; *h_out = h;
    return 1;
}

/* RdeBoundDegree for a primitive (RT_LOG) / base monomial: leading-degree bound on
 * deg_tau(q) for a Dq + b q = c (Bronstein §6.3, RdeBoundDegreeBase/Prim; for a
 * primitive tau, delta(tau) = 0, so the base formula applies).  A proven UPPER
 * bound, so an over-estimate only wastes SPDE/NoCancel iterations. */
static long rde_bound_degree_prim(Expr* a, Expr* b, Expr* c, RdeCtx* C) {
    long da = rc_deg(a, C), db = rc_deg(b, C), dc = rc_deg(c, C);
    long mx = (db > da - 1) ? db : (da - 1);
    long n = dc - mx; if (n < 0) n = 0;
    if (db == da - 1 && da >= 1) {                  /* cancellation sub-case */
        Expr* lcb = rc_lc(b, C); Expr* lca = rc_lc(a, C);
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
    if (!s) { if (eg) expr_free(eg); expr_free(a1); expr_free(b1); expr_free(c1); return 0; }
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
        expr_free(cc); cc = t2;
        expr_free(dp); expr_free(bp); expr_free(t1); expr_free(p);
        n = m - 1;
    }
    expr_free(cc); expr_free(lcb);
    if (fail) { expr_free(q); return NULL; }
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
    if (C->T->kind[C->m] != RT_LOG) return NULL;    /* EXP top -> 1b */

    Expr *a, *b, *c, *h;
    if (!rde_normal_denominator_field(f, g, C, &a, &b, &c, &h)) return NULL;
    Expr* result = NULL;
    if (rc_ispoly(a, C) && rc_ispoly(b, C) && rc_ispoly(c, C)) {
        long n = rde_bound_degree_prim(a, b, c, C);
        RdeSpde sp;
        if (rde_spde_field(a, b, c, C, n, &sp)) {
            Expr* H = NULL;
            if (rt_is_zero(sp.b)) {
                /* D H = sp.c.  sp.c == 0 -> H = 0 (so q = beta); a nonzero sp.c is
                 * primitive antidifferentiation, which needs LimitedIntegrate
                 * (Gap 2) — decline. */
                if (rt_is_zero(sp.c)) H = expr_new_integer(0);
            } else if (rc_deg(sp.b, C) >= 1) {
                H = rde_polyrischde_nocancel_field(sp.b, sp.c, C, sp.m);
            }
            /* deg_tau(sp.b) == 0 (b in k -> PolyRischDECancelPrim, 1c) declines. */
            if (H) {
                Expr* aH = rde_mul(sp.alpha, H);
                Expr* q = rde_add(aH, sp.beta);
                expr_free(aH); expr_free(H);
                result = rde_divq(q, h);            /* y = q / h */
                expr_free(q);
            }
            rde_spde_free(&sp);
        }
    }
    expr_free(a); expr_free(b); expr_free(c); expr_free(h);

    /* Correct-by-construction gate: D_tower[y] + f y - g must be exactly 0. */
    if (result) {
        Expr* Dy = rc_D(result, C);
        Expr* fy = rde_mul(f, result);
        Expr* mg = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(g) }, 2));
        Expr* resid = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ Dy, fy, mg }, 3);
        Expr* tog = rt_eval1("Together", resid);
        Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;
        bool good = rnum && rt_is_zero(rnum);
        if (rnum) expr_free(rnum);
        if (!good) { expr_free(result); result = NULL; }
    }
    return result;
}
