/*
 * sum_converge.c -- convergence of an infinite series / product, shared by the
 * Sum and Product dispatchers (and, later, NSum / NProduct).
 *
 *   sum_series_convergence(f, k)  -> SUM_CONV_YES / _NO / _UNKNOWN
 *   product_convergence(f, k)     -> PROD_CONV_YES / _DIV / _ZERO / _UNKNOWN
 *
 * f is the EVALUATED term at a symbolic index k (the caller shadows k); the
 * question is about the tail k -> Infinity, so the lower bound never matters.
 *
 * The closed-form stages each sum an analytic continuation; none of them can
 * tell a divergent series from a convergent one by itself, and an answer like
 * Sum[2^k, {k, 0, Infinity}] = -1 is exactly what a continuation gives.  This
 * module is the one place that decides, so every infinite stage is gated by
 * the same tests and the reference's Sum::div / Product::div come from one
 * emitter.  UNKNOWN never blocks a stage: a symbolic parameter (x^k, k^-s)
 * closes unconditionally, as in the reference (1/(1 - x), Zeta[s]).
 *
 * Series tests, in order:
 *   1. A constant term: 0 converges, a numeric nonzero one diverges.
 *   2. A sum of terms: all convergent -> converges; one divergent among
 *      convergent ones -> diverges; otherwise unknown.
 *   3. Hypergeometric terms (rational in k, c^k, Factorial / Gamma / Binomial /
 *      Pochhammer of linear arguments, integer powers of these): the term
 *      ratio r(k) is built EXACTLY factor by factor (no Simplify -- the cost
 *      trap documented in sum_gosper.c), then the ratio test with Raabe's test
 *      on the unit circle: r ~ L (1 - a/k): |L| < 1 converges, |L| > 1
 *      diverges; L = 1 converges iff a > 1; another unimodular L iff a > 0
 *      (|t_k| ~ k^-a decays monotonically, Dirichlet / Leibniz).
 *   4. Algebraic-log terms (k^alpha, Log[P(k)]^m, (-1)^k, c^k, bounded
 *      Sin / Cos of a linear argument): the p-test / Bertrand's test on the
 *      growth k^ad Log[k]^lg; with an alternating or oscillating factor only
 *      the term test (terms not tending to 0 diverge) is decisive.
 *
 * Product tests: Prod f converges to a nonzero limit iff f -> 1 and Sum (f-1)
 * converges (with Sum (f-1)^2 when f - 1 alternates).  The verdicts follow
 * the reference: a product that tends to 0 IS evaluated (to 0, no message);
 * one that tends to Infinity or has no limit is Product::div.
 */

#include "sum_internal.h"
#include "sum.h"
#include "sym_names.h"
#include "eval.h"
#include "expr.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define SC_MAX_FACTORIAL_SLOPE 16     /* Factorial[a k + b] ratio: a factors */
#define SC_LIMIT_SECONDS 2            /* budget of one internal Limit call  */

static bool sc_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

static bool sc_is_int(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

/* Numeric value of a var-free expression as re + i im; false when not numeric. */
static bool sc_num(Expr* e, double* re, double* im) {
    Expr* n = sum_eval("N", (Expr*[]){ expr_copy(e) }, 1);
    bool ok = false;
    *re = 0; *im = 0;
    if (n && n->type == EXPR_REAL) { *re = n->data.real; ok = true; }
    else if (n && n->type == EXPR_INTEGER) { *re = (double)n->data.integer; ok = true; }
    else if (sc_head(n, SYM_Complex) && n->data.function.arg_count == 2) {
        Expr* a = n->data.function.args[0];
        Expr* b = n->data.function.args[1];
        if ((a->type == EXPR_REAL || a->type == EXPR_INTEGER)
            && (b->type == EXPR_REAL || b->type == EXPR_INTEGER)) {
            *re = a->type == EXPR_REAL ? a->data.real : (double)a->data.integer;
            *im = b->type == EXPR_REAL ? b->data.real : (double)b->data.integer;
            ok = true;
        }
    }
    if (n) expr_free(n);
    return ok && isfinite(*re) && isfinite(*im);
}

static bool sc_polyq(Expr* e, Expr* k) {
    Expr* r = sum_eval("PolynomialQ", (Expr*[]){ expr_copy(e), expr_copy(k) }, 2);
    bool yes = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True;
    if (r) expr_free(r);
    return yes;
}

/* Degree of a polynomial in k (-1 when not a machine integer). */
static int64_t sc_degree(Expr* p, Expr* k) {
    Expr* d = sum_eval("Exponent", (Expr*[]){ expr_copy(p), expr_copy(k) }, 2);
    int64_t v = (d && d->type == EXPR_INTEGER) ? d->data.integer : -1;
    if (d) expr_free(d);
    return v;
}

static Expr* sc_coeff(Expr* p, Expr* k, int64_t j) {
    Expr* ex = sum_eval("Expand", (Expr*[]){ expr_copy(p) }, 1);
    return sum_eval("Coefficient", (Expr*[]){ ex, expr_copy(k), sum_int(j) }, 3);
}

/* a k + b with a, b free of k: sets owned *a, *b. */
static bool sc_linear(Expr* e, Expr* k, Expr** a, Expr** b) {
    if (!sc_polyq(e, k) || sc_degree(e, k) != 1) return false;
    *a = sc_coeff(e, k, 1);
    *b = sc_coeff(e, k, 0);
    return true;
}

/* ------------------------------------------------------------------ */
/*  Exact term ratio of a hypergeometric term                          */
/* ------------------------------------------------------------------ */

static Expr* sc_shift(Expr* e, Expr* k) {
    Expr* kp1 = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(k), sum_int(1) }, 2);
    Expr* r = sum_subst(e, k, kp1);
    expr_free(kp1);
    return r;
}

static Expr* sc_div(Expr* a, Expr* b) {   /* consumes a, b */
    return sum_eval("Times", (Expr*[]){ a,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ b, sum_int(-1) }, 2) }, 2);
}

/* Prod_{j=j0}^{j1} (x + j) for owned x (consumed); j0 <= j1. */
static Expr* sc_rising(Expr* x, int64_t j0, int64_t j1) {
    size_t n = (size_t)(j1 - j0 + 1), w = 0;
    Expr** f = malloc(sizeof(Expr*) * n);
    for (int64_t j = j0; j <= j1; j++)
        f[w++] = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ expr_copy(x), sum_int(j) }, 2);
    expr_free(x);
    Expr* t = expr_new_function(expr_new_symbol(SYM_Times), f, w);
    free(f);
    return t;
}

static Expr* sc_ratio(Expr* f, Expr* k);

/* Gamma[x] with x = a k + b (a a positive machine integer): the ratio
 * Gamma[x + a]/Gamma[x] = Prod_{j=0}^{a-1} (x + j); `shift0` = 1 for a
 * Factorial (x! = Gamma[x + 1]). */
static Expr* sc_gamma_ratio(Expr* x, Expr* k, int64_t shift0) {
    Expr *a = NULL, *b = NULL;
    if (!sc_linear(x, k, &a, &b)) return NULL;
    bool ok = a->type == EXPR_INTEGER && a->data.integer >= 1
           && a->data.integer <= SC_MAX_FACTORIAL_SLOPE;
    int64_t av = ok ? a->data.integer : 0;
    expr_free(a); expr_free(b);
    if (!ok) return NULL;
    return sc_rising(expr_copy(x), shift0, shift0 + av - 1);
}

static Expr* sc_ratio_factor(Expr* e, Expr* k) {
    if (sum_free_of(e, k)) return sum_int(1);
    if (sc_polyq(e, k)) return sc_div(sc_shift(e, k), expr_copy(e));
    if (sc_head(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex = e->data.function.args[1];
        if (ex->type == EXPR_INTEGER) {
            Expr* rb = sc_ratio(base, k);
            if (!rb) return NULL;
            return sum_eval("Power", (Expr*[]){ rb, expr_copy(ex) }, 2);
        }
        if (sum_free_of(base, k)) {             /* c^(a k + b): ratio c^a */
            Expr *a = NULL, *b = NULL;
            if (!sc_linear(ex, k, &a, &b)) return NULL;
            expr_free(b);
            return sum_eval("Power", (Expr*[]){ expr_copy(base), a }, 2);
        }
        return NULL;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count < 1) return NULL;
    Expr* x = e->data.function.args[0];
    if (sc_head(e, SYM_Factorial) && e->data.function.arg_count == 1)
        return sc_gamma_ratio(x, k, 1);
    if (sc_head(e, SYM_Gamma) && e->data.function.arg_count == 1)
        return sc_gamma_ratio(x, k, 0);
    if (sc_head(e, SYM_Pochhammer) && e->data.function.arg_count == 2
        && sum_free_of(x, k)) {                 /* (c)_m = Gamma[c + m]/Gamma[c] */
        Expr* cm = sum_eval("Plus", (Expr*[]){ expr_copy(x), expr_copy(e->data.function.args[1]) }, 2);
        Expr* r = sc_gamma_ratio(cm, k, 0);
        expr_free(cm);
        return r;
    }
    if (sc_head(e, SYM_Binomial) && e->data.function.arg_count == 2) {
        Expr* y = e->data.function.args[1];
        Expr* xmy = sum_eval("Plus", (Expr*[]){ expr_copy(x),
                        expr_new_function(expr_new_symbol(SYM_Times),
                            (Expr*[]){ sum_int(-1), expr_copy(y) }, 2) }, 2);
        /* x!/(y! (x - y)!), each factorial with a k-free or rising argument */
        Expr* parts[3] = { x, y, xmy };
        Expr* acc = sum_int(1);
        for (int i = 0; i < 3 && acc; i++) {
            Expr* r = sum_free_of(parts[i], k) ? sum_int(1) : sc_gamma_ratio(parts[i], k, 1);
            if (!r) { expr_free(acc); acc = NULL; break; }
            acc = i == 0 ? sum_eval("Times", (Expr*[]){ acc, r }, 2) : sc_div(acc, r);
        }
        expr_free(xmy);
        return acc;
    }
    return NULL;
}

/* Exact rational term ratio f(k+1)/f(k), or NULL when f is not recognised as
 * a hypergeometric term. */
static Expr* sc_ratio(Expr* f, Expr* k) {
    if (sc_head(f, SYM_Times)) {
        Expr* acc = sum_int(1);
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            Expr* r = sc_ratio_factor(f->data.function.args[i], k);
            if (!r) { expr_free(acc); return NULL; }
            acc = sum_eval("Times", (Expr*[]){ acc, r }, 2);
        }
        return acc;
    }
    if (sc_head(f, SYM_Plus)) {
        /* a rational function written as a sum of fractions */
        Expr* tg = sum_eval("Together", (Expr*[]){ expr_copy(f) }, 1);
        Expr* num = sum_eval("Numerator", (Expr*[]){ expr_copy(tg) }, 1);
        Expr* den = sum_eval("Denominator", (Expr*[]){ tg }, 1);
        Expr* r = NULL;
        if (sc_polyq(num, k) && sc_polyq(den, k))
            r = sc_div(sc_div(sc_shift(num, k), expr_copy(num)),
                       sc_div(sc_shift(den, k), expr_copy(den)));
        expr_free(num); expr_free(den);
        return r;
    }
    return sc_ratio_factor(f, k);
}

/* Ratio test + Raabe on r = N/D. */
static SumConvergence sc_ratio_test(Expr* ratio, Expr* k) {
    Expr* tg = sum_eval("Together", (Expr*[]){ expr_copy(ratio) }, 1);
    Expr* num = sum_eval("Numerator", (Expr*[]){ expr_copy(tg) }, 1);
    Expr* den = sum_eval("Denominator", (Expr*[]){ tg }, 1);
    SumConvergence v = SUM_CONV_UNKNOWN;
    if (!sc_polyq(num, k) || !sc_polyq(den, k)) goto done;
    int64_t p = sc_degree(num, k), q = sc_degree(den, k);
    if (p < 0 || q < 0) goto done;
    if (p > q) { v = SUM_CONV_NO; goto done; }        /* |r| -> Infinity */
    if (p < q) { v = SUM_CONV_YES; goto done; }       /* r -> 0          */
    {
        Expr* lcN = sc_coeff(num, k, p);
        Expr* lcD = sc_coeff(den, k, q);
        Expr* L = sum_eval("Times", (Expr*[]){ expr_copy(lcN),
            expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(lcD), sum_int(-1) }, 2) }, 2);
        double lr, li;
        bool numeric = sc_num(L, &lr, &li);
        expr_free(L);
        if (numeric) {
            double m = hypot(lr, li);
            if (m < 1.0 - 1e-12) v = SUM_CONV_YES;
            else if (m > 1.0 + 1e-12) v = SUM_CONV_NO;
            else if (p >= 1) {
                /* r = L (1 - a/k + ...), a = sub_D/lc_D - sub_N/lc_N */
                Expr* sN = sc_coeff(num, k, p - 1);
                Expr* sD = sc_coeff(den, k, q - 1);
                Expr* a = sum_eval("Plus", (Expr*[]){
                    sc_div(sD, expr_copy(lcD)),
                    expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ sum_int(-1), sc_div(sN, expr_copy(lcN)) }, 2) }, 2);
                double ar, ai;
                if (sc_num(a, &ar, &ai) && fabs(ai) < 1e-12) {
                    bool l_is_one = fabs(lr - 1.0) < 1e-12 && fabs(li) < 1e-12;
                    double crit = l_is_one ? 1.0 : 0.0;
                    v = ar > crit + 1e-12 ? SUM_CONV_YES : SUM_CONV_NO;
                    /* a exactly at the threshold with L = 1 (Sum 1/(k log k)-
                     * like terms do not reach here: they are not
                     * hypergeometric), a = 1 is the harmonic series. */
                }
                expr_free(a);
            } else {
                /* constant ratio on the unit circle: |t_k| constant, no decay */
                v = SUM_CONV_NO;
            }
        }
        expr_free(lcN); expr_free(lcD);
    }
done:
    expr_free(num); expr_free(den);
    return v;
}

/* ------------------------------------------------------------------ */
/*  Algebraic-log growth: |t_k| ~ C g^k k^ad Log[k]^lg                 */
/* ------------------------------------------------------------------ */

typedef struct {
    double ad, lg;      /* algebraic and logarithmic orders     */
    double geo;         /* |c| of a geometric factor c^k (1: none) */
    bool alt;           /* a (-1)^k factor                       */
    bool osc;           /* a bounded non-decaying oscillation    */
} ScGrowth;

static bool sc_growth(Expr* e, Expr* k, ScGrowth* g);

static bool sc_growth_power(Expr* base, Expr* ex, Expr* k, ScGrowth* g) {
    if (sum_free_of(base, k)) {                 /* c^(a k + b) */
        Expr *a = NULL, *b = NULL;
        if (!sc_linear(ex, k, &a, &b)) return false;
        expr_free(b);
        Expr* ca = sum_eval("Power", (Expr*[]){ expr_copy(base), a }, 2);
        double re, im;
        bool ok = sc_num(ca, &re, &im);
        expr_free(ca);
        if (!ok) return false;
        double m = hypot(re, im);
        if (fabs(m - 1.0) > 1e-12) { g->geo *= m; return true; }
        if (fabs(re + 1.0) < 1e-12 && fabs(im) < 1e-12) g->alt = true;
        else if (!(fabs(re - 1.0) < 1e-12 && fabs(im) < 1e-12)) g->osc = true;
        return true;
    }
    double er, ei;
    if (!sum_free_of(ex, k) || !sc_num(ex, &er, &ei) || fabs(ei) > 1e-12) return false;
    ScGrowth gb = { 0, 0, 1, false, false };
    if (!sc_growth(base, k, &gb) || gb.alt || gb.osc || gb.geo != 1.0) return false;
    g->ad += er * gb.ad;
    g->lg += er * gb.lg;
    return true;
}

static bool sc_growth(Expr* e, Expr* k, ScGrowth* g) {
    if (sum_free_of(e, k)) return true;
    if (sc_polyq(e, k)) {
        int64_t d = sc_degree(e, k);
        if (d < 0) return false;
        g->ad += (double)d;
        return true;
    }
    if (sc_head(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!sc_growth(e->data.function.args[i], k, g)) return false;
        return true;
    }
    if (sc_head(e, SYM_Power) && e->data.function.arg_count == 2)
        return sc_growth_power(e->data.function.args[0], e->data.function.args[1], k, g);
    if (sc_head(e, SYM_Log) && e->data.function.arg_count == 1) {
        Expr* x = e->data.function.args[0];
        if (!sc_polyq(x, k) || sc_degree(x, k) < 1) return false;
        g->lg += 1.0;
        return true;
    }
    if ((sc_head(e, SYM_Sin) || sc_head(e, SYM_Cos)) && e->data.function.arg_count == 1) {
        Expr *a = NULL, *b = NULL;
        if (!sc_linear(e->data.function.args[0], k, &a, &b)) return false;
        double ar, ai;
        bool ok = sc_num(a, &ar, &ai) && fabs(ai) < 1e-12 && ar != 0.0;
        expr_free(a); expr_free(b);
        if (!ok) return false;
        g->osc = true;
        return true;
    }
    return false;
}

static SumConvergence sc_growth_test(Expr* f, Expr* k) {
    ScGrowth g = { 0, 0, 1, false, false };
    if (!sc_growth(f, k, &g)) return SUM_CONV_UNKNOWN;
    if (g.geo < 1.0 - 1e-12) return SUM_CONV_YES;     /* geometric decay dominates */
    if (g.geo > 1.0 + 1e-12) return SUM_CONV_NO;
    bool no_decay = g.ad > 1e-12 || (fabs(g.ad) <= 1e-12 && g.lg >= -1e-12);
    if (g.osc) return no_decay ? SUM_CONV_NO : SUM_CONV_UNKNOWN;
    if (g.alt) return no_decay ? SUM_CONV_NO : SUM_CONV_YES;
    if (g.ad < -1.0 - 1e-12) return SUM_CONV_YES;      /* p-test   */
    if (g.ad > -1.0 + 1e-12) return SUM_CONV_NO;
    return g.lg < -1.0 - 1e-12 ? SUM_CONV_YES : SUM_CONV_NO;   /* Bertrand */
}

/* ------------------------------------------------------------------ */
/*  Series                                                              */
/* ------------------------------------------------------------------ */

SumConvergence sum_series_convergence(Expr* f, Expr* k) {
    if (!f || k->type != EXPR_SYMBOL) return SUM_CONV_UNKNOWN;
    if (sum_free_of(f, k)) {
        if (sc_is_int(f, 0)) return SUM_CONV_YES;
        double re, im;
        if (sc_num(f, &re, &im) && (re != 0.0 || im != 0.0)) return SUM_CONV_NO;
        return SUM_CONV_UNKNOWN;
    }
    Expr* ratio = sc_ratio(f, k);
    if (ratio) {
        SumConvergence v = sc_ratio_test(ratio, k);
        expr_free(ratio);
        if (v != SUM_CONV_UNKNOWN) return v;
    }
    SumConvergence v = sc_growth_test(f, k);
    if (v != SUM_CONV_UNKNOWN || !sc_head(f, SYM_Plus)) return v;
    /* A sum of terms: all convergent -> converges; exactly one divergent
     * among convergent ones -> diverges (it cannot cancel). */
    size_t ndiv = 0;
    for (size_t i = 0; i < f->data.function.arg_count; i++) {
        Expr* t = f->data.function.args[i];
        SumConvergence vt = SUM_CONV_UNKNOWN;
        Expr* rt = sc_ratio(t, k);
        if (rt) { vt = sc_ratio_test(rt, k); expr_free(rt); }
        if (vt == SUM_CONV_UNKNOWN) vt = sc_growth_test(t, k);
        if (vt == SUM_CONV_UNKNOWN) return SUM_CONV_UNKNOWN;
        if (vt == SUM_CONV_NO) ndiv++;
    }
    return ndiv == 0 ? SUM_CONV_YES : ndiv == 1 ? SUM_CONV_NO : SUM_CONV_UNKNOWN;
}

/* ------------------------------------------------------------------ */
/*  Products                                                            */
/* ------------------------------------------------------------------ */

/* Heads under which Limit is a sound and bounded question. */
static bool sc_limit_safe(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return true;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    const char* n = h->data.symbol.name;
    if (!(n == SYM_Plus || n == SYM_Times || n == SYM_Power || n == SYM_Log
          || n == SYM_Sin || n == SYM_Cos || n == SYM_Tan || n == SYM_Sinh
          || n == SYM_Cosh || n == SYM_Tanh || n == SYM_ArcTan
          || n == SYM_Rational || n == SYM_Complex)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!sc_limit_safe(e->data.function.args[i])) return false;
    return true;
}

/* TimeConstrained[Limit[e, k -> Infinity], budget]; owned result or NULL. */
static Expr* sc_limit(Expr* e, Expr* k) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ expr_copy(k), expr_new_symbol(SYM_Infinity) }, 2);
    Expr* lim = expr_new_function(expr_new_symbol(SYM_Limit),
                    (Expr*[]){ expr_copy(e), rule }, 2);
    Expr* r = sum_eval("TimeConstrained", (Expr*[]){ lim, sum_int(SC_LIMIT_SECONDS),
                                                     expr_new_symbol(SYM_Indeterminate) }, 3);
    return r;
}

static bool sc_is_sym(const Expr* e, const char* s) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == s;
}

/* Classify a numeric limit L of the factor: 1 -> needs the next order. */
static ProductConvergence sc_classify_limit(Expr* L, bool* is_one) {
    *is_one = false;
    if (sc_is_sym(L, SYM_Infinity) || sc_is_sym(L, SYM_ComplexInfinity)
        || sc_head(L, SYM_DirectedInfinity)) return PROD_CONV_DIV;
    if (sc_head(L, SYM_Times) && L->data.function.arg_count == 2
        && sc_is_int(L->data.function.args[0], -1)
        && sc_is_sym(L->data.function.args[1], SYM_Infinity)) return PROD_CONV_DIV;
    double re, im;
    if (!sc_num(L, &re, &im)) return PROD_CONV_UNKNOWN;
    double m = hypot(re, im);
    if (m < 1.0 - 1e-12) return PROD_CONV_ZERO;
    if (fabs(re - 1.0) < 1e-12 && fabs(im) < 1e-12) { *is_one = true; return PROD_CONV_UNKNOWN; }
    return PROD_CONV_DIV;   /* |L| > 1, or unimodular L != 1: no limit */
}

/* Prod (1 + g) with g -> 0: decided by the series of g. */
static ProductConvergence sc_product_tail(Expr* f, Expr* k) {
    Expr* g = sum_eval("Plus", (Expr*[]){ expr_copy(f), sum_int(-1) }, 2);
    ProductConvergence v = PROD_CONV_UNKNOWN;
    SumConvergence sg = sum_series_convergence(g, k);
    if (sg == SUM_CONV_YES) {
        /* An alternating g also needs Sum g^2 (Prod (1 + (-1)^k/Sqrt[k]) -> 0). */
        Expr* g2 = sum_eval("Expand", (Expr*[]){ sum_eval("Power", (Expr*[]){ expr_copy(g), sum_int(2) }, 2) }, 1);
        SumConvergence s2 = sum_series_convergence(g2, k);
        expr_free(g2);
        v = s2 == SUM_CONV_YES ? PROD_CONV_YES
          : s2 == SUM_CONV_NO ? PROD_CONV_ZERO : PROD_CONV_UNKNOWN;
    }
    if (v == PROD_CONV_UNKNOWN && sc_limit_safe(g)) {
        /* g ~ lambda/k: lambda > 0 -> Infinity, lambda < 0 -> 0 */
        Expr* kg = sum_eval("Times", (Expr*[]){ expr_copy(k), expr_copy(g) }, 2);
        Expr* lam = sc_limit(kg, k);
        expr_free(kg);
        double re, im;
        if (lam && sc_num(lam, &re, &im) && fabs(im) < 1e-12 && re != 0.0)
            v = re > 0 ? PROD_CONV_DIV : PROD_CONV_ZERO;
        if (lam) expr_free(lam);
    }
    expr_free(g);
    return v;
}

ProductConvergence product_convergence(Expr* f, Expr* k) {
    if (!f || k->type != EXPR_SYMBOL) return PROD_CONV_UNKNOWN;
    if (sum_free_of(f, k)) {
        if (sc_is_int(f, 1)) return PROD_CONV_YES;
        double re, im;
        if (!sc_num(f, &re, &im)) return PROD_CONV_UNKNOWN;
        return hypot(re, im) < 1.0 - 1e-12 ? PROD_CONV_ZERO : PROD_CONV_DIV;
    }
    /* Rational in k with numeric leading behaviour: decided exactly. */
    Expr* tg = sum_eval("Together", (Expr*[]){ expr_copy(f) }, 1);
    Expr* num = sum_eval("Numerator", (Expr*[]){ expr_copy(tg) }, 1);
    Expr* den = sum_eval("Denominator", (Expr*[]){ tg }, 1);
    ProductConvergence v = PROD_CONV_UNKNOWN;
    if (sc_polyq(num, k) && sc_polyq(den, k)) {
        int64_t p = sc_degree(num, k), q = sc_degree(den, k);
        if (p >= 0 && q >= 0) {
            if (p > q) v = PROD_CONV_DIV;
            else if (p < q) v = PROD_CONV_ZERO;
            else {
                Expr* L = sc_div(sc_coeff(num, k, p), sc_coeff(den, k, q));
                bool one = false;
                v = sc_classify_limit(L, &one);
                expr_free(L);
                if (one) v = sc_product_tail(f, k);
            }
        }
        expr_free(num); expr_free(den);
        return v;
    }
    expr_free(num); expr_free(den);
    if (!sc_limit_safe(f)) return PROD_CONV_UNKNOWN;
    Expr* L = sc_limit(f, k);
    bool one = false;
    if (L) { v = sc_classify_limit(L, &one); expr_free(L); }
    if (one) v = sc_product_tail(f, k);
    return v;
}
