/*
 * series_numeric.c -- numeric Taylor expansion of an inexact analytic
 * expression (see series_numeric.h for the contract).
 *
 * A series is N coefficients t^0 .. t^(N-1) of the expansion in t = x - x0,
 * each an mpfr_t at the working precision wp, plus a double "majorant" m[k]:
 * a bound on the magnitude of every term that was summed into coefficient k
 * (so m[k] >= |c[k]|, and the rounding error of c[k] is a small multiple of
 * m[k] 2^-wp).  Every operation propagates it with the same formula run on
 * magnitudes.  The two classifications the engine needs read it:
 *
 *   zero      |c| <= m 2^-(wp - slack): indistinguishable from 0 -- the exact
 *             expansion's exact zero (Sin[x] - x at order 1, the constant term
 *             of 2^(z/2 + 1) - 2), returned as the numericalised 0 the exact
 *             path returns;
 *   resolved  |c| >= m 2^-(wp - slack - tb - 2): known to the answer's tb bits.
 *
 * Anything in between declines the whole expansion to the exact round trip,
 * as does a division by, a logarithm of, or a fractional power of a series
 * whose leading coefficient is zero or unresolved (a pole or branch point is a
 * Laurent / Puiseux expansion, which the exact engine owns).
 *
 * Every series here is a plain Taylor series (valuation >= 0): the covered
 * operations -- Plus, Times, Power, Exp, Log, Sin, Cos, Sinh, Cosh -- never
 * create a negative power once those leading-coefficient cases decline, so no
 * order is lost and no padding is needed.  Powers use J. C. P. Miller's
 * recurrence (O(N^2) for any exponent, where binary powering is
 * O(N^2 log k)); Exp / Log / Sin / Cos the standard ODE recurrences.
 */

#include "series_numeric.h"
#include "rationalize.h"   /* internal_force_rationalize_bits, internal_contains_inexact */
#include "numeric.h"       /* numericalize, NumericSpec */
#include "eval.h"          /* eval_and_free */
#include "sym_names.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <float.h>

#ifdef USE_MPFR
#include <mpfr.h>

/* Largest expansion order handled numerically (an order past it declines). */
#define TS_MAX_TERMS 4096

typedef struct {
    size_t      N;      /* number of coefficients */
    mpfr_prec_t wp;     /* working precision */
    mpfr_prec_t tb;     /* the answer's precision (53 for machine) */
    long        slack;  /* rounding-error headroom, bits */
    const Expr* x;      /* the expansion variable */
    mpfr_t      x0;     /* the expansion point */
    bool        ok;     /* cleared by any decline */
} TsCtx;

typedef struct { mpfr_t* c; double* m; } Ts;

static Ts* ts_new(TsCtx* C) {
    Ts* a = malloc(sizeof(Ts));
    if (!a) { C->ok = false; return NULL; }
    a->c = malloc(C->N * sizeof(mpfr_t));
    a->m = calloc(C->N, sizeof(double));
    if (!a->c || !a->m) { free(a->c); free(a->m); free(a); C->ok = false; return NULL; }
    for (size_t k = 0; k < C->N; k++) {
        mpfr_init2(a->c[k], C->wp);
        mpfr_set_zero(a->c[k], 1);
    }
    return a;
}

static void ts_free(const TsCtx* C, Ts* a) {
    if (!a) return;
    for (size_t k = 0; k < C->N; k++) mpfr_clear(a->c[k]);
    free(a->c); free(a->m); free(a);
}

static double ts_absd(const mpfr_t v) { return fabs(mpfr_get_d(v, MPFR_RNDU)); }

/* Is coefficient k indistinguishable from zero? */
static bool ts_zero(TsCtx* C, const Ts* a, size_t k) {
    if (mpfr_zero_p(a->c[k])) return true;
    double m = a->m[k];
    if (!isfinite(m)) { C->ok = false; return false; }
    if (m <= 0.0) return false;
    int me; frexp(m, &me);                          /* m < 2^me */
    return (long)mpfr_get_exp(a->c[k]) <= (long)me - ((long)C->wp - C->slack);
}

/* Is coefficient k known to the answer's precision? (false for a zero) */
static bool ts_resolved(TsCtx* C, const Ts* a, size_t k) {
    if (mpfr_zero_p(a->c[k])) return false;
    double m = a->m[k];
    if (!isfinite(m)) { C->ok = false; return false; }
    if (m <= 0.0) return true;
    int me; frexp(m, &me);
    return (long)mpfr_get_exp(a->c[k]) >
           (long)me - ((long)C->wp - C->slack) + (long)C->tb + 2;
}

static bool ts_free_of(const Expr* e, const Expr* x) {
    if (e->type == EXPR_SYMBOL) return !expr_eq((Expr*)e, (Expr*)x);
    if (e->type != EXPR_FUNCTION) return true;
    if (!ts_free_of(e->data.function.head, x)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!ts_free_of(e->data.function.args[i], x)) return false;
    return true;
}

static bool ts_head_is(const Expr* e, const char* sym) {
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

/* The real value of an x-free numeric expression at the working precision.
 * Exact integers and rationals convert directly; anything else goes through
 * N at wp bits and must come back a real number. */
static bool ts_numeric_value(TsCtx* C, const Expr* e, mpfr_t out) {
    if (e->type == EXPR_INTEGER) { mpfr_set_si(out, (long)e->data.integer, MPFR_RNDN); return true; }
    if (e->type == EXPR_BIGINT)  { mpfr_set_z(out, e->data.bigint, MPFR_RNDN); return true; }
    if (ts_head_is(e, SYM_Rational) && e->data.function.arg_count == 2) {
        const Expr* pe = e->data.function.args[0];
        const Expr* qe = e->data.function.args[1];
        if (!expr_is_integer_like((Expr*)pe) || !expr_is_integer_like((Expr*)qe)) return false;
        mpz_t p, q;
        mpz_init(p); mpz_init(q);
        expr_to_mpz(pe, p);
        expr_to_mpz(qe, q);
        bool ok = mpz_sgn(q) != 0;
        if (ok) {
            mpfr_set_z(out, p, MPFR_RNDN);
            mpfr_div_z(out, out, q, MPFR_RNDN);
        }
        mpz_clear(p); mpz_clear(q);
        return ok;
    }
    NumericSpec spec = numeric_machine_spec();
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = (long)C->wp;
    Expr* v = numericalize(e, spec);
    bool ok = false;
    if (v) {
        if (v->type == EXPR_MPFR)        { mpfr_set(out, v->data.mpfr, MPFR_RNDN); ok = true; }
        else if (v->type == EXPR_REAL)   { mpfr_set_d(out, v->data.real, MPFR_RNDN); ok = true; }
        else if (v->type == EXPR_INTEGER){ mpfr_set_si(out, (long)v->data.integer, MPFR_RNDN); ok = true; }
        else if (v->type == EXPR_BIGINT) { mpfr_set_z(out, v->data.bigint, MPFR_RNDN); ok = true; }
        expr_free(v);
    }
    return ok && mpfr_number_p(out);
}

/* ---- elementary series ---------------------------------------------------- */

static Ts* ts_const(TsCtx* C, const mpfr_t v) {
    Ts* r = ts_new(C);
    if (!r) return NULL;
    mpfr_set(r->c[0], v, MPFR_RNDN);
    r->m[0] = ts_absd(v);
    return r;
}

static Ts* ts_var(TsCtx* C) {
    Ts* r = ts_const(C, C->x0);
    if (r && C->N > 1) { mpfr_set_si(r->c[1], 1, MPFR_RNDN); r->m[1] = 1.0; }
    return r;
}

static Ts* ts_add(TsCtx* C, const Ts* a, const Ts* b) {
    Ts* r = ts_new(C);
    if (!r) return NULL;
    for (size_t k = 0; k < C->N; k++) {
        mpfr_add(r->c[k], a->c[k], b->c[k], MPFR_RNDN);
        r->m[k] = a->m[k] + b->m[k];
    }
    return r;
}

/* s a, for a constant s known to the working precision. */
static Ts* ts_scale(TsCtx* C, const Ts* a, const mpfr_t s) {
    Ts* r = ts_new(C);
    if (!r) return NULL;
    double ms = ts_absd(s);
    for (size_t k = 0; k < C->N; k++) {
        mpfr_mul(r->c[k], a->c[k], s, MPFR_RNDN);
        r->m[k] = a->m[k] * ms;
    }
    return r;
}

static Ts* ts_mul(TsCtx* C, const Ts* a, const Ts* b) {
    Ts* r = ts_new(C);
    if (!r) return NULL;
    mpfr_t t;
    mpfr_init2(t, C->wp);
    for (size_t i = 0; i < C->N; i++) {
        if (mpfr_zero_p(a->c[i])) {
            if (a->m[i] == 0.0) continue;
        }
        for (size_t j = 0; i + j < C->N; j++) {
            if (!mpfr_zero_p(a->c[i]) && !mpfr_zero_p(b->c[j])) {
                mpfr_mul(t, a->c[i], b->c[j], MPFR_RNDN);
                mpfr_add(r->c[i + j], r->c[i + j], t, MPFR_RNDN);
            }
            r->m[i + j] += a->m[i] * b->m[j];
        }
    }
    mpfr_clear(t);
    return r;
}

/* Index of the first coefficient that is not a zero (C->N when none). */
static size_t ts_valuation(TsCtx* C, const Ts* a) {
    for (size_t k = 0; k < C->N; k++)
        if (!ts_zero(C, a, k)) return k;
    return C->N;
}

/* a^r.  `is_int` with exponent k, else the real exponent r.  Miller's
 * recurrence on the unit-leading part: with a = t^v (A0 + A1 t + ...),
 *   b = A^r,  b_0 = A0^r,
 *   b_m = 1/(m A0) sum_{j=1..m} ((r + 1) j - m) A_j b_{m-j},
 * and a^r = t^(r v) b (v > 0 only for a positive integer exponent). */
static Ts* ts_pow_miller(TsCtx* C, const Ts* a, bool is_int, long k, const mpfr_t r) {
    size_t N = C->N;
    size_t v = ts_valuation(C, a);
    if (!C->ok) return NULL;
    if (is_int && k == 0) {
        Ts* one = ts_new(C);
        if (one) { mpfr_set_si(one->c[0], 1, MPFR_RNDN); one->m[0] = 1.0; }
        return one;
    }
    if (v == N) {                                   /* the zero series */
        if (is_int && k > 0) return ts_new(C);
        C->ok = false; return NULL;
    }
    if (v > 0 && !(is_int && k > 0)) { C->ok = false; return NULL; }
    if (is_int && k == 1) {                         /* a itself */
        Ts* id = ts_new(C);
        if (!id) return NULL;
        for (size_t i = 0; i < N; i++) { mpfr_set(id->c[i], a->c[i], MPFR_RNDN); id->m[i] = a->m[i]; }
        return id;
    }
    if (!ts_resolved(C, a, v)) { C->ok = false; return NULL; }
    if (!is_int && mpfr_sgn(a->c[v]) < 0) { C->ok = false; return NULL; }
    if (is_int && (k > (1L << 30) || k < -(1L << 30))) { C->ok = false; return NULL; }

    Ts* res = ts_new(C);
    if (!res) return NULL;
    if (is_int && k > 0 && (uint64_t)v * (uint64_t)k >= N) return res;   /* all beyond the order */
    size_t shift = is_int ? v * (size_t)k : 0;
    size_t L = N - shift;                           /* coefficients of b needed */
    size_t La = N - v;                              /* coefficients of A known */

    mpfr_t* b = res->c + shift;
    double* mb = res->m + shift;
    const mpfr_t* A = (const mpfr_t*)(a->c + v);
    const double* mA = a->m + v;
    double absA0 = ts_absd(A[0]);
    double ra = mA[0] / absA0;                      /* relative uncertainty of A0 */

    mpfr_t rp1, coef, t;
    mpfr_inits2(C->wp, rp1, coef, t, (mpfr_ptr)0);
    if (is_int) mpfr_pow_si(b[0], A[0], k, MPFR_RNDN);
    else        mpfr_pow(b[0], A[0], r, MPFR_RNDN);
    double absr = is_int ? fabs((double)k) : ts_absd(r);
    mb[0] = ts_absd(b[0]) * (1.0 + absr * ra);
    if (!is_int) mpfr_add_si(rp1, r, 1, MPFR_RNDN);

    for (size_t m = 1; m < L; m++) {
        double ms = 0.0;
        size_t jmax = m < La - 1 ? m : La - 1;
        for (size_t j = 1; j <= jmax; j++) {
            double cabs;
            if (is_int) {
                int64_t cj = (int64_t)(k + 1) * (int64_t)j - (int64_t)m;
                if (cj == 0) continue;
                cabs = fabs((double)cj);
                if (!mpfr_zero_p(A[j]) && !mpfr_zero_p(b[m - j])) {
                    mpfr_mul(t, A[j], b[m - j], MPFR_RNDN);
                    mpfr_mul_si(t, t, (long)cj, MPFR_RNDN);
                    mpfr_add(res->c[shift + m], res->c[shift + m], t, MPFR_RNDN);
                }
            } else {
                mpfr_mul_ui(coef, rp1, (unsigned long)j, MPFR_RNDN);
                mpfr_sub_ui(coef, coef, (unsigned long)m, MPFR_RNDN);
                cabs = ts_absd(coef);
                if (!mpfr_zero_p(A[j]) && !mpfr_zero_p(b[m - j])) {
                    mpfr_mul(t, A[j], b[m - j], MPFR_RNDN);
                    mpfr_mul(t, t, coef, MPFR_RNDN);
                    mpfr_add(res->c[shift + m], res->c[shift + m], t, MPFR_RNDN);
                }
            }
            ms += cabs * mA[j] * mb[m - j];
        }
        mpfr_div_ui(b[m], b[m], (unsigned long)m, MPFR_RNDN);
        mpfr_div(b[m], b[m], A[0], MPFR_RNDN);
        mb[m] = ms / ((double)m * absA0) + ts_absd(b[m]) * ra;
    }
    mpfr_clears(rp1, coef, t, (mpfr_ptr)0);
    return res;
}

/* Does every coefficient of a come out zero or resolved? */
static bool ts_all_settled(TsCtx* C, const Ts* a) {
    for (size_t k = 0; k < C->N; k++)
        if (!ts_zero(C, a, k) && !ts_resolved(C, a, k)) return false;
    return C->ok;
}

/* a^k for an integer k >= 2 by binary powering: convolutions only, so a
 * series of one sign has no cancellation at all.  The coefficients ahead of
 * the valuation are zeros by ts_zero and are cleared first, as the exact
 * expansion has them. */
static Ts* ts_pow_binary(TsCtx* C, const Ts* a, long k) {
    size_t v = ts_valuation(C, a);
    Ts* base = ts_new(C);
    if (!base) return NULL;
    for (size_t i = v; i < C->N; i++) {
        mpfr_set(base->c[i], a->c[i], MPFR_RNDN);
        base->m[i] = a->m[i];
    }
    Ts* acc = NULL;
    while (k > 0 && C->ok) {
        if (k & 1) {
            if (!acc) {
                acc = ts_new(C);
                if (!acc) break;
                for (size_t i = 0; i < C->N; i++) {
                    mpfr_set(acc->c[i], base->c[i], MPFR_RNDN);
                    acc->m[i] = base->m[i];
                }
            } else {
                Ts* p = ts_mul(C, acc, base);
                ts_free(C, acc);
                acc = p;
                if (!acc) break;
            }
        }
        k >>= 1;
        if (k > 0) {
            Ts* sq = ts_mul(C, base, base);
            ts_free(C, base);
            base = sq;
            if (!base) break;
        }
    }
    ts_free(C, base);
    if (!C->ok) { ts_free(C, acc); return NULL; }
    return acc;
}

/* a^r: Miller's recurrence, O(N^2) for any exponent.  Its weights
 * (r + 1) j - m change sign, so for a small integer exponent a high-order
 * coefficient can cancel below the working precision (the square of a
 * positive series at order 100); when it leaves any coefficient unresolved, an
 * integer power is recomputed by convolution instead. */
static Ts* ts_pow(TsCtx* C, const Ts* a, bool is_int, long k, const mpfr_t r) {
    Ts* res = ts_pow_miller(C, a, is_int, k, r);
    if (!res || !is_int || k < 2 || ts_all_settled(C, res)) return res;
    if (!C->ok) { ts_free(C, res); return NULL; }
    ts_free(C, res);
    return ts_pow_binary(C, a, k);
}

/* Exp[a]:  e' = a' e,  e_k = (1/k) sum_{j=1..k} j a_j e_{k-j}. */
static Ts* ts_exp(TsCtx* C, const Ts* a) {
    Ts* e = ts_new(C);
    if (!e) return NULL;
    mpfr_t t;
    mpfr_init2(t, C->wp);
    mpfr_exp(e->c[0], a->c[0], MPFR_RNDN);
    e->m[0] = ts_absd(e->c[0]) * (1.0 + a->m[0]);
    for (size_t k = 1; k < C->N; k++) {
        double ms = 0.0;
        for (size_t j = 1; j <= k; j++) {
            if (!mpfr_zero_p(a->c[j]) && !mpfr_zero_p(e->c[k - j])) {
                mpfr_mul(t, a->c[j], e->c[k - j], MPFR_RNDN);
                mpfr_mul_ui(t, t, (unsigned long)j, MPFR_RNDN);
                mpfr_add(e->c[k], e->c[k], t, MPFR_RNDN);
            }
            ms += (double)j * a->m[j] * e->m[k - j];
        }
        mpfr_div_ui(e->c[k], e->c[k], (unsigned long)k, MPFR_RNDN);
        e->m[k] = ms / (double)k;
    }
    mpfr_clear(t);
    return e;
}

/* Log[a], a_0 > 0:  a l' = a',  l_k = (a_k - (1/k) sum_{j=1..k-1} j l_j a_{k-j}) / a_0. */
static Ts* ts_log(TsCtx* C, const Ts* a) {
    if (ts_zero(C, a, 0) || !ts_resolved(C, a, 0) || mpfr_sgn(a->c[0]) <= 0) {
        C->ok = false; return NULL;
    }
    Ts* l = ts_new(C);
    if (!l) return NULL;
    mpfr_t t;
    mpfr_init2(t, C->wp);
    double absa0 = ts_absd(a->c[0]);
    double ra = a->m[0] / absa0;
    mpfr_log(l->c[0], a->c[0], MPFR_RNDN);
    l->m[0] = ts_absd(l->c[0]) + ra;
    for (size_t k = 1; k < C->N; k++) {
        double ms = 0.0;
        for (size_t j = 1; j < k; j++) {
            if (!mpfr_zero_p(l->c[j]) && !mpfr_zero_p(a->c[k - j])) {
                mpfr_mul(t, l->c[j], a->c[k - j], MPFR_RNDN);
                mpfr_mul_ui(t, t, (unsigned long)j, MPFR_RNDN);
                mpfr_add(l->c[k], l->c[k], t, MPFR_RNDN);
            }
            ms += (double)j * l->m[j] * a->m[k - j];
        }
        mpfr_div_ui(l->c[k], l->c[k], (unsigned long)k, MPFR_RNDN);
        mpfr_sub(l->c[k], a->c[k], l->c[k], MPFR_RNDN);
        mpfr_div(l->c[k], l->c[k], a->c[0], MPFR_RNDN);
        l->m[k] = (a->m[k] + ms / (double)k) / absa0 + ts_absd(l->c[k]) * ra;
    }
    mpfr_clear(t);
    return l;
}

/* Sin/Cos (hyp = false) or Sinh/Cosh (hyp = true) of a, both computed by
 *   s' = a' c,  c' = -/+ a' s;  returns the sine-like one when want_sin. */
static Ts* ts_sincos(TsCtx* C, const Ts* a, bool hyp, bool want_sin) {
    Ts* s = ts_new(C);
    Ts* c = ts_new(C);
    if (!s || !c) { ts_free(C, s); ts_free(C, c); return NULL; }
    mpfr_t t;
    mpfr_init2(t, C->wp);
    if (hyp) { mpfr_sinh(s->c[0], a->c[0], MPFR_RNDN); mpfr_cosh(c->c[0], a->c[0], MPFR_RNDN); }
    else     { mpfr_sin(s->c[0], a->c[0], MPFR_RNDN);  mpfr_cos(c->c[0], a->c[0], MPFR_RNDN); }
    double s0 = ts_absd(s->c[0]), c0 = ts_absd(c->c[0]);
    s->m[0] = s0 + c0 * a->m[0];
    c->m[0] = c0 + s0 * a->m[0];
    for (size_t k = 1; k < C->N; k++) {
        double mss = 0.0, msc = 0.0;
        for (size_t j = 1; j <= k; j++) {
            if (!mpfr_zero_p(a->c[j])) {
                if (!mpfr_zero_p(c->c[k - j])) {
                    mpfr_mul(t, a->c[j], c->c[k - j], MPFR_RNDN);
                    mpfr_mul_ui(t, t, (unsigned long)j, MPFR_RNDN);
                    mpfr_add(s->c[k], s->c[k], t, MPFR_RNDN);
                }
                if (!mpfr_zero_p(s->c[k - j])) {
                    mpfr_mul(t, a->c[j], s->c[k - j], MPFR_RNDN);
                    mpfr_mul_ui(t, t, (unsigned long)j, MPFR_RNDN);
                    mpfr_add(c->c[k], c->c[k], t, MPFR_RNDN);
                }
            }
            mss += (double)j * a->m[j] * c->m[k - j];
            msc += (double)j * a->m[j] * s->m[k - j];
        }
        mpfr_div_ui(s->c[k], s->c[k], (unsigned long)k, MPFR_RNDN);
        mpfr_div_ui(c->c[k], c->c[k], (unsigned long)k, MPFR_RNDN);
        if (!hyp) mpfr_neg(c->c[k], c->c[k], MPFR_RNDN);
        s->m[k] = mss / (double)k;
        c->m[k] = msc / (double)k;
    }
    mpfr_clear(t);
    if (want_sin) { ts_free(C, c); return s; }
    ts_free(C, s);
    return c;
}

/* ---- expression walk ------------------------------------------------------ */

static Ts* ts_build(TsCtx* C, const Expr* e);

/* x-free subexpression as a constant series. */
static Ts* ts_build_const(TsCtx* C, const Expr* e) {
    mpfr_t v;
    mpfr_init2(v, C->wp);
    Ts* r = NULL;
    if (ts_numeric_value(C, e, v)) r = ts_const(C, v);
    else C->ok = false;
    mpfr_clear(v);
    return r;
}

/* Exp[s Log[c]] for a constant c > 0 (c^s, the exponent depending on x). */
static Ts* ts_const_pow(TsCtx* C, const Expr* base, const Expr* ex) {
    mpfr_t c;
    mpfr_init2(c, C->wp);
    Ts* r = NULL;
    if (ts_numeric_value(C, base, c) && mpfr_sgn(c) > 0) {
        mpfr_log(c, c, MPFR_RNDN);
        Ts* s = ts_build(C, ex);
        Ts* sl = s ? ts_scale(C, s, c) : NULL;
        if (sl) r = ts_exp(C, sl);
        ts_free(C, s); ts_free(C, sl);
    } else C->ok = false;
    mpfr_clear(c);
    return r;
}

static Ts* ts_build_power(TsCtx* C, const Expr* e) {
    const Expr* base = e->data.function.args[0];
    const Expr* ex   = e->data.function.args[1];
    bool bfree = ts_free_of(base, C->x), efree = ts_free_of(ex, C->x);
    if (bfree) return ts_const_pow(C, base, ex);
    if (efree) {
        Ts* b = ts_build(C, base);
        if (!b) return NULL;
        Ts* r = NULL;
        mpfr_t rv;
        mpfr_init2(rv, C->wp);
        if (ex->type == EXPR_INTEGER) {
            r = ts_pow(C, b, true, (long)ex->data.integer, rv);
        } else if (ts_numeric_value(C, ex, rv)) {
            if (mpfr_integer_p(rv) && mpfr_fits_slong_p(rv, MPFR_RNDN))
                r = ts_pow(C, b, true, mpfr_get_si(rv, MPFR_RNDN), rv);
            else
                r = ts_pow(C, b, false, 0, rv);
        } else C->ok = false;
        mpfr_clear(rv);
        ts_free(C, b);
        return r;
    }
    /* base and exponent both depend on x: Exp[ex Log[base]]. */
    Ts* b = ts_build(C, base);
    Ts* lb = b ? ts_log(C, b) : NULL;
    Ts* s = lb ? ts_build(C, ex) : NULL;
    Ts* p = s ? ts_mul(C, s, lb) : NULL;
    Ts* r = p ? ts_exp(C, p) : NULL;
    ts_free(C, b); ts_free(C, lb); ts_free(C, s); ts_free(C, p);
    return r;
}

static Ts* ts_build(TsCtx* C, const Expr* e) {
    if (!C->ok) return NULL;
    if (ts_free_of(e, C->x)) return ts_build_const(C, e);
    if (e->type == EXPR_SYMBOL) return ts_var(C);           /* e is x itself */
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) {
        C->ok = false; return NULL;
    }
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;

    if (h == SYM_Plus && n >= 1) {
        Ts* acc = ts_build(C, e->data.function.args[0]);
        for (size_t i = 1; acc && i < n; i++) {
            Ts* t = ts_build(C, e->data.function.args[i]);
            Ts* s = t ? ts_add(C, acc, t) : NULL;
            ts_free(C, acc); ts_free(C, t);
            acc = s;
        }
        return acc;
    }
    if (h == SYM_Times && n >= 1) {
        /* Constant factors are multiplied into one scale; the rest convolve. */
        mpfr_t k, v;
        mpfr_inits2(C->wp, k, v, (mpfr_ptr)0);
        mpfr_set_si(k, 1, MPFR_RNDN);
        bool have_k = false;
        Ts* acc = NULL;
        for (size_t i = 0; C->ok && i < n; i++) {
            const Expr* f = e->data.function.args[i];
            if (ts_free_of(f, C->x)) {
                if (ts_numeric_value(C, f, v)) { mpfr_mul(k, k, v, MPFR_RNDN); have_k = true; }
                else C->ok = false;
                continue;
            }
            Ts* t = ts_build(C, f);
            if (!t) break;
            if (!acc) { acc = t; continue; }
            Ts* p = ts_mul(C, acc, t);
            ts_free(C, acc); ts_free(C, t);
            acc = p;
            if (!acc) break;
        }
        if (acc && C->ok && have_k) {
            Ts* s = ts_scale(C, acc, k);
            ts_free(C, acc);
            acc = s;
        }
        mpfr_clears(k, v, (mpfr_ptr)0);
        if (!C->ok) { ts_free(C, acc); return NULL; }
        return acc;
    }
    if (h == SYM_Power && n == 2) return ts_build_power(C, e);
    if (n == 1 && (h == SYM_Exp || h == SYM_Log || h == SYM_Sin || h == SYM_Cos ||
                   h == SYM_Sinh || h == SYM_Cosh)) {
        Ts* a = ts_build(C, e->data.function.args[0]);
        if (!a) return NULL;
        Ts* r;
        if (h == SYM_Exp)       r = ts_exp(C, a);
        else if (h == SYM_Log)  r = ts_log(C, a);
        else if (h == SYM_Sin)  r = ts_sincos(C, a, false, true);
        else if (h == SYM_Cos)  r = ts_sincos(C, a, false, false);
        else if (h == SYM_Sinh) r = ts_sincos(C, a, true, true);
        else                    r = ts_sincos(C, a, true, false);
        ts_free(C, a);
        return r;
    }
    C->ok = false;
    return NULL;
}

/* ---- entry ---------------------------------------------------------------- */

/* The round trip's precision, as internal_rationalize_then_numericalize
 * chooses it: the smallest MPFR precision over the inexact leaves when every
 * one is arbitrary-precision, -1 when any is a machine number. */
static long ts_input_bits(const Expr* e) {
    if (e->type == EXPR_REAL) return -1;
    if (e->type == EXPR_MPFR) return (long)mpfr_get_prec(e->data.mpfr);
    if (e->type != EXPR_FUNCTION) return 0;
    long best = ts_input_bits(e->data.function.head);
    if (best < 0) return -1;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        long b = ts_input_bits(e->data.function.args[i]);
        if (b < 0) return -1;
        if (b > 0 && (best == 0 || b < best)) best = b;
    }
    return best;
}

static bool ts_mentions_head(const Expr* e, const char* sym) {
    if (e->type != EXPR_FUNCTION) return false;
    if (ts_head_is(e, sym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ts_mentions_head(e->data.function.args[i], sym)) return true;
    return false;
}

Expr* series_inexact_taylor(const Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    const Expr* f = res->data.function.args[0];
    const Expr* spec = res->data.function.args[1];
    if (!ts_head_is(spec, SYM_List) || spec->data.function.arg_count != 3) return NULL;
    const Expr* x = spec->data.function.args[0];
    const Expr* x0 = spec->data.function.args[1];
    const Expr* ne = spec->data.function.args[2];
    if (x->type != EXPR_SYMBOL || ne->type != EXPR_INTEGER) return NULL;
    if (ne->data.integer < 0 || ne->data.integer >= TS_MAX_TERMS) return NULL;
    if (ts_head_is(f, SYM_List) || ts_mentions_head(f, SYM_SeriesData)) return NULL;
    if (ts_free_of(f, x)) return NULL;

    long b1 = ts_input_bits(f), b2 = ts_input_bits(x0);
    long bits = (b1 < 0 || b2 < 0) ? -1 : (b1 == 0 ? b2 : (b2 == 0 ? b1 : (b1 < b2 ? b1 : b2)));
    if (bits <= 53) bits = 0;                        /* machine round trip */

    Expr* fr = internal_force_rationalize_bits(f, bits);
    Expr* x0r = internal_force_rationalize_bits(x0, bits);
    if (!fr || !x0r || internal_contains_inexact(fr) || internal_contains_inexact(x0r)) {
        if (fr) expr_free(fr);
        if (x0r) expr_free(x0r);
        return NULL;
    }
    Expr* fe = eval_and_free(fr);
    Expr* x0e = eval_and_free(x0r);

    TsCtx C;
    C.N = (size_t)ne->data.integer + 1;
    C.tb = bits > 0 ? (mpfr_prec_t)bits : 53;
    long lg = 0;
    while (((size_t)1 << lg) < C.N + 1) lg++;
    C.slack = 16 + 2 * lg;
    C.wp = C.tb + 64 + (mpfr_prec_t)C.slack;
    C.x = x;
    C.ok = true;
    mpfr_init2(C.x0, C.wp);

    Expr* out = NULL;
    Ts* s = NULL;
    if (!ts_free_of(fe, x) && ts_numeric_value(&C, x0e, C.x0))
        s = ts_build(&C, fe);

    if (s && C.ok) {
        NumericSpec ospec = numeric_machine_spec();
        if (bits > 0) { ospec.mode = NUMERIC_MODE_MPFR; ospec.bits = bits; }
        Expr** coefs = calloc(C.N, sizeof(Expr*));
        Expr* zero0 = expr_new_integer(0);
        Expr* zero = numericalize(zero0, ospec);
        expr_free(zero0);
        bool good = coefs && zero;
        for (size_t k = 0; good && k < C.N; k++) {
            /* MPFR: an exact 0 and the full working-precision value, for the
             * whole-series numericalisation below. */
            if (ts_zero(&C, s, k)) {
                coefs[k] = bits > 0 ? expr_new_integer(0) : expr_copy(zero);
                continue;
            }
            if (!C.ok || !ts_resolved(&C, s, k)) { good = false; break; }
            if (bits > 0) {
                coefs[k] = expr_new_mpfr_copy(s->c[k]);
            } else {
                double d = mpfr_get_d(s->c[k], MPFR_RNDN);
                if (!isfinite(d) || fabs(d) < DBL_MIN) { good = false; break; }
                coefs[k] = expr_new_real(d);
            }
        }
        if (good && C.ok && bits > 0) {
            /* At arbitrary precision the zeros must be the round trip's zeros,
             * not a bare N[0]: numericalising the whole exact SeriesData (a
             * call) raises the working precision by the guard bits, so its
             * exact zeros come back as 0 with accuracy bits + guard, and the
             * nonzero values are rounded to `bits`.  A bare N[0, bits] has
             * accuracy `bits`, and downstream significance arithmetic pays for
             * the difference: the 110 x 110 Schroeder solve of hd20.9 kept 49
             * digits of 60 instead of 59.  So numericalise the same shape: the
             * exact zeros, the working-precision values, and the exact point. */
            Expr* list = expr_new_function(expr_new_symbol(SYM_List), coefs, C.N);
            Expr* args[6] = { expr_copy((Expr*)x), expr_copy(x0e), list,
                              expr_new_integer(0), expr_new_integer((int64_t)C.N),
                              expr_new_integer(1) };
            Expr* exactish = expr_new_function(expr_new_symbol(SYM_SeriesData), args, 6);
            out = numericalize(exactish, ospec);
            expr_free(exactish);
            if (out && ts_head_is(out, SYM_SeriesData) && out->data.function.arg_count == 6) {
                if (!internal_contains_inexact(x0)) {       /* the exact point back */
                    Expr* a2[6];
                    for (size_t i = 0; i < 6; i++)
                        a2[i] = (i == 1) ? expr_copy(x0e) : expr_copy(out->data.function.args[i]);
                    expr_free(out);
                    out = expr_new_function(expr_new_symbol(SYM_SeriesData), a2, 6);
                }
            } else if (out) {
                expr_free(out);
                out = NULL;
            }
        } else if (good && C.ok) {
            Expr* list = expr_new_function(expr_new_symbol(SYM_List), coefs, C.N);
            Expr* x0out = internal_contains_inexact(x0) ? numericalize(x0e, ospec)
                                                        : expr_copy(x0e);
            Expr* args[6] = { expr_copy((Expr*)x), x0out, list,
                              expr_new_integer(0), expr_new_integer((int64_t)C.N),
                              expr_new_integer(1) };
            out = expr_new_function(expr_new_symbol(SYM_SeriesData), args, 6);
        } else if (coefs) {
            for (size_t k = 0; k < C.N; k++) if (coefs[k]) expr_free(coefs[k]);
        }
        free(coefs);
        if (zero) expr_free(zero);
    }
    ts_free(&C, s);
    mpfr_clear(C.x0);
    expr_free(fe);
    expr_free(x0e);
    return out;
}

#else  /* !USE_MPFR */

Expr* series_inexact_taylor(const Expr* res) { (void)res; return NULL; }

#endif
