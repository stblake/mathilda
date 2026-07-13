/*
 * contfrac.c — ContinuedFraction[x] and ContinuedFraction[x, n].
 *
 * Computes the *simple* continued-fraction expansion of x:
 *
 *     {a1, a2, a3, ...}  <->  a1 + 1/(a2 + 1/(a3 + ...))
 *
 * Four input regimes are handled, each with the appropriate exact or
 * precision-aware algorithm:
 *
 *  1. Exact rationals (Integer, BigInt, Rational[p, q]).
 *     The Euclidean algorithm with floor division produces the canonical
 *     terminating expansion (last term >= 2, never the {..., k-1, 1} form).
 *     A finite rational yields finitely many terms, so ContinuedFraction[x, n]
 *     may return fewer than n terms.
 *
 *  2. Quadratic irrationals of the form Sqrt[D], D a positive non-square
 *     integer.  The classic periodic surd recurrence
 *         m_{i+1} = a_i q_i - m_i,
 *         q_{i+1} = (D - m_{i+1}^2) / q_i,
 *         a_{i+1} = floor((a0 + m_{i+1}) / q_{i+1})
 *     is purely periodic after a0 for an integer square root, so the period
 *     is detected when the state (m, q) first repeats.  Without n the result
 *     is {a0, {b1, ..., bk}} (the bracketed list is the repeating block);
 *     with n the periodic sequence is unrolled to exactly n terms.
 *     General quadratic irrationals (e.g. (1 + Sqrt[5])/2, Sqrt[2/3]) are not
 *     recognised symbolically here — supply an explicit n to use the numeric
 *     path instead.
 *
 *  3. Inexact reals (machine Real, or arbitrary-precision MPFR).  Terms are
 *     extracted by repeated reciprocation while tracking the absolute
 *     uncertainty of the running value (initially |x| * 2^-prec).  Expansion
 *     stops as soon as the integer part of the value is no longer determined
 *     by the available precision — i.e. ContinuedFraction stops when it runs
 *     out of precision.  The iteration carries 64 guard bits so arithmetic
 *     round-off stays far below the modelled input uncertainty.
 *
 *  4. Exact symbolic reals with an explicit n (Pi, Sqrt[E], Exp[Pi Sqrt[163]],
 *     ...).  The value is numericised via N[x, digits] at successively
 *     doubled precision until n reliable terms are obtained.
 *
 * Anything else (an exact non-rational, non-quadratic with no n, or a
 * non-real numeric value) is left unevaluated (the builtin returns NULL).
 *
 * Attributes: Listable, Protected.
 */

#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "attr.h"
#include "sym_names.h"
#include "contfrac.h"

#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>   /* DBL_MANT_DIG */

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ */
/* Growable vector of GMP integers (one continued-fraction term each). */
/* ------------------------------------------------------------------ */
typedef struct {
    mpz_t* v;
    size_t n;
    size_t cap;
} TermVec;

static void tv_init(TermVec* t) {
    t->v = NULL;
    t->n = 0;
    t->cap = 0;
}

static void tv_push(TermVec* t, const mpz_t z) {
    if (t->n == t->cap) {
        size_t ncap = t->cap ? t->cap * 2 : 16;
        mpz_t* nv = realloc(t->v, ncap * sizeof(mpz_t));
        if (!nv) return;            /* OOM: silently stop growing */
        t->v = nv;
        t->cap = ncap;
    }
    mpz_init_set(t->v[t->n], z);
    t->n++;
}

static void tv_clear(TermVec* t) {
    for (size_t i = 0; i < t->n; i++) mpz_clear(t->v[i]);
    free(t->v);
    t->v = NULL;
    t->n = t->cap = 0;
}

/* mpz -> Expr, demoting to a machine integer when it fits. */
static Expr* mpz_to_expr(const mpz_t z) {
    Expr* e = expr_new_bigint_from_mpz(z);
    return expr_bigint_normalize(e);
}

/* Build a flat List[...] expression from the accumulated terms. */
static Expr* tv_to_list(const TermVec* t) {
    size_t n = t->n;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) args[i] = mpz_to_expr(t->v[i]);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), args, n);
    free(args);
    return list;
}

/* ------------------------------------------------------------------ */
/* Input recognisers.                                                  */
/* ------------------------------------------------------------------ */

/* If e is an exact rational (Integer / BigInt / Rational[int, int]), set
 * p, q (q > 0 not enforced here) and return true.  On true return both p
 * and q are mpz_init'd and the caller must mpz_clear them; on false return
 * neither is initialised. */
static bool expr_get_rational(const Expr* e, mpz_t p, mpz_t q) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {
        expr_to_mpz(e, p);          /* initialises p */
        mpz_init_set_ui(q, 1);
        return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2) {
        Expr* a = e->data.function.args[0];
        Expr* b = e->data.function.args[1];
        if (expr_is_integer_like(a) && expr_is_integer_like(b)) {
            expr_to_mpz(a, p);
            expr_to_mpz(b, q);
            return true;
        }
    }
    return false;
}

/* If e is Sqrt[D] = Power[D, Rational[1, 2]] with D a positive integer-like
 * value, set D and return true (D is mpz_init'd, caller clears).  Otherwise
 * return false with D untouched. */
static bool expr_is_sqrt_int(const Expr* e, mpz_t D) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Power ||
        e->data.function.arg_count != 2)
        return false;
    Expr* base = e->data.function.args[0];
    Expr* ex   = e->data.function.args[1];
    if (!expr_is_integer_like(base)) return false;
    /* exponent must be exactly Rational[1, 2] */
    if (!(ex->type == EXPR_FUNCTION && ex->data.function.head &&
          ex->data.function.head->type == EXPR_SYMBOL &&
          ex->data.function.head->data.symbol.name == SYM_Rational &&
          ex->data.function.arg_count == 2))
        return false;
    Expr* en = ex->data.function.args[0];
    Expr* ed = ex->data.function.args[1];
    if (!(en->type == EXPR_INTEGER && en->data.integer == 1 &&
          ed->type == EXPR_INTEGER && ed->data.integer == 2))
        return false;
    expr_to_mpz(base, D);
    if (mpz_sgn(D) <= 0) { mpz_clear(D); return false; }
    return true;
}

/* ------------------------------------------------------------------ */
/* 1. Exact rational expansion (Euclidean algorithm).                  */
/* ------------------------------------------------------------------ */
static void cf_rational(mpz_t p, mpz_t q, long nlimit, TermVec* out) {
    if (mpz_sgn(q) < 0) { mpz_neg(p, p); mpz_neg(q, q); }   /* keep q > 0 */
    mpz_t a, r;
    mpz_inits(a, r, (mpz_ptr)0);
    long count = 0;
    while (mpz_sgn(q) != 0) {
        if (nlimit > 0 && count >= nlimit) break;
        mpz_fdiv_qr(a, r, p, q);     /* a = floor(p/q), r = p - a*q in [0, q) */
        tv_push(out, a);
        count++;
        mpz_set(p, q);
        mpz_set(q, r);
    }
    mpz_clears(a, r, (mpz_ptr)0);
}

/* ------------------------------------------------------------------ */
/* 2. Quadratic surd Sqrt[D] (D positive integer).                     */
/* ------------------------------------------------------------------ */

/* Upper bound on the period length we are willing to materialise for the
 * no-count form.  The period of Sqrt[D] can be O(sqrt(D) log D), which is
 * astronomically large for big D; past this bound we decline to evaluate. */
#define CF_PERIOD_CAP 2000000L

typedef enum { CFSQ_OK, CFSQ_PERFECT_SQUARE, CFSQ_TOO_LONG } CfSqStatus;

/* Directly generate the first `count` surd terms a0, a1, ... into out (flat).
 * No period detection — bounded by `count`, so it is cheap even for huge D. */
static void cf_sqrt_terms(const mpz_t D, long count, TermVec* out) {
    if (count <= 0) return;
    if (mpz_perfect_square_p(D)) {   /* defensive: really an integer */
        mpz_t s;
        mpz_init(s);
        mpz_sqrt(s, D);
        tv_push(out, s);
        mpz_clear(s);
        return;
    }
    mpz_t a0, m, q, a, mnext, qnext, tmp;
    mpz_inits(a0, m, q, a, mnext, qnext, tmp, (mpz_ptr)0);
    mpz_sqrt(a0, D);
    mpz_set_ui(m, 0);
    mpz_set_ui(q, 1);
    mpz_set(a, a0);
    tv_push(out, a0);                /* a_0 */
    for (long i = 1; i < count; i++) {
        mpz_mul(mnext, a, q);
        mpz_sub(mnext, mnext, m);
        mpz_mul(tmp, mnext, mnext);
        mpz_sub(tmp, D, tmp);
        mpz_divexact(qnext, tmp, q);
        mpz_add(tmp, a0, mnext);
        mpz_fdiv_q(a, tmp, qnext);
        mpz_set(m, mnext);
        mpz_set(q, qnext);
        tv_push(out, a);
    }
    mpz_clears(a0, m, q, a, mnext, qnext, tmp, (mpz_ptr)0);
}

/* Compute a0 = floor(sqrt(D)) and the repeating block.  The expansion of an
 * integer square root is purely periodic after a0, so the period closes when
 * the state (m, q) first repeats.  a0_out is mpz_init'd by the caller; period
 * is an initialised empty TermVec. */
static CfSqStatus cf_sqrt_period(const mpz_t D, mpz_t a0_out, TermVec* period) {
    if (mpz_perfect_square_p(D)) return CFSQ_PERFECT_SQUARE;

    mpz_t a0, m, q, m1, q1, a, mnext, qnext, tmp;
    mpz_inits(a0, m, q, m1, q1, a, mnext, qnext, tmp, (mpz_ptr)0);

    mpz_sqrt(a0, D);                 /* floor(sqrt(D)) */
    mpz_set(a0_out, a0);
    mpz_set_ui(m, 0);                /* m_0 */
    mpz_set_ui(q, 1);                /* q_0 */
    mpz_set(a, a0);                  /* a_0 */

    CfSqStatus st = CFSQ_TOO_LONG;
    int first = 1;
    for (long step = 0; step < CF_PERIOD_CAP; step++) {
        mpz_mul(mnext, a, q);        /* m_{i+1} = a_i q_i - m_i */
        mpz_sub(mnext, mnext, m);
        mpz_mul(tmp, mnext, mnext);  /* q_{i+1} = (D - m_{i+1}^2) / q_i */
        mpz_sub(tmp, D, tmp);
        mpz_divexact(qnext, tmp, q);
        mpz_add(tmp, a0, mnext);     /* a_{i+1} = floor((a0 + m_{i+1}) / q_{i+1}) */
        mpz_fdiv_q(a, tmp, qnext);
        mpz_set(m, mnext);
        mpz_set(q, qnext);

        if (first) {
            mpz_set(m1, m);
            mpz_set(q1, q);
            first = 0;
        } else if (mpz_cmp(m, m1) == 0 && mpz_cmp(q, q1) == 0) {
            st = CFSQ_OK;            /* state repeated -> period complete */
            break;
        }
        tv_push(period, a);
    }

    mpz_clears(a0, m, q, m1, q1, a, mnext, qnext, tmp, (mpz_ptr)0);
    return st;
}

static Expr* cf_sqrt_build(const mpz_t D, long nlimit) {
    /* With an explicit count, generate exactly that many terms directly. */
    if (nlimit > 0) {
        TermVec out;
        tv_init(&out);
        cf_sqrt_terms(D, nlimit, &out);
        Expr* r = tv_to_list(&out);
        tv_clear(&out);
        return r;
    }

    /* No count: emit {a0, {period...}}. */
    mpz_t a0;
    mpz_init(a0);
    TermVec period;
    tv_init(&period);
    CfSqStatus st = cf_sqrt_period(D, a0, &period);

    if (st == CFSQ_PERFECT_SQUARE) {
        mpz_sqrt(a0, D);
        TermVec one;
        tv_init(&one);
        tv_push(&one, a0);
        Expr* r = tv_to_list(&one);
        tv_clear(&one);
        tv_clear(&period);
        mpz_clear(a0);
        return r;
    }
    if (st == CFSQ_TOO_LONG) {        /* period too long to materialise */
        tv_clear(&period);
        mpz_clear(a0);
        return NULL;
    }

    Expr* inner = tv_to_list(&period);
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = mpz_to_expr(a0);
    args[1] = inner;
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), args, 2);
    free(args);

    tv_clear(&period);
    mpz_clear(a0);
    return result;
}

/* ------------------------------------------------------------------ */
/* 3. Inexact / numeric expansion.                                     */
/* ------------------------------------------------------------------ */
#ifdef USE_MPFR
/* Extract continued-fraction terms of x, whose absolute uncertainty is
 * |x| * 2^-pbits.  Stops at nlimit terms (nlimit <= 0 == unlimited) or when
 * the next integer part is no longer determined by the precision. */
static void cf_float_mpfr(const mpfr_t x, long pbits, long nlimit, TermVec* out) {
    if (mpfr_nan_p(x) || mpfr_inf_p(x)) return;

    mpfr_prec_t wp = (mpfr_prec_t)pbits + 64;   /* guard bits */
    mpfr_t v, frac, vlo, vhi, tmp, err;
    mpfr_inits2(wp, v, frac, vlo, vhi, tmp, (mpfr_ptr)0);
    mpfr_init2(err, 96);

    mpfr_set(v, x, MPFR_RNDN);
    mpfr_abs(tmp, x, MPFR_RNDN);
    mpfr_mul_2si(err, tmp, -pbits, MPFR_RNDU);  /* |x| * 2^-pbits */

    mpz_t a, alo, ahi;
    mpz_inits(a, alo, ahi, (mpz_ptr)0);

    /* Round-off scale at the working precision: a value this close to an
     * integer is integral up to arithmetic noise (e.g. 0.75 -> 1/(4/3-1)
     * lands at 3 + O(2^-wp)).  This is far tighter than `err`, so it only
     * fires on genuinely terminating quotients and never adds a spurious
     * term where precision is merely exhausted. */
    long count = 0;
    for (long i = 0; i < 1000000L; i++) {
        if (nlimit > 0 && count >= nlimit) break;
        /* Integral (within round-off) -> the expansion terminates here. */
        mpfr_get_z(a, v, MPFR_RNDN);            /* nearest integer */
        mpfr_sub_z(frac, v, a, MPFR_RNDN);      /* signed distance to it */
        mpfr_abs(frac, frac, MPFR_RNDN);
        mpfr_abs(tmp, v, MPFR_RNDN);
        mpfr_mul_2si(tmp, tmp, -(long)(wp - 4), MPFR_RNDU);  /* round-off scale */
        if (mpfr_cmp(frac, tmp) <= 0) {
            tv_push(out, a);
            count++;
            break;
        }
        /* Is floor(v) determined?  [v-err, v+err] must not straddle an int. */
        mpfr_sub(vlo, v, err, MPFR_RNDN);
        mpfr_add(vhi, v, err, MPFR_RNDN);
        mpfr_get_z(alo, vlo, MPFR_RNDD);
        mpfr_get_z(ahi, vhi, MPFR_RNDD);
        if (mpz_cmp(alo, ahi) != 0) break;      /* precision exhausted */

        mpfr_get_z(a, v, MPFR_RNDD);            /* a = floor(v) */
        tv_push(out, a);
        count++;

        mpfr_sub_z(frac, v, a, MPFR_RNDN);      /* frac in [0, 1) */
        if (mpfr_cmp(frac, err) <= 0) break;    /* terminates (within noise) */

        mpfr_ui_div(v, 1, frac, MPFR_RNDN);     /* v <- 1/frac */
        mpfr_div(err, err, frac, MPFR_RNDU);    /* err <- err / frac^2 */
        mpfr_div(err, err, frac, MPFR_RNDU);
    }

    mpz_clears(a, alo, ahi, (mpz_ptr)0);
    mpfr_clears(v, frac, vlo, vhi, tmp, (mpfr_ptr)0);
    mpfr_clear(err);
}

/* Run the float expansion on an Expr that is a Real or MPFR value. */
static void cf_float_expr(const Expr* x, long nlimit, TermVec* out) {
    if (x->type == EXPR_MPFR) {
        cf_float_mpfr(x->data.mpfr, (long)mpfr_get_prec(x->data.mpfr), nlimit, out);
    } else { /* EXPR_REAL */
        mpfr_t v;
        mpfr_init2(v, DBL_MANT_DIG);
        mpfr_set_d(v, x->data.real, MPFR_RNDN);
        cf_float_mpfr(v, DBL_MANT_DIG, nlimit, out);
        mpfr_clear(v);
    }
}

#else  /* !USE_MPFR : machine-double fallback (degraded precision) */

static void cf_float_double(double x, long nlimit, TermVec* out) {
    if (!isfinite(x)) return;
    double v = x;
    double err = fabs(x) * ldexp(1.0, -(DBL_MANT_DIG - 1));
    mpz_t a;
    mpz_init(a);
    long count = 0;
    for (long i = 0; i < 100000L; i++) {
        if (nlimit > 0 && count >= nlimit) break;
        double k = floor(v + 0.5);        /* nearest integer */
        if (fabs(v - k) <= fabs(v) * 8.0 * DBL_EPSILON) {  /* integral -> terminates */
            mpz_set_d(a, k);
            tv_push(out, a);
            count++;
            break;
        }
        if (floor(v - err) != floor(v + err)) break;
        double af = floor(v);
        mpz_set_d(a, af);
        tv_push(out, a);
        count++;
        double frac = v - af;
        if (frac <= err) break;
        v = 1.0 / frac;
        err = err / (frac * frac);
    }
    mpz_clear(a);
}

static void cf_float_expr(const Expr* x, long nlimit, TermVec* out) {
    /* Without MPFR only machine Real values are inexact. */
    cf_float_double(x->data.real, nlimit, out);
}
#endif

static Expr* cf_inexact(const Expr* x, long nlimit) {
    TermVec out;
    tv_init(&out);
    cf_float_expr(x, nlimit, &out);
    Expr* r = tv_to_list(&out);
    tv_clear(&out);
    return r;
}

#ifdef USE_MPFR
/* Do TermVecs a and b agree on their first n terms? */
static bool tv_prefix_eq(const TermVec* a, const TermVec* b, long n) {
    if ((long)a->n < n || (long)b->n < n) return false;
    for (long i = 0; i < n; i++)
        if (mpz_cmp(a->v[i], b->v[i]) != 0) return false;
    return true;
}
#endif

/* 4. Adaptive numeric expansion for an exact symbolic real with explicit n.
 *
 * The precision tag of N[x, digits] can overstate the value's true accuracy
 * (e.g. exp amplifies the error of its argument), so the term nearest the
 * precision boundary cannot be trusted on a single evaluation.  Instead we
 * raise the precision until two consecutive evaluations agree on all n
 * terms — convergence that no single-precision boundary artefact survives. */
static Expr* cf_exact_numeric(const Expr* x, long nlimit) {
#ifdef USE_MPFR
    long digits = nlimit > 24 ? nlimit : 24;
    TermVec prev;
    tv_init(&prev);
    bool have_prev = false;

    for (int tries = 0; tries < 48 && digits <= 400000; tries++, digits *= 2) {
        Expr* nexpr = expr_new_function(
            expr_new_symbol(SYM_N),
            (Expr*[]){ expr_copy((Expr*)x), expr_new_integer(digits) }, 2);
        Expr* val = eval_and_free(nexpr);
        if (!val) { tv_clear(&prev); return NULL; }

        if (val->type == EXPR_MPFR || val->type == EXPR_REAL) {
            if (val->type == EXPR_REAL && !isfinite(val->data.real)) {
                expr_free(val);
                tv_clear(&prev);
                return NULL;
            }
            TermVec out;
            tv_init(&out);
            cf_float_expr(val, nlimit, &out);
            expr_free(val);
            if ((long)out.n >= nlimit) {
                if (have_prev && tv_prefix_eq(&prev, &out, nlimit)) {
                    Expr* r = tv_to_list(&out);   /* exactly nlimit terms */
                    tv_clear(&out);
                    tv_clear(&prev);
                    return r;
                }
                /* remember this precision's terms and require confirmation */
                tv_clear(&prev);
                tv_init(&prev);
                for (size_t i = 0; i < out.n; i++) tv_push(&prev, out.v[i]);
                have_prev = true;
            }
            tv_clear(&out);
            continue;                /* not yet confirmed -> more precision */
        }

        /* N produced an exact integer/rational: expand exactly. */
        {
            mpz_t p, q;
            if (expr_get_rational(val, p, q)) {
                TermVec out;
                tv_init(&out);
                cf_rational(p, q, nlimit, &out);
                Expr* r = tv_to_list(&out);
                tv_clear(&out);
                mpz_clears(p, q, (mpz_ptr)0);
                expr_free(val);
                tv_clear(&prev);
                return r;
            }
        }
        /* Not a real number -> cannot expand. */
        expr_free(val);
        tv_clear(&prev);
        return NULL;
    }
    tv_clear(&prev);
    return NULL;
#else
    /* Machine-only build: a single N pass at machine precision. */
    Expr* nexpr = expr_new_function(expr_new_symbol(SYM_N),
                                    (Expr*[]){ expr_copy((Expr*)x) }, 1);
    Expr* val = eval_and_free(nexpr);
    if (val && val->type == EXPR_REAL && isfinite(val->data.real)) {
        Expr* r = cf_inexact(val, nlimit);
        expr_free(val);
        return r;
    }
    if (val) expr_free(val);
    return NULL;
#endif
}

/* ------------------------------------------------------------------ */
/* Dispatcher.                                                         */
/* ------------------------------------------------------------------ */
Expr* builtin_continued_fraction(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t ac = res->data.function.arg_count;
    if (ac < 1 || ac > 2) return NULL;

    Expr* x = res->data.function.args[0];

    long nlimit = 0;                 /* 0 == "all determinable terms" */
    if (ac == 2) {
        Expr* nn = res->data.function.args[1];
        if (nn->type != EXPR_INTEGER || nn->data.integer <= 0) return NULL;
        nlimit = (long)nn->data.integer;
    }

    /* 1. exact rational */
    {
        mpz_t p, q;
        if (expr_get_rational(x, p, q)) {
            TermVec out;
            tv_init(&out);
            cf_rational(p, q, nlimit, &out);
            Expr* r = tv_to_list(&out);
            tv_clear(&out);
            mpz_clears(p, q, (mpz_ptr)0);
            return r;
        }
    }

    /* 2. quadratic surd Sqrt[D] */
    {
        mpz_t D;
        if (expr_is_sqrt_int(x, D)) {
            Expr* r = cf_sqrt_build(D, nlimit);
            mpz_clear(D);
            return r;
        }
    }

    /* 3. inexact real */
    if (x->type == EXPR_REAL
#ifdef USE_MPFR
        || x->type == EXPR_MPFR
#endif
       ) {
        if (x->type == EXPR_REAL && !isfinite(x->data.real)) return NULL;
        return cf_inexact(x, nlimit);
    }

    /* 4. exact symbolic real with explicit n -> adaptive numeric */
    if (nlimit > 0) {
        return cf_exact_numeric(x, nlimit);
    }

    /* 5. cannot determine (e.g. ContinuedFraction[Pi] with no count) */
    return NULL;
}

/* ================================================================== */
/* FromContinuedFraction — the inverse of ContinuedFraction.           */
/* ================================================================== */
/*
 * FromContinuedFraction[{a1, a2, ..., an}] reconstructs the value
 *     a1 + 1/(a2 + 1/(a3 + ... + 1/an)).
 * The terms a_i may be symbolic; the result is the convergent h_n / k_n
 * built from the fundamental recurrence
 *     h_i = a_i h_{i-1} + h_{i-2},   k_i = a_i k_{i-1} + k_{i-2},
 * with h_{-1}=1, h_{-2}=0, k_{-1}=0, k_{-2}=1.  Building the numerator and
 * denominator separately (rather than collapsing the nested fraction) yields
 * the Horner-nested polynomial form, e.g.
 *     FromContinuedFraction[{a,b,c,d}]
 *         = (1 + a b + (a + (1 + a b) c) d) / (b + (1 + b c) d).
 *
 * FromContinuedFraction[{a1, ..., am, {b1, ..., bk}}] reconstructs the exact
 * quadratic irrational whose continued-fraction terms begin with the a_i and
 * then cycle through the b_i forever.  The purely periodic tail
 * x = [b1; b2, ..., bk, b1, ...] satisfies x = (h_{k-1} x + h_{k-2}) /
 * (k_{k-1} x + k_{k-2}) (h, k the period's convergents), i.e. the quadratic
 *     k_{k-1} x^2 + (k_{k-2} - h_{k-1}) x - h_{k-2} = 0,
 * whose positive root x = (P + Q sqrt(R)) / S is then pushed through the
 * leading terms by a Mobius transform (H x + H')/(K x + K') and rationalised
 * to a single (P + Q sqrt(R)) / S in lowest terms.  All leading and period
 * terms must be exact integers in this form.
 */

/* Build f[a]. */
static Expr* fcf_un(const char* head, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(head), args, 1);
}
/* Build f[a, b]. */
static Expr* fcf_bin(const char* head, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(head), args, 2);
}

/* Is e a List[...] expression? */
static bool fcf_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* --- Simple (possibly symbolic) convergent h_{n-1} / k_{n-1}. --- */
static Expr* fcf_simple(Expr** terms, size_t n) {
    if (n == 0) return expr_new_integer(0);

    Expr* hp2 = expr_new_integer(0);   /* h_{-2} */
    Expr* hp1 = expr_new_integer(1);   /* h_{-1} */
    Expr* kp2 = expr_new_integer(1);   /* k_{-2} */
    Expr* kp1 = expr_new_integer(0);   /* k_{-1} */

    for (size_t i = 0; i < n; i++) {
        Expr* t = terms[i];
        /* h_i = t h_{i-1} + h_{i-2}, evaluated each step to keep trees small
         * (numeric terms collapse; symbolic ones stay in convergent form). */
        Expr* h = eval_and_free(fcf_bin("Plus",
            fcf_bin("Times", expr_copy(t), expr_copy(hp1)), expr_copy(hp2)));
        Expr* k = eval_and_free(fcf_bin("Plus",
            fcf_bin("Times", expr_copy(t), expr_copy(kp1)), expr_copy(kp2)));
        expr_free(hp2); hp2 = hp1; hp1 = h;
        expr_free(kp2); kp2 = kp1; kp1 = k;
    }

    /* result = h_{n-1} / k_{n-1}; consumes hp1 and kp1. */
    Expr* result = eval_and_free(fcf_bin("Times", hp1,
        fcf_bin("Power", kp1, expr_new_integer(-1))));
    expr_free(hp2);
    expr_free(kp2);
    return result;
}

/* Largest-square extraction: in = f^2 * sf with sf squarefree (best effort).
 * f and sf are mpz_init'd by the caller.  Trial division handles every prime
 * factor up to FCF_SF_CAP; a residual that is itself a perfect square is then
 * folded into f.  (A residual of the form p^2 q with both p, q > FCF_SF_CAP is
 * left un-reduced — astronomically unlikely for reconstructed CF data.) */
#define FCF_SF_CAP 1000000UL
static void fcf_extract_square(const mpz_t in, mpz_t f, mpz_t sf) {
    mpz_t r, dd;
    mpz_init_set(r, in);
    mpz_init(dd);
    mpz_set_ui(f, 1);
    mpz_set_ui(sf, 1);

    for (unsigned long d = 2; d <= FCF_SF_CAP; d++) {
        mpz_set_ui(dd, d);
        mpz_mul(dd, dd, dd);            /* d^2 */
        if (mpz_cmp(dd, r) > 0) break;  /* d^2 > remaining -> done */
        if (mpz_divisible_ui_p(r, d)) {
            int e = 0;
            while (mpz_divisible_ui_p(r, d)) { mpz_divexact_ui(r, r, d); e++; }
            for (int j = 0; j < e / 2; j++) mpz_mul_ui(f, f, d);
            if (e & 1) mpz_mul_ui(sf, sf, d);
        }
    }
    if (mpz_cmp_ui(r, 1) > 0) {
        if (mpz_perfect_square_p(r)) {
            mpz_t s; mpz_init(s);
            mpz_sqrt(s, r);
            mpz_mul(f, f, s);
            mpz_clear(s);
        } else {
            mpz_mul(sf, sf, r);
        }
    }
    mpz_clear(r);
    mpz_clear(dd);
}

/* Emit the value (P + Q sqrt(R)) / S in lowest terms.  R is squarefree; the
 * value is rational when Q == 0 or R == 1.  Consumes nothing; reads P,Q,S,R
 * (and may mutate them as scratch).  Returns NULL only if S == 0. */
static Expr* fcf_qirr_to_expr(mpz_t P, mpz_t Q, mpz_t S, mpz_t R) {
    if (mpz_sgn(S) == 0) return NULL;
    if (mpz_sgn(S) < 0) { mpz_neg(P, P); mpz_neg(Q, Q); mpz_neg(S, S); }

    bool rational = (mpz_sgn(Q) == 0) || (mpz_cmp_ui(R, 1) == 0);
    if (rational) {
        if (mpz_cmp_ui(R, 1) == 0 && mpz_sgn(Q) != 0) {  /* sqrt(1) = 1 */
            mpz_add(P, P, Q);
            mpz_set_ui(Q, 0);
        }
        Expr* num = mpz_to_expr(P);
        if (mpz_cmp_ui(S, 1) == 0) return num;
        return eval_and_free(fcf_bin("Times", num,
            fcf_bin("Power", mpz_to_expr(S), expr_new_integer(-1))));
    }

    /* Irrational: divide P, Q, S by their common gcd (keeps S > 0). */
    mpz_t g;
    mpz_init(g);
    mpz_gcd(g, P, Q);
    mpz_gcd(g, g, S);
    if (mpz_cmp_ui(g, 1) > 0) {
        mpz_divexact(P, P, g);
        mpz_divexact(Q, Q, g);
        mpz_divexact(S, S, g);
    }
    mpz_clear(g);

    /* inner = P + Q Sqrt[R]; the evaluator drops a zero P and a unit Q. */
    Expr* inner = eval_and_free(fcf_bin("Plus", mpz_to_expr(P),
        fcf_bin("Times", mpz_to_expr(Q), fcf_un("Sqrt", mpz_to_expr(R)))));
    if (mpz_cmp_ui(S, 1) == 0) return inner;
    return eval_and_free(fcf_bin("Times",
        fcf_bin("Power", mpz_to_expr(S), expr_new_integer(-1)), inner));
}

/* --- Periodic (eventually cyclic) -> exact quadratic irrational. --- */
static Expr* fcf_periodic(Expr** lead, size_t m, Expr** per, size_t k) {
    /* Period convergents: after the loop hp1=h_{k-1}, hp2=h_{k-2}, etc. */
    mpz_t hp2, hp1, kp2, kp1, h, kk;
    mpz_inits(hp2, hp1, kp2, kp1, h, kk, (mpz_ptr)0);
    mpz_set_ui(hp2, 0); mpz_set_ui(hp1, 1);
    mpz_set_ui(kp2, 1); mpz_set_ui(kp1, 0);
    for (size_t i = 0; i < k; i++) {
        mpz_t t; expr_to_mpz(per[i], t);
        mpz_mul(h, t, hp1);  mpz_add(h, h, hp2);
        mpz_mul(kk, t, kp1); mpz_add(kk, kk, kp2);
        mpz_set(hp2, hp1); mpz_set(hp1, h);
        mpz_set(kp2, kp1); mpz_set(kp1, kk);
        mpz_clear(t);
    }

    /* Quadratic A x^2 + B x + C = 0 for the purely periodic tail. */
    mpz_t A, B, C, Delta, P, Q, S, R, f, sf, tmp;
    mpz_inits(A, B, C, Delta, P, Q, S, R, f, sf, tmp, (mpz_ptr)0);
    mpz_set(A, kp1);                 /* k_{k-1} */
    mpz_sub(B, kp2, hp1);            /* k_{k-2} - h_{k-1} */
    mpz_neg(C, hp2);                 /* -h_{k-2} */
    mpz_mul(Delta, B, B);
    mpz_mul(tmp, A, C); mpz_mul_ui(tmp, tmp, 4);
    mpz_sub(Delta, Delta, tmp);      /* B^2 - 4AC */

    Expr* result = NULL;
    if (mpz_sgn(A) != 0 && mpz_sgn(Delta) >= 0) {
        mpz_neg(P, B);               /* x = (-B + Q sqrt(R)) / S */
        mpz_mul_ui(S, A, 2);
        if (mpz_sgn(Delta) == 0) {
            mpz_set_ui(Q, 0); mpz_set_ui(R, 1);
        } else {
            fcf_extract_square(Delta, f, sf);
            mpz_set(Q, f); mpz_set(R, sf);
            if (mpz_cmp_ui(R, 1) == 0) {     /* perfect square -> rational */
                mpz_add(P, P, Q); mpz_set_ui(Q, 0);
            }
        }

        /* Leading-term convergents (Mobius (H x + H') / (K x + K')). */
        mpz_t Hp2, Hp1, Kp2, Kp1, th, tk;
        mpz_inits(Hp2, Hp1, Kp2, Kp1, th, tk, (mpz_ptr)0);
        mpz_set_ui(Hp2, 0); mpz_set_ui(Hp1, 1);
        mpz_set_ui(Kp2, 1); mpz_set_ui(Kp1, 0);
        for (size_t i = 0; i < m; i++) {
            mpz_t t; expr_to_mpz(lead[i], t);
            mpz_mul(th, t, Hp1); mpz_add(th, th, Hp2);
            mpz_mul(tk, t, Kp1); mpz_add(tk, tk, Kp2);
            mpz_set(Hp2, Hp1); mpz_set(Hp1, th);
            mpz_set(Kp2, Kp1); mpz_set(Kp1, tk);
            mpz_clear(t);
        }

        /* (a + b sqrt R) / (c + d sqrt R), rationalised. */
        mpz_t a, b, c, d, Pn, Qn, Sn, t1;
        mpz_inits(a, b, c, d, Pn, Qn, Sn, t1, (mpz_ptr)0);
        mpz_mul(a, Hp1, P); mpz_mul(t1, Hp2, S); mpz_add(a, a, t1);
        mpz_mul(b, Hp1, Q);
        mpz_mul(c, Kp1, P); mpz_mul(t1, Kp2, S); mpz_add(c, c, t1);
        mpz_mul(d, Kp1, Q);
        mpz_mul(Pn, a, c); mpz_mul(t1, b, d); mpz_mul(t1, t1, R); mpz_sub(Pn, Pn, t1);
        mpz_mul(Qn, b, c); mpz_mul(t1, a, d); mpz_sub(Qn, Qn, t1);
        mpz_mul(Sn, c, c); mpz_mul(t1, d, d); mpz_mul(t1, t1, R); mpz_sub(Sn, Sn, t1);
        mpz_set(P, Pn); mpz_set(Q, Qn); mpz_set(S, Sn);
        mpz_clears(a, b, c, d, Pn, Qn, Sn, t1, (mpz_ptr)0);
        mpz_clears(Hp2, Hp1, Kp2, Kp1, th, tk, (mpz_ptr)0);

        result = fcf_qirr_to_expr(P, Q, S, R);
    }

    mpz_clears(hp2, hp1, kp2, kp1, h, kk, (mpz_ptr)0);
    mpz_clears(A, B, C, Delta, P, Q, S, R, f, sf, tmp, (mpz_ptr)0);
    return result;
}

Expr* builtin_from_continued_fraction(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    if (!fcf_is_list(arg)) return NULL;

    Expr** terms = arg->data.function.args;
    size_t n = arg->data.function.arg_count;
    if (n == 0) return expr_new_integer(0);   /* {} -> 0 */

    /* A trailing sub-list marks the cyclic (period) block. */
    bool periodic = fcf_is_list(terms[n - 1]);
    for (size_t i = 0; i + 1 < n; i++)
        if (fcf_is_list(terms[i])) return NULL;   /* sub-list only allowed last */

    if (!periodic) return fcf_simple(terms, n);

    /* Periodic: leading terms[0..n-2], period = terms[n-1]. */
    Expr* block = terms[n - 1];
    Expr** per = block->data.function.args;
    size_t k = block->data.function.arg_count;
    if (k == 0) return NULL;                   /* empty period block */

    /* The quadratic-irrational path requires exact integer terms throughout. */
    for (size_t i = 0; i + 1 < n; i++)
        if (!expr_is_integer_like(terms[i])) return NULL;
    for (size_t i = 0; i < k; i++)
        if (!expr_is_integer_like(per[i])) return NULL;

    return fcf_periodic(terms, n - 1, per, k);
}

/* ------------------------------------------------------------------ */
/* Registration.                                                       */
/* ------------------------------------------------------------------ */
void contfrac_init(void) {
    symtab_add_builtin("ContinuedFraction", builtin_continued_fraction);
    symtab_get_def("ContinuedFraction")->attributes |=
        ATTR_LISTABLE | ATTR_PROTECTED;

    symtab_add_builtin("FromContinuedFraction", builtin_from_continued_fraction);
    symtab_get_def("FromContinuedFraction")->attributes |= ATTR_PROTECTED;
}
