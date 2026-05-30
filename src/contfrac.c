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
    Expr* list = expr_new_function(expr_new_symbol("List"), args, n);
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
        e->data.function.head->data.symbol == SYM_Rational &&
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
    if (!h || h->type != EXPR_SYMBOL || h->data.symbol != SYM_Power ||
        e->data.function.arg_count != 2)
        return false;
    Expr* base = e->data.function.args[0];
    Expr* ex   = e->data.function.args[1];
    if (!expr_is_integer_like(base)) return false;
    /* exponent must be exactly Rational[1, 2] */
    if (!(ex->type == EXPR_FUNCTION && ex->data.function.head &&
          ex->data.function.head->type == EXPR_SYMBOL &&
          ex->data.function.head->data.symbol == SYM_Rational &&
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
    Expr* result = expr_new_function(expr_new_symbol("List"), args, 2);
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
            expr_new_symbol("N"),
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
    Expr* nexpr = expr_new_function(expr_new_symbol("N"),
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

/* ------------------------------------------------------------------ */
/* Registration.                                                       */
/* ------------------------------------------------------------------ */
void contfrac_init(void) {
    symtab_add_builtin("ContinuedFraction", builtin_continued_fraction);
    symtab_get_def("ContinuedFraction")->attributes |=
        ATTR_LISTABLE | ATTR_PROTECTED;
}
