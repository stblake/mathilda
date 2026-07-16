/*
 * real.c
 *
 * RealDigits builtin -- positional-notation digit expansion.
 *
 *   RealDigits[x]               default base 10, length set by Precision[x].
 *   RealDigits[x, b]            base b, length set by Precision[x] / Log10[b].
 *   RealDigits[x, b, len]       exactly `len` digits, MSD-first.
 *   RealDigits[x, b, len, n]    `len` digits, first one = coefficient of b^n.
 *
 * Result form is `{ digits-list, exp }`.  The first element of digits-list is
 * the coefficient of b^(exp-1).  Sign of x is discarded.  Exact rationals
 * with non-terminating base-b expansions return a list ending in a nested
 * list of the recurring block.  Inexact (machine or MPFR) reals get
 * Indeterminate for any requested digit beyond the available precision.
 *
 *   x can be: Integer, BigInt, Rational[n,d], Real (machine), or
 *             EXPR_MPFR (arbitrary precision; USE_MPFR builds only).
 *
 * The general algorithm scales |x| by base^(-low) where `low` is the lowest
 * digit position we need, floors to an integer N, and reads off the base-b
 * digits of N (padding with leading zeros as needed).  This single GMP /
 * MPFR shift handles every numeric type uniformly.  For the special case of
 * an exact rational with no explicit `len`, a remainder-tracked long
 * division detects terminating vs recurring expansions and emits the
 * nested-list form.
 *
 * Implementation only supports integer bases b >= 2.  Non-integer bases
 * (e.g. GoldenRatio) emit a `::ibase` diagnostic and leave the call
 * unevaluated -- adding them requires MPFR floor-iteration and has been
 * deferred.
 */

#include "real.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "sym_names.h"
#include "numeric.h"
#include "arithmetic.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

#ifndef M_LN2
#define M_LN2 0.6931471805599453
#endif
#ifndef M_LN10
#define M_LN10 2.302585092994046
#endif

/* Machine double mantissa width in bits (IEEE-754 binary64). */
#define REAL_MACHINE_BITS 53

/* Hard caps to keep adversarial inputs from blowing memory.  The
 * digit-count cap covers every realistic Precision[x]/Log10[b] return; the
 * cycle-detection cap stops period-finding on pathological rationals whose
 * denominator has a very large multiplicative order. */
#define REAL_MAX_LEN          1000000UL
#define REAL_MAX_CYCLE_DIGITS 100000UL

/* -------------------------------------------------------------------------
 *  Diagnostics.  All emit a Mathematica-style `RealDigits::<tag>` message
 *  to stderr and return NULL so the evaluator leaves the call unevaluated.
 * ----------------------------------------------------------------------- */

static Expr* rd_emit_argb(size_t argc) {
    fprintf(stderr,
            "RealDigits::argb: RealDigits called with %zu argument%s; "
            "between 1 and 4 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

static Expr* rd_emit_realx(Expr* res) {
    char* s = expr_to_string(res);
    fprintf(stderr,
            "RealDigits::realx: The number %s is not a real number.\n",
            s ? s : "?");
    free(s);
    return NULL;
}

static Expr* rd_emit_ibase(Expr* b_expr, Expr* res) {
    char* bs = expr_to_string(b_expr);
    char* cs = expr_to_string(res);
    fprintf(stderr,
            "RealDigits::ibase: Base %s in %s is not a real number greater "
            "than 1.\n",
            bs ? bs : "?", cs ? cs : "?");
    free(bs);
    free(cs);
    return NULL;
}

static Expr* rd_emit_intnn(size_t pos, Expr* res) {
    char* cs = expr_to_string(res);
    fprintf(stderr,
            "RealDigits::intnn: Non-negative machine-sized integer expected "
            "at position %zu in %s.\n",
            pos, cs ? cs : "?");
    free(cs);
    return NULL;
}

static Expr* rd_emit_int(size_t pos, Expr* res) {
    char* cs = expr_to_string(res);
    fprintf(stderr,
            "RealDigits::int: Integer expected at position %zu in %s.\n",
            pos, cs ? cs : "?");
    free(cs);
    return NULL;
}

/* -------------------------------------------------------------------------
 *  Internal helpers
 * ----------------------------------------------------------------------- */

/* Classification of the first argument's numeric kind. */
typedef enum {
    RD_KIND_INT,        /* EXPR_INTEGER, EXPR_BIGINT, or Rational[n,d] with d=1 */
    RD_KIND_RAT,        /* Rational[n,d] with d != 1 */
    RD_KIND_REAL,       /* EXPR_REAL (machine double) */
    RD_KIND_MPFR,       /* EXPR_MPFR */
    RD_KIND_OTHER       /* not a real number */
} RdKind;

/* True for any Rational[n,d] with integer-like components.  Wider than
 * is_rational, which only accepts machine-sized parts. */
static bool rd_is_rational_bigint(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != SYM_Rational
        || e->data.function.arg_count != 2) {
        return false;
    }
    return expr_is_integer_like(e->data.function.args[0])
        && expr_is_integer_like(e->data.function.args[1]);
}

static RdKind rd_classify(const Expr* e) {
    if (!e) return RD_KIND_OTHER;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return RD_KIND_INT;
    if (e->type == EXPR_REAL) return RD_KIND_REAL;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return RD_KIND_MPFR;
#endif
    if (rd_is_rational_bigint(e)) {
        /* If denominator is exactly 1 we treat it as an integer. */
        Expr* d = e->data.function.args[1];
        if (d->type == EXPR_INTEGER && d->data.integer == 1) return RD_KIND_INT;
        if (d->type == EXPR_BIGINT && mpz_cmp_ui(d->data.bigint, 1) == 0) return RD_KIND_INT;
        return RD_KIND_RAT;
    }
    return RD_KIND_OTHER;
}

/* Number of significant base-b digits available in a numeric expression.
 *   - Integer / BigInt / exact Rational: returns 0 to signal "exact, no cap".
 *   - Machine Real: ~53 / log2(b).
 *   - MPFR Real:   prec_bits / log2(b).
 *
 * Rounded to the nearest non-negative integer.  The base is the integer
 * value `b_ulong`; callers already validated that the base fits in
 * unsigned long. */
static size_t rd_available_digits(const Expr* x, unsigned long b_ulong) {
    double log2b = log2((double)b_ulong);
    if (log2b <= 0.0) return 0; /* defensive */
    if (x->type == EXPR_REAL) {
        double d = (double)REAL_MACHINE_BITS / log2b;
        return (size_t)(d + 0.5);
    }
#ifdef USE_MPFR
    if (x->type == EXPR_MPFR) {
        double d = (double)mpfr_get_prec(x->data.mpfr) / log2b;
        return (size_t)(d + 0.5);
    }
#endif
    return 0;
}

/* Build an Expr from an mpz that is known to hold a base-b digit. */
static Expr* rd_digit_expr_from_mpz(const mpz_t d) {
    if (mpz_fits_slong_p(d)) return expr_new_integer((int64_t)mpz_get_si(d));
    return expr_new_bigint_from_mpz(d);
}

/* Extract |x| into the mpq `out`.  Caller must mpq_init `out` first.
 * Pre: rd_classify(x) is RD_KIND_INT or RD_KIND_RAT. */
static void rd_to_mpq(const Expr* x, mpq_t out) {
    if (x->type == EXPR_INTEGER) {
        int64_t v = x->data.integer;
        if (v < 0) v = -v;
        mpq_set_ui(out, (unsigned long)v, 1);
        return;
    }
    if (x->type == EXPR_BIGINT) {
        mpq_set_z(out, x->data.bigint);
        mpq_abs(out, out);
        return;
    }
    /* Rational[n, d] */
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    expr_to_mpz(x->data.function.args[0], n);
    expr_to_mpz(x->data.function.args[1], d);
    mpz_abs(n, n);
    mpz_abs(d, d);
    mpq_set_num(out, n);
    mpq_set_den(out, d);
    mpq_canonicalize(out);
    mpz_clear(n);
    mpz_clear(d);
}

/* Natural exponent of a positive rational q in base b: the (1+floor(log_b q))
 * such that b^(exp-1) <= q < b^exp.  Returns 0 for q == 0 (caller handles
 * the zero case separately).  Computed via integer arithmetic: exp =
 * IntegerLength(floor(q), b) if floor(q) > 0, else -shift where shift is
 * the number of times we multiplied numerator by b before
 * numerator >= denominator.  Caller passes q already in canonical form. */
static long rd_rational_natural_exp(const mpq_t q, const mpz_t base) {
    if (mpq_sgn(q) == 0) return 0;
    /* Try integer part. */
    mpz_t ip, num, den;
    mpz_init(ip);
    mpz_init(num);
    mpz_init(den);
    mpz_set(num, mpq_numref(q));
    mpz_set(den, mpq_denref(q));
    mpz_tdiv_q(ip, num, den);

    if (mpz_sgn(ip) > 0) {
        /* digits of ip in base b (count).  Use the same trick as int.c. */
        size_t s;
        if (mpz_fits_ulong_p(base) && mpz_get_ui(base) <= 62) {
            unsigned long b = mpz_get_ui(base);
            s = mpz_sizeinbase(ip, (int)b);
            if ((b & (b - 1UL)) != 0 && s > 0) {
                mpz_t threshold;
                mpz_init(threshold);
                mpz_ui_pow_ui(threshold, b, (unsigned long)(s - 1));
                if (mpz_cmp(ip, threshold) < 0) s -= 1;
                mpz_clear(threshold);
            }
        } else {
            /* Slow path for huge integer base: repeated division. */
            mpz_t q2, ip_copy;
            mpz_init(q2);
            mpz_init_set(ip_copy, ip);
            s = 0;
            while (mpz_sgn(ip_copy) > 0) {
                mpz_tdiv_q(q2, ip_copy, base);
                mpz_swap(ip_copy, q2);
                s++;
            }
            mpz_clear(q2);
            mpz_clear(ip_copy);
        }
        mpz_clear(ip);
        mpz_clear(num);
        mpz_clear(den);
        return (long)s;
    }

    /* ip == 0, so 0 < q < 1.  Count how many times we multiply num by b
     * before num >= den; that count is the number of leading zeros after
     * the decimal point.  Exp = -(count). */
    long shift = 0;
    mpz_t cur;
    mpz_init_set(cur, num);
    while (mpz_cmp(cur, den) < 0) {
        mpz_mul(cur, cur, base);
        shift++;
        /* Cap absurd magnitudes.  log_b(q) cannot exceed e.g. 10^9 for
         * any input we care about; the cap protects against pathological
         * mpq with denominators 10^(very large). */
        if (shift > 10000000L) break;
    }
    mpz_clear(cur);
    mpz_clear(ip);
    mpz_clear(num);
    mpz_clear(den);
    return -(shift - 1);
}

/* Natural exponent of |x| in base b for inexact (machine or MPFR) x.
 * Returns LONG_MIN as a sentinel for "x is zero". */
static long rd_inexact_natural_exp(const Expr* x, unsigned long b_ulong) {
    double log_b = log((double)b_ulong);
    if (x->type == EXPR_REAL) {
        double v = fabs(x->data.real);
        if (v == 0.0) return LONG_MIN;
        return (long)floor(log(v) / log_b) + 1;
    }
#ifdef USE_MPFR
    if (x->type == EXPR_MPFR) {
        if (mpfr_zero_p(x->data.mpfr)) return LONG_MIN;
        /* mpfr_get_exp returns the exponent e such that x = m * 2^e with
         * 1/2 <= |m| < 1.  Equivalently, 2^(e-1) <= |x| < 2^e, so
         * log2(|x|) is in [e-1, e).  Converting to base b: log_b(|x|) =
         * log2(|x|) / log2(b).  We compute via mpfr_log which is exact
         * enough for the floor. */
        mpfr_t t;
        mpfr_init2(t, 64);
        mpfr_abs(t, x->data.mpfr, MPFR_RNDN);
        mpfr_log(t, t, MPFR_RNDN);
        double v = mpfr_get_d(t, MPFR_RNDN);
        mpfr_clear(t);
        return (long)floor(v / log_b) + 1;
    }
#endif
    return 0;
}

/* Generate digits of an exact rational q (q >= 0) in base b via long
 * division with remainder tracking.  Used only when no explicit `len` is
 * given AND x is exact.  Builds a flat List of digits, possibly ending
 * with a nested List of the recurring block; sets *out_exp to the
 * exponent.
 *
 * Algorithm:
 *   ip = q.num / q.den;  rem = q.num mod q.den;
 *   if (ip > 0): emit IntegerDigits[ip, base] then start tracking rem;
 *                exp = IntegerLength(ip, base).
 *   else if (rem == 0): emit {0}, exp = 0; (zero case)
 *   else: shift rem left in base b until rem * b >= den, counting leading
 *         zero positions to set exp negatively.  Then start tracking.
 *
 *   Tracking: while rem > 0
 *       if rem already seen at position p: cycle found; split frac into
 *           prefix [0..p) and cycle [p..end].
 *       record rem -> position p of next digit
 *       rem *= base
 *       digit = rem / den
 *       rem  %= den
 *       emit digit
 *
 *   If rem reaches 0: terminating; emit flat list.
 *   If position count exceeds REAL_MAX_CYCLE_DIGITS: bail; emit flat list
 *     truncated to the cap (with no Indeterminate -- the caller's explicit
 *     `len` mode is the way to get arbitrary length).
 */
static Expr* rd_rational_exact_form(const mpq_t q, const mpz_t base, long* out_exp) {
    if (mpq_sgn(q) == 0) {
        Expr** a = malloc(sizeof(Expr*));
        a[0] = expr_new_integer(0);
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), a, 1);
        free(a);
        *out_exp = 0;
        return list;
    }

    mpz_t num, den, ip, rem, scratch;
    mpz_init(num);
    mpz_init(den);
    mpz_init(ip);
    mpz_init(rem);
    mpz_init(scratch);
    mpz_set(num, mpq_numref(q));
    mpz_set(den, mpq_denref(q));
    mpz_tdiv_qr(ip, rem, num, den);

    /* Collect integer-part digits MSD-first using the same algorithm as
     * IntegerDigits. */
    Expr** int_digits = NULL;
    size_t int_count = 0;
    long exp_val;

    if (mpz_sgn(ip) > 0) {
        /* divmod loop, LSD-first into a stack. */
        size_t cap = 16;
        mpz_t* stack = malloc(sizeof(mpz_t) * cap);
        size_t n = 0;
        mpz_t q2, r2, ip_copy;
        mpz_init(q2);
        mpz_init(r2);
        mpz_init_set(ip_copy, ip);
        while (mpz_sgn(ip_copy) > 0) {
            if (n >= cap) {
                cap *= 2;
                stack = realloc(stack, sizeof(mpz_t) * cap);
            }
            mpz_tdiv_qr(q2, r2, ip_copy, base);
            mpz_init_set(stack[n++], r2);
            mpz_swap(ip_copy, q2);
        }
        mpz_clear(q2);
        mpz_clear(r2);
        mpz_clear(ip_copy);
        int_count = n;
        int_digits = malloc(sizeof(Expr*) * n);
        /* MSD-first.  Cleanup happens in a separate pass so we do not
         * read from an already-cleared mpz slot. */
        for (size_t i = 0; i < n; i++) {
            int_digits[i] = rd_digit_expr_from_mpz(stack[n - 1 - i]);
        }
        for (size_t i = 0; i < n; i++) mpz_clear(stack[i]);
        free(stack);
        exp_val = (long)int_count;
    } else if (mpz_sgn(rem) == 0) {
        /* q == 0 already handled. */
        exp_val = 0;
    } else {
        /* 0 < q < 1.  Shift rem until rem >= den, counting positions. */
        exp_val = 0;
        while (mpz_cmp(rem, den) < 0) {
            mpz_mul(rem, rem, base);
            exp_val--;
            if (exp_val < -10000000L) break; /* sanity */
        }
        /* `rem` is now numerator scaled so the next digit-extraction step
         * yields a non-zero digit; we re-do the extraction in the loop
         * below, so divide rem back by base to keep the invariant
         * "loop will extract digit by rem *= base; digit = rem/den". */
        mpz_tdiv_q(rem, rem, base);
        /* exp_val currently is one too many shifts (we shifted into the
         * first significant position then backed off by one).  The first
         * digit emitted by the fractional loop is the coefficient of
         * b^(exp_val), which after the back-off matches the natural exp. */
        exp_val += 1;
        /* exp_val now equals the exponent (first significant digit is
         * coefficient of b^(exp_val - 1)).  This matches the integer-part
         * path where exp_val = int_count. */
    }

    /* Fractional-loop with remainder tracking. */
    typedef struct {
        mpz_t r;
        size_t pos;
    } RemEntry;
    RemEntry* seen = NULL;
    size_t seen_count = 0, seen_cap = 0;

    Expr** frac_digits = NULL;
    size_t frac_cap = 0, frac_count = 0;
    size_t cycle_start = 0;
    bool cycle = false;

    while (mpz_sgn(rem) > 0 && frac_count < REAL_MAX_CYCLE_DIGITS) {
        /* Look up rem in `seen`. */
        bool found = false;
        size_t hit = 0;
        for (size_t i = 0; i < seen_count; i++) {
            if (mpz_cmp(seen[i].r, rem) == 0) { found = true; hit = i; break; }
        }
        if (found) {
            cycle_start = seen[hit].pos;
            cycle = true;
            break;
        }
        if (seen_count >= seen_cap) {
            seen_cap = seen_cap ? seen_cap * 2 : 16;
            seen = realloc(seen, sizeof(RemEntry) * seen_cap);
        }
        mpz_init_set(seen[seen_count].r, rem);
        seen[seen_count].pos = frac_count;
        seen_count++;

        mpz_mul(rem, rem, base);
        mpz_tdiv_qr(scratch, rem, rem, den);

        if (frac_count >= frac_cap) {
            frac_cap = frac_cap ? frac_cap * 2 : 16;
            frac_digits = realloc(frac_digits, sizeof(Expr*) * frac_cap);
        }
        frac_digits[frac_count++] = rd_digit_expr_from_mpz(scratch);
    }

    /* Build the resulting digit list.  Three cases:
     *   1. terminating (rem reached 0): flat list of int_digits + frac_digits.
     *   2. cycle detected: int_digits + frac_digits[0..cycle_start) +
     *      [List(frac_digits[cycle_start..end])].
     *   3. cap reached without cycle: degrade to flat int_digits + frac_digits.
     */
    Expr* list_head = expr_new_symbol(SYM_List);
    Expr** out_args;
    size_t out_len;

    if (cycle) {
        /* int_count + cycle_start non-cycle frac digits + 1 sublist. */
        out_len = int_count + cycle_start + 1;
        out_args = malloc(sizeof(Expr*) * out_len);
        for (size_t i = 0; i < int_count; i++) out_args[i] = int_digits[i];
        for (size_t i = 0; i < cycle_start; i++)
            out_args[int_count + i] = frac_digits[i];
        size_t cyc_len = frac_count - cycle_start;
        Expr** cyc_args = malloc(sizeof(Expr*) * cyc_len);
        for (size_t i = 0; i < cyc_len; i++)
            cyc_args[i] = frac_digits[cycle_start + i];
        Expr* sub = expr_new_function(expr_new_symbol(SYM_List), cyc_args, cyc_len);
        free(cyc_args);
        out_args[int_count + cycle_start] = sub;
    } else {
        out_len = int_count + frac_count;
        out_args = malloc(sizeof(Expr*) * (out_len ? out_len : 1));
        for (size_t i = 0; i < int_count; i++) out_args[i] = int_digits[i];
        for (size_t i = 0; i < frac_count; i++)
            out_args[int_count + i] = frac_digits[i];
        if (out_len == 0) {
            /* Shouldn't happen since q != 0 case is the only one here, but
             * be defensive. */
            out_args[0] = expr_new_integer(0);
            out_len = 1;
        }
    }

    Expr* digits_list = expr_new_function(list_head, out_args, out_len);
    free(out_args);

    /* Cleanup. */
    free(int_digits);
    free(frac_digits);
    for (size_t i = 0; i < seen_count; i++) mpz_clear(seen[i].r);
    free(seen);
    mpz_clear(num); mpz_clear(den); mpz_clear(ip); mpz_clear(rem);
    mpz_clear(scratch);

    *out_exp = exp_val;
    return digits_list;
}

/* Exact-only scaled floor: computes floor(|x| * base^(-low)) into `out`
 * when x is Integer, BigInt, or Rational[n, d].  Returns true on success;
 * false for any other input kind.  The exact path uses pure mpz / mpq
 * arithmetic, so the answer is bit-exact regardless of |low|. */
static bool rd_scaled_floor_exact(const Expr* x, const mpz_t base, long low,
                                   mpz_t out) {
    if (x->type != EXPR_INTEGER && x->type != EXPR_BIGINT
        && !rd_is_rational_bigint(x)) {
        return false;
    }
    mpq_t q;
    mpq_init(q);
    rd_to_mpq(x, q);
    if (low >= 0) {
        mpz_t scale, den;
        mpz_init(scale);
        mpz_init(den);
        mpz_ui_pow_ui(scale, mpz_get_ui(base), (unsigned long)low);
        mpz_mul(den, mpq_denref(q), scale);
        mpz_tdiv_q(out, mpq_numref(q), den);
        mpz_clear(scale);
        mpz_clear(den);
    } else {
        mpz_t scale, num;
        mpz_init(scale);
        mpz_init(num);
        mpz_ui_pow_ui(scale, mpz_get_ui(base), (unsigned long)(-low));
        mpz_mul(num, mpq_numref(q), scale);
        mpz_tdiv_q(out, num, mpq_denref(q));
        mpz_clear(scale);
        mpz_clear(num);
    }
    mpq_clear(q);
    return true;
}

#ifdef USE_MPFR
/* Canonical rounded base-b digits for an inexact value.  Returns a
 * heap-allocated array of `count` digit values (each in [0, base)), MSD-
 * first; *out_exp_p is set so the digit at index i is the coefficient of
 * base^(*out_exp_p - 1 - i), i.e. natural_exp = *out_exp_p.
 *
 * The number of digits `count` equals the available precision in base b
 * (= prec_bits / log2(b), rounded).
 *
 * Built on mpfr_get_str, which yields round-to-nearest digits and is the
 * canonical way to recover Mathematica's "shortest decimal representation"
 * for machine reals (e.g. 0.1d -> "1" not "999...").  Caller frees the
 * returned digit array with free(). */
static unsigned long* rd_inexact_canonical_digits(
        const mpfr_t v_abs, mpfr_prec_t prec, unsigned long b_ulong,
        size_t count, long* out_exp_p) {
    if (count == 0) {
        *out_exp_p = 0;
        return NULL;
    }
    if (mpfr_zero_p(v_abs)) {
        unsigned long* arr = calloc(count, sizeof(unsigned long));
        *out_exp_p = 0;
        return arr;
    }
    mpfr_exp_t ep = 0;
    /* mpfr_get_str(NULL, &exp, base, n, op, rnd):
     *   - allocates string of n digits (no leading zeros, no sign).
     *   - implicit radix point left of first digit, so value = 0.ddd...
     *     * base^exp. */
    char* s = mpfr_get_str(NULL, &ep, (int)b_ulong, count, v_abs, MPFR_RNDN);
    if (!s) {
        *out_exp_p = 0;
        unsigned long* arr = calloc(count, sizeof(unsigned long));
        return arr;
    }
    unsigned long* arr = calloc(count, sizeof(unsigned long));
    /* Skip a leading '-' should never happen since we passed v_abs, but
     * guard anyway. */
    const char* p = s;
    if (*p == '-') p++;
    for (size_t i = 0; i < count && p[i]; i++) {
        unsigned char c = (unsigned char)p[i];
        unsigned long d;
        if (c >= '0' && c <= '9') d = (unsigned long)(c - '0');
        else if (c >= 'a' && c <= 'z') d = 10 + (unsigned long)(c - 'a');
        else if (c >= 'A' && c <= 'Z') d = 10 + (unsigned long)(c - 'A');
        else d = 0;
        arr[i] = d;
    }
    mpfr_free_str(s);
    *out_exp_p = (long)ep;
    (void)prec;
    return arr;
}
#endif

/* If x is symbolic (Pi, E, EulerGamma, GoldenRatio, Catalan, Degree, ...)
 * or a function that should be numericalized for digit extraction, try to
 * convert it to an inexact numeric (EXPR_MPFR if available, else EXPR_REAL)
 * at sufficient precision to compute `digits_b` base-b digits.  Returns
 * the new Expr (caller-owned) on success, or NULL if x cannot be
 * numericalized.  Used only when the caller passes explicit `len` so we
 * know how much precision to request. */
static Expr* rd_try_numericalize(const Expr* x, unsigned long b_ulong,
                                  size_t digits_b) {
    NumericSpec spec;
#ifdef USE_MPFR
    /* Bits needed = digits_b * log2(b) + safety margin.  Use 32 guard
     * bits to absorb rounding noise. */
    double log2b = log2((double)b_ulong);
    long bits = (long)((double)digits_b * log2b) + 32;
    if (bits < 64) bits = 64;
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = bits;
#else
    spec = numeric_machine_spec();
#endif
    Expr* n = numericalize(x, spec);
    if (!n) return NULL;
    RdKind k = rd_classify(n);
    if (k == RD_KIND_OTHER) {
        expr_free(n);
        return NULL;
    }
    return n;
}

/* Compute the base-b digits of N, LSD-first, into a freshly-allocated mpz
 * array of length *out_count.  Caller frees the array and clears each mpz. */
static void rd_digits_lsd(const mpz_t N, const mpz_t base,
                           mpz_t** out_arr, size_t* out_count) {
    if (mpz_sgn(N) == 0) {
        mpz_t* arr = malloc(sizeof(mpz_t) * 1);
        mpz_init_set_ui(arr[0], 0);
        *out_arr = arr;
        *out_count = 1;
        return;
    }
    size_t cap = 16, n = 0;
    mpz_t* arr = malloc(sizeof(mpz_t) * cap);
    mpz_t q, r, cur;
    mpz_init(q);
    mpz_init(r);
    mpz_init_set(cur, N);
    while (mpz_sgn(cur) > 0) {
        if (n >= cap) {
            cap *= 2;
            arr = realloc(arr, sizeof(mpz_t) * cap);
        }
        mpz_tdiv_qr(q, r, cur, base);
        mpz_init_set(arr[n++], r);
        mpz_swap(cur, q);
    }
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(cur);
    *out_arr = arr;
    *out_count = n;
}

/* -------------------------------------------------------------------------
 *  Main builtin
 * ----------------------------------------------------------------------- */

Expr* builtin_realdigits(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 4) return rd_emit_argb(argc);

    Expr* x_expr = res->data.function.args[0];

    /* --- Classify x ------------------------------------------------- */
    /* `x_owned`: when we numericalize a symbolic input (Pi, E, …), we
     * stash the freshly allocated Expr here so it is freed before we
     * return. Otherwise it stays NULL and we reference the caller's
     * `x_expr` directly. */
    Expr* x_owned = NULL;
    RdKind kind = rd_classify(x_expr);
    if (kind == RD_KIND_OTHER) {
        /* Concrete-but-non-real (Complex with non-zero imaginary part)
         * emits ::realx.  Anything else (named constants, function
         * expressions, plain symbols) is given to numericalize when the
         * caller has provided enough precision context to make that
         * meaningful (i.e. has_len). */
        Expr *re_p, *im_p;
        if (x_expr->type == EXPR_FUNCTION
            && x_expr->data.function.head
            && x_expr->data.function.head->type == EXPR_SYMBOL
            && x_expr->data.function.head->data.symbol.name == SYM_Complex
            && is_complex(x_expr, &re_p, &im_p)) {
            /* True complex with imaginary != 0 -> not a real number. */
            return rd_emit_realx(res);
        }
        if (expr_is_numeric_like(x_expr)) return rd_emit_realx(res);
        /* Symbolic: only numericalize if we know how many digits to
         * request.  Without that we cannot pick a precision. */
        /* The numericalize precision depends on the (yet-unparsed) base
         * and len arguments; we defer this until after we have parsed
         * them.  Mark kind as OTHER for now and revisit. */
    }

    /* --- Base argument --------------------------------------------- */
    mpz_t base;
    if (argc >= 2) {
        Expr* b_expr = res->data.function.args[1];
        if (!expr_is_integer_like(b_expr)) {
            /* Non-integer base (Real, Rational, etc.) -> ::ibase. */
            if (expr_is_numeric_like(b_expr)) return rd_emit_ibase(b_expr, res);
            return NULL;
        }
        expr_to_mpz(b_expr, base);
        if (mpz_cmp_ui(base, 2) < 0) {
            rd_emit_ibase(b_expr, res);
            mpz_clear(base);
            return NULL;
        }
        if (!mpz_fits_ulong_p(base)) {
            /* Soft cap: bases that don't fit in unsigned long would
             * still work for `floor` but bog down the scaled-floor MPFR
             * shift.  Reject to keep the code paths uniform. */
            fprintf(stderr,
                "RealDigits::ibase: Base too large; use a base <= ULONG_MAX.\n");
            mpz_clear(base);
            return NULL;
        }
    } else {
        mpz_init_set_ui(base, 10);
    }

    /* --- len argument ---------------------------------------------- */
    bool has_len = false;
    bool len_automatic = false;
    size_t target_len = 0;
    if (argc >= 3) {
        Expr* l_expr = res->data.function.args[2];
        if (l_expr->type == EXPR_SYMBOL && l_expr->data.symbol.name == SYM_Automatic) {
            len_automatic = true;
            has_len = true;
        } else if (!expr_is_integer_like(l_expr)) {
            mpz_clear(base);
            if (expr_is_numeric_like(l_expr)) return rd_emit_intnn(3, res);
            return NULL;
        } else {
            mpz_t l;
            expr_to_mpz(l_expr, l);
            if (mpz_sgn(l) < 0 || !mpz_fits_ulong_p(l)
                || mpz_get_ui(l) > REAL_MAX_LEN) {
                rd_emit_intnn(3, res);
                mpz_clear(l);
                mpz_clear(base);
                return NULL;
            }
            target_len = (size_t)mpz_get_ui(l);
            mpz_clear(l);
            has_len = true;
        }
    }

    /* --- start offset (4th arg) ------------------------------------ */
    bool has_n = false;
    long start_n = 0;
    if (argc >= 4) {
        Expr* n_expr = res->data.function.args[3];
        if (!expr_is_integer_like(n_expr)) {
            mpz_clear(base);
            if (expr_is_numeric_like(n_expr)) return rd_emit_int(4, res);
            return NULL;
        }
        mpz_t nz;
        expr_to_mpz(n_expr, nz);
        if (!mpz_fits_slong_p(nz)) {
            fprintf(stderr,
                "RealDigits::int: Starting-position argument too large.\n");
            mpz_clear(nz);
            mpz_clear(base);
            return NULL;
        }
        start_n = mpz_get_si(nz);
        mpz_clear(nz);
        has_n = true;
    }

    unsigned long b_ulong = mpz_get_ui(base);

    /* --- Late numericalize for symbolic input -----------------------
     * If we couldn't classify the first arg as a recognised numeric
     * kind but the caller has given an explicit `len`, ask the
     * `N` machinery to produce an inexact approximation at the
     * matching precision.  This is what lets `RealDigits[Pi, 10, 25]`
     * succeed: a symbolic Pi is numericalized to enough MPFR bits to
     * yield 25 base-10 digits. */
    if (kind == RD_KIND_OTHER) {
        if (!has_len) {
            /* Cannot pick a precision.  Leave call unevaluated; no
             * diagnostic -- this is the documented Mathematica
             * behaviour for `RealDigits[Pi]`. */
            mpz_clear(base);
            return NULL;
        }
        size_t want = len_automatic ? (size_t)32 : target_len;
        if (has_n) {
            /* Need digits down to position start_n - want + 1.  Add the
             * absolute magnitude of start_n so we keep enough leading
             * precision to land in the requested window. */
            long mag = start_n < 0 ? -start_n : start_n;
            want += (size_t)mag;
        }
        x_owned = rd_try_numericalize(x_expr, b_ulong, want);
        if (!x_owned) {
            mpz_clear(base);
            return NULL;
        }
        x_expr = x_owned;
        kind = rd_classify(x_expr);
        if (kind == RD_KIND_OTHER) {
            /* Numericalize handed back something we still don't
             * recognise (e.g. a Plus[..] that didn't fold).  Punt. */
            expr_free(x_owned);
            mpz_clear(base);
            return NULL;
        }
    }

    /* --- Path A: exact (Integer/Rational) with no explicit `len`.
     *
     * Use the recurring-aware long division so terminating expansions
     * are flat and non-terminating ones acquire a nested-list cycle
     * block.  This path ignores the 4-arg form (start_n) since the
     * natural exponent IS the result. */
    if ((kind == RD_KIND_INT || kind == RD_KIND_RAT) && !has_len && !has_n) {
        mpq_t q;
        mpq_init(q);
        rd_to_mpq(x_expr, q);
        long exp_val;
        Expr* digits = rd_rational_exact_form(q, base, &exp_val);
        mpq_clear(q);
        mpz_clear(base);
        if (x_owned) expr_free(x_owned);

        Expr** pair = malloc(sizeof(Expr*) * 2);
        pair[0] = digits;
        pair[1] = expr_new_integer((int64_t)exp_val);
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
        free(pair);
        return result;
    }

    /* --- Path B: general digit extraction via scaled floor.
     *
     * Compute `start_pos` (the position of the first digit) and `len`,
     * then read off len base-b digits of floor(|x| * b^(-low)). */

    /* Natural exponent of |x|.  For zero we follow Mathematica's
     * documented behaviour:
     *   RealDigits[0]  -> {{0}, 0}
     *   RealDigits[0.] -> {{0}, -Floor[Accuracy[0.]]}.
     */
    long natural_exp = 0;
    bool x_is_zero = false;
    if (kind == RD_KIND_INT || kind == RD_KIND_RAT) {
        mpq_t q;
        mpq_init(q);
        rd_to_mpq(x_expr, q);
        if (mpq_sgn(q) == 0) x_is_zero = true;
        else natural_exp = rd_rational_natural_exp(q, base);
        mpq_clear(q);
    } else {
        long e = rd_inexact_natural_exp(x_expr, b_ulong);
        if (e == LONG_MIN) x_is_zero = true;
        else natural_exp = e;
    }

    /* --- Zero special case ----------------------------------------- */
    if (x_is_zero) {
        Expr* zero_list;
        long out_exp;
        if (kind == RD_KIND_INT || kind == RD_KIND_RAT) {
            /* Exact zero: {{0}, 0}. */
            Expr** za = malloc(sizeof(Expr*));
            za[0] = expr_new_integer(0);
            zero_list = expr_new_function(expr_new_symbol(SYM_List), za, 1);
            free(za);
            out_exp = has_n ? start_n + 1 : 0;
            /* With explicit len, pad to len. */
            if (has_len && !len_automatic) {
                expr_free(zero_list);
                Expr** za2 = malloc(sizeof(Expr*) * (target_len ? target_len : 1));
                for (size_t i = 0; i < target_len; i++) za2[i] = expr_new_integer(0);
                zero_list = expr_new_function(expr_new_symbol(SYM_List), za2,
                                               target_len);
                free(za2);
            }
        } else {
            /* Inexact zero: digits are all 0 within accuracy.  We follow
             * Mathematica's form {{0}, -Floor[Accuracy[0.]]}.  Accuracy of
             * machine 0. is MachinePrecision + (-log10($MinMachineNumber))
             * ≈ 323.607, giving -323.  Accuracy of MPFR 0 with p_bits is
             * p_bits / log2(10).  Keep this in sync with accuracy_of() in
             * src/precision.c. */
            double acc_digits;
#ifdef USE_MPFR
            if (x_expr->type == EXPR_MPFR) {
                acc_digits = (double)mpfr_get_prec(x_expr->data.mpfr) / log2(10.0);
            } else
#endif
            {
                acc_digits = (double)REAL_MACHINE_BITS / log2(10.0) - log10(DBL_MIN);
            }
            long acc_floor = (long)floor(acc_digits);
            Expr** za = malloc(sizeof(Expr*));
            za[0] = expr_new_integer(0);
            zero_list = expr_new_function(expr_new_symbol(SYM_List), za, 1);
            free(za);
            out_exp = -acc_floor;
            if (has_n) out_exp = start_n + 1;
            if (has_len && !len_automatic) {
                expr_free(zero_list);
                Expr** za2 = malloc(sizeof(Expr*) * (target_len ? target_len : 1));
                for (size_t i = 0; i < target_len; i++) za2[i] = expr_new_integer(0);
                zero_list = expr_new_function(expr_new_symbol(SYM_List), za2,
                                               target_len);
                free(za2);
            }
        }
        mpz_clear(base);
        if (x_owned) expr_free(x_owned);
        Expr** pair = malloc(sizeof(Expr*) * 2);
        pair[0] = zero_list;
        pair[1] = expr_new_integer((int64_t)out_exp);
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
        free(pair);
        return result;
    }

    /* --- Determine start_pos and len ------------------------------- */
    /* start_pos = exponent of the first digit (i.e. first digit is
     * coefficient of b^start_pos).  Default is natural_exp - 1. */
    long start_pos = has_n ? start_n : (natural_exp - 1);

    /* Effective length:
     *   - explicit non-Automatic len: target_len.
     *   - Automatic (4-arg only meaningful for inexact): clamp to available.
     *   - no len: default to Round[Precision[x]/Log10[b]] for inexact, or
     *     -- for exact --  the natural length from start_pos down to 0
     *     (integer part length).  We only reach this branch when has_n is
     *     true for exact x (otherwise Path A handled it).  Use a default
     *     of natural_exp digits in that case.
     */
    size_t avail = rd_available_digits(x_expr, b_ulong);
    size_t out_len;

    if (has_len && !len_automatic) {
        out_len = target_len;
    } else if (len_automatic) {
        /* Automatic: only emit digits within available precision.  The
         * window is [start_pos .. start_pos - avail + 1].  If the natural
         * exponent's most significant digit is at natural_exp - 1, then
         * the available range in absolute positions is
         *   [natural_exp - 1 .. natural_exp - avail].
         * We clamp [start_pos .. start_pos - out_len + 1] to that range. */
        long lo_avail = natural_exp - (long)avail;
        if (start_pos < lo_avail) {
            /* Entire requested window is below precision -- empty list. */
            out_len = 0;
        } else {
            long count = start_pos - lo_avail + 1;
            if (count < 0) count = 0;
            out_len = (size_t)count;
        }
    } else if (kind == RD_KIND_INT || kind == RD_KIND_RAT) {
        /* has_n true on an exact input but no len: emit natural_exp
         * digits from start_pos down to start_pos - natural_exp + 1.
         * If natural_exp <= 0, emit at least one digit. */
        long count = natural_exp > 0 ? natural_exp : 1;
        out_len = (size_t)count;
    } else {
        out_len = avail;
    }

    /* Cap. */
    if (out_len > REAL_MAX_LEN) {
        fprintf(stderr,
            "RealDigits::ovfl: Requested length %zu exceeds cap %lu.\n",
            out_len, (unsigned long)REAL_MAX_LEN);
        mpz_clear(base);
        if (x_owned) expr_free(x_owned);
        return NULL;
    }

    /* --- Extract digits ------------------------------------------- */
    long low = start_pos - (long)out_len + 1;

    /* Determine the precision-limited window for inexact inputs.  Any
     * requested digit at position < precision_low becomes Indeterminate. */
    bool inexact = (kind == RD_KIND_REAL || kind == RD_KIND_MPFR);

    /* Build the digits.  Special case out_len == 0. */
    Expr** out_args = NULL;
    if (out_len > 0) {
        out_args = malloc(sizeof(Expr*) * out_len);

        if (!inexact) {
            /* Exact (Integer / Rational) path.  Compute N = floor(|x| *
             * b^(-low)) with bit-exact mpz / mpq arithmetic and read off
             * the base-b digits LSD-first. */
            mpz_t N;
            mpz_init(N);
            if (!rd_scaled_floor_exact(x_expr, base, low, N)) {
                mpz_clear(N);
                free(out_args);
                mpz_clear(base);
                if (x_owned) expr_free(x_owned);
                fprintf(stderr,
                    "RealDigits::nrep: Internal error scaling input for "
                    "digit extraction.\n");
                return NULL;
            }
            mpz_t* lsd = NULL;
            size_t lsd_count = 0;
            rd_digits_lsd(N, base, &lsd, &lsd_count);
            /* The MSD-first ordering puts the lsd_count'th digit
             * (positions low + lsd_count - 1) leftmost.  If the window
             * is wider than lsd_count, pad with leading zeros.  If
             * narrower, truncate the high end -- those bits would be
             * outside the user-requested window. */
            size_t use = lsd_count > out_len ? out_len : lsd_count;
            size_t leading_zeros = out_len - use;
            for (size_t i = 0; i < leading_zeros; i++) {
                out_args[i] = expr_new_integer(0);
            }
            for (size_t i = 0; i < use; i++) {
                out_args[leading_zeros + i] =
                    rd_digit_expr_from_mpz(lsd[use - 1 - i]);
            }
            for (size_t i = 0; i < lsd_count; i++) mpz_clear(lsd[i]);
            free(lsd);
            mpz_clear(N);
        } else {
            /* Inexact (Real / MPFR) path.  Use mpfr_get_str to obtain
             * the canonical round-to-nearest digits at the available
             * precision; place them at positions
             *   [exp_p - 1 .. exp_p - avail]
             * in the output window.  Positions above exp_p - 1 become
             * 0 (since x < b^exp_p), positions below exp_p - avail
             * become Indeterminate.  This matches Mathematica's
             * "shortest round-trip" rendering: 0.1 + 0.2 displayed as
             * 0.30000000000000004, not the long binary tail. */
#ifdef USE_MPFR
            /* Build an absolute-value MPFR copy of x. */
            mpfr_t v;
            mpfr_prec_t prec_bits;
            if (x_expr->type == EXPR_REAL) {
                prec_bits = REAL_MACHINE_BITS;
                mpfr_init2(v, prec_bits);
                mpfr_set_d(v, fabs(x_expr->data.real), MPFR_RNDN);
            } else {
                prec_bits = mpfr_get_prec(x_expr->data.mpfr);
                mpfr_init2(v, prec_bits);
                mpfr_abs(v, x_expr->data.mpfr, MPFR_RNDN);
            }
            long exp_p = 0;
            unsigned long* canon = rd_inexact_canonical_digits(
                v, prec_bits, b_ulong, avail, &exp_p);
            mpfr_clear(v);
            /* Map output position -> digit.  out_args[i] is at
             * position start_pos - i. */
            long hi = exp_p - 1;
            long lo_p = exp_p - (long)avail;
            for (size_t i = 0; i < out_len; i++) {
                long pos = start_pos - (long)i;
                if (pos > hi) {
                    out_args[i] = expr_new_integer(0);
                } else if (pos < lo_p) {
                    out_args[i] = expr_new_symbol(SYM_Indeterminate);
                } else {
                    size_t idx = (size_t)(hi - pos);
                    out_args[i] = expr_new_integer(
                        (int64_t)(canon ? canon[idx] : 0));
                }
            }
            free(canon);
#else
            /* No MPFR: best-effort floor extraction via doubles. */
            double v = fabs(x_expr->data.real);
            for (size_t i = 0; i < out_len; i++) {
                long pos = start_pos - (long)i;
                if (pos < (natural_exp - (long)avail)) {
                    out_args[i] = expr_new_symbol(SYM_Indeterminate);
                } else if (pos > natural_exp - 1) {
                    out_args[i] = expr_new_integer(0);
                } else {
                    double scale = pow((double)b_ulong, -(double)pos);
                    unsigned long d = (unsigned long)floor(v * scale)
                                       % b_ulong;
                    out_args[i] = expr_new_integer((int64_t)d);
                }
            }
#endif
        }
    }

    /* --- Assemble output ------------------------------------------ */
    long out_exp = has_n ? (start_n + 1) : natural_exp;

    /* Special-case for inexact x when no `n` was given: the canonical
     * exponent comes from mpfr_get_str (`exp_p`) rather than the rough
     * log_b estimate in natural_exp.  We already used exp_p above when
     * placing digits; reflect it in out_exp too so the result is
     * self-consistent.  This is harmless when natural_exp == exp_p
     * (the typical case) and corrects the rare boundary case where the
     * MPFR rounding pushes the value into the next decade. */
    if (inexact && !has_n) {
#ifdef USE_MPFR
        /* Recompute exp_p without redoing the full mpfr_get_str by
         * asking for a single digit at default rounding.  This is
         * O(log avail) effectively a no-op for typical inputs. */
        mpfr_t v;
        mpfr_prec_t pb;
        if (x_expr->type == EXPR_REAL) {
            pb = REAL_MACHINE_BITS;
            mpfr_init2(v, pb);
            mpfr_set_d(v, fabs(x_expr->data.real), MPFR_RNDN);
        } else {
            pb = mpfr_get_prec(x_expr->data.mpfr);
            mpfr_init2(v, pb);
            mpfr_abs(v, x_expr->data.mpfr, MPFR_RNDN);
        }
        if (!mpfr_zero_p(v)) {
            mpfr_exp_t ep = 0;
            char* s = mpfr_get_str(NULL, &ep, (int)b_ulong, 2, v, MPFR_RNDN);
            if (s) {
                out_exp = (long)ep;
                mpfr_free_str(s);
            }
        }
        mpfr_clear(v);
#endif
    }

    Expr* digits_list = expr_new_function(expr_new_symbol(SYM_List),
                                            out_args, out_len);
    free(out_args);
    mpz_clear(base);
    if (x_owned) expr_free(x_owned);

    Expr** pair = malloc(sizeof(Expr*) * 2);
    pair[0] = digits_list;
    pair[1] = expr_new_integer((int64_t)out_exp);
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
    free(pair);
    return result;
}

/* -------------------------------------------------------------------------
 *  MantissaExponent
 *
 *  MantissaExponent[x]      -> {m, e} such that x = m * 10^e and
 *                              1/10 <= |m| < 1  (or m == 0).
 *  MantissaExponent[x, b]   -> base-b mantissa / exponent.
 *
 *  Numeric kinds accepted: Integer, BigInt, Rational, Real, MPFR.  Complex
 *  values emit `MantissaExponent::realx`; non-integer bases leave the call
 *  unevaluated (current implementation supports integer bases >= 2 only).
 *
 *  For exact inputs the mantissa is built as an exact Rational with the
 *  GMP-canonical form sign(x) * |x_num| / (x_den * b^e).  For inexact
 *  inputs (machine Real or MPFR) the mantissa is computed in the same
 *  numeric type and at the same precision as the input.
 * ----------------------------------------------------------------------- */

static Expr* me_emit_argt(size_t argc) {
    fprintf(stderr,
            "MantissaExponent::argt: MantissaExponent called with %zu "
            "argument%s; 1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

static Expr* me_emit_realx(Expr* x_expr) {
    char* s = expr_to_string(x_expr);
    fprintf(stderr,
            "MantissaExponent::realx: The value %s is not a real number.\n",
            s ? s : "?");
    free(s);
    return NULL;
}

static Expr* me_emit_ibase(Expr* b_expr) {
    char* bs = expr_to_string(b_expr);
    fprintf(stderr,
            "MantissaExponent::ibase: Base %s is not an integer greater "
            "than 1.\n",
            bs ? bs : "?");
    free(bs);
    return NULL;
}

/* Build a List[m, e] result; takes ownership of `m_expr`. */
static Expr* me_make_pair(Expr* m_expr, long e_val) {
    Expr** pair = malloc(sizeof(Expr*) * 2);
    pair[0] = m_expr;
    pair[1] = expr_new_integer((int64_t)e_val);
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
    free(pair);
    return result;
}

/* Convert a canonical mpq to an Expr.  Returns an Integer / BigInt when the
 * denominator is 1, otherwise a `Rational[n, d]` function node with
 * appropriately-typed integer leaves. */
static Expr* me_mpq_to_expr(const mpq_t q) {
    if (mpz_cmp_ui(mpq_denref(q), 1) == 0) {
        if (mpz_fits_slong_p(mpq_numref(q))) {
            return expr_new_integer((int64_t)mpz_get_si(mpq_numref(q)));
        }
        return expr_new_bigint_from_mpz(mpq_numref(q));
    }
    Expr** args = malloc(sizeof(Expr*) * 2);
    if (mpz_fits_slong_p(mpq_numref(q))) {
        args[0] = expr_new_integer((int64_t)mpz_get_si(mpq_numref(q)));
    } else {
        args[0] = expr_new_bigint_from_mpz(mpq_numref(q));
    }
    if (mpz_fits_slong_p(mpq_denref(q))) {
        args[1] = expr_new_integer((int64_t)mpz_get_si(mpq_denref(q)));
    } else {
        args[1] = expr_new_bigint_from_mpz(mpq_denref(q));
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
    free(args);
    return r;
}

/* Load |x|'s absolute value (sign discarded) into `out_q` and write the
 * original sign into `*sign_out` (+1, -1, or 0 for zero).  Caller must
 * mpq_init `out_q` first.  Pre: classified RD_KIND_INT or RD_KIND_RAT. */
static void me_to_signed_mpq(const Expr* x, mpq_t out_q, int* sign_out) {
    if (x->type == EXPR_INTEGER) {
        int64_t v = x->data.integer;
        if (v == 0)      { *sign_out = 0;  mpq_set_ui(out_q, 0, 1); return; }
        if (v < 0)       { *sign_out = -1; v = -v; }
        else             { *sign_out = +1; }
        mpq_set_ui(out_q, (unsigned long)v, 1);
        return;
    }
    if (x->type == EXPR_BIGINT) {
        int s = mpz_sgn(x->data.bigint);
        *sign_out = s;
        if (s == 0) { mpq_set_ui(out_q, 0, 1); return; }
        mpq_set_z(out_q, x->data.bigint);
        mpq_abs(out_q, out_q);
        return;
    }
    /* Rational[n, d] */
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    expr_to_mpz(x->data.function.args[0], n);
    expr_to_mpz(x->data.function.args[1], d);
    int sn = mpz_sgn(n);
    int sd = mpz_sgn(d);
    *sign_out = sn * sd;
    mpz_abs(n, n);
    mpz_abs(d, d);
    mpq_set_num(out_q, n);
    mpq_set_den(out_q, d);
    mpq_canonicalize(out_q);
    mpz_clear(n);
    mpz_clear(d);
}

Expr* builtin_mantissa_exponent(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return me_emit_argt(argc);

    Expr* x_expr = res->data.function.args[0];

    /* --- Reject true complex inputs ------------------------------------ */
    {
        Expr *re_p, *im_p;
        if (x_expr->type == EXPR_FUNCTION
            && x_expr->data.function.head
            && x_expr->data.function.head->type == EXPR_SYMBOL
            && x_expr->data.function.head->data.symbol.name == SYM_Complex
            && is_complex(x_expr, &re_p, &im_p)) {
            return me_emit_realx(x_expr);
        }
    }

    /* --- Base argument ------------------------------------------------- */
    mpz_t base;
    if (argc == 2) {
        Expr* b_expr = res->data.function.args[1];
        if (!expr_is_integer_like(b_expr)) {
            /* Non-integer bases (Real, Rational, symbolic E, …) are not
             * yet implemented.  Leave the call unevaluated rather than
             * raising a misleading diagnostic. */
            return NULL;
        }
        expr_to_mpz(b_expr, base);
        if (mpz_cmp_ui(base, 2) < 0) {
            me_emit_ibase(b_expr);
            mpz_clear(base);
            return NULL;
        }
        if (!mpz_fits_ulong_p(base)) {
            fprintf(stderr,
                "MantissaExponent::ibase: Base too large; use a base "
                "<= ULONG_MAX.\n");
            mpz_clear(base);
            return NULL;
        }
    } else {
        mpz_init_set_ui(base, 10);
    }
    unsigned long b_ulong = mpz_get_ui(base);

    /* --- Classify the input ------------------------------------------- */
    RdKind kind = rd_classify(x_expr);
    if (kind == RD_KIND_OTHER) {
        /* Symbolic, named constant, or unhandled numeric form -- leave
         * the call unevaluated, no diagnostic.  This mirrors the
         * Mathematica behaviour for `MantissaExponent[x]`. */
        mpz_clear(base);
        return NULL;
    }

    /* ====================================================================
     *  EXACT path (Integer / BigInt / Rational)
     * ================================================================== */
    if (kind == RD_KIND_INT || kind == RD_KIND_RAT) {
        mpq_t q;
        mpq_init(q);
        int sign = 0;
        me_to_signed_mpq(x_expr, q, &sign);

        if (sign == 0) {
            /* MantissaExponent[0] = {0, 0}. */
            mpq_clear(q);
            mpz_clear(base);
            return me_make_pair(expr_new_integer(0), 0);
        }

        long e_val = rd_rational_natural_exp(q, base);

        /* mantissa = sign * |q| / b^e_val.  Express as a single mpq by
         * multiplying the denominator by base^e_val.  All arithmetic on
         * non-negative |q|, sign re-applied at the end. */
        mpz_t scale;
        mpz_init(scale);
        if (e_val >= 0) {
            mpz_ui_pow_ui(scale, b_ulong, (unsigned long)e_val);
            mpz_mul(mpq_denref(q), mpq_denref(q), scale);
        } else {
            /* For exact rationals with |x| < 1 the natural exponent is
             * non-positive, so we multiply the numerator by b^(-e). */
            mpz_ui_pow_ui(scale, b_ulong, (unsigned long)(-e_val));
            mpz_mul(mpq_numref(q), mpq_numref(q), scale);
        }
        mpq_canonicalize(q);
        if (sign < 0) mpz_neg(mpq_numref(q), mpq_numref(q));

        Expr* m_expr = me_mpq_to_expr(q);

        mpz_clear(scale);
        mpq_clear(q);
        mpz_clear(base);
        return me_make_pair(m_expr, e_val);
    }

    /* ====================================================================
     *  INEXACT path -- machine double
     * ================================================================== */
    if (kind == RD_KIND_REAL) {
        double v = x_expr->data.real;
        if (v == 0.0) {
            mpz_clear(base);
            return me_make_pair(expr_new_real(0.0), 0);
        }
        double absv = fabs(v);
        double log_b = log((double)b_ulong);
        long e_val = (long)floor(log(absv) / log_b) + 1;

        /* Compute scale = b^e_val.  Off-by-one corrections cover the rare
         * case where floor(log / log) rounds across a power-of-b boundary
         * because of double-rounding in the log itself. */
        double scale = pow((double)b_ulong, (double)e_val);
        double m = v / scale;
        if (fabs(m) >= 1.0) {
            m /= (double)b_ulong;
            e_val++;
        } else if (fabs(m) * (double)b_ulong < 1.0) {
            m *= (double)b_ulong;
            e_val--;
        }
        mpz_clear(base);
        return me_make_pair(expr_new_real(m), e_val);
    }

#ifdef USE_MPFR
    /* ====================================================================
     *  INEXACT path -- MPFR
     * ================================================================== */
    if (kind == RD_KIND_MPFR) {
        mpfr_srcptr xv = x_expr->data.mpfr;
        mpfr_prec_t prec = mpfr_get_prec(xv);
        if (mpfr_zero_p(xv)) {
            mpz_clear(base);
            return me_make_pair(expr_new_mpfr_bits(prec), 0);
        }

        /* Determine the natural exponent with a comfortable margin so
         * floor() can never land on the wrong side of a power-of-b
         * boundary. */
        mpfr_prec_t work = prec + 64;
        mpfr_t lx, lb;
        mpfr_init2(lx, work);
        mpfr_init2(lb, work);
        mpfr_abs(lx, xv, MPFR_RNDN);
        mpfr_log(lx, lx, MPFR_RNDN);
        mpfr_set_ui(lb, b_ulong, MPFR_RNDN);
        mpfr_log(lb, lb, MPFR_RNDN);
        mpfr_div(lx, lx, lb, MPFR_RNDN);
        mpfr_floor(lx, lx);
        long e_val = mpfr_get_si(lx, MPFR_RNDD) + 1;
        mpfr_clear(lx);
        mpfr_clear(lb);

        /* mantissa = x * b^(-e_val), kept at the input precision. */
        mpfr_t m, bp;
        mpfr_init2(m, prec);
        mpfr_init2(bp, prec + 32);
        mpfr_set_ui(bp, b_ulong, MPFR_RNDN);
        mpfr_pow_si(bp, bp, e_val, MPFR_RNDN);
        mpfr_div(m, xv, bp, MPFR_RNDN);
        mpfr_clear(bp);

        /* Off-by-one correction: |m| should be in [1/b, 1). */
        mpfr_t am;
        mpfr_init2(am, prec);
        mpfr_abs(am, m, MPFR_RNDN);
        if (mpfr_cmp_ui(am, 1) >= 0) {
            mpfr_div_ui(m, m, b_ulong, MPFR_RNDN);
            e_val++;
        } else {
            mpfr_t one_over_b;
            mpfr_init2(one_over_b, work);
            mpfr_set_ui(one_over_b, b_ulong, MPFR_RNDN);
            mpfr_ui_div(one_over_b, 1, one_over_b, MPFR_RNDN);
            if (mpfr_cmp(am, one_over_b) < 0) {
                mpfr_mul_ui(m, m, b_ulong, MPFR_RNDN);
                e_val--;
            }
            mpfr_clear(one_over_b);
        }
        mpfr_clear(am);

        Expr* m_expr = expr_new_mpfr_copy(m);
        mpfr_clear(m);
        mpz_clear(base);
        return me_make_pair(m_expr, e_val);
    }
#endif

    mpz_clear(base);
    return NULL;
}

/* -------------------------------------------------------------------------
 *  RealExponent
 *
 *  RealExponent[x]     -> Log[10, |x|]    (always inexact).
 *  RealExponent[x, b]  -> Log[b,  |x|].
 *
 *  Numeric kinds accepted for x and b: Integer, BigInt, Rational, Real,
 *  and (USE_MPFR) MPFR.  Symbolic numeric inputs -- bare constants like
 *  `Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Degree`, and any
 *  numeric-valued composite such as `Pi^Pi` or `1/Pi` -- are routed
 *  through `numericalize` at the combined working precision and then
 *  treated as their numeric equivalents.  Plain symbols (no numeric
 *  value) leave the call unevaluated with no diagnostic.
 *
 *  Output kind:
 *    - x is MPFR: result is MPFR at the higher of x's and b's bit
 *      precisions (b is promoted from machine / exact to that precision).
 *    - b is MPFR but x is not: result is MPFR at b's precision.
 *    - x is Real or b is Real (but neither is MPFR): result is a machine
 *      Real.
 *    - both exact: result is a machine Real.
 *
 *  Zero handling (Mathematica-compatible):
 *    - Exact zero (Integer 0, BigInt 0, Rational 0/n): RealExponent[0]
 *      = -Infinity.
 *    - Machine zero (0.):  log_b($MinMachineNumber).  Unbased call
 *      (b = 10) returns log10(DBL_MIN) ≈ -307.65.
 *    - MPFR zero @ p bits: -(p / log2(10)) / log10(b).  Unbased call
 *      returns -(p / log2(10)) digits, matching MMA's RealExponent[0``p] = -p.
 *
 *  Diagnostics:
 *    - argc not in [1, 2]                  -> ::argt
 *    - Complex x with non-zero imaginary   -> ::realx
 *    - Base is Complex, or base <= 1       -> ::ibase
 * ----------------------------------------------------------------------- */

static Expr* re_emit_argt(size_t argc) {
    fprintf(stderr,
            "RealExponent::argt: RealExponent called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

static Expr* re_emit_realx(Expr* x_expr) {
    char* s = expr_to_string(x_expr);
    fprintf(stderr,
            "RealExponent::realx: The value %s is not a real number.\n",
            s ? s : "?");
    free(s);
    return NULL;
}

static Expr* re_emit_ibase(Expr* b_expr) {
    char* bs = expr_to_string(b_expr);
    fprintf(stderr,
            "RealExponent::ibase: Base %s is not a real number greater "
            "than 1.\n",
            bs ? bs : "?");
    free(bs);
    return NULL;
}

/* -Infinity = Times[-1, Infinity].  The evaluator's Times canonicalisation
 * does not collapse this further; the printer renders it as `-Infinity`. */
static Expr* re_make_minus_infinity(void) {
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_new_integer(-1);
    args[1] = expr_new_symbol(SYM_Infinity);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), args, 2);
    free(args);
    return r;
}

/* True if a recognized numeric expression represents zero. */
static bool re_is_zero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return mpfr_zero_p(e->data.mpfr);
#endif
    if (rd_is_rational_bigint(e)) {
        Expr* n = e->data.function.args[0];
        if (n->type == EXPR_INTEGER) return n->data.integer == 0;
        if (n->type == EXPR_BIGINT)  return mpz_sgn(n->data.bigint) == 0;
    }
    return false;
}

/* Sign of a recognized numeric expression: 1, 0, or -1. */
static int re_sign(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer > 0) return 1;
        if (e->data.integer < 0) return -1;
        return 0;
    }
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    if (e->type == EXPR_REAL) {
        if (e->data.real > 0.0) return 1;
        if (e->data.real < 0.0) return -1;
        return 0;
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_sgn(e->data.mpfr);
#endif
    if (rd_is_rational_bigint(e)) {
        mpz_t n, d;
        expr_to_mpz(e->data.function.args[0], n);
        expr_to_mpz(e->data.function.args[1], d);
        int s = mpz_sgn(n) * mpz_sgn(d);
        mpz_clear(n); mpz_clear(d);
        return s;
    }
    return 0;
}

#ifdef USE_MPFR
/* Set `out` (already mpfr_init2'd) to |e| where e is a recognized numeric
 * kind.  Returns true on success. */
static bool re_mpfr_set_abs(mpfr_t out, const Expr* e) {
    if (e->type == EXPR_INTEGER) {
        long v = (long)e->data.integer;
        if (v < 0) v = -v;
        mpfr_set_si(out, v, MPFR_RNDN);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_abs(tmp, e->data.bigint);
        mpfr_set_z(out, tmp, MPFR_RNDN);
        mpz_clear(tmp);
        return true;
    }
    if (e->type == EXPR_REAL) {
        mpfr_set_d(out, fabs(e->data.real), MPFR_RNDN);
        return true;
    }
    if (e->type == EXPR_MPFR) {
        mpfr_abs(out, e->data.mpfr, MPFR_RNDN);
        return true;
    }
    if (rd_is_rational_bigint(e)) {
        mpz_t n, d;
        expr_to_mpz(e->data.function.args[0], n);
        expr_to_mpz(e->data.function.args[1], d);
        mpz_abs(n, n);
        mpz_abs(d, d);
        mpq_t q;
        mpq_init(q);
        mpq_set_num(q, n);
        mpq_set_den(q, d);
        mpq_canonicalize(q);
        mpfr_set_q(out, q, MPFR_RNDN);
        mpq_clear(q);
        mpz_clear(n); mpz_clear(d);
        return true;
    }
    return false;
}
#else
/* Best-effort plain-double absolute value for the no-MPFR build path. */
static double re_to_double_abs(const Expr* e) {
    if (e->type == EXPR_INTEGER) {
        int64_t v = e->data.integer;
        if (v < 0) v = -v;
        return (double)v;
    }
    if (e->type == EXPR_BIGINT) return fabs(mpz_get_d(e->data.bigint));
    if (e->type == EXPR_REAL)   return fabs(e->data.real);
    if (rd_is_rational_bigint(e)) {
        mpz_t n, d;
        expr_to_mpz(e->data.function.args[0], n);
        expr_to_mpz(e->data.function.args[1], d);
        double v = fabs(mpz_get_d(n) / mpz_get_d(d));
        mpz_clear(n); mpz_clear(d);
        return v;
    }
    return 0.0;
}
#endif

Expr* builtin_real_exponent(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return re_emit_argt(argc);

    Expr* x_in = res->data.function.args[0];
    Expr* b_in = (argc == 2) ? res->data.function.args[1] : NULL;

    /* --- Reject true complex inputs ----------------------------------- */
    {
        Expr *re_p, *im_p;
        if (x_in->type == EXPR_FUNCTION
            && x_in->data.function.head
            && x_in->data.function.head->type == EXPR_SYMBOL
            && x_in->data.function.head->data.symbol.name == SYM_Complex
            && is_complex(x_in, &re_p, &im_p)) {
            return re_emit_realx(x_in);
        }
        if (b_in && b_in->type == EXPR_FUNCTION
            && b_in->data.function.head
            && b_in->data.function.head->type == EXPR_SYMBOL
            && b_in->data.function.head->data.symbol.name == SYM_Complex
            && is_complex(b_in, &re_p, &im_p)) {
            return re_emit_ibase(b_in);
        }
    }

    /* --- Working precision selection --------------------------------- *
     * Numericalize uses these bits when promoting symbolic constants
     * (Pi, E, ...) into MPFR values.  When either x or b is already MPFR
     * we lift the working precision so that downstream Log preserves it. */
#ifdef USE_MPFR
    long working_bits = 64;
    if (x_in->type == EXPR_MPFR) {
        long bp = (long)mpfr_get_prec(x_in->data.mpfr) + 32;
        if (bp > working_bits) working_bits = bp;
    }
    if (b_in && b_in->type == EXPR_MPFR) {
        long bp = (long)mpfr_get_prec(b_in->data.mpfr) + 32;
        if (bp > working_bits) working_bits = bp;
    }
#endif

    /* --- Numericalize symbolic x / b to recognized numeric kinds ----- */
    Expr* x_owned = NULL;
    Expr* b_owned = NULL;
    const Expr* x_use = x_in;
    const Expr* b_use = b_in;

    NumericSpec spec;
#ifdef USE_MPFR
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = working_bits;
#else
    spec = numeric_machine_spec();
#endif

    RdKind kind_x = rd_classify(x_in);
    if (kind_x == RD_KIND_OTHER) {
        x_owned = numericalize(x_in, spec);
        if (!x_owned) return NULL;
        kind_x = rd_classify(x_owned);
        if (kind_x == RD_KIND_OTHER) {
            expr_free(x_owned);
            return NULL;
        }
        x_use = x_owned;
    }

    RdKind kind_b = b_in ? rd_classify(b_in) : RD_KIND_INT;
    if (b_in && kind_b == RD_KIND_OTHER) {
        b_owned = numericalize(b_in, spec);
        if (!b_owned) {
            if (x_owned) expr_free(x_owned);
            return NULL;
        }
        kind_b = rd_classify(b_owned);
        if (kind_b == RD_KIND_OTHER) {
            expr_free(b_owned);
            if (x_owned) expr_free(x_owned);
            return NULL;
        }
        b_use = b_owned;
    }

    /* --- Base must be a real value strictly greater than 1 ----------- */
    if (b_in) {
        if (re_sign(b_use) <= 0) {
            if (x_owned) expr_free(x_owned);
            if (b_owned) expr_free(b_owned);
            return re_emit_ibase(b_in);
        }
        bool too_small = false;
        if (b_use->type == EXPR_INTEGER) {
            if (b_use->data.integer <= 1) too_small = true;
        } else if (b_use->type == EXPR_BIGINT) {
            if (mpz_cmp_ui(b_use->data.bigint, 1) <= 0) too_small = true;
        } else if (b_use->type == EXPR_REAL) {
            if (b_use->data.real <= 1.0) too_small = true;
#ifdef USE_MPFR
        } else if (b_use->type == EXPR_MPFR) {
            if (mpfr_cmp_ui(b_use->data.mpfr, 1) <= 0) too_small = true;
#endif
        } else if (rd_is_rational_bigint(b_use)) {
            /* |b| > 1  iff  |num| > |den|. */
            mpz_t n, d;
            expr_to_mpz(b_use->data.function.args[0], n);
            expr_to_mpz(b_use->data.function.args[1], d);
            mpz_abs(n, n); mpz_abs(d, d);
            if (mpz_cmp(n, d) <= 0) too_small = true;
            mpz_clear(n); mpz_clear(d);
        }
        if (too_small) {
            if (x_owned) expr_free(x_owned);
            if (b_owned) expr_free(b_owned);
            return re_emit_ibase(b_in);
        }
    }

    /* --- Zero case (exact only) --------------------------------------
     * Exact zero -> -Infinity.  Inexact zero (machine 0., MPFR 0) falls
     * through and is handled after output-kind selection so we can match
     * the result's representation to the input precision. */
    bool x_inexact_zero = false;
    if (re_is_zero(x_use)) {
        bool x_is_exact = (kind_x == RD_KIND_INT || kind_x == RD_KIND_RAT);
        if (x_is_exact) {
            if (x_owned) expr_free(x_owned);
            if (b_owned) expr_free(b_owned);
            return re_make_minus_infinity();
        }
        x_inexact_zero = true;
    }

    /* --- Output kind selection ---------------------------------------
     * Only the *original* input matters here.  Numericalize may have
     * promoted a symbolic Pi or E to a 64-bit MPFR purely for computation,
     * but if the user never asked for arbitrary precision the result
     * should come back at MachinePrecision -- matching Mathematica's
     * `RealExponent[Pi^Pi]` giving a machine-precision answer. */
#ifdef USE_MPFR
    long out_bits = 53;
    bool out_is_mpfr = false;
    if (x_in->type == EXPR_MPFR) {
        out_bits = (long)mpfr_get_prec(x_in->data.mpfr);
        out_is_mpfr = true;
    }
    if (b_in && b_in->type == EXPR_MPFR) {
        long bp = (long)mpfr_get_prec(b_in->data.mpfr);
        if (!out_is_mpfr || bp > out_bits) out_bits = bp;
        out_is_mpfr = true;
    }

    /* --- Inexact zero case ------------------------------------------- *
     * Mathematica's convention is RealExponent[0_inexact] = log_b(MinPositive)
     * where MinPositive is the smallest representable nonzero magnitude at
     * the input's precision.  For machine 0. that is DBL_MIN (≈ 2.225e-308),
     * giving log10 ≈ -307.65.  For MPFR 0 @ p bits the analog is 10^(-p_d)
     * where p_d = p/log2(10), giving -p_d.  Division by log10(b) converts
     * to base b. */
    if (x_inexact_zero) {
        double log10_min;
        if (x_use->type == EXPR_MPFR) {
            log10_min = -(double)mpfr_get_prec(x_use->data.mpfr) / log2(10.0);
        } else {
            log10_min = log10(DBL_MIN);
        }
        double log10_b = 1.0;
        if (b_use) {
            mpfr_t bt;
            mpfr_init2(bt, 64);
            re_mpfr_set_abs(bt, b_use);
            log10_b = log10(mpfr_get_d(bt, MPFR_RNDN));
            mpfr_clear(bt);
        }
        double result_val = log10_min / log10_b;
        Expr* result;
        if (out_is_mpfr) {
            mpfr_t out;
            mpfr_init2(out, out_bits);
            mpfr_set_d(out, result_val, MPFR_RNDN);
            result = expr_new_mpfr_move(out);
        } else {
            result = expr_new_real(result_val);
        }
        if (x_owned) expr_free(x_owned);
        if (b_owned) expr_free(b_owned);
        return result;
    }

    long work_bits = (out_is_mpfr ? out_bits : 53) + 32;
    if (work_bits < 64) work_bits = 64;

    mpfr_t mx, mb_v, lx, lb;
    mpfr_init2(mx,   work_bits);
    mpfr_init2(mb_v, work_bits);
    mpfr_init2(lx,   work_bits);
    mpfr_init2(lb,   work_bits);

    re_mpfr_set_abs(mx, x_use);
    if (b_use) re_mpfr_set_abs(mb_v, b_use);
    else       mpfr_set_ui(mb_v, 10, MPFR_RNDN);

    mpfr_log(lx, mx, MPFR_RNDN);
    mpfr_log(lb, mb_v, MPFR_RNDN);
    mpfr_div(lx, lx, lb, MPFR_RNDN);

    Expr* result;
    if (out_is_mpfr) {
        mpfr_t out;
        mpfr_init2(out, out_bits);
        mpfr_set(out, lx, MPFR_RNDN);
        result = expr_new_mpfr_move(out);
    } else {
        result = expr_new_real(mpfr_get_d(lx, MPFR_RNDN));
    }

    mpfr_clear(mx); mpfr_clear(mb_v); mpfr_clear(lx); mpfr_clear(lb);
#else
    /* No MPFR build: machine doubles only.  Inexact zero handled
     * via the same MinPositive convention as the MPFR path. */
    if (x_inexact_zero) {
        double bv = b_use ? re_to_double_abs(b_use) : 10.0;
        double result_val = log10(DBL_MIN) / log10(bv);
        Expr* result = expr_new_real(result_val);
        if (x_owned) expr_free(x_owned);
        if (b_owned) expr_free(b_owned);
        return result;
    }
    double xv = re_to_double_abs(x_use);
    double bv = b_use ? re_to_double_abs(b_use) : 10.0;
    Expr* result = expr_new_real(log(xv) / log(bv));
#endif

    if (x_owned) expr_free(x_owned);
    if (b_owned) expr_free(b_owned);
    return result;
}

/* -------------------------------------------------------------------------
 *  Module init
 * ----------------------------------------------------------------------- */

void real_init(void) {
    symtab_add_builtin("RealDigits", builtin_realdigits);
    symtab_get_def("RealDigits")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("MantissaExponent", builtin_mantissa_exponent);
    symtab_get_def("MantissaExponent")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("RealExponent", builtin_real_exponent);
    symtab_get_def("RealExponent")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);
}
