/*
 * int.c
 *
 * Integer-digit / integer-length / digit-count builtins.
 *
 *   IntegerDigits[n]            decimal digits of |n|, most significant first.
 *   IntegerDigits[n, b]         base-b digits of |n|.
 *   IntegerDigits[n, b, len]    pads with leading zeros to length len; if
 *                               |n| has more than len base-b digits the
 *                               result is the len least-significant digits.
 *
 *   IntegerLength[n]            number of base-10 digits in |n|.
 *   IntegerLength[n, b]         number of base-b digits in |n|, b >= 2.
 *
 *   DigitCount[n]               list {c(1), c(2), ..., c(9), c(0)} of
 *                               digit counts in the base-10 expansion of |n|.
 *   DigitCount[n, b]            same, but in base b, list of length b
 *                               ordered as {c(1), ..., c(b-1), c(0)}.
 *   DigitCount[n, b, d]         number of times digit d appears in the
 *                               base-b expansion of |n|.  Scalar result.
 *
 * Sign of n is discarded.  IntegerDigits[0] -> {0}; IntegerLength[0] -> 0;
 * DigitCount[0] -> {0,0,...,0} of length b -- zero has no significant
 * digits (matching Mathematica).  IntegerDigits / IntegerLength thread
 * automatically via the Listable attribute; DigitCount is intentionally
 * NOT Listable (only Protected), so DigitCount[{1,2,3}] is left
 * unevaluated rather than threading.
 *
 * All arithmetic is done in GMP, so both machine integers and arbitrary-
 * precision bignums are supported uniformly.
 */

#include "int.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

void int_init(void) {
    symtab_add_builtin("IntegerDigits", builtin_integerdigits);
    symtab_get_def("IntegerDigits")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("IntegerLength", builtin_integerlength);
    symtab_get_def("IntegerLength")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("DigitCount", builtin_digitcount);
    symtab_get_def("DigitCount")->attributes |= ATTR_PROTECTED;
}

/* Build an Expr from an mpz_t digit, demoting to EXPR_INTEGER if it fits. */
static Expr* expr_from_mpz_digit(const mpz_t d) {
    if (mpz_fits_slong_p(d)) {
        return expr_new_integer((int64_t)mpz_get_si(d));
    }
    return expr_new_bigint_from_mpz(d);
}

/* Print `IntegerDigits::argb: IntegerDigits called with N argument(s);
 * between 1 and 3 arguments are expected.` to stderr.  Return NULL so the
 * evaluator leaves the call unevaluated, matching Mathematica's surface
 * behaviour. */
static Expr* int_emit_argb(size_t argc) {
    fprintf(stderr,
            "IntegerDigits::argb: IntegerDigits called with %zu argument%s; "
            "between 1 and 3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Print `IntegerDigits::int: Integer expected at position <pos> in <call>.`
 * to stderr.  `res` is the original call expression so the user sees the
 * exact form they typed echoed back. */
static Expr* int_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "IntegerDigits::int: Integer expected at position %zu in %s.\n",
            pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

Expr* builtin_integerdigits(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return int_emit_argb(argc);

    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Symbolic n flows through silently so downstream rewrites can
         * still rewire the call; concrete-but-non-integer (Real, Rational,
         * Complex, ...) gets a Mathematica-compatible IntegerDigits::int
         * diagnostic before the call is left unevaluated. */
        if (expr_is_numeric_like(n_expr)) return int_emit_int(1, res);
        return NULL;
    }

    /* --- base argument --------------------------------------------- */
    mpz_t base;
    if (argc >= 2) {
        Expr* b_expr = res->data.function.args[1];
        if (!expr_is_integer_like(b_expr)) {
            if (expr_is_numeric_like(b_expr)) return int_emit_int(2, res);
            return NULL;
        }
        expr_to_mpz(b_expr, base);
        if (mpz_cmp_ui(base, 2) < 0) {
            fprintf(stderr,
                "IntegerDigits::ibase: Base argument must be an integer >= 2.\n");
            mpz_clear(base);
            return NULL;
        }
    } else {
        mpz_init_set_ui(base, 10);
    }

    /* --- length argument ------------------------------------------- */
    bool has_len = false;
    size_t target_len = 0;
    if (argc >= 3) {
        Expr* l_expr = res->data.function.args[2];
        if (!expr_is_integer_like(l_expr)) {
            mpz_clear(base);
            if (expr_is_numeric_like(l_expr)) return int_emit_int(3, res);
            return NULL;
        }
        mpz_t l;
        expr_to_mpz(l_expr, l);
        if (mpz_sgn(l) < 0 || !mpz_fits_ulong_p(l)) {
            fprintf(stderr,
                "IntegerDigits::intnn: Non-negative machine-sized integer "
                "expected at position 3.\n");
            mpz_clear(l);
            mpz_clear(base);
            return NULL;
        }
        target_len = (size_t)mpz_get_ui(l);
        mpz_clear(l);
        has_len = true;
    }

    /* --- n = |n| ---------------------------------------------------- */
    mpz_t n;
    expr_to_mpz(n_expr, n);
    mpz_abs(n, n);

    /* --- divmod loop, least-significant digit first ----------------- */
    /* `digits[0]` is the least-significant digit (n mod base), `digits[1]`
     * is the next, etc.  Length-bounded by O(log_base |n|) so a single
     * geometric realloc keeps amortised cost linear. */
    mpz_t* digits = NULL;
    size_t cap = 0, count = 0;

    if (mpz_sgn(n) == 0) {
        cap = 1;
        digits = malloc(sizeof(mpz_t) * cap);
        mpz_init_set_ui(digits[0], 0);
        count = 1;
    } else {
        mpz_t q, r;
        mpz_init(q);
        mpz_init(r);
        while (mpz_sgn(n) > 0) {
            if (count >= cap) {
                size_t new_cap = (cap == 0) ? 16 : cap * 2;
                mpz_t* tmp = realloc(digits, sizeof(mpz_t) * new_cap);
                if (!tmp) {
                    /* Out-of-memory: clean up and bail. */
                    for (size_t i = 0; i < count; i++) mpz_clear(digits[i]);
                    free(digits);
                    mpz_clear(q);
                    mpz_clear(r);
                    mpz_clear(n);
                    mpz_clear(base);
                    return NULL;
                }
                digits = tmp;
                cap = new_cap;
            }
            mpz_tdiv_qr(q, r, n, base);
            mpz_init_set(digits[count++], r);
            mpz_set(n, q);
        }
        mpz_clear(q);
        mpz_clear(r);
    }

    /* --- assemble output -------------------------------------------- */
    /* Without len:    output is exactly `count` digits, MSD-first.
     * With len >= count: prepend (len - count) zeros.
     * With len < count: keep only the `len` least-significant digits
     *                   (i.e. the first `len` entries of `digits`).
     */
    size_t output_len;
    size_t leading_zeros;
    size_t digits_to_use;

    if (has_len) {
        output_len = target_len;
        if (count >= target_len) {
            leading_zeros = 0;
            digits_to_use = target_len;
        } else {
            leading_zeros = target_len - count;
            digits_to_use = count;
        }
    } else {
        output_len = count;
        leading_zeros = 0;
        digits_to_use = count;
    }

    Expr** out_args = NULL;
    if (output_len > 0) {
        out_args = malloc(sizeof(Expr*) * output_len);
        for (size_t i = 0; i < leading_zeros; i++) {
            out_args[i] = expr_new_integer(0);
        }
        /* Emit MSD-first by iterating from `digits_to_use - 1` down to 0. */
        for (size_t i = 0; i < digits_to_use; i++) {
            const mpz_t* d = (const mpz_t*) &digits[digits_to_use - 1 - i];
            out_args[leading_zeros + i] = expr_from_mpz_digit(*d);
        }
    }

    /* --- free working state ---------------------------------------- */
    for (size_t i = 0; i < count; i++) mpz_clear(digits[i]);
    free(digits);
    mpz_clear(n);
    mpz_clear(base);

    Expr* list_head = expr_new_symbol("List");
    Expr* result = expr_new_function(list_head, out_args, output_len);
    free(out_args);
    return result;
}

/* =====================================================================
 * IntegerLength
 * ===================================================================*/

/* Print `IntegerLength::argt: IntegerLength called with N arguments;
 * 1 or 2 arguments are expected.` to stderr.  Mathematica uses `argt`
 * for the variable-arity case (1 OR 2 expected). */
static Expr* intlen_emit_argt(size_t argc) {
    fprintf(stderr,
            "IntegerLength::argt: IntegerLength called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Print `IntegerLength::int: Integer expected at position <pos> in <call>.`
 * matching Mathematica's surface diagnostic. */
static Expr* intlen_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "IntegerLength::int: Integer expected at position %zu in %s.\n",
            pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* Compute the number of base-b digits of |n|.  Returns 0 when n == 0
 * (matching Mathematica: IntegerLength[0] is 0).
 *
 * Strategy:
 *   - Fast path for bases 2..62: mpz_sizeinbase is exact for power-of-2
 *     bases and at most 1 too big otherwise.  When it might be 1 too big
 *     we verify by comparing |n| with base^(size-1), which costs one
 *     pow + one cmp -- O(log size) ulimb operations -- independent of
 *     the number of digits.
 *   - Fallback for arbitrary-precision bases (>= 63): repeated division.
 *     This is O(d^2) bit-ops but only kicks in for the rare arbitrary-
 *     base case.
 *
 * Caller owns `n` and `base`; this routine does not modify or clear them. */
static size_t intlen_count_digits(const mpz_t n, const mpz_t base) {
    if (mpz_sgn(n) == 0) return 0;

    if (mpz_cmp_ui(base, 62) <= 0) {
        unsigned long b = mpz_get_ui(base);
        size_t s = mpz_sizeinbase(n, (int)b);
        /* Power-of-two bases: mpz_sizeinbase is exact. */
        if ((b & (b - 1UL)) == 0) return s;
        /* Otherwise s is either exact or 1 too big; verify. */
        if (s == 0) return 0;  /* defensive: shouldn't happen for n != 0 */
        mpz_t threshold, abs_n;
        mpz_init(threshold);
        mpz_init(abs_n);
        mpz_abs(abs_n, n);
        mpz_ui_pow_ui(threshold, b, (unsigned long)(s - 1));
        int cmp = mpz_cmp(abs_n, threshold);
        mpz_clear(threshold);
        mpz_clear(abs_n);
        return (cmp < 0) ? s - 1 : s;
    }

    /* Arbitrary-precision base: repeated division. */
    mpz_t abs_n, q;
    mpz_init(abs_n);
    mpz_init(q);
    mpz_abs(abs_n, n);
    size_t count = 0;
    while (mpz_sgn(abs_n) > 0) {
        mpz_tdiv_q(q, abs_n, base);
        mpz_swap(abs_n, q);
        count++;
    }
    mpz_clear(q);
    mpz_clear(abs_n);
    return count;
}

Expr* builtin_integerlength(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return intlen_emit_argt(argc);

    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Symbolic n flows through silently; concrete non-integer types
         * (Real, Rational, Complex, ...) trigger the ::int diagnostic
         * before the call is left unevaluated. */
        if (expr_is_numeric_like(n_expr)) return intlen_emit_int(1, res);
        return NULL;
    }

    mpz_t base;
    if (argc == 2) {
        Expr* b_expr = res->data.function.args[1];
        if (!expr_is_integer_like(b_expr)) {
            if (expr_is_numeric_like(b_expr)) return intlen_emit_int(2, res);
            return NULL;
        }
        expr_to_mpz(b_expr, base);
        if (mpz_cmp_ui(base, 2) < 0) {
            fprintf(stderr,
                "IntegerLength::ibase: Base argument must be an integer >= 2.\n");
            mpz_clear(base);
            return NULL;
        }
    } else {
        mpz_init_set_ui(base, 10);
    }

    mpz_t n;
    expr_to_mpz(n_expr, n);
    size_t len = intlen_count_digits(n, base);
    mpz_clear(n);
    mpz_clear(base);

    /* `len` cannot exceed mpz bit-size, which fits in size_t.  size_t can
     * exceed INT64_MAX only on theoretical 128-bit platforms; on all real
     * targets the cast is lossless and the value lands well within
     * EXPR_INTEGER range. */
    return expr_new_integer((int64_t)len);
}

/* =====================================================================
 * DigitCount
 * ===================================================================*/

/* Print `DigitCount::argb: DigitCount called with N argument(s);
 * between 1 and 3 arguments are expected.` matching Mathematica's
 * surface diagnostic exactly. */
static Expr* dc_emit_argb(size_t argc) {
    fprintf(stderr,
            "DigitCount::argb: DigitCount called with %zu argument%s; "
            "between 1 and 3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Print `DigitCount::int: Integer expected at position <pos> in <call>.`
 * Used for non-integer numeric first arg.  Pure symbolic inputs flow
 * through silently (returning NULL) without touching this path. */
static Expr* dc_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "DigitCount::int: Integer expected at position %zu in %s.\n",
            pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* Print `DigitCount::base: The base <b> at position 2 of <call> should
 * be an integer greater than 1.` matching Mathematica's surface
 * diagnostic verbatim.  Used for non-integer base (e.g. 5/2) AND for
 * integer base < 2 -- Mathematica collapses both into the same message. */
static Expr* dc_emit_base(Expr* b_expr, Expr* res) {
    char* b_str = expr_to_string(b_expr);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "DigitCount::base: The base %s at position 2 of %s should be "
            "an integer greater than 1.\n",
            b_str ? b_str : "?",
            call_str ? call_str : "?");
    free(b_str);
    free(call_str);
    return NULL;
}

/* Print `DigitCount::digit: The digit <d> at position 3 of <call>
 * should be a non-negative integer less than the base.`  Used for
 * non-integer digit OR integer digit out of [0, base). */
static Expr* dc_emit_digit(Expr* d_expr, Expr* res) {
    char* d_str = expr_to_string(d_expr);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "DigitCount::digit: The digit %s at position 3 of %s should "
            "be a non-negative integer less than the base.\n",
            d_str ? d_str : "?",
            call_str ? call_str : "?");
    free(d_str);
    free(call_str);
    return NULL;
}

/* Hard cap on the base for the list-returning form (1 / 2-arg).
 * The output list has size b, so an unbounded b would let a malicious
 * input request gigabytes of integer slots.  2^20 = ~1M entries is
 * already excessive but covers every realistic use; the 3-arg form
 * with a bignum digit handles the >2^20 case anyway. */
#define DC_MAX_BASE_FOR_LIST ((unsigned long)1 << 20)

/* Build a per-digit histogram of |n| in base b into `counts[0..b-1]`.
 * `counts` must be zero-initialised by the caller.  Uses GMP's
 * `mpz_tdiv_q_ui` which returns the remainder for ulong divisors --
 * roughly twice as fast as the full mpz_tdiv_qr path when b is small.
 * Bignum bases are not used here (the 1/2-arg form caps b at
 * DC_MAX_BASE_FOR_LIST, which always fits in unsigned long). */
static void dc_build_histogram(const mpz_t n, unsigned long b_val,
                                int64_t* counts) {
    if (mpz_sgn(n) == 0) return;
    mpz_t abs_n, q;
    mpz_init(q);
    mpz_init(abs_n);
    mpz_abs(abs_n, n);
    while (mpz_sgn(abs_n) > 0) {
        unsigned long rem = mpz_tdiv_q_ui(q, abs_n, b_val);
        counts[rem]++;
        mpz_swap(abs_n, q);
    }
    mpz_clear(q);
    mpz_clear(abs_n);
}

/* Count occurrences of `target` in the base-b expansion of |n|.
 * Works for bignum base AND bignum target, so the 3-arg form imposes
 * no upper bound on either.  Returns the count via *out_count, which
 * must already be mpz_init'd by the caller. */
static void dc_count_one_digit(const mpz_t n, const mpz_t base,
                                const mpz_t target, mpz_t out_count) {
    mpz_set_ui(out_count, 0);
    if (mpz_sgn(n) == 0) return;

    /* Fast path: base fits in unsigned long AND target fits in unsigned
     * long.  Avoid an mpz_t for the per-digit remainder. */
    if (mpz_fits_ulong_p(base) && mpz_fits_ulong_p(target)) {
        unsigned long b = mpz_get_ui(base);
        unsigned long t = mpz_get_ui(target);
        mpz_t abs_n, q;
        mpz_init(q);
        mpz_init(abs_n);
        mpz_abs(abs_n, n);
        while (mpz_sgn(abs_n) > 0) {
            unsigned long rem = mpz_tdiv_q_ui(q, abs_n, b);
            if (rem == t) mpz_add_ui(out_count, out_count, 1);
            mpz_swap(abs_n, q);
        }
        mpz_clear(q);
        mpz_clear(abs_n);
        return;
    }

    /* Slow path: bignum base (or bignum target). */
    mpz_t abs_n, q, r;
    mpz_init(q);
    mpz_init(r);
    mpz_init(abs_n);
    mpz_abs(abs_n, n);
    while (mpz_sgn(abs_n) > 0) {
        mpz_tdiv_qr(q, r, abs_n, base);
        if (mpz_cmp(r, target) == 0) {
            mpz_add_ui(out_count, out_count, 1);
        }
        mpz_swap(abs_n, q);
    }
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(abs_n);
}

Expr* builtin_digitcount(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return dc_emit_argb(argc);

    /* --- n -------------------------------------------------------- */
    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Symbolic n flows through silently; concrete non-integer
         * numerics (Real / Rational / Complex) trip the ::int
         * diagnostic.  Matches IntegerLength / IntegerDigits. */
        if (expr_is_numeric_like(n_expr)) return dc_emit_int(1, res);
        return NULL;
    }

    /* --- base ----------------------------------------------------- */
    mpz_t base;
    if (argc >= 2) {
        Expr* b_expr = res->data.function.args[1];
        if (!expr_is_integer_like(b_expr)) {
            /* Mathematica collapses fractional / real / complex bases
             * into ::base, not ::int.  Symbolic bases stay silent. */
            if (expr_is_numeric_like(b_expr)) return dc_emit_base(b_expr, res);
            return NULL;
        }
        expr_to_mpz(b_expr, base);
        if (mpz_cmp_ui(base, 2) < 0) {
            dc_emit_base(b_expr, res);
            mpz_clear(base);
            return NULL;
        }
    } else {
        mpz_init_set_ui(base, 10);
    }

    /* --- 3-arg form: scalar count -------------------------------- */
    if (argc == 3) {
        Expr* d_expr = res->data.function.args[2];
        if (!expr_is_integer_like(d_expr)) {
            mpz_clear(base);
            if (expr_is_numeric_like(d_expr)) return dc_emit_digit(d_expr, res);
            return NULL;
        }
        mpz_t target;
        expr_to_mpz(d_expr, target);
        /* Digit must satisfy 0 <= d < base. */
        if (mpz_sgn(target) < 0 || mpz_cmp(target, base) >= 0) {
            dc_emit_digit(d_expr, res);
            mpz_clear(target);
            mpz_clear(base);
            return NULL;
        }

        mpz_t n;
        expr_to_mpz(n_expr, n);
        mpz_t count;
        mpz_init(count);
        dc_count_one_digit(n, base, target, count);
        mpz_clear(n);
        mpz_clear(target);
        mpz_clear(base);

        Expr* result;
        if (mpz_fits_slong_p(count)) {
            result = expr_new_integer((int64_t)mpz_get_si(count));
        } else {
            result = expr_new_bigint_from_mpz(count);
        }
        mpz_clear(count);
        return result;
    }

    /* --- 1/2-arg form: per-digit histogram ----------------------- */
    /* Need to allocate a length-`base` output list -- enforce the
     * soft cap to avoid OOM on absurd bases.  We still emit a
     * Mathematica-style diagnostic before bailing. */
    if (!mpz_fits_ulong_p(base) || mpz_get_ui(base) > DC_MAX_BASE_FOR_LIST) {
        char* call_str = expr_to_string(res);
        fprintf(stderr,
                "DigitCount::ovfl: Base too large for the list-returning "
                "form in %s; use DigitCount[n, b, d] for a single digit "
                "instead.\n",
                call_str ? call_str : "?");
        free(call_str);
        mpz_clear(base);
        return NULL;
    }

    unsigned long b_val = mpz_get_ui(base);

    /* calloc gives zero-initialised counts, which is exactly what
     * dc_build_histogram expects.  Each count is bounded by the
     * total number of digits ~ log_b(|n|) which fits in int64 for
     * every plausible bignum. */
    int64_t* counts = calloc((size_t)b_val, sizeof(int64_t));
    if (!counts) {
        mpz_clear(base);
        return NULL;
    }

    mpz_t n;
    expr_to_mpz(n_expr, n);
    dc_build_histogram(n, b_val, counts);
    mpz_clear(n);
    mpz_clear(base);

    /* Reorder to {c(1), c(2), ..., c(b-1), c(0)}: digit 0 comes LAST,
     * matching Mathematica's surface convention. */
    Expr** out_args = malloc(sizeof(Expr*) * (size_t)b_val);
    for (unsigned long i = 1; i < b_val; i++) {
        out_args[i - 1] = expr_new_integer(counts[i]);
    }
    out_args[b_val - 1] = expr_new_integer(counts[0]);
    free(counts);

    Expr* list_head = expr_new_symbol("List");
    Expr* result = expr_new_function(list_head, out_args, (size_t)b_val);
    free(out_args);
    return result;
}
