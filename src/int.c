/*
 * int.c
 *
 * Integer-digit / integer-length / digit-count / digit-assembly builtins.
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
 *   FromDigits[list]            reconstructs a base-10 integer from a list
 *                               of digits, most-significant first.
 *   FromDigits[list, b]         base-b digit assembly.  Digits and base may
 *                               be arbitrary expressions (symbolic, Real,
 *                               negative, larger than b); the all-integer
 *                               case is computed via GMP, otherwise the
 *                               call expands to the Horner polynomial
 *                               d[0]*b^(n-1) + ... + d[n-1] for the
 *                               evaluator to simplify.
 *   FromDigits["string"]        decodes characters in base 10, where the
 *                               characters '0'-'9', 'a'-'z'/'A'-'Z'
 *                               represent digit values 0-9, 10-35.  As with
 *                               the list form, character "digits" can
 *                               exceed the base and are "carried" through
 *                               Horner.
 *   FromDigits["string", b]     decode characters in base b.
 *
 * Sign of n is discarded.  IntegerDigits[0] -> {0}; IntegerLength[0] -> 0;
 * DigitCount[0] -> {0,0,...,0} of length b -- zero has no significant
 * digits (matching Mathematica).  IntegerDigits / IntegerLength thread
 * automatically via the Listable attribute; DigitCount and FromDigits are
 * intentionally NOT Listable (only Protected): the natural first argument
 * of each is itself a list, so threading would be wrong.
 *
 * All arithmetic is done in GMP, so both machine integers and arbitrary-
 * precision bignums are supported uniformly.
 */

#include "int.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "sym_names.h"

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

    symtab_add_builtin("FromDigits", builtin_fromdigits);
    symtab_get_def("FromDigits")->attributes |= ATTR_PROTECTED;
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

/* =====================================================================
 * FromDigits
 *
 * The inverse of IntegerDigits / IntegerString.  Reconstructs a number
 * from a list (or string) of digits in MSD-first order:
 *
 *   FromDigits[{d0, d1, ..., d_{n-1}}, b]
 *     == d0 * b^(n-1) + d1 * b^(n-2) + ... + d_{n-1}
 *
 * Default base is 10.  When every digit AND the base are integer-like,
 * the result is computed via Horner in GMP and demoted to EXPR_INTEGER
 * when it fits in int64.  When any input is symbolic / Real / Rational /
 * Complex, the call expands to the same Horner sum as an Expr tree and
 * is handed back to the evaluator -- one code path handles symbolic
 * bases, symbolic digits, mixed digits, and inexact bases uniformly.
 *
 * Digits in the list (or characters in the string) need not be in
 * [0, b); they are "carried" through the Horner accumulation, matching
 * Mathematica's behaviour (e.g. FromDigits[{7,11,0,0,0,122}] == 810122
 * and FromDigits["1A3C"] == 2042).
 * ===================================================================*/

/* `FromDigits::argb: FromDigits called with N argument(s); 1 or 2
 * arguments are expected.`  Print to stderr and return NULL so the
 * call is left unevaluated, matching Mathematica's surface behaviour. */
static Expr* fd_emit_argb(size_t argc) {
    fprintf(stderr,
            "FromDigits::argb: FromDigits called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* `FromDigits::nlst: Argument <arg> at position 1 in <call> should be a
 * list or a string.`  Fired when the first argument is not a List and
 * not a String AND is concrete enough that leaving it unevaluated would
 * surprise the user (numbers).  Pure symbolic first args (e.g.
 * FromDigits[x]) flow through silently. */
static Expr* fd_emit_nlst(Expr* arg, Expr* res) {
    char* arg_str = expr_to_string(arg);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "FromDigits::nlst: Argument %s at position 1 in %s should "
            "be a list or a string.\n",
            arg_str ? arg_str : "?",
            call_str ? call_str : "?");
    free(arg_str);
    free(call_str);
    return NULL;
}

/* `FromDigits::char: Invalid digit character '<c>' in string at
 * position <pos> in <call>.`  Fired for the string form when a
 * character outside [0-9a-zA-Z] is encountered. */
static Expr* fd_emit_char(int c, size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    /* Render the offending byte readably regardless of its value. */
    char buf[8];
    if (c >= 0x20 && c < 0x7F) {
        snprintf(buf, sizeof(buf), "'%c'", (char)c);
    } else {
        snprintf(buf, sizeof(buf), "\\x%02X", (unsigned)c & 0xFFu);
    }
    fprintf(stderr,
            "FromDigits::char: Invalid digit character %s in string at "
            "position %zu in %s.\n",
            buf, pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* Map a single character to its digit value (0..35).  Returns -1 for
 * any character outside [0-9a-zA-Z].  Uses explicit ranges instead of
 * isdigit / isalpha so the result is locale-independent and matches
 * Mathematica's surface convention exactly. */
static int fd_char_to_digit(int c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'z') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'Z') return 10 + (c - 'A');
    return -1;
}

/* True iff every element of `list` is integer-like (EXPR_INTEGER or
 * EXPR_BIGINT).  `list` must be a List[]-headed function. */
static bool fd_all_integer_digits(Expr* list) {
    size_t n = list->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        if (!expr_is_integer_like(list->data.function.args[i])) return false;
    }
    return true;
}

/* Compute the integer FromDigits result via Horner in GMP, given an
 * already-validated all-integer digit list and an mpz base.  Returns a
 * fresh Expr (demoted to EXPR_INTEGER when it fits in int64, else
 * EXPR_BIGINT).  Empty list -> 0.
 *
 * Note `expr_to_mpz` does an mpz_init_set_* into its output, so the
 * per-iteration `digit` is paired with a matching `mpz_clear` at the
 * bottom of the loop -- skipping that clear leaks an mp_limb_t per
 * digit. */
static Expr* fd_compute_integer_list(Expr* list, const mpz_t base) {
    size_t n = list->data.function.arg_count;
    mpz_t acc, tmp;
    mpz_init_set_ui(acc, 0);
    mpz_init(tmp);
    for (size_t i = 0; i < n; i++) {
        mpz_t digit;
        expr_to_mpz(list->data.function.args[i], digit);
        mpz_mul(tmp, acc, base);
        mpz_add(acc, tmp, digit);
        mpz_clear(digit);
    }
    mpz_clear(tmp);

    Expr* out;
    if (mpz_fits_slong_p(acc)) {
        out = expr_new_integer((int64_t)mpz_get_si(acc));
    } else {
        out = expr_new_bigint_from_mpz(acc);
    }
    mpz_clear(acc);
    return out;
}

/* Compute the integer FromDigits result for a string of characters in
 * an integer base.  Returns a fresh Expr on success; on invalid char
 * emits the `::char` diagnostic and returns NULL with `res` left
 * unevaluated.  Empty string -> 0. */
static Expr* fd_compute_integer_string(const char* s, const mpz_t base,
                                        Expr* res) {
    mpz_t acc, tmp;
    mpz_init_set_ui(acc, 0);
    mpz_init(tmp);
    size_t pos = 0;
    for (const char* p = s; *p; p++, pos++) {
        int d = fd_char_to_digit((unsigned char)*p);
        if (d < 0) {
            mpz_clear(acc);
            mpz_clear(tmp);
            return fd_emit_char((unsigned char)*p, pos + 1, res);
        }
        mpz_mul(tmp, acc, base);
        mpz_add_ui(acc, tmp, (unsigned long)d);
    }
    mpz_clear(tmp);

    Expr* out;
    if (mpz_fits_slong_p(acc)) {
        out = expr_new_integer((int64_t)mpz_get_si(acc));
    } else {
        out = expr_new_bigint_from_mpz(acc);
    }
    mpz_clear(acc);
    return out;
}

/* Build the symbolic Horner expansion
 *   Plus[ Times[d[0], Power[base, n-1]],
 *         Times[d[1], Power[base, n-2]],
 *         ...,
 *         Times[d[n-2], base],
 *         d[n-1] ]
 *
 * The result is left for the evaluator to simplify on the next pass,
 * so digits-larger-than-base, negative digits, Real or symbolic bases,
 * etc. all collapse onto a single code path.  Caller owns `list` and
 * `base_expr`; we deep-copy what we need into the new tree. */
static Expr* fd_build_symbolic(Expr* list, Expr* base_expr) {
    size_t n = list->data.function.arg_count;
    if (n == 0) return expr_new_integer(0);

    Expr** plus_args = malloc(sizeof(Expr*) * n);

    for (size_t i = 0; i < n; i++) {
        Expr* digit = expr_copy(list->data.function.args[i]);
        size_t exponent = n - 1 - i;

        if (exponent == 0) {
            /* Units term: just the digit itself. */
            plus_args[i] = digit;
            continue;
        }

        /* Power[base, exponent].  `exponent` is at most n-1 where n is
         * an arg_count, comfortably within int64 range. */
        Expr** power_args = malloc(sizeof(Expr*) * 2);
        power_args[0] = expr_copy(base_expr);
        power_args[1] = expr_new_integer((int64_t)exponent);
        Expr* power_node = expr_new_function(expr_new_symbol("Power"),
                                              power_args, 2);
        free(power_args);

        /* Times[digit, Power[base, exponent]]. */
        Expr** times_args = malloc(sizeof(Expr*) * 2);
        times_args[0] = digit;
        times_args[1] = power_node;
        Expr* times_node = expr_new_function(expr_new_symbol("Times"),
                                              times_args, 2);
        free(times_args);
        plus_args[i] = times_node;
    }

    /* Single-element list collapses to Plus[x] which the evaluator
     * simplifies to x via the OneIdentity attribute -- still correct
     * even though we could short-circuit it here. */
    Expr* plus_node = expr_new_function(expr_new_symbol("Plus"),
                                         plus_args, n);
    free(plus_args);
    return plus_node;
}

Expr* builtin_fromdigits(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return fd_emit_argb(argc);

    Expr* first = res->data.function.args[0];

    /* --- base argument --------------------------------------------- */
    /* If no explicit base, default to the integer 10.  We don't allocate
     * a fresh Expr for the default -- only the mpz fast path needs the
     * value, and the symbolic path uses the original argument expression
     * (or a fresh expr_new_integer(10) when it falls back). */
    Expr* base_expr = NULL;
    if (argc == 2) base_expr = res->data.function.args[1];

    /* --- Case 1: first argument is a string ------------------------ */
    if (first->type == EXPR_STRING) {
        /* String form: every "digit" is a character with value in
         * [0, 36).  We currently require an integer base.  Symbolic /
         * non-integer base over a string would expand to a polynomial
         * in characters, which is not a meaningful Mathematica form
         * (FromDigits["abc", x] is left unevaluated). */
        mpz_t base;
        bool base_owned = false;
        if (base_expr == NULL) {
            mpz_init_set_ui(base, 10);
            base_owned = true;
        } else if (expr_is_integer_like(base_expr)) {
            expr_to_mpz(base_expr, base);
            base_owned = true;
            if (mpz_cmp_ui(base, 2) < 0) {
                fprintf(stderr,
                    "FromDigits::ibase: Base argument must be an integer >= 2.\n");
                mpz_clear(base);
                return NULL;
            }
        } else {
            /* Symbolic / non-integer base over a string: leave
             * unevaluated.  No diagnostic -- the user is likely
             * deliberately holding the base for later substitution. */
            return NULL;
        }

        Expr* out = fd_compute_integer_string(first->data.string, base, res);
        if (base_owned) mpz_clear(base);
        return out;
    }

    /* --- Case 2: first argument is a List -------------------------- */
    if (first->type == EXPR_FUNCTION
        && first->data.function.head
        && first->data.function.head->type == EXPR_SYMBOL
        && first->data.function.head->data.symbol == SYM_List) {

        /* Decide between integer fast path and symbolic Horner. */
        bool all_int_digits = fd_all_integer_digits(first);
        bool int_base = (base_expr == NULL) || expr_is_integer_like(base_expr);

        if (all_int_digits && int_base) {
            /* Pure-integer path: GMP Horner.  Reject base < 2 with a
             * Mathematica-style diagnostic. */
            mpz_t base;
            if (base_expr == NULL) {
                mpz_init_set_ui(base, 10);
            } else {
                expr_to_mpz(base_expr, base);
                if (mpz_cmp_ui(base, 2) < 0) {
                    fprintf(stderr,
                        "FromDigits::ibase: Base argument must be an "
                        "integer >= 2.\n");
                    mpz_clear(base);
                    return NULL;
                }
            }
            Expr* out = fd_compute_integer_list(first, base);
            mpz_clear(base);
            return out;
        }

        /* Symbolic / non-integer path: build the Horner polynomial and
         * hand it back to the evaluator. */
        Expr* base_for_poly;
        bool base_for_poly_owned = false;
        if (base_expr == NULL) {
            base_for_poly = expr_new_integer(10);
            base_for_poly_owned = true;
        } else {
            base_for_poly = base_expr;  /* aliased; fd_build_symbolic deep-copies */
        }
        Expr* out = fd_build_symbolic(first, base_for_poly);
        if (base_for_poly_owned) expr_free(base_for_poly);
        return out;
    }

    /* --- Case 3: anything else ------------------------------------- */
    /* Concrete-but-not-list-or-string (Real, Rational, Complex, plain
     * Integer) fires the ::nlst diagnostic.  Symbolic first args flow
     * through silently so downstream rewrites can still reach the call. */
    if (expr_is_numeric_like(first)) return fd_emit_nlst(first, res);
    return NULL;
}
