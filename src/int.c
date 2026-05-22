/*
 * int.c
 *
 * Integer-digit / integer-length builtins.
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
 * Sign of n is discarded.  IntegerDigits[0] -> {0}; IntegerLength[0] -> 0
 * (matching Mathematica: zero has no significant digits).  Threading over
 * list arguments is handled automatically by the evaluator (Listable).
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
#include <gmp.h>

void int_init(void) {
    symtab_add_builtin("IntegerDigits", builtin_integerdigits);
    symtab_get_def("IntegerDigits")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("IntegerLength", builtin_integerlength);
    symtab_get_def("IntegerLength")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);
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
