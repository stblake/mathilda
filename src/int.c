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
 *   DigitSum[n]                 sum of the base-10 digits of |n|.
 *   DigitSum[n, b]              sum of the base-b digits of |n|.
 *
 * Sign of n is discarded.  IntegerDigits[0] -> {0}; IntegerLength[0] -> 0;
 * DigitCount[0] -> {0,0,...,0} of length b -- zero has no significant
 * digits (matching Mathematica).  DigitSum[0] -> 0.
 * IntegerDigits / IntegerLength / DigitSum thread automatically via the
 * Listable attribute; DigitCount and FromDigits are intentionally NOT
 * Listable (only Protected): the natural first argument of each is itself
 * a list, so threading would be wrong.
 *
 * All arithmetic is done in GMP, so both machine integers and arbitrary-
 * precision bignums are supported uniformly.
 */

#include "linkedlist.h"
#include "int.h"
#include "ndarray.h"      /* is_ndarray, ndarray_delist_and_reeval */
#include "pack.h"         /* ndbuild_open_i64 */
#include "ndinteger.h"   /* NDArray buffer fast path */
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "sym_names.h"
#include "eval.h"           /* eval_and_free */
#include "common.h"         /* head_is */
#include "arithmetic.h"     /* expr_numeric_sign */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */

void int_init(void) {
    symtab_add_builtin("IntegerDigits", builtin_integerdigits);
    symtab_get_def("IntegerDigits")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("IntegerLength", builtin_integerlength);
    symtab_get_def("IntegerLength")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("IntegerExponent", builtin_integerexponent);
    symtab_get_def("IntegerExponent")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("DigitCount", builtin_digitcount);
    symtab_get_def("DigitCount")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FromDigits", builtin_fromdigits);
    symtab_get_def("FromDigits")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("IntegerString", builtin_integerstring);
    symtab_get_def("IntegerString")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("DigitSum", builtin_digitsum);
    symtab_get_def("DigitSum")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE | ATTR_NUMERICFUNCTION);
}

/* Build an Expr from an mpz_t digit, demoting to EXPR_INTEGER if it fits. */
static Expr* expr_from_mpz_digit(const mpz_t d) {
    int64_t v;
    if (expr_mpz_get_int64(d, &v)) {
        return expr_new_integer(v);
    }
    return expr_new_bigint_from_mpz(d);
}

/* Print `IntegerDigits::argb: IntegerDigits called with N argument(s);
 * between 1 and 3 arguments are expected.` to stderr.  Return NULL so the
 * evaluator leaves the call unevaluated, matching Mathematica's surface
 * behaviour. */
static Expr* int_emit_argb(size_t argc) {
    mth_message("IntegerDigits", "argb", "IntegerDigits called with %zu argument%s; between 1 and 3 arguments are expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Print `IntegerDigits::int: Integer expected at position <pos> in <call>.`
 * to stderr.  `res` is the original call expression so the user sees the
 * exact form they typed echoed back. */
static Expr* int_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    mth_message("IntegerDigits", "int", "Integer expected at position %zu in %s.", pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* IntegerDigits[n, MixedRadix[{b1, ..., bk}]] and the length form: the
 * digits of |n| in the mixed radix, least significant against bk.  Each
 * radix is an integer > 1.  What is left after b1 is one more, unbounded,
 * leading digit (IntegerDigits[10^6, MixedRadix[{24, 60, 60}]] is
 * {11, 13, 46, 40}); leading zeros are dropped, and a length pads with zeros
 * or keeps the low digits, as for a plain base.  NULL (unevaluated) for a
 * radix list that is not all integers > 1. */
static Expr* int_digits_mixed_radix(Expr* res, Expr* n_expr) {
    size_t argc = res->data.function.arg_count;
    Expr* mr = res->data.function.args[1];
    if (mr->data.function.arg_count != 1 || !head_is(mr->data.function.args[0], SYM_List))
        return NULL;
    Expr* radices = mr->data.function.args[0];
    size_t k = radices->data.function.arg_count;
    for (size_t i = 0; i < k; i++) {
        Expr* b = radices->data.function.args[i];
        if (!expr_is_integer_like(b) || expr_numeric_sign(b) <= 0 ||
            (b->type == EXPR_INTEGER && b->data.integer < 2))
            return NULL;
    }
    size_t target_len = 0;
    bool has_len = false;
    if (argc >= 3) {
        Expr* l = res->data.function.args[2];
        if (l->type != EXPR_INTEGER || l->data.integer < 0) return NULL;
        target_len = (size_t)l->data.integer;
        has_len = true;
    }
    mpz_t n, b, d;
    mpz_inits(n, b, d, NULL);
    { mpz_t t; expr_to_mpz(n_expr, t); mpz_abs(n, t); mpz_clear(t); }
    Expr** dig = malloc(sizeof(Expr*) * (k + 1));
    size_t nd = 0;                              /* digits, least significant first */
    for (size_t i = k; i-- > 0 && mpz_sgn(n) != 0;) {
        mpz_t bt;
        expr_to_mpz(radices->data.function.args[i], bt);
        mpz_fdiv_qr(n, d, n, bt);
        mpz_clear(bt);
        dig[nd++] = expr_bigint_normalize(expr_new_bigint_from_mpz(d));
    }
    if (mpz_sgn(n) != 0) dig[nd++] = expr_bigint_normalize(expr_new_bigint_from_mpz(n));
    if (nd == 0) dig[nd++] = expr_new_integer(0);   /* IntegerDigits[0, ...] is {0} */
    mpz_clears(n, b, d, NULL);
    size_t out_n = has_len ? target_len : nd;
    Expr** out = malloc(sizeof(Expr*) * (out_n ? out_n : 1));
    for (size_t j = 0; j < out_n; j++) {
        size_t from_low = out_n - 1 - j;            /* position from the least significant */
        out[j] = from_low < nd ? expr_copy(dig[from_low]) : expr_new_integer(0);
    }
    for (size_t j = 0; j < nd; j++) expr_free(dig[j]);
    free(dig);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), out, out_n);
    free(out);
    return r;
}

Expr* builtin_integerdigits(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* Base-10 digits of a rank-1 int64 buffer, read element by element out of
     * the buffer instead of materialising it first. A base argument, a higher
     * rank or a non-integer dtype degrades explicitly -- expr_is_integer_like
     * below is false for an NDArray. */
    if (res->data.function.arg_count >= 1 &&
        is_ndarray(res->data.function.args[0])) {
        Expr* nd = ndint_integerdigits(res);
        return nd ? nd : ndarray_delist_and_reeval(res);
    }
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

    /* --- IntegerDigits[n, MixedRadix[{b1, ..., bk}]] (, len) --------- */
    if (argc >= 2 && head_is(res->data.function.args[1], SYM_MixedRadix))
        return int_digits_mixed_radix(res, n_expr);

    /* --- base argument --------------------------------------------- */
    mpz_t base;
    if (argc >= 2) {
        Expr* b_expr = res->data.function.args[1];
        bool ok_base = expr_is_integer_like(b_expr);
        if (ok_base) {
            expr_to_mpz(b_expr, base);
            if (mpz_cmp_ui(base, 2) < 0) { mpz_clear(base); ok_base = false; }
        }
        if (!ok_base) {
            /* IntegerDigits::ibase, for any number that is
             * not an integer > 1; a symbolic base stays unevaluated quietly. */
            if (expr_is_numeric_like(b_expr)) {
                char* bs = expr_to_string(b_expr);
                mth_message("IntegerDigits", "ibase", "Base %s is not an integer greater than 1.",
                            bs ? bs : "?");
                free(bs);
            }
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
        bool ok_len = false;
        if (expr_is_integer_like(l_expr)) {
            mpz_t l;
            expr_to_mpz(l_expr, l);
            ok_len = mpz_sgn(l) >= 0 && mpz_fits_ulong_p(l);
            if (ok_len) target_len = (size_t)mpz_get_ui(l);
            mpz_clear(l);
        }
        if (!ok_len) {
            mpz_clear(base);
            if (expr_is_numeric_like(l_expr)) {
                char* call = expr_to_string(res);
                mth_message("IntegerDigits", "intnm",
                            "Non-negative machine-sized integer expected at position 3 in %s.",
                            call ? call : "?");
                free(call);
            }
            return NULL;
        }
        has_len = true;
    }

    /* --- machine fast path: int64 n, base up to 2^31 ----------------- */
    if (n_expr->type == EXPR_INTEGER && mpz_cmp_ui(base, 2147483647UL) <= 0 &&
        n_expr->data.integer != INT64_MIN) {
        uint64_t b = mpz_get_ui(base);
        uint64_t v = (uint64_t)(n_expr->data.integer < 0 ? -n_expr->data.integer : n_expr->data.integer);
        int64_t dg[64];
        size_t cnt = 0;
        do { dg[cnt++] = (int64_t)(v % b); v /= b; } while (v);
        size_t outn = has_len ? target_len : cnt;
        size_t use = (has_len && target_len < cnt) ? target_len : cnt;
        Expr** out_args = outn ? malloc(sizeof(Expr*) * outn) : NULL;
        size_t lead = outn - use;
        for (size_t i = 0; i < lead; i++) out_args[i] = expr_new_integer(0);
        for (size_t i = 0; i < use; i++) out_args[lead + i] = expr_new_integer(dg[use - 1 - i]);
        mpz_clear(base);
        Expr* r = expr_new_function(expr_new_symbol(SYM_List), out_args, outn);
        free(out_args);
        return r;
    }

    /* --- bignum, base 2..62: GMP's subquadratic radix conversion ------- */
    if (mpz_cmp_ui(base, 62) <= 0) {
        mpz_t n;
        expr_to_mpz(n_expr, n);
        mpz_abs(n, n);
        int b = (int)mpz_get_ui(base);
        char* str = mpz_get_str(NULL, b, n);
        size_t cnt = strlen(str);
        size_t outn = has_len ? target_len : cnt;
        size_t use = (has_len && target_len < cnt) ? target_len : cnt;
        size_t lead = outn - use;
        const char* src = str + (cnt - use);          /* the `use` least-significant */
        /* The answer is a flat list of machine integers, so it is written
         * straight into a packed int64 buffer (the canonical form):
         * 95k boxed digits cost 2.5x the radix conversion itself. NULL means
         * packing is off or the list is under threshold -> boxed digits. */
        int64_t* pbuf = NULL;
        Expr* packed = ndbuild_open_i64((int64_t)outn, &pbuf);
        Expr** out_args = (!packed && outn) ? malloc(sizeof(Expr*) * outn) : NULL;
        for (size_t i = 0; i < outn; i++) {
            int d = 0;
            if (i >= lead) {
                int c = (unsigned char)src[i - lead];
                d = (c >= '0' && c <= '9') ? c - '0'
                  : (c >= 'a' && c <= 'z') ? (b <= 36 ? 10 + c - 'a' : 36 + c - 'a')
                  : 10 + c - 'A';                            /* GMP: base > 36 uses A-Z then a-z */
            }
            if (packed) pbuf[i] = d;
            else out_args[i] = expr_new_integer(d);
        }
        void (*gmp_free)(void*, size_t);
        mp_get_memory_functions(NULL, NULL, &gmp_free);
        gmp_free(str, cnt + 1);
        mpz_clear(n);
        mpz_clear(base);
        if (packed) return packed;
        Expr* r = expr_new_function(expr_new_symbol(SYM_List), out_args, outn);
        free(out_args);
        return r;
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

    Expr* list_head = expr_new_symbol(SYM_List);
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
    mth_message("IntegerLength", "argt", "IntegerLength called with %zu argument%s; 1 or 2 arguments are expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Print `IntegerLength::int: Integer expected at position <pos> in <call>.`
 * matching Mathematica's surface diagnostic. */
static Expr* intlen_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    mth_message("IntegerLength", "int", "Integer expected at position %zu in %s.", pos, call_str ? call_str : "?");
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
            mth_message("IntegerLength", "ibase", "Base argument must be an integer >= 2.");
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
 * IntegerExponent
 *
 *   IntegerExponent[n]     == IntegerExponent[n, 10]
 *   IntegerExponent[n, b]  highest k such that b^k divides n.
 *
 * Equivalent to the number of trailing zero digits of |n| in base b.
 * Sign of n is discarded.  IntegerExponent[0, b] == Infinity (every
 * power of b divides 0, matching Mathematica's convention).
 * ===================================================================*/

/* `IntegerExponent::argt: IntegerExponent called with N argument(s);
 * 1 or 2 arguments are expected.` */
static Expr* intexp_emit_argt(size_t argc) {
    mth_message("IntegerExponent", "argt", "IntegerExponent called with %zu argument%s; 1 or 2 arguments are expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* `IntegerExponent::int: Integer expected at position <pos> in <call>.` */
static Expr* intexp_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    mth_message("IntegerExponent", "int", "Integer expected at position %zu in %s.", pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* Return the largest k such that base^k divides |n|.  Pre: n != 0,
 * base >= 2.  Both `n` and `base` are read-only and not cleared.
 *
 * Fast path for base == 2: `mpz_scan1(n, 0)` returns the position of
 * the lowest 1-bit, which equals the 2-adic valuation of |n|.  GMP's
 * scan1 acts on the twos-complement representation but treats trailing
 * zeros identically for ±n, so we can pass `n` directly without taking
 * the absolute value.
 *
 * General path: GMP's `mpz_remove(q, |n|, base)` divides |n| by base
 * as many times as possible into q and returns that count -- one
 * library call, O((log n) * M(d log b)) bit-ops. */
static size_t intexp_count(const mpz_t n, const mpz_t base) {
    if (mpz_cmp_ui(base, 2) == 0) {
        return (size_t)mpz_scan1(n, 0);
    }
    mpz_t abs_n, q;
    mpz_init(q);
    mpz_init(abs_n);
    mpz_abs(abs_n, n);
    mp_bitcnt_t cnt = mpz_remove(q, abs_n, base);
    mpz_clear(q);
    mpz_clear(abs_n);
    return (size_t)cnt;
}

Expr* builtin_integerexponent(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return intexp_emit_argt(argc);

    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Symbolic n flows through silently; concrete-but-non-integer
         * (Real, Rational, Complex) gets the ::int diagnostic. */
        if (expr_is_numeric_like(n_expr)) return intexp_emit_int(1, res);
        return NULL;
    }

    mpz_t base;
    if (argc == 2) {
        Expr* b_expr = res->data.function.args[1];
        bool ok_base = expr_is_integer_like(b_expr);
        if (ok_base) {
            expr_to_mpz(b_expr, base);
            if (mpz_cmp_ui(base, 2) < 0) { mpz_clear(base); ok_base = false; }
        }
        if (!ok_base) {
            /* IntegerExponent::ibase for any number that is not an integer
             * > 1; a symbolic base stays unevaluated quietly. */
            if (expr_is_numeric_like(b_expr)) {
                char* bs = expr_to_string(b_expr);
                mth_message("IntegerExponent", "ibase", "Base %s is not an integer greater than 1.",
                            bs ? bs : "?");
                free(bs);
            }
            return NULL;
        }
    } else {
        mpz_init_set_ui(base, 10);
    }

    mpz_t n;
    expr_to_mpz(n_expr, n);

    /* IntegerExponent[0, b] -> Infinity (every power of b divides 0). */
    if (mpz_sgn(n) == 0) {
        mpz_clear(n);
        mpz_clear(base);
        return expr_new_symbol(SYM_Infinity);
    }

    size_t k = intexp_count(n, base);
    mpz_clear(n);
    mpz_clear(base);

    /* `k` is bounded by IntegerLength[n, base], which fits in size_t
     * and -- for any plausible bignum -- comfortably within int64. */
    return expr_new_integer((int64_t)k);
}

/* =====================================================================
 * DigitCount
 * ===================================================================*/

/* Print `DigitCount::argb: DigitCount called with N argument(s);
 * between 1 and 3 arguments are expected.` matching Mathematica's
 * surface diagnostic exactly. */
static Expr* dc_emit_argb(size_t argc) {
    mth_message("DigitCount", "argb", "DigitCount called with %zu argument%s; between 1 and 3 arguments are expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Print `DigitCount::int: Integer expected at position <pos> in <call>.`
 * Used for non-integer numeric first arg.  Pure symbolic inputs flow
 * through silently (returning NULL) without touching this path. */
static Expr* dc_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    mth_message("DigitCount", "int", "Integer expected at position %zu in %s.", pos, call_str ? call_str : "?");
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
    mth_message("DigitCount", "base", "The base %s at position 2 of %s should be an integer greater than 1.", b_str ? b_str : "?", call_str ? call_str : "?");
    free(b_str);
    free(call_str);
    return NULL;
}

/* Print `DigitCount::digit: The digit <d> at position 3 of <call>
 * should be a non-negative integer less than the base.`  Used for
 * non-integer digit OR integer digit out of [0, base). */
static Expr* dc_emit_digit(Expr* d_expr, Expr* res) {
    (void)d_expr;
    char* call_str = expr_to_string(res);
    char* b_str = NULL;
    if (res->data.function.arg_count >= 2) b_str = expr_to_string(res->data.function.args[1]);
    mth_message("DigitCount", "digp", "Positive integer or list of non-negative integers less than the base %s expected at position 3 in %s.",
                b_str ? b_str : "10", call_str ? call_str : "?");
    free(b_str);
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
    /* 0 is the single digit 0 -- except in base 2, which counts
     * no digits at all (DigitCount[0] is {0, ..., 0, 1}, DigitCount[0, 2]
     * is {0, 0}). */
    if (mpz_sgn(n) == 0) { if (b_val != 2) counts[0] = 1; return; }
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
    if (mpz_sgn(n) == 0) {                     /* see dc_build_histogram */
        if (mpz_sgn(target) == 0 && mpz_cmp_ui(base, 2) != 0) mpz_set_ui(out_count, 1);
        return;
    }

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

    Expr** dargs = res->data.function.args;

    /* A LinkedList of numbers threads as the List of them does (below): the
     * same ::int check on every element, then one count per element. */
    if (dargs[0]->type == EXPR_LINKEDLIST) {
        const LinkedList* nl = dargs[0]->data.llist;
        bool numeric = true;
        LLKind nk = ll_kind(nl);
        if (nk == LL_BOXED) {               /* typed slots are numbers already */
            for (int64_t i = 0, n = ll_length(nl); i < n && numeric; i++) {
                const Expr* e = ll_peek_boxed(nl, i);
                numeric = expr_is_integer_like(e) || e->type == EXPR_REAL
#ifdef USE_MPFR
                    || e->type == EXPR_MPFR
#endif
                    ;
            }
        } else if (nk == LL_CPLX) {
            numeric = false;
        }
        if (!numeric) return dc_emit_int(1, res);
        return ll_map_arg(res, 0);
    }
    if (argc >= 2 && dargs[1]->type == EXPR_LINKEDLIST) return dc_emit_base(dargs[1], res);
    if (argc == 3 && dargs[2]->type == EXPR_LINKEDLIST) return ll_map_arg(res, 2);

    /* DigitCount[{n1, n2, ...}, ...] threads over the first argument only
     * (it is not Listable): one count per
     * number, integer or real.  A non-numeric element is ::int for the
     * whole call; a list as the base is ::base. */
    if (dargs[0]->type == EXPR_FUNCTION &&
        dargs[0]->data.function.head->type == EXPR_SYMBOL &&
        dargs[0]->data.function.head->data.symbol.name == SYM_List) {
        size_t nn = dargs[0]->data.function.arg_count;
        for (size_t i = 0; i < nn; i++) {
            Expr* e = dargs[0]->data.function.args[i];
            if (!expr_is_integer_like(e) && e->type != EXPR_REAL
#ifdef USE_MPFR
                && e->type != EXPR_MPFR
#endif
                ) return dc_emit_int(1, res);
        }
        Expr** items = malloc((nn ? nn : 1) * sizeof(Expr*));
        for (size_t i = 0; i < nn; i++) {
            Expr* ca[3];
            ca[0] = expr_copy(dargs[0]->data.function.args[i]);
            for (size_t j = 1; j < argc; j++) ca[j] = expr_copy(dargs[j]);
            items[i] = expr_new_function(expr_new_symbol("DigitCount"), ca, argc);
        }
        Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_List), items, nn));
        free(items);
        return out;
    }
    if (argc >= 2 && dargs[1]->type == EXPR_FUNCTION &&
        dargs[1]->data.function.head->type == EXPR_SYMBOL &&
        dargs[1]->data.function.head->data.symbol.name == SYM_List)
        return dc_emit_base(dargs[1], res);

    /* DigitCount[n, b, {d1, d2, ...}]: a list of counts. */
    if (argc == 3 && dargs[2]->type == EXPR_FUNCTION &&
        dargs[2]->data.function.head->type == EXPR_SYMBOL &&
        dargs[2]->data.function.head->data.symbol.name == SYM_List) {
        size_t nd = dargs[2]->data.function.arg_count;
        Expr** items = malloc((nd ? nd : 1) * sizeof(Expr*));
        for (size_t i = 0; i < nd; i++) {
            Expr* ca[3] = { expr_copy(dargs[0]), expr_copy(dargs[1]),
                            expr_copy(dargs[2]->data.function.args[i]) };
            items[i] = expr_new_function(expr_new_symbol("DigitCount"), ca, 3);
        }
        Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_List), items, nd));
        free(items);
        return out;
    }

    /* A real n counts the digits RealDigits[n, b] gives (DigitCount[2.5, 2]
     * is {2, 51}). */
    if (dargs[0]->type == EXPR_REAL
#ifdef USE_MPFR
        || dargs[0]->type == EXPR_MPFR
#endif
        ) {
        Expr* b_e = argc >= 2 ? dargs[1] : NULL;
        int64_t b = 10;
        if (b_e) {
            if (b_e->type != EXPR_INTEGER || b_e->data.integer < 2)
                return expr_is_numeric_like(b_e) ? dc_emit_base(b_e, res) : NULL;
            b = b_e->data.integer;
        }
        if ((uint64_t)b > DC_MAX_BASE_FOR_LIST) return NULL;
        int64_t target = -1;
        if (argc == 3) {
            Expr* d = dargs[2];
            if (d->type != EXPR_INTEGER || d->data.integer < 0 || d->data.integer >= b)
                return expr_is_numeric_like(d) ? dc_emit_digit(d, res) : NULL;
            target = d->data.integer;
        }
        Expr* rd_args[2] = { expr_copy(dargs[0]), expr_new_integer(b) };
        Expr* rd = eval_and_free(expr_new_function(expr_new_symbol("RealDigits"), rd_args, 2));
        Expr* digs = (rd->type == EXPR_FUNCTION && rd->data.function.arg_count == 2)
                         ? rd->data.function.args[0] : NULL;
        if (!digs || digs->type != EXPR_FUNCTION) { expr_free(rd); return NULL; }
        int64_t* counts = calloc((size_t)b, sizeof(int64_t));
        for (size_t i = 0; i < digs->data.function.arg_count; i++) {
            Expr* d = digs->data.function.args[i];
            if (d->type == EXPR_INTEGER && d->data.integer >= 0 && d->data.integer < b)
                counts[d->data.integer]++;
        }
        expr_free(rd);
        Expr* out;
        if (target >= 0) {
            out = expr_new_integer(counts[target]);
        } else {
            Expr** oa = malloc((size_t)b * sizeof(Expr*));
            for (int64_t i = 1; i < b; i++) oa[i - 1] = expr_new_integer(counts[i]);
            oa[b - 1] = expr_new_integer(counts[0]);
            out = expr_new_function(expr_new_symbol(SYM_List), oa, (size_t)b);
            free(oa);
        }
        free(counts);
        return out;
    }

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
        int64_t cv;
        if (expr_mpz_get_int64(count, &cv)) {
            result = expr_new_integer(cv);
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
        mth_message("DigitCount", "ovfl", "Base too large for the list-returning form in %s; use DigitCount[n, b, d] for a single digit instead.", call_str ? call_str : "?");
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

    Expr* list_head = expr_new_symbol(SYM_List);
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
    /* The expected tag for "1 or 2 arguments" is argt. */
    mth_message("FromDigits", "argt", "FromDigits called with %zu arguments; 1 or 2 arguments are expected.", argc);
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
    mth_message("FromDigits", "nlst", "Argument %s at position 1 in %s should be a list or a string.", arg_str ? arg_str : "?", call_str ? call_str : "?");
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
    mth_message("FromDigits", "char", "Invalid digit character %s in string at position %zu in %s.", buf, pos, call_str ? call_str : "?");
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
/* Digit i of a List (its argument, borrowed) or a LinkedList (its slot, held
 * in *hold until the next read; the caller frees *hold): digits are read one
 * at a time, never gathered. */
static Expr* fd_at(const Expr* seq, size_t i, Expr** hold) {
    if (seq->type != EXPR_LINKEDLIST) return seq->data.function.args[i];
    expr_free(*hold);
    *hold = ll_get(seq->data.llist, (int64_t)i);
    return *hold;
}

static bool fd_all_integer_digits(const Expr* seq, size_t n) {
    Expr* hold = NULL;
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) ok = expr_is_integer_like(fd_at(seq, i, &hold));
    expr_free(hold);
    return ok;
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
static Expr* fd_compute_integer_list(const Expr* seq, size_t n, const mpz_t base) {
    Expr* hold = NULL;
    mpz_t acc;
    mpz_init_set_ui(acc, 0);

    /* Every digit in [0, base) with 2 <= base <= 62: GMP's subquadratic radix
     * conversion from a digit string, instead of an O(n^2) Horner loop --
     * FromDigits[IntegerDigits[3^100000]] in a few ms rather than 0.27 s. */
    bool simple = mpz_cmp_ui(base, 2) >= 0 && mpz_cmp_ui(base, 62) <= 0 && n > 0;
    unsigned long b = simple ? mpz_get_ui(base) : 0;
    for (size_t i = 0; i < n && simple; i++) {
        Expr* d = fd_at(seq, i, &hold);
        simple = d->type == EXPR_INTEGER && d->data.integer >= 0 &&
                 (uint64_t)d->data.integer < b;
    }
    if (simple) {
        char* str = malloc(n + 1);
        for (size_t i = 0; i < n; i++) {
            int64_t d = fd_at(seq, i, &hold)->data.integer;
            str[i] = (char)(d < 10 ? '0' + d : (b <= 36 ? 'a' + (d - 10)
                                                : (d < 36 ? 'A' + (d - 10) : 'a' + (d - 36))));
        }
        str[n] = 0;
        mpz_set_str(acc, str, (int)b);
        free(str);
    } else {
        mpz_t tmp;
        mpz_init(tmp);
        for (size_t i = 0; i < n; i++) {
            mpz_t digit;
            expr_to_mpz(fd_at(seq, i, &hold), digit);
            mpz_mul(tmp, acc, base);
            mpz_add(acc, tmp, digit);
            mpz_clear(digit);
        }
        mpz_clear(tmp);
    }
    expr_free(hold);
    Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(acc));
    mpz_clear(acc);
    return out;
}

/* Compute the integer FromDigits result for a string of characters in
 * an integer base.  Returns a fresh Expr on success; on invalid char
 * emits the `::char` diagnostic and returns NULL with `res` left
 * unevaluated.  Empty string -> 0. */
static Expr* fd_compute_integer_string(const char* s, const mpz_t base,
                                        Expr* res) {
    /* Every character a digit below the base (base <= 36, case-insensitive):
     * GMP's subquadratic conversion. Otherwise ("12a" in base 10 is 130)
     * the character values go through Horner. */
    if (mpz_cmp_ui(base, 36) <= 0 && *s) {
        unsigned long b = mpz_get_ui(base);
        bool simple = true;
        for (const char* p = s; *p && simple; p++) {
            int d = fd_char_to_digit((unsigned char)*p);
            simple = d >= 0 && (unsigned long)d < b;
        }
        if (simple) {
            mpz_t v;
            mpz_init_set_str(v, s, (int)b);
            Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(v));
            mpz_clear(v);
            return out;
        }
    }
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
    int64_t av;
    if (expr_mpz_get_int64(acc, &av)) {
        out = expr_new_integer(av);
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
static Expr* fd_build_symbolic(const Expr* seq, size_t n, Expr* base_expr) {
    Expr* hold = NULL;
    if (n == 0) return expr_new_integer(0);

    Expr** plus_args = malloc(sizeof(Expr*) * n);

    for (size_t i = 0; i < n; i++) {
        Expr* digit = expr_copy(fd_at(seq, i, &hold));
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
        Expr* power_node = expr_new_function(expr_new_symbol(SYM_Power),
                                              power_args, 2);
        free(power_args);

        /* Times[digit, Power[base, exponent]]. */
        Expr** times_args = malloc(sizeof(Expr*) * 2);
        times_args[0] = digit;
        times_args[1] = power_node;
        Expr* times_node = expr_new_function(expr_new_symbol(SYM_Times),
                                              times_args, 2);
        free(times_args);
        plus_args[i] = times_node;
    }

    /* Single-element list collapses to Plus[x] which the evaluator
     * simplifies to x via the OneIdentity attribute -- still correct
     * even though we could short-circuit it here. */
    Expr* plus_node = expr_new_function(expr_new_symbol(SYM_Plus),
                                         plus_args, n);
    free(plus_args);
    expr_free(hold);
    return plus_node;
}

/* Nested Horner form b (b (d0) + d1) + d2 ..., the form
 * given for symbolic digits in a numeric base: FromDigits[{a, b, c}] is
 * 10 (10 a + b) + c. */
static Expr* fd_build_horner(const Expr* seq, size_t n, Expr* base_expr) {
    Expr* hold = NULL;
    if (n == 0) return expr_new_integer(0);
    Expr* acc = expr_copy(fd_at(seq, 0, &hold));
    for (size_t i = 1; i < n; i++) {
        Expr* t_args[2] = { expr_copy(base_expr), acc };
        Expr* t = expr_new_function(expr_new_symbol(SYM_Times), t_args, 2);
        Expr* p_args[2] = { t, expr_copy(fd_at(seq, i, &hold)) };
        acc = expr_new_function(expr_new_symbol(SYM_Plus), p_args, 2);
    }
    expr_free(hold);
    return acc;
}

static bool fd_is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* FromDigits[{{d1, ..., dm}, e}, b] and FromDigits[{{d1, ..., dm, {r1, ...}}, e}, b]:
 * the digits d1 d2 ... with the "point" placed so the value is
 * sum_i d_i b^(e - i), a trailing sublist repeating forever (the
 * RealDigits form): (FromDigits[d] + R / (b^p - 1)) b^(e - m) with R the
 * repeating block's value and p its length. Built as an expression and
 * evaluated, so it is exact for integer b and e and symbolic otherwise. */
static Expr* fd_digits(const Expr* seq, size_t fn, Expr* base_expr);

/* `digits` is a List or a LinkedList (its digits read one at a time). */
static Expr* fd_mantissa(Expr* digits, Expr* e, Expr* base_expr) {
    size_t n = (size_t)ll_seq_len(digits);
    Expr* hold = NULL;
    Expr* last = n > 0 ? fd_at(digits, n - 1, &hold) : NULL;
    Expr* rec = (last && ll_is_seq(last)) ? expr_copy(last) : NULL;
    expr_free(hold);
    size_t m = rec ? n - 1 : n;
    Expr* mant = fd_digits(digits, m, base_expr);
    if (!mant) { expr_free(rec); return NULL; }
    if (rec) {
        /* + FromDigits[rec, b] / (b^p - 1) */
        size_t p = rec->type == EXPR_LINKEDLIST ? (size_t)ll_length(rec->data.llist)
                                                : rec->data.function.arg_count;
        Expr* r_args[2] = { rec, expr_copy(base_expr) };
        Expr* r = expr_new_function(expr_new_symbol("FromDigits"), r_args, 2);
        Expr* bp_args[2] = { expr_copy(base_expr), expr_new_integer((int64_t)p) };
        Expr* bp = expr_new_function(expr_new_symbol(SYM_Power), bp_args, 2);
        Expr* den_args[2] = { bp, expr_new_integer(-1) };
        Expr* den = expr_new_function(expr_new_symbol(SYM_Plus), den_args, 2);
        Expr* inv_args[2] = { den, expr_new_integer(-1) };
        Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), inv_args, 2);
        Expr* q_args[2] = { r, inv };
        Expr* q = expr_new_function(expr_new_symbol(SYM_Times), q_args, 2);
        Expr* s_args[2] = { mant, q };
        mant = expr_new_function(expr_new_symbol(SYM_Plus), s_args, 2);
    }
    /* * b^(e - m) */
    Expr* ex_args[2] = { expr_copy(e), expr_new_integer(-(int64_t)m) };
    Expr* ex = expr_new_function(expr_new_symbol(SYM_Plus), ex_args, 2);
    Expr* pw_args[2] = { expr_copy(base_expr), ex };
    Expr* pw = expr_new_function(expr_new_symbol(SYM_Power), pw_args, 2);
    Expr* out_args[2] = { mant, pw };
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), out_args, 2));
}

/* FromDigits[s, "Roman"]: a Roman numeral, case-insensitive, read additively
 * with the subtractive pairs (a smaller numeral before a larger one is
 * subtracted) -- "XLII" is 42, "IIII" is 4.
 * Anything else draws FromRomanNumeral::nrom and gives $Failed. */
/* The Roman numeral system, shared by FromDigits[s, "Roman"] (reading) and
 * IntegerString[n, "Roman"] (writing): the seven numerals and, for writing,
 * the six subtractive pairs, largest first. */
static const struct { int64_t value; const char* text; } roman_table[] = {
    {1000, "M"}, {900, "CM"}, {500, "D"}, {400, "CD"}, {100, "C"}, {90, "XC"},
    {50, "L"}, {40, "XL"}, {10, "X"}, {9, "IX"}, {5, "V"}, {4, "IV"}, {1, "I"}
};
#define ROMAN_TABLE_LEN (sizeof(roman_table) / sizeof(roman_table[0]))

/* Value of one numeral character, case-insensitive; 0 if it is not one. */
static int64_t roman_numeral_value(char ch) {
    int c = ch >= 'a' && ch <= 'z' ? ch - 32 : ch;
    for (size_t k = 0; k < ROMAN_TABLE_LEN; k++)
        if (roman_table[k].text[1] == '\0' && roman_table[k].text[0] == c)
            return roman_table[k].value;
    return 0;
}

static Expr* fd_roman(const char* s) {
    int64_t total = 0;
    size_t len = strlen(s);
    bool ok = len > 0;
    for (size_t i = 0; i < len && ok; i++) {
        int64_t v = roman_numeral_value(s[i]);
        if (!v) { ok = false; break; }
        int64_t nx = (i + 1 < len) ? roman_numeral_value(s[i + 1]) : 0;
        total += (nx > v) ? -v : v;
    }
    if (!ok) {
        mth_message("FromRomanNumeral", "nrom",
                    "String \"%s\" does not represent a valid roman numeral.", s);
        return expr_new_symbol("$Failed");
    }
    return expr_new_integer(total);
}

/* FromDigits of the digit sequence fv[0..fn) -- a List's elements, or a
 * LinkedList's read from its slots -- in base base_expr (NULL: 10). */
static Expr* fd_digits(const Expr* seq, size_t fn, Expr* base_expr) {
    /* the first two elements, for the mantissa forms (new references) */
    Expr* f0 = NULL;
    Expr* f1 = NULL;
    if (fn == 2) {
        bool o;
        Expr* v = ll_seq_view(seq, 0, &o); f0 = o ? v : expr_copy(v);
        v = ll_seq_view(seq, 1, &o);       f1 = o ? v : expr_copy(v);
    }

    Expr* base_for_poly = base_expr ? expr_copy(base_expr) : expr_new_integer(10);

    /* {{digits...}, exponent}: the mantissa-exponent (RealDigits) form. */
    if (fn == 2 && ll_is_seq(f0) &&
        !ll_is_seq(f1)) {
        Expr* out = fd_mantissa(f0, f1,
                                base_for_poly);
        expr_free(base_for_poly);
        expr_free(f0); expr_free(f1);
        return out;
    }
    /* {{digits...}, {e1, e2, ...}} threads over the exponents:
     * FromDigits[{{1, 2}, {3, 4}}] is {120, 1200}. */
    if (fn == 2 && ll_is_seq(f0) &&
        ll_is_seq(f1)) {
        size_t ne = (size_t)ll_seq_len(f1);
        Expr* eh = NULL;
        bool flat = true;
        for (size_t i = 0; i < ne && flat; i++) flat = !ll_is_seq(fd_at(f1, i, &eh));
        if (flat) {
            Expr** items = malloc((ne ? ne : 1) * sizeof(Expr*));
            size_t got = 0;
            for (; got < ne; got++) {
                items[got] = fd_mantissa(f0, fd_at(f1, got, &eh), base_for_poly);
                if (!items[got]) break;
            }
            if (got < ne) {
                for (size_t i = 0; i < got; i++) expr_free(items[i]);
                free(items);
                expr_free(base_for_poly);
                expr_free(eh); expr_free(f0); expr_free(f1);
                return NULL;
            }
            /* a LinkedList of exponents answers with a LinkedList, built directly */
            Expr* out;
            if (f1->type == EXPR_LINKEDLIST) {
                LinkedList* ll = ll_new();
                for (size_t i = 0; i < ne; i++) ll_push_back(ll, items[i]);
                out = expr_new_linkedlist(ll);
            } else {
                out = expr_new_function(expr_new_symbol(SYM_List), items, ne);
            }
            free(items);
            expr_free(base_for_poly);
            expr_free(eh); expr_free(f0); expr_free(f1);
            return out;
        }
        expr_free(eh);
    }

    /* Integer digits in an integer base (any base: 0, 1, -10 are all just
     * Horner): exact in GMP. */
    if (fd_all_integer_digits(seq, fn) && expr_is_integer_like(base_for_poly)) {
        mpz_t base;
        expr_to_mpz(base_for_poly, base);
        Expr* out = fd_compute_integer_list(seq, fn, base);
        mpz_clear(base);
        expr_free(base_for_poly);
        expr_free(f0); expr_free(f1);
        return out;
    }

    /* Otherwise an expression: nested Horner in a numeric base
     * (10 (10 a + b) + c), the expanded sum in a symbolic one
     * (3 + 2 b + b^2). */
    Expr* out = expr_is_numeric_like(base_for_poly)
                    ? fd_build_horner(seq, fn, base_for_poly)
                    : fd_build_symbolic(seq, fn, base_for_poly);
    expr_free(base_for_poly);
    expr_free(f0); expr_free(f1);
    return out;
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
    if (first->type == EXPR_STRING && base_expr && base_expr->type == EXPR_STRING &&
        strcmp(base_expr->data.string, "Roman") == 0)
        return fd_roman(first->data.string);

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
                mth_message("FromDigits", "ibase", "Base argument must be an integer >= 2.");
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

    /* --- Case 1b: a packed / NDArray int64 digit buffer -------------- */
    /* FromDigits is on pack.c's AWARE list, so a packed digit list arrives
     * as its buffer instead of 10^5 materialised Integers (which cost 5x the
     * conversion itself). Digits in [0, b), b <= 62, go straight into
     * mpz_set_str; anything else is delisted and takes the List path. */
    if (is_ndarray(first)) {
        if (first->data.ndarray.rank == 1 && first->data.ndarray.dtype == NDT_INT64 &&
            (!base_expr || (base_expr->type == EXPR_INTEGER &&
                            base_expr->data.integer >= 2 && base_expr->data.integer <= 62))) {
            size_t n = ndarray_size(first);
            const void* buf = first->data.ndarray.data;
            int64_t b = base_expr ? base_expr->data.integer : 10;
            char* str = malloc(n + 1);
            size_t i = 0;
            for (; i < n; i++) {
                int64_t d = ndt_get_i(buf, i, NDT_INT64);
                if (d < 0 || d >= b) break;
                str[i] = (char)(d < 10 ? '0' + d : (b <= 36 ? 'a' + (d - 10)
                                                    : (d < 36 ? 'A' + (d - 10) : 'a' + (d - 36))));
            }
            if (i == n && n > 0) {
                str[n] = 0;
                mpz_t acc;
                mpz_init(acc);
                mpz_set_str(acc, str, (int)b);
                free(str);
                Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(acc));
                mpz_clear(acc);
                return out;
            }
            free(str);
        }
        return ndarray_delist_and_reeval(res);
    }

    /* --- Case 2: first argument is a List (or a LinkedList) -------- */
    if (fd_is_list(first) || first->type == EXPR_LINKEDLIST)
        return fd_digits(first, (size_t)ll_seq_len(first), base_expr);

    /* --- Case 3: anything else ------------------------------------- */
    /* Concrete-but-not-list-or-string (Real, Rational, Complex, plain
     * Integer) fires the ::nlst diagnostic.  Symbolic first args flow
     * through silently so downstream rewrites can still reach the call. */
    if (expr_is_numeric_like(first)) return fd_emit_nlst(first, res);
    return NULL;
}

/* =====================================================================
 * IntegerString
 *
 *   IntegerString[n]            decimal-digit string of |n|.
 *   IntegerString[n, b]         base-b digit string, 2 <= b <= 36.
 *                               Digits 10..35 use lowercase 'a'..'z'.
 *   IntegerString[n, b, len]    pads on the left with '0' to length len.
 *                               If |n| has more than len base-b digits the
 *                               result is the len least-significant digits.
 *
 * Sign of n is discarded; IntegerString[0] -> "0".  Attributes are
 * Listable | Protected so the function threads over List args in any
 * position element-wise.  Built on GMP's `mpz_get_str`, which renders
 * arbitrary-precision integers in any base 2..62 in one library call;
 * we cap the base at 36 to match Mathematica's surface contract.
 * ===================================================================*/

/* `IntegerString::argb: IntegerString called with N argument(s); between
 * 1 and 3 arguments are expected.` */
static Expr* is_emit_argb(size_t argc) {
    mth_message("IntegerString", "argb", "IntegerString called with %zu argument%s; between 1 and 3 arguments are expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* `IntegerString::int: Integer expected at position <pos> in <call>.` */
static Expr* is_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    mth_message("IntegerString", "int", "Integer expected at position %zu in %s.", pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* `IntegerString::basf: <b> is not a valid base for IntegerString;
 * the base must be an integer between 2 and 36.` */
static Expr* is_emit_basf(Expr* b_expr, Expr* res) {
    char* b_str = expr_to_string(b_expr);
    char* call_str = expr_to_string(res);
    mth_message("IntegerString", "basf", "%s is not a valid base for IntegerString in %s; the base must be an integer between 2 and 36.", b_str ? b_str : "?", call_str ? call_str : "?");
    free(b_str);
    free(call_str);
    return NULL;
}

/* `IntegerString::intnn: Non-negative machine-sized integer expected at
 * position 3.` */
static Expr* is_emit_intnn(void) {
    mth_message("IntegerString", "intnn", "Non-negative machine-sized integer expected at position 3.");
    return NULL;
}

/* Maximum allowed base for IntegerString.  Mathematica caps this at 36
 * because the surface convention only defines digit characters up to
 * 'z' (10 numerals + 26 lowercase letters).  Higher bases are
 * representable via GMP (up to 62) but produce strings users cannot
 * decode without consulting our own conventions. */
#define IS_MAX_BASE 36

/* ASCII case-insensitive strcmp (strcasecmp is POSIX, hidden under -std=c99). */
static int strcasecmp_ascii(const char* a, const char* b) {
    for (;; a++, b++) {
        int ca = (*a >= 'A' && *a <= 'Z') ? *a + 32 : *a;
        int cb = (*b >= 'A' && *b <= 'Z') ? *b + 32 : *b;
        if (ca != cb || !ca) return ca - cb;
    }
}

/* IntegerString[n, "Roman"]: |n| as a Roman numeral, greedy over the shared
 * roman_table (so 4 is "IV" and 4000 is "MMMM"); 0 is "N" (nulla). The
 * reference renders 5000 and above with overbarred (x1000) numerals, which is
 * a typeset box, not a string -- those stay unevaluated. A length argument
 * does not pad or truncate a Roman numeral. */
static Expr* is_roman(Expr* n_expr) {
    if (n_expr->type != EXPR_INTEGER) return NULL;   /* bigint: >= 5000 */
    int64_t v = n_expr->data.integer;
    if (v <= -5000 || v >= 5000) return NULL;
    if (v < 0) v = -v;                                /* |n| */
    if (v == 0) return expr_new_string("N");
    char buf[32];                                     /* "MMMMDCCCLXXXVIII" is 16 */
    size_t o = 0;
    for (size_t k = 0; k < ROMAN_TABLE_LEN; k++)
        while (v >= roman_table[k].value) {
            for (const char* t = roman_table[k].text; *t; t++) buf[o++] = *t;
            v -= roman_table[k].value;
        }
    buf[o] = '\0';
    return expr_new_string(buf);
}

Expr* builtin_integerstring(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return is_emit_argb(argc);

    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Symbolic n flows through silently so downstream rewrites can
         * still reach the call; concrete-but-non-integer (Real, Rational,
         * Complex, ...) gets the ::int diagnostic before the call is left
         * unevaluated. */
        if (expr_is_numeric_like(n_expr)) return is_emit_int(1, res);
        return NULL;
    }

    /* --- IntegerString[n, "Roman"] -------------------------------- */
    if (argc >= 2 && res->data.function.args[1]->type == EXPR_STRING &&
        strcasecmp_ascii(res->data.function.args[1]->data.string, "Roman") == 0)
        return is_roman(n_expr);

    /* --- base argument --------------------------------------------- */
    unsigned long base = 10;
    if (argc >= 2) {
        Expr* b_expr = res->data.function.args[1];
        if (!expr_is_integer_like(b_expr)) {
            /* Mathematica collapses fractional / real / complex bases
             * into ::basf, and so is a string other than "Roman"
             * (IntegerString[10, "Foo"]).  Symbolic bases stay silent. */
            if (expr_is_numeric_like(b_expr) || b_expr->type == EXPR_STRING)
                return is_emit_basf(b_expr, res);
            return NULL;
        }
        mpz_t b;
        expr_to_mpz(b_expr, b);
        if (mpz_cmp_ui(b, 2) < 0 || mpz_cmp_ui(b, IS_MAX_BASE) > 0) {
            is_emit_basf(b_expr, res);
            mpz_clear(b);
            return NULL;
        }
        base = mpz_get_ui(b);
        mpz_clear(b);
    }

    /* --- length argument ------------------------------------------- */
    bool has_len = false;
    size_t target_len = 0;
    if (argc == 3) {
        Expr* l_expr = res->data.function.args[2];
        if (!expr_is_integer_like(l_expr)) {
            if (expr_is_numeric_like(l_expr)) return is_emit_int(3, res);
            return NULL;
        }
        mpz_t l;
        expr_to_mpz(l_expr, l);
        if (mpz_sgn(l) < 0 || !mpz_fits_ulong_p(l)) {
            mpz_clear(l);
            return is_emit_intnn();
        }
        target_len = (size_t)mpz_get_ui(l);
        mpz_clear(l);
        has_len = true;
    }

    /* --- render |n| via mpz_get_str -------------------------------- */
    /* GMP's mpz_get_str renders arbitrary-precision integers in any
     * base 2..62 in one library call.  For positive bases in [2, 36] it
     * emits '0'..'9','a'..'z' (lowercase), exactly matching Mathematica's
     * IntegerString convention.  The required buffer is bounded by
     * `mpz_sizeinbase(n, base) + 2` (one byte for an optional sign,
     * one for the nul terminator); since we take |n| before rendering
     * the sign byte is unused but kept in the allocation for safety. */
    mpz_t n;
    expr_to_mpz(n_expr, n);
    mpz_abs(n, n);

    size_t need = mpz_sizeinbase(n, (int)base) + 2;
    char* digits = malloc(need);
    if (!digits) {
        mpz_clear(n);
        return NULL;
    }
    mpz_get_str(digits, (int)base, n);
    mpz_clear(n);

    size_t dlen = strlen(digits);

    /* --- pad / truncate to target_len ------------------------------- */
    char* out;
    if (has_len) {
        out = malloc(target_len + 1);
        if (!out) {
            free(digits);
            return NULL;
        }
        if (dlen >= target_len) {
            /* Keep len least-significant digits = the last target_len
             * characters.  When target_len == 0 this is the empty string,
             * matching Mathematica's IntegerString[n, b, 0] -> "". */
            if (target_len > 0) {
                memcpy(out, digits + (dlen - target_len), target_len);
            }
            out[target_len] = '\0';
        } else {
            size_t pad = target_len - dlen;
            memset(out, '0', pad);
            memcpy(out + pad, digits, dlen);
            out[target_len] = '\0';
        }
        free(digits);
    } else {
        out = digits;
    }

    Expr* result = expr_new_string(out);
    free(out);
    return result;
}

/* =====================================================================
 * DigitSum
 *
 *   DigitSum[n]      sum of the base-10 digits of |n|.
 *   DigitSum[n, b]   sum of the base-b digits of |n|, b >= 2.
 *
 * DigitSum[0] -> 0.  Sign is discarded.  Both machine integers and
 * arbitrary-precision bignums are supported via GMP.
 *
 * Fast path: ulong base  — uses mpz_tdiv_q_ui (returns remainder
 *   without a separate mpz_t for r).
 * Slow path: bignum base — uses mpz_tdiv_qr.
 * ===================================================================*/

Expr* builtin_digitsum(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) {
        mth_message("DigitSum", "argb", "DigitSum called with %zu argument%s; 1 or 2 arguments are expected.", argc, argc == 1 ? "" : "s");
        return NULL;
    }

    /* --- n -------------------------------------------------------- */
    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        if (expr_is_numeric_like(n_expr)) {
            char* call_str = expr_to_string(res);
            mth_message("DigitSum", "int", "Integer expected at position 1 in %s.", call_str ? call_str : "?");
            free(call_str);
            return NULL;
        }
        return NULL; /* symbolic — leave unevaluated silently */
    }

    /* --- base ----------------------------------------------------- */
    mpz_t base;
    int bignum_base = 0;
    if (argc == 2) {
        Expr* b_expr = res->data.function.args[1];
        if (!expr_is_integer_like(b_expr)) {
            if (expr_is_numeric_like(b_expr)) {
                char* b_str   = expr_to_string(b_expr);
                char* call_str = expr_to_string(res);
                mth_message("DigitSum", "base", "The base %s at position 2 of %s should be an integer greater than 1.", b_str ? b_str : "?", call_str ? call_str : "?");
                free(b_str);
                free(call_str);
                return NULL;
            }
            return NULL; /* symbolic base — leave unevaluated */
        }
        expr_to_mpz(b_expr, base);
        if (mpz_cmp_ui(base, 2) < 0) {
            char* b_str   = expr_to_string(b_expr);
            char* call_str = expr_to_string(res);
            mth_message("DigitSum", "base", "The base %s at position 2 of %s should be an integer greater than 1.", b_str ? b_str : "?", call_str ? call_str : "?");
            free(b_str);
            free(call_str);
            mpz_clear(base);
            return NULL;
        }
        bignum_base = !mpz_fits_ulong_p(base);
    } else {
        mpz_init_set_ui(base, 10);
    }

    /* --- compute digit sum ---------------------------------------- */
    mpz_t abs_n;
    expr_to_mpz(n_expr, abs_n);
    mpz_abs(abs_n, abs_n); /* discard sign */

    mpz_t sum, q;
    mpz_init_set_ui(sum, 0);
    mpz_init(q);

    if (!bignum_base) {
        /* Fast path: base fits in unsigned long. */
        unsigned long b_val = mpz_get_ui(base);
        while (mpz_sgn(abs_n) > 0) {
            unsigned long rem = mpz_tdiv_q_ui(q, abs_n, b_val);
            mpz_add_ui(sum, sum, rem);
            mpz_swap(abs_n, q);
        }
    } else {
        /* Slow path: bignum base — full mpz_tdiv_qr. */
        mpz_t r;
        mpz_init(r);
        while (mpz_sgn(abs_n) > 0) {
            mpz_tdiv_qr(q, r, abs_n, base);
            mpz_add(sum, sum, r);
            mpz_swap(abs_n, q);
        }
        mpz_clear(r);
    }

    mpz_clear(abs_n);
    mpz_clear(q);
    mpz_clear(base);

    Expr* result;
    int64_t sv;
    if (expr_mpz_get_int64(sum, &sv)) {
        result = expr_new_integer(sv);
    } else {
        result = expr_new_bigint_from_mpz(sum);
    }
    mpz_clear(sum);
    return result;
}
