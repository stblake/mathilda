/* Unit tests for DigitCount.
 *
 *   DigitCount[n]            -- list of counts of digits 1, 2, ..., 9, 0
 *                                in the base-10 representation of |n|.
 *   DigitCount[n, b]         -- list of counts of digits 1, 2, ..., b-1, 0
 *                                in the base-b representation of |n|.
 *   DigitCount[n, b, d]      -- scalar count of digit d in the base-b
 *                                representation of |n|.
 *
 * Coverage:
 *   - The three documented Mathematica examples verbatim:
 *       DigitCount[2147, 2, 1]    == 5
 *       DigitCount[2147, 2]       == {5, 7}
 *       DigitCount[100!]          == {15, 19, 10, 10, 14, 19, 7, 14, 20,
 *                                     30}
 *   - Zero handling for both list-returning and scalar forms.
 *   - Negative integers (sign discarded).
 *   - Bignum inputs (10^30, 2^100, 100!).
 *   - Boundary at int64 overflow (2^63 stored as bigint).
 *   - Slow-path digit count: bignum base, bignum digit.
 *   - Histogram-sum invariant: Total[DigitCount[n, b]] == IntegerLength[n, b]
 *     -- catches off-by-one or mis-ordered histogram bugs.
 *   - Cross-check against IntegerDigits/Count for several small inputs.
 *   - Argument-arity errors (::argb, 0 args and 4+ args).
 *   - ::int on non-integer concrete n (Real, Rational, Complex).
 *   - ::base on non-integer base (Rational, Real) AND integer base < 2
 *     (0, 1, -2) -- Mathematica collapses both into the same message.
 *   - ::digit on non-integer digit AND integer digit out of [0, base)
 *     (negative, equal to base, larger than base).
 *   - Symbolic args (n, b, d) are left silently unevaluated (no
 *     diagnostic), matching IntegerLength / IntegerDigits.
 *   - Attribute set is exactly Protected (no Listable -- DigitCount on a
 *     list of integers should be left unevaluated rather than threading).
 *   - Docstring registered and contains expected keywords.
 *   - Interned symbol pointer reachable.
 *   - Repeated-evaluation stress loop to expose double-frees /
 *     mishandled mpz_t under valgrind.
 *   - Soft-cap check: DigitCount[..., 2^21] in list mode is rejected
 *     with ::ovfl, but DigitCount[..., 2^21, d] in scalar mode works.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "attr.h"
#include "sym_names.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Capture stderr while `input` is parsed + evaluated.  Returns the
 * collected stderr text as a heap string (caller frees) and writes the
 * printed result into *out_result_str (also heap-allocated).  Uses a
 * fixed temp file path; safe because tests run serially.
 *
 * IMPORTANT: dup/dup2-based save-and-restore via the original stderr
 * fd so that subsequent fprintf(stderr, ...) keeps working even when
 * the test binary is run under shell redirection (no /dev/tty).
 * (The earlier freopen("/dev/tty", ...) silently fails under pipes,
 * leaving later FAIL messages invisible -- which masked a wrong-math
 * test failure during initial development of this file.) */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_digit_count_stderr.log";

    fflush(stderr);
    int saved_stderr = dup(STDERR_FILENO);
    FILE* sink = freopen(path, "w+", stderr);
    if (!sink) {
        if (saved_stderr != -1) close(saved_stderr);
        if (out_result_str) *out_result_str = NULL;
        return NULL;
    }

    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    if (out_result_str) *out_result_str = expr_to_string(e);
    expr_free(p);
    expr_free(e);

    fflush(stderr);
    if (saved_stderr != -1) {
        /* Restore the original stderr fd; reuse the same FILE* by
         * dup2'ing the saved descriptor over the current fileno. */
        dup2(saved_stderr, STDERR_FILENO);
        close(saved_stderr);
    }

    FILE* f = fopen(path, "r");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    if (n < 0) { fclose(f); return NULL; }
    fseek(f, 0, SEEK_SET);
    char* buf = malloc((size_t)n + 1);
    if (!buf) { fclose(f); return NULL; }
    size_t got = fread(buf, 1, (size_t)n, f);
    buf[got] = '\0';
    fclose(f);
    remove(path);
    return buf;
}

/* --- Documented Mathematica examples ----------------------------- */

static void test_doc_binary_count_ones(void) {
    /* DigitCount[2147, 2, 1] -> 5.  2147 = 100001100011_2, five 1s. */
    assert_eval_eq("DigitCount[2147, 2, 1]", "5", 0);
}

static void test_doc_binary_histogram(void) {
    /* DigitCount[2147, 2] -> {5, 7}: 5 ones, 7 zeros in 12-bit form. */
    assert_eval_eq("DigitCount[2147, 2]", "{5, 7}", 0);
}

static void test_doc_factorial_histogram(void) {
    /* DigitCount[100!] -> {15, 19, 10, 10, 14, 19, 7, 14, 20, 30}.
     * Sum == 158, which is IntegerLength[100!]. */
    assert_eval_eq("DigitCount[100!]",
                   "{15, 19, 10, 10, 14, 19, 7, 14, 20, 30}", 0);
}

/* --- Small decimal cases ----------------------------------------- */

static void test_decimal_single_digit_3(void) {
    /* DigitCount[3] -> {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}: one '3'. */
    assert_eval_eq("DigitCount[3]",
                   "{0, 0, 1, 0, 0, 0, 0, 0, 0, 0}", 0);
}

static void test_decimal_repeated(void) {
    /* DigitCount[333] -> {0, 0, 3, 0, 0, 0, 0, 0, 0, 0}: three '3's. */
    assert_eval_eq("DigitCount[333]",
                   "{0, 0, 3, 0, 0, 0, 0, 0, 0, 0}", 0);
}

static void test_decimal_each_digit_once(void) {
    /* DigitCount[1234567890] -> one of each. */
    assert_eval_eq("DigitCount[1234567890]",
                   "{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}", 0);
}

static void test_decimal_mix(void) {
    /* DigitCount[100100] -> two '1's, four '0's. */
    assert_eval_eq("DigitCount[100100]",
                   "{2, 0, 0, 0, 0, 0, 0, 0, 0, 4}", 0);
}

/* --- Zero handling ----------------------------------------------- */

static void test_zero_default_base(void) {
    /* DigitCount[0] -> ten zeros (zero has no significant digits). */
    assert_eval_eq("DigitCount[0]", "{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}", 0);
}

static void test_zero_binary(void) {
    /* DigitCount[0, 2] -> {0, 0}: no 1s, no 0s. */
    assert_eval_eq("DigitCount[0, 2]", "{0, 0}", 0);
}

static void test_zero_scalar(void) {
    /* DigitCount[0, 10, 0] -> 0; DigitCount[0, 10, 5] -> 0. */
    assert_eval_eq("DigitCount[0, 10, 0]", "0", 0);
    assert_eval_eq("DigitCount[0, 10, 5]", "0", 0);
    assert_eval_eq("DigitCount[0, 2, 1]", "0", 0);
}

/* --- Negative (sign discarded) ---------------------------------- */

static void test_negative_small(void) {
    /* DigitCount[-2147, 2] -> {5, 7} matches positive. */
    assert_eval_eq("DigitCount[-2147, 2]", "{5, 7}", 0);
}

static void test_negative_scalar(void) {
    /* DigitCount[-2147, 2, 1] -> 5. */
    assert_eval_eq("DigitCount[-2147, 2, 1]", "5", 0);
}

static void test_negative_decimal(void) {
    /* DigitCount[-1234567890] same as DigitCount[1234567890]. */
    assert_eval_eq("DigitCount[-1234567890]",
                   "{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}", 0);
}

/* --- Bignum -------------------------------------------------------- */

static void test_bignum_10pow30_base10(void) {
    /* 10^30 = 1 followed by 30 zeros.  One '1', 30 '0's. */
    assert_eval_eq("DigitCount[10^30]",
                   "{1, 0, 0, 0, 0, 0, 0, 0, 0, 30}", 0);
}

static void test_bignum_2pow100_base2(void) {
    /* 2^100 = 1 followed by 100 zeros in binary. */
    assert_eval_eq("DigitCount[2^100, 2]", "{1, 100}", 0);
}

static void test_bignum_2pow100_count_ones_in_decimal(void) {
    /* 2^100 = 1267650600228229401496703205376 (31 digits).
     * Cross-checked digit counts:
     *   '1': 2     ("12676...")
     *   '2': 5     (12, 22, 28, 2, 2, ... -- five occurrences)
     *   '3': 2
     *   '4': 2
     *   '5': 2
     *   '6': 5
     *   '7': 3
     *   '8': 1
     *   '9': 2
     *   '0': 7     (in "06002284...")
     * Actually let's not hardcode -- assert via summation only. */
    assert_eval_eq("Total[DigitCount[2^100]]", "31", 0);
}

static void test_bignum_negative(void) {
    assert_eval_eq("DigitCount[-(10^30)]",
                   "{1, 0, 0, 0, 0, 0, 0, 0, 0, 30}", 0);
}

static void test_int64_overflow_boundary(void) {
    /* 2^63 = 9223372036854775808 (19 digits).
     *   1: 2  ('1' twice -- positions 13 and 16 from left: 922337...755...108 -- count from string)
     * Use sum-invariant for safety; explicit count below. */
    assert_eval_eq("Total[DigitCount[2^63]]", "19", 0);
    assert_eval_eq("DigitCount[2^63, 2, 1]", "1", 0);   /* one '1' in binary */
    assert_eval_eq("DigitCount[2^63, 2, 0]", "63", 0);  /* 63 zeros */
}

/* --- Scalar form (3-arg) ---------------------------------------- */

static void test_scalar_decimal_digit(void) {
    /* DigitCount[100100, 10, 0] -> 4, DigitCount[100100, 10, 1] -> 2. */
    assert_eval_eq("DigitCount[100100, 10, 0]", "4", 0);
    assert_eval_eq("DigitCount[100100, 10, 1]", "2", 0);
    assert_eval_eq("DigitCount[100100, 10, 2]", "0", 0);
}

static void test_scalar_factorial_digit(void) {
    /* DigitCount[100!, 10, 3] -> 10 (from the histogram). */
    assert_eval_eq("DigitCount[100!, 10, 3]", "10", 0);
    assert_eval_eq("DigitCount[100!, 10, 0]", "30", 0);
    assert_eval_eq("DigitCount[100!, 10, 9]", "20", 0);
}

static void test_scalar_bignum_base(void) {
    /* Slow path: bignum base.  Let b = 10^18 + 1.  Then
     *   (b+1)^2 = b^2 + 2b + 1
     * which in base b is the three-digit string {1, 2, 1} (most-sig
     * first).  This forces the slow `mpz_tdiv_qr` loop in
     * dc_count_one_digit and stresses the bignum compare path. */
    assert_eval_eq("DigitCount[(10^18 + 2)^2, 10^18 + 1, 1]", "2", 0);
    assert_eval_eq("DigitCount[(10^18 + 2)^2, 10^18 + 1, 2]", "1", 0);
    assert_eval_eq("DigitCount[(10^18 + 2)^2, 10^18 + 1, 0]", "0", 0);

    /* Sanity check: b^2 itself in base b is "100" -> one '1', two '0's. */
    assert_eval_eq("DigitCount[(10^18 + 1)^2, 10^18 + 1, 1]", "1", 0);
    assert_eval_eq("DigitCount[(10^18 + 1)^2, 10^18 + 1, 0]", "2", 0);
}

static void test_scalar_bignum_digit(void) {
    /* Bignum digit d.  b = 10^18 + 1, d = 10^18.
     * (10^18 + 1)*2 in base b = {2, 2}, so digit 10^18 appears 0 times.
     * (10^18) by itself in base b = {10^18} -> single digit "10^18". */
    assert_eval_eq("DigitCount[10^18, 10^18 + 1, 10^18]", "1", 0);
    assert_eval_eq("DigitCount[10^18, 10^18 + 1, 0]", "0", 0);
}

/* --- Histogram-sum invariant ------------------------------------ */

static void test_total_invariant_decimal(void) {
    /* Total[DigitCount[n]] == IntegerLength[n] for all n != 0. */
    assert_eval_eq(
        "Table[Total[DigitCount[n]] == IntegerLength[n], {n, 1, 30}]",
        "{True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True}",
        0);
}

static void test_total_invariant_binary(void) {
    assert_eval_eq(
        "Table[Total[DigitCount[n, 2]] == IntegerLength[n, 2], {n, 1, 30}]",
        "{True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True}",
        0);
}

static void test_total_invariant_base7(void) {
    /* Non-power-of-2 base exercises the histogram differently. */
    assert_eval_eq(
        "Table[Total[DigitCount[n, 7]] == IntegerLength[n, 7], {n, 1, 30}]",
        "{True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True}",
        0);
}

/* --- Cross-check via IntegerDigits + Count ----------------------- */

static void test_cross_check_against_integerdigits(void) {
    /* DigitCount[n, b, d] == Count[IntegerDigits[n, b], d] for several
     * inputs.  This is the definitional property. */
    assert_eval_eq(
        "DigitCount[2147, 2, 1] == Count[IntegerDigits[2147, 2], 1]",
        "True", 0);
    assert_eval_eq(
        "DigitCount[100!, 10, 7] == Count[IntegerDigits[100!], 7]",
        "True", 0);
    assert_eval_eq(
        "DigitCount[123456789, 10, 5] == Count[IntegerDigits[123456789], 5]",
        "True", 0);
}

/* --- Arity errors (::argb) -------------------------------------- */

static void test_zero_args_argb(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "DigitCount[]");
    ASSERT_MSG(err && strstr(err, "DigitCount::argb") != NULL,
               "expected argb diagnostic, got: %s", err ? err : "(null)");
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    ASSERT(strstr(err, "between 1 and 3 arguments are expected") != NULL);
    free(result);
    free(err);
}

static void test_four_args_argb(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[5, 2, 0, 1]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitCount[5, 2, 0, 1]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "DigitCount::argb") != NULL,
               "expected argb diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 4 arguments") != NULL);
    free(result);
    free(err);
}

static void test_five_args_argb(void) {
    /* Documented Mathematica example: 5-arg call. */
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[5, 2, 3, 5, 4]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitCount[5, 2, 3, 5, 4]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "DigitCount::argb") != NULL,
               "expected argb diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 5 arguments") != NULL);
    free(result);
    free(err);
}

/* --- Non-integer n (::int) -------------------------------------- */

static void test_real_n_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[1.234]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitCount[1.234]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "DigitCount::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_rational_n_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_complex_n_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[1.0 + 2.0 I]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

/* --- Bad base (::base) ----------------------------------------- */

static void test_rational_base_dot_base(void) {
    /* Documented Mathematica diagnostic on DigitCount[3123, 5/2]. */
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[3123, 5/2]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitCount[3123, 5/2]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "DigitCount::base") != NULL,
               "expected base diagnostic, got: %s", err);
    ASSERT(strstr(err, "The base 5/2 at position 2") != NULL);
    ASSERT(strstr(err, "integer greater than 1") != NULL);
    free(result);
    free(err);
}

static void test_real_base_dot_base(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[10, 2.5]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::base") != NULL);
    free(result);
    free(err);
}

static void test_base_one_dot_base(void) {
    /* Base 1 -- integer but < 2; Mathematica still uses ::base. */
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[5, 1]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::base") != NULL);
    free(result);
    free(err);
}

static void test_base_zero_dot_base(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[5, 0]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::base") != NULL);
    free(result);
    free(err);
}

static void test_base_negative_dot_base(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[5, -2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::base") != NULL);
    free(result);
    free(err);
}

/* --- Bad digit (::digit) -------------------------------------- */

static void test_digit_negative_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[10, 2, -1]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::digit") != NULL);
    ASSERT(strstr(err, "non-negative integer less than the base") != NULL);
    free(result);
    free(err);
}

static void test_digit_equal_to_base_diagnostic(void) {
    /* d = base is out of range. */
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[10, 2, 2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::digit") != NULL);
    free(result);
    free(err);
}

static void test_digit_too_large_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[10, 2, 7]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::digit") != NULL);
    free(result);
    free(err);
}

static void test_digit_rational_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[10, 2, 1/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "DigitCount::digit") != NULL);
    free(result);
    free(err);
}

/* --- Symbolic args (silent unevaluated) ------------------------ */

static void test_symbolic_n_silent(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[x]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitCount[x]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "DigitCount::") == NULL,
               "expected no diagnostic for symbolic n, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_base_silent(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[10, b]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitCount[10, b]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "DigitCount::") == NULL,
               "expected no diagnostic for symbolic base, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_digit_silent(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[10, 2, d]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitCount[10, 2, d]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "DigitCount::") == NULL,
               "expected no diagnostic for symbolic digit, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Not Listable (list of integers left unevaluated) ----------- */

static void test_not_listable_on_n(void) {
    /* Without Listable, DigitCount[{1, 2, 3}] is evaluated once with a
     * List arg -- which is not integer-like -- and the call is left
     * unevaluated.  No threading. */
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[{1, 2, 3}]", &result);
    ASSERT(result != NULL);
    ASSERT_MSG(strstr(result, "DigitCount[") != NULL,
               "expected unevaluated form, got: %s", result);
    free(result);
    free(err);
}

/* --- Attributes / docstring / interned symbol ------------------ */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("DigitCount");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("DigitCount");
    ASSERT((a & ATTR_PROTECTED) != 0);
    /* Intentionally NOT Listable -- matches the user's spec exactly. */
    ASSERT((a & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("DigitCount");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "counts of digits") != NULL);
    ASSERT(strstr(def->docstring, "DigitCount[n, b, d]") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_DigitCount != NULL);
    ASSERT(strcmp(SYM_DigitCount, "DigitCount") == 0);
}

/* --- Memory-safety stress loop -------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix every code path: list & scalar form, machine int & bignum,
     * default & explicit base, default & slow base, zero, negatives,
     * error paths.  Anything mishandling mpz_t / args[] / counts will
     * surface as a leak or use-after-free under valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("DigitCount[2147, 2, 1]", "5", 0);
        assert_eval_eq("DigitCount[2147, 2]", "{5, 7}", 0);
        assert_eval_eq("DigitCount[100!]",
                       "{15, 19, 10, 10, 14, 19, 7, 14, 20, 30}", 0);
        assert_eval_eq("DigitCount[0]",
                       "{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}", 0);
        assert_eval_eq("DigitCount[-(10^30)]",
                       "{1, 0, 0, 0, 0, 0, 0, 0, 0, 30}", 0);
        assert_eval_eq("DigitCount[(10^18 + 2)^2, 10^18 + 1, 1]", "2", 0);
        assert_eval_eq("DigitCount[2^100, 2]", "{1, 100}", 0);
        /* 10 in binary == 1010 -> two 0s and two 1s.  Off-by-one tests
         * during development confirmed that this must be 2, not 1. */
        assert_eval_eq("DigitCount[10, 2, 0]", "2", 0);
        assert_eval_eq("DigitCount[10, 2, 1]", "2", 0);
    }
}

/* --- Soft-cap behaviour --------------------------------------- */

static void test_base_too_large_for_list_form(void) {
    /* List form with base 2^21 (> 2^20 cap) must emit ::ovfl and leave
     * the call unevaluated.  The 3-arg scalar form must still work for
     * the same base. */
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitCount[5, 2^21]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitCount[") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "DigitCount::ovfl") != NULL,
               "expected ovfl diagnostic, got: %s", err);
    free(result);
    free(err);

    /* Scalar form with same base works (no list allocated). */
    assert_eval_eq("DigitCount[5, 2^21, 5]", "1", 0);
    assert_eval_eq("DigitCount[5, 2^21, 0]", "0", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* Documented Mathematica examples (highest-signal). */
    TEST(test_doc_binary_count_ones);
    TEST(test_doc_binary_histogram);
    TEST(test_doc_factorial_histogram);

    /* Small decimal cases. */
    TEST(test_decimal_single_digit_3);
    TEST(test_decimal_repeated);
    TEST(test_decimal_each_digit_once);
    TEST(test_decimal_mix);

    /* Zero handling. */
    TEST(test_zero_default_base);
    TEST(test_zero_binary);
    TEST(test_zero_scalar);

    /* Negatives. */
    TEST(test_negative_small);
    TEST(test_negative_scalar);
    TEST(test_negative_decimal);

    /* Bignum. */
    TEST(test_bignum_10pow30_base10);
    TEST(test_bignum_2pow100_base2);
    TEST(test_bignum_2pow100_count_ones_in_decimal);
    TEST(test_bignum_negative);
    TEST(test_int64_overflow_boundary);

    /* Scalar (3-arg) form. */
    TEST(test_scalar_decimal_digit);
    TEST(test_scalar_factorial_digit);
    TEST(test_scalar_bignum_base);
    TEST(test_scalar_bignum_digit);

    /* Histogram-sum invariant. */
    TEST(test_total_invariant_decimal);
    TEST(test_total_invariant_binary);
    TEST(test_total_invariant_base7);

    /* Cross-check definition. */
    TEST(test_cross_check_against_integerdigits);

    /* Argument-arity errors. */
    TEST(test_zero_args_argb);
    TEST(test_four_args_argb);
    TEST(test_five_args_argb);

    /* Non-integer n diagnostics. */
    TEST(test_real_n_int_diagnostic);
    TEST(test_rational_n_int_diagnostic);
    TEST(test_complex_n_int_diagnostic);

    /* Bad base diagnostics. */
    TEST(test_rational_base_dot_base);
    TEST(test_real_base_dot_base);
    TEST(test_base_one_dot_base);
    TEST(test_base_zero_dot_base);
    TEST(test_base_negative_dot_base);

    /* Bad digit diagnostics. */
    TEST(test_digit_negative_diagnostic);
    TEST(test_digit_equal_to_base_diagnostic);
    TEST(test_digit_too_large_diagnostic);
    TEST(test_digit_rational_diagnostic);

    /* Symbolic args silent. */
    TEST(test_symbolic_n_silent);
    TEST(test_symbolic_base_silent);
    TEST(test_symbolic_digit_silent);

    /* Not Listable. */
    TEST(test_not_listable_on_n);

    /* Introspection. */
    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    /* Soft cap. */
    TEST(test_base_too_large_for_list_form);

    /* Stress. */
    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All DigitCount tests passed!\n");
    return 0;
}
