/* Unit tests for IntegerLength.
 *
 *   IntegerLength[n]            -- number of decimal digits in |n|.
 *   IntegerLength[n, b]         -- number of base-b digits in |n|, b >= 2.
 *
 * Coverage:
 *   - Decimal-digit count: small, single-digit, 10, 999, multi-digit.
 *   - Binary, base-8, hex (base 16), base-256 (mpz_sizeinbase path) and
 *     base-100 (non-power-of-2 verification path).
 *   - Zero special-case (IntegerLength[0] == 0 in any base).
 *   - Negative integers (sign discarded).
 *   - Bignum inputs: 10^30, 2^100, 100! including the Mathematica
 *     documented example IntegerLength[100!, 2] == 525, and the
 *     table {IntegerLength[100!, n], n=2..20}.
 *   - Boundary at int64 overflow (2^63 -> 19, 2^63 - 1 -> 19).
 *   - Arbitrary-precision base (base = 10^18 + 1, larger than 62, so
 *     forces the slow division path).
 *   - Listable threading on n, on b, and on both.
 *   - Error paths: 0 args, 5 args (::argt), real / rational / complex
 *     first arg (::int), real / negative / zero / one base (::int /
 *     ::ibase), symbolic args (silent unevaluated).
 *   - Attribute, docstring, and interned-symbol introspection.
 *   - Floor[Log[b, n]] + 1 equivalence as a property test.
 *   - Repeated-evaluation stress loop to catch double-frees under
 *     valgrind.
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
 * fixed temp file path; safe because tests run serially. */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_int_length_stderr.log";
    fflush(stderr);
    if (!freopen(path, "w+", stderr)) {
        if (out_result_str) *out_result_str = NULL;
        return NULL;
    }

    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    if (out_result_str) *out_result_str = expr_to_string(e);
    expr_free(p);
    expr_free(e);

    fflush(stderr);
    freopen("/dev/tty", "w", stderr);

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

/* --- Decimal --------------------------------------------------------- */

static void test_decimal_basic(void) {
    assert_eval_eq("IntegerLength[123456789]", "9", 0);
}

static void test_decimal_single_digit(void) {
    assert_eval_eq("IntegerLength[1]", "1", 0);
    assert_eval_eq("IntegerLength[9]", "1", 0);
}

static void test_decimal_ten(void) {
    assert_eval_eq("IntegerLength[10]", "2", 0);
}

static void test_decimal_nine_nines(void) {
    assert_eval_eq("IntegerLength[999]", "3", 0);
}

static void test_decimal_powers_of_ten(void) {
    /* 10^k has k+1 digits in base 10. */
    assert_eval_eq("IntegerLength[100]", "3", 0);
    assert_eval_eq("IntegerLength[1000]", "4", 0);
    assert_eval_eq("IntegerLength[100000]", "6", 0);
}

static void test_decimal_one_less_than_power(void) {
    /* 10^k - 1 has k digits.  Stresses the off-by-one verification
     * inside intlen_count_digits. */
    assert_eval_eq("IntegerLength[99]", "2", 0);
    assert_eval_eq("IntegerLength[9999]", "4", 0);
    assert_eval_eq("IntegerLength[99999999]", "8", 0);
}

/* --- Other bases ----------------------------------------------------- */

static void test_base2_small(void) {
    /* 7 in binary = 111 -> 3 digits. */
    assert_eval_eq("IntegerLength[7, 2]", "3", 0);
}

static void test_base2_at_boundary(void) {
    /* 2^k requires k+1 bits; 2^k - 1 requires k bits. */
    assert_eval_eq("IntegerLength[8, 2]", "4", 0);
    assert_eval_eq("IntegerLength[7, 2]", "3", 0);
    assert_eval_eq("IntegerLength[1024, 2]", "11", 0);
    assert_eval_eq("IntegerLength[1023, 2]", "10", 0);
}

static void test_base8(void) {
    /* 58127 in octal = 161417, 6 digits. */
    assert_eval_eq("IntegerLength[58127, 8]", "6", 0);
}

static void test_base16(void) {
    /* 58127 in hex = E30F, 4 digits. */
    assert_eval_eq("IntegerLength[58127, 16]", "4", 0);
}

static void test_base16_one_less(void) {
    /* 15 = F (one hex digit); 16 = 10 (two hex digits). */
    assert_eval_eq("IntegerLength[15, 16]", "1", 0);
    assert_eval_eq("IntegerLength[16, 16]", "2", 0);
}

static void test_base256(void) {
    /* 58127 in base 256: 0xE30F -> {227, 15}, 2 digits. */
    assert_eval_eq("IntegerLength[58127, 256]", "2", 0);
}

static void test_base_nonpower_of_two(void) {
    /* base = 100 is not a power of 2 -- exercises the verification
     * branch (mpz_sizeinbase could be 1 too big for non-power-of-2
     * bases <= 62, but here base > 62 so it falls into the slow path). */
    assert_eval_eq("IntegerLength[100, 100]", "2", 0);
    assert_eval_eq("IntegerLength[9999, 100]", "2", 0);
    assert_eval_eq("IntegerLength[10000, 100]", "3", 0);
}

static void test_base_max_supported_fast_path(void) {
    /* base = 62 is the upper limit of mpz_sizeinbase's "fast" range.
     * Also a non-power-of-2 -> verification triggers. */
    assert_eval_eq("IntegerLength[61, 62]", "1", 0);
    assert_eval_eq("IntegerLength[62, 62]", "2", 0);
    assert_eval_eq("IntegerLength[62*62 - 1, 62]", "2", 0);
    assert_eval_eq("IntegerLength[62*62, 62]", "3", 0);
}

static void test_base_three_verification(void) {
    /* base 3 hammers the non-power-of-2 verification path with a small
     * base.  3^k requires k+1 ternary digits; 3^k - 1 requires k. */
    assert_eval_eq("IntegerLength[3, 3]", "2", 0);
    assert_eval_eq("IntegerLength[2, 3]", "1", 0);
    assert_eval_eq("IntegerLength[9, 3]", "3", 0);
    assert_eval_eq("IntegerLength[8, 3]", "2", 0);
    assert_eval_eq("IntegerLength[27, 3]", "4", 0);
    assert_eval_eq("IntegerLength[26, 3]", "3", 0);
    assert_eval_eq("IntegerLength[243, 3]", "6", 0);
    assert_eval_eq("IntegerLength[242, 3]", "5", 0);
}

/* --- Zero ----------------------------------------------------------- */

static void test_zero_decimal(void) {
    /* Mathematica: IntegerLength[0] == 0. */
    assert_eval_eq("IntegerLength[0]", "0", 0);
}

static void test_zero_any_base(void) {
    assert_eval_eq("IntegerLength[0, 2]", "0", 0);
    assert_eval_eq("IntegerLength[0, 16]", "0", 0);
    assert_eval_eq("IntegerLength[0, 1000]", "0", 0);
}

/* --- Sign discarded ------------------------------------------------- */

static void test_negative(void) {
    assert_eval_eq("IntegerLength[-123456789]", "9", 0);
    assert_eval_eq("IntegerLength[-1]", "1", 0);
    assert_eval_eq("IntegerLength[-7, 2]", "3", 0);
}

/* --- Bignum -------------------------------------------------------- */

static void test_bignum_10pow30(void) {
    /* 10^30 has 31 decimal digits. */
    assert_eval_eq("IntegerLength[10^30]", "31", 0);
}

static void test_bignum_10pow30_minus_1(void) {
    /* 10^30 - 1 has 30 decimal digits. */
    assert_eval_eq("IntegerLength[10^30 - 1]", "30", 0);
}

static void test_bignum_2pow100_base2(void) {
    /* 2^100 in base 2 is 1 followed by 100 zeros -> 101 digits. */
    assert_eval_eq("IntegerLength[2^100, 2]", "101", 0);
}

static void test_bignum_2pow100_base10(void) {
    /* 2^100 = 1267650600228229401496703205376, which has 31 digits. */
    assert_eval_eq("IntegerLength[2^100]", "31", 0);
}

static void test_bignum_100_factorial_base10(void) {
    /* 100! has 158 decimal digits. */
    assert_eval_eq("IntegerLength[100!]", "158", 0);
}

static void test_bignum_100_factorial_base2(void) {
    /* Documented Mathematica example: IntegerLength[100!, 2] -> 525. */
    assert_eval_eq("IntegerLength[100!, 2]", "525", 0);
}

static void test_bignum_100_factorial_table(void) {
    /* Documented Mathematica example:
     *   Table[IntegerLength[100!, n], {n, 2, 20}]
     *   -> {525, 332, 263, 227, 204, 187, 175, 166, 158, 152,
     *       147, 142, 138, 135, 132, 129, 126, 124, 122} */
    assert_eval_eq(
        "Table[IntegerLength[100!, n], {n, 2, 20}]",
        "{525, 332, 263, 227, 204, 187, 175, 166, 158, 152,"
        " 147, 142, 138, 135, 132, 129, 126, 124, 122}",
        0);
}

static void test_bignum_negative(void) {
    /* Sign discarded on a bignum input. */
    assert_eval_eq("IntegerLength[-(10^30)]", "31", 0);
}

static void test_int64_overflow_boundary(void) {
    /* 2^63 is the first value past INT64_MAX -> stored as bigint.
     * Its decimal expansion is 9223372036854775808 -> 19 digits.
     * 2^63 - 1 = 9223372036854775807 is also 19 digits, but is a
     * machine int (INT64_MAX). */
    assert_eval_eq("IntegerLength[2^63]", "19", 0);
    assert_eval_eq("IntegerLength[2^63 - 1]", "19", 0);
    assert_eval_eq("IntegerLength[2^63, 2]", "64", 0);
    assert_eval_eq("IntegerLength[2^63 - 1, 2]", "63", 0);
}

static void test_bignum_base(void) {
    /* Slow-path: base > 62 forces the repeated-division loop.
     * b = 10^18 + 1; b^2 - 1 < b^2, b^2 has 3 digits in base b. */
    assert_eval_eq("IntegerLength[(10^18 + 1)^2, 10^18 + 1]", "3", 0);
    assert_eval_eq("IntegerLength[(10^18 + 1)^2 - 1, 10^18 + 1]", "2", 0);
    assert_eval_eq("IntegerLength[10^18 + 1, 10^18 + 1]", "2", 0);
    assert_eval_eq("IntegerLength[10^18, 10^18 + 1]", "1", 0);
}

/* --- Listable threading ------------------------------------------- */

static void test_listable_on_n(void) {
    assert_eval_eq("IntegerLength[{1, 10, 100, 1000}]", "{1, 2, 3, 4}", 0);
}

static void test_listable_on_base(void) {
    /* IntegerLength[8, {2, 3, 4, 8, 16}]: 8 has 4 binary digits (1000),
     * 3 ternary digits (22; 8 = 2*3 + 2 = 22_3 -> 2 digits actually).
     *   8 base 2: 1000      -> 4
     *   8 base 3: 22        -> 2
     *   8 base 4: 20        -> 2
     *   8 base 8: 10        -> 2
     *   8 base 16: 8        -> 1
     */
    assert_eval_eq("IntegerLength[8, {2, 3, 4, 8, 16}]",
                   "{4, 2, 2, 2, 1}", 0);
}

static void test_listable_on_both(void) {
    /* {2, 3, 4} paired element-wise with {2, 3, 4}:
     *   IntegerLength[2, 2] = 2
     *   IntegerLength[3, 3] = 2
     *   IntegerLength[4, 4] = 2
     */
    assert_eval_eq("IntegerLength[{2, 3, 4}, {2, 3, 4}]",
                   "{2, 2, 2}", 0);
}

static void test_listable_nested(void) {
    /* Threading reaches into nested lists. */
    assert_eval_eq("IntegerLength[{{1, 10}, {100, 1000}}]",
                   "{{1, 2}, {3, 4}}", 0);
}

/* --- Floor[Log[b, n]] + 1 property test -------------------------- */

static void test_floor_log_identity(void) {
    /* IntegerLength[n, b] == Floor[Log[b, n]] + 1 for n > 0, b > 1. */
    assert_eval_eq(
        "Table[IntegerLength[n, 2] == Floor[Log[2, n]] + 1, {n, 1, 30}]",
        "{True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True}",
        0);
    assert_eval_eq(
        "Table[IntegerLength[n, 10] == Floor[Log[10, n]] + 1, {n, 1, 30}]",
        "{True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True}",
        0);
}

/* --- Error paths --------------------------------------------------- */

static void test_zero_args_argt(void) {
    /* Documented Mathematica diagnostic:
     *   IntegerLength::argt: IntegerLength called with 0 arguments;
     *   1 or 2 arguments are expected. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerLength[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerLength[]");
    ASSERT_MSG(err && strstr(err, "IntegerLength::argt") != NULL,
               "expected argt diagnostic, got: %s", err ? err : "(null)");
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    ASSERT(strstr(err, "1 or 2 arguments are expected") != NULL);
    free(result);
    free(err);
}

static void test_three_args_argt(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerLength[1, 2, 3]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerLength") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerLength::argt") != NULL,
               "expected argt diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 3 arguments") != NULL);
    free(result);
    free(err);
}

static void test_five_args_argt(void) {
    /* Documented Mathematica example: 5-arg call. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerLength[1, 2, 3, 4, 5]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerLength[1, 2, 3, 4, 5]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerLength::argt") != NULL,
               "expected argt diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 5 arguments") != NULL);
    free(result);
    free(err);
}

static void test_real_n_int_diagnostic(void) {
    /* Documented Mathematica diagnostic:
     *   IntegerLength::int: Integer expected at position 1 in
     *   IntegerLength[1.1234]. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerLength[1.1234]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerLength[1.1234]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerLength::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    ASSERT(strstr(err, "IntegerLength[1.1234]") != NULL);
    free(result);
    free(err);
}

static void test_complex_n_int_diagnostic(void) {
    /* Documented Mathematica example: IntegerLength[1.1234 - 9 I]. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerLength[1.1234 - 9 I]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerLength::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_rational_n_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerLength[3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerLength::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_real_base_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerLength[10, 2.5]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerLength::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 2") != NULL);
    free(result);
    free(err);
}

static void test_bad_base_one(void) {
    /* Base 1 prints IntegerLength::ibase and leaves the call unevaluated. */
    const char* in = "IntegerLength[5, 1]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerLength") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_bad_base_zero(void) {
    const char* in = "IntegerLength[5, 0]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerLength") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_bad_base_negative(void) {
    const char* in = "IntegerLength[5, -2]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerLength") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_symbolic_n_left_unevaluated(void) {
    /* Symbolic first arg: silent (no diagnostic), call retained. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerLength[x]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerLength[x]") != NULL);
    /* No diagnostic printed for symbolic input. */
    ASSERT_MSG(err == NULL || strstr(err, "IntegerLength::") == NULL,
               "expected no diagnostic for symbolic input, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_base_left_unevaluated(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerLength[10, b]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerLength[10, b]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "IntegerLength::") == NULL,
               "expected no diagnostic for symbolic base, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Attribute / docstring / interned-symbol introspection -------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("IntegerLength");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("IntegerLength");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("IntegerLength");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "decimal digits") != NULL);
    ASSERT(strstr(def->docstring, "base b") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_IntegerLength != NULL);
    ASSERT(strcmp(SYM_IntegerLength, "IntegerLength") == 0);
}

/* --- Memory-safety stress loop ------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix machine-int / bignum / listable / power-of-2 / non-power-of-2 /
     * slow-base / error paths so any mishandled mpz_t or args[] slot
     * surfaces under valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("IntegerLength[123456789]", "9", 0);
        assert_eval_eq("IntegerLength[58127, 16]", "4", 0);
        assert_eval_eq("IntegerLength[0]", "0", 0);
        assert_eval_eq("IntegerLength[-7, 2]", "3", 0);
        assert_eval_eq("IntegerLength[10^30]", "31", 0);
        assert_eval_eq("IntegerLength[2^63]", "19", 0);
        assert_eval_eq("IntegerLength[{1, 10, 100, 1000}]", "{1, 2, 3, 4}", 0);
        assert_eval_eq("IntegerLength[100, 100]", "2", 0);
        assert_eval_eq("IntegerLength[(10^18 + 1)^2, 10^18 + 1]", "3", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_decimal_basic);
    TEST(test_decimal_single_digit);
    TEST(test_decimal_ten);
    TEST(test_decimal_nine_nines);
    TEST(test_decimal_powers_of_ten);
    TEST(test_decimal_one_less_than_power);

    TEST(test_base2_small);
    TEST(test_base2_at_boundary);
    TEST(test_base8);
    TEST(test_base16);
    TEST(test_base16_one_less);
    TEST(test_base256);
    TEST(test_base_nonpower_of_two);
    TEST(test_base_max_supported_fast_path);
    TEST(test_base_three_verification);

    TEST(test_zero_decimal);
    TEST(test_zero_any_base);

    TEST(test_negative);

    TEST(test_bignum_10pow30);
    TEST(test_bignum_10pow30_minus_1);
    TEST(test_bignum_2pow100_base2);
    TEST(test_bignum_2pow100_base10);
    TEST(test_bignum_100_factorial_base10);
    TEST(test_bignum_100_factorial_base2);
    TEST(test_bignum_100_factorial_table);
    TEST(test_bignum_negative);
    TEST(test_int64_overflow_boundary);
    TEST(test_bignum_base);

    TEST(test_listable_on_n);
    TEST(test_listable_on_base);
    TEST(test_listable_on_both);
    TEST(test_listable_nested);

    TEST(test_floor_log_identity);

    TEST(test_zero_args_argt);
    TEST(test_three_args_argt);
    TEST(test_five_args_argt);
    TEST(test_real_n_int_diagnostic);
    TEST(test_complex_n_int_diagnostic);
    TEST(test_rational_n_int_diagnostic);
    TEST(test_real_base_int_diagnostic);
    TEST(test_bad_base_one);
    TEST(test_bad_base_zero);
    TEST(test_bad_base_negative);
    TEST(test_symbolic_n_left_unevaluated);
    TEST(test_symbolic_base_left_unevaluated);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All IntegerLength tests passed!\n");
    return 0;
}
