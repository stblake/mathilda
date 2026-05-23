/* Unit tests for IntegerString.
 *
 *   IntegerString[n]            -- decimal-digit string of |n|.
 *   IntegerString[n, b]         -- base-b digit string (2 <= b <= 36).
 *                                  Digit values 10..35 use 'a'..'z'.
 *   IntegerString[n, b, len]    -- left-pad with '0' to length len; if
 *                                  |n| has more than len base-b digits,
 *                                  the len least-significant digits are
 *                                  returned.
 *
 * Coverage:
 *   - Documented Mathematica examples verbatim
 *     (17651 in base 2; 50! in base 16 and base 36).
 *   - Default decimal form, single-digit, multi-digit, multi-line.
 *   - Base 2, 8, 10, 16, 36 (boundary), and mid-range bases.
 *   - Padding with `len > digits`, exact match `len == digits`,
 *     truncation `len < digits`, and `len == 0` -> "".
 *   - Negative inputs (sign discarded) in the 1/2/3-arg forms.
 *   - Zero special cases under every form and base.
 *   - Bignum inputs: 100!, 2^100, -10^30, and IntegerString[10^k, b, len]
 *     padding boundary.
 *   - Listable threading on n, on base, on both, and nested lists.
 *   - Round-trip with FromDigits / IntegerDigits.
 *   - Error paths: 0 args, 4 args (::argb); Real / Rational / Complex n
 *     (::int); Real / negative / zero / one / 37 / 50 base (::basf);
 *     negative len (::intnn); symbolic args (silent unevaluated).
 *   - Attribute, docstring, and interned-symbol introspection.
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
    const char* path = "/tmp/mathilda_int_string_stderr.log";
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

/* --- Documented Mathematica examples ------------------------------- */

static void test_documented_binary(void) {
    /* Spec: IntegerString[17651, 2] -> "100010011110011" */
    assert_eval_eq("IntegerString[17651, 2]", "\"100010011110011\"", 0);
}

static void test_documented_50_factorial_base16(void) {
    /* Spec: IntegerString[50!, 16] ->
     * "49eebc961ed279b02b1ef4f28d19a84f5973a1d2c7800000000000" */
    assert_eval_eq(
        "IntegerString[50!, 16]",
        "\"49eebc961ed279b02b1ef4f28d19a84f5973a1d2c7800000000000\"", 0);
}

static void test_documented_50_factorial_base36(void) {
    /* Spec: IntegerString[50!, 36] ->
     * "4q7eyp9zizmtqt0648txt4fm720cc1s00000000000" */
    assert_eval_eq("IntegerString[50!, 36]",
                   "\"4q7eyp9zizmtqt0648txt4fm720cc1s00000000000\"", 0);
}

/* --- Default-base (10) -------------------------------------------- */

static void test_decimal_single_digit(void) {
    assert_eval_eq("IntegerString[0]", "\"0\"", 0);
    assert_eval_eq("IntegerString[1]", "\"1\"", 0);
    assert_eval_eq("IntegerString[9]", "\"9\"", 0);
}

static void test_decimal_multi_digit(void) {
    assert_eval_eq("IntegerString[10]", "\"10\"", 0);
    assert_eval_eq("IntegerString[12345]", "\"12345\"", 0);
    assert_eval_eq("IntegerString[1000000]", "\"1000000\"", 0);
}

static void test_decimal_negative_discards_sign(void) {
    /* Spec: IntegerString[n] discards the sign of n. */
    assert_eval_eq("IntegerString[-1]", "\"1\"", 0);
    assert_eval_eq("IntegerString[-12345]", "\"12345\"", 0);
    assert_eval_eq("IntegerString[-1000000]", "\"1000000\"", 0);
}

/* --- Base 2 -------------------------------------------------------- */

static void test_binary_small(void) {
    assert_eval_eq("IntegerString[0, 2]", "\"0\"", 0);
    assert_eval_eq("IntegerString[1, 2]", "\"1\"", 0);
    assert_eval_eq("IntegerString[2, 2]", "\"10\"", 0);
    assert_eval_eq("IntegerString[7, 2]", "\"111\"", 0);
    assert_eval_eq("IntegerString[8, 2]", "\"1000\"", 0);
    assert_eval_eq("IntegerString[255, 2]", "\"11111111\"", 0);
    assert_eval_eq("IntegerString[256, 2]", "\"100000000\"", 0);
}

static void test_binary_negative(void) {
    assert_eval_eq("IntegerString[-255, 2]", "\"11111111\"", 0);
}

/* --- Bases 8 and 16 ----------------------------------------------- */

static void test_octal(void) {
    assert_eval_eq("IntegerString[8, 8]", "\"10\"", 0);
    assert_eval_eq("IntegerString[63, 8]", "\"77\"", 0);
    assert_eval_eq("IntegerString[64, 8]", "\"100\"", 0);
    assert_eval_eq("IntegerString[511, 8]", "\"777\"", 0);
}

static void test_hex_basic(void) {
    /* Spec: "In base 16, a through f are used as digits". */
    assert_eval_eq("IntegerString[10, 16]", "\"a\"", 0);
    assert_eval_eq("IntegerString[15, 16]", "\"f\"", 0);
    assert_eval_eq("IntegerString[16, 16]", "\"10\"", 0);
    assert_eval_eq("IntegerString[255, 16]", "\"ff\"", 0);
    assert_eval_eq("IntegerString[256, 16]", "\"100\"", 0);
    assert_eval_eq("IntegerString[4095, 16]", "\"fff\"", 0);
    assert_eval_eq("IntegerString[65535, 16]", "\"ffff\"", 0);
}

/* --- Base 36 (maximum) -------------------------------------------- */

static void test_base36_boundary(void) {
    /* Digit 35 prints as 'z' in base 36; 36 wraps to "10". */
    assert_eval_eq("IntegerString[35, 36]", "\"z\"", 0);
    assert_eval_eq("IntegerString[36, 36]", "\"10\"", 0);
    /* 36^2 - 1 = 1295 -> "zz" */
    assert_eval_eq("IntegerString[1295, 36]", "\"zz\"", 0);
    /* 36^2 = 1296 -> "100" */
    assert_eval_eq("IntegerString[1296, 36]", "\"100\"", 0);
}

/* --- Mid-range bases ---------------------------------------------- */

static void test_base3(void) {
    /* 0,1,2,10,11,12,20,21,22,100 in base 3 */
    assert_eval_eq("IntegerString[0, 3]", "\"0\"", 0);
    assert_eval_eq("IntegerString[3, 3]", "\"10\"", 0);
    assert_eval_eq("IntegerString[9, 3]", "\"100\"", 0);
    assert_eval_eq("IntegerString[26, 3]", "\"222\"", 0);
    assert_eval_eq("IntegerString[27, 3]", "\"1000\"", 0);
}

static void test_base7(void) {
    /* 7^3 = 343 -> "1000"; 7^3 + 1 -> "1001"; 7^3 - 1 -> "666" */
    assert_eval_eq("IntegerString[343, 7]", "\"1000\"", 0);
    assert_eval_eq("IntegerString[342, 7]", "\"666\"", 0);
}

static void test_base11_through_15(void) {
    /* Digits 10..14 in base 11..15: each prints as 'a'..'e'. */
    assert_eval_eq("IntegerString[10, 11]", "\"a\"", 0);
    assert_eval_eq("IntegerString[11, 12]", "\"b\"", 0);
    assert_eval_eq("IntegerString[12, 13]", "\"c\"", 0);
    assert_eval_eq("IntegerString[13, 14]", "\"d\"", 0);
    assert_eval_eq("IntegerString[14, 15]", "\"e\"", 0);
}

static void test_base25(void) {
    /* 24 in base 25 -> "o" (10:a..23:n,24:o). */
    assert_eval_eq("IntegerString[24, 25]", "\"o\"", 0);
    assert_eval_eq("IntegerString[25, 25]", "\"10\"", 0);
}

/* --- Padding (3-arg form) ----------------------------------------- */

static void test_pad_decimal(void) {
    /* len > digits: pad with leading zeros. */
    assert_eval_eq("IntegerString[1, 10, 5]", "\"00001\"", 0);
    assert_eval_eq("IntegerString[42, 10, 6]", "\"000042\"", 0);
    /* len == digits: identity (no padding, no truncation). */
    assert_eval_eq("IntegerString[12345, 10, 5]", "\"12345\"", 0);
}

static void test_pad_binary(void) {
    /* 5 = 101 in binary; pad to length 8 -> "00000101". */
    assert_eval_eq("IntegerString[5, 2, 8]", "\"00000101\"", 0);
    /* 255 = 11111111 in binary; pad to length 16 -> "0000000011111111" */
    assert_eval_eq("IntegerString[255, 2, 16]", "\"0000000011111111\"", 0);
}

static void test_pad_hex(void) {
    assert_eval_eq("IntegerString[10, 16, 4]", "\"000a\"", 0);
    assert_eval_eq("IntegerString[255, 16, 4]", "\"00ff\"", 0);
    /* Width 8 hex byte string. */
    assert_eval_eq("IntegerString[0, 16, 8]", "\"00000000\"", 0);
}

static void test_pad_zero_value(void) {
    /* IntegerString[0, b, len] is len zeros. */
    assert_eval_eq("IntegerString[0, 10, 4]", "\"0000\"", 0);
    assert_eval_eq("IntegerString[0, 2, 8]", "\"00000000\"", 0);
    assert_eval_eq("IntegerString[0, 36, 3]", "\"000\"", 0);
}

static void test_pad_negative_input(void) {
    /* Sign discarded; padding applied to absolute value. */
    assert_eval_eq("IntegerString[-42, 10, 6]", "\"000042\"", 0);
    assert_eval_eq("IntegerString[-5, 2, 8]", "\"00000101\"", 0);
}

/* --- Truncation (len < digits) ------------------------------------ */

static void test_trunc_keeps_least_significant_decimal(void) {
    /* Spec: "If len is less than the number of digits in n, then the
     * len least significant digits are returned." */
    assert_eval_eq("IntegerString[12345, 10, 3]", "\"345\"", 0);
    assert_eval_eq("IntegerString[12345, 10, 1]", "\"5\"", 0);
    assert_eval_eq("IntegerString[1000000, 10, 4]", "\"0000\"", 0);
}

static void test_trunc_keeps_least_significant_binary(void) {
    /* 13 = "1101"; keep last 2 -> "01". */
    assert_eval_eq("IntegerString[13, 2, 2]", "\"01\"", 0);
    /* 255 = "11111111"; keep last 4 -> "1111". */
    assert_eval_eq("IntegerString[255, 2, 4]", "\"1111\"", 0);
}

static void test_len_zero_returns_empty(void) {
    /* IntegerString[n, b, 0] -> "" for any n. */
    assert_eval_eq("IntegerString[12345, 10, 0]", "\"\"", 0);
    assert_eval_eq("IntegerString[0, 10, 0]", "\"\"", 0);
    assert_eval_eq("IntegerString[-7, 2, 0]", "\"\"", 0);
}

/* --- Zero special cases ------------------------------------------- */

static void test_zero_default(void) {
    assert_eval_eq("IntegerString[0]", "\"0\"", 0);
}

static void test_zero_each_base(void) {
    assert_eval_eq("IntegerString[0, 2]", "\"0\"", 0);
    assert_eval_eq("IntegerString[0, 7]", "\"0\"", 0);
    assert_eval_eq("IntegerString[0, 16]", "\"0\"", 0);
    assert_eval_eq("IntegerString[0, 36]", "\"0\"", 0);
}

/* --- Bignum inputs ------------------------------------------------- */

static void test_bignum_100_factorial(void) {
    /* 100! has 158 decimal digits and exactly 24 trailing zeros
     * (Legendre: exponent of 5 in 100! is 25 + 5 + 1 = 24, so the
     * trailing-zero count is min(97, 24) = 24). */
    assert_eval_eq(
        "IntegerString[100!]",
        "\"93326215443944152681699238856266700490715968264381621468592963"
        "8952175999932299156089414639761565182862536979208272237582511852"
        "10916864000000000000000000000000\"", 0);
}

static void test_bignum_2pow100_binary(void) {
    /* 2^100 -> "1" followed by 100 zeros in binary (101 characters).
     * Defining property: StringLength[IntegerString[2^k, 2]] == k + 1. */
    assert_eval_eq("StringLength[IntegerString[2^100, 2]]", "101", 0);
    /* And it should match the concatenation IntegerString[1, 2]
     * <> StringRepeat-like padding via 3-arg form. */
    assert_eval_eq("IntegerString[2^100, 2] == "
                   "IntegerString[2^100, 2, 101]", "True", 0);
}

static void test_bignum_2pow100_hex(void) {
    /* 2^100 = 0x1 followed by 25 hex zeros (since 100 = 4 * 25). */
    assert_eval_eq("IntegerString[2^100, 16]",
                   "\"10000000000000000000000000\"", 0);
}

static void test_bignum_negative_10pow30(void) {
    /* Sign discarded for bignums too. */
    assert_eval_eq("IntegerString[-10^30]",
                   "\"1000000000000000000000000000000\"", 0);
}

static void test_bignum_pad_boundary(void) {
    /* 10^10 has 11 decimal digits.  Pad to length 15 -> 4 leading zeros. */
    assert_eval_eq("IntegerString[10^10, 10, 15]",
                   "\"000010000000000\"", 0);
}

/* --- Listable threading -------------------------------------------- */

static void test_listable_on_n(void) {
    assert_eval_eq("IntegerString[{0, 1, 12, 255}]",
                   "{\"0\", \"1\", \"12\", \"255\"}", 0);
}

static void test_listable_on_base(void) {
    /* IntegerString[10, {2, 3, 4, 16}]:
     *   base 2: "1010", base 3: "101", base 4: "22", base 16: "a". */
    assert_eval_eq("IntegerString[10, {2, 3, 4, 16}]",
                   "{\"1010\", \"101\", \"22\", \"a\"}", 0);
}

static void test_listable_on_both(void) {
    /* Element-wise threading. */
    assert_eval_eq("IntegerString[{10, 100, 1000}, {2, 8, 16}]",
                   "{\"1010\", \"144\", \"3e8\"}", 0);
}

static void test_listable_nested(void) {
    assert_eval_eq("IntegerString[{{1, 2}, {3, 4}}, 2]",
                   "{{\"1\", \"10\"}, {\"11\", \"100\"}}", 0);
}

static void test_listable_with_padding(void) {
    /* IntegerString[Range[0, 7], 2, 3]:
     *   {0->"000", 1->"001", 2->"010", ..., 7->"111"} */
    assert_eval_eq("IntegerString[Range[0, 7], 2, 3]",
                   "{\"000\", \"001\", \"010\", \"011\", "
                   "\"100\", \"101\", \"110\", \"111\"}", 0);
}

/* --- Round-trip properties ---------------------------------------- */

static void test_roundtrip_fromdigits_decimal(void) {
    /* FromDigits[IntegerString[n]] == n for any non-negative n.  Use
     * `==` rather than printed-string compare so the bignum side (which
     * prints in full) is normalised against the algebraic form. */
    assert_eval_eq("FromDigits[IntegerString[12345]]", "12345", 0);
    assert_eval_eq("FromDigits[IntegerString[100!]] == 100!", "True", 0);
}

static void test_roundtrip_fromdigits_binary(void) {
    /* FromDigits[IntegerString[n, 2], 2] == n. */
    assert_eval_eq("FromDigits[IntegerString[12345, 2], 2]", "12345", 0);
    assert_eval_eq("FromDigits[IntegerString[2^100, 2], 2] == 2^100",
                   "True", 0);
}

static void test_roundtrip_fromdigits_base36(void) {
    assert_eval_eq("FromDigits[IntegerString[123456789, 36], 36]",
                   "123456789", 0);
}

static void test_roundtrip_with_integerdigits(void) {
    /* StringLength[IntegerString[n, b]] == IntegerLength[n, b] (n != 0). */
    assert_eval_eq("StringLength[IntegerString[12345]] == "
                   "IntegerLength[12345]", "True", 0);
    assert_eval_eq("StringLength[IntegerString[100!, 16]] == "
                   "IntegerLength[100!, 16]", "True", 0);
}

/* --- Error paths --------------------------------------------------- */

static void test_zero_args_argb(void) {
    /* Documented:
     *   IntegerString[]
     *   IntegerString::argb: IntegerString called with 0 arguments;
     *   between 1 and 3 arguments are expected.
     *   IntegerString[] */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerString[]");
    ASSERT_MSG(err && strstr(err, "IntegerString::argb") != NULL,
               "expected argb diagnostic, got: %s", err ? err : "(null)");
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    ASSERT(strstr(err, "between 1 and 3 arguments are expected") != NULL);
    free(result);
    free(err);
}

static void test_four_args_argb(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[1, 2, 3, 4]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerString[1, 2, 3, 4]");
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerString::argb") != NULL);
    ASSERT(strstr(err, "called with 4 arguments") != NULL);
    free(result);
    free(err);
}

static void test_real_n_int_diagnostic(void) {
    /* Documented:
     *   IntegerString[11.3423]
     *   IntegerString::int: Integer expected at position 1 in
     *   IntegerString[11.3423].
     *   IntegerString[11.3423] */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[11.3423]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerString[11.3423]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerString::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    ASSERT(strstr(err, "IntegerString[11.3423]") != NULL);
    free(result);
    free(err);
}

static void test_rational_n_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerString::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_complex_n_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[2 + 3 I]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerString::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_real_base_basf_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[10, 2.5]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerString::basf") != NULL,
               "expected basf diagnostic, got: %s", err);
    free(result);
    free(err);
}

static void test_bad_base_one(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[10, 1]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerString[10, 1]") != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerString::basf") != NULL);
    free(result);
    free(err);
}

static void test_bad_base_zero(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[10, 0]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerString[10, 0]") != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerString::basf") != NULL);
    free(result);
    free(err);
}

static void test_bad_base_negative(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[10, -3]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerString::basf") != NULL);
    free(result);
    free(err);
}

static void test_bad_base_37(void) {
    /* Spec: max base is 36. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[10, 37]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerString[10, 37]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerString::basf") != NULL,
               "expected basf for base 37, got: %s", err);
    free(result);
    free(err);
}

static void test_bad_base_50(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[10, 50]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerString::basf") != NULL);
    free(result);
    free(err);
}

static void test_negative_len_intnn(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[12, 10, -3]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerString[12, 10, -3]") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerString::intnn") != NULL,
               "expected intnn diagnostic, got: %s", err);
    free(result);
    free(err);
}

static void test_real_len_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[12, 10, 2.5]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerString::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 3") != NULL);
    free(result);
    free(err);
}

static void test_symbolic_n_left_unevaluated(void) {
    /* Symbolic n: silent, call retained. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[x]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerString[x]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "IntegerString::") == NULL,
               "expected no diagnostic for symbolic input, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_base_left_unevaluated(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[12, b]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerString[12, b]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "IntegerString::") == NULL,
               "expected no diagnostic for symbolic base, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_len_left_unevaluated(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerString[12, 10, k]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerString[12, 10, k]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "IntegerString::") == NULL,
               "expected no diagnostic for symbolic len, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Attribute / docstring / interned-symbol introspection -------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("IntegerString");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("IntegerString");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("IntegerString");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "decimal digits") != NULL);
    ASSERT(strstr(def->docstring, "base-b digits") != NULL);
    ASSERT(strstr(def->docstring, "maximum allowed base is 36") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_IntegerString != NULL);
    ASSERT(strcmp(SYM_IntegerString, "IntegerString") == 0);
}

/* --- Memory-safety stress loop ------------------------------------ */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix machine-int / bignum / 1/2/3-arg / pad / truncate / negative /
     * zero / listable / fast and slow base paths. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("IntegerString[17651, 2]", "\"100010011110011\"", 0);
        assert_eval_eq("IntegerString[12345]", "\"12345\"", 0);
        assert_eval_eq("IntegerString[12345, 10, 3]", "\"345\"", 0);
        assert_eval_eq("IntegerString[-12345, 10, 7]", "\"0012345\"", 0);
        assert_eval_eq("IntegerString[0, 16, 5]", "\"00000\"", 0);
        assert_eval_eq("IntegerString[2^100, 16]",
                       "\"10000000000000000000000000\"", 0);
        assert_eval_eq("IntegerString[{0, 1, 12, 255}]",
                       "{\"0\", \"1\", \"12\", \"255\"}", 0);
        assert_eval_eq(
            "IntegerString[50!, 36]",
            "\"4q7eyp9zizmtqt0648txt4fm720cc1s00000000000\"", 0);
        assert_eval_eq("IntegerString[Range[0, 3], 2, 2]",
                       "{\"00\", \"01\", \"10\", \"11\"}", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_documented_binary);
    TEST(test_documented_50_factorial_base16);
    TEST(test_documented_50_factorial_base36);

    TEST(test_decimal_single_digit);
    TEST(test_decimal_multi_digit);
    TEST(test_decimal_negative_discards_sign);

    TEST(test_binary_small);
    TEST(test_binary_negative);

    TEST(test_octal);
    TEST(test_hex_basic);

    TEST(test_base36_boundary);

    TEST(test_base3);
    TEST(test_base7);
    TEST(test_base11_through_15);
    TEST(test_base25);

    TEST(test_pad_decimal);
    TEST(test_pad_binary);
    TEST(test_pad_hex);
    TEST(test_pad_zero_value);
    TEST(test_pad_negative_input);

    TEST(test_trunc_keeps_least_significant_decimal);
    TEST(test_trunc_keeps_least_significant_binary);
    TEST(test_len_zero_returns_empty);

    TEST(test_zero_default);
    TEST(test_zero_each_base);

    TEST(test_bignum_100_factorial);
    TEST(test_bignum_2pow100_binary);
    TEST(test_bignum_2pow100_hex);
    TEST(test_bignum_negative_10pow30);
    TEST(test_bignum_pad_boundary);

    TEST(test_listable_on_n);
    TEST(test_listable_on_base);
    TEST(test_listable_on_both);
    TEST(test_listable_nested);
    TEST(test_listable_with_padding);

    TEST(test_roundtrip_fromdigits_decimal);
    TEST(test_roundtrip_fromdigits_binary);
    TEST(test_roundtrip_fromdigits_base36);
    TEST(test_roundtrip_with_integerdigits);

    TEST(test_zero_args_argb);
    TEST(test_four_args_argb);
    TEST(test_real_n_int_diagnostic);
    TEST(test_rational_n_int_diagnostic);
    TEST(test_complex_n_int_diagnostic);
    TEST(test_real_base_basf_diagnostic);
    TEST(test_bad_base_one);
    TEST(test_bad_base_zero);
    TEST(test_bad_base_negative);
    TEST(test_bad_base_37);
    TEST(test_bad_base_50);
    TEST(test_negative_len_intnn);
    TEST(test_real_len_int_diagnostic);
    TEST(test_symbolic_n_left_unevaluated);
    TEST(test_symbolic_base_left_unevaluated);
    TEST(test_symbolic_len_left_unevaluated);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All IntegerString tests passed!\n");
    return 0;
}
