/* Unit tests for FromDigits.
 *
 *   FromDigits[list]            -- decimal digits, MSD-first.
 *   FromDigits[list, b]         -- base-b digits.
 *   FromDigits["string"]        -- string of digits in base 10
 *                                  (letters a-z / A-Z denote 10-35).
 *   FromDigits["string", b]     -- string of digits in base b.
 *
 * Coverage:
 *   - Documentation examples verbatim: FromDigits[{5,1,2,8}],
 *     FromDigits[{1,0,1,1,0,1,1},2], FromDigits["1923"],
 *     FromDigits["1011011",2], FromDigits[{7,11,0,0,0,122}],
 *     FromDigits["1A3C"], FromDigits[{a,b,c,d,e},x].
 *   - Inverse-of-IntegerDigits round-trip property:
 *     FromDigits[IntegerDigits[n]] == Abs[n] for both positive and
 *     negative n, across machine and bignum sizes.
 *   - Bignum results: 2^100, 10^30 + 12345, 100!.
 *   - "Carrying" of digits >= base (FromDigits[{7,11,0,0,0,122}],
 *     FromDigits[{1,2,30}, 10] == 150).
 *   - Negative digits within the list.
 *   - Symbolic digits or base produce the Horner polynomial.
 *   - Real / Rational base produce concrete numeric results.
 *   - Hex / base-16 string with mixed-case letters.
 *   - Multi-character string with no explicit base ("A1", "1A3C",
 *     "7", "D" -> single chars 7..13).
 *   - Edge cases: empty list, empty string, single zero, zero list,
 *     leading zeros.
 *   - Error paths: 0 args, 3 args (::argb), numeric / unevaluated
 *     first arg (::nlst), base < 2 over list / string (::ibase),
 *     invalid character in string (::char), symbolic base over a
 *     string (silent unevaluated), symbolic first arg
 *     (silent unevaluated).
 *   - Attribute, docstring, and interned-symbol introspection.
 *   - Repeated-evaluation stress loop to catch double-frees / leaks
 *     under valgrind.
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
    const char* path = "/tmp/mathilda_from_digits_stderr.log";
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

/* --- Documentation examples ----------------------------------------- */

static void test_doc_decimal_list(void) {
    /* FromDigits[{5,1,2,8}] -> 5128 */
    assert_eval_eq("FromDigits[{5,1,2,8}]", "5128", 0);
}

static void test_doc_binary_list(void) {
    /* FromDigits[{1,0,1,1,0,1,1},2] -> 91 */
    assert_eval_eq("FromDigits[{1,0,1,1,0,1,1},2]", "91", 0);
}

static void test_doc_decimal_string(void) {
    /* FromDigits["1923"] -> 1923 */
    assert_eval_eq("FromDigits[\"1923\"]", "1923", 0);
}

static void test_doc_binary_string(void) {
    /* FromDigits["1011011",2] -> 91 */
    assert_eval_eq("FromDigits[\"1011011\",2]", "91", 0);
}

static void test_doc_carrying(void) {
    /* FromDigits[{7,11,0,0,0,122}] -> 810122 */
    assert_eval_eq("FromDigits[{7,11,0,0,0,122}]", "810122", 0);
}

static void test_doc_string_base10_letters(void) {
    /* FromDigits["1A3C"] -> 2042 = 1*1000 + 10*100 + 3*10 + 12 */
    assert_eval_eq("FromDigits[\"1A3C\"]", "2042", 0);
}

static void test_doc_symbolic_polynomial(void) {
    /* FromDigits[{a,b,c,d,e},x] -> e + d x + c x^2 + b x^3 + a x^4
     * (Mathematica documented example, after the implicit Plus sort.) */
    assert_eval_eq("FromDigits[{a,b,c,d,e},x]",
                   "e + d x + c x^2 + b x^3 + a x^4", 0);
}

static void test_doc_single_chars(void) {
    /* From the docs: {"7","8","9","A","B","C","D"} -> {7,8,9,10,11,12,13} */
    assert_eval_eq("FromDigits /@ {\"7\",\"8\",\"9\",\"A\",\"B\",\"C\",\"D\"}",
                   "{7, 8, 9, 10, 11, 12, 13}", 0);
}

/* --- Inverse-of-IntegerDigits round trip ---------------------------- */

static void test_inverse_positive(void) {
    /* For non-negative n, FromDigits[IntegerDigits[n]] == n exactly. */
    assert_eval_eq("FromDigits[IntegerDigits[0]]", "0", 0);
    assert_eval_eq("FromDigits[IntegerDigits[1]]", "1", 0);
    assert_eval_eq("FromDigits[IntegerDigits[58127]]", "58127", 0);
    assert_eval_eq("FromDigits[IntegerDigits[1234567890]]", "1234567890", 0);
}

static void test_inverse_negative_drops_sign(void) {
    /* IntegerDigits discards the sign, so FromDigits[IntegerDigits[-n]]
     * is Abs[-n] = n, NOT -n.  Mathematica documents this verbatim. */
    assert_eval_eq("FromDigits /@ IntegerDigits[Range[-10, 10]]",
                   "{10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,"
                   " 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}", 0);
}

static void test_inverse_arbitrary_base(void) {
    /* Round-trip in several non-default bases. */
    assert_eval_eq("FromDigits[IntegerDigits[58127, 16], 16]", "58127", 0);
    assert_eval_eq("FromDigits[IntegerDigits[58127, 2], 2]", "58127", 0);
    assert_eval_eq("FromDigits[IntegerDigits[1000000, 256], 256]", "1000000", 0);
}

/* --- Bignum results ------------------------------------------------- */

static void test_bignum_2_to_100(void) {
    /* IntegerDigits[2^100] -> 31 decimal digits.  Round-trip via
     * FromDigits must produce the exact bignum value back. */
    assert_eval_eq("FromDigits[IntegerDigits[2^100]]",
                   "1267650600228229401496703205376", 0);
}

static void test_bignum_2_to_100_base2(void) {
    /* In base 2 the digit list has 101 entries, all-1 = 0 except the
     * leading 1; round-trip still exact. */
    assert_eval_eq("FromDigits[IntegerDigits[2^100, 2], 2]",
                   "1267650600228229401496703205376", 0);
}

static void test_bignum_10_to_30_plus(void) {
    /* Inputs that overflow int64 but only mildly. */
    assert_eval_eq("FromDigits[IntegerDigits[10^30 + 12345]]",
                   "1000000000000000000000000012345", 0);
}

static void test_bignum_factorial(void) {
    /* 100! has 158 decimal digits.  Round-trip via FromDigits should
     * land on the exact same bignum. */
    assert_eval_eq("FromDigits[IntegerDigits[100!]] == 100!", "True", 0);
}

static void test_bignum_negative_digit_in_list(void) {
    /* {1, 0, -1, 5} in base 10 -> 1000 + 0 - 10 + 5 == 995.
     * Validates that the integer fast path handles signed digits. */
    assert_eval_eq("FromDigits[{1, 0, -1, 5}]", "995", 0);
}

/* --- Carrying behaviour -------------------------------------------- */

static void test_carry_two_digit(void) {
    /* {1, 30} in base 10 -> 10 + 30 == 40. */
    assert_eval_eq("FromDigits[{1, 30}]", "40", 0);
}

static void test_carry_three_digit(void) {
    /* {1, 2, 30} -> 100 + 20 + 30 == 150. */
    assert_eval_eq("FromDigits[{1, 2, 30}]", "150", 0);
}

static void test_carry_large(void) {
    /* {10, 10, 10} in base 2 -> 10*4 + 10*2 + 10 == 70.
     * Stress the carry when *every* digit exceeds the base. */
    assert_eval_eq("FromDigits[{10, 10, 10}, 2]", "70", 0);
}

/* --- Symbolic / non-integer base ----------------------------------- */

static void test_symbolic_base(void) {
    /* Integer digits but symbolic base produce a polynomial in b. */
    assert_eval_eq("FromDigits[{1, 2, 3}, b]", "3 + 2 b + b^2", 0);
}

static void test_symbolic_digit_in_list(void) {
    /* Mixed integer / symbolic digit list with concrete base. */
    assert_eval_eq("FromDigits[{1, a, 3}, 10]", "103 + 10 a", 0);
}

static void test_symbolic_both(void) {
    /* All-symbolic list + symbolic base. */
    assert_eval_eq("FromDigits[{a, b}, x]", "b + a x", 0);
}

static void test_real_base(void) {
    /* Real base -> Real result via the polynomial path. */
    assert_eval_eq("FromDigits[{2, 3, 4}, 2.5]", "24.0", 0);
}

static void test_rational_base(void) {
    /* Rational base also lands on the polynomial path. */
    /* {1,2,3} in base 1/2 -> 1*(1/4) + 2*(1/2) + 3 == 1/4 + 1 + 3 == 17/4. */
    assert_eval_eq("FromDigits[{1, 2, 3}, 1/2]", "17/4", 0);
}

/* --- String forms --------------------------------------------------- */

static void test_string_hex(void) {
    /* "A1" in base 16 -> 10*16 + 1 == 161. */
    assert_eval_eq("FromDigits[\"A1\", 16]", "161", 0);
}

static void test_string_hex_mixed_case(void) {
    /* "abc" and "ABC" must agree (case-insensitive). */
    assert_eval_eq("FromDigits[\"ABC\", 16]", "2748", 0);
    assert_eval_eq("FromDigits[\"abc\", 16]", "2748", 0);
}

static void test_string_base_36(void) {
    /* "Z" in base 36 -> 35. */
    assert_eval_eq("FromDigits[\"Z\", 36]", "35", 0);
    assert_eval_eq("FromDigits[\"10\", 36]", "36", 0);
}

static void test_string_bignum(void) {
    /* A 50-digit string -> bignum result. */
    assert_eval_eq("FromDigits[\"12345678901234567890123456789012345678901234567890\"]",
                   "12345678901234567890123456789012345678901234567890", 0);
}

/* --- Edge cases ----------------------------------------------------- */

static void test_empty_list(void) {
    /* Empty digit list -> 0. */
    assert_eval_eq("FromDigits[{}]", "0", 0);
}

static void test_empty_list_with_base(void) {
    assert_eval_eq("FromDigits[{}, 2]", "0", 0);
}

static void test_empty_string(void) {
    /* Empty string of digits -> 0. */
    assert_eval_eq("FromDigits[\"\"]", "0", 0);
}

static void test_zero_list(void) {
    assert_eval_eq("FromDigits[{0}]", "0", 0);
    assert_eval_eq("FromDigits[{0, 0, 0}]", "0", 0);
}

static void test_leading_zeros(void) {
    /* Leading zeros are part of the polynomial but contribute nothing. */
    assert_eval_eq("FromDigits[{0, 0, 1, 2, 3}]", "123", 0);
}

static void test_single_digit_list(void) {
    /* Single-digit list collapses to that digit. */
    assert_eval_eq("FromDigits[{42}]", "42", 0);
    assert_eval_eq("FromDigits[{7}, 2]", "7", 0);
}

/* --- Error paths ---------------------------------------------------- */

static void test_zero_args_argb(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("FromDigits[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "FromDigits[]");
    ASSERT_MSG(err && strstr(err, "FromDigits::argb") != NULL,
               "expected argb diagnostic, got: %s", err ? err : "(null)");
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    ASSERT(strstr(err, "1 or 2 arguments are expected") != NULL);
    free(result);
    free(err);
}

static void test_three_args_argb(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("FromDigits[{1,2}, 10, 5]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "FromDigits") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "FromDigits::argb") != NULL,
               "expected argb diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 3 arguments") != NULL);
    free(result);
    free(err);
}

static void test_numeric_first_arg_nlst(void) {
    /* A plain integer is neither a list nor a string. */
    char* result = NULL;
    char* err = eval_capturing_stderr("FromDigits[5]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "FromDigits") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "FromDigits::nlst") != NULL,
               "expected nlst diagnostic, got: %s", err);
    ASSERT(strstr(err, "position 1") != NULL);
    free(result);
    free(err);
}

static void test_real_first_arg_nlst(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("FromDigits[1.5]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "FromDigits::nlst") != NULL);
    free(result);
    free(err);
}

static void test_bad_base_zero_list(void) {
    /* Base 0 over a list -> ::ibase, call left unevaluated. */
    char* result = NULL;
    char* err = eval_capturing_stderr("FromDigits[{1,2,3}, 0]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "FromDigits") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "FromDigits::ibase") != NULL,
               "expected ibase diagnostic, got: %s", err);
    free(result);
    free(err);
}

static void test_bad_base_negative_string(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("FromDigits[\"123\", -2]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "FromDigits") != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "FromDigits::ibase") != NULL);
    free(result);
    free(err);
}

static void test_invalid_string_char(void) {
    /* '$' is outside [0-9a-zA-Z] -> ::char diagnostic, unevaluated. */
    char* result = NULL;
    char* err = eval_capturing_stderr("FromDigits[\"1$3\"]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "FromDigits") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "FromDigits::char") != NULL,
               "expected char diagnostic, got: %s", err);
    ASSERT(strstr(err, "position 2") != NULL);
    free(result);
    free(err);
}

static void test_symbolic_base_over_string_silent(void) {
    /* Symbolic base over a string is left unevaluated WITHOUT any
     * diagnostic -- the user may be deliberately holding the base for
     * later substitution. */
    char* result = NULL;
    char* err = eval_capturing_stderr("FromDigits[\"123\", b]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "FromDigits") != NULL);
    ASSERT(strstr(result, "\"123\"") != NULL);
    ASSERT_MSG(err && err[0] == '\0',
               "expected silent unevaluated, got stderr: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_first_arg_silent(void) {
    /* Pure symbolic first argument flows through silently. */
    char* result = NULL;
    char* err = eval_capturing_stderr("FromDigits[x]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "FromDigits[x]");
    ASSERT_MSG(err && err[0] == '\0',
               "expected silent unevaluated, got stderr: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Attribute / docstring / interned-symbol introspection ---------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("FromDigits");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("FromDigits");
    ASSERT((a & ATTR_PROTECTED) != 0);
    /* FromDigits must NOT be Listable -- the first argument IS a list. */
    ASSERT((a & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("FromDigits");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "list of decimal digits") != NULL);
    ASSERT(strstr(def->docstring, "base b") != NULL);
    ASSERT(strstr(def->docstring, "string") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_FromDigits != NULL);
    ASSERT(strcmp(SYM_FromDigits, "FromDigits") == 0);
}

/* --- Memory-safety stress loop -------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix integer list / string / bignum / symbolic / error paths so
     * any mishandled mpz_t or args[] slot surfaces under valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("FromDigits[{5,1,2,8}]", "5128", 0);
        assert_eval_eq("FromDigits[\"1923\"]", "1923", 0);
        assert_eval_eq("FromDigits[{1,0,1,1,0,1,1}, 2]", "91", 0);
        assert_eval_eq("FromDigits[{a,b,c,d}, x]",
                       "d + c x + b x^2 + a x^3", 0);
        assert_eval_eq("FromDigits[IntegerDigits[2^100]]",
                       "1267650600228229401496703205376", 0);
        assert_eval_eq("FromDigits[{}]", "0", 0);
        assert_eval_eq("FromDigits[\"ABC\", 16]", "2748", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_doc_decimal_list);
    TEST(test_doc_binary_list);
    TEST(test_doc_decimal_string);
    TEST(test_doc_binary_string);
    TEST(test_doc_carrying);
    TEST(test_doc_string_base10_letters);
    TEST(test_doc_symbolic_polynomial);
    TEST(test_doc_single_chars);

    TEST(test_inverse_positive);
    TEST(test_inverse_negative_drops_sign);
    TEST(test_inverse_arbitrary_base);

    TEST(test_bignum_2_to_100);
    TEST(test_bignum_2_to_100_base2);
    TEST(test_bignum_10_to_30_plus);
    TEST(test_bignum_factorial);
    TEST(test_bignum_negative_digit_in_list);

    TEST(test_carry_two_digit);
    TEST(test_carry_three_digit);
    TEST(test_carry_large);

    TEST(test_symbolic_base);
    TEST(test_symbolic_digit_in_list);
    TEST(test_symbolic_both);
    TEST(test_real_base);
    TEST(test_rational_base);

    TEST(test_string_hex);
    TEST(test_string_hex_mixed_case);
    TEST(test_string_base_36);
    TEST(test_string_bignum);

    TEST(test_empty_list);
    TEST(test_empty_list_with_base);
    TEST(test_empty_string);
    TEST(test_zero_list);
    TEST(test_leading_zeros);
    TEST(test_single_digit_list);

    TEST(test_zero_args_argb);
    TEST(test_three_args_argb);
    TEST(test_numeric_first_arg_nlst);
    TEST(test_real_first_arg_nlst);
    TEST(test_bad_base_zero_list);
    TEST(test_bad_base_negative_string);
    TEST(test_invalid_string_char);
    TEST(test_symbolic_base_over_string_silent);
    TEST(test_symbolic_first_arg_silent);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All FromDigits tests passed!\n");
    return 0;
}
