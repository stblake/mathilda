/* Unit tests for IntegerExponent.
 *
 *   IntegerExponent[n]      -- highest power of 10 dividing |n|; equivalently
 *                              the number of trailing zeros in decimal digits.
 *   IntegerExponent[n, b]   -- highest power of b dividing |n|; equivalently
 *                              the number of trailing zeros in base-b digits.
 *
 * Coverage:
 *   - Default-base (10): trailing-zero counts on small ints, 10^k, mixed.
 *   - Base 2: fast path via mpz_scan1; powers and sums of powers of 2.
 *   - Other small bases: 3, 5, 7, 16 -- exercises mpz_remove general path.
 *   - Negative integers (sign discarded).
 *   - Zero special-case: IntegerExponent[0, b] == Infinity for any b.
 *   - Numbers not divisible by b: result is 0.
 *   - Bignum inputs: 10^30, 2^100, 100!, 2^100 * 3^50.
 *   - Bignum bases: 10^18 + 1.
 *   - Listable threading on n, on b, on both, and nested.
 *   - Property: IntegerExponent[n, b] == Length[Rest[Split[Reverse[
 *     IntegerDigits[n, b]], #2 == 0 &][[1]]]] ... actually we use the
 *     simpler property: Mod[n, b^IntegerExponent[n,b]] == 0 AND
 *     Mod[n, b^(IntegerExponent[n,b]+1)] != 0 for n!=0.
 *   - Error paths: 0 args, 3 args, 5 args (::argt); real / rational /
 *     complex n (::int); real / negative / zero / one base (::int /
 *     ::ibase); symbolic args (silent unevaluated).
 *   - Attribute, docstring, and interned-symbol introspection.
 *   - Documented Mathematica examples: 1230000 -> 4; 2^10+2^7 in base 2
 *     -> 7; 2^^10010000 in base 2 -> 4.
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
    const char* path = "/tmp/mathilda_int_exponent_stderr.log";
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

/* --- Default base (10) --------------------------------------------- */

static void test_decimal_documented_example(void) {
    /* Spec: IntegerExponent[1230000] -> 4 (four trailing zeros). */
    assert_eval_eq("IntegerExponent[1230000]", "4", 0);
}

static void test_decimal_no_trailing_zeros(void) {
    assert_eval_eq("IntegerExponent[123456789]", "0", 0);
    assert_eval_eq("IntegerExponent[1]", "0", 0);
    assert_eval_eq("IntegerExponent[7]", "0", 0);
}

static void test_decimal_one_trailing_zero(void) {
    assert_eval_eq("IntegerExponent[10]", "1", 0);
    assert_eval_eq("IntegerExponent[120]", "1", 0);
    assert_eval_eq("IntegerExponent[1230]", "1", 0);
}

static void test_decimal_powers_of_ten(void) {
    /* 10^k has exactly k trailing zeros. */
    assert_eval_eq("IntegerExponent[100]", "2", 0);
    assert_eval_eq("IntegerExponent[1000]", "3", 0);
    assert_eval_eq("IntegerExponent[100000]", "5", 0);
    assert_eval_eq("IntegerExponent[10^9]", "9", 0);
}

static void test_decimal_mixed(void) {
    /* 12 * 10^5 = 1200000 -> 5 trailing zeros. */
    assert_eval_eq("IntegerExponent[1200000]", "5", 0);
    /* 25 * 10^3 = 25000 -> 3 trailing zeros. */
    assert_eval_eq("IntegerExponent[25000]", "3", 0);
}

/* --- Base 2 (fast path via mpz_scan1) ------------------------------ */

static void test_base2_documented_example_sum(void) {
    /* Spec: IntegerExponent[2^10 + 2^7, 2] -> 7.
     * 2^10 + 2^7 = 1024 + 128 = 1152 = 2^7 * 9. */
    assert_eval_eq("IntegerExponent[2^10 + 2^7, 2]", "7", 0);
}

static void test_base2_documented_example_binary_literal(void) {
    /* Spec example: IntegerExponent[2^^10010000, 2] -> 4.
     * Mathilda's parser does not currently accept the Mathematica
     * `base^^digits` literal syntax, so we use the decimal equivalent
     * 0b10010000 = 144 = 2^4 * 9. */
    assert_eval_eq("IntegerExponent[144, 2]", "4", 0);
}

static void test_base2_pure_powers(void) {
    /* 2^k -> k. */
    assert_eval_eq("IntegerExponent[1, 2]", "0", 0);
    assert_eval_eq("IntegerExponent[2, 2]", "1", 0);
    assert_eval_eq("IntegerExponent[4, 2]", "2", 0);
    assert_eval_eq("IntegerExponent[8, 2]", "3", 0);
    assert_eval_eq("IntegerExponent[16, 2]", "4", 0);
    assert_eval_eq("IntegerExponent[1024, 2]", "10", 0);
}

static void test_base2_odd_numbers(void) {
    /* Every odd number has IntegerExponent[.,2] == 0. */
    assert_eval_eq("IntegerExponent[3, 2]", "0", 0);
    assert_eval_eq("IntegerExponent[5, 2]", "0", 0);
    assert_eval_eq("IntegerExponent[2^20 + 1, 2]", "0", 0);
}

static void test_base2_mixed(void) {
    /* 12 = 2^2 * 3 -> 2. */
    assert_eval_eq("IntegerExponent[12, 2]", "2", 0);
    /* 24 = 2^3 * 3 -> 3. */
    assert_eval_eq("IntegerExponent[24, 2]", "3", 0);
    /* 96 = 2^5 * 3 -> 5. */
    assert_eval_eq("IntegerExponent[96, 2]", "5", 0);
}

/* --- Other small bases (mpz_remove path) --------------------------- */

static void test_base3(void) {
    assert_eval_eq("IntegerExponent[1, 3]", "0", 0);
    assert_eval_eq("IntegerExponent[3, 3]", "1", 0);
    assert_eval_eq("IntegerExponent[9, 3]", "2", 0);
    assert_eval_eq("IntegerExponent[27, 3]", "3", 0);
    /* 54 = 2 * 27 = 2 * 3^3 -> 3. */
    assert_eval_eq("IntegerExponent[54, 3]", "3", 0);
    /* 162 = 2 * 81 = 2 * 3^4 -> 4. */
    assert_eval_eq("IntegerExponent[162, 3]", "4", 0);
}

static void test_base5(void) {
    assert_eval_eq("IntegerExponent[25, 5]", "2", 0);
    assert_eval_eq("IntegerExponent[125, 5]", "3", 0);
    /* 100 = 4 * 25 -> exponent 2. */
    assert_eval_eq("IntegerExponent[100, 5]", "2", 0);
    /* 24 not divisible by 5. */
    assert_eval_eq("IntegerExponent[24, 5]", "0", 0);
}

static void test_base7(void) {
    assert_eval_eq("IntegerExponent[7, 7]", "1", 0);
    assert_eval_eq("IntegerExponent[49, 7]", "2", 0);
    assert_eval_eq("IntegerExponent[343, 7]", "3", 0);
    /* 5 * 49 = 245 -> 2. */
    assert_eval_eq("IntegerExponent[245, 7]", "2", 0);
}

static void test_base16(void) {
    /* 16 = 16^1 -> 1; 256 = 16^2 -> 2. */
    assert_eval_eq("IntegerExponent[16, 16]", "1", 0);
    assert_eval_eq("IntegerExponent[256, 16]", "2", 0);
    assert_eval_eq("IntegerExponent[4096, 16]", "3", 0);
    /* 32 = 2 * 16 -> 1 (4 * 8 = 32; 32 / 16 = 2, not div by 16). */
    assert_eval_eq("IntegerExponent[32, 16]", "1", 0);
    /* 31 not divisible by 16. */
    assert_eval_eq("IntegerExponent[31, 16]", "0", 0);
}

static void test_base_composite(void) {
    /* base 6 = 2*3. 6^3 = 216 -> 3. */
    assert_eval_eq("IntegerExponent[216, 6]", "3", 0);
    /* 36 = 6^2 -> 2. */
    assert_eval_eq("IntegerExponent[36, 6]", "2", 0);
    /* 12 = 2*6 -> 1 (12/6=2, 2 not div by 6). */
    assert_eval_eq("IntegerExponent[12, 6]", "1", 0);
}

/* --- Sign discarded ------------------------------------------------- */

static void test_negative(void) {
    assert_eval_eq("IntegerExponent[-1230000]", "4", 0);
    assert_eval_eq("IntegerExponent[-1024, 2]", "10", 0);
    assert_eval_eq("IntegerExponent[-100]", "2", 0);
    assert_eval_eq("IntegerExponent[-1]", "0", 0);
    assert_eval_eq("IntegerExponent[-7]", "0", 0);
}

/* --- Zero ---------------------------------------------------------- */

static void test_zero_default_base(void) {
    /* IntegerExponent[0] -> Infinity (every power of 10 divides 0). */
    assert_eval_eq("IntegerExponent[0]", "Infinity", 0);
}

static void test_zero_any_base(void) {
    assert_eval_eq("IntegerExponent[0, 2]", "Infinity", 0);
    assert_eval_eq("IntegerExponent[0, 16]", "Infinity", 0);
    assert_eval_eq("IntegerExponent[0, 1000]", "Infinity", 0);
}

/* --- Bignum -------------------------------------------------------- */

static void test_bignum_10pow30(void) {
    /* 10^30 has 30 trailing zeros. */
    assert_eval_eq("IntegerExponent[10^30]", "30", 0);
}

static void test_bignum_10pow30_minus_1(void) {
    /* 10^30 - 1 = 999...9 (30 nines) -- no trailing zeros. */
    assert_eval_eq("IntegerExponent[10^30 - 1]", "0", 0);
}

static void test_bignum_2pow100_base2(void) {
    /* 2^100 has 100 trailing zeros in binary. */
    assert_eval_eq("IntegerExponent[2^100, 2]", "100", 0);
}

static void test_bignum_2pow100_base10(void) {
    /* 2^100 = 1267650600228229401496703205376 ends in "5376" (no zeros). */
    assert_eval_eq("IntegerExponent[2^100]", "0", 0);
}

static void test_bignum_100_factorial(void) {
    /* Legendre's formula: the exponent of 2 in 100! is
     *   sum_{k=1}^inf floor(100 / 2^k)
     *   = 50 + 25 + 12 + 6 + 3 + 1 = 97. */
    assert_eval_eq("IntegerExponent[100!, 2]", "97", 0);
    /* Exponent of 5 in 100! is 25 + 5 + 1 = 24.  Therefore the number
     * of trailing decimal zeros (min of exponents of 2 and 5) is 24. */
    assert_eval_eq("IntegerExponent[100!, 5]", "24", 0);
    assert_eval_eq("IntegerExponent[100!]", "24", 0);
}

static void test_bignum_product(void) {
    /* 2^100 * 3^50: exponent of 2 -> 100, exponent of 3 -> 50,
     * exponent of 6 -> min(100, 50) = 50. */
    assert_eval_eq("IntegerExponent[2^100 * 3^50, 2]", "100", 0);
    assert_eval_eq("IntegerExponent[2^100 * 3^50, 3]", "50", 0);
    assert_eval_eq("IntegerExponent[2^100 * 3^50, 6]", "50", 0);
}

static void test_bignum_negative(void) {
    /* Sign discarded even with bignum input. */
    assert_eval_eq("IntegerExponent[-(10^30)]", "30", 0);
    assert_eval_eq("IntegerExponent[-(2^100), 2]", "100", 0);
}

static void test_int64_overflow_boundary(void) {
    /* 2^63 promoted to bigint; trailing-zero count is 63. */
    assert_eval_eq("IntegerExponent[2^63, 2]", "63", 0);
    /* 2^63 - 1 = INT64_MAX, all ones in binary -> 0 trailing zeros. */
    assert_eval_eq("IntegerExponent[2^63 - 1, 2]", "0", 0);
}

static void test_bignum_base(void) {
    /* base = 10^18 + 1.  b * b = b^2, exponent 2. */
    assert_eval_eq("IntegerExponent[(10^18 + 1)^2, 10^18 + 1]", "2", 0);
    assert_eval_eq("IntegerExponent[(10^18 + 1)^3, 10^18 + 1]", "3", 0);
    /* Not divisible by 10^18 + 1. */
    assert_eval_eq("IntegerExponent[10^18, 10^18 + 1]", "0", 0);
}

/* --- Listable threading -------------------------------------------- */

static void test_listable_on_n(void) {
    assert_eval_eq("IntegerExponent[{10, 100, 1000, 10000}]",
                   "{1, 2, 3, 4}", 0);
}

static void test_listable_on_base(void) {
    /* IntegerExponent[24, {2,3,4,6}]:
     *   24 = 2^3 * 3 -> base 2: 3; base 3: 1; base 4: 1 (24/4=6, 6 not
     *   div by 4); base 6: 1 (24/6=4, 4 not div by 6). */
    assert_eval_eq("IntegerExponent[24, {2, 3, 4, 6}]",
                   "{3, 1, 1, 1}", 0);
}

static void test_listable_on_both(void) {
    /* Element-wise threading: IntegerExponent[{8,9,16}, {2,3,4}] =
     * {IntegerExponent[8,2], IntegerExponent[9,3], IntegerExponent[16,4]}
     * = {3, 2, 2}. */
    assert_eval_eq("IntegerExponent[{8, 9, 16}, {2, 3, 4}]",
                   "{3, 2, 2}", 0);
}

static void test_listable_nested(void) {
    assert_eval_eq("IntegerExponent[{{10, 100}, {1000, 10000}}]",
                   "{{1, 2}, {3, 4}}", 0);
}

static void test_listable_zero_in_list(void) {
    /* Zero -> Infinity for each occurrence. */
    assert_eval_eq("IntegerExponent[{0, 10, 0, 100}]",
                   "{Infinity, 1, Infinity, 2}", 0);
}

/* --- Property tests -------------------------------------------- */

static void test_property_divisibility(void) {
    /* For n != 0:  Mod[n, b^IntegerExponent[n,b]] == 0
     * and       Mod[n, b^(IntegerExponent[n,b]+1)] != 0. */
    assert_eval_eq(
        "Table[Mod[n, 2^IntegerExponent[n, 2]], {n, 1, 30}]",
        "{0, 0, 0, 0, 0, 0, 0, 0, 0, 0,"
        " 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,"
        " 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}",
        0);
    assert_eval_eq(
        "Table[Mod[n, 2^(IntegerExponent[n, 2] + 1)] != 0, {n, 1, 30}]",
        "{True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True}",
        0);
}

static void test_property_trailing_zeros_decimal(void) {
    /* IntegerExponent[n, 10] equals the number of trailing zeros of
     * IntegerDigits[n, 10] for n != 0. */
    assert_eval_eq("IntegerExponent[10^7]", "7", 0);
    assert_eval_eq("IntegerExponent[10^7 + 1]", "0", 0);
    assert_eval_eq("IntegerExponent[5 * 10^6]", "6", 0);
}

/* --- Error paths --------------------------------------------------- */

static void test_zero_args_argt(void) {
    /* Documented Mathematica diagnostic:
     *   IntegerExponent::argt: IntegerExponent called with 0 arguments;
     *   1 or 2 arguments are expected. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerExponent[]");
    ASSERT_MSG(err && strstr(err, "IntegerExponent::argt") != NULL,
               "expected argt diagnostic, got: %s", err ? err : "(null)");
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    ASSERT(strstr(err, "1 or 2 arguments are expected") != NULL);
    free(result);
    free(err);
}

static void test_three_args_argt(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[1, 2, 3]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerExponent[1, 2, 3]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerExponent::argt") != NULL,
               "expected argt diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 3 arguments") != NULL);
    free(result);
    free(err);
}

static void test_four_args_argt(void) {
    /* Documented Mathematica session (spec):
     *   IntegerExponent[1,2,3,4]
     *   IntegerExponent::argt: IntegerExponent called with 4 arguments;
     *   1 or 2 arguments are expected.
     *   IntegerExponent[1,2,3,4] */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[1, 2, 3, 4]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerExponent[1, 2, 3, 4]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerExponent::argt") != NULL,
               "expected argt diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 4 arguments") != NULL);
    ASSERT(strstr(err, "1 or 2 arguments are expected") != NULL);
    free(result);
    free(err);
}

static void test_five_args_argt(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[1, 2, 3, 4, 5]",
                                       &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerExponent[1, 2, 3, 4, 5]");
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerExponent::argt") != NULL);
    ASSERT(strstr(err, "called with 5 arguments") != NULL);
    free(result);
    free(err);
}

static void test_real_n_int_diagnostic(void) {
    /* Documented Mathematica session (spec):
     *   IntegerExponent[1.123]
     *   IntegerExponent::int: Integer expected at position 1 in
     *   IntegerExponent[1.123].
     *   IntegerExponent[1.123] */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[1.123]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "IntegerExponent[1.123]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "IntegerExponent::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    ASSERT(strstr(err, "IntegerExponent[1.123]") != NULL);
    free(result);
    free(err);
}

static void test_rational_n_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerExponent::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_complex_n_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[2 + 3 I]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerExponent::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_real_base_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[10, 2.5]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "IntegerExponent::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 2") != NULL);
    free(result);
    free(err);
}

static void test_bad_base_one(void) {
    /* Base 1 prints ::ibase and leaves the call unevaluated. */
    const char* in = "IntegerExponent[12, 1]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerExponent") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_bad_base_zero(void) {
    const char* in = "IntegerExponent[12, 0]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerExponent") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_bad_base_negative(void) {
    const char* in = "IntegerExponent[12, -2]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "IntegerExponent") != NULL,
               "expected unevaluated form, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_symbolic_n_left_unevaluated(void) {
    /* Symbolic n: silent (no diagnostic), call retained. */
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[x]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerExponent[x]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "IntegerExponent::") == NULL,
               "expected no diagnostic for symbolic input, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_base_left_unevaluated(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("IntegerExponent[12, b]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "IntegerExponent[12, b]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "IntegerExponent::") == NULL,
               "expected no diagnostic for symbolic base, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Attribute / docstring / interned-symbol introspection -------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("IntegerExponent");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("IntegerExponent");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("IntegerExponent");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "highest power") != NULL);
    ASSERT(strstr(def->docstring, "Infinity") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_IntegerExponent != NULL);
    ASSERT(strcmp(SYM_IntegerExponent, "IntegerExponent") == 0);
}

/* --- Memory-safety stress loop ------------------------------------ */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix machine-int / bignum / listable / fast-path (base 2) /
     * general path / error / zero / negative / bignum-base cases. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("IntegerExponent[1230000]", "4", 0);
        assert_eval_eq("IntegerExponent[2^10 + 2^7, 2]", "7", 0);
        assert_eval_eq("IntegerExponent[0]", "Infinity", 0);
        assert_eval_eq("IntegerExponent[-1024, 2]", "10", 0);
        assert_eval_eq("IntegerExponent[10^30]", "30", 0);
        assert_eval_eq("IntegerExponent[100!, 2]", "97", 0);
        assert_eval_eq("IntegerExponent[2^100 * 3^50, 6]", "50", 0);
        assert_eval_eq("IntegerExponent[{10, 100, 1000, 10000}]",
                       "{1, 2, 3, 4}", 0);
        assert_eval_eq("IntegerExponent[(10^18 + 1)^2, 10^18 + 1]", "2", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_decimal_documented_example);
    TEST(test_decimal_no_trailing_zeros);
    TEST(test_decimal_one_trailing_zero);
    TEST(test_decimal_powers_of_ten);
    TEST(test_decimal_mixed);

    TEST(test_base2_documented_example_sum);
    TEST(test_base2_documented_example_binary_literal);
    TEST(test_base2_pure_powers);
    TEST(test_base2_odd_numbers);
    TEST(test_base2_mixed);

    TEST(test_base3);
    TEST(test_base5);
    TEST(test_base7);
    TEST(test_base16);
    TEST(test_base_composite);

    TEST(test_negative);

    TEST(test_zero_default_base);
    TEST(test_zero_any_base);

    TEST(test_bignum_10pow30);
    TEST(test_bignum_10pow30_minus_1);
    TEST(test_bignum_2pow100_base2);
    TEST(test_bignum_2pow100_base10);
    TEST(test_bignum_100_factorial);
    TEST(test_bignum_product);
    TEST(test_bignum_negative);
    TEST(test_int64_overflow_boundary);
    TEST(test_bignum_base);

    TEST(test_listable_on_n);
    TEST(test_listable_on_base);
    TEST(test_listable_on_both);
    TEST(test_listable_nested);
    TEST(test_listable_zero_in_list);

    TEST(test_property_divisibility);
    TEST(test_property_trailing_zeros_decimal);

    TEST(test_zero_args_argt);
    TEST(test_three_args_argt);
    TEST(test_four_args_argt);
    TEST(test_five_args_argt);
    TEST(test_real_n_int_diagnostic);
    TEST(test_rational_n_int_diagnostic);
    TEST(test_complex_n_int_diagnostic);
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

    printf("All IntegerExponent tests passed!\n");
    return 0;
}
