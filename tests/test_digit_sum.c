/* Unit tests for DigitSum.
 *
 *   DigitSum[n]      sum of the base-10 digits of |n|.
 *   DigitSum[n, b]   sum of the base-b digits of |n|.
 *
 * Coverage:
 *   - Basic decimal cases: DigitSum[1234] == 10
 *   - Zero: DigitSum[0] -> 0
 *   - Negative (sign discarded): DigitSum[-1234] == DigitSum[1234]
 *   - Two-arg base form: DigitSum[255, 16] == 30 (F+F = 15+15)
 *   - Cross-check invariant: DigitSum[n] == Total[IntegerDigits[n]]
 *   - Bignum inputs (10^30, 100!)
 *   - Arity errors: 0-arg and 3+ args leave unevaluated
 *   - Non-integer n (Real, Rational) left unevaluated or diagnostics
 *   - Symbolic args left unevaluated silently
 *   - Bad base (< 2, non-integer) left unevaluated with diagnostic
 *   - Attributes: PROTECTED | LISTABLE | NUMERICFUNCTION
 *   - Docstring registered and contains expected keywords
 *   - SYM_DigitSum interned symbol pointer
 *   - Stress loop for valgrind hygiene
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
#include <unistd.h>

/* Capture stderr while `input` is parsed + evaluated.  Returns the
 * collected stderr text as a heap string (caller frees) and writes the
 * printed result into *out_result_str (also heap-allocated). */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_digit_sum_stderr.log";

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

/* --- Basic decimal cases ------------------------------------------ */

static void test_basic_1234(void) {
    /* 1+2+3+4 = 10 */
    assert_eval_eq("DigitSum[1234]", "10", 0);
}

static void test_basic_single_digit(void) {
    assert_eval_eq("DigitSum[7]", "7", 0);
}

static void test_basic_repunit(void) {
    /* 111 -> 3 */
    assert_eval_eq("DigitSum[111]", "3", 0);
}

static void test_basic_9999(void) {
    /* 9+9+9+9 = 36 */
    assert_eval_eq("DigitSum[9999]", "36", 0);
}

static void test_basic_with_zeros(void) {
    /* 1001 -> 1+0+0+1 = 2 */
    assert_eval_eq("DigitSum[1001]", "2", 0);
}

/* --- Zero ---------------------------------------------------------- */

static void test_zero(void) {
    /* Zero has no significant digits; sum is 0. */
    assert_eval_eq("DigitSum[0]", "0", 0);
}

static void test_zero_base2(void) {
    assert_eval_eq("DigitSum[0, 2]", "0", 0);
}

/* --- Negative (sign discarded) ------------------------------------- */

static void test_negative_small(void) {
    assert_eval_eq("DigitSum[-1234]", "10", 0);
}

static void test_negative_matches_positive(void) {
    assert_eval_eq("DigitSum[-9999]", "36", 0);
}

/* --- Two-arg base form --------------------------------------------- */

static void test_base2_popcount(void) {
    /* 1234 in binary = 10011010010; popcount = 5 */
    assert_eval_eq("DigitSum[1234, 2]", "5", 0);
}

static void test_base16_ff(void) {
    /* 255 = 0xFF; digits are {15, 15}; sum = 30 */
    assert_eval_eq("DigitSum[255, 16]", "30", 0);
}

static void test_base8(void) {
    /* 64 = 100 in base 8; digits {1,0,0}; sum = 1 */
    assert_eval_eq("DigitSum[64, 8]", "1", 0);
}

static void test_base10_explicit(void) {
    /* Explicit base 10 matches default. */
    assert_eval_eq("DigitSum[1234, 10]", "10", 0);
}

/* --- Cross-check invariant DigitSum[n] == Total[IntegerDigits[n]] -- */

static void test_cross_check_small(void) {
    assert_eval_eq(
        "Table[DigitSum[n] == Total[IntegerDigits[n]], {n, 1, 20}]",
        "{True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True}",
        0);
}

static void test_cross_check_base3(void) {
    assert_eval_eq(
        "Table[DigitSum[n, 3] == Total[IntegerDigits[n, 3]], {n, 1, 20}]",
        "{True, True, True, True, True, True, True, True, True, True,"
        " True, True, True, True, True, True, True, True, True, True}",
        0);
}

/* --- Bignum -------------------------------------------------------- */

static void test_bignum_10pow30(void) {
    /* 10^30 = 1 followed by 30 zeros; digit sum = 1 */
    assert_eval_eq("DigitSum[10^30]", "1", 0);
}

static void test_bignum_2pow100_base2(void) {
    /* 2^100 in binary is 1 followed by 100 zeros; popcount = 1 */
    assert_eval_eq("DigitSum[2^100, 2]", "1", 0);
}

static void test_bignum_factorial(void) {
    /* DigitSum[100!] == Total[IntegerDigits[100!]] (sum of digit values) */
    assert_eval_eq(
        "DigitSum[100!] == Total[IntegerDigits[100!]]",
        "True", 0);
}

static void test_bignum_negative(void) {
    assert_eval_eq("DigitSum[-(10^30)]", "1", 0);
}

/* --- Arity errors (leave unevaluated) ------------------------------ */

static void test_zero_args(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitSum[]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitSum[]") != NULL);
    ASSERT_MSG(err && strstr(err, "DigitSum::argb") != NULL,
               "expected argb diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_three_args(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitSum[5, 2, 1]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitSum[") != NULL);
    ASSERT_MSG(err && strstr(err, "DigitSum::argb") != NULL,
               "expected argb diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Non-integer n ------------------------------------------------- */

static void test_real_n_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitSum[1.5]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitSum[1.5]") != NULL);
    ASSERT_MSG(err && strstr(err, "DigitSum::int") != NULL,
               "expected int diagnostic, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_rational_n_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitSum[3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err && strstr(err, "DigitSum::int") != NULL);
    free(result);
    free(err);
}

/* --- Symbolic args (silent unevaluated) ---------------------------- */

static void test_symbolic_n_silent(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitSum[x]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitSum[x]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "DigitSum::") == NULL,
               "expected no diagnostic for symbolic n, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_base_silent(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitSum[10, b]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "DigitSum[10, b]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "DigitSum::") == NULL,
               "expected no diagnostic for symbolic base, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Bad base ------------------------------------------------------ */

static void test_base_one_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitSum[5, 1]", &result);
    ASSERT(result != NULL);
    ASSERT(err && strstr(err, "DigitSum::base") != NULL);
    free(result);
    free(err);
}

static void test_base_zero_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitSum[5, 0]", &result);
    ASSERT(result != NULL);
    ASSERT(err && strstr(err, "DigitSum::base") != NULL);
    free(result);
    free(err);
}

static void test_rational_base_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("DigitSum[5, 3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err && strstr(err, "DigitSum::base") != NULL);
    free(result);
    free(err);
}

/* --- Listable threading -------------------------------------------- */

static void test_listable_threads(void) {
    /* Listable attribute threads over a list of integers. */
    assert_eval_eq("DigitSum[{1234, 0, 99}]", "{10, 0, 18}", 0);
}

/* --- Attributes / docstring / interned symbol ---------------------- */

static void test_attributes(void) {
    uint32_t a = get_attributes("DigitSum");
    ASSERT((a & ATTR_PROTECTED)        != 0);
    ASSERT((a & ATTR_LISTABLE)         != 0);
    ASSERT((a & ATTR_NUMERICFUNCTION)  != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("DigitSum");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "sum") != NULL);
    ASSERT(strstr(def->docstring, "DigitSum[n, b]") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_DigitSum != NULL);
    ASSERT(strcmp(SYM_DigitSum, "DigitSum") == 0);
}

/* --- Memory-safety stress loop ------------------------------------- */

static void test_stress_loop(void) {
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("DigitSum[1234]",    "10", 0);
        assert_eval_eq("DigitSum[0]",       "0",  0);
        assert_eval_eq("DigitSum[-9999]",   "36", 0);
        assert_eval_eq("DigitSum[255, 16]", "30", 0);
        assert_eval_eq("DigitSum[1234, 2]", "5",  0);
        assert_eval_eq("DigitSum[10^30]",   "1",  0);
        assert_eval_eq("DigitSum[{1, 99}]", "{1, 18}", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic_1234);
    TEST(test_basic_single_digit);
    TEST(test_basic_repunit);
    TEST(test_basic_9999);
    TEST(test_basic_with_zeros);

    TEST(test_zero);
    TEST(test_zero_base2);

    TEST(test_negative_small);
    TEST(test_negative_matches_positive);

    TEST(test_base2_popcount);
    TEST(test_base16_ff);
    TEST(test_base8);
    TEST(test_base10_explicit);

    TEST(test_cross_check_small);
    TEST(test_cross_check_base3);

    TEST(test_bignum_10pow30);
    TEST(test_bignum_2pow100_base2);
    TEST(test_bignum_factorial);
    TEST(test_bignum_negative);

    TEST(test_zero_args);
    TEST(test_three_args);

    TEST(test_real_n_diagnostic);
    TEST(test_rational_n_diagnostic);

    TEST(test_symbolic_n_silent);
    TEST(test_symbolic_base_silent);

    TEST(test_base_one_diagnostic);
    TEST(test_base_zero_diagnostic);
    TEST(test_rational_base_diagnostic);

    TEST(test_listable_threads);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_stress_loop);

    printf("All DigitSum tests passed!\n");
    return 0;
}
