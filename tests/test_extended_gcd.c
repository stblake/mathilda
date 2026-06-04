/* Unit tests for ExtendedGCD.
 *
 *   ExtendedGCD[n1, n2, ...] -> {g, {r1, r2, ...}}
 *     where g == GCD[n1, ...] and g == r1 n1 + r2 n2 + ....
 *
 * Coverage:
 *   - Documented Mathematica examples ({1,{-1,1}}, {3,{-2,1,0}}, ...).
 *   - The leading element equals GCD of the arguments.
 *   - The Bezout identity g == Sum[ri ni] holds (checked via assignment).
 *   - Edge cases: ExtendedGCD[] -> {0,{}}; single arg; zeros; negatives.
 *   - BigInt arguments and mixed machine/BigInt arguments (the cofactors
 *     and gcd are computed in GMP and demoted back to machine integers
 *     when they fit).
 *   - Listable threading on first arg, second arg, and both.
 *   - Error paths: inexact (Real) argument -> ExtendedGCD::exact; exact
 *     non-integer (Rational) argument -> ExtendedGCD::egcd; symbolic
 *     argument -> silent, left unevaluated.
 *   - Attribute (Protected + Listable, NOT Flat/Orderless), docstring,
 *     and interned-symbol introspection.
 *   - Repeated-evaluation stress loop to catch double-frees / leaks under
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
    const char* path = "/tmp/mathilda_extended_gcd_stderr.log";
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

/* --- Documented examples ------------------------------------------- */

static void test_doc_two_args(void) {
    /* Spec: ExtendedGCD[2, 3] -> {1, {-1, 1}}; -1*2 + 1*3 == 1. */
    assert_eval_eq("ExtendedGCD[2, 3]", "{1, {-1, 1}}", 0);
}

static void test_doc_three_args(void) {
    /* Spec: ExtendedGCD[6, 15, 30] -> {3, {-2, 1, 0}}. */
    assert_eval_eq("ExtendedGCD[6, 15, 30]", "{3, {-2, 1, 0}}", 0);
}

static void test_doc_six_twentyone(void) {
    /* Spec: ExtendedGCD[6, 21] -> {3, {-3, 1}}; GCD[6,21] == 3. */
    assert_eval_eq("ExtendedGCD[6, 21]", "{3, {-3, 1}}", 0);
    assert_eval_eq("GCD[6, 21]", "3", 0);
}

static void test_doc_threading_example(void) {
    /* Spec: ExtendedGCD[3, {5, 15}] -> {{1, {2, -1}}, {3, {1, 0}}}. */
    assert_eval_eq("ExtendedGCD[3, {5, 15}]",
                   "{{1, {2, -1}}, {3, {1, 0}}}", 0);
}

/* --- More small integer cases (authoritative values) --------------- */

static void test_small_cases(void) {
    assert_eval_eq("ExtendedGCD[100, 35]", "{5, {-1, 3}}", 0);
    assert_eval_eq("ExtendedGCD[35, 100]", "{5, {3, -1}}", 0);
    assert_eval_eq("ExtendedGCD[12, 18]", "{6, {-1, 1}}", 0);
    assert_eval_eq("ExtendedGCD[17, 5]", "{1, {-2, 7}}", 0);
    assert_eval_eq("ExtendedGCD[8, 12, 20]", "{4, {-1, 1, 0}}", 0);
}

/* --- The leading element is the GCD -------------------------------- */

static void test_first_element_is_gcd(void) {
    assert_eval_eq("ExtendedGCD[6, 21][[1]]", "3", 0);
    assert_eval_eq("ExtendedGCD[6, 21][[1]] == GCD[6, 21]", "True", 0);
    assert_eval_eq("ExtendedGCD[12, 18, 30][[1]] == GCD[12, 18, 30]",
                   "True", 0);
}

/* --- Bezout identity g == Sum[ri ni] ------------------------------- */

static void test_bezout_two(void) {
    /* {g,{a,b}} = ExtendedGCD[2,3]; verify 2 a + 3 b == g. */
    assert_eval_eq("{g, {a, b}} = ExtendedGCD[2, 3]; 2 a + 3 b == g",
                   "True", 0);
}

static void test_bezout_three(void) {
    assert_eval_eq(
        "{g, {a, b, c}} = ExtendedGCD[6, 15, 30]; 6 a + 15 b + 30 c == g",
        "True", 0);
}

static void test_bezout_dot_property(void) {
    /* For r = ExtendedGCD[ns], r[[1]] == r[[2]] . ns. */
    assert_eval_eq(
        "r = ExtendedGCD[100, 35, 14]; r[[1]] == r[[2]] . {100, 35, 14}",
        "True", 0);
}

/* --- Edge cases ---------------------------------------------------- */

static void test_empty(void) {
    /* Spec: ExtendedGCD[] -> {0, {}}. */
    assert_eval_eq("ExtendedGCD[]", "{0, {}}", 0);
}

static void test_single_arg(void) {
    assert_eval_eq("ExtendedGCD[12]", "{12, {1}}", 0);
    assert_eval_eq("ExtendedGCD[1]", "{1, {1}}", 0);
}

static void test_single_arg_negative(void) {
    /* g is non-negative; cofactor flips sign: -1 * -6 == 6. */
    assert_eval_eq("ExtendedGCD[-6]", "{6, {-1}}", 0);
}

static void test_zeros(void) {
    assert_eval_eq("ExtendedGCD[0, 5]", "{5, {0, 1}}", 0);
    assert_eval_eq("ExtendedGCD[5, 0]", "{5, {1, 0}}", 0);
    assert_eval_eq("ExtendedGCD[0, 0]", "{0, {0, 0}}", 0);
}

static void test_negatives(void) {
    /* GCD is normalised non-negative; identity still holds. */
    assert_eval_eq("ExtendedGCD[-4, 6]", "{2, {1, 1}}", 0);
    assert_eval_eq("ExtendedGCD[4, -6]", "{2, {-1, -1}}", 0);
    assert_eval_eq(
        "{g, {a, b}} = ExtendedGCD[-4, 6]; (-4) a + 6 b == g", "True", 0);
}

/* --- BigInt and mixed machine/BigInt arguments --------------------- */

static void test_bigint_args(void) {
    /* gcd(2^100, 2^100 + 2) = gcd(2^100, 2) = 2. */
    assert_eval_eq("ExtendedGCD[2^100, 2^100 + 2]", "{2, {-1, 1}}", 0);
    /* gcd(10^20, 10^20 + 5^20) = 5^20 = 95367431640625. */
    assert_eval_eq("ExtendedGCD[10^20, 10^20 + 5^20]",
                   "{95367431640625, {-1, 1}}", 0);
}

static void test_bigint_bezout(void) {
    assert_eval_eq(
        "{g, {a, b}} = ExtendedGCD[2^100, 2^100 + 2]; "
        "a 2^100 + b (2^100 + 2) == g",
        "True", 0);
    assert_eval_eq(
        "{g, {a, b}} = ExtendedGCD[10^20, 10^20 + 5^20]; "
        "a 10^20 + b (10^20 + 5^20) == g",
        "True", 0);
}

static void test_mixed_machine_bigint(void) {
    /* One machine int, one bigint: result demotes to machine ints. */
    assert_eval_eq("ExtendedGCD[6, 10^30][[1]] == GCD[6, 10^30]", "True", 0);
    assert_eval_eq(
        "{g, {a, b}} = ExtendedGCD[6, 10^30]; 6 a + 10^30 b == g", "True", 0);
}

/* --- Listable threading -------------------------------------------- */

static void test_listable_first_arg(void) {
    /* Broadcast scalar second arg. */
    assert_eval_eq("ExtendedGCD[{12, 18}, {8, 20}]",
                   "{{4, {1, -1}}, {2, {-1, 1}}}", 0);
}

static void test_listable_scalar_broadcast(void) {
    assert_eval_eq("ExtendedGCD[3, {5, 15}]",
                   "{{1, {2, -1}}, {3, {1, 0}}}", 0);
}

static void test_listable_nested(void) {
    assert_eval_eq("ExtendedGCD[{{2, 3}}, {{4, 9}}]",
                   "{{{2, {1, 0}}, {3, {1, 0}}}}", 0);
}

/* --- Error paths --------------------------------------------------- */

static void test_inexact_real_diagnostic(void) {
    /* Spec: ExtendedGCD[4.5, 3.1]
     *   ExtendedGCD::exact: Argument 4.5 in ExtendedGCD[4.5, 3.1] is
     *   not an exact number.  Call left unevaluated. */
    char* result = NULL;
    char* err = eval_capturing_stderr("ExtendedGCD[4.5, 3.1]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "ExtendedGCD") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "ExtendedGCD::exact") != NULL,
               "expected exact diagnostic, got: %s", err);
    ASSERT(strstr(err, "is not an exact number") != NULL);
    free(result);
    free(err);
}

static void test_inexact_single_real(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("ExtendedGCD[6, 2.5]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "ExtendedGCD::exact") != NULL);
    free(result);
    free(err);
}

static void test_rational_egcd_diagnostic(void) {
    /* Spec: ExtendedGCD[4/5, 3/2]
     *   ExtendedGCD::egcd: Arguments in ExtendedGCD[4/5, 3/2] should be
     *   integers.  Call left unevaluated. */
    char* result = NULL;
    char* err = eval_capturing_stderr("ExtendedGCD[4/5, 3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "ExtendedGCD") != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "ExtendedGCD::egcd") != NULL,
               "expected egcd diagnostic, got: %s", err);
    ASSERT(strstr(err, "should be integers") != NULL);
    free(result);
    free(err);
}

static void test_rational_mixed_egcd(void) {
    /* An integer together with a rational still triggers ::egcd. */
    char* result = NULL;
    char* err = eval_capturing_stderr("ExtendedGCD[6, 3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "ExtendedGCD::egcd") != NULL);
    free(result);
    free(err);
}

static void test_symbolic_silent(void) {
    /* Symbolic argument: no diagnostic, call retained. */
    char* result = NULL;
    char* err = eval_capturing_stderr("ExtendedGCD[x, 2]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "ExtendedGCD[x, 2]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "ExtendedGCD::") == NULL,
               "expected no diagnostic for symbolic input, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

static void test_symbolic_only(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("ExtendedGCD[x, y]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "ExtendedGCD[x, y]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "ExtendedGCD::") == NULL,
               "expected no diagnostic for symbolic input, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Attribute / docstring / interned-symbol introspection -------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("ExtendedGCD");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("ExtendedGCD");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
    /* ExtendedGCD is order-sensitive: must NOT be Flat / Orderless /
     * OneIdentity (unlike GCD), or the cofactor list would be scrambled. */
    ASSERT((a & ATTR_FLAT) == 0);
    ASSERT((a & ATTR_ORDERLESS) == 0);
    ASSERT((a & ATTR_ONEIDENTITY) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("ExtendedGCD");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "extended GCD") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_ExtendedGCD != NULL);
    ASSERT(strcmp(SYM_ExtendedGCD, "ExtendedGCD") == 0);
}

/* --- Memory-safety stress loop ------------------------------------ */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix machine-int / bignum / listable / single / zero / negative /
     * error cases to catch double-frees and leaks under valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("ExtendedGCD[2, 3]", "{1, {-1, 1}}", 0);
        assert_eval_eq("ExtendedGCD[6, 15, 30]", "{3, {-2, 1, 0}}", 0);
        assert_eval_eq("ExtendedGCD[]", "{0, {}}", 0);
        assert_eval_eq("ExtendedGCD[-6]", "{6, {-1}}", 0);
        assert_eval_eq("ExtendedGCD[2^100, 2^100 + 2]", "{2, {-1, 1}}", 0);
        assert_eval_eq("ExtendedGCD[3, {5, 15}]",
                       "{{1, {2, -1}}, {3, {1, 0}}}", 0);
        char* result = NULL;
        char* err = eval_capturing_stderr("ExtendedGCD[4.5, 3.1]", &result);
        free(result);
        free(err);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_doc_two_args);
    TEST(test_doc_three_args);
    TEST(test_doc_six_twentyone);
    TEST(test_doc_threading_example);

    TEST(test_small_cases);

    TEST(test_first_element_is_gcd);

    TEST(test_bezout_two);
    TEST(test_bezout_three);
    TEST(test_bezout_dot_property);

    TEST(test_empty);
    TEST(test_single_arg);
    TEST(test_single_arg_negative);
    TEST(test_zeros);
    TEST(test_negatives);

    TEST(test_bigint_args);
    TEST(test_bigint_bezout);
    TEST(test_mixed_machine_bigint);

    TEST(test_listable_first_arg);
    TEST(test_listable_scalar_broadcast);
    TEST(test_listable_nested);

    TEST(test_inexact_real_diagnostic);
    TEST(test_inexact_single_real);
    TEST(test_rational_egcd_diagnostic);
    TEST(test_rational_mixed_egcd);
    TEST(test_symbolic_silent);
    TEST(test_symbolic_only);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All ExtendedGCD tests passed!\n");
    return 0;
}
