/* Unit tests for CharacteristicPolynomial.
 *
 *   CharacteristicPolynomial[m, x]      == Det[m - x I]
 *   CharacteristicPolynomial[{m, a}, x] == Det[m - x a]
 *
 * Covers integer / rational / machine-real / complex / symbolic entries, both
 * the ordinary and generalized forms, degree-drop (infinite generalized
 * eigenvalue), numeric-variable substitution, arity errors, and a repeated-call
 * leak tripwire.
 *
 * The comparison helpers use the ASSERT_STR_EQ macro (hard exit(1)) rather than
 * test_utils.h's assert_eval_eq, whose libc assert() is compiled out under the
 * Release build's -DNDEBUG. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Evaluate `input`, compare its OutputForm printout to `expected`. */
static void check_out(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* r = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0)
        fprintf(stderr, "  input: %s\n", input);
    ASSERT_STR_EQ(s, expected);
    printf("  PASS: %s -> %s\n", input, s);
    free(s);
    expr_free(r);
}

/* Evaluate `input`, compare its FullForm printout to `expected`. */
static void check_ff(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* r = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string_fullform(r);
    if (strcmp(s, expected) != 0)
        fprintf(stderr, "  input: %s\n", input);
    ASSERT_STR_EQ(s, expected);
    printf("  PASS: %s -> %s\n", input, s);
    free(s);
    expr_free(r);
}

/* ---- Ordinary case: exact entries ---- */

static void test_integer_2x2(void) {
    check_out("CharacteristicPolynomial[{{1,2},{3,4}},x]", "-2 - 5 x + x^2");
    check_ff("CharacteristicPolynomial[{{1,2},{3,4}},x]",
             "Plus[-2, Times[-5, x], Power[x, 2]]");
}

static void test_symbolic_2x2(void) {
    check_out("CharacteristicPolynomial[{{a,b},{c,d}},x]",
              "-b c + a d - a x - d x + x^2");
    /* Strongest check: equals the definition Det[m - x I], expanded. */
    check_out("Simplify[CharacteristicPolynomial[{{a,b},{c,d}},x]"
              " - Expand[Det[{{a,b},{c,d}} - x*IdentityMatrix[2]]]]", "0");
}

static void test_identity_3x3(void) {
    /* (1 - L)^3 = 1 - 3 L + 3 L^2 - L^3 ; odd n carries a leading -L^3. */
    check_out("CharacteristicPolynomial[IdentityMatrix[3],L]",
              "1 - 3 L + 3 L^2 - L^3");
}

static void test_zero_matrix(void) {
    check_out("CharacteristicPolynomial[ConstantArray[0,{4,4}],x]", "x^4");
}

static void test_1x1(void) {
    check_out("CharacteristicPolynomial[{{a}},x]", "a - x");
    check_ff("CharacteristicPolynomial[{{a}},x]", "Plus[a, Times[-1, x]]");
    check_out("CharacteristicPolynomial[{{7}},x]", "7 - x");
}

static void test_rational_3x3(void) {
    check_out("CharacteristicPolynomial[{{1/3,1/2,3/5},{1/2,4/5,1},{3/5,1,9/7}},x]",
              "1/10500 - 239/2100 x + 254/105 x^2 - x^3");
}

/* ---- Ordinary case: inexact entries (checked by degree + Chop identity) ---- */

static void test_machine_3x3(void) {
    const char* m = "{{1.1,2.2,3.25},{0.76,4.6,5},{0.1,0.1,6.1}}";
    char buf[512];
    snprintf(buf, sizeof buf,
             "Exponent[CharacteristicPolynomial[%s,x],x]", m);
    check_out(buf, "3");
    /* Constant term ~= 19.9688 -> rounds to 20. */
    snprintf(buf, sizeof buf,
             "Round[CharacteristicPolynomial[%s,x] /. x -> 0]", m);
    check_out(buf, "20");
    /* Agrees with Det[m - x I] to machine precision (Chop clears round-off). */
    snprintf(buf, sizeof buf,
             "Chop[CharacteristicPolynomial[%s,x]"
             " - Expand[Det[%s - x*IdentityMatrix[3]]]]", m, m);
    check_out(buf, "0");
}

static void test_complex_3x3(void) {
    check_out("Exponent[CharacteristicPolynomial["
              "{{1.2+I,3-2 I,3 Pi},{-0.2,5I,2},{1,2.3,E}},x],x]", "3");
}

/* ---- Generalized case: CharacteristicPolynomial[{m, a}, x] == Det[m - x a] ---- */

static void test_generalized_integer(void) {
    check_out("CharacteristicPolynomial[{{{1,2},{5,4}},{{4,3},{6,4}}},x]",
              "-6 + 7 x - 2 x^2");
    check_ff("CharacteristicPolynomial[{{{1,2},{5,4}},{{4,3},{6,4}}},x]",
             "Plus[-6, Times[7, x], Times[-2, Power[x, 2]]]");
    check_out("Simplify[CharacteristicPolynomial[{{{1,2},{5,4}},{{4,3},{6,4}}},x]"
              " - Expand[Det[{{1,2},{5,4}} - x*{{4,3},{6,4}}]]]", "0");
}

static void test_generalized_machine(void) {
    check_out("Exponent[CharacteristicPolynomial["
              "{{{1.,1.5,2.},{3.1,2.,2.9},{3.,2.,1.}},"
              "{{1.3,.5,1.1},{0.,1.5,2.3},{1.,0.,1.}}},x],x]", "3");
}

static void test_generalized_degree_drop(void) {
    /* An infinite generalized eigenvalue drops the leading x^3 term:
     * det here is quadratic even though the matrices are 3x3. */
    check_out("CharacteristicPolynomial[{{{1,1,1},{1,0,1},{0,0,1}},"
              "{{0,1,1},{0,1,1},{1,0,0}}},x]", "-1 - x + x^2");
}

static void test_generalized_symbolic(void) {
    /* Matrices carry a free symbol x; the polynomial variable is y. */
    check_out("CharacteristicPolynomial[{{{x,1+x},{1-x,x}},{{1,1},{1,2x}}},y]",
              "-1 + 2 x^2 + 2 y - x y - 2 x^2 y - y^2 + 2 x y^2");
}

/* ---- Packed / NDArray input (CharacteristicPolynomial is on the AWARE list) ---- */

static void test_packed_and_ndarray(void) {
    /* A large integer matrix auto-packs; faddeev pack_unpacks it internally. */
    check_out("Exponent[CharacteristicPolynomial[IdentityMatrix[20],x],x]", "20");
    check_out("CharacteristicPolynomial[2 IdentityMatrix[3],x]",
              "8 - 12 x + 6 x^2 - x^3");
    /* A visible NDArray argument is delisted and re-evaluated. */
    check_out("Exponent[CharacteristicPolynomial[NDArray[{{1,2},{3,4}}],x],x]",
              "2");
}

/* ---- Variable forms ---- */

static void test_numeric_variable(void) {
    /* Second argument a number: evaluate the polynomial at that point. */
    check_out("CharacteristicPolynomial[{{1,2},{3,4}},5]", "-2");
    check_out("CharacteristicPolynomial[{{1,2},{3,4}},0]", "-2");
}

/* ---- Arity + shape errors: leave the call unevaluated ---- */

static void test_arity_errors(void) {
    check_out("CharacteristicPolynomial[]", "CharacteristicPolynomial[]");
    check_out("CharacteristicPolynomial[{{1,2},{3,4}}]",
              "CharacteristicPolynomial[{{1, 2}, {3, 4}}]");
    check_out("CharacteristicPolynomial[{{1,2},{3,4}},x,y]",
              "CharacteristicPolynomial[{{1, 2}, {3, 4}}, x, y]");
}

static void test_non_square(void) {
    check_out("CharacteristicPolynomial[{{1,2,3},{4,5,6}},x]",
              "CharacteristicPolynomial[{{1, 2, 3}, {4, 5, 6}}, x]");
}

/* ---- Leak tripwire: many calls, each result freed ---- */

static void test_repeated_no_leak(void) {
    const char* inputs[] = {
        "CharacteristicPolynomial[{{1,2},{3,4}},x]",
        "CharacteristicPolynomial[{{a,b},{c,d}},x]",
        "CharacteristicPolynomial[{{1.1,2.2},{3.3,4.4}},x]",
        "CharacteristicPolynomial[{{{1,2},{5,4}},{{4,3},{6,4}}},x]",
        "CharacteristicPolynomial[]",
    };
    for (int i = 0; i < 25; i++) {
        const char* in = inputs[i % 5];
        Expr* p = parse_expression(in);
        Expr* r = evaluate(p);
        expr_free(p);
        expr_free(r);
    }
    printf("  PASS: 125 repeated calls freed cleanly\n");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running CharacteristicPolynomial tests...\n");
    TEST(test_integer_2x2);
    TEST(test_symbolic_2x2);
    TEST(test_identity_3x3);
    TEST(test_zero_matrix);
    TEST(test_1x1);
    TEST(test_rational_3x3);
    TEST(test_machine_3x3);
    TEST(test_complex_3x3);
    TEST(test_generalized_integer);
    TEST(test_generalized_machine);
    TEST(test_generalized_degree_drop);
    TEST(test_generalized_symbolic);
    TEST(test_packed_and_ndarray);
    TEST(test_numeric_variable);
    TEST(test_arity_errors);
    TEST(test_non_square);
    TEST(test_repeated_no_leak);
    printf("All CharacteristicPolynomial tests passed!\n");

    symtab_clear();
    return 0;
}
