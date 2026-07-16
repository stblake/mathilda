/* test_polynomialsqrt.c — PolynomialSqrt[p] / PolynomialSqrt[p, x].
 *
 * PolynomialSqrt returns s with s^2 == p when p is a perfect square (every
 * non-constant irreducible factor to an even multiplicity; numeric content
 * carried through Sqrt), else $Failed.  This is Cherry-substrate refactor R5 —
 * the polynomial square-root primitive the completing-square Erf-argument step
 * needs (r_i = Sqrt[p + beta q]).  Accepted only after the exact
 * Expand[s^2 - p] == 0 certificate, so a wrong root can never be returned.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Assert PolynomialSqrt gives a genuine root: Expand[s^2 - p] == 0 and s != $Failed. */
static void assert_is_sqrt(const char* p) {
    char buf[512];
    snprintf(buf, sizeof buf, "Expand[PolynomialSqrt[%s]^2 - (%s)]", p, p);
    run_test(buf, "0");
    snprintf(buf, sizeof buf, "PolynomialSqrt[%s] === $Failed", p);
    run_test(buf, "False");
}

static void assert_failed(const char* p) {
    char buf[512];
    snprintf(buf, sizeof buf, "PolynomialSqrt[%s] === $Failed", p);
    run_test(buf, "True");
}

static void test_perfect_squares(void) {
    run_test("PolynomialSqrt[(x + 1)^2]", "Plus[1, x]");
    run_test("PolynomialSqrt[x^2 + 2 x + 1]", "Plus[1, x]");
    run_test("PolynomialSqrt[0]", "0");
    /* Multi-factor, higher even multiplicities. */
    assert_is_sqrt("(x^2 + 1)^2 (x - 3)^4");
    assert_is_sqrt("x^4 + 2 x^2 + 1");
    assert_is_sqrt("(x - 1)^2 (x + 2)^2");
    assert_is_sqrt("(x^2 + x + 1)^2");
    /* Cherry-style completing-square target: q = s^2 with s = x. */
    assert_is_sqrt("x^2");
    /* Numeric-content cases: content carried through Sqrt but still a valid root. */
    assert_is_sqrt("4 (x + 1)^2");
    assert_is_sqrt("2 (x + 1)^2");
    assert_is_sqrt("9 x^2");
}

static void test_content(void) {
    run_test("PolynomialSqrt[4 (x + 1)^2]", "Times[2, Plus[1, x]]");
    run_test("PolynomialSqrt[9 x^2]", "Times[3, x]");
    run_test("PolynomialSqrt[2 (x + 1)^2]", "Times[Power[2, Rational[1, 2]], Plus[1, x]]");
}

static void test_non_squares(void) {
    assert_failed("(x + 1) (x + 2)");           /* distinct linear factors */
    assert_failed("x^2 + 1");                    /* irreducible, multiplicity 1 */
    assert_failed("(x + 1)^2 (x + 2)");          /* one odd multiplicity */
    assert_failed("(x + 1)^3");                  /* odd multiplicity */
    assert_failed("x^3");                         /* odd multiplicity */
    assert_failed("x^2 + x + 1");                 /* irreducible */
}

/* Constants: perfect-square constant -> its root; else Sqrt form (still a root). */
static void test_constants(void) {
    run_test("PolynomialSqrt[4]", "2");
    run_test("Expand[PolynomialSqrt[7]^2 - 7]", "0");
}

int main(void) {
    core_init();

    TEST(test_perfect_squares);
    TEST(test_content);
    TEST(test_non_squares);
    TEST(test_constants);

    printf("All PolynomialSqrt tests passed.\n");
    return 0;
}
