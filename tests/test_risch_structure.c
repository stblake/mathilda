/* test_risch_structure.c — the Risch structure theorems (Bronstein §9.3).
 *
 * Verifies:
 *   - Risch`RationalSpan: the rational Q-linear membership decision;
 *   - Risch`LogReducible / Risch`ExpReducible: the structure-theorem decision of
 *     whether a proposed new logarithm/exponential is reducible over the tower
 *     (returning the rational coefficients) or a genuinely new monomial (False).
 *
 * Classic reducibility facts under test: log(x^2)=2log(x), log(2x)=log2+log(x),
 * exp(2x)=exp(x)^2, exp(x+log x)=x e^x; and genuine-new cases log(x+1),
 * exp(x^2), exp(x log x)=x^x.
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
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    }
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Assert that the rational expression `diff` is identically zero. */
static void run_zero(const char* diff) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Expand[Together[%s]]", diff);
    run_test(buf, "0");
}

/* ---- The rational Q-span decision ------------------------------------- */
static void test_rational_span(void) {
    run_test("Risch`RationalSpan[1/x, {1/x}, {x}]", "List[1]");
    run_test("Risch`RationalSpan[2/x, {1/x}, {x}]", "List[2]");
    run_test("Risch`RationalSpan[1/(2 x), {1/x}, {x}]", "List[Rational[1, 2]]");
    run_test("Risch`RationalSpan[3/x + 2/(x + 1), {1/x, 1/(x + 1)}, {x}]", "List[3, 2]");
    run_test("Risch`RationalSpan[1/(x + 1), {1/x}, {x}]", "False");
    run_test("Risch`RationalSpan[0, {1/x}, {x}]", "List[0]");
    /* Dependent generators: a solution exists (non-unique); check it reconstructs. */
    run_zero("{1/x, 2/x} . Risch`RationalSpan[1/x, {1/x, 2/x}, {x}] - 1/x");
}

/* ---- Structure theorem: logarithms ------------------------------------ */
static void test_log_reducible(void) {
    /* log(x^2) = 2 log(x): Da/a = 2/x = 2 (1/x). */
    run_test("Risch`LogReducible[x^2, x, {{t1, \"Log\", 1/x}}]", "List[2]");
    /* log(2x) = log 2 + log(x): Da/a = 1/x. */
    run_test("Risch`LogReducible[2 x, x, {{t1, \"Log\", 1/x}}]", "List[1]");
    /* log(x) itself is (trivially) reducible. */
    run_test("Risch`LogReducible[x, x, {{t1, \"Log\", 1/x}}]", "List[1]");
    /* log(x+1) is a genuinely new monomial. */
    run_test("Risch`LogReducible[x + 1, x, {{t1, \"Log\", 1/x}}]", "False");
}

/* ---- Structure theorem: exponentials ---------------------------------- */
static void test_exp_reducible(void) {
    /* exp(2x) = exp(x)^2: Db = 2 = 2 (Dt/t) with Dt/t = 1. */
    run_test("Risch`ExpReducible[2 x, x, {{t1, \"Exp\", t1}}]", "List[2]");
    /* exp(3x): Db = 3. */
    run_test("Risch`ExpReducible[3 x, x, {{t1, \"Exp\", t1}}]", "List[3]");
    /* exp(x) itself. */
    run_test("Risch`ExpReducible[x, x, {{t1, \"Exp\", t1}}]", "List[1]");
    /* exp(x^2) is a genuinely new monomial: Db = 2x is not constant. */
    run_test("Risch`ExpReducible[x^2, x, {{t1, \"Exp\", t1}}]", "False");
}

/* ---- Mixed tower: log(x) and exp(x) together -------------------------- */
static void test_mixed_tower(void) {
    /* tower {t1 = log(x), Dt1 = 1/x ; t2 = exp(x), Dt2 = t2}. */
#define MIX "{{t1, \"Log\", 1/x}, {t2, \"Exp\", t2}}"
    /* exp(x + log x) = x e^x is reducible: Db = 1 + 1/x = 1(1) + 1(1/x). */
    run_test("Risch`ExpReducible[x + t1, x, " MIX "]", "List[1, 1]");
    /* exp(x log x) = x^x is a new monomial: Db = t1 + 1 has a t1 term. */
    run_test("Risch`ExpReducible[x t1, x, " MIX "]", "False");
    /* log(x) reduces onto the log generator only: {1, 0}. */
    run_test("Risch`LogReducible[x, x, " MIX "]", "List[1, 0]");
    /* exp(x) reduces onto the exp generator only: {0, 1}. */
    run_test("Risch`ExpReducible[x, x, " MIX "]", "List[0, 1]");
#undef MIX
}

int main(void) {
    core_init();

    TEST(test_rational_span);
    TEST(test_log_reducible);
    TEST(test_exp_reducible);
    TEST(test_mixed_tower);

    printf("All risch_structure tests passed.\n");
    return 0;
}
