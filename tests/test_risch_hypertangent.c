/* test_risch_hypertangent.c — the hypertangent case (Bronstein §5.10).
 *
 * Verifies:
 *   - Risch`PolynomialReduce (§5.4): reducing a polynomial modulo derivatives in
 *     a nonlinear monomial to a remainder of degree < deg_t(Dt);
 *   - Risch`IntegrateHypertangentPolynomial (§5.10): the polynomial-part tangent
 *     integrator, returning {q, c} with p - D[q] - c D(t^2+1)/(t^2+1) in k.
 *
 * The monomial t = tan(x) has Dt = 1 + t^2 (a = 1); t = tan(2x) has Dt = 2(1+t^2)
 * (a = 2).  Worked example 5.10.1: integral of (tan^2 x + x tan x + 1).
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

static void run_zero(const char* diff) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Expand[Together[%s]]", diff);
    run_test(buf, "0");
}

/* t = tan(x): Dt = 1 + t^2. */
#define HT "{x -> 1, t -> 1 + t^2}"

/* ---- PolynomialReduce (§5.4) ------------------------------------------ */
static void test_polynomial_reduce(void) {
    /* Example 5.10.1: reduce t^2 + x t + 1 -> q = t, r = x t. */
    run_zero("Risch`PolynomialReduce[t^2 + x t + 1, t, " HT "][[1]] - t");
    run_zero("Risch`PolynomialReduce[t^2 + x t + 1, t, " HT "][[2]] - x t");
    /* Reconstruction p = D[q] + r for a higher-degree p. */
    run_zero("Risch`Derivation[Risch`PolynomialReduce[t^4 + t + x, t, " HT "][[1]], " HT "] "
             "+ Risch`PolynomialReduce[t^4 + t + x, t, " HT "][[2]] - (t^4 + t + x)");
    /* The remainder has degree < deg_t(Dt) = 2. */
    run_zero("Coefficient[Risch`PolynomialReduce[t^4 + t + x, t, " HT "][[2]], t, 2]");
    run_zero("Coefficient[Risch`PolynomialReduce[t^4 + t + x, t, " HT "][[2]], t, 3]");
}

/* ---- IntegrateHypertangentPolynomial (§5.10) -------------------------- */
static void test_integrate_hypertangent_poly(void) {
#define IHP "Risch`IntegrateHypertangentPolynomial"
    /* Example 5.10.1: {q, c} = {t, x/2}. */
    run_zero(IHP "[t^2 + x t + 1, t, " HT "][[1]] - t");
    run_zero(IHP "[t^2 + x t + 1, t, " HT "][[2]] - x/2");
    /* D[c] = 1/2 != 0, so (the x tan x part of) the integral is not elementary. */
    run_test("D[" IHP "[t^2 + x t + 1, t, " HT "][[2]], x]", "Rational[1, 2]");
    /* The §5.10 certificate: p - D[q] - c D[t^2+1]/(t^2+1) is in k (here 0). */
    run_zero("(t^2 + x t + 1) - Risch`Derivation[" IHP "[t^2 + x t + 1, t, " HT "][[1]], " HT "] "
             "- " IHP "[t^2 + x t + 1, t, " HT "][[2]] Risch`Derivation[t^2 + 1, " HT "]/(t^2 + 1)");

    /* Integral of tan(x): {0, 1/2}; D[c] = 0, so it IS elementary
     * (∫ tan x = (1/2) Log[tan^2 x + 1] = -Log Cos x). */
    run_test(IHP "[t, t, " HT "]", "List[0, Rational[1, 2]]");
    run_test("D[" IHP "[t, t, " HT "][[2]], x]", "0");

    /* tan(2x): Dt = 2(1 + t^2), a = 2, so c = 1/(2 a) = 1/4. */
    run_test(IHP "[t, t, {x -> 1, t -> 2 (1 + t^2)}]", "List[0, Rational[1, 4]]");
#undef IHP
}

int main(void) {
    core_init();

    TEST(test_polynomial_reduce);
    TEST(test_integrate_hypertangent_poly);

    printf("All risch_hypertangent tests passed.\n");
    return 0;
}
