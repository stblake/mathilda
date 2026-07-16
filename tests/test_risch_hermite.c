/* test_risch_hermite.c — the Hermite reduction over a monomial extension
 * (Bronstein §5.3, Theorem 5.3.1 / the quadratic HermiteReduce, p.139).
 *
 * Verifies:
 *   - Risch`HermiteReduce returns {g, h, r} with the exact invariant
 *     f == D[g] + h + r (the strongest correctness certificate);
 *   - Bronstein Example 5.3.1 exact values;
 *   - h is simple (squarefree denominator), r is reduced (no normal poles);
 *   - correctness with arbitrary rational k = C(x) coefficients — the case the
 *     earlier polynomial-ansatz Hermite could not reach;
 *   - across log / exponential / hypertangent monomials and higher pole orders.
 *
 * Derivations:  t = tan(x): Dt = 1 + t^2;  t = log(x): Dt = 1/x;  t = e^x: Dt = t.
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

static void run_zero(const char* diff) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Expand[Together[%s]]", diff);
    run_test(buf, "0");
}

#define TAN "{x -> 1, t -> 1 + t^2}"
#define LOG "{x -> 1, t -> 1/x}"
#define EXP "{x -> 1, t -> t}"
#define HR  "Risch`HermiteReduce"

/* The defining invariant: f == D[g] + h + r. */
static void run_invariant(const char* f, const char* deriv) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "With[{r = " HR "[%s, t, %s]}, "
        "Expand[Together[(%s) - (Risch`Derivation[r[[1]], %s] + r[[2]] + r[[3]])]]]",
        f, deriv, f, deriv);
    run_test(buf, "0");
}

/* h (the simple part = r[[2]]) has a SQUAREFREE denominator in t: gcd with its
 * t-derivative is a unit (degree 0 in t). */
static void run_h_simple(const char* f, const char* deriv) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "With[{h = " HR "[%s, t, %s][[2]]}, "
        "With[{d = Denominator[Together[h]]}, "
        "FreeQ[PolynomialGCD[d, D[d, t]], t]]]",
        f, deriv);
    run_test(buf, "True");
}

/* ---- Bronstein Example 5.3.1 (exact values) --------------------------- */
static void test_example_531(void) {
    /* f = (x - t)/t^2, t = tan(x): HermiteReduce -> (-x/t, 0, -x). */
    run_invariant("(x - t)/t^2", TAN);
    /* The remaining reduced integrand is -x (so the integral is -x/tan(x) - x^2/2). */
    run_zero(HR "[(x - t)/t^2, t, " TAN "][[3]] + x");
    run_zero(HR "[(x - t)/t^2, t, " TAN "][[1]] + x/t");
    run_test(HR "[(x - t)/t^2, t, " TAN "][[2]]", "0");
}

/* ---- Invariant across monomial kinds ---------------------------------- */
static void test_invariant_log(void) {
    run_invariant("1/t^2", LOG);                 /* 1/log(x)^2 -> -x/log + li part */
    run_invariant("1/t^3", LOG);
    run_invariant("(x + t)/t^2", LOG);
    run_invariant("1/(t^2 (t + 1))", LOG);       /* normal repeated + normal simple */
    run_invariant("x/((x + 1) t^2)", LOG);       /* rational-in-x coefficient */
}

static void test_invariant_exp(void) {
    run_invariant("1/(t - 1)^2", EXP);           /* t = e^x; t-1 normal */
    run_invariant("1/(t - 1)^3", EXP);
    run_invariant("(x^2 t + 1)/(t + x)^2", EXP);
    run_invariant("1/((t - 1)^2 (t - 2))", EXP);
    run_invariant("t/((t - x)^2)", EXP);         /* x-dependent normal pole */
}

static void test_invariant_tan(void) {
    /* For tan, t^2+1 is SPECIAL: it must NOT be Hermite-reduced (stays in r). */
    run_invariant("(x - t)/t^2", TAN);
    run_invariant("1/(t - x)^2", TAN);           /* t-x is normal for tan */
    run_invariant("1/((t^2 + 1)^2)", TAN);       /* special denom -> entirely in r */
    run_invariant("(t + 1)/((t - 1)^2 (t^2 + 1))", TAN);  /* mixed normal + special */
}

/* ---- Rational k = C(x) coefficients (the ansatz-Hermite gap) ---------- */
static void test_rational_coefficients(void) {
    /* Numerator/denominator coefficients are genuine rational functions of x;
     * the earlier polynomial-degree-bounded ansatz would decline these. */
    run_invariant("1/((x^2 + 1) t^2)", LOG);
    run_invariant("x/((x^2 - 2) (t - 1)^2)", EXP);
    run_invariant("(x/(x + 1) + t)/(t^2 + x)^2", EXP);
    run_invariant("1/((x^3 + x + 1) t^3)", LOG);
}

/* ---- h simple, r reduced --------------------------------------------- */
static void test_h_simple_r_reduced(void) {
    run_h_simple("1/t^3", LOG);
    run_h_simple("(x + t)/t^2", LOG);
    run_h_simple("1/(t - 1)^3", EXP);
    run_h_simple("1/((t - 1)^2 (t - 2)^2)", EXP);
    /* r (the reduced part) for tan carries the special (t^2+1) poles unchanged. */
    run_test("FreeQ[" HR "[1/(t^2 + 1)^2, t, " TAN "][[1]], t]", "True");  /* g = 0 */
    run_test(HR "[1/(t^2 + 1)^2, t, " TAN "][[2]]", "0");                  /* h = 0 */
}

/* ---- Already-simple / polynomial inputs: g = 0 ------------------------ */
static void test_no_normal_pole(void) {
    /* A polynomial in t: g = 0, h = 0, r = the polynomial. */
    run_test(HR "[t^2 + x t + 1, t, " TAN "][[1]]", "0");
    run_test(HR "[t^2 + x t + 1, t, " TAN "][[2]]", "0");
    run_zero(HR "[t^2 + x t + 1, t, " TAN "][[3]] - (t^2 + x t + 1)");
    /* A simple (squarefree-denominator) fraction: g = 0, it is already h. */
    run_test(HR "[1/(t - x), t, " EXP "][[1]]", "0");
    run_invariant("1/(t - x)", EXP);
    run_invariant("(3 t + x)/(t^2 - x)", EXP);   /* squarefree denom */
    /* Zero. */
    run_test(HR "[0, t, " TAN "]", "List[0, 0, 0]");
}

/* ---- Higher multiplicities + multiple factors ------------------------- */
static void test_higher_multiplicity(void) {
    run_invariant("1/t^4", LOG);
    run_invariant("1/t^5", LOG);
    run_invariant("1/((t - 1)^3 (t - 2)^2)", EXP);
    run_invariant("(t^2 + 1)/((t - 3)^3)", EXP);   /* constant pole location */
    run_invariant("1/((t - 1)^2 (t - 2)^3 (t - 3))", EXP);
}

/* ---- Robustness ------------------------------------------------------- */
static void test_robustness(void) {
    /* t not a monomial variable of the derivation -> unevaluated. */
    run_test("Head[" HR "[1/s^2, s, " TAN "]] === " HR, "True");
    /* Wrong arity / malformed derivation -> unevaluated. */
    run_test("Head[" HR "[1/t^2, t]] === " HR, "True");
    run_test("Head[" HR "[1/t^2, t, {x, t}]] === " HR, "True");
}

int main(void) {
    core_init();

    TEST(test_example_531);
    TEST(test_invariant_log);
    TEST(test_invariant_exp);
    TEST(test_invariant_tan);
    TEST(test_rational_coefficients);
    TEST(test_h_simple_r_reduced);
    TEST(test_no_normal_pole);
    TEST(test_higher_multiplicity);
    TEST(test_robustness);

    printf("All risch_hermite tests passed.\n");
    return 0;
}
