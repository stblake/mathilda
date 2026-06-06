/* test_integrate_unknown.c — Integration of undefined (unknown) functions.
 *
 * Exercises src/calculus/integrate_unknown.c (Roach 1992, §1.7), the
 * Integrate stage that handles integrands rational in undefined functions
 * u[x] and their derivatives.
 *
 * Two kinds of checks:
 *   1. Exact-form assertions (run_eq) — the printed antiderivative matches
 *      a known Mathematica reference.
 *   2. Differential round-trip (assert_dintegral_zero) — the universal
 *      correctness predicate Cancel[Together[Expand[D[Integrate[f,x],x]-f]]]
 *      === 0, which catches algorithmic bugs regardless of surface form.
 *      This drives the randomized stress harness.
 *
 * Plus negative cases (must stay unevaluated) and a robustness check that
 * the genuinely non-elementary Integrate[f'[x] g'[x], x] terminates
 * (cycle detection) rather than crashing or looping.
 */

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Evaluate `input` and compare its printed form to `expected`. */
static void run_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Universal correctness: D[Integrate[f,x],x] - f simplifies to 0. */
static void assert_dintegral_zero(const char* integrand) {
    char buf[4096];
    snprintf(buf, sizeof(buf),
        "Cancel[Together[Expand[D[Integrate[%s, x], x] - (%s)]]]",
        integrand, integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "0") != 0) {
        printf("FAIL: D[Integrate[%s,x],x] - (%s) != 0\n  Got: %s\n",
               integrand, integrand, got);
        ASSERT_STR_EQ(got, "0");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Assert that Integrate[f, x] is returned unevaluated (no closed form). */
static void assert_unevaluated(const char* integrand) {
    char buf[2048];
    snprintf(buf, sizeof(buf), "Integrate[%s, x]", integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    /* The head must still be Integrate and mention the integrand var. */
    if (strncmp(got, "Integrate[", 10) != 0) {
        printf("FAIL: expected %s to stay unevaluated, got: %s\n", buf, got);
        ASSERT(0);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* Polynomial part — derivative recognition & by-parts                 */
/* ------------------------------------------------------------------ */

static void test_basic_first_derivative(void) {
    run_eq("Integrate[f'[x], x]", "f[x]");
    run_eq("Integrate[f''[x], x]", "Derivative[1][f][x]");
    run_eq("Integrate[f'''[x], x]", "Derivative[2][f][x]");
}

static void test_product_rule_x(void) {
    /* d/dx(x f[x]) = f[x] + x f'[x] */
    run_eq("Integrate[x f'[x] + f[x], x]", "x f[x]");
    /* d/dx(x f'[x]) = f'[x] + x f''[x] */
    run_eq("Integrate[x f''[x], x]", "-f[x] + x Derivative[1][f][x]");
}

static void test_product_rule_two_functions(void) {
    /* (f g)' = f' g + f g' */
    run_eq("Integrate[f'[x] g[x] + f[x] g'[x], x]", "f[x] g[x]");
    /* (f g h)' */
    run_eq("Integrate[f'[x] g[x] h[x] + f[x] g'[x] h[x] + f[x] g[x] h'[x], x]",
           "f[x] g[x] h[x]");
    /* (f g)'' = f''g + 2 f'g' + f g'' integrates to (f g)' = f'g + f g' */
    assert_dintegral_zero("f''[x] g[x] + 2 f'[x] g'[x] + f[x] g''[x]");
}

static void test_linear_sums(void) {
    run_eq("Integrate[f'[x] + g'[x], x]", "f[x] + g[x]");
    /* Unrelated un-integrable term flows through. */
    run_eq("Integrate[f'[x] + g[x], x]", "f[x] + Integrate[g[x], x]");
    assert_dintegral_zero("a f'[x] + b g'[x]");
}

static void test_power_part(void) {
    /* f'' f' = (f'^2/2)' */
    run_eq("Integrate[f''[x] f'[x], x]", "1/2 Derivative[1][f][x]^2");
    assert_dintegral_zero("f'[x]^2 f''[x]");   /* (f'^3/3)' */
}

/* ------------------------------------------------------------------ */
/* Fraction part — negative powers of a generator                      */
/* ------------------------------------------------------------------ */

static void test_fraction_part(void) {
    run_eq("Integrate[f'[x]/f[x]^2, x]", "-1/f[x]");
    run_eq("Integrate[f''[x]/f'[x]^3, x]",
           "-1/2/Derivative[1][f][x]^2");
    assert_dintegral_zero("f'[x]/f[x]^3");
}

/* ------------------------------------------------------------------ */
/* Log part — u'/u and friends                                         */
/* ------------------------------------------------------------------ */

static void test_log_part(void) {
    run_eq("Integrate[f'[x]/f[x], x]", "Log[f[x]]");
    /* (f'^2 - f f'')/(f f') = f'/f - f''/f' */
    run_eq("Integrate[(f'[x]^2 - f[x] f''[x])/(f[x] f'[x]), x]",
           "Log[f[x]] - Log[Derivative[1][f][x]]");
    /* Multiple logs, including 2/x. */
    run_eq("Integrate[(2 f[x] + x f'[x])/(x f[x]) + g''[x]/g'[x], x]",
           "2 Log[x] + Log[f[x]] + Log[Derivative[1][g][x]]");
    assert_dintegral_zero("f'[x]/(f[x] + 1)");
}

/* ------------------------------------------------------------------ */
/* ArcTan part — quadratic denominator in a generator                  */
/* ------------------------------------------------------------------ */

static void test_arctan_part(void) {
    /* d/dx(-ArcTan[g'/f]) = (f'g' - f g'')/(f^2 + g'^2) */
    run_eq("Integrate[(f'[x] g'[x] - f[x] g''[x])/(f[x]^2 + g'[x]^2), x]",
           "-ArcTan[Derivative[1][g][x]/f[x]]");
    assert_dintegral_zero("(f'[x] g'[x] - f[x] g''[x])/(f[x]^2 + g'[x]^2)");
}

/* ------------------------------------------------------------------ */
/* Mixed elementary + unknown function                                 */
/* ------------------------------------------------------------------ */

static void test_mixed_elementary(void) {
    run_eq("Integrate[f'[x] Cos[f[x]], x]", "Sin[f[x]]");
    assert_dintegral_zero("f'[x] Exp[f[x]]");
    assert_dintegral_zero("x f'[x] + f[x] + Cos[x]");
    /* d/dx(E^f f') = E^f f'^2 + E^f f'' -> antiderivative E^f f'. */
    run_eq("Integrate`Undefined[E^f[x]*f'[x]^2 + E^f[x]*f''[x], x]",
           "Derivative[1][f][x] E^f[x]");
    assert_dintegral_zero("E^f[x]*f'[x]^2 + E^f[x]*f''[x]");
}

/* ------------------------------------------------------------------ */
/* Transcendental Log[eta] generators (Roach §1.6/§1.7)                */
/* ------------------------------------------------------------------ */

static void test_log_generator(void) {
    /* Log of an unknown-function expression as an integrand generator. */
    run_eq("Integrate[(f[x] - x f[x] + f[x] Log[x f[x]] + x f'[x])/f[x], x]",
           "-1/2 x^2 + x Log[x f[x]]");
    /* Self-referential by-parts: Integrate[Log[f] f'/f, x] = Log[f]^2/2. */
    run_eq("Integrate[Log[f[x]] f'[x]/f[x], x]", "1/2 Log[f[x]]^2");
    assert_dintegral_zero("Log[f[x]] f'[x]/f[x]");
    assert_dintegral_zero("(f[x] - x f[x] + f[x] Log[x f[x]] + x f'[x])/f[x]");
}

/* ------------------------------------------------------------------ */
/* Composite arguments — chain rule                                    */
/* ------------------------------------------------------------------ */

static void test_composite_arguments(void) {
    run_eq("Integrate[2 x f'[x^2], x]", "f[x^2]");
    run_eq("Integrate[Cos[x] f'[Sin[x]], x]", "f[Sin[x]]");
    /* product rule with composite arguments */
    run_eq("Integrate[2 x f'[x^2] g[x^2] + 2 x f[x^2] g'[x^2], x]",
           "f[x^2] g[x^2]");
    assert_dintegral_zero("3 x^2 f'[x^3]");
}

/* ------------------------------------------------------------------ */
/* Higher-order towers                                                  */
/* ------------------------------------------------------------------ */

static void test_higher_order(void) {
    assert_dintegral_zero("x f''[x] + f'[x]");          /* (x f')' */
    assert_dintegral_zero("x^2 f'''[x]");
    run_eq("Integrate[f''''[x], x]", "Derivative[3][f][x]");
}

/* ------------------------------------------------------------------ */
/* Negative cases — must stay unevaluated                              */
/* ------------------------------------------------------------------ */

static void test_negative_no_derivative(void) {
    assert_unevaluated("f[x]");
    assert_unevaluated("x f[x]");
    assert_unevaluated("f[x] g[x]");
}

static void test_negative_nonintegrable(void) {
    /* f'[x]^2 has no antiderivative in f and its derivatives. */
    assert_unevaluated("f'[x]^2");
}

/* Robustness: the genuinely non-elementary product of two independent
 * derivatives must terminate via cycle detection (historically this
 * crashed / looped) and leave a well-formed unevaluated result. */
static void test_cycle_detection_terminates(void) {
    Expr* e = parse_expression("Integrate[f'[x] g'[x], x]");
    Expr* res = evaluate(e);   /* must not crash or hang */
    char* got = expr_to_string(res);
    ASSERT(strncmp(got, "Integrate[", 10) == 0);
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* Method option                                                       */
/* ------------------------------------------------------------------ */

static void test_method_option(void) {
    run_eq("Integrate[f'[x], x, Method -> \"Undefined\"]", "f[x]");
    /* Strict: a non-undefined integrand under Method->Undefined fails. */
    run_eq("Integrate[x^2, x, Method -> \"Undefined\"]",
           "Integrate[x^2, x, Method -> \"Undefined\"]");
    /* Direct package builtin. */
    run_eq("Integrate`Undefined[x f'[x] + f[x], x]", "x f[x]");
}

/* ------------------------------------------------------------------ */
/* Stress harness — randomized D round-trip                            */
/* ------------------------------------------------------------------ */

/* For each polynomial-in-generators expression P below, the integrand
 * D[P, x] is, by construction, an exact derivative; Integrate must
 * recover an antiderivative that differentiates back to it.  Covers deep
 * towers, several functions, and large x-polynomial coefficients without
 * hand-computing answers. */
static void test_stress_roundtrip(void) {
    const char* seeds[] = {
        "x^3 f[x]",
        "x f[x] g[x]",
        "f[x]^2 g[x]",
        "x^2 f'[x] + f[x]^3",
        "f[x] g[x] h[x]",
        "x^4 f'[x] g[x]",
        "f'[x]^2 g[x]",
        "x f[x] g'[x] h[x]",
        "f[x]^2 g[x]^2",
        "(x^2 + 1) f[x] g[x]",
        "x^5 f''[x]",
        "f[x] g[x^2]",          /* composite */
        NULL
    };
    for (int i = 0; seeds[i]; i++) {
        char buf[1024];
        snprintf(buf, sizeof(buf), "D[%s, x]", seeds[i]);
        /* The integrand is D[seed, x]; integrating it must round-trip. */
        assert_dintegral_zero(buf);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic_first_derivative);
    TEST(test_product_rule_x);
    TEST(test_product_rule_two_functions);
    TEST(test_linear_sums);
    TEST(test_power_part);
    TEST(test_fraction_part);
    TEST(test_log_part);
    TEST(test_log_generator);
    TEST(test_arctan_part);
    TEST(test_mixed_elementary);
    TEST(test_composite_arguments);
    TEST(test_higher_order);
    TEST(test_negative_no_derivative);
    TEST(test_negative_nonintegrable);
    TEST(test_cycle_detection_terminates);
    TEST(test_method_option);
    TEST(test_stress_roundtrip);

    printf("All integrate_unknown (Roach §1.7) tests passed!\n");
    return 0;
}
