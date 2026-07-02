/* test_extension_options.c — Phase 0 of the Integrate plan.
 *
 * Coverage: the `Extension -> α` option on PolynomialGCD, PolynomialLCM,
 * Cancel, Together, and Apart.  These tests pin the behavior the
 * BronsteinRational pipeline relies on: GCD / cancellation closed over
 * Q(α), with Extension -> None / no option preserving the existing
 * (pre-Phase-0) behavior bit-for-bit so unrelated test suites don't
 * regress.
 *
 * Out of scope for Phase 0:
 *   - Extension -> Automatic (auto-detection of algebraic numbers in
 *     input — currently treated as Extension -> None).
 *   - Extension -> {α₁, α₂, ...} tower form (currently falls back to
 *     no-extension on the GCD/LCM/Cancel/Together path; Factor already
 *     handles this case end-to-end via qa_factor_with_extension_tower).
 *   - Multivariate polynomials with an extension (the QAUPoly substrate
 *     is single-variable; multivariate-with-extension falls back to the
 *     standard path).
 *   - Together-with-Extension on inputs whose summands have Sqrt[α]-
 *     laden denominators: the underlying PolynomialQuotient does not
 *     yet accept Extension, so the no-α multivariate-GCD path stalls.
 *     Phase 0.5 will add Extension to PolynomialQuotient / Remainder.
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

/* Compare evaluator output (default form) against an `expected` string. */
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

/* ------------------------------------------------------------------ */
/* Option parsing                                                      */
/* ------------------------------------------------------------------ */

void test_parsing(void) {
    /* Extension -> None is consumed but treated as no option. */
    run_eq("Cancel[(x^2 - 4)/(x - 2), Extension -> None]", "2 + x");
    run_eq("Together[1/x + 1/x, Extension -> None]", "2/x");
    run_eq("PolynomialGCD[x^2 - 1, x - 1, Extension -> None]", "-1 + x");
    run_eq("PolynomialLCM[x - 1, x + 1, Extension -> None]", "(-1 + x) (1 + x)");

    /* Extension -> Automatic is currently identical to None (auto-
     * detection deferred). */
    run_eq("Cancel[(x^2 - 4)/(x - 2), Extension -> Automatic]", "2 + x");

    /* Unrecognized options leave the call unevaluated (mirroring
     * Mathematica's behavior of not silently swallowing typos). */
    run_eq("Cancel[expr, Foo -> Bar]", "Cancel[expr, Foo -> Bar]");
}

/* ------------------------------------------------------------------ */
/* PolynomialGCD with Extension -> α                                   */
/* ------------------------------------------------------------------ */

void test_polygcd_extension(void) {
    /* Sqrt[2]: x^2 - 2 = (x - Sqrt[2])(x + Sqrt[2]) over Q[Sqrt[2]]. */
    run_eq("PolynomialGCD[x^2 - 2, x - Sqrt[2], Extension -> Sqrt[2]]",
           "-Sqrt[2] + x");
    run_eq("PolynomialGCD[x^2 - 2, x + Sqrt[2], Extension -> Sqrt[2]]",
           "Sqrt[2] + x");

    /* Sqrt[3] in higher-degree polynomial. */
    run_eq("PolynomialGCD[x^4 - 9, x^2 - 3, Extension -> Sqrt[3]]",
           "-3 + x^2");

    /* Cube root: x^3 - 2 has only one rational root in Q[2^(1/3)],
     * namely x = 2^(1/3).  GCD with x - 2^(1/3) is x - 2^(1/3). */
    run_eq("PolynomialGCD[x^3 - 2, x - 2^(1/3), Extension -> 2^(1/3)]",
           "-2^(1/3) + x");

    /* I (imaginary unit): x^2 + 1 = (x - I)(x + I) over Q[I]. */
    run_eq("PolynomialGCD[x^2 + 1, x - I, Extension -> I]", "-I + x");

    /* Three-arg form folds left-to-right. */
    run_eq("PolynomialGCD[x^4 - 4, x^2 - 2, x - Sqrt[2], "
           "Extension -> Sqrt[2]]", "-Sqrt[2] + x");
}

/* ------------------------------------------------------------------ */
/* PolynomialLCM with Extension -> α                                   */
/* ------------------------------------------------------------------ */

void test_polylcm_extension(void) {
    /* lcm(x - Sqrt[2], x + Sqrt[2]) = x^2 - 2 over Q[Sqrt[2]]. The
     * extension path returns the monic, expanded form (qaupoly_to_expr
     * renders the polynomial in canonical Plus form). */
    run_eq("PolynomialLCM[x - Sqrt[2], x + Sqrt[2], "
           "Extension -> Sqrt[2]]",
           "-2 + x^2");

    /* When one input divides the other, the lcm is the larger one. */
    run_eq("PolynomialLCM[x^2 - 2, x - Sqrt[2], Extension -> Sqrt[2]]",
           "-2 + x^2");
}

/* ------------------------------------------------------------------ */
/* Cancel with Extension -> α                                          */
/* ------------------------------------------------------------------ */

void test_cancel_extension(void) {
    /* (x^2 - 2)/(x - Sqrt[2]) = x + Sqrt[2] over Q[Sqrt[2]]. */
    run_eq("Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]",
           "Sqrt[2] + x");

    /* (x^4 - 4)/(x^2 - 2) = x^2 + 2 (works without extension too, but
     * verifies the extension path doesn't break easy cases). */
    run_eq("Cancel[(x^4 - 4)/(x^2 - 2), Extension -> Sqrt[2]]",
           "2 + x^2");

    /* Without an explicit extension option: the FLINT engine auto-detects the
     * field Q(Sqrt[2]) from the operands and cancels rigorously,
     * (x^2 - 2)/(x - Sqrt[2]) = x + Sqrt[2]. Without FLINT, Sqrt[2] is opaque to
     * the classical path and the fraction is returned essentially unchanged. */
#ifdef USE_FLINT
    run_eq("Cancel[(x^2 - 2)/(x - Sqrt[2])]", "Sqrt[2] + x");
#else
    run_eq("Cancel[(x^2 - 2)/(x - Sqrt[2])]", "(-2 + x^2)/(-Sqrt[2] + x)");
#endif

    /* Cube-root extension: (x^3 - 2)/(x - 2^(1/3)) = x^2 + 2^(1/3) x + 2^(2/3). */
    run_eq("Cancel[(x^3 - 2)/(x - 2^(1/3)), Extension -> 2^(1/3)]",
           "2^(2/3) + 2^(1/3) x + x^2");

    /* Imaginary unit. */
    run_eq("Cancel[(x^2 + 1)/(x - I), Extension -> I]", "I + x");
}

/* ------------------------------------------------------------------ */
/* Together with Extension -> α                                        */
/* ------------------------------------------------------------------ */

void test_together_extension(void) {
    /* 1/(x-Sqrt[2]) + 1/(x+Sqrt[2]) = 2x/(x^2 - 2). */
    run_eq("Together[1/(x - Sqrt[2]) + 1/(x + Sqrt[2]), "
           "Extension -> Sqrt[2]]",
           "(2 x)/(-2 + x^2)");

    /* 1/(x - I) + 1/(x + I) = 2x/(x^2 + 1) over Q[I].  Q[I] = Q(ζ_4) is
     * the cyclotomic field for the primitive 4th root of unity, so the
     * extension-aware Together now routes through the cyclotomic lift
     * (rootofunity.c): the Complex[0,1] coefficients lift to the ζ_4
     * basis and the denominator reduces to the canonical x^2 + 1,
     * matching the Sqrt[2] case above and Mathematica. */
    run_eq("Together[1/(x - I) + 1/(x + I), Extension -> I]",
           "(2 x)/(1 + x^2)");

    /* Cyclotomic Q(ζ_6): a = (-1)^(1/3).  1/(x-a) + 1/(x+a) = 2x/(x^2-a^2),
     * and a^2 = (-1)^(2/3) reduces modulo Φ_6 (ζ^2 = ζ - 1) to
     * (-1)^(1/3) - 1, so the denominator collapses to x^2 + 1 - (-1)^(1/3).
     * This is the reduce-mod-minimal-polynomial behavior that makes
     * cyclotomic Together fast instead of blowing up the Q-path. */
    run_eq("Together[1/(x - (-1)^(1/3)) + 1/(x + (-1)^(1/3)), "
           "Extension -> (-1)^(1/3)]",
           "(2 x)/(1 - (-1)^(1/3) + x^2)");

    /* Same field auto-detected (base -1 collected by extension_autodetect)
     * and cancelled: (x^2 - (-1)^(2/3))/(x - (-1)^(1/3)) = x + (-1)^(1/3). */
    run_eq("Cancel[(x^2 - (-1)^(2/3))/(x - (-1)^(1/3)), Extension -> Automatic]",
           "(-1)^(1/3) + x");

    /* Auto-detected cyclotomic Together with a = ζ_3 = (-1)^(2/3): a^2 =
     * (-1)^(4/3) reduces modulo Φ_6 to -(-1)^(1/3), so x^2 - a^2 collapses
     * to x^2 + (-1)^(1/3). */
    run_eq("Together[1/(x - (-1)^(2/3)) + 1/(x + (-1)^(2/3)), "
           "Extension -> Automatic]",
           "(2 x)/((-1)^(1/3) + x^2)");

    /* Higher order Q(ζ_10), a = (-1)^(1/5): φ(10) = 4, so (-1)^(2/5) is a
     * basis element and stays — the denominator is x^2 - (-1)^(2/5). */
    run_eq("Together[1/(x - (-1)^(1/5)) + 1/(x + (-1)^(1/5)), "
           "Extension -> Automatic]",
           "(2 x)/(-(-1)^(2/5) + x^2)");

    /* Cancel over Q[I] with explicit Extension -> I. */
    run_eq("Cancel[(x^2 + 1)/(x - I), Extension -> I]", "I + x");

    /* No-extension default: the existing structural combine still works
     * for Sqrt[3] inputs because the LCM happens to factor cleanly. */
    run_eq("Together[1/(x - Sqrt[3]) + 1/(x + Sqrt[3])]",
           "(2 x)/(-3 + x^2)");
}

/* ------------------------------------------------------------------ */
/* Apart with Extension -> α                                            */
/* ------------------------------------------------------------------ */

void test_apart_extension(void) {
    /* Apart[1/(x^2 - 2), x, Extension -> Sqrt[2]] splits over Q[Sqrt[2]].
     * Mathilda's canonical output for this is the partial-fraction
     * decomposition with Sqrt[2] coefficients.  Pin the exact form so
     * we catch any future regressions in the algebraic-extension Apart
     * path. With the radical canonicalization in builtin_times, the
     * coefficient 1/2 * 1/Sqrt[2] (and Sqrt[2]/4) collapse to the unified
     * 1/2^(3/2) form. */
    run_eq("Apart[1/(x^2 - 2), x, Extension -> Sqrt[2]]",
           "-1/2 1/(Sqrt[2] (Sqrt[2] + x)) + 1/2 1/(Sqrt[2] (-Sqrt[2] + x))");

    /* No-extension default: Apart[1/(x^2 - 2), x] returns the input
     * unchanged because the denominator is irreducible over Q. */
    run_eq("Apart[1/(x^2 - 2), x]", "1/(-2 + x^2)");

    /* Apart over Q[I] splits 1/(x^2 + 1).  Mathilda auto-evaluates 1/I
     * to -I, producing the I*… coefficient form rather than a 1/(I*…)
     * fraction. */
    run_eq("Apart[1/(x^2 + 1), x, Extension -> I]",
           "(-1/2*I)/(-I + x) + (1/2*I)/(I + x)");

    /* Composite: numerator/denominator share an algebraic factor; Apart
     * cancels it via Together-with-Extension and returns the cleaned
     * partial fraction. */
    run_eq("Apart[(x^4 - 4)/((x - 1) (x^2 - 2)), x, "
           "Extension -> Sqrt[2]]", "1 + x + 3/(-1 + x)");
}

/* ------------------------------------------------------------------ */
/* PolynomialQuotient / PolynomialRemainder with Extension -> α        */
/* ------------------------------------------------------------------ */

void test_polydivrem_extension(void) {
    /* Sqrt[2]: x^2 - 2 = (x - Sqrt[2])(x + Sqrt[2]) over Q[Sqrt[2]]. */
    run_eq("PolynomialQuotient[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]",
           "Sqrt[2] + x");
    run_eq("PolynomialRemainder[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]",
           "0");

    /* Cube root: x^3 - 2 = (x - 2^(1/3))(x^2 + 2^(1/3) x + 2^(2/3)). */
    run_eq("PolynomialQuotient[x^3 - 2, x - 2^(1/3), x, Extension -> 2^(1/3)]",
           "2^(2/3) + 2^(1/3) x + x^2");
    run_eq("PolynomialRemainder[x^3 - 2, x - 2^(1/3), x, Extension -> 2^(1/3)]",
           "0");

    /* Imaginary unit. */
    run_eq("PolynomialQuotient[x^2 + 1, x - I, x, Extension -> I]", "I + x");
    run_eq("PolynomialRemainder[x^2 + 1, x - I, x, Extension -> I]", "0");

    /* Higher-degree case. */
    run_eq("PolynomialQuotient[x^4 - 4, x^2 - 2, x, Extension -> Sqrt[2]]",
           "2 + x^2");
    run_eq("PolynomialRemainder[x^4 - 4, x^2 - 2, x, Extension -> Sqrt[2]]",
           "0");

    /* Non-trivial remainder. */
    run_eq("PolynomialRemainder[x^2 - 3, x - Sqrt[2], x, Extension -> Sqrt[2]]",
           "-1");
    run_eq("PolynomialQuotient[x^2 - 3, x - Sqrt[2], x, Extension -> Sqrt[2]]",
           "Sqrt[2] + x");

    /* Extension -> None / Automatic / no option all behave the same as
     * the standard path: Sqrt[2] is treated as opaque, but division of
     * x^2 - 2 by x - Sqrt[2] still succeeds because the quotient has
     * Sqrt[2] in its constant term (no rational arithmetic with Sqrt[2]
     * needed to carry out the steps). */
    run_eq("PolynomialQuotient[x^2 - 2, x - Sqrt[2], x]", "Sqrt[2] + x");
    run_eq("PolynomialQuotient[x^2 - 2, x - Sqrt[2], x, Extension -> None]",
           "Sqrt[2] + x");
    run_eq("PolynomialQuotient[x^2 - 2, x - Sqrt[2], x, Extension -> Automatic]",
           "Sqrt[2] + x");

    /* Multivariate fallback: extension path requires univariate input
     * (single live variable other than the alpha render symbol).  When
     * a second variable is present the lift fails and we fall through
     * to the standard path. */
    run_eq("PolynomialQuotient[x^2 + y, x - 1, x, Extension -> Sqrt[2]]",
           "1 + x");
}

/* ------------------------------------------------------------------ */
/* Backward compatibility: no-option calls match pre-Phase-0 output    */
/* ------------------------------------------------------------------ */

void test_backward_compat(void) {
    /* These should match the existing rat_tests / parfrac_tests
     * expectations exactly — no behavioral drift on the default path. */
    run_eq("Cancel[(x^2 - 4)/(x - 2)]", "2 + x");
    run_eq("Together[1/x + 1/(x + 1)]", "(1 + 2 x)/(x + x^2)");
    run_eq("PolynomialGCD[x^2 - 1, x - 1]", "-1 + x");
    run_eq("PolynomialLCM[x - 1, x + 1]", "(-1 + x) (1 + x)");
    run_eq("Apart[1/((x - 1) (x - 2)), x]", "1/(-2 + x) - 1/(-1 + x)");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_parsing);
    TEST(test_polygcd_extension);
    TEST(test_polylcm_extension);
    TEST(test_cancel_extension);
    TEST(test_together_extension);
    TEST(test_apart_extension);
    TEST(test_polydivrem_extension);
    TEST(test_backward_compat);

    printf("All Extension-option tests passed!\n");
    return 0;
}
