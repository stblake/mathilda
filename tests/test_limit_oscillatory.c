/*
 * test_limit_oscillatory.c -- the oscillatory normal form (src/calculus/
 * limit_osc.c), reached from Limit at +/-Infinity.
 *
 * Four groups, one per decision rule, plus the abstentions. The rules and
 * their proofs are documented at the top of limit_osc.c:
 *
 *   R1  every oscillation decays          -> the limit of the rest
 *   R3  a non-oscillatory part dominates  -> +/-Infinity
 *   R0  one oscillation strictly dominates-> Indeterminate (IVT)
 *   R2  polynomial phases, Weyl means     -> Indeterminate
 *
 * An abstention (the echoed Limit[...] expression) is a *correct* answer
 * whenever a hypothesis is not verifiable, and several tests pin exactly
 * that so a future relaxation has to be deliberate.
 */

#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", input);
    Expr* v = evaluate(e);
    char* got = expr_to_string(v);
    ASSERT_MSG(strcmp(got, expected) == 0,
               "Limit mismatch for %s:\n    expected: %s\n    got:      %s",
               input, expected, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* ----------------------------------------------------------------- */
/* The motivating case.                                              */
/*                                                                    */
/*   (Cos[x^2]/x^2 - Cos[(x+1)^2]/(x+1)^2) x^3                        */
/*     = 2 x Sin[x^2 + x + 1/2] Sin[x + 1/2] + O(1)                   */
/*                                                                    */
/* Unbounded, and with no dominant summand: the four exponential      */
/* groups (phases +/-x^2 and +/-(x^2 + 2x)) all have modulus ~ x/2.   */
/* Only the Weyl mean-square rule settles it.                        */
/* ----------------------------------------------------------------- */
static void test_motivating_case(void) {
    check("Limit[(Cos[x^2]/(x^2) - Cos[(x+1)^2]/((x+1)^2))/(1/(x^3)), "
          "x -> Infinity]", "Indeterminate");
    /* Same shape one power lower: now the amplitudes vanish and the very
     * same normal form yields an honest 0 by the squeeze rule. */
    check("Limit[(Cos[x^2]/(x^2) - Cos[(x+1)^2]/((x+1)^2))/(1/x), "
          "x -> Infinity]", "0");
}

/* ----------------------------------------------------------------- */
/* R1 -- every oscillatory amplitude decays; the limit is that of the */
/* non-oscillatory remainder.                                         */
/* ----------------------------------------------------------------- */
static void test_r1_decaying_oscillation(void) {
    check("Limit[Sin[x]/x, x -> Infinity]", "0");
    check("Limit[Sin[x]/x + 1/x, x -> Infinity]", "0");
    check("Limit[2 + Cos[x^2]/x, x -> Infinity]", "2");
    check("Limit[x Sin[x]/(1 + x^2), x -> Infinity]", "0");
    check("Limit[3 + (Cos[x] + Sin[7 x^3])/Sqrt[x], x -> Infinity]", "3");
    /* The oscillation cancels outright -- no group survives. */
    check("Limit[Sin[x]^2 + Cos[x]^2, x -> Infinity]", "1");
    check("Limit[Sin[x + 2 Pi] - Sin[x], x -> Infinity]", "0");
}

/* ----------------------------------------------------------------- */
/* R3 -- SUM |c_j| stays strictly below |c_0|, which diverges.        */
/* ----------------------------------------------------------------- */
static void test_r3_dominant_smooth_part(void) {
    check("Limit[x + Cos[x], x -> Infinity]", "Infinity");
    check("Limit[x^2 + x Sin[x], x -> Infinity]", "Infinity");
    /* Envelope ratio 1/2 < 1: a term-by-term o(c_0) test would abstain. */
    check("Limit[x^2 (2 + Cos[x]), x -> Infinity]", "Infinity");
    check("Limit[-x^2 (2 + Cos[x]), x -> Infinity]", "-Infinity");
    check("Limit[3 x + 2 Sin[x] + Cos[5 x], x -> Infinity]", "Infinity");
}

/* ----------------------------------------------------------------- */
/* R0 -- one oscillation strictly dominates; the IVT produces two     */
/* distinct accumulation values. No polynomial restriction on the     */
/* phase, so exponential and logarithmic phases are covered.          */
/* ----------------------------------------------------------------- */
static void test_r0_dominant_oscillation(void) {
    check("Limit[x Sin[x], x -> Infinity]", "Indeterminate");
    check("Limit[x Cos[x], x -> -Infinity]", "Indeterminate");
    check("Limit[E^x Cos[x], x -> Infinity]", "Indeterminate");
    check("Limit[x^5 Cos[x], x -> Infinity]", "Indeterminate");
    check("Limit[Sin[Log[x]], x -> Infinity]", "Indeterminate");
    check("Limit[ArcTan[x] Sin[x], x -> Infinity]", "Indeterminate");
}

/* ----------------------------------------------------------------- */
/* R2 -- distinct polynomial phases, Weyl mean / mean-square.         */
/* ----------------------------------------------------------------- */
static void test_r2_polynomial_phases(void) {
    check("Limit[Sin[x], x -> Infinity]", "Indeterminate");
    check("Limit[Cos[x^2], x -> Infinity]", "Indeterminate");
    check("Limit[Sin[x] + Cos[x], x -> Infinity]", "Indeterminate");
    check("Limit[Sin[x^2] + Cos[x], x -> Infinity]", "Indeterminate");
    check("Limit[Cos[x] + Cos[Pi x], x -> Infinity]", "Indeterminate");
    /* Sin[x]^2 = 1/2 - E^(2 I x)/4 - E^(-2 I x)/4: a non-zero constant
     * part plus a surviving oscillation. */
    check("Limit[Sin[x]^2, x -> Infinity]", "Indeterminate");
    /* Phases +/-x +/- x^2, four distinct quadratics. */
    check("Limit[Sin[x] Sin[x^2], x -> Infinity]", "Indeterminate");
    /* A constant phase offset is absorbed into the amplitude, so these
     * two share the phase x and combine to (1 - E^I)/2 != 0. */
    check("Limit[Cos[x] - Cos[x + 1], x -> Infinity]", "Indeterminate");
    /* Real polynomial phase, unbounded amplitude, |c| = O(x^deg). */
    check("Limit[x (Cos[x^2] + Sin[3 x^2]), x -> Infinity]", "Indeterminate");
}

/* ----------------------------------------------------------------- */
/* Honest abstentions.                                                */
/* ----------------------------------------------------------------- */
static void test_abstentions(void) {
    /* A symbolic amplitude: a = 0 gives the limit 0, so no verdict. */
    check("Limit[a Sin[x], x -> Infinity]",
          "Limit[a Sin[x], x -> Infinity]");
    /* Envelope exactly equal to the smooth part -- x^2 (1 + Cos[x]) is
     * >= 0 with zeros at every odd multiple of Pi, so it has no limit,
     * but neither the triangle-inequality bound nor the dominant-
     * oscillation IVT can see that. */
    check("Limit[x^2 (1 + Cos[x]), x -> Infinity]",
          "Limit[x^2 (1 + Cos[x]), x -> Infinity]");
    /* Tan leaves an exponential in a denominator, so the amplitude is
     * not oscillation-free and the normal form is incomplete. */
    check("Limit[Tan[x], x -> Infinity]", "Limit[Tan[x], x -> Infinity]");
}

/* ----------------------------------------------------------------- */
/* A finite limit point reduces to +Infinity through x = a +/- 1/t,   */
/* so an oscillation *at a point* gets the identical normal form. A   */
/* two-sided limit needs both sides to agree.                         */
/* ----------------------------------------------------------------- */
static void test_finite_point(void) {
    check("Limit[Sin[1/x], x -> 0]", "Indeterminate");
    check("Limit[Sin[1/x]/x, x -> 0]", "Indeterminate");
    check("Limit[Sin[1/x] + Cos[1/x], x -> 0]", "Indeterminate");
    check("Limit[Sin[1/x]^2, x -> 0]", "Indeterminate");
    /* Plain substitution used to fold both terms to Cos[ComplexInfinity]
     * and cancel them to a confident 0. The function oscillates over
     * [-2 Sin[1/2], 2 Sin[1/2]]. */
    check("Limit[Cos[1/x] - Cos[1/x + 1], x -> 0]", "Indeterminate");
    check("Limit[Sin[1/x]/x, x -> 0, Direction -> \"FromAbove\"]",
          "Indeterminate");
    /* The decaying envelope still squeezes to 0 on this path. */
    check("Limit[x Sin[1/x], x -> 0]", "0");
    check("Limit[x^2 Sin[1/x^3], x -> 0]", "0");
    /* Away from the singularity nothing changes. */
    check("Limit[Sin[1/x], x -> 1]", "Sin[1]");
}

/* ----------------------------------------------------------------- */
/* A pure imaginary exponent. The f^g -> Exp[g Log f] reduction used  */
/* to fold Limit[I x] = DirectedInfinity[I] straight back through Exp */
/* and hand back the unfolded `E^DirectedInfinity[I]` -- swallowing   */
/* the 1/x in the second case, which is a wrong answer, not just an   */
/* ugly one. Both now reach the oscillatory layer.                    */
/* ----------------------------------------------------------------- */
static void test_imaginary_exponent(void) {
    check("Limit[E^(I x), x -> Infinity]", "Indeterminate");
    check("Limit[E^(I x)/x, x -> Infinity]", "0");
    check("Limit[E^(I x) + E^(-I x), x -> Infinity]", "Indeterminate");
    /* Unmated and unbounded: |f| -> Infinity with a rotating direction, so
     * ComplexInfinity is arguably the answer and the layer abstains rather
     * than pick one. */
    check("Limit[x E^(I x), x -> Infinity]",
          "Limit[x E^(I x), x -> Infinity]");
}

/* ----------------------------------------------------------------- */
/* Method -> "Oscillatory" selects this layer alone.                  */
/* ----------------------------------------------------------------- */
static void test_method_selection(void) {
    check("Limit[Sin[x] + Cos[x], x -> Infinity, Method -> \"Oscillatory\"]",
          "Indeterminate");
    check("Limit[x^2 (2 + Cos[x]), x -> Infinity, Method -> \"Oscillatory\"]",
          "Infinity");
    /* No oscillation at all -> the layer abstains and, with the cascade
     * restricted to it, the whole Limit stays unevaluated. */
    check("Limit[x^2, x -> Infinity, Method -> \"Oscillatory\"]",
          "Limit[x^2, x -> Infinity, Method -> \"Oscillatory\"]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_motivating_case);
    TEST(test_r1_decaying_oscillation);
    TEST(test_r3_dominant_smooth_part);
    TEST(test_r0_dominant_oscillation);
    TEST(test_r2_polynomial_phases);
    TEST(test_abstentions);
    TEST(test_finite_point);
    TEST(test_imaginary_exponent);
    TEST(test_method_selection);

    printf("All oscillatory-limit tests passed.\n");
    return 0;
}
