/* test_integrate_symmetry.c
 *
 * Tests for definite integration by origin-symmetry reduction
 * (src/calculus/integrate_symmetry.c):
 *   Integrate[f, {x, -c, c}]                         (auto-dispatch, after residue)
 *   Integrate[f, {x, -c, c}, Method -> "Symmetry"]   (strict: no fallback)
 *   Integrate`Symmetry[f, {x, -c, c}]                (explicit entry point)
 *
 * Coverage: odd -> 0 (finite, symbolic-bound, and convergent infinite), even
 * -> 2*half (finite and infinite), the dispatch/regression guarantees that the
 * residue families keep their clean closed forms, and the negative controls
 * that MUST stay unevaluated (asymmetric interval, no definite parity, and --
 * critically -- an odd integrand whose half diverges, whose principal value is
 * 0 but whose integral does not converge).
 *
 * Closed forms that Mathilda leaves in an equivalent-but-unsimplified surface
 * form are pinned numerically via Chop[N[value - reference]] == 0.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* -------------------------------------------------------------------------
 * Unit: the explicit Integrate`Symmetry entry point in isolation.
 * ---------------------------------------------------------------------- */
static void test_unit_odd(void) {
    /* Odd polynomials over a finite symmetric interval -> 0. */
    check_eq("Integrate`Symmetry[x^3, {x, -1, 1}]", "0");
    check_eq("Integrate`Symmetry[x^5 - 3 x^3 + x, {x, -2, 2}]", "0");
    /* Odd rational, symbolic bound -> 0 (half is finite in the bound). */
    check_eq("Integrate`Symmetry[x^3/(1+x^2), {x, -a, a}]", "0");
    check_eq("Integrate`Symmetry[x/(x^2 + 4), {x, -b, b}]", "0");
    /* Odd, convergent infinite interval -> 0. */
    check_eq("Integrate`Symmetry[x Exp[-x^2], {x, -Infinity, Infinity}]", "0");
    /* Odd trig over a full period -> 0 (half Integrate[Sin,{x,0,Pi}] = 2). */
    check_eq("Integrate`Symmetry[Sin[x], {x, -Pi, Pi}]", "0");
    /* Odd = even * odd: x Cos[x] (half computes by parts) -> 0. */
    check_eq("Integrate`Symmetry[x Cos[x], {x, -Pi, Pi}]", "0");
}

static void test_unit_even(void) {
    /* Even polynomial -> 2*half. */
    check_eq("Integrate`Symmetry[x^2, {x, -1, 1}]", "2/3");
    check_eq("Integrate`Symmetry[x^2, {x, -3, 3}]", "18");
    /* Even rational, finite interval -> 2 ArcTan[1] = Pi/2. */
    check_eq("Integrate`Symmetry[1/(1+x^2), {x, -1, 1}]", "1/2 Pi");
    /* Even rational, whole line: 2*half unlocks the half-line closed form. */
    check_eq("Integrate`Symmetry[1/(1+x^4), {x, -Infinity, Infinity}]", "Pi/Sqrt[2]");
    /* Cos is even: 2*Integrate[Cos, {x,0,Pi}] = 2 Sin[Pi] = 0. */
    check_eq("Integrate`Symmetry[Cos[x], {x, -Pi, Pi}]", "0");
}

/* -------------------------------------------------------------------------
 * Unit negative controls: MUST stay unevaluated (return NULL from the builtin).
 * ---------------------------------------------------------------------- */
static void test_unit_negative(void) {
    /* Asymmetric interval. */
    check_eq("Integrate`Symmetry[x^2, {x, 0, 1}]",
             "Integrate`Symmetry[x^2, {x, 0, 1}]");
    /* Neither odd nor even. */
    check_eq("Integrate`Symmetry[x^2 + x, {x, -1, 1}]",
             "Integrate`Symmetry[x + x^2, {x, -1, 1}]");
    /* Odd but the half diverges: PV is 0, the integral is not -- do NOT say 0. */
    check_eq("Integrate`Symmetry[1/x, {x, -1, 1}]",
             "Integrate`Symmetry[1/x, {x, -1, 1}]");
    check_eq("Integrate`Symmetry[x/(1+x^2), {x, -Infinity, Infinity}]",
             "Integrate`Symmetry[x/(1 + x^2), {x, -Infinity, Infinity}]");
}

/* -------------------------------------------------------------------------
 * Auto-dispatch regression: symmetry closes new cases while the residue
 * families keep their exact closed forms (symmetry runs after residue).
 * ---------------------------------------------------------------------- */
static void test_dispatch_regression(void) {
    /* New coverage via symmetry. */
    check_eq("Integrate[x^3/(1+x^2), {x, -a, a}]", "0");
    check_eq("Integrate[Sin[x]^3, {x, -Pi, Pi}]", "0");
    check_eq("Integrate[x^2, {x, -1, 1}]", "2/3");
    /* Even Gaussian on the whole line: symmetry -> 2*half unlocks Ramanujan. */
    check_eq("Chop[N[Integrate[x^2 Exp[-x^2], {x, -Infinity, Infinity}] "
             "- Sqrt[Pi]/2]]", "0");
    /* Residue families unchanged (must NOT regress to an uglier form). */
    check_eq("Integrate[1/(1+x^2), {x, -Infinity, Infinity}]", "Pi");
    check_eq("Integrate[1/(1+x^4), {x, -Infinity, Infinity}]", "Pi/Sqrt[2]");
    /* Divergent odd stays divergent (unevaluated), not 0. */
    check_eq("Integrate[x/(1+x^2), {x, -Infinity, Infinity}]",
             "Integrate[x/(1 + x^2), {x, -Infinity, Infinity}]");
}

/* -------------------------------------------------------------------------
 * Method -> "Symmetry" pins the mechanism (strict: no fallback).
 * ---------------------------------------------------------------------- */
static void test_method_option(void) {
    check_eq("Integrate[x^3, {x, -a, a}, Method -> \"Symmetry\"]", "0");
    check_eq("Integrate[x^5 - x, {x, -2, 2}, Method -> \"Symmetry\"]", "0");
    /* Strict: an asymmetric interval is not this mechanism's job -> unevaluated. */
    check_eq("Integrate[x^2, {x, 0, 1}, Method -> \"Symmetry\"]",
             "Integrate[x^2, {x, 0, 1}, Method -> \"Symmetry\"]");
}

/* -------------------------------------------------------------------------
 * Stress: high odd/even powers, mixed odd*even, symbolic parameters.
 * ---------------------------------------------------------------------- */
static void test_stress(void) {
    check_eq("Integrate`Symmetry[x^11 - 7 x^9 + 2 x^3 - x, {x, -3, 3}]", "0");
    check_eq("Integrate`Symmetry[x^2 Sin[x], {x, -Pi, Pi}]", "0"); /* even*odd = odd */
    check_eq("Integrate`Symmetry[x^4 + x^2 + 1, {x, -2, 2}]", "332/15");
    /* Odd, convergent infinite (rational, deg gap 3) -> 0. */
    check_eq("Integrate`Symmetry[x^3/(1+x^6), {x, -Infinity, Infinity}]", "0");
    /* Even Gaussian on the whole line -> 2*half via Ramanujan (numeric pin). */
    check_eq("Chop[N[Integrate`Symmetry[Exp[-x^2], {x, -Infinity, Infinity}] "
             "- Sqrt[Pi]]]", "0");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_unit_odd);
    TEST(test_unit_even);
    TEST(test_unit_negative);
    TEST(test_dispatch_regression);
    TEST(test_method_option);
    TEST(test_stress);

    printf("All integrate_symmetry tests passed.\n");
    return 0;
}
