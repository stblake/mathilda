/*
 * test_factor_baseline.c
 * ----------------------
 * Performance and correctness baseline for the multivariate factoring
 * pipeline.  Each case runs a Factor or Simplify, verifies the result
 * matches the expected form, and asserts a generous wall-clock budget.
 *
 * This file is the integration point for the cumulative Factor work
 * (Phases 0-6 + tactical fixes).  When the budgets here become too
 * tight, treat that as a regression signal: either the algorithm
 * changed shape (in which case update the budget after auditing) or
 * something genuinely slowed down (investigate).
 *
 * Budgets are intentionally generous (3-5x typical observed timings)
 * to remain stable across machines.  Reference timings on a 2020 x86
 * Mac (gcc -O3): see the comment by each test.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERTION FAILED: %s\n  at %s:%d\n", \
                #cond, __FILE__, __LINE__); \
        abort(); \
    } \
} while (0)

/* Evaluate `input`, time it, and check that the FullForm of the result
 * matches `expected`.  Returns the elapsed milliseconds. */
static double eval_check(const char* input, const char* expected,
                         double budget_ms) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    clock_t t0 = clock();
    Expr* result = evaluate(parsed);
    clock_t t1 = clock();
    double ms = 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC;

    char* got = expr_to_string_fullform(result);
    if (strcmp(got, expected) != 0) {
        fprintf(stderr, "FAIL [correctness]: %s\n"
                        "  expected: %s\n"
                        "  got:      %s\n  (%.2f ms)\n",
                input, expected, got, ms);
        free(got); expr_free(result); expr_free(parsed);
        ASSERT(0);
    }
    if (ms > budget_ms) {
        fprintf(stderr, "FAIL [perf]: %s\n"
                        "  budget: %.2f ms, actual: %.2f ms\n",
                input, budget_ms, ms);
        free(got); expr_free(result); expr_free(parsed);
        ASSERT(0);
    }
    free(got); expr_free(result); expr_free(parsed);
    return ms;
}

/* ====================================================================== */
/*  (1) Direct Factor calls -- correctness + perf                         */
/* ====================================================================== */

/* Univariate Factor: should be sub-millisecond. */
static void test_univariate_factors(void) {
    eval_check("Factor[x^2 - 1]",
               "Times[Plus[-1, x], Plus[1, x]]", 50.0);
    eval_check("Factor[x^3 - 6 x^2 + 11 x - 6]",
               "Times[Plus[-3, x], Plus[-2, x], Plus[-1, x]]", 50.0);
    eval_check("Factor[x^10 - 1]",
               "Times[Plus[-1, x], Plus[1, x], "
               "Plus[1, Times[-1, x], Power[x, 2], Times[-1, Power[x, 3]], Power[x, 4]], "
               "Plus[1, x, Power[x, 2], Power[x, 3], Power[x, 4]]]", 100.0);
}

/* Regression: an UNEXPANDED Power[Plus[...], n] with numeric coefficients must not
 * lose its leading term.  flint_univariate_factor previously computed the degree of
 * the raw (un-expanded) argument, and get_degree_poly counts Power[base,k] only when
 * base IS the variable — so (5-3s)^2 reported degree 0, capping the whole degree at 1
 * and silently dropping the s^2 coefficient (Factor[(5-3s)^2 - 30 s] -> -5(-5+12s),
 * a wrong lower-degree answer).  Factor now expands before computing the degree. */
static void test_unexpanded_power_input(void) {
    eval_check("Factor[(5 - 3 s)^2 - 30 s]",                 /* irreducible over Q */
               "Plus[25, Times[-60, s], Times[9, Power[s, 2]]]", 50.0);
    eval_check("Factor[(1 - s)^2 - (1 - s)]",
               "Times[Plus[-1, s], s]", 50.0);
    eval_check("Factor[(1 + 2 s)^2 - (1 + 2 s)]",
               "Times[2, s, Plus[1, Times[2, s]]]", 50.0);
}

/* Bivariate via factor_binomial / monomial-content (Phase 0). */
static void test_bivariate_simple(void) {
    eval_check("Factor[x^2 - 4 y^2]",
               "Times[Plus[x, Times[-2, y]], Plus[x, Times[2, y]]]", 50.0);
    eval_check("Factor[x^3 + y^3]",
               "Times[Plus[x, y], Plus[Power[x, 2], Times[-1, x, y], Power[y, 2]]]", 50.0);
    eval_check("Factor[x^2 y + x y^2]",
               "Times[x, y, Plus[x, y]]", 50.0);
}

/* Bivariate via Hensel lift (Phase 5a) -- previously stalled. */
static void test_bivariate_hensel(void) {
    eval_check("Factor[(x + y - 1)(x + y + 1)]",
               "Times[Plus[-1, x, y], Plus[1, x, y]]", 100.0);
    eval_check("Factor[(x^2 + y - 1)(x + y + 2)]",
               "Times[Plus[2, x, y], Plus[-1, Power[x, 2], y]]", 100.0);
    eval_check("Factor[(x - 1)(x + 1)(x + y)]",
               "Times[Plus[-1, x], Plus[1, x], Plus[x, y]]", 100.0);
}

/* Multivariate via factor_roots fallback (still works after Phase 0+1). */
static void test_multivariate_via_factor_roots(void) {
    /* = (a - x)(a + x)(3 a^2 - 2 x y).  Wolfram normalises each irreducible
     * factor to a positive leading coefficient in the FIRST (alphabetically
     * smallest) variable it contains — a for (a - x), a^2 for (3 a^2 - 2 x y)
     * — giving a-leading factors, not the x-leading associates. */
    eval_check("Factor[2 x^3 y - 2 a^2 x y - 3 a^2 x^2 + 3 a^4]",
               "Times[Plus[a, Times[-1, x]], Plus[a, x], "
               "Plus[Times[3, Power[a, 2]], Times[-2, x, y]]]", 100.0);
}

/* Irreducibility short-circuit (Phase 1) -- 60x speedup vs main. */
static void test_irreducibility_short_circuit(void) {
    eval_check("Factor[x^2 + y^2 + 1]",
               "Plus[1, Power[x, 2], Power[y, 2]]", 50.0);
    eval_check("Factor[x^4 + y^4 + 1]",
               "Plus[1, Power[x, 4], Power[y, 4]]", 50.0);
    eval_check("Factor[3 x^2 - 3 - y^2]",
               "Plus[-3, Times[3, Power[x, 2]], Times[-1, Power[y, 2]]]", 50.0);
}

/* Threading over comparison heads (Phase 5d, 2nd attempt). */
static void test_factor_threads_over_logic_heads(void) {
    /* Cubic inside Less: factor each comparison argument. */
    eval_check("Factor[1 < 1 + 2 x + x^2 + 1/(1+x) < 2]",
               "Inequality[1, Less, "
               "Times[Power[Plus[1, x], -1], Plus[2, x], Plus[1, x, Power[x, 2]]], "
               "Less, 2]",
               150.0);
}

/* Denominator factoring with extra variables (Phase 5d, context-aware). */
static void test_factor_denominator_with_extra_vars(void) {
    eval_check("Factor[(x^3+2 x^2)/(x^2-4 y^2) - (x+2)/(x^2-4 y^2)]",
               "Times[Plus[-1, x], Plus[1, x], Plus[2, x], "
                     "Power[Plus[x, Times[-2, y]], -1], "
                     "Power[Plus[x, Times[2, y]], -1]]",
               150.0);
}

/* A PRODUCT MUST NOT COME BACK AS A SUM.
 *
 * Factor normalises through Together, and Together expands a multivariate
 * product: Together[E^x (x^2 - 1)] is -E^x + x^2 E^x. That expansion is not a
 * polynomial in the one variable collect_variables finds (x) -- its coefficients
 * carry E^x -- so there was nothing to factor and the expansion came straight
 * back: Factor[E^x (x^2 - 1)] returned -E^x + x^2 E^x, LESS factored than its
 * own input. Reported from the notebook as `Integrate[x^5 E^x, x] // Factor`.
 *
 * Factoring is multiplicative, so a denominator-free product is factored factor
 * by factor. The denominator case must NOT take that route, which the last case
 * pins: Factor[(x^2 - 1)/(x - 1)] is legitimately the sum 1 + x, and threading
 * over that Times would hand back an uncancelled (x - 1)(x + 1)/(x - 1). */
static void test_factor_of_a_product_stays_a_product(void) {
    eval_check("Factor[E^x (x^2 - 1)]",
               "Times[Power[E, x], Plus[-1, x], Plus[1, x]]", 200.0);
    /* The reported case: the quintic is irreducible over Q, so the answer is the
       input -- but it must be the input, not the input multiplied out. */
    eval_check("Factor[E^x (x^5 - 5 x^4 + 20 x^3 - 60 x^2 + 120 x - 120)]",
               "Times[Power[E, x], Plus[-120, Times[120, x], Times[-60, Power[x, 2]], "
               "Times[20, Power[x, 3]], Times[-5, Power[x, 4]], Power[x, 5]]]", 200.0);
    eval_check("Factor[E^x (x^2 + 2 x + 1)]",
               "Times[Power[E, x], Power[Plus[1, x], 2]]", 200.0);
    /* Not exponential-specific: any factor the polynomial path cannot see. */
    eval_check("Factor[Log[x] (x^2 - 1)]",
               "Times[Plus[-1, x], Plus[1, x], Log[x]]", 200.0);
    eval_check("Factor[y (x^2 - 1)]",
               "Times[Plus[-1, x], Plus[1, x], y]", 200.0);
    /* A denominator is the excluded case: a sum here is the right answer. */
    eval_check("Factor[(x^2 - 1)/(x - 1)]", "Plus[1, x]", 200.0);
}

/* ====================================================================== */
/*  (2) Simplify integration -- the user-reported workloads               */
/* ====================================================================== */

/* The user's primary case.  Reference: ~478 ms.  Pre-Phase-0 it was
 * ~1158 ms (2.4x improvement). */
static void test_simplify_user_case(void) {
    double ms = eval_check("Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]",
                           "Times[-3, Power[Sin[x], 3]]",
                           1500.0);
    fprintf(stderr, "  user case: %.0f ms\n", ms);
}

/* Sensitive Tan case: rolled back v1 of 5d caused a 26x slowdown here.
 * After context-aware fix (Phase 5d v2) this stays around 2.5 s,
 * which is faster than pristine main (3.3 s) due to Factor memo. */
static void test_simplify_tan_double_angle(void) {
    double ms = eval_check("Simplify[Tan[2 x] - 2 Tan[x]/(1 - Tan[x]^2)]",
                           "0", 5000.0);
    fprintf(stderr, "  Tan double-angle: %.0f ms\n", ms);
}

/* Pythagorean to-the-fourth.  Reference: ~80 ms. */
static void test_simplify_pythagorean_power(void) {
    double ms = eval_check("Simplify[(Sin[x]^2 + Cos[x]^2)^4]",
                           "1", 500.0);
    fprintf(stderr, "  Pythag^4: %.0f ms\n", ms);
}

/* Hyperbolic Pythagorean.  Should be near-instant. */
static void test_simplify_hyperbolic_pythag(void) {
    double ms = eval_check("Simplify[Sinh[x]^2 - Cosh[x]^2]",
                           "-1", 200.0);
    fprintf(stderr, "  Sinh^2 - Cosh^2: %.0f ms\n", ms);
}

/* Identity that needs Factor working on the cubic 4 Cos^3 - 3 Cos. */
static void test_simplify_chebyshev_cubic(void) {
    double ms = eval_check("Simplify[Cos[3 x]/(4 Cos[x]^3 - 3 Cos[x])]",
                           "1", 500.0);
    fprintf(stderr, "  Cos[3x]/Cheb: %.0f ms\n", ms);
}

/* Compound expansion: (x+y)^3 - (x-y)^3 = 6 x^2 y + 2 y^3.
 * Reference: ~46 ms. */
static void test_simplify_cubic_expansion(void) {
    double ms = eval_check("Simplify[(x + y)^3 - (x - y)^3]",
                           "Times[2, y, Plus[Times[3, Power[x, 2]], Power[y, 2]]]",
                           300.0);
    fprintf(stderr, "  (x+y)^3 - (x-y)^3: %.0f ms\n", ms);
}

/* Sin double angle through division. */
static void test_simplify_sin2x_over_sinx(void) {
    double ms = eval_check("Simplify[Sin[2 x] / Sin[x]]",
                           "Times[2, Cos[x]]", 500.0);
    fprintf(stderr, "  Sin[2x]/Sin[x]: %.0f ms\n", ms);
}

/* ====================================================================== */
/*  (3) Sanity: previously slow / known-tricky cases                     */
/* ====================================================================== */

/* The user's reported pre-Phase-0 case.  Was returning input
 * unfactored AND in 364 ms; now returns the correct factorisation in
 * single-digit ms. */
static void test_factor_monomial_cube(void) {
    double ms = eval_check(
        "Factor[3 a^2 b - 3 b - b^3]",
        "Times[b, Plus[-3, Times[3, Power[a, 2]], Times[-1, Power[b, 2]]]]",
        100.0);
    fprintf(stderr, "  3a^2 b - 3b - b^3: %.0f ms\n", ms);
}

/* Trig form of the same shape, derived from Simplify intermediate. */
static void test_factor_trig_monomial(void) {
    double ms = eval_check(
        "Factor[3 Cos[x]^2 Sin[x] - 3 Sin[x] - Sin[x]^3]",
        "Times[Sin[x], Plus[-3, Times[3, Power[Cos[x], 2]], "
                          "Times[-1, Power[Sin[x], 2]]]]",
        100.0);
    fprintf(stderr, "  trig monomial: %.0f ms\n", ms);
}

/* Opaque kernels are independent generators.  The >=3-generator z-split
 * step specialised x -> 0 by ReplaceAll, which also rewrote the x INSIDE
 * BesselJ[-2/3, x] (= ComplexInfinity at 0) and leaked Infinity::indet out
 * of Factor and Simplify; f[x] was silently specialised to f[0].  Factor
 * now maps kernels to fresh symbols first.  Each case: Check returns the
 * value (no message) and the result is equal to the input. */
static void factor_kernel_case(const char* poly) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Check[Expand[Factor[%s] - (%s)], $Failed]", poly, poly);
    assert_eval_eq(buf, "0", 0);
}

/* Gaussian-integer coefficients: the content is taken over Z[i] (its
 * first-quadrant associate), times I when the rest is purely imaginary, and
 * the primitive part's sign goes into it -- the reference's forms. */
static void test_factor_gaussian_content(void) {
    assert_eval_eq("{Factor[4 I a - 4 I b], Factor[2 I x + 2 I], Factor[-I a + I b], "
                   "Factor[(1 + I) a + 2 b], Factor[2 I a + 4 b], Factor[I x^2 - I], "
                   "Factor[-2 I a - 2 I b], Factor[(2 + 2 I) a + (4 + 4 I) b], Factor[I a + b], "
                   "Factor[56 I - 30 I d + 3 I d^2], Factor[I a b - I c d], FactorList[4 I a - 4 I b]}",
                   "{(4*I) (a - b), (2*I) (1 + x), -I (a - b), (1 + I) (a + (1 - I) b), "
                   "2 (I a + 2 b), I (-1 + x) (1 + x), -(2*I) (a + b), (2 + 2*I) (a + 2 b), I a + b, "
                   "I (56 - 30 d + 3 d^2), I (a b - c d), {{4*I, 1}, {a - b, 1}}}", 0);
}

/* Issue #114: a multivariate polynomial with Gaussian coefficients is factored
 * over Q(I) (FLINT declines it), and the classical path took 63 s on this
 * irreducible 28-term, 30-kernel input.  Degree one in e[...] with coprime
 * coefficients settles irreducibility at once (Gauss's lemma). */
static const char* ISSUE114_POLY =
        "I e[l[i2], l[i3], l[i6], l[i7]] p[l[i1, n], l[i5, n]] + 2 p[l[i2], l[i7]] p[l[i3], l[i6]] "
        "p[l[i1, n], l[i5, n]] - 2 p[l[i2], l[i6]] p[l[i3], l[i7]] p[l[i1, n], l[i5, n]] + 2 i "
        "e[l[i2], l[i3], l[i5], l[i7]] p[l[i1, n], l[i6, n]] - 2 p[l[i2], l[i7]] p[l[i3], l[i5]] "
        "p[l[i1, n], l[i6, n]] + 2 p[l[i2], l[i5]] p[l[i3], l[i7]] p[l[i1, n], l[i6, n]] + I "
        "e[l[i2], l[i3], l[i5], l[i6]] p[l[i1, n], l[i7, n]] + 2 p[l[i2], l[i6]] p[l[i3], l[i5]] "
        "p[l[i1, n], l[i7, n]] - 2 p[l[i2], l[i5]] p[l[i3], l[i6]] p[l[i1, n], l[i7, n]] + I "
        "e[l[i1], l[i3], l[i6], l[i7]] p[l[i2, n], l[i5, n]] - 2 p[l[i1], l[i7]] p[l[i3], l[i6]] "
        "p[l[i2, n], l[i5, n]] + 2 p[l[i1], l[i6]] p[l[i3], l[i7]] p[l[i2, n], l[i5, n]] - 2 "
        "p[l[i6], l[i7]] p[l[i1, 4 + n], l[i3, 4 + n]] p[l[i2, n], l[i5, n]] + I e[l[i1], l[i3], "
        "l[i5], l[i7]] p[l[i2, n], l[i6, n]] + 2 p[l[i1], l[i7]] p[l[i3], l[i5]] p[l[i2, n], l[i6, "
        "n]] - 2 p[l[i1], l[i5]] p[l[i3], l[i7]] p[l[i2, n], l[i6, n]] + 2 p[l[i5], l[i7]] p[l[i1, "
        "4 + n], l[i3, 4 + n]] p[l[i2, n], l[i6, n]] + I e[l[i1], l[i3], l[i5], l[i6]] p[l[i2, n], "
        "l[i7, n]] - 2 p[l[i1], l[i6]] p[l[i3], l[i5]] p[l[i2, n], l[i7, n]] + 2 p[l[i1], l[i5]] "
        "p[l[i3], l[i6]] p[l[i2, n], l[i7, n]] - 2 p[l[i5], l[i6]] p[l[i1, 4 + n], l[i3, 4 + n]] "
        "p[l[i2, n], l[i7, n]] + 2 p[l[i6], l[i7]] p[l[i1, n], l[i5, n]] p[l[i2, 4 + n], l[i3, 4 + "
        "n]] - 2 p[l[i5], l[i7]] p[l[i1, n], l[i6, n]] p[l[i2, 4 + n], l[i3, 4 + n]] + 2 p[l[i5], "
        "l[i6]] p[l[i1, n], l[i7, n]] p[l[i2, 4 + n], l[i3, 4 + n]] + I e[l[i1], l[i2], l[i6], "
        "l[i7]] p[l[i3, n], l[i5, n]] + 2 p[l[i1], l[i7]] p[l[i2], l[i6]] p[l[i3, n], l[i5, n]] - 2"
        " p[l[i1], l[i6]] p[l[i2], l[i7]] p[l[i3, n], l[i5, n]] + 2 p[l[i6], l[i7]] p[l[i1, 4 + n],"
        " l[i2, 4 + n]] p[l[i3, n], l[i5, n]]";

static void test_factor_gaussian_multivariate(void) {
    char buf[8192];
    snprintf(buf, sizeof buf, "Factor[%s] === %s", ISSUE114_POLY, ISSUE114_POLY);
    eval_check(buf, "True", 2000.0);
    /* factorable Gaussian products still factor over Q(I) */
    assert_eval_eq("Factor[Expand[(x + I)^2 (y + 1)]]", "(I + x)^2 (1 + y)", 0);
    assert_eval_eq("Factor[Expand[(x + I y) (x - 2 I z + y)]]", "(x + I y) (x + y - (2*I) z)", 0);
    assert_eval_eq("Factor[Expand[(a + I b + f[c]) (a - I b + 2 g[d])]]",
                   "(a + I b + f[c]) (a - I b + 2 g[d])", 0);
    /* an algebraic constant is not an indeterminate: no degree-one shortcut */
    assert_eval_eq("Factor[x^2 + Sqrt[3] x - 6]", "-6 + Sqrt[3] x + x^2", 0);
}

static void test_factor_kernels_are_opaque(void) {
    factor_kernel_case("x^2 BesselJ[-2/3, x]^2 + x BesselJ[1/3, x]^2 + 1");
    factor_kernel_case("x^2 Log[x]^2 + x Log[2 x]^2 + 1");
    factor_kernel_case("32 BesselJ[1/3, x]^3 + 81 x^4 BesselJ[1/3, x]^3 + 324 x^3 BesselJ[4/3, x]^3"
        " - 360 x^2 BesselJ[1/3, x] BesselJ[4/3, x]^2 - 162 x^4 BesselJ[7/3, x] BesselJ[4/3, x]^2"
        " - 108 x^3 BesselJ[-2/3, x] BesselJ[1/3, x]^2 - 108 x^3 BesselJ[-2/3, x] BesselJ[4/3, x]^2"
        " - 81 x^4 BesselJ[7/3, x] BesselJ[1/3, x]^2 - 72 x^2 BesselJ[7/3, x] BesselJ[1/3, x]^2"
        " - 48 x BesselJ[-2/3, x] BesselJ[1/3, x]^2 + 48 x BesselJ[4/3, x] BesselJ[1/3, x]^2"
        " + 216 x^3 BesselJ[4/3, x] BesselJ[1/3, x]^2 + 144 x^2 BesselJ[-2/3, x] BesselJ[1/3, x] BesselJ[4/3, x]"
        " + 162 x^4 BesselJ[-2/3, x] BesselJ[1/3, x] BesselJ[4/3, x]"
        " + 216 x^3 BesselJ[1/3, x] BesselJ[4/3, x] BesselJ[7/3, x]");
    /* genuinely factorable over kernels */
    assert_eval_eq("Factor[x^2 a[x]^2 - b[x]^2]", "(x a[x] - b[x]) (x a[x] + b[x])", 0);
    assert_eval_eq("Factor[x^2 Sin[x] + 2 x Sin[x] + Sin[x]]", "(1 + x)^2 Sin[x]", 0);
    assert_eval_eq("Factor[E^x (x^2 - 1)]", "E^x (-1 + x) (1 + x)", 0);
}

/* User-reported case (2026-04): a trivariate 18-term polynomial whose
 * FactorSquareFree pre-check ran a univariate gcd over a degree-31
 * image.  The original Expr-level `poly_gcd_internal` (Knuth-style
 * primitive PRS) suffered exponential coefficient growth and took
 * >120s to terminate.  The fix routes `univariate_squarefree`'s
 * internal gcd through `zupoly_gcd` (subresultant PRS), bringing the
 * full Factor call to under 1.5 s.  Mathematica solves it in 2 ms,
 * Mathics in 100 ms; we're not at parity but we're no longer hung.
 *
 * Expected: x^2 * (z^13 - x^12) * (17 - 5y - z^14) * (z^4 + 3x^9 - y^13).
 * Each irreducible factor is normalised to a positive lexicographic leading
 * term (variables in canonical order, x first), so (17 - 5y - z^14) is
 * emitted as -17 + 5y + z^14 and (z^13 - x^12) as x^12 - z^13, the two sign
 * flips cancelling overall. */
static void test_factor_trivariate_high_degree_squarefree(void) {
    double ms = eval_check(
        "Factor[Expand[x^2 (z^13 - x^12) (z^4 + 3 x^9 - y^13) (17 - 5 y - z^14)]]",
        "Times[Power[x, 2], "
        "Plus[Times[3, Power[x, 9]], Times[-1, Power[y, 13]], Power[z, 4]], "
        "Plus[Power[x, 12], Times[-1, Power[z, 13]]], "
        "Plus[-17, Times[5, y], Power[z, 14]]]",
        5000.0);
    fprintf(stderr, "  trivariate user case: %.0f ms\n", ms);
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running factor baseline tests...\n");

    /* Direct Factor */
    TEST(test_univariate_factors);
    TEST(test_unexpanded_power_input);
    TEST(test_bivariate_simple);
    TEST(test_bivariate_hensel);
    TEST(test_multivariate_via_factor_roots);
    TEST(test_irreducibility_short_circuit);
    TEST(test_factor_threads_over_logic_heads);
    TEST(test_factor_denominator_with_extra_vars);
    TEST(test_factor_of_a_product_stays_a_product);

    /* Simplify integration */
    TEST(test_simplify_user_case);
    TEST(test_simplify_tan_double_angle);
    TEST(test_simplify_pythagorean_power);
    TEST(test_simplify_hyperbolic_pythag);
    TEST(test_simplify_chebyshev_cubic);
    TEST(test_simplify_cubic_expansion);
    TEST(test_simplify_sin2x_over_sinx);

    /* Direct Factor regression cases */
    TEST(test_factor_monomial_cube);
    TEST(test_factor_trig_monomial);
    TEST(test_factor_trivariate_high_degree_squarefree);
    TEST(test_factor_kernels_are_opaque);
    TEST(test_factor_gaussian_multivariate);
    TEST(test_factor_gaussian_content);

    printf("All baseline tests passed!\n");
    return 0;
}
