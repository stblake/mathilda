/* Unit tests for Minimize / Maximize (src/calculus/minimize.c).
 *
 * Minimize is the EXACT (symbolic) global optimizer. Exact input yields exact
 * output; inexact input delegates to NMinimize. The engine is sound by
 * construction: when a case cannot be decided exactly it leaves the expression
 * unevaluated rather than guessing, so several tests assert a *decline*
 * (Head[...] === Minimize) as the correct behaviour.
 *
 * Coverage:
 *   - Result shape / typing; {f_opt, {x -> x_opt, ...}}.
 *   - Univariate polynomial: quadratic, quartic, constant, high-degree, the
 *     degree-6 symmetric showcase, unbounded (odd degree / negative leading).
 *   - Maximize via objective negation + min/max duality.
 *   - Multivariate unconstrained polynomial with isolated minima.
 *   - Constrained polynomial over the Reals (disk, equality+inequality, LP).
 *   - Infeasible region -> {Infinity, {x -> Indeterminate, ...}}.
 *   - Inexact input -> numeric NMinimize fallback.
 *   - Sound declines: transcendental exact input, positive-dimensional minima,
 *     unbounded-constrained, bad arity.
 *   - Robustness & scope expansion (v0.305-v0.310): internal time budget
 *     (torus solves, no hang), compact-region + radical-cleared certificates
 *     (irrational-algebraic optima), radical/fractional-power objectives, Abs /
 *     piecewise objectives, bounded Integers-domain optimization, TimeConstraint
 *     option, and the Minimize::natt message routing (Check[]-catchable).
 *   - Memory-hygiene smoke loop.
 *
 * Run binary directly: ./minimize_tests */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* Parse, evaluate, FullForm-compare. */
static void check_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected) != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Parse+evaluate a boolean predicate and require it to be True. */
static void check_true(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, "True") != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: True\n  got:      %s\n",
                input, got);
        ASSERT_STR_EQ(got, "True");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* 1. Result shape                                                     */
/* ------------------------------------------------------------------ */

static void test_result_shape(void) {
    check_eq("Head[Minimize[x^2, x]]", "List");
    check_eq("Length[Minimize[x^2, x]]", "2");
    check_eq("Head[Last[Minimize[x^2, x]]]", "List");
    check_eq("Head[First[Last[Minimize[x^2, x]]]]", "Rule");
}

/* ------------------------------------------------------------------ */
/* 2. Univariate unconstrained polynomial                              */
/* ------------------------------------------------------------------ */

static void test_quadratic(void) {
    check_eq("Minimize[2 x^2 - 3 x + 5, x]",
             "List[Rational[31, 8], List[Rule[x, Rational[3, 4]]]]");
}

static void test_quadratic_simple(void) {
    check_eq("Minimize[x^2, x]", "List[0, List[Rule[x, 0]]]");
}

static void test_shifted_parabola(void) {
    check_eq("Minimize[(x - 3)^2 + 1, x]", "List[1, List[Rule[x, 3]]]");
}

static void test_constant(void) {
    check_eq("Minimize[5, x]", "List[5, List[Rule[x, 0]]]");
}

static void test_quartic_value(void) {
    /* x^4 - 2 x^2 -> global min -1 at x = ±1 (one point returned). */
    check_true("First[Minimize[x^4 - 2 x^2, x]] == -1");
    check_true("(x /. Last[Minimize[x^4 - 2 x^2, x]])^2 == 1");
}

static void test_degree6_showcase(void) {
    /* The classic exact-global-min showcase; minimum -35721/64 at x = 11/2. */
    check_eq("Minimize[Expand[(x-1)(x-2)(x-4)(x-7)(x-9)(x-10)], x]",
             "List[Rational[-35721, 64], List[Rule[x, Rational[11, 2]]]]");
}

static void test_unbounded_odd(void) {
    check_true("Minimize[x^3, x] === {-Infinity, {x -> Indeterminate}}");
}

static void test_unbounded_negative_leading(void) {
    check_true("Minimize[-x^2, x] === {-Infinity, {x -> Indeterminate}}");
}

/* ------------------------------------------------------------------ */
/* 3. Maximize                                                         */
/* ------------------------------------------------------------------ */

static void test_maximize_parabola(void) {
    check_eq("Maximize[-x^2, x]", "List[0, List[Rule[x, 0]]]");
}

static void test_maximize_quadratic(void) {
    /* max of -(2x^2-3x+5) = -(min) at the same point. */
    check_eq("Maximize[-(2 x^2 - 3 x + 5), x]",
             "List[Rational[-31, 8], List[Rule[x, Rational[3, 4]]]]");
}

static void test_maximize_unbounded(void) {
    check_true("Maximize[x^2, x] === {Infinity, {x -> Indeterminate}}");
}

static void test_min_max_duality(void) {
    check_true("First[Maximize[-(x - 2)^2 - 1, x]] == -1");
    check_true("(x /. Last[Maximize[-(x - 2)^2 - 1, x]]) == 2");
}

/* ------------------------------------------------------------------ */
/* 4. Multivariate unconstrained polynomial                            */
/* ------------------------------------------------------------------ */

static void test_bowl_2d(void) {
    check_eq("Minimize[x^2 + y^2, {x, y}]",
             "List[0, List[Rule[x, 0], Rule[y, 0]]]");
}

static void test_shifted_bowl_2d(void) {
    check_eq("Minimize[(x - 1)^2 + (y + 2)^2, {x, y}]",
             "List[0, List[Rule[x, 1], Rule[y, -2]]]");
}

static void test_positive_dimensional_solves(void) {
    /* Minimizer set is the curve x^2 - 2y = 1/2 (positive-dimensional): the
     * critical-point method finds no isolated candidate (Solve::nsdim), so the
     * infimum is read off real QE (Reduce[ForAll[...]]) and a minimiser is
     * realised with FindInstance. f = u^2 - u - 1 with u = x^2 - 2y, min -5/4
     * at u = 1/2. Assert the value and that the reported point attains it. */
    check_true("First[Minimize[(x^2 - 2 y)^2 - x^2 + 2 y - 1, {x, y}]] == -5/4");
    check_true("(((x^2 - 2 y)^2 - x^2 + 2 y - 1) /. "
               "Last[Minimize[(x^2 - 2 y)^2 - x^2 + 2 y - 1, {x, y}]]) == -5/4");
}

static void test_flat_valley_hyperbola(void) {
    /* The originally-reported case: f = (x y - 3)^2 + 1 is 1 everywhere on the
     * hyperbola x y == 3 (a positive-dimensional minimum). Value 1, attained. */
    check_true("First[Minimize[(x y - 3)^2 + 1, {x, y}]] == 1");
    check_true("(((x y - 3)^2 + 1) /. Last[Minimize[(x y - 3)^2 + 1, {x, y}]]) == 1");
    /* Line valley x + y == 2: min 0, attained. */
    check_true("First[Minimize[(x + y - 2)^2, {x, y}]] == 0");
    check_true("(((x + y - 2)^2) /. Last[Minimize[(x + y - 2)^2, {x, y}]]) == 0");
}

static void test_flat_valley_algebraic(void) {
    /* Quartic flat valley: f = (x y - 3)^4 - x y + 1 is minimised along the
     * hyperbola x y == 3 + 4^(-1/3). The infimum is an algebraic Root, so the
     * witness is realised by minimizing a univariate slice (FindInstance can't
     * instantiate the level set at an algebraic value). Assert the value
     * numerically against -2 - (3/4) 2^(-2/3) and that the reported point
     * attains it exactly. */
    check_true("Module[{f = (x y - 3)^4 - x y + 1, r}, "
               "r = Minimize[f, {x, y}]; "
               "Simplify[(f /. Last[r]) - First[r]] == 0 && "
               "Abs[N[First[r]] - (-2 - (3/4) 2^(-2/3))] < 10^-9]");
    /* A sextic flat valley with a rational infimum (FindInstance path). */
    check_eq("Minimize[(x y - 2)^6 + 3, {x, y}]",
             "List[3, List[Rule[x, -1], Rule[y, -2]]]");
}

static void test_positive_dim_unbounded(void) {
    /* Saddle x^2 - y^2: unbounded below along x==0 -> -Infinity, not attained.
     * (The critical point (0,0) is a saddle, so mz_exact_poly declines and the
     * QE fallback detects Reduce[ForAll[...]] === False.) */
    check_eq("Minimize[x^2 - y^2, {x, y}]",
             "List[Times[-1, Infinity], List[Rule[x, Indeterminate], "
             "Rule[y, Indeterminate]]]");
    /* Maximize of a flat valley that is unbounded above -> +Infinity. */
    check_eq("Maximize[(x y - 3)^2 + 1, {x, y}]",
             "List[Infinity, List[Rule[x, Indeterminate], Rule[y, Indeterminate]]]");
    /* Maximize mirror of the hyperbola valley: finite max -1, attained. */
    check_true("First[Maximize[-(x y - 3)^2 - 1, {x, y}]] == -1");
}

/* ------------------------------------------------------------------ */
/* 5. Constrained polynomial over the Reals                            */
/* ------------------------------------------------------------------ */

static void test_disk_linear(void) {
    /* {x + y, x^2 + y^2 <= 1} -> -Sqrt[2] at (-1/Sqrt2, -1/Sqrt2). */
    check_true("Abs[N[First[Minimize[{x + y, x^2 + y^2 <= 1}, {x, y}]]] - N[-Sqrt[2]]] < 1*^-9");
    check_true("Abs[N[(x /. Last[Minimize[{x + y, x^2 + y^2 <= 1}, {x, y}]])] + N[1/Sqrt[2]]] < 1*^-9");
}

static void test_equality_inequality(void) {
    /* {x + 2y, x^2 + 2y^2 <= 3 && x + y == 2 && x >= 1} -> 7/3 at (5/3, 1/3). */
    check_eq("Minimize[{x + 2 y, x^2 + 2 y^2 <= 3 && x + y == 2 && x >= 1}, {x, y}]",
             "List[Rational[7, 3], List[Rule[x, Rational[5, 3]], Rule[y, Rational[1, 3]]]]");
}

static void test_linear_program(void) {
    /* {2x + 3y - z, 1<=x+y+z<=2 && 1<=x-y+z<=2 && x-y-z==3} -> 3. */
    check_true("First[Minimize[{2 x + 3 y - z, "
               "1 <= x + y + z <= 2 && 1 <= x - y + z <= 2 && x - y - z == 3}, "
               "{x, y, z}]] == 3");
}

static void test_infeasible(void) {
    check_true("Minimize[{x + y, x^2 < -1}, {x, y}] === "
               "{Infinity, {x -> Indeterminate, y -> Indeterminate}}");
}

static void test_univariate_constrained(void) {
    /* {3x^2 - x + 9, 2x^3 + 5x - 7 >= 0} -> 11 at x = 1. */
    check_eq("Minimize[{3 x^2 - x + 9, 2 x^3 + 5 x - 7 >= 0}, x]",
             "List[11, List[Rule[x, 1]]]");
}

/* ------------------------------------------------------------------ */
/* 6. Inexact input -> NMinimize fallback                              */
/* ------------------------------------------------------------------ */

static void test_inexact_fallback(void) {
    /* 2.5 x^2 - 3 x -> min -0.9 at x = 0.6 (numeric). */
    check_true("Abs[First[Minimize[2.5 x^2 - 3 x, x]] - (-0.9)] < 1*^-3");
    check_true("Head[First[Minimize[2.5 x^2 - 3 x, x]]] === Real");
}

/* ------------------------------------------------------------------ */
/* 7. Sound declines                                                   */
/* ------------------------------------------------------------------ */

static void test_transcendental_declines(void) {
    check_eq("Head[Minimize[Sin[x] + Cos[x], x]]", "Minimize");
}

static void test_unbounded_constrained_declines(void) {
    /* Unbounded below over the feasible region: the constrained QE infimum path
     * (Reduce[ForAll[{x,y}, x<=y^2 => x+y>=b]] === False) now settles this
     * exactly as -Infinity rather than declining. */
    check_true("Minimize[{x + y, x <= y^2}, {x, y}] === "
               "{-Infinity, {x -> Indeterminate, y -> Indeterminate}}");
}

static void test_bad_arity(void) {
    check_eq("Head[Minimize[x^2]]", "Minimize");
}

/* ------------------------------------------------------------------ */
/* 9. Robustness & scope expansion (v0.305 - v0.310), driven by a       */
/*    20-problem stress suite. Each guards one new capability or one     */
/*    deliberate decline against silently regressing.                    */
/* ------------------------------------------------------------------ */

/* v0.306 compact-region shortcut: an irrational-algebraic optimum on a circle
 * (the Reduce certificate alone cannot compare against 14 - 2 Sqrt[13]). */
static void test_circle_irrational(void) {
    check_true("Simplify[First[Minimize[{x^2 + y^2, (x-2)^2 + (y-3)^2 == 1}, "
               "{x, y}]] == 14 - 2 Sqrt[13]]");
}

/* v0.305 internal time budget: the torus case that used to hang 55 s+ now SOLVES
 * (returns a finite {value, rules}) in bounded time. */
static void test_torus_solves(void) {
    check_eq("Head[Minimize[{x + y + z, (x^2 + y^2 + z^2 + 3)^2 == 16 (x^2 + y^2)}, "
             "{x, y, z}]]", "List");
}

/* v0.307 radical-cleared certificate: irrational optimum under an inequality. */
static void test_disk_inequality_irrational(void) {
    check_true("Simplify[First[Minimize[{x^2 + y^2, (x-2)^2 + (y-3)^2 <= 1}, "
               "{x, y}]] == 14 - 2 Sqrt[13]]");
}

/* v0.308 radical objective (well-formed {f, cons}, vars syntax), nested roots. */
static void test_radical_objective(void) {
    check_true("Simplify[First[Minimize[{Sqrt[x + Sqrt[x]] + Sqrt[x - Sqrt[x]], "
               "x >= 1}, x]] == Sqrt[2]]");
}

/* v0.308 limitation: a fractional power whose substitution raises the constraint
 * degree past the CAD (x^(2/3) with x^4 -> a^12) declines, never guesses. */
static void test_fractional_power_declines(void) {
    check_eq("Head[Minimize[{x^(2/3) + y^(2/3), x^4 + y^4 <= 1}, {x, y}, "
             "TimeConstraint -> 2]]", "Minimize");
}

/* v0.309 univariate Abs / piecewise: the sum-of-abs weighted median. */
static void test_abs_sum_median(void) {
    check_eq("First[Minimize[Sum[Abs[x - i^2], {i, 1, 10}], x]]", "275");
}

/* v0.309 Abs piecewise with an interior stationary point of a piece. */
static void test_abs_interior_min(void) {
    check_eq("First[Minimize[x^2 + Abs[x - 2], x]]", "Rational[7, 4]");
}

/* v0.309 Abs objective unbounded below (read off an end piece). */
static void test_abs_unbounded(void) {
    check_true("Minimize[2 x + Abs[x], x] === {-Infinity, {x -> Indeterminate}}");
}

/* v0.310 equality (Diophantine) integer optimum over x^2 + y^2 == 25. */
static void test_integer_equality(void) {
    check_eq("First[Minimize[{x + y, x^2 + y^2 == 25}, {x, y}, Integers]]", "-7");
}

/* v0.310 bounded inequality integer program via box enumeration. */
static void test_integer_box(void) {
    check_eq("First[Minimize[{x^2 + y^2, x + y >= 3 && 0 <= x <= 5 && 0 <= y <= 5}, "
             "{x, y}, Integers]]", "5");
}

/* v0.310 honest scope: an unbounded / hard Diophantine integer region declines. */
static void test_integer_hard_declines(void) {
    check_eq("Head[Minimize[{x^2 + y^2 + z^2, x^3 + y^3 + z^3 == 33}, {x, y, z}, "
             "Integers, TimeConstraint -> 2]]", "Minimize");
}

/* Message routing: Minimize::natt is catchable by Check[] (hence suppressed by
 * Quiet[]). A raw stderr write would make Check[] take the wrong branch. */
static void test_message_routing_check(void) {
    check_true("Check[Minimize[x^3, x], \"caught\"] === \"caught\"");
}

/* The TimeConstraint option is accepted and leaves an easy solve intact. */
static void test_timeconstraint_option(void) {
    check_eq("First[Minimize[x^2, x, TimeConstraint -> 5]]", "0");
}

/* ------------------------------------------------------------------ */
/* 9b. Campaign II (21-40 stress suite): new scope, each with a decline */
/*     sibling that must stay unevaluated (soundness pin).              */
/* ------------------------------------------------------------------ */

/* Separable / additive decomposition (variable-disjoint unconstrained): the
 * global min equals the sum of the univariate block minima. */
static void test_separable_decomposition(void) {
    check_true("Head[Minimize[(x^4 - 16 x^2 + 5 x) + (y^4 - 16 y^2 + 5 y), {x, y}]] === List");
    check_true("Abs[N[First[Minimize[(x^4-16x^2+5x)+(y^4-16y^2+5y),{x,y}]]] "
               "- 2 N[First[Minimize[x^4-16x^2+5x,x]]]] < 10^-6");
}

/* General compact-region shortcut: a rational optimum on a high-degree region
 * (sextic over the simplex) that the Reduce certificate cannot decide. */
static void test_compact_simplex_sextic(void) {
    check_true("First[Minimize[{x^6 + y^6 + z^6 - x y z, x + y + z == 1, "
               "x >= 0, y >= 0, z >= 0}, {x, y, z}]] == -8/243");
}

/* General compact-region shortcut: an irrational Root optimum on a compact
 * sphere-cap-cylinder curve (two equalities). */
static void test_compact_viviani(void) {
    check_true("Head[Minimize[{x + y + z, x^2 + y^2 + z^2 == 4 && "
               "(x - 1)^2 + y^2 == 1}, {x, y, z}]] === List");
    check_true("Abs[N[First[Minimize[{x + y + z, x^2 + y^2 + z^2 == 4 && "
               "(x - 1)^2 + y^2 == 1}, {x, y, z}]]] + 2.3009759965] < 10^-6");
}

/* Equality-constraint variable elimination: x==t, y==t^2, z==t^3 collapses to a
 * univariate problem in t. */
static void test_equality_elimination(void) {
    check_true("Head[Minimize[{(x-1)^2 + (y-2)^2 + (z-3)^2, "
               "x == t && y == t^2 && z == t^3}, {x, y, z, t}]] === List");
    check_true("Abs[N[First[Minimize[{(x-1)^2+(y-2)^2+(z-3)^2, "
               "x==t && y==t^2 && z==t^3}, {x,y,z,t}]]] - 0.1924713154] < 10^-6");
}

/* Rational-function objective: a positive-definite denominator reduces to a
 * polynomial problem in w = p/q. */
static void test_rational_objective(void) {
    check_true("First[Minimize[{x^2 y^2 / (x^2 + y^2 + 1), x^2 + y^2 >= 1}, {x, y}]] == 0");
    check_true("First[Minimize[x^2/(x^2 + 1), x]] == 0");
}

/* Rational objective whose denominator has indefinite sign on the closure
 * (boundary poles): a sound decline. */
static void test_rational_poles_decline(void) {
    check_eq("Head[Minimize[{1/x + 1/y + 1/z, x + y + z == 1, x > 0, y > 0, z > 0}, "
             "{x, y, z}, TimeConstraint -> 3]]", "Minimize");
}

/* Integers domain via the Element[{vars}, Integers] spelling AND an infinite
 * parametric Diophantine family with a coercive objective. */
static void test_integer_parametric(void) {
    check_true("First[Minimize[{x^2 + y^2 + z^2, 5 x + 7 y + 11 z == 13, "
               "Element[{x, y, z}, Integers]}, {x, y, z}]] == 3");
}

/* Mixed integer/continuous (only some variables declared integer) is out of
 * scope and must decline, not be silently treated as all-integer. */
static void test_mixed_integer_declines(void) {
    check_eq("Head[Minimize[{x y - z w, x + y + z + w == 10, 0 <= x <= 5, "
             "0 <= y <= 5, Element[{z, w}, Integers]}, {x, y, z, w}, "
             "TimeConstraint -> 3]]", "Minimize");
}

/* ------------------------------------------------------------------ */
/* 8. Memory smoke                                                     */
/* ------------------------------------------------------------------ */

static void test_memory_smoke(void) {
    for (int i = 0; i < 50; i++) {
        Expr* e = parse_expression("Minimize[2 x^2 - 3 x + 5, x]");
        Expr* r = evaluate(e);
        expr_free(e);
        expr_free(r);
    }
    for (int i = 0; i < 20; i++) {
        Expr* e = parse_expression("Minimize[{x + y, x^2 + y^2 <= 1}, {x, y}]");
        Expr* r = evaluate(e);
        expr_free(e);
        expr_free(r);
    }
}

int main(void) {
    symtab_init();
    core_init();

    /* Intended-error / decline tests write diagnostics to stderr; keep output
     * clean. Successful exact solves are silent. */
    freopen("/dev/null", "w", stderr);

    /* 1. Shape */
    TEST(test_result_shape);

    /* 2. Univariate */
    TEST(test_quadratic);
    TEST(test_quadratic_simple);
    TEST(test_shifted_parabola);
    TEST(test_constant);
    TEST(test_quartic_value);
    TEST(test_degree6_showcase);
    TEST(test_unbounded_odd);
    TEST(test_unbounded_negative_leading);

    /* 3. Maximize */
    TEST(test_maximize_parabola);
    TEST(test_maximize_quadratic);
    TEST(test_maximize_unbounded);
    TEST(test_min_max_duality);

    /* 4. Multivariate unconstrained */
    TEST(test_bowl_2d);
    TEST(test_shifted_bowl_2d);
    TEST(test_positive_dimensional_solves);
    TEST(test_flat_valley_hyperbola);
    TEST(test_flat_valley_algebraic);
    TEST(test_positive_dim_unbounded);

    /* 5. Constrained */
    TEST(test_disk_linear);
    TEST(test_equality_inequality);
    TEST(test_linear_program);
    TEST(test_infeasible);
    TEST(test_univariate_constrained);

    /* 6. Inexact fallback */
    TEST(test_inexact_fallback);

    /* 7. Declines */
    TEST(test_transcendental_declines);
    TEST(test_unbounded_constrained_declines);
    TEST(test_bad_arity);

    /* 9. Robustness & scope expansion (v0.305 - v0.310) */
    TEST(test_circle_irrational);
    TEST(test_torus_solves);
    TEST(test_disk_inequality_irrational);
    TEST(test_radical_objective);
    TEST(test_fractional_power_declines);
    TEST(test_abs_sum_median);
    TEST(test_abs_interior_min);
    TEST(test_abs_unbounded);
    TEST(test_integer_equality);
    TEST(test_integer_box);
    TEST(test_integer_hard_declines);
    TEST(test_message_routing_check);
    TEST(test_timeconstraint_option);

    /* 9b. Campaign II: 21-40 stress-suite scope + soundness pins */
    TEST(test_separable_decomposition);
    TEST(test_compact_simplex_sextic);
    TEST(test_compact_viviani);
    TEST(test_equality_elimination);
    TEST(test_rational_objective);
    TEST(test_rational_poles_decline);
    TEST(test_integer_parametric);
    TEST(test_mixed_integer_declines);

    /* 8. Memory */
    TEST(test_memory_smoke);

    printf("All Minimize tests passed.\n");
    return 0;
}
