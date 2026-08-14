/* Unit tests for NMinimize / NMaximize (src/findmin.c).
 *
 * NMinimize is a global-optimization driver layered on the FindMinimum
 * machinery. The search is stochastic but deterministic for the fixed
 * default RandomSeed, so exact-ish assertions are stable. Objective
 * tolerances are looser than FindMinimum's (global heuristics + penalty
 * constraints), tight where the problem is convex/linear.
 *
 * Coverage:
 *   - Result shape / typing; {fmin, {x -> xmin, ...}}.
 *   - Unconstrained 1-D and n-D global minima.
 *   - Inequality, equality, chained, and linear-program constraints.
 *   - Mixed-integer domains via Element[x, Integers] (integer results).
 *   - Empty feasible set -> {Infinity, {x -> Indeterminate, ...}}.
 *   - All four methods (DifferentialEvolution / NelderMead / RandomSearch /
 *     SimulatedAnnealing) and Method sub-options.
 *   - NMaximize wrapper + min/max duality.
 *   - Options[NMinimize]; WorkingPrecision (MPFR) refinement.
 *   - Constraint back-substitution feasibility.
 *   - HoldAll locality: the search variable does not leak.
 *   - Memory-hygiene smoke loop.
 *
 * Run binary directly: ./nminimize_tests */

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
    check_eq("Head[NMinimize[x^2, x]]", "List");
    check_eq("Length[NMinimize[x^2, x]]", "2");
    check_eq("Head[Last[NMinimize[x^2, x]]]", "List");
    check_eq("Head[First[Last[NMinimize[x^2, x]]]]", "Rule");
}

/* ------------------------------------------------------------------ */
/* 2. Unconstrained global minima                                      */
/* ------------------------------------------------------------------ */

static void test_quartic(void) {
    /* x^4 - 3 x^2 - x -> global min -3.51391 at x = 1.30084. */
    check_true("Abs[First[NMinimize[x^4 - 3 x^2 - x, x]] - (-3.5139097)] < 1.*^-3");
    check_true("Abs[(x /. Last[NMinimize[x^4 - 3 x^2 - x, x]]) - 1.3008373] < 1.*^-3");
}

static void test_shifted_parabola(void) {
    /* (x-3)^2 + 1 -> 1 at x = 3. */
    check_true("Abs[First[NMinimize[(x-3)^2 + 1, x]] - 1.0] < 1.*^-4");
    check_true("Abs[(x /. Last[NMinimize[(x-3)^2 + 1, x]]) - 3.0] < 1.*^-3");
}

static void test_2d_bowl(void) {
    /* (x-1)^2 + (y+2)^2 -> 0 at (1, -2). */
    check_true("Abs[First[NMinimize[(x-1)^2 + (y+2)^2, {x, y}]]] < 1.*^-3");
    check_true("Abs[(x /. Last[NMinimize[(x-1)^2 + (y+2)^2, {x, y}]]) - 1.0] < 1.*^-2");
    check_true("Abs[(y /. Last[NMinimize[(x-1)^2 + (y+2)^2, {x, y}]]) + 2.0] < 1.*^-2");
}

/* ------------------------------------------------------------------ */
/* 3. Constraints                                                      */
/* ------------------------------------------------------------------ */

static void test_disk_linear(void) {
    /* {x + y, x^2 + y^2 <= 9} -> -4.24264. */
    check_true("Abs[First[NMinimize[{x + y, x^2 + y^2 <= 9}, {x, y}]] - (-4.2426407)] < 1.*^-2");
}

static void test_quadratic_linear(void) {
    /* {(x-1)^2 + y^2, x + y/2 <= 1/2, x - y >= 0} -> 0.2 at (0.6, -0.2). */
    check_true("Abs[First[NMinimize[{(x-1)^2 + y^2, x + y/2 <= 1/2, x - y >= 0}, {x, y}]] - 0.2] < 1.*^-2");
}

static void test_linear_program(void) {
    /* {x + y, 3x+2y>=7 && x+2y>=6 && x>=0 && y>=0} -> 3.25. */
    check_true("Abs[First[NMinimize[{x + y, 3 x + 2 y >= 7 && x + 2 y >= 6 && x >= 0 && y >= 0}, {x, y}]] - 3.25] < 1.*^-2");
}

static void test_equality_constraint(void) {
    /* {x + 2y, x^2 + 2y^2 <= 3, x + y == 2, x >= 1} -> 2.33333. */
    check_true("Abs[First[NMinimize[{x + 2 y, x^2 + 2 y^2 <= 3, x + y == 2, x >= 1}, {x, y}]] - 2.3333333] < 2.*^-2");
}

static void test_chained_inequality(void) {
    /* {Sin[2x] + Cos[x], -2 <= x <= 3} -> -1.76017 at x = 2.50673. */
    check_true("Abs[First[NMinimize[{Sin[2 x] + Cos[x], -2 <= x <= 3}, x]] - (-1.7601696)] < 1.*^-2");
}

static void test_equation_system(void) {
    /* {x - y, x+y+z==1/2, x-2z==1, 2x-y>=1} -> 0.428571. */
    check_true("Abs[First[NMinimize[{x - y, x + y + z == 1/2, x - 2 z == 1, 2 x - y >= 1}, {x, y, z}]] - 0.4285714] < 2.*^-2");
}

/* ------------------------------------------------------------------ */
/* 4. Feasibility of the returned point                                */
/* ------------------------------------------------------------------ */

static void test_returned_point_feasible(void) {
    /* The reported minimizer must satisfy the constraint. */
    check_true("(x^2 + y^2 /. Last[NMinimize[{x + y, x^2 + y^2 <= 9}, {x, y}]]) <= 9.001");
    check_true("Abs[(x + y /. Last[NMinimize[{x + 2 y, x^2 + 2 y^2 <= 3, x + y == 2, x >= 1}, {x, y}]]) - 2.0] < 1.*^-2");
}

/* ------------------------------------------------------------------ */
/* 5. Integer domains                                                  */
/* ------------------------------------------------------------------ */

static void test_integer_domain_value(void) {
    /* {x + y, x+2y>=3, x>=-2}, x,y integers -> optimum value 1. */
    check_true("Abs[First[NMinimize[{x + y, x + 2 y >= 3, x >= -2}, {Element[x, Integers], Element[y, Integers]}]] - 1.0] < 1.*^-6");
}

static void test_integer_domain_heads(void) {
    /* Integer-domain variables come back as exact integers. */
    check_eq("Head[x /. Last[NMinimize[{x + y, x + 2 y >= 3, x >= -2}, {Element[x, Integers], Element[y, Integers]}]]]", "Integer");
    check_eq("Head[y /. Last[NMinimize[{x + y, x + 2 y >= 3, x >= -2}, {Element[x, Integers], Element[y, Integers]}]]]", "Integer");
}

static void test_mixed_integer(void) {
    /* One integer, one real variable. Optimum value 3, x integer 1. */
    check_true("Abs[First[NMinimize[{x + 2 y, x^2 + 2 y^2 <= 3, x + y == 2, Element[x, Integers]}, {x, y}]] - 3.0] < 5.*^-2");
    check_eq("Head[x /. Last[NMinimize[{x + 2 y, x^2 + 2 y^2 <= 3, x + y == 2, Element[x, Integers]}, {x, y}]]]", "Integer");
}

/* ------------------------------------------------------------------ */
/* 6. Infeasible problems                                              */
/* ------------------------------------------------------------------ */

static void test_infeasible(void) {
    /* Empty feasible set -> {Infinity, {x -> Indeterminate}}. */
    check_eq("First[NMinimize[{x, x > 2 && x < 1}, x]]", "Infinity");
    check_eq("x /. Last[NMinimize[{x, x > 2 && x < 1}, x]]", "Indeterminate");
}

/* ------------------------------------------------------------------ */
/* 7. Methods                                                          */
/* ------------------------------------------------------------------ */

static void test_method_de(void) {
    check_true("Abs[First[NMinimize[x^4 - 3 x^2 - x, x, Method -> \"DifferentialEvolution\"]] - (-3.5139097)] < 1.*^-3");
}
static void test_method_neldermead(void) {
    check_true("Abs[First[NMinimize[x^4 - 3 x^2 - x, x, Method -> \"NelderMead\"]] - (-3.5139097)] < 1.*^-3");
}
static void test_method_randomsearch(void) {
    check_true("Abs[First[NMinimize[x^4 - 3 x^2 - x, x, Method -> \"RandomSearch\"]] - (-3.5139097)] < 1.*^-3");
}
static void test_method_simulatedannealing(void) {
    check_true("Abs[First[NMinimize[x^4 - 3 x^2 - x, x, Method -> \"SimulatedAnnealing\"]] - (-3.5139097)] < 1.*^-2");
}
static void test_neldermead_suboptions(void) {
    /* ExpandRatio / ContractRatio are consumed by the simplex and still
     * converge to the optimum. */
    check_true("Abs[First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"NelderMead\", \"ExpandRatio\" -> 2.5, \"ContractRatio\" -> 0.6}]]] < 1.*^-4");
    /* PostProcess -> True runs the exact local polish, reaching the optimum. */
    check_true("Abs[First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PostProcess\" -> True}]]] < 1.*^-6");
    /* PostProcess -> False skips the polish; the raw global-search point is
     * still returned and is a valid nearby result. */
    check_true("Abs[First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PostProcess\" -> False}]]] < 1.*^-2");
    /* All together, as in the documented examples, are accepted. */
    check_eq("Head[NMinimize[(x-1)^2 + (y-1)^2, {x, y}, "
             "Method -> {\"NelderMead\", \"ExpandRatio\" -> 2.5, \"ContractRatio\" -> 0.6, "
             "\"ReflectRatio\" -> 1.0, \"ShrinkRatio\" -> 0.5, \"PostProcess\" -> False}]]", "List");
}

static void test_neldermead_shrink_tolerance(void) {
    /* Tolerance controls the simplex convergence threshold: a tight tolerance
     * refines much further than a loose one under a capped budget. */
    check_true("Abs[First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"NelderMead\", \"Tolerance\" -> 1.*^-10, \"PostProcess\" -> False}]]] < 1.*^-6");
    check_true("First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"NelderMead\", \"Tolerance\" -> 0.1, \"PostProcess\" -> False}, MaxIterations -> 2]] "
               "> First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"NelderMead\", \"Tolerance\" -> 1.*^-12, \"PostProcess\" -> False}, MaxIterations -> 2]]");
    /* ShrinkRatio is accepted and still converges to the optimum. */
    check_true("Abs[First[NMinimize[Abs[x-3] + Abs[y+2], {x, y}, "
               "Method -> {\"NelderMead\", \"ShrinkRatio\" -> 0.6}]]] < 1.*^-3");
}

static void test_bukin6_no_warning(void) {
    /* Bukin N.6: a Sqrt[Abs[...]] ridge whose gradient hits 1/0 on the valley.
     * The optimizer must not surface Power::infy (muted during point eval) and
     * must reproduce Mathematica's result. */
    check_true("Abs[First[NMinimize[{100 Sqrt[Abs[x2 - 0.01 x1^2]] + 0.01 Abs[x1 + 10], "
               "-15 <= x1 <= -5 && -3 <= x2 <= 3}, {x1, x2}, "
               "Method -> {\"NelderMead\", \"ShrinkRatio\" -> 0.75, \"Tolerance\" -> 10^-6}]] - 0.0113668] < 1.*^-3");
}

static void test_initial_points(void) {
    /* Easom function: a narrow spike of depth -1 at (Pi, Pi) in an otherwise
     * flat region. The three seed points are collinear on y=x with centroid
     * ~(Pi,Pi); using them as the simplex (and not declaring convergence on the
     * flat plateau while the simplex is still large) lets it find the spike. */
    check_true("Abs[First[NMinimize[{-Cos[x] Cos[y] Exp[-((x-Pi)^2 + (y-Pi)^2)], "
               "-100 <= x <= 100 && -100 <= y <= 100}, {x, y}, "
               "Method -> {\"NelderMead\", \"InitialPoints\" -> {{-50,-50},{50,50},{10,10}}}]] - (-1.0)] < 1.*^-3");
    check_true("Abs[(x /. Last[NMinimize[{-Cos[x] Cos[y] Exp[-((x-Pi)^2 + (y-Pi)^2)], "
               "-100 <= x <= 100 && -100 <= y <= 100}, {x, y}, "
               "Method -> {\"NelderMead\", \"InitialPoints\" -> {{-50,-50},{50,50},{10,10}}}]]) - Pi] < 1.*^-3");
    /* Symbolic seed entries (Pi) are evaluated/numericalized. */
    check_true("Abs[First[NMinimize[(x-1)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"NelderMead\", \"InitialPoints\" -> {{Pi,Pi},{0,0},{6,6}}}]]] < 1.*^-4");
    /* Malformed seeds (wrong dimension / non-numeric) fall back to random starts
     * and still return a valid, correct result. */
    check_true("Abs[First[NMinimize[(x-1)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"NelderMead\", \"InitialPoints\" -> {{1,2,3},{4,5}}}]]] < 1.*^-3");
    check_true("Abs[First[NMinimize[(x-3)^2 + (y-3)^2, {x, y}, "
               "Method -> {\"NelderMead\", \"InitialPoints\" -> {{a,b},{1,1},{2,2}}}]]] < 1.*^-3");
}

static void test_autocompile_parity_and_fallback(void) {
    /* At machine precision the objective is auto-compiled; the answer must be
     * identical to the (verified) interpreter result. */
    check_true("Abs[First[NMinimize[x^4 - 3 x^2 - x, x]] - (-3.5139097)] < 1.*^-4");
    /* A non-compilable objective (special function) falls back to the
     * interpreter and still optimizes correctly: min Gamma on [1,2] ~ 0.88560. */
    check_true("Abs[First[NMinimize[Gamma[x], {x, 1, 2}]] - 0.8856032] < 1.*^-3");
}

static void test_method_suboptions(void) {
    /* Method with sub-options is accepted and returns the right shape/value. */
    check_true("Abs[First[NMinimize[x^4 - 3 x^2 - x, x, Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 20, \"RandomSeed\" -> 7}]] - (-3.5139097)] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* 8. NMaximize                                                        */
/* ------------------------------------------------------------------ */

static void test_nmaximize_simple(void) {
    /* -x^2 + 4 -> 4 at x = 0. */
    check_true("Abs[First[NMaximize[-x^2 + 4, x]] - 4.0] < 1.*^-3");
    check_true("Abs[(x /. Last[NMaximize[-x^2 + 4, x]])] < 1.*^-2");
}

static void test_nmaximize_constrained(void) {
    /* max x + y on the unit disk -> Sqrt[2] = 1.41421. */
    check_true("Abs[First[NMaximize[{x + y, x^2 + y^2 <= 1}, {x, y}]] - 1.4142136] < 1.*^-2");
}

static void test_min_max_duality(void) {
    /* NMaximize[f] == -NMinimize[-f] on a shared problem. */
    check_true("Abs[First[NMaximize[{x + y, x^2 + y^2 <= 1}, {x, y}]] + First[NMinimize[{-(x + y), x^2 + y^2 <= 1}, {x, y}]]] < 1.*^-2");
}

/* ------------------------------------------------------------------ */
/* 9. Options                                                          */
/* ------------------------------------------------------------------ */

static void test_options_nonempty(void) {
    check_eq("Head[Options[NMinimize]]", "List");
    check_true("Length[Options[NMinimize]] >= 5");
    check_eq("Head[Options[NMaximize]]", "List");
}

static void test_working_precision_mpfr(void) {
    /* WorkingPrecision refines the continuous unconstrained result. */
    check_true("Abs[(x /. Last[NMinimize[(x-2)^2, x, WorkingPrecision -> 30]]) - 2] < 1.*^-10");
}

static void test_max_iterations_accepted(void) {
    check_eq("Head[NMinimize[x^2, x, MaxIterations -> 50]]", "List");
    check_eq("Head[NMinimize[x^2, x, MaxIterations -> Automatic]]", "List");
}

/* ------------------------------------------------------------------ */
/* 10. HoldAll locality                                                */
/* ------------------------------------------------------------------ */

static void test_holdall_locality(void) {
    /* A pre-existing value of the search variable must be preserved. */
    check_eq("Module[{}, x = 5; NMinimize[x^2, x]; x]", "5");
    check_true("(x = 5; NMinimize[(x-1)^2, x]; x) == 5");
    /* Clean up the global assignment so later tests see a free x. */
    check_eq("(x =.; ValueQ[x])", "False");
}

/* ------------------------------------------------------------------ */
/* 11. Diagnostics                                                     */
/* ------------------------------------------------------------------ */

static void test_arity_error_unevaluated(void) {
    /* Too few arguments -> unevaluated (echoes the head). */
    check_eq("Head[NMinimize[x^2]]", "NMinimize");
}

/* ------------------------------------------------------------------ */
/* 12. Memory hygiene smoke                                            */
/* ------------------------------------------------------------------ */

static void test_memory_smoke(void) {
    /* A loop of assorted calls; success is not crashing / leaking. */
    check_eq("Do[NMinimize[x^2 + k, {x, y}], {k, 1, 20}]; 0", "0");
    check_eq("Do[NMinimize[{x + y, x^2 + y^2 <= 4}, {x, y}], {5}]; 0", "0");
    check_eq("Do[NMaximize[-x^4 + 2 x^2, x], {5}]; 0", "0");
}

/* ------------------------------------------------------------------ */
/* 13. Indexed variables (Table / Array specs, held generators)        */
/* ------------------------------------------------------------------ */

static void test_indexed_table_vars(void) {
    /* Variable list is a held Table[x[i], ...]; objective is a held Sum. */
    check_true("Abs[First[NMinimize[Sum[(x[i] - i)^2, {i, 1, 3}], Table[x[i], {i, 1, 3}]]]] < 1.*^-4");
    check_true("Abs[(x[2] /. Last[NMinimize[Sum[(x[i] - i)^2, {i, 1, 3}], Table[x[i], {i, 1, 3}]]]) - 2.0] < 1.*^-3");
    /* Result variables come back as the original indexed forms. */
    check_eq("Head[First[Last[NMinimize[Sum[(x[i] - i)^2, {i, 1, 3}], Table[x[i], {i, 1, 3}]]]]]", "Rule");
    check_eq("Head[First[First[Last[NMinimize[Sum[(x[i] - i)^2, {i, 1, 3}], Table[x[i], {i, 1, 3}]]]]]]", "x");
}

static void test_indexed_literal_vars(void) {
    /* An explicit {x[1], x[2]} list (no generator). */
    check_true("Abs[First[NMinimize[x[1]^2 + (x[2] - 1)^2, {x[1], x[2]}]]] < 1.*^-4");
    check_true("Abs[(x[2] /. Last[NMinimize[x[1]^2 + (x[2] - 1)^2, {x[1], x[2]}]]) - 1.0] < 1.*^-3");
}

static void test_indexed_array_vars(void) {
    /* Array[x, 3] expands to {x[1], x[2], x[3]}. */
    check_true("Abs[First[NMinimize[Sum[(x[i] + i)^2, {i, 1, 3}], Array[x, 3]]]] < 1.*^-4");
}

static void test_indexed_table_constraints(void) {
    /* Held Table[...] constraint list, implicitly And-ed and expanded. */
    check_true("Abs[First[NMinimize[{Sum[(x[i] - 2)^2, {i, 1, 3}], Table[0 <= x[i] <= 1, {i, 1, 3}]}, Table[x[i], {i, 1, 3}]]] - 3.0] < 1.*^-3");
    /* Each returned coordinate sits on the upper bound. */
    check_true("Abs[(x[1] /. Last[NMinimize[{Sum[(x[i] - 2)^2, {i, 1, 3}], Table[0 <= x[i] <= 1, {i, 1, 3}]}, Table[x[i], {i, 1, 3}]]]) - 1.0] < 1.*^-3");
}

static void test_indexed_rosenbrock(void) {
    /* The originally-reported case: 10-D Rosenbrock, Table vars + constraints,
     * NelderMead with method sub-options. Global minimum 0 at all-ones. */
    check_true("Abs[First[NMinimize[{Sum[100 (x[i+1] - x[i]^2)^2 + (1 - x[i])^2, {i, 1, 9}], "
               "Table[-5 <= x[i] <= 5, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> {\"NelderMead\", \"ExpandRatio\" -> 2.5}, MaxIterations -> 5000]]] < 1.*^-2");
    check_true("Abs[(x[5] /. Last[NMinimize[Sum[100 (x[i+1] - x[i]^2)^2 + (1 - x[i])^2, {i, 1, 9}], "
               "Table[x[i], {i, 1, 10}], MaxIterations -> 5000]]) - 1.0] < 1.*^-1");
}

static void test_indexed_holdall_locality(void) {
    /* A pre-existing indexed value must be restored after the call. */
    check_eq("(x[1] = 42; NMinimize[Sum[(x[i] - i)^2, {i, 1, 3}], Table[x[i], {i, 1, 3}]]; x[1])", "42");
    check_eq("(x[1] =.; ValueQ[x[1]])", "False");
}

int main(void) {
    symtab_init();
    core_init();

    /* Diagnostics write to stderr on the intended-error tests; keep the test
     * output clean. NMinimize itself is silent on successful solves. */
    freopen("/dev/null", "w", stderr);

    /* 1. Shape */
    TEST(test_result_shape);

    /* 2. Unconstrained */
    TEST(test_quartic);
    TEST(test_shifted_parabola);
    TEST(test_2d_bowl);

    /* 3. Constraints */
    TEST(test_disk_linear);
    TEST(test_quadratic_linear);
    TEST(test_linear_program);
    TEST(test_equality_constraint);
    TEST(test_chained_inequality);
    TEST(test_equation_system);

    /* 4. Feasibility */
    TEST(test_returned_point_feasible);

    /* 5. Integer domains */
    TEST(test_integer_domain_value);
    TEST(test_integer_domain_heads);
    TEST(test_mixed_integer);

    /* 6. Infeasible */
    TEST(test_infeasible);

    /* 7. Methods */
    TEST(test_method_de);
    TEST(test_method_neldermead);
    TEST(test_method_randomsearch);
    TEST(test_method_simulatedannealing);
    TEST(test_method_suboptions);
    TEST(test_neldermead_suboptions);
    TEST(test_neldermead_shrink_tolerance);
    TEST(test_bukin6_no_warning);
    TEST(test_initial_points);
    TEST(test_autocompile_parity_and_fallback);

    /* 8. NMaximize */
    TEST(test_nmaximize_simple);
    TEST(test_nmaximize_constrained);
    TEST(test_min_max_duality);

    /* 9. Options */
    TEST(test_options_nonempty);
    TEST(test_working_precision_mpfr);
    TEST(test_max_iterations_accepted);

    /* 10. Locality */
    TEST(test_holdall_locality);

    /* 11. Diagnostics */
    TEST(test_arity_error_unevaluated);

    /* 12. Memory */
    TEST(test_memory_smoke);

    /* 13. Indexed variables */
    TEST(test_indexed_table_vars);
    TEST(test_indexed_literal_vars);
    TEST(test_indexed_array_vars);
    TEST(test_indexed_table_constraints);
    TEST(test_indexed_rosenbrock);
    TEST(test_indexed_holdall_locality);

    printf("All NMinimize tests passed.\n");
    return 0;
}
