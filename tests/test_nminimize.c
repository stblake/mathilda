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
 *   - Variable locality (Protected, not HoldAll): the search variable does not
 *     leak, and an assigned optimization variable makes the call unevaluated.
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

static void test_integer_domain_alternatives(void) {
    /* Element[x|y, Integers] (Alternatives) declares BOTH x and y integer,
     * exactly like a pair of single declarations. Same LP as
     * test_integer_domain_value -> optimum 1, both variables integer. */
    check_true("Abs[First[NMinimize[{x + y, x + 2 y >= 3, x >= -2, Element[x | y, Integers]}, {x, y}]] - 1.0] < 1.*^-6");
    check_eq("Head[x /. Last[NMinimize[{x + y, x + 2 y >= 3, x >= -2, Element[x | y, Integers]}, {x, y}]]]", "Integer");
    check_eq("Head[y /. Last[NMinimize[{x + y, x + 2 y >= 3, x >= -2, Element[x | y, Integers]}, {x, y}]]]", "Integer");
}

static void test_integer_domain_list(void) {
    /* Element[{x, y}, Integers] (List) is the same multi-variable declaration. */
    check_true("Abs[First[NMinimize[{x + y, x + 2 y >= 3, x >= -2, Element[{x, y}, Integers]}, {x, y}]] - 1.0] < 1.*^-6");
    check_eq("Head[y /. Last[NMinimize[{x + y, x + 2 y >= 3, x >= -2, Element[{x, y}, Integers]}, {x, y}]]]", "Integer");
}

static void test_region_expansion_rescue(void) {
    /* Feasible region x + y >= 80 lies entirely OUTSIDE the default +-10 DE
     * sampling span (there x + y <= 20). Adaptive region expansion must grow
     * the fully-unbounded coordinates until the feasible basin is reached and
     * return the true optimum (50, 40) with value 0, not {Infinity, ...}. */
    check_true("First[NMinimize[{(x-50)^2 + (y-40)^2, x + y >= 80}, {x, y}]] < 1.*^-2");
    check_true("Abs[(x /. Last[NMinimize[{(x-50)^2 + (y-40)^2, x + y >= 80}, {x, y}]]) - 50.0] < 1.*^-1");
    /* The same rescue with an integer coordinate. */
    check_true("First[NMinimize[{(x-50)^2 + (y-40)^2, x + y >= 80, Element[y, Integers]}, {x, y}]] < 1.*^-2");
    /* A genuinely infeasible unbounded problem still returns Infinity (the
     * expansion exhausts without ever finding a feasible point). */
    check_eq("First[NMinimize[{x, x^2 + 1 <= 0}, x]]", "Infinity");
}

static void test_mixed_integer_outside_region(void) {
    /* The continuous optimum (x = 15) lies well outside the +-10 default DE
     * sampling span; y is integer with optimum 3. The mixed-integer polish's
     * continuous-relaxation step must recover it: value 0 at (15, 3). Without
     * it the search would be stranded at the region wall (x ~ 10, value ~ 25). */
    check_true("Abs[First[NMinimize[{(x - 15)^2 + (y - 3)^2, Element[y, Integers]}, {x, y}]]] < 1.*^-2");
    check_true("Abs[(x /. Last[NMinimize[{(x - 15)^2 + (y - 3)^2, Element[y, Integers]}, {x, y}]]) - 15.0] < 1.*^-2");
    check_eq("Head[y /. Last[NMinimize[{(x - 15)^2 + (y - 3)^2, Element[y, Integers]}, {x, y}]]]", "Integer");
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

static void test_postprocess_values(void) {
    /* Full PostProcess value set. On a smooth problem the polish reaches the
     * optimum exactly, so "polish on" values all give ~0 and "polish off"
     * values return the raw (nonzero) SimulatedAnnealing point.
     * Polish ON: True, Automatic, and any named local method (string). */
    check_true("Abs[First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PostProcess\" -> Automatic}]]] < 1.*^-6");
    check_true("Abs[First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PostProcess\" -> \"InteriorPoint\"}]]] < 1.*^-6");
    check_true("Abs[First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PostProcess\" -> \"FindMinimum\"}]]] < 1.*^-6");
    check_true("Abs[First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PostProcess\" -> {\"InteriorPoint\", \"Tolerance\" -> 1.*^-8}}]]] < 1.*^-6");
    /* Polish OFF: None behaves like False (raw point, > 0 on this seed). */
    check_true("First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PostProcess\" -> None}]] > 1.*^-9");
    check_true("First[NMinimize[(x-3)^2 + (y+2)^2, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PostProcess\" -> False}]] > 1.*^-9");
    /* The reported ReflectRatio + string-PostProcess example is accepted. */
    check_eq("Head[NMinimize[{x^2 + y^2 + 0.1 Sin[1000 x] + 0.1 Sin[1000 y], -5 <= x <= 5 && -5 <= y <= 5}, "
             "{x, y}, Method -> {\"NelderMead\", \"ReflectRatio\" -> 1.5, \"PostProcess\" -> \"InteriorPoint\"}]]", "List");
}

static void test_sa_suboptions(void) {
    /* The three "SimulatedAnnealing" sub-options are honored — nm_sa used to
     * do (void)nc and ignore all of them, and the parser silently dropped
     * "PerturbationScale" / "BoltzmannExponent". Deterministic under the fixed
     * default seed. */

    /* "SearchPoints" -> K runs K independent annealing restarts and keeps the
     * global best. On a rugged multimodal surface, 30 restarts reach a far
     * better raw (PostProcess -> False) point than a single chain. */
    check_true("First[NMinimize[{x^2 + y^2 + 10 Sin[3 x]^2 + 10 Sin[3 y]^2, "
               "-5 <= x <= 5 && -5 <= y <= 5}, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"SearchPoints\" -> 30, \"PostProcess\" -> False}]] < "
               "First[NMinimize[{x^2 + y^2 + 10 Sin[3 x]^2 + 10 Sin[3 y]^2, "
               "-5 <= x <= 5 && -5 <= y <= 5}, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"SearchPoints\" -> 1, \"PostProcess\" -> False}]]");

    /* "PerturbationScale" scales the trial-step size. A tiny scale keeps the
     * raw walk near its random start, far from the (3, -2) optimum of this wide
     * box; the default scale explores and reaches it. */
    check_true("First[NMinimize[{(x-3)^2 + (y+2)^2, -50 <= x <= 50 && -50 <= y <= 50}, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PerturbationScale\" -> 0.002, \"PostProcess\" -> False}]] > 10");
    check_true("First[NMinimize[{(x-3)^2 + (y+2)^2, -50 <= x <= 50 && -50 <= y <= 50}, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"PerturbationScale\" -> 1.0, \"PostProcess\" -> False}]] < 1");

    /* "BoltzmannExponent" -> f supplies the acceptance-probability exponent
     * f[i, df, f0]; a proper (negative-for-uphill) exponent still anneals and
     * the polish reaches the optimum. */
    check_true("Abs[First[NMinimize[{(x-3)^2 + (y+2)^2, -50 <= x <= 50 && -50 <= y <= 50}, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"BoltzmannExponent\" -> (-#2/(#1 + 1) &), "
               "\"PostProcess\" -> True}]]] < 1.*^-6");

    /* Invalid sub-option values warn (NMinimize::sopt / ::bexp) and fall back
     * to the defaults rather than failing: the solve still returns a List. */
    check_eq("Head[NMinimize[{(x-3)^2 + (y+2)^2, -50 <= x <= 50 && -50 <= y <= 50}, {x, y}, "
             "Method -> {\"SimulatedAnnealing\", \"PerturbationScale\" -> -3}]]", "List");
    check_eq("Head[NMinimize[{(x-3)^2 + (y+2)^2, -50 <= x <= 50 && -50 <= y <= 50}, {x, y}, "
             "Method -> {\"SimulatedAnnealing\", \"BoltzmannExponent\" -> \"nope\"}]]", "List");
}

static void test_griewank_simulatedannealing(void) {
    /* Griewank-10 on [-600, 600]^10 is strongly multimodal. The default
     * SimulatedAnnealing runs Min[2 n, 50] = 20 independent chains and polishes
     * each chain's best into its basin minimum before ranking (as RandomSearch
     * does per restart), so it reaches a far deeper basin than a single walk can
     * — ~0.015 here, well below Mathematica's ~0.175 on the same problem, and a
     * large improvement on the old single-chain ~0.31. Deterministic under the
     * fixed default seed. */
    check_true("Module[{r = NMinimize[{Sum[x[i]^2/4000, {i, 1, 10}] "
               "- Product[Cos[x[i]/Sqrt[i]], {i, 1, 10}] + 1, "
               "Table[-600 <= x[i] <= 600, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> \"SimulatedAnnealing\"]}, 0 <= First[r] < 0.1]");

    /* The exact reported invocation — all three sub-options together — is
     * accepted and returns a valid, finite feasible minimum. The degenerate
     * "BoltzmannExponent" -> (1/# &) makes acceptance certain (Exp[1/i] > 1),
     * turning each chain into a wide random walk that, with 2x perturbation,
     * settles at a higher local minimum than the default Metropolis schedule;
     * the point of the test is that every sub-option is parsed and applied
     * without error and the result stays finite/feasible. */
    check_true("Module[{r = NMinimize[{Sum[x[i]^2/4000, {i, 1, 10}] "
               "- Product[Cos[x[i]/Sqrt[i]], {i, 1, 10}] + 1, "
               "Table[-600 <= x[i] <= 600, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> {\"SimulatedAnnealing\", \"PerturbationScale\" -> 2.0, "
               "\"BoltzmannExponent\" -> (1/# &), \"SearchPoints\" -> 50}]}, "
               "Head[r] === List && 0 <= First[r] < 200]");
}

static void test_griewank_differentialevolution(void) {
    /* Griewank-10 under explicit "DifferentialEvolution". The final population is
     * spread across basins; polishing the best Min[2 n, 50] distinct members and
     * keeping the deepest local minimum (not the raw global best's basin) reaches
     * the global 0 here — far below Mathematica's ~0.175. Deterministic under the
     * fixed default seed. */
    check_true("Module[{r = NMinimize[{Sum[x[i]^2/4000, {i, 1, 10}] "
               "- Product[Cos[x[i]/Sqrt[i]], {i, 1, 10}] + 1, "
               "Table[-600 <= x[i] <= 600, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> \"DifferentialEvolution\"]}, 0 <= First[r] < 0.1]");

    /* Monotonicity: enlarging "SearchPoints" leaves the population less converged
     * per generation, which used to STRAND the single polished best in a shallow
     * basin and report a worse optimum (SP 100 -> 0.197). With per-member basin
     * polishing a larger population is no worse — it stays well under 0.1. */
    check_true("Module[{r = NMinimize[{Sum[x[i]^2/4000, {i, 1, 10}] "
               "- Product[Cos[x[i]/Sqrt[i]], {i, 1, 10}] + 1, "
               "Table[-600 <= x[i] <= 600, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 100}]}, "
               "0 <= First[r] < 0.1]");
}

static void test_griewank_neldermead(void) {
    /* NelderMead now polishes each restart's converged vertex into its basin
     * minimum before ranking the restarts, rather than ranking raw simplex
     * vertices and polishing only the winner. The default runs Min[2 n, 20]
     * restarts; on Griewank-10 it stays well below Mathematica's ~0.175, and 40
     * restarts reach the global 0. Deterministic under the fixed default seed. */
    check_true("Module[{r = NMinimize[{Sum[x[i]^2/4000, {i, 1, 10}] "
               "- Product[Cos[x[i]/Sqrt[i]], {i, 1, 10}] + 1, "
               "Table[-600 <= x[i] <= 600, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> \"NelderMead\"]}, 0 <= First[r] < 0.1]");
    check_true("Module[{r = NMinimize[{Sum[x[i]^2/4000, {i, 1, 10}] "
               "- Product[Cos[x[i]/Sqrt[i]], {i, 1, 10}] + 1, "
               "Table[-600 <= x[i] <= 600, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> {\"NelderMead\", \"SearchPoints\" -> 40}]}, 0 <= First[r] < 0.01]");
}

static void test_randomsearch_searchpoints_verbatim(void) {
    /* RandomSearch "SearchPoints" was silently capped at 40, making any larger
     * value a no-op. It is now honored verbatim: on a Griewank box narrow enough
     * that random starts can reach the good basins, 200 starts beat 10 — which is
     * only possible if the 200 are actually run. Deterministic (fixed seed). */
    check_true("First[NMinimize[{Sum[x[i]^2/4000, {i, 1, 5}] "
               "- Product[Cos[x[i]/Sqrt[i]], {i, 1, 5}] + 1, "
               "Table[-15 <= x[i] <= 15, {i, 1, 5}]}, Table[x[i], {i, 1, 5}], "
               "Method -> {\"RandomSearch\", \"SearchPoints\" -> 200}]] < "
               "First[NMinimize[{Sum[x[i]^2/4000, {i, 1, 5}] "
               "- Product[Cos[x[i]/Sqrt[i]], {i, 1, 5}] + 1, "
               "Table[-15 <= x[i] <= 15, {i, 1, 5}]}, Table[x[i], {i, 1, 5}], "
               "Method -> {\"RandomSearch\", \"SearchPoints\" -> 10}]]");
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

static void test_penalty_function(void) {
    /* "PenaltyFunction" scores infeasible points during the global search;
     * Automatic is the built-in squared penalty. On a genuinely constrained
     * problem (min x+y on the unit disk -> -Sqrt[2]) every well-formed penalty
     * choice must still converge to the same feasible optimum. */
    check_true("Abs[First[NMinimize[{x + y, x^2 + y^2 <= 1}, {x, y}, "
               "Method -> {\"DifferentialEvolution\", \"PenaltyFunction\" -> Automatic}]] + 1.4142136] < 1.*^-2");
    /* #^2 & is exactly the Automatic penalty. */
    check_true("Abs[First[NMinimize[{x + y, x^2 + y^2 <= 1}, {x, y}, "
               "Method -> {\"DifferentialEvolution\", \"PenaltyFunction\" -> (#^2 &)}]] + 1.4142136] < 1.*^-2");
    /* A different (still monotone nonnegative) penalty also converges. */
    check_true("Abs[First[NMinimize[{x + y, x^2 + y^2 <= 1}, {x, y}, "
               "Method -> {\"DifferentialEvolution\", \"PenaltyFunction\" -> (10 # &)}]] + 1.4142136] < 1.*^-2");
    check_true("Abs[First[NMinimize[{x + y, x^2 + y^2 <= 1}, {x, y}, "
               "Method -> {\"NelderMead\", \"PenaltyFunction\" -> Sqrt}]] + 1.4142136] < 1.*^-2");
    /* An invalid value (a number) is rejected with NMinimize::penf and the run
     * falls back to Automatic, still returning the correct shape and answer. */
    check_true("Abs[First[NMinimize[{x + y, x^2 + y^2 <= 1}, {x, y}, "
               "Method -> {\"DifferentialEvolution\", \"PenaltyFunction\" -> 5}]] + 1.4142136] < 1.*^-2");
    /* On a box-only problem PenaltyFunction is inert (no general constraints);
     * the user's scaled-quadratic example still reaches ~0. */
    check_true("First[NMinimize[{Sum[10^(6 (i-1)/9) x[i]^2, {i, 1, 10}], "
               "Table[-10 <= x[i] <= 10, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> {\"NelderMead\", \"Tolerance\" -> 10^-8, \"PenaltyFunction\" -> Automatic}]] < 1.*^-6");
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

static void test_symbol_indirection(void) {
    /* The problem and the variable spec are held (HoldAll). A problem or a
     * variable list passed via a bound symbol — prob = {f, cons}; NMinimize[prob,
     * vars] — must resolve to its value, not be treated as a single opaque
     * variable. Regression: the {f, cons} list-via-symbol returned the infeasible
     * sentinel, and a vars-list-via-symbol produced one rule keyed by the whole
     * list. Module scopes the bound symbols so no global state leaks. */
    check_true("Module[{p = {x^2 + y^2, x + y >= 1}}, "
               "Abs[First[NMinimize[p, {x, y}]] - 0.5] < 1.*^-2]");           /* objective via symbol */
    check_true("Module[{v = {x, y}}, "
               "Abs[First[NMinimize[{x^2 + y^2, x + y >= 1}, v]] - 0.5] < 1.*^-2]"); /* vars via symbol */
    check_true("Module[{p = {x^2 + y^2, x + y >= 1}, v = {x, y}}, "
               "Length[Last[NMinimize[p, v]]] == 2]");                        /* both, correct arity */
    check_true("Module[{p = {Sum[x[i]^2, {i, 1, 3}], Table[-5 <= x[i] <= 5, {i, 1, 3}]}, "
               "v = Table[x[i], {i, 1, 3}]}, Abs[First[NMinimize[p, v]]] < 1.*^-2]"); /* indexed via symbol */
    /* An unbound bare symbol is still a single optimization variable. */
    check_true("Abs[First[NMinimize[x^4 - 3 x^2 - x, x]] - (-3.5139097)] < 1.*^-4");
}

static void test_search_points_honored(void) {
    /* Regression: an explicit DifferentialEvolution "SearchPoints" must be
     * honored, not silently capped at the automatic ceiling (was 40). Two runs
     * that differ only in SearchPoints, both above the old cap, once returned
     * byte-identical results because the cap swallowed the value — the RNG
     * trajectory is population-size dependent, so on a multimodal surface a
     * genuinely different NP lands on a different basin. Deterministic under the
     * fixed default seed. */
    check_true("First[NMinimize[{Sum[x[i]^2 - 10 Cos[2 Pi x[i]], {i, 1, 8}], "
               "Table[-5.12 <= x[i] <= 5.12, {i, 1, 8}]}, Table[x[i], {i, 1, 8}], "
               "Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 250}]] != "
               "First[NMinimize[{Sum[x[i]^2 - 10 Cos[2 Pi x[i]], {i, 1, 8}], "
               "Table[-5.12 <= x[i] <= 5.12, {i, 1, 8}]}, Table[x[i], {i, 1, 8}], "
               "Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 60}]]");
    /* A large explicit population still solves a convex bowl to the optimum. */
    check_true("Abs[First[NMinimize[{Sum[x[i]^2, {i, 1, 3}], Table[-5 <= x[i] <= 5, {i, 1, 3}]}, "
               "Table[x[i], {i, 1, 3}], Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 120}]]] < 1.*^-3");
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
/* 10. Variable locality (Protected, not HoldAll)                      */
/* ------------------------------------------------------------------ */

static void test_variable_locality(void) {
    /* Attributes[NMinimize] == {Protected} (matching Mathematica), NOT HoldAll.
     * The Block-style snapshot/restore of a search variable's OwnValue during the
     * numeric search must not leak: an unbound variable stays unbound. */
    check_eq("(ClearAll[x]; NMinimize[{x^2, -5 <= x <= 5}, x]; ValueQ[x])", "False");
    /* Because arguments are evaluated (not held), an optimization variable that
     * already carries a value resolves before the call, the variable spec is then
     * malformed, and NMinimize returns unevaluated — exactly as in Mathematica —
     * without mutating the pre-existing value. */
    check_eq("(x = 5; Head[NMinimize[x^2, x]])", "NMinimize");
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

static void test_indexed_variable_locality(void) {
    /* Indexed analogue. An unbound indexed variable stays unbound after the
     * search (indexed vars are rewritten to fresh synthetic symbols, so x itself
     * is never bound). And a pre-existing indexed value makes the now-evaluated
     * variable spec malformed, so the call returns unevaluated and leaves the
     * value intact — matching Mathematica's non-HoldAll NMinimize. */
    check_eq("(ClearAll[x]; NMinimize[Sum[(x[i]-i)^2, {i,1,3}], Table[x[i], {i,1,3}]]; ValueQ[x[1]])", "False");
    check_eq("(x[1] = 42; Head[NMinimize[Sum[(x[i] - i)^2, {i, 1, 3}], Table[x[i], {i, 1, 3}]]])", "NMinimize");
    check_eq("(x[1] = 42; NMinimize[Sum[(x[i] - i)^2, {i, 1, 3}], Table[x[i], {i, 1, 3}]]; x[1])", "42");
    check_eq("(x[1] =.; ValueQ[x[1]])", "False");
}

static void test_indexed_real_coefficient(void) {
    /* Regression: a Real coefficient multiplying an indexed term used to
     * numericalize the index itself (x[i] -> x[i.]) when the held objective was
     * expanded, leaving an unbound x[1.] in the compiled body. Every trial
     * point then scored non-numeric and NMinimize returned the infeasible
     * 1e+300 sentinel. The fix stops inexact contagion threading N into a
     * non-numeric head's arguments (see test_inexact_contagion in
     * test_numeric.c). 3-D Ackley (global min 0 at the origin, inside the box)
     * must now solve to ~0; deterministic under the fixed seed (~9.5e-10). */
    check_true("First[NMinimize[{-20 Exp[-0.2 Sqrt[1/3 Sum[x[i]^2, {i, 1, 3}]]] "
               "- Exp[1/3 Sum[Cos[2 Pi x[i]], {i, 1, 3}]] + 20 + E, "
               "Table[-32 <= x[i] <= 32, {i, 1, 3}]}, Table[x[i], {i, 1, 3}]]] < 0.01");
    /* Minimal reproducer: real coefficient, indexed vars, no Sum at all. */
    check_true("First[NMinimize[{Sqrt[0.333 (x[1]^2 + x[2]^2 + x[3]^2)], "
               "Table[-5 <= x[i] <= 5, {i, 1, 3}]}, Table[x[i], {i, 1, 3}]]] < 0.001");
}

static void test_de_boundary_no_stagnation(void) {
    /* Regression: DifferentialEvolution used to CLAMP an out-of-range mutant to
     * the box boundary. On the 10-D Schwefel function (global min ~0 at
     * x_i = 420.9687, well inside [-500, 500]) that stranded the search: once
     * several members shared the exact boundary value for a coordinate, their
     * mutation differentials for it collapsed to zero and it froze on the wall.
     * The reported run returned 1758.88 with five coordinates pinned to +/-500,
     * worse than Mathematica's 1221.25. Bounce-back reinitialisation (nm_de)
     * fixed it: the default-seed objective now solves to ~963, and -- the
     * platform-independent invariant this test asserts -- NO coordinate sits on
     * the boundary. The objective threshold is loose (buggy 1758 vs fixed ~963)
     * to tolerate cross-platform floating-point drift in the search path. */
    check_true("First[NMinimize[{418.9829*10 - Sum[x[i] Sin[Sqrt[Abs[x[i]]]], {i, 1, 10}], "
               "Table[-500 <= x[i] <= 500, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 100, "
               "\"PostProcess\" -> False}]] < 1400");
    /* No coordinate parked on the +/-500 boundary (the stagnation signature). */
    check_true("Count[Table[x[i], {i, 1, 10}] /. Last[NMinimize[{418.9829*10 "
               "- Sum[x[i] Sin[Sqrt[Abs[x[i]]]], {i, 1, 10}], "
               "Table[-500 <= x[i] <= 500, {i, 1, 10}]}, Table[x[i], {i, 1, 10}], "
               "Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 100, "
               "\"PostProcess\" -> False}]], v_ /; Abs[Abs[v] - 500] < 0.001] == 0");
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
    TEST(test_integer_domain_alternatives);
    TEST(test_integer_domain_list);
    TEST(test_region_expansion_rescue);
    TEST(test_mixed_integer_outside_region);

    /* 6. Infeasible */
    TEST(test_infeasible);

    /* 7. Methods */
    TEST(test_method_de);
    TEST(test_method_neldermead);
    TEST(test_method_randomsearch);
    TEST(test_method_simulatedannealing);
    TEST(test_method_suboptions);
    TEST(test_search_points_honored);
    TEST(test_symbol_indirection);
    TEST(test_neldermead_suboptions);
    TEST(test_postprocess_values);
    TEST(test_sa_suboptions);
    TEST(test_griewank_simulatedannealing);
    TEST(test_griewank_differentialevolution);
    TEST(test_griewank_neldermead);
    TEST(test_randomsearch_searchpoints_verbatim);
    TEST(test_neldermead_shrink_tolerance);
    TEST(test_bukin6_no_warning);
    TEST(test_initial_points);
    TEST(test_penalty_function);
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
    TEST(test_variable_locality);

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
    TEST(test_indexed_variable_locality);
    TEST(test_indexed_real_coefficient);
    TEST(test_de_boundary_no_stagnation);

    printf("All NMinimize tests passed.\n");
    return 0;
}
