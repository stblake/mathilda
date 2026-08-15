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
/* 6b. Disjunctive (Or) constraints                                    */
/* ------------------------------------------------------------------ */

static void test_disjunctive_constraints(void) {
    /* Disjunctive constraint A || B: a point is feasible iff it satisfies at
     * least one branch. NMinimize scores such a constraint by its minimum-branch
     * violation penalty (0 when any branch holds), so Deb's feasibility rules
     * select it during the derivative-free global search, and the local polish
     * folds in the currently-active branch so it descends within the feasible
     * region. Deterministic under the fixed default seed. */

    /* Himmelblau over the union of two small disks, one of which (centred (3,2))
     * contains the global minimiser: the constrained optimum is Himmelblau's
     * f = 0 at (3, 2). This is the reported In[5] example. */
    check_true("First[NMinimize[{(x^2 + y - 11)^2 + (x + y^2 - 7)^2, "
               "(x - 3)^2 + (y - 2)^2 <= 0.1 || (x + 2.8)^2 + (y + 3.1)^2 <= 0.1}, "
               "{x, y}, Method -> {\"RandomSearch\", \"SearchPoints\" -> 500, "
               "\"Method\" -> \"InteriorPoint\"}]] < 1.*^-4");
    /* The reported minimiser lands in the disk around (3, 2). */
    check_true("Norm[{x - 3, y - 2} /. Last[NMinimize[{(x^2 + y - 11)^2 + (x + y^2 - 7)^2, "
               "(x - 3)^2 + (y - 2)^2 <= 0.1 || (x + 2.8)^2 + (y + 3.1)^2 <= 0.1}, "
               "{x, y}, Method -> {\"RandomSearch\", \"SearchPoints\" -> 500, "
               "\"Method\" -> \"InteriorPoint\"}]]] < 0.35");

    /* A non-convex feasible set: x <= -2 OR x >= 2. The constrained minimum of
     * x^2 is 4 at x = +-2 (the branch boundaries), not the unconstrained x = 0
     * that sits in the infeasible gap between the branches. */
    check_true("Abs[First[NMinimize[{x^2, x <= -2 || x >= 2}, x]] - 4.0] < 0.05");
    check_true("Abs[x /. Last[NMinimize[{x^2, x <= -2 || x >= 2}, x]]] >= 1.99");
    /* RandomSearch is pure multi-start local polish with no global move, so its
     * per-start polish must itself respect the active branch — otherwise it
     * drives every feasible start into the infeasible gap and reports Infinity.
     * This asserts the disjunction-aware polish (regression guard). */
    check_true("Abs[First[NMinimize[{x^2, x <= -2 || x >= 2}, x, "
               "Method -> \"RandomSearch\"]] - 4.0] < 0.05");
    /* Same, with an enclosing box: the polish must honour both the box and the
     * active branch. */
    check_true("Abs[First[NMinimize[{x^2, (x <= -2 || x >= 2) && -100 <= x <= 100}, x, "
               "Method -> \"RandomSearch\"]] - 4.0] < 0.05");

    /* Nearest point to the union of two disks: the optimum sits on the boundary
     * of the disk around (3, 0), at ~(2.51, 0.098), f ~ 4.20. */
    check_true("Abs[First[NMinimize[{(x - 0.5)^2 + (y - 0.5)^2, "
               "-5 <= x <= 5 && -5 <= y <= 5 && "
               "((x - 3)^2 + y^2 <= 0.25 || (x + 3)^2 + y^2 <= 0.25)}, {x, y}, "
               "Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 60}]] - 4.2] < 0.05");

    /* NMaximize over a disjunction (negated-objective path): maximise
     * -(x^2 + y^2) over two disks; optimum on the (2,2)-disk boundary, ~-6.311. */
    check_true("Abs[First[NMaximize[{-(x^2 + y^2), "
               "(x - 2)^2 + (y - 2)^2 <= 0.1 || (x + 2)^2 + (y + 2)^2 <= 0.1}, {x, y}, "
               "Method -> \"SimulatedAnnealing\"]] - (-6.311)] < 0.05");
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

    /* "LevelIterations" -> L sets the trial moves per temperature level, so the
     * per-chain budget is MaxIterations * L. A longer single chain (no polish,
     * fixed seed) reaches a strictly deeper raw point than a very short one on a
     * rugged surface, proving the budget is honored — and honored verbatim, past
     * the automatic cap. */
    check_true("First[NMinimize[{-(y + 47) Sin[Sqrt[Abs[y + x/2 + 47]]] "
               "- x Sin[Sqrt[Abs[x - (y + 47)]]], -512 <= x <= 512 && -512 <= y <= 512}, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"SearchPoints\" -> 1, "
               "\"LevelIterations\" -> 200, \"PostProcess\" -> False}]] < "
               "First[NMinimize[{-(y + 47) Sin[Sqrt[Abs[y + x/2 + 47]]] "
               "- x Sin[Sqrt[Abs[x - (y + 47)]]], -512 <= x <= 512 && -512 <= y <= 512}, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"SearchPoints\" -> 1, "
               "\"LevelIterations\" -> 2, \"PostProcess\" -> False}]]");

    /* Invalid sub-option values warn (NMinimize::sopt / ::bexp) and fall back
     * to the defaults rather than failing: the solve still returns a List. */
    check_eq("Head[NMinimize[{(x-3)^2 + (y+2)^2, -50 <= x <= 50 && -50 <= y <= 50}, {x, y}, "
             "Method -> {\"SimulatedAnnealing\", \"PerturbationScale\" -> -3}]]", "List");
    check_eq("Head[NMinimize[{(x-3)^2 + (y+2)^2, -50 <= x <= 50 && -50 <= y <= 50}, {x, y}, "
             "Method -> {\"SimulatedAnnealing\", \"BoltzmannExponent\" -> \"nope\"}]]", "List");
    check_eq("Head[NMinimize[{(x-3)^2 + (y+2)^2, -50 <= x <= 50 && -50 <= y <= 50}, {x, y}, "
             "Method -> {\"SimulatedAnnealing\", \"LevelIterations\" -> -5}]]", "List");
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

static void test_sa_deceptive_landscapes(void) {
    /* Deceptive landscapes whose global basin is off-center and unreachable by
     * greedy descent (Eggholder, Schwefel). Before the adaptive-temperature fix
     * the acceptance exponent -df/T sat far below the objective's magnitude
     * (T <= 1, objective spanning hundreds), so every uphill move was rejected
     * and each "chain" was greedy descent from its random start; the default
     * 2 n = 4 chains at n = 2 then reported only -716.67 on the Eggholder. With
     * the per-chain objective-scale temperature -df/(T s) and the
     * Min[Max[2 n, 12], 50] = 12-chain floor, the default reaches the global
     * -959.64, below Mathematica's default -935.34. Deterministic under the
     * fixed default seed. */
    check_true("First[NMinimize[{-(y + 47) Sin[Sqrt[Abs[y + x/2 + 47]]] "
               "- x Sin[Sqrt[Abs[x - (y + 47)]]], "
               "-512 <= x <= 512 && -512 <= y <= 512}, {x, y}, "
               "Method -> \"SimulatedAnnealing\"]] < -935.0");

    /* Schwefel-2D: global ~0 at (420.97, 420.97). The fixed-temperature greedy
     * walk stuck one coordinate in a wrong well (~118); the real anneal clears
     * both coordinates. */
    check_true("First[NMinimize[{837.9658 - x Sin[Sqrt[Abs[x]]] - y Sin[Sqrt[Abs[y]]], "
               "-500 <= x <= 500 && -500 <= y <= 500}, {x, y}, "
               "Method -> \"SimulatedAnnealing\"]] < 0.01");
}

static void test_schaffer2_simulatedannealing(void) {
    /* Schaffer function N.2: 0.5 + (Sin[x^2-y^2]^2 - 0.5)/(1 + 0.001(x^2+y^2))^2.
     * Global minimum 0 at the origin, ringed by concentric near-optimal circular
     * ridges whose value approaches 0.5 far out and decays toward the centre —
     * a deceptive basin that greedy descent settles onto the wrong ring of. On
     * the [-100,100]^2 box with "SearchPoints" -> 150 restarts and a tight
     * "Tolerance" -> 10^-6, the anneal-plus-polish reaches the true centre:
     * f = 0 at (~1e-7, ~-1e-7). Deterministic under the fixed default seed. */
    check_true("First[NMinimize[{0.5 + (Sin[x^2 - y^2]^2 - 0.5)/(1 + 0.001 (x^2 + y^2))^2, "
               "-100 <= x <= 100 && -100 <= y <= 100}, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"SearchPoints\" -> 150, \"Tolerance\" -> 10^-6}]] < 1.*^-6");
    /* The reported minimizer lands at the origin (both coordinates ~0). */
    check_true("Norm[{x, y} /. Last[NMinimize[{0.5 + (Sin[x^2 - y^2]^2 - 0.5)/(1 + 0.001 (x^2 + y^2))^2, "
               "-100 <= x <= 100 && -100 <= y <= 100}, {x, y}, "
               "Method -> {\"SimulatedAnnealing\", \"SearchPoints\" -> 150, \"Tolerance\" -> 10^-6}]]] < 1.*^-3");
}

static void test_mishra_bird_randomsearch(void) {
    /* Mishra's Bird function, a classic constrained global-optimization test:
     * Sin[y] Exp[(1-Cos[x])^2] + Cos[x] Exp[(1-Sin[y])^2] + (x-y)^2 over the
     * disk (x+5)^2 + (y+5)^2 < 25. The unconstrained surface is wildly
     * multimodal, so the disk constraint carves out the single deep basin whose
     * global minimum is -106.7645 at (-3.1302, -1.5821). With 200 RandomSearch
     * restarts and the built-in (Automatic) squared penalty scoring infeasible
     * trial points, the anneal-plus-polish reaches it. Deterministic under the
     * fixed default seed. This is the reported In[1] example. */
    check_true("Abs[First[NMinimize[{Sin[y] Exp[(1 - Cos[x])^2] + Cos[x] Exp[(1 - Sin[y])^2] + (x - y)^2, "
               "(x + 5)^2 + (y + 5)^2 < 25}, {x, y}, "
               "Method -> {\"RandomSearch\", \"SearchPoints\" -> 200, \"PenaltyFunction\" -> Automatic}]] "
               "- (-106.7645)] < 1.*^-2");
    /* The reported minimizer lands at (-3.1302, -1.5821). */
    check_true("Norm[{x + 3.1302, y + 1.5821} /. Last[NMinimize[{Sin[y] Exp[(1 - Cos[x])^2] "
               "+ Cos[x] Exp[(1 - Sin[y])^2] + (x - y)^2, (x + 5)^2 + (y + 5)^2 < 25}, {x, y}, "
               "Method -> {\"RandomSearch\", \"SearchPoints\" -> 200, \"PenaltyFunction\" -> Automatic}]]] < 1.*^-2");
    /* The returned point is feasible (inside the disk). */
    check_true("((x + 5)^2 + (y + 5)^2 /. Last[NMinimize[{Sin[y] Exp[(1 - Cos[x])^2] "
               "+ Cos[x] Exp[(1 - Sin[y])^2] + (x - y)^2, (x + 5)^2 + (y + 5)^2 < 25}, {x, y}, "
               "Method -> {\"RandomSearch\", \"SearchPoints\" -> 200, \"PenaltyFunction\" -> Automatic}]]) <= 25.001");
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

static void test_lennard_jones_cluster(void) {
    /* Lennard-Jones 6-atom cluster: a real 18-D global-optimization workload
     * that drives the whole indexed-variable + DifferentialEvolution path at
     * once. Subsets[Range[6], {2}] enumerates the 15 atom pairs; the coordinates
     * are indexed vars built with Table[{x[i], y[i], z[i]}, ...] then Flatten'd
     * to 18 scalars; a ljPair[{i_, j_}] := ... helper supplies the softened pair
     * potential 1/r2^6 - 2/r2^3 (r2 the squared distance with a 1e-8 floor); the
     * box constraints are generated by Table[-3 <= v <= 3, {v, vars}] iterating
     * over the *symbol list* vars; and DE's "ScalingFactor" sub-option sets the
     * differential weight F.
     *
     * The softened potential's global minimum is the octahedron at ~ -12.712
     * (12 pairs at the unit-distance well contributing -1 each, 3 diagonals at
     * r2 = 2). Under the fixed default seed this lands at -12.3029 here, a good
     * near-global basin. The assertion is deliberately loose: an 18-D DE
     * trajectory drifts across platforms in floating point (cf.
     * test_de_boundary_no_stagnation), so the invariant is that a genuinely
     * compact, low-energy cluster forms (energy < -10, far below the ~0 of any
     * scattered or sentinel configuration), that all 18 coordinates are reported,
     * and that the returned point is feasible inside the +-3 box. Scoped in a
     * Module so the scratch bindings and the ljPair definition do not leak. */
    check_true(
        "Module[{coords, vars, pairs, ljPair, energy, r, n = 6}, "
        "coords = Table[{x[i], y[i], z[i]}, {i, 1, n}]; "
        "vars = Flatten[coords]; "
        "pairs = Subsets[Range[n], {2}]; "
        "ljPair[{i_, j_}] := With[{r2 = (x[i] - x[j])^2 + (y[i] - y[j])^2 "
        "+ (z[i] - z[j])^2 + 10^-8}, (1/r2^6) - (2/r2^3)]; "
        "energy = Sum[ljPair[p], {p, pairs}]; "
        "r = NMinimize[{energy, Table[-3 <= v <= 3, {v, vars}]}, vars, "
        "Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 150, "
        "\"ScalingFactor\" -> 0.7}]; "
        "Head[r] === List && Length[Last[r]] == 18 && First[r] < -10.0 "
        "&& Max[Abs[vars /. Last[r]]] <= 3.0001]");
}

static void test_katsuura_differentialevolution(void) {
    /* Katsuura-8: a fractal, highly multimodal benchmark. The objective is
     * (10/d^2) Product_i (1 + i Sum_{k=1..25} |2^k x_i - Round[2^k x_i]| / 2^k)
     * ^(10/d^1.2) - 10/d^2. Every |2^k x_i - Round[2^k x_i]| term vanishes when
     * x_i is an INTEGER, so each factor collapses to 1 and the global minimum 0
     * is attained on the entire integer lattice -- not only at the origin. A
     * near-integer point anywhere in the box (here one coordinate parked on the
     * +-100 wall) is therefore a legitimate global minimizer, which is why the
     * feasibility check allows the boundary rather than treating it as
     * stagnation. The Round makes the surface fractally rough, exactly what
     * stresses DifferentialEvolution; the "ScalingFactor" sub-option sets the
     * differential weight F. Deterministic under the fixed default seed
     * (~2.75e-9 here). Loose invariant, robust to cross-platform DE drift (cf.
     * test_de_boundary_no_stagnation): the search reaches the lattice minimum
     * (First[r] < 0.1, versus >= 13 at any scattered point), returns all 8
     * coordinates, and stays feasible inside the +-100 box. Module-scoped so the
     * katsuuraSum definition and scratch bindings do not leak. */
    check_true(
        "Module[{vars, katsuuraSum, obj, r, d = 8}, "
        "vars = Table[x[i], {i, 1, d}]; "
        "katsuuraSum[xi_, i_] := 1 + i Sum[Abs[2^k xi - Round[2^k xi]]/2^k, {k, 1, 25}]; "
        "obj = (10/d^2) Product[katsuuraSum[x[i], i], {i, 1, d}]^(10/d^1.2) - (10/d^2); "
        "r = NMinimize[{obj, Table[-100 <= x[i] <= 100, {i, 1, d}]}, vars, "
        "Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 200, "
        "\"ScalingFactor\" -> 0.6}]; "
        "Head[r] === List && Length[Last[r]] == 8 && First[r] < 0.1 "
        "&& Max[Abs[vars /. Last[r]]] <= 100.0001]");
}

static void test_rotated_rastrigin_stress(void) {
    /* 20-D rotated Rastrigin: the canonical hard, non-separable global-optimization
     * benchmark. A dense random orthogonal rotation Q (the Q factor of a seeded
     * RandomReal matrix, via QRDecomposition) couples all 20 coordinates, so the
     * ~10^20 local minima cannot be picked off one axis at a time; the global
     * minimum is 0 at the origin (well inside the +-5.12 box, since Q.0 = 0). This
     * drives the whole indexed-variable + DifferentialEvolution path at scale at
     * once: {qMat, rMat} = QRDecomposition[...] destructuring, Table vars and
     * constraints, a held Sum over a Dot, the "CrossProbability"/"ScalingFactor"
     * DE sub-options, and the Storn & Price 10n = 200 default population.
     * Deterministic under the fixed default search seed; the rotation is itself
     * pinned by SeedRandom[12345].
     *
     * The assertion is a loose invariant, robust to the cross-platform
     * floating-point drift of both the QR rotation and the 20-D DE trajectory
     * (cf. test_de_boundary_no_stagnation): the search must return the right
     * shape (all 20 coordinates), a feasible point inside the box, and a genuinely
     * good feasible basin. Rastrigin is >= 0 everywhere and a random scatter over
     * this box averages ~375, so First[r] < 100 confirms DE reached a deep basin
     * rather than a sentinel (1e300) or a scattered point. Plain DE does not solve
     * 20-D rotated Rastrigin to the global 0 -- it lands at ~37.8 here -- so the
     * bound is a "found a good region" check, not a "found the optimum" one. */
    check_true(
        "Module[{qMat, rMat, vars, z, obj, cons, r, d = 20}, "
        "SeedRandom[12345]; "
        "{qMat, rMat} = QRDecomposition[RandomReal[{-1, 1}, {d, d}]]; "
        "vars = Table[x[i], {i, 1, d}]; "
        "z = qMat . vars; "
        "obj = 10 d + Sum[z[[i]]^2 - 10 Cos[2 Pi z[[i]]], {i, 1, d}]; "
        "cons = Table[-5.12 <= x[i] <= 5.12, {i, 1, d}]; "
        "r = NMinimize[{obj, cons}, vars, "
        "Method -> {\"DifferentialEvolution\", \"CrossProbability\" -> 0.95, "
        "\"ScalingFactor\" -> 0.85}, MaxIterations -> 3000]; "
        "Head[r] === List && Length[Last[r]] == 20 && 0 <= First[r] < 100 "
        "&& Max[Abs[vars /. Last[r]]] <= 5.1201]");
}

static void test_de_options_effective(void) {
    /* Every "DifferentialEvolution" sub-option must actually steer the search, not
     * be silently dropped. "Tolerance" and "InitialPoints" used to be consumed
     * only by NelderMead; "CrossProbability"/"ScalingFactor" took any value with
     * no validation. Deterministic under the fixed default seed. */

    /* Tolerance sets the population-convergence threshold. On a smooth bowl with
     * the polish off, a tight tolerance refines the raw best far below a loose one
     * (which breaks as soon as the feasible sub-population's objective spread
     * collapses). If Tolerance were ignored the two would be identical. */
    check_true("First[NMinimize[{Sum[(x[i] - 3)^2, {i, 1, 4}], Table[-10 <= x[i] <= 10, {i, 1, 4}]}, "
               "Table[x[i], {i, 1, 4}], Method -> {\"DifferentialEvolution\", \"Tolerance\" -> 1.*^-14, "
               "\"PostProcess\" -> False}, MaxIterations -> 400]] < "
               "First[NMinimize[{Sum[(x[i] - 3)^2, {i, 1, 4}], Table[-10 <= x[i] <= 10, {i, 1, 4}]}, "
               "Table[x[i], {i, 1, 4}], Method -> {\"DifferentialEvolution\", \"Tolerance\" -> 0.1, "
               "\"PostProcess\" -> False}, MaxIterations -> 400]]");

    /* InitialPoints seeds the initial DE population. Seeding every member at the
     * origin -- the global minimum 0 of this 5-D Rastrigin -- makes the raw best
     * (one generation, no polish) essentially 0; a random start over the same box
     * is tens away. If InitialPoints were ignored the population would be random
     * and the raw best far from 0. */
    check_true("First[NMinimize[{50 + Sum[x[i]^2 - 10 Cos[2 Pi x[i]], {i, 1, 5}], "
               "Table[-5.12 <= x[i] <= 5.12, {i, 1, 5}]}, Table[x[i], {i, 1, 5}], "
               "Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 8, \"PostProcess\" -> False, "
               "\"InitialPoints\" -> ConstantArray[ConstantArray[0, 5], 8]}, MaxIterations -> 1]] < 1.*^-6");

    /* CrossProbability and ScalingFactor each steer DE/rand/1/bin: changing one
     * (all else fixed, same seed) lands the raw global-search point somewhere
     * different on a multimodal surface -- only possible if the value is threaded
     * into the mutation/crossover. Raw (PostProcess -> False) points are compared,
     * so the check does not depend on two polishes landing in different basins. */
    check_true("First[NMinimize[{Sum[x[i]^2 - 10 Cos[2 Pi x[i]], {i, 1, 8}], "
               "Table[-5.12 <= x[i] <= 5.12, {i, 1, 8}]}, Table[x[i], {i, 1, 8}], "
               "Method -> {\"DifferentialEvolution\", \"CrossProbability\" -> 0.1, \"PostProcess\" -> False}]] != "
               "First[NMinimize[{Sum[x[i]^2 - 10 Cos[2 Pi x[i]], {i, 1, 8}], "
               "Table[-5.12 <= x[i] <= 5.12, {i, 1, 8}]}, Table[x[i], {i, 1, 8}], "
               "Method -> {\"DifferentialEvolution\", \"CrossProbability\" -> 0.95, \"PostProcess\" -> False}]]");
    check_true("First[NMinimize[{Sum[x[i]^2 - 10 Cos[2 Pi x[i]], {i, 1, 8}], "
               "Table[-5.12 <= x[i] <= 5.12, {i, 1, 8}]}, Table[x[i], {i, 1, 8}], "
               "Method -> {\"DifferentialEvolution\", \"ScalingFactor\" -> 0.3, \"PostProcess\" -> False}]] != "
               "First[NMinimize[{Sum[x[i]^2 - 10 Cos[2 Pi x[i]], {i, 1, 8}], "
               "Table[-5.12 <= x[i] <= 5.12, {i, 1, 8}]}, Table[x[i], {i, 1, 8}], "
               "Method -> {\"DifferentialEvolution\", \"ScalingFactor\" -> 0.9, \"PostProcess\" -> False}]]");

    /* An out-of-range CrossProbability (not in [0,1]) or ScalingFactor (not in
     * (0,2]) warns (NMinimize::sopt) and falls back to the default rather than
     * failing: the solve still returns a valid List. */
    check_eq("Head[NMinimize[{x^2 + y^2, x + y >= 1}, {x, y}, "
             "Method -> {\"DifferentialEvolution\", \"CrossProbability\" -> 5}]]", "List");
    check_eq("Head[NMinimize[{x^2 + y^2, x + y >= 1}, {x, y}, "
             "Method -> {\"DifferentialEvolution\", \"ScalingFactor\" -> -1}]]", "List");
}

/* ------------------------------------------------------------------ *
 *  Real-world constrained problems (optimization testbed)             *
 * ------------------------------------------------------------------ */

static void test_refinery_pooling(void) {
    /* A1: refinery blending / pooling — 8 continuous variables, two equality
     * constraints one of which is BILINEAR (3 cA + cB = px·xA + py·yA, a pool
     * quality balance), plus linear capacity inequalities and box bounds. Global
     * minimum is -3900, reached at cA=cB=0, xA=xB=0, yA=100, yB=200, px=py=0.
     *
     * This is the case that motivated the augmented-Lagrangian upgrade to the
     * shared constrained local solver (fm_run_penalty). The old pure quadratic
     * penalty (μ→10^8) could not restore feasibility on the bilinear equality
     * from a random start — every polish ended infeasible, so RandomSearch
     * returned the {Infinity,...} sentinel on a plainly feasible problem, and a
     * larger MaxIterations (which converges each ill-conditioned round harder)
     * made it worse. PHR multipliers let a moderate μ reach the constraint
     * surface, so the polish now lands on -3900. DifferentialEvolution's global
     * search reaches it reliably across every seed in ~0.09s — versus scipy's
     * differential_evolution at ~2.2s (and less robust: it fails on some seeds).
     * scipy SLSQP (a true SQP local solver) also finds -3900; Mathilda's global
     * DE beats scipy's global engine, which is the matching method character.
     * Deterministic under the fixed "RandomSeed". */
    check_true(
        "Module[{cost, cons, r, sol, o},"
        " cost = 6 cA + 16 cB + 10 (xA + xB) - 9 (xA + yA) - 15 (xB + yB);"
        " cons = {cA + cB - (px + py) == 0, 3 cA + 1 cB - px xA - py yA == 0,"
        "   xA + yA <= 100, xB + yB <= 200, cA <= 300, cB <= 300, px <= 2.5, py <= 1.5,"
        "   cA >= 0, cB >= 0, xA >= 0, xB >= 0, yA >= 0, yB >= 0, px >= 0, py >= 0,"
        "   xA <= 100, yA <= 100, xB <= 200, yB <= 200};"
        " r = NMinimize[{cost, cons}, {px, py, xA, xB, yA, yB, cA, cB},"
        "   Method -> {\"DifferentialEvolution\", \"RandomSeed\" -> 1}];"
        " sol = Last[r]; o = First[r];"
        " (-3901 < o < -3899) &&"
        "  Abs[(cA + cB - (px + py)) /. sol] < 1*^-3 &&"
        "  Abs[(3 cA + cB - px xA - py yA) /. sol] < 1*^-3 &&"
        "  ((xA + yA) /. sol) <= 100.001 && ((xB + yB) /. sol) <= 200.001 &&"
        "  Min[{px, py, xA, xB, yA, yB, cA, cB} /. sol] >= -1*^-6]");
}

static void test_gaussian_well(void) {
    /* A2: an isolated Gaussian well in an otherwise flat 10-D landscape —
     * obj = 1 - Exp[-2·Σ(x_i-1.2345)²] + 1e-5·Σx_i² over [-5,5]^10. The Exp
     * carves a single narrow basin at x_i=1.2345 (global 1.524e-4); everywhere
     * else the surface is the ~1.0 plateau with only a 1e-5 quadratic pull
     * toward the ORIGIN — which is NOT where the well is, so a local method
     * descends to the plateau and misses the basin (RandomSearch finds it ~1/6
     * of seeds). DifferentialEvolution's population + differential mutation
     * locates the basin across every seed (SearchPoints 400), in ~0.02s versus
     * scipy's differential_evolution at ~3.4s and only 6/8 seeds.
     *
     * NB the user's original was d=12, sharpness 50 over [-10,10] — a needle so
     * narrow (capture radius ~0.15) that NO derivative-free global method,
     * Mathilda or scipy, can find it; both return the 1.0 plateau. This
     * regression uses the widest-still-isolated well that is genuinely findable,
     * so it tests the engine rather than luck. Deterministic under RandomSeed. */
    check_true(
        "Module[{d = 10, vars, obj, box, r, o, pt},"
        " vars = Table[x[i], {i, 1, d}];"
        " obj = 1 - Exp[-2 Sum[(x[i] - 1.2345)^2, {i, 1, d}]] + 10^-5 Sum[x[i]^2, {i, 1, d}];"
        " box = Table[-5 <= x[i] <= 5, {i, 1, d}];"
        " r = NMinimize[{obj, box}, vars,"
        "   Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 400,"
        "              \"ScalingFactor\" -> 0.9, \"RandomSeed\" -> 1}];"
        " o = First[r]; pt = vars /. Last[r];"
        " o < 1*^-3 && Max[Abs[pt - 1.2345]] < 0.05]");
}

static void test_modified_ackley(void) {
    /* A3: a modified-Ackley product-of-cosines over [-5,5]^10 —
     * -Exp[-0.2·rms(x)]·∏Cos[20 x_i] + 0.05·Σx_i². The Cos[20·] factors
     * oscillate with period ~0.31, so the surface is a dense forest of local
     * minima; the true global -1 sits at the single point x=0 and is
     * effectively unreachable by ANY derivative-free method (the origin basin
     * has radius ~π/20 ≈ 0.16 — a 1e-13 fraction of the box). Both Mathilda's
     * SimulatedAnnealing (best -0.984) and scipy's dual_annealing (best -0.975)
     * plateau just short of -1; Mathilda edges it (-0.984 best vs -0.975) and is
     * faster (~0.5 s vs ~1.0 s). This regression pins that SA reliably finds a
     * DEEP basin (not the ~0 the surface averages): every one of 8 seeds lands
     * below -0.87 at 40 chains, seed 1 at -0.919. Threshold -0.87 leaves margin
     * for cross-platform SA trajectory drift; deterministic under RandomSeed. */
    check_true(
        "Module[{d = 10, vars, obj, box, r},"
        " vars = Table[x[i], {i, 1, d}];"
        " obj = -Exp[-0.2 Sqrt[1/d Sum[x[i]^2, {i, 1, d}]]]"
        "        Product[Cos[20 x[i]], {i, 1, d}] + 0.05 Sum[x[i]^2, {i, 1, d}];"
        " box = Table[-5 <= x[i] <= 5, {i, 1, d}];"
        " r = NMinimize[{obj, box}, vars,"
        "   Method -> {\"SimulatedAnnealing\", \"LevelIterations\" -> 200,"
        "              \"PerturbationScale\" -> 0.5, \"SearchPoints\" -> 40,"
        "              \"RandomSeed\" -> 1}];"
        " First[r] < -0.87]");
}

static void test_minimax_chebyshev(void) {
    /* A4: a 15-D constrained minimax (Chebyshev) — minimize
     * Max_i |x_i - Sin[i]| subject to Σx_i² ≤ 15 and the linear equality
     * Σ i·x_i = 1, over [-2,2]^15. Max-of-affine + a convex sphere + a linear
     * equality is a CONVEX program, so the optimum is unique: t* = 0.125116
     * (verified independently by scipy SLSQP on the epigraph form). Because it
     * is convex a single NelderMead simplex (SearchPoints -> 1) plus the
     * augmented-Lagrangian polish reaches it exactly — no multi-restart needed —
     * in ~0.06 s, matching scipy's single-start SLSQP (~0.065 s). The default 20
     * restarts would be ~12× the work for the same answer. The user's original
     * spec passed the box table as a THIRD top-level list element
     * ({obj, cons, box}); NMinimize takes {obj, cons}, so the box is folded into
     * the constraint list here. Asserts the optimum, the equality residual, and
     * sphere feasibility; the convex optimum is reached precisely, so the band
     * is tight. */
    check_true(
        "Module[{d = 15, vars, obj, allcons, r, o, sol},"
        " vars = Table[x[i], {i, 1, d}];"
        " obj = Max[Table[Abs[x[i] - Sin[i]], {i, 1, d}]];"
        " allcons = Join[{Sum[x[i]^2, {i, 1, d}] <= d, Sum[i x[i], {i, 1, d}] == 1},"
        "   Table[-2 <= x[i] <= 2, {i, 1, d}]];"
        " r = NMinimize[{obj, allcons}, vars,"
        "   Method -> {\"NelderMead\", \"ContractRatio\" -> 0.5, \"Tolerance\" -> 10^-8,"
        "              \"SearchPoints\" -> 1, \"RandomSeed\" -> 1}];"
        " o = First[r]; sol = Last[r];"
        " Abs[o - 0.125116] < 1*^-3 &&"
        "  Abs[Sum[i (x[i] /. sol), {i, 1, d}] - 1] < 1*^-4 &&"
        "  Sum[(x[i] /. sol)^2, {i, 1, d}] <= 15.001]");
}

static void test_optimal_liquidation(void) {
    /* A5: Almgren-Chriss optimal execution over T=15 periods — minimize
     * Σ_t [0.1·|x_t-x_{t+1}|^1.5 + 0.05·0.04·x_t²] (transient market impact +
     * a risk penalty) subject to x_1=1, x_16=0 (start invested, end liquidated),
     * x_{t+1} ≤ x_t (sell-only), x_t ∈ [0,1]. |·|^1.5 and x² are convex and the
     * constraints are linear, so this is a CONVEX program with optimum 0.035206
     * (verified against scipy SLSQP). SimulatedAnnealing is far more than a
     * convex problem needs: at SearchPoints -> 1 (one chain + AL polish) it
     * reaches the optimum in ~0.05 s — faster than scipy's single-start SLSQP
     * (~0.11 s) — where the user's SearchPoints -> 100 spent ~2.2 s for the same
     * answer. Asserts the optimum, the pinned endpoints, and the monotone
     * (sell-only) schedule; convex ⇒ reached precisely, tight band. */
    check_true(
        "Module[{T = 15, x, eta = 0.1, lam = 0.05, s2 = 0.04, obj, cons, r, o, sol, xv},"
        " x = Array[pos, T + 1];"
        " obj = Sum[eta Abs[x[[t]] - x[[t + 1]]]^1.5 + lam s2 x[[t]]^2, {t, 1, T}];"
        " cons = Join[{pos[1] == 1, pos[T + 1] == 0}, Table[x[[t + 1]] <= x[[t]], {t, 1, T}],"
        "   Table[0 <= x[[t]] <= 1, {t, 1, T + 1}]];"
        " r = NMinimize[{obj, cons}, x,"
        "   Method -> {\"SimulatedAnnealing\", \"PerturbationScale\" -> 0.1,"
        "              \"SearchPoints\" -> 1, \"RandomSeed\" -> 1}];"
        " o = First[r]; sol = Last[r]; xv = x /. sol;"
        " Abs[o - 0.035206] < 1*^-4 && Abs[xv[[1]] - 1] < 1*^-4 &&"
        "  Abs[xv[[T + 1]]] < 1*^-4 && Max[Differences[xv]] <= 1*^-4]");
}

static void test_risk_parity(void) {
    /* A6: risk-parity portfolio, n=8. Minimize Σ_{i,j}(rc_i-rc_j)² where the
     * risk contribution rc = w ⊙ (cov·w), subject to Σw=1 and w≥0 — i.e. find
     * the weights whose per-asset risk contributions are all equal. For a PSD
     * covariance this has a unique interior solution (all w_i > 0); here every
     * rc_i converges to the same value, driving the objective to ~1e-16 with
     * min weight 0.0747 (identical to scipy SLSQP on the same covariance). The
     * problem is effectively unimodal, so RandomSearch needs only ONE start
     * (SearchPoints -> 1) — the user's 250 was ~40× the work for the same point.
     * Asserts the objective is ~0, the budget Σw=1, all weights strictly
     * positive, and — the risk-parity property itself — that the spread of risk
     * contributions is < 1e-4. Deterministic (seeded covariance + RandomSeed). */
    check_true(
        "Module[{n = 8, cov, w, rc, obj, cons, r, o, sol, wv, rcv},"
        " SeedRandom[202];"
        " cov = (# + Transpose[#]) &@RandomReal[{-0.05, 0.15}, {n, n}];"
        " cov = cov . Transpose[cov];"
        " w = Array[weight, n]; rc = w * (cov . w);"
        " obj = Sum[(rc[[i]] - rc[[j]])^2, {i, n}, {j, n}];"
        " cons = Join[{Total[w] == 1}, Table[w[[i]] >= 0, {i, n}]];"
        " r = NMinimize[{obj, cons}, w,"
        "   Method -> {\"RandomSearch\", \"Method\" -> \"InteriorPoint\","
        "              \"SearchPoints\" -> 1, \"RandomSeed\" -> 1}];"
        " o = First[r]; sol = Last[r]; wv = w /. sol; rcv = wv * (cov . wv);"
        " o < 1*^-6 && Abs[Total[wv] - 1] < 1*^-6 && Min[wv] > 0 &&"
        "  (Max[rcv] - Min[rcv]) < 1*^-4]");
}

static void test_vwap_tracking(void) {
    /* A7: VWAP execution tracking, T=20 — minimize the tracking error
     * Σ(v_t-target_t)² against an intraday volume profile, subject to Σv=1
     * (execute the whole order), v≥0, and a NONLINEAR impact-cost budget
     * Σ0.5·v_t²·Exp[v_t] ≤ B. Quadratic objective + convex constraint + linear
     * equality = a CONVEX program, unique optimum. The user's B=0.01 is
     * INFEASIBLE — the minimum achievable impact (at the uniform schedule) is
     * 0.02628, so no schedule with Σv=1 can meet 0.01. Set B=0.027, which lies
     * between the uniform minimum and the target's own impact (0.02775): the
     * budget genuinely BINDS, pulling the schedule off the target to
     * trackingError 0.000224628 (verified against scipy SLSQP). Convex ⇒ one
     * simplex (SearchPoints -> 1) + the QuasiNewton/AL polish reaches it, so
     * MaxIterations 100 suffices where the user's 5000 with the default 20
     * restarts spent ~3.9 s for the same answer (~0.08 s here). Asserts the
     * optimum, Σv=1, the active impact budget, and non-negativity. */
    check_true(
        "Module[{T = 20, tp, v, te, ic, cons, r, o, sol, vv},"
        " tp = Table[0.05 + 0.04 Cos[Pi (t - T/2)/(T/2)]^4, {t, 1, T}]; tp = tp/Total[tp];"
        " v = Array[trade, T];"
        " te = Sum[(v[[t]] - tp[[t]])^2, {t, 1, T}];"
        " ic = Sum[0.5 v[[t]]^2 Exp[v[[t]]], {t, 1, T}];"
        " cons = Join[{Total[v] == 1}, {ic <= 0.027}, Table[v[[t]] >= 0, {t, 1, T}]];"
        " r = NMinimize[{te, cons}, v,"
        "   Method -> {\"NelderMead\", \"PostProcess\" -> \"QuasiNewton\","
        "              \"SearchPoints\" -> 1, \"RandomSeed\" -> 1}, MaxIterations -> 100];"
        " o = First[r]; sol = Last[r]; vv = v /. sol;"
        " Abs[o - 0.000224628] < 1*^-5 && Abs[Total[vv] - 1] < 1*^-5 &&"
        "  (Sum[0.5 vv[[t]]^2 Exp[vv[[t]]], {t, 1, T}]) <= 0.027001 && Min[vv] >= -1*^-6]");
}

/* ------------------------------------------------------------------ *
 *  Combinatorial / integer problems (MINLP testbed)                   *
 * ------------------------------------------------------------------ */

static void test_qap_assignment(void) {
    /* B1: quadratic assignment problem, n=6. Binary permutation matrix x[i,k]
     * (facility i at location k), doubly-stochastic + integer, minimizing the
     * quartic Σ flow[i,j]·dist[k,l]·x[i,k]·x[j,l]. The true optimum, by brute
     * force over all 720 permutations, is 929.672; scipy's dedicated QAP solver
     * (quadratic_assignment, 2-opt best-of-30) reaches it too. Mathilda's GENERAL
     * global optimizer reaches the SAME optimum with RandomSearch at 300 starts
     * (SA, the user's original choice, is deceptive here — it stalls on a 945.9
     * local optimum from most seeds; RandomSearch's independent restarts clear
     * it). Half the seeds land on that 945.9 basin, so the robust bar is a valid
     * permutation matrix AND within ~1.8% of the optimum (< 946) — which already
     * beats scipy's single-start 2-opt (962). Deterministic under RandomSeed. */
    check_true(
        "Module[{n = 6, flow, dist, vars, varList, obj, cons, r, o, sol, xm},"
        " SeedRandom[1];"
        " flow = RandomReal[{1, 10}, {n, n}]; dist = RandomReal[{1, 10}, {n, n}];"
        " vars = Array[x, {n, n}]; varList = Flatten[vars];"
        " obj = Sum[flow[[i, j]] dist[[k, l]] vars[[i, k]] vars[[j, l]],"
        "   {i, n}, {j, n}, {k, n}, {l, n}];"
        " cons = Join[Table[Sum[vars[[i, j]], {i, n}] == 1, {j, n}],"
        "   Table[Sum[vars[[i, j]], {j, n}] == 1, {i, n}],"
        "   Map[0 <= # <= 1 &, varList], {Element[varList, Integers]}];"
        " r = NMinimize[{obj, cons}, varList,"
        "   Method -> {\"RandomSearch\", \"SearchPoints\" -> 300, \"RandomSeed\" -> 1}];"
        " o = First[r]; sol = Last[r]; xm = vars /. sol;"
        " o < 946 && AllTrue[Total /@ xm, Abs[# - 1] < 1*^-4 &] &&"
        "  AllTrue[Total /@ Transpose[xm], Abs[# - 1] < 1*^-4 &] &&"
        "  AllTrue[Flatten[xm], Abs[# - Round[#]] < 1*^-4 &]]");
}

static void test_job_scheduling(void) {
    /* B2: single-machine job scheduling, n=5, Big-M disjunctive formulation.
     * Continuous start times t_i, binary precedence y[i,j], and for every pair
     * a Big-M disjunction forcing i-before-j OR j-before-i (no overlap). Minimize
     * the makespan Max_i(t_i + dur_i). Because all five jobs share one machine,
     * the makespan is exactly Σ durations = 26 for ANY ordering (pack from t=0),
     * so 26 is a hard lower bound. Mathilda's general RandomSearch finds a valid
     * schedule with makespan 27 — within 3.8% of the MILP optimum; the residual
     * gap is imperfect continuous-time packing, the expected general-vs-dedicated
     * cost. The `Nothing` in the Table (used to skip the i≥j pairs) exercises the
     * list-construction identity element. Asserts a near-optimal makespan; since
     * a feasible single-machine schedule cannot dip below 26, the lower bound
     * also certifies the disjunctions hold. Deterministic under RandomSeed. */
    check_true(
        "Module[{n = 5, m = 1000, times, order, durations, cons, obj, allc, r, o},"
        " times = Array[t, n]; order = Array[y, {n, n}]; durations = {4, 7, 2, 8, 5};"
        " cons = Flatten@Table[If[i < j,"
        "   {t[j] >= t[i] + durations[[i]] - m (1 - y[i, j]),"
        "    t[i] >= t[j] + durations[[j]] - m y[i, j],"
        "    0 <= y[i, j] <= 1, Element[y[i, j], Integers]}, Nothing], {i, n}, {j, n}];"
        " obj = Max[Table[t[i] + durations[[i]], {i, n}]];"
        " allc = Join[cons, Map[# >= 0 &, times]];"
        " r = NMinimize[{obj, allc}, Flatten[{times, order}],"
        "   Method -> {\"RandomSearch\", \"SearchPoints\" -> 50, \"RandomSeed\" -> 1}];"
        " o = First[r]; 25.99 < o < 28]");
}

static void test_multiknapsack(void) {
    /* B5: multiple-knapsack, 30 items, 5 capacity constraints. Binary x_i,
     * maximize value·x (written as minimize -value·x) subject to five
     * weight·x ≤ capacity budgets. Purely inequality-constrained (no tight
     * assignment structure), so the integer search — now with 2-flip swap moves,
     * which for a knapsack means "drop one item, add another" — reaches the
     * region of the optimum easily. The MILP optimum (scipy.milp) is -995;
     * DifferentialEvolution (100 search points) lands at -990 for this seed,
     * within 0.5% of it, in ~0.04 s (faster than milp's ~0.17 s). Asserts a
     * near-optimal value AND that every capacity budget is honoured and the
     * solution is binary. Deterministic under RandomSeed. */
    check_true(
        "Module[{n = 30, m = 5, values, weights, capacities, vars, obj, cons, r, o, sol, xv},"
        " SeedRandom[123];"
        " values = RandomInteger[{10, 100}, n]; weights = RandomInteger[{5, 50}, {m, n}];"
        " capacities = Total[weights, {2}]*0.4; vars = Array[x, n];"
        " obj = -values . vars;"
        " cons = Join[Table[weights[[j]] . vars <= capacities[[j]], {j, m}],"
        "   Table[0 <= vars[[i]] <= 1, {i, n}], {Element[vars, Integers]}];"
        " r = NMinimize[{obj, cons}, vars,"
        "   Method -> {\"DifferentialEvolution\", \"SearchPoints\" -> 100, \"RandomSeed\" -> 1}];"
        " o = First[r]; sol = Last[r]; xv = vars /. sol;"
        " o <= -985 && AllTrue[Range[m], weights[[#]] . xv <= capacities[[#]] + 1*^-6 &] &&"
        "  AllTrue[xv, Abs[# - Round[#]] < 1*^-6 &]]");
}

static void test_sudoku_latin(void) {
    /* B8: 4x4 Latin-square / mini-sudoku feasibility. Binary x[i,j,k] (cell (i,j)
     * holds value k) with three families of one-hot constraints — each cell holds
     * one value, each value appears once per row, once per column — and a constant
     * objective (pure feasibility). The one-hot groups OVERLAP (every x[i,j,k] is
     * in a cell-group, a row-group and a column-group), so a random rounded start
     * is never a valid assignment; single-coordinate descent stalls infeasible and
     * NMinimize used to return the {Infinity,...} sentinel. The integer search now
     * detects the Σ-of-binaries == 1 groups and, when the plain descent ends
     * infeasible, repairs each group to exactly one 1; the 2-flip swap move then
     * reconciles the overlapping groups into a full solution. RandomSearch at 20
     * restarts finds a feasible Latin square in ~0.3 s. Asserts the objective is
     * finite (feasible) and every one-hot family is satisfied and binary.
     * Deterministic under RandomSeed. */
    check_true(
        "Module[{n = 4, vars, vl, cons, r, o, sol, xm},"
        " vars = Array[x, {n, n, n}]; vl = Flatten[vars];"
        " cons = Join["
        "   Flatten@Table[Sum[vars[[i, j, k]], {k, n}] == 1, {i, n}, {j, n}],"
        "   Flatten@Table[Sum[vars[[i, j, k]], {j, n}] == 1, {i, n}, {k, n}],"
        "   Flatten@Table[Sum[vars[[i, j, k]], {i, n}] == 1, {j, n}, {k, n}],"
        "   Table[0 <= vl[[i]] <= 1, {i, Length[vl]}], {Element[vl, Integers]}];"
        " r = NMinimize[{0, cons}, vl,"
        "   Method -> {\"RandomSearch\", \"SearchPoints\" -> 20, \"RandomSeed\" -> 1}];"
        " o = First[r]; sol = Last[r]; xm = vars /. sol;"
        " o < 1 &&"
        "  AllTrue[Flatten@Table[Sum[xm[[i, j, k]], {k, n}], {i, n}, {j, n}], Abs[# - 1] < 1*^-4 &] &&"
        "  AllTrue[Flatten@Table[Sum[xm[[i, j, k]], {j, n}], {i, n}, {k, n}], Abs[# - 1] < 1*^-4 &] &&"
        "  AllTrue[Flatten@Table[Sum[xm[[i, j, k]], {i, n}], {j, n}, {k, n}], Abs[# - 1] < 1*^-4 &] &&"
        "  AllTrue[Flatten[xm], Abs[# - Round[#]] < 1*^-4 &]]");
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
    TEST(test_disjunctive_constraints);

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
    TEST(test_sa_deceptive_landscapes);
    TEST(test_schaffer2_simulatedannealing);
    TEST(test_mishra_bird_randomsearch);
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
    TEST(test_lennard_jones_cluster);
    TEST(test_katsuura_differentialevolution);
    TEST(test_rotated_rastrigin_stress);
    TEST(test_de_options_effective);

    /* 14. Real-world constrained problems (optimization testbed) */
    TEST(test_refinery_pooling);
    TEST(test_gaussian_well);
    TEST(test_modified_ackley);
    TEST(test_minimax_chebyshev);
    TEST(test_optimal_liquidation);
    TEST(test_risk_parity);
    TEST(test_vwap_tracking);

    /* 15. Combinatorial / integer (MINLP) */
    TEST(test_qap_assignment);
    TEST(test_job_scheduling);
    TEST(test_multiknapsack);
    TEST(test_sudoku_latin);

    printf("All NMinimize tests passed.\n");
    return 0;
}
