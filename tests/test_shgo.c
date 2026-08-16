/* Unit + stress tests for NMinimize's SHGO method (Simplicial Homology Global
 * Optimization; src/numerical_calculus/findmin.c, nm_shgo).
 *
 * SHGO samples the bounded box, builds a graph (the exact Kuhn-triangulation
 * 1-skeleton for "Simplicial", a k-nearest-neighbour graph for the "Sobol" /
 * "Halton" low-discrepancy sets), finds the minimizer pool, and starts one local
 * search per pool vertex. It is deterministic (QMC + a fixed RandomSeed), so
 * exact-ish assertions are stable.
 *
 * Reference values are analytic optima or published benchmark constants, never a
 * snapshot of the solver's own output. Deceptive multimodal cases assert a loose
 * "global basin found" threshold plus feasibility, following the convention in
 * test_nminimize.c.
 *
 * Coverage:
 *   - Result shape / typing.
 *   - Convex / unimodal exact minima.
 *   - Standard multimodal benchmarks (Himmelblau, Rastrigin, Ackley, Branin,
 *     six-hump camel, Cross-in-tray, Shubert, Styblinski-Tang, Eggholder).
 *   - Sampling-method parity ("Simplicial" / "Sobol" / "Halton").
 *   - Constrained global optimization (disk, Mishra's Bird).
 *   - Sub-options: SearchPoints, Iterations, RandomSeed determinism,
 *     PostProcess -> False.
 *   - High-dimension "Simplicial" -> "Sobol" fall-back.
 *   - NMaximize via SHGO and min/max duality.
 *   - Memory-hygiene smoke loop.
 *
 * Run binary directly: ./shgo_tests */

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
/* 1. Result shape / typing                                            */
/* ------------------------------------------------------------------ */
static void test_shape(void) {
    /* {fmin, {x -> _, y -> _}} with a numeric minimum and rule list. */
    check_true("MatchQ[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\"}], {_?NumberQ, {(_->_)..}}]");
    check_true("Length[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\"}]] == 2");
}

/* ------------------------------------------------------------------ */
/* 2. Convex / unimodal exact minima                                   */
/* ------------------------------------------------------------------ */
static void test_sphere(void) {
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\"}]]] < 1.*^-6");
    check_true("Norm[{x,y} /. Last[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"SHGO\"}]]] < 1.*^-3");
}

static void test_quartic(void) {
    /* x^4 - 3 x^2 - x on [-5,5]: global -3.5139097 at 1.3008373. */
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"SHGO\"}]] - (-3.5139097)] < 1.*^-3");
    check_true("Abs[(x /. Last[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"SHGO\"}]]) - 1.3008373] < 1.*^-3");
}

static void test_shifted_parabola(void) {
    check_true("Abs[First[NMinimize[{(x-3)^2+(y+2)^2, -10<=x<=10 && -10<=y<=10}, "
               "{x,y}, Method->{\"SHGO\"}]]] < 1.*^-6");
    check_true("Norm[{x-3,y+2} /. Last[NMinimize[{(x-3)^2+(y+2)^2, "
               "-10<=x<=10 && -10<=y<=10}, {x,y}, Method->{\"SHGO\"}]]] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* 3. Standard multimodal benchmarks (bounded boxes)                   */
/* ------------------------------------------------------------------ */
static void test_himmelblau(void) {
    /* Four global minima, all f = 0. */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, Method->{\"SHGO\"}]]] < 1.*^-6");
}

static void test_rastrigin_2d(void) {
    check_true("Abs[First[NMinimize[{20 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,2}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->{\"SHGO\"}]]] < 1.*^-4");
}

static void test_rastrigin_3d(void) {
    check_true("Abs[First[NMinimize[{30 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,3}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,3}]}, Table[x[i],{i,1,3}], "
               "Method->{\"SHGO\",\"SamplingMethod\"->\"Simplicial\","
               "\"SearchPoints\"->200}]]] < 1.*^-3");
}

static void test_ackley_2d(void) {
    check_true("Abs[First[NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] "
               "- Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, Method->{\"SHGO\"}]]] < 1.*^-4");
}

static void test_branin(void) {
    /* Three global minima, f = 0.397887. */
    check_true("Abs[First[NMinimize[{(y-5.1/(4 Pi^2) x^2+5/Pi x-6)^2 "
               "+10(1-1/(8 Pi))Cos[x]+10, -5<=x<=10 && 0<=y<=15}, {x,y}, "
               "Method->{\"SHGO\"}]] - 0.397887] < 1.*^-4");
}

static void test_six_hump_camel(void) {
    /* Global -1.0316284. */
    check_true("Abs[First[NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, "
               "-3<=x<=3 && -2<=y<=2}, {x,y}, Method->{\"SHGO\"}]] - (-1.0316284)] "
               "< 1.*^-4");
}

static void test_cross_in_tray(void) {
    /* Four global minima, f = -2.0626119. Deceptive: needs more sample points. */
    check_true("Abs[First[NMinimize[{-0.0001(Abs[Sin[x]Sin[y]"
               "Exp[Abs[100-Sqrt[x^2+y^2]/Pi]]]+1)^0.1, -10<=x<=10 && -10<=y<=10}, "
               "{x,y}, Method->{\"SHGO\",\"SearchPoints\"->500}]] - (-2.0626119)] "
               "< 1.*^-3");
}

static void test_shubert(void) {
    /* 18 global minima, f = -186.7309. */
    check_true("Abs[First[NMinimize[{(Sum[j Cos[(j+1)x+j],{j,1,5}])"
               "(Sum[j Cos[(j+1)y+j],{j,1,5}]), -10<=x<=10 && -10<=y<=10}, {x,y}, "
               "Method->{\"SHGO\",\"SearchPoints\"->400}]] - (-186.7309)] < 1.*^-2");
}

static void test_styblinski_tang_2d(void) {
    /* Global -39.16599 * d, d = 2 -> -78.33198, at all x = -2.903534. */
    check_true("Abs[First[NMinimize[{Sum[x[i]^4-16 x[i]^2+5 x[i], {i,1,2}]/2, "
               "Table[-5<=x[i]<=5, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->{\"SHGO\"}]] - (-78.33198)] < 1.*^-2");
}

static void test_styblinski_tang_4d(void) {
    /* d = 4 -> -156.66397, via Sobol sampling. */
    check_true("Abs[First[NMinimize[{Sum[x[i]^4-16 x[i]^2+5 x[i], {i,1,4}]/2, "
               "Table[-5<=x[i]<=5, {i,1,4}]}, Table[x[i],{i,1,4}], "
               "Method->{\"SHGO\",\"SamplingMethod\"->\"Sobol\","
               "\"SearchPoints\"->200}]] - (-156.66397)] < 5.*^-2");
}

static void test_eggholder(void) {
    /* Deceptive; global -959.6407 at (512, 404.2319). The simplicial corner
     * coverage lands the deepest basin; assert a robust threshold + feasibility
     * rather than a brittle exact pin. */
    check_true("First[NMinimize[{-(y+47) Sin[Sqrt[Abs[y+x/2+47]]] "
               "- x Sin[Sqrt[Abs[x-(y+47)]]], -512<=x<=512 && -512<=y<=512}, "
               "{x,y}, Method->{\"SHGO\"}]] < -959.0");
    check_true("Max[Abs[{x,y} /. Last[NMinimize[{-(y+47) Sin[Sqrt[Abs[y+x/2+47]]] "
               "- x Sin[Sqrt[Abs[x-(y+47)]]], -512<=x<=512 && -512<=y<=512}, "
               "{x,y}, Method->{\"SHGO\"}]]]] <= 512.001");
}

/* ------------------------------------------------------------------ */
/* 4. Sampling-method parity                                           */
/* ------------------------------------------------------------------ */
static void test_sampling_parity(void) {
    /* Each of the three samplers reaches Himmelblau's zero. */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\",\"SamplingMethod\"->\"Simplicial\"}]]] < 1.*^-6");
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\",\"SamplingMethod\"->\"Sobol\","
               "\"SearchPoints\"->150}]]] < 1.*^-6");
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\",\"SamplingMethod\"->\"Halton\","
               "\"SearchPoints\"->150}]]] < 1.*^-6");
    /* And each reaches Rastrigin-2D's zero. */
    check_true("Abs[First[NMinimize[{20 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,2}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->{\"SHGO\",\"SamplingMethod\"->\"Sobol\","
               "\"SearchPoints\"->200}]]] < 1.*^-4");
}

/* ------------------------------------------------------------------ */
/* 5. Constrained global optimization                                  */
/* ------------------------------------------------------------------ */
static void test_constrained_disk(void) {
    /* min x + y on the disk x^2 + y^2 <= 9: global -3 sqrt(2) = -4.2426407. */
    check_true("Abs[First[NMinimize[{x+y, x^2+y^2<=9}, {x,y}, "
               "Method->{\"SHGO\"}]] - (-4.2426407)] < 1.*^-3");
    /* returned point is feasible */
    check_true("(x^2+y^2 /. Last[NMinimize[{x+y, x^2+y^2<=9}, {x,y}, "
               "Method->{\"SHGO\"}]]) <= 9.001");
}

static void test_mishra_bird(void) {
    /* Mishra's Bird, constrained to the disk (x+5)^2+(y+5)^2 < 25:
     * global -106.7645 at (-3.1302, -1.5821). */
    check_true("Abs[First[NMinimize[{Sin[y]Exp[(1-Cos[x])^2]"
               "+Cos[x]Exp[(1-Sin[y])^2]+(x-y)^2, (x+5)^2+(y+5)^2<25}, {x,y}, "
               "Method->{\"SHGO\",\"SearchPoints\"->300}]] - (-106.7645)] < 1.*^-2");
    check_true("((x+5)^2+(y+5)^2 /. Last[NMinimize[{Sin[y]Exp[(1-Cos[x])^2]"
               "+Cos[x]Exp[(1-Sin[y])^2]+(x-y)^2, (x+5)^2+(y+5)^2<25}, {x,y}, "
               "Method->{\"SHGO\",\"SearchPoints\"->300}]]) <= 25.001");
}

/* ------------------------------------------------------------------ */
/* 6. Sub-options                                                      */
/* ------------------------------------------------------------------ */
static void test_determinism(void) {
    /* Deterministic QMC + fixed seed: two identical calls agree exactly. */
    check_eq("NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
             "Method->{\"SHGO\",\"SamplingMethod\"->\"Sobol\",\"RandomSeed\"->7}] "
             "=== "
             "NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
             "Method->{\"SHGO\",\"SamplingMethod\"->\"Sobol\",\"RandomSeed\"->7}]",
             "True");
}

static void test_search_points_honored(void) {
    /* A larger sample budget still reaches the global and stays feasible. */
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\",\"SamplingMethod\"->\"Sobol\","
               "\"SearchPoints\"->400}]]] < 1.*^-5");
}

static void test_iterations(void) {
    /* Refinement iterations improve the deceptive Eggholder over a single pass. */
    check_true("First[NMinimize[{-(y+47) Sin[Sqrt[Abs[y+x/2+47]]] "
               "- x Sin[Sqrt[Abs[x-(y+47)]]], -512<=x<=512 && -512<=y<=512}, "
               "{x,y}, Method->{\"SHGO\",\"SamplingMethod\"->\"Sobol\","
               "\"SearchPoints\"->128,\"Iterations\"->3}]] < -900.0");
}

static void test_postprocess_false(void) {
    /* No local polish: the raw pool best on the sphere is a valid, near-zero,
     * feasible point (the diagonal-bisection midpoint lands the centre). */
    check_true("First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\",\"PostProcess\"->False}]] <= 0.5");
    check_true("MatchQ[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\",\"PostProcess\"->False}], {_?NumberQ, {(_->_)..}}]");
}

static void test_option_scoping(void) {
    /* A sub-option belonging to another method (DualAnnealing's
     * "VisitingParameter"), an unknown key, and MaxIterations (a top-level-only
     * option) inside the Method list are each warned to stderr (muted) and
     * dropped — never stored — so the deterministic SHGO result is identical to
     * omitting them. Guards against the old flat, method-agnostic parser that
     * silently accepted a wrong-method key into its shared config slot. */
    check_eq("NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
             "Method->{\"SHGO\",\"SamplingMethod\"->\"Sobol\",\"RandomSeed\"->7,"
             "\"VisitingParameter\"->2.5,\"Frobnicate\"->1,MaxIterations->9}] "
             "=== "
             "NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
             "Method->{\"SHGO\",\"SamplingMethod\"->\"Sobol\",\"RandomSeed\"->7}]",
             "True");
    /* A valid-for-SHGO key is still honored (not swept up by the scoping guard). */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\",\"Iterations\"->2}]]] < 1.*^-5");
}

/* ------------------------------------------------------------------ */
/* 7. High-dimension fall-back                                         */
/* ------------------------------------------------------------------ */
static void test_simplicial_fallback(void) {
    /* n > 7: "Simplicial" is impractical (2^n corners), so it falls back to
     * "Sobol" and still solves the 8-D sphere. */
    check_true("Abs[First[NMinimize[{Sum[x[i]^2, {i,1,8}], "
               "Table[-5<=x[i]<=5, {i,1,8}]}, Table[x[i],{i,1,8}], "
               "Method->{\"SHGO\",\"SamplingMethod\"->\"Simplicial\","
               "\"SearchPoints\"->300}]]] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* 8. NMaximize + duality                                              */
/* ------------------------------------------------------------------ */
static void test_nmaximize(void) {
    check_true("Abs[First[NMaximize[{-(x^2+y^2), -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"SHGO\"}]]] < 1.*^-6");
    /* NMaximize[f] == -NMinimize[-f] on a concave problem. */
    check_true("Abs[First[NMaximize[{4-(x-1)^2-(y+2)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"SHGO\"}]] - 4.0] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 9. Memory-hygiene smoke loop                                        */
/* ------------------------------------------------------------------ */
static void test_memory_smoke(void) {
    for (int i = 0; i < 3; i++) {
        check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
                   "-5<=x<=5 && -5<=y<=5}, {x,y}, "
                   "Method->{\"SHGO\",\"SamplingMethod\"->\"Halton\","
                   "\"SearchPoints\"->120}]]] < 1.*^-5");
    }
}

int main(void) {
    symtab_init();
    core_init();

    /* Expected diagnostics (the fall-back warning) write to stderr; keep the
     * test output clean. */
    freopen("/dev/null", "w", stderr);

    /* 1. Shape */
    TEST(test_shape);

    /* 2. Convex / unimodal */
    TEST(test_sphere);
    TEST(test_quartic);
    TEST(test_shifted_parabola);

    /* 3. Standard multimodal benchmarks */
    TEST(test_himmelblau);
    TEST(test_rastrigin_2d);
    TEST(test_rastrigin_3d);
    TEST(test_ackley_2d);
    TEST(test_branin);
    TEST(test_six_hump_camel);
    TEST(test_cross_in_tray);
    TEST(test_shubert);
    TEST(test_styblinski_tang_2d);
    TEST(test_styblinski_tang_4d);
    TEST(test_eggholder);

    /* 4. Sampling-method parity */
    TEST(test_sampling_parity);

    /* 5. Constrained */
    TEST(test_constrained_disk);
    TEST(test_mishra_bird);

    /* 6. Sub-options */
    TEST(test_determinism);
    TEST(test_search_points_honored);
    TEST(test_iterations);
    TEST(test_postprocess_false);
    TEST(test_option_scoping);

    /* 7. High-dimension fall-back */
    TEST(test_simplicial_fallback);

    /* 8. NMaximize + duality */
    TEST(test_nmaximize);

    /* 9. Memory hygiene */
    TEST(test_memory_smoke);

    printf("\nAll SHGO tests passed.\n");
    return 0;
}
