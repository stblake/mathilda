/* Unit + stress tests for NMinimize's DualAnnealing method (Generalized
 * Simulated Annealing; src/numerical_calculus/findmin.c, nm_dual_annealing).
 *
 * Dual Annealing draws each trial jump from a Tsallis q-generalized heavy-tailed
 * visiting distribution, accepts uphill moves by a generalized Metropolis rule,
 * reanneals the temperature, and runs a local search after each Markov chain. It
 * mirrors scipy.optimize.dual_annealing and is deterministic for a fixed
 * RandomSeed.
 *
 * Reference values are analytic optima or published benchmark constants, never a
 * snapshot of the solver's own output. Multimodal cases that default dual
 * annealing does not reliably solve to the exact global on every RNG stream
 * (Rastrigin is the classic one) either fix a RandomSeed known to reach it or
 * assert a loose "reached a deep basin" threshold plus feasibility, following the
 * convention in test_nminimize.c and test_shgo.c.
 *
 * Coverage:
 *   - Result shape / typing.
 *   - Convex / unimodal exact minima.
 *   - Standard multimodal benchmarks (Himmelblau, Rastrigin, Ackley, Schwefel,
 *     Styblinski-Tang, six-hump camel, Booth, Beale, Eggholder).
 *   - Constrained global optimization (disk, linear constraint).
 *   - Mixed-integer domains.
 *   - Sub-options: VisitingParameter, AcceptanceParameter, InitialTemperature,
 *     RestartTemperatureRatio, LocalSearch, SearchPoints, RandomSeed determinism,
 *     PostProcess -> False, out-of-range validation.
 *   - NMaximize via DualAnnealing and min/max duality.
 *   - High-dimension stress.
 *   - Memory-hygiene smoke loop.
 *
 * Run binary directly: ./dual_annealing_tests */

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
    check_true("MatchQ[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DualAnnealing\"}], {_?NumberQ, {(_->_)..}}]");
    check_true("Length[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DualAnnealing\"}]] == 2");
}

/* ------------------------------------------------------------------ */
/* 2. Convex / unimodal exact minima                                   */
/* ------------------------------------------------------------------ */
static void test_sphere(void) {
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DualAnnealing\"}]]] < 1.*^-6");
    check_true("Norm[{x,y} /. Last[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"DualAnnealing\"}]]] < 1.*^-3");
}

static void test_quartic(void) {
    /* x^4 - 3 x^2 - x on [-5,5]: global -3.5139097 at 1.3008373. */
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"DualAnnealing\"}]] - (-3.5139097)] < 1.*^-3");
    check_true("Abs[(x /. Last[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"DualAnnealing\"}]]) - 1.3008373] < 1.*^-3");
}

static void test_shifted_parabola(void) {
    check_true("Abs[First[NMinimize[{(x-3)^2+(y+2)^2, -10<=x<=10 && -10<=y<=10}, "
               "{x,y}, Method->{\"DualAnnealing\"}]]] < 1.*^-6");
    check_true("Norm[{x-3,y+2} /. Last[NMinimize[{(x-3)^2+(y+2)^2, "
               "-10<=x<=10 && -10<=y<=10}, {x,y}, Method->{\"DualAnnealing\"}]]] "
               "< 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* 3. Standard multimodal benchmarks (bounded boxes)                   */
/* ------------------------------------------------------------------ */
static void test_himmelblau(void) {
    /* Four global minima, all f = 0. */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, Method->{\"DualAnnealing\"}]]] "
               "< 1.*^-5");
}

static void test_booth(void) {
    /* Global f = 0 at (1, 3). */
    check_true("Abs[First[NMinimize[{(x+2y-7)^2+(2x+y-5)^2, "
               "-10<=x<=10 && -10<=y<=10}, {x,y}, Method->{\"DualAnnealing\"}]]] "
               "< 1.*^-5");
}

static void test_beale(void) {
    /* Global f = 0 at (3, 0.5). */
    check_true("Abs[First[NMinimize[{(1.5-x+x y)^2+(2.25-x+x y^2)^2"
               "+(2.625-x+x y^3)^2, -4.5<=x<=4.5 && -4.5<=y<=4.5}, {x,y}, "
               "Method->{\"DualAnnealing\"}]]] < 1.*^-5");
}

static void test_rastrigin_2d(void) {
    /* Rastrigin's tiny central basin: default dual annealing reaches it on many
     * but not all RNG streams, so pin a RandomSeed known to land the global 0
     * (RandomSeed determinism keeps this stable). */
    check_true("Abs[First[NMinimize[{20 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,2}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->{\"DualAnnealing\",\"RandomSeed\"->1}]]] < 1.*^-4");
}

static void test_rastrigin_3d(void) {
    check_true("Abs[First[NMinimize[{30 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,3}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,3}]}, Table[x[i],{i,1,3}], "
               "Method->{\"DualAnnealing\",\"RandomSeed\"->1}]]] < 1.*^-3");
}

static void test_ackley_2d(void) {
    check_true("Abs[First[NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] "
               "- Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, Method->{\"DualAnnealing\"}]]] "
               "< 1.*^-4");
}

static void test_ackley_4d(void) {
    check_true("Abs[First[NMinimize[{-20 Exp[-0.2 Sqrt[Sum[x[i]^2,{i,1,4}]/4]] "
               "- Exp[Sum[Cos[2 Pi x[i]],{i,1,4}]/4] + 20 + E, "
               "Table[-5<=x[i]<=5,{i,1,4}]}, Table[x[i],{i,1,4}], "
               "Method->{\"DualAnnealing\"}]]] < 1.*^-4");
}

static void test_schwefel_2d(void) {
    /* Deceptive; global ~0 at (420.97, 420.97) near the box boundary. The
     * heavy-tailed visiting distribution reaches it. */
    check_true("First[NMinimize[{837.9658 - x Sin[Sqrt[Abs[x]]] "
               "- y Sin[Sqrt[Abs[y]]], -500<=x<=500 && -500<=y<=500}, {x,y}, "
               "Method->{\"DualAnnealing\"}]] < 0.01");
}

static void test_styblinski_tang_2d(void) {
    /* Global -39.16599 * d, d = 2 -> -78.33198, at all x = -2.903534. */
    check_true("Abs[First[NMinimize[{Sum[x[i]^4-16 x[i]^2+5 x[i], {i,1,2}]/2, "
               "Table[-5<=x[i]<=5, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->{\"DualAnnealing\"}]] - (-78.33198)] < 1.*^-2");
}

static void test_six_hump_camel(void) {
    /* Global -1.0316284. */
    check_true("Abs[First[NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, "
               "-3<=x<=3 && -2<=y<=2}, {x,y}, Method->{\"DualAnnealing\"}]] "
               "- (-1.0316284)] < 1.*^-3");
}

static void test_eggholder(void) {
    /* Deceptive; the -959.64 global is missed by default dual annealing (scipy's
     * dual_annealing lands ~-894.6 too). Assert a robust "reached a deep basin"
     * threshold plus feasibility rather than the unreachable exact global. */
    check_true("First[NMinimize[{-(y+47) Sin[Sqrt[Abs[y+x/2+47]]] "
               "- x Sin[Sqrt[Abs[x-(y+47)]]], -512<=x<=512 && -512<=y<=512}, "
               "{x,y}, Method->{\"DualAnnealing\"}]] < -800.0");
    check_true("Max[Abs[{x,y} /. Last[NMinimize[{-(y+47) Sin[Sqrt[Abs[y+x/2+47]]] "
               "- x Sin[Sqrt[Abs[x-(y+47)]]], -512<=x<=512 && -512<=y<=512}, "
               "{x,y}, Method->{\"DualAnnealing\"}]]]] <= 512.001");
}

/* ------------------------------------------------------------------ */
/* 4. Constrained global optimization                                  */
/* ------------------------------------------------------------------ */
static void test_constrained_disk(void) {
    /* min x + y on the disk x^2 + y^2 <= 9: global -3 sqrt(2) = -4.2426407. */
    check_true("Abs[First[NMinimize[{x+y, x^2+y^2<=9}, {x,y}, "
               "Method->{\"DualAnnealing\"}]] - (-4.2426407)] < 1.*^-3");
    check_true("(x^2+y^2 /. Last[NMinimize[{x+y, x^2+y^2<=9}, {x,y}, "
               "Method->{\"DualAnnealing\"}]]) <= 9.001");
}

static void test_constrained_linear(void) {
    /* min x^2 + y^2 s.t. x + 2 y >= 4: global 16/5 = 3.2 at (0.8, 1.6). */
    check_true("Abs[First[NMinimize[{x^2+y^2, x+2y>=4 && -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"DualAnnealing\"}]] - 3.2] < 1.*^-3");
    check_true("(x+2y /. Last[NMinimize[{x^2+y^2, x+2y>=4 && -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"DualAnnealing\"}]]) >= 3.999");
}

/* ------------------------------------------------------------------ */
/* 5. Mixed-integer domains                                            */
/* ------------------------------------------------------------------ */
static void test_integer_domain(void) {
    /* min (x-3)^2 + (y+2)^2 over integers: global 0 at (3, -2). */
    check_true("Abs[First[NMinimize[{(x-3)^2+(y+2)^2, "
               "-10<=x<=10 && -10<=y<=10 && Element[{x,y},Integers]}, {x,y}, "
               "Method->{\"DualAnnealing\"}]]] < 1.*^-6");
    check_true("Norm[({x,y} /. Last[NMinimize[{(x-3)^2+(y+2)^2, "
               "-10<=x<=10 && -10<=y<=10 && Element[{x,y},Integers]}, {x,y}, "
               "Method->{\"DualAnnealing\"}]]) - {3,-2}] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 6. Sub-options                                                      */
/* ------------------------------------------------------------------ */
static void test_visiting_parameter(void) {
    /* q_v = 2 (Cauchy / fast SA) still solves the quartic. */
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"DualAnnealing\",\"VisitingParameter\"->2.0}]] "
               "- (-3.5139097)] < 1.*^-3");
}

static void test_acceptance_parameter(void) {
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"DualAnnealing\",\"AcceptanceParameter\"->-8.0}]] "
               "- (-3.5139097)] < 1.*^-3");
}

static void test_initial_temperature(void) {
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"DualAnnealing\",\"InitialTemperature\"->1000.0}]] "
               "- (-3.5139097)] < 1.*^-3");
}

static void test_restart_ratio(void) {
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"DualAnnealing\",\"RestartTemperatureRatio\"->0.0001}]] "
               "- (-3.5139097)] < 1.*^-3");
}

static void test_local_search_off(void) {
    /* Pure Generalized Simulated Annealing (no per-chain local search); the final
     * driver PostProcess polish still applies, so it reaches the basin. */
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"DualAnnealing\",\"LocalSearch\"->False}]] "
               "- (-3.5139097)] < 1.*^-2");
}

static void test_search_points_chains(void) {
    /* Multiple independent chains keep the Deb-best. */
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"DualAnnealing\",\"SearchPoints\"->4}]] - (-3.5139097)] "
               "< 1.*^-3");
}

static void test_determinism(void) {
    /* Deterministic PRNG + fixed seed: two identical calls agree exactly. */
    check_eq("NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
             "Method->{\"DualAnnealing\",\"RandomSeed\"->5}] === "
             "NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
             "Method->{\"DualAnnealing\",\"RandomSeed\"->5}]",
             "True");
}

static void test_postprocess_false(void) {
    /* No local polish: the raw GSA global-search best on the sphere is a valid,
     * near-zero, feasible point. */
    check_true("First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DualAnnealing\",\"PostProcess\"->False}]] <= 0.5");
    check_true("MatchQ[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DualAnnealing\",\"PostProcess\"->False}], "
               "{_?NumberQ, {(_->_)..}}]");
}

static void test_option_validation(void) {
    /* Out-of-range sub-option values warn (to stderr, muted) and fall back to the
     * default, so the solve still succeeds. */
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"DualAnnealing\",\"VisitingParameter\"->5.0,"
               "\"AcceptanceParameter\"->-2.0}]] - (-3.5139097)] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* 7. NMaximize + duality                                              */
/* ------------------------------------------------------------------ */
static void test_nmaximize(void) {
    check_true("Abs[First[NMaximize[{-(x^2+y^2), -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DualAnnealing\"}]]] < 1.*^-6");
    check_true("Abs[First[NMaximize[{10-x^2-y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DualAnnealing\"}]] - 10.0] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 8. High-dimension stress                                            */
/* ------------------------------------------------------------------ */
static void test_ackley_5d(void) {
    check_true("First[NMinimize[{-20 Exp[-0.2 Sqrt[Sum[x[i]^2,{i,1,5}]/5]] "
               "- Exp[Sum[Cos[2 Pi x[i]],{i,1,5}]/5] + 20 + E, "
               "Table[-5<=x[i]<=5,{i,1,5}]}, Table[x[i],{i,1,5}], "
               "Method->{\"DualAnnealing\"}]] < 0.01");
}

static void test_rastrigin_5d_stress(void) {
    /* Higher-dimension Rastrigin: assert a finite, feasible result in a loose
     * near-global interval (default dual annealing does not guarantee the exact
     * global here on every stream, matching scipy's behaviour). */
    check_true("Module[{r = NMinimize[{50 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,5}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,5}]}, Table[x[i],{i,1,5}], "
               "Method->{\"DualAnnealing\",\"RandomSeed\"->1}]}, "
               "Head[r] === List && 0 <= First[r] < 2.0]");
}

/* ------------------------------------------------------------------ */
/* 9. Memory-hygiene smoke loop                                        */
/* ------------------------------------------------------------------ */
static void test_memory_smoke(void) {
    for (int i = 0; i < 3; i++) {
        check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
                   "-5<=x<=5 && -5<=y<=5}, {x,y}, "
                   "Method->{\"DualAnnealing\",\"SearchPoints\"->2}]]] < 1.*^-4");
    }
}

int main(void) {
    symtab_init();
    core_init();

    /* Expected diagnostics (out-of-range sub-option warnings) write to stderr;
     * keep the test output clean. */
    freopen("/dev/null", "w", stderr);

    /* 1. Shape */
    TEST(test_shape);

    /* 2. Convex / unimodal */
    TEST(test_sphere);
    TEST(test_quartic);
    TEST(test_shifted_parabola);

    /* 3. Standard multimodal benchmarks */
    TEST(test_himmelblau);
    TEST(test_booth);
    TEST(test_beale);
    TEST(test_rastrigin_2d);
    TEST(test_rastrigin_3d);
    TEST(test_ackley_2d);
    TEST(test_ackley_4d);
    TEST(test_schwefel_2d);
    TEST(test_styblinski_tang_2d);
    TEST(test_six_hump_camel);
    TEST(test_eggholder);

    /* 4. Constrained */
    TEST(test_constrained_disk);
    TEST(test_constrained_linear);

    /* 5. Mixed-integer */
    TEST(test_integer_domain);

    /* 6. Sub-options */
    TEST(test_visiting_parameter);
    TEST(test_acceptance_parameter);
    TEST(test_initial_temperature);
    TEST(test_restart_ratio);
    TEST(test_local_search_off);
    TEST(test_search_points_chains);
    TEST(test_determinism);
    TEST(test_postprocess_false);
    TEST(test_option_validation);

    /* 7. NMaximize + duality */
    TEST(test_nmaximize);

    /* 8. High-dimension stress */
    TEST(test_ackley_5d);
    TEST(test_rastrigin_5d_stress);

    /* 9. Memory hygiene */
    TEST(test_memory_smoke);

    printf("\nAll DualAnnealing tests passed.\n");
    return 0;
}
