/* Unit + stress tests for NMinimize's BasinHopping method (Monte-Carlo
 * minimization; src/numerical_calculus/nm_basin_hopping.c).
 *
 * Basin Hopping (Wales & Doye 1997) perturbs the current point by a random
 * displacement, LOCALLY MINIMIZES the perturbed point ("quench"), then applies a
 * Metropolis accept/reject to the two locally-minimized energies, with an
 * adaptive step size that targets a fixed acceptance rate. It mirrors
 * scipy.optimize.basinhopping and is deterministic for a fixed RandomSeed.
 *
 * Reference values are analytic optima or published benchmark constants, never a
 * snapshot of the solver's own output. Because Mathilda quenches with its own
 * well-behaved local minimizer (which, unlike scipy's L-BFGS-B, does not overshoot
 * across a basin boundary on its first step), a single-run Basin Hopping walk
 * crosses widely-separated basins purely through the random displacement, so
 * genuinely multi-basin cases (Styblinski-Tang) use "SearchPoints" -> K
 * multi-start for robustness, following the seed/threshold convention in
 * test_dual_annealing.c and test_shgo.c.
 *
 * Coverage:
 *   - Result shape / typing.
 *   - Convex / unimodal exact minima.
 *   - Standard multimodal benchmarks (Himmelblau, Booth, Beale, Rosenbrock,
 *     six-hump camel, Ackley 2-D/4-D, Rastrigin 2-D/3-D, Styblinski-Tang).
 *   - Constrained global optimization (disk, linear constraint).
 *   - Mixed-integer domains.
 *   - Sub-options: Temperature, StepSize, StepInterval, TargetAcceptanceRate,
 *     StepFactor, SuccessIterations, SearchPoints, RandomSeed determinism,
 *     PostProcess -> False, out-of-range validation, per-method option scoping.
 *   - NMaximize via BasinHopping and min/max duality.
 *   - High-dimension stress.
 *   - Memory-hygiene smoke loop.
 *
 * Run binary directly: ./basin_hopping_tests */

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
               "Method->{\"BasinHopping\"}], {_?NumberQ, {(_->_)..}}]");
    check_true("Length[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"BasinHopping\"}]] == 2");
}

/* ------------------------------------------------------------------ */
/* 2. Convex / unimodal exact minima                                   */
/* ------------------------------------------------------------------ */
static void test_sphere(void) {
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"BasinHopping\"}]]] < 1.*^-6");
    check_true("Norm[{x,y} /. Last[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"BasinHopping\"}]]] < 1.*^-3");
}

static void test_quartic(void) {
    /* x^4 - 3 x^2 - x on [-5,5]: global -3.5139097 at 1.3008373. Two basins;
     * SearchPoints -> 4 makes the multi-start robust to the starting basin. */
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"BasinHopping\",\"SearchPoints\"->4}]] - (-3.5139097)] < 1.*^-3");
    check_true("Abs[(x /. Last[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"BasinHopping\",\"SearchPoints\"->4}]]) - 1.3008373] < 1.*^-3");
}

static void test_shifted_parabola(void) {
    check_true("Abs[First[NMinimize[{(x-3)^2+(y+2)^2, -10<=x<=10 && -10<=y<=10}, "
               "{x,y}, Method->{\"BasinHopping\"}]]] < 1.*^-6");
    check_true("Norm[{x-3,y+2} /. Last[NMinimize[{(x-3)^2+(y+2)^2, "
               "-10<=x<=10 && -10<=y<=10}, {x,y}, Method->{\"BasinHopping\"}]]] "
               "< 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* 3. Standard multimodal benchmarks (bounded boxes)                   */
/* ------------------------------------------------------------------ */
static void test_himmelblau(void) {
    /* Four global minima, all f = 0. */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, Method->{\"BasinHopping\"}]]] "
               "< 1.*^-5");
}

static void test_booth(void) {
    /* Global f = 0 at (1, 3). */
    check_true("Abs[First[NMinimize[{(x+2y-7)^2+(2x+y-5)^2, "
               "-10<=x<=10 && -10<=y<=10}, {x,y}, Method->{\"BasinHopping\"}]]] "
               "< 1.*^-5");
}

static void test_beale(void) {
    /* Global f = 0 at (3, 0.5). */
    check_true("Abs[First[NMinimize[{(1.5-x+x y)^2+(2.25-x+x y^2)^2"
               "+(2.625-x+x y^3)^2, -4.5<=x<=4.5 && -4.5<=y<=4.5}, {x,y}, "
               "Method->{\"BasinHopping\"}]]] < 1.*^-5");
}

static void test_rosenbrock(void) {
    /* Global f = 0 at (1, 1); the quench walks straight down the valley. */
    check_true("Abs[First[NMinimize[{100(y-x^2)^2+(1-x)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"BasinHopping\"}]]] < 1.*^-6");
}

static void test_rastrigin_2d(void) {
    /* Rastrigin's tiny central basin: the adaptive-step hops walk in toward the
     * origin global 0. Pin a RandomSeed for a stable, reproducible run. */
    check_true("Abs[First[NMinimize[{20 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,2}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->{\"BasinHopping\",\"RandomSeed\"->1}]]] < 1.*^-4");
}

static void test_rastrigin_3d(void) {
    check_true("Abs[First[NMinimize[{30 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,3}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,3}]}, Table[x[i],{i,1,3}], "
               "Method->{\"BasinHopping\",\"RandomSeed\"->1}]]] < 1.*^-3");
}

static void test_ackley_2d(void) {
    check_true("Abs[First[NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] "
               "- Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, Method->{\"BasinHopping\"}]]] "
               "< 1.*^-4");
}

static void test_ackley_4d(void) {
    check_true("Abs[First[NMinimize[{-20 Exp[-0.2 Sqrt[Sum[x[i]^2,{i,1,4}]/4]] "
               "- Exp[Sum[Cos[2 Pi x[i]],{i,1,4}]/4] + 20 + E, "
               "Table[-5<=x[i]<=5,{i,1,4}]}, Table[x[i],{i,1,4}], "
               "Method->{\"BasinHopping\"}]]] < 1.*^-4");
}

static void test_six_hump_camel(void) {
    /* Global -1.0316284. */
    check_true("Abs[First[NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, "
               "-3<=x<=3 && -2<=y<=2}, {x,y}, Method->{\"BasinHopping\"}]] "
               "- (-1.0316284)] < 1.*^-3");
}

static void test_styblinski_tang_2d(void) {
    /* Global -39.16599 * d, d = 2 -> -78.33198, at all x = -2.903534. Four
     * well-separated basins; multi-start (SearchPoints -> 6) makes the single
     * conservative-quench walk robust to its starting basin. */
    check_true("Abs[First[NMinimize[{Sum[x[i]^4-16 x[i]^2+5 x[i], {i,1,2}]/2, "
               "Table[-5<=x[i]<=5, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->{\"BasinHopping\",\"SearchPoints\"->6}]] - (-78.33198)] < 1.*^-2");
}

/* ------------------------------------------------------------------ */
/* 4. Constrained global optimization                                  */
/* ------------------------------------------------------------------ */
static void test_constrained_disk(void) {
    /* min x + y on the disk x^2 + y^2 <= 9: global -3 sqrt(2) = -4.2426407. */
    check_true("Abs[First[NMinimize[{x+y, x^2+y^2<=9}, {x,y}, "
               "Method->{\"BasinHopping\"}]] - (-4.2426407)] < 1.*^-3");
    check_true("(x^2+y^2 /. Last[NMinimize[{x+y, x^2+y^2<=9}, {x,y}, "
               "Method->{\"BasinHopping\"}]]) <= 9.001");
}

static void test_constrained_linear(void) {
    /* min x^2 + y^2 s.t. x + 2 y >= 4: global 16/5 = 3.2 at (0.8, 1.6). */
    check_true("Abs[First[NMinimize[{x^2+y^2, x+2y>=4 && -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"BasinHopping\"}]] - 3.2] < 1.*^-3");
    check_true("(x+2y /. Last[NMinimize[{x^2+y^2, x+2y>=4 && -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"BasinHopping\"}]]) >= 3.999");
}

/* ------------------------------------------------------------------ */
/* 5. Mixed-integer domains                                            */
/* ------------------------------------------------------------------ */
static void test_integer_domain(void) {
    /* min (x-3)^2 + (y+2)^2 over integers: global 0 at (3, -2). */
    check_true("Abs[First[NMinimize[{(x-3)^2+(y+2)^2, "
               "-10<=x<=10 && -10<=y<=10 && Element[{x,y},Integers]}, {x,y}, "
               "Method->{\"BasinHopping\"}]]] < 1.*^-6");
    check_true("Norm[({x,y} /. Last[NMinimize[{(x-3)^2+(y+2)^2, "
               "-10<=x<=10 && -10<=y<=10 && Element[{x,y},Integers]}, {x,y}, "
               "Method->{\"BasinHopping\"}]]) - {3,-2}] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 6. Sub-options                                                      */
/* ------------------------------------------------------------------ */
static void test_temperature(void) {
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"BasinHopping\",\"Temperature\"->2.5}]]] < 1.*^-6");
}

static void test_step_size(void) {
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"BasinHopping\",\"StepSize\"->1.0}]]] < 1.*^-5");
}

static void test_step_interval(void) {
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"BasinHopping\",\"StepInterval\"->10}]]] < 1.*^-5");
}

static void test_target_accept(void) {
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"BasinHopping\",\"TargetAcceptanceRate\"->0.3}]]] < 1.*^-5");
}

static void test_step_factor(void) {
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"BasinHopping\",\"StepFactor\"->0.8}]]] < 1.*^-5");
}

static void test_success_iterations(void) {
    /* Early-stop once the best stalls for N hops; the sphere is found long before. */
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"BasinHopping\",\"SuccessIterations\"->20}]]] < 1.*^-6");
}

static void test_search_points(void) {
    /* Multiple independent multi-start runs keep the Deb-best. */
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->{\"BasinHopping\",\"SearchPoints\"->6}]] - (-3.5139097)] "
               "< 1.*^-3");
}

static void test_determinism(void) {
    /* Deterministic PRNG + fixed seed: two identical calls agree exactly. */
    check_eq("NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
             "Method->{\"BasinHopping\",\"RandomSeed\"->5}] === "
             "NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
             "Method->{\"BasinHopping\",\"RandomSeed\"->5}]",
             "True");
}

static void test_postprocess_false(void) {
    /* No final driver polish: each hop already quenches, so the raw best on the
     * sphere is a valid, near-zero, feasible point. */
    check_true("First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"BasinHopping\",\"PostProcess\"->False}]] <= 0.5");
    check_true("MatchQ[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"BasinHopping\",\"PostProcess\"->False}], "
               "{_?NumberQ, {(_->_)..}}]");
}

static void test_option_validation(void) {
    /* Out-of-range sub-option values warn (to stderr, muted) and fall back to the
     * default, so the solve still succeeds. */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"BasinHopping\",\"Temperature\"->-1.0,"
               "\"StepSize\"->0,\"TargetAcceptanceRate\"->2.0,\"StepFactor\"->1.5}]]] "
               "< 1.*^-5");
}

static void test_option_scoping(void) {
    /* A sub-option belonging to another method (DualAnnealing's "VisitingParameter"),
     * an unknown key, and MaxIterations (a top-level-only option) inside the Method
     * list are each warned to stderr (muted) and dropped — never stored — so the
     * seed-deterministic result is identical to omitting them. */
    check_eq("NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
             "Method->{\"BasinHopping\",\"RandomSeed\"->5,"
             "\"VisitingParameter\"->2.5,\"Frobnicate\"->1,MaxIterations->9}] "
             "=== "
             "NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
             "Method->{\"BasinHopping\",\"RandomSeed\"->5}]",
             "True");
    /* A valid-for-BasinHopping key is still honored (not swept up by the guard). */
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"BasinHopping\",\"RandomSeed\"->5,"
               "\"StepSize\"->0.7}]]] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 7. NMaximize + duality                                              */
/* ------------------------------------------------------------------ */
static void test_nmaximize(void) {
    check_true("Abs[First[NMaximize[{-(x^2+y^2), -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"BasinHopping\"}]]] < 1.*^-6");
    check_true("Abs[First[NMaximize[{10-x^2-y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"BasinHopping\"}]] - 10.0] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 8. High-dimension stress                                            */
/* ------------------------------------------------------------------ */
static void test_ackley_5d(void) {
    check_true("First[NMinimize[{-20 Exp[-0.2 Sqrt[Sum[x[i]^2,{i,1,5}]/5]] "
               "- Exp[Sum[Cos[2 Pi x[i]],{i,1,5}]/5] + 20 + E, "
               "Table[-5<=x[i]<=5,{i,1,5}]}, Table[x[i],{i,1,5}], "
               "Method->{\"BasinHopping\"}]] < 0.01");
}

static void test_rastrigin_5d_stress(void) {
    /* Higher-dimension Rastrigin: assert a finite, feasible result in a loose
     * near-global interval (a single conservative-quench walk does not guarantee
     * the exact global here on every stream). */
    check_true("Module[{r = NMinimize[{50 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,5}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,5}]}, Table[x[i],{i,1,5}], "
               "Method->{\"BasinHopping\",\"RandomSeed\"->1}]}, "
               "Head[r] === List && 0 <= First[r] < 2.0]");
}

/* ------------------------------------------------------------------ */
/* 9. Memory-hygiene smoke loop                                        */
/* ------------------------------------------------------------------ */
static void test_memory_smoke(void) {
    for (int i = 0; i < 3; i++) {
        check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
                   "-5<=x<=5 && -5<=y<=5}, {x,y}, "
                   "Method->{\"BasinHopping\",\"SearchPoints\"->2}]]] < 1.*^-4");
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
    TEST(test_rosenbrock);
    TEST(test_rastrigin_2d);
    TEST(test_rastrigin_3d);
    TEST(test_ackley_2d);
    TEST(test_ackley_4d);
    TEST(test_six_hump_camel);
    TEST(test_styblinski_tang_2d);

    /* 4. Constrained */
    TEST(test_constrained_disk);
    TEST(test_constrained_linear);

    /* 5. Mixed-integer */
    TEST(test_integer_domain);

    /* 6. Sub-options */
    TEST(test_temperature);
    TEST(test_step_size);
    TEST(test_step_interval);
    TEST(test_target_accept);
    TEST(test_step_factor);
    TEST(test_success_iterations);
    TEST(test_search_points);
    TEST(test_determinism);
    TEST(test_postprocess_false);
    TEST(test_option_validation);
    TEST(test_option_scoping);

    /* 7. NMaximize + duality */
    TEST(test_nmaximize);

    /* 8. High-dimension stress */
    TEST(test_ackley_5d);
    TEST(test_rastrigin_5d_stress);

    /* 9. Memory hygiene */
    TEST(test_memory_smoke);

    printf("\nAll BasinHopping tests passed.\n");
    return 0;
}
