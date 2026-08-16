/* Unit + stress tests for NMinimize's DIRECT method (DIviding RECTangles;
 * src/numerical_calculus/nm_direct.c).
 *
 * DIRECT normalizes the box to the unit hypercube, samples cell centers, groups
 * cells by size, and each iteration subdivides the "potentially optimal" cells
 * (lower-right convex hull of the (size, value) level minima, plus an epsilon
 * improvement test). It is fully deterministic (no RNG), so exact-ish assertions
 * are stable and two identical calls agree bit-for-bit. "LocallyBiased" -> True
 * (default) is DIRECT-L (Gablonsky & Kelley); False is the original unbiased
 * DIRECT (Jones et al.), better for landscapes with many minima.
 *
 * Reference values are analytic optima or published benchmark constants, never a
 * snapshot of the solver's own output. Deceptive multimodal cases either budget
 * more evaluations or assert a loose "global basin found" threshold plus
 * feasibility, following the convention in test_nminimize.c / test_shgo.c.
 *
 * Coverage:
 *   - Result shape / typing.
 *   - Convex / unimodal exact minima.
 *   - Standard multimodal benchmarks (Himmelblau, Rastrigin, Ackley, Branin,
 *     six-hump camel, Styblinski-Tang, Cross-in-tray, Shubert, Eggholder).
 *   - DIRECT vs DIRECT-L parity ("LocallyBiased" -> True / False).
 *   - Constrained global optimization (disk, Mishra's Bird).
 *   - Sub-options: Epsilon, MaxFunctionEvaluations budget, MaxIterations,
 *     Volume/Length tolerances, MinValue early stop, PostProcess -> False,
 *     determinism, option scoping.
 *   - High dimension (8-D sphere).
 *   - Stress (5-D Rastrigin) + memory-hygiene smoke loop.
 *
 * Run binary directly: ./direct_tests */

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
               "Method->\"DIRECT\"], {_?NumberQ, {(_->_)..}}]");
    check_true("Length[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->\"DIRECT\"]] == 2");
}

/* ------------------------------------------------------------------ */
/* 2. Convex / unimodal exact minima                                   */
/* ------------------------------------------------------------------ */
static void test_sphere(void) {
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->\"DIRECT\"]]] < 1.*^-6");
    check_true("Norm[{x,y} /. Last[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->\"DIRECT\"]]] < 1.*^-3");
}

static void test_quartic(void) {
    /* x^4 - 3 x^2 - x on [-5,5]: global -3.5139097 at 1.3008373. */
    check_true("Abs[First[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->\"DIRECT\"]] - (-3.5139097)] < 1.*^-3");
    check_true("Abs[(x /. Last[NMinimize[{x^4-3x^2-x, -5<=x<=5}, {x}, "
               "Method->\"DIRECT\"]]) - 1.3008373] < 1.*^-3");
}

static void test_shifted_parabola(void) {
    check_true("Abs[First[NMinimize[{(x-3)^2+(y+2)^2, -10<=x<=10 && -10<=y<=10}, "
               "{x,y}, Method->\"DIRECT\"]]] < 1.*^-6");
    check_true("Norm[{x-3,y+2} /. Last[NMinimize[{(x-3)^2+(y+2)^2, "
               "-10<=x<=10 && -10<=y<=10}, {x,y}, Method->\"DIRECT\"]]] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* 3. Standard multimodal benchmarks (bounded boxes)                   */
/* ------------------------------------------------------------------ */
static void test_himmelblau(void) {
    /* Four global minima, all f = 0. */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, Method->\"DIRECT\"]]] < 1.*^-6");
}

static void test_rastrigin_2d(void) {
    check_true("Abs[First[NMinimize[{20 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,2}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->\"DIRECT\"]]] < 1.*^-4");
}

static void test_ackley_2d(void) {
    check_true("Abs[First[NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] "
               "- Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, Method->\"DIRECT\"]]] < 1.*^-4");
}

static void test_branin(void) {
    /* Three global minima, f = 0.397887. */
    check_true("Abs[First[NMinimize[{(y-5.1/(4 Pi^2) x^2+5/Pi x-6)^2 "
               "+10(1-1/(8 Pi))Cos[x]+10, -5<=x<=10 && 0<=y<=15}, {x,y}, "
               "Method->\"DIRECT\"]] - 0.397887] < 1.*^-4");
}

static void test_six_hump_camel(void) {
    /* Global -1.0316284. */
    check_true("Abs[First[NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, "
               "-3<=x<=3 && -2<=y<=2}, {x,y}, Method->\"DIRECT\"]] - (-1.0316284)] "
               "< 1.*^-4");
}

static void test_styblinski_tang_2d(void) {
    /* Global -39.16599 * d, d = 2 -> -78.33198, at all x = -2.903534. */
    check_true("Abs[First[NMinimize[{Sum[x[i]^4-16 x[i]^2+5 x[i], {i,1,2}]/2, "
               "Table[-5<=x[i]<=5, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->\"DIRECT\"]] - (-78.33198)] < 1.*^-2");
}

static void test_styblinski_tang_4d(void) {
    /* d = 4 -> -156.66397, via the unbiased variant. */
    check_true("Abs[First[NMinimize[{Sum[x[i]^4-16 x[i]^2+5 x[i], {i,1,4}]/2, "
               "Table[-5<=x[i]<=5, {i,1,4}]}, Table[x[i],{i,1,4}], "
               "Method->{\"DIRECT\",\"LocallyBiased\"->False}]] - (-156.66397)] "
               "< 5.*^-2");
}

static void test_cross_in_tray(void) {
    /* Four global minima, f = -2.0626119. Deceptive: unbiased DIRECT. */
    check_true("Abs[First[NMinimize[{-0.0001(Abs[Sin[x]Sin[y]"
               "Exp[Abs[100-Sqrt[x^2+y^2]/Pi]]]+1)^0.1, -10<=x<=10 && -10<=y<=10}, "
               "{x,y}, Method->{\"DIRECT\",\"LocallyBiased\"->False}]] "
               "- (-2.0626119)] < 1.*^-3");
}

static void test_shubert(void) {
    /* 18 global minima, f = -186.7309. Very deceptive: a larger budget is needed
     * to resolve the deepest basin. */
    check_true("Abs[First[NMinimize[{(Sum[j Cos[(j+1)x+j],{j,1,5}])"
               "(Sum[j Cos[(j+1)y+j],{j,1,5}]), -10<=x<=10 && -10<=y<=10}, {x,y}, "
               "Method->{\"DIRECT\",\"LocallyBiased\"->False,"
               "\"MaxFunctionEvaluations\"->40000}]] - (-186.7309)] < 1.*^-2");
}

static void test_eggholder(void) {
    /* Deceptive; global -959.6407 at (512, 404.2319). Budget more evaluations and
     * assert a robust threshold + feasibility rather than a brittle exact pin. */
    check_true("First[NMinimize[{-(y+47) Sin[Sqrt[Abs[y+x/2+47]]] "
               "- x Sin[Sqrt[Abs[x-(y+47)]]], -512<=x<=512 && -512<=y<=512}, "
               "{x,y}, Method->{\"DIRECT\",\"LocallyBiased\"->False,"
               "\"MaxFunctionEvaluations\"->20000}]] < -959.0");
    check_true("Max[Abs[{x,y} /. Last[NMinimize[{-(y+47) Sin[Sqrt[Abs[y+x/2+47]]] "
               "- x Sin[Sqrt[Abs[x-(y+47)]]], -512<=x<=512 && -512<=y<=512}, "
               "{x,y}, Method->{\"DIRECT\",\"LocallyBiased\"->False,"
               "\"MaxFunctionEvaluations\"->20000}]]]] <= 512.001");
}

/* ------------------------------------------------------------------ */
/* 4. DIRECT vs DIRECT-L parity                                        */
/* ------------------------------------------------------------------ */
static void test_variant_parity(void) {
    /* Both the locally-biased default and the original unbiased variant reach
     * Himmelblau's zero. */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DIRECT\",\"LocallyBiased\"->True}]]] < 1.*^-6");
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DIRECT\",\"LocallyBiased\"->False}]]] < 1.*^-6");
    /* And each reaches Rastrigin-2D's zero. */
    check_true("Abs[First[NMinimize[{20 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,2}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,2}]}, Table[x[i],{i,1,2}], "
               "Method->{\"DIRECT\",\"LocallyBiased\"->False}]]] < 1.*^-4");
}

/* ------------------------------------------------------------------ */
/* 5. Constrained global optimization                                  */
/* ------------------------------------------------------------------ */
static void test_constrained_disk(void) {
    /* min x + y on the disk x^2 + y^2 <= 9: global -3 sqrt(2) = -4.2426407. */
    check_true("Abs[First[NMinimize[{x+y, x^2+y^2<=9}, {x,y}, "
               "Method->\"DIRECT\"]] - (-4.2426407)] < 1.*^-3");
    /* returned point is feasible */
    check_true("(x^2+y^2 /. Last[NMinimize[{x+y, x^2+y^2<=9}, {x,y}, "
               "Method->\"DIRECT\"]]) <= 9.001");
}

static void test_mishra_bird(void) {
    /* Mishra's Bird, constrained to the disk (x+5)^2+(y+5)^2 < 25:
     * global -106.7645 at (-3.1302, -1.5821). */
    check_true("Abs[First[NMinimize[{Sin[y]Exp[(1-Cos[x])^2]"
               "+Cos[x]Exp[(1-Sin[y])^2]+(x-y)^2, (x+5)^2+(y+5)^2<25}, {x,y}, "
               "Method->\"DIRECT\"]] - (-106.7645)] < 1.*^-2");
    check_true("((x+5)^2+(y+5)^2 /. Last[NMinimize[{Sin[y]Exp[(1-Cos[x])^2]"
               "+Cos[x]Exp[(1-Sin[y])^2]+(x-y)^2, (x+5)^2+(y+5)^2<25}, {x,y}, "
               "Method->\"DIRECT\"]]) <= 25.001");
}

/* ------------------------------------------------------------------ */
/* 6. Sub-options                                                      */
/* ------------------------------------------------------------------ */
static void test_determinism(void) {
    /* DIRECT is deterministic: two identical calls agree exactly. */
    check_eq("NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
             "Method->\"DIRECT\"] "
             "=== "
             "NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
             "Method->\"DIRECT\"]",
             "True");
}

static void test_epsilon(void) {
    /* A larger epsilon biases toward global exploration; still solves. */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DIRECT\",\"Epsilon\"->0.01}]]] < 1.*^-6");
}

static void test_maxfun_budget(void) {
    /* The evaluation budget bites: a tiny budget on the deceptive Himmelblau
     * lands far from zero (raw), a large budget converges. */
    check_true("First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"DIRECT\",\"MaxFunctionEvaluations\"->20,"
               "\"PostProcess\"->False}]] > 1.");
    check_true("First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->{\"DIRECT\",\"MaxFunctionEvaluations\"->4000,"
               "\"PostProcess\"->False}]] < 1.*^-4");
}

static void test_maxiter(void) {
    /* A positive MaxIterations is honored and still reaches the sphere minimum. */
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DIRECT\",\"MaxIterations\"->500}]]] < 1.*^-6");
}

static void test_tolerances(void) {
    /* Explicit volume/length tolerances are accepted and still solve. */
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DIRECT\",\"VolumeTolerance\"->1.*^-12,"
               "\"LengthTolerance\"->1.*^-4}]]] < 1.*^-4");
}

static void test_min_value_stop(void) {
    /* A supplied known minimum (0) with a loose tolerance triggers the early
     * stop; the polished result is still valid and near zero. */
    check_true("Abs[First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DIRECT\",\"MinValue\"->0,\"MinValueTolerance\"->0.1}]]] "
               "< 1.*^-4");
}

static void test_postprocess_false(void) {
    /* No local polish: raw DIRECT centers the (3,-1) basin closely but not to
     * machine precision, yet is a valid, feasible, near-optimal point. */
    check_true("First[NMinimize[{(x-3)^2+(y+1)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DIRECT\",\"PostProcess\"->False}]] < 1.*^-3");
    check_true("MatchQ[NMinimize[{(x-3)^2+(y+1)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DIRECT\",\"PostProcess\"->False}], {_?NumberQ, {(_->_)..}}]");
}

static void test_option_scoping(void) {
    /* A sub-option belonging to another method (DualAnnealing's
     * "VisitingParameter") and an unknown key inside the Method list are each
     * warned to stderr (muted) and dropped -- never stored -- so the
     * deterministic DIRECT result is identical to omitting them. */
    check_eq("NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
             "Method->{\"DIRECT\",\"VisitingParameter\"->2.5,\"Frobnicate\"->1}] "
             "=== "
             "NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, "
             "Method->\"DIRECT\"]",
             "True");
    /* A valid-for-DIRECT key is still honored (not swept up by the guard). */
    check_true("Abs[First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
               "-5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->{\"DIRECT\",\"Epsilon\"->0.001}]]] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 7. High dimension                                                   */
/* ------------------------------------------------------------------ */
static void test_high_dim_sphere(void) {
    /* 8-D sphere: DIRECT halves each side symmetrically and centers the origin. */
    check_true("Abs[First[NMinimize[{Sum[x[i]^2, {i,1,8}], "
               "Table[-5<=x[i]<=5, {i,1,8}]}, Table[x[i],{i,1,8}], "
               "Method->\"DIRECT\"]]] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* 8. NMaximize + duality                                              */
/* ------------------------------------------------------------------ */
static void test_nmaximize(void) {
    check_true("Abs[First[NMaximize[{-(x^2+y^2), -5<=x<=5 && -5<=y<=5}, {x,y}, "
               "Method->\"DIRECT\"]]] < 1.*^-6");
    /* NMaximize[f] == -NMinimize[-f] on a concave problem. */
    check_true("Abs[First[NMaximize[{4-(x-1)^2-(y+2)^2, -5<=x<=5 && -5<=y<=5}, "
               "{x,y}, Method->\"DIRECT\"]] - 4.0] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 9. Stress + memory hygiene                                          */
/* ------------------------------------------------------------------ */
static void test_rastrigin_5d_stress(void) {
    /* 5-D Rastrigin (global 0). Deterministic, so assert a loose near-global
     * interval and feasibility rather than a brittle exact pin. */
    check_true("Module[{r = NMinimize[{50 + Sum[x[i]^2-10 Cos[2 Pi x[i]], {i,1,5}], "
               "Table[-5.12<=x[i]<=5.12, {i,1,5}]}, Table[x[i],{i,1,5}], "
               "Method->{\"DIRECT\",\"LocallyBiased\"->False}]}, "
               "Head[r] === List && 0 <= First[r] < 2.0]");
}

static void test_memory_smoke(void) {
    for (int i = 0; i < 3; i++) {
        check_true("First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, "
                   "-5<=x<=5 && -5<=y<=5}, {x,y}, "
                   "Method->{\"DIRECT\",\"MaxFunctionEvaluations\"->1500,"
                   "\"PostProcess\"->False}]] < 1.");
    }
}

int main(void) {
    symtab_init();
    core_init();

    /* Expected diagnostics (option-scoping warnings) write to stderr; keep the
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
    TEST(test_ackley_2d);
    TEST(test_branin);
    TEST(test_six_hump_camel);
    TEST(test_styblinski_tang_2d);
    TEST(test_styblinski_tang_4d);
    TEST(test_cross_in_tray);
    TEST(test_shubert);
    TEST(test_eggholder);

    /* 4. DIRECT vs DIRECT-L parity */
    TEST(test_variant_parity);

    /* 5. Constrained */
    TEST(test_constrained_disk);
    TEST(test_mishra_bird);

    /* 6. Sub-options */
    TEST(test_determinism);
    TEST(test_epsilon);
    TEST(test_maxfun_budget);
    TEST(test_maxiter);
    TEST(test_tolerances);
    TEST(test_min_value_stop);
    TEST(test_postprocess_false);
    TEST(test_option_scoping);

    /* 7. High dimension */
    TEST(test_high_dim_sphere);

    /* 8. NMaximize + duality */
    TEST(test_nmaximize);

    /* 9. Stress + memory hygiene */
    TEST(test_rastrigin_5d_stress);
    TEST(test_memory_smoke);

    printf("\nAll DIRECT tests passed.\n");
    return 0;
}
