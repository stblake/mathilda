/* test_risch_coupled.c — coupled differential systems + the hypertangent
 * reduced case (Bronstein §5.10, Chapter 8).
 *
 * Verifies:
 *   - Risch`CoupledDESystem (§8.1): the coupled 2x2 real system solved via the
 *     reduction to a single Risch DE over C(i)(x).  Book worked values
 *     (Examples 5.10.3, 8.4.1) plus defining-equation round-trips.
 *   - Risch`IntegrateHypertangentReduced (§5.10, p.169): peeling t^2+1 poles
 *     one multiplicity at a time; the elementary round-trip p - D[q] in k[t],
 *     and the non-elementary certificate (Example 5.10.2, ∫ sin x / x).
 *
 * The base field is C(x), D = d/dx.  A hypertangent monomial t = tan(x) has
 * Dt = 1 + t^2; the half-angle t = tan(x/2) has Dt = (1 + t^2)/2.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

static void run_zero(const char* diff) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Expand[Together[%s]]", diff);
    run_test(buf, "0");
}

/* t = tan(x): Dt = 1 + t^2. */
#define HT "{x -> 1, t -> 1 + t^2}"
#define CDS "Risch`CoupledDESystem"
#define IHR "Risch`IntegrateHypertangentReduced"

/* ---- CoupledDESystem: book worked values (§5.10, §8.4) ---------------- */
static void test_coupled_book_values(void) {
    /* Example 5.10.3 reduced sub-solve: (c, d) = (x/18 + 1/6, 1/108 - x^2/6). */
    run_zero(CDS "[0, 6, x^2, 1, x][[1]] - (x/18 + 1/6)");
    run_zero(CDS "[0, 6, x^2, 1, x][[2]] - (1/108 - x^2/6)");
    /* Example 8.4.1 step (s1, s2) = (-1, 2x + 1). */
    run_zero(CDS "[0, 4 x - 2, 2 - 8 x^2, 4 - 4 x, x][[1]] + 1");
    run_zero(CDS "[0, 4 x - 2, 2 - 8 x^2, 4 - 4 x, x][[2]] - (2 x + 1)");
    /* Example 8.4.1 inner recursion (s1, s2) = (1, 0). */
    run_zero(CDS "[0, 4 x, 0, 4 x, x][[1]] - 1");
    run_test(CDS "[0, 4 x, 0, 4 x, x][[2]]", "0");
    /* Homogeneous RHS -> trivial solution. */
    run_test(CDS "[x, 1/x, 0, 0, x]", "List[0, 0]");
}

/* ---- CoupledDESystem: defining-equation round-trips ------------------- */
/* {y1, y2} must satisfy  D[y1] + f1 y1 - f2 y2 = g1  and
 *                        D[y2] + f2 y1 + f1 y2 = g2  over C(x). */
static void run_coupled_rt(const char* f1, const char* f2,
                           const char* g1, const char* g2) {
    char buf[4096];
    /* Solvable by assumption: assert a solution WAS returned, and that it
     * satisfies both defining equations. */
    snprintf(buf, sizeof buf,
        "With[{s = " CDS "[%s, %s, %s, %s, x]}, "
        "Head[s] =!= " CDS " && "
        "Expand[Together[D[s[[1]], x] + (%s) s[[1]] - (%s) s[[2]] - (%s)]] === 0 && "
        "Expand[Together[D[s[[2]], x] + (%s) s[[1]] + (%s) s[[2]] - (%s)]] === 0]",
        f1, f2, g1, g2, f1, f2, g1, f2, f1, g2);
    run_test(buf, "True");
}

static void test_coupled_roundtrips(void) {
    run_coupled_rt("0", "6", "x^2", "1");                 /* polynomial solution */
    run_coupled_rt("0", "4 x - 2", "2 - 8 x^2", "4 - 4 x");
    run_coupled_rt("0", "1", "1", "0");                   /* constant coefficients */
    run_coupled_rt("1/x", "0", "1", "x");                 /* primitive-type f1 */
    run_coupled_rt("0", "2", "x", "x^2");
}

/* ---- CoupledDESystem stress: constructed-solution round-trips --------- */
/* Given ANY chosen (y1, y2, f1, f2), build g1 = D[y1] + f1 y1 - f2 y2 and
 * g2 = D[y2] + f2 y1 + f1 y2, so the system is solvable by construction; then
 * the solver must return a {s1, s2} that satisfies BOTH defining equations
 * (not necessarily y1, y2 themselves — homogeneous solutions are allowed). */
static void run_coupled_constructed(const char* y1, const char* y2,
                                    const char* f1, const char* f2) {
    char buf[6000];
    snprintf(buf, sizeof buf,
        "With[{g1 = D[%s, x] + (%s)(%s) - (%s)(%s), g2 = D[%s, x] + (%s)(%s) + (%s)(%s)}, "
        "With[{s = " CDS "[%s, %s, g1, g2, x]}, "
        "Head[s] =!= " CDS " && "
        "Expand[Together[D[s[[1]], x] + (%s) s[[1]] - (%s) s[[2]] - g1]] === 0 && "
        "Expand[Together[D[s[[2]], x] + (%s) s[[1]] + (%s) s[[2]] - g2]] === 0]]",
        y1, f1, y1, f2, y2,  y2, f2, y1, f1, y2,
        f1, f2,  f1, f2,  f2, f1);
    run_test(buf, "True");
}

static void test_coupled_constructed_stress(void) {
    /* Polynomial solutions, varied degree, f2 constant and non-constant. */
    run_coupled_constructed("x", "x^2", "0", "1");
    run_coupled_constructed("x^3 - x", "x^2 + 1", "0", "3");
    run_coupled_constructed("x^4", "x", "0", "2 x");
    run_coupled_constructed("x^2 + 2 x + 1", "x^3", "0", "x^2");
    /* f1 != 0 (the general Chapter-8 system, primitive-type real part). */
    run_coupled_constructed("x^3", "1", "1", "x");
    run_coupled_constructed("x^2 + 1", "x - 1", "x", "1");
    run_coupled_constructed("x", "x^2 + x", "2", "3");
    /* Rational solutions (poles off the real axis and on it). */
    run_coupled_constructed("1/(x^2 + 1)", "x", "0", "2");
    run_coupled_constructed("1/x", "1/x^2", "0", "x");
    run_coupled_constructed("(x - 1)/(x + 2)", "1/(x^2 + 4)", "0", "1");
    run_coupled_constructed("x/(x^2 + 1)", "1/(x^2 + 1)", "0", "4");
    /* Both coefficients non-trivial, rational data. */
    run_coupled_constructed("x^2", "1/(x + 1)", "1/x", "0");
    run_coupled_constructed("x + 1", "x^2 - 1", "0", "x + 1");
}

/* ---- CoupledDESystem: no-solution + robustness ------------------------ */
static void test_coupled_no_solution(void) {
    /* Example 5.10.2 kernel: y' + i y = 2/x has no rational solution
     * (this is exactly why ∫ sin x / x is not elementary — the ExpIntegral
     * obstruction).  Simple-pole right-hand sides c/x are unsolvable. */
    run_test("Head[" CDS "[0, 1, 2/x, 0, x]] === " CDS, "True");
    run_test("Head[" CDS "[0, 1, 1/x, 0, x]] === " CDS, "True");
    run_test("Head[" CDS "[0, 1, 0, 1/x, x]] === " CDS, "True");
    run_test("Head[" CDS "[0, 1, 1/x^2, 0, x]] === " CDS, "True");
    /* Malformed: non-symbol variable / wrong arity -> unevaluated. */
    run_test("Head[" CDS "[0, 1, 1, 0, 5]] === " CDS, "True");
    run_test("Head[" CDS "[0, 1, 1, 0]] === " CDS, "True");
}

/* ---- IntegrateHypertangentReduced: elementary round-trips (§5.10) ----- */
/* q removes all t^2+1 poles: p - D[q] must be a polynomial in t (denominator
 * free of t), and beta must be True. */
static void run_reduced_elementary(const char* p, const char* deriv) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "With[{r = " IHR "[%s, t, %s]}, "
        "r[[2]] && FreeQ[Denominator[Together[(%s) - Risch`Derivation[r[[1]], %s]]], t]]",
        p, deriv, p, deriv);
    run_test(buf, "True");
}

static void test_reduced_elementary(void) {
    run_reduced_elementary("1/(t^2 + 1)", HT);
    run_reduced_elementary("1/(t^2 + 1)^2", HT);
    run_reduced_elementary("(2 t + 3)/(t^2 + 1)", HT);
    run_reduced_elementary("t/(t^2 + 1)^2", HT);
    /* Example 5.10.3: (t^3 + t^2 + x^2 t + 1)/(t^2 + 1)^2. */
    run_reduced_elementary("(t^3 + t^2 + x^2 t + 1)/(t^2 + 1)^2", HT);
    /* Deeper pole. */
    run_reduced_elementary("(t + x)/(t^2 + 1)^3", HT);
    /* Half-angle monomial t = tan(x/2), Dt = (1 + t^2)/2. */
    run_reduced_elementary("t/(t^2 + 1)", "{x -> 1, t -> (1 + t^2)/2}");

    /* No pole at t^2+1 (already a polynomial): q = 0, beta = True. */
    run_test(IHR "[t^2 + x t + 1, t, " HT "]", "List[0, True]");
    run_test(IHR "[3 t + x, t, " HT "]", "List[0, True]");
    run_test(IHR "[0, t, " HT "]", "List[0, True]");
}

/* ---- IntegrateHypertangentReduced stress: high multiplicity + scaling -- */
/* t = tan(2x): Dt = 2(1+t^2);  t = tan(3x): Dt = 3(1+t^2). */
#define HT2 "{x -> 1, t -> 2 (1 + t^2)}"
#define HT3 "{x -> 1, t -> 3 (1 + t^2)}"
static void test_reduced_stress(void) {
    /* Deep poles at t^2+1 (a=1). */
    run_reduced_elementary("t/(t^2 + 1)^3", HT);
    run_reduced_elementary("(t + x)/(t^2 + 1)^4", HT);
    run_reduced_elementary("(t^2 + x t + 1)/(t^2 + 1)^5", HT);
    run_reduced_elementary("(x^3 t + x)/(t^2 + 1)^4", HT);
    run_reduced_elementary("(t^3 + t)/(t^2 + 1)^3", HT);
    /* Scaled monomials tan(2x), tan(3x). */
    run_reduced_elementary("t/(t^2 + 1)^2", HT2);
    run_reduced_elementary("1/(t^2 + 1)^2", HT2);
    run_reduced_elementary("(x^2 t + x)/(t^2 + 1)^3", HT3);
    run_reduced_elementary("(t + 1)/(t^2 + 1)^2", HT2);
    /* Half-angle monomial tan(x/2), higher multiplicity. */
    run_reduced_elementary("t/(t^2 + 1)^2", "{x -> 1, t -> (1 + t^2)/2}");
    run_reduced_elementary("(t + x)/(t^2 + 1)^3", "{x -> 1, t -> (1 + t^2)/2}");
    /* Concrete value: tan(2x), ∫ over t/(t^2+1)^2 -> q = -1/(8(t^2+1)^2). */
    run_zero(IHR "[t/(t^2 + 1)^2, t, " HT2 "][[1]] + 1/(8 (t^2 + 1)^2)");
}

/* ---- IntegrateHypertangentReduced: non-elementarity certificate ------- */
static void test_reduced_non_elementary(void) {
    /* Example 5.10.2: ∫ sin x / x, t = tan(x/2), p = 2t/(x(t^2+1)).  The pole
     * peel needs CoupledDESystem[0,1,2/x,0] which has no solution -> beta False. */
    run_test(IHR "[2 t/(x (t^2 + 1)), t, {x -> 1, t -> (1 + t^2)/2}]", "List[0, False]");
}

/* ---- Robustness: non-hypertangent monomial + malformed input ---------- */
static void test_reduced_robustness(void) {
    /* Exponential monomial (Dt = t) is not hypertangent -> unevaluated. */
    run_test("Head[" IHR "[1/(t^2 + 1), t, {x -> 1, t -> t}]] === " IHR, "True");
    /* Log monomial (Dt = 1/x) is not hypertangent -> unevaluated. */
    run_test("Head[" IHR "[1/(t^2 + 1), t, {x -> 1, t -> 1/x}]] === " IHR, "True");
    /* Wrong arity / malformed derivation. */
    run_test("Head[" IHR "[1/(t^2 + 1), t]] === " IHR, "True");
    run_test("Head[" IHR "[1/(t^2 + 1), t, {x, t}]] === " IHR, "True");
}

/* ---- CoupledDECancelTan (Bronstein §8.4): book value + diff-back ------ */
/* t = tan(2x): Dt = 2(1+t^2), eta = 2;  t = tan(3x): Dt = 3(1+t^2), eta = 3. */
#define HT2 "{x -> 1, t -> 2 (1 + t^2)}"
#define HT3 "{x -> 1, t -> 3 (1 + t^2)}"
#define HALF "{x -> 1, t -> (1 + t^2)/2}"
#define CCT "Risch`CoupledDECancelTan"

static void test_cancel_tan_book(void) {
    /* Example 8.4.1: CoupledDECancelTan[0, 4x, -t^2+2t-8x^2+1, 2(1-2x), t, D, 2]
     * has solution (q1, q2) = (t - 1, 2 x). */
    run_zero(CCT "[0, 4 x, -t^2 + 2 t - 8 x^2 + 1, 2 (1 - 2 x), t, " HT ", 2][[1]] - (t - 1)");
    run_zero(CCT "[0, 4 x, -t^2 + 2 t - 8 x^2 + 1, 2 (1 - 2 x), t, " HT ", 2][[2]] - 2 x");
    /* Inner recursion CoupledDECancelTan(0, 4x+1, -t, 4x+1, D, 1) = (1, 0). */
    run_test(CCT "[0, 4 x + 1, -t, 4 x + 1, t, " HT ", 1][[1]]", "1");
    run_test(CCT "[0, 4 x + 1, -t, 4 x + 1, t, " HT ", 1][[2]]", "0");
    /* Deepest base call CoupledDECancelTan(0, 4x+2, 0, 0, D, 0) = (0, 0). */
    run_test(CCT "[0, 4 x + 2, 0, 0, t, " HT ", 0]", "List[0, 0]");
}

/* ---- CoupledDECancelTan stress: constructed-solution diff-back -------- */
/* Given a chosen (q1, q2) in k[t] of degree <= n and b0, b2 in k, build the RHS
 * of the SYMMETRIC tangent-cancellation system (diagonal b0 - n eta t, forced by
 * the complex-scalar structure of Bronstein eq. 8.14):
 *     c1 = D[q1] + (b0 - n eta t) q1 - b2 q2
 *     c2 = D[q2] + b2 q1 + (b0 - n eta t) q2
 * so the system is solvable by construction; the solver must return SOME
 * {s1, s2} satisfying BOTH equations (homogeneous solutions allowed). */
static void run_cancel_tan_ct(const char* q1, const char* q2,
                              const char* b0, const char* b2,
                              int n, const char* deriv, const char* eta) {
    char buf[8192];
    snprintf(buf, sizeof buf,
      "With[{c1 = Risch`Derivation[%s, %s] + ((%s) - %d (%s) t)(%s) - (%s)(%s), "
            "c2 = Risch`Derivation[%s, %s] + (%s)(%s) + ((%s) - %d (%s) t)(%s)}, "
      "With[{s = " CCT "[%s, %s, c1, c2, t, %s, %d]}, "
      "Head[s] =!= " CCT " && "
      "Expand[Together[Risch`Derivation[s[[1]], %s] + ((%s) - %d (%s) t) s[[1]] - (%s) s[[2]] - c1]] === 0 && "
      "Expand[Together[Risch`Derivation[s[[2]], %s] + (%s) s[[1]] + ((%s) - %d (%s) t) s[[2]] - c2]] === 0]]",
      q1, deriv, b0, n, eta, q1, b2, q2,
      q2, deriv, b2, q1, b0, n, eta, q2,
      b0, b2, deriv, n,
      deriv, b0, n, eta, b2,
      deriv, b2, b0, n, eta);
    run_test(buf, "True");
}

static void test_cancel_tan_constructed(void) {
    /* Example 8.4.1 reconstructed from its own solution. */
    run_cancel_tan_ct("t - 1", "2 x", "0", "4 x", 2, HT, "1");
    /* n = 0 base case: (q1, q2) in k = C(x). */
    run_cancel_tan_ct("x", "1", "0", "1", 0, HT, "1");
    run_cancel_tan_ct("x^2", "x", "0", "2", 0, HT, "1");
    run_cancel_tan_ct("x", "x^2", "0", "x", 0, HT, "1");
    /* n = 1, 2, 3 over tan(x). */
    run_cancel_tan_ct("t + x", "1", "0", "x", 1, HT, "1");
    run_cancel_tan_ct("x t", "t - 1", "0", "2", 1, HT, "1");
    run_cancel_tan_ct("t^2 + x", "t", "0", "1", 2, HT, "1");
    run_cancel_tan_ct("t^2 - t", "t^2 + 1", "0", "x^2", 2, HT, "1");
    run_cancel_tan_ct("t^3", "x t", "0", "x", 3, HT, "1");
    /* Nonzero b0 (primitive-type real part, constant in k). */
    run_cancel_tan_ct("t", "1", "1", "1", 1, HT, "1");
    run_cancel_tan_ct("t^2 + t", "x", "2", "x", 2, HT, "1");
}

static void test_cancel_tan_stress(void) {
    /* Scaled monomial tan(2x): eta = 2. */
    run_cancel_tan_ct("t - 1", "x", "0", "1", 1, HT2, "2");
    run_cancel_tan_ct("t^2", "t", "0", "2 x", 2, HT2, "2");
    run_cancel_tan_ct("x t + 1", "t", "0", "x", 2, HT2, "2");
    /* tan(3x): eta = 3. */
    run_cancel_tan_ct("t^2 + x t", "t - x", "0", "x^2", 2, HT3, "3");
    run_cancel_tan_ct("t^3 + t", "x", "0", "1", 3, HT3, "3");
    /* Half-angle tan(x/2): eta = 1/2. */
    run_cancel_tan_ct("t + 1", "1", "0", "1", 1, HALF, "1/2");
    run_cancel_tan_ct("t^2 + x", "t", "0", "x", 2, HALF, "1/2");
    /* Deeper degree bounds. */
    run_cancel_tan_ct("t^4", "t^2 + 1", "0", "x", 4, HT, "1");
    run_cancel_tan_ct("x^2 t^3 + t", "t^2 - x", "0", "2 x", 3, HT, "1");
}

/* ---- CoupledDECancelTan: no-solution + robustness -------------------- */
static void test_cancel_tan_no_solution(void) {
    /* Base CoupledDESystem(0,1,2/x,0) is the Example 5.10.2 kernel (∫ sin x/x
     * obstruction) — no rational solution, so CoupledDECancelTan declines. */
    run_test("Head[" CCT "[0, 1, 2/x, 0, t, " HT ", 0]] === " CCT, "True");
    run_test("Head[" CCT "[0, 1, 1/x, 0, t, " HT ", 0]] === " CCT, "True");
    /* n = 1: the simple-pole obstruction survives the peel to a decoupled
     * integration D s1 = 2/x with no rational primitive. */
    run_test("Head[" CCT "[0, 1, 2/x, 0, t, " HT ", 1]] === " CCT, "True");
}

static void test_cancel_tan_robustness(void) {
    /* Wrong arity. */
    run_test("Head[" CCT "[0, 1, 0, 0, t, " HT "]] === " CCT, "True");
    /* Non-symbol monomial. */
    run_test("Head[" CCT "[0, 1, 0, 0, 5, " HT ", 1]] === " CCT, "True");
    /* Negative or non-integer degree bound. */
    run_test("Head[" CCT "[0, 1, 0, 0, t, " HT ", -1]] === " CCT, "True");
    run_test("Head[" CCT "[0, 1, 0, 0, t, " HT ", x]] === " CCT, "True");
    /* Non-hypertangent monomial: exponential (Dt = t) / logarithmic (Dt = 1/x). */
    run_test("Head[" CCT "[0, 1, 0, 0, t, {x -> 1, t -> t}, 1]] === " CCT, "True");
    run_test("Head[" CCT "[0, 1, 0, 0, t, {x -> 1, t -> 1/x}, 1]] === " CCT, "True");
}

int main(void) {
    core_init();

    TEST(test_coupled_book_values);
    TEST(test_coupled_roundtrips);
    TEST(test_coupled_constructed_stress);
    TEST(test_coupled_no_solution);
    TEST(test_reduced_elementary);
    TEST(test_reduced_stress);
    TEST(test_reduced_non_elementary);
    TEST(test_reduced_robustness);
    TEST(test_cancel_tan_book);
    TEST(test_cancel_tan_constructed);
    TEST(test_cancel_tan_stress);
    TEST(test_cancel_tan_no_solution);
    TEST(test_cancel_tan_robustness);

    printf("All risch_coupled tests passed.\n");
    return 0;
}
