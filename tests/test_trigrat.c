/*
 * test_trigrat.c -- Tests for the algebraic-normal-form trig rational
 * fast path in src/trigrat.c.
 *
 * Coverage:
 *   - Pythagorean reductions (trig and hyperbolic).
 *   - Tan/Sec/Csc/Cot rationalisation.
 *   - Multi-argument cases (within budget).
 *   - The two user-reported derivative-of-Risch-Norman cases as
 *     correctness + performance smoke tests.
 *   - Non-applicable inputs: the leaf-count gate must not let the
 *     fast path regress short, identity-recognisable inputs.
 */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <time.h>

#define TEST(func) do { printf("Running test: %s\n", #func); func(); } while(0)

/* ---- Pythagorean reductions ---- */

void test_trigrat_pythagorean(void) {
    assert_eval_eq("Simplify[Sin[x]^2 + Cos[x]^2]", "1", 0);
    assert_eval_eq("Simplify[(Sin[x]^2 + Cos[x]^2)^5]", "1", 0);
}

void test_trigrat_hyperbolic_pythagorean(void) {
    assert_eval_eq("Simplify[Cosh[x]^2 - Sinh[x]^2]", "1", 0);
}

void test_trigrat_sec_minus_tan(void) {
    assert_eval_eq("Simplify[Sec[x]^2 - Tan[x]^2]", "1", 0);
}

/* ---- Tan/Cot/Sec/Csc rationalisation ---- */

void test_trigrat_tan_cos_product(void) {
    assert_eval_eq("Simplify[Tan[x] Cos[x]]", "Sin[x]", 0);
}

void test_trigrat_sin_over_cos_equals_tan(void) {
    assert_eval_eq("Simplify[Sin[x]/Cos[x] - Tan[x]]", "0", 0);
}

/* ---- Multi-argument (within budget) ---- */

void test_trigrat_angle_addition_sin(void) {
    /* This is recognised by simp_search's TanAddition / TrigReduce; the
     * trig-rational fast path does not see it (under the leaf-count
     * floor), but the result must remain correct. */
    assert_eval_eq("Simplify[Sin[a] Cos[b] + Cos[a] Sin[b]]",
                   "Sin[a + b]", 0);
}

/* ---- Large user-reported cases (correctness + performance) ---- */

/* Wall-clock helper. */
static double seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

void test_trigrat_risch_tan_squared_derivative(void) {
    /* D[Integrate`RischNorman[Tan[x]^2 + Tan[x] + 1, x], x] // Simplify
     * should reduce to the integrand. Without the trig-rational fast
     * path, Simplify hangs indefinitely. We allow <= 30 s wall-clock
     * and accept any equivalent compact form (Sec[x]^2 + Tan[x], or the
     * factored Sec[x]^2 (1 + Cos[x] Sin[x]), etc.). */
    struct Expr* parsed = parse_expression(
        "Simplify[D[Integrate`RischNorman[Tan[x]^2 + Tan[x] + 1, x], x]]");
    assert(parsed != NULL);
    double t0 = seconds();
    struct Expr* r = evaluate(parsed);
    double dt = seconds() - t0;
    expr_free(parsed);
    char* s = expr_to_string(r);
    /* The result should be small (under ~20 leaves) and not contain any
     * leftover Log[...] -- if the Log terms didn't cancel, the
     * algorithm failed. */
    int has_log = (strstr(s, "Log[") != NULL);
    if (has_log) {
        fprintf(stderr, "FAIL: trigrat_risch_tan_squared: Log not cancelled. Got: %s\n", s);
    }
    if (dt > 30.0) {
        fprintf(stderr, "FAIL: trigrat_risch_tan_squared: took %.2f s (> 30 s budget)\n", dt);
    }
    assert(!has_log);
    assert(dt <= 30.0);
    free(s);
    expr_free(r);
}

void test_trigrat_risch_xn_sin_cos_derivative(void) {
    /* D[Integrate`RischNorman[x^4 Sin[x] Cos[x], x], x] // Simplify
     * should give back x^4 Sin[x] Cos[x] (or an equivalent compact
     * polynomial-times-trig product). Performance budget: < 5 s. */
    struct Expr* parsed = parse_expression(
        "Simplify[D[Integrate`RischNorman[x^4 Sin[x] Cos[x], x], x]]");
    assert(parsed != NULL);
    double t0 = seconds();
    struct Expr* r = evaluate(parsed);
    double dt = seconds() - t0;
    expr_free(parsed);
    char* s = expr_to_string(r);
    /* The simplified form contains x^4, Sin, Cos and is short. We
     * assert the result is *not* large (under 50 leaves printed). */
    size_t L = strlen(s);
    if (L > 200) {
        fprintf(stderr, "FAIL: trigrat_risch_xn_sin_cos: result too large (len=%zu). Got: %s\n",
                L, s);
    }
    if (dt > 5.0) {
        fprintf(stderr, "FAIL: trigrat_risch_xn_sin_cos: took %.2f s (> 5 s budget)\n", dt);
    }
    assert(L <= 200);
    assert(dt <= 5.0);
    free(s);
    expr_free(r);
}

/* ---- Non-applicable inputs: must not regress simp_search wins ---- */

void test_trigrat_double_angle_preserved(void) {
    /* 2 Tan[x]/(1+Tan[x]^2) -> Sin[2x]. The fast path does NOT fire
     * (below leaf floor). simp_search must still produce Sin[2x]. */
    assert_eval_eq("Simplify[2 Tan[x]/(1+Tan[x]^2)]", "Sin[2 x]", 0);
}

void test_trigrat_half_angle_preserved(void) {
    /* Sin[x]/(2 (Cos[x]+1)) -> Tan[x/2]/2. Below leaf floor; relies
     * on simp_search's HalfAngle transform. */
    assert_eval_eq("Simplify[Sin[x]/(2 (Cos[x] + 1))]",
                   "1/2 Tan[1/2 x]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_trigrat_pythagorean);
    TEST(test_trigrat_hyperbolic_pythagorean);
    TEST(test_trigrat_sec_minus_tan);
    TEST(test_trigrat_tan_cos_product);
    TEST(test_trigrat_sin_over_cos_equals_tan);
    TEST(test_trigrat_angle_addition_sin);
    TEST(test_trigrat_risch_tan_squared_derivative);
    TEST(test_trigrat_risch_xn_sin_cos_derivative);
    TEST(test_trigrat_double_angle_preserved);
    TEST(test_trigrat_half_angle_preserved);

    printf("All trigrat tests passed!\n");
    return 0;
}
