#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/* Local strong-assert wrapper. The shared test_utils.h `assert_eval_eq`
 * uses assert(), which the cmake build silences via NDEBUG; so a
 * mismatch only prints FAIL and the run exits 0. We need a hard signal,
 * so this file routes through a wrapper that calls exit(1) on
 * mismatch. */
static int g_failures = 0;
static void check_eval_eq(const char* input, const char* expected) {
    struct Expr* parsed = parse_expression(input);
    if (!parsed) {
        fprintf(stderr, "FAIL: parse failure for: %s\n", input);
        g_failures++;
        return;
    }
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    if (strcmp(str, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n",
                input, expected, str);
        g_failures++;
    }
    free(str);
    expr_free(evaluated);
}
/* Override the shared macro so existing call sites benefit from the
 * stronger check without changing other test files. */
#define assert_eval_eq(in, exp, ff) check_eval_eq((in), (exp))

/*
 * Tests for algorithmic radical simplification.
 *
 * Three phases of work, each toggled by a #define below. The full set
 * of 10 user-supplied cases lives in this file from day one: cases not
 * yet implemented are stubbed via SKIP_PHASE so the suite stays green.
 *
 *   Phase 1 — Sqrt-of-Sqrt denesting via the half-sum identity
 *               sqrt(A + sqrt(B)) = sqrt((A+s)/2) + sqrt((A-s)/2),
 *             where s^2 = A^2 - B simplifies to a clean square.
 *             Closes cases 1, 2, 3, 6, 7.
 *
 *   Phase 2 — Multi-generator algebraic-number reduction; rationalises
 *             denominators containing radicals via extended Euclidean
 *             in Q[x_1,...,x_k]/(m_1,...,m_k). Closes cases 8, 10.
 *
 *   Phase 3 — Cube-root denesting and sum-of-conjugate-cube-roots
 *             pattern recognisers. Closes cases 4, 5.
 *
 * Case 9 (zero-test for 1/(1+Sqrt[2]+Sqrt[3]) - (2+Sqrt[2]-Sqrt[6])/4)
 * is already handled by the existing simp_algebraic seed and therefore
 * runs unconditionally as a regression sentinel.
 */

#define PHASE_1
#define PHASE_2
#define PHASE_3

#define SKIP_PHASE(reason) do { \
    printf("  [SKIP] %s (phase not enabled): %s\n", __func__, reason); \
} while (0)

/* ---- Case 9 — regression sentinel (already passes) ---- */

void test_radical_case9_zero_test_already_passes(void) {
    /* 1/(1 + Sqrt[2] + Sqrt[3]) - (2 + Sqrt[2] - Sqrt[6])/4 == 0
     * Existing simp_algebraic handles the single composite extension
     * Q(sqrt(6) = sqrt(2)*sqrt(3)). This test guards against a phase-2
     * regression in that pipeline. */
    assert_eval_eq(
        "Simplify[1/(1 + Sqrt[2] + Sqrt[3]) - (2 + Sqrt[2] - Sqrt[6])/4]",
        "0", 0);
}

/* ---- Phase 1 — Square-root denesting ---- */

void test_radical_case1_sqrt_3_plus_2sqrt2(void) {
    /* sqrt(3 + 2 sqrt(2)) = 1 + sqrt(2). */
#ifdef PHASE_1
    assert_eval_eq("Simplify[Sqrt[3 + 2 Sqrt[2]] - (1 + Sqrt[2])]", "0", 0);
#else
    SKIP_PHASE("Phase 1: sqrt-of-sqrt denesting");
#endif
}

void test_radical_case2_sqrt_17_minus_12sqrt2(void) {
    /* sqrt(17 - 12 sqrt(2)) = 3 - 2 sqrt(2).
     * Note the sign: 3 - 2*sqrt(2) ~= 3 - 2.83 > 0, so this is the
     * principal-branch square root. (1 - sqrt(2)) - sqrt(...) should
     * NOT come out as the simplified form. */
#ifdef PHASE_1
    assert_eval_eq("Simplify[Sqrt[17 - 12 Sqrt[2]] - (3 - 2 Sqrt[2])]", "0", 0);
#else
    SKIP_PHASE("Phase 1: sqrt-of-sqrt denesting");
#endif
}

void test_radical_case3_sqrt_2_plus_sqrt3(void) {
    /* sqrt(2 + sqrt(3)) = (sqrt(6) + sqrt(2))/2.
     * The half-sum identity gives sqrt(3/2) + sqrt(1/2), which the
     * existing radical canonicaliser must collapse to (sqrt(6)+sqrt(2))/2. */
#ifdef PHASE_1
    assert_eval_eq("Simplify[Sqrt[2 + Sqrt[3]] - (Sqrt[6] + Sqrt[2])/2]", "0", 0);
#else
    SKIP_PHASE("Phase 1: sqrt-of-sqrt denesting");
#endif
}

void test_radical_case6_symbolic_two_term(void) {
    /* sqrt(x + y + 2 sqrt(x y)) = sqrt(x) + sqrt(y) for x, y > 0.
     * Discriminant: A^2 - B = (x+y)^2 - 4xy = (x-y)^2. */
#ifdef PHASE_1
    assert_eval_eq(
        "Assuming[x > 0 && y > 0, Simplify[Sqrt[x + y + 2 Sqrt[x y]] - (Sqrt[x] + Sqrt[y])]]",
        "0", 0);
#else
    SKIP_PHASE("Phase 1: sqrt-of-sqrt denesting");
#endif
}

void test_radical_case7_symbolic_half_angle(void) {
    /* sqrt(x + sqrt(x^2 - y^2)) = sqrt((x+y)/2) + sqrt((x-y)/2) for x>y>0.
     * Discriminant: A^2 - B = x^2 - (x^2 - y^2) = y^2. */
#ifdef PHASE_1
    assert_eval_eq(
        "Assuming[x > y && y > 0, Simplify[Sqrt[x + Sqrt[x^2 - y^2]] - (Sqrt[(x+y)/2] + Sqrt[(x-y)/2])]]",
        "0", 0);
#else
    SKIP_PHASE("Phase 1: sqrt-of-sqrt denesting");
#endif
}

/* ---- Phase 1 — branch-soundness subtests ---- */

void test_radical_phase1_no_assumption_no_denest_symbolic(void) {
    /* Without an Assuming, we cannot prove that (A+s)/2 and (A-s)/2 are
     * both nonneg for symbolic x, y, so the denesting must NOT fire on
     * Sqrt[x + Sqrt[y]]: the input must be returned unchanged
     * (modulo trivial canonicalisation). */
#ifdef PHASE_1
    assert_eval_eq(
        "Simplify[Sqrt[x + Sqrt[y]]]",
        "Sqrt[x + Sqrt[y]]", 0);
#else
    SKIP_PHASE("Phase 1: branch-soundness scaffold");
#endif
}

void test_radical_phase1_negative_discriminant_no_fire(void) {
    /* Sqrt[1 + Sqrt[2]]: A=1, B=2, D=A^2-B=-1. Negative discriminant
     * means the algorithm cannot produce a real denesting; it must
     * leave the expression alone. */
#ifdef PHASE_1
    assert_eval_eq(
        "Simplify[Sqrt[1 + Sqrt[2]]]",
        "Sqrt[1 + Sqrt[2]]", 0);
#else
    SKIP_PHASE("Phase 1: branch-soundness scaffold");
#endif
}

void test_radical_phase1_denest_minus_branch(void) {
    /* sqrt(3 - 2 sqrt(2)) = sqrt(2) - 1 (NOT 1 - sqrt(2), which is
     * negative). The principal-branch sqrt of a positive real is
     * positive. */
#ifdef PHASE_1
    assert_eval_eq("Simplify[Sqrt[3 - 2 Sqrt[2]] - (Sqrt[2] - 1)]", "0", 0);
#else
    SKIP_PHASE("Phase 1: branch-soundness scaffold");
#endif
}

void test_radical_phase1_idempotent_already_denested(void) {
    /* The denested form must be a fixed point under further Simplify. */
#ifdef PHASE_1
    assert_eval_eq("Simplify[1 + Sqrt[2]]", "1 + Sqrt[2]", 0);
#else
    SKIP_PHASE("Phase 1: branch-soundness scaffold");
#endif
}

/* ---- Phase 2 — Algebraic-number rationalisation ---- */

void test_radical_case8_cube_root_minus_one_inverse(void) {
    /* 1/(2^(1/3) - 1) = 4^(1/3) + 2^(1/3) + 1.
     * Algebraic identity: (a-1)(a^2+a+1) = a^3 - 1, with a = 2^(1/3),
     * a^3 = 2, so (a-1)(a^2+a+1) = 1. */
#ifdef PHASE_2
    assert_eval_eq(
        "Simplify[1/(2^(1/3) - 1) - (4^(1/3) + 2^(1/3) + 1)]",
        "0", 0);
#else
    SKIP_PHASE("Phase 2: algebraic-number rationalisation");
#endif
}

void test_radical_case10_mixed_extension_inverse(void) {
    /* 1/(sqrt(2) + 2^(1/3))
     *
     * The user's original transcript provided an r6 with sign errors
     * — it is NOT equal to 1/(sqrt(2) + 2^(1/3)) (verifiable
     * numerically). The correct rationalisation, derived via
     * extended Euclidean in Q[α]/(α^6 - 2) with α = 2^(1/6) the
     * primitive element of Q(sqrt(2), 2^(1/3)):
     *
     *   1/(sqrt(2) + 2^(1/3))
     *     = (2 sqrt(2) + sqrt(2)*2^(1/3) + sqrt(2)*4^(1/3)
     *        - 2*2^(1/3) - 4^(1/3) - 2) / 2.
     *
     * Mixed extension Q(sqrt(2), 2^(1/3)), dim 6 over Q. */
#ifdef PHASE_2
    assert_eval_eq(
        "Simplify[1/(Sqrt[2] + 2^(1/3)) - "
        "(2 Sqrt[2] + Sqrt[2]*2^(1/3) + Sqrt[2]*4^(1/3) "
        " - 2*2^(1/3) - 4^(1/3) - 2)/2]",
        "0", 0);
#else
    SKIP_PHASE("Phase 2: algebraic-number rationalisation");
#endif
}

/* ---- Phase 3 — Cube-root denesting ---- */

void test_radical_case4_cuberoot_denest(void) {
    /* (sqrt(5) + 2)^(1/3) = (sqrt(5) + 1)/2.
     * Borodin–Fagin–Hopcroft–Tompa: cubing (a + b sqrt(5))/2 with
     * a = b = 1 gives sqrt(5) + 2. */
#ifdef PHASE_3
    assert_eval_eq(
        "Simplify[(Sqrt[5] + 2)^(1/3) - (Sqrt[5] + 1)/2]",
        "0", 0);
#else
    SKIP_PHASE("Phase 3: cube-root denesting");
#endif
}

void test_radical_case5_sum_of_cuberoot_conjugates(void) {
    /* (2 + sqrt(5))^(1/3) + (2 - sqrt(5))^(1/3) = 1.
     * Cardano: let s be the sum. Then s^3 = 4 + 3((2+sqrt(5))(2-sqrt(5)))^(1/3) s
     * = 4 + 3*(-1)^(1/3) s = 4 - 3 s, so s^3 + 3 s - 4 = 0,
     * (s-1)(s^2+s+4) = 0, real root s = 1. */
#ifdef PHASE_3
    assert_eval_eq(
        "Simplify[(2 + Sqrt[5])^(1/3) + (2 - Sqrt[5])^(1/3) - 1]",
        "0", 0);
#else
    SKIP_PHASE("Phase 3: cube-root sum-of-conjugates");
#endif
}

int main(void) {
    symtab_init();
    core_init();

    /* Regression sentinel — runs in every phase. */
    TEST(test_radical_case9_zero_test_already_passes);

    /* Phase 1. */
    TEST(test_radical_case1_sqrt_3_plus_2sqrt2);
    TEST(test_radical_case2_sqrt_17_minus_12sqrt2);
    TEST(test_radical_case3_sqrt_2_plus_sqrt3);
    TEST(test_radical_case6_symbolic_two_term);
    TEST(test_radical_case7_symbolic_half_angle);
    TEST(test_radical_phase1_no_assumption_no_denest_symbolic);
    TEST(test_radical_phase1_negative_discriminant_no_fire);
    TEST(test_radical_phase1_denest_minus_branch);
    TEST(test_radical_phase1_idempotent_already_denested);

    /* Phase 2. */
    TEST(test_radical_case8_cube_root_minus_one_inverse);
    TEST(test_radical_case10_mixed_extension_inverse);

    /* Phase 3. */
    TEST(test_radical_case4_cuberoot_denest);
    TEST(test_radical_case5_sum_of_cuberoot_conjugates);

    if (g_failures > 0) {
        fprintf(stderr, "\n%d radical_simplify test(s) FAILED.\n", g_failures);
        return 1;
    }
    printf("All radical_simplify tests passed!\n");
    return 0;
}
