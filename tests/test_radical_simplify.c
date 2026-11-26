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

/* ---- Regression: single-radical fraction integrand recovery ---- */

void test_radical_integrand_recovery_sqrt_quadratic(void) {
    /* D[Integrate[Sqrt[x^2-x-1]/x, x], x] // Simplify must recover the
     * integrand exactly.  Previously Simplify returned a value 4x too
     * large because Factor's bivariate-Hensel path dropped the integer
     * content of the generator-polynomial it factored (see the facpoly
     * test_factor_bivariate_content_regression).  Pin both the value and
     * the canonical single-radical form. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[x^2 - x - 1]/x, x], x]]",
        "Sqrt[-1 - x + x^2]/x", 0);
    /* The factor-of-4 corruption in isolation: Simplify of the combined
     * fraction must preserve value. */
    assert_eval_eq(
        "Simplify[Together[D[Integrate[Sqrt[x^2 - x - 1]/x, x], x]]]",
        "Sqrt[-1 - x + x^2]/x", 0);
}

/* Algebraic-function-field canonical form (simp_algebraic_field): Simplify of a
 * re-differentiated algebraic antiderivative returns the integrand. */
static void test_radical_function_field_roundtrips(void) {
    check_eval_eq("Simplify[D[Integrate[1/(x^4 + 1)^(1/4), x, Method -> \"Trager\"], x]]",
                  "1/(1 + x^4)^(1/4)");
    check_eval_eq("Simplify[D[Integrate[1/(x^4 + 1)^(1/4), x, Method -> \"ParallelMixedTower\"], x]]",
                  "1/(1 + x^4)^(1/4)");
    check_eval_eq("Simplify[D[Integrate[x^(4/3)/(a Sqrt[x] + b)^(2/3), x], x]]",
                  "x^(4/3)/(b + a Sqrt[x])^(2/3)");
    check_eval_eq("Simplify[D[Integrate[1/(x^3 (a + b x)^(1/3)), x], x]]",
                  "1/(x^3 (a + b x)^(1/3))");
    /* Non-elementary (Chebychev): since v0.6763 the antiderivative is a 2F1
     * (x^2/2 2F1(1/3, 2/3; 5/3; -x^3)), which Simplify does not yet reduce
     * back to the radical.  Its derivative must still be the integrand, and
     * Simplify must keep it finite: a degenerate contiguous shift (c = b) made
     * it Indeterminate. */
    check_eval_eq("PossibleZeroQ[Simplify[D[Integrate[x/(1 + x^3)^(1/3), x], x]] - x/(1 + x^3)^(1/3)]",
                  "True");
    check_eval_eq("Simplify[D[Integrate[Sqrt[x]/(1 + x^(1/3)), x], x]]",
                  "Sqrt[x]/(1 + x^(1/3))");
    check_eval_eq("Simplify[D[Integrate[Sqrt[1 + Sqrt[1 - x]], x], x]]",
                  "Sqrt[1 + Sqrt[1 - x]]");
    check_eval_eq("Simplify[(x - a)/(x^(1/3) - a^(1/3))]",
                  "a^(2/3) + a^(1/3) x^(1/3) + x^(2/3)");
    /* soundness: branch-dependent non-identities are left alone */
    check_eval_eq("Simplify[Sqrt[x y] - Sqrt[x] Sqrt[y]]", "-Sqrt[x] Sqrt[y] + Sqrt[x y]");
    check_eval_eq("Simplify[Sqrt[x^2] - x]", "-x + Sqrt[x^2]");
    /* assumptions still apply to the canonical form */
    check_eval_eq("Simplify[Sqrt[x^2] + (x^2 - 1)/(x - 1), x > 0]", "1 + 2 x");
}

/* The field over transcendental atoms: E^(k x), b^x, Log[x], Sin[x], ...
 * enter Simplify's algebraic-function-field canonical form as parameters
 * (exponentials as integer powers of one generator, trig over E^(I x), the
 * imaginary unit as the root of unity omega). Out-of-sample sweep:
 * tools/simplify_atoms/gen.m. */
static void test_radical_function_field_transcendental_atoms(void) {
    /* the reported case: Sqrt[a] and Sqrt[a + E^(3x)] over Q(a, E^x) */
    check_eval_eq("Simplify[-1/Sqrt[a] - E^(3 x)/(a Sqrt[a + E^(3 x)] - Sqrt[a] (a + E^(3 x)))]",
                  "1/Sqrt[a + E^(3 x)]");
    check_eval_eq("Simplify[D[Integrate[1/Sqrt[E^(3 x) + a], x, Method -> \"Bronstein\"], x]]",
                  "1/Sqrt[a + E^(3 x)]");
    /* commensurate exponentials share one generator */
    check_eval_eq("Simplify[D[Integrate[E^x/Sqrt[a + E^(2 x)], x, Method -> \"Bronstein\"], x]]",
                  "E^x/Sqrt[a + E^(2 x)]");
    check_eval_eq("Simplify[D[Integrate[E^(2 x)/Sqrt[b + E^x], x, Method -> \"Bronstein\"], x]]",
                  "E^(2 x)/Sqrt[b + E^x]");
    check_eval_eq("Simplify[D[Integrate[2^x/Sqrt[2 + 2^x], x], x]]", "2^x/Sqrt[2 + 2^x]");
    /* kernels: Log[x] is an opaque parameter */
    check_eval_eq("Simplify[-1/Sqrt[a] - Log[x]/(a Sqrt[a + Log[x]] - Sqrt[a] (a + Log[x]))]",
                  "1/Sqrt[a + Log[x]]");
    check_eval_eq("Simplify[D[Integrate[1/(x Log[x] Sqrt[a + Log[x]]), x, Method -> \"Bronstein\"], x]]",
                  "1/(x Log[x] Sqrt[a + Log[x]])");
    /* trig over E^(I x): Bronstein's radicand is 3 + Sin[x] in disguise */
    check_eval_eq("Simplify[D[Integrate[Cos[x]/Sqrt[3 + Sin[x]], x, Method -> \"Bronstein\"], x]]",
                  "Cos[x]/Sqrt[3 + Sin[x]]");
    /* the imaginary unit as omega, ordered below the radical it sits in
     * (with omega leading this was a wrong answer) */
    check_eval_eq("Simplify[-1/Sqrt[a] - I E^x/(a Sqrt[a + I E^x] - Sqrt[a] (a + I E^x))]",
                  "1/Sqrt[a + I E^x]");
    /* no atom in any radicand: the layer stays off, the general search keeps it */
    check_eval_eq("Simplify[-1/Sqrt[a] - I y/(a Sqrt[a + I y] - Sqrt[a] (a + I y))]",
                  "-1/Sqrt[a] - (I y)/(a Sqrt[a + I y] - Sqrt[a] (a + I y))");
    check_eval_eq("DSolve[{D[u[y,t],t] == nu D[u[y,t],{y,2}], u[y,0] == 0, u[0,t] == U}, u[y,t], {y,t}]",
                  "{{u[y, t] -> U - U Erf[y/(2 Sqrt[nu] Sqrt[t])]}}");
    check_eval_eq("Simplify[(Sqrt[1 + (-1)^(1/3) y] - 1) (Sqrt[1 + (-1)^(1/3) y] + 1)]",
                  "(-1)^(1/3) y");
    /* an atom beside an ordinary variable */
    check_eval_eq("Simplify[(Sqrt[x + E^x] - x) (Sqrt[x + E^x] + x) - E^x]", "(1 - x) x");
    /* no-ops and declines */
    check_eval_eq("Simplify[Sqrt[Sin[x]]]", "Sqrt[Sin[x]]");
    check_eval_eq("Simplify[1/(Sqrt[x] + I)]", "1/(I + Sqrt[x])");
    /* RootReduce keeps the plain field (atoms off) */
    check_eval_eq("RootReduce[(Sqrt[a + E^x] - Sqrt[a]) (Sqrt[a + E^x] + Sqrt[a])]",
                  "(-Sqrt[a] + Sqrt[a + E^x]) (Sqrt[a] + Sqrt[a + E^x])");
}

int main(void) {
    symtab_init();
    core_init();

    /* Regression sentinel — runs in every phase. */
    TEST(test_radical_case9_zero_test_already_passes);
    TEST(test_radical_integrand_recovery_sqrt_quadratic);

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

    /* Algebraic function field. */
    TEST(test_radical_function_field_roundtrips);
    TEST(test_radical_function_field_transcendental_atoms);

    if (g_failures > 0) {
        fprintf(stderr, "\n%d radical_simplify test(s) FAILED.\n", g_failures);
        return 1;
    }
    printf("All radical_simplify tests passed!\n");
    return 0;
}
