#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Tests for phase-2 multi-extension nested-radical denesting:
 *
 *   Sqrt[A + b*Sqrt[C]]  with A in Q(Sqrt[k]) and C also potentially
 *   carrying nested radicals.
 *
 * Phase 1 (the original denester in simp.c) handles single-extension
 * inputs where A is a rational and there is exactly one sqrt-bearing
 * term. Phase 2 extends both directions:
 *
 *   (a) the splitter (split_plus_into_a_plus_b_sqrt_c) accepts
 *       radicands with multiple sqrt-bearing terms and iterates over
 *       candidates for the "outer" sqrt, ranked deepest-first;
 *   (b) sqrt_if_clean_square recurses into the denester (with a
 *       depth budget) when its discriminant D is itself a Q(Sqrt[k])
 *       element.
 *
 * The motivating expression is from the user's bug report:
 *
 *   Sqrt[16 - 2 Sqrt[29] + 2 Sqrt[55 - 10 Sqrt[29]]]
 *     -> Sqrt[5] + Sqrt[11 - 2 Sqrt[29]]
 *
 * which Mathematica simplifies in ~100ms but earlier picocas left
 * untouched. This file exercises that case plus controls that confirm
 * (a) phase-1 inputs still denest with byte-identical output,
 * (b) non-denestable multi-sqrt inputs remain unchanged, and
 * (c) the recursion budget terminates.
 */

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

/* ---------------------------------------------------------------- */
/* The motivating case from the user's bug report.                  */
/* ---------------------------------------------------------------- */

void test_user_report_multi_extension_difference(void) {
    /* The exact expression the user filed:
     *   e1 = Sqrt[16 - 2 Sqrt[29] + 2 Sqrt[55 - 10 Sqrt[29]]]
     *   e2 = Sqrt[5] + Sqrt[11 - 2 Sqrt[29]]
     *   Simplify[e1 - e2] must reduce to 0. */
    check_eval_eq(
        "Simplify[Sqrt[16 - 2 Sqrt[29] + 2 Sqrt[55 - 10 Sqrt[29]]]"
        " - Sqrt[5] - Sqrt[11 - 2 Sqrt[29]]]",
        "0");
}

void test_motivating_radicand_denests_to_pair(void) {
    /* Direct denesting of the motivating radicand. */
    check_eval_eq(
        "Simplify[Sqrt[16 - 2 Sqrt[29] + 2 Sqrt[55 - 10 Sqrt[29]]]]",
        "Sqrt[5] + Sqrt[11 - 2 Sqrt[29]]");
}

/* ---------------------------------------------------------------- */
/* Inner step: the discriminant 152 - 24 Sqrt[29] -> 2 Sqrt[29] - 6. */
/* This must keep working under the single-extension path; phase-2  */
/* would not be triggered here.                                     */
/* ---------------------------------------------------------------- */

void test_phase1_inner_step_still_works(void) {
    check_eval_eq(
        "Simplify[Sqrt[152 - 24 Sqrt[29]]]",
        "-6 + 2 Sqrt[29]");
}

/* ---------------------------------------------------------------- */
/* Phase-1 sanity: single-extension denesting unaffected by phase-2. */
/* ---------------------------------------------------------------- */

void test_phase1_5_plus_2sqrt6(void) {
    check_eval_eq("Simplify[Sqrt[5 + 2 Sqrt[6]]]", "Sqrt[2] + Sqrt[3]");
}

void test_phase1_7_plus_4sqrt3(void) {
    check_eval_eq("Simplify[Sqrt[7 + 4 Sqrt[3]]]", "2 + Sqrt[3]");
}

void test_phase1_3_plus_2sqrt2(void) {
    check_eval_eq("Simplify[Sqrt[3 + 2 Sqrt[2]]]", "1 + Sqrt[2]");
}

/* ---------------------------------------------------------------- */
/* Non-denestable phase-2 inputs must remain unchanged (no wrong   */
/* branch, no crash).                                              */
/* ---------------------------------------------------------------- */

void test_unfireable_multi_sqrt_stays_put(void) {
    /* 1 + Sqrt[2] + Sqrt[3] has no half-sum denesting -- two
     * sqrt-bearing terms with independent radicals, no closed-form
     * discriminant in Q(Sqrt[2]) or Q(Sqrt[3]) reduces to a clean
     * square. The denester must leave it alone. */
    check_eval_eq(
        "Simplify[Sqrt[1 + Sqrt[2] + Sqrt[3]]]",
        "Sqrt[1 + Sqrt[2] + Sqrt[3]]");
}

void test_negative_discriminant_no_fire(void) {
    /* Single-sqrt case where A^2 - b^2 C < 0; the denester must
     * refuse rather than emit a complex result. */
    check_eval_eq(
        "Simplify[Sqrt[1 + Sqrt[5]]]",
        "Sqrt[1 + Sqrt[5]]");
}

/* ---------------------------------------------------------------- */
/* Idempotency: a phase-2 denested result should be unchanged under */
/* a second Simplify pass.                                          */
/* ---------------------------------------------------------------- */

void test_idempotent_phase2_result(void) {
    /* The denested form Sqrt[5] + Sqrt[11 - 2 Sqrt[29]] is already
     * minimal under picocas's complexity measure; Simplify on it
     * must not re-introduce the nested radicand. */
    check_eval_eq(
        "Simplify[Sqrt[5] + Sqrt[11 - 2 Sqrt[29]]]",
        "Sqrt[5] + Sqrt[11 - 2 Sqrt[29]]");
}

int main(void) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    symtab_init();
    core_init();

    TEST(test_user_report_multi_extension_difference);
    TEST(test_motivating_radicand_denests_to_pair);
    TEST(test_phase1_inner_step_still_works);
    TEST(test_phase1_5_plus_2sqrt6);
    TEST(test_phase1_7_plus_4sqrt3);
    TEST(test_phase1_3_plus_2sqrt2);
    TEST(test_unfireable_multi_sqrt_stays_put);
    TEST(test_negative_discriminant_no_fire);
    TEST(test_idempotent_phase2_result);

    if (g_failures > 0) {
        fprintf(stderr, "\n*** %d test(s) failed ***\n", g_failures);
        return 1;
    }
    printf("\nAll phase-2 denesting tests passed.\n");
    return 0;
}
