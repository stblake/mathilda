/*
 * tests/test_nf_rowreduce.c
 *
 * Unit tests for the native number-field RowReduce fast path
 * (flint_qqbar_nf_mat_rref, wired into builtin_rowreduce): a matrix whose
 * entries are rationals and AlgebraicNumber[theta, ...] over ONE common theta
 * is reduced with FLINT antic nf_elem arithmetic instead of per-element
 * evaluator dispatch.
 *
 * The reduced row echelon form is canonical, so the native result must be
 * IDENTICAL to the classical AlgebraicNumber RowReduce (any Method) -- these
 * expected values were captured from the classical path (MATHILDA_NO_NF_RREF=1)
 * and must match with the fast path on (the default here).
 *
 * The randomized byte-identical differential (native vs MATHILDA_NO_NF_RREF over
 * 120 matrices across Q(sqrt2)/Q(i)/Q(2^(1/3))/Q(cubic)) is a separate shell
 * gate; this suite pins fixed cases plus RREF self-consistency (idempotence).
 *
 * Memory: assert_eval_eq frees its trees; the binary is valgrind-clean.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Rationals are handled by the fmpq_mat path; the field fast path must not
 * perturb them. */
static void test_rational_unchanged(void) {
    assert_eval_eq("RowReduce[{{1, 2}, {3, 4}}]", "{{1, 0}, {0, 1}}", 0);
    assert_eval_eq("RowReduce[{{2, 3}, {4, 6}}]", "{{1, 3/2}, {0, 0}}", 0);
    /* an explicit Method must still reach the same canonical RREF */
    assert_eval_eq("RowReduce[{{2, 3}, {4, 6}}, Method -> \"OneStepRowReduction\"]",
                   "{{1, 3/2}, {0, 0}}", 0);
}

/* Q(sqrt2): rank-deficient 2x2, pivot scaled by a field inverse. */
static void test_q_sqrt2(void) {
    assert_eval_eq(
        "RowReduce[{{AlgebraicNumber[Sqrt[2], {0, 1}], 1}, "
        "{2, AlgebraicNumber[Sqrt[2], {0, 1}]}}]",
        "{{1, AlgebraicNumber[Sqrt[2], {0, 1/2}]}, {0, 0}}", 0);
    /* proportional rows -> rank 1 */
    assert_eval_eq(
        "RowReduce[{{AlgebraicNumber[Sqrt[2], {0, 1}], 2}, "
        "{AlgebraicNumber[Sqrt[2], {0, 2}], 4}}]",
        "{{1, AlgebraicNumber[Sqrt[2], {0, 1}]}, {0, 0}}", 0);
}

/* Q(i): a full-rank Gaussian matrix reduces to the identity. */
static void test_q_i(void) {
    assert_eval_eq("RowReduce[{{I, 1, 2}, {1, I, 3}, {2, 3, I}}]",
                   "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}", 0);
    /* Q(sqrt(-5)) full rank -> identity */
    assert_eval_eq(
        "RowReduce[{{AlgebraicNumber[Sqrt[-5], {0, 1}], 4}, "
        "{2, AlgebraicNumber[Sqrt[-5], {0, 2}]}}]",
        "{{1, 0}, {0, 1}}", 0);
}

/* A cubic field Q(2^(1/3)): full rank -> identity (exercises degree > 2). */
static void test_cubic_field(void) {
    assert_eval_eq(
        "RowReduce[{{AlgebraicNumber[2^(1/3), {0, 1, 0}], 1}, "
        "{1, AlgebraicNumber[2^(1/3), {0, 0, 1}]}}]",
        "{{1, 0}, {0, 1}}", 0);
}

/* RREF is idempotent: reducing the reduced form is a no-op. A strong
 * self-consistency check independent of the pinned expected strings. */
static void test_idempotent(void) {
    assert_eval_eq(
        "With[{m = {{AlgebraicNumber[Sqrt[2], {1, 1}], 2, 3}, "
        "{1, AlgebraicNumber[Sqrt[2], {0, 1}], 5}, "
        "{2, 3, AlgebraicNumber[Sqrt[2], {1, 0}]}}}, "
        "RowReduce[m] === RowReduce[RowReduce[m]]]",
        "True", 0);
    /* mixed rational + field entries, non-square */
    assert_eval_eq(
        "With[{m = {{AlgebraicNumber[Sqrt[3], {0, 1}], 2, 1, 0}, "
        "{4, AlgebraicNumber[Sqrt[3], {1, 1}], 0, 1}}}, "
        "RowReduce[m] === RowReduce[RowReduce[m]]]",
        "True", 0);
}

/* RowReduce[m, ZeroTest -> f]: f a predicate (RootReduce-based body or a head
 * such as PossibleZeroQ), consulted for every zero decision.  Shares the
 * machinery with NullSpace (matsol_rref_with_zerotest). */
static void test_rowreduce_zerotest(void) {
    assert_eval_eq("RowReduce[{{2, 3}, {4, 6}}, ZeroTest -> PossibleZeroQ]",
                   "{{1, 3/2}, {0, 0}}", 0);
    assert_eval_eq(
        "RowReduce[{{2, 3}, {4, 6}}, ZeroTest -> (RootReduce[Together[#]] === 0 &)]",
        "{{1, 3/2}, {0, 0}}", 0);
    /* algebraic rank deficiency (row2 = Sqrt[2]*row1) */
    assert_eval_eq(
        "RowReduce[{{1, Sqrt[2]}, {Sqrt[2], 2}}, ZeroTest -> PossibleZeroQ]",
        "{{1, Sqrt[2]}, {0, 0}}", 0);
    /* Method and ZeroTest combined */
    assert_eval_eq(
        "RowReduce[{{2, 3}, {4, 6}}, Method -> \"OneStepRowReduction\", "
        "ZeroTest -> PossibleZeroQ]",
        "{{1, 3/2}, {0, 0}}", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_rational_unchanged);
    TEST(test_rowreduce_zerotest);
    TEST(test_q_sqrt2);
    TEST(test_q_i);
    TEST(test_cubic_field);
    TEST(test_idempotent);

    printf("All native number-field RowReduce tests passed!\n");
    return 0;
}
