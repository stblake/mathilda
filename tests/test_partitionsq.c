/* Unit tests for PartitionsQ (src/partitions.c).
 *
 * PartitionsQ[n] is q(n), the number of partitions of n into distinct parts
 * (= into odd parts; OEIS A000009). Two engines: an exact pentagonal-
 * convolution recurrence and the Hardy-Ramanujan-Rademacher / Hagis series.
 * These tests exercise the documented examples, big-integer correctness,
 * edge cases (zero/negative), Listable threading, symbolic passthrough, the
 * arg-count error, and an exhaustive HRR-equals-recurrence cross-check.
 *
 * We use the hard-failing ASSERT* macros (not libc assert) so failures abort
 * with a non-zero exit code even under CMake's -DNDEBUG Release build. */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "partitions.h"
#include "test_utils.h"
#include <stdio.h>
#include <gmp.h>

/* Evaluate `input`, assert the printed result equals `expected`. */
static void check(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(result);
    ASSERT_STR_EQ(str, expected);
    free(str);
    expr_free(result);
}

/* ----- documented examples (exact WL output) --------------------------- */

static void test_base_table(void) {
    /* Table[PartitionsQ[k], {k, 0, 20}] */
    check("Table[PartitionsQ[k], {k, 0, 20}]",
          "{1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 27, 32, 38, "
          "46, 54, 64}");
}

static void test_small_scalars(void) {
    check("PartitionsQ[0]", "1");
    check("PartitionsQ[1]", "1");
    check("PartitionsQ[2]", "1");
    check("PartitionsQ[5]", "3");
    check("PartitionsQ[10]", "10");
    check("PartitionsQ[100]", "444793");
}

/* Large values: q(1000) and q(100000) exercise the HRR path through the
 * builtin (threshold is 1000), and must be exact big integers. */
static void test_bigint_values(void) {
    check("PartitionsQ[1000]", "8635565795744155161506");
    check("PartitionsQ[100000]",
          "42494159403332317292526619504218136903700576932083624292980870857"
          "93661601651601912151502208964867232719338338068057175972722741603"
          "68211837446740514571940417111414290856263711241960579022839958369"
          "76239181670821800000403741232325992196887134172550");
}

/* q(1000) overflows int64, so it must really be EXPR_BIGINT. */
static void test_bigint_type(void) {
    Expr* e = parse_expression("PartitionsQ[1000]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_BIGINT);
    expr_free(e);
    expr_free(r);
}

/* ----- edge cases ------------------------------------------------------ */

static void test_negative(void) {
    /* q(n) = 0 for n < 0 */
    check("PartitionsQ[-1]", "0");
    check("PartitionsQ[-5]", "0");
    check("PartitionsQ[-100]", "0");
}

/* ----- Listable threading ---------------------------------------------- */

static void test_listable(void) {
    check("PartitionsQ[{2, 4, 6}]", "{1, 2, 4}");
    /* nested lists thread to all depths */
    check("PartitionsQ[{{2, 4}, {6, 10}}]", "{{1, 2}, {4, 10}}");
}

/* ----- symbolic / non-integer passthrough ------------------------------ */

static void test_symbolic(void) {
    check("PartitionsQ[x]", "PartitionsQ[x]");
    check("PartitionsQ[1/2]", "PartitionsQ[1/2]");
    check("PartitionsQ[12.1]", "PartitionsQ[12.1]");
}

/* ----- wrong arg count: emits PartitionsQ::argx, stays unevaluated ------ */

static void test_argx(void) {
    check("PartitionsQ[]", "PartitionsQ[]");
    check("PartitionsQ[3, 4]", "PartitionsQ[3, 4]");
}

/* ----- HRR engine cross-validation against the recurrence --------------
 * The pentagonal-convolution recurrence (engine 1) is an independent, exact
 * integer computation. We use it as an oracle to validate the Hagis /
 * Hardy-Ramanujan-Rademacher engine (engine 2): for every n the two must
 * agree exactly. This is the correctness guarantee for the HRR path, which
 * has no clean closed-form remainder bound. */
static void check_hrr_eq_recurrence(unsigned long n) {
    mpz_t a, b;
    mpz_init(a);
    mpz_init(b);
    ASSERT(partitionsq_recurrence(n, a) == 0);
    ASSERT(partitionsq_hrr(n, b) == 0);
    if (mpz_cmp(a, b) != 0) {
        gmp_fprintf(stderr,
            "FAIL: HRR != recurrence at n=%lu\n  recurrence: %Zd\n  hrr:        %Zd\n",
            n, a, b);
        ASSERT(0);
    }
    mpz_clear(a);
    mpz_clear(b);
}

static void test_hrr_matches_recurrence(void) {
    /* exhaustive small range (HRR is valid for n >= 2) */
    for (unsigned long n = 2; n <= 1500; n++) check_hrr_eq_recurrence(n);
    /* strided sample further out */
    for (unsigned long n = 1500; n <= 6000; n += 7) check_hrr_eq_recurrence(n);
    /* a scattering of larger values (recurrence still cheap here) */
    unsigned long big[] = {7777, 10000, 20000, 50000, 100000};
    for (size_t i = 0; i < sizeof(big) / sizeof(big[0]); i++)
        check_hrr_eq_recurrence(big[i]);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_base_table);
    TEST(test_small_scalars);
    TEST(test_bigint_values);
    TEST(test_bigint_type);
    TEST(test_negative);
    TEST(test_listable);
    TEST(test_symbolic);
    TEST(test_argx);
    TEST(test_hrr_matches_recurrence);

    printf("All PartitionsQ tests passed!\n");
    return 0;
}
