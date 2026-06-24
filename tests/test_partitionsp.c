/* Unit tests for PartitionsP (src/partitions.c).
 *
 * PartitionsP[n] is the partition-counting function p(n), computed with
 * Euler's pentagonal-number recurrence. These tests exercise the documented
 * examples, big-integer correctness, edge cases (zero/negative), Listable
 * threading, symbolic passthrough and the arg-count error.
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
    /* Table[PartitionsP[k], {k, 0, 12}] */
    check("Table[PartitionsP[k], {k, 0, 12}]",
          "{1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77}");
}

static void test_small_scalars(void) {
    check("PartitionsP[0]", "1");
    check("PartitionsP[1]", "1");
    check("PartitionsP[5]", "7");
    check("PartitionsP[10]", "42");
    check("PartitionsP[100]", "190569292");
}

/* The Table[PartitionsP[2^k], {k,0,12}] spec example, spot-checked at the
 * large end where the result is a big integer. */
static void test_bigint_values(void) {
    check("PartitionsP[256]", "365749566870782");
    check("PartitionsP[512]", "4453575699570940947378");
    check("PartitionsP[1024]",
          "61847822068260244309086870983975");
    check("PartitionsP[2048]",
          "18116048323611252751541173214616030020513022685");
    check("PartitionsP[4096]",
          "6927233917602120527467409170319882882996950147283323368445315320451");
}

/* The largest documented values must really be EXPR_BIGINT, not silently
 * truncated to a machine integer. */
static void test_bigint_type(void) {
    Expr* e = parse_expression("PartitionsP[4096]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_BIGINT);
    expr_free(e);
    expr_free(r);

    /* p(512) overflows int64, so it is a bigint too. (p(256) ~ 3.7e14 still
     * fits in int64 and normalizes to EXPR_INTEGER.) */
    Expr* e2 = parse_expression("PartitionsP[512]");
    Expr* r2 = evaluate(e2);
    ASSERT(r2->type == EXPR_BIGINT);
    expr_free(e2);
    expr_free(r2);
}

/* ----- edge cases ------------------------------------------------------ */

static void test_negative(void) {
    /* p(n) = 0 for n < 0 */
    check("PartitionsP[-1]", "0");
    check("PartitionsP[-5]", "0");
    check("PartitionsP[-100]", "0");
}

/* ----- Listable threading ---------------------------------------------- */

static void test_listable(void) {
    check("PartitionsP[{2, 4, 6}]", "{2, 5, 11}");
    /* nested lists thread to all depths */
    check("PartitionsP[{{1, 2}, {3, 4}}]", "{{1, 2}, {3, 5}}");
}

/* ----- symbolic / non-integer passthrough ------------------------------ */

static void test_symbolic(void) {
    check("PartitionsP[x]", "PartitionsP[x]");
    check("PartitionsP[1/2]", "PartitionsP[1/2]");
    check("PartitionsP[2.5]", "PartitionsP[2.5]");
}

/* ----- composite spec example: # of Abelian groups of order n ---------- */

static void test_abelian_groups(void) {
    check("Table[Times @@ PartitionsP[Last /@ FactorInteger[n]], {n, 12}]",
          "{1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2}");
}

/* ----- wrong arg count: emits PartitionsP::argx, stays unevaluated ------ */

static void test_argx(void) {
    check("PartitionsP[]", "PartitionsP[]");
    check("PartitionsP[3, 4]", "PartitionsP[3, 4]");
}

/* ----- HRR engine cross-validation against the recurrence --------------
 * The pentagonal recurrence (engine 1) is an independent, exact integer
 * computation. We use it as an oracle to validate the Hardy-Ramanujan-
 * Rademacher engine (engine 2): for every n the two must agree exactly.
 * This exhaustively exercises the HRR path far below its production
 * threshold and at a few large values. */
static void check_hrr_eq_recurrence(unsigned long n) {
    mpz_t a, b;
    mpz_init(a);
    mpz_init(b);
    ASSERT(partitionsp_recurrence(n, a) == 0);
    ASSERT(partitionsp_hrr(n, b) == 0);
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
    for (unsigned long n = 2; n <= 3000; n++) check_hrr_eq_recurrence(n);
    /* a scattering of larger values (recurrence still cheap here) */
    unsigned long big[] = {5000, 7777, 10000, 20000, 50000, 100000};
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
    TEST(test_abelian_groups);
    TEST(test_argx);
    TEST(test_hrr_matches_recurrence);

    printf("All PartitionsP tests passed!\n");
    return 0;
}
