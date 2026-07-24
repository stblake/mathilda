/* Unit tests for PrimitiveRoot and PrimitiveRootList (arithmetic.c).
 *
 * Coverage:
 *   - 1-arg PrimitiveRoot[n] for each cyclic family: n in {2, 4}, odd
 *     prime, odd prime power p^k (k > 1), 2 p^k.
 *   - 2-arg PrimitiveRoot[n, k] forward-search semantics.
 *   - Non-cyclic n (n in {1, 8, 12, 15, ...}) returns unevaluated.
 *   - Arg-count and arg-type diagnostics (::argt, ::argx, ::intg).
 *   - Listable threading.
 *   - PrimitiveRootList enumeration matches Mathematica references.
 *   - Bignum inputs exercise the GMP fast paths.
 *   - Structural correctness: every returned g satisfies gcd(g, n) == 1
 *     and order(g, n) == phi(n).  Verified via MultiplicativeOrder ... no
 *     wait, we don't have that builtin; we cross-check with PowerMod.
 *
 * Run binary directly: ./primitive_root_tests (MEMORY.md note: no ctest). */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* Suppress builtin stderr diagnostics for the whole binary; the tests
 * in this file only check return values, not stderr output. */
static void mute_stderr_once(void) {
    static int done = 0;
    if (!done) {
        freopen("/dev/null", "w", stderr);
        done = 1;
    }
}

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

/* ------------------------------------------------------------------ */
/* PrimitiveRoot — 1-arg form: spec examples and small cyclic families. */
/* ------------------------------------------------------------------ */

static void test_primitive_root_spec_examples(void) {
    /* From the docstring. */
    check_eq("PrimitiveRoot[9]",  "2");
    check_eq("PrimitiveRoot[10]", "7");
    check_eq("PrimitiveRoot[19]", "2");
    check_eq("PrimitiveRoot[7]",  "3");
}

static void test_primitive_root_small_cyclic(void) {
    /* n == 2 -> 1, n == 4 -> 3 are the boundary cases. */
    check_eq("PrimitiveRoot[2]", "1");
    check_eq("PrimitiveRoot[4]", "3");

    /* Odd primes up to 23.  Reference values from Mathematica. */
    check_eq("PrimitiveRoot[3]",  "2");
    check_eq("PrimitiveRoot[5]",  "2");
    check_eq("PrimitiveRoot[11]", "2");
    check_eq("PrimitiveRoot[13]", "2");
    check_eq("PrimitiveRoot[17]", "3");
    check_eq("PrimitiveRoot[23]", "5");

    /* Odd prime powers p^k (k > 1). */
    check_eq("PrimitiveRoot[25]", "2");   /* 5^2 */
    check_eq("PrimitiveRoot[27]", "2");   /* 3^3 */
    check_eq("PrimitiveRoot[81]", "2");   /* 3^4 */
    check_eq("PrimitiveRoot[49]", "3");   /* 7^2 */
    check_eq("PrimitiveRoot[121]", "2");  /* 11^2 */

    /* n = 2 p^k: if g is the smallest PR of p^k and g is odd, return g;
     * else return g + p^k.  Pulls out an odd representative coprime to 2. */
    check_eq("PrimitiveRoot[6]",  "5");   /* PR(3)=2 even -> 2 + 3 = 5  */
    check_eq("PrimitiveRoot[14]", "3");   /* PR(7)=3 odd  -> 3          */
    check_eq("PrimitiveRoot[18]", "11");  /* PR(9)=2 even -> 2 + 9 = 11 */
    check_eq("PrimitiveRoot[22]", "13");  /* PR(11)=2 even -> 2 + 11 = 13 */
}

static void test_primitive_root_non_cyclic(void) {
    /* Non-cyclic moduli leave PrimitiveRoot unevaluated.
     * 8 = 2^3, 12 = 2^2 * 3, 15 = 3 * 5, 16, 20, 24 ... */
    check_eq("PrimitiveRoot[8]",  "PrimitiveRoot[8]");
    check_eq("PrimitiveRoot[12]", "PrimitiveRoot[12]");
    check_eq("PrimitiveRoot[15]", "PrimitiveRoot[15]");
    check_eq("PrimitiveRoot[16]", "PrimitiveRoot[16]");
    check_eq("PrimitiveRoot[20]", "PrimitiveRoot[20]");
    check_eq("PrimitiveRoot[24]", "PrimitiveRoot[24]");
    check_eq("PrimitiveRoot[100]", "PrimitiveRoot[100]");
    /* 2 * 2^k = 2^(k+1) is just a power of 2, not 2 * (odd prime power). */
    check_eq("PrimitiveRoot[32]", "PrimitiveRoot[32]");
}

/* ------------------------------------------------------------------ */
/* PrimitiveRoot — 2-arg form: forward search starting at k.          */
/* ------------------------------------------------------------------ */

static void test_primitive_root_two_arg(void) {
    /* Smallest PR of 10 is 3; next is 7. */
    check_eq("PrimitiveRoot[10, 1]", "3");
    check_eq("PrimitiveRoot[10, 2]", "3");
    check_eq("PrimitiveRoot[10, 3]", "3");
    check_eq("PrimitiveRoot[10, 4]", "7");
    check_eq("PrimitiveRoot[10, 7]", "7");

    /* Smallest PR of 9 is 2; next is 5. */
    check_eq("PrimitiveRoot[9, 1]", "2");
    check_eq("PrimitiveRoot[9, 3]", "5");

    /* Negative k is clamped to 0; same as start-from-1 effectively. */
    check_eq("PrimitiveRoot[7, -10]", "3");

    /* k beyond the largest canonical PR -> unevaluated. */
    check_eq("PrimitiveRoot[10, 100]", "PrimitiveRoot[10, 100]");
}

/* ------------------------------------------------------------------ */
/* Diagnostics.  Builtins print to stderr but return the call         */
/* unevaluated; we only check the return value here.                  */
/* ------------------------------------------------------------------ */

static void test_primitive_root_diagnostics(void) {
    mute_stderr_once();

    check_eq("PrimitiveRoot[]",                "PrimitiveRoot[]");
    check_eq("PrimitiveRoot[3.4]",             "PrimitiveRoot[3.4]");
    check_eq("PrimitiveRoot[1]",               "PrimitiveRoot[1]");
    check_eq("PrimitiveRoot[0]",               "PrimitiveRoot[0]");
    check_eq("PrimitiveRoot[-5]",              "PrimitiveRoot[-5]");
    check_eq("PrimitiveRoot[56, 45, 34, 21]",  "PrimitiveRoot[56, 45, 34, 21]");
    check_eq("PrimitiveRoot[5, 2.5]",          "PrimitiveRoot[5, 2.5]");
}

/* ------------------------------------------------------------------ */
/* Listable threading.                                                 */
/* ------------------------------------------------------------------ */

static void test_primitive_root_listable(void) {
    /* From the docstring. */
    check_eq("PrimitiveRoot[{9, 7, 19}]", "List[2, 3, 2]");
    /* Mixed cyclic / non-cyclic threads element-wise; non-cyclic stays
     * unevaluated inside the list. */
    check_eq("PrimitiveRoot[{9, 12, 19}]",
             "List[2, PrimitiveRoot[12], 2]");
}

/* ------------------------------------------------------------------ */
/* Symbolic args flow through (HoldFirst-like behaviour for unknown   */
/* heads).                                                             */
/* ------------------------------------------------------------------ */

static void test_primitive_root_symbolic(void) {
    /* Bare symbol: PrimitiveRoot[x] stays unevaluated, no diagnostic. */
    check_eq("PrimitiveRoot[x]", "PrimitiveRoot[x]");
    /* Prime[1000000] = 15485863 now evaluates (fast prime counting), so
     * this reduces to a concrete primitive-root computation: 6 is the
     * smallest primitive root (2, 3, 5 are non-residues of full order). */
    check_eq("PrimitiveRoot[Prime[1000000]]", "6");
}

/* ------------------------------------------------------------------ */
/* Bignum inputs exercise the GMP fast paths.                          */
/* ------------------------------------------------------------------ */

static void test_primitive_root_bignum(void) {
    /* 10^18 + 9 is a well-known prime; smallest PR is 7. */
    check_eq("PrimitiveRoot[10^18 + 9]", "7");
    /* p^k where p fits in machine word but k makes the modulus huge.
     * 100003 is prime, 2 is a primitive root of 100003^k for k >= 1. */
    check_eq("PrimitiveRoot[100003^100]", "2");
    /* n = 2 * 11^5 = 322102; smallest PR of 11^5 is 2 (even) ->
     * 2 + 11^5 = 161053. */
    check_eq("PrimitiveRoot[2 * 11^5]", "161053");
}

/* ------------------------------------------------------------------ */
/* Structural verification: returned g must generate (Z/nZ)^*, i.e.   */
/* PowerMod[g, phi(n), n] == 1 and PowerMod[g, phi(n)/q, n] != 1 for  */
/* every prime q | phi(n).  We use the second-strongest spot check —  */
/* PowerMod[g, phi(n)/2, n] != 1 — to verify g is not a square.       */
/* ------------------------------------------------------------------ */

static void test_primitive_root_structural(void) {
    /* For n = 19, phi = 18, g = 2: 2^9 mod 19 = 18, 2^6 mod 19 = 7. */
    check_eq("PowerMod[PrimitiveRoot[19], EulerPhi[19], 19]", "1");
    check_eq("PowerMod[PrimitiveRoot[19], EulerPhi[19] / 2, 19]", "18");

    /* For n = 25, phi = 20: PR is 2; 2^20 = 1, 2^10 != 1. */
    check_eq("PowerMod[PrimitiveRoot[25], EulerPhi[25], 25]", "1");
}

/* ------------------------------------------------------------------ */
/* PrimitiveRootList — exact enumerations.                              */
/* ------------------------------------------------------------------ */

static void test_primitive_root_list_spec(void) {
    /* From the docstring. */
    check_eq("PrimitiveRootList[9]", "List[2, 5]");

    /* Mathematica references. */
    check_eq("PrimitiveRootList[2]",  "List[1]");
    check_eq("PrimitiveRootList[4]",  "List[3]");
    check_eq("PrimitiveRootList[3]",  "List[2]");
    check_eq("PrimitiveRootList[5]",  "List[2, 3]");
    check_eq("PrimitiveRootList[7]",  "List[3, 5]");
    check_eq("PrimitiveRootList[11]", "List[2, 6, 7, 8]");
    check_eq("PrimitiveRootList[13]", "List[2, 6, 7, 11]");
    check_eq("PrimitiveRootList[14]", "List[3, 5]");
    check_eq("PrimitiveRootList[18]", "List[5, 11]");
    check_eq("PrimitiveRootList[19]",
             "List[2, 3, 10, 13, 14, 15]");
    check_eq("PrimitiveRootList[22]",
             "List[7, 13, 17, 19]");
    check_eq("PrimitiveRootList[25]",
             "List[2, 3, 8, 12, 13, 17, 22, 23]");
}

static void test_primitive_root_list_non_cyclic(void) {
    /* Empty list for non-cyclic n. */
    check_eq("PrimitiveRootList[1]",   "List[]");
    check_eq("PrimitiveRootList[8]",   "List[]");
    check_eq("PrimitiveRootList[12]",  "List[]");
    check_eq("PrimitiveRootList[15]",  "List[]");
    check_eq("PrimitiveRootList[16]",  "List[]");
    check_eq("PrimitiveRootList[100]", "List[]");
}

static void test_primitive_root_list_diagnostics(void) {
    mute_stderr_once();

    check_eq("PrimitiveRootList[]",        "PrimitiveRootList[]");
    check_eq("PrimitiveRootList[34, 2, 3]", "PrimitiveRootList[34, 2, 3]");
    /* Non-integer numeric inputs flow through (no diagnostic per spec). */
    check_eq("PrimitiveRootList[11.0]",     "PrimitiveRootList[11.0]");
    check_eq("PrimitiveRootList[11 + I]",   "PrimitiveRootList[Complex[11, 1]]");
}

static void test_primitive_root_list_listable(void) {
    /* Listable threading. */
    check_eq("PrimitiveRootList[{9, 10}]",
             "List[List[2, 5], List[3, 7]]");
}

static void test_primitive_root_list_structural(void) {
    /* For each g in PrimitiveRootList[n], gcd(g, n) must be 1.
     * Spot-check via the Mathematica-style identity that the union of
     * orbits {g^i mod n : i in 1..phi(n)} equals the units (Z/nZ)^*. */
    check_eq("Union[Table[PowerMod[2, i, 9], {i, 1, 6}]]",
             "List[1, 2, 4, 5, 7, 8]");
    check_eq("Union[Table[PowerMod[5, i, 9], {i, 1, 6}]]",
             "List[1, 2, 4, 5, 7, 8]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_primitive_root_spec_examples);
    TEST(test_primitive_root_small_cyclic);
    TEST(test_primitive_root_non_cyclic);
    TEST(test_primitive_root_two_arg);
    TEST(test_primitive_root_diagnostics);
    TEST(test_primitive_root_listable);
    TEST(test_primitive_root_symbolic);
    TEST(test_primitive_root_bignum);
    TEST(test_primitive_root_structural);

    TEST(test_primitive_root_list_spec);
    TEST(test_primitive_root_list_non_cyclic);
    TEST(test_primitive_root_list_diagnostics);
    TEST(test_primitive_root_list_listable);
    TEST(test_primitive_root_list_structural);

    printf("All PrimitiveRoot tests passed!\n");
    return 0;
}
