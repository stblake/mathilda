/* Unit tests for MultiplicativeOrder (arithmetic.c).
 *
 * Coverage:
 *   - 2-arg form: smallest positive m with k^m == 1 (mod n).
 *   - 3-arg form: discrete-log search over a residue list.
 *   - Edge cases: n == 1, |n|, negative or out-of-range k, non-coprime
 *     gcd, empty list, list with non-integer / out-of-group elements.
 *   - Arg-count diagnostic (::argt).
 *   - Non-integer numeric inputs flow through unevaluated silently.
 *   - GMP / bignum k or n exercised via 10^N exponents.
 *   - Structural identities cross-checked against PowerMod and EulerPhi.
 *
 * Run binary directly: ./multiplicative_order_tests
 * (per MEMORY.md note: ctest is not configured in tests/CMakeLists.txt). */

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
/* 2-arg form — spec examples and small moduli.                        */
/* ------------------------------------------------------------------ */

static void test_mo_spec_examples(void) {
    /* From the docstring. */
    check_eq("MultiplicativeOrder[5, 8]",  "2");
    check_eq("MultiplicativeOrder[5, 7]",  "6");
    check_eq("MultiplicativeOrder[-5, 7]", "3");      /* -5 == 2 mod 7 */
    check_eq("MultiplicativeOrder[-5, 3]", "1");      /* -5 == 1 mod 3 */
    check_eq("MultiplicativeOrder[5, 3]",  "2");
    check_eq("MultiplicativeOrder[10, 21]", "6");
    check_eq("MultiplicativeOrder[21, 10]", "1");     /* 21 == 1 mod 10 */
}

static void test_mo_trivial_modulus(void) {
    /* n == 1: trivial group, order is 1 for every k. */
    check_eq("MultiplicativeOrder[5, 1]", "1");
    check_eq("MultiplicativeOrder[0, 1]", "1");
    check_eq("MultiplicativeOrder[1, 1]", "1");
    /* n == 2: phi == 1, only unit is 1. */
    check_eq("MultiplicativeOrder[1, 2]", "1");
    check_eq("MultiplicativeOrder[3, 2]", "1");
    /* k == 1 is always order 1 (1^1 == 1 (mod n)). */
    check_eq("MultiplicativeOrder[1, 7]", "1");
    check_eq("MultiplicativeOrder[1, 100]", "1");
}

static void test_mo_negative_n(void) {
    /* Negative n is treated as |n|; same group structure. */
    check_eq("MultiplicativeOrder[5, -8]", "2");
    check_eq("MultiplicativeOrder[5, -7]", "6");
}

static void test_mo_small_primes(void) {
    /* Mathematica reference values. */
    check_eq("MultiplicativeOrder[2, 5]",  "4");
    check_eq("MultiplicativeOrder[2, 7]",  "3");
    check_eq("MultiplicativeOrder[3, 7]",  "6");      /* 3 is a PR of 7 */
    check_eq("MultiplicativeOrder[2, 11]", "10");     /* 2 is a PR of 11 */
    check_eq("MultiplicativeOrder[2, 13]", "12");
    check_eq("MultiplicativeOrder[3, 13]", "3");
    check_eq("MultiplicativeOrder[5, 17]", "16");
    check_eq("MultiplicativeOrder[2, 23]", "11");
}

static void test_mo_composite_n(void) {
    /* n = 9: phi(9) = 6.  Orders divide 6. */
    check_eq("MultiplicativeOrder[2, 9]", "6");       /* primitive */
    check_eq("MultiplicativeOrder[4, 9]", "3");
    check_eq("MultiplicativeOrder[7, 9]", "3");
    check_eq("MultiplicativeOrder[8, 9]", "2");
    /* n = 15: phi(15) = 8. */
    check_eq("MultiplicativeOrder[2, 15]", "4");
    check_eq("MultiplicativeOrder[4, 15]", "2");
    check_eq("MultiplicativeOrder[7, 15]", "4");
    check_eq("MultiplicativeOrder[11, 15]", "2");
    /* n = 16: phi(16) = 8, but the group is non-cyclic; max order is 4. */
    check_eq("MultiplicativeOrder[3, 16]", "4");
    check_eq("MultiplicativeOrder[5, 16]", "4");
    check_eq("MultiplicativeOrder[7, 16]", "2");
    check_eq("MultiplicativeOrder[9, 16]", "2");
    check_eq("MultiplicativeOrder[15, 16]", "2");
}

/* ------------------------------------------------------------------ */
/* Non-coprime k and n leave the call unevaluated.                     */
/* ------------------------------------------------------------------ */

static void test_mo_non_coprime(void) {
    /* CoprimeQ[10, 22] == False; gcd is 2. */
    check_eq("MultiplicativeOrder[10, 22]", "MultiplicativeOrder[10, 22]");
    check_eq("MultiplicativeOrder[22, 10]", "MultiplicativeOrder[22, 10]");
    check_eq("MultiplicativeOrder[6, 9]",   "MultiplicativeOrder[6, 9]");
    check_eq("MultiplicativeOrder[0, 5]",   "MultiplicativeOrder[0, 5]");
    /* n == 0 has no group at all. */
    check_eq("MultiplicativeOrder[5, 0]",   "MultiplicativeOrder[5, 0]");
}

/* ------------------------------------------------------------------ */
/* 3-arg form — discrete-log search over a residue list.               */
/* ------------------------------------------------------------------ */

static void test_mo_three_arg_spec(void) {
    /* From the docstring. */
    check_eq("MultiplicativeOrder[5, 7, {3, 11}]", "2");  /* 5^2 == 4 == 11 (mod 7) */
    check_eq("MultiplicativeOrder[5, 7, {4}]",     "2");
    check_eq("MultiplicativeOrder[5, 7, {2, 3, 4}]", "2");
}

static void test_mo_three_arg_basic(void) {
    /* {1} alone reduces to the 2-arg order. */
    check_eq("MultiplicativeOrder[5, 7, {1}]", "6");
    check_eq("MultiplicativeOrder[5, 8, {1}]", "2");
    /* Single target equal to k itself => m == 1. */
    check_eq("MultiplicativeOrder[5, 7, {5}]", "1");
    /* Target reduces mod n. */
    check_eq("MultiplicativeOrder[5, 7, {12}]", "1");     /* 12 == 5 (mod 7) */
    /* Negative target — reduced into [0, n-1]. */
    check_eq("MultiplicativeOrder[5, 7, {-2}]", "1");     /* -2 == 5 (mod 7) */
}

static void test_mo_three_arg_no_hit(void) {
    /* No power of 5 (mod 7) is even, so {2} only matches when 5^m == 2.
     * Walk: 5,4,6,2,3,1 — m == 4 hits 2. */
    check_eq("MultiplicativeOrder[5, 7, {2}]", "4");
    /* Residue not in the orbit at all: gcd(k,n)==1, but the cycle is
     * {5,4,6,2,3,1} so 0 (a non-unit) cannot be reached.  FullForm
     * print represents {x, ...} as List[x, ...]. */
    check_eq("MultiplicativeOrder[5, 7, {0}]",
             "MultiplicativeOrder[5, 7, List[0]]");
}

static void test_mo_three_arg_edge_cases(void) {
    /* Empty target list: unevaluated. */
    check_eq("MultiplicativeOrder[5, 7, {}]",
             "MultiplicativeOrder[5, 7, List[]]");
    /* Non-list third argument: unevaluated. */
    check_eq("MultiplicativeOrder[5, 7, 3]",
             "MultiplicativeOrder[5, 7, 3]");
    /* Non-integer element in the list: unevaluated. */
    check_eq("MultiplicativeOrder[5, 7, {1.5}]",
             "MultiplicativeOrder[5, 7, List[1.5]]");
    /* Non-coprime k propagates through to unevaluated even with 3 args. */
    check_eq("MultiplicativeOrder[10, 22, {1}]",
             "MultiplicativeOrder[10, 22, List[1]]");
}

/* ------------------------------------------------------------------ */
/* Diagnostics — argument count.                                       */
/* ------------------------------------------------------------------ */

static void test_mo_diagnostics(void) {
    mute_stderr_once();

    check_eq("MultiplicativeOrder[]",            "MultiplicativeOrder[]");
    check_eq("MultiplicativeOrder[1]",           "MultiplicativeOrder[1]");
    check_eq("MultiplicativeOrder[1, 2, 3, 4]",  "MultiplicativeOrder[1, 2, 3, 4]");
    /* Non-integer numerics flow through silently (no diagnostic). */
    check_eq("MultiplicativeOrder[10., 21]",     "MultiplicativeOrder[10.0, 21]");
    check_eq("MultiplicativeOrder[10., 21.]",    "MultiplicativeOrder[10.0, 21.0]");
    /* Symbolic args also stay unevaluated. */
    check_eq("MultiplicativeOrder[x, 7]",        "MultiplicativeOrder[x, 7]");
    check_eq("MultiplicativeOrder[2, y]",        "MultiplicativeOrder[2, y]");
}

/* ------------------------------------------------------------------ */
/* Bignum inputs exercise the GMP fast paths.                          */
/* ------------------------------------------------------------------ */

static void test_mo_bignum(void) {
    /* From the docstring: 7919 is prime, phi = 7918 = 2 * 37 * 107.
     * 10 mod 7919 has some order dividing 7918; 10^10000 mod 7919 has
     * the same residue class as 10^(10000 mod 7918) since 10 and 7919
     * are coprime.  The spec asserts the value is 3959. */
    check_eq("MultiplicativeOrder[10^10000, 7919]", "3959");
    /* Bignum n: phi(10^18 + 9) = 10^18 + 8 (since it's prime).
     * Order of 2 modulo (10^18 + 9) is (10^18 + 8) / 4 — i.e. 2 lies in
     * the index-4 subgroup of the unit group. */
    check_eq("MultiplicativeOrder[2, 10^18 + 9]", "250000000000000002");
    /* Sanity cross-check: 2^order == 1 (mod 10^18 + 9). */
    check_eq("PowerMod[2, MultiplicativeOrder[2, 10^18 + 9], 10^18 + 9]", "1");
    /* Pure bignum k. */
    check_eq("MultiplicativeOrder[2^100, 7]", "3");
}

/* ------------------------------------------------------------------ */
/* Structural identities cross-checked against PowerMod / EulerPhi.    */
/* ------------------------------------------------------------------ */

static void test_mo_structural(void) {
    /* k^order(k, n) == 1 (mod n). */
    check_eq("PowerMod[3, MultiplicativeOrder[3, 16], 16]", "1");
    check_eq("PowerMod[7, MultiplicativeOrder[7, 100], 100]", "1");
    check_eq("PowerMod[2, MultiplicativeOrder[2, 17], 17]", "1");

    /* For a primitive root g of n, order(g, n) == phi(n). */
    check_eq("MultiplicativeOrder[PrimitiveRoot[109], 109] == EulerPhi[109]", "True");
    check_eq("MultiplicativeOrder[PrimitiveRoot[97], 97] == EulerPhi[97]", "True");

    /* The order always divides phi(n). */
    check_eq("Mod[EulerPhi[100], MultiplicativeOrder[7, 100]]", "0");
    check_eq("Mod[EulerPhi[81], MultiplicativeOrder[5, 81]]", "0");

    /* PrimitiveRootList agreement: the order of a PR equals phi(n). */
    check_eq("MultiplicativeOrder[3, 7] == EulerPhi[7]", "True");
    check_eq("MultiplicativeOrder[5, 7] == EulerPhi[7]", "True");
    check_eq("MultiplicativeOrder[6, 11] == EulerPhi[11]", "True");
    check_eq("MultiplicativeOrder[7, 11] == EulerPhi[11]", "True");
}

static void test_mo_three_arg_consistency(void) {
    /* For {r} where r is in the orbit, k^MultiplicativeOrder[k,n,{r}] mod n
     * must equal r mod n. */
    check_eq("PowerMod[5, MultiplicativeOrder[5, 7, {4}], 7]", "4");
    check_eq("PowerMod[5, MultiplicativeOrder[5, 7, {3}], 7]", "3");
    check_eq("PowerMod[5, MultiplicativeOrder[5, 7, {6}], 7]", "6");
    /* And of course PowerMod with the result m == order matches 1. */
    check_eq("PowerMod[2, MultiplicativeOrder[2, 13, {1}], 13]", "1");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_mo_spec_examples);
    TEST(test_mo_trivial_modulus);
    TEST(test_mo_negative_n);
    TEST(test_mo_small_primes);
    TEST(test_mo_composite_n);
    TEST(test_mo_non_coprime);

    TEST(test_mo_three_arg_spec);
    TEST(test_mo_three_arg_basic);
    TEST(test_mo_three_arg_no_hit);
    TEST(test_mo_three_arg_edge_cases);

    TEST(test_mo_diagnostics);
    TEST(test_mo_bignum);

    TEST(test_mo_structural);
    TEST(test_mo_three_arg_consistency);

    printf("All MultiplicativeOrder tests passed!\n");
    return 0;
}
