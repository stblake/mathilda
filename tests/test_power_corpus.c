/*
 * test_power_corpus.c -- Regression tests for Power algorithmic fixes.
 *
 * Bug 3 -- (I*Pi)^(1/2) -> (-1)^(1/4) Sqrt[Pi]
 *
 *   builtin_power now distributes a rational exponent over a Times-base
 *   when the base contains a pure-imaginary factor (Complex[0, k]) and
 *   every other factor is known strictly positive (Pi, E, positive
 *   integer/rational). The branch falls through unchanged when the
 *   factor list is mixed (Sqrt[I x] etc.) so we don't cross a branch
 *   cut for unknown-sign inputs.
 *
 * Bug 4 -- Power[2, -5/3] -> 1/(2 * 2^(2/3))
 *
 *   The integer-part extraction in builtin_power had a guard `p >= q`
 *   which excluded negative `p`. We relax to `|p| >= q` so the C99
 *   truncating divide produces a_int = -1, b_rem = -2 for p=-5, q=3,
 *   yielding the canonical Times[Rational[1,2], Power[2, -2/3]].
 *
 *   The Times-side canonicalizer (src/times.c radical-canonical pass)
 *   was simultaneously updated to NOT pull factors of `b` from den
 *   into the exponent past the proper-fraction boundary; otherwise the
 *   1/2 * 2^(-2/3) form would round-trip back into Power[2, -5/3] and
 *   loop. This test pins the new fixed-point.
 *
 * Bug 5 -- Power[0, 0] -> Indeterminate
 *
 *   Previously returned 1 because Power[_, 0] -> 1 fired before the
 *   0^anything dispatch. Now 0^0 is intercepted first and emits the
 *   `Power::indet` warning.
 */

#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Bug 3: pure-imaginary times positive base, rational exponent. */
static void test_bug3_imag_times_positive(void) {
    assert_eval_eq("(I*Pi)^(1/2)",
                   "(-1)^(1/4) Sqrt[Pi]", 0);
    /* Sqrt[I*Pi] = (I*Pi)^(1/2). */
    assert_eval_eq("Sqrt[I*Pi]",
                   "(-1)^(1/4) Sqrt[Pi]", 0);
    /* Multiple positive factors. */
    assert_eval_eq("(I*Pi*2)^(1/2)",
                   "(-1)^(1/4) Sqrt[2] Sqrt[Pi]", 0);
    /* -I (negative imaginary). Power[-1, -1/4] canonicalises (floor of the
     * exponent into [0,1)) to -(-1)^(3/4), matching Mathematica:
     * Sqrt[Pi]/(-1)^(1/4) == (-1)^(-1/4) Sqrt[Pi] == -(-1)^(3/4) Sqrt[Pi]. */
    assert_eval_eq("(-I*Pi)^(1/2)",
                   "-(-1)^(3/4) Sqrt[Pi]", 0);
    /* Cube root case. */
    assert_eval_eq("(I*Pi)^(1/3)",
                   "(-1)^(1/6) Pi^(1/3)", 0);
    /* Negative branch: x of unknown sign -- must NOT distribute. */
    assert_eval_eq("Sqrt[I*x]",
                   "Sqrt[I x]", 0);
    /* No imaginary factor: leave as-is. */
    assert_eval_eq("Sqrt[2*Pi]",
                   "Sqrt[2 Pi]", 0);
}

/* Bug 4: integer-part extraction for negative rational exponents. */
static void test_bug4_negative_integer_part(void) {
    /* Print form: 1/2/2^(2/3) (Times-print of Rational[1,2] * Power[2,-2/3]). */
    assert_eval_eq("Power[2, -5/3]",       "1/2/2^(2/3)", 0);
    /* FullForm pinpoints the Times[Rational[1,2], Power[2, Rational[-2,3]]]. */
    assert_eval_eq("FullForm[Power[2, -5/3]]",
                   "Times[Rational[1, 2], Power[2, Rational[-2, 3]]]", 0);
    /* Negative-half exponent: 2^(-3/2) = 1/2 * 2^(-1/2). */
    assert_eval_eq("Power[2, -3/2]", "1/2/Sqrt[2]", 0);
    /* Sqrt[2]/4 must canonicalise to the same proper-fraction form. */
    assert_eval_eq("Sqrt[2]/4", "1/2/Sqrt[2]", 0);
    /* Positive side still works (regression). */
    assert_eval_eq("Power[2, 5/3]", "2 2^(2/3)", 0);
    /* Proper-fraction stays as-is. */
    assert_eval_eq("Power[2, -1/3]", "1/2^(1/3)", 0);
    /* Pulling factors from numerator coefficient still works. */
    assert_eval_eq("8 * Power[2, -1/3]", "4 2^(2/3)", 0);
}

/* Bug 5: 0^0 is Indeterminate. */
static void test_bug5_zero_to_zero(void) {
    /* 0^0 must produce the symbol Indeterminate. */
    assert_eval_eq("Power[0, 0]", "Indeterminate", 0);
    assert_eval_eq("Power[0., 0]", "Indeterminate", 0);
    assert_eval_eq("Power[0, 0.]", "Indeterminate", 0);
    /* 0^positive integer still 0. */
    assert_eval_eq("Power[0, 3]", "0", 0);
    /* Positive base ^ 0 still 1 (regression). */
    assert_eval_eq("Power[5, 0]", "1", 0);
    /* 0^negative still ComplexInfinity (regression). */
    assert_eval_eq("Power[0, -1]", "ComplexInfinity", 0);
}

/* Rational-base + rational-exponent canonicalisation (2026-05-12).
 * Powers of Rational[n, d] with non-integer exponent now distribute
 * to Power[n, p/q] * Power[d, -p/q] when at least one piece has a
 * perfect q-th-power extraction (or |n| == 1). Cases without any
 * extraction (Sqrt[2/3], (4/9)^(2/3)) stay unevaluated. */
static void test_rational_base_rational_exp(void) {
    /* Numerator-only trigger (|n| == 1). */
    assert_eval_eq("Power[1/54, 2/3]", "1/9/2^(2/3)", 0);
    assert_eval_eq("Power[1/8, 2/3]",  "1/4", 0);
    assert_eval_eq("Power[1/4, 1/2]",  "1/2", 0);

    /* Both numerator and denominator extract. */
    assert_eval_eq("Power[8/27, 2/3]", "4/9", 0);
    assert_eval_eq("Power[4/9, 1/2]",  "2/3", 0);

    /* Numerator-only perfect extraction with non-perfect denominator. */
    assert_eval_eq("Power[8/3, 1/2]", "2 Sqrt[2/3]", 0);

    /* No extraction possible -- left unevaluated. */
    assert_eval_eq("Power[2/3, 1/2]",  "Sqrt[2/3]", 0);
    assert_eval_eq("Power[4/9, 2/3]",  "(4/9)^(2/3)", 0);

    /* Negative numerator routes through the integer negative-base path. */
    assert_eval_eq("Power[-1/8, 1/3]", "1/2 (-1)^(1/3)", 0);
    assert_eval_eq("Power[-1/4, 1/2]", "1/2*I", 0);
}

/* Power[Times[positive_factors], p/q] distribution (2026-05-12).
 * Distributes when every factor is known positive AND at least one
 * factor cleanly reduces to a rational coefficient under Power[_, p/q]. */
static void test_power_of_times_positives(void) {
    /* Integer perfect-square factor extracts cleanly. */
    assert_eval_eq("Sqrt[4 Pi]",   "2 Sqrt[Pi]", 0);
    assert_eval_eq("(Pi/8)^(2/3)", "1/4 Pi^(2/3)", 0);

    /* Conservative when no factor cleanly reduces. */
    assert_eval_eq("Sqrt[2 Pi]",       "Sqrt[2 Pi]", 0);
    assert_eval_eq("Power[4 Pi, 2/3]", "(4 Pi)^(2/3)", 0);
    assert_eval_eq("Sqrt[2 Sqrt[3]]",  "Sqrt[2 Sqrt[3]]", 0);

    /* Chains with rational-base distribution: Sqrt[(1/54)^(2/3)] cascades
     * (1/54)^(2/3) -> 1/9 * 2^(-2/3), then Sqrt of that distributes. */
    assert_eval_eq("Sqrt[Power[1/54, 2/3]]", "1/3/2^(1/3)", 0);

    /* Nested clean extraction. */
    assert_eval_eq("Sqrt[4 Sqrt[3]]", "2 3^(1/4)", 0);
}

/* Power[Times[positives], p/q] distribution when the trigger is a
 * Power[positive, p'/q'] factor that composes cleanly (composed
 * reduced denominator <= q'), not a perfect-power rational factor
 * (2026-05-13).  Closes Sqrt[k/2^(2/3)] -> Sqrt[k]/2^(1/3) for
 * arbitrary positive k. */
static void test_power_of_times_power_factor_composes(void) {
    /* Power factor triggers distribution; rational coefficient need
     * not factor cleanly under the outer p/q. */
    assert_eval_eq("Sqrt[3/2^(2/3)]",      "Sqrt[3]/2^(1/3)",       0);
    assert_eval_eq("Sqrt[1/3/2^(2/3)]",    "1/(2^(1/3) Sqrt[3])",   0);
    assert_eval_eq("Sqrt[Pi/2^(2/3)]",     "Sqrt[Pi]/2^(1/3)",      0);

    /* Pre-existing perfect-square cases still work. */
    assert_eval_eq("Sqrt[1/2^(2/3)]",      "1/2^(1/3)",             0);
    assert_eval_eq("Sqrt[9/2^(2/3)]",      "3/2^(1/3)",             0);

    /* Adjacency cases must remain unchanged: composition would
     * introduce a new (larger-denominator) radical. */
    assert_eval_eq("Sqrt[2 Sqrt[3]]",      "Sqrt[2 Sqrt[3]]",       0);
    assert_eval_eq("Sqrt[2 Pi]",           "Sqrt[2 Pi]",            0);
    assert_eval_eq("Power[4 Pi, 2/3]",     "(4 Pi)^(2/3)",          0);
}

/* Power-of-Power composition with positive base + rational outer (2026-05-12). */
static void test_power_of_power_positive_base(void) {
    /* Positive base with a *real* (rational) inner exponent merges:
     * (2^(2/3))^(1/2) = 2^(1/3), since 2^(2/3) is a positive real. */
    assert_eval_eq("Sqrt[Power[2, 2/3]]", "2^(1/3)", 0);

    /* Positive base but *symbolic* inner exponent must NOT merge: a may be
     * complex, and Sqrt[2^a] != 2^(a/2) on the principal branch (e.g.
     * a = 3 Pi I / Log[2] gives Sqrt[2^a] = I but 2^(a/2) = -I). Matches
     * Mathematica, which leaves Sqrt[2^a] unevaluated -- PowerExpand is the
     * documented route to the merged form. */
    assert_eval_eq("Sqrt[Power[2, a]]", "Sqrt[2^a]", 0);

    /* Pi positive -- Sqrt[Pi^2] reduces to Pi (not |Pi|). */
    assert_eval_eq("Sqrt[Power[Pi, 2]]", "Pi", 0);

    /* Unknown-sign symbol -- must NOT reduce (Sqrt[x^4] != x^2 for x<0). */
    assert_eval_eq("Sqrt[Power[x, 4]]", "Sqrt[x^4]", 0);

    /* Integer outer exponent on negative base still works (pre-existing rule). */
    assert_eval_eq("Power[Power[-1, 2], 3]", "1", 0);

    /* (E^x)^y with symbolic exponents must stay nested -- x may be complex,
     * so the core evaluator must not fold it to E^(x y). PowerExpand handles
     * the merge when the user asserts it is safe. */
    assert_eval_eq("(E^x)^y", "E^x^y", 0);
    assert_eval_eq("(E^x)^y", "Power[Power[E, x], y]", 1);
}

/* Perfect-power extraction on BigInt bases: previously skipped because
 * the gate insisted on EXPR_INTEGER. */
static void test_bigint_perfect_root(void) {
    assert_eval_eq("Sqrt[10^20]",  "10000000000", 0);
    assert_eval_eq("Sqrt[10^100]", "100000000000000000000000000000000000000000000000000", 0);
    assert_eval_eq("Sqrt[2^64]",   "4294967296", 0);
    assert_eval_eq("Sqrt[2^100]",  "1125899906842624", 0);
    assert_eval_eq("Power[10^60, 1/3]",  "100000000000000000000", 0);
    assert_eval_eq("Power[10^100, 1/2]", "100000000000000000000000000000000000000000000000000", 0);
    assert_eval_eq("Power[2^100, 1/4]",  "33554432", 0);
    assert_eval_eq("Power[10^100, 1/5]", "100000000000000000000", 0);
    assert_eval_eq("Power[10^21, 1/3]",  "10000000", 0);
    /* Mixed exponent: 10^20 = 10^(20), so (10^20)^(2/3) = 10^(40/3) which
     * factors as 10^13 * 10^(1/3). */
    assert_eval_eq("Power[10^20, 2/3]",  "10000000000000 10^(1/3)", 0);

    /* Negative BigInt base with q == 2: Sqrt[-bigint] = I * Sqrt[|bigint|]
     * via the I^p * |n|^(p/2) identity. Previously left unreduced because
     * the negative-base path gated on EXPR_INTEGER only. */
    assert_eval_eq("Sqrt[-(10^100)]", "100000000000000000000000000000000000000000000000000*I", 0);
    assert_eval_eq("Sqrt[-(10^50)]",  "10000000000000000000000000*I", 0);
    assert_eval_eq("Sqrt[-(2^64)]",   "4294967296*I", 0);
    /* (-bigint)^(3/2): I^3 = -I. */
    assert_eval_eq("Power[-(10^100), 3/2]",
        "-1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*I", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_bug3_imag_times_positive);
    TEST(test_bug4_negative_integer_part);
    TEST(test_bug5_zero_to_zero);
    TEST(test_rational_base_rational_exp);
    TEST(test_power_of_times_positives);
    TEST(test_power_of_power_positive_base);
    TEST(test_power_of_times_power_factor_composes);
    TEST(test_bigint_perfect_root);

    printf("All power_corpus tests passed!\n");
    return 0;
}
