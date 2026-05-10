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
    /* -I (negative imaginary). */
    assert_eval_eq("(-I*Pi)^(1/2)",
                   "Sqrt[Pi]/(-1)^(1/4)", 0);
    /* Cube root case. */
    assert_eval_eq("(I*Pi)^(1/3)",
                   "(-1)^(1/6) Pi^(1/3)", 0);
    /* Negative branch: x of unknown sign -- must NOT distribute. */
    assert_eval_eq("Sqrt[I*x]",
                   "Sqrt[(I) x]", 0);
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

int main(void) {
    symtab_init();
    core_init();

    TEST(test_bug3_imag_times_positive);
    TEST(test_bug4_negative_integer_part);
    TEST(test_bug5_zero_to_zero);

    printf("All power_corpus tests passed!\n");
    return 0;
}
