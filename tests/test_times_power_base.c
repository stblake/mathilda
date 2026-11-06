/* test_times_power_base.c
 *
 * Times merges a factor with a power of that same factor when the factor is
 * itself a power, as the reference does: m^2 (m^2)^(-1 + D/2) is
 * (m^2)^(D/2) (FeynCalc's FCFeynmanParametrize-ID13).  The base grouping
 * keyed m^2 as (m, 2) and (m^2)^a as (m^2, a), so the two never met.  Only
 * exact equality merges.  An inexact zero absorbs symbolic factors
 * (0. d is 0.).  Expected values computed by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_power_base_merge(void) {
    assert_eval_eq("m^2 (m^2)^(-1 + D/2)", "(m^2)^(1/2 D)", 0);
    assert_eval_eq("(m^2)^a m^2", "(m^2)^(1 + a)", 0);
    assert_eval_eq("Sqrt[m^2] m^2", "(m^2)^(3/2)", 0);
    assert_eval_eq("x^2 (x^2)^a (x^2)^b", "(x^2)^(1 + a + b)", 0);
    /* not the same base: no merge */
    assert_eval_eq("m^4 Sqrt[m^2]", "m^4 Sqrt[m^2]", 0);
    assert_eval_eq("x^2 (x^2)^a x", "x^3 (x^2)^a", 0);
    assert_eval_eq("m^2 m^3", "m^5", 0);
}

static void test_inexact_zero_absorbs(void) {
    assert_eval_eq("0. d", "0.0", 0);
    assert_eval_eq("0. (d y)", "0.0", 0);
    assert_eval_eq("0. Sin[x]", "0.0", 0);
    assert_eval_eq("Quiet[0. Infinity]", "Indeterminate", 0);
}

/* A product whose collected form exposes a new merge is not final: the
 * collector's own output (a^-1 a, from (a^-1)^(1/3) (a^-1)^(2/3) a) is
 * re-collected.  Exponents of one base add: (1/a)^(1/3) (1/a)^(2/3) = 1/a. */
static void test_collected_product_is_recollected(void) {
    assert_eval_eq("Times[-1, Power[Power[ca, -1], 1/3], Power[Power[ca, -1], 2/3], ca]", "-1", 0);
    assert_eval_eq("Power[Power[cb, -1], 1/2] Power[Power[cb, -1], 1/2] cb^2", "cb", 0);
    /* the folded coefficient and untouched factors: already final */
    assert_eval_eq("-1/(3 Pi)", "-1/(3 Pi)", 0);
    assert_eval_eq("FullForm[-1/(3 Pi)]", "Times[Rational[-1, 3], Power[Pi, -1]]", 0);
    assert_eval_eq("-Cos[6/5] Sin[6/5]^4", "-Cos[6/5] Sin[6/5]^4", 0);
    /* Sin[k Pi/2]/(k Pi): Sin = 1, 0, -1, 0 at k = 1..4 */
    assert_eval_eq("Table[Sin[Pi k/2]/(Pi k), {k, 4}]", "{1/Pi, 0, -1/(3 Pi), 0}", 0);
    assert_eval_eq("Table[1/(i + j - 1), {i, 3}, {j, 3}]",
                   "{{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}", 0);
}

/* Exact powers with a machine-integer exponent (the fast path ahead of the
 * shape tests).  Values by hand: 3^40 = 12157665459056928801,
 * 5^40 = 9094947017729282379150390625, 2^70 = 1180591620717411303424. */
static void test_exact_integer_powers(void) {
    assert_eval_eq("(-3)^-3", "-1/27", 0);
    assert_eval_eq("(2/3)^-2", "9/4", 0);
    assert_eval_eq("(-2/3)^3", "-8/27", 0);
    assert_eval_eq("(-1)^-3", "-1", 0);
    assert_eval_eq("1^-5", "1", 0);
    assert_eval_eq("2^-70", "1/1180591620717411303424", 0);
    assert_eval_eq("(3/5)^40", "12157665459056928801/9094947017729282379150390625", 0);
    assert_eval_eq("(5/3)^-40", "12157665459056928801/9094947017729282379150390625", 0);
    /* a bignum Rational base: the reduction is skipped only for a coprime base */
    assert_eval_eq("(2^70/3)^2", "1393796574908163946345982392040522594123776/9", 0);
    assert_eval_eq("Quiet[0^-1]", "ComplexInfinity", 0);
    assert_eval_eq("Rational[6, 4]", "3/2", 0);
    assert_eval_eq("Rational[2^40 3, 2^35 9]", "32/3", 0);
}

/* Machine-integer comparisons, canonical order of calls with different heads
 * (alphabetical by head, case-blind), and trig ratio canonicalisation. */
static void test_compare_order_trig(void) {
    assert_eval_eq("{3 < 5, 5 == 5, 7 > 9, -4 <= -4, 2^62 > -2^62, 4 != 4}",
                   "{True, True, False, True, True, False}", 0);
    assert_eval_eq("Sort[{Sin[x], f[x], Cos[x], h[x, y], g[x], h[x]}]",
                   "{Cos[x], f[x], g[x], h[x], h[x, y], Sin[x]}", 0);
    assert_eval_eq("Tan[x] Cos[x]", "Sin[x]", 0);
    assert_eval_eq("Sin[x]^2 Csc[x]", "Sin[x]", 0);
    assert_eval_eq("Sec[x] Sin[x]", "Tan[x]", 0);
    assert_eval_eq("Cos[x] Sin[x]^3", "Cos[x] Sin[x]^3", 0);
    assert_eval_eq("Cosh[y]^2 Sinh[y]", "Cosh[y]^2 Sinh[y]", 0);
    assert_eval_eq("Sin[ArcSin[z]] + Cos[ArcCos[z]] + Sinh[ArcSinh[z]]", "3 z", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_power_base_merge);
    TEST(test_inexact_zero_absorbs);
    TEST(test_collected_product_is_recollected);
    TEST(test_exact_integer_powers);
    TEST(test_compare_order_trig);
    printf("All Times power-base tests passed.\n");
    return 0;
}
