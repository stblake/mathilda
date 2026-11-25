/* Builtin-audit regressions: power series (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* Series#1 (critical): Series[x, {x, 0, 2^62}] sized its dense coefficient
 * buffer unchecked and segfaulted.  An order no machine can hold aborts with
 * General::nomem; ordinary orders are unchanged. */
static void test_series_huge_order(void) {
    audit_chk("Series#1", "Quiet[CheckAbort[Series[x, {x, 0, 2^62}], aborted]]", "aborted");
    audit_chk("Series#1",
              "Quiet[Check[CheckAbort[Series[Sin[x], {x, 0, 2^62}], 0], fired, General::nomem]]",
              "fired");
    audit_chk("Series#1", "Quiet[CheckAbort[Series[x, {x, 0, 9223372036854775807}], aborted]]",
              "aborted");
    audit_true("Series#1", "Normal[Series[Sin[x], {x, 0, 5}]] === x - x^3/6 + x^5/120");
    /* A constant is returned verbatim whatever the order (no buffer). */
    audit_chk("Series#1", "Series[y, {x, 0, 2^62}]", "y");
}

/* SeriesData#1 (high): a one-argument numeric function of a SeriesData
 * composes (only Power/Exp did). Expected values are the Taylor series. */
static void test_seriesdata_function_composition(void) {
    audit_true("SeriesData#1",
               "Sin[Series[x, {x, 0, 3}]] === SeriesData[x, 0, {1, 0, -1/6}, 1, 4, 1]");
    audit_true("SeriesData#1",
               "Cos[Series[x, {x, 0, 3}]] === SeriesData[x, 0, {1, 0, -1/2, 0}, 0, 4, 1]");
    audit_true("SeriesData#1",
               "Log[Series[1 + x, {x, 0, 3}]] === SeriesData[x, 0, {1, -1/2, 1/3}, 1, 4, 1]");
    audit_true("SeriesData#1",
               "Sinh[Series[x, {x, 0, 3}]] === SeriesData[x, 0, {1, 0, 1/6}, 1, 4, 1]");
    audit_true("SeriesData#1",
               "ArcTan[Series[x, {x, 0, 5}]] === SeriesData[x, 0, {1, 0, -1/3, 0, 1/5}, 1, 6, 1]");
    audit_true("SeriesData#1",
               "Normal[Sin[Series[x, {x, 0, 3}]]] === Normal[Series[Sin[x], {x, 0, 3}]]");
    audit_true("SeriesData#1",
               "Exp[Series[x, {x, 0, 3}]] === SeriesData[x, 0, {1, 1, 1/2, 1/6}, 0, 4, 1]");
}

/* SeriesData#2 (high): Exp of a series u with nmin not dividing the order took
 * floor(order/nmin) kernel terms, one short: Exp[x^2 + O[x]^3] lost its x^2,
 * and so did base^x = Exp[x Log[base]]. (1+x)^x = Exp[x^2 - x^3/2 + x^4/3 ...]
 * = 1 + x^2 - x^3/2 + 5/6 x^4 + O[x]^5. */
static void test_seriesdata_exp_kernel_terms(void) {
    audit_true("SeriesData#2",
               "SeriesData[x, 0, {1, 1, 0}, 0, 3, 1]^x === SeriesData[x, 0, {1, 0, 1}, 0, 3, 1]");
    audit_true("SeriesData#2",
               "SeriesData[x, 0, {1, 1, 0}, 0, 3, 1]^x === Series[(1 + x)^x, {x, 0, 2}]");
    audit_true("SeriesData#2",
               "SeriesData[x, 0, {1, 1, 0, 0, 0}, 0, 5, 1]^x === "
               "SeriesData[x, 0, {1, 0, 1, -1/2, 5/6}, 0, 5, 1]");
    audit_true("SeriesData#2",
               "Exp[SeriesData[x, 0, {1}, 2, 3, 1]] === SeriesData[x, 0, {1, 0, 1}, 0, 3, 1]");
    audit_true("SeriesData#2",
               "SeriesData[x, 0, {1, 1, 0, 0}, 0, 4, 1]^x === SeriesData[x, 0, {1, 0, 1, -1/2}, 0, 4, 1]");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* ArcCsch#1: ArcCsch[0] */
    audit_true("ArcCsch#1", "ArcCsch[0] === ComplexInfinity");
    /* ArcSech#2: ArcSech[0] */
    audit_true("ArcSech#2", "ArcSech[0] === Infinity");
    audit_true("ArcSech#2", "ArcSech[1] === 0");
    /* ArcSinh#1: ArcSinh along the imaginary axis */
    audit_true("ArcSinh#1", "ArcSinh[I Infinity] === Infinity");
    audit_true("ArcSinh#1", "ArcSinh[-I Infinity] === -Infinity");
    audit_true("ArcSinh#1", "Re[ArcSinh[10.^8 I]] > 19");
    /* ArcCosh#3: ArcCosh of -Infinity and I Infinity */
    audit_true("ArcCosh#3", "ArcCosh[-Infinity] === Infinity");
    audit_true("ArcCosh#3", "ArcCosh[I Infinity] === Infinity");
    audit_true("ArcCosh#3", "Re[ArcCosh[-10.^8]] > 19");
    /* ArcCsch#3: ArcCsch[Indeterminate] */
    audit_true("ArcCsch#3", "ArcCsch[Indeterminate] === Indeterminate");
    /* ArcSinh#2: ArcSinh[Indeterminate] */
    audit_true("ArcSinh#2", "ArcSinh[Indeterminate] === Indeterminate");
    /* SeriesData#4: a series about Infinity is in powers of 1/x */
    audit_true("SeriesData#4", "Check[Normal[SeriesData[x, Infinity, {1, 2}, 1, 3, 1]], \"m\"] === 1/x + 2/x^2");
    audit_true("SeriesData#4", "Normal[SeriesData[x, Infinity, {3}, 0, 2, 1]] === 3");
}

/* Wave 4 (G5): Series of relations, nested Normal, multivariate
 * SeriesCoefficient, ArcCosh on [-1, 1). */
static void test_wave4_series(void) {
    audit_chk("Series#2", "Series[x == Sin[x], {x, 0, 3}]", "x + O[x]^4 == x - 1/6 x^3 + O[x]^4");
    audit_chk("Normal#1", "Normal[Series[Exp[x y], {x, 0, 2}, {y, 0, 2}]]", "1 + x y + 1/2 x^2 y^2");
    audit_chk("SeriesCoefficient#1", "SeriesCoefficient[Exp[x y], {x, 0, 2}, {y, 0, 2}]", "1/2");
    audit_chk("ArcCsch#2", "With[{F = Integrate[ArcCsch[x], x]}, Chop[N[(D[F, x] - ArcCsch[x]) /. x -> -2]]]", "0");
    audit_chk("ArcCsch#2", "With[{F = Integrate[x ArcCsch[x], x]}, Chop[N[(D[F, x] - x ArcCsch[x]) /. x -> -2]]]", "0");
    audit_chk("ArcCosh#1", "ArcCosh[0]", "(1/2*I) Pi");
    audit_chk("ArcCosh#1", "ArcCosh[-1]", "I Pi");
    audit_chk("ArcCosh#1", "ArcCosh[1/2]", "(1/3*I) Pi");
}

int main(void) {
    audit_init();
    TEST(test_series_huge_order);
    TEST(test_seriesdata_function_composition);
    TEST(test_seriesdata_exp_kernel_terms);
    TEST(test_s11_special_values);
    TEST(test_wave4_series);
    return audit_finish("audit_power_series");
}
