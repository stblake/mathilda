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

int main(void) {
    audit_init();
    TEST(test_series_huge_order);
    TEST(test_seriesdata_function_composition);
    TEST(test_seriesdata_exp_kernel_terms);
    return audit_finish("audit_power_series");
}
