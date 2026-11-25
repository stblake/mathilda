/* Builtin-audit regressions: comparisons (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* Inequality#1 (high): "Residual chain malformed (even arg count, missing
 * operator) when a decided middle pair is dropped". */
static void test_inequality_residual_wellformed(void) {
    audit_run("ClearAll[x, y, z, r]");
    audit_chk("Inequality#1", "Inequality[x, Less, 1, Less, 2, Less, y]", "x < 1 && 2 < y");
    audit_chk("Inequality#1", "x < 1 < 2 < y", "x < 1 && 2 < y");
    audit_chk("Inequality#1", "Inequality[x, Less, 1, Less, 2, Less, y] /. {x -> 0, y -> 3}", "True");
    audit_chk("Inequality#1", "Inequality[x, Less, 1, Less, 2, Less, y] /. {x -> 5, y -> 3}", "False");
    audit_chk("Inequality#1", "Inequality[x, Less, 1, Unequal, 2, Less, y]", "x < 1 && 2 < y");
    /* sanity: adjacent undecided pairs stay one chain */
    audit_chk("Inequality#1", "x < y < 2 < z", "x < y < 2 < z");
    audit_chk("Inequality#1", "x < y < z", "x < y < z");
    audit_chk("Inequality#1", "1 < 2 < 3", "True");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* SameQ#2: DirectedInfinity[d] evaluates to the canonical infinity */
    audit_true("SameQ#2", "(DirectedInfinity[1] === Infinity) === True");
    audit_true("SameQ#2", "(DirectedInfinity[] === ComplexInfinity) === True");
    audit_true("SameQ#2", "(DirectedInfinity[-1] === -Infinity) === True");
    audit_true("SameQ#2", "(DirectedInfinity[2 I] === I Infinity) === True");
    audit_true("SameQ#2", "DirectedInfinity[1] == Infinity");
    /* UnsameQ#2: Infinity =!= DirectedInfinity[1] is False */
    audit_true("UnsameQ#2", "(Infinity =!= DirectedInfinity[1]) === False");
}

static void test_wave4_variadic_inequalities(void) {
    /* Variadic inequalities decide each adjacent pair; a False pair wins,
     * decided True pairs drop out (src/comparisons.c evaluate_inequality). */
    audit_chk("GreaterEqual#1", "GreaterEqual[x,4,5]", "False");
    audit_chk("GreaterEqual#1", "GreaterEqual[x,y,3,4]", "False");
    audit_chk("LessEqual#1", "LessEqual[x, 3, 2]", "False");
    audit_chk("LessEqual#1", "LessEqual[1, 2, x, 3]", "2 <= x <= 3");
    audit_chk("LessEqual#1", "Less[x, 1, 2, 3, y]", "x < 1 < 3 < y");
}

/* Wave 5 G7 (arith-misc): low-severity findings. Message counts are read
 * off $MessageList inside a Block (each diagnostic once per evaluation). */
static void test_wave5_g7_comparisons(void) {
    audit_chk("Equal#2", "True == False", "False");
    audit_chk("Equal#2", "True == True", "True");
    audit_chk("Unequal#1", "True != False", "True");
    audit_chk("Greater#4", "Quiet[Check[1 + I > 0, fired, Greater::nord]]", "fired");
    audit_chk("Greater#4", "Quiet[1 + I > 0]", "1 + I > 0");
    audit_chk("GreaterEqual#2", "Quiet[Check[I >= 0, fired, GreaterEqual::nord]]", "fired");
    audit_chk("Less#3", "Quiet[Check[I < 1, fired, Less::nord]]", "fired");
    audit_chk("Less#3", "Quiet[I < 1]", "I < 1");
    audit_chk("LessEqual#4", "Quiet[Check[I <= 2, fired, LessEqual::nord]]", "fired");
    audit_chk("LessEqual#3", "Interval[{1, 2}] <= 2", "True");
    audit_chk("LessEqual#3", "Interval[{1, 2}] <= 3", "True");
    audit_chk("LessEqual#3", "2 >= Interval[{1, 2}]", "True");
    audit_chk("LessEqual#3", "Interval[{1, 2}] >= 2", "Interval[{1, 2}] >= 2");
}

int main(void) {
    audit_init();
    TEST(test_inequality_residual_wellformed);
    TEST(test_s11_special_values);
    TEST(test_wave4_variadic_inequalities);
    TEST(test_wave5_g7_comparisons);
    return audit_finish("audit_comparisons");
}
