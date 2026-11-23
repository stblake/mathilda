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

int main(void) {
    audit_init();
    TEST(test_inequality_residual_wellformed);
    return audit_finish("audit_comparisons");
}
