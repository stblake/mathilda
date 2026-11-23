/* Builtin-audit regressions: control flow (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* For#1 (high, hang): "Throw in For's start is swallowed, so
 * Catch[For[Throw[7], True, i++]] loops forever".  The start's value was
 * freed with no flow-control check.  The same holds for the test (For#2) and
 * for the evaluated iterator bounds of Do / Table / Sum / Product. */
static void test_for_start_throw(void) {
    audit_run("ClearAll[fi]");
    audit_chk("For#1", "Catch[For[Throw[7], True, fi++]]", "7");
    audit_chk("For#1", "Catch[For[fi = 0, Throw[5], fi++]]", "5");
    audit_chk("For#1", "Catch[For[fi = 0, If[fi == 2, Throw[fi], fi < 5], fi++]]", "2");
    audit_chk("For#1", "Quiet[CheckAbort[For[Abort[], True, fi++], aborted]]", "aborted");
    audit_chk("For#1", "Catch[Do[1, {fi, Throw[3]}]]", "3");
    audit_chk("For#1", "Catch[Table[1, {fi, 1, Throw[4]}]]", "4");
    audit_chk("For#1", "Catch[Sum[fi, {fi, Throw[6]}]]", "6");
    audit_chk("For#1", "Catch[Product[fi, {fi, Throw[8]}]]", "8");
    /* sanity: ordinary loops */
    audit_chk("For#1", "Module[{s = 0}, For[fi = 1, fi <= 4, fi++, s += fi]; s]", "10");
    audit_chk("For#1", "For[fi = 0, fi < 3, fi++]", "Null");
    audit_chk("For#1", "f0[] := For[fi = 0, True, fi++, If[fi == 3, Return[fi]]]; f0[]", "3");
    audit_chk("For#1", "Table[fi, {fi, 3}]", "{1, 2, 3}");
}

int main(void) {
    audit_init();
    TEST(test_for_start_throw);
    return audit_finish("audit_control_flow");
}
