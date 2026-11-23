/* Builtin-audit regressions: pattern matching (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* PatternSequence#1 (high): "Condition on a PatternSequence never matches". */
static void test_condition_on_patternsequence(void) {
    audit_run("ClearAll[psh, x, y, s]");
    audit_chk("PatternSequence#1", "MatchQ[f[1, 2], f[PatternSequence[x_, y_] /; x < y]]", "True");
    audit_chk("PatternSequence#1", "MatchQ[f[2, 1], f[PatternSequence[x_, y_] /; x < y]]", "False");
    audit_run("psh[s : PatternSequence[x_, y_] /; x < y] := {s}");
    audit_chk("PatternSequence#1", "{psh[1, 2], psh[2, 1]}", "{{1, 2}, psh[2, 1]}");
    audit_chk("PatternSequence#1", "f[1, 2, 3] /. f[PatternSequence[x_, y_] /; x > 1, z_] :> {x, y, z}", "f[1, 2, 3]");
    audit_chk("PatternSequence#1", "f[3, 2, 1] /. f[PatternSequence[x_, y_] /; x > 1, z_] :> {x, y, z}", "{3, 2, 1}");
    /* sanity */
    audit_chk("PatternSequence#1", "MatchQ[f[1, 2], f[PatternSequence[x_, y_]]]", "True");
    audit_chk("PatternSequence#1", "MatchQ[f[1, 2], f[x__ /; Length[{x}] == 2]]", "True");
}

int main(void) {
    audit_init();
    TEST(test_condition_on_patternsequence);
    return audit_finish("audit_pattern_matching");
}
