/* Builtin-audit regressions: machine learning (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* Classify#1 (high, crash): "Classify[{Rule[]}] segfaults: the first rule's
 * lhs is read without checking its argument count".  Malformed training data
 * leaves the call unevaluated; Predict's reader had the same read. */
static void test_classify_malformed_rules(void) {
    audit_chk("Classify#1", "Classify[{Rule[]}]", "Classify[{Rule[]}]");
    audit_chk("Classify#1", "Classify[{Rule[1]}]", "Classify[{Rule[1]}]");
    audit_chk("Classify#1", "Classify[{Rule[1, 2, 3]}] === Unevaluated[Classify[{Rule[1, 2, 3]}]]", "True");
    audit_chk("Classify#1", "Classify[{RuleDelayed[]}]", "Classify[{RuleDelayed[]}]");
    audit_chk("Classify#1", "Classify[{}]", "Classify[{}]");
    audit_chk("Classify#1", "Predict[{Rule[]}]", "Predict[{Rule[]}]");
    audit_chk("Classify#1", "Predict[{Rule[1]}]", "Predict[{Rule[1]}]");
    /* sanity: well-formed data still trains */
    audit_chk("Classify#1",
              "Classify[{1. -> \"a\", 2. -> \"a\", 10. -> \"b\", 11. -> \"b\"}][10.5]", "\"b\"");
}

int main(void) {
    audit_init();
    TEST(test_classify_malformed_rules);
    return audit_finish("audit_machine_learning");
}
