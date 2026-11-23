/* Builtin-audit regressions: expression information (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* ByteCount#1 (high): "ByteCount segfaults on a deeply nested expression
 * (unbounded C recursion)".  The walk is iterative, and its total is the sum
 * over nodes, so f-nesting adds a fixed amount per level. */
static void test_bytecount_deep_nesting(void) {
    audit_run("auditDeep = Nest[f, x, 200000]");
    audit_true("ByteCount#1", "IntegerQ[ByteCount[auditDeep]]");
    audit_true("ByteCount#1",
               "ByteCount[auditDeep] - ByteCount[Nest[f, x, 199999]] == "
               "ByteCount[f[x]] - ByteCount[x]");
    audit_chk("ByteCount#1", "LeafCount[auditDeep]", "200001");
    audit_run("Clear[auditDeep]");
    /* Shallow values: a list's count is its own node plus its parts'. */
    audit_true("ByteCount#1",
               "ByteCount[{a, \"bc\"}] > ByteCount[a] + ByteCount[\"bc\"] + ByteCount[List]");
    audit_true("ByteCount#1", "ByteCount[2^1000] > ByteCount[2]");
}

int main(void) {
    audit_init();
    TEST(test_bytecount_deep_nesting);
    return audit_finish("audit_expression_information");
}
