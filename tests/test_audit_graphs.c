/* Builtin-audit regressions: graphs (review/issues/<Head>.md; tests/test_audit.h). */
#include "test_audit.h"

/* CompleteGraph#1 (critical): CompleteGraph[10^7] crashed (unchecked calloc
 * of n(n-1)/2 edge slots).  The documented resource limit (10^8 vertices or
 * 5x10^7 edges) leaves such a family unevaluated. */
static void test_complete_graph_resource_limit(void) {
    audit_chk("CompleteGraph#1", "CompleteGraph[10^7]", "CompleteGraph[10000000]");
    audit_chk("CompleteGraph#1", "CompleteGraph[2^62]", "CompleteGraph[4611686018427387904]");
    audit_chk("CompleteGraph#1", "EdgeCount[CompleteGraph[5]]", "10");
    audit_chk("CompleteGraph#1", "EdgeCount[CompleteGraph[1]]", "0");
    audit_chk("CompleteGraph#1", "VertexCount[CompleteGraph[0]]", "0");
}

/* StarGraph#1 (high): StarGraph[10^15] crashed (unchecked calloc). */
static void test_star_graph_resource_limit(void) {
    audit_chk("StarGraph#1", "StarGraph[10^15]", "StarGraph[1000000000000000]");
    audit_chk("StarGraph#1", "EdgeCount[StarGraph[6]]", "5");
    audit_chk("StarGraph#1", "CycleGraph[10^15]", "CycleGraph[1000000000000000]");
    audit_chk("StarGraph#1", "PathGraph[10^15]", "PathGraph[1000000000000000]");
    audit_chk("StarGraph#1", "RandomGraph[{2^40, 1}]", "RandomGraph[{1099511627776, 1}]");
    audit_chk("StarGraph#1", "EdgeCount[CycleGraph[4]]", "4");
}

int main(void) {
    audit_init();
    TEST(test_complete_graph_resource_limit);
    TEST(test_star_graph_resource_limit);
    return audit_finish("audit_graphs");
}
