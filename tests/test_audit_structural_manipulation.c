/* Builtin-audit regressions: structural manipulation (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* ArrayFlatten#1 (critical): "ArrayFlatten segfaults when a block is nested
 * more than 64 levels deep (stack buffer overflow in af_tensor_dims)".  A
 * block needs only its first r dimensions; deeper levels are content. */
static void test_array_flatten_deep_block(void) {
    audit_true("ArrayFlatten#1",
               "ArrayFlatten[{{Nest[List, 1, 100]}}] === Nest[List, 1, 100]");
    audit_true("ArrayFlatten#1",
               "ArrayFlatten[{{Nest[List, 1, 5000]}}] === Nest[List, 1, 5000]");
    audit_true("ArrayFlatten#1",
               "ArrayFlatten[{{{{Nest[List, 1, 200]}}}}, 1] === {{{Nest[List, 1, 200]}}}");
    /* Content jagged below the r flattened levels is carried along. */
    audit_chk("ArrayFlatten#1", "ArrayFlatten[{{{{1, {2}}, {3, 4}}}}]", "{{1, {2}}, {3, 4}}");
    /* Ordinary block matrices and scalar padding are unchanged. */
    audit_chk("ArrayFlatten#1", "ArrayFlatten[{{{{1, 2}, {3, 4}}, 0}, {0, {{5}}}}]",
              "{{1, 2, 0}, {3, 4, 0}, {0, 0, 5}}");
    audit_chk("ArrayFlatten#1", "ArrayFlatten[{{{{1, 2}}, {{3}, {4}}}}]",
              "ArrayFlatten[{{{{1, 2}}, {{3}, {4}}}}]");
}

/* Depth#1 (high): "Depth segfaults on a malformed Rule option (Rule[] or
 * Rule[Heads])".  Now Depth::nonopt and the call stays unevaluated. */
static void test_depth_malformed_option(void) {
    audit_chk("Depth#1", "Check[Depth[x, Rule[]], fired, Depth::nonopt]", "fired");
    audit_chk("Depth#1", "Quiet[Depth[x, Rule[]]]", "Depth[x, Rule[]]");
    audit_chk("Depth#1", "Check[Depth[x, Rule[Heads]], fired, Depth::nonopt]", "fired");
    audit_chk("Depth#1", "Quiet[Depth[x, y]]", "Depth[x, y]");
    audit_chk("Depth#1", "Depth[h[{{{a}}}][x, y], Heads -> True]", "6");
    audit_chk("Depth#1", "Depth[h[{{{a}}}][x, y], Heads -> False]", "2");
    audit_chk("Depth#1", "Depth[a -> b]", "2");
}

/* Depth#2 (medium): "Depth crashes (stack overflow) on deeply nested
 * expressions".  The walk is iterative; LeafCount crashed the same way. */
static void test_depth_deep_nesting(void) {
    audit_run("auditDeep = Nest[List, 1, 200000]");
    audit_chk("Depth#2", "Depth[auditDeep]", "200001");
    audit_chk("Depth#2", "LeafCount[auditDeep]", "200001");
    audit_chk("Depth#2", "LeafCount[auditDeep, Heads -> False]", "1");
    audit_run("auditDeepH = Nest[g[#][x] &, y, 100000]");
    audit_chk("Depth#2", "Depth[auditDeepH]", "2");
    audit_chk("Depth#2", "Depth[auditDeepH, Heads -> True]", "200001");
    audit_run("Clear[auditDeep, auditDeepH]");
    audit_chk("Depth#2", "Depth[{{1, 2}, {3, {4}}}]", "4");
    audit_chk("Depth#2", "Depth[<|a -> <|b -> 1|>|>]", "3");
    audit_chk("Depth#2", "Depth[{}]", "2");
    audit_chk("Depth#2", "Depth[{1/2, 1 + I}]", "2");
    audit_chk("Depth#2", "LeafCount[1 + a + b^2]", "6");
}

int main(void) {
    audit_init();
    TEST(test_array_flatten_deep_block);
    TEST(test_depth_malformed_option);
    TEST(test_depth_deep_nesting);
    return audit_finish("audit_structural_manipulation");
}
