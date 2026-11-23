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

/* An impossible size aborts with General::nomem: CheckAbort sees the abort and
 * Check sees the (routed) message. */
static void nomem_chk(const char* id, const char* expr) {
    char buf[512];
    snprintf(buf, sizeof buf, "Quiet[CheckAbort[%s, aborted]]", expr);
    audit_chk(id, buf, "aborted");
    snprintf(buf, sizeof buf, "Quiet[Check[CheckAbort[%s, 0], fired, General::nomem]]", expr);
    audit_chk(id, buf, "fired");
}

/* PadLeft#1 (high): a huge machine-integer length crashed the kernel
 * (unchecked malloc; -n overflowed at INT64_MIN).  PadRight shares the
 * driver, the packed rank-1 path and ArrayPad have their own buffers. */
static void test_pad_huge_length(void) {
    nomem_chk("PadLeft#1", "PadLeft[{1}, 10^15]");
    nomem_chk("PadLeft#1", "PadLeft[{1}, 2^62]");
    nomem_chk("PadLeft#1", "PadLeft[{1}, -2^63]");
    nomem_chk("PadLeft#1", "PadRight[{1}, 10^15]");
    nomem_chk("PadLeft#1", "PadRight[{{1}}, {10^10, 10^10}]");
    nomem_chk("PadLeft#1", "PadLeft[Range[100], 2^62]");          /* packed */
    nomem_chk("PadLeft#1", "PadRight[N[Range[100]], 2^62, 0.]");  /* packed */
    nomem_chk("PadLeft#1", "ArrayPad[{1}, 10^15]");
    nomem_chk("PadLeft#1", "ArrayPad[{1}, {0, 2^62}]");
    nomem_chk("PadLeft#1", "ArrayPad[Range[100], {0, 2^62}]");     /* packed */
    audit_chk("PadLeft#1", "PadLeft[{1, 2}, 4]", "{0, 0, 1, 2}");
    audit_chk("PadLeft#1", "PadRight[{1, 2}, -4]", "{0, 0, 1, 2}");
    audit_chk("PadLeft#1", "PadRight[{{1}, {2, 3}}, {2, 3}, x]", "{{1, x, x}, {2, 3, x}}");
    audit_chk("PadLeft#1", "ArrayPad[{1, 2}, 1]", "{0, 1, 2, 0}");
    audit_chk("PadLeft#1", "ArrayPad[Range[5], {-2, 1}]", "{3, 4, 5, 0}");
    audit_chk("PadLeft#1", "ArrayPad[Range[3], -2^63]", "{}");
}

/* Drop#1 (high): Drop[{a,b,c},{1,2^63-1}] segfaulted -- the index count came
 * from the raw spec before the range check (also Take, and the packed path). */
static void test_drop_huge_sequence_spec(void) {
    audit_chk("Drop#1", "Drop[{a, b, c}, {1, 9223372036854775807}]",
              "Drop[{a, b, c}, {1, 9223372036854775807}]");
    audit_chk("Drop#1", "Drop[{a, b, c}, {1, 2^62}]", "Drop[{a, b, c}, {1, 4611686018427387904}]");
    audit_chk("Drop#1", "Take[{a, b, c}, {1, 9223372036854775807}]",
              "Take[{a, b, c}, {1, 9223372036854775807}]");
    audit_chk("Drop#1", "Take[{a, b, c}, {3, -9223372036854775807 - 1, -1}]",
              "Take[{a, b, c}, {3, -9223372036854775808, -1}]");
    audit_chk("Drop#1", "Take[{a, b, c}, {1, -1, -9223372036854775807 - 1}]", "{}");
    audit_chk("Drop#1", "Head[Take[Range[10], {1, 2^63 - 1}]]", "Take");
    audit_chk("Drop#1", "Drop[{a, b, c, d, e}, {1, -1, 2}]", "{b, d}");
    audit_chk("Drop#1", "Take[{a, b, c, d, e}, {5, 1, -2}]", "{e, c, a}");
    audit_chk("Drop#1", "Take[{a, b, c, d, e}, {1, 4, 2}]", "{a, c}");
    audit_chk("Drop#1", "Take[Range[10], {2, 9, 3}]", "{2, 5, 8}");
}

int main(void) {
    audit_init();
    TEST(test_array_flatten_deep_block);
    TEST(test_depth_malformed_option);
    TEST(test_depth_deep_nesting);
    TEST(test_pad_huge_length);
    TEST(test_drop_huge_sequence_spec);
    return audit_finish("audit_structural_manipulation");
}
