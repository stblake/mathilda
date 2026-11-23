/* Builtin-audit regressions: string operations (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* StringRiffle#1 (medium): "StringRiffle segfaults on a deeply nested list
 * (unbounded recursion)".  Both walks are iterative now. */
static void test_stringriffle_deep_nesting(void) {
    audit_run("auditDeep = Nest[List, \"a\", 100000]");
    audit_chk("StringRiffle#1", "StringLength[StringRiffle[auditDeep]]", "1");
    audit_true("StringRiffle#1", "StringRiffle[auditDeep] === \"a\"");
    audit_run("Clear[auditDeep]");
    /* Deep and wide: one separator per level, the default scheme. */
    audit_true("StringRiffle#1",
               "StringRiffle[Nest[{#, #} &, \"a\", 3]] === "
               "\"a a\\na a\\n\\na a\\na a\"");
    /* Every documented form is unchanged. */
    audit_true("StringRiffle#1", "StringRiffle[{\"a\", \"b\", \"c\"}] === \"a b c\"");
    audit_true("StringRiffle#1", "StringRiffle[{{1, 2}, {3, 4}}] === \"1 2\\n3 4\"");
    audit_true("StringRiffle#1",
               "StringRiffle[{\"a\", \"b\"}, {\"(\", \", \", \")\"}] === \"(a, b)\"");
    audit_true("StringRiffle#1",
               "StringRiffle[{{a, b}, {c}}, {\"[\", \";\", \"]\"}, \"-\"] === \"[a-b;c]\"");
    audit_true("StringRiffle#1", "StringRiffle[{}] === \"\"");
    audit_true("StringRiffle#1", "StringRiffle[{{}, {x}}] === \"\\nx\"");
}

int main(void) {
    audit_init();
    TEST(test_stringriffle_deep_nesting);
    return audit_finish("audit_string_operations");
}
