/* Builtin-audit regressions: scoping constructs (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* Block#1 (critical): "Duplicate local variable in Block overwrites the
 * global value on exit".  x = 5; Block[{x = 1, x = 2}, x] left x = 1: the
 * frame saved x twice and restored the stale second slot last.  A duplicate
 * local is rejected with ::dup and the construct stays unevaluated, for
 * Module and With as for Block; the global is untouched. */
static void test_block_duplicate_local(void) {
    audit_run("ClearAll[bdx]; bdx = 5");
    audit_chk("Block#1", "Check[Block[{bdx = 1, bdx = 2}, bdx], fired, Block::dup]", "fired");
    audit_chk("Block#1", "bdx", "5");
    audit_chk("Block#1", "Quiet[Block[{bdx = 1, bdx = 2}, bdx]]", "Block[{bdx = 1, bdx = 2}, bdx]");
    audit_chk("Block#1", "bdx", "5");
    audit_chk("Block#1", "Check[Block[{bdx, bdx = 2}, bdx], fired, Block::dup]", "fired");
    audit_chk("Block#1", "bdx", "5");
    audit_chk("Block#1", "Check[Module[{bdy, bdy}, 1], fired, Module::dup]", "fired");
    audit_chk("Block#1", "Quiet[Module[{bdy, bdy = 2}, 1]]", "Module[{bdy, bdy = 2}, 1]");
    audit_chk("Block#1", "Check[With[{bdz = 1, bdz = 2}, bdz], fired, With::dup]", "fired");
    audit_chk("Block#1", "Check[With[{bdz = 1, bdz := 2}, bdz], fired, With::dup]", "fired");
    /* sanity: distinct locals still work and restore */
    audit_chk("Block#1", "Block[{bdx = 1, bdw = 2}, bdx + bdw]", "3");
    audit_chk("Block#1", "bdx", "5");
    audit_chk("Block#1", "Module[{p = 1, q = 2}, p + q]", "3");
    audit_chk("Block#1", "With[{p = 1, q = 2}, p + q]", "3");
}

/* With#2 (medium): "With substitutes into pattern names, producing invalid
 * patterns like Pattern[5, _]".  A rule or definition scopes its own pattern
 * variables: one named like the With variable is renamed to x$, not given
 * the value. */
static void test_with_pattern_names(void) {
    audit_run("ClearAll[wpx, wpf, wpg]");
    audit_chk("With#2", "With[{wpx = 5}, wpx_ -> wpx]", "wpx$_ -> wpx$");
    audit_chk("With#2", "With[{wpx = 2}, wpx_Integer :> wpx]", "wpx$_Integer :> wpx$");
    audit_run("With[{wpx = 5}, wpf[wpx_] := wpx]");
    audit_chk("With#2", "wpf[2]", "2");
    audit_chk("With#2", "3 /. With[{wpx = 5}, wpx_ -> wpx]", "3");
    audit_chk("With#2", "Module[{wpx}, wpx_ -> wpx]", "wpx$_ -> wpx$");
    /* sanity: a variable used, not bound, by the rule is still substituted */
    audit_chk("With#2", "With[{wpx = 5}, y_ -> y + wpx]", "y_ -> 5 + y");
    audit_run("With[{wpx = 5}, wpg[y_] := y + wpx]");
    audit_chk("With#2", "wpg[1]", "6");
}

int main(void) {
    audit_init();
    TEST(test_block_duplicate_local);
    TEST(test_with_pattern_names);
    return audit_finish("audit_scoping_constructs");
}
