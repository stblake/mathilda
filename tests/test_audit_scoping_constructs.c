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

/* With#1 (high): "With captures free variables in nested Function/Module/
 * With (no renaming)". */
static void test_with_capture(void) {
    audit_run("ClearAll[x, y]");
    audit_true("With#1", "With[{y = x}, Function[x, x + y][10]] === 10 + x");
    audit_true("With#1", "With[{y = x}, Module[{x = 2}, x + y]] === 2 + x");
    audit_true("With#1", "With[{y = x}, With[{x = 2}, x + y]] === 2 + x");
    audit_true("With#1", "(With[{y = x}, Function[x, x + y]] /. x -> 1)[10] === 11");
    /* sanity: no clash, no rename; init values still see the outer value */
    audit_chk("With#1", "With[{y = 3}, Function[x, x + y][10]]", "13");
    audit_chk("With#1", "With[{q = 12}, With[{k = q}, k]]", "12");
    audit_chk("With#1", "With[{y = 1}, Module[{z = y}, z + 1]]", "2");
}

/* Wave 4 G4: malformed local specs (Block/Module::lvsym, With::lvset), long
 * Module names, nested With lists, Unique's name check and attributes. */
static void test_wave4_g4_scoping(void) {
    audit_chk("Block#2", "Check[Block[{1}, 2], fired, Block::lvsym]", "fired");
    audit_chk("Block#2", "Quiet[Block[{1}, 2]]", "Block[{1}, 2]");
    audit_chk("Block#2", "Check[Module[{f[1]}, 2], fired, Module::lvsym]", "fired");
    audit_chk("Block#2", "Block[{bx = 1, by}, {bx, by}]", "{1, by}");
    audit_run("r1g4 = Module[{vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv = 1}, Hold[vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv]]; r2g4 = Module[{vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv = 2}, Hold[vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv]]");
    audit_chk("Module#1", "{ReleaseHold[r1g4], ReleaseHold[r2g4]}", "{1, 2}");
    audit_chk("With#3", "Check[With[{x}, x], fired, With::lvset]", "fired");
    audit_chk("With#3", "Quiet[{With[{x}, x], With[{1 = 2}, 3], With[{f[1] = 2}, 3]}]", "{With[{x}, x], With[{1 = 2}, 3], With[{f[1] = 2}, 3]}");
    audit_chk("With#4", "With[{x = 1}, {y = x + 1}, y]", "2");
    audit_chk("With#4", "With[{x = 1}, {y = x + 1}, {z = y + x}, {x, y, z}]", "{1, 2, 3}");
    audit_chk("Unique#1", "Check[Unique[\"1abc\"], fired, Unique::usym]", "fired");
    audit_chk("Unique#1", "Quiet[{Unique[\"a b\"], Unique[\"\"]}]", "{Unique[\"a b\"], Unique[\"\"]}");
    audit_chk("Unique#2", "Attributes[Evaluate[Unique[yg4, {Flat, Orderless}]]]", "{Flat, Orderless, Temporary}");
    audit_chk("Unique#2", "Attributes[Evaluate[Unique[yg4, Temporary]]]", "{Temporary}");
}

int main(void) {
    audit_init();
    TEST(test_block_duplicate_local);
    TEST(test_with_pattern_names);
    TEST(test_with_capture);
    TEST(test_wave4_g4_scoping);
    return audit_finish("audit_scoping_constructs");
}
