/* Unit tests for Sow / Reap (Wolfram-Language dynamic-scope result collection).
 *
 * Covers: basic reaping from a CompoundExpression, tag grouping and ordering,
 * tag lists (with repeats), the {p1,...} pattern-list form, the f[tag,{e...}]
 * handler form, tag-selective reaping, empty reaps, nested Reap routing, Sow's
 * return value, and interaction with a Throw crossing a Reap. See
 * src/sowreap.c and the registration in src/core.c.
 *
 * NOTE on Sum: Mathilda's Sum evaluates its summand symbolically once and
 * closed-forms the result (a pre-existing Sum property that affects every side
 * effect — Print, Sow alike), so Sow inside Sum fires only once. The faithful
 * iteration primitive is Do, which we use here for the "reap across a loop"
 * case (WL: Reap[Sum[Sow[i^2]+1,{i,10}]] == {395,{{1,4,...,100}}}).
 */
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ----- Basic reaping from a sequence of expressions ----- */
void test_basic() {
    assert_eval_eq("Reap[Sow[a]; b; Sow[c]; Sow[d]; e]", "{e, {{a, c, d}}}", 0);
    assert_eval_eq("Reap[x]", "{x, {}}", 0);            /* nothing sown */
    assert_eval_eq("Reap[5]", "{5, {}}", 0);
    /* Reap collects in exactly the order values are sown (MapAll traversal). */
    assert_eval_eq("Reap[MapAll[Sow, (a + b) (c + x^2)];]",
                   "{Null, {{a, b, a + b, c, x, 2, x^2, c + x^2, (a + b) (c + x^2)}}}", 0);
}

/* ----- Reap across an iteration produces per-step numeric values ----- */
void test_iteration() {
    assert_eval_eq("Reap[Do[Sow[i^2], {i, 10}]]",
                   "{Null, {{1, 4, 9, 16, 25, 36, 49, 64, 81, 100}}}", 0);
    /* Sum computes the value; the loop equivalent above proves the reap. */
    assert_eval_eq("Reap[Module[{s = 0}, Do[s = s + Sow[i^2], {i, 3}]; s]]",
                   "{14, {{1, 4, 9}}}", 0);
    /* Table returns its list as the value, with the sown values alongside. */
    assert_eval_eq("Reap[Table[Sow[i^2], {i, 3}]]", "{{1, 4, 9}, {{1, 4, 9}}}", 0);
}

/* ----- Tag grouping: distinct tags -> distinct sublists, first-encounter order ----- */
void test_tag_grouping() {
    assert_eval_eq("Reap[Sow[1, x]; Sow[2, y]; Sow[3, x]; Sow[4, y]]",
                   "{4, {{1, 3}, {2, 4}}}", 0);
    /* The list for the first tag encountered is given first. */
    assert_eval_eq("Reap[Sow[1, y]; Sow[2, x]; Sow[3, y]]", "{3, {{1, 3}, {2}}}", 0);
}

/* ----- Tag lists on Sow (repeats make a value appear multiple times) ----- */
void test_tag_lists() {
    assert_eval_eq("Reap[Sow[1, {a, a, a}], _, Rule]", "{1, {a -> {1, 1, 1}}}", 0);
    /* Sow[e, {{tag}}] sows once under the single tag {tag}. */
    assert_eval_eq("Reap[Sow[9, {{q}}], _, Rule]", "{9, {{q} -> {9}}}", 0);
}

/* ----- Reap[expr, {p1,...}] makes one slot per pattern ----- */
void test_pattern_list() {
    assert_eval_eq("Reap[Sow[1, {x, x}]; Sow[2, y]; Sow[3, x], {x, x, y}]",
                   "{3, {{{1, 1, 3}}, {{1, 1, 3}}, {{2}}}}", 0);
    /* An empty list argument reaps nothing into zero slots. */
    assert_eval_eq("Reap[Sow[1], {}]", "{1, {}}", 0);
    /* A pattern that matches no tag gets an empty slot. */
    assert_eval_eq("Reap[Sow[1, a]; Sow[2, b], {a, c}]", "{2, {{{1}}, {}}}", 0);
}

/* ----- Tag-selective reaping ----- */
void test_selective() {
    assert_eval_eq("Reap[Sow[1, x]; Sow[2, y]; Sow[3, x]; Sow[4, y], x]",
                   "{4, {{1, 3}}}", 0);
}

/* ----- The f[tag, {e...}] handler form ----- */
void test_handler_form() {
    assert_eval_eq("Reap[Sow[1, {x, x}]; Sow[2, y]; Sow[3, x], _, f]",
                   "{3, {f[x, {1, 1, 3}], f[y, {2}]}}", 0);
    assert_eval_eq("Reap[Sow[1, {x, x}]; Sow[2, y]; Sow[3, x], _, Rule]",
                   "{3, {x -> {1, 1, 3}, y -> {2}}}", 0);
    /* Handler with a pure function (tag -> total of values). */
    assert_eval_eq("Reap[Sow[3, a]; Sow[4, a]; Sow[5, b], _, #1 -> Total[#2] &]",
                   "{5, {a -> 7, b -> 5}}", 0);
}

/* ----- Sow returns its first argument (even outside any Reap) ----- */
void test_sow_returns() {
    assert_eval_eq("Sow[7]", "7", 0);
    assert_eval_eq("Sow[7, mytag]", "7", 0);
    assert_eval_eq("Reap[Sow[42]]", "{42, {{42}}}", 0);
    /* Nested Sow: inner sows+returns, then outer sows the returned value. */
    assert_eval_eq("Reap[Sow[Sow[1]]]", "{1, {{1, 1}}}", 0);
}

/* ----- Nested Reap: a non-matching inner tag propagates to the outer Reap ----- */
void test_nested() {
    assert_eval_eq("Reap[Sow[1, x]; Reap[Sow[2, y]; Sow[3, x], y]]",
                   "{{3, {{2}}}, {{1, 3}}}", 0);
    /* Inner Reap consumes Sow[2]; Sow[1]/Sow[3] reach the outer Reap. The inner
     * Reap is a discarded middle statement, so the value is Sow[3] -> 3. */
    assert_eval_eq("Reap[Sow[1]; Reap[Sow[2]]; Sow[3]]",
                   "{3, {{1, 3}}}", 0);
}

/* ----- A Throw crossing a Reap unwinds cleanly (sown values discarded) ----- */
void test_throw_across_reap() {
    assert_eval_eq("Catch[Reap[Sow[1]; Throw[q]; Sow[2]]]", "q", 0);
    /* Reap completes normally if the Throw is caught inside it. */
    assert_eval_eq("Reap[Catch[Sow[1]; Throw[q]; Sow[2]]]", "{q, {{1}}}", 0);
}

/* ----- Attributes and arity errors ----- */
void test_attributes_and_errors() {
    assert_eval_eq("Attributes[Sow]", "{Protected}", 0);
    assert_eval_eq("Attributes[Reap]", "{HoldFirst, Protected}", 0);
    /* Wrong arity leaves the expression unevaluated (emits ::argt on stderr). */
    assert_eval_eq("Sow[]", "Sow[]", 0);
    assert_eval_eq("Reap[]", "Reap[]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_basic);
    TEST(test_iteration);
    TEST(test_tag_grouping);
    TEST(test_tag_lists);
    TEST(test_pattern_list);
    TEST(test_selective);
    TEST(test_handler_form);
    TEST(test_sow_returns);
    TEST(test_nested);
    TEST(test_throw_across_reap);
    TEST(test_attributes_and_errors);

    printf("All sow/reap tests passed!\n");
    symtab_clear();
    return 0;
}
