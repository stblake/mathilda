/* NDArray support for the selection and structural heads.
 *
 * An NDArray is an ATOMIC value, so a builtin that walks
 * `arg->data.function.args` looks straight past one and returns the call
 * UNEVALUATED — which is exactly what Select, TakeWhile, SelectFirst,
 * AllTrue/AnyTrue/NoneTrue, Join, Differences, First/Last/Most/Rest,
 * RotateLeft/RotateRight, Riffle and Partition all did, while the identical
 * List call worked fine.
 *
 * The contract these tests enforce is the one the rest of the ND layer keeps
 * (src/ndstruct.h): the packed answer must be identical to the equivalent List
 * answer. Every case is therefore written as a comparison against the List
 * call itself rather than against a hand-written expected value — a literal
 * would only re-state the List implementation, and would not notice if BOTH
 * paths changed together.
 *
 * `Normal[...] === <List call>` is the assertion: SameQ, so element values,
 * length and structure must all match, and Normal unpacks whichever side came
 * back packed so a kind difference does not mask a value difference. Kinds are
 * asserted separately below. */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void run(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* The packed call and the List call must agree, for every head. */
static void agrees(const char* packed, const char* listed) {
    char buf[512];
    snprintf(buf, sizeof buf, "Normal[%s] === (%s)", packed, listed);
    Expr* e = parse_expression(buf);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, "True") == 0,
               "packed/List disagree: %s  vs  %s  (got %s)", packed, listed, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void test_select_family(void) {
    agrees("Select[NDArray[{1., 2., 3.}], # > 1 &]", "Select[{1., 2., 3.}, # > 1 &]");
    agrees("Select[NDArray[{1., 2., 3.}], # > 0 &]", "Select[{1., 2., 3.}, # > 0 &]");
    /* No survivors: an empty result has no packed form, so it stays a List —
     * which is what the List call gives too. */
    agrees("Select[NDArray[{1., 2.}], # > 9 &]", "Select[{1., 2.}, # > 9 &]");
    run("Select[NDArray[{1., 2.}], # > 9 &]", "{}");
    /* The n-limited form. */
    agrees("Select[NDArray[{1., 2., 3., 4.}], # > 1 &, 2]",
           "Select[{1., 2., 3., 4.}, # > 1 &, 2]");
}

static void test_takewhile_lengthwhile_selectfirst(void) {
    agrees("TakeWhile[NDArray[{1., 2., 9., 3.}], # < 5 &]",
           "TakeWhile[{1., 2., 9., 3.}, # < 5 &]");
    agrees("TakeWhile[NDArray[{9., 1.}], # < 5 &]", "TakeWhile[{9., 1.}, # < 5 &]");
    agrees("LengthWhile[NDArray[{1., 2., 9.}], # < 5 &]",
           "LengthWhile[{1., 2., 9.}, # < 5 &]");
    agrees("SelectFirst[NDArray[{1., 7., 3.}], # > 5 &]",
           "SelectFirst[{1., 7., 3.}, # > 5 &]");
    /* No match: Missing["NotFound"] from both sides. */
    agrees("SelectFirst[NDArray[{1., 2.}], # > 9 &]", "SelectFirst[{1., 2.}, # > 9 &]");
    agrees("SelectFirst[NDArray[{1., 2.}], # > 9 &, 0.]",
           "SelectFirst[{1., 2.}, # > 9 &, 0.]");
}

static void test_all_any_none(void) {
    agrees("AllTrue[NDArray[{1., 2.}], # > 0 &]",  "AllTrue[{1., 2.}, # > 0 &]");
    agrees("AllTrue[NDArray[{1., 2.}], # > 1 &]",  "AllTrue[{1., 2.}, # > 1 &]");
    agrees("AnyTrue[NDArray[{1., 2.}], # > 1 &]",  "AnyTrue[{1., 2.}, # > 1 &]");
    agrees("AnyTrue[NDArray[{1., 2.}], # > 9 &]",  "AnyTrue[{1., 2.}, # > 9 &]");
    agrees("NoneTrue[NDArray[{1., 2.}], # > 9 &]", "NoneTrue[{1., 2.}, # > 9 &]");
    agrees("NoneTrue[NDArray[{1., 2.}], # > 1 &]", "NoneTrue[{1., 2.}, # > 1 &]");
}

static void test_join_and_differences(void) {
    agrees("Join[NDArray[{1., 2.}], NDArray[{3.}]]", "Join[{1., 2.}, {3.}]");
    agrees("Join[NDArray[{1., 2.}], {3.}]",          "Join[{1., 2.}, {3.}]");
    agrees("Join[{1., 2.}, NDArray[{3.}]]",          "Join[{1., 2.}, {3.}]");
    agrees("Join[NDArray[{1.}], NDArray[{2.}], NDArray[{3.}]]",
           "Join[{1.}, {2.}, {3.}]");
    agrees("Differences[NDArray[{1., 4., 9.}]]", "Differences[{1., 4., 9.}]");
    agrees("Differences[NDArray[{1., 4., 9., 16.}], 2]",
           "Differences[{1., 4., 9., 16.}, 2]");
}

static void test_first_last_most_rest(void) {
    agrees("First[NDArray[{1., 2.}]]", "First[{1., 2.}]");
    agrees("Last[NDArray[{1., 2.}]]",  "Last[{1., 2.}]");
    agrees("Most[NDArray[{1., 2., 3.}]]", "Most[{1., 2., 3.}]");
    agrees("Rest[NDArray[{1., 2., 3.}]]", "Rest[{1., 2., 3.}]");
    /* At rank 2 First/Last give a ROW, which packs back to a rank-1 array. */
    agrees("First[NDArray[{{1., 2.}, {3., 4.}}]]", "First[{{1., 2.}, {3., 4.}}]");
    agrees("Rest[NDArray[{{1., 2.}, {3., 4.}}]]",  "Rest[{{1., 2.}, {3., 4.}}]");
    /* The default form still answers when there is nothing to take. */
    agrees("First[NDArray[{1.}], 0.]", "First[{1.}, 0.]");
}

static void test_rotate_riffle_partition(void) {
    agrees("RotateLeft[NDArray[{1., 2., 3.}], 1]",  "RotateLeft[{1., 2., 3.}, 1]");
    agrees("RotateLeft[NDArray[{1., 2., 3.}]]",     "RotateLeft[{1., 2., 3.}]");
    agrees("RotateRight[NDArray[{1., 2., 3.}], 1]", "RotateRight[{1., 2., 3.}, 1]");
    agrees("RotateRight[NDArray[{1., 2., 3.}]]",    "RotateRight[{1., 2., 3.}]");
    agrees("Riffle[NDArray[{1., 2., 3.}], 0.]",     "Riffle[{1., 2., 3.}, 0.]");
    agrees("Partition[NDArray[{1., 2., 3., 4.}], 2]", "Partition[{1., 2., 3., 4.}, 2]");
    agrees("Partition[NDArray[{1., 2., 3., 4., 5.}], 2, 1]",
           "Partition[{1., 2., 3., 4., 5.}, 2, 1]");
}

static void test_mapthread_repacks(void) {
    agrees("MapThread[Plus, {NDArray[{1., 2.}], NDArray[{3., 4.}]}]",
           "MapThread[Plus, {{1., 2.}, {3., 4.}}]");
    agrees("MapThread[Plus, {NDArray[{1., 2.}], {3., 4.}}]",
           "MapThread[Plus, {{1., 2.}, {3., 4.}}]");
}

/* Packed in, packed out — the convention Map and the ndstruct_* ops keep. A
 * head that unpacked its result would force a re-pack on every pipeline stage. */
static void test_result_stays_packed(void) {
    run("Head[Select[NDArray[{1., 2., 3.}], # > 1 &]]", "NDArray");
    run("Head[Join[NDArray[{1., 2.}], NDArray[{3.}]]]", "NDArray");
    run("Head[Rest[NDArray[{1., 2., 3.}]]]",            "NDArray");
    run("Head[RotateLeft[NDArray[{1., 2.}], 1]]",       "NDArray");
    run("Head[Partition[NDArray[{1., 2., 3., 4.}], 2]]", "NDArray");
    /* ...but a scalar result is a scalar, not a length-1 array. */
    run("Head[First[NDArray[{1., 2.}]]]",               "Real");
    run("Head[LengthWhile[NDArray[{1., 2.}], # < 5 &]]", "Integer");
    run("Head[AllTrue[NDArray[{1., 2.}], # > 0 &]]",    "Symbol");
    /* A result that cannot be packed stays exactly as the List path built it. */
    run("Head[Riffle[NDArray[{1., 2.}], zz]]",          "List");
}

/* The argument is borrowed: none of these may mutate what the caller passed. */
static void test_argument_not_mutated(void) {
    run("nds = NDArray[{1., 2., 3.}]; Select[nds, # > 1 &]; Normal[nds]",
        "{1.0, 2.0, 3.0}");
    run("ndj = NDArray[{1., 2.}]; Join[ndj, ndj]; Normal[ndj]", "{1.0, 2.0}");
    run("ndr = NDArray[{1., 2., 3.}]; RotateLeft[ndr, 1]; Normal[ndr]",
        "{1.0, 2.0, 3.0}");
}

/* Plain Lists must be untouched by any of this. */
static void test_plain_lists_unregressed(void) {
    run("Select[{1, 2, 3}, # > 1 &]", "{2, 3}");
    run("Join[{1, 2}, {3}]",          "{1, 2, 3}");
    run("First[{a, b}]",              "a");
    run("Rest[{a, b, c}]",            "{b, c}");
    run("Riffle[{a, b}, x]",          "{a, x, b}");
    run("Partition[{1, 2, 3, 4}, 2]", "{{1, 2}, {3, 4}}");
    run("RotateLeft[{1, 2, 3}]",      "{2, 3, 1}");
    /* An atom is still not a collection. */
    run("First[5]",                   "First[5]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_select_family);
    TEST(test_takewhile_lengthwhile_selectfirst);
    TEST(test_all_any_none);
    TEST(test_join_and_differences);
    TEST(test_first_last_most_rest);
    TEST(test_rotate_riffle_partition);
    TEST(test_mapthread_repacks);
    TEST(test_result_stays_packed);
    TEST(test_argument_not_mutated);
    TEST(test_plain_lists_unregressed);

    printf("All NDArray selection/structural tests passed!\n");
    return 0;
}
