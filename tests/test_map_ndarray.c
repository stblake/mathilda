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

/* NDArray support for the Map family (Map, MapIndexed, MapThread, MapAll).
 * An NDArray is an atomic value, so before these changes each of these was a
 * silent no-op over one. Assertions avoid bare machine reals (whose printed
 * form is formatting-sensitive) by mapping Head over elements, checking a Part
 * via equality, or counting _Real leaves. */
static void run(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Map ---------- */

static void test_map_ndarray_rank1(void) {
    /* f is applied to each scalar element (boxed as Real). */
    run("Map[Head, NDArray[Range[4]]]", "{Real, Real, Real, Real}");
    /* An explicit {1} level spec takes the same fast path. */
    run("Map[Head, NDArray[Range[4]], {1}]", "{Real, Real, Real, Real}");
}

static void test_map_ndarray_rank2_rows(void) {
    /* Level-1 parts of a rank-2 array are sub-NDArray rows. */
    run("Map[Head, NDArray[Table[i+j,{i,2},{j,3}]]]", "{NDArray, NDArray}");
}

static void test_map_ndarray_numeric_repacks(void) {
    /* A numeric result repacks into an NDArray (packed in -> packed out). */
    run("Head[Map[#^2 &, NDArray[Range[4]]]]", "NDArray");
    run("Map[#^2 &, NDArray[Range[4]]][[3]] == 9", "True");
}

static void test_map_ndarray_symbolic_is_list(void) {
    /* A symbolic result stays a List. */
    run("Head[Map[f, NDArray[Range[4]]]]", "List");
}

static void test_map_ndarray_deep_level(void) {
    /* A non-default level spec materializes and maps generically. */
    run("Map[Head, NDArray[Table[i+j,{i,2},{j,3}]], {2}]",
        "{{Real, Real, Real}, {Real, Real, Real}}");
}

/* ---------- MapIndexed ---------- */

static void test_mapindexed_ndarray_positions(void) {
    run("MapIndexed[#2 &, NDArray[Range[3]]]", "{{1}, {2}, {3}}");
    run("MapIndexed[#2 &, NDArray[Table[i+j,{i,2},{j,2}]]]", "{{1}, {2}}");
}

static void test_mapindexed_ndarray_values(void) {
    /* The value passed to f is the leading-axis part: scalar for rank 1, a
     * sub-NDArray row for rank 2. */
    run("MapIndexed[Head[#1] &, NDArray[Range[3]]]", "{Real, Real, Real}");
    run("MapIndexed[Head[#1] &, NDArray[Table[i+j,{i,2},{j,2}]]]", "{NDArray, NDArray}");
}

/* ---------- MapThread ---------- */

static void test_mapthread_ndarray_with_list(void) {
    /* An NDArray entry threads element-wise against a plain List. */
    run("MapThread[#2 &, {NDArray[Range[3]], {a, b, c}}]", "{a, b, c}");
    run("MapThread[#1 &, {NDArray[Range[3]], {a, b, c}}][[2]] == 2", "True");
}

static void test_mapthread_ndarray_pair(void) {
    run("MapThread[Plus, {NDArray[{1,2,3}], NDArray[{10,20,30}]}][[3]] == 33", "True");
}

/* ---------- MapAll ---------- */

static void test_mapall_ndarray_wraps_all_levels(void) {
    /* f wraps every level including the whole array at level 0 (head becomes f)
     * and each scalar leaf. */
    run("Head[MapAll[f, NDArray[Range[3]]]]", "f");
    run("Count[MapAll[f, NDArray[Range[3]]], _Real, Infinity]", "3");
    run("Count[MapAll[f, NDArray[Table[i+j,{i,2},{j,2}]]], _Real, Infinity]", "4");
}

/* ---------- Regression: plain lists still behave ---------- */

static void test_plain_list_unregressed(void) {
    run("Map[f, {a, b, c}]", "{f[a], f[b], f[c]}");
    run("Map[#^2 &, {1, 2, 3}]", "{1, 4, 9}");
    run("MapIndexed[f, {a, b, c}]", "{f[a, {1}], f[b, {2}], f[c, {3}]}");
    run("MapThread[f, {{1, 2}, {3, 4}}]", "{f[1, 3], f[2, 4]}");
    run("MapAll[f, {a, {b}}]", "f[{f[a], f[{f[b]}]}]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_map_ndarray_rank1);
    TEST(test_map_ndarray_rank2_rows);
    TEST(test_map_ndarray_numeric_repacks);
    TEST(test_map_ndarray_symbolic_is_list);
    TEST(test_map_ndarray_deep_level);
    TEST(test_mapindexed_ndarray_positions);
    TEST(test_mapindexed_ndarray_values);
    TEST(test_mapthread_ndarray_with_list);
    TEST(test_mapthread_ndarray_pair);
    TEST(test_mapall_ndarray_wraps_all_levels);
    TEST(test_plain_list_unregressed);

    printf("All Map-family NDArray tests passed!\n");
    return 0;
}
