/* Read-side Part and ordering functions on associations.
 *
 * Every expected value below was confirmed against Mathematica 15 (via
 * wolframscript), except that an unevaluated Part prints here as Part[...]
 * where Mathematica's InputForm writes expr[[...]].
 *
 * Unlike the abort-on-first-failure assert_eval_eq, this harness runs every
 * case, reports each failure, and exits nonzero if any failed. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

static int g_run = 0, g_failed = 0;

static void check(const char* input, const char* expected) {
    g_run++;
    Expr* parsed = parse_expression(input);
    if (!parsed) {
        fprintf(stderr, "FAIL (parse): %s\n", input);
        g_failed++;
        return;
    }
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    if (strcmp(str, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, str);
        g_failed++;
    }
    free(str);
    expr_free(evaluated);
}

#define A4 "<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3, \"d\" -> 4|>"

/* ---------- Part: sub-associations from All, Span and index lists ---------- */
static void test_part_span(void) {
    check(A4 "[[2;;3]]", "<|\"b\" -> 2, \"c\" -> 3|>");
    check(A4 "[[1;;1]]", "<|\"a\" -> 1|>");
    check(A4 "[[2;;]]", "<|\"b\" -> 2, \"c\" -> 3, \"d\" -> 4|>");
    check(A4 "[[;;;;2]]", "<|\"a\" -> 1, \"c\" -> 3|>");
    check(A4 "[[-2;;-1]]", "<|\"c\" -> 3, \"d\" -> 4|>");
    check(A4 "[[-3;;]]", "<|\"b\" -> 2, \"c\" -> 3, \"d\" -> 4|>");
    check(A4 "[[3;;2]]", "<||>");
    check(A4 "[[5;;4]]", "<||>");
    check(A4 "[[3;;1;;-1]]", "<|\"c\" -> 3, \"b\" -> 2, \"a\" -> 1|>");
    check(A4 "[[4;;1;;-2]]", "<|\"d\" -> 4, \"b\" -> 2|>");
    /* Off either end: unevaluated (Part::take). */
    check(A4 "[[2;;7]]", "Part[" A4 ", Span[2, 7]]");
    check(A4 "[[1;;0;;-1]]", "Part[" A4 ", Span[1, 0, -1]]");
}

static void test_part_lists(void) {
    check(A4 "[[{1,3}]]", "<|\"a\" -> 1, \"c\" -> 3|>");
    check(A4 "[[{-1,1}]]", "<|\"d\" -> 4, \"a\" -> 1|>");
    check(A4 "[[{\"a\",\"c\"}]]", "<|\"a\" -> 1, \"c\" -> 3|>");
    check("<|a -> 1, b -> 2|>[[{Key[b]}]]", "<|b -> 2|>");
    check(A4 "[[{\"c\", Key[\"a\"]}]]", "<|\"c\" -> 3, \"a\" -> 1|>");
    check(A4 "[[{}]]", "<||>");
    check(A4 "[[{1,1}]]", "<|\"a\" -> 1|>");
    check(A4 "[[{\"a\",\"z\"}]]", "<|\"a\" -> 1, \"z\" -> Missing[\"KeyAbsent\", \"z\"]|>");
    check("<|a -> 1|>[[{Key[z]}]]", "<|z -> Missing[\"KeyAbsent\", Key[z]]|>");
    /* Out of range, or positions mixed with keys: unevaluated. */
    check(A4 "[[{1,7}]]", "Part[" A4 ", {1, 7}]");
    check(A4 "[[{Key[\"b\"], 1}]]", "Part[" A4 ", {Key[\"b\"], 1}]");
}

static void test_part_all_and_deep(void) {
    check(A4 "[[All]]", A4);
    check("<|\"x\" -> {1, 2}, \"y\" -> {3, 4}|>[[All, 1]]", "<|\"x\" -> 1, \"y\" -> 3|>");
    check("<|\"x\" -> {1, 2}, \"y\" -> {3, 4}|>[[All, 2;;]]", "<|\"x\" -> {2}, \"y\" -> {4}|>");
    check("<|\"p\" -> <|\"x\" -> 1|>, \"q\" -> <|\"x\" -> 3|>|>[[All, \"x\"]]",
          "<|\"p\" -> 1, \"q\" -> 3|>");
    check("{<|a -> 1, b -> 2|>, <|a -> 3|>}[[All, Key[a]]]", "{1, 3}");
    check("{<|x -> 1|>, <|y -> 2|>}[[All, Key[x]]]", "{1, Missing[\"KeyAbsent\", Key[x]]}");
    check("<|\"x\" -> {1, 2}, \"y\" -> {3, 4}|>[[{2}, -1]]", "<|\"y\" -> 4|>");
    check("<|\"a\" -> {1, 2}|>[[1;;1, {2, 1}]]", "<|\"a\" -> {2, 1}|>");
    check(A4 "[[{1, -1}, 0]]", "<|\"a\" -> Integer, \"d\" -> Integer|>");
    check("<|\"a\" -> <|\"b\" -> {5}|>|>[[\"a\", \"b\", 1]]", "5");
    check("<|\"a\" -> {1, 2}|>[[{\"a\", \"z\"}, 1]]",
          "<|\"a\" -> 1, \"z\" -> Missing[\"KeyAbsent\", \"z\"]|>");
    /* A deeper index that fails on an atomic value: unevaluated (Part::partd). */
    check(A4 "[[All, 1]]", "Part[" A4 ", All, 1]");
    check(A4 "[[1, 1]]", "Part[" A4 ", 1, 1]");
}

static void test_part_single(void) {
    check(A4 "[[2]]", "2");
    check(A4 "[[-1]]", "4");
    check(A4 "[[\"c\"]]", "3");
    check(A4 "[[Key[\"c\"]]]", "3");
    check(A4 "[[0]]", "Association");
    check("<|a -> 1|>[[Key[b]]]", "Missing[\"KeyAbsent\", Key[b]]");
    check(A4 "[[\"z\"]]", "Missing[\"KeyAbsent\", \"z\"]");
    check(A4 "[[Key[\"z\"], 1]]", "Missing[\"KeyAbsent\", Key[\"z\"]]");
    check(A4 "[[5]]", "Part[" A4 ", 5]");
    check(A4 "[[-5]]", "Part[" A4 ", -5]");
    /* A bare symbol or real is not a part spec (Part::pkspec1). */
    check("<|x -> 1|>[[x]]", "Part[<|x -> 1|>, x]");
    check(A4 "[[1.5]]", "Part[" A4 ", 1.5]");
    /* RuleDelayed entries stay delayed in a sub-association. */
    check("<|\"a\" :> 1 + 1|>[[{1}]]", "<|\"a\" :> 1 + 1|>");
    check("<|\"a\" :> 1 + 1|>[[1]]", "2");
}

/* ---------- The accessor takes a literal key ---------- */
static void test_accessor_literal_key(void) {
    check("<|\"a\" -> 1|>[Key[\"a\"]]", "Missing[\"KeyAbsent\", Key[\"a\"]]");
    check("<|Key[\"a\"] -> 1|>[Key[\"a\"]]", "1");
    check("<|\"a\" -> <|\"b\" -> 1|>|>[\"a\", \"b\"]", "1");
    check("<|\"a\" -> <|\"b\" -> 1|>|>[\"z\", \"b\"]", "Missing[\"KeyAbsent\", \"z\"]");
    check("Key[\"a\"][<|\"a\" -> 1|>]", "1");
}

/* ---------- Sort / SortBy / Ordering / ReverseSort with ordering functions ---------- */
#define L4 "{{2, \"a\"}, {2, \"b\"}, {1, \"c\"}, {2, \"d\"}}"
static void test_sort_ordering_function(void) {
    check("Sort[<|\"a\" -> 2, \"b\" -> 3, \"c\" -> 1|>, Greater]",
          "<|\"b\" -> 3, \"a\" -> 2, \"c\" -> 1|>");
    check("Sort[<|\"a\" -> 2, \"b\" -> 3, \"c\" -> 1|>, #1 > #2 &]",
          "<|\"b\" -> 3, \"a\" -> 2, \"c\" -> 1|>");
    check("Sort[<|\"a\" -> 2, \"b\" -> 2, \"c\" -> 1, \"d\" -> 2|>, Greater]",
          "<|\"d\" -> 2, \"b\" -> 2, \"a\" -> 2, \"c\" -> 1|>");
    check("Sort[<|\"a\" -> 2, \"b\" -> 2, \"c\" -> 1, \"d\" -> 2|>, GreaterEqual]",
          "<|\"a\" -> 2, \"b\" -> 2, \"d\" -> 2, \"c\" -> 1|>");
    check("Sort[" L4 ", #1[[1]] > #2[[1]] &]", "{{2, \"d\"}, {2, \"b\"}, {2, \"a\"}, {1, \"c\"}}");
    check("Sort[" L4 ", -1 &]", "{{2, \"d\"}, {1, \"c\"}, {2, \"b\"}, {2, \"a\"}}");
    check("Sort[Range[10], Mod[#1, 3] > Mod[#2, 3] &]", "{8, 5, 2, 10, 7, 4, 1, 9, 6, 3}");
    check("Sort[Range[13], Mod[#1, 3] >= Mod[#2, 5] &]",
          "{13, 12, 9, 8, 7, 4, 3, 2, 1, 5, 6, 10, 11}");
    check("Ordering[<|\"a\" -> 2, \"b\" -> 2, \"c\" -> 1|>, All, Greater]", "{2, 1, 3}");
    check("Ordering[{1, 1, 0, 1}, All, Greater]", "{4, 2, 1, 3}");
    check("Ordering[<|\"a\" -> 2, \"b\" -> 3, \"c\" -> 1|>]", "{3, 1, 2}");
    check("Ordering[<|\"a\" -> 2, \"b\" -> 3, \"c\" -> 1|>, 2]", "{3, 1}");
    check("ReverseSort[<|\"a\" -> 2, \"b\" -> 2, \"c\" -> 1, \"d\" -> 2|>]",
          "<|\"a\" -> 2, \"b\" -> 2, \"d\" -> 2, \"c\" -> 1|>");
    check("ReverseSort[<|\"a\" -> 2, \"b\" -> 2, \"c\" -> 1, \"d\" -> 2|>, Less]",
          "<|\"a\" -> 2, \"b\" -> 2, \"d\" -> 2, \"c\" -> 1|>");
    check("ReverseSort[" L4 ", #1[[1]] < #2[[1]] &]",
          "{{2, \"a\"}, {2, \"b\"}, {2, \"d\"}, {1, \"c\"}}");
}

static void test_sortby(void) {
    check("SortBy[<|\"a\" -> {1, 2}, \"b\" -> {0, 5}, \"c\" -> {1, 1}|>, First, Greater]",
          "<|\"a\" -> {1, 2}, \"c\" -> {1, 1}, \"b\" -> {0, 5}|>");
    check("SortBy[<|\"a\" -> {1, 2}, \"b\" -> {0, 5}, \"c\" -> {1, 1}|>, First, Less]",
          "<|\"b\" -> {0, 5}, \"a\" -> {1, 2}, \"c\" -> {1, 1}|>");
    check("SortBy[{{1, \"a\"}, {1, \"b\"}, {0, \"c\"}, {1, \"d\"}}, First, Greater]",
          "{{1, \"d\"}, {1, \"b\"}, {1, \"a\"}, {0, \"c\"}}");
    check("SortBy[Range[12], Mod[#, 3] &, Greater]", "{11, 8, 5, 2, 10, 7, 4, 1, 12, 9, 6, 3}");
    check("SortBy[{{1, \"b\"}, {1, \"a\"}}, First]", "{{1, \"a\"}, {1, \"b\"}}");
    check("SortBy[{{1, \"b\"}, {1, \"a\"}}, {First}]", "{{1, \"b\"}, {1, \"a\"}}");
    check("SortBy[<|\"x\" -> {1, \"b\"}, \"y\" -> {1, \"a\"}|>, First]",
          "<|\"x\" -> {1, \"b\"}, \"y\" -> {1, \"a\"}|>");
    check("SortBy[<|\"a\" -> {1, \"z\"}, \"b\" -> {1, \"y\"}|>, First, Greater]",
          "<|\"a\" -> {1, \"z\"}, \"b\" -> {1, \"y\"}|>");
    check("ReverseSortBy[" L4 ", First, Less]", "{{2, \"d\"}, {2, \"b\"}, {2, \"a\"}, {1, \"c\"}}");
    check("ReverseSortBy[<|\"x\" -> {1, \"b\"}, \"y\" -> {1, \"a\"}, \"z\" -> {0, \"c\"}|>, First]",
          "<|\"y\" -> {1, \"a\"}, \"x\" -> {1, \"b\"}, \"z\" -> {0, \"c\"}|>");
}

/* ---------- KeySort / KeyTake ---------- */
static void test_keysort_keytake(void) {
    check("KeySort[<|\"c\" -> 1, \"a\" -> 2, \"b\" -> 3|>, Order]",
          "<|\"a\" -> 2, \"b\" -> 3, \"c\" -> 1|>");
    check("KeySort[<|\"c\" -> 1, \"a\" -> 2, \"b\" -> 3|>, -Order[#1, #2] &]",
          "<|\"c\" -> 1, \"b\" -> 3, \"a\" -> 2|>");
    check("KeySort[<|3 -> 1, 1 -> 2, 2 -> 3|>, Greater]", "<|3 -> 1, 2 -> 3, 1 -> 2|>");
    check("KeySort[<|\"c\" -> 1, \"a\" -> 2, \"b\" -> 3|>, Greater]",
          "<|\"c\" -> 1, \"a\" -> 2, \"b\" -> 3|>");
    check("KeyTake[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, {\"c\", \"a\"}]",
          "<|\"c\" -> 3, \"a\" -> 1|>");
    check("KeyTake[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, {\"c\", \"a\", \"z\"}]",
          "<|\"c\" -> 3, \"a\" -> 1|>");
    check("KeyTake[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, {\"c\", \"a\", \"c\"}]",
          "<|\"a\" -> 1, \"c\" -> 3|>");
    check("KeyTake[{<|\"a\" -> 1, \"b\" -> 2|>, <|\"b\" -> 3|>}, {\"b\", \"a\"}]",
          "{<|\"b\" -> 2, \"a\" -> 1|>, <|\"b\" -> 3|>}");
}

/* ---------- Prepend / Append / PrependTo / AppendTo ---------- */
static void test_prepend_append(void) {
    check("Prepend[<|\"a\" -> 1, \"b\" -> 2|>, \"b\" -> 9]", "<|\"b\" -> 9, \"a\" -> 1|>");
    check("Append[<|\"a\" -> 1, \"b\" -> 2|>, \"a\" -> 9]", "<|\"b\" -> 2, \"a\" -> 9|>");
    check("Prepend[<|\"a\" -> 1, \"b\" -> 2|>, {\"c\" -> 9, \"a\" -> 0}]",
          "<|\"c\" -> 9, \"a\" -> 0, \"b\" -> 2|>");
    check("Append[<|\"a\" -> 1, \"b\" -> 2|>, <|\"c\" -> 9, \"a\" -> 0|>]",
          "<|\"b\" -> 2, \"c\" -> 9, \"a\" -> 0|>");
    check("Prepend[<|\"a\" -> 1|>, 5]", "Prepend[<|\"a\" -> 1|>, 5]");
    check("Module[{x = <|\"a\" -> 1, \"b\" -> 2|>}, PrependTo[x, \"b\" -> 9]; x]",
          "<|\"b\" -> 9, \"a\" -> 1|>");
    check("Module[{x = <|\"a\" -> 1, \"b\" -> 2|>}, AppendTo[x, \"a\" -> 9]; x]",
          "<|\"b\" -> 2, \"a\" -> 9|>");
    check("Module[{x = <|\"a\" -> 1, \"b\" -> 2|>}, AppendTo[x, \"c\" -> 3]]",
          "<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>");
}

/* ---------- Catenate / Partition / Insert / Pick ---------- */
static void test_catenate_partition(void) {
    check("Catenate[{<|a -> 1|>, <|b -> 2|>}]", "{1, 2}");
    check("Catenate[{<|a -> 1, b -> 2|>, <|a -> 3|>}]", "{1, 2, 3}");
    check("Catenate[<|x -> <|a -> 1|>, y -> <|b -> 2|>|>]", "{1, 2}");
    check("Catenate[<|x -> {1, 2}, y -> {3}|>]", "{1, 2, 3}");
    check("Catenate[{<|a -> 1|>, {2, 3}}]", "{1, 2, 3}");
    check("Partition[" A4 ", 2]", "Partition[" A4 ", 2]");
    check("Partition[" A4 ", 2, 1]", "Partition[" A4 ", 2, 1]");
}

#define A3 "<|a -> 1, b -> 2, c -> 3|>"
static void test_insert(void) {
    check("Insert[" A3 ", d -> 3, Key[b]]", "<|a -> 1, d -> 3, b -> 2, c -> 3|>");
    check("Insert[" A3 ", d -> 3, 2]", "<|a -> 1, d -> 3, b -> 2, c -> 3|>");
    check("Insert[" A3 ", d -> 3, -1]", "<|a -> 1, b -> 2, c -> 3, d -> 3|>");
    check("Insert[" A3 ", d -> 3, -2]", "<|a -> 1, b -> 2, d -> 3, c -> 3|>");
    check("Insert[" A3 ", a -> 9, 3]", "<|b -> 2, a -> 9, c -> 3|>");
    check("Insert[" A3 ", a -> 3, Key[c]]", "<|b -> 2, a -> 3, c -> 3|>");
    check("Insert[" A3 ", d -> 3, 5]", "Insert[" A3 ", d -> 3, 5]");
    check("Insert[" A3 ", d -> 3, Key[z]]", "Insert[" A3 ", d -> 3, Key[z]]");
    check("Insert[" A3 ", 5, 2]", A3);
    check("Insert[" A3 ", {d -> 3, e -> 4}, 2]", "<|a -> 1, d -> 3, e -> 4, b -> 2, c -> 3|>");
    check("Insert[" A3 ", d -> 3, {{1}, {3}}]", "<|d -> 3, a -> 1, b -> 2, c -> 3|>");
    check("Insert[<|a -> {1, 2}|>, 9, {1, 1}]", "<|a -> {9, 1, 2}|>");
}

static void test_pick(void) {
    check("Pick[<|a -> 1, b -> 2|>, {True, False}]", "<|a -> 1|>");
    check("Pick[<|a -> 1, b -> 2, c -> 3|>, {1, 0, 1}, 1]", "<|a -> 1, c -> 3|>");
    check("Pick[<|a -> {1, 2}, b -> {3, 4}|>, {True, {False, True}}]",
          "<|a -> {1, 2}, b -> {4}|>");
    check("{Pick[<|a -> 1, b -> 2|>, <|a -> True, b -> False|>]}", "{}");
    check("Pick[<|a -> 1|>, True]", "<|a -> 1|>");
}

/* ---------- GroupBy with a list of functions ---------- */
static void test_groupby_multilevel(void) {
    check("GroupBy[{1, 2, 3, 4, 5, 6}, {EvenQ, # > 3 &}]",
          "<|False -> <|False -> {1, 3}, True -> {5}|>, True -> <|False -> {2}, True -> {4, 6}|>|>");
    check("GroupBy[{1, 2, 3, 4, 5, 6}, {EvenQ, # > 3 &}, Total]",
          "<|False -> <|False -> 4, True -> 5|>, True -> <|False -> 2, True -> 10|>|>");
    check("GroupBy[{1, 2, 3, 4, 5, 6}, {EvenQ}]", "<|False -> {1, 3, 5}, True -> {2, 4, 6}|>");
    check("GroupBy[{{1, x}, {2, y}, {1, z}}, {First, Last}, Length]",
          "<|1 -> <|x -> 1, z -> 1|>, 2 -> <|y -> 1|>|>");
    check("GroupBy[<|a -> 1, b -> 2, c -> 3|>, {OddQ, # > 1 &}]",
          "<|True -> <|False -> <|a -> 1|>, True -> <|c -> 3|>|>, False -> <|True -> <|b -> 2|>|>|>");
}

/* ---------- Structural heads on associations ---------- */
static void test_structural(void) {
    check("First[" A4 "]", "1");
    check("Last[" A4 "]", "4");
    check("Most[" A4 "]", "<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>");
    check("Rest[" A4 "]", "<|\"b\" -> 2, \"c\" -> 3, \"d\" -> 4|>");
    check("Take[" A4 ", {2, 3}]", "<|\"b\" -> 2, \"c\" -> 3|>");
    check("Drop[" A4 ", -1]", "<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>");
    check("Reverse[" A4 "]", "<|\"d\" -> 4, \"c\" -> 3, \"b\" -> 2, \"a\" -> 1|>");
    check("Reverse[<|\"x\" -> {1, 2}, \"y\" -> {3, 4}|>, 2]", "<|\"x\" -> {2, 1}, \"y\" -> {4, 3}|>");
    check("RotateLeft[" A4 ", 2]", "<|\"c\" -> 3, \"d\" -> 4, \"a\" -> 1, \"b\" -> 2|>");
    check("RotateRight[" A4 "]", "<|\"d\" -> 4, \"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>");
    check("TakeWhile[" A4 ", # < 3 &]", "<|\"a\" -> 1, \"b\" -> 2|>");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_part_span);
    TEST(test_part_lists);
    TEST(test_part_all_and_deep);
    TEST(test_part_single);
    TEST(test_accessor_literal_key);
    TEST(test_sort_ordering_function);
    TEST(test_sortby);
    TEST(test_keysort_keytake);
    TEST(test_prepend_append);
    TEST(test_catenate_partition);
    TEST(test_insert);
    TEST(test_pick);
    TEST(test_groupby_multilevel);
    TEST(test_structural);

    printf("assoc_read: %d/%d assertions passed\n", g_run - g_failed, g_run);
    return g_failed ? 1 : 0;
}
