/* End-to-end tests for the Association (<| ... |>) data structure.
 *
 * Each assert_eval_eq drives the full pipeline: parse -> evaluate -> print,
 * so these exercise the parser (<| |> syntax), the evaluator/builtins, and
 * the printer together. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ---------- Literal syntax, construction, canonicalisation ---------- */

void test_assoc_literal_basic() {
    assert_eval_eq("<|\"a\" -> 1, \"b\" -> 2|>",
                   "<|\"a\" -> 1, \"b\" -> 2|>", 0);
}

void test_assoc_empty() {
    assert_eval_eq("<||>", "<||>", 0);
    assert_eval_eq("Association[]", "<||>", 0);
}

void test_assoc_head_form_equals_literal() {
    assert_eval_eq("Association[\"a\" -> 1, \"b\" -> 2]",
                   "<|\"a\" -> 1, \"b\" -> 2|>", 0);
}

void test_assoc_dedup_last_wins() {
    /* Duplicate keys collapse to the last value, first-occurrence order. */
    assert_eval_eq("<|\"a\" -> 1, \"b\" -> 2, \"a\" -> 99|>",
                   "<|\"a\" -> 99, \"b\" -> 2|>", 0);
}

void test_assoc_dedup_order_preserved() {
    assert_eval_eq("<|1 -> 10, 2 -> 20, 1 -> 30, 3 -> 40|>",
                   "<|1 -> 30, 2 -> 20, 3 -> 40|>", 0);
}

void test_assoc_splice_list_of_rules() {
    assert_eval_eq("Association[{\"a\" -> 1, \"b\" -> 2}]",
                   "<|\"a\" -> 1, \"b\" -> 2|>", 0);
}

void test_assoc_splice_nested_association() {
    assert_eval_eq("Association[<|\"a\" -> 1|>, <|\"b\" -> 2|>]",
                   "<|\"a\" -> 1, \"b\" -> 2|>", 0);
}

void test_assoc_symbol_keys() {
    assert_eval_eq("<|a -> 1, b -> 2|>", "<|a -> 1, b -> 2|>", 0);
}

void test_assoc_invalid_unevaluated() {
    /* A non-rule argument leaves Association unevaluated. */
    assert_eval_eq("Association[x]", "Association[x]", 0);
}

/* ---------- Keys / Values / Normal ---------- */

void test_keys_basic() {
    assert_eval_eq("Keys[<|\"a\" -> 1, \"b\" -> 2|>]", "{\"a\", \"b\"}", 0);
}

void test_values_basic() {
    assert_eval_eq("Values[<|\"a\" -> 1, \"b\" -> 2|>]", "{1, 2}", 0);
}

void test_keys_empty() {
    assert_eval_eq("Keys[<||>]", "{}", 0);
    assert_eval_eq("Values[<||>]", "{}", 0);
}

void test_normal_roundtrip() {
    assert_eval_eq("Normal[<|\"a\" -> 1, \"b\" -> 2|>]",
                   "{\"a\" -> 1, \"b\" -> 2}", 0);
    assert_eval_eq("Association[Normal[<|a -> 1, b -> 2|>]]",
                   "<|a -> 1, b -> 2|>", 0);
}

/* ---------- Lookup ---------- */

void test_lookup_present() {
    assert_eval_eq("Lookup[<|\"a\" -> 1, \"b\" -> 2|>, \"b\"]", "2", 0);
}

void test_lookup_missing() {
    assert_eval_eq("Lookup[<|\"a\" -> 1|>, \"z\"]",
                   "Missing[\"KeyAbsent\", \"z\"]", 0);
}

void test_lookup_default() {
    assert_eval_eq("Lookup[<|\"a\" -> 1|>, \"z\", 0]", "0", 0);
}

void test_lookup_list_of_keys() {
    assert_eval_eq("Lookup[<|\"a\" -> 1, \"b\" -> 2|>, {\"a\", \"b\", \"c\"}, -1]",
                   "{1, 2, -1}", 0);
}

/* ---------- KeyExistsQ / AssociationQ ---------- */

void test_keyexistsq() {
    assert_eval_eq("KeyExistsQ[<|\"a\" -> 1|>, \"a\"]", "True", 0);
    assert_eval_eq("KeyExistsQ[<|\"a\" -> 1|>, \"b\"]", "False", 0);
}

void test_associationq() {
    assert_eval_eq("AssociationQ[<|\"a\" -> 1|>]", "True", 0);
    assert_eval_eq("AssociationQ[<||>]", "True", 0);
    assert_eval_eq("AssociationQ[{1, 2, 3}]", "False", 0);
    assert_eval_eq("AssociationQ[42]", "False", 0);
}

/* ---------- Part-style access ---------- */

void test_part_string_key() {
    assert_eval_eq("<|\"a\" -> 10, \"b\" -> 20|>[[\"a\"]]", "10", 0);
}

void test_part_key_wrapper() {
    assert_eval_eq("<|a -> 10, b -> 20|>[[Key[b]]]", "20", 0);
}

void test_part_positional() {
    assert_eval_eq("<|\"a\" -> 10, \"b\" -> 20|>[[2]]", "20", 0);
    assert_eval_eq("<|\"a\" -> 10, \"b\" -> 20|>[[-1]]", "20", 0);
}

void test_part_missing_key() {
    assert_eval_eq("<|\"a\" -> 10|>[[\"z\"]]",
                   "Missing[\"KeyAbsent\", \"z\"]", 0);
}

void test_part_zero_gives_head() {
    /* assoc[[0]] is the head (as for any expression); an integer key needs Key[]. */
    assert_eval_eq("<|7 -> 5|>[[0]]", "Association", 0);
    assert_eval_eq("<|7 -> 5|>[[Key[7]]]", "5", 0);
}

/* ---------- KeyDrop / KeyTake ---------- */

void test_keydrop_single() {
    assert_eval_eq("KeyDrop[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, \"b\"]",
                   "<|\"a\" -> 1, \"c\" -> 3|>", 0);
}

void test_keydrop_list() {
    assert_eval_eq("KeyDrop[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, {\"a\", \"c\"}]",
                   "<|\"b\" -> 2|>", 0);
}

void test_keytake_list() {
    assert_eval_eq("KeyTake[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, {\"c\", \"a\"}]",
                   "<|\"a\" -> 1, \"c\" -> 3|>", 0);  /* association order preserved */
}

/* ---------- KeyValueMap ---------- */

void test_keyvaluemap() {
    assert_eval_eq("KeyValueMap[f, <|a -> 1, b -> 2|>]",
                   "{f[a, 1], f[b, 2]}", 0);
}

void test_keyvaluemap_plus() {
    assert_eval_eq("KeyValueMap[Plus, <|1 -> 10, 2 -> 20|>]",
                   "{11, 22}", 0);
}

/* ---------- AssociationThread ---------- */

void test_associationthread_two_lists() {
    assert_eval_eq("AssociationThread[{\"a\", \"b\"}, {1, 2}]",
                   "<|\"a\" -> 1, \"b\" -> 2|>", 0);
}

void test_associationthread_rule() {
    assert_eval_eq("AssociationThread[{\"a\", \"b\"} -> {1, 2}]",
                   "<|\"a\" -> 1, \"b\" -> 2|>", 0);
}

/* ---------- Counts ---------- */

void test_counts_basic() {
    assert_eval_eq("Counts[{1, 2, 2, 3, 3, 3}]",
                   "<|1 -> 1, 2 -> 2, 3 -> 3|>", 0);
}

void test_counts_strings() {
    assert_eval_eq("Counts[{\"a\", \"b\", \"a\"}]",
                   "<|\"a\" -> 2, \"b\" -> 1|>", 0);
}

void test_counts_empty() {
    assert_eval_eq("Counts[{}]", "<||>", 0);
}

/* ---------- GroupBy ---------- */

void test_groupby_parity() {
    assert_eval_eq("GroupBy[{1, 2, 3, 4, 5, 6}, EvenQ]",
                   "<|False -> {1, 3, 5}, True -> {2, 4, 6}|>", 0);
}

/* ---------- Merge ---------- */

void test_merge_total() {
    assert_eval_eq("Merge[{<|\"a\" -> 1|>, <|\"a\" -> 2, \"b\" -> 3|>}, Total]",
                   "<|\"a\" -> 3, \"b\" -> 3|>", 0);
}

void test_merge_list() {
    assert_eval_eq("Merge[{<|\"a\" -> 1|>, <|\"a\" -> 2|>}, Identity]",
                   "<|\"a\" -> {1, 2}|>", 0);
}

/* ---------- AssociateTo (in-place) ---------- */

void test_associate_to_add() {
    assert_eval_eq("(asc = <|\"a\" -> 1|>; AssociateTo[asc, \"b\" -> 2]; asc)",
                   "<|\"a\" -> 1, \"b\" -> 2|>", 0);
}

void test_associate_to_update() {
    assert_eval_eq("(asc2 = <|\"a\" -> 1, \"b\" -> 2|>; AssociateTo[asc2, \"a\" -> 99]; asc2)",
                   "<|\"a\" -> 99, \"b\" -> 2|>", 0);
}

/* ---------- Interaction with generic tools ---------- */

void test_length_of_association() {
    assert_eval_eq("Length[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>]", "3", 0);
}

void test_association_equality() {
    assert_eval_eq("<|\"a\" -> 1, \"b\" -> 2|> === <|\"a\" -> 1, \"b\" -> 2|>",
                   "True", 0);
}

void test_nested_association() {
    assert_eval_eq("<|\"outer\" -> <|\"inner\" -> 5|>|>[[\"outer\"]][[\"inner\"]]",
                   "5", 0);
}

/* ---------- Map / Select thread over values (Wolfram semantics) ---------- */

void test_map_over_values() {
    assert_eval_eq("Map[f, <|\"a\" -> 1, \"b\" -> 2|>]",
                   "<|\"a\" -> f[1], \"b\" -> f[2]|>", 0);
}

void test_map_square_values() {
    assert_eval_eq("Map[#^2 &, <|\"x\" -> 3, \"y\" -> 4|>]",
                   "<|\"x\" -> 9, \"y\" -> 16|>", 0);
}

void test_select_by_value() {
    assert_eval_eq("Select[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, # > 1 &]",
                   "<|\"b\" -> 2, \"c\" -> 3|>", 0);
}

/* ---------- KeySort / KeySortBy ---------- */

void test_keysort() {
    assert_eval_eq("KeySort[<|\"c\" -> 3, \"a\" -> 1, \"b\" -> 2|>]",
                   "<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>", 0);
}

void test_keysortby_length() {
    assert_eval_eq("KeySortBy[<|\"bbb\" -> 1, \"a\" -> 2, \"cc\" -> 3|>, StringLength]",
                   "<|\"a\" -> 2, \"cc\" -> 3, \"bbb\" -> 1|>", 0);
}

/* ---------- KeyMap / KeySelect ---------- */

void test_keymap() {
    assert_eval_eq("KeyMap[f, <|1 -> 10, 2 -> 20|>]",
                   "<|f[1] -> 10, f[2] -> 20|>", 0);
}

void test_keyselect() {
    assert_eval_eq("KeySelect[<|1 -> 10, 2 -> 20, 3 -> 30|>, EvenQ]",
                   "<|2 -> 20|>", 0);
}

/* ---------- CountsBy / PositionIndex / AssociationMap ---------- */

void test_countsby() {
    assert_eval_eq("CountsBy[Range[10], EvenQ]",
                   "<|False -> 5, True -> 5|>", 0);
}

void test_positionindex() {
    assert_eval_eq("PositionIndex[{a, b, a, c, a, b}]",
                   "<|a -> {1, 3, 5}, b -> {2, 6}, c -> {4}|>", 0);
}

void test_associationmap() {
    assert_eval_eq("AssociationMap[#^2 &, {1, 2, 3, 4}]",
                   "<|1 -> 1, 2 -> 4, 3 -> 9, 4 -> 16|>", 0);
}

/* ---------- Sort / Total / Min / Max thread over values ---------- */

void test_sort_by_value() {
    assert_eval_eq("Sort[<|\"a\" -> 3, \"b\" -> 1, \"c\" -> 2|>]",
                   "<|\"b\" -> 1, \"c\" -> 2, \"a\" -> 3|>", 0);
}

void test_total_of_values() {
    assert_eval_eq("Total[<|\"a\" -> 3, \"b\" -> 1, \"c\" -> 2|>]", "6", 0);
}

void test_total_exact_rationals() {
    assert_eval_eq("Total[<|\"a\" -> 1/2, \"b\" -> 1/3|>]", "5/6", 0);
}

void test_min_max_of_values() {
    assert_eval_eq("Min[<|\"a\" -> 3, \"b\" -> 1, \"c\" -> 2|>]", "1", 0);
    assert_eval_eq("Max[<|\"a\" -> 3, \"b\" -> 1, \"c\" -> 2|>]", "3", 0);
}

void test_join_merges_associations() {
    assert_eval_eq("Join[<|\"a\" -> 1, \"b\" -> 2|>, <|\"b\" -> 3, \"c\" -> 4|>]",
                   "<|\"a\" -> 1, \"b\" -> 3, \"c\" -> 4|>", 0);
}

/* ---------- Part assignment: assoc[[key]] = val ---------- */

void test_part_assign_update() {
    assert_eval_eq("(pa = <|\"x\" -> 1, \"y\" -> 2|>; pa[[\"x\"]] = 99; pa)",
                   "<|\"x\" -> 99, \"y\" -> 2|>", 0);
}

void test_part_assign_add_key() {
    assert_eval_eq("(pb = <|\"x\" -> 1|>; pb[[\"z\"]] = 5; pb)",
                   "<|\"x\" -> 1, \"z\" -> 5|>", 0);
}

void test_part_assign_key_wrapper() {
    assert_eval_eq("(pc = <|1 -> 10, 2 -> 20|>; pc[[Key[2]]] = 200; pc)",
                   "<|1 -> 10, 2 -> 200|>", 0);
}

void test_part_assign_positional() {
    assert_eval_eq("(pd = <|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>; pd[[2]] = 20; pd)",
                   "<|\"a\" -> 1, \"b\" -> 20, \"c\" -> 3|>", 0);
}

void test_part_assign_read_modify_write() {
    assert_eval_eq("(pe = <|\"x\" -> 1|>; pe[[\"x\"]] = pe[[\"x\"]] + 1; pe)",
                   "<|\"x\" -> 2|>", 0);
}

void test_part_assign_nested_multi_index() {
    assert_eval_eq("(pf = <|\"p\" -> <|\"q\" -> 1|>|>; pf[[\"p\", \"q\"]] = 42; pf)",
                   "<|\"p\" -> <|\"q\" -> 42|>|>", 0);
}

void test_part_assign_into_list_value() {
    assert_eval_eq("(pg = <|\"a\" -> {1, 2, 3}|>; pg[[\"a\", 2]] = 20; pg)",
                   "<|\"a\" -> {1, 20, 3}|>", 0);
}

/* ---------- Mean over values ---------- */

void test_mean_of_values() {
    assert_eval_eq("Mean[<|\"a\" -> 2, \"b\" -> 4, \"c\" -> 6|>]", "4", 0);
}

/* ---------- Cases / Count / DeleteCases over values ---------- */

void test_cases_over_values() {
    assert_eval_eq("Cases[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, x_ /; x > 1]",
                   "{2, 3}", 0);
}

void test_cases_pattern_head() {
    assert_eval_eq("Cases[<|\"a\" -> 1, \"b\" -> 2|>, _Integer]", "{1, 2}", 0);
}

void test_count_over_values() {
    assert_eval_eq("Count[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, x_ /; x > 1]", "2", 0);
}

void test_deletecases_over_values() {
    assert_eval_eq("DeleteCases[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, x_ /; x > 1]",
                   "<|\"a\" -> 1|>", 0);
}

/* ---------- AllTrue / AnyTrue / NoneTrue (lists and association values) ---------- */

void test_alltrue_list() {
    assert_eval_eq("AllTrue[{2, 4, 6}, EvenQ]", "True", 0);
    assert_eval_eq("AllTrue[{2, 3, 6}, EvenQ]", "False", 0);
    assert_eval_eq("AllTrue[{}, EvenQ]", "True", 0);
}

void test_anytrue_list() {
    assert_eval_eq("AnyTrue[{1, 3, 4}, EvenQ]", "True", 0);
    assert_eval_eq("AnyTrue[{1, 3, 5}, EvenQ]", "False", 0);
}

void test_nonetrue_list() {
    assert_eval_eq("NoneTrue[{1, 3, 5}, EvenQ]", "True", 0);
    assert_eval_eq("NoneTrue[{1, 2, 5}, EvenQ]", "False", 0);
}

void test_alltrue_over_values() {
    assert_eval_eq("AllTrue[<|\"a\" -> 2, \"b\" -> 4|>, EvenQ]", "True", 0);
    assert_eval_eq("AnyTrue[<|\"a\" -> 2, \"b\" -> 3|>, OddQ]", "True", 0);
    assert_eval_eq("NoneTrue[<|\"a\" -> 1, \"b\" -> 3|>, EvenQ]", "True", 0);
}

void test_alltrue_indeterminate_unevaluated() {
    assert_eval_eq("AllTrue[{2, x, 6}, # > 0 &]", "AllTrue[{2, x, 6}, #1 > 0 &]", 0);
}

void test_memberq_over_values() {
    assert_eval_eq("MemberQ[<|\"a\" -> 1, \"b\" -> 2|>, 2]", "True", 0);
    assert_eval_eq("MemberQ[<|\"a\" -> 1, \"b\" -> 2|>, 5]", "False", 0);
}

/* ---------- SortBy (lists and associations) ---------- */

void test_sortby_list_identity() {
    assert_eval_eq("SortBy[{5, 1, 3, 2}, Identity]", "{1, 2, 3, 5}", 0);
}

void test_sortby_list_abs() {
    assert_eval_eq("SortBy[{-3, 1, -2, 4}, Abs]", "{1, -2, -3, 4}", 0);
}

void test_sortby_list_stringlength() {
    assert_eval_eq("SortBy[{\"ccc\", \"a\", \"bb\"}, StringLength]",
                   "{\"a\", \"bb\", \"ccc\"}", 0);
}

void test_sortby_association() {
    assert_eval_eq("SortBy[<|\"a\" -> 3, \"b\" -> 1, \"c\" -> 2|>, Identity]",
                   "<|\"b\" -> 1, \"c\" -> 2, \"a\" -> 3|>", 0);
}

void test_sortby_association_by_first() {
    assert_eval_eq("SortBy[<|\"a\" -> {9}, \"b\" -> {1}|>, First]",
                   "<|\"b\" -> {1}, \"a\" -> {9}|>", 0);
}

void test_sortby_operator_form() {
    assert_eval_eq("SortBy[Abs][{-3, 1, -2}]", "{1, -2, -3}", 0);
}

void test_sortby_multi_criteria() {
    /* sort by First, break ties by Last */
    assert_eval_eq("SortBy[{{1, 3}, {1, 1}, {2, 0}, {1, 2}}, {First, Last}]",
                   "{{1, 1}, {1, 2}, {1, 3}, {2, 0}}", 0);
}

void test_sortby_multi_criteria_strings() {
    /* by length, then alphabetically */
    assert_eval_eq("SortBy[{\"bb\", \"a\", \"cc\", \"d\"}, {StringLength, Identity}]",
                   "{\"a\", \"d\", \"bb\", \"cc\"}", 0);
}

/* ---------- MaximalBy / MinimalBy (lists and associations) ---------- */

void test_maximalby_list_ties() {
    assert_eval_eq("MaximalBy[{1, -5, 3, -5, 2}, Abs]", "{-5, -5}", 0);
}

void test_minimalby_list_ties() {
    assert_eval_eq("MinimalBy[{3, 1, 4, 1, 5}, Identity]", "{1, 1}", 0);
}

void test_maximalby_association() {
    assert_eval_eq("MaximalBy[<|\"a\" -> 1, \"b\" -> 3, \"c\" -> 3|>, Identity]",
                   "<|\"b\" -> 3, \"c\" -> 3|>", 0);
}

void test_minimalby_association() {
    assert_eval_eq("MinimalBy[<|\"a\" -> 1, \"b\" -> 3, \"c\" -> 2|>, Identity]",
                   "<|\"a\" -> 1|>", 0);
}

void test_maximalby_operator_form() {
    assert_eval_eq("MaximalBy[Abs][{1, -5, 3}]", "{-5}", 0);
}

/* ---------- TakeLargest / TakeSmallest (+By), lists and associations ---------- */

void test_takelargest_list() {
    assert_eval_eq("TakeLargest[{3, 1, 4, 1, 5, 9, 2, 6}, 3]", "{9, 6, 5}", 0);
}

void test_takesmallest_list() {
    assert_eval_eq("TakeSmallest[{3, 1, 4, 1, 5, 9, 2, 6}, 3]", "{1, 1, 2}", 0);
}

void test_takelargest_association() {
    assert_eval_eq("TakeLargest[<|\"a\" -> 3, \"b\" -> 9, \"c\" -> 1, \"d\" -> 6|>, 2]",
                   "<|\"b\" -> 9, \"d\" -> 6|>", 0);
}

void test_takesmallest_association() {
    assert_eval_eq("TakeSmallest[<|\"a\" -> 3, \"b\" -> 9, \"c\" -> 1|>, 1]",
                   "<|\"c\" -> 1|>", 0);
}

void test_takelargestby_list() {
    assert_eval_eq("TakeLargestBy[{-9, 2, -3, 5}, Abs, 2]", "{-9, 5}", 0);
}

void test_takelargestby_association() {
    assert_eval_eq("TakeLargestBy[<|\"a\" -> {1}, \"b\" -> {9}, \"c\" -> {4}|>, First, 2]",
                   "<|\"b\" -> {9}, \"c\" -> {4}|>", 0);
}

void test_takelargest_more_than_length() {
    assert_eval_eq("TakeLargest[{1, 2}, 5]", "{2, 1}", 0);
}

/* ---------- GroupBy reducer form + GatherBy ---------- */

void test_groupby_reducer_total() {
    assert_eval_eq("GroupBy[Range[10], EvenQ, Total]",
                   "<|False -> 25, True -> 30|>", 0);
}

void test_groupby_reducer_length() {
    assert_eval_eq("GroupBy[{1, 2, 3, 4, 5, 6}, Mod[#, 3] &, Length]",
                   "<|1 -> 2, 2 -> 2, 0 -> 2|>", 0);
}

void test_groupby_reducer_mean() {
    assert_eval_eq("GroupBy[{1, 2, 3, 4}, EvenQ, Mean]",
                   "<|False -> 2, True -> 3|>", 0);
}

void test_groupby_key_value_transform() {
    assert_eval_eq("GroupBy[{{\"x\", 1}, {\"y\", 2}, {\"x\", 3}}, First -> Last]",
                   "<|\"x\" -> {1, 3}, \"y\" -> {2}|>", 0);
}

void test_groupby_key_value_reduce() {
    assert_eval_eq("GroupBy[{{\"x\", 1}, {\"y\", 2}, {\"x\", 3}}, First -> Last, Total]",
                   "<|\"x\" -> 4, \"y\" -> 2|>", 0);
}

void test_groupby_key_value_mean() {
    assert_eval_eq("GroupBy[{{\"a\", 5}, {\"a\", 7}, {\"b\", 2}}, First -> Last, Mean]",
                   "<|\"a\" -> 6, \"b\" -> 2|>", 0);
}

void test_gatherby_parity() {
    assert_eval_eq("GatherBy[{1, 2, 3, 4, 5, 6}, EvenQ]",
                   "{{1, 3, 5}, {2, 4, 6}}", 0);
}

void test_gatherby_stringlength() {
    assert_eval_eq("GatherBy[{\"aa\", \"b\", \"cc\", \"d\"}, StringLength]",
                   "{{\"aa\", \"cc\"}, {\"b\", \"d\"}}", 0);
}

/* ---------- ReverseSort / ReverseSortBy ---------- */

void test_reversesort_list() {
    assert_eval_eq("ReverseSort[{3, 1, 4, 1, 5, 9, 2}]", "{9, 5, 4, 3, 2, 1, 1}", 0);
}

void test_reversesort_association() {
    assert_eval_eq("ReverseSort[<|\"a\" -> 3, \"b\" -> 1, \"c\" -> 2|>]",
                   "<|\"a\" -> 3, \"c\" -> 2, \"b\" -> 1|>", 0);
}

void test_reversesortby_list() {
    assert_eval_eq("ReverseSortBy[{-9, 2, -3, 5}, Abs]", "{-9, 5, -3, 2}", 0);
}

void test_reversesortby_association() {
    assert_eval_eq("ReverseSortBy[<|\"a\" -> 1, \"b\" -> 9, \"c\" -> 4|>, Identity]",
                   "<|\"b\" -> 9, \"c\" -> 4, \"a\" -> 1|>", 0);
}

/* ---------- SelectFirst / FirstCase (lists and association values) ---------- */

void test_selectfirst_found() {
    assert_eval_eq("SelectFirst[{1, 3, 4, 5, 6}, EvenQ]", "4", 0);
}

void test_selectfirst_missing() {
    assert_eval_eq("SelectFirst[{1, 3, 5}, EvenQ]", "Missing[\"NotFound\"]", 0);
}

void test_selectfirst_default() {
    assert_eval_eq("SelectFirst[{1, 3, 5}, EvenQ, none]", "none", 0);
}

void test_selectfirst_association() {
    assert_eval_eq("SelectFirst[<|\"a\" -> 1, \"b\" -> 4, \"c\" -> 6|>, EvenQ]", "4", 0);
}

void test_firstcase_found() {
    assert_eval_eq("FirstCase[{1, 2, 3, 4}, _?EvenQ]", "2", 0);
}

void test_firstcase_default() {
    assert_eval_eq("FirstCase[{1, 3, 5}, _?EvenQ, none]", "none", 0);
}

void test_firstcase_association() {
    assert_eval_eq("FirstCase[<|\"a\" -> 1, \"b\" -> 2|>, _?EvenQ]", "2", 0);
}

void test_firstcase_missing() {
    assert_eval_eq("FirstCase[{a, b, c}, _Integer]", "Missing[\"NotFound\"]", 0);
}

/* ---------- DeleteMissing ---------- */

void test_deletemissing_list() {
    assert_eval_eq("DeleteMissing[{1, Missing[\"KeyAbsent\", \"x\"], 3, Missing[]}]",
                   "{1, 3}", 0);
}

void test_deletemissing_lookup_result() {
    assert_eval_eq("DeleteMissing[Lookup[<|\"a\" -> 1, \"b\" -> 2|>, {\"a\", \"z\", \"b\"}]]",
                   "{1, 2}", 0);
}

void test_deletemissing_association() {
    assert_eval_eq("DeleteMissing[<|\"a\" -> 1, \"b\" -> Missing[], \"c\" -> 3|>]",
                   "<|\"a\" -> 1, \"c\" -> 3|>", 0);
}

void test_deletemissing_noop() {
    assert_eval_eq("DeleteMissing[{1, 2, 3}]", "{1, 2, 3}", 0);
}

/* ---------- Integration: multi-builtin pipelines ---------- */

void test_pipeline_counts_takelargest() {
    /* word frequency -> top 2 */
    assert_eval_eq("TakeLargest[Counts[{a, b, a, c, a, b}], 2]",
                   "<|a -> 3, b -> 2|>", 0);
}

void test_pipeline_select_keys() {
    assert_eval_eq("Keys[Select[<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>, # > 1 &]]",
                   "{\"b\", \"c\"}", 0);
}

void test_pipeline_map_total() {
    assert_eval_eq("Total[Map[#^2 &, <|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>]]", "14", 0);
}

void test_pipeline_merge_counts() {
    assert_eval_eq("Merge[{Counts[{a, b}], Counts[{b, c}]}, Total]",
                   "<|a -> 1, b -> 2, c -> 1|>", 0);
}

void test_pipeline_groupby_reduce_rank() {
    /* group transactions by key, sum amounts, rank descending */
    assert_eval_eq("ReverseSort[GroupBy[{{\"x\", 1}, {\"y\", 2}, {\"x\", 3}}, First, Total[#[[All, 2]]] &]]",
                   "<|\"x\" -> 4, \"y\" -> 2|>", 0);
}

void test_pipeline_lookup_deletemissing_total() {
    assert_eval_eq("Total[DeleteMissing[Lookup[<|\"a\" -> 10, \"b\" -> 20|>, {\"a\", \"z\", \"b\"}]]]",
                   "30", 0);
}

/* ---------- Edge cases: empty-collection reductions match list behaviour ---------- */

void test_edge_empty_map_sort_keysort() {
    assert_eval_eq("Map[f, <||>]", "<||>", 0);
    assert_eval_eq("KeySort[<||>]", "<||>", 0);
    assert_eval_eq("GroupBy[{}, EvenQ]", "<||>", 0);
    assert_eval_eq("Merge[{}, Total]", "<||>", 0);
}

/* ---------- Append / Prepend on associations ---------- */

void test_append_association() {
    assert_eval_eq("Append[<|\"a\" -> 1, \"b\" -> 2|>, \"c\" -> 3]",
                   "<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>", 0);
}

void test_append_association_updates() {
    assert_eval_eq("Append[<|\"a\" -> 1|>, \"a\" -> 99]", "<|\"a\" -> 99|>", 0);
}

void test_prepend_association() {
    assert_eval_eq("Prepend[<|\"a\" -> 1|>, \"b\" -> 2]", "<|\"b\" -> 2, \"a\" -> 1|>", 0);
}

/* ---------- KeyValuePattern: destructuring associations in patterns ---------- */

void test_kvp_key_present() {
    assert_eval_eq("MatchQ[<|\"a\" -> 1, \"b\" -> 2|>, KeyValuePattern[{\"a\" -> _}]]", "True", 0);
}

void test_kvp_multiple_keys() {
    assert_eval_eq("MatchQ[<|\"a\" -> 1, \"b\" -> 2|>, KeyValuePattern[{\"a\" -> _, \"b\" -> _}]]", "True", 0);
}

void test_kvp_absent_key() {
    assert_eval_eq("MatchQ[<|\"a\" -> 1|>, KeyValuePattern[{\"z\" -> _}]]", "False", 0);
}

void test_kvp_value_constraint() {
    assert_eval_eq("MatchQ[<|\"a\" -> 1, \"b\" -> 2|>, KeyValuePattern[{\"b\" -> 2}]]", "True", 0);
    assert_eval_eq("MatchQ[<|\"a\" -> 1, \"b\" -> 2|>, KeyValuePattern[{\"b\" -> 3}]]", "False", 0);
}

void test_kvp_single_rule_form() {
    assert_eval_eq("MatchQ[<|\"a\" -> 1|>, KeyValuePattern[\"a\" -> _Integer]]", "True", 0);
}

void test_kvp_binding_extraction() {
    assert_eval_eq("Replace[<|\"a\" -> 5, \"b\" -> 2|>, KeyValuePattern[{\"a\" -> v_}] :> v]", "5", 0);
}

void test_kvp_on_list_of_rules() {
    assert_eval_eq("MatchQ[{\"a\" -> 1, \"b\" -> 2}, KeyValuePattern[{\"a\" -> _}]]", "True", 0);
}

void test_kvp_cases_record_filter() {
    assert_eval_eq("Cases[{<|\"t\" -> 1|>, <|\"t\" -> 2|>, <|\"x\" -> 3|>}, KeyValuePattern[{\"t\" -> _}]]",
                   "{<|\"t\" -> 1|>, <|\"t\" -> 2|>}", 0);
}

/* ---------- Association as accessor: assoc[key] ---------- */

void test_accessor_present() {
    assert_eval_eq("<|\"a\" -> 1, \"b\" -> 2|>[\"a\"]", "1", 0);
}

void test_accessor_missing() {
    assert_eval_eq("<|\"a\" -> 1|>[\"z\"]", "Missing[\"KeyAbsent\", \"z\"]", 0);
}

void test_accessor_key_wrapper() {
    assert_eval_eq("<|1 -> 10, 2 -> 20|>[Key[2]]", "20", 0);
}

void test_accessor_bound_symbol() {
    assert_eval_eq("(acc = <|\"x\" -> 42|>; acc[\"x\"])", "42", 0);
}

void test_accessor_in_map() {
    assert_eval_eq("Map[#[\"a\"] &, {<|\"a\" -> 1|>, <|\"a\" -> 2|>, <|\"a\" -> 3|>}]",
                   "{1, 2, 3}", 0);
}

void test_accessor_nested_multikey() {
    assert_eval_eq("<|\"a\" -> <|\"b\" -> 5|>|>[\"a\", \"b\"]", "5", 0);
}

void test_accessor_nested_first_key_only() {
    assert_eval_eq("<|\"a\" -> <|\"b\" -> 5|>|>[\"a\"]", "<|\"b\" -> 5|>", 0);
}

void test_accessor_nested_missing_first() {
    assert_eval_eq("<|\"a\" -> <|\"b\" -> 5|>|>[\"z\", \"b\"]",
                   "Missing[\"KeyAbsent\", \"z\"]", 0);
}

void test_accessor_tabular() {
    assert_eval_eq("(tab = <|\"r1\" -> <|\"x\" -> 1, \"y\" -> 2|>, \"r2\" -> <|\"x\" -> 3, \"y\" -> 4|>|>; tab[\"r2\", \"y\"])",
                   "4", 0);
}

/* ---------- KeyValuePattern in function definitions (DownValues) ---------- */

void test_kvp_downvalue_binding() {
    assert_eval_eq("(pdv[KeyValuePattern[{\"x\" -> x_}]] := got[x]; pdv[<|\"x\" -> 9, \"y\" -> 1|>])",
                   "got[9]", 0);
}

void test_kvp_downvalue_multikey() {
    assert_eval_eq("(area[KeyValuePattern[{\"w\" -> w_, \"h\" -> h_}]] := w h; area[<|\"w\" -> 3, \"h\" -> 4|>])",
                   "12", 0);
}

void test_kvp_downvalue_no_match() {
    assert_eval_eq("(pdv2[KeyValuePattern[{\"x\" -> _}]] := yes; pdv2[<|\"y\" -> 1|>])",
                   "pdv2[<|\"y\" -> 1|>]", 0);
}

void test_except_downvalue() {
    /* The same dispatch-filter fix also repairs Except in a DownValue LHS. */
    assert_eval_eq("(qdv[Except[0]] := nonzero; {qdv[5], qdv[0]})",
                   "{nonzero, qdv[0]}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_assoc_literal_basic);
    TEST(test_assoc_empty);
    TEST(test_assoc_head_form_equals_literal);
    TEST(test_assoc_dedup_last_wins);
    TEST(test_assoc_dedup_order_preserved);
    TEST(test_assoc_splice_list_of_rules);
    TEST(test_assoc_splice_nested_association);
    TEST(test_assoc_symbol_keys);
    TEST(test_assoc_invalid_unevaluated);

    TEST(test_keys_basic);
    TEST(test_values_basic);
    TEST(test_keys_empty);
    TEST(test_normal_roundtrip);

    TEST(test_lookup_present);
    TEST(test_lookup_missing);
    TEST(test_lookup_default);
    TEST(test_lookup_list_of_keys);

    TEST(test_keyexistsq);
    TEST(test_associationq);

    TEST(test_part_string_key);
    TEST(test_part_key_wrapper);
    TEST(test_part_positional);
    TEST(test_part_missing_key);
    TEST(test_part_zero_gives_head);

    TEST(test_keydrop_single);
    TEST(test_keydrop_list);
    TEST(test_keytake_list);

    TEST(test_keyvaluemap);
    TEST(test_keyvaluemap_plus);

    TEST(test_associationthread_two_lists);
    TEST(test_associationthread_rule);

    TEST(test_counts_basic);
    TEST(test_counts_strings);
    TEST(test_counts_empty);

    TEST(test_groupby_parity);

    TEST(test_merge_total);
    TEST(test_merge_list);

    TEST(test_associate_to_add);
    TEST(test_associate_to_update);

    TEST(test_length_of_association);
    TEST(test_association_equality);
    TEST(test_nested_association);

    TEST(test_map_over_values);
    TEST(test_map_square_values);
    TEST(test_select_by_value);
    TEST(test_keysort);
    TEST(test_keysortby_length);
    TEST(test_keymap);
    TEST(test_keyselect);
    TEST(test_countsby);
    TEST(test_positionindex);
    TEST(test_associationmap);

    TEST(test_sort_by_value);
    TEST(test_total_of_values);
    TEST(test_total_exact_rationals);
    TEST(test_min_max_of_values);
    TEST(test_join_merges_associations);

    TEST(test_part_assign_update);
    TEST(test_part_assign_add_key);
    TEST(test_part_assign_key_wrapper);
    TEST(test_part_assign_positional);
    TEST(test_part_assign_read_modify_write);
    TEST(test_part_assign_nested_multi_index);
    TEST(test_part_assign_into_list_value);
    TEST(test_mean_of_values);

    TEST(test_cases_over_values);
    TEST(test_cases_pattern_head);
    TEST(test_count_over_values);
    TEST(test_deletecases_over_values);

    TEST(test_alltrue_list);
    TEST(test_anytrue_list);
    TEST(test_nonetrue_list);
    TEST(test_alltrue_over_values);
    TEST(test_alltrue_indeterminate_unevaluated);
    TEST(test_memberq_over_values);

    TEST(test_sortby_list_identity);
    TEST(test_sortby_list_abs);
    TEST(test_sortby_list_stringlength);
    TEST(test_sortby_association);
    TEST(test_sortby_association_by_first);
    TEST(test_sortby_operator_form);
    TEST(test_sortby_multi_criteria);
    TEST(test_sortby_multi_criteria_strings);

    TEST(test_maximalby_list_ties);
    TEST(test_minimalby_list_ties);
    TEST(test_maximalby_association);
    TEST(test_minimalby_association);
    TEST(test_maximalby_operator_form);

    TEST(test_takelargest_list);
    TEST(test_takesmallest_list);
    TEST(test_takelargest_association);
    TEST(test_takesmallest_association);
    TEST(test_takelargestby_list);
    TEST(test_takelargestby_association);
    TEST(test_takelargest_more_than_length);

    TEST(test_groupby_reducer_total);
    TEST(test_groupby_reducer_length);
    TEST(test_groupby_reducer_mean);
    TEST(test_groupby_key_value_transform);
    TEST(test_groupby_key_value_reduce);
    TEST(test_groupby_key_value_mean);
    TEST(test_gatherby_parity);
    TEST(test_gatherby_stringlength);

    TEST(test_reversesort_list);
    TEST(test_reversesort_association);
    TEST(test_reversesortby_list);
    TEST(test_reversesortby_association);

    TEST(test_selectfirst_found);
    TEST(test_selectfirst_missing);
    TEST(test_selectfirst_default);
    TEST(test_selectfirst_association);
    TEST(test_firstcase_found);
    TEST(test_firstcase_default);
    TEST(test_firstcase_association);
    TEST(test_firstcase_missing);

    TEST(test_deletemissing_list);
    TEST(test_deletemissing_lookup_result);
    TEST(test_deletemissing_association);
    TEST(test_deletemissing_noop);

    TEST(test_pipeline_counts_takelargest);
    TEST(test_pipeline_select_keys);
    TEST(test_pipeline_map_total);
    TEST(test_pipeline_merge_counts);
    TEST(test_pipeline_groupby_reduce_rank);
    TEST(test_pipeline_lookup_deletemissing_total);
    TEST(test_edge_empty_map_sort_keysort);

    TEST(test_append_association);
    TEST(test_append_association_updates);
    TEST(test_prepend_association);
    TEST(test_kvp_key_present);
    TEST(test_kvp_multiple_keys);
    TEST(test_kvp_absent_key);
    TEST(test_kvp_value_constraint);
    TEST(test_kvp_single_rule_form);
    TEST(test_kvp_binding_extraction);
    TEST(test_kvp_on_list_of_rules);
    TEST(test_kvp_cases_record_filter);

    TEST(test_accessor_present);
    TEST(test_accessor_missing);
    TEST(test_accessor_key_wrapper);
    TEST(test_accessor_bound_symbol);
    TEST(test_accessor_in_map);
    TEST(test_accessor_nested_multikey);
    TEST(test_accessor_nested_first_key_only);
    TEST(test_accessor_nested_missing_first);
    TEST(test_accessor_tabular);

    TEST(test_kvp_downvalue_binding);
    TEST(test_kvp_downvalue_multikey);
    TEST(test_kvp_downvalue_no_match);
    TEST(test_except_downvalue);

    printf("All Association tests passed.\n");
    return 0;
}
