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

    printf("All Association tests passed.\n");
    return 0;
}
