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

    printf("All Association tests passed.\n");
    return 0;
}
