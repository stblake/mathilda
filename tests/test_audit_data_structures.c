/* Builtin-audit regressions: data structures (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* GatherBy#1 (high): "GatherBy[list, {f1, f2, ...}] nested gathering is not
 * supported; returns singleton groups". */
static void test_gatherby_nested(void) {
    audit_chk("GatherBy#1", "GatherBy[Range[6], {EvenQ}]", "{{1, 3, 5}, {2, 4, 6}}");
    audit_chk("GatherBy#1", "GatherBy[Range[12], {EvenQ, Mod[#, 3] &}]",
              "{{{1, 7}, {3, 9}, {5, 11}}, {{2, 8}, {4, 10}, {6, 12}}}");
    audit_true("GatherBy#1",
               "GatherBy[Range[24], {EvenQ, Mod[#, 3] &, Mod[#, 4] &}] === "
               "Map[GatherBy[#, Mod[#, 4] &] &, GatherBy[Range[24], {EvenQ, Mod[#, 3] &}], {2}]");
    audit_chk("GatherBy#1", "GatherBy[{}, {EvenQ, OddQ}]", "{}");
    audit_chk("GatherBy#1", "GatherBy[<|a -> 1, b -> 2, c -> 3|>, {OddQ}]", "{<|a -> 1, c -> 3|>, <|b -> 2|>}");
    /* the single-function form is unchanged */
    audit_chk("GatherBy#1", "GatherBy[Range[6], EvenQ]", "{{1, 3, 5}, {2, 4, 6}}");
}

/* Key#1 (high): "assoc[Key[k]] does not strip Key: returns
 * Missing["KeyAbsent", Key[k]] for a present key".  NOT a defect: the
 * accessor's argument is the literal key, as in the reference
 * (<|"b" -> 2|>[Key["b"]] is Missing["KeyAbsent", Key["b"]]); the spec page
 * claimed otherwise and was corrected.  Pin the contract and the unwrapping
 * forms. */
static void test_key_accessor_contract(void) {
    audit_chk("Key#1", "<|\"a\" -> 1, \"b\" -> 2|>[Key[\"b\"]]", "Missing[\"KeyAbsent\", Key[\"b\"]]");
    audit_chk("Key#1", "<|Key[\"b\"] -> 5|>[Key[\"b\"]]", "5");
    audit_chk("Key#1", "<|\"a\" -> 1, \"b\" -> 2|>[[Key[\"b\"]]]", "2");
    audit_chk("Key#1", "Lookup[<|{1, 2} -> \"pair\"|>, Key[{1, 2}]]", "\"pair\"");
    audit_chk("Key#1", "Key[\"b\"][<|\"a\" -> 1, \"b\" -> 2|>]", "2");
}

/* KeySort#1 (high): "KeySort[assoc, p] stays unevaluated for any association
 * with 4 or more keys" (Ordering's packed result was not read). */
static void test_keysort_ordering_packed(void) {
    audit_chk("KeySort#1", "KeySort[<|4 -> 1, 1 -> 2, 3 -> 3, 2 -> 4|>, Greater]",
              "<|4 -> 1, 3 -> 3, 2 -> 4, 1 -> 2|>");
    audit_chk("KeySort#1", "KeySort[<|5 -> a, 1 -> b, 4 -> c, 2 -> d, 3 -> e|>, Order]",
              "<|1 -> b, 2 -> d, 3 -> e, 4 -> c, 5 -> a|>");
    audit_true("KeySort#1",
               "Keys[KeySort[AssociationThread[Range[50] -> Range[50]], Greater]] === Range[50, 1, -1]");
    audit_chk("KeySort#1",
              "KeySort[<|\"bb\" -> 1, \"a\" -> 2, \"dddd\" -> 3, \"ccc\" -> 4|>, StringLength[#1] <= StringLength[#2] &]",
              "<|\"a\" -> 2, \"bb\" -> 1, \"ccc\" -> 4, \"dddd\" -> 3|>");
    audit_chk("KeySort#1", "KeySort[<|3 -> a, 1 -> b, 2 -> c|>, #1 > #2 &]", "<|3 -> a, 2 -> c, 1 -> b|>");
}

/* ArrayQ#1 / DeleteMissing#1 (low): no depth cap; a BigInt level. */
static void test_s14_arrayq_deletemissing(void) {
    audit_chk("ArrayQ#1", "ArrayQ[Nest[List, 1, 33]]", "True");
    audit_chk("ArrayQ#1", "ArrayQ[Nest[List, 1, 40]]", "True");
    audit_chk("ArrayQ#1", "ArrayQ[Nest[List, 1, 40], 40]", "True");
    audit_chk("ArrayQ#1", "ArrayQ[{{1, 2}, {3}}]", "False");
    audit_chk("ArrayQ#1", "ArrayQ[{1, {2}}]", "False");
    audit_true("DeleteMissing#1", "DeleteMissing[{1, Missing[]}, 2^70] === {1}");
    audit_true("DeleteMissing#1", "DeleteMissing[{{1, Missing[]}, Missing[]}, 2^70] === {{1}}");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Rule#3", "Quiet[Check[Rule[a], fired, Rule::argr]]", "fired");
    audit_chk("Rule#3", "Quiet[Rule[a]]", "Rule[a]");
    audit_chk("ApplyTo#3", "Quiet[Check[ApplyTo[x], fired, ApplyTo::argr]]", "fired");
    audit_chk("ApplyTo#3", "Quiet[ApplyTo[x]]", "ApplyTo[x]");
    audit_chk("CountDistinct#2", "Quiet[Check[CountDistinct[], fired, CountDistinct::argt]]", "fired");
    audit_chk("CountDistinct#2", "Quiet[CountDistinct[]]", "CountDistinct[]");
    audit_chk("CountsBy#1", "Quiet[Check[CountsBy[{1, 2}, f, g], fired, CountsBy::argt]]", "fired");
    audit_chk("CountsBy#1", "Quiet[CountsBy[{1, 2}, f, g]]", "CountsBy[{1, 2}, f, g]");
    audit_chk("GroupBy#1", "Quiet[Check[GroupBy[{1, 2, 3}, EvenQ, Total, 4], fired, GroupBy::argb]]", "fired");
    audit_chk("GroupBy#1", "Quiet[GroupBy[{1, 2, 3}, EvenQ, Total, 4]]", "GroupBy[{1, 2, 3}, EvenQ, Total, 4]");
    audit_chk("KeyComplement#1", "Quiet[Check[KeyComplement[], fired, KeyComplement::argx]]", "fired");
    audit_chk("KeyComplement#1", "Quiet[KeyComplement[]]", "KeyComplement[]");
    audit_chk("KeyValueMap#1", "Quiet[Check[KeyValueMap[f, <|a -> 1|>, g], fired, KeyValueMap::argt]]", "fired");
    audit_chk("Merge#1", "Quiet[Check[Merge[{<|a -> 1|>}, f, g], fired, Merge::argt]]", "fired");
    audit_chk("MissingQ#1", "Quiet[Check[MissingQ[Sequence[1, 2]], fired, MissingQ::argx]]", "fired");
    audit_chk("MissingQ#1", "Quiet[MissingQ[Sequence[1, 2]]]", "MissingQ[1, 2]");
    audit_chk("PositionIndex#1", "Quiet[Check[PositionIndex[{a, b}, {c}], fired, PositionIndex::argx]]", "fired");
    audit_chk("PositionIndex#1", "Quiet[PositionIndex[{a, b}, {c}]]", "PositionIndex[{a, b}, {c}]");
}

int main(void) {
    audit_init();
    TEST(test_gatherby_nested);
    TEST(test_key_accessor_contract);
    TEST(test_keysort_ordering_packed);
    TEST(test_s14_arrayq_deletemissing);
    TEST(test_wave3_arity_options);
    return audit_finish("audit_data_structures");
}
