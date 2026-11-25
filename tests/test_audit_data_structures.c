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

/* ================================================================== */
/* Wave 4, G2-structural (data structures): infix Rule printing, a    */
/* Throw from a callback ending the loop, RuleDelayed entries kept    */
/* delayed, rule lists for KeyTake/KeyDrop, entry patterns.           */
/* ================================================================== */

static void test_wave4_data_structures(void) {
    /* x //= f is ApplyTo */
    audit_chk("ApplyTo#1", "x4a = 5; x4a //= f", "f[5]");
    audit_chk("ApplyTo#1", "x4a", "f[5]");
    audit_chk("ApplyTo#1", "FullForm[Hold[y //= g; z]]", "Hold[CompoundExpression[ApplyTo[y, g], z]]");
    /* Rule prints with the parentheses its associativity needs */
    audit_chk("Rule#1", "ToExpression[ToString[InputForm[(a -> b) -> c]]] === ((a -> b) -> c)", "True");
    audit_chk("Rule#1", "ToString[InputForm[(a -> b) -> c]]", "\"(a -> b) -> c\"");
    audit_chk("Rule#1", "ToString[InputForm[a -> (b -> c)]]", "\"a -> b -> c\"");
    audit_chk("Rule#1", "ToString[InputForm[Hold[(a = b) = c]]]", "\"Hold[(a = b) = c]\"");
    audit_chk("Rule#1", "ToString[InputForm[Equal[Equal[a, b], c]]]", "\"(a == b) == c\"");
    audit_chk("Rule#1", "ToString[InputForm[Hold[(a /; b) /; c]]]", "\"Hold[a /; b /; c]\"");
    /* Rule with three arguments prints in functional form */
    audit_chk("Rule#2", "Quiet[ToString[InputForm[Rule[a, b, c]]]]", "\"Rule[a, b, c]\"");
    audit_chk("Rule#2", "Quiet[ToString[Rule[a, b, c]]]", "\"Rule[a, b, c]\"");
    /* AssociateTo merges an association */
    audit_chk("AssociateTo#1", "a4s = <|\"x\" -> 1|>; AssociateTo[a4s, <|\"y\" -> 2|>]", "<|\"x\" -> 1, \"y\" -> 2|>");
    audit_chk("AssociateTo#1", "a4s", "<|\"x\" -> 1, \"y\" -> 2|>");
    /* Delete[assoc, Key[k]] */
    audit_chk("Association#1", "Delete[<|a -> 1, b -> 2|>, Key[a]]", "<|b -> 2|>");
    audit_chk("Key#2", "Delete[<|\"a\"->1, \"b\"->2|>, Key[\"b\"]]", "<|\"a\" -> 1|>");
    audit_chk("Key#2", "Delete[<|\"a\"->1, \"b\"->2|>, {Key[\"b\"]}]", "<|\"a\" -> 1|>");
    /* RuleDelayed entries stay delayed */
    audit_chk("AssociationMap#1", "AssociationMap[Identity, <|a :> 1+1|>]", "<|a :> 1 + 1|>");
    audit_chk("AssociationMap#1", "AssociationMap[Reverse, <|a -> 1, b -> 2|>]", "<|1 -> a, 2 -> b|>");
    audit_chk("KeyIntersection#1", "KeyIntersection[{<|a :> 1 + 1, b -> 2|>, <|a -> 3|>}]", "{<|a :> 1 + 1|>, <|a -> 3|>}");
    audit_chk("KeyMap#1", "KeyMap[f, <|x :> 1 + 1|>]", "<|f[x] :> 1 + 1|>");
    audit_chk("KeyUnion#1", "KeyUnion[{<|a :> 1 + 1|>, <|b -> 2|>}]",
              "{<|a :> 1 + 1, b -> Missing[\"KeyAbsent\", b]|>, <|a -> Missing[\"KeyAbsent\", a], b -> 2|>}");
    /* CountDistinct with a test */
    audit_chk("CountDistinct#1", "CountDistinct[{1, 1.}, Equal]", "1");
    audit_chk("CountDistinct#1", "CountDistinct[{1, 2, 1}]", "2");
    /* a Throw from the key function ends the call at once */
    audit_chk("CountDistinctBy#1", "k4c = 0; Catch[CountDistinctBy[{1,2,3}, (k4c++; Throw[#])&]]", "1");
    audit_chk("CountDistinctBy#1", "k4c", "1");
    audit_chk("CountDistinctBy#1", "Catch[CountDistinctBy[{1,2,3}, Throw]]", "1");
    audit_chk("GatherBy#2", "k4g = 0; Catch[GatherBy[{1,2,3,4}, (k4g++; If[# > 1, Throw[#], 0])&]]", "2");
    audit_chk("GatherBy#2", "k4g", "2");
    /* DeleteDuplicatesBy is linear in the number of distinct keys */
    audit_chk("DeleteDuplicatesBy#1", "TimeConstrained[Length[DeleteDuplicatesBy[Range[100000], Identity]], 10]", "100000");
    audit_chk("DeleteDuplicatesBy#1", "DeleteDuplicatesBy[{1, 2, 3, 4, 5}, EvenQ]", "{1, 2}");
    /* KeyDrop / KeyTake over a list of rules */
    audit_chk("KeyDrop#1", "KeyDrop[{\"a\"->1,\"b\"->2},\"a\"]", "{\"b\" -> 2}");
    audit_chk("KeyTake#1", "KeyTake[{\"a\"->1, \"b\"->2, \"c\"->3}, {\"a\",\"c\"}]", "{\"a\" -> 1, \"c\" -> 3}");
    /* KeyMemberQ with sequence patterns */
    audit_chk("KeyMemberQ#1", "KeyMemberQ[<|f[1, 2] -> 0|>, f[PatternSequence[1, 2]]]", "True");
    audit_chk("KeyMemberQ#1", "KeyMemberQ[<|f[1, 2] -> 0|>, f[OrderlessPatternSequence[2, 1]]]", "True");
    audit_chk("KeyMemberQ#1", "KeyMemberQ[<|3 -> 0|>, Shortest[3]]", "True");
    audit_chk("KeyMemberQ#1", "KeyMemberQ[<|3 -> 0|>, 4]", "False");
    /* KeyValuePattern element patterns that are not literal rules */
    audit_chk("KeyValuePattern#1", "MatchQ[<|\"a\" -> 1|>, KeyValuePattern[{_Rule}]]", "True");
    audit_chk("KeyValuePattern#1", "MatchQ[<|\"a\" -> 1|>, KeyValuePattern[{(\"a\" -> 1) | (\"b\" -> 2)}]]", "True");
    audit_chk("KeyValuePattern#1", "MatchQ[<|\"a\" -> 1|>, KeyValuePattern[{(x_ -> y_) /; y > 0}]]", "True");
    audit_chk("KeyValuePattern#1", "MatchQ[<|\"a\" -> 1|>, KeyValuePattern[{(x_ -> y_) /; y > 5}]]", "False");
}

int main(void) {
    audit_init();
    TEST(test_gatherby_nested);
    TEST(test_key_accessor_contract);
    TEST(test_keysort_ordering_packed);
    TEST(test_s14_arrayq_deletemissing);
    TEST(test_wave3_arity_options);
    TEST(test_wave4_data_structures);
    return audit_finish("audit_data_structures");
}
