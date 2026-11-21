/* LinkedList through the sorting, set, grouping/counting, Association and
 * restructuring heads (src/sort.c, src/list/*, src/assoc*.c, src/replace.c,
 * src/random.c), each run on the slots -- through an element-array view
 * (src/linkedlist_seq.h) or a typed list's raw-buffer kernel -- never on a
 * List built from the list.
 *
 * The core check, llchk[h, l] (as in test_linkedlist.c): h on a LinkedList of
 * l agrees with h[Normal[ll]] -- held, so an unevaluated h[ll] is never
 * rescued by evaluating h[{...}] afterwards -- without materialising the list
 * (only Normal may) and without changing it. Normal of a long typed list is
 * PACKED, so the long machine fixtures pin the answers bit-for-bit to the
 * packed List's (the raw kernels are the ones the packed path runs).
 */
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "parse.h"
#include "eval.h"
#include "linkedlist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char* LLCHK =
    "llchk[h_, l0_] := Module[{ll = LinkedList[l0], l, c0, r, c1}, l = Normal[ll]; "
    "c0 = Internal`LinkedListConversionCount[]; r = h[ll]; "
    "c1 = Internal`LinkedListConversionCount[]; "
    "{c1 == c0, (Hold[Evaluate[r]] //. x_LinkedList :> With[{v = Normal[x]}, v /; True]) === "
    "Hold[Evaluate[h[l]]], Normal[ll] === l}]; ";

/* llchk with nested LinkedLists in the fixture: the reference is h of the
 * fully Normal nest, and both answers are compared with every LinkedList
 * replaced by its elements. */
static const char* LLCHK2 =
    "dnll[e_] := e //. x_LinkedList :> With[{v = Normal[x]}, v /; True]; "
    "llchk2[h_, l0_] := Module[{ll = LinkedList[l0], c0, r, c1}, "
    "c0 = Internal`LinkedListConversionCount[]; r = h[ll]; "
    "c1 = Internal`LinkedListConversionCount[]; "
    "{c1 == c0, dnll[Hold[Evaluate[r]]] === dnll[Hold[Evaluate[h[dnll[l0]]]]]}]; ";

static void expect_llchk2(const char* heads, const char* fixtures) {
    char buf[8192];
    snprintf(buf, sizeof buf,
             "%sUnion[Flatten[Quiet[Table[llchk2[h, l], {h, %s}, {l, %s}]]]]", LLCHK2, heads, fixtures);
    assert_eval_eq(buf, "{True}", 0);
}

static void expect_llchk(const char* heads, const char* fixtures) {
    char buf[8192];
    snprintf(buf, sizeof buf,
             "%sUnion[Flatten[Quiet[Table[llchk[h, l], {h, %s}, {l, %s}]]]]", LLCHK, heads, fixtures);
    assert_eval_eq(buf, "{True}", 0);
}

/* Small mixed fixtures, plus long typed ones whose Normal is packed. */
#define FX_SMALL "{{3, 1, 2, 2, -7}, {1.5, -2., 3., 0.5, 3.}, {s1a, s1b, s1a, s1c}, " \
                 "{\"b\", \"a\", \"c\", \"a\"}, {3, s1a, 2.5, s1a, 1/2}, {}}"
#define FX_LONG  "{Mod[Range[3000]^2, 97] - 40, N[Sin[Range[3000]]], Range[2000, 1, -1]}"

/* Sorting and ordering: typed lists through the packed Sort / Ordering kernels. */
void test_sorting() {
    const char* heads =
        "{Sort, Sort[#, Greater] &, ReverseSort, SortBy[#, Minus] &, SortBy[#, {Sign, Abs}] &, "
        "SortBy[#, Abs, Greater] &, ReverseSortBy[#, Abs] &, ReverseSortBy[#, Abs, Less] &, Ordering, "
        "Ordering[#, 2] &, Ordering[#, -2] &, Ordering[#, All, Greater] &, MaximalBy[#, Abs] &, "
        "MinimalBy[#, Abs] &, TakeLargest[#, 2] &, TakeSmallest[#, 2] &, RankedMax[#, 2] &, "
        "RankedMin[#, 1] &}";
    expect_llchk(heads, FX_SMALL);
    expect_llchk(heads, FX_LONG);
}

/* Set operations: LinkedList operands mix with Lists (the answer is a
 * LinkedList when any operand is one). */
void test_sets() {
    const char* heads =
        "{Union, Union[#, #] &, Union[#, {2, s1a}] &, Union[#, SameTest -> (Abs[#1] == Abs[#2] &)] &, "
        "Intersection[#, {2, 3., s1a, \"a\"}] &, Intersection[#, #] &, "
        "Complement[#, {2, s1b}] &, Complement[{1, 2, 3, s1a}, #] &, IntersectingQ[#, {2, s1c}] &, "
        "SubsetQ[#, {2}] &, SubsetQ[{1, 2, 3}, #] &, DeleteDuplicates, "
        "DeleteDuplicates[#, Abs[#1] == Abs[#2] &] &, DeleteDuplicatesBy[#, Abs] &, DuplicateFreeQ, "
        "DuplicateFreeQ[#, Abs[#1] == Abs[#2] &] &, Tally, Tally[#, Abs[#1] == Abs[#2] &] &, "
        "Commonest, Commonest[#, 2] &, CountDistinct, CountDistinctBy[#, Abs] &}";
    expect_llchk(heads, FX_SMALL);
    expect_llchk("{Union, Intersection[#, Range[50]] &, Complement[#, Range[-10, 10]] &, DeleteDuplicates, "
                 "Tally, Commonest, CountDistinct}", FX_LONG);
    /* Mixing a LinkedList with a List operand. */
    assert_eval_eq("Normal[Union[LinkedList[{3, 1}], {2, 1}]]", "{1, 2, 3}", 0);
    assert_eval_eq("Head[Complement[{1, 2, 3}, LinkedList[{2}]]]", "LinkedList", 0);
    assert_eval_eq("IntersectingQ[LinkedList[{s2a, s2b}], LinkedList[{s2b}]]", "True", 0);
}

/* Grouping, counting and the Association family. */
void test_grouping_and_assoc() {
    expect_llchk("{Gather, Gather[#, Abs[#1] == Abs[#2] &] &, GatherBy[#, Abs] &, Split, "
                 "Split[#, Abs[#1] == Abs[#2] &] &, SplitBy[#, Abs] &, SplitBy[#, {Sign, Abs}] &, "
                 "Counts, CountsBy[#, Abs] &, PositionIndex, GroupBy[#, Abs] &, GroupBy[#, Abs, Length] &, "
                 "GroupBy[#, {Sign, Abs}] &, GroupBy[#, Abs -> Minus] &, AssociationThread[#, #] &, "
                 "AssociationThread[# -> 1] &, AssociationMap[s3f, #] &}", FX_SMALL);
    /* A LinkedList of rules is a rule list to the key heads and to Association. */
    assert_eval_eq("s3r = LinkedList[{s3a -> 1, s3b -> 2}]; {Lookup[s3r, s3b], Lookup[s3r, s3c, 0], "
                   "Lookup[s3r, {s3a, s3c}], KeyExistsQ[s3r, s3a], KeyMemberQ[s3r, _], KeyFreeQ[s3r, s3c], "
                   "Association[s3r], Association[s3r, s3c -> 3]}",
                   "{2, 0, {1, Missing[\"KeyAbsent\", s3c]}, True, True, True, <|s3a -> 1, s3b -> 2|>, "
                   "<|s3a -> 1, s3b -> 2, s3c -> 3|>}", 0);
    /* Lookup threads over a LinkedList of associations, and of keys. */
    assert_eval_eq("Normal[Lookup[LinkedList[{<|s3k -> 1|>, <|s3k -> 2|>}], s3k]]", "{1, 2}", 0);
    assert_eval_eq("Normal[Lookup[<|s3k -> 1|>, LinkedList[{s3k, s3j}]]]", "{1, Missing[\"KeyAbsent\", s3j]}", 0);
    /* GroupBy values and GatherBy groups are LinkedLists for a LinkedList. */
    assert_eval_eq("Head /@ Values[GroupBy[LinkedList[{1, 2, 3}], OddQ]]", "{LinkedList, LinkedList}", 0);
}

/* Restructuring: sub-sequences are LinkedLists inside an outer LinkedList. */
void test_restructuring() {
    expect_llchk("{Partition[#, 2] &, Partition[#, 2, 1] &, Partition[#, 3, 1, 1] &, Partition[#, UpTo[2]] &, "
                 "Partition[#, 2, 1, {1, 1}, {}] &, Subsets, Subsets[#, {2}] &, Subsets[#, All, 3] &, "
                 "Permutations[#, 2] &, Tuples[#, 2] &, Tuples[#, {1, 2}] &, Riffle[#, 0] &, Riffle[#, {s4x, s4y}] &, "
                 "RotateLeft, RotateRight[#, 2] &, PadLeft[#, 7] &, PadRight[#, 7, s4p] &, PadRight[#, 2] &, "
                 "PadLeft[#, {7}, {s4p, s4q}] &, Flatten, ArrayPad[#, 2] &, ArrayPad[#, {1, 2}, s4p] &, "
                 "ArrayReshape[#, {2, 2}] &, ArrayRules, Accumulate, Differences, Differences[#, 2] &, "
                 "Differences[#, 1, 2] &, Ratios, Ratios[#, 2] &, Max, Min, MinMax, Max[#, 2] &, "
                 "Nearest[#, 2] &, ComposeList[{s4f, s4g}, #] &}", FX_SMALL);
    expect_llchk("{Partition[#, 7] &, Accumulate, Differences, Differences[#, 3] &, Ratios, Max, Min, MinMax, "
                 "Flatten, RotateLeft[#, 5] &}", FX_LONG);
    /* Nested LinkedLists: rows are read from their slots, at every level. */
    expect_llchk2("{Flatten, Flatten[#, 1] &, Catenate, ArrayFlatten[#, 1] &, ArrayReshape[#, {3}] &, "
                 "PadRight, PadLeft[#, {3, 3}] &, ArrayPad[#, 1] &, Tuples, Max, Partition[#, 1] &, "
                 "Tally}",
                 "{{LinkedList[{1, 2}], LinkedList[{3}]}, {LinkedList[{s4a}], {s4b, s4c}}}");
    assert_eval_eq("Normal[ComposeList[LinkedList[{s4f, s4g}], s4x]]", "{s4x, s4f[s4x], s4g[s4f[s4x]]}", 0);
    assert_eval_eq("Head /@ Normal[Partition[LinkedList[{1, 2, 3, 4}], 2]]", "{LinkedList, LinkedList}", 0);
    /* Sub-sequences are small LinkedLists: header plus one chunk, on the
     * order of a List node (the chunk directory is sized to the list). */
    assert_eval_eq("s4s = Partition[LinkedList[Range[10]], 2][[1]]; "
                   "{Head[s4s], Normal[s4s], ByteCount[s4s] < 400}", "{LinkedList, {1, 2}, True}", 0);
    /* Fixed-shape records stay Lists. */
    assert_eval_eq("{Head[MinMax[LinkedList[{3, 1}]]], Head[First[Tally[LinkedList[{3, 1}]]]]}",
                   "{List, List}", 0);
}

/* Machine answers agree bit-for-bit with the packed List; arithmetic that
 * leaves the machine range takes the element path, as the packed List does. */
void test_machine_bit_identity() {
    assert_eval_eq("s5v = N[Sin[Range[5000]]] 10^3; s5l = LinkedList[s5v]; "
                   "{Normal[Accumulate[s5l]] === Accumulate[s5v], Normal[Sort[s5l]] === Sort[s5v], "
                   "Normal[Differences[s5l]] === Differences[s5v], Normal[Ratios[s5l]] === Ratios[s5v], "
                   "Max[s5l] === Max[s5v], Min[s5l] === Min[s5v], Normal[Ordering[s5l]] === Ordering[s5v], "
                   "Normal[ReverseSort[s5l]] === ReverseSort[s5v]}",
                   "{True, True, True, True, True, True, True, True}", 0);
    assert_eval_eq("s5i = LinkedList[{2^62, 2^62, -3}]; {Normal[Accumulate[s5i]], Normal[Differences[s5i]]}",
                   "{{4611686018427387904, 9223372036854775808, 9223372036854775805}, "
                   "{0, -4611686018427387907}}", 0);
    assert_eval_eq("Normal[Accumulate[LinkedList[{1.*^308, 1.*^308}]]] === Accumulate[{1.*^308, 1.*^308}]",
                   "True", 0);
}

/* Callbacks may change the list (the elements come from a snapshot), and a
 * Throw leaves the head. */
void test_callbacks() {
    assert_eval_eq("s6x = LinkedList[{3, 1, 2}]; {Normal[SortBy[s6x, (AppendTo[s6x, 0]; -#) &]], Length[s6x]}",
                   "{{3, 2, 1}, 6}", 0);
    assert_eval_eq("s6y = LinkedList[{1, 2, 3, 4}]; Normal /@ Normal[GatherBy[s6y, (PopFirst[s6y]; OddQ[#]) &]]",
                   "{{1, 3}, {2, 4}}", 0);
    assert_eval_eq("Catch[GroupBy[LinkedList[{1, 2, 3}], If[# == 2, Throw[s6t], #] &]]", "s6t", 0);
}

/* Replacement with a LinkedList of rules (or of rule sets). */
void test_replace_with_linked_rules() {
    assert_eval_eq("s7r = LinkedList[{s7a -> 1, s7b :> 2 + 3}]; c0 = Internal`LinkedListConversionCount[]; "
                   "{{s7a, s7b, s7c} /. s7r, Replace[s7a, s7r], {s7a, s7b} //. LinkedList[{s7a -> s7b, s7b -> s7c}], "
                   "ReplaceList[s7a, LinkedList[{s7a -> 1, s7a -> 2}]], "
                   "Internal`LinkedListConversionCount[] == c0}",
                   "{{1, 5, s7c}, 1, {s7c, s7c}, {1, 2}, True}", 0);
    assert_eval_eq("Normal[{s7a, s7b} /. LinkedList[{{s7a -> 1}, {s7b -> 2}}]]", "{{1, s7b}, {s7a, 2}}", 0);
    /* A LinkedList of rule sets reaches each set's fixed point (it used to grow
     * without bound: the per-scan answer was a new LinkedList every time). */
    assert_eval_eq("{Normal[{s7a, s7b} //. LinkedList[{{s7a -> 1}, {s7b -> 2}}]], "
                   "Normal[ReplaceRepeated[s7f, LinkedList[{3, 1}]]]}",
                   "{{{1, s7b}, {s7a, 2}}, {s7f, s7f}}", 0);
}

/* RandomChoice / RandomSample draw the same stream as from the List. */
void test_random() {
    assert_eval_eq("s8l = LinkedList[{1., 2., 3.}]; "
                   "{SeedRandom[3]; Normal[RandomChoice[s8l, 6]], SeedRandom[3]; RandomChoice[{1., 2., 3.}, 6], "
                   "SeedRandom[3]; Normal /@ Normal[RandomChoice[s8l, {2, 2}]], "
                   "SeedRandom[3]; RandomChoice[{1., 2., 3.}, {2, 2}], "
                   "SeedRandom[5]; Normal[RandomChoice[LinkedList[{1, 2}] -> LinkedList[{s8a, s8b}], 5]], "
                   "SeedRandom[5]; RandomChoice[{1, 2} -> {s8a, s8b}, 5], "
                   "SeedRandom[7]; Normal[RandomSample[LinkedList[Range[6]], 3]], "
                   "SeedRandom[7]; RandomSample[Range[6], 3], "
                   "SeedRandom[7]; RandomChoice[s8l], SeedRandom[7]; RandomChoice[{1., 2., 3.}]} // "
                   "Partition[#, 2] & // Map[SameQ @@ # &]",
                   "{True, True, True, True, True}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_sorting);
    TEST(test_sets);
    TEST(test_grouping_and_assoc);
    TEST(test_restructuring);
    TEST(test_machine_bit_identity);
    TEST(test_callbacks);
    TEST(test_replace_with_linked_rules);
    TEST(test_random);

    printf("All LinkedList set/grouping tests passed!\n");
    symtab_clear();
    return 0;
}
