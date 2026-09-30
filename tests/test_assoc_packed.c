/* Unit tests for the Association builtins on the packed-array and visible
 * NDArray surfaces (src/assoc_packed.c, the with_list_args routing in
 * src/assoc.c, and the Values/Keys packing).
 *
 * The contract under test, for every association head that takes a list:
 *   1. SAME ANSWER on all three representations -- a plain List, a packed List
 *      and a visible NDArray[...] -- because a representation may never change
 *      a value (src/pack.h). The plain twin is FromNDArray[packed], which is
 *      the materialised list, so each check compares the head against its own
 *      List path rather than a hand-written constant.
 *   2. A VISIBLE NDArray is never left unevaluated (CLAUDE.md: that is a wrong
 *      answer, not a slow one).
 *   3. The NO-NESTING INVARIANT holds on the way out: nothing packed ends up
 *      inside an association (a structural pattern over the result matches,
 *      which it could not if a buffer were nested in it).
 *   4. The exactness gate of the compiled key function: a key function the VM
 *      would answer differently (x/2 on integers, a tie in Round, overflow)
 *      still gives the interpreter's keys.
 *   5. Values / Keys / a bulk Lookup of machine numbers come back PACKED, and
 *      a mixed or short column does not.
 *   6. Repeated evaluation, for leaks and double frees (run under leaks/valgrind).
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

static void run(const char* input) {
    Expr* r = eval_and_free(parse_expression(input));
    expr_free(r);
}

/* Integer data with repeats, and its Real twin; both over the 250-element
 * packing threshold so the buffer paths are live. */
static void setup(void) {
    run("tv = Mod[Range[2000]*7919, 37] - 18;");
    run("tr = N[tv]/4;");
    run("ta = AssociationThread[Range[1000], N[Range[1000]]/2];");
}

static void test_inputs_are_packed(void) {
    assert_true("PackedArrayQ[tv] && PackedArrayQ[tr]");
    assert_true("!PackedArrayQ[FromNDArray[tv]]");
}

/* ---- 1. packed == plain -------------------------------------------------- */

static void test_positionindex(void) {
    assert_true("PositionIndex[tv] === PositionIndex[FromNDArray[tv]]");
    assert_true("PositionIndex[tr] === PositionIndex[FromNDArray[tr]]");
    assert_true("Keys[PositionIndex[tv]] === DeleteDuplicates[FromNDArray[tv]]");
    assert_true("Total[Length /@ Values[PositionIndex[tv]]] === 2000");
    /* 0. and -0. stay two keys, exactly as expr_eq (and Tally) keep them */
    assert_true("PositionIndex[ToNDArray[Join[{0., -0.}, Range[1., 300.]]]] === "
                "PositionIndex[Join[{0., -0.}, Range[1., 300.]]]");
}

static void test_groupby(void) {
    assert_true("GroupBy[tv, Mod[#, 5] &] === GroupBy[FromNDArray[tv], Mod[#, 5] &]");
    assert_true("GroupBy[tv, EvenQ] === GroupBy[FromNDArray[tv], EvenQ]");
    assert_true("GroupBy[tv, # > 0 &] === GroupBy[FromNDArray[tv], # > 0 &]");
    assert_true("GroupBy[tv, Abs[#] + 1 &, Total] === GroupBy[FromNDArray[tv], Abs[#] + 1 &, Total]");
    assert_true("GroupBy[tv, Mod[#, 3] &, Mean] === GroupBy[FromNDArray[tv], Mod[#, 3] &, Mean]");
    assert_true("GroupBy[tr, Floor, Length] === GroupBy[FromNDArray[tr], Floor, Length]");
    assert_true("GroupBy[tr, Floor[#/2.] &] === GroupBy[FromNDArray[tr], Floor[#/2.] &]");
    assert_true("GroupBy[tv, Sign] === GroupBy[FromNDArray[tv], Sign]");
    /* a key function outside the compiled subset takes the List path */
    assert_true("GroupBy[tv, If[# > 0, \"pos\", \"neg\"] &] === "
                "GroupBy[FromNDArray[tv], If[# > 0, \"pos\", \"neg\"] &]");
    /* the key -> value transform form */
    assert_true("GroupBy[tv, EvenQ -> (#^2 &)] === GroupBy[FromNDArray[tv], EvenQ -> (#^2 &)]");
}

static void test_groupby_exactness_gate(void) {
    /* x/2 on integers is a Rational to the interpreter and a double to the VM */
    assert_true("GroupBy[tv, #/2 &] === GroupBy[FromNDArray[tv], #/2 &]");
    assert_true("Head[First[Keys[GroupBy[tv, #/2 &]]]] === Rational || "
                "IntegerQ[First[Keys[GroupBy[tv, #/2 &]]]]");
    /* Round goes to the even neighbour on a tie */
    assert_true("GroupBy[tr, Round] === GroupBy[FromNDArray[tr], Round]");
    /* integer overflow leaves the machine range: the VM bails, bignum keys */
    assert_true("GroupBy[tv, # 10^18 &] === GroupBy[FromNDArray[tv], # 10^18 &]");
    /* Equal on Reals is tolerant in the interpreter; never compiled */
    assert_true("GroupBy[tr, # == 1. &] === GroupBy[FromNDArray[tr], # == 1. &]");
}

static void test_gatherby_countsby(void) {
    assert_true("GatherBy[tv, Mod[#, 4] &] === GatherBy[FromNDArray[tv], Mod[#, 4] &]");
    assert_true("CountsBy[tv, Mod[#, 5] &] === CountsBy[FromNDArray[tv], Mod[#, 5] &]");
    assert_true("CountsBy[tr, Positive] === CountsBy[FromNDArray[tr], Positive]");
}

static void test_thread_and_map(void) {
    assert_true("AssociationThread[tv, Range[2000]] === AssociationThread[FromNDArray[tv], Range[2000]]");
    assert_true("AssociationThread[tv -> tr] === AssociationThread[FromNDArray[tv] -> FromNDArray[tr]]");
    assert_true("AssociationThread[tr, tv] === AssociationThread[FromNDArray[tr], FromNDArray[tv]]");
    assert_true("AssociationMap[#^2 &, tv] === AssociationMap[#^2 &, FromNDArray[tv]]");
    assert_true("AssociationMap[Sqrt, Range[300]] === AssociationMap[Sqrt, FromNDArray[Range[300]]]");
    assert_true("AssociationMap[f, tv] === AssociationMap[f, FromNDArray[tv]]");
}

static void test_lookup_buffer_keys(void) {
    run("tk = Mod[Range[600]*31, 1200] + 1;");
    assert_true("PackedArrayQ[tk]");
    assert_true("Lookup[ta, tk] === Lookup[ta, FromNDArray[tk]]");
    assert_true("Lookup[ta, tk, 0.] === Lookup[ta, FromNDArray[tk], 0.]");
    /* a packed DEFAULT is copied into the answer, so it is materialised */
    assert_true("Lookup[ta, {5000}, Range[300]] === {Range[300]}");
    assert_true("Lookup[ta, N[tk]] === Lookup[ta, FromNDArray[N[tk]]]");
}

/* ---- 2. visible NDArray -------------------------------------------------- */

static void test_visible_ndarray(void) {
    assert_eval_eq("AssociationThread[{a, b, c}, NDArray[{1., 2., 3.}]]",
                   "<|a -> 1.0, b -> 2.0, c -> 3.0|>", 0);
    assert_eval_eq("AssociationThread[NDArray[{1, 2, 1}, DataType -> \"int64\"], {x, y, z}]",
                   "<|1 -> z, 2 -> y|>", 0);
    assert_eval_eq("AssociationMap[#^2 &, NDArray[{1., 2., 3.}]]",
                   "<|1.0 -> 1.0, 2.0 -> 4.0, 3.0 -> 9.0|>", 0);
    assert_eval_eq("GroupBy[NDArray[{1., 2., 3.}], # > 1 &]",
                   "<|False -> {1.0}, True -> {2.0, 3.0}|>", 0);
    assert_eval_eq("GroupBy[NDArray[{1., 2., 3.}], # > 1 &, Total]",
                   "<|False -> 1.0, True -> 5.0|>", 0);
    assert_eval_eq("GroupBy[NDArray[{{1, 2}, {3, 4}, {1, 5}}, DataType -> \"int64\"], First]",
                   "<|1 -> {{1, 2}, {1, 5}}, 3 -> {{3, 4}}|>", 0);
    assert_eval_eq("GatherBy[NDArray[{1, 2, 3, 4}, DataType -> \"int64\"], EvenQ]",
                   "{{1, 3}, {2, 4}}", 0);
    assert_eval_eq("CountsBy[NDArray[{1.5, 2.5, 3.5}], Floor]", "<|1 -> 1, 2 -> 1, 3 -> 1|>", 0);
    assert_eval_eq("PositionIndex[NDArray[{1, 2, 1}, DataType -> \"int64\"]]",
                   "<|1 -> {1, 3}, 2 -> {2}|>", 0);
    /* complex: no word key, so the materialised List path answers */
    assert_true("PositionIndex[NDArray[{1 + I, 2, 1 + I}, DataType -> \"complex64\"]] === "
                "PositionIndex[{1. + 1. I, 2., 1. + 1. I}]");
    assert_true("Counts[NDArray[{1 + I, 2, 1 + I}, DataType -> \"complex64\"]] === "
                "Counts[{1. + 1. I, 2., 1. + 1. I}]");
    assert_eval_eq("Lookup[<|1 -> x, 2 -> y, 3 -> z|>, NDArray[{1, 3}, DataType -> \"int64\"]]",
                   "{x, z}", 0);
    assert_eval_eq("KeyDrop[<|1 -> x, 2 -> y, 3 -> z|>, NDArray[{1, 3}, DataType -> \"int64\"]]",
                   "<|2 -> y|>", 0);
    assert_eval_eq("KeyTake[<|1 -> x, 2 -> y, 3 -> z|>, NDArray[{1, 3}, DataType -> \"int64\"]]",
                   "<|1 -> x, 3 -> z|>", 0);
}

/* ---- 3. nothing packed nested in a result -------------------------------- */

static void test_no_nesting(void) {
    /* A structural pattern walks args[]: a nested buffer would not match. */
    assert_true("MatchQ[PositionIndex[tv], <|(_Integer -> {__Integer}) ..|>]");
    assert_true("MatchQ[GroupBy[tv, Mod[#, 3] &], <|(_Integer -> {__Integer}) ..|>]");
    assert_true("MatchQ[GroupBy[tv, Mod[#, 3] &, Identity], <|(_Integer -> {__Integer}) ..|>]");
    assert_true("MatchQ[GroupBy[tv, Mod[#, 3] &, Sort], <|(_Integer -> {__Integer}) ..|>]");
    /* a packed list stored as a VALUE is boxed at Rule (the gate), by design */
    run("tb = <|\"x\" -> tr|>;");
    assert_true("MatchQ[tb, <|\"x\" -> {__Real}|>]");
    assert_true("!PackedArrayQ[tb[\"x\"]] && tb[\"x\"] === tr");
}

/* ---- 5. Values / Keys / Lookup packing ----------------------------------- */

static void test_values_keys_packed(void) {
    assert_true("PackedArrayQ[Values[ta]] && PackedArrayQ[Keys[ta]]");
    assert_true("Values[ta] === FromNDArray[Values[ta]] && Head[Values[ta]] === List");
    assert_true("Total[Values[ta]] === 250250.");
    assert_true("Keys[ta] === Range[1000]");
    /* integer values pack to int64 and stay exact */
    run("tc = AssociationThread[Range[400], 2^62 + Range[400]];");
    assert_true("PackedArrayQ[Values[tc]] && Last[Values[tc]] === 2^62 + 400");
    /* Booleans pack to a bool buffer and read back as True/False */
    assert_true("PackedArrayQ[Values[AssociationThread[Range[300], EvenQ /@ Range[300]]]]");
    assert_true("Count[Values[AssociationThread[Range[300], EvenQ /@ Range[300]]], True] === 150");
    /* a mixed column does not pack; nor does a short one */
    assert_true("!PackedArrayQ[Values[Append[ta, 5000 -> 1]]]");
    assert_true("!PackedArrayQ[Values[AssociationThread[Range[10], N[Range[10]]]]]");
    assert_true("!PackedArrayQ[Keys[AssociationThread[ToString /@ Range[300], Range[300]]]]");
    /* the rule-list form packs the same way */
    assert_true("PackedArrayQ[Values[Normal[ta]]]");
    /* bulk Lookup of machine values */
    assert_true("PackedArrayQ[Lookup[ta, Range[500]]]");
    assert_true("!PackedArrayQ[Lookup[ta, Range[995, 1005]]]");
}

/* ---- 6. repeated evaluation ---------------------------------------------- */

static void test_repeated_no_leak(void) {
    for (int i = 0; i < 20; i++) {
        run("PositionIndex[tv]; GroupBy[tv, Mod[#, 5] &]; GroupBy[tv, Mod[#, 5] &, Total];"
            "GatherBy[tr, Floor]; CountsBy[tv, EvenQ]; AssociationThread[tv, tr];"
            "AssociationMap[# + 1 &, tv]; Lookup[ta, Range[400]]; Values[ta]; Keys[ta];"
            "GroupBy[NDArray[{1., 2., 3.}], # > 1 &]; GroupBy[tv, #/2 &];");
    }
    assert_true("Length[PositionIndex[tv]] === 37");
}

int main(void) {
    symtab_init();
    core_init();
    setup();

    TEST(test_inputs_are_packed);
    TEST(test_positionindex);
    TEST(test_groupby);
    TEST(test_groupby_exactness_gate);
    TEST(test_gatherby_countsby);
    TEST(test_thread_and_map);
    TEST(test_lookup_buffer_keys);
    TEST(test_visible_ndarray);
    TEST(test_no_nesting);
    TEST(test_values_keys_packed);
    TEST(test_repeated_no_leak);

    printf("All Association packed/NDArray tests passed!\n");
    return 0;
}
