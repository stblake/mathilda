/* test_sparsearray.c -- SparseArray (src/sparse/): the EXPR_SPARSE raw object,
 * its constructors, canonical form, structural queries, Part, arithmetic and
 * element-wise functions.
 *
 * Element-wise operations are checked against the invariant they promise:
 * f applied to SparseArrays equals f applied to their Normal[...] forms.
 *
 * Failures are counted rather than asserted: the CMake Release build passes
 * -DNDEBUG, which would no-op assert(). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

static int failures = 0;

/* Evaluate `input`; fail unless it prints exactly `expected`. */
static void chk(const char* input, const char* expected) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, s);
        failures++;
    }
    free(s);
    expr_free(e);
}

/* The InputForm (canonical form) of a SparseArray-valued input. */
static void chk_if(const char* input, const char* expected) {
    char buf[4096];
    snprintf(buf, sizeof buf, "InputForm[%s]", input);
    chk(buf, expected);
}

static void chk_true(const char* input) { chk(input, "True"); }

/* ------------------------------------------------------------------------ */

static void test_constructors(void) {
    /* explicit positions; the first rule for a position wins */
    chk_if("SparseArray[{{1,1}->1,{2,3}->x,{1,1}->5},{3,3}]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 2}, {{1}, {3}}}, {1, x}}]");
    chk_if("SparseArray[{1 -> a, 1 -> b}]", "SparseArray[Automatic, {1}, 0, {1, {{0, 1}, {{1}}}, {a}}]");
    chk_if("SparseArray[{{1}, {1}} -> {1, 2}]", "SparseArray[Automatic, {1}, 0, {1, {{0, 1}, {{1}}}, {1}}]");
    /* dims inferred from |index|; negative indices count from the end */
    chk_if("SparseArray[{-1 -> 9}, 4]", "SparseArray[Automatic, {4}, 0, {1, {{0, 1}, {{4}}}, {9}}]");
    chk_if("SparseArray[{-1 -> 9}]", "SparseArray[Automatic, {1}, 0, {1, {{0, 1}, {{1}}}, {9}}]");
    chk_if("SparseArray[{{1, 1} -> 1}, 3]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 1, 1}, {{1}}}, {1}}]");
    chk_if("SparseArray[{1 -> a}, Automatic]", "SparseArray[Automatic, {1}, 0, {1, {{0, 1}, {{1}}}, {a}}]");
    /* background-valued entries are left implicit at construction */
    chk_if("SparseArray[{1 -> 0, 2 -> 1}]", "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{2}}}, {1}}]");
    chk_if("SparseArray[{1 -> 0., 2 -> 1}]", "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{2}}}, {1}}]");
    chk_if("SparseArray[{1 -> 0.}, 2, 0]", "SparseArray[Automatic, {2}, 0, {1, {{0, 0}, {}}, {}}]");
    /* {pos..} -> {vals..} and {pos..} -> v */
    chk_if("SparseArray[{{1}, {2}} -> {a, b}, 3]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 2}, {{1}, {2}}}, {a, b}}]");
    chk_if("SparseArray[{{1, 1}, {2, 2}} -> a, {2, 2}]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 2}, {{1}, {2}}}, {a, a}}]");
    /* patterns */
    chk_if("SparseArray[{i_, j_} /; i > j :> i + j, {3, 3}]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 0, 1, 3}, {{1}, {1}, {2}}}, {3, 4, 5}}]");
    chk_if("SparseArray[{{i_, i_} -> i}, {3, 3}]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {2}, {3}}}, {1, 2, 3}}]");
    chk_if("SparseArray[{x_ /; x > 2 -> 1}, 4]",
           "SparseArray[Automatic, {4}, 0, {1, {{0, 2}, {{3}, {4}}}, {1, 1}}]");
    chk("Normal[SparseArray[{i_, j_} :> i + 10 j, {2, 3}]]", "{{11, 21, 31}, {12, 22, 32}}");
    chk("Normal[SparseArray[{i_} :> i^2, 4]]", "{1, 4, 9, 16}");
    /* an all-matching rule names the background (and shadows later rules) */
    chk_if("SparseArray[{x_ -> 1}, 3]", "SparseArray[Automatic, {3}, 1, {1, {{0, 0}, {}}, {}}]");
    chk_if("SparseArray[{{1, 1} -> 1, {_, _} -> 5}]",
           "SparseArray[Automatic, {1, 1}, 5, {1, {{0, 1}, {{1}}}, {1}}]");
    chk_if("SparseArray[{{1, 1} -> 1, {_, _} -> 5}, {2, 2}]",
           "SparseArray[Automatic, {2, 2}, 5, {1, {{0, 1, 1}, {{1}}}, {1}}]");
    chk_if("SparseArray[{{1, 1} -> 1, {_, _} -> 5}, {2, 2}, 7]",
           "SparseArray[Automatic, {2, 2}, 7, {1, {{0, 2, 4}, {{1}, {2}, {1}, {2}}}, {1, 5, 5, 5}}]");
    /* Band: scalar, vector, subblocks, steps, inferred dims */
    chk_if("SparseArray[Band[{1, 1}] -> {1, 2, 3}]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {2}, {3}}}, {1, 2, 3}}]");
    chk_if("SparseArray[Band[{1, 1}] -> 1]", "SparseArray[Automatic, {1, 1}, 0, {1, {{0, 1}, {{1}}}, {1}}]");
    chk("Normal[SparseArray[Band[{1, 1}] -> {{a, b}, {c, d}}, {4, 4}]]",
        "{{a, b, 0, 0}, {c, d, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}");
    chk("Normal[SparseArray[Band[{1, 1}] -> {{{a, b}, {c, d}}, {{e, f}, {g, h}}}, {4, 4}]]",
        "{{a, b, 0, 0}, {c, d, 0, 0}, {0, 0, e, f}, {0, 0, g, h}}");
    chk("Normal[SparseArray[Band[{2, 1}] -> x, {3, 3}]]", "{{0, 0, 0}, {x, 0, 0}, {0, x, 0}}");
    chk("Normal[SparseArray[Band[{1, 1}, {2, 2}] -> x, {3, 3}]]", "{{x, 0, 0}, {0, x, 0}, {0, 0, 0}}");
    chk("Normal[SparseArray[Band[{-1, -1}] -> x, {3, 3}]]", "{{0, 0, 0}, {0, 0, 0}, {0, 0, x}}");
    chk("Normal[SparseArray[Band[{1, 3}, {3, 1}, {1, -1}] -> x, {3, 3}]]",
        "{{0, 0, x}, {0, x, 0}, {x, 0, 0}}");
    chk("Normal[SparseArray[{Band[{1, 1}] -> 1, Band[{1, 2}] -> 2}, {3, 3}]]",
        "{{1, 2, 0}, {0, 1, 2}, {0, 0, 1}}");
    chk("Normal[SparseArray[{Band[{1, 1}] -> {a, b}, {1, 1} -> z}, {2, 2}]]", "{{a, 0}, {0, b}}");
    chk("Normal[SparseArray[{Band[{1, 1}] -> {1, 2}, {3, 3} -> 5}]]",
        "{{1, 0, 0}, {0, 2, 0}, {0, 0, 5}}");
    /* dense arrays: the zero is 0. for an all-inexact array */
    chk_if("SparseArray[{1., 0., 2.}]",
           "SparseArray[Automatic, {3}, 0.0, {1, {{0, 2}, {{1}, {3}}}, {1.0, 2.0}}]");
    chk_if("SparseArray[{1., 0, 2.}]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 2}, {{1}, {3}}}, {1.0, 2.0}}]");
    chk_if("SparseArray[{{1, 0}, {0, 2}}]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 2}, {{1}, {2}}}, {1, 2}}]");
    chk_if("SparseArray[{1. + I, 0.}]",
           "SparseArray[Automatic, {2}, 0.0, {1, {{0, 1}, {{1}}}, {1.0 + 1.0*I}}]");
    chk_if("SparseArray[{a, 0, b}, 3, 7]",
           "SparseArray[Automatic, {3}, 7, {1, {{0, 3}, {{1}, {2}, {3}}}, {a, 0, b}}]");
    chk_if("SparseArray[{1, 0, 2}, 3, 0.]",
           "SparseArray[Automatic, {3}, 0.0, {1, {{0, 2}, {{1}, {3}}}, {1.0, 2.0}}]");
    chk_if("SparseArray[{{1, 2}, {3, 4}}, {3, 3}]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 2, 4, 4}, {{1}, {2}, {1}, {2}}}, {1, 2, 3, 4}}]");
    chk_if("SparseArray[{{1, 2}, {3, 4}}, {1, 1}]",
           "SparseArray[Automatic, {1, 1}, 0, {1, {{0, 1}, {{1}}}, {1}}]");
    chk_if("SparseArray[{True, False}]",
           "SparseArray[Automatic, {2}, 0, {1, {{0, 2}, {{1}, {2}}}, {True, False}}]");
    chk_if("SparseArray[{\"a\", 0}]", "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{1}}}, {\"a\"}}]");
    chk_if("SparseArray[{}, 3]", "SparseArray[Automatic, {3}, 0, {1, {{0, 0}, {}}, {}}]");
    /* packed input: read straight from the buffer */
    chk_if("SparseArray[Range[5.] - 3.]",
           "SparseArray[Automatic, {5}, 0.0, {1, {{0, 4}, {{1}, {2}, {4}, {5}}}, {-2.0, -1.0, 1.0, 2.0}}]");
    chk_if("SparseArray[Range[6] - 3]",
           "SparseArray[Automatic, {6}, 0, {1, {{0, 5}, {{1}, {2}, {4}, {5}, {6}}}, {-2, -1, 1, 2, 3}}]");
    chk_if("SparseArray[Transpose[{Range[4], Range[4]}] -> N[Range[4]], {4, 4}]",
           "SparseArray[Automatic, {4, 4}, 0, {1, {{0, 1, 2, 3, 4}, {{1}, {2}, {3}, {4}}}, {1.0, 2.0, 3.0, 4.0}}]");
    chk_if("SparseArray[NDArray[{0., 3.}]]",
           "SparseArray[Automatic, {2}, 0.0, {1, {{0, 1}, {{2}}}, {3.0}}]");
    /* existing arrays: minimise, resize, rebackground */
    chk_if("SparseArray[SparseArray[{1, 1, 1, 0}, 4]]",
           "SparseArray[Automatic, {4}, 0, {1, {{0, 3}, {{1}, {2}, {3}}}, {1, 1, 1}}]");
    chk_if("SparseArray[SparseArray[{1 -> 1, 2 -> -1}, 2] + SparseArray[{1 -> -1}, 2]]",
           "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{2}}}, {-1}}]");
    chk_if("SparseArray[SparseArray[{{1,1}->1,{2,3}->x},{3,3}], {4, 4}]",
           "SparseArray[Automatic, {4, 4}, 0, {1, {{0, 1, 2, 2, 2}, {{1}, {3}}}, {1, x}}]");
    chk_if("SparseArray[SparseArray[{{1,1}->1,{2,3}->x},{3,3}], {2, 2}]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 1}, {{1}}}, {1}}]");
    chk_if("SparseArray[SparseArray[{{1,1}->1,{2,3}->x},{3,3}], {3, 3}, 5]",
           "SparseArray[Automatic, {3, 3}, 5, {1, {{0, 3, 6, 9}, {{1}, {2}, {3}, {1}, {2}, {3}, "
           "{1}, {2}, {3}}}, {1, 0, 0, 0, 0, x, 0, 0, 0}}]");
    /* the canonical form reads back; explicit entries are kept as given */
    chk("Normal[SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {2}, {3}}}, {1, 1, 1}}]]",
        "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}");
    chk_if("SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{1}}}, {0}}]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{1}}}, {0}}]");
    chk_if("SparseArray[Automatic, {3}, 0, {1, {{0, 2}, {{1}, {1}}}, {a, b}}]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{1}}}, {a}}]");
    chk_true("With[{s = SparseArray[{{1,1}->1,{2,3}->x},{3,3}]}, ToExpression[ToString[s, InputForm]] === s]");
    chk_true("With[{s = SparseArray[{{1,2,3}->a,{2,1,1}->2.5},{2,2,3}]}, ToExpression[ToString[s, InputForm]] === s]");
    /* rank 3 */
    chk_if("SparseArray[{{1, 2, 3} -> a, {2, 1, 1} -> b}, {2, 2, 3}]",
           "SparseArray[Automatic, {2, 2, 3}, 0, {1, {{0, 1, 2}, {{2, 3}, {1, 1}}}, {a, b}}]");
    /* bigint values are boxed */
    chk_if("SparseArray[{1 -> 2^70}, 2]",
           "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{1}}}, {1180591620717411303424}}]");
}

static void test_messages(void) {
    chk("Quiet[Check[SparseArray[{5 -> 9}, 4], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{0 -> 9}, 4], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{a -> 1}, 2], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{{1, 2} -> 1, 3 -> 2}], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{1 -> 1}, {2, -1}], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{{1, 1}, {2, 2}} -> {a}, {2, 2}], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[x], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{{1, 2}, {3}}], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{{}}], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{}], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{1 -> {1, 2}}, 2], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{1 -> 1}, 3, {1}], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{{1, 2}, {3, 4}}, {2}], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{{i_, i_} -> 1}], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{1 -> SparseArray[{1 -> 1}]}, 2], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{1 -> 1}, 2] + SparseArray[{1 -> 1}, 3], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{1 -> 1}, 2][\"foo\"], $Failed]]", "$Failed");
    chk("Quiet[Check[SparseArray[{1 -> 1}, 2][[3]], $Failed]]", "$Failed");
    /* a malformed spec stays unevaluated */
    chk("Quiet[SparseArray[x]]", "SparseArray[x]");
    chk("Quiet[SparseArray[{5 -> 9}, 4]]", "SparseArray[{5 -> 9}, 4]");
    /* Band with too many values warns but still builds */
    chk("Quiet[Normal[SparseArray[Band[{1, 1}] -> {1, 2, 3, 4, 5}, {3, 3}]]]",
        "{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}");
}

static void test_raw_object(void) {
    chk("s = SparseArray[{{1,1}->1,{2,3}->x},{3,3}]", "SparseArray[<2>, {3, 3}]");
    chk("{s, 1}", "{SparseArray[<2>, {3, 3}], 1}");
    chk("SparseArray[{}, {2, 2}]", "SparseArray[<0>, {2, 2}]");
    chk("FullForm[SparseArray[{1 -> a}, 2]]",
        "SparseArray[Automatic, List[2], 0, List[1, List[List[0, 1], List[List[1]]], List[a]]]");
    chk("AtomQ[s]", "True");
    chk("Head[s]", "SparseArray");
    chk("MatchQ[s, _SparseArray]", "True");
    chk("Quiet[MatchQ[s, SparseArray[___]]]", "False");
    chk("s /. SparseArray[a_, b_, c_, d_] :> {b}", "SparseArray[<2>, {3, 3}]");
    chk("Cases[{1, s, x}, _SparseArray]", "{SparseArray[<2>, {3, 3}]}");
    chk("NumericQ[s]", "False");
    /* values are never re-evaluated */
    chk("x = 7; r = InputForm[s]; Clear[x]; r",
        "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 2}, {{1}, {3}}}, {1, x}}]");
    /* SameQ is equality of the canonical form */
    chk("SparseArray[{1 -> 2}, 2] === SparseArray[{1 -> 2, 2 -> 0}, 2]", "True");
    chk("SparseArray[{1 -> 2}, 2] === SparseArray[{1 -> 2.}, 2]", "False");
    chk("SparseArray[{1 -> 1}, 2] === {1, 0}", "False");
    chk("SparseArray[{1 -> 1}, 2] == SparseArray[{1 -> 1}, 2]", "True");
    chk("Length[Union[{SparseArray[{1 -> 2}, 2], SparseArray[{1 -> 2, 2 -> 0}, 2], SparseArray[{2 -> 2}, 2]}]]", "2");
    chk("OrderedQ[{SparseArray[{1 -> 1}, 2], SparseArray[{1 -> 2}, 2]}]", "True");
    chk("Sort[{SparseArray[{1 -> 2}, 2], SparseArray[{1 -> 1}, 3], 1, x}]",
        "{1, SparseArray[<1>, {2}], SparseArray[<1>, {3}], x}");
    chk("TeXForm[SparseArray[{1 -> 1}, 2]]", "\\{1,0\\}");
}

static void test_structure(void) {
    chk("s = SparseArray[{{1,1}->1,{2,3}->x},{3,3}]; t = SparseArray[{{1,2,3}->a,{2,1,1}->b},{2,2,3}];", "Null");
    chk("Normal[s]", "{{1, 0, 0}, {0, 0, x}, {0, 0, 0}}");
    chk("Normal[{SparseArray[{2 -> 1}], 3}]", "{{0, 1}, 3}");
    chk("Normal[SparseArray[{2 -> 1.5}, 3, 0.]]", "{0.0, 1.5, 0.0}");
    chk("Normal[SparseArray[{1 -> 1.}, 3]]", "{1.0, 0, 0}");
    chk("Developer`PackedArrayQ[Normal[SparseArray[{1 -> 1.}, 30, 0.]]]", "True");
    chk("Normal[SparseArray[{}, {0, 3}]]", "{}");
    chk("Normal[SparseArray[{}, {2, 0}]]", "{{}, {}}");
    chk("Length[s]", "3");
    chk("Length[t]", "2");
    chk("Dimensions[t]", "{2, 2, 3}");
    chk("Dimensions[t, 2]", "{2, 2}");
    chk("ArrayDepth[t]", "3");
    chk("MatrixQ[s]", "True");
    chk("MatrixQ[t]", "False");
    chk("VectorQ[SparseArray[{1 -> 1}, 2]]", "True");
    chk("VectorQ[SparseArray[{1 -> 1}, 2], NumberQ]", "True");
    chk("VectorQ[SparseArray[{1 -> x}, 2], NumberQ]", "False");
    chk("MatrixQ[s, NumberQ]", "False");
    chk("ArrayQ[s]", "True");
    chk("ArrayQ[t, 3]", "True");
    chk("ArrayQ[t, 2]", "False");
    chk("ArrayRules[s]", "{{1, 1} -> 1, {2, 3} -> x, {_, _} -> 0}");
    chk("ArrayRules[t]", "{{1, 2, 3} -> a, {2, 1, 1} -> b, {_, _, _} -> 0}");
    chk("ArrayRules[{0, 1, 0, 2}]", "{{2} -> 1, {4} -> 2, {_} -> 0}");
    chk("ArrayRules[s, 1]",
        "{{1, 2} -> 0, {1, 3} -> 0, {2, 1} -> 0, {2, 2} -> 0, {2, 3} -> x, {3, 1} -> 0, {3, 2} -> 0, "
        "{3, 3} -> 0, {_, _} -> 1}");
    chk("SparseArray[ArrayRules[s]] === SparseArray[s, {3, 3}] || Normal[SparseArray[ArrayRules[s], {3, 3}]] === Normal[s]", "True");
    /* properties */
    chk("s[\"NonzeroPositions\"]", "{{1, 1}, {2, 3}}");
    chk("s[\"NonzeroValues\"]", "{1, x}");
    chk("s[\"Background\"]", "0");
    chk("s[\"ColumnIndices\"]", "{{1}, {3}}");
    chk("s[\"RowPointers\"]", "{0, 1, 2, 2}");
    chk("s[\"ExplicitLength\"]", "2");
    chk("s[\"Density\"] == 2/9", "True");
    chk("SparseArray[{{1,1}->1,{3,2}->2},{3,3}][\"BandWidth\"]", "{1, 0}");
    chk("SparseArray[{{1,1}->1,{3,2}->2},{3,3}][\"AdjacencyLists\"]", "{{1}, {}, {2}}");
    /* ArrayQ / ArrayDepth on ordinary arrays */
    chk("ArrayQ[{{1, 2}, {3, 4}}]", "True");
    chk("ArrayQ[{1, {2}}]", "False");
    chk("ArrayQ[{1, 2}, 1]", "True");
    chk("ArrayQ[{1, 2}, _, IntegerQ]", "True");
    chk("ArrayQ[{1, 2.}, _, IntegerQ]", "False");
    chk("ArrayQ[x]", "False");
    chk("ArrayQ[Range[10.]]", "True");
    chk("ArrayDepth[{{1, 2}, {3, 4}}]", "2");
    chk("ArrayDepth[{{1, 2}, {3}}]", "1");
    chk("ArrayDepth[x]", "0");
}

static void test_part(void) {
    chk("s = SparseArray[{{1,1}->1,{2,3}->x},{3,3}]; t = SparseArray[{{1,2,3}->a,{2,1,1}->b},{2,2,3}];", "Null");
    chk_if("s[[2]]", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{3}}}, {x}}]");
    chk("s[[2, 3]]", "x");
    chk("s[[2, -1]]", "x");
    chk("s[[3, 3]]", "0");
    chk_if("s[[All, 3]]", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{2}}}, {x}}]");
    chk_if("s[[1 ;; 2]]", "SparseArray[Automatic, {2, 3}, 0, {1, {{0, 1, 2}, {{1}, {3}}}, {1, x}}]");
    chk_if("s[[{2, 1}, 3]]", "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{1}}}, {x}}]");
    chk_if("SparseArray[{1, 2, 3}][[{1, 1, 3}]]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 3}, {{1}, {2}, {3}}}, {1, 1, 3}}]");
    chk("SparseArray[{{1, 0}, {0, 2}}][[2 ;; 1]]", "{}");
    chk("s[[0]]", "List");
    chk("s[[1, 0]]", "List");
    chk_if("t[[2]]", "SparseArray[Automatic, {2, 3}, 0, {1, {{0, 1, 1}, {{1}}}, {b}}]");
    chk("t[[1, 2, 3]]", "a");
    chk("Quiet[s[[4]]]", "Part[SparseArray[<2>, {3, 3}], 4]");
    chk("Normal[s[[{3, 2}, {3, 1}]]] === Normal[s][[{3, 2}, {3, 1}]]", "True");
    /* a SparseArray reached below the top level: the rest of the spec is its
     * own Part ({c0, sa}[[2, 1]], as CoefficientArrays results are indexed) */
    chk("{1, s}[[2, 2, 3]]", "x");
    chk_if("{1, s}[[2, 2]]", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{3}}}, {x}}]");
    chk_if("{{s}}[[1, 1, All, 3]]", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{2}}}, {x}}]");
    chk("Extract[{1, s}, {2, 1, 1}]", "1");
    chk("{1, s}[[2, 0]]", "List");
    chk("Head[{1, s}[[2]]]", "SparseArray");
    /* element assignment */
    chk_if("u = SparseArray[{1 -> 1}, 3]; u[[2]] = 5; u",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 2}, {{1}, {2}}}, {1, 5}}]");
    chk_if("u[[1]] = 0; u", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{2}}}, {5}}]");
    chk_if("v = s; v[[3, 3]] = q; v",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {3}, {3}}}, {1, x, q}}]");
    chk_if("s", "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 2}, {{1}, {3}}}, {1, x}}]");
    chk_if("w = SparseArray[{1 -> 1.5}, 3]; w[[3]] = 2.5; w",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 2}, {{1}, {3}}}, {1.5, 2.5}}]");
    /* a non-element assignment goes through the dense array */
    chk("m = SparseArray[{{1,1}->1},{2,2}]; m[[All, 2]] = 7; Normal[m]", "{{1, 7}, {0, 7}}");
    chk("Head[m]", "SparseArray");
}

static void test_arithmetic(void) {
    chk("s = SparseArray[{{1,1}->1,{2,3}->x},{3,3}];", "Null");
    chk_if("s + 1", "SparseArray[Automatic, {3, 3}, 1, {1, {{0, 1, 2, 2}, {{1}, {3}}}, {2, 1 + x}}]");
    chk("s + {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}", "{{2, 2, 3}, {4, 5, 6 + x}, {7, 8, 9}}");
    chk("Range[4] + SparseArray[{1 -> 1, 3 -> 2}, 4]", "{2, 2, 5, 4}");
    chk_if("SparseArray[{1 -> 1}, 2] + x",
           "SparseArray[Automatic, {2}, x, {1, {{0, 1}, {{1}}}, {1 + x}}]");
    chk_if("SparseArray[{1 -> 1}, 2] + 1 + x",
           "SparseArray[Automatic, {2}, 1 + x, {1, {{0, 1}, {{1}}}, {2 + x}}]");
    chk_if("SparseArray[{1 -> a, 2 -> b}, 3] * SparseArray[{2 -> c, 3 -> d}, 3]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{2}}}, {b*c}}]");
    chk_if("SparseArray[{1 -> a}, 3] * SparseArray[{2 -> c}, 3, 1]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 2}, {{1}, {2}}}, {a, 0}}]");
    chk_if("SparseArray[{1 -> a, 2 -> b}, 3] * SparseArray[{2 -> c, 3 -> d}, 3, 1]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 3}, {{1}, {2}, {3}}}, {a, b*c, 0}}]");
    chk_if("SparseArray[{1 -> 1, 2 -> -1}, 2] + SparseArray[{1 -> -1}, 2]",
           "SparseArray[Automatic, {2}, 0, {1, {{0, 2}, {{1}, {2}}}, {0, -1}}]");
    chk_if("SparseArray[{1 -> 1}, 2] + SparseArray[{2 -> 1}, 2, 1]",
           "SparseArray[Automatic, {2}, 1, {1, {{0, 1}, {{1}}}, {2}}]");
    chk_if("SparseArray[{1 -> 2}, 3] * SparseArray[{1 -> 3, 2 -> 4}, 3] * SparseArray[{2 -> 5}, 3]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 0}, {}}, {}}]");
    chk_if("SparseArray[{1 -> 1}, 3] + SparseArray[{2 -> 1}, 3] + SparseArray[{3 -> 1}, 3]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 3}, {{1}, {2}, {3}}}, {1, 1, 1}}]");
    chk_if("0 * SparseArray[{1 -> 1, 2 -> x}, 3]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 2}, {{1}, {2}}}, {0, 0}}]");
    chk_if("2 * SparseArray[{1 -> 1}, 3]", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{1}}}, {2}}]");
    chk_if("SparseArray[{1 -> 1}, 3] + 0.5",
           "SparseArray[Automatic, {3}, 0.5, {1, {{0, 1}, {{1}}}, {1.5}}]");
    chk_if("SparseArray[{1 -> 1, 3 -> 2}, 4] / 2",
           "SparseArray[Automatic, {4}, 0, {1, {{0, 2}, {{1}, {3}}}, {1/2, 1}}]");
    chk_if("-SparseArray[{1 -> 1, 3 -> 2}, 4]",
           "SparseArray[Automatic, {4}, 0, {1, {{0, 2}, {{1}, {3}}}, {-1, -2}}]");
    chk_if("2.5 SparseArray[{1 -> 1, 3 -> 2}, 4]",
           "SparseArray[Automatic, {4}, 0.0, {1, {{0, 2}, {{1}, {3}}}, {2.5, 5.0}}]");
    chk_if("SparseArray[{1 -> 2}, 3]^2", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{1}}}, {4}}]");
    chk_if("2^SparseArray[{1 -> 2}, 3]", "SparseArray[Automatic, {3}, 1, {1, {{0, 1}, {{1}}}, {4}}]");
    chk_if("Quiet[SparseArray[{1 -> a}, 3]^0]",
           "SparseArray[Automatic, {3}, Indeterminate, {1, {{0, 1}, {{1}}}, {1}}]");
    chk_if("Quiet[SparseArray[{1 -> 2}, 3]^SparseArray[{2 -> 2}, 3]]",
           "SparseArray[Automatic, {3}, Indeterminate, {1, {{0, 2}, {{1}, {2}}}, {1, 0}}]");
    chk_if("Quiet[SparseArray[{1 -> 1}, 2] / SparseArray[{2 -> 1}, 2]]",
           "SparseArray[Automatic, {2}, Indeterminate, {1, {{0, 2}, {{1}, {2}}}, {ComplexInfinity, 0}}]");
    /* typed fast paths: overflow goes exact, machine values stay machine */
    chk_if("SparseArray[{1 -> 2^62}, 2] * 4",
           "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{1}}}, {18446744073709551616}}]");
    chk_if("SparseArray[{1 -> 2^62, 2 -> 1}, 2] + SparseArray[{1 -> 2^62}, 2]",
           "SparseArray[Automatic, {2}, 0, {1, {{0, 2}, {{1}, {2}}}, {9223372036854775808, 1}}]");
    chk_if("SparseArray[{{1,1}->1.,{2,2}->2.,{3,1}->4.},{3,3}] + SparseArray[{{1,2}->1.,{2,2}->3.,{3,3}->5.},{3,3}]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 2, 3, 5}, {{1}, {2}, {2}, {1}, {3}}}, {1.0, 1.0, 5.0, 4.0, 5.0}}]");
    chk_if("SparseArray[{{1,1}->1.,{2,2}->2.,{3,1}->4.},{3,3}] * SparseArray[{{1,2}->1.,{2,2}->3.,{3,3}->5.},{3,3}]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 0, 1, 1}, {{2}}}, {6.0}}]");
    chk_if("SparseArray[{1 -> 1.5}, 3] + SparseArray[{1 -> -1.5}, 3, 2]",
           "SparseArray[Automatic, {3}, 2, {1, {{0, 1}, {{1}}}, {0.0}}]");
    chk_if("(1. + 2. I) SparseArray[{1 -> 1, 3 -> 2}, 4]",
           "SparseArray[Automatic, {4}, 0.0, {1, {{0, 2}, {{1}, {3}}}, {1.0 + 2.0*I, 2.0 + 4.0*I}}]");
    chk_if("SparseArray[{1 -> 1. + 2. I}, 2] + (-2. I)",
           "SparseArray[Automatic, {2}, 0.0 - 2.0*I, {1, {{0, 1}, {{1}}}, {1.0}}]");
    /* the invariant, over exact, symbolic and machine operands */
    chk("a1 = SparseArray[{1 -> 1, 3 -> 2}, 5]; b1 = SparseArray[{2 -> 3/2, 3 -> -2, 5 -> y}, 5, 1];", "Null");
    chk_true("Normal[a1 + b1] === Normal[a1] + Normal[b1]");
    chk_true("Normal[a1 b1] === Normal[a1] Normal[b1]");
    chk_true("Normal[a1 - 2 b1 + 3] === Normal[a1] - 2 Normal[b1] + 3");
    chk_true("Normal[a1^2 + b1^3] === Normal[a1]^2 + Normal[b1]^3");
    chk("c1 = SparseArray[{{1,1}->1.5,{2,3}->-2.,{4,4}->3.}, {4,4}]; "
        "d1 = SparseArray[{{1,1}->2.,{3,3}->1.,{4,1}->-1.}, {4,4}];", "Null");
    chk_true("Normal[c1 + d1] === Normal[c1] + Normal[d1]");
    chk_true("Normal[c1 - 0.5 d1] === Normal[c1] - 0.5 Normal[d1]");
    chk_true("Normal[c1 + 1.] === Normal[c1] + 1.");
    chk_true("Normal[c1 d1] == Normal[c1] Normal[d1]");
    /* pattern sharing does not alias: the operand is unchanged */
    chk_true("c2 = c1; c3 = 2 c1; c2 === c1 && Normal[c3] === 2 Normal[c1]");
}

static void test_functions(void) {
    chk_if("Exp[SparseArray[{1 -> 2.}, 3]]",
           "SparseArray[Automatic, {3}, 1, {1, {{0, 1}, {{1}}}, {7.38905609893065}}]");
    chk_if("Sin[SparseArray[{1 -> 1.}, 2]]",
           "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{1}}}, {0.8414709848078965}}]");
    chk_if("Cos[SparseArray[{1 -> 1.}, 2]]",
           "SparseArray[Automatic, {2}, 1, {1, {{0, 1}, {{1}}}, {0.5403023058681398}}]");
    chk_if("Sin[SparseArray[{1 -> Pi}, 3]]", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{1}}}, {0}}]");
    chk_if("Sin[SparseArray[{1 -> Pi}, 3, 1]]",
           "SparseArray[Automatic, {3}, Sin[1], {1, {{0, 1}, {{1}}}, {0}}]");
    chk_if("Exp[SparseArray[{1 -> x}, 2]]", "SparseArray[Automatic, {2}, 1, {1, {{0, 1}, {{1}}}, {E^x}}]");
    chk_if("Log[SparseArray[{1 -> 2}, 2]]",
           "SparseArray[Automatic, {2}, -Infinity, {1, {{0, 1}, {{1}}}, {Log[2]}}]");
    chk_if("Abs[SparseArray[{1 -> -2}, 2]]", "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{1}}}, {2}}]");
    chk_if("Re[SparseArray[{1 -> 1 + 2 I}, 2]]", "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{1}}}, {1}}]");
    chk_if("Floor[SparseArray[{1 -> 2.5}, 2]]", "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{1}}}, {2}}]");
    chk_if("Mod[SparseArray[{1 -> 7}, 2], 3]", "SparseArray[Automatic, {2}, 0, {1, {{0, 1}, {{1}}}, {1}}]");
    chk_if("N[SparseArray[{1 -> 1/2}, 2]]",
           "SparseArray[Automatic, {2}, 0.0, {1, {{0, 1}, {{1}}}, {0.5}}]");
    chk_if("N[SparseArray[{1 -> x}, 2]]", "SparseArray[Automatic, {2}, 0.0, {1, {{0, 1}, {{1}}}, {x}}]");
    chk_if("f[SparseArray[{1 -> a}, 3], SparseArray[{2 -> c}, 3]] /. f -> g",
           "g[SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{1}}}, {a}}], "
           "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{2}}}, {c}}]]");
    chk("SetAttributes[lf, Listable]; Null", "Null");
    chk_if("lf[SparseArray[{1 -> a}, 3], SparseArray[{2 -> c}, 3]]",
           "SparseArray[Automatic, {3}, lf[0, 0], {1, {{0, 2}, {{1}, {2}}}, {lf[a, 0], lf[0, c]}}]");
    chk_if("lf[SparseArray[{1 -> a}, 3], 2]",
           "SparseArray[Automatic, {3}, lf[0, 2], {1, {{0, 1}, {{1}}}, {lf[a, 2]}}]");
    chk("zf[_] := 0; SetAttributes[zf, Listable]; Null", "Null");
    chk_if("zf[SparseArray[{1 -> 1}, 3]]", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{1}}}, {0}}]");
    /* the kernel tier agrees with the boxed per-element evaluation */
    chk("v0 = SparseArray[{1 -> 0.5, 4 -> 2.5, 7 -> -1.25}, 9];", "Null");
    const char* fns[] = { "Sin", "Cos", "Exp", "Tanh", "ArcTan", "Erf", "Gamma", "Sqrt", "Log",
                          "BesselJ[0, #]&", "Floor", "Abs", NULL };
    for (int i = 0; fns[i]; i++) {
        char buf[256];
        snprintf(buf, sizeof buf, "With[{p = Normal[(%s)[v0]], q = (%s) /@ Normal[v0]}, p === q || p == q]",
                 fns[i], fns[i]);
        chk_true(buf);
    }
    chk_true("Normal[ArcTan[v0, 2.]] == ArcTan[Normal[v0], 2.]");
    chk_true("Normal[Sqrt[SparseArray[{1 -> -4., 2 -> 4.}, 3]]] == Sqrt[{-4., 4., 0}]");
}

/* The machine fast paths (packed construction through the radix sort, the
 * fused typed merges, the packed dense scan) against the boxed paths that
 * define the semantics, on pseudo-random data with repeated positions. */
static void test_fast_paths(void) {
    chk("k = Range[3000]; pi = Mod[Mod[k 48271, 2147483647], 50] + 1; "
        "pj = Mod[Mod[k 69621, 2147483647], 40] + 1; pk = Mod[k, 3] + 1; "
        "pv = N[Mod[k, 97]]/10. - 4.; pw = Mod[k, 13] - 6;", "Null");
    chk_true("Developer`PackedArrayQ[Transpose[{pi, pj}]] && Developer`PackedArrayQ[pv]");
    /* packed {pos} -> {vals} (radix, first wins) == the boxed rule list */
    chk_true("SparseArray[Transpose[{pi, pj}] -> pv, {50, 40}] === "
             "SparseArray[MapThread[Rule, {Transpose[{pi, pj}], pv}] /. r_Rule :> r, {50, 40}]");
    chk_true("SparseArray[Transpose[{pi, pj}] -> pw, {50, 40}] === "
             "SparseArray[Thread[Rule[Normal[Transpose[{pi, pj}]], Normal[pw]]], {50, 40}]");
    chk_true("SparseArray[Transpose[{pi, pj, pk}] -> pv] === "
             "SparseArray[Thread[Rule[Normal[Transpose[{pi, pj, pk}]], Normal[pv]]]]");
    chk_true("SparseArray[Transpose[{-pi, pj}] -> pv, {50, 40}] === "
             "SparseArray[Transpose[{51 - pi, pj}] -> pv, {50, 40}]");
    /* one long row: a vector built out of order */
    chk_true("SparseArray[Transpose[{Mod[k 7919, 100000] + 1}] -> pv, 100000] === "
             "SparseArray[Thread[Rule[Mod[Normal[k] 7919, 100000] + 1, Normal[pv]]], 100000]");
    chk_true("Quiet[Check[SparseArray[Transpose[{pi, pj}] -> pv, {10, 40}], $Failed]] === $Failed");
    /* fused merges against the dense arithmetic */
    chk("fa = SparseArray[Transpose[{pi, pj}] -> pv, {50, 40}]; "
        "fb = SparseArray[Transpose[{pj, pi}] -> pv, {40, 50}][[All, 1 ;; 40]] ; "
        "fb = SparseArray[fb, {50, 40}];", "Null");
    chk_true("Normal[fa + fb] === Normal[fa] + Normal[fb]");
    chk_true("Normal[fa - 2.5 fb] === Normal[fa] - 2.5 Normal[fb]");
    chk_true("Normal[fa fb] == Normal[fa] Normal[fb]");
    chk_true("Normal[(fa + 1.) (fb - 2.)] == (Normal[fa] + 1.) (Normal[fb] - 2.)");
    chk("ia = SparseArray[Transpose[{pi, pj}] -> pw, {50, 40}]; "
        "ib = SparseArray[Transpose[{pj, pi}] -> 3 pw, {40, 50}][[All, 1 ;; 40]]; "
        "ib = SparseArray[ib, {50, 40}];", "Null");
    chk_true("Normal[ia + ib] === Normal[ia] + Normal[ib]");
    chk_true("Normal[ia ib] == Normal[ia] Normal[ib]");
    chk_true("Normal[ia - 7] === Normal[ia] - 7");
    chk_true("Normal[ia + fa] == Normal[ia] + Normal[fa]");
    /* the packed dense scan agrees with the boxed one */
    chk("dm = N[Table[If[Mod[i j, 7] == 0, i - j, 0], {i, 30}, {j, 30}]];", "Null");
    chk_true("Developer`PackedArrayQ[dm] && SparseArray[dm] === SparseArray[Developer`FromPackedArray[dm]]");
    chk_true("Normal[SparseArray[dm]] === dm");
    chk_true("Developer`PackedArrayQ[Normal[SparseArray[dm]]]");
    chk_true("Developer`PackedArrayQ[SparseArray[dm][\"NonzeroPositions\"]]");
    chk_true("SparseArray[dm][\"NonzeroPositions\"] === Position[dm, x_ /; x != 0]");
    /* elementwise functions on a large typed array match the dense surface */
    chk_true("Normal[Sin[fa]] == Sin[Normal[fa]]");
    chk_true("Normal[Exp[fa]] == Exp[Normal[fa]]");
    chk_true("Normal[fa^2] == Normal[fa]^2");
    chk_true("Normal[2.^fa] == 2.^Normal[fa]");
}

/* The transparency gate (sparse_gate.c): heads that answer on Normal[s], heads
 * whose array answer comes back sparse, level-1 walkers, and the pattern heads
 * that keep treating the array as an atom. */
static void test_gate(void) {
    chk("S = SparseArray[{{1,1}->1,{2,3}->2,{3,2}->-1,{3,3}->4},{3,3}]; "
        "V = SparseArray[{1->3,4->-2},{5}]; Null", "Null");
    /* SPG_DENSE */
    chk("Total[S]", "{1, -1, 6}");
    chk("Total[S, 2]", "6");
    chk("Max[S]", "4");
    chk("Min[V]", "-2");
    chk("Tr[S]", "5");
    chk("Mean[V]", "1/5");
    chk("Accumulate[V]", "{3, 3, 3, 1, 1}");
    chk("Det[S]", "2");
    chk("Inverse[S]", "{{1, 0, 0}, {0, 2, -1}, {0, 1/2, 0}}");
    chk("MatrixRank[S]", "3");
    chk("S . {1, 2, 3}", "{1, 6, 10}");
    chk("{1, 2, 3} . S", "{1, -3, 16}");
    chk_true("S == Normal[S]");
    chk("Equal[S, {{1, 0, 2}, {0, 3, 0}, {4, 0, 5}}]", "False");
    /* SPG_RESPARSE: sparse when the array leads */
    chk_true("Head[Transpose[S]] === SparseArray && Normal[Transpose[S]] === Transpose[Normal[S]]");
    chk_true("Head[Reverse[V]] === SparseArray && Normal[Reverse[V]] === {0, -2, 0, 0, 3}");
    chk_true("Normal[Take[S, 2]] === {{1, 0, 0}, {0, 0, 2}}");
    chk_true("Normal[Flatten[S]] === {1, 0, 0, 0, 0, 2, 0, -1, 4}");
    chk_true("Head[First[S]] === SparseArray && Normal[Rest[V]] === {0, 0, -2, 0}");
    chk_true("Head[Join[S, S]] === SparseArray && Normal[Join[S, S]] === Join[Normal[S], Normal[S]]");
    chk("Head[Join[S, {{1, 0, 2}}]]", "List");
    chk_true("Head[Outer[Times, V, {1, 2}]] === SparseArray");
    chk_true("Head[ArrayFlatten[{{S, S}}]] === SparseArray && Dimensions[ArrayFlatten[{{S, S}}]] === {3, 6}");
    chk_true("Normal[Diagonal[S]] === {1, 0, 4}");
    /* Map: rows of a matrix, elements of a vector (and the background maps) */
    chk("Map[g, S] /. a_SparseArray :> Normal[a]", "{g[{1, 0, 0}], g[{0, 0, 2}], g[{0, -1, 4}]}");
    chk("Map[g, V][\"Background\"]", "g[0]");
    chk_true("Normal[Map[g, V]] === {g[3], g[0], g[0], g[-2], g[0]}");
    chk("Apply[Plus, S] // Normal", "{1, -1, 6}");
    /* pattern heads keep the atom */
    chk("Count[S, 1]", "0");
    chk("Position[S, 2]", "{}");
    chk("MemberQ[S, 2]", "False");
    chk("AtomQ[S]", "True");
    chk("Select[V, Positive]", "{3}");
    chk_true("Head[First[Select[S, VectorQ]]] === SparseArray");
    /* a declined pos -> vals rule with only ONE packed side keeps the other
     * side (it was replaced by NULL, and the next step crashed) */
    chk("Quiet[SparseArray[g[1] -> N[Range[3]], {4}]] /. SparseArray[r_, __] :> r",
        "g[1] -> {1.0, 2.0, 3.0}");
    chk("Quiet[SparseArray[Transpose[{Range[3], Range[3]}] -> {h, k}, {4, 4}]] /. SparseArray[r_, __] :> First[r]",
        "{{1, 1}, {2, 2}, {3, 3}}");
    /* a too-large implicit densification is refused, not attempted */
    chk("Quiet[Median[SparseArray[{{1, 1} -> 1}, {10^5, 10^5}]]] // Head", "Median");
    chk("Check[Median[SparseArray[{{1, 1} -> 1}, {10^5, 10^5}]], fail]", "fail");
}

int main(void) {
    symtab_init();
    core_init();
    test_constructors();
    test_messages();
    test_raw_object();
    test_structure();
    test_part();
    test_arithmetic();
    test_functions();
    test_fast_paths();
    test_gate();
    if (failures) {
        fprintf(stderr, "\n%d SparseArray test(s) FAILED\n", failures);
        return 1;
    }
    printf("All SparseArray tests passed.\n");
    return 0;
}
