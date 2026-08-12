/* Packed lists: the invisibility contract.
 *
 * A packed list (ToNDArray[...], and in due course anything the system packs on
 * its own) must be indistinguishable from the nested List it stands for. Most
 * of that falls out of the evaluator's transparency gate, but three functions
 * sit BELOW the evaluator and have to be got right by hand: expr_eq, expr_hash
 * and expr_compare. Each walks the buffer and reproduces what
 * ndarray_to_nested_list would have built, WITHOUT building it -- so a drift
 * between the two is silent, shows up only in a hash table or a Sort, and no
 * other test in the tree would catch it.
 *
 * So the assertions here are deliberately differential: they compare the packed
 * value against the REAL materialised form, computed at run time, rather than
 * against a hand-written expectation that could be wrong in the same direction
 * as the code. */

#include "test_utils.h"
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "ndarray.h"
#include "print_latex.h"
#include "pack.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static Expr* ev(const char* src) {
    Expr* p = parse_expression(src);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    expr_free(p);
    return r;
}

/* The core differential check: for a source expression that evaluates to a
 * packed list, every observable must agree with its own materialised form. */
static void check_invisible(const char* src) {
    Expr* packed = ev(src);
    ASSERT(is_packed_list(packed));

    Expr* plain = ndarray_to_nested_list(packed);
    ASSERT(plain != NULL);

    /* Identity: equal both ways, and hashing to the same bucket. A packed list
     * that hashes differently cannot be looked up in an Association keyed by
     * the plain form, and Union splits one value into two. */
    ASSERT(expr_eq(packed, plain));
    ASSERT(expr_eq(plain, packed));
    ASSERT(expr_hash(packed) == expr_hash(plain));
    ASSERT(expr_compare(packed, plain) == 0);
    ASSERT(expr_compare(plain, packed) == 0);

    /* Printing, in every form. */
    char* sp = expr_to_string(packed);
    char* sl = expr_to_string(plain);
    if (strcmp(sp, sl) != 0)
        fprintf(stderr, "FAIL: print mismatch for %s\n  packed: %s\n  plain:  %s\n",
                src, sp, sl);
    ASSERT(strcmp(sp, sl) == 0);
    free(sp); free(sl);

    char* fp = expr_to_string_fullform(packed);
    char* fl = expr_to_string_fullform(plain);
    if (strcmp(fp, fl) != 0)
        fprintf(stderr, "FAIL: fullform mismatch for %s\n  packed: %s\n  plain:  %s\n",
                src, fp, fl);
    ASSERT(strcmp(fp, fl) == 0);
    free(fp); free(fl);

    expr_free(plain);
    expr_free(packed);
}

void test_invisible_across_ranks_and_dtypes(void) {
    /* Rank 1-3, float64 and int64, including the values whose bit patterns are
     * easiest to get wrong (negative zero, denormals, large exact integers). */
    check_invisible("ToNDArray[{1., 2., 3.}]");
    check_invisible("ToNDArray[{1, 2, 3}]");
    check_invisible("ToNDArray[{{1., 2.}, {3., 4.}}]");
    check_invisible("ToNDArray[{{1, 2}, {3, 4}}]");
    check_invisible("ToNDArray[{{{1., 2.}, {3., 4.}}, {{5., 6.}, {7., 8.}}}]");
    check_invisible("ToNDArray[{-0., 0., 1.5, -2.25}]");
    check_invisible("ToNDArray[{9007199254740993, -9007199254740993}]");
    check_invisible("ToNDArray[Table[N[k], {k, 300}]]");
    check_invisible("ToNDArray[Table[k, {k, 300}]]");
}

/* ---------- What the user sees ---------- */

void test_presents_as_list(void) {
    assert_eval_eq("ToNDArray[{1., 2., 3.}]", "{1.0, 2.0, 3.0}", 0);
    assert_eval_eq("Head[ToNDArray[{1., 2., 3.}]]", "List", 0);
    assert_eval_eq("AtomQ[ToNDArray[{1., 2., 3.}]]", "False", 0);
    assert_eval_eq("NDArrayQ[ToNDArray[{1., 2., 3.}]]", "True", 0);
    assert_eval_eq("FullForm[ToNDArray[{1., 2.}]]", "List[1.0, 2.0]", 0);
    assert_eval_eq("InputForm[ToNDArray[{1., 2.}]]", "{1.0, 2.0}", 0);
    /* SameQ against the plain List, in both argument orders. */
    assert_eval_eq("ToNDArray[{1., 2., 3.}] === {1., 2., 3.}", "True", 0);
    assert_eval_eq("{1., 2., 3.} === ToNDArray[{1., 2., 3.}]", "True", 0);
}

void test_visible_ndarray_unchanged(void) {
    /* The explicit NDArray[...] head keeps every bit of its old behaviour: it
     * is a DIFFERENT value from the List, with a different Head, and it is
     * atomic. Packing must not have leaked into it. */
    assert_eval_eq("NDArray[{1., 2., 3.}]", "NDArray[{1.0, 2.0, 3.0}]", 0);
    assert_eval_eq("Head[NDArray[{1., 2., 3.}]]", "NDArray", 0);
    assert_eval_eq("AtomQ[NDArray[{1., 2., 3.}]]", "True", 0);
    assert_eval_eq("NDArray[{1., 2., 3.}] === {1., 2., 3.}", "False", 0);
    assert_eval_eq("NDArray[{1., 2., 3.}] === ToNDArray[{1., 2., 3.}]", "False", 0);
}

/* ---------- Exactness ---------- */

void test_exactness_preserved(void) {
    /* An all-Integer list packs to an int64 buffer and its elements come back
     * as Integers. Packing may not turn {1, 2, 3} into {1., 2., 3.}. */
    assert_eval_eq("DataType[ToNDArray[{1, 2, 3}]]", "\"int64\"", 0);
    assert_eval_eq("Head[ToNDArray[{1, 2, 3}][[2]]]", "Integer", 0);
    assert_eval_eq("ToNDArray[{1, 2, 3}] === {1, 2, 3}", "True", 0);
    assert_eval_eq("ToNDArray[{1, 2, 3}] === {1., 2., 3.}", "False", 0);
    /* Exact integers beyond 2^53, where the generic ndt_get/ndt_set pair (which
     * routes through double) would round. */
    assert_eval_eq("ToNDArray[{9007199254740993}][[1]]", "9007199254740993", 0);
}

/* An int64 buffer must give the SAME answer as the ordinary integer list, for
 * every head -- value, head and all.
 *
 * This is the one place packing came closest to shipping wrong answers. Most of
 * the ND layer reads elements through ndt_get, which routes via `double` and is
 * exact only to 2^53; that was safe while only Compile[] could create an
 * integer buffer, and stopped being safe the moment an all-Integer list packed
 * to one. Before the int64 gate, `Total[{1,2,3}]` came back as 6. instead of 6
 * and `Sin[{1,2,3}]` as {0,0,0} instead of symbolic.
 *
 * Written as a differential rather than as fixed expectations: each expression
 * is evaluated over the packed list and over the plain one, and the two printed
 * results must be identical. That way it keeps testing the right thing even if
 * an answer legitimately changes for unrelated reasons. */
static void same_as_plain(const char* fmt_packed, const char* fmt_plain) {
    Expr* a = ev(fmt_packed);
    Expr* b = ev(fmt_plain);
    char* sa = expr_to_string(a);
    char* sb = expr_to_string(b);
    if (strcmp(sa, sb) != 0)
        fprintf(stderr, "FAIL: packed and plain disagree\n  %s -> %s\n  %s -> %s\n",
                fmt_packed, sa, fmt_plain, sb);
    ASSERT(strcmp(sa, sb) == 0);
    free(sa); free(sb);
    expr_free(a); expr_free(b);
}

void test_int64_matches_plain_integer_lists(void) {
    static const char* const EXPRS[] = {
        "Total[%s]",  "Head[Total[%s]]", "Mean[%s]",   "Median[%s]",
        "Max[%s]",    "Min[%s]",         "Variance[%s]",
        "%s . %s",    "%s + 1",          "%s * 2",     "%s * 5/2",
        "Sin[%s]",    "Sqrt[%s]",        "Exp[%s]",    "Abs[%s]",
        "Sort[%s]",   "Reverse[%s]",     "Accumulate[%s]", "Differences[%s]",
        "Precision[%s]", "Length[%s]",   "Dimensions[%s]", "Depth[%s]",
        "%s[[2]]",    "Head[%s[[2]]]",   "Map[#^2 &, %s]", "Total[%s^3]",
        "Tally[%s]",  "Head[Tally[%s][[1,1]]]",
    };
    const char* packed_src = "ToNDArray[{1, 2, 3, 4}]";
    const char* plain_src  = "{1, 2, 3, 4}";
    char bufp[256], bufl[256];
    for (size_t i = 0; i < sizeof(EXPRS) / sizeof(EXPRS[0]); i++) {
        /* The one two-slot pattern is `%s . %s`; snprintf with a repeated
         * argument covers both shapes. */
        snprintf(bufp, sizeof(bufp), EXPRS[i], packed_src, packed_src);
        snprintf(bufl, sizeof(bufl), EXPRS[i], plain_src, plain_src);
        same_as_plain(bufp, bufl);
    }
    /* Exactness that a double buffer could not hold: the sum overflows past
     * what 2^53 can represent, and the interpreter promotes to a bigint. */
    same_as_plain("Total[ToNDArray[Table[k, {k, 1000000}]]]",
                  "Total[Table[k, {k, 1000000}]]");
    same_as_plain("Total[ToNDArray[{1000000000, 1000000000, 1000000000}]^3]",
                  "Total[{1000000000, 1000000000, 1000000000}^3]");
    /* Both operands int64 buffers -- ndarray_elementwise_power's exact path,
     * which the %s-templated cases above never reach (their exponent is either
     * a scalar or a Real list). The second pair overflows int64, which must
     * abandon to the List path and reach GMP rather than wrap. */
    same_as_plain("ToNDArray[{2, 3, 4, 5}]^ToNDArray[{3, 2, 5, 4}]",
                  "{2, 3, 4, 5}^{3, 2, 5, 4}");
    same_as_plain("ToNDArray[{1000000, 1000000}]^ToNDArray[{4, 4}]",
                  "{1000000, 1000000}^{4, 4}");
}

void test_bool_packing_matches_plain(void) {
    /* An all-True/False list packs to a one-byte bool buffer and is otherwise
     * indistinguishable from the plain list. */
    assert_eval_eq("DataType[ToNDArray[{True, False, True}]]", "\"bool\"", 0);
    assert_eval_eq("Head[ToNDArray[{True, False}]]", "List", 0);
    assert_eval_eq("ToNDArray[{True, False, True}] === {True, False, True}", "True", 0);
    /* A mixed bool/number list cannot form a uniform buffer, so it declines. */
    assert_eval_eq("NDArrayQ[ToNDArray[{True, 1}]]", "False", 0);
    assert_eval_eq("NDArrayQ[ToNDArray[{True, x}]]", "False", 0);

    /* Every observable of a packed bool list must equal the plain list's. The
     * producers (sign predicates) turn a numeric buffer INTO a bool one; the
     * consumers (Boole, the quantifiers) turn one back or reduce it. */
    static const char* const BOOL_EXPRS[] = {
        "%s",             "Head[%s]",       "Length[%s]",   "%s[[2]]",
        "Boole[%s]",      "Normal[%s]",     "Reverse[%s]",  "First[%s]",
        "AllTrue[%s, TrueQ]", "AnyTrue[%s, TrueQ]", "NoneTrue[%s, TrueQ]",
        "And @@ %s",      "Or @@ %s",       "Count[%s, True]",
    };
    const char* packed_src = "ToNDArray[{True, False, True}]";
    const char* plain_src  = "{True, False, True}";
    char bufp[256], bufl[256];
    for (size_t i = 0; i < sizeof(BOOL_EXPRS) / sizeof(BOOL_EXPRS[0]); i++) {
        snprintf(bufp, sizeof(bufp), BOOL_EXPRS[i], packed_src);
        snprintf(bufl, sizeof(bufl), BOOL_EXPRS[i], plain_src);
        same_as_plain(bufp, bufl);
    }
    /* Sign predicates: a numeric buffer in, a bool buffer out, matching the List
     * path element for element on both the packed and visible surfaces. */
    same_as_plain("Positive[ToNDArray[{-2, 0, 3, -1}]]", "Positive[{-2, 0, 3, -1}]");
    same_as_plain("NonNegative[ToNDArray[{-1., 0., 2.}]]", "NonNegative[{-1., 0., 2.}]");
    same_as_plain("Boole[Positive[ToNDArray[{-1, 2, -3}]]]", "Boole[Positive[{-1, 2, -3}]]");
}

void test_declines_what_it_cannot_represent(void) {
    /* Each of these must come back as the ORIGINAL list, unpacked. A mixed
     * exact/inexact list is the important one: a uniform buffer cannot hold an
     * Integer head on one element and a Real head on another, so widening it
     * would silently change {1, 2.5} into {1., 2.5}. */
    assert_eval_eq("NDArrayQ[ToNDArray[{1, 2.5}]]", "False", 0);
    assert_eval_eq("ToNDArray[{1, 2.5}]", "{1, 2.5}", 0);
    assert_eval_eq("NDArrayQ[ToNDArray[{1, 2, x}]]", "False", 0);
    assert_eval_eq("NDArrayQ[ToNDArray[{1/2, 1/3}]]", "False", 0);
    assert_eval_eq("NDArrayQ[ToNDArray[{{1., 2.}, {3.}}]]", "False", 0);   /* ragged */
    assert_eval_eq("NDArrayQ[ToNDArray[{}]]", "False", 0);
    assert_eval_eq("NDArrayQ[ToNDArray[{2^70, 1}]]", "False", 0);          /* BigInt */
    /* An explicit DataType may widen an exact list, but never round an inexact
     * one into an integer buffer -- that would change values, not storage. */
    assert_eval_eq("NDArrayQ[ToNDArray[{1, 2, 3}, DataType -> \"float64\"]]", "True", 0);
    assert_eval_eq("DataType[ToNDArray[{1, 2, 3}, DataType -> \"float64\"]]", "\"float64\"", 0);
    assert_eval_eq("NDArrayQ[ToNDArray[{1., 2.5}, DataType -> \"int64\"]]", "False", 0);
}

/* ---------- Ordering ---------- */

/* Sorting a packed list against a plain one must give the same answer as
 * sorting two plain ones. The buffer fast path is only valid when the shapes
 * match exactly: element-wise, {2., 0.} sorts AFTER {1., 9., 9.} on its first
 * element but BEFORE it on length, and List order settles length first. */
/* Tally hashes the machine word rather than the Expr (src/ndreduce.c), which is
 * the one fast path here that builds its OWN keys instead of moving elements.
 * Everything that can differ between keying on a word and keying on an
 * expression is checked against the plain list: the first-appearance ordering,
 * the element heads of the reconstructed values, integers past 2^53 (a float64
 * gather would merge them), the two zeros, and the custom-test form -- which
 * must not take the fast path at all, since a user test has to see expressions.
 */
void test_tally_matches_plain_lists(void) {
    static const char* const CASES[][2] = {
        /* Repeats, and first-appearance order -- NOT sorted order. */
        { "Tally[ToNDArray[{3, 1, 3, 2, 1, 3}]]",   "Tally[{3, 1, 3, 2, 1, 3}]" },
        { "Tally[ToNDArray[{5, 5, 5, 5}]]",         "Tally[{5, 5, 5, 5}]" },
        { "Tally[ToNDArray[{-2, 0, -2, 7}]]",       "Tally[{-2, 0, -2, 7}]" },
        /* Every value distinct, and every value identical: the two ends of the
         * table-growth path. */
        { "Tally[ToNDArray[Range[500]]]",           "Tally[Range[500]]" },
        { "Tally[ToNDArray[ConstantArray[4, 500]]]","Tally[ConstantArray[4, 500]]" },
        { "Length[Tally[ToNDArray[Mod[Range[1000]^2, 37]]]]",
          "Length[Tally[Mod[Range[1000]^2, 37]]]" },
        /* Past 2^53: 9007199254740993 and 9007199254740992 are the same double
         * and must stay two tallies. */
        { "Tally[ToNDArray[{9007199254740993, 9007199254740992, 9007199254740993}]]",
          "Tally[{9007199254740993, 9007199254740992, 9007199254740993}]" },
        /* Reals, including the two zeros, which are ONE value in both paths. */
        { "Tally[ToNDArray[{3.5, 1.5, 3.5}]]",      "Tally[{3.5, 1.5, 3.5}]" },
        { "Tally[ToNDArray[{0., -0., 1.}]]",        "Tally[{0., -0., 1.}]" },
        { "Tally[ToNDArray[{1., 2., 1., 2., 3.}]]", "Tally[{1., 2., 1., 2., 3.}]" },
        /* Heads of the rebuilt values and of the counts. */
        { "Head[Tally[ToNDArray[{1, 2, 1}]][[1, 1]]]",   "Head[Tally[{1, 2, 1}][[1, 1]]]" },
        { "Head[Tally[ToNDArray[{1., 2., 1.}]][[1, 1]]]","Head[Tally[{1., 2., 1.}][[1, 1]]]" },
        { "Head[Tally[ToNDArray[{1, 2, 1}]][[1, 2]]]",   "Head[Tally[{1, 2, 1}][[1, 2]]]" },
        /* A custom test must reach the ordinary implementation. Handing an
         * NDArray to it unchanged answered {} -- the generic code indexes its
         * argument as a List and an NDArray is not one. */
        { "Tally[ToNDArray[{1, 2, 3}], Divisible[#2, #1] &]",
          "Tally[{1, 2, 3}, Divisible[#2, #1] &]" },
        { "Tally[ToNDArray[{1., 2., 4.}], Abs[#1 - #2] < 3 &]",
          "Tally[{1., 2., 4.}, Abs[#1 - #2] < 3 &]" },
        /* Rank 2 tallies ROWS, which are not machine words. */
        { "Tally[ToNDArray[{{1, 2}, {1, 2}, {3, 4}}]]", "Tally[{{1, 2}, {1, 2}, {3, 4}}]" },
        /* The direct-index path (src/ndreduce.c) is chosen on the value RANGE,
         * so both sides of that decision need pinning. It is taken when the
         * range is no wider than the input and declines to the hash otherwise;
         * the two must agree, including on the exact boundary. */
        { "Tally[ToNDArray[{0, 4, 2, 4, 0}]]",  "Tally[{0, 4, 2, 4, 0}]" },   /* range == n */
        { "Tally[ToNDArray[{0, 100, 200}]]",    "Tally[{0, 100, 200}]" },     /* range > n  */
        { "Tally[ToNDArray[{-1000000000000, 0, 1000000000000, 0}]]",
          "Tally[{-1000000000000, 0, 1000000000000, 0}]" },                   /* hash path  */
        /* A negative minimum: the offset is computed in unsigned arithmetic
         * precisely so this cannot wrap. */
        { "Tally[ToNDArray[{-5, -3, -5, -4}]]", "Tally[{-5, -3, -5, -4}]" },
    };
    for (size_t i = 0; i < sizeof(CASES) / sizeof(CASES[0]); i++)
        same_as_plain(CASES[i][0], CASES[i][1]);
}

/* DeleteDuplicates / Union / Intersection / Complement over packed int64
 * (src/list/setops.c). Like Tally these build their own answer rather than
 * moving elements, and each is chosen by the value RANGE: a bounded range is
 * direct-indexed, anything wider falls back (a hash set for DeleteDuplicates, a
 * sort-merge for the rest). Both sides of every such decision are pinned here
 * against the plain-list result.
 *
 * The custom-SameTest rows are the important ones. Putting DeleteDuplicates on
 * pack.c's AWARE list stops the transparency gate materialising for EVERY call,
 * not just the ones the fast path handles -- and the generic implementation
 * tests `type != EXPR_FUNCTION`, which an NDArray is not, so it returned its
 * argument unchanged and deduplicated nothing. Silent, and no gate can see it. */
void test_setops_match_plain_lists(void) {
    static const char* const CASES[][2] = {
        /* First-appearance order, which is DeleteDuplicates' contract. */
        { "DeleteDuplicates[ToNDArray[{5, 1, 5, 3, 1, 9}]]", "DeleteDuplicates[{5, 1, 5, 3, 1, 9}]" },
        { "DeleteDuplicates[ToNDArray[{4, 4, 4, 4}]]",       "DeleteDuplicates[{4, 4, 4, 4}]" },
        { "DeleteDuplicates[ToNDArray[{-2, -5, -2, 0}]]",    "DeleteDuplicates[{-2, -5, -2, 0}]" },
        { "DeleteDuplicates[ToNDArray[Range[500]]]",         "DeleteDuplicates[Range[500]]" },
        /* Range too wide to index -- the inline-key hash set. */
        { "DeleteDuplicates[ToNDArray[{0, 1000000000000, 0, -1000000000000}]]",
          "DeleteDuplicates[{0, 1000000000000, 0, -1000000000000}]" },
        /* Past 2^53: must stay distinct, as for Tally. */
        { "DeleteDuplicates[ToNDArray[{9007199254740993, 9007199254740992, 9007199254740993}]]",
          "DeleteDuplicates[{9007199254740993, 9007199254740992, 9007199254740993}]" },
        /* A custom SameTest must reach the ordinary implementation. */
        { "DeleteDuplicates[ToNDArray[{1, 2, 3, 4}], Mod[#1 - #2, 2] == 0 &]",
          "DeleteDuplicates[{1, 2, 3, 4}, Mod[#1 - #2, 2] == 0 &]" },
        { "DeleteDuplicates[ToNDArray[{2, 4, 5}], Divisible[#2, #1] &]",
          "DeleteDuplicates[{2, 4, 5}, Divisible[#2, #1] &]" },
        /* Union of one list is sorted-unique; the range walk emits it ascending
         * without sorting, so the ORDER is the thing to check. */
        { "Union[ToNDArray[{5, 1, 5, 3, 1, 9}]]", "Union[{5, 1, 5, 3, 1, 9}]" },
        { "Union[ToNDArray[{-3, 7, -3}]]",        "Union[{-3, 7, -3}]" },
        { "Union[ToNDArray[{0, 1000000000000}]]", "Union[{0, 1000000000000}]" },
        /* Two- and three-list forms: one bitmask per value, same walk, different
         * predicate. Disjoint, identical, and empty operands are the edges. */
        { "Union[ToNDArray[{3, 1}], ToNDArray[{2, 1}]]", "Union[{3, 1}, {2, 1}]" },
        { "Union[ToNDArray[{3}], ToNDArray[{2}], ToNDArray[{1}]]", "Union[{3}, {2}, {1}]" },
        { "Intersection[ToNDArray[{3, 1, 2}], ToNDArray[{2, 1}]]", "Intersection[{3, 1, 2}, {2, 1}]" },
        { "Intersection[ToNDArray[{1, 2}], ToNDArray[{5, 6}]]",    "Intersection[{1, 2}, {5, 6}]" },
        { "Intersection[ToNDArray[{1, 2, 3}], ToNDArray[{2, 3}], ToNDArray[{3}]]",
          "Intersection[{1, 2, 3}, {2, 3}, {3}]" },
        { "Complement[ToNDArray[{3, 1, 2}], ToNDArray[{2}]]",   "Complement[{3, 1, 2}, {2}]" },
        { "Complement[ToNDArray[{1, 2}], ToNDArray[{1, 2}]]",   "Complement[{1, 2}, {1, 2}]" },
        { "Complement[ToNDArray[{1, 2, 3}], ToNDArray[{2}], ToNDArray[{3}]]",
          "Complement[{1, 2, 3}, {2}, {3}]" },
        /* An empty operand on either side. */
        { "Union[ToNDArray[{3, 1}], ToNDArray[{}]]",        "Union[{3, 1}, {}]" },
        { "Intersection[ToNDArray[{3, 1}], ToNDArray[{}]]", "Intersection[{3, 1}, {}]" },
        { "Complement[ToNDArray[{3, 1}], ToNDArray[{}]]",   "Complement[{3, 1}, {}]" },
        { "Complement[ToNDArray[{}], ToNDArray[{3}]]",      "Complement[{}, {3}]" },
        /* Element heads survive the round trip. */
        { "Head[DeleteDuplicates[ToNDArray[{1, 2, 1}]][[1]]]", "Head[DeleteDuplicates[{1, 2, 1}][[1]]]" },
        { "Head[Union[ToNDArray[{2, 1, 2}]][[1]]]",            "Head[Union[{2, 1, 2}][[1]]]" },
    };
    for (size_t i = 0; i < sizeof(CASES) / sizeof(CASES[0]); i++)
        same_as_plain(CASES[i][0], CASES[i][1]);
}

/* Subtract and Divide became packed-aware on 2026-08-02, and a symbolic operand
 * beside an invisible packed array stopped being a wrong answer at the same
 * time.  Both were found by the gate pass of tools/nd_fastpath_sweep.py.
 *
 * The two are one test because the second was uncovered by checking the first:
 * Subtract and Divide rewrite to Plus/Times/Power, so marking them aware makes
 * them inherit whatever those heads do with an array — and what those heads did
 * with `packedArray + x` was WARN AND LEAVE IT UNEVALUATED, where the same short
 * list threads.  Packing keys on length, so one expression answered two ways:
 *
 *     Range[1., 3.]   + x  ->  {1. + x, 2. + x, 3. + x}
 *     Range[1., 300.] + x  ->  unevaluated, with NDArray::sym
 *
 * A visible NDArray[...] still declines, which is its documented contract; only
 * the invisible surface materialises and re-runs. */
void test_arithmetic_rewrites_match_plain_lists(void) {
    static const char* const CASES[][2] = {
        /* Subtract / Divide: the head the gate used to fire at is one rewrite
         * away from Plus and Times, which have had buffer paths all along. */
        { "Take[Range[1., 300.] - 1., 3]",        "Range[1., 3.] - 1." },
        { "Take[1. - Range[1., 300.], 3]",        "1. - Range[1., 3.]" },
        { "Take[Range[1., 300.]/2., 3]",          "Range[1., 3.]/2." },
        { "Take[2./Range[1., 300.], 3]",          "2./Range[1., 3.]" },
        { "Take[Subtract[Range[1., 300.], 1.], 3]", "Subtract[Range[1., 3.], 1.]" },
        /* THE WRONG ANSWER.  builtin_divide's Real branch fires when EITHER
         * operand is Real and then reads the other through scalar type tests
         * that an array fails all of, defaulting to 0.0 -- so this answered
         * `0.` for the whole call once the gate stopped materialising. */
        { "Take[Divide[Range[1., 300.], 2.], 3]", "Divide[Range[1., 3.], 2.]" },
        { "Head[Divide[Range[1., 300.], 2.]]",    "Head[Divide[Range[1., 3.], 2.]]" },
        /* Integer buffers: Subtract stays exact, Divide gives Rationals and so
         * must NOT keep the buffer (it is deliberately absent from INT64_OK). */
        { "Take[Range[300] - 1, 3]",              "Range[3] - 1" },
        { "Take[Range[300]/2, 3]",                "Range[3]/2" },
        { "Head[First[Range[300] - 1]]",          "Head[First[Range[3] - 1]]" },
        { "Head[First[Range[300]/2]]",            "Head[First[Range[3]/2]]" },
        /* Array against array, and the degenerate divisor. */
        { "Take[Range[1., 300.] - Range[1., 300.], 3]", "Range[1., 3.] - Range[1., 3.]" },
        { "Take[Range[1., 300.]/Range[1., 300.], 3]",   "Range[1., 3.]/Range[1., 3.]" },
        { "Take[Range[1., 300.]/0, 3]",           "Range[1., 3.]/0" },
        /* A SYMBOLIC operand must thread, not decline — for every head that
         * reaches ndarray_symbolic_delist_retry. */
        { "Take[Range[1., 300.] + x, 3]",         "Range[1., 3.] + x" },
        { "Take[Range[1., 300.] x, 3]",           "Range[1., 3.] x" },
        { "Take[Range[1., 300.] - x, 3]",         "Range[1., 3.] - x" },
        { "Take[Range[1., 300.]/x, 3]",           "Range[1., 3.]/x" },
        { "Take[Range[1., 300.]^x, 3]",           "Range[1., 3.]^x" },
        { "Take[x^Range[1., 300.], 3]",           "x^Range[1., 3.]" },
        { "Take[Range[300] + x, 3]",              "Range[3] + x" },
        /* An exact irrational is numeric, not symbolic: it must still fold. */
        { "Take[Range[1., 300.] + Sqrt[2], 3]",   "Range[1., 3.] + Sqrt[2]" },
        /* Rank 2, nesting, a compound symbolic term, and two symbols — the
         * materialise-and-re-run has to be the ORDINARY evaluation, not a
         * special case that only handles a bare symbol at rank 1. */
        { "Table[i*1., {i, 20}, {j, 20}][[1, 1]] + x", "1. + x" },
        { "Take[(Range[1., 300.] + x)*2, 3]",     "(Range[1., 3.] + x)*2" },
        { "Take[Range[1., 300.] + f[y], 3]",      "Range[1., 3.] + f[y]" },
        { "Take[Range[1., 300.] + x + y, 3]",     "Range[1., 3.] + x + y" },
        { "Take[Sin[Range[1., 300.] + x], 3]",    "Sin[Range[1., 3.] + x]" },
    };
    for (size_t i = 0; i < sizeof(CASES) / sizeof(CASES[0]); i++)
        same_as_plain(CASES[i][0], CASES[i][1]);

    /* The results still pack — the point of the change was speed, and a fix
     * that quietly stopped packing would pass every assertion above. */
    assert_eval_eq("{NDArrayQ[Range[1., 300.] - 1.], NDArrayQ[Range[1., 300.]/2.], "
                   "NDArrayQ[Range[300] - 1]}", "{True, True, True}", 0);

    /* A VISIBLE NDArray keeps the numeric-only contract: naming the head is
     * asking for a purely numeric object, and the warning is the answer.  One
     * visible operand anywhere in the call vetoes the retry, so marking
     * Subtract and Divide aware cannot leak the new behaviour into a mixed
     * call -- the alternative was tried, and it moved the divergence instead
     * of removing it. */
    assert_eval_eq("Head[NDArray[{1., 2., 3.}] + x]", "Plus", 0);
    assert_eval_eq("Head[NDArray[{1., 2., 3.}] - x]", "Plus", 0);
    assert_eval_eq("Head[NDArray[{1., 2., 3.}]/x]", "Times", 0);

    /* KNOWN OPEN, and pre-existing: mixing a VISIBLE NDArray with a list and a
     * symbol still answers two ways, because Plus is Listable and threads over
     * a plain List while a packed one is an EXPR_NDARRAY that apply_listable
     * does not recognise as threadable:
     *
     *   NDArray[{1., 2.}] + {1., 2.}                 + x  ->  threads
     *   NDArray[{1., 2.}] + Take[Range[1., 300.], 2] + x  ->  unevaluated
     *
     * Deliberately NOT asserted either way: closing it means deciding what
     * `NDArray + List + symbol` should mean, which is a semantic question and
     * not this change's to answer.  Recorded so the next sweep finds it
     * already known rather than rediscovering it. */
}

/* Commonest shares Tally's counting routine (ndred_tally / ndred_commonest in
 * src/ndreduce.c), so what needs pinning is the SELECTION: which distinct values
 * a count tie picks, and in what order they come back. Both are decided by
 * first-appearance index, and both surfaces must decide them identically. */
void test_commonest_matches_plain_lists(void) {
    static const char* const CASES[][2] = {
        /* A tie at the top count: BOTH survivors, in first-appearance order. */
        { "Commonest[ToNDArray[{5, 1, 5, 3, 1, 9}]]", "Commonest[{5, 1, 5, 3, 1, 9}]" },
        { "Commonest[ToNDArray[{3, 1, 1, 3, 2}]]",    "Commonest[{3, 1, 1, 3, 2}]" },
        /* No tie, all distinct, and a single repeated value. */
        { "Commonest[ToNDArray[{1, 2, 3}]]",          "Commonest[{1, 2, 3}]" },
        { "Commonest[ToNDArray[{4, 4, 4, 4}]]",       "Commonest[{4, 4, 4, 4}]" },
        { "Commonest[ToNDArray[{-2, -5, -2, 0}]]",    "Commonest[{-2, -5, -2, 0}]" },
        { "Commonest[ToNDArray[{7}]]",                "Commonest[{7}]" },
        /* The count forms. n past the distinct count is clamped (and warns);
         * UpTo asks for "at most" and is silent; zero or less selects nothing. */
        { "Commonest[ToNDArray[{5, 1, 5, 3, 1, 9}], 2]", "Commonest[{5, 1, 5, 3, 1, 9}, 2]" },
        { "Commonest[ToNDArray[{5, 1, 5, 3, 1, 9}], 4]", "Commonest[{5, 1, 5, 3, 1, 9}, 4]" },
        { "Commonest[ToNDArray[{1, 1, 2}], UpTo[7]]",    "Commonest[{1, 1, 2}, UpTo[7]]" },
        { "Commonest[ToNDArray[{1, 1, 2}], UpTo[1]]",    "Commonest[{1, 1, 2}, UpTo[1]]" },
        { "Commonest[ToNDArray[{1, 1, 2}], 0]",          "Commonest[{1, 1, 2}, 0]" },
        { "Commonest[ToNDArray[{1, 1, 2}], -1]",         "Commonest[{1, 1, 2}, -1]" },
        { "Commonest[ToNDArray[{1, 1, 2}], -2]",         "Commonest[{1, 1, 2}, -2]" },
        /* A count that is neither an Integer nor UpTo[Integer] is not the fast
         * path's to interpret: both surfaces leave the call unevaluated. */
        { "Commonest[ToNDArray[{1, 1, 2}], nn]",         "Commonest[{1, 1, 2}, nn]" },
        { "Commonest[ToNDArray[{1, 1, 2}], 2.5]",        "Commonest[{1, 1, 2}, 2.5]" },
        /* Range too wide to direct-index -- the inline-key hash. */
        { "Commonest[ToNDArray[{0, 1000000000000, 0, -1000000000000}]]",
          "Commonest[{0, 1000000000000, 0, -1000000000000}]" },
        /* Past 2^53: keyed on the raw word, so these stay three distinct values
         * and the answer is the one that really does appear twice. */
        { "Commonest[ToNDArray[{9007199254740993, 9007199254740992, 9007199254740993}]]",
          "Commonest[{9007199254740993, 9007199254740992, 9007199254740993}]" },
        /* Reals. 0. and -0. print differently and expr_eq holds them apart, so
         * the bit-pattern key must NOT normalise them -- the same trap Tally
         * has. Here it is visible in the ANSWER, not just in the count. */
        { "Commonest[ToNDArray[{1.5, 2.5, 1.5}]]",  "Commonest[{1.5, 2.5, 1.5}]" },
        { "Commonest[ToNDArray[{0., -0., 1.}]]",    "Commonest[{0., -0., 1.}]" },
        { "Commonest[ToNDArray[{1., 2., 3.}], 2]",  "Commonest[{1., 2., 3.}, 2]" },
        /* A non-finite element cannot be keyed on its bits and stay faithful, so
         * the whole call is handed back. */
        { "Commonest[ToNDArray[Exp[{1000., 1000., 2.}]]]", "Commonest[Exp[{1000., 1000., 2.}]]" },
        /* Rank 2: the ELEMENTS are rows, which are not machine words. */
        { "Commonest[ToNDArray[{{1, 2}, {1, 2}, {3, 4}}]]", "Commonest[{{1, 2}, {1, 2}, {3, 4}}]" },
        /* Element heads survive the round trip, for both dtypes. */
        { "Head[Commonest[ToNDArray[{1, 2, 1}]][[1]]]",    "Head[Commonest[{1, 2, 1}][[1]]]" },
        { "Head[Commonest[ToNDArray[{1., 2., 1.}]][[1]]]", "Head[Commonest[{1., 2., 1.}][[1]]]" },
    };
    for (size_t i = 0; i < sizeof(CASES) / sizeof(CASES[0]); i++)
        same_as_plain(CASES[i][0], CASES[i][1]);
}

void test_ordering_matches_plain_lists(void) {
    assert_eval_eq("Sort[{ToNDArray[{2., 1.}], ToNDArray[{1., 9.}]}]",
                   "{{1.0, 9.0}, {2.0, 1.0}}", 0);
    assert_eval_eq("Sort[{{2., 1.}, {1., 9.}}]",
                   "{{1.0, 9.0}, {2.0, 1.0}}", 0);
    /* Differing shapes: the fast path must decline and materialise. */
    assert_eval_eq("Sort[{ToNDArray[{2., 0.}], {1., 9., 9.}}]",
                   "{{2.0, 0.0}, {1.0, 9.0, 9.0}}", 0);
    assert_eval_eq("Sort[{{2., 0.}, {1., 9., 9.}}]",
                   "{{2.0, 0.0}, {1.0, 9.0, 9.0}}", 0);
    /* An int64 buffer must order exactly above 2^53, where comparing through
     * double would call two distinct values equal and break the total order
     * that Sort and Union depend on. */
    assert_eval_eq("Sort[{ToNDArray[{9007199254740993}], ToNDArray[{9007199254740992}]}]",
                   "{{9007199254740992}, {9007199254740993}}", 0);
}

/* The Ordering BUILTIN (distinct from test_ordering_matches_plain_lists above,
 * which pins Sort's canonical order of packed sublists). Ordering[a] argsorts the
 * buffer and returns an int64 permutation; every case must agree with the plain
 * List Ordering, and the int64 argsort must be exact past 2^53. */
void test_ordering_builtin_matches_plain(void) {
    same_as_plain("Ordering[ToNDArray[{2, 6, 1, 9, 1, 2, 3}]]",
                  "Ordering[{2, 6, 1, 9, 1, 2, 3}]");
    same_as_plain("Ordering[ToNDArray[{2., 6., 1., 9., 1., 2., 3.}]]",
                  "Ordering[{2., 6., 1., 9., 1., 2., 3.}]");
    /* The 2nd-argument Take-spec forms, on the buffer. */
    same_as_plain("Ordering[ToNDArray[{2, 6, 1, 9, 1, 2, 3}], 4]",
                  "Ordering[{2, 6, 1, 9, 1, 2, 3}, 4]");
    same_as_plain("Ordering[ToNDArray[{2, 6, 1, 9, 1, 2, 3}], -2]",
                  "Ordering[{2, 6, 1, 9, 1, 2, 3}, -2]");
    same_as_plain("Ordering[ToNDArray[{2, 6, 1, 9, 1, 2, 3}], {4, -1}]",
                  "Ordering[{2, 6, 1, 9, 1, 2, 3}, {4, -1}]");
    same_as_plain("Ordering[ToNDArray[{2., 6., 1., 9., 2.}], UpTo[6]]",
                  "Ordering[{2., 6., 1., 9., 2.}, UpTo[6]]");
    /* Result element is Integer regardless of the input dtype. */
    assert_eval_eq("Head[Ordering[ToNDArray[{1., 2., 1.}]][[1]]]", "Integer", 0);
    assert_eval_eq("Head[Ordering[ToNDArray[{1, 2, 1}]][[1]]]", "Integer", 0);
    /* int64 exactness past 2^53: comparing through a double would call these two
     * equal and stably keep {1, 2}; the int64 argsort orders them {2, 1}. */
    assert_eval_eq("Ordering[ToNDArray[{9007199254740993, 9007199254740992}]]",
                   "{2, 1}", 0);
    /* Ordering is an aware head: its packed input yields a packed permutation. */
    {
        Expr* r = ev("Ordering[ToNDArray[Range[300]]]");
        ASSERT(is_packed_list(r));
        expr_free(r);
    }
    /* Rank 2 declines to the List path (rows order by canonical Expr order), and
     * must still agree with the plain form. */
    same_as_plain("Ordering[ToNDArray[{{2, 1}, {1, 9}, {1, 2}}]]",
                  "Ordering[{{2, 1}, {1, 9}, {1, 2}}]");
}

void test_hash_consumers(void) {
    /* These are the reason expr_hash has to agree bit for bit. */
    assert_eval_eq("Union[{ToNDArray[{1., 2.}], {1., 2.}}]", "{{1.0, 2.0}}", 0);
    assert_eval_eq("DeleteDuplicates[{ToNDArray[{1., 2.}], {1., 2.}}]", "{{1.0, 2.0}}", 0);
    assert_eval_eq("<|ToNDArray[{1., 2.}] -> \"x\"|>[{1., 2.}]", "\"x\"", 0);
    assert_eval_eq("<|{1., 2.} -> \"x\"|>[ToNDArray[{1., 2.}]]", "\"x\"", 0);
}

/* ---------- Storage semantics ---------- */

void test_value_semantics_on_assignment(void) {
    /* Copy-on-write: a packed list assigned to a second symbol and then mutated
     * through that symbol must not disturb the first. The buffer is shared by
     * refcount until written, so a missing unshare here would be invisible
     * until exactly this test. */
    assert_eval_eq("Module[{p, q}, p = ToNDArray[{1., 2., 3.}]; q = p; "
                   "q[[1]] = 9.; {p, q}]",
                   "{{1.0, 2.0, 3.0}, {9.0, 2.0, 3.0}}", 0);
}

/* Writing into a packed list must give the same answer as writing into the
 * plain one, INCLUDING each element's head. The buffer path coerces, which is
 * right for a visible NDArray[...] but would be a silent answer change here:
 * p[[1]] = 1 into a float64 buffer reads back as 1., where the ordinary list
 * keeps the exact Integer. So a right-hand side that cannot be stored with its
 * head intact unpacks first.
 *
 * Each case is asserted against the plain-list answer written out in full,
 * because that is exactly what a user would compare against. */
void test_part_assignment_preserves_heads(void) {
    /* exact Integer into a float64 buffer */
    assert_eval_eq("Module[{p}, p = ToNDArray[{1., 2., 3.}]; p[[1]] = 1; p]",
                   "{1, 2.0, 3.0}", 0);
    assert_eval_eq("Module[{p}, p = {1., 2., 3.}; p[[1]] = 1; p]",
                   "{1, 2.0, 3.0}", 0);
    /* Rational into a float64 buffer */
    assert_eval_eq("Module[{p}, p = ToNDArray[{1., 2., 3.}]; p[[1]] = 1/3; p]",
                   "{1/3, 2.0, 3.0}", 0);
    /* Real into an int64 buffer */
    assert_eval_eq("Module[{p}, p = ToNDArray[{1, 2, 3}]; p[[1]] = 2.5; p]",
                   "{2.5, 2, 3}", 0);
    /* BigInt into an int64 buffer */
    assert_eval_eq("Module[{p}, p = ToNDArray[{1, 2, 3}]; p[[1]] = 2^70; p]",
                   "{1180591620717411303424, 2, 3}", 0);
    /* A symbolic value unpacks too. */
    assert_eval_eq("Module[{p}, p = ToNDArray[{1., 2., 3.}]; p[[1]] = zz; p]",
                   "{zz, 2.0, 3.0}", 0);
    /* The representable case still takes the buffer path and stays packed. */
    assert_eval_eq("Module[{p}, p = ToNDArray[{1., 2., 3.}]; p[[1]] = 9.; p]",
                   "{9.0, 2.0, 3.0}", 0);
    assert_eval_eq("Module[{p}, p = ToNDArray[{1., 2., 3.}]; p[[1]] = 9.; NDArrayQ[p]]",
                   "True", 0);
    /* A visible NDArray[...] keeps its coercing behaviour: the user asked for a
     * machine buffer, and its dtype is part of the value. */
    assert_eval_eq("Module[{p}, p = NDArray[{1., 2., 3.}]; p[[1]] = 1; p]",
                   "NDArray[{1.0, 2.0, 3.0}]", 0);
}

void test_part_inherits_presentation(void) {
    /* A sub-array carved out of a packed list is itself a packed LIST, not a
     * visible NDArray -- this is the expr_new_ndarray_like propagation. If it
     * regresses, transparency breaks at the first Part of a matrix. */
    assert_eval_eq("ToNDArray[{{1., 2.}, {3., 4.}}][[2]]", "{3.0, 4.0}", 0);
    assert_eval_eq("Head[ToNDArray[{{1., 2.}, {3., 4.}}][[2]]]", "List", 0);
    assert_eval_eq("NDArrayQ[ToNDArray[{{1., 2.}, {3., 4.}}][[2]]]", "True", 0);
    assert_eval_eq("ToNDArray[{{1., 2.}, {3., 4.}}][[2, 1]]", "3.0", 0);
}

void test_roundtrip_builtins(void) {
    assert_eval_eq("FromNDArray[ToNDArray[{1., 2., 3.}]]", "{1.0, 2.0, 3.0}", 0);
    assert_eval_eq("NDArrayQ[FromNDArray[ToNDArray[{1., 2., 3.}]]]", "False", 0);
    assert_eval_eq("Head[FromNDArray[ToNDArray[{1., 2., 3.}]]]", "List", 0);
    /* FromNDArray also undoes an explicit NDArray[...]. */
    assert_eval_eq("FromNDArray[NDArray[{1., 2.}]]", "{1.0, 2.0}", 0);
    /* And leaves anything else alone. */
    assert_eval_eq("FromNDArray[{1, 2, x}]", "{1, 2, x}", 0);
    assert_eval_eq("FromNDArray[7]", "7", 0);
    /* ToNDArray of an already-packed list is a no-op, not a double wrap. */
    assert_eval_eq("NDArrayQ[ToNDArray[ToNDArray[{1., 2., 3.}]]]", "True", 0);
    assert_eval_eq("ToNDArray[ToNDArray[{1., 2., 3.}]]", "{1.0, 2.0, 3.0}", 0);
    /* ToNDArray of a visible NDArray restates it as a List. */
    assert_eval_eq("Head[ToNDArray[NDArray[{1., 2.}]]]", "List", 0);
}

void test_kill_switch(void) {
    /* pack_set_enabled(false) must stop automatic packing without changing any
     * answer -- that property is what the differential suite rests on. It does
     * NOT disable ToNDArray, which is an explicit request. */
    pack_set_enabled(false);
    assert_eval_eq("NDArrayQ[ToNDArray[{1., 2., 3.}]]", "True", 0);
    assert_eval_eq("ToNDArray[{1., 2., 3.}]", "{1.0, 2.0, 3.0}", 0);
    pack_set_enabled(true);
}

/* ---------- The transparency gate ---------- */

/* Heads that know nothing about packing must still be right. They get a
 * materialised List from the gate in evaluate_step, so none of them needed a
 * line of packing-specific code -- which is the whole point, and also why this
 * has to be tested: if the gate regresses, every one of these silently returns
 * a confident wrong answer (Count would say 0) rather than crashing. */
void test_unaware_heads_are_correct(void) {
    assert_eval_eq("Count[ToNDArray[{1., 2., 3., 4.}], _Real]", "4", 0);
    assert_eval_eq("Cases[ToNDArray[{1., 2., 3., 4.}], x_ /; x > 2.]", "{3.0, 4.0}", 0);
    assert_eval_eq("Position[ToNDArray[{1., 2., 3.}], 2.]", "{{2}}", 0);
    assert_eval_eq("Level[ToNDArray[{1., 2.}], {1}]", "{1.0, 2.0}", 0);
    assert_eval_eq("ToNDArray[{1., 2., 3.}] /. 2. -> 9.", "{1.0, 9.0, 3.0}", 0);
    assert_eval_eq("LeafCount[ToNDArray[{1., 2., 3.}]]", "4", 0);
    assert_eval_eq("Insert[ToNDArray[{1., 2.}], 9., 2]", "{1.0, 9.0, 2.0}", 0);
    assert_eval_eq("Delete[ToNDArray[{1., 2., 3.}], 1]", "{2.0, 3.0}", 0);
    assert_eval_eq("ReplacePart[ToNDArray[{1., 2.}], 1 -> 7.]", "{7.0, 2.0}", 0);
    assert_eval_eq("Append[ToNDArray[{1., 2.}], 3.]", "{1.0, 2.0, 3.0}", 0);
    assert_eval_eq("Prepend[ToNDArray[{1., 2.}], 0.]", "{0.0, 1.0, 2.0}", 0);
    assert_eval_eq("ListQ[ToNDArray[{1., 2.}]]", "True", 0);
    assert_eval_eq("VectorQ[ToNDArray[{1., 2.}]]", "True", 0);
    assert_eval_eq("MatrixQ[ToNDArray[{{1., 2.}, {3., 4.}}]]", "True", 0);
}

void test_pattern_matching(void) {
    /* The matcher cannot descend a buffer, so these all rely on the gate
     * running BEFORE DownValues are tried. */
    assert_eval_eq("MatchQ[ToNDArray[{1., 2.}], {__Real}]", "True", 0);
    assert_eval_eq("MatchQ[ToNDArray[{1., 2.}], _List]", "True", 0);
    assert_eval_eq("MatchQ[ToNDArray[{1., 2.}], {_, _}]", "True", 0);
    assert_eval_eq("Module[{f}, f[x_List] := Length[x]; f[ToNDArray[{1., 2., 3.}]]]",
                   "3", 0);
    assert_eval_eq("Module[{g}, g[{a_, b_}] := a + b; g[ToNDArray[{1., 2.}]]]",
                   "3.0", 0);
}

/* The third HPC sweep's fast paths, each as a differential against the identical
 * plain list. Every one of these was added because the packed form was SLOWER
 * than the plain form, or because a small unpacked operand disabled the buffer
 * path for a large packed one -- so "same answer" is the whole acceptance
 * criterion, and the sizes are above PACK_MIN_ELEMENTS so the fast paths
 * actually run. */
void test_third_sweep_fast_paths(void) {
    /* Dot with one operand a plain List: the small side gets packed, and the
     * result must equal the fully-plain product exactly. */
    same_as_plain("ToNDArray[Table[N[i + j], {i, 40}, {j, 8}]] . {1., 2., 3., 4., 5., 6., 7., 8.}",
                  "Table[N[i + j], {i, 40}, {j, 8}] . {1., 2., 3., 4., 5., 6., 7., 8.}");
    same_as_plain("ToNDArray[Table[i + j, {i, 40}, {j, 8}]] . {1, 2, 3, 4, 5, 6, 7, 8}",
                  "Table[i + j, {i, 40}, {j, 8}] . {1, 2, 3, 4, 5, 6, 7, 8}");
    /* ... and the exact/symbolic operands it must NOT touch. */
    same_as_plain("ToNDArray[{{1., 2.}, {3., 4.}}] . {1/2, 1/3}",
                  "{{1., 2.}, {3., 4.}} . {1/2, 1/3}");
    same_as_plain("ToNDArray[{{1., 2.}, {3., 4.}}] . {aa, bb}",
                  "{{1., 2.}, {3., 4.}} . {aa, bb}");

    /* Total with a level SPEC, which used to fall through to the List path. */
    static const char* const TOT[] = {
        "Total[%s]", "Total[%s, {1}]", "Total[%s, {2}]", "Total[%s, {1,2}]",
        "Total[%s, {2,Infinity}]", "Total[%s, 2]", "Total[%s, {0}]",
        "Total[%s, {3}]", "Total[%s, {-1}]", "Head[Total[%s, {2}][[1]]]",
    };
    char bp[320], bl[320];
    for (size_t i = 0; i < sizeof(TOT) / sizeof(TOT[0]); i++) {
        snprintf(bp, sizeof(bp), TOT[i], "ToNDArray[Table[i*10 + j, {i, 30}, {j, 20}]]");
        snprintf(bl, sizeof(bl), TOT[i], "Table[i*10 + j, {i, 30}, {j, 20}]");
        same_as_plain(bp, bl);
    }

    /* Outer over two real vectors, and the forms nd_outer2 declines. */
    static const char* const OUT[] = {
        "Outer[Subtract, %s, %s]", "Outer[Plus, %s, %s]", "Outer[Times, %s, %s]",
        "Outer[Min, %s, %s]", "Outer[Max, %s, %s]",
        "Outer[Divide, %s, %s]", "Outer[List, %s, %s]", "Outer[ff, %s, %s]",
        "Dimensions[Outer[Times, %s, %s]]",
    };
    for (size_t i = 0; i < sizeof(OUT) / sizeof(OUT[0]); i++) {
        snprintf(bp, sizeof(bp), OUT[i], "ToNDArray[Range[1., 30.]]",
                 "ToNDArray[Range[1., 20.]]");
        snprintf(bl, sizeof(bl), OUT[i], "Range[1., 30.]", "Range[1., 20.]");
        same_as_plain(bp, bl);
    }
    /* Integer operands keep the tree path: exactness, not speed, decides. */
    same_as_plain("Outer[Times, ToNDArray[Range[30]], ToNDArray[Range[20]]]",
                  "Outer[Times, Range[30], Range[20]]");

    /* Leading-axis broadcast: rank-2 against rank-1, both orders, both packed
     * and with the smaller side left as a plain List. */
    static const char* const BC[] = {
        "%s + {10., 20., 30.}", "%s - {10., 20., 30.}", "%s * {10., 20., 30.}",
        "{10., 20., 30.} + %s", "Total[(%s - {10., 20., 30.})^2]",
        "%s + {10, 20, 30}", "%s ^ {2., 3., 2.}",
        "%s + {10., 20.}",                       /* length mismatch: stays put */
    };
    for (size_t i = 0; i < sizeof(BC) / sizeof(BC[0]); i++) {
        snprintf(bp, sizeof(bp), BC[i], "ToNDArray[Table[N[i*100 + j], {i, 3}, {j, 100}]]");
        snprintf(bl, sizeof(bl), BC[i], "Table[N[i*100 + j], {i, 3}, {j, 100}]");
        same_as_plain(bp, bl);
    }

    /* A Listable head with one packed and one plain operand of the SAME shape:
     * the List is lifted rather than the buffer materialised. Includes the
     * symbolic case, where lifting must decline and threading must still win. */
    same_as_plain("ToNDArray[Range[1., 300.]] + Normal[Range[1., 300.]]",
                  "Normal[Range[1., 300.]] + Normal[Range[1., 300.]]");
    same_as_plain("ToNDArray[Range[1., 300.]] * Normal[Range[1., 300.]]",
                  "Normal[Range[1., 300.]] * Normal[Range[1., 300.]]");
    same_as_plain("Total[ToNDArray[Range[1., 300.]] + Normal[Range[1., 300.]] + zz]",
                  "Total[Normal[Range[1., 300.]] + Normal[Range[1., 300.]] + zz]");
    same_as_plain("ToNDArray[Range[300]] + Normal[Range[300]]",
                  "Normal[Range[300]] + Normal[Range[300]]");
    same_as_plain("ToNDArray[Range[300]] + Table[1/2, {300}]",
                  "Normal[Range[300]] + Table[1/2, {300}]");

    /* MapThread's column fold, and the heads it declines. */
    static const char* const MT[] = {
        "MapThread[Min, %s]", "MapThread[Max, %s]", "MapThread[Plus, %s]",
        "MapThread[Times, %s]", "MapThread[List, %s]", "MapThread[Subtract, %s]",
        "MapThread[ff, %s]", "MapThread[Plus, %s, 1]",
    };
    for (size_t i = 0; i < sizeof(MT) / sizeof(MT[0]); i++) {
        snprintf(bp, sizeof(bp), MT[i], "ToNDArray[Table[N[Mod[i*7 + j, 11]], {i, 4}, {j, 80}]]");
        snprintf(bl, sizeof(bl), MT[i], "Table[N[Mod[i*7 + j, 11]], {i, 4}, {j, 80}]");
        same_as_plain(bp, bl);
    }
    same_as_plain("MapThread[Plus, ToNDArray[Table[Mod[i*7 + j, 11], {i, 4}, {j, 80}]]]",
                  "MapThread[Plus, Table[Mod[i*7 + j, 11], {i, 4}, {j, 80}]]");

    /* A List of packed rows absorbs into one array -- ragged and mixed
     * exact/inexact must decline, and every one must equal the plain form. */
    same_as_plain("{Range[1., 300.], Range[1., 300.]}",
                  "{Normal[Range[1., 300.]], Normal[Range[1., 300.]]}");
    same_as_plain("Total[{Range[1., 300.], Range[1., 300.]}, {2}]",
                  "Total[{Normal[Range[1., 300.]], Normal[Range[1., 300.]]}, {2}]");
    same_as_plain("{Range[1., 300.], Range[1., 299.]}",
                  "{Normal[Range[1., 300.]], Normal[Range[1., 299.]]}");
    same_as_plain("{Range[300], Range[1., 300.]}",
                  "{Normal[Range[300]], Normal[Range[1., 300.]]}");
    same_as_plain("{Range[300], Range[300]}",
                  "{Normal[Range[300]], Normal[Range[300]]}");
    same_as_plain("Head[{Range[300], Range[300]}[[1, 1]]]",
                  "Head[{Normal[Range[300]], Normal[Range[300]]}[[1, 1]]]");
}

/* Fourth sweep (2026-07-31) -- the structural, scan and convolution paths.
 *
 * Same acceptance criterion as the third sweep above: the packed form must give
 * the byte-identical answer the plain form does, at sizes ABOVE
 * PACK_MIN_ELEMENTS so the fast paths actually run. (A test written below the
 * threshold tests nothing -- two stale claims survived the last sweep exactly
 * that way.) Each group names the defect it guards.
 */
void test_fourth_sweep_fast_paths(void) {
    /* First / Last / Most / Rest read the buffer instead of materialising it.
     * First and Last were asymptotically wrong: an O(1) element read cost O(n). */
    static const char* const HEADTAIL[] = {
        "First[%s]", "Last[%s]", "Most[%s]", "Rest[%s]",
        "Head[First[%s]]", "Head[Last[%s]]",
        "Length[Rest[%s]]", "Length[Most[%s]]",
        "First[%s, zz]", "Last[%s, zz]",
    };
    static const char* const HT_SRC[][2] = {
        /* real rank 1, int rank 1 (exactness), rank 2 (the axis must DROP) */
        { "ToNDArray[Table[N[i], {i, 300}]]",          "Table[N[i], {i, 300}]" },
        { "ToNDArray[Table[i, {i, 300}]]",             "Table[i, {i, 300}]" },
        { "ToNDArray[Table[N[i j], {i, 60}, {j, 5}]]", "Table[N[i j], {i, 60}, {j, 5}]" },
        { "ToNDArray[Table[i j, {i, 60}, {j, 5}]]",    "Table[i j, {i, 60}, {j, 5}]" },
    };
    for (size_t e = 0; e < sizeof(HEADTAIL) / sizeof(HEADTAIL[0]); e++)
        for (size_t k = 0; k < sizeof(HT_SRC) / sizeof(HT_SRC[0]); k++) {
            char pk[512], pl[512];
            snprintf(pk, sizeof(pk), HEADTAIL[e], HT_SRC[k][0]);
            snprintf(pl, sizeof(pl), HEADTAIL[e], HT_SRC[k][1]);
            same_as_plain(pk, pl);
        }
    /* Rest/Most of a SINGLE row is {}, which no buffer shape holds -- it must
     * degrade to the List path rather than answer with something else. */
    same_as_plain("Rest[ToNDArray[Table[N[i], {i, 300}]][[1 ;; 1]]]", "Rest[{1.}]");
    same_as_plain("Most[ToNDArray[Table[N[i], {i, 300}]][[1 ;; 1]]]", "Most[{1.}]");

    /* Clip is back on the buffer, gated on the BOUNDS: Real bounds are uniform,
     * an exact bound is safe only when nothing is clipped. */
    static const char* const CLIPS[] = {
        "Clip[%s, {1.2, 1.8}]",     /* Real bounds -> buffer */
        "Clip[%s, {1, 2}]",         /* exact bounds, clipping -> exact Integers */
        "Clip[%s, {-99, 99}]",      /* exact bounds, nothing clipped -> input */
        "Clip[%s]",                 /* 1-arg: default bounds are exact */
        "Clip[%s, {1.2, 1.8}, {0., 9.}]",   /* 3-arg replacement form */
        "Head[First[Clip[%s, {1, 2}]]]",
        "Head[Last[Clip[%s, {1, 2}]]]",
    };
    static const char* const CLIP_SRC[][2] = {
        { "ToNDArray[Table[N[i]/100., {i, 300}]]", "Table[N[i]/100., {i, 300}]" },
        { "ToNDArray[Table[i, {i, 300}]]",         "Table[i, {i, 300}]" },
    };
    for (size_t e = 0; e < sizeof(CLIPS) / sizeof(CLIPS[0]); e++)
        for (size_t k = 0; k < sizeof(CLIP_SRC) / sizeof(CLIP_SRC[0]); k++) {
            char pk[512], pl[512];
            snprintf(pk, sizeof(pk), CLIPS[e], CLIP_SRC[k][0]);
            snprintf(pl, sizeof(pl), CLIPS[e], CLIP_SRC[k][1]);
            same_as_plain(pk, pl);
        }

    /* Scans. Both spellings of the operator must agree -- and BOTH must agree
     * with the plain list, which is what separates "fast" from "right". */
    static const char* const SCANS[] = {
        "FoldList[Max, 0., %s]",  "FoldList[Min, 9999., %s]",
        "FoldList[Plus, 0., %s]", "FoldList[Times, 1., %s]",
        "Fold[Max, 0., %s]",      "Fold[Plus, 0., %s]",
        "FoldList[Max, %s]",      "FoldList[Plus, %s]",          /* seedless */
        "FoldList[Max[#1, #2] &, 0., %s]",                        /* pure-function */
        "FoldList[#1 + #2 &, 0., %s]",
        "FoldList[Function[{p, q}, p + q^2], 0., %s]",            /* general body */
        "FoldList[Subtract, 0., %s]",                             /* NOT a scan op */
        "FoldList[ff, 0., %s]",                                   /* symbolic f */
        "Head[First[FoldList[Max, 0, %s]]]",
        "Length[FoldList[Max, 0., %s]]",
    };
    static const char* const SCAN_SRC[][2] = {
        { "ToNDArray[Table[N[Mod[7 i, 13]], {i, 300}]]", "Table[N[Mod[7 i, 13]], {i, 300}]" },
        { "ToNDArray[Table[Mod[7 i, 13], {i, 300}]]",    "Table[Mod[7 i, 13], {i, 300}]" },
    };
    for (size_t e = 0; e < sizeof(SCANS) / sizeof(SCANS[0]); e++)
        for (size_t k = 0; k < sizeof(SCAN_SRC) / sizeof(SCAN_SRC[0]); k++) {
            char pk[512], pl[512];
            snprintf(pk, sizeof(pk), SCANS[e], SCAN_SRC[k][0]);
            snprintf(pl, sizeof(pl), SCANS[e], SCAN_SRC[k][1]);
            same_as_plain(pk, pl);
        }
    /* An exact seed beside a Real buffer: Max hands back one of its ARGUMENTS,
     * so the exact 0 can survive into the answer and the scan must decline. */
    same_as_plain("FoldList[Max, 0, ToNDArray[Table[N[i]/1000., {i, 300}]]]",
                  "FoldList[Max, 0, Table[N[i]/1000., {i, 300}]]");
    /* Times over an int64 buffer overflows int64; the List answer is a bigint,
     * which no buffer holds, so the whole scan is abandoned. */
    same_as_plain("FoldList[Times, 1, ToNDArray[Table[10^3, {i, 300}]]]",
                  "FoldList[Times, 1, Table[10^3, {i, 300}]]");

    /* PRESENTATION PARITY. A visible NDArray[...] argument gives a visible
     * NDArray[...] result, whichever internal path answered. The compiled-VM
     * scan lost this and only a cross-spelling comparison caught it: the SAME
     * FoldList answered with head NDArray for Plus and head List for the
     * equivalent lambda. */
    same_as_plain("Head[FoldList[Plus, 0., NDArray[Table[N[i], {i, 300}]]]]",
                  "Head[FoldList[Function[{p, q}, p + q], 0., NDArray[Table[N[i], {i, 300}]]]]");
    same_as_plain("Head[Map[Sqrt, NDArray[Table[N[i], {i, 300}]]]]",
                  "Head[FoldList[Function[{p, q}, p + q^2], 0., NDArray[Table[N[i], {i, 300}]]]]");

    /* ListConvolve / ListCorrelate on the buffer, including every form that
     * must degrade: exact data, symbolic kernel, custom g/h. */
    static const char* const CONVS[] = {
        "ListCorrelate[%s, %s]", "ListConvolve[%s, %s]",
        "ListCorrelate[%s, %s, 1]", "ListCorrelate[%s, %s, -1]",
        "ListCorrelate[%s, %s, {-1, 1}]",
        "ListCorrelate[%s, %s, {1, -1}, 0.]",
        "ListCorrelate[%s, %s, {-1, 1}, {}]",
        "ListCorrelate[%s, %s, {-1, 1}, {}, Times, Plus]",
        "ListCorrelate[%s, %s, {-1, 1}, {}, List, Plus]",
    };
    static const char* const CV_K[2] = { "ToNDArray[{1., 2., 3.}]", "{1., 2., 3.}" };
    static const char* const CV_L[2] = { "ToNDArray[Table[N[Mod[7 i, 11]], {i, 300}]]",
                                         "Table[N[Mod[7 i, 11]], {i, 300}]" };
    for (size_t e = 0; e < sizeof(CONVS) / sizeof(CONVS[0]); e++) {
        char pk[512], pl[512];
        snprintf(pk, sizeof(pk), CONVS[e], CV_K[0], CV_L[0]);
        snprintf(pl, sizeof(pl), CONVS[e], CV_K[1], CV_L[1]);
        same_as_plain(pk, pl);
        /* ...and each side packed ALONE: a mixed pair must not change the answer. */
        snprintf(pk, sizeof(pk), CONVS[e], CV_K[0], CV_L[1]);
        same_as_plain(pk, pl);
        snprintf(pk, sizeof(pk), CONVS[e], CV_K[1], CV_L[0]);
        same_as_plain(pk, pl);
    }
    /* Exact data stays exact -- an int64 buffer must NOT come back as Reals. */
    same_as_plain("ListCorrelate[ToNDArray[{1, 2, 3}], ToNDArray[Table[Mod[7 i, 11], {i, 300}]]]",
                  "ListCorrelate[{1, 2, 3}, Table[Mod[7 i, 11], {i, 300}]]");
    same_as_plain("Head[First[ListCorrelate[ToNDArray[{1, 2, 3}], ToNDArray[Table[Mod[7 i, 11], {i, 300}]]]]]",
                  "Head[First[ListCorrelate[{1, 2, 3}, Table[Mod[7 i, 11], {i, 300}]]]]");
    /* Symbolic kernel: no buffer anywhere in the answer. */
    same_as_plain("ListCorrelate[{aa, bb}, ToNDArray[Table[N[Mod[7 i, 11]], {i, 300}]]]",
                  "ListCorrelate[{aa, bb}, Table[N[Mod[7 i, 11]], {i, 300}]]");
    /* Rank 2, both engines: small takes the direct loop, large takes the FFT. */
    same_as_plain("ListCorrelate[ToNDArray[Table[N[Mod[i + j, 3]], {i, 3}, {j, 3}]], "
                  "ToNDArray[Table[N[Mod[i j, 7]], {i, 40}, {j, 40}]]]",
                  "ListCorrelate[Table[N[Mod[i + j, 3]], {i, 3}, {j, 3}], "
                  "Table[N[Mod[i j, 7]], {i, 40}, {j, 40}]]");
    same_as_plain("ListCorrelate[ToNDArray[Table[N[Mod[i, 5]], {i, 300}]], "
                  "ToNDArray[Table[N[Mod[3 i, 13]], {i, 2000}]]]",
                  "ListCorrelate[Table[N[Mod[i, 5]], {i, 300}], "
                  "Table[N[Mod[3 i, 13]], {i, 2000}]]");

    /* Accumulate / Differences float64 arms, and the exact arms beside them. */
    static const char* const REDUCE[] = {
        "Accumulate[%s]", "Differences[%s]",
        "Head[First[Accumulate[%s]]]", "Head[First[Differences[%s]]]",
    };
    for (size_t e = 0; e < sizeof(REDUCE) / sizeof(REDUCE[0]); e++)
        for (size_t k = 0; k < sizeof(SCAN_SRC) / sizeof(SCAN_SRC[0]); k++) {
            char pk[512], pl[512];
            snprintf(pk, sizeof(pk), REDUCE[e], SCAN_SRC[k][0]);
            snprintf(pl, sizeof(pl), REDUCE[e], SCAN_SRC[k][1]);
            same_as_plain(pk, pl);
        }
    same_as_plain("Accumulate[ToNDArray[Table[N[i j], {i, 60}, {j, 5}]]]",
                  "Accumulate[Table[N[i j], {i, 60}, {j, 5}]]");

    /* Part: the contiguous-span block copy, and every selector that is NOT
     * contiguous and must keep the general gather (steps, reversals, position
     * lists, out-of-order position lists). */
    static const char* const PARTS[] = {
        "%s[[2 ;; -1]]", "%s[[1 ;; 100]]", "%s[[-50 ;; -1]]",
        "%s[[1 ;; -1 ;; 3]]", "%s[[300 ;; 1 ;; -1]]",
        "%s[[{1, 5, 9}]]", "%s[[{9, 1, 5}]]", "%s[[{3, 3, 3}]]",
    };
    for (size_t e = 0; e < sizeof(PARTS) / sizeof(PARTS[0]); e++)
        for (size_t k = 0; k < sizeof(SCAN_SRC) / sizeof(SCAN_SRC[0]); k++) {
            char pk[512], pl[512];
            snprintf(pk, sizeof(pk), PARTS[e], SCAN_SRC[k][0]);
            snprintf(pl, sizeof(pl), PARTS[e], SCAN_SRC[k][1]);
            same_as_plain(pk, pl);
        }
    static const char* const PARTS2[] = {
        "%s[[3 ;; 7, 2 ;; 4]]", "%s[[All, 3 ;; 5]]", "%s[[2 ;; 4]]",
        "%s[[All, 4]]", "%s[[5, All]]", "%s[[{1, 3, 5}, All]]",
        "%s[[All, {2, 4}]]", "%s[[{3, 1, 2}, All]]",
    };
    for (size_t e = 0; e < sizeof(PARTS2) / sizeof(PARTS2[0]); e++) {
        char pk[512], pl[512];
        snprintf(pk, sizeof(pk), PARTS2[e], "ToNDArray[Table[N[i j], {i, 60}, {j, 5}]]");
        snprintf(pl, sizeof(pl), PARTS2[e], "Table[N[i j], {i, 60}, {j, 5}]");
        same_as_plain(pk, pl);
    }

    /* Outer with ONE operand below the packing threshold -- the small-operand
     * trap again, this time costing 849 ms on a 100000 x 16 product. */
    same_as_plain("Outer[Plus, ToNDArray[Table[N[i], {i, 300}]], {1., 2., 3.}]",
                  "Outer[Plus, Table[N[i], {i, 300}], {1., 2., 3.}]");
    same_as_plain("Outer[Times, {1., 2., 3.}, ToNDArray[Table[N[i], {i, 300}]]]",
                  "Outer[Times, {1., 2., 3.}, Table[N[i], {i, 300}]]");
    /* ...and the operands it must still decline. */
    same_as_plain("Outer[Plus, ToNDArray[Table[N[i], {i, 300}]], {1, 2, 3}]",
                  "Outer[Plus, Table[N[i], {i, 300}], {1, 2, 3}]");
    same_as_plain("Outer[Plus, ToNDArray[Table[N[i], {i, 300}]], {aa, bb}]",
                  "Outer[Plus, Table[N[i], {i, 300}], {aa, bb}]");
    same_as_plain("Outer[ff, ToNDArray[Table[N[i], {i, 300}]], {1., 2.}]",
                  "Outer[ff, Table[N[i], {i, 300}], {1., 2.}]");
}

/* ---------------------------------------------------------------------------
 *  Fifth sweep (2026-07-31): the packed index, the integer scan and thread,
 *  Ramp, the infinite Clip bound, and the packed set operations.
 *
 *  Every case below is the SAME source evaluated twice -- once with automatic
 *  packing on, once with it off -- so what is compared is the new buffer path
 *  against the interpreter's List path on identical input. That is stronger
 *  than comparing against a written-down expectation, which can be wrong in the
 *  same direction as the code, and it is the form the third and fourth sweeps
 *  arrived at after three separate tests turned out to assert nothing.
 *
 *  Every source is ABOVE PACK_MIN_ELEMENTS. A differential case built on 3
 *  elements never runs the path it claims to test (experiment 10, experiment 11).
 * ------------------------------------------------------------------------- */
static void both_ways(const char* src) {
    Expr* on = ev(src);
    pack_set_enabled(false);
    Expr* off = ev(src);
    pack_set_enabled(true);
    char* sa = expr_to_string(on);
    char* sb = expr_to_string(off);
    if (strcmp(sa, sb) != 0)
        fprintf(stderr, "FAIL: packed and unpacked disagree\n  %s\n"
                        "  packed:   %s\n  unpacked: %s\n", src, sa, sb);
    ASSERT(strcmp(sa, sb) == 0);
    free(sa); free(sb);
    expr_free(on); expr_free(off);
}

void test_fifth_sweep_fast_paths(void) {
    /* ---- 1. Part with a PACKED index list -- the gather. ------------------
     * build_axis_selector accepted a List of Integer positions but not a packed
     * one, and every producer of an index list (Flatten, Range, RandomInteger,
     * arithmetic on any of them) now hands back a packed one -- so x[[idx]]
     * degraded the whole Part and materialised BOTH arrays. */
    static const char* const GATHER[] = {
        "Table[N[i]/7., {i, 300}][[Table[Mod[7 i, 300] + 1, {i, 300}]]]",
        "Table[i, {i, 300}][[Table[Mod[7 i, 300] + 1, {i, 300}]]]",
        "Table[N[i]/7., {i, 300}][[Range[300]]]",
        "Table[N[i]/7., {i, 300}][[-Range[300]]]",
        "Table[N[i], {i, 300}][[Table[1, {300}]]]",
        "Table[N[i j], {i, 60}, {j, 5}][[Table[Mod[3 i, 60] + 1, {i, 300}]]]",
        "Table[N[i j], {i, 60}, {j, 5}][[Table[Mod[3 i, 60] + 1, {i, 300}], 2]]",
        "Table[N[i j], {i, 60}, {j, 5}][[All, Table[Mod[i, 5] + 1, {i, 300}]]]",
        "Table[i j, {i, 60}, {j, 5}][[Table[Mod[3 i, 60] + 1, {i, 300}]]]",
        /* The answer's ELEMENT HEADS, not just its values: an int64 source must
         * still gather exact Integers. */
        "Head[Table[i, {i, 300}][[Range[300]]][[1]]]",
        "Head[Table[N[i], {i, 300}][[Range[300]]][[1]]]",
        "Total[Table[N[i], {i, 300}][[Table[Mod[11 i, 300] + 1, {i, 300}]]]]",
        "Length[Table[N[i], {i, 300}][[Table[Mod[11 i, 300] + 1, {i, 300}]]]]",
    };
    for (size_t i = 0; i < sizeof(GATHER) / sizeof(GATHER[0]); i++)
        both_ways(GATHER[i]);

    /* Forms the gather must DECLINE, each for its own reason. All four have to
     * answer exactly what the List path answers, which for three of them is an
     * error or an unevaluated Part. */
    both_ways("Table[N[i], {i, 300}][[Table[1., {300}]]]");        /* Real positions */
    both_ways("Table[N[i], {i, 300}][[Table[400, {300}]]]");       /* out of range */
    both_ways("Table[N[i], {i, 300}][[Table[-400, {300}]]]");      /* out of range, negative */
    both_ways("Table[N[i], {i, 300}][[{}]]");                      /* empty spec */
    both_ways("Table[N[i], {i, 300}][[Table[Mod[i, 300] + 1, {i, 60}, {j, 5}]]]");  /* rank-2 spec */

    /* ---- 2. MapThread over an INTEGER buffer. ----------------------------
     * nd_mapthread2 was float64-only and MapThread was absent from INT64_OK, so
     * an integer pair materialised twice over. */
    static const char* const MT[] = {
        "MapThread[Max, {Table[Mod[7 i, 97], {i, 300}], Table[Mod[11 i, 89], {i, 300}]}]",
        "MapThread[Min, {Table[Mod[7 i, 97], {i, 300}], Table[Mod[11 i, 89], {i, 300}]}]",
        "MapThread[Plus, {Table[Mod[7 i, 97], {i, 300}], Table[Mod[11 i, 89], {i, 300}]}]",
        "MapThread[Times, {Table[Mod[7 i, 97], {i, 300}], Table[Mod[11 i, 89], {i, 300}]}]",
        "MapThread[Max, {Table[Mod[7 i, 97], {i, 300}], Table[Mod[11 i, 89], {i, 300}],"
        " Table[Mod[13 i, 83], {i, 300}]}]",
        /* Element heads: an integer thread must not come back Real. */
        "Head[MapThread[Max, {Table[Mod[7 i, 97], {i, 300}],"
        " Table[Mod[11 i, 89], {i, 300}]}][[1]]]",
        /* The float64 control, unchanged by this sweep. */
        "MapThread[Max, {Table[N[Mod[7 i, 97]], {i, 300}], Table[N[Mod[11 i, 89]], {i, 300}]}]",
        /* Mixed exactness never packs in the first place, so this exercises the
         * List path on both sides and must still agree. */
        "MapThread[Max, {Table[Mod[7 i, 97], {i, 300}], Table[N[Mod[11 i, 89]], {i, 300}]}]",
        /* OVERFLOW abandons the whole array so GMP answers exactly -- the one
         * behaviour a wrapping int64 loop would get silently wrong. */
        "MapThread[Times, {Table[10^18, {300}], Table[10^18, {300}]}][[1]]",
        "MapThread[Plus, {Table[9223372036854775807, {300}], Table[1, {300}]}][[1]]",
    };
    for (size_t i = 0; i < sizeof(MT) / sizeof(MT[0]); i++) both_ways(MT[i]);

    /* ---- 3. Fold / FoldList over an INTEGER buffer -- the scan. ----------- */
    static const char* const SCAN[] = {
        "FoldList[Max, 0, Table[Mod[7 i, 97], {i, 300}]]",
        "FoldList[Min, 999, Table[Mod[7 i, 97], {i, 300}]]",
        "FoldList[Plus, 0, Table[Mod[7 i, 97], {i, 300}]]",
        "Fold[Max, 0, Table[Mod[7 i, 97], {i, 300}]]",
        "Fold[Plus, 0, Table[Mod[7 i, 97], {i, 300}]]",
        "FoldList[Max, Table[Mod[7 i, 97], {i, 300}]]",           /* seedless */
        "Head[FoldList[Max, 0, Table[Mod[7 i, 97], {i, 300}]][[1]]]",
        "Head[FoldList[Plus, 0, Table[Mod[7 i, 97], {i, 300}]][[-1]]]",
        /* A REAL seed over an integer buffer is the exactness trap: Max hands
         * back one of its arguments, so the answer is Real exactly where the
         * seed won and exact everywhere else -- mixed, and unpackable. */
        "FoldList[Max, 0., Table[Mod[7 i, 97], {i, 300}]]",
        "FoldList[Max, 1000., Table[Mod[7 i, 97], {i, 300}]]",
        /* Overflow, again: the product leaves int64 and must reach GMP. */
        "Fold[Times, 1, Table[If[i <= 5, 10^15, 1], {i, 300}]]",
        /* The pure-function spelling must agree with the bare symbol -- the
         * presentation-parity trap of experiment 11. */
        "FoldList[Max[#1, #2] &, 0, Table[Mod[7 i, 97], {i, 300}]]",
        "Head[FoldList[Max[#1, #2] &, 0, Table[Mod[7 i, 97], {i, 300}]]]",
        "Head[FoldList[Max, 0, Table[Mod[7 i, 97], {i, 300}]]]",
    };
    for (size_t i = 0; i < sizeof(SCAN) / sizeof(SCAN[0]); i++) both_ways(SCAN[i]);

    /* ---- 4. Ramp -- a new builtin, and its buffer kernel. ----------------- */
    /* Against the Mathematica 14.0 reference, read off directly. */
    assert_eval_eq("Ramp[{-1., 0., 2.5}]", "{0.0, 0.0, 2.5}", 0);
    assert_eval_eq("Head /@ Ramp[{-1., 0., 2.5}]", "{Real, Real, Real}", 0);
    assert_eval_eq("Ramp[{-3, 0, 4}]", "{0, 0, 4}", 0);
    assert_eval_eq("Ramp[{-1/2, 3/4}]", "{0, 3/4}", 0);
    assert_eval_eq("Ramp[x]", "Ramp[x]", 0);
    assert_eval_eq("Ramp[1. + 2. I]", "Ramp[1.0 + 2.0*I]", 0);
    assert_eval_eq("Attributes[Ramp]", "{Listable, NumericFunction, Protected}", 0);
    /* An exact symbolic argument whose sign IS decidable resolves; Ramp is not
     * allowed to give up on it just because it is not a machine number. */
    assert_eval_eq("Ramp[Sqrt[2] - 1]", "-1 + Sqrt[2]", 0);
    assert_eval_eq("Ramp[1 - Sqrt[2]]", "0", 0);

    static const char* const RAMP[] = {
        "Ramp[Table[N[i] - 150., {i, 300}]]",
        "Ramp[Table[i - 150, {i, 300}]]",
        "Ramp[Table[N[i j] - 150., {i, 60}, {j, 5}]]",
        "Head[Ramp[Table[N[i] - 150., {i, 300}]][[1]]]",
        "Head[Ramp[Table[i - 150, {i, 300}]][[1]]]",
        "Total[Ramp[Table[N[i] - 150., {i, 300}]]]",
        "Total[Ramp[Table[i - 150, {i, 300}]]]",
        "Ramp[Table[i/7, {i, 300}] - 20]",          /* Rationals never pack */
        "Ramp[Table[If[i == 7, xx, N[i] - 150.], {i, 300}]]",   /* one symbol: no pack */
    };
    for (size_t i = 0; i < sizeof(RAMP) / sizeof(RAMP[0]); i++) both_ways(RAMP[i]);

    /* ---- 5. Clip with an INFINITE bound. ---------------------------------
     * Clip[x, {0., Infinity}] is how the positive part is spelled, and the
     * bound simply failed to parse: every such call returned unevaluated. */
    assert_eval_eq("Clip[{-2., 0.5, 3.}, {0., Infinity}]", "{0.0, 0.5, 3.0}", 0);
    assert_eval_eq("Clip[{-2., 0.5, 3.}, {-Infinity, 1.}]", "{-2.0, 0.5, 1.0}", 0);
    assert_eval_eq("Clip[{-2., 0.5, 3.}, {-Infinity, Infinity}]", "{-2.0, 0.5, 3.0}", 0);
    assert_eval_eq("Clip[{-2, 5}, {0, Infinity}]", "{0, 5}", 0);
    assert_eval_eq("Head /@ Clip[{-2., 0.5, 3.}, {0., Infinity}]", "{Real, Real, Real}", 0);

    static const char* const CLIPINF[] = {
        "Clip[Table[N[i] - 150., {i, 300}], {0., Infinity}]",
        "Clip[Table[N[i] - 150., {i, 300}], {-Infinity, 0.}]",
        "Clip[Table[N[i] - 150., {i, 300}], {-Infinity, Infinity}]",
        /* An EXACT bound beside Real data still has to put an exact head at
         * every clipped position, infinite partner or not. */
        "Clip[Table[N[i] - 150., {i, 300}], {0, Infinity}]",
        "Head[Clip[Table[N[i] - 150., {i, 300}], {0, Infinity}][[1]]]",
        "Head[Clip[Table[N[i] - 150., {i, 300}], {0., Infinity}][[1]]]",
        "Clip[Table[i - 150, {i, 300}], {0, Infinity}]",
        "Total[Clip[Table[N[i] - 150., {i, 300}], {0., Infinity}]]",
    };
    for (size_t i = 0; i < sizeof(CLIPINF) / sizeof(CLIPINF[0]); i++) both_ways(CLIPINF[i]);

    /* ---- 6. The packed set operations. ------------------------------------
     * Union / Intersection / Complement over integer LABELS is what a graph
     * traversal is made of; the generic path allocates one Expr per element and
     * sorts through expr_compare. */
    static const char* const SETOPS[] = {
        "Union[Table[Mod[7 i, 97], {i, 300}]]",
        "Union[Table[Mod[7 i, 97], {i, 300}], Table[Mod[11 i, 89], {i, 300}]]",
        "Union[Table[Mod[7 i, 97], {i, 300}], Table[Mod[11 i, 89], {i, 300}],"
        " Table[Mod[13 i, 83], {i, 300}]]",
        "Intersection[Table[Mod[7 i, 97], {i, 300}], Table[Mod[11 i, 89], {i, 300}]]",
        "Complement[Table[Mod[7 i, 97], {i, 300}], Table[Mod[11 i, 89], {i, 300}]]",
        "Complement[Table[i, {i, 300}], Table[i, {i, 300}]]",           /* empty */
        "Complement[Table[i, {i, 300}], Table[2 i, {i, 300}]]",
        "Intersection[Table[i, {i, 300}], Table[i + 1000, {i, 300}]]",  /* empty */
        "Length[Union[Table[Mod[7 i, 97], {i, 300}]]]",
        "Head[Union[Table[Mod[7 i, 97], {i, 300}]][[1]]]",
        /* Forms that must degrade: Reals (0. and -0. compare equal and print
         * differently, so which representative survives would differ), a
         * SameTest, and a mix of packed and plain operands. */
        "Union[Table[N[Mod[7 i, 97]], {i, 300}]]",
        "Union[Table[Mod[7 i, 97], {i, 300}], SameTest -> (Abs[#1 - #2] < 2 &)]",
        "Union[Table[Mod[7 i, 97], {i, 300}], {1000, 2000}]",
        "Complement[Table[Mod[7 i, 97], {i, 300}], {1, 2, 3}]",
        "Union[Table[Mod[7 i, 97], {i, 300}], {a, b}]",                 /* symbolic */
    };
    for (size_t i = 0; i < sizeof(SETOPS) / sizeof(SETOPS[0]); i++) both_ways(SETOPS[i]);

    /* ---- 6b. Round two: the int64 arms the first pass missed, Join's small
     * operands, and a gather whose SOURCE is below the threshold. Each was found
     * by re-profiling a kernel that had not moved, not by reading the code. */

    /* Abs of an exact integer is an exact integer -- a projection kernel writes
     * a real dtype, so this arm has to exist for the buffer to answer the same
     * thing the List path does. |INT64_MIN| abandons to GMP. */
    assert_eval_eq("Abs[{-3, 2, -1}]", "{3, 2, 1}", 0);
    assert_eval_eq("Head /@ Abs[{-3, 2, -1}]", "{Integer, Integer, Integer}", 0);
    assert_eval_eq("Abs[{-3., 2.}]", "{3.0, 2.0}", 0);
    assert_eval_eq("Abs[3. + 4. I]", "5.0", 0);

    static const char* const R2[] = {
        /* Abs: integer, real, and the overflow that must reach a bignum. */
        "Abs[Table[i - 150, {i, 300}]]",
        "Head[Abs[Table[i - 150, {i, 300}]][[1]]]",
        "Total[Abs[Table[i - 150, {i, 300}]]]",
        "Abs[Table[N[i] - 150., {i, 300}]]",
        "Abs[Table[If[i == 3, -9223372036854775808, -i], {i, 300}]][[3]]",
        "Head[Abs[Table[If[i == 3, -9223372036854775808, -i], {i, 300}]][[1]]]",
        /* First / Last / Most / Rest on an INTEGER buffer. */
        "Most[Table[i, {i, 300}]]", "Rest[Table[i, {i, 300}]]",
        "First[Table[i, {i, 300}]]", "Last[Table[i, {i, 300}]]",
        "Head[Most[Table[i, {i, 300}]][[1]]]", "Head[First[Table[i, {i, 300}]]]",
        "Most[Table[i j, {i, 60}, {j, 5}]]", "First[Table[i j, {i, 60}, {j, 5}]]",
        /* Join with small plain operands beside a buffer -- the explicit
         * finite-difference boundary. The lifted operand must SNIFF to the same
         * dtype, so an exact 1 beside Real data still gives a mixed answer. */
        "Join[{0.}, Table[N[i], {i, 300}], {9.}]",
        "Head[Join[{0.}, Table[N[i], {i, 300}], {9.}][[1]]]",
        "Length[Join[{0.}, Table[N[i], {i, 300}], {9.}]]",
        "Join[{1}, Table[N[i], {i, 300}]]",
        "Head[Join[{1}, Table[N[i], {i, 300}]][[1]]]",
        "Join[{1}, Table[i, {i, 300}]]",
        "Join[Table[N[i], {i, 300}], {a, b}]",              /* symbolic: declines */
        "Join[{{1., 2.}}, Table[{N[i], N[i]}, {i, 300}]]",  /* rank 2 */
        "Join[{{1., 2.}}, {{3., 4.}}, 2]",                  /* the level form */
        "Join[Table[N[i], {i, 300}], Table[i, {i, 300}]]",  /* dtype mismatch */
        /* Partition with an OFFSET -- the sliding window. The default d == k
         * tiles and is one memcpy; any other d overlaps and is `rows` strided
         * copies. Only complete rows, matching the List path; the padded 4-arg
         * form is declined. */
        "Partition[Table[N[i], {i, 300}], 5, 3]",
        "Partition[Table[N[i], {i, 300}], 5]",
        "Partition[Table[N[i], {i, 300}], 7, 4]",
        "Partition[Table[N[i], {i, 300}], 5, 7]",       /* d > k: rows do not overlap */
        "Partition[Table[N[i], {i, 300}], 300, 1]",     /* exactly one row */
        "Partition[Table[N[i], {i, 300}], 301, 1]",     /* shorter than k: {} */
        "Partition[Table[N[i], {i, 300}], 5, 3, 1]",    /* padded: declines */
        "Partition[Table[i, {i, 300}], 5, 3]",
        "Head[Partition[Table[i, {i, 300}], 5, 3][[1, 1]]]",
        "Dimensions[Partition[Table[N[i], {i, 300}], 12, 1]]",
        "Total[Partition[Table[N[i], {i, 300}], 12, 1], 2]",
        "Partition[Table[N[i j], {i, 60}, {j, 5}], 7, 3]",   /* rank 2: declines */
        /* A gather whose SOURCE is below the packing threshold and whose INDEX
         * is above it. Lifting must not change the answer, an out-of-range
         * diagnostic included. */
        "Table[N[i], {i, 65}][[Table[Mod[7 i, 65] + 1, {i, 300}]]]",
        "Total[Table[N[i], {i, 65}][[Table[Mod[7 i, 65] + 1, {i, 300}]]]]",
        "Table[i, {i, 65}][[Table[Mod[7 i, 65] + 1, {i, 300}]]]",
        "Head[Table[i, {i, 65}][[Table[Mod[7 i, 65] + 1, {i, 300}]]][[1]]]",
        "Table[N[i], {i, 65}][[Table[100, {300}]]]",        /* out of range */
        "Table[i/7, {i, 65}][[Table[Mod[7 i, 65] + 1, {i, 300}]]]",  /* Rationals */
        "{a, b, c}[[Table[Mod[i, 3] + 1, {i, 300}]]]",      /* symbolic source */
    };
    for (size_t i = 0; i < sizeof(R2) / sizeof(R2[0]); i++) both_ways(R2[i]);

    /* ---- 7. The compositions the sweep was actually built out of. ---------
     * A fast path can be right in isolation and wrong where it meets the next
     * one; these are the four inner loops from experiments 12, 14 and 19. */
    both_ways("Module[{p = Table[1./300., {300}], f = Table[Mod[7 i, 300] + 1, {i, 1200}]},"
              " Total[Total[Partition[p[[f]], 4], {2}]]]");                  /* SpMV */
    both_ways("Module[{v = Table[Mod[3 i, 61], {i, 300}], w = Table[Mod[5 i, 47], {i, 300}]},"
              " Total[FoldList[Max, 0, MapThread[Max, {v, w}] + Range[300]]]]");  /* NW row */
    both_ways("Module[{a = Table[Mod[7 i, 97] + 1, {i, 300}]},"
              " Length[Complement[Union[a], Table[i, {i, 50}]]]]");          /* BFS level */
    both_ways("Module[{t = Table[N[i] - 150., {i, 300}], c = Table[Mod[i, 40] + 1, {i, 300}]},"
              " Total[Ramp[t] Table[N[i]/7., {i, 40}][[c]]]]");              /* shade + gather */
    both_ways("Module[{v = Table[N[i]/300., {i, 300}], p = Table[N[i]/600., {i, 300}], k = 0},"
              " While[k < 5, v = Join[{0.}, 0.25 Most[Most[v]] + 0.5 Take[v, {2, -2}]"
              " + 0.25 Rest[Rest[v]], {1.}]; v = MapThread[Max, {v, p}]; k = k + 1];"
              " Total[v]]");                                                  /* explicit FD sweep */
    both_ways("Module[{s = Table[Mod[7 i, 4], {i, 300}], c = 2},"
              " Total[2 - 3 UnitStep[Abs[s - c] - 1]]]");                      /* NW substitution */
    both_ways("Module[{c = Table[Mod[i^2 + 3 i, 4], {i, 300}], p = Table[4^(6 - j), {j, 6}]},"
              " Length[Union[Partition[c, 6, 1] . p]]]");                      /* k-mer count */
}

/* The seventh round: the coverage register's fast paths (2026-08-01).
 *
 * What makes this set worth its own function is that FOUR of the heads in it
 * were already on pack.c's AWARE list and were slow anyway -- Extract,
 * MatrixPower, PseudoInverse and LeastSquares each declined the buffer on their
 * own first line. No static check sees that, so the only guard is a
 * differential one. */
void test_seventh_round_fast_paths(void) {
    /* ---- 1. Extract reads the buffer. -----------------------------------
     * expr_part cannot index an NDArray (is_atomic is true for one), so Extract
     * returned NULL and the POST-gate materialised on the way to rest. The
     * declines matter as much as the hits: an out-of-range position and [[0]]
     * must answer exactly what the List path answers. */
    static const char* const EXTRACT[] = {
        "Extract[Table[N[i], {i, 300}], {{1}, {2}, {300}}]",
        "Extract[Table[i, {i, 300}], {{1}, {2}, {300}}]",
        "Head[Extract[Table[i, {i, 300}], {{7}}][[1]]]",       /* exact Integer out */
        "Head[Extract[Table[N[i], {i, 300}], {{7}}][[1]]]",
        "Extract[Table[N[i], {i, 300}], 5]",                   /* scalar position */
        "Extract[Table[N[i], {i, 300}], {5}]",
        "Extract[Table[N[i], {i, 300}], {-1}]",
        "Extract[Table[N[i j], {i, 60}, {j, 5}], {{3, 4}, {7, 1}}]",
        "Extract[Table[N[i], {i, 300}], {0}]",                 /* head extraction */
        "Extract[Table[N[i], {i, 300}], {{999}}]",             /* out of range */
        "Extract[Table[N[i], {i, 300}], {{1}}, Hold]",
    };
    for (size_t i = 0; i < sizeof(EXTRACT) / sizeof(EXTRACT[0]); i++)
        both_ways(EXTRACT[i]);

    /* ---- 2. MatrixPower keeps the buffer. --------------------------------
     * It used to delist and run square-and-multiply over boxed matrices. An
     * INTEGER matrix must still materialise and stay exact -- MatrixPower is
     * deliberately not on INT64_OK, because the exact fourth power of an
     * integer matrix does not fit a float64 buffer past 2^53. */
    static const char* const MATPOW[] = {
        "Total[Flatten[MatrixPower[Table[N[1/(i + j)], {i, 20}, {j, 20}], 3]]]",
        "Head[MatrixPower[Table[i + j, {i, 20}, {j, 20}], 3][[1, 1]]]",
        "MatrixPower[Table[i + j, {i, 20}, {j, 20}], 0][[1, 1]]",
        "Head[MatrixPower[Table[N[i + j], {i, 20}, {j, 20}], 2][[1, 1]]]",
        "MatrixPower[Table[N[i + j], {i, 20}, {j, 20}], nn]",  /* symbolic: unevaluated */
    };
    for (size_t i = 0; i < sizeof(MATPOW) / sizeof(MATPOW[0]); i++)
        both_ways(MATPOW[i]);

    /* ---- 3. The SVD pair. -------------------------------------------------
     * PseudoInverse and LeastSquares on a machine matrix take one gesdd, and on
     * an EXACT one keep the rationalised pipeline. Only the second is testable
     * both ways bit-for-bit: the machine path is a different algorithm from the
     * exact one, so the assertions on it are the Moore-Penrose identities. */
    both_ways("PseudoInverse[Table[If[i == j, 2, 1], {i, 20}, {j, 20}]][[1, 1]]");
    both_ways("Head[PseudoInverse[Table[If[i == j, 2, 1], {i, 20}, {j, 20}]][[1, 1]]]");
    assert_eval_eq("Module[{a = Table[N[Mod[7 i + 3 j, 11]] + If[i == j, 30., 0.],"
                   "                    {i, 20}, {j, 20}]},"
                   " Max[Abs[Flatten[a . PseudoInverse[a] . a - a]]] < 1.*^-9]", "True", 0);
    assert_eval_eq("Module[{a = Table[N[Mod[7 i + 3 j, 11]], {i, 30}, {j, 12}],"
                   "        b = Table[N[Mod[5 i, 7]], {i, 30}]},"
                   " Max[Abs[Transpose[a] . (a . LeastSquares[a, b] - b)]] < 1.*^-9]", "True", 0);
    /* An exact system stays exact and is unaffected by the SVD path. */
    assert_eval_eq("LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {1, 2, 3}]", "{0, 1}", 0);

    /* ---- 4. The producers. ------------------------------------------------
     * Each writes a buffer only where every element it would have boxed has the
     * same head, so the mixed cases below must come back UNPACKED and identical
     * to the List path -- that is the whole of the rule. */
    static const char* const PRODUCE[] = {
        "Subdivide[0., 1., 300]",
        "Subdivide[3., 1., 300]",              /* descending */
        "Subdivide[5., 5., 300]",              /* degenerate */
        "Subdivide[0, 1., 300]",               /* one Real endpoint: all Real */
        "Subdivide[3., 1, 300]",
        "Subdivide[1/2, 1., 300]",             /* the Rational numericalises */
        "Subdivide[10, 300]",                  /* exact: Rationals, no pack */
        "Subdivide[a, b, 4]",                  /* symbolic */
        "Head[Subdivide[0., 1., 300][[1]]]",
        "Head[Subdivide[0, 1., 300][[1]]]",
        "IdentityMatrix[20]",
        "Head[IdentityMatrix[20][[1, 1]]]",
        "IdentityMatrix[{3, 5}]",
        "DiagonalMatrix[Table[i, {i, 20}]]",
        "DiagonalMatrix[Table[N[i], {i, 20}]]",
        "Head[DiagonalMatrix[Table[N[i], {i, 20}]][[1, 2]]]",
        "DiagonalMatrix[Table[If[i == 1, 1., i], {i, 20}]]",   /* one Real widens all */
        "Head[DiagonalMatrix[Table[If[i == 1, 1., i], {i, 20}]][[5, 5]]]",
        "DiagonalMatrix[Table[i, {i, 20}], 2]",
        "DiagonalMatrix[Table[i, {i, 20}], -2]",
        "UnitVector[300, 7]",
        "Head[UnitVector[300, 7][[1]]]",
        "UnitVector[300, 7, WorkingPrecision -> MachinePrecision]",
        "Rescale[Table[N[i], {i, 300}]]",
        "Rescale[Table[i, {i, 300}]]",            /* exact: Rationals */
        "Rescale[Table[N[i], {i, 300}], {0., 10.}]",
        "Rescale[Table[N[i], {i, 300}], {0., 10.}, {1., 2.}]",
        "Rescale[Table[i, {i, 20}, {j, 3}]]",     /* nested threading */
        "Head[Rescale[Table[N[i], {i, 300}]][[1]]]",
    };
    for (size_t i = 0; i < sizeof(PRODUCE) / sizeof(PRODUCE[0]); i++)
        both_ways(PRODUCE[i]);

    /* MACHINE-REAL CONTAGION. One machine Real makes the WHOLE result Real --
     * the zeros a producer invents and the exact entries it was given alike.
     * Every expected value below is Mathematica 14.0's, checked against it
     * directly; each was WRONG here before 2026-08-01, writing an exact Integer
     * 0 (or keeping an exact endpoint) and thereby making the result
     * two-headed, which is what stopped it packing. Only MACHINE Real is
     * contagious: the MPFR and exact rows must keep their exact zeros. */
    assert_eval_eq("DiagonalMatrix[{1., 2., 3.}]",
                   "{{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}}", 0);
    assert_eval_eq("DiagonalMatrix[{1, 2, 3.}]",
                   "{{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}}", 0);
    assert_eval_eq("DiagonalMatrix[{1/2, 1.}]", "{{0.5, 0.0}, {0.0, 1.0}}", 0);
    assert_eval_eq("DiagonalMatrix[{a, 1.}]", "{{a, 0.0}, {0.0, 1.0}}", 0);
    assert_eval_eq("DiagonalMatrix[{1, 2, 3}]",
                   "{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}", 0);
    assert_eval_eq("DiagonalMatrix[{1/2, 3/4}]", "{{1/2, 0}, {0, 3/4}}", 0);
    assert_eval_eq("DiagonalMatrix[{a, b}]", "{{a, 0}, {0, b}}", 0);
    assert_eval_eq("Subdivide[0, 1., 4]", "{0.0, 0.25, 0.5, 0.75, 1.0}", 0);
    assert_eval_eq("Subdivide[3., 1, 4]", "{3.0, 2.5, 2.0, 1.5, 1.0}", 0);
    assert_eval_eq("Subdivide[1/2, 1., 4]", "{0.5, 0.625, 0.75, 0.875, 1.0}", 0);
    assert_eval_eq("Subdivide[10, 4]", "{0, 5/2, 5, 15/2, 10}", 0);
    /* An MPFR endpoint is NOT machine-real, so the exact 0 survives. */
    assert_eval_eq("Head[Subdivide[0, 1.`30, 4][[1]]]", "Integer", 0);
    assert_eval_eq("Head[DiagonalMatrix[{1.`30, 2}][[1, 2]]]", "Integer", 0);
    /* And it is the exactness, not the packing, that decides: a Real diagonal
     * now packs precisely because it is one dtype. */
    assert_eval_eq("NDArrayQ[DiagonalMatrix[Table[N[i], {i, 300}]]]", "True", 0);
    assert_eval_eq("NDArrayQ[Subdivide[0, 1., 300]]", "True", 0);

    /* Subdivide's interior points are min + i*step, which is what both
     * Mathematica and numpy.linspace compute -- checked bit-for-bit against
     * linspace. The fourth point of {0., 1., 10} is the discriminating one. */
    /* Not spelled with ==: compare_numeric treats two inexact operands as equal
     * when they agree to a relative 2^-46, so 0.30000000000000004 == 0.3 is
     * True. Scaling to an exact Integer is the only way to see the last bit. */
    assert_eval_eq("Round[10^17 Subdivide[0., 1., 10][[4]]]", "30000000000000004", 0);
    assert_eval_eq("Round[10^17 Subdivide[0., 1., 10][[7]]]", "60000000000000008", 0);
    assert_eval_eq("Subdivide[0., 1., 10][[1]] === 0. && Subdivide[0., 1., 10][[11]] === 1.",
                   "True", 0);

    /* ---- 5. ConjugateTranspose. -------------------------------------------
     * Its NDArray path handed a rank-1 buffer to Transpose, which declines one,
     * so the composite came back unevaluated -- which is why the head was
     * EXEMPT from the packed-aware audit rather than aware. */
    static const char* const CT[] = {
        "ConjugateTranspose[Table[N[i j], {i, 20}, {j, 20}]] == "
        "Transpose[Table[N[i j], {i, 20}, {j, 20}]]",
        "ConjugateTranspose[Table[N[i], {i, 300}]]",       /* rank 1: shape kept */
        "ConjugateTranspose[Table[i, {i, 300}]]",
        "Head[ConjugateTranspose[Table[i, {i, 20}, {j, 20}]][[1, 1]]]",
        "ConjugateTranspose[Table[i j, {i, 20}, {j, 20}], {2, 1}]",
        "ConjugateTranspose[Table[N[i j], {i, 20}, {j, 20}], bad]",   /* bad spec */
    };
    for (size_t i = 0; i < sizeof(CT) / sizeof(CT[0]); i++) both_ways(CT[i]);

    /* ---- 6. ProductLog's machine kernel was WRONG just above x = 1. -------
     * The asymptotic seed log x - log log x + log log x / log x was used for
     * every x >= 1; it is log(0)/0 at x = 1 and diverges above it, and Halley
     * cannot recover. ProductLog[1.01] returned -338.392 for an answer of
     * 0.5707. Only the compiled and packed paths used that kernel, so the
     * assertion has to name them: the interpreter went to MPFR and was right. */
    assert_eval_eq("Module[{f = Compile[{x}, ProductLog[x]]},"
                   " Abs[f[1.01] - 0.570752468] < 1.*^-6]", "True", 0);
    static const char* const PLOG[] = {
        "Total[ProductLog[Table[1. + N[i]/1000., {i, 300}]]]",
        "Max[Abs[Module[{w = ProductLog[Table[1. + N[i]/1000., {i, 300}]]},"
        " w Exp[w] - Table[1. + N[i]/1000., {i, 300}]]]] < 1.*^-12",
        "Total[ProductLog[Table[N[i]/100., {i, 300}]]]",
    };
    for (size_t i = 0; i < sizeof(PLOG) / sizeof(PLOG[0]); i++) both_ways(PLOG[i]);
}

/* NO TWO-HEADED RESULT FROM A MACHINE INPUT (2026-08-01).
 *
 * The rule: a routine handed a machine array answers with a machine array. The
 * elements a routine INVENTS -- a pad zero, an identity's 1, a pivot's 1, a
 * free variable's 1, a Vandermonde x^0 -- take the input's exactness, they do
 * not default to exact.
 *
 * Every expected value here was checked against Mathematica 14.0 except the
 * RowReduce block, which is a stated divergence and says so. `Head /@ Flatten`
 * rather than the printed form, because the printed form of a two-headed
 * matrix is easy to read past: {{1, 0.}, {0, 1}} looks fine.
 *
 * tools/check_array_exactness.py is the sweep this pins; these are the rows it
 * found. */
static void test_no_two_headed_results(void) {
    /* --- constructors: the elements are invented outright ---------------- */
    assert_eval_eq("Union[Head /@ Flatten[HankelMatrix[{1, 2, 3.}]]]", "{Real}", 0);
    assert_eval_eq("HankelMatrix[{1, 2, 3.}]",
                   "{{1.0, 2.0, 3.0}, {2.0, 3.0, 0.0}, {3.0, 0.0, 0.0}}", 0);
    assert_eval_eq("HankelMatrix[{1, 2, 3}]",
                   "{{1, 2, 3}, {2, 3, 0}, {3, 0, 0}}", 0);
    assert_eval_eq("Union[Head /@ Flatten[ToeplitzMatrix[{1, 2, 3.}]]]", "{Real}", 0);
    assert_eval_eq("Union[Head /@ Flatten[VandermondeMatrix[{1., 2., 3.}]]]", "{Real}", 0);
    assert_eval_eq("VandermondeMatrix[{1., 2., 3.}]",
                   "{{1.0, 1.0, 1.0}, {1.0, 2.0, 4.0}, {1.0, 3.0, 9.0}}", 0);
    assert_eval_eq("VandermondeMatrix[{1, 2, 3}]",
                   "{{1, 1, 1}, {1, 2, 4}, {1, 3, 9}}", 0);

    /* --- MatrixPower[m, 0] is the identity at m's exactness -------------- */
    assert_eval_eq("MatrixPower[{{2., 0.}, {0., 2.}}, 0]",
                   "{{1.0, 0.0}, {0.0, 1.0}}", 0);
    assert_eval_eq("MatrixPower[{{2, 0}, {0, 2}}, 0]", "{{1, 0}, {0, 1}}", 0);

    /* --- NullSpace's free-variable slot ---------------------------------- */
    assert_eval_eq("Union[Head /@ Flatten[NullSpace[{{1., 2.}, {2., 4.}}]]]", "{Real}", 0);
    assert_eval_eq("NullSpace[{{1, 2}, {2, 4}}]", "{{-2, 1}}", 0);

    /* --- RowReduce: a DELIBERATE divergence ------------------------------
     * Mathematica writes the exact 1 and 0 regardless of the input, so its
     * RowReduce[{{2., 0.}, {0., 4.}}] has heads {Integer, Real, Integer,
     * Integer} and PackedArrayQ False. Mathilda follows the project rule: a
     * machine matrix in, a machine matrix out. Values agree; only the heads of
     * the structural entries differ, and only for an inexact input. */
    assert_eval_eq("Union[Head /@ Flatten[RowReduce[{{2., 0.}, {0., 4.}}]]]", "{Real}", 0);
    assert_eval_eq("RowReduce[{{1., 2.}, {3., 4.}}]", "{{1.0, 0.0}, {0.0, 1.0}}", 0);
    assert_eval_eq("RowReduce[{{1, 2}, {3, 4}}]", "{{1, 0}, {0, 1}}", 0);

    /* --- and the mixtures that are CORRECT, so the rule is not overreached
     * into a coercion. Both match Mathematica. */
    assert_eval_eq("Chop[{1., 1.*^-12, 2.}]", "{1.0, 0, 2.0}", 0);
    assert_eval_eq("PadRight[{1., 2.}, 4]", "{1.0, 2.0, 0, 0}", 0);
    assert_eval_eq("Riffle[{1., 2.}, 0]", "{1.0, 0, 2.0}", 0);

    /* --- and MPFR does NOT trigger the contagion ------------------------- */
    assert_eval_eq("Head[HankelMatrix[{1.`30, 2}][[2, 2]]]", "Integer", 0);

    /* THE HAZARD THAT COMES WITH MAKING A HEAD PACK.
     *
     * Once RowReduce answers with a BUFFER, every internal consumer that walks
     * its result structurally breaks -- silently, and only above the
     * 250-element packing threshold. get_tensor_dims returns 0 for an NDArray,
     * so NullSpace of any machine matrix at 16x16 or larger came back
     * UNEVALUATED while 14x14 was fine. Six call sites (NullSpace x2,
     * MatrixRank x2, Inverse's mat_rref, parfrac, eigen_common) now use
     * pack_eval_plain, which is what pack.h prescribes for exactly this.
     *
     * The sizes straddle the threshold on purpose: a test at one size only
     * proves whichever side it landed on. */
    static const char* const THRESHOLD[] = {
        /* {rank-1 matrix, so the null space has dimension n-1} */
        "Length[NullSpace[Table[1. i j, {i, 14}, {j, 14}]]] == 13",
        "Length[NullSpace[Table[1. i j, {i, 16}, {j, 16}]]] == 15",
        "Length[NullSpace[Table[1. i j, {i, 20}, {j, 20}]]] == 19",
        "MatrixRank[Table[1. i j, {i, 14}, {j, 14}]] == 1",
        "MatrixRank[Table[1. i j, {i, 20}, {j, 20}]] == 1",
        "MatrixRank[Table[If[i == j, 20., 1./(1. + Abs[i - j])], {i, 20}, {j, 20}]] == 20",
        "Head[PseudoInverse[Table[1. i j, {i, 20}, {j, 20}]]] === List",
        "Dimensions[Inverse[Table[If[i == j, 20., 1./(1. + Abs[i - j])], {i, 20}, {j, 20}]]] == {20, 20}",
        "Length[Eigenvalues[Table[If[i == j, 20., 1./(1. + Abs[i - j])], {i, 20}, {j, 20}]]] == 20",
    };
    for (size_t i = 0; i < sizeof(THRESHOLD) / sizeof(THRESHOLD[0]); i++)
        assert_eval_eq(THRESHOLD[i], "True", 0);
    /* Apart goes through parfrac.c's RowReduce. */
    assert_eval_eq("Apart[1/(x^2 - 1)]", "-1/2/(1 + x) + 1/2/(-1 + x)", 0);

    /* The point of all of the above: single-headed means packable. */
    static const char* const PACKS[] = {
        "NDArrayQ[HankelMatrix[Table[N[i], {i, 30}]]]",
        "NDArrayQ[ToeplitzMatrix[Table[N[i], {i, 30}]]]",
        "NDArrayQ[VandermondeMatrix[Table[N[i], {i, 30}]]]",
        "NDArrayQ[MatrixPower[Table[N[1./(1. + Abs[i - j])], {i, 20}, {j, 20}], 0]]",
        "NDArrayQ[RowReduce[Table[If[i == j, 20., 1./(1. + Abs[i - j])], {i, 20}, {j, 20}]]]",
    };
    for (size_t i = 0; i < sizeof(PACKS) / sizeof(PACKS[0]); i++)
        assert_eval_eq(PACKS[i], "True", 0);
}

void test_no_nesting_invariant(void) {
    /* THE INVARIANT: a packed node may never come to rest inside a plain
     * EXPR_FUNCTION tree, where an unaware recursive walker would meet it. That
     * is what keeps the gate's top-level scan complete. If any of these reports
     * True the gate has become unsound and needs to be deep (or the flag-bit
     * design in the plan).
     *
     * There are two ways to honour it, and List uses BOTH. When the rows can
     * form one array they are absorbed into a single rank-2 buffer (see
     * test_list_of_packed_rows_absorbs); when they cannot -- a mixed list like
     * this one -- they are materialised as the List is built. What is never
     * allowed is the third outcome, a plain List holding buffers. */
    assert_eval_eq("NDArrayQ[{ToNDArray[{1., 2.}], 5}[[1]]]", "False", 0);
    assert_eval_eq("Head[{ToNDArray[{1., 2.}], 5}[[1]]]", "List", 0);
    assert_eval_eq("{ToNDArray[{1., 2.}], 5}", "{{1.0, 2.0}, 5}", 0);
    /* Absorbed rather than materialised, and STILL not a nested buffer: the
     * element of a packed rank-2 is a fresh packed row, not an interior node of
     * a plain List. */
    {
        Expr* r = ev("{Range[1., 300.], Range[1., 300.]}");
        ASSERT(is_packed_list(r));
        ASSERT(r->data.ndarray.rank == 2);
        expr_free(r);
    }

    /* The hole this nearly shipped with. An unevaluated application like
     * gg[xx][p] has a non-symbol head and no builtin to run, so it comes to
     * REST -- and if the gate exempts it as though it were a pure function, a
     * packed node stays nested inside a plain function node where the shallow
     * scan at the enclosing level never looks. Every one of these returned a
     * confident wrong answer while that exemption was too broad. */
    assert_eval_eq("Count[gg[xx][ToNDArray[{1., 2., 3., 4.}]], _Real, 2]", "4", 0);
    assert_eval_eq("Count[gg[xx][{1., 2., 3., 4.}], _Real, 2]", "4", 0);
    assert_eval_eq("LeafCount[gg[xx][ToNDArray[{1., 2., 3., 4.}]]]", "7", 0);
    assert_eval_eq("LeafCount[gg[xx][{1., 2., 3., 4.}]]", "7", 0);
    assert_eval_eq("Cases[gg[xx][ToNDArray[{1., 2.}]], _Real, 2]", "{1.0, 2.0}", 0);
    assert_eval_eq("gg[xx][ToNDArray[{1., 2.}]] /. 2. -> 9.", "gg[xx][{1.0, 9.0}]", 0);
    assert_eval_eq("NDArrayQ[gg[xx][ToNDArray[{1., 2.}]][[1]]]", "False", 0);

    /* A genuine pure Function IS exempt, and must stay so: substitution drops
     * the value into the body where every head is gated on the next pass, so
     * nothing comes to rest nested and the fast path is kept. */
    assert_eval_eq("(# + 1. &)[ToNDArray[{1., 2.}]]", "{2.0, 3.0}", 0);
    assert_eval_eq("NDArrayQ[(# + 1. &)[ToNDArray[{1., 2.}]]]", "True", 0);
}

/* Packed in, packed out. Every one of these took the materialise-reuse-repack
 * detour at some point in its implementation, and each would silently hand back
 * a VISIBLE NDArray[...] if the repack forgot to inherit the presentation. */
void test_aware_heads_stay_packed(void) {
    static const char* const SRCS[] = {
        "Map[#^2 &, ToNDArray[{1., 2., 3., 4.}]]",
        "Select[ToNDArray[{1., 2., 3., 4.}], # > 2. &]",
        "Rest[ToNDArray[{1., 2., 3.}]]",
        "Most[ToNDArray[{1., 2., 3.}]]",
        "Join[ToNDArray[{1., 2.}], {3.}]",
        "Partition[ToNDArray[{1., 2., 3., 4.}], 2]",
        "RotateLeft[ToNDArray[{1., 2., 3.}], 1]",
        "Differences[ToNDArray[{1., 2., 4.}]]",
        "Riffle[ToNDArray[{1., 2.}], 0.]",
        "TakeWhile[ToNDArray[{1., 2., 3.}], # < 3. &]",
        "FoldList[Plus, 0., ToNDArray[{1., 2.}]]",
        "Take[ToNDArray[{1., 2., 3.}], 2]",
        "Drop[ToNDArray[{1., 2., 3.}], 1]",
        "Accumulate[ToNDArray[{1., 2.}]]",
        "Reverse[ToNDArray[{1., 2.}]]",
        "Sort[ToNDArray[{2., 1.}]]",
        "Transpose[ToNDArray[{{1., 2.}, {3., 4.}}]]",
        "ToNDArray[{1., 2.}] + 1.",
        "ToNDArray[{1., 2.}] * 2.",
        "Sin[ToNDArray[{1., 2.}]]",
        /* Every element Commonest returns was copied out of the buffer, so the
         * answer is one dtype and stays packed for whatever reads it next. */
        "Commonest[ToNDArray[{1., 2., 2., 3.}]]",
        "Commonest[ToNDArray[{1., 2., 2., 3., 3.}], 2]",
    };
    for (size_t i = 0; i < sizeof(SRCS) / sizeof(SRCS[0]); i++) {
        Expr* r = ev(SRCS[i]);
        if (!is_packed_list(r)) {
            char* s = expr_to_string(r);
            fprintf(stderr, "FAIL: %s\n  did not stay packed; got %s\n", SRCS[i], s);
            free(s);
        }
        ASSERT(is_packed_list(r));
        expr_free(r);
    }

    /* Clip USED to be on this list and had to come off: Clip clamps to its
     * EXACT bounds, so a clipped element comes back as the Integer 1 where a
     * float64 buffer can only hold 1.. It is correct but unpacked now -- the
     * value is what matters, and the value has to be the plain list's. */
    same_as_plain("Clip[ToNDArray[{1., 5.}], {2., 3.}]", "Clip[{1., 5.}, {2., 3.}]");
    same_as_plain("Clip[ToNDArray[{1., 2., 3.}]]", "Clip[{1., 2., 3.}]");
}

/* The known, DELIBERATE cost of the no-nesting invariant.
 *
 * A head whose argument is a List *containing* arrays -- MapThread and
 * Transpose over several vectors are the two that matter -- never sees packed
 * elements at all: List is not packed-aware, so building {p, q} materialises
 * both. The answers stay right, they are just computed the ordinary way.
 *
 * This is the price of keeping the transparency gate a top-level O(argc) scan
 * instead of an O(tree) walk on every unaware head, and it is pinned here so
 * that a future change to the invariant shows up as a test to update rather
 * than as an unexplained performance cliff. Lifting it needs the transitive
 * "contains a packed node" bit described in docs/design/packed_arrays.md. */
void test_nesting_limitation_is_correct_but_unpacked(void) {
    /* Below PACK_MIN_ELEMENTS the two rows are materialised rather than
     * absorbed, so MapThread gets an ordinary List of ordinary lists and
     * answers from the generic path. Correct value, ordinary storage.
     *
     * The threshold is PINNED here rather than assumed. This arm is about what
     * happens BELOW it, so it has to name a size that is below it -- and the
     * default moved (250 -> 4 on 2026-08-02, so a 4-element result began
     * packing and this assertion began failing on a value that had not
     * changed). pack_set_min_elements exists for exactly this; 0 restores the
     * compiled-in default. */
    pack_set_min_elements(1000);
    assert_eval_eq("MapThread[Plus, {ToNDArray[{1., 2.}], ToNDArray[{3., 4.}]}]",
                   "{4.0, 6.0}", 0);
    Expr* r = ev("MapThread[Plus, {ToNDArray[{1., 2.}], ToNDArray[{3., 4.}]}]");
    ASSERT(!is_packed_list(r));     /* correct value, ordinary List storage */
    expr_free(r);
    pack_set_min_elements(0);

    /* And the SAME expression packs once the threshold is below its size --
     * which is the real content of the arm above: storage follows the
     * threshold, the answer does not. */
    {
        Expr* small = ev("MapThread[Plus, {ToNDArray[{1., 2.}], ToNDArray[{3., 4.}]}]");
        char* ss = expr_to_string(small);
        ASSERT(strcmp(ss, "{4.0, 6.0}") == 0);
        free(ss);
        expr_free(small);
    }

    /* Above it the same expression absorbs into a rank-2 buffer and MapThread's
     * column fold answers instead -- same value, packed storage. The two arms
     * together are the point: which one runs is a size decision, and it is not
     * observable in the answer. */
    {
        Expr* big = ev("MapThread[Plus, {Range[1., 300.], Range[1., 300.]}]");
        Expr* plain = ev("MapThread[Plus, {Normal[Range[1., 300.]], "
                         "Normal[Range[1., 300.]]}]");
        char* sa = expr_to_string(big);
        char* sb = expr_to_string(plain);
        ASSERT(strcmp(sa, sb) == 0);
        ASSERT(is_packed_list(big));
        free(sa); free(sb);
        expr_free(big); expr_free(plain);
    }

    /* MapIndexed is a separate, PRE-EXISTING gap, not a packing one: it never
     * repacked, on either surface. MapIndexed[f, NDArray[{1.,2.}]] has always
     * come back as a plain List while Map[f, NDArray[{1.,2.}]] comes back as an
     * NDArray -- mi_ndarray_axis returns unevaluated f[part, {i}] nodes for the
     * evaluator to reduce and never repacks the reduced result. Recorded here
     * so the inconsistency is visible; fixing it belongs with MapIndexed, not
     * with packing. */
    assert_eval_eq("MapIndexed[#1 &, ToNDArray[{1., 2.}]]", "{1.0, 2.0}", 0);
    Expr* mi = ev("MapIndexed[#1 &, ToNDArray[{1., 2.}]]");
    ASSERT(!is_packed_list(mi));
    expr_free(mi);

    /* Whereas a DIRECT packed argument does keep its packing. */
    assert_eval_eq("Join[ToNDArray[{1., 2.}], {3.}]", "{1.0, 2.0, 3.0}", 0);
    Expr* j = ev("Join[ToNDArray[{1., 2.}], {3.}]");
    ASSERT(is_packed_list(j));
    expr_free(j);
}

/* The mirror image: a VISIBLE NDArray[...] must not become a List just because
 * the same repack machinery now knows about presentations. */
void test_visible_stays_visible(void) {
    assert_eval_eq("Map[# + 1. &, NDArray[{1., 2.}]]", "NDArray[{2.0, 3.0}]", 0);
    assert_eval_eq("Select[NDArray[{1., 2., 3.}], # > 1. &]", "NDArray[{2.0, 3.0}]", 0);
    assert_eval_eq("Rest[NDArray[{1., 2., 3.}]]", "NDArray[{2.0, 3.0}]", 0);
    assert_eval_eq("Sort[NDArray[{2., 1.}]]", "NDArray[{1.0, 2.0}]", 0);
    /* A visible operand dominates a packed one in a binary op: the user asked
     * for an NDArray on one side, so they get one back. */
    assert_eval_eq("NDArray[{1., 2.}] + ToNDArray[{1., 2.}]", "NDArray[{2.0, 4.0}]", 0);
    assert_eval_eq("ToNDArray[{1., 2.}] + ToNDArray[{1., 2.}]", "{2.0, 4.0}", 0);
}


/* ---------------------------------------------------------------- Phase 4:
 * automatic packing at the producers.
 *
 * Everything above exercises packing through the explicit ToNDArray. These test
 * that the producers pack on their own, that the 250-element threshold is not
 * observable, and that the specific defects automatic packing exposed stay
 * fixed. Each of the last group was a WRONG ANSWER reachable from code that
 * never mentions an array. */

void test_producers_pack(void) {
    /* Direct construction: the buffer is written, the elements never built. */
    assert_eval_eq("NDArrayQ[Range[1., 300.]]", "True", 0);
    assert_eval_eq("NDArrayQ[Range[300]]", "True", 0);
    assert_eval_eq("DataType[Range[300]]", "\"int64\"", 0);
    assert_eval_eq("DataType[Range[1., 300.]]", "\"float64\"", 0);
    assert_eval_eq("NDArrayQ[ConstantArray[1., 300]]", "True", 0);
    assert_eval_eq("NDArrayQ[ConstantArray[7, {3, 100}]]", "True", 0);
    assert_eval_eq("NDArrayQ[RandomReal[1, 300]]", "True", 0);
    assert_eval_eq("NDArrayQ[Table[i^2, {i, 1., 300.}]]", "True", 0);
    /* Offered after building. */
    assert_eval_eq("NDArrayQ[Table[i^2, {i, 1, 300}]]", "True", 0);
    assert_eval_eq("NDArrayQ[Table[i, {i, 300}]]", "True", 0);
    assert_eval_eq("NDArrayQ[RandomInteger[10, 300]]", "True", 0);
    assert_eval_eq("NDArrayQ[Sort[RandomReal[1, 300]]]", "True", 0);
    assert_eval_eq("NDArrayQ[Select[Range[1., 300.], # > 0. &]]", "True", 0);
    assert_eval_eq("NDArrayQ[NestList[# + 1. &, 0., 300]]", "True", 0);
    assert_eval_eq("NDArrayQ[FoldList[Plus, 0., Range[1., 300.]]]", "True", 0);
    /* A body whose result is an Integer must NOT take the float64 direct path:
     * Table[1, {i, 1., 300.}] is a list of Integers, and the compiled result
     * type is what gates it. */
    assert_eval_eq("Head[Table[1, {i, 1., 300.}][[1]]]", "Integer", 0);
    assert_eval_eq("DataType[Table[1, {i, 1., 300.}]]", "\"int64\"", 0);
    /* Nested producers give a real rank-2 array, not a list of packed rows --
     * the packer absorbs rows that are already buffers. */
    assert_eval_eq("NDArrayQ[Table[i j, {i, 300}, {j, 300}]]", "True", 0);
    assert_eval_eq("Dimensions[Table[i j, {i, 300}, {j, 300}]]", "{300, 300}", 0);
    assert_eval_eq("Total[Table[i j, {i, 20}, {j, 20}], 2]", "44100", 0);
    /* Nothing symbolic or mixed packs, at any size. */
    assert_eval_eq("NDArrayQ[Table[x, {300}]]", "False", 0);
    assert_eval_eq("NDArrayQ[Table[k/2, {k, 300}]]", "False", 0);
    assert_eval_eq("NDArrayQ[Array[Sqrt[#] &, 300]]", "False", 0);
}

/* N over a packed integer array widens the whole buffer to a packed float64
 * array (numericalize_rec's EXPR_NDARRAY case + N's packed_int64_ok claim). The
 * regression: N[Range[10^6]] used to materialise to a list of boxed reals -- the
 * gate handed the int64 buffer to N as a plain list, and the element-by-element
 * rebuild dropped packing -- so every consumer downstream fell off the buffer. */
void test_n_over_integer_packs(void) {
    assert_eval_eq("NDArrayQ[N[Range[300]]]", "True", 0);
    assert_eval_eq("DataType[N[Range[300]]]", "\"float64\"", 0);
    assert_eval_eq("Head[N[Range[300]][[1]]]", "Real", 0);
    /* Values agree with the arithmetic int->real route, which always packed. */
    assert_eval_eq("N[Range[300]] === Range[1., 300.]", "True", 0);
    assert_eval_eq("Total[N[Range[300]]] == Total[Range[1., 300.]]", "True", 0);
    /* A visible integer NDArray converts to a visible float64 NDArray. */
    assert_eval_eq("DataType[N[NDArray[{1, 2, 3}, DataType -> \"int64\"]]]",
                   "\"float64\"", 0);
    assert_eval_eq("Head[N[NDArray[{1, 2, 3}, DataType -> \"int64\"]]]", "NDArray", 0);
    /* A real buffer is already machine-precision and passes through packed. */
    assert_eval_eq("NDArrayQ[N[Range[1., 300.]]]", "True", 0);
    /* With packing off the flag is respected: the source is unpacked, so is N's
     * result. (The restore runs inside the evaluated expression regardless.) */
    assert_eval_eq("Module[{r}, $AutoArrayPacking = False; "
                   "r = NDArrayQ[N[Range[300]]]; $AutoArrayPacking = True; r]",
                   "False", 0);
}

void test_threshold_is_not_observable(void) {
    /* 249 is under the threshold and 250 over it, so the two run through
     * different representations. Every observable except NDArrayQ must agree.
     *
     * The boundary is PINNED here with the test-only override rather than read
     * off PACK_MIN_ELEMENTS, because this test is about the threshold being
     * unobservable, not about where it sits: when the constant moved 250 -> 4
     * on 2026-08-02 (so a small matrix reaches LAPACK) these three assertions
     * were the only thing in the suite that noticed, and they were measuring
     * the constant instead of the invariant. */
    pack_set_min_elements(250);
    assert_eval_eq("NDArrayQ[Range[249]]", "False", 0);
    assert_eval_eq("NDArrayQ[Range[250]]", "True", 0);
    assert_eval_eq("NDArrayQ[Range[251]]", "True", 0);
    for (int n = 249; n <= 251; n++) {
        char a[192], b[192];
        static const char* const EXPRS[] = {
            "Total[Range[%d]]", "Head[Range[%d][[1]]]", "Mean[Range[%d]]",
            "Length[Range[%d]]", "Head[Range[%d]]", "Depth[Range[%d]]",
            "Total[Range[1., %d.]]", "Count[Range[%d], _Integer]",
            "MatchQ[Range[%d], {__Integer}]", "Last[Sort[Reverse[Range[%d]]]]",
            "LeafCount[Range[%d]]", "Position[Range[%d], 7]",
            "Range[%d] === Table[k, {k, %d}]",
        };
        for (size_t i = 0; i < sizeof(EXPRS) / sizeof(EXPRS[0]); i++) {
            snprintf(a, sizeof(a), EXPRS[i], n, n);
            /* The same expression with packing off is the reference. */
            pack_set_enabled(false);
            Expr* off = ev(a);
            pack_set_enabled(true);
            Expr* on = ev(a);
            char* so = expr_to_string(off);
            char* sn = expr_to_string(on);
            if (strcmp(so, sn) != 0)
                fprintf(stderr, "FAIL: threshold observable at n=%d\n  %s\n"
                                "  packed: %s\n  plain:  %s\n", n, a, sn, so);
            ASSERT(strcmp(so, sn) == 0);
            free(so); free(sn);
            expr_free(off); expr_free(on);
        }
        snprintf(b, sizeof(b), "Range[%d]", n);
        (void)b;
    }
    pack_set_min_elements(0);   /* restore PACK_MIN_ELEMENTS */
}

void test_exact_int64_arithmetic_on_a_buffer(void) {
    /* Every one of these arrives as an int64 buffer now that Range packs, and
     * each head below is on the verified int64 list -- so these run ON the
     * buffer rather than materialising. The answers must still be the exact ones.
     * Before the exact paths, Total[Range[10^6]] was 1.55x SLOWER packed than
     * plain, because the gate had to materialise first. */
    assert_eval_eq("Total[Range[300]]", "45150", 0);
    assert_eval_eq("Head[Total[Range[300]]]", "Integer", 0);
    assert_eval_eq("Mean[Range[300]]", "301/2", 0);
    assert_eval_eq("Median[Range[300]]", "301/2", 0);
    assert_eval_eq("Median[Range[299]]", "150", 0);
    assert_eval_eq("Max[Range[300]]", "300", 0);
    assert_eval_eq("Min[Range[300]]", "1", 0);
    assert_eval_eq("Last[Accumulate[Range[300]]]", "45150", 0);
    assert_eval_eq("Take[Range[300] + 1, 2]", "{2, 3}", 0);
    assert_eval_eq("Take[2 Range[300], 2]", "{2, 4}", 0);
    assert_eval_eq("Take[Range[300]^2, 2]", "{1, 4}", 0);
    assert_eval_eq("Take[Sort[-Range[300]], 2]", "{-300, -299}", 0);
    /* Overflow must promote, never wrap. */
    same_as_plain("Total[Range[1000000]^3]",
                  "Total[Table[k, {k, 1000000}]^3]");
    same_as_plain("Range[1000000] . Range[1000000]",
                  "Table[k, {k, 1000000}] . Table[k, {k, 1000000}]");
    /* A Real scalar widens; anything exact-but-not-integer goes to the List
     * path, where the answer is Rationals, radicals or exact Complex. */
    assert_eval_eq("Take[Range[300] 2.5, 2]", "{2.5, 5.0}", 0);
    assert_eval_eq("Take[Range[300] 5/2, 2]", "{5/2, 5}", 0);
    assert_eval_eq("Take[Range[300]^-1, 2]", "{1, 1/2}", 0);
    assert_eval_eq("Take[Range[300]^(1/2), 2]", "{1, Sqrt[2]}", 0);
    assert_eval_eq("Take[Range[300] + I, 2]", "{1 + I, 2 + I}", 0);
    assert_eval_eq("Take[2^Range[300], 3]", "{2, 4, 8}", 0);
    /* Sort in int64: past 2^53 two integers compare equal as doubles, so a
     * double sort would reorder them and round every element. */
    assert_eval_eq("Sort[ToNDArray[{9007199254740995, 9007199254740993, "
                   "9007199254740994}]]",
                   "{9007199254740993, 9007199254740994, 9007199254740995}", 0);
    /* Deliberately degraded: the exact answer is a Rational or a radical. */
    assert_eval_eq("Variance[Range[300]]", "7525", 0);
    assert_eval_eq("Quartiles[Range[300]]", "{151/2, 301/2, 451/2}", 0);
    assert_eval_eq("Take[Clip[Range[300], {3, 7}], 4]", "{3, 3, 3, 4}", 0);
    /* Precision is Listable in Mathilda, so a list answers per element -- the
     * point is that the packed list agrees with the plain one, not that it
     * collapses to a scalar. Its buffer fast path DID collapse it, which is why
     * Precision and Accuracy are on pack.c's NOT_AWARE list. */
    assert_eval_eq("First[Precision[Range[300]]]", "Infinity", 0);
    assert_eval_eq("Length[Precision[Range[300]]]", "300", 0);
    same_as_plain("Precision[Range[1., 300.]]", "Precision[Table[1. k, {k, 300}]]");
    same_as_plain("Accuracy[Range[1., 300.]]", "Accuracy[Table[1. k, {k, 300}]]");
    /* Same class, found by the same sweep: a real-closed kernel keeps the
     * float64 dtype, so Floor and friends wrote 1.0 where the list gives the
     * exact Integer 1, and Clip clamped to a Real where the list uses its exact
     * bound. */
    assert_eval_eq("Take[Floor[Range[1., 300.]/4], 4]", "{0, 0, 0, 1}", 0);
    assert_eval_eq("Take[Ceiling[Range[1., 300.]/4], 4]", "{1, 1, 1, 1}", 0);
    assert_eval_eq("Take[Round[Range[1., 300.]/4], 4]", "{0, 0, 1, 1}", 0);
    assert_eval_eq("Take[IntegerPart[Range[1., 300.]/4], 4]", "{0, 0, 0, 1}", 0);
    assert_eval_eq("Take[Sign[Range[1., 300.]], 2]", "{1, 1}", 0);
    assert_eval_eq("Take[Im[Range[1., 300.]], 2]", "{0, 0}", 0);
    assert_eval_eq("Take[Clip[Range[1., 300.]], 3]", "{1.0, 1, 1}", 0);
    /* An aware Listable head skips threading so its kernel can fire; when no
     * kernel matches the arity, the post-gate materialises the buffer so the
     * resting expression is the one the plain list gives. */
    same_as_plain("Mod[ToNDArray[{1., 2., 3.}]]", "Mod[{1., 2., 3.}]");
}

void test_regressions_automatic_packing_exposed(void) {
    /* (1) The evaluator's fixed-point test used to discard the gate's work.
     * expr_eq is blind to packing by design, so a step whose only effect was
     * materialising looked like no progress and the packed form was kept -- and
     * a nested Table came back as a List of 300 packed rows, so Dimensions was
     * {300} and Total[m, 2] a list instead of a number. */
    assert_eval_eq("Dimensions[Table[i j, {i, 300}, {j, 300}]]", "{300, 300}", 0);
    assert_eval_eq("Head[Table[i j, {i, 300}, {j, 300}][[1]]]", "List", 0);
    assert_eval_eq("Total[Table[i j, {i, 300}, {j, 300}], 2]", "2038522500", 0);

    /* (2) A Listable head with BOTH a plain list and a buffer threaded on the
     * plain one only and broadcast the buffer, giving an outer product. */
    assert_eval_eq("ToNDArray[{1., 2., 3.}] * {10, 20, 30}", "{10.0, 40.0, 90.0}", 0);
    assert_eval_eq("ToNDArray[{1., 2., 3.}] + {10, 20, 30}", "{11.0, 22.0, 33.0}", 0);
    assert_eval_eq("{10, 20, 30} * ToNDArray[{1., 2., 3.}]", "{10.0, 40.0, 90.0}", 0);
    /* Both packed still runs on the buffers. */
    assert_eval_eq("NDArrayQ[ToNDArray[Range[1., 300.]] * "
                   "ToNDArray[Range[1., 300.]]]", "True", 0);

    /* (3) A repack coerced elements the operation had introduced: Join and
     * Riffle repacked at the SOURCE dtype, turning an appended exact 1 into 1.. */
    same_as_plain("Take[Join[Range[1., 300.], {1}], -2]",
                  "Take[Join[Table[1. k, {k, 300}], {1}], -2]");
    assert_eval_eq("Take[Join[Range[1., 300.], {1}], -2]", "{300.0, 1}", 0);
    assert_eval_eq("Take[Riffle[Range[1., 300.], 1], 4]", "{1.0, 1, 2.0, 1}", 0);
    assert_eval_eq("Take[Map[If[# > 150., #, 1] &, Range[1., 300.]], 2]",
                   "{1, 1}", 0);

    /* (4) Range's exact branch was driven by the double shadow of its bounds,
     * and one ulp at 10^18 is 128 -- so val += 1 never advanced and this ran to
     * the 10^6-element safety cap. */
    assert_eval_eq("Range[1000000000000000000, 1000000000000000003]",
                   "{1000000000000000000, 1000000000000000001, "
                   "1000000000000000002, 1000000000000000003}", 0);
    assert_eval_eq("Length[Range[1000000000000000000, 1000000000000000003]]", "4", 0);

    /* (5) An internal evaluate() bypasses the gate: Median builds Sort[data],
     * evaluates it, then indexes the result's args -- and once Sort packed, that
     * result was a buffer, so Median[Range[300]] came back UNEVALUATED. */
    assert_eval_eq("Median[Table[k, {k, 300}]]", "301/2", 0);
    assert_eval_eq("Quartiles[Table[k, {k, 300}]]", "{151/2, 301/2, 451/2}", 0);

    /* (6) expr_to_latex had no arm for a buffer, so every packed result came
     * back with an empty latex field over the NDJSON pipe. */
    {
        Expr* p = ev("Range[1., 300.]");
        ASSERT(is_packed_list(p));
        Expr* plain = ndarray_to_nested_list(p);
        char* lp = expr_to_latex(p);
        char* ll = expr_to_latex(plain);
        ASSERT(lp != NULL && ll != NULL);
        ASSERT(lp[0] != '\0');
        if (strcmp(lp, ll) != 0)
            fprintf(stderr, "FAIL: latex mismatch\n  packed: %.80s\n  plain:  %.80s\n",
                    lp, ll);
        ASSERT(strcmp(lp, ll) == 0);
        free(lp); free(ll);
        expr_free(plain); expr_free(p);
    }
}


/* ------------------------------------------------------- $-flags and aliases
 *
 * $AutoCompilation and $AutoArrayPacking each switch off an optimisation that is
 * invisible by construction, so the assertions below are about three things:
 * that the flag actually reaches the C-side switch, that reading the symbol back
 * never lies about which path is running, and that turning one off changes speed
 * and nothing else. ToPackedArray / FromPackedArray are Mathematica's names for
 * ToNDArray / FromNDArray and must be the same function, not a near-copy. */

void test_dollar_autoarraypacking(void) {
    assert_eval_eq("$AutoArrayPacking", "True", 0);
    assert_eval_eq("NDArrayQ[Range[1., 300.]]", "True", 0);

    /* Off: nothing packs, from any producer. */
    assert_eval_eq("$AutoArrayPacking = False", "False", 0);
    assert_eval_eq("$AutoArrayPacking", "False", 0);
    assert_eval_eq("NDArrayQ[Range[1., 300.]]", "False", 0);
    assert_eval_eq("NDArrayQ[Range[300]]", "False", 0);
    assert_eval_eq("NDArrayQ[ConstantArray[1., 300]]", "False", 0);
    assert_eval_eq("NDArrayQ[RandomReal[1, 300]]", "False", 0);
    assert_eval_eq("NDArrayQ[Table[i^2, {i, 1., 300.}]]", "False", 0);
    assert_eval_eq("NDArrayQ[Sort[RandomReal[1, 300]]]", "False", 0);
    assert_eval_eq("NDArrayQ[NestList[# + 1. &, 0., 300]]", "False", 0);
    /* But the answers are unchanged, and the explicit requests still work --
     * ToNDArray is the user asking, not the system guessing. */
    assert_eval_eq("Total[Range[300]]", "45150", 0);
    assert_eval_eq("Mean[Range[300]]", "301/2", 0);
    assert_eval_eq("NDArrayQ[ToNDArray[{1., 2., 3.}]]", "True", 0);
    assert_eval_eq("NDArrayQ[ToPackedArray[{1., 2., 3.}]]", "True", 0);
    assert_eval_eq("Head[NDArray[{1., 2.}]]", "NDArray", 0);

    assert_eval_eq("$AutoArrayPacking = True", "True", 0);
    assert_eval_eq("$AutoArrayPacking", "True", 0);
    assert_eval_eq("NDArrayQ[Range[1., 300.]]", "True", 0);

    /* Anything but True/False is refused, and the symbol is rolled back to the
     * LIVE state rather than keeping a value that does not describe it. */
    assert_eval_eq("$AutoArrayPacking = 7", "7", 0);          /* Set returns its rhs */
    assert_eval_eq("$AutoArrayPacking", "True", 0);           /* not 7 */
    assert_eval_eq("NDArrayQ[Range[1., 300.]]", "True", 0);   /* still on */
    assert_eval_eq("$AutoArrayPacking = \"no\"", "\"no\"", 0);
    assert_eval_eq("$AutoArrayPacking", "True", 0);
}

void test_dollar_autocompilation(void) {
    assert_eval_eq("$AutoCompilation", "True", 0);

    /* The contract is that the compiled and interpreted paths agree, so the
     * assertion is a differential: the same expressions with the flag off must
     * give exactly what they give with it on. Covers the autocompile adapter
     * (Table over an inexact iterator) and numloop (Map, Nest, Do bodies). */
    static const char* const EXPRS[] = {
        "Take[Table[i^2, {i, 1., 20.}], 4]",
        "Take[Table[Sin[t] Exp[-t], {t, 0., 2., 0.1}], 3]",
        "Take[Table[1, {i, 1., 20.}], 3]",
        "Head[Table[1, {i, 1., 20.}][[1]]]",
        "Take[Map[#^2 &, {1., 2., 3., 4.}], 2]",
        "Nest[3.5 # (1 - #) &, 0.31, 40]",
        "Take[NestList[# + 1. &, 0., 20], 3]",
        "Fold[Plus, 0., Table[1. k, {k, 20}]]",
        "q1 = 1.; Do[q1 = 2. q1, {5}]; q1",
        "Total[Table[i^2, {i, 1., 300.}]]",
        "Take[Table[i j, {i, 1., 5.}, {j, 1., 5.}], 2]",
    };
    for (size_t i = 0; i < sizeof(EXPRS) / sizeof(EXPRS[0]); i++) {
        Expr* pr = parse_expression("$AutoCompilation = True");
        expr_free(evaluate(pr)); expr_free(pr);
        Expr* on = ev(EXPRS[i]);
        pr = parse_expression("$AutoCompilation = False");
        expr_free(evaluate(pr)); expr_free(pr);
        Expr* off = ev(EXPRS[i]);
        char* so = expr_to_string(on);
        char* sf = expr_to_string(off);
        if (strcmp(so, sf) != 0)
            fprintf(stderr, "FAIL: $AutoCompilation changed an answer\n  %s\n"
                            "  True:  %s\n  False: %s\n", EXPRS[i], so, sf);
        ASSERT(strcmp(so, sf) == 0);
        free(so); free(sf);
        expr_free(on); expr_free(off);
    }
    assert_eval_eq("$AutoCompilation = False", "False", 0);
    assert_eval_eq("$AutoCompilation", "False", 0);
    /* An explicitly built CompiledFunction is NOT affected -- the user asked for
     * that one. Only the compilation nobody asked for is switched off. */
    assert_eval_eq("cf = Compile[{{u, _Real}}, u^2 + 1.]; cf[3.]", "10.0", 0);
    /* (Head[] of a CompiledFunction object comes back unevaluated -- a
     * pre-existing gap in builtin_head for EXPR_COMPILED, unrelated to these
     * flags, so the check above is that the object still WORKS.) */
    assert_eval_eq("$AutoCompilation = True", "True", 0);
    assert_eval_eq("$AutoCompilation", "True", 0);

    assert_eval_eq("$AutoCompilation = 0", "0", 0);
    assert_eval_eq("$AutoCompilation", "True", 0);
}

void test_flags_are_independent(void) {
    /* The two switches must not be wired to each other: Table's direct packed
     * path needs BOTH (it is gated on the compiled result type), so a bug that
     * conflated them would still look right on that one producer. Range needs
     * only packing, and Nest needs only compilation. */
    assert_eval_eq("$AutoCompilation = False; $AutoArrayPacking = True; "
                   "NDArrayQ[Range[1., 300.]]", "True", 0);
    assert_eval_eq("$AutoCompilation = True; $AutoArrayPacking = False; "
                   "NDArrayQ[Range[1., 300.]]", "False", 0);
    /* Both off, then both on: the answers never move. */
    assert_eval_eq("$AutoCompilation = False; $AutoArrayPacking = False; "
                   "Total[Table[i^2, {i, 1., 300.}]]", "9.04505e+06", 0);
    assert_eval_eq("$AutoCompilation = True; $AutoArrayPacking = True; "
                   "Total[Table[i^2, {i, 1., 300.}]]", "9.04505e+06", 0);
    assert_eval_eq("{$AutoCompilation, $AutoArrayPacking}", "{True, True}", 0);
}

void test_sysflag_hook_touches_only_the_flags(void) {
    /* The assignment hook is reached by EVERY $-prefixed symbol, and it builds a
     * probe value -- evaluating the right-hand side when the assignment is
     * delayed -- so that it can validate the flag. That must happen for the two
     * flags and for nothing else: the REPL hooks ($Pre, $PreRead, $Post,
     * $PrePrint, $Epilog) live in the same namespace with deliberately HELD
     * right-hand sides, and evaluating one of those probes broke them
     * (repl_hooks_tests caught it). The name is now checked before any probe is
     * built.
     *
     * The observable: a delayed assignment to a $-symbol that is not a flag must
     * not evaluate its right-hand side even once. A counter in the RHS makes that
     * visible. */
    assert_eval_eq("hookCount = 0; $NotAFlagXyz := (hookCount = hookCount + 1; 5); "
                   "hookCount", "0", 0);
    /* Reading it once evaluates it once, not twice. */
    assert_eval_eq("$NotAFlagXyz", "5", 0);
    assert_eval_eq("hookCount", "1", 0);
    /* And the flags themselves still validate through the same hook. */
    assert_eval_eq("$AutoArrayPacking = False; $AutoArrayPacking", "False", 0);
    assert_eval_eq("$AutoArrayPacking = True; $AutoArrayPacking", "True", 0);
    /* A DELAYED assignment to a flag is validated from its evaluated value, so
     * the symbol still describes the live state rather than holding a rule. */
    assert_eval_eq("$AutoArrayPacking := False; NDArrayQ[Range[1., 300.]]", "False", 0);
    assert_eval_eq("$AutoArrayPacking = True; NDArrayQ[Range[1., 300.]]", "True", 0);
}

void test_packedarray_aliases(void) {
    /* Same function under Mathematica's name, so every property has to match --
     * including the ones an approximate re-implementation would get wrong. */
    assert_eval_eq("ToPackedArray[{1., 2., 3.}] === ToNDArray[{1., 2., 3.}]", "True", 0);
    assert_eval_eq("NDArrayQ[ToPackedArray[{1., 2., 3.}]]", "True", 0);
    assert_eval_eq("Head[ToPackedArray[{1., 2., 3.}]]", "List", 0);
    assert_eval_eq("ToPackedArray[{1., 2., 3.}]", "{1.0, 2.0, 3.0}", 0);
    assert_eval_eq("DataType[ToPackedArray[{1, 2, 3}]]", "\"int64\"", 0);
    assert_eval_eq("DataType[ToPackedArray[{1., 2.}, DataType -> \"float32\"]]",
                   "\"float32\"", 0);
    assert_eval_eq("Dimensions[ToPackedArray[{{1., 2.}, {3., 4.}}]]", "{2, 2}", 0);
    /* Declines identically, rather than throwing or half-packing. */
    assert_eval_eq("ToPackedArray[{1, 2.5}]", "{1, 2.5}", 0);
    assert_eval_eq("NDArrayQ[ToPackedArray[{1, 2.5}]]", "False", 0);
    assert_eval_eq("NDArrayQ[ToPackedArray[{1, 2, x}]]", "False", 0);
    assert_eval_eq("NDArrayQ[ToPackedArray[{}]]", "False", 0);
    /* Idempotent, and inverse of the other alias in both directions. */
    assert_eval_eq("NDArrayQ[ToPackedArray[ToPackedArray[{1., 2., 3.}]]]", "True", 0);
    assert_eval_eq("FromPackedArray[ToPackedArray[{1., 2., 3.}]]", "{1.0, 2.0, 3.0}", 0);
    assert_eval_eq("NDArrayQ[FromPackedArray[ToPackedArray[{1., 2., 3.}]]]", "False", 0);
    assert_eval_eq("FromNDArray[ToPackedArray[{1., 2.}]] === "
                   "FromPackedArray[ToNDArray[{1., 2.}]]", "True", 0);
    /* FromPackedArray also un-packs the VISIBLE head, like FromNDArray. */
    assert_eval_eq("FromPackedArray[NDArray[{1., 2.}]]", "{1.0, 2.0}", 0);
    assert_eval_eq("Head[FromPackedArray[NDArray[{1., 2.}]]]", "List", 0);
    /* Anything that is not an array at all comes back untouched. */
    assert_eval_eq("FromPackedArray[{1, 2, x}]", "{1, 2, x}", 0);
    assert_eval_eq("FromPackedArray[7]", "7", 0);
    assert_eval_eq("ToPackedArray[7]", "7", 0);
    /* Attributes and docstrings, which the project requires of every builtin. */
    assert_eval_eq("Attributes[ToPackedArray]", "{Protected}", 0);
    assert_eval_eq("Attributes[FromPackedArray]", "{Protected}", 0);
    ASSERT(symtab_get_docstring("ToPackedArray") != NULL);
    ASSERT(symtab_get_docstring("FromPackedArray") != NULL);
    ASSERT(symtab_get_docstring("$AutoCompilation") != NULL);
    ASSERT(symtab_get_docstring("$AutoArrayPacking") != NULL);
    /* The whole packed surface works through the alias, not just its identity. */
    assert_eval_eq("Total[ToPackedArray[{1, 2, 3}]]", "6", 0);
    assert_eval_eq("Head[Total[ToPackedArray[{1, 2, 3}]]]", "Integer", 0);
    assert_eval_eq("Take[ToPackedArray[Table[1. k, {k, 300}]], 2]", "{1.0, 2.0}", 0);
}


/* --------------------------------------------------------- the Compile[] boundary
 *
 * A CompiledFunction must answer with the head the INTERPRETER would give for the
 * same input, and automatic packing added a third possible input where there were
 * two. Before this, a packed argument had its buffer materialised by the gate
 * before compiled_function_apply ever saw it, so f[Range[1., 200000.]] measured
 * ~75x slower than f[NDArray[Range[1., 200000.]]] -- values differing only in
 * `present_as`. The assertions below are about the head rules and the exactness
 * rules; the speed is in the changelog. */

void test_compile_boundary_presentation(void) {
    assert_eval_eq("fr = Compile[{{u, _Real, 1}}, u^2 + 1.]; 0", "0", 0);
    /* Packed in -> packed out, at any size: a DERIVED array inherits its source's
     * presentation, exactly as Sin[packed] does, with no threshold re-applied. */
    assert_eval_eq("NDArrayQ[fr[ToNDArray[{1., 2., 3.}]]]", "True", 0);
    assert_eval_eq("Head[fr[ToNDArray[{1., 2., 3.}]]]", "List", 0);
    assert_eval_eq("fr[ToNDArray[{1., 2., 3.}]]", "{2.0, 5.0, 10.0}", 0);
    assert_eval_eq("NDArrayQ[fr[Range[1., 300.]]]", "True", 0);
    /* Plain List in -> plain List out, because that is what the interpreter's own
     * threading gives. */
    assert_eval_eq("NDArrayQ[fr[{1., 2., 3.}]]", "False", 0);
    assert_eval_eq("fr[{1., 2., 3.}]", "{2.0, 5.0, 10.0}", 0);
    /* The visible head is untouched, and dominates in a mixed call. */
    assert_eval_eq("fr[NDArray[{1., 2., 3.}]]", "NDArray[{2.0, 5.0, 10.0}]", 0);
    assert_eval_eq("Head[fr[NDArray[{1., 2., 3.}]]]", "NDArray", 0);

    /* A BUILT result has no presentation to inherit, so it follows the PRODUCER
     * rule -- the threshold and $AutoArrayPacking, exactly like ConstantArray
     * itself. Below the threshold that still yields a plain List, which is why
     * the rule is pack_offer and not a flag flip. */
    assert_eval_eq("fb = Compile[{{n, _Integer}}, ConstantArray[1., n]]; "
                   "{NDArrayQ[fb[300]], Length[fb[300]]}", "{True, 300}", 0);
    /* Length 3, not 5: PACK_MIN_ELEMENTS moved 250 -> 4 on 2026-08-02 and this
     * assertion, which is about the producer rule and not about where the
     * threshold sits, went on naming a length that is now ABOVE it. */
    assert_eval_eq("{NDArrayQ[fb[3]], fb[3]}", "{False, {1.0, 1.0, 1.0}}", 0);

    /* A COMPLEX result never wears the packed presentation. A zero imaginary part
     * materialises as Complex[re, 0.], which the evaluator never produces (it
     * folds to a Real), so a packed complex list would not round-trip: the same
     * call printed {1. + 0.*I, 4. + 0.*I} against the plain List's {1., 4.}. */
    assert_eval_eq("fc = Compile[{{u, _Complex, 1}}, u^2]; fc[ToNDArray[{1., 2.}]]",
                   "{1.0, 4.0}", 0);
    assert_eval_eq("NDArrayQ[fc[ToNDArray[{1., 2.}]]]", "False", 0);
    same_as_plain("fc[ToNDArray[{1., 2., 3.}]]", "fc[{1., 2., 3.}]");
}

void test_compile_boundary_exactness(void) {
    /* An int64 buffer reaching an _Integer parameter stays exact, elements and
     * heads. This combination could not arise while only the compiler made int64
     * buffers; Range[n] makes them now. */
    assert_eval_eq("fi = Compile[{{u, _Integer, 1}}, u * 2 + 1]; "
                   "Take[fi[Range[300]], 3]", "{3, 5, 7}", 0);
    assert_eval_eq("Head[fi[Range[300]][[1]]]", "Integer", 0);
    same_as_plain("Take[fi[Range[300]], 4]", "Take[fi[Table[k, {k, 300}]], 4]");
    /* Exact past 2^53, where a double round trip would lose the low bits. */
    assert_eval_eq("fi[ToNDArray[{9007199254740993, 9007199254740995}]]",
                   "{18014398509481987, 18014398509481991}", 0);
    /* Overflow inside a compiled body abandons and the interpreter promotes,
     * the same contract the scalar path has. */
    same_as_plain("Compile[{{u, _Integer, 1}}, u^3][Range[3000000, 3000002]]",
                  "Range[3000000, 3000002]^3");

    /* An int64 buffer at a _Real parameter is CAST in one O(n) pass rather than
     * declining the whole call to the interpreter (which is what the VM's own
     * dtype guard did, at ~31x the cost). */
    assert_eval_eq("Take[fr[Range[300]], 3]", "{2.0, 5.0, 10.0}", 0);
    same_as_plain("Take[fr[Range[300]], 3]", "Take[fr[Table[1. k, {k, 300}]], 3]");
    /* But the cast DECLINES rather than rounds when there is no exact double, so
     * the interpreter gives the answer instead. */
    same_as_plain("fr[ToNDArray[{9007199254740993, 9007199254740995}]]",
                  "fr[{9007199254740993, 9007199254740995}]");
    /* Narrowing a float buffer into an _Integer slot is a value change, not a
     * representation change: declined, and the interpreter answers. */
    same_as_plain("fi[ToNDArray[{1., 2., 3.}]]", "fi[{1., 2., 3.}]");
}

void test_compile_boundary_listable(void) {
    /* RuntimeAttributes -> Listable threads by RANK and repacks by re-sniffing the
     * dtype. Forcing float64 (which the repack used to do) turned an
     * integer-valued body over a packed integer list into a list of Reals. */
    assert_eval_eq("fL = Compile[{{x, _Real}}, x^2 + 1., "
                   "RuntimeAttributes -> {Listable}]; "
                   "NDArrayQ[fL[Range[1., 300.]]]", "True", 0);
    assert_eval_eq("Take[fL[Range[1., 300.]], 3]", "{2.0, 5.0, 10.0}", 0);
    assert_eval_eq("fLi = Compile[{{x, _Integer}}, x + 1, "
                   "RuntimeAttributes -> {Listable}]; "
                   "Take[fLi[Range[300]], 3]", "{2, 3, 4}", 0);
    assert_eval_eq("Head[fLi[Range[300]][[1]]]", "Integer", 0);
    /* Plain List and visible NDArray keep their own answers. */
    assert_eval_eq("fL[{1., 2.}]", "{2.0, 5.0}", 0);
    assert_eval_eq("fL[NDArray[{1., 2.}]]", "NDArray[{2.0, 5.0}]", 0);
    /* Rank-2 threading over a packed matrix. */
    assert_eval_eq("Dimensions[fL[Table[1. i j, {i, 20}, {j, 20}]]]", "{20, 20}", 0);
    same_as_plain("fL[ToNDArray[{{1., 2.}, {3., 4.}}]]", "fL[{{1., 2.}, {3., 4.}}]");
}

void test_compile_integer_array_scalar_exactness(void) {
    /* Found while wiring the boundary, and older than packing: a scalar operand of
     * an array op could only be Real or Complex, so the literal in
     * Compile[{{u, _Integer, 1}}, u * 2] was loaded into a REAL register and boxed
     * as a Real -- and the interpreter's elementwise path then correctly widened
     * the int64 buffer, answering {2., 4., 6.} where Range[3] * 2 is {2, 4, 6}.
     * A wrong COMPILED answer, which the engine may not produce. Fixed with an
     * AK_INT operand kind. */
    assert_eval_eq("Compile[{{u, _Integer, 1}}, u * 2][{1, 2, 3}]", "{2, 4, 6}", 0);
    assert_eval_eq("Compile[{{u, _Integer, 1}}, u + 1][{1, 2, 3}]", "{2, 3, 4}", 0);
    assert_eval_eq("Compile[{{u, _Integer, 1}}, u * 2 + 3][{1, 2, 3}]", "{5, 7, 9}", 0);
    same_as_plain("Compile[{{u, _Integer, 1}}, u * 2][{1, 2, 3}]", "{1, 2, 3} * 2");
    /* A Part assignment writes an exact right-hand side into an int64 buffer for
     * the same reason. */
    assert_eval_eq("Compile[{{u, _Integer, 1}}, Module[{v = u}, v[[1]] = 5; v]]"
                   "[{1, 2, 3}]", "{5, 2, 3}", 0);
    /* A Real element type still makes the whole operation inexact, which is the
     * interpreter's answer too. */
    assert_eval_eq("Compile[{{u, _Real, 1}}, u * 2][{1., 2., 3.}]", "{2.0, 4.0, 6.0}", 0);
    assert_eval_eq("Compile[{{u, _Integer, 1}}, u * 2.][{1, 2, 3}]", "{2.0, 4.0, 6.0}", 0);
}

void test_map_over_a_packed_list_stays_compiled(void) {
    /* Map over a packed list must reach numloop's COMPILED loop, not the ndarray
     * leading-axis walk that applies f through the interpreter per element.
     * Measured at 10^6: Map[#^2 &, x] was 424 ms packed against 222 ms plain and
     * Map[Sin[#] Exp[-#] &, x] 1120 ms against 180 ms -- automatic packing made
     * the most-used functional head up to 6x SLOWER until numloop_map learned to
     * read a buffer. */
    assert_eval_eq("NDArrayQ[Map[#^2 &, ToNDArray[{1., 2., 3.}]]]", "True", 0);
    assert_eval_eq("Map[#^2 &, ToNDArray[{1., 2., 3.}]]", "{1.0, 4.0, 9.0}", 0);
    /* Derived, so packed at ANY size -- four elements are under the producer
     * threshold and must not be materialised by it. */
    assert_eval_eq("NDArrayQ[Map[#^2 &, ToNDArray[{1., 2., 3., 4.}]]]", "True", 0);
    same_as_plain("Map[#^2 &, ToNDArray[{1., 2., 3.}]]", "Map[#^2 &, {1., 2., 3.}]");
    same_as_plain("Take[Map[Sin[#] Exp[-#] &, Range[1., 300.]], 3]",
                  "Take[Map[Sin[#] Exp[-#] &, Table[1. k, {k, 300}]], 3]");
    /* An exact int64 source: the dtype answers "every element is exact" in O(1),
     * and the all-inexact rule then needs the body to carry a Real -- so
     * Map[# &, packedInts] passes them through exactly, as the plain List does. */
    same_as_plain("Take[Map[# &, Range[300]], 3]", "Take[Map[# &, Table[k, {k, 300}]], 3]");
    same_as_plain("Take[Map[#^2. &, Range[300]], 3]",
                  "Take[Map[#^2. &, Table[k, {k, 300}]], 3]");
    /* A non-numeric body still falls through to the general mapper. */
    same_as_plain("Take[Map[f, Range[1., 300.]], 2]",
                  "Take[Map[f, Table[1. k, {k, 300}]], 2]");
    /* Rank 2 maps over ROWS, which numloop's scalar body cannot do -- it must
     * decline and leave that to the ndarray path. */
    same_as_plain("Map[Total, ToNDArray[{{1., 2.}, {3., 4.}}]]",
                  "Map[Total, {{1., 2.}, {3., 4.}}]");
    /* Scan over a packed list answers Null and still stops where the interpreter
     * stops. */
    assert_eval_eq("Scan[#^2 &, ToNDArray[{1., 2., 3.}]]", "Null", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_invisible_across_ranks_and_dtypes);
    TEST(test_presents_as_list);
    TEST(test_visible_ndarray_unchanged);
    TEST(test_exactness_preserved);
    TEST(test_int64_matches_plain_integer_lists);
    TEST(test_bool_packing_matches_plain);
    TEST(test_tally_matches_plain_lists);
    TEST(test_setops_match_plain_lists);
    TEST(test_arithmetic_rewrites_match_plain_lists);
    TEST(test_commonest_matches_plain_lists);
    TEST(test_declines_what_it_cannot_represent);
    TEST(test_ordering_matches_plain_lists);
    TEST(test_ordering_builtin_matches_plain);
    TEST(test_hash_consumers);
    TEST(test_value_semantics_on_assignment);
    TEST(test_part_assignment_preserves_heads);
    TEST(test_part_inherits_presentation);
    TEST(test_roundtrip_builtins);
    TEST(test_unaware_heads_are_correct);
    TEST(test_pattern_matching);
    TEST(test_third_sweep_fast_paths);
    TEST(test_fourth_sweep_fast_paths);
    TEST(test_fifth_sweep_fast_paths);
    TEST(test_seventh_round_fast_paths);
    TEST(test_no_two_headed_results);
    TEST(test_no_nesting_invariant);
    TEST(test_aware_heads_stay_packed);
    TEST(test_nesting_limitation_is_correct_but_unpacked);
    TEST(test_visible_stays_visible);
    TEST(test_kill_switch);
    TEST(test_producers_pack);
    TEST(test_n_over_integer_packs);
    TEST(test_threshold_is_not_observable);
    TEST(test_exact_int64_arithmetic_on_a_buffer);
    TEST(test_regressions_automatic_packing_exposed);
    TEST(test_dollar_autoarraypacking);
    TEST(test_dollar_autocompilation);
    TEST(test_flags_are_independent);
    TEST(test_sysflag_hook_touches_only_the_flags);
    TEST(test_packedarray_aliases);
    TEST(test_compile_boundary_presentation);
    TEST(test_compile_boundary_exactness);
    TEST(test_compile_boundary_listable);
    TEST(test_compile_integer_array_scalar_exactness);
    TEST(test_map_over_a_packed_list_stays_compiled);

    printf("All packed-list tests passed.\n");
    return 0;
}
