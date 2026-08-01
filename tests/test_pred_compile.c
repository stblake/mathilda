/*
 * test_pred_compile.c -- the compiled boolean predicate path (src/numloop.c,
 * src/funcprog.c).
 *
 * Seven heads apply a test to every element -- Select, AllTrue, AnyTrue,
 * NoneTrue, TakeWhile, LengthWhile, SelectFirst. Each now compiles an ordinary
 * comparison predicate to bytecode and runs it over the float64 buffer instead
 * of materialising the buffer and calling the test through the interpreter.
 *
 * WHY EVERY ASSERTION HERE IS DIFFERENTIAL. The fast path and the interpreter
 * must answer identically, and the only trustworthy statement of that is to
 * evaluate THE SAME SOURCE TWICE -- once with numloop_set_enabled(true), once
 * with it false -- and compare. A written-down expectation can be wrong in the
 * same direction as the code, which docs/design/performance.md §13 records
 * happening three separate times in this tree.
 *
 * THE CASE THAT MATTERS MOST is tolerance_is_reproduced() below. Mathilda does
 * not compare machine reals with C's `<`: compare_numeric in src/comparisons.c
 * treats two inexact operands as EQUAL when they agree to a relative 2^-46, so
 * Less[1., 1. + 2.^-47] is False where a naive compiled comparison says True.
 * That difference is invisible on ordinary data -- every other case in this file
 * passes with or without the tolerance -- and would have been a silent wrong
 * answer on exactly the near-tie data a numerical program produces.
 *
 * The DECLINE cases are equally load-bearing: an int64 buffer must NOT reach
 * this path, because two exact Integers compare through GMP with no tolerance
 * at all. A head that quietly accepted one would answer a different question.
 *
 * Every source below is above PACK_MIN_ELEMENTS (250). A shorter one is not
 * packed, the buffer path never runs, and the test would pass while testing
 * nothing -- the third appearance of that lesson in this tree.
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "numloop.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Evaluate `src` and return its printed form (caller frees). */
static char* eval_str(const char* src) {
    Expr* p = parse_expression(src);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(r);
    expr_free(r);
    return s;
}

/* THE differential assertion: the same source, compiled and interpreted, must
 * print identically. */
static void same_both_ways(const char* src) {
    numloop_set_enabled(true);
    char* fast = eval_str(src);
    numloop_set_enabled(false);
    char* slow = eval_str(src);
    numloop_set_enabled(true);

    if (strcmp(fast, slow) != 0)
        fprintf(stderr, "FAIL: %s\n  compiled:    %s\n  interpreted: %s\n",
                src, fast, slow);
    ASSERT(strcmp(fast, slow) == 0);
    free(fast);
    free(slow);
}

/* As same_both_ways, and additionally pin the answer, for the cases where
 * "both paths agree" would still be satisfied by both being wrong. */
static void same_and_equals(const char* src, const char* expect) {
    same_both_ways(src);
    char* got = eval_str(src);
    if (strcmp(got, expect) != 0)
        fprintf(stderr, "FAIL: %s\n  got:      %s\n  expected: %s\n",
                src, got, expect);
    ASSERT(strcmp(got, expect) == 0);
    free(got);
}

/* ---------------------------------------------------------------- */

static void quantifiers(void) {
    same_and_equals("AllTrue[Range[1., 300.], # > 0. &]", "True");
    same_and_equals("AllTrue[Range[1., 300.], # > 100. &]", "False");
    same_and_equals("AnyTrue[Range[1., 300.], # > 299. &]", "True");
    same_and_equals("AnyTrue[Range[1., 300.], # > 1000. &]", "False");
    same_and_equals("NoneTrue[Range[1., 300.], # > 1000. &]", "True");
    same_and_equals("NoneTrue[Range[1., 300.], # > 299. &]", "False");

    /* The short-circuit must not change the answer when the deciding element is
     * first, last, or absent. */
    same_and_equals("AllTrue[Range[1., 300.], # > 1. &]", "False");   /* first */
    same_and_equals("AllTrue[Range[1., 300.], # < 300. &]", "False"); /* last  */
}

static void selection(void) {
    same_both_ways("Select[Range[1., 300.], # > 150. &]");
    same_both_ways("Select[Range[1., 300.], # > 1000. &]");   /* empty result */
    same_both_ways("Select[Range[1., 300.], # > 0. &]");      /* everything   */
    same_both_ways("TakeWhile[Range[1., 300.], # < 50. &]");
    same_and_equals("LengthWhile[Range[1., 300.], # < 50. &]", "49");
    same_and_equals("LengthWhile[Range[1., 300.], # > 1000. &]", "0");
    same_and_equals("LengthWhile[Range[1., 300.], # > 0. &]", "300");
    same_and_equals("SelectFirst[Range[1., 300.], # > 250. &]", "251.0");
    same_both_ways("SelectFirst[Range[1., 300.], # > 1000. &]");
    same_and_equals("SelectFirst[Range[1., 300.], # > 1000. &, -1.]", "-1.0");
}

/* Chained comparisons. `a < # < b` is the ordinary way to write a band, and it
 * does NOT parse as a three-argument Less -- it is
 * Inequality[a, Less, Slot[1], Less, b], operands and operator symbols
 * alternating. The first version of the compiler generalised Less to argc >= 2
 * and changed nothing at all, because this head never reached that branch. */
static void chained_comparisons(void) {
    same_both_ways("Select[Range[1., 300.], 100. < # < 200. &]");
    same_both_ways("Select[Range[1., 300.], 100. <= # <= 200. &]");
    same_both_ways("Select[Range[1., 300.], 200. > # > 100. &]");
    same_both_ways("Select[Range[1., 300.], 100. < # < 150. < 200. &]");
    same_and_equals("LengthWhile[Range[1., 300.], 0. < # < 50. &]", "49");
    same_both_ways("AllTrue[Range[1., 300.], 0. < # < 400. &]");
    same_both_ways("AnyTrue[Range[1., 300.], 100. < # < 200. &]");
    /* A MIXED chain, where the two operators differ -- each consecutive pair
     * must take its own operator, not the first one. */
    same_both_ways("Select[Range[1., 300.], 100. < # <= 200. &]");
    same_both_ways("Select[Range[1., 300.], 100. <= # < 200. &]");
    /* An operator outside the order set anywhere in the chain must decline the
     * whole predicate rather than the pair. */
    same_both_ways("Select[Range[1., 300.], 100. < # != 200. &]");
}

/* And / Or / Not, and comparison operands that are themselves expressions. */
static void predicate_shapes(void) {
    same_both_ways("Select[Range[1., 300.], 100. < # < 200. &]");
    same_both_ways("Select[Range[1., 300.], # < 10. || # > 290. &]");
    same_both_ways("Select[Range[1., 300.], Not[# > 150.] &]");
    same_both_ways("Select[Range[1., 300.], Abs[# - 150.] < 10. &]");
    same_both_ways("Select[Range[1., 300.], #^2. > 10000. &]");
    same_both_ways("Select[Range[1., 300.], Sqrt[#] > 10. &]");
    same_both_ways("Select[Range[1., 300.], 2. # + 1. > 100. &]");
    same_both_ways("Select[Range[1., 300.], 150. > # &]");        /* var on right */
    same_both_ways("Select[Range[1., 300.], # >= 150. &]");
    same_both_ways("Select[Range[1., 300.], # <= 150. &]");
    /* Named-parameter form, not just Slot. */
    same_both_ways("Select[Range[1., 300.], Function[x, x > 150.]]");
    same_both_ways("AllTrue[Range[1., 300.], Function[x, x > 0.]]");
}

/* THE tolerance case. 2^-47 is about 7.1e-15 relative at 1.0, inside
 * compare_numeric's 2^-46 (1.42e-14) band, so the interpreter says
 * 1. + 2.^-47 is NOT greater than 1. -- the two compare EQUAL. A compiled `>`
 * without the tolerance would answer the opposite, and only here.
 *
 * The list is built so the near-tie sits among ordinary values and the result
 * is above the packing threshold either way. */
static void tolerance_is_reproduced(void) {
    /* Sanity: the scalar comparison really does answer False, so the case below
     * is discriminating rather than vacuous. */
    char* s = eval_str("1. + 2.^-47 > 1.");
    if (strcmp(s, "False") != 0)
        fprintf(stderr, "FAIL: tolerance premise changed -- `1. + 2.^-47 > 1.` "
                        "is %s, so the differential case below tests nothing\n", s);
    ASSERT(strcmp(s, "False") == 0);
    free(s);

    same_both_ways("Select[Join[Range[2., 300.], {1., 1. + 2.^-47}], # > 1. &]");
    same_both_ways("Count[Select[Join[Range[2., 300.], {1., 1. + 2.^-47}],"
                   " # > 1. &], _]");
    same_and_equals("AllTrue[Join[Range[2., 300.], {1. + 2.^-47}], # > 1. &]",
                    "False");
    same_both_ways("AnyTrue[{1. + 2.^-47}, # > 1. &]");

    /* A difference just OUTSIDE the band must still compare greater, so the
     * tolerance is not swallowing everything. */
    same_and_equals("1. + 2.^-40 > 1.", "True");
}

/* Shapes the fast path must decline. Each must still give the interpreter's
 * answer -- the point is that declining is silent and correct, not that it
 * errors. */
static void declines(void) {
    /* An int64 buffer: exact Integers compare with no tolerance, a different
     * rule, so pred_open rejects the dtype outright. */
    same_both_ways("Select[Range[300], # > 150 &]");
    same_both_ways("AllTrue[Range[300], # > 0 &]");
    same_both_ways("Select[Range[300], EvenQ]");

    /* Bodies outside the compilable subset. */
    same_both_ways("Select[Range[1., 300.], PrimeQ[Round[#]] &]");
    same_both_ways("Select[Range[1., 300.], # == 150. &]");
    same_both_ways("AllTrue[Range[1., 300.], NumberQ[#] &]");
    same_both_ways("Select[Range[1., 300.], StringQ[#] &]");

    /* A predicate that is not a pure Function at all. */
    same_both_ways("Select[Range[1., 300.], Positive]");

    /* Non-finite intermediates: Log of a non-positive goes complex in the
     * interpreter, so the compiled scan must abandon and let it. */
    same_both_ways("Select[Range[-150., 149.], Log[#] > 0. &]");
    same_both_ways("Select[Range[1., 300.], 1./(# - 150.) > 0. &]");

    /* Not a rank-1 numeric vector. */
    same_both_ways("Select[Partition[Range[1., 300.], 3], Length[#] > 2 &]");
    same_both_ways("Select[Range[1., 300.] + I, Re[#] > 150. &]");
}

/* Select's counted form takes the old path (the fast one has no n limit), and
 * must be unaffected. */
static void counted_select(void) {
    same_both_ways("Select[Range[1., 300.], # > 150. &, 5]");
    same_both_ways("Select[Range[1., 300.], # > 1000. &, 5]");
}

int main(void) {
    symtab_init();
    core_init();

    quantifiers();
    selection();
    chained_comparisons();
    predicate_shapes();
    tolerance_is_reproduced();
    declines();
    counted_select();

    symtab_clear();
    printf("test_pred_compile: all checks passed\n");
    return 0;
}
