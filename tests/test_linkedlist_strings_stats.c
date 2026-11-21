/* LinkedList everywhere a List works: strings, statistics, distances, dates.
 *
 * Every check runs a head on a LinkedList and on Normal of the same list (a
 * typed list's Normal is PACKED, so a machine statistic must match the packed
 * kernel to the bit) and requires the same answer with every LinkedList in the
 * result replaced by its elements, no materialisation of the argument
 * (Internal`LinkedListConversionCount[] does not move) and the argument
 * unchanged. Fixture symbols are prefixed (s9...) so the other suites' bindings
 * of a, b, x, ... cannot leak in.
 */
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "parse.h"
#include "eval.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* llchk[h, l]: h on LinkedList[l] against h on its Normal, held -- an
 * unevaluated h[ll] is never rescued by evaluating h[{...}] afterwards. */
static const char* LLCHK =
    "llchk[h_, l0_] := Module[{ll = LinkedList[l0], l, c0, r, c1}, l = Normal[ll]; "
    "c0 = Internal`LinkedListConversionCount[]; r = h[ll]; "
    "c1 = Internal`LinkedListConversionCount[]; "
    "{c1 == c0, (Hold[Evaluate[r]] //. x_LinkedList :> With[{v = Normal[x]}, v /; True]) === "
    "Hold[Evaluate[h[l]]], Normal[ll] === l}]; ";

/* llmchk[h, m]: the matrix shapes -- a List of LinkedList rows, and a
 * LinkedList of List rows -- against h on the matrix of their Normals. */
static const char* LLMCHK =
    "llmchk[h_, m0_] := Module[{rows = LinkedList /@ m0, outer = LinkedList[m0], ref, c0, r1, r2, c1, "
    "nl = (Hold[#] //. x_LinkedList :> With[{v = Normal[x]}, v /; True]) &}, "
    "ref = Normal /@ rows; c0 = Internal`LinkedListConversionCount[]; "
    "r1 = h[rows]; r2 = h[outer]; c1 = Internal`LinkedListConversionCount[]; "
    "{c1 == c0, nl[r1] === Hold[Evaluate[h[ref]]], nl[r2] === Hold[Evaluate[h[Normal[outer]]]], "
    "(Normal /@ rows) === ref}]; ";

static void expect_llchk(const char* heads, const char* fixtures) {
    size_t len = strlen(LLCHK) + strlen(heads) + strlen(fixtures) + 128;
    char* buf = malloc(len);
    snprintf(buf, len, "%sUnion[Flatten[Quiet[Table[llchk[h, l], {h, %s}, {l, %s}]]]]",
             LLCHK, heads, fixtures);
    assert_eval_eq(buf, "{True}", 0);
    free(buf);
}

static void expect_llmchk(const char* heads, const char* fixtures) {
    size_t len = strlen(LLMCHK) + strlen(heads) + strlen(fixtures) + 128;
    char* buf = malloc(len);
    snprintf(buf, len, "%sUnion[Flatten[Quiet[Table[llmchk[h, m], {h, %s}, {m, %s}]]]]",
             LLMCHK, heads, fixtures);
    assert_eval_eq(buf, "{True}", 0);
    free(buf);
}

/* ------------------------------------------------------------ statistics */

/* Machine-real vectors longer than the pairwise block (128) and the thread
 * split, so the summation tree is exercised; ints and rationals (exact),
 * symbols, complexes, and the empty list. */
static const char* STAT_FIXTURES =
    "{N[Sin[Range[1000]]], N[Range[-40, 300]]/7, Range[-5, 14], Range[7]/3, "
    "{s9a, s9b, s9a, 2}, N[Range[300]] (1. + 0.5 I), {1.5, -2., 3., 0.5}, {}}";

void test_stats_vectors() {
    expect_llchk("{Mean, Variance, StandardDeviation, RootMeanSquare, Skewness, Kurtosis, "
                 "Median, Quartiles, InterquartileRange, MeanDeviation, MedianDeviation, "
                 "Moment[#, 3] &, Moment[#, 2.5] &, CentralMoment[#, 2] &, CentralMoment[#, 3] &}",
                 STAT_FIXTURES);
}

void test_stats_moving_and_scaling() {
    expect_llchk("{MovingAverage[#, 3] &, MovingAverage[#, {1, 2, 3}] &, MovingMedian[#, 3] &, "
                 "ExponentialMovingAverage[#, 0.3] &, ExponentialMovingAverage[#, 1/3] &, "
                 "Standardize, Rescale, Rescale[#, {-1, 1}] &, Rescale[#, {0, 2}, {5, 6}] &}",
                 STAT_FIXTURES);
    /* a LinkedList of weights where the List form takes a List */
    assert_eval_eq("s9w = N[Range[20]]^2; {MovingAverage[s9w, LinkedList[{1, 2, 1}]] === "
                   "MovingAverage[s9w, {1, 2, 1}], Head[MovingAverage[LinkedList[s9w], 4]]}",
                   "{True, LinkedList}", 0);
}

void test_stats_two_vectors() {
    expect_llchk("{Covariance[#, #] &, Correlation[#, #] &, Covariance[#, # + 1] &, "
                 "EuclideanDistance[#, #] &, SquaredEuclideanDistance[#, # + 2] &, "
                 "ManhattanDistance[#, 0 #] &, CosineDistance[#, # + 1] &}", STAT_FIXTURES);
    /* a LinkedList against a List */
    assert_eval_eq("s9v = N[Range[200]]/9; {Covariance[LinkedList[s9v], Reverse[s9v]] === "
                   "Covariance[s9v, Reverse[s9v]], EuclideanDistance[LinkedList[{1, 2}], {4, 6}]}",
                   "{True, 5}", 0);
}

/* Columnwise: a List of LinkedList rows and a LinkedList of List rows. */
void test_stats_matrices() {
    expect_llmchk("{Mean, Variance, StandardDeviation, RootMeanSquare, Skewness, Kurtosis, Median, "
                  "Quartiles, InterquartileRange, MeanDeviation, MedianDeviation, Moment[#, 2] &, "
                  "CentralMoment[#, 3] &, Covariance, Correlation, Standardize, Rescale}",
                  "{Partition[N[Sin[Range[60]]], 3], Partition[Range[12], 3], "
                  "{{s9a, 1}, {s9b, 2}}, Partition[N[Cos[Range[400]]], 2]}");
    /* a matrix answer comes back in the image of its rows */
    assert_eval_eq("{Head[Mean[{LinkedList[{1, 2}], LinkedList[{3, 5}]}]], "
                   "Head[Mean[LinkedList[{{1, 2}, {3, 5}}]]], Normal[Mean[{LinkedList[{1, 2}], LinkedList[{3, 5}]}]]}",
                   "{LinkedList, List, {2, 7/2}}", 0);
}

/* The machine reductions read the chunks of a list grown at both ends (the
 * front run is stored downward) and agree with the packed kernel to the bit. */
void test_stats_both_runs_bit_parity() {
    assert_eval_eq("s9q = LinkedList[N[Sin[Range[1500]]]]; Do[PrependTo[s9q, N[Cos[s9i]]], {s9i, 1500}]; "
                   "s9p = Normal[s9q]; "
                   "Union[Table[(h[s9q] /. s9l_LinkedList :> Normal[s9l]) === h[s9p], "
                   "{h, {Mean, Variance, StandardDeviation, RootMeanSquare, Skewness, Kurtosis, Median, "
                   "Quartiles, CentralMoment[#, 5] &, Moment[#, 4] &, MovingAverage[#, 7] &, "
                   "MovingMedian[#, 4] &, ExponentialMovingAverage[#, 0.2] &}}]]",
                   "{True}", 0);
}

/* --------------------------------------------------------------- strings */

static const char* STR_FIXTURES =
    "{{\"b\", \"a\", \"c\", \"a\"}, {\"hello world\", \"  foo bar \", \"x,y,z\"}, {}}";

void test_string_subjects() {
    expect_llchk("{StringJoin, StringRiffle, StringRiffle[#, \",\"] &, StringRiffle[#, {\"(\", \"|\", \")\"}] &, "
                 "StringSplit, StringSplit[#, \",\"] &, StringCases[#, \"o\" | \"a\"] &, "
                 "StringContainsQ[#, \"a\"] &, StringFreeQ[#, \"a\"] &, StringStartsQ[#, \"h\"] &, "
                 "StringEndsQ[#, \"a\"] &, StringMatchQ[#, \"a\"] &, StringCount[#, \"o\"] &, "
                 "StringPosition[#, \"o\"] &, StringDelete[#, \"o\"] &, StringReplace[#, \"a\" -> \"A\"] &, "
                 "StringDrop[#, 1] &, StringTake[#, 1] &, StringTake[#, {1, -1, 2}] &, StringPart[#, 1] &, "
                 "StringPadLeft[#, 6] &, StringPadRight, StringTrim, StringExtract[#, 1] &, "
                 "Capitalize, Decapitalize, ToCharacterCode, StringJoin[#, #] &}", STR_FIXTURES);
}

void test_string_pattern_and_spec_lists() {
    /* a LinkedList of patterns, separators or specs where the List form takes a List */
    expect_llchk("{StringContainsQ[{\"abc\", \"xyz\"}, #] &, StringStartsQ[{\"abc\", \"xyz\"}, #] &, "
                 "StringEndsQ[\"abc\", #] &, StringFreeQ[\"abc\", #] &, StringMatchQ[\"b\", #] &, "
                 "StringCases[\"abcab\", #] &, StringCount[\"abcab\", #] &, StringPosition[\"abcab\", #] &, "
                 "StringDelete[\"abcab\", #] &, StringSplit[\"a b,c\", #] &, StringTrim[\"ab ba\", #] &, "
                 "StringMatchQ[\"xa\", \"x\" ~~ #] &, StringRiffle[{\"p\", \"q\"}, #] &}",
                 "{{\"a\", \"b\"}, {\" \", \",\"}, {\"(\", \"-\", \")\"}}");
    expect_llchk("{StringTake[\"abcdef\", #] &, StringDrop[\"abcdef\", #] &, StringPart[\"abcdef\", #] &, "
                 "StringExtract[\"a b c d\", #] &, FromCharacterCode}",
                 "{{2, 4}, {1, -1, 2}, {3}, {104, 105}}");
    expect_llchk("{EditDistance[#, {\"a\", \"b\"}] &, HammingDistance[#, Reverse[#]] &, EditDistance[#, #] &}",
                 STR_FIXTURES);
    /* code lists: a LinkedList of code LinkedLists, and a List of them */
    assert_eval_eq("{FromCharacterCode[LinkedList[{LinkedList[{104, 105}], {33}}]] // Normal, "
                   "FromCharacterCode[{LinkedList[{104, 105}], LinkedList[{33}]}]}",
                   "{{\"hi\", \"!\"}, {\"hi\", \"!\"}}", 0);
    /* a mutated pattern list is never served from the compiled-rule cache */
    assert_eval_eq("s9p1 = LinkedList[{\"a\"}]; s9r1 = StringCount[\"abab\", s9p1]; AppendTo[s9p1, \"b\"]; "
                   "{s9r1, StringCount[\"abab\", s9p1]}", "{2, 4}", 0);
}

/* ----------------------------------------------------------------- dates */

void test_dates() {
    expect_llchk("{DateList, AbsoluteTime, UnixTime}",
                 "{{2024, 3, 1}, {2024, 3, 1, 12, 30, 15.5}, {1999}, {3, 1, 2, 2}}");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_stats_vectors);
    TEST(test_stats_moving_and_scaling);
    TEST(test_stats_two_vectors);
    TEST(test_stats_matrices);
    TEST(test_stats_both_runs_bit_parity);
    TEST(test_string_subjects);
    TEST(test_string_pattern_and_spec_lists);
    TEST(test_dates);

    printf("All LinkedList strings/statistics tests passed!\n");
    symtab_clear();
    return 0;
}
