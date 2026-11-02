/* test_sparse_reduce.c -- native SparseArray reductions and comparisons
 * (src/sparse/sparse_reduce.c): Total, Mean, Variance, StandardDeviation, Tr,
 * Max, Min, MinMax, Norm, Accumulate, Equal, Unequal.
 *
 * Values are checked against fixed expected answers (and the answer KIND:
 * dense List, scalar, or SparseArray), and against the invariant
 * f[s] == f[Normal[s]] where the two must agree.
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

/* The InputForm of a (possibly SparseArray-valued) input. */
static void chk_if(const char* input, const char* expected) {
    char buf[4096];
    snprintf(buf, sizeof buf, "InputForm[%s]", input);
    chk(buf, expected);
}

static void chk_true(const char* input) { chk(input, "True"); }

/* ------------------------------------------------------------------------ */

static void setup(void) {
    chk("S = SparseArray[{{1,1}->1,{2,3}->2,{3,2}->-1,{3,3}->4},{3,3}]; "
        "S3 = SparseArray[{{1,1,1}->1,{2,3,1}->2,{2,2,2}->5},{2,3,2}]; "
        "S3b = SparseArray[{{1,1,1}->1,{2,3,1}->2,{2,2,2}->5},{2,3,2}, 1]; "
        "V = SparseArray[{1->3,4->-2},{5}]; "
        "R = SparseArray[{{1,1}->1.5,{2,1}->-2.5,{2,3}->0.25},{3,4}]; Null", "Null");
}

static void test_total(void) {
    chk("Total[S]", "{1, -1, 6}");
    chk("Total[S, {2}]", "{1, 2, 3}");
    chk("Total[S, 2]", "6");
    chk("Total[S, -1]", "6");
    chk("Total[S, Infinity]", "6");
    chk("Total[S, {1, 2}]", "6");
    chk("Total[V]", "1");
    chk("Head[Total[S]]", "List");
    /* rank >= 2 results over an exact 0 background stay sparse */
    chk_if("Total[S3]", "SparseArray[Automatic, {3, 2}, 0, {1, {{0, 1, 2, 3}, {{1}, {2}, {1}}}, {1, 5, 2}}]");
    chk_if("Total[S3, {2}]", "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 3}, {{1}, {1}, {2}}}, {1, 2, 5}}]");
    chk_if("Total[S3, {3}]", "SparseArray[Automatic, {2, 3}, 0, {1, {{0, 1, 3}, {{1}, {2}, {3}}}, {1, 5, 2}}]");
    chk("Total[S3, 2]", "{3, 5}");
    chk("Total[S3, {2, 3}]", "{1, 7}");
    chk("Total[S3, -1]", "8");
    chk_true("Total[S3, 0] === S3");
    /* ... but dense over any other background */
    chk("Total[S3b]", "{{2, 2}, {2, 6}, {3, 2}}");
    chk("Total[S3b, {1, 3}]", "17");
    chk("Total[S3b, {1, 2}]", "{7, 10}");
    /* explicit cancelling entries are not re-minimised */
    chk_if("Total[SparseArray[{{1,1,1}->1, {2,1,1}->-1},{2,2,2}]]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 1}, {{1}}}, {0}}]");
    chk_if("Total[SparseArray[{{1,1,1}->1.5},{2,2,2}]]",
           "SparseArray[Automatic, {2, 2}, 0.0, {1, {{0, 1, 1}, {{1}}}, {1.5}}]");
    /* exactness: exact where every summand is exact */
    chk_if("Total[SparseArray[{{1,1}->1.5},{2,2},3]]", "{4.5, 6}");
    chk_if("Total[SparseArray[{{1,1}->1, {2,1}->2},{2,2},0.5]]", "{3, 1.0}");
    chk_if("Total[SparseArray[{{1,1}->1.0},{2,2}]]", "{1.0, 0.0}");
    chk_if("Total[SparseArray[{{1,1}->1},{2,2},0.0]]", "{1.0, 0.0}");
    chk_if("Total[SparseArray[{{1,1}->1/2},{2,2}]]", "{1/2, 0}");
    chk("Total[SparseArray[{{1,1}->1},{2,2},x]]", "{1 + x, 2 x}");
    chk("Total[SparseArray[{{1,1,1}->1, {1,1,2}->2},{2,2,2}, x]]", "{{1 + x, 2 + x}, {2 x, 2 x}}");
    /* symbolic leaves are leaves for negative / Infinity levels */
    chk("Total[SparseArray[{1->b, 2->a+1},{2}], {-1}]", "1 + a + b");
    chk("Total[SparseArray[{z^2, 1}], Infinity]", "1 + z^2");
    chk("Total[SparseArray[{{z^2, 2},{3,4}}], {-2}]", "{3 + z^2, 6}");
    /* ... as on the dense form (Total counts List depth only) */
    chk("Total[{z^2, 1}, Infinity]", "1 + z^2");
    chk("Total[{{z^2, 2}, {3, 4}}, {-2}]", "{3 + z^2, 6}");
    chk("Total[{f[1, 2], 3}, -1]", "3 + f[1, 2]");
    /* integer overflow falls back to exact big sums */
    chk("Total[SparseArray[{1 -> 2^62, 2 -> 2^62, 3 -> 2^62}]]", "13835058055282163712");
    chk("Total[SparseArray[{{1,1} -> 2^62, {2,1} -> 2^62}, {2, 2}]]", "{9223372036854775808, 0}");
    /* agreement with the dense form */
    chk_true("Total[R] == Total[Normal[R]] && Total[R, {2}] == Total[Normal[R], {2}]");
    chk_true("Total[R, 2] == Total[Normal[R], 2]");
}

static void test_mean_variance(void) {
    chk("Mean[S]", "{1/3, -1/3, 2}");
    chk("Mean[V]", "1/5");
    chk_if("Mean[S3]", "SparseArray[Automatic, {3, 2}, 0, {1, {{0, 1, 2, 3}, {{1}, {2}, {1}}}, {1/2, 5/2, 1}}]");
    chk("Mean[S3b]", "{{1, 1}, {1, 3}, {3/2, 1}}");
    chk_if("Mean[SparseArray[{{1,1}->1.5},{2,2},3]]", "{2.25, 3}");
    chk("Variance[V]", "16/5");
    chk("Variance[S]", "{1/3, 1/3, 4}");
    chk("StandardDeviation[S]", "{1/Sqrt[3], 1/Sqrt[3], 2}");
    chk_if("Variance[SparseArray[{1->1.5},4]]", "0.5625");
    chk_if("Variance[SparseArray[{{1,1}->1,{2,2}->1.0},{3,2}]]", "{1/3, 0.33333333333333337}");
    chk("StandardDeviation[SparseArray[{1->1},4]]", "1/2");
    chk_true("Abs[StandardDeviation[N[V]] - StandardDeviation[Normal[N[V]]]] < 10^-15");
}

static void test_extremes(void) {
    chk("Max[S]", "4");
    chk("Min[S]", "-1");
    chk("Min[V]", "-2");
    chk("MinMax[S]", "{-1, 4}");
    chk("Max[S, 10]", "10");
    chk("Max[S, y]", "Max[4, y]");
    chk("Max[S, {7, 8}]", "8");
    chk("Max[S, S3]", "5");
    chk("Max[SparseArray[{{1,1}->-1},{1,1}]]", "-1");
    chk("Max[SparseArray[{1->x},2]]", "Max[0, x]");
    chk("Max[SparseArray[{},{3}]]", "0");
    chk_if("Max[SparseArray[{1->2.0},3]]", "2.0");
    chk_if("MinMax[SparseArray[{1->-2.0},3]]", "{-2.0, 0}");
    chk_if("Max[SparseArray[{1->-2.0},3]]", "0");
    chk("Min[SparseArray[{1->-2, 2->3, 3->1},3]]", "-2");
}

static void test_tr(void) {
    chk("Tr[S]", "5");
    chk("Tr[S, List]", "{1, 0, 4}");
    chk("Tr[S3]", "6");
    chk("Tr[SparseArray[{{1,1}->1},{2,3},x]]", "1 + x");
    chk_if("Tr[SparseArray[{{1,1}->1.5},{3,3},2]]", "5.5");
    chk("Tr[V]", "1");
}

static void test_norm(void) {
    chk("Norm[V]", "Sqrt[13]");
    chk("Norm[V, 1]", "5");
    chk("Norm[V, Infinity]", "3");
    chk("Norm[V, 3]", "35^(1/3)");
    chk_true("Abs[Norm[N[V], 3] - 3.2710663101885897] < 10^-14");
    chk("Norm[SparseArray[{1->3, 2->4},5, 1]]", "2 Sqrt[7]");
    chk_if("Norm[SparseArray[{1->3, 2->4},3], 2.0]", "5.0");
    chk("Norm[S, 1]", "6");
    chk("Norm[S, Infinity]", "5");
    chk("Norm[S, \"Frobenius\"]", "Sqrt[22]");
    chk("Norm[SparseArray[{{1,1}->3, {2,1}->-4},{2,2}, 1], 1]", "7");
    chk("Norm[SparseArray[{{1,1}->3, {2,1}->-4},{2,2}, 1], \"Frobenius\"]", "3 Sqrt[3]");
    chk_if("Norm[SparseArray[{{1,1}->3.0, {2,1}->-4},{2,2}], Infinity]", "4.0");
    /* the dense induced 1- and Infinity-norms of exact matrices agree */
    chk("Norm[{{1, 2}, {-3, 4}}, 1]", "6");
    chk("Norm[{{1, 2}, {-3, 4}}, Infinity]", "7");
    chk("Norm[{{0}}, 1]", "0");
    chk("Norm[{{a, 0}, {0, b}}, 1]", "Max[Abs[a], Abs[b]]");
    /* the machine 2-norm: Lanczos agrees with the dense SVD */
    chk_true("Abs[Norm[N[S]] - Norm[Normal[N[S]]]] < 10^-14");
    chk_true("SeedRandom[5]; M = SparseArray[RandomInteger[{1, 40}, {300, 2}] -> RandomReal[{-1, 1}, 300], {40, 60}]; "
             "Abs[Norm[M] - Norm[Normal[M]]] < 10^-13 Norm[M]");
    chk_true("Abs[Norm[Transpose[M]] - Norm[Normal[M]]] < 10^-13 Norm[M]");
    chk_if("Norm[SparseArray[{}, {3, 3}, 0.0]]", "0.0");
    /* an exact matrix 2-norm is not the machine kernel's to answer */
    chk("Head[Norm[S]]", "Norm");
}

static void test_accumulate(void) {
    chk("Accumulate[V]", "{3, 3, 3, 1, 1}");
    chk("Accumulate[S]", "{{1, 0, 0}, {1, 0, 2}, {1, -1, 6}}");
    chk("Accumulate[S3]", "{{{1, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 5}, {2, 0}}}");
    chk_if("Accumulate[SparseArray[{1->1.0},3]]", "{1.0, 1.0, 1.0}");
    chk_if("Accumulate[SparseArray[{2->1.5},3, 0.5]]", "{0.5, 2.0, 2.5}");
    chk("Accumulate[SparseArray[{2->5},3, 2]]", "{2, 7, 9}");
    chk("Head[Accumulate[V]]", "List");
    chk_true("Accumulate[N[S3b]] == Accumulate[Normal[N[S3b]]]");
}

static void test_equal(void) {
    chk_true("S == Normal[S]");
    chk_true("S == S");
    chk("S != S", "False");
    chk_true("S == N[S]");
    chk("S == Transpose[S]", "False");
    chk("S == {1, 2}", "False");
    chk("S == S3", "False");
    chk("SparseArray[{1->a},2] == SparseArray[{1->a},2]", "True");
    chk("Head[SparseArray[{1->a},2] == SparseArray[{1->b},2]]", "Equal");
    chk_true("SparseArray[{1->1},2] == SparseArray[{1->1.0},2]");
    chk_true("SparseArray[{1->1},2] == SparseArray[{2->0},2, 1]");
    chk_true("Equal[SparseArray[{1->1},2], SparseArray[{2->0},2, 1], {1, 0}]");
    chk_true("SparseArray[{1->1},2] != SparseArray[{1->2},2]");
    chk("SparseArray[{1->1},2] != {1, 0}", "False");
    chk("Head[SparseArray[{1->1},2] == 1]", "Equal");
    chk("SparseArray[{1->1},2] == {1, 0, 0}", "False");
    chk("SparseArray[{1->a},2] == {b, 1}", "False");
    chk_true("SparseArray[{1 -> 0.1 + 0.2}, 3] == {0.3, 0, 0}");
    chk("SparseArray[{1 -> 1.0}, 2] == SparseArray[{1 -> 1.001}, 2]", "False");
    /* a packed dense operand */
    chk_true("SparseArray[{5 -> 2}, 300] == ReplacePart[ConstantArray[0, 300], 5 -> 2]");
    chk("SparseArray[{5 -> 2}, 300] == ConstantArray[0, 300]", "False");
}

static void test_large(void) {
    /* sizes the dense gate refuses: answered natively */
    chk("L = SparseArray[{{1, 1} -> 2.0, {7, 99999} -> -1.0, {100000, 3} -> 5.0}, {10^5, 10^5}]; Null", "Null");
    chk("Total[L, 2]", "6.0");
    chk("Total[L][[{1, 3, 99999}]]", "{2.0, 5.0, -1.0}");
    chk("Total[L, {2}][[{1, 7, 100000}]]", "{2.0, -1.0, 5.0}");
    chk("Max[L]", "5.0");
    chk("Min[L]", "-1.0");
    chk("Norm[L, \"Frobenius\"]", "5.47723");
    chk("Norm[L]", "5.0");
    chk("Norm[L, 1]", "5.0");
    chk("Tr[L]", "2.0");
    chk_true("L == L");
    chk("L == L + 1.0", "False");
    chk_true("Mean[L][[3]] == 0.00005");
    chk("Accumulate[SparseArray[{10 -> 1}, 10^6]][[{9, 10, 10^6}]]", "{0, 1, 1}");
}

/* Arrays large enough for the partitioned (threaded) kernels. */
static void test_parallel(void) {
    chk("SeedRandom[11]; P = SparseArray[Transpose[{RandomInteger[{1, 2000}, 10^6], "
        "RandomInteger[{1, 400}, 10^6]}] -> RandomInteger[{-9, 9}, 10^6], {2000, 400}]; "
        "PN = Normal[P]; Length[P[\"NonzeroValues\"]] > 2^18", "True");
    chk_true("Total[P] === Total[PN]");
    chk_true("Total[P, {2}] === Total[PN, {2}]");
    chk_true("Total[P, 2] === Total[PN, 2]");
    chk_true("Max[P] === Max[PN] && Min[P] === Min[PN]");
    chk_true("Total[P + 7] === Total[PN + 7]");
    chk_true("Max[Abs[Total[N[P] + 0.5] - Total[N[PN] + 0.5]]] < 10^-9");
    chk_true("Max[Abs[Total[N[P], {2}] - Total[N[PN], {2}]]] < 10^-9");
    chk_true("Abs[Total[N[P], 2] - Total[N[PN], 2]] < 10^-8");
    chk_true("Max[N[P]] == Max[N[PN]] && Min[N[P]] == Min[N[PN]]");
    chk_true("Abs[Norm[N[P], \"Frobenius\"] - Norm[N[PN], \"Frobenius\"]] < 10^-9");
    chk_true("Norm[P, 1] === Norm[PN, 1] && Norm[P, Infinity] === Norm[PN, Infinity]");
    chk_true("Mean[P] === Mean[PN]");
    chk_true("P == PN && P == P");
}

int main(void) {
    symtab_init();
    core_init();
    setup();
    test_total();
    test_mean_variance();
    test_extremes();
    test_tr();
    test_norm();
    test_accumulate();
    test_equal();
    test_large();
    test_parallel();
    if (failures) {
        fprintf(stderr, "\n%d SparseArray reduction test(s) FAILED\n", failures);
        return 1;
    }
    printf("All SparseArray reduction tests passed.\n");
    return 0;
}
