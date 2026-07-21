/* test_ndarray_reduce.c — reduction, order-statistic, structural and moving-
 * statistic fast paths over NDArray objects (src/ndreduce.c, src/ndstruct.c).
 *
 * The correctness property throughout is "matches the List path": for machine
 * data D, f[NDArray[D]] must equal f[D]. NDArray-returning results are compared
 * elementwise after Normal via
 *     Max[Abs[Flatten[Normal[f[NDArray[D]]] - f[D]]]] < eps
 * (robust to last-ULP differences between the C buffer loops and the scalar
 * builtins, and to the exact-int vs real distinction the List path sometimes
 * keeps); scalar results via Abs[... - ...] < eps. Structural results (result
 * stays a packed NDArray) and the faithful degrade at unsupported specs are
 * asserted directly.
 *
 * exit(1) on failure (the CMake Release build passes -DNDEBUG, no-oping assert). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

static int failures = 0;

/* Evaluate `input`; fail unless it prints exactly `expected`. */
static void chk_eq(const char* input, const char* expected) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n",
                input, expected, s);
        failures++;
    }
    free(s);
    expr_free(e);
}

/* Assert f[NDArray[data]] matches f[data] (elementwise, within 1e-9). `fn` is
 * the whole call with %s placeholders for the argument spelling. */
static void chk_scalar(const char* call_fmt, const char* nd, const char* lst) {
    char callnd[512], calll[512], buf[1152];
    snprintf(callnd, sizeof(callnd), call_fmt, nd);
    snprintf(calll, sizeof(calll), call_fmt, lst);
    snprintf(buf, sizeof(buf), "Abs[(%s) - (%s)] < 1/1000000000", callnd, calll);
    chk_eq(buf, "True");
}

/* Assert Normal[f[NDArray[data]]] matches f[data] elementwise within eps. */
static void chk_array(const char* call_fmt, const char* nd, const char* lst) {
    char callnd[512], calll[512], buf[1280];
    snprintf(callnd, sizeof(callnd), call_fmt, nd);
    snprintf(calll, sizeof(calll), call_fmt, lst);
    snprintf(buf, sizeof(buf),
             "Max[Abs[Flatten[Normal[%s] - (%s)]]] < 1/1000000000", callnd, calll);
    chk_eq(buf, "True");
}

/* ------------------------------------------------------------- reductions */

static void test_total(void) {
    const char* v = "{1.5, 2.5, 3.5, 4.5}";
    const char* m = "{{1.,2.,3.},{4.,5.,6.}}";
    char nd[128], lst[128];
    /* vector -> scalar */
    snprintf(nd, sizeof(nd), "NDArray[%s]", v);
    chk_scalar("Total[%s]", nd, v);
    /* matrix default (leading axis) -> vector */
    snprintf(nd, sizeof(nd), "NDArray[%s]", m);
    chk_array("Total[%s]", nd, m);
    chk_eq("NDArrayQ[Total[NDArray[{{1.,2.},{3.,4.}}]]]", "True");
    /* Total[m, 2] and Total[m, Infinity] -> scalar full sum */
    chk_eq("Total[NDArray[{{1.,2.},{3.,4.}}], 2] == Total[{{1.,2.},{3.,4.}}, 2]", "True");
    chk_eq("Total[NDArray[{{1.,2.},{3.,4.}}], Infinity] == 10.", "True");
    /* {k} level spec degrades to the exact List result */
    chk_eq("Total[NDArray[{{1.,2.},{3.,4.}}], {2}] == Total[{{1.,2.},{3.,4.}}, {2}]", "True");
}

static void test_mean_variance(void) {
    const char* v = "{2., 4., 4., 4., 5., 5., 7., 9.}";
    const char* m = "{{1.,2.},{3.,4.},{5.,9.}}";
    char nd[128];
    snprintf(nd, sizeof(nd), "NDArray[%s]", v);
    chk_scalar("Mean[%s]", nd, v);
    chk_scalar("Variance[%s]", nd, v);
    chk_scalar("StandardDeviation[%s]", nd, v);
    chk_scalar("RootMeanSquare[%s]", nd, v);
    /* matrix columnwise */
    snprintf(nd, sizeof(nd), "NDArray[%s]", m);
    chk_array("Mean[%s]", nd, m);
    chk_array("Variance[%s]", nd, m);
    chk_array("StandardDeviation[%s]", nd, m);
    /* complex Variance yields a real */
    chk_eq("Variance[NDArray[{1.+2.*I,3.-1.*I,0.+0.*I}, DataType->\"complex64\"]] == "
           "Variance[{1.+2.*I,3.-1.*I,0.+0.*I}]", "True");
    /* n < 2 degrades (List gives the symbolic/edge result, so no crash) */
    chk_eq("Variance[NDArray[{5.}]] == Variance[{5.}]", "True");
}

static void test_maxmin_accumulate(void) {
    const char* v = "{3., 1., 9., 2., 7.}";
    char nd[128];
    snprintf(nd, sizeof(nd), "NDArray[%s]", v);
    chk_scalar("Max[%s]", nd, v);
    chk_scalar("Min[%s]", nd, v);
    /* Max/Min flatten a matrix fully */
    chk_eq("Max[NDArray[{{1.,9.},{3.,4.}}]] == 9.", "True");
    chk_array("Accumulate[%s]", nd, v);
    chk_array("Accumulate[%s]", "NDArray[{{1.,2.},{3.,4.},{5.,6.}}]",
              "{{1.,2.},{3.,4.},{5.,6.}}");
    /* Max on complex degrades (stays symbolic, matching the List path) */
    chk_eq("Max[NDArray[{1.+2.*I,3.+4.*I}, DataType->\"complex64\"]] == "
           "Max[{1.+2.*I,3.+4.*I}]", "True");
}

static void test_order_stats(void) {
    chk_scalar("Median[%s]", "NDArray[{3.,1.,2.,5.,4.}]", "{3.,1.,2.,5.,4.}");   /* odd */
    chk_scalar("Median[%s]", "NDArray[{3.,1.,2.,4.}]", "{3.,1.,2.,4.}");         /* even */
    chk_array("Median[%s]", "NDArray[{{1.,2.},{3.,4.},{5.,9.}}]",
              "{{1.,2.},{3.,4.},{5.,9.}}");                                       /* columnwise */
    chk_array("Quartiles[%s]", "NDArray[Range[1.,10.]]", "Range[1.,10.]");
    chk_array("Quartiles[%s]", "NDArray[Range[1.,20.]]", "Range[1.,20.]");
}

/* ------------------------------------------------------------- structural */

static void test_structural(void) {
    chk_array("Sort[%s]", "NDArray[{3.,1.,2.,5.,4.}]", "{3.,1.,2.,5.,4.}");
    chk_eq("NDArrayQ[Sort[NDArray[{3.,1.,2.}]]]", "True");
    chk_eq("DataType[Sort[NDArray[{3.,1.,2.}, DataType->\"float32\"]]]", "\"float32\"");
    chk_array("Reverse[%s]", "NDArray[{1.,2.,3.,4.}]", "{1.,2.,3.,4.}");
    chk_array("Reverse[%s]", "NDArray[{{1.,2.},{3.,4.},{5.,6.}}]",
              "{{1.,2.},{3.,4.},{5.,6.}}");
    chk_array("Transpose[%s]", "NDArray[{{1.,2.,3.},{4.,5.,6.}}]",
              "{{1.,2.,3.},{4.,5.,6.}}");
    chk_eq("Dimensions[Transpose[NDArray[{{1.,2.,3.},{4.,5.,6.}}]]] == {3, 2}", "True");
    chk_array("Flatten[%s]", "NDArray[{{1.,2.},{3.,4.}}]", "{{1.,2.},{3.,4.}}");
    /* Take / Drop leading-axis specs */
    chk_array("Take[%s, 3]", "NDArray[{1.,2.,3.,4.,5.}]", "{1.,2.,3.,4.,5.}");
    chk_array("Take[%s, -2]", "NDArray[{1.,2.,3.,4.,5.}]", "{1.,2.,3.,4.,5.}");
    chk_array("Take[%s, {2,4}]", "NDArray[{1.,2.,3.,4.,5.}]", "{1.,2.,3.,4.,5.}");
    chk_array("Drop[%s, 2]", "NDArray[{1.,2.,3.,4.,5.}]", "{1.,2.,3.,4.,5.}");
    chk_array("Drop[%s, -1]", "NDArray[{1.,2.,3.,4.,5.}]", "{1.,2.,3.,4.,5.}");
    chk_array("Drop[%s, {2,4}]", "NDArray[{1.,2.,3.,4.,5.}]", "{1.,2.,3.,4.,5.}");
    /* Clip default [-1,1] and explicit bounds */
    chk_array("Clip[%s]", "NDArray[{-2.,-0.5,0.5,2.}]", "{-2.,-0.5,0.5,2.}");
    chk_array("Clip[%s, {2.,8.}]", "NDArray[{1.,5.,10.}]", "{1.,5.,10.}");
    chk_eq("NDArrayQ[Clip[NDArray[{-2.,0.,2.}]]]", "True");
}

/* ------------------------------------------------------------ moving stats */

static void test_moving(void) {
    chk_array("MovingAverage[%s, 3]", "NDArray[{1.,2.,3.,4.,5.}]", "{1.,2.,3.,4.,5.}");
    chk_array("MovingMedian[%s, 3]", "NDArray[{5.,1.,3.,2.,4.,6.}]", "{5.,1.,3.,2.,4.,6.}");
    chk_array("ExponentialMovingAverage[%s, 0.5]", "NDArray[{1.,2.,3.,4.}]",
              "{1.,2.,3.,4.}");
    chk_eq("NDArrayQ[MovingAverage[NDArray[{1.,2.,3.,4.,5.}], 2]]", "True");
}

int main(void) {
    symtab_init();
    core_init();

    test_total();
    test_mean_variance();
    test_maxmin_accumulate();
    test_order_stats();
    test_structural();
    test_moving();

    if (failures) {
        fprintf(stderr, "\n%d NDArray reduction/structural test(s) FAILED\n", failures);
        return 1;
    }
    printf("All NDArray reduction/structural tests passed.\n");
    return 0;
}
