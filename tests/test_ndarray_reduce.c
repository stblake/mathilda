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
    /* {k} and {n1,n2} level specs are computed on the buffer, not degraded to
     * the List path. This assertion used to read `... == Total[list, {2}]` and
     * pass because the spec fell through and materialised; a contiguous axis
     * range is now summed in place, so the result stays an NDArray and is
     * compared the same way every other supported spec is. (Total[m, {1}] and
     * Total[m] are the same operation by definition and used to differ by 190x.) */
    snprintf(nd, sizeof(nd), "NDArray[%s]", m);
    chk_array("Total[%s, {1}]", nd, m);
    chk_array("Total[%s, {2}]", nd, m);
    chk_array("Total[%s, {1,2}]", nd, m);
    chk_eq("NDArrayQ[Total[NDArray[{{1.,2.},{3.,4.}}], {2}]]", "True");
    /* rank 3: every contiguous range of levels, against the List path */
    {
        const char* t = "{{{1.,2.},{3.,4.}},{{5.,6.},{7.,8.}},{{9.,10.},{11.,12.}}}";
        char nd3[256];
        snprintf(nd3, sizeof(nd3), "NDArray[%s]", t);
        chk_array("Total[%s, {1}]", nd3, t);
        chk_array("Total[%s, {2}]", nd3, t);
        chk_array("Total[%s, {3}]", nd3, t);
        chk_array("Total[%s, {1,2}]", nd3, t);
        chk_array("Total[%s, {2,3}]", nd3, t);
        chk_array("Total[%s, {2,Infinity}]", nd3, t);
        chk_scalar("Total[%s, {1,3}]", nd3, t);
    }
    /* Exactness survives: an int64 buffer keeps Integer heads, and a sum that
     * overflows int64 abandons the buffer rather than wrapping.
     *
     * ToNDArray, not NDArray[...]: the bare constructor defaults to float64, so
     * NDArray[{{1,2},{3,4}}] holds 1., 2., ... and answering Real to Head is
     * correct there (Total[..., 2] has always done so). Only the dtype-inferring
     * spelling gives an int64 buffer, which is the one with something to prove. */
    chk_eq("Head[Total[ToNDArray[{{1,2},{3,4}}], {2}][[1]]]", "Integer");
    chk_eq("Total[ToNDArray[{{4611686018427387904, 4611686018427387904},{1,2}}], {2}] == "
           "Total[{{4611686018427387904, 4611686018427387904},{1,2}}, {2}]", "True");
    /* Specs outside 1..rank still degrade faithfully to the List path. */
    chk_eq("Total[NDArray[{{1.,2.},{3.,4.}}], {0}] == Total[{{1.,2.},{3.,4.}}, {0}]", "True");
    chk_eq("Total[NDArray[{{1.,2.},{3.,4.}}], {3}] == Total[{{1.,2.},{3.,4.}}, {3}]", "True");
    chk_eq("Total[NDArray[{{1.,2.},{3.,4.}}], {-1}] == Total[{{1.,2.},{3.,4.}}, {-1}]", "True");
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
    /* Clip default [-1,1] and explicit bounds.
     *
     * Clip returns the BOUND at every clipped position, with the bound's own
     * head -- so on a Real array an EXACT bound produces exact Integers and the
     * result is not uniform:
     *
     *     Clip[{-2., 0., 2.}]           ->  {-1, 0., 1}     (default bounds are exact)
     *     Clip[{-2., 0., 2.}, {-1., 1.}] ->  {-1., 0., 1.}   (Real bounds)
     *
     * This used to assert NDArrayQ -> True for the FIRST of those, i.e. that
     * Clip answered {-1., 0., 1.} where the List path answered {-1, 0., 1}. The
     * chk_array checks beside it could not catch that: they compare numeric
     * distance, and the two differ only in their element HEADS. The invariant is
     * that the packed answer equals the List answer, so that is what is asserted
     * now -- through a comparison that is sensitive to exactness. */
    chk_array("Clip[%s]", "NDArray[{-2.,-0.5,0.5,2.}]", "{-2.,-0.5,0.5,2.}");
    chk_array("Clip[%s, {2.,8.}]", "NDArray[{1.,5.,10.}]", "{1.,5.,10.}");
    chk_eq("Clip[NDArray[{-2.,0.,2.}]] === Clip[{-2.,0.,2.}]", "True");
    chk_eq("Map[Head, Clip[NDArray[{-2.,0.,2.}]]] === {Integer, Real, Integer}", "True");
    chk_eq("Clip[NDArray[{-2.,0.,2.}], {-1,1}] === Clip[{-2.,0.,2.}, {-1,1}]", "True");
    /* Real bounds ARE uniform, so those stay on the buffer... */
    chk_eq("NDArrayQ[Clip[NDArray[{-2.,0.,2.}], {-1.,1.}]]", "True");
    /* Normal[] on the left because a VISIBLE NDArray[...] argument gives a
     * visible NDArray[...] result -- that presentation is the point of naming
     * the head, and it is not SameQ to a List. The VALUES are what must agree. */
    chk_eq("Normal[Clip[NDArray[{-2.,0.,2.}], {-1.,1.}]] === Clip[{-2.,0.,2.}, {-1.,1.}]", "True");
    /* ...and so does an exact bound that nothing reaches, where the answer is
     * just the input. */
    chk_eq("NDArrayQ[Clip[NDArray[{-2.,0.,2.}], {-99,99}]]", "True");
    chk_eq("Normal[Clip[NDArray[{-2.,0.,2.}], {-99,99}]] === Clip[{-2.,0.,2.}, {-99,99}]", "True");
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
