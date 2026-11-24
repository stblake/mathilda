/* Builtin-audit regressions: graphics (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* Plot3D#1 (high): a huge PlotPoints crashed the kernel -- build_grid's
 * n*n GridPt malloc wrapped (2^32 squared is 0 mod 2^64) and was written
 * through.  The grid is checked first: past memory aborts with General::nomem. */
static void test_plot3d_huge_plot_points(void) {
    audit_chk("Plot3D#1",
              "Quiet[CheckAbort[Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2^32]], aborted]]",
              "aborted");
    audit_chk("Plot3D#1",
              "Quiet[Check[CheckAbort[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2^31], 0], "
              "fired, General::nomem]]",
              "fired");
    audit_chk("Plot3D#1", "Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 10]]", "Graphics3D");
}

/* Histogram#1 (high): with an explicit {min, max, dw} bin spec, data outside the bins
 * was clamped into the edge bins (15 counted in [5, 10), 4 in [3, 4)).  The bins are
 * left-closed [min + k dw, min + (k+1) dw) and data in none of them is not counted. */
static void test_histogram_explicit_range_excludes(void) {
    audit_chk("Histogram#1",
              "Cases[Histogram[{1, 5, 9, 15}, {0, 10, 5}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity]",
              "{{0.0, 5.0, 1.0}, {5.0, 10.0, 2.0}}");
    audit_chk("Histogram#1",
              "Cases[Histogram[{1, 2, 3, 4}, {0, 4, 1}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity]",
              "{{1.0, 2.0, 1.0}, {2.0, 3.0, 1.0}, {3.0, 4.0, 1.0}}");
    audit_chk("Histogram#1",
              "Cases[Histogram[{-3, 0, 2, 2.5, 7}, {0, 5, 2.5}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity]",
              "{{0.0, 2.5, 2.0}, {2.5, 5.0, 1.0}}");
    /* the counted total is the number of data in [min, max) */
    audit_true("Histogram#1",
               "SeedRandom[9]; Module[{d = RandomReal[{-1, 2}, 500]},"
               " Total[Cases[Histogram[d, {0, 1, 0.1}], Rectangle[_, {_, h_}] :> h, Infinity]] =="
               " Count[d, v_ /; 0 <= v < 1]]");
    /* sanity: a derived range still counts every datum, the maximum included */
    audit_true("Histogram#1",
               "Total[Cases[Histogram[{1, 2, 3, 4}], Rectangle[_, {_, h_}] :> h, Infinity]] == 4");
    audit_true("Histogram#1",
               "Total[Cases[Histogram[{1, 2, 3, 4, 9}, 3], Rectangle[_, {_, h_}] :> h, Infinity]] == 5");
}

int main(void) {
    audit_init();
    TEST(test_plot3d_huge_plot_points);
    TEST(test_histogram_explicit_range_excludes);
    return audit_finish("audit_graphics");
}
