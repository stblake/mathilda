/* Pure-numeric tests for src/graphics/sampling.c -- no Expr/symtab/Raylib
 * involved, so these run headless on any build (USE_GRAPHICS on or off). */
#include "sampling.h"
#include "test_utils.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static bool sin_fn(double x, void* ctx, double* y) { (void)ctx; *y = sin(x); return true; }
static bool inv_fn(double x, void* ctx, double* y) {
    (void)ctx;
    if (fabs(x) < 1e-9) return false;
    *y = 1.0 / x;
    return true;
}
static bool sharp_fn(double x, void* ctx, double* y) { (void)ctx; *y = (x < 0) ? -1.0 : 1.0; return true; }

void test_sampling_no_recursion_returns_exact_grid(void) {
    size_t n;
    PlotPoint* pts = plot_sample_adaptive(sin_fn, NULL, 0.0, 1.0, 5, 0, -1, &n);
    ASSERT(pts != NULL);
    ASSERT_MSG(n == 5, "expected exactly 5 points with MaxRecursion=0, got %zu", n);
    plot_points_free(pts);
}

void test_sampling_refines_curved_function(void) {
    size_t n;
    PlotPoint* pts = plot_sample_adaptive(sin_fn, NULL, 0.0, 2 * M_PI, 10, 6, -1, &n);
    ASSERT(pts != NULL);
    ASSERT_MSG(n > 10, "expected adaptive refinement beyond the initial PlotPoints grid, got %zu", n);
    for (size_t i = 0; i < n; i++) ASSERT(pts[i].valid);
    plot_points_free(pts);
}

void test_sampling_detects_singularity(void) {
    size_t n;
    PlotPoint* pts = plot_sample_adaptive(inv_fn, NULL, -1.0, 1.0, 9, 6, -1, &n);
    ASSERT(pts != NULL);
    bool found_break = false;
    for (size_t i = 0; i < n; i++) if (pts[i].break_before) found_break = true;
    ASSERT_MSG(found_break, "expected a break_before around the 1/x singularity at x=0");
    plot_points_free(pts);
}

void test_sampling_respects_max_plot_points(void) {
    size_t n;
    PlotPoint* pts = plot_sample_adaptive(sharp_fn, NULL, -1.0, 1.0, 5, 12, 20, &n);
    ASSERT(pts != NULL);
    ASSERT_MSG(n <= 30, "expected point count to stay near the MaxPlotPoints budget (20), got %zu", n);
    plot_points_free(pts);
}

void test_sampling_rejects_degenerate_range(void) {
    size_t n;
    PlotPoint* pts = plot_sample_adaptive(sin_fn, NULL, 1.0, 1.0, 10, 6, -1, &n);
    ASSERT(pts == NULL);
    ASSERT(n == 0);
}

void test_sampling_rejects_too_few_plot_points(void) {
    size_t n;
    PlotPoint* pts = plot_sample_adaptive(sin_fn, NULL, 0.0, 1.0, 1, 6, -1, &n);
    ASSERT(pts == NULL);
    ASSERT(n == 0);
}

/* ---- plot_robust_yrange: spike-resistant vertical auto-range ---- */

void test_yrange_clips_asymptote_spikes(void) {
    /* Tan-like: a dense bulk in [-1, 1] plus a few enormous spikes (the
     * sampler's climb toward an asymptote). The band must ignore the spikes. */
    double ys[40];
    for (int i = 0; i < 36; i++) ys[i] = -1.0 + 2.0 * (i / 35.0);
    ys[36] = 1e6; ys[37] = -1e6; ys[38] = 3e5; ys[39] = -2e5;
    double lo, hi;
    plot_robust_yrange(ys, 40, &lo, &hi);
    ASSERT_MSG(hi < 100.0 && lo > -100.0,
               "expected spikes clipped to a small band, got [%g, %g]", lo, hi);
    ASSERT_MSG(hi > 0.5 && lo < -0.5, "expected the bulk kept, got [%g, %g]", lo, hi);
}

void test_yrange_keeps_legitimate_extrema(void) {
    /* A smooth ramp 0..9 with no spikes -- nothing should be clipped. */
    double ys[30];
    for (int i = 0; i < 30; i++) ys[i] = 9.0 * (i / 29.0);
    double lo, hi;
    plot_robust_yrange(ys, 30, &lo, &hi);
    ASSERT_MSG(lo == 0.0 && hi == 9.0, "expected [0, 9] unchanged, got [%g, %g]", lo, hi);
}

void test_yrange_few_points_uses_full_extent(void) {
    /* Under 8 points: too little to estimate a robust spread, show all. */
    double ys[5] = { -1.0, 0.0, 2.0, 1e6, 3.0 };
    double lo, hi;
    plot_robust_yrange(ys, 5, &lo, &hi);
    ASSERT_MSG(lo == -1.0 && hi == 1e6, "expected exact min/max fallback, got [%g, %g]", lo, hi);
}

void test_yrange_symmetric_spikes_clipped_both_sides(void) {
    /* 1/x-like: bulk near 0, large spikes both directions -- both clamped and
     * the band stays roughly symmetric and small. */
    double ys[40];
    for (int i = 0; i < 32; i++) ys[i] = -0.5 + 1.0 * (i / 31.0);
    ys[32] = 1e6; ys[33] = -1e6; ys[34] = 5e5; ys[35] = -5e5;
    ys[36] = 8e5; ys[37] = -8e5; ys[38] = 2e5; ys[39] = -2e5;
    double lo, hi;
    plot_robust_yrange(ys, 40, &lo, &hi);
    ASSERT_MSG(hi < 100.0 && lo > -100.0,
               "expected both sides clipped, got [%g, %g]", lo, hi);
    ASSERT_MSG(fabs(hi + lo) < 0.5 * (hi - lo),
               "expected a roughly symmetric band, got [%g, %g]", lo, hi);
}

int main(void) {
    TEST(test_sampling_no_recursion_returns_exact_grid);
    TEST(test_sampling_refines_curved_function);
    TEST(test_sampling_detects_singularity);
    TEST(test_sampling_respects_max_plot_points);
    TEST(test_sampling_rejects_degenerate_range);
    TEST(test_sampling_rejects_too_few_plot_points);
    TEST(test_yrange_clips_asymptote_spikes);
    TEST(test_yrange_keeps_legitimate_extrema);
    TEST(test_yrange_few_points_uses_full_extent);
    TEST(test_yrange_symmetric_spikes_clipped_both_sides);
    printf("All graphics sampling tests passed!\n");
    return 0;
}
