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

int main(void) {
    TEST(test_sampling_no_recursion_returns_exact_grid);
    TEST(test_sampling_refines_curved_function);
    TEST(test_sampling_detects_singularity);
    TEST(test_sampling_respects_max_plot_points);
    TEST(test_sampling_rejects_degenerate_range);
    TEST(test_sampling_rejects_too_few_plot_points);
    printf("All graphics sampling tests passed!\n");
    return 0;
}
