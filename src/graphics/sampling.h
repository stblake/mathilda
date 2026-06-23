/* sampling.h — pure-numeric adaptive curve sampler.
 *
 * No Expr/symtab/Raylib dependency by design: this is the one algorithm
 * `Plot` (and, later, `ListPlot`/`ParametricPlot`) all sample through, so
 * PlotPoints/MaxRecursion/MaxPlotPoints behave identically everywhere they
 * apply. Being pure C-doubles-in, C-doubles-out also makes it unit
 * testable with plain functions like sin() with no display required. */
#ifndef MATHILDA_GRAPHICS_SAMPLING_H
#define MATHILDA_GRAPHICS_SAMPLING_H

#include <stdbool.h>
#include <stddef.h>

/* Evaluate the sampled function at `x`. Returns false when the result is
 * not a finite real (singularity, complex, undefined) -- the sampler
 * treats that x as a gap rather than erroring. */
typedef bool (*PlotSampleFn)(double x, void* ctx, double* y_out);

typedef struct {
    double x, y;
    bool   valid;        /* always true for stored points; gaps are dropped, not stored */
    bool   break_before; /* true if this point starts a new segment after a gap */
} PlotPoint;

/*
 * Adaptively sample `fn` across [xmin, xmax].
 *
 *   plot_points    initial uniform grid size (>= 2; Mathematica's PlotPoints)
 *   max_recursion  max bisection depth per initial interval (MaxRecursion)
 *   max_plot_points  overall cap on stored points; <= 0 means unbounded
 *                    (MaxPlotPoints)
 *
 * Returns a malloc'd array of *out_count points (caller frees with
 * plot_points_free), or NULL if xmin >= xmax or plot_points < 2.
 */
PlotPoint* plot_sample_adaptive(PlotSampleFn fn, void* ctx,
                                 double xmin, double xmax,
                                 long plot_points, int max_recursion,
                                 long max_plot_points, size_t* out_count);

void plot_points_free(PlotPoint* pts);

/*
 * Robust [lo, hi] band for the vertical axis of a plot.
 *
 * Functions with asymptotes (Tan, 1/x, ...) produce a sparse heavy tail of
 * finite-but-huge sampled y-values as the adaptive sampler climbs toward each
 * singularity. Taking the raw min/max would let those spikes dominate the
 * range and flatten the visible curve. This computes a median/MAD fence and
 * intersects it with the true data extent: genuine extrema (parabola peaks,
 * exponential growth) sit inside the fence and survive untouched, while sparse
 * spikes are clamped to the fence so they run off-screen instead.
 *
 * `ys` is read-only and may be in any order. Falls back to the exact
 * [min, max] when n < 8 or the bulk is flat (MAD == 0) -- too little signal to
 * clip safely. With n == 0, returns lo > hi (an empty range) so the caller can
 * detect "nothing to bound".
 */
void plot_robust_yrange(const double* ys, size_t n, double* out_lo, double* out_hi);

#endif /* MATHILDA_GRAPHICS_SAMPLING_H */
