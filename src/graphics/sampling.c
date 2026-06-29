/* sampling.c — see sampling.h. */

#include "sampling.h"
#include <math.h>
#include <stdlib.h>

typedef struct {
    PlotPoint* pts;
    size_t     len;
    size_t     cap;
    long       max_points;     /* <= 0 means unbounded */
    bool       pending_break;  /* a gap is open; next pushed point gets break_before */
} SampleBuf;

static void buf_init(SampleBuf* buf, long max_points) {
    buf->cap = 64;
    buf->pts = malloc(sizeof(PlotPoint) * buf->cap);
    buf->len = 0;
    buf->max_points = max_points;
    buf->pending_break = false;
}

static void buf_push(SampleBuf* buf, double x, double y) {
    if (buf->len == buf->cap) {
        buf->cap *= 2;
        buf->pts = realloc(buf->pts, sizeof(PlotPoint) * buf->cap);
    }
    buf->pts[buf->len].x = x;
    buf->pts[buf->len].y = y;
    buf->pts[buf->len].valid = true;
    buf->pts[buf->len].break_before = buf->pending_break;
    buf->pending_break = false;
    buf->len++;
}

static bool buf_over_budget(const SampleBuf* buf) {
    return buf->max_points > 0 && (long)buf->len >= buf->max_points;
}

/* Maximum vertical deviation of the curve from the straight chord a->b that
 * the renderer would draw, as a fraction of the displayed y-extent. Below this
 * the segment is accepted as visually flat. This is the PRIMARY refinement
 * driver: because deviation grows with the segment's curvature (the sagitta of
 * a chord is ~|f''|*h^2/8), a threshold this small makes the sample density
 * track curvature -- dense where the curve bends (the peaks of Sin, the ends of
 * x^3, the body of a Runge spike), sparse where it is straight (the steep but
 * locally linear zero-crossings of Sin, a sloped line). At ~1/1666 of the
 * y-range it sits comfortably under a pixel on a normal window, so accepted
 * segments read as smooth while the density stays curvature-proportional.
 * (A looser value lets ordinary smooth curves pass as "flat" at the initial
 * grid scale, leaving the chord cap below as the only refiner -- which keys on
 * slope, not curvature, and inverts the density.) */
#define FLAT_TOL 0.0006

/* Upper bound on the screen-space length of any accepted leaf segment, as a
 * fraction of the displayed frame (Euclidean over the x- and y-normalized
 * axes). This is a LOOSE BACKSTOP, not the primary driver: the deviation test
 * above keys on curvature, but it is blind to *steepness*, so where a curve is
 * steep but locally straight (e.g. Log[x] near 0, or any near-vertical approach
 * to an asymptote) the chord can be sub-pixel-accurate yet the sample points --
 * and the Mesh dots drawn at them -- spread far apart, thinning abruptly at a
 * grid boundary. Capping the chord length bounds that gap in BOTH axes at once.
 * It must stay generous (a large fraction of the frame): tightening it makes the
 * cap out-compete the curvature test on ordinary moderately-sloped curves and
 * pour samples into steep-but-straight stretches (Sin's zero-crossings) instead
 * of the curvy parts (Sin's peaks), inverting the intended density. At ~1/12 of
 * the frame it only ever binds where a single segment would otherwise leave a
 * conspicuous on-screen gap. */
#define MAX_CHORD_FRAC 0.08

static bool sample_one(PlotSampleFn fn, void* ctx, double x, double* y) {
    return fn(x, ctx, y) && isfinite(*y);
}

/* Clamp a sampled y into the displayed band [lo, hi] before the flatness test.
 * A degenerate band (lo >= hi) means "no clip" and returns y unchanged. */
static double clamp_band(double y, double lo, double hi) {
    if (!(lo < hi)) return y;
    if (y < lo) return lo;
    if (y > hi) return hi;
    return y;
}

/* Recursively refine the interval [a, b] (endpoints already classified as
 * valid/invalid by the caller) and append every point on the *right* side
 * of each accepted leaf segment to `buf` (the left endpoint of the very
 * first interval is pushed once by the top-level driver below).
 *
 * Flatness is judged by the *vertical* gap between the curve and the chord at
 * three interior probes (the 1/4, 1/2 and 3/4 points), normalized by the
 * displayed y-extent `yspan`. Probing more than the midpoint alone is what
 * makes the sampler robust to oscillation: a single-midpoint test resonates
 * with periodic functions (e.g. Sin[22 x] on a coarse grid) -- three samples
 * land collinear while the curve wiggles between them, so the interval looks
 * flat and refinement never starts, no matter how high MaxRecursion is. With
 * three probes that false "flat" requires four aligned samples, which the
 * off-period grid does not provide, so the wiggle is caught and refined.
 * Vertical (not perpendicular) deviation is exactly the on-screen error for a
 * y = f(x) curve, and normalizing by the spike-clamped yspan keeps the test in
 * the coordinates the curve is actually displayed in.
 *
 * `yclip_lo`/`yclip_hi` are the displayed y-band (from an explicit PlotRange,
 * or the current zoom's visible band). Every probe is clamped into it before
 * the deviation is measured, so the test judges only the part of the curve the
 * window actually shows. The payoff is for functions that dive far outside the
 * frame -- a truncated Taylor series, a steep asymptote -- where the unclamped
 * test would pour refinement into the off-screen plunge (both chord and probes
 * agree it is "curved") and starve the visible crossing. Clamped, the off-band
 * stretch collapses onto the clip line and reads as flat, redirecting the
 * recursion budget to the on-screen body. A degenerate band (lo >= hi) is the
 * no-clip default and leaves the original behaviour untouched. */
static void subdivide(PlotSampleFn fn, void* ctx, SampleBuf* buf,
                       double ax, double ay, bool a_valid,
                       double bx, double by, bool b_valid,
                       int depth, int max_recursion, double yspan, double xspan,
                       double yclip_lo, double yclip_hi) {
    if (buf_over_budget(buf) || depth >= max_recursion) {
        if (b_valid) buf_push(buf, bx, by);
        else buf->pending_break = true;
        return;
    }

    double mx = (ax + bx) / 2.0, my = 0.0;
    bool m_valid = sample_one(fn, ctx, mx, &my);

    if (a_valid && b_valid && m_valid) {
        /* Probe the two quarter points as well; reject (refine) if any of the
         * three interior samples strays from the chord, or if a probe falls in
         * a gap (a singularity hides between the grid points). Deviation is
         * measured in the clamped (displayed) frame; see the header comment. */
        double q1x = (ax + mx) / 2.0, q3x = (mx + bx) / 2.0, q1y = 0.0, q3y = 0.0;
        bool q1_valid = sample_one(fn, ctx, q1x, &q1y);
        bool q3_valid = sample_one(fn, ctx, q3x, &q3y);
        if (q1_valid && q3_valid) {
            double cay = clamp_band(ay, yclip_lo, yclip_hi);
            double cby = clamp_band(by, yclip_lo, yclip_hi);
            double d1 = fabs(clamp_band(q1y, yclip_lo, yclip_hi) - (cay + 0.25 * (cby - cay)));
            double d2 = fabs(clamp_band(my,  yclip_lo, yclip_hi) - (cay + 0.50 * (cby - cay)));
            double d3 = fabs(clamp_band(q3y, yclip_lo, yclip_hi) - (cay + 0.75 * (cby - cay)));
            double dev = d1 > d2 ? d1 : d2;
            if (d3 > dev) dev = d3;
            /* Screen-space length of the chord the renderer would draw, with
             * each axis normalized by its displayed extent and y measured in
             * the clamped frame (so an off-band plunge collapses to zero length
             * and is not chased). See MAX_CHORD_FRAC. */
            double cdx = (bx - ax) / xspan;
            double cdy = (cby - cay) / yspan;
            double chord = sqrt(cdx * cdx + cdy * cdy);
            if (dev < FLAT_TOL * yspan && chord < MAX_CHORD_FRAC) {
                buf_push(buf, bx, by);
                return;
            }
        }
        /* else: a probe is invalid -- fall through and bisect to localize the
         * singularity, exactly as the invalid-endpoint case below does. */
    }

    subdivide(fn, ctx, buf, ax, ay, a_valid, mx, my, m_valid, depth + 1, max_recursion, yspan, xspan, yclip_lo, yclip_hi);
    subdivide(fn, ctx, buf, mx, my, m_valid, bx, by, b_valid, depth + 1, max_recursion, yspan, xspan, yclip_lo, yclip_hi);
}

PlotPoint* plot_sample_adaptive(PlotSampleFn fn, void* ctx,
                                 double xmin, double xmax,
                                 long plot_points, int max_recursion,
                                 long max_plot_points,
                                 double yclip_lo, double yclip_hi,
                                 size_t* out_count) {
    if (out_count) *out_count = 0;
    if (!fn || xmin >= xmax || plot_points < 2) return NULL;
    if (max_recursion < 0) max_recursion = 0;

    /* Pass 1: the uniform grid. We sample it up front (rather than streaming
     * one point per iteration) so we can estimate the vertical extent the
     * curve will be *displayed* over before refining. The curvature test in
     * subdivide() judges flatness in that normalized frame; see
     * segment_is_flat(). */
    double* gx = malloc(sizeof(double) * (size_t)plot_points);
    double* gy = malloc(sizeof(double) * (size_t)plot_points);
    bool*   gv = malloc(sizeof(bool)   * (size_t)plot_points);
    double* yvals = malloc(sizeof(double) * (size_t)plot_points);
    if (!gx || !gy || !gv || !yvals) {
        free(gx); free(gy); free(gv); free(yvals);
        return NULL;
    }

    size_t nvalid = 0;
    for (long i = 0; i < plot_points; i++) {
        double t = (double)i / (double)(plot_points - 1);
        double x = xmin + t * (xmax - xmin);
        if (i == plot_points - 1) x = xmax; /* avoid float drift off the edge */
        double y = 0.0;
        bool valid = sample_one(fn, ctx, x, &y);
        gx[i] = x; gy[i] = y; gv[i] = valid;
        /* Clamp into the displayed band so an off-screen excursion can't
         * inflate the robust y-extent and loosen the flatness tolerance. */
        if (valid) yvals[nvalid++] = clamp_band(y, yclip_lo, yclip_hi);
    }

    /* The flatness test normalizes the curve-to-chord gap by the *robust*
     * y-extent -- the same spike-clamped band render uses for the visible axis
     * -- so refinement tracks the displayed curve, and a lone asymptote spike
     * cannot inflate the scale and re-flatten everything else. */
    double yspan = xmax - xmin; /* square-ish fallback when nothing to estimate */
    if (nvalid > 0) {
        double ylo, yhi;
        plot_robust_yrange(yvals, nvalid, &ylo, &yhi);
        if (yhi > ylo) yspan = yhi - ylo;
    }
    free(yvals);
    if (!(yspan > 0.0)) yspan = 1.0;

    /* Y scale that normalizes both the flatness gap and the chord length. When
     * an explicit PlotRange band is in force this is the *displayed* height --
     * the curve is drawn over exactly that span, so a fraction of it is the
     * honest on-screen fraction. (The robust `yspan` above is the right scale
     * only for the auto-range case, where it both feeds the visible axis and
     * keeps a lone asymptote spike from loosening the test.) Using the robust
     * extent when the band is wider over-weights Δy and over-subdivides steep
     * stretches; using the band keeps screen-space spacing uniform -- and is the
     * frame ParametricPlot will need, where the curve is not a graph of x. */
    double ynorm = (yclip_hi > yclip_lo) ? (yclip_hi - yclip_lo) : yspan;
    if (!(ynorm > 0.0)) ynorm = 1.0;

    SampleBuf buf;
    buf_init(&buf, max_plot_points);

    if (gv[0]) buf_push(&buf, gx[0], gy[0]);
    else buf.pending_break = true;

    for (long i = 1; i < plot_points; i++) {
        subdivide(fn, ctx, &buf, gx[i - 1], gy[i - 1], gv[i - 1],
                  gx[i], gy[i], gv[i], 0, max_recursion, ynorm, (xmax - xmin), yclip_lo, yclip_hi);
    }

    free(gx); free(gy); free(gv);

    if (out_count) *out_count = buf.len;
    if (buf.len == 0) {
        free(buf.pts);
        return NULL;
    }
    return buf.pts;
}

void plot_points_free(PlotPoint* pts) {
    free(pts);
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Median of an already-sorted ascending array (n >= 1). */
static double sorted_median(const double* s, size_t n) {
    return (n & 1) ? s[n / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2]);
}

void plot_robust_yrange(const double* ys, size_t n, double* out_lo, double* out_hi) {
    if (n == 0) { *out_lo = 1.0; *out_hi = -1.0; return; }

    double* s = malloc(sizeof(double) * n);
    if (!s) { /* degrade to a linear min/max scan */
        double lo = ys[0], hi = ys[0];
        for (size_t i = 1; i < n; i++) { if (ys[i] < lo) lo = ys[i]; if (ys[i] > hi) hi = ys[i]; }
        *out_lo = lo; *out_hi = hi; return;
    }
    for (size_t i = 0; i < n; i++) s[i] = ys[i];
    qsort(s, n, sizeof(double), cmp_double);

    double data_min = s[0], data_max = s[n - 1];

    /* Too few points to estimate a robust spread -- show everything. */
    if (n < 8) { *out_lo = data_min; *out_hi = data_max; free(s); return; }

    double m = sorted_median(s, n);

    /* MAD = median(|y - m|). Reuse `s` for the deviations, then re-sort. */
    for (size_t i = 0; i < n; i++) s[i] = fabs(s[i] - m);
    qsort(s, n, sizeof(double), cmp_double);
    double mad = sorted_median(s, n);
    free(s);

    /* Flat bulk: no robust scale to clip against -- show everything. */
    if (mad <= 0.0) { *out_lo = data_min; *out_hi = data_max; return; }

    /* 1.4826 makes the scale match a standard deviation for normal data; K is
     * generous so smooth curves (Sin, polynomials, Exp growth) stay inside the
     * fence and only sparse asymptote spikes get clamped. */
    const double K = 4.0;
    double scale = 1.4826 * mad;
    double fence_lo = m - K * scale;
    double fence_hi = m + K * scale;

    /* Intersect with the true extent: this IS the spike guard. A genuine
     * extreme that lies within the fence survives (min/max kept); one that
     * shoots far past is clamped to the fence. */
    *out_lo = (data_min > fence_lo) ? data_min : fence_lo;
    *out_hi = (data_max < fence_hi) ? data_max : fence_hi;
}
