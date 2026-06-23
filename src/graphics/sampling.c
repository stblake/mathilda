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
 * the segment is accepted as visually flat. ~1/400 of the y-range is well
 * under a pixel on a normal window, so accepted segments read as smooth. */
#define FLAT_TOL 0.0025

static bool sample_one(PlotSampleFn fn, void* ctx, double x, double* y) {
    return fn(x, ctx, y) && isfinite(*y);
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
 * the coordinates the curve is actually displayed in. */
static void subdivide(PlotSampleFn fn, void* ctx, SampleBuf* buf,
                       double ax, double ay, bool a_valid,
                       double bx, double by, bool b_valid,
                       int depth, int max_recursion, double yspan) {
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
         * a gap (a singularity hides between the grid points). */
        double q1x = (ax + mx) / 2.0, q3x = (mx + bx) / 2.0, q1y = 0.0, q3y = 0.0;
        bool q1_valid = sample_one(fn, ctx, q1x, &q1y);
        bool q3_valid = sample_one(fn, ctx, q3x, &q3y);
        if (q1_valid && q3_valid) {
            double d1 = fabs(q1y - (ay + 0.25 * (by - ay)));
            double d2 = fabs(my  - (ay + 0.50 * (by - ay)));
            double d3 = fabs(q3y - (ay + 0.75 * (by - ay)));
            double dev = d1 > d2 ? d1 : d2;
            if (d3 > dev) dev = d3;
            if (dev < FLAT_TOL * yspan) {
                buf_push(buf, bx, by);
                return;
            }
        }
        /* else: a probe is invalid -- fall through and bisect to localize the
         * singularity, exactly as the invalid-endpoint case below does. */
    }

    subdivide(fn, ctx, buf, ax, ay, a_valid, mx, my, m_valid, depth + 1, max_recursion, yspan);
    subdivide(fn, ctx, buf, mx, my, m_valid, bx, by, b_valid, depth + 1, max_recursion, yspan);
}

PlotPoint* plot_sample_adaptive(PlotSampleFn fn, void* ctx,
                                 double xmin, double xmax,
                                 long plot_points, int max_recursion,
                                 long max_plot_points, size_t* out_count) {
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
        if (valid) yvals[nvalid++] = y;
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

    SampleBuf buf;
    buf_init(&buf, max_plot_points);

    if (gv[0]) buf_push(&buf, gx[0], gy[0]);
    else buf.pending_break = true;

    for (long i = 1; i < plot_points; i++) {
        subdivide(fn, ctx, &buf, gx[i - 1], gy[i - 1], gv[i - 1],
                  gx[i], gy[i], gv[i], 0, max_recursion, yspan);
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
