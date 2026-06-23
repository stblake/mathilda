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

/* Curvature test on the (a, mid, b) triple: angle between the a->mid and
 * mid->b segment vectors. 0 = perfectly straight. Degenerate (near-zero
 * length) segments are treated as flat so the recursion terminates. */
static bool segment_is_flat(double ax, double ay, double mx, double my,
                             double bx, double by) {
    double v1x = mx - ax, v1y = my - ay;
    double v2x = bx - mx, v2y = by - my;
    double len1 = hypot(v1x, v1y), len2 = hypot(v2x, v2y);
    if (len1 < 1e-12 || len2 < 1e-12) return true;
    double cosang = (v1x * v2x + v1y * v2y) / (len1 * len2);
    if (cosang > 1.0) cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;
    const double ANGLE_THRESHOLD = 0.06; /* radians, ~3.4 degrees */
    return acos(cosang) < ANGLE_THRESHOLD;
}

static bool sample_one(PlotSampleFn fn, void* ctx, double x, double* y) {
    return fn(x, ctx, y) && isfinite(*y);
}

/* Recursively refine the interval [a, b] (endpoints already classified as
 * valid/invalid by the caller) and append every point on the *right* side
 * of each accepted leaf segment to `buf` (the left endpoint of the very
 * first interval is pushed once by the top-level driver below). */
static void subdivide(PlotSampleFn fn, void* ctx, SampleBuf* buf,
                       double ax, double ay, bool a_valid,
                       double bx, double by, bool b_valid,
                       int depth, int max_recursion) {
    if (buf_over_budget(buf) || depth >= max_recursion) {
        if (b_valid) buf_push(buf, bx, by);
        else buf->pending_break = true;
        return;
    }

    double mx = (ax + bx) / 2.0, my = 0.0;
    bool m_valid = sample_one(fn, ctx, mx, &my);

    if (!a_valid || !b_valid || !m_valid) {
        /* A singularity is somewhere in this interval; keep bisecting to
         * localize it instead of running the (cheaper) curvature test,
         * which needs all three points to be meaningful. */
        subdivide(fn, ctx, buf, ax, ay, a_valid, mx, my, m_valid, depth + 1, max_recursion);
        subdivide(fn, ctx, buf, mx, my, m_valid, bx, by, b_valid, depth + 1, max_recursion);
        return;
    }

    if (segment_is_flat(ax, ay, mx, my, bx, by)) {
        buf_push(buf, bx, by);
        return;
    }

    subdivide(fn, ctx, buf, ax, ay, true, mx, my, true, depth + 1, max_recursion);
    subdivide(fn, ctx, buf, mx, my, true, bx, by, true, depth + 1, max_recursion);
}

PlotPoint* plot_sample_adaptive(PlotSampleFn fn, void* ctx,
                                 double xmin, double xmax,
                                 long plot_points, int max_recursion,
                                 long max_plot_points, size_t* out_count) {
    if (out_count) *out_count = 0;
    if (!fn || xmin >= xmax || plot_points < 2) return NULL;
    if (max_recursion < 0) max_recursion = 0;

    SampleBuf buf;
    buf_init(&buf, max_plot_points);

    double prev_x = xmin, prev_y = 0.0;
    bool prev_valid = sample_one(fn, ctx, prev_x, &prev_y);
    if (prev_valid) buf_push(&buf, prev_x, prev_y);
    else buf.pending_break = true;

    for (long i = 1; i < plot_points; i++) {
        double t = (double)i / (double)(plot_points - 1);
        double x = xmin + t * (xmax - xmin);
        if (i == plot_points - 1) x = xmax; /* avoid float drift off the edge */
        double y = 0.0;
        bool valid = sample_one(fn, ctx, x, &y);

        subdivide(fn, ctx, &buf, prev_x, prev_y, prev_valid, x, y, valid, 0, max_recursion);

        prev_x = x;
        prev_y = y;
        prev_valid = valid;
    }

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
