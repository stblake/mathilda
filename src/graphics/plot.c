/* plot.c — Plot[f, {x, xmin, xmax}, opts...].
 *
 * HoldAll, like Table/Do: f and the iterator spec must not be
 * pre-evaluated (x has no value yet). Splits trailing options into the
 * sampler's own (PlotPoints/MaxRecursion/MaxPlotPoints/Mesh/RegionFunction/
 * Exclusions/ColorFunction/Filling/..., consumed here) and everything else
 * (PlotRange/AspectRatio/PlotStyle/Axes/.../ImageSize, copied through onto
 * the resulting Graphics[...] unevaluated -- render.c is the single place
 * that interprets those, whether reached via Plot's auto-display or a
 * later Show[]). */

#include "plot.h"
#include "plot_common.h"
#include "show.h"
#include "sampling.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "print.h"
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* True if `e` is a 2-element List of numericizable values, returning them. */
static bool list2_nums(Expr* e, double* a, double* b) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List
        && e->data.function.arg_count == 2
        && numericize_bound(e->data.function.args[0], a)
        && numericize_bound(e->data.function.args[1], b);
}

/* Extract an explicit numeric y-band from a PlotRange value, mirroring the
 * forms render.c honours: {{xmin,xmax},{ymin,ymax}} and {xspec,{ymin,ymax}}
 * fix y via the second sublist, while a bare {ymin,ymax} pins y directly.
 * Anything with a non-numeric (Automatic/All) y side leaves the outputs
 * untouched and returns false. The band is fed to the sampler as a clip frame;
 * see PlotSampleOpts. */
static bool plotrange_yband(Expr* rhs, double* lo, double* hi) {
    Expr* v = evaluate(expr_copy(rhs));
    double a, b;
    bool ok = false;
    if (v->type == EXPR_FUNCTION && v->data.function.head
        && v->data.function.head->type == EXPR_SYMBOL
        && v->data.function.head->data.symbol == SYM_List
        && v->data.function.arg_count == 2) {
        Expr* a1 = v->data.function.args[1];
        if (list2_nums(a1, &a, &b)) ok = true;                 /* {xspec, {ymin,ymax}} */
        else ok = numericize_bound(v->data.function.args[0], &a)
               && numericize_bound(a1, &b);                    /* {ymin, ymax} */
    }
    expr_free(v);
    if (ok && a > b) { double t = a; a = b; b = t; }
    if (ok && a < b) { *lo = a; *hi = b; return true; }
    return false;
}

typedef struct {
    Expr* var;             /* iterator symbol, borrowed */
    Expr* body;             /* f, borrowed */
    Expr* region_function;  /* borrowed; NULL = none */
} PlotEvalCtx;

static bool plot_eval_fn(double x, void* ctx_, double* y_out) {
    PlotEvalCtx* ctx = (PlotEvalCtx*)ctx_;
    Expr* xval = expr_new_real(x);
    symtab_add_own_value(ctx->var->data.symbol, ctx->var, xval);
    Expr* result = evaluate(ctx->body);

    double y;
    bool ok = expr_to_real_double(result, &y) && isfinite(y);
    expr_free(result);

    if (ok && ctx->region_function && !eval_region(ctx->region_function, x, y)) ok = false;

    expr_free(xval);
    if (ok) *y_out = y;
    return ok;
}

/* ---- Option parsing ---- */

typedef struct {
    long plot_points;
    int  max_recursion;
    long max_plot_points; /* <= 0 means unbounded */
    bool mesh;            /* Mesh -> All: also emit the sample points as dots */
    Expr* region_function; /* borrowed; held; NULL = none */
    Expr* exclusions;      /* borrowed; held; NULL = none */
    Expr* color_function;  /* borrowed; held (function, or string "Rainbow"); NULL = none */
    bool  color_function_scaling; /* default true */
    Expr* filling;         /* borrowed; held (Axis/Bottom/Top symbol, or a number); NULL = none */
    Expr* filling_style;   /* borrowed; held color; NULL = default (curve colour @ Opacity[0.3]) */
    /* Displayed y-band from an explicit numeric PlotRange, fed to the sampler
     * so a curve that dives off-screen (a divergent Taylor tail, a steep
     * asymptote) doesn't starve refinement of its on-screen body. A degenerate
     * band (lo >= hi) means "no explicit PlotRange y" -> sample full extent. */
    double yclip_lo, yclip_hi;
} PlotSampleOpts;

/* Splits res's trailing Rule args (starting at index 2) into the sampler
 * options above and a passthrough list of borrowed-then-copied Rule
 * expressions destined for the Graphics[...] result. `single_color_out`
 * receives a freshly owned copy of the resolved PlotStyle color (or the
 * default blue, if PlotStyle wasn't given) for single-curve Filling/
 * restore-after-fill use; caller must expr_free it. Returns false on a
 * malformed trailing argument (not a Rule, or a known option with a
 * badly-typed value) so the caller can decline evaluation entirely. */
static bool split_options(Expr* res, PlotSampleOpts* sopts,
                           Expr*** passthrough_out, size_t* passthrough_count_out,
                           Expr** single_color_out) {
    sopts->plot_points = 50;
    sopts->max_recursion = 6;
    sopts->max_plot_points = -1;
    sopts->mesh = false;
    sopts->region_function = NULL;
    sopts->exclusions = NULL;
    sopts->color_function = NULL;
    sopts->color_function_scaling = true;
    sopts->filling = NULL;
    sopts->filling_style = NULL;
    sopts->yclip_lo = 0.0;
    sopts->yclip_hi = -1.0; /* degenerate => no clip until an explicit PlotRange y is seen */
    *single_color_out = NULL;

    size_t argc = res->data.function.arg_count;
    /* +3 headroom for the Axes/AspectRatio/PlotStyle defaults potentially
     * appended below when the caller didn't already supply them. */
    size_t cap = (argc > 2 ? argc - 2 : 0) + 3;
    Expr** passthrough = malloc(sizeof(Expr*) * cap);
    size_t n = 0;

    bool have_axes = false, have_aspect = false, have_style = false, have_frame = false;

#define FAIL_CLEANUP() do { free(passthrough); expr_free(*single_color_out); *single_color_out = NULL; return false; } while (0)

    for (size_t i = 2; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) FAIL_CLEANUP();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol : NULL;

        if (name == SYM_PlotPoints) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 2) FAIL_CLEANUP();
            sopts->plot_points = v;
        } else if (name == SYM_MaxRecursion) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 0) FAIL_CLEANUP();
            sopts->max_recursion = (int)v;
        } else if (name == SYM_MaxPlotPoints) {
            Expr* v = evaluate(expr_copy(rhs));
            bool is_inf = (v->type == EXPR_SYMBOL && v->data.symbol == SYM_Infinity);
            long lv = -1;
            bool ok = is_inf;
            if (!ok && v->type == EXPR_INTEGER && v->data.integer > 0) { lv = (long)v->data.integer; ok = true; }
            expr_free(v);
            if (!ok) FAIL_CLEANUP();
            sopts->max_plot_points = lv;
        } else if (name == SYM_Mesh) {
            /* Mesh -> All (or True) overlays the evaluation points as dots;
             * None/False (the default) draws the line only. Consumed here, so
             * it never reaches the Graphics[...] result. */
            Expr* v = evaluate(expr_copy(rhs));
            bool on = (v->type == EXPR_SYMBOL
                       && (v->data.symbol == SYM_All || v->data.symbol == SYM_True));
            bool off = (v->type == EXPR_SYMBOL
                        && (v->data.symbol == SYM_None || v->data.symbol == SYM_False));
            expr_free(v);
            if (!on && !off) FAIL_CLEANUP();
            sopts->mesh = on;
        } else if (name == SYM_RegionFunction) {
            sopts->region_function = rhs;
        } else if (name == SYM_Exclusions) {
            if (!(rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_None)) sopts->exclusions = rhs;
        } else if (name == SYM_ColorFunction) {
            sopts->color_function = rhs;
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = evaluate(expr_copy(rhs));
            sopts->color_function_scaling = !(v->type == EXPR_SYMBOL && v->data.symbol == SYM_False);
            expr_free(v);
        } else if (name == SYM_Filling) {
            if (!(rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_None)) sopts->filling = rhs;
        } else if (name == SYM_FillingStyle) {
            sopts->filling_style = rhs;
        } else if (name == SYM_PlotStyle) {
            have_style = true;
            if (*single_color_out) expr_free(*single_color_out);
            *single_color_out = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), expr_copy(*single_color_out) };
            passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        } else if (name == SYM_AspectRatio) {
            /* The renderer has no evaluator, so resolve the value here. The
             * symbolic settings Automatic and Full are interpreted downstream
             * and pass through verbatim; anything else (a plain ratio, or a
             * symbolic-numeric one like 1/GoldenRatio) is numericalized to a
             * real so render.c receives a machine number. A value that won't
             * reduce to a positive real falls back to Automatic. */
            have_aspect = true;
            if (rhs->type == EXPR_SYMBOL
                && (rhs->data.symbol == SYM_Automatic || rhs->data.symbol == SYM_Full)) {
                passthrough[n++] = expr_copy(arg);
            } else {
                double v;
                Expr* val = (numericize_bound(rhs, &v) && v > 0)
                    ? expr_new_real(v) : expr_new_symbol(SYM_Automatic);
                Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), val };
                passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
            }
        } else {
            if (name == SYM_Axes) have_axes = true;
            else if (name == SYM_Frame) {
                /* A frame (Frame -> True or a per-edge spec) takes the place of
                 * the interior Axes cross, so suppress Plot's Axes -> True
                 * default below. Frame -> False/None leaves the axes default. */
                if (!(rhs->type == EXPR_SYMBOL
                      && (rhs->data.symbol == SYM_False || rhs->data.symbol == SYM_None)))
                    have_frame = true;
            }
            else if (name == SYM_PlotRange) {
                /* An explicit numeric y-range becomes the sampler's clip frame
                 * so refinement concentrates on the on-screen curve rather than
                 * an off-screen plunge. The value still passes through below to
                 * the Graphics[...] result, which clips the display to it. */
                plotrange_yband(rhs, &sopts->yclip_lo, &sopts->yclip_hi);
            }
            /* Plot is HoldAll, so this option's value would otherwise reach
             * the Graphics[...] result completely unevaluated -- e.g. a
             * named color constant like Red (an OwnValue -> RGBColor[...])
             * would never resolve. Evaluate it here, once, exactly as a
             * non-held Graphics[]'s own arguments already would. */
            Expr* val = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    /* Plot-specific defaults (distinct from bare Graphics[]'s defaults of
     * Axes->False / AspectRatio->Automatic): inject only what the caller
     * didn't already specify. */
    if (!have_axes && !have_frame) {
        Expr* a[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_True) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_aspect) {
        const double inv_phi = 2.0 / (1.0 + sqrt(5.0)); /* 1/GoldenRatio */
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_real(inv_phi) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_style) {
        Expr* rgb_args[3] = { expr_new_real(0.2), expr_new_real(0.4), expr_new_real(0.8) };
        Expr* rgb = expr_new_function(expr_new_symbol(SYM_RGBColor), rgb_args, 3);
        Expr* a[2] = { expr_new_symbol(SYM_PlotStyle), expr_copy(rgb) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        *single_color_out = rgb;
    }

    *passthrough_out = passthrough;
    *passthrough_count_out = n;
    return true;
}
#undef FAIL_CLEANUP

/* ---- Exclusions: explicit forced discontinuities ---- */

typedef struct { double lo, hi; } Range1D;

/* Numericizes each element of `exclusions` (a number, an equation x == a,
 * or a List of either) into a sorted array of x-values. Unparseable
 * elements are silently skipped. */
static double* parse_exclusions(Expr* exclusions, Expr* var, size_t* out_count) {
    *out_count = 0;
    if (!exclusions) return NULL;

    Expr** items;
    size_t n;
    bool is_list = (exclusions->type == EXPR_FUNCTION && exclusions->data.function.head->type == EXPR_SYMBOL
                    && exclusions->data.function.head->data.symbol == SYM_List);
    if (is_list) { items = exclusions->data.function.args; n = exclusions->data.function.arg_count; }
    else { items = &exclusions; n = 1; }
    if (n == 0) return NULL;

    double* out = malloc(sizeof(double) * n);
    size_t cnt = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* it = items[i];
        Expr* value_expr = it;
        /* x == a (Equal[x, a]): take the side that isn't the iterator var. */
        if (it->type == EXPR_FUNCTION && it->data.function.head->type == EXPR_SYMBOL
            && it->data.function.head->data.symbol == SYM_Equal && it->data.function.arg_count == 2) {
            Expr* lhs = it->data.function.args[0];
            Expr* rhsv = it->data.function.args[1];
            value_expr = (lhs->type == EXPR_SYMBOL && lhs->data.symbol == var->data.symbol) ? rhsv : lhs;
        }
        double v;
        if (numericize_bound(value_expr, &v)) out[cnt++] = v;
    }
    /* Small insertion sort -- exclusion lists are always tiny. */
    for (size_t i = 1; i < cnt; i++) {
        double key = out[i];
        size_t j = i;
        while (j > 0 && out[j - 1] > key) { out[j] = out[j - 1]; j--; }
        out[j] = key;
    }
    *out_count = cnt;
    return out;
}

/* Splits [xmin,xmax] into sub-ranges at each in-range exclusion x, leaving a
 * tiny epsilon gap at each split so the curve doesn't visually rejoin. */
static Range1D* split_at_exclusions(double xmin, double xmax, const double* excl, size_t nexcl,
                                     size_t* out_count) {
    double eps = (xmax - xmin) * 1e-6;
    Range1D* ranges = malloc(sizeof(Range1D) * (nexcl + 1));
    size_t rc = 0;
    double lo = xmin;
    for (size_t i = 0; i < nexcl; i++) {
        if (excl[i] > lo + eps && excl[i] < xmax - eps) {
            ranges[rc++] = (Range1D){ lo, excl[i] - eps };
            lo = excl[i] + eps;
        }
    }
    ranges[rc++] = (Range1D){ lo, xmax };
    *out_count = rc;
    return ranges;
}

/* ---- Filling: baseline resolution + per-run fill polygon ---- */

static double filling_baseline(Expr* filling, double run_ymin, double run_ymax) {
    if (!filling) return 0.0;
    if (filling->type == EXPR_SYMBOL) {
        if (filling->data.symbol == SYM_Axis) return 0.0;
        if (filling->data.symbol == SYM_Bottom) return run_ymin;
        if (filling->data.symbol == SYM_Top) return run_ymax;
    }
    double v;
    if (numericize_bound(filling, &v)) return v;
    return 0.0;
}

/* One small quad Polygon[] per consecutive point pair in the run, each
 * (p_i, p_{i+1}, (x_{i+1},baseline), (x_i,baseline)) -- rather than one
 * big polygon tracing the whole run then back along the baseline.
 *
 * render.c's Polygon fill is a triangle *fan* from the first vertex, which
 * only renders correctly for a star-shaped-from-vertex-0 outline (convex
 * is the easy common case). A fill spanning an entire wavy run (e.g. a
 * full period of Sin[x]) is nowhere close to that and fans out into a
 * wedge of garbage triangles instead of hugging the curve. Each individual
 * quad between two adjacent (and thus close-together) sample points is
 * always simple regardless of the curve's overall shape, so this sidesteps
 * the fan's star-shaped requirement entirely. Returns the quad count via
 * *out_count (0 if run_len < 2); caller owns the returned array and Exprs. */
/* Appends a 4-vertex Polygon[] (x0,y0)-(x1,y1)-(x1,baseline)-(x0,baseline)
 * to `shapes` at *n, advancing *n by 1. */
static void push_fill_quad(Expr** shapes, size_t* n, double x0, double y0,
                            double x1, double y1, double baseline) {
    Expr* c0[2] = { expr_new_real(x0), expr_new_real(y0) };
    Expr* c1[2] = { expr_new_real(x1), expr_new_real(y1) };
    Expr* c2[2] = { expr_new_real(x1), expr_new_real(baseline) };
    Expr* c3[2] = { expr_new_real(x0), expr_new_real(baseline) };
    Expr* verts[4] = {
        expr_new_function(expr_new_symbol(SYM_List), c0, 2),
        expr_new_function(expr_new_symbol(SYM_List), c1, 2),
        expr_new_function(expr_new_symbol(SYM_List), c2, 2),
        expr_new_function(expr_new_symbol(SYM_List), c3, 2),
    };
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, 4);
    Expr* poly_args[1] = { vlist };
    shapes[(*n)++] = expr_new_function(expr_new_symbol(SYM_Polygon), poly_args, 1);
}

/* Appends a 3-vertex Polygon[] (triangle) to `shapes` at *n. */
static void push_fill_triangle(Expr** shapes, size_t* n, double x0, double y0,
                                double x1, double y1, double x2, double y2) {
    Expr* a0[2] = { expr_new_real(x0), expr_new_real(y0) };
    Expr* a1[2] = { expr_new_real(x1), expr_new_real(y1) };
    Expr* a2[2] = { expr_new_real(x2), expr_new_real(y2) };
    Expr* verts[3] = {
        expr_new_function(expr_new_symbol(SYM_List), a0, 2),
        expr_new_function(expr_new_symbol(SYM_List), a1, 2),
        expr_new_function(expr_new_symbol(SYM_List), a2, 2),
    };
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, 3);
    Expr* poly_args[1] = { vlist };
    shapes[(*n)++] = expr_new_function(expr_new_symbol(SYM_Polygon), poly_args, 1);
}

static Expr** build_fill_quads(const PlotPoint* pts, size_t run_start, size_t run_len,
                                double baseline, size_t* out_count) {
    if (run_len < 2) { *out_count = 0; return NULL; }
    size_t nseg = run_len - 1;
    /* Worst case: every segment crosses the baseline and needs 2 triangles
     * instead of 1 quad. */
    Expr** shapes = malloc(sizeof(Expr*) * nseg * 2);
    size_t n = 0;
    for (size_t j = 0; j < nseg; j++) {
        double x0 = pts[run_start + j].x,     y0 = pts[run_start + j].y;
        double x1 = pts[run_start + j + 1].x, y1 = pts[run_start + j + 1].y;
        double d0 = y0 - baseline, d1 = y1 - baseline;

        if ((d0 > 0 && d1 < 0) || (d0 < 0 && d1 > 0)) {
            /* The segment crosses the baseline -- one quad here would be a
             * self-intersecting "bowtie" (its closing edge, lying exactly
             * along the baseline, crosses the curve edge at the very point
             * the curve crosses zero), which a triangle fan turns into a
             * stray sliver right at the crossing. Split into two simple,
             * always-fan-correct triangles instead, one on each side. */
            double t = d0 / (d0 - d1);
            double xc = x0 + t * (x1 - x0);
            push_fill_triangle(shapes, &n, x0, y0, xc, baseline, x0, baseline);
            push_fill_triangle(shapes, &n, xc, baseline, x1, y1, x1, baseline);
        } else {
            push_fill_quad(shapes, &n, x0, y0, x1, y1, baseline);
        }
    }
    *out_count = n;
    return shapes;
}

/* Sample one function `body` over [xmin,xmax] (split at any in-range
 * Exclusions, each sub-range subject to RegionFunction) and return its
 * primitives: an optional Filling Polygon[] per run, then either one plain
 * Line[run] (no ColorFunction) or one color-directive + 2-point Line[] pair
 * per segment (ColorFunction set), plus the optional Mesh dot overlay.
 * `curve_color` (borrowed, may be NULL) is the color this curve is drawn in
 * -- used only to restore it after a Filling polygon when ColorFunction
 * isn't also overriding every segment's color anyway. Returns NULL with
 * *out_count == 0 when nothing is plottable. The caller must already have
 * shadowed the iterator variable. */
static Expr** sample_lines(Expr* body, Expr* var, double xmin, double xmax,
                            const PlotSampleOpts* sopts, Expr* curve_color,
                            size_t* out_count) {
    *out_count = 0;

    size_t nexcl = 0;
    double* excl = parse_exclusions(sopts->exclusions, var, &nexcl);
    size_t nranges = 0;
    Range1D* ranges = split_at_exclusions(xmin, xmax, excl, nexcl, &nranges);
    free(excl);

    PlotPoint* pts = NULL;
    size_t npts = 0, cap = 0;
    for (size_t r = 0; r < nranges; r++) {
        if (!(ranges[r].lo < ranges[r].hi)) continue;
        PlotEvalCtx ctx = { .var = var, .body = body, .region_function = sopts->region_function };
        size_t rn;
        PlotPoint* rpts = plot_sample_adaptive(plot_eval_fn, &ctx, ranges[r].lo, ranges[r].hi,
                                                sopts->plot_points, sopts->max_recursion,
                                                sopts->max_plot_points,
                                                sopts->yclip_lo, sopts->yclip_hi, &rn);
        if (!rpts || rn == 0) { plot_points_free(rpts); continue; }
        for (size_t i = 0; i < rn; i++) {
            if (npts == cap) { cap = cap ? cap * 2 : 64; pts = realloc(pts, sizeof(PlotPoint) * cap); }
            PlotPoint p = rpts[i];
            if (i == 0 && npts > 0) p.break_before = true; /* force a gap between sub-ranges */
            pts[npts++] = p;
        }
        plot_points_free(rpts);
    }
    free(ranges);

    if (!pts || npts == 0) { free(pts); return NULL; }

    /* Worst case: every point becomes its own colored 2-point segment (2
     * primitives each, ~2*npts) plus Filling's strip (up to 2*(npts-1)
     * shapes, since a baseline-crossing segment becomes 2 triangles, plus
     * 2 Opacity brackets and a color restore per run) plus the Mesh
     * overlay's PointSize/Point pair. Generous on purpose -- a malloc
     * over-allocation here is harmless, unlike under-allocating. */
    Expr** prims = malloc(sizeof(Expr*) * (npts * 6 + 24));
    size_t prim_count = 0;
    size_t run_start = 0;
    for (size_t i = 1; i <= npts; i++) {
        bool end_of_run = (i == npts) || pts[i].break_before;
        if (end_of_run) {
            size_t run_len = i - run_start;
            if (run_len >= 2) {
                if (sopts->filling) {
                    double run_ymin = pts[run_start].y, run_ymax = pts[run_start].y;
                    for (size_t j = 1; j < run_len; j++) {
                        double yy = pts[run_start + j].y;
                        if (yy < run_ymin) run_ymin = yy;
                        if (yy > run_ymax) run_ymax = yy;
                    }
                    double baseline = filling_baseline(sopts->filling, run_ymin, run_ymax);
                    if (sopts->filling_style) {
                        prims[prim_count++] = evaluate(expr_copy(sopts->filling_style));
                    } else {
                        Expr* op_arg[1] = { expr_new_real(0.3) };
                        prims[prim_count++] = expr_new_function(expr_new_symbol(SYM_Opacity), op_arg, 1);
                    }
                    size_t nquads;
                    Expr** quads = build_fill_quads(pts, run_start, run_len, baseline, &nquads);
                    for (size_t qi = 0; qi < nquads; qi++) prims[prim_count++] = quads[qi];
                    free(quads);
                    if (!sopts->filling_style) {
                        Expr* op_arg2[1] = { expr_new_integer(1) };
                        prims[prim_count++] = expr_new_function(expr_new_symbol(SYM_Opacity), op_arg2, 1);
                    }
                    if (curve_color && !sopts->color_function) prims[prim_count++] = expr_copy(curve_color);
                }

                if (sopts->color_function) {
                    for (size_t j = 0; j + 1 < run_len; j++) {
                        size_t a = run_start + j, b = a + 1;
                        double midx = (pts[a].x + pts[b].x) / 2.0;
                        double midy = (pts[a].y + pts[b].y) / 2.0;
                        prims[prim_count++] = eval_color_function(sopts->color_function, midx, midy,
                                                                   xmin, xmax, sopts->color_function_scaling);
                        Expr* xyA[2] = { expr_new_real(pts[a].x), expr_new_real(pts[a].y) };
                        Expr* xyB[2] = { expr_new_real(pts[b].x), expr_new_real(pts[b].y) };
                        Expr* seg_pts[2] = {
                            expr_new_function(expr_new_symbol(SYM_List), xyA, 2),
                            expr_new_function(expr_new_symbol(SYM_List), xyB, 2),
                        };
                        Expr* seg_list = expr_new_function(expr_new_symbol(SYM_List), seg_pts, 2);
                        Expr* line_args[1] = { seg_list };
                        prims[prim_count++] = expr_new_function(expr_new_symbol(SYM_Line), line_args, 1);
                    }
                } else {
                    Expr** line_pts = malloc(sizeof(Expr*) * run_len);
                    for (size_t j = 0; j < run_len; j++) {
                        Expr* xy[2] = { expr_new_real(pts[run_start + j].x), expr_new_real(pts[run_start + j].y) };
                        line_pts[j] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
                    }
                    Expr* pts_list2 = expr_new_function(expr_new_symbol(SYM_List), line_pts, run_len);
                    free(line_pts);
                    Expr* line_args[1] = { pts_list2 };
                    prims[prim_count++] = expr_new_function(expr_new_symbol(SYM_Line), line_args, 1);
                }
            }
            run_start = i;
        }
    }

    /* Mesh -> All: overlay one dot per evaluation point, in the curve's colour
     * (the directives precede these so the dots inherit it). Sized as a small
     * fraction of the x-span, which render scales by the same zoom as the
     * curve, so dots stay ~constant on screen regardless of the plot range. */
    if (sopts->mesh && npts > 0) {
        double msize = (xmax - xmin) * 0.0025;
        Expr* ps_arg[1] = { expr_new_real(msize > 0 ? msize : 0.005) };
        prims[prim_count++] = expr_new_function(expr_new_symbol(SYM_PointSize), ps_arg, 1);

        Expr** dots = malloc(sizeof(Expr*) * npts);
        for (size_t j = 0; j < npts; j++) {
            Expr* xy[2] = { expr_new_real(pts[j].x), expr_new_real(pts[j].y) };
            dots[j] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
        }
        Expr* dot_list = expr_new_function(expr_new_symbol(SYM_List), dots, npts);
        free(dots);
        Expr* pt_args[1] = { dot_list };
        prims[prim_count++] = expr_new_function(expr_new_symbol(SYM_Point), pt_args, 1);
    }

    plot_points_free(pts);
    if (prim_count == 0) { free(prims); return NULL; }
    *out_count = prim_count;
    return prims;
}

/* Sample every body over [xmin,xmax] and assemble the flat primitive
 * List[...] that becomes the first argument of Graphics[...]: palette colour
 * directives interleaved for multi-curve plots, plain Line runs for a single
 * curve (which inherits the PlotStyle colour at render time). Shadows `var`
 * internally. `single_color` (borrowed, may be NULL) is the resolved
 * PlotStyle color used for single-curve Filling/colour-restore. Returns NULL
 * when nothing is plottable. Shared by builtin_plot (initial draw) and
 * plot_resample (re-draw on zoom). */
static Expr* build_plot_primitives(Expr** bodies, size_t nfun, Expr* var,
                                    double xmin, double xmax,
                                    const PlotSampleOpts* sopts, Expr* single_color) {
    if (nfun == 0 || !(xmin < xmax)) return NULL;

    Rule* old_own = iter_spec_shadow(var);
    Expr*** per = malloc(sizeof(Expr**) * nfun);
    size_t* per_count = malloc(sizeof(size_t) * nfun);
    size_t total = 0;
    bool any = false;
    bool multi = (nfun > 1);
    for (size_t fi = 0; fi < nfun; fi++) {
        Expr* curve_color = multi ? palette_color(fi) : single_color;
        per[fi] = sample_lines(bodies[fi], var, xmin, xmax, sopts, curve_color, &per_count[fi]);
        if (multi) expr_free(curve_color);
        total += per_count[fi];
        if (per_count[fi] > 0) any = true;
    }
    iter_spec_restore(var, old_own);

    if (!any) { free(per); free(per_count); return NULL; }

    size_t cap = total + (multi ? nfun : 0);
    Expr** prims = malloc(sizeof(Expr*) * (cap > 0 ? cap : 1));
    size_t prim_count = 0;
    for (size_t fi = 0; fi < nfun; fi++) {
        if (multi) prims[prim_count++] = palette_color(fi);
        for (size_t j = 0; j < per_count[fi]; j++) prims[prim_count++] = per[fi][j];
        free(per[fi]);
    }
    free(per); free(per_count);

    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, prim_count);
    free(prims);
    return prim_list;
}

/* PlotLegends -> Automatic | "Expressions" | {labels...}: builds the
 * internal $PlotLegendData[{color1,label1}, ...] metadata, read by
 * render.c's draw_legend(). `legends` is the (already-evaluated) option
 * value; bodies/nfun supply per-curve colors and "Expressions" labels.
 * Returns NULL if PlotLegends wasn't given (or resolves to None). */
static Expr* build_legend_meta(Expr* legends, Expr** bodies, size_t nfun, Expr* single_color) {
    if (!legends) return NULL;
    if (legends->type == EXPR_SYMBOL && legends->data.symbol == SYM_None) return NULL;

    bool explicit_list = (legends->type == EXPR_FUNCTION && legends->data.function.head->type == EXPR_SYMBOL
                           && legends->data.function.head->data.symbol == SYM_List);
    bool multi = (nfun > 1);

    Expr** entries = malloc(sizeof(Expr*) * (nfun > 0 ? nfun : 1));
    for (size_t i = 0; i < nfun; i++) {
        Expr* color = multi ? palette_color(i) : (single_color ? expr_copy(single_color) : palette_color(0));
        Expr* label;
        if (explicit_list && i < legends->data.function.arg_count) {
            label = expr_copy(legends->data.function.args[i]);
        } else {
            char* s = expr_to_string(bodies[i]);
            label = expr_new_string(s ? s : "");
            free(s);
        }
        Expr* a[2] = { color, label };
        entries[i] = expr_new_function(expr_new_symbol(SYM_List), a, 2);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_PlotLegendData), entries, nfun);
    free(entries);
    return result;
}

/* Build the hidden $PlotResample[var, {bodies}, {opts...}] node embedded in
 * the returned Graphics[...], capturing everything plot_resample needs to
 * re-sample later -- including ColorFunction/Filling/RegionFunction/
 * Exclusions, so zooming doesn't silently drop them. $PlotResample is
 * HoldAll, so the bodies and held option values survive the Graphics
 * re-evaluation unevaluated. */
static Expr* build_resample_meta(Expr** bodies, size_t nfun, Expr* var,
                                  const PlotSampleOpts* sopts) {
    Expr** bcopies = malloc(sizeof(Expr*) * (nfun > 0 ? nfun : 1));
    for (size_t i = 0; i < nfun; i++) bcopies[i] = expr_copy(bodies[i]);
    Expr* blist = expr_new_function(expr_new_symbol(SYM_List), bcopies, nfun);
    free(bcopies);

    Expr* none_expr = expr_new_symbol(SYM_None);
    Expr* oargs[10] = {
        expr_new_integer(sopts->plot_points),
        expr_new_integer(sopts->max_recursion),
        expr_new_integer(sopts->max_plot_points),
        expr_new_integer(sopts->mesh ? 1 : 0),
        sopts->region_function ? expr_copy(sopts->region_function) : expr_copy(none_expr),
        sopts->exclusions ? expr_copy(sopts->exclusions) : expr_copy(none_expr),
        sopts->color_function ? expr_copy(sopts->color_function) : expr_copy(none_expr),
        expr_new_integer(sopts->color_function_scaling ? 1 : 0),
        sopts->filling ? expr_copy(sopts->filling) : expr_copy(none_expr),
        sopts->filling_style ? expr_copy(sopts->filling_style) : expr_copy(none_expr),
    };
    expr_free(none_expr);
    Expr* opts_list = expr_new_function(expr_new_symbol(SYM_List), oargs, 10);

    Expr* margs[3] = { expr_copy(var), blist, opts_list };
    return expr_new_function(expr_new_symbol(SYM_PlotResample), margs, 3);
}

static Expr* opt_or_none(Expr* e) {
    return (e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_None) ? NULL : e;
}

Expr* plot_resample(const Expr* graphics_expr, double xmin, double xmax,
                    double yclip_lo, double yclip_hi) {
    if (!graphics_expr || graphics_expr->type != EXPR_FUNCTION || !(xmin < xmax)) return NULL;

    /* Locate the $PlotResample[var, {bodies}, {opts...}] metadata arg. */
    const Expr* meta = NULL;
    size_t argc = graphics_expr->data.function.arg_count;
    for (size_t i = 1; i < argc; i++) {
        const Expr* a = graphics_expr->data.function.args[i];
        if (a && a->type == EXPR_FUNCTION && a->data.function.head
            && a->data.function.head->type == EXPR_SYMBOL
            && a->data.function.head->data.symbol == SYM_PlotResample
            && a->data.function.arg_count == 3) { meta = a; break; }
    }
    if (!meta) return NULL;

    Expr* var = meta->data.function.args[0];
    Expr* blist = meta->data.function.args[1];
    Expr* olist = meta->data.function.args[2];
    if (!var || var->type != EXPR_SYMBOL) return NULL;
    if (!blist || blist->type != EXPR_FUNCTION || !blist->data.function.head
        || blist->data.function.head->type != EXPR_SYMBOL
        || blist->data.function.head->data.symbol != SYM_List) return NULL;
    if (!olist || olist->type != EXPR_FUNCTION || olist->data.function.arg_count != 10) return NULL;

    PlotSampleOpts sopts = {
        .plot_points = 50, .max_recursion = 6, .max_plot_points = -1, .mesh = false,
        .region_function = NULL, .exclusions = NULL, .color_function = NULL,
        .color_function_scaling = true, .filling = NULL, .filling_style = NULL,
        /* Re-sample against the current zoom's visible band so detail tracks
         * what's on screen (a degenerate band from the caller disables it). */
        .yclip_lo = yclip_lo, .yclip_hi = yclip_hi,
    };
    Expr** o = olist->data.function.args;
    if (o[0]->type == EXPR_INTEGER && o[0]->data.integer >= 2) sopts.plot_points = (long)o[0]->data.integer;
    if (o[1]->type == EXPR_INTEGER && o[1]->data.integer >= 0) sopts.max_recursion = (int)o[1]->data.integer;
    if (o[2]->type == EXPR_INTEGER) sopts.max_plot_points = (long)o[2]->data.integer;
    if (o[3]->type == EXPR_INTEGER) sopts.mesh = (o[3]->data.integer != 0);
    sopts.region_function = opt_or_none(o[4]);
    sopts.exclusions = opt_or_none(o[5]);
    sopts.color_function = opt_or_none(o[6]);
    if (o[7]->type == EXPR_INTEGER) sopts.color_function_scaling = (o[7]->data.integer != 0);
    sopts.filling = opt_or_none(o[8]);
    sopts.filling_style = opt_or_none(o[9]);

    /* Re-sampling has no single-curve PlotStyle handy here; NULL just skips
     * the Filling restore-color step (the original draw's colour directive,
     * already in the static primitive list's defaults, stays in force). */
    return build_plot_primitives(blist->data.function.args, blist->data.function.arg_count,
                                 var, xmin, xmax, &sopts, NULL);
}

Expr* builtin_plot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    IterSpec ispec;
    if (!iter_spec_parse(spec, &ispec)) return NULL;
    if (ispec.kind != ITER_KIND_RANGE || !ispec.var) { iter_spec_free(&ispec); return NULL; }

    double xmin, xmax;
    if (!numericize_bound(ispec.imin, &xmin) || !numericize_bound(ispec.imax, &xmax) || !(xmin < xmax)) {
        iter_spec_free(&ispec);
        return NULL;
    }

    /* PlotLegends is consumed directly here (not via split_options/sopts)
     * since it only affects the separate $PlotLegendData metadata, not
     * sampling. Scan for it (and validate/evaluate it) before split_options
     * runs so a malformed value can still decline the whole call cleanly. */
    Expr* legends = NULL;
    for (size_t i = 2; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) continue;
        Expr* lhs = arg->data.function.args[0];
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol == SYM_PlotLegends) {
            legends = evaluate(expr_copy(arg->data.function.args[1]));
        }
    }

    PlotSampleOpts sopts;
    Expr** passthrough = NULL;
    size_t passthrough_count = 0;
    Expr* single_color = NULL;
    if (!split_options(res, &sopts, &passthrough, &passthrough_count, &single_color)) {
        iter_spec_free(&ispec);
        expr_free(legends);
        return NULL;
    }

    /* One curve (Plot[f, …]) or several (Plot[{f1, f2, …}, …]). A List head
     * on the first argument is the multi-curve form. */
    Expr* single = f;
    Expr** bodies;
    size_t nfun;
    bool is_list = (f->type == EXPR_FUNCTION && f->data.function.head
                    && f->data.function.head->type == EXPR_SYMBOL
                    && f->data.function.head->data.symbol == SYM_List);
    if (is_list) {
        nfun = f->data.function.arg_count;
        bodies = f->data.function.args; /* borrowed */
    } else {
        nfun = 1;
        bodies = &single;               /* borrowed */
    }

    if (nfun == 0) { /* Plot[{}, …] — nothing to draw */
        iter_spec_free(&ispec);
        expr_free(legends);
        expr_free(single_color);
        for (size_t i = 0; i < passthrough_count; i++) expr_free(passthrough[i]);
        free(passthrough);
        return NULL;
    }

    /* Sample every curve (shadowing the iterator variable) into the initial
     * primitive list shown at the home view. */
    Expr* prim_list = build_plot_primitives(bodies, nfun, ispec.var, xmin, xmax, &sopts, single_color);
    if (!prim_list) {
        iter_spec_free(&ispec);
        expr_free(legends);
        expr_free(single_color);
        for (size_t i = 0; i < passthrough_count; i++) expr_free(passthrough[i]);
        free(passthrough);
        return NULL;
    }

    Expr* legend_meta = build_legend_meta(legends, bodies, nfun, single_color);
    expr_free(legends);
    expr_free(single_color);

    /* Capture f/var/options so the renderer can re-sample at the current
     * zoom (see plot_resample) -- this is what keeps a magnified Sin[1/x^2]
     * smooth rather than revealing the original coarse grid. */
    Expr* meta = build_resample_meta(bodies, nfun, ispec.var, &sopts);
    iter_spec_free(&ispec);

    size_t gargc = 1 + passthrough_count + 1 + (legend_meta ? 1 : 0);
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = prim_list;
    for (size_t i = 0; i < passthrough_count; i++) gargs[1 + i] = passthrough[i];
    gargs[1 + passthrough_count] = meta;
    if (legend_meta) gargs[gargc - 1] = legend_meta;
    free(passthrough);

    Expr* graphics = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
    free(gargs);

    /* The REPL front end renders any top-level Graphics[...] result, so
     * we just return the object (no graphics_show here). See repl.c. */
    return graphics;
}
