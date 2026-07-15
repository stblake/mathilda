/* streamplot.c — StreamPlot[{vx,vy}, {x,xmin,xmax}, {y,ymin,ymax}, opts...]
 *
 * Traces streamlines of a 2-D vector field by RK4 integration from a grid
 * of seed points, emitting one Arrow[...] primitive per stream.  The result
 * is returned as a Graphics[...] object (auto-displayed by the REPL). */

#include "streamplot.h"
#include "plot_common.h"
#include "show.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Forward declaration so we can use it inside option parsing          */
/* ------------------------------------------------------------------ */

static bool numericize(Expr* e, double* out) { return numericize_bound(e, out); }

/* ------------------------------------------------------------------ */
/* Vector field evaluation                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    Expr* xvar;  /* borrowed */
    Expr* yvar;  /* borrowed */
    Expr* vx;    /* borrowed */
    Expr* vy;    /* borrowed */
} FieldCtx;

/* Evaluate the vector field at (x, y), storing results in (*vx_out, *vy_out).
 * Returns false if either component is non-finite or fails to evaluate. */
static bool eval_field(const FieldCtx* ctx, double x, double y,
                       double* vx_out, double* vy_out) {
    Expr* xval = expr_new_real(x);
    Expr* yval = expr_new_real(y);
    symtab_add_own_value(ctx->xvar->data.symbol, ctx->xvar, xval);
    symtab_add_own_value(ctx->yvar->data.symbol, ctx->yvar, yval);

    Expr* rx = evaluate(ctx->vx);
    Expr* ry = evaluate(ctx->vy);

    double dvx, dvy;
    bool ok = expr_to_real_double(rx, &dvx) && isfinite(dvx)
           && expr_to_real_double(ry, &dvy) && isfinite(dvy);
    expr_free(rx);
    expr_free(ry);
    expr_free(xval);
    expr_free(yval);

    if (ok) { *vx_out = dvx; *vy_out = dvy; }
    return ok;
}

/* ------------------------------------------------------------------ */
/* Color function evaluation                                           */
/* ------------------------------------------------------------------ */

/* Evaluate StreamColorFunction at the stream's midpoint.
 * Tries: f[x,y,vx,vy,speed], f[x,y,vx,vy], f[x,y], f[speed],
 * and "Rainbow" (hue = scaled speed). Returns NULL if not set. */
static Expr* eval_stream_color(Expr* cfn,
                               double x, double y,
                               double vx, double vy,
                               double speed,
                               double speed_min, double speed_max) {
    if (!cfn) return NULL;

    if (cfn->type == EXPR_STRING) {
        double t = (speed_max > speed_min)
            ? (speed - speed_min) / (speed_max - speed_min) : 0.5;
        Expr* c = named_color_ramp(cfn->data.string, t);
        if (c) return c;
    }

    /* Try multi-arg forms, falling back to fewer args.
     * f[x,y,vx,vy,speed] → f[x,y,vx,vy] → f[x,y] → f[speed] */
    static const int arities[] = { 5, 4, 2, 1 };
    double vals[5] = { x, y, vx, vy, speed };
    for (size_t ai = 0; ai < 4; ai++) {
        int arity = arities[ai];
        Expr** fargs = malloc(sizeof(Expr*) * (size_t)arity);
        for (int j = 0; j < arity; j++) fargs[j] = expr_new_real(vals[j]);
        Expr* fcall = expr_new_function(expr_copy(cfn), fargs, (size_t)arity);
        free(fargs);
        Expr* result = evaluate(fcall);
        bool is_color = (result->type == EXPR_FUNCTION
                         && result->data.function.head
                         && result->data.function.head->type == EXPR_SYMBOL
                         && (result->data.function.head->data.symbol == SYM_RGBColor
                          || result->data.function.head->data.symbol == SYM_GrayLevel
                          || result->data.function.head->data.symbol == SYM_Hue));
        if (is_color) return result;
        expr_free(result);
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* RK4 streamline integration                                          */
/* ------------------------------------------------------------------ */

typedef struct { double x, y; } Point2;

/* Integrate one streamline forward from (x0, y0) for at most max_steps
 * steps of size h, staying within [xmin,xmax] x [ymin,ymax].
 * Reverses direction when forward==false.
 * Returns a malloc'd array of Point2; *n_out is the count (including seed).
 * The array always starts at {x0, y0}. */
static Point2* rk4_integrate(const FieldCtx* ctx,
                              double x0, double y0,
                              double xmin, double xmax,
                              double ymin, double ymax,
                              double h,
                              int max_steps,
                              bool forward,
                              Expr* region_fn,
                              size_t* n_out) {
    size_t cap = (size_t)max_steps + 2;
    Point2* pts = malloc(sizeof(Point2) * cap);
    size_t n = 0;
    pts[n++] = (Point2){ x0, y0 };

    double sign = forward ? 1.0 : -1.0;
    double x = x0, y = y0;

    for (int step = 0; step < max_steps; step++) {
        double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;
        if (!eval_field(ctx, x, y, &k1x, &k1y)) break;
        if (!eval_field(ctx, x + 0.5 * h * sign * k1x,
                              y + 0.5 * h * sign * k1y, &k2x, &k2y)) break;
        if (!eval_field(ctx, x + 0.5 * h * sign * k2x,
                              y + 0.5 * h * sign * k2y, &k3x, &k3y)) break;
        if (!eval_field(ctx, x + h * sign * k3x,
                              y + h * sign * k3y, &k4x, &k4y)) break;

        double nx = x + sign * h * (k1x + 2*k2x + 2*k3x + k4x) / 6.0;
        double ny = y + sign * h * (k1y + 2*k2y + 2*k3y + k4y) / 6.0;

        /* Clamp: stop if we leave the domain. */
        if (nx < xmin || nx > xmax || ny < ymin || ny > ymax) break;

        /* RegionFunction check. */
        if (region_fn && !eval_region(region_fn, nx, ny)) break;

        x = nx; y = ny;
        if (n >= cap) { cap *= 2; pts = realloc(pts, sizeof(Point2) * cap); }
        pts[n++] = (Point2){ x, y };
    }

    *n_out = n;
    return pts;
}

/* ------------------------------------------------------------------ */
/* Seed placement and proximity check                                  */
/* ------------------------------------------------------------------ */

typedef struct {
    double* xs;
    double* ys;
    size_t  n;
    size_t  cap;
} SeedList;

static void seed_push(SeedList* sl, double x, double y) {
    if (sl->n >= sl->cap) {
        sl->cap = sl->cap ? sl->cap * 2 : 64;
        sl->xs = realloc(sl->xs, sizeof(double) * sl->cap);
        sl->ys = realloc(sl->ys, sizeof(double) * sl->cap);
    }
    sl->xs[sl->n] = x;
    sl->ys[sl->n] = y;
    sl->n++;
}

/* ------------------------------------------------------------------ */
/* Option parsing                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    int    stream_points;     /* seeds per axis; default 15 */
    double stream_scale;      /* max arc-length / domain diag; <=0 = None (run full) */
    Expr*  stream_style;      /* borrowed list of directives, or NULL */
    Expr*  color_function;    /* borrowed; NULL = none */
    Expr*  region_function;   /* borrowed; NULL = none */
    Expr*  legends;           /* borrowed; NULL = none */
    ScaleFnType sf_x, sf_y;   /* ScalingFunctions per-axis */
} StreamOpts;

static bool split_stream_options(Expr* res, StreamOpts* so,
                                  Expr*** passthrough_out, size_t* pt_count_out) {
    so->stream_points = 15;
    so->stream_scale  = 0.08; /* 8% of domain diagonal by default */
    so->stream_style  = NULL;
    so->color_function = NULL;
    so->region_function = NULL;
    so->legends = NULL;
    so->sf_x = SF_NONE;
    so->sf_y = SF_NONE;

    size_t argc = res->data.function.arg_count;
    size_t cap = (argc > 3 ? argc - 3 : 0) + 4;
    Expr** passthrough = malloc(sizeof(Expr*) * cap);
    size_t n = 0;

    bool have_axes = false, have_aspect = false, have_frame = false;

#define SP_FAIL() do { free(passthrough); return false; } while(0)

    for (size_t i = 3; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) SP_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol : NULL;

        if (name == SYM_StreamPoints || name == SYM_PlotPoints) {
            Expr* v = evaluate(expr_copy(rhs));
            if (v->type == EXPR_INTEGER && v->data.integer >= 1) {
                so->stream_points = (int)v->data.integer;
                expr_free(v);
            } else {
                /* Automatic / Fine / Coarse keywords */
                expr_free(v);
                /* keep default */
            }
        } else if (name == SYM_StreamScale) {
            Expr* v = evaluate(expr_copy(rhs));
            if (v->type == EXPR_SYMBOL && v->data.symbol == SYM_None) {
                so->stream_scale = -1.0; /* run until boundary */
            } else if (v->type == EXPR_SYMBOL && v->data.symbol == SYM_Automatic) {
                so->stream_scale = 0.08;
            } else {
                double sv;
                if (numericize(v, &sv) && sv > 0) so->stream_scale = sv;
            }
            expr_free(v);
        } else if (name == SYM_StreamStyle) {
            so->stream_style = rhs; /* borrowed */
        } else if (name == SYM_StreamColorFunction || name == SYM_ColorFunction) {
            so->color_function = rhs; /* borrowed */
        } else if (name == SYM_RegionFunction) {
            so->region_function = rhs; /* borrowed */
        } else if (name == SYM_PlotLegends) {
            so->legends = rhs; /* borrowed */
        } else if (name == SYM_ScalingFunctions) {
            Expr* v = evaluate(expr_copy(rhs));
            parse_scaling_functions(v, &so->sf_x, &so->sf_y);
            expr_free(v);
        } else {
            /* Pass through to Graphics[...] */
            if (name == SYM_Axes)         have_axes   = true;
            else if (name == SYM_AspectRatio) have_aspect = true;
            else if (name == SYM_Frame) {
                if (!(rhs->type == EXPR_SYMBOL
                      && (rhs->data.symbol == SYM_False || rhs->data.symbol == SYM_None)))
                    have_frame = true;
            }
            Expr* val = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    /* StreamPlot defaults: Frame->True, Axes->False, square aspect. */
    if (!have_axes && !have_frame) {
        Expr* fa[2] = { expr_new_symbol(SYM_Frame), expr_new_symbol(SYM_True) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), fa, 2);
        Expr* aa[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_False) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), aa, 2);
    } else if (!have_axes) {
        Expr* aa[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_False) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), aa, 2);
    }
    if (!have_aspect) {
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_integer(1) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }

    *passthrough_out = passthrough;
    *pt_count_out = n;
    return true;

#undef SP_FAIL
}

/* ------------------------------------------------------------------ */
/* Build Line[...] and Arrow[...] Exprs from a Point2 array            */
/* ------------------------------------------------------------------ */

/* Build Arrow[{{x1,y1},...,{xn,yn}}] from the full point array. */
static Expr* make_arrow(const Point2* pts, size_t n) {
    if (n < 2) return NULL;
    Expr** coord_exprs = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* xy[2] = { expr_new_real(pts[i].x), expr_new_real(pts[i].y) };
        coord_exprs[i] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
    }
    Expr* pt_list = expr_new_function(expr_new_symbol(SYM_List), coord_exprs, n);
    free(coord_exprs);
    Expr* args[1] = { pt_list };
    return expr_new_function(expr_new_symbol(SYM_Arrow), args, 1);
}

/* Build Line[{{x1,y1},...,{xn,yn}}] — the streamline shaft. */
static Expr* make_line(const Point2* pts, size_t n) {
    if (n < 2) return NULL;
    Expr** coord_exprs = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* xy[2] = { expr_new_real(pts[i].x), expr_new_real(pts[i].y) };
        coord_exprs[i] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
    }
    Expr* pt_list = expr_new_function(expr_new_symbol(SYM_List), coord_exprs, n);
    free(coord_exprs);
    Expr* args[1] = { pt_list };
    return expr_new_function(expr_new_symbol(SYM_Line), args, 1);
}

/* Build a fixed-length Arrow[{p1, p2}] centered at the stream midpoint,
 * oriented along the local flow direction.  `len` is in plot coordinates and
 * is the same for every stream, so arrows near tight attractors (where the
 * stream is very short) are the same size as arrows in open regions. */
static Expr* make_fixed_arrow(const Point2* pts, size_t n, double len) {
    if (n < 2 || len <= 0) return NULL;
    size_t mid = n / 2;
    if (mid == 0) mid = 1;

    /* Local flow direction at the midpoint. */
    double dx = pts[mid].x - pts[mid - 1].x;
    double dy = pts[mid].y - pts[mid - 1].y;
    double d  = sqrt(dx * dx + dy * dy);
    if (d < 1e-30) return NULL;
    double ux = dx / d, uy = dy / d;

    double cx = pts[mid].x, cy = pts[mid].y;
    double half = len * 0.5;
    Point2 seg[2] = {
        { cx - half * ux, cy - half * uy },
        { cx + half * ux, cy + half * uy },
    };
    return make_arrow(seg, 2);
}

/* Default speed-to-color: dark blue-purple (low) → yellow (high), matching
 * Mathematica's default StreamPlot thermal palette. */
static Expr* default_stream_color(double spd, double spd_min, double spd_max) {
    double t = (spd_max > spd_min + 1e-30)
        ? (spd - spd_min) / (spd_max - spd_min) : 0.5;
    double rv, gv, bv;
    thermal_rgb(t, &rv, &gv, &bv);
    Expr* args[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), args, 3);
}

/* ------------------------------------------------------------------ */
/* Main builtin                                                         */
/* ------------------------------------------------------------------ */

Expr* builtin_streamplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 3) return NULL;

    /* Arg 0: {vx, vy} — must be a 2-element List (held, not yet evaluated). */
    Expr* field_arg = res->data.function.args[0];
    if (!field_arg || field_arg->type != EXPR_FUNCTION
        || !field_arg->data.function.head
        || field_arg->data.function.head->type != EXPR_SYMBOL
        || field_arg->data.function.head->data.symbol != SYM_List
        || field_arg->data.function.arg_count != 2)
        return NULL;

    Expr* vx_body = field_arg->data.function.args[0]; /* borrowed */
    Expr* vy_body = field_arg->data.function.args[1]; /* borrowed */

    /* Arg 1: {x, xmin, xmax} */
    IterSpec xspec;
    if (!iter_spec_parse(res->data.function.args[1], &xspec)) return NULL;
    if (xspec.kind != ITER_KIND_RANGE || !xspec.var) { iter_spec_free(&xspec); return NULL; }
    double xmin, xmax;
    if (!numericize(xspec.imin, &xmin) || !numericize(xspec.imax, &xmax) || !(xmin < xmax)) {
        iter_spec_free(&xspec); return NULL;
    }

    /* Arg 2: {y, ymin, ymax} */
    IterSpec yspec;
    if (!iter_spec_parse(res->data.function.args[2], &yspec)) { iter_spec_free(&xspec); return NULL; }
    if (yspec.kind != ITER_KIND_RANGE || !yspec.var) { iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL; }
    double ymin, ymax;
    if (!numericize(yspec.imin, &ymin) || !numericize(yspec.imax, &ymax) || !(ymin < ymax)) {
        iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL;
    }

    /* Options */
    StreamOpts so;
    Expr** passthrough = NULL;
    size_t pt_count = 0;
    if (!split_stream_options(res, &so, &passthrough, &pt_count)) {
        iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL;
    }

    /* World-space (scaled) axis bounds — integration happens in data space,
     * but seed placement and PlotRange use world space. */
    double u_xmin = scale_apply(so.sf_x, xmin), u_xmax = scale_apply(so.sf_x, xmax);
    double u_ymin = scale_apply(so.sf_y, ymin), u_ymax = scale_apply(so.sf_y, ymax);

    /* Always embed an explicit PlotRange (world coords) so the renderer uses
     * the user's domain.  Skip if the user already passed their own PlotRange. */
    {
        bool have_pr = false;
        for (size_t i = 0; i < pt_count; i++) {
            const Expr* e = passthrough[i];
            if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
                && e->data.function.args[0]->type == EXPR_SYMBOL
                && e->data.function.args[0]->data.symbol == SYM_PlotRange)
                { have_pr = true; break; }
        }
        if (!have_pr) {
            passthrough = realloc(passthrough, sizeof(Expr*) * (pt_count + 1));
            Expr* xr[2] = { expr_new_real(u_xmin), expr_new_real(u_xmax) };
            Expr* yr[2] = { expr_new_real(u_ymin), expr_new_real(u_ymax) };
            Expr* xrange = expr_new_function(expr_new_symbol(SYM_List), xr, 2);
            Expr* yrange = expr_new_function(expr_new_symbol(SYM_List), yr, 2);
            Expr* rlist[2] = { xrange, yrange };
            Expr* pr = expr_new_function(expr_new_symbol(SYM_List), rlist, 2);
            Expr* ra[2] = { expr_new_symbol(SYM_PlotRange), pr };
            passthrough[pt_count++] =
                expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
    }

    emit_scaling_meta(so.sf_x, so.sf_y, &passthrough, &pt_count);

    /* Shadow both iterator variables so the field sees no stray OwnValues. */
    Rule* old_x = iter_spec_shadow(xspec.var);
    Rule* old_y = iter_spec_shadow(yspec.var);

    FieldCtx ctx = { .xvar = xspec.var, .yvar = yspec.var,
                     .vx = vx_body, .vy = vy_body };

    /* Integration parameters in world space (for consistent arrow sizing). */
    double w_dx = u_xmax - u_xmin, w_dy = u_ymax - u_ymin;
    double diag = sqrt(w_dx * w_dx + w_dy * w_dy);
    double h = diag / (so.stream_points * 20.0);
    if (h <= 0) h = 1e-4;

    /* Uniform arrow indicator length: ~70% of seed spacing so every direction
     * chevron is the same size regardless of local stream arc length. */
    double seed_spacing = diag / so.stream_points;
    double arrow_len = seed_spacing * 0.70;

    /* max_steps: enough to cross the whole domain from corner to corner,
     * unless StreamScale limits arc-length. */
    double max_arc = (so.stream_scale > 0) ? so.stream_scale * diag : diag * 4.0;
    int max_steps = (int)(max_arc / h) + 1;
    if (max_steps > 4000) max_steps = 4000;

    /* Build seed list from a uniform grid in world space.  Convert seeds to
     * data space for field evaluation; RegionFunction uses data coords. */
    int np = so.stream_points;
    double du_x = (u_xmax - u_xmin) / np;
    double du_y = (u_ymax - u_ymin) / np;
    SeedList seeds = { NULL, NULL, 0, 0 };
    for (int ix = 0; ix < np; ix++) {
        double sx = scale_invert(so.sf_x, u_xmin + (ix + 0.5) * du_x);
        for (int iy = 0; iy < np; iy++) {
            double sy = scale_invert(so.sf_y, u_ymin + (iy + 0.5) * du_y);
            if (so.region_function && !eval_region(so.region_function, sx, sy)) continue;
            /* Quick sanity: skip seeds where the field is zero (or fails). */
            double tvx, tvy;
            if (!eval_field(&ctx, sx, sy, &tvx, &tvy)) continue;
            seed_push(&seeds, sx, sy);
        }
    }

    /* Collect speed samples to scale ColorFunction if requested. */
    /* We do two passes: first integrate all streams, then colorize. */
    size_t max_streams = seeds.n;
    Point2** all_streams = malloc(sizeof(Point2*) * (max_streams > 0 ? max_streams : 1));
    size_t*  all_lengths = malloc(sizeof(size_t)  * (max_streams > 0 ? max_streams : 1));
    double*  all_speed   = malloc(sizeof(double)  * (max_streams > 0 ? max_streams : 1));
    size_t nstreams = 0;

    for (size_t si = 0; si < seeds.n; si++) {
        size_t n_pts;
        Point2* pts = rk4_integrate(&ctx, seeds.xs[si], seeds.ys[si],
                                    xmin, xmax, ymin, ymax, h, max_steps,
                                    true, so.region_function, &n_pts);
        if (n_pts < 2) { free(pts); continue; }

        /* Midpoint speed for colorization (in data space, before scaling). */
        size_t mid = n_pts / 2;
        double svx, svy;
        double spd = 0.0;
        if (eval_field(&ctx, pts[mid].x, pts[mid].y, &svx, &svy))
            spd = sqrt(svx * svx + svy * svy);

        /* Transform stream points from data space to world space. */
        if (so.sf_x != SF_NONE || so.sf_y != SF_NONE) {
            for (size_t k = 0; k < n_pts; k++) {
                pts[k].x = scale_apply(so.sf_x, pts[k].x);
                pts[k].y = scale_apply(so.sf_y, pts[k].y);
            }
        }

        all_streams[nstreams] = pts;
        all_lengths[nstreams] = n_pts;
        all_speed[nstreams]   = spd;
        nstreams++;
    }

    free(seeds.xs); free(seeds.ys);

    /* Restore iterator variable bindings. */
    iter_spec_restore(xspec.var, old_x);
    iter_spec_restore(yspec.var, old_y);

    if (nstreams == 0) {
        free(all_streams); free(all_lengths); free(all_speed);
        iter_spec_free(&xspec); iter_spec_free(&yspec);
        for (size_t i = 0; i < pt_count; i++) expr_free(passthrough[i]);
        free(passthrough);
        return NULL;
    }

    /* Speed range for ColorFunction scaling. */
    double spd_min = all_speed[0], spd_max = all_speed[0];
    for (size_t i = 1; i < nstreams; i++) {
        if (all_speed[i] < spd_min) spd_min = all_speed[i];
        if (all_speed[i] > spd_max) spd_max = all_speed[i];
    }

    /* ---- Build primitive list ---- */
    /* Per-stream layout: [color?, Line[all_pts], Arrow[mid_2pt]]
     *
     * color source (in priority order):
     *   1. StreamColorFunction / ColorFunction  — explicit user function
     *   2. speed-based default (blue-violet → orange) — when no StreamStyle
     *   3. StreamStyle directive              — single global style, no per-stream color
     *
     * The Arrow is a tiny 2-point segment at 40–60% of the stream arc so it
     * acts as a small direction indicator mid-stream rather than a large cap. */
    bool have_style = (so.stream_style != NULL);
    bool have_cfn   = (so.color_function != NULL);
    bool use_speed_color = !have_style && !have_cfn;

    /* 4 primitives per stream (color + Line + Arrow + possible re-evaluation overhead). */
    size_t prim_cap = nstreams * 4 + 8;
    Expr** prims = malloc(sizeof(Expr*) * prim_cap);
    size_t nprim = 0;

    /* Global thickness directive (thin line weight for all streams). */
    {
        Expr* thick_arg[1] = { expr_new_real(0.006) };
        prims[nprim++] = expr_new_function(expr_new_symbol(SYM_Thickness), thick_arg, 1);
    }

    /* If StreamStyle given, apply it once before all streams. */
    if (have_style) {
        prims[nprim++] = evaluate(expr_copy(so.stream_style));
    }

    for (size_t i = 0; i < nstreams; i++) {
        Point2* pts   = all_streams[i];
        size_t  n_pts = all_lengths[i];
        double  spd   = all_speed[i];

        /* Per-stream color. */
        Expr* col = NULL;
        if (have_cfn) {
            size_t mid = n_pts / 2;
            double svx = 0, svy = 0;
            Rule* sx2 = iter_spec_shadow(xspec.var);
            Rule* sy2 = iter_spec_shadow(yspec.var);
            eval_field(&ctx, pts[mid].x, pts[mid].y, &svx, &svy);
            iter_spec_restore(xspec.var, sx2);
            iter_spec_restore(yspec.var, sy2);
            col = eval_stream_color(so.color_function,
                                    pts[mid].x, pts[mid].y,
                                    svx, svy, spd, spd_min, spd_max);
        } else if (use_speed_color) {
            col = default_stream_color(spd, spd_min, spd_max);
        }

        if (col) {
            if (nprim + 3 >= prim_cap) { prim_cap *= 2; prims = realloc(prims, sizeof(Expr*) * prim_cap); }
            prims[nprim++] = col;
        }

        /* Full streamline as a smooth Line[...]. */
        Expr* shaft = make_line(pts, n_pts);
        if (shaft) {
            if (nprim + 2 >= prim_cap) { prim_cap *= 2; prims = realloc(prims, sizeof(Expr*) * prim_cap); }
            prims[nprim++] = shaft;
        }

        /* Fixed-length direction Arrow at the stream midpoint. */
        Expr* mid_arrow = make_fixed_arrow(pts, n_pts, arrow_len);
        free(pts);
        if (mid_arrow) {
            if (nprim + 1 >= prim_cap) { prim_cap *= 2; prims = realloc(prims, sizeof(Expr*) * prim_cap); }
            prims[nprim++] = mid_arrow;
        }
    }
    free(all_streams); free(all_lengths); free(all_speed);

    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, nprim);
    free(prims);

    /* Color-scale bar metadata: $StreamColorBar[spd_min, spd_max].
     * Only emitted when PlotLegends is set and speed-based coloring is active.
     * Replaces the discrete swatch legend used by Plot — StreamPlot's legend
     * is a continuous gradient bar, not a labeled-line list. */
    Expr* color_bar_meta = NULL;
    if (so.legends && use_speed_color) {
        Expr* eval_legends = evaluate(expr_copy(so.legends));
        bool want_legend = !(eval_legends->type == EXPR_SYMBOL
                             && (eval_legends->data.symbol == SYM_None
                              || eval_legends->data.symbol == SYM_False));
        expr_free(eval_legends);
        if (want_legend) {
            Expr* cfn_copy = so.color_function
                             ? expr_copy(so.color_function)
                             : expr_new_symbol(SYM_Automatic);
            Expr* cb_args[3] = { expr_new_real(spd_min), expr_new_real(spd_max), cfn_copy };
            color_bar_meta = expr_new_function(
                expr_new_symbol(SYM_StreamColorBar), cb_args, 3);
        }
    }

    /* Assemble Graphics[primitives, opts..., colorBarMeta?]. */
    size_t gargc = 1 + pt_count + (color_bar_meta ? 1 : 0);
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = prim_list;
    for (size_t i = 0; i < pt_count; i++) gargs[1 + i] = passthrough[i];
    if (color_bar_meta) gargs[gargc - 1] = color_bar_meta;
    free(passthrough);

    Expr* graphics = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
    free(gargs);

    iter_spec_free(&xspec);
    iter_spec_free(&yspec);

    return graphics;
}
