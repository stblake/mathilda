/* parametricplot.c — ParametricPlot[body, {t, tmin, tmax}, opts...]
 *                  — ParametricPlot[body, {t, tmin, tmax}, {r, rmin, rmax}, opts...]
 *
 * HoldAll: the body and all iterator specs are unevaluated when received.
 *
 * One-iterator form:
 *   body must evaluate to {x, y} for each t. Single body or
 *   {body1, body2,...} (list of bodies) for multi-curve.
 *   Adaptive 2D sampling with a Euclidean chord-deviation flatness test.
 *   Output: Graphics[{Line[...], ...}, opts].
 *
 * Two-iterator form:
 *   body evaluates to {x, y} for each (t, r) pair.
 *   Samples a PlotPoints x PlotPoints grid and builds Polygon[] quads,
 *   like Plot3D but mapped back into the xy-plane.
 *   Output: Graphics[{Polygon[...], ...}, opts].
 *
 * In both forms the body can be any expression that evaluates to a
 * 2-element numeric list: a literal {fx, fy}, or a computed form such as
 * r^2 * {Sqrt[t] Cos[t], Sin[t]}, as long as the result has head List and
 * exactly two finite-real elements. */

#include "parametricplot.h"
#include "plot_common.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include <stdlib.h>
#include <math.h>

/* ---- Shared evaluation context ---- */

typedef struct {
    Expr* var1;             /* first (or only) iterator variable */
    Expr* var2;             /* second iterator variable; NULL for 1-iterator form */
    Expr* body;             /* expression that evaluates to {x, y} when vars are set */
    Expr* region_function;  /* borrowed; NULL = no filter */
} ParamEvalCtx;

/* Evaluate ctx->body (with vars already set in the symbol table) and
 * extract the resulting {x, y} pair.  Returns true on success.
 * The evaluator owns the result; we free it here. */
static bool eval_body_xy(ParamEvalCtx* ctx, double* x_out, double* y_out) {
    Expr* result = evaluate(ctx->body);
    bool ok = false;
    if (result && result->type == EXPR_FUNCTION
        && result->data.function.head->type == EXPR_SYMBOL
        && result->data.function.head->data.symbol.name == SYM_List
        && result->data.function.arg_count == 2) {
        double x, y;
        if (expr_to_real_double(result->data.function.args[0], &x) && isfinite(x)
            && expr_to_real_double(result->data.function.args[1], &y) && isfinite(y)) {
            if (!ctx->region_function || eval_region(ctx->region_function, x, y)) {
                *x_out = x; *y_out = y; ok = true;
            }
        }
    }
    expr_free(result);
    return ok;
}

/* Set var1 = t and evaluate the body (1-iterator form). */
static bool param_eval(double t, ParamEvalCtx* ctx, double* x_out, double* y_out) {
    symtab_add_own_value(ctx->var1->data.symbol.name, ctx->var1, expr_new_real(t));
    return eval_body_xy(ctx, x_out, y_out);
}

/* ---- 2D adaptive sampler (1-iterator form) ---- */

typedef struct {
    double t, x, y;
    bool   valid;
    bool   break_before;
} ParamPt;

#define PARAM_FLAT_TOL 0.0025

typedef struct {
    ParamPt* pts;
    size_t   len;
    size_t   cap;
    long     max_points;
    bool     pending_break;
} ParamBuf;

static void pbuf_push(ParamBuf* buf, double t, double x, double y) {
    if (buf->len == buf->cap) {
        buf->cap *= 2;
        buf->pts = realloc(buf->pts, sizeof(ParamPt) * buf->cap);
    }
    ParamPt* p = &buf->pts[buf->len++];
    p->t = t; p->x = x; p->y = y;
    p->valid = true;
    p->break_before = buf->pending_break;
    buf->pending_break = false;
}

static bool pbuf_capped(const ParamBuf* buf) {
    return buf->max_points > 0 && (long)buf->len >= buf->max_points;
}

static void param_subdivide(ParamEvalCtx* ctx, ParamBuf* buf,
                             double at, double ax, double ay, bool a_ok,
                             double bt, double bx, double by, bool b_ok,
                             int depth, int max_rec, double diag) {
    if (pbuf_capped(buf) || depth >= max_rec) {
        if (b_ok) pbuf_push(buf, bt, bx, by);
        else buf->pending_break = true;
        return;
    }

    double mt = (at + bt) / 2.0, mx = 0.0, my = 0.0;
    bool m_ok = param_eval(mt, ctx, &mx, &my);

    if (a_ok && b_ok && m_ok) {
        double lx = (ax + bx) / 2.0, ly = (ay + by) / 2.0;
        double dev = sqrt((mx-lx)*(mx-lx) + (my-ly)*(my-ly));

        /* Three-probe anti-aliasing (mirrors sampling.c). */
        double q1t = (at + mt) / 2.0, q3t = (mt + bt) / 2.0;
        double q1x = 0.0, q1y = 0.0, q3x = 0.0, q3y = 0.0;
        bool q1_ok = param_eval(q1t, ctx, &q1x, &q1y);
        bool q3_ok = param_eval(q3t, ctx, &q3x, &q3y);
        if (q1_ok && q3_ok) {
            double d1 = sqrt((q1x-(ax+0.25*(bx-ax)))*(q1x-(ax+0.25*(bx-ax)))
                           + (q1y-(ay+0.25*(by-ay)))*(q1y-(ay+0.25*(by-ay))));
            double d3 = sqrt((q3x-(ax+0.75*(bx-ax)))*(q3x-(ax+0.75*(bx-ax)))
                           + (q3y-(ay+0.75*(by-ay)))*(q3y-(ay+0.75*(by-ay))));
            if (d1 > dev) dev = d1;
            if (d3 > dev) dev = d3;
        }
        if (dev < PARAM_FLAT_TOL * diag) {
            pbuf_push(buf, bt, bx, by);
            return;
        }
    }

    param_subdivide(ctx, buf, at, ax, ay, a_ok, mt, mx, my, m_ok, depth+1, max_rec, diag);
    param_subdivide(ctx, buf, mt, mx, my, m_ok, bt, bx, by, b_ok, depth+1, max_rec, diag);
}

static ParamPt* param_sample(ParamEvalCtx* ctx,
                               double tmin, double tmax,
                               long plot_points, int max_rec, long max_pts,
                               size_t* out_count) {
    *out_count = 0;
    if (!ctx || !(tmin < tmax) || plot_points < 2) return NULL;

    double* gt  = malloc(sizeof(double) * (size_t)plot_points);
    double* gx  = malloc(sizeof(double) * (size_t)plot_points);
    double* gy  = malloc(sizeof(double) * (size_t)plot_points);
    bool*   gok = malloc(sizeof(bool)   * (size_t)plot_points);

    double xlo = 1e300, xhi = -1e300, ylo = 1e300, yhi = -1e300;
    for (long i = 0; i < plot_points; i++) {
        double t = tmin + (tmax - tmin) * (double)i / (double)(plot_points - 1);
        if (i == plot_points - 1) t = tmax;
        gt[i] = t;
        gok[i] = param_eval(t, ctx, &gx[i], &gy[i]);
        if (gok[i]) {
            if (gx[i] < xlo) xlo = gx[i];
            if (gx[i] > xhi) xhi = gx[i];
            if (gy[i] < ylo) ylo = gy[i];
            if (gy[i] > yhi) yhi = gy[i];
        }
    }
    double xspan = (xhi > xlo) ? (xhi - xlo) : 1.0;
    double yspan = (yhi > ylo) ? (yhi - ylo) : 1.0;
    double diag  = sqrt(xspan*xspan + yspan*yspan);
    if (diag < 1e-12) diag = 1.0;

    ParamBuf buf;
    buf.cap   = (size_t)plot_points * 4;
    buf.pts   = malloc(sizeof(ParamPt) * buf.cap);
    buf.len   = 0;
    buf.max_points   = max_pts;
    buf.pending_break = false;

    if (gok[0]) pbuf_push(&buf, gt[0], gx[0], gy[0]);
    else buf.pending_break = true;

    for (long i = 0; i + 1 < plot_points && !pbuf_capped(&buf); i++) {
        param_subdivide(ctx, &buf,
                         gt[i],   gx[i],   gy[i],   gok[i],
                         gt[i+1], gx[i+1], gy[i+1], gok[i+1],
                         0, max_rec, diag);
    }

    free(gt); free(gx); free(gy); free(gok);

    if (buf.len == 0) { free(buf.pts); return NULL; }
    *out_count = buf.len;
    return buf.pts;
}

/* ---- Option parsing ---- */

typedef struct {
    long  plot_points;
    int   max_recursion;
    long  max_plot_points;
    bool  mesh;
    Expr* region_function;  /* borrowed */
    Expr* color_function;   /* borrowed */
    bool  color_function_scaling;
} ParamSampleOpts;

/* `opts_start`: index of the first trailing Rule[] arg in res.
 * `default_plot_points`: 25 for 1-iterator (adaptive), 75 for 2-iterator (uniform grid).
 * For 1-iterator form: 2.  For 2-iterator form: 3. */
static bool split_options_param(Expr* res, size_t opts_start, long default_plot_points,
                                  ParamSampleOpts* sopts,
                                  Expr*** pass_out, size_t* pass_count_out,
                                  Expr** single_color_out) {
    sopts->plot_points            = default_plot_points;
    sopts->max_recursion          = 6;
    sopts->max_plot_points        = -1;
    sopts->mesh                   = false;
    sopts->region_function        = NULL;
    sopts->color_function         = NULL;
    sopts->color_function_scaling = true;
    *single_color_out = NULL;

    size_t argc = res->data.function.arg_count;
    size_t cap  = (argc > opts_start ? argc - opts_start : 0) + 4;
    Expr** pass = malloc(sizeof(Expr*) * cap);
    size_t n    = 0;

    bool have_axes = false, have_aspect = false, have_style = false, have_frame = false;

#define PFAIL() do { free(pass); expr_free(*single_color_out); *single_color_out = NULL; return false; } while (0)

    for (size_t i = opts_start; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) PFAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol.name : NULL;

        if (name == SYM_PlotPoints) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 2) PFAIL();
            sopts->plot_points = v;
        } else if (name == SYM_MaxRecursion) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 0) PFAIL();
            sopts->max_recursion = (int)v;
        } else if (name == SYM_MaxPlotPoints) {
            Expr* v = evaluate(expr_copy(rhs));
            bool is_inf = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Infinity);
            long lv = -1;
            bool ok = is_inf;
            if (!ok && v->type == EXPR_INTEGER && v->data.integer > 0) { lv = v->data.integer; ok = true; }
            expr_free(v);
            if (!ok) PFAIL();
            sopts->max_plot_points = lv;
        } else if (name == SYM_Mesh) {
            Expr* v = evaluate(expr_copy(rhs));
            bool on  = (v->type == EXPR_SYMBOL && (v->data.symbol.name == SYM_All  || v->data.symbol.name == SYM_True));
            bool off = (v->type == EXPR_SYMBOL && (v->data.symbol.name == SYM_None || v->data.symbol.name == SYM_False));
            expr_free(v);
            if (!on && !off) PFAIL();
            sopts->mesh = on;
        } else if (name == SYM_RegionFunction) {
            sopts->region_function = rhs;
        } else if (name == SYM_ColorFunction) {
            sopts->color_function = rhs;
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = evaluate(expr_copy(rhs));
            sopts->color_function_scaling = !(v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_False);
            expr_free(v);
        } else if (name == SYM_PlotStyle) {
            have_style = true;
            if (*single_color_out) expr_free(*single_color_out);
            *single_color_out = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), expr_copy(*single_color_out) };
            pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        } else if (name == SYM_AspectRatio) {
            have_aspect = true;
            if (rhs->type == EXPR_SYMBOL
                && (rhs->data.symbol.name == SYM_Automatic || rhs->data.symbol.name == SYM_Full)) {
                pass[n++] = expr_copy(arg);
            } else {
                double v;
                Expr* val = (numericize_bound(rhs, &v) && v > 0)
                    ? expr_new_real(v) : expr_new_symbol(SYM_Automatic);
                Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), val };
                pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
            }
        } else {
            if (name == SYM_Axes)  have_axes  = true;
            if (name == SYM_Frame) {
                if (!(rhs->type == EXPR_SYMBOL
                      && (rhs->data.symbol.name == SYM_False || rhs->data.symbol.name == SYM_None)))
                    have_frame = true;
            }
            Expr* val = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    if (!have_axes && !have_frame) {
        Expr* a[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_True) };
        pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_aspect) {
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_real(1.0) };
        pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_style) {
        Expr* rgb_a[3] = { expr_new_real(0.2), expr_new_real(0.4), expr_new_real(0.8) };
        Expr* rgb = expr_new_function(expr_new_symbol(SYM_RGBColor), rgb_a, 3);
        Expr* a[2] = { expr_new_symbol(SYM_PlotStyle), expr_copy(rgb) };
        pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        *single_color_out = rgb;
    }

    *pass_out       = pass;
    *pass_count_out = n;
    return true;
#undef PFAIL
}

/* ---- Primitive assembly: 1-iterator (Line) form ---- */

static Expr** build_param_curve(Expr* body, Expr* var1,
                                  double tmin, double tmax,
                                  const ParamSampleOpts* sopts,
                                  Expr* curve_color,
                                  size_t* out_count) {
    *out_count = 0;

    ParamEvalCtx ctx = { .var1 = var1, .var2 = NULL, .body = body,
                          .region_function = sopts->region_function };

    size_t npts;
    ParamPt* pts = param_sample(&ctx, tmin, tmax,
                                  sopts->plot_points, sopts->max_recursion,
                                  sopts->max_plot_points, &npts);
    if (!pts || npts == 0) { free(pts); return NULL; }

    Expr** prims = malloc(sizeof(Expr*) * (npts * 4 + 4));
    size_t pc = 0;
    size_t run_start = 0;

    for (size_t i = 1; i <= npts; i++) {
        bool end_of_run = (i == npts) || pts[i].break_before;
        if (!end_of_run) continue;

        size_t run_len = i - run_start;
        if (run_len >= 2) {
            if (sopts->color_function) {
                for (size_t j = 0; j + 1 < run_len; j++) {
                    size_t a = run_start + j, b = a + 1;
                    double tmid = (pts[a].t + pts[b].t) / 2.0;
                    prims[pc++] = eval_color_function(sopts->color_function,
                                                       tmid, 0.0,
                                                       tmin, tmax,
                                                       sopts->color_function_scaling);
                    Expr* pa[2] = { expr_new_real(pts[a].x), expr_new_real(pts[a].y) };
                    Expr* pb[2] = { expr_new_real(pts[b].x), expr_new_real(pts[b].y) };
                    Expr* seg[2] = {
                        expr_new_function(expr_new_symbol(SYM_List), pa, 2),
                        expr_new_function(expr_new_symbol(SYM_List), pb, 2),
                    };
                    Expr* seg_list = expr_new_function(expr_new_symbol(SYM_List), seg, 2);
                    Expr* largs[1] = { seg_list };
                    prims[pc++] = expr_new_function(expr_new_symbol(SYM_Line), largs, 1);
                }
            } else {
                Expr** line_pts = malloc(sizeof(Expr*) * run_len);
                for (size_t j = 0; j < run_len; j++) {
                    Expr* xy[2] = { expr_new_real(pts[run_start + j].x),
                                    expr_new_real(pts[run_start + j].y) };
                    line_pts[j] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
                }
                Expr* pts_list = expr_new_function(expr_new_symbol(SYM_List), line_pts, run_len);
                free(line_pts);
                Expr* largs[1] = { pts_list };
                prims[pc++] = expr_new_function(expr_new_symbol(SYM_Line), largs, 1);
            }
        }
        run_start = i;
    }

    if (sopts->mesh && npts > 0) {
        if (curve_color && !sopts->color_function) prims[pc++] = expr_copy(curve_color);
        double xlo = pts[0].x, xhi = pts[0].x, ylo = pts[0].y, yhi = pts[0].y;
        for (size_t j = 1; j < npts; j++) {
            if (pts[j].x < xlo) xlo = pts[j].x;
            if (pts[j].x > xhi) xhi = pts[j].x;
            if (pts[j].y < ylo) ylo = pts[j].y;
            if (pts[j].y > yhi) yhi = pts[j].y;
        }
        double diag = sqrt((xhi-xlo)*(xhi-xlo) + (yhi-ylo)*(yhi-ylo));
        double msize = (diag > 0) ? diag * 0.0025 : 0.005;
        Expr* ps_arg[1] = { expr_new_real(msize) };
        prims[pc++] = expr_new_function(expr_new_symbol(SYM_PointSize), ps_arg, 1);
        Expr** dots = malloc(sizeof(Expr*) * npts);
        for (size_t j = 0; j < npts; j++) {
            Expr* xy[2] = { expr_new_real(pts[j].x), expr_new_real(pts[j].y) };
            dots[j] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
        }
        Expr* dot_list = expr_new_function(expr_new_symbol(SYM_List), dots, npts);
        free(dots);
        Expr* pt_args[1] = { dot_list };
        prims[pc++] = expr_new_function(expr_new_symbol(SYM_Point), pt_args, 1);
    }

    free(pts);
    if (pc == 0) { free(prims); return NULL; }
    *out_count = pc;
    return prims;
}

/* ---- Primitive assembly: 2-iterator (Polygon) form ---- */

/* Grid point in (x,y)-space. */
typedef struct { double x, y; bool valid; } GridPt2D;

/* Helper: make a 2-element {x,y} List. */
static Expr* xy_pair(double x, double y) {
    Expr* a[2] = { expr_new_real(x), expr_new_real(y) };
    return expr_new_function(expr_new_symbol(SYM_List), a, 2);
}

/* Helper: make a Line[{{xa,ya},{xb,yb}}] between two grid points. */
static Expr* grid_line_2d(const GridPt2D* a, const GridPt2D* b) {
    Expr* pts[2] = { xy_pair(a->x, a->y), xy_pair(b->x, b->y) };
    Expr* pts_list = expr_new_function(expr_new_symbol(SYM_List), pts, 2);
    Expr* largs[1] = { pts_list };
    return expr_new_function(expr_new_symbol(SYM_Line), largs, 1);
}

/* Build a Polygon/Line primitive list from a 2-iterator ParametricPlot.
 * Samples an n×n grid in (t,r)-space, maps each point to (x,y) via body,
 * and emits one Polygon[] per valid 2×2 cell. */
static Expr** build_param_region(Expr* body, Expr* var1, Expr* var2,
                                   double t1min, double t1max,
                                   double t2min, double t2max,
                                   const ParamSampleOpts* sopts,
                                   Expr* color,
                                   size_t* out_count) {
    *out_count = 0;
    long n = sopts->plot_points;
    if (n < 2) return NULL;

    size_t ntotal = (size_t)(n * n);
    GridPt2D* grid = malloc(sizeof(GridPt2D) * ntotal);

    ParamEvalCtx ctx = { .var1 = var1, .var2 = var2, .body = body,
                          .region_function = sopts->region_function };

    for (long i = 0; i < n; i++) {
        double t = (n > 1) ? t1min + (t1max - t1min) * (double)i / (double)(n-1) : t1min;
        if (i == n-1) t = t1max;
        /* Set var1 once per row. */
        symtab_add_own_value(ctx.var1->data.symbol.name, ctx.var1, expr_new_real(t));
        for (long j = 0; j < n; j++) {
            double r = (n > 1) ? t2min + (t2max - t2min) * (double)j / (double)(n-1) : t2min;
            if (j == n-1) r = t2max;
            symtab_add_own_value(ctx.var2->data.symbol.name, ctx.var2, expr_new_real(r));
            GridPt2D* p = &grid[i * n + j];
            p->valid = eval_body_xy(&ctx, &p->x, &p->y);
        }
    }

    /* Worst-case allocation: one Polygon per cell + color + opacity + mesh. */
    long ncells = (n-1) * (n-1);
    long nmesh  = sopts->mesh ? 2 * (n-1) * n : 0;
    Expr** prims = malloc(sizeof(Expr*) * (size_t)(ncells + nmesh + 4));
    size_t pc = 0;

    /* Color directive up front (solid by default; pass PlotStyle -> {color, Opacity[a]}
     * to get a transparent region fill). */
    if (color && !sopts->color_function)
        prims[pc++] = expr_copy(color);

    /* Polygon quads. */
    for (long i = 0; i < n-1; i++) {
        for (long j = 0; j < n-1; j++) {
            GridPt2D* p00 = &grid[i * n + j];
            GridPt2D* p10 = &grid[(i+1) * n + j];
            GridPt2D* p11 = &grid[(i+1) * n + (j+1)];
            GridPt2D* p01 = &grid[i * n + (j+1)];

            if (!p00->valid || !p10->valid || !p11->valid || !p01->valid) continue;

            if (sopts->color_function) {
                /* Color by (t_mid, r_mid) passed as (x_scaled, y_scaled). */
                double tmid = t1min + (t1max - t1min) * (i + 0.5) / (double)(n-1);
                double rmid = t2min + (t2max - t2min) * (j + 0.5) / (double)(n-1);
                prims[pc++] = eval_color_function(sopts->color_function,
                                                   tmid, rmid,
                                                   t1min, t1max,
                                                   sopts->color_function_scaling);
            }

            Expr* verts[4] = {
                xy_pair(p00->x, p00->y),
                xy_pair(p10->x, p10->y),
                xy_pair(p11->x, p11->y),
                xy_pair(p01->x, p01->y),
            };
            Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, 4);
            Expr* pargs[1] = { vlist };
            prims[pc++] = expr_new_function(expr_new_symbol(SYM_Polygon), pargs, 1);
        }
    }

    /* Optional mesh: draw interior grid edges as Line[]. */
    if (sopts->mesh) {
        for (long i = 0; i < n-1; i++) {
            for (long j = 0; j < n; j++) {
                GridPt2D* pa = &grid[i * n + j];
                GridPt2D* pb = &grid[(i+1) * n + j];
                if (pa->valid && pb->valid) prims[pc++] = grid_line_2d(pa, pb);
            }
        }
        for (long i = 0; i < n; i++) {
            for (long j = 0; j < n-1; j++) {
                GridPt2D* pa = &grid[i * n + j];
                GridPt2D* pb = &grid[i * n + (j+1)];
                if (pa->valid && pb->valid) prims[pc++] = grid_line_2d(pa, pb);
            }
        }
    }

    free(grid);
    if (pc == 0) { free(prims); return NULL; }
    *out_count = pc;
    return prims;
}

/* ---- Helper: is this a 3-element iterator spec {var, min, max}? ---- */
static bool is_iterator(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List
        && e->data.function.arg_count == 3
        && e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* ---- Helpers ---- */

/* Pre-scan trailing Rule args starting at `opts_start` for PlotLegends.
 * Returns the evaluated value (caller owns), or NULL if absent. */
static Expr* prescan_plotlegends(Expr* res, size_t opts_start) {
    size_t argc = res->data.function.arg_count;
    for (size_t i = opts_start; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) continue;
        Expr* lhs = arg->data.function.args[0];
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_PlotLegends)
            return evaluate(expr_copy(arg->data.function.args[1]));
    }
    return NULL;
}

/* Assemble a Graphics[] expression from a primitives list, pass-through
 * options, and an optional legend metadata node. */
static Expr* make_graphics(Expr* prim_list, Expr** pass, size_t npass,
                             Expr* legend_meta) {
    size_t g_argc = 1 + npass + (legend_meta ? 1 : 0);
    Expr** g_args = malloc(sizeof(Expr*) * g_argc);
    g_args[0] = prim_list;
    for (size_t i = 0; i < npass; i++) g_args[1 + i] = pass[i];
    if (legend_meta) g_args[g_argc - 1] = legend_meta;
    Expr* g = expr_new_function(expr_new_symbol(SYM_Graphics), g_args, g_argc);
    free(g_args);
    return g;
}

/* ---- Entry point ---- */

Expr* builtin_parametricplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* arg 1 must be a {var, min, max} iterator. */
    Expr* iter1 = res->data.function.args[1];
    if (!is_iterator(iter1)) return NULL;

    Expr* var1 = iter1->data.function.args[0];
    double t1min, t1max;
    if (!numericize_bound(iter1->data.function.args[1], &t1min)
        || !numericize_bound(iter1->data.function.args[2], &t1max)
        || !(t1min < t1max)) return NULL;

    /* Detect 2-iterator form: arg 2 is also an iterator. */
    bool two_iter = (argc >= 3) && is_iterator(res->data.function.args[2]);
    size_t opts_start = two_iter ? 3 : 2;

    /* Extract PlotLegends before split_options_param (same pattern as Plot). */
    Expr* legends = prescan_plotlegends(res, opts_start);

    if (two_iter) {
        /* ---- 2-iterator form: filled region(s) ---- */
        Expr* iter2 = res->data.function.args[2];
        Expr* var2  = iter2->data.function.args[0];
        double t2min, t2max;
        if (!numericize_bound(iter2->data.function.args[1], &t2min)
            || !numericize_bound(iter2->data.function.args[2], &t2max)
            || !(t2min < t2max)) { expr_free(legends); return NULL; }

        Expr* body = res->data.function.args[0];

        /* Multi-surface: body is a List whose first element is also a List,
         * e.g. {{2r Cos[t], r Sin[t]}, {r Cos[t], 2r Sin[t]}}. Each sub-body
         * is rendered as a separate region in a distinct palette colour.
         * Single-surface: everything else ({fx,fy} or a computed expression). */
        bool body_is_list = (body->type == EXPR_FUNCTION
                             && body->data.function.head->type == EXPR_SYMBOL
                             && body->data.function.head->data.symbol.name == SYM_List);
        bool multi2 = (body_is_list
                       && body->data.function.arg_count >= 1
                       && body->data.function.args[0]->type == EXPR_FUNCTION
                       && body->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                       && body->data.function.args[0]->data.function.head->data.symbol.name == SYM_List);

        Expr** sub_bodies;
        size_t nsub;
        if (multi2) {
            sub_bodies = body->data.function.args;
            nsub       = body->data.function.arg_count;
        } else {
            sub_bodies = &body;
            nsub       = 1;
        }

        ParamSampleOpts sopts;
        Expr** pass;
        size_t npass;
        Expr* color = NULL;
        if (!split_options_param(res, 3, 75, &sopts, &pass, &npass, &color)) {
            expr_free(legends); return NULL;
        }

        Rule* old1 = iter_spec_shadow(var1);
        Rule* old2 = iter_spec_shadow(var2);

        size_t total = 0;
        Expr*** per_prims  = malloc(sizeof(Expr**) * nsub);
        size_t* per_counts = malloc(sizeof(size_t)  * nsub);
        bool any = false;

        for (size_t si = 0; si < nsub; si++) {
            /* build_param_region prepends the color directive internally, so
             * just pass palette_color(si) for multi or the single color for
             * single-surface -- no external prepend needed. */
            Expr* sc = multi2 ? palette_color(si) : color;
            per_prims[si] = build_param_region(sub_bodies[si], var1, var2,
                                                t1min, t1max, t2min, t2max,
                                                &sopts, sc, &per_counts[si]);
            if (multi2) expr_free(sc);
            if (per_counts[si] > 0) any = true;
            total += per_counts[si];
        }

        iter_spec_restore(var2, old2);
        iter_spec_restore(var1, old1);

        if (!any) {
            for (size_t si = 0; si < nsub; si++) free(per_prims[si]);
            free(per_prims); free(per_counts);
            free(pass); expr_free(color); expr_free(legends);
            return NULL;
        }

        Expr** prims = malloc(sizeof(Expr*) * total);
        size_t pc = 0;
        for (size_t si = 0; si < nsub; si++) {
            for (size_t j = 0; j < per_counts[si]; j++) prims[pc++] = per_prims[si][j];
            free(per_prims[si]);
        }
        free(per_prims); free(per_counts);

        Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, pc);
        free(prims);

        Expr* legend_meta = build_legend_meta(legends, sub_bodies, nsub, color);
        expr_free(legends);

        /* legend_meta ownership transfers into make_graphics → Graphics[...] */
        Expr* g = make_graphics(prim_list, pass, npass, legend_meta);
        free(pass);
        expr_free(color);
        return g;

    } else {
        /* ---- 1-iterator form: parametric curve(s) ---- */
        Expr* spec = res->data.function.args[0];

        /* Determine single-body vs. multi-body.
         * Multi: spec is a List whose first element is also a List.
         * Single: everything else (a literal {fx,fy}, or a computed body). */
        bool is_list_head = (spec->type == EXPR_FUNCTION
                             && spec->data.function.head->type == EXPR_SYMBOL
                             && spec->data.function.head->data.symbol.name == SYM_List);

        bool multi = (is_list_head
                      && spec->data.function.arg_count >= 1
                      && spec->data.function.args[0]->type == EXPR_FUNCTION
                      && spec->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                      && spec->data.function.args[0]->data.function.head->data.symbol.name == SYM_List);

        Expr** bodies;
        size_t nbodies;
        if (multi) {
            bodies  = spec->data.function.args;
            nbodies = spec->data.function.arg_count;
        } else {
            bodies  = &spec;    /* single body: the whole spec expression */
            nbodies = 1;
        }

        ParamSampleOpts sopts;
        Expr** pass;
        size_t npass;
        Expr* single_color = NULL;
        if (!split_options_param(res, 2, 25, &sopts, &pass, &npass, &single_color)) {
            expr_free(legends); return NULL;
        }

        Rule* old1 = iter_spec_shadow(var1);

        size_t total = 0;
        Expr*** per_prims  = malloc(sizeof(Expr**) * nbodies);
        size_t* per_counts = malloc(sizeof(size_t)  * nbodies);
        bool any = false;

        for (size_t ci = 0; ci < nbodies; ci++) {
            Expr* cc = multi ? palette_color(ci) : single_color;
            per_prims[ci] = build_param_curve(bodies[ci], var1, t1min, t1max,
                                               &sopts, cc, &per_counts[ci]);
            if (multi) expr_free(cc);
            if (per_counts[ci] > 0) any = true;
            total += per_counts[ci];
        }

        iter_spec_restore(var1, old1);

        if (!any) {
            for (size_t ci = 0; ci < nbodies; ci++) free(per_prims[ci]);
            free(per_prims); free(per_counts);
            free(pass); expr_free(single_color); expr_free(legends);
            return NULL;
        }

        size_t cap = total + (multi ? nbodies : 0) + 1;
        Expr** prims = malloc(sizeof(Expr*) * cap);
        size_t pc = 0;
        for (size_t ci = 0; ci < nbodies; ci++) {
            if (multi) prims[pc++] = palette_color(ci);
            for (size_t j = 0; j < per_counts[ci]; j++) prims[pc++] = per_prims[ci][j];
            free(per_prims[ci]);
        }
        free(per_prims); free(per_counts);

        Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, pc);
        free(prims);

        Expr* legend_meta = build_legend_meta(legends, bodies, nbodies, single_color);
        expr_free(legends);

        /* legend_meta ownership transfers into make_graphics → Graphics[...] */
        Expr* g = make_graphics(prim_list, pass, npass, legend_meta);
        free(pass);
        expr_free(single_color);
        return g;
    }
}
