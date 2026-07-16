/* parametricplot3d.c — ParametricPlot3D[body, {t, tmin, tmax}, opts...]
 *                     — ParametricPlot3D[body, {t, tmin, tmax}, {u, umin, umax}, opts...]
 *
 * HoldAll: the body and all iterator specs are unevaluated when received.
 *
 * One-iterator form:
 *   body must evaluate to {x, y, z} for each t. Single body or
 *   {body1, body2,...} (list of bodies) for multi-curve.
 *   Adaptive 3D sampling with a Euclidean chord-deviation test in (x,y,z) space.
 *   Three interior probes per interval prevent aliasing against periodic curves.
 *   Output: Graphics3D[{Line[...], ...}, opts].
 *
 * Two-iterator form:
 *   body evaluates to {x, y, z} for each (t, u) pair.
 *   Samples a PlotPoints x PlotPoints grid and builds Polygon[] quads,
 *   like Plot3D but driven by a parametric mapping instead of f(x,y).
 *   Output: Graphics3D[{Polygon[...], ...}, opts].
 *
 * In both forms the body can be any expression that evaluates to a
 * 3-element numeric list: a literal {fx, fy, fz}, or a computed form,
 * as long as the result has head List and exactly three finite-real elements.
 * ColorFunction receives scaled (x,y,z) coordinates (spatial, not parameter),
 * identical to Plot3D's convention. "Rainbow" sweeps hue over the z-extent. */

#include "parametricplot3d.h"
#include "plot_common.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include <stdlib.h>
#include <math.h>

/* ---- 3D region function: f[x,y,z] tried first, then f[x,y] fallback ---- */

static bool eval_region3d(Expr* rfn, double x, double y, double z) {
    Expr* args[3] = { expr_new_real(x), expr_new_real(y), expr_new_real(z) };
    Expr* call = expr_new_function(expr_copy(rfn), args, 3);
    Expr* r = evaluate(call);
    expr_free(call);
    bool is_true  = (r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True);
    bool is_false = (r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_False);
    expr_free(r);
    if (is_true)  return true;
    if (is_false) return false;
    return eval_region(rfn, x, y);
}

/* ---- Body evaluation context ---- */

typedef struct {
    Expr* var1;             /* first (or only) iterator variable */
    Expr* var2;             /* second iterator variable; NULL for 1-iterator form */
    Expr* body;             /* expression that evaluates to {x, y, z} */
    Expr* region_function;  /* borrowed; NULL = no filter */
} Param3DEvalCtx;

/* Evaluate ctx->body (with vars already set in symtab) and extract {x,y,z}.
 * Returns true on success and stores the coordinates in xo, yo, zo. */
static bool eval_body_xyz(Param3DEvalCtx* ctx, double* xo, double* yo, double* zo) {
    Expr* result = evaluate(ctx->body);
    bool ok = false;
    if (result && result->type == EXPR_FUNCTION
        && result->data.function.head->type == EXPR_SYMBOL
        && result->data.function.head->data.symbol.name == SYM_List
        && result->data.function.arg_count == 3) {
        double x, y, z;
        if (expr_to_real_double(result->data.function.args[0], &x) && isfinite(x)
            && expr_to_real_double(result->data.function.args[1], &y) && isfinite(y)
            && expr_to_real_double(result->data.function.args[2], &z) && isfinite(z)) {
            if (!ctx->region_function || eval_region3d(ctx->region_function, x, y, z)) {
                *xo = x; *yo = y; *zo = z; ok = true;
            }
        }
    }
    expr_free(result);
    return ok;
}

/* Set var1 = t and evaluate the body (1-iterator form). */
static bool param3d_eval(double t, Param3DEvalCtx* ctx,
                          double* xo, double* yo, double* zo) {
    symtab_add_own_value(ctx->var1->data.symbol.name, ctx->var1, expr_new_real(t));
    return eval_body_xyz(ctx, xo, yo, zo);
}

/* ---- 3D adaptive sampler (1-iterator form) ---- */

typedef struct {
    double t, x, y, z;
    bool   valid;
    bool   break_before;
} ParamPt3D;

#define PARAM3D_FLAT_TOL 0.0025

typedef struct {
    ParamPt3D* pts;
    size_t     len;
    size_t     cap;
    long       max_points;
    bool       pending_break;
} ParamBuf3D;

static void pbuf3d_push(ParamBuf3D* buf, double t, double x, double y, double z) {
    if (buf->len == buf->cap) {
        buf->cap *= 2;
        buf->pts = realloc(buf->pts, sizeof(ParamPt3D) * buf->cap);
    }
    ParamPt3D* p = &buf->pts[buf->len++];
    p->t = t; p->x = x; p->y = y; p->z = z;
    p->valid = true;
    p->break_before = buf->pending_break;
    buf->pending_break = false;
}

static bool pbuf3d_capped(const ParamBuf3D* buf) {
    return buf->max_points > 0 && (long)buf->len >= buf->max_points;
}

static void param3d_subdivide(Param3DEvalCtx* ctx, ParamBuf3D* buf,
                               double at, double ax, double ay, double az, bool a_ok,
                               double bt, double bx, double by, double bz, bool b_ok,
                               int depth, int max_rec, double diag) {
    if (pbuf3d_capped(buf) || depth >= max_rec) {
        if (b_ok) pbuf3d_push(buf, bt, bx, by, bz);
        else buf->pending_break = true;
        return;
    }

    double mt = (at + bt) / 2.0, mx = 0.0, my = 0.0, mz = 0.0;
    bool m_ok = param3d_eval(mt, ctx, &mx, &my, &mz);

    if (a_ok && b_ok && m_ok) {
        double lx = (ax + bx) / 2.0, ly = (ay + by) / 2.0, lz = (az + bz) / 2.0;
        double dev = sqrt((mx-lx)*(mx-lx) + (my-ly)*(my-ly) + (mz-lz)*(mz-lz));

        /* Three-probe anti-aliasing: prevents resonance with periodic curves. */
        double q1t = (at + mt) / 2.0, q3t = (mt + bt) / 2.0;
        double q1x = 0.0, q1y = 0.0, q1z = 0.0;
        double q3x = 0.0, q3y = 0.0, q3z = 0.0;
        bool q1_ok = param3d_eval(q1t, ctx, &q1x, &q1y, &q1z);
        bool q3_ok = param3d_eval(q3t, ctx, &q3x, &q3y, &q3z);
        if (q1_ok && q3_ok) {
            double ix1 = ax + 0.25*(bx-ax), iy1 = ay + 0.25*(by-ay), iz1 = az + 0.25*(bz-az);
            double ix3 = ax + 0.75*(bx-ax), iy3 = ay + 0.75*(by-ay), iz3 = az + 0.75*(bz-az);
            double d1 = sqrt((q1x-ix1)*(q1x-ix1) + (q1y-iy1)*(q1y-iy1) + (q1z-iz1)*(q1z-iz1));
            double d3 = sqrt((q3x-ix3)*(q3x-ix3) + (q3y-iy3)*(q3y-iy3) + (q3z-iz3)*(q3z-iz3));
            if (d1 > dev) dev = d1;
            if (d3 > dev) dev = d3;
        }
        if (dev < PARAM3D_FLAT_TOL * diag) {
            pbuf3d_push(buf, bt, bx, by, bz);
            return;
        }
    }

    param3d_subdivide(ctx, buf,
                       at, ax, ay, az, a_ok,
                       mt, mx, my, mz, m_ok,
                       depth+1, max_rec, diag);
    param3d_subdivide(ctx, buf,
                       mt, mx, my, mz, m_ok,
                       bt, bx, by, bz, b_ok,
                       depth+1, max_rec, diag);
}

/* Sample the parametric curve over [tmin, tmax].
 * Returns a heap-allocated ParamPt3D array (caller frees), or NULL. */
static ParamPt3D* param3d_sample(Param3DEvalCtx* ctx,
                                   double tmin, double tmax,
                                   long plot_points, int max_rec, long max_pts,
                                   size_t* out_count) {
    *out_count = 0;
    if (!ctx || !(tmin < tmax) || plot_points < 2) return NULL;

    double* gt  = malloc(sizeof(double) * (size_t)plot_points);
    double* gx  = malloc(sizeof(double) * (size_t)plot_points);
    double* gy  = malloc(sizeof(double) * (size_t)plot_points);
    double* gz  = malloc(sizeof(double) * (size_t)plot_points);
    bool*   gok = malloc(sizeof(bool)   * (size_t)plot_points);

    double xlo = 1e300, xhi = -1e300;
    double ylo = 1e300, yhi = -1e300;
    double zlo = 1e300, zhi = -1e300;
    for (long i = 0; i < plot_points; i++) {
        double t = tmin + (tmax - tmin) * (double)i / (double)(plot_points - 1);
        if (i == plot_points - 1) t = tmax;
        gt[i] = t;
        gok[i] = param3d_eval(t, ctx, &gx[i], &gy[i], &gz[i]);
        if (gok[i]) {
            if (gx[i] < xlo) xlo = gx[i];
            if (gx[i] > xhi) xhi = gx[i];
            if (gy[i] < ylo) ylo = gy[i];
            if (gy[i] > yhi) yhi = gy[i];
            if (gz[i] < zlo) zlo = gz[i];
            if (gz[i] > zhi) zhi = gz[i];
        }
    }
    double xspan = (xhi > xlo) ? (xhi - xlo) : 1.0;
    double yspan = (yhi > ylo) ? (yhi - ylo) : 1.0;
    double zspan = (zhi > zlo) ? (zhi - zlo) : 1.0;
    double diag  = sqrt(xspan*xspan + yspan*yspan + zspan*zspan);
    if (diag < 1e-12) diag = 1.0;

    ParamBuf3D buf;
    buf.cap          = (size_t)plot_points * 4;
    buf.pts          = malloc(sizeof(ParamPt3D) * buf.cap);
    buf.len          = 0;
    buf.max_points   = max_pts;
    buf.pending_break = false;

    if (gok[0]) pbuf3d_push(&buf, gt[0], gx[0], gy[0], gz[0]);
    else buf.pending_break = true;

    for (long i = 0; i + 1 < plot_points && !pbuf3d_capped(&buf); i++) {
        param3d_subdivide(ctx, &buf,
                           gt[i],   gx[i],   gy[i],   gz[i],   gok[i],
                           gt[i+1], gx[i+1], gy[i+1], gz[i+1], gok[i+1],
                           0, max_rec, diag);
    }

    free(gt); free(gx); free(gy); free(gz); free(gok);

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
    Expr* region_function;   /* borrowed */
    Expr* color_function;    /* borrowed */
    bool  color_function_scaling;
} Param3DSampleOpts;

/* Parses trailing Rule args starting at opts_start.
 * Outputs the evaluated pass-through option list (for Graphics3D) and the
 * resolved single PlotStyle color. AspectRatio is silently dropped (orbit
 * camera has no ratio). */
static bool split_options_param3d(Expr* res, size_t opts_start, long default_plot_points,
                                    Param3DSampleOpts* sopts,
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

    bool have_axes = false, have_style = false;

#define P3D_FAIL() do { \
    free(pass); \
    expr_free(*single_color_out); \
    *single_color_out = NULL; \
    return false; \
} while (0)

    for (size_t i = opts_start; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) P3D_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol.name : NULL;

        if (name == SYM_PlotPoints) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 2) P3D_FAIL();
            sopts->plot_points = v;
        } else if (name == SYM_MaxRecursion) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 0) P3D_FAIL();
            sopts->max_recursion = (int)v;
        } else if (name == SYM_MaxPlotPoints) {
            Expr* v = evaluate(expr_copy(rhs));
            bool is_inf = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Infinity);
            long lv = -1;
            bool ok = is_inf;
            if (!ok && v->type == EXPR_INTEGER && v->data.integer > 0) { lv = v->data.integer; ok = true; }
            expr_free(v);
            if (!ok) P3D_FAIL();
            sopts->max_plot_points = lv;
        } else if (name == SYM_Mesh) {
            Expr* v = evaluate(expr_copy(rhs));
            bool on  = (v->type == EXPR_SYMBOL
                        && (v->data.symbol.name == SYM_All || v->data.symbol.name == SYM_True));
            bool off = (v->type == EXPR_SYMBOL
                        && (v->data.symbol.name == SYM_None || v->data.symbol.name == SYM_False));
            expr_free(v);
            if (!on && !off) P3D_FAIL();
            sopts->mesh = on;
        } else if (name == SYM_RegionFunction) {
            sopts->region_function = rhs;
        } else if (name == SYM_ColorFunction) {
            sopts->color_function = rhs;
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = evaluate(expr_copy(rhs));
            sopts->color_function_scaling =
                !(v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_False);
            expr_free(v);
        } else if (name == SYM_PlotStyle) {
            have_style = true;
            if (*single_color_out) expr_free(*single_color_out);
            *single_color_out = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), expr_copy(*single_color_out) };
            pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        } else {
            /* AspectRatio: silently drop — no meaning for the orbit camera. */
            if (name == SYM_AspectRatio) continue;
            if (name == SYM_Axes) have_axes = true;
            Expr* val = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    if (!have_axes) {
        Expr* a[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_True) };
        pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_style) {
        Expr* rgb_a[3] = { expr_new_real(0.4), expr_new_real(0.7), expr_new_real(1.0) };
        Expr* rgb = expr_new_function(expr_new_symbol(SYM_RGBColor), rgb_a, 3);
        Expr* a[2] = { expr_new_symbol(SYM_PlotStyle), expr_copy(rgb) };
        pass[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        *single_color_out = rgb;
    }

    *pass_out       = pass;
    *pass_count_out = n;
    return true;
#undef P3D_FAIL
}

/* ---- Primitive helpers ---- */

static Expr* xyz_point(double x, double y, double z) {
    Expr* a[3] = { expr_new_real(x), expr_new_real(y), expr_new_real(z) };
    return expr_new_function(expr_new_symbol(SYM_List), a, 3);
}

static Expr* line3d_two_pts(double ax, double ay, double az,
                             double bx, double by, double bz) {
    Expr* pts[2] = { xyz_point(ax, ay, az), xyz_point(bx, by, bz) };
    Expr* plist  = expr_new_function(expr_new_symbol(SYM_List), pts, 2);
    Expr* largs[1] = { plist };
    return expr_new_function(expr_new_symbol(SYM_Line), largs, 1);
}

/* ---- Primitive assembly: 1-iterator (Line) form ---- */

static Expr** build_param3d_curve(Expr* body, Expr* var1,
                                    double tmin, double tmax,
                                    const Param3DSampleOpts* sopts,
                                    Expr* curve_color,
                                    size_t* out_count) {
    *out_count = 0;

    Param3DEvalCtx ctx = { .var1 = var1, .var2 = NULL, .body = body,
                             .region_function = sopts->region_function };

    size_t npts;
    ParamPt3D* pts = param3d_sample(&ctx, tmin, tmax,
                                      sopts->plot_points, sopts->max_recursion,
                                      sopts->max_plot_points, &npts);
    if (!pts || npts == 0) { free(pts); return NULL; }

    /* Compute xyz range for ColorFunction spatial scaling. */
    double xlo = pts[0].x, xhi = pts[0].x;
    double ylo = pts[0].y, yhi = pts[0].y;
    double zlo = pts[0].z, zhi = pts[0].z;
    for (size_t j = 1; j < npts; j++) {
        if (pts[j].x < xlo) xlo = pts[j].x;
        if (pts[j].x > xhi) xhi = pts[j].x;
        if (pts[j].y < ylo) ylo = pts[j].y;
        if (pts[j].y > yhi) yhi = pts[j].y;
        if (pts[j].z < zlo) zlo = pts[j].z;
        if (pts[j].z > zhi) zhi = pts[j].z;
    }

    /* Generous capacity: one color + one Line per segment at most, plus mesh. */
    Expr** prims = malloc(sizeof(Expr*) * (npts * 4 + 4));
    size_t pc = 0;
    size_t run_start = 0;

    for (size_t i = 1; i <= npts; i++) {
        bool end_of_run = (i == npts) || pts[i].break_before;
        if (!end_of_run) continue;

        size_t run_len = i - run_start;
        if (run_len >= 2) {
            if (sopts->color_function) {
                /* Per-segment color using spatial midpoint coordinates. */
                for (size_t j = 0; j + 1 < run_len; j++) {
                    size_t a = run_start + j, b = a + 1;
                    double cx = (pts[a].x + pts[b].x) / 2.0;
                    double cy = (pts[a].y + pts[b].y) / 2.0;
                    double cz = (pts[a].z + pts[b].z) / 2.0;
                    prims[pc++] = eval_color_function3(sopts->color_function,
                                                        cx, cy, cz,
                                                        xlo, xhi, ylo, yhi, zlo, zhi,
                                                        sopts->color_function_scaling);
                    prims[pc++] = line3d_two_pts(pts[a].x, pts[a].y, pts[a].z,
                                                  pts[b].x, pts[b].y, pts[b].z);
                }
            } else {
                Expr** line_pts = malloc(sizeof(Expr*) * run_len);
                for (size_t j = 0; j < run_len; j++) {
                    line_pts[j] = xyz_point(pts[run_start + j].x,
                                             pts[run_start + j].y,
                                             pts[run_start + j].z);
                }
                Expr* pts_list = expr_new_function(expr_new_symbol(SYM_List),
                                                    line_pts, run_len);
                free(line_pts);
                Expr* largs[1] = { pts_list };
                prims[pc++] = expr_new_function(expr_new_symbol(SYM_Line), largs, 1);
            }
        }
        run_start = i;
    }

    if (sopts->mesh && npts > 0) {
        if (curve_color && !sopts->color_function) prims[pc++] = expr_copy(curve_color);
        double diag = sqrt((xhi-xlo)*(xhi-xlo) + (yhi-ylo)*(yhi-ylo) + (zhi-zlo)*(zhi-zlo));
        double msize = (diag > 0) ? diag * 0.0025 : 0.005;
        Expr* ps_arg[1] = { expr_new_real(msize) };
        prims[pc++] = expr_new_function(expr_new_symbol(SYM_PointSize), ps_arg, 1);
        Expr** dots = malloc(sizeof(Expr*) * npts);
        for (size_t j = 0; j < npts; j++)
            dots[j] = xyz_point(pts[j].x, pts[j].y, pts[j].z);
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

/* ---- Primitive assembly: 2-iterator (Polygon surface) form ---- */

typedef struct { double x, y, z; bool valid; } GridPt3DParam;

static Expr** build_param3d_surface(Expr* body, Expr* var1, Expr* var2,
                                     double t1min, double t1max,
                                     double t2min, double t2max,
                                     const Param3DSampleOpts* sopts,
                                     Expr* color,
                                     size_t* out_count) {
    *out_count = 0;
    long n = sopts->plot_points;
    if (n < 2) return NULL;

    size_t ntotal = (size_t)(n * n);
    GridPt3DParam* grid = malloc(sizeof(GridPt3DParam) * ntotal);

    Param3DEvalCtx ctx = { .var1 = var1, .var2 = var2, .body = body,
                             .region_function = sopts->region_function };

    /* Sample the grid and compute xyz bbox for ColorFunction scaling. */
    double xlo = 1e300, xhi = -1e300;
    double ylo = 1e300, yhi = -1e300;
    double zlo = 1e300, zhi = -1e300;

    for (long i = 0; i < n; i++) {
        double t = (n > 1) ? t1min + (t1max - t1min) * (double)i / (double)(n-1) : t1min;
        if (i == n-1) t = t1max;
        symtab_add_own_value(ctx.var1->data.symbol.name, ctx.var1, expr_new_real(t));
        for (long j = 0; j < n; j++) {
            double u = (n > 1) ? t2min + (t2max - t2min) * (double)j / (double)(n-1) : t2min;
            if (j == n-1) u = t2max;
            symtab_add_own_value(ctx.var2->data.symbol.name, ctx.var2, expr_new_real(u));
            GridPt3DParam* p = &grid[i * n + j];
            p->valid = eval_body_xyz(&ctx, &p->x, &p->y, &p->z);
            if (p->valid) {
                if (p->x < xlo) xlo = p->x;
                if (p->x > xhi) xhi = p->x;
                if (p->y < ylo) ylo = p->y;
                if (p->y > yhi) yhi = p->y;
                if (p->z < zlo) zlo = p->z;
                if (p->z > zhi) zhi = p->z;
            }
        }
    }

    /* Worst-case capacity: 2 per cell (color+polygon) + mesh lines. */
    long ncells = (n-1) * (n-1);
    long nmesh  = sopts->mesh ? 2 * (n-1) * n : 0;
    Expr** prims = malloc(sizeof(Expr*) * (size_t)(2 * ncells + nmesh + 4));
    size_t pc = 0;

    /* Color directive up front (solid by default; pass PlotStyle -> {color, Opacity[a]}
     * to get a transparent surface). */
    if (color && !sopts->color_function)
        prims[pc++] = expr_copy(color);

    /* Polygon quads. */
    for (long i = 0; i < n-1; i++) {
        for (long j = 0; j < n-1; j++) {
            GridPt3DParam* p00 = &grid[i * n + j];
            GridPt3DParam* p10 = &grid[(i+1) * n + j];
            GridPt3DParam* p11 = &grid[(i+1) * n + (j+1)];
            GridPt3DParam* p01 = &grid[i * n + (j+1)];

            if (!p00->valid || !p10->valid || !p11->valid || !p01->valid) continue;

            if (sopts->color_function) {
                double cx = (p00->x + p10->x + p11->x + p01->x) / 4.0;
                double cy = (p00->y + p10->y + p11->y + p01->y) / 4.0;
                double cz = (p00->z + p10->z + p11->z + p01->z) / 4.0;
                prims[pc++] = eval_color_function3(sopts->color_function,
                                                    cx, cy, cz,
                                                    xlo, xhi, ylo, yhi, zlo, zhi,
                                                    sopts->color_function_scaling);
            }

            Expr* verts[4] = {
                xyz_point(p00->x, p00->y, p00->z),
                xyz_point(p10->x, p10->y, p10->z),
                xyz_point(p11->x, p11->y, p11->z),
                xyz_point(p01->x, p01->y, p01->z),
            };
            Expr* vlist  = expr_new_function(expr_new_symbol(SYM_List), verts, 4);
            Expr* pargs[1] = { vlist };
            prims[pc++] = expr_new_function(expr_new_symbol(SYM_Polygon), pargs, 1);
        }
    }

    /* Optional mesh: interior grid lines (all edges since we have no
     * neighbour-validity filter; the surface is parametric, not over a
     * Cartesian domain, so all cells are equally "interior"). */
    if (sopts->mesh) {
        for (long i = 0; i < n-1; i++) {
            for (long j = 0; j < n; j++) {
                GridPt3DParam* pa = &grid[i * n + j];
                GridPt3DParam* pb = &grid[(i+1) * n + j];
                if (pa->valid && pb->valid)
                    prims[pc++] = line3d_two_pts(pa->x, pa->y, pa->z,
                                                  pb->x, pb->y, pb->z);
            }
        }
        for (long i = 0; i < n; i++) {
            for (long j = 0; j < n-1; j++) {
                GridPt3DParam* pa = &grid[i * n + j];
                GridPt3DParam* pb = &grid[i * n + (j+1)];
                if (pa->valid && pb->valid)
                    prims[pc++] = line3d_two_pts(pa->x, pa->y, pa->z,
                                                  pb->x, pb->y, pb->z);
            }
        }
    }

    free(grid);
    if (pc == 0) { free(prims); return NULL; }
    *out_count = pc;
    return prims;
}

/* ---- Helper: is `e` a valid {var, min, max} iterator spec? ---- */
static bool is_iterator_3d(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List
        && e->data.function.arg_count == 3
        && e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* ---- Pre-scan for PlotLegends before option parsing ---- */
static Expr* prescan_plotlegends3d(Expr* res, size_t opts_start) {
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

/* ---- Assemble Graphics3D[prim_list, opts..., legend_meta?] ---- */
static Expr* make_graphics3d(Expr* prim_list, Expr** pass, size_t npass,
                               Expr* legend_meta) {
    size_t g_argc = 1 + npass + (legend_meta ? 1 : 0);
    Expr** g_args = malloc(sizeof(Expr*) * g_argc);
    g_args[0] = prim_list;
    for (size_t i = 0; i < npass; i++) g_args[1 + i] = pass[i];
    if (legend_meta) g_args[g_argc - 1] = legend_meta;
    Expr* g = expr_new_function(expr_new_symbol(SYM_Graphics3D), g_args, g_argc);
    free(g_args);
    return g;
}

/* ---- Entry point ---- */

Expr* builtin_parametricplot3d(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* First positional arg after the body must be a {var, min, max} iterator. */
    Expr* iter1 = res->data.function.args[1];
    if (!is_iterator_3d(iter1)) return NULL;

    Expr* var1 = iter1->data.function.args[0];
    double t1min, t1max;
    if (!numericize_bound(iter1->data.function.args[1], &t1min)
        || !numericize_bound(iter1->data.function.args[2], &t1max)
        || !(t1min < t1max)) return NULL;

    /* Check for optional second iterator → 2-iterator (surface) form. */
    bool two_iter = (argc >= 3) && is_iterator_3d(res->data.function.args[2]);
    size_t opts_start = two_iter ? 3 : 2;

    Expr* legends = prescan_plotlegends3d(res, opts_start);

    if (two_iter) {
        /* ---- 2-iterator form: parametric 3D surface patch ---- */
        Expr* iter2 = res->data.function.args[2];
        Expr* var2  = iter2->data.function.args[0];
        double t2min, t2max;
        if (!numericize_bound(iter2->data.function.args[1], &t2min)
            || !numericize_bound(iter2->data.function.args[2], &t2max)
            || !(t2min < t2max)) { expr_free(legends); return NULL; }

        Expr* body = res->data.function.args[0];

        /* Multi-surface detection: body is a List whose first element is also
         * a List (i.e. {{fx1,fy1,fz1}, {fx2,fy2,fz2}}), not a single
         * {fx,fy,fz} triple. */
        bool body_is_list = (body->type == EXPR_FUNCTION
                              && body->data.function.head->type == EXPR_SYMBOL
                              && body->data.function.head->data.symbol.name == SYM_List);
        bool multi2 = (body_is_list
                       && body->data.function.arg_count >= 1
                       && body->data.function.args[0]->type == EXPR_FUNCTION
                       && body->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                       && body->data.function.args[0]->data.function.head->data.symbol.name
                          == SYM_List);

        Expr** sub_bodies;
        size_t nsub;
        if (multi2) {
            sub_bodies = body->data.function.args;
            nsub       = body->data.function.arg_count;
        } else {
            sub_bodies = &body;
            nsub       = 1;
        }

        Param3DSampleOpts sopts;
        Expr** pass;
        size_t npass;
        Expr* color = NULL;
        if (!split_options_param3d(res, 3, 25, &sopts, &pass, &npass, &color)) {
            expr_free(legends); return NULL;
        }

        Rule* old1 = iter_spec_shadow(var1);
        Rule* old2 = iter_spec_shadow(var2);

        size_t total = 0;
        Expr*** per_prims  = malloc(sizeof(Expr**) * nsub);
        size_t* per_counts = malloc(sizeof(size_t)  * nsub);
        bool any = false;

        for (size_t si = 0; si < nsub; si++) {
            Expr* sc = multi2 ? palette_color(si) : color;
            per_prims[si] = build_param3d_surface(sub_bodies[si], var1, var2,
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

        Expr* g = make_graphics3d(prim_list, pass, npass, legend_meta);
        free(pass);
        expr_free(color);
        return g;

    } else {
        /* ---- 1-iterator form: 3D parametric space curve(s) ---- */
        Expr* spec = res->data.function.args[0];

        /* Multi-curve: spec is a List whose first element is also a List. */
        bool is_list_head = (spec->type == EXPR_FUNCTION
                              && spec->data.function.head->type == EXPR_SYMBOL
                              && spec->data.function.head->data.symbol.name == SYM_List);
        bool multi = (is_list_head
                      && spec->data.function.arg_count >= 1
                      && spec->data.function.args[0]->type == EXPR_FUNCTION
                      && spec->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                      && spec->data.function.args[0]->data.function.head->data.symbol.name
                         == SYM_List);

        Expr** bodies;
        size_t nbodies;
        if (multi) {
            bodies  = spec->data.function.args;
            nbodies = spec->data.function.arg_count;
        } else {
            bodies  = &spec;
            nbodies = 1;
        }

        Param3DSampleOpts sopts;
        Expr** pass;
        size_t npass;
        Expr* single_color = NULL;
        if (!split_options_param3d(res, 2, 25, &sopts, &pass, &npass, &single_color)) {
            expr_free(legends); return NULL;
        }

        Rule* old1 = iter_spec_shadow(var1);

        size_t total = 0;
        Expr*** per_prims  = malloc(sizeof(Expr**) * nbodies);
        size_t* per_counts = malloc(sizeof(size_t)  * nbodies);
        bool any = false;

        for (size_t ci = 0; ci < nbodies; ci++) {
            Expr* cc = multi ? palette_color(ci) : single_color;
            per_prims[ci] = build_param3d_curve(bodies[ci], var1, t1min, t1max,
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

        Expr* g = make_graphics3d(prim_list, pass, npass, legend_meta);
        free(pass);
        expr_free(single_color);
        return g;
    }
}
