/* complexplot.c — ComplexPlot and ComplexPlot3D.
 *
 * Both functions share the same domain-parsing logic and option set.
 * The only structural difference is what they build from the evaluated
 * grid: ComplexPlot emits Rectangle primitives into Graphics[], and
 * ComplexPlot3D emits Polygon quads (height = |w|, colour = arg(w))
 * into Graphics3D[].
 *
 * Coloring convention: the default phase-to-color mapping uses the same
 * thermal_rgb ramp that DensityPlot and Plot3D use, keyed to the
 * normalised argument: t = (atan2(im, re) + π) / (2π) ∈ [0, 1].
 * This keeps the palette consistent across the whole graphics engine.
 * A custom ColorFunction receives (re, im) as a 2-arg call; with
 * ColorFunctionScaling→True (default) the arguments are first scaled
 * so that re ∈ [0,1] and im ∈ [0,1] across the sampled domain.
 *
 * Both are HoldAll: the body and the iterator spec are held unevaluated
 * until z is bound to Complex[x, y] at each grid point.
 *
 * Domain spec:
 *   {z, zmin, zmax}  — z is the complex iterator variable; zmin and zmax
 *   are evaluated and their Re/Im parts define the rectangular plotting
 *   domain: xmin=Re(zmin), xmax=Re(zmax), ymin=Im(zmin), ymax=Im(zmax).
 *   Both endpoints may be real (imaginary part = 0). */

#include "complexplot.h"
#include "plot_common.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "arithmetic.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Option struct shared by both 2D and 3D                               */
/* ------------------------------------------------------------------ */

typedef struct {
    int   plot_points;
    Expr* color_function;         /* borrowed; NULL = thermal default */
    bool  color_function_scaling;
    Expr* region_function;        /* borrowed; NULL = none */
    bool  show_legend;            /* PlotLegends -> Automatic / True */
} CPlotOpts;

/* Parse trailing Rule args starting at index `first_opt_idx`.
 * Fills `o`, builds a passthrough list for Graphics/Graphics3D opts.
 * Returns true on success; on failure frees *pt_out and returns false. */
static bool split_cplot_options(Expr* res, size_t first_opt_idx,
                                 CPlotOpts* o,
                                 Expr*** pt_out, size_t* pt_n_out) {
    o->plot_points            = 200;
    o->color_function         = NULL;
    o->color_function_scaling = true;
    o->region_function        = NULL;
    o->show_legend            = false;

    size_t argc = res->data.function.arg_count;
    size_t extra = (argc > first_opt_idx ? argc - first_opt_idx : 0);
    /* +4 headroom: AspectRatio, Axes, PlotRange, $StreamColorBar defaults */
    Expr** pt = malloc(sizeof(Expr*) * (extra + 5));
    size_t n = 0;

    bool have_axes   = false;
    bool have_aspect = false;
    bool have_frame  = false;

#define CP_FAIL() do { free(pt); return false; } while (0)

    for (size_t i = first_opt_idx; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) CP_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol : NULL;

        if (name == SYM_PlotPoints) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 2) CP_FAIL();
            o->plot_points = (int)v;
        } else if (name == SYM_ColorFunction) {
            o->color_function = rhs; /* borrowed */
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = evaluate(expr_copy(rhs));
            o->color_function_scaling = !(v->type == EXPR_SYMBOL
                                           && v->data.symbol == SYM_False);
            expr_free(v);
        } else if (name == SYM_RegionFunction) {
            o->region_function = rhs; /* borrowed */
        } else if (name == SYM_PlotLegends) {
            Expr* v = evaluate(expr_copy(rhs));
            o->show_legend = !(v->type == EXPR_SYMBOL
                               && (v->data.symbol == SYM_None
                                   || v->data.symbol == SYM_False));
            expr_free(v);
        } else {
            if      (name == SYM_Axes)        have_axes   = true;
            else if (name == SYM_AspectRatio) have_aspect = true;
            else if (name == SYM_Frame) {
                if (!(rhs->type == EXPR_SYMBOL
                      && (rhs->data.symbol == SYM_False
                          || rhs->data.symbol == SYM_None)))
                    have_frame = true;
            }
            Expr* val  = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    /* Inject defaults that match the existing field-plot conventions. */
    if (!have_axes && !have_frame) {
        Expr* a[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_True) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_aspect) {
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_integer(1) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }

    *pt_out   = pt;
    *pt_n_out = n;
    return true;
#undef CP_FAIL
}

/* ------------------------------------------------------------------ */
/* Domain parsing: {z, zmin, zmax}                                      */
/* ------------------------------------------------------------------ */

/* Extract (re, im) from an already-evaluated expression.
 * Accepts Complex[a,b], a pure real, or an integer; returns false otherwise. */
static bool bound_to_complex(Expr* e, double* re_out, double* im_out) {
    double r;
    if (expr_to_real_double(e, &r) && isfinite(r)) {
        *re_out = r;
        *im_out = 0.0;
        return true;
    }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double rv, iv;
        if (expr_to_real_double(re, &rv) && isfinite(rv)
            && expr_to_real_double(im, &iv) && isfinite(iv)) {
            *re_out = rv;
            *im_out = iv;
            return true;
        }
    }
    return false;
}

/* Parse iterator arg `{z, zmin, zmax}` and extract:
 *   zvar   — borrowed symbol expression from the iterator
 *   xmin/xmax — Re(zmin), Re(zmax)
 *   ymin/ymax — Im(zmin), Im(zmax)
 * Returns false on any parse / numeric failure. */
static bool parse_complex_iterator(Expr* iter, Expr** zvar_out,
                                    double* xmin, double* xmax,
                                    double* ymin, double* ymax) {
    if (!iter || iter->type != EXPR_FUNCTION
        || !iter->data.function.head
        || iter->data.function.head->type != EXPR_SYMBOL
        || iter->data.function.head->data.symbol != SYM_List
        || iter->data.function.arg_count != 3)
        return false;

    Expr* zvar = iter->data.function.args[0];
    if (!zvar || zvar->type != EXPR_SYMBOL) return false;

    Expr* e_zmin = evaluate(expr_copy(iter->data.function.args[1]));
    Expr* e_zmax = evaluate(expr_copy(iter->data.function.args[2]));

    double x0, y0, x1, y1;
    bool ok = bound_to_complex(e_zmin, &x0, &y0)
           && bound_to_complex(e_zmax, &x1, &y1);
    expr_free(e_zmin);
    expr_free(e_zmax);

    if (!ok) return false;
    if (x0 == x1 || y0 == y1) return false; /* degenerate domain */

    *zvar_out = zvar;   /* borrowed from res — do not free */
    *xmin = (x0 < x1 ? x0 : x1);
    *xmax = (x0 < x1 ? x1 : x0);
    *ymin = (y0 < y1 ? y0 : y1);
    *ymax = (y0 < y1 ? y1 : y0);
    return true;
}

/* ------------------------------------------------------------------ */
/* Grid evaluation                                                      */
/* ------------------------------------------------------------------ */

/* Bind z = Complex[x, y], evaluate body, extract (re, im) of result.
 * Returns true and sets *re_out, *im_out on success; returns false on
 * failure (unevaluated, non-numeric, or infinite result). */
static bool cp_eval(Expr* zvar, Expr* body, double x, double y,
                    double* re_out, double* im_out) {
    Expr* ra[2] = { expr_new_real(x), expr_new_real(y) };
    Expr* zval  = expr_new_function(expr_new_symbol(SYM_Complex), ra, 2);
    symtab_add_own_value(zvar->data.symbol, zvar, zval);

    Expr* result = evaluate(expr_copy(body));
    bool ok = false;
    double r = 0.0, i = 0.0;

    double rv;
    Expr *rp, *ip;
    if (expr_to_real_double(result, &rv) && isfinite(rv)) {
        r = rv; i = 0.0; ok = true;
    } else if (is_complex(result, &rp, &ip)) {
        double rv2, iv;
        if (expr_to_real_double(rp, &rv2) && isfinite(rv2)
            && expr_to_real_double(ip, &iv) && isfinite(iv)) {
            r = rv2; i = iv; ok = true;
        }
    }
    expr_free(result);
    expr_free(zval);

    if (ok) { *re_out = r; *im_out = i; }
    return ok;
}

/* ------------------------------------------------------------------ */
/* Coloring                                                             */
/* ------------------------------------------------------------------ */

/* Default: thermal_rgb keyed to normalised argument.
 * t = (atan2(im, re) + π) / (2π) — wraps continuously in [0, 1].
 * Multiplied by brightness b = |w|/(1+|w|) to fade near the origin. */
static Expr* cp_default_color(double re, double im) {
    double arg = atan2(im, re);
    double t   = (arg + M_PI) / (2.0 * M_PI);
    /* Clamp floating-point edge cases (atan2 range is (-π, π]). */
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    /* Modulus-based brightness: maps [0,∞) → [0,1) via x/(1+x). */
    double mod   = sqrt(re * re + im * im);
    double bright = mod / (1.0 + mod);
    double r, g, b;
    thermal_rgb(t, &r, &g, &b);
    /* Scale by brightness so the origin fades to black. */
    r *= bright; g *= bright; b *= bright;
    Expr* a[3] = { expr_new_real(r), expr_new_real(g), expr_new_real(b) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
}

static bool is_color_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    return h == SYM_RGBColor || h == SYM_GrayLevel || h == SYM_Hue || h == SYM_CMYKColor;
}

/* Resolve color for one grid cell from a custom ColorFunction or the default.
 * re/im: raw result values; re_sc/im_sc: values scaled to [0,1] for
 * ColorFunctionScaling→True. */
static Expr* cp_color(Expr* cfn, bool scaling,
                       double re, double im, double re_sc, double im_sc) {
    if (!cfn) return cp_default_color(re, im);

    if (cfn->type == EXPR_STRING) {
        /* "PhaseRings" needs both re and im — intercept before the 1-D ramp path. */
        if (strcmp(cfn->data.string, "PhaseRings") == 0) {
            double rv, gv, bv;
            phase_rings_rgb(re, im, &rv, &gv, &bv);
            Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
            return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
        }
        /* Other named ramps: key to normalised argument ∈ [0, 1]. */
        double arg = atan2(im, re);
        double t   = (arg + M_PI) / (2.0 * M_PI);
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        Expr* c = named_color_ramp(cfn->data.string, t);
        if (c) return c;
    }

    /* Custom function: try f[re, im] (2-arg) then f[re] (1-arg). */
    double u = scaling ? re_sc : re;
    double v = scaling ? im_sc : im;

    Expr* a2[2] = { expr_new_real(u), expr_new_real(v) };
    Expr* call2 = expr_new_function(expr_copy(cfn), a2, 2);
    Expr* r2    = evaluate(call2);
    expr_free(call2);
    if (is_color_head(r2)) return r2;
    expr_free(r2);

    Expr* a1[1] = { expr_new_real(u) };
    Expr* call1 = expr_new_function(expr_copy(cfn), a1, 1);
    Expr* r1    = evaluate(call1);
    expr_free(call1);
    if (is_color_head(r1)) return r1;
    expr_free(r1);

    /* Fallback: neutral mid-gray */
    Expr* ga[1] = { expr_new_real(0.5) };
    return expr_new_function(expr_new_symbol(SYM_GrayLevel), ga, 1);
}

/* ------------------------------------------------------------------ */
/* Grid storage                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    double re, im;  /* evaluated complex result */
    bool   valid;   /* evaluation succeeded and RegionFunction passed */
} CGrid;

/* Height-cap for ComplexPlot3D: sort valid |f(z)| values and return the
 * 95th-percentile so that pole spikes are clipped to a sensible ceiling
 * (matching Mathematica's automatic PlotRange clipping behaviour).
 * A cap of exactly 0 is replaced by 1 to avoid a degenerate surface. */
static int cmp_double(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

static double compute_height_cap(const CGrid* grid, size_t total) {
    /* Collect valid moduli */
    double* hs = malloc(sizeof(double) * total);
    if (!hs) return 1.0;
    size_t nc = 0;
    for (size_t k = 0; k < total; k++) {
        if (!grid[k].valid) continue;
        double h = sqrt(grid[k].re * grid[k].re + grid[k].im * grid[k].im);
        if (isfinite(h)) hs[nc++] = h;
    }
    if (nc == 0) { free(hs); return 1.0; }
    qsort(hs, nc, sizeof(double), cmp_double);
    /* 95th percentile clips the top 5% of heights (pole spikes) */
    size_t idx = (size_t)(nc * 0.95);
    if (idx >= nc) idx = nc - 1;
    double cap = hs[idx];
    free(hs);
    return (cap > 1e-10) ? cap : 1.0;
}

static CGrid* build_cgrid(Expr* zvar, Expr* body, Expr* region_fn,
                           double xmin, double xmax,
                           double ymin, double ymax, int N) {
    CGrid* grid = malloc(sizeof(CGrid) * (size_t)(N + 1) * (size_t)(N + 1));
    if (!grid) return NULL;

    double dx = (xmax - xmin) / N;
    double dy = (ymax - ymin) / N;

    for (int iy = 0; iy <= N; iy++) {
        double y = ymin + iy * dy;
        if (iy == N) y = ymax;
        for (int ix = 0; ix <= N; ix++) {
            double x = xmin + ix * dx;
            if (ix == N) x = xmax;
            CGrid* p = &grid[iy * (N + 1) + ix];
            double re, im;
            bool ok = cp_eval(zvar, body, x, y, &re, &im);
            if (ok && region_fn && !eval_region(region_fn, x, y)) ok = false;
            p->re    = ok ? re : 0.0;
            p->im    = ok ? im : 0.0;
            p->valid = ok;
        }
    }
    return grid;
}

/* Compute the range of re and im values across valid grid cells.
 * Used for ColorFunctionScaling normalization. */
static void grid_rerange(const CGrid* grid, size_t total,
                          double* re_min, double* re_max,
                          double* im_min, double* im_max) {
    *re_min =  1e300; *re_max = -1e300;
    *im_min =  1e300; *im_max = -1e300;
    for (size_t k = 0; k < total; k++) {
        if (!grid[k].valid) continue;
        if (grid[k].re < *re_min) *re_min = grid[k].re;
        if (grid[k].re > *re_max) *re_max = grid[k].re;
        if (grid[k].im < *im_min) *im_min = grid[k].im;
        if (grid[k].im > *im_max) *im_max = grid[k].im;
    }
    if (*re_min > *re_max) { *re_min = 0.0; *re_max = 1.0; }
    if (*im_min > *im_max) { *im_min = 0.0; *im_max = 1.0; }
}

/* ------------------------------------------------------------------ */
/* Embed PlotRange matching the domain                                  */
/* ------------------------------------------------------------------ */

static void embed_plot_range(double xmin, double xmax, double ymin, double ymax,
                              Expr*** pt, size_t* pt_n) {
    /* Check if caller already supplied a PlotRange rule */
    for (size_t i = 0; i < *pt_n; i++) {
        const Expr* e = (*pt)[i];
        if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
            && e->data.function.args[0]->type == EXPR_SYMBOL
            && e->data.function.args[0]->data.symbol == SYM_PlotRange)
            return;
    }
    *pt = realloc(*pt, sizeof(Expr*) * (*pt_n + 1));
    Expr* xr[2] = { expr_new_real(xmin), expr_new_real(xmax) };
    Expr* yr[2] = { expr_new_real(ymin), expr_new_real(ymax) };
    Expr* rl[2] = { expr_new_function(expr_new_symbol(SYM_List), xr, 2),
                    expr_new_function(expr_new_symbol(SYM_List), yr, 2) };
    Expr* pr    = expr_new_function(expr_new_symbol(SYM_List), rl, 2);
    Expr* a[2]  = { expr_new_symbol(SYM_PlotRange), pr };
    (*pt)[(*pt_n)++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

/* 3D variant: x/y domain + explicit z-range so the renderer frames the
 * surface correctly when pole heights are clamped. */
static void embed_plot_range3(double xmin, double xmax,
                               double ymin, double ymax,
                               double zmin, double zmax,
                               Expr*** pt, size_t* pt_n) {
    for (size_t i = 0; i < *pt_n; i++) {
        const Expr* e = (*pt)[i];
        if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
            && e->data.function.args[0]->type == EXPR_SYMBOL
            && e->data.function.args[0]->data.symbol == SYM_PlotRange)
            return;
    }
    *pt = realloc(*pt, sizeof(Expr*) * (*pt_n + 1));
    Expr* xr[2] = { expr_new_real(xmin), expr_new_real(xmax) };
    Expr* yr[2] = { expr_new_real(ymin), expr_new_real(ymax) };
    Expr* zr[2] = { expr_new_real(zmin), expr_new_real(zmax) };
    Expr* rl[3] = { expr_new_function(expr_new_symbol(SYM_List), xr, 2),
                    expr_new_function(expr_new_symbol(SYM_List), yr, 2),
                    expr_new_function(expr_new_symbol(SYM_List), zr, 2) };
    Expr* pr   = expr_new_function(expr_new_symbol(SYM_List), rl, 3);
    Expr* a[2] = { expr_new_symbol(SYM_PlotRange), pr };
    (*pt)[(*pt_n)++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

/* Append $StreamColorBar[-π, π, cfn_or_Automatic] so the renderer draws a
 * vertical phase-angle color scale.  The bar parameter t ∈ [0,1] maps to
 * arg = -π + t·2π, which is identical to our normalization
 * t = (arg+π)/(2π) — so the renderer's thermal_rgb(t) or named-ramp
 * call produces exactly the same colors as cp_default_color / cp_color. */
static void emit_phase_color_bar(Expr* cfn, Expr*** pt, size_t* pt_n) {
    *pt = realloc(*pt, sizeof(Expr*) * (*pt_n + 1));
    Expr* cfn_copy = cfn ? expr_copy(cfn) : expr_new_symbol(SYM_Automatic);
    Expr* cb_args[3] = { expr_new_real(-M_PI), expr_new_real(M_PI), cfn_copy };
    (*pt)[(*pt_n)++] = expr_new_function(expr_new_symbol(SYM_StreamColorBar), cb_args, 3);
}

/* ------------------------------------------------------------------ */
/* builtin_complexplot                                                  */
/* ------------------------------------------------------------------ */

Expr* builtin_complexplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* body = res->data.function.args[0]; /* held */

    /* Parse {z, zmin, zmax} */
    Expr*  zvar;
    double xmin, xmax, ymin, ymax;
    if (!parse_complex_iterator(res->data.function.args[1],
                                 &zvar, &xmin, &xmax, &ymin, &ymax))
        return NULL;

    CPlotOpts opts;
    Expr** pt = NULL; size_t pt_n = 0;
    if (!split_cplot_options(res, 2, &opts, &pt, &pt_n)) return NULL;

    int N = opts.plot_points;
    size_t stride = (size_t)(N + 1);

    Rule* old_z = iter_spec_shadow(zvar);
    CGrid* grid = build_cgrid(zvar, body, opts.region_function,
                               xmin, xmax, ymin, ymax, N);
    iter_spec_restore(zvar, old_z);

    if (!grid) {
        for (size_t i = 0; i < pt_n; i++) expr_free(pt[i]);
        free(pt);
        return NULL;
    }

    /* Compute re/im ranges for ColorFunctionScaling */
    double re_min, re_max, im_min, im_max;
    grid_rerange(grid, stride * stride, &re_min, &re_max, &im_min, &im_max);
    double re_span = (re_max > re_min) ? (re_max - re_min) : 1.0;
    double im_span = (im_max > im_min) ? (im_max - im_min) : 1.0;

    double dx = (xmax - xmin) / N;
    double dy = (ymax - ymin) / N;

    /* 2 primitives per cell: color directive + Rectangle */
    size_t cap   = (size_t)N * (size_t)N * 2 + 2;
    Expr** prims = malloc(sizeof(Expr*) * cap);
    size_t np    = 0;

    for (int iy = 0; iy < N; iy++) {
        double y0 = ymin + iy * dy;
        for (int ix = 0; ix < N; ix++) {
            double x0 = xmin + ix * dx;

            /* Average the four corner values for smooth appearance */
            const CGrid* p00 = &grid[ iy      * (int)stride + ix    ];
            const CGrid* p10 = &grid[ iy      * (int)stride + ix + 1];
            const CGrid* p11 = &grid[(iy + 1) * (int)stride + ix + 1];
            const CGrid* p01 = &grid[(iy + 1) * (int)stride + ix    ];

            if (!p00->valid || !p10->valid || !p11->valid || !p01->valid)
                continue;

            double re_avg = (p00->re + p10->re + p11->re + p01->re) * 0.25;
            double im_avg = (p00->im + p10->im + p11->im + p01->im) * 0.25;

            double re_sc = (re_avg - re_min) / re_span;
            double im_sc = (im_avg - im_min) / im_span;
            if (re_sc < 0.0) re_sc = 0.0; if (re_sc > 1.0) re_sc = 1.0;
            if (im_sc < 0.0) im_sc = 0.0; if (im_sc > 1.0) im_sc = 1.0;

            prims[np++] = cp_color(opts.color_function, opts.color_function_scaling,
                                    re_avg, im_avg, re_sc, im_sc);

            /* Rectangle in plot coordinates (x = Re axis, y = Im axis) */
            Expr* p1[2] = { expr_new_real(x0),      expr_new_real(y0) };
            Expr* p2[2] = { expr_new_real(x0 + dx), expr_new_real(y0 + dy) };
            Expr* ra[2] = { expr_new_function(expr_new_symbol(SYM_List), p1, 2),
                             expr_new_function(expr_new_symbol(SYM_List), p2, 2) };
            prims[np++] = expr_new_function(expr_new_symbol(SYM_Rectangle), ra, 2);
        }
    }

    free(grid);

    embed_plot_range(xmin, xmax, ymin, ymax, &pt, &pt_n);
    if (opts.show_legend) emit_phase_color_bar(opts.color_function, &pt, &pt_n);

    Expr* plist  = expr_new_function(expr_new_symbol(SYM_List), prims, np);
    free(prims);

    size_t gargc = 1 + pt_n;
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = plist;
    for (size_t i = 0; i < pt_n; i++) gargs[1 + i] = pt[i];
    free(pt);

    return expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
}

/* ------------------------------------------------------------------ */
/* builtin_complexplot3d                                                */
/* ------------------------------------------------------------------ */

static Expr* point3d(double x, double y, double z) {
    Expr* a[3] = { expr_new_real(x), expr_new_real(y), expr_new_real(z) };
    return expr_new_function(expr_new_symbol(SYM_List), a, 3);
}

Expr* builtin_complexplot3d(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* body = res->data.function.args[0]; /* held */

    Expr*  zvar;
    double xmin, xmax, ymin, ymax;
    if (!parse_complex_iterator(res->data.function.args[1],
                                 &zvar, &xmin, &xmax, &ymin, &ymax))
        return NULL;

    CPlotOpts opts;
    Expr** pt = NULL; size_t pt_n = 0;
    /* 3D: default passthrough gets Axes → True (reusing split_cplot_options) */
    if (!split_cplot_options(res, 2, &opts, &pt, &pt_n)) return NULL;

    int N = opts.plot_points;
    size_t stride = (size_t)(N + 1);

    Rule* old_z = iter_spec_shadow(zvar);
    CGrid* grid = build_cgrid(zvar, body, opts.region_function,
                               xmin, xmax, ymin, ymax, N);
    iter_spec_restore(zvar, old_z);

    if (!grid) {
        for (size_t i = 0; i < pt_n; i++) expr_free(pt[i]);
        free(pt);
        return NULL;
    }

    /* Re/im ranges for ColorFunctionScaling */
    double re_min, re_max, im_min, im_max;
    grid_rerange(grid, stride * stride, &re_min, &re_max, &im_min, &im_max);
    double re_span = (re_max > re_min) ? (re_max - re_min) : 1.0;
    double im_span = (im_max > im_min) ? (im_max - im_min) : 1.0;

    /* Clip heights at the 95th-percentile of |f(z)| so poles appear as
     * flat-topped cylinders rather than infinite spikes (Mathematica style). */
    double hcap = compute_height_cap(grid, stride * stride);

    /* 2 primitives per cell: color + Polygon */
    size_t cap   = (size_t)N * (size_t)N * 2 + 2;
    Expr** prims = malloc(sizeof(Expr*) * cap);
    size_t np    = 0;

    for (int iy = 0; iy < N; iy++) {
        for (int ix = 0; ix < N; ix++) {
            const CGrid* p00 = &grid[ iy      * (int)stride + ix    ];
            const CGrid* p10 = &grid[ iy      * (int)stride + ix + 1];
            const CGrid* p11 = &grid[(iy + 1) * (int)stride + ix + 1];
            const CGrid* p01 = &grid[(iy + 1) * (int)stride + ix    ];

            if (!p00->valid || !p10->valid || !p11->valid || !p01->valid)
                continue;

            /* Heights: |f(z)| clamped to hcap so poles become flat cylinders. */
#define HCLAMP(v) ((v) > hcap ? hcap : (v))
            double h00 = HCLAMP(sqrt(p00->re*p00->re + p00->im*p00->im));
            double h10 = HCLAMP(sqrt(p10->re*p10->re + p10->im*p10->im));
            double h11 = HCLAMP(sqrt(p11->re*p11->re + p11->im*p11->im));
            double h01 = HCLAMP(sqrt(p01->re*p01->re + p01->im*p01->im));
#undef HCLAMP

            /* Color at cell center */
            double re_avg = (p00->re + p10->re + p11->re + p01->re) * 0.25;
            double im_avg = (p00->im + p10->im + p11->im + p01->im) * 0.25;

            double re_sc = (re_avg - re_min) / re_span;
            double im_sc = (im_avg - im_min) / im_span;
            if (re_sc < 0.0) re_sc = 0.0; if (re_sc > 1.0) re_sc = 1.0;
            if (im_sc < 0.0) im_sc = 0.0; if (im_sc > 1.0) im_sc = 1.0;

            /* For 3D, the default color uses arg of the center but without
             * modulus-brightness attenuation (same visual weight as Plot3D). */
            Expr* color;
            if (!opts.color_function) {
                double arg = atan2(im_avg, re_avg);
                double t   = (arg + M_PI) / (2.0 * M_PI);
                if (t < 0.0) t = 0.0;
                if (t > 1.0) t = 1.0;
                double r, g, b;
                thermal_rgb(t, &r, &g, &b);
                Expr* ca[3] = { expr_new_real(r), expr_new_real(g), expr_new_real(b) };
                color = expr_new_function(expr_new_symbol(SYM_RGBColor), ca, 3);
            } else {
                color = cp_color(opts.color_function, opts.color_function_scaling,
                                  re_avg, im_avg, re_sc, im_sc);
            }
            prims[np++] = color;

            /* x axis = Re(z), y axis = Im(z), z axis = |f(z)| */
            double x0 = xmin + ix       * (xmax - xmin) / N;
            double x1 = xmin + (ix + 1) * (xmax - xmin) / N;
            double y0 = ymin + iy       * (ymax - ymin) / N;
            double y1 = ymin + (iy + 1) * (ymax - ymin) / N;

            Expr* verts[4] = { point3d(x0, y0, h00),
                               point3d(x1, y0, h10),
                               point3d(x1, y1, h11),
                               point3d(x0, y1, h01) };
            Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, 4);
            Expr* pargs[1] = { vlist };
            prims[np++] = expr_new_function(expr_new_symbol(SYM_Polygon), pargs, 1);
        }
    }

    free(grid);

    /* PlotRange for 3D: supply explicit z-extent so the renderer frames
     * the surface at the clamped ceiling rather than the raw pole height. */
    embed_plot_range3(xmin, xmax, ymin, ymax, 0.0, hcap, &pt, &pt_n);
    if (opts.show_legend) emit_phase_color_bar(opts.color_function, &pt, &pt_n);

    Expr* plist  = expr_new_function(expr_new_symbol(SYM_List), prims, np);
    free(prims);

    size_t gargc = 1 + pt_n;
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = plist;
    for (size_t i = 0; i < pt_n; i++) gargs[1 + i] = pt[i];
    free(pt);

    return expr_new_function(expr_new_symbol(SYM_Graphics3D), gargs, gargc);
}
