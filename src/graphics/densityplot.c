/* densityplot.c — DensityPlot[f, {x,xmin,xmax}, {y,ymin,ymax}, opts...]
 *
 * Renders f(x,y) as a heatmap: each grid cell is coloured by its average
 * function value via ColorFunction (default: thermal blue→yellow ramp from
 * plot_common's thermal_rgb). Returns Graphics[...] auto-displayed by the REPL.
 *
 * DensityPlot is HoldAll: f and iterator specs are held unevaluated until x
 * and y get numeric values, matching ContourPlot's semantics.
 *
 * Options:
 *   PlotPoints            grid resolution per axis (default 50)
 *   ColorFunction         f[t] → color, or "Rainbow" / "Temperature"
 *   ColorFunctionScaling  True (default): normalise z to [0,1]; False: raw z
 *   RegionFunction        f[x,y] mask; excluded cells are not drawn
 *   PlotLegends           Automatic: attach a $StreamColorBar colour scale
 *   Standard Graphics options pass through (Axes, AspectRatio→1, Frame, …) */

#include "densityplot.h"
#include "plot_common.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- options ---- */

typedef struct {
    int         plot_points;
    Expr*       color_function;         /* borrowed; NULL = thermal default */
    bool        color_function_scaling;
    Expr*       region_function;        /* borrowed; NULL = none */
    bool        show_legend;
    ScaleFnType sf_x, sf_y;
} DensityOpts;

static bool split_density_options(Expr* res, DensityOpts* o,
                                   Expr*** pt_out, size_t* pt_n_out) {
    o->plot_points            = 50;
    o->color_function         = NULL;
    o->color_function_scaling = true;
    o->region_function        = NULL;
    o->show_legend            = false;
    o->sf_x                   = SF_NONE;
    o->sf_y                   = SF_NONE;

    size_t argc = res->data.function.arg_count;
    size_t cap  = (argc > 3 ? argc - 3 : 0) + 4;
    Expr** pt   = malloc(sizeof(Expr*) * cap);
    size_t n    = 0;
    bool have_axes = false, have_aspect = false, have_frame = false;

#define DP_FAIL() do { free(pt); return false; } while(0)

    for (size_t i = 3; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) DP_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol : NULL;

        if (name == SYM_PlotPoints) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 2) DP_FAIL();
            o->plot_points = (int)v;
        } else if (name == SYM_ColorFunction) {
            o->color_function = rhs;          /* borrowed */
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = evaluate(expr_copy(rhs));
            o->color_function_scaling = !(v->type == EXPR_SYMBOL
                                          && v->data.symbol == SYM_False);
            expr_free(v);
        } else if (name == SYM_RegionFunction) {
            o->region_function = rhs;          /* borrowed */
        } else if (name == SYM_PlotLegends) {
            Expr* v = evaluate(expr_copy(rhs));
            o->show_legend = !(v->type == EXPR_SYMBOL
                               && (v->data.symbol == SYM_None
                                   || v->data.symbol == SYM_False));
            expr_free(v);
        } else if (name == SYM_ScalingFunctions) {
            Expr* v = evaluate(expr_copy(rhs));
            parse_scaling_functions(v, &o->sf_x, &o->sf_y);
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
            Expr* val    = evaluate(expr_copy(rhs));
            Expr* a[2]   = { expr_copy(lhs), val };
            pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    if (!have_axes && !have_frame) {
        Expr* fa[2] = { expr_new_symbol(SYM_Frame), expr_new_symbol(SYM_True) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), fa, 2);
        Expr* aa[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_False) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), aa, 2);
    } else if (!have_axes) {
        Expr* aa[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_False) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), aa, 2);
    }
    if (!have_aspect) {
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_integer(1) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }

    *pt_out   = pt;
    *pt_n_out = n;
    return true;
#undef DP_FAIL
}

/* ---- grid evaluation ---- */

static double dp_eval(Expr* xvar, Expr* yvar, Expr* body, double x, double y) {
    Expr* xv = expr_new_real(x);
    Expr* yv = expr_new_real(y);
    symtab_add_own_value(xvar->data.symbol, xvar, xv);
    symtab_add_own_value(yvar->data.symbol, yvar, yv);
    Expr* r = evaluate(body);
    double v;
    bool ok = expr_to_real_double(r, &v) && isfinite(v);
    expr_free(r); expr_free(xv); expr_free(yv);
    return ok ? v : NAN;
}

/* ---- colour resolution ---- *
 * t is already in [0,1] (normalised by the caller).
 * For a custom ColorFunction: try f[t] (1-arg), then fall back.
 * For "Rainbow": Hue sweep.  Default: thermal_rgb. */

static bool is_color_head_dp(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    return h == SYM_RGBColor || h == SYM_GrayLevel || h == SYM_Hue || h == SYM_CMYKColor;
}

static Expr* dp_color(Expr* cfn, double t) {
    if (!cfn) {
        double rv, gv, bv;
        thermal_rgb(t, &rv, &gv, &bv);
        Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
        return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
    }
    if (cfn->type == EXPR_STRING) {
        Expr* c = named_color_ramp(cfn->data.string, t);
        if (c) return c;
    }
    /* Custom function: try f[t]. */
    Expr* a1[1] = { expr_new_real(t) };
    Expr* call  = expr_new_function(expr_copy(cfn), a1, 1);
    Expr* r     = evaluate(call);
    expr_free(call);
    if (is_color_head_dp(r)) return r;
    expr_free(r);
    /* Fallback: neutral mid-gray. */
    Expr* ga[1] = { expr_new_real(0.5) };
    return expr_new_function(expr_new_symbol(SYM_GrayLevel), ga, 1);
}

/* ---- main builtin ---- */

Expr* builtin_densityplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 3) return NULL;

    Expr* body = res->data.function.args[0]; /* held */

    IterSpec xspec;
    if (!iter_spec_parse(res->data.function.args[1], &xspec)) return NULL;
    if (xspec.kind != ITER_KIND_RANGE || !xspec.var) { iter_spec_free(&xspec); return NULL; }
    double xmin, xmax;
    if (!numericize_bound(xspec.imin, &xmin) || !numericize_bound(xspec.imax, &xmax)
        || !(xmin < xmax)) { iter_spec_free(&xspec); return NULL; }

    IterSpec yspec;
    if (!iter_spec_parse(res->data.function.args[2], &yspec)) {
        iter_spec_free(&xspec); return NULL; }
    if (yspec.kind != ITER_KIND_RANGE || !yspec.var) {
        iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL; }
    double ymin, ymax;
    if (!numericize_bound(yspec.imin, &ymin) || !numericize_bound(yspec.imax, &ymax)
        || !(ymin < ymax)) { iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL; }

    DensityOpts opts;
    Expr** pt = NULL; size_t pt_n = 0;
    if (!split_density_options(res, &opts, &pt, &pt_n)) {
        iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL; }

    int    N  = opts.plot_points;
    /* Sample uniformly in world (scaled) space */
    double u_xmin = scale_apply(opts.sf_x, xmin), u_xmax = scale_apply(opts.sf_x, xmax);
    double u_ymin = scale_apply(opts.sf_y, ymin), u_ymax = scale_apply(opts.sf_y, ymax);
    double du_x = (u_xmax - u_xmin) / N;
    double du_y = (u_ymax - u_ymin) / N;

    double* grid = malloc(sizeof(double) * (size_t)(N + 1) * (size_t)(N + 1));
    if (!grid) {
        iter_spec_free(&xspec); iter_spec_free(&yspec);
        for (size_t i = 0; i < pt_n; i++) expr_free(pt[i]); free(pt);
        return NULL;
    }

    Rule* old_x = iter_spec_shadow(xspec.var);
    Rule* old_y = iter_spec_shadow(yspec.var);

    double zmin =  1e300, zmax = -1e300;
    for (int iy = 0; iy <= N; iy++) {
        double uy = u_ymin + iy * du_y;
        double y  = scale_invert(opts.sf_y, uy);
        for (int ix = 0; ix <= N; ix++) {
            double ux = u_xmin + ix * du_x;
            double x  = scale_invert(opts.sf_x, ux);
            double v  = dp_eval(xspec.var, yspec.var, body, x, y);
            grid[iy * (N + 1) + ix] = v;
            if (isfinite(v)) {
                if (v < zmin) zmin = v;
                if (v > zmax) zmax = v;
            }
        }
    }

    iter_spec_restore(xspec.var, old_x);
    iter_spec_restore(yspec.var, old_y);

    double zspan = (zmax > zmin) ? (zmax - zmin) : 1.0;

    /* 2 primitives per cell (colour directive + Rectangle). */
    size_t cap   = (size_t)N * N * 2 + 4;
    Expr** prims = malloc(sizeof(Expr*) * cap);
    size_t np    = 0;

    for (int iy = 0; iy < N; iy++) {
        double uy0 = u_ymin + iy * du_y;
        for (int ix = 0; ix < N; ix++) {
            double ux0 = u_xmin + ix * du_x;
            double v00 = grid[ iy      * (N + 1) + ix    ];
            double v10 = grid[ iy      * (N + 1) + ix + 1];
            double v11 = grid[(iy + 1) * (N + 1) + ix + 1];
            double v01 = grid[(iy + 1) * (N + 1) + ix    ];

            if (!isfinite(v00) || !isfinite(v10) || !isfinite(v11) || !isfinite(v01))
                continue;

            /* Region check uses original data coordinates at cell center */
            double cx = scale_invert(opts.sf_x, ux0 + du_x * 0.5);
            double cy = scale_invert(opts.sf_y, uy0 + du_y * 0.5);
            if (opts.region_function && !eval_region(opts.region_function, cx, cy))
                continue;

            double avg = (v00 + v10 + v11 + v01) * 0.25;
            double t   = opts.color_function_scaling
                       ? (avg - zmin) / zspan : avg;
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;

            prims[np++] = dp_color(opts.color_function, t);

            /* Rectangle in world coordinates */
            Expr* p1[2] = { expr_new_real(ux0),         expr_new_real(uy0) };
            Expr* p2[2] = { expr_new_real(ux0 + du_x),  expr_new_real(uy0 + du_y) };
            Expr* ra[2] = { expr_new_function(expr_new_symbol(SYM_List), p1, 2),
                             expr_new_function(expr_new_symbol(SYM_List), p2, 2) };
            prims[np++] = expr_new_function(expr_new_symbol(SYM_Rectangle), ra, 2);
        }
    }

    free(grid);

    /* Embed PlotRange matching the domain (renderer uses it for axis bounds). */
    {
        bool have_pr = false;
        for (size_t i = 0; i < pt_n; i++) {
            const Expr* e = pt[i];
            if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
                && e->data.function.args[0]->type == EXPR_SYMBOL
                && e->data.function.args[0]->data.symbol == SYM_PlotRange)
                { have_pr = true; break; }
        }
        if (!have_pr) {
            pt = realloc(pt, sizeof(Expr*) * (pt_n + 1));
            Expr* xr[2]  = { expr_new_real(u_xmin), expr_new_real(u_xmax) };
            Expr* yr[2]  = { expr_new_real(u_ymin), expr_new_real(u_ymax) };
            Expr* xrng   = expr_new_function(expr_new_symbol(SYM_List), xr, 2);
            Expr* yrng   = expr_new_function(expr_new_symbol(SYM_List), yr, 2);
            Expr* rl[2]  = { xrng, yrng };
            Expr* pr     = expr_new_function(expr_new_symbol(SYM_List), rl, 2);
            Expr* ra[2]  = { expr_new_symbol(SYM_PlotRange), pr };
            pt[pt_n++]   = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
    }

    emit_scaling_meta(opts.sf_x, opts.sf_y, &pt, &pt_n);

    /* PlotLegends -> Automatic: attach $StreamColorBar so the renderer draws
     * a vertical colour scale on the right of the plot. */
    if (opts.show_legend) {
        double bar_lo = isfinite(zmin) ? zmin : 0.0;
        double bar_hi = isfinite(zmax) ? zmax : 1.0;
        if (bar_lo == bar_hi) bar_hi = bar_lo + 1.0;
        pt = realloc(pt, sizeof(Expr*) * (pt_n + 1));
        Expr* cfn_copy = opts.color_function
                         ? expr_copy(opts.color_function)
                         : expr_new_symbol(SYM_Automatic);
        Expr* cb_args[3] = { expr_new_real(bar_lo), expr_new_real(bar_hi), cfn_copy };
        pt[pt_n++] = expr_new_function(expr_new_symbol(SYM_StreamColorBar), cb_args, 3);
    }

    Expr* plist = expr_new_function(expr_new_symbol(SYM_List), prims, np);
    free(prims);

    size_t  gargc = 1 + pt_n;
    Expr**  gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = plist;
    for (size_t i = 0; i < pt_n; i++) gargs[1 + i] = pt[i];
    free(pt);

    Expr* g = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
    free(gargs);

    iter_spec_free(&xspec);
    iter_spec_free(&yspec);
    return g;
}
