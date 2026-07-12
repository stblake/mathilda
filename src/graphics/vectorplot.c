/* vectorplot.c — VectorPlot[{vx, vy}, {x,xmin,xmax}, {y,ymin,ymax}, opts...]
 *
 * Draws a grid of arrows showing the direction (and optionally magnitude) of
 * the 2-D vector field {vx, vy}. Returns a Graphics[...] object auto-
 * displayed by the REPL.
 *
 * VectorPlot is HoldAll: vx, vy, and the iterator specs are held unevaluated
 * until x and y are bound to numeric values (same semantics as ContourPlot).
 *
 * Options:
 *   VectorPoints    integer n → n × n grid (default 15); Automatic = 15
 *   VectorScale     Automatic: normalise all arrows to equal display length
 *                   None:      proportional to magnitude
 *                   real f:    arrow length = f × grid_spacing
 *   VectorStyle     style directive(s) applied to all arrows
 *   ColorFunction   f[vx,vy,speed,x,y] (or fewer args) → color, or "Rainbow"
 *   ColorFunctionScaling  True (default): normalise speed to [0,1]
 *   RegionFunction  f[x,y] mask: skip grid points outside the region
 *   Standard Graphics options pass through (Axes, AspectRatio→1, Frame, …) */

#include "vectorplot.h"
#include "plot_common.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- options ---- */

typedef enum { VS_AUTOMATIC, VS_NONE, VS_FRACTION } VecScaleMode;

typedef struct {
    int          vector_points;
    VecScaleMode scale_mode;
    double       scale_frac;       /* VS_FRACTION only */
    Expr*        vector_style;     /* borrowed; NULL = default palette */
    Expr*        color_function;   /* borrowed; NULL = none */
    bool         color_function_scaling;
    Expr*        region_function;  /* borrowed; NULL = none */
} VecOpts;

static bool split_vector_options(Expr* res, VecOpts* o,
                                  Expr*** pt_out, size_t* pt_n_out) {
    o->vector_points          = 15;
    o->scale_mode             = VS_AUTOMATIC;
    o->scale_frac             = 0.5;
    o->vector_style           = NULL;
    o->color_function         = NULL;
    o->color_function_scaling = true;
    o->region_function        = NULL;

    size_t argc = res->data.function.arg_count;
    size_t cap  = (argc > 3 ? argc - 3 : 0) + 4;
    Expr** pt   = malloc(sizeof(Expr*) * cap);
    size_t n    = 0;
    bool have_axes = false, have_aspect = false;

#define VP_FAIL() do { free(pt); return false; } while(0)

    for (size_t i = 3; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) VP_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol : NULL;

        if (name == SYM_VectorPoints) {
            Expr* v = evaluate(expr_copy(rhs));
            if (v->type == EXPR_SYMBOL && v->data.symbol == SYM_Automatic) {
                o->vector_points = 15; expr_free(v);
            } else {
                long lv;
                expr_free(v);
                if (!parse_long_value(rhs, &lv) || lv < 1) VP_FAIL();
                o->vector_points = (int)lv;
            }
        } else if (name == SYM_VectorScale) {
            Expr* v = evaluate(expr_copy(rhs));
            if (v->type == EXPR_SYMBOL && v->data.symbol == SYM_Automatic) {
                o->scale_mode = VS_AUTOMATIC;
            } else if (v->type == EXPR_SYMBOL && v->data.symbol == SYM_None) {
                o->scale_mode = VS_NONE;
            } else {
                double f;
                if (expr_to_real_double(v, &f) && f > 0.0) {
                    o->scale_mode  = VS_FRACTION;
                    o->scale_frac  = f;
                } else {
                    o->scale_mode = VS_AUTOMATIC;
                }
            }
            expr_free(v);
        } else if (name == SYM_VectorStyle) {
            o->vector_style = rhs;   /* borrowed */
        } else if (name == SYM_ColorFunction) {
            o->color_function = rhs; /* borrowed */
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = evaluate(expr_copy(rhs));
            o->color_function_scaling = !(v->type == EXPR_SYMBOL
                                          && v->data.symbol == SYM_False);
            expr_free(v);
        } else if (name == SYM_RegionFunction) {
            o->region_function = rhs; /* borrowed */
        } else {
            if (name == SYM_Axes)        have_axes   = true;
            if (name == SYM_AspectRatio) have_aspect = true;
            Expr* val  = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    if (!have_axes) {
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
#undef VP_FAIL
}

/* ---- vector field evaluation at (x,y) ---- */

typedef struct { double vx, vy; bool ok; } Vec2;

static Vec2 vp_eval(Expr* xvar, Expr* yvar,
                     Expr* vx_body, Expr* vy_body,
                     double x, double y) {
    Expr* xv = expr_new_real(x);
    Expr* yv = expr_new_real(y);
    symtab_add_own_value(xvar->data.symbol, xvar, xv);
    symtab_add_own_value(yvar->data.symbol, yvar, yv);

    Expr* rvx = evaluate(vx_body);
    Expr* rvy = evaluate(vy_body);
    double vx_val = 0.0, vy_val = 0.0;
    bool ok = expr_to_real_double(rvx, &vx_val) && isfinite(vx_val)
           && expr_to_real_double(rvy, &vy_val) && isfinite(vy_val);
    expr_free(rvx); expr_free(rvy);
    expr_free(xv); expr_free(yv);
    return (Vec2){ vx_val, vy_val, ok };
}

/* ---- arrow color ---- */

static bool is_color_head_vp(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    return h == SYM_RGBColor || h == SYM_GrayLevel || h == SYM_Hue || h == SYM_CMYKColor;
}

/* speed_norm is pre-normalised [0,1]; raw_vx/vy/speed are used when the
 * ColorFunction wants physical values. */
static Expr* vp_color(Expr* cfn, double speed_norm,
                       double vx, double vy, double speed) {
    if (!cfn) {
        /* Default: thermal ramp based on speed. */
        double r, g, b;
        thermal_rgb(speed_norm, &r, &g, &b);
        Expr* a[3] = { expr_new_real(r), expr_new_real(g), expr_new_real(b) };
        return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
    }
    if (cfn->type == EXPR_STRING) {
        Expr* c = named_color_ramp(cfn->data.string, speed_norm);
        if (c) return c;
    }
    /* Try f[vx,vy,speed], f[speed], f[speed_norm] in sequence. */
    static const int arities[] = { 3, 1 };
    for (size_t ai = 0; ai < 2; ai++) {
        int ar = arities[ai];
        Expr** fargs = malloc(sizeof(Expr*) * (size_t)ar);
        if (ar == 3) {
            fargs[0] = expr_new_real(vx);
            fargs[1] = expr_new_real(vy);
            fargs[2] = expr_new_real(speed);
        } else {
            fargs[0] = expr_new_real(speed_norm);
        }
        Expr* call = expr_new_function(expr_copy(cfn), fargs, (size_t)ar);
        free(fargs);
        Expr* r = evaluate(call);
        expr_free(call);
        if (is_color_head_vp(r)) return r;
        expr_free(r);
    }
    /* Fallback: mid-gray. */
    Expr* ga[1] = { expr_new_real(0.5) };
    return expr_new_function(expr_new_symbol(SYM_GrayLevel), ga, 1);
}

/* ---- main builtin ---- */

Expr* builtin_vectorplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 4) return NULL;

    /* args[0] must be a List[vx_body, vy_body] */
    Expr* field = res->data.function.args[0];
    if (!field || field->type != EXPR_FUNCTION
        || field->data.function.head->type != EXPR_SYMBOL
        || field->data.function.head->data.symbol != SYM_List
        || field->data.function.arg_count != 2)
        return NULL;

    Expr* vx_body = field->data.function.args[0];  /* held */
    Expr* vy_body = field->data.function.args[1];  /* held */

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

    VecOpts opts;
    Expr** pt = NULL; size_t pt_n = 0;
    if (!split_vector_options(res, &opts, &pt, &pt_n)) {
        iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL; }

    int    N   = opts.vector_points;
    double dx  = (xmax - xmin) / (N - 1 > 0 ? N - 1 : 1);
    double dy  = (ymax - ymin) / (N - 1 > 0 ? N - 1 : 1);

    /* Evaluate the field on an N×N grid; collect vectors. */
    typedef struct { double x, y, vx, vy; } Arrow;
    Arrow* arrows  = malloc(sizeof(Arrow) * (size_t)N * N);
    size_t n_arr   = 0;
    double spd_max = 0.0;

    Rule* old_x = iter_spec_shadow(xspec.var);
    Rule* old_y = iter_spec_shadow(yspec.var);

    for (int iy = 0; iy < N; iy++) {
        double y = ymin + iy * dy;
        for (int ix = 0; ix < N; ix++) {
            double x = xmin + ix * dx;
            if (opts.region_function && !eval_region(opts.region_function, x, y)) continue;
            Vec2 v = vp_eval(xspec.var, yspec.var, vx_body, vy_body, x, y);
            if (!v.ok) continue;
            double spd = sqrt(v.vx * v.vx + v.vy * v.vy);
            if (spd > spd_max) spd_max = spd;
            arrows[n_arr++] = (Arrow){ x, y, v.vx, v.vy };
        }
    }

    iter_spec_restore(xspec.var, old_x);
    iter_spec_restore(yspec.var, old_y);

    /* Arrow length scale: uniform grid spacing used as reference. */
    double grid_step = (dx < dy ? dx : dy);

    /* 3 prims per arrow: color directive + Thickness + Arrow. */
    size_t cap   = n_arr * 4 + 4;
    Expr** prims = malloc(sizeof(Expr*) * cap);
    size_t np    = 0;

    /* Apply VectorStyle as a leading directive if set. */
    if (opts.vector_style) {
        Expr* vs = evaluate(expr_copy(opts.vector_style));
        prims[np++] = vs;
    } else {
        /* Default thin line weight so arrows don't look blocky. */
        Expr* th_a[1] = { expr_new_real(0.005) };
        prims[np++] = expr_new_function(expr_new_symbol(SYM_Thickness), th_a, 1);
    }

    double spd_range = (spd_max > 0.0) ? spd_max : 1.0;

    for (size_t k = 0; k < n_arr; k++) {
        Arrow* a   = &arrows[k];
        double spd = sqrt(a->vx * a->vx + a->vy * a->vy);
        double snorm = opts.color_function_scaling
                     ? (spd / spd_range) : spd;
        if (snorm < 0.0) snorm = 0.0;
        if (snorm > 1.0) snorm = 1.0;

        prims[np++] = vp_color(opts.color_function, snorm, a->vx, a->vy, spd);

        /* Arrow half-length in plot coordinates. */
        double len;
        if (opts.scale_mode == VS_NONE) {
            /* Proportional: normalise by spd_max then scale to grid_step/2. */
            len = (spd_max > 0.0) ? (spd / spd_max) * grid_step * 0.5 : 0.0;
        } else if (opts.scale_mode == VS_FRACTION) {
            len = opts.scale_frac * grid_step;
        } else {
            /* Automatic: uniform length = grid_step * 0.5; show direction only. */
            len = grid_step * 0.5;
        }

        double nx = 0.0, ny = 0.0;
        if (spd > 0.0) { nx = a->vx / spd; ny = a->vy / spd; }
        else           { len = 0.0; }

        /* Centre the arrow on the grid point: tail at (x - len/2 * n),
         * head at (x + len/2 * n). */
        double tx = a->x - nx * len * 0.5, ty = a->y - ny * len * 0.5;
        double hx = a->x + nx * len * 0.5, hy = a->y + ny * len * 0.5;

        Expr* tail_pt[2] = { expr_new_real(tx), expr_new_real(ty) };
        Expr* head_pt[2] = { expr_new_real(hx), expr_new_real(hy) };
        Expr* tail = expr_new_function(expr_new_symbol(SYM_List), tail_pt, 2);
        Expr* head = expr_new_function(expr_new_symbol(SYM_List), head_pt, 2);
        Expr* seg[2] = { tail, head };
        Expr* seg_list = expr_new_function(expr_new_symbol(SYM_List), seg, 2);
        Expr* arr_args[1] = { seg_list };
        prims[np++] = expr_new_function(expr_new_symbol(SYM_Arrow), arr_args, 1);
    }

    free(arrows);

    /* Embed explicit PlotRange. */
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
            Expr* xr[2] = { expr_new_real(xmin), expr_new_real(xmax) };
            Expr* yr[2] = { expr_new_real(ymin), expr_new_real(ymax) };
            Expr* xrng  = expr_new_function(expr_new_symbol(SYM_List), xr, 2);
            Expr* yrng  = expr_new_function(expr_new_symbol(SYM_List), yr, 2);
            Expr* rl[2] = { xrng, yrng };
            Expr* pr    = expr_new_function(expr_new_symbol(SYM_List), rl, 2);
            Expr* ra[2] = { expr_new_symbol(SYM_PlotRange), pr };
            pt[pt_n++]  = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
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
