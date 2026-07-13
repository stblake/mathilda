/* plot_common.c — see plot_common.h. Extracted from plot.c so Plot3D
 * (plot3d.c) reuses the exact same option/evaluation idioms as Plot. */

#include "plot_common.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "print.h"
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <string.h>
#include <math.h>

bool expr_to_real_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)   { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational) {
        Expr* n = e->data.function.args[0];
        Expr* d = e->data.function.args[1];
        if (n->type == EXPR_INTEGER && d->type == EXPR_INTEGER && d->data.integer != 0) {
            *out = (double)n->data.integer / (double)d->data.integer;
            return true;
        }
    }
    return false;
}

bool numericize_bound(Expr* e, double* out) {
    Expr* n_arg[1] = { expr_copy(e) };
    Expr* n_call = expr_new_function(expr_new_symbol("N"), n_arg, 1);
    Expr* result = evaluate(n_call);
    expr_free(n_call);
    bool ok = expr_to_real_double(result, out) && isfinite(*out);
    expr_free(result);
    return ok;
}

bool is_rule_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    return (h->data.symbol.name == SYM_Rule || h->data.symbol.name == SYM_RuleDelayed)
        && e->data.function.arg_count == 2;
}

bool parse_long_value(Expr* rhs, long* out) {
    Expr* v = evaluate(expr_copy(rhs));
    bool ok = (v->type == EXPR_INTEGER);
    if (ok) *out = (long)v->data.integer;
    expr_free(v);
    return ok;
}

Expr* palette_color(size_t i) {
    static const double pal[][3] = {
        { 0.368417, 0.506779, 0.709798 },
        { 0.880722, 0.611041, 0.142051 },
        { 0.560181, 0.691569, 0.194885 },
        { 0.922526, 0.385626, 0.209179 },
        { 0.528488, 0.470624, 0.701351 },
        { 0.772079, 0.431554, 0.102387 },
        { 0.363898, 0.618501, 0.782349 },
        { 1.000000, 0.750000, 0.000000 },
        { 0.647624, 0.378160, 0.614037 },
        { 0.571589, 0.586483, 0.000000 },
    };
    size_t k = i % (sizeof(pal) / sizeof(pal[0]));
    Expr* a[3] = { expr_new_real(pal[k][0]), expr_new_real(pal[k][1]), expr_new_real(pal[k][2]) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
}

bool eval_region(Expr* region_fn, double x, double y) {
    Expr* args2[2] = { expr_new_real(x), expr_new_real(y) };
    Expr* call2 = expr_new_function(expr_copy(region_fn), args2, 2);
    Expr* r2 = evaluate(call2);
    expr_free(call2); /* evaluate() borrows its argument; the call node is ours to free */
    bool true2 = (r2->type == EXPR_SYMBOL && r2->data.symbol.name == SYM_True);
    bool false2 = (r2->type == EXPR_SYMBOL && r2->data.symbol.name == SYM_False);
    expr_free(r2);
    if (true2) return true;
    if (false2) return false;

    /* The 2-arg call didn't resolve to a boolean (likely a 1-arg function,
     * e.g. Function[x, x > 0]) -- retry with just x. */
    Expr* args1[1] = { expr_new_real(x) };
    Expr* call1 = expr_new_function(expr_copy(region_fn), args1, 1);
    Expr* r1 = evaluate(call1);
    expr_free(call1);
    bool ok = (r1->type == EXPR_SYMBOL && r1->data.symbol.name == SYM_True);
    expr_free(r1);
    return ok;
}

static bool is_color_head(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_RGBColor
            || e->data.function.head->data.symbol.name == SYM_GrayLevel
            || e->data.function.head->data.symbol.name == SYM_Hue);
}

Expr* eval_color_function(Expr* color_fn, double x, double y,
                           double xmin, double xmax, bool scaling) {
    if (color_fn->type == EXPR_STRING && strcmp(color_fn->data.string, "Rainbow") == 0) {
        double t = (xmax > xmin) ? (x - xmin) / (xmax - xmin) : 0.0;
        Expr* a[1] = { expr_new_real(t * 0.8) }; /* stop short of wrapping back to red */
        return expr_new_function(expr_new_symbol(SYM_Hue), a, 1);
    }

    double cx = (scaling && xmax > xmin) ? (x - xmin) / (xmax - xmin) : x;

    Expr* args2[2] = { expr_new_real(cx), expr_new_real(y) };
    Expr* call2 = expr_new_function(expr_copy(color_fn), args2, 2);
    Expr* r2 = evaluate(call2);
    expr_free(call2);
    if (is_color_head(r2)) return r2;
    expr_free(r2);

    Expr* args1[1] = { expr_new_real(cx) };
    Expr* call1 = expr_new_function(expr_copy(color_fn), args1, 1);
    Expr* r1 = evaluate(call1);
    expr_free(call1);
    if (is_color_head(r1)) return r1;
    expr_free(r1);

    Expr* a[1] = { expr_new_real(0.5) };
    return expr_new_function(expr_new_symbol(SYM_GrayLevel), a, 1);
}

Expr* eval_color_function3(Expr* color_fn,
                            double x,    double y,    double z,
                            double xmin, double xmax,
                            double ymin, double ymax,
                            double zmin, double zmax,
                            bool scaling) {
    if (color_fn->type == EXPR_STRING && strcmp(color_fn->data.string, "Rainbow") == 0) {
        /* Height-based: cold blue → hot red, capped at 0.8 hue to avoid
         * wrapping back to red from the violet end. */
        double t = (zmax > zmin) ? (z - zmin) / (zmax - zmin) : 0.0;
        Expr* a[1] = { expr_new_real((1.0 - t) * 0.8) };
        return expr_new_function(expr_new_symbol(SYM_Hue), a, 1);
    }

    double xs = (scaling && xmax > xmin) ? (x - xmin) / (xmax - xmin) : x;
    double ys = (scaling && ymax > ymin) ? (y - ymin) / (ymax - ymin) : y;
    double zs = (scaling && zmax > zmin) ? (z - zmin) / (zmax - zmin) : z;

    /* Try f[xs, ys, zs] first (Mathematica's Plot3D ColorFunction convention). */
    Expr* args3[3] = { expr_new_real(xs), expr_new_real(ys), expr_new_real(zs) };
    Expr* call3 = expr_new_function(expr_copy(color_fn), args3, 3);
    Expr* r3 = evaluate(call3);
    expr_free(call3);
    if (is_color_head(r3)) return r3;
    expr_free(r3);

    /* Fall back to f[xs, zs] — common for height-coloured ramps. */
    Expr* args2[2] = { expr_new_real(xs), expr_new_real(zs) };
    Expr* call2 = expr_new_function(expr_copy(color_fn), args2, 2);
    Expr* r2 = evaluate(call2);
    expr_free(call2);
    if (is_color_head(r2)) return r2;
    expr_free(r2);

    /* Last resort: f[zs] — a univariate colour ramp indexed by height. */
    Expr* args1[1] = { expr_new_real(zs) };
    Expr* call1 = expr_new_function(expr_copy(color_fn), args1, 1);
    Expr* r1 = evaluate(call1);
    expr_free(call1);
    if (is_color_head(r1)) return r1;
    expr_free(r1);

    Expr* a[1] = { expr_new_real(0.5) };
    return expr_new_function(expr_new_symbol(SYM_GrayLevel), a, 1);
}

/* build_legend_meta — shared by Plot, ParametricPlot, and Plot3D.
 *
 * `legends`: already-evaluated PlotLegends value (Automatic, "Expressions",
 *   or an explicit {label1, label2, ...} List).  If NULL or resolves to None,
 *   returns NULL (no legend).
 * `bodies`: the nfun body/curve-spec expressions used as auto-labels.
 * `single_color`: the resolved PlotStyle color for single-curve plots
 *   (NULL falls back to palette_color(0)).
 *
 * Returns a fresh $PlotLegendData[{color1,label1}, ...] node, or NULL.
 * The renderer (render.c: draw_legend) reads it at display time. */
Expr* build_legend_meta(Expr* legends, Expr** bodies, size_t nfun, Expr* single_color) {
    if (!legends) return NULL;
    if (legends->type == EXPR_SYMBOL && legends->data.symbol.name == SYM_None) return NULL;
    if (nfun == 0) return NULL;

    bool explicit_list = (legends->type == EXPR_FUNCTION
                          && legends->data.function.head->type == EXPR_SYMBOL
                          && legends->data.function.head->data.symbol.name == SYM_List);
    bool multi = (nfun > 1);

    Expr** entries = malloc(sizeof(Expr*) * nfun);
    for (size_t i = 0; i < nfun; i++) {
        Expr* color = multi ? palette_color(i)
                            : (single_color ? expr_copy(single_color) : palette_color(0));
        Expr* label;
        if (explicit_list && i < legends->data.function.arg_count) {
            label = expr_copy(legends->data.function.args[i]);
        } else {
            char* s = (bodies && bodies[i]) ? expr_to_string(bodies[i]) : NULL;
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

/* thermal_rgb — 5-stop RGB interpolation approximating Mathematica's default
 * StreamPlot speed colormap: dark blue-purple at t=0, bright yellow at t=1. */
void thermal_rgb(double t, double* r, double* g, double* b) {
    static const double stops[5][3] = {
        { 0.04, 0.00, 0.30 }, /* dark blue-purple */
        { 0.40, 0.00, 0.60 }, /* purple           */
        { 0.80, 0.10, 0.20 }, /* red              */
        { 1.00, 0.55, 0.00 }, /* orange           */
        { 1.00, 0.95, 0.20 }, /* bright yellow    */
    };
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    double idx = t * 4.0;
    int    i   = (int)idx;
    if (i > 3) i = 3;
    double f = idx - (double)i;
    *r = stops[i][0] + f * (stops[i + 1][0] - stops[i][0]);
    *g = stops[i][1] + f * (stops[i + 1][1] - stops[i][1]);
    *b = stops[i][2] + f * (stops[i + 1][2] - stops[i][2]);
}
