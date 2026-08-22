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
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
            || e->data.function.head->data.symbol.name == SYM_Hue
            || e->data.function.head->data.symbol.name == SYM_CMYKColor);
}

Expr* eval_color_function(Expr* color_fn, double x, double y,
                           double xmin, double xmax, bool scaling) {
    if (color_fn->type == EXPR_STRING) {
        double t = (xmax > xmin) ? (x - xmin) / (xmax - xmin) : 0.0;
        Expr* c = named_color_ramp(color_fn->data.string, t);
        if (c) return c;
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
    if (color_fn->type == EXPR_STRING) {
        double t = (zmax > zmin) ? (z - zmin) / (zmax - zmin) : 0.0;
        if (strcmp(color_fn->data.string, "Rainbow") == 0) {
            /* 3D Rainbow: inverted height-based sweep (cold blue at top → red at bottom). */
            Expr* a[1] = { expr_new_real((1.0 - t) * 0.8) };
            return expr_new_function(expr_new_symbol(SYM_Hue), a, 1);
        }
        Expr* c = named_color_ramp(color_fn->data.string, t);
        if (c) return c;
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

/* cool_tones_rgb — near-white ice blue (t=0) → deep navy/indigo (t=1). */
void cool_tones_rgb(double t, double* r, double* g, double* b) {
    static const double stops[5][3] = {
        { 0.91, 0.96, 1.00 }, /* near-white, ice-blue tint */
        { 0.53, 0.78, 0.96 }, /* light sky blue            */
        { 0.18, 0.50, 0.83 }, /* cornflower blue           */
        { 0.07, 0.23, 0.60 }, /* royal/cobalt blue         */
        { 0.02, 0.08, 0.35 }, /* deep navy/indigo          */
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

/* warm_tones_rgb — pale cream (t=0) → amber → orange → deep crimson (t=1). */
void warm_tones_rgb(double t, double* r, double* g, double* b) {
    static const double stops[5][3] = {
        { 1.00, 0.97, 0.80 }, /* pale cream/yellow */
        { 1.00, 0.80, 0.38 }, /* warm amber        */
        { 0.95, 0.48, 0.08 }, /* orange            */
        { 0.78, 0.14, 0.06 }, /* red-orange        */
        { 0.45, 0.03, 0.03 }, /* deep crimson      */
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

/* ---------------------------------------------------------------------- */
/* Perceptually-uniform colormaps (matplotlib/viridisLite) + "Haze".       */
/*                                                                          */
/* Each of viridis/magma/plasma/inferno/cividis is a 32-stop resampling of  */
/* the authoritative 256-entry matplotlib table (linear-interp error is     */
/* sub-perceptual). "Haze" is Sam Blake's palette, reproduced exactly from  */
/* its LinearSegmentedColormap definition: 7 anchors at non-uniform         */
/* positions (white→pink→blue→green→yellow→orange→red).                     */
/* ---------------------------------------------------------------------- */

static const double viridis_stops[32][3] = {
    { 0.267004, 0.004874, 0.329415 }, { 0.277018, 0.050344, 0.375715 },
    { 0.282327, 0.094955, 0.417331 }, { 0.282623, 0.140926, 0.457517 },
    { 0.278012, 0.180367, 0.486697 }, { 0.269308, 0.218818, 0.509577 },
    { 0.257322, 0.256130, 0.526563 }, { 0.241237, 0.296485, 0.539709 },
    { 0.225863, 0.330805, 0.547314 }, { 0.210503, 0.363727, 0.552206 },
    { 0.195860, 0.395433, 0.555276 }, { 0.182256, 0.426184, 0.557120 },
    { 0.168126, 0.459988, 0.558082 }, { 0.156270, 0.489624, 0.557936 },
    { 0.144759, 0.519093, 0.556572 }, { 0.133743, 0.548535, 0.553541 },
    { 0.123463, 0.581687, 0.547445 }, { 0.119423, 0.611141, 0.538982 },
    { 0.124780, 0.640461, 0.527068 }, { 0.143303, 0.669459, 0.511215 },
    { 0.180653, 0.701402, 0.488189 }, { 0.226397, 0.728888, 0.462789 },
    { 0.281477, 0.755203, 0.432552 }, { 0.344074, 0.780029, 0.397381 },
    { 0.412913, 0.803041, 0.357269 }, { 0.496615, 0.826376, 0.306377 },
    { 0.575563, 0.844566, 0.256415 }, { 0.657642, 0.860219, 0.203082 },
    { 0.741388, 0.873449, 0.149561 }, { 0.835270, 0.886029, 0.102646 },
    { 0.916242, 0.896091, 0.100717 }, { 0.993248, 0.906157, 0.143936 },
};

static const double magma_stops[32][3] = {
    { 0.001462, 0.000466, 0.013866 }, { 0.013708, 0.011771, 0.068667 },
    { 0.039608, 0.031090, 0.133515 }, { 0.078815, 0.054184, 0.211667 },
    { 0.118405, 0.066479, 0.286321 }, { 0.165308, 0.067911, 0.361816 },
    { 0.218512, 0.061158, 0.425392 }, { 0.278493, 0.061978, 0.469190 },
    { 0.329114, 0.075972, 0.489287 }, { 0.378211, 0.095332, 0.500067 },
    { 0.426877, 0.115395, 0.505714 }, { 0.475780, 0.134577, 0.507921 },
    { 0.531507, 0.154739, 0.506895 }, { 0.581819, 0.171596, 0.502777 },
    { 0.632805, 0.187893, 0.495332 }, { 0.684224, 0.204286, 0.484219 },
    { 0.742004, 0.224025, 0.467018 }, { 0.792427, 0.244242, 0.447543 },
    { 0.840636, 0.268953, 0.424666 }, { 0.884651, 0.300530, 0.400047 },
    { 0.925937, 0.346844, 0.374959 }, { 0.953099, 0.397563, 0.361438 },
    { 0.971582, 0.454210, 0.361030 }, { 0.983485, 0.513280, 0.374198 },
    { 0.990871, 0.572706, 0.398714 }, { 0.995480, 0.639027, 0.436607 },
    { 0.997186, 0.697349, 0.477182 }, { 0.997138, 0.755190, 0.522806 },
    { 0.995680, 0.812706, 0.572645 }, { 0.992831, 0.877168, 0.633109 },
    { 0.989815, 0.934329, 0.690198 }, { 0.987053, 0.991438, 0.749504 },
};

static const double plasma_stops[32][3] = {
    { 0.050383, 0.029803, 0.527975 }, { 0.132381, 0.022258, 0.563250 },
    { 0.193374, 0.018354, 0.590330 }, { 0.254627, 0.013882, 0.615419 },
    { 0.306210, 0.008902, 0.633694 }, { 0.356359, 0.003798, 0.647810 },
    { 0.405503, 0.000678, 0.656977 }, { 0.459623, 0.003574, 0.660277 },
    { 0.506454, 0.016333, 0.656202 }, { 0.551715, 0.043136, 0.645277 },
    { 0.595011, 0.077190, 0.627917 }, { 0.636008, 0.112092, 0.605205 },
    { 0.679160, 0.151848, 0.575189 }, { 0.714883, 0.187299, 0.546338 },
    { 0.748289, 0.222711, 0.516834 }, { 0.779604, 0.258078, 0.487539 },
    { 0.812612, 0.297928, 0.455338 }, { 0.840155, 0.333580, 0.427455 },
    { 0.866078, 0.369660, 0.400126 }, { 0.890340, 0.406398, 0.373130 },
    { 0.915471, 0.448807, 0.342890 }, { 0.935630, 0.487712, 0.315952 },
    { 0.953428, 0.527960, 0.288883 }, { 0.968526, 0.569700, 0.261721 },
    { 0.980556, 0.613039, 0.234646 }, { 0.989935, 0.663787, 0.204859 },
    { 0.994103, 0.710698, 0.180097 }, { 0.993851, 0.759304, 0.159092 },
    { 0.988648, 0.809579, 0.145357 }, { 0.976265, 0.868016, 0.143351 },
    { 0.959276, 0.921407, 0.151566 }, { 0.940015, 0.975158, 0.131326 },
};

static const double inferno_stops[32][3] = {
    { 0.001462, 0.000466, 0.013866 }, { 0.013995, 0.011225, 0.071862 },
    { 0.042253, 0.028139, 0.141141 }, { 0.087411, 0.044556, 0.224813 },
    { 0.135778, 0.046856, 0.299776 }, { 0.190367, 0.039309, 0.361447 },
    { 0.244967, 0.037055, 0.400007 }, { 0.303568, 0.049396, 0.422182 },
    { 0.354032, 0.066925, 0.430906 }, { 0.403894, 0.085580, 0.433179 },
    { 0.453651, 0.103848, 0.430498 }, { 0.503493, 0.121575, 0.423356 },
    { 0.559624, 0.141346, 0.410078 }, { 0.609330, 0.159474, 0.393589 },
    { 0.658463, 0.178962, 0.372748 }, { 0.706500, 0.200728, 0.347777 },
    { 0.758422, 0.229097, 0.315266 }, { 0.801871, 0.258674, 0.283099 },
    { 0.841969, 0.292933, 0.248564 }, { 0.878001, 0.332060, 0.212268 },
    { 0.912966, 0.381636, 0.169755 }, { 0.938675, 0.430091, 0.130438 },
    { 0.959114, 0.482014, 0.089499 }, { 0.974176, 0.536780, 0.048392 },
    { 0.983779, 0.593849, 0.023770 }, { 0.987926, 0.660250, 0.051750 },
    { 0.985566, 0.720782, 0.112229 }, { 0.977497, 0.782258, 0.185923 },
    { 0.964394, 0.843848, 0.273391 }, { 0.948683, 0.910473, 0.395289 },
    { 0.951740, 0.960587, 0.524203 }, { 0.988362, 0.998364, 0.644924 },
};

static const double cividis_stops[32][3] = {
    { 0.000000, 0.135112, 0.304751 }, { 0.000000, 0.157932, 0.357521 },
    { 0.000000, 0.178802, 0.414764 }, { 0.032110, 0.201199, 0.440785 },
    { 0.110658, 0.223170, 0.435067 }, { 0.159733, 0.245221, 0.429528 },
    { 0.199764, 0.267099, 0.425497 }, { 0.239312, 0.291562, 0.423167 },
    { 0.271639, 0.313253, 0.422837 }, { 0.302169, 0.334963, 0.424213 },
    { 0.331474, 0.356744, 0.427144 }, { 0.359916, 0.378641, 0.431501 },
    { 0.391151, 0.403464, 0.438096 }, { 0.418383, 0.425733, 0.445560 },
    { 0.445148, 0.448226, 0.454885 }, { 0.471501, 0.470960, 0.466357 },
    { 0.503185, 0.496851, 0.472305 }, { 0.532829, 0.520135, 0.472401 },
    { 0.562972, 0.543741, 0.470488 }, { 0.593622, 0.567697, 0.466401 },
    { 0.628576, 0.595104, 0.459641 }, { 0.660082, 0.619904, 0.451534 },
    { 0.691971, 0.645145, 0.441491 }, { 0.724274, 0.670859, 0.429194 },
    { 0.756975, 0.697071, 0.414659 }, { 0.794298, 0.727190, 0.395016 },
    { 0.827959, 0.754553, 0.374292 }, { 0.862105, 0.782494, 0.350011 },
    { 0.896818, 0.811030, 0.320982 }, { 0.936660, 0.843841, 0.280876 },
    { 0.973114, 0.873550, 0.234677 }, { 0.995737, 0.909344, 0.217772 },
};

/* Haze — 7 anchors at non-uniform positions (LinearSegmentedColormap). */
static const double haze_stops[7][3] = {
    { 1.000, 1.000, 1.000 },  /* white  */
    { 0.954, 0.752, 0.958 },  /* pink   */
    { 0.316, 0.713, 0.938 },  /* blue   */
    { 0.260, 0.818, 0.484 },  /* green  */
    { 0.944, 0.955, 0.277 },  /* yellow */
    { 0.991, 0.544, 0.051 },  /* orange */
    { 0.914, 0.000, 0.097 },  /* red    */
};
static const double haze_pos[7] = { 0.0, 0.158, 0.438, 0.543, 0.657, 0.845, 1.0 };

/* ramp_lerp — linear interpolation over an n-stop table. When `pos` is NULL the
 * stops are treated as evenly spaced; otherwise pos[] gives each stop's
 * position in [0,1] (must be ascending, pos[0]=0, pos[n-1]=1). */
static void ramp_lerp(double t, const double (*stops)[3], const double* pos,
                      int n, double* r, double* g, double* b) {
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    int i; double f;
    if (pos) {
        i = 0;
        while (i < n - 2 && t > pos[i + 1]) i++;
        double span = pos[i + 1] - pos[i];
        f = (span > 0.0) ? (t - pos[i]) / span : 0.0;
    } else {
        double idx = t * (double)(n - 1);
        i = (int)idx;
        if (i > n - 2) i = n - 2;
        f = idx - (double)i;
    }
    *r = stops[i][0] + f * (stops[i + 1][0] - stops[i][0]);
    *g = stops[i][1] + f * (stops[i + 1][1] - stops[i][1]);
    *b = stops[i][2] + f * (stops[i + 1][2] - stops[i][2]);
}

void viridis_rgb(double t, double* r, double* g, double* b) {
    ramp_lerp(t, viridis_stops, NULL, 32, r, g, b);
}

/* The single place the system default colormap is chosen. */
void default_ramp_rgb(double t, double* r, double* g, double* b) {
    viridis_rgb(t, r, g, b);
}

/* matplotlib_family_rgb — resolve one of the perceptually-uniform maps (plus
 * Haze) by name to raw RGB. Returns 1 on a hit, 0 otherwise. Shared by both
 * named_color_ramp (Expr) and resolve_ramp_to_rgb (raw) so they cannot drift. */
static int matplotlib_family_rgb(const char* name, double t,
                                 double* r, double* g, double* b) {
    if (strcmp(name, "Viridis") == 0) { ramp_lerp(t, viridis_stops, NULL, 32, r, g, b); return 1; }
    if (strcmp(name, "Magma")   == 0) { ramp_lerp(t, magma_stops,   NULL, 32, r, g, b); return 1; }
    if (strcmp(name, "Plasma")  == 0) { ramp_lerp(t, plasma_stops,  NULL, 32, r, g, b); return 1; }
    if (strcmp(name, "Inferno") == 0) { ramp_lerp(t, inferno_stops, NULL, 32, r, g, b); return 1; }
    if (strcmp(name, "Cividis") == 0) { ramp_lerp(t, cividis_stops, NULL, 32, r, g, b); return 1; }
    if (strcmp(name, "Haze")    == 0) { ramp_lerp(t, haze_stops,    haze_pos, 7, r, g, b); return 1; }
    return 0;
}

/* named_color_ramp — resolve a ColorFunction name string and a normalised
 * parameter t ∈ [0,1] to a concrete color expression (caller owns).
 * Returns NULL when the name is not recognised.
 *
 * Recognised names:
 *   "Rainbow"               — Hue sweep (red→violet, stops short of wrap)
 *   "Temperature"/"Thermal" — thermal blue-purple→yellow ramp
 *   "CoolTones"/"Cool"      — ice-white → sky-blue → deep navy
 *   "WarmTones"/"Warm"      — cream → amber → orange → crimson
 *   "Greyscale"/"Grayscale" — white (t=0) → black (t=1)
 */
Expr* named_color_ramp(const char* name, double t) {
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    {   /* Viridis/Magma/Plasma/Inferno/Cividis/Haze — all RGBColor ramps. */
        double rv, gv, bv;
        if (matplotlib_family_rgb(name, t, &rv, &gv, &bv)) {
            Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
            return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
        }
    }

    if (strcmp(name, "Rainbow") == 0) {
        Expr* a[1] = { expr_new_real(t * 0.8) };
        return expr_new_function(expr_new_symbol(SYM_Hue), a, 1);
    }
    if (strcmp(name, "Temperature") == 0 || strcmp(name, "Thermal") == 0) {
        double rv, gv, bv; thermal_rgb(t, &rv, &gv, &bv);
        Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
        return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
    }
    if (strcmp(name, "CoolTones") == 0 || strcmp(name, "Cool") == 0) {
        double rv, gv, bv; cool_tones_rgb(t, &rv, &gv, &bv);
        Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
        return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
    }
    if (strcmp(name, "WarmTones") == 0 || strcmp(name, "Warm") == 0) {
        double rv, gv, bv; warm_tones_rgb(t, &rv, &gv, &bv);
        Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
        return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
    }
    if (strcmp(name, "Greyscale") == 0 || strcmp(name, "Grayscale") == 0
        || strcmp(name, "GrayScale") == 0 || strcmp(name, "GreyScale") == 0
        || strcmp(name, "Grey") == 0 || strcmp(name, "Gray") == 0) {
        Expr* a[1] = { expr_new_real(1.0 - t) }; /* white at t=0, black at t=1 */
        return expr_new_function(expr_new_symbol(SYM_GrayLevel), a, 1);
    }
    /* "PhaseRings" — 1-D bar variant: pure hue sweep at full brightness.
     * The modulus-ring brightness oscillation cannot be shown in 1-D; the bar
     * label range (−π to π) already conveys the phase mapping. */
    if (strcmp(name, "PhaseRings") == 0) {
        Expr* a[1] = { expr_new_real(t) };
        return expr_new_function(expr_new_symbol(SYM_Hue), a, 1);
    }
    return NULL;
}

/* phase_rings_rgb — domain-colouring ramp for complex functions.
 *
 * Hue = Arg(re + i·im) / (2π)  maps phase continuously around the color wheel.
 * Value (brightness) = 0.1 + 0.9 · (1 + cos(2π · log|w|)) / 2  creates one
 * bright/dark ring per e-fold of |w|, compressing near poles (|w|→∞) and zeros
 * (|w|→0) to make them instantly visible as concentric ring clusters.
 *
 * Used by ComplexPlot when ColorFunction → "PhaseRings"; the color bar path
 * (which only has a 1-D t parameter) uses hue_to_rgb instead, showing the phase
 * sweep at full brightness so the bar tick labels (−π to π) are legible. */
void phase_rings_rgb(double re, double im, double* r, double* g, double* b) {
    double arg   = atan2(im, re);
    double hue   = (arg + M_PI) / (2.0 * M_PI);   /* phase → [0, 1] */
    double mod   = sqrt(re * re + im * im);
    double bright;
    if (mod < 1e-300) {
        bright = 0.0;
    } else {
        double lmod = log(mod);
        bright = (1.0 + cos(2.0 * M_PI * lmod)) * 0.5;
        bright = 0.1 + 0.9 * bright;               /* floor at 0.1 */
    }
    /* HSV → RGB: S = 1, V = bright */
    double h6 = hue * 6.0;
    h6 -= floor(h6 / 6.0) * 6.0;                  /* wrap to [0, 6) */
    int    sec = (int)h6;
    double f   = h6 - (double)sec;
    double q   = bright * (1.0 - f);
    double tv  = bright * f;
    switch (sec % 6) {
        case 0: *r = bright; *g = tv;    *b = 0.0;   break;
        case 1: *r = q;      *g = bright; *b = 0.0;  break;
        case 2: *r = 0.0;    *g = bright; *b = tv;   break;
        case 3: *r = 0.0;    *g = q;     *b = bright; break;
        case 4: *r = tv;     *g = 0.0;   *b = bright; break;
        default: *r = bright; *g = 0.0;  *b = q;     break;
    }
}

/* hue_to_rgb — HSV → RGB with s=v=1 (pure saturated hue sweep). */
static void hue_to_rgb(double h, double* r, double* g, double* b) {
    h = h - floor(h);
    double hh = h * 6.0;
    int    i  = (int)floor(hh);
    double f  = hh - (double)i;
    double q  = 1.0 - f;
    switch (((i % 6) + 6) % 6) {
        case 0: *r = 1.0; *g = f;   *b = 0.0; break;
        case 1: *r = q;   *g = 1.0; *b = 0.0; break;
        case 2: *r = 0.0; *g = 1.0; *b = f;   break;
        case 3: *r = 0.0; *g = q;   *b = 1.0; break;
        case 4: *r = f;   *g = 0.0; *b = 1.0; break;
        default: *r = 1.0; *g = 0.0; *b = q;  break;
    }
}

/* resolve_ramp_to_rgb — same lookup table as named_color_ramp but writes raw
 * RGB doubles instead of constructing an Expr.  Returns 1 on success, 0 if
 * the name is not recognised. */
int resolve_ramp_to_rgb(const char* name, double t, double* r, double* g, double* b) {
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    if (matplotlib_family_rgb(name, t, r, g, b)) return 1;
    if (strcmp(name, "Rainbow") == 0) {
        hue_to_rgb(t * 0.8, r, g, b); return 1;
    }
    if (strcmp(name, "Temperature") == 0 || strcmp(name, "Thermal") == 0) {
        thermal_rgb(t, r, g, b); return 1;
    }
    if (strcmp(name, "CoolTones") == 0 || strcmp(name, "Cool") == 0) {
        cool_tones_rgb(t, r, g, b); return 1;
    }
    if (strcmp(name, "WarmTones") == 0 || strcmp(name, "Warm") == 0) {
        warm_tones_rgb(t, r, g, b); return 1;
    }
    if (strcmp(name, "Greyscale") == 0 || strcmp(name, "Grayscale") == 0
        || strcmp(name, "GrayScale") == 0 || strcmp(name, "GreyScale") == 0
        || strcmp(name, "Grey") == 0 || strcmp(name, "Gray") == 0) {
        *r = *g = *b = 1.0 - t; return 1;
    }
    if (strcmp(name, "PhaseRings") == 0) {
        hue_to_rgb(t, r, g, b); return 1;
    }
    return 0;
}

/* ---------------------------------------------------------------------- */
/* Axis scaling helpers (ScalingFunctions)                                 */
/* ---------------------------------------------------------------------- */

double scale_apply(ScaleFnType sf, double x) {
    switch (sf) {
        case SF_LOG:     return (x > 0.0) ? log(x)   : -1e30;
        case SF_LOG2:    return (x > 0.0) ? log2(x)  : -1e30;
        case SF_LOG10:   return (x > 0.0) ? log10(x) : -1e30;
        case SF_REVERSE: return -x;
        default:         return x;
    }
}

double scale_invert(ScaleFnType sf, double w) {
    switch (sf) {
        case SF_LOG:     return exp(w);
        case SF_LOG2:    return pow(2.0, w);
        case SF_LOG10:   return pow(10.0, w);
        case SF_REVERSE: return -w;
        default:         return w;
    }
}

ScaleFnType parse_scale_fn(Expr* e) {
    if (!e) return SF_NONE;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == SYM_None || e->data.symbol.name == SYM_Automatic) return SF_NONE;
    }
    if (e->type == EXPR_STRING) {
        const char* s = e->data.string;
        if (strcmp(s, "Log")     == 0) return SF_LOG;
        if (strcmp(s, "Log2")    == 0) return SF_LOG2;
        if (strcmp(s, "Log10")   == 0) return SF_LOG10;
        if (strcmp(s, "Reverse") == 0 || strcmp(s, "Reversed") == 0) return SF_REVERSE;
    }
    return SF_NONE;
}

void parse_scaling_functions(Expr* rhs, ScaleFnType* sf_x, ScaleFnType* sf_y) {
    *sf_x = SF_NONE; *sf_y = SF_NONE;
    if (!rhs) return;
    /* {sfx, sfy} two-element List */
    if (rhs->type == EXPR_FUNCTION
        && rhs->data.function.head->type == EXPR_SYMBOL
        && rhs->data.function.head->data.symbol.name == SYM_List
        && rhs->data.function.arg_count == 2) {
        *sf_x = parse_scale_fn(rhs->data.function.args[0]);
        *sf_y = parse_scale_fn(rhs->data.function.args[1]);
    } else {
        /* Single spec applies to both axes */
        *sf_x = *sf_y = parse_scale_fn(rhs);
    }
}

void emit_scaling_meta(ScaleFnType sf_x, ScaleFnType sf_y,
                       Expr*** pt, size_t* pt_n) {
    if (sf_x == SF_NONE && sf_y == SF_NONE) return;
    *pt = realloc(*pt, sizeof(Expr*) * (*pt_n + 1));
    Expr* sm_args[2] = { expr_new_integer((int64_t)sf_x),
                         expr_new_integer((int64_t)sf_y) };
    (*pt)[(*pt_n)++] = expr_new_function(expr_new_symbol(SYM_ScalingMeta), sm_args, 2);
}
