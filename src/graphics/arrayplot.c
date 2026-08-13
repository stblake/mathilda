/* arrayplot.c — ArrayPlot[array, opts...]
 *
 * Renders a 2D array of values as a grid of coloured cells (a discrete
 * heatmap, no interpolation): row 0 of the array is drawn at the top,
 * column 0 at the left, matching Mathematica's ArrayPlot orientation.
 * Each numeric cell's value is normalised to [0,1] (ColorFunctionScaling,
 * default True) and mapped through ColorFunction. The default ramp is the
 * shared "Greyscale" ramp from plot_common (white at the minimum, black at
 * the maximum) -- the same ramp DensityPlot/ComplexPlot expose by name, and
 * Mathematica's own ArrayPlot default.
 *
 * A cell that is already a colour literal (RGBColor/GrayLevel/Hue/
 * CMYKColor) is painted directly instead, and numeric and colour cells
 * freely mix within the same array -- ArrayPlot[{{1, 0, Pink}, {0, 1, Red}}]
 * calls out two cells explicitly while the rest still follow the normal
 * heatmap, and ArrayPlot[{{Red, Blue}, {Blue, Red}}] (every cell a colour)
 * lets ArrayPlot double as a raw pixel-grid renderer.
 *
 * Unlike the HoldAll function-sampling plots (DensityPlot, Plot, ...),
 * ArrayPlot is not HoldAll: its argument is an ordinary already-evaluated
 * nested List or NDArray. An NDArray input (a dense numeric buffer by
 * construction) goes through na_load_matrix (src/linalg/numarray.h); a
 * nested List goes through arrayplot_load below, which classifies each
 * cell individually. */

#include "arrayplot.h"
#include "plot_common.h"
#include "numarray.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>
#include <math.h>

/* evaluate() borrows its argument (it copies internally rather than
 * consuming what it's handed), so `evaluate(expr_copy(x))` with the copy
 * inlined leaks that temporary -- the same bug class fixed for
 * guard_all_true in src/match.c. Every option value evaluated below routes
 * through this helper so the input copy is always freed. */
static Expr* eval_borrowed(Expr* x) {
    Expr* copy   = expr_copy(x);
    Expr* result = evaluate(copy);
    expr_free(copy);
    return result;
}

/* ---- options ---- */

typedef struct {
    Expr* color_function;          /* borrowed; NULL = default greyscale ramp */
    bool  color_function_scaling;  /* default true */
    bool  mesh;                    /* draw grid lines between cells */
    bool  show_legend;
    Expr* color_rules;             /* owned (evaluated); NULL = none */
} ArrayPlotOpts;

static bool split_arrayplot_options(Expr* res, ArrayPlotOpts* o,
                                     Expr*** pt_out, size_t* pt_n_out) {
    o->color_function         = NULL;
    o->color_function_scaling = true;
    o->mesh                   = false;
    o->show_legend            = false;
    o->color_rules            = NULL;

    size_t argc = res->data.function.arg_count;
    size_t cap  = (argc > 1 ? argc - 1 : 0) + 4;
    Expr** pt   = malloc(sizeof(Expr*) * cap);
    size_t n    = 0;
    bool have_axes = false, have_frame = false;

#define AP_FAIL() do { free(pt); if (o->color_rules) expr_free(o->color_rules); return false; } while(0)

    for (size_t i = 1; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) AP_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol.name : NULL;

        if (name == SYM_ColorFunction) {
            o->color_function = rhs;          /* borrowed */
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = eval_borrowed(rhs);
            o->color_function_scaling = !(v->type == EXPR_SYMBOL
                                          && v->data.symbol.name == SYM_False);
            expr_free(v);
        } else if (name == SYM_Mesh) {
            Expr* v = eval_borrowed(rhs);
            o->mesh = (v->type == EXPR_SYMBOL
                      && (v->data.symbol.name == SYM_All || v->data.symbol.name == SYM_True));
            expr_free(v);
        } else if (name == SYM_PlotLegends) {
            Expr* v = eval_borrowed(rhs);
            o->show_legend = !(v->type == EXPR_SYMBOL
                               && (v->data.symbol.name == SYM_None
                                   || v->data.symbol.name == SYM_False));
            expr_free(v);
        } else if (name == SYM_ColorRules) {
            Expr* v = eval_borrowed(rhs);
            if (v->type == EXPR_SYMBOL
                && (v->data.symbol.name == SYM_None || v->data.symbol.name == SYM_Automatic)) {
                expr_free(v);
            } else {
                o->color_rules = v;   /* owned; freed by builtin_arrayplot */
            }
        } else {
            if (name == SYM_Axes) have_axes = true;
            else if (name == SYM_Frame) {
                if (!(rhs->type == EXPR_SYMBOL
                      && (rhs->data.symbol.name == SYM_False
                          || rhs->data.symbol.name == SYM_None)))
                    have_frame = true;
            }
            Expr* val  = eval_borrowed(rhs);
            Expr* a[2] = { expr_copy(lhs), val };
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

    *pt_out   = pt;
    *pt_n_out = n;
    return true;
#undef AP_FAIL
}

/* ---- colour resolution ---- */

static bool is_color_head_ap(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return h == SYM_RGBColor || h == SYM_GrayLevel || h == SYM_Hue || h == SYM_CMYKColor;
}

/* t is already in [0,1]. Custom ColorFunction: try f[t]; named ramp string:
 * look up in plot_common's shared table; no ColorFunction: the "Greyscale"
 * ramp (white at t=0, black at t=1) -- Mathematica's own ArrayPlot default. */
static Expr* ap_color(Expr* cfn, double t) {
    if (!cfn) return named_color_ramp("Greyscale", t);
    if (cfn->type == EXPR_STRING) {
        Expr* c = named_color_ramp(cfn->data.string, t);
        if (c) return c;
    }
    Expr* a1[1] = { expr_new_real(t) };
    Expr* call  = expr_new_function(expr_copy(cfn), a1, 1);
    Expr* r     = evaluate(call);
    expr_free(call);
    if (is_color_head_ap(r)) return r;
    expr_free(r);
    Expr* ga[1] = { expr_new_real(0.5) };
    return expr_new_function(expr_new_symbol(SYM_GrayLevel), ga, 1);
}

/* ColorRules -> {key1 -> color1, ...} (or a bare Rule): an explicit colour
 * for particular array values, checked before the continuous ColorFunction
 * ramp. `rules` is already evaluated (see split_arrayplot_options), so keys
 * are literal numbers and colours are resolved literals (e.g. Pink's
 * OwnValue has already become RGBColor[...]). Returns an owned colour Expr
 * on an exact match, NULL if no rule's key equals v -- the cell then falls
 * through to the normal scaled ColorFunction, matching Mathematica's
 * ArrayPlot: ColorRules overrides only the values it names. */
static Expr* color_rules_lookup(const Expr* rules, double v) {
    if (!rules) return NULL;
    bool is_list = (rules->type == EXPR_FUNCTION && rules->data.function.head->type == EXPR_SYMBOL
                    && rules->data.function.head->data.symbol.name == SYM_List);
    size_t n = is_list ? rules->data.function.arg_count : 1;
    for (size_t i = 0; i < n; i++) {
        const Expr* rule = is_list ? rules->data.function.args[i] : rules;
        if (rule->type != EXPR_FUNCTION || rule->data.function.head->type != EXPR_SYMBOL
            || rule->data.function.arg_count != 2) continue;
        const char* h = rule->data.function.head->data.symbol.name;
        if (h != SYM_Rule && h != SYM_RuleDelayed) continue;
        double kd;
        if (!expr_to_real_double(rule->data.function.args[0], &kd)) continue;
        if (kd == v) return expr_copy(rule->data.function.args[1]);
    }
    return NULL;
}

/* Walks `data` as a rectangular List-of-Lists, classifying each cell as
 * either a colour literal (RGBColor/GrayLevel/Hue/CMYKColor -- painted
 * directly) or a plain number (coloured via ColorRules/ColorFunction as
 * usual). On success, rows and cols (via the out-parameters) hold the
 * dimensions and cell_colors/buf (each rows*cols, malloc'd, caller-owned)
 * hold one or the other per cell:
 * (*cell_colors)[k] is a borrowed pointer to the literal when the cell is a
 * colour, NULL when it's numeric (in which case (*buf)[k] holds the value,
 * otherwise unspecified). Returns false (nothing allocated) for a ragged or
 * non-List-of-Lists shape, or a cell that is neither a colour literal nor a
 * plain number. */
static bool arrayplot_load(const Expr* data, size_t* rows, size_t* cols,
                           const Expr*** cell_colors, double** buf) {
    if (data->type != EXPR_FUNCTION || data->data.function.head->type != EXPR_SYMBOL
        || data->data.function.head->data.symbol.name != SYM_List) return false;
    size_t nrows = data->data.function.arg_count;
    if (nrows == 0) return false;
    const Expr* row0 = data->data.function.args[0];
    if (row0->type != EXPR_FUNCTION || row0->data.function.head->type != EXPR_SYMBOL
        || row0->data.function.head->data.symbol.name != SYM_List) return false;
    size_t ncols = row0->data.function.arg_count;
    if (ncols == 0) return false;

    const Expr** colors = malloc(sizeof(Expr*) * nrows * ncols);
    double*      values = malloc(sizeof(double) * nrows * ncols);

    for (size_t i = 0; i < nrows; i++) {
        const Expr* row = data->data.function.args[i];
        if (row->type != EXPR_FUNCTION || row->data.function.head->type != EXPR_SYMBOL
            || row->data.function.head->data.symbol.name != SYM_List
            || row->data.function.arg_count != ncols) {
            free(colors); free(values); return false;
        }
        for (size_t j = 0; j < ncols; j++) {
            const Expr* cell = row->data.function.args[j];
            size_t k = i * ncols + j;
            if (is_color_head_ap(cell)) {
                colors[k] = cell; values[k] = 0.0;
            } else if (expr_to_real_double(cell, &values[k])) {
                colors[k] = NULL;
            } else {
                free(colors); free(values); return false;
            }
        }
    }
    *rows = nrows; *cols = ncols;
    *cell_colors = colors; *buf = values;
    return true;
}

/* ---- main builtin ---- */

Expr* builtin_arrayplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    Expr* data = res->data.function.args[0];

    size_t rows = 0, cols = 0;
    const Expr** cell_colors = NULL;   /* rows*cols; borrowed; NULL entry = numeric cell */
    double* buf = NULL;                /* rows*cols; valid where cell_colors[k] == NULL */

    if (data->type == EXPR_NDARRAY) {
        /* A dense numeric buffer by construction -- no cell can be a colour
         * literal, so every cell is numeric. */
        int nrows, ncols;
        if (!na_load_matrix(data, false, false, &nrows, &ncols, &buf)) return NULL;
        if (nrows < 1 || ncols < 1) { free(buf); return NULL; }
        rows = (size_t)nrows; cols = (size_t)ncols;
        cell_colors = calloc(rows * cols, sizeof(Expr*));
    } else {
        if (!arrayplot_load(data, &rows, &cols, &cell_colors, &buf)) return NULL;
    }

    ArrayPlotOpts opts;
    Expr** pt = NULL; size_t pt_n = 0;
    if (!split_arrayplot_options(res, &opts, &pt, &pt_n)) {
        free(buf); free((void*)cell_colors); return NULL;
    }

    double zmin =  1e300, zmax = -1e300;
    bool any_numeric = false;
    for (size_t k = 0; k < rows * cols; k++) {
        if (cell_colors[k]) continue;
        double v = buf[k];
        if (isfinite(v)) { any_numeric = true; if (v < zmin) zmin = v; if (v > zmax) zmax = v; }
    }
    double zspan = (zmax > zmin) ? (zmax - zmin) : 1.0;

    /* 2 primitives per cell (colour directive + Rectangle), plus mesh lines. */
    size_t cap   = rows * cols * 2 + (opts.mesh ? (rows + cols + 3) : 0) + 4;
    Expr** prims = malloc(sizeof(Expr*) * cap);
    size_t np    = 0;

    for (size_t i = 0; i < rows; i++) {
        double y0 = (double)(rows - 1 - i);   /* row 0 at the top */
        for (size_t j = 0; j < cols; j++) {
            double x0 = (double)j;
            size_t k = i * cols + j;
            Expr* color;
            if (cell_colors[k]) {
                color = expr_copy((Expr*)cell_colors[k]);
            } else {
                double v = buf[k];
                if (!isfinite(v)) continue;
                color = color_rules_lookup(opts.color_rules, v);
                if (!color) {
                    double t = v;
                    if (opts.color_function_scaling) {
                        t = (v - zmin) / zspan;
                        if (t < 0.0) t = 0.0;
                        if (t > 1.0) t = 1.0;
                    }
                    color = ap_color(opts.color_function, t);
                }
            }
            prims[np++] = color;

            Expr* p1[2] = { expr_new_real(x0),        expr_new_real(y0) };
            Expr* p2[2] = { expr_new_real(x0 + 1.0),   expr_new_real(y0 + 1.0) };
            Expr* ra[2] = { expr_new_function(expr_new_symbol(SYM_List), p1, 2),
                             expr_new_function(expr_new_symbol(SYM_List), p2, 2) };
            prims[np++] = expr_new_function(expr_new_symbol(SYM_Rectangle), ra, 2);
        }
    }

    free(buf);
    free((void*)cell_colors);
    if (opts.color_rules) expr_free(opts.color_rules);

    if (opts.mesh) {
        /* Same mid-grey already used for GridLines (render.c's default
         * grid_color ~ GrayLevel[0.82]) sits too close to white for a
         * cell-boundary line; 0.5 is the shared neutral-grey used
         * elsewhere in this module for a ColorFunction fallback. */
        Expr* mc[1] = { expr_new_real(0.5) };
        prims[np++] = expr_new_function(expr_new_symbol(SYM_GrayLevel), mc, 1);
        for (size_t i = 0; i <= rows; i++) {
            Expr* p1[2]  = { expr_new_real(0.0),        expr_new_real((double)i) };
            Expr* p2[2]  = { expr_new_real((double)cols), expr_new_real((double)i) };
            Expr* pts[2] = { expr_new_function(expr_new_symbol(SYM_List), p1, 2),
                              expr_new_function(expr_new_symbol(SYM_List), p2, 2) };
            Expr* la[1]  = { expr_new_function(expr_new_symbol(SYM_List), pts, 2) };
            prims[np++]  = expr_new_function(expr_new_symbol(SYM_Line), la, 1);
        }
        for (size_t j = 0; j <= cols; j++) {
            Expr* p1[2]  = { expr_new_real((double)j), expr_new_real(0.0) };
            Expr* p2[2]  = { expr_new_real((double)j), expr_new_real((double)rows) };
            Expr* pts[2] = { expr_new_function(expr_new_symbol(SYM_List), p1, 2),
                              expr_new_function(expr_new_symbol(SYM_List), p2, 2) };
            Expr* la[1]  = { expr_new_function(expr_new_symbol(SYM_List), pts, 2) };
            prims[np++]  = expr_new_function(expr_new_symbol(SYM_Line), la, 1);
        }
    }

    /* AspectRatio: match the array's row/column ratio unless overridden, so
     * cells render square rather than Graphics' default 1/GoldenRatio. */
    {
        bool have_aspect = false;
        for (size_t i = 0; i < pt_n; i++) {
            const Expr* e = pt[i];
            if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
                && e->data.function.args[0]->type == EXPR_SYMBOL
                && e->data.function.args[0]->data.symbol.name == SYM_AspectRatio)
                { have_aspect = true; break; }
        }
        if (!have_aspect) {
            pt = realloc(pt, sizeof(Expr*) * (pt_n + 1));
            Expr* a[2] = { expr_new_symbol(SYM_AspectRatio),
                           expr_new_real((double)rows / (double)cols) };
            pt[pt_n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    /* PlotRange: exact array extent (renderer uses it for axis bounds). */
    {
        bool have_pr = false;
        for (size_t i = 0; i < pt_n; i++) {
            const Expr* e = pt[i];
            if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
                && e->data.function.args[0]->type == EXPR_SYMBOL
                && e->data.function.args[0]->data.symbol.name == SYM_PlotRange)
                { have_pr = true; break; }
        }
        if (!have_pr) {
            pt = realloc(pt, sizeof(Expr*) * (pt_n + 1));
            Expr* xr[2] = { expr_new_real(0.0), expr_new_real((double)cols) };
            Expr* yr[2] = { expr_new_real(0.0), expr_new_real((double)rows) };
            Expr* xrng  = expr_new_function(expr_new_symbol(SYM_List), xr, 2);
            Expr* yrng  = expr_new_function(expr_new_symbol(SYM_List), yr, 2);
            Expr* rl[2] = { xrng, yrng };
            Expr* pr    = expr_new_function(expr_new_symbol(SYM_List), rl, 2);
            Expr* ra[2] = { expr_new_symbol(SYM_PlotRange), pr };
            pt[pt_n++]  = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
    }

    /* PlotLegends -> Automatic: attach $StreamColorBar -- only meaningful
     * when at least one cell is numeric (an all-colour array has no scalar
     * scale to show). */
    if (opts.show_legend && any_numeric) {
        double bar_lo = isfinite(zmin) ? zmin : 0.0;
        double bar_hi = isfinite(zmax) ? zmax : 1.0;
        if (bar_lo == bar_hi) bar_hi = bar_lo + 1.0;
        pt = realloc(pt, sizeof(Expr*) * (pt_n + 1));
        Expr* cfn_copy = opts.color_function
                         ? expr_copy(opts.color_function)
                         : expr_new_string("Greyscale");
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
    return g;
}
