/* barchart.c — BarChart[data, opts...] and Histogram[data, opts...]
 *
 * BarChart renders a vertical bar chart from explicit heights.
 * Histogram bins numeric data and renders a frequency histogram.
 * Both return Graphics[...] objects auto-displayed by the REPL.
 *
 * BarChart[{v1,...,vn}, opts...]
 *   n bars at x = 1..n with heights v1..vn, coloured via ChartStyle or
 *   the default palette.
 * BarChart[{{v1,...}, {w1,...}, ...}, opts...]
 *   Multiple grouped datasets, each dataset in a distinct palette colour.
 *
 * Histogram[data, opts...]
 * Histogram[data, k, opts...]          k equal-width bins
 * Histogram[data, {step}, opts...]     bins of given width
 * Histogram[data, {min,max,step}, opts...] explicit range + width
 *
 * Options (both):
 *   ChartStyle   — color/style list cycling through bars
 *   ChartLabels  — label list for x-axis ticks
 *   BarSpacing   — gap fraction of bar width (default 0.2)
 *   PlotLabel, standard Graphics options pass through */

#include "barchart.h"
#include "plot_common.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Shared option struct                                                 */
/* ------------------------------------------------------------------ */

typedef struct {
    double  bar_spacing;     /* default layout: 0.2 of each unit slot is gap */
    bool    bar_spacing_set; /* BarSpacing -> s given: the gap is s times the bar width */
    Expr*   chart_style;     /* borrowed List[...] or single directive; NULL = palette */
    Expr*   chart_labels;    /* borrowed List[...] of label expressions; NULL = none */
} ChartOpts;

/* Parse options shared by BarChart and Histogram; extra args start at
 * `start_idx`.  Unrecognised Rule args are collected into *pt_out. */
static bool split_chart_options(Expr* res, size_t start_idx,
                                 ChartOpts* co,
                                 Expr*** pt_out, size_t* pt_n_out) {
    co->bar_spacing  = 0.2;
    co->bar_spacing_set = false;
    co->chart_style  = NULL;
    co->chart_labels = NULL;

    size_t argc = res->data.function.arg_count;
    size_t cap  = (argc > start_idx ? argc - start_idx : 0) + 4;
    Expr** pt   = malloc(sizeof(Expr*) * cap);
    size_t n    = 0;
    bool have_axes = false, have_aspect = false;

#define BC_FAIL() do { free(pt); return false; } while(0)

    for (size_t i = start_idx; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) BC_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol.name : NULL;

        if (name == SYM_BarSpacing) {
            double v;
            if (!numericize_bound(rhs, &v) || v < 0.0) BC_FAIL();
            co->bar_spacing = v;
            co->bar_spacing_set = true;
        } else if (name == SYM_ChartStyle) {
            co->chart_style = rhs;  /* borrowed */
        } else if (name == SYM_ChartLabels) {
            co->chart_labels = rhs;  /* borrowed */
        } else {
            if (name == SYM_Axes)        have_axes   = true;
            if (name == SYM_AspectRatio) have_aspect = true;
            Expr* val  = evaluate(rhs);
            Expr* a[2] = { expr_copy(lhs), val };
            pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    if (!have_axes) {
        Expr* a[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_True) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_aspect) {
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_real(0.618) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }

    *pt_out   = pt;
    *pt_n_out = n;
    return true;
#undef BC_FAIL
}

/* ------------------------------------------------------------------ */
/* Color for bar i in dataset d                                        */
/* ------------------------------------------------------------------ */

static Expr* bar_color(Expr* chart_style, size_t i) {
    if (chart_style) {
        Expr* cs = evaluate(chart_style);
        Expr* color = NULL;
        if (cs->type == EXPR_FUNCTION
            && cs->data.function.head->type == EXPR_SYMBOL
            && cs->data.function.head->data.symbol.name == SYM_List
            && cs->data.function.arg_count > 0) {
            size_t idx = i % cs->data.function.arg_count;
            color = evaluate(cs->data.function.args[idx]);
        } else {
            /* Single directive: use for all bars. */
            color = cs;
            cs = NULL;
        }
        if (cs) expr_free(cs);
        if (color) return color;
    }
    return palette_color(i);
}

/* ------------------------------------------------------------------ */
/* Build a filled Rectangle + outline from corner coords.              */
/* Emits: color_directive, Rectangle, GrayLevel[0] (outline), Line   */
/* Returns number of primitives written to prims[np..]. */
/* ------------------------------------------------------------------ */

static size_t emit_bar(Expr** prims, size_t np,
                        double x0, double y0, double x1, double y1,
                        Expr* color) {
    prims[np++] = color;  /* fill color — takes ownership */

    Expr* p1[2] = { expr_new_real(x0), expr_new_real(y0) };
    Expr* p2[2] = { expr_new_real(x1), expr_new_real(y1) };
    Expr* ra[2] = { expr_new_function(expr_new_symbol(SYM_List), p1, 2),
                     expr_new_function(expr_new_symbol(SYM_List), p2, 2) };
    prims[np++] = expr_new_function(expr_new_symbol(SYM_Rectangle), ra, 2);

    /* Thin dark outline for visual separation. */
    Expr* gl_a[1] = { expr_new_real(0.15) };
    prims[np++] = expr_new_function(expr_new_symbol(SYM_GrayLevel), gl_a, 1);

    Expr* th_a[1] = { expr_new_real(0.004) };
    prims[np++] = expr_new_function(expr_new_symbol(SYM_Thickness), th_a, 1);

    /* Outline path: bottom-left → top-left → top-right → bottom-right → close.
     * Line[{{x0,y0},...}] — single List-of-points argument, as the renderer
     * expects. */
    double xs[5] = { x0, x0, x1, x1, x0 };
    double ys[5] = { y0, y1, y1, y0, y0 };
    Expr** pts = malloc(sizeof(Expr*) * 5);
    for (int k = 0; k < 5; k++) {
        Expr* pp[2] = { expr_new_real(xs[k]), expr_new_real(ys[k]) };
        pts[k] = expr_new_function(expr_new_symbol(SYM_List), pp, 2);
    }
    Expr* pts_list = expr_new_function(expr_new_symbol(SYM_List), pts, 5);
    free(pts);
    Expr* line_args[1] = { pts_list };
    prims[np++] = expr_new_function(expr_new_symbol(SYM_Line), line_args, 1);

    return np;
}

/* ------------------------------------------------------------------ */
/* Flatten a possibly nested dataset list                              */
/* Returns false if vals is not a list of numeric values.             */
/* ------------------------------------------------------------------ */

/* keep_slots: a non-numeric entry (a symbol, Missing[], a complex, Infinity) keeps its slot as
 * NaN, so the bars after it stay at their positions (audit BarChart#1); otherwise it is dropped.
 * Exact numeric constants (Pi, E) are numericized (audit Histogram#3). */
static bool extract_reals(Expr* list, double** vals_out, size_t* n_out, bool keep_slots) {
    if (!list || list->type != EXPR_FUNCTION
        || list->data.function.head->type != EXPR_SYMBOL
        || list->data.function.head->data.symbol.name != SYM_List)
        return false;

    size_t n = list->data.function.arg_count;
    double* v = malloc(sizeof(double) * (n > 0 ? n : 1));
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* ev = evaluate(list->data.function.args[i]);
        double x;
        if ((expr_to_real_double(ev, &x) || numericize_bound(ev, &x)) && isfinite(x)) v[k++] = x;
        else if (keep_slots) v[k++] = NAN;
        expr_free(ev);
    }
    *vals_out = v;
    *n_out    = k;
    return k > 0;
}

/* ------------------------------------------------------------------ */
/* BarChart                                                             */
/* ------------------------------------------------------------------ */

Expr* builtin_barchart(Expr* res) {
    if (res->data.function.arg_count < 1) return NULL;

    Expr* data_arg = evaluate(res->data.function.args[0]);
    if (!data_arg || data_arg->type != EXPR_FUNCTION
        || data_arg->data.function.head->type != EXPR_SYMBOL
        || data_arg->data.function.head->data.symbol.name != SYM_List) {
        expr_free(data_arg);
        return NULL;
    }

    /* Determine if this is a multi-dataset call:
     * args[0] is a List of Lists (each inner list is one dataset). */
    bool multi = false;
    size_t outer_n = data_arg->data.function.arg_count;
    if (outer_n > 0) {
        Expr* first = data_arg->data.function.args[0];
        if (first->type == EXPR_FUNCTION
            && first->data.function.head->type == EXPR_SYMBOL
            && first->data.function.head->data.symbol.name == SYM_List)
            multi = true;
    }

    ChartOpts co;
    Expr** pt = NULL; size_t pt_n = 0;
    if (!split_chart_options(res, 1, &co, &pt, &pt_n)) {
        expr_free(data_arg); return NULL; }

    /* Build arrays: datasets[d] = array of heights, lengths[d] = count. */
    size_t  n_datasets = multi ? outer_n : 1;
    double** datasets  = malloc(sizeof(double*) * n_datasets);
    size_t*  lengths   = malloc(sizeof(size_t)  * n_datasets);

    if (multi) {
        for (size_t d = 0; d < n_datasets; d++) {
            if (!extract_reals(data_arg->data.function.args[d],
                                &datasets[d], &lengths[d], true)) {
                /* Fill with zeros on bad input. */
                datasets[d] = malloc(sizeof(double));
                datasets[d][0] = 0.0; lengths[d] = 1;
            }
        }
    } else {
        if (!extract_reals(data_arg, &datasets[0], &lengths[0], true)) {
            expr_free(data_arg);
            free(datasets); free(lengths);
            for (size_t i = 0; i < pt_n; i++) expr_free(pt[i]);
            free(pt);
            return NULL;
        }
    }

    expr_free(data_arg);

    /* Determine n_bars (max length across datasets). */
    size_t n_bars = 0;
    for (size_t d = 0; d < n_datasets; d++)
        if (lengths[d] > n_bars) n_bars = lengths[d];

    /* Layout. One dataset: bar i centred at x = i. Several datasets: one GROUP per dataset,
     * centred at x = d, its elements side by side and coloured by position within the group
     * (audit BarChart#3: grouping by element index transposed the reference's layout).
     * The bar width: 1 - 0.2 by default; with BarSpacing -> s the gap between groups is s times
     * the group width, so the group is 1/(1 + s) wide (audit BarChart#2: s >= 1 gave empty or
     * inverted bars). */
    size_t n_groups = multi ? n_datasets : n_bars;
    size_t per_group = multi ? n_bars : 1;
    double group_w  = co.bar_spacing_set ? 1.0 / (1.0 + co.bar_spacing) : 1.0 - co.bar_spacing;
    double sub_w    = group_w / (double)per_group;

    /* ymin for PlotRange: 0 unless negative values exist. */
    double ymin_all = 0.0, ymax_all = 0.0;
    for (size_t d = 0; d < n_datasets; d++)
        for (size_t i = 0; i < lengths[d]; i++) {
            double v = datasets[d][i];
            if (isnan(v)) continue;
            if (v < ymin_all) ymin_all = v;
            if (v > ymax_all) ymax_all = v;
        }
    if (ymax_all == ymin_all) ymax_all = ymin_all + 1.0;

    /* Capacity: at most (n_datasets * n_bars) bars × 5 prims each + labels. */
    size_t cap   = n_datasets * n_bars * 6 + n_bars * 2 + n_datasets * 2 + 4;
    Expr** prims = malloc(sizeof(Expr*) * cap);
    size_t np    = 0;

    for (size_t g = 0; g < n_groups; g++) {
        double group_left = (double)(g + 1) - group_w * 0.5;
        for (size_t j = 0; j < per_group; j++) {
            size_t d = multi ? g : 0;
            size_t i = multi ? j : g;
            if (i >= lengths[d]) continue;
            double v   = datasets[d][i];
            if (isnan(v)) continue;                     /* a missing entry: an empty slot */
            double x0  = group_left + (double)j * sub_w;
            double x1  = x0 + sub_w;
            double y_lo = (v < 0.0) ? v : 0.0;
            double y_hi = (v > 0.0) ? v : 0.0;
            Expr*  col  = bar_color(co.chart_style, i);
            np = emit_bar(prims, np, x0, y_lo, x1, y_hi, col);
        }
    }
    n_bars = n_groups;                                  /* the label / range slots below */

    /* ChartLabels: packaged as $BarChartLabels[{x1,"label1"}, {x2,"label2"}, ...]
     * so the renderer can draw them in screen space (crisp, bold, fixed size)
     * below the x-axis, matching the style of axis tick numbers. */
    if (co.chart_labels) {
        Expr* labels_ev = evaluate(co.chart_labels);
        if (labels_ev->type == EXPR_FUNCTION
            && labels_ev->data.function.head->type == EXPR_SYMBOL
            && labels_ev->data.function.head->data.symbol.name == SYM_List) {
            size_t nl = labels_ev->data.function.arg_count;
            size_t count = (nl < n_bars) ? nl : n_bars;
            Expr** pairs = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                Expr* lbl = expr_copy(labels_ev->data.function.args[i]);
                Expr* pair[2] = { expr_new_real((double)(i + 1)), lbl };
                pairs[i] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
            }
            pt = realloc(pt, sizeof(Expr*) * (pt_n + 1));
            pt[pt_n++] = expr_new_function(expr_new_symbol(SYM_BarChartLabels), pairs, count);
            free(pairs);
        }
        expr_free(labels_ev);
    }

    /* Free dataset arrays. */
    for (size_t d = 0; d < n_datasets; d++) free(datasets[d]);
    free(datasets); free(lengths);

    /* PlotRange: explicit so the renderer knows where the axes are. */
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
            double pad_y = (ymax_all - ymin_all) * 0.05 + 0.01;
            Expr* xr[2] = { expr_new_real(0.5), expr_new_real((double)n_bars + 0.5) };
            Expr* yr[2] = { expr_new_real(ymin_all - pad_y),
                            expr_new_real(ymax_all + pad_y) };
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

    /* expr_new_function copies the argument array (adopting its elements), so the array
     * itself is ours to free. */
    Expr* g = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
    free(gargs);
    return g;
}

/* ------------------------------------------------------------------ */
/* Histogram                                                           */
/* ------------------------------------------------------------------ */

/* Parse the optional bin specification at args[1] (if it exists and is
 * numeric/list).  Returns the first option-arg index. */
static size_t histogram_parse_bins(Expr* res,
                                    double data_min, double data_max,
                                    double* bin_min_out, double* bin_max_out,
                                    int*    n_bins_out, bool* explicit_out) {
    size_t argc = res->data.function.arg_count;
    *explicit_out = false;

    /* Defaults: Sturges' rule, range = [data_min, data_max]. */
    *bin_min_out = data_min;
    *bin_max_out = data_max;

    /* From the context summary: Histogram[data, opts...] or
     * Histogram[data, k, opts...] or Histogram[data, {step}, opts...] or
     * Histogram[data, {min,max,step}, opts...]. */
    if (argc < 2) goto use_default;

    Expr* spec = res->data.function.args[1];
    if (is_rule_arg(spec)) goto use_default;  /* opts start here */

    Expr* ev = evaluate(spec);
    if (ev->type == EXPR_INTEGER) {
        int k = (int)ev->data.integer;
        expr_free(ev);
        if (k >= 1) { *n_bins_out = k; return 2; }
        goto use_default;
    }
    double vd;
    if (expr_to_real_double(ev, &vd) && vd > 0.0) {
        expr_free(ev);
        /* Interpret as number of bins (Real). */
        *n_bins_out = (int)vd;
        return 2;
    }
    if (ev->type == EXPR_FUNCTION
        && ev->data.function.head->type == EXPR_SYMBOL
        && ev->data.function.head->data.symbol.name == SYM_List) {
        size_t lc = ev->data.function.arg_count;
        if (lc == 1) {
            /* {step} */
            double step;
            if (numericize_bound(ev->data.function.args[0], &step) && step > 0.0) {
                /* Bins of width exactly `step`, on multiples of it, covering the data (audit
                 * Histogram#2: the range used to be stretched to [min, max]). */
                double lo = floor(data_min / step) * step;
                int nb = (int)floor((data_max - lo) / step) + 1;
                if (nb < 1) nb = 1;
                *bin_min_out = lo;
                *bin_max_out = lo + nb * step;
                *n_bins_out = nb;
                expr_free(ev);
                return 2;
            }
        } else if (lc == 3) {
            /* {min, max, step} */
            double lo, hi, step;
            if (numericize_bound(ev->data.function.args[0], &lo)
                && numericize_bound(ev->data.function.args[1], &hi)
                && numericize_bound(ev->data.function.args[2], &step)
                && step > 0.0 && hi > lo) {
                /* width exactly `step`: the last bin reaches past hi when step does not
                 * divide the range (audit Histogram#2) */
                *explicit_out = true;
                int nb = (int)ceil((hi - lo) / step - 1e-12);
                if (nb < 1) nb = 1;
                *bin_min_out = lo;
                *bin_max_out = lo + nb * step;
                *n_bins_out = nb;
                expr_free(ev);
                return 2;
            }
        }
    }
    expr_free(ev);

use_default:;
    /* Sturges' rule: k = ceil(log2(n)) + 1, but we don't know n here; caller
     * will override n_bins_out after counting data. Signal with -1. */
    *n_bins_out = -1;
    return (argc >= 2 && !is_rule_arg(res->data.function.args[1])) ? 2 : 1;
}

Expr* builtin_histogram(Expr* res) {
    if (res->data.function.arg_count < 1) return NULL;

    Expr* data_arg = evaluate(res->data.function.args[0]);
    double* vals   = NULL;
    size_t  n_vals = 0;
    if (!extract_reals(data_arg, &vals, &n_vals, false) || n_vals == 0) {
        expr_free(data_arg); free(vals);
        return NULL;
    }
    expr_free(data_arg);

    double dmin = vals[0], dmax = vals[0];
    for (size_t i = 1; i < n_vals; i++) {
        if (vals[i] < dmin) dmin = vals[i];
        if (vals[i] > dmax) dmax = vals[i];
    }
    if (dmin == dmax) dmax = dmin + 1.0;

    double bin_min, bin_max;
    int    n_bins;
    bool   explicit_range;
    size_t opts_start = histogram_parse_bins(res, dmin, dmax,
                                              &bin_min, &bin_max, &n_bins,
                                              &explicit_range);

    /* Apply Sturges' rule if not specified. */
    if (n_bins < 1) {
        n_bins = (int)ceil(log2((double)n_vals)) + 1;
        if (n_bins < 2) n_bins = 2;
        if (n_bins > 50) n_bins = 50;
    }

    /* Bin the data. */
    int* counts = calloc((size_t)n_bins, sizeof(int));
    double bin_w = (bin_max - bin_min) / n_bins;
    for (size_t i = 0; i < n_vals; i++) {
        /* An EXPLICIT {min, max, dw} range counts only the data inside its left-closed bins
         * [min + k dw, min + (k+1) dw); a value below min, or at or above the last edge, is not in
         * any bin and is not counted. Clamping it into the edge bin (as before) made 15 count in
         * [5, 10). A derived range spans exactly [data min, data max], so there the maximum itself
         * belongs to the last bin and the clamp is only absorbing it. */
        double q = (vals[i] - bin_min) / bin_w;
        if (explicit_range && !(q >= 0.0 && q < (double)n_bins)) continue;
        int b = (int)q;
        if (b < 0) b = 0;
        if (b >= n_bins) b = n_bins - 1;
        counts[b]++;
    }
    free(vals);

    int max_count = 0;
    for (int b = 0; b < n_bins; b++)
        if (counts[b] > max_count) max_count = counts[b];

    ChartOpts co;
    Expr** pt = NULL; size_t pt_n = 0;
    if (!split_chart_options(res, opts_start, &co, &pt, &pt_n)) {
        free(counts); return NULL; }

    /* 5 prims per bar (color + rect + outline color + thickness + Line). */
    size_t cap   = (size_t)n_bins * 6 + 4;
    Expr** prims = malloc(sizeof(Expr*) * cap);
    size_t np    = 0;

    /* Histogram uses a single palette colour (palette_color(0)) by default. */
    for (int b = 0; b < n_bins; b++) {
        if (counts[b] == 0) continue;
        double x0 = bin_min + b * bin_w;
        double x1 = x0 + bin_w;
        Expr*  col = bar_color(co.chart_style, (size_t)b);
        np = emit_bar(prims, np, x0, 0.0, x1, (double)counts[b], col);
    }
    free(counts);

    /* PlotRange covering full x domain and [0, max_count * 1.05]. */
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
            Expr* xr[2] = { expr_new_real(bin_min), expr_new_real(bin_max) };
            Expr* yr[2] = { expr_new_real(0.0),
                            expr_new_real((double)max_count * 1.06 + 0.5) };
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
    free(gargs);        /* the array is copied; its elements are adopted */
    return g;
}
