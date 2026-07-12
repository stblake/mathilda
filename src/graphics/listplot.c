/* listplot.c — ListPlot[data, opts...].
 *
 * Unlike Plot (which is HoldAll because its function body must stay symbolic
 * while x is unbound), ListPlot's data is concrete and must be evaluated
 * (so ListPlot[Table[i^2, {i, 5}]] / ListPlot[Range[10]] work), so ListPlot
 * is a plain protected builtin. Its arguments — the data and the option
 * values — therefore arrive already evaluated (named colours like Red are
 * RGBColor[...] by the time we see them, exactly as a bare Graphics[]'s own
 * arguments would be).
 *
 * The work is purely constructive: classify the data into one or more
 * datasets, turn each into Point[...] (or Line[...] when Joined) primitives,
 * and wrap them in a Graphics[...] carrying the passthrough options. The
 * existing renderer (render.c) interprets PlotRange/Axes/AspectRatio/Frame/
 * PlotStyle/PlotLegends, so ListPlot inherits all of them for free. Sampler
 * helpers numericize_bound/palette_color are shared from plot.c (see plot.h). */

#include "listplot.h"
#include "plot.h"
#include "plot_common.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- small structural predicates (kept local; trivial) ---- */

/* is_rule_arg provided by plot_common.h */

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List;
}

/* Coerce an element to a double. Fast path for literal numbers; falls back
 * to N[] (via numericize_bound) for symbolic-but-numeric values like Pi. */
static bool to_double(Expr* e, double* v) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *v = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *v = e->data.real; return true; }
    return numericize_bound(e, v);
}

/* True if `e` is a 2-element List of numericizable values, returning them. */
static bool point2(Expr* e, double* x, double* y) {
    return is_list(e) && e->data.function.arg_count == 2
        && to_double(e->data.function.args[0], x)
        && to_double(e->data.function.args[1], y);
}

/* ---- data model ---- */

typedef struct { double* xs; double* ys; size_t n; } PointSet;

static void pointset_free(PointSet* p) { free(p->xs); free(p->ys); p->xs = p->ys = NULL; p->n = 0; }

/* DataRange option, resolved. `force_multi` (DataRange -> All) makes a flat
 * list of pairs be read as several height-datasets rather than as points. */
typedef struct { bool has_range; double xmin, xmax; bool force_multi; } DataRange;

/* Turn one dataset list into a PointSet. A list whose every element is a
 * numeric 2-pair is read as explicit points; otherwise it is a list of
 * heights {y_i}, plotted at x running over the DataRange (default 1..n),
 * with non-numeric (missing) entries skipped but still occupying their
 * index slot. */
static PointSet extract_dataset(Expr* ds, const DataRange* dr) {
    PointSet ps = { NULL, NULL, 0 };
    if (!is_list(ds)) return ps;
    size_t n = ds->data.function.arg_count;
    if (n == 0) return ps;
    Expr** el = ds->data.function.args;

    bool pairs = true;
    for (size_t i = 0; i < n; i++) {
        double a, b;
        if (!point2(el[i], &a, &b)) { pairs = false; break; }
    }

    ps.xs = malloc(sizeof(double) * n);
    ps.ys = malloc(sizeof(double) * n);
    if (pairs) {
        for (size_t i = 0; i < n; i++) {
            double a, b;
            point2(el[i], &a, &b);
            ps.xs[ps.n] = a; ps.ys[ps.n] = b; ps.n++;
        }
    } else {
        double xmin = dr->has_range ? dr->xmin : 1.0;
        double xmax = dr->has_range ? dr->xmax : (double)n;
        for (size_t i = 0; i < n; i++) {
            double y;
            if (!to_double(el[i], &y)) continue; /* missing */
            double x = (n == 1) ? xmin : xmin + (double)i * (xmax - xmin) / (double)(n - 1);
            ps.xs[ps.n] = x; ps.ys[ps.n] = y; ps.n++;
        }
    }
    if (ps.n == 0) pointset_free(&ps);
    return ps;
}

/* Classify the top-level data into datasets. Returns a malloc'd array of
 * borrowed dataset Expr* (caller frees the array, not the Exprs), or NULL on
 * nothing classifiable. *nsets is the count. The single-dataset forms
 * ({y...} heights and {{x,y}...} points, unless DataRange -> All) return the
 * data list itself as the one dataset. */
static Expr** classify(Expr* data, const DataRange* dr, size_t* nsets) {
    *nsets = 0;
    if (!is_list(data)) return NULL;
    size_t nd = data->data.function.arg_count;
    if (nd == 0) return NULL;
    Expr** el = data->data.function.args;

    /* A list with no sublists is a flat list of heights (non-numeric entries
     * are missing); a list whose every element is a numeric 2-pair is a list
     * of points; anything else (some elements are sublists) is several
     * datasets. DataRange -> All forces the multi reading of a pair list. */
    bool has_list = false;
    for (size_t i = 0; i < nd; i++)
        if (is_list(el[i])) { has_list = true; break; }

    bool all_pairs = !dr->force_multi;
    if (all_pairs) {
        for (size_t i = 0; i < nd; i++) {
            double a, b;
            if (!point2(el[i], &a, &b)) { all_pairs = false; break; }
        }
    }

    Expr** sets;
    if (all_pairs || !has_list) {
        sets = malloc(sizeof(Expr*));
        sets[0] = data;
        *nsets = 1;
    } else {
        sets = malloc(sizeof(Expr*) * nd);
        for (size_t i = 0; i < nd; i++)
            if (is_list(el[i])) sets[(*nsets)++] = el[i]; /* non-list elements are missing */
        if (*nsets == 0) { free(sets); return NULL; }
    }
    return sets;
}

/* ---- option splitting (mirrors plot.c's split_options) ---- */

typedef struct {
    bool   joined;          /* Joined -> True */
    Expr*  filling;         /* borrowed; held (Axis/Bottom/Top/number); NULL = None */
    Expr*  filling_style;   /* borrowed; held colour; NULL = default */
    bool   have_markers;    /* PlotMarkers given (accepted; glyphs not yet drawn) */
    DataRange dr;
    ScaleFnType sf_x, sf_y; /* ScalingFunctions per-axis */
} ListPlotOpts;

static double stem_baseline(Expr* filling, double ymin, double ymax) {
    if (!filling) return 0.0;
    if (filling->type == EXPR_SYMBOL) {
        if (filling->data.symbol == SYM_Axis)   return 0.0;
        if (filling->data.symbol == SYM_Bottom) return ymin;
        if (filling->data.symbol == SYM_Top)    return ymax;
    }
    double v;
    if (numericize_bound(filling, &v)) return v;
    return 0.0;
}

/* Split res's trailing Rule args into ListPlot's consumed options and a
 * passthrough Rule list bound for the Graphics[...] result. *legends_out
 * receives the (owned, evaluated) PlotLegends value or NULL; *single_color_out
 * a freshly owned copy of the resolved PlotStyle colour (default blue) used to
 * restore the curve colour after a Filling block on a single-curve plot.
 * Returns false on a malformed trailing argument. */
static bool split_options(Expr* res, ListPlotOpts* o,
                          Expr*** passthrough_out, size_t* passthrough_count_out,
                          Expr** single_color_out, Expr** legends_out) {
    o->joined = false;
    o->filling = NULL;
    o->filling_style = NULL;
    o->have_markers = false;
    o->dr.has_range = false;
    o->dr.xmin = o->dr.xmax = 0.0;
    o->dr.force_multi = false;
    o->sf_x = SF_NONE;
    o->sf_y = SF_NONE;
    *single_color_out = NULL;
    *legends_out = NULL;

    size_t argc = res->data.function.arg_count;
    /* +3 headroom for injected Axes/AspectRatio/PlotStyle defaults. */
    size_t cap = (argc > 1 ? argc - 1 : 0) + 3;
    Expr** passthrough = malloc(sizeof(Expr*) * cap);
    size_t n = 0;
    bool have_axes = false, have_aspect = false, have_style = false, have_frame = false;

#define FAIL_CLEANUP() do { \
        for (size_t k = 0; k < n; k++) expr_free(passthrough[k]); \
        free(passthrough); \
        expr_free(*single_color_out); *single_color_out = NULL; \
        expr_free(*legends_out); *legends_out = NULL; \
        return false; \
    } while (0)

    for (size_t i = 1; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) FAIL_CLEANUP();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol : NULL;

        if (name == SYM_Joined) {
            o->joined = (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_True);
        } else if (name == SYM_Filling) {
            if (!(rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_None)) o->filling = rhs;
        } else if (name == SYM_FillingStyle) {
            if (!(rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic)) o->filling_style = rhs;
        } else if (name == SYM_PlotMarkers) {
            o->have_markers = !(rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_None);
        } else if (name == SYM_DataRange) {
            if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_All) {
                o->dr.force_multi = true;
            } else if (point2(rhs, &o->dr.xmin, &o->dr.xmax)) {
                o->dr.has_range = (o->dr.xmin < o->dr.xmax);
            } else if (!(rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic)) {
                FAIL_CLEANUP();
            }
        } else if (name == SYM_ScalingFunctions) {
            Expr* v = evaluate(expr_copy(rhs));
            parse_scaling_functions(v, &o->sf_x, &o->sf_y);
            expr_free(v);
        } else if (name == SYM_PlotLegends) {
            if (*legends_out) expr_free(*legends_out);
            *legends_out = expr_copy(rhs);
        } else if (name == SYM_PlotStyle) {
            have_style = true;
            if (*single_color_out) expr_free(*single_color_out);
            *single_color_out = expr_copy(rhs);
            Expr* a[2] = { expr_copy(lhs), expr_copy(*single_color_out) };
            passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        } else if (name == SYM_AspectRatio) {
            /* render.c has no evaluator: resolve here. Automatic/Full pass
             * through verbatim; anything else numericalizes to a real. */
            have_aspect = true;
            if (rhs->type == EXPR_SYMBOL
                && (rhs->data.symbol == SYM_Automatic || rhs->data.symbol == SYM_Full)) {
                passthrough[n++] = expr_copy(arg);
            } else {
                double v;
                Expr* val = (numericize_bound(rhs, &v) && v > 0)
                    ? expr_new_real(v) : expr_new_symbol(SYM_Automatic);
                Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), val };
                passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
            }
        } else {
            if (name == SYM_Axes) have_axes = true;
            else if (name == SYM_Frame) {
                if (!(rhs->type == EXPR_SYMBOL
                      && (rhs->data.symbol == SYM_False || rhs->data.symbol == SYM_None)))
                    have_frame = true;
            }
            passthrough[n++] = expr_copy(arg);
        }
    }

    /* ListPlot-specific defaults: Axes -> True, AspectRatio -> 1/GoldenRatio,
     * PlotStyle -> blue, injected only when the caller didn't specify them. */
    if (!have_axes && !have_frame) {
        Expr* a[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_True) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_aspect) {
        const double inv_phi = 2.0 / (1.0 + sqrt(5.0)); /* 1/GoldenRatio */
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_real(inv_phi) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_style) {
        Expr* rgb_args[3] = { expr_new_real(0.2), expr_new_real(0.4), expr_new_real(0.8) };
        Expr* rgb = expr_new_function(expr_new_symbol(SYM_RGBColor), rgb_args, 3);
        Expr* a[2] = { expr_new_symbol(SYM_PlotStyle), expr_copy(rgb) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        *single_color_out = rgb;
    }

    *passthrough_out = passthrough;
    *passthrough_count_out = n;
    return true;
}
#undef FAIL_CLEANUP

/* ---- primitive construction ---- */

/* Build a List[{x1,y1}, {x2,y2}, ...] from a PointSet. */
static Expr* coords_list(const PointSet* ps) {
    Expr** pts = malloc(sizeof(Expr*) * (ps->n ? ps->n : 1));
    for (size_t j = 0; j < ps->n; j++) {
        Expr* xy[2] = { expr_new_real(ps->xs[j]), expr_new_real(ps->ys[j]) };
        pts[j] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), pts, ps->n);
    free(pts);
    return list;
}

/* Append this dataset's primitives to prims[*pc]: an optional palette colour
 * (multi-dataset only), optional Filling stems, then either a Line polyline
 * (Joined) or a Point[...] cloud (markers use render.c's data-scaled default
 * point size, matching every other Point primitive). `single_color`
 * (borrowed, may be NULL) restores the curve colour after a single-curve
 * Filling block. */
static void emit_dataset(Expr** prims, size_t* pc, const PointSet* ps, size_t idx,
                         bool multi, const ListPlotOpts* o, Expr* single_color) {
    if (multi) prims[(*pc)++] = palette_color(idx);

    if (o->filling && ps->n > 0) {
        double ymin = ps->ys[0], ymax = ps->ys[0];
        for (size_t j = 1; j < ps->n; j++) {
            if (ps->ys[j] < ymin) ymin = ps->ys[j];
            if (ps->ys[j] > ymax) ymax = ps->ys[j];
        }
        double base = stem_baseline(o->filling, ymin, ymax);

        if (o->filling_style) {
            prims[(*pc)++] = expr_copy(o->filling_style);
        } else {
            Expr* op_arg[1] = { expr_new_real(0.3) };
            prims[(*pc)++] = expr_new_function(expr_new_symbol(SYM_Opacity), op_arg, 1);
        }
        if (o->joined) {
            /* The points are drawn as a connected polyline, so fill the whole
             * continuous region between that curve and the baseline (one quad
             * per segment, shared with Plot's filling), not isolated stems. */
            size_t nq = 0;
            Expr** quads = gfx_build_fill_quads(ps->xs, ps->ys, ps->n, base, &nq);
            for (size_t j = 0; j < nq; j++) prims[(*pc)++] = quads[j];
            free(quads);
        } else {
            /* Unconnected points: a vertical stem from the baseline to each. */
            for (size_t j = 0; j < ps->n; j++) {
                Expr* lo[2] = { expr_new_real(ps->xs[j]), expr_new_real(base) };
                Expr* hi[2] = { expr_new_real(ps->xs[j]), expr_new_real(ps->ys[j]) };
                Expr* seg[2] = {
                    expr_new_function(expr_new_symbol(SYM_List), lo, 2),
                    expr_new_function(expr_new_symbol(SYM_List), hi, 2),
                };
                Expr* sl = expr_new_function(expr_new_symbol(SYM_List), seg, 2);
                Expr* la[1] = { sl };
                prims[(*pc)++] = expr_new_function(expr_new_symbol(SYM_Line), la, 1);
            }
        }
        if (!o->filling_style) {
            Expr* op2[1] = { expr_new_integer(1) };
            prims[(*pc)++] = expr_new_function(expr_new_symbol(SYM_Opacity), op2, 1);
        }
        /* Restore the curve colour the markers/line should draw in. */
        if (multi) prims[(*pc)++] = palette_color(idx);
        else if (single_color) prims[(*pc)++] = expr_copy(single_color);
    }

    Expr* coords = coords_list(ps);
    Expr* ca[1] = { coords };
    prims[(*pc)++] = expr_new_function(
        expr_new_symbol(o->joined ? SYM_Line : SYM_Point), ca, 1);
}

/* Build the internal $PlotLegendData[{color,label}, ...] metadata for the
 * renderer's draw_legend(), mirroring plot.c. `legends` (the evaluated option
 * value) is an explicit {labels...} list or Automatic; Automatic labels each
 * dataset with its 1-based index. Returns NULL if legends is None/absent. */
static Expr* build_listplot_legend_meta(Expr* legends, size_t nsets, Expr* single_color) {
    if (!legends) return NULL;
    if (legends->type == EXPR_SYMBOL && legends->data.symbol == SYM_None) return NULL;

    bool explicit_list = is_list(legends);
    bool multi = (nsets > 1);

    Expr** entries = malloc(sizeof(Expr*) * (nsets ? nsets : 1));
    for (size_t i = 0; i < nsets; i++) {
        Expr* color = multi ? palette_color(i)
                     : (single_color ? expr_copy(single_color) : palette_color(0));
        Expr* label;
        if (explicit_list && i < legends->data.function.arg_count) {
            label = expr_copy(legends->data.function.args[i]);
        } else {
            char buf[32];
            snprintf(buf, sizeof(buf), "%zu", i + 1);
            label = expr_new_string(buf);
        }
        Expr* a[2] = { color, label };
        entries[i] = expr_new_function(expr_new_symbol(SYM_List), a, 2);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_PlotLegendData), entries, nsets);
    free(entries);
    return result;
}

Expr* builtin_listplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;
    Expr* data = res->data.function.args[0];
    if (!is_list(data)) return NULL;

    ListPlotOpts o;
    Expr** passthrough = NULL;
    size_t passthrough_count = 0;
    Expr* single_color = NULL;
    Expr* legends = NULL;
    if (!split_options(res, &o, &passthrough, &passthrough_count, &single_color, &legends))
        return NULL;

    size_t nsets = 0;
    Expr** sets = classify(data, &o.dr, &nsets);
    if (!sets) goto fail_opts;

    PointSet* psets = malloc(sizeof(PointSet) * nsets);
    size_t live = 0;
    for (size_t i = 0; i < nsets; i++) {
        PointSet p = extract_dataset(sets[i], &o.dr);
        if (p.n == 0) { pointset_free(&p); continue; }
        psets[live++] = p;
    }
    free(sets);

    if (live == 0) { free(psets); goto fail_opts; }

    /* Apply ScalingFunctions: transform coords to world space in-place. */
    if (o.sf_x != SF_NONE || o.sf_y != SF_NONE) {
        for (size_t i = 0; i < live; i++) {
            for (size_t j = 0; j < psets[i].n; j++) {
                psets[i].xs[j] = scale_apply(o.sf_x, psets[i].xs[j]);
                psets[i].ys[j] = scale_apply(o.sf_y, psets[i].ys[j]);
            }
        }
    }

    bool multi = (live > 1);
    /* Worst case per dataset: 1 palette colour + (filling: style + fill shapes
     * + close + restore) + 1 Point/Line. A Joined fill emits up to 2*(n-1)
     * polygons (every segment crossing the baseline splits into 2 triangles),
     * which dominates the n stems of the unjoined case. Generous. */
    size_t cap = 0;
    for (size_t i = 0; i < live; i++) cap += 2 * psets[i].n + 8;
    Expr** prims = malloc(sizeof(Expr*) * (cap ? cap : 1));
    size_t pc = 0;
    for (size_t i = 0; i < live; i++)
        emit_dataset(prims, &pc, &psets[i], i, multi, &o, single_color);

    for (size_t i = 0; i < live; i++) pointset_free(&psets[i]);
    free(psets);

    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, pc);
    free(prims);

    emit_scaling_meta(o.sf_x, o.sf_y, &passthrough, &passthrough_count);

    Expr* legend_meta = build_listplot_legend_meta(legends, live, single_color);
    expr_free(legends);
    expr_free(single_color);

    size_t gargc = 1 + passthrough_count + (legend_meta ? 1 : 0);
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = prim_list;
    for (size_t i = 0; i < passthrough_count; i++) gargs[1 + i] = passthrough[i];
    if (legend_meta) gargs[gargc - 1] = legend_meta;
    free(passthrough);

    Expr* graphics = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
    free(gargs);
    /* The REPL renders any top-level Graphics[...] result (see repl.c). */
    return graphics;

fail_opts:
    for (size_t i = 0; i < passthrough_count; i++) expr_free(passthrough[i]);
    free(passthrough);
    expr_free(single_color);
    expr_free(legends);
    return NULL;
}
