/* plot.c — Plot[f, {x, xmin, xmax}, opts...].
 *
 * HoldAll, like Table/Do: f and the iterator spec must not be
 * pre-evaluated (x has no value yet). Splits trailing options into the
 * sampler's own (PlotPoints/MaxRecursion/MaxPlotPoints, consumed here)
 * and everything else (PlotRange/AspectRatio/PlotStyle/Axes/.../ImageSize,
 * copied through onto the resulting Graphics[...] unevaluated -- render.c
 * is the single place that interprets those, whether reached via Plot's
 * auto-display or a later Show[]). */

#include "plot.h"
#include "show.h"
#include "sampling.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <stdlib.h>
#include <math.h>

static bool expr_to_real_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)   { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational) {
        Expr* n = e->data.function.args[0];
        Expr* d = e->data.function.args[1];
        if (n->type == EXPR_INTEGER && d->type == EXPR_INTEGER && d->data.integer != 0) {
            *out = (double)n->data.integer / (double)d->data.integer;
            return true;
        }
    }
    return false;
}

/* Iterator bounds are frequently symbolic-but-numeric (2 Pi, E, Sqrt[2],
 * ...) rather than a literal machine number, so route them through N[]
 * (exactly as a user typing N[2 Pi] would) before extracting a double. */
static bool numericize_bound(Expr* e, double* out) {
    Expr* n_arg[1] = { expr_copy(e) };
    Expr* n_call = expr_new_function(expr_new_symbol("N"), n_arg, 1);
    Expr* result = evaluate(n_call);
    expr_free(n_call);
    bool ok = expr_to_real_double(result, out) && isfinite(*out);
    expr_free(result);
    return ok;
}

typedef struct {
    Expr* var;  /* iterator symbol, borrowed */
    Expr* body; /* f, borrowed */
} PlotEvalCtx;

static bool plot_eval_fn(double x, void* ctx_, double* y_out) {
    PlotEvalCtx* ctx = (PlotEvalCtx*)ctx_;
    Expr* xval = expr_new_real(x);
    symtab_add_own_value(ctx->var->data.symbol, ctx->var, xval);
    Expr* result = evaluate(ctx->body);
    expr_free(xval);

    double y;
    bool ok = expr_to_real_double(result, &y) && isfinite(y);
    expr_free(result);
    if (ok) *y_out = y;
    return ok;
}

/* ---- Option parsing: PlotPoints / MaxRecursion / MaxPlotPoints ---- */

typedef struct {
    long plot_points;
    int  max_recursion;
    long max_plot_points; /* <= 0 means unbounded */
} PlotSampleOpts;

static bool is_rule_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    return (h->data.symbol == SYM_Rule || h->data.symbol == SYM_RuleDelayed)
        && e->data.function.arg_count == 2;
}

static bool parse_long_value(Expr* rhs, long* out) {
    Expr* v = evaluate(expr_copy(rhs));
    bool ok = (v->type == EXPR_INTEGER);
    if (ok) *out = (long)v->data.integer;
    expr_free(v);
    return ok;
}

/* Splits res's trailing Rule args (starting at index 2) into the sampler
 * options above and a passthrough list of borrowed-then-copied Rule
 * expressions destined for the Graphics[...] result. Returns false on a
 * malformed trailing argument (not a Rule, or a known option with a
 * badly-typed value) so the caller can decline evaluation entirely. */
static bool split_options(Expr* res, PlotSampleOpts* sopts,
                           Expr*** passthrough_out, size_t* passthrough_count_out) {
    sopts->plot_points = 25;
    sopts->max_recursion = 6;
    sopts->max_plot_points = -1;

    size_t argc = res->data.function.arg_count;
    /* +3 headroom for the Axes/AspectRatio/PlotStyle defaults potentially
     * appended below when the caller didn't already supply them. */
    size_t cap = (argc > 2 ? argc - 2 : 0) + 3;
    Expr** passthrough = malloc(sizeof(Expr*) * cap);
    size_t n = 0;

    bool have_axes = false, have_aspect = false, have_style = false;

    for (size_t i = 2; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) { free(passthrough); return false; }
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol : NULL;

        if (name == SYM_PlotPoints) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 2) { free(passthrough); return false; }
            sopts->plot_points = v;
        } else if (name == SYM_MaxRecursion) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 0) { free(passthrough); return false; }
            sopts->max_recursion = (int)v;
        } else if (name == SYM_MaxPlotPoints) {
            Expr* v = evaluate(expr_copy(rhs));
            bool is_inf = (v->type == EXPR_SYMBOL && v->data.symbol == SYM_Infinity);
            long lv = -1;
            bool ok = is_inf;
            if (!ok && v->type == EXPR_INTEGER && v->data.integer > 0) { lv = (long)v->data.integer; ok = true; }
            expr_free(v);
            if (!ok) { free(passthrough); return false; }
            sopts->max_plot_points = lv;
        } else {
            if (name == SYM_Axes) have_axes = true;
            else if (name == SYM_AspectRatio) have_aspect = true;
            else if (name == SYM_PlotStyle) have_style = true;
            passthrough[n++] = expr_copy(arg);
        }
    }

    /* Plot-specific defaults (distinct from bare Graphics[]'s defaults of
     * Axes->False / AspectRatio->Automatic): inject only what the caller
     * didn't already specify. */
    if (!have_axes) {
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
        Expr* a[2] = { expr_new_symbol(SYM_PlotStyle), rgb };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }

    *passthrough_out = passthrough;
    *passthrough_count_out = n;
    return true;
}

/* Distinct, harmonious per-curve colours for multi-function plots. This is
 * Mathematica's default plot palette (ColorData[97], "the 97 colours"),
 * hand-tuned so any prefix of the list reads well together; we cycle it for
 * plots with more curves than entries. */
static Expr* palette_color(size_t i) {
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

/* Sample one function `body` over [xmin,xmax] and return its run-split
 * Line[...] primitives (a run boundary is a recorded gap/singularity).
 * Returns NULL with *out_count == 0 when nothing is plottable. The caller
 * must already have shadowed the iterator variable. */
static Expr** sample_lines(Expr* body, Expr* var, double xmin, double xmax,
                           const PlotSampleOpts* sopts, size_t* out_count) {
    *out_count = 0;
    PlotEvalCtx ctx = { .var = var, .body = body };

    size_t npts;
    PlotPoint* pts = plot_sample_adaptive(plot_eval_fn, &ctx, xmin, xmax,
                                           sopts->plot_points, sopts->max_recursion,
                                           sopts->max_plot_points, &npts);
    if (!pts || npts == 0) { plot_points_free(pts); return NULL; }

    Expr** prims = malloc(sizeof(Expr*) * npts);
    size_t prim_count = 0;
    size_t run_start = 0;
    for (size_t i = 1; i <= npts; i++) {
        bool end_of_run = (i == npts) || pts[i].break_before;
        if (end_of_run) {
            size_t run_len = i - run_start;
            if (run_len >= 2) {
                Expr** line_pts = malloc(sizeof(Expr*) * run_len);
                for (size_t j = 0; j < run_len; j++) {
                    Expr* xy[2] = { expr_new_real(pts[run_start + j].x), expr_new_real(pts[run_start + j].y) };
                    line_pts[j] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
                }
                Expr* pts_list = expr_new_function(expr_new_symbol(SYM_List), line_pts, run_len);
                free(line_pts);
                Expr* line_args[1] = { pts_list };
                prims[prim_count++] = expr_new_function(expr_new_symbol(SYM_Line), line_args, 1);
            }
            run_start = i;
        }
    }
    plot_points_free(pts);
    if (prim_count == 0) { free(prims); return NULL; }
    *out_count = prim_count;
    return prims;
}

Expr* builtin_plot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    IterSpec ispec;
    if (!iter_spec_parse(spec, &ispec)) return NULL;
    if (ispec.kind != ITER_KIND_RANGE || !ispec.var) { iter_spec_free(&ispec); return NULL; }

    double xmin, xmax;
    if (!numericize_bound(ispec.imin, &xmin) || !numericize_bound(ispec.imax, &xmax) || !(xmin < xmax)) {
        iter_spec_free(&ispec);
        return NULL;
    }

    PlotSampleOpts sopts;
    Expr** passthrough = NULL;
    size_t passthrough_count = 0;
    if (!split_options(res, &sopts, &passthrough, &passthrough_count)) {
        iter_spec_free(&ispec);
        return NULL;
    }

    /* One curve (Plot[f, …]) or several (Plot[{f1, f2, …}, …]). A List head
     * on the first argument is the multi-curve form. */
    Expr* single = f;
    Expr** bodies;
    size_t nfun;
    bool is_list = (f->type == EXPR_FUNCTION && f->data.function.head
                    && f->data.function.head->type == EXPR_SYMBOL
                    && f->data.function.head->data.symbol == SYM_List);
    if (is_list) {
        nfun = f->data.function.arg_count;
        bodies = f->data.function.args; /* borrowed */
    } else {
        nfun = 1;
        bodies = &single;               /* borrowed */
    }

    if (nfun == 0) { /* Plot[{}, …] — nothing to draw */
        iter_spec_free(&ispec);
        for (size_t i = 0; i < passthrough_count; i++) expr_free(passthrough[i]);
        free(passthrough);
        return NULL;
    }

    /* Sample every curve under a single shadow of the iterator variable. */
    Rule* old_own = iter_spec_shadow(ispec.var);
    Expr*** per = malloc(sizeof(Expr**) * nfun);
    size_t* per_count = malloc(sizeof(size_t) * nfun);
    size_t total = 0;
    bool any = false;
    for (size_t fi = 0; fi < nfun; fi++) {
        per[fi] = sample_lines(bodies[fi], ispec.var, xmin, xmax, &sopts, &per_count[fi]);
        total += per_count[fi];
        if (per_count[fi] > 0) any = true;
    }
    iter_spec_restore(ispec.var, old_own);
    iter_spec_free(&ispec);

    if (!any) {
        free(per); free(per_count); /* all per[fi] are NULL when nothing plotted */
        for (size_t i = 0; i < passthrough_count; i++) expr_free(passthrough[i]);
        free(passthrough);
        return NULL;
    }

    /* Assemble the primitive list. With more than one curve, each curve's
     * runs are introduced by a palette colour directive so the series read
     * as distinct yet harmonious; a single curve keeps the PlotStyle colour
     * that split_options supplied. */
    bool multi = (nfun > 1);
    size_t cap = total + (multi ? nfun : 0);
    Expr** prims = malloc(sizeof(Expr*) * (cap > 0 ? cap : 1));
    size_t prim_count = 0;
    for (size_t fi = 0; fi < nfun; fi++) {
        if (multi) prims[prim_count++] = palette_color(fi);
        for (size_t j = 0; j < per_count[fi]; j++) prims[prim_count++] = per[fi][j];
        free(per[fi]);
    }
    free(per); free(per_count);

    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, prim_count);
    free(prims);

    size_t gargc = 1 + passthrough_count;
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = prim_list;
    for (size_t i = 0; i < passthrough_count; i++) gargs[1 + i] = passthrough[i];
    free(passthrough);

    Expr* graphics = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
    free(gargs);

    graphics_show(graphics);
    return graphics;
}
