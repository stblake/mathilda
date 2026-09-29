/* show.c — Show[graphics, opts...] / Show[g1, g2, ..., opts...] plus the
 * no-Raylib graphics_show() stub. When USE_GRAPHICS is compiled in, render.c
 * provides the real graphics_show(); this file's stub is excluded then (see
 * the #ifndef below) so there's exactly one definition either way.
 *
 * Show follows Mathematica:
 *
 *   Show[g, opts]           g with opts overriding its own options.
 *   Show[g1, g2, ..., opts] the graphics overlaid. Each input's primitives
 *                           become one List -- a directive scope -- so a
 *                           colour or Dashing in g1 never restyles g2. An
 *                           input's own PlotStyle (the style a single-curve
 *                           Plot/ListPlot is drawn in) is baked into its
 *                           scope, because the combined object keeps only
 *                           g1's options. Options come from g1 unless opts
 *                           overrides them; PlotRange is the union of the
 *                           inputs' ranges (a graphic without an explicit
 *                           range contributes its primitives' extent), and
 *                           is left automatic when no input fixes one.
 *   Show[{g1, g2, ...}, ..] the same for (nested) lists of graphics.
 *
 * Plot's $PlotResample metadata re-samples only its own curves, so it is
 * dropped from a combination (the static primitives stay); the inputs'
 * $PlotLegendData entries are concatenated into one legend. 2D and 3D
 * graphics cannot be mixed (the call stays unevaluated). */

#include "show.h"
#include "render3d.h"
#include "sym_names.h"
#include "plot_common.h"
#include <stdlib.h>
#include <string.h>

static bool head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

static bool show_is_rule(const Expr* e) {
    return (head_is(e, SYM_Rule) || head_is(e, SYM_RuleDelayed))
        && e->data.function.arg_count == 2;
}

/* Show[] accepts either container head -- Graphics[...] (2D) or
 * Graphics3D[...] (Plot3D's output) -- and preserves whichever one it got. */
static bool is_graphics(const Expr* e) {
    return (head_is(e, SYM_Graphics) || head_is(e, SYM_Graphics3D))
        && e->data.function.arg_count >= 1;
}

/* ------------------------------------------------------ argument split --- */

typedef struct { const Expr** v; size_t n, cap; } PtrVec;

static void pv_push(PtrVec* p, const Expr* e) {
    if (p->n == p->cap) {
        p->cap = p->cap ? p->cap * 2 : 8;
        p->v = realloc(p->v, sizeof(Expr*) * p->cap);
    }
    p->v[p->n++] = e;
}

/* Collect graphics from `e`, descending into Lists. False on anything that
 * is neither a graphics object nor a List of them. */
static bool collect_graphics(const Expr* e, PtrVec* gs) {
    if (is_graphics(e)) { pv_push(gs, e); return true; }
    if (head_is(e, SYM_List)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!collect_graphics(e->data.function.args[i], gs)) return false;
        return true;
    }
    return false;
}

/* The value of option `sym` among a graphic's own options, or NULL. */
static const Expr* graphic_option(const Expr* g, const char* sym) {
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        const Expr* a = g->data.function.args[i];
        if (show_is_rule(a) && a->data.function.args[0]->type == EXPR_SYMBOL
            && a->data.function.args[0]->data.symbol.name == sym)
            return a->data.function.args[1];
    }
    return NULL;
}

static bool opts_have(const PtrVec* opts, const char* sym) {
    for (size_t i = 0; i < opts->n; i++) {
        const Expr* lhs = opts->v[i]->data.function.args[0];
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == sym) return true;
    }
    return false;
}

/* ----------------------------------------------------- one-graphic Show --- */

/* Show[g, opts]: g's options with each of opts overriding by name or
 * appending, matching Mathematica's later-wins option semantics. */
static Expr* show_single(const Expr* g, const PtrVec* opts) {
    if (opts->n == 0) return expr_copy((Expr*)g);

    size_t gargc = g->data.function.arg_count;
    size_t max_opts = (gargc - 1) + opts->n;
    Expr** out = malloc(sizeof(Expr*) * (max_opts > 0 ? max_opts : 1));
    size_t n = 0;
    for (size_t i = 1; i < gargc; i++) out[n++] = expr_copy(g->data.function.args[i]);

    for (size_t i = 0; i < opts->n; i++) {
        const Expr* newopt = opts->v[i];
        const Expr* newlhs = newopt->data.function.args[0];
        bool replaced = false;
        for (size_t j = 0; j < n; j++) {
            if (!show_is_rule(out[j])) continue;          /* metadata node */
            const Expr* lhs = out[j]->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL && newlhs->type == EXPR_SYMBOL
                && lhs->data.symbol.name == newlhs->data.symbol.name) {
                expr_free(out[j]);
                out[j] = expr_copy((Expr*)newopt);
                replaced = true;
                break;
            }
        }
        if (!replaced) out[n++] = expr_copy((Expr*)newopt);
    }

    Expr** gargs = malloc(sizeof(Expr*) * (1 + n));
    gargs[0] = expr_copy(g->data.function.args[0]);
    for (size_t i = 0; i < n; i++) gargs[1 + i] = out[i];
    free(out);
    Expr* merged = expr_new_function(expr_copy(g->data.function.head), gargs, 1 + n);
    free(gargs);
    return merged;
}

/* ------------------------------------------------------ range handling --- */

#define SHOW_MAX_DIM 3
typedef struct { double lo[SHOW_MAX_DIM], hi[SHOW_MAX_DIM]; bool any; } Box;

static void box_add(Box* b, const double* p, size_t d) {
    for (size_t k = 0; k < d; k++) {
        if (!b->any || p[k] < b->lo[k]) b->lo[k] = p[k];
        if (!b->any || p[k] > b->hi[k]) b->hi[k] = p[k];
    }
    b->any = true;
}

/* A List of exactly d numbers, read as a point. */
static bool read_point(const Expr* e, size_t d, double* p) {
    if (!head_is(e, SYM_List) || e->data.function.arg_count != d) return false;
    for (size_t k = 0; k < d; k++)
        if (!gfx_coerce_double(e->data.function.args[k], &p[k])) return false;
    return true;
}

/* Every point inside a coordinate argument (a point, or nested Lists). */
static void coords_extent(const Expr* e, size_t d, Box* b) {
    double p[SHOW_MAX_DIM];
    if (read_point(e, d, p)) { box_add(b, p, d); return; }
    if (head_is(e, SYM_List))
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            coords_extent(e->data.function.args[i], d, b);
}

/* Extent of a primitive tree: coordinates of the geometric primitives only,
 * so a directive such as Dashing[{0.02, 0.02}] is never mistaken for a point. */
static void prims_extent(const Expr* e, size_t d, Box* b) {
    if (!e || e->type != EXPR_FUNCTION) return;
    size_t n = e->data.function.arg_count;
    if (head_is(e, SYM_List)) {
        for (size_t i = 0; i < n; i++) prims_extent(e->data.function.args[i], d, b);
        return;
    }
    if ((head_is(e, SYM_Disk) || head_is(e, SYM_Circle)) && d == 2) {
        /* Circle[] / Circle[c] / Circle[c, r]: unit circle at the origin by default. */
        double c[2] = { 0.0, 0.0 }, r = 1.0;
        if (n >= 1 && !read_point(e->data.function.args[0], 2, c)) return;
        if (n >= 2) gfx_coerce_double(e->data.function.args[1], &r);
        double lo[2] = { c[0] - r, c[1] - r }, hi[2] = { c[0] + r, c[1] + r };
        box_add(b, lo, 2); box_add(b, hi, 2);
        return;
    }
    if (n == 0) return;
    if (head_is(e, SYM_Line) || head_is(e, SYM_Point) || head_is(e, SYM_Polygon)
        || head_is(e, SYM_Arrow)) {
        coords_extent(e->data.function.args[0], d, b);
    } else if (head_is(e, SYM_Rectangle)) {
        for (size_t i = 0; i < n && i < 2; i++) coords_extent(e->data.function.args[i], d, b);
    } else if (head_is(e, SYM_Text) && n >= 2) {
        coords_extent(e->data.function.args[1], d, b);
    }
}

/* A graphic's range: its explicit numeric PlotRange where given (the 2D
 * y-only form {ymin, ymax} keeps the x extent of the primitives), else the
 * extent of its primitives. *explicit_out reports whether any axis was fixed. */
static void graphic_range(const Expr* g, size_t d, Box* b, bool* explicit_out) {
    memset(b, 0, sizeof(*b));
    *explicit_out = false;
    prims_extent(g->data.function.args[0], d, b);
    const Expr* pr = graphic_option(g, SYM_PlotRange);
    if (!head_is(pr, SYM_List)) return;
    double lo, hi;
    if (pr->data.function.arg_count == d) {
        bool all = true;
        double rl[SHOW_MAX_DIM], rh[SHOW_MAX_DIM];
        for (size_t k = 0; k < d && all; k++) {
            const Expr* ax = pr->data.function.args[k];
            all = head_is(ax, SYM_List) && ax->data.function.arg_count == 2
                && gfx_coerce_double(ax->data.function.args[0], &rl[k])
                && gfx_coerce_double(ax->data.function.args[1], &rh[k]);
        }
        if (all) {
            for (size_t k = 0; k < d; k++) { b->lo[k] = rl[k]; b->hi[k] = rh[k]; }
            b->any = true;
            *explicit_out = true;
            return;
        }
    }
    if (d == 2 && pr->data.function.arg_count == 2
        && gfx_coerce_double(pr->data.function.args[0], &lo)
        && gfx_coerce_double(pr->data.function.args[1], &hi)) {
        if (!b->any) { b->lo[0] = 0.0; b->hi[0] = 1.0; b->any = true; }
        b->lo[1] = lo; b->hi[1] = hi;
        *explicit_out = true;
    }
}

/* PlotRange -> union of the inputs' ranges, or NULL when none fixes one. */
static Expr* merged_plot_range(const PtrVec* gs, size_t d) {
    Box u;
    memset(&u, 0, sizeof(u));
    bool any_explicit = false;
    for (size_t i = 0; i < gs->n; i++) {
        Box b; bool ex;
        graphic_range(gs->v[i], d, &b, &ex);
        any_explicit |= ex;
        if (!b.any) continue;
        box_add(&u, b.lo, d);
        box_add(&u, b.hi, d);
    }
    if (!any_explicit || !u.any) return NULL;
    Expr* axes[SHOW_MAX_DIM];
    for (size_t k = 0; k < d; k++) {
        Expr* ab[2] = { expr_new_real(u.lo[k]), expr_new_real(u.hi[k]) };
        axes[k] = expr_new_function(expr_new_symbol(SYM_List), ab, 2);
    }
    Expr* val = expr_new_function(expr_new_symbol(SYM_List), axes, d);
    Expr* ra[2] = { expr_new_symbol(SYM_PlotRange), val };
    return expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
}

/* --------------------------------------------------- multi-graphic Show --- */

/* One input's primitives as a directive scope: {style..., prims}. The style
 * is the input's PlotStyle as the back ends apply it to a single-curve plot
 * (a List of styles is per curve, so its first entry). */
static Expr* scoped_prims(const Expr* g) {
    const Expr* prims = g->data.function.args[0];
    const Expr* ps = graphic_option(g, SYM_PlotStyle);
    if (head_is(ps, SYM_List))
        ps = ps->data.function.arg_count ? ps->data.function.args[0] : NULL;
    if (ps && ps->type == EXPR_SYMBOL
        && (ps->data.symbol.name == SYM_None || ps->data.symbol.name == SYM_Automatic))
        ps = NULL;

    size_t ns = !ps ? 0 : head_is(ps, SYM_List) ? ps->data.function.arg_count : 1;
    if (ns == 0 && head_is(prims, SYM_List)) return expr_copy((Expr*)prims);

    Expr** items = malloc(sizeof(Expr*) * (ns + 1));
    for (size_t k = 0; k < ns; k++)
        items[k] = expr_copy(head_is(ps, SYM_List) ? ps->data.function.args[k] : (Expr*)ps);
    items[ns] = expr_copy((Expr*)prims);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, ns + 1);
    free(items);
    return out;
}

static Expr* show_multi(const PtrVec* gs, const PtrVec* opts) {
    const Expr* g1 = gs->v[0];
    const char* head = g1->data.function.head->data.symbol.name;
    size_t d = head == SYM_Graphics3D ? 3 : 2;

    /* Overlaid primitives, one scope per input. */
    Expr** scopes = malloc(sizeof(Expr*) * gs->n);
    /* (The 3D renderer takes a surface's style from the PlotStyle option
     * rather than from directives, so 3D inputs are overlaid as they are
     * and keep g1's PlotStyle.) */
    for (size_t i = 0; i < gs->n; i++)
        scopes[i] = d == 3 ? expr_copy(gs->v[i]->data.function.args[0])
                           : scoped_prims(gs->v[i]);
    Expr* prims = expr_new_function(expr_new_symbol(SYM_List), scopes, gs->n);
    free(scopes);

    /* Concatenated legend entries. */
    size_t nleg = 0;
    for (size_t i = 0; i < gs->n; i++)
        for (size_t j = 1; j < gs->v[i]->data.function.arg_count; j++)
            if (head_is(gs->v[i]->data.function.args[j], SYM_PlotLegendData))
                nleg += gs->v[i]->data.function.args[j]->data.function.arg_count;

    size_t cap = opts->n + g1->data.function.arg_count + 3;
    Expr** out = malloc(sizeof(Expr*) * cap);
    size_t n = 0;

    /* Show's own options win. */
    for (size_t i = 0; i < opts->n; i++) out[n++] = expr_copy((Expr*)opts->v[i]);

    /* Then g1's, minus what Show overrides, what was baked into the scopes
     * (PlotStyle), what is merged (PlotRange, legends) and what no longer
     * describes the combination ($PlotResample). */
    for (size_t i = 1; i < g1->data.function.arg_count; i++) {
        const Expr* a = g1->data.function.args[i];
        if (show_is_rule(a)) {
            const Expr* lhs = a->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL
                && ((lhs->data.symbol.name == SYM_PlotStyle && d == 2)
                    || lhs->data.symbol.name == SYM_PlotRange
                    || opts_have(opts, lhs->data.symbol.name)))
                continue;
        } else if (head_is(a, SYM_PlotResample) || head_is(a, SYM_PlotLegendData)) {
            continue;
        }
        out[n++] = expr_copy((Expr*)a);
    }

    if (!opts_have(opts, SYM_PlotRange)) {
        Expr* pr = merged_plot_range(gs, d);
        if (pr) out[n++] = pr;
    }

    if (nleg > 0) {
        Expr** ents = malloc(sizeof(Expr*) * nleg);
        size_t k = 0;
        for (size_t i = 0; i < gs->n; i++)
            for (size_t j = 1; j < gs->v[i]->data.function.arg_count; j++) {
                const Expr* ld = gs->v[i]->data.function.args[j];
                if (!head_is(ld, SYM_PlotLegendData)) continue;
                for (size_t m = 0; m < ld->data.function.arg_count; m++)
                    ents[k++] = expr_copy(ld->data.function.args[m]);
            }
        out[n++] = expr_new_function(expr_new_symbol(SYM_PlotLegendData), ents, nleg);
        free(ents);
    }

    Expr** gargs = malloc(sizeof(Expr*) * (1 + n));
    gargs[0] = prims;
    for (size_t i = 0; i < n; i++) gargs[1 + i] = out[i];
    free(out);
    Expr* combined = expr_new_function(expr_new_symbol(head), gargs, 1 + n);
    free(gargs);
    return combined;
}

Expr* builtin_show(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    /* Graphics (possibly in Lists) first, then trailing options. */
    PtrVec gs = { NULL, 0, 0 }, opts = { NULL, 0, 0 };
    bool ok = true;
    for (size_t i = 0; i < argc && ok; i++) {
        const Expr* a = res->data.function.args[i];
        if (show_is_rule(a)) {
            if (gs.n == 0) ok = false;              /* Show[opts] alone */
            else pv_push(&opts, a);
        } else if (opts.n > 0) {
            ok = false;                              /* graphics after options */
        } else {
            ok = collect_graphics(a, &gs);
        }
    }
    if (ok && gs.n == 0) ok = false;
    /* All inputs must share one container head: 2D and 3D do not mix. */
    for (size_t i = 1; ok && i < gs.n; i++)
        if (gs.v[i]->data.function.head->data.symbol.name
            != gs.v[0]->data.function.head->data.symbol.name) ok = false;

    /* Rendering is owned by the front end (the REPL renders any top-level
     * Graphics[...] result); Show[] just builds the object and returns it.
     * See the auto-display block in repl.c::process_input. */
    Expr* out = NULL;
    if (ok) out = gs.n == 1 ? show_single(gs.v[0], &opts) : show_multi(&gs, &opts);
    free(gs.v);
    free(opts.v);
    return out;
}

#ifndef USE_GRAPHICS
#include <stdio.h>
void graphics_show(const Expr* graphics_expr) {
    (void)graphics_expr;
    printf("Graphics: not rendered -- graphics support not compiled in "
           "(install raylib and rebuild).\n");
}
void graphics3d_show(const Expr* graphics3d_expr) {
    (void)graphics3d_expr;
    printf("Graphics3D: not rendered -- graphics support not compiled in "
           "(install raylib and rebuild).\n");
}
#endif
