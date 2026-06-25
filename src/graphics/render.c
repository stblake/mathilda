/* render.c — Raylib backend: graphics_show() (see show.h).
 *
 * Coordinate convention used throughout this file: math-space is y-up;
 * Raylib draw-space is y-down. Every world point's y is negated exactly
 * once, right where it's converted to a Vector2 for a draw call -- there
 * are no other sign flips anywhere else in this file. */

#include "show.h"
#include "plot.h"
#include "render.h"
#include "sampling.h"
#include "hershey_font.h"
#include "primitives.h"
#include "sym_names.h"
#include "print.h"
#include "eval.h"
#include <raylib.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define HERSHEY_CAP_HEIGHT 7.0

/* ---------------- Expr coercion ---------------- */

static bool expr_to_d(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational) {
        double n, d;
        if (expr_to_d(e->data.function.args[0], &n) && expr_to_d(e->data.function.args[1], &d) && d != 0) {
            *out = n / d;
            return true;
        }
    }
    /* A symbolic-but-numeric coordinate (Pi/2, Sqrt[2], E, GoldenRatio, or
     * arithmetic of these) isn't a plain literal -- plain evaluate() alone
     * never reduces Pi to a number, only N[] does (e.g. Point[{Pi/2, 1}]
     * would otherwise silently fail to draw at all). Numericize it the same
     * way plot.c's numericize_bound resolves Plot's iterator bounds. Hand-
     * authored Graphics[]/Epilog/Prolog content is typically small, so the
     * extra evaluate() call here (vs. the sampler's already-literal curve
     * points, which never reach this fallback) is not a hot-path concern. */
    if (e->type == EXPR_SYMBOL || e->type == EXPR_FUNCTION) {
        Expr* n_arg[1] = { expr_copy((Expr*)e) };
        Expr* n_call = expr_new_function(expr_new_symbol("N"), n_arg, 1);
        Expr* result = evaluate(n_call);
        expr_free(n_call);
        bool ok = false;
        if (result->type == EXPR_REAL) { *out = result->data.real; ok = true; }
        else if (result->type == EXPR_INTEGER) { *out = (double)result->data.integer; ok = true; }
#ifdef USE_MPFR
        else if (result->type == EXPR_MPFR) { *out = mpfr_get_d(result->data.mpfr, MPFR_RNDN); ok = true; }
#endif
        expr_free(result);
        if (ok && isfinite(*out)) return true;
    }
    return false;
}

static bool expr_point(const Expr* e, double* x, double* y) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    return expr_to_d(e->data.function.args[0], x) && expr_to_d(e->data.function.args[1], y);
}

static bool looks_like_point(const Expr* e) {
    double x, y;
    return expr_point(e, &x, &y);
}

/* Resolve a Circle[...]/Disk[...] node's centre and radius, applying
 * Mathematica's defaults: Circle[] is the unit circle at the origin,
 * Circle[{x,y}] takes radius 1, Circle[{x,y}, r] is fully explicit.
 * Returns false only when a supplied argument is present but unreadable. */
static bool circle_params(const Expr* e, double* cx, double* cy, double* r) {
    size_t n = e->data.function.arg_count;
    *cx = 0; *cy = 0; *r = 1;
    if (n >= 1 && !expr_point(e->data.function.args[0], cx, cy)) return false;
    if (n >= 2 && !expr_to_d(e->data.function.args[1], r)) return false;
    return true;
}

static bool expr_is_sym(const Expr* e, const char* sym) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == sym;
}

/* ---------------- Color helpers ---------------- */

static RGBA8 rgba_from_rgbcolor(const Expr* e) {
    double r = 0, g = 0, b = 0, a = 1;
    size_t n = e->data.function.arg_count;
    if (n >= 3) {
        expr_to_d(e->data.function.args[0], &r);
        expr_to_d(e->data.function.args[1], &g);
        expr_to_d(e->data.function.args[2], &b);
    }
    if (n >= 4) expr_to_d(e->data.function.args[3], &a);
    RGBA8 c = { (unsigned char)(r * 255), (unsigned char)(g * 255), (unsigned char)(b * 255), (unsigned char)(a * 255) };
    return c;
}

static RGBA8 rgba_from_graylevel(const Expr* e) {
    double g = 0, a = 1;
    size_t n = e->data.function.arg_count;
    if (n >= 1) expr_to_d(e->data.function.args[0], &g);
    if (n >= 2) expr_to_d(e->data.function.args[1], &a);
    RGBA8 c = { (unsigned char)(g * 255), (unsigned char)(g * 255), (unsigned char)(g * 255), (unsigned char)(a * 255) };
    return c;
}

/* Hue[h] (s=v=1), Hue[h,s,v], Hue[h,s,v,a] -- standard HSB/HSV to RGB. */
static RGBA8 rgba_from_hue(const Expr* e) {
    double h = 0, s = 1, v = 1, a = 1;
    size_t n = e->data.function.arg_count;
    if (n >= 1) expr_to_d(e->data.function.args[0], &h);
    if (n >= 3) { expr_to_d(e->data.function.args[1], &s); expr_to_d(e->data.function.args[2], &v); }
    if (n >= 4) expr_to_d(e->data.function.args[3], &a);
    h = h - floor(h); /* wrap to [0,1), matching Mathematica's Hue */
    if (s < 0) s = 0; if (s > 1) s = 1;
    if (v < 0) v = 0; if (v > 1) v = 1;
    double r, g, b;
    if (s <= 0.0) {
        r = g = b = v;
    } else {
        double hh = h * 6.0;
        int i = (int)floor(hh);
        double f = hh - i;
        double p = v * (1.0 - s);
        double q = v * (1.0 - f * s);
        double t = v * (1.0 - (1.0 - f) * s);
        switch (((i % 6) + 6) % 6) {
            case 0: r = v; g = t; b = p; break;
            case 1: r = q; g = v; b = p; break;
            case 2: r = p; g = v; b = t; break;
            case 3: r = p; g = q; b = v; break;
            case 4: r = t; g = p; b = v; break;
            default: r = v; g = p; b = q; break;
        }
    }
    RGBA8 c = { (unsigned char)(r * 255), (unsigned char)(g * 255), (unsigned char)(b * 255), (unsigned char)(a * 255) };
    return c;
}

static Color to_raylib(RGBA8 c) { Color out = { c.r, c.g, c.b, c.a }; return out; }

/* Resolves any recognized color-literal head (RGBColor, GrayLevel, Hue) to
 * an RGBA8, leaving *out untouched (returning false) for anything else --
 * the single place every color-bearing option/directive in this file goes
 * through, so adding a future color form means touching only this. */
static bool resolve_color(const Expr* e, RGBA8* out) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h == SYM_RGBColor)  { *out = rgba_from_rgbcolor(e); return true; }
    if (h == SYM_GrayLevel) { *out = rgba_from_graylevel(e); return true; }
    if (h == SYM_Hue)       { *out = rgba_from_hue(e); return true; }
    return false;
}

/* ---------------- Option parsing ---------------- */

/* Frame edge indices, used for the per-edge frame_edge[]/frame_ticks[]
 * arrays. The order matches Mathematica's Frame -> {{left, right},
 * {bottom, top}} once flattened. */
enum { FR_LEFT = 0, FR_RIGHT = 1, FR_BOTTOM = 2, FR_TOP = 3 };

typedef struct {
    bool axes;
    bool frame;             /* any frame edge drawn at all */
    bool frame_edge[4];     /* per edge: draw this frame line (L,R,B,T) */
    bool frame_ticks[4];    /* per edge: draw ticks/labels (gated by frame_edge) */
    RGBA8 frame_color;      /* FrameStyle colour for the box, ticks and labels */
    bool x_auto;        /* x extent auto-computed from the data */
    bool y_auto;        /* y extent auto-computed from the data */
    bool clip_outliers; /* on the auto y-axis, clamp sparse asymptote spikes */
    PlotRange2D range;  /* explicit bounds; only the non-auto axes are read */
    double aspect_ratio; /* height/width; <= 0 means Automatic (true geometry) */
    bool aspect_full;    /* AspectRatio -> Full: stretch to fill the window */
    bool height_pinned;  /* ImageSize -> {w,h}: height is fixed, AspectRatio
                          * shapes the data within rather than the window */
    RGBA8 style_color;
    RGBA8 background;
    long width, height;
    const Expr* axes_label; /* borrowed */
    const Expr* plot_label; /* borrowed */

    bool axes_origin_set;     /* AxesOrigin given explicitly */
    double axes_origin_x, axes_origin_y;
    RGBA8 axes_color;         /* AxesStyle: axis lines + tick marks */
    RGBA8 ticks_color;        /* TicksStyle: tick label text */
    const Expr* frame_label;  /* borrowed; FrameLabel -> {xlabel, ylabel} */
    bool rotate_label;        /* RotateLabel: FrameLabel y-label orientation */
    double pad_x_frac, pad_y_frac; /* PlotRangePadding (auto-fit only) */
    bool grid_x_on, grid_y_on;     /* GridLines: draw vertical/horizontal lines at all */
    const Expr* grid_x_list;  /* borrowed explicit x positions; NULL = Automatic majors */
    const Expr* grid_y_list;  /* borrowed explicit y positions; NULL = Automatic majors */
    RGBA8 grid_color;         /* GridLinesStyle */
    const Expr* prolog;       /* borrowed; drawn first, in data space */
    const Expr* epilog;       /* borrowed; drawn last, in data space */
} GfxOptions;

/* True when `e` is the symbol True or All (the "on" forms for Frame and
 * FrameTicks edge settings); False/None/anything else reads as "off". */
static bool frame_edge_on(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && (e->data.symbol == SYM_True || e->data.symbol == SYM_All);
}

/* True for the "show ticks" forms of a FrameTicks edge setting. Automatic and
 * True request ticks; None and False suppress them. */
static bool frame_ticks_on(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && (e->data.symbol == SYM_Automatic || e->data.symbol == SYM_True);
}

/* Read a nested {{left, right}, {bottom, top}} option value into per-edge
 * results via `pred`. Returns false if the value isn't that 2x2 List shape,
 * leaving `out` untouched so the caller can keep its default. */
static bool parse_edge_pairs(const Expr* rhs, bool (*pred)(const Expr*), bool out[4]) {
    if (!rhs || rhs->type != EXPR_FUNCTION || rhs->data.function.head->type != EXPR_SYMBOL
        || rhs->data.function.head->data.symbol != SYM_List || rhs->data.function.arg_count != 2)
        return false;
    const Expr* lr = rhs->data.function.args[0];
    const Expr* bt = rhs->data.function.args[1];
    if (lr->type != EXPR_FUNCTION || lr->data.function.arg_count != 2
        || bt->type != EXPR_FUNCTION || bt->data.function.arg_count != 2)
        return false;
    out[FR_LEFT]   = pred(lr->data.function.args[0]);
    out[FR_RIGHT]  = pred(lr->data.function.args[1]);
    out[FR_BOTTOM] = pred(bt->data.function.args[0]);
    out[FR_TOP]    = pred(bt->data.function.args[1]);
    return true;
}

static void gfx_options_parse(const Expr* graphics, GfxOptions* o) {
    o->axes = false;
    o->frame = false;
    for (int e = 0; e < 4; e++) { o->frame_edge[e] = false; o->frame_ticks[e] = true; }
    o->frame_color = (RGBA8){ 90, 90, 90, 255 }; /* matches DARKGRAY-ish axis tone */
    o->x_auto = true;
    o->y_auto = true;
    o->clip_outliers = true;
    o->aspect_ratio = -1.0;
    o->aspect_full = false;
    o->height_pinned = false;
    o->style_color = (RGBA8){ 30, 80, 180, 255 };
    o->background = (RGBA8){ 255, 255, 255, 255 };
    o->width = 800;
    o->height = 600;
    o->axes_label = NULL;
    o->plot_label = NULL;

    o->axes_origin_set = false;
    o->axes_origin_x = 0.0; o->axes_origin_y = 0.0;
    o->axes_color = (RGBA8){ 90, 90, 90, 255 };
    o->ticks_color = (RGBA8){ 90, 90, 90, 255 };
    o->frame_label = NULL;
    o->rotate_label = true;
    o->pad_x_frac = 0.08; o->pad_y_frac = 0.08;
    o->grid_x_on = false; o->grid_y_on = false;
    o->grid_x_list = NULL; o->grid_y_list = NULL;
    o->grid_color = (RGBA8){ 210, 210, 210, 255 };
    o->prolog = NULL;
    o->epilog = NULL;

    size_t argc = graphics->data.function.arg_count;

    /* LabelStyle -> color seeds the axis/ticks/frame text-and-line defaults
     * before the main pass below runs, so any of AxesStyle/TicksStyle/
     * FrameStyle the caller also gives still overrides it, regardless of
     * each option's position in the argument list. */
    for (size_t i = 1; i < argc; i++) {
        const Expr* opt = graphics->data.function.args[i];
        if (!opt || opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) continue;
        const Expr* h0 = opt->data.function.head;
        if (!h0 || h0->type != EXPR_SYMBOL || (h0->data.symbol != SYM_Rule && h0->data.symbol != SYM_RuleDelayed)) continue;
        const Expr* lhs0 = opt->data.function.args[0];
        if (lhs0->type == EXPR_SYMBOL && lhs0->data.symbol == SYM_LabelStyle) {
            RGBA8 c;
            if (resolve_color(opt->data.function.args[1], &c)) {
                o->axes_color = c; o->ticks_color = c; o->frame_color = c;
            }
        }
    }

    for (size_t i = 1; i < argc; i++) {
        const Expr* opt = graphics->data.function.args[i];
        if (!opt || opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) continue;
        const Expr* h = opt->data.function.head;
        if (!h || h->type != EXPR_SYMBOL) continue;
        if (h->data.symbol != SYM_Rule && h->data.symbol != SYM_RuleDelayed) continue;
        const Expr* lhs = opt->data.function.args[0];
        const Expr* rhs = opt->data.function.args[1];
        if (lhs->type != EXPR_SYMBOL) continue;
        const char* name = lhs->data.symbol;

        if (name == SYM_Axes) {
            o->axes = expr_is_sym(rhs, SYM_True);
        } else if (name == SYM_Frame) {
            /* Frame -> True       : box on all four edges.
             * Frame -> False/None : no frame.
             * Frame -> {{l,r},{b,t}} : per-edge, each True/False. */
            if (frame_edge_on(rhs)) {
                for (int e = 0; e < 4; e++) o->frame_edge[e] = true;
            } else if (rhs->type == EXPR_SYMBOL
                       && (rhs->data.symbol == SYM_False || rhs->data.symbol == SYM_None)) {
                for (int e = 0; e < 4; e++) o->frame_edge[e] = false;
            } else {
                parse_edge_pairs(rhs, frame_edge_on, o->frame_edge);
            }
            o->frame = o->frame_edge[0] || o->frame_edge[1] || o->frame_edge[2] || o->frame_edge[3];
        } else if (name == SYM_FrameTicks) {
            /* FrameTicks -> Automatic/True : ticks on every drawn edge (default).
             * FrameTicks -> None/False     : frame box but no ticks or labels.
             * FrameTicks -> {{l,r},{b,t}}  : per-edge, each Automatic/None. */
            if (rhs->type == EXPR_SYMBOL
                && (rhs->data.symbol == SYM_None || rhs->data.symbol == SYM_False)) {
                for (int e = 0; e < 4; e++) o->frame_ticks[e] = false;
            } else if (frame_ticks_on(rhs)) {
                for (int e = 0; e < 4; e++) o->frame_ticks[e] = true;
            } else {
                parse_edge_pairs(rhs, frame_ticks_on, o->frame_ticks);
            }
        } else if (name == SYM_FrameStyle) {
            resolve_color(rhs, &o->frame_color);
        } else if (name == SYM_AspectRatio) {
            /* Automatic (default): keep the <=0 sentinel -> true geometry.
             * Full: fill the window; resolved to height/width after the loop
             * since ImageSize may follow AspectRatio in the option list.
             * Numeric a>0: explicit height-to-width ratio. */
            if (expr_is_sym(rhs, SYM_Automatic)) {
                o->aspect_ratio = -1.0; o->aspect_full = false;
            } else if (expr_is_sym(rhs, SYM_Full)) {
                o->aspect_full = true;
            } else {
                double v;
                if (expr_to_d(rhs, &v) && v > 0) { o->aspect_ratio = v; o->aspect_full = false; }
            }
        } else if (name == SYM_PlotStyle) {
            resolve_color(rhs, &o->style_color);
        } else if (name == SYM_Background) {
            resolve_color(rhs, &o->background);
        } else if (name == SYM_ImageSize) {
            if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL
                && rhs->data.function.head->data.symbol == SYM_List && rhs->data.function.arg_count == 2) {
                double w, hh;
                if (expr_to_d(rhs->data.function.args[0], &w)) o->width = (long)w;
                if (expr_to_d(rhs->data.function.args[1], &hh)) { o->height = (long)hh; o->height_pinned = true; }
            } else {
                double s;
                if (expr_to_d(rhs, &s) && s > 0) { o->width = (long)s; o->height = (long)(s * 0.75); }
            }
        } else if (name == SYM_PlotRange) {
            /* Symbol forms first, mirroring Wolfram:
             *   Automatic -- data-driven y with spike clipping (the default)
             *   All       -- show every sampled point, no spike clipping */
            if (rhs->type == EXPR_SYMBOL
                && (rhs->data.symbol == SYM_Automatic || rhs->data.symbol == SYM_All)) {
                o->clip_outliers = (rhs->data.symbol != SYM_All);
            }
            /* Two numeric forms:
             *   {{xmin, xmax}, {ymin, ymax}}  -- fix both axes
             *   {ymin, ymax}                  -- fix y only, x stays automatic */
            else if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL
                && rhs->data.function.head->data.symbol == SYM_List
                && rhs->data.function.arg_count == 2) {
                const Expr* a0 = rhs->data.function.args[0];
                const Expr* a1 = rhs->data.function.args[1];
                double xmin, xmax, ymin, ymax;
                if (a0->type == EXPR_FUNCTION && a0->data.function.arg_count == 2
                    && a1->type == EXPR_FUNCTION && a1->data.function.arg_count == 2
                    && expr_to_d(a0->data.function.args[0], &xmin) && expr_to_d(a0->data.function.args[1], &xmax)
                    && expr_to_d(a1->data.function.args[0], &ymin) && expr_to_d(a1->data.function.args[1], &ymax)) {
                    o->x_auto = false; o->y_auto = false;
                    o->range.xmin = xmin; o->range.xmax = xmax;
                    o->range.ymin = ymin; o->range.ymax = ymax;
                } else if (expr_to_d(a0, &ymin) && expr_to_d(a1, &ymax)) {
                    o->y_auto = false;
                    o->range.ymin = ymin; o->range.ymax = ymax;
                }
            }
        } else if (name == SYM_AxesLabel) {
            o->axes_label = rhs;
        } else if (name == SYM_PlotLabel) {
            o->plot_label = rhs;
        } else if (name == SYM_AxesOrigin) {
            double ox, oy;
            if (expr_point(rhs, &ox, &oy)) {
                o->axes_origin_set = true;
                o->axes_origin_x = ox; o->axes_origin_y = oy;
            }
        } else if (name == SYM_AxesStyle) {
            resolve_color(rhs, &o->axes_color);
        } else if (name == SYM_TicksStyle) {
            resolve_color(rhs, &o->ticks_color);
        } else if (name == SYM_FrameLabel) {
            o->frame_label = rhs;
        } else if (name == SYM_RotateLabel) {
            if (expr_is_sym(rhs, SYM_False)) o->rotate_label = false;
            else if (expr_is_sym(rhs, SYM_True)) o->rotate_label = true;
        } else if (name == SYM_PlotRangePadding) {
            if (expr_is_sym(rhs, SYM_None)) {
                o->pad_x_frac = 0.0; o->pad_y_frac = 0.0;
            } else if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL
                       && rhs->data.function.head->data.symbol == SYM_List && rhs->data.function.arg_count == 2) {
                double px, py;
                if (expr_to_d(rhs->data.function.args[0], &px)) o->pad_x_frac = px;
                if (expr_to_d(rhs->data.function.args[1], &py)) o->pad_y_frac = py;
            } else {
                double p;
                if (expr_to_d(rhs, &p)) { o->pad_x_frac = p; o->pad_y_frac = p; }
            }
        } else if (name == SYM_GridLines) {
            /* None (default) : no grid.
             * Automatic       : both axes, at the same major ticks as Axes/Frame.
             * {xspec, yspec}  : independently None/Automatic/an explicit List[...]
             *                   of positions per axis. */
            if (expr_is_sym(rhs, SYM_None)) {
                o->grid_x_on = false; o->grid_y_on = false;
            } else if (expr_is_sym(rhs, SYM_Automatic)) {
                o->grid_x_on = true; o->grid_y_on = true;
                o->grid_x_list = NULL; o->grid_y_list = NULL;
            } else if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL
                       && rhs->data.function.head->data.symbol == SYM_List && rhs->data.function.arg_count == 2) {
                const Expr* xs = rhs->data.function.args[0];
                const Expr* ys = rhs->data.function.args[1];
                if (expr_is_sym(xs, SYM_None)) o->grid_x_on = false;
                else if (expr_is_sym(xs, SYM_Automatic)) { o->grid_x_on = true; o->grid_x_list = NULL; }
                else if (xs->type == EXPR_FUNCTION && xs->data.function.head->type == EXPR_SYMBOL
                         && xs->data.function.head->data.symbol == SYM_List) { o->grid_x_on = true; o->grid_x_list = xs; }
                if (expr_is_sym(ys, SYM_None)) o->grid_y_on = false;
                else if (expr_is_sym(ys, SYM_Automatic)) { o->grid_y_on = true; o->grid_y_list = NULL; }
                else if (ys->type == EXPR_FUNCTION && ys->data.function.head->type == EXPR_SYMBOL
                         && ys->data.function.head->data.symbol == SYM_List) { o->grid_y_on = true; o->grid_y_list = ys; }
            }
        } else if (name == SYM_GridLinesStyle) {
            resolve_color(rhs, &o->grid_color);
        } else if (name == SYM_Prolog) {
            o->prolog = rhs;
        } else if (name == SYM_Epilog) {
            o->epilog = rhs;
        }
    }

    if (o->width < 100) o->width = 100;
    if (o->height < 100) o->height = 100;

    /* AspectRatio -> Full: a height-to-width ratio equal to the window's own
     * ratio makes the data fill it edge to edge (fit_by_width == fit_by_height
     * in the scaling below). Resolved here so it tracks the final ImageSize. */
    if (o->aspect_full) o->aspect_ratio = (double)o->height / (double)o->width;
}

long gfx_window_height(long width, long height, double aspect_ratio,
                       bool aspect_full, bool height_pinned,
                       double data_w, double data_h) {
    if (height_pinned || aspect_full) return height; /* keep the ImageSize box */
    double a = aspect_ratio > 0 ? aspect_ratio : (data_h / data_w);
    if (!isfinite(a) || a <= 0) return height;
    long h = (long)(width * a + 0.5);
    if (h < 100)  h = 100;
    if (h > 2000) h = 2000;
    return h;
}

/* ---------------- Bounding box ---------------- */

static void update_bbox(PlotRange2D* bb, double x, double y) {
    if (x < bb->xmin) bb->xmin = x;
    if (x > bb->xmax) bb->xmax = x;
    if (y < bb->ymin) bb->ymin = y;
    if (y > bb->ymax) bb->ymax = y;
}

static void compute_bbox(const Expr* node, PlotRange2D* bb) {
    if (!node || node->type != EXPR_FUNCTION) return;
    const Expr* h = node->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return;
    const char* name = h->data.symbol;
    size_t n = node->data.function.arg_count;

    if (name == SYM_List) {
        double x, y;
        if (n == 2 && expr_to_d(node->data.function.args[0], &x) && expr_to_d(node->data.function.args[1], &y)) {
            update_bbox(bb, x, y);
            return;
        }
        for (size_t i = 0; i < n; i++) compute_bbox(node->data.function.args[i], bb);
    } else if ((name == SYM_Point || name == SYM_Line || name == SYM_Polygon) && n >= 1) {
        compute_bbox(node->data.function.args[0], bb);
    } else if (name == SYM_Rectangle && n >= 2) {
        compute_bbox(node->data.function.args[0], bb);
        compute_bbox(node->data.function.args[1], bb);
    } else if (name == SYM_Circle || name == SYM_Disk) {
        double cx, cy, r;
        if (circle_params(node, &cx, &cy, &r)) {
            update_bbox(bb, cx - r, cy - r);
            update_bbox(bb, cx + r, cy + r);
        }
    } else if (name == SYM_Text && n >= 2) {
        compute_bbox(node->data.function.args[1], bb);
    }
}

/* Append a y-value to a growable buffer, doubling capacity as needed. */
static void ybuf_push(double y, double** buf, size_t* len, size_t* cap) {
    if (*len == *cap) {
        *cap = *cap ? *cap * 2 : 64;
        *buf = realloc(*buf, sizeof(double) * *cap);
    }
    (*buf)[(*len)++] = y;
}

static int cmp_double_asc(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Append the absolute slope |dy/dx| of each segment of one polyline. */
static void polyline_slopes(const Expr* pts, double** buf, size_t* len, size_t* cap) {
    size_t n = pts->data.function.arg_count;
    double px = 0, py = 0; bool have = false;
    for (size_t i = 0; i < n; i++) {
        double x, y;
        if (!expr_point(pts->data.function.args[i], &x, &y)) { have = false; continue; }
        if (have) { double dx = x - px; if (dx != 0.0) ybuf_push(fabs((y - py) / dx), buf, len, cap); }
        px = x; py = y; have = true;
    }
}

/* Gather every Line segment's |slope| (recursing through List wrappers and
 * both Line point-list shapes), mirroring draw_polyline's structure walk. */
static void collect_slopes(const Expr* node, double** buf, size_t* len, size_t* cap) {
    if (!node || node->type != EXPR_FUNCTION) return;
    const Expr* h = node->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return;
    const char* name = h->data.symbol;
    size_t n = node->data.function.arg_count;

    if (name == SYM_Line && n >= 1 && node->data.function.args[0]->type == EXPR_FUNCTION) {
        const Expr* arg = node->data.function.args[0];
        size_t m = arg->data.function.arg_count;
        if (m > 0 && looks_like_point(arg->data.function.args[0])) {
            polyline_slopes(arg, buf, len, cap);
        } else {
            for (size_t i = 0; i < m; i++)
                if (arg->data.function.args[i]->type == EXPR_FUNCTION)
                    polyline_slopes(arg->data.function.args[i], buf, len, cap);
        }
    } else if (name == SYM_List) {
        for (size_t i = 0; i < n; i++) collect_slopes(node->data.function.args[i], buf, len, cap);
    }
}

/* True when a curve in `prims` contains a near-vertical runaway -- the
 * signature of an asymptote. tan(x) never lands exactly on a pole, so it has
 * no break to key off; instead its climb produces a segment whose |slope| is
 * astronomically larger than the curve's typical slope. A max-to-median slope
 * ratio above ~1000 cleanly separates poles (Tan ~1e17, 1/x ~1e4, Gamma ~1e5)
 * from merely steep smooth curves (x^5 ~200, x^3 ~16, Exp ~4), which the
 * adaptive sampler keeps slope-bounded by refining. Only such curves get the
 * spike-clipping band; everything smooth keeps its full, legitimate extent. */
static bool prims_have_runaway(const Expr* prims) {
    double* sl = NULL; size_t n = 0, cap = 0;
    collect_slopes(prims, &sl, &n, &cap);
    bool runaway = false;
    if (n >= 4) {
        qsort(sl, n, sizeof(double), cmp_double_asc);
        double med = sl[n / 2], mx = sl[n - 1];
        if (med > 0.0 && mx > 1000.0 * med) runaway = true;
    }
    free(sl);
    return runaway;
}

/* Collect every primitive y-value into `buf` (same node set compute_bbox
 * visits), so plot_robust_yrange can pick a spike-resistant vertical band.
 * Mirrors compute_bbox's traversal exactly to stay consistent with it. */
static void gather_ys(const Expr* node, double** buf, size_t* len, size_t* cap) {
    if (!node || node->type != EXPR_FUNCTION) return;
    const Expr* h = node->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return;
    const char* name = h->data.symbol;
    size_t n = node->data.function.arg_count;

    if (name == SYM_List) {
        double x, y;
        if (n == 2 && expr_to_d(node->data.function.args[0], &x) && expr_to_d(node->data.function.args[1], &y)) {
            ybuf_push(y, buf, len, cap);
            return;
        }
        for (size_t i = 0; i < n; i++) gather_ys(node->data.function.args[i], buf, len, cap);
    } else if ((name == SYM_Point || name == SYM_Line || name == SYM_Polygon) && n >= 1) {
        gather_ys(node->data.function.args[0], buf, len, cap);
    } else if (name == SYM_Rectangle && n >= 2) {
        gather_ys(node->data.function.args[0], buf, len, cap);
        gather_ys(node->data.function.args[1], buf, len, cap);
    } else if (name == SYM_Circle || name == SYM_Disk) {
        double cx, cy, r;
        if (circle_params(node, &cx, &cy, &r)) { ybuf_push(cy - r, buf, len, cap); ybuf_push(cy + r, buf, len, cap); }
    } else if (name == SYM_Text && n >= 2) {
        gather_ys(node->data.function.args[1], buf, len, cap);
    }
}

/* ---------------- Drawing ---------------- */

typedef struct {
    Color color;
    float thickness;  /* world units; <= 0 means hairline */
    float point_size; /* world units (radius) */
    float text_scale; /* world units per Hershey grid unit, for Text[] */
    float yscale;     /* data-y -> render-y factor (non-uniform PlotRange/AspectRatio) */
} DrawState;

static void draw_polyline(const Expr* pts_list, const DrawState* state) {
    size_t n = pts_list->data.function.arg_count;
    Vector2 prev = { 0, 0 };
    bool have_prev = false;
    for (size_t i = 0; i < n; i++) {
        double x, y;
        if (!expr_point(pts_list->data.function.args[i], &x, &y)) { have_prev = false; continue; }
        Vector2 cur = { (float)x, (float)(-y * state->yscale) };
        if (have_prev) {
            if (state->thickness > 0.0001f) DrawLineEx(prev, cur, state->thickness, state->color);
            else DrawLineV(prev, cur, state->color);
        }
        prev = cur;
        have_prev = true;
    }
}

static void draw_primitive(const Expr* node, DrawState* state) {
    if (!node || node->type != EXPR_FUNCTION) return;
    const Expr* h = node->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return;
    const char* name = h->data.symbol;
    size_t n = node->data.function.arg_count;

    if (name == SYM_List) {
        for (size_t i = 0; i < n; i++) draw_primitive(node->data.function.args[i], state);
        return;
    }
    {
        RGBA8 c;
        if (resolve_color(node, &c)) { state->color = to_raylib(c); return; }
    }
    if (name == SYM_Opacity) {
        double a;
        if (n >= 1 && expr_to_d(node->data.function.args[0], &a)) state->color.a = (unsigned char)(a * 255);
        return;
    }
    if (name == SYM_Thickness) {
        double t;
        if (n >= 1 && expr_to_d(node->data.function.args[0], &t)) state->thickness = (float)t;
        return;
    }
    if (name == SYM_PointSize) {
        double s;
        if (n >= 1 && expr_to_d(node->data.function.args[0], &s)) state->point_size = (float)s;
        return;
    }
    if (name == SYM_Point && n >= 1) {
        const Expr* arg = node->data.function.args[0];
        double x, y;
        if (expr_point(arg, &x, &y)) {
            DrawCircleV((Vector2){ (float)x, (float)(-y * state->yscale) }, state->point_size, state->color);
        } else if (arg->type == EXPR_FUNCTION) {
            for (size_t i = 0; i < arg->data.function.arg_count; i++) {
                double px, py;
                if (expr_point(arg->data.function.args[i], &px, &py))
                    DrawCircleV((Vector2){ (float)px, (float)(-py * state->yscale) }, state->point_size, state->color);
            }
        }
        return;
    }
    if (name == SYM_Line && n >= 1 && node->data.function.args[0]->type == EXPR_FUNCTION) {
        const Expr* arg = node->data.function.args[0];
        size_t m = arg->data.function.arg_count;
        if (m > 0 && looks_like_point(arg->data.function.args[0])) {
            draw_polyline(arg, state);
        } else {
            for (size_t i = 0; i < m; i++) {
                const Expr* sub = arg->data.function.args[i];
                if (sub->type == EXPR_FUNCTION) draw_polyline(sub, state);
            }
        }
        return;
    }
    if (name == SYM_Rectangle && n >= 2) {
        double x1, y1, x2, y2;
        if (expr_point(node->data.function.args[0], &x1, &y1) && expr_point(node->data.function.args[1], &x2, &y2)) {
            double xmin = x1 < x2 ? x1 : x2, xmax = x1 < x2 ? x2 : x1;
            double ymin = y1 < y2 ? y1 : y2, ymax = y1 < y2 ? y2 : y1;
            Rectangle r = { (float)xmin, (float)(-ymax * state->yscale),
                            (float)(xmax - xmin), (float)((ymax - ymin) * state->yscale) };
            DrawRectangleRec(r, state->color);
        }
        return;
    }
    if (name == SYM_Circle || name == SYM_Disk) {
        double cx, cy, r;
        if (circle_params(node, &cx, &cy, &r)) {
            Vector2 center = { (float)cx, (float)(-cy * state->yscale) };
            if (name == SYM_Disk) DrawCircleV(center, (float)r, state->color);
            else DrawCircleLinesV(center, (float)r, state->color);
        }
        return;
    }
    if (name == SYM_Polygon && n >= 1 && node->data.function.args[0]->type == EXPR_FUNCTION) {
        const Expr* arg = node->data.function.args[0];
        size_t m = arg->data.function.arg_count;
        if (m >= 3) {
            Vector2* poly = malloc(sizeof(Vector2) * m);
            size_t cnt = 0;
            for (size_t i = 0; i < m; i++) {
                double x, y;
                if (expr_point(arg->data.function.args[i], &x, &y)) poly[cnt++] = (Vector2){ (float)x, (float)(-y * state->yscale) };
            }
            if (cnt >= 3) DrawTriangleFan(poly, (int)cnt, state->color);
            free(poly);
        }
        return;
    }
    if (name == SYM_Text && n >= 2) {
        double x, y;
        if (!expr_point(node->data.function.args[1], &x, &y)) return;
        const Expr* content = node->data.function.args[0];
        char* owned = NULL;
        const char* label;
        if (content->type == EXPR_STRING) {
            label = content->data.string;
        } else {
            owned = expr_to_string((Expr*)content);
            label = owned ? owned : "";
        }
        float w = hershey_text_width(label, state->text_scale);
        hershey_draw_text(label, (float)x - w / 2.0f, (float)(-y * state->yscale), state->text_scale, 0.0f, state->color);
        if (owned) free(owned);
        return;
    }
}

/* ---------------- Axes ---------------- */

static double nice_step(double range, int target_ticks) {
    if (range <= 0) return 1.0;
    double raw = range / target_ticks;
    double mag = pow(10.0, floor(log10(raw)));
    double norm = raw / mag;
    double step;
    if (norm < 1.5) step = 1.0;
    else if (norm < 3.0) step = 2.0;
    else if (norm < 7.0) step = 5.0;
    else step = 10.0;
    return step * mag;
}

/* The portion of world-space currently visible through `camera` in a
 * `win_w` x `win_h` window -- axes are drawn against this (recomputed
 * every frame) rather than the data's original bounding box, so they
 * always reach the window edges and their tick spacing adapts as you
 * zoom, instead of staying pinned to a fixed world-space length. */
static PlotRange2D compute_visible_range(Camera2D camera, int win_w, int win_h) {
    Vector2 corners[4] = {
        GetScreenToWorld2D((Vector2){ 0, 0 }, camera),
        GetScreenToWorld2D((Vector2){ (float)win_w, 0 }, camera),
        GetScreenToWorld2D((Vector2){ 0, (float)win_h }, camera),
        GetScreenToWorld2D((Vector2){ (float)win_w, (float)win_h }, camera),
    };
    PlotRange2D r = { DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX };
    for (int i = 0; i < 4; i++) {
        double x = corners[i].x, y = -corners[i].y; /* undo the y-flip convention */
        if (x < r.xmin) r.xmin = x;
        if (x > r.xmax) r.xmax = x;
        if (y < r.ymin) r.ymin = y;
        if (y > r.ymax) r.ymax = y;
    }
    return r;
}

/* Resolve where the axes cross: the explicit AxesOrigin override (its
 * data-space y converted to render-space via ysc), or the default (0 if in
 * range, else clamp to the near edge) when unset. */
static void axes_origin(const PlotRange2D* range, double ysc, const GfxOptions* o,
                         double* ox, double* oy) {
    if (o->axes_origin_set) {
        *ox = o->axes_origin_x;
        *oy = o->axes_origin_y * ysc;
        return;
    }
    *ox = (range->xmin <= 0 && range->xmax >= 0) ? 0.0 : range->xmin;
    *oy = (range->ymin <= 0 && range->ymax >= 0) ? 0.0 : range->ymin;
}

/* World-space axis/tick lines -- call inside BeginMode2D so they pan and
 * zoom with the data. */
/* `range` is in render space (post y-scale). `ysc` maps data-y -> render-y;
 * y ticks are placed at nice *data* values (data*ysc), so labels read true.
 * `zoom` is the camera zoom: the strokes are drawn 1.5 px wide on screen (to
 * match the frame's weight) by expressing the world-space width as 1.5/zoom. */
static void draw_axes_lines(const PlotRange2D* range, double ysc, float zoom, const GfxOptions* o) {
    double ox, oy;
    axes_origin(range, ysc, o, &ox, &oy);
    float w = (zoom > 0.0f) ? 1.5f / zoom : 1.5f; /* world units -> 1.5 px on screen */
    Color col = to_raylib(o->axes_color);

    DrawLineEx((Vector2){ (float)range->xmin, (float)-oy }, (Vector2){ (float)range->xmax, (float)-oy }, w, col);
    DrawLineEx((Vector2){ (float)ox, (float)-range->ymin }, (Vector2){ (float)ox, (float)-range->ymax }, w, col);

    double xstep = nice_step(range->xmax - range->xmin, 8);
    double xtick = (range->ymax - range->ymin) * 0.015;
    double ytick = (range->xmax - range->xmin) * 0.015;

    for (double tx = ceil(range->xmin / xstep) * xstep; tx <= range->xmax + 1e-9; tx += xstep) {
        DrawLineEx((Vector2){ (float)tx, (float)-(oy - xtick) }, (Vector2){ (float)tx, (float)-(oy + xtick) }, w, col);
    }
    double dymin = range->ymin / ysc, dymax = range->ymax / ysc;
    double ystep = nice_step(dymax - dymin, 6);
    for (double ty = ceil(dymin / ystep) * ystep; ty <= dymax + 1e-9; ty += ystep) {
        double ry = ty * ysc;
        DrawLineEx((Vector2){ (float)(ox - ytick), (float)-ry }, (Vector2){ (float)(ox + ytick), (float)-ry }, w, col);
    }
}

/* Light grid lines at major tick positions (GridLines -> Automatic) or
 * explicit numeric positions (GridLines -> {xlist, ylist}). Call inside
 * BeginMode2D, before the axes/curve, so they sit underneath everything. */
static void draw_gridlines(const PlotRange2D* range, double ysc, float zoom, const GfxOptions* o) {
    if (!o->grid_x_on && !o->grid_y_on) return;
    float w = (zoom > 0.0f) ? 1.0f / zoom : 1.0f;
    Color col = to_raylib(o->grid_color);

    if (o->grid_x_on) {
        if (o->grid_x_list) {
            size_t n = o->grid_x_list->data.function.arg_count;
            for (size_t i = 0; i < n; i++) {
                double tx;
                if (expr_to_d(o->grid_x_list->data.function.args[i], &tx))
                    DrawLineEx((Vector2){ (float)tx, (float)-range->ymin }, (Vector2){ (float)tx, (float)-range->ymax }, w, col);
            }
        } else {
            double xstep = nice_step(range->xmax - range->xmin, 8);
            for (double tx = ceil(range->xmin / xstep) * xstep; tx <= range->xmax + 1e-9; tx += xstep)
                DrawLineEx((Vector2){ (float)tx, (float)-range->ymin }, (Vector2){ (float)tx, (float)-range->ymax }, w, col);
        }
    }
    if (o->grid_y_on) {
        if (o->grid_y_list) {
            size_t n = o->grid_y_list->data.function.arg_count;
            for (size_t i = 0; i < n; i++) {
                double ty;
                if (expr_to_d(o->grid_y_list->data.function.args[i], &ty)) {
                    double ry = ty * ysc;
                    DrawLineEx((Vector2){ (float)range->xmin, (float)-ry }, (Vector2){ (float)range->xmax, (float)-ry }, w, col);
                }
            }
        } else {
            double dymin = range->ymin / ysc, dymax = range->ymax / ysc;
            double ystep = nice_step(dymax - dymin, 6);
            for (double ty = ceil(dymin / ystep) * ystep; ty <= dymax + 1e-9; ty += ystep) {
                double ry = ty * ysc;
                DrawLineEx((Vector2){ (float)range->xmin, (float)-ry }, (Vector2){ (float)range->xmax, (float)-ry }, w, col);
            }
        }
    }
}

/* Screen-space tick labels -- call after EndMode2D with the same camera
 * so labels stay a crisp, fixed pixel size regardless of zoom.
 *
 * Each label is anchored to the *outer* end of its tick mark (not the axis
 * line) plus a small screen-pixel gap, so the digits sit clear of the ticks
 * instead of overprinting them. The Hershey baseline grows upward by
 * `cap = HERSHEY_CAP_HEIGHT * scale` px, which the offsets account for:
 * x-labels drop the baseline a full cap-height below the tick so the whole
 * glyph clears it; y-labels recentre vertically on the tick by half a cap. */
static void draw_axes_labels(const PlotRange2D* range, Camera2D camera, double ysc, const GfxOptions* o) {
    double ox, oy;
    axes_origin(range, ysc, o, &ox, &oy);
    Color col = to_raylib(o->ticks_color);
    double xstep = nice_step(range->xmax - range->xmin, 8);
    /* y ticks are stepped in data space; render position is ty*ysc below. */
    double dymin = range->ymin / ysc, dymax = range->ymax / ysc;
    double ystep = nice_step(dymax - dymin, 6);
    /* World-space tick half-lengths, matching draw_axes_lines exactly. */
    double xtick = (range->ymax - range->ymin) * 0.015;
    double ytick = (range->xmax - range->xmin) * 0.015;
    const float scale = 1.5f;
    const float cap = HERSHEY_CAP_HEIGHT * scale; /* glyph height in px */
    const float gap = 5.0f;                       /* clearance past the tick end */
    char buf[64];

    for (double tx = ceil(range->xmin / xstep) * xstep; tx <= range->xmax + 1e-9; tx += xstep) {
        if (fabs(tx) < 1e-9) tx = 0.0;
        snprintf(buf, sizeof(buf), "%g", tx);
        /* Screen position of the tick's lower (below-axis) end. */
        Vector2 tip = GetWorldToScreen2D((Vector2){ (float)tx, (float)-(oy - xtick) }, camera);
        float w = hershey_text_width(buf, scale);
        hershey_draw_text_ex(buf, tip.x - w / 2.0f, tip.y + gap + cap, scale, 0.0f, col, 1.5f);
    }
    for (double ty = ceil(dymin / ystep) * ystep; ty <= dymax + 1e-9; ty += ystep) {
        if (fabs(ty) < 1e-9) ty = 0.0;
        snprintf(buf, sizeof(buf), "%g", ty);
        /* Screen position of the tick's left (outside-axis) end. */
        Vector2 tip = GetWorldToScreen2D((Vector2){ (float)(ox - ytick), (float)(-ty * ysc) }, camera);
        float w = hershey_text_width(buf, scale);
        hershey_draw_text_ex(buf, tip.x - w - gap, tip.y + cap / 2.0f, scale, 0.0f, col, 1.5f);
    }
}

/* ---------------- Frame ----------------
 *
 * A Frame (Frame -> True) rules the plot with a rectangle along its edges,
 * replacing the through-the-origin cross of Axes. Unlike the axes, the frame
 * is a *fixed screen-space rectangle* inset from the window by a reserved
 * margin (see graphics_show): the data is fitted to and clipped against this
 * interior, so the curve never spills past the frame. The frame stays put
 * while the data pans/zooms within it, and the tick values -- read from the
 * data coordinates at the current frame edges -- update live, like a ruler.
 *
 * Ticks point inward; numeric labels for the major ticks sit just *outside*
 * the frame, in the reserved margin (bottom & left edges by default, falling
 * back to top & right). Major ticks land on the same "nice" values as the
 * axes; minor (sub-)ticks subdivide each major interval, the count chosen from
 * the leading digit of the major step so they always fall on round values: a
 * step of 1 splits into 5 (minors every 0.2), 2 into 4 (every 0.5), 5 into 5
 * (every 1). Minor ticks are half-length and unlabeled. */

/* Minor-tick subdivisions per major interval, from the step's leading digit.
 * Declared in render.h so the policy is unit-testable headless. */
int frame_minor_divs(double step) {
    if (step <= 0) return 5;
    double mag = pow(10.0, floor(log10(step)));
    double lead = step / mag;            /* ~1, 2 or 5 from nice_step */
    if (lead < 1.5) return 5;            /* 1 -> 0.2 minors */
    if (lead < 3.5) return 4;            /* 2 -> 0.5 minors */
    return 5;                            /* 5 -> 1   minors */
}

/* Draw the frame box, ticks and labels in screen space. The frame occupies the
 * fixed rectangle [rx, ry, rw, rh] (the reserved interior); `camera`/`ysc` map
 * data coordinates to that rectangle. Call after EndMode2D so the numbers stay
 * a crisp fixed size and the labels land in the margin outside the box. */
static void draw_frame(float rx, float ry, float rw, float rh,
                       Camera2D camera, double ysc, const GfxOptions* o) {
    Color col = to_raylib(o->frame_color);
    float L = rx, R = rx + rw, T = ry, B = ry + rh;
    /* 1.5x the 1px hairline default used by the axes, applied uniformly to the
     * frame box, its ticks and its labels. */
    const float lw = 1.5f;

    if (o->frame_edge[FR_BOTTOM]) DrawLineEx((Vector2){ L, B }, (Vector2){ R, B }, lw, col);
    if (o->frame_edge[FR_TOP])    DrawLineEx((Vector2){ L, T }, (Vector2){ R, T }, lw, col);
    if (o->frame_edge[FR_LEFT])   DrawLineEx((Vector2){ L, T }, (Vector2){ L, B }, lw, col);
    if (o->frame_edge[FR_RIGHT])  DrawLineEx((Vector2){ R, T }, (Vector2){ R, B }, lw, col);

    /* Data coordinates currently shown at the frame corners. wTL/wBR are in
     * draw space (y-down, negated render-y); undo that to get data-y. */
    Vector2 wTL = GetScreenToWorld2D((Vector2){ L, T }, camera);
    Vector2 wBR = GetScreenToWorld2D((Vector2){ R, B }, camera);
    double xL = wTL.x, xR = wBR.x;
    double dymax = (-wTL.y) / ysc, dymin = (-wBR.y) / ysc;
    if (!(xR > xL) || !(dymax > dymin) || !isfinite(xL) || !isfinite(dymax)) return;

    const float scale = 1.5f;
    const float cap = HERSHEY_CAP_HEIGHT * scale;
    const float gap = 5.0f;
    const float maj = 6.0f, minr = 3.0f; /* inward tick pixel lengths */
    char buf[64];

    /* Which edge carries labels: bottom/left preferred, else top/right. */
    int xlab = (o->frame_edge[FR_BOTTOM] && o->frame_ticks[FR_BOTTOM]) ? FR_BOTTOM
             : (o->frame_edge[FR_TOP] && o->frame_ticks[FR_TOP]) ? FR_TOP : -1;
    int ylab = (o->frame_edge[FR_LEFT] && o->frame_ticks[FR_LEFT]) ? FR_LEFT
             : (o->frame_edge[FR_RIGHT] && o->frame_ticks[FR_RIGHT]) ? FR_RIGHT : -1;

    /* X ticks along the bottom and/or top edges. */
    double xstep = nice_step(xR - xL, 8);
    int xsub = frame_minor_divs(xstep);
    double xm = xstep / xsub;
    if (xm > 0) {
        long i0 = (long)ceil(xL / xm - 1e-9), i1 = (long)floor(xR / xm + 1e-9);
        for (long i = i0; i <= i1; i++) {
            double tx = i * xm;
            bool major = (i % xsub == 0);
            float len = major ? maj : minr;
            float sx = GetWorldToScreen2D((Vector2){ (float)tx, 0.0f }, camera).x;
            if (sx < L - 0.5f || sx > R + 0.5f) continue;
            if (o->frame_edge[FR_BOTTOM] && o->frame_ticks[FR_BOTTOM])
                DrawLineEx((Vector2){ sx, B }, (Vector2){ sx, B - len }, lw, col);
            if (o->frame_edge[FR_TOP] && o->frame_ticks[FR_TOP])
                DrawLineEx((Vector2){ sx, T }, (Vector2){ sx, T + len }, lw, col);
            if (major && xlab >= 0) {
                snprintf(buf, sizeof(buf), "%g", fabs(tx) < 1e-9 ? 0.0 : tx);
                float w = hershey_text_width(buf, scale);
                /* Outside the frame: below the bottom edge, or above the top. */
                float baseline = (xlab == FR_BOTTOM) ? (B + gap + cap) : (T - gap);
                hershey_draw_text_ex(buf, sx - w / 2.0f, baseline, scale, 0.0f, col, lw);
            }
        }
    }

    /* Y ticks along the left and/or right edges. */
    double ystep = nice_step(dymax - dymin, 6);
    int ysub = frame_minor_divs(ystep);
    double ym = ystep / ysub;
    if (ym > 0) {
        long i0 = (long)ceil(dymin / ym - 1e-9), i1 = (long)floor(dymax / ym + 1e-9);
        for (long i = i0; i <= i1; i++) {
            double ty = i * ym;
            bool major = (i % ysub == 0);
            float len = major ? maj : minr;
            float sy = GetWorldToScreen2D((Vector2){ 0.0f, (float)(-ty * ysc) }, camera).y;
            if (sy < T - 0.5f || sy > B + 0.5f) continue;
            if (o->frame_edge[FR_LEFT] && o->frame_ticks[FR_LEFT])
                DrawLineEx((Vector2){ L, sy }, (Vector2){ L + len, sy }, lw, col);
            if (o->frame_edge[FR_RIGHT] && o->frame_ticks[FR_RIGHT])
                DrawLineEx((Vector2){ R, sy }, (Vector2){ R - len, sy }, lw, col);
            if (major && ylab >= 0) {
                snprintf(buf, sizeof(buf), "%g", fabs(ty) < 1e-9 ? 0.0 : ty);
                float w = hershey_text_width(buf, scale);
                /* Outside the frame: left of the left edge, or right of the right. */
                float tx = (ylab == FR_LEFT) ? (L - gap - w) : (R + gap);
                hershey_draw_text_ex(buf, tx, sy + cap / 2.0f, scale, 0.0f, col, lw);
            }
        }
    }
}

/* Frame-mode equivalent of AxesLabel: xlabel centered below the frame's
 * bottom edge, ylabel beside the left edge -- rotated 90 deg unless
 * RotateLabel -> False. Call alongside draw_frame, in screen space. */
static void draw_frame_label(float rx, float ry, float rw, float rh, const GfxOptions* o) {
    if (!o->frame_label || o->frame_label->type != EXPR_FUNCTION || o->frame_label->data.function.arg_count != 2) return;
    const Expr* xl = o->frame_label->data.function.args[0];
    const Expr* yl = o->frame_label->data.function.args[1];

    char* ox = NULL;
    const char* xlabel = (xl->type == EXPR_STRING) ? xl->data.string : (ox = expr_to_string((Expr*)xl));
    if (xlabel) {
        float w = hershey_text_width(xlabel, 1.8f);
        hershey_draw_text(xlabel, rx + rw / 2.0f - w / 2.0f, ry + rh + 40.0f, 1.8f, 0.0f, BLACK);
    }
    free(ox);

    char* oy = NULL;
    const char* ylabel = (yl->type == EXPR_STRING) ? yl->data.string : (oy = expr_to_string((Expr*)yl));
    if (ylabel) {
        if (o->rotate_label) {
            hershey_draw_text(ylabel, rx - 32.0f, ry + rh / 2.0f, 1.8f, 90.0f, BLACK);
        } else {
            float w = hershey_text_width(ylabel, 1.8f);
            hershey_draw_text(ylabel, rx - 16.0f - w, ry + rh / 2.0f, 1.8f, 0.0f, BLACK);
        }
    }
    free(oy);
}

static void draw_extra_labels(const GfxOptions* opts, int win_w, int win_h) {
    if (opts->plot_label) {
        char* owned = NULL;
        const char* label = NULL;
        if (opts->plot_label->type == EXPR_STRING) label = opts->plot_label->data.string;
        else { owned = expr_to_string((Expr*)opts->plot_label); label = owned; }
        if (label) {
            float w = hershey_text_width(label, 2.0f);
            hershey_draw_text(label, win_w / 2.0f - w / 2.0f, 28.0f, 2.0f, 0.0f, BLACK);
        }
        free(owned);
    }
    if (opts->axes_label && opts->axes_label->type == EXPR_FUNCTION && opts->axes_label->data.function.arg_count == 2) {
        const Expr* xl = opts->axes_label->data.function.args[0];
        const Expr* yl = opts->axes_label->data.function.args[1];

        char* ox = NULL;
        const char* xlabel = (xl->type == EXPR_STRING) ? xl->data.string : (ox = expr_to_string((Expr*)xl));
        if (xlabel) {
            float w = hershey_text_width(xlabel, 1.8f);
            hershey_draw_text(xlabel, win_w / 2.0f - w / 2.0f, win_h - 36.0f, 1.8f, 0.0f, BLACK);
        }
        free(ox);

        char* oy = NULL;
        const char* ylabel = (yl->type == EXPR_STRING) ? yl->data.string : (oy = expr_to_string((Expr*)yl));
        if (ylabel) hershey_draw_text(ylabel, 18.0f, (float)win_h / 2.0f, 1.8f, 90.0f, BLACK);
        free(oy);
    }
}

/* ---------------- Toolbar (Plotly-style modebar) ----------------
 *
 * A row of icon buttons in the top-right corner whose actions mirror
 * Plotly's 2D modebar. Everything here is screen-space UI chrome drawn
 * after EndMode2D; the icons are hand-drawn vector glyphs (no image
 * assets), in keeping with the Hershey vector font used elsewhere.
 *
 * Plotly's Box-Select / Lasso-Select are intentionally omitted: Mathilda
 * renders continuous primitives, not a discrete dataset, so there is
 * nothing to "select". */

typedef enum { TOOL_PAN = 0, TOOL_ZOOM = 1 } ToolMode;

typedef enum {
    TB_SAVE = 0, TB_ZOOMBOX, TB_PAN, TB_ZOOMIN, TB_ZOOMOUT,
    TB_AUTOSCALE, TB_RESET, TB_CLOSE, TB_COUNT
} ToolButton;

#define TB_BTN    30.0f   /* button edge (px)            */
#define TB_GAP     4.0f   /* gap between buttons (px)    */
#define TB_MARGIN 10.0f   /* inset from window top/right */

/* Screen rect of button i, laid left-to-right ending TB_MARGIN from the
 * window's right edge. */
static Rectangle tb_rect(int i, int win_w) {
    float total = TB_COUNT * TB_BTN + (TB_COUNT - 1) * TB_GAP;
    float x0 = (float)win_w - TB_MARGIN - total;
    return (Rectangle){ x0 + i * (TB_BTN + TB_GAP), TB_MARGIN, TB_BTN, TB_BTN };
}

static int tb_hit(Vector2 m, int win_w) {
    for (int i = 0; i < TB_COUNT; i++)
        if (CheckCollisionPointRec(m, tb_rect(i, win_w))) return i;
    return -1;
}

static const char* tb_tip(int k) {
    switch (k) {
        case TB_SAVE:      return "Download plot as PNG";
        case TB_ZOOMBOX:   return "Zoom (drag a box)";
        case TB_PAN:       return "Pan (drag to move)";
        case TB_ZOOMIN:    return "Zoom in";
        case TB_ZOOMOUT:   return "Zoom out";
        case TB_AUTOSCALE: return "Autoscale to fit";
        case TB_RESET:     return "Reset axes";
        case TB_CLOSE:     return "Close window";
        default:           return "";
    }
}

/* Keep the camera's zoom within the same bounds the scroll-wheel uses. */
static void clamp_zoom(Camera2D* c, float base) {
    if (c->zoom < base * 0.05f) c->zoom = base * 0.05f;
    if (c->zoom > base * 50.0f) c->zoom = base * 50.0f;
}

/* --- glyph drawing helpers --------------------------------------------- *
 * Raylib's DrawLineEx draws a bare quad with butt caps, so polylines show
 * notches where segments meet and bare ends look chopped. These helpers
 * give every glyph rounded caps (a disc at each endpoint closes the joints)
 * and solid arrowheads, for a clean modebar-style finish. */

#define TB_LW 2.2f   /* shared stroke width for all glyphs */

/* A line segment with rounded end-caps. */
static void stroke(Vector2 a, Vector2 b, Color c) {
    DrawLineEx(a, b, TB_LW, c);
    DrawCircleV(a, TB_LW * 0.5f, c);
    DrawCircleV(b, TB_LW * 0.5f, c);
}

/* A filled triangle, winding-agnostic (raylib back-face culls one order, so
 * drawing both guarantees it shows regardless of vertex orientation). */
static void tri(Vector2 a, Vector2 b, Vector2 c, Color col) {
    DrawTriangle(a, b, c, col);
    DrawTriangle(a, c, b, col);
}

/* A solid arrowhead whose point sits at `tip`, opening back along unit
 * direction `dir` (pointing toward the tip): `len` deep, `half` to a side. */
static void arrowhead(Vector2 tip, Vector2 dir, float len, float half, Color c) {
    Vector2 back = { tip.x - dir.x * len, tip.y - dir.y * len };
    Vector2 perp = { -dir.y, dir.x };
    Vector2 l = { back.x + perp.x * half, back.y + perp.y * half };
    Vector2 r = { back.x - perp.x * half, back.y - perp.y * half };
    tri(tip, l, r, c);
}

/* --- icon glyphs: each fills the inset content box `b` with colour `c` --- */

static void icon_save(Rectangle b, Color c) {            /* camera */
    DrawRectangleRoundedLinesEx((Rectangle){ b.x, b.y + b.height * 0.30f, b.width, b.height * 0.52f },
                                0.25f, 8, TB_LW, c);
    DrawRectangleRoundedLinesEx((Rectangle){ b.x + b.width * 0.32f, b.y + b.height * 0.12f,
                                             b.width * 0.30f, b.height * 0.20f }, 0.4f, 6, TB_LW, c);
    float lr = b.width * 0.16f;
    DrawRing((Vector2){ b.x + b.width * 0.5f, b.y + b.height * 0.56f },
             lr - TB_LW * 0.5f, lr + TB_LW * 0.5f, 0.0f, 360.0f, 32, c);
}

static void icon_magnifier(Rectangle b, Color c, int sign) { /* glass; sign: 0 plain, + in, - out */
    Vector2 ctr = { b.x + b.width * 0.40f, b.y + b.height * 0.40f };
    float r = b.width * 0.30f;
    /* a filled annulus, half the stroke width to a side -- same visual
     * weight as the strokes used by every other glyph */
    DrawRing(ctr, r - TB_LW * 0.5f, r + TB_LW * 0.5f, 0.0f, 360.0f, 48, c);
    stroke((Vector2){ ctr.x + r * 0.72f, ctr.y + r * 0.72f },
           (Vector2){ b.x + b.width * 0.96f, b.y + b.height * 0.96f }, c);
    if (sign != 0) {
        float s = r * 0.5f;
        stroke((Vector2){ ctr.x - s, ctr.y }, (Vector2){ ctr.x + s, ctr.y }, c);
        if (sign > 0) stroke((Vector2){ ctr.x, ctr.y - s }, (Vector2){ ctr.x, ctr.y + s }, c);
    }
}

static void icon_move(Rectangle b, Color c) {            /* four-way arrows */
    float cx = b.x + b.width * 0.5f, cy = b.y + b.height * 0.5f;
    float L = b.width * 0.48f;
    float hl = b.width * 0.22f, hw = b.width * 0.15f;   /* head length / half-width */
    /* shafts stop just short of each tip so the solid heads cap them */
    float s = L - hl * 0.5f;
    stroke((Vector2){ cx - s, cy }, (Vector2){ cx + s, cy }, c);
    stroke((Vector2){ cx, cy - s }, (Vector2){ cx, cy + s }, c);
    arrowhead((Vector2){ cx - L, cy }, (Vector2){ -1, 0 }, hl, hw, c);
    arrowhead((Vector2){ cx + L, cy }, (Vector2){ +1, 0 }, hl, hw, c);
    arrowhead((Vector2){ cx, cy - L }, (Vector2){ 0, -1 }, hl, hw, c);
    arrowhead((Vector2){ cx, cy + L }, (Vector2){ 0, +1 }, hl, hw, c);
}

static void icon_expand(Rectangle b, Color c) {          /* outward diagonal arrows */
    float cx = b.x + b.width * 0.5f, cy = b.y + b.height * 0.5f;
    float hl = b.width * 0.24f, hw = b.width * 0.15f;
    Vector2 corner[4] = {
        { b.x,           b.y           },
        { b.x + b.width, b.y           },
        { b.x,           b.y + b.height },
        { b.x + b.width, b.y + b.height },
    };
    for (int i = 0; i < 4; i++) {
        Vector2 e = corner[i];
        Vector2 d = { e.x - cx, e.y - cy };
        float len = sqrtf(d.x * d.x + d.y * d.y);
        if (len > 0) { d.x /= len; d.y /= len; }
        stroke((Vector2){ cx, cy },
               (Vector2){ e.x - d.x * hl * 0.5f, e.y - d.y * hl * 0.5f }, c);
        arrowhead(e, d, hl, hw, c);
    }
}

static void icon_home(Rectangle b, Color c) {            /* house */
    float cx = b.x + b.width * 0.5f;
    float roofY = b.y + b.height * 0.12f, eaveY = b.y + b.height * 0.45f, baseY = b.y + b.height * 0.88f;
    float lx = b.x + b.width * 0.16f, rx = b.x + b.width * 0.84f;
    stroke((Vector2){ lx, eaveY }, (Vector2){ cx, roofY }, c);
    stroke((Vector2){ cx, roofY }, (Vector2){ rx, eaveY }, c);
    float wlx = b.x + b.width * 0.26f, wrx = b.x + b.width * 0.74f;
    stroke((Vector2){ wlx, eaveY }, (Vector2){ wlx, baseY }, c);
    stroke((Vector2){ wrx, eaveY }, (Vector2){ wrx, baseY }, c);
    stroke((Vector2){ wlx, baseY }, (Vector2){ wrx, baseY }, c);
    float dlx = b.x + b.width * 0.44f, drx = b.x + b.width * 0.56f, dY = b.y + b.height * 0.62f;
    stroke((Vector2){ dlx, baseY }, (Vector2){ dlx, dY }, c);
    stroke((Vector2){ drx, baseY }, (Vector2){ drx, dY }, c);
    stroke((Vector2){ dlx, dY }, (Vector2){ drx, dY }, c);
}

static void icon_close(Rectangle b, Color c) {           /* an X */
    float m = b.width * 0.18f;
    stroke((Vector2){ b.x + m, b.y + m },
           (Vector2){ b.x + b.width - m, b.y + b.height - m }, c);
    stroke((Vector2){ b.x + b.width - m, b.y + m },
           (Vector2){ b.x + m, b.y + b.height - m }, c);
}

/* Draw the whole toolbar. `tool` highlights the active mode button;
 * `hover` (a button index or -1) gets a hover background + tooltip. */
static void draw_toolbar(int win_w, int tool, int hover) {
    float total = TB_COUNT * TB_BTN + (TB_COUNT - 1) * TB_GAP;
    Rectangle panel = { (float)win_w - TB_MARGIN - total - 5.0f, TB_MARGIN - 5.0f,
                        total + 10.0f, TB_BTN + 10.0f };
    DrawRectangleRounded(panel, 0.3f, 6, (Color){ 248, 248, 248, 225 });

    /* All glyphs share one gray; the active mode is shown by the button's
     * background tint, not by recolouring its icon. */
    const Color icol = { 90, 90, 90, 255 };
    for (int i = 0; i < TB_COUNT; i++) {
        Rectangle r = tb_rect(i, win_w);
        bool active = (i == TB_PAN && tool == TOOL_PAN) || (i == TB_ZOOMBOX && tool == TOOL_ZOOM);
        /* Close gets a red hover tint to flag its destructive action; every
         * other button shares the neutral blue-gray highlight. */
        if (i == hover) {
            Color hb = (i == TB_CLOSE) ? (Color){ 240, 205, 205, 255 } : (Color){ 220, 227, 236, 255 };
            DrawRectangleRounded(r, 0.3f, 6, hb);
        }
        else if (active)     DrawRectangleRounded(r, 0.3f, 6, (Color){ 205, 222, 245, 255 });
        Rectangle ic = { r.x + 6, r.y + 6, r.width - 12, r.height - 12 };
        switch (i) {
            case TB_SAVE:      icon_save(ic, icol);          break;
            case TB_ZOOMBOX:   icon_magnifier(ic, icol, 0);  break;
            case TB_PAN:       icon_move(ic, icol);          break;
            case TB_ZOOMIN:    icon_magnifier(ic, icol, +1); break;
            case TB_ZOOMOUT:   icon_magnifier(ic, icol, -1); break;
            case TB_AUTOSCALE: icon_expand(ic, icol);        break;
            case TB_RESET:     icon_home(ic, icol);          break;
            case TB_CLOSE:     icon_close(ic, (i == hover) ? (Color){ 190, 60, 60, 255 } : icol); break;
            default: break;
        }
    }

    if (hover >= 0) {
        const char* t = tb_tip(hover);
        int tw = MeasureText(t, 12);
        Rectangle r = tb_rect(hover, win_w);
        float tx = r.x + r.width * 0.5f - tw * 0.5f;
        if (tx + tw + 6 > win_w) tx = (float)win_w - tw - 6;
        if (tx < 4) tx = 4;
        float ty = r.y + r.height + 7;
        DrawRectangle((int)tx - 5, (int)ty - 3, tw + 10, 19, (Color){ 40, 40, 40, 235 });
        DrawText(t, (int)tx, (int)ty, 12, RAYWHITE);
    }
}

/* Finds the $PlotLegendData[{color,label}, ...] metadata arg on
 * graphics_expr's option list, if present (built by plot.c's
 * build_legend_meta when PlotLegends was given). Borrowed; NULL if absent. */
static const Expr* find_legend_data(const Expr* graphics_expr) {
    size_t argc = graphics_expr->data.function.arg_count;
    for (size_t i = 1; i < argc; i++) {
        const Expr* a = graphics_expr->data.function.args[i];
        if (a && a->type == EXPR_FUNCTION && a->data.function.head
            && a->data.function.head->type == EXPR_SYMBOL
            && a->data.function.head->data.symbol == SYM_PlotLegendData) return a;
    }
    return NULL;
}

/* Screen-space legend box: one swatch+label row per curve. Anchored
 * top-right, like the toolbar, but below it (the toolbar already owns
 * that corner) so the two never collide. */
static void draw_legend(const Expr* legend_data, int win_w) {
    size_t n = legend_data->data.function.arg_count;
    if (n == 0) return;

    const float row_h = 22.0f, swatch_w = 22.0f, pad = 8.0f, scale = 1.5f;
    const float top = TB_MARGIN + TB_BTN + 10.0f;

    float max_label_w = 0.0f;
    for (size_t i = 0; i < n; i++) {
        const Expr* entry = legend_data->data.function.args[i];
        if (entry->type != EXPR_FUNCTION || entry->data.function.arg_count != 2) continue;
        const Expr* label = entry->data.function.args[1];
        if (label->type != EXPR_STRING) continue;
        float w = hershey_text_width(label->data.string, scale);
        if (w > max_label_w) max_label_w = w;
    }
    float box_w = swatch_w + pad * 3 + max_label_w;
    float box_h = pad * 2 + n * row_h;
    float box_x = (float)win_w - TB_MARGIN - box_w;
    float box_y = top;

    DrawRectangleRec((Rectangle){ box_x, box_y, box_w, box_h }, (Color){ 255, 255, 255, 230 });
    DrawRectangleLinesEx((Rectangle){ box_x, box_y, box_w, box_h }, 1.0f, (Color){ 150, 150, 150, 255 });

    for (size_t i = 0; i < n; i++) {
        const Expr* entry = legend_data->data.function.args[i];
        if (entry->type != EXPR_FUNCTION || entry->data.function.arg_count != 2) continue;
        RGBA8 c = { 0, 0, 0, 255 };
        resolve_color(entry->data.function.args[0], &c);

        float ry = box_y + pad + (float)i * row_h;
        DrawRectangleRec((Rectangle){ box_x + pad, ry + 3.0f, swatch_w, row_h - 8.0f }, to_raylib(c));

        const Expr* label = entry->data.function.args[1];
        if (label->type == EXPR_STRING) {
            hershey_draw_text(label->data.string, box_x + pad * 2 + swatch_w, ry + row_h - 6.0f,
                               scale, 0.0f, BLACK);
        }
    }
}

/* ---------------- Main entry point ---------------- */

void graphics_show(const Expr* graphics_expr) {
    if (!graphics_expr || graphics_expr->type != EXPR_FUNCTION || graphics_expr->data.function.arg_count < 1) return;

    /* Escape hatch for headless test/CI runs: even when USE_GRAPHICS is
     * compiled in, a window must never pop up and block during automated
     * test execution. tests/test_graphics.c sets this before evaluating
     * anything. */
    const char* no_window = getenv("MATHILDA_NO_GRAPHICS_WINDOW");
    if (no_window && no_window[0] != '\0') {
        printf("Show: graphics window suppressed (MATHILDA_NO_GRAPHICS_WINDOW set).\n");
        return;
    }

    GfxOptions opts;
    gfx_options_parse(graphics_expr, &opts);
    const Expr* prims = graphics_expr->data.function.args[0];
    const Expr* legend_data = find_legend_data(graphics_expr);

    /* Live re-sampling: a Plot-produced Graphics carries its function so we
     * can re-sample at the current zoom (plot_resample), keeping curves like
     * Sin[1/x^2] smooth when magnified instead of exposing the home grid.
     * draw_prims is what the loop draws; it starts as the static `prims` and
     * is swapped for freshly sampled ones, which dyn_prims owns. */
    const Expr* draw_prims = prims;
    Expr* dyn_prims = NULL;
    bool resample_ok = true;          /* cleared if this isn't a Plot object */
    /* The x-span currently sampled (with margin) and the visible width at the
     * last re-sample -- the loop re-samples once the view leaves either. */
    double cov_lo = 0.0, cov_hi = 0.0, ref_vw = -1.0;

    PlotRange2D range = { DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX };
    if (!opts.x_auto && !opts.y_auto) {
        range = opts.range;
    } else {
        compute_bbox(prims, &range);
        if (range.xmin > range.xmax) { range.xmin = -1; range.xmax = 1; }
        if (range.ymin > range.ymax) { range.ymin = -1; range.ymax = 1; }
        /* Replace the raw y-extent with a spike-resistant band so a curve like
         * Tan[x] (whose adaptive samples climb to ~1e16 near its asymptotes)
         * frames its visible body instead of collapsing to a flat line. Gated
         * on prims_have_runaway: only curves with a near-vertical asymptote are
         * clipped, so smooth steep curves (x^3, Exp) keep their full, legitimate
         * extent. Skipped entirely under PlotRange -> All. */
        if (opts.y_auto && opts.clip_outliers && prims_have_runaway(prims)) {
            double* ys = NULL; size_t yn = 0, ycap = 0;
            gather_ys(prims, &ys, &yn, &ycap);
            if (yn > 0) {
                double lo, hi;
                plot_robust_yrange(ys, yn, &lo, &hi);
                if (lo <= hi) { range.ymin = lo; range.ymax = hi; }
            }
            free(ys);
        }
        double xpad = (range.xmax - range.xmin) * opts.pad_x_frac;
        double ypad = (range.ymax - range.ymin) * opts.pad_y_frac;
        /* A zero-width range still needs *some* padding to be visible -- but
         * only when padding is actually wanted (PlotRangePadding -> None
         * means exactly 0, always). */
        if (xpad <= 0 && opts.pad_x_frac > 0) xpad = 1.0;
        if (ypad <= 0 && opts.pad_y_frac > 0) ypad = 1.0;
        range.xmin -= xpad; range.xmax += xpad;
        range.ymin -= ypad; range.ymax += ypad;
        /* Override the explicitly fixed axis (e.g. PlotRange -> {ymin, ymax}
         * pins y while x remains data-driven). */
        if (!opts.x_auto) { range.xmin = opts.range.xmin; range.xmax = opts.range.xmax; }
        if (!opts.y_auto) { range.ymin = opts.range.ymin; range.ymax = opts.range.ymax; }
    }

    double data_w = range.xmax - range.xmin;
    double data_h = range.ymax - range.ymin;
    if (data_w <= 0) data_w = 1;
    if (data_h <= 0) data_h = 1;

    /* AspectRatio sets the *window's* height-to-width ratio, not merely the
     * internal y-scale: a wide-screen Plot (default 1/GoldenRatio) opens a
     * short, wide window with the curve filling it edge to edge, rather than
     * a fixed box with the curve letterboxed inside. See gfx_window_height. */
    opts.height = gfx_window_height(opts.width, opts.height, opts.aspect_ratio,
                                    opts.aspect_full, opts.height_pinned,
                                    data_w, data_h);

    /* 4x MSAA smooths every vector stroke -- the plot curves, the axes, and
     * the hand-drawn toolbar glyphs all read as crisp anti-aliased lines
     * rather than the stair-stepped default. Must precede InitWindow. */
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow((int)opts.width, (int)opts.height, "Mathilda");
    SetTargetFPS(60);

    /* Plot region: the rectangle the data is fitted to and clipped against.
     * Without a frame it is the whole window (unchanged behaviour). A frame
     * reserves a margin (~5% of each window dimension) around the region so the
     * box, its outward tick labels and the bottom help line sit outside the
     * data; the bottom/left margins are floored a little larger so multi-digit
     * labels and the help text always fit. The data fits *inside* the region,
     * so the curve never spills past the frame. */
    float reg_x = 0.0f, reg_y = 0.0f;
    float reg_w = (float)opts.width, reg_h = (float)opts.height;
    if (opts.frame) {
        /* FrameLabel adds one more text line outside the tick labels, so it
         * needs extra room past the tick-label minimums above. */
        float labelL = opts.frame_label ? 26.0f : 0.0f;
        float labelB = opts.frame_label ? 26.0f : 0.0f;
        float mL = (float)opts.width  * 0.05f; if (mL < 50.0f + labelL) mL = 50.0f + labelL;
        float mR = (float)opts.width  * 0.05f; if (mR < 20.0f) mR = 20.0f;
        float mT = (float)opts.height * 0.05f; if (mT < 20.0f) mT = 20.0f;
        float mB = (float)opts.height * 0.05f; if (mB < 48.0f + labelB) mB = 48.0f + labelB;
        /* Never let the margins swallow the whole window. */
        if (mL + mR < (float)opts.width  - 40.0f && mT + mB < (float)opts.height - 40.0f) {
            reg_x = mL; reg_y = mT;
            reg_w = (float)opts.width - mL - mR;
            reg_h = (float)opts.height - mT - mB;
        }
    }

    double aspect = opts.aspect_ratio > 0 ? opts.aspect_ratio : (data_h / data_w);
    double fit_by_width  = reg_w / data_w;
    double fit_by_height = reg_h / (data_w * aspect);
    float base_zoom = (float)(fit_by_width < fit_by_height ? fit_by_width : fit_by_height);
    if (base_zoom <= 0 || !isfinite(base_zoom)) base_zoom = 1.0f;

    /* Non-uniform vertical scale: stretch data-y into "render space" so the
     * x-range fills the width and the y-range fills its aspect-determined share
     * of the height independently. With AspectRatio Automatic this collapses to
     * 1 (true geometry); a forced AspectRatio (e.g. Plot's 1/GoldenRatio) makes
     * PlotRange -> {ymin,ymax} actually frame the requested band, which a single
     * uniform zoom cannot. The camera, mouse, and re-sampler all live in render
     * space; only primitive y-coords and y-tick labels convert. */
    double ysc = aspect * data_w / data_h;
    if (!isfinite(ysc) || ysc <= 0) ysc = 1.0;

    Camera2D camera = { 0 };
    /* Anchor the camera on the region centre (window centre when unframed) so
     * the data fills the region rather than the whole window. */
    camera.offset = (Vector2){ reg_x + reg_w / 2.0f, reg_y + reg_h / 2.0f };
    camera.target = (Vector2){ (float)((range.xmin + range.xmax) / 2.0), (float)(-(range.ymin + range.ymax) / 2.0 * ysc) };
    camera.rotation = 0.0f;
    camera.zoom = base_zoom;
    const Camera2D home = camera;

    DrawState init_state;
    init_state.color = to_raylib(opts.style_color);
    /* Default stroke renders 1.5 px wide at the home zoom (1.5x the former
     * 1-px hairline). Expressed in world units (= screen px / zoom) so the
     * line scales with zoom like point_size, rather than staying pinned to
     * a single pixel. base_zoom is guaranteed positive and finite above. */
    init_state.thickness = 1.5f / base_zoom;
    init_state.point_size = (float)(fmax(data_w, data_h) * 0.006);
    init_state.text_scale = (float)(fmax(data_w, data_h) * 0.03 / HERSHEY_CAP_HEIGHT);
    init_state.yscale = (float)ysc;

    int tool = TOOL_PAN;          /* active left-drag tool (toolbar-selected) */
    bool selecting = false;       /* mid box-zoom drag */
    bool left_drag_canvas = false;/* current left drag began on the canvas, not the toolbar */
    Vector2 sel_start = { 0, 0 };
    bool shot = false;            /* this frame renders chrome-free for a capture */
    int toast = 0;                /* frames left to flash the "saved" confirmation */

    /* Seed the sampled span from the curve's actual x-extent. */
    {
        PlotRange2D dbb = { DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX };
        compute_bbox(prims, &dbb);
        if (dbb.xmin < dbb.xmax) { cov_lo = dbb.xmin; cov_hi = dbb.xmax; }
        else { cov_lo = range.xmin; cov_hi = range.xmax; }
    }

    while (!WindowShouldClose()) {
        Vector2 mouse = GetMousePosition();
        int hover = tb_hit(mouse, (int)opts.width);

        /* Scroll zoom, except while pointing at the toolbar. Zoom *about the
         * cursor*: the world point under the mouse stays fixed, so scrolling
         * magnifies wherever the user points instead of always diving toward
         * the plot centre. This matters for any plot whose centre is the worst
         * place to land -- e.g. Plot[Sin[1/x^2], {x, -2 Pi, 2 Pi}], where the
         * centre is the x=0 essential singularity that no sampling can resolve;
         * cursor-anchored zoom lets the user reach the resolvable detail the
         * live re-sampler then refines. */
        float wheel = GetMouseWheelMove();
        if (wheel != 0.0f && hover < 0) {
            Vector2 before = GetScreenToWorld2D(mouse, camera);
            camera.zoom *= (wheel > 0) ? 1.1f : (1.0f / 1.1f);
            clamp_zoom(&camera, base_zoom);
            Vector2 after = GetScreenToWorld2D(mouse, camera);
            camera.target.x += before.x - after.x;
            camera.target.y += before.y - after.y;
        }

        /* Left-press: a toolbar button fires its action; otherwise the press
         * begins a canvas gesture (pan or box-zoom, per the active tool). */
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            if (hover >= 0) {
                left_drag_canvas = false;
                switch (hover) {
                    case TB_SAVE:      shot = true; break;
                    case TB_ZOOMBOX:   tool = TOOL_ZOOM; break;
                    case TB_PAN:       tool = TOOL_PAN; break;
                    case TB_ZOOMIN:    camera.zoom *= 1.3f;          clamp_zoom(&camera, base_zoom); break;
                    case TB_ZOOMOUT:   camera.zoom *= (1.0f / 1.3f); clamp_zoom(&camera, base_zoom); break;
                    case TB_AUTOSCALE: camera = home; break;
                    case TB_RESET:     camera = home; break;
                    case TB_CLOSE:     goto close_window;
                    default: break;
                }
            } else {
                left_drag_canvas = true;
                if (tool == TOOL_ZOOM) { selecting = true; sel_start = mouse; }
            }
        }

        /* Left-release: finish a box-zoom by fitting the camera to the drawn
         * rectangle (ignored if it was really just a click). */
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) {
            if (selecting) {
                if (fabsf(mouse.x - sel_start.x) > 4.0f && fabsf(mouse.y - sel_start.y) > 4.0f) {
                    Vector2 a = GetScreenToWorld2D(sel_start, camera);
                    Vector2 b = GetScreenToWorld2D(mouse, camera);
                    float rw = fabsf(b.x - a.x), rh = fabsf(b.y - a.y);
                    if (rw > 1e-9f && rh > 1e-9f) {
                        camera.target = (Vector2){ (a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f };
                        float zx = reg_w / rw, zy = reg_h / rh;
                        camera.zoom = zx < zy ? zx : zy;
                        clamp_zoom(&camera, base_zoom);
                    }
                }
                selecting = false;
            }
            left_drag_canvas = false;
        }

        /* Pan: right/middle drag always; left drag only in Pan mode and only
         * when the drag started on the canvas. */
        bool panning = IsMouseButtonDown(MOUSE_BUTTON_RIGHT) || IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)
                    || (tool == TOOL_PAN && left_drag_canvas && IsMouseButtonDown(MOUSE_BUTTON_LEFT));
        if (panning) {
            Vector2 delta = GetMouseDelta();
            float ang = -camera.rotation * (float)M_PI / 180.0f;
            float dx = delta.x * cosf(ang) - delta.y * sinf(ang);
            float dy = delta.x * sinf(ang) + delta.y * cosf(ang);
            camera.target.x -= dx / camera.zoom;
            camera.target.y -= dy / camera.zoom;
        }

        if (IsKeyDown(KEY_Q)) camera.rotation -= 60.0f * GetFrameTime();
        if (IsKeyDown(KEY_E)) camera.rotation += 60.0f * GetFrameTime();
        if (IsKeyPressed(KEY_R)) camera = home;
        if (IsKeyPressed(KEY_ESCAPE)) break;

        PlotRange2D visible = compute_visible_range(camera, (int)opts.width, (int)opts.height);

        /* Re-sample the curve for the current x-window when the view has
         * moved enough -- zoomed past ~30% (resolution change) or panned to
         * the edge of what we last sampled (new territory). Held off while a
         * drag is in flight so gestures stay smooth: the resample lands on
         * release; discrete scroll/zoom-button steps resample next frame.
         * plot_resample returns NULL for a non-Plot Graphics, after which we
         * stop trying and keep the static primitives. */
        if (ref_vw < 0.0) ref_vw = visible.xmax - visible.xmin;
        if (resample_ok
            && !IsMouseButtonDown(MOUSE_BUTTON_LEFT)
            && !IsMouseButtonDown(MOUSE_BUTTON_RIGHT)
            && !IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
            double vw = visible.xmax - visible.xmin;
            bool zoom_changed = vw < ref_vw * 0.7 || vw > ref_vw * 1.4;
            double edge = vw * 0.02;
            bool panned_off = visible.xmin < cov_lo + edge || visible.xmax > cov_hi - edge;
            if (vw > 0.0 && isfinite(vw) && (zoom_changed || panned_off)) {
                double margin = vw * 0.25;
                double nx0 = visible.xmin - margin, nx1 = visible.xmax + margin;
                /* Clip refinement to the on-screen y-band so a curve that
                 * dives out of frame doesn't waste detail off-screen. */
                Expr* np = plot_resample(graphics_expr, nx0, nx1, visible.ymin, visible.ymax);
                if (np) {
                    if (dyn_prims) expr_free(dyn_prims);
                    dyn_prims = np;
                    draw_prims = dyn_prims;
                    cov_lo = nx0; cov_hi = nx1; ref_vw = vw;
                } else {
                    resample_ok = false;
                }
            }
        }

        BeginDrawing();
        ClearBackground(to_raylib(opts.background));

        /* With a frame, clip the data (and axes) to the region so nothing
         * spills past the box; the frame, labels and chrome draw unclipped. */
        if (opts.frame) BeginScissorMode((int)reg_x, (int)reg_y, (int)reg_w, (int)reg_h);
        BeginMode2D(camera);
        if (opts.prolog) { DrawState ps = init_state; draw_primitive(opts.prolog, &ps); }
        draw_gridlines(&visible, ysc, camera.zoom, &opts);
        if (opts.axes) draw_axes_lines(&visible, ysc, camera.zoom, &opts);
        DrawState state = init_state;
        draw_primitive(draw_prims, &state);
        if (opts.epilog) { DrawState es = init_state; draw_primitive(opts.epilog, &es); }
        EndMode2D();
        if (opts.frame) EndScissorMode();

        if (opts.axes) draw_axes_labels(&visible, camera, ysc, &opts);
        if (opts.frame) draw_frame(reg_x, reg_y, reg_w, reg_h, camera, ysc, &opts);
        if (opts.frame) draw_frame_label(reg_x, reg_y, reg_w, reg_h, &opts);
        draw_extra_labels(&opts, (int)opts.width, (int)opts.height);
        if (legend_data) draw_legend(legend_data, (int)opts.width);

        /* On a capture frame suppress every bit of UI chrome so the saved
         * PNG holds only the plot; the capture happens just after the swap. */
        if (!shot) {
            if (selecting) {
                Rectangle sel = { fminf(sel_start.x, mouse.x), fminf(sel_start.y, mouse.y),
                                  fabsf(mouse.x - sel_start.x), fabsf(mouse.y - sel_start.y) };
                DrawRectangleRec(sel, (Color){ 30, 80, 180, 40 });
                DrawRectangleLinesEx(sel, 1.0f, (Color){ 30, 80, 180, 180 });
            }
            draw_toolbar((int)opts.width, tool, hover);
            DrawText("drag: pan/zoom per tool   scroll: zoom   right-drag: pan   Q/E: rotate   R: reset   Esc: close",
                     10, (int)opts.height - 22, 14, GRAY);
            if (toast > 0) {
                const char* msg = "Saved mathilda_plot.png";
                int tw = MeasureText(msg, 16);
                DrawRectangle(10, 10, tw + 16, 26, (Color){ 40, 40, 40, 220 });
                DrawText(msg, 18, 15, 16, RAYWHITE);
                toast--;
            }
        }

        EndDrawing();

        /* EndDrawing has swapped buffers, so the chrome-free frame is now
         * presented -- capture it here (raylib's own F12 capture uses this
         * exact point), then clear the flag and flash a confirmation. */
        if (shot) {
            TakeScreenshot("mathilda_plot.png");
            shot = false;
            toast = 120;
        }
    }

close_window:
    if (dyn_prims) expr_free(dyn_prims);
    CloseWindow();
}
