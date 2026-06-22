/* render.c — Raylib backend: graphics_show() (see show.h).
 *
 * Coordinate convention used throughout this file: math-space is y-up;
 * Raylib draw-space is y-down. Every world point's y is negated exactly
 * once, right where it's converted to a Vector2 for a draw call -- there
 * are no other sign flips anywhere else in this file. */

#include "show.h"
#include "hershey_font.h"
#include "primitives.h"
#include "sym_names.h"
#include "print.h"
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

static Color to_raylib(RGBA8 c) { Color out = { c.r, c.g, c.b, c.a }; return out; }

/* ---------------- Option parsing ---------------- */

typedef struct {
    bool axes;
    bool range_auto;
    PlotRange2D range;
    double aspect_ratio; /* <= 0 means Automatic */
    RGBA8 style_color;
    RGBA8 background;
    long width, height;
    const Expr* axes_label; /* borrowed */
    const Expr* plot_label; /* borrowed */
} GfxOptions;

static void gfx_options_parse(const Expr* graphics, GfxOptions* o) {
    o->axes = false;
    o->range_auto = true;
    o->aspect_ratio = -1.0;
    o->style_color = (RGBA8){ 30, 80, 180, 255 };
    o->background = (RGBA8){ 255, 255, 255, 255 };
    o->width = 800;
    o->height = 600;
    o->axes_label = NULL;
    o->plot_label = NULL;

    size_t argc = graphics->data.function.arg_count;
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
        } else if (name == SYM_AspectRatio) {
            double v;
            if (expr_to_d(rhs, &v) && v > 0) o->aspect_ratio = v;
        } else if (name == SYM_PlotStyle && rhs->type == EXPR_FUNCTION
                   && rhs->data.function.head->type == EXPR_SYMBOL) {
            if (rhs->data.function.head->data.symbol == SYM_RGBColor) o->style_color = rgba_from_rgbcolor(rhs);
            else if (rhs->data.function.head->data.symbol == SYM_GrayLevel) o->style_color = rgba_from_graylevel(rhs);
        } else if (name == SYM_Background && rhs->type == EXPR_FUNCTION
                   && rhs->data.function.head->type == EXPR_SYMBOL) {
            if (rhs->data.function.head->data.symbol == SYM_RGBColor) o->background = rgba_from_rgbcolor(rhs);
            else if (rhs->data.function.head->data.symbol == SYM_GrayLevel) o->background = rgba_from_graylevel(rhs);
        } else if (name == SYM_ImageSize) {
            if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL
                && rhs->data.function.head->data.symbol == SYM_List && rhs->data.function.arg_count == 2) {
                double w, hh;
                if (expr_to_d(rhs->data.function.args[0], &w)) o->width = (long)w;
                if (expr_to_d(rhs->data.function.args[1], &hh)) o->height = (long)hh;
            } else {
                double s;
                if (expr_to_d(rhs, &s) && s > 0) { o->width = (long)s; o->height = (long)(s * 0.75); }
            }
        } else if (name == SYM_PlotRange) {
            if (rhs->type == EXPR_FUNCTION && rhs->data.function.arg_count == 2) {
                const Expr* xr = rhs->data.function.args[0];
                const Expr* yr = rhs->data.function.args[1];
                double xmin, xmax, ymin, ymax;
                if (xr->type == EXPR_FUNCTION && xr->data.function.arg_count == 2
                    && yr->type == EXPR_FUNCTION && yr->data.function.arg_count == 2
                    && expr_to_d(xr->data.function.args[0], &xmin) && expr_to_d(xr->data.function.args[1], &xmax)
                    && expr_to_d(yr->data.function.args[0], &ymin) && expr_to_d(yr->data.function.args[1], &ymax)) {
                    o->range_auto = false;
                    o->range.xmin = xmin; o->range.xmax = xmax;
                    o->range.ymin = ymin; o->range.ymax = ymax;
                }
            }
        } else if (name == SYM_AxesLabel) {
            o->axes_label = rhs;
        } else if (name == SYM_PlotLabel) {
            o->plot_label = rhs;
        }
    }

    if (o->width < 100) o->width = 100;
    if (o->height < 100) o->height = 100;
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
    } else if ((name == SYM_Circle || name == SYM_Disk) && n >= 2) {
        double cx, cy, r;
        if (expr_point(node->data.function.args[0], &cx, &cy) && expr_to_d(node->data.function.args[1], &r)) {
            update_bbox(bb, cx - r, cy - r);
            update_bbox(bb, cx + r, cy + r);
        }
    } else if (name == SYM_Text && n >= 2) {
        compute_bbox(node->data.function.args[1], bb);
    }
}

/* ---------------- Drawing ---------------- */

typedef struct {
    Color color;
    float thickness;  /* world units; <= 0 means hairline */
    float point_size; /* world units (radius) */
    float text_scale; /* world units per Hershey grid unit, for Text[] */
} DrawState;

static void draw_polyline(const Expr* pts_list, const DrawState* state) {
    size_t n = pts_list->data.function.arg_count;
    Vector2 prev = { 0, 0 };
    bool have_prev = false;
    for (size_t i = 0; i < n; i++) {
        double x, y;
        if (!expr_point(pts_list->data.function.args[i], &x, &y)) { have_prev = false; continue; }
        Vector2 cur = { (float)x, (float)-y };
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
    if (name == SYM_RGBColor) { state->color = to_raylib(rgba_from_rgbcolor(node)); return; }
    if (name == SYM_GrayLevel) { state->color = to_raylib(rgba_from_graylevel(node)); return; }
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
            DrawCircleV((Vector2){ (float)x, (float)-y }, state->point_size, state->color);
        } else if (arg->type == EXPR_FUNCTION) {
            for (size_t i = 0; i < arg->data.function.arg_count; i++) {
                double px, py;
                if (expr_point(arg->data.function.args[i], &px, &py))
                    DrawCircleV((Vector2){ (float)px, (float)-py }, state->point_size, state->color);
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
            Rectangle r = { (float)xmin, (float)-ymax, (float)(xmax - xmin), (float)(ymax - ymin) };
            DrawRectangleRec(r, state->color);
        }
        return;
    }
    if ((name == SYM_Circle || name == SYM_Disk) && n >= 2) {
        double cx, cy, r;
        if (expr_point(node->data.function.args[0], &cx, &cy) && expr_to_d(node->data.function.args[1], &r)) {
            Vector2 center = { (float)cx, (float)-cy };
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
                if (expr_point(arg->data.function.args[i], &x, &y)) poly[cnt++] = (Vector2){ (float)x, (float)-y };
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
        hershey_draw_text(label, (float)x - w / 2.0f, (float)-y, state->text_scale, 0.0f, state->color);
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

/* World-space axis/tick lines -- call inside BeginMode2D so they pan and
 * zoom with the data. */
static void draw_axes_lines(const PlotRange2D* range) {
    double ox = (range->xmin <= 0 && range->xmax >= 0) ? 0.0 : range->xmin;
    double oy = (range->ymin <= 0 && range->ymax >= 0) ? 0.0 : range->ymin;

    DrawLineV((Vector2){ (float)range->xmin, (float)-oy }, (Vector2){ (float)range->xmax, (float)-oy }, DARKGRAY);
    DrawLineV((Vector2){ (float)ox, (float)-range->ymin }, (Vector2){ (float)ox, (float)-range->ymax }, DARKGRAY);

    double xstep = nice_step(range->xmax - range->xmin, 8);
    double ystep = nice_step(range->ymax - range->ymin, 6);
    double xtick = (range->ymax - range->ymin) * 0.015;
    double ytick = (range->xmax - range->xmin) * 0.015;

    for (double tx = ceil(range->xmin / xstep) * xstep; tx <= range->xmax + 1e-9; tx += xstep) {
        DrawLineV((Vector2){ (float)tx, (float)-(oy - xtick) }, (Vector2){ (float)tx, (float)-(oy + xtick) }, DARKGRAY);
    }
    for (double ty = ceil(range->ymin / ystep) * ystep; ty <= range->ymax + 1e-9; ty += ystep) {
        DrawLineV((Vector2){ (float)(ox - ytick), (float)-ty }, (Vector2){ (float)(ox + ytick), (float)-ty }, DARKGRAY);
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
static void draw_axes_labels(const PlotRange2D* range, Camera2D camera) {
    double ox = (range->xmin <= 0 && range->xmax >= 0) ? 0.0 : range->xmin;
    double oy = (range->ymin <= 0 && range->ymax >= 0) ? 0.0 : range->ymin;
    double xstep = nice_step(range->xmax - range->xmin, 8);
    double ystep = nice_step(range->ymax - range->ymin, 6);
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
        hershey_draw_text(buf, tip.x - w / 2.0f, tip.y + gap + cap, scale, 0.0f, DARKGRAY);
    }
    for (double ty = ceil(range->ymin / ystep) * ystep; ty <= range->ymax + 1e-9; ty += ystep) {
        if (fabs(ty) < 1e-9) ty = 0.0;
        snprintf(buf, sizeof(buf), "%g", ty);
        /* Screen position of the tick's left (outside-axis) end. */
        Vector2 tip = GetWorldToScreen2D((Vector2){ (float)(ox - ytick), (float)-ty }, camera);
        float w = hershey_text_width(buf, scale);
        hershey_draw_text(buf, tip.x - w - gap, tip.y + cap / 2.0f, scale, 0.0f, DARKGRAY);
    }
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

    PlotRange2D range = { DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX };
    if (!opts.range_auto) {
        range = opts.range;
    } else {
        compute_bbox(prims, &range);
        if (range.xmin > range.xmax) { range.xmin = -1; range.xmax = 1; }
        if (range.ymin > range.ymax) { range.ymin = -1; range.ymax = 1; }
        double xpad = (range.xmax - range.xmin) * 0.08;
        double ypad = (range.ymax - range.ymin) * 0.08;
        if (xpad <= 0) xpad = 1.0;
        if (ypad <= 0) ypad = 1.0;
        range.xmin -= xpad; range.xmax += xpad;
        range.ymin -= ypad; range.ymax += ypad;
    }

    double data_w = range.xmax - range.xmin;
    double data_h = range.ymax - range.ymin;
    if (data_w <= 0) data_w = 1;
    if (data_h <= 0) data_h = 1;

    InitWindow((int)opts.width, (int)opts.height, "Mathilda");
    SetTargetFPS(60);

    double aspect = opts.aspect_ratio > 0 ? opts.aspect_ratio : (data_h / data_w);
    double fit_by_width  = opts.width / data_w;
    double fit_by_height = opts.height / (data_w * aspect);
    float base_zoom = (float)(fit_by_width < fit_by_height ? fit_by_width : fit_by_height);
    if (base_zoom <= 0 || !isfinite(base_zoom)) base_zoom = 1.0f;

    Camera2D camera = { 0 };
    camera.offset = (Vector2){ opts.width / 2.0f, opts.height / 2.0f };
    camera.target = (Vector2){ (float)((range.xmin + range.xmax) / 2.0), (float)(-(range.ymin + range.ymax) / 2.0) };
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

    while (!WindowShouldClose()) {
        float wheel = GetMouseWheelMove();
        if (wheel != 0.0f) {
            camera.zoom *= (wheel > 0) ? 1.1f : (1.0f / 1.1f);
            if (camera.zoom < base_zoom * 0.05f) camera.zoom = base_zoom * 0.05f;
            if (camera.zoom > base_zoom * 50.0f) camera.zoom = base_zoom * 50.0f;
        }
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) || IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
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

        BeginDrawing();
        ClearBackground(to_raylib(opts.background));

        BeginMode2D(camera);
        if (opts.axes) draw_axes_lines(&visible);
        DrawState state = init_state;
        draw_primitive(prims, &state);
        EndMode2D();

        if (opts.axes) draw_axes_labels(&visible, camera);
        draw_extra_labels(&opts, (int)opts.width, (int)opts.height);

        DrawText("scroll: zoom   right-drag: pan   Q/E: rotate   R: reset view   Esc: close",
                 10, (int)opts.height - 22, 14, GRAY);

        EndDrawing();
    }

    CloseWindow();
}
