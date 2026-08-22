/* graphics_export.c -- vector PDF export for 2D Graphics[...] expressions.
 *
 * This is the dependency-free, headless half of graphics export (the raster
 * half lives in render.c behind USE_GRAPHICS). It walks the same primitive
 * tree the Raylib renderer and graphics_json.c walk -- see the shape summary
 * in graphics_json.c -- and emits a single-page PDF whose content stream maps
 * one-to-one onto the primitives:
 *
 *     Line       -> polyline stroke (m/l/S)
 *     Polygon    -> filled path      (m/l/h/f)
 *     Disk       -> filled circle    (4 Beziers, f)
 *     Circle     -> stroked circle   (4 Beziers, S)
 *     Point      -> filled dot
 *     Rectangle  -> filled + stroked box
 *     Arrow      -> polyline + solid arrowhead
 *     Text       -> Helvetica string (BT .. Tj .. ET)
 *
 * plus the RGBColor/GrayLevel/Hue/CMYKColor/Opacity/Thickness/PointSize
 * directives. PDF's coordinate system is y-up with the origin at the lower
 * left, which is exactly the mathematical convention, so world y needs no
 * flip. Colour, opacity and the "nice" tick policy mirror the renderer so the
 * PDF and the on-screen/PNG output agree.
 *
 * The writer needs no third-party code: text uses the PDF base-14 Helvetica
 * font (present in every conforming viewer), so no font file is embedded.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "expr.h"
#include "sym_names.h"
#include "graphics_export.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Bezier "magic" constant for approximating a quarter circle. */
#define KAPPA 0.5522847498307936

/* ------------------------------------------------------------------ Buf --- */
/* A growable byte buffer. PDF content is ASCII, so a plain char buffer is
 * enough; we never embed binary. */
typedef struct { char* p; size_t len, cap; } Buf;

static void buf_init(Buf* b) { b->p = NULL; b->len = 0; b->cap = 0; }
static void buf_free(Buf* b) { free(b->p); b->p = NULL; b->len = b->cap = 0; }

static int buf_reserve(Buf* b, size_t extra) {
    if (b->len + extra + 1 <= b->cap) return 1;
    size_t ncap = b->cap ? b->cap * 2 : 1024;
    while (ncap < b->len + extra + 1) ncap *= 2;
    char* np = (char*)realloc(b->p, ncap);
    if (!np) return 0;
    b->p = np; b->cap = ncap;
    return 1;
}

static void buf_cat(Buf* b, const char* s) {
    size_t n = strlen(s);
    if (!buf_reserve(b, n)) return;
    memcpy(b->p + b->len, s, n);
    b->len += n; b->p[b->len] = '\0';
}

/* printf into the buffer (bounded scratch; content lines are short). */
static void buf_catf(Buf* b, const char* fmt, ...) {
    char tmp[512];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(tmp, sizeof(tmp), fmt, ap);
    va_end(ap);
    buf_cat(b, tmp);
}

/* ------------------------------------------------------- Expr helpers ----- */

static int head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

static int to_double(const Expr* e, double* out) {
    if (!e) return 0;
    if (e->type == EXPR_REAL)    { *out = e->data.real;              return 1; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer;   return 1; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return 1; }
    return 0;
}

static double clip01(double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); }

/* Resolve a colour directive to RGBA in [0,1]. Mirrors graphics_json.c's
 * resolve_color_rgb so PDF colours match the renderer for all four forms. */
static int resolve_color(const Expr* e, double* r, double* g, double* b, double* a) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) return 0;
    const char*  h    = e->data.function.head->data.symbol.name;
    const Expr** args = (const Expr**)e->data.function.args;
    size_t       n    = e->data.function.arg_count;
    double rr = 0, gg = 0, bb = 0, aa = 1;

    if (h == SYM_RGBColor) {
        if (n < 3) return 0;
        if (!to_double(args[0], &rr) || !to_double(args[1], &gg) || !to_double(args[2], &bb))
            return 0;
        if (n >= 4) to_double(args[3], &aa);
    } else if (h == SYM_GrayLevel) {
        if (n < 1 || !to_double(args[0], &rr)) return 0;
        gg = bb = rr;
        if (n >= 2) to_double(args[1], &aa);
    } else if (h == SYM_Hue) {
        double hh = 0, s = 1, v = 1;
        if (n >= 1) to_double(args[0], &hh);
        if (n >= 3) { to_double(args[1], &s); to_double(args[2], &v); }
        if (n >= 4) to_double(args[3], &aa);
        hh -= floor(hh); s = clip01(s); v = clip01(v);
        if (s <= 0.0) { rr = gg = bb = v; }
        else {
            double h6 = hh * 6.0; int qi = (int)floor(h6); double f = h6 - qi;
            double p = v*(1-s), q = v*(1-f*s), t = v*(1-(1-f)*s);
            switch (((qi % 6) + 6) % 6) {
                case 0: rr=v;  gg=t;  bb=p;  break;
                case 1: rr=q;  gg=v;  bb=p;  break;
                case 2: rr=p;  gg=v;  bb=t;  break;
                case 3: rr=p;  gg=q;  bb=v;  break;
                case 4: rr=t;  gg=p;  bb=v;  break;
                default: rr=v; gg=p;  bb=q;  break;
            }
        }
    } else if (h == SYM_CMYKColor) {
        const Expr** ca = args; size_t cn = n;
        if (cn == 1 && head_is(args[0], SYM_List)) {
            cn = args[0]->data.function.arg_count;
            ca = (const Expr**)args[0]->data.function.args;
        }
        double c=0,m=0,y=0,k=0;
        if (cn >= 3) { to_double(ca[0],&c); to_double(ca[1],&m); to_double(ca[2],&y); }
        if (cn >= 4) to_double(ca[3],&k);
        if (cn >= 5) to_double(ca[4],&aa);
        double w = 1.0 - clip01(k);
        rr = (1-clip01(c))*w; gg = (1-clip01(m))*w; bb = (1-clip01(y))*w;
    } else {
        return 0;
    }
    *r = clip01(rr); *g = clip01(gg); *b = clip01(bb); *a = clip01(aa);
    return 1;
}

/* {cx,cy} from a 2-element numeric List. */
static int get_pt(const Expr* e, double* x, double* y) {
    if (!head_is(e, SYM_List) || e->data.function.arg_count != 2) return 0;
    return to_double(e->data.function.args[0], x)
        && to_double(e->data.function.args[1], y);
}

/* --------------------------------------------------- Bounding box pass ---- */

typedef struct { double xmin, xmax, ymin, ymax; } Range;

static void bb_pt(Range* r, double x, double y) {
    if (x < r->xmin) r->xmin = x;
    if (x > r->xmax) r->xmax = x;
    if (y < r->ymin) r->ymin = y;
    if (y > r->ymax) r->ymax = y;
}

/* Accumulate all point coordinates that contribute to the plot extent. */
static void bbox_walk(const Expr* p, Range* r) {
    if (!p) return;
    if (head_is(p, SYM_List)) {
        for (size_t i = 0; i < p->data.function.arg_count; i++)
            bbox_walk(p->data.function.args[i], r);
        return;
    }
    if ((head_is(p, SYM_Line) || head_is(p, SYM_Polygon) || head_is(p, SYM_Arrow))
        && p->data.function.arg_count >= 1) {
        const Expr* pts = p->data.function.args[0];
        if (head_is(pts, SYM_List))
            for (size_t i = 0; i < pts->data.function.arg_count; i++) {
                double x, y;
                if (get_pt(pts->data.function.args[i], &x, &y)) bb_pt(r, x, y);
            }
        return;
    }
    if (head_is(p, SYM_Point) && p->data.function.arg_count >= 1) {
        const Expr* a = p->data.function.args[0];
        double x, y;
        if (get_pt(a, &x, &y)) { bb_pt(r, x, y); return; }
        if (head_is(a, SYM_List))
            for (size_t i = 0; i < a->data.function.arg_count; i++)
                if (get_pt(a->data.function.args[i], &x, &y)) bb_pt(r, x, y);
        return;
    }
    if ((head_is(p, SYM_Disk) || head_is(p, SYM_Circle))
        && p->data.function.arg_count >= 1) {
        double cx = 0, cy = 0, rad = 1;
        get_pt(p->data.function.args[0], &cx, &cy);
        if (p->data.function.arg_count >= 2) to_double(p->data.function.args[1], &rad);
        bb_pt(r, cx - rad, cy - rad); bb_pt(r, cx + rad, cy + rad);
        return;
    }
    if (head_is(p, SYM_Rectangle) && p->data.function.arg_count >= 2) {
        double x1, y1, x2, y2;
        if (get_pt(p->data.function.args[0], &x1, &y1)
            && get_pt(p->data.function.args[1], &x2, &y2)) {
            bb_pt(r, x1, y1); bb_pt(r, x2, y2);
        }
        return;
    }
    if (head_is(p, SYM_Text) && p->data.function.arg_count >= 2) {
        double x, y;
        if (get_pt(p->data.function.args[1], &x, &y)) bb_pt(r, x, y);
        return;
    }
}

/* -------------------------------------------------------- Options parse --- */

typedef struct {
    Range range;          /* plot range in world coords */
    int   have_range;     /* PlotRange gave explicit numeric bounds */
    int   frame;          /* Frame -> True */
    int   axes;           /* Axes -> True (default for Plot) */
    int   axes_set;       /* Axes option was present */
    double bg_r, bg_g, bg_b; int have_bg;
    double width, height; /* page size in points (from ImageSize) */
} Opts;

/* Look for option `sym -> value` among the Graphics args (index >= 1). */
static const Expr* find_option(const Expr* g, const char* sym) {
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        const Expr* a = g->data.function.args[i];
        if (head_is(a, SYM_Rule) && a->data.function.arg_count == 2) {
            const Expr* lhs = a->data.function.args[0];
            if (lhs && lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == sym)
                return a->data.function.args[1];
        }
    }
    return NULL;
}

static int opt_is_true(const Expr* v) {
    return v && v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True;
}

/* Parse PlotRange -> {{x0,x1},{y0,y1}} when fully numeric. */
static void parse_range(const Expr* v, Opts* o) {
    if (!head_is(v, SYM_List) || v->data.function.arg_count != 2) return;
    const Expr* xr = v->data.function.args[0];
    const Expr* yr = v->data.function.args[1];
    double x0, x1, y0, y1;
    if (head_is(xr, SYM_List) && xr->data.function.arg_count == 2
        && head_is(yr, SYM_List) && yr->data.function.arg_count == 2
        && to_double(xr->data.function.args[0], &x0)
        && to_double(xr->data.function.args[1], &x1)
        && to_double(yr->data.function.args[0], &y0)
        && to_double(yr->data.function.args[1], &y1)) {
        o->range.xmin = x0; o->range.xmax = x1;
        o->range.ymin = y0; o->range.ymax = y1;
        o->have_range = 1;
    }
}

static void parse_opts(const Expr* g, Opts* o) {
    o->have_range = 0; o->frame = 0; o->axes = 1; o->axes_set = 0;
    o->have_bg = 0; o->width = 504.0; o->height = 360.0; /* 7in x 5in default */

    const Expr* v;
    if ((v = find_option(g, SYM_PlotRange))) parse_range(v, o);
    if ((v = find_option(g, SYM_Frame)))     o->frame = opt_is_true(v);
    if ((v = find_option(g, SYM_Axes)))      { o->axes = opt_is_true(v); o->axes_set = 1; }
    if ((v = find_option(g, SYM_Background))) {
        double r, gg, b, a;
        if (resolve_color(v, &r, &gg, &b, &a)) { o->bg_r=r; o->bg_g=gg; o->bg_b=b; o->have_bg=1; }
    }
    if ((v = find_option(g, SYM_ImageSize))) {
        double w, h;
        if (to_double(v, &w) && w > 0) { o->width = w; o->height = w * 0.7; }
        else if (head_is(v, SYM_List) && v->data.function.arg_count == 2
                 && to_double(v->data.function.args[0], &w)
                 && to_double(v->data.function.args[1], &h) && w > 0 && h > 0) {
            o->width = w; o->height = h;
        }
    }
}

/* ---------------------------------------------------------- tick policy --- */

/* "Nice" tick step (1/2/5 x 10^k) for a range and a target count. */
static double nice_step(double range, int target) {
    if (range <= 0 || target < 1) return 1.0;
    double raw = range / target;
    double mag = pow(10.0, floor(log10(raw)));
    double n = raw / mag;
    double nice = (n < 1.5) ? 1 : (n < 3) ? 2 : (n < 7) ? 5 : 10;
    return nice * mag;
}

/* Format a tick label without trailing zero noise. */
static void fmt_tick(char* out, size_t sz, double v, double step) {
    if (fabs(v) < step * 1e-6) v = 0.0;               /* kill -0 / fuzz */
    int dec = 0; double s = step;
    while (s < 1.0 - 1e-9 && dec < 8) { s *= 10.0; dec++; }
    snprintf(out, sz, "%.*f", dec, v);
}

/* ---------------------------------------------------- content emission ---- */

/* Opacity is applied through PDF ExtGState objects, one per distinct alpha<1,
 * collected here and referenced from the content as /GSk inside a q..Q pair. */
#define MAX_ALPHA 16
typedef struct {
    Buf*   c;                     /* content stream */
    double alphas[MAX_ALPHA];
    int    nalpha;
    int    used_text;             /* did we emit any text? (font resource) */
    /* current graphics state */
    double r, g, b, opacity, thickness, point_size;
    /* world->page transform */
    double ox, oy, sx, sy;
} Emit;

static double X(const Emit* e, double wx) { return e->ox + (wx) * e->sx; }
static double Y(const Emit* e, double wy) { return e->oy + (wy) * e->sy; }

static int gs_for_alpha(Emit* e, double a) {
    for (int i = 0; i < e->nalpha; i++)
        if (fabs(e->alphas[i] - a) < 1e-6) return i;
    if (e->nalpha < MAX_ALPHA) { e->alphas[e->nalpha] = a; return e->nalpha++; }
    return -1;
}

static void set_fill(Emit* e)   { buf_catf(e->c, "%.4f %.4f %.4f rg\n", e->r, e->g, e->b); }
static void set_stroke(Emit* e) { buf_catf(e->c, "%.4f %.4f %.4f RG\n", e->r, e->g, e->b); }

/* Stroke width in points: Thickness is a fraction of the plot width. */
static double stroke_w(const Emit* e, double plot_w) {
    double w = e->thickness > 0 ? e->thickness * plot_w : 1.0;
    if (w < 0.4) w = 0.4;
    return w;
}

/* Emit a circle centred at page (cx,cy), radius rp, as four Beziers. */
static void path_circle(Emit* e, double cx, double cy, double rp) {
    double k = KAPPA * rp;
    buf_catf(e->c, "%.3f %.3f m\n", cx + rp, cy);
    buf_catf(e->c, "%.3f %.3f %.3f %.3f %.3f %.3f c\n", cx+rp, cy+k, cx+k, cy+rp, cx, cy+rp);
    buf_catf(e->c, "%.3f %.3f %.3f %.3f %.3f %.3f c\n", cx-k, cy+rp, cx-rp, cy+k, cx-rp, cy);
    buf_catf(e->c, "%.3f %.3f %.3f %.3f %.3f %.3f c\n", cx-rp, cy-k, cx-k, cy-rp, cx, cy-rp);
    buf_catf(e->c, "%.3f %.3f %.3f %.3f %.3f %.3f c\n", cx+k, cy-rp, cx+rp, cy-k, cx+rp, cy);
}

/* Escape a string for a PDF literal ( ). */
static void emit_pdf_string(Buf* c, const char* s) {
    buf_cat(c, "(");
    for (const char* p = s; *p; p++) {
        if (*p == '(' || *p == ')' || *p == '\\') { char t[3] = {'\\', *p, 0}; buf_cat(c, t); }
        else if (*p == '\n') buf_cat(c, "\\n");
        else if (*p == '\r') buf_cat(c, "\\r");
        else { char t[2] = { *p, 0 }; buf_cat(c, t); }
    }
    buf_cat(c, ")");
}

/* Text content -> a display string (caller frees), or NULL if unrenderable. */
static char* text_string(const Expr* e) {
    if (!e) return NULL;
    if (e->type == EXPR_STRING) {
        size_t n = strlen(e->data.string) + 1;
        char* s = (char*)malloc(n); if (s) memcpy(s, e->data.string, n); return s;
    }
    char tmp[64];
    if (e->type == EXPR_INTEGER) { snprintf(tmp, sizeof(tmp), "%lld", (long long)e->data.integer); }
    else if (e->type == EXPR_REAL) { snprintf(tmp, sizeof(tmp), "%g", e->data.real); }
    else return NULL;
    size_t n = strlen(tmp) + 1; char* s = (char*)malloc(n); if (s) memcpy(s, tmp, n); return s;
}

static void draw_prim(Emit* e, const Expr* p, double plot_w);

static void draw_list(Emit* e, const Expr* lst, double plot_w) {
    if (!head_is(lst, SYM_List)) return;
    for (size_t i = 0; i < lst->data.function.arg_count; i++)
        draw_prim(e, lst->data.function.args[i], plot_w);
}

static void draw_prim(Emit* e, const Expr* p, double plot_w) {
    if (!p) return;

    if (head_is(p, SYM_List)) { draw_list(e, p, plot_w); return; }

    /* Directives -------------------------------------------------------- */
    { double r, g, b, a;
      if (resolve_color(p, &r, &g, &b, &a)) {
          e->r = r; e->g = g; e->b = b;
          if (a < 1.0) e->opacity = a;    /* RGBColor's 4th arg sets opacity */
          return;
      } }
    if (head_is(p, SYM_Opacity) && p->data.function.arg_count >= 1) {
        double o; if (to_double(p->data.function.args[0], &o)) e->opacity = clip01(o);
        return;
    }
    if (head_is(p, SYM_Thickness) && p->data.function.arg_count >= 1) {
        double t; if (to_double(p->data.function.args[0], &t)) e->thickness = t;
        return;
    }
    if (head_is(p, SYM_PointSize) && p->data.function.arg_count >= 1) {
        double s; if (to_double(p->data.function.args[0], &s)) e->point_size = s;
        return;
    }

    /* Line -------------------------------------------------------------- */
    if (head_is(p, SYM_Line) && p->data.function.arg_count >= 1) {
        const Expr* pts = p->data.function.args[0];
        if (!head_is(pts, SYM_List)) return;
        set_stroke(e);
        buf_catf(e->c, "%.3f w\n", stroke_w(e, plot_w));
        int first = 1;
        for (size_t i = 0; i < pts->data.function.arg_count; i++) {
            double x, y;
            if (!get_pt(pts->data.function.args[i], &x, &y)) continue;
            buf_catf(e->c, "%.3f %.3f %s\n", X(e, x), Y(e, y), first ? "m" : "l");
            first = 0;
        }
        if (!first) buf_cat(e->c, "S\n");
        return;
    }

    /* Polygon ----------------------------------------------------------- */
    if (head_is(p, SYM_Polygon) && p->data.function.arg_count >= 1) {
        const Expr* pts = p->data.function.args[0];
        if (!head_is(pts, SYM_List)) return;
        int gsi = e->opacity < 1.0 ? gs_for_alpha(e, e->opacity) : -1;
        if (gsi >= 0) buf_catf(e->c, "q /GS%d gs\n", gsi);
        set_fill(e);
        int first = 1;
        for (size_t i = 0; i < pts->data.function.arg_count; i++) {
            double x, y;
            if (!get_pt(pts->data.function.args[i], &x, &y)) continue;
            buf_catf(e->c, "%.3f %.3f %s\n", X(e, x), Y(e, y), first ? "m" : "l");
            first = 0;
        }
        if (!first) buf_cat(e->c, "h f\n");
        if (gsi >= 0) buf_cat(e->c, "Q\n");
        return;
    }

    /* Disk / Circle ----------------------------------------------------- */
    if ((head_is(p, SYM_Disk) || head_is(p, SYM_Circle))
        && p->data.function.arg_count >= 1) {
        double cx = 0, cy = 0, rad = 1;
        get_pt(p->data.function.args[0], &cx, &cy);
        if (p->data.function.arg_count >= 2) to_double(p->data.function.args[1], &rad);
        double rp = rad * fabs(e->sx); /* radius in page units (x scale) */
        int filled = head_is(p, SYM_Disk);
        int gsi = (filled && e->opacity < 1.0) ? gs_for_alpha(e, e->opacity) : -1;
        if (gsi >= 0) buf_catf(e->c, "q /GS%d gs\n", gsi);
        if (filled) set_fill(e); else { set_stroke(e); buf_catf(e->c, "%.3f w\n", stroke_w(e, plot_w)); }
        path_circle(e, X(e, cx), Y(e, cy), rp);
        buf_cat(e->c, filled ? "f\n" : "S\n");
        if (gsi >= 0) buf_cat(e->c, "Q\n");
        return;
    }

    /* Point(s) ---------------------------------------------------------- */
    if (head_is(p, SYM_Point) && p->data.function.arg_count >= 1) {
        const Expr* a = p->data.function.args[0];
        double ps = (e->point_size > 0 ? e->point_size : 0.008) * plot_w;
        if (ps < 1.0) ps = 1.0;
        set_fill(e);
        double x, y;
        if (get_pt(a, &x, &y)) { path_circle(e, X(e, x), Y(e, y), ps); buf_cat(e->c, "f\n"); return; }
        if (head_is(a, SYM_List))
            for (size_t i = 0; i < a->data.function.arg_count; i++)
                if (get_pt(a->data.function.args[i], &x, &y)) {
                    path_circle(e, X(e, x), Y(e, y), ps); buf_cat(e->c, "f\n");
                }
        return;
    }

    /* Rectangle --------------------------------------------------------- */
    if (head_is(p, SYM_Rectangle) && p->data.function.arg_count >= 2) {
        double x1, y1, x2, y2;
        if (!get_pt(p->data.function.args[0], &x1, &y1)
            || !get_pt(p->data.function.args[1], &x2, &y2)) return;
        int gsi = e->opacity < 1.0 ? gs_for_alpha(e, e->opacity) : -1;
        if (gsi >= 0) buf_catf(e->c, "q /GS%d gs\n", gsi);
        set_fill(e);
        double px = X(e, x1), py = Y(e, y1);
        buf_catf(e->c, "%.3f %.3f %.3f %.3f re f\n", px, py, X(e, x2) - px, Y(e, y2) - py);
        if (gsi >= 0) buf_cat(e->c, "Q\n");
        return;
    }

    /* Arrow ------------------------------------------------------------- */
    if (head_is(p, SYM_Arrow) && p->data.function.arg_count >= 1) {
        const Expr* pts = p->data.function.args[0];
        if (!head_is(pts, SYM_List) || pts->data.function.arg_count < 2) return;
        set_stroke(e);
        buf_catf(e->c, "%.3f w\n", stroke_w(e, plot_w));
        double lx = 0, ly = 0, px = 0, py = 0; int first = 1;
        for (size_t i = 0; i < pts->data.function.arg_count; i++) {
            double x, y;
            if (!get_pt(pts->data.function.args[i], &x, &y)) continue;
            px = lx; py = ly; lx = X(e, x); ly = Y(e, y);
            buf_catf(e->c, "%.3f %.3f %s\n", lx, ly, first ? "m" : "l");
            first = 0;
        }
        if (!first) buf_cat(e->c, "S\n");
        /* Solid arrowhead on the final segment. */
        double dx = lx - px, dy = ly - py, len = sqrt(dx*dx + dy*dy);
        if (len > 1e-6) {
            double ux = dx/len, uy = dy/len, hl = 8.0, hw = 3.0;
            double bx = lx - ux*hl, by = ly - uy*hl;
            set_fill(e);
            buf_catf(e->c, "%.3f %.3f m %.3f %.3f l %.3f %.3f l h f\n",
                     lx, ly, bx - uy*hw, by + ux*hw, bx + uy*hw, by - ux*hw);
        }
        return;
    }

    /* Text -------------------------------------------------------------- */
    if (head_is(p, SYM_Text) && p->data.function.arg_count >= 2) {
        double x, y;
        if (!get_pt(p->data.function.args[1], &x, &y)) return;
        char* s = text_string(p->data.function.args[0]);
        if (!s) return;
        double fs = 10.0;
        double tx = X(e, x) - 0.25 * fs * (double)strlen(s); /* rough centring */
        double ty = Y(e, y) - 0.35 * fs;
        set_fill(e);
        buf_cat(e->c, "BT /F1 ");
        buf_catf(e->c, "%.1f Tf %.3f %.3f Td ", fs, tx, ty);
        emit_pdf_string(e->c, s);
        buf_cat(e->c, " Tj ET\n");
        e->used_text = 1;
        free(s);
        return;
    }
}

/* Axis frame + ticks + numeric labels, drawn in page space. */
static void draw_axes_frame(Emit* e, const Range* r, double rx, double ry,
                            double rw, double rh) {
    /* Frame box. */
    buf_cat(e->c, "0.2 0.2 0.2 RG\n0.8 w\n");
    buf_catf(e->c, "%.3f %.3f %.3f %.3f re S\n", rx, ry, rw, rh);

    double xstep = nice_step(r->xmax - r->xmin, 6);
    double ystep = nice_step(r->ymax - r->ymin, 6);
    double tick = 5.0; /* tick length in points */

    /* X ticks + labels along the bottom. */
    double xv = ceil(r->xmin / xstep) * xstep;
    for (; xv <= r->xmax + xstep * 1e-6; xv += xstep) {
        double px = X(e, xv);
        buf_catf(e->c, "%.3f %.3f m %.3f %.3f l S\n", px, ry, px, ry + tick);
        buf_catf(e->c, "%.3f %.3f m %.3f %.3f l S\n", px, ry + rh, px, ry + rh - tick);
        char lbl[32]; fmt_tick(lbl, sizeof(lbl), xv, xstep);
        buf_cat(e->c, "0 0 0 rg\nBT /F1 9 Tf ");
        buf_catf(e->c, "%.3f %.3f Td ", px - 3.0 * (double)strlen(lbl), ry - 12.0);
        emit_pdf_string(e->c, lbl);
        buf_cat(e->c, " Tj ET\n0.2 0.2 0.2 RG\n");
        e->used_text = 1;
    }
    /* Y ticks + labels along the left. */
    double yv = ceil(r->ymin / ystep) * ystep;
    for (; yv <= r->ymax + ystep * 1e-6; yv += ystep) {
        double py = Y(e, yv);
        buf_catf(e->c, "%.3f %.3f m %.3f %.3f l S\n", rx, py, rx + tick, py);
        buf_catf(e->c, "%.3f %.3f m %.3f %.3f l S\n", rx + rw, py, rx + rw - tick, py);
        char lbl[32]; fmt_tick(lbl, sizeof(lbl), yv, ystep);
        buf_cat(e->c, "0 0 0 rg\nBT /F1 9 Tf ");
        buf_catf(e->c, "%.3f %.3f Td ", rx - 6.0 * (double)strlen(lbl) - 4.0, py - 3.0);
        emit_pdf_string(e->c, lbl);
        buf_cat(e->c, " Tj ET\n0.2 0.2 0.2 RG\n");
        e->used_text = 1;
    }
}

/* ------------------------------------------------------- PDF assembly ----- */

int graphics_export_pdf(const Expr* g, const char* path) {
    if (!g || g->type != EXPR_FUNCTION || g->data.function.arg_count < 1) return 0;
    if (!head_is(g, SYM_Graphics)) return 0;   /* 2D only; Graphics3D unsupported */

    const Expr* prims = g->data.function.args[0];

    Opts o; parse_opts(g, &o);

    /* Resolve the plot range. */
    Range r;
    if (o.have_range) {
        r = o.range;
    } else {
        r.xmin = r.ymin = 1e300; r.xmax = r.ymax = -1e300;
        bbox_walk(prims, &r);
        if (r.xmin > r.xmax) { r.xmin = -1; r.xmax = 1; }
        if (r.ymin > r.ymax) { r.ymin = -1; r.ymax = 1; }
        double xpad = (r.xmax - r.xmin) * 0.02, ypad = (r.ymax - r.ymin) * 0.04;
        if (xpad <= 0) xpad = 1;
        if (ypad <= 0) ypad = 1;
        r.xmin -= xpad; r.xmax += xpad; r.ymin -= ypad; r.ymax += ypad;
    }
    double dw = r.xmax - r.xmin, dh = r.ymax - r.ymin;
    if (dw <= 0) dw = 1;
    if (dh <= 0) dh = 1;

    int draw_axes = o.frame || o.axes;   /* Plot defaults to Axes->True */

    /* Page + plot region. */
    double W = o.width, H = o.height;
    double mL = draw_axes ? 44.0 : 8.0, mR = 12.0;
    double mT = 12.0, mB = draw_axes ? 30.0 : 8.0;
    double rx = mL, ry = mB, rw = W - mL - mR, rh = H - mT - mB;
    if (rw < 20) rw = 20;
    if (rh < 20) rh = 20;

    /* Build the content stream. */
    Buf content; buf_init(&content);
    Emit e;
    memset(&e, 0, sizeof(e));
    e.c = &content; e.opacity = 1.0; e.thickness = 0; e.point_size = 0;
    e.r = 0.368; e.g = 0.507; e.b = 0.71; /* Mathematica default line blue */
    e.sx = rw / dw; e.sy = rh / dh;
    e.ox = rx - r.xmin * e.sx; e.oy = ry - r.ymin * e.sy;

    /* Background fill (page). */
    if (o.have_bg)
        buf_catf(&content, "%.4f %.4f %.4f rg 0 0 %.3f %.3f re f\n", o.bg_r, o.bg_g, o.bg_b, W, H);

    /* Clip the primitives to the plot region so nothing spills into margins. */
    buf_catf(&content, "q %.3f %.3f %.3f %.3f re W n\n", rx, ry, rw, rh);
    draw_prim(&e, prims, rw);
    buf_cat(&content, "Q\n");

    if (draw_axes) draw_axes_frame(&e, &r, rx, ry, rw, rh);

    /* ---- Assemble the PDF objects. --------------------------------------
     * 1 Catalog, 2 Pages, 3 Page, 4 Font, 5..(4+nalpha) ExtGState,
     * last = Contents. */
    int nobj = 4 + e.nalpha + 1;
    int content_obj = nobj;                 /* last object number */
    long* off = (long*)calloc((size_t)nobj + 1, sizeof(long));
    if (!off) { buf_free(&content); return 0; }

    Buf pdf; buf_init(&pdf);
    buf_cat(&pdf, "%PDF-1.4\n%\xE2\xE3\xCF\xD3\n");

    #define OBJ(n) do { off[n] = (long)pdf.len; buf_catf(&pdf, "%d 0 obj\n", (n)); } while (0)

    OBJ(1); buf_cat(&pdf, "<< /Type /Catalog /Pages 2 0 R >>\nendobj\n");
    OBJ(2); buf_cat(&pdf, "<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n");

    OBJ(3);
    buf_catf(&pdf, "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.2f %.2f]\n", W, H);
    buf_cat(&pdf, "   /Resources << /Font << /F1 4 0 R >>");
    if (e.nalpha > 0) {
        buf_cat(&pdf, " /ExtGState <<");
        for (int i = 0; i < e.nalpha; i++) buf_catf(&pdf, " /GS%d %d 0 R", i, 5 + i);
        buf_cat(&pdf, " >>");
    }
    buf_catf(&pdf, " >>\n   /Contents %d 0 R >>\nendobj\n", content_obj);

    OBJ(4);
    buf_cat(&pdf, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>\nendobj\n");

    for (int i = 0; i < e.nalpha; i++) {
        OBJ(5 + i);
        buf_catf(&pdf, "<< /Type /ExtGState /ca %.4f /CA %.4f >>\nendobj\n",
                 e.alphas[i], e.alphas[i]);
    }

    OBJ(content_obj);
    buf_catf(&pdf, "<< /Length %zu >>\nstream\n", content.len);
    buf_cat(&pdf, content.p ? content.p : "");
    buf_cat(&pdf, "endstream\nendobj\n");

    long xref = (long)pdf.len;
    buf_catf(&pdf, "xref\n0 %d\n", nobj + 1);
    buf_cat(&pdf, "0000000000 65535 f \n");
    for (int n = 1; n <= nobj; n++) buf_catf(&pdf, "%010ld 00000 n \n", off[n]);
    buf_catf(&pdf, "trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%ld\n%%%%EOF\n",
             nobj + 1, xref);

    #undef OBJ

    /* Write it out. */
    FILE* f = fopen(path, "wb");
    int ok = 0;
    if (f) {
        ok = (fwrite(pdf.p, 1, pdf.len, f) == pdf.len);
        fclose(f);
    }
    free(off);
    buf_free(&content);
    buf_free(&pdf);
    return ok;
}
