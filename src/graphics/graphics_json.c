/*
 * graphics_json.c
 *
 * Convert a Mathilda Graphics[...] expression to a Plotly JSON string
 * suitable for sending through the NDJSON pipe protocol to the notebook
 * frontend.
 *
 * The Graphics Expr produced by Plot[] looks like:
 *
 *   Graphics[
 *     List[
 *       Opacity[0.3],
 *       Polygon[List[List[x,y], ...]],
 *       RGBColor[r, g, b],
 *       Line[List[List[x,y], ...]],
 *       ...  (more curves repeat the pattern)
 *     ],
 *     Rule[Axes, True],
 *     Rule[AspectRatio, 0.618034],
 *     Rule[PlotStyle, RGBColor[r,g,b]],
 *     $PlotResample[...]
 *   ]
 *
 * GraphPlot[] instead emits node-link diagrams built from Line (edges),
 * Disk (vertices), and Text (labels). When Disk or Text primitives are
 * present the converter switches to "diagram mode": vertices become marker
 * points, labels become annotations, and the layout uses an equal-aspect
 * (scaleanchor) square with hidden axes -- so a graph draws as a graph, not
 * as an XY function plot.
 *
 * Output JSON (Plotly format):
 *
 *   {
 *     "data": [
 *       {"type":"scatter","x":[...],"y":[...],"mode":"lines",
 *        "line":{"color":"rgba(51,102,204,1.0)","width":1.5}},
 *       ...
 *     ],
 *     "layout": { ... }
 *   }
 */

#include "graphics_json.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* -----------------------------------------------------------------------
 * Dynamic string buffer
 * --------------------------------------------------------------------- */

typedef struct {
    char*  buf;
    size_t len;
    size_t cap;
} Buf;

static int buf_init(Buf* b, size_t initial) {
    b->buf = malloc(initial);
    if (!b->buf) return 0;
    b->buf[0] = '\0';
    b->len = 0;
    b->cap = initial;
    return 1;
}

static int buf_grow(Buf* b, size_t need) {
    if (b->len + need + 1 <= b->cap) return 1;
    size_t new_cap = (b->cap + need) * 2;
    char* p = realloc(b->buf, new_cap);
    if (!p) return 0;
    b->buf = p;
    b->cap = new_cap;
    return 1;
}

static int buf_cat(Buf* b, const char* s) {
    size_t n = strlen(s);
    if (!buf_grow(b, n)) return 0;
    memcpy(b->buf + b->len, s, n + 1);
    b->len += n;
    return 1;
}

/* Append a double value — use %g to keep it compact. Snap sub-noise
 * magnitudes to 0 so floating-point dust (e.g. sin(pi) ~ 1e-16) does not
 * blow up an auto-ranged axis. */
static int buf_catd(Buf* b, double v) {
    char tmp[64];
    if (fabs(v) < 1e-12) v = 0.0;
    snprintf(tmp, sizeof(tmp), "%.7g", v);
    return buf_cat(b, tmp);
}

static void buf_free(Buf* b) { free(b->buf); b->buf = NULL; b->len = b->cap = 0; }

/* -----------------------------------------------------------------------
 * Expr helpers
 * --------------------------------------------------------------------- */

static int head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

static int expr_to_double(const Expr* e, double* out) {
    if (!e) return 0;
    if (e->type == EXPR_REAL)    { *out = e->data.real;               return 1; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer;    return 1; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint);  return 1; }
    return 0;
}


/* -----------------------------------------------------------------------
 * Coordinate-list serialiser
 * Appends "[x1,x2,...]" or "[y1,y2,...]" to buf from
 * List[List[x,y], ...] using component index (0=x, 1=y).
 * --------------------------------------------------------------------- */

static int append_coord_array(Buf* b, const Expr* pts, int component) {
    if (!head_is(pts, SYM_List)) return 0;
    buf_cat(b, "[");
    int first = 1;
    for (size_t i = 0; i < pts->data.function.arg_count; i++) {
        const Expr* pair = pts->data.function.args[i];
        if (!head_is(pair, SYM_List) || pair->data.function.arg_count < 2) continue;
        double v;
        if (!expr_to_double(pair->data.function.args[component], &v)) continue;
        if (!first) buf_cat(b, ",");
        buf_catd(b, v);
        first = 0;
    }
    return buf_cat(b, "]");
}

/* -----------------------------------------------------------------------
 * rgba() string helper
 * --------------------------------------------------------------------- */

static void rgba_str(char* out, size_t outsz, double r, double g, double b, double a) {
    int ri = (int)(r * 255.0 + 0.5);
    int gi = (int)(g * 255.0 + 0.5);
    int bi = (int)(b * 255.0 + 0.5);
    if (ri < 0) ri = 0;
    if (ri > 255) ri = 255;
    if (gi < 0) gi = 0;
    if (gi > 255) gi = 255;
    if (bi < 0) bi = 0;
    if (bi > 255) bi = 255;
    snprintf(out, outsz, "rgba(%d,%d,%d,%.3f)", ri, gi, bi, a);
}

/* -----------------------------------------------------------------------
 * Color resolution — mirrors render.c's resolve_color() so the JSON
 * output is consistent with the Raylib renderer for all four color forms.
 * Returns 1 and writes (r,g,b,a) in [0,1] on success; 0 otherwise.
 * --------------------------------------------------------------------- */

static double clip01(double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); }

static int resolve_color_rgb(const Expr* e,
                              double* r_out, double* g_out,
                              double* b_out, double* a_out) {
    if (!e || e->type != EXPR_FUNCTION
        || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL)
        return 0;

    const char*   h    = e->data.function.head->data.symbol.name;
    const Expr**  args = (const Expr**)e->data.function.args;
    size_t        n    = e->data.function.arg_count;
    double r = 0, g = 0, b = 0, a = 1;

    if (h == SYM_RGBColor) {
        if (n < 3) return 0;
        if (!expr_to_double(args[0], &r) ||
            !expr_to_double(args[1], &g) ||
            !expr_to_double(args[2], &b)) return 0;
        if (n >= 4) expr_to_double(args[3], &a);

    } else if (h == SYM_GrayLevel) {
        if (n < 1) return 0;
        double lv = 0;
        if (!expr_to_double(args[0], &lv)) return 0;
        r = g = b = lv;
        if (n >= 2) expr_to_double(args[1], &a);

    } else if (h == SYM_Hue) {
        /* Hue[h], Hue[h,s,v], Hue[h,s,v,a] — HSB/HSV to RGB. */
        double hh = 0, s = 1, v = 1;
        if (n >= 1) expr_to_double(args[0], &hh);
        if (n >= 3) { expr_to_double(args[1], &s); expr_to_double(args[2], &v); }
        if (n >= 4) expr_to_double(args[3], &a);
        hh = hh - floor(hh);
        s = clip01(s); v = clip01(v);
        if (s <= 0.0) {
            r = g = b = v;
        } else {
            double hh6 = hh * 6.0;
            int   qi   = (int)floor(hh6);
            double f   = hh6 - qi;
            double p   = v * (1.0 - s);
            double q_  = v * (1.0 - f * s);
            double t_  = v * (1.0 - (1.0 - f) * s);
            switch (((qi % 6) + 6) % 6) {
                case 0: r = v;  g = t_; b = p;  break;
                case 1: r = q_; g = v;  b = p;  break;
                case 2: r = p;  g = v;  b = t_; break;
                case 3: r = p;  g = q_; b = v;  break;
                case 4: r = t_; g = p;  b = v;  break;
                default: r = v; g = p;  b = q_; break;
            }
        }

    } else if (h == SYM_CMYKColor) {
        /* CMYKColor[c,m,y], [c,m,y,k], [c,m,y,k,a], or list form. */
        const Expr** ca = args;
        size_t cn = n;
        if (cn == 1 && args[0] && args[0]->type == EXPR_FUNCTION
            && args[0]->data.function.head->type == EXPR_SYMBOL
            && args[0]->data.function.head->data.symbol.name == SYM_List) {
            cn = args[0]->data.function.arg_count;
            ca = (const Expr**)args[0]->data.function.args;
        }
        double c = 0, m = 0, y = 0, k = 0;
        if (cn >= 3) { expr_to_double(ca[0], &c); expr_to_double(ca[1], &m); expr_to_double(ca[2], &y); }
        if (cn >= 4) expr_to_double(ca[3], &k);
        if (cn >= 5) expr_to_double(ca[4], &a);
        double w = 1.0 - clip01(k);
        r = (1.0 - clip01(c)) * w;
        g = (1.0 - clip01(m)) * w;
        b = (1.0 - clip01(y)) * w;
        a = clip01(a);

    } else {
        return 0;
    }

    *r_out = clip01(r); *g_out = clip01(g);
    *b_out = clip01(b); *a_out = clip01(a);
    return 1;
}

/* -----------------------------------------------------------------------
 * Geometry helpers shared by 2D primitives
 * --------------------------------------------------------------------- */

/* Extract {cx, cy} from a 2-element List; defaults to {0,0}. */
static int get_center(const Expr* e, double* cx, double* cy) {
    *cx = 0; *cy = 0;
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    return expr_to_double(e->data.function.args[0], cx)
        && expr_to_double(e->data.function.args[1], cy);
}

/* Parse Disk[]/Disk[{cx,cy}]/Disk[{cx,cy},r] — same defaults as render.c. */
static int disk_params(const Expr* e, double* cx, double* cy, double* r) {
    size_t n = e->data.function.arg_count;
    *cx = 0; *cy = 0; *r = 1;
    if (n >= 1 && !get_center(e->data.function.args[0], cx, cy)) return 0;
    if (n >= 2 && !expr_to_double(e->data.function.args[1], r)) return 0;
    return 1;
}

/* Append a Plotly scatter trace representing a circle/disk approximated by
 * NSEG segments. is_filled=1 → Disk (filled), 0 → Circle (outline). */
#define DISK_NSEG 64
static void append_disk_trace(Buf* b, double cx, double cy, double r,
                               double col_r, double col_g, double col_b,
                               double opacity, int is_filled, int* first_trace) {
    char color_str[64];
    rgba_str(color_str, sizeof(color_str), col_r, col_g, col_b, 1.0);

    if (!(*first_trace)) buf_cat(b, ",");
    *first_trace = 0;

    buf_cat(b, "{\"type\":\"scatter\",\"mode\":\"lines\",");
    if (is_filled) {
        char fill_str[64];
        rgba_str(fill_str, sizeof(fill_str), col_r, col_g, col_b, opacity);
        buf_cat(b, "\"fill\":\"toself\",\"fillcolor\":\"");
        buf_cat(b, fill_str);
        buf_cat(b, "\",");
    }
    buf_cat(b, "\"x\":[");
    for (int i = 0; i <= DISK_NSEG; i++) {
        double angle = 2.0 * M_PI * i / DISK_NSEG;
        if (i) buf_cat(b, ",");
        buf_catd(b, cx + r * cos(angle));
    }
    buf_cat(b, "],\"y\":[");
    for (int i = 0; i <= DISK_NSEG; i++) {
        double angle = 2.0 * M_PI * i / DISK_NSEG;
        if (i) buf_cat(b, ",");
        buf_catd(b, cy + r * sin(angle));
    }
    buf_cat(b, "],\"line\":{\"color\":\"");
    buf_cat(b, color_str);
    buf_cat(b, "\",\"width\":");
    buf_cat(b, is_filled ? "0" : "1.5");
    buf_cat(b, "},\"showlegend\":false}");
}

/* -----------------------------------------------------------------------
 * 2D drawing state (shared across the recursive primitive walker)
 * --------------------------------------------------------------------- */

typedef struct {
    Buf*   buf;
    double cur_r, cur_g, cur_b;
    double cur_opacity;
    int    first_trace;
} State2D;

/* Forward declaration for mutual recursion (List recurses into draw_prim_2d). */
static void draw_prim_2d(State2D* s, const Expr* p);

static void draw_prim_list_2d(State2D* s, const Expr* lst) {
    if (!head_is(lst, SYM_List)) return;
    for (size_t i = 0; i < lst->data.function.arg_count; i++)
        draw_prim_2d(s, lst->data.function.args[i]);
}

static void draw_prim_2d(State2D* s, const Expr* p) {
    if (!p) return;

    /* ----- Nested List — recurse (handles Table[{color,prim},...]) --- */
    if (head_is(p, SYM_List)) {
        draw_prim_list_2d(s, p);
        return;
    }

    /* ----- Color directive (RGBColor / GrayLevel / Hue / CMYKColor) --- */
    {
        double r, g_c, b, a;
        if (resolve_color_rgb(p, &r, &g_c, &b, &a)) {
            s->cur_r = r; s->cur_g = g_c; s->cur_b = b;
            return;
        }
    }

    /* ----- Opacity[n] — update fill opacity ----------------------- */
    if (head_is(p, SYM_Opacity) && p->data.function.arg_count >= 1) {
        double op;
        if (expr_to_double(p->data.function.args[0], &op))
            s->cur_opacity = op;
        return;
    }

    /* ----- Line[List[pts...]] — line trace ----------------------- */
    if (head_is(p, SYM_Line) && p->data.function.arg_count >= 1) {
        const Expr* pts = p->data.function.args[0];
        if (!head_is(pts, SYM_List)) return;

        char color_str[64];
        rgba_str(color_str, sizeof(color_str), s->cur_r, s->cur_g, s->cur_b, 1.0);

        if (!s->first_trace) buf_cat(s->buf, ",");
        s->first_trace = 0;

        buf_cat(s->buf, "{\"type\":\"scatter\",\"mode\":\"lines\",");
        buf_cat(s->buf, "\"x\":");
        append_coord_array(s->buf, pts, 0);
        buf_cat(s->buf, ",\"y\":");
        append_coord_array(s->buf, pts, 1);
        buf_cat(s->buf, ",\"line\":{\"color\":\"");
        buf_cat(s->buf, color_str);
        buf_cat(s->buf, "\",\"width\":1.5},\"showlegend\":false}");
        return;
    }

    /* ----- Polygon[List[pts...]] — filled area ------------------- */
    if (head_is(p, SYM_Polygon) && p->data.function.arg_count >= 1) {
        const Expr* pts = p->data.function.args[0];
        if (!head_is(pts, SYM_List)) return;

        char fill_str[64];
        rgba_str(fill_str, sizeof(fill_str), s->cur_r, s->cur_g, s->cur_b, s->cur_opacity);

        if (!s->first_trace) buf_cat(s->buf, ",");
        s->first_trace = 0;

        buf_cat(s->buf, "{\"type\":\"scatter\",\"mode\":\"lines\",");
        buf_cat(s->buf, "\"fill\":\"toself\",\"fillcolor\":\"");
        buf_cat(s->buf, fill_str);
        buf_cat(s->buf, "\",\"x\":");
        append_coord_array(s->buf, pts, 0);
        buf_cat(s->buf, ",\"y\":");
        append_coord_array(s->buf, pts, 1);
        buf_cat(s->buf, ",\"line\":{\"color\":\"transparent\"},\"showlegend\":false}");
        return;
    }

    /* ----- Disk[{cx,cy}, r] / Circle[{cx,cy}, r] ---------------- */
    if ((head_is(p, SYM_Disk) || head_is(p, SYM_Circle))
            && p->data.function.arg_count <= 2) {
        double cx, cy, r;
        if (disk_params(p, &cx, &cy, &r)) {
            int is_filled = head_is(p, SYM_Disk);
            append_disk_trace(s->buf, cx, cy, r,
                              s->cur_r, s->cur_g, s->cur_b,
                              s->cur_opacity, is_filled, &s->first_trace);
        }
        return;
    }

    /* ----- Point[{x,y}] / Point[{{x1,y1},...}] — markers -------- */
    if (head_is(p, SYM_Point) && p->data.function.arg_count >= 1) {
        const Expr* arg = p->data.function.args[0];
        if (!head_is(arg, SYM_List)) return;

        char color_str[64];
        rgba_str(color_str, sizeof(color_str), s->cur_r, s->cur_g, s->cur_b, 1.0);

        /* Discriminate: {x,y} with numeric elements is a single point;
         * otherwise treat as a list of points. */
        double px, py;
        if (get_center(arg, &px, &py)) {
            if (!s->first_trace) buf_cat(s->buf, ",");
            s->first_trace = 0;
            buf_cat(s->buf, "{\"type\":\"scatter\",\"mode\":\"markers\",");
            char pt_buf[256];
            snprintf(pt_buf, sizeof(pt_buf),
                     "\"x\":[%.10g],\"y\":[%.10g],"
                     "\"marker\":{\"color\":\"%s\",\"size\":6},"
                     "\"showlegend\":false}",
                     px, py, color_str);
            buf_cat(s->buf, pt_buf);
        } else {
            /* list of points */
            if (!s->first_trace) buf_cat(s->buf, ",");
            s->first_trace = 0;
            buf_cat(s->buf, "{\"type\":\"scatter\",\"mode\":\"markers\",\"x\":[");
            int fp = 1;
            for (size_t k = 0; k < arg->data.function.arg_count; k++) {
                double qx, qy;
                if (!get_center(arg->data.function.args[k], &qx, &qy)) continue;
                if (!fp) buf_cat(s->buf, ",");
                fp = 0; buf_catd(s->buf, qx);
            }
            buf_cat(s->buf, "],\"y\":["); fp = 1;
            for (size_t k = 0; k < arg->data.function.arg_count; k++) {
                double qx, qy;
                if (!get_center(arg->data.function.args[k], &qx, &qy)) continue;
                if (!fp) buf_cat(s->buf, ",");
                fp = 0; buf_catd(s->buf, qy);
            }
            char ms[128];
            snprintf(ms, sizeof(ms),
                     "],\"marker\":{\"color\":\"%s\",\"size\":6},"
                     "\"showlegend\":false}", color_str);
            buf_cat(s->buf, ms);
        }
        return;
    }

    /* ----- Rectangle[{x1,y1}, {x2,y2}] — filled box ------------- */
    if (head_is(p, SYM_Rectangle) && p->data.function.arg_count >= 2) {
        double x1, y1, x2, y2;
        if (get_center(p->data.function.args[0], &x1, &y1)
                && get_center(p->data.function.args[1], &x2, &y2)) {
            char fill_str[64];
            rgba_str(fill_str, sizeof(fill_str), s->cur_r, s->cur_g, s->cur_b, s->cur_opacity);
            char col_str[64];
            rgba_str(col_str, sizeof(col_str), s->cur_r, s->cur_g, s->cur_b, 1.0);
            if (!s->first_trace) buf_cat(s->buf, ",");
            s->first_trace = 0;
            char rect_buf[512];
            snprintf(rect_buf, sizeof(rect_buf),
                     "{\"type\":\"scatter\",\"mode\":\"lines\","
                     "\"fill\":\"toself\",\"fillcolor\":\"%s\","
                     "\"x\":[%.10g,%.10g,%.10g,%.10g,%.10g],"
                     "\"y\":[%.10g,%.10g,%.10g,%.10g,%.10g],"
                     "\"line\":{\"color\":\"%s\",\"width\":0},"
                     "\"showlegend\":false}",
                     fill_str,
                     x1, x2, x2, x1, x1,
                     y1, y1, y2, y2, y1,
                     col_str);
            buf_cat(s->buf, rect_buf);
        }
        return;
    }

    /* All other primitives (Rule, $PlotResample, etc.) — skip. */
}

/* -----------------------------------------------------------------------
 * Main converter
 * --------------------------------------------------------------------- */

char* graphics_to_plotly_json(const Expr* g) {
    if (!g || !head_is(g, SYM_Graphics) || g->data.function.arg_count < 1)
        return NULL;

    /* The first argument is the primitive list. */
    const Expr* prim_list = g->data.function.args[0];
    if (!head_is(prim_list, SYM_List)) return NULL;

    size_t n = prim_list->data.function.arg_count;

    /* Diagram mode: a node-link picture (GraphPlot) rather than a function
     * plot. Signalled by the presence of Disk (vertices) or Text (labels),
     * which Plot[] never emits. */
    int diagram_mode = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* p = prim_list->data.function.args[i];
        if (head_is(p, SYM_Disk) || head_is(p, SYM_Text)) { diagram_mode = 1; break; }
    }

    Buf data_buf, anno_buf;
    if (!buf_init(&data_buf, 65536)) return NULL;
    if (!buf_init(&anno_buf, 4096)) { buf_free(&data_buf); return NULL; }
    buf_cat(&anno_buf, "[");

    State2D state = {
        .buf         = &data_buf,
        .cur_r       = 0.2, .cur_g = 0.4, .cur_b = 0.8, /* Mathematica default blue */
        .cur_opacity = 0.3,
        .first_trace = 1,
    };

    buf_cat(&data_buf, "[");
    draw_prim_list_2d(&state, prim_list);
    buf_cat(&data_buf, "]");
    buf_cat(&anno_buf, "]");

    /* Build the layout object. */
    Buf layout;
    if (!buf_init(&layout, 1024)) { buf_free(&data_buf); buf_free(&anno_buf); return NULL; }
    if (diagram_mode) {
        /* Square, axis-free canvas with equal x/y scaling so a circular
         * vertex layout stays circular. */
        buf_cat(&layout,
            "{\"xaxis\":{\"visible\":false,\"showgrid\":false,\"zeroline\":false},"
            "\"yaxis\":{\"visible\":false,\"showgrid\":false,\"zeroline\":false,"
            "\"scaleanchor\":\"x\",\"scaleratio\":1},"
            "\"margin\":{\"l\":20,\"r\":20,\"t\":20,\"b\":20},"
            "\"height\":400,\"showlegend\":false,"
            "\"plot_bgcolor\":\"#fff\",\"paper_bgcolor\":\"#fff\",");
        buf_cat(&layout, "\"annotations\":");
        buf_cat(&layout, anno_buf.buf);
        buf_cat(&layout, "}");
    } else {
        buf_cat(&layout,
            "{\"xaxis\":{\"showgrid\":true,\"zeroline\":true,"
            "\"zerolinecolor\":\"#aaa\",\"zerolinewidth\":1},"
            "\"yaxis\":{\"showgrid\":true,\"zeroline\":true,"
            "\"zerolinecolor\":\"#aaa\",\"zerolinewidth\":1},"
            "\"margin\":{\"l\":50,\"r\":20,\"t\":20,\"b\":50},"
            "\"height\":350,"
            "\"plot_bgcolor\":\"#fff\",\"paper_bgcolor\":\"#fff\"}");
    }

    /* Assemble the final JSON object. */
    Buf out;
    if (!buf_init(&out, data_buf.len + layout.len + 64)) {
        buf_free(&data_buf); buf_free(&anno_buf); buf_free(&layout);
        return NULL;
    }
    buf_cat(&out, "{\"data\":");
    buf_cat(&out, data_buf.buf);
    buf_cat(&out, ",\"layout\":");
    buf_cat(&out, layout.buf);
    buf_cat(&out, "}");

    buf_free(&data_buf);
    buf_free(&anno_buf);
    buf_free(&layout);
    return out.buf; /* caller frees */
}

/* -----------------------------------------------------------------------
 * 3D coordinate helpers
 * --------------------------------------------------------------------- */

/* Append one numeric component of a 3-element point list. */
static int append_coord3_array(Buf* b, const Expr* pts, int component) {
    if (!head_is(pts, SYM_List)) return 0;
    buf_cat(b, "[");
    int first = 1;
    for (size_t i = 0; i < pts->data.function.arg_count; i++) {
        const Expr* triple = pts->data.function.args[i];
        if (!head_is(triple, SYM_List) || triple->data.function.arg_count < 3) continue;
        double v;
        if (!expr_to_double(triple->data.function.args[component], &v)) continue;
        if (!first) buf_cat(b, ",");
        buf_catd(b, v);
        first = 0;
    }
    return buf_cat(b, "]");
}

/* -----------------------------------------------------------------------
 * graphics3d_to_plotly_json
 *
 * Converts Graphics3D[...] produced by ParametricPlot3D / Plot3D into
 * Plotly scatter3d (space curves from Line[]) and mesh3d (surfaces from
 * Polygon[] quads) traces.
 *
 * scatter3d carries the curve color via line.color; mesh3d carries it via
 * the color field.  Both respect RGBColor[]/Opacity[] directives in the
 * primitive list exactly as the 2D converter does.
 * --------------------------------------------------------------------- */

char* graphics3d_to_plotly_json(const Expr* g) {
    if (!g || !head_is(g, SYM_Graphics3D) || g->data.function.arg_count < 1)
        return NULL;

    const Expr* prim_list = g->data.function.args[0];
    if (!head_is(prim_list, SYM_List)) return NULL;

    Buf data_buf;
    if (!buf_init(&data_buf, 65536)) return NULL;

    double cur_r = 0.4, cur_g_val = 0.7, cur_b = 1.0; /* default: light blue */
    double cur_opacity = 0.85;
    int first_trace = 1;

    /* mesh3d accumulates all polygon vertices into one trace per color block.
     * We flush at every color change and at the end. */
    Buf mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k;
    buf_init(&mesh_x, 4096); buf_init(&mesh_y, 4096); buf_init(&mesh_z, 4096);
    buf_init(&mesh_i, 2048); buf_init(&mesh_j, 2048); buf_init(&mesh_k, 2048);
    size_t mesh_vcount = 0;  /* total vertices so far in current mesh block */
    size_t mesh_qcount = 0;  /* total quads so far in current mesh block */
    double mesh_r = cur_r, mesh_g_val = cur_g_val, mesh_b_color = cur_b;
    double mesh_opacity = cur_opacity;

    /* Emit the accumulated mesh3d trace and reset buffers. */
#define FLUSH_MESH() do { \
    if (mesh_qcount > 0) { \
        char mc[64]; \
        rgba_str(mc, sizeof(mc), mesh_r, mesh_g_val, mesh_b_color, mesh_opacity); \
        if (!first_trace) buf_cat(&data_buf, ","); \
        first_trace = 0; \
        buf_cat(&data_buf, "{\"type\":\"mesh3d\","); \
        buf_cat(&data_buf, "\"x\":"); buf_cat(&data_buf, mesh_x.buf); \
        buf_cat(&data_buf, ",\"y\":"); buf_cat(&data_buf, mesh_y.buf); \
        buf_cat(&data_buf, ",\"z\":"); buf_cat(&data_buf, mesh_z.buf); \
        buf_cat(&data_buf, ",\"i\":"); buf_cat(&data_buf, mesh_i.buf); \
        buf_cat(&data_buf, ",\"j\":"); buf_cat(&data_buf, mesh_j.buf); \
        buf_cat(&data_buf, ",\"k\":"); buf_cat(&data_buf, mesh_k.buf); \
        buf_cat(&data_buf, ",\"color\":\""); buf_cat(&data_buf, mc); \
        buf_cat(&data_buf, "\",\"flatshading\":true,\"showscale\":false}"); \
        /* Reset mesh buffers. */ \
        buf_free(&mesh_x); buf_free(&mesh_y); buf_free(&mesh_z); \
        buf_free(&mesh_i); buf_free(&mesh_j); buf_free(&mesh_k); \
        buf_init(&mesh_x, 4096); buf_init(&mesh_y, 4096); buf_init(&mesh_z, 4096); \
        buf_init(&mesh_i, 2048); buf_init(&mesh_j, 2048); buf_init(&mesh_k, 2048); \
        mesh_vcount = 0; mesh_qcount = 0; \
        buf_cat(&mesh_x, "["); buf_cat(&mesh_y, "["); buf_cat(&mesh_z, "["); \
        buf_cat(&mesh_i, "["); buf_cat(&mesh_j, "["); buf_cat(&mesh_k, "["); \
    } \
} while (0)

    /* Prime the mesh buffers with opening brackets. */
    buf_cat(&mesh_x, "["); buf_cat(&mesh_y, "["); buf_cat(&mesh_z, "[");
    buf_cat(&mesh_i, "["); buf_cat(&mesh_j, "["); buf_cat(&mesh_k, "[");

    buf_cat(&data_buf, "[");

    size_t n = prim_list->data.function.arg_count;
    for (size_t pi = 0; pi < n; pi++) {
        const Expr* p = prim_list->data.function.args[pi];
        if (!p) continue;

        /* Color directive (RGBColor / GrayLevel / Hue / CMYKColor) → flush mesh, update color. */
        {
            double r, g_c, b, a;
            if (resolve_color_rgb(p, &r, &g_c, &b, &a)) {
                FLUSH_MESH();
                cur_r = r; cur_g_val = g_c; cur_b = b;
                mesh_r = r; mesh_g_val = g_c; mesh_b_color = b;
                continue;
            }
        }

        /* Opacity → update fill opacity. */
        if (head_is(p, SYM_Opacity) && p->data.function.arg_count >= 1) {
            double op;
            if (expr_to_double(p->data.function.args[0], &op)) {
                cur_opacity = op;
                mesh_opacity = op;
            }
            continue;
        }

        /* Line[List[List[x,y,z], ...]] → scatter3d */
        if (head_is(p, SYM_Line) && p->data.function.arg_count >= 1) {
            const Expr* pts = p->data.function.args[0];
            if (!head_is(pts, SYM_List) || pts->data.function.arg_count < 2) continue;
            /* Confirm it's 3D by checking the first point has 3 elements. */
            const Expr* first_pt = pts->data.function.args[0];
            if (!head_is(first_pt, SYM_List) || first_pt->data.function.arg_count < 3) continue;

            char color_str[64];
            rgba_str(color_str, sizeof(color_str), cur_r, cur_g_val, cur_b, 1.0);

            if (!first_trace) buf_cat(&data_buf, ",");
            first_trace = 0;

            buf_cat(&data_buf, "{\"type\":\"scatter3d\",\"mode\":\"lines\",");
            buf_cat(&data_buf, "\"x\":"); append_coord3_array(&data_buf, pts, 0);
            buf_cat(&data_buf, ",\"y\":"); append_coord3_array(&data_buf, pts, 1);
            buf_cat(&data_buf, ",\"z\":"); append_coord3_array(&data_buf, pts, 2);
            buf_cat(&data_buf, ",\"line\":{\"color\":\"");
            buf_cat(&data_buf, color_str);
            buf_cat(&data_buf, "\",\"width\":3},\"showlegend\":false}");
            continue;
        }

        /* Polygon[List[pts...]] → accumulate into mesh3d.
         * Each quad (4 verts) is split into 2 triangles: (0,1,2) and (0,2,3).
         * Triangles (3 verts) are emitted as-is: (0,1,2). */
        if (head_is(p, SYM_Polygon) && p->data.function.arg_count >= 1) {
            const Expr* pts = p->data.function.args[0];
            if (!head_is(pts, SYM_List)) continue;
            size_t nv = pts->data.function.arg_count;
            if (nv < 3) continue;
            const Expr* fp = pts->data.function.args[0];
            if (!head_is(fp, SYM_List) || fp->data.function.arg_count < 3) continue;

            /* Append vertices. */
            for (size_t vi = 0; vi < nv; vi++) {
                const Expr* vp = pts->data.function.args[vi];
                if (!head_is(vp, SYM_List) || vp->data.function.arg_count < 3) continue;
                double vx, vy, vz;
                if (!expr_to_double(vp->data.function.args[0], &vx)) continue;
                if (!expr_to_double(vp->data.function.args[1], &vy)) continue;
                if (!expr_to_double(vp->data.function.args[2], &vz)) continue;
                if (mesh_vcount > 0) {
                    buf_cat(&mesh_x, ","); buf_cat(&mesh_y, ","); buf_cat(&mesh_z, ",");
                }
                buf_catd(&mesh_x, vx); buf_catd(&mesh_y, vy); buf_catd(&mesh_z, vz);
                mesh_vcount++;
            }

            /* Triangle indices. base = start of this polygon's vertices. */
            size_t base = mesh_vcount - nv;
            /* Fan: (0, k, k+1) for k in 1..nv-2. */
            for (size_t k = 1; k + 1 < nv; k++) {
                if (mesh_qcount > 0) {
                    buf_cat(&mesh_i, ","); buf_cat(&mesh_j, ","); buf_cat(&mesh_k, ",");
                }
                char ibuf[32];
                snprintf(ibuf, sizeof(ibuf), "%zu", base);     buf_cat(&mesh_i, ibuf);
                snprintf(ibuf, sizeof(ibuf), "%zu", base + k); buf_cat(&mesh_j, ibuf);
                snprintf(ibuf, sizeof(ibuf), "%zu", base + k + 1); buf_cat(&mesh_k, ibuf);
                mesh_qcount++;
            }
            continue;
        }
    }

    /* Flush any remaining mesh. */
    buf_cat(&mesh_x, "]"); buf_cat(&mesh_y, "]"); buf_cat(&mesh_z, "]");
    buf_cat(&mesh_i, "]"); buf_cat(&mesh_j, "]"); buf_cat(&mesh_k, "]");
    FLUSH_MESH();
#undef FLUSH_MESH

    buf_free(&mesh_x); buf_free(&mesh_y); buf_free(&mesh_z);
    buf_free(&mesh_i); buf_free(&mesh_j); buf_free(&mesh_k);

    buf_cat(&data_buf, "]");

    /* 3D layout: scene with equal-ish axis ranges, no legend, tight margins. */
    const char* layout =
        "{\"scene\":{"
          "\"xaxis\":{\"showgrid\":true,\"zeroline\":false},"
          "\"yaxis\":{\"showgrid\":true,\"zeroline\":false},"
          "\"zaxis\":{\"showgrid\":true,\"zeroline\":false},"
          "\"aspectmode\":\"data\""
        "},"
        "\"margin\":{\"l\":0,\"r\":0,\"t\":20,\"b\":0},"
        "\"height\":420,"
        "\"showlegend\":false,"
        "\"paper_bgcolor\":\"#fff\"}";

    Buf out;
    if (!buf_init(&out, data_buf.len + 1024)) {
        buf_free(&data_buf);
        return NULL;
    }
    buf_cat(&out, "{\"data\":");
    buf_cat(&out, data_buf.buf);
    buf_cat(&out, ",\"layout\":");
    buf_cat(&out, layout);
    buf_cat(&out, "}");

    buf_free(&data_buf);
    return out.buf;
}
