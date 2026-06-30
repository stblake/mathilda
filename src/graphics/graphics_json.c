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
 * Output JSON (Plotly format):
 *
 *   {
 *     "data": [
 *       {"type":"scatter","x":[...],"y":[...],"mode":"lines",
 *        "line":{"color":"rgba(51,102,204,1.0)","width":1.5}},
 *       ...
 *     ],
 *     "layout": {
 *       "xaxis":{"showgrid":true,"zeroline":true},
 *       "yaxis":{"showgrid":true,"zeroline":true},
 *       "margin":{"l":50,"r":20,"t":20,"b":50},
 *       "height":350,
 *       "plot_bgcolor":"#fff","paper_bgcolor":"#fff"
 *     }
 *   }
 */

#include "graphics_json.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/* Append a double value — use %g to keep it compact. */
static int buf_catd(Buf* b, double v) {
    char tmp[64];
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
        && e->data.function.head->data.symbol == sym;
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
    if (ri < 0) ri = 0; if (ri > 255) ri = 255;
    if (gi < 0) gi = 0; if (gi > 255) gi = 255;
    if (bi < 0) bi = 0; if (bi > 255) bi = 255;
    snprintf(out, outsz, "rgba(%d,%d,%d,%.3f)", ri, gi, bi, a);
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

    Buf data_buf;
    if (!buf_init(&data_buf, 65536)) return NULL;

    /* Track current drawing state. */
    double cur_r = 0.2, cur_g = 0.4, cur_b = 0.8; /* Mathematica default blue */
    double cur_opacity = 0.3;                        /* fill opacity */
    int first_trace = 1;

    buf_cat(&data_buf, "[");

    size_t n = prim_list->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        const Expr* p = prim_list->data.function.args[i];
        if (!p) continue;

        /* ----- RGBColor[r, g, b] — update current colour ------------- */
        if (head_is(p, SYM_RGBColor) && p->data.function.arg_count >= 3) {
            double r, gg, b;
            if (expr_to_double(p->data.function.args[0], &r) &&
                expr_to_double(p->data.function.args[1], &gg) &&
                expr_to_double(p->data.function.args[2], &b)) {
                cur_r = r; cur_g = gg; cur_b = b;
            }
            continue;
        }

        /* ----- Opacity[n] — update fill opacity ----------------------- */
        if (head_is(p, SYM_Opacity) && p->data.function.arg_count >= 1) {
            double op;
            if (expr_to_double(p->data.function.args[0], &op))
                cur_opacity = op;
            continue;
        }

        /* ----- Line[List[pts...]] — line trace ----------------------- */
        if (head_is(p, SYM_Line) && p->data.function.arg_count >= 1) {
            const Expr* pts = p->data.function.args[0];
            if (!head_is(pts, SYM_List)) continue;

            char color_str[64];
            rgba_str(color_str, sizeof(color_str), cur_r, cur_g, cur_b, 1.0);

            if (!first_trace) buf_cat(&data_buf, ",");
            first_trace = 0;

            buf_cat(&data_buf, "{\"type\":\"scatter\",\"mode\":\"lines\",");
            buf_cat(&data_buf, "\"x\":");
            append_coord_array(&data_buf, pts, 0);
            buf_cat(&data_buf, ",\"y\":");
            append_coord_array(&data_buf, pts, 1);
            buf_cat(&data_buf, ",\"line\":{\"color\":\"");
            buf_cat(&data_buf, color_str);
            buf_cat(&data_buf, "\",\"width\":1.5},\"showlegend\":false}");
            continue;
        }

        /* ----- Polygon[List[pts...]] — filled area ------------------- */
        if (head_is(p, SYM_Polygon) && p->data.function.arg_count >= 1) {
            const Expr* pts = p->data.function.args[0];
            if (!head_is(pts, SYM_List)) continue;

            char fill_str[64];
            rgba_str(fill_str, sizeof(fill_str), cur_r, cur_g, cur_b, cur_opacity);

            if (!first_trace) buf_cat(&data_buf, ",");
            first_trace = 0;

            buf_cat(&data_buf, "{\"type\":\"scatter\",\"mode\":\"lines\",");
            buf_cat(&data_buf, "\"fill\":\"toself\",");
            buf_cat(&data_buf, "\"fillcolor\":\"");
            buf_cat(&data_buf, fill_str);
            buf_cat(&data_buf, "\",\"x\":");
            append_coord_array(&data_buf, pts, 0);
            buf_cat(&data_buf, ",\"y\":");
            append_coord_array(&data_buf, pts, 1);
            buf_cat(&data_buf, ",\"line\":{\"color\":\"transparent\"},\"showlegend\":false}");
            continue;
        }

        /* All other primitives (Rule, $PlotResample, etc.) — skip. */
    }

    buf_cat(&data_buf, "]");

    /* Build the layout object. */
    char layout[512];
    snprintf(layout, sizeof(layout),
             "{\"xaxis\":{\"showgrid\":true,\"zeroline\":true,"
             "\"zerolinecolor\":\"#aaa\",\"zerolinewidth\":1},"
             "\"yaxis\":{\"showgrid\":true,\"zeroline\":true,"
             "\"zerolinecolor\":\"#aaa\",\"zerolinewidth\":1},"
             "\"margin\":{\"l\":50,\"r\":20,\"t\":20,\"b\":50},"
             "\"height\":350,"
             "\"plot_bgcolor\":\"#fff\","
             "\"paper_bgcolor\":\"#fff\"}");

    /* Assemble the final JSON object. */
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
    return out.buf; /* caller frees */
}
