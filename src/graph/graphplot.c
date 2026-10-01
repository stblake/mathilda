/* graphplot.c - GraphPlot[g, opts]: draw a graph as a Graphics[...] expression,
 * plus the drawing helpers GraphPlot shares with HypergraphPlot (glayout.h).
 *
 * PIPELINE
 *   1. Layout (glayout.c): GraphLayout -> Automatic (tidy trees for branching
 *      forests, layered DAGs, stress majorization otherwise), or an explicit
 *      embedding; VertexCoordinates overrides any subset of the positions.
 *   2. Size: vertex disks get a radius scaled to the median edge length and the
 *      layout extent, so a 10-vertex and a 1000-vertex graph both read.
 *   3. Frame (gd_frame): labels are placed beside their vertex on the side
 *      with the widest angular gap between incident edges, and the world box
 *      and page size are solved together so labels are never clipped.
 *   4. Emit: edges (Line, or Arrow shortened to stop at the target disk, with
 *      Arrowheads sized to the disks; mutual pairs u->v, v->u offset apart),
 *      then vertex disks with a thin darker rim, then labels.
 *
 * The output uses only primitives every renderer draws -- Line, Arrow, Disk,
 * Circle, Text, colour and Thickness/Arrowheads directives -- with
 * AspectRatio -> Automatic, Axes -> False and an explicit PlotRange/ImageSize.
 * Colours: vertices RGBColor[0.368417, 0.506779, 0.709798] (ColorData[97][1]),
 * edges a medium grey-blue, GraphHighlight in red and thicker, as Mathematica.
 *
 * Memory (SPEC section 4): returns a freshly-allocated Graphics tree; the
 * evaluator frees res. Nothing borrowed from res outlives the call.
 */

#include "graph.h"
#include "glayout.h"
#include "graphics_export.h"
#include "expr.h"
#include "eval.h"
#include "print.h"
#include "sym_names.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Page geometry of the PDF exporter with Axes -> False (see
 * graphics_export.c): 8pt left/bottom and 12pt right/top margins. */
#define PAGE_MARGIN_X 20.0
#define PAGE_MARGIN_Y 20.0
#define LABEL_FS 10.0          /* exporter's default Text size, in points     */
#define EDGE_R 0.571589        /* default edge colour: a medium grey-blue     */
#define EDGE_G 0.586483
#define EDGE_B 0.699215

/* =========================================================== shared helpers */

static bool is_head(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

static bool is_rule(const Expr* e) {
    return (is_head(e, SYM_Rule) || is_head(e, SYM_RuleDelayed))
        && e->data.function.arg_count == 2;
}

static bool num(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (is_head(e, SYM_Rational) && e->data.function.arg_count == 2) {
        double a, b;
        if (num(e->data.function.args[0], &a) && num(e->data.function.args[1], &b) && b != 0) {
            *out = a / b; return true;
        }
    }
    return false;
}

static bool pair(const Expr* e, double* x, double* y) {
    return is_head(e, SYM_List) && e->data.function.arg_count == 2
        && num(e->data.function.args[0], x) && num(e->data.function.args[1], y);
}

const Expr* gd_option(const Expr* res, size_t first, const char* name) {
    for (size_t i = first; i < res->data.function.arg_count; i++) {
        const Expr* a = res->data.function.args[i];
        if (!is_rule(a)) continue;
        const Expr* lhs = a->data.function.args[0];
        if (lhs->type == EXPR_SYMBOL && strcmp(lhs->data.symbol.name, name) == 0)
            return a->data.function.args[1];
    }
    return NULL;
}

int gd_apply_vertex_coordinates(const Expr* spec, int n, double* xy,
                                int (*vindex)(const void* ctx, const Expr* v),
                                const void* ctx) {
    if (!is_head(spec, SYM_List)) return 0;
    size_t k = spec->data.function.arg_count;
    if (k == 0) return 0;
    double x, y;
    if (!is_rule(spec->data.function.args[0])) {
        /* {{x,y}, ...}: exactly one pair per vertex, all numeric. */
        if ((int)k != n) return 0;
        for (size_t i = 0; i < k; i++)
            if (!pair(spec->data.function.args[i], &x, &y)) return 0;
        for (size_t i = 0; i < k; i++) {
            pair(spec->data.function.args[i], &x, &y);
            xy[2 * i] = x; xy[2 * i + 1] = y;
        }
        return n;
    }
    int set = 0;
    for (size_t i = 0; i < k; i++) {
        const Expr* r = spec->data.function.args[i];
        if (!is_rule(r) || !pair(r->data.function.args[1], &x, &y)) continue;
        int v = vindex(ctx, r->data.function.args[0]);
        if (v < 0 || v >= n) continue;
        xy[2 * v] = x; xy[2 * v + 1] = y; set++;
    }
    return set;
}

void gd_push(GDPrims* P, Expr* e) {
    if (!e) { P->oom = 1; return; }
    if (P->n == P->cap) {
        size_t nc = P->cap ? 2 * P->cap : 64;
        Expr** np = realloc(P->p, nc * sizeof(Expr*));
        if (!np) { expr_free(e); P->oom = 1; return; }
        P->p = np; P->cap = nc;
    }
    P->p[P->n++] = e;
}

void gd_prims_free(GDPrims* P) {
    for (size_t i = 0; i < P->n; i++) expr_free(P->p[i]);
    free(P->p);
    P->p = NULL; P->n = P->cap = 0;
}

/* Rounds to 6 significant decimals so the Graphics tree prints compactly and
 * identically across runs. */
static double tidy(double v) {
    if (v == 0 || !isfinite(v)) return 0;
    double a = fabs(v);
    double sc = pow(10.0, 5 - (int)floor(log10(a)));
    double r = floor(v * sc + 0.5) / sc;
    return r == 0 ? 0 : r;
}

Expr* gd_pt(double x, double y) {
    Expr* a[2] = { expr_new_real(tidy(x)), expr_new_real(tidy(y)) };
    return expr_new_function(expr_new_symbol(SYM_List), a, 2);
}

Expr* gd_rgb(double r, double g, double b) {
    Expr* a[3] = { expr_new_real(tidy(r)), expr_new_real(tidy(g)), expr_new_real(tidy(b)) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
}

Expr* gd_head1(const char* head, Expr* a) {
    Expr* v[1] = { a };
    return expr_new_function(expr_new_symbol(head), v, 1);
}

Expr* gd_head2(const char* head, Expr* a, Expr* b) {
    Expr* v[2] = { a, b };
    return expr_new_function(expr_new_symbol(head), v, 2);
}

const Expr* gd_style_for(const Expr* spec, const void* item,
                         int (*match)(const Expr* key, const void* item)) {
    if (!spec) return NULL;
    if (is_rule(spec))
        return match(spec->data.function.args[0], item) ? spec->data.function.args[1] : NULL;
    if (is_head(spec, SYM_List)) {
        for (size_t i = 0; i < spec->data.function.arg_count; i++) {
            const Expr* r = spec->data.function.args[i];
            if (is_rule(r) && match(r->data.function.args[0], item))
                return r->data.function.args[1];
        }
        /* A List with no rule-like element is a list of directives: a global
         * style. Anything else (a rule list that did not match, or a
         * malformed one) gives no style. */
        for (size_t i = 0; i < spec->data.function.arg_count; i++) {
            const Expr* a = spec->data.function.args[i];
            if (is_rule(a) || is_head(a, SYM_TwoWayRule) || is_head(a, SYM_UndirectedEdge)
                || is_head(a, SYM_DirectedEdge))
                return NULL;
        }
        return spec;
    }
    if (spec->type == EXPR_SYMBOL
        && (spec->data.symbol.name == SYM_Automatic || spec->data.symbol.name == SYM_None))
        return NULL;
    return spec;
}

char* gd_label_text(const Expr* e) {
    if (e->type == EXPR_STRING) {
        size_t n = strlen(e->data.string) + 1;
        char* s = malloc(n);
        if (s) memcpy(s, e->data.string, n);
        return s;
    }
    return expr_to_string((Expr*)e);
}

GDLabelMode gd_label_mode(const Expr* v) {
    if (!v) return GD_LBL_NONE;
    if (v->type == EXPR_SYMBOL) {
        const char* n = v->data.symbol.name;
        if (n == SYM_True || n == SYM_All || n == SYM_Automatic) return GD_LBL_NAME;
        return GD_LBL_NONE;
    }
    if (v->type == EXPR_STRING) return GD_LBL_NAME;       /* "Name", "Index"... */
    if (is_rule(v) || is_head(v, SYM_List)) return GD_LBL_RULES;
    return GD_LBL_NONE;
}

const Expr* gd_rule_lookup(const Expr* rules, const Expr* key) {
    if (is_rule(rules))
        return expr_eq(rules->data.function.args[0], key) ? rules->data.function.args[1] : NULL;
    if (!is_head(rules, SYM_List)) return NULL;
    for (size_t i = 0; i < rules->data.function.arg_count; i++) {
        const Expr* r = rules->data.function.args[i];
        if (is_rule(r) && expr_eq(r->data.function.args[0], key))
            return r->data.function.args[1];
    }
    return NULL;
}

void gd_palette(int k, double* r, double* g, double* b) {
    static const double P[10][3] = {
        {0.368417, 0.506779, 0.709798}, {0.880722, 0.611041, 0.142051},
        {0.560181, 0.691569, 0.194885}, {0.922526, 0.385626, 0.209179},
        {0.528488, 0.470624, 0.701351}, {0.772079, 0.431554, 0.102387},
        {0.363898, 0.618501, 0.782349}, {1.0, 0.75, 0.0},
        {0.647624, 0.37816, 0.614037},  {0.571589, 0.586483, 0.699215},
    };
    int i = ((k % 10) + 10) % 10;
    *r = P[i][0]; *g = P[i][1]; *b = P[i][2];
}

/* ---- label placement and framing ---------------------------------------- */

/* Candidate label directions in order of preference (upper right first). */
static const double LBL_DIRS[8] = {
    M_PI / 4, 3 * M_PI / 4, -M_PI / 4, -3 * M_PI / 4, M_PI / 2, 0.0, M_PI, -M_PI / 2
};

static double ang_dist(double a, double b) {
    double d = fmod(fabs(a - b), 2 * M_PI);
    return d > M_PI ? 2 * M_PI - d : d;
}

/* Placement of one label: anchor point, Text offset, and world box. */
typedef struct { double ax, ay, ox, oy, box[4]; } LblPlace;

static void place_label(const char* text, double x, double y, double r,
                        const double* avoid, int na, double pt2w, LblPlace* L) {
    int best = 0;
    if (na > 0) {
        double bestscore = -1;
        for (int c = 0; c < 8; c++) {
            double sc = 1e9;
            for (int a = 0; a < na; a++) {
                double d = ang_dist(LBL_DIRS[c], avoid[a]);
                if (d < sc) sc = d;
            }
            if (sc > bestscore + 1e-6) { bestscore = sc; best = c; }
        }
    }
    double dx = cos(LBL_DIRS[best]), dy = sin(LBL_DIRS[best]);
    if (fabs(dx) < 1e-12) dx = 0;
    if (fabs(dy) < 1e-12) dy = 0;
    double gap = r + 1.5 * pt2w;
    L->ax = x + gap * dx; L->ay = y + gap * dy;
    double mx = fabs(dx) > fabs(dy) ? fabs(dx) : fabs(dy);
    L->ox = -dx / mx; L->oy = -dy / mx;
    double tw = graphics_helvetica_width(text) * LABEL_FS * pt2w;
    double th = 0.70 * LABEL_FS * pt2w;
    /* box centre = anchor - offset * half size */
    double cx = L->ax - L->ox * tw / 2, cy = L->ay - L->oy * th / 2;
    L->box[0] = cx - tw / 2; L->box[1] = cx + tw / 2;
    L->box[2] = cy - th / 2 - 0.25 * th;   /* allow for descenders */
    L->box[3] = cy + th / 2;
}

/* Page size and points-per-world-unit for world box W x H. */
static void page_for(double W, double H, double iw, double ih, double D,
                     double* pw, double* ph, double* k) {
    if (W < 1e-9) W = 1e-9;
    if (H < 1e-9) H = 1e-9;
    if (iw > 0 && ih > 0) {
        double kx = (iw - PAGE_MARGIN_X) / W, ky = (ih - PAGE_MARGIN_Y) / H;
        *k = kx < ky ? kx : ky; *pw = iw; *ph = ih;
    } else if (iw > 0) {
        *k = (iw - PAGE_MARGIN_X) / W; *pw = iw; *ph = H * *k + PAGE_MARGIN_Y;
    } else if (W >= H) {
        *pw = D; *k = (D - PAGE_MARGIN_X) / W; *ph = H * *k + PAGE_MARGIN_Y;
        if (*ph < 90) { *ph = 90; }
    } else {
        *ph = D; *k = (D - PAGE_MARGIN_Y) / H; *pw = W * *k + PAGE_MARGIN_X;
        if (*pw < 90) { *pw = 90; }
    }
}

void gd_frame(double* bb, const Expr* res, size_t first, int n,
              const double* xy, double r, char* const* texts,
              const int* aoff, const double* ang, GDPrims* L,
              double* pw, double* ph) {
    double iw = 0, ih = 0, v, w2;
    const Expr* is = gd_option(res, first, "ImageSize");
    if (is && num(is, &v) && v > 0) iw = v;
    else if (is && pair(is, &v, &w2) && v > 0 && w2 > 0) { iw = v; ih = w2; }
    /* Default size: 300pt, growing gently with the vertex count. */
    double D = 300.0 * sqrt(n > 40 ? n / 40.0 : 1.0);
    if (D > 600) D = 600;

    double base[4] = { bb[0], bb[1], bb[2], bb[3] };
    if (!(base[0] <= base[1])) { base[0] = -1; base[1] = 1; }
    if (!(base[2] <= base[3])) { base[2] = -1; base[3] = 1; }
    double cur[4];
    memcpy(cur, base, sizeof(cur));
    double k = 1, pwv = 0, phv = 0;
    LblPlace lp;
    for (int iter = 0; iter < 6; iter++) {
        page_for(cur[1] - cur[0], cur[3] - cur[2], iw, ih, D, &pwv, &phv, &k);
        double pt2w = 1.0 / k;
        double nb[4];
        memcpy(nb, base, sizeof(nb));
        for (int i = 0; texts && i < n; i++) {
            if (!texts[i]) continue;
            place_label(texts[i], xy[2 * i], xy[2 * i + 1], r,
                        ang ? ang + aoff[i] : NULL, ang ? aoff[i + 1] - aoff[i] : 0, pt2w, &lp);
            if (lp.box[0] < nb[0]) nb[0] = lp.box[0];
            if (lp.box[1] > nb[1]) nb[1] = lp.box[1];
            if (lp.box[2] < nb[2]) nb[2] = lp.box[2];
            if (lp.box[3] > nb[3]) nb[3] = lp.box[3];
        }
        /* Margin: 4pt plus 2% of the extent. */
        double ext = (nb[1] - nb[0]) > (nb[3] - nb[2]) ? nb[1] - nb[0] : nb[3] - nb[2];
        double mg = 4.0 * pt2w + 0.02 * ext;
        if (ext <= 1e-12) mg = 1.0;
        nb[0] -= mg; nb[1] += mg; nb[2] -= mg; nb[3] += mg;
        memcpy(cur, nb, sizeof(cur));
    }
    page_for(cur[1] - cur[0], cur[3] - cur[2], iw, ih, D, &pwv, &phv, &k);
    /* Emit the labels at the final scale. */
    if (texts) {
        int any = 0;
        for (int i = 0; i < n; i++) {
            if (!texts[i]) continue;
            if (!any) { gd_push(L, gd_head1(SYM_GrayLevel, expr_new_real(0.0))); any = 1; }
            place_label(texts[i], xy[2 * i], xy[2 * i + 1], r,
                        ang ? ang + aoff[i] : NULL, ang ? aoff[i + 1] - aoff[i] : 0, 1.0 / k, &lp);
            Expr* ta[3] = { expr_new_string(texts[i]), gd_pt(lp.ax, lp.ay), gd_pt(lp.ox, lp.oy) };
            gd_push(L, expr_new_function(expr_new_symbol(SYM_Text), ta, 3));
        }
    }
    memcpy(bb, cur, sizeof(cur));
    *pw = floor(pwv + 0.5); *ph = floor(phv + 0.5);
}

Expr* gd_finish(GDPrims* P, const double* bb, double pw, double ph,
                const Expr* res, size_t first, const char* const* consumed) {
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), P->p, P->n);
    free(P->p); P->p = NULL; P->n = P->cap = 0;
    size_t argc = res->data.function.arg_count;
    Expr** ga = malloc(sizeof(Expr*) * (argc + 6));
    if (!ga) { expr_free(list); return NULL; }
    size_t k = 0;
    ga[k++] = list;
    /* Pass-through options first: the exporter honours the first occurrence,
     * so a user's PlotRange/AspectRatio/... beats the computed one. */
    for (size_t i = first; i < argc; i++) {
        const Expr* a = res->data.function.args[i];
        if (!is_rule(a) || a->data.function.args[0]->type != EXPR_SYMBOL) continue;
        const char* nm = a->data.function.args[0]->data.symbol.name;
        bool skip = strcmp(nm, "ImageSize") == 0;
        for (const char* const* c = consumed; !skip && c && *c; c++)
            if (strcmp(nm, *c) == 0) skip = true;
        if (!skip) ga[k++] = expr_copy((Expr*)a);
    }
    ga[k++] = gd_head2(SYM_Rule, expr_new_symbol(SYM_PlotRange),
        gd_head2(SYM_List, gd_head2(SYM_List, expr_new_real(tidy(bb[0])), expr_new_real(tidy(bb[1]))),
                           gd_head2(SYM_List, expr_new_real(tidy(bb[2])), expr_new_real(tidy(bb[3])))));
    ga[k++] = gd_head2(SYM_Rule, expr_new_symbol(SYM_AspectRatio), expr_new_symbol(SYM_Automatic));
    ga[k++] = gd_head2(SYM_Rule, expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_False));
    ga[k++] = gd_head2(SYM_Rule, expr_new_symbol(SYM_ImageSize),
                       gd_head2(SYM_List, expr_new_integer((int64_t)pw), expr_new_integer((int64_t)ph)));
    Expr* g = expr_new_function(expr_new_symbol(SYM_Graphics), ga, k);
    free(ga);
    return g;
}

void gd_min_extent(double* bb, double unit) {
    if (bb[1] - bb[0] >= unit || bb[3] - bb[2] >= unit) return;
    double cx = (bb[0] + bb[1]) / 2, cy = (bb[2] + bb[3]) / 2;
    bb[0] = cx - unit / 2; bb[1] = cx + unit / 2;
    bb[2] = cy - unit / 2; bb[3] = cy + unit / 2;
}

/* ================================================================ GraphPlot */

typedef struct { const Expr* g; } VCtx;

static int graph_vindex(const void* ctx, const Expr* v) {
    return graph_vertex_position(((const VCtx*)ctx)->g, v);
}

int gd_vertex_match(const Expr* key, const void* item) {
    return expr_eq(key, (const Expr*)item);
}

/* An edge as seen by the style / highlight matchers. */
typedef struct { const Expr* u; const Expr* v; int directed; } EdgeRef;

/* Does the edge spec `key` (UndirectedEdge/TwoWayRule/DirectedEdge/Rule)
 * denote edge e? An undirected spec matches an undirected edge either way
 * round; a directed spec matches that directed edge, and (leniently) an
 * undirected edge in either orientation. */
static int edge_match(const Expr* key, const void* item) {
    const EdgeRef* e = (const EdgeRef*)item;
    if (!key || key->type != EXPR_FUNCTION || key->data.function.arg_count != 2) return 0;
    bool und = is_head(key, SYM_UndirectedEdge) || is_head(key, SYM_TwoWayRule);
    bool dir = is_head(key, SYM_DirectedEdge) || is_head(key, SYM_Rule);
    if (!und && !dir) return 0;
    const Expr* a = key->data.function.args[0];
    const Expr* b = key->data.function.args[1];
    bool fwd = expr_eq(a, e->u) && expr_eq(b, e->v);
    bool rev = expr_eq(a, e->v) && expr_eq(b, e->u);
    if (e->directed) return dir ? fwd : 0;
    return fwd || rev;
}

/* Is `item` in the GraphHighlight spec (a List of items, or one item)? */
int gd_in_highlight(const Expr* hl, const void* item,
                         int (*match)(const Expr*, const void*)) {
    if (!hl) return false;
    if (is_head(hl, SYM_List)) {
        for (size_t i = 0; i < hl->data.function.arg_count; i++)
            if (match(hl->data.function.args[i], item)) return true;
        return false;
    }
    return match(hl, item) != 0;
}

/* RGB of a colour directive, if it is one. */
static bool rgb_of(const Expr* c, double* r, double* g, double* b) {
    if (is_head(c, SYM_RGBColor) && c->data.function.arg_count >= 3)
        return num(c->data.function.args[0], r) && num(c->data.function.args[1], g)
            && num(c->data.function.args[2], b);
    if (is_head(c, SYM_GrayLevel) && c->data.function.arg_count >= 1 && num(c->data.function.args[0], r)) {
        *g = *b = *r; return true;
    }
    return false;
}

/* Pushes `style` (a directive, or a List of directives) unless it equals the
 * current one. */
static void push_style(GDPrims* P, const Expr* style, const Expr** cur) {
    if (*cur && expr_eq(*cur, style)) return;
    if (is_head(style, SYM_List))
        for (size_t i = 0; i < style->data.function.arg_count; i++)
            gd_push(P, expr_copy(style->data.function.args[i]));
    else gd_push(P, expr_copy((Expr*)style));
    *cur = style;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

/* Draws n vertices at xy with radius r: disk in its style, thin darker rim.
 * vstyle[i] may be NULL (default colour); hl[i] marks highlighted vertices. */
void gd_emit_vertices(GDPrims* P, int n, const double* xy, double r,
                      const Expr* const* vstyle, const unsigned char* hl,
                      double pt2page) {
    const Expr* cur = NULL;
    Expr* defc = gd_rgb(GD_VERTEX_R, GD_VERTEX_G, GD_VERTEX_B);
    Expr* red = gd_rgb(1, 0, 0);
    double rim = 0.6 / pt2page;              /* 0.6pt rim as a Thickness */
    gd_push(P, gd_head1(SYM_Thickness, expr_new_real(tidy(rim))));
    for (int i = 0; i < n; i++) {
        const Expr* st = (hl && hl[i]) ? red : (vstyle && vstyle[i]) ? vstyle[i] : defc;
        double rr = (hl && hl[i]) ? r * 1.15 : r;
        push_style(P, st, &cur);
        gd_push(P, gd_head2(SYM_Disk, gd_pt(xy[2 * i], xy[2 * i + 1]), expr_new_real(tidy(rr))));
        double cr, cg, cb;
        const Expr* colour = st;
        if (is_head(st, SYM_List)) {
            colour = NULL;
            for (size_t t = 0; t < st->data.function.arg_count; t++)
                if (rgb_of(st->data.function.args[t], &cr, &cg, &cb)) colour = st->data.function.args[t];
        }
        Expr* rimc = (colour && rgb_of(colour, &cr, &cg, &cb))
                   ? gd_rgb(cr * 0.6, cg * 0.6, cb * 0.6) : gd_head1(SYM_GrayLevel, expr_new_real(0.3));
        gd_push(P, rimc);
        gd_push(P, gd_head2(SYM_Circle, gd_pt(xy[2 * i], xy[2 * i + 1]), expr_new_real(tidy(rr))));
        cur = NULL;                               /* the rim colour replaced it */
    }
    expr_free(defc); expr_free(red);
}

/* Typical spacing of a drawing: median edge length, else median nearest-
 * neighbour distance, else 1. */
/* Median nearest-neighbour distance (0 when unknown / too many vertices). */
static double nn_median(int n, const double* xy) {
    if (n < 2 || n > 3000) return 0;
    double* L = malloc(sizeof(double) * (size_t)n);
    if (!L) return 0;
    int t = 0;
    for (int i = 0; i < n; i++) {
        double best = 1e300;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            double dx = xy[2 * i] - xy[2 * j], dy = xy[2 * i + 1] - xy[2 * j + 1];
            double d = dx * dx + dy * dy;
            if (d > 1e-24 && d < best) best = d;
        }
        if (best < 1e300) L[t++] = sqrt(best);
    }
    double med = 0;
    if (t > 0) { qsort(L, (size_t)t, sizeof(double), cmp_double); med = L[t / 2]; }
    free(L);
    return med;
}

double gd_drawing_unit(int n, int m, const int* eu, const int* ev, const double* xy) {
    double unit = 0;
    if (m > 0) {
        double* L = malloc(sizeof(double) * (size_t)m);
        if (L) {
            int t = 0;
            for (int k = 0; k < m; k++) {
                double dx = xy[2 * eu[k]] - xy[2 * ev[k]], dy = xy[2 * eu[k] + 1] - xy[2 * ev[k] + 1];
                double d = sqrt(dx * dx + dy * dy);
                if (d > 1e-12) L[t++] = d;
            }
            if (t > 0) { qsort(L, (size_t)t, sizeof(double), cmp_double); unit = L[t / 2]; }
            free(L);
        }
    }
    if (unit <= 0) unit = nn_median(n, xy);
    return unit > 0 ? unit : 1.0;
}

/* Vertex radius from the unit spacing and the layout extent. */
double gd_vertex_radius(int n, const double* xy, double unit, int nncap, double* bb) {
    bb[0] = bb[2] = 1e300; bb[1] = bb[3] = -1e300;
    for (int i = 0; i < n; i++) {
        if (xy[2 * i] < bb[0]) bb[0] = xy[2 * i];
        if (xy[2 * i] > bb[1]) bb[1] = xy[2 * i];
        if (xy[2 * i + 1] < bb[2]) bb[2] = xy[2 * i + 1];
        if (xy[2 * i + 1] > bb[3]) bb[3] = xy[2 * i + 1];
    }
    double ext = n ? ((bb[1] - bb[0]) > (bb[3] - bb[2]) ? bb[1] - bb[0] : bb[3] - bb[2]) : 0;
    if (ext < unit) ext = unit;
    double r = 0.075 * unit;
    if (r > 0.02 * ext) r = 0.02 * ext;
    if (r < 0.004 * ext) r = 0.004 * ext;
    /* Vertices may sit closer than edges are long (the leaves of a wide
     * tree): with nncap, never let typical neighbours' disks touch. A
     * hairball's disks may overlap, as they do in any force layout. */
    double nn = nncap ? nn_median(n, xy) : 0;
    if (nn > 0 && r > 0.3 * nn) r = 0.3 * nn;
    return r;
}

/* Builds the per-vertex label texts per VertexLabels (NULL when none). */
char** gd_vertex_texts(const Expr* spec, int n, const Expr* verts) {
    GDLabelMode mode = gd_label_mode(spec);
    if (mode == GD_LBL_NONE || n == 0) return NULL;
    char** t = calloc((size_t)n, sizeof(char*));
    if (!t) return NULL;
    for (int i = 0; i < n; i++) {
        const Expr* v = verts->data.function.args[i];
        if (mode == GD_LBL_NAME) { t[i] = gd_label_text(v); continue; }
        const Expr* l = gd_rule_lookup(spec, v);
        if (!l) continue;
        if (l->type == EXPR_STRING && strcmp(l->data.string, "Name") == 0) t[i] = gd_label_text(v);
        else if (!(l->type == EXPR_SYMBOL && l->data.symbol.name == SYM_None)) t[i] = gd_label_text(l);
    }
    return t;
}

void gd_free_texts(char** t, int n) {
    if (!t) return;
    for (int i = 0; i < n; i++) free(t[i]);
    free(t);
}

/* Options GraphPlot interprets itself (everything else passes to Graphics). */
static const char* const GP_CONSUMED[] = {
    "GraphLayout", "VertexCoordinates", "VertexLabels", "GraphHighlight",
    "VertexStyle", "EdgeStyle", "EdgeLabels", "VertexSize", "Method", NULL
};

Expr* builtin_graph_plot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;
    const Expr* g = res->data.function.args[0];
    Expr* built = NULL;
    if (is_head(g, SYM_List)) {
        /* GraphPlot[{u -> v, ...}]: the Graph of those rules. */
        Expr* call = gd_head1(SYM_Graph, expr_copy((Expr*)g));
        built = evaluate(call);
        expr_free(call);
        g = built;
    }
    if (!graph_is_valid(g)) { if (built) expr_free(built); return NULL; }
    /* Options must be Rules. */
    for (size_t i = 1; i < argc; i++)
        if (!is_rule(res->data.function.args[i])) { if (built) expr_free(built); return NULL; }

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    int m = (int)edges->data.function.arg_count;

    const int *beu = NULL, *bev = NULL; const unsigned char* bdir = NULL;
    graph_edge_indices(g, &beu, &bev, &bdir);
    int* eu = malloc(sizeof(int) * (size_t)(m + 1));
    int* ev = malloc(sizeof(int) * (size_t)(m + 1));
    unsigned char* dir = malloc((size_t)m + 1);
    double* xy = calloc(2 * (size_t)n + 2, sizeof(double));
    if (!eu || !ev || !dir || !xy) {
        free(eu); free(ev); free(dir); free(xy); if (built) expr_free(built); return NULL;
    }
    for (int k = 0; k < m; k++) { eu[k] = beu[k]; ev[k] = bev[k]; dir[k] = bdir[k]; }

    /* ---- layout -------------------------------------------------------- */
    GLMethod method = GL_AUTOMATIC;
    const Expr* lo = gd_option(res, 1, "GraphLayout");
    if (lo && !glayout_parse_method(lo, &method)) method = GL_AUTOMATIC;
    const Expr* vc = gd_option(res, 1, "VertexCoordinates");
    VCtx vctx = { g };
    int given = 0, used = GL_AUTOMATIC;
    if (vc) {
        /* Probe: a complete spec skips the layout entirely. */
        given = gd_apply_vertex_coordinates(vc, n, xy, graph_vindex, &vctx);
    }
    if (given < n) {
        used = glayout_compute(n, m, eu, ev, dir, method, xy);
        if (used < 0) {
            free(eu); free(ev); free(dir); free(xy); if (built) expr_free(built); return NULL;
        }
        if (vc) gd_apply_vertex_coordinates(vc, n, xy, graph_vindex, &vctx);
    }

    /* ---- sizes ---------------------------------------------------------- */
    double unit = gd_drawing_unit(n, m, eu, ev, xy);
    double bb[4];
    double r = gd_vertex_radius(n, xy, unit, used == GL_LAYERED || used == GL_GRID, bb);
    const Expr* vs_opt = gd_option(res, 1, "VertexSize");
    double vsz;
    if (vs_opt && num(vs_opt, &vsz) && vsz > 0) r = vsz * unit / 2;  /* diameter, in units */
    if (n > 0) { bb[0] -= 1.15 * r; bb[1] += 1.15 * r; bb[2] -= 1.15 * r; bb[3] += 1.15 * r; }
    else { bb[0] = bb[2] = -1; bb[1] = bb[3] = 1; }
    gd_min_extent(bb, unit);

    /* ---- labels & frame ------------------------------------------------- */
    char** texts = gd_vertex_texts(gd_option(res, 1, "VertexLabels"), n, verts);
    int* aoff = NULL; double* ang = NULL;
    if (texts) {
        aoff = calloc((size_t)n + 1, sizeof(int));
        ang = malloc(sizeof(double) * (size_t)(2 * m + 1));
        if (aoff && ang) {
            for (int k = 0; k < m; k++) { aoff[eu[k] + 1]++; aoff[ev[k] + 1]++; }
            for (int i = 0; i < n; i++) aoff[i + 1] += aoff[i];
            int* fill = malloc(sizeof(int) * (size_t)(n + 1));
            if (fill) {
                memcpy(fill, aoff, sizeof(int) * (size_t)n);
                for (int k = 0; k < m; k++) {
                    double dx = xy[2 * ev[k]] - xy[2 * eu[k]], dy = xy[2 * ev[k] + 1] - xy[2 * eu[k] + 1];
                    ang[fill[eu[k]]++] = atan2(dy, dx);
                    ang[fill[ev[k]]++] = atan2(-dy, -dx);
                }
                free(fill);
            } else { free(aoff); free(ang); aoff = NULL; ang = NULL; }
        } else { free(aoff); free(ang); aoff = NULL; ang = NULL; }
    }
    GDPrims labels = {0};
    double pw, ph;
    gd_frame(bb, res, 1, n, xy, r, texts, aoff, ang, &labels, &pw, &ph);
    double wpt = (pw - PAGE_MARGIN_X) / (bb[1] - bb[0]);       /* points per unit */
    double hpt = (ph - PAGE_MARGIN_Y) / (bb[3] - bb[2]);
    if (hpt < wpt) wpt = hpt;
    double plot_w_pt = (bb[1] - bb[0]) * wpt;                  /* exporter's plot_w */

    /* ---- edges ---------------------------------------------------------- */
    GDPrims P = {0};
    const Expr* hl = gd_option(res, 1, "GraphHighlight");
    const Expr* es = gd_option(res, 1, "EdgeStyle");
    const Expr* vsopt = gd_option(res, 1, "VertexStyle");
    Expr* edge_def = gd_rgb(EDGE_R, EDGE_G, EDGE_B);
    Expr* red = gd_rgb(1, 0, 0);
    double thin = 1.1 / plot_w_pt, thick = 2.8 / plot_w_pt;
    double head_len = 2.2 * r;                                 /* world units */
    if (head_len * wpt < 5.0) head_len = 5.0 / wpt;           /* >= 5pt */
    if (head_len * wpt > 12.0) head_len = 12.0 / wpt;         /* <= 12pt */
    unsigned char* ehl = calloc((size_t)m + 1, 1);
    for (int pass = 0; pass < 2 && ehl; pass++) {
        /* pass 0: ordinary edges; pass 1: highlighted edges, drawn on top */
        const Expr* cur = NULL;
        int started = 0;
        for (int k = 0; k < m; k++) {
            const Expr* e = edges->data.function.args[k];
            EdgeRef ref = { e->data.function.args[0], e->data.function.args[1], dir[k] };
            if (pass == 0) ehl[k] = gd_in_highlight(hl, &ref, edge_match);
            if ((int)ehl[k] != pass) continue;
            if (!started) {
                gd_push(&P, gd_head1(SYM_Thickness, expr_new_real(tidy(pass ? thick : thin))));
                gd_push(&P, gd_head1(SYM_Arrowheads,
                                     expr_new_real(tidy(head_len * (pass ? 1.25 : 1.0) * wpt / plot_w_pt))));
                started = 1;
            }
            const Expr* st = pass ? red : gd_style_for(es, &ref, edge_match);
            push_style(&P, st ? st : edge_def, &cur);
            double x0 = xy[2 * eu[k]], y0 = xy[2 * eu[k] + 1];
            double x1 = xy[2 * ev[k]], y1 = xy[2 * ev[k] + 1];
            if (!dir[k]) {
                gd_push(&P, gd_head1(SYM_Line, gd_head2(SYM_List, gd_pt(x0, y0), gd_pt(x1, y1))));
                continue;
            }
            double dx = x1 - x0, dy = y1 - y0, len = sqrt(dx * dx + dy * dy);
            if (len < 1e-12) continue;
            double ux = dx / len, uy = dy / len;
            /* A mutual pair u->v, v->u: shift each to its own right. */
            if (graph_has_edge(g, e->data.function.args[1], e->data.function.args[0], 1) == 1) {
                double off = 0.8 * r;
                x0 += uy * off; y0 -= ux * off; x1 += uy * off; y1 -= ux * off;
            }
            /* Tip stops just short of the target disk (highlighted disks are
             * 15% larger, so clear that radius). */
            double re = r * 1.15 + 0.5 / wpt;
            gd_push(&P, gd_head1(SYM_Arrow, gd_head2(SYM_List,
                gd_pt(x0 + ux * r, y0 + uy * r), gd_pt(x1 - ux * re, y1 - uy * re))));
        }
    }

    /* ---- edge labels ---------------------------------------------------- */
    const Expr* el = gd_option(res, 1, "EdgeLabels");
    if (el && m > 0 && !(el->type == EXPR_SYMBOL && (el->data.symbol.name == SYM_None
                                                     || el->data.symbol.name == SYM_False))) {
        Expr* weights = NULL;
        if (el->type == EXPR_STRING && strcmp(el->data.string, "EdgeWeight") == 0)
            weights = graph_resolve_edge_weights(g);
        double cx = 0, cy = 0;
        for (int i = 0; i < n; i++) { cx += xy[2 * i]; cy += xy[2 * i + 1]; }
        if (n) { cx /= n; cy /= n; }
        int any = 0;
        for (int k = 0; k < m; k++) {
            const Expr* e = edges->data.function.args[k];
            EdgeRef ref = { e->data.function.args[0], e->data.function.args[1], dir[k] };
            const Expr* lab = NULL;
            if (weights && is_head(weights, SYM_List) && (size_t)k < weights->data.function.arg_count)
                lab = weights->data.function.args[k];
            else if (!weights) lab = gd_style_for(el, &ref, edge_match);
            if (!lab || (el->type == EXPR_STRING && !weights)) continue;
            char* s = gd_label_text(lab);
            if (!s) continue;
            double x0 = xy[2 * eu[k]], y0 = xy[2 * eu[k] + 1];
            double x1 = xy[2 * ev[k]], y1 = xy[2 * ev[k] + 1];
            double mx = (x0 + x1) / 2, my = (y0 + y1) / 2;
            double dx = x1 - x0, dy = y1 - y0, len = sqrt(dx * dx + dy * dy);
            double px = len > 0 ? -dy / len : 0, py = len > 0 ? dx / len : 1;
            if ((mx - cx) * px + (my - cy) * py > 0) { px = -px; py = -py; }   /* inward */
            if (!any) { gd_push(&P, gd_head1(SYM_GrayLevel, expr_new_real(0.25))); any = 1; }
            double off = 7.0 / wpt;
            Expr* ta[3] = { expr_new_string(s), gd_pt(mx + px * off, my + py * off), gd_pt(0, 0) };
            gd_push(&P, expr_new_function(expr_new_symbol(SYM_Text), ta, 3));
            free(s);
        }
        if (weights) expr_free(weights);
    }

    /* ---- vertices ------------------------------------------------------- */
    const Expr** vstyle = n ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    unsigned char* vhl = n ? calloc((size_t)n, 1) : NULL;
    for (int i = 0; i < n && vstyle && vhl; i++) {
        const Expr* v = verts->data.function.args[i];
        vstyle[i] = gd_style_for(vsopt, v, gd_vertex_match);
        vhl[i] = gd_in_highlight(hl, v, gd_vertex_match);
    }
    if (n) gd_emit_vertices(&P, n, xy, r, vstyle, vhl, plot_w_pt);
    for (size_t i = 0; i < labels.n; i++) gd_push(&P, labels.p[i]);
    free(labels.p);

    Expr* out = NULL;
    if (!P.oom) out = gd_finish(&P, bb, pw, ph, res, 1, GP_CONSUMED);
    gd_prims_free(&P);
    expr_free(edge_def); expr_free(red);
    free((void*)vstyle); free(vhl); free(ehl);
    gd_free_texts(texts, n); free(aoff); free(ang);
    free(eu); free(ev); free(dir); free(xy);
    if (built) expr_free(built);
    return out;
}
