/* render3d.c — Raylib backend for Graphics3D[...]/Plot3D[...]: graphics3d_show().
 *
 * Deliberately much smaller than render.c's 2D path: there's no Frame/
 * GridLines/Prolog/Epilog/legend/box-zoom-tool/live-resample-on-zoom here,
 * because none of them have a meaningful 3D analogue in this engine (see
 * plot3d.c's header comment) -- an orbiting camera never changes which
 * (x,y) domain is sampled, so there's nothing to re-sample. What *is*
 * shared with the 2D renderer (color resolution, numeric coercion, tick
 * spacing) comes from render_common.h instead of being re-implemented.
 *
 * Coordinate convention: the data is z-up (z is the plotted function
 * value, as in Mathematica), but Raylib's Camera3D is y-up. Every world
 * point's y and z are swapped exactly once, in to_v3() -- there is no other
 * axis remap anywhere else in this file.
 *
 * Performance note: the expression tree is converted to a flat C mesh
 * (BakedMesh) ONCE before the render loop.  This eliminates per-frame
 * malloc/free for every polygon and the per-frame recursive tree walk,
 * both of which scaled badly for dense multi-surface plots. */

#include "render3d.h"
#include "render_common.h"
#include "sym_names.h"
#include "print.h"
#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------- bounding box ---------------- */

typedef struct { double xmin, xmax, ymin, ymax, zmin, zmax; } Box3D;

static void update_bbox3(Box3D* bb, double x, double y, double z) {
    if (x < bb->xmin) bb->xmin = x;
    if (x > bb->xmax) bb->xmax = x;
    if (y < bb->ymin) bb->ymin = y;
    if (y > bb->ymax) bb->ymax = y;
    if (z < bb->zmin) bb->zmin = z;
    if (z > bb->zmax) bb->zmax = z;
}

/* Mirrors render.c's compute_bbox, walking List/Point/Line/Polygon to find
 * every {x,y,z} triple. */
static void compute_bbox3(const Expr* node, Box3D* bb) {
    if (!node || node->type != EXPR_FUNCTION) return;
    const Expr* h = node->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return;
    const char* name = h->data.symbol;
    size_t n = node->data.function.arg_count;

    if (name == SYM_List) {
        double x, y, z;
        if (n == 3 && expr_to_d(node->data.function.args[0], &x)
            && expr_to_d(node->data.function.args[1], &y)
            && expr_to_d(node->data.function.args[2], &z)) {
            update_bbox3(bb, x, y, z);
            return;
        }
        for (size_t i = 0; i < n; i++) compute_bbox3(node->data.function.args[i], bb);
    } else if ((name == SYM_Point || name == SYM_Line || name == SYM_Polygon) && n >= 1) {
        compute_bbox3(node->data.function.args[0], bb);
    }
}

/* ---------------- options ---------------- */

typedef struct {
    bool axes;
    bool lighting;          /* Lambertian shading on Polygon faces; default true */
    RGBA8 style_color;
    RGBA8 background;
    long width, height;
    const Expr* plot_label; /* borrowed */
} Gfx3DOptions;

static void gfx3d_options_parse(const Expr* graphics3d, Gfx3DOptions* o) {
    o->axes = false;
    o->lighting = true;
    o->style_color = (RGBA8){ 30, 80, 180, 255 };
    o->background = (RGBA8){ 255, 255, 255, 255 };
    o->width = 800; o->height = 600;
    o->plot_label = NULL;

    size_t argc = graphics3d->data.function.arg_count;
    for (size_t i = 1; i < argc; i++) {
        const Expr* opt = graphics3d->data.function.args[i];
        if (!opt || opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) continue;
        const Expr* h = opt->data.function.head;
        if (!h || h->type != EXPR_SYMBOL) continue;
        if (h->data.symbol != SYM_Rule && h->data.symbol != SYM_RuleDelayed) continue;
        const Expr* lhs = opt->data.function.args[0];
        const Expr* rhs = opt->data.function.args[1];
        if (lhs->type != EXPR_SYMBOL) continue;
        const char* name = lhs->data.symbol;

        if (name == SYM_Lighting) {
            /* Lighting -> None/False disables shading; anything else keeps it on. */
            o->lighting = !(rhs->type == EXPR_SYMBOL
                            && (rhs->data.symbol == SYM_None
                                || rhs->data.symbol == SYM_False));
        } else if (name == SYM_Axes) {
            o->axes = (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_True);
        } else if (name == SYM_PlotStyle) {
            resolve_color(rhs, &o->style_color);
        } else if (name == SYM_Background) {
            resolve_color(rhs, &o->background);
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
        } else if (name == SYM_PlotLabel) {
            o->plot_label = rhs;
        }
    }

    if (o->width < 100) o->width = 100;
    if (o->height < 100) o->height = 100;
}

/* ---------------- drawing ---------------- */

static Vector3 to_v3(double x, double y, double z) {
    return (Vector3){ (float)x, (float)z, (float)y };
}

/* Lambertian intensity: ambient + diffuse * |n·L|.
 * |n·L| (absolute value) gives two-sided lighting so back faces are also shaded
 * rather than receiving only the flat ambient term — important since backface
 * culling is disabled for mathematical surfaces.
 * `light_dir` is computed per-frame from the camera's local axes so that the
 * lit face always appears consistently relative to the viewer as they orbit. */
static Color shade_color(Color base, Vector3 face_normal, Vector3 light_dir) {
    const float Ka = 0.30f;  /* ambient */
    const float Kd = 0.70f;  /* diffuse */
    float ndotl = fabsf(Vector3DotProduct(face_normal, light_dir));
    float intensity = Ka + Kd * ndotl;
    return (Color){
        (unsigned char)(base.r * intensity),
        (unsigned char)(base.g * intensity),
        (unsigned char)(base.b * intensity),
        base.a,
    };
}

/* Reads a 3-coordinate {x,y,z} List into a Vector3 (already axis-remapped). */
static bool expr_point3(const Expr* e, Vector3* out) {
    double x, y, z;
    if (e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 3
        && expr_to_d(e->data.function.args[0], &x)
        && expr_to_d(e->data.function.args[1], &y)
        && expr_to_d(e->data.function.args[2], &z)) {
        *out = to_v3(x, y, z);
        return true;
    }
    return false;
}

/* ---------------- baked mesh ---------------- *
 *
 * The expression tree is walked ONCE before the render loop to produce a pair
 * of flat C arrays: triangles (BakedTri) and line segments (BakedSeg).  Each
 * frame only iterates those arrays -- no malloc, no recursive tree walk.
 *
 * Normals are pre-computed at bake time; per-frame shading only needs a dot
 * product and a multiply per triangle, not a cross-product + normalize. */

typedef struct {
    Vector3 v0, v1, v2;
    Vector3 normal;     /* pre-computed face normal */
    Color   base_color; /* unshaded colour */
} BakedTri;

typedef struct {
    Vector3 a, b;
    Color   color;
} BakedSeg;

typedef struct {
    BakedTri* tris;  size_t n_tris, tri_cap;
    BakedSeg* segs;  size_t n_segs, seg_cap;
} BakedMesh;

static void bm_push_tri(BakedMesh* bm, Vector3 v0, Vector3 v1, Vector3 v2, Color c) {
    if (bm->n_tris == bm->tri_cap) {
        bm->tri_cap = bm->tri_cap ? bm->tri_cap * 2 : 512;
        bm->tris = (BakedTri*)realloc(bm->tris, sizeof(BakedTri) * bm->tri_cap);
    }
    Vector3 e1 = Vector3Subtract(v1, v0);
    Vector3 e2 = Vector3Subtract(v2, v0);
    bm->tris[bm->n_tris++] = (BakedTri){
        v0, v1, v2,
        Vector3Normalize(Vector3CrossProduct(e1, e2)),
        c
    };
}

static void bm_push_seg(BakedMesh* bm, Vector3 a, Vector3 b, Color c) {
    if (bm->n_segs == bm->seg_cap) {
        bm->seg_cap = bm->seg_cap ? bm->seg_cap * 2 : 512;
        bm->segs = (BakedSeg*)realloc(bm->segs, sizeof(BakedSeg) * bm->seg_cap);
    }
    bm->segs[bm->n_segs++] = (BakedSeg){ a, b, c };
}

static void bm_free(BakedMesh* bm) {
    free(bm->tris); bm->tris = NULL; bm->n_tris = bm->tri_cap = 0;
    free(bm->segs); bm->segs = NULL; bm->n_segs = bm->seg_cap = 0;
}

/* Walk the expression tree and append primitives to `bm`.  Mirrors
 * draw_primitive3's logic but outputs to flat arrays instead of Raylib calls.
 * `cur_color` carries the mutable colour state across recursive calls. */
static void bake_node(const Expr* node, Color* cur_color, BakedMesh* bm) {
    if (!node || node->type != EXPR_FUNCTION) return;
    const Expr* h = node->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return;
    const char* name = h->data.symbol;
    size_t n = node->data.function.arg_count;

    if (name == SYM_List) {
        for (size_t i = 0; i < n; i++) bake_node(node->data.function.args[i], cur_color, bm);
        return;
    }
    {
        RGBA8 c;
        if (resolve_color(node, &c)) { *cur_color = to_raylib(c); return; }
    }
    if (name == SYM_Opacity && n >= 1) {
        double a;
        if (expr_to_d(node->data.function.args[0], &a)) cur_color->a = (unsigned char)(a * 255);
        return;
    }
    if (name == SYM_Polygon && n >= 1 && node->data.function.args[0]->type == EXPR_FUNCTION) {
        const Expr* arg = node->data.function.args[0];
        size_t m = arg->data.function.arg_count;
        if (m >= 3) {
            Vector3* v = (Vector3*)malloc(sizeof(Vector3) * m);
            size_t cnt = 0;
            for (size_t i = 0; i < m; i++) {
                if (expr_point3(arg->data.function.args[i], &v[cnt])) cnt++;
            }
            /* Fan triangulation from v[0] — valid because every Polygon
             * Plot3D/ParametricPlot3D emits is a convex quad. */
            for (size_t i = 1; i + 1 < cnt; i++)
                bm_push_tri(bm, v[0], v[i], v[i + 1], *cur_color);
            free(v);
        }
        return;
    }
    if (name == SYM_Line && n >= 1 && node->data.function.args[0]->type == EXPR_FUNCTION) {
        const Expr* arg = node->data.function.args[0];
        size_t m = arg->data.function.arg_count;
        Vector3 prev = { 0, 0, 0 };
        bool have_prev = false;
        for (size_t i = 0; i < m; i++) {
            Vector3 cur;
            if (expr_point3(arg->data.function.args[i], &cur)) {
                if (have_prev) bm_push_seg(bm, prev, cur, *cur_color);
                prev = cur;
                have_prev = true;
            }
        }
        return;
    }
}

/* Convert the primitive expression to a flat BakedMesh. Call once, before
 * the render loop; render with render_baked() every frame. */
static BakedMesh bake_mesh(const Expr* prims, RGBA8 default_color) {
    BakedMesh bm = { NULL, 0, 0, NULL, 0, 0 };
    Color c = to_raylib(default_color);
    bake_node(prims, &c, &bm);
    return bm;
}

/* Per-frame render: iterate the flat arrays, applying Lambertian shading when
 * lighting is on. No allocations, no tree traversal. */
static void render_baked(const BakedMesh* bm, bool lighting, Vector3 light_dir) {
    for (size_t i = 0; i < bm->n_tris; i++) {
        const BakedTri* t = &bm->tris[i];
        Color c = lighting ? shade_color(t->base_color, t->normal, light_dir) : t->base_color;
        DrawTriangle3D(t->v0, t->v1, t->v2, c);
    }
    for (size_t i = 0; i < bm->n_segs; i++) {
        DrawLine3D(bm->segs[i].a, bm->segs[i].b, bm->segs[i].color);
    }
}

static void draw_box3(const Box3D* bb, Color col) {
    Vector3 c[8] = {
        to_v3(bb->xmin, bb->ymin, bb->zmin), to_v3(bb->xmax, bb->ymin, bb->zmin),
        to_v3(bb->xmax, bb->ymax, bb->zmin), to_v3(bb->xmin, bb->ymax, bb->zmin),
        to_v3(bb->xmin, bb->ymin, bb->zmax), to_v3(bb->xmax, bb->ymin, bb->zmax),
        to_v3(bb->xmax, bb->ymax, bb->zmax), to_v3(bb->xmin, bb->ymax, bb->zmax),
    };
    static const int bot[4] = { 0, 1, 2, 3 }, top[4] = { 4, 5, 6, 7 };
    for (int i = 0; i < 4; i++) DrawLine3D(c[bot[i]], c[bot[(i + 1) % 4]], col);
    for (int i = 0; i < 4; i++) DrawLine3D(c[top[i]], c[top[(i + 1) % 4]], col);
    for (int i = 0; i < 4; i++) DrawLine3D(c[bot[i]], c[top[i]], col);
}

static void draw_tick_label3(Vector3 world_pos, Camera cam, const char* text, int win_w, int win_h, Color color) {
    Vector2 s = GetWorldToScreenEx(world_pos, cam, win_w, win_h);
    int tw = MeasureText(text, 14);
    DrawText(text, (int)s.x - tw / 2, (int)s.y, 14, color);
}

/* One tick label every nice_step() interval along each of the three box
 * edges meeting at (xmin,ymin,zmin) -- reusing render.c's tick-spacing
 * policy instead of re-deriving it. */
static void draw_box_ticks(const Box3D* bb, Camera cam, int win_w, int win_h, Color col) {
    char buf[32];
    double xstep = nice_step(bb->xmax - bb->xmin, 5);
    double ystep = nice_step(bb->ymax - bb->ymin, 5);
    double zstep = nice_step(bb->zmax - bb->zmin, 5);

    for (double x = ceil(bb->xmin / xstep) * xstep; x <= bb->xmax + 1e-9; x += xstep) {
        snprintf(buf, sizeof(buf), "%g", x);
        draw_tick_label3(to_v3(x, bb->ymin, bb->zmin), cam, buf, win_w, win_h, col);
    }
    for (double y = ceil(bb->ymin / ystep) * ystep; y <= bb->ymax + 1e-9; y += ystep) {
        snprintf(buf, sizeof(buf), "%g", y);
        draw_tick_label3(to_v3(bb->xmax, y, bb->zmin), cam, buf, win_w, win_h, col);
    }
    for (double z = ceil(bb->zmin / zstep) * zstep; z <= bb->zmax + 1e-9; z += zstep) {
        snprintf(buf, sizeof(buf), "%g", z);
        draw_tick_label3(to_v3(bb->xmin, bb->ymin, z), cam, buf, win_w, win_h, col);
    }
}

/* ---------------- 3D toolbar ---------------- *
 * Three buttons: Save / Reset / Close.  Icon glyphs are the same hand-drawn
 * vector shapes as render.c's 2D toolbar; the helpers are duplicated here
 * because they are static in render.c and the two renderers are separate
 * translation units. */

#define TB3_BTN    30.0f
#define TB3_GAP     4.0f
#define TB3_MARGIN 10.0f
#define TB3_LW      2.2f
#define TB3_COUNT   3

typedef enum { TB3_SAVE = 0, TB3_RESET, TB3_CLOSE } Tb3Btn;

static Rectangle tb3_rect(int i, int win_w) {
    float total = TB3_COUNT * TB3_BTN + (TB3_COUNT - 1) * TB3_GAP;
    float x0 = (float)win_w - TB3_MARGIN - total;
    return (Rectangle){ x0 + i * (TB3_BTN + TB3_GAP), TB3_MARGIN, TB3_BTN, TB3_BTN };
}

static int tb3_hit(Vector2 m, int win_w) {
    for (int i = 0; i < TB3_COUNT; i++)
        if (CheckCollisionPointRec(m, tb3_rect(i, win_w))) return i;
    return -1;
}

static const char* tb3_tip(int k) {
    switch (k) {
        case TB3_SAVE:  return "Save as PNG";
        case TB3_RESET: return "Reset view";
        case TB3_CLOSE: return "Close window";
        default:        return "";
    }
}

static void tb3_stroke(Vector2 a, Vector2 b, Color c) {
    DrawLineEx(a, b, TB3_LW, c);
    DrawCircleV(a, TB3_LW * 0.5f, c);
    DrawCircleV(b, TB3_LW * 0.5f, c);
}

static void tb3_icon_save(Rectangle b, Color c) {
    DrawRectangleRoundedLinesEx((Rectangle){ b.x, b.y + b.height * 0.30f,
                                             b.width, b.height * 0.52f }, 0.25f, 8, TB3_LW, c);
    DrawRectangleRoundedLinesEx((Rectangle){ b.x + b.width * 0.32f, b.y + b.height * 0.12f,
                                             b.width * 0.30f, b.height * 0.20f }, 0.4f, 6, TB3_LW, c);
    float lr = b.width * 0.16f;
    DrawRing((Vector2){ b.x + b.width * 0.5f, b.y + b.height * 0.56f },
             lr - TB3_LW * 0.5f, lr + TB3_LW * 0.5f, 0.0f, 360.0f, 32, c);
}

static void tb3_icon_home(Rectangle b, Color c) {
    float cx = b.x + b.width * 0.5f;
    float roofY = b.y + b.height * 0.12f, eaveY = b.y + b.height * 0.45f, baseY = b.y + b.height * 0.88f;
    float lx = b.x + b.width * 0.16f, rx = b.x + b.width * 0.84f;
    tb3_stroke((Vector2){ lx, eaveY }, (Vector2){ cx, roofY }, c);
    tb3_stroke((Vector2){ cx, roofY }, (Vector2){ rx, eaveY }, c);
    float wlx = b.x + b.width * 0.26f, wrx = b.x + b.width * 0.74f;
    tb3_stroke((Vector2){ wlx, eaveY }, (Vector2){ wlx, baseY }, c);
    tb3_stroke((Vector2){ wrx, eaveY }, (Vector2){ wrx, baseY }, c);
    tb3_stroke((Vector2){ wlx, baseY }, (Vector2){ wrx, baseY }, c);
    float dlx = b.x + b.width * 0.44f, drx = b.x + b.width * 0.56f, dY = b.y + b.height * 0.62f;
    tb3_stroke((Vector2){ dlx, baseY }, (Vector2){ dlx, dY }, c);
    tb3_stroke((Vector2){ drx, baseY }, (Vector2){ drx, dY }, c);
    tb3_stroke((Vector2){ dlx, dY }, (Vector2){ drx, dY }, c);
}

static void tb3_icon_close(Rectangle b, Color c) {
    float m = b.width * 0.18f;
    tb3_stroke((Vector2){ b.x + m, b.y + m },
               (Vector2){ b.x + b.width - m, b.y + b.height - m }, c);
    tb3_stroke((Vector2){ b.x + b.width - m, b.y + m },
               (Vector2){ b.x + m, b.y + b.height - m }, c);
}

static void draw_toolbar3(int win_w, int hover) {
    float total = TB3_COUNT * TB3_BTN + (TB3_COUNT - 1) * TB3_GAP;
    Rectangle panel = { (float)win_w - TB3_MARGIN - total - 5.0f, TB3_MARGIN - 5.0f,
                        total + 10.0f, TB3_BTN + 10.0f };
    DrawRectangleRounded(panel, 0.3f, 6, (Color){ 248, 248, 248, 225 });

    const Color icol = { 90, 90, 90, 255 };
    for (int i = 0; i < TB3_COUNT; i++) {
        Rectangle r = tb3_rect(i, win_w);
        if (i == hover) {
            Color hb = (i == TB3_CLOSE) ? (Color){ 240, 205, 205, 255 } : (Color){ 220, 227, 236, 255 };
            DrawRectangleRounded(r, 0.3f, 6, hb);
        }
        Rectangle ic = { r.x + 6, r.y + 6, r.width - 12, r.height - 12 };
        switch (i) {
            case TB3_SAVE:  tb3_icon_save(ic, icol); break;
            case TB3_RESET: tb3_icon_home(ic, icol); break;
            case TB3_CLOSE: tb3_icon_close(ic, (i == hover) ? (Color){ 190, 60, 60, 255 } : icol); break;
            default: break;
        }
    }

    if (hover >= 0) {
        const char* t = tb3_tip(hover);
        int tw = MeasureText(t, 12);
        Rectangle r = tb3_rect(hover, win_w);
        float tx = r.x + r.width * 0.5f - (float)tw * 0.5f;
        if (tx + (float)tw + 6 > (float)win_w) tx = (float)win_w - (float)tw - 6;
        if (tx < 4) tx = 4;
        float ty = r.y + r.height + 7;
        DrawRectangle((int)tx - 5, (int)ty - 3, tw + 10, 19, (Color){ 40, 40, 40, 235 });
        DrawText(t, (int)tx, (int)ty, 12, RAYWHITE);
    }
}

/* ---------------- main loop ---------------- */

void graphics3d_show(const Expr* graphics3d_expr) {
    if (!graphics3d_expr || graphics3d_expr->type != EXPR_FUNCTION
        || graphics3d_expr->data.function.arg_count < 1) return;

    /* Same headless-CI escape hatch as graphics_show (render.c). */
    const char* no_window = getenv("MATHILDA_NO_GRAPHICS_WINDOW");
    if (no_window && no_window[0] != '\0') {
        printf("Show: graphics window suppressed (MATHILDA_NO_GRAPHICS_WINDOW set).\n");
        return;
    }

    Gfx3DOptions opts;
    gfx3d_options_parse(graphics3d_expr, &opts);
    const Expr* prims = graphics3d_expr->data.function.args[0];

    Box3D bb = { DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX };
    compute_bbox3(prims, &bb);
    if (bb.xmin > bb.xmax) { bb.xmin = -1; bb.xmax = 1; }
    if (bb.ymin > bb.ymax) { bb.ymin = -1; bb.ymax = 1; }
    if (bb.zmin > bb.zmax) { bb.zmin = -1; bb.zmax = 1; }

    Vector3 center = to_v3((bb.xmin + bb.xmax) / 2.0, (bb.ymin + bb.ymax) / 2.0, (bb.zmin + bb.zmax) / 2.0);
    double diag = sqrt(pow(bb.xmax - bb.xmin, 2) + pow(bb.ymax - bb.ymin, 2) + pow(bb.zmax - bb.zmin, 2));
    if (!(diag > 0.0)) diag = 2.0;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow((int)opts.width, (int)opts.height, "Mathilda");
    SetTargetFPS(60);

    /* Pre-bake the expression tree into flat C arrays before the render loop.
     * This converts every Polygon into (n-2) BakedTri entries and every Line
     * into BakedSeg entries — one-time O(n) work that eliminates per-frame
     * malloc/free and recursive tree traversal. */
    BakedMesh mesh = bake_mesh(prims, opts.style_color);

    Camera3D camera = { 0 };
    camera.up = (Vector3){ 0, 1, 0 };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    camera.target = center;

    /* Orbit state: spherical coordinates around `camera.target`. A modest
     * default elevation/azimuth gives a recognizable 3/4 view instead of
     * staring straight down an axis. */
    double azimuth = -60.0 * M_PI / 180.0;
    double elevation = 25.0 * M_PI / 180.0;
    double distance = diag * 1.6;
    const double home_az = azimuth, home_el = elevation, home_dist = distance;
    const Vector3 home_target = center;

    const Color axes_color = { 90, 90, 90, 255 };

    int tb3_hover = -1;

    while (!WindowShouldClose()) {
        /* Use the actual current window size every frame so toolbar hit
         * detection and 2D overlays stay correct if the OS resizes the window
         * (e.g. window manager placement, Mission Control animations). */
        int win_w = GetScreenWidth();
        int win_h = GetScreenHeight();

        Vector2 mouse = GetMousePosition();
        Vector2 mdelta = GetMouseDelta();
        tb3_hover = tb3_hit(mouse, win_w);

        /* Toolbar clicks: only act on a release inside the button so that
         * a drag starting on the toolbar doesn't orbit the camera. */
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) && tb3_hover >= 0) {
            switch (tb3_hover) {
                case TB3_SAVE:  TakeScreenshot("mathilda_plot.png"); break;
                case TB3_RESET:
                    azimuth = home_az; elevation = home_el;
                    distance = home_dist; camera.target = home_target;
                    break;
                case TB3_CLOSE: goto done3d;
            }
        }

        /* Only orbit/pan when the mouse is NOT over the toolbar. */
        if (tb3_hover < 0) {
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                azimuth   -= mdelta.x * 0.005;
                elevation += mdelta.y * 0.005;
                const double lim = 89.0 * M_PI / 180.0;
                if (elevation > lim) elevation = lim;
                if (elevation < -lim) elevation = -lim;
            } else if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) || IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
                /* Pan the orbit target along the camera's screen-aligned
                 * right/up axes, scaled by distance so the pan speed feels
                 * consistent regardless of zoom. */
                Vector3 fwd = Vector3Normalize(Vector3Subtract(camera.target, camera.position));
                Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, camera.up));
                Vector3 trueUp = Vector3CrossProduct(right, fwd);
                float k = (float)(distance * 0.0015);
                camera.target = Vector3Add(camera.target,
                    Vector3Add(Vector3Scale(right, -mdelta.x * k), Vector3Scale(trueUp, mdelta.y * k)));
            }
        }

        float wheel = GetMouseWheelMove();
        if (wheel != 0.0f && tb3_hover < 0) {
            distance *= (wheel > 0) ? (1.0 / 1.1) : 1.1;
            if (distance < diag * 0.05) distance = diag * 0.05;
            if (distance > diag * 20.0) distance = diag * 20.0;
        }

        if (IsKeyPressed(KEY_R)) { azimuth = home_az; elevation = home_el; distance = home_dist; camera.target = home_target; }
        if (IsKeyPressed(KEY_S)) TakeScreenshot("mathilda_plot.png");
        if (IsKeyPressed(KEY_ESCAPE)) break;

        camera.position = (Vector3){
            (float)(camera.target.x + distance * cos(elevation) * cos(azimuth)),
            (float)(camera.target.y + distance * sin(elevation)),
            (float)(camera.target.z + distance * cos(elevation) * sin(azimuth)),
        };

        BeginDrawing();
        ClearBackground(to_raylib(opts.background));

        /* Camera-relative light direction: fixed in view space so the shading
         * updates as the user orbits.  Light is at upper-right-front in camera
         * space ({0.3, 1.0, 0.5} right/up/forward), transformed to world space
         * via the camera's local axes, then normalised. */
        Vector3 fwd   = Vector3Normalize(Vector3Subtract(camera.target, camera.position));
        Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, camera.up));
        Vector3 cam_up = Vector3CrossProduct(right, fwd);
        Vector3 light_dir = Vector3Normalize(
            Vector3Add(Vector3Add(Vector3Scale(right,  0.3f),
                                  Vector3Scale(cam_up,  1.0f)),
                                  Vector3Scale(fwd,     0.5f)));

        BeginMode3D(camera);
        rlDisableBackfaceCulling(); /* surfaces are visible from both sides */
        render_baked(&mesh, opts.lighting, light_dir);
        if (opts.axes) draw_box3(&bb, axes_color);
        EndMode3D(); /* EndMode3D flushes the render batch; rlDisableBackfaceCulling
                      * must stay in effect until here so the flush renders both faces */
        rlEnableBackfaceCulling();

        if (opts.axes) draw_box_ticks(&bb, camera, win_w, win_h, axes_color);
        if (opts.plot_label) {
            char* s = expr_to_string((Expr*)opts.plot_label);
            if (s) {
                int tw = MeasureText(s, 18);
                DrawText(s, (win_w - tw) / 2, 10, 18, BLACK);
                free(s);
            }
        }

        draw_toolbar3(win_w, tb3_hover);

        DrawText("drag: rotate   scroll: zoom   right-drag: pan",
                 10, win_h - 22, 14, GRAY);

        EndDrawing();
    }
    done3d:;

    bm_free(&mesh);
    CloseWindow();
}
