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
#include "label_font.h"
#include "sym_names.h"
#include "print.h"
#include "plot_common.h"   /* raylib_verbose_enabled -- $RaylibVerbose */
#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
    const char* name = h->data.symbol.name;
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
    bool have_box_ratios;   /* BoxRatios -> {rx,ry,rz} given? (else Automatic/true-scale) */
    double box_ratios[3];   /* the {rx,ry,rz} display side-length ratios */
} Gfx3DOptions;

static void gfx3d_options_parse(const Expr* graphics3d, Gfx3DOptions* o) {
    o->axes = false;
    o->lighting = true;
    o->style_color = (RGBA8){ 30, 80, 180, 255 };
    o->background = (RGBA8){ 255, 255, 255, 255 };
    o->width = 800; o->height = 600;
    o->plot_label = NULL;
    o->have_box_ratios = false;
    o->box_ratios[0] = o->box_ratios[1] = o->box_ratios[2] = 1.0;

    size_t argc = graphics3d->data.function.arg_count;
    for (size_t i = 1; i < argc; i++) {
        const Expr* opt = graphics3d->data.function.args[i];
        if (!opt || opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) continue;
        const Expr* h = opt->data.function.head;
        if (!h || h->type != EXPR_SYMBOL) continue;
        if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed) continue;
        const Expr* lhs = opt->data.function.args[0];
        const Expr* rhs = opt->data.function.args[1];
        if (lhs->type != EXPR_SYMBOL) continue;
        const char* name = lhs->data.symbol.name;

        if (name == SYM_Lighting) {
            /* Lighting -> None/False disables shading; anything else keeps it on. */
            o->lighting = !(rhs->type == EXPR_SYMBOL
                            && (rhs->data.symbol.name == SYM_None
                                || rhs->data.symbol.name == SYM_False));
        } else if (name == SYM_Axes) {
            o->axes = (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_True);
        } else if (name == SYM_PlotStyle) {
            resolve_color(rhs, &o->style_color);
        } else if (name == SYM_Background) {
            resolve_color(rhs, &o->background);
        } else if (name == SYM_ImageSize) {
            if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL
                && rhs->data.function.head->data.symbol.name == SYM_List && rhs->data.function.arg_count == 2) {
                double w, hh;
                if (expr_to_d(rhs->data.function.args[0], &w)) o->width = (long)w;
                if (expr_to_d(rhs->data.function.args[1], &hh)) o->height = (long)hh;
            } else {
                double s;
                if (expr_to_d(rhs, &s) && s > 0) { o->width = (long)s; o->height = (long)(s * 0.75); }
            }
        } else if (name == SYM_PlotLabel) {
            o->plot_label = rhs;
        } else if (name == SYM_BoxRatios) {
            /* BoxRatios -> {rx, ry, rz}: display box side-length ratios.
             * Automatic (or any non-3-list) leaves have_box_ratios false, i.e.
             * true-scale rendering. All three must be positive reals. */
            if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->type == EXPR_SYMBOL
                && rhs->data.function.head->data.symbol.name == SYM_List
                && rhs->data.function.arg_count == 3) {
                double rx, ry, rz;
                if (expr_to_d(rhs->data.function.args[0], &rx)
                    && expr_to_d(rhs->data.function.args[1], &ry)
                    && expr_to_d(rhs->data.function.args[2], &rz)
                    && rx > 0.0 && ry > 0.0 && rz > 0.0) {
                    o->box_ratios[0] = rx;
                    o->box_ratios[1] = ry;
                    o->box_ratios[2] = rz;
                    o->have_box_ratios = true;
                }
            }
        }
    }

    if (o->width < 100) o->width = 100;
    if (o->height < 100) o->height = 100;
}

/* ---------------- drawing ---------------- */

/* Display transform: maps world (data) coordinates to box-ratio-normalised
 * display coordinates. Set once per render (view3d_configure) from the data
 * bounding box and the BoxRatios option, then applied by *every* to_v3() call,
 * so the surface vertices, the bounding box, the axis ticks and the camera
 * target all share one mapping and stay consistent. Tick *labels* still read
 * true data values — they are formatted from the untransformed bbox and only
 * their screen positions pass through to_v3.
 *
 * Identity (c=0, s=1) reproduces the pre-BoxRatios behaviour exactly, so raw
 * Graphics3D[...] with no BoxRatios renders at true scale as before. */
typedef struct { double cx, cy, cz, sx, sy, sz; } View3DXform;
static View3DXform g_view3d = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };

static Vector3 to_v3(double x, double y, double z) {
    double dx = (x - g_view3d.cx) * g_view3d.sx;
    double dy = (y - g_view3d.cy) * g_view3d.sy;
    double dz = (z - g_view3d.cz) * g_view3d.sz;
    /* y and z are swapped exactly once here: data is z-up, Raylib is y-up. */
    return (Vector3){ (float)dx, (float)dz, (float)dy };
}

/* Configure the shared display transform from the data bbox and BoxRatios.
 * `ratios` (length 3) requests a box whose side lengths are exactly
 * {rx, ry, rz} — each axis is normalised to its data span then scaled by its
 * ratio (Mathematica's BoxRatios semantics). Pass NULL for Automatic / true
 * scale (identity, centred at the origin-preserving map). Writes back the
 * transformed box centre and diagonal for camera framing, so the camera
 * distance tracks whatever box the transform actually produced. */
static void view3d_configure(const Box3D* bb, const double* ratios,
                             Vector3* center_out, double* diag_out) {
    double spanx = bb->xmax - bb->xmin;
    double spany = bb->ymax - bb->ymin;
    double spanz = bb->zmax - bb->zmin;
    if (ratios) {
        g_view3d.cx = (bb->xmin + bb->xmax) / 2.0;
        g_view3d.cy = (bb->ymin + bb->ymax) / 2.0;
        g_view3d.cz = (bb->zmin + bb->zmax) / 2.0;
        g_view3d.sx = (spanx > 1e-12) ? ratios[0] / spanx : 1.0;
        g_view3d.sy = (spany > 1e-12) ? ratios[1] / spany : 1.0;
        g_view3d.sz = (spanz > 1e-12) ? ratios[2] / spanz : 1.0;
    } else {
        g_view3d.cx = g_view3d.cy = g_view3d.cz = 0.0;
        g_view3d.sx = g_view3d.sy = g_view3d.sz = 1.0;
    }
    *center_out = to_v3((bb->xmin + bb->xmax) / 2.0,
                        (bb->ymin + bb->ymax) / 2.0,
                        (bb->zmin + bb->zmax) / 2.0);
    /* Transformed extents: the map is monotone per axis, so the transformed
     * span is simply span * scale. */
    double tx = spanx * g_view3d.sx;
    double ty = spany * g_view3d.sy;
    double tz = spanz * g_view3d.sz;
    *diag_out = sqrt(tx * tx + ty * ty + tz * tz);
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

/* scene_key_light — the default Plot3D/Graphics3D key-light direction (the
 * surface->light vector used by shade_color). Dominated by world-up (render +y,
 * which to_v3 maps from the plotted height z), so upward-facing surfaces — the
 * tops the viewer looks down on — are the brightest at EVERY camera angle,
 * including looking straight down, where a camera-up-based light collapses and
 * darkens the tops. A small camera-relative lateral (`right`) and toward-viewer
 * (`-fwd`) tilt adds surface form without darkening the tops. Shared by both
 * render paths so their lighting cannot drift apart. */
static Vector3 scene_key_light(Vector3 right, Vector3 fwd) {
    Vector3 world_up = (Vector3){ 0.0f, 1.0f, 0.0f };
    return Vector3Normalize(
        Vector3Add(Vector3Add(world_up, Vector3Scale(right, 0.30f)),
                   Vector3Scale(fwd, -0.25f)));
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
    const char* name = h->data.symbol.name;
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
    int tw = label_font_measure_px(text, 14);
    label_font_draw_px(text, (int)s.x - tw / 2, (int)s.y, 14, color);
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
    /* y = ymin is the point (xmax, ymin, zmin) -- the same 3D point (and
     * therefore the same projected screen position) as the x-loop's
     * x = xmax tick above. Drawing both would overwrite one label with
     * the other; skip the y-axis's copy and let the x-axis's tick stand
     * for that corner. */
    for (double y = ceil(bb->ymin / ystep) * ystep; y <= bb->ymax + 1e-9; y += ystep) {
        if (fabs(y - bb->ymin) < 1e-9) continue;
        snprintf(buf, sizeof(buf), "%g", y);
        draw_tick_label3(to_v3(bb->xmax, y, bb->zmin), cam, buf, win_w, win_h, col);
    }
    /* z = zmin is the point (xmin, ymin, zmin) -- the same corner as the
     * x-loop's x = xmin tick; skip it for the same reason as above. */
    for (double z = ceil(bb->zmin / zstep) * zstep; z <= bb->zmax + 1e-9; z += zstep) {
        if (fabs(z - bb->zmin) < 1e-9) continue;
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
#define TB3_COUNT   4

typedef enum { TB3_SAVE = 0, TB3_SLICE, TB3_RESET, TB3_CLOSE } Tb3Btn;

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
        case TB3_SLICE: return "Toggle cross-section slicing";
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

static void tb3_icon_slice(Rectangle b, Color c) {       /* box cut by a plane */
    Rectangle box = { b.x + b.width * 0.16f, b.y + b.height * 0.16f,
                       b.width * 0.68f, b.height * 0.68f };
    DrawRectangleLinesEx(box, TB3_LW, c);
    tb3_stroke((Vector2){ b.x, b.y + b.height * 0.5f },
               (Vector2){ b.x + b.width, b.y + b.height * 0.5f }, c);
}

static void tb3_icon_close(Rectangle b, Color c) {
    float m = b.width * 0.18f;
    tb3_stroke((Vector2){ b.x + m, b.y + m },
               (Vector2){ b.x + b.width - m, b.y + b.height - m }, c);
    tb3_stroke((Vector2){ b.x + b.width - m, b.y + m },
               (Vector2){ b.x + m, b.y + b.height - m }, c);
}

static void draw_toolbar3(int win_w, int hover, bool slicing) {
    float total = TB3_COUNT * TB3_BTN + (TB3_COUNT - 1) * TB3_GAP;
    Rectangle panel = { (float)win_w - TB3_MARGIN - total - 5.0f, TB3_MARGIN - 5.0f,
                        total + 10.0f, TB3_BTN + 10.0f };
    DrawRectangleRounded(panel, 0.3f, 6, (Color){ 248, 248, 248, 225 });

    const Color icol = { 90, 90, 90, 255 };
    for (int i = 0; i < TB3_COUNT; i++) {
        Rectangle r = tb3_rect(i, win_w);
        bool active = (i == TB3_SLICE && slicing);
        if (i == hover) {
            Color hb = (i == TB3_CLOSE) ? (Color){ 240, 205, 205, 255 } : (Color){ 220, 227, 236, 255 };
            DrawRectangleRounded(r, 0.3f, 6, hb);
        }
        else if (active) DrawRectangleRounded(r, 0.3f, 6, (Color){ 205, 222, 245, 255 });
        Rectangle ic = { r.x + 6, r.y + 6, r.width - 12, r.height - 12 };
        switch (i) {
            case TB3_SAVE:  tb3_icon_save(ic, icol); break;
            case TB3_SLICE: tb3_icon_slice(ic, icol); break;
            case TB3_RESET: tb3_icon_home(ic, icol); break;
            case TB3_CLOSE: tb3_icon_close(ic, (i == hover) ? (Color){ 190, 60, 60, 255 } : icol); break;
            default: break;
        }
    }

    if (hover >= 0) {
        const char* t = tb3_tip(hover);
        int tw = label_font_measure_px(t, 12);
        Rectangle r = tb3_rect(hover, win_w);
        float tx = r.x + r.width * 0.5f - (float)tw * 0.5f;
        if (tx + (float)tw + 6 > (float)win_w) tx = (float)win_w - (float)tw - 6;
        if (tx < 4) tx = 4;
        float ty = r.y + r.height + 7;
        DrawRectangle((int)tx - 5, (int)ty - 3, tw + 10, 19, (Color){ 40, 40, 40, 235 });
        label_font_draw_px(t, (int)tx, (int)ty, 12, RAYWHITE);
    }
}

/* ---------------- hover trace (3D) ---------------- *
 *
 * Brute-force closest ray/triangle intersection over the baked mesh --
 * there is no spatial index (render_baked() itself is a flat linear scan,
 * see above), but a default-resolution Plot3D surface bakes to roughly a
 * thousand triangles, cheap enough for a per-frame scan at interactive
 * rates. GetRayCollisionTriangle is raylib's own routine, so no hand-rolled
 * Moller-Trumbore is needed here. */
static bool pick_baked_mesh(const BakedMesh* bm, Ray ray, Vector3* hit_point) {
    bool found = false;
    float best_dist = 0.0f;
    for (size_t i = 0; i < bm->n_tris; i++) {
        RayCollision rc = GetRayCollisionTriangle(ray, bm->tris[i].v0, bm->tris[i].v1, bm->tris[i].v2);
        if (rc.hit && (!found || rc.distance < best_dist)) {
            found = true;
            best_dist = rc.distance;
            *hit_point = rc.point;
        }
    }
    return found;
}

/* ---------------- cross-section slicing (3D, standalone show only) ---------------- *
 *
 * An axis-aligned plane (X, Y, or Z in data space) cuts the baked mesh.
 * There is no shortcut via Plot3D's original sample grid: build_surface_primitives
 * only has that (i,j) structure while it's building the Graphics3D[...] result --
 * by the time it reaches a flat BakedTri array here, all that's left is an
 * unstructured triangle soup, so the intersection is genuine "marching
 * triangles": per triangle, test each vertex's signed distance to the plane,
 * and where exactly two of the three edges cross zero, linearly interpolate
 * one segment. Segments are drawn independently, with no graph-stitching
 * into a single ordered curve -- adjacent triangles share edges, so they
 * still read as one continuous outline. */

typedef struct { Vector3 a, b; } SliceSeg;

/* Data axis -> the Raylib component to_v3 puts it in (to_v3 maps data
 * (x,y,z) to Raylib (x,z,y): data x stays x, data y becomes Raylib z, data
 * z becomes Raylib y). */
static float axis_component(Vector3 v, int axis) {
    switch (axis) {
        case 0: return v.x;
        case 1: return v.z;
        default: return v.y;
    }
}

static Vector3 lerp3(Vector3 a, Vector3 b, float t) {
    return (Vector3){ a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t };
}

static void axis_range(const Box3D* bb, int axis, double* lo, double* hi) {
    switch (axis) {
        case 0: *lo = bb->xmin; *hi = bb->xmax; break;
        case 1: *lo = bb->ymin; *hi = bb->ymax; break;
        default: *lo = bb->zmin; *hi = bb->zmax; break;
    }
}

static void slice_seg_push(SliceSeg** segs, size_t* len, size_t* cap, Vector3 a, Vector3 b) {
    if (*len == *cap) {
        *cap = *cap ? *cap * 2 : 256;
        *segs = (SliceSeg*)realloc(*segs, sizeof(SliceSeg) * (*cap));
    }
    (*segs)[(*len)++] = (SliceSeg){ a, b };
}

static void slice_mesh(const BakedMesh* bm, int axis, float pos,
                        SliceSeg** segs, size_t* len, size_t* cap) {
    for (size_t i = 0; i < bm->n_tris; i++) {
        Vector3 v0 = bm->tris[i].v0, v1 = bm->tris[i].v1, v2 = bm->tris[i].v2;
        float d0 = axis_component(v0, axis) - pos;
        float d1 = axis_component(v1, axis) - pos;
        float d2 = axis_component(v2, axis) - pos;
        Vector3 pts[2];
        int n = 0;
        if ((d0 < 0.0f) != (d1 < 0.0f) && d0 != d1) pts[n++] = lerp3(v0, v1, d0 / (d0 - d1));
        if (n < 2 && (d1 < 0.0f) != (d2 < 0.0f) && d1 != d2) pts[n++] = lerp3(v1, v2, d1 / (d1 - d2));
        if (n < 2 && (d2 < 0.0f) != (d0 < 0.0f) && d2 != d0) pts[n++] = lerp3(v2, v0, d2 / (d2 - d0));
        if (n == 2) slice_seg_push(segs, len, cap, pts[0], pts[1]);
    }
}

/* Translucent quad spanning the bbox's other two axes at the current plane
 * position -- backface culling is already off for this whole render pass
 * (mathematical surfaces are drawn double-sided), so one winding suffices. */
static void draw_slice_plane(const Box3D* bb, int axis, float pos, Color col) {
    Vector3 c[4];
    switch (axis) {
        case 0:
            c[0] = to_v3(pos, bb->ymin, bb->zmin); c[1] = to_v3(pos, bb->ymax, bb->zmin);
            c[2] = to_v3(pos, bb->ymax, bb->zmax); c[3] = to_v3(pos, bb->ymin, bb->zmax);
            break;
        case 1:
            c[0] = to_v3(bb->xmin, pos, bb->zmin); c[1] = to_v3(bb->xmax, pos, bb->zmin);
            c[2] = to_v3(bb->xmax, pos, bb->zmax); c[3] = to_v3(bb->xmin, pos, bb->zmax);
            break;
        default:
            c[0] = to_v3(bb->xmin, bb->ymin, pos); c[1] = to_v3(bb->xmax, bb->ymin, pos);
            c[2] = to_v3(bb->xmax, bb->ymax, pos); c[3] = to_v3(bb->xmin, bb->ymax, pos);
            break;
    }
    DrawTriangle3D(c[0], c[1], c[2], col);
    DrawTriangle3D(c[0], c[2], c[3], col);
}

/* Project a Raylib-space point onto the cutting plane's own 2D basis (the
 * two data axes *not* being sliced), for the flattened inset panel. */
static void free_axes_2d(Vector3 v, int slice_axis, float* fx, float* fy) {
    switch (slice_axis) {
        case 0: *fx = v.z; *fy = v.y; break; /* data (y, z) */
        case 1: *fx = v.x; *fy = v.y; break; /* data (x, z) */
        default: *fx = v.x; *fy = v.z; break; /* data (x, y) */
    }
}

#define SLICE_INSET_W 170.0f
#define SLICE_INSET_H 170.0f
#define SLICE_INSET_MARGIN 14.0f

/* The literal "CT-scan" view: the cross-section flattened into 2D, drawn as
 * a small fixed-size inset panel, autoscaled to fit. Segments are drawn as
 * disconnected lines (same reasoning as the 3D pass -- no stitching needed
 * for them to read as one continuous curve). */
static void draw_slice_inset(const SliceSeg* segs, size_t n_segs, int slice_axis,
                              int win_w, float floor_y) {
    if (n_segs == 0) return;
    float minx = 1e30f, maxx = -1e30f, miny = 1e30f, maxy = -1e30f;
    for (size_t i = 0; i < n_segs; i++) {
        float ax, ay, bx, by;
        free_axes_2d(segs[i].a, slice_axis, &ax, &ay);
        free_axes_2d(segs[i].b, slice_axis, &bx, &by);
        if (ax < minx) minx = ax;
        if (ax > maxx) maxx = ax;
        if (bx < minx) minx = bx;
        if (bx > maxx) maxx = bx;
        if (ay < miny) miny = ay;
        if (ay > maxy) maxy = ay;
        if (by < miny) miny = by;
        if (by > maxy) maxy = by;
    }
    float dw = maxx - minx, dh = maxy - miny;
    if (dw < 1e-6f) dw = 1.0f;
    if (dh < 1e-6f) dh = 1.0f;

    float box_x = (float)win_w - SLICE_INSET_MARGIN - SLICE_INSET_W;
    float box_y = floor_y - SLICE_INSET_MARGIN - SLICE_INSET_H;
    if (box_y < 40.0f) box_y = 40.0f;
    DrawRectangle((int)box_x - 4, (int)box_y - 4, (int)SLICE_INSET_W + 8, (int)SLICE_INSET_H + 8,
                  (Color){ 248, 248, 248, 235 });
    DrawRectangleLines((int)box_x - 4, (int)box_y - 4, (int)SLICE_INSET_W + 8, (int)SLICE_INSET_H + 8,
                        (Color){ 140, 140, 150, 255 });

    float pad = 12.0f;
    float sx = (SLICE_INSET_W - 2.0f * pad) / dw, sy = (SLICE_INSET_H - 2.0f * pad) / dh;
    float s = sx < sy ? sx : sy;
    float ox = box_x + SLICE_INSET_W * 0.5f - (minx + maxx) * 0.5f * s;
    float oy = box_y + SLICE_INSET_H * 0.5f + (miny + maxy) * 0.5f * s;

    for (size_t i = 0; i < n_segs; i++) {
        float ax, ay, bx, by;
        free_axes_2d(segs[i].a, slice_axis, &ax, &ay);
        free_axes_2d(segs[i].b, slice_axis, &bx, &by);
        Vector2 pa = { ox + ax * s, oy - ay * s };
        Vector2 pb = { ox + bx * s, oy - by * s };
        DrawLineEx(pa, pb, 2.0f, (Color){ 220, 60, 40, 255 });
    }
    const char* lbl = slice_axis == 0 ? "y-z cross-section"
                     : slice_axis == 1 ? "x-z cross-section"
                                       : "x-y cross-section";
    label_font_draw_px(lbl, (int)box_x, (int)(box_y - 16), 11, (Color){ 80, 80, 90, 255 });
}

/* ---------------- slice control row (axis select + position slider) ---------------- */

#define SLICE_ROW_H     30.0f
#define SLICE_PAD       10.0f
#define SLICE_AX0       54.0f
#define SLICE_AXBTN_W   22.0f
#define SLICE_AXBTN_GAP  4.0f
#define SLICE_VALUE_W   70.0f

static float slice_track_x0(void) { return SLICE_AX0 + 3.0f * (SLICE_AXBTN_W + SLICE_AXBTN_GAP) + 12.0f; }
static float slice_track_x1(int win_w) { return (float)win_w - SLICE_VALUE_W - SLICE_PAD; }

static int slice_axis_hit(Vector2 m, float row_y) {
    for (int k = 0; k < 3; k++) {
        Rectangle r = { SLICE_AX0 + k * (SLICE_AXBTN_W + SLICE_AXBTN_GAP),
                        row_y + (SLICE_ROW_H - SLICE_AXBTN_W) * 0.5f, SLICE_AXBTN_W, SLICE_AXBTN_W };
        if (CheckCollisionPointRec(m, r)) return k;
    }
    return -1;
}

static bool slice_track_hit(Vector2 m, float row_y, int win_w) {
    float x0 = slice_track_x0(), x1 = slice_track_x1(win_w);
    float mid_y = row_y + SLICE_ROW_H * 0.5f;
    return m.x >= x0 - 14.0f && m.x <= x1 + 14.0f && m.y >= mid_y - 12.0f && m.y <= mid_y + 12.0f;
}

static double slice_frac_from_mouse(float mx, int win_w) {
    float x0 = slice_track_x0(), x1 = slice_track_x1(win_w);
    float frac = (mx - x0) / (x1 - x0);
    if (frac < 0.0f) frac = 0.0f;
    if (frac > 1.0f) frac = 1.0f;
    return (double)frac;
}

static void draw_slice_row(int win_w, float row_y, int slice_axis, double slice_pos,
                            double lo, double hi, Vector2 mouse, bool dragging) {
    DrawRectangle(0, (int)row_y, win_w, (int)SLICE_ROW_H, (Color){ 235, 235, 240, 255 });
    DrawRectangle(0, (int)(row_y + SLICE_ROW_H - 1), win_w, 1, (Color){ 200, 200, 210, 255 });

    float mid_y = row_y + SLICE_ROW_H * 0.5f;
    label_font_draw_px("Slice", 6, (int)(mid_y - 7), 13, (Color){ 60, 60, 100, 255 });

    const char* names[3] = { "X", "Y", "Z" };
    for (int k = 0; k < 3; k++) {
        Rectangle r = { SLICE_AX0 + k * (SLICE_AXBTN_W + SLICE_AXBTN_GAP),
                        row_y + (SLICE_ROW_H - SLICE_AXBTN_W) * 0.5f, SLICE_AXBTN_W, SLICE_AXBTN_W };
        bool active = (k == slice_axis);
        bool hov = CheckCollisionPointRec(mouse, r);
        Color bg = active ? (Color){ 50, 100, 220, 255 }
                           : (hov ? (Color){ 205, 205, 215, 255 } : (Color){ 222, 222, 228, 255 });
        DrawRectangleRec(r, bg);
        DrawRectangleLinesEx(r, 1.0f, (Color){ 150, 150, 165, 255 });
        int tw = label_font_measure_px(names[k], 12);
        label_font_draw_px(names[k], (int)(r.x + r.width * 0.5f - (float)tw * 0.5f), (int)(r.y + r.height * 0.5f - 6),
                 12, active ? WHITE : (Color){ 40, 40, 60, 255 });
    }

    float track_x0 = slice_track_x0(), track_x1 = slice_track_x1(win_w);
    float track_len = track_x1 - track_x0;
    DrawRectangle((int)track_x0, (int)(mid_y - 2), (int)track_len + 1, 4, (Color){ 185, 185, 195, 255 });
    double span = hi - lo;
    double frac = (span > 0) ? (slice_pos - lo) / span : 0.0;
    if (frac < 0.0) frac = 0.0;
    if (frac > 1.0) frac = 1.0;
    float fill = (float)(frac * track_len);
    if (fill > 0.0f)
        DrawRectangle((int)track_x0, (int)(mid_y - 2), (int)fill + 1, 4, (Color){ 220, 80, 60, 255 });
    float hx = track_x0 + fill;
    bool hover_handle = fabsf(mouse.x - hx) < 12.0f && mouse.y >= row_y && mouse.y < row_y + SLICE_ROW_H;
    Color hc = (hover_handle || dragging) ? (Color){ 190, 50, 30, 255 } : (Color){ 220, 80, 60, 255 };
    DrawCircle((int)hx, (int)mid_y, 7.0f, hc);
    DrawCircleLines((int)hx, (int)mid_y, 7.0f, (Color){ 140, 30, 15, 255 });

    char val[32];
    snprintf(val, sizeof(val), "%.4g", slice_pos);
    int vw = label_font_measure_px(val, 12);
    label_font_draw_px(val, (int)((float)win_w - SLICE_VALUE_W / 2.0f - (float)vw / 2.0f), (int)(mid_y - 6), 12,
             (Color){ 40, 40, 80, 255 });
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
    const Expr* color_bar_data = find_color_bar(graphics3d_expr);

    Box3D bb = { DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX };
    compute_bbox3(prims, &bb);
    if (bb.xmin > bb.xmax) { bb.xmin = -1; bb.xmax = 1; }
    if (bb.ymin > bb.ymax) { bb.ymin = -1; bb.ymax = 1; }
    if (bb.zmin > bb.zmax) { bb.zmin = -1; bb.zmax = 1; }

    /* Configure the box-ratio display transform BEFORE baking the mesh (bake
     * runs to_v3 on every vertex). center/diag come back in transformed space
     * so the camera framing matches whatever box the ratios produce. */
    Vector3 center;
    double diag;
    view3d_configure(&bb, opts.have_box_ratios ? opts.box_ratios : NULL, &center, &diag);
    if (!(diag > 0.0)) diag = 2.0;

    SetTraceLogLevel(raylib_verbose_enabled() ? LOG_ALL : LOG_NONE); /* $RaylibVerbose */
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow((int)opts.width, (int)opts.height, "Mathilda");
    SetTargetFPS(60);
    label_font_load();

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

    /* Cross-section slicing state (toggled via the toolbar). */
    bool slicing = false;
    int slice_axis = 0;                          /* 0=X, 1=Y, 2=Z (data space) */
    double slice_pos = (bb.xmin + bb.xmax) / 2.0; /* current plane position */
    bool slice_dragging = false;

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
                case TB3_SLICE: slicing = !slicing; break;
                case TB3_RESET:
                    azimuth = home_az; elevation = home_el;
                    distance = home_dist; camera.target = home_target;
                    break;
                case TB3_CLOSE: goto done3d;
            }
        }

        /* Slice control row: only present/interactive while slicing is on.
         * Positioned above the hint-text line at the very bottom. */
        float slice_row_y = (float)win_h - 22.0f - SLICE_ROW_H - 4.0f;
        bool slice_row_hover = slicing && mouse.y >= slice_row_y && mouse.y < slice_row_y + SLICE_ROW_H;
        if (slicing && IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && slice_row_hover) {
            int ax_hit = slice_axis_hit(mouse, slice_row_y);
            if (ax_hit >= 0) {
                slice_axis = ax_hit;
                double lo, hi; axis_range(&bb, slice_axis, &lo, &hi);
                slice_pos = (lo + hi) * 0.5;
            } else if (slice_track_hit(mouse, slice_row_y, win_w)) {
                slice_dragging = true;
                double lo, hi; axis_range(&bb, slice_axis, &lo, &hi);
                slice_pos = lo + slice_frac_from_mouse(mouse.x, win_w) * (hi - lo);
            }
        }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) slice_dragging = false;
        if (slice_dragging) {
            double lo, hi; axis_range(&bb, slice_axis, &lo, &hi);
            slice_pos = lo + slice_frac_from_mouse(mouse.x, win_w) * (hi - lo);
        }

        /* Only orbit/pan when the mouse is NOT over the toolbar or the slice row. */
        if (tb3_hover < 0 && !slice_row_hover && !slice_dragging) {
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
        if (wheel != 0.0f && tb3_hover < 0 && !slice_row_hover) {
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

        /* Hover trace: pick the closest triangle under the cursor, gated to
         * "not over the toolbar" and "not currently orbiting/panning" so it
         * doesn't fight the drag gestures. */
        bool hv_found = false;
        Vector3 hv_pt = { 0, 0, 0 };
        if (tb3_hover < 0 && !IsMouseButtonDown(MOUSE_BUTTON_LEFT)
            && !IsMouseButtonDown(MOUSE_BUTTON_RIGHT) && !IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
            Ray ray = GetScreenToWorldRay(mouse, camera);
            hv_found = pick_baked_mesh(&mesh, ray, &hv_pt);
        }

        /* Cross-section slicing: re-intersect the mesh fresh every frame --
         * the mesh itself is static here (baked once before the loop), but
         * this keeps the same "recompute, don't cache" shape as the rest of
         * the interactive overlays and is cheap at these triangle counts. */
        SliceSeg* slice_segs = NULL;
        size_t n_slice_segs = 0, slice_cap = 0;
        if (slicing) slice_mesh(&mesh, slice_axis, (float)slice_pos, &slice_segs, &n_slice_segs, &slice_cap);

        BeginDrawing();
        ClearBackground(to_raylib(opts.background));

        /* Overhead key light (see scene_key_light): world-up dominated so the
         * surface tops stay brightest, with a camera-relative tilt (`right`,
         * `-fwd`) so the shading still gives form and updates as the user
         * orbits. */
        Vector3 fwd   = Vector3Normalize(Vector3Subtract(camera.target, camera.position));
        Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, camera.up));
        Vector3 light_dir = scene_key_light(right, fwd);

        BeginMode3D(camera);
        rlDisableBackfaceCulling(); /* surfaces are visible from both sides */
        render_baked(&mesh, opts.lighting, light_dir);
        if (opts.axes) draw_box3(&bb, axes_color);
        if (hv_found) DrawSphere(hv_pt, (float)(diag * 0.012), (Color){ 230, 60, 60, 255 });
        if (slicing) {
            draw_slice_plane(&bb, slice_axis, (float)slice_pos, (Color){ 80, 140, 220, 70 });
            for (size_t i = 0; i < n_slice_segs; i++)
                DrawLine3D(slice_segs[i].a, slice_segs[i].b, (Color){ 220, 60, 40, 255 });
        }
        EndMode3D(); /* EndMode3D flushes the render batch; rlDisableBackfaceCulling
                      * must stay in effect until here so the flush renders both faces */
        rlEnableBackfaceCulling();

        if (hv_found) {
            /* Undo to_v3's axis swap: raylib (x,y,z) -> data (x,z,y). */
            Vector2 spos = GetWorldToScreenEx(hv_pt, camera, win_w, win_h);
            char lbl[80];
            snprintf(lbl, sizeof(lbl), "(%.4g, %.4g, %.4g)", (double)hv_pt.x, (double)hv_pt.z, (double)hv_pt.y);
            int tw = label_font_measure_px(lbl, 13);
            float lx = spos.x + 10.0f, ly = spos.y - 18.0f;
            DrawRectangle((int)lx - 4, (int)ly - 3, tw + 8, 19, (Color){ 40, 40, 40, 215 });
            label_font_draw_px(lbl, (int)lx, (int)ly, 13, RAYWHITE);
        }

        if (opts.axes) draw_box_ticks(&bb, camera, win_w, win_h, axes_color);
        if (opts.plot_label) {
            char* s = expr_to_string((Expr*)opts.plot_label);
            if (s) {
                int tw = label_font_measure_px(s, 18);
                label_font_draw_px(s, (win_w - tw) / 2, 10, 18, BLACK);
                free(s);
            }
        }

        /* Vertical phase color-scale bar (ComplexPlot3D PlotLegends). */
        if (color_bar_data && color_bar_data->data.function.arg_count >= 2) {
            double cb_min = 0.0, cb_max = 1.0;
            const Expr* a0 = color_bar_data->data.function.args[0];
            const Expr* a1 = color_bar_data->data.function.args[1];
            if (a0->type == EXPR_REAL) cb_min = a0->data.real;
            if (a1->type == EXPR_REAL) cb_max = a1->data.real;
            const Expr* cb_cfn = (color_bar_data->data.function.arg_count >= 3)
                                 ? color_bar_data->data.function.args[2] : NULL;
            const float cb_margin = 50.0f;
            const float cb_w = 18.0f;
            float cb_h = (float)win_h - 2.0f * cb_margin;
            if (cb_h < 20.0f) cb_h = 20.0f;
            float cb_x = (float)win_w - 70.0f;
            float cb_y = cb_margin;
            draw_color_bar(cb_x, cb_y, cb_w, cb_h, cb_min, cb_max, cb_cfn);
        }

        if (slicing) {
            draw_slice_inset(slice_segs, n_slice_segs, slice_axis, win_w, slice_row_y);
            double lo, hi; axis_range(&bb, slice_axis, &lo, &hi);
            draw_slice_row(win_w, slice_row_y, slice_axis, slice_pos, lo, hi, mouse, slice_dragging);
        }

        draw_toolbar3(win_w, tb3_hover, slicing);

        label_font_draw_px(slicing
                 ? "drag: rotate   scroll: zoom   right-drag: pan   drag slider: move slice"
                 : "drag: rotate   scroll: zoom   right-drag: pan",
                 10, win_h - 22, 14, GRAY);

        EndDrawing();
        free(slice_segs);
    }
    done3d:;

    bm_free(&mesh);
    label_font_unload();
    CloseWindow();
}

/* ---------------- embedded region rendering (Animate/Manipulate) ---------------- *
 *
 * Unlike graphics3d_show(), the caller owns the window and re-evaluates its
 * body -- and therefore hands us a freshly-baked Expr* -- every frame. The
 * only thing that must survive across calls is the orbit camera; everything
 * else (bbox, baked mesh) is cheap to recompute and has no stable identity
 * to cache against (mirrors graphics_render_in_region's 2D behaviour). */

struct Graphics3DEmbedState {
    bool     initialized;
    double   azimuth, elevation, distance;
    Vector3  target;
    double   home_az, home_el, home_dist;
    Vector3  home_target;
    bool     dragging_orbit, dragging_pan;
    RenderTexture2D rt;
    bool     rt_loaded;
    int      rt_w, rt_h;
};

Graphics3DEmbedState* graphics3d_embed_state_new(void) {
    return (Graphics3DEmbedState*)calloc(1, sizeof(Graphics3DEmbedState));
}

void graphics3d_embed_state_free(Graphics3DEmbedState* state) {
    if (!state) return;
    if (state->rt_loaded) UnloadRenderTexture(state->rt);
    free(state);
}

void graphics3d_embed_state_reset_view(Graphics3DEmbedState* state) {
    if (!state || !state->initialized) return;
    state->azimuth   = state->home_az;
    state->elevation = state->home_el;
    state->distance  = state->home_dist;
    state->target    = state->home_target;
}

void graphics3d_render_in_region(const Expr* graphics3d_expr,
                                  float rx, float ry, float rw, float rh,
                                  Graphics3DEmbedState* state) {
    if (!state || !graphics3d_expr || graphics3d_expr->type != EXPR_FUNCTION
        || graphics3d_expr->data.function.arg_count < 1) return;
    if (rw <= 0.0f || rh <= 0.0f) return;

    Gfx3DOptions opts;
    gfx3d_options_parse(graphics3d_expr, &opts);
    const Expr* prims = graphics3d_expr->data.function.args[0];
    const Expr* color_bar_data = find_color_bar(graphics3d_expr);

    Box3D bb = { DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX };
    compute_bbox3(prims, &bb);
    if (bb.xmin > bb.xmax) { bb.xmin = -1; bb.xmax = 1; }
    if (bb.ymin > bb.ymax) { bb.ymin = -1; bb.ymax = 1; }
    if (bb.zmin > bb.zmax) { bb.zmin = -1; bb.zmax = 1; }

    /* Configure the box-ratio display transform BEFORE baking the mesh (bake
     * runs to_v3 on every vertex). center/diag come back in transformed space
     * so the camera framing matches whatever box the ratios produce. */
    Vector3 center;
    double diag;
    view3d_configure(&bb, opts.have_box_ratios ? opts.box_ratios : NULL, &center, &diag);
    if (!(diag > 0.0)) diag = 2.0;

    /* Establish the home view once; later calls keep whatever the user
     * orbited to even as re-evaluation reshapes the bounding box. */
    if (!state->initialized) {
        state->azimuth  = state->home_az   = -60.0 * M_PI / 180.0;
        state->elevation = state->home_el  = 25.0 * M_PI / 180.0;
        state->distance  = state->home_dist = diag * 1.6;
        state->target    = state->home_target = center;
        state->initialized = true;
    }

    BakedMesh mesh = bake_mesh(prims, opts.style_color);

    int want_w = (int)rw, want_h = (int)rh;
    if (!state->rt_loaded || state->rt_w != want_w || state->rt_h != want_h) {
        if (state->rt_loaded) UnloadRenderTexture(state->rt);
        state->rt = LoadRenderTexture(want_w, want_h);
        state->rt_loaded = true;
        state->rt_w = want_w;
        state->rt_h = want_h;
    }

    /* ---- input: orbit / pan / zoom, gated to this region ---- */
    Vector2 mouse  = GetMousePosition();
    Vector2 mdelta = GetMouseDelta();
    bool inside = mouse.x >= rx && mouse.x <= rx + rw && mouse.y >= ry && mouse.y <= ry + rh;

    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && inside)              state->dragging_orbit = true;
    if ((IsMouseButtonPressed(MOUSE_BUTTON_RIGHT)
         || IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE)) && inside)       state->dragging_pan = true;
    if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))                       state->dragging_orbit = false;
    if (IsMouseButtonReleased(MOUSE_BUTTON_RIGHT)
        || IsMouseButtonReleased(MOUSE_BUTTON_MIDDLE))                  state->dragging_pan = false;

    Camera3D camera = { 0 };
    camera.up = (Vector3){ 0, 1, 0 };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    camera.target = state->target;

    if (state->dragging_orbit) {
        state->azimuth   -= mdelta.x * 0.005;
        state->elevation += mdelta.y * 0.005;
        const double lim = 89.0 * M_PI / 180.0;
        if (state->elevation > lim) state->elevation = lim;
        if (state->elevation < -lim) state->elevation = -lim;
    } else if (state->dragging_pan) {
        Vector3 cur_pos = (Vector3){
            (float)(camera.target.x + state->distance * cos(state->elevation) * cos(state->azimuth)),
            (float)(camera.target.y + state->distance * sin(state->elevation)),
            (float)(camera.target.z + state->distance * cos(state->elevation) * sin(state->azimuth)),
        };
        Vector3 fwd    = Vector3Normalize(Vector3Subtract(camera.target, cur_pos));
        Vector3 right  = Vector3Normalize(Vector3CrossProduct(fwd, camera.up));
        Vector3 trueUp = Vector3CrossProduct(right, fwd);
        float k = (float)(state->distance * 0.0015);
        state->target = Vector3Add(state->target,
            Vector3Add(Vector3Scale(right, -mdelta.x * k), Vector3Scale(trueUp, mdelta.y * k)));
        camera.target = state->target;
    }

    float wheel = GetMouseWheelMove();
    if (wheel != 0.0f && inside) {
        state->distance *= (wheel > 0) ? (1.0 / 1.1) : 1.1;
        if (state->distance < diag * 0.05) state->distance = diag * 0.05;
        if (state->distance > diag * 20.0) state->distance = diag * 20.0;
    }

    camera.position = (Vector3){
        (float)(camera.target.x + state->distance * cos(state->elevation) * cos(state->azimuth)),
        (float)(camera.target.y + state->distance * sin(state->elevation)),
        (float)(camera.target.z + state->distance * cos(state->elevation) * sin(state->azimuth)),
    };

    Vector3 fwd    = Vector3Normalize(Vector3Subtract(camera.target, camera.position));
    Vector3 right  = Vector3Normalize(Vector3CrossProduct(fwd, camera.up));
    Vector3 light_dir = scene_key_light(right, fwd);  /* overhead key light */

    /* Hover trace: gated to the region and to "no drag in progress" (own or
     * the caller's, per the Left/Right/Middle state also gating
     * dragging_orbit/dragging_pan above). The ray is cast in the render
     * texture's own local coordinate space (mouse - region origin), since
     * that's the viewport BeginMode3D actually renders into here. */
    bool hv_found = false;
    Vector3 hv_pt = { 0, 0, 0 };
    if (inside && !state->dragging_orbit && !state->dragging_pan
        && !IsMouseButtonDown(MOUSE_BUTTON_LEFT) && !IsMouseButtonDown(MOUSE_BUTTON_RIGHT)
        && !IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
        Vector2 local_mouse = { mouse.x - rx, mouse.y - ry };
        Ray ray = GetScreenToWorldRayEx(local_mouse, camera, want_w, want_h);
        hv_found = pick_baked_mesh(&mesh, ray, &hv_pt);
    }

    /* ---- draw into an offscreen target sized exactly rw x rh, so
     * BeginMode3D's projection gets the region's own aspect ratio rather
     * than the full window's ---- */
    BeginTextureMode(state->rt);
    ClearBackground(to_raylib(opts.background));

    BeginMode3D(camera);
    rlDisableBackfaceCulling();
    render_baked(&mesh, opts.lighting, light_dir);
    if (opts.axes) draw_box3(&bb, (Color){ 90, 90, 90, 255 });
    if (hv_found) DrawSphere(hv_pt, (float)(diag * 0.012), (Color){ 230, 60, 60, 255 });
    EndMode3D();
    rlEnableBackfaceCulling();

    if (hv_found) {
        Vector2 spos = GetWorldToScreenEx(hv_pt, camera, want_w, want_h);
        char lbl[80];
        snprintf(lbl, sizeof(lbl), "(%.4g, %.4g, %.4g)", (double)hv_pt.x, (double)hv_pt.z, (double)hv_pt.y);
        int tw = label_font_measure_px(lbl, 13);
        float lx = spos.x + 10.0f, ly = spos.y - 18.0f;
        DrawRectangle((int)lx - 4, (int)ly - 3, tw + 8, 19, (Color){ 40, 40, 40, 215 });
        label_font_draw_px(lbl, (int)lx, (int)ly, 13, RAYWHITE);
    }

    if (opts.axes) draw_box_ticks(&bb, camera, want_w, want_h, (Color){ 90, 90, 90, 255 });
    if (opts.plot_label) {
        char* s = expr_to_string((Expr*)opts.plot_label);
        if (s) {
            int tw = label_font_measure_px(s, 18);
            label_font_draw_px(s, (want_w - tw) / 2, 10, 18, BLACK);
            free(s);
        }
    }
    if (color_bar_data && color_bar_data->data.function.arg_count >= 2) {
        double cb_min = 0.0, cb_max = 1.0;
        const Expr* a0 = color_bar_data->data.function.args[0];
        const Expr* a1 = color_bar_data->data.function.args[1];
        if (a0->type == EXPR_REAL) cb_min = a0->data.real;
        if (a1->type == EXPR_REAL) cb_max = a1->data.real;
        const Expr* cb_cfn = (color_bar_data->data.function.arg_count >= 3)
                             ? color_bar_data->data.function.args[2] : NULL;
        const float cb_margin = 50.0f;
        const float cb_w = 18.0f;
        float cb_h = (float)want_h - 2.0f * cb_margin;
        if (cb_h < 20.0f) cb_h = 20.0f;
        float cb_x = (float)want_w - 70.0f;
        float cb_y = cb_margin;
        draw_color_bar(cb_x, cb_y, cb_w, cb_h, cb_min, cb_max, cb_cfn);
    }

    EndTextureMode();
    bm_free(&mesh);

    /* Blit the render texture into the caller's window rect. Render
     * textures are stored bottom-up, hence the negative source height. */
    Rectangle src = { 0, 0, (float)state->rt.texture.width, -(float)state->rt.texture.height };
    Rectangle dst = { rx, ry, rw, rh };
    DrawTexturePro(state->rt.texture, src, dst, (Vector2){ 0, 0 }, 0.0f, WHITE);
}

/* Raster (PNG/JPEG) export of a Graphics3D[...] via the offscreen 3D renderer.
 * Mirrors render.c's graphics_render_rgba but drives graphics3d_render_in_region,
 * so the file matches the on-screen 3D view. Returns a fresh RGBA8 buffer
 * (caller frees with free()) or NULL when there is no usable GUI session. */
unsigned char* graphics3d_render_rgba(const Expr* g, int w, int h,
                                       int* out_w, int* out_h) {
    if (!g || w <= 0 || h <= 0) return NULL;
    if (w > 8000) w = 8000;
    if (h > 8000) h = 8000;

    bool opened = false;
    if (!IsWindowReady()) {
        if (!gui_session_available()) return NULL;   /* headless: fail, don't crash */
        SetTraceLogLevel(raylib_verbose_enabled() ? LOG_ALL : LOG_NONE); /* $RaylibVerbose */
        SetConfigFlags(FLAG_MSAA_4X_HINT);
        InitWindow(w, h, "mathilda-export");
        if (!IsWindowReady()) return NULL;
        SetWindowState(FLAG_WINDOW_HIDDEN);
        opened = true;
    }

    label_font_load();
    Graphics3DEmbedState* st = graphics3d_embed_state_new();

    /* graphics3d_render_in_region renders into st->rt then blits it to the
     * active target; give the blit a real BeginDrawing target (discarded), then
     * read st->rt back. This mirrors how Animate/Manipulate drive it. */
    BeginDrawing();
    ClearBackground(RAYWHITE);
    graphics3d_render_in_region(g, 0.0f, 0.0f, (float)w, (float)h, st);
    EndDrawing();

    unsigned char* out = NULL;
    if (st->rt_loaded) {
        Image img = LoadImageFromTexture(st->rt.texture);
        ImageFlipVertical(&img);                       /* render textures are y-flipped */
        ImageFormat(&img, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
        size_t n = (size_t)img.width * (size_t)img.height * 4;
        if (img.data && n > 0) {
            out = (unsigned char*)malloc(n);
            if (out) {
                memcpy(out, img.data, n);
                *out_w = img.width;
                *out_h = img.height;
            }
        }
        UnloadImage(img);
    }

    graphics3d_embed_state_free(st);                    /* unloads st->rt */
    label_font_unload();
    if (opened) CloseWindow();
    return out;
}
