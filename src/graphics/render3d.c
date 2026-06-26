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
 * axis remap anywhere else in this file. */

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
    RGBA8 style_color;
    RGBA8 background;
    long width, height;
    const Expr* plot_label; /* borrowed */
} Gfx3DOptions;

static void gfx3d_options_parse(const Expr* graphics3d, Gfx3DOptions* o) {
    o->axes = false;
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

        if (name == SYM_Axes) {
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

typedef struct { Color color; } DrawState3D;

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

/* Mirrors render.c's draw_primitive, restricted to the primitives Plot3D
 * actually emits (List, color directives, Opacity, Polygon, Line). A
 * hand-built Graphics3D[] using a 2D-only primitive (Point, Rectangle,
 * Circle, Text, ...) simply isn't drawn -- the same "unrecognized node is
 * silently skipped" contract render.c's draw_primitive already has. */
static void draw_primitive3(const Expr* node, DrawState3D* state) {
    if (!node || node->type != EXPR_FUNCTION) return;
    const Expr* h = node->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return;
    const char* name = h->data.symbol;
    size_t n = node->data.function.arg_count;

    if (name == SYM_List) {
        for (size_t i = 0; i < n; i++) draw_primitive3(node->data.function.args[i], state);
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
    if (name == SYM_Polygon && n >= 1 && node->data.function.args[0]->type == EXPR_FUNCTION) {
        const Expr* arg = node->data.function.args[0];
        size_t m = arg->data.function.arg_count;
        if (m >= 3) {
            Vector3* v = malloc(sizeof(Vector3) * m);
            size_t cnt = 0;
            for (size_t i = 0; i < m; i++) {
                if (expr_point3(arg->data.function.args[i], &v[cnt])) cnt++;
            }
            /* Fan from vertex 0, exactly Polygon's 2D contract (render.c) --
             * valid here because every Polygon Plot3D builds is a convex
             * quad. Backface culling is disabled for the whole 3D draw
             * pass (see graphics3d_show), so winding doesn't matter for
             * visibility. */
            for (size_t i = 1; i + 1 < cnt; i++) DrawTriangle3D(v[0], v[i], v[i + 1], state->color);
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
                if (have_prev) DrawLine3D(prev, cur, state->color);
                prev = cur;
                have_prev = true;
            }
        }
        return;
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

    while (!WindowShouldClose()) {
        Vector2 mdelta = GetMouseDelta();

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

        float wheel = GetMouseWheelMove();
        if (wheel != 0.0f) {
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

        BeginMode3D(camera);
        rlDisableBackfaceCulling(); /* surfaces are visible from both sides */
        DrawState3D state = { .color = to_raylib(opts.style_color) };
        draw_primitive3(prims, &state);
        if (opts.axes) draw_box3(&bb, axes_color);
        EndMode3D(); /* EndMode3D flushes the render batch; rlDisableBackfaceCulling
                      * must stay in effect until here so the flush renders both faces */
        rlEnableBackfaceCulling();

        if (opts.axes) draw_box_ticks(&bb, camera, (int)opts.width, (int)opts.height, axes_color);
        if (opts.plot_label) {
            char* s = expr_to_string((Expr*)opts.plot_label);
            if (s) {
                int tw = MeasureText(s, 18);
                DrawText(s, ((int)opts.width - tw) / 2, 10, 18, BLACK);
                free(s);
            }
        }
        DrawText("drag: rotate   scroll: zoom   right-drag: pan   R: reset   S: save   Esc: close",
                 10, (int)opts.height - 22, 14, GRAY);

        EndDrawing();
    }

    CloseWindow();
}
