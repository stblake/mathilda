/* manipulate.c — Manipulate[expr, {u, umin, umax}, ...]
 *
 * Opens a Raylib window (when USE_GRAPHICS is compiled in) with one
 * control row per variable, stacked above a thin footer:
 *
 *   ┌─────────────────────────────────────────────────────────────┐
 *   │                                                               │
 *   │                        (rendered content)                    │
 *   │                                                               │
 *   ├─────────────────────────────────────────────────────────────┤
 *   │ u  │ 0 ─────────────●────────── 10              │  3.42     │
 *   │ f  │ [Sin[x]] [Cos[x]] [Tan[x]]                              │
 *   ├─────────────────────────────────────────────────────────────┤
 *   │ [Reset]  Esc: close                                          │
 *   └─────────────────────────────────────────────────────────────┘
 *
 * Each control is independently user-driven from the first frame — there
 * is no animation phase or playback transport (see Animate for that).
 * A control is either:
 *   - a continuous drag slider:  {u, umin, umax}, {u, umin, umax, du},
 *     {{u, u0}, umin, umax}, {{u, u0}, umin, umax, du}
 *   - a discrete button set:     {u, {v1, v2, ...}}, {{u, u0}, {v1, ...}}
 */

#include "manipulate.h"
#include "show.h"
#include "render.h"
#include "render3d.h"
#include "sym_names.h"
#include "symtab.h"
#include "iter.h"
#include "eval.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Per-control state                                                    */
/* ------------------------------------------------------------------ */

#define MAX_CTRLS            8
#define MAX_DISCRETE_VALUES  32

typedef enum { MCTRL_RANGE, MCTRL_DISCRETE } ManipCtrlKind;

typedef struct {
    ManipCtrlKind kind;
    Expr*  var;              /* owned symbol Expr */

    /* MCTRL_RANGE */
    double vmin, vmax, step; /* step == 0 means no snapping */
    double value;
    double default_value;

    /* MCTRL_DISCRETE */
    Expr** values;           /* owned array of owned Expr* */
    size_t n_values;
    size_t selected_idx;
    size_t default_idx;
} ManipCtrl;

static void manip_ctrl_free(ManipCtrl* c) {
    if (!c) return;
    if (c->var) expr_free(c->var);
    if (c->kind == MCTRL_DISCRETE && c->values) {
        for (size_t i = 0; i < c->n_values; i++)
            if (c->values[i]) expr_free(c->values[i]);
        free(c->values);
    }
}

/* ------------------------------------------------------------------ */
/* Control-spec parsing (independent of USE_GRAPHICS)                   */
/* ------------------------------------------------------------------ */

static bool numericize(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    Expr* a[1] = { expr_copy((Expr*)e) };
    Expr* nc   = expr_new_function(expr_new_symbol("N"), a, 1);
    Expr* r    = evaluate(nc);
    bool ok = false;
    if (r->type == EXPR_REAL)    { *out = r->data.real;            ok = true; }
    if (r->type == EXPR_INTEGER) { *out = (double)r->data.integer; ok = true; }
    expr_free(r);
    return ok;
}

static bool is_list_expr(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* First element of a control spec: either a bare symbol `var`, or an
 * explicit-default pair `{var, default}`. On success sets *var_sym and,
 * for the pair form, evaluates and returns the default in *default_out
 * (caller owns it; NULL means "no explicit default"). */
static bool parse_var_part(const Expr* first, const char** var_sym,
                            Expr** default_out) {
    *default_out = NULL;
    if (!first) return false;
    if (first->type == EXPR_SYMBOL) {
        *var_sym = first->data.symbol.name;
        return true;
    }
    if (is_list_expr(first) && first->data.function.arg_count == 2) {
        const Expr* v = first->data.function.args[0];
        if (!v || v->type != EXPR_SYMBOL) return false;
        *var_sym = v->data.symbol.name;
        *default_out = evaluate(expr_copy((Expr*)first->data.function.args[1]));
        return true;
    }
    return false;
}

/* Parse one already-evaluated control-spec List into `out`. Returns false
 * (leaving *out untouched) if `e` isn't a recognized control spec. */
static bool parse_control_spec(const Expr* e, ManipCtrl* out) {
    if (!is_list_expr(e) || e->data.function.arg_count < 2) return false;

    const char* var_sym = NULL;
    Expr* default_expr = NULL;
    if (!parse_var_part(e->data.function.args[0], &var_sym, &default_expr))
        return false;

    size_t rest_n = e->data.function.arg_count - 1;
    Expr** rest = e->data.function.args + 1;

    /* Discrete: exactly one remaining arg, itself a nonempty List. */
    if (rest_n == 1 && is_list_expr(rest[0])) {
        size_t nv = rest[0]->data.function.arg_count;
        if (nv == 0 || nv > MAX_DISCRETE_VALUES) {
            if (default_expr) expr_free(default_expr);
            return false;
        }
        out->kind = MCTRL_DISCRETE;
        out->var  = expr_new_symbol(var_sym);
        out->values = calloc(nv, sizeof(Expr*));
        out->n_values = nv;
        for (size_t i = 0; i < nv; i++)
            out->values[i] = expr_copy((Expr*)rest[0]->data.function.args[i]);
        out->selected_idx = 0;
        out->default_idx  = 0;
        if (default_expr) {
            for (size_t i = 0; i < nv; i++) {
                if (expr_eq(out->values[i], default_expr)) {
                    out->selected_idx = out->default_idx = i;
                    break;
                }
            }
            expr_free(default_expr);
        }
        return true;
    }

    /* Continuous range: 2 (min,max) or 3 (min,max,step) numeric args. */
    if (rest_n == 2 || rest_n == 3) {
        double lo, hi, step = 0.0;
        if (!numericize(rest[0], &lo) || !numericize(rest[1], &hi)) {
            if (default_expr) expr_free(default_expr);
            return false;
        }
        if (rest_n == 3 && !numericize(rest[2], &step)) {
            if (default_expr) expr_free(default_expr);
            return false;
        }
        if (lo > hi) { double t = lo; lo = hi; hi = t; }
        double def = lo;
        if (default_expr) {
            double dv;
            if (numericize(default_expr, &dv)) def = dv;
            expr_free(default_expr);
        }
        if (def < lo) def = lo;
        if (def > hi) def = hi;

        out->kind  = MCTRL_RANGE;
        out->var   = expr_new_symbol(var_sym);
        out->vmin  = lo;
        out->vmax  = hi;
        out->step  = (step > 0.0) ? step : 0.0;
        out->value = def;
        out->default_value = def;
        return true;
    }

    if (default_expr) expr_free(default_expr);
    return false;
}

/* ------------------------------------------------------------------ */
/* Rendering (Raylib)                                                    */
/* ------------------------------------------------------------------ */

#ifdef USE_GRAPHICS
#include <raylib.h>
#include "label_font.h"

#define ROW_H        30.0f   /* height of each control row                 */
#define FOOTER_H     28.0f   /* Reset button + Esc hint row                */
#define LABEL_W      56.0f   /* pixels reserved for the var-name label     */
#define VALUE_W      54.0f   /* pixels reserved for the current-value text */
#define PAD          10.0f   /* horizontal gap around the track            */
#define BTN_H        22.0f   /* discrete button height                     */
#define BTN_PAD_X     8.0f   /* horizontal padding inside a discrete button*/
#define BTN_GAP       6.0f   /* gap between discrete buttons                */
#define RESET_BTN_W  52.0f
#define RESET_BTN_H  20.0f

static bool is_graphics2d_expr(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Graphics
        && e->data.function.arg_count >= 1;
}

static bool is_graphics3d_expr(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Graphics3D
        && e->data.function.arg_count >= 1;
}

/* ---- continuous range row ---- */

static void draw_range_row(const ManipCtrl* c, int win_w, float row_y,
                            Vector2 mouse, bool dragging) {
    float track_x0 = LABEL_W + PAD;
    float track_x1 = (float)win_w - VALUE_W - PAD;
    float track_len = track_x1 - track_x0;
    float mid_y = row_y + ROW_H * 0.5f;

    DrawRectangle(0, (int)row_y, win_w, (int)ROW_H, (Color){235,235,240,255});
    DrawRectangle(0, (int)(row_y + ROW_H - 1), win_w, 1, (Color){200,200,210,255});

    label_font_draw_px(c->var->data.symbol.name, 6, (int)(mid_y - 7), 14, (Color){60,60,100,255});

    char lo[24], hi[24];
    snprintf(lo, sizeof(lo), "%.4g", c->vmin);
    snprintf(hi, sizeof(hi), "%.4g", c->vmax);
    int lo_w = label_font_measure_px(lo, 10);
    label_font_draw_px(lo, (int)(track_x0 - lo_w - 3), (int)(mid_y - 5), 10, DARKGRAY);
    label_font_draw_px(hi, (int)(track_x1 + 3),         (int)(mid_y - 5), 10, DARKGRAY);

    DrawRectangle((int)track_x0, (int)(mid_y - 2), (int)track_len + 1, 4,
                  (Color){185,185,195,255});

    double span = c->vmax - c->vmin;
    double frac = (span > 0) ? (c->value - c->vmin) / span : 0.0;
    if (frac < 0.0) frac = 0.0;
    if (frac > 1.0) frac = 1.0;

    float fill = (float)(frac * track_len);
    if (fill > 0.0f)
        DrawRectangle((int)track_x0, (int)(mid_y - 2), (int)fill + 1, 4,
                      (Color){60,110,215,255});

    float hx = track_x0 + fill;
    bool hover_handle = fabsf(mouse.x - hx) < 12.0f
                     && mouse.y >= row_y && mouse.y < row_y + ROW_H;
    Color handle_col = (hover_handle || dragging)
                       ? (Color){30,80,200,255}
                       : (Color){50,100,220,255};
    DrawCircle((int)hx, (int)mid_y, 7.5f, handle_col);
    DrawCircleLines((int)hx, (int)mid_y, 7.5f, (Color){20,60,160,255});

    char val[32];
    snprintf(val, sizeof(val), "%.5g", c->value);
    int val_w = label_font_measure_px(val, 12);
    label_font_draw_px(val, (int)((float)win_w - VALUE_W / 2.0f - (float)val_w / 2.0f),
             (int)(mid_y - 6), 12, (Color){40,40,80,255});
}

static bool range_row_hit(Vector2 m, float row_y, int win_w) {
    float track_x0 = LABEL_W + PAD;
    float track_x1 = (float)win_w - VALUE_W - PAD;
    float mid_y    = row_y + ROW_H * 0.5f;
    return m.x >= track_x0 - 14.0f && m.x <= track_x1 + 14.0f
        && m.y >= mid_y - 12.0f && m.y <= mid_y + 12.0f;
}

static double range_frac_from_mouse(float mx, int win_w) {
    float track_x0 = LABEL_W + PAD;
    float track_x1 = (float)win_w - VALUE_W - PAD;
    float frac = (mx - track_x0) / (track_x1 - track_x0);
    if (frac < 0.0f) frac = 0.0f;
    if (frac > 1.0f) frac = 1.0f;
    return (double)frac;
}

static double range_snap(const ManipCtrl* c, double v) {
    if (c->step > 0.0) {
        double n = round((v - c->vmin) / c->step);
        v = c->vmin + n * c->step;
    }
    if (v < c->vmin) v = c->vmin;
    if (v > c->vmax) v = c->vmax;
    return v;
}

/* ---- discrete button-set row ----
 * Layout is deterministic given the (immutable) option texts, so drawing
 * and hit-testing recompute the same geometry independently. */

static int discrete_row_layout(const ManipCtrl* c, float* out_x, float* out_w) {
    float x = LABEL_W + PAD;
    int n = 0;
    for (size_t i = 0; i < c->n_values; i++) {
        char* s = expr_to_string(c->values[i]);
        int tw = s ? label_font_measure_px(s, 12) : 8;
        if (s) free(s);
        float w = (float)tw + 2.0f * BTN_PAD_X;
        out_x[n] = x;
        out_w[n] = w;
        n++;
        x += w + BTN_GAP;
    }
    return n;
}

static void draw_discrete_row(const ManipCtrl* c, int win_w, float row_y,
                               Vector2 mouse) {
    float mid_y = row_y + ROW_H * 0.5f;
    DrawRectangle(0, (int)row_y, win_w, (int)ROW_H, (Color){235,235,240,255});
    DrawRectangle(0, (int)(row_y + ROW_H - 1), win_w, 1, (Color){200,200,210,255});
    label_font_draw_px(c->var->data.symbol.name, 6, (int)(mid_y - 7), 14, (Color){60,60,100,255});

    float bx[MAX_DISCRETE_VALUES], bw[MAX_DISCRETE_VALUES];
    int n = discrete_row_layout(c, bx, bw);
    float by = row_y + (ROW_H - BTN_H) * 0.5f;
    for (int i = 0; i < n; i++) {
        char* s = expr_to_string(c->values[i]);
        const char* text = s ? s : "?";
        bool selected = ((size_t)i == c->selected_idx);
        bool hovered = mouse.x >= bx[i] && mouse.x <= bx[i] + bw[i]
                     && mouse.y >= by && mouse.y <= by + BTN_H;
        Color bg = selected ? (Color){50,100,220,255}
                             : (hovered ? (Color){205,205,215,255} : (Color){222,222,228,255});
        DrawRectangleRec((Rectangle){bx[i], by, bw[i], BTN_H}, bg);
        DrawRectangleLinesEx((Rectangle){bx[i], by, bw[i], BTN_H}, 1.0f,
                              (Color){150,150,165,255});
        Color txt_col = selected ? WHITE : (Color){40,40,60,255};
        label_font_draw_px(text, (int)(bx[i] + BTN_PAD_X), (int)(mid_y - 6), 12, txt_col);
        if (s) free(s);
    }
}

static int discrete_row_hit(const ManipCtrl* c, Vector2 m, float row_y) {
    float by = row_y + (ROW_H - BTN_H) * 0.5f;
    if (m.y < by || m.y > by + BTN_H) return -1;
    float bx[MAX_DISCRETE_VALUES], bw[MAX_DISCRETE_VALUES];
    int n = discrete_row_layout(c, bx, bw);
    for (int i = 0; i < n; i++)
        if (m.x >= bx[i] && m.x <= bx[i] + bw[i]) return i;
    return -1;
}

/* ------------------------------------------------------------------ */
/* Main loop                                                            */
/* ------------------------------------------------------------------ */

static void graphics_manipulate(const Expr* body, ManipCtrl* ctrls, int n_ctrls) {
    const char* no_window = getenv("MATHILDA_NO_GRAPHICS_WINDOW");
    if (no_window && no_window[0] != '\0') {
        printf("Manipulate: suppressed (MATHILDA_NO_GRAPHICS_WINDOW set).\n");
        return;
    }

    int win_w = 800, win_h = 500;
    float ctrl_h    = (float)n_ctrls * ROW_H + FOOTER_H;
    float content_h = (float)win_h - ctrl_h;
    float ctrl_y    = content_h;
    float footer_y  = ctrl_y + (float)n_ctrls * ROW_H;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(win_w, win_h, "Mathilda - Manipulate");
    SetTargetFPS(60);
    label_font_load();

    int   drag_ctrl  = -1;
    Expr* frame_expr = NULL;
    Graphics3DEmbedState* cam3d = graphics3d_embed_state_new();

    float reset_x = PAD;
    float reset_y = footer_y + (FOOTER_H - RESET_BTN_H) * 0.5f;

    while (!WindowShouldClose()) {
        /* ---- Evaluate body at current control values ---- */
        Rule* saved[MAX_CTRLS];
        for (int k = 0; k < n_ctrls; k++) {
            Expr* val = (ctrls[k].kind == MCTRL_RANGE)
                      ? expr_new_real(ctrls[k].value)
                      : expr_copy(ctrls[k].values[ctrls[k].selected_idx]);
            saved[k] = iter_spec_shadow(ctrls[k].var);
            Expr* vpat = expr_new_symbol(ctrls[k].var->data.symbol.name);
            symtab_add_own_value(ctrls[k].var->data.symbol.name, vpat, val);
            expr_free(vpat);
        }
        Expr* result = evaluate(expr_copy((Expr*)body));
        for (int k = n_ctrls - 1; k >= 0; k--)
            iter_spec_restore(ctrls[k].var, saved[k]);
        if (frame_expr) expr_free(frame_expr);
        frame_expr = result;

        /* ---- Handle input ---- */
        Vector2 mouse = GetMousePosition();

        if (IsKeyPressed(KEY_ESCAPE)) break;

        bool reset_hover = mouse.x >= reset_x && mouse.x <= reset_x + RESET_BTN_W
                         && mouse.y >= reset_y && mouse.y <= reset_y + RESET_BTN_H;

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            drag_ctrl = -1;
            bool handled = false;
            for (int k = 0; k < n_ctrls && !handled; k++) {
                float row_y = ctrl_y + (float)k * ROW_H;
                if (ctrls[k].kind == MCTRL_RANGE) {
                    if (range_row_hit(mouse, row_y, win_w)) {
                        drag_ctrl = k;
                        double frac = range_frac_from_mouse(mouse.x, win_w);
                        double v = ctrls[k].vmin + frac * (ctrls[k].vmax - ctrls[k].vmin);
                        ctrls[k].value = range_snap(&ctrls[k], v);
                        handled = true;
                    }
                } else {
                    int hit = discrete_row_hit(&ctrls[k], mouse, row_y);
                    if (hit >= 0) {
                        ctrls[k].selected_idx = (size_t)hit;
                        handled = true;
                    }
                }
            }
            if (!handled && reset_hover) {
                for (int k = 0; k < n_ctrls; k++) {
                    if (ctrls[k].kind == MCTRL_RANGE) ctrls[k].value = ctrls[k].default_value;
                    else                              ctrls[k].selected_idx = ctrls[k].default_idx;
                }
                graphics3d_embed_state_reset_view(cam3d);
            }
        }

        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) drag_ctrl = -1;

        if (drag_ctrl >= 0 && drag_ctrl < n_ctrls && ctrls[drag_ctrl].kind == MCTRL_RANGE) {
            double frac = range_frac_from_mouse(mouse.x, win_w);
            double v = ctrls[drag_ctrl].vmin + frac * (ctrls[drag_ctrl].vmax - ctrls[drag_ctrl].vmin);
            ctrls[drag_ctrl].value = range_snap(&ctrls[drag_ctrl], v);
        }

        /* ---- Render ---- */
        BeginDrawing();
        ClearBackground((Color){240,240,245,255});

        bool is_3d = frame_expr && is_graphics3d_expr(frame_expr);
        if (frame_expr && is_graphics2d_expr(frame_expr)) {
            graphics_render_in_region(frame_expr, 0.0f, 0.0f, (float)win_w, content_h);
        } else if (is_3d) {
            graphics3d_render_in_region(frame_expr, 0.0f, 0.0f, (float)win_w, content_h, cam3d);
        } else if (frame_expr) {
            char* s = expr_to_string(frame_expr);
            if (s) { label_font_draw_px(s, 20, 20, 16, DARKGRAY); free(s); }
        }

        DrawRectangle(0, (int)content_h, win_w, 1, (Color){170,170,185,255});

        for (int k = 0; k < n_ctrls; k++) {
            float row_y = ctrl_y + (float)k * ROW_H;
            if (ctrls[k].kind == MCTRL_RANGE)
                draw_range_row(&ctrls[k], win_w, row_y, mouse, drag_ctrl == k);
            else
                draw_discrete_row(&ctrls[k], win_w, row_y, mouse);
        }

        DrawRectangle(0, (int)footer_y, win_w, (int)FOOTER_H, (Color){225,225,230,255});
        Color reset_bg = reset_hover ? (Color){195,195,195,255} : (Color){215,215,215,255};
        DrawRectangleRec((Rectangle){reset_x, reset_y, RESET_BTN_W, RESET_BTN_H}, reset_bg);
        DrawRectangleLinesEx((Rectangle){reset_x, reset_y, RESET_BTN_W, RESET_BTN_H}, 1.0f,
                              (Color){145,145,145,255});
        label_font_draw_px("Reset", (int)(reset_x + 6), (int)(reset_y + 4), 11, BLACK);
        const char* hint = is_3d ? "drag: rotate  scroll: zoom  right-drag: pan   Esc: close"
                                  : "Esc: close";
        label_font_draw_px(hint, (int)(reset_x + RESET_BTN_W + 12), (int)(reset_y + 4), 11,
                 (Color){130,130,140,255});

        EndDrawing();
    }

    if (frame_expr) { expr_free(frame_expr); frame_expr = NULL; }
    graphics3d_embed_state_free(cam3d);
    label_font_unload();
    CloseWindow();
}

#else /* !USE_GRAPHICS */

static void graphics_manipulate(const Expr* body, ManipCtrl* ctrls, int n_ctrls) {
    (void)body; (void)ctrls; (void)n_ctrls;
    printf("Manipulate: not rendered -- graphics support not compiled in.\n");
}

#endif /* USE_GRAPHICS */

/* ------------------------------------------------------------------ */
/* builtin_manipulate                                                   */
/* ------------------------------------------------------------------ */

Expr* builtin_manipulate(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    const Expr* body = res->data.function.args[0];  /* held */

    /* Collect all positional control specs, evaluating each candidate arg
     * first (Manipulate is HoldAll, so nothing has been evaluated yet).
     * Stop at the first one that doesn't parse as a control spec. */
    ManipCtrl ctrls[MAX_CTRLS];
    int n_ctrls = 0;

    for (size_t i = 1; i < argc && n_ctrls < MAX_CTRLS; i++) {
        Expr* arg_ev = evaluate(expr_copy(res->data.function.args[i]));
        ManipCtrl c;
        memset(&c, 0, sizeof(c));
        bool ok = parse_control_spec(arg_ev, &c);
        expr_free(arg_ev);
        if (!ok) break;
        ctrls[n_ctrls++] = c;
    }

    if (n_ctrls == 0) return NULL;

    graphics_manipulate(body, ctrls, n_ctrls);

    for (int i = 0; i < n_ctrls; i++) manip_ctrl_free(&ctrls[i]);

    return expr_new_symbol("Null");
}
