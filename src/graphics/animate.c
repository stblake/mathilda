/* animate.c — Animate[expr, {t, tmin, tmax}, opts...]
 *
 * Opens a Raylib window (when USE_GRAPHICS is compiled in) with a
 * Manipulate-style control bar:
 *
 *   ┌─────────────────────────────────────────────────────────────┐
 *   │ a  │ 0 ─────────────●────────── 5  │  2.31               │
 *   ├─────────────────────────────────────────────────────────────┤
 *   │ [R][|<][▶][>|][dir]                              1x [-][+] │
 *   │ Space:play/pause  ←/→:step  R:reset  Esc:close            │
 *   └─────────────────────────────────────────────────────────────┘
 *
 * Each positional iterator argument {var, min, max} adds a labeled
 * slider row.  All iterators share a common animation phase so they
 * advance together.  Additional iterator arguments beyond the first
 * work as simultaneously animated variables (not static sliders).
 *
 * Options supported:
 *   AnimationDirection    Forward | Backward | ForwardBackward |
 *                         BackwardForward (default Forward)
 *   AnimationRate         parameter units per second of the FIRST iterator
 *                         (overrides DefaultDuration)
 *   AnimationRepetitions  integer or Infinity (default Infinity)
 *   AnimationRunning      True | False (default True)
 *   AppearanceElements    All | None | {"PlayPauseButton", ...}
 *   DefaultDuration       seconds for one full pass (default 1.0)
 *   ControlPlacement      Bottom | Top (default Bottom)
 *   RefreshRate           target display FPS (default 60)
 */

#include "animate.h"
#include "show.h"
#include "render.h"
#include "sym_names.h"
#include "symtab.h"
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
/* Per-iterator state                                                   */
/* ------------------------------------------------------------------ */

#define MAX_ITERS 8

typedef struct {
    const char* var_sym;
    double      tmin, tmax;
    double      t;          /* current value (= tmin + phase*(tmax-tmin)) */
} AnimIter;

/* ------------------------------------------------------------------ */
/* Option structs                                                       */
/* ------------------------------------------------------------------ */

typedef enum {
    ANIM_DIR_FORWARD = 0,
    ANIM_DIR_BACKWARD,
    ANIM_DIR_FORWARD_BACKWARD,
    ANIM_DIR_BACKWARD_FORWARD
} AnimDirection;

typedef enum { CTRL_BOTTOM = 0, CTRL_TOP } CtrlPlacement;

typedef struct {
    AnimDirection  direction;
    double         rate;           /* first-iter units/sec; 0 = use duration */
    int            repetitions;    /* -1 = Infinity */
    bool           running;
    double         duration;       /* seconds for one full pass (default 1.0) */
    CtrlPlacement  placement;
    double         refresh_rate;
    /* AppearanceElements flags */
    bool           show_play_pause;
    bool           show_step_buttons;
    bool           show_direction_button;
    bool           show_speed_controls;
    bool           show_reset_button;
} AnimateOpts;

static void animate_opts_defaults(AnimateOpts* o) {
    o->direction             = ANIM_DIR_FORWARD;
    o->rate                  = 0.0;
    o->repetitions           = -1;
    o->running               = true;
    o->duration              = 1.0;
    o->placement             = CTRL_BOTTOM;
    o->refresh_rate          = 60.0;
    o->show_play_pause       = true;
    o->show_step_buttons     = true;
    o->show_direction_button = true;
    o->show_speed_controls   = true;
    o->show_reset_button     = true;
}

static bool is_rule_arg(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_Rule
            || e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}

static bool is_sym(const Expr* e, const char* sym) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == sym;
}

static void parse_animate_opts(Expr** args, size_t start, size_t argc,
                                AnimateOpts* o) {
    for (size_t i = start; i < argc; i++) {
        Expr* arg = args[i];
        if (!is_rule_arg(arg)) continue;
        const Expr* lhs = arg->data.function.args[0];
        const Expr* rhs = arg->data.function.args[1];
        if (!lhs || lhs->type != EXPR_SYMBOL) continue;
        const char* key = lhs->data.symbol.name;

        if (key == SYM_AnimationDirection) {
            if      (is_sym(rhs, SYM_Forward))         o->direction = ANIM_DIR_FORWARD;
            else if (is_sym(rhs, SYM_Backward))        o->direction = ANIM_DIR_BACKWARD;
            else if (is_sym(rhs, SYM_ForwardBackward)) o->direction = ANIM_DIR_FORWARD_BACKWARD;
            else if (is_sym(rhs, SYM_BackwardForward)) o->direction = ANIM_DIR_BACKWARD_FORWARD;
        } else if (key == SYM_AnimationRate) {
            Expr* ev = evaluate(expr_copy((Expr*)rhs));
            if (ev->type == EXPR_REAL)         o->rate = ev->data.real;
            else if (ev->type == EXPR_INTEGER) o->rate = (double)ev->data.integer;
            expr_free(ev);
        } else if (key == SYM_AnimationRepetitions) {
            if (is_sym(rhs, SYM_Infinity)) {
                o->repetitions = -1;
            } else {
                Expr* ev = evaluate(expr_copy((Expr*)rhs));
                if (ev->type == EXPR_INTEGER) o->repetitions = (int)ev->data.integer;
                expr_free(ev);
            }
        } else if (key == SYM_AnimationRunning) {
            if      (is_sym(rhs, SYM_True))  o->running = true;
            else if (is_sym(rhs, SYM_False)) o->running = false;
        } else if (key == SYM_DefaultDuration) {
            Expr* ev = evaluate(expr_copy((Expr*)rhs));
            if (ev->type == EXPR_REAL)         o->duration = ev->data.real;
            else if (ev->type == EXPR_INTEGER) o->duration = (double)ev->data.integer;
            expr_free(ev);
        } else if (key == SYM_ControlPlacement) {
            if      (is_sym(rhs, SYM_Bottom)) o->placement = CTRL_BOTTOM;
            else if (is_sym(rhs, SYM_Top))    o->placement = CTRL_TOP;
        } else if (key == SYM_RefreshRate) {
            Expr* ev = evaluate(expr_copy((Expr*)rhs));
            if (ev->type == EXPR_REAL)         o->refresh_rate = ev->data.real;
            else if (ev->type == EXPR_INTEGER) o->refresh_rate = (double)ev->data.integer;
            expr_free(ev);
        } else if (key == SYM_AppearanceElements) {
            if (is_sym(rhs, SYM_None)) {
                o->show_play_pause       = false;
                o->show_step_buttons     = false;
                o->show_direction_button = false;
                o->show_speed_controls   = false;
                o->show_reset_button     = false;
            } else if (rhs->type == EXPR_FUNCTION
                       && rhs->data.function.head->type == EXPR_SYMBOL
                       && rhs->data.function.head->data.symbol.name == SYM_List) {
                o->show_play_pause       = false;
                o->show_step_buttons     = false;
                o->show_direction_button = false;
                o->show_speed_controls   = false;
                o->show_reset_button     = false;
                for (size_t j = 0; j < rhs->data.function.arg_count; j++) {
                    const Expr* el = rhs->data.function.args[j];
                    if (!el || el->type != EXPR_STRING) continue;
                    const char* s = el->data.string;
                    if (strcmp(s, "PlayPauseButton") == 0)     o->show_play_pause = true;
                    if (strcmp(s, "StepLeftButton") == 0 ||
                        strcmp(s, "StepRightButton") == 0)     o->show_step_buttons = true;
                    if (strcmp(s, "DirectionButton") == 0)     o->show_direction_button = true;
                    if (strcmp(s, "FasterSlowerButtons") == 0) o->show_speed_controls = true;
                    if (strcmp(s, "ResetButton") == 0)         o->show_reset_button = true;
                }
            }
            /* All: defaults already set */
        }
    }
}

/* ------------------------------------------------------------------ */
/* Numericize an iterator bound.                                        */
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

/* Try to parse {var, tmin, tmax} from a List expression.
 * Returns true on success; `var_sym`, `tmin`, `tmax` are set. */
static bool parse_iter_spec(const Expr* e, const char** var_sym,
                             double* tmin, double* tmax) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (e->data.function.head->data.symbol.name != SYM_List) return false;
    if (e->data.function.arg_count != 3) return false;
    const Expr* ve = e->data.function.args[0];
    if (!ve || ve->type != EXPR_SYMBOL) return false;
    *var_sym = ve->data.symbol.name;
    return numericize(e->data.function.args[1], tmin)
        && numericize(e->data.function.args[2], tmax);
}

/* ------------------------------------------------------------------ */
/* Rendering helpers (Raylib)                                           */
/* ------------------------------------------------------------------ */

#ifdef USE_GRAPHICS
#include <raylib.h>

#define SLIDER_ROW_H   28.0f   /* height of each labeled-variable slider row */
#define BTNS_ROW_H     34.0f   /* height of the playback button row          */
#define HELP_ROW_H     14.0f   /* help text at the very bottom               */
#define LABEL_W        48.0f   /* pixels reserved for the var-name label     */
#define VALUE_W        54.0f   /* pixels reserved for the current-value text */
#define SLIDER_PAD     10.0f   /* horizontal gap around the track            */
#define BTN_SIZE       26.0f
#define BTN_GAP         7.0f

static float ctrl_total_h(int n_iters) {
    return (float)n_iters * SLIDER_ROW_H + BTNS_ROW_H + HELP_ROW_H;
}

/* ---- Icon drawing ---- */

static void tri(Vector2 a, Vector2 b, Vector2 c, Color col) {
    DrawTriangle(a, b, c, col);
}

static void draw_play_pause_icon(float x, float y, float sz, bool playing) {
    float cx = x + sz * 0.5f, cy = y + sz * 0.5f, m = sz * 0.27f;
    if (playing) {
        DrawRectangle((int)(cx - m - 2), (int)(cy - m), 4, (int)(m * 2), BLACK);
        DrawRectangle((int)(cx + 1),     (int)(cy - m), 4, (int)(m * 2), BLACK);
    } else {
        tri((Vector2){ cx - m * 0.8f, cy - m },
            (Vector2){ cx - m * 0.8f, cy + m },
            (Vector2){ cx + m,         cy     }, BLACK);
    }
}

static void draw_step_left_icon(float x, float y, float sz) {
    float cx = x + sz * 0.5f, cy = y + sz * 0.5f, m = sz * 0.24f;
    DrawRectangle((int)(cx - m - 3), (int)(cy - m), 3, (int)(m * 2), BLACK);
    tri((Vector2){ cx + m * 0.6f, cy - m },
        (Vector2){ cx - m * 0.4f, cy     },
        (Vector2){ cx + m * 0.6f, cy + m }, BLACK);
}

static void draw_step_right_icon(float x, float y, float sz) {
    float cx = x + sz * 0.5f, cy = y + sz * 0.5f, m = sz * 0.24f;
    /* CW winding (screen y-down): bottom-left, right-center, top-left */
    tri((Vector2){ cx - m * 0.6f, cy + m },
        (Vector2){ cx + m * 0.4f, cy     },
        (Vector2){ cx - m * 0.6f, cy - m }, BLACK);
    DrawRectangle((int)(cx + m * 0.4f), (int)(cy - m), 3, (int)(m * 2), BLACK);
}

/* Direction button: text label is clearer than small triangles. */
static void draw_direction_label(float x, float y, float sz, AnimDirection dir) {
    const char* lbl;
    Color col;
    switch (dir) {
        case ANIM_DIR_FORWARD:          lbl = "FWD"; col = (Color){40,80,200,255};  break;
        case ANIM_DIR_BACKWARD:         lbl = "BWD"; col = (Color){190,40,40,255};  break;
        case ANIM_DIR_FORWARD_BACKWARD: lbl = "F<>"; col = (Color){120,40,160,255}; break;
        case ANIM_DIR_BACKWARD_FORWARD: lbl = "<>F"; col = (Color){120,40,160,255}; break;
        default:                         lbl = "?";   col = DARKGRAY;                break;
    }
    int tw = MeasureText(lbl, 10);
    DrawText(lbl, (int)(x + sz * 0.5f - (float)tw * 0.5f),
                  (int)(y + sz * 0.5f - 5), 10, col);
}

static void draw_btn_bg(float x, float y, float sz, bool hovered) {
    Color bg = hovered ? (Color){195,195,195,255} : (Color){215,215,215,255};
    DrawRectangleRec((Rectangle){x, y, sz, sz}, bg);
    DrawRectangleLinesEx((Rectangle){x, y, sz, sz}, 1.0f, (Color){145,145,145,255});
}

static bool btn_hit(Vector2 m, float bx, float by, float sz) {
    return m.x >= bx && m.x <= bx + sz && m.y >= by && m.y <= by + sz;
}

/* ---- Per-iterator labeled slider row ----
 *
 * Each row occupies [row_y, row_y + SLIDER_ROW_H].
 *
 *  | LABEL_W | SLIDER_PAD | track ... | SLIDER_PAD | VALUE_W |
 *
 * Returns true if the user is clicking/dragging inside this row.
 * Sets *drag_frac to [0,1] when dragging; caller updates iter->t. */
static void draw_iter_slider(const AnimIter* it, int win_w,
                              float row_y, Vector2 mouse,
                              bool dragging) {
    float track_x0 = LABEL_W + SLIDER_PAD;
    float track_x1 = (float)win_w - VALUE_W - SLIDER_PAD;
    float track_len = track_x1 - track_x0;
    float mid_y     = row_y + SLIDER_ROW_H * 0.5f;

    /* Background stripe */
    DrawRectangle(0, (int)row_y, win_w, (int)SLIDER_ROW_H, (Color){235,235,240,255});
    DrawRectangle(0, (int)(row_y + SLIDER_ROW_H - 1), win_w, 1, (Color){200,200,210,255});

    /* Variable name label */
    DrawText(it->var_sym, 6, (int)(mid_y - 7), 14, (Color){60,60,100,255});

    /* Min / max hints */
    char lo[24], hi[24];
    snprintf(lo, sizeof(lo), "%.4g", it->tmin);
    snprintf(hi, sizeof(hi), "%.4g", it->tmax);
    int lo_w = MeasureText(lo, 10);
    DrawText(lo, (int)(track_x0 - lo_w - 3), (int)(mid_y - 5), 10, DARKGRAY);
    DrawText(hi, (int)(track_x1 + 3),         (int)(mid_y - 5), 10, DARKGRAY);

    /* Track */
    DrawRectangle((int)track_x0, (int)(mid_y - 2), (int)track_len + 1, 4,
                  (Color){185,185,195,255});

    double span = it->tmax - it->tmin;
    double frac = (span > 0) ? (it->t - it->tmin) / span : 0.0;
    if (frac < 0.0) frac = 0.0;
    if (frac > 1.0) frac = 1.0;

    /* Filled portion */
    float fill = (float)(frac * track_len);
    if (fill > 0.0f)
        DrawRectangle((int)track_x0, (int)(mid_y - 2), (int)fill + 1, 4,
                      (Color){60,110,215,255});

    /* Handle */
    float hx = track_x0 + fill;
    bool hover_handle = fabsf(mouse.x - hx) < 12.0f
                     && mouse.y >= row_y && mouse.y < row_y + SLIDER_ROW_H;
    Color handle_col = (hover_handle || dragging)
                       ? (Color){30,80,200,255}
                       : (Color){50,100,220,255};
    DrawCircle((int)hx, (int)mid_y, 7.5f, handle_col);
    DrawCircleLines((int)hx, (int)mid_y, 7.5f, (Color){20,60,160,255});

    /* Current value */
    char val[32];
    snprintf(val, sizeof(val), "%.5g", it->t);
    int val_w = MeasureText(val, 12);
    DrawText(val, (int)((float)win_w - VALUE_W / 2.0f - (float)val_w / 2.0f),
             (int)(mid_y - 6), 12, (Color){40,40,80,255});
}

/* Check whether a click/drag in (mx, my) targets the slider for row `i`. */
static bool slider_row_hit(Vector2 m, int i, int win_w) {
    float row_y   = (float)i * SLIDER_ROW_H;
    float track_x0 = LABEL_W + SLIDER_PAD;
    float track_x1 = (float)win_w - VALUE_W - SLIDER_PAD;
    float mid_y    = row_y + SLIDER_ROW_H * 0.5f;
    return m.x >= track_x0 - 14.0f && m.x <= track_x1 + 14.0f
        && m.y >= mid_y - 12.0f && m.y <= mid_y + 12.0f;
}

/* Compute the fractional slider position from mouse x for row i. */
static double slider_frac_from_mouse(float mx, int win_w) {
    float track_x0 = LABEL_W + SLIDER_PAD;
    float track_x1 = (float)win_w - VALUE_W - SLIDER_PAD;
    float frac = (mx - track_x0) / (track_x1 - track_x0);
    if (frac < 0.0f) frac = 0.0f;
    if (frac > 1.0f) frac = 1.0f;
    return (double)frac;
}

/* ------------------------------------------------------------------ */
/* Main animation loop                                                  */
/* ------------------------------------------------------------------ */

static bool is_graphics_expr(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_Graphics
            || e->data.function.head->data.symbol.name == SYM_Graphics3D)
        && e->data.function.arg_count >= 1;
}

static void graphics_animate(const Expr* body,
                               AnimIter* iters, int n_iters,
                               const AnimateOpts* in_opts) {
    const char* no_window = getenv("MATHILDA_NO_GRAPHICS_WINDOW");
    if (no_window && no_window[0] != '\0') {
        printf("Animate: suppressed (MATHILDA_NO_GRAPHICS_WINDOW set).\n");
        return;
    }

    AnimateOpts opts = *in_opts;

    /* Base rate from first iterator's span. */
    double span0   = iters[0].tmax - iters[0].tmin;
    if (span0 <= 0.0) span0 = 1.0;
    double base_rate = opts.rate;
    if (base_rate <= 0.0 && opts.duration > 0.0)
        base_rate = span0 / opts.duration;
    if (base_rate <= 0.0) base_rate = span0;

    double speed_mult = 1.0;

    int win_w = 800, win_h = 500;
    float ctrl_h    = ctrl_total_h(n_iters);
    /* Place slider rows first, then buttons, regardless of ControlPlacement. */
    float sliders_y = (opts.placement == CTRL_TOP) ? 0.0f : ((float)win_h - ctrl_h);
    float btns_y    = sliders_y + (float)n_iters * SLIDER_ROW_H;
    float content_y = (opts.placement == CTRL_TOP) ? ctrl_h : 0.0f;
    float content_h = (float)win_h - ctrl_h;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(win_w, win_h, "Mathilda - Animate");
    SetTargetFPS((int)opts.refresh_rate);

    /* Animation phase in [0,1]; drives all iterators simultaneously. */
    double phi      = 0.0;      /* 0 = tmin end, 1 = tmax end */
    bool going_fwd  = (opts.direction != ANIM_DIR_BACKWARD &&
                       opts.direction != ANIM_DIR_BACKWARD_FORWARD);
    if (!going_fwd) phi = 1.0;
    bool running    = opts.running;
    int  reps_done  = 0;

    /* Per-slider drag state. */
    int  drag_slider = -1;   /* which slider is being dragged (-1 = none) */

    Expr* frame_expr = NULL;

    /* Button positions (recomputed each frame but constant during a frame). */
    float btn_y    = 0.0f;
    int   n_slots  = 0;
    float reset_x = 0.0f, step_l_x = 0.0f, pp_x = 0.0f;
    float step_r_x = 0.0f, dir_x = 0.0f;
    float slower_x = 0.0f, faster_x = 0.0f;

    /* Sync iterators from current phase. */
    for (int k = 0; k < n_iters; k++) {
        double span = iters[k].tmax - iters[k].tmin;
        iters[k].t = iters[k].tmin + phi * span;
    }

    while (!WindowShouldClose()) {
        double dt = (double)GetFrameTime();

        /* ---- Update phase ---- */
        if (running) {
            /* Convert base_rate (first-iter units/sec) to phase/sec. */
            double phase_rate = base_rate * speed_mult / span0;
            if (going_fwd) phi += phase_rate * dt;
            else           phi -= phase_rate * dt;

            if (going_fwd && phi >= 1.0) {
                switch (opts.direction) {
                    case ANIM_DIR_FORWARD:
                        if (opts.repetitions < 0 || reps_done + 1 < opts.repetitions) {
                            phi = 0.0; reps_done++;
                        } else { phi = 1.0; running = false; }
                        break;
                    case ANIM_DIR_FORWARD_BACKWARD:
                        phi = 1.0; going_fwd = false; reps_done++;
                        if (opts.repetitions >= 0 && reps_done >= opts.repetitions)
                            running = false;
                        break;
                    case ANIM_DIR_BACKWARD_FORWARD:
                        phi = 1.0; going_fwd = false; break;
                    default:
                        phi = 1.0; running = false; break;
                }
            } else if (!going_fwd && phi <= 0.0) {
                switch (opts.direction) {
                    case ANIM_DIR_BACKWARD:
                        if (opts.repetitions < 0 || reps_done + 1 < opts.repetitions) {
                            phi = 1.0; reps_done++;
                        } else { phi = 0.0; running = false; }
                        break;
                    case ANIM_DIR_FORWARD_BACKWARD:
                        phi = 0.0; going_fwd = true; break;
                    case ANIM_DIR_BACKWARD_FORWARD:
                        phi = 0.0; going_fwd = true; reps_done++;
                        if (opts.repetitions >= 0 && reps_done >= opts.repetitions)
                            running = false;
                        break;
                    default:
                        phi = 0.0; running = false; break;
                }
            }
            if (phi < 0.0) phi = 0.0;
            if (phi > 1.0) phi = 1.0;
            /* Sync iterators. */
            for (int k = 0; k < n_iters; k++) {
                double span = iters[k].tmax - iters[k].tmin;
                iters[k].t = iters[k].tmin + phi * span;
            }
        }

        /* ---- Evaluate body at current iterator values ---- */
        for (int k = 0; k < n_iters; k++) {
            Expr* tval = expr_new_real(iters[k].t);
            Expr* vpat = expr_new_symbol(iters[k].var_sym);
            symtab_add_own_value(iters[k].var_sym, vpat, expr_copy(tval));
            expr_free(vpat);
            expr_free(tval);
        }
        Expr* result = evaluate(expr_copy((Expr*)body));
        for (int k = 0; k < n_iters; k++)
            symtab_clear_symbol(iters[k].var_sym);
        if (frame_expr) expr_free(frame_expr);
        frame_expr = result;

        /* ---- Button layout ---- */
        btn_y = btns_y + (BTNS_ROW_H - BTN_SIZE) * 0.5f;
        #define SLOT(bx) bx = BTN_GAP + (float)(n_slots++) * (BTN_SIZE + BTN_GAP)
        reset_x = step_l_x = pp_x = step_r_x = dir_x = 0.0f;
        slower_x = faster_x = 0.0f;
        n_slots = 0;
        if (opts.show_reset_button)     { SLOT(reset_x); }
        if (opts.show_step_buttons)     { SLOT(step_l_x); }
        if (opts.show_play_pause)       { SLOT(pp_x); }
        if (opts.show_step_buttons)     { SLOT(step_r_x); }
        if (opts.show_direction_button) { SLOT(dir_x); }
        #undef SLOT

        /* Speed controls on right side. */
        float right_edge = (float)win_w - BTN_GAP;
        if (opts.show_speed_controls) {
            faster_x = right_edge - BTN_SIZE;
            slower_x = faster_x - BTN_GAP - BTN_SIZE;
            right_edge = slower_x - BTN_GAP - 42.0f; /* speed-mult label */
        }
        (void)right_edge;

        /* ---- Handle input ---- */
        Vector2 mouse = GetMousePosition();

        if (IsKeyPressed(KEY_SPACE)) running = !running;
        if (IsKeyPressed(KEY_ESCAPE)) break;
        if (IsKeyPressed(KEY_R)) {
            phi = (opts.direction == ANIM_DIR_BACKWARD ||
                   opts.direction == ANIM_DIR_BACKWARD_FORWARD) ? 1.0 : 0.0;
            going_fwd = (opts.direction != ANIM_DIR_BACKWARD &&
                         opts.direction != ANIM_DIR_BACKWARD_FORWARD);
            reps_done = 0;
            for (int k = 0; k < n_iters; k++) {
                double span = iters[k].tmax - iters[k].tmin;
                iters[k].t  = iters[k].tmin + phi * span;
            }
        }
        if (IsKeyPressed(KEY_LEFT)) {
            running = false;
            phi -= 0.02;
            if (phi < 0.0) phi = 0.0;
            for (int k = 0; k < n_iters; k++) {
                double span = iters[k].tmax - iters[k].tmin;
                iters[k].t  = iters[k].tmin + phi * span;
            }
        }
        if (IsKeyPressed(KEY_RIGHT)) {
            running = false;
            phi += 0.02;
            if (phi > 1.0) phi = 1.0;
            for (int k = 0; k < n_iters; k++) {
                double span = iters[k].tmax - iters[k].tmin;
                iters[k].t  = iters[k].tmin + phi * span;
            }
        }

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            /* Check slider rows first. */
            drag_slider = -1;
            for (int k = 0; k < n_iters; k++) {
                /* slider_row_hit expects y relative to slider area origin */
                if (slider_row_hit((Vector2){mouse.x, mouse.y - sliders_y}, k, win_w)) {
                    drag_slider = k;
                    running = false;
                    double frac = slider_frac_from_mouse(mouse.x, win_w);
                    double span  = iters[k].tmax - iters[k].tmin;
                    iters[k].t   = iters[k].tmin + frac * span;
                    /* Independent: only update this slider, not the others. */
                    break;
                }
            }

            /* Button row. */
            if (drag_slider < 0) {
                if (opts.show_reset_button && btn_hit(mouse, reset_x, btn_y, BTN_SIZE)) {
                    phi = (opts.direction == ANIM_DIR_BACKWARD ||
                           opts.direction == ANIM_DIR_BACKWARD_FORWARD) ? 1.0 : 0.0;
                    going_fwd = (opts.direction != ANIM_DIR_BACKWARD &&
                                 opts.direction != ANIM_DIR_BACKWARD_FORWARD);
                    reps_done = 0;
                    for (int k = 0; k < n_iters; k++) {
                        double span = iters[k].tmax - iters[k].tmin;
                        iters[k].t  = iters[k].tmin + phi * span;
                    }
                }
                if (opts.show_step_buttons && btn_hit(mouse, step_l_x, btn_y, BTN_SIZE)) {
                    running = false;
                    phi -= 0.05; if (phi < 0.0) phi = 0.0;
                    for (int k = 0; k < n_iters; k++) {
                        double span = iters[k].tmax - iters[k].tmin;
                        iters[k].t  = iters[k].tmin + phi * span;
                    }
                }
                if (opts.show_play_pause && btn_hit(mouse, pp_x, btn_y, BTN_SIZE))
                    running = !running;
                if (opts.show_step_buttons && btn_hit(mouse, step_r_x, btn_y, BTN_SIZE)) {
                    running = false;
                    phi += 0.05; if (phi > 1.0) phi = 1.0;
                    for (int k = 0; k < n_iters; k++) {
                        double span = iters[k].tmax - iters[k].tmin;
                        iters[k].t  = iters[k].tmin + phi * span;
                    }
                }
                if (opts.show_direction_button && btn_hit(mouse, dir_x, btn_y, BTN_SIZE)) {
                    switch (opts.direction) {
                        case ANIM_DIR_FORWARD:          opts.direction = ANIM_DIR_BACKWARD;         break;
                        case ANIM_DIR_BACKWARD:         opts.direction = ANIM_DIR_FORWARD_BACKWARD; break;
                        case ANIM_DIR_FORWARD_BACKWARD: opts.direction = ANIM_DIR_BACKWARD_FORWARD; break;
                        case ANIM_DIR_BACKWARD_FORWARD: opts.direction = ANIM_DIR_FORWARD;          break;
                    }
                    going_fwd = (opts.direction != ANIM_DIR_BACKWARD &&
                                 opts.direction != ANIM_DIR_BACKWARD_FORWARD);
                }
                if (opts.show_speed_controls) {
                    if (btn_hit(mouse, slower_x, btn_y, BTN_SIZE)) {
                        speed_mult /= 2.0;
                        if (speed_mult < 0.125) speed_mult = 0.125;
                    }
                    if (btn_hit(mouse, faster_x, btn_y, BTN_SIZE)) {
                        speed_mult *= 2.0;
                        if (speed_mult > 8.0) speed_mult = 8.0;
                    }
                }
            }
        }

        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) drag_slider = -1;

        /* Drag update — only the dragged slider moves. */
        if (drag_slider >= 0 && drag_slider < n_iters) {
            int k = drag_slider;
            double frac = slider_frac_from_mouse(mouse.x, win_w);
            double span  = iters[k].tmax - iters[k].tmin;
            iters[k].t   = iters[k].tmin + frac * span;
        }

        /* ---- Render ---- */
        BeginDrawing();
        ClearBackground((Color){240,240,245,255});

        /* Content area */
        if (frame_expr && is_graphics_expr(frame_expr)) {
            graphics_render_in_region(frame_expr, 0.0f, content_y,
                                       (float)win_w, content_h);
        } else if (frame_expr) {
            char* s = expr_to_string(frame_expr);
            if (s) { DrawText(s, 20, (int)(content_y + 20), 16, DARKGRAY); free(s); }
        }

        /* Separator line between content and control bar. */
        float sep_y = (opts.placement == CTRL_TOP)
                      ? ctrl_h - 1.0f
                      : (float)win_h - ctrl_h;
        DrawRectangle(0, (int)sep_y, win_w, 1, (Color){170,170,185,255});

        /* Slider rows */
        for (int k = 0; k < n_iters; k++) {
            float row_y = sliders_y + (float)k * SLIDER_ROW_H;
            draw_iter_slider(&iters[k], win_w, row_y, mouse,
                              drag_slider == k);
        }

        /* Button row background */
        DrawRectangle(0, (int)btns_y, win_w,
                      (int)(BTNS_ROW_H + HELP_ROW_H), (Color){225,225,230,255});

        /* Buttons */
        if (opts.show_reset_button) {
            draw_btn_bg(reset_x, btn_y, BTN_SIZE, btn_hit(mouse, reset_x, btn_y, BTN_SIZE));
            DrawText("R", (int)(reset_x + BTN_SIZE * 0.5f - 4),
                     (int)(btn_y + BTN_SIZE * 0.5f - 7), 14, BLACK);
        }
        if (opts.show_step_buttons) {
            draw_btn_bg(step_l_x, btn_y, BTN_SIZE, btn_hit(mouse, step_l_x, btn_y, BTN_SIZE));
            draw_step_left_icon(step_l_x, btn_y, BTN_SIZE);
            draw_btn_bg(step_r_x, btn_y, BTN_SIZE, btn_hit(mouse, step_r_x, btn_y, BTN_SIZE));
            draw_step_right_icon(step_r_x, btn_y, BTN_SIZE);
        }
        if (opts.show_play_pause) {
            draw_btn_bg(pp_x, btn_y, BTN_SIZE, btn_hit(mouse, pp_x, btn_y, BTN_SIZE));
            draw_play_pause_icon(pp_x, btn_y, BTN_SIZE, running);
        }
        if (opts.show_direction_button) {
            draw_btn_bg(dir_x, btn_y, BTN_SIZE, btn_hit(mouse, dir_x, btn_y, BTN_SIZE));
            draw_direction_label(dir_x, btn_y, BTN_SIZE, opts.direction);
        }
        if (opts.show_speed_controls) {
            char spd[16];
            snprintf(spd, sizeof(spd), "%.3gx", speed_mult);
            DrawText(spd, (int)(slower_x - 40.0f), (int)(btn_y + BTN_SIZE * 0.5f - 6),
                     11, DARKGRAY);
            draw_btn_bg(slower_x, btn_y, BTN_SIZE, btn_hit(mouse, slower_x, btn_y, BTN_SIZE));
            DrawText("-", (int)(slower_x + BTN_SIZE * 0.5f - 3),
                     (int)(btn_y + BTN_SIZE * 0.5f - 8), 16, BLACK);
            draw_btn_bg(faster_x, btn_y, BTN_SIZE, btn_hit(mouse, faster_x, btn_y, BTN_SIZE));
            DrawText("+", (int)(faster_x + BTN_SIZE * 0.5f - 4),
                     (int)(btn_y + BTN_SIZE * 0.5f - 8), 16, BLACK);
        }

        /* Help text */
        float help_y = btns_y + BTNS_ROW_H + 1.0f;
        DrawText("Space: play/pause   \xE2\x86\x90/\xE2\x86\x92: step 2%   R: reset   Esc: close",
                 4, (int)help_y, 10, (Color){130,130,140,255});

        EndDrawing();
    }

    if (frame_expr) { expr_free(frame_expr); frame_expr = NULL; }
    CloseWindow();
}

#else /* !USE_GRAPHICS */

static void graphics_animate(const Expr* body,
                               AnimIter* iters, int n_iters,
                               const AnimateOpts* opts) {
    (void)body; (void)iters; (void)n_iters; (void)opts;
    printf("Animate: not rendered -- graphics support not compiled in.\n");
}

#endif /* USE_GRAPHICS */

/* ------------------------------------------------------------------ */
/* builtin_animate                                                      */
/* ------------------------------------------------------------------ */


Expr* builtin_animate(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    const Expr* body = res->data.function.args[0];  /* held */

    /* Collect all positional iterator specs {var, min, max}
     * starting from arg 1. Stop at the first Rule[...] or non-list. */
    AnimIter iters[MAX_ITERS];
    int n_iters = 0;
    size_t opts_start = 1;

    for (size_t i = 1; i < argc && n_iters < MAX_ITERS; i++) {
        Expr* arg_ev = evaluate(expr_copy(res->data.function.args[i]));
        const char* vs = NULL;
        double lo = 0.0, hi = 1.0;
        bool ok = parse_iter_spec(arg_ev, &vs, &lo, &hi);
        expr_free(arg_ev);
        if (!ok) { opts_start = i; break; }
        if (lo > hi) { double tmp = lo; lo = hi; hi = tmp; }
        iters[n_iters].var_sym = vs;
        iters[n_iters].tmin   = lo;
        iters[n_iters].tmax   = hi;
        iters[n_iters].t      = lo;
        n_iters++;
        opts_start = i + 1;
    }

    if (n_iters == 0) return NULL;

    AnimateOpts opts;
    animate_opts_defaults(&opts);
    parse_animate_opts(res->data.function.args, opts_start, argc, &opts);

    graphics_animate(body, iters, n_iters, &opts);

    return expr_new_symbol("Null");
}
