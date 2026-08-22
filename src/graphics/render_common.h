/* render_common.h — small Raylib-backend helpers shared between the 2D
 * renderer (render.c) and the 3D renderer (render3d.c).
 *
 * Everything here is already dimensionless (a color, a coerced double, a
 * tick-spacing policy) -- none of it is specific to Camera2D vs Camera3D,
 * so render3d.c reuses these instead of re-implementing them. */
#ifndef MATHILDA_GRAPHICS_RENDER_COMMON_H
#define MATHILDA_GRAPHICS_RENDER_COMMON_H

#include "expr.h"
#include "primitives.h"
#include <raylib.h>
#include <stdbool.h>

/* Coerce a literal or symbolic-but-numeric (Pi/2, Sqrt[2], ...) Expr to a
 * double, the latter via N[] -- see render.c for the full rationale. */
bool expr_to_d(const Expr* e, double* out);

/* Resolves any recognized color-literal head (RGBColor, GrayLevel, Hue) to
 * an RGBA8; returns false (leaving *out untouched) for anything else. */
bool resolve_color(const Expr* e, RGBA8* out);

Color to_raylib(RGBA8 c);

/* "Nice" axis tick step (1/2/5 x 10^k) for a given data range and a target
 * tick count. */
double nice_step(double range, int target_ticks);

/* True if this process can open a GUI window (Aqua session on macOS, X11/
 * Wayland display on Linux). The offscreen export paths screen for this before
 * InitWindow, which aborts inside GLFW rather than failing when no monitor is
 * available. Defined in render.c; used by both renderers. */
bool gui_session_available(void);

/* Find the $StreamColorBar[lo, hi, cfn] metadata node embedded by the
 * plotter in the Graphics/Graphics3D option list.  Returns NULL if absent. */
const Expr* find_color_bar(const Expr* graphics_expr);

/* Draw a vertical color scale bar in screen space.  Usable from both the
 * 2D and 3D render paths (called after EndMode3D so drawing is 2D). */
void draw_color_bar(float bar_x, float bar_y, float bar_w, float bar_h,
                    double spd_min, double spd_max, const Expr* cfn);

#endif /* MATHILDA_GRAPHICS_RENDER_COMMON_H */
