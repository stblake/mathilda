/* render.h — Raylib backend for the graphics engine.
 *
 * render.c implements graphics_show() (declared in show.h) when
 * USE_GRAPHICS is compiled in. The only other entry point is the pure
 * window-sizing policy below, factored out so it can be unit-tested
 * without a live display. */
#ifndef MATHILDA_GRAPHICS_RENDER_H
#define MATHILDA_GRAPHICS_RENDER_H

#include <stdbool.h>
#include "expr.h"

/* Resolve the on-screen window height (px) from the AspectRatio setting.
 *
 * AspectRatio is the height-to-width ratio of the plot, so the window is
 * reshaped to width * ratio rather than letterboxing the curve inside a
 * fixed box. aspect_ratio > 0 is an explicit ratio; aspect_ratio <= 0 means
 * Automatic (use the data extent ratio data_h/data_w). The given height is
 * kept unchanged when height_pinned (ImageSize -> {w,h}) or aspect_full
 * (AspectRatio -> Full fills the ImageSize box). The result is clamped to a
 * screen-friendly [100, 2000] band so an extreme Automatic ratio letterboxes
 * gracefully instead of opening a giant window. */
long gfx_window_height(long width, long height, double aspect_ratio,
                       bool aspect_full, bool height_pinned,
                       double data_w, double data_h);

/* Frame/Axes margin policy: mL/mR depend only on window width, mT/mB only on
 * window height (partly floored in px), matching the layout render.c uses
 * for the actual plot region. Factored out so gfx_window_height_fit_region
 * below and the renderer's own region layout can never drift apart, and so
 * both are unit-testable headless. */
void gfx_horizontal_margins(float width, bool frame, bool axes,
                            bool frame_label, float* mL, float* mR);
void gfx_vertical_margins(float height, bool frame, bool axes,
                          bool frame_label, float* mT, float* mB);

/* Refines gfx_window_height's result against the margin-reduced DATA REGION
 * rather than the outer window. gfx_window_height alone sizes the window so
 * *it* matches aspect_ratio, but the Frame/Axes margins above are fixed
 * pixel amounts (partly floored) rather than proportional, so the region
 * left after subtracting them lands at a different ratio than the window --
 * an AspectRatio-driven Frame plot (ArrayPlot's rows/cols default,
 * DensityPlot's 1, ...) then letterboxes inside its own frame instead of
 * filling it edge-to-edge. Solves for the height whose region hits the
 * target ratio via fixed-point iteration (a contraction, since each margin
 * is either height-independent or a small fraction of height); a no-op
 * (returns height unchanged) in every case gfx_window_height itself is a
 * no-op, when neither frame nor axes is drawn, or when there isn't room for
 * both margins to fit. */
long gfx_window_height_fit_region(long width, long height, double aspect_ratio,
                                  bool aspect_full, bool height_pinned,
                                  double data_w, double data_h,
                                  bool frame, bool axes, bool frame_label);

/* Pixel dimensions for a rasterised 2D Graphics[...] export (PNG/JPEG). Sized
 * exactly the way graphics_show sizes its on-screen window: the width from
 * ImageSize (default 800), the height derived from AspectRatio via
 * gfx_window_height + gfx_window_height_fit_region so the data region fills the
 * frame in both directions. Passing a fixed 4:3 canvas instead letterboxed
 * every aspect-driven plot (ArrayPlot, DensityPlot, ContourPlot, ...) with wide
 * horizontal padding. Writes the chosen size to out_w and out_h (clamped to
 * [16, 8000]); falls back to 800x600 when g is not a Graphics[...]. */
void graphics_raster_dims(const Expr* g, int* out_w, int* out_h);

/* Convert a CMYK color to RGB using the standard (profile-free) subtractive
 * model: r = (1-c)(1-k), g = (1-m)(1-k), b = (1-y)(1-k). Each input is first
 * clipped to [0,1]; outputs land in [0,1]. Factored out of the renderer so the
 * conversion is unit-testable headless. */
void cmyk_to_rgb(double c, double m, double y, double k,
                 double* r, double* g, double* b);

/* Render a 2D Graphics[...] expression into a sub-region of an already-open
 * Raylib window. The region is [rx, ry] origin with size [rw x rh] in screen
 * pixels. Must be called between BeginDrawing() / EndDrawing(). The background
 * color of the region is painted from the graphics object's Background option.
 * No-op when USE_GRAPHICS is not compiled in (this file is excluded then). */
void graphics_render_in_region(const Expr* graphics_expr,
                                float rx, float ry, float rw, float rh);

/* Number of minor (sub-)tick intervals per major frame-tick interval, chosen
 * from the leading digit of `step` so minor ticks land on round values: a
 * "nice" step of 1 splits into 5, 2 into 4, 5 into 5. Factored out of the
 * frame renderer so the tick-subdivision policy can be unit-tested headless. */
int frame_minor_divs(double step);

/* Signed area (shoelace formula) of a closed polygon given as parallel x/y
 * arrays of length count. Positive in a y-down (screen-like) convention
 * means the vertices are wound clockwise as drawn; Polygon[]'s fill
 * reverses the vertex order whenever this is positive, since the
 * underlying triangle-fan fill requires counter-clockwise winding but
 * Polygon[] itself (like Mathematica's) imposes no winding convention on
 * the caller. Factored out so the winding-detection policy is
 * unit-testable headless. */
double polygon_signed_area(const double* xs, const double* ys, int count);

#endif /* MATHILDA_GRAPHICS_RENDER_H */
