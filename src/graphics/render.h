/* render.h — Raylib backend for the graphics engine.
 *
 * render.c implements graphics_show() (declared in show.h) when
 * USE_GRAPHICS is compiled in. The only other entry point is the pure
 * window-sizing policy below, factored out so it can be unit-tested
 * without a live display. */
#ifndef MATHILDA_GRAPHICS_RENDER_H
#define MATHILDA_GRAPHICS_RENDER_H

#include <stdbool.h>

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
