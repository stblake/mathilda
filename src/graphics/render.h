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

#endif /* MATHILDA_GRAPHICS_RENDER_H */
