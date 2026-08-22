/* graphics_export.h -- Export a Graphics[...] expression to an image file.
 *
 * Two backends, chosen by the requested format:
 *
 *   * graphics_export_pdf  -- a dependency-free vector PDF writer. It walks the
 *     Graphics primitive tree directly (no raylib), so it is compiled in every
 *     configuration and works headless. This is the high-quality path for
 *     print/publication and is what `Export["f.pdf", Plot[...]]` uses.
 *
 *   * graphics_export_raster -- PNG/JPEG via the existing Raylib renderer drawn
 *     offscreen into a RenderTexture, so the pixels match the on-screen plot
 *     exactly (text, axes, ticks, everything). Only meaningful when USE_GRAPHICS
 *     is compiled in AND a GL context is available; returns 0 otherwise (the
 *     caller then reports a clear error while PDF export still works).
 *
 * Both return 1 on success, 0 on failure. Only 2D Graphics[...] is handled;
 * Graphics3D returns 0 (a deliberate, documented limit -- see the docstring).
 */
#ifndef MATHILDA_GRAPHICS_EXPORT_H
#define MATHILDA_GRAPHICS_EXPORT_H

#include "expr.h"

/* Vector PDF. Dependency-free, headless. `path` is the output file. */
int graphics_export_pdf(const Expr* graphics_expr, const char* path);

/* Render a 2D Graphics[...] to a freshly malloc'd RGBA8 pixel buffer via the
 * Raylib renderer offscreen (w*h*4 bytes, top row first). Returns NULL on
 * failure (no GUI session / no GL context / Graphics3D). On success writes the
 * actual pixel dimensions to *out_w,*out_h; the caller owns the buffer and frees
 * it with free(), then encodes it (PNG/JPEG) with the vendored stb writers, so
 * encoding does not depend on which formats this Raylib build supports.
 * Defined in render.c under USE_GRAPHICS; guard the call with #ifdef
 * USE_GRAPHICS so the reference vanishes when the renderer is not built. */
unsigned char* graphics_render_rgba(const Expr* graphics_expr, int w, int h,
                                    int* out_w, int* out_h);

#endif /* MATHILDA_GRAPHICS_EXPORT_H */
