#ifndef MATHILDA_GRAPHICS_RENDER3D_H
#define MATHILDA_GRAPHICS_RENDER3D_H

#include "expr.h"

/* Render `graphics3d_expr` (a Graphics3D[...] expression, as built by
 * Plot3D) in a blocking Raylib window with an orbit camera; control
 * returns to the caller once the window is closed. Mirrors show.h's
 * graphics_show() for the 2D case -- see render3d.c. Does not take
 * ownership of `graphics3d_expr`. */
void graphics3d_show(const Expr* graphics3d_expr);

/* Embedded 3D content, for Animate/Manipulate's per-frame content region.
 * Opaque (no Raylib types here) so callers can include this header and hold
 * a Graphics3DEmbedState* regardless of USE_GRAPHICS. */
typedef struct Graphics3DEmbedState Graphics3DEmbedState;

Graphics3DEmbedState* graphics3d_embed_state_new(void);
void graphics3d_embed_state_free(Graphics3DEmbedState* state);

/* Snaps the orbit camera back to its home view (the view derived from the
 * content's bounding box the first time this state was rendered). No-op if
 * nothing has been rendered with this state yet. */
void graphics3d_embed_state_reset_view(Graphics3DEmbedState* state);

/* Renders `graphics3d_expr` (a Graphics3D[...] expression) into the
 * sub-rectangle (rx,ry,rw,rh) of the currently open window. Must be called
 * inside BeginDrawing()/EndDrawing(). Re-parses and re-bakes the expression
 * every call (the caller re-evaluates its body fresh each frame, so there is
 * no stable expression identity to cache against -- mirrors
 * graphics_render_in_region's lack of caching in the 2D case). The orbit
 * camera (drag to rotate, scroll to zoom, right/middle-drag to pan) persists
 * in `state` across calls and is only gated to (rx,ry,rw,rh) for *starting*
 * a drag; an in-progress drag continues regardless of mouse position, same
 * as the caller's own slider-drag state. Does not take ownership of
 * `graphics3d_expr` or `state`. */
void graphics3d_render_in_region(const Expr* graphics3d_expr,
                                  float rx, float ry, float rw, float rh,
                                  Graphics3DEmbedState* state);

#endif /* MATHILDA_GRAPHICS_RENDER3D_H */
