#ifndef MATHILDA_GRAPHICS_RENDER3D_H
#define MATHILDA_GRAPHICS_RENDER3D_H

#include "expr.h"

/* Render `graphics3d_expr` (a Graphics3D[...] expression, as built by
 * Plot3D) in a blocking Raylib window with an orbit camera; control
 * returns to the caller once the window is closed. Mirrors show.h's
 * graphics_show() for the 2D case -- see render3d.c. Does not take
 * ownership of `graphics3d_expr`. */
void graphics3d_show(const Expr* graphics3d_expr);

#endif /* MATHILDA_GRAPHICS_RENDER3D_H */
