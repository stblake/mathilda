#ifndef MATHILDA_GRAPHICS_PLOT3D_H
#define MATHILDA_GRAPHICS_PLOT3D_H

#include "expr.h"

/* Plot3D[f, {x,xmin,xmax}, {y,ymin,ymax}, opts...] -- HoldAll. Samples f
 * over a uniform (x,y) grid (adaptively refined by MaxRecursion; see
 * plot3d.c), builds a Graphics3D[{Polygon[...], ...}, opts] expression,
 * and returns it. The REPL front end renders any top-level Graphics3D[...]
 * result, exactly as it does for Graphics[...] (see repl.c). */
Expr* builtin_plot3d(Expr* res);

#endif /* MATHILDA_GRAPHICS_PLOT3D_H */
