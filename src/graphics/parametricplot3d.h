#ifndef MATHILDA_GRAPHICS_PARAMETRICPLOT3D_H
#define MATHILDA_GRAPHICS_PARAMETRICPLOT3D_H

#include "expr.h"

/* ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, opts...]
 * ParametricPlot3D[{{fx1,fy1,fz1}, ...}, {t, tmin, tmax}, opts...]
 * ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, {u, umin, umax}, opts...]
 *
 * HoldAll. One-iterator form: adaptively samples the parametric 3D space
 * curve(s) over [tmin,tmax], builds a Graphics3D[{Line[...], ...}, opts]
 * expression, auto-displays it, and returns it.
 * Two-iterator form: samples a PlotPoints x PlotPoints grid of (t,u) pairs,
 * maps each to (x,y,z) via the body, and emits Polygon[] quads, producing
 * Graphics3D[{Polygon[...], ...}, opts]. */
Expr* builtin_parametricplot3d(Expr* res);

#endif /* MATHILDA_GRAPHICS_PARAMETRICPLOT3D_H */
