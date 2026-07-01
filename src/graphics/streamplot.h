#ifndef MATHILDA_GRAPHICS_STREAMPLOT_H
#define MATHILDA_GRAPHICS_STREAMPLOT_H

#include "expr.h"

/* StreamPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
 *
 * HoldAll: the vector field components and iterator specs are not
 * pre-evaluated (x and y have no values at call time).  Integrates
 * streamlines of the 2-D vector field using RK4 from a grid of seed
 * points, builds a Graphics[{Arrow[...], ...}, opts] expression, and
 * returns it (auto-displayed by the REPL front-end).
 *
 * Supported options
 * -----------------
 *  StreamPoints          – Integer n (n×n seed grid) or {pts} explicit
 *                          seeds (default Automatic = 15×15 grid).
 *  StreamScale           – Automatic (moderate length), None (let streams
 *                          run until they leave the domain), or a positive
 *                          real giving the maximum arc-length fraction of
 *                          the domain diagonal (default Automatic ≈ 0.08).
 *  StreamStyle           – A style directive or list applied to every stream
 *                          (e.g. Thickness[0.002], RGBColor[...]).
 *  StreamColorFunction   – A function f[x,y,vx,vy,speed] (or fewer args)
 *                          that returns a color directive, applied per stream
 *                          at its midpoint.
 *  ColorFunction         – Alias for StreamColorFunction.
 *  RegionFunction        – f[x,y] mask; seeds and integration steps outside
 *                          the region are skipped.
 *  PlotLegends           – Automatic / "Expressions" / explicit list.
 *  PlotPoints            – Alias for StreamPoints (integer only).
 *  Standard Graphics options (PlotRange, Axes, AspectRatio, Frame, …)
 *                          pass through to the returned Graphics[...]. */
Expr* builtin_streamplot(Expr* res);

#endif /* MATHILDA_GRAPHICS_STREAMPLOT_H */
