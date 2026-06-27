#ifndef MATHILDA_GRAPHICS_LISTPLOT_H
#define MATHILDA_GRAPHICS_LISTPLOT_H

#include "expr.h"

/* ListPlot[data, opts...] -- HoldAll. Plots explicit data as a point
 * (scatter) plot, or a connected line when Joined -> True, returning a
 * Graphics[...] object rendered by the existing engine (render.c).
 *
 * Data forms:
 *   {y1, ..., yn}             -> points {i, yi} (x from DataRange).
 *   {{x1,y1}, ..., {xn,yn}}   -> the given points (scatter).
 *   {data1, data2, ...}       -> several datasets, one palette colour each.
 *
 * Elements not of a plottable form are treated as missing and skipped. */
Expr* builtin_listplot(Expr* res);

#endif /* MATHILDA_GRAPHICS_LISTPLOT_H */
