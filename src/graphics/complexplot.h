#ifndef MATHILDA_GRAPHICS_COMPLEXPLOT_H
#define MATHILDA_GRAPHICS_COMPLEXPLOT_H

#include "expr.h"

/* ComplexPlot[f, {z, zmin, zmax}, opts...]  -- HoldAll.
 *
 * Domain-colouring plot of the complex function f over the rectangular
 * region in the complex plane with corners zmin and zmax.  Each grid cell
 * is coloured by Arg(f(z)) via the system thermal ramp (matching DensityPlot
 * and Plot3D defaults), with brightness proportional to |f(z)|/(1+|f(z)|).
 * When ColorFunction is supplied it receives f[re, im] and must return a
 * color directive; ColorFunctionScaling (default True) normalises the
 * argument to [0,1] before calling it.
 *
 * Returns a Graphics[...] expression (auto-displayed by the REPL). */
Expr* builtin_complexplot(Expr* res);

/* ComplexPlot3D[f, {z, zmin, zmax}, opts...]  -- HoldAll.
 *
 * Three-dimensional surface over the same rectangular complex domain:
 * height = |f(z)|, colour = Arg(f(z)) via the thermal ramp (matching Plot3D
 * defaults).  Returns a Graphics3D[...] expression (auto-displayed). */
Expr* builtin_complexplot3d(Expr* res);

#endif /* MATHILDA_GRAPHICS_COMPLEXPLOT_H */
