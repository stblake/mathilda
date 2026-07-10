#ifndef MATHILDA_GRAPHICS_POLARPLOT_H
#define MATHILDA_GRAPHICS_POLARPLOT_H

#include "expr.h"

/* PolarPlot[r, {theta, tmin, tmax}, opts...]
 * PolarPlot[{r1, r2, ...}, {theta, tmin, tmax}, opts...]
 *
 * HoldAll: r and the iterator spec are held unevaluated on entry.
 * Converts r(theta) -> {r*Cos[theta], r*Sin[theta]} then delegates to
 * builtin_parametricplot, so all ParametricPlot sampling options apply.
 * Returns a Graphics[{Line[...], ...}, opts] expression. */
Expr* builtin_polarplot(Expr* res);

#endif /* MATHILDA_GRAPHICS_POLARPLOT_H */
