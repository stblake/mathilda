#ifndef MATHILDA_GRAPHICS_PARAMETRICPLOT_H
#define MATHILDA_GRAPHICS_PARAMETRICPLOT_H

#include "expr.h"

/* ParametricPlot[{fx, fy}, {t, tmin, tmax}, opts...]
 * ParametricPlot[{{fx1,fy1}, {fx2,fy2}, ...}, {t, tmin, tmax}, opts...]
 *
 * HoldAll.  Adaptively samples the parametric curve(s) over [tmin,tmax],
 * builds a Graphics[{Line[...], ...}, opts] expression, auto-displays it,
 * and returns it. */
Expr* builtin_parametricplot(Expr* res);

#endif /* MATHILDA_GRAPHICS_PARAMETRICPLOT_H */
