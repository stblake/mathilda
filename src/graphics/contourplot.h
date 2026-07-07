/* contourplot.h — ContourPlot[f, {x,xmin,xmax}, {y,ymin,ymax}, opts...] */
#ifndef MATHILDA_GRAPHICS_CONTOURPLOT_H
#define MATHILDA_GRAPHICS_CONTOURPLOT_H

#include "expr.h"

/* HoldAll builtin: f and the iterator specs must not be pre-evaluated. */
Expr* builtin_contourplot(Expr* res);

#endif /* MATHILDA_GRAPHICS_CONTOURPLOT_H */
