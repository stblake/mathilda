#ifndef MATHILDA_GRAPHICS_PLOT_H
#define MATHILDA_GRAPHICS_PLOT_H

#include "expr.h"

/* Plot[f, {x, xmin, xmax}, opts...] -- HoldAll. Adaptively samples f over
 * [xmin, xmax] (see sampling.h), builds a Graphics[{Line[...], ...}, opts]
 * expression, auto-displays it (graphics_show, see show.h), and returns
 * the Graphics[...] expression. */
Expr* builtin_plot(Expr* res);

/* Re-sample a Plot-produced Graphics object over a new x-window. The
 * renderer calls this on zoom/pan so curves like Sin[1/x^2] stay smooth
 * when magnified, instead of revealing the original coarse sampling.
 *
 * `graphics_expr` must be a Graphics[...] carrying the $PlotResample
 * metadata that builtin_plot embeds; returns a freshly sampled primitive
 * List[...] (caller owns, free with expr_free) for [xmin, xmax], or NULL
 * if the object has no such metadata (e.g. a hand-built Graphics) or
 * nothing is plottable over the window. */
Expr* plot_resample(const Expr* graphics_expr, double xmin, double xmax);

#endif /* MATHILDA_GRAPHICS_PLOT_H */
