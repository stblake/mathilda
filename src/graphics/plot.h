#ifndef MATHILDA_GRAPHICS_PLOT_H
#define MATHILDA_GRAPHICS_PLOT_H

#include "expr.h"

/* Plot[f, {x, xmin, xmax}, opts...] -- HoldAll. Adaptively samples f over
 * [xmin, xmax] (see sampling.h), builds a Graphics[{Line[...], ...}, opts]
 * expression, auto-displays it (graphics_show, see show.h), and returns
 * the Graphics[...] expression. */
Expr* builtin_plot(Expr* res);

#endif /* MATHILDA_GRAPHICS_PLOT_H */
