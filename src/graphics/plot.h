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
 * nothing is plottable over the window.
 *
 * [yclip_lo, yclip_hi] is the visible y-band the adaptive refinement should
 * track (the current zoom's extent); pass a degenerate band (lo >= hi) to
 * sample over the curve's full extent. */
Expr* plot_resample(const Expr* graphics_expr, double xmin, double xmax,
                    double yclip_lo, double yclip_hi);

/* Shared with listplot.c. Routes a (possibly symbolic-but-numeric, e.g.
 * 2 Pi / Sqrt[2]) expression through N[] and extracts a finite double.
 * Returns false if the result isn't a finite real. */
bool numericize_bound(Expr* e, double* out);

/* Shared with listplot.c. Returns a freshly owned RGBColor[...] directive
 * from Mathematica's default plot palette (ColorData[97]), cycling for
 * i beyond the table length. */
Expr* palette_color(size_t i);

/* Shared with listplot.c. Builds a continuous Filling strip between the
 * polyline through (xs[i], ys[i]) and the horizontal `baseline`: one quad per
 * segment, split into two triangles where a segment crosses the baseline (so
 * render.c's triangle-fan Polygon fill stays correct). Returns a malloc'd
 * Expr* array (caller owns the array and the Exprs); *out_count is the shape
 * count, 0 (and NULL returned) when n < 2. */
Expr** gfx_build_fill_quads(const double* xs, const double* ys, size_t n,
                            double baseline, size_t* out_count);

#endif /* MATHILDA_GRAPHICS_PLOT_H */
