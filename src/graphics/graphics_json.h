#ifndef MATHILDA_GRAPHICS_JSON_H
#define MATHILDA_GRAPHICS_JSON_H

#include "expr.h"

/*
 * graphics_to_plotly_json(g)
 *
 * Convert a Graphics[...] Expr to a Plotly-compatible JSON string.
 * Returns a heap-allocated string (caller must free) on success,
 * or NULL if the expression cannot be converted.
 *
 * Handled primitives:
 *   Line[List[List[x,y], ...]]  — scatter trace, mode="lines"
 *   Polygon[List[pts...]]       — scatter trace, fill="toself"
 *   RGBColor[r, g, b]           — sets the current stroke/fill colour
 *   Opacity[n]                  — sets the current fill opacity
 *
 * Skipped: $PlotResample, Rule[...], unknown heads.
 */
char* graphics_to_plotly_json(const Expr* g);

#endif /* MATHILDA_GRAPHICS_JSON_H */
