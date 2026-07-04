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

/*
 * graphics3d_to_plotly_json(g)
 *
 * Convert a Graphics3D[...] Expr (as produced by Graph3D) to a Plotly
 * scatter3d scene JSON string. Handled primitives:
 *   Line[List[List[x,y,z], List[x,y,z]]] — a 3D line segment (edge)
 *   Point[List[List[x,y,z], ...]]        — 3D markers (vertices)
 *   RGBColor[r, g, b]                    — current colour
 * Returns a heap-allocated string (caller frees) or NULL.
 */
char* graphics3d_to_plotly_json(const Expr* g);

#endif /* MATHILDA_GRAPHICS_JSON_H */
