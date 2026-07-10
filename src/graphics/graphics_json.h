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
 *   GrayLevel[g]                — sets colour via gray level (0=black, 1=white)
 *   Hue[h], Hue[h,s,v,a]       — sets colour via HSB/HSV
 *   CMYKColor[c,m,y,k]         — sets colour via CMYK
 *   Opacity[n]                  — sets the current fill opacity
 *
 * Skipped: $PlotResample, Rule[...], unknown heads.
 */
char* graphics_to_plotly_json(const Expr* g);

/*
 * graphics3d_to_plotly_json(g)
 *
 * Convert a Graphics3D[...] Expr to a Plotly-compatible JSON string.
 * Returns a heap-allocated string (caller must free) on success,
 * or NULL if the expression cannot be converted.
 *
 * Handled primitives:
 *   Line[List[List[x,y,z], ...]]  — scatter3d trace, mode="lines"
 *   Polygon[List[pts...]]         — mesh3d trace (quads triangulated)
 *   RGBColor[r, g, b]             — sets the current stroke/fill colour
 *   GrayLevel[g]                  — sets colour via gray level (0=black, 1=white)
 *   Hue[h], Hue[h,s,v,a]         — sets colour via HSB/HSV
 *   CMYKColor[c,m,y,k]           — sets colour via CMYK
 *   Opacity[n]                    — sets the current fill opacity
 *
 * Skipped: Rule[...], unknown heads.
 */
char* graphics3d_to_plotly_json(const Expr* g);

#endif /* MATHILDA_GRAPHICS_JSON_H */
