/*
 * geometry.h -- core computational-geometry builtins (GEO-1).
 *
 * Five heads over 2D data, in the spirit of Wolfram Language 12 core geometry:
 *
 *   Area[Polygon[{{x1,y1},...}]]            shoelace area (exact or machine)
 *   Perimeter[Polygon[{{x1,y1},...}]]       edge-length sum (symbolic-exact or machine)
 *   RegionCentroid[Polygon[{{x1,y1},...}]]  area centroid (exact or machine)
 *   RegionMember[Polygon[...], {x,y}]       boundary-inclusive point-in-polygon
 *   ConvexHullRegion[{{x1,y1},...}]         monotone-chain hull -> Polygon/Line/Point
 *
 * Simple (non-self-intersecting) polygons only; anything out of scope declines
 * (returns NULL) and stays unevaluated. See docs/spec/builtins/geometry.md.
 */
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "expr.h"

void geometry_init(void);

Expr* builtin_area(Expr* res);
Expr* builtin_perimeter(Expr* res);
Expr* builtin_region_centroid(Expr* res);
Expr* builtin_region_member(Expr* res);
Expr* builtin_convex_hull_region(Expr* res);

#endif /* GEOMETRY_H */
