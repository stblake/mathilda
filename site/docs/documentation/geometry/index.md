# Computational geometry

5 built-in function(s) in this category.

- [`Area`](Area.md) — Area[Polygon[{{x1, y1}, ...}]] gives the area of a simple 2D polygon. Exact (Integer/Rational) coordinates give an exact result; any Real coordinate gives a machine-precision result. A polygon with fewer than 3 distinct vertices has Undefined area.  _(Stable)_
- [`ConvexHullRegion`](ConvexHullRegion.md) — ConvexHullRegion[{{x1, y1}, ...}] gives the convex hull of a set of 2D points: a Polygon with the hull vertices in counterclockwise order, a Line for collinear input, or a Point for a single point.  _(Stable)_
- [`Perimeter`](Perimeter.md) — Perimeter[Polygon[{{x1, y1}, ...}]] gives the perimeter of a simple 2D polygon: the sum of its edge lengths, including the closing edge. Exact coordinates give an exact (possibly symbolic, e.g. 2 + Sqrt[2]) result; Real coordinates give a machine-precision result.  _(Stable)_
- [`RegionCentroid`](RegionCentroid.md) — RegionCentroid[Polygon[{{x1, y1}, ...}]] gives the centroid {cx, cy} of a simple 2D polygon with nonzero area. Exact coordinates give an exact result; Real coordinates give a machine-precision result.  _(Stable)_
- [`RegionMember`](RegionMember.md) — RegionMember[Polygon[{{x1, y1}, ...}], {x, y}] gives True if the point lies inside or on the boundary of the simple 2D polygon, and False otherwise.  _(Stable)_
