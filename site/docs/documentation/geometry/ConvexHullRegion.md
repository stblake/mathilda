# ConvexHullRegion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConvexHullRegion[{{x1, y1}, ...}] gives the convex hull of a set of 2D points: a Polygon with the hull vertices in counterclockwise order, a Line for collinear input, or a Point for a single point.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)
