# RegionCentroid

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RegionCentroid[Polygon[{{x1, y1}, ...}]] gives the centroid {cx, cy} of a simple 2D polygon with nonzero area. Exact coordinates give an exact result; Real coordinates give a machine-precision result.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)
