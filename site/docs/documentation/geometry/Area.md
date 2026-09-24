# Area

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Area[Polygon[{{x1, y1}, ...}]] gives the area of a simple 2D polygon. Exact (Integer/Rational) coordinates give an exact result; any Real coordinate gives a machine-precision result. A polygon with fewer than 3 distinct vertices has Undefined area.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)
