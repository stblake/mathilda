# RegionMember

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RegionMember[Polygon[{{x1, y1}, ...}], {x, y}] gives True if the point lies inside or on the boundary of the simple 2D polygon, and False otherwise.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)
