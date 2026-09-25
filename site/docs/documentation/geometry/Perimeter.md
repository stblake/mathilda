# Perimeter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Perimeter[Polygon[{{x1, y1}, ...}]] gives the perimeter of a simple 2D polygon: the sum of its edge lengths, including the closing edge. Exact coordinates give an exact (possibly symbolic, e.g. 2 + Sqrt[2]) result; Real coordinates give a machine-precision result.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Sqrt](../../arithmetic/Sqrt/)

- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)
