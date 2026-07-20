# PolynomialSqrt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolynomialSqrt[p] gives a polynomial s with s^2 == p when p is a perfect square (every non-constant irreducible factor has even multiplicity; the numeric content is carried through Sqrt), and $Failed otherwise. PolynomialSqrt[p, x] treats p as a polynomial in x.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/facpoly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/facpoly.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
