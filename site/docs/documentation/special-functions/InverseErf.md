# InverseErf

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
InverseErf[s]
    gives the inverse error function: the z solving s = Erf[z].
InverseErf[z0, s]
    gives the inverse of the generalized error function Erf[z0, z].
InverseErf[0] = 0, InverseErf[1] = Infinity, InverseErf[-1] = -Infinity.
Odd in s. Numerical values are given only for real s in [-1, 1], at
machine or arbitrary (MPFR) precision; D[InverseErf[z], z] =
(Sqrt[Pi]/2) E^(InverseErf[z]^2). Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `InverseErf[0] = 0`, `InverseErf[1] = Infinity`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
