# Erf

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Erf[z]
    gives the error function erf(z) = (2/Sqrt[Pi]) Integral_0^z e^(-t^2) dt.
Erf[z0, z1]
    gives the generalized error function erf(z1) - erf(z0).
Erf[0] = 0, Erf[Infinity] = 1, Erf[-Infinity] = -1. An entire function,
odd in z. Real and complex inputs evaluate numerically at machine or
arbitrary (MPFR) precision; D[Erf[z], z] = (2/Sqrt[Pi]) E^(-z^2). Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `Erf[0] = 0`, `Erf[Infinity] = 1`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
