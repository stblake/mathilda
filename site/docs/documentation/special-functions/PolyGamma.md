# PolyGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolyGamma[z]
    gives the digamma function psi(z) (rewritten as PolyGamma[0, z]).
PolyGamma[n, z]
    gives the n-th derivative of the digamma function, psi^(n)(z).
Positive-integer arguments reduce to exact values: psi(m) to a rational
minus EulerGamma, and psi^(n)(m) for odd n to a rational plus a rational
multiple of Pi^(n+1); even orders stay symbolic. Non-positive integer
arguments give ComplexInfinity. Inexact real and complex arguments evaluate
numerically at machine or arbitrary (MPFR) precision. PolyGamma[-1, z] gives
LogGamma[z]. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
