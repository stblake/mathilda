# HarmonicNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HarmonicNumber[n]
    gives the n-th harmonic number H_n = Sum_{i=1}^n 1/i.
HarmonicNumber[n, r]
    gives the order-r harmonic number H_n^(r) = Sum_{i=1}^n 1/i^r.
Non-negative integer n expands to the exact finite sum (a rational for
integer r, an explicit sum for symbolic r); HarmonicNumber[Infinity, r] is
Zeta[r]; a non-positive integer order r gives the Faulhaber polynomial in n.
Inexact arguments evaluate numerically at machine or arbitrary (MPFR)
precision, including complex order, via Zeta[r] - Zeta[r, n+1] (and the
digamma form for r = 1). Listable.
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
