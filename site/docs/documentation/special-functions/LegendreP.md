# LegendreP

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LegendreP[n, x]
    gives the Legendre polynomial P_n(x).
LegendreP[n, m, x] gives the associated Legendre function P_n^m(x).
LegendreP[n, m, a, x] gives the Legendre function of type a (a in
{1, 2, 3}, default 1). Integer n yields the explicit polynomial; a
non-integer order with an inexact argument evaluates numerically at
machine or arbitrary (MPFR) precision, real or complex. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= LegendreP[3, x]
Out[1]= -3/2 x + 5/2 x^3

In[2]:= LegendreP[10, 2, x]
Out[2]= (1 - x^2) (3465/128 - 45045/32 x^2 + 675675/64 x^4 - 765765/32 x^6 + 2078505/128 x^8)
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
