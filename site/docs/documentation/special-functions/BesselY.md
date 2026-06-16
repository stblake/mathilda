# BesselY

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
BesselY[n, z]
    gives the Bessel function of the second kind Y_n(z), the solution of
    z^2 y'' + z y' + (z^2 - n^2) y = 0 singular at the origin.
Y_0(0) = -Infinity, Y_n(0) = ComplexInfinity for integer n != 0; Y_n has
a logarithmic branch point at 0 and a branch cut along the negative real
z axis, with Y_{-n} = (-1)^n Y_n for integer n. Real and complex order
and argument evaluate numerically at machine or arbitrary (MPFR)
precision; D[BesselY[n, z], z] = (BesselY[n-1, z] - BesselY[n+1, z])/2.
Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= BesselY[0, 2.5]
Out[1]= 0.49807

In[2]:= D[BesselY[n, x], x]
Out[2]= 1/2 (BesselY[-1 + n, x] - BesselY[1 + n, x])
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
