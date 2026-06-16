# BesselJ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
BesselJ[n, z]
    gives the Bessel function of the first kind J_n(z), a solution of
    z^2 y'' + z y' + (z^2 - n^2) y = 0 regular at the origin.
J_0(0) = 1, J_n(0) = 0 for integer n != 0. Has a branch cut along the
negative real z axis for non-integer n. Real and complex order and
argument evaluate numerically at machine or arbitrary (MPFR) precision;
D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= BesselJ[0, 5.2]
Out[1]= -0.11029

In[2]:= D[BesselJ[n, x], x]
Out[2]= 1/2 (BesselJ[-1 + n, x] - BesselJ[1 + n, x])
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
